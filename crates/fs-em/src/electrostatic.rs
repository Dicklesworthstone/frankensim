//! Spatial linear dielectric electrostatics from the same material data as
//! circuit capacitance. Solve -div(epsilon grad(V)) = rho_free on a planar
//! P1 mesh using the existing material Hodge and checked-CG owners.
//!
//! Prescribed exterior voltage traces supply electrode charge. Explicit normal
//! displacement traces, internal free sheet charge and fixed impressed
//! polarization use the same weak operator; other exterior edges insulate.
//! Fields are not averaged across dielectric interfaces. Charge reactions and
//! the induced dielectric storage integral E.epsilon.E/2 are
//! checked before publication. These are discrete residual checks, not spatial
//! error bounds. This static solve adds no conduction, dielectric dispersion,
//! breakdown or geometry evolution. [`relaxation`] composes its native storage
//! with Ohmic transport on a retained charge history at constant permittivity.

mod force;
mod sources;
pub use sources::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, PlanarElectrostaticSources,
};
/// Multiple source-resolved electrodes and passive native circuit reduction.
pub mod network;
/// Reversible nonlinear dielectrics, mixed electrode controls and bias tangents.
pub mod nonlinear;
/// Retained dielectric charge, native Ohmic current and physical relaxation heat.
pub mod relaxation;
/// Three-dimensional dielectric weak forms and polarization loads on native P1 tetrahedra.
pub mod tetrahedral;

use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError,
        planar::{PlanarDielectricLaw, PlanarPermittivity},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_phs::circuit::Branch;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, SolveReport};
use fs_sparse::Csr;

use crate::planar::{MagneticError, gradients};
use crate::scalar_potential::{self, PotentialError};

/// Independent algebraic, charge and energy acceptance budgets.
#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticControls {
    pub linear: CheckedCgConfig,
    /// Absolute allowed net charge and free-node charge residual [C].
    pub charge_tolerance_c: f64,
    /// Absolute allowed field/source energy-identity residual [J].
    pub energy_tolerance_j: f64,
}

/// Source, geometry, numerical or physical-balance refusal.
#[derive(Debug)]
pub enum ElectrostaticError {
    Material(FieldStorageError),
    Assembly(WeightedError),
    Invalid(&'static str),
    Cancelled,
    Nonconvergence(SolveReport),
    Balance {
        charge_residual_c: f64,
        maximum_free_charge_residual_c: f64,
        energy_residual_j: f64,
    },
}
impl core::fmt::Display for ElectrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar electrostatics: {self:?}")
    }
}
impl std::error::Error for ElectrostaticError {}

impl From<PotentialError> for ElectrostaticError {
    fn from(error: PotentialError) -> Self {
        match error {
            PotentialError::Invalid(reason) => Self::Invalid(reason),
            PotentialError::Cancelled => Self::Cancelled,
            PotentialError::Nonconvergence(report) => Self::Nonconvergence(report),
        }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ElectrostaticError> {
    if cancelled() {
        Err(ElectrostaticError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, ElectrostaticError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ElectrostaticError::Invalid(
            "electric field arithmetic is not representable",
        ))
    }
}
fn accumulate(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), ElectrostaticError> {
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn geometry(error: MagneticError) -> ElectrostaticError {
    // The existing scalar mesh-gradient and exterior-trace helpers perform
    // geometry/gauge admission only; no magnetic law or material is substituted.
    match error {
        MagneticError::Cancelled => ElectrostaticError::Cancelled,
        MagneticError::Invalid(reason) => ElectrostaticError::Invalid(reason),
        _ => ElectrostaticError::Invalid("invalid scalar-potential mesh"),
    }
}

/// Immutable geometry, complete per-triangle material states and epsilon Hodge.
pub struct PlanarElectrostatics {
    mesh: TriComplex2,
    materials: Vec<ResolvedMaterialStatePoint>,
    permittivity_laws: Vec<PlanarDielectricLaw>,
    permittivity: Vec<PlanarPermittivity>,
    scalar_permittivity_f_m: Option<Vec<f64>>,
    stiffness: Csr,
}

struct ElectricCellResponse {
    electric_field_v_m: Vec<[f64; 2]>,
    displacement_c_m2: Vec<[f64; 2]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
}
impl PlanarElectrostatics {
    /// Assemble d0^T M(epsilon) d0 from the SAME scalar absolute SI
    /// permittivity law admitted by ResolvedCapacitor.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectrostaticError> {
        Self::assemble(mesh, materials, None, limits, cancelled)
    }

    /// Assemble the SAME native dielectric operator with explicit scalar or
    /// oriented source laws per cell. Every oriented law names one common mesh
    /// frame; its material axes need not coincide with any conductivity axes.
    /// Materials retain every source requirement and state coordinate.
    pub fn new_with_laws(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        laws: Vec<PlanarDielectricLaw>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectrostaticError> {
        Self::assemble(mesh, materials, Some(laws), limits, cancelled)
    }

    fn assemble(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        laws: Option<Vec<PlanarDielectricLaw>>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectrostaticError> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. }) {
            return Err(ElectrostaticError::Invalid(
                "this dielectric owner requires a planar metric",
            ));
        }
        if materials.len() != mesh.faces().len()
            || laws
                .as_ref()
                .is_some_and(|values| values.len() != materials.len())
            || mesh.vertices().len() > limits.max_dofs
            || mesh.edges().len() > limits.max_dofs
            || mesh.faces().len() > limits.max_cells
            || mesh
                .faces()
                .len()
                .checked_mul(9)
                .is_none_or(|n| n > limits.max_triplets)
        {
            return Err(ElectrostaticError::Invalid(
                "one material per triangle and bounded assembly extents required",
            ));
        }
        let permittivity_laws =
            laws.unwrap_or_else(|| vec![PlanarDielectricLaw::Isotropic; materials.len()]);
        let mut permittivity = Vec::with_capacity(materials.len());
        let mut weights = Vec::with_capacity(materials.len());
        let mut geometry_frame = None;
        for (material, law) in materials.iter().zip(&permittivity_laws) {
            poll(cancelled)?;
            if let Some(frame) = law.frame() {
                if geometry_frame.is_some_and(|id| id != frame.geometry_frame()) {
                    return Err(ElectrostaticError::Invalid(
                        "all oriented dielectric laws must name one mesh-coordinate frame",
                    ));
                }
                geometry_frame = Some(frame.geometry_frame());
            }
            let epsilon = law
                .resolve(material)
                .map_err(ElectrostaticError::Material)?;
            // Tensor SPD admission stays with the original outward-interval
            // weighted FEEC owner; the solver never substitutes an average.
            weights.push(match epsilon {
                PlanarPermittivity::Scalar(value) => CellWeight2::Scalar(value),
                PlanarPermittivity::Tensor(matrix) => CellWeight2::Tensor(matrix),
            });
            permittivity.push(epsilon);
        }
        let scalar_permittivity_f_m = permittivity
            .iter()
            .map(|value| match value {
                PlanarPermittivity::Scalar(value) => Some(*value),
                PlanarPermittivity::Tensor(_) => None,
            })
            .collect();
        let mass = weighted_mass_matrix_2d(
            &mesh,
            1,
            &weights,
            ELECTRIC_PERMITTIVITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => ElectrostaticError::Cancelled,
            error => ElectrostaticError::Assembly(error),
        })?;
        poll(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&mesh.d0()), mass.matrix());
        for row in 0..stiffness.nrows() {
            poll(cancelled)?;
            if stiffness.row(row).1.iter().any(|v| !v.is_finite()) {
                return Err(ElectrostaticError::Invalid(
                    "unrepresentable dielectric stiffness",
                ));
            }
        }
        Ok(Self {
            mesh,
            materials,
            permittivity_laws,
            permittivity,
            scalar_permittivity_f_m,
            stiffness,
        })
    }
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        &self.materials
    }
    /// Exact original source-law choices, including material/mesh frames.
    pub fn permittivity_laws(&self) -> &[PlanarDielectricLaw] {
        &self.permittivity_laws
    }
    /// Actual scalar/tensor coefficients in original mesh coordinates [F/m].
    pub fn permittivity_coefficients(&self) -> &[PlanarPermittivity] {
        &self.permittivity
    }
    /// Scalar view only when EVERY cell declares an isotropic law. A mixed or
    /// tensor field has no scalar surrogate; use the complete coefficient view.
    pub fn permittivity_f_m(&self) -> Option<&[f64]> {
        self.scalar_permittivity_f_m.as_deref()
    }

    /// Solve prescribed free charge density [C/m3] and exterior voltage [V].
    /// Every connected component requires a positive-length prescribed trace;
    /// a point gauge is not a physical electrode. Unprescribed edges are insulated.
    /// The cell charge is integrated by the original P1 consistent load.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectricField<'_>, ElectrostaticError> {
        self.solve_with_sources(
            PlanarElectrostaticSources::volume_charge(free_charge_density_c_m3),
            boundary_voltage_v,
            controls,
            cancelled,
        )
    }

    /// Solve with explicit volume/sheet free charge, outward displacement and
    /// fixed impressed polarization on the original mesh. The same sourced
    /// epsilon, scalar solve, electrode reactions and energy gates are reused.
    /// Every component still needs a complete prescribed voltage trace.
    pub fn solve_with_sources(
        &self,
        sources: PlanarElectrostaticSources<'_>,
        boundary_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectricField<'_>, ElectrostaticError> {
        poll(cancelled)?;
        if boundary_voltage_v.len() > self.mesh.vertices().len()
            || !positive(controls.linear.tolerance)
            || controls.linear.tolerance >= 1.0
            || !positive(controls.charge_tolerance_c)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(ElectrostaticError::Invalid(
                "complete finite charges and positive solve tolerances required",
            ));
        }
        let prescribed =
            scalar_potential::prescribed_boundary(&self.mesh, boundary_voltage_v, cancelled)?;
        let loads = sources::assemble(&self.mesh, sources, &prescribed, cancelled)?;
        let solved = scalar_potential::solve_relative(
            &self.stiffness,
            &loads.rhs,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        let mut field = self.recover_field(
            solved.relative,
            solved.references,
            sources.free_charge_density_c_m3,
            &loads.rhs,
            &prescribed,
            finite(loads.volume_charge_c + loads.sheet_charge_c - loads.outward_displacement_c)?,
            controls.charge_tolerance_c,
            controls.energy_tolerance_j,
            solved.report,
            cancelled,
        )?;
        // The reaction calculation above subtracts ALL authored source loads.
        // Polarization changes D, but its externally maintained source is not
        // assigned the induced dielectric energy E.epsilon.E/2.
        sources::retain(&mut field, sources, loads, cancelled)?;
        poll(cancelled)?;
        Ok(field)
    }

    // Share actual E, D and storage with the retained charge owner without
    // imposing a separate prescribed-volume-charge equilibrium on that history.
    fn cell_response(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectricCellResponse, ElectrostaticError> {
        if voltage.len() != self.mesh.vertices().len() || voltage.iter().any(|v| !v.is_finite()) {
            return Err(ElectrostaticError::Invalid(
                "complete finite native potentials required",
            ));
        }
        let mut electric_field_v_m = Vec::with_capacity(self.materials.len());
        let mut displacement_c_m2 = Vec::with_capacity(self.materials.len());
        let mut cell_energy_j = Vec::with_capacity(self.materials.len());
        let mut energy_j = 0.0;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let gradients =
                gradients(face.map(|v| self.mesh.vertices()[v as usize])).map_err(geometry)?;
            let mut e = [0.0; 2];
            for local in 1..3 {
                // A constant voltage offset must not create an electric field.
                let delta = finite(voltage[face[local] as usize] - voltage[face[0] as usize])?;
                for axis in 0..2 {
                    e[axis] = finite((-delta).mul_add(gradients[local][axis], e[axis]))?;
                }
            }
            let displacement = self.permittivity[cell].displacement_c_m2(e);
            let d = [finite(displacement[0])?, finite(displacement[1])?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ElectrostaticError::Invalid(
                    "missing dielectric cell measure",
                ))?;
            let energy = finite(0.5 * volume * e[0].mul_add(d[0], e[1] * d[1]))?;
            if energy < 0.0 || (e != [0.0; 2] && energy == 0.0) {
                return Err(ElectrostaticError::Invalid(
                    "field energy is not representably positive",
                ));
            }
            energy_j = finite(energy_j + energy)?;
            electric_field_v_m.push(e);
            displacement_c_m2.push(d);
            cell_energy_j.push(energy);
        }
        Ok(ElectricCellResponse {
            electric_field_v_m,
            displacement_c_m2,
            cell_energy_j,
            energy_j,
        })
    }

    // Native solves and admitted electrode-basis superpositions use this same
    // gradient, source-reaction and physical-energy recovery.
    #[allow(clippy::too_many_arguments)]
    fn recover_field(
        &self,
        relative_voltage: Vec<f64>,
        reference_voltage: Vec<f64>,
        free_charge_density_c_m3: &[f64],
        rhs: &[f64],
        prescribed: &[Option<f64>],
        volume_charge_c: f64,
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectricField<'_>, ElectrostaticError> {
        let n = self.mesh.vertices().len();
        let ElectricCellResponse {
            electric_field_v_m,
            displacement_c_m2,
            cell_energy_j,
            energy_j,
        } = self.cell_response(&relative_voltage, cancelled)?;
        let mut nodal_boundary_charge_c = vec![None; n];
        let mut boundary_charge_c = 0.0;
        let mut boundary_charge_correction = 0.0;
        let mut maximum_free_charge_residual_c: f64 = 0.0;
        let mut source_pairing_j = 0.0;
        let mut source_pairing_correction = 0.0;
        for v in 0..n {
            poll(cancelled)?;
            let (cols, vals) = self.stiffness.row(v);
            let mut applied = 0.0;
            let mut applied_correction = 0.0;
            for (&w, &value) in cols.iter().zip(vals) {
                // The full scalar operator annihilates constants. Apply it to
                // edge differences before multiplication, as the field-gradient
                // owner does, rather than subtracting enormous biased charges.
                let difference = finite(relative_voltage[w] - relative_voltage[v])?;
                accumulate(
                    &mut applied,
                    &mut applied_correction,
                    finite(value * difference)?,
                )?;
            }
            let applied = finite(applied + applied_correction)?;
            let reaction = finite(applied - rhs[v])?;
            let supplied = if prescribed[v].is_some() {
                nodal_boundary_charge_c[v] = Some(reaction);
                accumulate(
                    &mut boundary_charge_c,
                    &mut boundary_charge_correction,
                    reaction,
                )?;
                // rhs + reaction is exactly the applied operator; do not lose
                // that sum by subtracting and then re-adding a large cell load.
                applied
            } else {
                maximum_free_charge_residual_c = maximum_free_charge_residual_c.max(reaction.abs());
                rhs[v]
            };
            // Gauge balance is independent on each disconnected component.
            // Pair with its own voltage differences; an arbitrary common bias
            // must not amplify roundoff in charge balance into spurious energy.
            accumulate(
                &mut source_pairing_j,
                &mut source_pairing_correction,
                finite(relative_voltage[v] * supplied)?,
            )?;
        }
        let boundary_charge_c = finite(boundary_charge_c + boundary_charge_correction)?;
        let source_pairing_j = finite(source_pairing_j + source_pairing_correction)?;
        let charge_balance_residual_c = finite(volume_charge_c + boundary_charge_c)?;
        let energy_balance_residual_j = finite(energy_j - 0.5 * source_pairing_j)?;
        if charge_balance_residual_c.abs() > charge_tolerance_c
            || maximum_free_charge_residual_c > charge_tolerance_c
            || energy_balance_residual_j.abs() > energy_tolerance_j
        {
            return Err(ElectrostaticError::Balance {
                charge_residual_c: charge_balance_residual_c,
                maximum_free_charge_residual_c,
                energy_residual_j: energy_balance_residual_j,
            });
        }
        // Absolute voltages are an observational view only. Reconstruct them
        // after every physical calculation, preserving authored traces exactly.
        let mut voltage = Vec::with_capacity(n);
        for v in 0..n {
            poll(cancelled)?;
            voltage.push(match prescribed[v] {
                Some(value) => value,
                None => finite(reference_voltage[v] + relative_voltage[v])?,
            });
        }
        poll(cancelled)?;
        Ok(ElectricField {
            problem: self,
            free_charge_density_c_m3: free_charge_density_c_m3.to_vec(),
            impressed_polarization_c_m2: None,
            free_charge_sheets: Vec::new(),
            displacement_traces: Vec::new(),
            source_loaded_cells: Vec::new(),
            sheet_charge_c: 0.0,
            outward_displacement_c: 0.0,
            polarization_pairing_j: 0.0,
            voltage_v: voltage,
            relative_voltage_v: relative_voltage,
            reference_voltage_v: reference_voltage,
            electric_field_v_m,
            displacement_c_m2,
            cell_energy_j,
            energy_j,
            nodal_boundary_charge_c,
            volume_charge_c,
            boundary_charge_c,
            charge_balance_residual_c,
            maximum_free_charge_residual_c,
            energy_balance_residual_j,
            report,
        })
    }

    /// Extract a two-terminal capacitor from a one-volt, zero-volume-charge solve.
    /// Each terminal is a nonempty union of complete exterior mesh edges.
    /// Remaining exterior edges are insulated, so their physical placement is
    /// part of the capacitor geometry. Fringing inside this domain is solved;
    /// exterior/end fields absent from the mesh are not guessed.
    pub fn capacitance(
        &self,
        positive_terminal: &[usize],
        negative_terminal: &[usize],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarCapacitance<'_>, ElectrostaticError> {
        poll(cancelled)?;
        let excitation_v = 1.0;
        let boundary = self.terminal_boundary(
            &[positive_terminal, negative_terminal],
            &[excitation_v, 0.0],
            cancelled,
        )?;
        let field = self.solve(
            &vec![0.0; self.materials.len()],
            &boundary,
            controls,
            cancelled,
        )?;
        let mut terminal_charge = 0.0;
        for &v in positive_terminal {
            poll(cancelled)?;
            terminal_charge = finite(terminal_charge + field.nodal_boundary_charge_c[v].unwrap())?;
        }
        // C=2U/V^2, independently checked against Q/V in the same units.
        let capacitance_f = finite(2.0 * field.energy_j / (excitation_v * excitation_v))?;
        let charge_capacitance_f = finite(terminal_charge / excitation_v)?;
        let allowance_f = controls.charge_tolerance_c / excitation_v
            + 2.0 * controls.energy_tolerance_j / (excitation_v * excitation_v);
        if !positive(capacitance_f)
            || !positive(capacitance_f.recip())
            || terminal_charge <= 0.0
            || (charge_capacitance_f - capacitance_f).abs() > allowance_f
        {
            return Err(ElectrostaticError::Invalid(
                "nonpositive or unbalanced two-terminal capacitance",
            ));
        }
        poll(cancelled)?;
        Ok(PlanarCapacitance {
            field,
            capacitance_f,
        })
    }

    // The scalar and multi-electrode reductions admit exactly the same traces.
    fn terminal_boundary(
        &self,
        terminals: &[&[usize]],
        voltages_v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<(usize, f64)>, ElectrostaticError> {
        scalar_potential::terminal_boundary(&self.mesh, terminals, voltages_v, cancelled)
            .map_err(Into::into)
    }
}

/// Accepted electrostatic field tied to its exact geometry and material states.
pub struct ElectricField<'a> {
    problem: &'a PlanarElectrostatics,
    free_charge_density_c_m3: Vec<f64>,
    impressed_polarization_c_m2: Option<Vec<[f64; 2]>>,
    free_charge_sheets: Vec<ElectrostaticChargeSheet>,
    displacement_traces: Vec<ElectrostaticDisplacementTrace>,
    source_loaded_cells: Vec<bool>,
    sheet_charge_c: f64,
    outward_displacement_c: f64,
    polarization_pairing_j: f64,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    reference_voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 2]>,
    displacement_c_m2: Vec<[f64; 2]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    nodal_boundary_charge_c: Vec<Option<f64>>,
    volume_charge_c: f64,
    boundary_charge_c: f64,
    charge_balance_residual_c: f64,
    maximum_free_charge_residual_c: f64,
    energy_balance_residual_j: f64,
    report: Option<SolveReport>,
}
impl ElectricField<'_> {
    pub fn problem(&self) -> &PlanarElectrostatics {
        self.problem
    }
    /// Absolute nodal voltage view. Prescribed traces are retained exactly;
    /// small induced voltages can round away under a large common reference.
    /// All physical fields and balances use `relative_voltage_v` instead.
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    /// Actual solved potentials relative to each component's prescribed gauge.
    /// These retain charge-induced voltages below the spacing of absolute V.
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    /// The component's prescribed reference at each original vertex [V].
    /// Absolute voltage is reference plus relative voltage, subject to rounding.
    pub fn reference_voltage_v(&self) -> &[f64] {
        &self.reference_voltage_v
    }
    /// Original prescribed volume-charge source per cell [C/m3], retained even
    /// when positive and negative charges cancel in the total inventory.
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    /// Cell E=-grad(V) [V/m], in the original planar frame.
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.electric_field_v_m
    }
    /// Cell D=epsilon E+P_impressed [C/m2], without interface averaging.
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        &self.displacement_c_m2
    }
    /// Induced dielectric storage E.epsilon.E/2 per cell [J]. Fixed impressed
    /// polarization has no inferred source energy, hysteresis or work law.
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Weak reactions K*V-f [C] on prescribed nodes, positive for charge
    /// supplied to the electrode. Sum a complete terminal's nodes for its charge.
    /// Free nodes carry None, not a fictitious physical charge.
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        &self.nodal_boundary_charge_c
    }
    pub fn volume_charge_c(&self) -> f64 {
        self.volume_charge_c
    }
    pub fn boundary_charge_c(&self) -> f64 {
        self.boundary_charge_c
    }
    pub fn charge_balance_residual_c(&self) -> f64 {
        self.charge_balance_residual_c
    }
    pub fn maximum_free_charge_residual_c(&self) -> f64 {
        self.maximum_free_charge_residual_c
    }
    /// U - V.(volume + sheet - outward displacement + polarization load +
    /// electrode charge)/2 [J], using each component's prescribed reference.
    /// The polarization load pairs as -integral(E.P_impressed), not free charge.
    pub fn energy_balance_residual_j(&self) -> f64 {
        self.energy_balance_residual_j
    }
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}

/// Source-bound unit-voltage field and its energy/charge-checked capacitance.
pub struct PlanarCapacitance<'a> {
    field: ElectricField<'a>,
    capacitance_f: f64,
}
impl PlanarCapacitance<'_> {
    pub fn unit_voltage_field(&self) -> &ElectricField<'_> {
        &self.field
    }
    pub fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
    /// Lower into the existing Kirchhoff/pHS circuit owner. Keep this object
    /// alongside the branch to retain geometry/source identity. This reduction
    /// assumes a fixed dielectric snapshot; changing its state in a live circuit
    /// needs explicit parameter-work accounting, not replacement of one number.
    pub fn circuit_branch(&self) -> Branch {
        Branch::Capacitor {
            farads: self.capacitance_f,
        }
    }
}
