//! Axisymmetric electrostatics on the original meridional `(r,z)` complex.
//!
//! `E = -grad(V)` and `div(D) = rho_free`, with `D = epsilon E + P`, are
//! integrated against `dV = angular_span * r dr dz`. The existing weighted
//! Whitney Hodge integrates that radial factor exactly; scalar elimination,
//! checked CG, material resolution and field/charge/energy recovery are shared
//! with the planar owner. Cell and surface charge moments also integrate `r`
//! exactly rather than dividing their total equally between vertices.
//!
//! All data are independent of azimuth. The electric/displacement vectors are
//! local `[radial, axial]` components, with no azimuthal component. Integrated
//! charge, capacitance and energy describe the supplied angular span; use TAU
//! for a complete body of revolution. A smaller span is a symmetry sector,
//! not a wedge with independently driven azimuthal faces. The zero-area axis
//! carries the natural regularity condition and cannot act as an electrode.
//!
//! References: COMSOL's [electrostatic equations][1] and [axis boundary][2].
//! This is a fixed-geometry linear dielectric model. It does not infer exterior
//! fields absent from the mesh, nonlinear/dispersive response, conduction,
//! electrostriction, force or mechanical motion from a static permittivity.
//!
//! [1]: https://doc.comsol.com/6.4/doc/com.comsol.help.edis/edis_ug_acdc.7.71.html
//! [2]: https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_acdc.26.29.html

/// Multiple physical electrodes, floating charge and native circuit storage.
pub mod network;

use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    field_storage::planar::{PlanarDielectricLaw, PlanarPermittivity},
    state_point::ResolvedMaterialStatePoint,
};
use fs_phs::circuit::Branch;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::SolveReport;

use super::{
    ElectricField, ElectrostaticChargeSheet, ElectrostaticControls, ElectrostaticDisplacementTrace,
    ElectrostaticError, PlanarCapacitance, PlanarElectrostaticSources, PlanarElectrostatics, poll,
};

/// Axisymmetric sources in original cell/edge order and absolute SI units.
/// Shared edge descriptors denote surfaces of revolution, not planar sheets.
/// Impressed P is in the local `(r,z)` frame and excludes polarization already
/// contained in the supplied material permittivity.
#[derive(Clone, Copy, Debug)]
pub struct AxisymmetricElectrostaticSources<'a> {
    pub free_charge_density_c_m3: &'a [f64],
    pub impressed_polarization_c_m2: Option<&'a [[f64; 2]]>,
    pub charge_sheets: &'a [ElectrostaticChargeSheet],
    pub displacement_traces: &'a [ElectrostaticDisplacementTrace],
}

impl<'a> AxisymmetricElectrostaticSources<'a> {
    /// Volume free charge with zero impressed polarization and homogeneous
    /// outward displacement on every unprescribed physical boundary.
    pub fn volume_charge(free_charge_density_c_m3: &'a [f64]) -> Self {
        Self {
            free_charge_density_c_m3,
            impressed_polarization_c_m2: None,
            charge_sheets: &[],
            displacement_traces: &[],
        }
    }
}

/// Prepared radial Hodge, original geometry and complete dielectric sources.
/// Its private shared scalar core cannot escape through a planar force API.
pub struct AxisymmetricElectrostatics {
    core: PlanarElectrostatics,
}

impl AxisymmetricElectrostatics {
    /// Assemble source-backed isotropic epsilon on an axisymmetric mesh.
    /// Negative radii and nonpositive physical cell measures are refused by
    /// the original `TriComplex2` geometry owner.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElectrostaticError> {
        Self::assemble(mesh, materials, None, limits, cancelled)
    }

    /// Admit scalar or explicitly oriented SPD meridional laws per cell.
    /// The existing two-coordinate material law is interpreted in local
    /// `(r,z)` coordinates: its axes rotate with azimuth. The omitted azimuthal
    /// direction must be decoupled; a fixed Cartesian anisotropy is not silently
    /// revolved into an axisymmetric one. All laws name one meridional frame.
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
        if !matches!(mesh.metric(), Metric2::Axisymmetric { .. }) {
            return Err(ElectrostaticError::Invalid(
                "axisymmetric electrostatics requires an explicit radial metric",
            ));
        }
        Ok(Self {
            core: PlanarElectrostatics::assemble_metric(mesh, materials, laws, limits, cancelled)?,
        })
    }

    /// Original mesh with metre-valued `[radius, axial_coordinate]` vertices.
    pub fn mesh(&self) -> &TriComplex2 {
        self.core.mesh()
    }

    /// Complete retained material source/state bundles in original cell order.
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        self.core.materials()
    }

    /// Original source laws, including the meridional material-axis frames.
    pub fn permittivity_laws(&self) -> &[PlanarDielectricLaw] {
        self.core.permittivity_laws()
    }

    /// Absolute SI permittivity in each cell's local `(r,z)` frame [F/m].
    pub fn permittivity_coefficients(&self) -> &[PlanarPermittivity] {
        self.core.permittivity_coefficients()
    }

    /// Solve volume charge and exterior voltage traces. Every connected
    /// component needs a positive-area voltage trace; an axis line is not one.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricElectricField<'_>, ElectrostaticError> {
        self.solve_with_sources(
            AxisymmetricElectrostaticSources::volume_charge(free_charge_density_c_m3),
            boundary_voltage_v,
            controls,
            cancelled,
        )
    }

    /// Solve fixed volume/surface charge, outward D.n and impressed P with
    /// their exact radial P1 moments. No source is converted into a fictitious
    /// conductivity, conductor constraint or additional dielectric storage.
    pub fn solve_with_sources(
        &self,
        sources: AxisymmetricElectrostaticSources<'_>,
        boundary_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricElectricField<'_>, ElectrostaticError> {
        let field = self.core.solve_with_sources(
            PlanarElectrostaticSources {
                free_charge_density_c_m3: sources.free_charge_density_c_m3,
                impressed_polarization_c_m2: sources.impressed_polarization_c_m2,
                charge_sheets: sources.charge_sheets,
                displacement_traces: sources.displacement_traces,
            },
            boundary_voltage_v,
            controls,
            cancelled,
        )?;
        Ok(AxisymmetricElectricField {
            problem: self,
            field,
        })
    }

    /// Extract the capacitance of two disjoint equipotential exterior traces,
    /// including the original angular span. The unit-voltage solution checks
    /// Q/V against 2U/V² through the same charge and energy admission gates.
    pub fn capacitance(
        &self,
        positive_terminal: &[usize],
        negative_terminal: &[usize],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricCapacitance<'_>, ElectrostaticError> {
        admit_terminal_traces(
            self.mesh(),
            &[positive_terminal, negative_terminal],
            cancelled,
        )?;
        let PlanarCapacitance {
            field,
            capacitance_f,
        } = self
            .core
            .capacitance(positive_terminal, negative_terminal, controls, cancelled)?;
        Ok(AxisymmetricCapacitance {
            field: AxisymmetricElectricField {
                problem: self,
                field,
            },
            capacitance_f,
        })
    }
}

// Admit each conductor separately. A zero-area axis terminal cannot borrow
// physical support from adjacent edges belonging to a different conductor.
// The shared scalar owner still checks disjoint, complete exterior traces.
pub(super) fn admit_terminal_traces(
    mesh: &TriComplex2,
    terminals: &[&[usize]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ElectrostaticError> {
    poll(cancelled)?;
    if terminals.len() < 2 {
        return Err(ElectrostaticError::Invalid(
            "at least two bounded axisymmetric electrode traces required",
        ));
    }
    admit_electrode_traces(mesh, terminals, cancelled)
}

/// General radial Gauss fields may have one grounded physical conductor.
/// Capacitance/network callers retain their separate two-terminal admission.
pub(super) fn admit_electrode_traces(
    mesh: &TriComplex2,
    terminals: &[&[usize]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ElectrostaticError> {
    poll(cancelled)?;
    let n = mesh.vertices().len();
    crate::scalar_potential::electrode_boundary(
        mesh,
        terminals,
        &vec![0.0; terminals.len()],
        cancelled,
    )?;
    let mut owner = vec![None; n];
    for (terminal, vertices) in terminals.iter().enumerate() {
        for &v in *vertices {
            poll(cancelled)?;
            owner[v] = Some(terminal);
        }
    }
    let trace = mesh
        .boundary_trace()
        .map_err(|_| ElectrostaticError::Invalid("invalid axisymmetric electrode trace"))?;
    let mut supported = vec![false; n];
    for edge in &trace.edges {
        poll(cancelled)?;
        let [a, b] = edge.oriented_vertices.map(|v| v as usize);
        if owner[a].is_some()
            && owner[a] == owner[b]
            && mesh
                .edge_measure(edge.global_edge)
                .is_some_and(|area| area > 0.0)
        {
            supported[a] = true;
            supported[b] = true;
        }
    }
    if owner
        .iter()
        .zip(supported)
        .any(|(terminal, supported)| terminal.is_some() && !supported)
    {
        return Err(ElectrostaticError::Invalid(
            "each axisymmetric electrode needs its own physical surface area, not the symmetry axis",
        ));
    }
    poll(cancelled)
}

// A meridional axis edge has positive line length but zero physical surface
// area. It supplies neither an electrode nor a physical potential anchor.
pub(super) fn admit_prescribed_traces(
    mesh: &TriComplex2,
    prescribed: &[Option<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ElectrostaticError> {
    let trace = mesh
        .boundary_trace()
        .map_err(|_| ElectrostaticError::Invalid("invalid axisymmetric voltage boundary trace"))?;
    let mut supported = vec![false; prescribed.len()];
    for edge in &trace.edges {
        poll(cancelled)?;
        let [a, b] = edge.oriented_vertices.map(|v| v as usize);
        if prescribed[a].is_some()
            && prescribed[b].is_some()
            && mesh
                .edge_measure(edge.global_edge)
                .is_some_and(|area| area > 0.0)
        {
            supported[a] = true;
            supported[b] = true;
        }
    }
    if prescribed
        .iter()
        .zip(supported)
        .any(|(value, supported)| value.is_some() && !supported)
    {
        return Err(ElectrostaticError::Invalid(
            "axisymmetric voltage traces require physical surface area, not the symmetry axis",
        ));
    }
    poll(cancelled)
}

/// Accepted axisymmetric field tied to its original radial geometry/sources.
/// No planar mechanical load method is exposed on this distinct field type.
pub struct AxisymmetricElectricField<'a> {
    problem: &'a AxisymmetricElectrostatics,
    field: ElectricField<'a>,
}

impl AxisymmetricElectricField<'_> {
    pub fn problem(&self) -> &AxisymmetricElectrostatics {
        self.problem
    }
    /// Observational absolute potentials [V]; physical recovery uses differences.
    pub fn voltage_v(&self) -> &[f64] {
        self.field.voltage_v()
    }
    /// Potentials relative to each connected component's prescribed gauge [V].
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.field.relative_voltage_v()
    }
    /// Original gauge attached to every vertex [V].
    pub fn reference_voltage_v(&self) -> &[f64] {
        self.field.reference_voltage_v()
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        self.field.free_charge_density_c_m3()
    }
    /// Original cellwise `[P_r, P_z]`, without induced polarization [C/m²].
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 2]]> {
        self.field.impressed_polarization_c_m2()
    }
    /// Cellwise `[E_r, E_z]` [V/m]; the azimuthal component is zero.
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        self.field.electric_field_v_m()
    }
    /// Cellwise `[D_r, D_z]` [C/m²], with original interface discontinuities.
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        self.field.displacement_c_m2()
    }
    /// Induced dielectric storage integrated over each revolved cell [J].
    pub fn cell_energy_j(&self) -> &[f64] {
        self.field.cell_energy_j()
    }
    /// Total storage for the authored angular span [J], not per radian.
    pub fn energy_j(&self) -> f64 {
        self.field.energy_j()
    }
    /// Weak electrode reactions for the authored angular span [C].
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        self.field.nodal_boundary_charge_c()
    }
    pub fn volume_charge_c(&self) -> f64 {
        self.field.volume_charge_c()
    }
    pub fn sheet_charge_c(&self) -> f64 {
        self.field.sheet_charge_c()
    }
    pub fn boundary_charge_c(&self) -> f64 {
        self.field.boundary_charge_c()
    }
    pub fn outward_displacement_c(&self) -> f64 {
        self.field.outward_displacement_c()
    }
    /// Integral E.P for the fixed impressed source [J], not stored energy.
    pub fn polarization_pairing_j(&self) -> f64 {
        self.field.polarization_pairing_j()
    }
    pub fn charge_balance_residual_c(&self) -> f64 {
        self.field.charge_balance_residual_c()
    }
    pub fn maximum_free_charge_residual_c(&self) -> f64 {
        self.field.maximum_free_charge_residual_c()
    }
    pub fn energy_balance_residual_j(&self) -> f64 {
        self.field.energy_balance_residual_j()
    }
    pub fn report(&self) -> Option<&SolveReport> {
        self.field.report()
    }
}

/// Source-resolved axisymmetric capacitor on its original angular span.
pub struct AxisymmetricCapacitance<'a> {
    field: AxisymmetricElectricField<'a>,
    capacitance_f: f64,
}

impl AxisymmetricCapacitance<'_> {
    pub fn unit_voltage_field(&self) -> &AxisymmetricElectricField<'_> {
        &self.field
    }
    pub fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
    /// Native circuit reduction at fixed geometry and material state. Circuit
    /// capacitor energy and this same field storage must be counted only once.
    pub fn circuit_branch(&self) -> Branch {
        Branch::Capacitor {
            farads: self.capacitance_f,
        }
    }
}
