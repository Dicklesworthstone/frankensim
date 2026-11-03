//! Three-dimensional source-bound electrostatics on the original P1 complex.
//!
//! Solve -div(epsilon grad(V))=rho, with E=-grad(V) and D=epsilon E.
//! The tetrahedral Whitney Hodge, geometry, trace integration and checked CG
//! are the SAME owners used by the other spatial fields. Permittivity is an
//! absolute SI dielectric source, never a conductivity or a relative constant.
//!
//! Prescribed D.n is positive OUT of the dielectric; its weak load is negative.
//! Dirichlet reactions are free electrode charge. A pure-Neumann component
//! must balance its volume charge against its outward displacement before a
//! deterministic numerical gauge is introduced. That gauge supplies no charge.
//!
//! This is a fixed-geometry, linear, nondispersive snapshot. Unlisted exterior
//! faces have D.n=0; a finite enclosing box is not an infinite-domain model.
//! Cell fields retain material jumps. Reported weak residuals and raw interface
//! jumps are not continuum error bounds or an equilibrated H(div) certificate.
//! Formulation reference: <https://mfem.org/electromagnetics/#electrostatics>.

mod force;
mod terminals;
pub use terminals::SpatialCapacitance;

use super::{
    ScalarBoundary, SpatialGeometry, SpatialMagneticError, dot, prepare_geometry,
    prepare_scalar_boundary_traces_with_sources, solve_scalar_system,
};
use crate::electrostatic::ElectrostaticControls;
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{CellWeight, WeightedAssemblyLimits, WeightedError, weighted_mass_matrix},
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError,
        spatial::{SpatialDielectricLaw, SpatialPermittivity},
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::TetComplex;
use fs_solver::SolveReport;
use fs_sparse::Csr;
use std::sync::Arc;

/// P1 voltage trace on an original exterior face, in sorted face-vertex order.
/// Shared vertices must receive the same potential from every incident trace.
#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticVoltageFace {
    pub face: usize,
    pub voltage_v: [f64; 3],
}

/// Uniform OUTWARD electric displacement D.n [C/m2] on an exterior face.
/// This is a distributed Neumann condition, not a floating equipotential metal.
#[derive(Clone, Copy, Debug)]
pub struct ElectrostaticDisplacementFace {
    pub face: usize,
    pub outward_displacement_c_m2: f64,
}

/// Unlisted faces have zero normal displacement. Voltage and displacement
/// traces cannot occupy the same face. Disconnected components retain their
/// own reference potential and charge-compatibility check.
#[derive(Clone, Debug, Default)]
pub struct SpatialElectrostaticBoundary {
    pub voltage: Vec<ElectrostaticVoltageFace>,
    pub displacement: Vec<ElectrostaticDisplacementFace>,
}

/// No partially admitted field returns on a material, solve or balance error.
#[derive(Debug)]
pub enum SpatialElectrostaticError {
    /// Shared volume geometry, boundary, compatibility or checked-solver error.
    Native(SpatialMagneticError),
    Material(FieldStorageError),
    Assembly(WeightedError),
    Invalid(&'static str),
    Cancelled,
    Balance {
        charge_residual_c: f64,
        maximum_component_charge_residual_c: f64,
        maximum_free_charge_residual_c: f64,
        energy_residual_j: f64,
        hodge_energy_residual_j: f64,
    },
}
impl From<SpatialMagneticError> for SpatialElectrostaticError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Native(other),
        }
    }
}
impl core::fmt::Display for SpatialElectrostaticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral electrostatics: {self:?}")
    }
}
impl std::error::Error for SpatialElectrostaticError {}
type Result<T> = core::result::Result<T, SpatialElectrostaticError>;

fn invalid(reason: &'static str) -> SpatialElectrostaticError {
    SpatialElectrostaticError::Invalid(reason)
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialElectrostaticError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable electrostatic arithmetic"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let adjusted = finite(value - *correction)?;
    let next = finite(*sum + adjusted)?;
    *correction = finite((next - *sum) - adjusted)?;
    *sum = next;
    Ok(())
}

/// Immutable complete material states and their original tetrahedral operator.
/// Geometry and element gradients are prepared once and reused across sources
/// and electrode excitations. Scalar and oriented tensor cells may coexist.
pub struct SpatialElectrostatics {
    geometry: Arc<SpatialGeometry>,
    materials: Vec<ResolvedMaterialStatePoint>,
    laws: Vec<SpatialDielectricLaw>,
    permittivity: Vec<SpatialPermittivity>,
    stiffness: Csr,
}

impl SpatialElectrostatics {
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if tetrahedra.len() != materials.len() || tetrahedra.len() > limits.max_cells {
            return Err(invalid(
                "one material per tetrahedron within the explicit cell budget required",
            ));
        }
        let laws = vec![SpatialDielectricLaw::Isotropic; materials.len()];
        Self::new_with_laws(positions_m, tetrahedra, materials, laws, limits, cancelled)
    }

    pub fn new_with_laws(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedMaterialStatePoint>,
        laws: Vec<SpatialDielectricLaw>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if tetrahedra.len() != materials.len() || laws.len() != materials.len() {
            return Err(invalid(
                "one complete dielectric source and law per tetrahedron required",
            ));
        }
        let geometry = Arc::new(prepare_geometry(
            positions_m,
            tetrahedra,
            limits,
            cancelled,
        )?);
        let mut permittivity = Vec::with_capacity(materials.len());
        let mut weights = Vec::with_capacity(materials.len());
        for (material, law) in materials.iter().zip(&laws) {
            check(cancelled)?;
            let epsilon = law
                .resolve(material)
                .map_err(SpatialElectrostaticError::Material)?;
            weights.push(match epsilon {
                SpatialPermittivity::Scalar(value) => CellWeight::Scalar(value),
                SpatialPermittivity::Tensor(matrix) => CellWeight::Tensor(matrix),
            });
            permittivity.push(epsilon);
        }
        let mass = weighted_mass_matrix(
            &geometry.mesh,
            &geometry.geometry,
            1,
            &weights,
            ELECTRIC_PERMITTIVITY_DIMS,
            limits,
            cancelled,
        )
        .map_err(|error| match error {
            WeightedError::Cancelled => SpatialElectrostaticError::Cancelled,
            other => SpatialElectrostaticError::Assembly(other),
        })?;
        check(cancelled)?;
        let stiffness = stiffness(&incidence_to_csr(&geometry.mesh.d0()), mass.matrix());
        for row in 0..stiffness.nrows() {
            check(cancelled)?;
            for &value in stiffness.row(row).1 {
                finite(value)?;
            }
        }
        Ok(Self {
            geometry,
            materials,
            laws,
            permittivity,
            stiffness,
        })
    }

    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.positions_m
    }
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        &self.materials
    }
    pub fn permittivity_laws(&self) -> &[SpatialDielectricLaw] {
        &self.laws
    }
    pub fn permittivity_coefficients(&self) -> &[SpatialPermittivity] {
        &self.permittivity
    }
    pub fn cell_volume_m3(&self, cell: usize) -> Option<f64> {
        self.geometry.geometry.vol_signed.get(cell).map(|v| v.abs())
    }
    pub fn exterior_faces(&self) -> impl Iterator<Item = usize> + '_ {
        self.geometry
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, face)| (face.cells.len() == 1).then_some(i))
    }
    pub fn exterior_normal(&self, face: usize) -> Option<[f64; 3]> {
        self.geometry
            .faces
            .get(face)
            .filter(|f| f.cells.len() == 1)
            .map(|f| f.outward_normal)
    }
    pub fn face_area_m2(&self, face: usize) -> Option<f64> {
        self.geometry.faces.get(face).map(|f| f.area_m2)
    }

    /// Solve fixed per-cell free volume charge and the explicitly supplied
    /// face traces. P1 volume and face loads use exact constant-source moments.
    /// Field recovery keeps component-relative potentials throughout; adding
    /// a large common voltage does not pollute charge reactions or stored energy.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'_>> {
        let (volume_load, traces) =
            self.prepare_boundary(free_charge_density_c_m3, boundary, controls, cancelled)?;
        let rhs = volume_load
            .iter()
            .zip(&traces.flux_load_wb)
            .map(|(&volume, &surface)| finite(volume + surface))
            .collect::<Result<Vec<_>>>()?;
        let (relative, report) = solve_scalar_system(
            &self.stiffness,
            &rhs,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?;
        self.recover_field(
            relative,
            free_charge_density_c_m3,
            volume_load,
            traces,
            boundary,
            controls,
            report,
            cancelled,
        )
    }

    fn prepare_boundary(
        &self,
        rho: &[f64],
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, ScalarBoundary)> {
        check(cancelled)?;
        if rho.len() != self.mesh().tets.len()
            || rho.iter().any(|v| !v.is_finite())
            || boundary.voltage.len() > self.geometry.faces.len()
            || boundary.displacement.len() > self.geometry.faces.len()
            || !positive(controls.linear.tolerance)
            || controls.linear.tolerance >= 1.0
            || !positive(controls.charge_tolerance_c)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(invalid(
                "complete finite volume sources, bounded traces and positive tolerances required",
            ));
        }
        let g = &self.geometry;
        let mut volume_load = vec![0.0; g.mesh.vertex_count];
        let mut corrections = vec![0.0; g.mesh.vertex_count];
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let load = finite(rho[cell] * g.geometry.vol_signed[cell].abs() / 4.0)?;
            for &v in tet {
                add(
                    &mut volume_load[v as usize],
                    &mut corrections[v as usize],
                    load,
                )?;
            }
        }
        let traces = prepare_scalar_boundary_traces_with_sources(
            &g.mesh,
            &g.faces,
            &g.components,
            boundary.voltage.iter().map(|b| (b.face, b.voltage_v)),
            boundary
                .displacement
                .iter()
                .map(|b| (b.face, -b.outward_displacement_c_m2)),
            Some(&volume_load),
            controls.charge_tolerance_c,
            cancelled,
        )?;
        Ok((volume_load, traces))
    }

    #[allow(clippy::too_many_arguments)]
    fn recover_field(
        &self,
        relative: Vec<f64>,
        rho: &[f64],
        volume_load: Vec<f64>,
        traces: ScalarBoundary,
        boundary: &SpatialElectrostaticBoundary,
        controls: ElectrostaticControls,
        report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'_>> {
        let ScalarBoundary {
            physical_prescribed,
            flux_load_wb: surface_load,
            gauges,
            references,
            ..
        } = traces;
        let g = &self.geometry;
        let n = g.mesh.vertex_count;
        let mut nodal_boundary_charge_c = vec![None; n];
        let mut component_charge = vec![0.0; n];
        let mut component_correction = vec![0.0; n];
        let (mut volume_charge_c, mut volume_correction) = (0.0, 0.0);
        let (mut boundary_charge_c, mut boundary_correction) = (0.0, 0.0);
        let (mut pairing_j, mut pairing_correction) = (0.0, 0.0);
        let (mut hodge_j, mut hodge_correction) = (0.0, 0.0);
        let mut maximum_free_charge_residual_c = 0.0_f64;
        for v in 0..n {
            check(cancelled)?;
            let (columns, values) = self.stiffness.row(v);
            let (mut applied, mut correction) = (0.0, 0.0);
            for (&w, &coefficient) in columns.iter().zip(values) {
                add(
                    &mut applied,
                    &mut correction,
                    finite(coefficient * finite(relative[w] - relative[v])?)?,
                )?;
            }
            add(
                &mut hodge_j,
                &mut hodge_correction,
                finite(0.5 * relative[v] * applied)?,
            )?;
            let charge = if physical_prescribed[v] {
                // Includes prescribed displacement on adjacent faces ONCE.
                Some(finite(applied - volume_load[v])?)
            } else {
                let residual = finite(applied - finite(volume_load[v] + surface_load[v])?)?;
                maximum_free_charge_residual_c = maximum_free_charge_residual_c.max(residual.abs());
                (surface_load[v] != 0.0).then_some(surface_load[v])
            };
            let boundary_charge = charge.unwrap_or(0.0);
            let supplied = if physical_prescribed[v] {
                // Keep the exact operator reaction instead of subtracting then
                // re-adding a potentially much larger prescribed volume load.
                applied
            } else {
                finite(volume_load[v] + boundary_charge)?
            };
            add(
                &mut pairing_j,
                &mut pairing_correction,
                finite(0.5 * relative[v] * supplied)?,
            )?;
            add(&mut volume_charge_c, &mut volume_correction, volume_load[v])?;
            add(
                &mut boundary_charge_c,
                &mut boundary_correction,
                boundary_charge,
            )?;
            let component = g.components[v];
            add(
                &mut component_charge[component],
                &mut component_correction[component],
                supplied,
            )?;
            nodal_boundary_charge_c[v] = charge;
        }
        let mut electric_field_v_m = Vec::with_capacity(self.materials.len());
        let mut displacement_c_m2 = Vec::with_capacity(self.materials.len());
        let mut cell_energy_j = Vec::with_capacity(self.materials.len());
        let (mut energy_j, mut energy_correction) = (0.0, 0.0);
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            let mut e = [0.0; 3];
            for local in 1..4 {
                let difference = finite(relative[tet[local] as usize] - relative[tet[0] as usize])?;
                for (axis, out) in e.iter_mut().enumerate() {
                    *out =
                        finite((-difference).mul_add(g.geometry.grads[cell][local][axis], *out))?;
                }
            }
            let d = self.permittivity[cell].displacement_c_m2(e);
            for value in d {
                finite(value)?;
            }
            let energy = finite(0.5 * dot(e, d) * g.geometry.vol_signed[cell].abs())?;
            if energy < 0.0 || (e != [0.0; 3] && energy == 0.0) {
                return Err(invalid(
                    "nonzero dielectric field lost representably positive storage",
                ));
            }
            add(&mut energy_j, &mut energy_correction, energy)?;
            electric_field_v_m.push(e);
            displacement_c_m2.push(d);
            cell_energy_j.push(energy);
        }
        let charge_balance_residual_c = finite(volume_charge_c + boundary_charge_c)?;
        let maximum_component_charge_residual_c =
            component_charge.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let energy_balance_residual_j = finite(energy_j - pairing_j)?;
        let hodge_energy_residual_j = finite(energy_j - hodge_j)?;
        if charge_balance_residual_c.abs() > controls.charge_tolerance_c
            || maximum_component_charge_residual_c > controls.charge_tolerance_c
            || maximum_free_charge_residual_c > controls.charge_tolerance_c
            || energy_balance_residual_j.abs() > controls.energy_tolerance_j
            || hodge_energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(SpatialElectrostaticError::Balance {
                charge_residual_c: charge_balance_residual_c,
                maximum_component_charge_residual_c,
                maximum_free_charge_residual_c,
                energy_residual_j: energy_balance_residual_j,
                hodge_energy_residual_j,
            });
        }
        let mut interface_displacement_jumps_c = Vec::new();
        for (face, data) in g.faces.iter().enumerate() {
            check(cancelled)?;
            if let [left, right] = data.cells.as_slice() {
                let jump = std::array::from_fn(|axis| {
                    displacement_c_m2[*left][axis] - displacement_c_m2[*right][axis]
                });
                interface_displacement_jumps_c
                    .push((face, finite(dot(jump, data.outward_normal) * data.area_m2)?));
            }
        }
        let mut voltage_v = relative
            .iter()
            .enumerate()
            .map(|(v, &value)| finite(value + references[g.components[v]]))
            .collect::<Result<Vec<_>>>()?;
        // Preserve caller-supplied physical traces bit-for-bit for reporting;
        // recovered fields continue to use the unrounded relative solution.
        for trace in &boundary.voltage {
            for (&v, &value) in g.mesh.faces[trace.face].iter().zip(&trace.voltage_v) {
                voltage_v[v as usize] = value;
            }
        }
        check(cancelled)?;
        Ok(SpatialElectricField {
            problem: self,
            free_charge_density_c_m3: rho.to_vec(),
            voltage_v,
            relative_voltage_v: relative,
            electric_field_v_m,
            displacement_c_m2,
            cell_energy_j,
            energy_j,
            nodal_boundary_charge_c,
            nodal_volume_charge_c: volume_load,
            nodal_applied_surface_charge_c: surface_load,
            volume_charge_c,
            boundary_charge_c,
            charge_balance_residual_c,
            maximum_component_charge_residual_c,
            maximum_free_charge_residual_c,
            energy_balance_residual_j,
            hodge_energy_residual_j,
            interface_displacement_jumps_c,
            gauges,
            report,
        })
    }
}

/// Accepted field retaining its exact source/geometry owner. Relative voltages
/// preserve a small physical signal even when an absolute output rounds under
/// a very large common bias. Field, charge and energy use the relative solution.
pub struct SpatialElectricField<'a> {
    problem: &'a SpatialElectrostatics,
    free_charge_density_c_m3: Vec<f64>,
    voltage_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    electric_field_v_m: Vec<[f64; 3]>,
    displacement_c_m2: Vec<[f64; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    nodal_boundary_charge_c: Vec<Option<f64>>,
    nodal_volume_charge_c: Vec<f64>,
    nodal_applied_surface_charge_c: Vec<f64>,
    volume_charge_c: f64,
    boundary_charge_c: f64,
    charge_balance_residual_c: f64,
    maximum_component_charge_residual_c: f64,
    maximum_free_charge_residual_c: f64,
    energy_balance_residual_j: f64,
    hodge_energy_residual_j: f64,
    interface_displacement_jumps_c: Vec<(usize, f64)>,
    gauges: Vec<usize>,
    report: Option<SolveReport>,
}
impl SpatialElectricField<'_> {
    pub fn problem(&self) -> &SpatialElectrostatics {
        self.problem
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.electric_field_v_m
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 3]] {
        &self.displacement_c_m2
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// Free electrode charge, including any supplied face-displacement charge
    /// at shared trace vertices. Free unforced nodes carry None, including gauges.
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        &self.nodal_boundary_charge_c
    }
    pub fn nodal_volume_charge_c(&self) -> &[f64] {
        &self.nodal_volume_charge_c
    }
    /// Signed free surface load = -integral(N_i D.n dA) [C].
    pub fn nodal_applied_surface_charge_c(&self) -> &[f64] {
        &self.nodal_applied_surface_charge_c
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
    pub fn maximum_component_charge_residual_c(&self) -> f64 {
        self.maximum_component_charge_residual_c
    }
    pub fn maximum_free_charge_residual_c(&self) -> f64 {
        self.maximum_free_charge_residual_c
    }
    pub fn energy_balance_residual_j(&self) -> f64 {
        self.energy_balance_residual_j
    }
    pub fn hodge_energy_residual_j(&self) -> f64 {
        self.hodge_energy_residual_j
    }
    /// Raw P1 normal-D discontinuity integrated on each interior face [C].
    pub fn interface_displacement_jumps_c(&self) -> &[(usize, f64)] {
        &self.interface_displacement_jumps_c
    }
    pub fn gauge_vertices(&self) -> &[usize] {
        &self.gauges
    }
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}
