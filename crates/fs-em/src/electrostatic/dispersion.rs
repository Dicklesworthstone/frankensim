//! Spatial multipole dielectric memory on planar and axisymmetric P1 meshes.
//!
//! Triangle and tetrahedron owners use one shared evolution kernel. Geometry
//! adapters preserve their original native Hodge, physical trace and metric.
//! Local polarization is eliminated before the potential solve; equal-sized
//! steps reuse the prepared operator. Geometry and material conditions are fixed.
//! The opt-in `resolved-leaky-dielectric` lane includes independently sourced
//! Ohmic current in charge continuity and recovers separate Joule/polarization
//! heat. No conductivity floor, hidden heat reservoir or induction is supplied.
//! See `kernel` for physical charge, energy and loss recovery.
//! Total-current formulation: <https://www.comsol.com/multiphysics/electroquasistatics-theory>.

use super::{
    ElectrostaticControls, ElectrostaticError, accumulate, finite, geometry, poll, positive,
};
use crate::{planar::gradients, scalar_potential};
use fs_feec::{
    assembly::{incidence_to_csr, stiffness},
    weighted::{WeightedAssemblyLimits, WeightedError},
    weighted2::{CellWeight2, weighted_mass_matrix_2d},
};
use fs_material::field_storage::dielectric::bulk::ResolvedDebyeMaterial;
#[cfg(feature = "resolved-leaky-dielectric")]
use fs_material::field_storage::dielectric::leaky::ResolvedLeakyDielectric;
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::SolveReport;
use fs_sparse::Csr;
use kernel::{DebyeGeometry, Field, Prepared};
use std::sync::Arc;

pub(crate) mod kernel;
/// Inertial polarization on the same native potential/charge evolution.
pub mod lorentz;
mod source;
pub use source::DielectricMaterialState;
pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Original immutable field identity and physical cell accounts used by
/// receiving multiphysics owners. Implementations are native admitted fields.
pub trait DielectricFieldState: sealed::Sealed + Clone {
    type Material: DielectricMaterialState;
    fn materials(&self) -> &[Self::Material];
    fn cell_volumes_m3(&self) -> &[f64];
    fn conductivity_s_m(&self) -> &[f64];
    /// None denotes the explicitly insulating constructor; Some retains the
    /// exact source name, including a sourced zero conductivity.
    fn conductivity_properties(&self) -> &[Option<String>];
    fn stored_energy_j(&self) -> f64;
    fn time_s(&self) -> f64;
    fn same_problem(&self, other: &Self) -> bool;
}
type Result<T> = core::result::Result<T, ElectrostaticError>;
fn invalid(reason: &'static str) -> ElectrostaticError {
    ElectrostaticError::Invalid(reason)
}

struct TriangleGeometry {
    mesh: TriComplex2,
    gradients: Vec<[[f64; 2]; 3]>,
    volumes_m3: Vec<f64>,
}
impl DebyeGeometry<2> for TriangleGeometry {
    fn vertex_count(&self) -> usize {
        self.mesh.vertices().len()
    }
    fn cell_count(&self) -> usize {
        self.mesh.faces().len()
    }
    fn volumes_m3(&self) -> &[f64] {
        &self.volumes_m3
    }
    fn assemble(
        &self,
        coefficients: &[f64],
        dimensions: Dims,
        limits: WeightedAssemblyLimits,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Csr> {
        let weights: Vec<_> = coefficients
            .iter()
            .copied()
            .map(CellWeight2::Scalar)
            .collect();
        let mass =
            weighted_mass_matrix_2d(&self.mesh, 1, &weights, dimensions, limits, &mut || {
                cancelled()
            })
            .map_err(|error| match error {
                WeightedError::Cancelled => ElectrostaticError::Cancelled,
                other => ElectrostaticError::Assembly(other),
            })?;
        poll(&mut || cancelled())?;
        let matrix = stiffness(&incidence_to_csr(&self.mesh.d0()), mass.matrix());
        for row in 0..matrix.nrows() {
            poll(&mut || cancelled())?;
            if matrix.row(row).1.iter().any(|value| !value.is_finite()) {
                return Err(invalid("dielectric field operator is not representable"));
            }
        }
        Ok(matrix)
    }
    fn fields(
        &self,
        voltage: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<[f64; 2]>> {
        let mut fields = Vec::with_capacity(self.mesh.faces().len());
        for (face, g) in self.mesh.faces().iter().zip(&self.gradients) {
            poll(&mut || cancelled())?;
            let mut e = [0.0; 2];
            for local in 1..3 {
                let difference = finite(voltage[face[local] as usize] - voltage[face[0] as usize])?;
                for axis in 0..2 {
                    e[axis] = finite((-difference).mul_add(g[local][axis], e[axis]))?;
                }
            }
            fields.push(e);
        }
        Ok(fields)
    }

    // Same constant-cell P1 moment as the impressed-polarization weak load.
    // The first local entry closes partition of unity exactly before reduction.
    fn polarization_load(
        &self,
        polarization: &[[f64; 2]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        let mut rhs = vec![0.0; self.mesh.vertices().len()];
        let mut corrections = vec![0.0; rhs.len()];
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(&mut || cancelled())?;
            let g = self.gradients[cell];
            let p = polarization[cell];
            let volume = self.volumes_m3[cell];
            let b1 = finite(volume * finite(g[1][0].mul_add(p[0], g[1][1] * p[1]))?)?;
            let b2 = finite(volume * finite(g[2][0].mul_add(p[0], g[2][1] * p[1]))?)?;
            for (&vertex, load) in face.iter().zip([-finite(b1 + b2)?, b1, b2]) {
                accumulate(
                    &mut rhs[vertex as usize],
                    &mut corrections[vertex as usize],
                    load,
                )?;
            }
        }
        for (value, correction) in rhs.iter_mut().zip(corrections) {
            *value = finite(*value + correction)?;
        }
        Ok(rhs)
    }
}

/// Fixed-geometry dielectric field in an explicit meridional metric. Use the
/// named planar/axisymmetric aliases; they cannot exchange prepared operators.
#[derive(Clone)]
pub struct MeridionalDebyeDielectric<const AXISYMMETRIC: bool> {
    geometry: Arc<TriangleGeometry>,
    core: Field<2>,
}
/// Extruded Cartesian P1 dielectric field; original source-compatible owner.
pub type PlanarDebyeDielectric = MeridionalDebyeDielectric<false>;
/// Revolved meridional P1 dielectric field with local radial/axial vectors.
/// All integrated storage, charge and heat include the authored angular span.
pub type AxisymmetricDebyeDielectric = MeridionalDebyeDielectric<true>;

/// Condensed endpoint operator for one immutable geometry/material owner.
pub struct PreparedDebyeStep<const AXISYMMETRIC: bool = false> {
    core: Prepared<2>,
}
pub type PreparedAxisymmetricDebyeStep = PreparedDebyeStep<true>;
impl<const AXISYMMETRIC: bool> PreparedDebyeStep<AXISYMMETRIC> {
    pub fn time_step_s(&self) -> f64 {
        self.core.time_step_s()
    }
}

/// Physical exchanges on one accepted-candidate interval. Source-supplied
/// charge is positive into the device and includes both stored electrode
/// charge change and conduction. Local heat can enter a thermal owner without
/// treating stored electric/polarization energy as heat.
#[derive(Clone, Debug)]
pub struct DebyeFieldStep<const D: usize = 2> {
    pub interval_s: [f64; 2],
    /// Change of stored weak electrode charge only, excluding transported charge.
    pub nodal_boundary_charge_increment_c: Vec<Option<f64>>,
    /// Total supplied charge Delta Q + integrated conduction outflow. This is
    /// the charge paired with midpoint voltage to calculate source work.
    pub nodal_boundary_supplied_charge_c: Vec<Option<f64>>,
    /// Integrated weak conduction outflow at all nodes. Its negative is the
    /// free-charge change at an unprescribed node.
    pub nodal_conduction_outflow_c: Vec<f64>,
    pub polarization_current_density_a_m2: Vec<[f64; D]>,
    pub conduction_current_density_a_m2: Vec<[f64; D]>,
    pub cell_polarization_heat_j: Vec<f64>,
    pub cell_joule_heat_j: Vec<f64>,
    /// Total local polarization plus Joule heat, each included exactly once.
    pub cell_heat_j: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Sum of local endpoint storage increments, without subtracting two large
    /// total energies. This is the change used in the physical balance.
    pub stored_change_j: f64,
    pub polarization_heat_j: f64,
    pub joule_heat_j: f64,
    pub heat_j: f64,
    /// Work components in the same per-component voltage references. Their
    /// sum is physical source work; the separated terms retain small dielectric
    /// exchanges under large DC through-power and are not separate gauge-free
    /// energy accounts.
    pub supplied_displacement_work_j: f64,
    pub supplied_conduction_work_j: f64,
    pub supplied_work_j: f64,
    pub maximum_free_charge_residual_c: f64,
    pub boundary_charge_residual_c: f64,
    pub energy_residual_j: f64,
    pub linear: Option<SolveReport>,
}

/// Shared physical charge, storage and heat accounts for dielectric memory.
/// Lorentz and Debye field candidates use this same account type.
pub use DebyeFieldStep as DielectricFieldStep;

/// Transactional field candidate with independent physical charge/work checks.
pub struct DielectricFieldCandidate<F, const D: usize> {
    pub(crate) field: F,
    pub(crate) step: DebyeFieldStep<D>,
}
pub type DebyeFieldCandidate<const AXISYMMETRIC: bool = false> =
    DielectricFieldCandidate<MeridionalDebyeDielectric<AXISYMMETRIC>, 2>;
pub type AxisymmetricDebyeFieldCandidate = DebyeFieldCandidate<true>;
impl<F, const D: usize> DielectricFieldCandidate<F, D> {
    pub fn field(&self) -> &F {
        &self.field
    }
    pub fn step(&self) -> &DebyeFieldStep<D> {
        &self.step
    }
    pub fn into_parts(self) -> (F, DebyeFieldStep<D>) {
        (self.field, self.step)
    }
}

fn triangle_geometry<const AXISYMMETRIC: bool>(
    mesh: TriComplex2,
    limits: WeightedAssemblyLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Arc<TriangleGeometry>> {
    poll(cancelled)?;
    let correct_metric = if AXISYMMETRIC {
        matches!(mesh.metric(), Metric2::Axisymmetric { .. })
    } else {
        matches!(mesh.metric(), Metric2::Planar { .. })
    };
    if !correct_metric
        || mesh.faces().len() > limits.max_cells
        || mesh.vertices().len() > limits.max_dofs
        || mesh.edges().len() > limits.max_dofs
        || mesh
            .faces()
            .len()
            .checked_mul(9)
            .is_none_or(|n| n > limits.max_triplets)
    {
        return Err(invalid(
            "bounded dielectric mesh must match its explicit planar/axisymmetric owner",
        ));
    }
    let mut cell_gradients = Vec::with_capacity(mesh.faces().len());
    let mut volumes_m3 = Vec::with_capacity(mesh.faces().len());
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        cell_gradients
            .push(gradients(face.map(|v| mesh.vertices()[v as usize])).map_err(geometry)?);
        volumes_m3.push(
            mesh.face_measure(cell)
                .ok_or_else(|| invalid("dielectric cell lacks its physical volume"))?,
        );
    }
    Ok(Arc::new(TriangleGeometry {
        mesh,
        gradients: cell_gradients,
        volumes_m3,
    }))
}
fn boundary<const AXISYMMETRIC: bool>(
    mesh: &TriComplex2,
    voltages: &[(usize, f64)],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Option<f64>>> {
    let prescribed = scalar_potential::prescribed_boundary(mesh, voltages, cancelled)?;
    if AXISYMMETRIC {
        super::axisymmetric::admit_prescribed_traces(mesh, &prescribed, cancelled)?;
    }
    Ok(prescribed)
}
impl<const AXISYMMETRIC: bool> MeridionalDebyeDielectric<AXISYMMETRIC> {
    /// Initialize retained polarization and the corresponding weak free charge.
    /// Every component requires a physical voltage trace. No conductivity is
    /// inferred by this explicitly insulating constructor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedDebyeMaterial>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; 2]>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let geometry = triangle_geometry::<AXISYMMETRIC>(mesh, limits, cancelled)?;
        if voltage_v.len() != geometry.mesh.vertices().len()
            || voltage_vertices.len() > voltage_v.len()
            || voltage_vertices.iter().any(|&v| v >= voltage_v.len())
        {
            return Err(invalid(
                "complete dielectric voltage state and valid trace vertices required",
            ));
        }
        let traces: Vec<_> = voltage_vertices
            .iter()
            .map(|&v| (v, voltage_v[v]))
            .collect();
        let prescribed = boundary::<AXISYMMETRIC>(&geometry.mesh, &traces, cancelled)?;
        let core = Field::new(
            Arc::clone(&geometry) as Arc<dyn DebyeGeometry<2>>,
            materials,
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )?;
        Ok(Self { geometry, core })
    }
    /// Initialize a leaky dielectric from one complete polarization/conductivity
    /// source bundle per cell. Free charge subsequently follows total-current
    /// continuity; a voltage source supplies both displacement and conduction.
    #[cfg(feature = "resolved-leaky-dielectric")]
    #[allow(clippy::too_many_arguments)]
    pub fn new_leaky(
        mesh: TriComplex2,
        materials: Vec<ResolvedLeakyDielectric>,
        voltage_v: Vec<f64>,
        polarization_c_m2: Vec<Vec<[f64; 2]>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let geometry = triangle_geometry::<AXISYMMETRIC>(mesh, limits, cancelled)?;
        if voltage_v.len() != geometry.mesh.vertices().len()
            || voltage_vertices.len() > voltage_v.len()
            || voltage_vertices.iter().any(|&v| v >= voltage_v.len())
        {
            return Err(invalid(
                "complete dielectric voltage state and valid trace vertices required",
            ));
        }
        let traces: Vec<_> = voltage_vertices
            .iter()
            .map(|&v| (v, voltage_v[v]))
            .collect();
        let prescribed = boundary::<AXISYMMETRIC>(&geometry.mesh, &traces, cancelled)?;
        let core = Field::new_leaky(
            Arc::clone(&geometry) as Arc<dyn DebyeGeometry<2>>,
            materials,
            voltage_v,
            polarization_c_m2,
            prescribed,
            time_s,
            limits,
            controls,
            cancelled,
        )?;
        Ok(Self { geometry, core })
    }
    pub fn mesh(&self) -> &TriComplex2 {
        &self.geometry.mesh
    }
    pub fn materials(&self) -> &[ResolvedDebyeMaterial] {
        self.core.materials()
    }
    pub fn voltage_v(&self) -> &[f64] {
        self.core.voltage_v()
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.core.relative_voltage_v()
    }
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        self.core.electric_field_v_m()
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        self.core.displacement_c_m2()
    }
    pub fn polarization_c_m2(&self) -> &[Vec<[f64; 2]>] {
        self.core.states()
    }
    pub fn nodal_charge_c(&self) -> &[f64] {
        self.core.nodal_charge_c()
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        self.core.cell_energy_j()
    }
    pub fn cell_volumes_m3(&self) -> &[f64] {
        self.core.cell_volumes_m3()
    }
    pub fn conductivity_s_m(&self) -> &[f64] {
        self.core.conductivity_s_m()
    }
    pub fn conductivity_properties(&self) -> &[Option<String>] {
        self.core.conductivity_properties()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.core.stored_energy_j()
    }
    pub fn time_s(&self) -> f64 {
        self.core.time_s()
    }
    pub fn same_problem(&self, other: &Self) -> bool {
        self.core.same_problem(&other.core)
    }
    /// Prepare once and reuse at the same dt without repeating Hodge assembly.
    pub fn prepare_step(
        &self,
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedDebyeStep<AXISYMMETRIC>> {
        Ok(PreparedDebyeStep {
            core: self.core.prepare_step(dt_s, cancelled)?,
        })
    }
    /// Advance with linearly ramped voltage traces. No accepted state mutates.
    pub fn step_candidate(
        &self,
        prepared: &PreparedDebyeStep<AXISYMMETRIC>,
        endpoint_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DebyeFieldCandidate<AXISYMMETRIC>> {
        let prescribed =
            boundary::<AXISYMMETRIC>(&self.geometry.mesh, endpoint_voltage_v, cancelled)?;
        let (core, step) =
            self.core
                .step_candidate(&prepared.core, prescribed, controls, cancelled)?;
        Ok(DielectricFieldCandidate {
            field: Self {
                geometry: Arc::clone(&self.geometry),
                core,
            },
            step,
        })
    }
}

impl<const AXISYMMETRIC: bool> sealed::Sealed for MeridionalDebyeDielectric<AXISYMMETRIC> {}
impl<const AXISYMMETRIC: bool> DielectricFieldState for MeridionalDebyeDielectric<AXISYMMETRIC> {
    type Material = ResolvedDebyeMaterial;
    fn materials(&self) -> &[ResolvedDebyeMaterial] {
        self.core.materials()
    }
    fn cell_volumes_m3(&self) -> &[f64] {
        self.core.cell_volumes_m3()
    }
    fn conductivity_s_m(&self) -> &[f64] {
        self.core.conductivity_s_m()
    }
    fn conductivity_properties(&self) -> &[Option<String>] {
        self.core.conductivity_properties()
    }
    fn stored_energy_j(&self) -> f64 {
        self.core.stored_energy_j()
    }
    fn time_s(&self) -> f64 {
        self.core.time_s()
    }
    fn same_problem(&self, other: &Self) -> bool {
        self.core.same_problem(&other.core)
    }
}
