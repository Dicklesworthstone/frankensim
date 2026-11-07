//! Self-consistent inertial polarization on planar and axisymmetric fields.
//!
//! Each original material pole retains polarization and polarization current.
//! Its passive midpoint condensation enters the existing native potential
//! operator; recovery keeps both restoring and kinetic energy. The field and
//! poles exchange work through the same discrete Gauss/current balance.
//! No local sinusoidal steady-state or positive real frequency-domain epsilon
//! assumption replaces the time-domain state.
//!
//! This is fixed-geometry, isotropic electroquasistatics: magnetic induction,
//! propagating Maxwell waves and temperature-dependent coefficient work are
//! outside this owner. It is explicitly insulating; pole friction is its only
//! material heat. Midpoint conserves lossless discrete energy, but its phase
//! error requires time refinement through the highest participating resonance.

use super::{
    DielectricFieldCandidate, DielectricFieldState, ElectrostaticControls, Result,
    TriangleGeometry, boundary, invalid,
    kernel::{DebyeGeometry, Field, Prepared},
    sealed, triangle_geometry,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::field_storage::dielectric::lorentz::{LorentzState, ResolvedLorentzMaterial};
use fs_rep_mesh::TriComplex2;
use std::sync::Arc;

/// Original meridional geometry and all inertial material states.
#[derive(Clone)]
pub struct MeridionalLorentzDielectric<const AXISYMMETRIC: bool> {
    geometry: Arc<TriangleGeometry>,
    core: Field<2, ResolvedLorentzMaterial>,
}
pub type PlanarLorentzDielectric = MeridionalLorentzDielectric<false>;
/// Radial/axial polarization vectors on the original revolved metric.
pub type AxisymmetricLorentzDielectric = MeridionalLorentzDielectric<true>;

/// Positive endpoint operator bound to one original geometry/material owner.
/// Reusing it at fixed dt does not repeat native Hodge assembly.
pub struct PreparedLorentzFieldStep<const AXISYMMETRIC: bool = false> {
    core: Prepared<2, ResolvedLorentzMaterial>,
}
pub type PreparedAxisymmetricLorentzFieldStep = PreparedLorentzFieldStep<true>;
impl<const AXISYMMETRIC: bool> PreparedLorentzFieldStep<AXISYMMETRIC> {
    pub fn time_step_s(&self) -> f64 {
        self.core.time_step_s()
    }
}

pub type LorentzFieldCandidate<const AXISYMMETRIC: bool = false> =
    DielectricFieldCandidate<MeridionalLorentzDielectric<AXISYMMETRIC>, 2>;
pub type AxisymmetricLorentzFieldCandidate = LorentzFieldCandidate<true>;

impl<const AXISYMMETRIC: bool> MeridionalLorentzDielectric<AXISYMMETRIC> {
    /// Retain the original field, pole positions/currents and weak free charge.
    /// Each connected component needs an original physical voltage trace.
    /// Unprescribed nodal free charge is retained throughout the trajectory.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedLorentzMaterial>,
        voltage_v: Vec<f64>,
        states: Vec<Vec<LorentzState<2>>>,
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
            states,
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
    pub fn materials(&self) -> &[ResolvedLorentzMaterial] {
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
    /// Original cell/pole order; current is an independent inertial state.
    pub fn polarization_states(&self) -> &[Vec<LorentzState<2>>] {
        self.core.states()
    }
    pub fn nodal_charge_c(&self) -> &[f64] {
        self.core.nodal_charge_c()
    }
    /// Electric plus all original restoring and kinetic polarization storage.
    pub fn cell_energy_j(&self) -> &[f64] {
        self.core.cell_energy_j()
    }
    pub fn cell_volumes_m3(&self) -> &[f64] {
        self.core.cell_volumes_m3()
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
    pub fn prepare_step(
        &self,
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedLorentzFieldStep<AXISYMMETRIC>> {
        Ok(PreparedLorentzFieldStep {
            core: self.core.prepare_step(dt_s, cancelled)?,
        })
    }
    /// Immutable trial under linearly ramped physical voltage traces. Charge,
    /// work and energy checks complete before publishing any endpoint or heat.
    pub fn step_candidate(
        &self,
        prepared: &PreparedLorentzFieldStep<AXISYMMETRIC>,
        endpoint_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LorentzFieldCandidate<AXISYMMETRIC>> {
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

impl<const AXISYMMETRIC: bool> sealed::Sealed for MeridionalLorentzDielectric<AXISYMMETRIC> {}
impl<const AXISYMMETRIC: bool> DielectricFieldState for MeridionalLorentzDielectric<AXISYMMETRIC> {
    type Material = ResolvedLorentzMaterial;
    fn materials(&self) -> &[ResolvedLorentzMaterial] {
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
