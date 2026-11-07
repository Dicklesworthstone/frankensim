//! Three-dimensional self-consistent Lorentz polarization and charge.
//!
//! Reuses the original tetrahedral Whitney Hodge, physical electrode traces,
//! potential solver and shared charge/work recovery. Each local pole keeps
//! its original polarization and current vectors; only scalar-potential DOFs
//! enter the global solve. Complete restoring/kinetic storage and frictional
//! heat are recovered separately. See the planar Lorentz owner for the fixed
//! geometry, electroquasistatic and time-resolution boundaries.

use super::{
    DebyeGeometry, DielectricFieldCandidate, ElectrostaticControls, Field, Prepared, Result,
    TetrahedralGeometry,
};
use crate::electrostatic::dispersion::{DielectricFieldState, sealed};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::field_storage::dielectric::lorentz::{LorentzState, ResolvedLorentzMaterial};
use fs_rep_mesh::TetComplex;
use std::sync::Arc;

#[derive(Clone)]
pub struct SpatialLorentzDielectric {
    geometry: Arc<TetrahedralGeometry>,
    core: Field<3, ResolvedLorentzMaterial>,
}
pub struct PreparedSpatialLorentzStep {
    core: Prepared<3, ResolvedLorentzMaterial>,
}
impl PreparedSpatialLorentzStep {
    pub fn time_step_s(&self) -> f64 {
        self.core.time_step_s()
    }
}
pub type SpatialLorentzFieldCandidate = DielectricFieldCandidate<SpatialLorentzDielectric, 3>;

impl SpatialLorentzDielectric {
    /// Initialize an insulating field with original polarization and current.
    /// Every component needs a positive-area exterior voltage trace. Free
    /// nodal charge follows the supplied field and remains fixed thereafter.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        materials: Vec<ResolvedLorentzMaterial>,
        voltage_v: Vec<f64>,
        states: Vec<Vec<LorentzState<3>>>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let geometry = TetrahedralGeometry::new(positions_m, tetrahedra, limits, cancelled)?;
        let prescribed = geometry.initial_boundary(&voltage_v, voltage_vertices, cancelled)?;
        let core = Field::new(
            Arc::clone(&geometry) as Arc<dyn DebyeGeometry<3>>,
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
    pub fn mesh(&self) -> &TetComplex {
        &self.geometry.native.mesh
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.geometry.native.positions_m
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
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        self.core.electric_field_v_m()
    }
    pub fn displacement_c_m2(&self) -> &[[f64; 3]] {
        self.core.displacement_c_m2()
    }
    pub fn polarization_states(&self) -> &[Vec<LorentzState<3>>] {
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
    ) -> Result<PreparedSpatialLorentzStep> {
        Ok(PreparedSpatialLorentzStep {
            core: self.core.prepare_step(dt_s, cancelled)?,
        })
    }
    /// Preserve the accepted field on cancellation, a failed solve or a
    /// physical charge/energy refusal. Prepared operators are owner-bound.
    pub fn step_candidate(
        &self,
        prepared: &PreparedSpatialLorentzStep,
        endpoint_voltage_v: &[(usize, f64)],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzFieldCandidate> {
        let prescribed = self.geometry.boundary(endpoint_voltage_v, cancelled)?;
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

impl sealed::Sealed for SpatialLorentzDielectric {}
impl DielectricFieldState for SpatialLorentzDielectric {
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
