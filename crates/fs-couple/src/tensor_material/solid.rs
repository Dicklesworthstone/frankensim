//! Original 3D viscoelastic solids drive the existing full-memory mechanics.
//! Geometry, material quadrature, consistent mass and force/velocity mapping
//! all come from one retained source basis. No strain reconstruction, damping
//! fit, second eigensolve or alternate time integrator is introduced.

/// Physical element loss deposited on the original finite thermal timeline.
pub mod thermal;
/// Exact source queries and memory-preserving operator rebuilding.
pub mod source;

use super::*;
use fs_solid::linear3::visco::{ProjectedViscoSolid, ViscoElementError};

#[derive(Debug)]
pub enum SolidError {
    Source(ViscoElementError),
    Dynamics(TensorMaterialError),
    MaterialSource(fs_material::visco_tensor::source::TensorSourceError),
}
impl core::fmt::Display for SolidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => write!(f, "tensor solid source: {e}"),
            Self::Dynamics(e) => write!(f, "tensor solid dynamics: {e}"),
            Self::MaterialSource(e) => write!(f, "tensor solid material source: {e}"),
        }
    }
}
impl std::error::Error for SolidError {}
impl From<TensorMaterialError> for SolidError {
    fn from(e: TensorMaterialError) -> Self { Self::Dynamics(e) }
}
impl From<ViscoElementError> for SolidError {
    fn from(e: ViscoElementError) -> Self { Self::Source(e) }
}
impl From<fs_material::visco_tensor::source::TensorSourceError> for SolidError {
    fn from(e: fs_material::visco_tensor::source::TensorSourceError) -> Self { Self::MaterialSource(e) }
}
type SolidResult<T> = core::result::Result<T, SolidError>;

/// One fixed source solid/basis and its existing full-memory dynamics owner.
/// All original material factor rows survive; the dense capacity never removes
/// a point, mode or branch. Coefficients are fixed during each native interval;
/// source-bound continuation can rebuild them at accepted interval boundaries.
pub struct TensorSolidBody<'s> {
    projection: ProjectedViscoSolid<'s>,
    body: TensorMaterialBody,
    sources: Option<Vec<source::SolidMaterialSource>>,
}
impl<'s> TensorSolidBody<'s> {
    pub fn new(
        projection: ProjectedViscoSolid<'s>, budget: TensorMaterialBudget, gate: &CancelGate,
    ) -> SolidResult<Self> {
        poll(gate)?;
        // Preflight before copying immutable material contributions.
        full_3d_extent(projection.coordinates(), projection.material_points(), budget)?;
        let body = TensorMaterialBody::from_mass_matrix(
            projection.coordinates(), projection.mass_matrix_kg().to_vec(),
            projection.material_points().to_vec(), budget, gate,
        )?;
        Ok(Self { projection, body, sources:None })
    }
    pub fn projection(&self) -> &ProjectedViscoSolid<'s> { &self.projection }
    pub fn body(&self) -> &TensorMaterialBody { &self.body }

    /// Coordinate displacement/velocity and constitutive prehistory are explicit.
    /// The original source basis maps them to physical nodal motion.
    pub fn history<'a>(
        &'a self, displacement_m: &[f64], velocity_m_s: &[f64],
        memory: InitialMemory, gate: &CancelGate,
    ) -> SolidResult<TensorSolidHistory<'a, 's>> {
        let inner = self.body.history(displacement_m, velocity_m_s, memory, gate)?;
        Ok(TensorSolidHistory { solid: self, inner })
    }
}

#[derive(Clone, Debug)]
pub struct SolidStep {
    /// Unmodified exact-ZOH record in the retained coordinate order.
    pub dynamics: StepRecord,
    /// Independently integrated Kelvin+Maxwell heat per original element [J].
    pub element_heat_j: Vec<f64>,
    pub heat_attribution_residual_j: f64,
    /// Step-mean nodal velocity from the SAME force-conjugate port output.
    /// dt * sum(f_node dot mean_velocity_node) equals the reported port work,
    /// up to the original propagation/projection arithmetic.
    pub mean_nodal_velocity_m_s: Vec<[f64; 3]>,
}

#[derive(Clone)]
pub struct TensorSolidHistory<'a, 's> {
    solid: &'a TensorSolidBody<'s>,
    inner: TensorMaterialHistory<'a>,
}
impl TensorSolidHistory<'_, '_> {
    pub fn material_history(&self) -> &TensorMaterialHistory<'_> { &self.inner }
    pub fn nodal_displacement_m(&self, gate: &CancelGate) -> SolidResult<Vec<[f64; 3]>> {
        Ok(self.solid.projection.expand_nodal(self.inner.displacement_m(), gate)?)
    }
    pub fn nodal_velocity_m_s(&self, gate: &CancelGate) -> SolidResult<Vec<[f64; 3]>> {
        Ok(self.solid.projection.expand_nodal(&self.inner.velocity_m_s(), gate)?)
    }

    /// Apply held nodal forces through the SAME source basis. Material heat
    /// is attributed to its actual element, not to output modes or a residual.
    /// The explicit absolute heat allowance gates the actual element handoff.
    /// Failed source/heat/expansion/cancellation checks leave the whole history.
    pub fn advance(
        &mut self, nodal_force_n: &[[f64; 3]], dt_s: f64,
        absolute_heat_attribution_tolerance_j: f64, gate: &CancelGate,
    ) -> SolidResult<SolidStep> {
        poll(gate)?;
        let force = self.solid.projection.project_nodal_force(nodal_force_n, gate)?;
        let mut candidate = self.inner.clone();
        let step = candidate.advance_with_material_heat(
            &force, dt_s, absolute_heat_attribution_tolerance_j, gate,
        )?;
        let mut element_heat = vec![0.0; self.solid.projection.element_count()];
        for (&element, &heat) in self.solid.projection.point_elements().iter().zip(&step.point_heat_j) {
            element_heat[element] += heat;
        }
        if !finite(&element_heat) { return Err(bad("solid element heat accumulation overflow").into()); }
        let heat_residual = element_heat.iter().sum::<f64>()-step.dynamics.dissipated;
        if !heat_residual.is_finite() || heat_residual.abs() > absolute_heat_attribution_tolerance_j {
            return Err(bad("solid element heat attribution exceeds its absolute allowance").into());
        }
        let velocity = self.solid.projection.expand_nodal(&step.dynamics.y, gate)?;
        poll(gate)?;
        self.inner = candidate;
        Ok(SolidStep {
            dynamics: step.dynamics, element_heat_j: element_heat,
            heat_attribution_residual_j: heat_residual,
            mean_nodal_velocity_m_s: velocity,
        })
    }
}
