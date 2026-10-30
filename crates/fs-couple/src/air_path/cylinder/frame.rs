//! Source-defined elastic frames in the original causal cylinder fluid.
//!
//! Dry frame modes/normalization come from fs-solid. The same source section's
//! outer radius and the native cubic transverse trace define its fluid load.
//! Four positive Gauss points integrate cubic-trace products on each original
//! member or constant-fluid subspan. The reduced inertia is FULL:
//! M = I + sum(w*m_added*b*b^T), not a list of retuned modal frequencies.
//! Each fluid branch retains its own momentum, with resistance
//! w*c*(b^T*v-v_fluid)^2. Original material viscosity keeps all cross terms.
//!
//! This is frozen, small, planar vibration at fixed supports in a stationary
//! Newtonian exterior. Circular prismatic members, local 2-D Stokes flow and
//! explicit modal/spectral truncation apply. No axial/rotary fluid drag,
//! inter-member or inter-region hydrodynamic interactions, thermal evolution,
//! moving supports, cavitation, nonlinear wake or acoustic radiation is inferred.
//! Regional boundaries are prescribed reference data, NOT resolved free surfaces:
//! waves, menisci and cross-section transitions are not modeled. Hollow members
//! receive EXTERIOR loading only; internal fluid needs its own owner.

use super::{CylinderFlow, CylinderFlowResolution};
use crate::acoustic_realize::AcousticRealizeError;
use crate::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget, FrameVibrationError};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_la::factor::{FactorError, cholesky};
use fs_material::fluid_transport::FluidTransport;
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};
use fs_solid::frame2d::Frame2d;

mod mechanics;
mod regions;
pub use regions::FrameFluidRegion;

/// Original member index and original exterior material source; never a radius,
/// mass or stiffness override. Different members may have different fluids.
pub struct FrameFluidBinding<F> {
    pub element: usize,
    pub fluid: F,
}

/// Complete preparation limits, in addition to the native dry-frame budget.
#[derive(Clone, Copy, Debug)]
pub struct FrameFluidBudget {
    /// 2*retained_modes + 4*loaded_regions*spectral_intervals.
    pub max_states: usize,
    /// Maximum entries in EACH full dense state matrix; not total byte usage.
    /// Native factorizations and time steppers also require workspace.
    pub max_dense_entries: usize,
    /// Six scalar projection terms per retained mode per fluid station or
    /// original material resistance port. Checked before projected allocation.
    pub max_projection_terms: usize,
    /// Infinity norm allowance for M*M_inverse-I, finite and in (0,1).
    /// Arithmetic admission only, not a continuum-error certificate.
    pub mass_inverse_tolerance: f64,
}

#[derive(Debug)]
pub enum FrameFluidError {
    Invalid(&'static str),
    Frame(FrameVibrationError),
    Fluid(AcousticRealizeError),
    Factor(FactorError),
    Native(PhsError),
    Balance { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for FrameFluidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material frame/fluid: {self:?}")
    }
}
impl std::error::Error for FrameFluidError {}
impl From<FrameVibrationError> for FrameFluidError {
    fn from(error: FrameVibrationError) -> Self {
        match error { FrameVibrationError::Cancelled => Self::Cancelled, other => Self::Frame(other) }
    }
}
impl From<fs_solid::SolidError> for FrameFluidError {
    fn from(error: fs_solid::SolidError) -> Self { Self::Frame(error.into()) }
}
impl From<PhsError> for FrameFluidError {
    fn from(error: PhsError) -> Self { Self::Native(error) }
}
type Result<T> = core::result::Result<T, FrameFluidError>;
fn invalid(message: &'static str) -> FrameFluidError { FrameFluidError::Invalid(message) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(FrameFluidError::Cancelled) } else { Ok(()) }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() { Ok(x) } else { Err(invalid("unrepresentable frame/fluid coefficient or state")) }
}
fn positive(x: f64) -> Result<f64> {
    finite(x)?;
    if x > 0.0 { Ok(x) } else { Err(invalid("positive frame/fluid coefficient required")) }
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() { return Err(invalid("frame/fluid projection dimensions differ")); }
    finite(a.iter().zip(b).map(|(a,b)| a*b).sum())
}

/// Loaded portion of an unchanged member and its original source. Several
/// entries may refer to disjoint spans of the same original structural element.
pub struct LoadedFrameMember<F> {
    element: usize,
    span: [f64; 2],
    length_m: f64,
    flow: CylinderFlow<F>,
}
impl<F> LoadedFrameMember<F> {
    pub const fn element(&self) -> usize { self.element }
    /// Oriented reference interval in the original member's xi coordinates.
    pub const fn span(&self) -> [f64; 2] { self.span }
    /// Loaded reference length [m], not the whole member length.
    pub const fn length_m(&self) -> f64 { self.length_m }
    pub const fn flow(&self) -> &CylinderFlow<F> { &self.flow }
}

struct Station {
    member: usize,
    weight_m: f64,
    projection: Vec<f64>,
    first_memory: usize,
}
struct MaterialPort {
    projection: Vec<f64>,
    resistance: f64,
}

/// Frozen complete physical component, retaining the original frame and fluids.
/// State order: [dry modal q, LOADED modal momentum p=M*qdot, fluid momenta].
/// Fluid coordinates follow region declaration, then Gauss-point, then original
/// spectral-branch order. They never reset when an interval is accepted.
pub struct MaterialFluidFrame<'a, F> {
    frame: &'a Frame2d,
    basis: FrameVibrationBasis,
    members: Vec<LoadedFrameMember<F>>,
    stations: Vec<Station>,
    material_ports: Vec<MaterialPort>,
    modal_mass: Vec<f64>,
    memory_masses: Vec<f64>,
    inverse_residual: f64,
    system: PortHamiltonian,
}

/// Native proposal plus separate physical dissipation accounts. These losses
/// have NOT been deposited into any thermal owner. `native.x` is the complete
/// next state; publication remains the caller's joint-transaction decision.
pub struct FrameFluidStep {
    pub native: StepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Same order as basis().material_damping().ports(), or empty for elastic frames.
    pub material_port_dissipation_j: Vec<f64>,
    /// One entry per members() region in declaration order, including low tail
    /// and memory-slip loss. Repeated element IDs do not merge these heats.
    pub fluid_member_dissipation_j: Vec<f64>,
    /// Held nodal forces/moments dotted with actual reconstructed displacement change.
    pub nodal_work_j: f64,
    pub energy_residual_j: f64,
}

const GAUSS: [(f64, f64); 4] = [
    (0.06943184420297371, 0.17392742256872693),
    (0.33000947820757187, 0.32607257743127307),
    (0.6699905217924281, 0.32607257743127307),
    (0.9305681557970262, 0.17392742256872693),
];

impl<F> MaterialFluidFrame<'_, F> {
    pub const fn frame(&self) -> &Frame2d { self.frame }
    pub const fn basis(&self) -> &FrameVibrationBasis { &self.basis }
    /// Loaded regions in declaration order. Original element IDs may repeat;
    /// span, source and per-region losses remain separate rather than averaged.
    pub fn members(&self) -> &[LoadedFrameMember<F>] { &self.members }
    /// Full loaded modal inertia in the ORIGINAL dry mass-normalized basis.
    pub fn modal_mass_matrix(&self) -> &[f64] { &self.modal_mass }
    pub fn memory_masses_kg(&self) -> &[f64] { &self.memory_masses }
    pub const fn mass_inverse_residual(&self) -> f64 { self.inverse_residual }
    /// Native force/velocity system. Ports are dry-basis generalized forces;
    /// use the retained basis for physical force/moment projection. Raw native
    /// stepping bypasses step_candidate's source-clock and independent work gates.
    pub const fn system(&self) -> &PortHamiltonian { &self.system }
}

fn project(basis: &FrameVibrationBasis, dofs: &[usize; 6], weights: &[f64; 6]) -> Result<Vec<f64>> {
    basis.modes().iter().map(|mode| finite(dofs.iter().zip(weights)
        .map(|(&dof,&weight)| weight*mode.displacement[dof]).sum())).collect()
}
fn add_outer(matrix: &mut [f64], n: usize, offset: usize, shape: &[f64], coefficient: f64) -> Result<()> {
    for i in 0..shape.len() {
        for k in 0..=i {
            let value = finite(coefficient*shape[i]*shape[k])?;
            let a = (offset+i)*n+offset+k;
            matrix[a] = finite(matrix[a]+value)?;
            if k != i { matrix[(offset+k)*n+offset+i] = matrix[a]; }
        }
    }
    Ok(())
}
fn inverse_mass(mass: &[f64], n: usize, tolerance: f64, gate: &CancelGate) -> Result<(Vec<f64>, f64)> {
    let factor = cholesky(mass, n).map_err(FrameFluidError::Factor)?;
    let mut inverse = vec![0.0; n*n];
    for column in 0..n {
        poll(gate)?;
        let mut rhs = vec![0.0; n];
        rhs[column] = 1.0;
        factor.solve(&mut rhs);
        for row in 0..n { inverse[row*n+column] = finite(rhs[row])?; }
    }
    // Retain an exactly symmetric Hamiltonian; check the actual averaged
    // inverse against the ORIGINAL full mass, rather than assuming symmetry
    // or silently accepting a poorly resolved factorization.
    for i in 0..n {
        for k in 0..i {
            let value = 0.5*inverse[i*n+k]+0.5*inverse[k*n+i];
            inverse[i*n+k] = value; inverse[k*n+i] = value;
        }
    }
    let mut residual = 0.0_f64;
    for i in 0..n {
        poll(gate)?;
        let mut row_sum = 0.0;
        for k in 0..n {
            let product: f64 = (0..n).map(|j| mass[i*n+j]*inverse[j*n+k]).sum();
            row_sum += finite(product-if i == k { 1.0 } else { 0.0 })?.abs();
        }
        residual = residual.max(finite(row_sum)?);
    }
    if residual > tolerance { return Err(invalid("loaded modal mass inverse exceeds the declared residual allowance")); }
    Ok((inverse, residual))
}

#[cfg(test)]
mod tests;

/// Prepared native stepping and source-pinned, fixed-inertia viscosity feedback.
pub mod viscosity;

/// Exact loaded-inertia coordinates and original fluid memory for shared assemblies.
pub mod modal;
