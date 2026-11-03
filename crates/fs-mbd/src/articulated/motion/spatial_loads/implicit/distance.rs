//! Rigid distance and material-point coincidence loops on the original tree.
//!
//! An objective initial constraint impulse and the original implicit midpoint are
//! solved together for endpoint position closure. A second constraint impulse at
//! that endpoint enforces tangency with the original free-floating mass action.
//! Both impulses retain their actual quadratic work; no position is snapped and
//! no finite-step constraint work is labelled heat or silently discarded.
//!
//! This uses the SHAKE/RATTLE position/hidden-constraint construction described
//! by McLachlan et al., <https://arxiv.org/abs/1207.3367>. The underlying native
//! Lie midpoint is not asserted symplectic, so this composition makes no such
//! claim. Only regular independent constraint rows are admitted. Point closure
//! uses `g = R_A^T (x_B-x_A)`, not three zero-length distance constraints. Its
//! exact differential supplies the A-side geometric couple away from closure.

use super::super::MaterialPoint;
use super::super::bushing::coupled::CoupledSpatialBushing;
use super::super::bushing::{self, BushingLaw, SpatialBushing};
use super::{
    ArticulatedModel, Cell, EffortLaw, FallibleFiniteDifference, FiniteDifferenceError,
    FreeFloatingBaseState, FreeFloatingMotionError, FreeFloatingMotionState, ImplicitSpatialConfig,
    ImplicitSpatialError, MAX_FREE_FLOATING_STEP_LINKS, MidpointEquations, MidpointInitial,
    MomentumMap, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, RefCell, SpatialMotionLoads,
    SpatialMotionStep, Vec3, advance, finish_midpoint, poll, validate_config,
    validate_external_wrenches, validate_forward_dynamics_header, validate_slice,
};
use crate::articulated::motion::reconfigure::constraints::{
    ConstraintError, JointConstraint, admit_constraints, apply_rows,
};
use crate::articulated::motion::sum_products;
use crate::articulated::{
    ArticulatedError, BaseState, FreeFloatingImpulseResponse, JointKind, Kinematics,
    apply_free_floating_impulses, forward_kinematics, free_floating_impulse_delta,
    solve_floating_base_system,
};
use fs_ga::{GaError, Se3, Twist, Wrench};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;

/// Bound on simultaneous independent scalar rows: a point closure consumes three.
pub const MAX_BILATERAL_ROWS: usize = 32;

/// Bound on simultaneous independently reacting distance loops.
pub const MAX_DISTANCE_LOOPS: usize = MAX_BILATERAL_ROWS;

/// A positive rigid separation between two distinct original link material points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceLoop {
    /// First endpoint, unchanged in its link frame.
    pub a: MaterialPoint,
    /// Second endpoint on a different link.
    pub b: MaterialPoint,
    /// Required finite positive distance in metres.
    pub length_m: f64,
}

/// Coincidence of two material points, with all relative rotations unrestricted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointCoincidenceLoop {
    /// First endpoint and the material frame in which the three rows are expressed.
    pub a: MaterialPoint,
    /// Second endpoint on a different original link.
    pub b: MaterialPoint,
}

/// One rigid loop, in caller order, participating in the same simultaneous solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BilateralLoop {
    /// One positive scalar separation constraint.
    Distance(DistanceLoop),
    /// Three objective material-frame components of a spherical/ball closure.
    PointCoincidence(PointCoincidenceLoop),
}

impl BilateralLoop {
    fn rows(self) -> usize {
        match self {
            Self::Distance(_) => 1,
            Self::PointCoincidence(_) => 3,
        }
    }

    fn endpoints(self) -> [MaterialPoint; 2] {
        match self {
            Self::Distance(c) => [c.a, c.b],
            Self::PointCoincidence(c) => [c.a, c.b],
        }
    }
}

/// Interpretation of one flattened residual and reaction-impulse coordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BilateralRow {
    /// Original caller loop index.
    pub constraint: usize,
    /// `None` for distance; `Some(0..=2)` for the A-material-frame x/y/z axis.
    pub component: Option<usize>,
}

/// Both objective coordinates and actual physical endpoint observations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointClosureResidual {
    /// Original caller loop index.
    pub constraint: usize,
    /// `R_A^T (x_B-x_A)`, in metres.
    pub separation_a_m: Vec3,
    /// `R_A^T (v_B-v_A-omega_A cross (x_B-x_A))`, in m/s.
    pub separation_rate_a_mps: Vec3,
    /// Actual `x_B-x_A` in world coordinates, computed without common translation.
    pub separation_world_m: Vec3,
    /// Actual `v_B-v_A` in world coordinates, including common angular motion.
    pub relative_velocity_world_mps: Vec3,
}

/// Nonlinear work, physical constraint tolerances and regularity admission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BilateralLoopConfig {
    /// Existing scaled-rate solve policy. Its evaluation cap includes the
    /// coupled solve's stencils and its actual initial-impulse verification.
    pub nonlinear: ImplicitSpatialConfig,
    /// Maximum initial/final distance error or point separation norm, in metres.
    pub position_tolerance_m: f64,
    /// Maximum distance rate, point-coordinate rate norm AND actual point
    /// relative-velocity norm at the initial/final configuration, in m/s.
    pub velocity_tolerance_mps: f64,
    /// Minimum mobility/free-body ratio and squared Cholesky pivot of the
    /// unit-diagonal mobility, in [4096*EPSILON,1).
    /// This is a floating regularity admission, not a certified condition bound.
    pub minimum_scaled_pivot: f64,
}

impl Default for BilateralLoopConfig {
    fn default() -> Self {
        Self {
            nonlinear: ImplicitSpatialConfig::default(),
            position_tolerance_m: 1.0e-9,
            velocity_tolerance_mps: 1.0e-9,
            minimum_scaled_pivot: 1.0e-10,
        }
    }
}

/// Controls for the distance-only specialization of the common loop engine.
pub type DistanceLoopConfig = BilateralLoopConfig;

/// One position/velocity-admitted constrained step, including both impulse events.
#[derive(Debug, Clone, PartialEq)]
pub struct BilateralLoopStep {
    /// Canonical final constrained state; use this for the next step.
    pub state_after: FreeFloatingMotionState,
    /// Native smooth part between the initial and final constraint impulses.
    /// Its state is the intermediate endpoint BEFORE final tangent correction.
    pub midpoint: SpatialMotionStep,
    /// Actual native initial constraint impulse response, with quadratic work.
    pub initial_reaction: FreeFloatingImpulseResponse,
    /// Actual native final tangent impulse response, with quadratic work.
    pub final_reaction: FreeFloatingImpulseResponse,
    /// Flattened row interpretation for all residual and scalar impulse vectors.
    pub rows: Vec<BilateralRow>,
    /// Initial impulse per row, N s. Positive distance impulse pushes apart;
    /// point components give the impulse on B in A's initial material frame.
    pub initial_impulse_ns: Vec<f64>,
    /// Endpoint impulse per row, N s, using A's ENDPOINT material frame.
    pub final_impulse_ns: Vec<f64>,
    /// Actual initial signed length or A-frame point-component residuals, metres.
    pub position_residual_before_m: Vec<f64>,
    /// Actual final row residuals, metres; positions are never projected.
    pub position_residual_after_m: Vec<f64>,
    /// Actual initial objective row rates, m/s.
    pub velocity_residual_before_mps: Vec<f64>,
    /// Actual final objective row rates, m/s.
    pub velocity_residual_after_mps: Vec<f64>,
    /// Initial objective and physical observations, in point-loop caller order.
    pub point_residual_before: Vec<PointClosureResidual>,
    /// Corresponding endpoint observations; both physical and objective rates pass.
    pub point_residual_after: Vec<PointClosureResidual>,
    /// Initial A-side `-(x_B-x_A) cross impulse_on_B`, in world N m s,
    /// in point-loop caller order. This is part of the constraint impulse.
    pub initial_geometric_couple_impulse_world: Vec<Vec3>,
    /// Corresponding endpoint geometric couples in the endpoint configuration.
    pub final_geometric_couple_impulse_world: Vec<Vec3>,
    /// Initial geometric-couple work at actual midpoint event angular velocities,
    /// J; already included in `constraint_impulse_work`, not additional work.
    pub initial_geometric_couple_work: f64,
    /// Corresponding endpoint geometric-couple work, J.
    pub final_geometric_couple_work: f64,
    /// Minimum admitted unit-diagonal mobility pivot at the initial configuration.
    pub initial_scaled_pivot: f64,
    /// The corresponding endpoint admission; no constraint may silently disappear.
    pub final_scaled_pivot: f64,
    /// Coupled position/midpoint solve report, with the actual impulse rechecked.
    pub nonlinear: NewtonReport,
    /// All nonlinear stage evaluations including directional-difference probes.
    pub residual_evaluations: usize,
    /// Sum of the actual initial/final impulse work, J. Ideal continuous
    /// constraints do zero work, but this finite-step defect is not assumed zero.
    pub constraint_impulse_work: f64,
    /// Actual change of kinetic+finite-spring storage plus damping loss minus
    /// actuator/external work, J. Includes constraint-work and integration defects.
    pub mechanical_work_residual: f64,
    /// Actual net linear momentum change minus the external smooth impulse, N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual orbital-plus-spin change minus the external smooth impulse, N m s.
    pub angular_momentum_residual_world: Vec3,
    /// Actual native summed reaction linear impulse, including floating defects.
    pub reaction_linear_impulse_world: Vec3,
    /// Actual native summed reaction angular impulse about the world origin.
    /// Objective constraints give zero in exact arithmetic, also off manifold.
    pub reaction_angular_impulse_world: Vec3,
    /// Linear momentum defect after also subtracting the actual reaction impulse.
    pub linear_momentum_balance_residual_world: Vec3,
    /// Angular momentum defect after also subtracting the actual reaction impulse.
    pub angular_momentum_balance_residual_world: Vec3,
}

/// One simultaneous body-loop/constant-joint step.
///
/// Body-loop vectors retain their original metre/N s coordinates. These extra
/// vectors follow the supplied joint-row order and its declared Angle/Length
/// units. Both native impulse receipts and their total work in `bilateral`
/// already include the joint reactions and every body/joint cross term.
#[derive(Debug, Clone, PartialEq)]
pub struct MixedBilateralLoopStep {
    pub bilateral: BilateralLoopStep,
    /// Initial A^T multiplier coordinates: N m s for Angle, N s for Length.
    pub initial_joint_multipliers: Vec<f64>,
    /// Corresponding final multipliers at the actual endpoint configuration.
    pub final_joint_multipliers: Vec<f64>,
    /// Actual summed initial A^T lambda in compact native joint order.
    pub initial_joint_impulse: Vec<f64>,
    /// Actual summed final A^T lambda in compact native joint order.
    pub final_joint_impulse: Vec<f64>,
    /// Actual A q - b at the two independently admitted configurations.
    pub joint_position_residuals_before: Vec<f64>,
    pub joint_position_residuals_after: Vec<f64>,
    /// Actual A qdot at the two independently admitted physical states.
    pub joint_velocity_residuals_before: Vec<f64>,
    pub joint_velocity_residuals_after: Vec<f64>,
    /// Joint impulse paired with the actual mean pre/post event rates, J.
    /// This is INCLUDED in bilateral.constraint_impulse_work, not extra work.
    pub initial_joint_reaction_work_j: f64,
    /// Corresponding endpoint joint work, also included in the total.
    pub final_joint_reaction_work_j: f64,
}

/// Result of the distance-only specialization; point-specific vectors are empty.
pub type DistanceLoopStep = BilateralLoopStep;

/// Atomic refusal: no state, arbitrary dependent reaction or clipped pose is returned.
#[derive(Debug, Clone, PartialEq)]
pub enum BilateralLoopError {
    /// Malformed loop, controls, count or nonfinite arithmetic.
    Invalid(&'static str),
    /// Original force/dynamics/nonlinear/budget/cancellation refusal.
    Motion(ImplicitSpatialError),
    /// Input or endpoint is outside the admitted constraint position tolerance.
    Position {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// Flattened row index; a point norm reports its first row.
        constraint: usize,
        /// Actual signed distance/component error or point separation norm, m.
        residual_m: f64,
    },
    /// Initial or final velocity is not tangent within the requested tolerance.
    Velocity {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// Flattened row index; a point norm reports its first row.
        constraint: usize,
        /// Actual objective rate component or point-coordinate rate norm, m/s.
        residual_mps: f64,
    },
    /// Objective tangency alone does not hide a physical finite-gap spin velocity.
    PointVelocity {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// Original caller loop index.
        constraint: usize,
        /// Actual world-space relative velocity, m/s.
        relative_velocity_world_mps: Vec3,
    },
    /// Dependent, singular, insufficiently regular or nonreciprocal mobility.
    /// This lane refuses instead of assigning an arbitrary redundant reaction.
    Mobility {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// First refused flattened row/pivot in the supplied order.
        constraint: usize,
        /// Concrete floating admission that failed.
        reason: &'static str,
    },
}

/// Refusals shared with the distance-only specialization.
pub type DistanceLoopError = BilateralLoopError;

impl core::fmt::Display for DistanceLoopError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated bilateral loop: {self:?}")
    }
}
impl std::error::Error for DistanceLoopError {}
impl From<ImplicitSpatialError> for DistanceLoopError {
    fn from(value: ImplicitSpatialError) -> Self {
        Self::Motion(value)
    }
}
impl From<FreeFloatingMotionError> for DistanceLoopError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<ArticulatedError> for DistanceLoopError {
    fn from(value: ArticulatedError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<ConstraintError> for DistanceLoopError {
    fn from(value: ConstraintError) -> Self {
        match value {
            ConstraintError::Cancelled => FreeFloatingMotionError::Cancelled.into(),
            other => Self::Motion(other.into()),
        }
    }
}
impl From<GaError> for DistanceLoopError {
    fn from(value: GaError) -> Self {
        Self::Motion(value.into())
    }
}

/// Advance a closed tree with regular internal rigid distance loops.
///
/// Initial positions AND velocities must satisfy the declared constraints. The
/// first central impulses and native midpoint rates are simultaneous unknowns;
/// the endpoint tangent solve uses the full free-floating mobility, including
/// base recoil, all joints, off-centre arms and complete inertia tensors. Both
/// impulse applications use the original owner. No mass is clamped or replaced
/// by a penalty spring. Forces, finite springs and dampers retain their separate
/// storage/work account from the unconstrained lane.
///
/// At most 32 independent distance loops and 4096 links are admitted. Cached
/// initial inverse-mass columns cost O(loops*links), with O(loops*joints) storage;
/// each coupled residual adds one original FK pass to the native midpoint work.
/// The endpoint solves only a loops-by-loops SPD mobility. No dense generalized
/// mass/Jacobian matrix or nested nonlinear integration is introduced.
///
/// Returned endpoint residuals are floating observations, not a state enclosure.
/// Intermediate drift states need not lie on the manifold. There is no collision
/// detection, general nonlinear joint/gear loops, nonholonomic rolling, automatic assembly,
/// redundancy resolution, branch certificate or adjoint. Numerical constraint
/// work remains visible; no general symplecticity or energy-exactness is claimed.
pub fn step_free_floating_distance_loops(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[DistanceLoop],
    duration: f64,
    config: DistanceLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<DistanceLoopStep, DistanceLoopError> {
    step_free_floating_distance_loops_with_bushings(
        model,
        state,
        loads,
        loops,
        &[],
        duration,
        config,
        cancelled,
    )
}

/// Advance the same rigid loops together with finite spatial bushings.
///
/// Both finite bushing strains/rates and rigid loop reactions use the original
/// tree at the simultaneous stage. All bushing storage, loss and internal work
/// is included in `midpoint` and the complete mechanical-work residual. A rigid
/// distance does not turn the bushing into a penalty constraint or erase its
/// tangential/rotational mechanics. An empty bushing slice is exactly the plain
/// distance-loop operation. All of that operation's limits and refusals apply.
pub fn step_free_floating_distance_loops_with_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[DistanceLoop],
    bushings: &[SpatialBushing],
    duration: f64,
    config: DistanceLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<DistanceLoopStep, DistanceLoopError> {
    if loops.len() > MAX_DISTANCE_LOOPS {
        return Err(DistanceLoopError::Invalid("distance loop count exceeds 32"));
    }
    let constraints: Vec<_> = loops.iter().copied().map(BilateralLoop::Distance).collect();
    step_free_floating_bilateral_loops_with_bushings(
        model,
        state,
        loads,
        &constraints,
        bushings,
        duration,
        config,
        cancelled,
    )
}

/// Advance a simultaneous mixture of rigid distance and spherical/ball loops.
///
/// Point closure constrains translation only. Its objective rows are
/// `g = R_A^T (x_B-x_A)`, with derivative
/// `gdot = R_A^T (v_B-v_A-omega_A cross (x_B-x_A))`. For an impulse `lambda`
/// conjugate to `gdot`, B receives `F = R_A lambda`, A receives `-F`, and A also
/// receives `-(x_B-x_A) cross F`. This geometric couple is the exact derivative
/// of the material frame, not an angular-momentum repair after the solve. The
/// same complete wrench enters mobility, both native impulses, and work.
///
/// The initial and endpoint point separation norm, objective rate norm, AND
/// actual world relative-velocity norm must pass the physical tolerances. The
/// last check matters for a small admitted gap spinning quickly. There is no
/// coordinate snapping, frozen base, zero-length distance normalization, large
/// penalty stiffness, or sequential independent-loop projection.
///
/// Up to 32 independent scalar rows are admitted; each point closure consumes
/// three. The original bounded midpoint/Krylov solve, full native impulse map,
/// regularity tests, explicit work defects and cancellation semantics apply.
/// Point reactions have three components in A's material frame at EACH event;
/// initial and final components generally use different world directions.
pub fn step_free_floating_bilateral_loops(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    duration: f64,
    config: BilateralLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralLoopStep, BilateralLoopError> {
    step_free_floating_bilateral_loops_with_bushings(
        model,
        state,
        loads,
        loops,
        &[],
        duration,
        config,
        cancelled,
    )
}

/// The same simultaneous rigid loops composed with finite spatial bushings.
///
/// All finite forces are reevaluated at the simultaneous native stage, and their
/// storage/loss remains in the complete work account. An empty bushing slice is
/// exactly [`step_free_floating_bilateral_loops`].
pub fn step_free_floating_bilateral_loops_with_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    bushings: &[SpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralLoopStep, BilateralLoopError> {
    step_bilateral_evaluated(
        model,
        state,
        loads,
        loops,
        &[],
        bushings,
        duration,
        config,
        None,
        cancelled,
    )
    .map(|step| step.bilateral)
}

/// Advance rigid bilateral loops, finite bushings and state-dependent internal
/// joint effort in the same simultaneous native solve.
///
/// The callback receives elapsed time `duration / 2`, the complete reconstructed
/// trial midpoint, and an output slice reset to `loads.held.generalized_force`.
/// It may add to or replace those efforts. The trial includes the current
/// initial constraint impulse; it need not yet satisfy closure. Joint forces,
/// rigid reactions and bushing strains/rates therefore participate in one
/// residual and its numerical Jacobian, including the independent verification
/// with the actual summed native initial impulse.
///
/// Native ABA supplies the carrier reaction. The accepted midpoint effort
/// supplies actuator work; neither joint effort nor objective loop reaction is
/// an external momentum source. Both native constraint impulses retain their
/// actual work, and endpoint position, hidden-velocity and mobility admissions
/// are unchanged. The final `state_after`, including the endpoint reaction,
/// remains the canonical state for the next step.
///
/// The effort units, repeatability and trial-history requirements of
/// [`super::step_free_floating_spatial_implicit_with_effort`] apply. In
/// particular, callbacks must not advance accepted physical/controller history
/// during probes. Their side effects cannot be rolled back. General effort
/// laws need not be passive; actuator storage/loss remains caller-owned.
/// One O(joints) output buffer is reused throughout the bounded solve.
///
/// The existing midpoint stability, nonlinear geometry and refinement limits
/// remain. This does not add coupled electrical/controller state, eliminate
/// finite-step constraint work or identify a unique nonlinear branch.
///
/// # Errors
/// All original loop/motion/regularity/budget refusals, invalid or over-limit
/// trial effort, cancellation, or a callback error. No accepted native state
/// is published on refusal; no effort, constraint or position is clipped.
pub fn step_free_floating_bilateral_loops_with_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    bushings: &[SpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralLoopStep, BilateralLoopError> {
    step_bilateral_evaluated(
        model,
        state,
        loads,
        loops,
        &[],
        bushings,
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
    .map(|step| step.bilateral)
}

/// Advance rigid loops with the original full passive stiffness and damping laws.
///
/// Every entry of the prepared six-coordinate matrices acts at the same trial
/// midpoint as the rigid closures. Translation/rotation and cross-axis forces
/// use the existing reciprocal material-frame wrench lift; they cannot inject
/// external momentum. The original loop impulse work, bushing storage and
/// damping loss remain separate in the returned mechanics account.
///
/// An empty slice gives the original bilateral-loop motion. Diagonal matrices
/// have the same physical meaning as the diagonal-bushing entry point.
pub fn step_free_floating_bilateral_loops_with_coupled_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralLoopStep, BilateralLoopError> {
    step_bilateral_evaluated(
        model,
        state,
        loads,
        loops,
        &[],
        bushings,
        duration,
        config,
        None,
        cancelled,
    )
    .map(|step| step.bilateral)
}

/// Simultaneously solve rigid closures, full passive bushings and joint feedback.
///
/// The callback has the same repeatable-trial contract as
/// [`step_free_floating_bilateral_loops_with_effort`]. Its actual midpoint
/// effort is the actuator-work source; passive bushing reactions and the two
/// rigid-loop impulse events remain in their existing storage/loss/work owners.
/// Neither the matrix law nor the rigid constraints waive native effort,
/// speed, geometry, regularity, cancellation or nonlinear-budget admission.
#[allow(clippy::too_many_arguments)]
pub fn step_free_floating_bilateral_loops_with_coupled_bushings_and_effort(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    bushings: &[CoupledSpatialBushing],
    duration: f64,
    config: BilateralLoopConfig,
    mut effort: impl FnMut(
        f64,
        &FreeFloatingMotionState,
        &mut [f64],
    ) -> Result<(), FreeFloatingMotionError>,
    cancelled: impl FnMut() -> bool,
) -> Result<BilateralLoopStep, BilateralLoopError> {
    step_bilateral_evaluated(
        model,
        state,
        loads,
        loops,
        &[],
        bushings,
        duration,
        config,
        Some(&mut effort),
        cancelled,
    )
    .map(|step| step.bilateral)
}

/// Constant ideal gears/locks combined with distance and spherical closures.
///
/// All rows participate in one native free-floating impulse mobility, including
/// their cross responses. The initial impulses and midpoint are solved
/// simultaneously; the endpoint tangent impulse enforces all rows together.
/// Each joint row retains its own units and absolute tolerances. Trial probes
/// may leave the manifold; accepted positions and rates are independently
/// admitted. No gear is enforced by a sequential velocity projection.
///
/// The combined system admits at most 32 scalar rows; each spherical closure
/// consumes three. Mobility refusal indices list the flattened body-loop rows
/// first, then the joint rows in caller order. Body-specific metre/N s receipt
/// vectors remain separate from the original-unit joint vectors.
///
/// Both impulse events retain their actual work, including cross terms. This
/// RATTLE-style composition is distinct from the continuous affine-reaction
/// midpoint lane: intermediate kicked rates need not be tangent, and finite
/// impulse work is not presumed zero. Activation, changing ratios, backlash,
/// nonholonomic relations and arbitrary redundant reactions are not modeled.
pub mod joint {
    use super::*;

    /// Advance body loops and ideal joint rows with diagonal passive bushings.
    /// An empty bushing slice gives the plain mixed-constraint operation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_bushings(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        loads: SpatialMotionLoads<'_>,
        loops: &[BilateralLoop],
        rows: &[JointConstraint],
        bushings: &[SpatialBushing],
        duration: f64,
        config: BilateralLoopConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MixedBilateralLoopStep, BilateralLoopError> {
        step_bilateral_evaluated(
            model, state, loads, loops, rows, bushings, duration, config, None, cancelled,
        )
    }

    /// The same solve with smooth, repeatable trial actuator feedback.
    /// The parent bilateral callback contract and physical effort limits apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_bushings_and_effort(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        loads: SpatialMotionLoads<'_>,
        loops: &[BilateralLoop],
        rows: &[JointConstraint],
        bushings: &[SpatialBushing],
        duration: f64,
        config: BilateralLoopConfig,
        mut effort: impl FnMut(
            f64,
            &FreeFloatingMotionState,
            &mut [f64],
        ) -> Result<(), FreeFloatingMotionError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MixedBilateralLoopStep, BilateralLoopError> {
        step_bilateral_evaluated(
            model,
            state,
            loads,
            loops,
            rows,
            bushings,
            duration,
            config,
            Some(&mut effort),
            cancelled,
        )
    }

    /// The same simultaneous rows with complete passive six-coordinate matrices.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_coupled_bushings(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        loads: SpatialMotionLoads<'_>,
        loops: &[BilateralLoop],
        rows: &[JointConstraint],
        bushings: &[CoupledSpatialBushing],
        duration: f64,
        config: BilateralLoopConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MixedBilateralLoopStep, BilateralLoopError> {
        step_bilateral_evaluated(
            model, state, loads, loops, rows, bushings, duration, config, None, cancelled,
        )
    }

    /// Jointly solve all rigid rows, full passive matrices and actuator feedback.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_coupled_bushings_and_effort(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        loads: SpatialMotionLoads<'_>,
        loops: &[BilateralLoop],
        rows: &[JointConstraint],
        bushings: &[CoupledSpatialBushing],
        duration: f64,
        config: BilateralLoopConfig,
        mut effort: impl FnMut(
            f64,
            &FreeFloatingMotionState,
            &mut [f64],
        ) -> Result<(), FreeFloatingMotionError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MixedBilateralLoopStep, BilateralLoopError> {
        step_bilateral_evaluated(
            model,
            state,
            loads,
            loops,
            rows,
            bushings,
            duration,
            config,
            Some(&mut effort),
            cancelled,
        )
    }
}

fn step_bilateral_evaluated<B: BushingLaw>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[BilateralLoop],
    joint_rows: &[JointConstraint],
    bushings: &[B],
    duration: f64,
    config: BilateralLoopConfig,
    effort: Option<&mut EffortLaw<'_>>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixedBilateralLoopStep, BilateralLoopError> {
    poll(&mut cancelled)?;
    validate_config(config.nonlinear)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration.into());
    }
    if loops.is_empty() && joint_rows.is_empty()
        || loops.len() > MAX_BILATERAL_ROWS
        || joint_rows.len() > MAX_BILATERAL_ROWS
        || loops.iter().map(|c| c.rows()).sum::<usize>() + joint_rows.len() > MAX_BILATERAL_ROWS
    {
        return Err(DistanceLoopError::Invalid(
            "independent bilateral row count must be in 1..=32",
        ));
    }
    if [
        config.position_tolerance_m,
        config.velocity_tolerance_mps,
        config.minimum_scaled_pivot,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || config.minimum_scaled_pivot >= 1.0
        || config.minimum_scaled_pivot < 4096. * f64::EPSILON
    {
        return Err(DistanceLoopError::Invalid(
            "finite positive constraint tolerances and resolved mobility pivot required",
        ));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    model.free_floating_complexity()?;
    validate_forward_dynamics_header(
        model,
        state.base.twist_body,
        &state.joint_position,
        &state.joint_velocity,
        loads.held.generalized_force,
        loads.held.gravity_world,
        loads.held.external_body_wrench,
    )?;
    validate_external_wrenches(loads.held.external_body_wrench)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    admit_loops(model, loops)?;
    let joint_before = admit_constraints(
        joint_rows,
        &state.joint_position,
        &state.joint_velocity,
        &mut cancelled,
    )?;
    let initial_map = MomentumMap::new(model, state, &mut cancelled)?;
    let initial_coordinates = initial_map.coordinates(state.base)?;
    let initial_load =
        bushing::evaluate(&initial_map, state.base, loads, bushings, &mut cancelled)?.load;
    let before = geometry(&initial_map.shape, state.base, loops)?;
    check_position(&before, config, "initial")?;
    check_velocity(&before, config, "initial")?;
    let initial_mobility = Mobility::new(
        model,
        state,
        &initial_map,
        &before,
        joint_rows,
        config,
        "initial",
        true,
        &mut cancelled,
    )?;
    let n = model.dof_count();
    let loop_count = before.gaps.len();
    let count = loop_count + joint_rows.len();
    let mut scales = vec![config.nonlinear.angular_velocity_scale; n + 3];
    for (link, part) in model.links().iter().enumerate() {
        if part.joint.kind() == JointKind::Prismatic
            && let Some(dof) = model.dof_index(link)
        {
            scales[dof] = config.nonlinear.linear_velocity_scale;
        }
    }
    // A coefficient carries row-units / joint-units. Scaling by the native
    // coordinate rates preserves these units and row-rescaling covariance.
    let mut row_scales = vec![config.nonlinear.linear_velocity_scale; loop_count];
    for row in joint_rows {
        let mut scale = 0.0_f64;
        for (&a, &v) in row.coefficients.iter().zip(&scales) {
            let term = a.abs() * v;
            if !term.is_finite() {
                return Err(DistanceLoopError::Invalid(
                    "unrepresentable joint row scaling",
                ));
            }
            scale = scale.max(term);
        }
        row_scales.push(scale);
    }
    let impulse_scales: Vec<_> = initial_mobility
        .diagonal
        .iter()
        .zip(&row_scales)
        .map(|(d, scale)| scale / (d * d))
        .collect();
    let gap_scales: Vec<_> = row_scales.iter().map(|v| duration * v).collect();
    let mut tolerances = vec![config.position_tolerance_m; loop_count];
    tolerances.extend(joint_rows.iter().map(|row| row.position_tolerance));
    let target = tolerances.iter().zip(&gap_scales).fold(
        config.nonlinear.residual_tolerance,
        |target, (tolerance, scale)| target.min(tolerance / scale),
    );
    if gap_scales.iter().any(|s| !s.is_finite() || *s <= 0.0)
        || !target.is_finite()
        || target <= 0.0
        || impulse_scales.iter().any(|s| !s.is_finite() || *s <= 0.0)
    {
        return Err(DistanceLoopError::Invalid(
            "unrepresentable constraint scaling",
        ));
    }
    let gravity_only_external = loads.world_points.is_empty()
        && loads
            .held
            .external_body_wrench
            .iter()
            .all(|w| w.force == Vec3::default() && w.torque == Vec3::default());
    let mut guess = state.joint_velocity.clone();
    let w = state.base.twist_body.angular;
    guess.extend([w.x, w.y, w.z]);
    for (x, s) in guess.iter_mut().zip(&scales) {
        *x /= s;
    }
    guess.resize(n + 3 + count, 0.0);
    let evaluations = Cell::new(0usize);
    let cancelled = RefCell::new(cancelled);
    let effort_values = if effort.is_some() {
        loads.held.generalized_force.to_vec()
    } else {
        Vec::new()
    };
    let effort = RefCell::new((effort, effort_values));
    let charge = || -> Result<(), DistanceLoopError> {
        let used = evaluations.get();
        if used >= config.nonlinear.maximum_evaluations {
            return Err(ImplicitSpatialError::EvaluationBudget {
                used,
                maximum: config.nonlinear.maximum_evaluations,
            }
            .into());
        }
        evaluations.set(used + 1);
        Ok(())
    };
    let residual = |x: &[f64]| -> Result<Vec<f64>, DistanceLoopError> {
        let mut cancel = cancelled.borrow_mut();
        poll(&mut *cancel)?;
        charge()?;
        let mut kicked = state.clone();
        let mut coordinates = initial_coordinates;
        for (j, column) in initial_mobility.columns.iter().enumerate() {
            let magnitude = x[n + 3 + j] * impulse_scales[j];
            for (v, dv) in kicked.joint_velocity.iter_mut().zip(&column.joint) {
                *v += magnitude * dv;
            }
            kicked.base.twist_body = kicked.base.twist_body.plus(column.base.scale(magnitude));
            // Objective rows have zero resultants in exact arithmetic. Retain
            // the actual rounded wrench resultants instead of assuming that.
            coordinates.linear_world =
                coordinates.linear_world + column.linear_world.scale(magnitude);
            coordinates.angular_com_world =
                coordinates.angular_com_world + column.angular_com_world.scale(magnitude);
        }
        let equations = MidpointEquations {
            model,
            state: &kicked,
            loads,
            duration,
            initial: coordinates,
            mass: initial_map.mass(),
            scales: &scales,
            gravity_only_external,
            bushings,
        };
        let midpoint = {
            let mut effort = effort.borrow_mut();
            let (law, values) = &mut *effort;
            if let Some(law) = law.as_deref_mut() {
                equations.evaluate_with_effort(
                    &x[..n + 3],
                    Some((law, values.as_mut_slice())),
                    &mut *cancel,
                )?
            } else {
                equations.evaluate(&x[..n + 3], &mut *cancel)?
            }
        };
        let q = advance(
            &state.joint_position,
            &midpoint.state.joint_velocity,
            duration,
        );
        poll(&mut *cancel)?;
        let shape = forward_kinematics(
            model,
            BaseState::stationary(Se3::identity()),
            &q,
            &vec![0.; n],
        )?;
        let end = geometry(&shape, state.base, loops)?;
        let mut residual = midpoint.residual;
        let mut gaps = end.gaps;
        gaps.extend(joint_gaps(joint_rows, &q, &mut *cancel)?);
        residual.extend(gaps.iter().zip(&gap_scales).map(|(g, s)| g / s));
        validate_slice(&residual, "distance_loop.equations")?;
        poll(&mut *cancel)?;
        Ok(residual)
    };
    let problem =
        FallibleFiniteDifference::new(n + 3 + count, config.nonlinear.difference_step, residual)
            .map_err(difference_error)?;
    let solver_config = NewtonKrylovConfig {
        absolute_tolerance: target,
        relative_tolerance: 0.,
        linear_restart: config.nonlinear.krylov_restart.min(n + 3 + count),
        max_linear_cycles: 8,
        forcing_minimum: 1e-3,
        forcing_maximum: 0.1,
        ..Default::default()
    };
    let created = NewtonKrylovState::new(&problem, guess, solver_config);
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = created.map_err(ImplicitSpatialError::Solver)?;
    let run = solver.run_controlled(&problem, config.nonlinear.maximum_iterations, &mut || {
        problem.has_failed() || (*cancelled.borrow_mut())()
    });
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    match run.progress {
        fs_solver::nonlinear::SolverRunProgress::Complete => {}
        fs_solver::nonlinear::SolverRunProgress::Paused => {
            return Err(FreeFloatingMotionError::Cancelled.into());
        }
        other => return Err(ImplicitSpatialError::SolverProgress(other).into()),
    }
    let mut report = run.report;
    if !report.converged {
        return Err(ImplicitSpatialError::Unresolved(report).into());
    }
    let initial_impulse: Vec<_> = solver.x[n + 3..]
        .iter()
        .zip(&impulse_scales)
        .map(|(x, s)| x * s)
        .collect();
    let mut cancel = cancelled.borrow_mut();
    poll(&mut *cancel)?;
    let (initial_reaction, initial_joint_impulse) = apply(
        model,
        state,
        &before,
        joint_rows,
        &initial_impulse,
        &mut *cancel,
    )?;
    let kicked = FreeFloatingMotionState {
        base: initial_reaction.base_after,
        joint_position: state.joint_position.clone(),
        joint_velocity: initial_reaction.joint_velocity_after.clone(),
    };
    // Re-evaluate with the ACTUAL native summed impulse, not cached unit-column
    // superposition or assumed zero floating momentum defect.
    charge()?;
    let kicked_map = MomentumMap::new(model, &kicked, &mut *cancel)?;
    let kicked_coordinates = kicked_map.coordinates(kicked.base)?;
    let equations = MidpointEquations {
        model,
        state: &kicked,
        loads,
        duration,
        initial: kicked_coordinates,
        mass: kicked_map.mass(),
        scales: &scales,
        gravity_only_external,
        bushings,
    };
    let actual_midpoint = {
        let mut effort = effort.borrow_mut();
        let (law, values) = &mut *effort;
        if let Some(law) = law.as_deref_mut() {
            equations.evaluate_with_effort(
                &solver.x[..n + 3],
                Some((law, values.as_mut_slice())),
                &mut *cancel,
            )?
        } else {
            equations.evaluate(&solver.x[..n + 3], &mut *cancel)?
        }
    };
    let final_q = advance(
        &kicked.joint_position,
        &actual_midpoint.state.joint_velocity,
        duration,
    );
    let shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &final_q,
        &vec![0.; n],
    )?;
    let final_geometry = geometry(&shape, state.base, loops)?;
    let mut actual_residual = actual_midpoint.residual.clone();
    let mut final_gaps = final_geometry.gaps;
    final_gaps.extend(joint_gaps(joint_rows, &final_q, &mut *cancel)?);
    actual_residual.extend(final_gaps.iter().zip(&gap_scales).map(|(g, s)| g / s));
    report.residual_norm = fs_solver::norm2(&actual_residual);
    if !report.residual_norm.is_finite() || report.residual_norm > target {
        report.converged = false;
        return Err(ImplicitSpatialError::Unresolved(report).into());
    }
    let (initial_geometric_couple_impulse_world, initial_geometric_couple_work) =
        geometric_work(&before, &kicked_map.shape, kicked.base, &initial_impulse)?;
    let midpoint = finish_midpoint(
        model,
        &kicked,
        loads,
        duration,
        MidpointInitial {
            map: kicked_map,
            coordinates: kicked_coordinates,
            load: initial_load,
            momentum: initial_reaction.momentum_after,
        },
        actual_midpoint,
        bushings,
        &mut *cancel,
    )?;
    let endpoint = &midpoint.motion.state_after;
    let final_map = MomentumMap::new(model, endpoint, &mut *cancel)?;
    let end = geometry(&final_map.shape, endpoint.base, loops)?;
    check_position(&end, config, "endpoint")?;
    let final_mobility = Mobility::new(
        model,
        endpoint,
        &final_map,
        &end,
        joint_rows,
        config,
        "endpoint",
        false,
        &mut *cancel,
    )?;
    let mut end_rates = end.rates.clone();
    end_rates.extend(apply_rows(
        joint_rows,
        &endpoint.joint_velocity,
        &mut *cancel,
    )?);
    let final_impulse = final_mobility.cancel_rates(&end_rates)?;
    poll(&mut *cancel)?;
    let (final_reaction, final_joint_impulse) = apply(
        model,
        endpoint,
        &end,
        joint_rows,
        &final_impulse,
        &mut *cancel,
    )?;
    let state_after = FreeFloatingMotionState {
        base: final_reaction.base_after,
        joint_position: endpoint.joint_position.clone(),
        joint_velocity: final_reaction.joint_velocity_after.clone(),
    };
    let final_shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    let actual = geometry(&final_shape, state_after.base, loops)?;
    check_position(&actual, config, "endpoint")?;
    check_velocity(&actual, config, "endpoint")?;
    let joint_after = admit_constraints(
        joint_rows,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut *cancel,
    )?;
    let initial_joint_reaction_work_j = joint_work(
        &initial_joint_impulse,
        &state.joint_velocity,
        &kicked.joint_velocity,
    )?;
    let final_joint_reaction_work_j = joint_work(
        &final_joint_impulse,
        &endpoint.joint_velocity,
        &state_after.joint_velocity,
    )?;
    let (final_geometric_couple_impulse_world, final_geometric_couple_work) =
        geometric_work(&end, &final_shape, state_after.base, &final_impulse)?;
    let initial_momentum = initial_reaction.momentum_before;
    let final_momentum = final_reaction.momentum_after;
    let constraint_impulse_work = initial_reaction.impulse_work + final_reaction.impulse_work;
    let mechanical_work_residual = final_momentum.kinetic_energy - initial_momentum.kinetic_energy
        + midpoint.stored_energy_after
        - midpoint.stored_energy_before
        + midpoint.dissipated_energy
        - midpoint.motion.external_work
        - midpoint.motion.actuator_work;
    let linear_momentum_residual_world = final_momentum.linear_momentum_world
        - initial_momentum.linear_momentum_world
        - midpoint.motion.external_linear_impulse_world;
    let angular_momentum_residual_world = final_momentum.angular_momentum_world
        - initial_momentum.angular_momentum_world
        - midpoint.motion.external_angular_impulse_world;
    let reaction_linear_impulse_world = initial_reaction.external_linear_impulse_world
        + final_reaction.external_linear_impulse_world;
    let reaction_angular_impulse_world = initial_reaction.external_angular_impulse_world
        + final_reaction.external_angular_impulse_world;
    let linear_momentum_balance_residual_world =
        linear_momentum_residual_world - reaction_linear_impulse_world;
    let angular_momentum_balance_residual_world =
        angular_momentum_residual_world - reaction_angular_impulse_world;
    validate_slice(
        &[
            constraint_impulse_work,
            mechanical_work_residual,
            initial_geometric_couple_work,
            final_geometric_couple_work,
        ],
        "distance_loop.work",
    )?;
    super::validate_vec3(
        linear_momentum_residual_world,
        "distance_loop.linear_momentum",
    )?;
    super::validate_vec3(
        angular_momentum_residual_world,
        "distance_loop.angular_momentum",
    )?;
    for value in [
        reaction_linear_impulse_world,
        reaction_angular_impulse_world,
        linear_momentum_balance_residual_world,
        angular_momentum_balance_residual_world,
    ] {
        super::validate_vec3(value, "bilateral_loop.reaction_momentum")?;
    }
    poll(&mut *cancel)?;
    Ok(MixedBilateralLoopStep {
        initial_joint_multipliers: initial_impulse[loop_count..].to_vec(),
        final_joint_multipliers: final_impulse[loop_count..].to_vec(),
        initial_joint_impulse,
        final_joint_impulse,
        joint_position_residuals_before: joint_before.position,
        joint_position_residuals_after: joint_after.position,
        joint_velocity_residuals_before: joint_before.velocity,
        joint_velocity_residuals_after: joint_after.velocity,
        initial_joint_reaction_work_j,
        final_joint_reaction_work_j,
        bilateral: DistanceLoopStep {
            state_after,
            midpoint,
            initial_reaction,
            final_reaction,
            rows: before.rows,
            initial_impulse_ns: initial_impulse[..loop_count].to_vec(),
            final_impulse_ns: final_impulse[..loop_count].to_vec(),
            position_residual_before_m: before.gaps,
            position_residual_after_m: actual.gaps,
            velocity_residual_before_mps: before.rates,
            velocity_residual_after_mps: actual.rates,
            point_residual_before: before.points.iter().map(|p| p.residual).collect(),
            point_residual_after: actual.points.iter().map(|p| p.residual).collect(),
            initial_geometric_couple_impulse_world,
            final_geometric_couple_impulse_world,
            initial_geometric_couple_work,
            final_geometric_couple_work,
            initial_scaled_pivot: initial_mobility.minimum_pivot,
            final_scaled_pivot: final_mobility.minimum_pivot,
            nonlinear: report,
            residual_evaluations: evaluations.get(),
            constraint_impulse_work,
            mechanical_work_residual,
            linear_momentum_residual_world,
            angular_momentum_residual_world,
            reaction_linear_impulse_world,
            reaction_angular_impulse_world,
            linear_momentum_balance_residual_world,
            angular_momentum_balance_residual_world,
        },
    })
}

fn joint_gaps(
    rows: &[JointConstraint],
    position: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, BilateralLoopError> {
    let mut gaps = apply_rows(rows, position, cancelled)?;
    for (gap, row) in gaps.iter_mut().zip(rows) {
        *gap -= row.position;
    }
    validate_slice(&gaps, "bilateral_loop.joint_gaps")?;
    Ok(gaps)
}

fn joint_work(impulse: &[f64], before: &[f64], after: &[f64]) -> Result<f64, BilateralLoopError> {
    let work = ((sum_products(impulse.iter().copied(), before.iter().copied())
        + sum_products(impulse.iter().copied(), after.iter().copied()))
        * Dd::from_f64(0.5))
    .to_f64();
    validate_slice(&[work], "bilateral_loop.joint_work")?;
    Ok(work)
}

struct Geometry {
    gaps: Vec<f64>,
    rates: Vec<f64>,
    terms: Vec<[(usize, Wrench); 2]>,
    rows: Vec<BilateralRow>,
    points: Vec<PointGeometry>,
}

struct PointGeometry {
    residual: PointClosureResidual,
    link_a: usize,
    row_start: usize,
    angular_a_world: Vec3,
    couple_world: [Vec3; 3],
}

fn admit_loops(model: &ArticulatedModel, loops: &[BilateralLoop]) -> Result<(), DistanceLoopError> {
    for c in loops {
        let [a, b] = c.endpoints();
        if a.link >= model.link_count() || b.link >= model.link_count() || a.link == b.link {
            return Err(DistanceLoopError::Invalid("bilateral material endpoints"));
        }
        if let BilateralLoop::Distance(d) = c
            && (!d.length_m.is_finite() || d.length_m <= 0.)
        {
            return Err(DistanceLoopError::Invalid(
                "finite positive distance required",
            ));
        }
        for p in [a, b] {
            super::validate_vec3(p.point_body, "bilateral_loop.material_point")?;
        }
    }
    Ok(())
}

fn geometry(
    shape: &Kinematics,
    base: FreeFloatingBaseState,
    loops: &[BilateralLoop],
) -> Result<Geometry, DistanceLoopError> {
    let count: usize = loops.iter().map(|c| c.rows()).sum();
    let mut gaps = Vec::with_capacity(count);
    let mut rates = Vec::with_capacity(count);
    let mut terms = Vec::with_capacity(count);
    let mut rows = Vec::with_capacity(count);
    let mut points = Vec::new();
    let rotation = base.world_from_base.rotation();
    for (constraint, c) in loops.iter().enumerate() {
        let [pa, pb] = c.endpoints();
        let a = shape.world_from_link[pa.link];
        let b = shape.world_from_link[pb.link];
        let delta = b.transform_point(pb.point_body)? - a.transform_point(pa.point_body)?;
        let va = shape.body_twist[pa.link];
        let vb = shape.body_twist[pb.link];
        // The shape has a stationary identity base. Common translation is never
        // subtracted, and common angular motion cancels in objective row rates.
        let relative = b
            .rotation()
            .rotate(vb.linear + vb.angular.cross(pb.point_body))?
            - a.rotation()
                .rotate(va.linear + va.angular.cross(pa.point_body))?;
        match c {
            BilateralLoop::Distance(c) => {
                let length = delta.norm();
                if !length.is_finite() || length <= 0. {
                    return Err(DistanceLoopError::Invalid(
                        "coincident/nonfinite distance endpoints",
                    ));
                }
                let normal = delta.scale(1. / length);
                let fa = a.rotation().inverse().rotate(normal.scale(-1.))?;
                let fb = b.rotation().inverse().rotate(normal)?;
                gaps.push(length - c.length_m);
                rates.push(normal.dot(relative));
                terms.push([
                    (pa.link, Wrench::new(pa.point_body.cross(fa), fa)),
                    (pb.link, Wrench::new(pb.point_body.cross(fb), fb)),
                ]);
                rows.push(BilateralRow {
                    constraint,
                    component: None,
                });
            }
            BilateralLoop::PointCoincidence(_) => {
                let gap = a.rotation().inverse().rotate(delta)?;
                let rate = a.rotation().inverse().rotate(relative)? - va.angular.cross(gap);
                let separation_world_m = rotation.rotate(delta)?;
                let relative_velocity_world_mps =
                    rotation.rotate(relative + base.twist_body.angular.cross(delta))?;
                let angular_a_world =
                    rotation.rotate(base.twist_body.angular + a.rotation().rotate(va.angular)?)?;
                for (value, field) in [
                    (gap, "point_loop.gap"),
                    (rate, "point_loop.objective_rate"),
                    (separation_world_m, "point_loop.world_gap"),
                    (relative_velocity_world_mps, "point_loop.world_rate"),
                    (angular_a_world, "point_loop.angular_rate"),
                ] {
                    super::validate_vec3(value, field)?;
                }
                let row_start = gaps.len();
                gaps.extend([gap.x, gap.y, gap.z]);
                rates.extend([rate.x, rate.y, rate.z]);
                let mut couple_world = [Vec3::default(); 3];
                for (component, direction) in [
                    Vec3::new(1., 0., 0.),
                    Vec3::new(0., 1., 0.),
                    Vec3::new(0., 0., 1.),
                ]
                .into_iter()
                .enumerate()
                {
                    let fa = direction.scale(-1.);
                    let fb = b
                        .rotation()
                        .inverse()
                        .rotate(a.rotation().rotate(direction)?)?;
                    // lambda dot gdot = F dot (vB-vA) - (d cross F) dot omegaA.
                    // Include this exact geometric couple in the same native
                    // wrench used by the mobility, application and work audit.
                    let couple_a = gap.cross(direction).scale(-1.);
                    couple_world[component] = rotation.rotate(a.rotation().rotate(couple_a)?)?;
                    terms.push([
                        (pa.link, Wrench::new(pa.point_body.cross(fa) + couple_a, fa)),
                        (pb.link, Wrench::new(pb.point_body.cross(fb), fb)),
                    ]);
                    rows.push(BilateralRow {
                        constraint,
                        component: Some(component),
                    });
                }
                points.push(PointGeometry {
                    residual: PointClosureResidual {
                        constraint,
                        separation_a_m: gap,
                        separation_rate_a_mps: rate,
                        separation_world_m,
                        relative_velocity_world_mps,
                    },
                    link_a: pa.link,
                    row_start,
                    angular_a_world,
                    couple_world,
                });
            }
        }
    }
    validate_slice(&gaps, "bilateral_loop.gaps")?;
    validate_slice(&rates, "bilateral_loop.rates")?;
    Ok(Geometry {
        gaps,
        rates,
        terms,
        rows,
        points,
    })
}

fn geometric_work(
    before: &Geometry,
    after_shape: &Kinematics,
    after_base: FreeFloatingBaseState,
    magnitude: &[f64],
) -> Result<(Vec<Vec3>, f64), DistanceLoopError> {
    let mut couples = Vec::with_capacity(before.points.len());
    let mut work = 0.;
    for p in &before.points {
        let mut couple = Vec3::default();
        for (axis, direction) in p.couple_world.iter().enumerate() {
            couple = couple + direction.scale(magnitude[p.row_start + axis]);
        }
        let angular_after = after_base.world_from_base.rotation().rotate(
            after_base.twist_body.angular
                + after_shape.world_from_link[p.link_a]
                    .rotation()
                    .rotate(after_shape.body_twist[p.link_a].angular)?,
        )?;
        work += couple.dot((p.angular_a_world + angular_after).scale(0.5));
        super::validate_vec3(couple, "point_loop.geometric_couple")?;
        couples.push(couple);
    }
    validate_slice(&[work], "point_loop.geometric_work")?;
    Ok((couples, work))
}

struct ImpulseColumn {
    joint: Vec<f64>,
    base: Twist,
    linear_world: Vec3,
    angular_com_world: Vec3,
}

struct Mobility {
    diagonal: Vec<f64>,
    factor: Cholesky,
    minimum_pivot: f64,
    columns: Vec<ImpulseColumn>,
}
impl Mobility {
    fn new(
        model: &ArticulatedModel,
        state: &FreeFloatingMotionState,
        map: &MomentumMap,
        geometry: &Geometry,
        joint_rows: &[JointConstraint],
        config: DistanceLoopConfig,
        stage: &'static str,
        retain_columns: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DistanceLoopError> {
        let loop_count = geometry.terms.len();
        let k = loop_count + joint_rows.len();
        let mut matrix = vec![0.; k * k];
        let mut columns = Vec::with_capacity(if retain_columns { k } else { 0 });
        let mut joints = vec![0.; model.dof_count()];
        let mut body = vec![Wrench::default(); model.link_count()];
        for column in 0..k {
            poll(cancelled)?;
            let terms = geometry.terms.get(column);
            if let Some(terms) = terms {
                for (link, wrench) in terms {
                    body[*link] = *wrench;
                }
            } else {
                joints.copy_from_slice(&joint_rows[column - loop_count].coefficients);
            }
            let response =
                free_floating_impulse_delta(model, &state.joint_position, &joints, &body)?;
            for (row, row_terms) in geometry.terms.iter().enumerate() {
                matrix[row * k + column] = row_terms
                    .iter()
                    .map(|(link, w)| w.pairing(response.body[*link]))
                    .sum();
            }
            for (row, value) in apply_rows(joint_rows, &response.joint, cancelled)?
                .into_iter()
                .enumerate()
            {
                matrix[(loop_count + row) * k + column] = value;
            }
            if retain_columns {
                let mut force_base = Vec3::default();
                let mut torque_com_base = Vec3::default();
                if let Some(terms) = terms {
                    for (link, wrench) in terms {
                        let pose = map.shape.world_from_link[*link];
                        let force = pose.rotation().rotate(wrench.force)?;
                        force_base = force_base + force;
                        torque_com_base = torque_com_base
                            + pose.rotation().rotate(wrench.torque)?
                            + (pose.translation() - map.center_base).cross(force);
                    }
                }
                // Joint covectors are internal and have no external resultant.
                let rotation = state.base.world_from_base.rotation();
                let linear_world = rotation.rotate(force_base)?;
                let angular_com_world = rotation.rotate(torque_com_base)?;
                super::validate_vec3(linear_world, "bilateral_loop.column_linear")?;
                super::validate_vec3(angular_com_world, "bilateral_loop.column_angular")?;
                columns.push(ImpulseColumn {
                    joint: response.joint,
                    base: response.base,
                    linear_world,
                    angular_com_world,
                });
            }
            if let Some(terms) = terms {
                for (link, _) in terms {
                    body[*link] = Wrench::default();
                }
            } else {
                joints.fill(0.0);
            }
        }
        validate_slice(&matrix, "distance_loop.mobility")?;
        let mut diagonal = Vec::with_capacity(k);
        for i in 0..k {
            let d = matrix[i * k + i];
            if d <= 0. || !d.is_finite() {
                return Err(DistanceLoopError::Mobility {
                    stage,
                    constraint: i,
                    reason: "constraint mobility is not finite and positive",
                });
            }
            if let Some(terms) = geometry.terms.get(i) {
                let mut free_mobility = 0.;
                for &(link, wrench) in terms {
                    let response =
                        solve_floating_base_system(&model.links()[link].inertia.matrix()?, wrench)?;
                    free_mobility += wrench.pairing(response);
                }
                if !free_mobility.is_finite()
                    || free_mobility <= 0.
                    || d / free_mobility < config.minimum_scaled_pivot
                {
                    return Err(DistanceLoopError::Mobility {
                        stage,
                        constraint: i,
                        reason: "constraint mobility vanishes relative to free endpoint bodies",
                    });
                }
            }
            diagonal.push(d.sqrt());
        }
        for i in 0..k {
            matrix[i * k + i] = 1.;
            for j in 0..i {
                let a = (matrix[i * k + j] / diagonal[i]) / diagonal[j];
                let b = (matrix[j * k + i] / diagonal[j]) / diagonal[i];
                if !a.is_finite() || !b.is_finite() || (a - b).abs() > 1e-10 {
                    return Err(DistanceLoopError::Mobility {
                        stage,
                        constraint: i,
                        reason: "constraint mobility reciprocity unresolved",
                    });
                }
                matrix[i * k + j] = 0.5 * (a + b);
                matrix[j * k + i] = matrix[i * k + j];
            }
        }
        poll(cancelled)?;
        let factor = cholesky(&matrix, k).map_err(|error| DistanceLoopError::Mobility {
            stage,
            constraint: match error {
                fs_la::factor::FactorError::NotSpd { index }
                | fs_la::factor::FactorError::Singular { index } => index,
            },
            reason: "dependent or indefinite constraint mobility",
        })?;
        let mut minimum_pivot = 1.0_f64;
        for i in 0..k {
            let pivot = factor.l(i, i) * factor.l(i, i);
            if !pivot.is_finite() || pivot < config.minimum_scaled_pivot {
                return Err(DistanceLoopError::Mobility {
                    stage,
                    constraint: i,
                    reason: "constraint mobility is insufficiently regular",
                });
            }
            minimum_pivot = minimum_pivot.min(pivot);
        }
        poll(cancelled)?;
        Ok(Self {
            diagonal,
            factor,
            minimum_pivot,
            columns,
        })
    }
    fn cancel_rates(&self, rates: &[f64]) -> Result<Vec<f64>, DistanceLoopError> {
        let mut impulse: Vec<_> = rates
            .iter()
            .zip(&self.diagonal)
            .map(|(v, d)| -v / d)
            .collect();
        self.factor.solve(&mut impulse);
        for (j, d) in impulse.iter_mut().zip(&self.diagonal) {
            *j /= d;
        }
        validate_slice(&impulse, "distance_loop.endpoint_impulse")?;
        Ok(impulse)
    }
}

fn apply(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    geometry: &Geometry,
    joint_rows: &[JointConstraint],
    magnitude: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(FreeFloatingImpulseResponse, Vec<f64>), DistanceLoopError> {
    let mut body = vec![Wrench::default(); model.link_count()];
    for (terms, scale) in geometry.terms.iter().zip(magnitude) {
        for (link, w) in terms {
            body[*link] = Wrench::new(
                body[*link].torque + w.torque.scale(*scale),
                body[*link].force + w.force.scale(*scale),
            );
        }
    }
    let mut joint_impulse = vec![0.0; model.dof_count()];
    let multipliers = &magnitude[geometry.terms.len()..];
    for (joint, impulse) in joint_impulse.iter_mut().enumerate() {
        if joint % 64 == 0 {
            poll(cancelled)?;
        }
        *impulse = sum_products(
            joint_rows.iter().map(|row| row.coefficients[joint]),
            multipliers.iter().copied(),
        )
        .to_f64();
    }
    validate_slice(&joint_impulse, "bilateral_loop.joint_impulse")?;
    let response = apply_free_floating_impulses(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &joint_impulse,
        &body,
    )?;
    Ok((response, joint_impulse))
}
fn check_position(
    geometry: &Geometry,
    config: DistanceLoopConfig,
    stage: &'static str,
) -> Result<(), DistanceLoopError> {
    for (constraint, &residual_m) in geometry.gaps.iter().enumerate() {
        if residual_m.abs() > config.position_tolerance_m {
            return Err(DistanceLoopError::Position {
                stage,
                constraint,
                residual_m,
            });
        }
    }
    for point in &geometry.points {
        let residual_m = point
            .residual
            .separation_a_m
            .norm()
            .max(point.residual.separation_world_m.norm());
        if residual_m > config.position_tolerance_m {
            return Err(DistanceLoopError::Position {
                stage,
                constraint: point.row_start,
                residual_m,
            });
        }
    }
    Ok(())
}
fn check_velocity(
    geometry: &Geometry,
    config: DistanceLoopConfig,
    stage: &'static str,
) -> Result<(), DistanceLoopError> {
    for (constraint, &residual_mps) in geometry.rates.iter().enumerate() {
        if residual_mps.abs() > config.velocity_tolerance_mps {
            return Err(DistanceLoopError::Velocity {
                stage,
                constraint,
                residual_mps,
            });
        }
    }
    for point in &geometry.points {
        let residual_mps = point.residual.separation_rate_a_mps.norm();
        if residual_mps > config.velocity_tolerance_mps {
            return Err(DistanceLoopError::Velocity {
                stage,
                constraint: point.row_start,
                residual_mps,
            });
        }
        let relative_velocity_world_mps = point.residual.relative_velocity_world_mps;
        if relative_velocity_world_mps.norm() > config.velocity_tolerance_mps {
            return Err(DistanceLoopError::PointVelocity {
                stage,
                constraint: point.residual.constraint,
                relative_velocity_world_mps,
            });
        }
    }
    Ok(())
}
fn difference_error(error: FiniteDifferenceError<DistanceLoopError>) -> DistanceLoopError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => ImplicitSpatialError::Difference(reason).into(),
    }
}

#[cfg(test)]
mod point_row_tests {
    use super::*;
    use fs_ga::{So3, So3Tangent};

    #[test]
    fn objective_point_rows_match_pose_differentials_and_internal_spatial_work() {
        let a = Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.4, -0.3, 0.2))).unwrap(),
            Vec3::new(0.2, 0.3, -0.1),
        )
        .unwrap();
        let b = Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(-0.2, 0.6, -0.1))).unwrap(),
            Vec3::new(0.6, -0.2, 0.4),
        )
        .unwrap();
        let constraint = BilateralLoop::PointCoincidence(PointCoincidenceLoop {
            a: MaterialPoint {
                link: 0,
                point_body: Vec3::new(0.2, -0.1, 0.3),
            },
            b: MaterialPoint {
                link: 1,
                point_body: Vec3::new(-0.3, 0.1, 0.2),
            },
        });
        let base = FreeFloatingBaseState::stationary(Se3::identity());
        let mut maximum_difference = 0.0_f64;
        let mut geometric_signal = 0.0_f64;
        for coordinate in 0..12 {
            let mut velocities = vec![Twist::zero(); 2];
            let mut axis = [0.; 6];
            axis[coordinate % 6] = 1.;
            velocities[coordinate / 6] = Twist::new(
                Vec3::new(axis[0], axis[1], axis[2]),
                Vec3::new(axis[3], axis[4], axis[5]),
            );
            let shape = Kinematics {
                world_from_link: vec![a, b],
                body_twist: velocities.clone(),
            };
            let rows = geometry(&shape, base, &[constraint]).unwrap();
            let perturbed = |h| {
                let poses = [a, b]
                    .into_iter()
                    .zip(&velocities)
                    .map(|(p, v)| p.compose(Se3::exp(v.scale(h)).unwrap()).unwrap())
                    .collect();
                geometry(
                    &Kinematics {
                        world_from_link: poses,
                        body_twist: velocities.clone(),
                    },
                    base,
                    &[constraint],
                )
                .unwrap()
            };
            let epsilon = 2e-6;
            let plus = perturbed(epsilon);
            let minus = perturbed(-epsilon);
            for row in 0..3 {
                let derivative = (plus.gaps[row] - minus.gaps[row]) / (2. * epsilon);
                let work: f64 = rows.terms[row]
                    .iter()
                    .map(|(link, w)| w.pairing(velocities[*link]))
                    .sum();
                maximum_difference = maximum_difference.max((derivative - work).abs());
                assert!(
                    (derivative - work).abs() < 2e-10,
                    "row={row}, coordinate={coordinate}: {derivative} != {work}"
                );
                assert!((rows.rates[row] - work).abs() < 5e-16);
                let mut force = Vec3::default();
                let mut torque = Vec3::default();
                for (link, w) in &rows.terms[row] {
                    let world = w.transform_by(&shape.world_from_link[*link]);
                    force = force + world.force;
                    torque = torque + world.torque;
                }
                assert!(force.norm() < 5e-16);
                assert!(torque.norm() < 8e-16);
                geometric_signal = geometric_signal.max(rows.points[0].couple_world[row].norm());
            }
            let moving_base = FreeFloatingBaseState::new(
                Se3::from_parts(a.rotation(), Vec3::new(1e3, -3e2, 2e2)).unwrap(),
                Twist::new(Vec3::new(100., -200., 300.), Vec3::new(1e12, -2e12, 3e12)),
            );
            let transformed = geometry(&shape, moving_base, &[constraint]).unwrap();
            assert_eq!(rows.gaps, transformed.gaps);
            assert_eq!(rows.rates, transformed.rates);
            assert_eq!(rows.terms, transformed.terms);
        }
        assert!(geometric_signal > 0.1);
        eprintln!(
            "objective point-row finite-difference maximum defect={maximum_difference:e}; required geometric-couple signal={geometric_signal:e}"
        );
    }
}
