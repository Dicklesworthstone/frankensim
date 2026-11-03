//! Finite motion of two native bodies driven by one prescribed internal wrench.
//!
//! Two accepted internal half-impulses surround the original native free rigid
//! drift. The opposite wrench is transported at each actual endpoint geometry,
//! so evolving attachments retain both linear and spin-plus-orbital momentum.
//! Full inertia, quaternion evolution and rotational nonlinear solve stay with
//! the original rigid-body owner; no pose or momentum is projected.
//!
//! This symmetric splitting is second order for the smooth held world wrench.
//! Its signed actuator work is the sum of the native conjugate impulse works.
//! Neither that discrete work identity nor a closed momentum budget establishes
//! continuous trajectory accuracy: resolve time and rotation by refinement.
//! This is active prescribed actuation, without a spring, contact or heat law.

use super::internal_impulse::{
    InternalWrenchImpulse, InternalWrenchImpulseControls, InternalWrenchImpulseError,
    InternalWrenchImpulseReceipt, apply_internal_wrench_impulse,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};
use fs_math::dd::Dd;

/// A constant world wrench prescribed on attachment B over one finite step.
///
/// Attachment A receives the transported opposite force/couple at each stage.
/// The command may supply or recover mechanical energy. Material arms are
/// COM-relative and use each body's original inertia/momentum coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrescribedInternalWrench {
    /// Original body-fixed material attachments [m].
    pub arms_body_m: [Vec3; 2],
    /// Force on B, fixed in world coordinates for the step [N].
    pub force_on_second_world_n: Vec3,
    /// Free couple on B, fixed in world coordinates for the step [N m].
    pub couple_on_second_world_nm: Vec3,
}

/// Accepted finite trajectory and actual whole-step physical accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchStep {
    /// Actual duration of the native free drift [s].
    pub duration_s: f64,
    /// Original native states with both poses and momenta advanced.
    pub states_after: [RigidBodyState; 2],
    /// Accepted initial and endpoint impulses, including their distinct
    /// material geometry, transported couples, and signed work.
    pub half_impulses: [InternalWrenchImpulseReceipt; 2],
    /// Fixed world origin used for both native and whole-step angular gates [m].
    pub angular_origin_world_m: Vec3,
    /// Actual whole-step quadratic kinetic-energy increment [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of accepted conjugate half-impulse works [J]. This is a discrete
    /// actuator-work quadrature; it is not an exact continuous power integral.
    pub actuator_work_j: f64,
    /// Signed kinetic-energy change in the original free drift [J]. The native
    /// implicit rotation keeps this at solve tolerance and roundoff. This defect
    /// remains in the total work gate and is never deposited as heat.
    pub native_drift_energy_defect_j: f64,
    /// Actual total linear-momentum change in this internally driven pair [N s].
    pub linear_momentum_residual_ns: Vec3,
    /// Actual total spin-plus-orbital change about the declared origin [N m s].
    pub angular_momentum_residual_nms: Vec3,
    /// Actual kinetic-energy change minus discrete actuator work [J].
    pub work_energy_residual_j: f64,
}

/// A failed trajectory returns no endpoint or partially accepted actuator work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InternalWrenchStepError {
    /// Original native state, full-inertia, duration or nonlinear-solve refusal.
    Dynamics(DynamicsError),
    /// The original impulse owner refused the supplied transfer or its accounts.
    Impulse(InternalWrenchImpulseError),
    /// The duration cannot be halved exactly, or a nonzero held wrench produces
    /// an underflowed/overflowed half-impulse. No source is silently discarded.
    UnrepresentableImpulse,
    /// Measured native free-drift or complete-trajectory budgets did not close.
    Balance {
        /// Some(i) identifies one body's free drift; None is the whole step.
        body: Option<usize>,
        /// Norm of actual linear change [N s].
        linear_residual_ns: f64,
        /// Norm of actual spin-plus-orbital change [N m s].
        angular_residual_nms: f64,
        /// Signed actual kinetic change minus supplied actuator work [J].
        work_residual_j: f64,
    },
    /// Cancellation was observed before complete publication.
    Cancelled,
}
impl From<DynamicsError> for InternalWrenchStepError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<InternalWrenchImpulseError> for InternalWrenchStepError {
    fn from(error: InternalWrenchImpulseError) -> Self {
        match error {
            InternalWrenchImpulseError::Cancelled => Self::Cancelled,
            error => Self::Impulse(error),
        }
    }
}
impl core::fmt::Display for InternalWrenchStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "internal wrench motion refused: {self:?}")
    }
}
impl std::error::Error for InternalWrenchStepError {}

/// Advance two independent finite native bodies through one held internal drive.
///
/// Half-kicks use the existing internal wrench owner; the intervening motion
/// uses RigidBodyIntegrator::step with zero external wrench and zero gravity.
/// The native implicit Euler-top solve and its exact constant-axis special cases
/// retain their original limits/refusals. Both endpoint poses feed wrench
/// transport: freezing the initial opposite couple over the drift is incorrect.
///
/// The original dimensional impulse controls also gate EACH free body's
/// momentum/energy change and the whole trajectory's momentum/work residual.
/// Endpoint differences, world rotations, orbital products and quadratic
/// kinetic increments are reduced in double-double precision. Free-drift
/// errors cannot cancel between bodies to evade their own acceptance gates.
///
/// Work and storage are constant and allocation-free: two native impulse calls
/// and two native free-body steps, without hidden substeps or a replacement
/// inertia/pose solver. Cancellation is checked between these bounded operations
/// and immediately before return. Inputs are values and remain unchanged on
/// every refusal. Repeating accepted calls gives a finite driven trajectory;
/// step refinement is still required for force/rotation phase accuracy.
///
/// No constraint, collision, field-dependent controller, elastic storage,
/// thermal allocation, external load or symplecticity claim is added.
pub fn step_internal_wrench(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    drive: PrescribedInternalWrench,
    duration_s: f64,
    controls: InternalWrenchImpulseControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InternalWrenchStep, InternalWrenchStepError> {
    poll(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 {
        return Err(DynamicsError::InvalidStepDuration.into());
    }
    let half = 0.5 * duration_s;
    if half == 0.0 || 2.0 * half != duration_s {
        return Err(InternalWrenchStepError::UnrepresentableImpulse);
    }
    let impulse = InternalWrenchImpulse {
        arms_body_m: drive.arms_body_m,
        impulse_on_second_world_ns: scaled_impulse(drive.force_on_second_world_n, half)?,
        couple_on_second_world_nms: scaled_impulse(drive.couple_on_second_world_nm, half)?,
    };
    let initial =
        apply_internal_wrench_impulse(states, properties, impulse, controls, &mut cancelled)?;
    let kicked = initial.bodies.map(|event| event.state_after);
    let mut drifted = kicked;
    let mut drift_energy = Dd::ZERO;
    let free = RigidBodyIntegrator::new(Gravity::ZERO);
    for i in 0..2 {
        poll(&mut cancelled)?;
        drifted[i] = free
            .step(kicked[i], properties[i], Wrench::ZERO, duration_s)?
            .state_after;
        poll(&mut cancelled)?;
        let change = state_change(
            kicked[i],
            drifted[i],
            properties[i],
            controls.angular_origin_world_m,
        )?;
        accept(change, Dd::ZERO, controls, Some(i))?;
        drift_energy = drift_energy + change.energy;
    }
    poll(&mut cancelled)?;
    let endpoint =
        apply_internal_wrench_impulse(drifted, properties, impulse, controls, &mut cancelled)?;
    let states_after = endpoint.bodies.map(|event| event.state_after);
    let changes = [
        state_change(
            states[0],
            states_after[0],
            properties[0],
            controls.angular_origin_world_m,
        )?,
        state_change(
            states[1],
            states_after[1],
            properties[1],
            controls.angular_origin_world_m,
        )?,
    ];
    let total = Change {
        linear: add(changes[0].linear, changes[1].linear),
        angular: add(changes[0].angular, changes[1].angular),
        energy: changes[0].energy + changes[1].energy,
    };
    let work = Dd::from_f64(initial.connector_work_j) + Dd::from_f64(endpoint.connector_work_j);
    accept(total, work, controls, None)?;
    let result = InternalWrenchStep {
        duration_s,
        states_after,
        half_impulses: [initial, endpoint],
        angular_origin_world_m: controls.angular_origin_world_m,
        kinetic_energy_change_j: finite(total.energy.to_f64())?,
        actuator_work_j: finite(work.to_f64())?,
        native_drift_energy_defect_j: finite(drift_energy.to_f64())?,
        linear_momentum_residual_ns: rounded(total.linear)?,
        angular_momentum_residual_nms: rounded(total.angular)?,
        work_energy_residual_j: finite((total.energy - work).to_f64())?,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

fn scaled_impulse(force: Vec3, half: f64) -> Result<Vec3, InternalWrenchStepError> {
    if !force.is_finite() {
        return Err(DynamicsError::NonFinite("internal_motion.wrench").into());
    }
    let impulse = force.scale(half);
    if !impulse.is_finite()
        || [force.x, force.y, force.z]
            .into_iter()
            .zip([impulse.x, impulse.y, impulse.z])
            .any(|(source, value)| source != 0.0 && value == 0.0)
    {
        return Err(InternalWrenchStepError::UnrepresentableImpulse);
    }
    Ok(impulse)
}

type WideVec = [Dd; 3];
#[derive(Clone, Copy)]
struct Change {
    linear: WideVec,
    angular: WideVec,
    energy: Dd,
}

fn state_change(
    before: RigidBodyState,
    after: RigidBodyState,
    properties: MassProperties,
    origin: Vec3,
) -> Result<Change, InternalWrenchStepError> {
    let dp = sub(
        dd(after.linear_momentum_world()),
        dd(before.linear_momentum_world()),
    );
    let dh = sub(
        dd(after.angular_momentum_body()),
        dd(before.angular_momentum_body()),
    );
    let energy = dot(
        midpoint(
            before.center_of_mass_velocity_world(properties)?,
            after.center_of_mass_velocity_world(properties)?,
        ),
        dp,
    ) + dot(
        midpoint(
            properties.angular_velocity_body_checked(before.angular_momentum_body())?,
            properties.angular_velocity_body_checked(after.angular_momentum_body())?,
        ),
        dh,
    );
    // Both orbital products use the original stored states and the same origin.
    // World spin uses native pose columns; no rounded total angular momenta are
    // subtracted, and a material-frame change cannot masquerade as an impulse.
    let angular = sub(
        add(
            cross(
                sub(dd(after.pose().position_world()), dd(origin)),
                dd(after.linear_momentum_world()),
            ),
            world_spin(after),
        ),
        add(
            cross(
                sub(dd(before.pose().position_world()), dd(origin)),
                dd(before.linear_momentum_world()),
            ),
            world_spin(before),
        ),
    );
    Ok(Change {
        linear: dp,
        angular,
        energy,
    })
}
fn accept(
    change: Change,
    work: Dd,
    controls: InternalWrenchImpulseControls,
    body: Option<usize>,
) -> Result<(), InternalWrenchStepError> {
    let linear = rounded(change.linear)?;
    let angular = rounded(change.angular)?;
    let energy = finite((change.energy - work).to_f64())?;
    let linear_error = linear.x.hypot(linear.y).hypot(linear.z);
    let angular_error = angular.x.hypot(angular.y).hypot(angular.z);
    if !linear_error.is_finite()
        || !angular_error.is_finite()
        || linear_error > controls.linear_tolerance_ns
        || angular_error > controls.angular_tolerance_nms
        || energy.abs() > controls.work_tolerance_j
    {
        return Err(InternalWrenchStepError::Balance {
            body,
            linear_residual_ns: linear_error,
            angular_residual_nms: angular_error,
            work_residual_j: energy,
        });
    }
    Ok(())
}

fn dd(v: Vec3) -> WideVec {
    [v.x, v.y, v.z].map(Dd::from_f64)
}
fn add(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[i] - b[i])
}
fn midpoint(a: Vec3, b: Vec3) -> WideVec {
    let (a, b, half) = (dd(a), dd(b), Dd::from_f64(0.5));
    core::array::from_fn(|i| a[i] * half + b[i] * half)
}
fn cross(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn dot(a: WideVec, b: WideVec) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn world_spin(state: RigidBodyState) -> WideVec {
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .map(|v| dd(state.pose().orientation().rotate_body_to_world(v)));
    let h = dd(state.angular_momentum_body());
    core::array::from_fn(|i| axes[0][i] * h[0] + axes[1][i] * h[1] + axes[2][i] * h[2])
}
fn finite(value: f64) -> Result<f64, InternalWrenchStepError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DynamicsError::NonFinite("internal_motion.account").into())
    }
}
fn rounded(value: WideVec) -> Result<Vec3, InternalWrenchStepError> {
    Ok(Vec3::new(
        finite(value[0].to_f64())?,
        finite(value[1].to_f64())?,
        finite(value[2].to_f64())?,
    ))
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), InternalWrenchStepError> {
    if cancelled() {
        Err(InternalWrenchStepError::Cancelled)
    } else {
        Ok(())
    }
}
