//! Internal spatial-wrench impulses between two original finite bodies.
//!
//! A force impulse `J` on attachment B and `-J` on A is not angularly internal
//! when the attachments differ: it supplies `(x_B-x_A) cross J`. A connector
//! transmitting couple `C_B` at B must also transmit
//! `C_A = -C_B - (x_B-x_A) cross J` at A. This is wrench transport, not a
//! contact law or an angular-momentum repair after the event.
//!
//! The original full-inertia momentum owner updates each body once. Measured
//! per-body increments, total spin plus orbital momentum and conjugate work
//! must close before publication. The connector may do positive or negative
//! work; internal actuation does not imply passivity or heat generation.

use crate::{DynamicsError, MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};
use fs_math::dd::Dd;

/// Force/couple impulse applied to the second attachment by an internal connector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchImpulse {
    /// Original COM-relative material arms, in each body's own coordinates [m].
    pub arms_body_m: [Vec3; 2],
    /// Force impulse on the second body; the first receives its negative [N s].
    pub impulse_on_second_world_ns: Vec3,
    /// Free couple at the second attachment, in world coordinates [N m s].
    /// The first receives the transported opposite wrench, not just `-couple`.
    pub couple_on_second_world_nms: Vec3,
}

/// Dimensional acceptance limits, without changing any admitted native state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchImpulseControls {
    /// Fixed world point about which total angular change is measured [m].
    pub angular_origin_world_m: Vec3,
    /// Maximum norm of EACH commanded/applied linear-impulse residual and
    /// the complete pair's net linear-momentum change [N s].
    pub linear_tolerance_ns: f64,
    /// Maximum norm of EACH original/applied world central angular residual
    /// (including force/couple summation and frame-conversion losses), and
    /// the complete pair's net orbital-plus-spin angular change [N m s].
    pub angular_tolerance_nms: f64,
    /// Maximum absolute actual kinetic-change minus connector-work defect [J].
    pub work_tolerance_j: f64,
}

/// Complete native event and its measured internal momentum/work accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchImpulseReceipt {
    /// Original native resultant updates; poses are retained exactly.
    pub bodies: [MomentumImpulseReceipt; 2],
    /// Attachment B minus attachment A in the world frame [m]. Computed from
    /// relative COM positions and original arms, without a large common offset.
    pub attachment_separation_world_m: Vec3,
    /// Rounded FREE-couple diagnostics at A and B in world coordinates [N m s].
    /// Work and acceptance retain the unrounded geometric moment transport.
    pub attachment_couples_world_nms: [Vec3; 2],
    /// Actual total final-minus-initial linear momentum [N s].
    pub linear_momentum_change_world_ns: Vec3,
    /// Actual total orbital-plus-spin change about the declared origin [N m s].
    pub angular_momentum_change_world_nms: Vec3,
    /// Sum of actual quadratic kinetic-energy increments, without subtracting
    /// large total energies [J]. Can have either sign.
    pub kinetic_energy_change_j: f64,
    /// Signed work from the supplied and transported impulses at the native
    /// midpoint velocities [J]. This is an instantaneous impulse work identity,
    /// not finite-displacement work or an inferred motor/heat model.
    pub connector_work_j: f64,
    /// Actual kinetic change minus conjugate connector work [J].
    pub work_energy_residual_j: f64,
}

/// Refusal publishes neither a partially changed body nor an accepted account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InternalWrenchImpulseError {
    /// Native state/inertia admission or derived finite arithmetic failed.
    Dynamics(DynamicsError),
    /// The audit origin or a nonnegative dimensional tolerance was invalid.
    InvalidControls,
    /// Maximum per-body/net momentum defects and the signed pair work defect.
    Unresolved {
        /// Maximum per-body or total linear-momentum error [N s].
        linear_residual_ns: f64,
        /// Maximum per-body or total angular-momentum error [N m s].
        angular_residual_nms: f64,
        /// Signed actual kinetic change minus supplied connector work [J].
        work_residual_j: f64,
    },
    /// Cancellation observed before complete publication.
    Cancelled,
}

impl From<DynamicsError> for InternalWrenchImpulseError {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}
impl core::fmt::Display for InternalWrenchImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "internal wrench impulse refused: {self:?}")
    }
}
impl std::error::Error for InternalWrenchImpulseError {}

/// Transfer an internal wrench between distinct native material attachments.
///
/// Both bodies are finite and retain their own full inertia, mass, pose, and
/// world/body momentum coordinates. The attachment arms may differ in position;
/// their free couples then differ by the physical moment-transport term. No
/// contact gap is inferred or repaired, and no joint constraint is imposed.
/// To model two independent external forces without that moment transport, use
/// the original point-impulse owner instead.
///
/// Constant-size double-double reductions retain small transfers against large
/// background momentum. Actual per-body increments are checked as well as net
/// balance: two swallowed opposite impulses cannot masquerade as a transfer.
/// All work is stack-only and bounded, with cancellation between native updates
/// and immediately before return. Inputs are values and never mutate on refusal.
pub fn apply_internal_wrench_impulse(
    states: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    impulse: InternalWrenchImpulse,
    controls: InternalWrenchImpulseControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InternalWrenchImpulseReceipt, InternalWrenchImpulseError> {
    poll(&mut cancelled)?;
    if !controls.angular_origin_world_m.is_finite()
        || [
            controls.linear_tolerance_ns,
            controls.angular_tolerance_nms,
            controls.work_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(InternalWrenchImpulseError::InvalidControls);
    }
    for vector in [
        impulse.impulse_on_second_world_ns,
        impulse.couple_on_second_world_nms,
    ] {
        finite_vector(vector)?;
    }
    // Admission and arm rotation stay with the original rigid-point owner.
    let points = [
        states[0].point_kinematics(properties[0], impulse.arms_body_m[0])?,
        states[1].point_kinematics(properties[1], impulse.arms_body_m[1])?,
    ];
    let separation = add(
        subtract(
            dd(states[1].pose().position_world()),
            dd(states[0].pose().position_world()),
        ),
        subtract(dd(points[1].arm_world), dd(points[0].arm_world)),
    );
    let force = dd(impulse.impulse_on_second_world_ns);
    let couple_b = dd(impulse.couple_on_second_world_nms);
    let couple_a = subtract([Dd::ZERO; 3], add(couple_b, cross(separation, force)));
    let attachment_couples_world_nms = [rounded(couple_a)?, rounded(couple_b)?];
    // Retain the original geometric wrench through work and per-body audit.
    // Rounding is only for the native momentum update: a small supplied couple
    // must not vanish behind a large r cross J before its acceptance check.
    let central = [
        subtract(couple_a, cross(dd(points[0].arm_world), force)),
        add(couple_b, cross(dd(points[1].arm_world), force)),
    ];
    poll(&mut cancelled)?;
    let first = states[0].apply_momentum_impulse(
        properties[0],
        impulse.impulse_on_second_world_ns.scale(-1.0),
        states[0]
            .pose()
            .orientation()
            .rotate_world_to_body(rounded(central[0])?),
    )?;
    poll(&mut cancelled)?;
    let second = states[1].apply_momentum_impulse(
        properties[1],
        impulse.impulse_on_second_world_ns,
        states[1]
            .pose()
            .orientation()
            .rotate_world_to_body(rounded(central[1])?),
    )?;
    let bodies = [first, second];
    let mut linear_change = [Dd::ZERO; 3];
    let mut angular_change = [Dd::ZERO; 3];
    let mut kinetic_change = Dd::ZERO;
    let mut connector_work = Dd::ZERO;
    let mut linear_error = 0.0_f64;
    let mut angular_error = 0.0_f64;
    for i in 0..2 {
        poll(&mut cancelled)?;
        let event = bodies[i];
        let after = event.state_after;
        let velocity_mid = midpoint(
            states[i].center_of_mass_velocity_world(properties[i])?,
            after.center_of_mass_velocity_world(properties[i])?,
        );
        let omega_mid = midpoint(
            properties[i].angular_velocity_body_checked(states[i].angular_momentum_body())?,
            properties[i].angular_velocity_body_checked(after.angular_momentum_body())?,
        );
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .map(|axis| states[i].pose().orientation().rotate_body_to_world(axis));
        let omega_world = rotate_wide(axes, omega_mid);
        // Subtract original endpoint bits in wide precision too. A large kick
        // may swallow the OLD momentum; rounding after-before first would hide
        // that loss just as readily as two swallowed small commanded kicks.
        let actual_linear = subtract(
            dd(after.linear_momentum_world()),
            dd(states[i].linear_momentum_world()),
        );
        let actual_angular = subtract(
            dd(after.angular_momentum_body()),
            dd(states[i].angular_momentum_body()),
        );
        let actual_world_angular = rotate_wide(axes, actual_angular);
        linear_change = add(linear_change, actual_linear);
        let arm = subtract(
            dd(states[i].pose().position_world()),
            dd(controls.angular_origin_world_m),
        );
        angular_change = add(
            angular_change,
            add(cross(arm, actual_linear), actual_world_angular),
        );
        kinetic_change =
            kinetic_change + dot(velocity_mid, actual_linear) + dot(omega_mid, actual_angular);
        connector_work = connector_work
            + dot(velocity_mid, dd(event.linear_impulse_world))
            + dot(omega_world, central[i]);
        linear_error = linear_error
            .max(norm(event.linear_residual_world))
            .max(norm(rounded(subtract(
                actual_linear,
                dd(event.linear_impulse_world),
            ))?));
        angular_error = angular_error
            .max(norm(event.angular_residual_body))
            .max(norm(rounded(subtract(actual_world_angular, central[i]))?));
    }
    let linear_momentum_change_world_ns = rounded(linear_change)?;
    let angular_momentum_change_world_nms = rounded(angular_change)?;
    linear_error = linear_error.max(norm(linear_momentum_change_world_ns));
    angular_error = angular_error.max(norm(angular_momentum_change_world_nms));
    let work_energy_residual_j = finite((kinetic_change - connector_work).to_f64())?;
    if !linear_error.is_finite()
        || !angular_error.is_finite()
        || linear_error > controls.linear_tolerance_ns
        || angular_error > controls.angular_tolerance_nms
        || work_energy_residual_j.abs() > controls.work_tolerance_j
    {
        return Err(InternalWrenchImpulseError::Unresolved {
            linear_residual_ns: linear_error,
            angular_residual_nms: angular_error,
            work_residual_j: work_energy_residual_j,
        });
    }
    let receipt = InternalWrenchImpulseReceipt {
        bodies,
        attachment_separation_world_m: rounded(separation)?,
        attachment_couples_world_nms,
        linear_momentum_change_world_ns,
        angular_momentum_change_world_nms,
        kinetic_energy_change_j: finite(kinetic_change.to_f64())?,
        connector_work_j: finite(connector_work.to_f64())?,
        work_energy_residual_j,
    };
    poll(&mut cancelled)?;
    Ok(receipt)
}

type WideVec = [Dd; 3];
fn dd(v: Vec3) -> WideVec {
    [v.x, v.y, v.z].map(Dd::from_f64)
}
fn add(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[i] + b[i])
}
fn subtract(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[i] - b[i])
}
fn midpoint(a: Vec3, b: Vec3) -> WideVec {
    let half = Dd::from_f64(0.5);
    core::array::from_fn(|i| dd(a)[i] * half + dd(b)[i] * half)
}
fn cross(a: WideVec, b: WideVec) -> WideVec {
    core::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn dot(a: WideVec, b: WideVec) -> Dd {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn rotate_wide(axes_world: [Vec3; 3], body: WideVec) -> WideVec {
    let axes = axes_world.map(dd);
    core::array::from_fn(|i| axes[0][i] * body[0] + axes[1][i] * body[1] + axes[2][i] * body[2])
}
fn rounded(v: WideVec) -> Result<Vec3, InternalWrenchImpulseError> {
    let result = Vec3::new(v[0].to_f64(), v[1].to_f64(), v[2].to_f64());
    finite_vector(result)?;
    Ok(result)
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn finite(value: f64) -> Result<f64, InternalWrenchImpulseError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DynamicsError::NonFinite("internal_wrench.account").into())
    }
}
fn finite_vector(value: Vec3) -> Result<(), InternalWrenchImpulseError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(DynamicsError::NonFinite("internal_wrench.vector").into())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), InternalWrenchImpulseError> {
    if cancelled() {
        Err(InternalWrenchImpulseError::Cancelled)
    } else {
        Ok(())
    }
}
