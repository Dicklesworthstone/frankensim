//! Rigid payload capture and release on articulated trees.
//!
//! The root API describes free motion; [`supported`] retains fixed or driven
//! support reactions and work during prescribed-base capture and release.
//!
//! An admitted attachment preserves the target link frame and every joint.
//! Full rotated centroidal inertia and COM orbital inertia enter the original
//! articulated-body mass solve; existing constant gear rows may react in that
//! same event. The supplied poses define the attachment geometry. Contact
//! discovery, finite-duration grasping, mass flux and a thermal destination for
//! the removed relative kinetic energy remain caller-owned physical models.

/// Simultaneous capture of multiple rigid payloads by one articulated tree.
#[path = "inventory/batch.rs"]
mod batch;
pub use batch::{PayloadBatchCapture, PayloadCaptureInput, capture_payloads};

/// Conservative payload release with optional internal separation impulses.
#[path = "inventory/release.rs"]
pub mod release;

/// Rigid payload capture with fixed or driven prescribed support momentum/work.
#[path = "inventory/supported.rs"]
pub mod supported;

/// Payload inventory events retaining distance, spherical and ideal joint rows.
#[path = "inventory/bilateral.rs"]
pub mod bilateral;

use super::reconfigure::constraints::{
    self, ConstraintControls, ConstraintError, JointConstraint, admit_constraints,
};
use super::reconfigure::{
    self, ReconfigurationControls, ReconfigurationError, audit_momentum, state_kinematics,
};
use super::spatial_loads::implicit::ImplicitSpatialError;
use super::spatial_loads::implicit::distance::BilateralLoopError;
use super::{
    FreeFloatingMotionError, FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS, SpatialSum,
};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, Kinematics, Link,
    SpatialInertia, inverse_dynamics, momentum_from_kinematics, validate_slice, validate_twist,
};
use core::fmt;
use fs_ga::{GaError, Mat3, Se3, Twist, Vec3, Wrench};
use fs_math::dd::Dd;

/// A free rigid payload, with inertia and twist expressed in its own body frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PayloadState {
    /// Physically admitted full inertia about an arbitrary payload-frame origin.
    pub inertia: SpatialInertia,
    /// Actual payload pose at the event; capture never snaps this pose.
    pub world_from_body: Se3,
    /// Actual angular and origin-linear velocities in payload coordinates.
    pub twist_body: Twist,
}

/// Momentum, energy and independent-constraint acceptance for an inventory event.
#[derive(Debug, Clone, Copy)]
pub struct PayloadControls<'a> {
    /// Relative componentwise momentum tolerance, finite and in `[0, 1e-8]`.
    pub relative_momentum_tolerance: f64,
    /// Positive absolute base/world angular momentum tolerance, N m s.
    pub angular_momentum_tolerance_nms: f64,
    /// Positive absolute base/world linear momentum tolerance, N s.
    pub linear_momentum_tolerance_ns: f64,
    /// Positive per-joint momentum allowances in native compact joint order.
    pub joint_momentum_tolerances: &'a [f64],
    /// Nonnegative absolute common-frame energy-balance residual allowance, J.
    pub energy_tolerance_j: f64,
    /// Nonnegative maximum capture loss or absolute supplied release work, J.
    pub maximum_energy_exchange_j: f64,
    /// Nonnegative absolute actual ideal-constraint midpoint-work allowance, J.
    pub maximum_constraint_work_j: f64,
    /// Normalized row-mobility reciprocal-condition floor in `[4096*EPSILON, 1]`.
    pub minimum_reciprocal_condition: f64,
}

/// Inventory identity created by capture; callers cannot forge its contents.
///
/// The exact loaded target inertia and its original kinematic signature bind
/// this handle. Further captures on that link must be released in reverse order.
#[derive(Debug, Clone, PartialEq)]
pub struct PayloadAttachment {
    target_link: usize,
    retained_link: Link,
    loaded_inertia: SpatialInertia,
    link_from_payload: Se3,
    payload_inertia: SpatialInertia,
}

impl PayloadAttachment {
    /// Target link index in the unchanged articulated topology.
    #[must_use]
    pub const fn target_link(&self) -> usize {
        self.target_link
    }

    /// Fixed transform from payload coordinates into the target link frame.
    #[must_use]
    pub const fn link_from_payload(&self) -> Se3 {
        self.link_from_payload
    }

    /// Original payload-frame inertia, retained without reconstruction.
    #[must_use]
    pub const fn payload_inertia(&self) -> SpatialInertia {
        self.payload_inertia
    }

    /// Exact target-link inertia to restore when this attachment is released.
    #[must_use]
    pub fn retained_inertia(&self) -> SpatialInertia {
        self.retained_link.inertia()
    }

    /// Combined inertia to which the attachment is bound.
    #[must_use]
    pub const fn loaded_inertia(&self) -> SpatialInertia {
        self.loaded_inertia
    }
}

/// Independent physical and native generalized-momentum accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct PayloadBalance {
    /// Actual complete-system world momentum and lab-frame energy before.
    pub momentum_before: ArticulatedMomentum,
    /// Actual complete-system world momentum and lab-frame energy after.
    pub momentum_after: ArticulatedMomentum,
    /// Original native `[base angular, base linear, joint]` projected momenta.
    pub generalized_momentum_before: Vec<f64>,
    /// Independently reconstructed momenta of the actual returned bodies.
    pub generalized_momentum_after: Vec<f64>,
    /// Native constraint covector in `[base angular, base linear, joint]` order.
    /// Gear-only events use `[0; A^T lambda]`; bilateral events also retain the
    /// original body-row projection, including its floating base resultants.
    pub generalized_constraint_impulse: Vec<f64>,
    /// Gear-row impulse multipliers in caller row order.
    pub multipliers: Vec<f64>,
    /// Actual constraint impulse paired with mean pre/post rates in the common
    /// translating accounting frame, J; never inferred heat.
    pub reaction_work_j: f64,
    /// Largest normalized physical or generalized momentum residual, at most one.
    pub maximum_normalized_momentum_residual: f64,
    /// Constant world velocity of the common inertial energy-accounting frame.
    pub accounting_boost_world: Vec3,
    /// Complete-system kinetic energy in that translating frame before, J.
    pub kinetic_before_j: f64,
    /// Complete-system kinetic energy in the same translating frame after, J.
    pub kinetic_after_j: f64,
    /// Actual signed common-frame energy defect; no numerical defect is heat.
    pub energy_residual_j: f64,
}

/// Complete accepted capture; all inputs and budgets remain caller-owned.
#[derive(Debug, Clone, PartialEq)]
pub struct PayloadCapture {
    /// Identical topology and joints with the target's loaded full inertia.
    pub model_after: ArticulatedModel,
    /// Unchanged coordinates and pose with simultaneous capture/reaction rates.
    pub state_after: FreeFloatingMotionState,
    /// Exact inventory identity required for subsequent release.
    pub attachment: PayloadAttachment,
    /// Actual rigidly attached payload pose and velocity in the returned state.
    pub payload_after: PayloadState,
    /// Nonnegative relative-motion quadratic loss; no thermal model is implied.
    pub removed_relative_energy_j: f64,
    /// Independent momentum and common-translating-frame energy audits.
    pub balance: PayloadBalance,
}

/// Atomic refusal without publishing a partial model, state or attachment.
#[derive(Debug, Clone, PartialEq)]
pub enum PayloadError {
    /// Invalid inventory, dimensions, controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Native geometry, full inertia, state, limit or articulated solve refusal.
    Dynamics(ArticulatedError),
    /// Native ideal-joint-row admission or bounded mobility solve refusal.
    Constraint(ConstraintError),
    /// Original mixed material-loop/joint impulse admission refused.
    Bilateral(BilateralLoopError),
    /// Shared momentum or state admission refusal.
    Reconfiguration(ReconfigurationError),
    /// Cancellation observed before publishing the operation.
    Cancelled,
    /// Capture loss or absolute supplied release work exceeds its allowance.
    EnergyExchange {
        /// Actual signed work or nonnegative removed energy, J.
        energy_j: f64,
        /// Declared absolute allowance, J.
        maximum_j: f64,
    },
    /// Actual common-frame energy defect exceeds its independent tolerance.
    EnergyResidual {
        /// Signed energy defect, J.
        residual_j: f64,
        /// Declared absolute tolerance, J.
        maximum_j: f64,
    },
    /// Actual finite ideal-reaction midpoint work exceeds its separate budget.
    ConstraintWork {
        /// Signed reaction work, J.
        work_j: f64,
        /// Declared absolute allowance, J.
        maximum_j: f64,
    },
}

impl fmt::Display for PayloadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "articulated payload inventory refused: {self:?}")
    }
}
impl std::error::Error for PayloadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Dynamics(source) => Some(source),
            Self::Constraint(source) => Some(source),
            Self::Bilateral(source) => Some(source),
            Self::Reconfiguration(source) => Some(source),
            _ => None,
        }
    }
}
impl From<ArticulatedError> for PayloadError {
    fn from(error: ArticulatedError) -> Self {
        Self::Dynamics(error)
    }
}
impl From<GaError> for PayloadError {
    fn from(error: GaError) -> Self {
        Self::Dynamics(ArticulatedError::Geometry(error))
    }
}
impl From<ReconfigurationError> for PayloadError {
    fn from(error: ReconfigurationError) -> Self {
        match error {
            ReconfigurationError::Cancelled => Self::Cancelled,
            other => Self::Reconfiguration(other),
        }
    }
}
impl From<ConstraintError> for PayloadError {
    fn from(error: ConstraintError) -> Self {
        match error {
            ConstraintError::Cancelled
            | ConstraintError::Reconfiguration(ReconfigurationError::Cancelled) => Self::Cancelled,
            other => Self::Constraint(other),
        }
    }
}

impl From<BilateralLoopError> for PayloadError {
    fn from(error: BilateralLoopError) -> Self {
        match error {
            BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::Cancelled,
            )) => Self::Cancelled,
            other => Self::Bilateral(other),
        }
    }
}

/// Capture a free payload rigidly at its actual pose on any target link.
///
/// This instantaneous, perfectly inelastic event conserves total world linear
/// and angular momentum. The original zero-bias articulated mass solve supplies
/// all base/joint rates, including recoil. Optional admitted constant gear rows
/// react on the loaded inertia before any final speed or energy admission.
/// Their initial positions and rates must already satisfy the declared rows.
///
/// Energy loss is evaluated as the sum of nonnegative old-body relative-motion
/// quadratics, never as a subtraction of large lab-frame kinetic energies. The
/// separate balance uses actual returned states in one inertial translating
/// frame and retains its finite residual. No geometry, input state or budget is
/// mutated. Native passes admit at most 4096 links; additional loops poll every
/// 64 entries, and cancellation is checked between every native tree pass.
pub fn capture_payload(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target_link: usize,
    payload: PayloadState,
    constraints: &[JointConstraint],
    controls: PayloadControls<'_>,
    cancelled: impl FnMut() -> bool,
) -> Result<PayloadCapture, PayloadError> {
    let result = capture_payloads(
        model,
        state,
        &[PayloadCaptureInput {
            target_link,
            payload,
        }],
        constraints,
        controls,
        cancelled,
    )?;
    Ok(single_capture(result))
}

fn single_capture(result: PayloadBatchCapture) -> PayloadCapture {
    // The singleton entry retains the same native batch solve and accounting;
    // only the one-element inventory/result containers are unwrapped.
    let PayloadBatchCapture {
        model_after,
        state_after,
        mut attachments,
        mut payloads_after,
        removed_relative_energy_j,
        balance,
    } = result;
    PayloadCapture {
        model_after,
        state_after,
        attachment: attachments.remove(0),
        payload_after: payloads_after.remove(0),
        removed_relative_energy_j,
        balance,
    }
}

fn reconfiguration_controls(controls: PayloadControls<'_>) -> ReconfigurationControls<'_> {
    ReconfigurationControls {
        relative_momentum_tolerance: controls.relative_momentum_tolerance,
        angular_momentum_tolerance_nms: controls.angular_momentum_tolerance_nms,
        linear_momentum_tolerance_ns: controls.linear_momentum_tolerance_ns,
        joint_momentum_tolerances: controls.joint_momentum_tolerances,
        maximum_parameter_work_j: controls.maximum_energy_exchange_j,
    }
}

fn validate_inputs(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    target_link: usize,
    rows: &[JointConstraint],
    controls: PayloadControls<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Kinematics, PayloadError> {
    poll(cancelled)?;
    if target_link >= model.link_count() {
        return Err(PayloadError::Invalid("payload target link index"));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(ReconfigurationError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    model.free_floating_complexity()?;
    reconfigure::validate_controls(model, reconfiguration_controls(controls))?;
    if !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j < 0.0 {
        return Err(PayloadError::Invalid("payload energy residual tolerance"));
    }
    constraints::validate_controls(ConstraintControls {
        reconfiguration: reconfiguration_controls(controls),
        maximum_reaction_work_j: controls.maximum_constraint_work_j,
        minimum_reciprocal_condition: controls.minimum_reciprocal_condition,
    })?;
    let kin = state_kinematics(model, state)?;
    poll(cancelled)?;
    admit_constraints(
        rows,
        &state.joint_position,
        &state.joint_velocity,
        cancelled,
    )?;
    Ok(kin)
}

fn replace_link_inertia(
    model: &ArticulatedModel,
    target: usize,
    inertia: SpatialInertia,
) -> Result<ArticulatedModel, PayloadError> {
    let mut links = model.links().to_vec();
    let old = links
        .get(target)
        .ok_or(PayloadError::Invalid("payload target link index"))?;
    links[target] = Link::new(
        old.name(),
        old.parent(),
        old.parent_from_child_zero(),
        old.joint(),
        inertia,
    );
    Ok(ArticulatedModel::new(links)?)
}

fn combine_inertias(
    retained: SpatialInertia,
    payload: SpatialInertia,
    link_from_payload: Se3,
) -> Result<SpatialInertia, PayloadError> {
    let mass = finite(retained.mass() + payload.mass(), "combined payload mass")?;
    if mass <= retained.mass() || mass <= payload.mass() {
        return Err(PayloadError::Invalid(
            "payload mass increment is not representable",
        ));
    }
    let delta =
        link_from_payload.transform_point(payload.center_of_mass())? - retained.center_of_mass();
    let center = retained.center_of_mass() + delta.scale(payload.mass() / mass);
    let reduced_mass = (retained.mass() / mass) * payload.mass();
    let u = delta.scale(reduced_mass.sqrt());
    let shift = [
        u.y * u.y + u.z * u.z,
        -u.x * u.y,
        -u.x * u.z,
        -u.x * u.y,
        u.x * u.x + u.z * u.z,
        -u.y * u.z,
        -u.x * u.z,
        -u.y * u.z,
        u.x * u.x + u.y * u.y,
    ];
    let rotation = link_from_payload.rotation().matrix();
    let rotated = rotation
        .compose(payload.inertia_com())
        .compose(rotation.transpose());
    let retained_tensor = retained.inertia_com();
    let mut tensor = Mat3 { m: [0.0; 9] };
    for i in 0..3 {
        for j in i..3 {
            // Accumulate each symmetric entry once. The reduced-mass form
            // avoids subtracting two large origin tensors to obtain the COM.
            let value = finite(
                retained_tensor.m[3 * i + j]
                    + 0.5 * rotated.m[3 * i + j]
                    + 0.5 * rotated.m[3 * j + i]
                    + shift[3 * i + j],
                "combined payload centroidal inertia",
            )?;
            tensor.m[3 * i + j] = value;
            tensor.m[3 * j + i] = value;
        }
    }
    Ok(SpatialInertia::new(mass, center, tensor)?)
}

fn payload_at_link(
    kin: &Kinematics,
    target: usize,
    link_from_payload: Se3,
    inertia: SpatialInertia,
) -> Result<PayloadState, PayloadError> {
    let world_from_body = kin.world_from_link[target].compose(link_from_payload)?;
    let twist_body = kin.body_twist[target].transform_by(&link_from_payload.inverse()?);
    validate_twist(twist_body, "attached payload twist")?;
    Ok(PayloadState {
        inertia,
        world_from_body,
        twist_body,
    })
}

fn account_boost(model: &ArticulatedModel, kin: &Kinematics) -> Result<Vec3, PayloadError> {
    let root_velocity = kin.body_twist[0];
    Ok(kin.world_from_link[0].rotation().rotate(
        root_velocity.linear
            + root_velocity
                .angular
                .cross(model.links()[0].inertia().center_of_mass()),
    )?)
}

struct SystemTotals {
    momentum: ArticulatedMomentum,
    generalized: Vec<f64>,
}

fn system_totals(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    kin: &Kinematics,
    extra: Option<(usize, &PayloadState)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SystemTotals, PayloadError> {
    system_totals_many(model, state, kin, extra, cancelled)
}

fn system_totals_many<'p>(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    kin: &Kinematics,
    extras: impl IntoIterator<Item = (usize, &'p PayloadState)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SystemTotals, PayloadError> {
    poll(cancelled)?;
    let momentum = momentum_from_kinematics(model, kin)?;
    let mut world = SpatialSum::default();
    world.add_world_wrench(
        momentum.linear_momentum_world,
        momentum.angular_momentum_world,
        &[],
    );
    let mut kinetic = Dd::from_f64(momentum.kinetic_energy);
    poll(cancelled)?;
    let mut sums = Vec::with_capacity(model.link_count());
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let value = link.inertia().momentum(kin.body_twist[index])?;
        let mut sum = SpatialSum::default();
        sum.add_world_wrench(value.force.scale(-1.0), value.torque.scale(-1.0), &[]);
        sums.push(sum);
    }
    for (index, (target, payload)) in extras.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let body = payload.inertia.momentum(payload.twist_body)?;
        let at_world = body.transform_by(&payload.world_from_body);
        world.add_world_wrench(at_world.force, at_world.torque, &[]);
        kinetic = kinetic + Dd::from_f64(payload.inertia.kinetic_energy(payload.twist_body)?);
        let link_from_payload = kin.world_from_link[target]
            .inverse()?
            .compose(payload.world_from_body)?;
        let at_link = body.transform_by(&link_from_payload);
        sums[target].add_world_wrench(at_link.force.scale(-1.0), at_link.torque.scale(-1.0), &[]);
    }
    let world = world.checked_wrench("payload world momentum")?;
    let momentum = ArticulatedMomentum {
        angular_momentum_world: world.torque,
        linear_momentum_world: world.force,
        kinetic_energy: finite(kinetic.to_f64(), "payload lab kinetic energy")?,
    };
    let mut external = Vec::with_capacity(sums.len());
    for (index, sum) in sums.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        external.push(sum.checked_wrench("payload projected link momentum")?);
    }
    poll(cancelled)?;
    // Zero-velocity RNEA projects independently evaluated body momenta through
    // the original joint subspaces. Each incoming body contributes once at its
    // target's virtual link; the native owner subtracts external wrenches.
    let zero = vec![0.0; model.dof_count()];
    let evaluated = inverse_dynamics(
        model,
        BaseState::prescribed(state.base.world_from_base, Twist::zero(), Twist::zero()),
        &state.joint_position,
        &zero,
        &zero,
        Vec3::default(),
        &external,
    )?;
    poll(cancelled)?;
    let base = evaluated.body_wrench[0].transform_by(&model.links()[0].parent_from_child_zero());
    let mut generalized = Vec::with_capacity(6 + model.dof_count());
    generalized.extend_from_slice(&base.to_array());
    generalized.extend_from_slice(&evaluated.generalized_force);
    validate_slice(&generalized, "payload generalized momentum")?;
    Ok(SystemTotals {
        momentum,
        generalized,
    })
}

fn twist_in_frame(twist: Twist, pose: Se3, boost: Vec3) -> Result<Twist, PayloadError> {
    Ok(twist.minus(Twist::new(
        Vec3::default(),
        pose.rotation().inverse().rotate(boost)?,
    )))
}

fn energy_in_frame(
    model: &ArticulatedModel,
    kin: &Kinematics,
    extra: Option<&PayloadState>,
    boost: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, PayloadError> {
    energy_in_frame_many(model, kin, extra, boost, cancelled)
}

fn energy_in_frame_many<'p>(
    model: &ArticulatedModel,
    kin: &Kinematics,
    extras: impl IntoIterator<Item = &'p PayloadState>,
    boost: Vec3,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, PayloadError> {
    let mut total = Dd::ZERO;
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        total = total
            + Dd::from_f64(link.inertia().kinetic_energy(twist_in_frame(
                kin.body_twist[index],
                kin.world_from_link[index],
                boost,
            )?)?);
    }
    for (index, payload) in extras.into_iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        total = total
            + Dd::from_f64(payload.inertia.kinetic_energy(twist_in_frame(
                payload.twist_body,
                payload.world_from_body,
                boost,
            )?)?);
    }
    finite(total.to_f64(), "payload common-frame kinetic energy")
}

#[allow(clippy::too_many_arguments)]
fn finish_common_frame_balance(
    before: SystemTotals,
    after: SystemTotals,
    impulse: Vec<f64>,
    multipliers: Vec<f64>,
    reaction_work_j: f64,
    boost: Vec3,
    kinetic_before_j: f64,
    kinetic_after_j: f64,
    energy_residual_j: f64,
    controls: PayloadControls<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PayloadBalance, PayloadError> {
    poll(cancelled)?;
    if before.generalized.len() != 6 + controls.joint_momentum_tolerances.len()
        || before.generalized.len() != after.generalized.len()
        || before.generalized.len() != impulse.len()
    {
        return Err(PayloadError::Invalid("payload momentum audit dimensions"));
    }
    finite(reaction_work_j, "payload reaction work")?;
    if reaction_work_j.abs() > controls.maximum_constraint_work_j {
        return Err(PayloadError::ConstraintWork {
            work_j: reaction_work_j,
            maximum_j: controls.maximum_constraint_work_j,
        });
    }
    let mut expected = Vec::with_capacity(before.generalized.len());
    for (index, (&old, &reaction)) in before.generalized.iter().zip(&impulse).enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        expected.push(finite(
            old + reaction,
            "payload reaction-adjusted momentum",
        )?);
    }
    let maximum_normalized_momentum_residual = audit_momentum(
        &expected,
        &after.generalized,
        before.momentum,
        after.momentum,
        reconfiguration_controls(controls),
    )?;
    finite(
        kinetic_before_j,
        "payload common-frame kinetic energy before",
    )?;
    finite(kinetic_after_j, "payload common-frame kinetic energy after")?;
    finite(energy_residual_j, "payload common-frame energy residual")?;
    if energy_residual_j.abs() > controls.energy_tolerance_j {
        return Err(PayloadError::EnergyResidual {
            residual_j: energy_residual_j,
            maximum_j: controls.energy_tolerance_j,
        });
    }
    poll(cancelled)?;
    Ok(PayloadBalance {
        momentum_before: before.momentum,
        momentum_after: after.momentum,
        generalized_momentum_before: before.generalized,
        generalized_momentum_after: after.generalized,
        generalized_constraint_impulse: impulse,
        multipliers,
        reaction_work_j,
        maximum_normalized_momentum_residual,
        accounting_boost_world: boost,
        kinetic_before_j,
        kinetic_after_j,
        energy_residual_j,
    })
}

fn add_wrench(left: Wrench, right: Wrench) -> Wrench {
    Wrench::new(left.torque + right.torque, left.force + right.force)
}

fn finite(value: f64, field: &'static str) -> Result<f64, PayloadError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PayloadError::Invalid(field))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PayloadError> {
    if cancelled() {
        Err(PayloadError::Cancelled)
    } else {
        Ok(())
    }
}
