//! Rigid payload capture and release on articulated trees.
//!
//! The root API describes free motion; [`supported`] retains fixed or driven
//! support reactions and work during a prescribed-base capture.
//!
//! An admitted attachment preserves the target link frame and every joint.
//! Full rotated centroidal inertia and COM orbital inertia enter the original
//! articulated-body mass solve; existing constant gear rows may react in that
//! same event. The supplied poses define the attachment geometry. Contact
//! discovery, finite-duration grasping, mass flux and a thermal destination for
//! the removed relative kinetic energy remain caller-owned physical models.

/// Conservative payload release with optional internal separation impulses.
#[path = "inventory/release.rs"]
pub mod release;

/// Rigid payload capture with fixed or driven prescribed support momentum/work.
#[path = "inventory/supported.rs"]
pub mod supported;

use super::reconfigure::constraints::{
    self, ConstraintControls, ConstraintError, JointConstraint, PreparedJointConstraints,
    admit_constraints, apply_rows,
};
use super::reconfigure::{
    self, ReconfigurationControls, ReconfigurationError, audit_momentum, state_kinematics,
};
use super::{FreeFloatingMotionState, MAX_FREE_FLOATING_STEP_LINKS};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, FreeFloatingBaseState,
    Kinematics, Link, SpatialInertia, forward_kinematics, free_floating_impulse_delta,
    inverse_dynamics, momentum_from_kinematics, validate_slice, validate_twist, validate_vec3,
};
use core::fmt;
use fs_ga::{GaError, Mat3, Se3, Twist, Vec3, Wrench};

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
    /// Simultaneous internal gear impulse `[0; A^T lambda]`.
    pub generalized_constraint_impulse: Vec<f64>,
    /// Gear-row impulse multipliers in caller row order.
    pub multipliers: Vec<f64>,
    /// Actual `sum lambda * (A v_before + A v_after)/2`, J; never inferred heat.
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
    mut cancelled: impl FnMut() -> bool,
) -> Result<PayloadCapture, PayloadError> {
    let before_kin = validate_inputs(
        model,
        state,
        target_link,
        constraints,
        controls,
        &mut cancelled,
    )?;
    validate_twist(payload.twist_body, "payload.twist_body")?;
    let link_from_payload = before_kin.world_from_link[target_link]
        .inverse()?
        .compose(payload.world_from_body)?;
    let loaded_inertia = combine_inertias(
        model.links()[target_link].inertia(),
        payload.inertia,
        link_from_payload,
    )?;
    let model_after = replace_link_inertia(model, target_link, loaded_inertia)?;
    poll(&mut cancelled)?;
    let before = system_totals(
        model,
        state,
        &before_kin,
        Some((target_link, &payload)),
        &mut cancelled,
    )?;
    let boost = account_boost(model, &before_kin)?;
    let boost_base = state
        .base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(boost)?;
    let centered = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state
                .base
                .twist_body
                .minus(Twist::new(Vec3::default(), boost_base)),
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )?;
    poll(&mut cancelled)?;
    let mut body_momenta = Vec::with_capacity(model.link_count());
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        body_momenta.push(link.inertia().momentum(centered.body_twist[index])?);
    }
    let incoming = payload
        .inertia
        .momentum(twist_in_frame(
            payload.twist_body,
            payload.world_from_body,
            boost,
        )?)?
        .transform_by(&link_from_payload);
    body_momenta[target_link] = add_wrench(body_momenta[target_link], incoming);
    poll(&mut cancelled)?;
    let solved = free_floating_impulse_delta(
        &model_after,
        &state.joint_position,
        &vec![0.0; model.dof_count()],
        &body_momenta,
    )?;
    poll(&mut cancelled)?;
    let mut state_after = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            state.base.world_from_base,
            solved.base.plus(Twist::new(Vec3::default(), boost_base)),
        ),
        joint_position: state.joint_position.clone(),
        joint_velocity: solved.joint,
    };
    let mut impulse = vec![0.0; 6 + model.dof_count()];
    let mut multipliers = Vec::new();
    if !constraints.is_empty() {
        let mut prepared = PreparedJointConstraints::new(
            &model_after,
            &state.joint_position,
            constraints,
            controls.minimum_reciprocal_condition,
            &mut cancelled,
        )?;
        let rhs: Vec<_> = apply_rows(constraints, &state_after.joint_velocity, &mut cancelled)?
            .into_iter()
            .map(|value| -value)
            .collect();
        let response = prepared.solve_response(&rhs, &mut cancelled)?;
        state_after.base.twist_body = state_after.base.twist_body.plus(response.base_delta);
        for (index, delta) in response.joint_delta.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            state_after.joint_velocity[index] = finite(
                state_after.joint_velocity[index] + delta,
                "payload final joint rate",
            )?;
        }
        impulse[6..].copy_from_slice(&response.joint_effort);
        multipliers = response.multipliers;
    }
    let final_rows = admit_constraints(
        constraints,
        &state_after.joint_position,
        &state_after.joint_velocity,
        &mut cancelled,
    )?;
    let after_kin = state_kinematics(&model_after, &state_after)?;
    poll(&mut cancelled)?;
    let payload_after =
        payload_at_link(&after_kin, target_link, link_from_payload, payload.inertia)?;
    let after = system_totals(&model_after, &state_after, &after_kin, None, &mut cancelled)?;

    let mut removed_relative_energy_j = 0.0;
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(&mut cancelled)?;
        }
        removed_relative_energy_j = finite(
            removed_relative_energy_j
                + link.inertia().kinetic_energy(
                    before_kin.body_twist[index].minus(after_kin.body_twist[index]),
                )?,
            "payload removed relative energy",
        )?;
    }
    removed_relative_energy_j = finite(
        removed_relative_energy_j
            + payload
                .inertia
                .kinetic_energy(payload.twist_body.minus(payload_after.twist_body))?,
        "payload removed relative energy",
    )?;
    if removed_relative_energy_j < 0.0 {
        return Err(PayloadError::Invalid("negative capture relative energy"));
    }
    if removed_relative_energy_j > controls.maximum_energy_exchange_j {
        return Err(PayloadError::EnergyExchange {
            energy_j: removed_relative_energy_j,
            maximum_j: controls.maximum_energy_exchange_j,
        });
    }
    let initial_rows = apply_rows(constraints, &state.joint_velocity, &mut cancelled)?;
    let mut reaction_work_j = 0.0;
    for (index, multiplier) in multipliers.iter().enumerate() {
        poll(&mut cancelled)?;
        reaction_work_j = finite(
            reaction_work_j
                + multiplier * (0.5 * initial_rows[index] + 0.5 * final_rows.velocity[index]),
            "payload reaction work",
        )?;
    }
    let kinetic_before_j =
        energy_in_frame(model, &before_kin, Some(&payload), boost, &mut cancelled)?;
    let kinetic_after_j = energy_in_frame(&model_after, &after_kin, None, boost, &mut cancelled)?;
    // The relative quadratic already describes the complete passive capture
    // projection. Finite gear work is independently retained and capped; it is
    // not relabelled as removed energy or subtracted to hide this defect.
    let energy_residual_j = finite(
        kinetic_after_j - kinetic_before_j + removed_relative_energy_j,
        "payload common-frame capture energy residual",
    )?;
    let balance = finish_common_frame_balance(
        before,
        after,
        impulse,
        multipliers,
        reaction_work_j,
        boost,
        kinetic_before_j,
        kinetic_after_j,
        energy_residual_j,
        controls,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    Ok(PayloadCapture {
        model_after,
        state_after,
        attachment: PayloadAttachment {
            target_link,
            retained_link: model.links()[target_link].clone(),
            loaded_inertia,
            link_from_payload,
            payload_inertia: payload.inertia,
        },
        payload_after,
        removed_relative_energy_j,
        balance,
    })
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
    poll(cancelled)?;
    let mut momentum = momentum_from_kinematics(model, kin)?;
    poll(cancelled)?;
    let mut external = Vec::with_capacity(model.link_count());
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        let value = link.inertia().momentum(kin.body_twist[index])?;
        external.push(Wrench::new(
            value.torque.scale(-1.0),
            value.force.scale(-1.0),
        ));
    }
    if let Some((target, payload)) = extra {
        let body = payload.inertia.momentum(payload.twist_body)?;
        let world = body.transform_by(&payload.world_from_body);
        momentum.linear_momentum_world = momentum.linear_momentum_world + world.force;
        momentum.angular_momentum_world = momentum.angular_momentum_world + world.torque;
        validate_vec3(
            momentum.linear_momentum_world,
            "payload total linear momentum",
        )?;
        validate_vec3(
            momentum.angular_momentum_world,
            "payload total angular momentum",
        )?;
        momentum.kinetic_energy = finite(
            momentum.kinetic_energy + payload.inertia.kinetic_energy(payload.twist_body)?,
            "payload lab kinetic energy",
        )?;
        let link_from_payload = kin.world_from_link[target]
            .inverse()?
            .compose(payload.world_from_body)?;
        let at_link = body.transform_by(&link_from_payload);
        external[target] = add_wrench(
            external[target],
            Wrench::new(at_link.torque.scale(-1.0), at_link.force.scale(-1.0)),
        );
    }
    poll(cancelled)?;
    // Zero-velocity RNEA projects the independently evaluated body momenta
    // through native joint subspaces. External wrenches are subtracted there.
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
    let mut total = 0.0;
    for (index, link) in model.links().iter().enumerate() {
        if index % 64 == 0 {
            poll(cancelled)?;
        }
        total = finite(
            total
                + link.inertia().kinetic_energy(twist_in_frame(
                    kin.body_twist[index],
                    kin.world_from_link[index],
                    boost,
                )?)?,
            "payload common-frame kinetic energy",
        )?;
    }
    if let Some(payload) = extra {
        total = finite(
            total
                + payload.inertia.kinetic_energy(twist_in_frame(
                    payload.twist_body,
                    payload.world_from_body,
                    boost,
                )?)?,
            "payload common-frame kinetic energy",
        )?;
    }
    Ok(total)
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
