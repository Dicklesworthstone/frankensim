//! Owner-composed Unitree G1 walking objective for the browser flagship.
//!
//! This L6 module is orchestration, not a robot-math island. The source-bound
//! model, 5,040-D feature map, residual policy, kinematics, and articulated-body
//! solve come from `fs-mbd`; `fs-ga` owns every pose/twist/wrench operation;
//! `fs-time` advances the base on SE(3); `fs-contact` supplies the compliant
//! normal law; and `fs-tribo` supplies dry friction. This module only declares
//! the experiment, composes those receipts in a fixed order, and scores the
//! resulting trajectory.

use core::fmt;

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchGeometry,
    NormalPatchIdentity, NormalPatchLaw, NormalPatchReceipt, NormalPatchRequest,
};
use fs_ga::{Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::{
    BaseState, FreeFloatingBaseState, forward_kinematics, free_floating_forward_dynamics,
};
use fs_mbd::robot_models::{
    CatalogRobotModel, G1_POLICY_ACTUATORS, G1PolicyObservation, G1ResidualPolicy,
    unitree_g1_lower_body_waist_15dof,
};
use fs_time::{RenormPolicy, se3_exp_step_renorm};
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialSlip,
};

/// Stable identity of the owner-composed walking experiment.
pub const G1_WALKING_MODEL_ID: &str = "fs-cmaes/g1-walking-owner-composition-v3";
/// Links retained by the source-bound lower-body-and-waist catalog.
pub const G1_LINK_COUNT: usize = 16;
/// Scalar pose words per link: translation xyz followed by quaternion wxyz.
pub const G1_LINK_POSE_WORDS: usize = 7;

const LEFT_FOOT_LINK: usize = 6;
const RIGHT_FOOT_LINK: usize = 12;
// The source catalog deliberately omits visual/collision meshes, so the
// experiment owns a small, explicit support footprint instead of pretending
// that one point under each ankle represents a foot. The four equal-height
// points form the only admitted ground patches; their moments are accumulated
// onto the source foot link before the articulated-body solve.
const FOOT_CONTACT_POINTS_BODY_M: [Vec3; 4] = [
    Vec3 {
        x: -0.045,
        y: -0.035,
        z: -0.032,
    },
    Vec3 {
        x: -0.045,
        y: 0.035,
        z: -0.032,
    },
    Vec3 {
        x: 0.095,
        y: -0.035,
        z: -0.032,
    },
    Vec3 {
        x: 0.095,
        y: 0.035,
        z: -0.032,
    },
];
const GRAVITY_WORLD_M_PER_S2: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: -9.806_65,
};
const TWO_PI: f64 = 2.0 * core::f64::consts::PI;
const MAX_CONTACT_INDENTATION_M: f64 = 0.035;
const FOOT_EFFECTIVE_RADIUS_M: f64 = 0.035;
const FOOT_REDUCED_MODULUS_PA: f64 = 2.0e6;
// Survival is lexicographically primary for this walking experiment. The
// secondary physical shaping score is smoothly bounded below half of one
// horizon-step charge, so one additional survived step dominates every
// possible difference between shaping scores.
const UNCOMPLETED_STEP_PENALTY: f64 = 1_000.0;
const SHAPING_SCORE_LIMIT: f64 = 400.0;
const SHAPING_SCORE_SCALE: f64 = 10_000.0;

/// Fixed, public experiment controls. They are intentionally not CMA search
/// coordinates: changing them defines a different black-box problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct G1WalkingConfig {
    /// Fixed integrator step [s].
    pub step_s: f64,
    /// Requested physical rollout duration [s].
    pub duration_s: f64,
    /// Commanded forward pelvis speed [m/s].
    pub target_forward_speed_m_per_s: f64,
    /// Open-loop phase frequency supplying the periodic policy basis [Hz].
    pub gait_frequency_hz: f64,
    /// State intervals between trace samples; objective-only runs retain none.
    pub trace_stride: usize,
}

impl Default for G1WalkingConfig {
    fn default() -> Self {
        Self {
            step_s: 1.0 / 480.0,
            duration_s: 1.5,
            target_forward_speed_m_per_s: 0.65,
            gait_frequency_hz: 1.55,
            trace_stride: 12,
        }
    }
}

/// One owner-derived render sample. No browser-side forward kinematics is
/// required or permitted by this contract.
#[derive(Debug, Clone, PartialEq)]
pub struct G1TraceSample {
    /// Physical sample time [s].
    pub time_s: f64,
    /// World-from-link poses in exact catalog link order, encoded xyz+wxyz.
    pub link_pose: [[f64; G1_LINK_POSE_WORDS]; G1_LINK_COUNT],
    /// Active left/right ground patches.
    pub foot_contact: [bool; 2],
}

/// Why a rollout ended. A constitutive-domain or state guard is an evaluated
/// black-box outcome, not a transport failure: CMA-ES must be able to rank the
/// remaining valid candidates instead of losing an entire generation.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G1TerminationReason {
    /// Every requested fixed step completed without a terminal guard.
    Horizon = 0,
    /// The pelvis crossed the declared minimum height.
    BaseHeight = 1,
    /// The pelvis tilt crossed the declared upright envelope.
    BaseTilt = 2,
    /// A foot exceeded the maximum admitted ground indentation.
    ContactIndentation = 3,
    /// A foot exceeded the admitted normal contact speed.
    ContactSpeed = 4,
    /// The normal law refused the candidate-dependent contact state.
    ContactModel = 5,
    /// A source joint crossed its hard position limit.
    JointPositionLimit = 6,
}

impl G1TerminationReason {
    /// Whether a guard ended the rollout instead of the requested horizon.
    #[must_use]
    pub const fn fell(self) -> bool {
        !matches!(self, Self::Horizon)
    }
}

/// Decomposed non-smooth objective and physical rollout diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct G1WalkingReceipt {
    /// Scalar minimized by CMA-ES.
    pub objective: f64,
    /// Forward pelvis displacement [m].
    pub distance_m: f64,
    /// Integral of squared forward-speed error [m²/s].
    pub speed_error_integral: f64,
    /// Absolute actuator work proxy `integral |tau*qdot| dt` [J].
    pub actuator_work_j: f64,
    /// Integrated tangent slip speed squared while in contact [m²/s].
    pub slip_integral: f64,
    /// Integrated squared tilt/height deviation.
    pub posture_integral: f64,
    /// Integrated squared normalized proximity to source hard limits.
    pub joint_limit_integral: f64,
    /// Integrated squared ground reaction, scaled to kilonewtons.
    pub impact_integral: f64,
    /// Number of fixed steps actually completed.
    pub completed_steps: usize,
    /// Exact horizon or terminal guard that ended the rollout.
    pub termination_reason: G1TerminationReason,
    /// Optional owner-derived trajectory for rendering.
    pub trace: Vec<G1TraceSample>,
}

/// Typed refusal surface for the composed experiment.
#[derive(Debug)]
pub enum G1WalkingError {
    /// A fixed experiment control is non-finite or outside its admitted range.
    InvalidConfig { field: &'static str },
    /// The multibody owner refused a malformed or numerically invalid state.
    Robot(fs_mbd::articulated::ArticulatedError),
    /// The `fs-mbd` policy owner refused the parameter vector or observation.
    Policy(fs_mbd::robot_models::G1PolicyError),
    /// The normal-contact owner refused an inapplicable contact state.
    Contact(fs_contact::normal_patch::NormalPatchError),
    /// The friction owner refused an invalid interface or slip state.
    Friction(fs_tribo::TriboError),
    /// The time-integration owner refused an invalid group step.
    Time(fs_time::Se3Error),
    /// The Lie-group owner refused invalid geometry.
    Geometry(fs_ga::GaError),
    /// The configured sphere/plane request unexpectedly returned a line receipt.
    UnexpectedContactReceipt,
    /// A completed rollout produced a non-finite score.
    NonFiniteObjective,
}

impl fmt::Display for G1WalkingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field } => write!(formatter, "invalid G1 rollout {field}"),
            Self::Robot(error) => write!(formatter, "G1 articulated owner refused: {error}"),
            Self::Policy(error) => write!(formatter, "G1 policy owner refused: {error}"),
            Self::Contact(error) => write!(formatter, "G1 normal-contact owner refused: {error}"),
            Self::Friction(error) => write!(formatter, "G1 friction owner refused: {error}"),
            Self::Time(error) => write!(formatter, "G1 time owner refused: {error}"),
            Self::Geometry(error) => write!(formatter, "G1 Lie owner refused: {error}"),
            Self::UnexpectedContactReceipt => {
                formatter.write_str("G1 sphere/plane contact returned a non-point receipt")
            }
            Self::NonFiniteObjective => {
                formatter.write_str("G1 rollout produced a non-finite objective")
            }
        }
    }
}

impl core::error::Error for G1WalkingError {}

impl From<fs_mbd::articulated::ArticulatedError> for G1WalkingError {
    fn from(value: fs_mbd::articulated::ArticulatedError) -> Self {
        Self::Robot(value)
    }
}

impl From<fs_mbd::robot_models::G1PolicyError> for G1WalkingError {
    fn from(value: fs_mbd::robot_models::G1PolicyError) -> Self {
        Self::Policy(value)
    }
}

impl From<fs_contact::normal_patch::NormalPatchError> for G1WalkingError {
    fn from(value: fs_contact::normal_patch::NormalPatchError) -> Self {
        Self::Contact(value)
    }
}

impl From<fs_tribo::TriboError> for G1WalkingError {
    fn from(value: fs_tribo::TriboError) -> Self {
        Self::Friction(value)
    }
}

impl From<fs_time::Se3Error> for G1WalkingError {
    fn from(value: fs_time::Se3Error) -> Self {
        Self::Time(value)
    }
}

impl From<fs_ga::GaError> for G1WalkingError {
    fn from(value: fs_ga::GaError) -> Self {
        Self::Geometry(value)
    }
}

/// Reusable single-threaded experiment object. The model and constitutive
/// interface identities are built once and reused across candidate rollouts.
#[derive(Debug, Clone)]
pub struct G1WalkingEvaluator {
    config: G1WalkingConfig,
    catalog: CatalogRobotModel,
    interface: InterfaceSystemRef,
    friction: FrictionLaw,
    contact_frame: ContactFrame,
    reference_position: [f64; G1_POLICY_ACTUATORS],
    initial_base_height_m: f64,
    step_count: usize,
}

impl G1WalkingEvaluator {
    /// Admit the fixed experiment and build the source-bound model once.
    pub fn new(config: G1WalkingConfig) -> Result<Self, G1WalkingError> {
        validate_config(config)?;
        let catalog = unitree_g1_lower_body_waist_15dof()?;
        let reference_position = reference_posture();
        let initial_base_height_m = initial_height(&catalog, &reference_position)?;
        let interface = InterfaceSystemRef::new(
            "g1-rubber-foot--rigid-dry-ground",
            "g1-walking-rollout-v1",
            "caller-declared-browser-demo-interface",
            InputAuthority::Estimated,
            InterfaceMedium::Dry,
        )?;
        let friction = FrictionLaw::Stribeck {
            static_mu: 0.85,
            kinetic_mu: 0.68,
            characteristic_speed: 0.08,
            viscous_per_speed: 0.015,
        };
        let contact_frame = ContactFrame::new([0.0, 0.0, 1.0])?;
        let step_count = rounded_step_count(config)?;
        Ok(Self {
            config,
            catalog,
            interface,
            friction,
            contact_frame,
            reference_position,
            initial_base_height_m,
            step_count,
        })
    }

    /// Fixed controls admitted by this evaluator.
    #[must_use]
    pub const fn config(&self) -> G1WalkingConfig {
        self.config
    }

    /// Evaluate one candidate without retaining a trajectory.
    pub fn evaluate(&self, parameters: &[f64]) -> Result<G1WalkingReceipt, G1WalkingError> {
        self.rollout(parameters, false)
    }

    /// Evaluate one candidate and retain decimated owner-derived link poses.
    pub fn trace(&self, parameters: &[f64]) -> Result<G1WalkingReceipt, G1WalkingError> {
        self.rollout(parameters, true)
    }

    fn rollout(
        &self,
        parameters: &[f64],
        retain_trace: bool,
    ) -> Result<G1WalkingReceipt, G1WalkingError> {
        let policy = G1ResidualPolicy::new(parameters)?;
        let model = self.catalog.model();
        let mut position = self.reference_position;
        let mut velocity = [0.0; G1_POLICY_ACTUATORS];
        let mut base = FreeFloatingBaseState::stationary(Se3::from_parts(
            fs_ga::So3::identity(),
            Vec3::new(0.0, 0.0, self.initial_base_height_m),
        )?);
        let initial_x = base.world_from_base.translation().x;
        let mut normal_request = normal_request(self.interface.clone(), self.config.step_s);
        let mut trace = if retain_trace {
            Vec::with_capacity(self.step_count / self.config.trace_stride + 2)
        } else {
            Vec::new()
        };
        // The symmetric reference posture begins at the analytically derived
        // static Hertz indentation. This supplies weight support on the first
        // step instead of injecting a timestep-dependent free-fall impact.
        let mut contact = [true; 2];
        let mut speed_error_integral = 0.0;
        let mut actuator_work_j = 0.0;
        let mut slip_integral = 0.0;
        let mut posture_integral = 0.0;
        let mut joint_limit_integral = 0.0;
        let mut impact_integral = 0.0;
        let mut completed_steps = 0;
        let mut termination_reason = G1TerminationReason::Horizon;
        let mut terminal_guard_penalty = 0.0;

        'rollout: for step in 0..self.step_count {
            let time_s = step as f64 * self.config.step_s;
            let kinematics = forward_kinematics(
                model,
                BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero()),
                &position,
                &velocity,
            )?;
            if retain_trace && step % self.config.trace_stride == 0 {
                trace.push(trace_sample(time_s, &kinematics, contact));
            }
            let rotation = base.world_from_base.rotation();
            let gravity_direction_body = rotation.inverse().rotate(Vec3::new(0.0, 0.0, -1.0))?;
            let target_velocity_body = rotation.inverse().rotate(Vec3::new(
                self.config.target_forward_speed_m_per_s,
                0.0,
                0.0,
            ))?;
            let observation = G1PolicyObservation {
                joint_position_rad: position,
                joint_velocity_rad_per_s: velocity,
                gravity_direction_body,
                angular_velocity_body_rad_per_s: base.twist_body.angular,
                target_velocity_error_body_m_per_s: target_velocity_body - base.twist_body.linear,
                foot_contact: contact,
                phase_rad: TWO_PI * self.config.gait_frequency_hz * time_s,
            };
            let residual = policy.evaluate(&observation)?;
            let mut external = [Wrench::default(); G1_LINK_COUNT];
            let mut next_contact = [false; 2];
            for (foot, link) in [LEFT_FOOT_LINK, RIGHT_FOOT_LINK].into_iter().enumerate() {
                let pose = kinematics.world_from_link[link];
                for point_body in FOOT_CONTACT_POINTS_BODY_M {
                    let point_world = pose.transform_point(point_body)?;
                    let point_velocity_body = kinematics.body_twist[link].linear
                        + kinematics.body_twist[link].angular.cross(point_body);
                    let point_velocity_world = pose.rotation().rotate(point_velocity_body)?;
                    let indentation_m = (-point_world.z).max(0.0);
                    if indentation_m == 0.0 {
                        continue;
                    }
                    if indentation_m > MAX_CONTACT_INDENTATION_M {
                        termination_reason = G1TerminationReason::ContactIndentation;
                        terminal_guard_penalty +=
                            220.0 + 120.0 * indentation_m / MAX_CONTACT_INDENTATION_M;
                        break 'rollout;
                    }
                    if point_velocity_world.z.abs() > 8.0 {
                        termination_reason = G1TerminationReason::ContactSpeed;
                        terminal_guard_penalty += 260.0 + 10.0 * point_velocity_world.z.abs();
                        break 'rollout;
                    }
                    next_contact[foot] = true;
                    normal_request.indentation_m = indentation_m;
                    normal_request.indentation_rate_m_per_s = -point_velocity_world.z;
                    let normal_force_n = match normal_request.evaluate() {
                        Ok(NormalPatchReceipt::Point(receipt)) => receipt.normal_force_n,
                        Ok(NormalPatchReceipt::Line(_)) => {
                            return Err(G1WalkingError::UnexpectedContactReceipt);
                        }
                        Err(_) => {
                            termination_reason = G1TerminationReason::ContactModel;
                            terminal_guard_penalty += 300.0;
                            break 'rollout;
                        }
                    };
                    let slip = TangentialSlip::new(
                        &self.contact_frame,
                        [point_velocity_world.x, point_velocity_world.y, 0.0],
                    )?;
                    let friction = self
                        .friction
                        .evaluate(&self.interface, normal_force_n, slip)?;
                    let traction = friction.traction_n();
                    let force_world = Vec3::new(traction[0], traction[1], normal_force_n);
                    let force_body = pose.rotation().inverse().rotate(force_world)?;
                    let previous = external[link];
                    external[link] = Wrench::new(
                        previous.torque + point_body.cross(force_body),
                        previous.force + force_body,
                    );
                    slip_integral += (point_velocity_world.x * point_velocity_world.x
                        + point_velocity_world.y * point_velocity_world.y)
                        * self.config.step_s;
                    impact_integral += (normal_force_n / 1_000.0).powi(2) * self.config.step_s;
                }
            }
            contact = next_contact;

            let generalized_force = controller_force(
                &self.catalog,
                &self.reference_position,
                &position,
                &velocity,
                &residual,
            );

            let dynamics = free_floating_forward_dynamics(
                model,
                base,
                &position,
                &velocity,
                &generalized_force,
                GRAVITY_WORLD_M_PER_S2,
                &external,
            )?;
            for actuator in 0..G1_POLICY_ACTUATORS {
                actuator_work_j +=
                    (generalized_force[actuator] * velocity[actuator]).abs() * self.config.step_s;
                let source = self.catalog.joints()[actuator];
                let next_velocity = velocity[actuator]
                    + dynamics.generalized_acceleration[actuator] * self.config.step_s;
                let velocity_limit = source.velocity_rad_per_second;
                if next_velocity.abs() > velocity_limit {
                    let normalized_overshoot =
                        (next_velocity.abs() - velocity_limit) / velocity_limit;
                    joint_limit_integral +=
                        25.0 * normalized_overshoot * normalized_overshoot * self.config.step_s;
                }
                velocity[actuator] = next_velocity.clamp(-velocity_limit, velocity_limit);
                let next_position = position[actuator] + velocity[actuator] * self.config.step_s;
                if next_position < source.lower_position_rad
                    || next_position > source.upper_position_rad
                {
                    let half_range = 0.5 * (source.upper_position_rad - source.lower_position_rad);
                    let overshoot = if next_position < source.lower_position_rad {
                        source.lower_position_rad - next_position
                    } else {
                        next_position - source.upper_position_rad
                    };
                    termination_reason = G1TerminationReason::JointPositionLimit;
                    terminal_guard_penalty += 200.0 + 100.0 * overshoot / half_range;
                    break 'rollout;
                }
                position[actuator] = next_position;
                let center = 0.5 * (source.lower_position_rad + source.upper_position_rad);
                let half_range = 0.5 * (source.upper_position_rad - source.lower_position_rad);
                let normalized = (position[actuator] - center) / half_range;
                joint_limit_integral += normalized.powi(8) * self.config.step_s;
            }
            base.twist_body = base.twist_body.plus(
                dynamics
                    .base_spatial_acceleration_body
                    .scale(self.config.step_s),
            );
            base.world_from_base = se3_exp_step_renorm(
                base.world_from_base,
                base.twist_body,
                self.config.step_s,
                &RenormPolicy::default(),
            )?
            .0;
            completed_steps = step + 1;

            let updated_rotation = base.world_from_base.rotation();
            let updated_gravity_direction_body = updated_rotation
                .inverse()
                .rotate(Vec3::new(0.0, 0.0, -1.0))?;
            let world_velocity = updated_rotation.rotate(base.twist_body.linear)?;
            let speed_error = world_velocity.x - self.config.target_forward_speed_m_per_s;
            speed_error_integral += speed_error * speed_error * self.config.step_s;
            let height_error = base.world_from_base.translation().z - self.initial_base_height_m;
            let tilt_error = 1.0 + updated_gravity_direction_body.z;
            posture_integral += (2.5 * height_error * height_error + 4.0 * tilt_error * tilt_error)
                * self.config.step_s;
            if base.world_from_base.translation().z < 0.32 {
                termination_reason = G1TerminationReason::BaseHeight;
                break;
            }
            if updated_gravity_direction_body.z > -0.35 {
                termination_reason = G1TerminationReason::BaseTilt;
                break;
            }
        }

        if retain_trace {
            let kinematics = forward_kinematics(
                model,
                BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero()),
                &position,
                &velocity,
            )?;
            trace.push(trace_sample(
                completed_steps as f64 * self.config.step_s,
                &kinematics,
                contact,
            ));
        }
        let distance_m = base.world_from_base.translation().x - initial_x;
        // Falling is one distinct failure in addition to every skipped step.
        // This keeps the ordering strict even at the horizon boundary: a fall
        // on the final step is worse than completing the horizon upright, and
        // a fall one step earlier is worse again.
        let failed_horizon_steps =
            survival_charge_steps(self.step_count, completed_steps, termination_reason.fell());
        let raw_shaping_score = -18.0 * distance_m
            + 12.0 * speed_error_integral
            + 0.008 * actuator_work_j
            + 16.0 * slip_integral
            + 30.0 * posture_integral
            + 2.0 * joint_limit_integral
            + 0.8 * impact_integral
            + terminal_guard_penalty;
        if !raw_shaping_score.is_finite() {
            return Err(G1WalkingError::NonFiniteObjective);
        }
        let bounded_shaping_score =
            SHAPING_SCORE_LIMIT * (raw_shaping_score / SHAPING_SCORE_SCALE).tanh();
        let objective =
            UNCOMPLETED_STEP_PENALTY * failed_horizon_steps as f64 + bounded_shaping_score;
        if !objective.is_finite() {
            return Err(G1WalkingError::NonFiniteObjective);
        }
        Ok(G1WalkingReceipt {
            objective,
            distance_m,
            speed_error_integral,
            actuator_work_j,
            slip_integral,
            posture_integral,
            joint_limit_integral,
            impact_integral,
            completed_steps,
            termination_reason,
            trace,
        })
    }
}

const fn survival_charge_steps(total_steps: usize, completed_steps: usize, fell: bool) -> usize {
    debug_assert!(completed_steps <= total_steps);
    total_steps - completed_steps + fell as usize
}

fn validate_config(config: G1WalkingConfig) -> Result<(), G1WalkingError> {
    for (value, field) in [
        (config.step_s, "step_s"),
        (config.duration_s, "duration_s"),
        (
            config.target_forward_speed_m_per_s,
            "target_forward_speed_m_per_s",
        ),
        (config.gait_frequency_hz, "gait_frequency_hz"),
    ] {
        if !value.is_finite() {
            return Err(G1WalkingError::InvalidConfig { field });
        }
    }
    if config.step_s <= 0.0 {
        return Err(G1WalkingError::InvalidConfig { field: "step_s" });
    }
    if config.duration_s <= 0.0 {
        return Err(G1WalkingError::InvalidConfig {
            field: "duration_s",
        });
    }
    if config.target_forward_speed_m_per_s < 0.0 {
        return Err(G1WalkingError::InvalidConfig {
            field: "target_forward_speed_m_per_s",
        });
    }
    if config.gait_frequency_hz <= 0.0 {
        return Err(G1WalkingError::InvalidConfig {
            field: "gait_frequency_hz",
        });
    }
    if config.trace_stride == 0 {
        return Err(G1WalkingError::InvalidConfig {
            field: "trace_stride",
        });
    }
    Ok(())
}

fn rounded_step_count(config: G1WalkingConfig) -> Result<usize, G1WalkingError> {
    let count = (config.duration_s / config.step_s).round();
    if !(count.is_finite() && (1.0..=10_000.0).contains(&count)) {
        return Err(G1WalkingError::InvalidConfig {
            field: "duration_s / step_s",
        });
    }
    Ok(count as usize)
}

const fn reference_posture() -> [f64; G1_POLICY_ACTUATORS] {
    [
        -0.20, 0.0, 0.0, 0.42, -0.22, 0.0, -0.20, 0.0, 0.0, 0.42, -0.22, 0.0, 0.0, 0.0, 0.0,
    ]
}

fn initial_height(
    catalog: &CatalogRobotModel,
    position: &[f64; G1_POLICY_ACTUATORS],
) -> Result<f64, G1WalkingError> {
    let kinematics = forward_kinematics(
        catalog.model(),
        BaseState::stationary(Se3::identity()),
        position,
        &[0.0; G1_POLICY_ACTUATORS],
    )?;
    let mut minimum_contact_z = f64::INFINITY;
    for link in [LEFT_FOOT_LINK, RIGHT_FOOT_LINK] {
        for point_body in FOOT_CONTACT_POINTS_BODY_M {
            let point = kinematics.world_from_link[link].transform_point(point_body)?;
            minimum_contact_z = minimum_contact_z.min(point.z);
        }
    }
    Ok(-minimum_contact_z - static_contact_indentation_m(catalog))
}

fn static_contact_indentation_m(catalog: &CatalogRobotModel) -> f64 {
    let total_mass_kg = catalog
        .model()
        .links()
        .iter()
        .map(|link| link.inertia().mass())
        .sum::<f64>();
    let contact_count = 2.0 * FOOT_CONTACT_POINTS_BODY_M.len() as f64;
    let force_per_patch_n = total_mass_kg * GRAVITY_WORLD_M_PER_S2.z.abs() / contact_count;
    let hertz_coefficient = (4.0 / 3.0) * FOOT_REDUCED_MODULUS_PA * FOOT_EFFECTIVE_RADIUS_M.sqrt();
    (force_per_patch_n / hertz_coefficient).powf(2.0 / 3.0)
}

fn controller_force(
    catalog: &CatalogRobotModel,
    reference: &[f64; G1_POLICY_ACTUATORS],
    position: &[f64; G1_POLICY_ACTUATORS],
    velocity: &[f64; G1_POLICY_ACTUATORS],
    residual: &[f64; G1_POLICY_ACTUATORS],
) -> [f64; G1_POLICY_ACTUATORS] {
    let mut force = [0.0; G1_POLICY_ACTUATORS];
    for actuator in 0..G1_POLICY_ACTUATORS {
        let effort_limit = catalog.joints()[actuator].effort_newton_metres;
        let proportional_gain = if matches!(actuator, 3 | 9) {
            95.0
        } else {
            58.0
        };
        let derivative_gain = if matches!(actuator, 3 | 9) { 3.8 } else { 2.5 };
        let posture = proportional_gain * (reference[actuator] - position[actuator])
            - derivative_gain * velocity[actuator];
        let residual_force = 0.32 * effort_limit * residual[actuator];
        force[actuator] = (posture + residual_force).clamp(-effort_limit, effort_limit);
    }
    force
}

fn normal_request(interface: InterfaceSystemRef, step_s: f64) -> NormalPatchRequest {
    NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: G1_WALKING_MODEL_ID.to_owned(),
            source_id: "estimated-rubber-foot-sphere".to_owned(),
            state_id: "g1-walking-rollout-contact-state".to_owned(),
        },
        interface,
        law: NormalPatchLaw::HuntCrossleySphere {
            effective_radius_m: FOOT_EFFECTIVE_RADIUS_M,
            reduced_modulus_pa: FOOT_REDUCED_MODULUS_PA,
            dissipation_s_per_m: 0.25,
        },
        geometry: NormalPatchGeometry::SpherePlane,
        indentation_m: 0.0,
        indentation_rate_m_per_s: 0.0,
        step_s,
        line_load_n_per_m: 0.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 0.12,
            layer_thickness_m: 0.06,
            yield_strength_pa: 20.0e6,
            characteristic_rate_m_per_s: 8.0,
            temperature_k: 293.15,
            adhesion_energy_j_per_m2: 0.0,
        },
        limits: ApplicabilityLimits {
            max_patch_to_radius: 1.1,
            max_strain: 1.1,
            max_patch_to_depth: 1.0,
            max_patch_to_layer: 1.0,
            max_pressure_to_yield: 1.0,
            max_rate_ratio: 1.0,
            min_temperature_k: 250.0,
            max_temperature_k: 330.0,
        },
        uncertainty: InputUncertainty {
            radius_relative: 0.20,
            modulus_relative: 0.50,
            load_relative: 0.0,
        },
    }
}

fn trace_sample(
    time_s: f64,
    kinematics: &fs_mbd::articulated::Kinematics,
    foot_contact: [bool; 2],
) -> G1TraceSample {
    let mut link_pose = [[0.0; G1_LINK_POSE_WORDS]; G1_LINK_COUNT];
    for (output, pose) in link_pose.iter_mut().zip(&kinematics.world_from_link) {
        let translation = pose.translation();
        let rotation = pose.rotation();
        let quaternion = rotation.as_quat();
        *output = [
            translation.x,
            translation.y,
            translation.z,
            quaternion.w,
            quaternion.x,
            quaternion.y,
            quaternion.z,
        ];
    }
    G1TraceSample {
        time_s,
        link_pose,
        foot_contact,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_policy_rollout_is_deterministic_and_owner_derived() -> Result<(), G1WalkingError> {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        let parameters = vec![0.0; fs_mbd::robot_models::G1_POLICY_DIMENSION];
        let first = evaluator.trace(&parameters)?;
        let second = evaluator.trace(&parameters)?;
        assert_eq!(first, second);
        assert!(first.trace.len() >= 5);
        assert_eq!(first.trace[0].link_pose.len(), G1_LINK_COUNT);
        assert!(first.objective.is_finite());
        assert!(
            first.completed_steps >= evaluator.step_count / 2,
            "zero policy terminated as {:?} after {} of {} steps",
            first.termination_reason,
            first.completed_steps,
            evaluator.step_count
        );
        Ok(())
    }

    #[test]
    fn evaluator_refuses_wrong_policy_shape_and_invalid_controls() -> Result<(), G1WalkingError> {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        assert!(matches!(
            evaluator.evaluate(&[]),
            Err(G1WalkingError::Policy(
                fs_mbd::robot_models::G1PolicyError::ParameterCount { .. }
            ))
        ));
        assert!(matches!(
            G1WalkingEvaluator::new(G1WalkingConfig {
                step_s: 0.0,
                ..G1WalkingConfig::default()
            }),
            Err(G1WalkingError::InvalidConfig { field: "step_s" })
        ));
        Ok(())
    }

    #[test]
    fn early_termination_cannot_win_by_escaping_the_remaining_horizon() -> Result<(), G1WalkingError>
    {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        let zero = evaluator.evaluate(&vec![0.0; fs_mbd::robot_models::G1_POLICY_DIMENSION])?;
        let aggressive =
            evaluator.evaluate(&vec![0.03; fs_mbd::robot_models::G1_POLICY_DIMENSION])?;
        assert!(
            aggressive.completed_steps < zero.completed_steps,
            "aggressive {aggressive:?}, zero {zero:?}"
        );
        assert!(
            aggressive.objective > zero.objective,
            "aggressive {aggressive:?}, zero {zero:?}"
        );
        Ok(())
    }

    #[test]
    fn survival_charge_is_strict_through_the_horizon_boundary() {
        assert_eq!(survival_charge_steps(180, 180, false), 0);
        assert_eq!(survival_charge_steps(180, 180, true), 1);
        assert_eq!(survival_charge_steps(180, 179, true), 2);
    }

    #[test]
    fn reference_feet_start_at_static_contact_preload() -> Result<(), G1WalkingError> {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        let kinematics = forward_kinematics(
            evaluator.catalog.model(),
            BaseState::stationary(Se3::from_parts(
                fs_ga::So3::identity(),
                Vec3::new(0.0, 0.0, evaluator.initial_base_height_m),
            )?),
            &evaluator.reference_position,
            &[0.0; G1_POLICY_ACTUATORS],
        )?;
        let expected_z = -static_contact_indentation_m(&evaluator.catalog);
        for link in [LEFT_FOOT_LINK, RIGHT_FOOT_LINK] {
            for point_body in FOOT_CONTACT_POINTS_BODY_M {
                let point = kinematics.world_from_link[link].transform_point(point_body)?;
                assert!((point.z - expected_z).abs() < 1.0e-10);
            }
        }
        Ok(())
    }
}
