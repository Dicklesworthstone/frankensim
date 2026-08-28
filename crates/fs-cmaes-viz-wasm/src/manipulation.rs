//! Owner-composed KUKA LBR iiwa household-manipulation objective.
//!
//! This L6 experiment declares a benchmark and composes existing owners. The
//! pinned iiwa topology, inertias, hard limits, Lie-group kinematics, inverse
//! dynamics, and articulated-body forward dynamics come from `fs-mbd` and
//! `fs-ga`; `fs-contact` supplies the compliant finger-pad normal response;
//! and `fs-tribo` supplies the dry-friction capacity. This module owns only the
//! task presets, a disclosed 128-coordinate joint/gripper knot trajectory, a
//! reduced bilateral-grasp state machine, integration, scoring, and receipts.
//!
//! The reduced grasp is intentionally narrower than a general rigid-contact
//! solve. An object is supported by a horizontal surface until two opposing
//! compliant pad contacts establish sufficient owner-reported static friction.
//! While that capacity remains available, a bilateral rigid grasp constraint
//! carries the object with the flange and its weight is applied to the source
//! arm as an external wrench. Release returns the object to ballistic motion
//! and a one-sided horizontal support. There is no mesh collision, self-
//! collision, impact impulse, grasp-planning, or hardware-transfer claim.

use core::fmt;

use fs_contact::normal_patch::{
    ApplicabilityInput, ApplicabilityLimits, InputUncertainty, NormalPatchGeometry,
    NormalPatchIdentity, NormalPatchLaw, NormalPatchReceipt, NormalPatchRequest,
};
use fs_ga::{Se3, Vec3, Wrench};
use fs_mbd::articulated::{BaseState, forward_dynamics, forward_kinematics, inverse_dynamics};
use fs_mbd::robot_models::{CatalogRobotModel, kuka_lbr_iiwa7_r800};
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialSlip,
};

/// Stable identity of the reduced owner-composed manipulation experiment.
pub const MANIPULATION_MODEL_ID: &str = "fs-cmaes/iiwa-household-manipulation-v1";
/// Source-bound joints in the KUKA LBR iiwa 7 R800 catalog.
pub const ARM_JOINT_COUNT: usize = 7;
/// Source-bound links, including fixed `iiwa_link_0`.
pub const ARM_LINK_COUNT: usize = 8;
/// Uniform trajectory knots per joint and for the reduced gripper actuator.
pub const ARM_POLICY_KNOTS: usize = 16;
/// Seven joint rows plus one gripper-width row, each with sixteen knots.
pub const ARM_POLICY_DIMENSION: usize = (ARM_JOINT_COUNT + 1) * ARM_POLICY_KNOTS;
/// Translation xyz followed by quaternion wxyz.
pub const ARM_LINK_POSE_WORDS: usize = 7;

const END_EFFECTOR_LINK: usize = ARM_LINK_COUNT - 1;
const GRAVITY_WORLD_M_PER_S2: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: -9.806_65,
};
/// Maximum/open finger separation admitted by the reduced gripper [m].
pub const OPEN_GRIPPER_WIDTH_M: f64 = 0.105;
/// Minimum finger separation admitted by the reduced gripper [m].
pub const MIN_GRIPPER_WIDTH_M: f64 = 0.020;
const MAX_GRIPPER_SPEED_M_PER_S: f64 = 0.18;
const MAX_PAD_INDENTATION_M: f64 = 0.004;
const PAD_EFFECTIVE_RADIUS_M: f64 = 0.012;
const PAD_REDUCED_MODULUS_PA: f64 = 1.0e6;
/// Maximum terminal object-centre error for a successful placement [m].
pub const PLACEMENT_TOLERANCE_M: f64 = 0.085;
/// Required maximum object-centre rise for successful transport [m].
pub const LIFT_TARGET_M: f64 = 0.09;
const MAX_POLICY_POPULATION: usize = 64;

/// The three household settings share one dynamics and policy contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
pub enum ManipulationTask {
    /// Move a benchmark stoneware mug between kitchen counter stations.
    KitchenMug = 0,
    /// Move a benchmark television remote into a living-room organizer.
    LivingRoomRemote = 1,
    /// Move a benchmark hand trowel into a backyard tool caddy.
    BackyardTrowel = 2,
}

/// Fixed workload controls for one evaluator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ManipulationConfig {
    /// Household scene/object preset.
    pub task: ManipulationTask,
    /// Fixed semi-implicit joint/object step [s].
    pub step_s: f64,
    /// Fixed horizon [s].
    pub duration_s: f64,
    /// Whole steps between retained render samples.
    pub trace_stride: usize,
}

impl Default for ManipulationConfig {
    fn default() -> Self {
        Self {
            task: ManipulationTask::KitchenMug,
            step_s: 1.0 / 90.0,
            duration_s: 4.0,
            trace_stride: 3,
        }
    }
}

/// Object and scene dimensions consumed by the renderer without inference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ManipulationScene {
    /// Benchmark object mass [kg]. These are declared demo estimates.
    pub object_mass_kg: f64,
    /// Full axis-aligned visual/support envelope [m].
    pub object_dimensions_m: Vec3,
    /// Object half-width opposed by the two finger pads [m].
    pub grasp_half_width_m: f64,
    /// Initial object centre in world coordinates [m].
    pub initial_object_position_m: Vec3,
    /// Desired released object centre in world coordinates [m].
    pub goal_object_position_m: Vec3,
    /// Horizontal support height [m].
    pub support_height_m: f64,
    /// Centre of one declared keep-out box [m].
    pub obstacle_center_m: Vec3,
    /// Half extents of the keep-out box [m].
    pub obstacle_half_extents_m: Vec3,
}

/// One decimated owner trajectory sample.
#[derive(Debug, Clone, PartialEq)]
pub struct ManipulationTraceSample {
    /// Simulation time [s].
    pub time_s: f64,
    /// Rate-limited finger separation [m].
    pub gripper_width_m: f64,
    /// Sum of the two compliant-pad normal forces [N].
    pub grip_normal_force_n: f64,
    /// Whether the reduced bilateral grasp constraint is active.
    pub grasped: bool,
    /// World object pose: translation xyz, quaternion wxyz.
    pub object_pose: [f64; ARM_LINK_POSE_WORDS],
    /// World-from-link poses in source catalog order.
    pub link_pose: [[f64; ARM_LINK_POSE_WORDS]; ARM_LINK_COUNT],
}

/// Decomposed objective receipt plus an optional render trace.
#[derive(Debug, Clone, PartialEq)]
pub struct ManipulationReceipt {
    /// Scalar minimized by CMA-ES.
    pub objective: f64,
    /// Final object-centre distance to the goal [m].
    pub final_object_error_m: f64,
    /// Closest pre-grasp flange-to-object distance [m].
    pub minimum_reach_error_m: f64,
    /// Maximum object-centre rise from its supported start [m].
    pub maximum_lift_m: f64,
    /// Absolute actuator work integral [J].
    pub actuator_work_j: f64,
    /// Integrated keep-out-box proximity/penetration proxy [m s].
    pub obstacle_integral: f64,
    /// Integrated proposed joint/gripper limit excess.
    pub control_limit_integral: f64,
    /// First successful bilateral grasp time, or horizon if never grasped [s].
    pub first_grasp_time_s: f64,
    /// Total time under the bilateral grasp constraint [s].
    pub grasp_duration_s: f64,
    /// Peak owner-reported summed finger normal force [N].
    pub peak_grip_force_n: f64,
    /// Whether a grasp was established at least once.
    pub ever_grasped: bool,
    /// Whether the object was released after being transported.
    pub released_after_transport: bool,
    /// Whether the terminal placement met all declared criteria.
    pub placed: bool,
    /// Whole completed fixed steps.
    pub completed_steps: usize,
    /// Trace samples; empty for objective-only evaluation.
    pub trace: Vec<ManipulationTraceSample>,
}

/// Typed refusal surface for configuration, owners, and numeric results.
#[derive(Debug)]
pub enum ManipulationError {
    /// A fixed workload control lies outside the browser envelope.
    InvalidConfig { field: &'static str },
    /// A policy has the wrong flat length.
    ParameterCount { expected: usize, actual: usize },
    /// A policy coordinate is not finite.
    NonFiniteParameter { index: usize },
    /// The source-bound articulated owner refused an input or derived value.
    Robot(fs_mbd::articulated::ArticulatedError),
    /// A canonical Lie-group operation refused.
    Geometry(fs_ga::GaError),
    /// The compliant normal owner refused a pad state.
    Contact(fs_contact::normal_patch::NormalPatchError),
    /// The dry-friction owner refused a pad state.
    Friction(fs_tribo::TriboError),
    /// A point-contact request returned an impossible line receipt.
    UnexpectedContactReceipt,
    /// The completed rollout produced a non-finite score or receipt.
    NonFiniteObjective,
}

impl fmt::Display for ManipulationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field } => {
                write!(formatter, "invalid manipulation config: {field}")
            }
            Self::ParameterCount { expected, actual } => {
                write!(
                    formatter,
                    "manipulation policy needs {expected} coordinates, received {actual}"
                )
            }
            Self::NonFiniteParameter { index } => {
                write!(
                    formatter,
                    "manipulation policy coordinate {index} is not finite"
                )
            }
            Self::Robot(error) => write!(formatter, "iiwa articulated owner refused: {error}"),
            Self::Geometry(error) => write!(formatter, "iiwa Lie owner refused: {error}"),
            Self::Contact(error) => write!(formatter, "finger contact owner refused: {error}"),
            Self::Friction(error) => write!(formatter, "finger friction owner refused: {error}"),
            Self::UnexpectedContactReceipt => {
                formatter.write_str("finger sphere/plane request returned a line receipt")
            }
            Self::NonFiniteObjective => {
                formatter.write_str("manipulation rollout produced a non-finite receipt")
            }
        }
    }
}

impl std::error::Error for ManipulationError {}

impl From<fs_mbd::articulated::ArticulatedError> for ManipulationError {
    fn from(value: fs_mbd::articulated::ArticulatedError) -> Self {
        Self::Robot(value)
    }
}

impl From<fs_ga::GaError> for ManipulationError {
    fn from(value: fs_ga::GaError) -> Self {
        Self::Geometry(value)
    }
}

impl From<fs_contact::normal_patch::NormalPatchError> for ManipulationError {
    fn from(value: fs_contact::normal_patch::NormalPatchError) -> Self {
        Self::Contact(value)
    }
}

impl From<fs_tribo::TriboError> for ManipulationError {
    fn from(value: fs_tribo::TriboError) -> Self {
        Self::Friction(value)
    }
}

#[derive(Debug, Clone, Copy)]
struct TaskDefinition {
    initial_yaw: f64,
    goal_yaw: f64,
    mass_kg: f64,
    dimensions_m: Vec3,
    grasp_half_width_m: f64,
    obstacle_center_m: Vec3,
    obstacle_half_extents_m: Vec3,
}

impl TaskDefinition {
    fn for_task(task: ManipulationTask) -> Self {
        match task {
            ManipulationTask::KitchenMug => Self {
                initial_yaw: -0.62,
                goal_yaw: 0.72,
                mass_kg: 0.34,
                dimensions_m: Vec3::new(0.092, 0.122, 0.104),
                grasp_half_width_m: 0.046,
                obstacle_center_m: Vec3::new(0.28, 0.02, 0.66),
                obstacle_half_extents_m: Vec3::new(0.10, 0.14, 0.20),
            },
            ManipulationTask::LivingRoomRemote => Self {
                initial_yaw: 0.48,
                goal_yaw: -0.78,
                mass_kg: 0.18,
                dimensions_m: Vec3::new(0.058, 0.188, 0.024),
                grasp_half_width_m: 0.029,
                obstacle_center_m: Vec3::new(0.22, -0.20, 0.62),
                obstacle_half_extents_m: Vec3::new(0.16, 0.09, 0.16),
            },
            ManipulationTask::BackyardTrowel => Self {
                initial_yaw: -0.82,
                goal_yaw: 0.88,
                mass_kg: 0.24,
                dimensions_m: Vec3::new(0.086, 0.302, 0.048),
                grasp_half_width_m: 0.021,
                obstacle_center_m: Vec3::new(-0.08, 0.24, 0.64),
                obstacle_half_extents_m: Vec3::new(0.13, 0.10, 0.18),
            },
        }
    }
}

/// Reusable evaluator that constructs the source model and constitutive owners once.
#[derive(Debug, Clone)]
pub struct ManipulationEvaluator {
    config: ManipulationConfig,
    catalog: CatalogRobotModel,
    scene: ManipulationScene,
    initial_object_pose: Se3,
    interface: InterfaceSystemRef,
    friction: FrictionLaw,
    contact_frame: ContactFrame,
    step_count: usize,
}

impl ManipulationEvaluator {
    /// Admit a complete fixed benchmark configuration.
    pub fn new(config: ManipulationConfig) -> Result<Self, ManipulationError> {
        validate_config(config)?;
        let catalog = kuka_lbr_iiwa7_r800()?;
        let definition = TaskDefinition::for_task(config.task);
        let initial_configuration = grasp_configuration(definition.initial_yaw);
        let goal_configuration = grasp_configuration(definition.goal_yaw);
        let initial_object_pose = end_effector_pose(&catalog, &initial_configuration)?;
        let goal_object_pose = end_effector_pose(&catalog, &goal_configuration)?;
        let initial_object_position_m = initial_object_pose.translation();
        let goal_object_position_m = goal_object_pose.translation();
        let support_height_m = initial_object_position_m.z - 0.5 * definition.dimensions_m.z;
        if !support_height_m.is_finite() || support_height_m <= 0.0 {
            return Err(ManipulationError::InvalidConfig {
                field: "derived support height",
            });
        }
        if (goal_object_position_m.z - initial_object_position_m.z).abs() > 1.0e-9 {
            return Err(ManipulationError::InvalidConfig {
                field: "source-derived grasp/place support heights",
            });
        }
        let scene = ManipulationScene {
            object_mass_kg: definition.mass_kg,
            object_dimensions_m: definition.dimensions_m,
            grasp_half_width_m: definition.grasp_half_width_m,
            initial_object_position_m,
            goal_object_position_m,
            support_height_m,
            obstacle_center_m: definition.obstacle_center_m,
            obstacle_half_extents_m: definition.obstacle_half_extents_m,
        };
        let interface = InterfaceSystemRef::new(
            "iiwa-polyurethane-pad--household-object",
            "iiwa-manipulation-rollout-v1",
            "caller-declared-browser-demo-interface",
            InputAuthority::Estimated,
            InterfaceMedium::Dry,
        )?;
        let friction = FrictionLaw::Coulomb {
            static_mu: 0.58,
            kinetic_mu: 0.46,
        };
        let contact_frame = ContactFrame::new([1.0, 0.0, 0.0])?;
        let step_count = rounded_step_count(config)?;
        Ok(Self {
            config,
            catalog,
            scene,
            initial_object_pose,
            interface,
            friction,
            contact_frame,
            step_count,
        })
    }

    /// Fixed admitted controls.
    #[must_use]
    pub const fn config(&self) -> ManipulationConfig {
        self.config
    }

    /// Source-derived scene geometry and declared object estimates.
    #[must_use]
    pub const fn scene(&self) -> ManipulationScene {
        self.scene
    }

    /// Evaluate one policy without retaining render poses.
    pub fn evaluate(&self, parameters: &[f64]) -> Result<ManipulationReceipt, ManipulationError> {
        self.rollout(parameters, false)
    }

    /// Evaluate one policy and retain decimated owner poses.
    pub fn trace(&self, parameters: &[f64]) -> Result<ManipulationReceipt, ManipulationError> {
        self.rollout(parameters, true)
    }

    fn rollout(
        &self,
        parameters: &[f64],
        retain_trace: bool,
    ) -> Result<ManipulationReceipt, ManipulationError> {
        validate_policy(parameters)?;
        let base = BaseState::stationary(Se3::identity());
        let mut joint_position = policy_joint_knots(parameters, 0);
        let mut joint_velocity = [0.0; ARM_JOINT_COUNT];
        project_initial_configuration(&self.catalog, &mut joint_position);
        let mut gripper_width_m = OPEN_GRIPPER_WIDTH_M;
        let mut object_pose = self.initial_object_pose;
        let mut object_velocity_world = Vec3::new(0.0, 0.0, 0.0);
        let mut grasped = false;
        let mut ever_grasped = false;
        let mut released_after_transport = false;
        let mut first_grasp_time_s = self.config.duration_s;
        let mut grasp_duration_s = 0.0;
        let mut peak_grip_force_n = 0.0_f64;
        let mut minimum_reach_error_m = f64::INFINITY;
        let mut maximum_lift_m = 0.0_f64;
        let mut actuator_work_j = 0.0_f64;
        let mut obstacle_integral = 0.0_f64;
        let mut control_limit_integral = 0.0_f64;
        let trace_capacity = if retain_trace {
            self.step_count / self.config.trace_stride + 2
        } else {
            0
        };
        let mut trace = Vec::with_capacity(trace_capacity);
        let mut previous_tool_position =
            end_effector_pose(&self.catalog, &joint_position)?.translation();

        for step in 0..self.step_count {
            let time_s = step as f64 * self.config.step_s;
            let progress = time_s / self.config.duration_s;
            let kinematics =
                forward_kinematics(self.catalog.model(), base, &joint_position, &joint_velocity)?;
            let tool_pose = kinematics.world_from_link[END_EFFECTOR_LINK];
            let tool_position = tool_pose.translation();
            let reach_error = vec_norm(object_pose.translation() - tool_position);
            if !ever_grasped {
                minimum_reach_error_m = minimum_reach_error_m.min(reach_error);
            }

            let (desired_position, desired_velocity, desired_gripper_width, proposed_excess) =
                desired_controls(
                    &self.catalog,
                    parameters,
                    progress,
                    self.config.step_s / self.config.duration_s,
                    self.config.step_s,
                );
            control_limit_integral += proposed_excess * self.config.step_s;
            let width_delta = (desired_gripper_width - gripper_width_m).clamp(
                -MAX_GRIPPER_SPEED_M_PER_S * self.config.step_s,
                MAX_GRIPPER_SPEED_M_PER_S * self.config.step_s,
            );
            gripper_width_m += width_delta;

            let mut external = [Wrench::default(); ARM_LINK_COUNT];
            if grasped {
                let force_world = GRAVITY_WORLD_M_PER_S2.scale(self.scene.object_mass_kg);
                let force_body = tool_pose.rotation().inverse().rotate(force_world)?;
                let application_body = tool_pose
                    .inverse()?
                    .transform_point(object_pose.translation())?;
                external[END_EFFECTOR_LINK] =
                    Wrench::new(application_body.cross(force_body), force_body);
            }
            // Use the articulated owner's inverse dynamics as a computed-torque
            // tracker. Feeding a bounded desired acceleration through the exact
            // source inertia/coriolis/gravity model avoids the joint-dependent
            // instability of applying one set of raw torque-space PD gains to
            // links whose reflected inertias differ by orders of magnitude.
            let mut desired_acceleration = [0.0; ARM_JOINT_COUNT];
            for joint in 0..ARM_JOINT_COUNT {
                desired_acceleration[joint] = (72.0
                    * (desired_position[joint] - joint_position[joint])
                    + 17.0 * (desired_velocity[joint] - joint_velocity[joint]))
                    .clamp(-24.0, 24.0);
            }
            let tracking_force = inverse_dynamics(
                self.catalog.model(),
                base,
                &joint_position,
                &joint_velocity,
                &desired_acceleration,
                GRAVITY_WORLD_M_PER_S2,
                &external,
            )?;
            let mut generalized_force = [0.0; ARM_JOINT_COUNT];
            for joint in 0..ARM_JOINT_COUNT {
                let metadata = self.catalog.joints()[joint];
                generalized_force[joint] = tracking_force.generalized_force[joint].clamp(
                    -metadata.effort_newton_metres,
                    metadata.effort_newton_metres,
                );
                actuator_work_j +=
                    (generalized_force[joint] * joint_velocity[joint]).abs() * self.config.step_s;
            }
            let dynamics = forward_dynamics(
                self.catalog.model(),
                base,
                &joint_position,
                &joint_velocity,
                &generalized_force,
                GRAVITY_WORLD_M_PER_S2,
                &external,
            )?;
            for joint in 0..ARM_JOINT_COUNT {
                let metadata = self.catalog.joints()[joint];
                joint_velocity[joint] +=
                    dynamics.generalized_acceleration[joint] * self.config.step_s;
                joint_velocity[joint] = joint_velocity[joint].clamp(
                    -metadata.velocity_rad_per_second,
                    metadata.velocity_rad_per_second,
                );
                let proposed = joint_position[joint] + joint_velocity[joint] * self.config.step_s;
                let projected =
                    proposed.clamp(metadata.lower_position_rad, metadata.upper_position_rad);
                control_limit_integral += (proposed - projected).abs();
                if projected != proposed {
                    joint_velocity[joint] = 0.0;
                }
                joint_position[joint] = projected;
            }

            let next_kinematics =
                forward_kinematics(self.catalog.model(), base, &joint_position, &joint_velocity)?;
            let next_tool_pose = next_kinematics.world_from_link[END_EFFECTOR_LINK];
            let next_tool_position = next_tool_pose.translation();
            let tool_velocity_world =
                (next_tool_position - previous_tool_position).scale(1.0 / self.config.step_s);
            previous_tool_position = next_tool_position;

            let grip = self.grip_state(
                next_tool_pose,
                object_pose,
                gripper_width_m,
                width_delta / self.config.step_s,
            )?;
            peak_grip_force_n = peak_grip_force_n.max(grip.normal_force_n);
            if !grasped && grip.can_capture {
                grasped = true;
                ever_grasped = true;
                first_grasp_time_s = first_grasp_time_s.min(time_s + self.config.step_s);
                object_pose = next_tool_pose;
                object_velocity_world = tool_velocity_world;
            } else if grasped && !grip.can_hold {
                grasped = false;
                released_after_transport |= maximum_lift_m >= 0.5 * LIFT_TARGET_M;
            }
            if grasped {
                grasp_duration_s += self.config.step_s;
                object_pose = next_tool_pose;
                object_velocity_world = tool_velocity_world;
            } else {
                object_velocity_world =
                    object_velocity_world + GRAVITY_WORLD_M_PER_S2.scale(self.config.step_s);
                let mut position =
                    object_pose.translation() + object_velocity_world.scale(self.config.step_s);
                let supported_z =
                    self.scene.support_height_m + 0.5 * self.scene.object_dimensions_m.z;
                if position.z < supported_z {
                    position.z = supported_z;
                    if object_velocity_world.z < 0.0 {
                        object_velocity_world.z = 0.0;
                    }
                    object_velocity_world.x *= 0.82;
                    object_velocity_world.y *= 0.82;
                }
                object_pose = Se3::from_parts(object_pose.rotation(), position)?;
            }
            maximum_lift_m = maximum_lift_m
                .max(object_pose.translation().z - self.scene.initial_object_position_m.z);
            obstacle_integral += obstacle_proximity(
                &next_kinematics.world_from_link,
                self.scene.obstacle_center_m,
                self.scene.obstacle_half_extents_m,
            ) * self.config.step_s;

            if retain_trace && step.is_multiple_of(self.config.trace_stride) {
                trace.push(trace_sample(
                    time_s + self.config.step_s,
                    &next_kinematics.world_from_link,
                    object_pose,
                    gripper_width_m,
                    grip.normal_force_n,
                    grasped,
                ));
            }
        }

        if retain_trace {
            let final_time = self.step_count as f64 * self.config.step_s;
            let needs_terminal = trace
                .last()
                .is_none_or(|sample| (sample.time_s - final_time).abs() > 1.0e-12);
            if needs_terminal {
                let kinematics = forward_kinematics(
                    self.catalog.model(),
                    base,
                    &joint_position,
                    &joint_velocity,
                )?;
                trace.push(trace_sample(
                    final_time,
                    &kinematics.world_from_link,
                    object_pose,
                    gripper_width_m,
                    0.0,
                    grasped,
                ));
            }
        }

        let final_object_error_m =
            vec_norm(object_pose.translation() - self.scene.goal_object_position_m);
        let placed = ever_grasped
            && released_after_transport
            && maximum_lift_m >= LIFT_TARGET_M
            && final_object_error_m <= PLACEMENT_TOLERANCE_M
            && !grasped;
        let grasp_penalty = if ever_grasped { 0.0 } else { 220.0 };
        let release_penalty = if released_after_transport { 0.0 } else { 100.0 };
        let lift_penalty = 160.0 * (1.0 - maximum_lift_m / LIFT_TARGET_M).clamp(0.0, 1.0);
        let placement_bonus = if placed { 180.0 } else { 0.0 };
        let objective = 320.0 * final_object_error_m.min(2.0)
            + 80.0 * minimum_reach_error_m.min(1.0)
            + lift_penalty
            + grasp_penalty
            + release_penalty
            + 45.0 * obstacle_integral
            + 18.0 * control_limit_integral
            + 0.0015 * actuator_work_j
            - placement_bonus;
        let receipt = ManipulationReceipt {
            objective,
            final_object_error_m,
            minimum_reach_error_m,
            maximum_lift_m,
            actuator_work_j,
            obstacle_integral,
            control_limit_integral,
            first_grasp_time_s,
            grasp_duration_s,
            peak_grip_force_n,
            ever_grasped,
            released_after_transport,
            placed,
            completed_steps: self.step_count,
            trace,
        };
        validate_receipt(&receipt)?;
        Ok(receipt)
    }

    fn grip_state(
        &self,
        tool_pose: Se3,
        object_pose: Se3,
        gripper_width_m: f64,
        gripper_speed_m_per_s: f64,
    ) -> Result<GripState, ManipulationError> {
        let tool_from_object = tool_pose.inverse()?.compose(object_pose)?;
        let relative = tool_from_object.translation();
        let angular_error = vec_norm(tool_from_object.log().angular);
        let centering_error = vec_norm(relative);
        let indentation_m = (self.scene.grasp_half_width_m - 0.5 * gripper_width_m)
            .clamp(0.0, MAX_PAD_INDENTATION_M);
        let aligned = centering_error <= 0.045 && angular_error <= 0.45;
        let normal_force_n = if aligned && indentation_m > 0.0 {
            2.0 * self.pad_force(indentation_m, -0.5 * gripper_speed_m_per_s)?
        } else {
            0.0
        };
        let static_capacity_n = if normal_force_n > 0.0 {
            let response = self.friction.evaluate(
                &self.interface,
                normal_force_n,
                TangentialSlip::new(&self.contact_frame, [0.0, 0.0, 0.0])?,
            )?;
            response.static_limit
        } else {
            0.0
        };
        let required_capacity_n = self.scene.object_mass_kg * GRAVITY_WORLD_M_PER_S2.z.abs() * 1.15;
        Ok(GripState {
            normal_force_n,
            can_capture: aligned
                && gripper_speed_m_per_s < -1.0e-5
                && static_capacity_n >= required_capacity_n,
            can_hold: aligned
                && gripper_width_m < OPEN_GRIPPER_WIDTH_M - 0.01
                && static_capacity_n >= required_capacity_n,
        })
    }

    fn pad_force(
        &self,
        indentation_m: f64,
        indentation_rate_m_per_s: f64,
    ) -> Result<f64, ManipulationError> {
        let mut request = pad_request(self.interface.clone(), self.config.step_s);
        request.indentation_m = indentation_m;
        request.indentation_rate_m_per_s = indentation_rate_m_per_s;
        match request.evaluate()? {
            NormalPatchReceipt::Point(receipt) => Ok(receipt.normal_force_n),
            NormalPatchReceipt::Line(_) => Err(ManipulationError::UnexpectedContactReceipt),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct GripState {
    normal_force_n: f64,
    can_capture: bool,
    can_hold: bool,
}

/// Disclosed source-feasible curriculum mean in the exact 128-coordinate layout.
#[must_use]
pub fn manipulation_curriculum_mean(task: ManipulationTask) -> [f64; ARM_POLICY_DIMENSION] {
    let definition = TaskDefinition::for_task(task);
    let mut mean = [0.0; ARM_POLICY_DIMENSION];
    for knot in 0..ARM_POLICY_KNOTS {
        let progress = knot as f64 / (ARM_POLICY_KNOTS - 1) as f64;
        let configuration = reference_configuration(definition, progress);
        for joint in 0..ARM_JOINT_COUNT {
            mean[joint * ARM_POLICY_KNOTS + knot] = configuration[joint];
        }
        mean[ARM_JOINT_COUNT * ARM_POLICY_KNOTS + knot] =
            reference_gripper_width(definition, progress);
    }
    mean
}

fn validate_config(config: ManipulationConfig) -> Result<(), ManipulationError> {
    if !config.step_s.is_finite() || !(1.0 / 240.0..=1.0 / 45.0).contains(&config.step_s) {
        return Err(ManipulationError::InvalidConfig { field: "step_s" });
    }
    if !config.duration_s.is_finite()
        || !(3.0..=6.0).contains(&config.duration_s)
        || config.duration_s < config.step_s
    {
        return Err(ManipulationError::InvalidConfig {
            field: "duration_s",
        });
    }
    if !(1..=1_000).contains(&config.trace_stride) {
        return Err(ManipulationError::InvalidConfig {
            field: "trace_stride",
        });
    }
    Ok(())
}

fn rounded_step_count(config: ManipulationConfig) -> Result<usize, ManipulationError> {
    let steps = (config.duration_s / config.step_s).round();
    if !steps.is_finite() || !(1.0..=1_440.0).contains(&steps) {
        return Err(ManipulationError::InvalidConfig {
            field: "rounded step count",
        });
    }
    Ok(steps as usize)
}

fn validate_policy(parameters: &[f64]) -> Result<(), ManipulationError> {
    if parameters.len() != ARM_POLICY_DIMENSION {
        return Err(ManipulationError::ParameterCount {
            expected: ARM_POLICY_DIMENSION,
            actual: parameters.len(),
        });
    }
    for (index, value) in parameters.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(ManipulationError::NonFiniteParameter { index });
        }
    }
    Ok(())
}

fn project_initial_configuration(
    catalog: &CatalogRobotModel,
    position: &mut [f64; ARM_JOINT_COUNT],
) {
    for (value, metadata) in position.iter_mut().zip(catalog.joints()) {
        *value = value.clamp(metadata.lower_position_rad, metadata.upper_position_rad);
    }
}

fn policy_joint_knots(parameters: &[f64], knot: usize) -> [f64; ARM_JOINT_COUNT] {
    core::array::from_fn(|joint| parameters[joint * ARM_POLICY_KNOTS + knot])
}

fn desired_controls(
    catalog: &CatalogRobotModel,
    parameters: &[f64],
    progress: f64,
    progress_step: f64,
    step_s: f64,
) -> ([f64; ARM_JOINT_COUNT], [f64; ARM_JOINT_COUNT], f64, f64) {
    let (raw_position, raw_width) = interpolate_policy(parameters, progress);
    let (next_position, _) = interpolate_policy(parameters, (progress + progress_step).min(1.0));
    let mut desired_position = [0.0; ARM_JOINT_COUNT];
    let mut desired_velocity = [0.0; ARM_JOINT_COUNT];
    let mut excess = 0.0;
    for joint in 0..ARM_JOINT_COUNT {
        let metadata = catalog.joints()[joint];
        desired_position[joint] =
            raw_position[joint].clamp(metadata.lower_position_rad, metadata.upper_position_rad);
        let next =
            next_position[joint].clamp(metadata.lower_position_rad, metadata.upper_position_rad);
        desired_velocity[joint] = ((next - desired_position[joint]) / step_s).clamp(
            -metadata.velocity_rad_per_second,
            metadata.velocity_rad_per_second,
        );
        excess += (raw_position[joint] - desired_position[joint]).abs();
    }
    let desired_width = raw_width.clamp(MIN_GRIPPER_WIDTH_M, OPEN_GRIPPER_WIDTH_M);
    excess += 20.0 * (raw_width - desired_width).abs();
    (desired_position, desired_velocity, desired_width, excess)
}

fn interpolate_policy(parameters: &[f64], progress: f64) -> ([f64; ARM_JOINT_COUNT], f64) {
    let coordinate = progress.clamp(0.0, 1.0) * (ARM_POLICY_KNOTS - 1) as f64;
    let left = (coordinate.floor() as usize).min(ARM_POLICY_KNOTS - 1);
    let right = (left + 1).min(ARM_POLICY_KNOTS - 1);
    let blend = smoothstep(coordinate - left as f64);
    let mut joints = [0.0; ARM_JOINT_COUNT];
    for (joint, value) in joints.iter_mut().enumerate() {
        let offset = joint * ARM_POLICY_KNOTS;
        *value = lerp(parameters[offset + left], parameters[offset + right], blend);
    }
    let offset = ARM_JOINT_COUNT * ARM_POLICY_KNOTS;
    let width = lerp(parameters[offset + left], parameters[offset + right], blend);
    (joints, width)
}

fn reference_configuration(definition: TaskDefinition, progress: f64) -> [f64; ARM_JOINT_COUNT] {
    let home = home_configuration(0.5 * (definition.initial_yaw + definition.goal_yaw));
    let pregrasp = pregrasp_configuration(definition.initial_yaw);
    let grasp = grasp_configuration(definition.initial_yaw);
    let lift_initial = lift_configuration(definition.initial_yaw);
    let lift_goal = lift_configuration(definition.goal_yaw);
    let place = grasp_configuration(definition.goal_yaw);
    if progress < 0.20 {
        blend_configuration(home, pregrasp, progress / 0.20)
    } else if progress < 0.34 {
        blend_configuration(pregrasp, grasp, (progress - 0.20) / 0.14)
    } else if progress < 0.47 {
        grasp
    } else if progress < 0.61 {
        blend_configuration(grasp, lift_initial, (progress - 0.47) / 0.14)
    } else if progress < 0.76 {
        blend_configuration(lift_initial, lift_goal, (progress - 0.61) / 0.15)
    } else if progress < 0.88 {
        blend_configuration(lift_goal, place, (progress - 0.76) / 0.12)
    } else if progress < 0.95 {
        place
    } else {
        blend_configuration(place, home, (progress - 0.95) / 0.05)
    }
}

fn reference_gripper_width(definition: TaskDefinition, progress: f64) -> f64 {
    let closed = (2.0 * definition.grasp_half_width_m - 0.002)
        .clamp(MIN_GRIPPER_WIDTH_M, OPEN_GRIPPER_WIDTH_M);
    if progress < 0.34 {
        OPEN_GRIPPER_WIDTH_M
    } else if progress < 0.46 {
        lerp(
            OPEN_GRIPPER_WIDTH_M,
            closed,
            smoothstep((progress - 0.34) / 0.12),
        )
    } else if progress < 0.88 {
        closed
    } else if progress < 0.94 {
        lerp(
            closed,
            OPEN_GRIPPER_WIDTH_M,
            smoothstep((progress - 0.88) / 0.06),
        )
    } else {
        OPEN_GRIPPER_WIDTH_M
    }
}

const fn home_configuration(yaw: f64) -> [f64; ARM_JOINT_COUNT] {
    [yaw, -0.42, 0.0, 0.82, 0.0, -0.40, 0.0]
}

const fn pregrasp_configuration(yaw: f64) -> [f64; ARM_JOINT_COUNT] {
    [yaw, -0.64, 0.0, 1.28, 0.0, -0.64, 0.0]
}

const fn grasp_configuration(yaw: f64) -> [f64; ARM_JOINT_COUNT] {
    [yaw, -0.76, 0.0, 1.52, 0.0, -0.76, 0.0]
}

const fn lift_configuration(yaw: f64) -> [f64; ARM_JOINT_COUNT] {
    [yaw, -0.48, 0.0, 0.96, 0.0, -0.48, 0.0]
}

fn blend_configuration(
    start: [f64; ARM_JOINT_COUNT],
    end: [f64; ARM_JOINT_COUNT],
    amount: f64,
) -> [f64; ARM_JOINT_COUNT] {
    let amount = smoothstep(amount.clamp(0.0, 1.0));
    core::array::from_fn(|joint| lerp(start[joint], end[joint], amount))
}

fn end_effector_pose(
    catalog: &CatalogRobotModel,
    position: &[f64; ARM_JOINT_COUNT],
) -> Result<Se3, ManipulationError> {
    let kinematics = forward_kinematics(
        catalog.model(),
        BaseState::stationary(Se3::identity()),
        position,
        &[0.0; ARM_JOINT_COUNT],
    )?;
    Ok(kinematics.world_from_link[END_EFFECTOR_LINK])
}

fn pad_request(interface: InterfaceSystemRef, step_s: f64) -> NormalPatchRequest {
    NormalPatchRequest {
        identity: NormalPatchIdentity {
            model_id: MANIPULATION_MODEL_ID.to_owned(),
            source_id: "estimated-polyurethane-finger-pad-sphere".to_owned(),
            state_id: "iiwa-household-grasp-state".to_owned(),
        },
        interface,
        law: NormalPatchLaw::HuntCrossleySphere {
            effective_radius_m: PAD_EFFECTIVE_RADIUS_M,
            reduced_modulus_pa: PAD_REDUCED_MODULUS_PA,
            dissipation_s_per_m: 0.20,
        },
        geometry: NormalPatchGeometry::SpherePlane,
        indentation_m: 0.0,
        indentation_rate_m_per_s: 0.0,
        step_s,
        line_load_n_per_m: 0.0,
        applicability: ApplicabilityInput {
            half_space_depth_m: 0.05,
            layer_thickness_m: 0.025,
            yield_strength_pa: 18.0e6,
            characteristic_rate_m_per_s: 2.0,
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
            radius_relative: 0.25,
            modulus_relative: 0.50,
            load_relative: 0.0,
        },
    }
}

fn obstacle_proximity(poses: &[Se3], center: Vec3, half_extents: Vec3) -> f64 {
    let margin = 0.075;
    let mut total = 0.0;
    for pose in poses.iter().skip(1) {
        let point = pose.translation();
        let dx = (point.x - center.x).abs() - half_extents.x;
        let dy = (point.y - center.y).abs() - half_extents.y;
        let dz = (point.z - center.z).abs() - half_extents.z;
        let outside = Vec3::new(dx.max(0.0), dy.max(0.0), dz.max(0.0));
        let distance = vec_norm(outside);
        if dx <= 0.0 && dy <= 0.0 && dz <= 0.0 {
            total += margin + (-dx).min((-dy).min(-dz));
        } else if distance < margin {
            total += margin - distance;
        }
    }
    total
}

fn trace_sample(
    time_s: f64,
    poses: &[Se3],
    object_pose: Se3,
    gripper_width_m: f64,
    grip_normal_force_n: f64,
    grasped: bool,
) -> ManipulationTraceSample {
    let mut link_pose = [[0.0; ARM_LINK_POSE_WORDS]; ARM_LINK_COUNT];
    for (output, pose) in link_pose.iter_mut().zip(poses) {
        *output = pose_words(*pose);
    }
    ManipulationTraceSample {
        time_s,
        gripper_width_m,
        grip_normal_force_n,
        grasped,
        object_pose: pose_words(object_pose),
        link_pose,
    }
}

fn pose_words(pose: Se3) -> [f64; ARM_LINK_POSE_WORDS] {
    let translation = pose.translation();
    let rotation = pose.rotation();
    let quaternion = rotation.as_quat();
    [
        translation.x,
        translation.y,
        translation.z,
        quaternion.w,
        quaternion.x,
        quaternion.y,
        quaternion.z,
    ]
}

fn validate_receipt(receipt: &ManipulationReceipt) -> Result<(), ManipulationError> {
    let scalars = [
        receipt.objective,
        receipt.final_object_error_m,
        receipt.minimum_reach_error_m,
        receipt.maximum_lift_m,
        receipt.actuator_work_j,
        receipt.obstacle_integral,
        receipt.control_limit_integral,
        receipt.first_grasp_time_s,
        receipt.grasp_duration_s,
        receipt.peak_grip_force_n,
    ];
    if scalars.iter().any(|value| !value.is_finite())
        || receipt.final_object_error_m < 0.0
        || receipt.minimum_reach_error_m < 0.0
        || receipt.maximum_lift_m < 0.0
        || receipt.actuator_work_j < 0.0
        || receipt.obstacle_integral < 0.0
        || receipt.control_limit_integral < 0.0
        || receipt.first_grasp_time_s < 0.0
        || receipt.grasp_duration_s < 0.0
        || receipt.peak_grip_force_n < 0.0
    {
        return Err(ManipulationError::NonFiniteObjective);
    }
    Ok(())
}

fn vec_norm(vector: Vec3) -> f64 {
    vector.dot(vector).sqrt()
}

fn smoothstep(value: f64) -> f64 {
    let value = value.clamp(0.0, 1.0);
    value * value * (3.0 - 2.0 * value)
}

fn lerp(start: f64, end: f64, amount: f64) -> f64 {
    start + amount * (end - start)
}

/// Maximum population accepted by the packed batch boundary.
pub const fn manipulation_max_population() -> usize {
    MAX_POLICY_POPULATION
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_feasible_waypoints_define_level_pick_and_place_stations()
    -> Result<(), ManipulationError> {
        for task in [
            ManipulationTask::KitchenMug,
            ManipulationTask::LivingRoomRemote,
            ManipulationTask::BackyardTrowel,
        ] {
            let evaluator = ManipulationEvaluator::new(ManipulationConfig {
                task,
                ..ManipulationConfig::default()
            })?;
            assert_eq!(evaluator.catalog.model().dof_count(), ARM_JOINT_COUNT);
            assert_eq!(evaluator.catalog.model().link_count(), ARM_LINK_COUNT);
            assert!(
                (evaluator.scene.initial_object_position_m.z
                    - evaluator.scene.goal_object_position_m.z)
                    .abs()
                    < 1.0e-10
            );
            let definition = TaskDefinition::for_task(task);
            let lifted = end_effector_pose(
                &evaluator.catalog,
                &lift_configuration(definition.initial_yaw),
            )?
            .translation();
            assert!(
                lifted.z - evaluator.scene.initial_object_position_m.z > LIFT_TARGET_M,
                "task={task:?}, initial={:?}, lifted={lifted:?}",
                evaluator.scene.initial_object_position_m
            );
        }
        Ok(())
    }

    #[test]
    fn curriculum_rollouts_are_deterministic_and_complete_all_three_tasks()
    -> Result<(), ManipulationError> {
        for task in [
            ManipulationTask::KitchenMug,
            ManipulationTask::LivingRoomRemote,
            ManipulationTask::BackyardTrowel,
        ] {
            let evaluator = ManipulationEvaluator::new(ManipulationConfig {
                task,
                ..ManipulationConfig::default()
            })?;
            let mean = manipulation_curriculum_mean(task);
            let first = evaluator.trace(&mean)?;
            let second = evaluator.trace(&mean)?;
            assert_eq!(first, second);
            assert_eq!(first.completed_steps, evaluator.step_count);
            assert!(first.trace.len() >= 50);
            let closest = first
                .trace
                .iter()
                .min_by(|left, right| {
                    let distance = |sample: &ManipulationTraceSample| {
                        let object = Vec3::new(
                            sample.object_pose[0],
                            sample.object_pose[1],
                            sample.object_pose[2],
                        );
                        let tool = Vec3::new(
                            sample.link_pose[END_EFFECTOR_LINK][0],
                            sample.link_pose[END_EFFECTOR_LINK][1],
                            sample.link_pose[END_EFFECTOR_LINK][2],
                        );
                        vec_norm(object - tool)
                    };
                    distance(left).total_cmp(&distance(right))
                })
                .expect("retained trace is non-empty");
            let snapshots = first
                .trace
                .iter()
                .step_by((first.trace.len() / 6).max(1))
                .map(|sample| {
                    (
                        sample.time_s,
                        &sample.link_pose[END_EFFECTOR_LINK][0..3],
                        sample.gripper_width_m,
                    )
                })
                .collect::<Vec<_>>();
            let summary = format!(
                "task={task:?}, objective={:.6}, final_error={:.6}, min_reach={:.6}, \
                 max_lift={:.6}, first_grasp={:.6}, grasp_duration={:.6}, peak_force={:.6}, \
                 ever_grasped={}, released={}, placed={}, closest_time={:.6}, \
                 closest_tool={:?}, closest_object={:?}, snapshots={snapshots:?}",
                first.objective,
                first.final_object_error_m,
                first.minimum_reach_error_m,
                first.maximum_lift_m,
                first.first_grasp_time_s,
                first.grasp_duration_s,
                first.peak_grip_force_n,
                first.ever_grasped,
                first.released_after_transport,
                first.placed,
                closest.time_s,
                &closest.link_pose[END_EFFECTOR_LINK][0..3],
                &closest.object_pose[0..3],
            );
            assert!(first.ever_grasped, "{summary}");
            assert!(first.maximum_lift_m >= LIFT_TARGET_M, "{summary}");
            assert!(first.released_after_transport, "{summary}");
            assert!(first.placed, "{summary}");
        }
        Ok(())
    }

    #[test]
    fn policy_validation_and_objective_are_candidate_sensitive() -> Result<(), ManipulationError> {
        let evaluator = ManipulationEvaluator::new(ManipulationConfig::default())?;
        let mean = manipulation_curriculum_mean(ManipulationTask::KitchenMug);
        let baseline = evaluator.evaluate(&mean)?;
        let mut changed = mean;
        changed[2 * ARM_POLICY_KNOTS + 5] += 0.18;
        let perturbed = evaluator.evaluate(&changed)?;
        assert_ne!(baseline.objective, perturbed.objective);
        let mut invalid = mean;
        invalid[17] = f64::NAN;
        assert!(matches!(
            evaluator.evaluate(&invalid),
            Err(ManipulationError::NonFiniteParameter { index: 17 })
        ));
        assert!(matches!(
            evaluator.evaluate(&mean[..ARM_POLICY_DIMENSION - 1]),
            Err(ManipulationError::ParameterCount { .. })
        ));
        Ok(())
    }
}
