//! Provenance-bound articulated robot model catalog.
//!
//! The builders in this module transcribe selected URDF/Xacro records into the
//! canonical [`crate::articulated`] and `fs-ga` types. They do not parse URDF at
//! runtime, load meshes, or introduce a second robot-math representation.

use crate::articulated::{
    ArticulatedError, ArticulatedModel, JointLimits, JointModel, Link, SpatialInertia,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3};
use fs_math::det;

/// Version of the typed in-source robot catalog layout.
pub const ROBOT_MODEL_CATALOG_VERSION: u32 = 1;

/// Number of actuators retained by the G1 lower-body-and-waist catalog.
pub const G1_POLICY_ACTUATORS: usize = 15;
/// Raw state and command signals supplied to the G1 residual policy.
pub const G1_POLICY_RAW_SIGNALS: usize = 42;
/// Periodic phase functions multiplying every raw signal.
pub const G1_POLICY_PHASE_BASIS: usize = 8;
/// Features seen by each G1 actuator.
pub const G1_POLICY_FEATURES_PER_ACTUATOR: usize = G1_POLICY_RAW_SIGNALS * G1_POLICY_PHASE_BASIS;
/// Exact search-space dimension of the flagship G1 residual policy.
pub const G1_POLICY_DIMENSION: usize = G1_POLICY_ACTUATORS * G1_POLICY_FEATURES_PER_ACTUATOR;

/// One complete observation for the catalog-owned 5,040-parameter G1 policy.
///
/// The signal ordering is deliberately closed and documented: bias, 15 joint
/// positions, 15 joint velocities, body-frame gravity direction, body angular
/// velocity, body-frame target-velocity error, and left/right contact state.
/// A browser or optimizer may choose the parameters, but it must not silently
/// invent a different feature map under the same dimensional label.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct G1PolicyObservation {
    /// Joint angles in compact catalog order [rad].
    pub joint_position_rad: [f64; G1_POLICY_ACTUATORS],
    /// Joint rates in compact catalog order [rad/s].
    pub joint_velocity_rad_per_s: [f64; G1_POLICY_ACTUATORS],
    /// Unit gravity direction expressed in the pelvis frame.
    pub gravity_direction_body: Vec3,
    /// Pelvis angular velocity expressed in the pelvis frame [rad/s].
    pub angular_velocity_body_rad_per_s: Vec3,
    /// Target minus current translational velocity in the pelvis frame [m/s].
    pub target_velocity_error_body_m_per_s: Vec3,
    /// Left and right foot contact-state indicators.
    pub foot_contact: [bool; 2],
    /// Gait phase in radians. Periodic wrapping is implicit in the basis.
    pub phase_rad: f64,
}

/// Refusal surface for the catalog-owned G1 policy representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G1PolicyError {
    /// An observation field contains NaN or infinity.
    NonFiniteObservation {
        /// Stable field name.
        field: &'static str,
    },
    /// The flat actuator-major matrix is not exactly 5,040 entries long.
    ParameterCount {
        /// Required parameter count.
        expected: usize,
        /// Supplied parameter count.
        actual: usize,
    },
    /// A policy weight contains NaN or infinity.
    NonFiniteParameter {
        /// Flat parameter index.
        index: usize,
    },
    /// A fixed-order activation produced a non-finite residual.
    NonFiniteOutput {
        /// Compact actuator index.
        actuator: usize,
    },
}

impl core::fmt::Display for G1PolicyError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFiniteObservation { field } => {
                write!(formatter, "non-finite G1 policy observation: {field}")
            }
            Self::ParameterCount { expected, actual } => write!(
                formatter,
                "G1 policy requires exactly {expected} parameters, got {actual}"
            ),
            Self::NonFiniteParameter { index } => {
                write!(formatter, "non-finite G1 policy parameter at index {index}")
            }
            Self::NonFiniteOutput { actuator } => {
                write!(
                    formatter,
                    "non-finite G1 policy output for actuator {actuator}"
                )
            }
        }
    }
}

impl std::error::Error for G1PolicyError {}

/// Expand one admitted observation into the exact 336-feature actuator basis.
///
/// Features are signal-major. Each of the 42 raw signals is multiplied by
/// `[1, sin(phi), cos(phi), sin(2phi), cos(2phi), sin(3phi), cos(3phi),
/// sin(4phi)]`. This makes the 5,040-D search space a typed part of the robot
/// owner rather than an undocumented browser convention.
pub fn g1_policy_features(
    observation: &G1PolicyObservation,
) -> Result<[f64; G1_POLICY_FEATURES_PER_ACTUATOR], G1PolicyError> {
    validate_g1_observation(observation)?;
    let basis = g1_policy_phase_basis(observation.phase_rad)?;
    let mut raw = [0.0; G1_POLICY_RAW_SIGNALS];
    raw[0] = 1.0;
    raw[1..16].copy_from_slice(&observation.joint_position_rad);
    raw[16..31].copy_from_slice(&observation.joint_velocity_rad_per_s);
    raw[31..34].copy_from_slice(&vec3_array(observation.gravity_direction_body));
    raw[34..37].copy_from_slice(&vec3_array(observation.angular_velocity_body_rad_per_s));
    raw[37..40].copy_from_slice(&vec3_array(observation.target_velocity_error_body_m_per_s));
    raw[40] = if observation.foot_contact[0] {
        1.0
    } else {
        0.0
    };
    raw[41] = if observation.foot_contact[1] {
        1.0
    } else {
        0.0
    };

    let mut features = [0.0; G1_POLICY_FEATURES_PER_ACTUATOR];
    for (signal_index, signal) in raw.iter().copied().enumerate() {
        let offset = signal_index * G1_POLICY_PHASE_BASIS;
        for (basis_index, multiplier) in basis.iter().copied().enumerate() {
            features[offset + basis_index] = signal * multiplier;
        }
    }
    Ok(features)
}

/// Evaluate the exact deterministic periodic basis shared by the G1 policy
/// and any owner-composed experiment that declares its phase schedule.
pub fn g1_policy_phase_basis(
    phase_rad: f64,
) -> Result<[f64; G1_POLICY_PHASE_BASIS], G1PolicyError> {
    if !phase_rad.is_finite() {
        return Err(G1PolicyError::NonFiniteObservation { field: "phase_rad" });
    }
    Ok([
        1.0,
        det::sin(phase_rad),
        det::cos(phase_rad),
        det::sin(2.0 * phase_rad),
        det::cos(2.0 * phase_rad),
        det::sin(3.0 * phase_rad),
        det::cos(3.0 * phase_rad),
        det::sin(4.0 * phase_rad),
    ])
}

/// An admitted, borrowed view of the exact 5,040-D G1 residual policy.
///
/// Construction validates the immutable parameter matrix once. A physical
/// rollout can then evaluate hundreds of observations without rescanning all
/// 5,040 weights at every integration step. The matrix is actuator-major
/// (`15 × 336`); each fixed-order dot product passes through deterministic
/// `tanh`, yielding a normalized residual in `[-1, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct G1ResidualPolicy<'a> {
    parameters: &'a [f64],
}

impl<'a> G1ResidualPolicy<'a> {
    /// Admit one finite parameter matrix in the catalog-owned layout.
    pub fn new(parameters: &'a [f64]) -> Result<Self, G1PolicyError> {
        if parameters.len() != G1_POLICY_DIMENSION {
            return Err(G1PolicyError::ParameterCount {
                expected: G1_POLICY_DIMENSION,
                actual: parameters.len(),
            });
        }
        for (index, value) in parameters.iter().copied().enumerate() {
            if !value.is_finite() {
                return Err(G1PolicyError::NonFiniteParameter { index });
            }
        }
        Ok(Self { parameters })
    }

    /// Evaluate one observation in the policy owner's exact feature order.
    pub fn evaluate(
        self,
        observation: &G1PolicyObservation,
    ) -> Result<[f64; G1_POLICY_ACTUATORS], G1PolicyError> {
        let features = g1_policy_features(observation)?;
        let mut output = [0.0; G1_POLICY_ACTUATORS];
        for (actuator, row) in self
            .parameters
            .chunks_exact(G1_POLICY_FEATURES_PER_ACTUATOR)
            .enumerate()
        {
            let mut activation = 0.0;
            for (weight, feature) in row.iter().copied().zip(features.iter().copied()) {
                activation += weight * feature;
            }
            output[actuator] = det::tanh(activation);
            if !output[actuator].is_finite() {
                return Err(G1PolicyError::NonFiniteOutput { actuator });
            }
        }
        Ok(output)
    }
}

fn validate_g1_observation(observation: &G1PolicyObservation) -> Result<(), G1PolicyError> {
    for (values, field) in [
        (
            observation.joint_position_rad.as_slice(),
            "joint_position_rad",
        ),
        (
            observation.joint_velocity_rad_per_s.as_slice(),
            "joint_velocity_rad_per_s",
        ),
    ] {
        if values.iter().any(|value| !value.is_finite()) {
            return Err(G1PolicyError::NonFiniteObservation { field });
        }
    }
    for (value, field) in [
        (observation.gravity_direction_body, "gravity_direction_body"),
        (
            observation.angular_velocity_body_rad_per_s,
            "angular_velocity_body_rad_per_s",
        ),
        (
            observation.target_velocity_error_body_m_per_s,
            "target_velocity_error_body_m_per_s",
        ),
    ] {
        if vec3_array(value)
            .iter()
            .any(|component| !component.is_finite())
        {
            return Err(G1PolicyError::NonFiniteObservation { field });
        }
    }
    if !observation.phase_rad.is_finite() {
        return Err(G1PolicyError::NonFiniteObservation { field: "phase_rad" });
    }
    Ok(())
}

const fn vec3_array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

/// Stable identity for one catalog model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobotModelId {
    /// Unitree G1 legacy 29-DoF source reduced to 12 leg and three waist joints.
    UnitreeG1LowerBodyWaist15,
    /// KUKA LBR iiwa 7 R800 description from the `iiwa_stack` iiwa7 Xacro.
    KukaLbrIiwa7R800,
}

/// One immutable upstream source record used by a catalog model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RobotModelSource {
    /// Purpose of this source within the transcription.
    pub role: &'static str,
    /// Content-stable URL pinned to the declared repository revision.
    pub url: &'static str,
    /// Full upstream Git revision containing the source.
    pub revision: &'static str,
    /// Git blob SHA-1 reported by the upstream repository for this file.
    pub git_blob_sha1: &'static str,
}

/// Honest provenance and reduction boundary for one built model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RobotModelMetadata {
    /// Typed catalog identity.
    pub id: RobotModelId,
    /// Stable human-readable catalog label.
    pub label: &'static str,
    /// Upstream organization or reference-model owner.
    pub authority: &'static str,
    /// Exact upstream source variant.
    pub source_variant: &'static str,
    /// Status of that source variant at the pinned revision.
    pub source_status: &'static str,
    /// Units used by all retained numeric records.
    pub units: &'static str,
    /// Exact derivation performed by this catalog builder.
    pub derivation: &'static str,
    /// Pinned files establishing topology, parameters, or source status.
    pub sources: &'static [RobotModelSource],
    /// Material omissions and approximation boundaries.
    pub omissions: &'static [&'static str],
}

/// Source joint label and hard-limit record in compact generalized order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobotJointMetadata {
    /// Stable source joint name.
    pub name: &'static str,
    /// Inclusive lower position limit in radians.
    pub lower_position_rad: f64,
    /// Inclusive upper position limit in radians.
    pub upper_position_rad: f64,
    /// Symmetric absolute velocity limit in radians per second.
    pub velocity_rad_per_second: f64,
    /// Symmetric absolute actuator-effort limit in newton-metres.
    pub effort_newton_metres: f64,
}

impl RobotJointMetadata {
    const fn new(
        name: &'static str,
        lower_position_rad: f64,
        upper_position_rad: f64,
        velocity_rad_per_second: f64,
        effort_newton_metres: f64,
    ) -> Self {
        Self {
            name,
            lower_position_rad,
            upper_position_rad,
            velocity_rad_per_second,
            effort_newton_metres,
        }
    }

    fn limits(self) -> Result<JointLimits, ArticulatedError> {
        JointLimits::new(
            self.lower_position_rad,
            self.upper_position_rad,
            self.velocity_rad_per_second,
            self.effort_newton_metres,
        )
    }
}

/// A validated articulated model with its stable source joint order and provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogRobotModel {
    model: ArticulatedModel,
    joints: Vec<RobotJointMetadata>,
    metadata: &'static RobotModelMetadata,
}

impl CatalogRobotModel {
    /// Validated articulated dynamics model.
    #[must_use]
    pub const fn model(&self) -> &ArticulatedModel {
        &self.model
    }

    /// Source joints in the exact compact generalized-coordinate order.
    #[must_use]
    pub fn joints(&self) -> &[RobotJointMetadata] {
        &self.joints
    }

    /// Pinned source and approximation metadata.
    #[must_use]
    pub const fn metadata(&self) -> &'static RobotModelMetadata {
        self.metadata
    }

    /// Consume the catalog wrapper and retain only the articulated model.
    #[must_use]
    pub fn into_model(self) -> ArticulatedModel {
        self.model
    }
}

const UNITREE_SOURCES: [RobotModelSource; 2] = [
    RobotModelSource {
        role: "selected link inertias, joint origins, axes, and hard limits",
        url: "https://raw.githubusercontent.com/unitreerobotics/unitree_ros/80b96438bcbc673379b5f35e767c40611f7e0af1/robots/g1_description/g1_29dof.urdf",
        revision: "80b96438bcbc673379b5f35e767c40611f7e0af1",
        git_blob_sha1: "49aeacd82072a6cb9847fac05ba59d9c35b0664e",
    },
    RobotModelSource {
        role: "variant inventory and upstream deprecation status",
        url: "https://raw.githubusercontent.com/unitreerobotics/unitree_ros/80b96438bcbc673379b5f35e767c40611f7e0af1/robots/g1_description/README.md",
        revision: "80b96438bcbc673379b5f35e767c40611f7e0af1",
        git_blob_sha1: "d1b8ce0f1859bf5f8ae8c4405533434f1d6c13be",
    },
];

const UNITREE_OMISSIONS: [&str; 4] = [
    "The 14 arm DoFs, arm links, fixed rubber hands, and their inertias are omitted rather than lumped into torso_link.",
    "Fixed pelvis-contour, logo, head, waist-support, and sensor attachments are omitted; their nonzero masses are not lumped into retained links.",
    "Visual/collision meshes, foot contact spheres, sensors, transmissions, and actuator/gear dynamics are not represented.",
    "The catalog itself prescribes no base boundary: callers may use fs-mbd's prescribed-base or unconstrained free-floating solve; neither adds the omitted upper-body mass nor ground contact.",
];

static UNITREE_METADATA: RobotModelMetadata = RobotModelMetadata {
    id: RobotModelId::UnitreeG1LowerBodyWaist15,
    label: "Unitree G1-inspired 15-DoF lower body and waist",
    authority: "Unitree Robotics official unitree_ros repository",
    source_variant: "g1_29dof.urdf, legacy mode_machine=2",
    source_status: "Deprecated by the pinned upstream README; retained as the explicitly requested source variant",
    units: "URDF SI: metres, kilograms, kilogram-metres squared, radians, radians/second, newton-metres",
    derivation: "Exact transcription of the pelvis, 12 leg, and three waist joint/link numeric records into fs-ga frames and fs-mbd articulated types",
    sources: &UNITREE_SOURCES,
    omissions: &UNITREE_OMISSIONS,
};

const IIWA_SOURCES: [RobotModelSource; 2] = [
    RobotModelSource {
        role: "link inertias, joint origins, axes, hard limits, and macro defaults",
        url: "https://raw.githubusercontent.com/IFL-CAMP/iiwa_stack/44f9d13c1b444d5dc9fd3e43ba60b7d3b2ea2bbb/iiwa_description/urdf/iiwa7.xacro",
        revision: "44f9d13c1b444d5dc9fd3e43ba60b7d3b2ea2bbb",
        git_blob_sha1: "44ebe189d96206503312183b84663bc255b6f39e",
    },
    RobotModelSource {
        role: "iiwa7 wrapper and default robot_name=iiwa",
        url: "https://raw.githubusercontent.com/IFL-CAMP/iiwa_stack/44f9d13c1b444d5dc9fd3e43ba60b7d3b2ea2bbb/iiwa_description/urdf/iiwa7.urdf.xacro",
        revision: "44f9d13c1b444d5dc9fd3e43ba60b7d3b2ea2bbb",
        git_blob_sha1: "6f2c4f87ffc98d2d12fb819cb2eed4d4d82d09a7",
    },
];

const IIWA_OMISSIONS: [&str; 4] = [
    "The wrapper's world link and fixed massless link_ee flange are omitted; iiwa_link_0 is attached to the prescribed BaseState.",
    "Visual/collision meshes, self-collision capsules, materials, Gazebo extensions, and transmissions are not represented.",
    "Joint damping and soft safety-controller limits are omitted; retained hard effort=300 N m and velocity=10 rad/s values are iiwa_stack macro defaults, not a KUKA-certified actuator envelope.",
    "No payload, tool, cable, motor/gear dynamics, contact, calibration, or manufacturer validation is inferred from the reference Xacro.",
];

static IIWA_METADATA: RobotModelMetadata = RobotModelMetadata {
    id: RobotModelId::KukaLbrIiwa7R800,
    label: "KUKA LBR iiwa 7 R800 reference arm",
    authority: "IFL-CAMP iiwa_stack reference description",
    source_variant: "iiwa7.xacro instantiated by iiwa7.urdf.xacro with default robot_name=iiwa",
    source_status: "Pinned community reference model; not a KUKA certification artifact",
    units: "Xacro/URDF SI: metres, kilograms, kilogram-metres squared, radians, radians/second, newton-metres",
    derivation: "Exact transcription of iiwa_link_0 through iiwa_link_7 and seven hard-limited revolute joints into fs-ga frames and fs-mbd articulated types",
    sources: &IIWA_SOURCES,
    omissions: &IIWA_OMISSIONS,
};

#[derive(Debug, Clone, Copy)]
enum SourceJoint {
    Fixed,
    Revolute {
        axis: [f64; 3],
        metadata: RobotJointMetadata,
    },
}

#[derive(Debug, Clone, Copy)]
struct SourceLink {
    name: &'static str,
    parent: Option<usize>,
    origin_xyz: [f64; 3],
    origin_rpy: [f64; 3],
    joint: SourceJoint,
    mass: f64,
    center_of_mass: [f64; 3],
    inertia: [f64; 6],
}

impl SourceLink {
    const fn fixed(
        name: &'static str,
        parent: Option<usize>,
        origin_xyz: [f64; 3],
        origin_rpy: [f64; 3],
        mass: f64,
        center_of_mass: [f64; 3],
        inertia: [f64; 6],
    ) -> Self {
        Self {
            name,
            parent,
            origin_xyz,
            origin_rpy,
            joint: SourceJoint::Fixed,
            mass,
            center_of_mass,
            inertia,
        }
    }

    // Keeping each source row self-contained makes transcription review against
    // the URDF materially less error-prone than splitting one record across
    // parallel tables.
    #[allow(clippy::too_many_arguments)]
    const fn revolute(
        name: &'static str,
        parent: usize,
        origin_xyz: [f64; 3],
        origin_rpy: [f64; 3],
        axis: [f64; 3],
        metadata: RobotJointMetadata,
        mass: f64,
        center_of_mass: [f64; 3],
        inertia: [f64; 6],
    ) -> Self {
        Self {
            name,
            parent: Some(parent),
            origin_xyz,
            origin_rpy,
            joint: SourceJoint::Revolute { axis, metadata },
            mass,
            center_of_mass,
            inertia,
        }
    }
}

// These decimals intentionally preserve the upstream URDF records instead of
// silently replacing them with nearby mathematical constants.
#[allow(clippy::approx_constant, clippy::unreadable_literal)]
const G1_LINKS: [SourceLink; 16] = [
    SourceLink::fixed(
        "pelvis",
        None,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        3.813,
        [0.0, 0.0, -0.07605],
        [0.010549, 0.0, 2.1e-6, 0.0093089, 0.0, 0.0079184],
    ),
    SourceLink::revolute(
        "left_hip_pitch_link",
        0,
        [0.0, 0.064452, -0.1027],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("left_hip_pitch_joint", -2.5307, 2.8798, 32.0, 88.0),
        1.35,
        [0.002741, 0.047791, -0.02606],
        [0.001811, 3.68e-5, -3.44e-5, 0.0014193, 0.000171, 0.0012812],
    ),
    SourceLink::revolute(
        "left_hip_roll_link",
        1,
        [0.0, 0.052, -0.030465],
        [0.0, -0.1749, 0.0],
        [1.0, 0.0, 0.0],
        RobotJointMetadata::new("left_hip_roll_joint", -0.5236, 2.9671, 32.0, 88.0),
        1.52,
        [0.029812, -0.001045, -0.087934],
        [
            0.0023773, -3.8e-6, -0.0003908, 0.0024123, 1.84e-5, 0.0016595,
        ],
    ),
    SourceLink::revolute(
        "left_hip_yaw_link",
        2,
        [0.025001, 0.0, -0.12412],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("left_hip_yaw_joint", -2.7576, 2.7576, 32.0, 88.0),
        1.702,
        [-0.057709, -0.010981, -0.15078],
        [
            0.0057774, -0.0005411, -0.0023948, 0.0076124, -0.0007072, 0.003149,
        ],
    ),
    SourceLink::revolute(
        "left_knee_link",
        3,
        [-0.078273, 0.0021489, -0.17734],
        [0.0, 0.1749, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("left_knee_joint", -0.087267, 2.8798, 20.0, 139.0),
        1.932,
        [0.005457, 0.003964, -0.12074],
        [0.011329, 4.82e-5, -4.49e-5, 0.011277, -0.0007146, 0.0015168],
    ),
    SourceLink::revolute(
        "left_ankle_pitch_link",
        4,
        [0.0, -9.4445e-5, -0.30001],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("left_ankle_pitch_joint", -0.87267, 0.5236, 30.0, 35.0),
        0.074,
        [-0.007269, 0.0, 0.011137],
        [8.4e-6, 0.0, -2.9e-6, 1.89e-5, 0.0, 1.26e-5],
    ),
    SourceLink::revolute(
        "left_ankle_roll_link",
        5,
        [0.0, 0.0, -0.017558],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        RobotJointMetadata::new("left_ankle_roll_joint", -0.2618, 0.2618, 30.0, 35.0),
        0.608,
        [0.026505, 0.0, -0.016425],
        [0.0002231, 2.0e-7, 8.91e-5, 0.0016161, -1.0e-7, 0.0016667],
    ),
    SourceLink::revolute(
        "right_hip_pitch_link",
        0,
        [0.0, -0.064452, -0.1027],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("right_hip_pitch_joint", -2.5307, 2.8798, 32.0, 88.0),
        1.35,
        [0.002741, -0.047791, -0.02606],
        [
            0.001811, -3.68e-5, -3.44e-5, 0.0014193, -0.000171, 0.0012812,
        ],
    ),
    SourceLink::revolute(
        "right_hip_roll_link",
        7,
        [0.0, -0.052, -0.030465],
        [0.0, -0.1749, 0.0],
        [1.0, 0.0, 0.0],
        RobotJointMetadata::new("right_hip_roll_joint", -2.9671, 0.5236, 32.0, 88.0),
        1.52,
        [0.029812, 0.001045, -0.087934],
        [
            0.0023773, 3.8e-6, -0.0003908, 0.0024123, -1.84e-5, 0.0016595,
        ],
    ),
    SourceLink::revolute(
        "right_hip_yaw_link",
        8,
        [0.025001, 0.0, -0.12412],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("right_hip_yaw_joint", -2.7576, 2.7576, 32.0, 88.0),
        1.702,
        [-0.057709, 0.010981, -0.15078],
        [
            0.0057774, 0.0005411, -0.0023948, 0.0076124, 0.0007072, 0.003149,
        ],
    ),
    SourceLink::revolute(
        "right_knee_link",
        9,
        [-0.078273, -0.0021489, -0.17734],
        [0.0, 0.1749, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("right_knee_joint", -0.087267, 2.8798, 20.0, 139.0),
        1.932,
        [0.005457, -0.003964, -0.12074],
        [0.011329, -4.82e-5, 4.49e-5, 0.011277, 0.0007146, 0.0015168],
    ),
    SourceLink::revolute(
        "right_ankle_pitch_link",
        10,
        [0.0, 9.4445e-5, -0.30001],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("right_ankle_pitch_joint", -0.87267, 0.5236, 30.0, 35.0),
        0.074,
        [-0.007269, 0.0, 0.011137],
        [8.4e-6, 0.0, -2.9e-6, 1.89e-5, 0.0, 1.26e-5],
    ),
    SourceLink::revolute(
        "right_ankle_roll_link",
        11,
        [0.0, 0.0, -0.017558],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        RobotJointMetadata::new("right_ankle_roll_joint", -0.2618, 0.2618, 30.0, 35.0),
        0.608,
        [0.026505, 0.0, -0.016425],
        [0.0002231, -2.0e-7, 8.91e-5, 0.0016161, 1.0e-7, 0.0016667],
    ),
    SourceLink::revolute(
        "waist_yaw_link",
        0,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("waist_yaw_joint", -2.618, 2.618, 32.0, 88.0),
        0.244,
        [0.003964, 0.0, 0.018769],
        [
            9.9587e-5, -1.833e-6, -1.2617e-5, 0.00012411, -1.18e-7, 0.00015586,
        ],
    ),
    SourceLink::revolute(
        "waist_roll_link",
        13,
        [-0.0039635, 0.0, 0.035],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        RobotJointMetadata::new("waist_roll_joint", -0.52, 0.52, 30.0, 35.0),
        0.047,
        [0.0, -0.000236, 0.010111],
        [7.515e-6, 0.0, 0.0, 6.398e-6, 9.9e-8, 3.988e-6],
    ),
    SourceLink::revolute(
        "torso_link",
        14,
        [0.0, 0.0, 0.019],
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        RobotJointMetadata::new("waist_pitch_joint", -0.52, 0.52, 30.0, 35.0),
        8.562,
        [0.002601, 0.000257, 0.153719],
        [
            0.065674966,
            -8.597e-5,
            -0.001737252,
            0.053535188,
            8.6899e-5,
            0.030808125,
        ],
    ),
];

const DEGREE: f64 = core::f64::consts::PI / 180.0;

// Literal grouping is intentionally identical to the Xacro decimal records.
#[allow(clippy::unreadable_literal)]
const IIWA_LINKS: [SourceLink; 8] = [
    SourceLink::fixed(
        "iiwa_link_0",
        None,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        5.0,
        [-0.1, 0.0, 0.07],
        [0.05, 0.0, 0.0, 0.06, 0.0, 0.03],
    ),
    SourceLink::revolute(
        "iiwa_link_1",
        0,
        [0.0, 0.0, 0.15],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_1", -170.0 * DEGREE, 170.0 * DEGREE, 10.0, 300.0),
        3.4525,
        [0.0, -0.03, 0.12],
        [0.02183, 0.0, 0.0, 0.007703, -0.003887, 0.02083],
    ),
    SourceLink::revolute(
        "iiwa_link_2",
        1,
        [0.0, 0.0, 0.19],
        [core::f64::consts::FRAC_PI_2, 0.0, core::f64::consts::PI],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_2", -120.0 * DEGREE, 120.0 * DEGREE, 10.0, 300.0),
        3.4821,
        [0.0003, 0.059, 0.042],
        [0.02076, 0.0, -0.003626, 0.02179, 0.0, 0.00779],
    ),
    SourceLink::revolute(
        "iiwa_link_3",
        2,
        [0.0, 0.21, 0.0],
        [core::f64::consts::FRAC_PI_2, 0.0, core::f64::consts::PI],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_3", -170.0 * DEGREE, 170.0 * DEGREE, 10.0, 300.0),
        4.05623,
        [0.0, 0.03, 0.13],
        [0.03204, 0.0, 0.0, 0.00972, 0.006227, 0.03042],
    ),
    SourceLink::revolute(
        "iiwa_link_4",
        3,
        [0.0, 0.0, 0.19],
        [core::f64::consts::FRAC_PI_2, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_4", -120.0 * DEGREE, 120.0 * DEGREE, 10.0, 300.0),
        3.4822,
        [0.0, 0.067, 0.034],
        [0.02178, 0.0, 0.0, 0.02075, -0.003625, 0.007785],
    ),
    SourceLink::revolute(
        "iiwa_link_5",
        4,
        [0.0, 0.21, 0.0],
        [-core::f64::consts::FRAC_PI_2, core::f64::consts::PI, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_5", -170.0 * DEGREE, 170.0 * DEGREE, 10.0, 300.0),
        2.1633,
        [0.0001, 0.021, 0.076],
        [0.01287, 0.0, 0.0, 0.005708, -0.003946, 0.01112],
    ),
    SourceLink::revolute(
        "iiwa_link_6",
        5,
        [0.0, 0.06070, 0.19],
        [core::f64::consts::FRAC_PI_2, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_6", -120.0 * DEGREE, 120.0 * DEGREE, 10.0, 300.0),
        2.3466,
        [0.0, 0.0006, 0.0004],
        [0.006509, 0.0, 0.0, 0.006259, 0.00031891, 0.004527],
    ),
    SourceLink::revolute(
        "iiwa_link_7",
        6,
        [0.0, 0.081, 0.06070],
        [-core::f64::consts::FRAC_PI_2, core::f64::consts::PI, 0.0],
        [0.0, 0.0, 1.0],
        RobotJointMetadata::new("iiwa_joint_7", -175.0 * DEGREE, 175.0 * DEGREE, 10.0, 300.0),
        3.129,
        [0.0, 0.0, 0.02],
        [0.01464, 0.0005912, 0.0, 0.01465, 0.0, 0.002872],
    ),
];

/// Build the requested Unitree G1-inspired 15-DoF lower-body and waist model.
///
/// # Errors
/// Returns the canonical articulated/Lie refusal if a transcribed source record
/// violates model, inertia, limit, or pose invariants.
pub fn unitree_g1_lower_body_waist_15dof() -> Result<CatalogRobotModel, ArticulatedError> {
    build_model(&G1_LINKS, &UNITREE_METADATA)
}

/// Build the seven-axis KUKA LBR iiwa 7 R800 reference model.
///
/// # Errors
/// Returns the canonical articulated/Lie refusal if a transcribed source record
/// violates model, inertia, limit, or pose invariants.
pub fn kuka_lbr_iiwa7_r800() -> Result<CatalogRobotModel, ArticulatedError> {
    build_model(&IIWA_LINKS, &IIWA_METADATA)
}

fn build_model(
    source_links: &[SourceLink],
    metadata: &'static RobotModelMetadata,
) -> Result<CatalogRobotModel, ArticulatedError> {
    let mut links = Vec::with_capacity(source_links.len());
    let mut joints = Vec::with_capacity(source_links.len().saturating_sub(1));
    for source in source_links {
        let joint = match source.joint {
            SourceJoint::Fixed => JointModel::FIXED,
            SourceJoint::Revolute { axis, metadata } => {
                joints.push(metadata);
                JointModel::revolute(vec3(axis), Some(metadata.limits()?))?
            }
        };
        links.push(Link::new(
            source.name,
            source.parent,
            pose(source.origin_xyz, source.origin_rpy)?,
            joint,
            inertia(source.mass, source.center_of_mass, source.inertia)?,
        ));
    }
    Ok(CatalogRobotModel {
        model: ArticulatedModel::new(links)?,
        joints,
        metadata,
    })
}

fn vec3(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn inertia(
    mass: f64,
    center_of_mass: [f64; 3],
    values: [f64; 6],
) -> Result<SpatialInertia, ArticulatedError> {
    let [ixx, ixy, ixz, iyy, iyz, izz] = values;
    SpatialInertia::new(
        mass,
        vec3(center_of_mass),
        Mat3 {
            m: [ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz],
        },
    )
}

fn pose(xyz: [f64; 3], rpy: [f64; 3]) -> Result<Se3, ArticulatedError> {
    let roll = So3::exp(So3Tangent::new(Vec3::new(rpy[0], 0.0, 0.0)))?;
    let pitch = So3::exp(So3Tangent::new(Vec3::new(0.0, rpy[1], 0.0)))?;
    let yaw = So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, rpy[2])))?;
    let rotation = yaw.compose(pitch)?.compose(roll)?;
    Ok(Se3::from_parts(rotation, vec3(xyz))?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::articulated::{BaseState, forward_dynamics, forward_kinematics};
    use fs_ga::Wrench;

    const TOLERANCE: f64 = 2.0e-10;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= TOLERANCE,
            "expected {expected:.16e}, got {actual:.16e}"
        );
    }

    fn neutral_kinematics(model: &CatalogRobotModel) -> crate::articulated::Kinematics {
        forward_kinematics(
            model.model(),
            BaseState::stationary(Se3::identity()),
            &vec![0.0; model.model().dof_count()],
            &vec![0.0; model.model().dof_count()],
        )
        .unwrap()
    }

    #[test]
    fn g1_catalog_has_stable_counts_and_source_order() {
        let catalog = unitree_g1_lower_body_waist_15dof().unwrap();
        assert_eq!(catalog.model().link_count(), 16);
        assert_eq!(catalog.model().dof_count(), 15);
        assert_eq!(
            catalog
                .model()
                .links()
                .iter()
                .map(Link::name)
                .collect::<Vec<_>>(),
            [
                "pelvis",
                "left_hip_pitch_link",
                "left_hip_roll_link",
                "left_hip_yaw_link",
                "left_knee_link",
                "left_ankle_pitch_link",
                "left_ankle_roll_link",
                "right_hip_pitch_link",
                "right_hip_roll_link",
                "right_hip_yaw_link",
                "right_knee_link",
                "right_ankle_pitch_link",
                "right_ankle_roll_link",
                "waist_yaw_link",
                "waist_roll_link",
                "torso_link",
            ]
        );
        assert_eq!(
            catalog
                .joints()
                .iter()
                .map(|joint| joint.name)
                .collect::<Vec<_>>(),
            [
                "left_hip_pitch_joint",
                "left_hip_roll_joint",
                "left_hip_yaw_joint",
                "left_knee_joint",
                "left_ankle_pitch_joint",
                "left_ankle_roll_joint",
                "right_hip_pitch_joint",
                "right_hip_roll_joint",
                "right_hip_yaw_joint",
                "right_knee_joint",
                "right_ankle_pitch_joint",
                "right_ankle_roll_joint",
                "waist_yaw_joint",
                "waist_roll_joint",
                "waist_pitch_joint",
            ]
        );
    }

    #[test]
    fn g1_policy_is_exactly_5040_dimensional_and_owner_ordered() -> Result<(), G1PolicyError> {
        assert_eq!(G1_POLICY_RAW_SIGNALS, 42);
        assert_eq!(G1_POLICY_FEATURES_PER_ACTUATOR, 336);
        assert_eq!(G1_POLICY_DIMENSION, 5_040);
        let mut observation = G1PolicyObservation {
            joint_position_rad: [0.0; G1_POLICY_ACTUATORS],
            joint_velocity_rad_per_s: [0.0; G1_POLICY_ACTUATORS],
            gravity_direction_body: Vec3::new(0.0, 0.0, -1.0),
            angular_velocity_body_rad_per_s: Vec3::new(0.0, 0.0, 0.0),
            target_velocity_error_body_m_per_s: Vec3::new(1.0, 0.0, 0.0),
            foot_contact: [true, false],
            phase_rad: 0.0,
        };
        observation.joint_position_rad[0] = 2.0;
        let features = g1_policy_features(&observation)?;
        assert_eq!(&features[0..8], &[1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]);
        assert_eq!(&features[8..16], &[2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0]);
        assert_eq!(features[40 * G1_POLICY_PHASE_BASIS], 1.0);
        assert_eq!(features[41 * G1_POLICY_PHASE_BASIS], 0.0);

        let zero_parameters = vec![0.0; G1_POLICY_DIMENSION];
        let policy = G1ResidualPolicy::new(&zero_parameters)?;
        assert_eq!(policy.evaluate(&observation)?, [0.0; G1_POLICY_ACTUATORS]);
        Ok(())
    }

    #[test]
    fn every_g1_policy_coordinate_reaches_exactly_one_actuator() -> Result<(), G1PolicyError> {
        let observation = G1PolicyObservation {
            joint_position_rad: core::array::from_fn(|index| 0.05 * (index + 1) as f64),
            joint_velocity_rad_per_s: core::array::from_fn(|index| -0.03 * (index + 1) as f64),
            gravity_direction_body: Vec3::new(0.20, -0.30, -0.93),
            angular_velocity_body_rad_per_s: Vec3::new(0.11, -0.17, 0.23),
            target_velocity_error_body_m_per_s: Vec3::new(0.41, -0.29, 0.13),
            foot_contact: [true, true],
            phase_rad: 0.37,
        };
        let features = g1_policy_features(&observation)?;
        assert!(features.iter().all(|feature| feature.abs() > 1.0e-9));
        let mut parameters = vec![0.0; G1_POLICY_DIMENSION];
        for coordinate in 0..G1_POLICY_DIMENSION {
            parameters[coordinate] = 0.25;
            let output = G1ResidualPolicy::new(&parameters)?.evaluate(&observation)?;
            let expected_actuator = coordinate / G1_POLICY_FEATURES_PER_ACTUATOR;
            let expected_feature = coordinate % G1_POLICY_FEATURES_PER_ACTUATOR;
            for (actuator, value) in output.into_iter().enumerate() {
                let expected = if actuator == expected_actuator {
                    det::tanh(0.25 * features[expected_feature])
                } else {
                    0.0
                };
                assert_eq!(
                    value.to_bits(),
                    expected.to_bits(),
                    "coordinate {coordinate}"
                );
            }
            parameters[coordinate] = 0.0;
        }
        Ok(())
    }

    #[test]
    fn g1_policy_refuses_shape_and_non_finite_inputs() {
        let observation = G1PolicyObservation {
            joint_position_rad: [0.0; G1_POLICY_ACTUATORS],
            joint_velocity_rad_per_s: [0.0; G1_POLICY_ACTUATORS],
            gravity_direction_body: Vec3::new(0.0, 0.0, -1.0),
            angular_velocity_body_rad_per_s: Vec3::new(0.0, 0.0, 0.0),
            target_velocity_error_body_m_per_s: Vec3::new(0.0, 0.0, 0.0),
            foot_contact: [false; 2],
            phase_rad: 0.0,
        };
        assert!(matches!(
            G1ResidualPolicy::new(&[]),
            Err(G1PolicyError::ParameterCount {
                expected: G1_POLICY_DIMENSION,
                actual: 0
            })
        ));
        let mut parameters = vec![0.0; G1_POLICY_DIMENSION];
        parameters[1_234] = f64::NAN;
        assert_eq!(
            G1ResidualPolicy::new(&parameters),
            Err(G1PolicyError::NonFiniteParameter { index: 1_234 })
        );
        assert_eq!(
            g1_policy_features(&observation).unwrap().len(),
            G1_POLICY_FEATURES_PER_ACTUATOR
        );
    }

    #[test]
    fn g1_neutral_forward_kinematics_preserves_bilateral_origin_symmetry() {
        let catalog = unitree_g1_lower_body_waist_15dof().unwrap();
        let kinematics = neutral_kinematics(&catalog);
        for (left, right) in (1..=6).zip(7..=12) {
            let left_origin = kinematics.world_from_link[left]
                .transform_point(Vec3::new(0.0, 0.0, 0.0))
                .unwrap();
            let right_origin = kinematics.world_from_link[right]
                .transform_point(Vec3::new(0.0, 0.0, 0.0))
                .unwrap();
            assert_close(left_origin.x, right_origin.x);
            assert_close(left_origin.y, -right_origin.y);
            assert_close(left_origin.z, right_origin.z);
        }
    }

    #[test]
    fn iiwa_catalog_has_stable_counts_order_and_neutral_endpoint() {
        let catalog = kuka_lbr_iiwa7_r800().unwrap();
        assert_eq!(catalog.model().link_count(), 8);
        assert_eq!(catalog.model().dof_count(), 7);
        assert_eq!(
            catalog
                .model()
                .links()
                .iter()
                .map(Link::name)
                .collect::<Vec<_>>(),
            [
                "iiwa_link_0",
                "iiwa_link_1",
                "iiwa_link_2",
                "iiwa_link_3",
                "iiwa_link_4",
                "iiwa_link_5",
                "iiwa_link_6",
                "iiwa_link_7",
            ]
        );
        assert_eq!(
            catalog
                .joints()
                .iter()
                .map(|joint| joint.name)
                .collect::<Vec<_>>(),
            [
                "iiwa_joint_1",
                "iiwa_joint_2",
                "iiwa_joint_3",
                "iiwa_joint_4",
                "iiwa_joint_5",
                "iiwa_joint_6",
                "iiwa_joint_7",
            ]
        );
        let kinematics = neutral_kinematics(&catalog);
        let endpoint = kinematics.world_from_link[7]
            .transform_point(Vec3::new(0.0, 0.0, 0.0))
            .unwrap();
        assert_close(endpoint.x, 0.0);
        assert_close(endpoint.y, 0.0);
        assert_close(endpoint.z, 1.2210);
    }

    #[test]
    fn catalog_limits_inertias_and_linear_complexity_are_admitted() {
        for catalog in [
            unitree_g1_lower_body_waist_15dof().unwrap(),
            kuka_lbr_iiwa7_r800().unwrap(),
        ] {
            assert_eq!(catalog.joints().len(), catalog.model().dof_count());
            for joint in catalog.joints() {
                assert!(joint.lower_position_rad <= 0.0);
                assert!(joint.upper_position_rad >= 0.0);
                assert!(joint.velocity_rad_per_second.is_finite());
                assert!(joint.velocity_rad_per_second > 0.0);
                assert!(joint.effort_newton_metres.is_finite());
                assert!(joint.effort_newton_metres > 0.0);
            }
            for link in catalog.model().links() {
                assert!(link.inertia().mass().is_finite());
                assert!(link.inertia().mass() > 0.0);
                assert!(link.inertia().matrix().is_ok());
            }
            let lower = catalog
                .joints()
                .iter()
                .map(|joint| joint.lower_position_rad)
                .collect::<Vec<_>>();
            let upper = catalog
                .joints()
                .iter()
                .map(|joint| joint.upper_position_rad)
                .collect::<Vec<_>>();
            let zero = vec![0.0; catalog.model().dof_count()];
            assert!(
                forward_kinematics(
                    catalog.model(),
                    BaseState::stationary(Se3::identity()),
                    &lower,
                    &zero,
                )
                .is_ok()
            );
            assert!(
                forward_kinematics(
                    catalog.model(),
                    BaseState::stationary(Se3::identity()),
                    &upper,
                    &zero,
                )
                .is_ok()
            );
            let external = vec![Wrench::default(); catalog.model().link_count()];
            let dynamics = forward_dynamics(
                catalog.model(),
                BaseState::stationary(Se3::identity()),
                &zero,
                &zero,
                &zero,
                Vec3::new(0.0, 0.0, 0.0),
                &external,
            )
            .unwrap();
            assert_eq!(dynamics.generalized_acceleration, zero);
            let complexity = catalog.model().complexity();
            assert_eq!(complexity.dense_generalized_matrix_entries, 0);
            assert_eq!(complexity.degrees_of_freedom, catalog.joints().len());
            assert_eq!(complexity.links, catalog.model().link_count());
        }
    }

    #[test]
    fn catalog_rebuild_is_deterministic_and_preserves_provenance() {
        let g1_first = unitree_g1_lower_body_waist_15dof().unwrap();
        let g1_second = unitree_g1_lower_body_waist_15dof().unwrap();
        assert_eq!(g1_first, g1_second);
        assert_eq!(g1_first.metadata().sources.len(), 2);
        assert!(g1_first.metadata().source_status.contains("Deprecated"));

        let iiwa_first = kuka_lbr_iiwa7_r800().unwrap();
        let iiwa_second = kuka_lbr_iiwa7_r800().unwrap();
        assert_eq!(iiwa_first, iiwa_second);
        assert_eq!(iiwa_first.metadata().sources.len(), 2);
        assert!(
            iiwa_first
                .metadata()
                .omissions
                .iter()
                .any(|item| item.contains("300 N m"))
        );
    }
}
