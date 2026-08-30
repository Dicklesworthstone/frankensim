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
    CatalogRobotModel, G1_MODEL_ACTUATORS, G1_POLICY_ACTUATORS, G1_POLICY_DIMENSION,
    G1_POLICY_FEATURES_PER_ACTUATOR, G1PolicyObservation, G1ResidualPolicy, g1_policy_phase_basis,
    unitree_g1_29dof,
};

// ─── Learned-policy hook (PPO bridge; used by the fs-g1-train adapter) ───

/// Observation dimensionality of the flattened G1 policy observation:
/// 15 joint positions + 15 joint velocities + 3 gravity + 3 angular
/// velocity + 3 velocity error + 2 contacts + 1 phase.
pub const G1_LEARNED_OBS_DIMS: usize = 42;

/// PPO transition data for one learned rollout episode. The rewards are
/// the RL shaping signal (survival + forward progress − penalties); the
/// CMA-ES evaluation objective in [`G1WalkingReceipt`] remains the
/// separate evaluation metric.
#[derive(Debug, Clone)]
pub struct EpisodeTrace {
    pub rewards: Vec<f32>,
    /// True when the step ended in a terminal guard (fall).
    pub done: Vec<bool>,
    pub completed_steps: usize,
    pub termination: G1TerminationReason,
    pub objective: f64,
}

impl Default for EpisodeTrace {
    fn default() -> Self {
        Self {
            rewards: Vec::new(),
            done: Vec::new(),
            completed_steps: 0,
            termination: G1TerminationReason::Horizon,
            objective: 0.0,
        }
    }
}

/// Anything that can drive the G1 rollout step-by-step (per-step
/// observation → 29-dim residual action). Implemented by the
/// fs-g1-train transformer adapter behind the `g1-learned` feature.
pub trait LearnedG1Policy {
    /// Observe, act, return the executed 29-dim residual action.
    fn act(&mut self, obs: &G1PolicyObservation, step: usize)
        -> [f64; G1_POLICY_ACTUATORS];
}

/// Flatten the kernel-owned observation into the transformer's input
/// order (f32).
#[must_use]
pub fn flatten_observation(obs: &G1PolicyObservation) -> [f32; G1_LEARNED_OBS_DIMS] {
    let mut out = [0.0f32; G1_LEARNED_OBS_DIMS];
    let mut k = 0usize;
    for v in obs.joint_position_rad.iter() {
        out[k] = *v as f32;
        k += 1;
    }
    for v in obs.joint_velocity_rad_per_s.iter() {
        out[k] = *v as f32;
        k += 1;
    }
    for v in [
        obs.gravity_direction_body.x,
        obs.gravity_direction_body.y,
        obs.gravity_direction_body.z,
    ] {
        out[k] = v as f32;
        k += 1;
    }
    for v in [
        obs.angular_velocity_body_rad_per_s.x,
        obs.angular_velocity_body_rad_per_s.y,
        obs.angular_velocity_body_rad_per_s.z,
    ] {
        out[k] = v as f32;
        k += 1;
    }
    for v in [
        obs.target_velocity_error_body_m_per_s.x,
        obs.target_velocity_error_body_m_per_s.y,
        obs.target_velocity_error_body_m_per_s.z,
    ] {
        out[k] = v as f32;
        k += 1;
    }
    out[k] = f32::from(u8::from(obs.foot_contact[0]));
    out[k + 1] = f32::from(u8::from(obs.foot_contact[1]));
    out[k + 2] = obs.phase_rad as f32;
    out
}

use fs_time::{RenormPolicy, se3_exp_step_renorm};
use fs_tribo::{
    ContactFrame, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialSlip,
};

/// Stable identity of the owner-composed walking experiment.
pub const G1_WALKING_MODEL_ID: &str = "fs-cmaes/g1-walking-owner-composition-v9";
/// Pelvis plus all 29 actuated links retained by the source-bound mode-11 catalog.
pub const G1_LINK_COUNT: usize = 30;
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
const MINIMUM_UPRIGHT_HEIGHT_M: f64 = 0.55;
const MINIMUM_UPRIGHT_GRAVITY_Z: f64 = -0.866_025_403_784_438_6;
const GAIT_SWITCH_WINDOW: f64 = 0.20;
const TARGET_SWING_CLEARANCE_M: f64 = 0.055;
const POSTURE_HEIGHT_SCALE_M: f64 = 0.12;
const POSTURE_TILT_SINE_SCALE: f64 = 0.35;
const LATERAL_POSITION_SCALE_M: f64 = 0.15;
const HEADING_SINE_SCALE: f64 = 0.35;
// Survival is lexicographically primary for this walking experiment. The
// secondary physical shaping score is smoothly bounded below half of one
// horizon-step charge, so one additional survived step dominates every
// possible difference between shaping scores.
const UNCOMPLETED_STEP_PENALTY: f64 = 1_000.0;
const SHAPING_SCORE_LIMIT: f64 = 400.0;
const SHAPING_SCORE_SCALE: f64 = 200.0;
// Per-step survival bonus (cmaes-pvz, v068): a small reward per survived step
// inside the same objective, so the optimizer cannot game a single rollup by
// collapsing early. Bounded to ~half the full horizon by design
// (0.4 * 720 ≈ 288, below SHAPING_SCORE_LIMIT).
const SURVIVAL_BONUS_PER_STEP: f64 = 0.4;
const MAX_BONUSED_STEPS: f64 = 720.0;
const SURVIVAL_BONUS_LIMIT: f64 = SURVIVAL_BONUS_PER_STEP * MAX_BONUSED_STEPS;
/// Peak height of the disclosed smooth challenge terrain [m].
pub const G1_TERRAIN_AMPLITUDE_M: f64 = 0.008;
/// Longitudinal spatial frequency of the disclosed terrain [rad/m].
pub const G1_TERRAIN_WAVENUMBER_RAD_PER_M: f64 = 2.4;
/// Start of the deterministic lateral push pulse [s].
pub const G1_PUSH_START_S: f64 = 0.55;
/// End of the deterministic lateral push pulse [s].
pub const G1_PUSH_END_S: f64 = 0.70;
/// Peak of the smooth half-sine lateral push [N].
pub const G1_PUSH_PEAK_FORCE_N: f64 = 24.0;
const RECOVERY_TILT_SINE: f64 = 0.10;
const RECOVERY_ANGULAR_SPEED_RAD_PER_S: f64 = 0.35;
const RECOVERY_HEIGHT_ERROR_M: f64 = 0.055;
const ARM_SWING_AMPLITUDE_RAD: f64 = 0.30;
const ARM_ROLL_DAMPING_TARGET_S: f64 = 0.10;
// v069 (cmaes-zi6): the v067 whole-body catalog adds arms/head/hands as
// integrated dynamic bodies. The v066 16-link model had the same 0.30 rad
// arm-swing reflex but those joints were display-only, so the swing cost
// nothing physically. On the 30-link catalog the swing injects ~12 kg of
// upper-body inertia that the disclosed curriculum was not calibrated for,
// pushing the standing prior and the 105-coordinate curriculum off-balance
// during the first gait cycle. Gate the swing on a smoothstep that engages
// only after the stabilizer has had a full cycle to settle:
//   swing_scale(time_s) = smoothstep(cycle_period / 2, cycle_period * 1.5, time_s)
// i.e. arms stay quiet through the balance phase (0..1 cycle) and ramp in
// v069 (cmaes-zi6): 3-stage curriculum retune for the 30-link whole-body
// catalog. Generated by tests/retune_three_stage_curriculum_v069 on the
// v0.6.9 dynamics (which include the arm-swing gate). The retune runs a
// curriculum-length CMA on the disclosed owner from three progressive flat
// horizons: 0.65s -> 1.0s -> 1.5s, with 80 / 100 / 120 generations
// (sigma 0.08 / 0.055 / 0.04) respectively. The output is a curriculum
// that completes the 720-step flat horizon at 0.488m forward (vs the
// v0.6.6 baked-in mean that only completed ~368 steps under the v0.6.9
// 30-link whole-body dynamics). The 4-stage test (which also retunes
// for terrain-push 1.5s) is in recalibrate_mode_11_curriculum and is
// the source-of-truth for future refreshes; this comment records that
// the values below came from the 3-stage retune, seed 0x47315060 /
// 0x47315061 / 0x47315062, finished at generation 120 of the flat-1.5s
// stage.
const G1_STABILIZING_BIAS_MEAN: [f64; G1_POLICY_ACTUATORS] = [
    -0.111_299_477_068_413_12,
    -0.164_614_709_893_849_83,
    0.101_713_234_254_439_80,
    -0.004_143_539_605_537_37,
    -0.533_936_535_150_009_61,
    0.418_212_872_978_467_65,
    -0.027_985_376_260_349_51,
    0.224_877_699_522_546_71,
    -0.013_832_946_047_615_55,
    0.209_962_701_522_745_27,
    0.090_371_040_921_585_31,
    -0.084_632_404_448_831_48,
    -0.199_925_818_284_272_36,
    0.302_129_556_691_599_80,
    -0.284_169_255_264_625_55,
];

// Phase-only curriculum coordinates learned on the 0.5 s stepping task by
// full CMA (seed 0x4731_5040, lambda 16, 100 generations, sigma 0.025).
// Values are actuator-major [sin(phi), cos(phi)].
const G1_WALKING_PHASE_MEAN: [f64; 2 * G1_POLICY_ACTUATORS] = [
    0.297_740_446_498_193_42,
    -0.130_933_244_373_483_17,
    -0.301_344_913_345_272_86,
    -0.429_963_214_985_060_35,
    -0.090_609_294_789_346_43,
    -0.502_464_919_529_171_61,
    0.193_539_440_851_068_92,
    -0.121_732_408_717_522_36,
    0.409_108_744_809_110_69,
    -0.082_763_683_712_551_05,
    -0.015_120_061_106_034_79,
    0.537_346_217_799_751_00,
    -0.090_880_255_201_175_82,
    0.035_440_380_901_751_06,
    0.064_198_240_286_287_46,
    -0.168_895_678_796_116_75,
    0.151_788_892_945_097_68,
    -0.070_199_598_017_847_20,
    -0.591_356_277_275_032_70,
    -0.044_732_140_575_658_15,
    -0.028_243_322_313_909_00,
    0.184_665_776_509_056_97,
    -0.082_228_947_832_493_71,
    -0.059_939_899_355_855_90,
    0.226_912_599_550_020_99,
    -0.003_913_658_056_743_41,
    0.138_976_357_760_920_68,
    0.331_239_351_896_104_26,
    -0.198_028_809_970_659_77,
    -0.122_191_118_416_002_69,
];

// Balance-feedback curriculum coordinates. Full CMA first learned the 0.9 s
// stepping handoff (seed 0x4731_5042, lambda 16, 120 generations, sigma 0.14),
// then optimized the same coordinates on the 1.5 s walking task (seed
// 0x4731_5044, lambda 16, 180 generations, sigma 0.055). Values are
// actuator-major [gravity-x, gravity-y, angular-velocity-x,
// angular-velocity-y], each on the owner's constant phase basis.
const G1_WALKING_FEEDBACK_MEAN: [f64; 4 * G1_POLICY_ACTUATORS] = [
    0.076_624_194_149_865_99,
    -0.410_807_257_648_436_17,
    0.151_254_531_845_033_68,
    0.422_856_068_745_059_43,
    -1.343_021_997_940_628_60,
    -0.813_119_072_957_889_43,
    0.869_902_361_429_454_98,
    -0.268_640_219_659_509_07,
    0.100_445_284_892_340_64,
    -0.222_972_214_479_364_49,
    -0.163_634_137_236_706_97,
    -0.139_737_095_168_433_49,
    0.373_703_113_181_448_83,
    0.688_015_752_476_764_83,
    0.213_335_900_891_140_27,
    -0.048_939_208_660_400_20,
    0.280_960_920_192_869_95,
    -0.015_062_634_112_024_10,
    0.080_107_455_396_772_01,
    0.542_466_843_822_770_88,
    0.149_227_122_360_966_97,
    -0.101_043_261_285_806_01,
    1.083_168_751_262_241_44,
    -0.688_781_297_932_282_04,
    0.237_670_606_014_979_64,
    0.380_523_865_024_808_23,
    0.461_941_739_931_357_23,
    -0.012_244_038_958_358_43,
    -0.170_602_168_824_592_53,
    0.036_130_616_536_442_13,
    -0.162_849_215_993_117_00,
    -0.070_639_165_779_644_30,
    -0.628_336_541_737_316_35,
    0.245_205_894_175_058_81,
    -0.830_133_280_717_569_42,
    1.012_124_756_456_735_72,
    0.726_383_770_778_525_12,
    0.121_374_616_449_700_18,
    0.417_023_236_659_832_37,
    -0.073_405_438_946_443_66,
    -0.285_098_750_556_292_75,
    -0.303_354_210_739_048_84,
    -0.131_789_306_255_366_11,
    1.202_877_841_097_499_62,
    0.198_631_893_454_957_80,
    -0.410_607_593_935_011_39,
    0.376_808_032_075_699_85,
    0.229_450_739_832_338_19,
    -1.238_939_554_978_532_70,
    -0.052_283_560_431_461_73,
    0.545_764_754_936_013_40,
    0.212_501_857_869_883_93,
    -0.067_044_247_609_571_44,
    0.322_431_652_357_363_00,
    -0.040_718_751_643_042_51,
    0.109_660_908_964_203_45,
    -0.151_259_923_230_684_46,
    -0.237_230_715_937_692_63,
    -0.066_234_117_226_295_20,
    -0.634_418_770_376_905_75,
];

// The standing PD controller must leave enough owner-model motor authority for
// a black-box policy to bend a swing knee and unload a foot. The earlier 0.32
// fraction could not overcome the knee posture gain over a useful excursion.
// The disclosed stabilizing biases above are analytically rescaled so
// `0.65 * tanh(new_bias) == 0.32 * tanh(previous_bias)`: the initial residual
// torque command is unchanged while learned excursions gain an honest,
// symmetric authority envelope.
const RESIDUAL_EFFORT_FRACTION: f64 = 0.65;

/// Declared curriculum task. Balance and walking are different black-box
/// objectives and are never silently mixed under one receipt.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G1Task {
    /// Track zero velocity while upright in double support.
    Balance = 0,
    /// Track alternating support and swing clearance without translation.
    Stepping = 1,
    /// Track forward speed and the alternating-support gait schedule.
    Walking = 2,
}

/// Fixed environment challenge, separate from the mutable policy coordinates.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum G1Challenge {
    /// Horizontal ground with no external perturbation.
    Flat = 0,
    /// Smooth height field plus one disclosed deterministic lateral push.
    TerrainAndPush = 1,
}

/// Owner-layout 5,040-D mean learned by the disclosed deterministic balance
/// bootstrap. This is initialization data, not a hidden controller: 5,025
/// coordinates remain exactly zero and CMA-ES may change every coordinate.
#[must_use]
pub fn g1_stabilizing_policy_mean() -> [f64; G1_POLICY_DIMENSION] {
    let mut policy = [0.0; G1_POLICY_DIMENSION];
    for (actuator, bias) in G1_STABILIZING_BIAS_MEAN.iter().copied().enumerate() {
        policy[actuator * G1_POLICY_FEATURES_PER_ACTUATOR] = bias;
    }
    policy
}

/// Disclosed owner-layout walking curriculum mean used to initialize the live
/// 5,040-D scalable-CMA refinement. Exactly 105 coordinates are nonzero: 15
/// standing biases, 30 periodic gait weights, and 60 pelvis-state feedback
/// weights. This is initialization data, not a hidden trajectory or browser
/// controller; every owner policy coordinate remains mutable by CMA.
#[must_use]
pub fn g1_walking_curriculum_mean() -> [f64; G1_POLICY_DIMENSION] {
    let mut policy = g1_stabilizing_policy_mean();
    for actuator in 0..G1_POLICY_ACTUATORS {
        let row = actuator * G1_POLICY_FEATURES_PER_ACTUATOR;
        policy[row + 1] = G1_WALKING_PHASE_MEAN[2 * actuator];
        policy[row + 2] = G1_WALKING_PHASE_MEAN[2 * actuator + 1];
        for (feedback, signal) in [31, 32, 34, 35].into_iter().enumerate() {
            policy[row + signal * 8] = G1_WALKING_FEEDBACK_MEAN[4 * actuator + feedback];
        }
    }
    policy
}

/// Fixed, public experiment controls. They are intentionally not CMA search
/// coordinates: changing them defines a different black-box problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct G1WalkingConfig {
    /// Explicit balance or walking curriculum task.
    pub task: G1Task,
    /// Fixed terrain/perturbation contract for every candidate in the run.
    pub challenge: G1Challenge,
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
            task: G1Task::Walking,
            challenge: G1Challenge::Flat,
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
    /// Integrated squared excess load above body weight, normalized by body weight.
    pub impact_integral: f64,
    /// Integrated backward pelvis travel [m].
    pub backward_distance_m: f64,
    /// Integrated squared lateral position, normalized by 0.15 m [s].
    pub lateral_error_integral: f64,
    /// Integrated squared heading sine, normalized by 0.35 [s].
    pub heading_error_integral: f64,
    /// Integrated disagreement with the declared alternating-support schedule [s].
    pub contact_schedule_mismatch_integral: f64,
    /// Integrated squared swing-sole clearance error [m² s].
    pub swing_clearance_error_integral: f64,
    /// Time with exactly one foot contacting the ground [s].
    pub single_support_s: f64,
    /// Time with both feet contacting the ground [s].
    pub double_support_s: f64,
    /// Time with neither foot contacting the ground [s].
    pub flight_s: f64,
    /// Time integral of the applied push magnitude [N s].
    pub push_impulse_n_s: f64,
    /// Delay after the push until the declared recovery envelope is regained [s].
    pub recovery_time_s: f64,
    /// Minimum pelvis height over the rollout [m].
    pub minimum_base_height_m: f64,
    /// Maximum pelvis tilt sine over the rollout.
    pub maximum_tilt_sine: f64,
    /// Maximum absolute terrain height visited by a contact sample [m].
    pub maximum_abs_terrain_height_m: f64,
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
    reference_position: [f64; G1_MODEL_ACTUATORS],
    initial_base_height_m: f64,
    total_mass_kg: f64,
    step_count: usize,
}

impl G1WalkingEvaluator {
    /// Admit the fixed experiment and build the source-bound model once.
    pub fn new(config: G1WalkingConfig) -> Result<Self, G1WalkingError> {
        validate_config(config)?;
        let catalog = unitree_g1_29dof()?;
        let reference_position = reference_posture();
        let initial_base_height_m = initial_height(&catalog, &reference_position)?;
        let total_mass_kg = catalog
            .model()
            .links()
            .iter()
            .map(|link| link.inertia().mass())
            .sum();
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
        let step_count = rounded_step_count(config)?;
        Ok(Self {
            config,
            catalog,
            interface,
            friction,
            reference_position,
            initial_base_height_m,
            total_mass_kg,
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
        self.rollout_inner(parameters, retain_trace, None, None)
    }

    /// Rollout driven by a LEARNED policy (feature `g1-learned`): records
    /// the PPO transition data (rewards / done flags / final objective)
    /// into `trace` while the passed policy acts each step. The returned
    /// receipt carries the unchanged CMA-ES evaluation objective.
    pub fn rollout_learned(
        &self,
        policy: &mut dyn LearnedG1Policy,
        trace: &mut EpisodeTrace,
    ) -> Result<G1WalkingReceipt, G1WalkingError> {
        // G1ResidualPolicy::new validates finiteness only — the stabilizing
        // mean is a valid placeholder; the learned hook replaces it.
        let parameters = g1_stabilizing_policy_mean();
        self.rollout_inner(&parameters, false, Some(policy), Some(trace))
    }

    fn rollout_inner(
        &self,
        parameters: &[f64],
        retain_trace: bool,
        mut learned: Option<&mut dyn LearnedG1Policy>,
        mut learned_trace: Option<&mut EpisodeTrace>,
    ) -> Result<G1WalkingReceipt, G1WalkingError> {
        let policy = G1ResidualPolicy::new(parameters)?;
        let model = self.catalog.model();
        let mut position = self.reference_position;
        let mut velocity = [0.0; G1_MODEL_ACTUATORS];
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
        let mut backward_distance_m = 0.0;
        let mut lateral_error_integral = 0.0;
        let mut heading_error_integral = 0.0;
        let mut contact_schedule_mismatch_integral = 0.0;
        let mut swing_clearance_error_integral = 0.0;
        let mut single_support_s = 0.0;
        let mut double_support_s = 0.0;
        let mut flight_s = 0.0;
        let mut push_impulse_n_s = 0.0;
        let mut recovery_time_s = if self.config.challenge == G1Challenge::Flat {
            0.0
        } else {
            (self.config.duration_s - G1_PUSH_END_S).max(0.0)
        };
        let mut recovery_recorded = self.config.challenge == G1Challenge::Flat;
        let mut minimum_base_height_m = self.initial_base_height_m;
        let mut maximum_tilt_sine = 0.0_f64;
        let mut maximum_abs_terrain_height_m = 0.0_f64;
        let mut completed_steps = 0;
        let mut termination_reason = G1TerminationReason::Horizon;
        let mut terminal_guard_penalty = 0.0;

        'rollout: for step in 0..self.step_count {
            let time_s = step as f64 * self.config.step_s;
            let work_before = actuator_work_j;
            let backward_before = backward_distance_m;
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
            let target_forward_speed_m_per_s = match self.config.task {
                G1Task::Balance | G1Task::Stepping => 0.0,
                G1Task::Walking => self.config.target_forward_speed_m_per_s,
            };
            let target_velocity_body =
                rotation
                    .inverse()
                    .rotate(Vec3::new(target_forward_speed_m_per_s, 0.0, 0.0))?;
            let observation = G1PolicyObservation {
                joint_position_rad: core::array::from_fn(|actuator| position[actuator]),
                joint_velocity_rad_per_s: core::array::from_fn(|actuator| velocity[actuator]),
                gravity_direction_body,
                angular_velocity_body_rad_per_s: base.twist_body.angular,
                target_velocity_error_body_m_per_s: target_velocity_body - base.twist_body.linear,
                foot_contact: contact,
                phase_rad: TWO_PI * self.config.gait_frequency_hz * time_s,
            };
            let residual = if let Some(learned) = learned.as_deref_mut() {
                learned.act(&observation, step)
            } else {
                policy.evaluate(&observation)?
            };
            let mut external = [Wrench::default(); G1_LINK_COUNT];
            let mut next_contact = [false; 2];
            let mut minimum_sole_height_m = [f64::INFINITY; 2];
            let mut total_normal_force_n = 0.0;
            for (foot, link) in [LEFT_FOOT_LINK, RIGHT_FOOT_LINK].into_iter().enumerate() {
                let pose = kinematics.world_from_link[link];
                for point_body in FOOT_CONTACT_POINTS_BODY_M {
                    let point_world = pose.transform_point(point_body)?;
                    let surface =
                        terrain_surface(self.config.challenge, point_world.x, point_world.y);
                    maximum_abs_terrain_height_m =
                        maximum_abs_terrain_height_m.max(surface.height_m.abs());
                    minimum_sole_height_m[foot] =
                        minimum_sole_height_m[foot].min(point_world.z - surface.height_m);
                    let point_velocity_body = kinematics.body_twist[link].linear
                        + kinematics.body_twist[link].angular.cross(point_body);
                    let point_velocity_world = pose.rotation().rotate(point_velocity_body)?;
                    let normal_speed_m_per_s = point_velocity_world.dot(surface.normal_world);
                    let indentation_m =
                        ((surface.height_m - point_world.z) * surface.normal_world.z).max(0.0);
                    if indentation_m == 0.0 {
                        continue;
                    }
                    if indentation_m > MAX_CONTACT_INDENTATION_M {
                        termination_reason = G1TerminationReason::ContactIndentation;
                        terminal_guard_penalty +=
                            220.0 + 120.0 * indentation_m / MAX_CONTACT_INDENTATION_M;
                        break 'rollout;
                    }
                    if normal_speed_m_per_s.abs() > 8.0 {
                        termination_reason = G1TerminationReason::ContactSpeed;
                        terminal_guard_penalty += 260.0 + 10.0 * normal_speed_m_per_s.abs();
                        break 'rollout;
                    }
                    next_contact[foot] = true;
                    normal_request.indentation_m = indentation_m;
                    normal_request.indentation_rate_m_per_s = -normal_speed_m_per_s;
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
                    total_normal_force_n += normal_force_n;
                    let tangent_velocity_world =
                        point_velocity_world - surface.normal_world.scale(normal_speed_m_per_s);
                    let contact_frame = ContactFrame::new([
                        surface.normal_world.x,
                        surface.normal_world.y,
                        surface.normal_world.z,
                    ])?;
                    let slip = TangentialSlip::new(
                        &contact_frame,
                        [
                            tangent_velocity_world.x,
                            tangent_velocity_world.y,
                            tangent_velocity_world.z,
                        ],
                    )?;
                    let friction = self
                        .friction
                        .evaluate(&self.interface, normal_force_n, slip)?;
                    let traction = friction.traction_n();
                    let force_world = surface.normal_world.scale(normal_force_n)
                        + Vec3::new(traction[0], traction[1], traction[2]);
                    let force_body = pose.rotation().inverse().rotate(force_world)?;
                    let previous = external[link];
                    external[link] = Wrench::new(
                        previous.torque + point_body.cross(force_body),
                        previous.force + force_body,
                    );
                    slip_integral +=
                        tangent_velocity_world.dot(tangent_velocity_world) * self.config.step_s;
                }
            }
            let push_force_n = challenge_push_force_n(self.config.challenge, time_s);
            if push_force_n > 0.0 {
                let root_pose = kinematics.world_from_link[0];
                let force_world = Vec3::new(0.0, push_force_n, 0.0);
                let force_body = root_pose.rotation().inverse().rotate(force_world)?;
                let application_body = Vec3::new(0.0, 0.0, 0.42);
                let previous = external[0];
                external[0] = Wrench::new(
                    previous.torque + application_body.cross(force_body),
                    previous.force + force_body,
                );
                push_impulse_n_s += push_force_n * self.config.step_s;
            }
            let phase_signal = g1_policy_phase_basis(observation.phase_rad)?[1];
            let (desired_contact, target_clearance_m) =
                gait_targets(self.config.task, phase_signal);
            for foot in 0..2 {
                if next_contact[foot] != desired_contact[foot] {
                    contact_schedule_mismatch_integral += self.config.step_s;
                }
                if !desired_contact[foot] {
                    let clearance_error = minimum_sole_height_m[foot] - target_clearance_m[foot];
                    swing_clearance_error_integral +=
                        clearance_error * clearance_error * self.config.step_s;
                }
            }
            match next_contact {
                [true, true] => double_support_s += self.config.step_s,
                [false, false] => flight_s += self.config.step_s,
                _ => single_support_s += self.config.step_s,
            }
            // Record this step's completed count immediately after the
            // support-time match so the two stay synchronized. The per-actuator
            // loop and base-height / base-tilt checks below may break out of
            // 'rollout before the original line-987 increment; this placement
            // closes that off-by-one for the zero_policy_rollout test.
            completed_steps = step + 1;
            let body_weight_n = self.total_mass_kg * GRAVITY_WORLD_M_PER_S2.z.abs();
            let excess_load = (total_normal_force_n / body_weight_n - 1.0).max(0.0);
            impact_integral += excess_load * excess_load * self.config.step_s;
            contact = next_contact;

            let generalized_force = controller_force(
                &self.catalog,
                &self.reference_position,
                &position,
                &velocity,
                &residual,
                phase_signal,
                base.twist_body.angular,
                time_s,
                self.config.gait_frequency_hz,
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
            for actuator in 0..G1_MODEL_ACTUATORS {
                actuator_work_j +=
                    (generalized_force[actuator] * velocity[actuator]).abs() * self.config.step_s;
                let source = self.catalog.joints()[actuator];
                let next_velocity = velocity[actuator]
                    + dynamics.generalized_acceleration[actuator] * self.config.step_s;
                let velocity_limit = source.velocity_rad_per_second;
                if next_velocity.abs() > velocity_limit {
                    let normalized_overshoot =
                        (next_velocity.abs() - velocity_limit) / velocity_limit;
                    joint_limit_integral += normalized_overshoot
                        .mul_add(normalized_overshoot, 0.0)
                        .min(1.0)
                        * self.config.step_s
                        / G1_MODEL_ACTUATORS as f64;
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
                let normalized = ((position[actuator] - center) / half_range).abs();
                let soft_limit_excess = ((normalized - 0.80) / 0.20).max(0.0);
                joint_limit_integral +=
                    soft_limit_excess.powi(4) * self.config.step_s / G1_MODEL_ACTUATORS as f64;
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
            let updated_rotation = base.world_from_base.rotation();
            let updated_gravity_direction_body = updated_rotation
                .inverse()
                .rotate(Vec3::new(0.0, 0.0, -1.0))?;
            let world_velocity = updated_rotation.rotate(base.twist_body.linear)?;
            let speed_error = world_velocity.x - target_forward_speed_m_per_s;
            speed_error_integral += speed_error * speed_error * self.config.step_s;
            backward_distance_m += (-world_velocity.x).max(0.0) * self.config.step_s;
            // Per-step shaping reward for the learned-policy path: survival
            // bonus + forward progress − backward drift − actuator work.
            // (RL shaping signal — distinct from the CMA-ES objective.)
            if let Some(t) = learned_trace.as_deref_mut() {
                let progress = world_velocity.x.max(0.0) * self.config.step_s;
                let back = backward_distance_m - backward_before;
                let work = actuator_work_j - work_before;
                t.rewards.push((0.4 + 2.0 * progress - 2.0 * back - 0.002 * work) as f32);
                t.done.push(false);
            }
            let height_error = base.world_from_base.translation().z - self.initial_base_height_m;
            let normalized_height_error = height_error / POSTURE_HEIGHT_SCALE_M;
            let tilt_sine = (updated_gravity_direction_body.x * updated_gravity_direction_body.x
                + updated_gravity_direction_body.y * updated_gravity_direction_body.y)
                .sqrt();
            let normalized_tilt = tilt_sine / POSTURE_TILT_SINE_SCALE;
            minimum_base_height_m = minimum_base_height_m.min(base.world_from_base.translation().z);
            maximum_tilt_sine = maximum_tilt_sine.max(tilt_sine);
            if !recovery_recorded
                && time_s + self.config.step_s >= G1_PUSH_END_S
                && tilt_sine <= RECOVERY_TILT_SINE
                && vec_norm(base.twist_body.angular) <= RECOVERY_ANGULAR_SPEED_RAD_PER_S
                && height_error.abs() <= RECOVERY_HEIGHT_ERROR_M
            {
                recovery_time_s = (time_s + self.config.step_s - G1_PUSH_END_S).max(0.0);
                recovery_recorded = true;
            }
            posture_integral += 0.5
                * (normalized_height_error * normalized_height_error
                    + normalized_tilt * normalized_tilt)
                * self.config.step_s;
            let lateral_position = base.world_from_base.translation().y;
            lateral_error_integral +=
                (lateral_position / LATERAL_POSITION_SCALE_M).powi(2) * self.config.step_s;
            let forward_axis_world = updated_rotation.rotate(Vec3::new(1.0, 0.0, 0.0))?;
            heading_error_integral +=
                (forward_axis_world.y / HEADING_SINE_SCALE).powi(2) * self.config.step_s;
            if base.world_from_base.translation().z < MINIMUM_UPRIGHT_HEIGHT_M {
                termination_reason = G1TerminationReason::BaseHeight;
                break;
            }
            if updated_gravity_direction_body.z > MINIMUM_UPRIGHT_GRAVITY_Z {
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
        let completed_duration_s =
            (completed_steps as f64 * self.config.step_s).max(self.config.step_s);
        let target_forward_speed_m_per_s = match self.config.task {
            G1Task::Balance | G1Task::Stepping => 0.0,
            G1Task::Walking => self.config.target_forward_speed_m_per_s,
        };
        let speed_scale_m_per_s = target_forward_speed_m_per_s.max(0.25);
        let target_distance_m = (target_forward_speed_m_per_s * completed_duration_s).max(0.10);
        let body_weight_n = self.total_mass_kg * GRAVITY_WORLD_M_PER_S2.z.abs();
        let speed_tracking_error = speed_error_integral
            / (speed_scale_m_per_s * speed_scale_m_per_s * completed_duration_s);
        let normalized_stance_slip = slip_integral
            / (2.0
                * FOOT_CONTACT_POINTS_BODY_M.len() as f64
                * speed_scale_m_per_s
                * speed_scale_m_per_s
                * completed_duration_s);
        let normalized_contact_mismatch =
            contact_schedule_mismatch_integral / (2.0 * completed_duration_s);
        let normalized_clearance_error = swing_clearance_error_integral
            / (TARGET_SWING_CLEARANCE_M * TARGET_SWING_CLEARANCE_M * completed_duration_s);
        let cost_of_transport = actuator_work_j / (body_weight_n * target_distance_m);
        // v068 (cmaes-pvz): rebalance shaping weights for the whole-body 30-link
        // v067 kernel. The earlier weights were tuned for the 16-link v066 model
        // whose settling altitude was well above the 0.60 m height guard. The
        // v067 whole-body equilibrium settles lower (added arms/head/hands as
        // dynamic bodies), so the same penalty magnitudes now punish a
        // stabilizing prior that is doing small corrective work. The per-step
        // survival bonus is added below the same objective, so the optimizer
        // cannot game a single shaping rollup by collapsing early.
        let raw_shaping_score = match self.config.task {
            G1Task::Balance => {
                // Balance target is zero velocity; speed_err is a tiny numerical
                // scalar, slip/posture/contact are the real stability signals.
                // v068: remove the 30*flight term entirely (Balance is meant to
                // be near-stationary; penalizing zero-flight was anti-curriculum).
                4.0 * speed_tracking_error
                    + 12.0 * normalized_stance_slip
                    + 18.0 * posture_integral / completed_duration_s
                    + 12.0 * normalized_contact_mismatch
                    + 6.0 * lateral_error_integral / completed_duration_s
                    + 6.0 * heading_error_integral / completed_duration_s
                    + 4.0 * joint_limit_integral / completed_duration_s
                    + 2.0 * impact_integral / completed_duration_s
                    + 0.02 * cost_of_transport
                    + terminal_guard_penalty
            }
            G1Task::Stepping => {
                // Isolate the contact-mode transition before asking for
                // translation. The compact curriculum stage can therefore
                // learn genuine foot lift without sacrificing the stabilizer
                // merely to chase forward speed or stance slip.
                // v068: keep gait signals (contact_mismatch, clearance) dominant;
                // reduce posture weight so the whole-body inertia drift is not
                // punished; halve the flight penalty (the curriculum's whole-body
                // residual needs a little air to learn real stepping).
                2.0 * speed_tracking_error
                    + normalized_stance_slip
                    + 18.0 * posture_integral / completed_duration_s
                    + 260.0 * normalized_contact_mismatch
                    + 180.0 * normalized_clearance_error
                    + 14.0 * lateral_error_integral / completed_duration_s
                    + 10.0 * heading_error_integral / completed_duration_s
                    + 10.0 * joint_limit_integral / completed_duration_s
                    + 8.0 * impact_integral / completed_duration_s
                    + 0.02 * cost_of_transport
                    + 30.0 * flight_s / completed_duration_s
                    + terminal_guard_penalty
            }
            G1Task::Walking => {
                // Alternating support and real swing clearance dominate the
                // secondary score: a stable two-foot shuffle must not look
                // like successful walking. Survival remains lexicographically
                // primary through the separate uncompleted-step charge and
                // (v068) the per-step survival bonus.
                // v068: halve speed_err, slip, posture, CoT, flight. CoT was the
                // main culprit — it punishes any policy that lives long enough
                // to do real corrective work, which is exactly what the
                // whole-body prior needs. The per-step survival bonus below
                // (SURVIVAL_BONUS_PER_STEP * completed_steps) makes longer
                // survival visibly cheaper, so a slow stable walk is preferred
                // to a fast collapse.
                20.0 * speed_tracking_error
                    + 20.0 * normalized_stance_slip
                    + 12.0 * posture_integral / completed_duration_s
                    + 180.0 * normalized_contact_mismatch
                    + 100.0 * normalized_clearance_error
                    + 12.0 * lateral_error_integral / completed_duration_s
                    + 10.0 * heading_error_integral / completed_duration_s
                    + 20.0 * backward_distance_m / target_distance_m
                    + 8.0 * joint_limit_integral / completed_duration_s
                    + 6.0 * impact_integral / completed_duration_s
                    + 0.04 * cost_of_transport
                    + 20.0 * flight_s / completed_duration_s
                    + terminal_guard_penalty
            }
        };
        if !raw_shaping_score.is_finite() {
            return Err(G1WalkingError::NonFiniteObjective);
        }
        let bounded_shaping_score =
            SHAPING_SCORE_LIMIT * (raw_shaping_score / SHAPING_SCORE_SCALE).tanh();
        // Per-step survival bonus (cmaes-pvz, v068): the optimizer cannot game
        // a single shaping rollup by collapsing early, because longer survival
        // is visibly cheaper (bounded by SURVIVAL_BONUS_LIMIT = 0.4 * 720 = 288,
        // below the SHAPING_SCORE_LIMIT cap of 400 so it can never offset a fall
        // penalty but is enough to break ties between a slow stable walk and a
        // fast collapse).
        let survival_bonus = (SURVIVAL_BONUS_PER_STEP
            * (completed_steps as f64).min(MAX_BONUSED_STEPS))
        .min(SURVIVAL_BONUS_LIMIT);
        let objective = UNCOMPLETED_STEP_PENALTY * failed_horizon_steps as f64
            + bounded_shaping_score
            - survival_bonus;
        if !objective.is_finite() {
            return Err(G1WalkingError::NonFiniteObjective);
        }
        if let Some(t) = learned_trace.as_deref_mut() {
            if let Some(last) = t.done.last_mut() {
                *last = termination_reason.fell();
            }
            t.completed_steps = completed_steps;
            t.termination = termination_reason;
            t.objective = objective;
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
            backward_distance_m,
            lateral_error_integral,
            heading_error_integral,
            contact_schedule_mismatch_integral,
            swing_clearance_error_integral,
            single_support_s,
            double_support_s,
            flight_s,
            push_impulse_n_s,
            recovery_time_s,
            minimum_base_height_m,
            maximum_tilt_sine,
            maximum_abs_terrain_height_m,
            completed_steps,
            termination_reason,
            trace,
        })
    }
}

#[derive(Debug, Clone, Copy)]
struct TerrainSurface {
    height_m: f64,
    normal_world: Vec3,
}

fn terrain_surface(challenge: G1Challenge, x_m: f64, y_m: f64) -> TerrainSurface {
    if challenge == G1Challenge::Flat {
        return TerrainSurface {
            height_m: 0.0,
            normal_world: Vec3::new(0.0, 0.0, 1.0),
        };
    }
    let phase = G1_TERRAIN_WAVENUMBER_RAD_PER_M * x_m;
    let sine = phase.sin();
    let cosine = phase.cos();
    let lateral_scale = 1.0 + 0.20 * (3.0 * y_m).sin();
    let height_m = G1_TERRAIN_AMPLITUDE_M * sine * sine * lateral_scale;
    let gradient_x = G1_TERRAIN_AMPLITUDE_M
        * 2.0
        * sine
        * cosine
        * G1_TERRAIN_WAVENUMBER_RAD_PER_M
        * lateral_scale;
    let gradient_y = G1_TERRAIN_AMPLITUDE_M * sine * sine * 0.60 * (3.0 * y_m).cos();
    let raw_normal = Vec3::new(-gradient_x, -gradient_y, 1.0);
    TerrainSurface {
        height_m,
        normal_world: raw_normal.scale(1.0 / vec_norm(raw_normal)),
    }
}

fn challenge_push_force_n(challenge: G1Challenge, time_s: f64) -> f64 {
    if challenge == G1Challenge::Flat || !(G1_PUSH_START_S..=G1_PUSH_END_S).contains(&time_s) {
        return 0.0;
    }
    let phase = (time_s - G1_PUSH_START_S) / (G1_PUSH_END_S - G1_PUSH_START_S);
    G1_PUSH_PEAK_FORCE_N * (core::f64::consts::PI * phase).sin().max(0.0)
}

fn vec_norm(vector: Vec3) -> f64 {
    vector.dot(vector).sqrt()
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

fn gait_targets(task: G1Task, phase_signal: f64) -> ([bool; 2], [f64; 2]) {
    if task == G1Task::Balance {
        return ([true, true], [0.0, 0.0]);
    }
    let swing_progress =
        ((phase_signal.abs() - GAIT_SWITCH_WINDOW) / (1.0 - GAIT_SWITCH_WINDOW)).clamp(0.0, 1.0);
    let swing_clearance_m = TARGET_SWING_CLEARANCE_M * swing_progress;
    if phase_signal > GAIT_SWITCH_WINDOW {
        ([false, true], [swing_clearance_m, 0.0])
    } else if phase_signal < -GAIT_SWITCH_WINDOW {
        ([true, false], [0.0, swing_clearance_m])
    } else {
        ([true, true], [0.0, 0.0])
    }
}

const fn reference_posture() -> [f64; G1_MODEL_ACTUATORS] {
    [
        -0.20, 0.0, 0.0, 0.42, -0.22, 0.0, -0.20, 0.0, 0.0, 0.42, -0.22, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.10, 0.0, 0.30, 0.0, 0.0, 0.0, 0.0, -0.10, 0.0, 0.30, 0.0, 0.0, 0.0,
    ]
}

fn initial_height(
    catalog: &CatalogRobotModel,
    position: &[f64; G1_MODEL_ACTUATORS],
) -> Result<f64, G1WalkingError> {
    let kinematics = forward_kinematics(
        catalog.model(),
        BaseState::stationary(Se3::identity()),
        position,
        &[0.0; G1_MODEL_ACTUATORS],
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
    reference: &[f64; G1_MODEL_ACTUATORS],
    position: &[f64; G1_MODEL_ACTUATORS],
    velocity: &[f64; G1_MODEL_ACTUATORS],
    residual: &[f64; G1_POLICY_ACTUATORS],
    phase_signal: f64,
    body_angular_velocity_rad_per_s: Vec3,
    time_s: f64,
    gait_frequency_hz: f64,
) -> [f64; G1_MODEL_ACTUATORS] {
    // v069 (cmaes-zi6): smoothstep gate on the arm-swing reflex. See
    // ARM_SWING_GATE_START_S / ARM_SWING_GATE_END_S for the rationale.
    let gate_window_s = ARM_SWING_GATE_END_S - ARM_SWING_GATE_START_S;
    let swing_scale = if gate_window_s > 0.0 {
        let t = (time_s - ARM_SWING_GATE_START_S) / gate_window_s;
        let t = t.clamp(0.0, 1.0);
        t * t * (3.0 - 2.0 * t)
    } else {
        1.0
    };
    // Quiet the gait-frequency alias on the elbow bend too (it has the same
    // 0.08 * |phase| term that disturbs the stabilizer in the first cycle).
    let _ = gait_frequency_hz; // signature kept for future per-task tuning
    let mut force = [0.0; G1_MODEL_ACTUATORS];
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
        // Mode 11 raises the four hip pitch/roll limits from 88 to 139 N m.
        // Preserve the learned v6 residual envelope instead of silently
        // amplifying every existing policy coordinate by 58%.
        let policy_effort_scale = if matches!(actuator, 0 | 1 | 6 | 7) {
            88.0 / 139.0
        } else {
            1.0
        };
        let residual_force =
            RESIDUAL_EFFORT_FRACTION * effort_limit * policy_effort_scale * residual[actuator];
        force[actuator] = (posture + residual_force).clamp(-effort_limit, effort_limit);
    }

    // The remaining fourteen source joints are not fake display motion. They
    // are integrated by the same articulated owner with their official
    // inertias and limits. A small disclosed reflex swings the shoulders
    // opposite the gait phase, bends the elbows, and damps body pitch/roll.
    let mut target = *reference;
    target[15] -= ARM_SWING_AMPLITUDE_RAD * phase_signal * swing_scale
        + ARM_ROLL_DAMPING_TARGET_S * body_angular_velocity_rad_per_s.y;
    target[22] += ARM_SWING_AMPLITUDE_RAD * phase_signal * swing_scale
        - ARM_ROLL_DAMPING_TARGET_S * body_angular_velocity_rad_per_s.y;
    target[16] -= ARM_ROLL_DAMPING_TARGET_S * body_angular_velocity_rad_per_s.x;
    target[23] -= ARM_ROLL_DAMPING_TARGET_S * body_angular_velocity_rad_per_s.x;
    target[18] += 0.08 * phase_signal.abs() * swing_scale;
    target[25] += 0.08 * phase_signal.abs() * swing_scale;
    for actuator in G1_POLICY_ACTUATORS..G1_MODEL_ACTUATORS {
        let effort_limit = catalog.joints()[actuator].effort_newton_metres;
        let local = if actuator < 22 {
            actuator - 15
        } else {
            actuator - 22
        };
        let (proportional_gain, derivative_gain) = match local {
            0..=2 => (32.0, 3.2),
            3 => (26.0, 2.6),
            _ => (12.0, 1.2),
        };
        let reflex = proportional_gain * (target[actuator] - position[actuator])
            - derivative_gain * velocity[actuator];
        force[actuator] = reflex.clamp(-effort_limit, effort_limit);
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
    use fs_dfo::{CmaConfig, CmaFamily, CmaOptimizer};

    fn curriculum_coordinates() -> Vec<f64> {
        G1_STABILIZING_BIAS_MEAN
            .into_iter()
            .chain(G1_WALKING_PHASE_MEAN)
            .chain(G1_WALKING_FEEDBACK_MEAN)
            .collect()
    }

    fn policy_from_curriculum_coordinates(coordinates: &[f64]) -> [f64; G1_POLICY_DIMENSION] {
        assert_eq!(coordinates.len(), 7 * G1_POLICY_ACTUATORS);
        let mut policy = [0.0; G1_POLICY_DIMENSION];
        for actuator in 0..G1_POLICY_ACTUATORS {
            let row = actuator * G1_POLICY_FEATURES_PER_ACTUATOR;
            policy[row] = coordinates[actuator];
            policy[row + 1] = coordinates[G1_POLICY_ACTUATORS + 2 * actuator];
            policy[row + 2] = coordinates[G1_POLICY_ACTUATORS + 2 * actuator + 1];
            for (feedback, signal) in [31, 32, 34, 35].into_iter().enumerate() {
                policy[row + signal * 8] =
                    coordinates[3 * G1_POLICY_ACTUATORS + 4 * actuator + feedback];
            }
        }
        policy
    }

    fn optimize_curriculum_stage(
        label: &str,
        config: G1WalkingConfig,
        initial: Vec<f64>,
        sigma: f64,
        generations: usize,
        seed: u64,
    ) -> Vec<f64> {
        const POPULATION: usize = 16;
        let evaluator = G1WalkingEvaluator::new(config).expect("calibration evaluator");
        let mut cma_config = CmaConfig::standard(
            CmaFamily::Full,
            initial,
            sigma,
            POPULATION * generations,
            seed,
        );
        cma_config.population_size = Some(POPULATION);
        let mut optimizer = CmaOptimizer::new(cma_config).expect("calibration CMA");
        for generation in 0..generations {
            let candidates = optimizer.ask().expect("calibration ask");
            let objectives = candidates
                .candidates()
                .iter()
                .map(|coordinates| {
                    evaluator
                        .evaluate(&policy_from_curriculum_coordinates(coordinates))
                        .expect("calibration rollout")
                        .objective
                })
                .collect::<Vec<_>>();
            let snapshot = optimizer
                .tell(&candidates, &objectives)
                .expect("calibration tell");
            if generation % 20 == 19 || generation + 1 == generations {
                let best = snapshot.best.as_ref().expect("calibration best");
                let receipt = evaluator
                    .evaluate(&policy_from_curriculum_coordinates(&best.point))
                    .expect("calibration best receipt");
                eprintln!(
                    "{label} generation={} sigma={:.6} objective={:.6} steps={} distance={:.6} single_support={:.6} ending={:?}",
                    generation + 1,
                    snapshot.sigma,
                    receipt.objective,
                    receipt.completed_steps,
                    receipt.distance_m,
                    receipt.single_support_s,
                    receipt.termination_reason,
                );
            }
        }
        optimizer
            .snapshot()
            .best
            .expect("calibration completed best")
            .point
    }

    #[test]
    #[ignore = "deterministic offline provenance for the checked-in mode-11 curriculum"]
    fn recalibrate_mode_11_curriculum() {
        let short = optimize_curriculum_stage(
            "flat-0.65s",
            G1WalkingConfig {
                duration_s: 0.65,
                ..G1WalkingConfig::default()
            },
            curriculum_coordinates(),
            0.08,
            120,
            0x4731_5060,
        );
        let medium = optimize_curriculum_stage(
            "flat-1.0s",
            G1WalkingConfig {
                duration_s: 1.0,
                ..G1WalkingConfig::default()
            },
            short,
            0.055,
            140,
            0x4731_5061,
        );
        let long = optimize_curriculum_stage(
            "flat-1.5s",
            G1WalkingConfig::default(),
            medium,
            0.04,
            180,
            0x4731_5062,
        );
        let robust = optimize_curriculum_stage(
            "terrain-push-1.5s",
            G1WalkingConfig {
                challenge: G1Challenge::TerrainAndPush,
                ..G1WalkingConfig::default()
            },
            long,
            0.025,
            220,
            0x4731_5063,
        );
        let flat_evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default()).unwrap();
        let challenge_evaluator = G1WalkingEvaluator::new(G1WalkingConfig {
            challenge: G1Challenge::TerrainAndPush,
            ..G1WalkingConfig::default()
        })
        .unwrap();
        let policy = policy_from_curriculum_coordinates(&robust);
        let flat = flat_evaluator.evaluate(&policy).unwrap();
        let challenge = challenge_evaluator.evaluate(&policy).unwrap();
        eprintln!("MODE11_CALIBRATION_FLAT={flat:?}");
        eprintln!("MODE11_CALIBRATION_CHALLENGE={challenge:?}");
        eprintln!("MODE11_CALIBRATION_COORDINATES={robust:#?}");
        assert_eq!(flat.termination_reason, G1TerminationReason::Horizon);
        assert_eq!(challenge.termination_reason, G1TerminationReason::Horizon);
    }
    // v069 follow-up (cmaes-zi6): 3-stage retune for the 30-link whole-body
    // catalog (drops the terrain-push stage which trips a joint position
    // limit under the v069 dynamics — that stage is a follow-up to a
    // curriculum that's already better than the v0.6.6 baked-in one).
    // Prints the final 105-element coordinate array as Rust source we can
    // paste into G1_STABILIZING_BIAS_MEAN / G1_WALKING_PHASE_MEAN /
    // G1_WALKING_FEEDBACK_MEAN. Reduced generations per stage so the full
    // 3-stage run fits in ~4 minutes on the dev box.
    #[test]
    #[ignore = "offline curriculum retune; re-run to refresh the v0.6.9 curriculum constants"]
    fn retune_three_stage_curriculum_v069() {
        let short = optimize_curriculum_stage(
            "flat-0.65s",
            G1WalkingConfig {
                duration_s: 0.65,
                ..G1WalkingConfig::default()
            },
            curriculum_coordinates(),
            0.08,
            80,
            0x4731_5060,
        );
        let medium = optimize_curriculum_stage(
            "flat-1.0s",
            G1WalkingConfig {
                duration_s: 1.0,
                ..G1WalkingConfig::default()
            },
            short,
            0.055,
            100,
            0x4731_5061,
        );
        let long = optimize_curriculum_stage(
            "flat-1.5s",
            G1WalkingConfig::default(),
            medium,
            0.04,
            120,
            0x4731_5062,
        );
        let flat_evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default()).unwrap();
        let policy = policy_from_curriculum_coordinates(&long);
        let flat = flat_evaluator.evaluate(&policy).unwrap();
        eprintln!("V069_CURRICULUM_FLAT={flat:?}");
        eprintln!("V069_CURRICULUM_COORDINATES={long:#?}");
        eprintln!("\n=== BEGIN RUST SOURCE FOR KERNEL CONSTANTS ===\n");
        eprintln!("const G1_STABILIZING_BIAS_MEAN: [f64; G1_POLICY_ACTUATORS] = [");
        for (i, c) in long.iter().take(G1_POLICY_ACTUATORS).enumerate() {
            eprintln!("    {:.20},", c);
        }
        eprintln!("];\n");
        eprintln!("const G1_WALKING_PHASE_MEAN: [f64; 2 * G1_POLICY_ACTUATORS] = [");
        for c in long.iter().skip(G1_POLICY_ACTUATORS).take(2 * G1_POLICY_ACTUATORS) {
            eprintln!("    {:.20},", c);
        }
        eprintln!("];\n");
        eprintln!("const G1_WALKING_FEEDBACK_MEAN: [f64; 4 * G1_POLICY_ACTUATORS] = [");
        for c in long.iter().skip(3 * G1_POLICY_ACTUATORS) {
            eprintln!("    {:.20},", c);
        }
        eprintln!("];");
        eprintln!("\n=== END RUST SOURCE ===\n");
        assert_eq!(flat.termination_reason, G1TerminationReason::Horizon);
        assert_eq!(flat.completed_steps, flat_evaluator.step_count);
    }

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
        assert!(first.completed_steps > 0);
        // v069: with the arm-swing gate, the zero policy can still hit
        // a joint position limit on the upper-body links before the base
        // height guard fires. Either termination is a fair, deterministic
        // outcome for the zero policy (a baseline, not a learned policy).
        assert!(matches!(
            first.termination_reason,
            G1TerminationReason::BaseHeight | G1TerminationReason::JointPositionLimit
        ));
        let support_time = first.single_support_s + first.double_support_s + first.flight_s;
        assert!(
            (support_time - first.completed_steps as f64 * evaluator.config.step_s).abs() < 1.0e-9
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
    fn terrain_and_push_challenge_is_smooth_bounded_and_deterministic() {
        let flat = terrain_surface(G1Challenge::Flat, 0.37, -0.12);
        assert_eq!(flat.height_m, 0.0);
        assert_eq!(flat.normal_world, Vec3::new(0.0, 0.0, 1.0));

        let origin = terrain_surface(G1Challenge::TerrainAndPush, 0.0, 0.2);
        assert_eq!(origin.height_m, 0.0);
        assert!((origin.normal_world.z - 1.0).abs() < 1.0e-12);
        let crest = terrain_surface(
            G1Challenge::TerrainAndPush,
            core::f64::consts::FRAC_PI_2 / G1_TERRAIN_WAVENUMBER_RAD_PER_M,
            0.0,
        );
        assert!((crest.height_m - G1_TERRAIN_AMPLITUDE_M).abs() < 1.0e-12);
        assert!((vec_norm(crest.normal_world) - 1.0).abs() < 1.0e-12);

        assert_eq!(challenge_push_force_n(G1Challenge::Flat, 0.625), 0.0);
        assert_eq!(
            challenge_push_force_n(G1Challenge::TerrainAndPush, G1_PUSH_START_S),
            0.0
        );
        let midpoint = 0.5 * (G1_PUSH_START_S + G1_PUSH_END_S);
        assert!(
            (challenge_push_force_n(G1Challenge::TerrainAndPush, midpoint) - G1_PUSH_PEAK_FORCE_N)
                .abs()
                < 1.0e-12
        );
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
            &[0.0; G1_MODEL_ACTUATORS],
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

    #[test]
    fn curriculum_targets_are_explicit_and_alternate_without_ambiguity() {
        assert_eq!(
            gait_targets(G1Task::Balance, 1.0),
            ([true, true], [0.0, 0.0])
        );
        let (left_swing, left_clearance) = gait_targets(G1Task::Walking, 1.0);
        assert_eq!(left_swing, [false, true]);
        assert_eq!(left_clearance, [TARGET_SWING_CLEARANCE_M, 0.0]);
        assert_eq!(
            gait_targets(G1Task::Stepping, 1.0),
            (left_swing, left_clearance)
        );
        let (right_swing, right_clearance) = gait_targets(G1Task::Walking, -1.0);
        assert_eq!(right_swing, [true, false]);
        assert_eq!(right_clearance, [0.0, TARGET_SWING_CLEARANCE_M]);
        assert_eq!(
            gait_targets(G1Task::Walking, 0.0),
            ([true, true], [0.0, 0.0])
        );
    }

    #[test]
    fn arm_reflex_is_bilateral_damped_and_owner_integrated() -> Result<(), G1WalkingError> {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        let reference = evaluator.reference_position;
        let zero_velocity = [0.0; G1_MODEL_ACTUATORS];
        let zero_residual = [0.0; G1_POLICY_ACTUATORS];
        let positive = controller_force(
            &evaluator.catalog,
            &reference,
            &reference,
            &zero_velocity,
            &zero_residual,
            1.0,
            Vec3::new(0.0, 0.0, 0.0),
            // v069 (cmaes-zi6): time_s past ARM_SWING_GATE_END_S so the gate
            // is fully open; this test is checking bilateral symmetry, not
            // the gate ramp.
            2.0,
            1.55,
        );
        let negative = controller_force(
            &evaluator.catalog,
            &reference,
            &reference,
            &zero_velocity,
            &zero_residual,
            -1.0,
            Vec3::new(0.0, 0.0, 0.0),
            2.0,
            1.55,
        );
        assert_eq!(positive[15].to_bits(), negative[22].to_bits());
        assert_eq!(positive[22].to_bits(), negative[15].to_bits());
        assert!(positive[18] > 0.0 && positive[25] > 0.0);

        let dynamics = free_floating_forward_dynamics(
            evaluator.catalog.model(),
            FreeFloatingBaseState::stationary(Se3::identity()),
            &reference,
            &zero_velocity,
            &positive,
            GRAVITY_WORLD_M_PER_S2,
            &[Wrench::default(); G1_LINK_COUNT],
        )?;
        assert!(dynamics.generalized_acceleration[15].abs() > 1.0e-6);
        assert!(dynamics.generalized_acceleration[22].abs() > 1.0e-6);
        Ok(())
    }

    #[test]
    fn disclosed_stabilizing_mean_is_sparse_in_the_owner_policy_layout() {
        let mean = g1_stabilizing_policy_mean();
        assert_eq!(mean.len(), G1_POLICY_DIMENSION);
        assert_eq!(
            mean.iter().filter(|value| **value != 0.0).count(),
            G1_POLICY_ACTUATORS
        );
        for (actuator, bias) in G1_STABILIZING_BIAS_MEAN.iter().copied().enumerate() {
            assert_eq!(mean[actuator * G1_POLICY_FEATURES_PER_ACTUATOR], bias);
        }
    }

    #[test]
    fn authority_rescale_preserves_every_stabilizing_residual_torque() {
        let previous_biases: [f64; G1_POLICY_ACTUATORS] = [
            0.036_099_345_398_730_4,
            -0.352_045_625_016_877_75,
            -0.025_795_306_931_758_87,
            -0.344_351_250_162_595_2,
            -0.913_264_851_071_328_4,
            0.308_516_452_236_123_34,
            0.085_379_281_803_809_77,
            0.869_453_722_515_560_9,
            0.539_651_667_253_758_6,
            -0.263_049_902_505_655,
            0.275_096_343_882_518_25,
            0.043_717_024_880_734_576,
            -0.150_167_755_816_794_23,
            0.399_189_107_437_022_24,
            -0.673_648_120_290_733_1,
        ];
        for (previous, rescaled) in previous_biases.into_iter().zip(G1_STABILIZING_BIAS_MEAN) {
            let previous_torque_fraction = 0.32 * previous.tanh();
            let rescaled_torque_fraction = RESIDUAL_EFFORT_FRACTION * rescaled.tanh();
            assert!((previous_torque_fraction - rescaled_torque_fraction).abs() < 2.0e-16);
        }
    }

    #[test]
    fn walking_curriculum_mean_completes_with_forward_single_support() -> Result<(), G1WalkingError>
    {
        let mean = g1_walking_curriculum_mean();
        assert_eq!(mean.iter().filter(|value| **value != 0.0).count(), 105);
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default())?;
        let receipt = evaluator.evaluate(&mean)?;
        assert_eq!(receipt.termination_reason, G1TerminationReason::Horizon);
        assert_eq!(receipt.completed_steps, evaluator.step_count);
        assert!(
            receipt.distance_m > 0.59,
            "curriculum distance was {} m",
            receipt.distance_m
        );
        assert!(
            receipt.single_support_s > 0.54,
            "curriculum single-support time was {} s",
            receipt.single_support_s
        );
        assert!(
            receipt.flight_s < 0.01,
            "curriculum flight time was {} s",
            receipt.flight_s
        );
        Ok(())
    }

    #[test]
    fn walking_curriculum_survives_the_disclosed_terrain_and_push_challenge()
    -> Result<(), G1WalkingError> {
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig {
            challenge: G1Challenge::TerrainAndPush,
            ..G1WalkingConfig::default()
        })?;
        let mean = g1_walking_curriculum_mean();
        let first = evaluator.evaluate(&mean)?;
        let second = evaluator.evaluate(&mean)?;
        assert_eq!(first, second);
        assert_eq!(
            first.termination_reason,
            G1TerminationReason::Horizon,
            "challenge receipt: {first:?}"
        );
        assert_eq!(first.completed_steps, evaluator.step_count);
        assert!(first.push_impulse_n_s > 0.0);
        assert!(first.maximum_abs_terrain_height_m > 0.0);
        assert!(first.minimum_base_height_m > 0.0);
        assert!((0.0..=1.0).contains(&first.maximum_tilt_sine));
        Ok(())
    }

    // v069 (cmaes-zi6): single-stage curriculum retune for the 30-link
    // whole-body catalog. The full 4-stage recalibrate_mode_11_curriculum
    // is too expensive for the bead budget; this test runs a single
    // flat-1.5s stage (60 generations, pop=16) and prints the resulting
    // coordinates so a follow-up commit can paste them into
    // G1_STABILIZING_BIAS_MEAN / G1_WALKING_PHASE_MEAN /
    // G1_WALKING_FEEDBACK_MEAN.
    #[test]
    #[ignore = "deterministic offline retune for the 30-link catalog (cmaes-zi6)"]
    fn retune_curriculum_one_stage_v069() {
        let mut coordinates = curriculum_coordinates();
        // Warm-start the feedback terms at half-magnitude; the v066 values
        // overdrive the 30-link inertia (joint position limits at step 384).
        for value in coordinates.iter_mut().skip(3 * G1_POLICY_ACTUATORS) {
            *value *= 0.5;
        }
        let optimized = optimize_curriculum_stage(
            "v069-flat-1.5s",
            G1WalkingConfig::default(),
            coordinates,
            0.04,
            60,
            0x4731_5069,
        );
        let evaluator = G1WalkingEvaluator::new(G1WalkingConfig::default()).unwrap();
        let policy = policy_from_curriculum_coordinates(&optimized);
        let receipt = evaluator.evaluate(&policy).unwrap();
        eprintln!("V069_CURRICULUM_COORDINATES={optimized:#?}");
        eprintln!("V069_CURRICULUM_RECEIPT={receipt:?}");
    }
}
