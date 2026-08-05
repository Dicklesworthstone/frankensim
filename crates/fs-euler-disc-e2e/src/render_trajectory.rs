//! Validated animation-grade trajectory semantics for the Euler-disc pipeline.
//!
//! This module freezes the accepted public state needed by rendering and sound.
//! It deliberately does not encode or decode an artifact; canonical transport,
//! content identity, and replay belong to the later trajectory-codec layer.

use core::fmt;

use fs_blake3::ContentHash;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

use crate::coupled_runner::{
    ChannelOwnership, ChannelWrench, ContactTransitionKind, CoupledContactBranch,
    CoupledNumericalRefusalReason, CoupledRun, CoupledSample, CoupledTerminal, qois,
};

/// Exact schema version admitted by [`RenderTrajectory::try_new`].
pub const EULER_RENDER_TRAJECTORY_SCHEMA_VERSION: u16 = 1;
/// Resource ceiling for retained samples in one in-memory trajectory.
pub const MAX_RENDER_TRAJECTORY_SAMPLES: usize = 10_000_000;
/// Resource ceiling for localized transitions attached to one sample.
pub const MAX_RENDER_TRANSITIONS_PER_SAMPLE: usize = 64;
/// Resource ceiling for mandatory no-claim declarations.
pub const MAX_RENDER_TRAJECTORY_NO_CLAIMS: usize = 64;

const UNIT_TOLERANCE: f64 = 1.0e-12;
const DERIVED_QOI_TOLERANCE: f64 = 1.0e-9;
const MAX_TEXT_BYTES: usize = 1024;
const INTERVAL_END_ULP_TOLERANCE: u64 = 32;

/// Frozen v1 world-frame convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderWorldFrame {
    /// Right-handed Cartesian world coordinates with `+z` opposing gravity.
    RightHandedZUp,
    /// Deliberately non-v1 convention, useful for explicit mismatch refusal.
    RightHandedYUp,
}

/// Frozen v1 unit convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderUnitSystem {
    /// Metres, kilograms, seconds, newtons, joules, and radians.
    SiRadians,
    /// Deliberately non-v1 angular convention.
    SiDegrees,
}

/// Authority carried by this artifact class.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderTrajectoryAuthority {
    /// Accepted state and diagnostics from a declared simulation model only.
    SimulationEvidence,
}

/// Binding between exact mass properties and their upstream identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderMassProperties {
    /// Content identity of the resolved mass-property artifact.
    pub identity: ContentHash,
    /// Exact properties used to interpret body-frame angular momentum.
    pub properties: MassProperties,
}

/// Declares which reduced-model channel payloads are present in every sample.
///
/// A present channel may legitimately be all zero. An unavailable channel must
/// be all zero, so consumers never have to guess whether zero means quiescent
/// physics or missing data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RenderChannelAvailability {
    /// Gravity wrench/work payloads are present.
    pub gravity: bool,
    /// Aggregate contact wrench/work payloads are present.
    pub contact: bool,
    /// Reduced rolling-resistance wrench/work payloads are present.
    pub rolling: bool,
    /// Reduced-base damping-channel payloads are present.
    pub base: bool,
    /// Exterior-gas body wrench/work payloads are present.
    pub gas: bool,
}

/// Nominal rigid frame of the reduced base plane.
///
/// The one-mode displacement translates this frame along its local `+z`
/// without changing orientation. V1 requires that axis to coincide with world
/// `+z`, while retaining yaw and origin so local contact coordinates survive
/// admissible horizontal rigid transforms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderBaseFrame {
    /// Nominal origin of the undisplaced base frame in world coordinates [m].
    pub origin_world_m: Vec3,
    /// Nominal base-to-world orientation.
    pub orientation_base_to_world: UnitQuaternion,
}

impl RenderChannelAvailability {
    /// Availability emitted by the current coupled runner, which evaluates
    /// every declared channel even when a configured coefficient is zero.
    pub const ALL_AVAILABLE: Self = Self {
        gravity: true,
        contact: true,
        rolling: true,
        base: true,
        gas: true,
    };

    /// Explicitly unavailable channel set for import/refusal tests.
    pub const NONE_AVAILABLE: Self = Self {
        gravity: false,
        contact: false,
        rolling: false,
        base: false,
        gas: false,
    };
}

/// Top-level metadata required to interpret every retained sample.
#[derive(Clone, Debug, PartialEq)]
pub struct RenderTrajectoryMetadata {
    /// Exact schema version.
    pub schema_version: u16,
    /// World frame repeated by every sample to prevent silent cross-wiring.
    pub world_frame: RenderWorldFrame,
    /// Unit system repeated by every sample to prevent silent cross-wiring.
    pub units: RenderUnitSystem,
    /// Resolved specimen/profile identity used for support geometry.
    pub specimen_profile_identity: ContentHash,
    /// Exact chart or contact-geometry identity.
    pub specimen_chart_identity: ContentHash,
    /// Mass properties and their identity.
    pub mass_properties: RenderMassProperties,
    /// Accepted initial rigid-body state.
    pub initial_state: RigidBodyState,
    /// Accepted initial one-mode base state.
    pub initial_base_mode: RenderBaseModeState,
    /// Base model/configuration identity.
    pub base_model_identity: ContentHash,
    /// Nominal coordinate frame of the reduced base plane.
    pub base_frame: RenderBaseFrame,
    /// Physics model identity.
    pub model_identity: ContentHash,
    /// Explicit availability of every interval wrench/work channel.
    pub channel_availability: RenderChannelAvailability,
    /// Full run configuration identity.
    pub configuration_identity: ContentHash,
    /// Existing reduced-run restart fingerprint, retained as an opaque audit aid.
    pub configuration_fingerprint: u64,
    /// Declared fixed macro timestep in seconds.
    pub timestep_s: f64,
    /// Producer implementation/version label.
    pub producer_version: String,
    /// Bounded Context-of-Use applicability statement.
    pub applicability: String,
    /// Explicit statements that constrain downstream interpretation.
    pub no_claims: Vec<String>,
    /// Authority ceiling of the trajectory.
    pub authority: RenderTrajectoryAuthority,
}

/// Derived Euler-disc quantities retained as auditable conveniences.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DerivedEulerQois {
    /// Angle between the symmetry axis and world `+z` [rad].
    pub inclination_rad: f64,
    /// Azimuthal symmetry-axis rate about world `+z` [rad/s].
    pub precession_rad_per_s: f64,
    /// Residual axial spin after removing precession [rad/s].
    pub spin_rad_per_s: f64,
    /// Finite-difference precession acceleration [rad/s^2].
    pub precession_acceleration_rad_per_s2: f64,
}

impl DerivedEulerQois {
    /// Derive the three intrinsic pose/twist QoIs from authoritative state.
    pub fn from_state(
        state: RigidBodyState,
        mass: MassProperties,
        precession_acceleration_rad_per_s2: f64,
    ) -> Result<Self, RenderTrajectoryError> {
        if !precession_acceleration_rad_per_s2.is_finite() {
            return Err(RenderTrajectoryError::NonFinite {
                sample: None,
                field: "qois.precession_acceleration_rad_per_s2",
            });
        }
        let (inclination_rad, precession_rad_per_s, spin_rad_per_s) = qois(state, mass)
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
        Ok(Self {
            inclination_rad,
            precession_rad_per_s,
            spin_rad_per_s,
            precession_acceleration_rad_per_s2,
        })
    }
}

/// Stable support-feature vocabulary for contact visualization.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderSupportFeature {
    /// Compatibility cylinder rim without a resolved profile index.
    CylinderRim,
    /// Feature index from the resolved profile support query.
    ProfileFeature(usize),
}

/// Contact point and normal associated with a closed contact branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderContactGeometry {
    /// World-space contact point [m].
    pub point_world_m: Vec3,
    /// Unit normal expressed in the world frame and directed base-to-disc.
    pub normal_world: Vec3,
    /// Geometry feature selected by the support query.
    pub support_feature: RenderSupportFeature,
}

/// Unilateral branch at the retained post-step state.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderContactBranch {
    /// Positive/open-gap branch.
    Open,
    /// Closed signed-gap branch; zero force remains possible at a localized root.
    Closed,
}

/// One localized branch transition retained inside the preceding interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderContactTransition {
    /// Opening or reimpact.
    pub kind: ContactTransitionKind,
    /// Localized event time [s].
    pub time_s: f64,
    /// Inclusive bracket start [s].
    pub bracket_start_s: f64,
    /// Inclusive bracket end [s].
    pub bracket_end_s: f64,
}

/// Complete one-mode base state required by the current reduced runner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderBaseModeState {
    /// Generalized vertical displacement [m].
    pub displacement_m: f64,
    /// Generalized vertical velocity [m/s].
    pub velocity_m_per_s: f64,
}

/// Localized inclination-threshold event attached only to a terminal sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RenderTerminalEvent {
    /// Localized threshold time [s].
    pub time_s: f64,
    /// Inclusive root bracket start [s].
    pub bracket_start_s: f64,
    /// Inclusive root bracket end [s].
    pub bracket_end_s: f64,
}

/// Why a retained sample ends the public trajectory.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderSampleDisposition {
    /// More accepted samples follow.
    Continue,
    /// The declared inclination threshold was localized.
    TerminalInclination,
    /// The declared time/step horizon was reached before a physical terminal.
    HorizonCensored,
    /// The numerical model refused to continue.
    NumericalRefusal(RenderNumericalRefusalReason),
}

/// Stable cross-backend reason carried by a numerical-refusal terminal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderNumericalRefusalReason {
    /// The declared separation/reimpact budget was exceeded.
    ReimpactLimitExceeded,
    /// A bracketed contact event could not make deterministic progress.
    ContactEventLocalizationFailed,
    /// Energy accounting or the reduced base state became non-finite.
    NonFiniteEnergyOrBaseState,
    /// Nonzero producer-specific code retained without widening its authority.
    BackendSpecific(u32),
}

/// Raw construction input. Validation converts its quaternion to the canonical
/// `q/-q` representative and retains a checked [`RigidBodyState`].
#[derive(Clone, Debug, PartialEq)]
pub struct RenderTrajectorySampleInput {
    /// Exact start of the accepted interval ending at `time_s` [s]. For a
    /// resumed segment this is the restart checkpoint time, not zero.
    pub interval_start_time_s: f64,
    /// Exact producer `f64` time value [s]; later encoding preserves its bits.
    pub time_s: f64,
    /// Repeated frame declaration.
    pub world_frame: RenderWorldFrame,
    /// Repeated unit declaration.
    pub units: RenderUnitSystem,
    /// Center-of-mass position in the world frame [m].
    pub center_of_mass_world_m: Vec3,
    /// `(w,x,y,z)` body-to-world quaternion; must already be unit length.
    pub orientation_body_to_world: [f64; 4],
    /// World-frame linear momentum [kg m/s].
    pub linear_momentum_world_kg_m_per_s: Vec3,
    /// Principal-body-frame angular momentum [kg m^2/s].
    pub angular_momentum_body_kg_m2_per_s: Vec3,
    /// Redundant symmetry axis in the world frame, checked against orientation.
    pub symmetry_axis_world: Vec3,
    /// Post-step unilateral branch.
    pub contact_branch: RenderContactBranch,
    /// Required exactly for a closed branch.
    pub contact_geometry: Option<RenderContactGeometry>,
    /// Signed support gap [m].
    pub signed_gap_m: f64,
    /// Whether any accepted subinterval used the closed contact branch. This
    /// is not inferred from force magnitude because a localized root may have
    /// zero penetration and zero force.
    pub interval_contact_active: bool,
    /// Normal load evaluated at the midpoint of the first accepted subinterval
    /// [N]. When an event splits the interval this is not its mean load; use the
    /// duration-weighted contact-channel force for interval controls.
    pub interval_normal_force_n: f64,
    /// Chronological localized transitions in the preceding interval.
    pub contact_transitions: Vec<RenderContactTransition>,
    /// Complete reduced-base state; omission is not silently interpreted as zero.
    pub base_mode: Option<RenderBaseModeState>,
    /// Per-channel wrench/work accounting.
    pub channels: ChannelOwnership,
    /// Total declared mechanical energy [J].
    pub mechanical_energy_j: f64,
    /// Declared energy-closure defect [J].
    pub energy_defect_j: f64,
    /// Redundant, independently checked Euler QoIs.
    pub qois: DerivedEulerQois,
    /// Terminal/censor state.
    pub disposition: RenderSampleDisposition,
    /// Required exactly when `disposition` is `TerminalInclination`.
    pub terminal_event: Option<RenderTerminalEvent>,
}

/// Admitted animation sample with canonical orientation and checked state.
#[derive(Clone, Debug, PartialEq)]
pub struct RenderTrajectorySample {
    input: RenderTrajectorySampleInput,
    state: RigidBodyState,
}

impl RenderTrajectorySample {
    /// Original semantic fields with canonicalized quaternion components.
    #[must_use]
    pub const fn input(&self) -> &RenderTrajectorySampleInput {
        &self.input
    }

    /// Checked accepted rigid-body state.
    #[must_use]
    pub const fn state(&self) -> RigidBodyState {
        self.state
    }
}

/// Validated public trajectory. It is simulation evidence, not calibrated truth.
#[derive(Clone, Debug, PartialEq)]
pub struct RenderTrajectory {
    metadata: RenderTrajectoryMetadata,
    samples: Vec<RenderTrajectorySample>,
}

impl RenderTrajectory {
    /// Convert the complete accepted states published by the reduced coupled
    /// runner. Provisional stages are never visible here. The final retained
    /// sample must match the runner checkpoint exactly; otherwise no public
    /// trajectory is minted.
    pub fn from_coupled_run(
        metadata: RenderTrajectoryMetadata,
        run: &CoupledRun,
    ) -> Result<Self, RenderTrajectoryError> {
        if metadata.configuration_fingerprint != run.checkpoint.configuration_fingerprint
            || metadata.mass_properties.properties != run.mass_properties
            || metadata.channel_availability != RenderChannelAvailability::ALL_AVAILABLE
            || metadata.base_frame.origin_world_m != Vec3::ZERO
            || metadata.base_frame.orientation_base_to_world != UnitQuaternion::IDENTITY
            || metadata.initial_state != run.configuration_initial_state
            || metadata.initial_base_mode.displacement_m.to_bits()
                != run.configuration_initial_base_deflection_m.to_bits()
            || metadata.initial_base_mode.velocity_m_per_s.to_bits()
                != run.configuration_initial_base_velocity_m_per_s.to_bits()
            || metadata.timestep_s.to_bits() != run.macro_timestep_s.to_bits()
        {
            return Err(RenderTrajectoryError::RunnerConfigurationMismatch);
        }
        let Some(last) = run.samples.last() else {
            return Err(RenderTrajectoryError::EmptyTrajectory);
        };
        if last.time_s.to_bits() != run.checkpoint.time_s.to_bits()
            || last.state != run.checkpoint.state
            || last.base_deflection_m.to_bits() != run.checkpoint.base_deflection_m.to_bits()
            || last.base_velocity_m_per_s.to_bits()
                != run.checkpoint.base_velocity_m_per_s.to_bits()
            || last.energy_defect_j.to_bits()
                != run.checkpoint.accumulated_energy_defect_j.to_bits()
        {
            return Err(RenderTrajectoryError::RunnerCheckpointMismatch);
        }
        for (index, sample) in run.samples.iter().enumerate() {
            let derived_velocity = sample
                .state
                .center_of_mass_velocity_world(run.mass_properties)
                .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
            if derived_velocity != sample.center_of_mass_velocity_world_m_per_s {
                return Err(RenderTrajectoryError::RunnerSampleMismatch(index));
            }
        }
        let final_index = run.samples.len() - 1;
        let inputs = run
            .samples
            .iter()
            .enumerate()
            .map(|(index, sample)| {
                coupled_sample_input(
                    sample,
                    if index == final_index {
                        render_disposition(run.terminal)
                    } else {
                        RenderSampleDisposition::Continue
                    },
                )
            })
            .collect();
        Self::try_new(metadata, inputs)
    }

    /// Validate metadata, all samples, cross-sample time/event rules, and the
    /// final terminal/censor placement.
    pub fn try_new(
        metadata: RenderTrajectoryMetadata,
        inputs: Vec<RenderTrajectorySampleInput>,
    ) -> Result<Self, RenderTrajectoryError> {
        validate_metadata(&metadata)?;
        if inputs.is_empty() {
            return Err(RenderTrajectoryError::EmptyTrajectory);
        }
        if inputs.len() > MAX_RENDER_TRAJECTORY_SAMPLES {
            return Err(RenderTrajectoryError::TooManySamples(inputs.len()));
        }

        let sample_count = inputs.len();
        let mut samples = Vec::new();
        samples.try_reserve_exact(sample_count).map_err(|_| {
            RenderTrajectoryError::Capacity {
                artifact: "render trajectory samples",
                requested: sample_count,
            }
        })?;
        let mut previous_time = None;
        let mut previous_branch = None;
        for (index, input) in inputs.into_iter().enumerate() {
            let sample = validate_sample(&metadata, input, index, previous_time)?;
            validate_transition_origin(&sample.input, previous_branch, index)?;
            previous_time = Some(sample.input.time_s);
            previous_branch = Some(sample.input.contact_branch);
            let is_last = index + 1 == sample_count;
            if !is_last && sample.input.disposition != RenderSampleDisposition::Continue {
                return Err(RenderTrajectoryError::TerminalBeforeFinalSample(index));
            }
            samples.push(sample);
        }
        if samples
            .last()
            .is_some_and(|sample| sample.input.disposition == RenderSampleDisposition::Continue)
        {
            return Err(RenderTrajectoryError::MissingFinalDisposition);
        }
        Ok(Self { metadata, samples })
    }

    /// Interpretation metadata.
    #[must_use]
    pub const fn metadata(&self) -> &RenderTrajectoryMetadata {
        &self.metadata
    }

    /// Strictly time-ordered admitted samples.
    #[must_use]
    pub fn samples(&self) -> &[RenderTrajectorySample] {
        &self.samples
    }
}

/// Structured refusal from trajectory admission.
#[derive(Clone, Debug, PartialEq)]
pub enum RenderTrajectoryError {
    /// Only the exact v1 schema is currently accepted.
    UnsupportedSchemaVersion(u16),
    /// Only right-handed `+z`-up coordinates are admitted by v1.
    UnsupportedWorldFrame(RenderWorldFrame),
    /// Only SI/radians are admitted by v1.
    UnsupportedUnits(RenderUnitSystem),
    /// The reduced base frame was non-finite or its local `+z` did not align
    /// with the frozen world `+z` axis.
    UnsupportedBaseFrame,
    /// A mandatory content identity was all zeroes.
    ZeroIdentity(&'static str),
    /// A bounded required text field was empty or too large.
    InvalidText(&'static str),
    /// The declared timestep was not finite and positive.
    InvalidTimestep,
    /// No samples were supplied.
    EmptyTrajectory,
    /// The retained sample ceiling was exceeded.
    TooManySamples(usize),
    /// A bounded trajectory allocation could not be admitted.
    Capacity {
        /// Stable artifact/component name.
        artifact: &'static str,
        /// Requested item count.
        requested: usize,
    },
    /// A sample repeated a different frame.
    FrameMismatch(usize),
    /// A sample repeated a different unit system.
    UnitMismatch(usize),
    /// A scalar or vector was non-finite.
    NonFinite {
        /// Sample index, or `None` for metadata/standalone derivation.
        sample: Option<usize>,
        /// Stable semantic field name.
        field: &'static str,
    },
    /// The raw quaternion was not already unit length.
    QuaternionNotUnit(usize),
    /// A checked rigid-body state could not be formed.
    InvalidRigidState(usize, String),
    /// Redundant symmetry axis disagreed with canonical orientation.
    SymmetryAxisMismatch(usize),
    /// Contact branch and optional geometry disagreed.
    ContactGeometryMismatch(usize),
    /// Closed-contact normal was not unit length.
    ContactNormalNotUnit(usize),
    /// Signed gap contradicted the declared contact branch.
    ContactGapMismatch(usize),
    /// Contact force was negative.
    NegativeNormalForce(usize),
    /// An interval declared no active contact but carried contact-only data.
    InactiveContactHasIntervalData(usize),
    /// A localized transition or its bracket was malformed.
    InvalidTransition {
        /// Sample containing the invalid transition.
        sample: usize,
        /// Transition index within the sample.
        transition: usize,
    },
    /// Localized transitions do not alternate or do not end on the declared branch.
    ContactTransitionBranchMismatch(usize),
    /// Required reduced-base displacement/velocity was absent.
    MissingBaseState(usize),
    /// Sample times were duplicated or decreased.
    NonMonotoneTime(usize),
    /// The exact interval start was non-finite, after its endpoint, or did not
    /// equal the preceding retained endpoint.
    InvalidIntervalStart(usize),
    /// A positive interval exceeded the declared fixed macro timestep beyond
    /// the binary64 comparison tolerance.
    IntervalExceedsDeclaredTimestep(usize),
    /// A zero-duration initial point carried interval-only data.
    ZeroDurationIntervalData(usize),
    /// A channel declared unavailable carried a nonzero payload.
    UnavailableChannelHasData {
        /// Sample containing the contradictory payload.
        sample: usize,
        /// Stable channel name.
        channel: &'static str,
    },
    /// Redundant QoIs did not agree with authoritative state/mass properties.
    DerivedQoiMismatch(usize),
    /// Authoritative-state QoI derivation refused.
    DerivedState(String),
    /// A terminal/censor flag appeared before the final sample.
    TerminalBeforeFinalSample(usize),
    /// Terminal disposition and localized threshold event disagreed.
    TerminalEventMismatch(usize),
    /// A backend-specific numerical refusal used the reserved zero code.
    InvalidNumericalRefusalCode(usize),
    /// The final sample did not explain why retention ended.
    MissingFinalDisposition,
    /// Trajectory metadata did not name the runner's admitted configuration.
    RunnerConfigurationMismatch,
    /// The final accepted sample did not exactly match the restart checkpoint.
    RunnerCheckpointMismatch,
    /// A runner sample's redundant accepted-state diagnostic was inconsistent.
    RunnerSampleMismatch(usize),
}

impl fmt::Display for RenderTrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RenderTrajectoryError {}

fn coupled_sample_input(
    sample: &CoupledSample,
    disposition: RenderSampleDisposition,
) -> RenderTrajectorySampleInput {
    let orientation = sample.state.pose().orientation();
    let contact_geometry =
        (sample.contact_branch == CoupledContactBranch::Closed).then_some(RenderContactGeometry {
            point_world_m: sample.endpoint_contact_geometry.point_world_m,
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            support_feature: sample.support_source_feature.map_or(
                RenderSupportFeature::CylinderRim,
                RenderSupportFeature::ProfileFeature,
            ),
        });
    RenderTrajectorySampleInput {
        interval_start_time_s: sample.interval_start_time_s,
        time_s: sample.time_s,
        world_frame: RenderWorldFrame::RightHandedZUp,
        units: RenderUnitSystem::SiRadians,
        center_of_mass_world_m: sample.state.pose().position_world(),
        orientation_body_to_world: orientation.components(),
        linear_momentum_world_kg_m_per_s: sample.state.linear_momentum_world(),
        angular_momentum_body_kg_m2_per_s: sample.state.angular_momentum_body(),
        symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        contact_branch: match sample.contact_branch {
            CoupledContactBranch::Open => RenderContactBranch::Open,
            CoupledContactBranch::Closed => RenderContactBranch::Closed,
        },
        contact_geometry,
        signed_gap_m: sample.endpoint_signed_gap_m,
        interval_contact_active: sample.contact_active,
        interval_normal_force_n: sample.interval_normal_force_n,
        contact_transitions: sample
            .contact_transitions
            .iter()
            .map(|transition| RenderContactTransition {
                kind: transition.kind,
                time_s: transition.time_s,
                bracket_start_s: transition.bracket_start_s,
                bracket_end_s: transition.bracket_end_s,
            })
            .collect(),
        base_mode: Some(RenderBaseModeState {
            displacement_m: sample.base_deflection_m,
            velocity_m_per_s: sample.base_velocity_m_per_s,
        }),
        channels: sample.channels,
        mechanical_energy_j: sample.mechanical_energy_j,
        energy_defect_j: sample.energy_defect_j,
        qois: DerivedEulerQois {
            inclination_rad: sample.inclination_rad,
            precession_rad_per_s: sample.precession_rad_per_s,
            spin_rad_per_s: sample.spin_rad_per_s,
            precession_acceleration_rad_per_s2: sample.precession_acceleration_rad_per_s2,
        },
        disposition,
        terminal_event: sample
            .terminal_inclination_event
            .map(|event| RenderTerminalEvent {
                time_s: event.time_s,
                bracket_start_s: event.bracket_start_s,
                bracket_end_s: event.bracket_end_s,
            }),
    }
}

const fn render_disposition(terminal: CoupledTerminal) -> RenderSampleDisposition {
    match terminal {
        CoupledTerminal::TerminalInclination => RenderSampleDisposition::TerminalInclination,
        CoupledTerminal::HorizonReached => RenderSampleDisposition::HorizonCensored,
        CoupledTerminal::NumericalRefusal { reason } => {
            RenderSampleDisposition::NumericalRefusal(match reason {
                CoupledNumericalRefusalReason::ReimpactLimitExceeded => {
                    RenderNumericalRefusalReason::ReimpactLimitExceeded
                }
                CoupledNumericalRefusalReason::ContactEventLocalizationFailed => {
                    RenderNumericalRefusalReason::ContactEventLocalizationFailed
                }
                CoupledNumericalRefusalReason::NonFiniteEnergyOrBaseState => {
                    RenderNumericalRefusalReason::NonFiniteEnergyOrBaseState
                }
            })
        }
    }
}

fn validate_metadata(metadata: &RenderTrajectoryMetadata) -> Result<(), RenderTrajectoryError> {
    if metadata.schema_version != EULER_RENDER_TRAJECTORY_SCHEMA_VERSION {
        return Err(RenderTrajectoryError::UnsupportedSchemaVersion(
            metadata.schema_version,
        ));
    }
    if metadata.world_frame != RenderWorldFrame::RightHandedZUp {
        return Err(RenderTrajectoryError::UnsupportedWorldFrame(
            metadata.world_frame,
        ));
    }
    if metadata.units != RenderUnitSystem::SiRadians {
        return Err(RenderTrajectoryError::UnsupportedUnits(metadata.units));
    }
    for (name, identity) in [
        (
            "specimen_profile_identity",
            metadata.specimen_profile_identity,
        ),
        ("specimen_chart_identity", metadata.specimen_chart_identity),
        (
            "mass_properties.identity",
            metadata.mass_properties.identity,
        ),
        ("base_model_identity", metadata.base_model_identity),
        ("model_identity", metadata.model_identity),
        ("configuration_identity", metadata.configuration_identity),
    ] {
        if identity.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(RenderTrajectoryError::ZeroIdentity(name));
        }
    }
    if !metadata.timestep_s.is_finite() || metadata.timestep_s <= 0.0 {
        return Err(RenderTrajectoryError::InvalidTimestep);
    }
    if !metadata.initial_base_mode.displacement_m.is_finite()
        || !metadata.initial_base_mode.velocity_m_per_s.is_finite()
    {
        return Err(RenderTrajectoryError::NonFinite {
            sample: None,
            field: "initial_base_mode",
        });
    }
    if !metadata.base_frame.origin_world_m.is_finite()
        || !vec_close(
            metadata
                .base_frame
                .orientation_base_to_world
                .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
            Vec3::new(0.0, 0.0, 1.0),
            UNIT_TOLERANCE,
        )
    {
        return Err(RenderTrajectoryError::UnsupportedBaseFrame);
    }
    metadata
        .initial_state
        .center_of_mass_velocity_world(metadata.mass_properties.properties)
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    metadata
        .mass_properties
        .properties
        .angular_velocity_body_checked(metadata.initial_state.angular_momentum_body())
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    validate_text("producer_version", &metadata.producer_version)?;
    validate_text("applicability", &metadata.applicability)?;
    if metadata.no_claims.is_empty() || metadata.no_claims.len() > MAX_RENDER_TRAJECTORY_NO_CLAIMS {
        return Err(RenderTrajectoryError::InvalidText("no_claims"));
    }
    for no_claim in &metadata.no_claims {
        validate_text("no_claim", no_claim)?;
    }
    Ok(())
}

fn validate_sample(
    metadata: &RenderTrajectoryMetadata,
    mut input: RenderTrajectorySampleInput,
    index: usize,
    previous_time: Option<f64>,
) -> Result<RenderTrajectorySample, RenderTrajectoryError> {
    if input.world_frame != metadata.world_frame {
        return Err(RenderTrajectoryError::FrameMismatch(index));
    }
    if input.units != metadata.units {
        return Err(RenderTrajectoryError::UnitMismatch(index));
    }
    finite_scalar(input.time_s, index, "time_s")?;
    if input.time_s < 0.0 || previous_time.is_some_and(|time| input.time_s <= time) {
        return Err(RenderTrajectoryError::NonMonotoneTime(index));
    }
    finite_scalar(input.interval_start_time_s, index, "interval_start_time_s")?;
    if input.interval_start_time_s < 0.0
        || input.interval_start_time_s > input.time_s
        || previous_time.is_some_and(|time| input.interval_start_time_s.to_bits() != time.to_bits())
    {
        return Err(RenderTrajectoryError::InvalidIntervalStart(index));
    }
    let declared_end_s = input.interval_start_time_s + metadata.timestep_s;
    let maximum_end_s = advance_nonnegative_ulps(declared_end_s, INTERVAL_END_ULP_TOLERANCE);
    if input.time_s > maximum_end_s {
        return Err(RenderTrajectoryError::IntervalExceedsDeclaredTimestep(
            index,
        ));
    }
    finite_vec(
        input.center_of_mass_world_m,
        index,
        "center_of_mass_world_m",
    )?;
    finite_vec(
        input.linear_momentum_world_kg_m_per_s,
        index,
        "linear_momentum_world_kg_m_per_s",
    )?;
    finite_vec(
        input.angular_momentum_body_kg_m2_per_s,
        index,
        "angular_momentum_body_kg_m2_per_s",
    )?;
    let orientation = checked_unit_quaternion(input.orientation_body_to_world, index)?;
    input.orientation_body_to_world = orientation.components();
    let pose = Pose::new(input.center_of_mass_world_m, orientation)
        .map_err(|error| RenderTrajectoryError::InvalidRigidState(index, error.to_string()))?;
    let state = RigidBodyState::new(
        pose,
        input.linear_momentum_world_kg_m_per_s,
        input.angular_momentum_body_kg_m2_per_s,
    )
    .map_err(|error| RenderTrajectoryError::InvalidRigidState(index, error.to_string()))?;

    finite_vec(input.symmetry_axis_world, index, "symmetry_axis_world")?;
    let expected_axis = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    if !unit_vector(input.symmetry_axis_world)
        || !vec_close(input.symmetry_axis_world, expected_axis, UNIT_TOLERANCE)
    {
        return Err(RenderTrajectoryError::SymmetryAxisMismatch(index));
    }
    validate_channels(&input.channels, metadata.channel_availability, index)?;
    validate_contact(metadata, &input, index)?;
    let base = input
        .base_mode
        .ok_or(RenderTrajectoryError::MissingBaseState(index))?;
    finite_scalar(base.displacement_m, index, "base_mode.displacement_m")?;
    finite_scalar(base.velocity_m_per_s, index, "base_mode.velocity_m_per_s")?;
    if input.interval_start_time_s.to_bits() == input.time_s.to_bits()
        && (input.interval_contact_active
            || input.interval_normal_force_n != 0.0
            || !input.contact_transitions.is_empty()
            || channel_ownership_has_data(&input.channels))
    {
        return Err(RenderTrajectoryError::ZeroDurationIntervalData(index));
    }
    finite_scalar(input.mechanical_energy_j, index, "mechanical_energy_j")?;
    finite_scalar(input.energy_defect_j, index, "energy_defect_j")?;
    validate_qois(metadata, &input, state, index)?;
    validate_transitions(
        &input.contact_transitions,
        input.contact_branch,
        input.disposition,
        input.interval_start_time_s,
        input.time_s,
        metadata.timestep_s,
        index,
    )?;
    validate_terminal_event(
        &input,
        input.interval_start_time_s,
        metadata.timestep_s,
        index,
    )?;

    Ok(RenderTrajectorySample { input, state })
}

fn checked_unit_quaternion(
    components: [f64; 4],
    index: usize,
) -> Result<UnitQuaternion, RenderTrajectoryError> {
    if components.iter().any(|component| !component.is_finite()) {
        return Err(RenderTrajectoryError::NonFinite {
            sample: Some(index),
            field: "orientation_body_to_world",
        });
    }
    let scale = components
        .iter()
        .fold(0.0_f64, |maximum, component| maximum.max(component.abs()));
    if scale == 0.0 {
        return Err(RenderTrajectoryError::QuaternionNotUnit(index));
    }
    let scaled_norm = components
        .iter()
        .map(|component| component / scale)
        .map(|component| component * component)
        .sum::<f64>()
        .sqrt();
    let norm = scale * scaled_norm;
    if !norm.is_finite() || (norm - 1.0).abs() > UNIT_TOLERANCE {
        return Err(RenderTrajectoryError::QuaternionNotUnit(index));
    }
    UnitQuaternion::new(components[0], components[1], components[2], components[3])
        .map_err(|error| RenderTrajectoryError::InvalidRigidState(index, error.to_string()))
}

fn advance_nonnegative_ulps(value: f64, ulps: u64) -> f64 {
    if value.is_infinite() {
        value
    } else {
        f64::from_bits(
            value
                .to_bits()
                .saturating_add(ulps)
                .min(f64::INFINITY.to_bits()),
        )
    }
}

fn validate_contact(
    metadata: &RenderTrajectoryMetadata,
    input: &RenderTrajectorySampleInput,
    index: usize,
) -> Result<(), RenderTrajectoryError> {
    finite_scalar(input.signed_gap_m, index, "signed_gap_m")?;
    finite_scalar(
        input.interval_normal_force_n,
        index,
        "interval_normal_force_n",
    )?;
    if input.interval_normal_force_n < 0.0 {
        return Err(RenderTrajectoryError::NegativeNormalForce(index));
    }
    if metadata.channel_availability.contact {
        let base_axis_world = metadata
            .base_frame
            .orientation_base_to_world
            .rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        if input.channels.contact.force_world_n.dot(base_axis_world) < 0.0 {
            return Err(RenderTrajectoryError::NegativeNormalForce(index));
        }
    }
    if !input.interval_contact_active
        && (input.interval_normal_force_n != 0.0
            || channel_has_data(input.channels.contact)
            || input
                .contact_transitions
                .iter()
                .any(|event| event.kind == ContactTransitionKind::Opening))
    {
        return Err(RenderTrajectoryError::InactiveContactHasIntervalData(index));
    }
    match (input.contact_branch, input.contact_geometry) {
        (RenderContactBranch::Open, None) => {
            if input.signed_gap_m < 0.0 {
                return Err(RenderTrajectoryError::ContactGapMismatch(index));
            }
        }
        (RenderContactBranch::Closed, Some(geometry)) => {
            if input.signed_gap_m > 0.0 {
                return Err(RenderTrajectoryError::ContactGapMismatch(index));
            }
            finite_vec(geometry.point_world_m, index, "contact.point_world_m")?;
            finite_vec(geometry.normal_world, index, "contact.normal_world")?;
            if !unit_vector(geometry.normal_world) {
                return Err(RenderTrajectoryError::ContactNormalNotUnit(index));
            }
        }
        _ => return Err(RenderTrajectoryError::ContactGeometryMismatch(index)),
    }
    let positive_duration = input.interval_start_time_s.to_bits() != input.time_s.to_bits();
    let inactive_closed_at_endpoint_reimpact = input.contact_branch == RenderContactBranch::Closed
        && input.contact_transitions.last().is_some_and(|event| {
            event.kind == ContactTransitionKind::Reimpact
                && event.time_s.to_bits() == input.time_s.to_bits()
        });
    if positive_duration
        && !input.interval_contact_active
        && input.contact_branch == RenderContactBranch::Closed
        && !inactive_closed_at_endpoint_reimpact
    {
        return Err(RenderTrajectoryError::InactiveContactHasIntervalData(index));
    }
    if positive_duration
        && input.interval_contact_active
        && input.contact_branch == RenderContactBranch::Open
        && !input
            .contact_transitions
            .iter()
            .any(|event| event.kind == ContactTransitionKind::Opening)
    {
        return Err(RenderTrajectoryError::ContactTransitionBranchMismatch(
            index,
        ));
    }
    Ok(())
}

fn validate_transitions(
    transitions: &[RenderContactTransition],
    final_branch: RenderContactBranch,
    disposition: RenderSampleDisposition,
    interval_start: f64,
    sample_time: f64,
    timestep_s: f64,
    sample: usize,
) -> Result<(), RenderTrajectoryError> {
    if transitions.len() > MAX_RENDER_TRANSITIONS_PER_SAMPLE {
        return Err(RenderTrajectoryError::InvalidTransition {
            sample,
            transition: MAX_RENDER_TRANSITIONS_PER_SAMPLE,
        });
    }
    let mut last_time = interval_start;
    let mut previous_kind = None;
    for (transition, event) in transitions.iter().enumerate() {
        let terminal_root_bracket = transition + 1 == transitions.len()
            && disposition
                == RenderSampleDisposition::NumericalRefusal(
                    RenderNumericalRefusalReason::ReimpactLimitExceeded,
                )
            && event.kind == ContactTransitionKind::Reimpact
            && event.time_s.to_bits() == sample_time.to_bits()
            && event.bracket_end_s - sample_time <= timestep_s;
        let valid = event.time_s.is_finite()
            && event.bracket_start_s.is_finite()
            && event.bracket_end_s.is_finite()
            && event.bracket_start_s <= event.time_s
            && event.time_s <= event.bracket_end_s
            && interval_start <= event.bracket_start_s
            && (event.bracket_end_s <= sample_time || terminal_root_bracket)
            && (transition == 0 && event.time_s >= last_time
                || transition > 0 && event.time_s > last_time)
            && previous_kind != Some(event.kind);
        if !valid {
            return Err(RenderTrajectoryError::InvalidTransition { sample, transition });
        }
        last_time = event.time_s;
        previous_kind = Some(event.kind);
    }
    if let Some(final_kind) = previous_kind {
        let final_branch_matches = matches!(
            (final_kind, final_branch),
            (ContactTransitionKind::Opening, RenderContactBranch::Open)
                | (ContactTransitionKind::Reimpact, RenderContactBranch::Closed)
        );
        if !final_branch_matches {
            return Err(RenderTrajectoryError::ContactTransitionBranchMismatch(
                sample,
            ));
        }
    }
    if disposition
        == RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        )
        && !transitions.last().is_some_and(|event| {
            event.kind == ContactTransitionKind::Reimpact
                && event.time_s.to_bits() == sample_time.to_bits()
        })
    {
        return Err(RenderTrajectoryError::InvalidTransition {
            sample,
            transition: transitions.len().saturating_sub(1),
        });
    }
    Ok(())
}

fn validate_transition_origin(
    input: &RenderTrajectorySampleInput,
    previous_branch: Option<RenderContactBranch>,
    sample: usize,
) -> Result<(), RenderTrajectoryError> {
    let Some(previous_branch) = previous_branch else {
        // The retained segment-start branch is unavailable for first-sample
        // preroll, so its transition origin cannot be reconstructed safely.
        return Ok(());
    };
    let origin_matches = input.contact_transitions.first().map_or(
        input.contact_branch == previous_branch,
        |transition| {
            matches!(
                (previous_branch, transition.kind),
                (RenderContactBranch::Open, ContactTransitionKind::Reimpact)
                    | (RenderContactBranch::Closed, ContactTransitionKind::Opening)
            )
        },
    );
    if origin_matches {
        Ok(())
    } else {
        Err(RenderTrajectoryError::ContactTransitionBranchMismatch(
            sample,
        ))
    }
}

fn validate_terminal_event(
    input: &RenderTrajectorySampleInput,
    interval_start: f64,
    timestep_s: f64,
    index: usize,
) -> Result<(), RenderTrajectoryError> {
    if matches!(
        input.disposition,
        RenderSampleDisposition::NumericalRefusal(RenderNumericalRefusalReason::BackendSpecific(0))
    ) {
        return Err(RenderTrajectoryError::InvalidNumericalRefusalCode(index));
    }
    match (input.disposition, input.terminal_event) {
        (RenderSampleDisposition::TerminalInclination, Some(event)) => {
            let retained_root_bracket = event.time_s.to_bits() == input.time_s.to_bits()
                && event.bracket_end_s - input.time_s <= timestep_s;
            if !event.time_s.is_finite()
                || !event.bracket_start_s.is_finite()
                || !event.bracket_end_s.is_finite()
                || event.time_s.to_bits() != input.time_s.to_bits()
                || event.bracket_start_s > event.time_s
                || event.time_s > event.bracket_end_s
                || event.bracket_start_s < interval_start
                || (event.bracket_end_s > input.time_s && !retained_root_bracket)
            {
                return Err(RenderTrajectoryError::TerminalEventMismatch(index));
            }
        }
        (RenderSampleDisposition::TerminalInclination, None)
        | (RenderSampleDisposition::Continue, Some(_))
        | (RenderSampleDisposition::HorizonCensored, Some(_))
        | (RenderSampleDisposition::NumericalRefusal(_), Some(_)) => {
            return Err(RenderTrajectoryError::TerminalEventMismatch(index));
        }
        _ => {}
    }
    Ok(())
}

fn validate_qois(
    metadata: &RenderTrajectoryMetadata,
    input: &RenderTrajectorySampleInput,
    state: RigidBodyState,
    index: usize,
) -> Result<(), RenderTrajectoryError> {
    for (field, value) in [
        ("qois.inclination_rad", input.qois.inclination_rad),
        ("qois.precession_rad_per_s", input.qois.precession_rad_per_s),
        ("qois.spin_rad_per_s", input.qois.spin_rad_per_s),
        (
            "qois.precession_acceleration_rad_per_s2",
            input.qois.precession_acceleration_rad_per_s2,
        ),
    ] {
        finite_scalar(value, index, field)?;
    }
    let derived = DerivedEulerQois::from_state(
        state,
        metadata.mass_properties.properties,
        input.qois.precession_acceleration_rad_per_s2,
    )?;
    if !scalar_close(input.qois.inclination_rad, derived.inclination_rad)
        || !scalar_close(
            input.qois.precession_rad_per_s,
            derived.precession_rad_per_s,
        )
        || !scalar_close(input.qois.spin_rad_per_s, derived.spin_rad_per_s)
    {
        return Err(RenderTrajectoryError::DerivedQoiMismatch(index));
    }
    Ok(())
}

fn validate_channels(
    channels: &ChannelOwnership,
    availability: RenderChannelAvailability,
    index: usize,
) -> Result<(), RenderTrajectoryError> {
    for (name, channel, available) in [
        ("gravity", channels.gravity, availability.gravity),
        ("contact", channels.contact, availability.contact),
        ("rolling", channels.rolling, availability.rolling),
        ("base", channels.base, availability.base),
        ("gas", channels.gas, availability.gas),
    ] {
        finite_vec(channel.force_world_n, index, "channels.force_world_n")?;
        finite_vec(channel.torque_world_nm, index, "channels.torque_world_nm")?;
        finite_scalar(channel.work_j, index, "channels.work_j")?;
        if !available && channel_has_data(channel) {
            return Err(RenderTrajectoryError::UnavailableChannelHasData {
                sample: index,
                channel: name,
            });
        }
    }
    Ok(())
}

fn channel_ownership_has_data(channels: &ChannelOwnership) -> bool {
    [
        channels.gravity,
        channels.contact,
        channels.rolling,
        channels.base,
        channels.gas,
    ]
    .into_iter()
    .any(channel_has_data)
}

fn channel_has_data(channel: ChannelWrench) -> bool {
    channel.force_world_n != Vec3::ZERO
        || channel.torque_world_nm != Vec3::ZERO
        || channel.work_j != 0.0
}

fn validate_text(field: &'static str, value: &str) -> Result<(), RenderTrajectoryError> {
    if value.trim().is_empty() || value.len() > MAX_TEXT_BYTES || value.contains('\0') {
        return Err(RenderTrajectoryError::InvalidText(field));
    }
    Ok(())
}

fn finite_scalar(
    value: f64,
    sample: usize,
    field: &'static str,
) -> Result<(), RenderTrajectoryError> {
    if !value.is_finite() {
        return Err(RenderTrajectoryError::NonFinite {
            sample: Some(sample),
            field,
        });
    }
    Ok(())
}

fn finite_vec(
    value: Vec3,
    sample: usize,
    field: &'static str,
) -> Result<(), RenderTrajectoryError> {
    if !value.is_finite() {
        return Err(RenderTrajectoryError::NonFinite {
            sample: Some(sample),
            field,
        });
    }
    Ok(())
}

fn unit_vector(value: Vec3) -> bool {
    value.is_finite() && (value.norm_squared() - 1.0).abs() <= UNIT_TOLERANCE
}

fn vec_close(left: Vec3, right: Vec3, tolerance: f64) -> bool {
    (left.x - right.x).abs() <= tolerance
        && (left.y - right.y).abs() <= tolerance
        && (left.z - right.z).abs() <= tolerance
}

fn scalar_close(left: f64, right: f64) -> bool {
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    (left - right).abs() <= DERIVED_QOI_TOLERANCE * scale
}
