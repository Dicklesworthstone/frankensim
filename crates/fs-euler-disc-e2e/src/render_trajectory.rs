//! Validated animation-grade trajectory semantics for the Euler-disc pipeline.
//!
//! This module freezes the accepted public state needed by rendering and sound.
//! It deliberately does not encode or decode an artifact; canonical transport,
//! content identity, and replay belong to the later trajectory-codec layer.

use core::{f64::consts::TAU, fmt};

use fs_blake3::{ContentHash, DomainHasher, hash_domain};
use fs_exec::Cx;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

use crate::contact_dynamics::{
    profile_contact_geometry, profile_mass_to_mbd, profile_state_at_ground_contact,
};
use crate::coupled_runner::{
    ChannelOwnership, ChannelWrench, ContactTransitionKind, CoupledContactBranch,
    CoupledNumericalRefusalReason, CoupledRun, CoupledSample, CoupledTerminal, qois,
};
use crate::production_coupling::{
    ProductionControlTrajectory, ProductionCouplingError, ProductionCouplingModel,
    ProductionCouplingReceipt, ProductionEventTrajectory, ProductionEventTrajectoryTermination,
    ProductionOpenFlightReceipt, ProductionTrajectoryBranch, ProductionTrajectoryStepReceipt,
    SmoothContactTrajectory, SmoothContactTrajectoryTermination,
};
use crate::reduced_decay::{
    BILDSTEN_PUBLISHED_POWER_COEFFICIENT, ChannelPowers, ChannelWork, REDUCED_DECAY_MODEL_ID,
    ReducedDecayRun, ReducedDecaySample, ReducedDecayTerminal,
    THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S, THORNE_2026_FITTED_ROLLING_COEFFICIENT,
    THORNE_2026_SOURCE_ID,
};
use crate::specimen::ResolvedDiscProfile;

/// Exact schema version admitted by [`RenderTrajectory::try_new`].
pub const EULER_RENDER_TRAJECTORY_SCHEMA_VERSION: u16 = 3;
/// Resource ceiling for retained samples in one in-memory trajectory.
pub const MAX_RENDER_TRAJECTORY_SAMPLES: usize = 10_000_000;
/// Resource ceiling for localized transitions attached to one sample.
pub const MAX_RENDER_TRANSITIONS_PER_SAMPLE: usize = 64;
/// Resource ceiling for mandatory no-claim declarations.
pub const MAX_RENDER_TRAJECTORY_NO_CLAIMS: usize = 64;
/// Exact deterministic convention for the reduced-decay render bridge.
pub const REDUCED_DECAY_RENDER_BRIDGE_VERSION: u32 = 3;
/// Exact producer version for the transactional production-prefix bridge.
pub const PRODUCTION_COUPLING_RENDER_BRIDGE_VERSION: u32 = 1;
/// Exact producer version for the event-aware production trajectory bridge.
pub const PRODUCTION_EVENT_RENDER_BRIDGE_VERSION: u32 = 1;
/// Cinematic tail retained by the default reduced-decay bridge [s].
///
/// The source run is slightly longer than eight seconds. Rebasing its final
/// eight seconds to `t = 0` keeps the positive-cutoff chirp inside an
/// eight-second film without changing coefficients or time-scaling motion.
pub const REDUCED_DECAY_RENDER_TAIL_HORIZON_S: f64 = 8.0;

const UNIT_TOLERANCE: f64 = 1.0e-12;
const DERIVED_QOI_TOLERANCE: f64 = 1.0e-9;
const MAX_TEXT_BYTES: usize = 1024;
const INTERVAL_END_ULP_TOLERANCE: u64 = 32;
// A tail bridge may rebase a producer clock by subtracting its source-time
// origin from each endpoint.  At sub-second output times that subtraction can
// leave a few femtoseconds of cancellation residue, far larger than 32 ULPs
// of a `1e-4 s` endpoint but still only bounded binary64 rounding noise.  Do
// not extend this allowance to large clocks: their existing endpoint-ULP
// bound remains the authority, including near `f64::MAX`.
const SUBSECOND_REBASED_CLOCK_TOLERANCE_S: f64 = INTERVAL_END_ULP_TOLERANCE as f64 * f64::EPSILON;
const TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES: usize = 1_024;
const REDUCED_DECAY_MODEL_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-disc.reduced-decay-render-model.v2";
const REDUCED_DECAY_CONFIGURATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-disc.reduced-decay-render-configuration.v2";
const REDUCED_DECAY_BASE_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-disc.reduced-decay-static-base.v1";
const REDUCED_DECAY_PHASE_CONVENTION: &str = "q=Rz(precession)*Ry(theta)*Rz(spin); phi_dot=Omega; psi_dot=-Omega*cos(theta); theta_dot=-Phi/dE_dtheta; trapezoidal-phase-v1";
const REDUCED_DECAY_CONTACT_REACTION_CONVENTION: &str = "positive-duration intervals carry N_bar=m*(g+(v_z1-v_z0)/dt) as a separately available base-normal load; full contact wrench/work remains unavailable; initial point carries no interval data; no resolved torque, tangential traction, subinterval force history, or contact patch";

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
    /// Exact sampling rule for the independently retained base-normal scalar.
    ///
    /// This is intentionally not a boolean: a first accepted-subinterval
    /// midpoint is useful contact diagnostics, but it is not a duration mean
    /// and must not be used as an acoustic force measure.
    pub normal_force_sampling: RenderNormalForceSampling,
    /// Reduced rolling-resistance wrench/work payloads are present.
    pub rolling: bool,
    /// Reduced-base damping-channel payloads are present.
    pub base: bool,
    /// Exterior-gas body wrench/work payloads are present.
    pub gas: bool,
}

/// Exact sampling rule for [`RenderTrajectorySampleInput::interval_normal_force_n`].
///
/// The explicit tag prevents consumers from guessing semantics from the
/// presence of full contact channels. In particular, a normal-only midpoint
/// can never be promoted to an interval mean for sound synthesis.
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RenderNormalForceSampling {
    /// No normal-load scalar is authoritative; its encoded value must be zero.
    Unavailable = 0,
    /// The scalar is the representative normal reaction at the first accepted
    /// contact subinterval midpoint. It is diagnostic-only unless a full
    /// contact wrench provides a separate duration mean.
    FirstAcceptedSubintervalMidpoint = 1,
    /// The scalar is the mean base-normal reaction over the accepted interval.
    IntervalMean = 2,
    /// The scalar is the normal-law evaluation applied as a constant world
    /// force throughout one accepted transactional mechanics substep. It is
    /// therefore the exact interval mean of the discretized zero-order-hold
    /// forcing, while remaining only a timestep-dependent approximation to the
    /// continuously evolving physical force.
    AppliedSubstepZeroOrderHold = 3,
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
        normal_force_sampling: RenderNormalForceSampling::FirstAcceptedSubintervalMidpoint,
        rolling: true,
        base: true,
        gas: true,
    };

    /// Explicitly unavailable channel set for import/refusal tests.
    pub const NONE_AVAILABLE: Self = Self {
        gravity: false,
        contact: false,
        normal_force_sampling: RenderNormalForceSampling::Unavailable,
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
    /// Precession acceleration derived from the admitted motion law [rad/s^2].
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
    /// Producer-declared normal load for the preceding interval [N]. Its exact
    /// sampling rule is `metadata.channel_availability.normal_force_sampling`.
    /// Consumers must never infer interval-mean authority from this scalar's
    /// magnitude or from unrelated channel availability.
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

    /// Convert an accepted smooth-contact prefix from the transactional
    /// production mechanics stack into the common render trajectory.
    ///
    /// Only fields established by that backend are published. The exact
    /// contact point and single-step normal-law evaluation are retained, while
    /// the render channel-work payloads remain explicitly unavailable because
    /// the current production receipt does not yet publish one common
    /// cross-channel work ledger. The mechanical-energy field is the complete
    /// disc-only `fs-mbd` diagnostic, and `energy_defect_j` is its constant-
    /// wrench midpoint work residual. Neither is mislabeled as full coupled
    /// disc/base/contact energy closure.
    pub fn from_production_coupling_prefix(
        model: &ProductionCouplingModel,
        source: &SmoothContactTrajectory,
        profile: &ResolvedDiscProfile,
        declared_maximum_timestep_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        cx.checkpoint()
            .map_err(|_| RenderTrajectoryError::Cancelled)?;
        if source.accepted_steps.is_empty() {
            return Err(RenderTrajectoryError::ProductionPrefixEmpty);
        }
        if !(declared_maximum_timestep_s.is_finite() && declared_maximum_timestep_s > 0.0) {
            return Err(RenderTrajectoryError::InvalidTimestep);
        }
        let profile_mass = profile_mass_to_mbd(profile.mass_properties)
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
        model
            .validate_checkpoint(&source.start_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        model
            .validate_checkpoint(&source.last_accepted_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        if model.disc_mass_properties != profile_mass {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }
        let expected_version = source
            .start_checkpoint
            .committed_version
            .checked_add(u64::try_from(source.accepted_steps.len()).map_err(|_| {
                RenderTrajectoryError::Capacity {
                    artifact: "production-prefix version",
                    requested: source.accepted_steps.len(),
                }
            })?)
            .ok_or(RenderTrajectoryError::Capacity {
                artifact: "production-prefix version",
                requested: source.accepted_steps.len(),
            })?;
        if source.last_accepted_checkpoint.committed_version != expected_version {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let mut inputs = Vec::new();
        inputs
            .try_reserve_exact(source.accepted_steps.len())
            .map_err(|_| RenderTrajectoryError::Capacity {
                artifact: "production render-prefix samples",
                requested: source.accepted_steps.len(),
            })?;
        let mut interval_start_time_s = source.start_checkpoint.elapsed_time_s();
        let mut previous_state = source.start_checkpoint.disc_state;
        let mut previous_base_version = source.start_checkpoint.committed_version;
        let mut last_base_displacement_m = source.start_checkpoint.base_displacement_m();
        let mut last_base_velocity_m_per_s = source.start_checkpoint.base_velocity_m_per_s();
        let final_index = source.accepted_steps.len() - 1;
        for (index, receipt) in source.accepted_steps.iter().enumerate() {
            if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
                cx.checkpoint()
                    .map_err(|_| RenderTrajectoryError::Cancelled)?;
            }
            let base = receipt.base.receipt();
            let duration_s = receipt.rigid_step.duration_seconds;
            if duration_s.to_bits() != base.timestep_s.to_bits()
                || duration_s > declared_maximum_timestep_s
                || receipt.rigid_step.state_before != previous_state
                || receipt.rigid_step.state_after != receipt.next_disc_state
                || base.parent_version != previous_base_version
                || previous_base_version.checked_add(1) != Some(base.next_version)
            {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "accepted state/clock lineage",
                });
            }
            let time_s = interval_start_time_s + duration_s;
            if !time_s.is_finite() {
                return Err(RenderTrajectoryError::NonFinite {
                    sample: Some(index),
                    field: "production-prefix time",
                });
            }
            let nominal_normal_force_n = match &receipt.normal.generic.receipt {
                fs_contact::normal_patch::NormalPatchReceipt::Point(point) => point.normal_force_n,
                fs_contact::normal_patch::NormalPatchReceipt::Line(_) => {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "line contact cannot produce a point-resultant trajectory",
                    });
                }
            };
            let expected_applied_normal_force_n = nominal_normal_force_n
                + receipt
                    .surface_excitation
                    .as_ref()
                    .map_or(0.0, |surface| surface.normal_force_perturbation_n);
            let normal_force_n = base.compressive_normal_force_on_base_n;
            if expected_applied_normal_force_n.to_bits() != normal_force_n.to_bits() {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "topography-perturbed normal-force lineage",
                });
            }
            // The force receipt is evaluated from the interval-start patch,
            // but render contact geometry is an endpoint field. Re-query the
            // actual profile support at the accepted post-step pose instead
            // of pairing a stale start point with the new rigid state.
            let endpoint_contact = profile_contact_geometry(
                &profile.chart,
                profile.mass_properties,
                receipt.next_disc_state.pose(),
                cx,
            )
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
            let normal_world = Vec3::new(0.0, 0.0, 1.0);
            let signed_gap_m = endpoint_contact.contact.gap_m - base.modal_displacement_end_m;
            if signed_gap_m > 0.0 {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "accepted smooth-contact step ends after an unlocalized opening",
                });
            }
            let before_qois = DerivedEulerQois::from_state(
                receipt.rigid_step.state_before,
                model.disc_mass_properties,
                0.0,
            )?;
            let after_qois_without_acceleration = DerivedEulerQois::from_state(
                receipt.rigid_step.state_after,
                model.disc_mass_properties,
                0.0,
            )?;
            let precession_acceleration_rad_per_s2 = (after_qois_without_acceleration
                .precession_rad_per_s
                - before_qois.precession_rad_per_s)
                / duration_s;
            let qois = DerivedEulerQois {
                precession_acceleration_rad_per_s2,
                ..after_qois_without_acceleration
            };
            let energy_defect_j = production_disc_work_residual_j(model, receipt, index)?;
            let orientation = receipt.next_disc_state.pose().orientation();
            inputs.push(RenderTrajectorySampleInput {
                interval_start_time_s,
                time_s,
                world_frame: RenderWorldFrame::RightHandedZUp,
                units: RenderUnitSystem::SiRadians,
                center_of_mass_world_m: receipt.next_disc_state.pose().position_world(),
                orientation_body_to_world: orientation.components(),
                linear_momentum_world_kg_m_per_s: receipt.next_disc_state.linear_momentum_world(),
                angular_momentum_body_kg_m2_per_s: receipt.next_disc_state.angular_momentum_body(),
                symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
                contact_branch: RenderContactBranch::Closed,
                contact_geometry: Some(RenderContactGeometry {
                    point_world_m: endpoint_contact.contact.point_world_m,
                    normal_world,
                    support_feature: RenderSupportFeature::ProfileFeature(
                        endpoint_contact.support_source_feature,
                    ),
                }),
                signed_gap_m,
                interval_contact_active: true,
                interval_normal_force_n: normal_force_n,
                contact_transitions: Vec::new(),
                base_mode: Some(RenderBaseModeState {
                    displacement_m: base.modal_displacement_end_m,
                    velocity_m_per_s: base.modal_velocity_end_m_per_s,
                }),
                channels: ChannelOwnership::default(),
                mechanical_energy_j: receipt.rigid_step.diagnostics_after.mechanical_energy,
                energy_defect_j,
                qois,
                disposition: if index == final_index {
                    production_prefix_disposition(&source.termination)
                } else {
                    RenderSampleDisposition::Continue
                },
                terminal_event: None,
            });
            interval_start_time_s = time_s;
            previous_state = receipt.next_disc_state;
            previous_base_version = base.next_version;
            last_base_displacement_m = base.modal_displacement_end_m;
            last_base_velocity_m_per_s = base.modal_velocity_end_m_per_s;
        }
        if previous_state != source.last_accepted_checkpoint.disc_state
            || previous_base_version != source.last_accepted_checkpoint.committed_version
            || interval_start_time_s.to_bits()
                != source.last_accepted_checkpoint.elapsed_time_s().to_bits()
            || last_base_displacement_m.to_bits()
                != source
                    .last_accepted_checkpoint
                    .base_displacement_m()
                    .to_bits()
            || last_base_velocity_m_per_s.to_bits()
                != source
                    .last_accepted_checkpoint
                    .base_velocity_m_per_s()
                    .to_bits()
        {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let identities = profile.content_identities();
        let model_identity = production_model_identity(model, identities.profile);
        let configuration_identity = production_configuration_identity(model, model_identity);
        let mut fingerprint = [0_u8; 8];
        fingerprint.copy_from_slice(&configuration_identity.as_bytes()[..8]);
        let (base_model_id, base_configuration_id) = model.base_port.identity_parts();
        let base_model_identity = hash_domain(
            "org.frankensim.euler-disc.production-base-model.v1",
            format!("{base_model_id}\0{base_configuration_id}").as_bytes(),
        );
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: identities.chart,
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: model.disc_mass_properties,
            },
            initial_state: source.start_checkpoint.disc_state,
            initial_base_mode: RenderBaseModeState {
                displacement_m: source.start_checkpoint.base_displacement_m(),
                velocity_m_per_s: source.start_checkpoint.base_velocity_m_per_s(),
            },
            base_model_identity,
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity,
            channel_availability: RenderChannelAvailability {
                gravity: false,
                contact: false,
                normal_force_sampling: RenderNormalForceSampling::AppliedSubstepZeroOrderHold,
                rolling: false,
                base: false,
                gas: false,
            },
            configuration_identity,
            configuration_fingerprint: u64::from_le_bytes(fingerprint),
            timestep_s: declared_maximum_timestep_s,
            producer_version: format!(
                "fs-euler-disc-e2e/production-coupling-render-bridge-v{}",
                PRODUCTION_COUPLING_RENDER_BRIDGE_VERSION
            ),
            applicability: "accepted smooth-contact prefix from the transactional finite-patch/partial-slip/rolling/gas/one-mode-base composition; ends at its explicit step budget or first typed mechanism refusal".to_owned(),
            no_claims: vec![
                "Estimate-authority simulation prefix; not experimental calibration or validated Euler-disc stopping-time prediction".to_owned(),
                "smooth-contact prefix only; separation, impact, reimpact, and terminal continuation are not synthesized".to_owned(),
                "normal force is the exact mean of the discretized applied zero-order hold; physical force bandwidth and convergence remain limited by the mechanics timestep".to_owned(),
                "render channel work is unavailable pending one shared cross-channel work ledger; zero payloads mean unavailable, not zero physical loss".to_owned(),
                "mechanical energy and defect are disc-only fs-mbd diagnostics, not full disc/base/contact/gas energy closure".to_owned(),
                "model identity binds the declared production configuration and exposed adapter/base identities; it does not independently introspect every private constituent law coefficient".to_owned(),
                "no thermal evolution, phase change, melting, plastic flow, structural-mode solve, or acoustic-radiation solve is implied".to_owned(),
            ],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        Self::try_new(metadata, inputs)
    }

    /// Convert an accepted event-aware production trajectory into the common
    /// render and sound state stream.
    ///
    /// The source interval branch controls whether a physical contact force was
    /// applied. A branch change detected at the following checkpoint is attached
    /// to the preceding interval endpoint with its full fixed-grid time bracket:
    /// contact-to-open therefore retains the contact interval force and ends
    /// open, while open-to-contact retains exactly zero normal force and ends
    /// closed. This is the state convention required by the common trajectory
    /// validator and prevents an impact-only audio impulse from being invented.
    pub fn from_production_event_trajectory(
        model: &ProductionCouplingModel,
        source: &ProductionEventTrajectory,
        profile: &ResolvedDiscProfile,
        declared_maximum_timestep_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        cx.checkpoint()
            .map_err(|_| RenderTrajectoryError::Cancelled)?;
        if source.accepted_steps.is_empty() {
            return Err(RenderTrajectoryError::ProductionPrefixEmpty);
        }
        if !(declared_maximum_timestep_s.is_finite() && declared_maximum_timestep_s > 0.0) {
            return Err(RenderTrajectoryError::InvalidTimestep);
        }
        let profile_mass = profile_mass_to_mbd(profile.mass_properties)
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
        model
            .validate_checkpoint(&source.start_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        model
            .validate_checkpoint(&source.last_accepted_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        if model.disc_mass_properties != profile_mass {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }
        let expected_version = source
            .start_checkpoint
            .committed_version
            .checked_add(u64::try_from(source.accepted_steps.len()).map_err(|_| {
                RenderTrajectoryError::Capacity {
                    artifact: "production-event version",
                    requested: source.accepted_steps.len(),
                }
            })?)
            .ok_or(RenderTrajectoryError::Capacity {
                artifact: "production-event version",
                requested: source.accepted_steps.len(),
            })?;
        if source.last_accepted_checkpoint.committed_version != expected_version {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let mut inputs = Vec::new();
        inputs
            .try_reserve_exact(source.accepted_steps.len())
            .map_err(|_| RenderTrajectoryError::Capacity {
                artifact: "production event render samples",
                requested: source.accepted_steps.len(),
            })?;
        let mut interval_start_time_s = source.start_checkpoint.elapsed_time_s();
        let mut previous_state = source.start_checkpoint.disc_state;
        let mut previous_base_version = source.start_checkpoint.committed_version;
        let mut last_base_displacement_m = source.start_checkpoint.base_displacement_m();
        let mut last_base_velocity_m_per_s = source.start_checkpoint.base_velocity_m_per_s();
        let mut transition_index = 0_usize;
        let final_index = source.accepted_steps.len() - 1;

        for (index, step) in source.accepted_steps.iter().enumerate() {
            if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
                cx.checkpoint()
                    .map_err(|_| RenderTrajectoryError::Cancelled)?;
            }
            if step.start_time_s.to_bits() != interval_start_time_s.to_bits() {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "event source start clock",
                });
            }
            let (rigid_step, base, contact_receipt) = match (&step.branch, &step.receipt) {
                (
                    ProductionTrajectoryBranch::CompliantContact,
                    ProductionTrajectoryStepReceipt::CompliantContact(receipt),
                ) => (&receipt.rigid_step, receipt.base.receipt(), Some(receipt)),
                (
                    ProductionTrajectoryBranch::OpenFlight,
                    ProductionTrajectoryStepReceipt::OpenFlight(receipt),
                ) => (&receipt.rigid_step, receipt.base.receipt(), None),
                _ => {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "event branch/receipt mismatch",
                    });
                }
            };
            let duration_s = rigid_step.duration_seconds;
            let time_s = interval_start_time_s + duration_s;
            if duration_s.to_bits() != base.timestep_s.to_bits()
                || duration_s > declared_maximum_timestep_s
                || step.end_time_s.to_bits() != time_s.to_bits()
                || rigid_step.state_before != previous_state
                || base.parent_version != previous_base_version
                || previous_base_version.checked_add(1) != Some(base.next_version)
            {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "event accepted state/clock lineage",
                });
            }

            let (next_disc_state, normal_force_n, energy_defect_j) = if let Some(receipt) =
                contact_receipt
            {
                if rigid_step.state_after != receipt.next_disc_state {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "contact endpoint state",
                    });
                }
                let nominal_normal_force_n = match &receipt.normal.generic.receipt {
                    fs_contact::normal_patch::NormalPatchReceipt::Point(point) => {
                        point.normal_force_n
                    }
                    fs_contact::normal_patch::NormalPatchReceipt::Line(_) => {
                        return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                            sample: index,
                            field: "line contact cannot produce a point-resultant trajectory",
                        });
                    }
                };
                let expected_force_n = nominal_normal_force_n
                    + receipt
                        .surface_excitation
                        .as_ref()
                        .map_or(0.0, |surface| surface.normal_force_perturbation_n);
                if expected_force_n.to_bits() != base.compressive_normal_force_on_base_n.to_bits() {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "event topography-perturbed normal-force lineage",
                    });
                }
                (
                    receipt.next_disc_state,
                    base.compressive_normal_force_on_base_n,
                    production_disc_work_residual_j(model, receipt, index)?,
                )
            } else {
                let ProductionTrajectoryStepReceipt::OpenFlight(receipt) = &step.receipt else {
                    unreachable!("branch/receipt match checked above")
                };
                if rigid_step.state_after != receipt.next_disc_state
                    || base.compressive_normal_force_on_base_n != 0.0
                {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "open endpoint state or zero contact load",
                    });
                }
                (
                    receipt.next_disc_state,
                    0.0,
                    production_open_disc_work_residual_j(model, receipt, index)?,
                )
            };

            let transition = source
                .transitions
                .get(transition_index)
                .filter(|transition| transition.bracket_end_s.to_bits() == time_s.to_bits());
            let endpoint_branch = transition.map_or(step.branch, |transition| transition.to);
            if let Some(transition) = transition {
                if transition.from != step.branch
                    || transition.bracket_start_s.to_bits() != interval_start_time_s.to_bits()
                    || !(transition.bracket_start_s <= transition.bracket_end_s)
                {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "event transition lineage",
                    });
                }
                transition_index += 1;
            }
            let endpoint_contact = profile_contact_geometry(
                &profile.chart,
                profile.mass_properties,
                next_disc_state.pose(),
                cx,
            )
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
            let signed_gap_m = endpoint_contact.contact.gap_m - base.modal_displacement_end_m;
            let contact_geometry = match endpoint_branch {
                ProductionTrajectoryBranch::CompliantContact => {
                    if signed_gap_m > 0.0 {
                        return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                            sample: index,
                            field: "closed event endpoint has positive gap",
                        });
                    }
                    Some(RenderContactGeometry {
                        point_world_m: endpoint_contact.contact.point_world_m,
                        normal_world: Vec3::new(0.0, 0.0, 1.0),
                        support_feature: RenderSupportFeature::ProfileFeature(
                            endpoint_contact.support_source_feature,
                        ),
                    })
                }
                ProductionTrajectoryBranch::OpenFlight => {
                    if signed_gap_m < 0.0 {
                        return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                            sample: index,
                            field: "open event endpoint has negative gap",
                        });
                    }
                    None
                }
            };
            let contact_transitions = transition
                .map(|transition| {
                    vec![RenderContactTransition {
                        kind: match (transition.from, transition.to) {
                            (
                                ProductionTrajectoryBranch::CompliantContact,
                                ProductionTrajectoryBranch::OpenFlight,
                            ) => ContactTransitionKind::Opening,
                            (
                                ProductionTrajectoryBranch::OpenFlight,
                                ProductionTrajectoryBranch::CompliantContact,
                            ) => ContactTransitionKind::Reimpact,
                            _ => unreachable!("driver records only branch changes"),
                        },
                        time_s,
                        bracket_start_s: transition.bracket_start_s,
                        bracket_end_s: transition.bracket_end_s,
                    }]
                })
                .unwrap_or_default();
            let before_qois = DerivedEulerQois::from_state(
                rigid_step.state_before,
                model.disc_mass_properties,
                0.0,
            )?;
            let after_qois_without_acceleration = DerivedEulerQois::from_state(
                rigid_step.state_after,
                model.disc_mass_properties,
                0.0,
            )?;
            let qois = DerivedEulerQois {
                precession_acceleration_rad_per_s2: (after_qois_without_acceleration
                    .precession_rad_per_s
                    - before_qois.precession_rad_per_s)
                    / duration_s,
                ..after_qois_without_acceleration
            };
            let orientation = next_disc_state.pose().orientation();
            inputs.push(RenderTrajectorySampleInput {
                interval_start_time_s,
                time_s,
                world_frame: RenderWorldFrame::RightHandedZUp,
                units: RenderUnitSystem::SiRadians,
                center_of_mass_world_m: next_disc_state.pose().position_world(),
                orientation_body_to_world: orientation.components(),
                linear_momentum_world_kg_m_per_s: next_disc_state.linear_momentum_world(),
                angular_momentum_body_kg_m2_per_s: next_disc_state.angular_momentum_body(),
                symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
                contact_branch: match endpoint_branch {
                    ProductionTrajectoryBranch::CompliantContact => RenderContactBranch::Closed,
                    ProductionTrajectoryBranch::OpenFlight => RenderContactBranch::Open,
                },
                contact_geometry,
                signed_gap_m,
                interval_contact_active: step.branch
                    == ProductionTrajectoryBranch::CompliantContact,
                interval_normal_force_n: normal_force_n,
                contact_transitions,
                base_mode: Some(RenderBaseModeState {
                    displacement_m: base.modal_displacement_end_m,
                    velocity_m_per_s: base.modal_velocity_end_m_per_s,
                }),
                channels: ChannelOwnership::default(),
                mechanical_energy_j: rigid_step.diagnostics_after.mechanical_energy,
                energy_defect_j,
                qois,
                disposition: if index == final_index {
                    production_event_disposition(&source.termination)
                } else {
                    RenderSampleDisposition::Continue
                },
                terminal_event: None,
            });
            interval_start_time_s = time_s;
            previous_state = next_disc_state;
            previous_base_version = base.next_version;
            last_base_displacement_m = base.modal_displacement_end_m;
            last_base_velocity_m_per_s = base.modal_velocity_end_m_per_s;
        }
        if transition_index != source.transitions.len()
            || previous_state != source.last_accepted_checkpoint.disc_state
            || previous_base_version != source.last_accepted_checkpoint.committed_version
            || interval_start_time_s.to_bits()
                != source.last_accepted_checkpoint.elapsed_time_s().to_bits()
            || last_base_displacement_m.to_bits()
                != source
                    .last_accepted_checkpoint
                    .base_displacement_m()
                    .to_bits()
            || last_base_velocity_m_per_s.to_bits()
                != source
                    .last_accepted_checkpoint
                    .base_velocity_m_per_s()
                    .to_bits()
        {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let identities = profile.content_identities();
        let model_identity = production_model_identity(model, identities.profile);
        let configuration_identity = production_configuration_identity(model, model_identity);
        let mut fingerprint = [0_u8; 8];
        fingerprint.copy_from_slice(&configuration_identity.as_bytes()[..8]);
        let (base_model_id, base_configuration_id) = model.base_port.identity_parts();
        let base_model_identity = hash_domain(
            "org.frankensim.euler-disc.production-base-model.v1",
            format!("{base_model_id}\0{base_configuration_id}").as_bytes(),
        );
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: identities.chart,
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: model.disc_mass_properties,
            },
            initial_state: source.start_checkpoint.disc_state,
            initial_base_mode: RenderBaseModeState {
                displacement_m: source.start_checkpoint.base_displacement_m(),
                velocity_m_per_s: source.start_checkpoint.base_velocity_m_per_s(),
            },
            base_model_identity,
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity,
            channel_availability: RenderChannelAvailability {
                gravity: false,
                contact: false,
                normal_force_sampling: RenderNormalForceSampling::AppliedSubstepZeroOrderHold,
                rolling: false,
                base: false,
                gas: false,
            },
            configuration_identity,
            configuration_fingerprint: u64::from_le_bytes(fingerprint),
            timestep_s: declared_maximum_timestep_s,
            producer_version: format!(
                "fs-euler-disc-e2e/production-event-render-bridge-v{}",
                PRODUCTION_EVENT_RENDER_BRIDGE_VERSION
            ),
            applicability: "accepted fixed-grid open/contact trajectory from the transactional finite-patch/partial-slip/rolling/gas/one-mode-base composition; branch times retain full timestep brackets".to_owned(),
            no_claims: vec![
                "Estimate-authority simulation trajectory; not experimental calibration or validated Euler-disc stopping-time prediction".to_owned(),
                "opening/reimpact times are fixed-grid brackets whose convergence must be demonstrated; no restitution impulse or exact event time is synthesized".to_owned(),
                "normal force is the exact mean of each discretized applied zero-order hold and exactly zero on open intervals; physical bandwidth remains timestep-dependent".to_owned(),
                "render channel work is unavailable pending one shared cross-channel work ledger; zero payloads mean unavailable, not zero physical loss".to_owned(),
                "mechanical energy and defect are disc-only fs-mbd diagnostics, not full disc/base/contact/gas energy closure".to_owned(),
                "no thermal evolution, phase change, melting, plastic flow, structural-mode solve, or acoustic-radiation solve is implied".to_owned(),
            ],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        Self::try_new(metadata, inputs)
    }

    /// Convert a bounded-memory production control trajectory into the common
    /// rendering and structural-acoustics stream.
    ///
    /// Each source interval is homogeneous in its contact branch. Its normal
    /// force is the exact mechanics impulse divided by interval duration, so a
    /// mechanics-to-control reduction preserves every published force-time
    /// cell rather than decimating point samples. Endpoint rigid/base state is
    /// the actual accepted state; no pose or contact event is synthesized.
    pub fn from_production_control_trajectory(
        model: &ProductionCouplingModel,
        source: &ProductionControlTrajectory,
        profile: &ResolvedDiscProfile,
        declared_maximum_timestep_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        cx.checkpoint()
            .map_err(|_| RenderTrajectoryError::Cancelled)?;
        if source.intervals.is_empty() || source.accepted_mechanics_steps == 0 {
            return Err(RenderTrajectoryError::ProductionPrefixEmpty);
        }
        if !(declared_maximum_timestep_s.is_finite() && declared_maximum_timestep_s > 0.0) {
            return Err(RenderTrajectoryError::InvalidTimestep);
        }
        let profile_mass = profile_mass_to_mbd(profile.mass_properties)
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
        model
            .validate_checkpoint(&source.start_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        model
            .validate_checkpoint(&source.last_accepted_checkpoint)
            .map_err(|_| RenderTrajectoryError::ProductionPrefixModelMismatch)?;
        if model.disc_mass_properties != profile_mass {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }
        let accepted_steps = u64::try_from(source.accepted_mechanics_steps).map_err(|_| {
            RenderTrajectoryError::Capacity {
                artifact: "production-control accepted mechanics steps",
                requested: source.accepted_mechanics_steps,
            }
        })?;
        if source
            .start_checkpoint
            .committed_version
            .checked_add(accepted_steps)
            != Some(source.last_accepted_checkpoint.committed_version)
            || source
                .intervals
                .iter()
                .map(|interval| interval.mechanics_substeps)
                .sum::<usize>()
                != source.accepted_mechanics_steps
        {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let mut inputs = Vec::new();
        inputs
            .try_reserve_exact(source.intervals.len())
            .map_err(|_| RenderTrajectoryError::Capacity {
                artifact: "production control render samples",
                requested: source.intervals.len(),
            })?;
        let time_origin_s = source.start_checkpoint.elapsed_time_s();
        let mut previous_time_s = 0.0_f64;
        let mut previous_state = source.start_checkpoint.disc_state;
        let mut previous_base_version = source.start_checkpoint.committed_version;
        let mut transition_index = 0_usize;
        let final_index = source.intervals.len() - 1;
        for (index, interval) in source.intervals.iter().enumerate() {
            if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
                cx.checkpoint()
                    .map_err(|_| RenderTrajectoryError::Cancelled)?;
            }
            let interval_start_time_s = interval.start_time_s - time_origin_s;
            let interval_end_time_s = interval.end_time_s - time_origin_s;
            let duration_s = interval_end_time_s - interval_start_time_s;
            if interval_start_time_s.to_bits() != previous_time_s.to_bits()
                || !(duration_s.is_finite()
                    && duration_s > 0.0
                    && duration_s
                        <= declared_maximum_timestep_s
                            + 64.0 * f64::EPSILON * declared_maximum_timestep_s.max(1.0))
                || interval.mechanics_substeps == 0
                || interval.state_before != previous_state
                || interval.base_parent_version != previous_base_version
                || interval.base_next_version
                    != previous_base_version
                        .checked_add(u64::try_from(interval.mechanics_substeps).map_err(|_| {
                            RenderTrajectoryError::Capacity {
                                artifact: "control interval mechanics substeps",
                                requested: interval.mechanics_substeps,
                            }
                        })?)
                        .ok_or(RenderTrajectoryError::Capacity {
                            artifact: "control interval base version",
                            requested: interval.mechanics_substeps,
                        })?
                || !(interval.mean_normal_force_n.is_finite()
                    && interval.mean_normal_force_n >= 0.0)
                || (interval.normal_impulse_n_s - interval.mean_normal_force_n * duration_s).abs()
                    > 64.0 * f64::EPSILON * interval.normal_impulse_n_s.abs().max(1.0)
                || (interval.branch == ProductionTrajectoryBranch::OpenFlight
                    && (interval.mean_normal_force_n != 0.0 || interval.normal_impulse_n_s != 0.0))
            {
                return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                    sample: index,
                    field: "control interval clock/state/impulse lineage",
                });
            }
            let transition = source
                .transitions
                .get(transition_index)
                .filter(|transition| {
                    transition.bracket_end_s.to_bits() == interval.end_time_s.to_bits()
                });
            let endpoint_branch = transition.map_or(interval.branch, |transition| transition.to);
            if let Some(transition) = transition {
                if transition.from != interval.branch
                    || transition.bracket_start_s < interval.start_time_s
                    || transition.bracket_start_s > transition.bracket_end_s
                {
                    return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                        sample: index,
                        field: "control event transition lineage",
                    });
                }
                transition_index += 1;
            }
            let endpoint_contact = profile_contact_geometry(
                &profile.chart,
                profile.mass_properties,
                interval.state_after.pose(),
                cx,
            )
            .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
            let signed_gap_m = endpoint_contact.contact.gap_m - interval.base_displacement_end_m;
            let contact_geometry = match endpoint_branch {
                ProductionTrajectoryBranch::CompliantContact => {
                    if signed_gap_m > 0.0 {
                        return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                            sample: index,
                            field: "closed control endpoint has positive bulk-surface gap",
                        });
                    }
                    Some(RenderContactGeometry {
                        point_world_m: endpoint_contact.contact.point_world_m,
                        normal_world: Vec3::new(0.0, 0.0, 1.0),
                        support_feature: RenderSupportFeature::ProfileFeature(
                            endpoint_contact.support_source_feature,
                        ),
                    })
                }
                ProductionTrajectoryBranch::OpenFlight => {
                    if signed_gap_m < 0.0 {
                        return Err(RenderTrajectoryError::ProductionPrefixSampleMismatch {
                            sample: index,
                            field: "open control endpoint has negative bulk-surface gap",
                        });
                    }
                    None
                }
            };
            let contact_transitions = transition
                .map(|transition| {
                    vec![RenderContactTransition {
                        kind: match (transition.from, transition.to) {
                            (
                                ProductionTrajectoryBranch::CompliantContact,
                                ProductionTrajectoryBranch::OpenFlight,
                            ) => ContactTransitionKind::Opening,
                            (
                                ProductionTrajectoryBranch::OpenFlight,
                                ProductionTrajectoryBranch::CompliantContact,
                            ) => ContactTransitionKind::Reimpact,
                            _ => unreachable!("only branch changes are retained"),
                        },
                        time_s: interval_end_time_s,
                        bracket_start_s: transition.bracket_start_s - time_origin_s,
                        bracket_end_s: transition.bracket_end_s - time_origin_s,
                    }]
                })
                .unwrap_or_default();
            let before_qois = DerivedEulerQois::from_state(
                interval.state_before,
                model.disc_mass_properties,
                0.0,
            )?;
            let after_qois_without_acceleration = DerivedEulerQois::from_state(
                interval.state_after,
                model.disc_mass_properties,
                0.0,
            )?;
            let qois = DerivedEulerQois {
                precession_acceleration_rad_per_s2: (after_qois_without_acceleration
                    .precession_rad_per_s
                    - before_qois.precession_rad_per_s)
                    / duration_s,
                ..after_qois_without_acceleration
            };
            let orientation = interval.state_after.pose().orientation();
            inputs.push(RenderTrajectorySampleInput {
                interval_start_time_s,
                time_s: interval_end_time_s,
                world_frame: RenderWorldFrame::RightHandedZUp,
                units: RenderUnitSystem::SiRadians,
                center_of_mass_world_m: interval.state_after.pose().position_world(),
                orientation_body_to_world: orientation.components(),
                linear_momentum_world_kg_m_per_s: interval.state_after.linear_momentum_world(),
                angular_momentum_body_kg_m2_per_s: interval.state_after.angular_momentum_body(),
                symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
                contact_branch: match endpoint_branch {
                    ProductionTrajectoryBranch::CompliantContact => RenderContactBranch::Closed,
                    ProductionTrajectoryBranch::OpenFlight => RenderContactBranch::Open,
                },
                contact_geometry,
                signed_gap_m,
                interval_contact_active: interval.branch
                    == ProductionTrajectoryBranch::CompliantContact,
                interval_normal_force_n: interval.mean_normal_force_n,
                contact_transitions,
                base_mode: Some(RenderBaseModeState {
                    displacement_m: interval.base_displacement_end_m,
                    velocity_m_per_s: interval.base_velocity_end_m_per_s,
                }),
                channels: ChannelOwnership::default(),
                mechanical_energy_j: interval.mechanical_energy_end_j,
                energy_defect_j: interval.disc_work_residual_j,
                qois,
                disposition: if index == final_index {
                    production_event_disposition(&source.termination)
                } else {
                    RenderSampleDisposition::Continue
                },
                terminal_event: None,
            });
            previous_time_s = interval_end_time_s;
            previous_state = interval.state_after;
            previous_base_version = interval.base_next_version;
        }
        if transition_index != source.transitions.len()
            || previous_time_s.to_bits()
                != (source.last_accepted_checkpoint.elapsed_time_s() - time_origin_s).to_bits()
            || previous_state != source.last_accepted_checkpoint.disc_state
            || previous_base_version != source.last_accepted_checkpoint.committed_version
        {
            return Err(RenderTrajectoryError::ProductionPrefixModelMismatch);
        }

        let identities = profile.content_identities();
        let model_identity = production_model_identity(model, identities.profile);
        let configuration_identity = production_configuration_identity(model, model_identity);
        let mut fingerprint = [0_u8; 8];
        fingerprint.copy_from_slice(&configuration_identity.as_bytes()[..8]);
        let (base_model_id, base_configuration_id) = model.base_port.identity_parts();
        let base_model_identity = hash_domain(
            "org.frankensim.euler-disc.production-base-model.v1",
            format!("{base_model_id}\0{base_configuration_id}").as_bytes(),
        );
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: identities.chart,
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: model.disc_mass_properties,
            },
            initial_state: source.start_checkpoint.disc_state,
            initial_base_mode: RenderBaseModeState {
                displacement_m: source.start_checkpoint.base_displacement_m(),
                velocity_m_per_s: source.start_checkpoint.base_velocity_m_per_s(),
            },
            base_model_identity,
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity,
            channel_availability: RenderChannelAvailability {
                gravity: false,
                contact: false,
                normal_force_sampling: RenderNormalForceSampling::AppliedSubstepZeroOrderHold,
                rolling: false,
                base: false,
                gas: false,
            },
            configuration_identity,
            configuration_fingerprint: u64::from_le_bytes(fingerprint),
            timestep_s: declared_maximum_timestep_s,
            producer_version: format!(
                "fs-euler-disc-e2e/production-control-render-bridge-v{}",
                PRODUCTION_EVENT_RENDER_BRIDGE_VERSION
            ),
            applicability: "bounded-memory homogeneous controls reduced from the transactional finite-patch/surface/partial-slip/rolling/gas/modal-base production composition; normal impulse is exact across reduction cells".to_owned(),
            no_claims: vec![
                "Estimate-authority simulation trajectory; not experimental calibration or validated Euler-disc stopping-time prediction".to_owned(),
                "mechanics-to-control reduction preserves accepted impulse and endpoints but does not establish mechanics timestep convergence".to_owned(),
                "opening/reimpact times retain mechanics-grid brackets; no restitution impulse or exact event time is synthesized".to_owned(),
                "mechanical energy and defect are disc-only fs-mbd diagnostics, not full disc/base/contact/gas energy closure".to_owned(),
                "no thermal evolution, phase change, melting, plastic flow, or phase-dependent remeshing is implied".to_owned(),
            ],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        Self::try_new(metadata, inputs)
    }

    /// Build a grounded animation trajectory from an admitted late-stage
    /// reduced-decay run and the exact resolved filleted specimen.
    ///
    /// The bridge integrates the declared precession and residual-spin phases
    /// with `phi_dot = Omega` and `psi_dot = -Omega cos(theta)`. The inclination
    /// rate comes from the retained dissipation divided by `dE/dtheta`; all
    /// three rates are included in the rigid twist before principal-body
    /// angular momentum is formed. Every pose is placed on the actual profile
    /// support through `profile_state_at_ground_contact`. The final
    /// [`REDUCED_DECAY_RENDER_TAIL_HORIZON_S`] seconds are rebased to `t = 0`;
    /// an exact crop-boundary sample follows the source integrator's
    /// left-boundary power law, so no physical coefficient or clock is changed.
    ///
    /// Rolling and Bildsten losses are retained only as interval work. No
    /// force/torque is invented for either energy-only closure. Each
    /// positive-duration interval carries the mean base-normal reaction
    /// required by exact endpoint vertical impulse balance. This is a
    /// kinematic reconstruction, not a resolved contact force history.
    /// Reaching the positive validity cutoff is published as horizon
    /// censoring, never as a `theta = 0`, contact-loss, or localized physical
    /// terminal event.
    pub fn from_reduced_decay_run(
        run: &ReducedDecayRun,
        profile: &ResolvedDiscProfile,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        Self::from_reduced_decay_run_tail(
            run,
            profile,
            REDUCED_DECAY_RENDER_TAIL_HORIZON_S,
            0.0,
            REDUCED_DECAY_RENDER_TAIL_HORIZON_S,
            None,
            "final eight-second tail of the Thorne et al. 2026 source-bound small-angle analytical steel-on-glass decay, rebased without time scaling and ending at a positive validity cutoff".to_owned(),
            cx,
        )
    }

    /// Build a longer source-bound tail for fixture-internal acoustic preroll.
    ///
    /// This keeps the public eight-second visual bridge unchanged. The caller
    /// must identify the longer tail as an internal warm-start source rather
    /// than presenting it as the picture trajectory.
    pub(crate) fn from_reduced_decay_run_with_tail_horizon(
        run: &ReducedDecayRun,
        profile: &ResolvedDiscProfile,
        tail_horizon_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        Self::from_reduced_decay_run_tail(
            run,
            profile,
            tail_horizon_s,
            0.0,
            tail_horizon_s,
            None,
            format!(
                "final {tail_horizon_s:.17e}-second internal acoustic-preroll tail of the Thorne et al. 2026 source-bound small-angle analytical steel-on-glass decay, rebased without time scaling and ending at a positive validity cutoff"
            ),
            cx,
        )
    }

    /// Build a causal warm-start source and the following picture/sound crop
    /// from one phase-continuous reduced-decay bridge.
    ///
    /// Both trajectories integrate the same exact inserted crop boundary. The
    /// source begins one preroll before the published trajectory; the published
    /// trajectory rebases that boundary to zero without resetting orientation,
    /// contact phase, translation, or any physical state.
    pub(crate) fn from_reduced_decay_run_with_causal_preroll(
        run: &ReducedDecayRun,
        profile: &ResolvedDiscProfile,
        preroll_s: f64,
        published_horizon_s: f64,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self), RenderTrajectoryError> {
        if !(preroll_s.is_finite()
            && preroll_s > 0.0
            && published_horizon_s.is_finite()
            && published_horizon_s > 0.0)
        {
            return Err(reduced_decay_bridge_refusal(
                "causal_preroll_window",
                "preroll and published horizons must be finite and positive",
            ));
        }
        let source_horizon_s = preroll_s + published_horizon_s;
        if !source_horizon_s.is_finite() {
            return Err(reduced_decay_bridge_refusal(
                "causal_preroll_window",
                "source horizon overflow",
            ));
        }
        let source = Self::from_reduced_decay_run_tail(
            run,
            profile,
            source_horizon_s,
            0.0,
            source_horizon_s,
            Some(preroll_s),
            format!(
                "final {source_horizon_s:.17e}-second causal picture-and-sound source tail with an exact crop boundary at {preroll_s:.17e} seconds"
            ),
            cx,
        )?;
        let published = Self::from_reduced_decay_run_tail(
            run,
            profile,
            source_horizon_s,
            preroll_s,
            published_horizon_s,
            Some(preroll_s),
            format!(
                "{published_horizon_s:.17e}-second published picture-and-sound crop following a {preroll_s:.17e}-second causal modal warm start; time is rebased without resetting pose, contact phase, translation, or mechanical state"
            ),
            cx,
        )?;
        Ok((source, published))
    }

    fn from_reduced_decay_run_tail(
        run: &ReducedDecayRun,
        profile: &ResolvedDiscProfile,
        tail_horizon_s: f64,
        publish_start_offset_s: f64,
        published_horizon_s: f64,
        required_boundary_offset_s: Option<f64>,
        applicability: String,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryError> {
        cx.checkpoint()
            .map_err(|_| RenderTrajectoryError::Cancelled)?;
        let mass = validate_reduced_decay_render_source(run, profile, cx)?;
        let inputs = reduced_decay_sample_inputs(
            run,
            profile,
            mass,
            tail_horizon_s,
            publish_start_offset_s,
            published_horizon_s,
            required_boundary_offset_s,
            cx,
        )?;
        let initial = inputs
            .first()
            .ok_or(RenderTrajectoryError::EmptyTrajectory)?;
        let initial_orientation =
            checked_unit_quaternion(initial.orientation_body_to_world, 0, false)?;
        let initial_state = RigidBodyState::new(
            Pose::new(initial.center_of_mass_world_m, initial_orientation).map_err(|error| {
                RenderTrajectoryError::ReducedDecayBridgeRefusal {
                    field: "initial_pose",
                    detail: error.to_string(),
                }
            })?,
            initial.linear_momentum_world_kg_m_per_s,
            initial.angular_momentum_body_kg_m2_per_s,
        )
        .map_err(|error| RenderTrajectoryError::ReducedDecayBridgeRefusal {
            field: "initial_state",
            detail: error.to_string(),
        })?;
        let identities = profile.content_identities();
        let model_identity = reduced_decay_model_identity(run);
        let source_configuration_identity = reduced_decay_configuration_identity(
            run,
            model_identity,
            identities.profile,
            identities.chart,
            identities.mass_properties,
        );
        let mut configuration_hasher =
            DomainHasher::new("org.frankensim.euler-disc.reduced-decay-render-window.v3");
        configuration_hasher.update(source_configuration_identity.as_bytes());
        configuration_hasher.update(&tail_horizon_s.to_bits().to_le_bytes());
        configuration_hasher.update(&publish_start_offset_s.to_bits().to_le_bytes());
        configuration_hasher.update(&published_horizon_s.to_bits().to_le_bytes());
        match required_boundary_offset_s {
            Some(offset_s) => {
                configuration_hasher.update(&[1]);
                configuration_hasher.update(&offset_s.to_bits().to_le_bytes());
            }
            None => configuration_hasher.update(&[0]),
        }
        let configuration_identity = configuration_hasher.finalize();
        let mut fingerprint_bytes = [0_u8; 8];
        fingerprint_bytes.copy_from_slice(&configuration_identity.as_bytes()[..8]);
        let metadata = RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identities.profile,
            specimen_chart_identity: identities.chart,
            mass_properties: RenderMassProperties {
                identity: identities.mass_properties,
                properties: mass,
            },
            initial_state,
            initial_base_mode: RenderBaseModeState {
                displacement_m: 0.0,
                velocity_m_per_s: 0.0,
            },
            base_model_identity: hash_domain(
                REDUCED_DECAY_BASE_IDENTITY_DOMAIN,
                b"static-rigid-z-up-ground-plane-no-base-dynamics",
            ),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity,
            channel_availability: RenderChannelAvailability {
                gravity: false,
                contact: false,
                normal_force_sampling: RenderNormalForceSampling::IntervalMean,
                rolling: true,
                base: false,
                gas: true,
            },
            configuration_identity,
            configuration_fingerprint: u64::from_le_bytes(fingerprint_bytes),
            timestep_s: run.parameters.timestep_s,
            producer_version: format!(
                "fs-euler-disc-e2e/reduced-decay-render-bridge-v{}",
                REDUCED_DECAY_RENDER_BRIDGE_VERSION
            ),
            applicability,
            no_claims: vec![
                "literature-calibrated analytical trajectory; not a raw measured trajectory"
                    .to_owned(),
                "Bildsten gas loss is energy-only; no aerodynamic wrench or exact full-FSI prefactor validation"
                    .to_owned(),
                "validity cutoff is horizon censoring; no theta-zero, loss-of-contact, or terminal-event claim"
                    .to_owned(),
                "ground contact is kinematic profile support; the separately available normal load is only the kinematically implied interval mean required by endpoint vertical impulse balance under gravity-plus-support closure, not a resolved subinterval force history or measured force"
                    .to_owned(),
                "the full contact wrench/work channel remains unavailable; no contact-torque, angular-impulse, tangential-traction, pressure-patch, deformation, or acoustic-radiation claim is made"
                    .to_owned(),
                "not a configuration, design, specimen, or target-ranking claim".to_owned(),
            ],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        };
        Self::try_new(metadata, inputs)
    }

    /// Validate metadata, all samples, cross-sample time/event rules, and the
    /// final terminal/censor placement.
    pub fn try_new(
        metadata: RenderTrajectoryMetadata,
        inputs: Vec<RenderTrajectorySampleInput>,
    ) -> Result<Self, RenderTrajectoryError> {
        Self::try_new_with_policy(metadata, inputs, false, &mut || Ok(()))
    }

    /// Re-admit canonical wire components without renormalizing quaternions,
    /// polling the supplied cancellation boundary at a fixed sample cadence.
    pub(crate) fn try_new_canonical(
        metadata: RenderTrajectoryMetadata,
        inputs: Vec<RenderTrajectorySampleInput>,
        checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryError>,
    ) -> Result<Self, RenderTrajectoryError> {
        Self::try_new_with_policy(metadata, inputs, true, checkpoint)
    }

    fn try_new_with_policy(
        metadata: RenderTrajectoryMetadata,
        inputs: Vec<RenderTrajectorySampleInput>,
        canonical_quaternions: bool,
        checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryError>,
    ) -> Result<Self, RenderTrajectoryError> {
        checkpoint()?;
        validate_metadata(&metadata)?;
        if inputs.is_empty() {
            return Err(RenderTrajectoryError::EmptyTrajectory);
        }
        if inputs.len() > MAX_RENDER_TRAJECTORY_SAMPLES {
            return Err(RenderTrajectoryError::TooManySamples(inputs.len()));
        }

        let sample_count = inputs.len();
        let mut samples = Vec::new();
        samples
            .try_reserve_exact(sample_count)
            .map_err(|_| RenderTrajectoryError::Capacity {
                artifact: "render trajectory samples",
                requested: sample_count,
            })?;
        let mut previous_time = None;
        let mut previous_branch = None;
        for (index, input) in inputs.into_iter().enumerate() {
            if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
                checkpoint()?;
            }
            let sample = validate_sample(
                &metadata,
                input,
                index,
                previous_time,
                canonical_quaternions,
            )?;
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
        checkpoint()?;
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
    /// A cancellation-aware admission caller requested an early stop.
    Cancelled,
    /// Only the exact v3 schema is currently accepted.
    UnsupportedSchemaVersion(u16),
    /// Only right-handed `+z`-up coordinates are admitted by v3.
    UnsupportedWorldFrame(RenderWorldFrame),
    /// Only SI/radians are admitted by v3.
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
    /// A production smooth-contact prefix accepted no mechanics substep.
    ProductionPrefixEmpty,
    /// Production model, specimen, start, or final checkpoint identities disagreed.
    ProductionPrefixModelMismatch,
    /// One accepted production receipt did not extend the preceding state exactly.
    ProductionPrefixSampleMismatch {
        /// Zero-based accepted-step index.
        sample: usize,
        /// Stable mismatch category.
        field: &'static str,
    },
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
    /// Full contact-wrench authority was declared without its required normal
    /// load diagnostic.
    InvalidChannelAvailability(&'static str),
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
    /// A source-bound reduced-decay run, profile, or derived bridge quantity
    /// contradicted the bridge contract.
    ReducedDecayBridgeRefusal {
        /// Stable semantic field or invariant.
        field: &'static str,
        /// Upstream or numerical detail.
        detail: String,
    },
}

impl fmt::Display for RenderTrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RenderTrajectoryError {}

fn reduced_decay_bridge_refusal(
    field: &'static str,
    detail: impl Into<String>,
) -> RenderTrajectoryError {
    RenderTrajectoryError::ReducedDecayBridgeRefusal {
        field,
        detail: detail.into(),
    }
}

fn validate_reduced_decay_render_source(
    run: &ReducedDecayRun,
    profile: &ResolvedDiscProfile,
    cx: &Cx<'_>,
) -> Result<MassProperties, RenderTrajectoryError> {
    let provenance = &run.provenance;
    if provenance.model_id != REDUCED_DECAY_MODEL_ID
        || provenance.small_angle_oracle_source_id != THORNE_2026_SOURCE_ID
        || provenance.model_authority != "literature-calibrated-analytical"
        || provenance.physical_validation != "no-raw-trajectory-or-full-fsi-validation-claimed"
        || provenance.published_rolling_source_id.as_deref() != Some(THORNE_2026_SOURCE_ID)
        || provenance.bildsten_source_id.as_deref() != Some(THORNE_2026_SOURCE_ID)
        || provenance
            .published_rolling_coefficient_mu
            .map(f64::to_bits)
            != Some(THORNE_2026_FITTED_ROLLING_COEFFICIENT.to_bits())
        || provenance.bildsten_dynamic_viscosity_pa_s.map(f64::to_bits)
            != Some(THORNE_2026_DECLARED_AIR_VISCOSITY_PA_S.to_bits())
        || provenance
            .bildsten_dimensionless_prefactor
            .map(f64::to_bits)
            != Some(1.0_f64.to_bits())
    {
        return Err(reduced_decay_bridge_refusal(
            "provenance",
            "run is not the exact source-bound Thorne 2026 analytical model",
        ));
    }
    let density_kg_per_m3 = provenance.bildsten_density_kg_per_m3.ok_or_else(|| {
        reduced_decay_bridge_refusal("provenance.bildsten_density_kg_per_m3", "missing")
    })?;
    if !(density_kg_per_m3.is_finite() && density_kg_per_m3 > 0.0) {
        return Err(reduced_decay_bridge_refusal(
            "provenance.bildsten_density_kg_per_m3",
            "must be finite and positive",
        ));
    }
    let specimen = provenance
        .literature_specimen
        .as_ref()
        .ok_or_else(|| reduced_decay_bridge_refusal("provenance.literature_specimen", "missing"))?;
    if specimen.source_id != THORNE_2026_SOURCE_ID || profile.spec != specimen.profile_spec() {
        return Err(reduced_decay_bridge_refusal(
            "resolved_profile.spec",
            "resolved profile does not match the source-bound filleted specimen",
        ));
    }
    if run.parameters.mass_kg.to_bits() != specimen.mass_kg.to_bits()
        || run.parameters.radius_m.to_bits() != (0.5 * specimen.diameter_m).to_bits()
        || profile.dimensions.outer_radius_m.to_bits() != run.parameters.radius_m.to_bits()
        || profile.dimensions.thickness_m.to_bits() != specimen.thickness_m.to_bits()
        || !scalar_close(profile.mass_properties.mass, specimen.mass_kg)
    {
        return Err(reduced_decay_bridge_refusal(
            "resolved_profile.mass_or_dimensions",
            "profile geometry/mass and reduced-decay parameters disagree",
        ));
    }
    for (field, value) in [
        ("parameters.mass_kg", run.parameters.mass_kg),
        ("parameters.radius_m", run.parameters.radius_m),
        (
            "parameters.gravity_m_per_s2",
            run.parameters.gravity_m_per_s2,
        ),
        (
            "parameters.initial_theta_rad",
            run.parameters.initial_theta_rad,
        ),
        (
            "parameters.validity_cutoff_theta_rad",
            run.parameters.validity_cutoff_theta_rad,
        ),
        ("parameters.timestep_s", run.parameters.timestep_s),
    ] {
        if !(value.is_finite() && value > 0.0) {
            return Err(reduced_decay_bridge_refusal(
                field,
                "must be finite and positive",
            ));
        }
    }
    if run.parameters.initial_theta_rad <= run.parameters.validity_cutoff_theta_rad
        || run.parameters.maximum_steps == 0
        || run.terminal != ReducedDecayTerminal::ValidityCutoff
        || run.samples.is_empty()
        || run.samples.len()
            > usize::try_from(run.parameters.maximum_steps)
                .unwrap_or(usize::MAX)
                .saturating_add(1)
    {
        return Err(reduced_decay_bridge_refusal(
            "run.bounds_or_terminal",
            "run must be a nonempty bounded validity-cutoff trajectory",
        ));
    }

    let energy_slope_j_per_rad =
        1.5 * run.parameters.mass_kg * run.parameters.gravity_m_per_s2 * run.parameters.radius_m;
    let rolling_coefficient = provenance.published_rolling_coefficient_mu.ok_or_else(|| {
        reduced_decay_bridge_refusal("provenance.published_rolling_coefficient_mu", "missing")
    })?;
    let viscosity_pa_s = provenance.bildsten_dynamic_viscosity_pa_s.ok_or_else(|| {
        reduced_decay_bridge_refusal("provenance.bildsten_dynamic_viscosity_pa_s", "missing")
    })?;
    let bildsten_prefactor = provenance.bildsten_dimensionless_prefactor.ok_or_else(|| {
        reduced_decay_bridge_refusal("provenance.bildsten_dimensionless_prefactor", "missing")
    })?;
    let mut previous: Option<&ReducedDecaySample> = None;
    for (index, sample) in run.samples.iter().enumerate() {
        if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
            cx.checkpoint()
                .map_err(|_| RenderTrajectoryError::Cancelled)?;
        }
        for (field, value) in [
            ("sample.time_s", sample.time_s),
            ("sample.theta_rad", sample.theta_rad),
            ("sample.omega_rad_s", sample.omega_rad_s),
            ("sample.energy_j", sample.energy_j),
            ("sample.powers.dry_contour_w", sample.powers.dry_contour_w),
            (
                "sample.powers.published_rolling_w",
                sample.powers.published_rolling_w,
            ),
            (
                "sample.powers.bildsten_boundary_layer_w",
                sample.powers.bildsten_boundary_layer_w,
            ),
            ("sample.work.dry_contour_j", sample.work.dry_contour_j),
            (
                "sample.work.published_rolling_j",
                sample.work.published_rolling_j,
            ),
            (
                "sample.work.bildsten_boundary_layer_j",
                sample.work.bildsten_boundary_layer_j,
            ),
        ] {
            if !value.is_finite() {
                return Err(reduced_decay_bridge_refusal(
                    field,
                    format!("sample {index} is non-finite"),
                ));
            }
        }
        if sample.time_s < 0.0
            || sample.theta_rad <= 0.0
            || sample.powers.dry_contour_w != 0.0
            || sample.work.dry_contour_j != 0.0
            || sample.powers.published_rolling_w <= 0.0
            || sample.powers.bildsten_boundary_layer_w <= 0.0
            || sample.work.published_rolling_j < 0.0
            || sample.work.bildsten_boundary_layer_j < 0.0
        {
            return Err(reduced_decay_bridge_refusal(
                "sample.channel_domain",
                format!("sample {index} contradicts the benchmark channel domains"),
            ));
        }
        let expected_omega = (4.0 * run.parameters.gravity_m_per_s2
            / (run.parameters.radius_m * sample.theta_rad.sin()))
        .sqrt();
        let expected_energy = energy_slope_j_per_rad * sample.theta_rad;
        let expected_rolling = rolling_coefficient
            * run.parameters.mass_kg
            * run.parameters.gravity_m_per_s2
            * run.parameters.radius_m
            * sample.theta_rad.cos()
            * sample.omega_rad_s;
        let expected_air = bildsten_prefactor
            * BILDSTEN_PUBLISHED_POWER_COEFFICIENT
            * (viscosity_pa_s * density_kg_per_m3).sqrt()
            * run.parameters.gravity_m_per_s2.powf(1.25)
            * run.parameters.radius_m.powf(2.75)
            * sample.theta_rad.powf(-1.25);
        if !scalar_close(sample.omega_rad_s, expected_omega)
            || !scalar_close(sample.energy_j, expected_energy)
            || !scalar_close(sample.powers.published_rolling_w, expected_rolling)
            || !scalar_close(sample.powers.bildsten_boundary_layer_w, expected_air)
        {
            return Err(reduced_decay_bridge_refusal(
                "sample.reduced_equations",
                format!("sample {index} does not satisfy the retained analytical equations"),
            ));
        }
        if let Some(previous) = previous {
            let dt_s = sample.time_s - previous.time_s;
            if !(dt_s.is_finite() && dt_s > 0.0) || sample.theta_rad >= previous.theta_rad {
                return Err(reduced_decay_bridge_refusal(
                    "sample.time_or_theta_order",
                    format!("sample {index} is not a bounded descending step"),
                ));
            }
            let time_to_cutoff_s = (previous.theta_rad - run.parameters.validity_cutoff_theta_rad)
                * energy_slope_j_per_rad
                / previous.powers.total_w();
            let expected_dt_s = run.parameters.timestep_s.min(time_to_cutoff_s);
            let expected_theta_rad = if time_to_cutoff_s <= run.parameters.timestep_s {
                run.parameters.validity_cutoff_theta_rad
            } else {
                previous.theta_rad
                    - previous.powers.total_w() * expected_dt_s / energy_slope_j_per_rad
            };
            // `time_s` is accumulated by repeated floating-point addition in
            // the producer. Subtracting adjacent late samples therefore has an
            // error measured in ULPs of the *absolute clock*, not ULPs of the
            // small fixed timestep. Retain a bounded endpoint-scale allowance
            // while still checking the exact left-boundary step model below.
            let time_difference_tolerance_s = 32.0
                * f64::EPSILON
                * previous
                    .time_s
                    .abs()
                    .max(sample.time_s.abs())
                    .max(expected_dt_s.abs());
            if (dt_s - expected_dt_s).abs() > time_difference_tolerance_s
                || !scalar_close(sample.theta_rad, expected_theta_rad)
            {
                return Err(reduced_decay_bridge_refusal(
                    "sample.integration_step",
                    format!("sample {index} does not follow the admitted left-boundary decay step"),
                ));
            }
            let rolling_increment =
                sample.work.published_rolling_j - previous.work.published_rolling_j;
            let gas_increment =
                sample.work.bildsten_boundary_layer_j - previous.work.bildsten_boundary_layer_j;
            if !scalar_close(
                rolling_increment,
                previous.powers.published_rolling_w * dt_s,
            ) || !scalar_close(
                gas_increment,
                previous.powers.bildsten_boundary_layer_w * dt_s,
            ) {
                return Err(reduced_decay_bridge_refusal(
                    "sample.interval_work",
                    format!("sample {index} work does not close against left-boundary power"),
                ));
            }
        } else if sample.time_s.to_bits() != 0.0_f64.to_bits()
            || sample.theta_rad.to_bits() != run.parameters.initial_theta_rad.to_bits()
            || sample.work.total_j() != 0.0
        {
            return Err(reduced_decay_bridge_refusal(
                "sample.initial",
                "initial sample must be the exact zero-time, zero-work initial state",
            ));
        }
        previous = Some(sample);
    }
    let final_sample = run
        .samples
        .last()
        .ok_or(RenderTrajectoryError::EmptyTrajectory)?;
    if final_sample.theta_rad.to_bits() != run.parameters.validity_cutoff_theta_rad.to_bits() {
        return Err(reduced_decay_bridge_refusal(
            "sample.final_theta_rad",
            "final sample does not equal the positive validity cutoff",
        ));
    }
    let closure = energy_slope_j_per_rad * run.parameters.initial_theta_rad
        - final_sample.energy_j
        - final_sample.work.total_j();
    if !scalar_close(closure, run.energy_closure_residual_j) {
        return Err(reduced_decay_bridge_refusal(
            "energy_closure_residual_j",
            "run-level energy residual does not match the retained samples",
        ));
    }
    profile_mass_to_mbd(profile.mass_properties).map_err(|error| {
        reduced_decay_bridge_refusal("resolved_profile.mass_properties", error.to_string())
    })
}

fn reduced_decay_tail_samples(
    run: &ReducedDecayRun,
    tail_horizon_s: f64,
    required_boundary_offset_s: Option<f64>,
) -> Result<(f64, Vec<ReducedDecaySample>), RenderTrajectoryError> {
    if !(tail_horizon_s.is_finite() && tail_horizon_s > 0.0) {
        return Err(reduced_decay_bridge_refusal(
            "tail_horizon_s",
            "must be finite and positive",
        ));
    }
    let final_sample = run
        .samples
        .last()
        .ok_or(RenderTrajectoryError::EmptyTrajectory)?;
    let time_origin_s = (final_sample.time_s - tail_horizon_s).max(0.0);
    let first_retained = run
        .samples
        .partition_point(|sample| sample.time_s < time_origin_s);
    if first_retained >= run.samples.len() {
        return Err(reduced_decay_bridge_refusal(
            "tail_horizon_s",
            "crop origin lies beyond the final retained sample",
        ));
    }
    let needs_boundary = first_retained > 0
        && run.samples[first_retained].time_s.to_bits() != time_origin_s.to_bits();
    let requested = run
        .samples
        .len()
        .saturating_sub(first_retained)
        .saturating_add(usize::from(needs_boundary));
    let mut samples = Vec::new();
    samples
        .try_reserve_exact(requested)
        .map_err(|_| RenderTrajectoryError::Capacity {
            artifact: "reduced-decay cropped source samples",
            requested,
        })?;
    if needs_boundary {
        let previous = &run.samples[first_retained - 1];
        let dt_s = time_origin_s - previous.time_s;
        let energy_slope_j_per_rad = 1.5
            * run.parameters.mass_kg
            * run.parameters.gravity_m_per_s2
            * run.parameters.radius_m;
        let theta_rad =
            previous.theta_rad - previous.powers.total_w() * dt_s / energy_slope_j_per_rad;
        let (omega_rad_s, powers) = reduced_decay_benchmark_powers(run, theta_rad)?;
        let work = ChannelWork {
            dry_contour_j: previous.work.dry_contour_j + previous.powers.dry_contour_w * dt_s,
            published_rolling_j: previous.work.published_rolling_j
                + previous.powers.published_rolling_w * dt_s,
            bildsten_boundary_layer_j: previous.work.bildsten_boundary_layer_j
                + previous.powers.bildsten_boundary_layer_w * dt_s,
        };
        samples.push(ReducedDecaySample {
            time_s: time_origin_s,
            theta_rad,
            omega_rad_s,
            energy_j: energy_slope_j_per_rad * theta_rad,
            powers,
            work,
        });
    }
    samples.extend_from_slice(&run.samples[first_retained..]);
    if let Some(offset_s) = required_boundary_offset_s {
        if !(offset_s.is_finite() && offset_s > 0.0 && offset_s < tail_horizon_s) {
            return Err(reduced_decay_bridge_refusal(
                "required_boundary_offset_s",
                "must be finite and strictly inside the retained tail",
            ));
        }
        let boundary_time_s = time_origin_s + offset_s;
        match samples.binary_search_by(|sample| sample.time_s.total_cmp(&boundary_time_s)) {
            Ok(_) => {}
            Err(index) => {
                let boundary = interpolated_reduced_decay_sample(run, boundary_time_s)?;
                samples
                    .try_reserve(1)
                    .map_err(|_| RenderTrajectoryError::Capacity {
                        artifact: "reduced-decay exact crop boundary",
                        requested: samples.len().saturating_add(1),
                    })?;
                samples.insert(index, boundary);
            }
        }
    }
    if samples.is_empty() {
        return Err(RenderTrajectoryError::EmptyTrajectory);
    }
    Ok((time_origin_s, samples))
}

fn interpolated_reduced_decay_sample(
    run: &ReducedDecayRun,
    time_s: f64,
) -> Result<ReducedDecaySample, RenderTrajectoryError> {
    if !time_s.is_finite() {
        return Err(reduced_decay_bridge_refusal(
            "crop.time_s",
            "must be finite",
        ));
    }
    let right = run.samples.partition_point(|sample| sample.time_s < time_s);
    if right < run.samples.len() && run.samples[right].time_s.to_bits() == time_s.to_bits() {
        return Ok(run.samples[right].clone());
    }
    if right == 0 || right >= run.samples.len() {
        return Err(reduced_decay_bridge_refusal(
            "crop.time_s",
            "lies outside the admitted reduced-decay run",
        ));
    }
    let previous = &run.samples[right - 1];
    let dt_s = time_s - previous.time_s;
    let energy_slope_j_per_rad =
        1.5 * run.parameters.mass_kg * run.parameters.gravity_m_per_s2 * run.parameters.radius_m;
    let theta_rad = previous.theta_rad - previous.powers.total_w() * dt_s / energy_slope_j_per_rad;
    let (omega_rad_s, powers) = reduced_decay_benchmark_powers(run, theta_rad)?;
    let work = ChannelWork {
        dry_contour_j: previous.work.dry_contour_j + previous.powers.dry_contour_w * dt_s,
        published_rolling_j: previous.work.published_rolling_j
            + previous.powers.published_rolling_w * dt_s,
        bildsten_boundary_layer_j: previous.work.bildsten_boundary_layer_j
            + previous.powers.bildsten_boundary_layer_w * dt_s,
    };
    Ok(ReducedDecaySample {
        time_s,
        theta_rad,
        omega_rad_s,
        energy_j: energy_slope_j_per_rad * theta_rad,
        powers,
        work,
    })
}

fn reduced_decay_benchmark_powers(
    run: &ReducedDecayRun,
    theta_rad: f64,
) -> Result<(f64, ChannelPowers), RenderTrajectoryError> {
    if !(theta_rad.is_finite() && theta_rad > 0.0) {
        return Err(reduced_decay_bridge_refusal(
            "crop.theta_rad",
            "interpolated inclination must be finite and positive",
        ));
    }
    let rolling_coefficient = run
        .provenance
        .published_rolling_coefficient_mu
        .ok_or_else(|| {
            reduced_decay_bridge_refusal("provenance.published_rolling_coefficient_mu", "missing")
        })?;
    let density_kg_per_m3 = run.provenance.bildsten_density_kg_per_m3.ok_or_else(|| {
        reduced_decay_bridge_refusal("provenance.bildsten_density_kg_per_m3", "missing")
    })?;
    let viscosity_pa_s = run
        .provenance
        .bildsten_dynamic_viscosity_pa_s
        .ok_or_else(|| {
            reduced_decay_bridge_refusal("provenance.bildsten_dynamic_viscosity_pa_s", "missing")
        })?;
    let bildsten_prefactor = run
        .provenance
        .bildsten_dimensionless_prefactor
        .ok_or_else(|| {
            reduced_decay_bridge_refusal("provenance.bildsten_dimensionless_prefactor", "missing")
        })?;
    let omega_rad_s = (4.0 * run.parameters.gravity_m_per_s2
        / (run.parameters.radius_m * theta_rad.sin()))
    .sqrt();
    let published_rolling_w = rolling_coefficient
        * run.parameters.mass_kg
        * run.parameters.gravity_m_per_s2
        * run.parameters.radius_m
        * theta_rad.cos()
        * omega_rad_s;
    let bildsten_boundary_layer_w = bildsten_prefactor
        * BILDSTEN_PUBLISHED_POWER_COEFFICIENT
        * (viscosity_pa_s * density_kg_per_m3).sqrt()
        * run.parameters.gravity_m_per_s2.powf(1.25)
        * run.parameters.radius_m.powf(2.75)
        * theta_rad.powf(-1.25);
    if !omega_rad_s.is_finite()
        || !published_rolling_w.is_finite()
        || !bildsten_boundary_layer_w.is_finite()
        || published_rolling_w <= 0.0
        || bildsten_boundary_layer_w <= 0.0
    {
        return Err(reduced_decay_bridge_refusal(
            "crop.powers",
            "interpolated source-bound rates are outside their finite positive domain",
        ));
    }
    Ok((
        omega_rad_s,
        ChannelPowers {
            dry_contour_w: 0.0,
            published_rolling_w,
            bildsten_boundary_layer_w,
        },
    ))
}

fn reduced_decay_sample_inputs(
    run: &ReducedDecayRun,
    profile: &ResolvedDiscProfile,
    mass: MassProperties,
    tail_horizon_s: f64,
    publish_start_offset_s: f64,
    published_horizon_s: f64,
    required_boundary_offset_s: Option<f64>,
    cx: &Cx<'_>,
) -> Result<Vec<RenderTrajectorySampleInput>, RenderTrajectoryError> {
    if !(publish_start_offset_s.is_finite()
        && publish_start_offset_s >= 0.0
        && publish_start_offset_s < tail_horizon_s)
    {
        return Err(reduced_decay_bridge_refusal(
            "publish_start_offset_s",
            "must be finite and inside the retained tail",
        ));
    }
    let reconstructed_horizon_s = tail_horizon_s - publish_start_offset_s;
    let horizon_tolerance_s = 32.0
        * f64::EPSILON
        * published_horizon_s
            .abs()
            .max(reconstructed_horizon_s.abs())
            .max(1.0);
    if !(published_horizon_s.is_finite()
        && published_horizon_s > 0.0
        && reconstructed_horizon_s.is_finite()
        && (published_horizon_s - reconstructed_horizon_s).abs() <= horizon_tolerance_s)
    {
        return Err(reduced_decay_bridge_refusal(
            "published_horizon_s",
            "must be finite, positive, and agree with the retained source window within bounded roundoff",
        ));
    }
    let (time_origin_s, source_samples) =
        reduced_decay_tail_samples(run, tail_horizon_s, required_boundary_offset_s)?;
    let publish_origin_s = time_origin_s + publish_start_offset_s;
    let required_boundary_time_s =
        required_boundary_offset_s.map(|offset_s| time_origin_s + offset_s);
    let mut inputs = Vec::new();
    inputs
        .try_reserve_exact(source_samples.len())
        .map_err(|_| RenderTrajectoryError::Capacity {
            artifact: "reduced-decay render trajectory samples",
            requested: source_samples.len(),
        })?;
    let energy_slope_j_per_rad =
        1.5 * run.parameters.mass_kg * run.parameters.gravity_m_per_s2 * run.parameters.radius_m;
    let inertia = mass.principal_inertia_body();
    let mut precession_phase_rad = 0.0_f64;
    let mut spin_phase_rad = 0.0_f64;
    let mut center_xy_m = (0.0_f64, 0.0_f64);
    let mut previous_velocity_world_m_per_s: Option<Vec3> = None;
    let mut previous_sample: Option<&ReducedDecaySample> = None;
    let mut previous_published_velocity_world_m_per_s: Option<Vec3> = None;
    for (index, sample) in source_samples.iter().enumerate() {
        if index % TRAJECTORY_ADMISSION_CHECKPOINT_SAMPLES == 0 {
            cx.checkpoint()
                .map_err(|_| RenderTrajectoryError::Cancelled)?;
        }
        let spin_rate_rad_per_s = -sample.omega_rad_s * sample.theta_rad.cos();
        if let Some(previous) = previous_sample {
            let dt_s = sample.time_s - previous.time_s;
            let previous_spin_rate = -previous.omega_rad_s * previous.theta_rad.cos();
            precession_phase_rad = wrapped_phase(
                precession_phase_rad + 0.5 * (previous.omega_rad_s + sample.omega_rad_s) * dt_s,
            );
            spin_phase_rad = wrapped_phase(
                spin_phase_rad + 0.5 * (previous_spin_rate + spin_rate_rad_per_s) * dt_s,
            );
        }
        let orientation =
            reduced_decay_orientation(precession_phase_rad, sample.theta_rad, spin_phase_rad)?;
        let theta_rate_rad_per_s = -sample.powers.total_w() / energy_slope_j_per_rad;
        let tilt_axis_world =
            Vec3::new(-precession_phase_rad.sin(), precession_phase_rad.cos(), 0.0);
        let symmetry_axis_world = orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let omega_world_rad_per_s = Vec3::new(0.0, 0.0, sample.omega_rad_s)
            .add(tilt_axis_world.scale(theta_rate_rad_per_s))
            .add(symmetry_axis_world.scale(spin_rate_rad_per_s));
        let omega_body_rad_per_s = orientation.rotate_world_to_body(omega_world_rad_per_s);
        let angular_momentum_body = Vec3::new(
            inertia.x * omega_body_rad_per_s.x,
            inertia.y * omega_body_rad_per_s.y,
            inertia.z * omega_body_rad_per_s.z,
        );
        let provisional = profile_state_at_ground_contact(
            &profile.chart,
            profile.mass_properties,
            orientation,
            Vec3::ZERO,
            angular_momentum_body,
            cx,
        )
        .map_err(|error| {
            reduced_decay_bridge_refusal("profile_state_at_ground_contact", error.to_string())
        })?;
        let provisional_contact = profile_contact_geometry(
            &profile.chart,
            profile.mass_properties,
            provisional.pose(),
            cx,
        )
        .map_err(|error| {
            reduced_decay_bridge_refusal("profile_contact_geometry", error.to_string())
        })?;
        let velocity_world_m_per_s = omega_world_rad_per_s
            .cross(provisional_contact.contact.radius_world_m)
            .scale(-1.0);
        if let (Some(previous), Some(previous_velocity)) =
            (previous_sample, previous_velocity_world_m_per_s)
        {
            let dt_s = sample.time_s - previous.time_s;
            center_xy_m.0 += 0.5 * (previous_velocity.x + velocity_world_m_per_s.x) * dt_s;
            center_xy_m.1 += 0.5 * (previous_velocity.y + velocity_world_m_per_s.y) * dt_s;
        }
        let linear_momentum_world = velocity_world_m_per_s.scale(mass.mass());
        let grounded = profile_state_at_ground_contact(
            &profile.chart,
            profile.mass_properties,
            orientation,
            linear_momentum_world,
            angular_momentum_body,
            cx,
        )
        .map_err(|error| {
            reduced_decay_bridge_refusal("profile_state_at_ground_contact", error.to_string())
        })?;
        let grounded_position = grounded.pose().position_world();
        let pose = Pose::new(
            Vec3::new(center_xy_m.0, center_xy_m.1, grounded_position.z),
            orientation,
        )
        .map_err(|error| {
            reduced_decay_bridge_refusal("translated_ground_pose", error.to_string())
        })?;
        let state = RigidBodyState::new(pose, linear_momentum_world, angular_momentum_body)
            .map_err(|error| reduced_decay_bridge_refusal("grounded_state", error.to_string()))?;
        let contact =
            profile_contact_geometry(&profile.chart, profile.mass_properties, state.pose(), cx)
                .map_err(|error| {
                    reduced_decay_bridge_refusal("profile_contact_geometry", error.to_string())
                })?;
        let precession_acceleration_rad_per_s2 = -0.5 * sample.omega_rad_s * sample.theta_rad.cos()
            / sample.theta_rad.sin()
            * theta_rate_rad_per_s;
        let qois = DerivedEulerQois {
            inclination_rad: sample.theta_rad,
            precession_rad_per_s: sample.omega_rad_s,
            spin_rad_per_s: spin_rate_rad_per_s,
            precession_acceleration_rad_per_s2,
        };
        let derived =
            DerivedEulerQois::from_state(state, mass, precession_acceleration_rad_per_s2)?;
        if !scalar_close(qois.inclination_rad, derived.inclination_rad)
            || !scalar_close(qois.precession_rad_per_s, derived.precession_rad_per_s)
            || !scalar_close(qois.spin_rad_per_s, derived.spin_rad_per_s)
        {
            return Err(reduced_decay_bridge_refusal(
                "phase_rate_qoi_consistency",
                format!("sample {index} pose/twist does not recover declared Euler rates"),
            ));
        }
        if sample.time_s < publish_origin_s {
            previous_velocity_world_m_per_s = Some(velocity_world_m_per_s);
            previous_sample = Some(sample);
            continue;
        }
        let first_published_sample = inputs.is_empty();
        let (rolling_work_j, gas_work_j, interval_start_time_s) = if first_published_sample {
            (0.0, 0.0, sample.time_s)
        } else if let Some(previous) = previous_sample {
            (
                -(sample.work.dry_contour_j + sample.work.published_rolling_j
                    - previous.work.dry_contour_j
                    - previous.work.published_rolling_j),
                -(sample.work.bildsten_boundary_layer_j - previous.work.bildsten_boundary_layer_j),
                previous.time_s,
            )
        } else {
            return Err(reduced_decay_bridge_refusal(
                "published_interval",
                "a noninitial published sample lost its predecessor",
            ));
        };
        let energy_defect_j = energy_slope_j_per_rad * run.parameters.initial_theta_rad
            - sample.energy_j
            - sample.work.total_j();
        let output_time_s = if first_published_sample {
            0.0
        } else if required_boundary_time_s
            .is_some_and(|boundary_s| sample.time_s.to_bits() == boundary_s.to_bits())
        {
            required_boundary_offset_s.expect("boundary time has an offset")
                - publish_start_offset_s
        } else if index + 1 == source_samples.len()
            && (time_origin_s > 0.0 || publish_start_offset_s > 0.0)
        {
            published_horizon_s
        } else {
            sample.time_s - publish_origin_s
        };
        let output_interval_start_time_s = if first_published_sample {
            0.0
        } else if required_boundary_time_s
            .is_some_and(|boundary_s| interval_start_time_s.to_bits() == boundary_s.to_bits())
        {
            required_boundary_offset_s.expect("boundary time has an offset")
                - publish_start_offset_s
        } else {
            interval_start_time_s - publish_origin_s
        };
        let retained_velocity_world_m_per_s =
            state.center_of_mass_velocity_world(mass).map_err(|error| {
                reduced_decay_bridge_refusal("grounded_velocity", error.to_string())
            })?;
        let mean_normal_reaction_n = if first_published_sample {
            0.0
        } else {
            let previous_velocity = previous_published_velocity_world_m_per_s.ok_or_else(|| {
                reduced_decay_bridge_refusal(
                    "kinematically_implied_mean_normal_reaction_n",
                    "a noninitial published sample lost its predecessor velocity",
                )
            })?;
            let dt_s = output_time_s - output_interval_start_time_s;
            if !dt_s.is_finite() || dt_s <= 0.0 {
                return Err(reduced_decay_bridge_refusal(
                    "kinematically_implied_mean_normal_reaction_n",
                    format!("sample {index} has invalid published interval duration {dt_s:.17e} s"),
                ));
            }
            let vertical_acceleration_m_per_s2 =
                (retained_velocity_world_m_per_s.z - previous_velocity.z) / dt_s;
            let reaction_n =
                mass.mass() * (run.parameters.gravity_m_per_s2 + vertical_acceleration_m_per_s2);
            if !reaction_n.is_finite() || reaction_n < 0.0 {
                return Err(reduced_decay_bridge_refusal(
                    "kinematically_implied_mean_normal_reaction_n",
                    format!(
                        "sample {index} requires inadmissible interval-mean reaction {reaction_n:.17e} N"
                    ),
                ));
            }
            reaction_n
        };
        let channels = ChannelOwnership {
            gravity: ChannelWrench::default(),
            contact: ChannelWrench::default(),
            rolling: ChannelWrench {
                force_world_n: Vec3::ZERO,
                torque_world_nm: Vec3::ZERO,
                work_j: rolling_work_j,
            },
            base: ChannelWrench::default(),
            gas: ChannelWrench {
                force_world_n: Vec3::ZERO,
                torque_world_nm: Vec3::ZERO,
                work_j: gas_work_j,
            },
        };
        inputs.push(RenderTrajectorySampleInput {
            interval_start_time_s: output_interval_start_time_s,
            time_s: output_time_s,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            center_of_mass_world_m: state.pose().position_world(),
            orientation_body_to_world: orientation.components(),
            linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
            angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
            symmetry_axis_world,
            contact_branch: RenderContactBranch::Closed,
            contact_geometry: Some(RenderContactGeometry {
                point_world_m: contact.contact.point_world_m,
                normal_world: Vec3::new(0.0, 0.0, 1.0),
                support_feature: RenderSupportFeature::ProfileFeature(
                    contact.support_source_feature,
                ),
            }),
            signed_gap_m: contact.contact.gap_m,
            interval_contact_active: !first_published_sample,
            interval_normal_force_n: mean_normal_reaction_n,
            contact_transitions: Vec::new(),
            base_mode: Some(RenderBaseModeState {
                displacement_m: 0.0,
                velocity_m_per_s: 0.0,
            }),
            channels,
            mechanical_energy_j: sample.energy_j,
            energy_defect_j,
            qois,
            disposition: if index + 1 == source_samples.len() {
                RenderSampleDisposition::HorizonCensored
            } else {
                RenderSampleDisposition::Continue
            },
            terminal_event: None,
        });
        previous_published_velocity_world_m_per_s = Some(retained_velocity_world_m_per_s);
        previous_velocity_world_m_per_s = Some(velocity_world_m_per_s);
        previous_sample = Some(sample);
    }
    cx.checkpoint()
        .map_err(|_| RenderTrajectoryError::Cancelled)?;
    Ok(inputs)
}

fn reduced_decay_orientation(
    precession_phase_rad: f64,
    theta_rad: f64,
    spin_phase_rad: f64,
) -> Result<UnitQuaternion, RenderTrajectoryError> {
    UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), precession_phase_rad)
        .and_then(|orientation| orientation.right_exp(Vec3::new(0.0, theta_rad, 0.0)))
        .and_then(|orientation| orientation.right_exp(Vec3::new(0.0, 0.0, spin_phase_rad)))
        .map_err(|error| reduced_decay_bridge_refusal("orientation", error.to_string()))
}

fn wrapped_phase(phase_rad: f64) -> f64 {
    phase_rad.rem_euclid(TAU)
}

fn reduced_decay_model_identity(run: &ReducedDecayRun) -> ContentHash {
    let mut hasher = DomainHasher::new(REDUCED_DECAY_MODEL_IDENTITY_DOMAIN);
    hasher.update(&REDUCED_DECAY_RENDER_BRIDGE_VERSION.to_le_bytes());
    hash_text(&mut hasher, REDUCED_DECAY_PHASE_CONVENTION);
    hash_text(&mut hasher, REDUCED_DECAY_CONTACT_REACTION_CONVENTION);
    hash_text(&mut hasher, run.provenance.model_id);
    hash_text(&mut hasher, &run.provenance.small_angle_oracle_source_id);
    hash_text(&mut hasher, run.provenance.model_authority);
    hash_text(&mut hasher, run.provenance.physical_validation);
    for value in [
        run.parameters.mass_kg,
        run.parameters.radius_m,
        run.parameters.gravity_m_per_s2,
        run.parameters.initial_theta_rad,
        run.parameters.validity_cutoff_theta_rad,
        run.parameters.timestep_s,
        run.provenance
            .published_rolling_coefficient_mu
            .unwrap_or(f64::NAN),
        run.provenance
            .bildsten_density_kg_per_m3
            .unwrap_or(f64::NAN),
        run.provenance
            .bildsten_dynamic_viscosity_pa_s
            .unwrap_or(f64::NAN),
        run.provenance
            .bildsten_dimensionless_prefactor
            .unwrap_or(f64::NAN),
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.update(&run.parameters.maximum_steps.to_le_bytes());
    if let Some(specimen) = &run.provenance.literature_specimen {
        hash_text(&mut hasher, &specimen.source_id);
        for value in [
            specimen.diameter_m,
            specimen.thickness_m,
            specimen.mass_kg,
            specimen.outer_fillet_radius_m,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn reduced_decay_configuration_identity(
    run: &ReducedDecayRun,
    model_identity: ContentHash,
    profile_identity: ContentHash,
    chart_identity: ContentHash,
    mass_identity: ContentHash,
) -> ContentHash {
    let mut hasher = DomainHasher::new(REDUCED_DECAY_CONFIGURATION_IDENTITY_DOMAIN);
    hasher.update(model_identity.as_bytes());
    hasher.update(profile_identity.as_bytes());
    hasher.update(chart_identity.as_bytes());
    hasher.update(mass_identity.as_bytes());
    hasher.update(
        &u64::try_from(run.samples.len())
            .unwrap_or(u64::MAX)
            .to_le_bytes(),
    );
    for sample in &run.samples {
        for value in [
            sample.time_s,
            sample.theta_rad,
            sample.omega_rad_s,
            sample.energy_j,
            sample.powers.dry_contour_w,
            sample.powers.published_rolling_w,
            sample.powers.bildsten_boundary_layer_w,
            sample.work.dry_contour_j,
            sample.work.published_rolling_j,
            sample.work.bildsten_boundary_layer_j,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.update(&run.energy_closure_residual_j.to_bits().to_le_bytes());
    hasher.update(&[match run.terminal {
        ReducedDecayTerminal::ValidityCutoff => 1,
        ReducedDecayTerminal::StepBudgetExhausted => 2,
    }]);
    hasher.finalize()
}

fn hash_text(hasher: &mut DomainHasher, value: &str) {
    hasher.update(&u64::try_from(value.len()).unwrap_or(u64::MAX).to_le_bytes());
    hasher.update(value.as_bytes());
}

fn production_model_identity(
    model: &ProductionCouplingModel,
    profile_identity: ContentHash,
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.euler-disc.production-render-model.v1");
    hash_text(&mut hasher, &model.identity.case_id);
    hash_text(&mut hasher, &model.identity.configuration_id);
    hash_text(&mut hasher, &model.identity.world_frame_id);
    hasher.update(profile_identity.as_bytes());
    let gravity = model.gravity.acceleration_world();
    for value in [gravity.x, gravity.y, gravity.z] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    let (base_model_id, base_configuration_id) = model.base_port.identity_parts();
    hash_text(&mut hasher, base_model_id);
    hash_text(&mut hasher, base_configuration_id);
    hash_text(&mut hasher, model.tangential_adapter.adapter_id());
    hash_text(&mut hasher, model.tangential_adapter.source_id());
    hasher.finalize()
}

fn production_configuration_identity(
    model: &ProductionCouplingModel,
    model_identity: ContentHash,
) -> ContentHash {
    let mut hasher =
        DomainHasher::new("org.frankensim.euler-disc.production-render-configuration.v1");
    hasher.update(model_identity.as_bytes());
    hash_text(&mut hasher, &model.identity.configuration_id);
    hasher.finalize()
}

fn production_disc_work_residual_j(
    model: &ProductionCouplingModel,
    receipt: &ProductionCouplingReceipt,
    sample: usize,
) -> Result<f64, RenderTrajectoryError> {
    production_wrench_work_residual_j(
        model,
        &receipt.rigid_step,
        receipt.total_force_world_n,
        receipt.total_moment_about_com_world_n_m,
        sample,
    )
}

fn production_open_disc_work_residual_j(
    model: &ProductionCouplingModel,
    receipt: &ProductionOpenFlightReceipt,
    sample: usize,
) -> Result<f64, RenderTrajectoryError> {
    production_wrench_work_residual_j(
        model,
        &receipt.rigid_step,
        receipt.total_force_world_n,
        receipt.total_moment_about_com_world_n_m,
        sample,
    )
}

fn production_wrench_work_residual_j(
    model: &ProductionCouplingModel,
    rigid_step: &fs_mbd::StepReceipt,
    force_world_n: Vec3,
    moment_about_com_world_n_m: Vec3,
    sample: usize,
) -> Result<f64, RenderTrajectoryError> {
    let before = rigid_step.state_before;
    let after = rigid_step.state_after;
    let duration_s = rigid_step.duration_seconds;
    let velocity_before = before
        .center_of_mass_velocity_world(model.disc_mass_properties)
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    let velocity_after = after
        .center_of_mass_velocity_world(model.disc_mass_properties)
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    let omega_before_body = model
        .disc_mass_properties
        .angular_velocity_body_checked(before.angular_momentum_body())
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    let omega_after_body = model
        .disc_mass_properties
        .angular_velocity_body_checked(after.angular_momentum_body())
        .map_err(|error| RenderTrajectoryError::DerivedState(error.to_string()))?;
    let velocity_mid = velocity_before.add(velocity_after).scale(0.5);
    let omega_mid_body = omega_before_body.add(omega_after_body).scale(0.5);
    // `fs-mbd::Wrench` holds torque constant in body axes, whereas force is
    // held constant in world axes. Reconstruct that exact discrete convention;
    // averaging two world-frame angular velocities would silently use a
    // different torque history as the body rotates.
    let torque_body = before
        .pose()
        .orientation()
        .rotate_world_to_body(moment_about_com_world_n_m);
    let wrench_work_j =
        duration_s * (force_world_n.dot(velocity_mid) + torque_body.dot(omega_mid_body));
    let residual = rigid_step.diagnostics_after.mechanical_energy
        - rigid_step.diagnostics_before.mechanical_energy
        - wrench_work_j;
    if residual.is_finite() {
        Ok(residual)
    } else {
        Err(RenderTrajectoryError::NonFinite {
            sample: Some(sample),
            field: "production-prefix disc work residual",
        })
    }
}

fn production_event_disposition(
    termination: &ProductionEventTrajectoryTermination,
) -> RenderSampleDisposition {
    match termination {
        ProductionEventTrajectoryTermination::StepLimitReached { .. } => {
            RenderSampleDisposition::HorizonCensored
        }
        ProductionEventTrajectoryTermination::Refused { error, .. } => {
            RenderSampleDisposition::NumericalRefusal(
                RenderNumericalRefusalReason::BackendSpecific(production_refusal_code(error)),
            )
        }
    }
}

fn production_prefix_disposition(
    termination: &SmoothContactTrajectoryTermination,
) -> RenderSampleDisposition {
    match termination {
        SmoothContactTrajectoryTermination::StepLimitReached { .. } => {
            RenderSampleDisposition::HorizonCensored
        }
        SmoothContactTrajectoryTermination::Refused { error, .. } => {
            RenderSampleDisposition::NumericalRefusal(
                RenderNumericalRefusalReason::BackendSpecific(production_refusal_code(error)),
            )
        }
    }
}

const fn production_refusal_code(error: &ProductionCouplingError) -> u32 {
    match error {
        ProductionCouplingError::CheckpointMismatch => 1,
        ProductionCouplingError::CheckpointVersionMismatch { .. } => 2,
        ProductionCouplingError::CheckpointIntegrityMismatch => 3,
        ProductionCouplingError::InputIdentityMismatch { .. } => 4,
        ProductionCouplingError::InvalidInput { .. } => 5,
        ProductionCouplingError::ProfileContact(_) => 6,
        ProductionCouplingError::ProfileModelMassMismatch => 7,
        ProductionCouplingError::ResolvedProfileMassMismatch => 8,
        ProductionCouplingError::ResolvedProfileIdentityMismatch { .. } => 9,
        ProductionCouplingError::Patch(_) => 10,
        ProductionCouplingError::CurvatureUnavailable => 11,
        ProductionCouplingError::UnsupportedMechanism { .. } => 12,
        ProductionCouplingError::UnsupportedLineNormalContact => 13,
        ProductionCouplingError::Normal(_) => 14,
        ProductionCouplingError::Tangential(_) => 15,
        ProductionCouplingError::Rolling(_) => 16,
        ProductionCouplingError::ExteriorAir(_) => 17,
        ProductionCouplingError::AirFilm(_) => 18,
        ProductionCouplingError::GasChannelMismatch => 19,
        ProductionCouplingError::ExteriorCandidateUnavailable => 20,
        ProductionCouplingError::Base(_) => 21,
        ProductionCouplingError::Dynamics(_) => 22,
        ProductionCouplingError::SurfaceExcitation(_) => 23,
        ProductionCouplingError::ModalBase(_) => 24,
        ProductionCouplingError::BaseBackendMismatch => 25,
        ProductionCouplingError::BaseStaticPreloadUnsupported => 26,
        ProductionCouplingError::StaticPreloadMismatch { .. } => 27,
    }
}

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
    if metadata.channel_availability.contact
        && metadata.channel_availability.normal_force_sampling
            == RenderNormalForceSampling::Unavailable
    {
        return Err(RenderTrajectoryError::InvalidChannelAvailability(
            "contact requires non-unavailable normal_force_sampling",
        ));
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
    canonical_quaternion: bool,
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
    let accepted_subsecond_rebase_residue = declared_end_s.is_finite()
        && declared_end_s.abs() < 1.0
        && input.time_s - declared_end_s <= SUBSECOND_REBASED_CLOCK_TOLERANCE_S;
    if input.time_s > maximum_end_s && !accepted_subsecond_rebase_residue {
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
    let orientation =
        checked_unit_quaternion(input.orientation_body_to_world, index, canonical_quaternion)?;
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
    canonical: bool,
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
    let admitted = if canonical {
        UnitQuaternion::from_canonical_components(components)
    } else {
        UnitQuaternion::new(components[0], components[1], components[2], components[3])
    };
    admitted.map_err(|error| RenderTrajectoryError::InvalidRigidState(index, error.to_string()))
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
    if metadata.channel_availability.normal_force_sampling == RenderNormalForceSampling::Unavailable
        && input.interval_normal_force_n != 0.0
    {
        return Err(RenderTrajectoryError::UnavailableChannelHasData {
            sample: index,
            channel: "normal_force_sampling",
        });
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
