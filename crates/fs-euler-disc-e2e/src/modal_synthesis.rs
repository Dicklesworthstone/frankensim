//! Deterministic damped-modal synthesis for the Euler-disc cinematic pipeline.
//!
//! Each mode is a physically normalized scalar coordinate governed by
//! `m q'' + 2 zeta m omega q' + m omega^2 q = F`. Generalized component forces
//! are held constant over one audio frame (zero-order hold); declared impulses
//! are applied as exact velocity jumps at the frame's left boundary. The
//! runtime emits dry disc/glass/base stems and internal modal-energy diagnostics.
//!
//! This is a physically informed auditory model, not calibrated structural
//! acoustics. Representative presets are not measured eigenfrequencies, modal
//! energy is not radiated energy, and cross-component participation is not a
//! coupled structural or fluid-structure solve.

use core::{f64::consts::TAU, fmt};

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::cinematic::DeclaredAcousticCalibrationReceipt;
use fs_evidence::cinematic_sound::{
    MAX_SOUND_MODE_DAMPING_RATIO, MAX_SOUND_MODE_PARTICIPATION, MAX_SOUND_MODE_RADIATION_GAIN,
    MAX_SOUND_MODES, MIN_SOUND_MODE_MASS_KG, SOUND_MASTER_SAMPLE_RATE_HZ,
    SOUND_MODE_NYQUIST_GUARD_FRACTION, SoundModalComponent, SoundMode, SoundModeParticipation,
    SoundSynthesisConfig,
};
use fs_exec::Cx;
use fs_math::{STRICT_CORE_GOLDEN_HASH, STRICT_CORE_SEMANTICS_VERSION, det};

/// Exact version of the sampled modal algorithm and checkpoint semantics.
pub const MODAL_SYNTHESIS_ALGORITHM_VERSION: u32 = 2;
/// Maximum frames admitted in one transactional chunk (about 21.8 s at 48 kHz).
pub const MAX_MODAL_SYNTHESIS_CHUNK_FRAMES: usize = 1_048_576;
/// Maximum total frames admitted in one checkpoint lineage (one hour at 48 kHz).
pub const MAX_MODAL_SYNTHESIS_TOTAL_FRAMES: u64 = 172_800_000;
/// Maximum magnitude of a caller-supplied normalized spatial participation factor.
pub const MAX_MODAL_SPATIAL_PARTICIPATION: f64 = 1.0;
/// Maximum delay between cancellation polls during synthesis (1.33 ms at 48 kHz).
pub const MODAL_CANCELLATION_POLL_FRAMES: usize = 64;

const MODEL_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.modal-synthesis-model.v2";
const PRESET_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.modal-preset.v1";
const PRESET_COMPONENT_DOMAIN: &str = "org.frankensim.euler-cinematic.modal-preset-component.v1";
const PARAMETER_SET_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.modal-parameter-set.v1";
/// Exact version of the admitted Euler modal-parameter-set identity.
pub const EULER_MODAL_PARAMETER_SET_VERSION: u32 = 1;
/// Maximum UTF-8 disclosure size accepted by one modal parameter set.
pub const MAX_EULER_MODAL_DISCLOSURE_BYTES: usize = 1_024;
const COEFFICIENT_TAYLOR_TERMS: usize = 18;
const SMALL_STEP_MATRIX_NORM_LIMIT: f64 = 0.125;

/// Explicit execution and state limits. Exceeding a limit refuses the entire
/// proposed chunk; state and samples are never silently clamped.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalSynthesisBudget {
    /// Maximum frames accepted over one resumable checkpoint lineage.
    pub maximum_total_sample_frames: u64,
    /// Maximum frames accepted by one transactional call.
    pub maximum_chunk_sample_frames: usize,
    /// Maximum absolute displacement at sample and post-impulse boundaries [m].
    pub maximum_abs_displacement_m: f64,
    /// Maximum absolute velocity at sample and post-impulse boundaries [m/s].
    pub maximum_abs_velocity_m_per_s: f64,
    /// Maximum energy of any mode at sample and post-impulse boundaries [J].
    pub maximum_mode_energy_j: f64,
    /// Maximum summed energy at sample and post-impulse boundaries [J].
    pub maximum_total_energy_j: f64,
    /// Maximum absolute instantaneous dry sample in digital-full-scale units.
    pub maximum_abs_output_fs: f64,
}

impl ModalSynthesisBudget {
    /// A bounded reference-film budget. The caller still chooses the exact
    /// lineage length; no duration is inferred from a trajectory.
    #[must_use]
    pub const fn reference_film(maximum_total_sample_frames: u64) -> Self {
        let maximum_chunk_sample_frames =
            if maximum_total_sample_frames < MAX_MODAL_SYNTHESIS_CHUNK_FRAMES as u64 {
                maximum_total_sample_frames as usize
            } else {
                MAX_MODAL_SYNTHESIS_CHUNK_FRAMES
            };
        Self {
            maximum_total_sample_frames,
            maximum_chunk_sample_frames,
            maximum_abs_displacement_m: 0.05,
            maximum_abs_velocity_m_per_s: 2_000.0,
            maximum_mode_energy_j: 100_000.0,
            maximum_total_energy_j: 500_000.0,
            maximum_abs_output_fs: 16.0,
        }
    }
}

/// Generalized force or impulse coordinates for the three emitted components.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ModalComponentValues {
    /// Disc-body coordinate.
    pub disc: f64,
    /// Glass-plate coordinate.
    pub glass_plate: f64,
    /// Base-assembly coordinate.
    pub base_assembly: f64,
}

impl ModalComponentValues {
    /// Exact zero value.
    pub const ZERO: Self = Self {
        disc: 0.0,
        glass_plate: 0.0,
        base_assembly: 0.0,
    };
}

/// One audio-frame drive with explicit localized and distributed source classes.
///
/// Both force fields are held over the entire frame, and both impulse fields are
/// applied at the frame's left boundary first. [`ModalSpatialParticipation`]
/// modulates only the localized fields (contact and rolling sources). Distributed
/// fields (base damping and exterior-gas sources) retain the modes' declared
/// participation without contact-location modulation.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ModalDriveFrame {
    /// Piecewise-constant localized generalized component forces [N].
    pub localized_generalized_force_n: ModalComponentValues,
    /// Piecewise-constant distributed generalized component forces [N].
    pub distributed_generalized_force_n: ModalComponentValues,
    /// Exact localized generalized component impulses at the left boundary [N s].
    pub localized_boundary_impulse_n_s: ModalComponentValues,
    /// Exact distributed generalized component impulses at the left boundary [N s].
    pub distributed_boundary_impulse_n_s: ModalComponentValues,
}

/// Optional caller-supplied source-location modulation. The row-major slice is
/// indexed by `(frame, canonical_mode)`. A caller may either multiply localized
/// component drive by a normalized factor or supply the already-participated,
/// already-filtered localized modal drive directly. Deriving either form from
/// contact geometry is owned by the excitation-mapping stage. Distributed force
/// and impulse never use this modulation.
#[derive(Debug, Clone, Copy)]
pub enum ModalSpatialParticipation<'a> {
    /// Use each mode's declared static participation unchanged.
    Declared,
    /// Multiply by one normalized signed factor per frame and canonical mode.
    PerFrameModeFactors(&'a [f64]),
    /// Use already-participated localized modal drive in row-major
    /// `(frame, canonical_mode)` order.
    ///
    /// This is the resampling-safe path when source location changes within the
    /// anti-alias filter support: the mapper filters the localized force-factor
    /// product, rather than filtering its operands independently. The frame's
    /// localized component fields must be exactly zero; distributed component
    /// fields remain active through each mode's declared static participation.
    PreparticipatedLocalizedDrive {
        /// Piecewise-constant localized generalized force per mode [N].
        generalized_force_n: &'a [f64],
        /// Localized left-boundary generalized impulse per mode [N s].
        boundary_impulse_n_s: &'a [f64],
    },
}

/// Complete model input. Modes may arrive in any order; construction sorts by
/// unique nonzero `mode_id` before hashing or arithmetic.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalSynthesisModelInput {
    /// Exact sample rate. The v2 cinematic master requires 48 kHz.
    pub sample_rate_hz: u32,
    /// Unordered modes to admit and canonicalize.
    pub modes: Vec<SoundMode>,
    /// Explicit execution and state limits bound into model identity.
    pub budget: ModalSynthesisBudget,
}

/// Whether component forcing is one-hot or includes declared cross-routing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ModalCouplingClass {
    /// Every mode participates only in its own component coordinate.
    Independent,
    /// At least one mode has nonzero off-component source participation.
    DeclaredCrossParticipation,
}

/// Immutable, content-identified bank of exact sampled modal transitions.
#[derive(Debug, Clone)]
pub struct ModalSynthesisModel {
    identity: ContentHash,
    sample_rate_hz: u32,
    sample_period_s: f64,
    modes: Vec<PreparedMode>,
    public_modes: Vec<SoundMode>,
    budget: ModalSynthesisBudget,
    coupling: ModalCouplingClass,
}

/// State of one canonical modal coordinate at an exact sample boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalModeState {
    /// Stable mode identifier.
    pub mode_id: u32,
    /// Displacement-normalized modal coordinate [m].
    pub displacement_m: f64,
    /// Modal velocity [m/s].
    pub velocity_m_per_s: f64,
}

/// Immutable resumable boundary. `next_sample_frame` is the index of the first
/// drive frame not yet consumed.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalSynthesisCheckpoint {
    model_identity: ContentHash,
    next_sample_frame: u64,
    states: Vec<ModalModeState>,
}

impl ModalSynthesisCheckpoint {
    /// Model identity to which this state is bound.
    #[must_use]
    pub const fn model_identity(&self) -> ContentHash {
        self.model_identity
    }

    /// Index of the next audio sample frame.
    #[must_use]
    pub const fn next_sample_frame(&self) -> u64 {
        self.next_sample_frame
    }

    /// Canonical per-mode states.
    #[must_use]
    pub fn states(&self) -> &[ModalModeState] {
        &self.states
    }
}

/// Dry output attributed to each model component at one sample boundary.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ModalStemFrame {
    /// Disc-body dry output [digital full scale].
    pub disc_fs: f64,
    /// Glass-plate dry output [digital full scale].
    pub glass_plate_fs: f64,
    /// Base-assembly dry output [digital full scale].
    pub base_assembly_fs: f64,
}

/// Energy of one mode at the chunk's successor boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModalModeEnergy {
    /// Stable mode identifier.
    pub mode_id: u32,
    /// Kinetic energy `0.5 m v^2` [J].
    pub kinetic_j: f64,
    /// Elastic energy `0.5 m omega^2 q^2` [J].
    pub elastic_j: f64,
    /// Kinetic plus elastic energy [J].
    pub total_j: f64,
}

/// Chunk-local diagnostics. Peak is the largest absolute instantaneous sample;
/// RMS is the population RMS over exactly `[start_sample_frame,end_sample_frame)`.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalSynthesisDiagnostics {
    /// First rendered sample index.
    pub start_sample_frame: u64,
    /// Exclusive successor sample index.
    pub end_sample_frame: u64,
    /// Absolute mixed-sample peak over this chunk [digital full scale].
    pub peak_abs_fs: f64,
    /// Population RMS of mixed end-boundary samples in this chunk.
    pub rms_fs: f64,
    /// Largest total internal modal energy at the start, post-impulse, or
    /// successor sample boundaries [J].
    pub maximum_total_modal_energy_j: f64,
    /// Successor-boundary energy of every canonical mode.
    pub final_mode_energies: Vec<ModalModeEnergy>,
    /// Declared routing class; this is not a structural-coupling claim.
    pub coupling: ModalCouplingClass,
}

/// Atomically published output and successor state for one bounded chunk.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalSynthesisChunk {
    /// Canonically summed mono dry signal [digital full scale].
    pub mixed_samples_fs: Vec<f64>,
    /// Component-attributed dry stems at each sample boundary.
    pub stem_frames: Vec<ModalStemFrame>,
    /// Total internal modal energy after each frame [J].
    pub total_modal_energy_j: Vec<f64>,
    /// Chunk-local peak/RMS/energy diagnostics.
    pub diagnostics: ModalSynthesisDiagnostics,
    /// State to use for the next transactional chunk.
    pub successor: ModalSynthesisCheckpoint,
}

/// Provenance tier attached to Euler auditory parameter sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ModalPresetAuthority {
    /// Plausible hand-authored values for rendering, not measurements or fits.
    RepresentativeUncalibrated,
    /// The caller declares that the parameters came from measurements.
    ///
    /// This is provenance, not authentication, calibration, validation, or a
    /// sound-pressure authority promotion.
    DeclaredMeasured,
}

/// Caller-supplied binding to a separately verified calibration declaration.
///
/// Admission checks only that `verification_identity` is nonzero and binds all
/// receipt fields into the parameter-set identity. This module does not execute
/// a verifier, authenticate either identity, inspect measurement bytes, or
/// establish that the calibration applies to the specimen and rig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CallerVerifiedCalibrationBinding {
    /// Structurally admitted external calibration declaration.
    pub receipt: DeclaredAcousticCalibrationReceipt,
    /// Caller-supplied identity of the separate verification result.
    pub verification_identity: ContentHash,
}

/// Complete input for one provenance-bearing Euler modal parameter set.
///
/// The embedded model input is admitted through
/// [`ModalSynthesisModel::try_new`], so mode canonicalization and every
/// numerical invariant stay owned by the production modal runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct EulerModalParameterSetInput {
    /// Parameter provenance; this never promotes sound authority by itself.
    pub authority: ModalPresetAuthority,
    /// Exact specimen artifact identity. Zero is invalid.
    pub specimen_identity: ContentHash,
    /// Exact support/rig artifact identity. Zero is invalid.
    pub rig_identity: ContentHash,
    /// Nonempty, bounded, single-line disclosure shown with the parameter set.
    pub disclosure: String,
    /// Optional caller assertion that an external calibration receipt was
    /// separately verified. Only `DeclaredMeasured` may carry it.
    pub calibration: Option<CallerVerifiedCalibrationBinding>,
    /// Complete sampled modal runtime input.
    pub model: ModalSynthesisModelInput,
}

/// Immutable admitted modal parameter set plus its prepared synthesis model.
///
/// The set identity binds provenance, specimen, rig, disclosure, the complete
/// prepared model identity, and every optional calibration-binding field.
#[derive(Debug, Clone)]
pub struct EulerModalParameterSet {
    identity: ContentHash,
    authority: ModalPresetAuthority,
    specimen_identity: ContentHash,
    rig_identity: ContentHash,
    disclosure: String,
    calibration: Option<CallerVerifiedCalibrationBinding>,
    model: ModalSynthesisModel,
}

impl EulerModalParameterSet {
    /// Validate provenance and prepare the embedded production modal model.
    pub fn try_admit(
        input: EulerModalParameterSetInput,
        cx: &Cx<'_>,
    ) -> Result<Self, EulerModalParameterSetError> {
        if is_zero_hash(input.specimen_identity) {
            return Err(EulerModalParameterSetError::InvalidIdentity("specimen"));
        }
        if is_zero_hash(input.rig_identity) {
            return Err(EulerModalParameterSetError::InvalidIdentity("rig"));
        }
        if input.disclosure.is_empty()
            || input.disclosure.trim() != input.disclosure
            || input.disclosure.len() > MAX_EULER_MODAL_DISCLOSURE_BYTES
            || input.disclosure.chars().any(char::is_control)
        {
            return Err(EulerModalParameterSetError::InvalidDisclosure);
        }
        if input
            .calibration
            .is_some_and(|binding| is_zero_hash(binding.verification_identity))
        {
            return Err(EulerModalParameterSetError::InvalidIdentity(
                "calibration-verification",
            ));
        }
        if input.authority == ModalPresetAuthority::RepresentativeUncalibrated
            && input.calibration.is_some()
        {
            return Err(EulerModalParameterSetError::UnexpectedCalibrationBinding);
        }

        let model = ModalSynthesisModel::try_new(input.model, cx)
            .map_err(EulerModalParameterSetError::InvalidModel)?;
        let identity = modal_parameter_set_identity(
            input.authority,
            input.specimen_identity,
            input.rig_identity,
            &input.disclosure,
            input.calibration,
            model.identity(),
        );
        Ok(Self {
            identity,
            authority: input.authority,
            specimen_identity: input.specimen_identity,
            rig_identity: input.rig_identity,
            disclosure: input.disclosure,
            calibration: input.calibration,
            model,
        })
    }

    /// Content identity of the complete admitted parameter set.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Declared parameter provenance, never an acoustic-authority promotion.
    #[must_use]
    pub const fn authority(&self) -> ModalPresetAuthority {
        self.authority
    }

    /// Exact specimen artifact identity.
    #[must_use]
    pub const fn specimen_identity(&self) -> ContentHash {
        self.specimen_identity
    }

    /// Exact support/rig artifact identity.
    #[must_use]
    pub const fn rig_identity(&self) -> ContentHash {
        self.rig_identity
    }

    /// Binding human-facing disclosure.
    #[must_use]
    pub fn disclosure(&self) -> &str {
        &self.disclosure
    }

    /// Optional caller-supplied external verification binding.
    #[must_use]
    pub const fn calibration(&self) -> Option<CallerVerifiedCalibrationBinding> {
        self.calibration
    }

    /// Prepared production modal model admitted from this exact set.
    #[must_use]
    pub const fn model(&self) -> &ModalSynthesisModel {
        &self.model
    }

    /// Consume the provenance wrapper and return its prepared model.
    #[must_use]
    pub fn into_model(self) -> ModalSynthesisModel {
        self.model
    }
}

/// Typed refusal from modal parameter-set provenance or model admission.
#[derive(Debug, Clone, PartialEq)]
pub enum EulerModalParameterSetError {
    /// A required content identity was all zeroes.
    InvalidIdentity(&'static str),
    /// Disclosure was empty, padded, multiline/control-bearing, or oversized.
    InvalidDisclosure,
    /// An uncalibrated representative set attempted to carry calibration.
    UnexpectedCalibrationBinding,
    /// The production modal runtime refused the sampled model.
    InvalidModel(ModalSynthesisError),
}

impl fmt::Display for EulerModalParameterSetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for EulerModalParameterSetError {}

/// Disc material choice for the representative Euler assembly preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RepresentativeDiscMaterial {
    /// Dense tungsten disc parameter set.
    Tungsten,
    /// Stainless-steel cone/disc parameter set.
    StainlessSteel,
}

/// Config-ready representative disc/glass/base modes with an explicit no-claim.
#[derive(Debug, Clone, PartialEq)]
pub struct RepresentativeModalPreset {
    identity: ContentHash,
    disc_material: RepresentativeDiscMaterial,
    modes: Vec<SoundMode>,
}

impl RepresentativeModalPreset {
    /// Content identity of the complete representative parameter set.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Disc material choice.
    #[must_use]
    pub const fn disc_material(&self) -> RepresentativeDiscMaterial {
        self.disc_material
    }

    /// Immutable canonical mode list.
    #[must_use]
    pub fn modes(&self) -> &[SoundMode] {
        &self.modes
    }

    /// Authority tier of every built-in value.
    #[must_use]
    pub const fn authority(&self) -> ModalPresetAuthority {
        ModalPresetAuthority::RepresentativeUncalibrated
    }

    /// Binding disclosure for UI, manifests, and reports.
    #[must_use]
    pub const fn disclosure(&self) -> &'static str {
        "representative auditory modes only; not measured eigenfrequencies, calibrated radiation, or SPL"
    }
}

/// Typed transactional refusal from model admission or chunk synthesis.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalSynthesisError {
    /// Execution scope requested cancellation before atomic publication.
    Cancelled,
    /// Only the exact reference-master rate is accepted by v2.
    InvalidSampleRate(u32),
    /// No mode was supplied.
    EmptyModes,
    /// The public mode ceiling was exceeded.
    TooManyModes(usize),
    /// One mode violates a named physical or numerical invariant.
    InvalidMode {
        /// Stable mode identifier, or zero when zero itself is invalid.
        mode_id: u32,
        /// Stable field/invariant name.
        field: &'static str,
    },
    /// Two modes use the same stable identifier.
    DuplicateModeId(u32),
    /// One explicit budget is invalid or exceeds a hard ceiling.
    InvalidBudget(&'static str),
    /// A zero-length synthesis transaction was requested.
    EmptyDrive,
    /// The drive exceeds the declared per-chunk frame budget.
    ChunkFrameBudgetExceeded {
        /// Requested drive-frame count.
        requested: usize,
        /// Admitted per-call maximum.
        limit: usize,
    },
    /// The successor would exceed the lineage frame budget.
    TotalFrameBudgetExceeded {
        /// Requested exclusive successor frame, or `u64::MAX` on overflow.
        requested: u64,
        /// Admitted lineage maximum.
        limit: u64,
    },
    /// Per-frame spatial factors do not match `frame_count * mode_count`.
    SpatialParticipationLength {
        /// Exact `frame_count * mode_count` length.
        expected: usize,
        /// Supplied factor count.
        actual: usize,
    },
    /// A drive coordinate was NaN or infinite.
    NonFiniteDrive {
        /// Chunk-local frame index.
        frame: usize,
        /// Stable coordinate name.
        field: &'static str,
    },
    /// One normalized position-participation factor is invalid.
    InvalidSpatialParticipation {
        /// Chunk-local frame index.
        frame: usize,
        /// Stable canonical mode ID.
        mode_id: u32,
    },
    /// Direct per-mode localized drive was combined with component-localized drive.
    ConflictingLocalizedDrive {
        /// Chunk-local frame index.
        frame: usize,
        /// Stable nonzero localized coordinate name.
        field: &'static str,
    },
    /// A checkpoint belongs to another complete modal model.
    CheckpointIdentityMismatch,
    /// Checkpoint mode ordering/identity or frame count is malformed.
    InvalidCheckpoint,
    /// Allocator refused an explicitly preflighted result/state capacity.
    Capacity {
        /// Result/state collection that could not reserve memory.
        artifact: &'static str,
        /// Exact requested element count.
        requested: usize,
    },
    /// An intermediate force, state, energy, or output became non-finite.
    NonFiniteResult {
        /// Absolute sample frame at which the refusal occurred.
        sample_frame: u64,
        /// Stable mode ID, or `None` for an aggregate result.
        mode_id: Option<u32>,
        /// Stable quantity name.
        field: &'static str,
    },
    /// A proposed state or output exceeded an explicit model limit.
    LimitExceeded {
        /// Absolute sample frame at which the refusal occurred.
        sample_frame: u64,
        /// Stable mode ID, or `None` for an aggregate result.
        mode_id: Option<u32>,
        /// Stable limit name.
        field: &'static str,
        /// Absolute proposed magnitude.
        magnitude: f64,
        /// Declared maximum magnitude.
        limit: f64,
    },
    /// An admitted sound configuration disagrees with this exact model.
    SoundConfigurationMismatch(&'static str),
}

impl fmt::Display for ModalSynthesisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ModalSynthesisError {}

#[derive(Debug, Clone, Copy)]
struct ModalTransition {
    a00: f64,
    a01: f64,
    a10: f64,
    a11: f64,
    /// Displacement increment per constant acceleration [s^2].
    gamma_q: f64,
    /// Velocity increment per constant acceleration [s].
    gamma_v: f64,
}

#[derive(Debug, Clone, Copy)]
struct PreparedMode {
    mode: SoundMode,
    stiffness_n_per_m: f64,
    transition: ModalTransition,
}

impl ModalSynthesisModel {
    /// Validate, canonicalize, prepare, and content-identify a modal model.
    pub fn try_new(
        mut input: ModalSynthesisModelInput,
        cx: &Cx<'_>,
    ) -> Result<Self, ModalSynthesisError> {
        checkpoint(cx)?;
        if input.sample_rate_hz != SOUND_MASTER_SAMPLE_RATE_HZ {
            return Err(ModalSynthesisError::InvalidSampleRate(input.sample_rate_hz));
        }
        validate_budget(input.budget)?;
        if input.modes.is_empty() {
            return Err(ModalSynthesisError::EmptyModes);
        }
        if input.modes.len() > MAX_SOUND_MODES {
            return Err(ModalSynthesisError::TooManyModes(input.modes.len()));
        }
        input.modes.sort_by_key(|mode| mode.mode_id);
        validate_modes(&input.modes, input.sample_rate_hz)?;
        let sample_period_s = 1.0 / f64::from(input.sample_rate_hz);
        let mut modes = Vec::new();
        modes
            .try_reserve_exact(input.modes.len())
            .map_err(|_| ModalSynthesisError::Capacity {
                artifact: "prepared modal modes",
                requested: input.modes.len(),
            })?;
        for mode in &input.modes {
            checkpoint(cx)?;
            let omega_rad_per_s = TAU * mode.frequency_hz;
            let stiffness_n_per_m = mode.modal_mass_kg * omega_rad_per_s * omega_rad_per_s;
            let transition = modal_transition(omega_rad_per_s, mode.damping_ratio, sample_period_s)
                .map_err(|error| match error {
                    ModalSynthesisError::InvalidMode { field, .. } => {
                        ModalSynthesisError::InvalidMode {
                            mode_id: mode.mode_id,
                            field,
                        }
                    }
                    other => other,
                })?;
            if !omega_rad_per_s.is_finite() || !stiffness_n_per_m.is_finite() {
                return Err(ModalSynthesisError::InvalidMode {
                    mode_id: mode.mode_id,
                    field: "derived stiffness",
                });
            }
            modes.push(PreparedMode {
                mode: *mode,
                stiffness_n_per_m,
                transition,
            });
        }
        let coupling = classify_coupling(&input.modes);
        let identity = model_identity(input.sample_rate_hz, &input.modes, input.budget, coupling);
        checkpoint(cx)?;
        Ok(Self {
            identity,
            sample_rate_hz: input.sample_rate_hz,
            sample_period_s,
            modes,
            public_modes: input.modes,
            budget: input.budget,
            coupling,
        })
    }

    /// Complete model identity, including deterministic-math semantics and limits.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact sample rate [Hz].
    #[must_use]
    pub const fn sample_rate_hz(&self) -> u32 {
        self.sample_rate_hz
    }

    /// Exact sample period [s].
    #[must_use]
    pub const fn sample_period_s(&self) -> f64 {
        self.sample_period_s
    }

    /// Canonically ordered public mode definitions.
    #[must_use]
    pub fn modes(&self) -> &[SoundMode] {
        &self.public_modes
    }

    /// Explicit state and execution budget bound into model identity.
    #[must_use]
    pub const fn budget(&self) -> ModalSynthesisBudget {
        self.budget
    }

    /// Declared component-force routing class.
    #[must_use]
    pub const fn coupling(&self) -> ModalCouplingClass {
        self.coupling
    }

    /// Derived quality factor `Q = 1/(2 zeta)`, or infinity for zero damping.
    #[must_use]
    pub fn quality_factor(&self, mode_id: u32) -> Option<f64> {
        let damping = self
            .public_modes
            .binary_search_by_key(&mode_id, |mode| mode.mode_id)
            .ok()
            .map(|index| self.public_modes[index].damping_ratio)?;
        Some(if damping == 0.0 {
            f64::INFINITY
        } else {
            0.5 / damping
        })
    }

    /// Verify that the admitted L6 sound configuration names this exact model,
    /// algorithm version, sample rate, and canonical mode summaries.
    pub fn validate_sound_configuration(
        &self,
        config: &SoundSynthesisConfig,
    ) -> Result<(), ModalSynthesisError> {
        if config.input().audio_clock.ticks_per_second_numerator() != self.sample_rate_hz
            || config.input().audio_clock.ticks_per_second_denominator() != 1
        {
            return Err(ModalSynthesisError::SoundConfigurationMismatch(
                "sample rate",
            ));
        }
        if config.input().sound_model.identity() != self.identity {
            return Err(ModalSynthesisError::SoundConfigurationMismatch(
                "sound model identity",
            ));
        }
        if config.input().sound_model.version() != MODAL_SYNTHESIS_ALGORITHM_VERSION {
            return Err(ModalSynthesisError::SoundConfigurationMismatch(
                "sound model version",
            ));
        }
        if config.input().modes.as_slice() != self.public_modes.as_slice() {
            return Err(ModalSynthesisError::SoundConfigurationMismatch(
                "canonical modes",
            ));
        }
        Ok(())
    }

    /// Construct the all-zero state at sample boundary zero.
    pub fn initial_checkpoint(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ModalSynthesisCheckpoint, ModalSynthesisError> {
        checkpoint(cx)?;
        let mut states = Vec::new();
        states
            .try_reserve_exact(self.modes.len())
            .map_err(|_| ModalSynthesisError::Capacity {
                artifact: "initial modal states",
                requested: self.modes.len(),
            })?;
        states.extend(self.modes.iter().map(|mode| ModalModeState {
            mode_id: mode.mode.mode_id,
            displacement_m: 0.0,
            velocity_m_per_s: 0.0,
        }));
        checkpoint(cx)?;
        Ok(ModalSynthesisCheckpoint {
            model_identity: self.identity,
            next_sample_frame: 0,
            states,
        })
    }

    /// Render one transactional, zero-order-held drive chunk. The input
    /// checkpoint is immutable; any refusal returns neither samples nor state.
    pub fn synthesize_chunk(
        &self,
        checkpoint_state: &ModalSynthesisCheckpoint,
        drive: &[ModalDriveFrame],
        spatial_participation: ModalSpatialParticipation<'_>,
        cx: &Cx<'_>,
    ) -> Result<ModalSynthesisChunk, ModalSynthesisError> {
        self.synthesize_chunk_with_checkpoint(
            checkpoint_state,
            drive,
            spatial_participation,
            &mut || checkpoint(cx),
        )
    }

    fn synthesize_chunk_with_checkpoint(
        &self,
        checkpoint_state: &ModalSynthesisCheckpoint,
        drive: &[ModalDriveFrame],
        spatial_participation: ModalSpatialParticipation<'_>,
        checkpoint_fn: &mut impl FnMut() -> Result<(), ModalSynthesisError>,
    ) -> Result<ModalSynthesisChunk, ModalSynthesisError> {
        checkpoint_fn()?;
        self.validate_checkpoint(checkpoint_state)?;
        if drive.is_empty() {
            return Err(ModalSynthesisError::EmptyDrive);
        }
        if drive.len() > self.budget.maximum_chunk_sample_frames {
            return Err(ModalSynthesisError::ChunkFrameBudgetExceeded {
                requested: drive.len(),
                limit: self.budget.maximum_chunk_sample_frames,
            });
        }
        let requested_u64 = u64::try_from(drive.len()).map_err(|_| {
            ModalSynthesisError::TotalFrameBudgetExceeded {
                requested: u64::MAX,
                limit: self.budget.maximum_total_sample_frames,
            }
        })?;
        let successor_frame = checkpoint_state
            .next_sample_frame
            .checked_add(requested_u64)
            .ok_or(ModalSynthesisError::TotalFrameBudgetExceeded {
                requested: u64::MAX,
                limit: self.budget.maximum_total_sample_frames,
            })?;
        if successor_frame > self.budget.maximum_total_sample_frames {
            return Err(ModalSynthesisError::TotalFrameBudgetExceeded {
                requested: successor_frame,
                limit: self.budget.maximum_total_sample_frames,
            });
        }
        let spatial_len = drive.len().checked_mul(self.modes.len()).ok_or(
            ModalSynthesisError::SpatialParticipationLength {
                expected: usize::MAX,
                actual: spatial_slice_len(spatial_participation),
            },
        )?;
        match spatial_participation {
            ModalSpatialParticipation::Declared => {}
            ModalSpatialParticipation::PerFrameModeFactors(factors) => {
                validate_spatial_length(spatial_len, factors.len())?;
            }
            ModalSpatialParticipation::PreparticipatedLocalizedDrive {
                generalized_force_n,
                boundary_impulse_n_s,
            } => {
                validate_spatial_length(spatial_len, generalized_force_n.len())?;
                validate_spatial_length(spatial_len, boundary_impulse_n_s.len())?;
            }
        }
        preflight_drive(drive, spatial_participation, &self.modes, checkpoint_fn)?;
        checkpoint_fn()?;

        let mut states = Vec::new();
        states
            .try_reserve_exact(checkpoint_state.states.len())
            .map_err(|_| ModalSynthesisError::Capacity {
                artifact: "successor modal states",
                requested: checkpoint_state.states.len(),
            })?;
        states.extend_from_slice(&checkpoint_state.states);
        let mut mixed_samples_fs = Vec::new();
        let mut stem_frames = Vec::new();
        let mut total_modal_energy_j = Vec::new();
        for (artifact, output) in [
            ("mixed modal samples", &mut mixed_samples_fs),
            ("modal energy trace", &mut total_modal_energy_j),
        ] {
            output
                .try_reserve_exact(drive.len())
                .map_err(|_| ModalSynthesisError::Capacity {
                    artifact,
                    requested: drive.len(),
                })?;
        }
        stem_frames
            .try_reserve_exact(drive.len())
            .map_err(|_| ModalSynthesisError::Capacity {
                artifact: "modal component stems",
                requested: drive.len(),
            })?;

        let initial_total = self.total_energy(&states, checkpoint_state.next_sample_frame)?;
        let mut maximum_total_modal_energy_j = initial_total;
        let mut peak_abs_fs = 0.0_f64;
        let mut rms = ScaledSumSquares::new();
        for (local_frame, frame) in drive.iter().enumerate() {
            if local_frame % MODAL_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let absolute_frame = checkpoint_state.next_sample_frame + local_frame as u64;
            let mut stem_sums = [NeumaierSum::new(); 3];
            let mut kicked_energy_sum = NeumaierSum::new();
            let mut energy_sum = NeumaierSum::new();
            for (mode_index, (mode, state)) in self.modes.iter().zip(&mut states).enumerate() {
                let modal_drive_index = local_frame * self.modes.len() + mode_index;
                let (localized_force_n, localized_impulse_n_s) = localized_modal_drive(
                    spatial_participation,
                    frame,
                    mode.mode.source_participation,
                    modal_drive_index,
                );
                let distributed_force_n = participation_dot(
                    mode.mode.source_participation,
                    frame.distributed_generalized_force_n,
                );
                let distributed_impulse_n_s = participation_dot(
                    mode.mode.source_participation,
                    frame.distributed_boundary_impulse_n_s,
                );
                let generalized_force_n = localized_force_n + distributed_force_n;
                let generalized_impulse_n_s = localized_impulse_n_s + distributed_impulse_n_s;
                if !generalized_force_n.is_finite() || !generalized_impulse_n_s.is_finite() {
                    return Err(ModalSynthesisError::NonFiniteResult {
                        sample_frame: absolute_frame,
                        mode_id: Some(mode.mode.mode_id),
                        field: "participated drive",
                    });
                }
                let kicked_velocity =
                    state.velocity_m_per_s + generalized_impulse_n_s / mode.mode.modal_mass_kg;
                let kicked_state = ModalModeState {
                    mode_id: state.mode_id,
                    displacement_m: state.displacement_m,
                    velocity_m_per_s: kicked_velocity,
                };
                validate_proposed_state(
                    (kicked_state.displacement_m, kicked_state.velocity_m_per_s),
                    mode,
                    self.budget,
                    absolute_frame,
                )?;
                let kicked_energy = mode_energy(mode, kicked_state, absolute_frame)?;
                if kicked_energy.total_j > self.budget.maximum_mode_energy_j {
                    return Err(ModalSynthesisError::LimitExceeded {
                        sample_frame: absolute_frame,
                        mode_id: Some(mode.mode.mode_id),
                        field: "mode energy after boundary impulse",
                        magnitude: kicked_energy.total_j,
                        limit: self.budget.maximum_mode_energy_j,
                    });
                }
                kicked_energy_sum.add(kicked_energy.total_j);
                let acceleration = generalized_force_n / mode.mode.modal_mass_kg;
                let proposed = advance_mode(
                    mode.transition,
                    state.displacement_m,
                    kicked_velocity,
                    acceleration,
                );
                validate_proposed_state(proposed, mode, self.budget, absolute_frame)?;
                state.displacement_m = proposed.0;
                state.velocity_m_per_s = proposed.1;
                let energy = mode_energy(mode, *state, absolute_frame)?;
                if energy.total_j > self.budget.maximum_mode_energy_j {
                    return Err(ModalSynthesisError::LimitExceeded {
                        sample_frame: absolute_frame,
                        mode_id: Some(mode.mode.mode_id),
                        field: "mode energy",
                        magnitude: energy.total_j,
                        limit: self.budget.maximum_mode_energy_j,
                    });
                }
                energy_sum.add(energy.total_j);
                let radiated = mode.mode.radiation_gain_fs_s_per_m * state.velocity_m_per_s;
                validate_output(
                    radiated,
                    self.budget.maximum_abs_output_fs,
                    absolute_frame,
                    Some(mode.mode.mode_id),
                    "mode radiation output",
                )?;
                stem_sums[component_index(mode.mode.component)].add(radiated);
            }
            let kicked_total_energy = kicked_energy_sum.total();
            validate_total_energy(kicked_total_energy, self.budget, absolute_frame)?;
            maximum_total_modal_energy_j = maximum_total_modal_energy_j.max(kicked_total_energy);
            let total_energy = energy_sum.total();
            validate_total_energy(total_energy, self.budget, absolute_frame)?;
            maximum_total_modal_energy_j = maximum_total_modal_energy_j.max(total_energy);
            let stems = ModalStemFrame {
                disc_fs: stem_sums[0].total(),
                glass_plate_fs: stem_sums[1].total(),
                base_assembly_fs: stem_sums[2].total(),
            };
            for (field, value) in [
                ("disc stem", stems.disc_fs),
                ("glass plate stem", stems.glass_plate_fs),
                ("base assembly stem", stems.base_assembly_fs),
            ] {
                validate_output(
                    value,
                    self.budget.maximum_abs_output_fs,
                    absolute_frame,
                    None,
                    field,
                )?;
            }
            let mut mixed = NeumaierSum::new();
            mixed.add(stems.disc_fs);
            mixed.add(stems.glass_plate_fs);
            mixed.add(stems.base_assembly_fs);
            let mixed = mixed.total();
            validate_output(
                mixed,
                self.budget.maximum_abs_output_fs,
                absolute_frame,
                None,
                "mixed output",
            )?;
            peak_abs_fs = peak_abs_fs.max(mixed.abs());
            rms.add(mixed);
            mixed_samples_fs.push(mixed);
            stem_frames.push(stems);
            total_modal_energy_j.push(total_energy);
        }
        checkpoint_fn()?;

        let final_mode_energies = self.mode_energies(&states, successor_frame)?;
        let rms_fs = rms.root_mean_square(drive.len());
        if !rms_fs.is_finite() {
            return Err(ModalSynthesisError::NonFiniteResult {
                sample_frame: successor_frame,
                mode_id: None,
                field: "chunk RMS",
            });
        }
        let successor = ModalSynthesisCheckpoint {
            model_identity: self.identity,
            next_sample_frame: successor_frame,
            states,
        };
        Ok(ModalSynthesisChunk {
            mixed_samples_fs,
            stem_frames,
            total_modal_energy_j,
            diagnostics: ModalSynthesisDiagnostics {
                start_sample_frame: checkpoint_state.next_sample_frame,
                end_sample_frame: successor_frame,
                peak_abs_fs,
                rms_fs,
                maximum_total_modal_energy_j,
                final_mode_energies,
                coupling: self.coupling,
            },
            successor,
        })
    }

    fn validate_checkpoint(
        &self,
        checkpoint_state: &ModalSynthesisCheckpoint,
    ) -> Result<(), ModalSynthesisError> {
        if checkpoint_state.model_identity != self.identity {
            return Err(ModalSynthesisError::CheckpointIdentityMismatch);
        }
        if checkpoint_state.next_sample_frame > self.budget.maximum_total_sample_frames
            || checkpoint_state.states.len() != self.modes.len()
        {
            return Err(ModalSynthesisError::InvalidCheckpoint);
        }
        for (mode, state) in self.modes.iter().zip(&checkpoint_state.states) {
            if state.mode_id != mode.mode.mode_id {
                return Err(ModalSynthesisError::InvalidCheckpoint);
            }
            validate_proposed_state(
                (state.displacement_m, state.velocity_m_per_s),
                mode,
                self.budget,
                checkpoint_state.next_sample_frame,
            )?;
            let energy = mode_energy(mode, *state, checkpoint_state.next_sample_frame)?;
            if energy.total_j > self.budget.maximum_mode_energy_j {
                return Err(ModalSynthesisError::InvalidCheckpoint);
            }
        }
        let total =
            self.total_energy(&checkpoint_state.states, checkpoint_state.next_sample_frame)?;
        if total > self.budget.maximum_total_energy_j {
            return Err(ModalSynthesisError::InvalidCheckpoint);
        }
        Ok(())
    }

    fn total_energy(
        &self,
        states: &[ModalModeState],
        sample_frame: u64,
    ) -> Result<f64, ModalSynthesisError> {
        let mut total = NeumaierSum::new();
        for (mode, state) in self.modes.iter().zip(states) {
            total.add(mode_energy(mode, *state, sample_frame)?.total_j);
        }
        let total = total.total();
        validate_total_energy(total, self.budget, sample_frame)?;
        Ok(total)
    }

    fn mode_energies(
        &self,
        states: &[ModalModeState],
        sample_frame: u64,
    ) -> Result<Vec<ModalModeEnergy>, ModalSynthesisError> {
        let mut energies = Vec::new();
        energies.try_reserve_exact(self.modes.len()).map_err(|_| {
            ModalSynthesisError::Capacity {
                artifact: "final modal energies",
                requested: self.modes.len(),
            }
        })?;
        for (mode, state) in self.modes.iter().zip(states) {
            energies.push(mode_energy(mode, *state, sample_frame)?);
        }
        Ok(energies)
    }
}

/// Build the representative disc/glass/base parameter set for one disc material.
#[must_use]
pub fn representative_modal_preset(
    disc_material: RepresentativeDiscMaterial,
) -> RepresentativeModalPreset {
    let base_identity = preset_component_identity("metal-pa12-three-foot-base");
    let glass_identity = preset_component_identity("thick-glass-plate");
    let disc_identity = match disc_material {
        RepresentativeDiscMaterial::Tungsten => preset_component_identity("tungsten-disc"),
        RepresentativeDiscMaterial::StainlessSteel => {
            preset_component_identity("stainless-steel-cone-disc")
        }
    };
    let disc_modes = match disc_material {
        RepresentativeDiscMaterial::Tungsten => [
            (1_560.0, 0.014, 0.020, 1.00, 0.020),
            (4_320.0, 0.018, 0.011, 0.64, -0.011),
            (9_100.0, 0.026, 0.006, 0.38, 0.005),
        ],
        RepresentativeDiscMaterial::StainlessSteel => [
            (2_180.0, 0.018, 0.011, 1.00, 0.024),
            (5_810.0, 0.024, 0.006, 0.61, -0.013),
            (12_100.0, 0.032, 0.0038, 0.34, 0.006),
        ],
    };
    let mut modes = Vec::with_capacity(8);
    for (mode_id, (frequency_hz, damping_ratio, modal_mass_kg, participation, gain)) in
        [1_u32, 2, 3].into_iter().zip(disc_modes)
    {
        modes.push(SoundMode {
            mode_id,
            component: SoundModalComponent::Disc,
            frequency_hz,
            damping_ratio,
            modal_mass_kg,
            source_participation: SoundModeParticipation {
                disc: participation,
                glass_plate: 0.035 * participation,
                base_assembly: 0.0,
            },
            radiation_gain_fs_s_per_m: gain,
            material_identity: disc_identity,
            base_identity,
        });
    }
    for (mode_id, frequency_hz, damping_ratio, modal_mass_kg, gain) in [
        (101, 780.0, 0.022, 0.45, 0.018),
        (102, 2_460.0, 0.030, 0.22, -0.010),
        (103, 6_120.0, 0.041, 0.09, 0.0045),
    ] {
        modes.push(SoundMode {
            mode_id,
            component: SoundModalComponent::GlassPlate,
            frequency_hz,
            damping_ratio,
            modal_mass_kg,
            source_participation: SoundModeParticipation {
                disc: 0.08,
                glass_plate: 1.0,
                base_assembly: 0.12,
            },
            radiation_gain_fs_s_per_m: gain,
            material_identity: glass_identity,
            base_identity,
        });
    }
    for (mode_id, frequency_hz, damping_ratio, modal_mass_kg, gain) in [
        (201, 210.0, 0.055, 1.8, 0.012),
        (202, 630.0, 0.072, 0.8, -0.006),
    ] {
        modes.push(SoundMode {
            mode_id,
            component: SoundModalComponent::BaseAssembly,
            frequency_hz,
            damping_ratio,
            modal_mass_kg,
            source_participation: SoundModeParticipation {
                disc: 0.02,
                glass_plate: 0.15,
                base_assembly: 1.0,
            },
            radiation_gain_fs_s_per_m: gain,
            material_identity: base_identity,
            base_identity,
        });
    }
    let identity = preset_identity(disc_material, &modes);
    RepresentativeModalPreset {
        identity,
        disc_material,
        modes,
    }
}

fn validate_budget(budget: ModalSynthesisBudget) -> Result<(), ModalSynthesisError> {
    if budget.maximum_total_sample_frames == 0
        || budget.maximum_total_sample_frames > MAX_MODAL_SYNTHESIS_TOTAL_FRAMES
    {
        return Err(ModalSynthesisError::InvalidBudget(
            "maximum_total_sample_frames",
        ));
    }
    if budget.maximum_chunk_sample_frames == 0
        || budget.maximum_chunk_sample_frames > MAX_MODAL_SYNTHESIS_CHUNK_FRAMES
        || u64::try_from(budget.maximum_chunk_sample_frames).ok()
            > Some(budget.maximum_total_sample_frames)
    {
        return Err(ModalSynthesisError::InvalidBudget(
            "maximum_chunk_sample_frames",
        ));
    }
    for (field, value) in [
        (
            "maximum_abs_displacement_m",
            budget.maximum_abs_displacement_m,
        ),
        (
            "maximum_abs_velocity_m_per_s",
            budget.maximum_abs_velocity_m_per_s,
        ),
        ("maximum_mode_energy_j", budget.maximum_mode_energy_j),
        ("maximum_total_energy_j", budget.maximum_total_energy_j),
        ("maximum_abs_output_fs", budget.maximum_abs_output_fs),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(ModalSynthesisError::InvalidBudget(field));
        }
    }
    Ok(())
}

fn validate_modes(modes: &[SoundMode], sample_rate_hz: u32) -> Result<(), ModalSynthesisError> {
    let frequency_limit = f64::from(sample_rate_hz) * 0.5 * SOUND_MODE_NYQUIST_GUARD_FRACTION;
    let mut previous = None;
    for mode in modes {
        if mode.mode_id == 0 {
            return Err(ModalSynthesisError::InvalidMode {
                mode_id: 0,
                field: "mode_id",
            });
        }
        if previous == Some(mode.mode_id) {
            return Err(ModalSynthesisError::DuplicateModeId(mode.mode_id));
        }
        previous = Some(mode.mode_id);
        let field = if !mode.frequency_hz.is_finite()
            || mode.frequency_hz <= 0.0
            || mode.frequency_hz >= frequency_limit
        {
            Some("frequency_hz")
        } else if !mode.damping_ratio.is_finite()
            || mode.damping_ratio < 0.0
            || mode.damping_ratio > MAX_SOUND_MODE_DAMPING_RATIO
        {
            Some("damping_ratio")
        } else if !mode.modal_mass_kg.is_finite() || mode.modal_mass_kg < MIN_SOUND_MODE_MASS_KG {
            Some("modal_mass_kg")
        } else if [
            mode.source_participation.disc,
            mode.source_participation.glass_plate,
            mode.source_participation.base_assembly,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > MAX_SOUND_MODE_PARTICIPATION)
        {
            Some("source_participation")
        } else if !mode.radiation_gain_fs_s_per_m.is_finite()
            || mode.radiation_gain_fs_s_per_m.abs() > MAX_SOUND_MODE_RADIATION_GAIN
        {
            Some("radiation_gain_fs_s_per_m")
        } else if is_zero_hash(mode.material_identity) {
            Some("material_identity")
        } else if is_zero_hash(mode.base_identity) {
            Some("base_identity")
        } else {
            None
        };
        if let Some(field) = field {
            return Err(ModalSynthesisError::InvalidMode {
                mode_id: mode.mode_id,
                field,
            });
        }
    }
    Ok(())
}

fn modal_transition(
    omega: f64,
    damping_ratio: f64,
    dt: f64,
) -> Result<ModalTransition, ModalSynthesisError> {
    let x = omega * dt;
    let matrix_norm_bound = (1.0 + 2.0 * damping_ratio) * x;
    let transition = if matrix_norm_bound <= SMALL_STEP_MATRIX_NORM_LIMIT {
        modal_transition_taylor(omega, damping_ratio, dt)
    } else if damping_ratio < 1.0 {
        let nu = omega * det::sqrt((1.0 - damping_ratio) * (1.0 + damping_ratio));
        let phase = nu * dt;
        let exponential = det::exp(-damping_ratio * omega * dt);
        let s = exponential * dt * sinc(phase);
        let c = exponential * det::cos(phase);
        let gamma = damping_ratio * omega;
        let a00 = c + gamma * s;
        ModalTransition {
            a00,
            a01: s,
            a10: -omega * omega * s,
            a11: c - gamma * s,
            gamma_q: (1.0 - a00) / (omega * omega),
            gamma_v: s,
        }
    } else if damping_ratio == 1.0 {
        let exponential = det::exp(-x);
        let s = exponential * dt;
        let a00 = exponential * (1.0 + x);
        ModalTransition {
            a00,
            a01: s,
            a10: -omega * omega * s,
            a11: exponential * (1.0 - x),
            gamma_q: (-det::expm1(-x) - x * exponential) / (omega * omega),
            gamma_v: s,
        }
    } else {
        modal_transition_overdamped(omega, damping_ratio, dt)
    };
    if [
        transition.a00,
        transition.a01,
        transition.a10,
        transition.a11,
        transition.gamma_q,
        transition.gamma_v,
    ]
    .iter()
    .any(|value| !value.is_finite())
    {
        return Err(ModalSynthesisError::InvalidMode {
            mode_id: 0,
            field: "sampled transition",
        });
    }
    Ok(transition)
}

fn modal_transition_taylor(omega: f64, damping_ratio: f64, dt: f64) -> ModalTransition {
    let x = omega * dt;
    let matrix = [[0.0, 1.0], [-1.0, -2.0 * damping_ratio]];
    let mut phi = [[1.0, 0.0], [0.0, 1.0]];
    let mut term = phi;
    for order in 1..=COEFFICIENT_TAYLOR_TERMS {
        term = scale_matrix(multiply_matrix(term, matrix), x / order as f64);
        add_matrix(&mut phi, term);
    }
    let mut g = [0.0, 0.0];
    let mut g_term = [0.0, x];
    add_vector(&mut g, g_term);
    for order in 1..COEFFICIENT_TAYLOR_TERMS {
        g_term = scale_vector(
            multiply_matrix_vector(matrix, g_term),
            x / (order + 1) as f64,
        );
        add_vector(&mut g, g_term);
    }
    ModalTransition {
        a00: phi[0][0],
        a01: phi[0][1] / omega,
        a10: phi[1][0] * omega,
        a11: phi[1][1],
        gamma_q: g[0] / (omega * omega),
        gamma_v: g[1] / omega,
    }
}

fn modal_transition_overdamped(omega: f64, damping_ratio: f64, dt: f64) -> ModalTransition {
    let root = det::sqrt((damping_ratio - 1.0) * (damping_ratio + 1.0));
    let sum = damping_ratio + root;
    let slow = omega / sum;
    let fast = omega * sum;
    let separation = fast - slow;
    let slow_decay = det::exp(-slow * dt);
    let s = slow_decay * (-det::expm1(-separation * dt)) / separation;
    let one_minus_slow = -det::expm1(-slow * dt);
    // These algebraically equivalent forms avoid subtracting nearly equal
    // pole-weighted exponentials as damping approaches the critical branch.
    let a00 = slow_decay + slow * s;
    let a11 = slow_decay - fast * s;
    ModalTransition {
        a00,
        a01: s,
        a10: -omega * omega * s,
        a11,
        gamma_q: (one_minus_slow - slow * s) / (omega * omega),
        gamma_v: s,
    }
}

fn sinc(value: f64) -> f64 {
    let squared = value * value;
    if squared <= 1.0e-8 {
        1.0 + squared * (-1.0 / 6.0 + squared * (1.0 / 120.0 - squared / 5_040.0))
    } else {
        det::sin(value) / value
    }
}

fn advance_mode(
    transition: ModalTransition,
    displacement: f64,
    velocity: f64,
    constant_acceleration: f64,
) -> (f64, f64) {
    (
        transition.a00 * displacement
            + transition.a01 * velocity
            + transition.gamma_q * constant_acceleration,
        transition.a10 * displacement
            + transition.a11 * velocity
            + transition.gamma_v * constant_acceleration,
    )
}

fn preflight_drive(
    drive: &[ModalDriveFrame],
    spatial: ModalSpatialParticipation<'_>,
    modes: &[PreparedMode],
    checkpoint_fn: &mut impl FnMut() -> Result<(), ModalSynthesisError>,
) -> Result<(), ModalSynthesisError> {
    for (frame, values) in drive.iter().enumerate() {
        if frame % MODAL_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        for (field, value) in [
            (
                "localized disc force",
                values.localized_generalized_force_n.disc,
            ),
            (
                "localized glass plate force",
                values.localized_generalized_force_n.glass_plate,
            ),
            (
                "localized base assembly force",
                values.localized_generalized_force_n.base_assembly,
            ),
            (
                "distributed disc force",
                values.distributed_generalized_force_n.disc,
            ),
            (
                "distributed glass plate force",
                values.distributed_generalized_force_n.glass_plate,
            ),
            (
                "distributed base assembly force",
                values.distributed_generalized_force_n.base_assembly,
            ),
            (
                "localized disc impulse",
                values.localized_boundary_impulse_n_s.disc,
            ),
            (
                "localized glass plate impulse",
                values.localized_boundary_impulse_n_s.glass_plate,
            ),
            (
                "localized base assembly impulse",
                values.localized_boundary_impulse_n_s.base_assembly,
            ),
            (
                "distributed disc impulse",
                values.distributed_boundary_impulse_n_s.disc,
            ),
            (
                "distributed glass plate impulse",
                values.distributed_boundary_impulse_n_s.glass_plate,
            ),
            (
                "distributed base assembly impulse",
                values.distributed_boundary_impulse_n_s.base_assembly,
            ),
        ] {
            if !value.is_finite() {
                return Err(ModalSynthesisError::NonFiniteDrive { frame, field });
            }
        }
        let row = frame * modes.len();
        match spatial {
            ModalSpatialParticipation::Declared => {}
            ModalSpatialParticipation::PerFrameModeFactors(factors) => {
                for (mode_index, mode) in modes.iter().enumerate() {
                    let factor = factors[row + mode_index];
                    if !factor.is_finite() || factor.abs() > MAX_MODAL_SPATIAL_PARTICIPATION {
                        return Err(ModalSynthesisError::InvalidSpatialParticipation {
                            frame,
                            mode_id: mode.mode.mode_id,
                        });
                    }
                }
            }
            ModalSpatialParticipation::PreparticipatedLocalizedDrive {
                generalized_force_n,
                boundary_impulse_n_s,
            } => {
                for mode_index in 0..modes.len() {
                    for (field, value) in [
                        (
                            "preparticipated localized force",
                            generalized_force_n[row + mode_index],
                        ),
                        (
                            "preparticipated localized impulse",
                            boundary_impulse_n_s[row + mode_index],
                        ),
                    ] {
                        if !value.is_finite() {
                            return Err(ModalSynthesisError::NonFiniteDrive { frame, field });
                        }
                    }
                }
            }
        }
        if matches!(
            spatial,
            ModalSpatialParticipation::PreparticipatedLocalizedDrive { .. }
        ) {
            for (field, value) in [
                (
                    "localized disc force",
                    values.localized_generalized_force_n.disc,
                ),
                (
                    "localized glass plate force",
                    values.localized_generalized_force_n.glass_plate,
                ),
                (
                    "localized base assembly force",
                    values.localized_generalized_force_n.base_assembly,
                ),
                (
                    "localized disc impulse",
                    values.localized_boundary_impulse_n_s.disc,
                ),
                (
                    "localized glass plate impulse",
                    values.localized_boundary_impulse_n_s.glass_plate,
                ),
                (
                    "localized base assembly impulse",
                    values.localized_boundary_impulse_n_s.base_assembly,
                ),
            ] {
                if value != 0.0 {
                    return Err(ModalSynthesisError::ConflictingLocalizedDrive { frame, field });
                }
            }
        }
    }
    Ok(())
}

fn validate_spatial_length(expected: usize, actual: usize) -> Result<(), ModalSynthesisError> {
    if actual != expected {
        return Err(ModalSynthesisError::SpatialParticipationLength { expected, actual });
    }
    Ok(())
}

fn localized_modal_drive(
    spatial: ModalSpatialParticipation<'_>,
    frame: &ModalDriveFrame,
    participation: SoundModeParticipation,
    row_major_index: usize,
) -> (f64, f64) {
    match spatial {
        ModalSpatialParticipation::Declared => (
            participation_dot(participation, frame.localized_generalized_force_n),
            participation_dot(participation, frame.localized_boundary_impulse_n_s),
        ),
        ModalSpatialParticipation::PerFrameModeFactors(factors) => {
            let factor = factors[row_major_index];
            (
                participation_dot(participation, frame.localized_generalized_force_n) * factor,
                participation_dot(participation, frame.localized_boundary_impulse_n_s) * factor,
            )
        }
        ModalSpatialParticipation::PreparticipatedLocalizedDrive {
            generalized_force_n,
            boundary_impulse_n_s,
        } => (
            generalized_force_n[row_major_index],
            boundary_impulse_n_s[row_major_index],
        ),
    }
}

fn spatial_slice_len(spatial: ModalSpatialParticipation<'_>) -> usize {
    match spatial {
        ModalSpatialParticipation::Declared => 0,
        ModalSpatialParticipation::PerFrameModeFactors(factors) => factors.len(),
        ModalSpatialParticipation::PreparticipatedLocalizedDrive {
            generalized_force_n,
            ..
        } => generalized_force_n.len(),
    }
}

fn participation_dot(participation: SoundModeParticipation, values: ModalComponentValues) -> f64 {
    let mut sum = NeumaierSum::new();
    sum.add(participation.disc * values.disc);
    sum.add(participation.glass_plate * values.glass_plate);
    sum.add(participation.base_assembly * values.base_assembly);
    sum.total()
}

fn classify_coupling(modes: &[SoundMode]) -> ModalCouplingClass {
    if modes.iter().any(|mode| {
        [
            (SoundModalComponent::Disc, mode.source_participation.disc),
            (
                SoundModalComponent::GlassPlate,
                mode.source_participation.glass_plate,
            ),
            (
                SoundModalComponent::BaseAssembly,
                mode.source_participation.base_assembly,
            ),
        ]
        .into_iter()
        .any(|(component, value)| component != mode.component && value != 0.0)
    }) {
        ModalCouplingClass::DeclaredCrossParticipation
    } else {
        ModalCouplingClass::Independent
    }
}

fn component_index(component: SoundModalComponent) -> usize {
    match component {
        SoundModalComponent::Disc => 0,
        SoundModalComponent::GlassPlate => 1,
        SoundModalComponent::BaseAssembly => 2,
    }
}

fn validate_proposed_state(
    proposed: (f64, f64),
    mode: &PreparedMode,
    budget: ModalSynthesisBudget,
    sample_frame: u64,
) -> Result<(), ModalSynthesisError> {
    for (field, value, limit) in [
        (
            "modal displacement",
            proposed.0,
            budget.maximum_abs_displacement_m,
        ),
        (
            "modal velocity",
            proposed.1,
            budget.maximum_abs_velocity_m_per_s,
        ),
    ] {
        if !value.is_finite() {
            return Err(ModalSynthesisError::NonFiniteResult {
                sample_frame,
                mode_id: Some(mode.mode.mode_id),
                field,
            });
        }
        if value.abs() > limit {
            return Err(ModalSynthesisError::LimitExceeded {
                sample_frame,
                mode_id: Some(mode.mode.mode_id),
                field,
                magnitude: value.abs(),
                limit,
            });
        }
    }
    Ok(())
}

fn mode_energy(
    mode: &PreparedMode,
    state: ModalModeState,
    sample_frame: u64,
) -> Result<ModalModeEnergy, ModalSynthesisError> {
    let kinetic_j = 0.5 * mode.mode.modal_mass_kg * state.velocity_m_per_s * state.velocity_m_per_s;
    let elastic_j = 0.5 * mode.stiffness_n_per_m * state.displacement_m * state.displacement_m;
    let total_j = kinetic_j + elastic_j;
    if !kinetic_j.is_finite()
        || !elastic_j.is_finite()
        || !total_j.is_finite()
        || kinetic_j < 0.0
        || elastic_j < 0.0
        || total_j < 0.0
    {
        return Err(ModalSynthesisError::NonFiniteResult {
            sample_frame,
            mode_id: Some(mode.mode.mode_id),
            field: "modal energy",
        });
    }
    Ok(ModalModeEnergy {
        mode_id: mode.mode.mode_id,
        kinetic_j,
        elastic_j,
        total_j,
    })
}

fn validate_total_energy(
    total: f64,
    budget: ModalSynthesisBudget,
    sample_frame: u64,
) -> Result<(), ModalSynthesisError> {
    if !total.is_finite() || total < 0.0 {
        return Err(ModalSynthesisError::NonFiniteResult {
            sample_frame,
            mode_id: None,
            field: "total modal energy",
        });
    }
    if total > budget.maximum_total_energy_j {
        return Err(ModalSynthesisError::LimitExceeded {
            sample_frame,
            mode_id: None,
            field: "total modal energy",
            magnitude: total,
            limit: budget.maximum_total_energy_j,
        });
    }
    Ok(())
}

fn validate_output(
    value: f64,
    limit: f64,
    sample_frame: u64,
    mode_id: Option<u32>,
    field: &'static str,
) -> Result<(), ModalSynthesisError> {
    if !value.is_finite() {
        return Err(ModalSynthesisError::NonFiniteResult {
            sample_frame,
            mode_id,
            field,
        });
    }
    if value.abs() > limit {
        return Err(ModalSynthesisError::LimitExceeded {
            sample_frame,
            mode_id,
            field,
            magnitude: value.abs(),
            limit,
        });
    }
    Ok(())
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ModalSynthesisError> {
    cx.checkpoint().map_err(|_| ModalSynthesisError::Cancelled)
}

#[derive(Clone, Copy)]
struct NeumaierSum {
    sum: f64,
    correction: f64,
}

impl NeumaierSum {
    const fn new() -> Self {
        Self {
            sum: 0.0,
            correction: 0.0,
        }
    }

    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        if self.sum.abs() >= value.abs() {
            self.correction += (self.sum - next) + value;
        } else {
            self.correction += (value - next) + self.sum;
        }
        self.sum = next;
    }

    fn total(self) -> f64 {
        self.sum + self.correction
    }
}

struct ScaledSumSquares {
    scale: f64,
    sum_squares: f64,
}

impl ScaledSumSquares {
    const fn new() -> Self {
        Self {
            scale: 0.0,
            sum_squares: 1.0,
        }
    }

    fn add(&mut self, value: f64) {
        let magnitude = value.abs();
        if magnitude == 0.0 {
            return;
        }
        if self.scale < magnitude {
            let ratio = self.scale / magnitude;
            self.sum_squares = 1.0 + self.sum_squares * ratio * ratio;
            self.scale = magnitude;
        } else {
            let ratio = magnitude / self.scale;
            self.sum_squares += ratio * ratio;
        }
    }

    fn root_mean_square(&self, count: usize) -> f64 {
        if self.scale == 0.0 {
            0.0
        } else {
            self.scale * det::sqrt(self.sum_squares / count as f64)
        }
    }
}

fn multiply_matrix(left: [[f64; 2]; 2], right: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [
        [
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ],
        [
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ],
    ]
}

fn multiply_matrix_vector(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> [f64; 2] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    ]
}

fn scale_matrix(mut matrix: [[f64; 2]; 2], scale: f64) -> [[f64; 2]; 2] {
    for row in &mut matrix {
        for value in row {
            *value *= scale;
        }
    }
    matrix
}

fn scale_vector(vector: [f64; 2], scale: f64) -> [f64; 2] {
    [vector[0] * scale, vector[1] * scale]
}

fn add_matrix(target: &mut [[f64; 2]; 2], addend: [[f64; 2]; 2]) {
    for (target_row, addend_row) in target.iter_mut().zip(addend) {
        for (target_value, addend_value) in target_row.iter_mut().zip(addend_row) {
            *target_value += addend_value;
        }
    }
}

fn add_vector(target: &mut [f64; 2], addend: [f64; 2]) {
    target[0] += addend[0];
    target[1] += addend[1];
}

fn model_identity(
    sample_rate_hz: u32,
    modes: &[SoundMode],
    budget: ModalSynthesisBudget,
    coupling: ModalCouplingClass,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(160 + modes.len() * 128);
    push_u32(&mut bytes, MODAL_SYNTHESIS_ALGORITHM_VERSION);
    push_u32(&mut bytes, STRICT_CORE_SEMANTICS_VERSION);
    push_u64(&mut bytes, STRICT_CORE_GOLDEN_HASH);
    push_u32(&mut bytes, sample_rate_hz);
    bytes.push(match coupling {
        ModalCouplingClass::Independent => 1,
        ModalCouplingClass::DeclaredCrossParticipation => 2,
    });
    bytes.extend_from_slice(b"zoh-left-impulse-end-velocity-output-chunk-population-rms-v1");
    push_u64(&mut bytes, budget.maximum_total_sample_frames);
    push_u64(&mut bytes, budget.maximum_chunk_sample_frames as u64);
    for value in [
        budget.maximum_abs_displacement_m,
        budget.maximum_abs_velocity_m_per_s,
        budget.maximum_mode_energy_j,
        budget.maximum_total_energy_j,
        budget.maximum_abs_output_fs,
    ] {
        push_f64(&mut bytes, value);
    }
    push_u32(&mut bytes, modes.len() as u32);
    encode_modes(&mut bytes, modes);
    hash_domain(MODEL_IDENTITY_DOMAIN, &bytes)
}

fn preset_identity(material: RepresentativeDiscMaterial, modes: &[SoundMode]) -> ContentHash {
    let mut bytes = Vec::with_capacity(16 + modes.len() * 128);
    bytes.push(match material {
        RepresentativeDiscMaterial::Tungsten => 1,
        RepresentativeDiscMaterial::StainlessSteel => 2,
    });
    bytes.extend_from_slice(b"representative-uncalibrated-v1");
    encode_modes(&mut bytes, modes);
    hash_domain(PRESET_IDENTITY_DOMAIN, &bytes)
}

fn modal_parameter_set_identity(
    authority: ModalPresetAuthority,
    specimen_identity: ContentHash,
    rig_identity: ContentHash,
    disclosure: &str,
    calibration: Option<CallerVerifiedCalibrationBinding>,
    model_identity: ContentHash,
) -> ContentHash {
    let mut bytes = Vec::with_capacity(256 + disclosure.len());
    push_u32(&mut bytes, EULER_MODAL_PARAMETER_SET_VERSION);
    bytes.push(match authority {
        ModalPresetAuthority::RepresentativeUncalibrated => 1,
        ModalPresetAuthority::DeclaredMeasured => 2,
    });
    bytes.extend_from_slice(specimen_identity.as_bytes());
    bytes.extend_from_slice(rig_identity.as_bytes());
    push_u64(&mut bytes, disclosure.len() as u64);
    bytes.extend_from_slice(disclosure.as_bytes());
    bytes.extend_from_slice(model_identity.as_bytes());
    match calibration {
        None => bytes.push(0),
        Some(binding) => {
            bytes.push(1);
            bytes.extend_from_slice(binding.receipt.dataset_identity().as_bytes());
            bytes.extend_from_slice(binding.receipt.method_identity().as_bytes());
            bytes.extend_from_slice(binding.receipt.validity_identity().as_bytes());
            push_u32(&mut bytes, binding.receipt.version());
            bytes.extend_from_slice(binding.verification_identity.as_bytes());
        }
    }
    hash_domain(PARAMETER_SET_IDENTITY_DOMAIN, &bytes)
}

fn encode_modes(bytes: &mut Vec<u8>, modes: &[SoundMode]) {
    for mode in modes {
        push_u32(bytes, mode.mode_id);
        bytes.push(mode.component as u8);
        push_f64(bytes, mode.frequency_hz);
        push_f64(bytes, mode.damping_ratio);
        push_f64(bytes, mode.modal_mass_kg);
        push_f64(bytes, mode.source_participation.disc);
        push_f64(bytes, mode.source_participation.glass_plate);
        push_f64(bytes, mode.source_participation.base_assembly);
        push_f64(bytes, mode.radiation_gain_fs_s_per_m);
        bytes.extend_from_slice(mode.material_identity.as_bytes());
        bytes.extend_from_slice(mode.base_identity.as_bytes());
    }
}

fn preset_component_identity(label: &str) -> ContentHash {
    hash_domain(PRESET_COMPONENT_DOMAIN, label.as_bytes())
}

fn is_zero_hash(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_f64(bytes: &mut Vec<u8>, value: f64) {
    bytes.extend_from_slice(&value.to_bits().to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_test_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x4d4f_4441_4c5f_554e,
                    kernel_id: 0x4555_4c45,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    #[test]
    fn coefficient_branches_are_continuous_around_critical_damping() {
        let omega = TAU * 1_000.0;
        let dt = 1.0 / 48_000.0;
        let critical = modal_transition(omega, 1.0, dt).unwrap();
        for damping in [1.0 - 1.0e-10, 1.0 + 1.0e-10] {
            let nearby = modal_transition(omega, damping, dt).unwrap();
            for (actual, expected) in [
                (nearby.a00, critical.a00),
                (nearby.a01, critical.a01),
                (nearby.a10, critical.a10),
                (nearby.a11, critical.a11),
                (nearby.gamma_q, critical.gamma_q),
                (nearby.gamma_v, critical.gamma_v),
            ] {
                assert!((actual - expected).abs() <= 2.0e-10 * expected.abs().max(1.0));
            }
        }
    }

    #[test]
    fn small_step_taylor_matches_closed_underdamped_form() {
        let omega = TAU * 10.0;
        let damping = 0.02;
        let dt = 1.0 / 192_000.0;
        let taylor = modal_transition_taylor(omega, damping, dt);
        let nu = omega * det::sqrt((1.0 - damping) * (1.0 + damping));
        let exponential = det::exp(-damping * omega * dt);
        let s = exponential * det::sin(nu * dt) / nu;
        let c = exponential * det::cos(nu * dt);
        let a00 = c + damping * omega * s;
        let closed = ModalTransition {
            a00,
            a01: s,
            a10: -omega * omega * s,
            a11: c - damping * omega * s,
            gamma_q: (1.0 - a00) / (omega * omega),
            gamma_v: s,
        };
        for (actual, expected) in [
            (taylor.a00, closed.a00),
            (taylor.a01, closed.a01),
            (taylor.a10, closed.a10),
            (taylor.a11, closed.a11),
            (taylor.gamma_v, closed.gamma_v),
        ] {
            assert!((actual - expected).abs() <= 2.0e-13 * expected.abs().max(1.0));
        }
        let gamma_q_series = 0.5 * dt * dt - damping * omega * dt * dt * dt / 3.0
            + (4.0 * damping * damping - 1.0) * omega * omega * dt.powi(4) / 24.0;
        assert!(
            (taylor.gamma_q - gamma_q_series).abs() <= 2.0e-10 * gamma_q_series.abs(),
            "gamma_q={} series={gamma_q_series}",
            taylor.gamma_q,
        );
    }

    #[test]
    fn injected_mid_synthesis_cancellation_publishes_no_successor() {
        let test_identity = hash_domain("org.frankensim.test.modal-cancellation.v1", b"disc-mode");
        let mode = SoundMode {
            mode_id: 1,
            component: SoundModalComponent::Disc,
            frequency_hz: 800.0,
            damping_ratio: 0.02,
            modal_mass_kg: 0.2,
            source_participation: SoundModeParticipation {
                disc: 1.0,
                glass_plate: 0.0,
                base_assembly: 0.0,
            },
            radiation_gain_fs_s_per_m: 0.1,
            material_identity: test_identity,
            base_identity: test_identity,
        };
        let budget = ModalSynthesisBudget {
            maximum_total_sample_frames: 256,
            maximum_chunk_sample_frames: 256,
            maximum_abs_displacement_m: 1.0,
            maximum_abs_velocity_m_per_s: 1_000.0,
            maximum_mode_energy_j: 1.0e6,
            maximum_total_energy_j: 1.0e6,
            maximum_abs_output_fs: 1.0e6,
        };
        let model = with_test_cx(|cx| {
            ModalSynthesisModel::try_new(
                ModalSynthesisModelInput {
                    sample_rate_hz: SOUND_MASTER_SAMPLE_RATE_HZ,
                    modes: vec![mode],
                    budget,
                },
                cx,
            )
            .unwrap()
        });
        let initial = with_test_cx(|cx| model.initial_checkpoint(cx).unwrap());
        let original = initial.clone();
        let drive = vec![ModalDriveFrame::default(); 256];
        let mut polls = 0_usize;
        let result = model.synthesize_chunk_with_checkpoint(
            &initial,
            &drive,
            ModalSpatialParticipation::Declared,
            &mut || {
                polls += 1;
                if polls == 9 {
                    Err(ModalSynthesisError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(polls, 9, "cancellation must be injected during synthesis");
        assert_eq!(result, Err(ModalSynthesisError::Cancelled));
        assert_eq!(initial, original);
    }
}
