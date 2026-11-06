//! Audible delivery of physical observer-pressure signals.
//!
//! [`PhysicalPressureSignal`](crate::structural_acoustics::PhysicalPressureSignal)
//! is the scientific product: its samples are pressure in pascals at a named
//! observer. A media player instead expects unitless digital-full-scale
//! samples. This module keeps that conversion explicit and replayable. It
//! never changes, normalizes, or relabels the source pressure signal; it emits
//! a separate listening master whose boundary window and gain are presentation
//! metadata. The window prevents a finite media cut through an already-ringing
//! pressure field from manufacturing a playback click at either endpoint.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_math::det;

use crate::audio_artifact::{
    AudioArtifactBudget, AudioArtifactError, AudioMeters, StereoSample, WavCodecReceipt,
    WavMetadata, WavSampleEncoding, encode_stereo_wav, measure_audio,
};
use fs_evidence::cinematic_sound::SOUND_MASTER_SAMPLE_RATE_HZ;
use crate::render_trajectory::RenderTrajectory;
use crate::structural_acoustics::{
    AcousticWorldReceiver, PhysicalPressureSignal, RetardedFarFieldObserverControls,
    RigidDiscBroadbandRadiationArtifact, StructuralBroadbandRadiationArtifact,
    StructuralBroadbandSourceStem, StructuralModalBasisError,
    StructuralResidualFlexibilityEstimateBasis, synthesize_retarded_far_field_world_observers,
    synthesize_rigid_disc_retarded_far_field_world_observers,
};

const PHYSICAL_LISTENING_MASTER_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-disc.physical-listening-master.v3";

/// Explicit policy for the nonphysical pressure-to-digital presentation step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureListeningMasterPolicy {
    /// Fixed pressure-to-digital calibration [FS / Pa]. It never depends on the
    /// rendered signal, so pressure differences remain audible differences.
    pub digital_gain_fs_per_pa: f64,
    /// Largest admitted four-times-oversampled peak in digital-full-scale units.
    /// Exceeding this ceiling refuses instead of silently normalizing or limiting.
    pub maximum_true_peak_fs: f64,
    /// Raised-cosine fade from exact zero at the media cut [sample frames].
    pub initial_fade_sample_frames: u32,
    /// Raised-cosine fade to exact zero at the media cut [sample frames].
    pub terminal_fade_sample_frames: u32,
}

impl PressureListeningMasterPolicy {
    /// Critique master using a fixed `512 FS / Pa` review calibration and a
    /// conventional -1 dBFS estimated true-peak ceiling. The power-of-two gain
    /// is exact, remains constant across runs, and can be overridden explicitly
    /// when a playback chain has a measured pressure-to-digital calibration.
    /// The short boundary windows exist only on the listening derivative; the
    /// separately identified physical pressure signals remain unchanged.
    pub const CRITIQUE: Self = Self {
        digital_gain_fs_per_pa: 512.0,
        maximum_true_peak_fs: 0.891_250_938_133_745_6,
        initial_fade_sample_frames: 960,
        terminal_fade_sample_frames: 240,
    };
}

/// Canonical audible derivative of one left/right physical observer pair.
#[derive(Clone, Debug, PartialEq)]
pub struct PhysicalPressureListeningMaster {
    wav_bytes: Vec<u8>,
    /// Canonical WAV codec receipt.
    pub wav: WavCodecReceipt,
    /// Decoded-sample meters after deterministic presentation gain.
    pub meters: AudioMeters,
    /// Left physical pressure-signal identity.
    pub left_pressure_identity: ContentHash,
    /// Right physical pressure-signal identity.
    pub right_pressure_identity: ContentHash,
    /// Applied unit conversion and presentation gain `[FS / Pa]`.
    pub digital_gain_fs_per_pa: f64,
    /// Applied gain in decibels relative to `1 FS / Pa`.
    pub digital_gain_db: f64,
    /// Largest absolute physical pressure over both observers [Pa].
    pub source_peak_abs_pressure_pa: f64,
    /// Identity binding both physical signals, policy, exact digital samples,
    /// codec receipt, and applied presentation gain.
    pub identity: ContentHash,
}

impl PhysicalPressureListeningMaster {
    /// Deterministically master and encode two physical observer signals.
    ///
    /// Pass the same signal twice for a valid dual-mono presentation. Distinct
    /// signals must share the same physical-field aggregate identities, start
    /// time, sample rate, and frame count; their radiation identities may
    /// differ because they represent different observer locations. An
    /// aggregate identity may represent one radiator or a deterministic SI
    /// pressure superposition of several simultaneously excited bodies.
    ///
    /// # Errors
    /// Refuses mismatched signals, silence, malformed policy, excessive gain,
    /// non-finite pressure, cancellation, metering, or codec failures.
    pub fn try_build(
        left: &PhysicalPressureSignal,
        right: &PhysicalPressureSignal,
        policy: PressureListeningMasterPolicy,
        metadata: &WavMetadata,
        budget: AudioArtifactBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, PhysicalListeningMasterError> {
        validate_pair(left, right)?;
        if !(policy.digital_gain_fs_per_pa > 0.0 && policy.digital_gain_fs_per_pa.is_finite())
            || !(policy.maximum_true_peak_fs > 0.0
                && policy.maximum_true_peak_fs < 1.0
                && policy.maximum_true_peak_fs.is_finite())
        {
            return Err(PhysicalListeningMasterError::InvalidPolicy);
        }
        validate_boundary_fades(policy, left.pressure_pa.len())?;

        let source_peak_abs_pressure_pa = left.peak_abs_pressure_pa.max(right.peak_abs_pressure_pa);
        if !(source_peak_abs_pressure_pa > 0.0 && source_peak_abs_pressure_pa.is_finite()) {
            return Err(PhysicalListeningMasterError::SilentPhysicalSignal);
        }

        let digital_gain_fs_per_pa = policy.digital_gain_fs_per_pa;
        let digital_gain_db = 20.0 * det::ln(digital_gain_fs_per_pa) / core::f64::consts::LN_10;
        if !digital_gain_db.is_finite() {
            return Err(PhysicalListeningMasterError::InvalidPolicy);
        }

        let samples = convert_pair(left, right, policy, digital_gain_fs_per_pa)?;
        let meters = measure_audio(&samples, budget, cx)?;
        let observed_peak = meters.sample_peak_fs.max(meters.true_peak_estimate_fs);
        let tolerance = 256.0 * f64::EPSILON * policy.maximum_true_peak_fs;
        if observed_peak > policy.maximum_true_peak_fs + tolerance {
            return Err(PhysicalListeningMasterError::PeakTargetExceeded {
                observed_fs: observed_peak,
                target_fs: policy.maximum_true_peak_fs,
            });
        }
        let (wav_bytes, wav) = encode_stereo_wav(
            &samples,
            left.sample_rate_hz,
            WavSampleEncoding::Pcm24,
            metadata,
            budget,
            cx,
        )?;
        let identity = listening_master_identity(
            left,
            right,
            policy,
            digital_gain_fs_per_pa,
            source_peak_abs_pressure_pa,
            wav,
            &samples,
        );
        Ok(Self {
            wav_bytes,
            wav,
            meters,
            left_pressure_identity: left.identity,
            right_pressure_identity: right.identity,
            digital_gain_fs_per_pa,
            digital_gain_db,
            source_peak_abs_pressure_pa,
            identity,
        })
    }

    /// Exact canonical PCM24 WAV bytes.
    #[must_use]
    pub fn wav_bytes(&self) -> &[u8] {
        &self.wav_bytes
    }
}

/// Existing broadband radiator consumed by the physical receiver audition.
#[derive(Clone, Copy, Debug)]
pub enum RetardedPressureSource<'a> {
    /// Flexible structural modes with their original radiation and damping.
    Structural(&'a StructuralBroadbandRadiationArtifact),
    /// Existing rigid-disc translational radiation estimate.
    RigidDisc(&'a RigidDiscBroadbandRadiationArtifact),
}

/// Physical inputs for one independently moving left/right receiver pair.
/// This composes the existing renderer; it never applies artistic panning.
pub struct RetardedPressureListeningInput<'a> {
    /// Radiation model that produced the supplied coefficient stem.
    pub source: RetardedPressureSource<'a>,
    /// Original body-frame source coefficients, before listener propagation.
    pub stem: &'a StructuralBroadbandSourceStem,
    /// Exact structural geometry identity shared with the source trajectory.
    pub basis: &'a StructuralResidualFlexibilityEstimateBasis,
    /// Immutable accepted source mechanics.
    pub trajectory: &'a RenderTrajectory,
    /// Left/right microphones, each fixed or independently moving.
    pub receivers: [AcousticWorldReceiver<'a>; 2],
    /// Explicit source/output bands, low-Mach source scope and work caps.
    /// The listening consumer accepts 48 kHz directly or 96 kHz with the
    /// existing antialias filter; the shared renderer supports other clocks.
    pub observation: RetardedFarFieldObserverControls,
}

/// Transactional physical moving-receiver audition and calibrated PCM24 WAV.
#[derive(Clone, Debug)]
pub struct RetardedPressureListeningArtifact {
    /// Common 48 kHz left/right pressure in pascals, before digital gain or fades.
    /// A 96 kHz observation is antialias-filtered onto these original master
    /// boundaries, retaining only its complete real filter-history interior.
    pub pressure: [PhysicalPressureSignal; 2],
    /// Existing fixed-calibration listening derivative of exactly these signals.
    pub listening_master: PhysicalPressureListeningMaster,
}

impl RetardedPressureListeningArtifact {
    /// Render the actual receiver field and encode its audible derivative.
    ///
    /// All source, receiver, stencil, Doppler-band and far-field admission is
    /// delegated to the shared retarded renderer. Moving receivers retain their
    /// finite histories and common reception clock. No pressure amplitude factor,
    /// gain normalization, extrapolated tail or microphone/head scattering is
    /// introduced. The result retains the original estimate-only authority.
    ///
    /// # Errors
    /// Refuses invalid physical inputs or budgets, unsupported listening clocks,
    /// incomplete filter interiors, cancellation, silence or listening overload.
    pub fn try_build(
        input: RetardedPressureListeningInput<'_>,
        policy: PressureListeningMasterPolicy,
        metadata: &WavMetadata,
        budget: AudioArtifactBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, RetardedPressureListeningError> {
        if input.observation.output_sample_rate_hz != SOUND_MASTER_SAMPLE_RATE_HZ
            && input.observation.output_sample_rate_hz != 2 * SOUND_MASTER_SAMPLE_RATE_HZ
        {
            return Err(PhysicalListeningMasterError::SignalMismatch {
                field: "retarded listening observation must be 48 or 96 kHz",
            }.into());
        }
        let signals = match input.source {
            RetardedPressureSource::Structural(source) => {
                synthesize_retarded_far_field_world_observers(
                    source, input.stem, input.basis, input.trajectory, &input.receivers,
                    input.observation, cx,
                )?
            }
            RetardedPressureSource::RigidDisc(source) => {
                synthesize_rigid_disc_retarded_far_field_world_observers(
                    source, input.stem, input.basis, input.trajectory, &input.receivers,
                    input.observation, cx,
                )?
            }
        };
        let [left, right]: [PhysicalPressureSignal; 2] = signals.try_into().map_err(|_| {
            StructuralModalBasisError::InvalidRequest {
                what: "retarded stereo renderer must return exactly two receivers",
            }
        })?;
        let pressure = if input.observation.output_sample_rate_hz == 2 * SOUND_MASTER_SAMPLE_RATE_HZ {
            [
                decimate_retarded_pressure_to_master(left, cx)?,
                decimate_retarded_pressure_to_master(right, cx)?,
            ]
        } else {
            [left, right]
        };
        let listening_master = PhysicalPressureListeningMaster::try_build(
            &pressure[0], &pressure[1], policy, metadata, budget, cx,
        )?;
        Ok(Self { pressure, listening_master })
    }
}

/// Typed refusal of either the physical radiation or its listening derivative.
#[derive(Debug)]
pub enum RetardedPressureListeningError {
    /// Original shared source/receiver radiation refusal.
    Radiation(StructuralModalBasisError),
    /// Antialias, calibration, metering or PCM24 publication refusal.
    Listening(PhysicalListeningMasterError),
}

impl From<StructuralModalBasisError> for RetardedPressureListeningError {
    fn from(error: StructuralModalBasisError) -> Self {
        Self::Radiation(error)
    }
}

impl From<PhysicalListeningMasterError> for RetardedPressureListeningError {
    fn from(error: PhysicalListeningMasterError) -> Self {
        Self::Listening(error)
    }
}

impl fmt::Display for RetardedPressureListeningError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Radiation(error) => write!(formatter, "retarded physical pressure: {error}"),
            Self::Listening(error) => write!(formatter, "retarded listening derivative: {error}"),
        }
    }
}

impl std::error::Error for RetardedPressureListeningError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Radiation(error) => Some(error),
            Self::Listening(error) => Some(error),
        }
    }
}

/// Typed refusal from physical-pressure listening-master construction.
#[derive(Clone, Debug, PartialEq)]
pub enum PhysicalListeningMasterError {
    /// The presentation policy is nonfinite or outside its safe domain.
    InvalidPolicy,
    /// Left and right signals disagree on a required shared physical/timeline
    /// field.
    SignalMismatch { field: &'static str },
    /// A pressure sample or cached source peak is nonfinite.
    NonFinitePressure { channel: &'static str, frame: usize },
    /// No nonzero physical pressure exists to audition.
    SilentPhysicalSignal,
    /// Audible presentation would require more gain than the caller admits.
    GainOutsidePolicy {
        required_db: f64,
        maximum_absolute_db: f64,
    },
    /// The fixed calibration exceeded the declared true-peak ceiling.
    PeakTargetExceeded { observed_fs: f64, target_fs: f64 },
    /// Original pressure antialias filter or its finite-history admission refused.
    PressureResampling(String),
    /// Canonical WAV/metering/cancellation refusal.
    Artifact(AudioArtifactError),
}

impl From<AudioArtifactError> for PhysicalListeningMasterError {
    fn from(error: AudioArtifactError) -> Self {
        Self::Artifact(error)
    }
}

impl fmt::Display for PhysicalListeningMasterError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PhysicalListeningMasterError {}

/// Antialias one 96 kHz physical receiver onto the original 48 kHz master grid.
///
/// Reuses the native half-band decimator, full-history interior and group-delay
/// compensation used by the cinematic fixed-microphone mix. Returned samples
/// remain pascals; their identity binds the original pressure, filter and clock.
///
/// # Errors
/// Refuses a different/off-grid clock, incomplete filter support, allocation,
/// original filter errors or cancellation. No missing history is padded.
pub fn decimate_retarded_pressure_to_master(
    mut signal: PhysicalPressureSignal,
    cx: &Cx<'_>,
) -> Result<PhysicalPressureSignal, PhysicalListeningMasterError> {
    let source_rate = f64::from(2 * SOUND_MASTER_SAMPLE_RATE_HZ);
    let coordinate = signal.start_time_s * source_rate;
    let rounded = coordinate.round();
    if signal.sample_rate_hz != (2 * SOUND_MASTER_SAMPLE_RATE_HZ)
        || !coordinate.is_finite()
        || rounded.abs() > 9_007_199_254_740_991.0
        || (coordinate - rounded).abs() > 1.0e-6
    {
        return Err(PhysicalListeningMasterError::PressureResampling(
            "retarded pressure must retain its explicit twice-master sample clock".into(),
        ));
    }
    // Each output is centered on an ORIGINAL master closing boundary. If the
    // common retarded horizon starts halfway through a master cell, retain
    // only the following complete cells rather than shifting their times.
    let skipped_input = rounded.rem_euclid(2.0) as usize;
    let input = signal.pressure_pa.get(skipped_input..).ok_or_else(|| {
        PhysicalListeningMasterError::PressureResampling("retarded pressure has no complete input group".into())
    })?;
    let mut filter = fs_couple::pcm_wav::decimate::Decimator::new(2, 1)
        .map_err(|error| PhysicalListeningMasterError::PressureResampling(error.to_string()))?;
    let delay = filter.delay_output_frames();
    let unsupported_prefix = (2.0 * delay) as usize;
    let retained = (input.len() / 2)
        .checked_sub(unsupported_prefix)
        .filter(|&count| count > 0)
        .ok_or_else(|| {
            PhysicalListeningMasterError::PressureResampling(
                "retarded pressure has no complete antialias-filter interior".into(),
            )
        })?;
    let start_time_s = signal.start_time_s
        + skipped_input as f64 / source_rate
        + (unsupported_prefix as f64 - delay) / f64::from(SOUND_MASTER_SAMPLE_RATE_HZ);
    let mut identity =
        DomainHasher::new("org.frankensim.euler-disc.retarded-pressure-decimation.v1");
    identity.update(signal.identity.as_bytes());
    identity.update(b"blackman-harris-integer-decimator-v1");
    identity.update(&SOUND_MASTER_SAMPLE_RATE_HZ.to_le_bytes());
    identity.update(&start_time_s.to_bits().to_le_bytes());
    identity.update(&delay.to_bits().to_le_bytes());
    identity.update(&(skipped_input as u64).to_le_bytes());
    identity.update(&(retained as u64).to_le_bytes());
    let mut output = Vec::new();
    output.try_reserve_exact(retained).map_err(|_| {
        PhysicalListeningMasterError::PressureResampling("retarded pressure decimation allocation refused".into())
    })?;
    let mut peak_abs_pressure_pa = 0.0_f64;
    for (frame, samples) in input.chunks_exact(2).enumerate() {
        if frame % 64 == 0 {
            cx.checkpoint()
                .map_err(|_| PhysicalListeningMasterError::Artifact(AudioArtifactError::Cancelled))?;
        }
        let value = filter.preview(samples).map_err(|error| PhysicalListeningMasterError::PressureResampling(error.to_string()))?[0];
        filter.commit();
        if frame >= unsupported_prefix {
            output.push(value);
            peak_abs_pressure_pa = peak_abs_pressure_pa.max(value.abs());
            identity.update(&value.to_bits().to_le_bytes());
        }
    }
    cx.checkpoint()
        .map_err(|_| PhysicalListeningMasterError::Artifact(AudioArtifactError::Cancelled))?;
    signal.start_time_s = start_time_s;
    signal.sample_rate_hz = SOUND_MASTER_SAMPLE_RATE_HZ;
    signal.pressure_pa = output;
    signal.peak_abs_pressure_pa = peak_abs_pressure_pa;
    signal.identity = identity.finalize();
    Ok(signal)
}


fn validate_pair(
    left: &PhysicalPressureSignal,
    right: &PhysicalPressureSignal,
) -> Result<(), PhysicalListeningMasterError> {
    for (field, matches) in [
        (
            "start_time_s",
            left.start_time_s.to_bits() == right.start_time_s.to_bits(),
        ),
        (
            "sample_rate_hz",
            left.sample_rate_hz == right.sample_rate_hz,
        ),
        (
            "sample_count",
            left.pressure_pa.len() == right.pressure_pa.len(),
        ),
        (
            "structural_basis_identity",
            left.structural_basis_identity == right.structural_basis_identity,
        ),
        (
            "damping_model_identity",
            left.damping_model_identity == right.damping_model_identity,
        ),
    ] {
        if !matches {
            return Err(PhysicalListeningMasterError::SignalMismatch { field });
        }
    }
    if left.pressure_pa.is_empty() {
        return Err(PhysicalListeningMasterError::SignalMismatch {
            field: "nonempty pressure signal",
        });
    }
    if left.sample_rate_hz != fs_evidence::cinematic_sound::SOUND_MASTER_SAMPLE_RATE_HZ {
        return Err(PhysicalListeningMasterError::SignalMismatch {
            field: "48 kHz cinematic sample rate",
        });
    }
    Ok(())
}

fn convert_pair(
    left: &PhysicalPressureSignal,
    right: &PhysicalPressureSignal,
    policy: PressureListeningMasterPolicy,
    gain_fs_per_pa: f64,
) -> Result<Vec<StereoSample>, PhysicalListeningMasterError> {
    let mut samples = Vec::new();
    samples
        .try_reserve_exact(left.pressure_pa.len())
        .map_err(|_| {
            PhysicalListeningMasterError::Artifact(AudioArtifactError::Capacity {
                artifact: "physical-pressure listening samples",
                requested: left.pressure_pa.len() as u64,
            })
        })?;
    for (frame, (&left_pa, &right_pa)) in
        left.pressure_pa.iter().zip(&right.pressure_pa).enumerate()
    {
        for (channel, pressure) in [("left", left_pa), ("right", right_pa)] {
            if !pressure.is_finite() {
                return Err(PhysicalListeningMasterError::NonFinitePressure { channel, frame });
            }
        }
        let window = boundary_window(frame, left.pressure_pa.len(), policy);
        let sample = if window.to_bits() == 0.0_f64.to_bits() {
            StereoSample {
                left_fs: 0.0,
                right_fs: 0.0,
            }
        } else {
            StereoSample {
                left_fs: left_pa * gain_fs_per_pa * window,
                right_fs: right_pa * gain_fs_per_pa * window,
            }
        };
        if !(sample.left_fs.is_finite() && sample.right_fs.is_finite()) {
            return Err(PhysicalListeningMasterError::GainOutsidePolicy {
                required_db: f64::INFINITY,
                maximum_absolute_db: 0.0,
            });
        }
        samples.push(sample);
    }
    Ok(samples)
}

fn validate_boundary_fades(
    policy: PressureListeningMasterPolicy,
    sample_frames: usize,
) -> Result<(), PhysicalListeningMasterError> {
    let initial = usize::try_from(policy.initial_fade_sample_frames)
        .map_err(|_| PhysicalListeningMasterError::InvalidPolicy)?;
    let terminal = usize::try_from(policy.terminal_fade_sample_frames)
        .map_err(|_| PhysicalListeningMasterError::InvalidPolicy)?;
    if initial < 2
        || terminal < 2
        || initial > sample_frames
        || terminal > sample_frames
        || initial
            .checked_add(terminal)
            .is_none_or(|sum| sum > sample_frames)
    {
        return Err(PhysicalListeningMasterError::InvalidPolicy);
    }
    Ok(())
}

fn boundary_window(
    frame: usize,
    sample_frames: usize,
    policy: PressureListeningMasterPolicy,
) -> f64 {
    let initial = policy.initial_fade_sample_frames as usize;
    if frame < initial {
        if frame == 0 {
            return 0.0;
        }
        if frame + 1 == initial {
            return 1.0;
        }
        let phase = core::f64::consts::PI * frame as f64 / (initial - 1) as f64;
        return 0.5 * (1.0 - det::cos(phase));
    }
    let terminal = policy.terminal_fade_sample_frames as usize;
    let terminal_start = sample_frames - terminal;
    if frame >= terminal_start {
        if frame == terminal_start {
            return 1.0;
        }
        if frame + 1 == sample_frames {
            return 0.0;
        }
        let phase = core::f64::consts::PI * (frame - terminal_start) as f64 / (terminal - 1) as f64;
        return 0.5 * (1.0 + det::cos(phase));
    }
    1.0
}

#[allow(clippy::too_many_arguments)]
fn listening_master_identity(
    left: &PhysicalPressureSignal,
    right: &PhysicalPressureSignal,
    policy: PressureListeningMasterPolicy,
    gain_fs_per_pa: f64,
    source_peak_abs_pressure_pa: f64,
    wav: WavCodecReceipt,
    samples: &[StereoSample],
) -> ContentHash {
    let mut hasher = DomainHasher::new(PHYSICAL_LISTENING_MASTER_IDENTITY_DOMAIN);
    hasher.update(left.identity.as_bytes());
    hasher.update(right.identity.as_bytes());
    for value in [
        policy.digital_gain_fs_per_pa,
        policy.maximum_true_peak_fs,
        gain_fs_per_pa,
        source_peak_abs_pressure_pa,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.update(&policy.initial_fade_sample_frames.to_le_bytes());
    hasher.update(&policy.terminal_fade_sample_frames.to_le_bytes());
    hasher.update(wav.wav_identity().as_bytes());
    hasher.update(&(samples.len() as u64).to_le_bytes());
    for sample in samples {
        hasher.update(&sample.left_fs.to_bits().to_le_bytes());
        hasher.update(&sample.right_fs.to_bits().to_le_bytes());
    }
    hasher.finalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x5052_4553_5357_4156,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn signal(scale: f64, radiation: u8) -> PhysicalPressureSignal {
        let pressure_pa = (0..48_000)
            .map(|sample| {
                let phase = 2.0 * core::f64::consts::PI * 440.0 * sample as f64 / 48_000.0;
                scale * det::sin(phase)
            })
            .collect::<Vec<_>>();
        let peak_abs_pressure_pa = pressure_pa
            .iter()
            .fold(0.0_f64, |peak, value| peak.max(value.abs()));
        PhysicalPressureSignal {
            start_time_s: 0.0,
            sample_rate_hz: 48_000,
            pressure_pa,
            peak_abs_pressure_pa,
            contact_force_sampling: crate::structural_acoustics::PhysicalContactForceSampling::IntervalMeanAtClosingElseOpeningEndpointZohV1,
            observer: crate::structural_acoustics::PhysicalPressureObserver::WorldFixed(
                crate::structural_acoustics::AcousticWorldObserver {
                    position_world_m: [f64::from(radiation), 0.0, 1.0],
                },
            ),
            structural_basis_identity: ContentHash([1; 32]),
            radiation_identity: ContentHash([radiation; 32]),
            damping_model_identity: ContentHash([2; 32]),
            identity: ContentHash([radiation.wrapping_add(10); 32]),
        }
    }

    #[test]
    fn moving_pressure_antialias_and_pcm24_keep_physical_time_and_receiver_identity() {
        use crate::structural_acoustics::PhysicalPressureObserver;

        for opening_index in [0, 1] {
            let mut input = signal(1.0e-4, 3);
            input.start_time_s = f64::from(opening_index) / 96_000.0;
            input.sample_rate_hz = 96_000;
            input.pressure_pa = (0..4_097).map(|frame| {
                let time = input.start_time_s + (frame + 1) as f64 / 96_000.0;
                1.0e-4 * (
                    0.8 * det::sin(core::f64::consts::TAU * 2_000.0 * time)
                    + 0.5 * det::sin(core::f64::consts::TAU * 30_000.0 * time)
                )
            }).collect();
            input.peak_abs_pressure_pa = input.pressure_pa.iter()
                .fold(0.0_f64, |peak, value| peak.max(value.abs()));
            input.observer = PhysicalPressureObserver::WorldMoving {
                trajectory_identity: ContentHash([17; 32]),
            };
            let observer = input.observer;
            let pressure = with_cx(|cx| decimate_retarded_pressure_to_master(input, cx)).unwrap();
            assert_eq!(pressure.observer, observer);
            let opening = pressure.start_time_s * 48_000.0;
            assert!((opening - opening.round()).abs() < 1.0e-10);
            for (frame, value) in pressure.pressure_pa.iter().enumerate() {
                let time = pressure.start_time_s + (frame + 1) as f64 / 48_000.0;
                let expected = 0.8e-4 * det::sin(core::f64::consts::TAU * 2_000.0 * time);
                assert!((value - expected).abs() < 2.0e-9);
            }
            let master = with_cx(|cx| PhysicalPressureListeningMaster::try_build(
                &pressure, &pressure, PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::default(), AudioArtifactBudget::DEFAULT, cx,
            )).unwrap();
            assert_eq!(master.left_pressure_identity, pressure.identity);
            assert_eq!(master.right_pressure_identity, pressure.identity);
            assert_eq!(master.wav.sample_frame_count(), pressure.pressure_pa.len() as u64);
            assert_eq!(master.wav.encoding(), WavSampleEncoding::Pcm24);
            assert!(master.meters.sample_peak_fs > 0.03);
            assert!(master.meters.sample_peak_fs < 0.05);
        }
    }

    #[test]
    fn g0_listening_master_uses_fixed_pressure_calibration_without_normalization() {
        let left = signal(2.0e-5, 3);
        let right = signal(1.0e-5, 4);
        let master = with_cx(|cx| {
            PhysicalPressureListeningMaster::try_build(
                &left,
                &right,
                PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::try_new(Some(
                    "physical Pa observers; fixed-calibration listening derivative".to_owned(),
                ))
                .unwrap(),
                AudioArtifactBudget::DEFAULT,
                cx,
            )
        })
        .unwrap();
        assert_eq!(master.left_pressure_identity, left.identity);
        assert_eq!(master.right_pressure_identity, right.identity);
        assert_eq!(master.wav.sample_frame_count(), 48_000);
        assert_eq!(master.digital_gain_fs_per_pa.to_bits(), 512.0_f64.to_bits());
        assert!(master.digital_gain_db > 54.0 && master.digital_gain_db < 54.3);
        assert!(master.meters.sample_peak_fs > 4.0e-3);
        assert!(master.meters.sample_peak_fs < 16.0e-3);
        assert!(
            master.meters.true_peak_estimate_fs
                <= PressureListeningMasterPolicy::CRITIQUE.maximum_true_peak_fs + 1.0e-12
        );
    }

    #[test]
    fn g0_fixed_pressure_calibration_refuses_instead_of_normalizing_overload() {
        let left = signal(1.0, 3);
        let right = signal(0.5, 4);
        let error = with_cx(|cx| {
            PhysicalPressureListeningMaster::try_build(
                &left,
                &right,
                PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            )
        })
        .unwrap_err();
        assert!(matches!(
            error,
            PhysicalListeningMasterError::PeakTargetExceeded { target_fs, .. }
                if target_fs.to_bits()
                    == PressureListeningMasterPolicy::CRITIQUE
                        .maximum_true_peak_fs
                        .to_bits()
        ));
    }

    #[test]
    fn g0_listening_boundary_window_removes_cut_click_without_mutating_pressure() {
        let mut left = signal(2.0e-5, 3);
        let mut right = signal(1.0e-5, 4);
        left.pressure_pa[0] = 2.0e-5;
        right.pressure_pa[0] = -1.0e-5;
        let left_before = left.pressure_pa.clone();
        let right_before = right.pressure_pa.clone();
        assert_ne!(left_before[0].to_bits(), 0.0_f64.to_bits());
        assert_ne!(right_before[0].to_bits(), 0.0_f64.to_bits());
        let samples =
            convert_pair(&left, &right, PressureListeningMasterPolicy::CRITIQUE, 1.0).unwrap();

        assert_eq!(
            samples.first().unwrap().left_fs.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            samples.first().unwrap().right_fs.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(samples.last().unwrap().left_fs.to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            samples.last().unwrap().right_fs.to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            boundary_window(
                PressureListeningMasterPolicy::CRITIQUE.initial_fade_sample_frames as usize - 1,
                samples.len(),
                PressureListeningMasterPolicy::CRITIQUE,
            )
            .to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(left.pressure_pa, left_before);
        assert_eq!(right.pressure_pa, right_before);
    }

    #[test]
    fn g0_silence_and_cross_basis_pairing_refuse() {
        let mut left = signal(0.0, 3);
        let right = signal(0.0, 4);
        assert!(matches!(
            with_cx(|cx| PhysicalPressureListeningMaster::try_build(
                &left,
                &right,
                PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            )),
            Err(PhysicalListeningMasterError::SilentPhysicalSignal)
        ));
        left.structural_basis_identity = ContentHash([9; 32]);
        let right = signal(1.0, 4);
        assert!(matches!(
            with_cx(|cx| PhysicalPressureListeningMaster::try_build(
                &left,
                &right,
                PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::default(),
                AudioArtifactBudget::DEFAULT,
                cx,
            )),
            Err(PhysicalListeningMasterError::SignalMismatch {
                field: "structural_basis_identity"
            })
        ));
    }

    #[test]
    fn g0_physical_pressure_crop_is_exact_and_rebased_without_processing() {
        let source = signal(2.0e-5, 3);
        let cropped = source.try_crop_rebased(2_000, 48_000, 0.0).unwrap();
        assert_eq!(cropped.start_time_s, 0.0);
        assert_eq!(cropped.pressure_pa, source.pressure_pa[2_000..48_000]);
        assert_eq!(cropped.sample_rate_hz, source.sample_rate_hz);
        assert_eq!(
            cropped.structural_basis_identity,
            source.structural_basis_identity
        );
        assert_ne!(cropped.identity, source.identity);
        assert!(source.try_crop_rebased(2_000, 2_000, 0.0).is_err());
        assert!(source.try_crop_rebased(0, 48_001, 0.0).is_err());
    }
}
