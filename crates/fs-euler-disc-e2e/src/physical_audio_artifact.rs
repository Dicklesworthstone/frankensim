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
use crate::structural_acoustics::PhysicalPressureSignal;

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
