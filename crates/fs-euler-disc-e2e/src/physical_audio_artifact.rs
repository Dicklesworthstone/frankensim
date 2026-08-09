//! Audible delivery of physical observer-pressure signals.
//!
//! [`PhysicalPressureSignal`](crate::structural_acoustics::PhysicalPressureSignal)
//! is the scientific product: its samples are pressure in pascals at a named
//! observer. A media player instead expects unitless digital-full-scale
//! samples. This module keeps that conversion explicit and replayable. It
//! never changes, normalizes, or relabels the source pressure signal; it emits
//! a separate listening master whose gain is presentation metadata.

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
    "org.frankensim.euler-disc.physical-listening-master.v1";

/// Explicit policy for the nonphysical pressure-to-digital presentation step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PressureListeningMasterPolicy {
    /// Desired four-times-oversampled peak in digital-full-scale units.
    pub target_true_peak_fs: f64,
    /// Largest admitted absolute pressure-to-digital gain adjustment [dB].
    pub maximum_absolute_gain_db: f64,
}

impl PressureListeningMasterPolicy {
    /// Conservative critique master: about -6.94 dBFS estimated true peak.
    pub const CRITIQUE: Self = Self {
        target_true_peak_fs: 0.45,
        maximum_absolute_gain_db: 180.0,
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
    /// signals must share the same structural basis, damping model, start
    /// time, sample rate, and frame count; their radiation identities may
    /// differ because they represent different observer locations.
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
        if !(policy.target_true_peak_fs > 0.0
            && policy.target_true_peak_fs < 1.0
            && policy.target_true_peak_fs.is_finite()
            && policy.maximum_absolute_gain_db > 0.0
            && policy.maximum_absolute_gain_db.is_finite())
        {
            return Err(PhysicalListeningMasterError::InvalidPolicy);
        }

        let source_peak_abs_pressure_pa = left.peak_abs_pressure_pa.max(right.peak_abs_pressure_pa);
        if !(source_peak_abs_pressure_pa > 0.0 && source_peak_abs_pressure_pa.is_finite()) {
            return Err(PhysicalListeningMasterError::SilentPhysicalSignal);
        }

        // First place the stored-sample peak at the requested ceiling, then
        // account for intersample overshoot using the same declared estimator
        // published by the WAV artifact layer.
        let mut digital_gain_fs_per_pa = policy.target_true_peak_fs / source_peak_abs_pressure_pa;
        let mut samples = convert_pair(left, right, digital_gain_fs_per_pa)?;
        let first_meters = measure_audio(&samples, budget, cx)?;
        let first_peak = first_meters
            .sample_peak_fs
            .max(first_meters.true_peak_estimate_fs);
        if !(first_peak > 0.0 && first_peak.is_finite()) {
            return Err(PhysicalListeningMasterError::SilentPhysicalSignal);
        }
        digital_gain_fs_per_pa *= policy.target_true_peak_fs / first_peak;
        let digital_gain_db = 20.0 * det::ln(digital_gain_fs_per_pa) / core::f64::consts::LN_10;
        if !digital_gain_db.is_finite() || digital_gain_db.abs() > policy.maximum_absolute_gain_db {
            return Err(PhysicalListeningMasterError::GainOutsidePolicy {
                required_db: digital_gain_db,
                maximum_absolute_db: policy.maximum_absolute_gain_db,
            });
        }

        samples = convert_pair(left, right, digital_gain_fs_per_pa)?;
        let meters = measure_audio(&samples, budget, cx)?;
        let observed_peak = meters.sample_peak_fs.max(meters.true_peak_estimate_fs);
        let tolerance = 256.0 * f64::EPSILON * policy.target_true_peak_fs;
        if observed_peak > policy.target_true_peak_fs + tolerance {
            return Err(PhysicalListeningMasterError::PeakTargetExceeded {
                observed_fs: observed_peak,
                target_fs: policy.target_true_peak_fs,
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
    /// The declared true-peak ceiling was not met after deterministic scaling.
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
        let sample = StereoSample {
            left_fs: left_pa * gain_fs_per_pa,
            right_fs: right_pa * gain_fs_per_pa,
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
        policy.target_true_peak_fs,
        policy.maximum_absolute_gain_db,
        gain_fs_per_pa,
        source_peak_abs_pressure_pa,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
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
            structural_basis_identity: ContentHash([1; 32]),
            radiation_identity: ContentHash([radiation; 32]),
            damping_model_identity: ContentHash([2; 32]),
            identity: ContentHash([radiation.wrapping_add(10); 32]),
        }
    }

    #[test]
    fn g0_physical_pressure_remains_bound_while_listening_master_is_audible() {
        let left = signal(2.0e-5, 3);
        let right = signal(1.0e-5, 4);
        let master = with_cx(|cx| {
            PhysicalPressureListeningMaster::try_build(
                &left,
                &right,
                PressureListeningMasterPolicy::CRITIQUE,
                &WavMetadata::try_new(Some(
                    "physical Pa observers; presentation-normalized listening derivative"
                        .to_owned(),
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
        assert!(master.digital_gain_db > 80.0);
        assert!(master.meters.sample_peak_fs > 0.4);
        assert!(master.meters.true_peak_estimate_fs <= 0.45 + 1.0e-12);
        assert!(master.meters.integrated_loudness_lufs.is_some());
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
}
