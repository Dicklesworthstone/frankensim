//! Shared physical trajectories without mixing synchronized sensor channels.

use super::{
    ChirpedToneFit, DampedChirpedToneFit, DampedTonePoleFit, HarmonicFit, PartialFitFailure,
    PartialFrame, PartialFrameEstimate, PartialFrameSchedule, PartialTrackingError,
    PartialTrackingModel, PartialTrackingSummary, ToneFrequencyFit, TrackedPartial, track_channels,
};
use core::ops::ControlFlow;
use fs_fft::harmonic::{
    ChirpedToneChannelsFit, DampedChirpedToneChannelsFit, DampedTonePoleChannelsFit,
    ToneFrequencyChannelsFit,
};

/// A shared physical partial law and separate original-channel observations.
///
/// No channel is averaged, normalized, reordered or selected as a reference.
/// A mode absent at one sensor can still be identified through other sensors.
/// Numerical conditioning and convergence refer to the whole joint fit.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialChannelsFrameEstimate {
    /// Shared independent frequencies and per-channel peak phasors.
    Stationary(ToneFrequencyChannelsFit),
    /// Shared oscillation/decay poles and per-channel first-sample phasors.
    Ringdown(DampedTonePoleChannelsFit),
    /// Shared first/last frequency endpoints and per-channel peak phasors.
    Chirped(ChirpedToneChannelsFit),
    /// Shared frequency endpoints and decay rates, with separate initial phasors.
    DampedChirped(DampedChirpedToneChannelsFit),
}

impl PartialChannelsFrameEstimate {
    /// Original-channel DC, peak phasors and unexplained sampled RMS.
    /// These retain the input units and each frame's first-sample phase origin.
    #[must_use]
    pub fn fits(&self) -> &[HarmonicFit] {
        match self {
            Self::Stationary(result) => &result.fits,
            Self::Ringdown(result) => &result.fits,
            Self::Chirped(result) => &result.fits,
            Self::DampedChirped(result) => &result.fits,
        }
    }

    /// Observe one original partial in one original channel. Shared frequencies
    /// and decay laws accompany that channel's own complex peak coefficient.
    /// Invalid channel or partial indices return `None` without reassignment.
    #[must_use]
    pub fn partial(&self, channel: usize, index: usize) -> Option<TrackedPartial> {
        let peak_phasor = *self.fits().get(channel)?.phasors.get(index)?;
        let (frequency_hz, end_frequency_hz, decay_rate_per_s) = match self {
            Self::Stationary(result) => (*result.frequencies_hz.get(index)?, None, None),
            Self::Ringdown(result) => {
                let pole = result.poles.get(index)?;
                (pole.frequency_hz, None, Some(pole.decay_rate_per_s))
            }
            Self::Chirped(result) => {
                let tone = result.tones.get(index)?;
                (tone.start_frequency_hz, Some(tone.end_frequency_hz), None)
            }
            Self::DampedChirped(result) => {
                let tone = result.tones.get(index)?;
                (
                    tone.start_frequency_hz,
                    Some(tone.end_frequency_hz),
                    Some(tone.decay_rate_per_s),
                )
            }
        };
        Some(TrackedPartial {
            frequency_hz,
            end_frequency_hz,
            decay_rate_per_s,
            peak_phasor,
        })
    }

    // Only called for the scalar facade, whose admitted channel count is one.
    // Move every retained vector; scalar compatibility adds no result clones.
    fn into_scalar(self) -> PartialFrameEstimate {
        match self {
            Self::Stationary(mut result) => PartialFrameEstimate::Stationary(ToneFrequencyFit {
                frequencies_hz: result.frequencies_hz,
                fit: result.fits.remove(0),
                evaluations: result.evaluations,
                last_corrections_hz: result.last_corrections_hz,
                frequency_reciprocal_condition: result.frequency_reciprocal_condition,
            }),
            Self::Ringdown(mut result) => PartialFrameEstimate::Ringdown(DampedTonePoleFit {
                poles: result.poles,
                fit: result.fits.remove(0),
                evaluations: result.evaluations,
                last_corrections: result.last_corrections,
                pole_reciprocal_condition: result.pole_reciprocal_condition,
            }),
            Self::Chirped(mut result) => PartialFrameEstimate::Chirped(ChirpedToneFit {
                tones: result.tones,
                fit: result.fits.remove(0),
                evaluations: result.evaluations,
                last_corrections: result.last_corrections,
                chirp_reciprocal_condition: result.chirp_reciprocal_condition,
            }),
            Self::DampedChirped(mut result) => {
                PartialFrameEstimate::DampedChirped(DampedChirpedToneFit {
                    tones: result.tones,
                    fit: result.fits.remove(0),
                    evaluations: result.evaluations,
                    last_corrections: result.last_corrections,
                    trajectory_reciprocal_condition: result.trajectory_reciprocal_condition,
                })
            }
        }
    }
}

/// One synchronized frame, retaining the entire joint-fit outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialChannelsFrame {
    /// Exact first sample on every channel's shared source clock.
    pub first_sample_index: u64,
    /// Number of original samples in EACH channel of this frame.
    pub sample_count: usize,
    /// Uniform physical sampling interval [s], common to all channels.
    pub timestep_s: f64,
    /// A failed shared fit publishes no partial channel spectra or zero fills.
    pub estimate: Result<PartialChannelsFrameEstimate, PartialFitFailure>,
}

impl PartialChannelsFrame {
    pub(super) fn into_scalar(self) -> PartialFrame {
        PartialFrame {
            first_sample_index: self.first_sample_index,
            sample_count: self.sample_count,
            timestep_s: self.timestep_s,
            estimate: self.estimate.map(PartialChannelsFrameEstimate::into_scalar),
        }
    }
}

/// Track synchronized, equal-length COMMON-UNIT records with shared partial laws.
///
/// All channels inform the same frequencies, decay rates and/or chirp endpoints;
/// each keeps its own DC, first-sample peak phasors and residual. The objective
/// is the sum of original-channel squared residuals. No averaging or individual
/// normalization can erase opposite-sign observations or alter physical vector
/// weights. Orthogonal rotations of common-unit vector components preserve this
/// objective; different units or sensor-noise weights require explicit caller
/// conversion. No delay alignment, resampling or microphone calibration occurs.
/// In particular, independently delayed chirps need not share one frequency law.
///
/// This uses [`super::track_partials`]' complete-record schedule, immutable
/// domains, successful-anchor warm starts, exact integer source clock and typed
/// unresolved-frame semantics. A channel at a modal node is permitted when the
/// full array identifies the shared trajectory; failure concerns the WHOLE
/// physical fit. No missing channel is filled or silently dropped. All channel
/// lengths and samples are admitted before the first observation callback.
///
/// Only one shared amplitude QR is prepared per candidate. With C channels,
/// N samples per frame, K partials and at most E candidates per frame, work is
/// O(E*C*(N*K² + K³)); numeric workspace is O(N*K + K² + C*K), plus the charged
/// warm seeds and C reusable slice views. One channel uses a stack view and
/// preserves the scalar owner's original storage budget. The callback borrows
/// one result at a time; caller-retained clones are outside that budget.
/// Cancellation is observed between admitted frame fits, after input admission,
/// and stops future frames when the callback returns `ControlFlow::Break(())`.
pub fn track_partial_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    schedule: PartialFrameSchedule,
    model: PartialTrackingModel<'_>,
    mut visit: impl FnMut(&PartialChannelsFrame) -> ControlFlow<()>,
) -> Result<PartialTrackingSummary, PartialTrackingError> {
    track_channels(channels, timestep_s, schedule, model, |frame| visit(&frame))
}
