//! Bounded physical partial trajectories over a complete sampled record.
//!
//! Every frame jointly fits DC and independently spaced partials through the
//! existing `fs_fft::harmonic` owners. Stationary frames refine frequencies;
//! ringdown frames refine both damped frequencies and amplitude decay rates.
//! Linear-chirp frames infer first- and last-sample frequencies, retaining
//! within-frame drift instead of charging it to a stationary fit's residual.
//! This avoids forcing stiff-string partials onto integer harmonic bins or
//! fitting one decay slope to an entire note. Equal-frequency, different-decay
//! poles remain separate when the supplied decay bounds distinguish them.
//!
//! Fixed disjoint search domains define partial identities. Successful fits
//! warm-start the next frame; unresolved frames retain their typed refusal and
//! do not change those seeds. No automatic discovery, reassignment, silence
//! substitution, uncertainty, noise-floor detection or instrument acceptance
//! is implied. Nonlinear drift, attacks, growth and other omitted dynamics can
//! bias a fit and remain in its residual; a frame-length ladder is needed when
//! the selected local model is only an approximation.

use core::{mem::size_of, ops::ControlFlow};
use fs_fft::C64;
use fs_fft::harmonic::{
    ChirpedTone, ChirpedToneBounds, ChirpedToneFit, ChirpedToneRefinementBudget,
    ChirpedToneRefinementError, DampedToneBounds, DampedTonePole, DampedTonePoleFit,
    DampedToneRefinementBudget, DampedToneRefinementError, HarmonicFit, HarmonicFitError,
    ToneFrequencyBudget, ToneFrequencyError, ToneFrequencyFit, refine_chirped_tones,
    refine_damped_tone_poles, refine_tone_frequencies,
};

/// Complete-record frame schedule and its explicit work cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PartialFrameSchedule {
    /// Number of uniformly sampled values in each fit.
    pub frame_samples: usize,
    /// Positive advance, at most `frame_samples`, so coverage has no gaps.
    pub hop_samples: usize,
    /// Refuse the entire analysis before fitting if more frames are required.
    pub max_frames: usize,
    /// Exact index of the record's first sample on the caller's source clock.
    /// Large epochs are retained as integers, never rounded to f64 seconds.
    pub first_sample_index: u64,
}

/// Caller-authored local model, partial identities and numerical budgets.
///
/// The fit's workspace cap also pays for the tracker's mutable seed array.
/// Chirps retain both a successful anchor and a reusable prediction array.
/// Input records, borrowed declarations and memory retained by the callback
/// are caller-owned and excluded, as is allocator metadata.
#[derive(Debug, Clone, Copy)]
pub enum PartialTrackingModel<'a> {
    /// Independent stationary tones; decay is not estimated or assumed zero.
    /// Frequency intervals must be ordered, positive, disjoint and sub-Nyquist.
    Stationary {
        /// Initial frequencies in the immutable partial-identity order.
        initial_frequencies_hz: &'a [f64],
        /// Fixed local search interval for each partial, in hertz.
        frequency_bounds_hz: &'a [[f64; 2]],
        /// Per-frame refinement, conditioning and total workspace limits.
        budget: ToneFrequencyBudget,
    },
    /// Independently decaying oscillatory partials. Disjoint decay intervals
    /// may distinguish poles whose frequency intervals overlap.
    Ringdown {
        /// Initial damped frequencies and amplitude decay rates.
        initial_poles: &'a [DampedTonePole],
        /// Fixed local search rectangle for each partial.
        bounds: &'a [DampedToneBounds],
        /// Per-frame refinement, conditioning and total workspace limits.
        budget: DampedToneRefinementBudget,
    },
    /// Independent linear frequency trajectories with constant per-frame
    /// amplitudes. Frequencies refer to the first and last sampled values;
    /// their separation in time is `(frame_samples - 1) * timestep_s`.
    Chirped {
        /// Initial endpoint frequencies for the first scheduled frame.
        initial_tones: &'a [ChirpedTone],
        /// Fixed endpoint rectangles, disjoint in at least one coordinate.
        /// They must admit every frame in the intended trajectory.
        bounds: &'a [ChirpedToneBounds],
        /// Per-frame refinement, conditioning and total workspace limits.
        budget: ChirpedToneRefinementBudget,
    },
}

/// One estimated partial, in the original declaration order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrackedPartial {
    /// Measured oscillation frequency [Hz], not an integer harmonic label.
    /// For a chirp this is the frequency at THIS frame's first sample.
    pub frequency_hz: f64,
    /// Last-sample instantaneous frequency [Hz] for a chirp. Other models
    /// return `None` because they fit one frequency throughout the frame.
    pub end_frequency_hz: Option<f64>,
    /// Measured amplitude decay [1/s] for a ringdown model. Other models
    /// return `None`: they do not measure decay. Zero is an estimated zero.
    pub decay_rate_per_s: Option<f64>,
    /// Peak `exp(-i phase(t))` coefficient in input units, with phase zero at
    /// THIS frame's first sample. A decaying coefficient is its initial amplitude,
    /// not its record-averaged amplitude. Frames have different phase origins.
    pub peak_phasor: C64,
}

/// Successful local fit, retaining the numerical owner's full diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialFrameEstimate {
    /// Joint stationary-tone identification.
    Stationary(ToneFrequencyFit),
    /// Joint oscillatory-decay identification.
    Ringdown(DampedTonePoleFit),
    /// Joint first/last-sample frequency-law identification.
    Chirped(ChirpedToneFit),
}

impl PartialFrameEstimate {
    /// Joint DC, peak phasors and unexplained sampled RMS in input units.
    #[must_use]
    pub fn fit(&self) -> &HarmonicFit {
        match self {
            Self::Stationary(result) => &result.fit,
            Self::Ringdown(result) => &result.fit,
            Self::Chirped(result) => &result.fit,
        }
    }

    /// Observe an original partial identity without sorting by fitted frequency.
    #[must_use]
    pub fn partial(&self, index: usize) -> Option<TrackedPartial> {
        let peak_phasor = *self.fit().phasors.get(index)?;
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
        };
        Some(TrackedPartial {
            frequency_hz,
            end_frequency_hz,
            decay_rate_per_s,
            peak_phasor,
        })
    }
}

/// The exact numerical owner's refusal, without a fabricated partial spectrum.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialFitFailure {
    /// Independent-tone refusal.
    Stationary(ToneFrequencyError),
    /// Damped-pole refusal.
    Ringdown(DampedToneRefinementError),
    /// Linear-chirp refusal.
    Chirped(ChirpedToneRefinementError),
}

impl PartialFitFailure {
    fn unresolved(&self) -> bool {
        matches!(
            self,
            Self::Stationary(
                ToneFrequencyError::Unidentifiable { .. }
                    | ToneFrequencyError::NotConverged { .. }
                    | ToneFrequencyError::Fit(HarmonicFitError::IllConditioned { .. })
            ) | Self::Ringdown(
                DampedToneRefinementError::Unidentifiable { .. }
                    | DampedToneRefinementError::NotConverged { .. }
                    | DampedToneRefinementError::Fit(HarmonicFitError::IllConditioned { .. })
            ) | Self::Chirped(
                ChirpedToneRefinementError::Unidentifiable { .. }
                    | ChirpedToneRefinementError::NotConverged { .. }
                    | ChirpedToneRefinementError::Fit(HarmonicFitError::IllConditioned { .. })
            )
        )
    }
}

impl core::fmt::Display for PartialFitFailure {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Stationary(error) => error.fmt(f),
            Self::Ringdown(error) => error.fmt(f),
            Self::Chirped(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for PartialFitFailure {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Stationary(error) => Some(error),
            Self::Ringdown(error) => Some(error),
            Self::Chirped(error) => Some(error),
        }
    }
}

/// One complete analysis window, including an explicit unresolved outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct PartialFrame {
    /// Exact first sample on the supplied source clock; local phasors refer here.
    pub first_sample_index: u64,
    /// Number of original samples in this frame.
    pub sample_count: usize,
    /// Uniform physical sampling interval [s].
    pub timestep_s: f64,
    /// No coefficients are published when local identification is unresolved.
    pub estimate: Result<PartialFrameEstimate, PartialFitFailure>,
}

/// Work admission or a fatal failure of the requested observation.
#[derive(Debug, Clone, PartialEq)]
pub enum PartialTrackingError {
    /// Invalid clock, coverage, record, index arithmetic or work extent.
    InvalidParameter {
        /// The refused condition.
        what: &'static str,
    },
    /// Complete coverage exceeds the caller's frame cap; no callback has run.
    FrameBudgetExceeded {
        /// Frames required, including a terminal overlapping frame if needed.
        required_frames: usize,
        /// Caller-supplied cap.
        max_frames: usize,
    },
    /// Invalid fit declarations, insufficient storage or nonfinite arithmetic.
    /// Previous callbacks, if any, remain observations already delivered.
    Fit {
        /// Exact starting sample of the refused frame.
        first_sample_index: u64,
        /// Original numerical refusal. Its workspace cap has the tracker seed
        /// storage subtracted, as documented by [`PartialTrackingModel`].
        failure: PartialFitFailure,
    },
}

impl core::fmt::Display for PartialTrackingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidParameter { what } => write!(f, "partial tracking: {what}"),
            Self::FrameBudgetExceeded {
                required_frames,
                max_frames,
            } => write!(
                f,
                "partial tracking needs {required_frames} frames, cap is {max_frames}"
            ),
            Self::Fit {
                first_sample_index,
                failure,
            } => write!(
                f,
                "partial tracking at sample {first_sample_index}: {failure}"
            ),
        }
    }
}

impl std::error::Error for PartialTrackingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit { failure, .. } => Some(failure),
            _ => None,
        }
    }
}

/// Completion and the conservative admitted work bound, including failures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PartialTrackingSummary {
    /// Frames needed for complete record coverage.
    pub planned_frames: usize,
    /// Frames delivered to the callback, whether resolved or unresolved.
    pub visited_frames: usize,
    /// Delivered frames with a converged physical fit.
    pub resolved_frames: usize,
    /// `planned_frames * per_frame_max_evaluations`, not a measured count.
    /// Rejected line-search candidates are included by the shared owner.
    pub maximum_candidate_evaluations: usize,
    /// False only when the callback stopped before the final frame.
    pub completed: bool,
}

enum Seeds {
    Stationary(Vec<f64>),
    Ringdown(Vec<DampedTonePole>),
    Chirped {
        anchor: Vec<ChirpedTone>,
        prediction: Vec<ChirpedTone>,
        anchor_offset: usize,
    },
}

fn invalid(what: &'static str) -> PartialTrackingError {
    PartialTrackingError::InvalidParameter { what }
}

/// Track a complete record with bounded per-frame fits and streaming results.
///
/// Regular windows start at 0, hop, 2*hop, ... . If the regular grid misses the
/// record's end, one final FULL window ends exactly at the last input sample;
/// no zero padding, dropped tail, default attack skip or minimum frequency
/// filter is applied. The last hop can therefore be shorter than requested.
///
/// Only the latest successful frequencies/poles are retained as seeds. All
/// search bounds and partial identities stay fixed, including across failures.
/// A chirp's two endpoints advance from its last successful frame by the ACTUAL
/// sample-offset gap, including shorter terminal hops and unresolved frames.
/// Predictions are only initial guesses and are clamped to the fixed domains;
/// fitted results must still satisfy the original untruncated-correction gates.
/// The original source epoch stays integer; no absolute timestamp enters the
/// predictor. Fixed rectangles do not provide arbitrary track reassignment.
/// Unidentifiable, ill-conditioned and unconverged frames are delivered as
/// `Err` observations; invalid declarations, storage refusals and nonfinite
/// arithmetic abort with [`PartialTrackingError`]. All input samples are checked
/// before the first callback. Stopping and restarting is not a continuation
/// protocol: a separate call starts from the caller's declarations again.
///
/// The callback receives one borrowed result at a time and may stop with
/// `ControlFlow::Break(())`. Cancellation is bounded by one admitted frame fit;
/// the numerical owners have no mid-factorization cancellation. Caller-retained
/// result clones are not charged to the working-memory budget. Work is at most
/// F*E candidate sets, with the existing joint QR cost per candidate; numeric
/// working storage is O(N*K + K*K) plus K warm seeds (two arrays for chirps),
/// independent of F. Chirped frames do not also identify amplitude decay.
pub fn track_partials(
    samples: &[f64],
    timestep_s: f64,
    schedule: PartialFrameSchedule,
    mut model: PartialTrackingModel<'_>,
    mut visit: impl FnMut(&PartialFrame) -> ControlFlow<()>,
) -> Result<PartialTrackingSummary, PartialTrackingError> {
    if !(timestep_s.is_finite() && timestep_s > 0.0)
        || schedule.frame_samples == 0
        || schedule.hop_samples == 0
        || schedule.hop_samples > schedule.frame_samples
        || samples.len() < schedule.frame_samples
    {
        return Err(invalid(
            "need a positive clock and full, gap-free analysis windows",
        ));
    }
    let last_offset = samples.len() - schedule.frame_samples;
    let planned_frames = (last_offset / schedule.hop_samples)
        .checked_add(1 + usize::from(last_offset % schedule.hop_samples != 0))
        .ok_or_else(|| invalid("frame count overflows"))?;
    if planned_frames > schedule.max_frames {
        return Err(PartialTrackingError::FrameBudgetExceeded {
            required_frames: planned_frames,
            max_frames: schedule.max_frames,
        });
    }
    let last_index = u64::try_from(samples.len() - 1)
        .ok()
        .and_then(|offset| schedule.first_sample_index.checked_add(offset));
    if last_index.is_none() || samples.iter().any(|sample| !sample.is_finite()) {
        return Err(invalid(
            "source sample indices must fit u64 and all samples must be finite",
        ));
    }
    let fit_error = |failure| PartialTrackingError::Fit {
        first_sample_index: schedule.first_sample_index,
        failure,
    };
    let (seeds, max_evaluations) = match &mut model {
        PartialTrackingModel::Stationary {
            initial_frequencies_hz,
            budget,
            ..
        } => {
            reserve_seeds(
                initial_frequencies_hz.len(),
                size_of::<f64>(),
                &mut budget.fit,
            )
            .map_err(|error| {
                fit_error(PartialFitFailure::Stationary(ToneFrequencyError::Fit(
                    error,
                )))
            })?;
            (
                Seeds::Stationary(initial_frequencies_hz.to_vec()),
                budget.max_evaluations,
            )
        }
        PartialTrackingModel::Ringdown {
            initial_poles,
            budget,
            ..
        } => {
            reserve_seeds(
                initial_poles.len(),
                size_of::<DampedTonePole>(),
                &mut budget.fit,
            )
            .map_err(|error| {
                fit_error(PartialFitFailure::Ringdown(DampedToneRefinementError::Fit(
                    error,
                )))
            })?;
            (
                Seeds::Ringdown(initial_poles.to_vec()),
                budget.max_evaluations,
            )
        }
        PartialTrackingModel::Chirped {
            initial_tones,
            budget,
            ..
        } => {
            reserve_seeds(
                initial_tones.len(),
                2 * size_of::<ChirpedTone>(),
                &mut budget.fit,
            )
            .map_err(|error| {
                fit_error(PartialFitFailure::Chirped(ChirpedToneRefinementError::Fit(
                    error,
                )))
            })?;
            (
                Seeds::Chirped {
                    anchor: initial_tones.to_vec(),
                    prediction: initial_tones.to_vec(),
                    anchor_offset: 0,
                },
                budget.max_evaluations,
            )
        }
    };
    let maximum_candidate_evaluations = planned_frames
        .checked_mul(max_evaluations)
        .ok_or_else(|| invalid("candidate evaluation count overflows"))?;
    let mut seeds = seeds;
    let mut summary = PartialTrackingSummary {
        planned_frames,
        visited_frames: 0,
        resolved_frames: 0,
        maximum_candidate_evaluations,
        completed: false,
    };
    let mut offset = 0;
    loop {
        let first_sample_index = schedule.first_sample_index + offset as u64;
        let block = &samples[offset..offset + schedule.frame_samples];
        let estimate = match (&model, &mut seeds) {
            (
                PartialTrackingModel::Stationary {
                    frequency_bounds_hz,
                    budget,
                    ..
                },
                Seeds::Stationary(seed),
            ) => refine_tone_frequencies(block, timestep_s, seed, frequency_bounds_hz, *budget)
                .map(|result| {
                    seed.copy_from_slice(&result.frequencies_hz);
                    PartialFrameEstimate::Stationary(result)
                })
                .map_err(PartialFitFailure::Stationary),
            (PartialTrackingModel::Ringdown { bounds, budget, .. }, Seeds::Ringdown(seed)) => {
                refine_damped_tone_poles(block, timestep_s, seed, bounds, *budget)
                    .map(|result| {
                        seed.copy_from_slice(&result.poles);
                        PartialFrameEstimate::Ringdown(result)
                    })
                    .map_err(PartialFitFailure::Ringdown)
            }
            (
                PartialTrackingModel::Chirped { bounds, budget, .. },
                Seeds::Chirped {
                    anchor,
                    prediction,
                    anchor_offset,
                },
            ) => {
                let gap = offset - *anchor_offset;
                if gap == 0 {
                    // Admit the original declarations through the numerical
                    // owner before any prediction may clamp search bounds.
                    prediction.copy_from_slice(anchor);
                } else {
                    let advance = gap as f64 / (schedule.frame_samples - 1) as f64;
                    for ((trial, accepted), domain) in
                        prediction.iter_mut().zip(anchor.iter()).zip(bounds.iter())
                    {
                        let change = accepted.end_frequency_hz - accepted.start_frequency_hz;
                        trial.start_frequency_hz = change
                            .mul_add(advance, accepted.start_frequency_hz)
                            .clamp(domain.start_frequency_hz[0], domain.start_frequency_hz[1]);
                        trial.end_frequency_hz = change
                            .mul_add(advance, accepted.end_frequency_hz)
                            .clamp(domain.end_frequency_hz[0], domain.end_frequency_hz[1]);
                    }
                }
                refine_chirped_tones(block, timestep_s, prediction, bounds, *budget)
                    .map(|result| {
                        anchor.copy_from_slice(&result.tones);
                        *anchor_offset = offset;
                        PartialFrameEstimate::Chirped(result)
                    })
                    .map_err(PartialFitFailure::Chirped)
            }
            _ => unreachable!("warm seeds retain their model kind"),
        };
        if let Err(failure) = &estimate {
            if !failure.unresolved() {
                return Err(PartialTrackingError::Fit {
                    first_sample_index,
                    failure: estimate.expect_err("matched fit failure"),
                });
            }
        } else {
            summary.resolved_frames += 1;
        }
        let frame = PartialFrame {
            first_sample_index,
            sample_count: schedule.frame_samples,
            timestep_s,
            estimate,
        };
        summary.visited_frames += 1;
        summary.completed = offset == last_offset;
        if visit(&frame).is_break() || summary.completed {
            return Ok(summary);
        }
        offset = offset.saturating_add(schedule.hop_samples).min(last_offset);
    }
}

fn reserve_seeds(
    count: usize,
    element_bytes: usize,
    budget: &mut fs_fft::harmonic::HarmonicFitBudget,
) -> Result<(), HarmonicFitError> {
    let required_bytes = count.checked_mul(element_bytes);
    let remaining = required_bytes.and_then(|bytes| budget.max_workspace_bytes.checked_sub(bytes));
    budget.max_workspace_bytes = remaining.ok_or(HarmonicFitError::BudgetExceeded {
        required_bytes: required_bytes.unwrap_or(usize::MAX),
        max_bytes: budget.max_workspace_bytes,
    })?;
    Ok(())
}
