//! Joint exponential amplitude decay and linear instantaneous-frequency drift.

use super::refinement::{self, Controls, Domain};
use super::{
    C64, ChirpedTone, FitDesign, HarmonicFit, HarmonicFitBudget, HarmonicFitError,
    HarmonicFitSummary,
};
use fs_math::{dd::Dd, det};

/// One decaying linear instantaneous-frequency trajectory.
///
/// With `T = (sample_count - 1) * timestep_s`, its real contribution is
/// `exp(-decay*t) * (a*cos(phase(t)) + b*sin(phase(t)))`, where phase in cycles
/// is `start_frequency_hz*t + (end_frequency_hz-start_frequency_hz)*t*t/(2*T)`.
/// Decay acts on amplitude. Both frequencies refer to actual sampled endpoints;
/// a changing frequency is not a stationary modal pole or a natural frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedChirpedTone {
    /// Positive instantaneous frequency at the first sample [Hz].
    pub start_frequency_hz: f64,
    /// Positive instantaneous frequency at the last sample [Hz].
    pub end_frequency_hz: f64,
    /// Nonnegative exponential amplitude decay rate [1/s].
    pub decay_rate_per_s: f64,
}

/// Reusable joint DC/amplitude fit at supplied decaying frequency trajectories.
///
/// Peak phasors `a + i*b` use `exp(-i phase)` and refer to the first sample,
/// before exponential attenuation. Tone order is retained; trajectories may
/// cross or share both endpoint frequencies if their decays distinguish them.
/// The complete sampled design must pass the declared conditioning floor.
///
/// This supplied-law plan does not infer frequency, decay, mode count or
/// uncertainty. Nonlinear drift, nonexponential envelopes and omitted tones
/// remain in the residual and may bias the fit. It shares the existing QR,
/// physical-unit scaling and allocation-free `fit_into` owner. Zero decay
/// preserves [`super::ChirpedToneFitPlan`] bits; zero drift preserves
/// [`super::DampedToneFitPlan`] bits for corresponding declarations and records.
#[derive(Debug, Clone)]
pub struct DampedChirpedToneFitPlan {
    pub(super) tones: Vec<DampedChirpedTone>,
    pub(super) design: FitDesign,
}

impl DampedChirpedToneFitPlan {
    /// Prepare distinct supplied trajectories on one uniform physical clock.
    /// Both frequency endpoints must be positive and strictly sub-Nyquist;
    /// linearity then admits the whole instantaneous-frequency path. Decays
    /// must be finite and nonnegative. At least `1 + 2*tones.len()` samples
    /// and a finite positive duration are required. Storage includes all three
    /// retained physical parameters per tone before allocating numeric arrays.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        tones: &[DampedChirpedTone],
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if tones.is_empty() || sample_count < 2 {
            return Err(super::invalid(
                "damped chirp fitting needs tones and at least two samples",
            ));
        }
        let duration = (sample_count - 1) as f64 * timestep_s;
        if !(timestep_s.is_finite() && timestep_s > 0.0 && duration.is_finite() && duration > 0.0) {
            return Err(super::invalid(
                "damped chirp clock and duration must be positive and finite",
            ));
        }
        for (k, tone) in tones.iter().enumerate() {
            if ![tone.start_frequency_hz, tone.end_frequency_hz]
                .into_iter()
                .all(|f| f.is_finite() && f > 0.0 && f * timestep_s > 0.0 && f * timestep_s < 0.5)
                || !(tone.decay_rate_per_s.is_finite() && tone.decay_rate_per_s >= 0.0)
            {
                return Err(super::invalid(
                    "damped chirps need finite positive sub-Nyquist endpoints and nonnegative decay",
                ));
            }
            if tones[..k].contains(tone) {
                return Err(super::invalid(
                    "damped chirp trajectories must be distinct in frequency or decay",
                ));
            }
        }
        // The shared owner already budgets one retained f64 per tone.
        let extra_bytes = tones
            .len()
            .checked_mul(2 * core::mem::size_of::<f64>())
            .ok_or_else(|| super::invalid("damped chirp parameter storage size overflows"))?;
        let inverse_intervals = Dd::ONE / Dd::from_f64((sample_count - 1) as f64);
        let design = FitDesign::with_basis(
            sample_count,
            timestep_s,
            tones.len(),
            |k| tones[k].start_frequency_hz.max(tones[k].end_frequency_hz),
            |j, k| {
                let (s, c) = phase_sin_cos(tones[k], timestep_s, inverse_intervals, j);
                let amplitude = envelope(tones[k].decay_rate_per_s, timestep_s, j);
                (s * amplitude, c * amplitude)
            },
            extra_bytes,
            budget,
        )?;
        Ok(Self {
            tones: tones.to_vec(),
            design,
        })
    }

    /// Supplied physical parameters in the same order as fitted phasors.
    #[must_use]
    pub fn tones(&self) -> &[DampedChirpedTone] {
        &self.tones
    }

    /// Required record and scratch length.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.design.samples
    }

    /// Uniform physical sample interval [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.design.timestep_s
    }

    /// Reciprocal 1-norm condition of the complete linear amplitude design.
    #[must_use]
    pub const fn reciprocal_condition(&self) -> f64 {
        self.design.reciprocal_condition
    }

    /// Conservative peak numeric-array storage admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.design.workspace_bytes
    }

    /// Fit one record, allocating scratch and first-sample peak phasors.
    pub fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        self.design.fit(samples)
    }

    /// Allocation-free fit; output phasors remain unchanged on refusal.
    pub fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        self.design.fit_into(samples, scratch, phasors)
    }
}

pub(super) fn phase_sin_cos(
    tone: DampedChirpedTone,
    dt: f64,
    inverse_intervals: Dd,
    index: usize,
) -> (f64, f64) {
    super::chirp::sample_sin_cos(
        ChirpedTone {
            start_frequency_hz: tone.start_frequency_hz,
            end_frequency_hz: tone.end_frequency_hz,
        },
        dt,
        inverse_intervals,
        index,
    )
}

pub(super) fn envelope(decay: f64, dt: f64, index: usize) -> f64 {
    if index == 0 || decay == 0.0 {
        1.0
    } else {
        // Preserve the damped-tone owner's operation order. Overflow of a
        // positive exponent gives zero amplitude; t=0 is handled above.
        det::exp(-((decay * dt) * index as f64))
    }
}

/// Closed search box for one decaying frequency trajectory.
/// Boxes must be strictly disjoint in at least one of the three coordinates;
/// distinct decay intervals may therefore separate identical frequency laws.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedChirpedToneBounds {
    /// Positive, increasing first-sample frequency interval [Hz].
    pub start_frequency_hz: [f64; 2],
    /// Positive, increasing last-sample frequency interval [Hz].
    pub end_frequency_hz: [f64; 2],
    /// Nonnegative, increasing amplitude-decay interval [1/s].
    pub decay_rate_per_s: [f64; 2],
}

/// Numerical resolution, work and storage for joint decay/drift identification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedChirpedToneRefinementBudget {
    /// Both untruncated endpoint corrections must meet this positive limit [Hz].
    pub absolute_frequency_tolerance_hz: f64,
    /// Every untruncated decay correction must meet this positive limit [1/s].
    pub absolute_decay_tolerance_per_s: f64,
    /// All candidate trajectory sets, including rejected backtracking trials.
    pub max_evaluations: usize,
    /// Minimum reciprocal 1-norm condition of both full and projected physical
    /// Jacobians in coordinates `[f_start*T, f_end*T, decay*T]`. Original
    /// common-unit channel amplitudes survive; columns are not independently
    /// normalized. Resolved projected sensitivity is also required, so residual
    /// curvature cannot identify an absent component. This is not uncertainty.
    pub min_trajectory_reciprocal_condition: f64,
    /// Linear-design conditioning and total peak numeric-array storage,
    /// including refinement scratch and both candidates' channel fits.
    pub fit: HarmonicFitBudget,
}

/// Signed, untruncated physical corrections, not confidence intervals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedChirpedToneCorrection {
    /// First-sample frequency correction [Hz].
    pub start_frequency_hz: f64,
    /// Last-sample frequency correction [Hz].
    pub end_frequency_hz: f64,
    /// Amplitude-decay correction [1/s]; may be negative.
    pub decay_rate_per_s: f64,
}

/// Jointly inferred frequency trajectories, decay rates and initial amplitudes.
#[derive(Debug, Clone, PartialEq)]
pub struct DampedChirpedToneFit {
    /// Physical parameters in the caller's original seed order.
    pub tones: Vec<DampedChirpedTone>,
    /// DC, first-sample peak phasors and unexplained sampled RMS in input units.
    pub fit: HarmonicFit,
    /// Evaluated candidate sets, including rejected trials.
    pub evaluations: usize,
    /// Remaining untruncated corrections in the same tone order.
    pub last_corrections: Vec<DampedChirpedToneCorrection>,
    /// Minimum condition of the full and projected physical Jacobians.
    pub trajectory_reciprocal_condition: f64,
}

/// Shared decaying frequency laws with independent simultaneous-channel phasors.
#[derive(Debug, Clone, PartialEq)]
pub struct DampedChirpedToneChannelsFit {
    /// Shared physical parameters in the original seed order.
    pub tones: Vec<DampedChirpedTone>,
    /// Original channel order, input units, DC and first-sample peak phasors.
    pub fits: Vec<HarmonicFit>,
    /// Evaluated candidate sets; each shares one amplitude QR across channels.
    pub evaluations: usize,
    /// Remaining untruncated corrections in the same tone order.
    pub last_corrections: Vec<DampedChirpedToneCorrection>,
    /// Minimum condition of the full and projected physical Jacobians.
    pub trajectory_reciprocal_condition: f64,
}

/// Refusal to identify unresolved decay/drift or publish unconverged amplitudes.
#[derive(Debug, Clone, PartialEq)]
pub enum DampedChirpedToneRefinementError {
    /// The original amplitude owner refused data, conditioning or storage.
    Fit(HarmonicFitError),
    /// Invalid clock, local domains, channel shapes or controls.
    InvalidParameter {
        /// The failed admission condition.
        what: &'static str,
    },
    /// Joint decay/drift derivatives lack resolved physical signal or rank.
    Unidentifiable {
        /// Refused full/projected Jacobian reciprocal 1-norm condition.
        reciprocal_condition: f64,
        /// Requested minimum reciprocal condition.
        minimum: f64,
    },
    /// Work exhaustion, active bounds or stagnation left excessive corrections.
    /// No partial fitted amplitudes are published.
    NotConverged {
        /// Best accepted physical trajectories, in seed order.
        tones: Vec<DampedChirpedTone>,
        /// Joint unexplained RMS in the channels' common physical units.
        residual_rms: f64,
        /// Evaluated candidate sets, including rejected trials.
        evaluations: usize,
        /// Remaining untruncated corrections at that candidate.
        last_corrections: Vec<DampedChirpedToneCorrection>,
    },
}

impl core::fmt::Display for DampedChirpedToneRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(f),
            Self::InvalidParameter { what } => write!(f, "damped chirped tones: {what}"),
            Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => write!(
                f,
                "damped chirped tones are unidentifiable: reciprocal condition \
                 {reciprocal_condition:e} is below {minimum:e}"
            ),
            Self::NotConverged { evaluations, .. } => write!(
                f,
                "damped chirped tones did not converge after {evaluations} candidate sets"
            ),
        }
    }
}

impl std::error::Error for DampedChirpedToneRefinementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            _ => None,
        }
    }
}

fn tones(parameters: &[f64]) -> Vec<DampedChirpedTone> {
    parameters
        .chunks_exact(3)
        .map(|p| DampedChirpedTone {
            start_frequency_hz: p[0],
            end_frequency_hz: p[1],
            decay_rate_per_s: p[2],
        })
        .collect()
}

fn corrections(parameters: &[f64]) -> Vec<DampedChirpedToneCorrection> {
    parameters
        .chunks_exact(3)
        .map(|p| DampedChirpedToneCorrection {
            start_frequency_hz: p[0],
            end_frequency_hz: p[1],
            decay_rate_per_s: p[2],
        })
        .collect()
}

impl From<refinement::Error> for DampedChirpedToneRefinementError {
    fn from(error: refinement::Error) -> Self {
        match error {
            refinement::Error::Fit(error) => Self::Fit(error),
            refinement::Error::InvalidParameter { what } => Self::InvalidParameter { what },
            refinement::Error::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            },
            refinement::Error::NotConverged {
                parameters,
                residual_rms,
                evaluations,
                corrections: delta,
            } => Self::NotConverged {
                tones: tones(&parameters),
                residual_rms,
                evaluations,
                last_corrections: corrections(&delta),
            },
        }
    }
}

fn invalid(what: &'static str) -> DampedChirpedToneRefinementError {
    DampedChirpedToneRefinementError::InvalidParameter { what }
}

/// Jointly refine locally seeded exponential decay and linear frequency drift.
///
/// All three physical parameters are inferred together. First/last endpoint
/// frequencies remain in hertz, amplitude decay in reciprocal seconds, and
/// phasors in first-sample input units. No stationary-frequency or constant-
/// amplitude approximation is inserted between frequency and decay fits.
/// The complete variable-projection derivative and independent projected-
/// sensitivity gate retain the shared QR owner's original physics and weights.
/// Active bounds, unresolved parameters and excessive corrections refuse.
///
/// This is local inference of a declared model order and law. It does not
/// discover modes, fit arbitrary envelopes/nonlinear drift, infer material
/// damping, certify uncertainty or remove aliasing. Omitted motion remains in
/// the residual and can bias the fit. A changing frequency is not a modal pole.
pub fn refine_damped_chirped_tones(
    samples: &[f64],
    timestep_s: f64,
    initial_tones: &[DampedChirpedTone],
    bounds: &[DampedChirpedToneBounds],
    budget: DampedChirpedToneRefinementBudget,
) -> Result<DampedChirpedToneFit, DampedChirpedToneRefinementError> {
    let mut result = refine_damped_chirped_tones_channels(
        &[samples],
        timestep_s,
        initial_tones,
        bounds,
        budget,
    )?;
    Ok(DampedChirpedToneFit {
        tones: result.tones,
        fit: result.fits.remove(0),
        evaluations: result.evaluations,
        last_corrections: result.last_corrections,
        trajectory_reciprocal_condition: result.trajectory_reciprocal_condition,
    })
}

/// Infer shared decaying frequency laws from simultaneous common-clock channels.
/// Each channel keeps its own DC and initial phasors. Original common-unit
/// channel weights survive orthogonal changes of vector components; callers
/// must convert unrelated units explicitly. Constant channels retain DC without
/// identifying motion. All three coordinates need resolved joint sensitivity.
///
/// One amplitude QR serves each candidate. Full and projected derivative blocks
/// fold through small QRs, without normal equations or a stacked channel matrix.
/// Work is O(E C (N K² + K³)); numeric storage is O(N K + K² + C K).
/// The offline owner has no mid-factorization cancellation or throughput claim.
pub fn refine_damped_chirped_tones_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_tones: &[DampedChirpedTone],
    bounds: &[DampedChirpedToneBounds],
    budget: DampedChirpedToneRefinementBudget,
) -> Result<DampedChirpedToneChannelsFit, DampedChirpedToneRefinementError> {
    if initial_tones.is_empty() || initial_tones.len() != bounds.len() {
        return Err(invalid(
            "each damped chirp needs one three-parameter search box",
        ));
    }
    if !(timestep_s.is_finite() && timestep_s > 0.0)
        || ![
            budget.absolute_frequency_tolerance_hz,
            budget.absolute_decay_tolerance_per_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        || budget.max_evaluations == 0
        || !(budget.min_trajectory_reciprocal_condition.is_finite()
            && budget.min_trajectory_reciprocal_condition > 0.0
            && budget.min_trajectory_reciprocal_condition <= 1.0)
    {
        return Err(invalid(
            "invalid damped chirp clock, tolerances, work or condition limit",
        ));
    }
    for (k, (seed, bound)) in initial_tones.iter().zip(bounds).enumerate() {
        for (value, [lo, hi]) in [
            (seed.start_frequency_hz, bound.start_frequency_hz),
            (seed.end_frequency_hz, bound.end_frequency_hz),
        ] {
            if ![value, lo, hi].into_iter().all(f64::is_finite)
                || !(lo > 0.0 && lo < hi && lo * timestep_s > 0.0 && hi * timestep_s < 0.5)
                || !(value >= lo && value <= hi)
            {
                return Err(invalid(
                    "damped chirp frequency boxes must contain seeds and stay positive and sub-Nyquist",
                ));
            }
        }
        let [lo, hi] = bound.decay_rate_per_s;
        if ![seed.decay_rate_per_s, lo, hi]
            .into_iter()
            .all(f64::is_finite)
            || !(lo >= 0.0 && lo < hi)
            || !(seed.decay_rate_per_s >= lo && seed.decay_rate_per_s <= hi)
        {
            return Err(invalid(
                "damped chirp decay intervals must contain seeds and be nonnegative",
            ));
        }
        if bounds[..k].iter().any(|other| {
            !(bound.start_frequency_hz[1] < other.start_frequency_hz[0]
                || other.start_frequency_hz[1] < bound.start_frequency_hz[0]
                || bound.end_frequency_hz[1] < other.end_frequency_hz[0]
                || other.end_frequency_hz[1] < bound.end_frequency_hz[0]
                || bound.decay_rate_per_s[1] < other.decay_rate_per_s[0]
                || other.decay_rate_per_s[1] < bound.decay_rate_per_s[0])
        }) {
            return Err(invalid(
                "damped chirp boxes must be strictly disjoint in frequency or decay",
            ));
        }
    }
    let result = refinement::refine_channels(
        channels,
        timestep_s,
        Domain::DampedChirped(initial_tones, bounds),
        Controls {
            tolerances: [
                budget.absolute_frequency_tolerance_hz,
                budget.absolute_decay_tolerance_per_s,
            ],
            max_evaluations: budget.max_evaluations,
            min_reciprocal_condition: budget.min_trajectory_reciprocal_condition,
            fit: budget.fit,
        },
    )?;
    Ok(DampedChirpedToneChannelsFit {
        tones: tones(&result.parameters),
        fits: result.fits,
        evaluations: result.evaluations,
        last_corrections: corrections(&result.corrections),
        trajectory_reciprocal_condition: result.reciprocal_condition,
    })
}
