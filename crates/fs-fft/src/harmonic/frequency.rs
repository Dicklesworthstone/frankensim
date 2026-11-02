//! Bounded local variable projection for a common harmonic fundamental.

use super::{C64, HarmonicFit, HarmonicFitBudget, HarmonicFitError, HarmonicFitPlan};
use fs_math::det;

/// Numerical limits for local fundamental-frequency refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HarmonicFrequencyBudget {
    /// Positive, increasing local search bounds, inclusive, in hertz. Every
    /// retained harmonic must remain strictly below Nyquist throughout them.
    pub frequency_bounds_hz: [f64; 2],
    /// Convergence requires the untruncated projected correction to be no
    /// larger than this positive tolerance. This is not frequency uncertainty.
    pub absolute_tolerance_hz: f64,
    /// Maximum QR fits, including the seed and rejected backtracking trials.
    pub max_evaluations: usize,
    /// Conditioning floor and conservative total numeric-array peak storage,
    /// including refinement scratch and every channel's two candidate fits.
    pub fit: HarmonicFitBudget,
}

/// A local stationary harmonic fit, referenced to the first input sample.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicFrequencyFit {
    /// Refined common fundamental in hertz, within the admitted search bounds.
    pub fundamental_hz: f64,
    /// Joint DC and peak phasors at that frequency, in physical input units.
    pub fit: HarmonicFit,
    /// Number of QR fits performed, including rejected candidate steps.
    pub evaluations: usize,
    /// Untruncated projected correction at the returned frequency, in hertz.
    pub last_correction_hz: f64,
}

/// One common fundamental with independently fitted DC and phasors for each
/// simultaneous channel, retaining the original channel and harmonic order.
#[derive(Debug, Clone, PartialEq)]
pub struct HarmonicFrequencyChannelsFit {
    /// Refined fundamental shared by every channel [Hz].
    pub fundamental_hz: f64,
    /// Per-channel fits in the caller's common physical units. A constant
    /// channel retains its DC and zero phasors without inventing frequency data.
    pub fits: Vec<HarmonicFit>,
    /// Candidate frequencies evaluated; all channels reuse each candidate QR.
    pub evaluations: usize,
    /// Remaining untruncated projected correction [Hz].
    pub last_correction_hz: f64,
}

/// Refusal to manufacture a refined frequency or partial spectrum.
#[derive(Debug, Clone, PartialEq)]
pub enum HarmonicFrequencyError {
    /// The shared joint-fit owner refused sampling, conditioning, or storage.
    Fit(HarmonicFitError),
    /// The local frequency interval, tolerance, or evaluation cap is invalid.
    InvalidParameter {
        /// Admission condition that failed.
        what: &'static str,
    },
    /// The record or projected harmonic derivative carries no numerically
    /// resolved frequency information. A constant record has no pitch.
    Unidentifiable,
    /// Budget exhaustion, an active bound, or floating-point stagnation left
    /// the frequency correction above the requested tolerance. These fields
    /// describe the best accepted iterate; no partial phasors are returned.
    NotConverged {
        /// Best accepted frequency in hertz.
        fundamental_hz: f64,
        /// Unexplained RMS over all scalar samples/channels at that iterate,
        /// in the channels' common physical input units.
        residual_rms: f64,
        /// QR fits performed, including rejected trials.
        evaluations: usize,
        /// Remaining untruncated projected correction in hertz.
        last_correction_hz: f64,
    },
}

impl From<HarmonicFitError> for HarmonicFrequencyError {
    fn from(error: HarmonicFitError) -> Self {
        Self::Fit(error)
    }
}

impl core::fmt::Display for HarmonicFrequencyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(f),
            Self::InvalidParameter { what } => write!(f, "harmonic frequency: {what}"),
            Self::Unidentifiable => write!(f, "harmonic frequency is not identifiable"),
            Self::NotConverged {
                fundamental_hz,
                evaluations,
                last_correction_hz,
                ..
            } => write!(
                f,
                "harmonic frequency did not converge after {evaluations} fits: \
                 {fundamental_hz} Hz, remaining correction {last_correction_hz:e} Hz"
            ),
        }
    }
}

impl std::error::Error for HarmonicFrequencyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            _ => None,
        }
    }
}

fn invalid(what: &'static str) -> HarmonicFrequencyError {
    HarmonicFrequencyError::InvalidParameter { what }
}

struct Candidate {
    frequency: f64,
    fits: Vec<HarmonicFit>,
    residual_rms: f64,
    normalized_residual_squared: f64,
    correction_hz: f64,
}

impl Candidate {
    fn not_converged(&self, evaluations: usize) -> HarmonicFrequencyError {
        HarmonicFrequencyError::NotConverged {
            fundamental_hz: self.frequency,
            residual_rms: self.residual_rms,
            evaluations,
            last_correction_hz: self.correction_hz,
        }
    }
}

/// Refine a supplied local pitch seed while refitting DC and all harmonic
/// amplitudes together. Missing fundamentals and weak upper harmonics are
/// allowed: frequency information comes from all declared harmonic columns.
///
/// Variable projection eliminates the linear coefficients with the same
/// Householder QR as [`HarmonicFitPlan`]. The analytic sinusoid derivative is
/// projected orthogonally to those columns. A scalar projected Gauss--Newton
/// step uses monotone backtracking and moves at most a quarter fundamental
/// cycle across the record per trial. No normal equations, FFT-bin rounding,
/// or finite-difference derivative are used. Every attempted QR is budgeted.
///
/// This is a local stationary-spectrum estimator, not a pitch detector, global
/// optimizer, octave selector, time-varying tracker, or uncertainty certificate.
/// The caller supplies an independent seed and physically meaningful interval.
/// Noise, inharmonicity and changing frequency remain in the sampled residual.
/// Bound clipping, exhausted work, or stagnation cannot count as convergence.
/// Preparation costs O(E (N K^2 + K^3)) with O(N K) peak storage for at most E
/// fits; it is an offline observation path, not an allocation-free audio step.
pub fn refine_harmonic_frequency(
    samples: &[f64],
    timestep_s: f64,
    initial_fundamental_hz: f64,
    harmonics: &[usize],
    budget: HarmonicFrequencyBudget,
) -> Result<HarmonicFrequencyFit, HarmonicFrequencyError> {
    let mut result = refine_harmonic_frequency_channels(
        &[samples],
        timestep_s,
        initial_fundamental_hz,
        harmonics,
        budget,
    )?;
    Ok(HarmonicFrequencyFit {
        fundamental_hz: result.fundamental_hz,
        fit: result.fits.remove(0),
        evaluations: result.evaluations,
        last_correction_hz: result.last_correction_hz,
    })
}

/// Refine one fundamental from simultaneous equal-length, common-clock
/// channels. Each channel has independent DC/phasors, while one QR per candidate
/// supplies every fit and projected derivative. This is the same numerical
/// owner as [`refine_harmonic_frequency`], not an average of separately inferred
/// pitches or a choice of one arbitrarily preferred channel.
///
/// The objective is the sum of squared residuals in the channels' COMMON units,
/// with equal weights. Relative channel strength is retained; independently
/// normalizing channels to unit amplitude would change the physical objective.
/// In particular, Cartesian vector components retain their Euclidean metric
/// under orthogonal frame changes. Unrelated units or sensor noise weights must
/// be converted explicitly by the caller before use. DC preloads do not set the
/// amplitude scale. Constant channels carry no frequency information but are
/// retained, and only an entirely unresolved record is unidentifiable.
///
/// Storage is O(N K + C K + N) for N samples, K orders and C channels. Work is
/// O(E (N K^2 + C N K + K^3)); the evaluation cap counts common QR candidates.
/// Channels stream through two reusable length-N buffers. The same local,
/// stationary, sub-Nyquist and non-certification boundaries apply as for the
/// scalar estimator. The memory budget includes every retained channel fit.
pub fn refine_harmonic_frequency_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_fundamental_hz: f64,
    harmonics: &[usize],
    budget: HarmonicFrequencyBudget,
) -> Result<HarmonicFrequencyChannelsFit, HarmonicFrequencyError> {
    let [lower, upper] = budget.frequency_bounds_hz;
    if !(lower.is_finite() && upper.is_finite() && lower > 0.0 && lower < upper) {
        return Err(invalid(
            "frequency bounds must be finite, positive and increasing",
        ));
    }
    if !initial_fundamental_hz.is_finite()
        || initial_fundamental_hz < lower
        || initial_fundamental_hz > upper
    {
        return Err(invalid(
            "initial frequency must lie within the search bounds",
        ));
    }
    if !(budget.absolute_tolerance_hz.is_finite() && budget.absolute_tolerance_hz > 0.0)
        || budget.max_evaluations == 0
    {
        return Err(invalid(
            "frequency tolerance and evaluation cap must be positive",
        ));
    }
    let Some(samples) = channels.first() else {
        return Err(invalid("at least one simultaneous channel is required"));
    };
    let sample_count = samples.len();
    if channels.iter().any(|channel| channel.len() != sample_count) {
        return Err(invalid(
            "simultaneous channels must have the same sample count",
        ));
    }
    if channels
        .iter()
        .any(|channel| channel.iter().any(|x| !x.is_finite()))
    {
        return Err(HarmonicFitError::NonFinite {
            what: "input sample",
        }
        .into());
    }
    let columns = harmonics
        .len()
        .checked_mul(2)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("harmonic design size overflows"))?;
    if sample_count <= columns {
        return Err(invalid(
            "frequency refinement needs more than 1 + 2 * harmonic_count samples",
        ));
    }
    let duration = (sample_count - 1) as f64 * timestep_s;
    if !(duration.is_finite() && duration > 0.0) {
        return Err(invalid("record duration must be positive and finite"));
    }
    if harmonics.iter().any(|&h| {
        let cycles = (upper * timestep_s) * h as f64;
        !(cycles > 0.0 && cycles < 0.5 && (upper * h as f64).is_finite())
    }) {
        return Err(invalid(
            "all harmonics must remain below Nyquist over the search interval",
        ));
    }
    // The fixed-fit owner budgets a plan and one fit. Reserve both persistent
    // sample vectors and both candidates' channel fits/scales in addition, so
    // even preparation of a trial stays inside the declared storage peak.
    let channel_bytes = harmonics
        .len()
        .checked_mul(4 * core::mem::size_of::<f64>())
        .and_then(|n| n.checked_add(2 * core::mem::size_of::<HarmonicFit>()))
        .and_then(|n| n.checked_add(core::mem::size_of::<f64>()))
        .and_then(|n| n.checked_mul(channels.len()))
        .ok_or_else(|| invalid("channel fit workspace size overflows"))?;
    let extra_bytes = sample_count
        .checked_mul(2 * core::mem::size_of::<f64>())
        .and_then(|n| n.checked_add(channel_bytes))
        .ok_or_else(|| invalid("frequency refinement workspace size overflows"))?;
    let required_bytes = super::fit_workspace_bytes(sample_count, harmonics.len())?
        .checked_add(extra_bytes)
        .ok_or_else(|| invalid("frequency refinement workspace size overflows"))?;
    if required_bytes > budget.fit.max_workspace_bytes {
        return Err(HarmonicFitError::BudgetExceeded {
            required_bytes,
            max_bytes: budget.fit.max_workspace_bytes,
        }
        .into());
    }
    let plan_budget = HarmonicFitBudget {
        max_workspace_bytes: budget.fit.max_workspace_bytes - extra_bytes,
        ..budget.fit
    };
    // Validate the seed's harmonic declarations/conditioning before allocating
    // scratch. Only one QR plan is retained at any time, including line search.
    let seed_plan = HarmonicFitPlan::new(
        sample_count,
        timestep_s,
        initial_fundamental_hz,
        harmonics,
        plan_budget,
    )?;
    let mut ratios: Vec<f64> = channels
        .iter()
        .map(|channel| {
            let low = channel.iter().copied().fold(f64::INFINITY, f64::min);
            let high = channel.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let center = low.midpoint(high);
            channel
                .iter()
                .map(|x| (x - center).abs())
                .fold(0.0, f64::max)
        })
        .collect();
    let physical_scale = ratios.iter().copied().fold(0.0, f64::max);
    if physical_scale == 0.0 {
        return Err(HarmonicFrequencyError::Unidentifiable);
    }
    for ratio in &mut ratios {
        *ratio /= physical_scale;
    }
    let mut scratch = vec![0.0; sample_count];
    let mut derivative = vec![0.0; sample_count];
    let mut current = evaluate(
        channels,
        &ratios,
        physical_scale,
        seed_plan,
        duration,
        &mut scratch,
        &mut derivative,
    )?;
    let mut evaluations = 1;
    loop {
        if current.correction_hz.abs() <= budget.absolute_tolerance_hz {
            return Ok(HarmonicFrequencyChannelsFit {
                fundamental_hz: current.frequency,
                fits: current.fits,
                evaluations,
                last_correction_hz: current.correction_hz,
            });
        }
        let max_step_hz = 0.25 / duration;
        let mut step = current.correction_hz.clamp(-max_step_hz, max_step_hz);
        loop {
            if evaluations == budget.max_evaluations {
                return Err(current.not_converged(evaluations));
            }
            let trial_frequency = (current.frequency + step).clamp(lower, upper);
            if trial_frequency == current.frequency {
                return Err(current.not_converged(evaluations));
            }
            let plan = HarmonicFitPlan::new(
                sample_count,
                timestep_s,
                trial_frequency,
                harmonics,
                plan_budget,
            )?;
            evaluations += 1;
            let trial = evaluate(
                channels,
                &ratios,
                physical_scale,
                plan,
                duration,
                &mut scratch,
                &mut derivative,
            )?;
            if trial.normalized_residual_squared < current.normalized_residual_squared {
                current = trial;
                break;
            }
            step *= 0.5;
        }
    }
}

// Consumes/drops each plan before the next is prepared. Scratch after fit_into
// holds normalized solved coefficients followed by the Q^T residual. Thus all
// derivative/line-search arithmetic is independent of the input's unit scale.
#[allow(clippy::too_many_arguments)]
fn evaluate(
    channels: &[&[f64]],
    ratios: &[f64],
    physical_scale: f64,
    plan: HarmonicFitPlan,
    duration: f64,
    scratch: &mut [f64],
    derivative: &mut [f64],
) -> Result<Candidate, HarmonicFrequencyError> {
    let columns = 2 * plan.harmonics.len() + 1;
    let cycles = plan.fundamental_hz * plan.timestep_s();
    let mut fits = Vec::with_capacity(channels.len());
    let (mut residual_squared, mut numerator, mut denominator, mut joint_scale) =
        (0.0, 0.0, 0.0, 0.0_f64);
    for (&samples, &ratio) in channels.iter().zip(ratios) {
        let mut phasors = vec![C64::default(); plan.harmonics.len()];
        let summary = plan.fit_into(samples, scratch, &mut phasors)?;
        fits.push(HarmonicFit {
            dc: summary.dc,
            phasors,
            residual_rms: summary.residual_rms,
        });
        // A constant channel's fit_into returns without filling scratch; it has
        // exactly zero residual and supplies no derivative. A ratio below f64
        // representation also cannot affect this common-scale objective.
        if ratio == 0.0 {
            continue;
        }
        let local_residual: f64 = scratch[columns..].iter().map(|r| r * r).sum();
        residual_squared += local_residual * ratio * ratio;
        for (j, out) in derivative.iter_mut().enumerate() {
            // Centering the derivative's time coordinate changes it only by a
            // fitted sinusoid column, which its orthogonal projection removes.
            let centered_time = j as f64 / (samples.len() - 1) as f64 - 0.5;
            let mut value = 0.0;
            for (k, &h) in plan.harmonics.iter().enumerate() {
                let theta = core::f64::consts::TAU * ((cycles * h as f64) * j as f64);
                let (s, c) = det::sin_cos(theta);
                let slope = (-scratch[2 * k + 1]).mul_add(s, scratch[2 * k + 2] * c);
                value = (h as f64 * slope).mul_add(centered_time, value);
            }
            *out = core::f64::consts::TAU * value;
        }
        if derivative.iter().any(|d| !d.is_finite()) {
            return Err(HarmonicFitError::NonFinite {
                what: "frequency derivative",
            }
            .into());
        }
        let derivative_scale = derivative.iter().map(|d| d.abs()).fold(0.0, f64::max);
        let weighted_scale = derivative_scale * ratio;
        if weighted_scale == 0.0 {
            continue;
        }
        for d in derivative.iter_mut() {
            *d /= derivative_scale;
        }
        plan.design.factor.apply_qt(derivative);
        let local_denominator: f64 = derivative[columns..].iter().map(|d| d * d).sum();
        if local_denominator <= samples.len() as f64 * (64.0 * f64::EPSILON).powi(2) {
            continue;
        }
        let local_numerator: f64 = derivative[columns..]
            .iter()
            .zip(&scratch[columns..])
            .map(|(d, r)| d * r)
            .sum();
        // Maintain one derivative scale while streaming channels. Both the
        // gradient and information keep their original amplitude weights; no
        // independently normalized channel can dominate solely by being weak.
        if weighted_scale > joint_scale {
            let change = joint_scale / weighted_scale;
            numerator *= change;
            denominator *= change * change;
            joint_scale = weighted_scale;
        }
        let weight = weighted_scale / joint_scale;
        denominator += local_denominator * weight * weight;
        numerator += local_numerator * weight * ratio;
    }
    if joint_scale == 0.0 || denominator == 0.0 {
        return Err(HarmonicFrequencyError::Unidentifiable);
    }
    let correction_hz = ((numerator / denominator) / joint_scale) / duration;
    let residual_rms =
        det::sqrt(residual_squared / scratch.len() as f64 / channels.len() as f64) * physical_scale;
    if !correction_hz.is_finite() || !residual_squared.is_finite() || !residual_rms.is_finite() {
        return Err(HarmonicFitError::NonFinite {
            what: "projected frequency correction",
        }
        .into());
    }
    Ok(Candidate {
        frequency: plan.fundamental_hz,
        fits,
        residual_rms,
        normalized_residual_squared: residual_squared,
        correction_hz,
    })
}
