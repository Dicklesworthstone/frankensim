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
    /// including refinement scratch and two candidate coefficient vectors.
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
        /// Unexplained RMS at that iterate, in physical input units.
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
    fit: HarmonicFit,
    normalized_residual_squared: f64,
    correction_hz: f64,
}

impl Candidate {
    fn not_converged(&self, evaluations: usize) -> HarmonicFrequencyError {
        HarmonicFrequencyError::NotConverged {
            fundamental_hz: self.frequency,
            residual_rms: self.fit.residual_rms,
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
    let columns = harmonics
        .len()
        .checked_mul(2)
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| invalid("harmonic design size overflows"))?;
    if samples.len() <= columns {
        return Err(invalid(
            "frequency refinement needs more than 1 + 2 * harmonic_count samples",
        ));
    }
    let duration = (samples.len() - 1) as f64 * timestep_s;
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
    // sample vectors and a previous candidate's phasors in addition, so even
    // preparation of a trial fits inside the declared numeric-array peak.
    let extra_bytes = samples
        .len()
        .checked_mul(2)
        .and_then(|n| {
            harmonics
                .len()
                .checked_mul(4)
                .and_then(|k| n.checked_add(k))
        })
        .and_then(|n| n.checked_mul(core::mem::size_of::<f64>()))
        .ok_or_else(|| invalid("frequency refinement workspace size overflows"))?;
    let required_bytes = super::fit_workspace_bytes(samples.len(), harmonics.len())?
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
        samples.len(),
        timestep_s,
        initial_fundamental_hz,
        harmonics,
        plan_budget,
    )?;
    if samples.iter().any(|x| !x.is_finite()) {
        return Err(HarmonicFitError::NonFinite {
            what: "input sample",
        }
        .into());
    }
    if samples.iter().all(|x| *x == samples[0]) {
        return Err(HarmonicFrequencyError::Unidentifiable);
    }
    let mut scratch = vec![0.0; samples.len()];
    let mut derivative = vec![0.0; samples.len()];
    let mut current = evaluate(samples, seed_plan, duration, &mut scratch, &mut derivative)?;
    let mut evaluations = 1;
    loop {
        if current.correction_hz.abs() <= budget.absolute_tolerance_hz {
            return Ok(HarmonicFrequencyFit {
                fundamental_hz: current.frequency,
                fit: current.fit,
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
                samples.len(),
                timestep_s,
                trial_frequency,
                harmonics,
                plan_budget,
            )?;
            evaluations += 1;
            let trial = evaluate(samples, plan, duration, &mut scratch, &mut derivative)?;
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
fn evaluate(
    samples: &[f64],
    plan: HarmonicFitPlan,
    duration: f64,
    scratch: &mut [f64],
    derivative: &mut [f64],
) -> Result<Candidate, HarmonicFrequencyError> {
    let mut phasors = vec![C64::default(); plan.harmonics.len()];
    let summary = plan.fit_into(samples, scratch, &mut phasors)?;
    let columns = 2 * plan.harmonics.len() + 1;
    let residual_squared: f64 = scratch[columns..].iter().map(|r| r * r).sum();
    let cycles = plan.fundamental_hz * plan.timestep_s();
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
    let derivative_scale = derivative.iter().map(|d| d.abs()).fold(0.0, f64::max);
    if derivative_scale == 0.0 {
        return Err(HarmonicFrequencyError::Unidentifiable);
    }
    if !derivative_scale.is_finite() {
        return Err(HarmonicFitError::NonFinite {
            what: "frequency derivative",
        }
        .into());
    }
    for d in derivative.iter_mut() {
        *d /= derivative_scale;
    }
    plan.design.factor.apply_qt(derivative);
    let denominator: f64 = derivative[columns..].iter().map(|d| d * d).sum();
    if denominator <= samples.len() as f64 * (64.0 * f64::EPSILON).powi(2) {
        return Err(HarmonicFrequencyError::Unidentifiable);
    }
    let numerator: f64 = derivative[columns..]
        .iter()
        .zip(&scratch[columns..])
        .map(|(d, r)| d * r)
        .sum();
    let correction_hz = ((numerator / denominator) / derivative_scale) / duration;
    if !correction_hz.is_finite() || !residual_squared.is_finite() {
        return Err(HarmonicFitError::NonFinite {
            what: "projected frequency correction",
        }
        .into());
    }
    Ok(Candidate {
        frequency: plan.fundamental_hz,
        fit: HarmonicFit {
            dc: summary.dc,
            phasors,
            residual_rms: summary.residual_rms,
        },
        normalized_residual_squared: residual_squared,
        correction_hz,
    })
}
