//! Joint local identification of independent stationary tone frequencies.

use super::{C64, HarmonicFit, HarmonicFitBudget, HarmonicFitError, ToneFitPlan};
use fs_la::factor::{Qr, qr};
use fs_math::det;

/// Numerical limits for independent-frequency variable projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ToneFrequencyBudget {
    /// Every untruncated frequency correction must meet this positive limit
    /// [Hz]. This is a convergence criterion, not measurement uncertainty.
    pub absolute_tolerance_hz: f64,
    /// Maximum candidate frequency sets, including rejected line-search trials.
    pub max_evaluations: usize,
    /// Minimum reciprocal 1-norm condition of the joint projected frequency
    /// Jacobian's R factor. All columns use the same physical amplitude scale:
    /// weak and inseparable frequencies cannot hide behind column normalization.
    pub min_frequency_reciprocal_condition: f64,
    /// Linear-fit conditioning and total conservative numeric storage cap,
    /// including all refinement scratch and both candidates' channel fits.
    pub fit: HarmonicFitBudget,
}

/// Locally refined independent tones, referenced to the first input sample.
#[derive(Debug, Clone, PartialEq)]
pub struct ToneFrequencyFit {
    /// Refined physical frequencies [Hz], in the supplied tone order.
    pub frequencies_hz: Vec<f64>,
    /// Joint DC, peak negative-time phasors and unexplained sampled RMS.
    pub fit: HarmonicFit,
    /// Number of candidate frequency sets evaluated.
    pub evaluations: usize,
    /// Remaining untruncated projected corrections [Hz].
    pub last_corrections_hz: Vec<f64>,
    /// Reciprocal condition of the final projected frequency Jacobian.
    pub frequency_reciprocal_condition: f64,
}

/// One independent-frequency set fitted to simultaneous common-unit channels.
#[derive(Debug, Clone, PartialEq)]
pub struct ToneFrequencyChannelsFit {
    /// Refined frequencies [Hz], shared by every channel.
    pub frequencies_hz: Vec<f64>,
    /// Per-channel DC and phasors in the original channel and tone order.
    pub fits: Vec<HarmonicFit>,
    /// Number of shared candidate frequency sets evaluated.
    pub evaluations: usize,
    /// Remaining untruncated projected corrections [Hz].
    pub last_corrections_hz: Vec<f64>,
    /// Reciprocal condition of the final projected frequency Jacobian.
    pub frequency_reciprocal_condition: f64,
}

/// Refusal to infer an unresolved or unconverged independent spectrum.
#[derive(Debug, Clone, PartialEq)]
pub enum ToneFrequencyError {
    /// The existing joint-fit owner refused data, conditioning or storage.
    Fit(HarmonicFitError),
    /// Invalid frequency intervals, clock, channel shape or numerical limits.
    InvalidParameter {
        /// Failed admission condition.
        what: &'static str,
    },
    /// The projected frequency derivatives have insufficient numerical rank.
    /// A constant record and unresolved derivatives report zero condition.
    Unidentifiable {
        /// Reciprocal 1-norm condition of the projected R factor.
        reciprocal_condition: f64,
        /// Requested minimum reciprocal condition.
        minimum: f64,
    },
    /// The best accepted iterate still has an excessive unconstrained
    /// correction. Active bounds and floating-point stagnation are not success.
    NotConverged {
        /// Best accepted physical frequencies [Hz].
        frequencies_hz: Vec<f64>,
        /// Joint unexplained RMS in the channels' common physical units.
        residual_rms: f64,
        /// Candidate sets evaluated, including rejected trials.
        evaluations: usize,
        /// Remaining untruncated projected corrections [Hz].
        last_corrections_hz: Vec<f64>,
    },
}

impl From<HarmonicFitError> for ToneFrequencyError {
    fn from(error: HarmonicFitError) -> Self {
        Self::Fit(error)
    }
}

impl core::fmt::Display for ToneFrequencyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(f),
            Self::InvalidParameter { what } => write!(f, "tone frequency: {what}"),
            Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => write!(
                f,
                "tone frequencies are unidentifiable: reciprocal condition \
                 {reciprocal_condition:e} is below {minimum:e}"
            ),
            Self::NotConverged { evaluations, .. } => write!(
                f,
                "tone frequencies did not converge after {evaluations} candidate sets"
            ),
        }
    }
}

impl std::error::Error for ToneFrequencyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            _ => None,
        }
    }
}

fn invalid(what: &'static str) -> ToneFrequencyError {
    ToneFrequencyError::InvalidParameter { what }
}

fn unidentifiable(condition: f64, budget: ToneFrequencyBudget) -> ToneFrequencyError {
    ToneFrequencyError::Unidentifiable {
        reciprocal_condition: condition,
        minimum: budget.min_frequency_reciprocal_condition,
    }
}

struct Candidate {
    frequencies_hz: Vec<f64>,
    fits: Vec<HarmonicFit>,
    residual_rms: f64,
    residual_squared: f64,
    corrections_hz: Vec<f64>,
    reciprocal_condition: f64,
}

impl Candidate {
    fn not_converged(self, evaluations: usize) -> ToneFrequencyError {
        ToneFrequencyError::NotConverged {
            frequencies_hz: self.frequencies_hz,
            residual_rms: self.residual_rms,
            evaluations,
            last_corrections_hz: self.corrections_hz,
        }
    }
}

/// Refine independently seeded stationary tone frequencies and their amplitudes
/// together. Inharmonic partials and beating unisons retain their actual spacing.
/// Every seed has a positive, finite, closed search interval; these intervals
/// must be strictly separated and ordered and stay strictly below Nyquist.
///
/// Variable projection eliminates DC and all sinusoid amplitudes through the
/// existing [`ToneFitPlan`]. Analytic derivatives are projected outside that
/// shared design, then solved JOINTLY with Householder QR; independent scalar
/// frequency corrections would discard the cross terms of close tones.
/// Monotone backtracking and a quarter-cycle step cap bound each local trial.
/// This is the projected/Kaufman Gauss--Newton approximation; it omits the
/// residual-dependent part of the full variable-projection Jacobian.
///
/// This identifies a supplied local stationary spectrum. It neither discovers
/// the number of tones nor certifies uniqueness, noise uncertainty, damping,
/// time-varying frequency or arbitrary sub-bin resolution. Every requested tone
/// must supply numerically resolved frequency information. Omitted physics and
/// tones remain in the returned sampled residual, not a model-error certificate.
pub fn refine_tone_frequencies(
    samples: &[f64],
    timestep_s: f64,
    initial_frequencies_hz: &[f64],
    frequency_bounds_hz: &[[f64; 2]],
    budget: ToneFrequencyBudget,
) -> Result<ToneFrequencyFit, ToneFrequencyError> {
    let mut result = refine_tone_frequencies_channels(
        &[samples],
        timestep_s,
        initial_frequencies_hz,
        frequency_bounds_hz,
        budget,
    )?;
    Ok(ToneFrequencyFit {
        frequencies_hz: result.frequencies_hz,
        fit: result.fits.remove(0),
        evaluations: result.evaluations,
        last_corrections_hz: result.last_corrections_hz,
        frequency_reciprocal_condition: result.frequency_reciprocal_condition,
    })
}

/// Refine the common independent-frequency set of simultaneous, equal-length
/// channels with one sample clock. Each channel keeps independent DC/phasors;
/// the objective is the sum of squared residuals in COMMON physical units.
/// Original amplitude weights are retained, including under rotations of vector
/// components. Unrelated units or sensor-noise metrics need caller conversion.
/// Constant channels retain their DC but supply no frequency information.
///
/// One sinusoid QR serves every channel at each candidate. Projected derivative
/// and residual blocks are compressed and combined with QR, never normal
/// equations. Streaming channels costs O(E C (N K^2 + K^3)) work and
/// O(N K + K^2 + C K) numeric storage for N samples, K tones, C channels and
/// at most E candidate sets. This is an offline observation owner, without
/// mid-factorization cancellation or new cross-ISA certification claims.
pub fn refine_tone_frequencies_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_frequencies_hz: &[f64],
    frequency_bounds_hz: &[[f64; 2]],
    budget: ToneFrequencyBudget,
) -> Result<ToneFrequencyChannelsFit, ToneFrequencyError> {
    let tones = initial_frequencies_hz.len();
    if tones == 0 || frequency_bounds_hz.len() != tones {
        return Err(invalid(
            "each independent tone needs one frequency interval",
        ));
    }
    if !(timestep_s.is_finite() && timestep_s > 0.0)
        || !(budget.absolute_tolerance_hz.is_finite() && budget.absolute_tolerance_hz > 0.0)
        || budget.max_evaluations == 0
        || !(budget.min_frequency_reciprocal_condition.is_finite()
            && budget.min_frequency_reciprocal_condition > 0.0
            && budget.min_frequency_reciprocal_condition <= 1.0)
    {
        return Err(invalid(
            "invalid clock, frequency tolerance, work or condition limit",
        ));
    }
    for (k, (&seed, &[lower, upper])) in initial_frequencies_hz
        .iter()
        .zip(frequency_bounds_hz)
        .enumerate()
    {
        if !(lower.is_finite() && upper.is_finite() && lower > 0.0 && lower < upper)
            || !(seed.is_finite() && seed >= lower && seed <= upper)
            || !(lower * timestep_s > 0.0 && upper * timestep_s < 0.5)
            || (k > 0 && frequency_bounds_hz[k - 1][1] >= lower)
        {
            return Err(invalid(
                "frequency intervals must contain their seeds, be positive, separated, ordered and sub-Nyquist",
            ));
        }
    }
    let Some(samples) = channels.first() else {
        return Err(invalid("at least one simultaneous channel is required"));
    };
    let n = samples.len();
    let columns = tones
        .checked_mul(2)
        .and_then(|k| k.checked_add(1))
        .ok_or_else(|| invalid("tone design size overflows"))?;
    if n <= columns || channels.iter().any(|channel| channel.len() != n) {
        return Err(invalid(
            "equal-length channels need more than 1 + 2 * tone_count samples",
        ));
    }
    if (n - columns)
        .checked_mul(channels.len())
        .is_none_or(|rows| rows < tones)
    {
        return Err(invalid(
            "insufficient joint residual dimensions to identify all frequencies",
        ));
    }
    let duration = (n - 1) as f64 * timestep_s;
    if !(duration.is_finite() && duration > 0.0) {
        return Err(invalid("record duration must be finite and positive"));
    }
    if channels.iter().any(|c| c.iter().any(|x| !x.is_finite())) {
        return Err(HarmonicFitError::NonFinite {
            what: "input sample",
        }
        .into());
    }
    let width = tones + 1;
    let rows = (n - columns).max(width);
    let extra_bytes = refinement_bytes(n, tones, channels.len(), rows)?;
    let required_bytes = super::fit_workspace_bytes(n, tones)?
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
    let plan = ToneFitPlan::new(n, timestep_s, initial_frequencies_hz, plan_budget)?;
    let mut ratios: Vec<f64> = channels
        .iter()
        .map(|channel| {
            let lo = channel.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = channel.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let center = lo.midpoint(hi);
            channel
                .iter()
                .map(|x| (x - center).abs())
                .fold(0.0, f64::max)
        })
        .collect();
    let physical_scale = ratios.iter().copied().fold(0.0, f64::max);
    if physical_scale == 0.0 {
        return Err(unidentifiable(0.0, budget));
    }
    for ratio in &mut ratios {
        *ratio /= physical_scale;
    }
    let mut workspace = Workspace {
        scratch: vec![0.0; n],
        derivative: vec![0.0; n],
        augmented: vec![0.0; rows * width],
        folded: vec![0.0; 2 * width * width],
    };
    let mut current = evaluate(
        channels,
        &ratios,
        physical_scale,
        plan,
        duration,
        budget,
        &mut workspace,
    )?;
    let mut evaluations = 1;
    loop {
        let largest = current
            .corrections_hz
            .iter()
            .map(|x| x.abs())
            .fold(0.0, f64::max);
        if largest <= budget.absolute_tolerance_hz {
            return Ok(ToneFrequencyChannelsFit {
                frequencies_hz: current.frequencies_hz,
                fits: current.fits,
                evaluations,
                last_corrections_hz: current.corrections_hz,
                frequency_reciprocal_condition: current.reciprocal_condition,
            });
        }
        let mut fraction = (0.25 / duration / largest).min(1.0);
        loop {
            if evaluations == budget.max_evaluations {
                return Err(current.not_converged(evaluations));
            }
            let trial_frequencies: Vec<f64> = current
                .frequencies_hz
                .iter()
                .zip(&current.corrections_hz)
                .zip(frequency_bounds_hz)
                .map(|((&f, &d), &[lo, hi])| (f + fraction * d).clamp(lo, hi))
                .collect();
            if trial_frequencies == current.frequencies_hz {
                return Err(current.not_converged(evaluations));
            }
            evaluations += 1;
            let trial = ToneFitPlan::new(n, timestep_s, &trial_frequencies, plan_budget)
                .map_err(ToneFrequencyError::from)
                .and_then(|plan| {
                    evaluate(
                        channels,
                        &ratios,
                        physical_scale,
                        plan,
                        duration,
                        budget,
                        &mut workspace,
                    )
                });
            match trial {
                Ok(candidate) if candidate.residual_squared < current.residual_squared => {
                    current = candidate;
                    break;
                }
                Ok(_)
                | Err(ToneFrequencyError::Fit(HarmonicFitError::IllConditioned { .. }))
                | Err(ToneFrequencyError::Unidentifiable { .. }) => fraction *= 0.5,
                Err(error) => return Err(error),
            }
        }
    }
}

// Sum simultaneously live numeric storage conservatively, including both
// candidate reports and the previous/current small QR during a channel fold.
fn refinement_bytes(
    n: usize,
    k: usize,
    channels: usize,
    rows: usize,
) -> Result<usize, ToneFrequencyError> {
    let s = k + 1;
    let scalars = n
        .checked_mul(2)
        .and_then(|v| rows.checked_mul(s)?.checked_mul(2)?.checked_add(v))
        .and_then(|v| s.checked_mul(s)?.checked_mul(6)?.checked_add(v))
        .and_then(|v| s.checked_mul(3)?.checked_add(v))
        .and_then(|v| k.checked_mul(8)?.checked_add(v))
        .and_then(|v| {
            k.checked_mul(4)?
                .checked_add(1)?
                .checked_mul(channels)?
                .checked_add(v)
        })
        .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()));
    scalars
        .and_then(|v| {
            channels
                .checked_mul(2 * core::mem::size_of::<HarmonicFit>())?
                .checked_add(v)
        })
        .ok_or_else(|| invalid("frequency refinement workspace size overflows"))
}

struct Workspace {
    scratch: Vec<f64>,
    derivative: Vec<f64>,
    augmented: Vec<f64>,
    folded: Vec<f64>,
}

// D contains d(model)/d(f*T), after projection outside the fitted DC/sinusoid
// columns. For each channel, QR([D, residual]) compresses all cross terms.
// QR([R_previous; R_channel]) preserves the common-unit least-squares problem
// without squaring its condition number or storing C length-N derivative blocks.
#[allow(clippy::too_many_arguments)]
fn evaluate(
    channels: &[&[f64]],
    ratios: &[f64],
    physical_scale: f64,
    plan: ToneFitPlan,
    duration: f64,
    budget: ToneFrequencyBudget,
    work: &mut Workspace,
) -> Result<Candidate, ToneFrequencyError> {
    let k = plan.frequencies_hz.len();
    let columns = 2 * k + 1;
    let width = k + 1;
    let rows = work.augmented.len() / width;
    work.folded.fill(0.0);
    let mut factor: Option<Qr> = None;
    let mut fits = Vec::with_capacity(channels.len());
    let (mut residual_squared, mut joint_scale, mut raw_scale) = (0.0, 0.0_f64, 0.0_f64);
    for (&samples, &ratio) in channels.iter().zip(ratios) {
        let mut phasors = vec![C64::default(); k];
        let summary = plan.fit_into(samples, &mut work.scratch, &mut phasors)?;
        fits.push(HarmonicFit {
            dc: summary.dc,
            phasors,
            residual_rms: summary.residual_rms,
        });
        if ratio == 0.0 {
            continue;
        }
        for (j, &residual) in work.scratch[columns..].iter().enumerate() {
            let weighted = residual * ratio;
            work.augmented[j * width + k] = weighted;
            residual_squared = weighted.mul_add(weighted, residual_squared);
        }
        let mut channel_scale = 0.0_f64;
        for (tone, &frequency) in plan.frequencies_hz.iter().enumerate() {
            let cycles = frequency * plan.timestep_s();
            let a = work.scratch[2 * tone + 1];
            let b = work.scratch[2 * tone + 2];
            for (j, out) in work.derivative.iter_mut().enumerate() {
                let centered_time = j as f64 / (samples.len() - 1) as f64 - 0.5;
                let (s, c) = det::sin_cos(core::f64::consts::TAU * (cycles * j as f64));
                *out = core::f64::consts::TAU * centered_time * (-a).mul_add(s, b * c);
                raw_scale = raw_scale.max(out.abs() * ratio);
            }
            plan.design.factor.apply_qt(&mut work.derivative);
            for (j, &derivative) in work.derivative[columns..].iter().enumerate() {
                let weighted = derivative * ratio;
                if !weighted.is_finite() {
                    return Err(HarmonicFitError::NonFinite {
                        what: "projected tone derivative",
                    }
                    .into());
                }
                work.augmented[j * width + tone] = weighted;
                channel_scale = channel_scale.max(weighted.abs());
            }
        }
        if channel_scale > joint_scale {
            let change = joint_scale / channel_scale;
            for row in work.folded[..width * width].chunks_exact_mut(width) {
                for value in &mut row[..k] {
                    *value *= change;
                }
            }
            joint_scale = channel_scale;
        }
        if joint_scale == 0.0 {
            continue;
        }
        for row in work.augmented.chunks_exact_mut(width) {
            for value in &mut row[..k] {
                *value /= joint_scale;
            }
        }
        let local = qr(&work.augmented, rows, width);
        for i in 0..width {
            for j in 0..width {
                work.folded[(i + width) * width + j] = local.r(i, j);
            }
        }
        drop(local);
        let combined = qr(&work.folded, 2 * width, width);
        for i in 0..width {
            for j in 0..width {
                work.folded[i * width + j] = combined.r(i, j);
            }
        }
        factor = Some(combined);
    }
    if joint_scale <= 64.0 * f64::EPSILON * raw_scale {
        return Err(unidentifiable(0.0, budget));
    }
    let factor = factor.ok_or_else(|| unidentifiable(0.0, budget))?;
    let condition = super::triangular_reciprocal_condition(&factor, k);
    if !condition.is_finite() || condition < budget.min_frequency_reciprocal_condition {
        return Err(unidentifiable(condition, budget));
    }
    let mut correction: Vec<f64> = (0..k).map(|i| factor.r(i, k)).collect();
    for i in (0..k).rev() {
        let mut value = correction[i];
        for (j, &x) in correction.iter().enumerate().skip(i + 1) {
            value = (-factor.r(i, j)).mul_add(x, value);
        }
        correction[i] = value / factor.r(i, i);
    }
    for value in &mut correction {
        *value = (*value / joint_scale) / duration;
    }
    let residual_rms =
        det::sqrt(residual_squared / work.scratch.len() as f64 / channels.len() as f64)
            * physical_scale;
    if correction.iter().any(|x| !x.is_finite())
        || !residual_squared.is_finite()
        || !residual_rms.is_finite()
    {
        return Err(HarmonicFitError::NonFinite {
            what: "joint frequency correction",
        }
        .into());
    }
    Ok(Candidate {
        frequencies_hz: plan.frequencies_hz,
        fits,
        residual_rms,
        residual_squared,
        corrections_hz: correction,
        reciprocal_condition: condition,
    })
}
