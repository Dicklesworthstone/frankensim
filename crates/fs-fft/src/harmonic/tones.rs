//! Joint local identification of independent stationary tone frequencies.

use super::refinement::{self, Controls, Domain};
use super::{HarmonicFit, HarmonicFitBudget, HarmonicFitError};

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

impl From<refinement::Error> for ToneFrequencyError {
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
                corrections,
            } => Self::NotConverged {
                frequencies_hz: parameters,
                residual_rms,
                evaluations,
                last_corrections_hz: corrections,
            },
        }
    }
}

/// Refine independently seeded stationary tone frequencies and their amplitudes
/// together. Inharmonic partials and beating unisons retain their actual spacing.
/// Every seed has a positive, finite, closed search interval; these intervals
/// must be strictly separated and ordered and stay strictly below Nyquist.
///
/// Variable projection eliminates DC and all sinusoid amplitudes through the
/// existing [`super::ToneFitPlan`]. Analytic derivatives are projected outside that
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
    let result = refinement::refine_channels(
        channels,
        timestep_s,
        Domain::Stationary(initial_frequencies_hz, frequency_bounds_hz),
        Controls {
            tolerances: [budget.absolute_tolerance_hz; 2],
            max_evaluations: budget.max_evaluations,
            min_reciprocal_condition: budget.min_frequency_reciprocal_condition,
            fit: budget.fit,
        },
    )?;
    Ok(ToneFrequencyChannelsFit {
        frequencies_hz: result.parameters,
        fits: result.fits,
        evaluations: result.evaluations,
        last_corrections_hz: result.corrections,
        frequency_reciprocal_condition: result.reciprocal_condition,
    })
}
