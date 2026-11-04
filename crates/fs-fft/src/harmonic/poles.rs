//! Local identification of physical oscillatory ringdown poles.

use super::refinement::{self, Controls, Domain};
use super::{DampedTonePole, HarmonicFit, HarmonicFitBudget, HarmonicFitError};

/// Closed local search rectangle for one oscillatory decay pole.
/// Rectangles must be strictly disjoint in at least one coordinate. Distinct
/// decay intervals can therefore share a frequency interval for double decay.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedToneBounds {
    /// Positive, strictly sub-Nyquist frequency interval [Hz], with lower < upper.
    pub frequency_hz: [f64; 2],
    /// Nonnegative amplitude-decay interval [1/s], with lower < upper.
    pub decay_rate_per_s: [f64; 2],
}

/// Work, storage and numerical resolution admitted for joint pole refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedToneRefinementBudget {
    /// Every untruncated frequency correction must meet this positive limit [Hz].
    pub absolute_frequency_tolerance_hz: f64,
    /// Every untruncated amplitude-decay correction must meet this limit [1/s].
    pub absolute_decay_tolerance_per_s: f64,
    /// Maximum candidate pole sets, including rejected line-search trials.
    pub max_evaluations: usize,
    /// Minimum reciprocal 1-norm condition of BOTH the joint profiled Jacobian
    /// and its projected physical sensitivity, using coordinates `f * duration`
    /// and `decay * duration` without independent column normalization. Residual
    /// curvature cannot identify an absent pole; projected sensitivities must
    /// also exceed a roundoff floor in globally normalized signal units.
    /// This is not uncertainty.
    pub min_pole_reciprocal_condition: f64,
    /// Linear-fit conditioning and total conservative numeric storage cap,
    /// including all refinement scratch and both candidates' channel fits.
    pub fit: HarmonicFitBudget,
}

/// Signed, untruncated local correction, not a confidence interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedToneCorrection {
    /// Frequency correction [Hz].
    pub frequency_hz: f64,
    /// Amplitude-decay correction [1/s]; may be negative.
    pub decay_rate_per_s: f64,
}

/// Jointly refined poles and first-sample phasors for one physical record.
#[derive(Debug, Clone, PartialEq)]
pub struct DampedTonePoleFit {
    /// Physical damped frequencies and amplitude decay rates, in seed order.
    pub poles: Vec<DampedTonePole>,
    /// Fitted DC, peak `exp(-i omega t)` phasors and unexplained sampled RMS.
    pub fit: HarmonicFit,
    /// Candidate pole sets evaluated.
    pub evaluations: usize,
    /// Remaining untruncated corrections, in the same pole order.
    pub last_corrections: Vec<DampedToneCorrection>,
    /// Minimum reciprocal condition of the full and projected pole Jacobians.
    pub pole_reciprocal_condition: f64,
}

/// One physical pole set fitted to simultaneous common-unit channels.
#[derive(Debug, Clone, PartialEq)]
pub struct DampedTonePoleChannelsFit {
    /// Shared physical poles in the supplied seed order.
    pub poles: Vec<DampedTonePole>,
    /// Per-channel DC, first-sample peak phasors and unexplained RMS.
    pub fits: Vec<HarmonicFit>,
    /// Candidate pole sets evaluated; each candidate shares one amplitude QR.
    pub evaluations: usize,
    /// Remaining untruncated corrections, in the same pole order.
    pub last_corrections: Vec<DampedToneCorrection>,
    /// Minimum reciprocal condition of the full and projected pole Jacobians.
    pub pole_reciprocal_condition: f64,
}

/// Refusal to infer unresolved poles or publish an unconverged ringdown fit.
#[derive(Debug, Clone, PartialEq)]
pub enum DampedToneRefinementError {
    /// The shared joint-fit owner refused data, conditioning or storage.
    Fit(HarmonicFitError),
    /// Invalid sampling, local pole domains, channel shapes or controls.
    InvalidParameter {
        /// The failed admission condition.
        what: &'static str,
    },
    /// Pole derivatives have insufficient numerical rank or signal.
    Unidentifiable {
        /// The refused full/projected Jacobian's reciprocal 1-norm condition.
        reciprocal_condition: f64,
        /// Requested minimum reciprocal condition.
        minimum: f64,
    },
    /// Work exhaustion, an active bound or stagnation left an excessive
    /// untruncated correction. No partial phasors are published.
    NotConverged {
        /// Best accepted physical poles, in supplied seed order.
        poles: Vec<DampedTonePole>,
        /// Joint unexplained RMS in the channels' common physical units.
        residual_rms: f64,
        /// Candidate pole sets evaluated, including rejected trials.
        evaluations: usize,
        /// Remaining untruncated corrections at that candidate.
        last_corrections: Vec<DampedToneCorrection>,
    },
}

impl core::fmt::Display for DampedToneRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(f),
            Self::InvalidParameter { what } => write!(f, "damped tone poles: {what}"),
            Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => write!(
                f,
                "damped tone poles are unidentifiable: reciprocal condition \
                 {reciprocal_condition:e} is below {minimum:e}"
            ),
            Self::NotConverged { evaluations, .. } => write!(
                f,
                "damped tone poles did not converge after {evaluations} candidate sets"
            ),
        }
    }
}

impl std::error::Error for DampedToneRefinementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            _ => None,
        }
    }
}

fn poles(parameters: &[f64]) -> Vec<DampedTonePole> {
    parameters
        .chunks_exact(2)
        .map(|p| DampedTonePole {
            frequency_hz: p[0],
            decay_rate_per_s: p[1],
        })
        .collect()
}

fn corrections(parameters: &[f64]) -> Vec<DampedToneCorrection> {
    parameters
        .chunks_exact(2)
        .map(|p| DampedToneCorrection {
            frequency_hz: p[0],
            decay_rate_per_s: p[1],
        })
        .collect()
}

impl From<refinement::Error> for DampedToneRefinementError {
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
                poles: poles(&parameters),
                residual_rms,
                evaluations,
                last_corrections: corrections(&delta),
            },
        }
    }
}

fn invalid(what: &'static str) -> DampedToneRefinementError {
    DampedToneRefinementError::InvalidParameter { what }
}

/// Refine a locally seeded oscillatory ringdown without imposing harmonic
/// spacing or a common decay rate. The exact sampled model is
/// `dc + sum exp(-decay*t) * (a*cos(TAU*f*t) + b*sin(TAU*f*t))`.
/// Frequencies are the damped oscillation frequencies; decays act on amplitude.
/// All phasors refer to the first input sample and preserve supplied pole order.
///
/// Amplitude/DC elimination uses [`super::DampedToneFitPlan`]. The shared
/// bounded QR optimizer uses the COMPLETE variable-projection derivative:
/// for `r=y-A*c`, `-dr/dp = P_perp*D_p*c + (A^+)^T*D_p^T*r`. Its second,
/// residual-dependent term is retained using the original QR, including away
/// from an exact fit. See O'Leary and Rust, *Variable Projection for Nonlinear
/// Least Squares Problems*, §2.3:
/// <https://www.cs.umd.edu/users/oleary/software/varpro.pdf>.
///
/// This is a bounded local estimate with explicit numerical conditioning and
/// convergence criteria. It does not discover modes, certify noisy-record
/// separability, estimate uncertainty, fit growth or chirps, prevent aliasing,
/// or certify the underlying material model. Unmodeled motion stays in the
/// residual and may bias the fitted poles. An active bound is not convergence.
pub fn refine_damped_tone_poles(
    samples: &[f64],
    timestep_s: f64,
    initial_poles: &[DampedTonePole],
    bounds: &[DampedToneBounds],
    budget: DampedToneRefinementBudget,
) -> Result<DampedTonePoleFit, DampedToneRefinementError> {
    let mut result =
        refine_damped_tone_poles_channels(&[samples], timestep_s, initial_poles, bounds, budget)?;
    Ok(DampedTonePoleFit {
        poles: result.poles,
        fit: result.fits.remove(0),
        evaluations: result.evaluations,
        last_corrections: result.last_corrections,
        pole_reciprocal_condition: result.pole_reciprocal_condition,
    })
}

/// Jointly refine shared physical poles from equal-length simultaneous channels
/// with one sample clock. Each channel keeps its own DC and peak phasors. The
/// objective retains original COMMON-unit amplitude weights, including under
/// orthogonal changes of vector components; caller conversion is needed for
/// unrelated units or a different noise metric. Constant channels supply no
/// pole information. Channel or pole reordering is never implicit.
///
/// One amplitude QR serves all channels of each candidate. Full derivative and
/// residual blocks are folded by QR without normal equations or a channel-stack
/// matrix. With N samples, K poles, C channels and E candidate sets, work is
/// O(E C (N K² + K³)) and numeric storage O(N K + K² + C K). This is an offline
/// observation owner without mid-factorization cancellation or throughput claims.
pub fn refine_damped_tone_poles_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_poles: &[DampedTonePole],
    bounds: &[DampedToneBounds],
    budget: DampedToneRefinementBudget,
) -> Result<DampedTonePoleChannelsFit, DampedToneRefinementError> {
    if initial_poles.is_empty() || initial_poles.len() != bounds.len() {
        return Err(invalid("each damped pole needs one search rectangle"));
    }
    if !(timestep_s.is_finite() && timestep_s > 0.0)
        || ![
            budget.absolute_frequency_tolerance_hz,
            budget.absolute_decay_tolerance_per_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        || budget.max_evaluations == 0
        || !(budget.min_pole_reciprocal_condition.is_finite()
            && budget.min_pole_reciprocal_condition > 0.0
            && budget.min_pole_reciprocal_condition <= 1.0)
    {
        return Err(invalid("invalid clock, pole tolerance, work or condition limit"));
    }
    for (k, (seed, bound)) in initial_poles.iter().zip(bounds).enumerate() {
        let [fl, fh] = bound.frequency_hz;
        let [dl, dh] = bound.decay_rate_per_s;
        if ![fl, fh, dl, dh, seed.frequency_hz, seed.decay_rate_per_s]
            .into_iter()
            .all(f64::is_finite)
            || !(fl > 0.0 && fl < fh && fl * timestep_s > 0.0 && fh * timestep_s < 0.5)
            || !(dl >= 0.0 && dl < dh)
            || !(seed.frequency_hz >= fl && seed.frequency_hz <= fh)
            || !(seed.decay_rate_per_s >= dl && seed.decay_rate_per_s <= dh)
        {
            return Err(invalid(
                "pole rectangles must contain their seeds, have positive sub-Nyquist frequencies and nonnegative decays",
            ));
        }
        if bounds[..k].iter().any(|other| {
            !(fh < other.frequency_hz[0]
                || other.frequency_hz[1] < fl
                || dh < other.decay_rate_per_s[0]
                || other.decay_rate_per_s[1] < dl)
        }) {
            return Err(invalid("pole rectangles must be strictly disjoint in frequency or decay"));
        }
    }
    let result = refinement::refine_channels(
        channels,
        timestep_s,
        Domain::Damped(initial_poles, bounds),
        Controls {
            tolerances: [
                budget.absolute_frequency_tolerance_hz,
                budget.absolute_decay_tolerance_per_s,
            ],
            max_evaluations: budget.max_evaluations,
            min_reciprocal_condition: budget.min_pole_reciprocal_condition,
            fit: budget.fit,
        },
    )?;
    Ok(DampedTonePoleChannelsFit {
        poles: poles(&result.parameters),
        fits: result.fits,
        evaluations: result.evaluations,
        last_corrections: corrections(&result.corrections),
        pole_reciprocal_condition: result.reciprocal_condition,
    })
}
