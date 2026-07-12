//! fs-robust — objective epistemics (plan addendum, Proposal F). Layer: L4.
//!
//! Proposal 3 gives the three colors to the physics; Proposal F gives them to
//! the GOAL ITSELF. In "cheapest earthquake-resistant building" the structural
//! solve is the easy, well-posed half; what determines the answer is the
//! ground-motion ensemble, the cost spreadsheet, and the loss function. If
//! those live OUTSIDE the system as fixed inputs, the optimizer optimizes a
//! FICTION with certified precision — the Proposal 3 failure mode relocated to
//! the objective, where it is worse because nobody is watching.
//!
//! So an objective is a functional over `(design, environment)` whose inputs
//! carry COLORS ([`fs_evidence::Color`]), and:
//!
//! - the default deliverable is a ROBUST optimum — [`cvar`] (Conditional
//!   Value at Risk) over the cost distribution, not the nominal mean;
//! - the headline number's color is the color of its WEAKEST input
//!   ([`ColoredObjective::headline_color`]) — a verified solve under an
//!   estimated hazard is an ESTIMATED answer, and the report says so;
//! - [`robust_optimum`] enforces the amended optimization contract: NO
//!   optimization may run against an un-colored objective;
//! - the seismic deliverable is a colored [`fragility_curve`], not a binary
//!   "safe".
//!
//! Deterministic; sample paths are supplied by the caller (common random
//! numbers live in fs-scenario). This crate is the coloring + risk algebra.

pub use fs_evidence::{Color, ColorRank};

/// A structured objective-epistemics failure.
#[derive(Debug, Clone, PartialEq)]
pub enum RobustError {
    /// No cost samples were supplied.
    EmptySamples,
    /// The CVaR confidence level is not in `(0, 1)`.
    BadAlpha {
        /// The offending value.
        alpha: f64,
    },
    /// A supplied risk or fragility value was not finite.
    BadSample {
        /// The offending value.
        value: f64,
    },
    /// An objective declares no input colors — the optimization contract
    /// forbids optimizing it (it would optimize a fiction).
    UncoloredObjective {
        /// The offending design.
        design: String,
    },
    /// No candidate designs.
    NoCandidates,
}

/// Conditional Value at Risk at confidence `alpha`: the expected loss in the
/// worst `(1 − alpha)` tail of `samples` (costs/losses, higher = worse). More
/// robust than the mean because it weights the tail.
///
/// # Errors
/// [`RobustError::EmptySamples`] / [`RobustError::BadAlpha`] /
/// [`RobustError::BadSample`].
pub fn cvar(samples: &[f64], alpha: f64) -> Result<f64, RobustError> {
    if samples.is_empty() {
        return Err(RobustError::EmptySamples);
    }
    if !(alpha > 0.0 && alpha < 1.0) {
        return Err(RobustError::BadAlpha { alpha });
    }
    reject_non_finite(samples)?;
    let mut sorted: Vec<f64> = samples.to_vec();
    sorted.sort_by(f64::total_cmp);
    // Standard finite-sample empirical CVaR: when n*alpha is not an
    // integer, the boundary order statistic contributes only the fractional
    // mass needed to make the tail measure exactly n*(1-alpha). Giving every
    // one of ceil(n*(1-alpha)) samples equal weight dilutes the upper tail and
    // is anti-conservative.
    #[allow(clippy::cast_precision_loss)]
    let n = sorted.len() as f64;
    #[allow(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss
    )]
    let boundary_rank = ((n * alpha).ceil() as usize).clamp(1, sorted.len());
    let tail_mass = n * (1.0 - alpha);
    #[allow(clippy::cast_precision_loss)]
    let boundary_weight = boundary_rank as f64 - n * alpha;
    let (at_or_below, above) = sorted.split_at(boundary_rank);
    let boundary = *at_or_below.last().ok_or(RobustError::EmptySamples)?;
    let strictly_above = above.iter().sum::<f64>();
    Ok((boundary_weight * boundary + strictly_above) / tail_mass)
}

/// The weakest (lowest-rank) color among the inputs — the reporting rule.
#[must_use]
pub fn weakest_color(colors: &[Color]) -> Option<Color> {
    colors.iter().min_by_key(|c| c.rank()).cloned()
}

/// A colored objective for one design: its realized-cost samples under the
/// environment ensemble, plus the colors of the objective's ingredients
/// (hazard model, cost model, loss function, physics solve…).
#[derive(Debug, Clone, PartialEq)]
pub struct ColoredObjective {
    /// The design id.
    pub design: String,
    /// Realized-cost samples (common random numbers across designs).
    pub cost_samples: Vec<f64>,
    /// The colors of the objective's inputs.
    pub input_colors: Vec<Color>,
}

impl ColoredObjective {
    /// A colored objective.
    #[must_use]
    pub fn new(
        design: impl Into<String>,
        cost_samples: Vec<f64>,
        input_colors: Vec<Color>,
    ) -> ColoredObjective {
        ColoredObjective {
            design: design.into(),
            cost_samples,
            input_colors,
        }
    }

    /// The robust value: CVaR of the cost distribution at `alpha`.
    ///
    /// # Errors
    /// [`RobustError`] on bad samples/alpha.
    pub fn robust_value(&self, alpha: f64) -> Result<f64, RobustError> {
        cvar(&self.cost_samples, alpha)
    }

    /// The nominal value: the mean cost (what naive optimization targets).
    ///
    /// # Errors
    /// [`RobustError::EmptySamples`] if there are no samples;
    /// [`RobustError::BadSample`] if any sample is non-finite.
    pub fn nominal_value(&self) -> Result<f64, RobustError> {
        if self.cost_samples.is_empty() {
            return Err(RobustError::EmptySamples);
        }
        reject_non_finite(&self.cost_samples)?;
        Ok(self.cost_samples.iter().sum::<f64>() / self.cost_samples.len() as f64)
    }

    /// The headline color = the WEAKEST input color. An objective with no
    /// declared input colors is un-colored (contract violation).
    ///
    /// # Errors
    /// [`RobustError::UncoloredObjective`] if no input colors are declared.
    pub fn headline_color(&self) -> Result<Color, RobustError> {
        weakest_color(&self.input_colors).ok_or_else(|| RobustError::UncoloredObjective {
            design: self.design.clone(),
        })
    }
}

/// A robust-optimization report.
#[derive(Debug, Clone, PartialEq)]
pub struct RobustReport {
    /// The winning design.
    pub design: String,
    /// Its robust value (CVaR) — the objective minimized.
    pub robust_value: f64,
    /// Its nominal value (mean) — for comparison.
    pub nominal_value: f64,
    /// The headline color = the weakest input color of the winner.
    pub headline_color: Color,
}

/// Choose the robust optimum: the design minimizing CVaR at `alpha`. Enforces
/// the amended optimization contract — every candidate MUST be colored, else
/// the optimization is refused (it would optimize a fiction).
///
/// # Errors
/// [`RobustError::NoCandidates`], [`RobustError::UncoloredObjective`] (contract),
/// or a sample/alpha error.
pub fn robust_optimum(
    candidates: &[ColoredObjective],
    alpha: f64,
) -> Result<RobustReport, RobustError> {
    if candidates.is_empty() {
        return Err(RobustError::NoCandidates);
    }
    // contract: no optimization against an un-colored objective.
    for c in candidates {
        if c.input_colors.is_empty() {
            return Err(RobustError::UncoloredObjective {
                design: c.design.clone(),
            });
        }
    }
    let mut candidates = candidates.iter();
    let first = candidates.next().ok_or(RobustError::NoCandidates)?;
    let mut best = (first, first.robust_value(alpha)?);
    for c in candidates {
        let rv = c.robust_value(alpha)?;
        if rv < best.1 {
            best = (c, rv);
        }
    }
    let (winner, robust_value) = best;
    Ok(RobustReport {
        design: winner.design.clone(),
        robust_value,
        nominal_value: winner.nominal_value()?,
        headline_color: winner.headline_color()?,
    })
}

/// The Proposal-F kill-criterion test: is the robust optimum DOMINATED by
/// nominal-optimum-plus-standard-safety-factor on realized cost (at equal
/// achieved safety)? If robust designs are consistently dominated, the
/// ambiguity sets are miscalibrated.
#[must_use]
pub fn dominated_by_nominal(
    robust_realized_cost: f64,
    nominal_plus_safety_realized_cost: f64,
) -> bool {
    nominal_plus_safety_realized_cost < robust_realized_cost
}

/// A fragility-curve point: the probability of failure at a hazard intensity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragilityPoint {
    /// The hazard intensity.
    pub intensity: f64,
    /// `P(failure)` = fraction of capacity samples below the intensity.
    pub prob_failure: f64,
}

/// A fragility curve with its epistemic color band.
#[derive(Debug, Clone, PartialEq)]
pub struct ColoredFragility {
    /// The curve points, in intensity order.
    pub curve: Vec<FragilityPoint>,
    /// The color of the curve (the weakest input color — hazard + capacity).
    pub color: Color,
}

/// Build a fragility curve: at each intensity, `P(failure)` is the fraction of
/// `capacity_samples` that fall below it (the structure fails when demand
/// exceeds capacity). The `color` is the honest confidence band.
///
/// # Errors
/// [`RobustError::EmptySamples`] if there are no capacity samples;
/// [`RobustError::BadSample`] if any capacity sample or intensity is non-finite.
pub fn fragility_curve(
    capacity_samples: &[f64],
    intensities: &[f64],
    color: Color,
) -> Result<ColoredFragility, RobustError> {
    if capacity_samples.is_empty() {
        return Err(RobustError::EmptySamples);
    }
    reject_non_finite(capacity_samples)?;
    reject_non_finite(intensities)?;
    let n = capacity_samples.len() as f64;
    let curve = intensities
        .iter()
        .map(|&intensity| {
            let failures = capacity_samples.iter().filter(|&&c| c < intensity).count();
            FragilityPoint {
                intensity,
                prob_failure: failures as f64 / n,
            }
        })
        .collect();
    Ok(ColoredFragility { curve, color })
}

fn reject_non_finite(values: &[f64]) -> Result<(), RobustError> {
    if let Some(&value) = values.iter().find(|value| !value.is_finite()) {
        return Err(RobustError::BadSample { value });
    }
    Ok(())
}
