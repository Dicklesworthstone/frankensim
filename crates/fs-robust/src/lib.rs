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

pub use fs_evidence::{AdmittedColor, Color, ColorPayloadError, ColorRank};
use fs_evidence::validate_color_payload;

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
    /// An input color's payload is structurally malformed (bead 6pf9): a
    /// headline can never be derived over structural garbage.
    MalformedInputColor {
        /// The offending design.
        design: String,
        /// The exact structural defect.
        error: ColorPayloadError,
    },
    /// A declared input color is not covered by an admitted counterpart, so
    /// no positive (admitted) headline exists for this objective.
    UnadmittedInput {
        /// The offending design.
        design: String,
    },
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
    #[allow(clippy::cast_precision_loss)]
    let boundary_weight = boundary_rank as f64 - n * alpha;
    let (at_or_below, above) = sorted.split_at(boundary_rank);
    let boundary = *at_or_below.last().ok_or(RobustError::EmptySamples)?;
    weighted_mean(
        core::iter::once((boundary, boundary_weight))
            .chain(above.iter().copied().map(|value| (value, 1.0))),
    )
    .ok_or(RobustError::EmptySamples)
}

/// The weakest (lowest-rank) DECLARED color among the inputs — the
/// reporting rule over unverified candidates. Rank ties break by canonical
/// payload bytes, so the result is PERMUTATION-INVARIANT: reordering equal
/// inputs can never change which payload the report carries (bead 6pf9).
/// This is a declaration-level API; positive-evidence reporting goes
/// through [`admitted_headline_for`].
#[must_use]
pub fn weakest_color(colors: &[Color]) -> Option<Color> {
    colors
        .iter()
        .min_by_key(|c| (c.rank(), c.canonical_bytes()))
        .cloned()
}

/// The weakest ADMITTED color among positive inputs, permutation-invariant:
/// rank ties break by canonical payload bytes, then by admission-receipt
/// node hash, never by input order (bead 6pf9).
#[must_use]
pub fn weakest_admitted_color(inputs: &[AdmittedColor]) -> Option<AdmittedColor> {
    inputs
        .iter()
        .min_by_key(|a| {
            (
                a.rank(),
                a.admitted_color().canonical_bytes(),
                *a.receipt().node_hash().as_bytes(),
            )
        })
        .cloned()
}

/// The admitted headline for one objective: the weakest admitted input,
/// available ONLY when every declared input color is covered (count-aware,
/// canonical bytes) by an admitted counterpart. An Estimated declared input
/// can never be covered — [`AdmittedColor`] is always positive — so an
/// objective with estimated ingredients keeps its declared-only headline,
/// exactly as the no-laundering law requires.
///
/// # Errors
/// [`RobustError::UncoloredObjective`] with no declared inputs,
/// [`RobustError::MalformedInputColor`] on structural garbage, and
/// [`RobustError::UnadmittedInput`] when any declared input lacks an
/// admitted counterpart.
pub fn admitted_headline_for(
    objective: &ColoredObjective,
    admitted: &[AdmittedColor],
) -> Result<AdmittedColor, RobustError> {
    if objective.input_colors.is_empty() {
        return Err(RobustError::UncoloredObjective {
            design: objective.design.clone(),
        });
    }
    for color in &objective.input_colors {
        validate_color_payload(color).map_err(|error| RobustError::MalformedInputColor {
            design: objective.design.clone(),
            error,
        })?;
    }
    // Count-aware coverage: each declared input consumes one admitted
    // counterpart with identical canonical bytes. Only CONSUMED counterparts
    // enter the headline — surplus admitted values a caller happens to hold
    // must not influence this objective's report.
    let mut available: Vec<Option<&AdmittedColor>> = admitted.iter().map(Some).collect();
    let mut consumed: Vec<AdmittedColor> = Vec::with_capacity(objective.input_colors.len());
    for color in &objective.input_colors {
        let declared_bytes = color.canonical_bytes();
        // Among ALL available counterparts with matching canonical bytes,
        // consume the CANONICAL one — minimal under the SAME total key
        // `weakest_admitted_color` uses — never the first in input order. A
        // same-color surplus (two admitted colors covering one declared input
        // but carrying DIFFERENT receipts) would otherwise make WHICH receipt
        // enters the headline depend on the order of `admitted`, breaking the
        // permutation-invariance the report promises.
        let chosen = available
            .iter()
            .enumerate()
            .filter(|(_, slot)| {
                slot.is_some_and(|a| a.admitted_color().canonical_bytes() == declared_bytes)
            })
            .min_by_key(|(_, slot)| {
                let a = slot.expect("filtered to occupied");
                (
                    a.rank(),
                    a.admitted_color().canonical_bytes(),
                    *a.receipt().node_hash().as_bytes(),
                )
            })
            .map(|(index, _)| index);
        match chosen {
            Some(index) => consumed.push(
                available[index]
                    .take()
                    .expect("chosen slot is occupied")
                    .clone(),
            ),
            None => {
                return Err(RobustError::UnadmittedInput {
                    design: objective.design.clone(),
                });
            }
        }
    }
    weakest_admitted_color(&consumed).ok_or(RobustError::UnadmittedInput {
        design: objective.design.clone(),
    })
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
        let mut ordered = self.cost_samples.clone();
        ordered.sort_by(f64::total_cmp);
        weighted_mean(ordered.into_iter().map(|value| (value, 1.0)))
            .ok_or(RobustError::EmptySamples)
    }

    /// The headline color = the WEAKEST input color. An objective with no
    /// declared input colors is un-colored (contract violation).
    ///
    /// # Errors
    /// [`RobustError::UncoloredObjective`] if no input colors are declared.
    pub fn headline_color(&self) -> Result<Color, RobustError> {
        for color in &self.input_colors {
            validate_color_payload(color).map_err(|error| RobustError::MalformedInputColor {
                design: self.design.clone(),
                error,
            })?;
        }
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

/// A robust-optimization report whose headline is ADMITTED scientific
/// evidence (bead 6pf9): holding this value means every candidate's declared
/// inputs were covered by admitted counterparts and the winner's headline is
/// the weakest admitted input, permutation-invariant.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedRobustReport {
    /// The winning design.
    pub design: String,
    /// Its robust value (CVaR) — the objective minimized.
    pub robust_value: f64,
    /// Its nominal value (mean) — for comparison.
    pub nominal_value: f64,
    /// The admitted headline: the weakest admitted input of the winner.
    pub headline: AdmittedColor,
}

/// [`robust_optimum`] under the admitted-evidence contract: EVERY candidate's
/// declared inputs must be covered by its paired admitted colors before the
/// optimization is allowed to claim a positive headline — an optimization run
/// publishing admitted evidence must be wholly admitted, not admitted only at
/// the winner (fail closed).
///
/// # Errors
/// Everything [`robust_optimum`] refuses, plus
/// [`RobustError::MalformedInputColor`] and [`RobustError::UnadmittedInput`]
/// from any candidate.
pub fn robust_optimum_admitted(
    candidates: &[(ColoredObjective, &[AdmittedColor])],
    alpha: f64,
) -> Result<AdmittedRobustReport, RobustError> {
    if candidates.is_empty() {
        return Err(RobustError::NoCandidates);
    }
    let mut headlines = Vec::with_capacity(candidates.len());
    for (objective, admitted) in candidates {
        headlines.push(admitted_headline_for(objective, admitted)?);
    }
    let mut best: Option<(usize, f64)> = None;
    for (index, (objective, _)) in candidates.iter().enumerate() {
        let rv = objective.robust_value(alpha)?;
        if best.is_none_or(|(_, best_rv)| rv < best_rv) {
            best = Some((index, rv));
        }
    }
    let (index, robust_value) = best.ok_or(RobustError::NoCandidates)?;
    let winner = &candidates[index].0;
    Ok(AdmittedRobustReport {
        design: winner.design.clone(),
        robust_value,
        nominal_value: winner.nominal_value()?,
        headline: headlines.swap_remove(index),
    })
}

/// The Proposal-F kill-criterion test: is the robust optimum DOMINATED by
/// nominal-optimum-plus-standard-safety-factor on realized cost (at equal
/// achieved safety)? If robust designs are consistently dominated, the
/// ambiguity sets are miscalibrated.
///
/// # Errors
/// [`RobustError::BadSample`] if either realized cost is non-finite.
pub fn dominated_by_nominal(
    robust_realized_cost: f64,
    nominal_plus_safety_realized_cost: f64,
) -> Result<bool, RobustError> {
    reject_non_finite(&[robust_realized_cost, nominal_plus_safety_realized_cost])?;
    Ok(nominal_plus_safety_realized_cost < robust_realized_cost)
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
/// exceeds capacity). Intensities are returned in deterministic ascending
/// order, so the empirical CDF is structurally monotone. The `color` is the
/// honest confidence band.
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
    let mut ordered_intensities = intensities.to_vec();
    ordered_intensities.sort_by(f64::total_cmp);
    let curve = ordered_intensities
        .into_iter()
        .map(|intensity| {
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

fn weighted_mean(values: impl IntoIterator<Item = (f64, f64)>) -> Option<f64> {
    let mut mean = 0.0_f64;
    let mut total_weight = 0.0_f64;
    let mut present = false;
    for (value, weight) in values {
        if weight <= 0.0 {
            continue;
        }
        let combined_weight = total_weight + weight;
        if present {
            let incoming_share = weight / combined_weight;
            mean = if mean.is_sign_positive() == value.is_sign_positive() {
                mean + (value - mean) * incoming_share
            } else {
                mean * (total_weight / combined_weight) + value * incoming_share
            };
        } else {
            mean = value;
            present = true;
        }
        total_weight = combined_weight;
    }
    present.then_some(mean)
}
