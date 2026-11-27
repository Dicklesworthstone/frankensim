//! One-step expected value of sample information (EVSI) for finite hypotheses.
//! Payoffs are utilities to maximize; a loss is supplied as its negative.
//! Each test outcome selects a new action. Its contribution is the probability
//! of that outcome times the advantage over keeping the prior-optimal action.

/// Version of the finite decision algebra, admission rules and tie policy.
pub const FINITE_EVSI_VERSION: u32 = 1;
/// Largest admitted number of hypotheses, actions, or outcomes on each axis.
pub const MAX_AXIS: usize = 256;
/// Absolute tolerance on the mass of a supplied probability distribution.
/// Admitted rows are divided by their measured sum before use.
pub const MASS_TOLERANCE: f64 = 1e-12;

/// A prior on hypotheses and one row of utilities per available action.
#[derive(Debug, Clone, Copy)]
pub struct FiniteDecision<'a> {
    /// Hypothesis probabilities in table-column order.
    pub prior: &'a [f64],
    /// Payoff indexed by action, then hypothesis.
    pub payoff: &'a [&'a [f64]],
}

/// A measurement and its supplied cost in the same units as the payoffs.
#[derive(Debug, Clone, Copy)]
pub struct FiniteTest<'a> {
    /// Conditional outcome probabilities, one distribution per hypothesis.
    pub likelihood: &'a [&'a [f64]],
    /// Finite, nonnegative test cost.
    pub cost: f64,
}

/// Decision value of the test, before and after its cost.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteEvsi {
    /// Gross expected utility gain from observing the test result.
    pub evsi: f64,
    /// Expected gain from observing the hypothesis itself.
    pub evpi: f64,
    /// Gross gain less the supplied test cost. Can be negative.
    pub net_value: f64,
    /// Best action for each outcome, or `None` for a zero-probability outcome.
    pub posterior_policy: Vec<Option<usize>>,
    /// Marginal outcome probabilities in the declared order.
    pub outcome_probability: Vec<f64>,
}

/// Refusals from finite decision admission or evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteError {
    /// An axis is empty, exceeds the work bound, or rows disagree on a dimension.
    InvalidShape,
    /// The prior has an invalid entry or total mass.
    InvalidPrior,
    /// The identified likelihood row has an invalid entry or total mass.
    InvalidLikelihood(usize),
    /// A payoff is infinite or NaN.
    NonFinitePayoff,
    /// Test cost is negative or nonfinite.
    InvalidCost,
    /// An intermediate probability or utility lies outside the representable range.
    NumericalRange,
}
impl core::fmt::Display for FiniteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite decision evaluation: {self:?}")
    }
}
impl std::error::Error for FiniteError {}

/// Neumaier summation retains small advantages when large signed terms cancel.
fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut total = 0.0_f64;
    let mut correction = 0.0;
    for x in values {
        let next = total + x;
        correction += if total.abs() >= x.abs() {
            (total - next) + x
        } else {
            (x - next) + total
        };
        total = next;
    }
    total + correction
}

fn mass(row: &[f64]) -> Option<f64> {
    if row.iter().any(|&p| !(0.0..=1.0).contains(&p)) {
        return None;
    }
    let total = sum(row.iter().copied());
    ((total - 1.0).abs() <= MASS_TOLERANCE).then_some(total)
}

fn finite(x: f64) -> Result<f64, FiniteError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(FiniteError::NumericalRange)
    }
}

fn prior(d: &FiniteDecision<'_>) -> Result<Vec<f64>, FiniteError> {
    if !(1..=MAX_AXIS).contains(&d.prior.len())
        || !(1..=MAX_AXIS).contains(&d.payoff.len())
        || d.payoff.iter().any(|r| r.len() != d.prior.len())
    {
        return Err(FiniteError::InvalidShape);
    }
    if d.payoff.iter().any(|r| r.iter().any(|x| !x.is_finite())) {
        return Err(FiniteError::NonFinitePayoff);
    }
    let total = mass(d.prior).ok_or(FiniteError::InvalidPrior)?;
    Ok(d.prior.iter().map(|p| p / total).collect())
}

fn weighted_gain(p: f64, a: f64, b: f64) -> f64 {
    let difference = a - b;
    if difference.is_finite() {
        p * difference
    } else {
        p * a - p * b
    }
}

fn advantage(p: &[f64], candidate: &[f64], incumbent: &[f64]) -> Result<f64, FiniteError> {
    finite(sum(
        (0..p.len()).map(|h| weighted_gain(p[h], candidate[h], incumbent[h]))
    ))
}

fn best_action(payoff: &[&[f64]], p: &[f64]) -> Result<usize, FiniteError> {
    let mut best = 0;
    for a in 1..payoff.len() {
        if advantage(p, payoff[a], payoff[best])? > 0.0 {
            best = a;
        }
    }
    Ok(best)
}

/// Best action before testing and its expected utility. Ties use the lowest index.
///
/// # Errors
/// Refuses malformed tables and unrepresentable arithmetic.
pub fn value_now(d: &FiniteDecision<'_>) -> Result<(usize, f64), FiniteError> {
    let p = prior(d)?;
    let best = best_action(d.payoff, &p)?;
    let utility = finite(sum(p.iter().zip(d.payoff[best]).map(|(p, u)| p * u)))?;
    Ok((best, utility))
}

/// Enumerate all declared outcomes. The work is O(actions × hypotheses × outcomes)
/// and auxiliary storage is O(hypotheses + outcomes).
///
/// # Errors
/// Refuses malformed tables, invalid costs and unrepresentable arithmetic.
pub fn evsi(d: &FiniteDecision<'_>, test: &FiniteTest<'_>) -> Result<FiniteEvsi, FiniteError> {
    let p = prior(d)?;
    let outcomes = test.likelihood.first().map_or(0, |row| row.len());
    if !(1..=MAX_AXIS).contains(&outcomes)
        || test.likelihood.len() != p.len()
        || test.likelihood.iter().any(|row| row.len() != outcomes)
    {
        return Err(FiniteError::InvalidShape);
    }
    if !test.cost.is_finite() || test.cost < 0.0 {
        return Err(FiniteError::InvalidCost);
    }
    let row_mass: Vec<_> = test
        .likelihood
        .iter()
        .enumerate()
        .map(|(h, row)| mass(row).ok_or(FiniteError::InvalidLikelihood(h)))
        .collect::<Result<_, _>>()?;
    let current = best_action(d.payoff, &p)?;
    let perfect: Vec<_> = (0..p.len())
        .map(|h| {
            d.payoff
                .iter()
                .map(|row| row[h])
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();
    let evpi = advantage(&p, &perfect, d.payoff[current])?;
    let mut weights = vec![0.0; p.len()];
    let mut posterior_policy = Vec::with_capacity(outcomes);
    let mut outcome_probability = Vec::with_capacity(outcomes);
    let mut gains = Vec::with_capacity(outcomes);
    let mut attains_perfect = true;
    for y in 0..outcomes {
        // Factor out the largest supported likelihood first; a rare outcome
        // must not turn a well-conditioned posterior into an accidental tie.
        let likelihood = |h: usize| test.likelihood[h][y] / row_mass[h];
        let mut maximum = 0.0_f64;
        let mut minimum = 1.0_f64;
        for h in (0..p.len()).filter(|&h| p[h] > 0.0) {
            maximum = maximum.max(likelihood(h));
            minimum = minimum.min(likelihood(h));
        }
        if maximum == 0.0 {
            posterior_policy.push(None);
            outcome_probability.push(0.0);
            continue;
        }
        // Use a power of two so scaling introduces no rounding in ordinary
        // binary fractions. The weights are proportional to the posterior;
        // dividing each by its sum would add rounding at exact policy ties.
        let exponent = maximum.max(f64::MIN_POSITIVE).to_bits() & 0x7ff0_0000_0000_0000;
        let scale = (2.0 * f64::from_bits(exponent)).min(1.0);
        for h in 0..p.len() {
            weights[h] = if p[h] == 0.0 {
                0.0
            } else {
                p[h] * (likelihood(h) / scale)
            };
            if p[h] > 0.0 && likelihood(h) > 0.0 && weights[h] == 0.0 {
                return Err(FiniteError::NumericalRange);
            }
        }
        let scaled_mass = sum(weights.iter().copied());
        let probability = finite(scale * scaled_mass)?;
        if probability == 0.0 {
            return Err(FiniteError::NumericalRange);
        }
        let action = if minimum == maximum {
            current
        } else {
            best_action(d.payoff, &weights)?
        };
        attains_perfect &=
            (0..p.len()).all(|h| weights[h] == 0.0 || d.payoff[action][h] == perfect[h]);
        posterior_policy.push(Some(action));
        outcome_probability.push(probability);
        let gain = advantage(&weights, d.payoff[action], d.payoff[current])?;
        gains.push(scale * gain.max(0.0));
    }
    let evsi = if attains_perfect {
        evpi
    } else {
        finite(sum(gains.into_iter()))?.min(evpi)
    };
    Ok(FiniteEvsi {
        evsi,
        evpi,
        net_value: evsi - test.cost,
        posterior_policy,
        outcome_probability,
    })
}
