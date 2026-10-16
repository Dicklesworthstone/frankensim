//! Select complete experiments by regularized local observation information.
//!
//! Greedy block D-design: maximize log det(I + sum J_g^T J_g / ridge).
//! Rows must already have their intended coordinate and residual scaling.
//! The ridge is an explicit numerical design preference, not an inferred prior
//! or noise precision. No likelihood, posterior, or global optimum is supplied.
use fs_la::factor::{FactorError, cholesky};
use fs_math::det;

/// One indivisible experiment: all its observation rows are selected together.
#[derive(Clone, Debug, PartialEq)]
pub struct ExperimentGroup {
    /// Distinct, nonempty label; input order breaks numerically unresolved ties.
    pub name: String,
    /// Row-major weighted observation derivatives in common parameter units.
    pub rows: Vec<Vec<f64>>,
}

/// Explicit cardinality, numerical regularization and scoring-work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExperimentDesignConfig {
    /// Common number of columns, in 1..=32.
    pub parameters: usize,
    /// Select at most this many groups, in 1..=number of groups (at most 64).
    pub maximum_selected: usize,
    /// Strictly positive finite precision added on every parameter diagonal.
    /// Its meaning changes if the caller changes coordinate or residual scales.
    pub ridge_precision: f64,
    /// Maximum dense score factorizations, in 1..=4096, including the baseline.
    /// A whole candidate round is reserved before scoring it; no partial-round
    /// winner can depend on input order merely because a budget ran out.
    pub maximum_factorizations: usize,
}

impl ExperimentDesignConfig {
    /// Validate dimensions and caps without derivatives, allocations or physics.
    /// The total row cap is 1024; an individual group must also be nonempty.
    pub fn admit(self, groups: usize, rows: usize) -> Result<(), ExperimentDesignError> {
        if !(1..=32).contains(&self.parameters) || !(1..=64).contains(&groups)
            || !(1..=1024).contains(&rows) || self.maximum_selected == 0
            || self.maximum_selected > groups || self.maximum_factorizations == 0
            || self.maximum_factorizations > 4096
        {
            return Err(ExperimentDesignError::Invalid("experiment dimensions or budget exceed admission"));
        }
        if !self.ridge_precision.is_finite() || self.ridge_precision <= 0.0 {
            return Err(ExperimentDesignError::Invalid("ridge precision must be finite and positive"));
        }
        Ok(())
    }
}

/// Reason selection ended; none asserts a globally optimal subset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExperimentDesignStop {
    /// The requested maximum number of whole experiments was selected.
    SelectionLimit,
    /// The remaining scoring budget cannot cover a whole greedy round.
    EvaluationLimit,
    /// No remaining score increase exceeds the stated floating-point screen.
    NoResolvedGain,
}

/// One accepted whole-experiment addition.
#[derive(Clone, Debug, PartialEq)]
pub struct ExperimentDesignStep {
    /// Index in the original group array (not a flattened observation row).
    pub group: usize,
    /// Increase from the preceding selected set.
    pub marginal_gain: f64,
    /// log det(I + sum selected J_g^T J_g / ridge).
    pub log_determinant_gain: f64,
}

/// Complete accepted selection, even when the scoring budget stops it early.
#[derive(Clone, Debug, PartialEq)]
pub struct ExperimentDesignReport {
    /// Group indices in greedy selection order; groups are never repeated.
    pub selected: Vec<usize>,
    /// Accepted additions, with one entry for each selected group.
    pub steps: Vec<ExperimentDesignStep>,
    /// Final gain relative to the caller's isotropic regularization alone.
    pub log_determinant_gain: f64,
    /// Actual score factorizations, including the identity baseline.
    pub factorizations: usize,
    /// Cardinality, numerical gain or scoring-budget termination.
    pub stop: ExperimentDesignStop,
}

/// Structural, numerical or cooperative-cancellation refusal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExperimentDesignError {
    /// Invalid dimensions, names, rows or controls.
    Invalid(&'static str),
    /// Finite input could not be scored in the admitted numerical range.
    NumericalRange,
    /// The existing Cholesky owner refused a regularized information matrix.
    Factor(FactorError),
    /// Cancellation observed before publishing the complete report.
    Cancelled,
}
impl core::fmt::Display for ExperimentDesignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(why) => write!(f, "experiment design: {why}"),
            Self::NumericalRange => write!(f, "experiment information exceeds its numerical range"),
            Self::Factor(error) => write!(f, "experiment information: {error}"),
            Self::Cancelled => write!(f, "experiment design cancelled"),
        }
    }
}
impl std::error::Error for ExperimentDesignError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ExperimentDesignError> {
    if cancelled() { Err(ExperimentDesignError::Cancelled) } else { Ok(()) }
}
fn finite(x: f64) -> Result<f64, ExperimentDesignError> {
    if x.is_finite() { Ok(x) } else { Err(ExperimentDesignError::NumericalRange) }
}
fn score(matrix: &[f64], n: usize, cancelled: &mut impl FnMut() -> bool)
    -> Result<f64, ExperimentDesignError>
{
    poll(cancelled)?;
    let factor = cholesky(matrix, n).map_err(ExperimentDesignError::Factor)?;
    // The existing bounded (n<=32) dense factorization is not internally polled.
    let mut result = 0.0;
    for i in 0..n { result = finite(result + 2.0 * det::ln(factor.l(i, i)))?; }
    poll(cancelled)?;
    Ok(result)
}
fn gain_screen(a: f64, b: f64, n: usize) -> f64 {
    128.0 * f64::EPSILON * (1.0 + a.abs() + b.abs() + n as f64)
}

/// Greedily select complementary groups using their complete local Jacobians.
///
/// All groups are admitted before any score. Each round scores every unselected
/// group against the same accepted set, then publishes only its best resolved
/// increase. Numerical ties retain the lower original index. Zero-information
/// groups do not fill the requested cardinality with fabricated improvement.
///
/// Input matrices are immutable. Cancellation or numerical failure returns no
/// partial report; a work-cap stop returns the last fully selected set and an
/// explicit `EvaluationLimit`. There are no callbacks to a simulator. Dense
/// scoring is capped at 32 columns, 64 groups, 1024 rows and 4096 factorizations.
///
/// This is a greedy local quadratic-design heuristic, not exhaustive subset
/// optimization, rank certification or experimental validation. Additive row
/// information does not model correlated measurement errors across groups.
/// All scales and the ridge must be fixed before comparing candidate subsets.
pub fn select_experiments(groups: &[ExperimentGroup], config: ExperimentDesignConfig,
    mut cancelled: impl FnMut() -> bool) -> Result<ExperimentDesignReport, ExperimentDesignError>
{
    poll(&mut cancelled)?;
    // Reject large outer counts before traversing caller-owned row arrays.
    config.admit(groups.len(), 1)?;
    let rows = groups.iter().try_fold(0usize, |sum, group| sum.checked_add(group.rows.len()))
        .ok_or(ExperimentDesignError::Invalid("experiment row count overflow"))?;
    config.admit(groups.len(), rows)?;
    let n = config.parameters;
    for (i, group) in groups.iter().enumerate() {
        poll(&mut cancelled)?;
        if group.name.trim().is_empty() || group.name.len() > 128
            || groups[..i].iter().any(|g| g.name == group.name) || group.rows.is_empty()
            || group.rows.iter().any(|r| r.len() != n || r.iter().any(|x| !x.is_finite()))
        {
            return Err(ExperimentDesignError::Invalid("groups need unique names and complete finite nonempty Jacobians"));
        }
    }
    let root = det::sqrt(config.ridge_precision);
    let mut grams = Vec::with_capacity(groups.len());
    for group in groups {
        let mut gram = vec![0.0; n*n];
        for row in &group.rows {
            poll(&mut cancelled)?;
            let scaled: Vec<f64> = row.iter().map(|&x| {
                let y = finite(x / root)?;
                if x != 0.0 && y == 0.0 { Err(ExperimentDesignError::NumericalRange) } else { Ok(y) }
            }).collect::<Result<_, _>>()?;
            for i in 0..n {
                for j in 0..=i {
                    gram[i*n+j] = finite(scaled[i].mul_add(scaled[j], gram[i*n+j]))?;
                }
            }
        }
        grams.push(gram);
    }
    let mut matrix = vec![0.0; n*n];
    for i in 0..n { matrix[i*n+i] = 1.0; }
    let baseline = score(&matrix, n, &mut cancelled)?;
    let mut report = ExperimentDesignReport { selected: Vec::new(), steps: Vec::new(),
        log_determinant_gain: baseline, factorizations: 1, stop: ExperimentDesignStop::SelectionLimit };
    let mut used = vec![false; groups.len()];
    while report.selected.len() < config.maximum_selected {
        poll(&mut cancelled)?;
        if groups.len() - report.selected.len() > config.maximum_factorizations - report.factorizations {
            report.stop = ExperimentDesignStop::EvaluationLimit;
            break;
        }
        let mut best: Option<(usize, f64, Vec<f64>)> = None;
        for (i, gram) in grams.iter().enumerate() {
            if used[i] { continue; }
            poll(&mut cancelled)?;
            let trial: Vec<f64> = matrix.iter().zip(gram).map(|(a,b)| finite(a+b)).collect::<Result<_,_>>()?;
            report.factorizations += 1;
            let candidate = score(&trial, n, &mut cancelled)?;
            let floor = gain_screen(candidate, report.log_determinant_gain, n);
            if candidate < report.log_determinant_gain - floor {
                return Err(ExperimentDesignError::NumericalRange);
            }
            if candidate <= report.log_determinant_gain + floor { continue; }
            if best.as_ref().is_none_or(|(_, value, _)| candidate > *value + gain_screen(candidate, *value, n)) {
                best = Some((i, candidate, trial));
            }
        }
        let Some((group, value, trial)) = best else {
            report.stop = ExperimentDesignStop::NoResolvedGain;
            break;
        };
        poll(&mut cancelled)?;
        report.steps.push(ExperimentDesignStep { group, marginal_gain: value - report.log_determinant_gain,
            log_determinant_gain: value });
        report.selected.push(group); used[group] = true; matrix = trial;
        report.log_determinant_gain = value;
    }
    poll(&mut cancelled)?;
    Ok(report)
}
