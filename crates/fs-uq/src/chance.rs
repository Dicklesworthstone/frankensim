//! CHANCE CONSTRAINTS (bead qlvf, lane b): `P(g(x, ξ) ≤ 0) ≥ 1 − α`
//! lowered through the Chance treatment with ANYTIME-VALID probability
//! estimates. The complete candidate-query family is declared before sampling;
//! each query receives a Bonferroni share of the familywise CS budget, and only
//! candidates whose simultaneous lower bound is feasible may be returned.

use crate::anytime::estimate_probability_anytime;

const GRID_INTERVALS: u64 = 48;
const GRID_POINTS: u64 = GRID_INTERVALS + 1;
const MAX_CHANCE_QUERIES: u64 = 3_136;
const MAX_CHANCE_SAMPLES: u64 = 10_000_000;

/// A chance-constrained design backed by one member of a simultaneous
/// confidence-sequence family.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChanceSolution {
    /// Selected scalar design.
    pub x: f64,
    /// Sample probability at the selected design.
    pub probability: f64,
    /// Simultaneous lower confidence bound used for admission.
    pub probability_lo: f64,
    /// Simultaneous upper confidence bound.
    pub probability_hi: f64,
    /// Total samples consumed over all candidate queries.
    pub samples: u64,
    /// Number of confidence sequences in the predeclared family.
    pub queries: u64,
    /// Familywise error budget divided across those queries.
    pub familywise_alpha: f64,
}

/// The bounded scan found no candidate with an admitted feasible lower bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoFeasibleChanceCandidate;

impl core::fmt::Display for NoFeasibleChanceCandidate {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_str("no scanned candidate has a simultaneously feasible lower bound")
    }
}

impl std::error::Error for NoFeasibleChanceCandidate {}

/// CHANCE-CONSTRAINED minimization of `objective` over a scalar box:
/// `min f(x)` subject to `P(g(x, ξ) ≤ 0) ≥ 1 − α`, with the
/// probability estimated ANYTIME-VALIDLY per query. `sample_g(x, i)` draws
/// `g(x, ξᵢ)` with the i-th deterministic germ. The predeclared family spends
/// `cs_alpha / (outer_iters * 49)` on each query, so adaptive selection of the
/// returned candidate preserves the familywise confidence level.
///
/// # Errors
/// [`NoFeasibleChanceCandidate`] when no scanned candidate has a simultaneous
/// lower bound at least `1 - alpha`.
///
/// # Panics
/// Malformed scalar controls, non-finite callback outputs, or a requested work
/// envelope above ten million samples are programmer-contract violations and
/// fail before publishing a solution.
#[allow(clippy::too_many_arguments, clippy::cast_precision_loss)]
pub fn chance_constrained_min(
    objective: impl Fn(f64) -> f64,
    mut sample_g: impl FnMut(f64, u64) -> f64,
    alpha: f64,
    bounds: (f64, f64),
    outer_iters: usize,
    cs_alpha: f64,
    cs_half_width: f64,
    max_samples: u64,
) -> Result<ChanceSolution, NoFeasibleChanceCandidate> {
    let (lo, hi) = bounds;
    assert!(
        alpha.is_finite() && alpha > 0.0 && alpha < 1.0,
        "chance level alpha must be finite and lie in (0,1)"
    );
    assert!(
        cs_alpha.is_finite() && cs_alpha > 0.0 && cs_alpha < 1.0,
        "familywise CS alpha must be finite and lie in (0,1)"
    );
    assert!(
        lo.is_finite() && hi.is_finite() && lo < hi,
        "chance-constraint bounds must be finite and ordered"
    );
    assert!(outer_iters > 0, "at least one chance scan is required");
    assert!(
        max_samples > 0,
        "each chance query needs a positive sample cap"
    );
    let queries = u64::try_from(outer_iters)
        .ok()
        .and_then(|count| count.checked_mul(GRID_POINTS))
        .expect("chance-query count overflow");
    assert!(
        queries <= MAX_CHANCE_QUERIES,
        "chance-query count {queries} exceeds cap {MAX_CHANCE_QUERIES}"
    );
    let requested_samples = queries
        .checked_mul(max_samples)
        .expect("chance sample budget overflow");
    assert!(
        requested_samples <= MAX_CHANCE_SAMPLES,
        "chance sample envelope {requested_samples} exceeds cap {MAX_CHANCE_SAMPLES}"
    );
    let per_query_alpha = cs_alpha / queries as f64;
    let required_probability = 1.0 - alpha;
    let mut spent_total = 0u64;
    let mut best: Option<(f64, ChanceSolution)> = None;
    for outer in 0..outer_iters {
        let outer = u64::try_from(outer).expect("bounded outer index fits u64");
        let sample_base = outer
            .checked_mul(max_samples)
            .expect("bounded sample offset cannot overflow");
        for k in 0..=GRID_INTERVALS {
            let cand = lo + (hi - lo) * k as f64 / GRID_INTERVALS as f64;
            // Anytime-stopped probability estimate at this candidate.
            let est = estimate_probability_anytime(
                |i| {
                    let sample_index = sample_base
                        .checked_add(i)
                        .expect("bounded sample index cannot overflow");
                    let v = sample_g(cand, sample_index);
                    assert!(v.is_finite(), "chance-constraint sample must be finite");
                    f64::from(u8::from(v <= 0.0))
                },
                per_query_alpha,
                cs_half_width,
                max_samples,
            );
            spent_total = spent_total
                .checked_add(est.n)
                .expect("bounded total sample count cannot overflow");
            if est.lo >= required_probability {
                let objective_value = objective(cand);
                assert!(
                    objective_value.is_finite(),
                    "chance-constraint objective must be finite"
                );
                let solution = ChanceSolution {
                    x: cand,
                    probability: est.mean,
                    probability_lo: est.lo,
                    probability_hi: est.hi,
                    samples: spent_total,
                    queries,
                    familywise_alpha: cs_alpha,
                };
                if best.is_none_or(|(best_value, _)| objective_value < best_value) {
                    best = Some((objective_value, solution));
                }
            }
        }
    }
    let Some((_, mut solution)) = best else {
        return Err(NoFeasibleChanceCandidate);
    };
    solution.samples = spent_total;
    Ok(solution)
}
