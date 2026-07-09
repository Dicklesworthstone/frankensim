//! The co-optimizing allocator (bead gp3.9 V2, feature
//! `moonshot-planner`, [M]): exact where the problem is discrete,
//! convex where the models are rate-based, CMA-ES on itself where
//! they are not.
//!
//! - [`optimize_exact`]: per-track multiple-choice-knapsack DP. The
//!   tropical wall-clock (max over track sums ≤ B) DECOMPOSES: every
//!   track independently obeys the budget, so each track solves its
//!   own knapsack exactly (up to the cost grid) and errors add.
//! - [`waterfill`]: the convex co-optimizer for rate-based models
//!   e_k(c) = a_k·c^{−p_k} — KKT water-filling by bisection on the
//!   Lagrange multiplier (deterministic, certificate-grade for convex
//!   curves).
//! - [`cma_continuous`]: CMA-ES (fs-dfo) on the continuous allocation
//!   simplex for models that are NOT convex — the "CMA-ES on itself"
//!   the plan names. The battery cross-checks it against water-filling
//!   on convex fixtures and shows it escaping the basin water-filling
//!   rounds into on a non-convex one.
//!
//! PROMOTION GATE (real, and not passed here): V2 promotes only after
//! the Gauntlet shows its plans beat hand allocation on all three
//! flagships (huq.15). The fixture-matrix [`Scoreboard`] in the
//! battery is necessary evidence, not sufficient — the feature ships
//! OFF by default.

use crate::alloc::{AllocProblem, Plan, plan_total_error, plan_wall_clock};

/// Exact (grid-discretized) multiple-choice knapsack per track.
/// `grid` = number of cost buckets per track (resolution B/grid).
///
/// # Panics
/// If `grid == 0`.
#[must_use]
pub fn optimize_exact(problem: &AllocProblem, grid: usize) -> Option<Plan> {
    assert!(grid > 0, "cost grid must be positive");
    let knobs = &problem.knobs;
    let ntracks = knobs.iter().map(|k| k.track).max().unwrap_or(0) + 1;
    let step = problem.budget_s / grid as f64;
    let mut choice = vec![0usize; knobs.len()];
    for t in 0..ntracks {
        let members: Vec<usize> = (0..knobs.len()).filter(|&i| knobs[i].track == t).collect();
        if members.is_empty() {
            continue;
        }
        // DP over budget buckets: dp[b] = (min error, back-pointers).
        let inf = f64::INFINITY;
        let mut dp = vec![0.0f64; grid + 1];
        let mut back: Vec<Vec<usize>> = vec![Vec::new(); grid + 1];
        let mut initialized = false;
        for &ki in &members {
            let k = &knobs[ki];
            let mut ndp = vec![inf; grid + 1];
            let mut nback: Vec<Vec<usize>> = vec![Vec::new(); grid + 1];
            for b in 0..=grid {
                if initialized && !dp[b].is_finite() {
                    continue;
                }
                let base_err = if initialized { dp[b] } else { 0.0 };
                for (si, s) in k.settings.iter().enumerate() {
                    // Ceil-bucket the cost (conservative: never
                    // understates wall-clock).
                    let buckets = (s.cost / step).ceil() as usize;
                    let nb = b + buckets;
                    if nb > grid {
                        continue;
                    }
                    let e = base_err + s.error;
                    if e < ndp[nb] {
                        ndp[nb] = e;
                        let mut path = back[b].clone();
                        path.push(si);
                        nback[nb] = path;
                    }
                }
                if !initialized {
                    break; // first knob: only the b = 0 row seeds
                }
            }
            dp = ndp;
            back = nback;
            initialized = true;
        }
        // Best bucket for this track.
        let mut best: Option<(f64, usize)> = None;
        for (b, &e) in dp.iter().enumerate() {
            if e.is_finite() && best.is_none_or(|(be, _)| e < be) {
                best = Some((e, b));
            }
        }
        let (_, bb) = best?;
        for (mi, &ki) in members.iter().enumerate() {
            choice[ki] = back[bb][mi];
        }
    }
    let total_error = plan_total_error(knobs, &choice);
    let wall = plan_wall_clock(knobs, &choice);
    if wall > problem.budget_s {
        return None;
    }
    Some(Plan {
        choice,
        total_error,
        wall_clock: wall,
        rationale: vec!["exact per-track knapsack DP (moonshot V2)".to_owned()],
    })
}

/// A rate-based continuous knob: error(c) = a·c^{−p} for cost c > 0.
#[derive(Debug, Clone, Copy)]
pub struct RateKnob {
    /// Error prefactor.
    pub a: f64,
    /// Convergence rate exponent (> 0).
    pub p: f64,
}

/// Convex water-filling: minimize Σ aᵏ·cₖ^{−pₖ} s.t. Σ cₖ = budget,
/// by bisection on the KKT multiplier λ (aₖpₖcₖ^{−pₖ−1} = λ).
/// Returns the per-knob cost allocation.
#[must_use]
pub fn waterfill(knobs: &[RateKnob], budget: f64) -> Vec<f64> {
    let c_of_lambda = |lam: f64| -> Vec<f64> {
        knobs
            .iter()
            .map(|k| (k.a * k.p / lam).powf(1.0 / (k.p + 1.0)))
            .collect()
    };
    let total = |lam: f64| -> f64 { c_of_lambda(lam).iter().sum() };
    // Bracket: λ large → cheap; λ small → expensive.
    let (mut lo, mut hi) = (1e-18f64, 1e18f64);
    for _ in 0..200 {
        let mid = fs_math_sqrt(lo * hi); // geometric bisection
        if total(mid) > budget {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    c_of_lambda(fs_math_sqrt(lo * hi))
}

fn fs_math_sqrt(x: f64) -> f64 {
    // std sqrt is IEEE-correctly-rounded; keep the helper local to
    // avoid a whole fs-math dependency for one op.
    f64::sqrt(x)
}

/// Total rate-model error of a continuous allocation.
#[must_use]
pub fn rate_error(knobs: &[RateKnob], alloc: &[f64]) -> f64 {
    knobs
        .iter()
        .zip(alloc)
        .map(|(k, &c)| k.a * c.powf(-k.p))
        .sum()
}

/// CMA-ES on the continuous allocation (fs-dfo): parameterize the
/// budget split through a softmax so every sample is feasible, and
/// minimize an arbitrary (possibly NON-convex) error model. Uses
/// BIPOP restarts — MEASURED rejection: single-run CMA-ES converges
/// AWAY from activation cliffs (the surrogate-threshold fixture ended
/// at spend 0.0 with error 0.58 vs the 0.16 basin behind the cliff:
/// once σ shrinks below the cliff distance no sample ever crosses it;
/// the restart schedule with large populations does). Deterministic
/// via the fixed seed.
#[must_use]
pub fn cma_continuous<F: FnMut(&[f64]) -> f64>(
    nknobs: usize,
    budget: f64,
    error_of_alloc: &mut F,
    seed: u64,
) -> Vec<f64> {
    let mut objective = |z: &[f64]| -> f64 {
        let alloc = softmax_alloc(z, budget);
        error_of_alloc(&alloc)
    };
    let x0 = vec![0.0f64; nknobs];
    let report =
        fs_dfo::cma::bipop_cmaes(&mut objective, &x0, 1.5, 20_000, f64::NEG_INFINITY, seed);
    softmax_alloc(&report.best.x_best, budget)
}

fn softmax_alloc(z: &[f64], budget: f64) -> Vec<f64> {
    let m = z.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let exps: Vec<f64> = z.iter().map(|&v| (v - m).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.iter().map(|e| e / sum * budget).collect()
}

/// One scoreboard row: a fixture's accuracy-per-second comparison —
/// the promotion-gate evidence shape (necessary, not sufficient).
#[derive(Debug, Clone)]
pub struct ScoreRow {
    /// Fixture label.
    pub fixture: String,
    /// Errors achieved within the shared budget.
    pub hand_error: f64,
    /// Greedy V1.
    pub v1_error: f64,
    /// Co-optimizer V2.
    pub v2_error: f64,
}

impl ScoreRow {
    /// Does V2 beat (or tie) both hand allocation and V1 here?
    #[must_use]
    pub fn v2_wins(&self) -> bool {
        self.v2_error <= self.v1_error && self.v2_error <= self.hand_error
    }
}
