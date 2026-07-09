//! CHANCE CONSTRAINTS (bead qlvf, lane b): `P(g(x, ξ) ≤ 0) ≥ 1 − α`
//! lowered through the Chance treatment with ANYTIME-VALID probability
//! estimates — the e-process CS stops each probability query the
//! moment it is decision-grade, and the outer augmented-penalty loop
//! consumes the stopped estimates, using the CS LOWER bound so the
//! constraint holds even at the pessimistic end (validity feeding
//! feasibility).

use crate::anytime::estimate_probability_anytime;

/// CHANCE-CONSTRAINED minimization of `objective` over a scalar box:
/// `min f(x)` subject to `P(g(x, ξ) ≤ 0) ≥ 1 − α`, with the
/// probability estimated ANYTIME-VALIDLY per outer iterate (the CS
/// stops each query at decision grade) and enforced by an augmented
/// penalty on the shortfall. `sample_g(x, i)` draws g(x, ξᵢ) with the
/// i-th deterministic germ.
#[allow(clippy::too_many_arguments)]
pub fn chance_constrained_min(
    objective: impl Fn(f64) -> f64,
    mut sample_g: impl FnMut(f64, u64) -> f64,
    alpha: f64,
    bounds: (f64, f64),
    outer_iters: usize,
    cs_alpha: f64,
    cs_half_width: f64,
    max_samples: u64,
) -> (f64, f64, u64) {
    let (lo, hi) = bounds;
    let mut x = f64::midpoint(lo, hi);
    let mut penalty = 10.0f64;
    let mut spent_total = 0u64;
    let mut p_est = 0.0f64;
    for outer in 0..outer_iters {
        // Golden-section-ish scan of the penalized objective on the
        // box (deterministic; the objective is 1-D by fixture design).
        let mut best = (f64::INFINITY, x);
        let grid = 48;
        for k in 0..=grid {
            let cand = lo + (hi - lo) * f64::from(k) / f64::from(grid);
            // Anytime-stopped probability estimate at this candidate.
            let est = estimate_probability_anytime(
                |i| {
                    let v = sample_g(cand, i + (outer as u64) * 1_000_003);
                    f64::from(u8::from(v <= 0.0))
                },
                cs_alpha,
                cs_half_width,
                max_samples,
            );
            spent_total += est.n;
            // Use the CS LOWER bound: the constraint must hold even at
            // the pessimistic end (validity feeding feasibility).
            let shortfall = ((1.0 - alpha) - est.lo).max(0.0);
            let val = objective(cand) + penalty * shortfall * shortfall;
            if val < best.0 {
                best = (val, cand);
                p_est = est.mean;
            }
        }
        x = best.1;
        penalty *= 3.0;
    }
    (x, p_est, spent_total)
}
