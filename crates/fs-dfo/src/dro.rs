//! Wasserstein distributionally-robust inner sup with DISCRETE
//! candidate support: sup over {Q : W(Q, P̂_n) ≤ ρ, supp(Q) ⊆ {z_j}}
//! of E_Q[ℓ] via the EXACT dual reformulation
//! (Mohajerin Esfahani–Kuhn form):
//!
//!   min_{λ ≥ 0}  λ·ρ + (1/n)·Σᵢ max_j (ℓ_j − λ·C_ij),
//!
//! convex piecewise-linear in λ — minimized by deterministic
//! golden-section over a bracketing interval derived from the data.
//! Strong duality is EXHIBITED in the battery against an exact
//! tiny-LP enumeration oracle, not assumed.

/// DRO inner-sup report.
#[derive(Debug, Clone)]
pub struct DroReport {
    /// The worst-case expectation sup E_Q[ℓ].
    pub worst_case: f64,
    /// Optimal dual multiplier λ*.
    pub lambda: f64,
    /// Worst-case distribution over the candidate support (recovered
    /// from the argmax transport structure; ties broken by index).
    pub q: Vec<f64>,
}

fn validate_inputs(losses: &[f64], costs: &[f64], n: usize, rho: f64) -> usize {
    assert!(n > 0, "DRO empirical sample count must be positive");
    assert!(
        !losses.is_empty() && losses.iter().all(|l| l.is_finite()),
        "DRO losses must be non-empty and finite"
    );
    assert!(
        rho.is_finite() && rho >= 0.0,
        "DRO radius must be finite and non-negative"
    );
    let m = losses.len();
    assert!(
        n <= usize::MAX / m,
        "DRO cost matrix dimensions overflow usize"
    );
    let expected = n * m;
    assert_eq!(
        costs.len(),
        expected,
        "DRO costs must be row-major n*m entries"
    );
    assert!(
        costs.iter().all(|c| c.is_finite() && *c >= 0.0),
        "DRO transport costs must be finite and non-negative"
    );
    m
}

/// Dual objective g(λ) = λρ + (1/n)Σᵢ max_j (ℓ_j − λC_ij).
fn dual(lambda: f64, rho: f64, losses: &[f64], costs: &[f64], n: usize, m: usize) -> f64 {
    let mut acc = 0.0f64;
    for i in 0..n {
        let mut best = f64::NEG_INFINITY;
        for (j, &lj) in losses.iter().enumerate() {
            best = best.max(lambda.mul_add(-costs[i * m + j], lj));
        }
        acc += best;
    }
    lambda.mul_add(rho, acc / n as f64)
}

fn cheapest_max_loss_distribution(
    losses: &[f64],
    costs: &[f64],
    n: usize,
    m: usize,
    lmax: f64,
) -> (Vec<f64>, f64) {
    let loss_scale = losses
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0f64, f64::max)
        .max(f64::MIN_POSITIVE);
    let max_loss_tol = 64.0 * f64::EPSILON * loss_scale;
    let mut q = vec![0.0f64; m];
    let inv_n = 1.0 / n as f64;
    let mut cost = 0.0f64;
    for i in 0..n {
        let row = &costs[i * m..(i + 1) * m];
        let mut best_j = 0usize;
        let mut best_cost = f64::INFINITY;
        for (j, &loss) in losses.iter().enumerate() {
            if (lmax - loss).abs() <= max_loss_tol {
                let candidate_cost = row[j];
                if candidate_cost < best_cost
                    || (candidate_cost.to_bits() == best_cost.to_bits() && j < best_j)
                {
                    best_j = j;
                    best_cost = candidate_cost;
                }
            }
        }
        q[best_j] += inv_n;
        cost += best_cost * inv_n;
    }
    (q, cost)
}

fn recover_distribution(
    lambda: f64,
    rho: f64,
    losses: &[f64],
    costs: &[f64],
    n: usize,
    m: usize,
) -> Vec<f64> {
    let mut q = vec![0.0f64; m];
    let inv_n = 1.0 / n as f64;

    let scale = losses
        .iter()
        .copied()
        .map(f64::abs)
        .chain(costs.iter().copied().map(|c| (lambda * c).abs()))
        .fold(0.0f64, f64::max)
        .max(f64::MIN_POSITIVE);
    let active_tol = 512.0 * f64::EPSILON * scale;
    let mut bases = Vec::with_capacity(n);
    let mut caps = Vec::new();
    let mut base_cost = 0.0f64;
    for i in 0..n {
        let row = &costs[i * m..(i + 1) * m];
        let mut best = f64::NEG_INFINITY;
        let mut best_idx = 0usize;
        for (j, (&loss, &cost)) in losses.iter().zip(row).enumerate() {
            let reduced = lambda.mul_add(-cost, loss);
            if reduced > best {
                best = reduced;
                best_idx = j;
            }
        }
        let mut lo = best_idx;
        let mut hi = best_idx;
        for (j, (&loss, &cost)) in losses.iter().zip(row).enumerate() {
            let reduced = lambda.mul_add(-cost, loss);
            if best - reduced <= active_tol {
                if cost < row[lo] || (cost.to_bits() == row[lo].to_bits() && j < lo) {
                    lo = j;
                }
                if cost > row[hi] || (cost.to_bits() == row[hi].to_bits() && j > hi) {
                    hi = j;
                }
            }
        }
        q[lo] += inv_n;
        bases.push((lo, hi));
        base_cost += row[lo] * inv_n;
        let cap = (row[hi] - row[lo]) * inv_n;
        if cap > 0.0 {
            caps.push((i, cap));
        }
    }
    let mut remaining = (rho - base_cost).max(0.0);
    for (i, cap) in caps {
        if remaining <= 1e-14 {
            break;
        }
        let take = remaining.min(cap);
        let alpha = take / cap;
        let (lo, hi) = bases[i];
        q[lo] -= alpha * inv_n;
        q[hi] += alpha * inv_n;
        remaining -= take;
    }
    q
}

/// Worst-case expectation of `losses[j]` at support points z_j, over
/// the Wasserstein ball of radius `rho` around the empirical measure
/// of `n` samples, with transport costs `costs` (row-major n×m,
/// C_ij = c(ξ_i, z_j); each row must contain a zero-cost entry —
/// "stay put" — or the ball may be infeasible at ρ = 0, asserted).
#[must_use]
pub fn wasserstein_worst_case(losses: &[f64], costs: &[f64], n: usize, rho: f64) -> DroReport {
    let m = validate_inputs(losses, costs, n, rho);
    for i in 0..n {
        let row_min = (0..m)
            .map(|j| costs[i * m + j])
            .fold(f64::INFINITY, f64::min);
        assert!(
            row_min < 1e-12,
            "each sample needs a zero-cost support option (row {i} min {row_min})"
        );
    }
    // λ upper bracket: beyond max ℓ-range / min positive cost, every
    // argmax is the zero-cost option and g'(λ) = ρ ≥ 0 — the minimum
    // lies below.
    let lmax = losses.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let lmin = losses.iter().copied().fold(f64::INFINITY, f64::min);
    let min_pos_cost = costs
        .iter()
        .copied()
        .filter(|&c| c > 1e-12)
        .fold(f64::INFINITY, f64::min);
    let hi = if min_pos_cost.is_finite() {
        ((lmax - lmin) / min_pos_cost).max(1.0) * 2.0
    } else {
        1.0
    };
    // Golden-section on the convex dual (deterministic iteration
    // count sized for 1e-12 relative brackets).
    let phi = 0.618_033_988_749_894_9f64;
    let (mut a, mut b) = (0.0f64, hi);
    let mut c = b - phi * (b - a);
    let mut d = a + phi * (b - a);
    let mut fc = dual(c, rho, losses, costs, n, m);
    let mut fd = dual(d, rho, losses, costs, n, m);
    for _ in 0..200 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = dual(c, rho, losses, costs, n, m);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = dual(d, rho, losses, costs, n, m);
        }
    }
    let lambda = f64::midpoint(a, b);
    let worst_case = dual(lambda, rho, losses, costs, n, m)
        .min(dual(0.0, rho, losses, costs, n, m))
        .min(lmax); // the primal can never exceed the max loss
    let (saturated_q, saturated_cost) = cheapest_max_loss_distribution(losses, costs, n, m, lmax);
    let q = if saturated_cost <= rho + 1e-12 {
        saturated_q
    } else {
        recover_distribution(lambda, rho, losses, costs, n, m)
    };
    DroReport {
        worst_case,
        lambda,
        q,
    }
}
