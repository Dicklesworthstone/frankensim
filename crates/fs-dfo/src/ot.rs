//! Entropic optimal transport: Sinkhorn iterations in the LOG DOMAIN
//! (logsumexp-stabilized — small regularization without underflow),
//! deterministic schedules, and marginal-residual stopping. The
//! substrate for Wasserstein DRO and distributional robustness; the
//! known-answer gates live on 1D quadratic costs where the
//! unregularized optimum is the monotone coupling in closed form.

/// Outcome of a Sinkhorn solve.
#[derive(Debug, Clone)]
pub struct OtReport {
    /// Transport cost ⟨P, C⟩ of the entropic plan.
    pub cost: f64,
    /// The plan (row-major n×m).
    pub plan: Vec<f64>,
    /// Worst marginal residual max(‖P·1 − a‖∞, ‖Pᵀ1 − b‖∞).
    pub marginal_residual: f64,
    /// Sinkhorn iterations run.
    pub iters: usize,
}

fn logsumexp(row: &[f64]) -> f64 {
    let m = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let s: f64 = row.iter().map(|v| fs_math::det::exp(v - m)).sum();
    m + fs_math::det::ln(s)
}

/// Entropic OT between discrete measures `a` (length n) and `b`
/// (length m) with cost matrix `c` (row-major n×m), regularization
/// `epsilon`. Log-domain Sinkhorn with potentials (f, g):
/// fᵢ ← ε·ln aᵢ − ε·LSE_j((g_j − C_ij)/ε), symmetric for g.
///
/// # Panics
/// On non-positive weights or mass mismatch beyond 1e−12 (measures
/// must be normalized by the caller — an explicit contract).
#[must_use]
pub fn sinkhorn(a: &[f64], b: &[f64], c: &[f64], epsilon: f64, max_iters: usize) -> OtReport {
    let n = a.len();
    let m = b.len();
    assert_eq!(c.len(), n * m);
    assert!(a.iter().all(|&w| w > 0.0) && b.iter().all(|&w| w > 0.0));
    // The Gibbs kernel is exp(-c/epsilon); a non-positive epsilon gives
    // exp(±inf) -> 0/inf and NaN potentials. Fail closed.
    assert!(
        epsilon > 0.0,
        "sinkhorn entropic regularization epsilon must be positive"
    );
    let (sa, sb): (f64, f64) = (a.iter().sum(), b.iter().sum());
    assert!(
        (sa - sb).abs() < 1e-12,
        "balanced OT needs equal masses: {sa} vs {sb}"
    );
    let ln_a: Vec<f64> = a.iter().map(|&w| fs_math::det::ln(w)).collect();
    let ln_b: Vec<f64> = b.iter().map(|&w| fs_math::det::ln(w)).collect();
    let mut f = vec![0.0f64; n];
    let mut g = vec![0.0f64; m];
    let mut iters = 0usize;
    let mut scratch = vec![0.0f64; n.max(m)];
    for it in 0..max_iters {
        iters = it + 1;
        // f-update.
        for i in 0..n {
            for j in 0..m {
                scratch[j] = (g[j] - c[i * m + j]) / epsilon;
            }
            f[i] = epsilon * (ln_a[i] - logsumexp(&scratch[..m]));
        }
        // g-update.
        for j in 0..m {
            for i in 0..n {
                scratch[i] = (f[i] - c[i * m + j]) / epsilon;
            }
            g[j] = epsilon * (ln_b[j] - logsumexp(&scratch[..n]));
        }
        // Marginal residual every 10 iterations (row marginals are
        // exact right after the f-update; the g-update perturbs them —
        // check the ROW residual, the binding one).
        if it % 10 == 9 {
            let res = row_residual(a, &f, &g, c, epsilon, n, m);
            if res < 1e-10 {
                break;
            }
        }
    }
    // Recover the plan P_ij = exp((f_i + g_j − C_ij)/ε).
    let mut plan = vec![0.0f64; n * m];
    let mut cost = 0.0f64;
    for i in 0..n {
        for j in 0..m {
            let p = fs_math::det::exp((f[i] + g[j] - c[i * m + j]) / epsilon);
            plan[i * m + j] = p;
            cost = p.mul_add(c[i * m + j], cost);
        }
    }
    // Final residuals on both marginals.
    let mut worst = 0.0f64;
    for i in 0..n {
        let row: f64 = plan[i * m..(i + 1) * m].iter().sum();
        worst = worst.max((row - a[i]).abs());
    }
    for j in 0..m {
        let col: f64 = (0..n).map(|i| plan[i * m + j]).sum();
        worst = worst.max((col - b[j]).abs());
    }
    OtReport {
        cost,
        plan,
        marginal_residual: worst,
        iters,
    }
}

fn row_residual(
    a: &[f64],
    f: &[f64],
    g: &[f64],
    c: &[f64],
    epsilon: f64,
    n: usize,
    m: usize,
) -> f64 {
    let mut worst = 0.0f64;
    for i in 0..n {
        let mut row = 0.0f64;
        for j in 0..m {
            row += fs_math::det::exp((f[i] + g[j] - c[i * m + j]) / epsilon);
        }
        worst = worst.max((row - a[i]).abs());
    }
    worst
}

/// Squared-distance cost matrix for 1D point clouds.
#[must_use]
pub fn cost_sq_1d(x: &[f64], y: &[f64]) -> Vec<f64> {
    let mut c = Vec::with_capacity(x.len() * y.len());
    for &xi in x {
        for &yj in y {
            c.push((xi - yj) * (xi - yj));
        }
    }
    c
}

/// CLOSED-FORM 1D quadratic-cost OT between equal-weight point sets
/// of the same size: the optimum is the monotone (sorted) coupling —
/// the known-answer oracle for the batteries.
#[must_use]
pub fn monotone_cost_1d(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len());
    let mut xs = x.to_vec();
    let mut ys = y.to_vec();
    xs.sort_by(f64::total_cmp);
    ys.sort_by(f64::total_cmp);
    let n = xs.len() as f64;
    xs.iter()
        .zip(&ys)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        / n
}
