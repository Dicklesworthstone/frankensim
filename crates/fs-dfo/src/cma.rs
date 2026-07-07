//! CMA-ES as natural-gradient IGO (plan §9.3, Bet 6): weighted
//! recombination with log-rank weights, rank-µ + rank-1 covariance
//! updates, cumulative step-size adaptation — the standard Hansen
//! couplings, which ARE the natural-gradient couplings on the Gaussian
//! statistical manifold. Rank-based selection makes the evolution
//! invariant to monotone transforms of the objective BY CONSTRUCTION
//! (property-tested bitwise, not cited).
//!
//! Determinism: sampling from a keyed Philox stream, `total_cmp` ranking
//! with lowest-index tie-breaks, fixed eigendecomposition cadence via
//! the landed cyclic Jacobi — the trajectory is a pure function of the
//! seed.

use fs_la::eigen::jacobi_eigh;
use fs_rand::StreamKey;

/// Kernel id for CMA sampling streams (stable registry).
const K_CMA: u32 = 0xD1F0;

/// Tunables (defaults follow Hansen's standard settings).
#[derive(Debug, Clone)]
pub struct CmaParams {
    /// Population size λ (default 4 + ⌊3·ln n⌋).
    pub lambda: usize,
    /// Initial step size σ₀.
    pub sigma0: f64,
    /// Evaluation budget.
    pub max_evals: usize,
    /// Target objective value (stop when reached).
    pub f_target: f64,
    /// Generations between eigendecompositions (SPD refresh cadence).
    pub eigen_interval: usize,
}

impl CmaParams {
    /// Standard defaults for dimension `n`.
    #[must_use]
    pub fn standard(n: usize, sigma0: f64, max_evals: usize, f_target: f64) -> CmaParams {
        let lambda = 4 + (3.0 * fs_math::det::ln(n as f64)).floor() as usize;
        CmaParams {
            lambda,
            sigma0,
            max_evals,
            f_target,
            eigen_interval: 1,
        }
    }
}

/// Run evidence.
#[derive(Debug, Clone)]
pub struct CmaReport {
    /// Best point found.
    pub x_best: Vec<f64>,
    /// Best objective value.
    pub f_best: f64,
    /// Objective evaluations consumed.
    pub evals: usize,
    /// Generations run.
    pub generations: usize,
    /// Whether `f_target` was reached.
    pub converged: bool,
    /// Final step size (diagnostic).
    pub sigma: f64,
}

/// Full-covariance CMA-ES from `x0`. Deterministic per `seed`.
#[must_use]
#[allow(clippy::too_many_lines)] // the algorithm is one coherent loop
pub fn cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    params: &CmaParams,
    seed: u64,
) -> CmaReport {
    let n = x0.len();
    assert!(n >= 1, "dimension must be positive");
    let lambda = params.lambda.max(4);
    let mu = lambda / 2;
    // Log-rank recombination weights (Hansen standard).
    let raw: Vec<f64> = (0..mu)
        .map(|i| {
            fs_math::det::ln(f64::midpoint(lambda as f64, 1.0)) - fs_math::det::ln(i as f64 + 1.0)
        })
        .collect();
    let wsum: f64 = raw.iter().sum();
    let weights: Vec<f64> = raw.iter().map(|w| w / wsum).collect();
    let mu_eff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let nf = n as f64;
    // Standard strategy parameters (the IGO/natural-gradient couplings).
    let cc = (4.0 + mu_eff / nf) / (nf + 4.0 + 2.0 * mu_eff / nf);
    let cs = (mu_eff + 2.0) / (nf + mu_eff + 5.0);
    let c1 = 2.0 / ((nf + 1.3) * (nf + 1.3) + mu_eff);
    let cmu =
        (1.0 - c1).min(2.0 * (mu_eff - 2.0 + 1.0 / mu_eff) / ((nf + 2.0) * (nf + 2.0) + mu_eff));
    let damps = 1.0 + 2.0 * (fs_math::det::sqrt((mu_eff - 1.0) / (nf + 1.0)) - 1.0).max(0.0) + cs;
    // E‖N(0,I)‖ (Hansen's approximation).
    let chi_n = fs_math::det::sqrt(nf) * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));

    let mut mean = x0.to_vec();
    let mut sigma = params.sigma0;
    let mut cov = vec![0.0f64; n * n];
    for i in 0..n {
        cov[i * n + i] = 1.0;
    }
    let mut p_c = vec![0.0f64; n];
    let mut p_s = vec![0.0f64; n];
    // Eigen state: C = B·diag(d²)·Bᵀ; sqrt factors refreshed on cadence.
    let mut b_mat = cov.clone();
    let mut d_sqrt = vec![1.0f64; n];
    let mut stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 0,
    }
    .stream();

    let mut x_best = mean.clone();
    let mut f_best = f(&mean);
    let mut evals = 1usize;
    let mut generations = 0usize;
    // TolFun stagnation: generations since a meaningful f_best improvement.
    let mut gens_since_improve = 0usize;

    let mut zs: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut ys: Vec<Vec<f64>> = vec![vec![0.0; n]; lambda];
    let mut fitness: Vec<f64> = vec![0.0; lambda];

    while evals + lambda <= params.max_evals {
        generations += 1;
        // Refresh eigendecomposition on cadence (SPD maintenance).
        if generations % params.eigen_interval.max(1) == 1 || params.eigen_interval <= 1 {
            // Symmetrize (roundoff hygiene) then eigh; floor eigenvalues.
            for i in 0..n {
                for j in i + 1..n {
                    let avg = f64::midpoint(cov[i * n + j], cov[j * n + i]);
                    cov[i * n + j] = avg;
                    cov[j * n + i] = avg;
                }
            }
            let (vals, vecs) = jacobi_eigh(&cov, n);
            let vmax = vals.last().copied().unwrap_or(1.0).max(f64::MIN_POSITIVE);
            for (k, &v) in vals.iter().enumerate() {
                d_sqrt[k] = fs_math::det::sqrt(v.max(1e-14 * vmax));
            }
            b_mat.copy_from_slice(&vecs);
        }
        // Sample λ candidates: x = m + σ·B·diag(d)·z.
        for (k, z) in zs.iter_mut().enumerate() {
            for zi in z.iter_mut() {
                *zi = stream.next_normal();
            }
            let y = &mut ys[k];
            for i in 0..n {
                let mut acc = 0.0f64;
                for j in 0..n {
                    acc = (b_mat[i * n + j] * d_sqrt[j]).mul_add(z[j], acc);
                }
                y[i] = acc;
            }
            let x: Vec<f64> = mean
                .iter()
                .zip(y.iter())
                .map(|(m, yi)| sigma.mul_add(*yi, *m))
                .collect();
            fitness[k] = f(&x);
            evals += 1;
            if fitness[k] < f_best {
                if f_best - fitness[k] > 1e-12 * (1.0 + f_best.abs()) {
                    gens_since_improve = 0;
                }
                f_best = fitness[k];
                x_best = x;
            }
        }
        gens_since_improve += 1;
        if f_best <= params.f_target {
            return CmaReport {
                x_best,
                f_best,
                evals,
                generations,
                converged: true,
                sigma,
            };
        }
        // Rank (total_cmp, lowest index on ties — P2).
        let mut order: Vec<usize> = (0..lambda).collect();
        order.sort_by(|&a, &b| fitness[a].total_cmp(&fitness[b]).then(a.cmp(&b)));
        // Weighted recombination in y-space.
        let mut y_w = vec![0.0f64; n];
        for (w, &idx) in weights.iter().zip(&order) {
            for i in 0..n {
                y_w[i] = w.mul_add(ys[idx][i], y_w[i]);
            }
        }
        // Mean update.
        for i in 0..n {
            mean[i] = sigma.mul_add(y_w[i], mean[i]);
        }
        // CSA path: p_s ← (1−cs)p_s + √(cs(2−cs)µeff)·C^{−1/2}·y_w,
        // with C^{−1/2} = B·diag(1/d)·Bᵀ.
        let mut c_inv_half_yw = vec![0.0f64; n];
        for i in 0..n {
            // t = Bᵀ y_w
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[j * n + i].mul_add(y_w[j], acc);
            }
            c_inv_half_yw[i] = acc / d_sqrt[i];
        }
        let mut tmp = vec![0.0f64; n];
        for i in 0..n {
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = b_mat[i * n + j].mul_add(c_inv_half_yw[j], acc);
            }
            tmp[i] = acc;
        }
        let csn = fs_math::det::sqrt(cs * (2.0 - cs) * mu_eff);
        for i in 0..n {
            p_s[i] = (1.0 - cs).mul_add(p_s[i], csn * tmp[i]);
        }
        let ps_norm = fs_math::det::sqrt(p_s.iter().map(|t| t * t).sum::<f64>());
        // Step-size update (the natural-gradient-consistent coupling).
        sigma *= fs_math::det::exp((cs / damps) * (ps_norm / chi_n - 1.0));
        // STAGNATION STOP: once the search distribution has collapsed
        // (σ·√λmax(C) negligible vs σ₀) the run is dead — keep sampling
        // and it just burns budget polishing whatever basin it's in.
        // BIPOP's restart ladder DEPENDS on dead runs terminating
        // (measured during bring-up: without this, failed runs consumed
        // their entire 120k budget at f ≈ 1 on rastrigin).
        let spread = sigma * d_sqrt.iter().fold(0.0f64, |m, &d| m.max(d));
        if spread < 1e-12 * params.sigma0 || gens_since_improve > 120 {
            // TolX OR TolFun: σ-collapse alone fires too slowly inside a
            // per-run budget (measured: a λ=150 local-basin run burned
            // 120k evals with f stalled for hundreds of generations) —
            // the f-stall criterion is what actually frees the budget.
            break;
        }
        // Rank-1 path with stall indicator h_σ.
        let h_sig = ps_norm
            / fs_math::det::sqrt(
                1.0 - (1.0 - cs)
                    .powi(2 * i32::try_from(generations.min(100_000)).expect("generation count")),
            )
            < (1.4 + 2.0 / (nf + 1.0)) * chi_n;
        let ccn = fs_math::det::sqrt(cc * (2.0 - cc) * mu_eff);
        for i in 0..n {
            let h = if h_sig { ccn * y_w[i] } else { 0.0 };
            p_c[i] = (1.0 - cc).mul_add(p_c[i], h);
        }
        // Covariance update: rank-1 + rank-µ.
        let delta_h = if h_sig { 0.0 } else { cc * (2.0 - cc) };
        for i in 0..n {
            for j in 0..n {
                let mut rank_mu = 0.0f64;
                for (w, &idx) in weights.iter().zip(&order) {
                    rank_mu = (w * ys[idx][i]).mul_add(ys[idx][j], rank_mu);
                }
                let rank1 = p_c[i] * p_c[j];
                cov[i * n + j] = (1.0 - c1 - cmu).mul_add(
                    cov[i * n + j],
                    c1.mul_add(rank1 + delta_h * cov[i * n + j], cmu * rank_mu),
                );
            }
        }
    }
    CmaReport {
        x_best,
        f_best,
        evals,
        generations,
        converged: false,
        sigma,
    }
}

/// BIPOP restart evidence.
#[derive(Debug, Clone)]
pub struct BipopReport {
    /// The best report across all restarts.
    pub best: CmaReport,
    /// Population size used in each restart (the schedule evidence).
    pub schedule: Vec<usize>,
    /// Total evaluations across restarts.
    pub total_evals: usize,
}

/// BIPOP-CMA-ES: alternating large-population (doubling) and
/// small-population restarts under a shared budget — the standard
/// multimodal regime. Deterministic: restart seeds derive from the base
/// seed by counter.
#[must_use]
pub fn bipop_cmaes<F: FnMut(&[f64]) -> f64>(
    f: &mut F,
    x0: &[f64],
    sigma0: f64,
    total_budget: usize,
    f_target: f64,
    seed: u64,
) -> BipopReport {
    let n = x0.len();
    let base_lambda = 4 + (3.0 * fs_math::det::ln(n as f64)).floor() as usize;
    let mut schedule = Vec::new();
    let mut total_evals = 0usize;
    let mut best: Option<CmaReport> = None;
    let mut large_runs = 0u32;
    let mut restart = 0u64;
    let mut small_budget_used = 0usize;
    let mut large_budget_used = 0usize;
    // Deterministic restart-start perturbations (a restart from the SAME
    // point with a tiny sigma is just a polish run and cannot escape a
    // local basin — measured during bring-up on rastrigin).
    let mut restart_stream = StreamKey {
        seed,
        kernel: K_CMA,
        tile: 1,
    }
    .stream();
    while total_evals < total_budget {
        // BIPOP rule: run LARGE next if its cumulative budget lags.
        let run_large = large_budget_used <= small_budget_used;
        let lambda = if run_large {
            base_lambda * (1usize << large_runs)
        } else {
            base_lambda
        };
        // Per-run budget scales with the population (≈250 generations):
        // handing a small-λ run half the TOTAL budget just polishes one
        // local minimum expensively — the doubling ladder must be reached
        // (measured during bring-up on rastrigin).
        let budget = (lambda * 250).min(total_budget - total_evals);
        let params = CmaParams {
            lambda,
            sigma0,
            max_evals: budget,
            f_target,
            eigen_interval: 1,
        };
        // Restarts after the first launch from a perturbed start point.
        let start: Vec<f64> = if restart == 0 {
            x0.to_vec()
        } else {
            x0.iter()
                .map(|&v| sigma0.mul_add(restart_stream.next_normal(), v))
                .collect()
        };
        let rep = cmaes(f, &start, &params, seed.wrapping_add(restart * 0x9E37_79B9));
        total_evals += rep.evals;
        if run_large {
            large_budget_used += rep.evals;
            large_runs += 1;
        } else {
            small_budget_used += rep.evals;
        }
        schedule.push(lambda);
        let is_better = best.as_ref().is_none_or(|b| rep.f_best < b.f_best);
        if is_better {
            best = Some(rep.clone());
        }
        if rep.converged {
            break;
        }
        restart += 1;
        if large_runs > 8 {
            // Cap the LADDER, not total restarts: small runs are cheap
            // and interleave freely; counting them against the cap
            // stalled the ladder at λ ≈ 64 (measured during bring-up).
            break;
        }
    }
    BipopReport {
        best: best.expect("at least one run"),
        schedule,
        total_evals,
    }
}
