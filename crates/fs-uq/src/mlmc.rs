//! Multilevel Monte Carlo: telescoping E[P_L] = E[P₀] + Σ E[P_ℓ −
//! P_{ℓ−1}] over a level ladder, with the optimal sample allocation
//! N_ℓ ∝ √(V_ℓ/C_ℓ) (Giles). The TELESCOPING identity and the
//! variance-per-cost win over single-level MC are AUDITED in the
//! battery, not asserted from theory.

/// Per-level report.
#[derive(Debug, Clone)]
pub struct LevelStats {
    /// Samples taken.
    pub samples: usize,
    /// Mean of the level correction Y_ℓ = P_ℓ − P_{ℓ−1} (Y₀ = P₀).
    pub mean: f64,
    /// Sample variance of Y_ℓ.
    pub variance: f64,
    /// Unit cost supplied for the level.
    pub cost: f64,
}

/// MLMC outcome.
#[derive(Debug, Clone)]
pub struct MlmcReport {
    /// The multilevel estimate Σ mean_ℓ.
    pub estimate: f64,
    /// Per-level statistics (the ledgered evidence).
    pub levels: Vec<LevelStats>,
    /// Total cost spent (Σ N_ℓ·C_ℓ).
    pub total_cost: f64,
    /// Estimator variance (Σ V_ℓ/N_ℓ).
    pub estimator_variance: f64,
}

/// Run MLMC with the Giles allocation for a target estimator
/// variance. `sampler(level, germ_index)` returns the level
/// CORRECTION sample Y_ℓ (callers couple coarse/fine internally —
/// the same germ must drive both, which is what makes V_ℓ decay);
/// `costs[l]` is the unit cost of one level-ℓ sample. A pilot of
/// `pilot` samples per level estimates variances, then the
/// allocation tops up. Deterministic: germ indices are sequential
/// per level.
pub fn mlmc_estimate(
    sampler: &mut dyn FnMut(usize, u64) -> f64,
    costs: &[f64],
    pilot: usize,
    target_variance: f64,
) -> MlmcReport {
    let nl = costs.len();
    let mut sums = vec![(0.0f64, 0.0f64, 0usize); nl]; // (Σy, Σy², n)
    for (l, s) in sums.iter_mut().enumerate() {
        for g in 0..pilot {
            let y = sampler(l, g as u64);
            s.0 += y;
            s.1 = y.mul_add(y, s.1);
            s.2 += 1;
        }
    }
    let var_of = |s: &(f64, f64, usize)| -> f64 {
        let n = s.2 as f64;
        let mean = s.0 / n;
        (s.1 / n - mean * mean).max(1e-30)
    };
    // Giles allocation: N_ℓ = ceil(ε⁻²·√(V_ℓ/C_ℓ)·Σ√(V_ℓC_ℓ)).
    let sum_vc: f64 = sums
        .iter()
        .zip(costs)
        .map(|(s, c)| fs_math::det::sqrt(var_of(s) * c))
        .sum();
    for l in 0..nl {
        let v = var_of(&sums[l]);
        let n_opt = ((fs_math::det::sqrt(v / costs[l]) * sum_vc) / target_variance).ceil() as usize;
        let s = &mut sums[l];
        while s.2 < n_opt {
            let y = sampler(l, s.2 as u64);
            s.0 += y;
            s.1 = y.mul_add(y, s.1);
            s.2 += 1;
        }
    }
    let mut estimate = 0.0f64;
    let mut total_cost = 0.0f64;
    let mut est_var = 0.0f64;
    let levels: Vec<LevelStats> = sums
        .iter()
        .zip(costs)
        .map(|(s, &c)| {
            let n = s.2 as f64;
            let mean = s.0 / n;
            let variance = var_of(s);
            estimate += mean;
            total_cost = (n).mul_add(c, total_cost);
            est_var += variance / n;
            LevelStats {
                samples: s.2,
                mean,
                variance,
                cost: c,
            }
        })
        .collect();
    MlmcReport {
        estimate,
        levels,
        total_cost,
        estimator_variance: est_var,
    }
}
