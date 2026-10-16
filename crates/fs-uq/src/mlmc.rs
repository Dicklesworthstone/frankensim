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
    /// Bessel-corrected sample variance of Y_ℓ; zero for a singleton.
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
    /// Estimator variance (Σ V_ℓ/N_ℓ), using the Bessel-corrected
    /// per-level sample variances.
    pub estimator_variance: f64,
}

/// Numerically stable one-pass moments for a level correction.
///
/// Keeping the centered sum of squares avoids the catastrophic cancellation in
/// `E[Y²] - E[Y]²` when corrections have a large common offset but small
/// spread. The reported variance uses Bessel's correction, matching
/// `adaptive_mlmc`'s definition of a sample variance.
#[derive(Debug, Clone, Copy, Default)]
struct RunningStats {
    samples: usize,
    mean: f64,
    centered_sum_squares: f64,
}

impl RunningStats {
    fn push(&mut self, sample: f64) {
        self.samples += 1;
        #[allow(clippy::cast_precision_loss)]
        let count = self.samples as f64;
        let delta = sample - self.mean;
        self.mean += delta / count;
        let centered_delta = sample - self.mean;
        self.centered_sum_squares = delta.mul_add(centered_delta, self.centered_sum_squares);
    }

    fn sample_variance(self) -> f64 {
        if self.samples < 2 {
            return 0.0;
        }
        #[allow(clippy::cast_precision_loss)]
        let degrees_of_freedom = (self.samples - 1) as f64;
        // Roundoff can only make a mathematically non-negative M2 slightly
        // negative. Clamp that representation artifact, not the uncertainty.
        (self.centered_sum_squares / degrees_of_freedom).max(0.0)
    }

    fn allocation_variance(self) -> f64 {
        // A zero pilot variance must not remove a level from the allocation:
        // later samples may vary. This floor is an allocation safeguard only;
        // reports retain the honest zero sample variance.
        self.sample_variance().max(1e-30)
    }
}

/// Run MLMC with the Giles allocation for a target estimator
/// variance. `sampler(level, germ_index)` returns the level
/// CORRECTION sample Y_ℓ (callers couple coarse/fine internally —
/// the same germ must drive both, which is what makes V_ℓ decay);
/// `costs[l]` is the unit cost of one level-ℓ sample. A pilot of
/// `pilot` samples per level estimates variances, then the
/// allocation tops up. Per-level variance is accumulated with stable centered
/// moments and Bessel's correction; a singleton level has zero reported
/// variance. Deterministic: germ indices are sequential per level.
pub fn mlmc_estimate(
    sampler: &mut dyn FnMut(usize, u64) -> f64,
    costs: &[f64],
    pilot: usize,
    target_variance: f64,
) -> MlmcReport {
    // Fail closed on degenerate inputs (mirrors `adaptive_mlmc`'s pilot guard).
    // pilot == 0 supplies neither a level mean nor variance evidence; treating
    // it like the deliberately defined singleton case would report a
    // fake-confident zero estimator variance. A non-positive `target_variance`
    // or `costs[l]` drives `v / costs[l]` or `.../ target_variance` to +∞, so
    // `n_opt = f64::INFINITY as usize = usize::MAX` and the top-up loop samples
    // essentially forever (a hang, not a wrong number). NaN costs are rejected
    // by the same `> 0.0` test.
    assert!(
        pilot > 0,
        "mlmc_estimate needs a nonzero pilot sample count"
    );
    assert!(
        target_variance > 0.0,
        "mlmc_estimate needs a positive target variance"
    );
    assert!(
        costs.iter().all(|&c| c > 0.0),
        "mlmc_estimate needs strictly positive per-level costs"
    );
    let nl = costs.len();
    let mut statistics = vec![RunningStats::default(); nl];
    for (l, stats) in statistics.iter_mut().enumerate() {
        for g in 0..pilot {
            stats.push(sampler(l, g as u64));
        }
    }
    // Giles allocation: N_ℓ = ceil(ε⁻²·√(V_ℓ/C_ℓ)·Σ√(V_ℓC_ℓ)).
    let sum_vc: f64 = statistics
        .iter()
        .zip(costs)
        .map(|(stats, c)| fs_math::det::sqrt(stats.allocation_variance() * c))
        .sum();
    for l in 0..nl {
        let v = statistics[l].allocation_variance();
        let n_opt = ((fs_math::det::sqrt(v / costs[l]) * sum_vc) / target_variance).ceil() as usize;
        let stats = &mut statistics[l];
        while stats.samples < n_opt {
            stats.push(sampler(l, stats.samples as u64));
        }
    }
    let mut estimate = 0.0f64;
    let mut total_cost = 0.0f64;
    let mut est_var = 0.0f64;
    let levels: Vec<LevelStats> = statistics
        .iter()
        .zip(costs)
        .map(|(stats, &c)| {
            #[allow(clippy::cast_precision_loss)]
            let n = stats.samples as f64;
            let variance = stats.sample_variance();
            estimate += stats.mean;
            total_cost = (n).mul_add(c, total_cost);
            est_var += variance / n;
            LevelStats {
                samples: stats.samples,
                mean: stats.mean,
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
