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
    /// Estimated finest-level expectation. Independent-level MLMC sums the
    /// correction means; a coupled summary retains the original finest values
    /// so cancellation between large level means cannot erase a small result.
    pub estimate: f64,
    /// Per-level statistics (the ledgered evidence).
    pub levels: Vec<LevelStats>,
    /// Total cost spent (Σ N_ℓ·C_ℓ).
    pub total_cost: f64,
    /// Sample estimate of variance of the mean. Independent-level MLMC uses
    /// Σ V_ℓ/N_ℓ; coupled summaries include within-row level covariance by
    /// accumulating the original finest-level values directly. This quantity
    /// is descriptive after data-dependent stopping, not an anytime bound.
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

/// Refusal while summarizing already-observed, coupled resolution levels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoupledMlmcError {
    /// A resolution ladder must contain at least one level.
    EmptyLevels,
    /// No physical observations were supplied.
    EmptySamples,
    /// An original level's unit cost was not finite and strictly positive.
    InvalidCost { level: usize },
    /// An original physical observation was not finite.
    NonFiniteSample { sample: usize, level: usize },
    /// A level correction or its centered moments could not be represented.
    UnrepresentableMoments { level: usize },
    /// The complete observed work could not be represented.
    UnrepresentableCost,
}

impl core::fmt::Display for CoupledMlmcError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled refinement summary: {self:?}")
    }
}
impl std::error::Error for CoupledMlmcError {}

/// Summarize one observed cohort evaluated on a common resolution ladder.
///
/// Each row contains the original physical values [P_0, ..., P_L] for ONE
/// member, with the same excitation and physical horizon at every level.
/// Level statistics describe P_0 and the differences P_l - P_(l-1).
/// Costs are the original per-level solve costs, counting reused solves once.
///
/// Rows are consumed once, without calling a model, allocating more samples,
/// or imposing an independence assumption on level means. The estimate and
/// its sample variance come directly from P_L: this includes all cross-level
/// covariances and preserves a small P_L even when large corrections cancel.
/// A singleton has zero descriptive sample variance, not certainty.
///
/// Work is O(rows * levels), with O(levels) state and one report allocation.
/// These are observed-sample diagnostics. Population uncertainty needs the
/// caller's sampling assumptions; adaptive stopping needs its own confidence
/// sequence. No bias, discretization, or anytime-coverage bound is supplied.
///
/// # Errors
/// Empty levels/samples, nonpositive or nonfinite costs, nonfinite observations,
/// or unrepresentable corrections, moments, or total work.
pub fn coupled_mlmc_summary<const LEVELS: usize>(
    values: impl IntoIterator<Item = [f64; LEVELS]>,
    costs: [f64; LEVELS],
) -> Result<MlmcReport, CoupledMlmcError> {
    if LEVELS == 0 {
        return Err(CoupledMlmcError::EmptyLevels);
    }
    for (level, &cost) in costs.iter().enumerate() {
        if !cost.is_finite() || cost <= 0.0 {
            return Err(CoupledMlmcError::InvalidCost { level });
        }
    }
    let mut statistics = [RunningStats::default(); LEVELS];
    let mut finest = RunningStats::default();
    for (sample, row) in values.into_iter().enumerate() {
        let mut previous = 0.0;
        for (level, (&value, stats)) in row.iter().zip(&mut statistics).enumerate() {
            if !value.is_finite() {
                return Err(CoupledMlmcError::NonFiniteSample { sample, level });
            }
            let correction = value - previous;
            if !correction.is_finite() {
                return Err(CoupledMlmcError::UnrepresentableMoments { level });
            }
            stats.push(correction);
            if !stats.mean.is_finite() || !stats.centered_sum_squares.is_finite() {
                return Err(CoupledMlmcError::UnrepresentableMoments { level });
            }
            previous = value;
        }
        finest.push(previous);
        if !finest.mean.is_finite() || !finest.centered_sum_squares.is_finite() {
            return Err(CoupledMlmcError::UnrepresentableMoments { level: LEVELS - 1 });
        }
    }
    if finest.samples == 0 {
        return Err(CoupledMlmcError::EmptySamples);
    }
    #[allow(clippy::cast_precision_loss)]
    let count = finest.samples as f64;
    let mut total_cost = 0.0;
    let levels = statistics
        .iter()
        .zip(costs)
        .map(|(stats, cost)| {
            total_cost = count.mul_add(cost, total_cost);
            LevelStats {
                samples: stats.samples,
                mean: stats.mean,
                variance: stats.sample_variance(),
                cost,
            }
        })
        .collect();
    if !total_cost.is_finite() {
        return Err(CoupledMlmcError::UnrepresentableCost);
    }
    Ok(MlmcReport {
        estimate: finest.mean,
        levels,
        total_cost,
        estimator_variance: finest.sample_variance() / count,
    })
}

#[cfg(test)]
mod coupled_tests {
    use super::{CoupledMlmcError, coupled_mlmc_summary};

    #[test]
    fn g1_coupled_refinement_includes_cross_level_covariance() {
        let report =
            coupled_mlmc_summary([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0], [4.0, 8.0]], [3.0, 6.0])
                .unwrap();
        assert_eq!(report.estimate, 5.0);
        assert_eq!(report.total_cost, 36.0);
        assert_eq!(report.levels[0].samples, 4);
        assert_eq!(report.levels[1].samples, 4);
        assert!((report.estimator_variance - 5.0 / 3.0).abs() < 1e-14);
        let independent = report
            .levels
            .iter()
            .map(|level| level.variance / 4.0)
            .sum::<f64>();
        assert!((report.estimator_variance - 2.0 * independent).abs() < 1e-14);

        let opposing =
            coupled_mlmc_summary([[100.0, 1.0], [200.0, 2.0], [300.0, 3.0]], [1.0; 2]).unwrap();
        assert_eq!(opposing.estimate, 2.0);
        assert!((opposing.estimator_variance - 1.0 / 3.0).abs() < 1e-14);
    }

    #[test]
    fn g3_coupled_summary_retains_original_finest_values_and_visits_once() {
        let report =
            coupled_mlmc_summary([[1e200, 1.0], [1e200, 2.0], [1e200, 3.0]], [1.0; 2]).unwrap();
        // Rounded coarse and correction means cancel to zero; the supplied
        // physical finest responses still have mean two and variance one.
        assert_eq!(report.levels[0].mean + report.levels[1].mean, 0.0);
        assert_eq!(report.estimate, 2.0);
        assert!((report.estimator_variance - 1.0 / 3.0).abs() < 1e-14);

        let visits = core::cell::Cell::new(0);
        let rows = [3.0, 5.0, 7.0].into_iter().map(|value| {
            visits.set(visits.get() + 1);
            [value]
        });
        let single = coupled_mlmc_summary(rows, [2.0]).unwrap();
        assert_eq!(visits.get(), 3);
        assert_eq!(single.estimate, 5.0);
        assert_eq!(single.levels[0].variance, 4.0);
        assert_eq!(single.estimator_variance, 4.0 / 3.0);
        assert_eq!(single.total_cost, 6.0);
    }

    #[test]
    fn g0_coupled_summary_refuses_invalid_or_unrepresentable_observations() {
        assert!(matches!(
            coupled_mlmc_summary([[0.0; 0]], []),
            Err(CoupledMlmcError::EmptyLevels)
        ));
        assert!(matches!(
            coupled_mlmc_summary(core::iter::empty::<[f64; 1]>(), [1.0]),
            Err(CoupledMlmcError::EmptySamples)
        ));
        assert!(matches!(
            coupled_mlmc_summary([[1.0]], [f64::INFINITY]),
            Err(CoupledMlmcError::InvalidCost { level: 0 })
        ));
        assert!(matches!(
            coupled_mlmc_summary([[1.0, f64::NAN]], [1.0; 2]),
            Err(CoupledMlmcError::NonFiniteSample {
                sample: 0,
                level: 1
            })
        ));
        assert!(matches!(
            coupled_mlmc_summary([[-f64::MAX, f64::MAX]], [1.0; 2]),
            Err(CoupledMlmcError::UnrepresentableMoments { level: 1 })
        ));
        assert!(matches!(
            coupled_mlmc_summary([[0.0], [0.0]], [f64::MAX]),
            Err(CoupledMlmcError::UnrepresentableCost)
        ));
    }
}
