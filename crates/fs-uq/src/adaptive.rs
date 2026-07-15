//! ADAPTIVE MLMC (bead o5kc): level addition driven by BIAS estimates
//! — the slice-1 ladder was caller-fixed; this one grows itself.
//! Weak (α) and strong (β) convergence rates are ESTIMATED from level
//! statistics, the remaining bias is extrapolated as
//! `|mean_L| / (2^α − 1)`, and levels are added until the bias fits
//! inside half the tolerance (the other half goes to variance via the
//! standard optimal sample allocation `n_l ∝ √(V_l / C_l)`).

/// One level's running statistics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveLevel {
    /// Level index.
    pub level: usize,
    /// Correction mean `E[Y_l]`.
    pub mean: f64,
    /// Correction variance `V[Y_l]`.
    pub var: f64,
    /// Samples taken.
    pub n: usize,
    /// Unit cost (caller-declared, e.g. cells).
    pub cost: f64,
}

/// The adaptive-MLMC outcome: the estimate, the audited level table,
/// and the fitted rates.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveReport {
    /// The telescoped estimate.
    pub estimate: f64,
    /// Level statistics (the audit trail).
    pub levels: Vec<AdaptiveLevel>,
    /// Fitted weak rate α (mean decay per level, log2).
    pub alpha: f64,
    /// Fitted strong rate β (variance decay per level, log2).
    pub beta: f64,
    /// The extrapolated remaining bias at the stop.
    pub bias_estimate: f64,
    /// Audited estimator variance `sum_l V_l / N_l` after allocation.
    pub estimator_variance: f64,
}

/// Stable centered moments for one correction level.
#[derive(Debug, Clone, Copy, Default)]
struct RunningStats {
    samples: usize,
    mean: f64,
    centered_sum_squares: f64,
}

impl RunningStats {
    fn push(&mut self, sample: f64) {
        assert!(sample.is_finite(), "adaptive MLMC samples must be finite");
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
        (self.centered_sum_squares / degrees_of_freedom).max(0.0)
    }

    fn allocation_variance(self) -> f64 {
        // A zero pilot variance is not evidence that all future corrections
        // are constant. Keep the safeguard out of the reported variance.
        self.sample_variance().max(1e-30)
    }
}

fn fitted_rates(statistics: &[RunningStats]) -> (f64, f64) {
    let corr = &statistics[1..];
    let fit_slope = |ys: &[f64]| -> f64 {
        let n = ys.len();
        #[allow(clippy::cast_precision_loss)]
        let xbar = (0..n).map(|i| i as f64).sum::<f64>() / n as f64;
        let ybar = ys.iter().sum::<f64>() / ys.len() as f64;
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for (i, y) in ys.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let dx = i as f64 - xbar;
            numerator += dx * (y - ybar);
            denominator += dx * dx;
        }
        if denominator > 0.0 {
            -numerator / denominator
        } else {
            1.0
        }
    };
    let log_means: Vec<f64> = corr
        .iter()
        .map(|stats| stats.mean.abs().max(1e-300).log2())
        .collect();
    let log_variances: Vec<f64> = corr
        .iter()
        .map(|stats| stats.sample_variance().max(1e-300).log2())
        .collect();
    (
        fit_slope(&log_means).max(0.1),
        fit_slope(&log_variances).max(0.1),
    )
}

fn allocation_targets(
    statistics: &[RunningStats],
    costs: &[f64],
    target_variance: f64,
) -> Vec<usize> {
    let sum_sqrt_vc: f64 = statistics
        .iter()
        .zip(costs)
        .map(|(stats, &cost)| fs_math::det::sqrt(stats.allocation_variance() * cost))
        .sum();
    assert!(
        sum_sqrt_vc.is_finite() && sum_sqrt_vc > 0.0,
        "adaptive MLMC allocation weights must be finite"
    );
    statistics
        .iter()
        .zip(costs)
        .map(|(stats, &cost)| {
            let raw = fs_math::det::sqrt(stats.allocation_variance() / cost) * sum_sqrt_vc
                / target_variance;
            #[allow(clippy::cast_precision_loss)]
            let representable = raw.is_finite() && raw <= usize::MAX as f64;
            assert!(
                representable,
                "adaptive MLMC sample allocation exceeds usize"
            );
            (raw.ceil() as usize).max(stats.samples)
        })
        .collect()
}

/// Run adaptive MLMC: `sampler(level, i)` returns the level-`l`
/// correction sample `Y_l(ω_i)` (level 0 = the coarse value itself;
/// ONE germ per index drives both rungs of a correction, the slice-1
/// contract). `cost_of(level)` prices a sample. Levels are added while
/// the extrapolated bias exceeds `tol / 2`, up to `max_level`.
pub fn adaptive_mlmc(
    mut sampler: impl FnMut(usize, usize) -> f64,
    cost_of: impl Fn(usize) -> f64,
    tol: f64,
    n_pilot: usize,
    max_level: usize,
) -> AdaptiveReport {
    assert!(
        tol.is_finite() && tol > 0.0,
        "tol must be positive and finite"
    );
    assert!(n_pilot > 0, "n_pilot must be nonzero");
    assert!(
        max_level >= 1,
        "adaptive MLMC needs at least levels 0 and 1"
    );
    let variance_radius = tol / 2.0;
    let target_variance = variance_radius * variance_radius;
    assert!(
        target_variance.is_finite() && target_variance > 0.0,
        "tol is too small to represent an adaptive MLMC variance target"
    );

    let mut statistics: Vec<RunningStats> = Vec::new();
    let mut costs = Vec::new();
    let add_level = |l: usize,
                     statistics: &mut Vec<RunningStats>,
                     costs: &mut Vec<f64>,
                     sampler: &mut dyn FnMut(usize, usize) -> f64| {
        let cost = cost_of(l);
        assert!(
            cost.is_finite() && cost > 0.0,
            "adaptive MLMC costs must be positive and finite"
        );
        let mut stats = RunningStats::default();
        for i in 0..n_pilot {
            stats.push(sampler(l, i));
        }
        statistics.push(stats);
        costs.push(cost);
    };
    add_level(0, &mut statistics, &mut costs, &mut sampler);
    add_level(1, &mut statistics, &mut costs, &mut sampler);

    loop {
        let (alpha, _) = fitted_rates(&statistics);
        let bias = statistics.last().expect("levels").mean.abs() / (2f64.powf(alpha) - 1.0);
        if bias > variance_radius && statistics.len() <= max_level {
            let next = statistics.len();
            add_level(next, &mut statistics, &mut costs, &mut sampler);
            continue;
        }

        // Re-estimate after each top-up. A changed sample variance can increase
        // a target, so stop only when the Giles allocation is a fixed point.
        // Fail closed rather than return under-allocated evidence if an
        // adversarial, nonstationary sampler prevents stabilization.
        let mut top_up_passes = 0usize;
        loop {
            let targets = allocation_targets(&statistics, &costs, target_variance);
            if targets
                .iter()
                .zip(&statistics)
                .all(|(&target, stats)| target <= stats.samples)
            {
                break;
            }
            top_up_passes += 1;
            assert!(
                top_up_passes <= 64,
                "adaptive MLMC allocation did not stabilize"
            );
            for (level, target) in targets.into_iter().enumerate() {
                while statistics[level].samples < target {
                    let sample_index = statistics[level].samples;
                    statistics[level].push(sampler(level, sample_index));
                }
            }
        }

        // Top-ups change fitted means as well as variances. Re-open the bias
        // loop when capacity remains so the final bias evidence describes the
        // same samples used by the returned estimate.
        let (alpha, beta) = fitted_rates(&statistics);
        let bias = statistics.last().expect("levels").mean.abs() / (2f64.powf(alpha) - 1.0);
        if bias > variance_radius && statistics.len() <= max_level {
            let next = statistics.len();
            add_level(next, &mut statistics, &mut costs, &mut sampler);
            continue;
        }

        let estimator_variance: f64 = statistics
            .iter()
            .map(|stats| {
                #[allow(clippy::cast_precision_loss)]
                let samples = stats.samples as f64;
                stats.sample_variance() / samples
            })
            .sum();
        assert!(
            estimator_variance <= target_variance * (1.0 + 64.0 * f64::EPSILON),
            "adaptive MLMC allocation failed its variance target"
        );
        let levels: Vec<AdaptiveLevel> = statistics
            .iter()
            .zip(&costs)
            .enumerate()
            .map(|(level, (stats, &cost))| AdaptiveLevel {
                level,
                mean: stats.mean,
                var: stats.sample_variance(),
                n: stats.samples,
                cost,
            })
            .collect();
        let estimate = levels.iter().map(|level| level.mean).sum();
        return AdaptiveReport {
            estimate,
            levels,
            alpha,
            beta,
            bias_estimate: bias,
            estimator_variance,
        };
    }
}
