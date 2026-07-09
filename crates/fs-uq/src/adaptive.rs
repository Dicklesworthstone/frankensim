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
    let mut levels: Vec<AdaptiveLevel> = Vec::new();
    let mut add_level = |l: usize,
                         levels: &mut Vec<AdaptiveLevel>,
                         sampler: &mut dyn FnMut(usize, usize) -> f64| {
        let samples: Vec<f64> = (0..n_pilot).map(|i| sampler(l, i)).collect();
        #[allow(clippy::cast_precision_loss)]
        let mean = samples.iter().sum::<f64>() / n_pilot as f64;
        #[allow(clippy::cast_precision_loss)]
        let var = samples.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>()
            / (n_pilot as f64 - 1.0).max(1.0);
        levels.push(AdaptiveLevel {
            level: l,
            mean,
            var,
            n: n_pilot,
            cost: cost_of(l),
        });
    };
    add_level(0, &mut levels, &mut sampler);
    add_level(1, &mut levels, &mut sampler);
    loop {
        // Fit rates from the correction levels (l >= 1): log2 |mean|
        // and log2 var against level index, least squares.
        let corr: Vec<&AdaptiveLevel> = levels.iter().skip(1).collect();
        let fit_slope = |ys: &[f64]| -> f64 {
            let n = ys.len();
            #[allow(clippy::cast_precision_loss)]
            let xbar = (0..n).map(|i| i as f64).sum::<f64>() / n as f64;
            let ybar = ys.iter().sum::<f64>() / ys.len() as f64;
            let mut num = 0.0;
            let mut den = 0.0;
            for (i, y) in ys.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let dx = i as f64 - xbar;
                num += dx * (y - ybar);
                den += dx * dx;
            }
            if den > 0.0 { -num / den } else { 1.0 }
        };
        let log_means: Vec<f64> = corr
            .iter()
            .map(|l| l.mean.abs().max(1e-300).log2())
            .collect();
        let log_vars: Vec<f64> = corr.iter().map(|l| l.var.max(1e-300).log2()).collect();
        let alpha = fit_slope(&log_means).max(0.1);
        let beta = fit_slope(&log_vars).max(0.1);
        let last = levels.last().expect("levels");
        let bias = last.mean.abs() / (2f64.powf(alpha) - 1.0);
        if bias <= tol / 2.0 || levels.len() > max_level {
            let estimate: f64 = levels.iter().map(|l| l.mean).sum();
            return AdaptiveReport {
                estimate,
                levels,
                alpha,
                beta,
                bias_estimate: bias,
            };
        }
        let next = levels.len();
        add_level(next, &mut levels, &mut sampler);
    }
}
