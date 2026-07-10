//! fs-eproc — anytime-valid inference (plan §9.6, Bet 5): betting
//! e-processes, mixture confidence sequences, e-value arithmetic, and e-BH.
//!
//! # Why this exists
//! Classical tests are valid only at a FIXED sample size; an optimizer that
//! peeks and stops when it likes destroys their guarantees. An e-process is
//! valid at EVERY stopping time simultaneously (Ville's inequality:
//! P(∃t: E_t ≥ 1/α) ≤ α under the null), so ASCENT may monitor continuously
//! and kill candidates the moment evidence crosses threshold — the
//! statistical heart of e-racing and anytime-valid UQ stopping.
//!
//! # Validity is structural
//! The betting construction multiplies wealth by (1 + λ_t·(x_t − m)) with a
//! PREDICTABLE λ_t in the admissible range — under H₀ (mean m), wealth is a
//! nonnegative supermartingale REGARDLESS of the betting strategy. Strategy
//! choice affects POWER only. Our plug-in strategy follows the
//! Waudby-Smith–Ramdas lineage (empirical mean/variance, clipped).
//!
//! # Determinism
//! e-trajectories use fs-math strict functions only — combined with fs-rand
//! logical-identity streams, every tournament is bit-replayable from its
//! seed (Bet 8's reproducible-racing requirement).

#[cfg(feature = "conformal-hardening")]
pub mod hardening;

use core::fmt;

use fs_math::det;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ---------------------------------------------------------------------------
// Betting e-process for bounded outcomes.
// ---------------------------------------------------------------------------

/// One-sided betting e-process testing H₀: mean(X) ≤ `null_mean` against
/// "mean is LARGER", for outcomes in [0, 1]. (Race two candidates by feeding
/// their declared-span-normalized difference with null mean 1/2 — see
/// [`PairwiseRace`].)
#[derive(Debug, Clone)]
pub struct BettingEProcess {
    null_mean: f64,
    /// log of the wealth (e-value) — log-space avoids overflow at huge e.
    log_wealth: f64,
    /// Running count and moments for the predictable plug-in bet.
    n: u64,
    sum: f64,
    sum_sq: f64,
    /// Fraction of the maximal admissible bet to use (default 0.5 — the
    /// "half-Kelly"-style hedge that trades a little power for robustness).
    aggressiveness: f64,
}

impl BettingEProcess {
    /// Create a fresh process for H₀: mean ≤ `null_mean`, outcomes ∈ [0,1].
    ///
    /// # Panics
    /// If `null_mean` is not strictly inside (0, 1).
    #[must_use]
    pub fn new(null_mean: f64) -> Self {
        assert!(
            null_mean > 0.0 && null_mean < 1.0,
            "null mean must lie strictly inside (0,1); got {null_mean}"
        );
        BettingEProcess {
            null_mean,
            log_wealth: 0.0,
            n: 0,
            sum: 0.0,
            sum_sq: 0.0,
            aggressiveness: 0.5,
        }
    }

    /// The predictable bet for the NEXT observation (computed from data seen
    /// so far only — predictability is what makes validity structural).
    fn next_lambda(&self) -> f64 {
        // Plug-in: bet proportional to (μ̂ − m)/σ̂², clipped inside the
        // admissible range scaled by `aggressiveness`. Regularized moments
        // (add a pseudo-observation at the null) keep early bets tame.
        let n = self.n as f64;
        let mu = (self.sum + self.null_mean) / (n + 1.0);
        let var = ((self.sum_sq + self.null_mean * self.null_mean) / (n + 1.0) - mu * mu).max(1e-4);
        let raw = (mu - self.null_mean) / var;
        let cap = self.aggressiveness / self.null_mean.max(1.0 - self.null_mean);
        raw.clamp(0.0, cap) // one-sided: never bet on "smaller"
    }

    /// Observe one outcome; returns the updated log e-value.
    ///
    /// # Panics
    /// If `x` is outside [0, 1] (the boundedness the guarantee rests on —
    /// feeding unbounded data here would VOID validity, so it is refused).
    pub fn observe(&mut self, x: f64) -> f64 {
        assert!(
            (0.0..=1.0).contains(&x),
            "outcome {x} outside [0,1] voids the guarantee"
        );
        let lambda = self.next_lambda();
        // Wealth *= 1 + λ(x − m); log-space via strict ln (argument is
        // ≥ 1 − aggressiveness > 0 by the λ cap, so ln is safe).
        self.log_wealth += det::ln(lambda.mul_add(x - self.null_mean, 1.0));
        self.n += 1;
        self.sum += x;
        self.sum_sq += x * x;
        self.log_wealth
    }

    /// Current e-value (wealth). May be +∞-adjacent for huge evidence; use
    /// [`Self::log_e_value`] in thresholds.
    #[must_use]
    pub fn e_value(&self) -> f64 {
        det::exp(self.log_wealth)
    }

    /// Current log e-value.
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.log_wealth
    }

    /// Has evidence crossed `1/alpha` (the Ville threshold for level α)?
    #[must_use]
    pub fn rejects_at(&self, alpha: f64) -> bool {
        self.log_wealth >= -det::ln(alpha)
    }

    /// Observations consumed.
    #[must_use]
    pub fn len(&self) -> u64 {
        self.n
    }

    /// True before any observation.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

// ---------------------------------------------------------------------------
// Pairwise racing (the e-racing primitive).
// ---------------------------------------------------------------------------

/// Race candidate A against candidate B on paired noisy scores where LOWER
/// is better (losses). A declared support `s` maps the raw difference to
/// `d = ((b - a) / s + 1) / 2` in `[0, 1]`, feeding a betting e-process
/// with null mean 1/2: evidence accumulates that A beats B. Out-of-support
/// observations are refused because clipping would change the estimand.
#[derive(Debug, Clone)]
pub struct PairwiseRace {
    proc: BettingEProcess,
    loss_span: LossSpan,
}

/// Finite positive support bound for an absolute paired-loss difference.
/// Construction is checked so a race cannot carry a malformed scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LossSpan(f64);

impl Eq for LossSpan {}

impl LossSpan {
    /// Unit span for already-normalized losses.
    pub const ONE: Self = Self(1.0);

    /// Validate a raw span.
    ///
    /// # Errors
    /// [`PairwiseInputError::InvalidLossSpan`] unless `span` is finite
    /// and strictly positive.
    pub fn new(span: f64) -> Result<Self, PairwiseInputError> {
        if !span.is_finite() || span <= 0.0 {
            return Err(PairwiseInputError::InvalidLossSpan {
                span_bits: span.to_bits(),
            });
        }
        Ok(Self(span))
    }

    /// The checked raw span.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// A paired-loss observation that cannot support the claimed e-process.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairwiseInputError {
    /// The declared maximum absolute paired-loss difference is invalid.
    InvalidLossSpan {
        /// IEEE-754 bits of the rejected span.
        span_bits: u64,
    },
    /// At least one loss is non-finite.
    NonFiniteLoss {
        /// IEEE-754 bits of candidate A's loss.
        loss_a_bits: u64,
        /// IEEE-754 bits of candidate B's loss.
        loss_b_bits: u64,
    },
    /// Subtracting two finite losses overflowed.
    NonFiniteDifference {
        /// IEEE-754 bits of candidate A's loss.
        loss_a_bits: u64,
        /// IEEE-754 bits of candidate B's loss.
        loss_b_bits: u64,
    },
    /// The observed paired difference exceeds its declared support.
    DifferenceOutOfRange {
        /// IEEE-754 bits of `loss_b - loss_a`.
        difference_bits: u64,
        /// IEEE-754 bits of the checked declared span.
        span_bits: u64,
    },
}

impl fmt::Display for PairwiseInputError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            PairwiseInputError::InvalidLossSpan { span_bits } => write!(
                f,
                "pairwise loss span must be finite and positive; got {}",
                f64::from_bits(span_bits)
            ),
            PairwiseInputError::NonFiniteLoss {
                loss_a_bits,
                loss_b_bits,
            } => write!(
                f,
                "pairwise losses must be finite; got ({}, {})",
                f64::from_bits(loss_a_bits),
                f64::from_bits(loss_b_bits)
            ),
            PairwiseInputError::NonFiniteDifference {
                loss_a_bits,
                loss_b_bits,
            } => write!(
                f,
                "finite pairwise losses overflowed during subtraction: ({}, {})",
                f64::from_bits(loss_a_bits),
                f64::from_bits(loss_b_bits)
            ),
            PairwiseInputError::DifferenceOutOfRange {
                difference_bits,
                span_bits,
            } => write!(
                f,
                "paired-loss difference {} exceeds declared span {}",
                f64::from_bits(difference_bits),
                f64::from_bits(span_bits)
            ),
        }
    }
}

impl std::error::Error for PairwiseInputError {}

impl PairwiseRace {
    /// Fresh race with a fixed, checked paired-loss support.
    #[must_use]
    pub fn new(loss_span: LossSpan) -> Self {
        PairwiseRace {
            proc: BettingEProcess::new(0.5),
            loss_span,
        }
    }

    /// Observe one paired `(loss_a, loss_b)` evaluation. The declared
    /// span maps the raw difference linearly to `[0, 1]`; values outside
    /// that support are refused rather than clipped into a different
    /// estimand.
    ///
    /// # Errors
    /// Non-finite losses, subtraction overflow, or a paired difference
    /// outside the declared span. On error, wealth is unchanged.
    pub fn observe(&mut self, loss_a: f64, loss_b: f64) -> Result<(), PairwiseInputError> {
        if !loss_a.is_finite() || !loss_b.is_finite() {
            return Err(PairwiseInputError::NonFiniteLoss {
                loss_a_bits: loss_a.to_bits(),
                loss_b_bits: loss_b.to_bits(),
            });
        }
        let difference = loss_b - loss_a;
        if !difference.is_finite() {
            return Err(PairwiseInputError::NonFiniteDifference {
                loss_a_bits: loss_a.to_bits(),
                loss_b_bits: loss_b.to_bits(),
            });
        }
        if difference.abs() > self.loss_span.get() {
            return Err(PairwiseInputError::DifferenceOutOfRange {
                difference_bits: difference.to_bits(),
                span_bits: self.loss_span.get().to_bits(),
            });
        }
        let d = f64::midpoint(difference / self.loss_span.get(), 1.0);
        let _ = self.proc.observe(d);
        Ok(())
    }

    /// Declared maximum absolute paired-loss difference.
    #[must_use]
    pub fn loss_span(&self) -> LossSpan {
        self.loss_span
    }

    /// Does the race declare "A beats B" at level α?
    #[must_use]
    pub fn a_beats_b(&self, alpha: f64) -> bool {
        self.proc.rejects_at(alpha)
    }

    /// Current log e-value of "A beats B".
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.proc.log_e_value()
    }
}

// ---------------------------------------------------------------------------
// Gaussian (sub-Gaussian) mixture confidence sequence.
// ---------------------------------------------------------------------------

/// Robbins' normal-mixture confidence sequence for the mean of sub-Gaussian
/// observations with parameter `sigma`: time-uniform coverage
/// P(∀t: μ ∈ CS_t) ≥ 1 − α. The closed-form radius is
/// √( (tσ² + ρ)/t² · ( ln((tσ² + ρ)/ρ) + 2 ln(1/α) ) ).
#[derive(Debug, Clone)]
pub struct GaussianMixtureCs {
    sigma: f64,
    rho: f64,
    alpha: f64,
    n: u64,
    sum: f64,
}

impl GaussianMixtureCs {
    /// New CS at level `alpha` for sub-Gaussian-`sigma` data; `rho > 0` tunes
    /// WHERE the boundary is tightest (≈ the sample size you care about ×σ²).
    ///
    /// # Panics
    /// On non-positive `sigma`/`rho` or `alpha` outside (0,1).
    #[must_use]
    pub fn new(sigma: f64, rho: f64, alpha: f64) -> Self {
        assert!(sigma > 0.0 && rho > 0.0, "sigma/rho must be positive");
        assert!(alpha > 0.0 && alpha < 1.0, "alpha in (0,1)");
        GaussianMixtureCs {
            sigma,
            rho,
            alpha,
            n: 0,
            sum: 0.0,
        }
    }

    /// Observe a value.
    pub fn observe(&mut self, x: f64) {
        self.n += 1;
        self.sum += x;
    }

    /// Current interval (center, radius); `None` before any data.
    #[must_use]
    pub fn interval(&self) -> Option<(f64, f64)> {
        if self.n == 0 {
            return None;
        }
        let t = self.n as f64;
        let v = t * self.sigma * self.sigma + self.rho;
        let radius = det::sqrt(v * (det::ln(v / self.rho) + 2.0 * det::ln(1.0 / self.alpha))) / t;
        Some((self.sum / t, radius))
    }

    /// The two-sided e-value for H₀: mean = `m` (the mixture martingale
    /// itself — usable in e-value arithmetic).
    #[must_use]
    pub fn e_value_for(&self, m: f64) -> f64 {
        if self.n == 0 {
            return 1.0;
        }
        let t = self.n as f64;
        let s = self.sum - m * t;
        let v = t * self.sigma * self.sigma + self.rho;
        det::sqrt(self.rho / v) * det::exp(s * s / (2.0 * v))
    }

    /// Observations consumed.
    #[must_use]
    pub fn len(&self) -> u64 {
        self.n
    }

    /// True before any observation.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

// ---------------------------------------------------------------------------
// e-value arithmetic and e-BH.
// ---------------------------------------------------------------------------

/// Combine INDEPENDENT e-values by product (evidence multiplies).
#[must_use]
pub fn combine_product(log_e: &[f64]) -> f64 {
    log_e.iter().sum()
}

/// Combine ARBITRARILY DEPENDENT e-values by averaging (always valid —
/// the mixture is an e-value whatever the dependence). Input/output in
/// log space; computed with a max-shift for stability.
#[must_use]
pub fn combine_average(log_e: &[f64]) -> f64 {
    if log_e.is_empty() {
        return 0.0;
    }
    let m = log_e.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let sum_exp: f64 = log_e.iter().map(|&l| det::exp(l - m)).sum();
    m + det::ln(sum_exp / log_e.len() as f64)
}

/// An e-value's implied p-value bound: p ≤ 1/e (Markov).
#[must_use]
pub fn e_to_p(log_e: f64) -> f64 {
    det::exp(-log_e).min(1.0)
}

/// e-BH (Wang–Ramdas): given m log-e-values, control FDR at level α by
/// rejecting the k̂ hypotheses with largest e-values, where
/// k̂ = max{k : e_(k) ≥ m/(α·k)} (e_(k) the k-th LARGEST). Valid under
/// ARBITRARY dependence between the e-values. Returns rejected indices.
#[must_use]
pub fn e_benjamini_hochberg(log_e: &[f64], alpha: f64) -> Vec<usize> {
    let m = log_e.len();
    if m == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..m).collect();
    // Deterministic tie-breaking: by (descending e, ascending index).
    order.sort_by(|&a, &b| {
        log_e[b]
            .partial_cmp(&log_e[a])
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(a.cmp(&b))
    });
    let mut k_hat = 0;
    for (rank0, &idx) in order.iter().enumerate() {
        let k = rank0 + 1;
        let threshold = det::ln(m as f64 / (alpha * k as f64));
        if log_e[idx] >= threshold {
            k_hat = k;
        }
    }
    order.truncate(k_hat);
    order.sort_unstable();
    order
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_rand::StreamKey;

    fn key(tile: u32) -> StreamKey {
        StreamKey {
            seed: 0xE9_0C ^ 0x5EED,
            kernel: 1,
            tile,
        }
    }

    /// THE validity test: under the null, even an ADVERSARY who stops at the
    /// supremum of wealth over the whole horizon must cross 1/α with
    /// probability ≤ α (Ville). This is what "anytime-valid" means.
    #[test]
    fn ville_validity_under_adversarial_stopping() {
        const SIMS: u32 = 4_000;
        const HORIZON: usize = 2_000;
        let alpha = 0.05;
        let mut crossings = 0u32;
        for sim in 0..SIMS {
            let mut s = key(sim).stream();
            let mut e = BettingEProcess::new(0.5);
            let mut crossed = false;
            for _ in 0..HORIZON {
                let x = s.next_f64(); // true mean 0.5 == null
                let _ = e.observe(x);
                if e.rejects_at(alpha) {
                    crossed = true;
                    break; // adversary stops the instant it looks significant
                }
            }
            crossings += u32::from(crossed);
        }
        let rate = f64::from(crossings) / f64::from(SIMS);
        // Binomial slack: sd ≈ sqrt(.05*.95/4000) ≈ 0.0034; allow 4σ.
        assert!(
            rate <= alpha + 0.014,
            "adversarial-stopping type-I rate {rate} exceeds alpha {alpha} — validity broken"
        );
        println!(
            "{{\"suite\":\"fs-eproc\",\"case\":\"ville-validity\",\"verdict\":\"pass\",\"detail\":\"adversarial type-I {rate:.4} <= {alpha}\"}}"
        );
    }

    #[test]
    fn power_detects_true_effects_and_scales_with_size() {
        // Shifted Bernoulli-ish outcomes: mean 0.5 + delta must be detected,
        // faster for larger delta.
        let alpha = 0.05;
        let mut stop_times = Vec::new();
        for (tile, delta) in [(1_000u32, 0.05f64), (2_000, 0.15)] {
            let mut times = Vec::new();
            for rep in 0..200 {
                let mut s = key(tile + rep).stream();
                let mut e = BettingEProcess::new(0.5);
                let mut t = 0u64;
                while !e.rejects_at(alpha) && t < 100_000 {
                    let x = (s.next_f64() + delta).clamp(0.0, 1.0);
                    let _ = e.observe(x);
                    t += 1;
                }
                times.push(t);
            }
            times.sort_unstable();
            stop_times.push(times[times.len() / 2]);
        }
        assert!(stop_times[0] < 100_000, "delta=0.05 never detected");
        assert!(
            stop_times[1] * 4 < stop_times[0],
            "bigger effects must stop much sooner: medians {stop_times:?}"
        );
        println!(
            "{{\"suite\":\"fs-eproc\",\"case\":\"power\",\"verdict\":\"pass\",\"detail\":\"median stop times {stop_times:?} for deltas [0.05, 0.15]\"}}"
        );
    }

    #[test]
    fn cs_time_uniform_coverage_and_shrinkage() {
        const SIMS: u32 = 1_500;
        const HORIZON: usize = 1_500;
        let (sigma, rho, alpha) = (1.0, 10.0, 0.05);
        let true_mean = 0.7;
        let mut ever_missed = 0u32;
        for sim in 0..SIMS {
            let mut s = key(50_000 + sim).stream();
            let mut cs = GaussianMixtureCs::new(sigma, rho, alpha);
            let mut missed = false;
            for _ in 0..HORIZON {
                cs.observe(true_mean + s.next_normal() * sigma);
                let (c, r) = cs.interval().expect("has data");
                if (c - true_mean).abs() > r {
                    missed = true;
                    break;
                }
            }
            ever_missed += u32::from(missed);
        }
        let miss = f64::from(ever_missed) / f64::from(SIMS);
        assert!(
            miss <= alpha + 0.02,
            "time-uniform miss rate {miss} > {alpha}"
        );
        // Shrinkage: radius at t=1500 well below radius at t=30.
        let mut s = key(99_999).stream();
        let mut cs = GaussianMixtureCs::new(sigma, rho, alpha);
        let mut r30 = f64::NAN;
        for t in 1..=1_500 {
            cs.observe(true_mean + s.next_normal());
            if t == 30 {
                r30 = cs.interval().expect("data").1;
            }
        }
        let r_end = cs.interval().expect("data").1;
        assert!(r_end < r30 / 4.0, "radius must shrink: {r30} -> {r_end}");
        println!(
            "{{\"suite\":\"fs-eproc\",\"case\":\"cs-coverage\",\"verdict\":\"pass\",\"detail\":\"miss {miss:.4}; radius {r30:.3}->{r_end:.3}\"}}"
        );
    }

    #[test]
    fn e_bh_controls_fdr_and_finds_signals() {
        // 40 hypotheses: 10 true effects (delta .2), 30 nulls; run to a fixed
        // horizon, apply e-BH, measure FDR over sims.
        const SIMS: u32 = 300;
        let alpha = 0.1;
        let mut fdp_sum = 0.0;
        let mut power_sum = 0.0;
        for sim in 0..SIMS {
            let mut log_es = Vec::new();
            for h in 0..40u32 {
                let delta = if h < 10 { 0.2 } else { 0.0 };
                let mut s = key(200_000 + sim * 64 + h).stream();
                let mut e = BettingEProcess::new(0.5);
                for _ in 0..400 {
                    let _ = e.observe((s.next_f64() + delta).clamp(0.0, 1.0));
                }
                log_es.push(e.log_e_value());
            }
            let rejected = e_benjamini_hochberg(&log_es, alpha);
            let false_r = rejected.iter().filter(|&&i| i >= 10).count() as f64;
            let true_r = rejected.iter().filter(|&&i| i < 10).count() as f64;
            fdp_sum += if rejected.is_empty() {
                0.0
            } else {
                false_r / rejected.len() as f64
            };
            power_sum += true_r / 10.0;
        }
        let fdr = fdp_sum / f64::from(SIMS);
        let power = power_sum / f64::from(SIMS);
        assert!(fdr <= alpha + 0.03, "FDR {fdr} exceeds {alpha}");
        assert!(
            power > 0.8,
            "power {power} too low — thresholds miscomputed?"
        );
        println!(
            "{{\"suite\":\"fs-eproc\",\"case\":\"e-bh\",\"verdict\":\"pass\",\"detail\":\"FDR {fdr:.3} <= {alpha}, power {power:.2}\"}}"
        );
    }

    #[test]
    fn racing_decides_and_is_bit_replayable() {
        // A genuinely better candidate must win; the full e-trajectory must
        // replay bit-identically from the same stream key (Bet 8).
        let run = || -> (bool, u64, u64) {
            let mut s = key(777).stream();
            let mut race = PairwiseRace::new(LossSpan::ONE);
            let mut t = 0u64;
            while !race.a_beats_b(0.05) && t < 50_000 {
                let a = 0.4 + 0.1 * s.next_f64(); // better (lower loss)
                let b = 0.5 + 0.1 * s.next_f64();
                race.observe(a, b).expect("difference lies in [-1, 1]");
                t += 1;
            }
            (race.a_beats_b(0.05), t, race.log_e_value().to_bits())
        };
        let (won1, t1, bits1) = run();
        let (won2, t2, bits2) = run();
        assert!(won1, "the better candidate must win the race");
        assert_eq!(
            (won1, t1, bits1),
            (won2, t2, bits2),
            "tournament must be bit-replayable"
        );
        println!(
            "{{\"suite\":\"fs-eproc\",\"case\":\"race-replay\",\"verdict\":\"pass\",\"detail\":\"decided at t={t1}, bitwise replayable\"}}"
        );
    }

    #[test]
    fn pairwise_scale_is_checked_before_wealth_changes() {
        let mut race = PairwiseRace::new(LossSpan::ONE);
        let initial = race.log_e_value().to_bits();
        let outside = f64::from_bits(1.0f64.to_bits() + 1);
        assert!(matches!(
            race.observe(0.0, outside),
            Err(PairwiseInputError::DifferenceOutOfRange { .. })
        ));
        assert_eq!(race.log_e_value().to_bits(), initial);
        race.observe(0.0, 1.0).expect("inclusive boundary");
        assert!(LossSpan::new(f64::NAN).is_err());
        assert!(LossSpan::new(0.0).is_err());
    }

    #[test]
    fn skew_equal_mean_losses_cannot_be_clipped_into_evidence() {
        // Candidate B is 4 with probability 3/4 and 0 with probability
        // 1/4; candidate A is always 3. Both means are exactly 3. The
        // old silent clamp changed B-A from {+1,-3} to {+1,-1}, whose
        // positive mean manufactured evidence that A beats equal-mean B.
        let mut race = PairwiseRace::new(LossSpan::ONE);
        for loss_b in [4.0, 4.0, 4.0] {
            race.observe(3.0, loss_b).expect("upper boundary");
        }
        let before = race.log_e_value().to_bits();
        assert!(matches!(
            race.observe(3.0, 0.0),
            Err(PairwiseInputError::DifferenceOutOfRange { .. })
        ));
        assert_eq!(race.log_e_value().to_bits(), before);
    }

    #[test]
    fn e_value_arithmetic_laws() {
        // Product of independent e-values is an e-value; average is valid
        // under dependence; e-to-p is Markov.
        let logs = [0.5f64, 1.2, -0.3];
        assert!((combine_product(&logs) - 1.4).abs() < 1e-12);
        let avg = combine_average(&logs);
        let direct = det::ln((det::exp(0.5) + det::exp(1.2) + det::exp(-0.3)) / 3.0);
        assert!((avg - direct).abs() < 1e-12, "{avg} vs {direct}");
        assert!((e_to_p(det::ln(20.0)) - 0.05).abs() < 1e-12);
        assert!(e_to_p(-5.0) <= 1.0);
        // Empty and single-element edges.
        assert_eq!((combine_product(&[]) + 0.0).to_bits(), 0.0f64.to_bits()); // +0 normalizes -0
        assert_eq!((combine_average(&[]) + 0.0).to_bits(), 0.0f64.to_bits());
        assert!(e_benjamini_hochberg(&[], 0.1).is_empty());
    }

    #[test]
    fn bounded_input_contract_is_enforced() {
        let mut e = BettingEProcess::new(0.5);
        let result = std::panic::catch_unwind(core::panic::AssertUnwindSafe(|| {
            let _ = e.observe(1.5);
        }));
        assert!(
            result.is_err(),
            "unbounded input must be refused (it voids validity)"
        );
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}
