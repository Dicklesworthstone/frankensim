//! CONFORMAL HARDENING (patch Rev M, bead 7tv.9; [F] — behind the
//! `conformal-hardening` feature): the anytime-valid layer's
//! assumptions made into OPERATIONAL CONTRACTS. The system must not
//! advertise distribution-free protection while quietly violating
//! exchangeability by changing design regimes.
//!
//! Three mechanisms:
//! - MONDRIAN buckets: split-conformal calibration PER regime/fidelity/
//!   family bucket — marginal coverage can hide systematic per-regime
//!   failure; a bucket that lacks calibration mass REFUSES rather than
//!   extrapolating another bucket's quantile.
//! - DRIFT e-tests: the training-vs-candidate two-sample question as a
//!   sequential, anytime-valid betting e-process on the probability
//!   integral transform (uniform under exchangeability) — eating our
//!   own dogfood; detection triggers escalation and validity-domain
//!   shrinkage.
//! - FALSE-COVERAGE-RATE budgets: thousands of simultaneous coverage
//!   claims each carry a miscoverage e-process; e-BH over those
//!   e-values spends an EXPLICIT statistical error budget (P4 extended
//!   to coverage), and admission math reserves the per-claim share up
//!   front.

use crate::{BettingEProcess, e_benjamini_hochberg};
use std::collections::BTreeMap;

/// A per-bucket split-conformal calibrator. Buckets are opaque keys
/// composed by the caller (regime class, fidelity, family, ...).
#[derive(Debug, Clone, Default)]
pub struct MondrianConformal {
    buckets: BTreeMap<String, Vec<f64>>,
    /// Minimum calibration count before a bucket may claim coverage
    /// (below it, `band` refuses — the honest failure).
    pub min_calibration: usize,
}

/// A per-bucket band, or the refusal that keeps the guarantee honest.
#[derive(Debug, Clone, PartialEq)]
pub enum BucketBand {
    /// The bucket's split-conformal half-width at the requested level.
    Calibrated {
        /// Half-width (symmetric residual band).
        half_width: f64,
        /// Calibration sample count backing the claim.
        n: usize,
    },
    /// Not enough calibration mass in this bucket: NO coverage claim.
    /// The teaching field says how many more samples are needed.
    Refused {
        /// Samples present.
        have: usize,
        /// Samples required.
        need: usize,
    },
}

impl MondrianConformal {
    /// New calibrator; `min_calibration` gates per-bucket claims.
    #[must_use]
    pub fn new(min_calibration: usize) -> MondrianConformal {
        MondrianConformal {
            buckets: BTreeMap::new(),
            min_calibration,
        }
    }

    /// Record one calibration residual (absolute) in a bucket.
    pub fn add(&mut self, bucket: &str, abs_residual: f64) {
        self.buckets
            .entry(bucket.to_string())
            .or_default()
            .push(abs_residual);
    }

    /// The bucket's band at miscoverage `alpha`, or an honest refusal.
    /// Quantile index follows the split-conformal `⌈(n+1)(1−α)⌉` rule.
    #[must_use]
    pub fn band(&self, bucket: &str, alpha: f64) -> BucketBand {
        let res = self
            .buckets
            .get(bucket)
            .map_or(&[] as &[f64], Vec::as_slice);
        if res.len() < self.min_calibration {
            return BucketBand::Refused {
                have: res.len(),
                need: self.min_calibration,
            };
        }
        let mut sorted = res.to_vec();
        sorted.sort_by(f64::total_cmp);
        let n = sorted.len();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (((n + 1) as f64) * (1.0 - alpha)).ceil() as usize;
        // Split conformal: half-width = the k-th smallest residual, k =
        // ⌈(n+1)(1−α)⌉. When k > n there are too few calibration points for this
        // α — the (n+1)(1−α) quantile of {residuals} ∪ {+∞} lands on +∞, so the
        // honest band is INFINITE (still ≥1−α coverage, trivially). Capping at
        // sorted[n−1] (the max) UNDER-COVERS at exactly n/(n+1) < 1−α (bead
        // q2tf); `k.min(n) − 1` also underflowed usize for n = 0.
        let half_width = if k > n { f64::INFINITY } else { sorted[k - 1] };
        BucketBand::Calibrated { half_width, n }
    }

    /// The MARGINAL band over all buckets pooled (what a non-Mondrian
    /// calibrator would report — kept for the comparison the bead
    /// demands).
    #[must_use]
    pub fn marginal_band(&self, alpha: f64) -> BucketBand {
        let mut all: Vec<f64> = self.buckets.values().flatten().copied().collect();
        if all.len() < self.min_calibration {
            return BucketBand::Refused {
                have: all.len(),
                need: self.min_calibration,
            };
        }
        all.sort_by(f64::total_cmp);
        let n = all.len();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (((n + 1) as f64) * (1.0 - alpha)).ceil() as usize;
        // Infinite band when k > n — see `band` (bead q2tf).
        let half_width = if k > n { f64::INFINITY } else { all[k - 1] };
        BucketBand::Calibrated { half_width, n }
    }

    /// Bucket names (deterministic order).
    #[must_use]
    pub fn bucket_names(&self) -> Vec<String> {
        self.buckets.keys().cloned().collect()
    }
}

/// The drift verdict and its escalation contract.
#[derive(Debug, Clone, PartialEq)]
pub struct DriftVerdict {
    /// Anytime-valid rejection of exchangeability at the monitor's α.
    pub drifted: bool,
    /// Samples consumed when detection fired (0 if not fired).
    pub samples_at_detection: u64,
    /// The validity-domain SHRINK factor in (0, 1]: 1 = full trust;
    /// decays with the e-value once drift is detected. Surrogate
    /// validity domains multiply by this (Evidence updates flow
    /// through).
    pub validity_scale: f64,
}

/// A sequential two-sample drift monitor: candidates are ranked
/// against the frozen training sample; under exchangeability the
/// normalized rank (PIT) is uniform ON AVERAGE, but conditionally on a
/// FINITE training draw it carries an O(1/√n) bias — so the null is
/// COMPOSITE: the betting pair tests mean > 1/2 + δ (and < 1/2 − δ)
/// with δ = max(1/√n_train, 0.02), the finite-calibration tolerance.
/// Real shifts move the PIT mean by far more than δ; the slack is what
/// keeps 'no drift' from false-firing on the training sample's own
/// sampling noise (found empirically: an unslacked monitor false-fired
/// on a 2000-sample null run).
#[derive(Debug, Clone)]
pub struct DriftMonitor {
    train: Vec<f64>,
    up: BettingEProcess,
    down: BettingEProcess,
    alpha: f64,
    fired_at: u64,
    seen: u64,
}

impl DriftMonitor {
    /// Freeze the training distribution and start monitoring at `alpha`.
    #[must_use]
    pub fn new(mut train: Vec<f64>, alpha: f64) -> DriftMonitor {
        train.sort_by(f64::total_cmp);
        #[allow(clippy::cast_precision_loss)]
        let delta = (1.0 / (train.len() as f64).sqrt()).max(0.02);
        DriftMonitor {
            train,
            up: BettingEProcess::new((0.5 + delta).min(0.9)),
            down: BettingEProcess::new((0.5 + delta).min(0.9)),
            alpha,
            fired_at: 0,
            seen: 0,
        }
    }

    /// Observe one candidate value; returns the current verdict.
    pub fn observe(&mut self, x: f64) -> DriftVerdict {
        self.seen += 1;
        // PIT: rank of x in the frozen training sample.
        let below = self.train.partition_point(|&t| t < x);
        #[allow(clippy::cast_precision_loss)]
        let u = (below as f64 + 0.5) / (self.train.len() as f64 + 1.0);
        let _ = self.up.observe(u);
        let _ = self.down.observe(1.0 - u);
        let log_e = self.up.log_e_value().max(self.down.log_e_value());
        let drifted = self.up.rejects_at(self.alpha) || self.down.rejects_at(self.alpha);
        if drifted && self.fired_at == 0 {
            self.fired_at = self.seen;
        }
        DriftVerdict {
            drifted,
            samples_at_detection: self.fired_at,
            validity_scale: if drifted {
                (1.0 / (1.0 + log_e)).clamp(0.05, 1.0)
            } else {
                1.0
            },
        }
    }

    /// Current max log e-value across the two one-sided monitors.
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.up.log_e_value().max(self.down.log_e_value())
    }
}

/// One monitored coverage claim: a miscoverage e-process betting that
/// the miss rate exceeds the advertised `alpha` (null: miss ≤ alpha).
#[derive(Debug, Clone)]
pub struct CoverageClaim {
    /// Claim identifier (ledger key).
    pub name: String,
    eproc: BettingEProcess,
    hits: u64,
    misses: u64,
}

impl CoverageClaim {
    /// New claim advertised at miscoverage `alpha`.
    #[must_use]
    pub fn new(name: &str, alpha: f64) -> CoverageClaim {
        CoverageClaim {
            name: name.to_string(),
            eproc: BettingEProcess::new(alpha),
            hits: 0,
            misses: 0,
        }
    }

    /// Record one prediction outcome (true = the band covered).
    pub fn observe(&mut self, covered: bool) {
        if covered {
            self.hits += 1;
        } else {
            self.misses += 1;
        }
        let _ = self.eproc.observe(if covered { 0.0 } else { 1.0 });
    }

    /// Log e-value against the advertised rate.
    #[must_use]
    pub fn log_e_value(&self) -> f64 {
        self.eproc.log_e_value()
    }

    /// Empirical miss rate.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn miss_rate(&self) -> f64 {
        self.misses as f64 / (self.hits + self.misses).max(1) as f64
    }
}

/// The FALSE-COVERAGE budget over simultaneous claims: e-BH at the
/// budget level flags the claims whose advertised coverage is broken,
/// with FDR control inherited from the e-BH guarantee.
#[must_use]
pub fn fcr_flag(claims: &[CoverageClaim], budget_alpha: f64) -> Vec<usize> {
    let log_e: Vec<f64> = claims.iter().map(CoverageClaim::log_e_value).collect();
    e_benjamini_hochberg(&log_e, budget_alpha)
}

/// ADMISSION math: with `k` simultaneous claims under a total
/// miscoverage budget, the per-claim advertised alpha that keeps the
/// UNION miscoverage within budget (Bonferroni reservation — the
/// conservative admission bound; e-BH monitoring then recovers power
/// at run time).
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn admission_alpha(total_budget: f64, k: usize) -> f64 {
    total_budget / k.max(1) as f64
}

/// The exchangeability MODEL CARD: which conformal claims assume what,
/// declared rather than implied (ledger-ready via `to_json`).
#[derive(Debug, Clone, PartialEq)]
pub struct ExchangeabilityCard {
    /// The bucket scheme in force (Mondrian keys).
    pub bucketing: String,
    /// The drift monitor's alpha (0 = no monitoring — declared!).
    pub drift_alpha: f64,
    /// The FCR budget over simultaneous claims.
    pub fcr_budget: f64,
    /// Refresh policy note (how optimization-induced shift is handled).
    pub refresh_policy: String,
}

impl ExchangeabilityCard {
    /// Structured declaration for the ledger.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"bucketing\":\"{}\",\"drift_alpha\":{},\"fcr_budget\":{},\
             \"refresh_policy\":\"{}\"}}",
            self.bucketing, self.drift_alpha, self.fcr_budget, self.refresh_policy
        )
    }
}
