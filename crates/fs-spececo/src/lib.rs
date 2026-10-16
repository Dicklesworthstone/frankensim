//! fs-spececo — certified-speculation accept/reject ECONOMICS (plan addendum,
//! Proposal 9). Layer: L6 (an orchestration/decision layer; no numerical
//! dependencies).
//!
//! Untrusted proposers (fs-verify's proposer zoo) generate candidate
//! solutions; a certified verifier bounds each candidate's error. This crate
//! is the DECISION-AND-LEARNING layer over that:
//!
//! - **The decision** ([`decide`]): accept a candidate OUTRIGHT iff its
//!   certified bound meets the query tolerance; otherwise warm-start the true
//!   solver from the candidate and measure iterations saved. Correctness never
//!   depends on the proposer — a garbage candidate (or a non-finite bound) can
//!   only trigger a warm-start, never a false accept.
//! - **The telemetry** ([`ProposerTelemetry`]): per-`(proposer, regime)`
//!   accept rates and warm-start savings — the [`SolveRecord`] fields
//!   (`proposer_id, accepted, bound, iterations_saved`) are exactly the
//!   Error-Ledger solve-node schema this feeds. This telemetry is
//!   simultaneously the surrogate training signal, the query planner's cost
//!   model, and the drift detector's input.
//! - **Drift demotion** ([`DriftDetector`]): an accept-rate collapse in a
//!   regime demotes the proposer THERE — but only after a minimum sample
//!   count, and with HYSTERESIS, so a single unlucky reject cannot
//!   demote-then-re-promote (no flapping).
//!
//! Everything here is deterministic (no RNG, no I/O): replaying the same
//! records reproduces the same telemetry, decisions, and demotions.

use std::collections::BTreeMap;

/// What to do with a proposer's candidate given its certified error bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Decision {
    /// The certified bound meets the query tolerance — accept outright, no
    /// true solve needed.
    AcceptOutright,
    /// The bound does not meet tolerance (or is not a finite, non-negative
    /// bound) — warm-start the true solver from the candidate.
    WarmStart,
}

/// Decide a candidate's fate. Accept outright iff `certified_bound` is a
/// finite, non-negative number that is `<= query_tolerance`; otherwise
/// warm-start. This fails safe: a non-finite or negative bound never accepts.
#[must_use]
pub fn decide(certified_bound: f64, query_tolerance: f64) -> Decision {
    if certified_bound.is_finite()
        && certified_bound >= 0.0
        && query_tolerance.is_finite()
        && certified_bound <= query_tolerance
    {
        Decision::AcceptOutright
    } else {
        Decision::WarmStart
    }
}

/// One recorded speculation outcome — the fs-ledger solve-node telemetry
/// fields (`proposer_id`, `accepted`, `bound`, `iterations_saved`).
#[derive(Debug, Clone, PartialEq)]
pub struct SolveRecord {
    /// Which proposer produced the candidate.
    pub proposer_id: String,
    /// The problem regime (e.g. dimensionless-group bucket).
    pub regime: String,
    /// Was the candidate accepted OUTRIGHT (no true solve)?
    pub accepted: bool,
    /// The certified error bound achieved.
    pub bound: f64,
    /// Iterations saved by warm-starting from the candidate. May be NEGATIVE
    /// (a warm-start that was WORSE than a cold start).
    pub iterations_saved: i64,
}

impl SolveRecord {
    /// A record.
    #[must_use]
    pub fn new(
        proposer_id: impl Into<String>,
        regime: impl Into<String>,
        accepted: bool,
        bound: f64,
        iterations_saved: i64,
    ) -> SolveRecord {
        SolveRecord {
            proposer_id: proposer_id.into(),
            regime: regime.into(),
            accepted,
            bound,
            iterations_saved,
        }
    }
}

/// Accumulated statistics for one `(proposer, regime)` pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Stats {
    /// Total candidates recorded.
    pub attempts: u64,
    /// Candidates accepted outright.
    pub accepts: u64,
    /// Warm-starts that actually saved iterations (`iterations_saved > 0`).
    pub positive_saves: u64,
    /// Net iterations saved (sum; may be negative if warm-starts hurt).
    pub net_iterations_saved: i64,
}

/// Per-`(proposer, regime)` telemetry accumulator.
#[derive(Debug, Clone, Default)]
pub struct ProposerTelemetry {
    stats: BTreeMap<(String, String), Stats>,
}

impl ProposerTelemetry {
    /// An empty accumulator.
    #[must_use]
    pub fn new() -> ProposerTelemetry {
        ProposerTelemetry {
            stats: BTreeMap::new(),
        }
    }

    /// Fold a record into the telemetry.
    pub fn record(&mut self, rec: &SolveRecord) {
        let s = self
            .stats
            .entry((rec.proposer_id.clone(), rec.regime.clone()))
            .or_default();
        s.attempts += 1;
        if rec.accepted {
            s.accepts += 1;
        }
        // A warm-start counts as a "win" only if it genuinely saved iterations;
        // a negative saving is never a win, but it does lower the net total.
        if rec.iterations_saved > 0 {
            s.positive_saves += 1;
        }
        s.net_iterations_saved = s.net_iterations_saved.saturating_add(rec.iterations_saved);
    }

    /// Stats for a pair, if any have been recorded.
    #[must_use]
    pub fn stats(&self, proposer: &str, regime: &str) -> Option<&Stats> {
        self.stats.get(&(proposer.to_string(), regime.to_string()))
    }

    /// The accept-outright rate for a pair. Zero (conservative) when there is
    /// no telemetry — never a divide-by-zero.
    #[must_use]
    pub fn accept_rate(&self, proposer: &str, regime: &str) -> f64 {
        match self.stats(proposer, regime) {
            Some(s) if s.attempts > 0 => s.accepts as f64 / s.attempts as f64,
            _ => 0.0,
        }
    }

    /// The mean iterations saved per candidate for a pair (0 when no
    /// telemetry). May be negative if warm-starts hurt on average.
    #[must_use]
    pub fn mean_iterations_saved(&self, proposer: &str, regime: &str) -> f64 {
        match self.stats(proposer, regime) {
            Some(s) if s.attempts > 0 => s.net_iterations_saved as f64 / s.attempts as f64,
            _ => 0.0,
        }
    }
}

/// A drift detector: demotes a proposer in a regime when its accept-rate
/// collapses, and restores it when the rate recovers — with HYSTERESIS
/// (`demote_below < restore_above`) and a minimum sample count, so a single
/// unlucky reject cannot demote-then-re-promote (no flapping).
#[derive(Debug, Clone)]
pub struct DriftDetector {
    min_samples: u64,
    demote_below: f64,
    restore_above: f64,
    demoted: BTreeMap<(String, String), bool>,
}

impl DriftDetector {
    /// A detector. `demote_below` is the accept-rate at/under which a proposer
    /// is demoted; `restore_above` (should be `> demote_below`) is the rate it
    /// must exceed to be restored; `min_samples` gates both against noise.
    ///
    /// # Panics
    /// If `restore_above <= demote_below` (the hysteresis band would be
    /// inverted OR zero-width — an equal-threshold band has no dead-zone, so a
    /// rate hovering at the shared threshold flaps demote↔restore on
    /// consecutive updates, defeating the whole point of hysteresis).
    #[must_use]
    pub fn new(min_samples: u64, demote_below: f64, restore_above: f64) -> DriftDetector {
        assert!(
            restore_above > demote_below,
            "hysteresis band inverted or degenerate: restore_above must be > demote_below (a \
             zero-width band flaps)"
        );
        DriftDetector {
            min_samples,
            demote_below,
            restore_above,
            demoted: BTreeMap::new(),
        }
    }

    /// Update the demotion state for a pair from current telemetry and return
    /// whether it is now demoted. Below `min_samples` the state is unchanged
    /// (noise-gated).
    pub fn update(&mut self, telem: &ProposerTelemetry, proposer: &str, regime: &str) -> bool {
        let key = (proposer.to_string(), regime.to_string());
        let attempts = telem.stats(proposer, regime).map_or(0, |s| s.attempts);
        if attempts < self.min_samples {
            // not enough evidence to change state.
            return *self.demoted.get(&key).unwrap_or(&false);
        }
        let rate = telem.accept_rate(proposer, regime);
        let currently = *self.demoted.get(&key).unwrap_or(&false);
        // `rate` is a finite accept-rate in [0, 1], so `<=` is total here.
        let next = if currently {
            // stay demoted until the rate clears the upper threshold.
            rate <= self.restore_above
        } else {
            // demote only when the rate drops to/under the lower threshold.
            rate <= self.demote_below
        };
        self.demoted.insert(key, next);
        next
    }

    /// Is the pair currently demoted? (Pure query; does not update state.)
    #[must_use]
    pub fn is_demoted(&self, proposer: &str, regime: &str) -> bool {
        *self
            .demoted
            .get(&(proposer.to_string(), regime.to_string()))
            .unwrap_or(&false)
    }
}
