//! fs-race — e-RACING (plan §9.6, Bet 8 [M]). Layer: L4 (ASCENT).
//!
//! Anytime-valid sequential tests DRIVE structured cancellation: within
//! a generation, per-candidate loss streams feed a full pairwise
//! fs-eproc race matrix; the moment a candidate's elimination evidence
//! crosses the e-BH threshold its kill-handle fires through fs-exec's
//! [`KillRegistry`], cancelling the candidate's whole evaluation tree.
//!
//! BIT-REPRODUCIBLE BY CONSTRUCTION: rounds are the only clock. Every
//! surviving candidate consumes exactly one observation per round in
//! canonical index order, and e-value crossings are evaluated ONLY at
//! round boundaries — the elimination sequence is a pure function of
//! (seed, logical stream identities), never of wall-clock arrival.
//!
//! The [M] discipline: the 2–5× payoff claim is MEASURED (evaluations
//! used vs the fixed-N budget) on separated and inseparable fields
//! alike, and the battery's calibration study checks that the true
//! best is eliminated no more often than α promises. If the payoff
//! were not to materialize on some field, the ledger would say so —
//! that is the point of carrying `fixed_n_equivalent` in the outcome.

use fs_eproc::{PairwiseRace, e_benjamini_hochberg};
use fs_exec::KillRegistry;

/// Racing controls.
#[derive(Debug, Clone, Copy)]
pub struct RaceSettings {
    /// Family-wise elimination level α (e-BH across the population).
    pub alpha: f64,
    /// Round budget (the fixed-N design would spend this per
    /// candidate).
    pub max_rounds: u32,
    /// Rounds before the first elimination check (e-processes need a
    /// few observations before crossings mean anything; checks before
    /// this are skipped, never peeked).
    pub min_rounds: u32,
}

impl Default for RaceSettings {
    fn default() -> Self {
        RaceSettings {
            alpha: 0.05,
            max_rounds: 400,
            min_rounds: 8,
        }
    }
}

/// The tournament record — the auditable ledger row.
#[derive(Debug, Clone)]
pub struct RaceOutcome {
    /// Surviving candidate indices (ascending).
    pub survivors: Vec<usize>,
    /// Elimination events `(round, candidate)` in occurrence order
    /// (within a round: ascending candidate index — deterministic).
    pub eliminated: Vec<(u32, usize)>,
    /// Winner: the surviving candidate with the lowest running mean
    /// loss (ties break by index).
    pub winner: usize,
    /// Loss evaluations actually consumed.
    pub evaluations_used: u64,
    /// What a fixed-N design (every candidate to the full budget)
    /// would have consumed.
    pub fixed_n_equivalent: u64,
    /// Rounds executed.
    pub rounds: u32,
}

impl RaceOutcome {
    /// Evaluations saved as a ratio (≥ 1; the falsifiable payoff).
    #[must_use]
    pub fn savings(&self) -> f64 {
        #[allow(clippy::cast_precision_loss)] // fixture-scale counters
        {
            self.fixed_n_equivalent as f64 / (self.evaluations_used as f64).max(1.0)
        }
    }
}

/// Race a field of candidates with e-BH family-wise elimination.
/// `loss(candidate, observation)` must be a PURE function of its
/// arguments (deterministic streams — the caller keys them by seed and
/// candidate id); eliminated candidates' gates fire through `kills`
/// (register ids `0..n` before racing to hold handles).
///
/// # Panics
/// If `n_candidates < 2` (a race needs a field).
#[must_use]
pub fn race_field(
    loss: &mut dyn FnMut(usize, u64) -> f64,
    n_candidates: usize,
    settings: RaceSettings,
    kills: &KillRegistry,
) -> RaceOutcome {
    assert!(n_candidates >= 2, "a race needs at least two candidates");
    let n = n_candidates;
    // Pairwise race matrix (i, j), i < j: PairwiseRace observing
    // (loss_i, loss_j); a_beats_b == "i dominates j".
    let mut races: Vec<PairwiseRace> = (0..n * n).map(|_| PairwiseRace::new()).collect();
    let mut alive: Vec<bool> = vec![true; n];
    let mut sums = vec![0.0f64; n];
    let mut counts = vec![0u64; n];
    let mut eliminated: Vec<(u32, usize)> = Vec::new();
    let mut evaluations_used = 0u64;
    let mut round = 0u32;
    while round < settings.max_rounds && alive.iter().filter(|&&a| a).count() > 1 {
        // One observation per survivor, canonical order.
        let obs: Vec<Option<f64>> = (0..n)
            .map(|i| {
                if alive[i] {
                    evaluations_used += 1;
                    let v = loss(i, u64::from(round));
                    sums[i] += v;
                    counts[i] += 1;
                    Some(v)
                } else {
                    None
                }
            })
            .collect();
        // Feed every live pair in BOTH directions: slot (i, j) tracks
        // the evidence that i beats j, slot (j, i) the reverse.
        for i in 0..n {
            for j in (i + 1)..n {
                if let (Some(a), Some(b)) = (obs[i], obs[j]) {
                    races[i * n + j].observe(a, b);
                    races[j * n + i].observe(b, a);
                }
            }
        }
        round += 1;
        if round < settings.min_rounds {
            continue;
        }
        // Elimination evidence per survivor: the strongest surviving
        // opponent's log e-value for "opponent beats me".
        let live: Vec<usize> = (0..n).filter(|&i| alive[i]).collect();
        let mut log_e: Vec<f64> = Vec::with_capacity(live.len());
        for &i in &live {
            let mut best = f64::NEG_INFINITY;
            for &j in &live {
                if j == i {
                    continue;
                }
                // Slot (j, i) tracks "j beats i" — both directions
                // are fed each round.
                best = best.max(races[j * n + i].log_e_value());
            }
            log_e.push(best);
        }
        let condemned = e_benjamini_hochberg(&log_e, settings.alpha);
        if !condemned.is_empty() {
            let ids: Vec<usize> = condemned.iter().map(|&k| live[k]).collect();
            for &i in &ids {
                alive[i] = false;
                eliminated.push((round, i));
                let _ = kills.kill(i as u64);
            }
        }
    }
    let survivors: Vec<usize> = (0..n).filter(|&i| alive[i]).collect();
    let winner = survivors
        .iter()
        .copied()
        .min_by(|&a, &b| {
            let ma = sums[a] / counts[a].max(1) as f64;
            let mb = sums[b] / counts[b].max(1) as f64;
            ma.partial_cmp(&mb).expect("finite losses").then(a.cmp(&b))
        })
        .expect("at least one survivor");
    RaceOutcome {
        survivors,
        eliminated,
        winner,
        evaluations_used,
        fixed_n_equivalent: n as u64 * u64::from(settings.max_rounds),
        rounds: round,
    }
}

/// Successive-halving bracket: at each budget milestone, the bottom
/// (1 − 1/eta) of survivors BY RUNNING MEAN are killed (rank-based —
/// the standard SH semantics, which does NOT carry the e-guarantee;
/// documented, ledgered per bracket).
#[derive(Debug, Clone)]
pub struct BracketLedger {
    /// (milestone round, survivors before, survivors after).
    pub brackets: Vec<(u32, usize, usize)>,
    /// The outcome fields shared with [`RaceOutcome`].
    pub winner: usize,
    /// Loss evaluations consumed.
    pub evaluations_used: u64,
    /// Fixed-N equivalent.
    pub fixed_n_equivalent: u64,
}

/// Run a successive-halving tournament with reduction factor `eta`.
///
/// # Panics
/// If `n_candidates < 2` or `eta < 2`.
#[must_use]
pub fn successive_halving(
    loss: &mut dyn FnMut(usize, u64) -> f64,
    n_candidates: usize,
    base_rounds: u32,
    eta: u32,
    kills: &KillRegistry,
) -> BracketLedger {
    assert!(n_candidates >= 2, "a bracket needs at least two candidates");
    assert!(eta >= 2, "eta must halve at least");
    let n = n_candidates;
    let mut alive: Vec<bool> = vec![true; n];
    let mut sums = vec![0.0f64; n];
    let mut counts = vec![0u64; n];
    let mut evaluations_used = 0u64;
    let mut brackets = Vec::new();
    let mut milestone = base_rounds;
    let mut round = 0u32;
    let mut total_budget = 0u64;
    while alive.iter().filter(|&&a| a).count() > 1 {
        while round < milestone {
            for i in 0..n {
                if alive[i] {
                    evaluations_used += 1;
                    sums[i] += loss(i, u64::from(round));
                    counts[i] += 1;
                }
            }
            round += 1;
        }
        let mut live: Vec<usize> = (0..n).filter(|&i| alive[i]).collect();
        let before = live.len();
        live.sort_by(|&a, &b| {
            let ma = sums[a] / counts[a].max(1) as f64;
            let mb = sums[b] / counts[b].max(1) as f64;
            ma.partial_cmp(&mb).expect("finite losses").then(a.cmp(&b))
        });
        let keep = (before as u32).div_ceil(eta).max(1) as usize;
        for &i in &live[keep..] {
            alive[i] = false;
            let _ = kills.kill(i as u64);
        }
        brackets.push((round, before, keep));
        total_budget = total_budget.max(u64::from(round));
        milestone *= eta;
    }
    let winner = (0..n)
        .filter(|&i| alive[i])
        .min_by(|&a, &b| {
            let ma = sums[a] / counts[a].max(1) as f64;
            let mb = sums[b] / counts[b].max(1) as f64;
            ma.partial_cmp(&mb).expect("finite losses").then(a.cmp(&b))
        })
        .expect("one survivor");
    BracketLedger {
        brackets,
        winner,
        evaluations_used,
        fixed_n_equivalent: n as u64 * u64::from(round),
    }
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
