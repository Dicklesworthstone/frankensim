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

use core::fmt;

pub use fs_eproc::LossSpan;
use fs_eproc::{PairwiseInputError, PairwiseRace, combine_average, e_benjamini_hochberg};
use fs_exec::KillRegistry;

/// Racing controls.
#[derive(Debug, Clone, Copy)]
pub struct RaceSettings {
    /// False-discovery-rate elimination level α (e-BH across the population).
    pub alpha: f64,
    /// Round budget (the fixed-N design would spend this per
    /// candidate).
    pub max_rounds: u32,
    /// Rounds before the first elimination check (e-processes need a
    /// few observations before crossings mean anything; checks before
    /// this are skipped, never peeked).
    pub min_rounds: u32,
    /// Finite positive bound on `abs(loss_a - loss_b)` for every pair
    /// in every round. This is part of the validity and replay identity,
    /// not an observed-data tuning parameter.
    pub loss_span: LossSpan,
}

impl RaceSettings {
    /// Standard alpha and round budgets with an explicit statistical scale.
    #[must_use]
    pub const fn new(loss_span: LossSpan) -> Self {
        RaceSettings {
            alpha: 0.05,
            max_rounds: 400,
            min_rounds: 8,
            loss_span,
        }
    }
}

/// A tournament that cannot issue a statistically valid outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RaceError {
    /// A race needs at least two candidates.
    TooFewCandidates {
        /// Number of supplied candidates.
        count: usize,
    },
    /// Alpha must be finite and strictly between zero and one.
    InvalidAlpha {
        /// IEEE-754 bits of the rejected alpha.
        alpha_bits: u64,
    },
    /// The round budget must be nonzero and include `min_rounds`.
    InvalidRoundBudget {
        /// Requested first round at which elimination may occur.
        min_rounds: u32,
        /// Requested total round budget.
        max_rounds: u32,
    },
    /// One pair violated the declared support, so no race evidence is valid.
    PairwiseInput {
        /// One-based round where the invalid pair appeared.
        round: u32,
        /// First candidate in the oriented pair.
        candidate_a: usize,
        /// Second candidate in the oriented pair.
        candidate_b: usize,
        /// Exact paired-input contract violation.
        source: PairwiseInputError,
    },
    /// A non-finite loss invalidated the optional-stopping construction.
    NonFiniteLoss {
        /// One-based round where the invalid value appeared.
        round: u32,
        /// Candidate that emitted the invalid value.
        candidate: usize,
        /// IEEE-754 bits of the invalid value.
        value_bits: u64,
    },
}

impl fmt::Display for RaceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RaceError::TooFewCandidates { count } => {
                write!(f, "a race needs at least two candidates; got {count}")
            }
            RaceError::InvalidAlpha { alpha_bits } => write!(
                f,
                "race alpha must be finite and in (0, 1); got {}",
                f64::from_bits(*alpha_bits)
            ),
            RaceError::InvalidRoundBudget {
                min_rounds,
                max_rounds,
            } => write!(
                f,
                "race round budget must satisfy 1 <= min_rounds <= max_rounds; got {min_rounds} and {max_rounds}"
            ),
            RaceError::PairwiseInput {
                round,
                candidate_a,
                candidate_b,
                source,
            } => write!(
                f,
                "race lost its validity claim at round {round}, pair ({candidate_a}, {candidate_b}): {source}"
            ),
            RaceError::NonFiniteLoss {
                round,
                candidate,
                value_bits,
            } => write!(
                f,
                "race lost its validity claim when candidate {candidate} emitted non-finite loss {} at round {round}",
                f64::from_bits(*value_bits)
            ),
        }
    }
}

impl std::error::Error for RaceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RaceError::PairwiseInput { source, .. } => Some(source),
            _ => None,
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
    /// Declared paired-loss support used by every e-process.
    pub loss_span: LossSpan,
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

fn validate_race_settings(n_candidates: usize, settings: RaceSettings) -> Result<(), RaceError> {
    if n_candidates < 2 {
        return Err(RaceError::TooFewCandidates {
            count: n_candidates,
        });
    }
    if !settings.alpha.is_finite() || settings.alpha <= 0.0 || settings.alpha >= 1.0 {
        return Err(RaceError::InvalidAlpha {
            alpha_bits: settings.alpha.to_bits(),
        });
    }
    if settings.min_rounds == 0
        || settings.max_rounds == 0
        || settings.min_rounds > settings.max_rounds
    {
        return Err(RaceError::InvalidRoundBudget {
            min_rounds: settings.min_rounds,
            max_rounds: settings.max_rounds,
        });
    }
    Ok(())
}

fn update_mean(mean: &mut f64, count: &mut u64, value: f64) {
    if *count == 0 {
        *mean = value;
        *count = 1;
        return;
    }
    let next = *count + 1;
    let n = next as f64;
    if mean.is_sign_positive() == value.is_sign_positive() {
        *mean += (value - *mean) / n;
    } else {
        *mean = (*mean).mul_add(1.0 - 1.0 / n, value / n);
    }
    *count = next;
}

/// Race a field of candidates with e-BH false-discovery-rate control.
/// `loss(candidate, observation)` must be a PURE function of its
/// arguments (deterministic streams — the caller keys them by seed and
/// candidate id). Every paired difference must lie inside
/// [`RaceSettings::loss_span`]; violations return [`RaceError`] before
/// an invalid outcome can escape. Eliminated candidates' gates fire
/// through `kills` (register ids `0..n` before racing to hold handles).
///
/// VALIDITY (bead 7tv.7.1 — derivation in CONTRACT.md): candidate i's
/// elimination evidence is the MIXTURE (average) of its pairwise
/// e-processes over the FIXED, predeclared opponent family — all n−1
/// original opponents, with a dead opponent's process frozen at its
/// elimination round (a stopped supermartingale is a supermartingale).
/// Under i's null ("i is not worse than any opponent") every term has
/// expectation ≤ 1 at every stopping time, so the mixture is itself an
/// anytime-valid e-process; e-BH then controls the elimination FDR
/// under arbitrary dependence (Wang–Ramdas). The former construction —
/// the MAXIMUM over currently-surviving opponents — is not an e-value
/// and inflates false elimination; the battery's certifier test
/// demonstrates the inflation and pins this one below it.
///
/// Non-finite losses abort the race with no verdict. Condemning only the
/// offending stream would select a freeze time using the next observation,
/// which is not a valid stopped-supermartingale construction.
///
/// # Errors
/// Invalid settings, paired-loss support violations, or a field with no
/// valid candidate. These errors carry no race verdict.
pub fn race_field(
    loss: &mut dyn FnMut(usize, u64) -> f64,
    n_candidates: usize,
    settings: RaceSettings,
    kills: &KillRegistry,
) -> Result<RaceOutcome, RaceError> {
    validate_race_settings(n_candidates, settings)?;
    let prototype = PairwiseRace::new(settings.loss_span);
    let n = n_candidates;
    // Pairwise race matrix (i, j), i < j: PairwiseRace observing
    // (loss_i, loss_j); a_beats_b == "i dominates j".
    let mut races = vec![prototype; n * n];
    let mut alive: Vec<bool> = vec![true; n];
    let mut means = vec![0.0f64; n];
    let mut counts = vec![0u64; n];
    let mut eliminated: Vec<(u32, usize)> = Vec::new();
    let mut evaluations_used = 0u64;
    let mut round = 0u32;
    while round < settings.max_rounds && alive.iter().filter(|&&a| a).count() > 1 {
        // One observation per survivor, canonical order. A non-finite
        // loss aborts the entire evidence path: selectively freezing a
        // component at tau-1 after observing tau is not optional stopping.
        let mut obs: Vec<Option<f64>> = vec![None; n];
        for i in 0..n {
            if !alive[i] {
                continue;
            }
            evaluations_used += 1;
            let v = loss(i, u64::from(round));
            if v.is_finite() {
                update_mean(&mut means[i], &mut counts[i], v);
                obs[i] = Some(v);
            } else {
                return Err(RaceError::NonFiniteLoss {
                    round: round + 1,
                    candidate: i,
                    value_bits: v.to_bits(),
                });
            }
        }
        // Feed every live pair in BOTH directions: slot (i, j) tracks
        // the evidence that i beats j, slot (j, i) the reverse.
        for i in 0..n {
            for j in (i + 1)..n {
                if let (Some(a), Some(b)) = (obs[i], obs[j]) {
                    races[i * n + j]
                        .observe(a, b)
                        .map_err(|source| RaceError::PairwiseInput {
                            round: round + 1,
                            candidate_a: i,
                            candidate_b: j,
                            source,
                        })?;
                    races[j * n + i]
                        .observe(b, a)
                        .map_err(|source| RaceError::PairwiseInput {
                            round: round + 1,
                            candidate_a: j,
                            candidate_b: i,
                            source,
                        })?;
                }
            }
        }
        round += 1;
        if round < settings.min_rounds {
            continue;
        }
        // Elimination evidence per survivor i: the mixture over the
        // FIXED family of all n−1 ORIGINAL opponents — slot (j, i)
        // tracks "j beats i"; a dead j's process is frozen, which is
        // a stopped supermartingale, still valid in the mixture.
        let live: Vec<usize> = (0..n).filter(|&i| alive[i]).collect();
        let mut log_e: Vec<f64> = Vec::with_capacity(live.len());
        let mut family: Vec<f64> = Vec::with_capacity(n - 1);
        for &i in &live {
            family.clear();
            for j in 0..n {
                if j != i {
                    family.push(races[j * n + i].log_e_value());
                }
            }
            log_e.push(combine_average(&family));
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
    // The overflow-safe online means stay finite even when long streams
    // contain values near the ends of the finite f64 range.
    let winner = survivors
        .iter()
        .copied()
        .min_by(|&a, &b| means[a].total_cmp(&means[b]).then(a.cmp(&b)))
        .expect("finite race inputs leave at least one survivor");
    Ok(RaceOutcome {
        survivors,
        eliminated,
        winner,
        evaluations_used,
        fixed_n_equivalent: n as u64 * u64::from(settings.max_rounds),
        rounds: round,
        loss_span: settings.loss_span,
    })
}

/// Successive-halving bracket: at each budget milestone, the bottom
/// (1 − 1/eta) of survivors BY RUNNING MEAN are killed (rank-based —
/// the standard SH semantics, which does NOT carry the e-guarantee;
/// documented, ledgered per bracket).
#[derive(Debug, Clone)]
pub struct BracketLedger {
    /// (milestone round, survivors before, survivors after).
    pub brackets: Vec<(u32, usize, usize)>,
    /// Candidates structurally rejected for non-finite losses
    /// (`(round, candidate)` — rank-based bracket semantics only).
    pub invalid: Vec<(u32, usize)>,
    /// The outcome fields shared with [`RaceOutcome`].
    pub winner: usize,
    /// Loss evaluations consumed.
    pub evaluations_used: u64,
    /// Fixed-N equivalent.
    pub fixed_n_equivalent: u64,
}

/// Run a successive-halving tournament with reduction factor `eta`.
/// Non-finite losses condemn their candidate structurally (fail
/// closed), exactly as in [`race_field`].
///
/// # Panics
/// If `n_candidates < 2`, `eta < 2`, or every candidate is
/// structurally invalid.
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
    let mut means = vec![0.0f64; n];
    let mut counts = vec![0u64; n];
    let mut invalid: Vec<(u32, usize)> = Vec::new();
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
                    let v = loss(i, u64::from(round));
                    if v.is_finite() {
                        update_mean(&mut means[i], &mut counts[i], v);
                    } else {
                        alive[i] = false;
                        invalid.push((round + 1, i));
                        let _ = kills.kill(i as u64);
                    }
                }
            }
            round += 1;
        }
        let mut live: Vec<usize> = (0..n).filter(|&i| alive[i]).collect();
        let before = live.len();
        // Overflow-safe online means make this a finite total order.
        live.sort_by(|&a, &b| means[a].total_cmp(&means[b]).then(a.cmp(&b)));
        let keep = (before as u32).div_ceil(eta).max(1) as usize;
        for &i in &live[keep.min(live.len())..] {
            alive[i] = false;
            let _ = kills.kill(i as u64);
        }
        brackets.push((round, before, keep.min(live.len())));
        total_budget = total_budget.max(u64::from(round));
        milestone *= eta;
    }
    let winner = (0..n)
        .filter(|&i| alive[i])
        .min_by(|&a, &b| means[a].total_cmp(&means[b]).then(a.cmp(&b)))
        .expect("no valid candidate survived: every loss stream produced non-finite values");
    BracketLedger {
        brackets,
        invalid,
        winner,
        evaluations_used,
        fixed_n_equivalent: n as u64 * u64::from(round),
    }
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
