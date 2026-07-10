//! fs-bisect — git-bisect for a wrong number (plan addendum, Proposal 10),
//! plus the Gauntlet failure-compounding workflow ([`compound`], bead 6nb.9).
//! Layer: L6 (a version-control/orchestration concern; no numerical deps).
//!
//! Given a linear COMMIT SEQUENCE and a monotone predicate ("QoI within spec",
//! "certificate passes"), binary-search for the first commit that broke a
//! design — `O(log n)` predicate evaluations instead of `O(n)`. Two economies
//! layer on top:
//!
//! - **Two-tier fidelity** ([`bisect_two_tier`]): run the inner search on a
//!   cheap LOW-fidelity oracle (whose localization is only *estimated*), then
//!   CONFIRM the culprit at FULL fidelity (*verified*). If the full-fidelity
//!   oracle rejects the low-fidelity candidate, the search is re-run at full
//!   fidelity rather than returning a confident wrong culprit.
//! - **Non-monotonicity is detected, not assumed** ([`bisect_checked`],
//!   [`verify_monotone`]): a predicate that goes Good→Bad→Good (a regression
//!   introduced, fixed, then reintroduced) violates the bisect precondition;
//!   rather than return garbage, it is flagged with a witness.
//!
//! A bisect is only sound if replaying commit `k` reproduces commit `k`'s
//! state — that determinism is the ledger's `at(t)`/ExecMode contract, not
//! this crate's; here the [`CommitOracle`] IS that replay-plus-predicate.
//! Everything here is deterministic and side-effect-free; the oracle's own
//! (possibly expensive) evaluation runs under the caller's cancellation scope.

pub mod compound;

/// A commit's verdict under the predicate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// The design is good at this commit (predicate holds).
    Good,
    /// The design is broken at this commit (regressed).
    Bad,
}

/// Evaluate the predicate at a commit index. Commit 0 is the oldest (an
/// assumed-good baseline); `len-1` is the newest.
pub trait CommitOracle {
    /// The verdict at `commit`.
    fn evaluate(&self, commit: usize) -> Verdict;
}

/// The outcome of a bisect.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BisectResult {
    /// The first BAD commit — the culprit. `confirmed` is true when it was
    /// re-checked at full fidelity (a *verified* localization) rather than
    /// only at the search fidelity (an *estimated* one).
    Culprit {
        /// The culprit commit index.
        index: usize,
        /// Was it confirmed at full fidelity?
        confirmed: bool,
    },
    /// Every commit is Good — no regression in the range.
    AllGood,
    /// Every commit is Bad — already broken at the base.
    AllBad,
    /// The predicate is NON-MONOTONE — a Bad commit is followed by a later
    /// Good one, so bisect's precondition is violated. The witness names both.
    NonMonotone {
        /// The index of a Bad commit...
        bad: usize,
        /// ...followed by this later Good commit.
        later_good: usize,
    },
    /// The range is empty.
    Empty,
}

/// A bisect run: its result plus the search path (for structured logging).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BisectRun {
    /// The outcome.
    pub result: BisectResult,
    /// The `(commit, verdict)` probes made, in order.
    pub probes: Vec<(usize, Verdict)>,
}

/// Bisect a commit sequence, ASSUMING the predicate is monotone (Good…Good,
/// then Bad…Bad). The returned culprit is `confirmed = false` (single
/// fidelity). Boundary cases: `len == 0` → `Empty`; base already Bad → `AllBad`;
/// newest still Good → `AllGood`.
#[must_use]
pub fn bisect<O: CommitOracle + ?Sized>(len: usize, oracle: &O) -> BisectRun {
    let mut probes = Vec::new();
    let mut probe = |i: usize| {
        let v = oracle.evaluate(i);
        probes.push((i, v));
        v
    };
    if len == 0 {
        return BisectRun {
            result: BisectResult::Empty,
            probes,
        };
    }
    if probe(0) == Verdict::Bad {
        return BisectRun {
            result: BisectResult::AllBad,
            probes,
        };
    }
    if probe(len - 1) == Verdict::Good {
        return BisectRun {
            result: BisectResult::AllGood,
            probes,
        };
    }
    // invariant: probe(lo) == Good, probe(hi) == Bad, lo < hi.
    let mut lo = 0usize;
    let mut hi = len - 1;
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if probe(mid) == Verdict::Bad {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    BisectRun {
        result: BisectResult::Culprit {
            index: hi,
            confirmed: false,
        },
        probes,
    }
}

/// Scan for a non-monotonicity witness: the first Bad commit that is followed
/// by a later Good one. `None` if the sequence is monotone. `O(n)`.
#[must_use]
pub fn verify_monotone<O: CommitOracle + ?Sized>(len: usize, oracle: &O) -> Option<(usize, usize)> {
    let mut first_bad: Option<usize> = None;
    for i in 0..len {
        match oracle.evaluate(i) {
            Verdict::Bad if first_bad.is_none() => first_bad = Some(i),
            Verdict::Good => {
                if let Some(b) = first_bad {
                    return Some((b, i));
                }
            }
            Verdict::Bad => {}
        }
    }
    None
}

/// Bisect, but first VERIFY monotonicity (`O(n)`): a non-monotone predicate is
/// reported as [`BisectResult::NonMonotone`] rather than mis-localized.
#[must_use]
pub fn bisect_checked<O: CommitOracle + ?Sized>(len: usize, oracle: &O) -> BisectRun {
    if let Some((bad, later_good)) = verify_monotone(len, oracle) {
        let probes = (0..len).map(|i| (i, oracle.evaluate(i))).collect();
        return BisectRun {
            result: BisectResult::NonMonotone { bad, later_good },
            probes,
        };
    }
    bisect(len, oracle)
}

/// Two-tier bisect: narrow with a cheap `low`-fidelity oracle, then CONFIRM the
/// culprit at `full` fidelity. If `full` rejects the low-fidelity candidate
/// (the localization was wrong), re-search entirely at full fidelity. The
/// returned culprit is always `confirmed = true`.
#[must_use]
pub fn bisect_two_tier<L, F>(len: usize, low: &L, full: &F) -> BisectRun
where
    L: CommitOracle + ?Sized,
    F: CommitOracle + ?Sized,
{
    let low_run = bisect(len, low);
    if let BisectResult::Culprit { index, .. } = low_run.result {
        // confirm the candidate locally at full fidelity: it must be Bad, and
        // its predecessor (if any) must be Good.
        let mut probes = low_run.probes;
        let at = full.evaluate(index);
        probes.push((index, at));
        let pred_good = if index == 0 {
            true
        } else {
            let p = full.evaluate(index - 1);
            probes.push((index - 1, p));
            p == Verdict::Good
        };
        if at == Verdict::Bad && pred_good {
            return BisectRun {
                result: BisectResult::Culprit {
                    index,
                    confirmed: true,
                },
                probes,
            };
        }
        // full fidelity disagrees — re-search authoritatively.
        let full_run = bisect(len, full);
        probes.extend(full_run.probes);
        return BisectRun {
            result: confirm(full_run.result),
            probes,
        };
    }
    // low fidelity found no culprit (all-good / all-bad / empty) — confirm the
    // verdict at full fidelity.
    let full_run = bisect(len, full);
    BisectRun {
        result: confirm(full_run.result),
        probes: full_run.probes,
    }
}

/// Mark a bisect result as full-fidelity-confirmed.
fn confirm(r: BisectResult) -> BisectResult {
    match r {
        BisectResult::Culprit { index, .. } => BisectResult::Culprit {
            index,
            confirmed: true,
        },
        other => other,
    }
}
