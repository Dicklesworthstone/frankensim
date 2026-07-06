//! Speculative races with loser-cancellation (plan §5.2 behavior 1):
//! "try Delaunay refinement AND octree CutFEM concurrently; the first to
//! meet the error budget wins; the loser's scope is cancelled and its
//! arenas reclaimed."
//!
//! Deterministic victory rule (P2): the winner is the LOWEST-INDEX branch
//! whose accepted result exists — a pure function of branch outcomes, never
//! of arrival order. Early kills remain aggressive: branch `j` is cancelled
//! the moment any branch `i < j` produces an accepted result (j can no
//! longer win), and once a decision is sealed every other branch drains.
//! Arrival timing moves WHEN kills land, never WHO wins.
//! `ExecMode::Fast` relaxes the rule to first-arrival-wins, and the mode
//! rides every report (provenance).
//!
//! Deterministic liveness caveat (the honest price of P2): a branch BELOW
//! the current leader cannot be killed until it terminates on its own —
//! its outcome could still displace the leader. Branches must therefore be
//! budget/poll-disciplined (bounded work per poll, own error budgets);
//! Fast mode has no such requirement (everyone dies on first acceptance).

use crate::cx::{CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use core::fmt;
use std::sync::Mutex;

/// Race configuration (arena config feeds per-branch scope arenas).
#[derive(Debug, Clone)]
pub struct RacerConfig {
    /// Study seed for branch stream keys.
    pub seed: u64,
    /// Victory rule: `Deterministic` = lowest accepted index; `Fast` =
    /// first accepted arrival (recorded).
    pub mode: ExecMode,
    /// Arena configuration for branch scopes.
    pub arena: fs_alloc::ArenaConfig,
}

impl RacerConfig {
    /// Deterministic-mode config with default arenas.
    #[must_use]
    pub fn new(seed: u64) -> Self {
        RacerConfig {
            seed,
            mode: ExecMode::Deterministic,
            arena: fs_alloc::ArenaConfig::default(),
        }
    }
}

/// One competing strategy. The closure runs under its own [`Cx`] (own
/// cancel gate, own scope arena, own stream key) and must poll
/// `cx.checkpoint()` at bounded strides, returning `Err(Cancelled)`
/// promptly when killed.
pub struct RaceBranch<'a, T> {
    name: &'static str,
    run: Box<dyn FnOnce(&Cx<'_>) -> Result<T, Cancelled> + Send + 'a>,
}

impl<'a, T> RaceBranch<'a, T> {
    /// Name a strategy.
    pub fn new(
        name: &'static str,
        run: impl FnOnce(&Cx<'_>) -> Result<T, Cancelled> + Send + 'a,
    ) -> Self {
        RaceBranch {
            name,
            run: Box::new(run),
        }
    }
}

/// What happened to one branch (the "per-branch progress at kill time" of
/// the bead's log requirement).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BranchOutcome {
    /// Produced the accepted result that won.
    Won,
    /// Produced an accepted result, but a lower-index (Deterministic) or
    /// earlier (Fast) branch won.
    Outraced,
    /// Completed, but the victory predicate rejected its result.
    Rejected,
    /// Observed its kill and drained cleanly.
    Cancelled,
    /// Panicked; contained with the message (siblings unaffected).
    Panicked {
        /// The panic payload's message.
        message: String,
    },
}

impl BranchOutcome {
    fn name(&self) -> &'static str {
        match self {
            BranchOutcome::Won => "won",
            BranchOutcome::Outraced => "outraced",
            BranchOutcome::Rejected => "rejected",
            BranchOutcome::Cancelled => "cancelled",
            BranchOutcome::Panicked { .. } => "panicked",
        }
    }
}

/// One branch's report row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BranchReport {
    /// Branch name.
    pub name: &'static str,
    /// Branch index (the deterministic tie-break key).
    pub index: usize,
    /// Outcome.
    pub outcome: BranchOutcome,
}

/// A finished race: the winner and every branch's fate.
#[derive(Debug)]
pub struct RaceRun<T> {
    /// Winning branch index.
    pub winner: usize,
    /// The winning value.
    pub value: T,
    /// Per-branch outcomes (index order).
    pub reports: Vec<BranchReport>,
    /// The victory rule that produced this result (provenance).
    pub mode: &'static str,
}

impl<T> RaceRun<T> {
    /// Canonical JSON summary (deterministic order) for race-outcome logs.
    #[must_use]
    pub fn to_json(&self) -> String {
        race_json(self.mode, Some(self.winner), &self.reports)
    }
}

/// Structured race failure: nobody produced an accepted result.
#[derive(Debug)]
pub struct NoWinner {
    /// Per-branch outcomes (index order) — the diagnosis.
    pub reports: Vec<BranchReport>,
    /// The victory rule in force.
    pub mode: &'static str,
}

impl fmt::Display for NoWinner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "race produced no winner: every branch was rejected, cancelled, or panicked — \
             {}; loosen the victory predicate, raise branch budgets, or add a fallback branch",
            race_json(self.mode, None, &self.reports)
        )
    }
}

impl core::error::Error for NoWinner {}

fn race_json(mode: &str, winner: Option<usize>, reports: &[BranchReport]) -> String {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(96);
    let _ = write!(
        s,
        "{{\"mode\":\"{mode}\",\"winner\":{},\"branches\":[",
        winner.map_or_else(|| "null".to_string(), |w| w.to_string())
    );
    for (i, r) in reports.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(
            s,
            "{{\"index\":{},\"name\":\"{}\",\"outcome\":\"{}\"}}",
            r.index,
            r.name,
            r.outcome.name()
        );
    }
    s.push_str("]}");
    s
}

/// Speculative-race executor. Owns the arena pool whose quiescence after a
/// race proves the losers were fully drained and reclaimed.
pub struct Racer {
    config: RacerConfig,
    arenas: fs_alloc::ArenaPool,
}

/// Shared decision state (drives early kills; the FINAL winner is
/// recomputed from completed outcomes so timing cannot leak in).
struct Decision {
    /// Lowest accepted index so far (Deterministic) or first accepted
    /// arrival (Fast).
    leader: Option<usize>,
}

impl Racer {
    /// Build a racer.
    #[must_use]
    pub fn new(config: RacerConfig) -> Self {
        let arenas = fs_alloc::ArenaPool::new(config.arena.clone());
        Racer { config, arenas }
    }

    /// The arena pool backing branch scopes (G4 leak oracle).
    #[must_use]
    pub fn arena_pool(&self) -> &fs_alloc::ArenaPool {
        &self.arenas
    }

    /// Race branches to the first ACCEPTED result under the configured
    /// victory rule; losers are killed and drained before this returns.
    ///
    /// # Errors
    /// [`NoWinner`] when every branch is rejected, cancelled, or panicked.
    pub fn race<T: Send>(
        &self,
        branches: Vec<RaceBranch<'_, T>>,
        accept: impl Fn(&T) -> bool + Sync,
    ) -> Result<RaceRun<T>, NoWinner> {
        self.race_with_gate(branches, accept, &CancelGate::new())
    }

    /// [`Racer::race`] under an external parent gate (a kill-handle from
    /// the statistical-preemption registry cancels the WHOLE race tree).
    ///
    /// # Errors
    /// [`NoWinner`] as in [`Racer::race`] — including when the parent gate
    /// cancels everything.
    pub fn race_with_gate<T: Send>(
        &self,
        branches: Vec<RaceBranch<'_, T>>,
        accept: impl Fn(&T) -> bool + Sync,
        parent: &CancelGate,
    ) -> Result<RaceRun<T>, NoWinner> {
        let n = branches.len();
        let mode = self.config.mode;
        let gates: Vec<CancelGate> = (0..n).map(|_| CancelGate::new()).collect();
        let slots: Vec<Mutex<Option<Result<T, BranchOutcome>>>> =
            (0..n).map(|_| Mutex::new(None)).collect();
        let names: Vec<&'static str> = branches.iter().map(|b| b.name).collect();
        let decision = Mutex::new(Decision { leader: None });
        let race_done = std::sync::atomic::AtomicBool::new(false);

        std::thread::scope(|s| {
            // Parent-tree kill propagation: a registry kill-handle must
            // cancel the WHOLE race mid-flight, so a watcher folds the
            // parent gate into every branch gate at a bounded stride.
            let watcher_gates = &gates;
            let watcher_done = &race_done;
            s.spawn(move || {
                while !watcher_done.load(std::sync::atomic::Ordering::Acquire) {
                    if parent.is_requested() {
                        for g in watcher_gates {
                            g.request();
                        }
                        return;
                    }
                    std::thread::sleep(std::time::Duration::from_micros(50));
                }
            });
            for (index, branch) in branches.into_iter().enumerate() {
                let gates = &gates;
                let slots = &slots;
                let decision = &decision;
                let accept = &accept;
                let arenas = &self.arenas;
                let seed = self.config.seed;
                let race_done_ref = &race_done;
                s.spawn(move || {
                    let gate = &gates[index];
                    // Parent-tree cancellation folds into the branch gate.
                    if parent.is_requested() {
                        gate.request();
                    }
                    let key = StreamKey {
                        seed,
                        kernel_id: fs_obs::fnv1a64(names[index].as_bytes()),
                        tile: index as u64,
                        iteration: 0,
                    };
                    let run = branch.run;
                    let outcome = arenas.scope(|arena| {
                        let cx =
                            Cx::new(gate, arena, key, asupersync::types::Budget::INFINITE, mode);
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(&cx)))
                    });
                    let entry = match outcome {
                        Ok(Ok(value)) if accept(&value) => {
                            // Accepted: contend for leadership and kill the
                            // branches that can no longer win.
                            let mut d = decision.lock().expect("race decision");
                            let leads = match (mode, d.leader) {
                                (_, None) => true,
                                (ExecMode::Deterministic, Some(l)) => index < l,
                                (ExecMode::Fast, Some(_)) => false,
                            };
                            if leads {
                                d.leader = Some(index);
                            }
                            let leader = d.leader.unwrap_or(index);
                            drop(d);
                            for (j, g) in gates.iter().enumerate() {
                                let doomed = match mode {
                                    // j > leader can never win the
                                    // lowest-index rule; j < leader may
                                    // still displace it.
                                    ExecMode::Deterministic => j > leader,
                                    ExecMode::Fast => j != leader,
                                };
                                if doomed {
                                    g.request();
                                }
                            }
                            Ok(value)
                        }
                        Ok(Ok(_)) => Err(BranchOutcome::Rejected),
                        Ok(Err(Cancelled)) => Err(BranchOutcome::Cancelled),
                        Err(payload) => {
                            let message = payload
                                .downcast_ref::<&str>()
                                .map(ToString::to_string)
                                .or_else(|| payload.downcast_ref::<String>().cloned())
                                .unwrap_or_else(|| "non-string panic payload".to_string());
                            Err(BranchOutcome::Panicked { message })
                        }
                    };
                    *slots[index].lock().expect("race slot") = Some(entry);
                    // Last branch out releases the parent watcher.
                    if slots.iter().all(|s| s.lock().expect("race slot").is_some()) {
                        race_done_ref.store(true, std::sync::atomic::Ordering::Release);
                    }
                });
            }
        });

        // Decide from OUTCOMES, not timing: the winner is the lowest-index
        // accepted result (Deterministic) or the recorded first arrival
        // (Fast). Scope join above guarantees every loser fully drained.
        let recorded_leader = decision.into_inner().expect("race decision").leader;
        let mut winner: Option<(usize, T)> = None;
        let mut reports: Vec<BranchReport> = Vec::with_capacity(n);
        for (index, slot) in slots.into_iter().enumerate() {
            let entry = slot
                .into_inner()
                .expect("race slot")
                .expect("every branch reports");
            match entry {
                Ok(value) => {
                    let take = match mode {
                        ExecMode::Deterministic => winner.is_none(),
                        ExecMode::Fast => recorded_leader == Some(index),
                    };
                    if take {
                        winner = Some((index, value));
                        reports.push(BranchReport {
                            name: names[index],
                            index,
                            outcome: BranchOutcome::Won,
                        });
                    } else {
                        reports.push(BranchReport {
                            name: names[index],
                            index,
                            outcome: BranchOutcome::Outraced,
                        });
                    }
                }
                Err(outcome) => reports.push(BranchReport {
                    name: names[index],
                    index,
                    outcome,
                }),
            }
        }
        match winner {
            Some((index, value)) => Ok(RaceRun {
                winner: index,
                value,
                reports,
                mode: mode.name(),
            }),
            None => Err(NoWinner {
                reports,
                mode: mode.name(),
            }),
        }
    }
}

impl fmt::Debug for Racer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Racer")
            .field("mode", &self.config.mode.name())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spin_until_cancelled(cx: &Cx<'_>) -> Result<u64, Cancelled> {
        loop {
            cx.checkpoint()?;
            std::hint::spin_loop();
        }
    }

    #[test]
    fn lowest_accepted_index_wins_and_losers_drain() {
        let racer = Racer::new(RacerConfig::new(0xACE));
        let run = racer
            .race(
                vec![
                    RaceBranch::new("slow-but-lowest", |cx| {
                        // Real work; finishes later than branch 2 usually.
                        let mut acc = 1u64;
                        for i in 0..200_000u64 {
                            if i % 4096 == 0 {
                                cx.checkpoint()?;
                            }
                            acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i);
                        }
                        Ok(acc | 1)
                    }),
                    RaceBranch::new("spins-forever", spin_until_cancelled),
                    RaceBranch::new("instant", |_cx| Ok(7)),
                ],
                |v| *v > 0,
            )
            .expect("race has a winner");
        assert_eq!(run.winner, 0, "deterministic rule: lowest accepted index");
        assert_eq!(run.reports[1].outcome, BranchOutcome::Cancelled);
        assert_eq!(run.reports[2].outcome, BranchOutcome::Outraced);
        assert!(racer.arena_pool().stats().quiescent(), "losers drained");
        assert!(run.to_json().contains("\"winner\":0"), "{}", run.to_json());
    }

    #[test]
    fn rejected_and_panicked_branches_yield_a_teaching_no_winner() {
        let racer = Racer::new(RacerConfig::new(0xACE));
        let err = racer
            .race(
                vec![
                    RaceBranch::new("rejected", |_cx| Ok(-1i64)),
                    RaceBranch::new("bomb", |_cx| -> Result<i64, Cancelled> {
                        panic!("strategy exploded");
                    }),
                ],
                |v| *v > 0,
            )
            .expect_err("no winner");
        assert_eq!(err.reports[0].outcome, BranchOutcome::Rejected);
        assert!(matches!(
            &err.reports[1].outcome,
            BranchOutcome::Panicked { message } if message.contains("exploded")
        ));
        assert!(err.to_string().contains("fallback"), "{err}");
        assert!(racer.arena_pool().stats().quiescent());
    }

    #[test]
    fn parent_gate_kills_the_whole_race_tree() {
        let racer = Racer::new(RacerConfig::new(0xACE));
        let parent = CancelGate::new();
        parent.request();
        let err = racer
            .race_with_gate(
                vec![
                    RaceBranch::new("a", spin_until_cancelled),
                    RaceBranch::new("b", spin_until_cancelled),
                ],
                |_| true,
                &parent,
            )
            .expect_err("pre-cancelled tree has no winner");
        assert!(
            err.reports
                .iter()
                .all(|r| r.outcome == BranchOutcome::Cancelled),
            "{err}"
        );
        assert!(racer.arena_pool().stats().quiescent());
    }
}
