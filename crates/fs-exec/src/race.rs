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

/// A boxed strategy body (named for clippy's sake; the shape is the
/// contract: consume a per-branch `Cx`, return the candidate or observe
/// the kill).
type BranchFn<'a, T> = Box<dyn FnOnce(&Cx<'_>) -> Result<T, Cancelled> + Send + 'a>;

/// One competing strategy. The closure runs under its own [`Cx`] (own
/// cancel gate, own scope arena, own stream key) and must poll
/// `cx.checkpoint()` at bounded strides, returning `Err(Cancelled)`
/// promptly when killed.
pub struct RaceBranch<'a, T> {
    name: &'static str,
    run: BranchFn<'a, T>,
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
        if self.reports.is_empty() {
            return write!(
                f,
                "race refused: no branches were supplied (an empty race can never \
                 produce a winner; the former behavior hung the parent watcher — wf9.8.1)"
            );
        }
        write!(
            f,
            "race produced no winner: every branch was rejected, cancelled, or panicked — \
             {}; loosen the victory predicate, raise branch budgets, or add a fallback branch",
            race_json(self.mode, None, &self.reports)
        )
    }
}

impl core::error::Error for NoWinner {}

/// Panic-total lock: a poisoned mutex degrades to its data, never to a
/// second panic that could re-hang the drain protocol (wf9.8.1).
fn relock<T>(m: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    m.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
}

fn race_json(mode: &str, winner: Option<usize>, reports: &[BranchReport]) -> String {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(96);
    s.push_str("{\"mode\":");
    push_json_string(&mut s, mode);
    let _ = write!(
        s,
        ",\"winner\":{},\"branches\":[",
        winner.map_or_else(|| "null".to_string(), |w| w.to_string())
    );
    for (i, r) in reports.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        let _ = write!(s, "{{\"index\":{},\"name\":", r.index);
        push_json_string(&mut s, r.name);
        let _ = write!(s, ",\"outcome\":\"{}\"}}", r.outcome.name());
    }
    s.push_str("]}");
    s
}

fn push_json_string(out: &mut String, value: &str) {
    use std::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            other => out.push(other),
        }
    }
    out.push('"');
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
    // One coherent protocol (watcher -> branches -> outcome-keyed decision);
    // splitting would scatter the drain/decision invariants exec-010 audits.
    #[allow(clippy::too_many_lines)]
    pub fn race_with_gate<T: Send>(
        &self,
        branches: Vec<RaceBranch<'_, T>>,
        accept: impl Fn(&T) -> bool + Sync,
        parent: &CancelGate,
    ) -> Result<RaceRun<T>, NoWinner> {
        let n = branches.len();
        let mode = self.config.mode;
        // Structured refusal, never a hung watcher (wf9.8.1): with zero
        // branches nothing would ever set the completion flag.
        if n == 0 {
            return Err(NoWinner {
                reports: Vec::new(),
                mode: mode.name(),
            });
        }
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
                let names = &names;
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
                    // The WHOLE fallible path — branch body AND the
                    // acceptance predicate — runs inside one unwind
                    // guard (wf9.8.1): an accept panic used to unwind
                    // the worker before its slot was filled, hanging
                    // the parent watcher forever.
                    let outcome = arenas.scope(|arena| {
                        let cx =
                            Cx::new(gate, arena, key, asupersync::types::Budget::INFINITE, mode);
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            run(&cx).map(|value| {
                                let accepted = accept(&value);
                                (value, accepted)
                            })
                        }))
                    });
                    let entry = match outcome {
                        Ok(Ok((value, true))) => {
                            // Accepted: contend for leadership and kill the
                            // branches that can no longer win.
                            let mut d = relock(decision);
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
                        Ok(Ok((_, false))) => Err(BranchOutcome::Rejected),
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
                    // Terminal-slot guarantee: this epilogue is panic-
                    // free (poison-tolerant locks, atomic store), so the
                    // watcher ALWAYS gets released.
                    *relock(&slots[index]) = Some(entry);
                    if slots.iter().all(|s| relock(s).is_some()) {
                        race_done_ref.store(true, std::sync::atomic::Ordering::Release);
                    }
                });
            }
        });

        // Decide from OUTCOMES, not timing: the winner is the lowest-index
        // accepted result (Deterministic) or the recorded first arrival
        // (Fast). Scope join above guarantees every loser fully drained.
        let recorded_leader = decision
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .leader;
        let mut winner: Option<(usize, T)> = None;
        let mut reports: Vec<BranchReport> = Vec::with_capacity(n);
        for (index, slot) in slots.into_iter().enumerate() {
            let entry = slot
                .into_inner()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .expect("terminal-slot guarantee: the guarded epilogue always records");
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
    fn race_json_escapes_hostile_branch_names() {
        let reports = [BranchReport {
            name: "quote\" slash\\ line\n tab\t control\u{0001}",
            index: 0,
            outcome: BranchOutcome::Won,
        }];
        assert_eq!(
            race_json("deterministic", Some(0), &reports),
            "{\"mode\":\"deterministic\",\"winner\":0,\"branches\":[{\"index\":0,\"name\":\"quote\\\" slash\\\\ line\\n tab\\t control\\u0001\",\"outcome\":\"won\"}]}"
        );
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

    /// Run a race on a helper thread with a hard time bound: the
    /// wf9.8.1 failure mode is a HANG, so every regression here must
    /// stay bounded even when it fails.
    fn bounded<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(f());
        });
        rx.recv_timeout(std::time::Duration::from_secs(20))
            .expect("race returned within the bound (wf9.8.1: no hangs)")
    }

    #[test]
    fn empty_race_is_refused_not_hung() {
        let err = bounded(|| {
            let racer = Racer::new(RacerConfig::new(0xACE));
            racer
                .race(Vec::<RaceBranch<'_, u64>>::new(), |_| true)
                .map(|r| r.winner)
                .expect_err("empty race refuses")
        });
        assert!(err.reports.is_empty());
        assert!(err.to_string().contains("no branches"), "{err}");
    }

    #[test]
    fn accept_panic_is_contained_and_the_race_terminates() {
        let (winner, reports) = bounded(|| {
            let racer = Racer::new(RacerConfig::new(0xACE));
            let run = racer
                .race(
                    vec![
                        RaceBranch::new("poison", |_cx| Ok(7u64)),
                        RaceBranch::new("clean", |_cx| Ok(8u64)),
                    ],
                    |v| {
                        assert!(*v != 7, "acceptance predicate bomb");
                        true
                    },
                )
                .expect("clean branch wins");
            assert!(racer.arena_pool().stats().quiescent());
            (run.winner, run.reports)
        });
        assert_eq!(winner, 1);
        assert!(matches!(
            &reports[0].outcome,
            BranchOutcome::Panicked { message } if message.contains("predicate bomb")
        ));
        assert_eq!(reports[1].outcome, BranchOutcome::Won);
    }

    #[test]
    fn simultaneous_instant_completions_pick_the_lowest_index() {
        for _ in 0..32 {
            let winner = bounded(|| {
                let racer = Racer::new(RacerConfig::new(0xACE));
                racer
                    .race(
                        vec![
                            RaceBranch::new("a", |_cx| Ok(1u64)),
                            RaceBranch::new("b", |_cx| Ok(2u64)),
                            RaceBranch::new("c", |_cx| Ok(3u64)),
                        ],
                        |_| true,
                    )
                    .expect("winner")
                    .winner
            });
            assert_eq!(winner, 0, "deterministic under simultaneous completion");
        }
    }

    #[test]
    fn mid_flight_parent_cancellation_drains_everything() {
        let reports = bounded(|| {
            let racer = Racer::new(RacerConfig::new(0xACE));
            let parent = std::sync::Arc::new(CancelGate::new());
            let killer = std::sync::Arc::clone(&parent);
            let killer_thread = std::thread::spawn(move || {
                std::thread::sleep(std::time::Duration::from_millis(20));
                killer.request();
            });
            let err = racer
                .race_with_gate(
                    vec![
                        RaceBranch::new("a", spin_until_cancelled),
                        RaceBranch::new("b", spin_until_cancelled),
                    ],
                    |_| true,
                    &parent,
                )
                .expect_err("cancelled tree has no winner");
            killer_thread.join().expect("killer joins");
            assert!(racer.arena_pool().stats().quiescent(), "losers drained");
            err.reports
        });
        assert!(
            reports
                .iter()
                .all(|r| r.outcome == BranchOutcome::Cancelled)
        );
    }

    /// G4 STORM (wf9.8.1 acceptance): races run under registry-owned
    /// candidate gates while eliminations storm in from outside; every
    /// kill lands on a REGISTERED handle (nonzero, structured), every
    /// race returns, and the arenas end quiescent.
    #[test]
    fn g4_kill_storm_hits_registered_gates_and_drains() {
        let out = bounded(|| {
            let registry = std::sync::Arc::new(crate::kill::KillRegistry::new());
            let racer = Racer::new(RacerConfig::new(0xACE));
            let mut landed = 0u32;
            for candidate in 0..12u64 {
                let gate = registry.register(candidate);
                let storm_reg = std::sync::Arc::clone(&registry);
                let storm = std::thread::spawn(move || {
                    std::thread::sleep(std::time::Duration::from_millis(5));
                    storm_reg
                        .kill_registered(candidate)
                        .expect("storm kills only registered candidates");
                });
                let err = racer
                    .race_with_gate(
                        vec![
                            RaceBranch::new("x", spin_until_cancelled),
                            RaceBranch::new("y", spin_until_cancelled),
                        ],
                        |_| true,
                        &gate,
                    )
                    .expect_err("stormed candidate has no winner");
                assert!(
                    err.reports
                        .iter()
                        .all(|r| r.outcome == BranchOutcome::Cancelled)
                );
                storm.join().expect("storm joins");
                landed += 1;
                assert!(registry.release(candidate));
            }
            assert!(racer.arena_pool().stats().quiescent(), "arenas quiescent");
            (landed, registry.live())
        });
        assert_eq!(out.0, 12, "nonzero registered kills, all landed");
        assert_eq!(out.1, 0, "registry drained");
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
