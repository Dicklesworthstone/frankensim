//! FLYWHEEL CLOSES (bead lmp4.18; [F] — behind the `flywheel-e2e`
//! feature until its Gauntlet tier is green): the whole-loop harness
//! testing the addendum's CENTRAL CLAIM — that speculation (9),
//! incremental recompute (2), the sheaf-adjudicated merge (10), and
//! tombstones (E) COMPOUND, not merely work in isolation.
//!
//! The workload is the benchmark corpus's design-iteration model on the
//! CHT wedge: two concurrent agents iterate a design whose per-edit DAG
//! sizes and certifiable skip sets come from the recorded edit traces
//! (`fs-benchmark`), with a seeded fraction of candidate designs being
//! π-equivalent re-visits of tombstoned dead ends. COSTS ARE MODELED
//! UNITS from the corpus (real API calls, modeled physics): the loop
//! MECHANICS are measured; wall-clock physics lands with the wedge's
//! real solvers.
//!
//! The measurement (review round 3): isolated speedups per proposal and
//! the composed loop, over N seeded replays, asserting composed >
//! max(isolated) by a stated margin with across-replay variance
//! reported — plus laundering-across-the-loop (an estimated speculation
//! result is never upgraded anywhere downstream), whole-loop
//! determinism (G5: identical trace hashes), and a cancellation storm
//! (G4: a mid-loop cancel leaves a consistent ledger and no residue).
#![cfg(feature = "flywheel-e2e")]

use std::collections::BTreeMap;

use fs_evidence::{Color, IntervalOp, compose};
use fs_geom::sheaf_merge::{BranchState, MergeOutcome, three_way_merge};
use fs_geom::sheaf_repair::SheafSkeleton;
use fs_ledger::hash_bytes;
use fs_ledger::tombstone::{Descriptor, ExplorationVerdict, TombstoneIndex};
use fs_qty::{Dims, QtyAny};
use fs_recompute::{NodeRecord, ParamValue, SkipDecision, Store};
use fs_spececo::{Decision, ProposerTelemetry, SolveRecord, decide};

/// Which proposals are switched ON for a run (a feature-toggle matrix
/// by design: the harness measures every on/off combination).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct LoopConfig {
    /// Proposal 9: certified speculation (warm starts).
    pub speculation: bool,
    /// Proposal 2: incremental recompute (skips).
    pub recompute: bool,
    /// Proposal 10: sheaf-adjudicated merge (parallel credit).
    pub merge: bool,
    /// Proposal E: tombstone gate (dead candidates blocked).
    pub tombstones: bool,
    /// Cancel after this many stage transitions (G4 storm), if any.
    pub cancel_after_stages: Option<usize>,
}

impl LoopConfig {
    /// Everything off: the baseline.
    #[must_use]
    pub fn baseline() -> LoopConfig {
        LoopConfig {
            speculation: false,
            recompute: false,
            merge: false,
            tombstones: false,
            cancel_after_stages: None,
        }
    }

    /// Everything on: the closed loop.
    #[must_use]
    pub fn composed() -> LoopConfig {
        LoopConfig {
            speculation: true,
            recompute: true,
            merge: true,
            tombstones: true,
            cancel_after_stages: None,
        }
    }
}

/// One run's report — the whole flywheel's telemetry in one trace.
#[derive(Debug, Clone, PartialEq)]
pub struct LoopReport {
    /// Total modeled cost (wall-analog units).
    pub total_cost: f64,
    /// Iterations completed (== requested unless cancelled).
    pub iterations: usize,
    /// True when a G4 cancel fired mid-loop.
    pub cancelled: bool,
    /// Structured stage events (the trace; hashable for G5).
    pub events: Vec<String>,
    /// Speculation accept rate (0 when off).
    pub accept_rate: f64,
    /// Ops skipped by recompute.
    pub skips: usize,
    /// Merges resolved / conflicted.
    pub merges: (usize, usize),
    /// Candidates blocked by the tombstone gate.
    pub tombstone_blocks: usize,
    /// The headline color carried end-to-end.
    pub headline: Color,
}

impl LoopReport {
    /// The G5 trace hash: every event row, in order.
    #[must_use]
    pub fn trace_hash(&self) -> String {
        let joined = self.events.join("\n");
        hash_bytes(joined.as_bytes()).to_hex()
    }
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64
}

fn wedge_descriptor(name: &str, velocity: f64, scale: f64) -> Descriptor {
    let mut params = BTreeMap::new();
    params.insert(
        "velocity".to_string(),
        QtyAny::new(velocity, Dims([1, 0, -1, 0, 0])),
    );
    params.insert(
        "length".to_string(),
        QtyAny::new(scale, Dims([1, 0, 0, 0, 0])),
    );
    params.insert(
        "viscosity".to_string(),
        QtyAny::new(1.8e-5, Dims([2, 0, -1, 0, 0])),
    );
    Descriptor {
        name: name.to_string(),
        params,
    }
}

/// A 3-patch merge skeleton (the wedge's chart layout stand-in).
fn merge_skeleton() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// Run the loop: `iterations` design steps by two concurrent agents on
/// the corpus's first edit trace, with the configured proposals live.
/// Deterministic in `seed` (G5).
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn run_loop(config: &LoopConfig, iterations: usize, seed: u64) -> LoopReport {
    let trace = fs_benchmark::edit_traces()[0];
    let ops_per_iter = trace.total_ops;
    let skippable = trace.correct_skips;
    let mut state = seed;
    let mut events = Vec::new();
    let mut total_cost = 0.0f64;
    let mut stages = 0usize;
    let mut cancelled = false;

    let mut tombstones = TombstoneIndex::new();
    // Pre-seed the graveyard: three dead designs at known velocities.
    for v in [10.0, 20.0, 40.0] {
        tombstones.record_falsification_kill(
            wedge_descriptor("cht-wedge bracket", v, 0.1),
            "{\"kind\":\"tombstone\"}",
            vec!["estimated".to_string()],
            50.0,
            "2026-07-08",
            "agent:corpus",
        );
    }
    let mut store = Store::new();
    let mut telemetry = ProposerTelemetry::new();
    let mut accepts = 0usize;
    let mut proposals = 0usize;
    let mut skips = 0usize;
    let mut merges_ok = 0usize;
    let mut merges_conflict = 0usize;
    let mut blocks = 0usize;
    // The color the pipeline carries: speculation results are ESTIMATED
    // and must stay so through cache -> merge -> query.
    let mut headline = Color::Verified { lo: 0.0, hi: 1e-9 };
    let skeleton = merge_skeleton();
    let mut done = 0usize;

    'outer: for iter in 0..iterations {
        // ---- Stage E: the tombstone gate over this iteration's candidate.
        stages += 1;
        if let Some(limit) = config.cancel_after_stages
            && stages > limit
        {
            cancelled = true;
            break 'outer;
        }
        // Every third candidate is a REVISIT of a dead design (same
        // descriptor family, same pi-neighborhood); fresh explorations
        // carry genuinely distinct descriptors and physics (a fin array
        // at a different scale — decades away in pi-space).
        let revisit = iter % 3 == 2;
        let velocity = if revisit {
            20.0 + lcg(&mut state) * 0.4
        } else {
            100.0 + 50.0 * lcg(&mut state)
        };
        let candidate = if revisit {
            wedge_descriptor("cht-wedge bracket", velocity, 0.1)
        } else {
            wedge_descriptor(
                &format!("cht-wedge fin-array rev{}", iter % 5),
                velocity,
                0.5,
            )
        };
        if config.tombstones {
            if let ExplorationVerdict::Blocked { .. } = tombstones.pre_exploration_check(&candidate)
            {
                blocks += 1;
                events.push(format!("iter={iter} stage=tombstone verdict=blocked"));
                continue; // the whole candidate's cost is saved
            }
            events.push(format!(
                "iter={iter} stage=tombstone verdict=clear v={velocity:.3}"
            ));
        } else if revisit {
            // Without the gate the dead candidate is fully re-solved by
            // BOTH agents (and then re-discovered dead).
            total_cost += 2.0 * ops_per_iter as f64;
            events.push(format!(
                "iter={iter} stage=dead-resolve cost={ops_per_iter}"
            ));
            continue;
        }

        // ---- Stages 2+9: each agent solves its branch's DAG.
        let mut branch_costs = [0.0f64; 2];
        for (agent, cost) in branch_costs.iter_mut().enumerate() {
            stages += 1;
            if let Some(limit) = config.cancel_after_stages
                && stages > limit
            {
                cancelled = true;
                break 'outer;
            }
            for op in 0..ops_per_iter {
                let record = NodeRecord {
                    op_id: format!("wedge-op-{op}"),
                    input_hashes: Vec::new(),
                    params: vec![(
                        "iter-group".to_string(),
                        // Skippable ops share params across iterations
                        // (the unchanged part of the design); the rest
                        // change every iteration.
                        ParamValue::f(if op < skippable {
                            0.0
                        } else {
                            #[allow(clippy::cast_precision_loss)]
                            {
                                (iter * 2 + agent) as f64
                            }
                        }),
                    )],
                    code_version_hash: hash_bytes(b"wedge-v1"),
                    rng_seed: 7,
                    achieved_error: 1e-8,
                    required_tolerance: 1e-6,
                };
                if config.recompute
                    && matches!(store.can_skip(&record, 1e-6), SkipDecision::Hit { .. })
                {
                    skips += 1;
                    continue; // certified skip: no cost
                }
                // Speculation: a proposer offers a warm start.
                let mut op_cost = 1.0;
                if config.speculation {
                    proposals += 1;
                    let bound = if lcg(&mut state) < 0.7 { 5e-7 } else { 1e-3 };
                    let decision = decide(bound, 1e-6);
                    let accepted = decision == Decision::AcceptOutright;
                    if accepted {
                        accepts += 1;
                        op_cost = 0.15; // verification-only cost
                        // The accepted result is ESTIMATED: compose it
                        // into the headline (weakest input wins).
                        headline = compose(
                            &headline,
                            &Color::Estimated {
                                estimator: "wedge-proposer-v1".to_string(),
                                dispersion: 0.1,
                            },
                            IntervalOp::Hull,
                        );
                    }
                    telemetry.record(&SolveRecord::new(
                        "wedge-proposer-v1",
                        "re-1e5",
                        accepted,
                        bound,
                        if accepted { 40 } else { -2 },
                    ));
                }
                *cost += op_cost;
                let _ = store.put(record, b"artifact");
            }
            events.push(format!(
                "iter={iter} stage=solve agent={agent} cost={cost:.2}"
            ));
        }

        // ---- Stage 10: merge the two branches.
        stages += 1;
        if let Some(limit) = config.cancel_after_stages
            && stages > limit
        {
            cancelled = true;
            break 'outer;
        }
        if config.merge {
            let base = vec![0.0; 3];
            // Gauge-style concurrent edits (occasionally a circulation
            // taint that genuinely conflicts).
            let taint = lcg(&mut state) < 0.1;
            let x = BranchState {
                provenance: format!("agent-x@{iter}"),
                mismatch: if taint {
                    skeleton.d1t(&[0.05])
                } else {
                    skeleton.d0(&[0.0, 0.01 * lcg(&mut state), 0.0])
                },
                assignments: BTreeMap::new(),
            };
            let y = BranchState {
                provenance: format!("agent-y@{iter}"),
                mismatch: skeleton.d0(&[0.0, 0.0, 0.01 * lcg(&mut state)]),
                assignments: BTreeMap::new(),
            };
            match three_way_merge(&skeleton, &base, &x, &y, None, 1e-6, 1e-6) {
                MergeOutcome::Resolved { .. } | MergeOutcome::Trivial { .. } => {
                    merges_ok += 1;
                    // Parallel credit: wall time is the max branch.
                    total_cost += branch_costs[0].max(branch_costs[1]);
                    events.push(format!("iter={iter} stage=merge verdict=resolved"));
                }
                _ => {
                    merges_conflict += 1;
                    // Conflict: serialize + redo the cheaper branch.
                    total_cost +=
                        branch_costs[0] + branch_costs[1] + branch_costs[0].min(branch_costs[1]);
                    events.push(format!("iter={iter} stage=merge verdict=conflict"));
                }
            }
        } else {
            // No merge machinery: agents serialize.
            total_cost += branch_costs[0] + branch_costs[1];
            events.push(format!("iter={iter} stage=serialize"));
        }
        done += 1;
    }

    #[allow(clippy::cast_precision_loss)]
    let accept_rate = if proposals == 0 {
        0.0
    } else {
        accepts as f64 / proposals as f64
    };
    LoopReport {
        total_cost,
        iterations: done,
        cancelled,
        events,
        accept_rate,
        skips,
        merges: (merges_ok, merges_conflict),
        tombstone_blocks: blocks,
        headline,
    }
}

/// Isolated + composed speedups over one seed (baseline_cost / cost).
#[must_use]
pub fn speedups(iterations: usize, seed: u64) -> (BTreeMap<&'static str, f64>, f64) {
    let base = run_loop(&LoopConfig::baseline(), iterations, seed).total_cost;
    let one = |f: fn(&mut LoopConfig)| {
        let mut c = LoopConfig::baseline();
        f(&mut c);
        base / run_loop(&c, iterations, seed).total_cost
    };
    let mut isolated = BTreeMap::new();
    isolated.insert("speculation", one(|c| c.speculation = true));
    isolated.insert("recompute", one(|c| c.recompute = true));
    isolated.insert("merge", one(|c| c.merge = true));
    isolated.insert("tombstones", one(|c| c.tombstones = true));
    let composed = base / run_loop(&LoopConfig::composed(), iterations, seed).total_cost;
    (isolated, composed)
}
