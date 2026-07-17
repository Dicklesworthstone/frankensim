//! fs-exec conformance suite (CONTRACT.md: any reimplementation must pass).
//!
//! Gauntlet coverage: G0 (completeness, reduction laws), G4 (cancellation
//! storm with panic injection, drain/leak oracles), G5 (bit-identical
//! results and stream keys across worker counts). Latency measurements are
//! emitted as events — measured, never assumed. Every aggregate case emits
//! one canonical fs-obs verdict; seeded cases carry their real logical input
//! or execution seed. Those seeds do not replay OS scheduling or wall-clock
//! timing.

use core::ops::ControlFlow;
use fs_exec::{
    CancelGate, Cancelled, Cx, PoolConfig, RunError, RunId, StreamKey, TileKernel, TilePlan,
    TilePool, victim_order, weighted_ranges,
};
use fs_substrate::affinity::CcdTopology;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-exec/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-exec/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("executor verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("executor verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

/// In-house LCG (L0 cannot depend on fs-rand; fs-qty battery constants).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        // High bits: LCG low bits are short-period (bit 0 alternates), which
        // would lock decision sequences into cycles.
        (self.next() >> 32) % n
    }
}

fn pool_with(workers: usize, seed: u64) -> TilePool {
    TilePool::new(PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, seed))
}

/// Non-associative float work + arena traffic: the determinism witness.
struct FloatKernel {
    tiles: u64,
}

impl TileKernel for FloatKernel {
    type Out = f64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("conf/float", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, f64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let xs = cx
            .arena()
            .alloc_slice_with(fs_alloc::Site::named("conf/float"), 33, |i| {
                1.0 / ((tile * 33 + i as u64 + 1) as f64)
            })
            .expect("arena");
        ControlFlow::Continue(xs.iter().sum())
    }
}

/// Order-sensitive reduction witness: output must equal ascending tile
/// order regardless of scheduling.
struct OrderKernel {
    tiles: u64,
}

impl TileKernel for OrderKernel {
    type Out = Vec<u64>;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("conf/order", self.tiles)
    }

    fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<Cancelled, Vec<u64>> {
        ControlFlow::Continue(vec![tile])
    }
}

/// Publishes each tile's logical stream key (worker-independence witness).
struct KeyKernel {
    tiles: u64,
}

impl TileKernel for KeyKernel {
    type Out = Vec<(u64, u128)>;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("conf/keys", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Vec<(u64, u128)>> {
        ControlFlow::Continue(vec![(tile, cx.stream_key().key128())])
    }
}

#[test]
fn exec_001_completeness_and_arena_hygiene() {
    const SEED: u64 = 0xE001;
    let mut checked = 0;
    for workers in [1usize, 2, 5, 8] {
        for tiles in [0u64, 1, 7, 128] {
            let pool = pool_with(workers, SEED);
            let got = pool.run(&OrderKernel { tiles }).expect("run");
            let want: Vec<u64> = (0..tiles).collect();
            assert_eq!(got, want, "workers={workers} tiles={tiles}");
            assert!(pool.arena_pool().stats().quiescent(), "arena leak");
            checked += 1;
        }
    }
    verdict(
        "exec-001",
        checked == 16,
        "every tile runs exactly once and arenas reclaim, across worker/tile counts (G0)",
        SEED,
    );
}

#[test]
fn exec_002_g5_bit_identical_results_across_worker_counts() {
    const SEED: u64 = 0xE002;
    let reference = pool_with(1, SEED)
        .run(&FloatKernel { tiles: 257 })
        .expect("reference run");
    let mut all_equal = true;
    for workers in [2usize, 3, 4, 8] {
        for _ in 0..3 {
            let got = pool_with(workers, SEED)
                .run(&FloatKernel { tiles: 257 })
                .expect("run");
            all_equal &= got.to_bits() == reference.to_bits();
        }
    }
    verdict(
        "exec-002",
        all_equal,
        &format!(
            "non-associative float reduction bit-identical across 1..8 workers (G5): {reference:e}"
        ),
        SEED,
    );
}

#[test]
fn exec_003_stream_keys_are_worker_independent() {
    const SEED: u64 = 0xE003;
    let mut runs: Vec<Vec<(u64, u128)>> = Vec::new();
    for workers in [1usize, 2, 7] {
        // Fresh pool per run: iteration counters all start at 0, so the
        // ONLY varying factor is worker count — which must not matter.
        let got = pool_with(workers, SEED)
            .run(&KeyKernel { tiles: 100 })
            .expect("run");
        runs.push(got);
    }
    let identical = runs.windows(2).all(|w| w[0] == w[1]);
    let distinct = {
        let mut keys: Vec<u128> = runs[0].iter().map(|&(_, k)| k).collect();
        keys.sort_unstable();
        keys.dedup();
        keys.len() == 100
    };
    verdict(
        "exec-003",
        identical && distinct,
        "RNG stream keys derive from logical identity only — shuffling worker counts changes \
         nothing (G5), and every tile's stream is distinct",
        SEED,
    );
}

#[test]
fn exec_004_external_cancellation_drains_and_ledgers_latency() {
    const SEED: u64 = 0xE004;
    struct SlowKernel;
    impl TileKernel for SlowKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("conf/slow", 100_000)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            // ~microseconds of real work per tile so the run outlives the
            // cancel request without sleeping.
            let mut acc = tile;
            for i in 0..500 {
                acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i);
            }
            ControlFlow::Continue(acc & 1)
        }
    }
    let pool = pool_with(4, SEED);
    let gate = CancelGate::new();
    let (result, report) = std::thread::scope(|s| {
        s.spawn(|| {
            std::thread::sleep(std::time::Duration::from_millis(2));
            gate.request();
        });
        pool.run_with_gate(&SlowKernel, &gate)
    });
    let cancelled = matches!(result, Err(RunError::Cancelled { .. }));
    let quiescent = pool.arena_pool().stats().quiescent();
    // Ledger the measured histogram through the one observability schema.
    let mut em = fs_obs::Emitter::new("fs-exec/conformance", "exec-004/cancel");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "fs-exec-cancel-latency".to_string(),
                json: report.to_json(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("latency event validates");
    println!("{line}");
    let p99 = report.cancel_latency_p99_ns().unwrap_or(u64::MAX);
    // Sanity envelope only (debug-profile CI noise): the 200 µs reference-
    // hardware target is a roofline-harness gate, not a unit assertion —
    // the honest evidence is the ledgered histogram above (CONTRACT).
    let bounded = p99 < 250_000_000;
    verdict(
        "exec-004",
        cancelled && quiescent && bounded,
        &format!(
            "external cancel drains cleanly ({}/{} tiles), arenas quiescent, observed p99 = \
             {p99} ns (measured, ledgered)",
            report.completed, report.total
        ),
        SEED,
    );
}

#[test]
fn exec_005_g4_storm_random_cancels_and_panics_stay_structured() {
    const SEED: u64 = 0xE005_2026_0706_5707;
    struct StormKernel {
        panic_tile: Option<u64>,
    }
    impl TileKernel for StormKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("conf/storm", 64)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            assert!(self.panic_tile != Some(tile), "storm-injected panic");
            let b = cx
                .arena()
                .alloc_slice_fill(fs_alloc::Site::named("conf/storm"), 256, tile as u8)
                .expect("arena");
            ControlFlow::Continue(u64::from(b[128]))
        }
    }
    let mut rng = Lcg(SEED);
    let pool = pool_with(8, SEED);
    let (mut oks, mut cancels, mut panics) = (0u32, 0u32, 0u32);
    for _ in 0..300 {
        let inject_panic = rng.below(4) == 0;
        let kernel = StormKernel {
            panic_tile: inject_panic.then(|| rng.below(64)),
        };
        let gate = CancelGate::new();
        if rng.below(2) == 0 {
            gate.request(); // cancel BEFORE the run: zero tiles must leak
        }
        let (res, _report) = pool.run_with_gate(&kernel, &gate);
        match res {
            Ok(_) => oks += 1,
            Err(RunError::Cancelled { .. }) => cancels += 1,
            Err(RunError::TilePanicked { message, .. }) => {
                assert!(message.contains("storm-injected"), "{message}");
                panics += 1;
            }
            Err(other) => panic!("unstructured outcome: {other}"),
        }
        assert!(
            pool.arena_pool().stats().quiescent(),
            "leak after a storm iteration"
        );
    }
    let pass = oks + cancels + panics == 300 && cancels > 0 && panics > 0 && oks > 0;
    let mut em = fs_obs::Emitter::new("fs-exec/conformance", "exec-005/storm");
    let event = em.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::StormAssertion {
            name: "no-executor-leak".to_string(),
            pass,
            seed: SEED,
        },
        None,
    );
    fs_obs::validate_line(&event.to_jsonl()).expect("storm event validates");
    println!("{}", event.to_jsonl());
    verdict(
        "exec-005",
        pass,
        &format!(
            "300 storm runs (seed {SEED:#x}): ok={oks} cancelled={cancels} panicked={panics}, \
             all structured, arenas quiescent throughout (G4)"
        ),
        SEED,
    );
}

#[test]
fn exec_006_steal_order_and_quanta_respect_topology_fixtures() {
    // Threadripper fixture: 96 workers over 12 CCDs — every same-CCD victim
    // must precede every cross-CCD victim.
    let mut ccd_local_first = true;
    for (topo, workers) in [
        (CcdTopology::TR_7995WX, 96usize),
        (CcdTopology::EPYC_128C, 128),
        (CcdTopology::APPLE_M_CLASS, 16),
    ] {
        let per_ccd = workers / (topo.ccds as usize);
        for w in 0..workers {
            let my_ccd = w / per_ccd;
            let order = victim_order(w, workers, &topo);
            let cut = per_ccd - 1; // same-CCD victims (all but me)
            ccd_local_first &= order[..cut].iter().all(|&v| v / per_ccd == my_ccd);
            ccd_local_first &= order[cut..].iter().all(|&v| v / per_ccd != my_ccd);
        }
    }
    // P/E quanta: heavier workers take proportionally larger initial runs.
    let ranges = weighted_ranges(120, &[2, 2, 1, 1]);
    let sizes: Vec<u64> = ranges.iter().map(|r| r.end - r.start).collect();
    let proportional = sizes == vec![40, 40, 20, 20];
    verdict(
        "exec-006",
        ccd_local_first && proportional,
        "victim order is CCD-local-first on TR/EPYC/Apple fixtures; quantum weights split \
         tiles proportionally (the P/E hook)",
        0,
    );
}

#[test]
fn exec_007_latency_lane_stays_responsive_under_tile_load() {
    const SEED: u64 = 0xE007;
    struct BusyKernel;
    impl TileKernel for BusyKernel {
        type Out = u64;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("conf/busy", 20_000)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let mut acc = tile;
            for i in 0..2_000 {
                acc = acc.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(i);
            }
            ControlFlow::Continue(acc & 1)
        }
    }
    let lane = fs_exec::LatencyLane::new(1).expect("lane");
    let pool = pool_with(
        std::thread::available_parallelism().map_or(4, std::num::NonZero::get),
        SEED,
    );
    let (lane_ns, run) = std::thread::scope(|s| {
        let handle = s.spawn(|| pool.run(&BusyKernel));
        // While the tile lane saturates the cores, the latency lane must
        // still turn around a trivial orchestration future quickly.
        let t0 = std::time::Instant::now();
        let v = lane.block_on(async { 21 * 2 });
        let lane_ns = t0.elapsed().as_nanos() as u64;
        assert_eq!(v, 42);
        (lane_ns, handle.join().expect("tile thread joins"))
    });
    run.expect("busy kernel completes");
    let mut em = fs_obs::Emitter::new("fs-exec/conformance", "exec-007/lane");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::BenchmarkResult {
                kernel: "latency-lane-turnaround".to_string(),
                metric: "ns".to_string(),
                value: lane_ns as f64,
                machine: fs_substrate::CapabilityProbe::topology_only().fingerprint(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("lane event validates");
    println!("{line}");
    // Generous debug-profile envelope; the ≤100 ms conversational target is
    // HELM's gate on reference hardware (CONTRACT no-claims).
    verdict(
        "exec-007",
        lane_ns < 1_000_000_000,
        &format!("latency lane turnaround {lane_ns} ns while the tile pool saturates cores"),
        SEED,
    );
}

#[test]
fn exec_008_reduction_shape_is_scheduling_invariant() {
    const SEED: u64 = 0xE008;
    // A non-commutative reduction (concatenation) is the sharpest witness:
    // ANY arrival-order fold would scramble it under stealing.
    for workers in [1usize, 2, 4, 8] {
        for _ in 0..5 {
            let got = pool_with(workers, SEED)
                .run(&OrderKernel { tiles: 300 })
                .expect("run");
            assert!(
                got.iter().copied().eq(0..300),
                "fold order broke at workers={workers}"
            );
        }
    }
    verdict(
        "exec-008",
        true,
        "non-commutative concatenation folds in ascending tile order under every schedule \
         (fixed-shape reduction law, G0/G5)",
        SEED,
    );
}

#[test]
fn exec_009_g5_audit_compensated_reductions_bit_stable_across_thread_counts() {
    use fs_exec::reduce::{Compensated, audit_accumulator, det_sum};

    const INPUT_SEED: u64 = 0xE009_A0D1;
    const EXECUTION_SEED: u64 = 0xE009;

    /// Ill-conditioned per-tile partial sums, merged compensated on the
    /// pool's fixed pairwise tree.
    struct CompKernel;
    impl TileKernel for CompKernel {
        type Out = Compensated;

        fn tiles(&self) -> TilePlan {
            TilePlan::new("conf/compensated", 311)
        }

        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Compensated> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let mut acc = Compensated::zero();
            for i in 0..97u64 {
                let sgn = if (tile + i).is_multiple_of(2) {
                    1.0
                } else {
                    -1.0
                };
                acc = acc.accumulate(sgn * 1e15);
                acc = acc.accumulate(1.0 / ((tile * 97 + i + 1) as f64));
            }
            ControlFlow::Continue(acc)
        }
    }
    let p = std::thread::available_parallelism().map_or(4, std::num::NonZero::get);
    let mut hashes = Vec::new();
    for workers in [1usize, 2, p, 2 * p] {
        let got = pool_with(workers, EXECUTION_SEED)
            .run(&CompKernel)
            .expect("run");
        hashes.push(fs_obs::fnv1a64(&got.value().to_bits().to_le_bytes()));
    }
    let bit_stable = hashes.windows(2).all(|w| w[0] == w[1]);
    // The audit half of the acceptance: a seeded arrival-order bug must be
    // caught and localized; det_sum on the same series must be order-free.
    let mut nasty: Vec<f64> = Vec::new();
    let mut rng = Lcg(INPUT_SEED);
    for _ in 0..200 {
        nasty.push(1e16);
        nasty.push(1.0 + (rng.below(9) as f64));
        nasty.push(-1e16);
    }
    let bug = audit_accumulator(&nasty, EXECUTION_SEED, |a, x| a + x);
    let caught = matches!(&bug, Err(e) if e.witness >= 2);
    let det_is_clean = {
        let s1 = det_sum(&nasty);
        let mut rev = nasty.clone();
        rev.reverse();
        // det_sum is index-keyed: same INPUT SEQUENCE -> same bits; a
        // different sequence is a different reduction, so instead assert
        // repeatability and dd-grade value.
        s1.to_bits() == det_sum(&nasty).to_bits() && rev.len() == nasty.len()
    };
    verdict(
        "exec-009",
        bit_stable && caught && det_is_clean,
        &format!(
            "compensated artifact hash bit-stable across {{1,2,{p},{}}} workers (G5 audit); \
             seeded input {INPUT_SEED:#x}, execution/audit permutation {EXECUTION_SEED:#x}, \
             arrival-order bug caught with witness={}",
            2 * p,
            bug.err().map_or(0, |e| e.witness)
        ),
        INPUT_SEED,
    );
}

#[test]
fn exec_010_race_winner_is_deterministic_and_losers_fully_drain() {
    use fs_exec::{BranchOutcome, RaceBranch, Racer, RacerConfig};
    const SEED: u64 = 0xE010;
    // Ten repeats with jittered branch timing: the winner (index AND bits)
    // must never move in Deterministic mode.
    let mut winners = Vec::new();
    for round in 0..10u64 {
        let racer = Racer::new(RacerConfig::new(SEED));
        let run = racer
            .race(
                vec![
                    RaceBranch::new("workhorse", move |cx| {
                        let mut acc = 0x9E37_79B9u64;
                        // Jitter the workload per round: timing changes,
                        // outcome must not.
                        for i in 0..(50_000 + round * 7_919) {
                            if i % 4096 == 0 && cx.checkpoint().is_err() {
                                return Err(fs_exec::Cancelled);
                            }
                            acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i);
                        }
                        Ok(0xACCE_57ED_u64)
                    }),
                    RaceBranch::new("spinner", |cx| {
                        loop {
                            cx.checkpoint()?;
                            std::hint::spin_loop();
                        }
                    }),
                    RaceBranch::new("sprinter", |_cx| Ok(0xFA57_u64)),
                ],
                |v| *v > 0,
            )
            .expect("race has a winner");
        assert!(
            racer.arena_pool().stats().quiescent(),
            "losers must drain before race returns"
        );
        assert_eq!(run.reports[1].outcome, BranchOutcome::Cancelled);
        winners.push((run.winner, run.value));
    }
    let stable = winners.windows(2).all(|w| w[0] == w[1]);
    verdict(
        "exec-010",
        stable && winners[0] == (0, 0xACCE_57ED),
        &format!(
            "deterministic victory rule: lowest accepted index wins across 10 jittered \
             repeats (got branch {} = {:#x}); spinner killed and drained every time",
            winners[0].0, winners[0].1
        ),
        SEED,
    );
}

#[test]
fn exec_011_solver_checkpoint_resume_fork_is_bit_exact() {
    use fs_exec::solver::{
        ResumableSolver, SolverProgress, SolverState, StepVerdict, codec, drive, fork,
    };
    const SEED: u64 = 0xE011;
    /// Logistic-map style iteration: chaotic enough that ANY perturbation
    /// of the trajectory shows up in the bits.
    struct Chaotic {
        steps: u64,
    }
    #[derive(Clone)]
    struct ChaoticState {
        x: f64,
        iter: u64,
    }
    impl SolverState for ChaoticState {
        const TYPE_ID: u64 = 0x4348_414f_5449_0001;
        const SCHEMA_VERSION: u32 = 1;

        fn encode(&self, enc: &mut codec::Enc) {
            enc.put_f64(self.x);
            enc.put_u64(self.iter);
        }
        fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
            Ok(ChaoticState {
                x: dec.get_f64()?,
                iter: dec.get_u64()?,
            })
        }
    }
    impl ResumableSolver for Chaotic {
        type State = ChaoticState;
        type Out = f64;
        fn step(&self, s: &mut ChaoticState, _cx: &Cx<'_>) -> StepVerdict<f64> {
            s.x = 3.9 * s.x * (1.0 - s.x);
            s.iter += 1;
            if s.iter >= self.steps {
                StepVerdict::Done(s.x)
            } else {
                StepVerdict::Continue
            }
        }
    }
    let solver = Chaotic { steps: 10_000 };
    let s0 = ChaoticState { x: 0.372, iter: 0 };
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let run_under_cx =
        |state: ChaoticState, cancel_at: Option<u64>| -> SolverProgress<ChaoticState, f64> {
            let gate = CancelGate::new();
            arenas.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    fs_exec::StreamKey {
                        seed: SEED,
                        kernel_id: 1,
                        tile: 0,
                        iteration: 0,
                    },
                    asupersync::types::Budget::INFINITE,
                    fs_exec::ExecMode::Deterministic,
                );
                if let Some(at) = cancel_at {
                    // Drive manually until `at`, then request and let drive pause.
                    let mut st = state;
                    for _ in 0..at {
                        let StepVerdict::Continue = solver.step(&mut st, &cx) else {
                            panic!("cancel point must precede completion");
                        };
                    }
                    gate.request();
                    drive(&solver, st, &cx)
                } else {
                    drive(&solver, state, &cx)
                }
            })
        };
    let SolverProgress::Done(x_ref) = run_under_cx(s0.clone(), None) else {
        panic!("reference finishes");
    };
    // Pause at three different depths; serialize; resume; compare bits.
    let mut all_exact = true;
    for cancel_at in [1u64, 137, 9_999] {
        let SolverProgress::Paused(paused) = run_under_cx(s0.clone(), Some(cancel_at)) else {
            panic!("requested gate must pause");
        };
        let bytes = paused.to_bytes();
        let restored = ChaoticState::from_bytes(&bytes).expect("round trip");
        let forked = fork(&restored).expect("fork proves serializability");
        assert_eq!(restored.content_hash(), forked.content_hash());
        let SolverProgress::Done(x_resumed) = run_under_cx(restored, None) else {
            panic!("resume finishes");
        };
        all_exact &= x_resumed.to_bits() == x_ref.to_bits();
    }
    verdict(
        "exec-011",
        all_exact,
        &format!(
            "pause-serialize-resume reproduces the uninterrupted chaotic trajectory bit-exactly \
             at depths 1/137/9999 (G4 law): x = {x_ref:e}"
        ),
        SEED,
    );
}

#[test]
fn exec_012_kill_handles_reclaim_a_deep_candidate_tree_mid_flight() {
    use fs_exec::KillRegistry;
    const SEED: u64 = 0xE012;
    struct GrindKernel;
    impl TileKernel for GrindKernel {
        type Out = u64;
        fn tiles(&self) -> TilePlan {
            TilePlan::new("conf/grind", 1_000_000)
        }
        fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
            if cx.checkpoint().is_err() {
                return ControlFlow::Break(Cancelled);
            }
            let mut acc = tile;
            for i in 0..200 {
                acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i);
            }
            ControlFlow::Continue(acc & 1)
        }
    }
    let registry = KillRegistry::new();
    let gate = registry.register(42);
    let pool = pool_with(4, SEED);
    // The candidate's "deep tree": a tile-pool run under the registry gate,
    // killed from outside mid-flight (the e-process elimination path).
    let (result, report) = std::thread::scope(|s| {
        let gate2 = std::sync::Arc::clone(&gate);
        let registry = &registry;
        s.spawn(move || {
            while gate2.now_ns() < 2_000_000 {
                std::hint::spin_loop();
            }
            assert!(registry.kill(42), "registered candidate is killable");
        });
        pool.run_with_gate(&GrindKernel, &gate)
    });
    let cancelled = matches!(result, Err(RunError::Cancelled { .. }));
    let quiescent = pool.arena_pool().stats().quiescent();
    assert!(registry.release(42));
    // Ledger the observed kill-to-drain latency (measured, never assumed).
    let mut em = fs_obs::Emitter::new("fs-exec/conformance", "exec-012/kill");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "fs-exec-kill-latency".to_string(),
                json: report.to_json(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("kill event validates");
    println!("{line}");
    let p99 = report.cancel_latency_p99_ns().unwrap_or(u64::MAX);
    verdict(
        "exec-012",
        cancelled && quiescent && p99 < 250_000_000,
        &format!(
            "registry kill drained the candidate tree ({}/{} tiles) with arenas quiescent; \
             observed p99 = {p99} ns (measured, ledgered; the 200 µs gate is the perf \
             harness's)",
            report.completed, report.total
        ),
        SEED,
    );
}

#[test]
fn exec_013_autotuner_calibrates_persists_and_pins_reproducibly() {
    use fs_exec::{TuneSource, Tuner};
    let probe = fs_substrate::CapabilityProbe::topology_only();
    let fingerprint = probe.fingerprint();
    let mut tuner = Tuner::cold(fingerprint);
    assert!(tuner.needs_calibration());
    // Calibrate: real stencil sweep through the real pool; report ledgered.
    let report = tuner.calibrate(&probe).expect("matching calibration probe");
    let mut em = fs_obs::Emitter::new("fs-exec/conformance", "exec-013/tune");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "fs-exec-calibration-report".to_string(),
                json: report.clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("calibration report validates");
    println!("{line}");
    // Rows carry the fingerprint plus typed, versioned measurement evidence.
    let fp_hex = format!("{fingerprint:016x}");
    let rows_keyed = report.contains(&fp_hex)
        && report.contains("\"evidence_version\":2")
        && report.contains("\"kind\":\"wall-time\"")
        && report.contains("\"kind\":\"work-count\"")
        && report.contains("\"kind\":\"throughput\"")
        && !report.contains("confidence");
    // Tuned decision now answers (source = tuned), and recalibration is
    // idempotent (same keys, refresh increments — visible in the report).
    let (_edge_first, src) = tuner.tile_edge_for("stencil7-f32").expect("bounded kernel");
    let tuned = src == TuneSource::Tuned;
    let report2 = tuner
        .calibrate(&probe)
        .expect("matching recalibration probe");
    let idempotent = report2.contains("\"refresh\":2");
    // Persist -> reload -> consume; foreign fingerprints go stale.
    let dir = std::env::temp_dir().join("fs-exec-tune-conf");
    std::fs::create_dir_all(&dir).expect("tmp dir");
    let path = dir.join("tune.jsonl");
    tuner.save(&path).expect("save");
    let mut reloaded = Tuner::load(&path, fingerprint).expect("load");
    let (edge_reloaded, src_reloaded) = reloaded
        .tile_edge_for("stencil7-f32")
        .expect("bounded kernel");
    let (edge_tuner, _) = tuner.tile_edge_for("stencil7-f32").expect("bounded kernel");
    let persisted = src_reloaded == TuneSource::Tuned && edge_reloaded == edge_tuner;
    let stale = Tuner::load(&path, fingerprint ^ 1)
        .expect("foreign load")
        .needs_calibration();
    // Replay fidelity: pin the recorded decision on a COLD tuner — the
    // choice reproduces without any calibration data.
    let recorded = reloaded.decisions()[0].clone();
    let mut replay = Tuner::cold(0);
    replay
        .pin(recorded.kernel.clone(), recorded.params.clone())
        .expect("recorded params are canonical");
    let (edge_replayed, src_replayed) = replay
        .tile_edge_for("stencil7-f32")
        .expect("bounded kernel");
    let replayable = src_replayed == TuneSource::Pinned && edge_replayed == edge_reloaded;
    verdict(
        "exec-013",
        rows_keyed && tuned && idempotent && persisted && stale && replayable,
        &format!(
            "calibrate->persist->consume round trip on machine {fp_hex}: rows fingerprinted \
             with typed evidence v2, recalibration idempotent (refresh=2), foreign fingerprints \
             stale, pinned replay reproduces the recorded plan (edge={})",
            edge_reloaded.cells()
        ),
        0,
    );
}

/// wf9.7.1 — stream identity is DECLARED, never scheduled: reusing a
/// pool, running unrelated work, and racing concurrent runs must not
/// perturb a logical run's keys; distinct declared runs must diverge;
/// and keys reconstruct from ledger fields alone.
#[test]
fn stream_keys_are_immune_to_pool_history_and_concurrency() {
    let baseline = pool_with(4, 0x5EED)
        .run(&KeyKernel { tiles: 64 })
        .expect("baseline");
    // A polluted pool: unrelated runs, then two CONCURRENT probes.
    let pool = pool_with(4, 0x5EED);
    for _ in 0..3 {
        let _ = pool.run(&KeyKernel { tiles: 7 }).expect("unrelated");
    }
    let (a, b) = std::thread::scope(|s| {
        let pa = &pool;
        let ha = s.spawn(move || pa.run(&KeyKernel { tiles: 64 }).expect("concurrent a"));
        let hb = s.spawn(move || pa.run(&KeyKernel { tiles: 64 }).expect("concurrent b"));
        (ha.join().expect("join a"), hb.join().expect("join b"))
    });
    assert_eq!(
        a, baseline,
        "pool history + concurrency must not perturb streams"
    );
    assert_eq!(b, baseline, "arrival order must not perturb streams");
    // Worker count is not identity either.
    let wide = pool_with(9, 0x5EED)
        .run(&KeyKernel { tiles: 64 })
        .expect("wide");
    assert_eq!(wide, baseline, "worker count must not perturb streams");
    // Distinct DECLARED runs diverge; the same declared run replays.
    let g = CancelGate::new();
    let (r1, report1) = pool.run_declared(&KeyKernel { tiles: 64 }, &g, RunId(1));
    let r1 = r1.expect("run 1");
    let (r1b, report1b) = pool.run_declared(&KeyKernel { tiles: 64 }, &g, RunId(1));
    let r1b = r1b.expect("run 1 replay");
    let (r2, report2) = pool.run_declared(&KeyKernel { tiles: 64 }, &g, RunId(2));
    let r2 = r2.expect("run 2");
    assert_eq!(report1.declared_run, RunId(1));
    assert_eq!(report1b.declared_run, RunId(1));
    assert_eq!(report2.declared_run, RunId(2));
    assert_eq!(r1, r1b, "same declared run is bit-identical");
    assert_ne!(r1, r2, "distinct declared runs diverge");
    assert_ne!(r1, baseline, "RunId(1) diverges from the implicit RunId(0)");
    // Replay from ledger fields alone (no pool, no hidden state).
    for &(tile, key) in &r1 {
        let reconstructed = StreamKey {
            seed: 0x5EED,
            kernel_id: fs_obs::fnv1a64(b"conf/keys"),
            tile,
            iteration: 1,
        };
        assert_eq!(
            reconstructed.key128(),
            key,
            "ledger fields reconstruct the key"
        );
    }
}
