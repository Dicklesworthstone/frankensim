//! fs-exec conformance suite (CONTRACT.md: any reimplementation must pass).
//!
//! Gauntlet coverage: G0 (completeness, reduction laws), G4 (cancellation
//! storm with panic injection, drain/leak oracles), G5 (bit-identical
//! results and stream keys across worker counts). Latency measurements are
//! emitted as events — measured, never assumed. Every case prints one
//! JSON-line verdict; seeded cases carry their seed.

use core::ops::ControlFlow;
use fs_exec::{
    CancelGate, Cancelled, Cx, PoolConfig, RunError, TileKernel, TilePlan, TilePool, victim_order,
    weighted_ranges,
};
use fs_substrate::affinity::CcdTopology;

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-exec/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
    let mut checked = 0;
    for workers in [1usize, 2, 5, 8] {
        for tiles in [0u64, 1, 7, 128] {
            let pool = pool_with(workers, 0xE001);
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
    );
}

#[test]
fn exec_002_g5_bit_identical_results_across_worker_counts() {
    let reference = pool_with(1, 0xE002)
        .run(&FloatKernel { tiles: 257 })
        .expect("reference run");
    let mut all_equal = true;
    for workers in [2usize, 3, 4, 8] {
        for _ in 0..3 {
            let got = pool_with(workers, 0xE002)
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
    );
}

#[test]
fn exec_003_stream_keys_are_worker_independent() {
    let mut runs: Vec<Vec<(u64, u128)>> = Vec::new();
    for workers in [1usize, 2, 7] {
        // Fresh pool per run: iteration counters all start at 0, so the
        // ONLY varying factor is worker count — which must not matter.
        let got = pool_with(workers, 0xE003)
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
    );
}

#[test]
fn exec_004_external_cancellation_drains_and_ledgers_latency() {
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
    let pool = pool_with(4, 0xE004);
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
    );
}

#[test]
fn exec_007_latency_lane_stays_responsive_under_tile_load() {
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
        0xE007,
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
    );
}

#[test]
fn exec_008_reduction_shape_is_scheduling_invariant() {
    // A non-commutative reduction (concatenation) is the sharpest witness:
    // ANY arrival-order fold would scramble it under stealing.
    for workers in [1usize, 2, 4, 8] {
        for _ in 0..5 {
            let got = pool_with(workers, 0xE008)
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
    );
}
