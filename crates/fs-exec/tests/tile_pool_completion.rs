//! End-to-end completion-witness journeys for the real TilePool protocol.
//!
//! Run with:
//! `cargo test -p fs-exec --test tile_pool_completion -- --nocapture`
//!
//! Every terminal journey emits one deterministic JSONL evidence record. The
//! records intentionally contain semantic executor evidence, not wall-clock
//! latency, steal timing, scheduler replay, external-thread discovery, or
//! application-publication claims.

use core::ops::ControlFlow;
use fs_exec::{
    Budget, CancelGate, Cancelled, Cx, PoolConfig, Reduce, RunError, RunId, TileFailure,
    TileKernel, TilePlan, TilePool,
};

const BUILD_IDENTITY: &str = concat!("fs-exec-", env!("CARGO_PKG_VERSION"));

macro_rules! assert_common_completion {
    ($witness:expr) => {{
        assert_eq!(
            $witness.verify(),
            Ok(()),
            "executor-minted witness must verify: {}",
            $witness.to_canonical_json()
        );
        assert_eq!(
            $witness.joined_workers(),
            $witness.launched_workers(),
            "every entered worker must have exited before witness sealing"
        );
        assert!(
            $witness.worker_admission_closed(),
            "no worker can enter after witness sealing"
        );
        assert_eq!($witness.live_worker_guards_at_seal(), 0);
        assert_eq!($witness.live_tile_scopes_at_seal(), 0);
        assert!(
            $witness.executor_transients_quiescent(),
            "run-local worker, tile-scope, and root-charge transients must be closed"
        );
    }};
}

macro_rules! witnessed_parts {
    ($run:expr) => {{
        let bundle = $run.expect("executor must seal a valid witnessed-run bundle");
        bundle
            .verify_bundle()
            .expect("fresh executor bundle must verify");
        bundle.into_parts()
    }};
}

macro_rules! log_witness {
    ($case:expr, $sequence:expr, $witness:expr) => {{
        let line = $witness.to_jsonl_with_reuse($case, $sequence, BUILD_IDENTITY, None);
        assert!(
            !line.contains('\n'),
            "one journey must produce exactly one JSONL record"
        );
        assert!(line.contains("\"schema\":\"fs-exec-tilepool-completion-e2e-v2\""));
        assert!(
            line.len() < 16 * 1024,
            "focused evidence line must stay bounded"
        );
        assert!(line.contains("\"reuse_verdict\":null"));
        assert!(line.contains("\"witness_root\":\""));
        assert!(line.contains("\"no_claims\":["));
        println!("{line}");
        line
    }};
    ($case:expr, $sequence:expr, $witness:expr, reuse = $reuse:expr) => {{
        let line = $witness.to_jsonl_with_reuse($case, $sequence, BUILD_IDENTITY, Some($reuse));
        assert!(
            !line.contains('\n'),
            "one journey must produce exactly one JSONL record"
        );
        assert!(line.contains("\"schema\":\"fs-exec-tilepool-completion-e2e-v2\""));
        assert!(
            line.len() < 16 * 1024,
            "focused evidence line must stay bounded"
        );
        assert!(line.contains(concat!("\"reuse_verdict\":", stringify!($reuse))));
        assert!(line.contains("\"witness_root\":\""));
        assert!(line.contains("\"no_claims\":["));
        println!("{line}");
        line
    }};
}

fn pool(workers: usize) -> TilePool {
    TilePool::new(PoolConfig::for_host(workers, 0xC0_4D_50_1E_7E))
}

struct UnitKernel {
    name: &'static str,
    tiles: u64,
}

impl UnitKernel {
    const fn new(name: &'static str, tiles: u64) -> Self {
        Self { name, tiles }
    }
}

impl TileKernel for UnitKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new(self.name, self.tiles)
    }

    fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() {
            ControlFlow::Break(Cancelled)
        } else {
            ControlFlow::Continue(1)
        }
    }
}

struct CancelAfterFirstWave<'a> {
    gate: &'a CancelGate,
    rendezvous: std::sync::Barrier,
    tiles: u64,
}

impl TileKernel for CancelAfterFirstWave<'_> {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("completion/mid-run-cancel", self.tiles)
    }

    fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        self.rendezvous.wait();
        if tile == 0 {
            self.gate.request();
        }
        self.rendezvous.wait();
        ControlFlow::Continue(1)
    }
}

struct FaultAtZero {
    tiles: u64,
}

impl TileKernel for FaultAtZero {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("completion/typed-fault", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if tile == 0 {
            ControlFlow::Break(cx.refuse(TileFailure::InjectedFault {
                plan_version: 1,
                plan_seed: 0xFA_17,
                tiles: self.tiles,
                touches_per_tile: 1,
                touch: 1,
            }))
        } else {
            ControlFlow::Continue(1)
        }
    }
}

struct PanicAtZero {
    tiles: u64,
}

impl TileKernel for PanicAtZero {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("completion/panic", self.tiles)
    }

    fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        assert!(tile != 0, "completion witness panic fixture");
        ControlFlow::Continue(1)
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ReductionBomb(u64);

impl Reduce for ReductionBomb {
    fn identity() -> Self {
        Self(0)
    }

    fn merge(self, _other: Self) -> Self {
        panic!("completion witness reduction fixture");
    }
}

struct ReductionPanic;

impl TileKernel for ReductionPanic {
    type Out = ReductionBomb;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("completion/reduction-panic", 2)
    }

    fn run(&self, _tile: u64, _cx: &Cx<'_>) -> ControlFlow<Cancelled, ReductionBomb> {
        ControlFlow::Continue(ReductionBomb(1))
    }
}

struct ArenaAllocationRefusal;

impl TileKernel for ArenaAllocationRefusal {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("completion/arena-allocation-refusal", 1)
    }

    fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        match cx.arena().alloc_slice_fill(
            fs_alloc::Site::named("completion/arena-allocation-refusal"),
            1,
            0_u8,
        ) {
            Ok(_) => ControlFlow::Continue(1),
            Err(error) => ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
        }
    }
}

#[test]
fn success_zero_one_many_and_deterministic_identity_are_fully_accounted() {
    for (sequence, workers, tiles) in [(0, 1, 0), (1, 4, 1), (2, 4, 257)] {
        let p = pool(workers);
        let kernel = UnitKernel::new("completion/success-shapes", tiles);
        let (result, report, witness) = witnessed_parts!(p.run_declared_witnessed(
            &kernel,
            &CancelGate::new(),
            RunId(sequence)
        ));

        assert_eq!(result, Ok(tiles));
        assert_eq!(report.completed, tiles);
        assert_eq!(report.total, tiles);
        assert_common_completion!(witness);
        assert_eq!(witness.disposition_name(), "completed");
        assert_eq!(witness.scope_kind(), "std-thread-scope");
        assert_eq!(witness.parent_region_id(), None);
        assert_eq!(witness.parent_task_id(), None);
        assert!(witness.admission_completed());
        assert_eq!(witness.admitted_tiles(), tiles);
        assert_eq!(witness.unadmitted_tiles(), 0);
        assert_eq!(witness.claimed_tiles(), tiles);
        assert_eq!(witness.completed_tiles(), tiles);
        assert_eq!(witness.break_tiles(), 0);
        assert_eq!(witness.panicked_tiles(), 0);
        assert_eq!(witness.request_phase().name(), "not-requested");
        assert_eq!(witness.cancellation_observed_workers(), 0);
        assert_eq!(witness.failed_tiles(), 0);
        assert_eq!(witness.cancelled_tiles(), 0);
        assert_eq!(witness.tile_scopes_opened(), tiles);
        assert!(witness.root_charge_admitted());
        assert!(witness.root_charge_released());
        if tiles == 0 {
            assert_eq!(witness.planned_workers(), 1);
            assert_eq!(witness.launched_workers(), 1);
        }
        log_witness!("success-shape", sequence, witness);
        assert!(p.arena_pool().stats().quiescent());
    }

    let p = pool(4);
    let kernel = UnitKernel::new("completion/deterministic-root", 257);
    let gate_a = CancelGate::new();
    let (_, _, first) = witnessed_parts!(p.run_declared_witnessed(&kernel, &gate_a, RunId(77)));
    let gate_b = CancelGate::new();
    let (_, _, second) = witnessed_parts!(p.run_declared_witnessed(&kernel, &gate_b, RunId(77)));
    assert_eq!(first.plan_root_bytes(), second.plan_root_bytes());
    assert_eq!(
        first.call_replay_root_bytes(),
        second.call_replay_root_bytes(),
        "standalone replay roots are intentionally reproducible, not unique"
    );
    assert_eq!(first.root_bytes(), second.root_bytes());
    assert_eq!(first.to_canonical_json(), second.to_canonical_json());
    assert_eq!(
        first,
        first.clone(),
        "retained evidence is immutable value data"
    );

    let (_, _, different_run) =
        witnessed_parts!(p.run_declared_witnessed(&kernel, &CancelGate::new(), RunId(78)));
    assert_ne!(first.plan_root_bytes(), different_run.plan_root_bytes());
    assert_ne!(first.root_bytes(), different_run.root_bytes());

    let different_seed_pool = TilePool::new(PoolConfig::for_host(4, 0xC0_4D_50_1E_7F));
    let (_, _, different_seed) = witnessed_parts!(different_seed_pool.run_declared_witnessed(
        &kernel,
        &CancelGate::new(),
        RunId(77),
    ));
    assert_eq!(
        first.pool_placement_identity(),
        different_seed.pool_placement_identity(),
        "logical stream seed is deliberately not a placement dimension"
    );
    assert_ne!(
        first.plan_root_bytes(),
        different_seed.plan_root_bytes(),
        "the completion plan must still bind the exact logical stream seed"
    );

    let different_placement_pool = pool(2);
    let (_, _, different_placement) = witnessed_parts!(
        different_placement_pool.run_declared_witnessed(&kernel, &CancelGate::new(), RunId(77),)
    );
    assert_ne!(
        first.pool_placement_identity(),
        different_placement.pool_placement_identity()
    );
    assert_ne!(
        first.plan_root_bytes(),
        different_placement.plan_root_bytes()
    );

    let canonical = first.to_canonical_json();
    assert!(!canonical.contains("\"cancel_latencies_ns\""));
    assert!(!canonical.contains("\"steals\""));
    assert!(!canonical.contains("\"elapsed\""));
    assert!(canonical.contains("\"application-state-publication\""));
    assert!(canonical.contains("\"cross-call-uniqueness-without-affine-invocation-permit\""));
    assert!(!first.has_affine_invocation_permit());
    let envelope_a = first.to_jsonl("envelope-a", 10, "build-a");
    let envelope_b = first.to_jsonl("envelope-b", 11, "build-b");
    assert_ne!(envelope_a, envelope_b);
    assert_eq!(
        first.root_bytes(),
        second.root_bytes(),
        "logging-envelope metadata cannot alter witness identity"
    );
    log_witness!("deterministic-root", 3, first, reuse = true);
}

#[test]
fn prelaunch_and_mid_run_cancellation_seal_then_the_pool_reuses() {
    let p = pool(4);

    let pre_cancelled = CancelGate::new_clock_free();
    pre_cancelled.request();
    let (result, report, witness) = witnessed_parts!(p.run_declared_witnessed(
        &UnitKernel::new("completion/pre-cancelled", 64),
        &pre_cancelled,
        RunId(100),
    ));
    assert!(matches!(
        result,
        Err(RunError::Cancelled {
            completed: 0,
            total: 64,
            ..
        })
    ));
    assert_eq!(report.completed, 0);
    assert_common_completion!(witness);
    assert_eq!(witness.disposition_name(), "cancelled");
    assert_eq!(witness.request_phase().name(), "before-entry");
    assert!(witness.admission_completed());
    assert!(witness.cancellation_requested());
    assert!(witness.cancellation_observed_workers() > 0);
    assert_eq!(witness.claimed_tiles(), 0);
    assert_eq!(witness.completed_tiles(), 0);
    assert_eq!(witness.cancelled_tiles(), 64);
    let pre_cancelled_witness = witness;

    const WORKERS: usize = 4;
    const TILES: u64 = 16_384;
    let gate = CancelGate::new_clock_free();
    let kernel = CancelAfterFirstWave {
        gate: &gate,
        rendezvous: std::sync::Barrier::new(WORKERS),
        tiles: TILES,
    };
    let (result, report, witness) =
        witnessed_parts!(p.run_declared_witnessed(&kernel, &gate, RunId(101)));
    match result {
        Err(RunError::Cancelled {
            completed, total, ..
        }) => {
            assert_eq!(completed, WORKERS as u64);
            assert_eq!(total, TILES);
        }
        other => panic!("expected deterministic mid-run cancellation, got {other:?}"),
    }
    assert_eq!(report.completed, WORKERS as u64);
    assert_common_completion!(witness);
    assert_eq!(witness.claimed_tiles(), WORKERS as u64);
    assert_eq!(witness.completed_tiles(), WORKERS as u64);
    assert_eq!(witness.break_tiles(), 0);
    assert_eq!(witness.request_phase().name(), "before-terminal-decision");
    assert!(witness.cancellation_observed_workers() > 0);
    assert_eq!(witness.tile_scopes_opened(), WORKERS as u64);
    assert_eq!(witness.cancelled_tiles(), TILES - WORKERS as u64);
    log_witness!("pre-cancelled", 0, pre_cancelled_witness, reuse = true);

    assert_eq!(
        p.run(&UnitKernel::new("completion/reuse-after-cancel", 32))
            .expect("pool must survive cancellation"),
        32
    );
    log_witness!("mid-run-cancel", 1, witness, reuse = true);
    assert!(p.arena_pool().stats().quiescent());
}

#[test]
fn typed_fault_tile_panic_and_reduction_panic_are_distinct_and_reusable() {
    let p = pool(1);

    let (fault_result, _, fault_witness) =
        witnessed_parts!(p.run_witnessed(&FaultAtZero { tiles: 128 }));
    assert!(matches!(
        fault_result,
        Err(RunError::TileFailed {
            tile: 0,
            failure: TileFailure::InjectedFault { .. },
            completed: 0,
            ..
        })
    ));
    assert_common_completion!(fault_witness);
    assert_eq!(fault_witness.disposition_name(), "tile-failed");
    assert_eq!(fault_witness.claimed_tiles(), 1);
    assert_eq!(fault_witness.break_tiles(), 1);
    assert_eq!(fault_witness.panicked_tiles(), 0);
    assert!(fault_witness.cancellation_observed_workers() > 0);
    assert_eq!(fault_witness.retained_refusal_tiles(), 1);
    assert_eq!(fault_witness.failed_tiles(), 1);
    assert_eq!(fault_witness.cancelled_tiles(), 127);
    assert_eq!(fault_witness.first_failure_kind(), Some("tile-failed"));
    assert_eq!(fault_witness.first_failure_tile(), Some(0));
    assert_eq!(
        p.run(&UnitKernel::new("completion/reuse-after-fault", 8))
            .expect("pool survives typed refusal"),
        8
    );
    log_witness!("typed-fault", 0, fault_witness, reuse = true);

    let (panic_result, _, panic_witness) =
        witnessed_parts!(p.run_witnessed(&PanicAtZero { tiles: 128 }));
    assert!(matches!(
        panic_result,
        Err(RunError::TilePanicked {
            tile: 0,
            completed: 0,
            ..
        })
    ));
    assert_common_completion!(panic_witness);
    assert_eq!(panic_witness.disposition_name(), "tile-panicked");
    assert_eq!(panic_witness.claimed_tiles(), 1);
    assert_eq!(panic_witness.break_tiles(), 0);
    assert_eq!(panic_witness.panicked_tiles(), 1);
    assert!(panic_witness.cancellation_observed_workers() > 0);
    assert_eq!(panic_witness.failed_tiles(), 1);
    assert_eq!(panic_witness.cancelled_tiles(), 127);
    assert_eq!(panic_witness.first_failure_kind(), Some("tile-panicked"));
    assert_eq!(panic_witness.first_failure_tile(), Some(0));
    assert_eq!(
        p.run(&UnitKernel::new("completion/reuse-after-panic", 8))
            .expect("pool survives contained tile panic"),
        8
    );
    log_witness!("tile-panic", 1, panic_witness, reuse = true);

    let (reduction_result, _, reduction_witness) =
        witnessed_parts!(p.run_witnessed(&ReductionPanic));
    assert!(matches!(
        reduction_result,
        Err(RunError::ReductionPanicked { .. })
    ));
    assert_common_completion!(reduction_witness);
    assert_eq!(reduction_witness.disposition_name(), "reduction-panicked");
    assert_eq!(reduction_witness.claimed_tiles(), 2);
    assert_eq!(reduction_witness.completed_tiles(), 2);
    assert_eq!(reduction_witness.failed_tiles(), 0);
    assert_eq!(reduction_witness.cancellation_observed_workers(), 0);
    assert_eq!(
        p.run(&UnitKernel::new("completion/reuse-after-reduction", 8))
            .expect("pool survives contained reduction panic"),
        8
    );
    log_witness!("reduction-panic", 2, reduction_witness, reuse = true);
    assert!(p.arena_pool().stats().quiescent());
}

#[test]
fn lease_and_tile_arena_allocation_refusals_are_sealed_and_reusable() {
    let p = pool(2);
    let lease = fs_alloc::OperationMemoryLease::bounded(0);
    let gate = CancelGate::new_clock_free();
    let (result, report, witness) = witnessed_parts!(p.run_declared_leased_budgeted_witnessed(
        &UnitKernel::new("completion/root-lease-refusal", 64),
        &gate,
        RunId(200),
        Budget::INFINITE,
        &lease,
    ));
    match result {
        Err(RunError::MemoryRefused {
            requested_bytes,
            used_bytes,
            limit_bytes,
            ..
        }) => {
            assert!(requested_bytes > 0);
            assert_eq!(used_bytes, 0);
            assert_eq!(limit_bytes, 0);
        }
        other => panic!("expected root lease refusal, got {other:?}"),
    }
    assert_eq!(report.completed, 0);
    assert_common_completion!(witness);
    assert_eq!(witness.disposition_name(), "memory-refused");
    assert!(!witness.admission_completed());
    assert_eq!(witness.admitted_tiles(), 0);
    assert_eq!(witness.unadmitted_tiles(), 64);
    assert_eq!(witness.claimed_tiles(), 0);
    assert_eq!(witness.launched_workers(), 0);
    assert_eq!(witness.cancellation_observed_workers(), 0);
    assert!(!witness.root_charge_admitted());
    assert!(!witness.root_charge_released());
    assert_eq!(witness.lease_used_before(), 0);
    assert_eq!(witness.lease_used_after(), 0);
    assert_eq!(
        witness.lease_refusals_after(),
        witness.lease_refusals_before() + 1
    );
    assert_eq!(
        p.run(&UnitKernel::new("completion/reuse-after-lease-refusal", 8))
            .expect("pool survives prelaunch lease refusal"),
        8
    );
    log_witness!("root-lease-refusal", 0, witness, reuse = true);

    let mut config = PoolConfig::for_host(1, 0xA11_0C);
    config.arena.limit_bytes = Some(0);
    config.arena.free_list_max_bytes = 0;
    let allocation_pool = TilePool::new(config);
    let (result, _, witness) =
        witnessed_parts!(allocation_pool.run_witnessed(&ArenaAllocationRefusal));
    assert!(matches!(
        result,
        Err(RunError::TileFailed {
            tile: 0,
            failure: TileFailure::Allocation(_),
            completed: 0,
            ..
        })
    ));
    assert_common_completion!(witness);
    assert_eq!(witness.disposition_name(), "tile-failed");
    assert_eq!(witness.claimed_tiles(), 1);
    assert_eq!(witness.break_tiles(), 1);
    assert!(witness.cancellation_observed_workers() > 0);
    assert_eq!(witness.retained_refusal_tiles(), 1);
    assert_eq!(witness.failed_tiles(), 1);
    assert_eq!(
        allocation_pool
            .run(&UnitKernel::new("completion/reuse-after-arena-refusal", 8,))
            .expect("allocation-free work remains usable"),
        8
    );
    log_witness!("tile-arena-allocation-refusal", 1, witness, reuse = true);
    assert!(allocation_pool.arena_pool().stats().quiescent());
}

#[test]
fn parked_crew_witness_uses_the_same_join_authority_and_crew_reuses() {
    let p = pool(4);
    p.with_parked_crew_local(|parked| {
        let (result, report, witness) = witnessed_parts!(
            parked.run_witnessed(&UnitKernel::new("completion/parked-success", 257))
        );
        assert_eq!(result, Ok(257));
        assert_eq!(report.completed, 257);
        assert_common_completion!(witness);
        assert_eq!(witness.disposition_name(), "completed");
        assert_eq!(witness.scope_kind(), "std-thread-parked-crew");
        assert_eq!(witness.planned_workers(), 4);
        assert_eq!(witness.launched_workers(), 4);
        assert_eq!(witness.joined_workers(), 4);
        assert_eq!(witness.planned_crew_callbacks(), 4);
        assert_eq!(witness.entered_crew_callbacks(), 4);
        assert_eq!(witness.exited_crew_callbacks(), 4);
        assert_eq!(witness.claimed_tiles(), 257);

        let (short_result, short_report, short_witness) = witnessed_parts!(
            parked.run_witnessed(&UnitKernel::new("completion/parked-short-plan", 1))
        );
        assert_eq!(short_result, Ok(1));
        assert_eq!(short_report.completed, 1);
        assert_common_completion!(short_witness);
        assert_eq!(short_witness.planned_workers(), 1);
        assert_eq!(short_witness.launched_workers(), 1);
        assert_eq!(short_witness.joined_workers(), 1);
        assert_eq!(
            short_witness.planned_crew_callbacks(),
            4,
            "the complete parked callback set is distinct from the active worker subset"
        );
        assert_eq!(short_witness.entered_crew_callbacks(), 4);
        assert_eq!(short_witness.exited_crew_callbacks(), 4);
        log_witness!("parked-short-plan", 1, short_witness, reuse = true);

        assert_eq!(
            parked
                .run(&UnitKernel::new("completion/parked-reuse", 32))
                .expect("parked crew remains reusable"),
            32
        );
        log_witness!("parked-success", 0, witness, reuse = true);
    });
    assert!(p.arena_pool().stats().quiescent());
}
