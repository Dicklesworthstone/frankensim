//! Operation-memory-lease battery (bead wf9.16): one run-scoped lease
//! covers TilePool root metadata and every tile arena chunk, with atomic
//! reserve/release, canonical accounting, and refusals that drain cleanly.
//! Pool-cache history and near-limit concurrent admission can change peaks;
//! thread stacks, allocator overhead, and arbitrary kernel/output-owned heap
//! are explicit no-claims.

use core::ops::ControlFlow;

use fs_alloc::{AllocError, OperationMemoryLease, Site};
use fs_exec::{
    Budget, CancelGate, Cancelled, Cx, PoolConfig, RunError, RunId, TileFailure, TileKernel,
    TilePlan, TilePool,
};
use fs_substrate::affinity::CcdTopology;

fn small_chunk_pool(workers: usize, seed: u64) -> TilePool {
    let mut config = PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, seed);
    // Small chunks so per-tile lease charges are visible and cheap.
    config.arena.chunk_bytes = 4096;
    TilePool::new(config)
}

/// One arena chunk per tile; output is the tile id (sum-reduced).
struct AllocatingKernel {
    tiles: u64,
}

impl TileKernel for AllocatingKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("lease/allocating", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        match cx
            .arena()
            .alloc_slice_fill(Site::named("lease-battery/tile"), 512, tile as u8)
        {
            Ok(_) => ControlFlow::Continue(tile),
            Err(error) => ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
        }
    }
}

/// No arena traffic at all: the lease sees exactly the root metadata.
struct RootOnlyKernel {
    tiles: u64,
}

impl TileKernel for RootOnlyKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("lease/root-only", self.tiles)
    }

    fn run(&self, tile: u64, _cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        ControlFlow::Continue(tile)
    }
}

struct PanickingKernel;

impl TileKernel for PanickingKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("lease/panicking", 2)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        let _ = cx
            .arena()
            .alloc_slice_fill(Site::named("lease-battery/panic"), 512, 0u8);
        assert!(tile != 0, "tile 0 panics with its chunk held");
        ControlFlow::Continue(tile)
    }
}

struct SelfCancellingKernel;

impl TileKernel for SelfCancellingKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("lease/self-cancel", 4)
    }

    fn run(&self, _tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        let _ = cx
            .arena()
            .alloc_slice_fill(Site::named("lease-battery/cancel"), 512, 0u8);
        ControlFlow::Break(Cancelled)
    }
}

fn run_leased(
    pool: &TilePool,
    kernel: &impl TileKernel<Out = u64>,
    lease: &OperationMemoryLease,
) -> (Result<u64, RunError>, fs_exec::RunReport) {
    pool.run_declared_leased_budgeted(
        kernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        lease,
    )
}

#[test]
fn exact_accounting_is_deterministic_and_returns_to_zero() {
    let pool = small_chunk_pool(2, 0x1EA5E);
    let kernel = RootOnlyKernel { tiles: 4 };
    let first = OperationMemoryLease::unbounded();
    let second = OperationMemoryLease::unbounded();
    assert_eq!(run_leased(&pool, &kernel, &first).0.expect("run 1"), 6);
    assert_eq!(run_leased(&pool, &kernel, &second).0.expect("run 2"), 6);
    let (a, b) = (first.receipt(), second.receipt());
    assert_eq!(a, b, "identical plans must produce identical receipts");
    assert!(a.requested_bytes > 0, "root metadata is charged");
    assert_eq!(a.used_bytes, 0, "every charge is released by run end");
    assert_eq!(a.refusals, 0);
    assert_eq!(a.peak_bytes, a.requested_bytes, "one root charge, no churn");
}

#[test]
fn root_metadata_refusal_precedes_launch_and_pool_stays_usable() {
    let pool = small_chunk_pool(2, 0x1EA5F);
    let kernel = AllocatingKernel { tiles: 4 };
    let lease = OperationMemoryLease::bounded(1);
    let (result, report) = run_leased(&pool, &kernel, &lease);
    match result {
        Err(RunError::MemoryRefused {
            what: "tilepool-root-metadata",
            used_bytes: 0,
            limit_bytes: 1,
            ..
        }) => {}
        other => panic!("expected the pre-launch refusal, got {other:?}"),
    }
    assert_eq!(report.completed, 0);
    assert_eq!(lease.receipt().used_bytes, 0);
    assert_eq!(lease.receipt().refusals, 1);
    // Nothing launched, nothing leaked; the pool remains fully usable.
    assert!(pool.arena_pool().stats().quiescent());
    let retry = OperationMemoryLease::unbounded();
    assert_eq!(run_leased(&pool, &kernel, &retry).0.expect("retry"), 6);
}

#[test]
fn mid_run_chunk_refusal_drains_and_releases_exactly() {
    let pool = small_chunk_pool(1, 0x1EA60);
    let kernel = AllocatingKernel { tiles: 2 };
    // The lease gates the CONCURRENT live set: with one worker, tile 0's
    // chunk is released back before tile 1 acquires, so the cumulative
    // demand admits sequentially. Starve the PEAK instead — one byte under
    // the concurrency high-water refuses the first chunk deterministically.
    let probe = OperationMemoryLease::unbounded();
    assert_eq!(run_leased(&pool, &kernel, &probe).0.expect("probe"), 1);
    let peak_demand = probe.receipt().peak_bytes;
    let starved = OperationMemoryLease::bounded(peak_demand - 1);
    let (result, report) = run_leased(&pool, &kernel, &starved);
    match result {
        Err(RunError::TileFailed {
            tile: 0,
            failure: TileFailure::Allocation(AllocError::LeaseExhausted { .. }),
            completed: 0,
            ..
        }) => {}
        other => panic!("expected the first tile's chunk to be refused, got {other:?}"),
    }
    assert_eq!(report.completed, 0);
    let receipt = starved.receipt();
    assert_eq!(receipt.used_bytes, 0, "drain must release every charge");
    assert_eq!(receipt.refusals, 1);
    assert_eq!(
        receipt.first_refusal.as_ref().expect("recorded").what,
        "arena-chunk"
    );
    assert!(pool.arena_pool().stats().quiescent());
}

#[test]
fn panic_releases_all_charges_and_contains_the_unwind() {
    let pool = small_chunk_pool(2, 0x1EA61);
    let lease = OperationMemoryLease::bounded(1 << 20);
    let (result, _report) = pool.run_declared_leased_budgeted(
        &PanickingKernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        &lease,
    );
    assert!(
        matches!(result, Err(RunError::TilePanicked { tile: 0, .. })),
        "got {result:?}"
    );
    assert_eq!(
        lease.receipt().used_bytes,
        0,
        "an unwinding tile must release its chunk charge and the root charge"
    );
    assert!(pool.arena_pool().stats().quiescent());
}

#[test]
fn cancellation_releases_all_charges() {
    let pool = small_chunk_pool(2, 0x1EA62);
    let lease = OperationMemoryLease::bounded(1 << 20);
    let (result, _report) = pool.run_declared_leased_budgeted(
        &SelfCancellingKernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        &lease,
    );
    assert!(
        matches!(result, Err(RunError::Cancelled { .. })),
        "got {result:?}"
    );
    let receipt = lease.receipt();
    assert_eq!(receipt.used_bytes, 0);
    assert_eq!(receipt.refusals, 0);
    assert!(pool.arena_pool().stats().quiescent());
}

#[test]
fn recycled_chunks_charge_each_operation_exactly_once() {
    let pool = small_chunk_pool(2, 0x1EA63);
    let kernel = AllocatingKernel { tiles: 4 };
    let first = OperationMemoryLease::unbounded();
    assert_eq!(run_leased(&pool, &kernel, &first).0.expect("run 1"), 6);
    let recycled_before = pool.arena_pool().stats().chunks_recycled;
    let second = OperationMemoryLease::unbounded();
    assert_eq!(run_leased(&pool, &kernel, &second).0.expect("run 2"), 6);
    let stats = pool.arena_pool().stats();
    // A later run may reach a higher concurrent live set and legitimately
    // grow the process-wide pool. The invariant under test is that this
    // operation actually acquired cached storage and received the same
    // logical lease charge as a fresh acquisition.
    assert!(
        stats.chunks_recycled > recycled_before,
        "run 2 must acquire at least one cached chunk: recycled before={recycled_before}, after={stats:?}"
    );
    assert_eq!(
        first.receipt().requested_bytes,
        second.receipt().requested_bytes,
        "a recycled chunk charges the acquiring operation exactly what a fresh one does"
    );
    assert_eq!(second.receipt().used_bytes, 0);
    assert!(
        stats.quiescent(),
        "both runs must drain every arena: {stats:?}"
    );
}

#[test]
fn concurrent_admission_stays_under_the_limit() {
    let pool = small_chunk_pool(4, 0x1EA64);
    let kernel = AllocatingKernel { tiles: 16 };
    // Generous but bounded: all concurrent workers' chunks must fit, and
    // the peak must respect the limit under CAS contention.
    let lease = OperationMemoryLease::bounded(1 << 22);
    let (result, report) = run_leased(&pool, &kernel, &lease);
    assert_eq!(result.expect("bounded run succeeds"), (0..16).sum::<u64>());
    assert_eq!(report.completed, 16);
    let receipt = lease.receipt();
    assert_eq!(receipt.used_bytes, 0);
    assert_eq!(receipt.refusals, 0);
    assert!(receipt.peak_bytes <= 1 << 22);
    assert!(receipt.peak_bytes > 0);
}

#[test]
fn varying_worker_counts_hold_the_lease_invariants() {
    for workers in [1usize, 2, 8] {
        let pool = small_chunk_pool(workers, 0x1EA65 + workers as u64);
        let kernel = AllocatingKernel { tiles: 8 };
        let lease = OperationMemoryLease::bounded(1 << 22);
        let (result, _report) = run_leased(&pool, &kernel, &lease);
        assert_eq!(
            result.expect("run"),
            (0..8).sum::<u64>(),
            "workers={workers}"
        );
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0, "workers={workers}");
        assert!(receipt.peak_bytes <= receipt.requested_bytes);
        assert!(pool.arena_pool().stats().quiescent());
    }
}

/// wf9.16.1: a list-shaped output flows through the leased path with its
/// payload admitted BEFORE allocation via `cx.lease()`, folds through
/// `Concat`'s detached identity, and every payload charge releases with
/// the results' drop — the full-run live-set ceiling now covers output
/// payloads, not just headers.
struct AdmittedListKernel {
    tiles: u64,
}

impl TileKernel for AdmittedListKernel {
    type Out = fs_exec::Concat<u64>;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("lease/admitted-list", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, Self::Out> {
        let lease = cx.lease().expect("pool runs carry the operation lease");
        let mut out = match fs_alloc::LeasedVec::with_capacity(lease, "lease-battery/list", 4) {
            Ok(out) => out,
            Err(error) => return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error))),
        };
        for i in 0..4u64 {
            if let Err(error) = out.push(tile * 4 + i) {
                return ControlFlow::Break(cx.refuse(TileFailure::Allocation(error)));
            }
        }
        ControlFlow::Continue(fs_exec::Concat(out))
    }
}

#[test]
fn admitted_list_outputs_carry_their_payload_charges_through_the_fold() {
    let pool = small_chunk_pool(2, 0x1EA66);
    let kernel = AdmittedListKernel { tiles: 4 };
    let lease = OperationMemoryLease::bounded(1 << 20);
    let (result, report) = pool.run_declared_leased_budgeted(
        &kernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        &lease,
    );
    let merged = result.expect("leased list run");
    assert_eq!(report.completed, 4);
    // Deterministic fold order: tiles ascending, each contributing 4
    // consecutive values.
    assert_eq!(merged.0.as_slice(), (0..16).collect::<Vec<_>>().as_slice());
    // The merged output's payload is still lease-charged while owned...
    assert!(lease.receipt().used_bytes >= 16 * 8);
    drop(merged);
    // ...and releases exactly when the caller drops it.
    assert_eq!(lease.receipt().used_bytes, 0);
    assert!(pool.arena_pool().stats().quiescent());
}

/// wf9.16.1: an output-payload admission refusal mid-run is a typed tile
/// failure that drains and releases, exactly like an arena-chunk refusal.
#[test]
fn output_payload_refusal_drains_and_releases() {
    let pool = small_chunk_pool(1, 0x1EA67);
    let kernel = AdmittedListKernel { tiles: 2 };
    // Fit the root metadata but starve the first tile's payload+chunk.
    let probe = OperationMemoryLease::unbounded();
    let _ = pool.run_declared_leased_budgeted(
        &kernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        &probe,
    );
    let starved = OperationMemoryLease::bounded(probe.receipt().peak_bytes - 1);
    let (result, _report) = pool.run_declared_leased_budgeted(
        &kernel,
        &CancelGate::new(),
        RunId(7),
        Budget::INFINITE,
        &starved,
    );
    // One byte under the peak refuses wherever the peak lives: a tile's
    // payload/chunk admission (typed TileFailed) or the fold's merge
    // re-admission (the documented ReductionPanicked containment). Both
    // drain and both release every charge.
    match result {
        Err(
            RunError::TileFailed {
                failure: TileFailure::Allocation(AllocError::LeaseExhausted { .. }),
                ..
            }
            | RunError::ReductionPanicked { .. },
        ) => {}
        other => panic!("expected a lease refusal on the starved run, got {other:?}"),
    }
    assert_eq!(
        starved.receipt().used_bytes,
        0,
        "every charge releases whether the refusal hit a tile or the fold"
    );
    assert!(pool.arena_pool().stats().quiescent());
}
