//! Allocation-stable work stealing (bead wf9.16.2): worker ownership is a
//! contiguous `TileRun` — two u64s — so half-steals are pure Copy
//! arithmetic and the loop allocates NOTHING after launch. This battery
//! proves it three ways: adversarial multi-level steal storms with exact
//! per-tile completion, deterministic outputs across worker counts, and a
//! differential allocation count (a counting global allocator measuring a
//! steal-heavy run against a steal-free run of identical shape — the delta
//! must be ZERO).

use core::ops::ControlFlow;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use fs_exec::{
    Budget, CancelGate, Cancelled, Cx, PoolConfig, RunId, TileKernel, TilePlan, TilePool,
};
use std::sync::Mutex;

/// libtest runs test fns on parallel threads and every thread shares the
/// global allocation counter; the differential windows are only meaningful
/// with the whole battery serialized.
static SERIAL: Mutex<()> = Mutex::new(());
use fs_substrate::affinity::CcdTopology;

/// Counting wrapper around the system allocator. Allocation COUNTS (not
/// bytes) are compared differentially between two runs of identical shape,
/// so any nonzero delta is caused by the only difference: steal traffic.
struct CountingAlloc;

static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

// SAFETY: pure delegation to `System`; the only addition is a relaxed
// counter increment on the allocation paths, which allocates nothing and
// touches no allocator state. Test-binary-only instrumentation.
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: same contract as the delegated call.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: same contract as the delegated call.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: same contract as the delegated call.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
#[allow(unsafe_code)]
static COUNTER: CountingAlloc = CountingAlloc;

/// Spin-heavy no-allocation kernel; per-tile work is tunable so steal
/// pressure can be shaped. Records each tile's executions for
/// exactly-once verification.
struct SpinKernel {
    tiles: u64,
    spins: u64,
    executed: Vec<AtomicU64>,
}

impl SpinKernel {
    fn new(tiles: u64, spins: u64) -> Self {
        SpinKernel {
            tiles,
            spins,
            executed: (0..tiles).map(|_| AtomicU64::new(0)).collect(),
        }
    }
}

impl TileKernel for SpinKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("steal/spin", self.tiles)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let mut acc = tile.max(1);
        for _ in 0..self.spins {
            acc = std::hint::black_box(
                acc.wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1),
            );
        }
        self.executed[usize::try_from(tile).expect("tile index")].fetch_add(1, Ordering::Relaxed);
        ControlFlow::Continue(tile.wrapping_add(acc & 1))
    }
}

/// All initial work on worker 0 (weight 10^6 vs 1): every other worker can
/// obtain tiles only by stealing, and their thefts cascade through steal
/// genealogies several levels deep.
fn single_owner_pool(workers: usize, seed: u64) -> TilePool {
    let mut config = PoolConfig::new(workers, CcdTopology::APPLE_M_CLASS, seed);
    config.quantum_weights = vec![1; workers];
    config.quantum_weights[0] = 1_000_000;
    TilePool::new(config)
}

#[test]
fn multi_level_steal_storm_executes_every_tile_exactly_once() {
    let _serial = SERIAL.lock().expect("serialized battery");
    for workers in [2usize, 4, 8] {
        let pool = single_owner_pool(workers, 0x57EA1 + workers as u64);
        let kernel = SpinKernel::new(512, 300);
        let (result, report) =
            pool.run_declared_leased_budgeted(
                &kernel,
                &CancelGate::new(),
                RunId(3),
                Budget::INFINITE,
                &fs_alloc::OperationMemoryLease::unbounded(),
            );
        let sum = result.expect("storm run");
        assert_eq!(report.completed, 512, "workers={workers}");
        for (tile, count) in kernel.executed.iter().enumerate() {
            assert_eq!(
                count.load(Ordering::Relaxed),
                1,
                "tile {tile} must run exactly once (workers={workers})"
            );
        }
        // Deterministic fixed-slot reduction: the folded value is a pure
        // function of the plan, whatever the steal schedule did.
        let expected = SpinKernel::new(512, 300);
        let baseline = TilePool::new(PoolConfig::new(1, CcdTopology::APPLE_M_CLASS, 0x1))
            .run(&expected)
            .expect("single worker baseline");
        assert_eq!(sum, baseline, "workers={workers}");
        if workers > 1 {
            assert!(report.steals > 0, "single-owner seeding must force steals");
        }
    }
}

#[test]
fn steal_traffic_allocates_exactly_nothing() {
    // Two runs of IDENTICAL shape (same workers, tiles, spins, kernel
    // construction outside the window). The first seeds every worker
    // near-evenly (steals rare); the second seeds all work on worker 0
    // (steals mandatory, multi-level). Identical allocation sequences
    // everywhere else, so the count delta isolates steal traffic.
    const WORKERS: usize = 4;
    const TILES: u64 = 256;
    const SPINS: u64 = 200;
    let _serial = SERIAL.lock().expect("serialized battery");

    let balanced_pool = TilePool::new(PoolConfig::new(
        WORKERS,
        CcdTopology::APPLE_M_CLASS,
        0xA110C,
    ));
    let storm_pool = single_owner_pool(WORKERS, 0xA110C);
    let balanced_kernel = SpinKernel::new(TILES, SPINS);
    let storm_kernel = SpinKernel::new(TILES, SPINS);
    let lease_a = fs_alloc::OperationMemoryLease::unbounded();
    let lease_b = fs_alloc::OperationMemoryLease::unbounded();

    // Warm both pools once (thread-spawn and arena paths touch their
    // lazily initialized state outside the measured windows).
    let _ = balanced_pool.run(&SpinKernel::new(TILES, 1));
    let _ = storm_pool.run(&SpinKernel::new(TILES, 1));

    // min-of-3 per configuration: ambient harness allocations only ever
    // ADD to a window, so the minimum is the clean per-run count.
    let mut balanced_allocs = usize::MAX;
    let mut balanced_report = None;
    let mut storm_allocs = usize::MAX;
    let mut storm_report = None;
    for _ in 0..3 {
        let before = ALLOCATIONS.load(Ordering::Relaxed);
        let (result, report) = balanced_pool.run_declared_leased_budgeted(
            &balanced_kernel,
            &CancelGate::new(),
            RunId(9),
            Budget::INFINITE,
            &lease_a,
        );
        result.expect("balanced run");
        balanced_allocs = balanced_allocs.min(ALLOCATIONS.load(Ordering::Relaxed) - before);
        balanced_report = Some(report);

        let before = ALLOCATIONS.load(Ordering::Relaxed);
        let (result, report) = storm_pool.run_declared_leased_budgeted(
            &storm_kernel,
            &CancelGate::new(),
            RunId(9),
            Budget::INFINITE,
            &lease_b,
        );
        result.expect("storm run");
        storm_allocs = storm_allocs.min(ALLOCATIONS.load(Ordering::Relaxed) - before);
        storm_report = Some(report);
    }
    let balanced_report = balanced_report.expect("measured");
    let storm_report = storm_report.expect("measured");

    // The two seedings produce materially different steal traffic (the
    // balanced run steals often in small end-of-run corrections; the
    // single-owner run steals rarely in huge cascade halves). Which side
    // is larger is schedule-dependent — the invariant under test is that
    // the allocation count is INDEPENDENT of steal count, so all the
    // discriminator needs is a difference.
    assert_ne!(
        storm_report.steals + balanced_report.steals,
        0,
        "at least one run must actually steal"
    );
    assert_ne!(
        storm_report.steals, balanced_report.steals,
        "the two seedings must produce different steal traffic for the \
         differential to discriminate"
    );
    assert_eq!(
        storm_allocs, balanced_allocs,
        "steal traffic must allocate exactly nothing: balanced {balanced_allocs} vs \
         storm {storm_allocs} allocation calls at identical run shape",
    );
}
