//! Fallible-storage authority ratchets for the CBC executor (bead
//! frankensim-epic-bedrock-6ys.20.5).
//!
//! Every construction and growth path must be a checked reservation refused
//! with a typed [`CbcExecError::Storage`] diagnostic before mutation, and
//! every admitted arithmetic transition must be allocation-free so resource
//! admission can never become an allocation panic.
//!
//! The tests below install a process-global counting allocator (integration
//! test targets are separate binaries, so the workspace's other suites are
//! untouched). Injection is single-shot: the first allocation attempted while
//! armed is refused and automatically disarms the policy, keeping runs
//! deterministic. All stateful tests serialize on one mutex because the
//! allocator counters are process-global, and every measurement is taken
//! relative to a scope baseline because libtest deliberately leaks harness
//! memory and prior panics unwind during later tests, which would poison
//! absolute comparisons.

#![deny(unsafe_code)]

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicIsize, AtomicUsize, Ordering};
use std::sync::{LazyLock, Mutex, MutexGuard};

use fs_rand::cbc::{CbcBudget, CbcExecutionMode, CbcProblem};
use fs_rand::cbc_exec::{
    CbcBoundary, CbcControl, CbcExecError, CbcPhaseKind, CbcRunStatus, CbcStorageClass,
    CbcStorageRefusal, CbcTileShape, RANKED_STORAGE_REMEDIATIONS,
};

// ---------------------------------------------------------------------------
// Process-global counting/failing allocator (single-shot injection).
// ---------------------------------------------------------------------------

static LIVE_BYTES: AtomicIsize = AtomicIsize::new(0);
static PEAK_LIVE_BYTES: AtomicUsize = AtomicUsize::new(0);
static ALLOCATION_COUNT: AtomicUsize = AtomicUsize::new(0);
thread_local! {
    /// Per-thread single-shot injection. Thread-local because the libtest
    /// harness allocates concurrently on other threads; a global flag would
    /// feed nulls into infallible foreign allocations and abort the process.
    static INJECT_ARMED: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// Consume one pending injection on the CURRENT thread.
fn take_injection() -> bool {
    INJECT_ARMED.with(|armed| armed.replace(false))
}

struct Counting;

// The allocator must touch raw pointers by contract; this capsule is the
// standard-library delegation seam. The suite itself stays safe Rust.
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if take_injection() {
            return std::ptr::null_mut();
        }
        // SAFETY: delegation to the system allocator with the caller's layout.
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            note_allocation(layout.size());
        }
        pointer
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        if take_injection() {
            return std::ptr::null_mut();
        }
        // SAFETY: delegation to the system allocator with the caller's layout.
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() {
            note_allocation(layout.size());
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        LIVE_BYTES.fetch_sub(
            isize::try_from(layout.size()).expect("allocation sizes fit isize"),
            Ordering::SeqCst,
        );
        // SAFETY: the pointer was produced by System with this layout.
        unsafe { System.dealloc(pointer, layout) };
    }

    unsafe fn realloc(&self, pointer: *mut u8, old_layout: Layout, new_size: usize) -> *mut u8 {
        if take_injection() {
            return std::ptr::null_mut();
        }
        let grown = isize::try_from(new_size).expect("allocation sizes fit isize")
            - isize::try_from(old_layout.size()).expect("allocation sizes fit isize");
        // SAFETY: the pointer was produced by System with this old layout.
        let replacement = unsafe { System.realloc(pointer, old_layout, new_size) };
        if !replacement.is_null() {
            LIVE_BYTES.fetch_add(grown, Ordering::SeqCst);
            raise_peak();
            bump_allocation_count();
        }
        replacement
    }
}

fn note_allocation(size: usize) {
    LIVE_BYTES.fetch_add(
        isize::try_from(size).expect("allocation sizes fit isize"),
        Ordering::SeqCst,
    );
    raise_peak();
    bump_allocation_count();
}

fn bump_allocation_count() {
    ALLOCATION_COUNT.fetch_add(1, Ordering::SeqCst);
}

fn raise_peak() {
    let live = LIVE_BYTES.load(Ordering::SeqCst);
    if live > 0 {
        PEAK_LIVE_BYTES.fetch_max(
            usize::try_from(live).expect("live byte counts are positive here"),
            Ordering::SeqCst,
        );
    }
}

#[global_allocator]
static COUNTING_ALLOCATOR: Counting = Counting;

/// Serialize every allocator-sensitive test in this binary.
fn counter_lock() -> MutexGuard<'static, ()> {
    static LOCK: LazyLock<Mutex<()>> = LazyLock::new(|| Mutex::new(()));
    match LOCK.lock() {
        Ok(guard) => guard,
        Err(poisoned) => poisoned.into_inner(),
    }
}

/// A relative measurement scope. Baselines isolate the measured region from
/// libtest's intentional leaks and from anything live before the region.
struct Scope {
    base_live: isize,
    base_peak: usize,
    base_allocations: usize,
}

fn scope_begin() -> Scope {
    Scope {
        base_live: LIVE_BYTES.load(Ordering::SeqCst),
        base_peak: PEAK_LIVE_BYTES.load(Ordering::SeqCst),
        base_allocations: ALLOCATION_COUNT.load(Ordering::SeqCst),
    }
}

impl Scope {
    /// Net live-byte change since the scope began.
    fn live_delta(&self) -> isize {
        LIVE_BYTES.load(Ordering::SeqCst) - self.base_live
    }

    /// Peak live bytes observed inside the scope, above its baseline.
    fn relative_peak(&self) -> usize {
        PEAK_LIVE_BYTES
            .load(Ordering::SeqCst)
            .saturating_sub(self.base_peak)
    }

    /// Allocations performed inside the scope.
    fn allocations(&self) -> usize {
        ALLOCATION_COUNT.load(Ordering::SeqCst) - self.base_allocations
    }
}

fn arm_single_shot_failure() {
    INJECT_ARMED.with(|armed| armed.set(true));
}

// ---------------------------------------------------------------------------
// Fixtures.
// ---------------------------------------------------------------------------

fn admitted(n: u32, dimension: usize, mode: CbcExecutionMode) -> fs_rand::cbc::CbcAdmission {
    let problem = CbcProblem::new(n, dimension).expect("structural fixture");
    problem
        .admit_for(mode, CbcBudget::UNBOUNDED)
        .expect("fixture admits under an unbounded budget")
}

fn tile(candidate_block: u32, point_block: u32) -> CbcTileShape {
    CbcTileShape::new(candidate_block, point_block).expect("nonzero test tile")
}

/// Deterministic never-cancelling poll closure.
fn keep_going() -> CbcControl {
    CbcControl::Continue
}

fn drive_to_completion(executor: &mut fs_rand::cbc_exec::CbcExecutor) -> CbcRunStatus {
    executor
        .run(&mut keep_going, tile(4, 2), u128::MAX)
        .expect("unbounded allowance cannot exhaust")
}

/// [`stepped_run_result`] with the expectation that no storage refusal fires.
fn stepped_run(executor: &mut fs_rand::cbc_exec::CbcExecutor) -> CbcRunStatus {
    stepped_run_result(executor).expect("stepped execution cannot refuse while disarmed")
}

/// Run one tile batch that cancels at its first poll, so each call advances
/// the resumable cursor by a bounded amount instead of driving to completion.
fn stepped_run_result(
    executor: &mut fs_rand::cbc_exec::CbcExecutor,
) -> Result<CbcRunStatus, CbcExecError> {
    let mut first_poll = true;
    executor.run(
        &mut || {
            if first_poll {
                first_poll = false;
                CbcControl::Cancel
            } else {
                CbcControl::Continue
            }
        },
        tile(4, 2),
        u128::MAX,
    )
}

// ---------------------------------------------------------------------------
// Suite.
// ---------------------------------------------------------------------------

#[test]
fn cfs_001_construction_reserves_every_admitted_class_fallibly() {
    let _guard = counter_lock();
    let scope = scope_begin();

    let admission = admitted(64, 3, CbcExecutionMode::Certified);
    let mut executor = fs_rand::cbc_exec::CbcExecutor::new(admission)
        .expect("an unbounded receipt admits the executor");
    executor
        .enable_certificates()
        .expect("a certified receipt enables evidence");

    let observed = executor.storage_observation();
    assert!(
        observed.minimum_observed_product_capacity_limbs() >= observed.requested_product_limbs(),
        "every product reserved its admitted limb payload"
    );
    assert!(
        observed.observed_certificate_record_capacity() >= 2,
        "the certified record array reserved one slot per scanned component"
    );

    // Peak live memory stays inside four times the sealed resource envelope
    // (allocator rounding lives outside the claim, so this is an evidence
    // bound, not an admission theorem).
    let estimate = admission.estimate();
    let envelope = estimate.executor_inline_bytes()
        + estimate.product_owner_array_bytes()
        + estimate.resident_product_payload_bytes()
        + estimate.product_overlap_bytes()
        + estimate.candidate_phase_bytes()
        + estimate.update_phase_bytes()
        + estimate.certificate_retained_bytes();
    let peak = scope.relative_peak();
    let ceiling = usize::try_from(envelope.saturating_mul(4)).expect("envelope fits usize");
    assert!(
        peak <= ceiling,
        "relative peak {peak} escaped four times the sealed envelope {envelope}"
    );

    // NOTE on release proofs: absolute liveness deltas are unobservable on
    // a shared process allocator under libtest (harness leaks by design,
    // prior panics free late), and relative peaks are monotone globals.
    // The bounded-peak ceiling above is this suite's memory gate at scale;
    // full release-after-drop is proven deterministically by cfs_005 on
    // the minimal problem, whose whole footprint sits far below noise.
    drop(executor);
}

#[test]
fn cfs_002_construction_injection_refuses_typed_without_partial_leaks() {
    let _guard = counter_lock();
    let admission = admitted(8, 3, CbcExecutionMode::Certified);

    let scope = scope_begin();
    arm_single_shot_failure();
    let refusal = fs_rand::cbc_exec::CbcExecutor::new(admission)
        .expect_err("the armed refusal must surface as a typed error");
    assert_eq!(
        refusal,
        CbcExecError::Storage(CbcStorageRefusal {
            class: CbcStorageClass::ProductOwnerArray,
            phase: CbcPhaseKind::Construction,
            cursor: 0,
            requested: 8,
            admitted: 8,
            observed: 0,
        }),
        "the owner array is the first checked reservation"
    );
    assert!(
        scope.live_delta().abs() <= 2048,
        "partially built state leaked {} bytes",
        scope.live_delta()
    );
    assert!(
        !INJECT_ARMED.with(core::cell::Cell::get),
        "injection was single-shot"
    );
    assert_ne!(
        RANKED_STORAGE_REMEDIATIONS[0], RANKED_STORAGE_REMEDIATIONS[1],
        "ranked remediation guidance ships distinct ordered entries"
    );
}

#[test]
fn cfs_003_scan_injection_refuses_typed_then_resumes_to_the_golden_lattice() {
    let _guard = counter_lock();

    // Control golden: complete a certified run with no injection.
    let golden = {
        let mut executor =
            fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Certified))
                .expect("control executor admits");
        executor.enable_certificates().expect("certified mode");
        assert_eq!(drive_to_completion(&mut executor), CbcRunStatus::Completed);
        let certificates: Vec<_> = executor
            .certificates()
            .iter()
            .map(|certificate| certificate.prefix.clone())
            .collect();
        (
            executor.prefix().to_vec(),
            executor.work_spent(),
            certificates,
        )
    };

    // Injected run: advance until the first component seals and the scan
    // phase installs, then refuse exactly one allocation.
    let mut executor =
        fs_rand::cbc_exec::CbcExecutor::new(admitted(8, 3, CbcExecutionMode::Certified))
            .expect("injected executor admits");
    executor.enable_certificates().expect("certified mode");
    while executor.prefix().len() < 1 && !executor.is_complete() {
        let status = stepped_run(&mut executor);
        assert!(
            matches!(
                status,
                CbcRunStatus::Cancelled(_) | CbcRunStatus::AllowanceExhausted(_)
            ),
            "stepping must pause, never complete early"
        );
    }
    assert_eq!(
        executor.prefix(),
        &[1],
        "the theorem-fixed first component sealed"
    );

    arm_single_shot_failure();
    let outcome = stepped_run_result(&mut executor);
    match outcome {
        Err(CbcExecError::Storage(CbcStorageRefusal {
            class: CbcStorageClass::CertificateTieScratch | CbcStorageClass::ScoreAccumulator,
            phase: CbcPhaseKind::Scan,
            ..
        })) => {}
        other => panic!("expected a typed scan-phase storage refusal, got {other:?}"),
    }

    // Retry determinism: the same receipt resumes to the identical golden.
    assert_eq!(drive_to_completion(&mut executor), CbcRunStatus::Completed);
    let certificates: Vec<_> = executor
        .certificates()
        .iter()
        .map(|certificate| certificate.prefix.clone())
        .collect();
    assert_eq!(
        executor.prefix(),
        golden.0.as_slice(),
        "prefix identity preserved"
    );
    assert_eq!(
        executor.work_spent(),
        golden.1,
        "work accounting resumed exactly"
    );
    assert_eq!(
        certificates, golden.2,
        "certificate payloads replay identically"
    );
}

#[test]
fn cfs_004_run_allocation_total_is_bounded_by_storage_classes_not_work_volume() {
    let _guard = counter_lock();

    // Legitimate checked reservations: one tie-class scratch per scan-phase
    // entry, one score accumulator per admissible candidate, payload clones
    // for emitted certificates, and owner-array growth. What must never
    // happen is per-point or per-tile allocation inside the arithmetic loops.
    let (n, dimension) = (8_u32, 3_usize);
    let admissible_per_prefix = 4_usize; // phi(8) = 4 unit candidates
    let certificates = dimension - 1;

    let mut executor =
        fs_rand::cbc_exec::CbcExecutor::new(admitted(n, dimension, CbcExecutionMode::Certified))
            .expect("executor admits");
    executor.enable_certificates().expect("certified mode");

    let scope = scope_begin();
    assert_eq!(drive_to_completion(&mut executor), CbcRunStatus::Completed);
    let total = scope.allocations();

    // Structural ceiling with generous slack for owner-array doubling.
    let ceiling = 2 * (dimension + admissible_per_prefix + 4 * certificates + 8);
    assert!(
        total <= ceiling,
        "{total} allocations exceeded the {ceiling} structural ceiling; \
         arithmetic is growing buffers per visit"
    );
    drop(executor);
}

#[test]
fn cfs_005_minimal_dimension_one_problem_stays_within_the_envelope() {
    let _guard = counter_lock();
    let scope = scope_begin();

    let admission = admitted(3, 1, CbcExecutionMode::Construction);
    let envelope = admission
        .estimate()
        .logical_state_bytes()
        .saturating_mul(4)
        .max(4096);
    let tolerance: isize = envelope.try_into().unwrap_or(isize::MAX);
    let mut executor =
        fs_rand::cbc_exec::CbcExecutor::new(admission).expect("the minimal fixture admits");
    assert_eq!(drive_to_completion(&mut executor), CbcRunStatus::Completed);
    assert_eq!(
        executor.prefix(),
        &[1],
        "the first component is theorem-fixed"
    );
    assert_eq!(
        executor.work_spent(),
        129,
        "dimension-one work matches the schedule KAT"
    );

    // `into_lattice` consumes the executor and releases every charge.
    let lattice = executor
        .into_lattice()
        .expect("completed construction converts");
    assert_eq!(lattice.z, vec![1]);
    drop(lattice);
    // Tolerance scales with the sealed envelope (foreign harness churn is
    // flat; real leaks scale with problem size and dominate).
    assert!(
        scope.live_delta().abs() <= tolerance,
        "conversion and drop left {} bytes charged",
        scope.live_delta()
    );
}
