//! Allocation probe for hostile quantity-literal refusal.
//!
//! This is a separate integration-test binary so its test-only global
//! allocator cannot mix measurements with the broader conformance battery.

use fs_qty::parse::{ParseBudget, ParseErrorKind, ParseResource, parse_qty_with_budget};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

struct PeakAllocation;

static TRACKING: AtomicBool = AtomicBool::new(false);
static LARGEST_REQUEST: AtomicUsize = AtomicUsize::new(0);

fn record_allocation(size: usize) {
    if TRACKING.load(Ordering::Relaxed) {
        LARGEST_REQUEST.fetch_max(size, Ordering::Relaxed);
    }
}

// SAFETY: this test-only allocator delegates every operation to `System`
// unchanged. The sole instrumentation is relaxed atomic bookkeeping, which
// neither allocates nor touches the allocation being delegated.
#[allow(unsafe_code)]
unsafe impl GlobalAlloc for PeakAllocation {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record_allocation(layout.size());
        // SAFETY: `GlobalAlloc::alloc` forwards the caller's valid layout.
        unsafe { System.alloc(layout) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record_allocation(layout.size());
        // SAFETY: `GlobalAlloc::alloc_zeroed` forwards the valid layout.
        unsafe { System.alloc_zeroed(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: pointer and layout came from the delegated system allocator.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        record_allocation(new_size);
        // SAFETY: pointer/layout/new-size preserve `GlobalAlloc::realloc`'s
        // caller-provided contract and are delegated unchanged.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: PeakAllocation = PeakAllocation;

struct TrackingWindow;

impl TrackingWindow {
    fn start() -> Self {
        LARGEST_REQUEST.store(0, Ordering::Relaxed);
        assert!(
            !TRACKING.swap(true, Ordering::Relaxed),
            "allocation tracking windows must not overlap"
        );
        Self
    }
}

impl Drop for TrackingWindow {
    fn drop(&mut self) {
        TRACKING.store(false, Ordering::Relaxed);
    }
}

#[test]
fn oversized_refusal_never_allocates_a_source_sized_buffer() {
    let input = format!("1{}", "x".repeat(1_000_000));
    let budget = ParseBudget::new(32, 1, 8, 48);

    let window = TrackingWindow::start();
    let error = parse_qty_with_budget(&input, budget).expect_err("byte admission must refuse");
    let largest_request = LARGEST_REQUEST.load(Ordering::Relaxed);
    drop(window);

    assert!(matches!(
        error.kind,
        ParseErrorKind::BudgetExceeded {
            resource: ParseResource::InputBytes,
            limit: 32,
            observed_at_least: 1_000_001,
        }
    ));
    assert_eq!(error.source_hash(), None);
    assert!(error.preview.len() <= budget.max_diagnostic_bytes());
    assert!(
        largest_request <= 1_024,
        "oversized refusal made a {largest_request}-byte allocation for a {}-byte source",
        input.len()
    );
}
