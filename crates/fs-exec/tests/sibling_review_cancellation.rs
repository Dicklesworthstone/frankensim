//! Independent sibling review: asupersync cancellation contract, drilled at the
//! FrankenSim boundary (bead `frankensim-extreal-program-f85xj.13.5`).
//!
//! Charter: `docs/SIBLING_REVIEW_ASUPERSYNC.md`.
//!
//! METHOD — these drills were written **contract-first**: each one names the
//! claim it attacks, taken from `asupersync_v4_formal_semantics.md` (the
//! sibling's own stated guarantees) and `crates/fs-exec/CONTRACT.md`
//! (FrankenSim's usage assumptions), *before* reading the implementation of
//! either. That ordering is the point. Common authorship across the
//! constellation means a shared assumption is invisible to tests written by
//! someone who already knows how the code works.
//!
//! SCOPE — this is the `fs-exec` adapter boundary, which is what FrankenSim
//! actually depends on and what is FrankenSim-runnable. It is NOT a review of
//! asupersync's internal scheduler, and passing drills here do not certify the
//! sibling's own invariants: they certify that the boundary FrankenSim consumes
//! behaves as the contract says at the observable surface.
//!
//! Each drill is deliberately cheap and deterministic. None spawns a runtime,
//! so none can flake on scheduling.

use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};

fn key() -> StreamKey {
    StreamKey {
        seed: 0xD13_5,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    }
}

/// D1 — `inv.cancel.idempotence`, and `CancelGate::request`'s explicit claim
/// that "the FIRST request's timestamp is the one latency histograms measure
/// from".
///
/// Attack: request cancellation repeatedly. If a later request overwrote the
/// retained timestamp, every cancel-latency histogram in the workspace would
/// silently measure from the wrong origin — a corrupted measurement that still
/// looks like a plausible number, which is worse than an obviously missing one.
#[test]
fn d1_repeated_cancellation_requests_do_not_move_the_retained_first_timestamp() {
    let gate = CancelGate::new();
    assert!(!gate.is_requested(), "a fresh gate must be un-requested");
    assert_eq!(
        gate.requested_at_ns(),
        None,
        "an un-requested gate must not report a request timestamp"
    );

    gate.request();
    let first = gate
        .requested_at_ns()
        .expect("an ordinary gate retains a request timestamp");

    // Spin long enough that a re-stamping regression would be visible against
    // the monotonic origin, then request many more times.
    let mut spin = 0u64;
    while gate.now_ns() <= first {
        spin = spin.wrapping_add(1);
        if spin > 50_000_000 {
            break;
        }
    }
    for _ in 0..64 {
        gate.request();
    }

    assert_eq!(
        gate.requested_at_ns(),
        Some(first),
        "request() is documented idempotent: the retained timestamp must stay the FIRST one"
    );
    assert!(gate.is_requested());
}

/// D2 — `rule.cancel.checkpoint_masked`: a checkpoint must observe a
/// cancellation that was requested through the gate the context was built on.
#[test]
fn d2_checkpoint_observes_a_request_made_through_the_bound_gate() {
    let gate = CancelGate::new();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            key(),
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(
            cx.checkpoint().is_ok(),
            "an un-requested gate must not report cancellation"
        );
        gate.request();
        assert!(
            cx.checkpoint().is_err(),
            "a checkpoint must observe cancellation requested through its bound gate"
        );
    });
}

/// D3 — "Cancel protocol: request → drain → finalize is **monotone** and
/// idempotent" (proof obligation 5).
///
/// Attack: poll the checkpoint many times after cancellation. Monotone means a
/// cancelled context can never report itself runnable again. A transient Ok
/// here would let a kernel resume work after its owner cancelled it — the
/// failure mode cancellation exists to prevent.
#[test]
fn d3_cancellation_is_monotone_and_never_reverts_across_repeated_polls() {
    let gate = CancelGate::new();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            key(),
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        gate.request();
        for poll in 0..10_000 {
            assert!(
                cx.checkpoint().is_err(),
                "cancellation reverted at poll {poll}: the protocol is documented monotone"
            );
        }
        // Idempotence across the boundary: further requests change nothing.
        gate.request();
        assert!(cx.checkpoint().is_err());
    });
}

/// D4 — `fs-exec/CONTRACT.md`: "clock-free manual gates produce an empty
/// latency sample set and explicitly make no latency claim."
///
/// Attack: a clock-free gate must refuse to supply a latency observation rather
/// than fabricate a plausible one. This is the certificate-integrity form of the
/// cancellation contract: the honest answer is *no sample*, and a zero or a
/// sentinel presented as a measurement would be a fabricated one.
#[test]
fn d4_a_clock_free_gate_makes_no_latency_claim_rather_than_a_fabricated_one() {
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        gate.now_ns(),
        0,
        "a clock-free gate has no timestamp domain and must report 0, not a synthesized reading"
    );

    gate.request();
    assert!(
        gate.is_requested(),
        "a clock-free gate must still carry cancellation state"
    );

    // The contract permits the internal sentinel, but the gate must not present
    // a *derived elapsed latency*. Whatever it exposes must be the documented
    // sentinel, never a growing measurement.
    let stamp = gate.requested_at_ns();
    let mut spin = 0u64;
    while spin < 5_000_000 {
        spin = spin.wrapping_add(1);
    }
    assert_eq!(
        gate.requested_at_ns(),
        stamp,
        "a clock-free gate must not accrue elapsed time: it has no clock to measure with"
    );
    assert_eq!(
        gate.now_ns(),
        0,
        "a clock-free gate must never begin reporting a nonzero 'now'"
    );
}

/// D5 — the seeded-regression drill required by the bead's drill-quality bar.
///
/// Historically fixed sibling-boundary defect (recorded in the repository
/// README under "Caller-owned cancellation gates"): a race, pause, or
/// memory-pressure response could manufacture *private* cancellation state that
/// the owner could not observe. The fix made gates caller-supplied.
///
/// Attack: the context must observe **the caller's** gate, not a detached copy.
/// Seeded regression that this drill fails on: any change where `Cx::new`
/// snapshots or clones the gate's state instead of borrowing the caller's
/// instance. Under that regression the caller's later `request()` becomes
/// invisible to the context, `checkpoint()` keeps returning `Ok`, and the owner
/// has no way to stop work it started — exactly the defect that was fixed.
#[test]
fn d5_the_context_observes_the_callers_own_gate_not_a_detached_copy() {
    let gate = CancelGate::new();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        // Build the context FIRST, then cancel through the caller's handle.
        // A snapshotting implementation passes the pre-cancel state and fails.
        let cx = Cx::new(
            &gate,
            arena,
            key(),
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(cx.checkpoint().is_ok());

        gate.request();

        assert!(
            cx.checkpoint().is_err(),
            "REGRESSION: the context did not observe cancellation requested through the \
             caller's own gate. Cancellation state has become private to the context, so a \
             session owner can no longer stop work it started (see README, caller-owned \
             cancellation gates)."
        );
        // And the caller retains its own observability of the same state.
        assert!(
            gate.is_requested(),
            "the caller must be able to observe the state it requested"
        );
    });
}

/// D5b — NEGATIVE CONTROL for D5, and the evidence that D5 discriminates.
///
/// A suite in which every drill passes proves nothing on its own: an assertion
/// that cannot fail is not a test. This control seeds the regression D5 exists
/// to catch — a context bound to a gate that is *not* the caller's handle,
/// which is precisely the observable of "private cancellation state the owner
/// cannot observe" — and shows the boundary behaves measurably differently.
///
/// If this control ever starts reporting cancellation, the two situations have
/// become indistinguishable and D5 has silently stopped discriminating.
#[test]
fn d5b_negative_control_a_detached_gate_is_observably_different() {
    let caller_gate = CancelGate::new();
    let detached_gate = CancelGate::new();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        // The seeded regression: the context is built on something other than
        // the gate the owner will cancel through.
        let cx = Cx::new(
            &detached_gate,
            arena,
            key(),
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );

        caller_gate.request();

        assert!(
            caller_gate.is_requested(),
            "the caller did request cancellation"
        );
        assert!(
            cx.checkpoint().is_ok(),
            "control invalid: a context on a detached gate must NOT see the caller's request; \
             if it does, D5 cannot distinguish caller-owned state from private state"
        );

        // And the drill's own condition does fire once the bound gate is used,
        // confirming the assertion in D5 is the discriminating one.
        detached_gate.request();
        assert!(cx.checkpoint().is_err());
    });
}

/// D6 — two contexts sharing one caller gate must both observe a single
/// request. Cancellation "propagates down" (`inv.cancel.propagates_down`) to
/// every consumer of the owner's gate, not just the first one built.
#[test]
fn d6_one_request_is_observed_by_every_context_sharing_the_gate() {
    let gate = CancelGate::new();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let a = Cx::new(
            &gate,
            arena,
            key(),
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let b = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xD13_5,
                kernel_id: 2,
                tile: 7,
                iteration: 3,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(a.checkpoint().is_ok() && b.checkpoint().is_ok());
        gate.request();
        assert!(
            a.checkpoint().is_err() && b.checkpoint().is_err(),
            "a single owner request must reach every context built on that gate"
        );
    });
}
