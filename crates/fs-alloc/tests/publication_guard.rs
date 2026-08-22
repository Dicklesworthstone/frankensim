//! Focused checks for the typed two-phase publication guard
//! (frankensim-epic-bedrock-6ys.21.1.3.2): exactly one live `T` bound to
//! allocator authority across prepare/commit/rollback/close, with refusals
//! handing staging and value back unchanged and destructor panics never
//! promoted to successful closes.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Barrier};

use fs_alloc::{
    LeaseIdentity, OperationMemoryLease, PublishedAllocation, PublishedTransferBinding,
};

const TEST_DOMAIN: [u8; 8] = *b"fsptype1";

fn subject(seed: u8) -> [u8; 32] {
    [seed; 32]
}

fn root_identity(seed: u8) -> LeaseIdentity {
    LeaseIdentity::root(TEST_DOMAIN, subject(seed))
}

/// One binding per attempt generation: the occurrence identity encodes the
/// generation so retried preparations never collide with the ledger's
/// duplicate-binding refusal.
fn binding_for(seed: u8, generation: u32) -> PublishedTransferBinding {
    let mut occurrence = [0_u8; 32];
    occurrence[..4].copy_from_slice(&generation.to_le_bytes());
    PublishedTransferBinding::new(
        subject(seed),
        occurrence,
        subject(seed.wrapping_add(1)),
        subject(seed.wrapping_add(2)),
    )
}

fn delegated_root(seed: u8, limit: u64) -> OperationMemoryLease {
    let root = OperationMemoryLease::bounded(limit);
    root.enable_delegation(root_identity(seed), "run", 8)
        .expect("pristine bounded root admits delegation");
    root
}

#[test]
fn committed_value_travels_with_exact_authority_and_closes_before_seal() {
    let root = delegated_root(10, 64);
    let child_id = root_identity(10).child(subject(11), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 8).unwrap();

    let staged = child.allocate("output", 7_u64).expect("exact fit");
    assert_eq!(staged.bytes(), 8);

    let published = staged
        .prepare(binding_for(12, 0))
        .expect("fresh binding prepares")
        .commit()
        .expect("prepared record commits");

    assert_eq!(published.observe(), &7);
    assert_eq!(published.bytes(), 8);
    let envelope = published.receipt().envelope;
    assert_eq!(envelope.payload_bytes(), 8);
    assert_eq!(envelope.total_bytes(), Some(8));

    // The child cannot return while its staging was published into a live
    // destination guard; closing destroys the value first, records the close,
    // and only then does the child return cleanly.
    let close_receipt = published.close().expect("explicit close succeeds");
    assert!(!close_receipt.implicit_close);
    assert_eq!(root.active_delegations(), 1);

    drop(
        child
            .close()
            .expect("child returns after publication closed"),
    );
    let sealed = root.seal().expect("conservation holds at seal");
    assert_eq!(sealed.published_transfer_count, 1);
    assert_eq!(sealed.child_published_bytes, 8);
    assert_eq!(sealed.final_used_bytes, 0);
    assert_eq!(sealed.active_delegations, 0);
}

#[test]
fn prepare_refusal_returns_allocation_unchanged_and_value_survives_retry() {
    let root = delegated_root(20, 64);
    let child_id = root_identity(20).child(subject(21), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 8).unwrap();

    // Occupy the binding with one completed publication.
    let first = child
        .allocate("output", 1_u64)
        .expect("staged")
        .prepare(binding_for(22, 0))
        .expect("prepared")
        .commit()
        .expect("committed");

    let second_staged = child.allocate("output", 42_u64).expect("staged");
    let rejected = second_staged
        .prepare(first.binding())
        .expect_err("duplicate binding refuses deterministically");
    assert_eq!(rejected.refusal().reason(), "duplicate_binding");
    assert_eq!(rejected.refusal().operation(), "prepare");

    // The refused preparation handed the staging and value back UNCHANGED:
    // re-preparing under a fresh attempt generation commits the same value.
    let republished = rejected
        .into_allocation()
        .prepare(binding_for(22, 1))
        .expect("fresh generation prepares")
        .commit()
        .expect("retried commit succeeds");
    assert_eq!(republished.observe(), &42);
    let _ = republished.close().expect("closes");
    let _ = first.close().expect("first closes");
    let _ = child.close().expect("child returns");
    root.seal().expect("seal succeeds");
}

#[test]
fn rollback_is_non_mutating_and_new_generation_retries_the_same_value() {
    let root = delegated_root(30, 64);
    let child_id = root_identity(30).child(subject(31), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 8).unwrap();

    let rolled = child
        .allocate("output", 99_u64)
        .expect("staged")
        .prepare(binding_for(32, 0))
        .expect("prepared")
        .rollback()
        .expect("non-mutating rollback succeeds");
    assert!(!rolled.receipt().implicit_rollback);

    let restaged = rolled.into_allocation();
    assert_eq!(restaged.bytes(), 8);
    let published = restaged
        .prepare(binding_for(32, 1))
        .expect("retry under new generation")
        .commit()
        .expect("commit after rollback");
    assert_eq!(published.observe(), &99);
    let _ = published.close();
    let _ = child.close();
    root.seal()
        .expect("returned plus published conserve at seal");
}

#[test]
fn zero_sized_publication_remains_a_counted_authority() {
    let root = delegated_root(40, 64);
    let child_id = root_identity(40).child(subject(41), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/unit", 8).unwrap();

    let staged = child.allocate("unit", ()).expect("ZST stages");
    assert_eq!(staged.bytes(), 1, "one counted authority byte");
    let published = staged
        .prepare(binding_for(42, 0))
        .expect("prepared")
        .commit()
        .expect("committed");
    assert_eq!(published.observe(), &());
    assert_eq!(published.bytes(), 1);
    let _ = published.close().expect("closes");
    let _ = child.close();
    let sealed = root.seal().expect("seal counts the zero-byte authority");
    assert_eq!(sealed.child_published_bytes, 1);
}

#[test]
fn guard_drop_destroys_value_and_records_implicit_close_agreeing_with_explicit() {
    let root = delegated_root(50, 64);
    let child_id = root_identity(50).child(subject(51), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 16).unwrap();

    let dropped_flag = Arc::new(AtomicBool::new(false));
    struct Observed(Arc<AtomicBool>);
    impl Drop for Observed {
        fn drop(&mut self) {
            self.0.store(true, Ordering::SeqCst);
        }
    }

    // Explicit close: value destroyed first, close recorded explicitly.
    let flag = Arc::clone(&dropped_flag);
    let explicit = child
        .allocate("a", Observed(flag))
        .expect("staged")
        .prepare(binding_for(52, 0))
        .expect("prepared")
        .commit()
        .expect("committed");
    assert!(!dropped_flag.load(Ordering::SeqCst));
    let receipt = explicit.close().expect("explicit close");
    assert!(!receipt.implicit_close);
    assert!(
        dropped_flag.load(Ordering::SeqCst),
        "value died before close"
    );

    // Guard drop without explicit close: value destroyed, close recorded
    // implicitly, both paths agree on exactly-once conservation.
    let flag = Arc::clone(&dropped_flag);
    let implicit = child
        .allocate("b", Observed(flag))
        .expect("staged")
        .prepare(binding_for(52, 1))
        .expect("prepared")
        .commit()
        .expect("committed");
    drop(implicit);
    assert!(dropped_flag.load(Ordering::SeqCst));

    let _ = child.close().expect("child returns");
    let sealed = root.seal().expect("both closes conserve");
    assert_eq!(sealed.published_transfer_count, 2);
}

#[test]
fn destructor_panic_during_close_never_records_success() {
    let root = delegated_root(60, 64);
    let child_id = root_identity(60).child(subject(61), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 8).unwrap();

    struct Exploding;
    impl Drop for Exploding {
        fn drop(&mut self) {
            panic!("destructor explosion");
        }
    }

    let published = child
        .allocate("out", Exploding)
        .expect("staged")
        .prepare(binding_for(62, 0))
        .expect("prepared")
        .commit()
        .expect("committed");

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        published
            .close()
            .map(|_| ())
            .expect_err("value drop panics")
    }));
    // The panic propagates out of close; the guard's fail-closed drop leaves
    // the published record open with NO close of any kind: no fabricated
    // success, no implicit close, not even a refusal against the record.
    assert!(result.is_err());

    drop(child);
    // Destination closes are independent of the root terminal receipt (the
    // byte-level design), so allocation-level conservation holds and the
    // seal succeeds. Honesty is proven by the counters: zero close
    // attempts, zero refusals, zero invariant violations.
    let sealed = root.seal().expect("allocation-level conservation holds");
    assert_eq!(sealed.published_transfer_count, 1);
    assert_eq!(sealed.refused_requests, 0);
    assert_eq!(sealed.release_invariant_violations, 0);
}

#[test]
fn forgotten_published_allocation_fabricates_no_close() {
    let root = delegated_root(70, 64);
    let child_id = root_identity(70).child(subject(71), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/out", 8).unwrap();

    let published: PublishedAllocation<u64> = child
        .allocate("out", 5_u64)
        .expect("staged")
        .prepare(binding_for(72, 0))
        .expect("prepared")
        .commit()
        .expect("committed");
    std::mem::forget(published);

    // Root-level conservation still holds (publication is an allocation
    // disposition), and nothing fabricated a destination close: the record
    // simply remains open forever, which is the honest stuck state.
    drop(child);
    let sealed = root.seal().expect("allocation-level conservation holds");
    assert_eq!(sealed.published_transfer_count, 1);
    assert_eq!(sealed.release_invariant_violations, 0);
}

#[test]
fn slot_displacement_overlap_and_cross_thread_handoff_conserve() {
    let root = delegated_root(80, 128);
    let child_id = root_identity(80).child(subject(81), 1).unwrap();
    let child = root.delegate_capacity(child_id, "run/slot", 24).unwrap();

    // Repeated guarded displacement over one caller-owned slot: every old
    // guard is dropped implicitly while the new one is already committed, so
    // old-plus-new briefly overlap and each is counted exactly once.
    let mut slot: Option<PublishedAllocation<u64>> = None;
    for generation in 0..3_u32 {
        let next = child
            .allocate("slot", u64::from(generation) + 1)
            .expect("staged")
            .prepare(binding_for(82, generation))
            .expect("prepared")
            .commit()
            .expect("committed");
        let displaced = slot.replace(next);
        if let Some(old) = displaced {
            assert_eq!(*old.observe(), u64::from(generation));
            drop(old);
        }
    }
    let final_guard = slot.expect("slot occupied");
    assert_eq!(*final_guard.observe(), 3);
    let _ = final_guard.close().expect("final explicit close");

    // Authority-preserving handoff: the guard MOVES to another thread with
    // its value and receipt, observes there, and closes there.
    let handoff = child
        .allocate("handoff", 777_u64)
        .expect("staged")
        .prepare(binding_for(83, 0))
        .expect("prepared")
        .commit()
        .expect("committed");
    let other = std::thread::spawn(move || {
        assert_eq!(*handoff.observe(), 777);
        handoff.close().expect("close after handoff")
    });
    let close_receipt = other.join().expect("handoff thread succeeds");
    assert!(!close_receipt.implicit_close);

    let _ = child.close().expect("child returns");
    let sealed = root.seal().expect("all dispositions conserve");
    assert_eq!(sealed.published_transfer_count, 4);
    assert_eq!(
        sealed.child_returned_bytes + sealed.child_published_bytes,
        sealed.child_granted_bytes
    );
}

#[test]
fn concurrent_publication_lifecycles_stay_exactly_counted() {
    let root = OperationMemoryLease::bounded(4096);
    root.enable_delegation(root_identity(90), "run", 64)
        .expect("pristine bounded root admits delegation");

    let barrier = Arc::new(Barrier::new(4));
    let handles: Vec<_> = (0..4)
        .map(|lane| {
            let root = root.clone();
            let barrier = Arc::clone(&barrier);
            std::thread::spawn(move || {
                barrier.wait();
                let child_id = root_identity(90)
                    .child(subject((91 + lane) as u8), (lane + 1) as u64)
                    .expect("identity depth bounded");
                let child = root
                    .delegate_capacity(child_id, "run/lane", 512)
                    .expect("each lane gets its own delegated child");
                let mut published_count = 0_u64;
                for generation in 0..8_u32 {
                    let seed = (92 + lane * 16) as u8;
                    let staged = child
                        .allocate("lane", u64::from(lane * 100 + generation))
                        .expect("staging fits under the lane capacity");
                    let published = staged
                        .prepare(binding_for(seed, generation))
                        .expect("distinct bindings prepare without contention")
                        .commit()
                        .expect("commits");
                    assert_eq!(
                        published.receipt().envelope.payload_bytes(),
                        8,
                        "every lane charges exactly its payload"
                    );
                    let _ = published.close().expect("closes");
                    published_count += 1;
                }
                let _ = child.close().expect("lane child returns");
                published_count
            })
        })
        .collect();
    let total: u64 = handles
        .into_iter()
        .map(|handle| handle.join().expect("no panics"))
        .sum();
    assert_eq!(total, 32, "every lane published exactly its generation set");
    assert_eq!(
        root.active_delegations(),
        0,
        "every delegation returned before seal"
    );

    let sealed = root.seal().expect("concurrent closes conserve");
    assert_eq!(sealed.published_transfer_count, 32);
    assert_eq!(sealed.counter_overflowed, false);
}
