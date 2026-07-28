//! Focused G0/G4 checks for serialized affine memory delegation.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Barrier};

use fs_alloc::{
    LeaseIdentity, OperationMemoryLease, PublishedTransferBinding, PublishedTransferEnvelope,
};

const TEST_DOMAIN: [u8; 8] = *b"fstest01";

fn subject(seed: u8) -> [u8; 32] {
    [seed; 32]
}

fn root_identity(seed: u8) -> LeaseIdentity {
    LeaseIdentity::root(TEST_DOMAIN, subject(seed))
}

fn child_identity(parent: LeaseIdentity, seed: u8, component: u64) -> LeaseIdentity {
    parent
        .child(subject(seed), component)
        .expect("test identity depth is bounded")
}

fn publication_binding(seed: u8) -> PublishedTransferBinding {
    PublishedTransferBinding::new(
        subject(seed),
        subject(seed.wrapping_add(1)),
        subject(seed.wrapping_add(2)),
        subject(seed.wrapping_add(3)),
    )
}

#[test]
fn typed_identity_encoding_is_fixed_versioned_and_path_bounded() {
    let root = root_identity(200);
    let first = child_identity(root, 201, 7);
    let subject_variant = child_identity(root, 202, 7);
    let path_variant = child_identity(root, 201, 8);

    assert_eq!(LeaseIdentity::SCHEMA_VERSION, 1);
    assert_eq!(LeaseIdentity::DOMAIN_BYTES, 8);
    assert_eq!(LeaseIdentity::SUBJECT_BYTES, 32);
    assert_eq!(LeaseIdentity::ENCODED_BYTES, 235);
    assert_eq!(root.canonical_bytes().len(), LeaseIdentity::ENCODED_BYTES);
    assert_eq!(root.canonical_bytes(), root.canonical_bytes());
    assert_ne!(root.canonical_bytes(), first.canonical_bytes());
    assert_ne!(first.canonical_bytes(), subject_variant.canonical_bytes());
    assert_ne!(first.canonical_bytes(), path_variant.canonical_bytes());
    assert_eq!(first.root_subject(), root.root_subject());
    assert_eq!(first.parent_subject(), root.owner_subject());
    assert_eq!(first.path(), &[7]);
    assert!(first.to_json().contains("\"schema_version\":1"));
    assert!(!first.to_json().contains("0x"));

    let mut deepest = root;
    for component in 0..LeaseIdentity::MAX_PATH_COMPONENTS {
        deepest = deepest
            .child(
                subject(203),
                u64::try_from(component).expect("fixed depth fits u64"),
            )
            .expect("sixteen fixed path components fit");
    }
    let refusal = deepest
        .child(subject(204), 16)
        .expect_err("the seventeenth component is refused without allocation");
    assert_eq!(
        refusal.maximum_components(),
        LeaseIdentity::MAX_PATH_COMPONENTS
    );
}

#[test]
fn publication_binding_encoding_keeps_all_four_authority_axes_distinct() {
    let binding = publication_binding(220);
    let plan_variant =
        PublishedTransferBinding::new(subject(221), subject(221), subject(222), subject(223));
    let occurrence_variant =
        PublishedTransferBinding::new(subject(220), subject(222), subject(222), subject(223));
    let output_variant =
        PublishedTransferBinding::new(subject(220), subject(221), subject(223), subject(223));
    let destination_variant =
        PublishedTransferBinding::new(subject(220), subject(221), subject(222), subject(224));

    assert_eq!(PublishedTransferBinding::SCHEMA_VERSION, 1);
    assert_eq!(PublishedTransferBinding::FIELD_BYTES, 32);
    assert_eq!(PublishedTransferBinding::ENCODED_BYTES, 130);
    assert_eq!(
        binding.canonical_bytes().len(),
        PublishedTransferBinding::ENCODED_BYTES
    );
    assert_ne!(binding.canonical_bytes(), plan_variant.canonical_bytes());
    assert_ne!(
        binding.canonical_bytes(),
        occurrence_variant.canonical_bytes()
    );
    assert_ne!(binding.canonical_bytes(), output_variant.canonical_bytes());
    assert_ne!(
        binding.canonical_bytes(),
        destination_variant.canonical_bytes()
    );
    assert_eq!(binding.plan_identity(), subject(220));
    assert_eq!(binding.occurrence_identity(), subject(221));
    assert_eq!(binding.output_identity(), subject(222));
    assert_eq!(binding.destination_identity(), subject(223));
    let json = binding.to_json();
    assert!(json.contains("\"schema\":\"fs-alloc-published-transfer-binding-v1\""));
    assert!(json.contains("\"plan_identity\":\""));
    assert!(json.contains("\"destination_identity\":\""));
    assert!(!json.contains("0x"));
}

#[test]
fn real_allocation_exact_under_and_over_charge_one_root_envelope() {
    let root = OperationMemoryLease::bounded(64);
    let root_id = root_identity(1);
    let child_id = child_identity(root_id, 2, 0);
    root.enable_delegation(root_id, "run-exact", 2)
        .expect("metadata is pre-admitted on a pristine root");
    let child = root
        .delegate_capacity(child_id, "run-exact/output", 64)
        .expect("exact child envelope fits");

    let charge = child
        .reserve("output/payload", 63)
        .expect("one-under allocation fits");
    let mut allocation = Vec::<u8>::new();
    allocation
        .try_reserve_exact(63)
        .expect("real backing allocation succeeds");
    assert_eq!(child.used_bytes(), 63);
    assert_eq!(root.receipt().used_bytes, 64);
    assert_eq!(
        root.receipt().requested_bytes,
        64,
        "child allocations consume the existing envelope"
    );

    let one_over = child
        .reserve("output/one-over", 2)
        .expect_err("63 + 2 exceeds a 64-byte child");
    assert_eq!(one_over.reason(), "capacity");
    assert_eq!(one_over.root_id(), "run-exact");
    assert_eq!(one_over.logical_path(), "run-exact/output");
    assert_eq!(one_over.site(), "output/one-over");
    assert_eq!(one_over.requested_bytes(), 2);
    assert_eq!(one_over.used_bytes(), 63);
    assert_eq!(one_over.limit_bytes(), 64);
    assert_eq!(one_over.available_bytes(), Some(1));

    let exact = child
        .reserve("output/final-byte", 1)
        .expect("the final byte exactly fills the child");
    assert_eq!(child.used_bytes(), 64);
    assert_eq!(child.peak_used_bytes(), 64);
    drop(exact);
    drop(allocation);
    drop(charge);

    let child_receipt = child.close().expect("drained child returns exactly once");
    child_receipt
        .verify_for(root_id, None, child_id)
        .expect("child receipt verifies");
    assert_eq!(child_receipt.allocation_granted_bytes, 64);
    assert_eq!(child_receipt.allocation_returned_bytes, 64);
    assert_eq!(child_receipt.refused_requests, 1);

    let root_receipt = root.seal().expect("quiescent root closes");
    root_receipt
        .verify_for(root_id)
        .expect("root receipt verifies");
    assert_eq!(root_receipt.delegated_bytes, 64);
    assert_eq!(root_receipt.returned_delegated_bytes, 64);
    assert_eq!(root_receipt.child_granted_bytes, 64);
    assert_eq!(root_receipt.child_returned_bytes, 64);
    assert_eq!(root_receipt.receipt().used_bytes, 0);
}

#[test]
fn zero_capacity_is_affine_while_invalid_duplicate_and_metadata_escape_refuse() {
    let root = OperationMemoryLease::bounded(4);
    let root_id = root_identity(3);
    let zero_id = child_identity(root_id, 4, 0);
    let first_id = child_identity(root_id, 5, 1);
    let second_id = child_identity(root_id, 6, 2);
    let third_id = child_identity(root_id, 7, 3);
    let wrong_root = root_identity(205);
    let wrong_authority_id = child_identity(wrong_root, 5, 1);
    let path_collision_id = child_identity(root_id, 210, 1);
    root.enable_delegation(root_id, "run-meta", 3)
        .expect("three records are pre-admitted");

    let zero = root
        .delegate_capacity(zero_id, "run-meta/zero", 0)
        .expect("zero capacity is still a real affine control authority");
    assert_eq!(zero.capacity_bytes(), 0);
    assert_eq!(zero.used_bytes(), 0);
    assert_eq!(root.active_delegations(), 1);
    zero.close().expect("zero authority returns exactly once");

    let invalid = root
        .delegate_capacity(first_id, "other/not-a-child", 1)
        .expect_err("the path must include its exact root");
    assert_eq!(invalid.reason(), "invalid_logical_path");

    let invalid_identity = root
        .delegate_capacity(wrong_authority_id, "run-meta/wrong-root", 1)
        .expect_err("a diagnostic label cannot override typed parent identity");
    assert_eq!(invalid_identity.reason(), "invalid_identity_relationship");
    assert_eq!(invalid_identity.root_identity(), Some(root_id));
    assert_eq!(invalid_identity.identity(), wrong_authority_id);

    let first = root
        .delegate_capacity(first_id, "run-meta/first", 1)
        .expect("first identity fits");
    first.close().expect("first identity returns");
    let duplicate = root
        .delegate_capacity(first_id, "run-meta/first", 1)
        .expect_err("returned identities cannot be reused");
    assert_eq!(duplicate.reason(), "duplicate_identity");
    let path_collision = root
        .delegate_capacity(path_collision_id, "run-meta/path-collision", 1)
        .expect_err("a sibling ordinal cannot be rebound to another subject");
    assert_eq!(path_collision.reason(), "duplicate_path");

    let second = root
        .delegate_capacity(second_id, "run-meta/second", 1)
        .expect("second retained record fits");
    second.close().expect("second identity returns");
    let exhausted = root
        .delegate_capacity(third_id, "run-meta/third", 1)
        .expect_err("record history is bounded and retained");
    assert_eq!(exhausted.reason(), "metadata_exhausted");

    let receipt = root.seal().expect("all accepted transfers returned");
    assert_eq!(receipt.delegation_count, 3);
    assert_eq!(receipt.metadata_limit, 3);
    assert_eq!(receipt.refused_requests, 5);
    assert_eq!(receipt.refused_bytes, 5);
    receipt.verify_for(root_id).expect("receipt verifies");
}

#[test]
fn live_zero_capacity_child_blocks_seal_until_exact_return() {
    let root = OperationMemoryLease::bounded(0);
    let root_id = root_identity(211);
    let zero_id = child_identity(root_id, 212, 0);
    root.enable_delegation(root_id, "run-zero-live", 1)
        .expect("metadata is pre-admitted despite a zero byte cap");
    let zero = root
        .delegate_capacity(zero_id, "run-zero-live/control", 0)
        .expect("zero control authority is admitted");
    let refusal = root
        .seal()
        .expect_err("zero used bytes cannot hide a live affine authority");
    assert_eq!(refusal.reason(), "live_capacity");
    assert_eq!(refusal.used_bytes(), 0);
    assert_eq!(refusal.active_delegations(), 1);
    zero.close().expect("zero control authority returns");
    let terminal = root.seal().expect("drained zero-capacity root closes");
    terminal.verify_for(root_id).expect("root verifies");
    assert_eq!(terminal.delegated_bytes, 0);
    assert_eq!(terminal.returned_delegated_bytes, 0);
}

#[test]
fn delegation_configuration_is_fallible_bounded_and_pristine_only() {
    let unbounded = OperationMemoryLease::unbounded();
    let unbounded_id = root_identity(8);
    let refusal = unbounded
        .enable_delegation(unbounded_id, "unbounded", 1)
        .expect_err("compatibility accounting cannot create bounded evidence");
    assert_eq!(refusal.reason(), "unbounded_root");

    let zero_metadata = OperationMemoryLease::bounded(1);
    let zero_metadata_id = root_identity(9);
    let refusal = zero_metadata
        .enable_delegation(zero_metadata_id, "zero-metadata", 0)
        .expect_err("zero metadata capacity is refused");
    assert_eq!(refusal.reason(), "metadata_limit");

    let excessive_metadata = OperationMemoryLease::bounded(1);
    let excessive_id = root_identity(10);
    let refusal = excessive_metadata
        .enable_delegation(excessive_id, "excessive", usize::MAX)
        .expect_err("an unbounded metadata request is refused before allocation");
    assert_eq!(refusal.reason(), "metadata_limit");

    let used = OperationMemoryLease::bounded(8);
    let used_id = root_identity(11);
    let charge = used.reserve("legacy-use", 1).expect("root reserve fits");
    drop(charge);
    let refusal = used
        .enable_delegation(used_id, "too-late", 1)
        .expect_err("authority identity is configured before any use");
    assert_eq!(refusal.reason(), "root_not_pristine");

    let configured = OperationMemoryLease::bounded(8);
    let configured_id = root_identity(12);
    configured
        .enable_delegation(configured_id, "configured", 1)
        .expect("first configuration succeeds");
    let refusal = configured
        .enable_delegation(root_identity(13), "replacement", 1)
        .expect_err("root identity is immutable");
    assert_eq!(refusal.reason(), "already_configured");

    let non_root = OperationMemoryLease::bounded(8);
    let root_id = root_identity(14);
    let child_id = child_identity(root_id, 15, 0);
    let refusal = non_root
        .enable_delegation(child_id, "not-a-root", 1)
        .expect_err("a descendant cannot masquerade as a root");
    assert_eq!(refusal.reason(), "invalid_root_identity");
}

#[test]
fn nested_and_sibling_transfers_conserve_each_parent_without_double_charge() {
    let root = OperationMemoryLease::bounded(100);
    let root_id = root_identity(16);
    let left_id = child_identity(root_id, 17, 0);
    let right_id = child_identity(root_id, 18, 1);
    let overflow_id = child_identity(root_id, 19, 2);
    let scratch_id = child_identity(left_id, 20, 0);
    let too_large_id = child_identity(left_id, 21, 1);
    root.enable_delegation(root_id, "run-tree", 4)
        .expect("four child records are pre-admitted");
    let left = root
        .delegate_capacity(left_id, "run-tree/left", 60)
        .expect("left fits");
    let right = root
        .delegate_capacity(right_id, "run-tree/right", 40)
        .expect("right exactly fills root");

    let one_over = root
        .delegate_capacity(overflow_id, "run-tree/overflow", 1)
        .expect_err("root is exactly occupied");
    assert_eq!(one_over.reason(), "capacity");
    assert_eq!(one_over.used_bytes(), 100);
    assert_eq!(one_over.limit_bytes(), Some(100));
    assert_eq!(one_over.available_bytes(), Some(0));

    let scratch = left
        .delegate_capacity(scratch_id, "run-tree/left/scratch", 24)
        .expect("nested child fits its exact parent");
    let nested_over = left
        .delegate_capacity(too_large_id, "run-tree/left/too-large", 37)
        .expect_err("24 + 37 exceeds the left envelope");
    assert_eq!(nested_over.reason(), "capacity");
    assert_eq!(nested_over.parent_path(), Some("run-tree/left"));
    assert_eq!(nested_over.parent_identity(), Some(left_id));
    assert_eq!(nested_over.available_bytes(), Some(36));

    let scratch_charge = scratch
        .reserve("scratch/payload", 24)
        .expect("nested real allocation fits");
    assert_eq!(scratch.used_bytes(), 24);
    assert_eq!(left.used_bytes(), 24);
    assert_eq!(
        root.receipt().used_bytes,
        100,
        "nested ownership stays inside the top-level envelope"
    );
    drop(scratch_charge);
    let scratch_receipt = scratch.close().expect("scratch returns");
    scratch_receipt
        .verify_for(root_id, Some(left_id), scratch_id)
        .expect("nested receipt verifies");
    let same_path_wrong_parent = child_identity(root_id, 209, 0);
    let mut forged_parent = scratch_receipt.clone();
    forged_parent.parent_identity = Some(same_path_wrong_parent);
    forged_parent.receipt_root = forged_parent.recompute_root();
    assert_eq!(
        forged_parent
            .verify_for(root_id, Some(same_path_wrong_parent), scratch_id)
            .expect_err("same-path parent subject substitution is rejected")
            .reason(),
        "identity_relationship"
    );
    assert_eq!(left.used_bytes(), 0);

    let left_receipt = left.close().expect("left returns");
    let right_receipt = right.close().expect("right returns");
    left_receipt
        .verify_for(root_id, None, left_id)
        .expect("left verifies");
    assert_eq!(left_receipt.refused_requests, 1);
    assert_eq!(left_receipt.refused_bytes, 37);
    right_receipt
        .verify_for(root_id, None, right_id)
        .expect("right verifies");
    let root_receipt = root.seal().expect("tree is conserved");
    root_receipt.verify_for(root_id).expect("root verifies");
    assert_eq!(root_receipt.delegated_bytes, 100);
    assert_eq!(root_receipt.returned_delegated_bytes, 100);
}

#[test]
fn deterministic_caller_paths_keep_concurrent_siblings_distinct() {
    let root = OperationMemoryLease::bounded(96);
    let root_id = root_identity(22);
    let left_id = child_identity(root_id, 23, 0);
    let middle_id = child_identity(root_id, 24, 1);
    let right_id = child_identity(root_id, 25, 2);
    root.enable_delegation(root_id, "run-parallel", 3)
        .expect("three records are pre-admitted");
    let left = root
        .delegate_capacity(left_id, "run-parallel/a", 32)
        .expect("a fits");
    let middle = root
        .delegate_capacity(middle_id, "run-parallel/b", 32)
        .expect("b fits");
    let right = root
        .delegate_capacity(right_id, "run-parallel/c", 32)
        .expect("c fits");
    let barrier = Arc::new(Barrier::new(3));

    let (a, b, c) = std::thread::scope(|scope| {
        let a_barrier = Arc::clone(&barrier);
        let a = scope.spawn(move || {
            let charge = left.reserve("parallel/a", 32).expect("a charge fits");
            a_barrier.wait();
            drop(charge);
            left.close().expect("a returns")
        });
        let b_barrier = Arc::clone(&barrier);
        let b = scope.spawn(move || {
            let charge = middle.reserve("parallel/b", 32).expect("b charge fits");
            b_barrier.wait();
            drop(charge);
            middle.close().expect("b returns")
        });
        let c_barrier = Arc::clone(&barrier);
        let c = scope.spawn(move || {
            let charge = right.reserve("parallel/c", 32).expect("c charge fits");
            c_barrier.wait();
            drop(charge);
            right.close().expect("c returns")
        });
        (
            a.join().expect("a thread"),
            b.join().expect("b thread"),
            c.join().expect("c thread"),
        )
    });

    assert_eq!(a.logical_path, "run-parallel/a");
    assert_eq!(b.logical_path, "run-parallel/b");
    assert_eq!(c.logical_path, "run-parallel/c");
    assert_ne!(a.receipt_root, b.receipt_root);
    assert_ne!(b.receipt_root, c.receipt_root);
    let terminal = root.seal().expect("all concurrent siblings returned");
    assert_eq!(terminal.delegation_count, 3);
    terminal
        .verify_for(root_id)
        .expect("parallel receipt verifies");
}

#[test]
fn seal_vs_root_reserve_has_one_serialized_winner() {
    for iteration in 0..64 {
        let root = OperationMemoryLease::bounded(8);
        let root_id = root_identity(26);
        root.enable_delegation(root_id, "run-reserve-race", 1)
            .expect("root configured");
        let start = Arc::new(Barrier::new(2));
        let seal_done = Arc::new(AtomicBool::new(false));
        let reserve_succeeded = Arc::new(AtomicBool::new(false));
        let seal_succeeded = Arc::new(AtomicBool::new(false));

        std::thread::scope(|scope| {
            let reserve_root = root.clone();
            let start_reserve = Arc::clone(&start);
            let seal_done_reserve = Arc::clone(&seal_done);
            let reserve_succeeded_out = Arc::clone(&reserve_succeeded);
            scope.spawn(move || {
                start_reserve.wait();
                match reserve_root.reserve("race/root", 8) {
                    Ok(charge) => {
                        reserve_succeeded_out.store(true, Ordering::Release);
                        while !seal_done_reserve.load(Ordering::Acquire) {
                            std::thread::yield_now();
                        }
                        drop(charge);
                    }
                    Err(refusal) => assert_eq!(refusal.reason(), "sealed"),
                }
            });

            let seal_root = root.clone();
            let start_seal = Arc::clone(&start);
            let seal_done_out = Arc::clone(&seal_done);
            let seal_succeeded_out = Arc::clone(&seal_succeeded);
            scope.spawn(move || {
                start_seal.wait();
                let result = seal_root.seal();
                seal_succeeded_out.store(result.is_ok(), Ordering::Release);
                if let Err(refusal) = result {
                    assert_eq!(refusal.reason(), "live_capacity");
                }
                seal_done_out.store(true, Ordering::Release);
            });
        });

        assert_ne!(
            reserve_succeeded.load(Ordering::Acquire),
            seal_succeeded.load(Ordering::Acquire),
            "iteration {iteration} must have exactly one serialized winner"
        );
        root.seal()
            .expect("drain followed by replay/final close succeeds")
            .verify_for(root_id)
            .expect("receipt verifies");
    }
}

#[test]
fn seal_vs_delegate_has_one_serialized_winner() {
    for iteration in 0..64 {
        let root = OperationMemoryLease::bounded(8);
        let root_id = root_identity(27);
        let child_id = child_identity(root_id, 28, 0);
        root.enable_delegation(root_id, "run-delegate-race", 1)
            .expect("root configured");
        let start = Arc::new(Barrier::new(2));
        let seal_done = Arc::new(AtomicBool::new(false));
        let delegate_succeeded = Arc::new(AtomicBool::new(false));
        let seal_succeeded = Arc::new(AtomicBool::new(false));

        std::thread::scope(|scope| {
            let delegate_root = root.clone();
            let start_delegate = Arc::clone(&start);
            let seal_done_delegate = Arc::clone(&seal_done);
            let delegate_succeeded_out = Arc::clone(&delegate_succeeded);
            scope.spawn(move || {
                start_delegate.wait();
                match delegate_root.delegate_capacity(child_id, "run-delegate-race/child", 8) {
                    Ok(child) => {
                        delegate_succeeded_out.store(true, Ordering::Release);
                        while !seal_done_delegate.load(Ordering::Acquire) {
                            std::thread::yield_now();
                        }
                        drop(child);
                    }
                    Err(refusal) => assert_eq!(refusal.reason(), "root_sealed"),
                }
            });

            let seal_root = root.clone();
            let start_seal = Arc::clone(&start);
            let seal_done_out = Arc::clone(&seal_done);
            let seal_succeeded_out = Arc::clone(&seal_succeeded);
            scope.spawn(move || {
                start_seal.wait();
                let result = seal_root.seal();
                seal_succeeded_out.store(result.is_ok(), Ordering::Release);
                if let Err(refusal) = result {
                    assert_eq!(refusal.reason(), "live_capacity");
                }
                seal_done_out.store(true, Ordering::Release);
            });
        });

        assert_ne!(
            delegate_succeeded.load(Ordering::Acquire),
            seal_succeeded.load(Ordering::Acquire),
            "iteration {iteration} must have exactly one serialized winner"
        );
        root.seal()
            .expect("returned winner permits final close")
            .verify_for(root_id)
            .expect("receipt verifies");
    }
}

#[test]
fn root_seal_freezes_child_admission_and_terminal_replay() {
    let root = OperationMemoryLease::bounded(32);
    let root_id = root_identity(29);
    let child_id = child_identity(root_id, 30, 0);
    let late_id = child_identity(root_id, 31, 1);
    root.enable_delegation(root_id, "run-seal", 2)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-seal/publication", 32)
        .expect("child fits");
    let live = child
        .reserve("publication/live", 16)
        .expect("existing work fits");

    let refusal = root.seal().expect_err("live ownership blocks close");
    assert_eq!(refusal.reason(), "live_capacity");
    assert!(root.is_sealed());
    assert_eq!(refusal.root_identity(), Some(root_id));
    assert_eq!(refusal.root_id(), Some("run-seal"));
    assert_eq!(refusal.used_bytes(), 32);
    assert_eq!(refusal.active_delegations(), 1);
    assert_eq!(refusal.release_invariant_violations(), 0);

    let child_refusal = child
        .reserve("publication/late", 1)
        .expect_err("the root cut freezes child admission too");
    assert_eq!(child_refusal.reason(), "root_sealed");
    drop(live);
    child
        .close()
        .expect("existing charge drains and child returns");

    let first = root.seal().expect("drained root closes");
    first.verify_for(root_id).expect("first verifies");
    assert_eq!(first.refused_requests, 2);
    assert_eq!(
        first.refused_bytes, 1,
        "a close observation is not a refused byte request"
    );
    let frozen_json = first.to_json();
    let root_refusal = root
        .reserve("late/root", 1)
        .expect_err("late root admission remains closed");
    assert_eq!(root_refusal.reason(), "sealed");
    let late_child = root
        .delegate_capacity(late_id, "run-seal/late", 1)
        .expect_err("late transfer remains closed");
    assert_eq!(late_child.reason(), "root_sealed");
    let replay = root.seal().expect("terminal close replays");
    assert_eq!(frozen_json, replay.to_json());
    assert_eq!(first.receipt_root, replay.receipt_root);
}

#[test]
fn replacement_growth_counts_old_plus_new_at_the_overlap() {
    let root = OperationMemoryLease::bounded(24);
    let root_id = root_identity(32);
    let child_id = child_identity(root_id, 33, 0);
    root.enable_delegation(root_id, "run-growth", 1)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-growth/vector", 24)
        .expect("child fits");
    let old = child.reserve("growth/old", 8).expect("old buffer fits");
    let replacement = child
        .reserve("growth/new", 16)
        .expect("old plus replacement exactly fits");
    assert_eq!(child.used_bytes(), 24);
    assert_eq!(child.peak_used_bytes(), 24);
    let one_over = child
        .reserve("growth/one-over", 1)
        .expect_err("replacement overlap cannot hide one byte");
    assert_eq!(one_over.reason(), "capacity");

    drop(old);
    assert_eq!(child.used_bytes(), 16);
    drop(replacement);
    let receipt = child.close().expect("replacement ownership returned");
    assert_eq!(receipt.allocation_granted_bytes, 24);
    assert_eq!(receipt.allocation_returned_bytes, 24);
    assert_eq!(receipt.peak_used_bytes, 24);
    root.seal()
        .expect("root closes")
        .verify_for(root_id)
        .expect("root verifies");
}

#[test]
fn empty_charge_and_implicit_child_drop_return_exactly_once() {
    let root = OperationMemoryLease::bounded(1);
    let root_id = root_identity(206);
    let child_id = child_identity(root_id, 207, 0);
    root.enable_delegation(root_id, "run-empty", 1)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-empty/child", 1)
        .expect("child fits");
    let empty = child
        .reserve("empty/allocation", 0)
        .expect("empty allocations retain typed ownership without capacity");
    assert_eq!(empty.bytes(), 0);
    assert_eq!(child.used_bytes(), 0);
    drop(empty);
    drop(child);

    assert_eq!(root.active_delegations(), 0);
    assert_eq!(root.receipt().used_bytes, 0);
    let receipt = root.seal().expect("implicit affine return is conserved");
    receipt.verify_for(root_id).expect("root verifies");
    assert_eq!(receipt.delegated_bytes, 1);
    assert_eq!(receipt.returned_delegated_bytes, 1);
}

#[test]
fn controllable_real_allocation_failure_returns_charge() {
    let root = OperationMemoryLease::bounded(u64::MAX);
    let root_id = root_identity(34);
    let child_id = child_identity(root_id, 35, 0);
    root.enable_delegation(root_id, "run-allocation-failure", 1)
        .expect("single metadata record pre-admitted");
    let child = root
        .delegate_capacity(child_id, "run-allocation-failure/vector", u64::MAX)
        .expect("maximum logical envelope does not allocate payload");
    let charge = child
        .reserve("vector/impossible", u64::MAX)
        .expect("logical charge is representable");
    let mut allocation = Vec::<u8>::new();
    assert!(
        allocation.try_reserve_exact(usize::MAX).is_err(),
        "the allocator refusal is controlled and reported by Vec"
    );
    drop(charge);
    assert_eq!(child.used_bytes(), 0);
    child
        .close()
        .expect("failed allocation handed its authority back");
    root.seal()
        .expect("root closes after allocation refusal")
        .verify_for(root_id)
        .expect("receipt verifies");
}

#[test]
fn cumulative_reuse_exceeds_u64_without_exceeding_live_capacity() {
    let root = OperationMemoryLease::bounded(u64::MAX);
    let first = root
        .reserve("cumulative/first", u64::MAX)
        .expect("maximum live charge fits");
    drop(first);
    let second = root
        .reserve("cumulative/second", 1)
        .expect("reused capacity admits one more byte");
    drop(second);
    assert_eq!(root.receipt().used_bytes, 0);
    assert_eq!(root.receipt().peak_bytes, u64::MAX);
    assert_eq!(
        root.receipt().requested_bytes,
        u128::from(u64::MAX) + 1,
        "the cumulative counter remains exact beyond u64"
    );
}

#[test]
fn cancellation_fault_and_unwind_return_every_authority() {
    let root = OperationMemoryLease::bounded(80);
    let root_id = root_identity(36);
    let child_id = child_identity(root_id, 37, 0);
    root.enable_delegation(root_id, "run-unwind", 1)
        .expect("root configured");
    let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let child = root
            .delegate_capacity(child_id, "run-unwind/staging", 80)
            .expect("child fits");
        let _charge = child.reserve("staging/payload", 64).expect("charge fits");
        let mut allocation = Vec::<u8>::new();
        allocation
            .try_reserve_exact(64)
            .expect("real allocation succeeds");
        panic!("exercise cancellation/fault unwind handoff");
    }));
    assert!(unwind.is_err());
    assert_eq!(root.receipt().used_bytes, 0);
    assert_eq!(root.active_delegations(), 0);
    assert_eq!(root.receipt().release_invariant_violations, 0);
    root.seal()
        .expect("unwind returned charge then owner")
        .verify_for(root_id)
        .expect("receipt verifies");
}

#[test]
fn prepared_output_can_publish_after_seal_cut_and_destination_closes_separately() {
    let root = OperationMemoryLease::bounded(24);
    let root_id = root_identity(70);
    let child_id = child_identity(root_id, 71, 0);
    let binding = publication_binding(72);
    root.enable_delegation(root_id, "run-publish", 2)
        .expect("delegation and publication metadata are pre-admitted");
    let child = root
        .delegate_capacity(child_id, "run-publish/output", 24)
        .expect("staging envelope fits");
    let prepared = child
        .reserve("output/payload", 24)
        .expect("staging bytes fit exactly")
        .prepare_published_transfer(binding)
        .expect("publication is prepared before the cut");
    assert_eq!(prepared.bytes(), 24);
    assert_eq!(prepared.binding(), binding);
    assert_eq!(child.used_bytes(), 24);

    let first_seal = root
        .seal()
        .expect_err("prepared staging remains live at the first close observation");
    assert_eq!(first_seal.reason(), "live_capacity");
    assert!(root.is_sealed());

    let published = prepared
        .publish()
        .expect("pre-admitted preparation may drain by publishing after the cut");
    let published_receipt = published.receipt().clone();
    published_receipt
        .verify_for(
            root_id,
            None,
            child_id,
            binding,
            PublishedTransferEnvelope::payload_only(24),
        )
        .expect("publication receipt binds every authority axis");
    assert_eq!(published.bytes(), 24);
    assert_eq!(published.binding(), binding);
    assert_eq!(child.used_bytes(), 0);

    let child_receipt = child
        .close()
        .expect("published staging no longer blocks child return");
    child_receipt
        .verify_for(root_id, None, child_id)
        .expect("child receipt verifies published conservation");
    assert_eq!(child_receipt.allocation_granted_bytes, 24);
    assert_eq!(child_receipt.allocation_returned_bytes, 0);
    assert_eq!(child_receipt.allocation_published_bytes, 24);
    assert_eq!(child_receipt.publication_record_count, 1);
    assert_eq!(child_receipt.published_transfer_count, 1);
    assert_eq!(child_receipt.rolled_back_transfer_count, 0);

    let terminal = root
        .seal()
        .expect("published output permits root close without destination drop");
    terminal
        .verify_for(root_id)
        .expect("root receipt verifies returned plus published conservation");
    assert_eq!(terminal.child_granted_bytes, 24);
    assert_eq!(terminal.child_returned_bytes, 0);
    assert_eq!(terminal.child_published_bytes, 24);
    assert_eq!(terminal.publication_record_count, 1);
    assert_eq!(terminal.published_transfer_count, 1);
    assert_eq!(terminal.rolled_back_transfer_count, 0);
    let frozen_terminal = terminal.to_json();

    let destination_close = published
        .close()
        .expect("destination ownership closes exactly once after root close");
    destination_close
        .verify_for(&published_receipt)
        .expect("destination close binds the publication receipt");
    assert!(!destination_close.implicit_close);
    assert_eq!(
        root.seal().expect("terminal receipt replays").to_json(),
        frozen_terminal,
        "later destination close cannot rewrite terminal root evidence"
    );
}

#[test]
fn seal_cut_refuses_late_preparation_but_returns_the_consumed_staging_charge() {
    let root = OperationMemoryLease::bounded(8);
    let root_id = root_identity(76);
    let child_id = child_identity(root_id, 77, 0);
    root.enable_delegation(root_id, "run-late-prepare", 1)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-late-prepare/output", 8)
        .expect("child fits");
    let charge = child
        .reserve("output/already-staged", 8)
        .expect("staging precedes the seal");
    root.seal()
        .expect_err("the live staging allocation blocks terminal close");

    let refusal = charge
        .prepare_published_transfer(publication_binding(78))
        .expect_err("preparation is new admission and cannot cross the seal cut");
    assert_eq!(refusal.operation(), "prepare");
    assert_eq!(refusal.reason(), "root_sealed");
    assert_eq!(refusal.bytes(), 8);
    assert_eq!(
        child.used_bytes(),
        0,
        "consuming the refused charge returns staging exactly once"
    );
    child.close().expect("child drains after late refusal");
    root.seal()
        .expect("root closes after refusal cleanup")
        .verify_for(root_id)
        .expect("root verifies");
}

#[test]
fn prepared_publication_explicit_and_drop_rollback_return_staging_exactly_once() {
    let root = OperationMemoryLease::bounded(16);
    let root_id = root_identity(80);
    let child_id = child_identity(root_id, 81, 0);
    let explicit_binding = publication_binding(82);
    let implicit_binding = publication_binding(86);
    let unwind_binding = publication_binding(90);
    root.enable_delegation(root_id, "run-rollback", 3)
        .expect("three publication records are pre-admitted");
    let child = root
        .delegate_capacity(child_id, "run-rollback/staging", 16)
        .expect("child fits");

    let explicit = child
        .reserve("staging/explicit", 7)
        .expect("first staging charge fits")
        .prepare_published_transfer(explicit_binding)
        .expect("first preparation fits");
    let rollback = explicit.rollback().expect("explicit rollback succeeds");
    rollback
        .verify_for(root_id, None, child_id, explicit_binding)
        .expect("rollback receipt verifies");
    assert!(!rollback.implicit_rollback);
    assert_eq!(child.used_bytes(), 0);

    let implicit = child
        .reserve("staging/drop", 9)
        .expect("second staging charge fits")
        .prepare_published_transfer(implicit_binding)
        .expect("second preparation fits");
    drop(implicit);
    assert_eq!(child.used_bytes(), 0);

    let unwind = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _prepared = child
            .reserve("staging/unwind", 8)
            .expect("unwind staging fits")
            .prepare_published_transfer(unwind_binding)
            .expect("unwind preparation fits");
        panic!("exercise prepared-transfer rollback during unwind");
    }));
    assert!(unwind.is_err());
    assert_eq!(
        child.used_bytes(),
        0,
        "unwind rolls prepared staging back exactly once"
    );

    let child_receipt = child.close().expect("all rollback paths drained");
    child_receipt
        .verify_for(root_id, None, child_id)
        .expect("returned-only child receipt verifies");
    assert_eq!(child_receipt.allocation_granted_bytes, 24);
    assert_eq!(child_receipt.allocation_returned_bytes, 24);
    assert_eq!(child_receipt.allocation_published_bytes, 0);
    assert_eq!(child_receipt.publication_record_count, 3);
    assert_eq!(child_receipt.published_transfer_count, 0);
    assert_eq!(child_receipt.rolled_back_transfer_count, 3);
    let terminal = root.seal().expect("rollback-only root closes");
    terminal.verify_for(root_id).expect("root verifies");
    assert_eq!(terminal.child_granted_bytes, 24);
    assert_eq!(terminal.child_returned_bytes, 24);
    assert_eq!(terminal.child_published_bytes, 0);
}

#[test]
fn publication_refuses_zero_duplicate_and_exhausted_metadata_without_leaking_bytes() {
    let zero_root = OperationMemoryLease::bounded(1);
    let zero_root_id = root_identity(90);
    let zero_child_id = child_identity(zero_root_id, 91, 0);
    zero_root
        .enable_delegation(zero_root_id, "run-publish-zero", 1)
        .expect("root configured");
    let zero_child = zero_root
        .delegate_capacity(zero_child_id, "run-publish-zero/output", 1)
        .expect("child fits");
    let zero = zero_child
        .reserve("output/empty", 0)
        .expect("ordinary zero-byte charges remain representable")
        .prepare_published_transfer(publication_binding(92))
        .expect_err("publication requires a nonzero ownership transfer");
    assert_eq!(zero.operation(), "prepare");
    assert_eq!(zero.reason(), "zero_bytes");
    assert_eq!(zero.bytes(), 0);
    assert_eq!(zero.child_identity(), zero_child_id);
    assert_eq!(zero_child.used_bytes(), 0);
    zero_child.close().expect("zero refusal leaked no charge");
    zero_root.seal().expect("zero-refusal root closes");

    let root = OperationMemoryLease::bounded(3);
    let root_id = root_identity(96);
    let child_id = child_identity(root_id, 97, 0);
    let first_binding = publication_binding(98);
    root.enable_delegation(root_id, "run-publish-refusal", 1)
        .expect("one publication record is pre-admitted");
    let child = root
        .delegate_capacity(child_id, "run-publish-refusal/output", 3)
        .expect("child fits");
    let published = child
        .reserve("output/first", 1)
        .expect("first staging byte fits")
        .prepare_published_transfer(first_binding)
        .expect("first binding is retained")
        .publish()
        .expect("first binding publishes");

    let duplicate = child
        .reserve("output/duplicate", 1)
        .expect("duplicate probe has an owned staging byte")
        .prepare_published_transfer(first_binding)
        .expect_err("a retained binding can never be reused");
    assert_eq!(duplicate.reason(), "duplicate_binding");
    assert_eq!(duplicate.binding(), first_binding);
    assert_eq!(child.used_bytes(), 0, "refused charge returned on consume");

    let exhausted_binding = publication_binding(102);
    let exhausted = child
        .reserve("output/exhausted", 1)
        .expect("metadata probe has an owned staging byte")
        .prepare_published_transfer(exhausted_binding)
        .expect_err("publication history cannot exceed its pre-admitted bound");
    assert_eq!(exhausted.reason(), "metadata_exhausted");
    assert_eq!(exhausted.binding(), exhausted_binding);
    assert_eq!(child.used_bytes(), 0);

    let child_receipt = child.close().expect("refusals left no staging live");
    assert_eq!(child_receipt.allocation_granted_bytes, 3);
    assert_eq!(child_receipt.allocation_returned_bytes, 2);
    assert_eq!(child_receipt.allocation_published_bytes, 1);
    let terminal = root.seal().expect("mixed returned/published root closes");
    terminal.verify_for(root_id).expect("root verifies");
    drop(published);
    assert_eq!(
        root.seal().expect("terminal replay remains stable"),
        terminal
    );
}

#[test]
fn publication_and_destination_receipts_reject_rehashed_identity_substitution() {
    let root = OperationMemoryLease::bounded(4);
    let root_id = root_identity(110);
    let child_id = child_identity(root_id, 111, 0);
    let binding = publication_binding(112);
    root.enable_delegation(root_id, "run-publish-receipt", 1)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-publish-receipt/output", 4)
        .expect("child fits");
    let published = child
        .reserve("output/payload", 4)
        .expect("staging fits")
        .prepare_published_transfer(binding)
        .expect("preparation fits")
        .publish()
        .expect("publication succeeds");
    let receipt = published.receipt().clone();
    receipt
        .verify_for(
            root_id,
            None,
            child_id,
            binding,
            PublishedTransferEnvelope::payload_only(4),
        )
        .expect("exact context verifies");

    let substituted_binding = publication_binding(116);
    let mut substituted = receipt.clone();
    substituted.binding = substituted_binding;
    substituted.receipt_root = substituted.recompute_root();
    assert_eq!(
        substituted
            .verify_for(
                root_id,
                None,
                child_id,
                binding,
                PublishedTransferEnvelope::payload_only(4),
            )
            .expect_err("external binding context rejects a rehashed substitution")
            .reason(),
        "identity"
    );
    let close = published.close().expect("destination closes");
    close.verify_for(&receipt).expect("exact close verifies");
    let mut false_close = close.clone();
    false_close.published_receipt_root = [0; 32];
    false_close.receipt_root = false_close.recompute_root();
    assert_eq!(
        false_close
            .verify_for(&receipt)
            .expect_err("close cannot detach from successful publication")
            .reason(),
        "identity"
    );

    child.close().expect("staging child closes");
    root.seal()
        .expect("root closes")
        .verify_for(root_id)
        .expect("root verifies");
}

#[test]
fn publication_binding_changes_child_and_root_ledger_roots() {
    fn run(binding: PublishedTransferBinding) -> ([u8; 32], [u8; 32], [u8; 32]) {
        let root = OperationMemoryLease::bounded(2);
        let root_id = root_identity(120);
        let child_id = child_identity(root_id, 121, 0);
        root.enable_delegation(root_id, "run-publish-root", 1)
            .expect("root configured");
        let child = root
            .delegate_capacity(child_id, "run-publish-root/output", 2)
            .expect("child fits");
        let published = child
            .reserve("output/payload", 2)
            .expect("staging fits")
            .prepare_published_transfer(binding)
            .expect("preparation fits")
            .publish()
            .expect("publication succeeds");
        let published_root = published.receipt().receipt_root;
        let child_root = child.close().expect("child closes").publication_root;
        let terminal = root.seal().expect("root closes");
        let root_publication_root = terminal.publication_root;
        drop(published);
        (published_root, child_root, root_publication_root)
    }

    let first = run(publication_binding(122));
    let replay = run(publication_binding(122));
    let destination_variant = run(PublishedTransferBinding::new(
        subject(122),
        subject(123),
        subject(124),
        subject(126),
    ));
    assert_eq!(first, replay, "identical identity tuples replay exactly");
    assert_ne!(first.0, destination_variant.0);
    assert_ne!(first.1, destination_variant.1);
    assert_ne!(first.2, destination_variant.2);
}

#[test]
fn root_and_child_receipts_reject_mismatch_mutation_and_rehashed_conservation_errors() {
    let root = OperationMemoryLease::bounded(16);
    let root_id = root_identity(38);
    let child_id = child_identity(root_id, 39, 7);
    root.enable_delegation(root_id, "run-receipt", 1)
        .expect("root configured");
    let child = root
        .delegate_capacity(child_id, "run-receipt/child", 16)
        .expect("child fits");
    let charge = child.reserve("receipt/payload", 16).expect("fits");
    drop(charge);
    let child_receipt = child.close().expect("child closes");
    child_receipt
        .verify_for(root_id, None, child_id)
        .expect("exact child verifies");

    let mut child_under = child_receipt.clone();
    child_under.allocation_returned_bytes -= 1;
    child_under.receipt_root = child_under.recompute_root();
    assert_eq!(
        child_under
            .verify_for(root_id, None, child_id)
            .expect_err("one-under rehashed mutant fails")
            .reason(),
        "conservation"
    );
    let mut child_over = child_receipt.clone();
    child_over.allocation_returned_bytes += 1;
    child_over.receipt_root = child_over.recompute_root();
    assert_eq!(
        child_over
            .verify_for(root_id, None, child_id)
            .expect_err("one-over rehashed mutant fails")
            .reason(),
        "conservation"
    );
    assert_eq!(
        child_receipt
            .verify_for(root_identity(40), None, child_id)
            .expect_err("root mismatch fails")
            .reason(),
        "identity"
    );
    let wrong_child_subject = child_identity(root_id, 41, 7);
    assert_eq!(
        child_receipt
            .verify_for(root_id, None, wrong_child_subject)
            .expect_err("child subject mismatch fails")
            .reason(),
        "identity"
    );
    let wrong_path = child_identity(root_id, 39, 8);
    assert_eq!(
        child_receipt
            .verify_for(root_id, None, wrong_path)
            .expect_err("path mismatch fails")
            .reason(),
        "identity"
    );
    let unrelated_parent = child_identity(root_id, 42, 99);
    assert_eq!(
        child_receipt
            .verify_for(root_id, Some(unrelated_parent), child_id)
            .expect_err("parent mismatch fails")
            .reason(),
        "identity"
    );
    let detached_root = root_identity(208);
    let mut detached = child_receipt.clone();
    detached.root_identity = detached_root;
    detached.receipt_root = detached.recompute_root();
    assert_eq!(
        detached
            .verify_for(detached_root, None, child_id)
            .expect_err("rehashed child cannot detach from its typed root")
            .reason(),
        "identity_relationship"
    );
    let mut false_parent = child_receipt.clone();
    false_parent.parent_identity = Some(unrelated_parent);
    false_parent.receipt_root = false_parent.recompute_root();
    assert_eq!(
        false_parent
            .verify_for(root_id, Some(unrelated_parent), child_id)
            .expect_err("rehashed parent substitution violates the typed path")
            .reason(),
        "identity_relationship"
    );

    let terminal = root.seal().expect("root closes");
    terminal.verify_for(root_id).expect("exact root verifies");
    let mut raw_mutant = terminal.clone();
    raw_mutant.delegated_bytes += 1;
    assert_eq!(
        raw_mutant
            .verify_for(root_id)
            .expect_err("unrehashed mutation fails")
            .reason(),
        "conservation"
    );
    let mut rehashed_mutant = terminal.clone();
    rehashed_mutant.returned_delegated_bytes -= 1;
    rehashed_mutant.receipt_root = rehashed_mutant.recompute_root();
    assert_eq!(
        rehashed_mutant
            .verify_for(root_id)
            .expect_err("rehashed mutation still fails conservation")
            .reason(),
        "conservation"
    );
    assert_eq!(
        terminal
            .verify_for(root_identity(43))
            .expect_err("expected root identity is external context")
            .reason(),
        "identity"
    );
    let mut non_root_terminal = terminal.clone();
    non_root_terminal.root_identity = child_id;
    non_root_terminal.receipt_root = non_root_terminal.recompute_root();
    assert_eq!(
        non_root_terminal
            .verify_for(child_id)
            .expect_err("a child identity cannot verify as the root")
            .reason(),
        "identity"
    );
}

#[test]
fn distinct_root_and_child_identities_have_distinct_receipt_roots() {
    fn run(root_id: LeaseIdentity, child_id: LeaseIdentity) -> ([u8; 32], [u8; 32]) {
        let root = OperationMemoryLease::bounded(1);
        root.enable_delegation(root_id, "run-collision", 1)
            .expect("root configured");
        let child = root
            .delegate_capacity(child_id, "run-collision/child", 1)
            .expect("child fits");
        let child_receipt = child.close().expect("child returns");
        let root_receipt = root.seal().expect("root closes");
        (child_receipt.receipt_root, root_receipt.receipt_root)
    }

    let first_root = root_identity(44);
    let second_root = root_identity(45);
    let first = run(first_root, child_identity(first_root, 46, 0));
    let second = run(second_root, child_identity(second_root, 46, 0));
    assert_ne!(first.0, second.0);
    assert_ne!(first.1, second.1);

    let child_subject_variant = run(first_root, child_identity(first_root, 47, 0));
    let path_variant = run(first_root, child_identity(first_root, 46, 1));
    assert_ne!(first.0, child_subject_variant.0);
    assert_ne!(first.0, path_variant.0);
    assert_ne!(first.1, child_subject_variant.1);
    assert_ne!(first.1, path_variant.1);
}

#[test]
fn deterministic_jsonl_model_replays_byte_for_byte() {
    fn run_model() -> String {
        let root = OperationMemoryLease::bounded(12);
        let root_id = root_identity(48);
        let child_id = child_identity(root_id, 49, 5);
        root.enable_delegation(root_id, "run-jsonl", 2)
            .expect("root configured");
        let child = root
            .delegate_capacity(child_id, "run-jsonl/output", 12)
            .expect("child fits");
        let charge = child.reserve("output/exact", 12).expect("fits");
        let refusal = child
            .reserve("output/one-over", 1)
            .expect_err("one over refuses");
        drop(charge);
        let child_receipt = child.close().expect("child closes");
        let root_receipt = root.seal().expect("root closes");
        format!(
            "{{\"schema\":\"fs-alloc-delegation-e2e-v1\",\"case\":\"exact-and-one-over\",\"sequence\":{},\"root_id\":\"run-jsonl\",\"child_path\":\"run-jsonl/output\",\"declared_capacity\":12,\"refusal_reason\":\"{}\",\"terminal_disposition\":\"returned\",\"no_claim\":\"allocator-overhead-and-addresses\"}}\n{}\n{}",
            refusal.sequence(),
            refusal.reason(),
            child_receipt.to_json(),
            root_receipt.to_json()
        )
    }

    let first = run_model();
    let second = run_model();
    assert_eq!(
        first, second,
        "identical logical traces serialize identically"
    );
    assert_eq!(first.lines().count(), 3);
    assert!(first.contains("\"schema\":\"fs-alloc-delegation-e2e-v1\""));
    assert!(first.contains("\"root_identity\":{"));
    assert!(first.contains("\"parent_subject\":\""));
    assert!(first.contains("\"owner_subject\":\""));
    assert!(first.contains("\"path\":[5]"));
    assert!(first.contains("\"receipt_root\":\""));
    assert!(!first.contains("/Users/"));
    assert!(!first.contains("\"pid\""));
    assert!(!first.contains("\"wall_time\""));
    assert!(!first.contains("0x"));
}

#[test]
fn published_transfer_e2e_log_replays_every_ownership_boundary_byte_for_byte() {
    fn run_model() -> String {
        let root = OperationMemoryLease::bounded(10);
        let root_id = root_identity(130);
        let child_id = child_identity(root_id, 131, 9);
        let binding = publication_binding(132);
        root.enable_delegation(root_id, "run-published-jsonl", 1)
            .expect("root configured");
        let child = root
            .delegate_capacity(child_id, "run-published-jsonl/output", 10)
            .expect("child fits");
        let prepared = child
            .reserve("output/payload", 10)
            .expect("staging fits")
            .prepare_published_transfer(binding)
            .expect("preparation fits");
        let prepared_sequence = prepared.prepared_sequence();
        let published = prepared.publish().expect("publication succeeds");
        let published_receipt = published.receipt().clone();
        let child_receipt = child.close().expect("child closes after publication");
        let root_receipt = root
            .seal()
            .expect("root closes while destination owns bytes");
        let destination_close = published.close().expect("destination closes separately");
        format!(
            "{{\"schema\":\"fs-alloc-published-transfer-e2e-v1\",\"case\":\"successful-output-publication\",\"transitions\":[\"reserve-staging\",\"prepare\",\"publish\",\"close-child\",\"seal-root\",\"close-destination\"],\"prepared_sequence\":{},\"published_sequence\":{},\"closed_sequence\":{},\"bytes\":10,\"conservation\":\"child_granted=child_returned+child_published\",\"equation\":\"10=0+10\",\"binding\":{},\"no_claim\":\"payload-address-or-allocator-overhead\"}}\n{}\n{}\n{}\n{}",
            prepared_sequence,
            published_receipt.published_sequence,
            destination_close.closed_sequence,
            binding.to_json(),
            published_receipt.to_json(),
            child_receipt.to_json(),
            root_receipt.to_json(),
            destination_close.to_json()
        )
    }

    let first = run_model();
    let second = run_model();
    assert_eq!(
        first, second,
        "the full publication lifecycle emits deterministic detailed evidence"
    );
    assert_eq!(first.lines().count(), 5);
    assert!(first.contains("\"transitions\":[\"reserve-staging\",\"prepare\",\"publish\""));
    assert!(first.contains("\"equation\":\"10=0+10\""));
    assert!(first.contains("\"plan_identity\":\""));
    assert!(first.contains("\"occurrence_identity\":\""));
    assert!(first.contains("\"output_identity\":\""));
    assert!(first.contains("\"destination_identity\":\""));
    assert!(first.contains("\"allocation_published_bytes\":10"));
    assert!(first.contains("\"child_published_bytes\":10"));
    assert!(!first.contains("/Users/"));
    assert!(!first.contains("\"pid\""));
    assert!(!first.contains("\"wall_time\""));
    assert!(!first.contains("0x"));
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ModelAction {
    DelegateA,
    DelegateB,
    ReserveA,
    ReleaseA,
    ReturnA,
    ReturnB,
    Seal,
}

#[derive(Clone, Copy, Debug, Default)]
struct Model {
    root_used: u8,
    a_used: u8,
    a_live: bool,
    b_live: bool,
    sealed: bool,
}

impl Model {
    fn apply(&mut self, action: ModelAction) -> bool {
        match action {
            ModelAction::DelegateA if !self.sealed && !self.a_live && self.root_used <= 1 => {
                self.a_live = true;
                self.root_used += 1;
                true
            }
            ModelAction::DelegateB if !self.sealed && !self.b_live && self.root_used <= 1 => {
                self.b_live = true;
                self.root_used += 1;
                true
            }
            ModelAction::ReserveA if !self.sealed && self.a_live && self.a_used == 0 => {
                self.a_used = 1;
                true
            }
            ModelAction::ReleaseA if self.a_used == 1 => {
                self.a_used = 0;
                true
            }
            ModelAction::ReturnA if self.a_live && self.a_used == 0 => {
                self.a_live = false;
                self.root_used -= 1;
                true
            }
            ModelAction::ReturnB if self.b_live => {
                self.b_live = false;
                self.root_used -= 1;
                true
            }
            ModelAction::Seal => {
                self.sealed = true;
                self.root_used == 0 && !self.a_live && !self.b_live
            }
            _ => false,
        }
    }
}

#[test]
fn pure_model_enumerates_short_two_child_interleavings() {
    let actions = [
        ModelAction::DelegateA,
        ModelAction::DelegateB,
        ModelAction::ReserveA,
        ModelAction::ReleaseA,
        ModelAction::ReturnA,
        ModelAction::ReturnB,
        ModelAction::Seal,
    ];
    let mut explored = 0_u64;
    for a in actions {
        for b in actions {
            for c in actions {
                for d in actions {
                    let mut model = Model::default();
                    for action in [a, b, c, d] {
                        let was_sealed = model.sealed;
                        let accepted = model.apply(action);
                        assert!(model.root_used <= 2, "root conservation after {action:?}");
                        assert!(model.a_used <= u8::from(model.a_live));
                        if was_sealed {
                            assert!(
                                !matches!(
                                    action,
                                    ModelAction::DelegateA
                                        | ModelAction::DelegateB
                                        | ModelAction::ReserveA
                                ) || !accepted,
                                "sealed admission stays closed"
                            );
                        }
                    }
                    explored += 1;
                }
            }
        }
    }
    assert_eq!(explored, 2401);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PublicationModelAction {
    Reserve,
    Prepare,
    Publish,
    Rollback,
    Release,
    CloseChild,
    Seal,
    CloseDestination,
}

#[derive(Clone, Copy, Debug)]
struct PublicationModel {
    child_live: bool,
    staging_live: bool,
    prepared: bool,
    destination_live: bool,
    sealed: bool,
    granted: u8,
    returned: u8,
    published: u8,
}

impl Default for PublicationModel {
    fn default() -> Self {
        Self {
            child_live: true,
            staging_live: false,
            prepared: false,
            destination_live: false,
            sealed: false,
            granted: 0,
            returned: 0,
            published: 0,
        }
    }
}

impl PublicationModel {
    fn apply(&mut self, action: PublicationModelAction) -> bool {
        match action {
            PublicationModelAction::Reserve
                if !self.sealed && self.child_live && !self.staging_live && !self.prepared =>
            {
                self.staging_live = true;
                self.granted += 1;
                true
            }
            PublicationModelAction::Prepare
                if !self.sealed && self.staging_live && !self.prepared =>
            {
                self.staging_live = false;
                self.prepared = true;
                true
            }
            PublicationModelAction::Publish if self.prepared => {
                self.prepared = false;
                self.destination_live = true;
                self.published += 1;
                true
            }
            PublicationModelAction::Rollback if self.prepared => {
                self.prepared = false;
                self.returned += 1;
                true
            }
            PublicationModelAction::Release if self.staging_live => {
                self.staging_live = false;
                self.returned += 1;
                true
            }
            PublicationModelAction::CloseChild
                if self.child_live && !self.staging_live && !self.prepared =>
            {
                self.child_live = false;
                true
            }
            PublicationModelAction::Seal => {
                self.sealed = true;
                !self.child_live
            }
            PublicationModelAction::CloseDestination if self.destination_live => {
                self.destination_live = false;
                true
            }
            _ => false,
        }
    }
}

#[test]
fn pure_model_enumerates_publish_rollback_seal_and_destination_close_interleavings() {
    let actions = [
        PublicationModelAction::Reserve,
        PublicationModelAction::Prepare,
        PublicationModelAction::Publish,
        PublicationModelAction::Rollback,
        PublicationModelAction::Release,
        PublicationModelAction::CloseChild,
        PublicationModelAction::Seal,
        PublicationModelAction::CloseDestination,
    ];
    let mut explored = 0_u64;
    for a in actions {
        for b in actions {
            for c in actions {
                for d in actions {
                    let mut model = PublicationModel::default();
                    for action in [a, b, c, d] {
                        let was_sealed = model.sealed;
                        let accepted = model.apply(action);
                        assert_eq!(
                            model.granted,
                            model.returned
                                + model.published
                                + u8::from(model.staging_live)
                                + u8::from(model.prepared),
                            "every granted byte has one exact disposition after {action:?}"
                        );
                        assert!(
                            !(model.staging_live && model.prepared),
                            "staging and prepared ownership are mutually exclusive"
                        );
                        if was_sealed
                            && matches!(
                                action,
                                PublicationModelAction::Reserve | PublicationModelAction::Prepare
                            )
                        {
                            assert!(!accepted, "the seal cut forbids new admission");
                        }
                    }
                    explored += 1;
                }
            }
        }
    }
    assert_eq!(explored, 4096);
}
