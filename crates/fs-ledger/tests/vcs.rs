//! Physics-VCS conformance (addendum Proposal 10 base verbs, the lmp4.9
//! bead). Acceptance: commits are reproducible Merkle roots; branches
//! share artifacts by hash (storage audit: N branches ≈ 1× + deltas);
//! checkout of a nearby branch costs a delta, not a full re-solve;
//! branch bookkeeping supports the diff/bisect/merge consumers; GC never
//! collects anything reachable from a live branch; boundaries (unknown
//! commit, empty repo) are structured.

use fs_ledger::vcs::Vcs;
use fs_ledger::{
    EdgeRole, EventRow, FiveExplicits, Ledger, LedgerError, OpOutcome, travel::ExecMode,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/vcs\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn open_ledger(tag: &str) -> (Ledger, std::path::PathBuf) {
    let dir = std::env::temp_dir().join(format!("fs-vcs-{tag}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let ledger = Ledger::open(dir.join("v.led").to_str().expect("utf8")).expect("ledger");
    (ledger, dir)
}

const EXPLICITS: FiveExplicits<'_> = FiveExplicits {
    seed: b"\x5e\xed",
    versions: "{\"constellation\":\"2026-07\"}",
    budget: "{\"wall_s\":10}",
    capability: "{\"ops\":[\"*\"]}",
};

/// Run one op on a branch: begin, store an artifact, link, finish.
fn run_op(ledger: &Ledger, branch: i64, ir: &str, payload: &[u8], t: i64) -> i64 {
    let op = ledger
        .begin_op_on(branch, ExecMode::Deterministic, None, ir, &EXPLICITS, t)
        .expect("begin");
    let receipt = ledger
        .put_artifact("result", payload, None)
        .expect("artifact");
    ledger
        .link(op, &receipt.hash, fs_ledger::EdgeRole::Out)
        .expect("link");
    ledger
        .finish_op(op, OpOutcome::Ok, None, t + 1)
        .expect("finish");
    op
}

#[test]
fn vc_001_commits_are_reproducible_merkle_roots() {
    let (ledger_a, dir_a) = open_ledger("repro-a");
    let (ledger_b, dir_b) = open_ledger("repro-b");
    // Two ledgers, identical LOGICAL histories, different wall times.
    for (ledger, t0) in [(&ledger_a, 100i64), (&ledger_b, 9_999i64)] {
        run_op(ledger, 1, "{\"op\":\"solve\",\"dof\":64}", b"field-v1", t0);
        run_op(
            ledger,
            1,
            "{\"op\":\"optimize\",\"steps\":3}",
            b"design-v1",
            t0 + 50,
        );
    }
    let mut vcs_a = Vcs::new();
    let mut vcs_b = Vcs::new();
    let ca = vcs_a.commit(&ledger_a, 1).expect("commit a");
    let cb = vcs_b.commit(&ledger_b, 1).expect("commit b");
    assert_eq!(
        ca.root, cb.root,
        "identical logical histories yield identical roots (wall times excluded)"
    );
    // Recommit without changes: the state identity is idempotent and must not
    // create a self-parent cycle or duplicate event.
    let ca2 = vcs_a.commit(&ledger_a, 1).expect("recommit");
    assert_eq!(ca2.root, ca.root, "same state, same root");
    assert_eq!(ca2, ca, "unchanged recommit is the same commit");
    assert_eq!(
        ledger_a.table_count("events").expect("event count"),
        2,
        "exactly one vcs-identity + one vcs-commit event: the idempotent \
         recommit appended neither a duplicate commit nor a new identity"
    );
    // A new op changes the root.
    run_op(
        &ledger_a,
        1,
        "{\"op\":\"solve\",\"dof\":128}",
        b"field-v2",
        300,
    );
    let ca3 = vcs_a.commit(&ledger_a, 1).expect("commit 3");
    assert_ne!(ca3.root, ca.root, "state change changes the root");
    assert_eq!(ca3.parent, Some(ca.root), "changed states chain");
    // Different artifact content changes the root even with the same IR.
    let (ledger_c, dir_c) = open_ledger("repro-c");
    run_op(
        &ledger_c,
        1,
        "{\"op\":\"solve\",\"dof\":64}",
        b"DIFFERENT",
        100,
    );
    run_op(
        &ledger_c,
        1,
        "{\"op\":\"optimize\",\"steps\":3}",
        b"design-v1",
        150,
    );
    let mut vcs_c = Vcs::new();
    let cc = vcs_c.commit(&ledger_c, 1).expect("commit c");
    assert_ne!(cc.root, ca.root, "artifact hashes are folded into leaves");
    for d in [dir_a, dir_b, dir_c] {
        let _ = std::fs::remove_dir_all(&d);
    }
    verdict(
        "vc-001",
        "roots reproduce across ledgers and runs; unchanged commits are idempotent; changed \
         states chain and move on any state or artifact change",
    );
}

#[test]
fn vc_002_branches_share_artifacts_by_hash() {
    let (ledger, dir) = open_ledger("share");
    // A common history of 6 artifacts…
    for k in 0..6 {
        run_op(
            &ledger,
            1,
            &format!("{{\"op\":{k}}}"),
            format!("shared-{k}").as_bytes(),
            10 + k,
        );
    }
    // …then 4 branches, each adding ONE delta op.
    let mut branch_ids = vec![1i64];
    for b in 0..4 {
        let id = ledger.fork(&format!("exp-{b}"), 1).expect("fork");
        run_op(
            &ledger,
            id,
            "{\"op\":\"delta\"}",
            format!("delta-{b}").as_bytes(),
            100 + b,
        );
        branch_ids.push(id);
    }
    let vcs = Vcs::new();
    let audit = vcs.storage_audit(&ledger).expect("audit");
    assert_eq!(audit.branches, 5);
    // Physical: 6 shared + 4 deltas = 10. Logical: 5 branches × 6 + 4.
    assert_eq!(audit.physical_artifacts, 10, "N branches ≈ 1x + deltas");
    assert_eq!(audit.logical_references, 5 * 6 + 4);
    assert!(
        audit.logical_references > 3 * audit.physical_artifacts,
        "sharing is real: {} logical over {} physical",
        audit.logical_references,
        audit.physical_artifacts
    );
    let _ = branch_ids;
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-002",
        "storage audit: 10 physical rows serve 34 logical references",
    );
}

#[test]
fn vc_003_checkout_and_nearby_delta() {
    let (ledger, dir) = open_ledger("checkout");
    for k in 0..20 {
        run_op(
            &ledger,
            1,
            &format!("{{\"base\":{k}}}"),
            format!("b{k}").as_bytes(),
            10 + k,
        );
    }
    let exp = ledger.fork("exp", 1).expect("fork");
    run_op(&ledger, exp, "{\"tweak\":1}", b"t1", 100);
    run_op(&ledger, exp, "{\"tweak\":2}", b"t2", 101);
    let mut vcs = Vcs::new();
    let main_head = vcs.commit(&ledger, 1).expect("main commit");
    let exp_head = vcs.commit(&ledger, exp).expect("exp commit");
    // Checkout materializes the committed view.
    let snap = vcs
        .checkout(&ledger, exp, &exp_head.root)
        .expect("checkout");
    assert_eq!(snap.ops.len(), 22, "20 shared + 2 branch ops");
    // The nearby delta is 2 ops, with 20 shared — a delta-solve frontier,
    // not a full re-solve.
    let delta = vcs
        .checkout_delta(&main_head.id(), &exp_head.id())
        .expect("delta");
    assert_eq!(delta.added.len(), 2, "delta-solve frontier");
    assert_eq!(delta.removed.len(), 0);
    assert_eq!(delta.shared, 20, "the bulk is hash-shared");
    assert!(
        delta.added.len() * 5 < delta.shared,
        "checkout cost is a small fraction of a full re-solve"
    );
    // Boundary: an unknown root is a structured error.
    let bogus = fs_ledger::hash_bytes(b"no-such-commit");
    let err = vcs.checkout(&ledger, 1, &bogus).expect_err("must refuse");
    assert!(format!("{err}").contains("no commit"), "teaches: {err}");
    // Boundary: empty repo commits cleanly (the well-defined empty root).
    let (empty, edir) = open_ledger("empty");
    let mut evcs = Vcs::new();
    let e = evcs.commit(&empty, 1).expect("empty commit");
    assert!(e.frontier_op.is_none());
    let esnap = evcs.checkout(&empty, 1, &e.root).expect("checkout empty");
    assert!(esnap.ops.is_empty());
    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&edir);
    verdict(
        "vc-003",
        "checkout materializes committed views; nearby delta = 2 ops over 20 shared; \
         unknown root refused; empty repo well-defined",
    );
}

#[test]
fn vc_004_merge_views_and_branch_independence() {
    let (ledger, dir) = open_ledger("merge");
    for k in 0..5 {
        run_op(
            &ledger,
            1,
            &format!("{{\"base\":{k}}}"),
            format!("b{k}").as_bytes(),
            10 + k,
        );
    }
    let a = ledger.fork("a", 1).expect("fork a");
    let b = ledger.fork("b", 1).expect("fork b");
    let op_a = run_op(&ledger, a, "{\"side\":\"a\"}", b"a1", 100);
    let op_b1 = run_op(&ledger, b, "{\"side\":\"b1\"}", b"b1", 110);
    let op_b2 = run_op(&ledger, b, "{\"side\":\"b2\"}", b"b2", 111);
    let vcs = Vcs::new();
    let views = vcs.merge_views(&ledger, a, b).expect("views");
    assert_eq!(views.base.len(), 5, "merge base is the shared history");
    assert_eq!(views.only_a, vec![op_a]);
    assert_eq!(views.only_b, vec![op_b1, op_b2]);
    // Branch independence (inherits fork-independence): A's op never
    // appears in B's view and vice versa.
    let ops_b = ledger.visible_op_ids(b, None).expect("b view");
    assert!(!ops_b.contains(&op_a), "A's work never leaks into B");
    let ops_a = ledger.visible_op_ids(a, None).expect("a view");
    assert!(!ops_a.contains(&op_b1) && !ops_a.contains(&op_b2));
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-004",
        "merge views split base/only-A/only-B exactly; branches independent",
    );
}

#[test]
fn vc_005_gc_never_collects_live_branch_artifacts() {
    let (ledger, dir) = open_ledger("gc");
    // Many-branch fixture: shared history + per-branch deltas.
    let mut expected_alive = Vec::new();
    for k in 0..4 {
        run_op(
            &ledger,
            1,
            &format!("{{\"base\":{k}}}"),
            format!("keep-b{k}").as_bytes(),
            10 + k,
        );
        expected_alive.push(fs_ledger::hash_bytes(format!("keep-b{k}").as_bytes()));
    }
    for br in 0..6 {
        let id = ledger.fork(&format!("live-{br}"), 1).expect("fork");
        run_op(
            &ledger,
            id,
            "{\"op\":\"delta\"}",
            format!("keep-d{br}").as_bytes(),
            100 + br,
        );
        expected_alive.push(fs_ledger::hash_bytes(format!("keep-d{br}").as_bytes()));
    }
    // An UNREFERENCED artifact (no op links it): the only legal victim.
    let orphan = ledger
        .put_artifact("orphan", b"collect-me", None)
        .expect("orphan");
    let report = ledger.gc_unreferenced_artifacts(false).expect("gc");
    assert!(report.deleted >= 1, "the orphan is collectable");
    // THE G0 DATA-LOSS INVARIANT: every artifact reachable from any live
    // branch survives.
    for h in &expected_alive {
        assert!(
            ledger.get_artifact(h).expect("query").is_some(),
            "live-branch artifact {} must survive GC",
            h.to_hex()
        );
    }
    assert!(
        ledger.get_artifact(&orphan.hash).expect("query").is_none(),
        "the orphan is gone"
    );
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-005",
        "GC collected only the unreferenced orphan; all 10 live-branch artifacts survive",
    );
}

fn identity_op(
    ledger: &Ledger,
    mode: ExecMode,
    seed: &[u8],
    edge: Option<(EdgeRole, &[u8])>,
) -> i64 {
    let explicits = FiveExplicits { seed, ..EXPLICITS };
    let op = ledger
        .begin_op_on(1, mode, None, "{\"op\":\"identity\"}", &explicits, 1)
        .expect("begin identity op");
    if let Some((role, payload)) = edge {
        let artifact = ledger
            .put_artifact("identity", payload, None)
            .expect("identity artifact");
        ledger
            .link(op, &artifact.hash, role)
            .expect("identity edge");
    }
    ledger
        .finish_op(op, OpOutcome::Ok, None, 2)
        .expect("finish identity op");
    op
}

#[test]
fn vc_006_commit_identity_is_framed_role_qualified_and_mode_bound() {
    let framed_a = Ledger::open(":memory:").expect("framed a");
    let framed_b = Ledger::open(":memory:").expect("framed b");
    let artifact_hash = fs_ledger::hash_bytes(b"edge-payload");
    let op_a = identity_op(
        &framed_a,
        ExecMode::Deterministic,
        b"seed",
        Some((EdgeRole::Out, b"edge-payload")),
    );
    let mut boundary_alias_seed = b"seed".to_vec();
    boundary_alias_seed.extend_from_slice(artifact_hash.as_bytes());
    let op_b = identity_op(
        &framed_b,
        ExecMode::Deterministic,
        &boundary_alias_seed,
        None,
    );
    assert_ne!(
        framed_a.commit_leaf(op_a).expect("framed leaf a"),
        framed_b.commit_leaf(op_b).expect("framed leaf b"),
        "seed/edge field boundaries must not alias"
    );

    let role_in = Ledger::open(":memory:").expect("role in");
    let role_out = Ledger::open(":memory:").expect("role out");
    let in_op = identity_op(
        &role_in,
        ExecMode::Deterministic,
        b"seed",
        Some((EdgeRole::In, b"same-artifact")),
    );
    let out_op = identity_op(
        &role_out,
        ExecMode::Deterministic,
        b"seed",
        Some((EdgeRole::Out, b"same-artifact")),
    );
    assert_ne!(
        role_in.commit_leaf(in_op).expect("input leaf"),
        role_out.commit_leaf(out_op).expect("output leaf"),
        "edge role is semantic"
    );

    let deterministic = Ledger::open(":memory:").expect("deterministic");
    let fast = Ledger::open(":memory:").expect("fast");
    let deterministic_op = identity_op(&deterministic, ExecMode::Deterministic, b"seed", None);
    let fast_op = identity_op(&fast, ExecMode::Fast, b"seed", None);
    assert_ne!(
        deterministic
            .commit_leaf(deterministic_op)
            .expect("deterministic leaf"),
        fast.commit_leaf(fast_op).expect("fast leaf"),
        "execution mode is semantic"
    );
    verdict(
        "vc-006",
        "commit identity separates frames, edge roles, and execution modes",
    );
}

#[test]
fn vc_007_checkout_is_exactly_the_frozen_commit_view() {
    let (ledger, dir) = open_ledger("frozen-checkout");
    let first_op = run_op(&ledger, 1, "{\"op\":\"first\"}", b"first", 10);
    let mut vcs = Vcs::new();
    let first_commit = vcs.commit(&ledger, 1).expect("first commit");

    let late_link = ledger
        .put_artifact("result", b"linked-after-commit", None)
        .expect("late artifact");
    ledger
        .link(first_op, &late_link.hash, EdgeRole::Out)
        .expect("late edge");
    run_op(&ledger, 1, "{\"op\":\"future\"}", b"future", 20);

    let frozen = vcs
        .checkout(&ledger, 1, &first_commit.root)
        .expect("frozen checkout");
    assert_eq!(frozen.ops.len(), 1, "post-commit op leaked into checkout");
    assert_eq!(
        frozen.artifacts,
        vec![fs_ledger::hash_bytes(b"first")],
        "post-commit op or late edge leaked a future artifact"
    );

    let in_flight = ledger
        .begin_op_on(1, ExecMode::Deterministic, None, "{}", &EXPLICITS, 30)
        .expect("begin in-flight op");
    let error = vcs
        .commit(&ledger, 1)
        .expect_err("in-flight state cannot be frozen reproducibly");
    assert!(
        error
            .to_string()
            .contains(&format!("op {in_flight} is still in flight")),
        "unexpected refusal: {error}"
    );
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-007",
        "checkout excludes later ops and later edges; in-flight commits fail closed",
    );
}

#[test]
fn vc_008_equal_roots_do_not_clobber_branch_envelopes() {
    // gp3.17: a fork with NO divergence reaches the SAME semantic root as
    // its parent branch — the registry must hold BOTH commit envelopes.
    let (ledger, dir) = open_ledger("same-root");
    run_op(&ledger, 1, "{\"op\":\"solve\"}", b"field", 10);
    run_op(&ledger, 1, "{\"op\":\"optimize\"}", b"design", 20);
    let twin = ledger.fork("twin", 1).expect("fork");
    let mut vcs = Vcs::new();
    let main_commit = vcs.commit(&ledger, 1).expect("main commit");
    let twin_commit = vcs.commit(&ledger, twin).expect("twin commit");
    assert_eq!(
        main_commit.root, twin_commit.root,
        "undiverged fork shares the semantic root"
    );
    assert_ne!(
        main_commit.id(),
        twin_commit.id(),
        "but the envelopes are distinct commits"
    );
    // COLLISION: both envelopes are retrievable with THEIR OWN branch
    // metadata (the root-keyed registry clobbered one of these).
    let main_lookup = vcs.lookup(&main_commit.id()).expect("main envelope");
    let twin_lookup = vcs.lookup(&twin_commit.id()).expect("twin envelope");
    assert_eq!(main_lookup.branch, 1);
    assert_eq!(twin_lookup.branch, twin);
    assert_eq!(vcs.lookup_semantic(&main_commit.root).len(), 2);
    // Heads are branch-scoped, both live.
    assert_eq!(vcs.head(&ledger, 1).expect("head"), Some(main_commit.root));
    assert_eq!(
        vcs.head(&ledger, twin).expect("head"),
        Some(twin_commit.root)
    );
    // Checkout is envelope-scoped and works for both.
    assert_eq!(
        vcs.checkout(&ledger, 1, &main_commit.root)
            .expect("m")
            .branch,
        1
    );
    assert_eq!(
        vcs.checkout(&ledger, twin, &twin_commit.root)
            .expect("t")
            .branch,
        twin
    );
    // REPLAY: a fresh registry rebuilds identical envelopes deterministically.
    let mut replay = Vcs::new();
    let main_replay = replay.commit(&ledger, 1).expect("replay main");
    let twin_replay = replay.commit(&ledger, twin).expect("replay twin");
    assert_eq!(main_replay.id(), main_commit.id());
    assert_eq!(twin_replay.id(), twin_commit.id());
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-008",
        "equal semantic roots coexist as distinct branch envelopes; replay deterministic",
    );
}

#[test]
fn vc_009_cross_ledger_deltas_use_semantic_leaves_not_row_ids() {
    // gp3.17: two INDEPENDENT ledgers holding the same logical ops in a
    // DIFFERENT insertion order have different local row ids (and, since
    // the root binds sequence, different roots) — but the semantic delta
    // between them is EMPTY. The old id-based delta reported phantom
    // adds/removes here.
    let (ledger_a, dir_a) = open_ledger("cross-a");
    let (ledger_b, dir_b) = open_ledger("cross-b");
    run_op(&ledger_a, 1, "{\"op\":\"solve\"}", b"field", 10);
    run_op(&ledger_a, 1, "{\"op\":\"optimize\"}", b"design", 20);
    // REORDERED IMPORT into b: same two logical ops, swapped row order.
    run_op(&ledger_b, 1, "{\"op\":\"optimize\"}", b"design", 500);
    run_op(&ledger_b, 1, "{\"op\":\"solve\"}", b"field", 510);
    assert_ne!(
        ledger_a.vcs_identity().expect("id a"),
        ledger_b.vcs_identity().expect("id b"),
        "independent databases carry distinct persisted identities"
    );
    let mut vcs = Vcs::new();
    let ca = vcs.commit(&ledger_a, 1).expect("commit a");
    let cb = vcs.commit(&ledger_b, 1).expect("commit b");
    assert_ne!(ca.root, cb.root, "the root binds the SEQUENCE of history");
    let delta = vcs.checkout_delta(&ca.id(), &cb.id()).expect("cross delta");
    assert!(
        delta.added.is_empty() && delta.removed.is_empty(),
        "the SET of semantic leaves is identical — no phantom frontier: {delta:?}"
    );
    assert_eq!(delta.shared, 2, "both ops are hash-shared semantically");
    // A REAL divergence reports the target-local op id with its leaf.
    let z = run_op(&ledger_b, 1, "{\"op\":\"probe\"}", b"probe", 600);
    let cb2 = vcs.commit(&ledger_b, 1).expect("commit b2");
    let delta = vcs.checkout_delta(&ca.id(), &cb2.id()).expect("delta 2");
    assert_eq!(delta.added.len(), 1);
    assert_eq!(delta.added[0].local_op, z, "target-local id, b's own");
    assert_eq!(
        delta.added[0].leaf,
        ledger_b.commit_leaf(z).expect("leaf"),
        "identified by the portable semantic leaf"
    );
    assert!(delta.removed.is_empty());
    let _ = std::fs::remove_dir_all(&dir_a);
    let _ = std::fs::remove_dir_all(&dir_b);
    verdict(
        "vc-009",
        "cross-ledger + reordered-import deltas: empty semantic frontier, no row-id lies",
    );
}

#[test]
fn vc_010_public_events_cannot_preempt_ledger_identity() {
    let (ledger, dir) = open_ledger("reserved-identity");
    let forged_payload = format!(
        "{{\"kind\":\"vcs-identity\",\"identity\":\"{}\"}}",
        "00".repeat(32)
    );
    let forged = EventRow {
        session: None,
        t: 0,
        kind: "vcs-identity",
        payload: Some(&forged_payload),
    };
    assert!(matches!(
        ledger.append_event(&forged),
        Err(LedgerError::Invalid { field, .. }) if field == "kind"
    ));

    let ordinary = EventRow {
        session: None,
        t: 1,
        kind: "observer-note",
        payload: Some("{}"),
    };
    assert!(ledger.append_events(&[ordinary, forged]).is_err());
    assert_eq!(
        ledger.table_count("events").expect("event count"),
        0,
        "a reserved kind must roll back the entire public batch"
    );

    let identity = ledger.vcs_identity().expect("internal identity mint");
    assert_ne!(identity.to_hex(), "00".repeat(32));
    assert_eq!(ledger.table_count("events").expect("event count"), 1);
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "vc-010",
        "reserved VCS identity cannot be injected through single or batched public events",
    );
}
