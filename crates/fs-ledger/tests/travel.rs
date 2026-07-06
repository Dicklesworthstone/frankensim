//! fs-ledger time-travel conformance suite (plan §11.2; the gp3.3 bead).
//!
//! Acceptance criteria covered: replay of a deterministic fixture
//! reproduces every artifact hash; forks share artifacts (storage audit);
//! explain() reconstructs known lineage completely and fails LOUDLY on
//! orphan inputs; at(t) views are consistent at arbitrary interior
//! instants; kill -9 during fork traffic recovers to a lint-clean state.

use std::sync::atomic::{AtomicU32, Ordering};

use fs_ledger::{
    ContentHash, EdgeRole, ExecMode, FiveExplicits, Ledger, MAIN_BRANCH, OpOutcome,
    SCHEMA_VERSION,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!("fs-ledger-travel-{tag}-{}-{n}.db", std::process::id()))
        .display()
        .to_string()
}

fn cleanup_db(path: &str) {
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/travel\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

const FX: FiveExplicits<'static> = FiveExplicits {
    seed: &[0x5E, 0xED],
    versions: r#"{"constellation":"f92683cc4572a198"}"#,
    budget: r#"{"wall_s":10}"#,
    capability: r#"{"ops":["travel.*"]}"#,
};

/// One complete op on `branch` producing one deterministic artifact.
fn unit(l: &Ledger, branch: i64, mode: ExecMode, tag: u64, t: i64) -> ContentHash {
    let op = l
        .begin_op_on(branch, mode, Some(b"travel".as_slice()), &format!("{{\"tag\":{tag}}}"), &FX, t)
        .expect("begin");
    let bytes: Vec<u8> = (0..256u64).map(|i| ((tag * 197 + i * 3) % 251) as u8).collect();
    let receipt = l.put_artifact("fixture", &bytes, None).expect("put");
    l.link(op, &receipt.hash, EdgeRole::Out).expect("link");
    l.finish_op(op, OpOutcome::Ok, None, t + 5).expect("finish");
    receipt.hash
}

#[test]
fn tt_001_v1_database_migrates_to_v2_with_history_intact() {
    let db = temp_db("migrate");
    {
        // Construct a GENUINE v1 database: v1 DDL only, one op, version 1.
        let raw = fsqlite::Connection::open(&db).expect("raw open");
        for ddl in fs_ledger::schema::V1 {
            raw.execute(ddl).expect("v1 ddl");
        }
        raw.execute(
            "INSERT INTO ops(session, ir, seed, versions, budget, capability, t_start, \
             t_end, outcome) VALUES (NULL, '{}', X'AB', '{}', '{}', '{}', 1, 2, 'ok')",
        )
        .expect("v1 op");
        raw.execute("PRAGMA user_version = 1").expect("set v1");
    }
    let l = Ledger::open(&db).expect("open migrates");
    assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
    // Pre-v2 history reads as main-branch deterministic ops.
    assert_eq!(l.visible_op_ids(MAIN_BRANCH, None).unwrap(), vec![1]);
    let main = l.branch(MAIN_BRANCH).unwrap().expect("main branch seeded");
    assert_eq!(main.name, "main");
    assert!(l.lint().unwrap().is_clean());
    // And the ledger keeps working post-migration.
    unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 7, 10);
    assert_eq!(l.visible_op_ids(MAIN_BRANCH, None).unwrap(), vec![1, 2]);
    drop(l);
    cleanup_db(&db);
    verdict("tt-001", "genuine v1 db migrated to v2; history on main; lint clean");
}

#[test]
fn tt_002_forks_share_artifacts_storage_audit() {
    let db = temp_db("forks");
    let l = Ledger::open(&db).expect("open");
    let mut mains = Vec::new();
    for tag in 0..5 {
        mains.push(unit(&l, MAIN_BRANCH, ExecMode::Deterministic, tag, 10 + tag as i64 * 10));
    }
    assert_eq!(l.table_count("artifacts").unwrap(), 5);
    // Three speculative forks; each consumes shared artifacts and adds ONE
    // new artifact. Storage must grow by exactly the deltas.
    for k in 0..3u64 {
        let branch = l.fork(&format!("candidate-{k}"), MAIN_BRANCH).unwrap();
        let op = l
            .begin_op_on(branch, ExecMode::Deterministic, None, "{\"fork\":true}", &FX, 100)
            .unwrap();
        l.link(op, &mains[k as usize], EdgeRole::In).unwrap();
        let out = l.put_artifact("fixture", format!("delta-{k}").as_bytes(), None).unwrap();
        l.link(op, &out.hash, EdgeRole::Out).unwrap();
        l.finish_op(op, OpOutcome::Ok, None, 105).unwrap();
    }
    assert_eq!(
        l.table_count("artifacts").unwrap(),
        8,
        "3 forks over a 5-artifact study must cost exactly 3 deltas (1x + deltas)"
    );
    // Fork independence: each fork sees the shared prefix + its own op only.
    let a = l.branch_by_name("candidate-0").unwrap().unwrap().id;
    let b = l.branch_by_name("candidate-1").unwrap().unwrap().id;
    let diff = l.branch_diff(a, b).unwrap();
    assert_eq!(diff.shared, 5);
    assert_eq!(diff.only_a.len(), 1);
    assert_eq!(diff.only_b.len(), 1);
    // Ops on A never affect B's view (property over further writes).
    let before_b = l.visible_op_ids(b, None).unwrap();
    for tag in 20..25 {
        unit(&l, a, ExecMode::Deterministic, tag, 200 + tag as i64);
    }
    assert_eq!(l.visible_op_ids(b, None).unwrap(), before_b, "branch A writes leaked into B");
    assert!(l.lint().unwrap().is_clean());
    drop(l);
    cleanup_db(&db);
    verdict("tt-002", "N forks cost 1x artifacts + deltas; branch views independent");
}

/// Deterministic fixture study: 6 units; optionally perturb one unit's
/// artifact content, optionally run one op in fast mode.
fn build_fixture(db: &str, perturb: Option<u64>, fast_at: Option<u64>) -> Ledger {
    let l = Ledger::open(db).expect("open fixture");
    for tag in 0..6u64 {
        let mode =
            if fast_at == Some(tag) { ExecMode::Fast } else { ExecMode::Deterministic };
        let op = l
            .begin_op_on(MAIN_BRANCH, mode, None, &format!("{{\"tag\":{tag}}}"), &FX, tag as i64)
            .expect("begin");
        let salt = if perturb == Some(tag) { 1u64 } else { 0 };
        let bytes: Vec<u8> =
            (0..256u64).map(|i| ((tag * 197 + i * 3 + salt) % 251) as u8).collect();
        let receipt = l.put_artifact("fixture", &bytes, None).expect("put");
        l.link(op, &receipt.hash, EdgeRole::Out).expect("link");
        l.finish_op(op, OpOutcome::Ok, None, tag as i64 + 1).expect("finish");
    }
    l
}

#[test]
fn tt_003_replay_audit_battery() {
    let (da, db_, dc, dd) = (temp_db("rp-a"), temp_db("rp-b"), temp_db("rp-c"), temp_db("rp-d"));
    let original = build_fixture(&da, None, None);
    // Faithful replay: every artifact hash reproduced exactly.
    let replay = build_fixture(&db_, None, None);
    let v = original.replay_verdict(MAIN_BRANCH, &replay, MAIN_BRANCH).unwrap();
    assert!(v.is_replay_clean(), "faithful replay flagged: {v:?}");
    assert_eq!(v.compared, 6);
    // A deterministic op producing different bytes is a replay FAILURE.
    let broken = build_fixture(&dc, Some(3), None);
    let v = original.replay_verdict(MAIN_BRANCH, &broken, MAIN_BRANCH).unwrap();
    assert!(!v.is_replay_clean());
    assert_eq!(v.deterministic_mismatches.len(), 1);
    assert_eq!(v.deterministic_mismatches[0].position, 3);
    // The same divergence on a FAST op is reported but does not fail.
    let original_fast = build_fixture(&dd, None, Some(3));
    let de = temp_db("rp-e");
    let fast_diverged = build_fixture(&de, Some(3), Some(3));
    let v = original_fast.replay_verdict(MAIN_BRANCH, &fast_diverged, MAIN_BRANCH).unwrap();
    assert!(v.is_replay_clean(), "fast divergence must not fail the audit");
    assert_eq!(v.fast_divergences.len(), 1);
    for d in [&da, &db_, &dc, &dd, &de] {
        cleanup_db(d);
    }
    verdict("tt-003", "replay clean/deterministic-failure/fast-divergence all classified");
}

#[test]
fn tt_004_explain_reconstructs_full_lineage() {
    let db = temp_db("explain");
    let l = Ledger::open(&db).expect("open");
    // Chain: op1 -> A; op2(A) -> B; op3(B) -> C.
    let a = unit(&l, MAIN_BRANCH, ExecMode::Deterministic, 1, 10);
    let op2 = l.begin_op_on(MAIN_BRANCH, ExecMode::Deterministic, None, "{}", &FX, 20).unwrap();
    l.link(op2, &a, EdgeRole::In).unwrap();
    let b = l.put_artifact("mid", b"artifact B", None).unwrap();
    l.link(op2, &b.hash, EdgeRole::Out).unwrap();
    l.finish_op(op2, OpOutcome::Ok, None, 25).unwrap();
    let op3 = l.begin_op_on(MAIN_BRANCH, ExecMode::Fast, None, "{}", &FX, 30).unwrap();
    l.link(op3, &b.hash, EdgeRole::In).unwrap();
    let c = l.put_artifact("final", b"artifact C", None).unwrap();
    l.link(op3, &c.hash, EdgeRole::Out).unwrap();
    l.finish_op(op3, OpOutcome::Ok, None, 35).unwrap();

    let tree = l.explain(&c.hash, 16).unwrap().expect("explainable");
    assert!(!tree.truncated);
    assert_eq!(tree.produced_by.len(), 1);
    let op3_node = &tree.produced_by[0];
    assert_eq!(op3_node.exec_mode, "fast");
    let b_node = &op3_node.inputs[0];
    let op2_node = &b_node.produced_by[0];
    let a_node = &op2_node.inputs[0];
    assert_eq!(a_node.hash_hex, a.to_hex(), "lineage reaches the root input");
    assert!(a_node.produced_by[0].inputs.is_empty(), "root op has no inputs");
    // Renderings carry the story.
    let json = tree.to_json();
    assert!(json.contains(&a.to_hex()) && json.contains("\"exec_mode\":\"fast\""));
    assert!(tree.render_text().contains("<- op"));
    // Depth limiting truncates gracefully instead of exploding.
    let shallow = l.explain(&c.hash, 1).unwrap().unwrap();
    assert!(shallow.produced_by[0].inputs[0].truncated);
    // Orphan input (artifact deleted behind the ledger's back) fails LOUDLY.
    {
        let raw = fsqlite::Connection::open(&db).expect("raw");
        raw.query("PRAGMA foreign_keys=OFF").expect("fk off");
        raw.execute_with_params(
            "DELETE FROM artifacts WHERE hash = ?1",
            &[fsqlite::SqliteValue::Blob(a.as_bytes().to_vec().into())],
        )
        .expect("corrupt");
    }
    let err = l.explain(&c.hash, 16).unwrap_err();
    assert_eq!(err.code(), "LedgerCorruption", "orphan input must be loud: {err}");
    drop(l);
    cleanup_db(&db);
    verdict("tt-004", "full causal tree reconstructed; depth-limits; orphan inputs loud");
}

#[test]
fn tt_005_at_time_views_consistent_at_interior_instants() {
    let db = temp_db("attime");
    let l = Ledger::open(&db).expect("open");
    for tag in 0..8u64 {
        unit(&l, MAIN_BRANCH, ExecMode::Deterministic, tag, (tag as i64) * 100);
    }
    // Sweep cutoffs including instants strictly inside op lifetimes
    // (t_start = 100k, t_end = 100k + 5).
    let mut prev_ops = 0usize;
    let mut prev_artifacts = 0usize;
    for cutoff in [-1i64, 0, 3, 99, 102, 350, 703, 9_999] {
        let view = l.at_time(MAIN_BRANCH, cutoff).unwrap();
        assert!(view.ops.len() >= prev_ops, "op visibility must be monotone in t");
        assert!(view.artifacts.len() >= prev_artifacts, "artifact visibility monotone");
        // Each fixture op produces exactly one artifact, so visible
        // artifacts never exceed FINISHED visible ops.
        assert!(view.artifacts.len() <= view.ops.len() - view.in_flight);
        // Internal consistency: every in-flight op is outcome-masked.
        for op in &view.ops {
            if op.t_end.is_none() {
                assert!(op.outcome.is_none(), "mid-sweep op leaked its future outcome");
            }
        }
        prev_ops = view.ops.len();
        prev_artifacts = view.artifacts.len();
    }
    // Interior instant: t=3 shows op1 in flight (t_start 0, t_end 5) with
    // no artifact yet.
    let mid = l.at_time(MAIN_BRANCH, 3).unwrap();
    assert_eq!(mid.ops.len(), 1);
    assert_eq!(mid.in_flight, 1);
    assert!(mid.artifacts.is_empty(), "unfinished op's output is the future");
    drop(l);
    cleanup_db(&db);
    verdict("tt-005", "at(t) monotone and internally consistent incl. mid-sweep instants");
}

/// Child-process entry for the fork kill -9 battery. No-op unless armed.
/// Each fork + first-op group is one transaction: a kill at any instant
/// must leave either the whole group or none of it (the consistency the
/// parent asserts after recovery).
#[test]
fn tt_006_crash_child_fork_writer() {
    let Ok(db) = std::env::var("FS_LEDGER_TRAVEL_CRASH_DB") else {
        return;
    };
    let l = Ledger::open(&db).expect("child open");
    let mut i = 0u64;
    loop {
        l.begin().expect("child txn");
        let branch = l.fork(&format!("storm-{i}"), MAIN_BRANCH).expect("fork");
        unit(&l, branch, ExecMode::Deterministic, i, i as i64);
        l.commit().expect("child commit");
        i += 1;
    }
}

#[test]
fn tt_006_crash_kill9_during_fork_traffic() {
    let exe = std::env::current_exe().expect("current exe");
    let seed: u64 = 0x5EED_F04B_0000_0006;
    let mut x = seed;
    let mut lcg = move || {
        x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        x
    };
    for round in 0..4 {
        let db = temp_db(&format!("crashfork{round}"));
        drop(Ledger::open(&db).expect("pre-create"));
        let kill_ms = 30 + (lcg() % 220);
        let mut child = std::process::Command::new(&exe)
            .args(["--exact", "tt_006_crash_child_fork_writer", "--nocapture"])
            .env("FS_LEDGER_TRAVEL_CRASH_DB", &db)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()
            .expect("spawn crash child");
        std::thread::sleep(std::time::Duration::from_millis(kill_ms));
        child.kill().expect("SIGKILL child");
        let _ = child.wait();

        let l = Ledger::open(&db).expect("recovery open");
        let lint = l.lint().unwrap();
        assert!(lint.is_clean(), "round {round}: post-recovery lint dirty: {lint:?}");
        // Every surviving branch is internally consistent: chain walkable,
        // visible ops complete, replayable views.
        let branches = l.branches().unwrap();
        for b in &branches {
            let ids = l.visible_op_ids(b.id, None).unwrap();
            for id in ids {
                let op = l.op(id).unwrap().expect("visible op exists");
                assert_eq!(op.outcome.as_deref(), Some("ok"), "partial op survived crash");
            }
        }
        println!(
            "{{\"suite\":\"fs-ledger/travel\",\"case\":\"tt-006\",\"round\":{round},\
             \"seed\":\"{seed:#x}\",\"kill_ms\":{kill_ms},\"branches\":{},\
             \"lint_clean\":true}}",
            branches.len()
        );
        drop(l);
        cleanup_db(&db);
    }
    verdict("tt-006", "4 kill -9 rounds during fork traffic recovered lint-clean");
}
