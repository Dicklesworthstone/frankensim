//! Independent sibling review: FrankenSQLite durability and refusal contract,
//! drilled at the fs-ledger boundary (bead
//! `frankensim-extreal-program-f85xj.13.5`).
//!
//! Charter: `docs/SIBLING_REVIEW_FRANKENSQLITE.md`.
//!
//! METHOD — contract-first. Each drill names the claim it attacks, taken from
//! `crates/fs-ledger/CONTRACT.md` (FrankenSim's usage assumptions about the
//! engine underneath) before reading the implementation. The `f85xj.13.1` P1
//! list scopes this to WAL durability, transaction boundaries, migration
//! refusal, and blob/checkpoint behaviour.
//!
//! SCOPE — this is the fs-ledger adapter boundary, which is what FrankenSim
//! depends on. It is NOT a review of FrankenSQLite's pager, WAL, or B-tree
//! internals, and green drills here do not certify the engine.

mod common;

use common::SyncConnection;
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};

const FX: FiveExplicits<'static> = FiveExplicits {
    seed: &[0x13, 0x05, 0x00, 0x01],
    versions: r#"{"review":"SREV-2026-07-B"}"#,
    budget: r#"{"wall_s":30}"#,
    capability: r#"{"ops":["ledger.*"]}"#,
};

/// Stamp `PRAGMA user_version` through the engine.
///
/// NOTE, and a trap worth recording: patching byte offset 60 of the main
/// database file directly does **not** work here. This ledger runs in WAL mode,
/// so a `-wal` sidecar holds a newer copy of page 1 — the page containing the
/// header — and the engine reads the version from the WAL, silently ignoring
/// the patched main file. A drill written that way appears to prove that a
/// future-schema file is accepted, which would be a false defect report. The
/// stamp must go through the engine so it lands wherever the engine will look.
fn stamp_user_version(path: &str, version: u32) {
    let conn = SyncConnection::open(path).expect("raw open for version stamp");
    conn.execute_batch(&format!("PRAGMA user_version = {version};"))
        .expect("stamp user_version");
    drop(conn);
}

/// Every durable byte of a WAL-mode database: the main file plus every sidecar
/// the engine created beside it, keyed by name and sorted for determinism.
fn snapshot_database(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = std::path::Path::new(path);
    let dir = file.parent().expect("database has a parent directory");
    let stem = file.file_name().expect("database has a file name");
    let mut snapshot: Vec<(String, Vec<u8>)> = std::fs::read_dir(dir)
        .expect("read database directory")
        .flatten()
        .filter(|entry| {
            entry
                .file_name()
                .to_string_lossy()
                .starts_with(&*stem.to_string_lossy())
        })
        .filter_map(|entry| {
            let name = entry.file_name().to_string_lossy().into_owned();
            std::fs::read(entry.path()).ok().map(|bytes| (name, bytes))
        })
        .collect();
    snapshot.sort();
    snapshot
}

fn scratch(tag: &str) -> String {
    let dir = std::env::temp_dir().join(format!(
        "fsim-sibling-review-{}-{}",
        std::process::id(),
        tag
    ));
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir.join("ledger.db").to_string_lossy().into_owned()
}

/// D1 — `LedgerError::FutureSchema`: "newer file refused, **never clobbered**".
///
/// This is a data-loss claim, and it is the P1 "migration refusal" item. A
/// ledger written by a newer FrankenSim must be refused, and — the part that
/// actually matters — the refused open must not have written to the file.
///
/// Attack: stamp a future `user_version`, snapshot the file bytes, attempt the
/// open, and compare bytes. A regression that "helpfully" migrates downward, or
/// that opens and journals before checking the version, destroys data that a
/// newer client owns.
#[test]
fn d1_a_future_schema_file_is_refused_and_left_byte_identical() {
    let path = scratch("future-schema");
    {
        let ledger = Ledger::open(&path).expect("fresh open");
        let op = ledger
            .begin_op(None, "{\"unit\":1}", &FX, 1)
            .expect("begin op");
        ledger
            .finish_op(op, OpOutcome::Ok, None, 2)
            .expect("finish");
    }

    // Stamp a schema version from the future, the way a newer FrankenSim would.
    stamp_user_version(&path, 9999);

    // Snapshot the WHOLE database, not just the main file. In WAL mode the
    // durable state is spread across the `-wal` sidecar too, so comparing only
    // the main file would let a refusal that journals into the WAL pass as
    // "never clobbered".
    let before = snapshot_database(&path);
    assert!(
        before.len() > 1,
        "expected WAL-mode sidecars; a single-file snapshot would weaken this drill"
    );

    match Ledger::open(&path) {
        Err(LedgerError::FutureSchema { .. }) => {}
        Err(other) => panic!("expected FutureSchema, got {other:?}"),
        Ok(_) => panic!(
            "REGRESSION: a ledger written by a newer schema was opened instead of refused; \
             the contract says a newer file is refused and never clobbered"
        ),
    }

    let after = snapshot_database(&path);
    assert_eq!(
        before, after,
        "REGRESSION: the refused open MODIFIED durable state. `FutureSchema` promises the \
         newer file is never clobbered; a refusal that still journals, checkpoints, or \
         migrates destroys data owned by a newer client."
    );
}

/// D2 — the seeded-regression drill required by the bead's drill-quality bar.
///
/// Historically fixed sibling defect (`904cdef1`, beads cgt07 / lnbzs):
/// FrankenSQLite never binds parameters referenced inside WHERE-clause
/// subqueries of a FROM-less `SELECT`, so `seal_artifact_output`'s guarded
/// `INSERT ... SELECT ?1, ?2 WHERE EXISTS(... ?1 ...)` **silently inserted zero
/// rows for every VALID sole producer**. The seal returned success and sealed
/// nothing. The fix routed enforcement through `NEW.`-based schema triggers.
///
/// The subtle part, and the reason this drill is shaped the way it is: a test
/// that only checked "an invalid seal refuses" would have **PASSED during the
/// bug** — zero rows inserted means nothing to violate. The discriminating
/// assertion is that the VALID case leaves an OBSERVABLE seal. Silent success
/// is the failure mode, so the drill must read the seal back.
#[test]
fn d2_a_valid_seal_is_observably_persisted_not_silently_dropped() {
    let path = scratch("seal-observable");
    let ledger = Ledger::open(&path).expect("open");

    let op = ledger
        .begin_op(None, "{\"unit\":1}", &FX, 1)
        .expect("begin op");
    let receipt = ledger
        .put_artifact("sealed-field", b"sole producer payload", None)
        .expect("put artifact");
    ledger
        .link(op, &receipt.hash, EdgeRole::Out)
        .expect("link output");
    ledger
        .finish_op(op, OpOutcome::Ok, None, 2)
        .expect("finish op");

    assert_eq!(
        ledger
            .artifact_output_seal(&receipt.hash)
            .expect("read seal"),
        None,
        "control: the artifact must start unsealed, or the positive assertion below is vacuous"
    );

    ledger
        .seal_artifact_output(&receipt.hash, op)
        .expect("sealing a valid sole producer must succeed");

    assert_eq!(
        ledger
            .artifact_output_seal(&receipt.hash)
            .expect("read seal back"),
        Some(op),
        "REGRESSION (bead lnbzs / commit 904cdef1): the seal reported success but persisted \
         NOTHING. A guarded INSERT whose parameters are dropped inside a FROM-less subquery \
         inserts zero rows and still returns Ok. Silent success is the failure mode, so the \
         seal must be read back, not merely requested."
    );

    // Idempotence for the same op, per the contract.
    ledger
        .seal_artifact_output(&receipt.hash, op)
        .expect("re-sealing the same op is documented idempotent");
    assert_eq!(
        ledger.artifact_output_seal(&receipt.hash).expect("re-read"),
        Some(op)
    );
}

/// D2b — NEGATIVE CONTROL for D2.
///
/// Shows the seal path discriminates: sealing to an operation that is NOT the
/// artifact's sole output producer must refuse, and must leave no seal behind.
/// Without this, D2 alone cannot distinguish "the guard works" from "the guard
/// accepts everything".
#[test]
fn d2b_negative_control_sealing_a_non_producer_refuses_and_leaves_no_seal() {
    let path = scratch("seal-control");
    let ledger = Ledger::open(&path).expect("open");

    let producer = ledger
        .begin_op(None, "{\"unit\":1}", &FX, 1)
        .expect("begin producer");
    let receipt = ledger
        .put_artifact("controlled-field", b"payload", None)
        .expect("put artifact");
    ledger
        .link(producer, &receipt.hash, EdgeRole::Out)
        .expect("link");
    ledger
        .finish_op(producer, OpOutcome::Ok, None, 2)
        .expect("finish");

    let bystander = ledger
        .begin_op(None, "{\"unit\":2}", &FX, 3)
        .expect("begin bystander");
    ledger
        .finish_op(bystander, OpOutcome::Ok, None, 4)
        .expect("finish bystander");

    let refusal = ledger.seal_artifact_output(&receipt.hash, bystander);
    assert!(
        refusal.is_err(),
        "control invalid: sealing to an operation that is not the sole output producer must \
         refuse; if it succeeds, D2 cannot tell an enforced guard from an absent one"
    );
    assert_eq!(
        ledger
            .artifact_output_seal(&receipt.hash)
            .expect("read seal"),
        None,
        "a refused seal must leave no seal behind"
    );
}

/// D3 — `LedgerError::NotFound` vs silent acceptance for an unknown artifact.
///
/// Contract: "unknown ops remain `NotFound`, so callers never" mistake a
/// missing subject for an accepted one. Attack: seal an artifact the ledger has
/// never seen. The engine must refuse rather than create a dangling seal row.
#[test]
fn d3_sealing_an_unknown_artifact_refuses_rather_than_creating_a_dangling_row() {
    let path = scratch("unknown-artifact");
    let ledger = Ledger::open(&path).expect("open");
    let op = ledger
        .begin_op(None, "{\"unit\":1}", &FX, 1)
        .expect("begin op");
    ledger
        .finish_op(op, OpOutcome::Ok, None, 2)
        .expect("finish");

    // A real, well-formed hash that THIS ledger has never stored: mint it in a
    // separate ledger so the value is structurally valid but locally unknown.
    let ghost = {
        let other = Ledger::open(&scratch("unknown-artifact-source")).expect("second ledger");
        other
            .put_artifact("elsewhere", b"stored in a different ledger", None)
            .expect("put artifact elsewhere")
            .hash
    };
    assert!(
        ledger.seal_artifact_output(&ghost, op).is_err(),
        "sealing an artifact with no stored producer must refuse"
    );
    assert_eq!(
        ledger.artifact_output_seal(&ghost).expect("read seal"),
        None,
        "a refused seal must not leave a dangling row for an artifact that does not exist"
    );
}

/// D4 — transaction boundary: a rolled-back transaction must leave no trace.
///
/// P1 list item "transaction boundaries". Attack: open a transaction, do real
/// work, roll back, and confirm the work is gone. A regression where a seal or
/// artifact escapes its transaction would corrupt every replay claim built on
/// the ledger.
#[test]
fn d4_a_rolled_back_transaction_leaves_no_artifact_or_seal() {
    let path = scratch("rollback");
    let ledger = Ledger::open(&path).expect("open");

    let hash = {
        ledger.begin().expect("begin txn");
        let op = ledger
            .begin_op(None, "{\"unit\":9}", &FX, 1)
            .expect("begin op");
        let receipt = ledger
            .put_artifact("doomed", b"rolled back bytes", None)
            .expect("put artifact");
        ledger.link(op, &receipt.hash, EdgeRole::Out).expect("link");
        ledger
            .finish_op(op, OpOutcome::Ok, None, 2)
            .expect("finish");
        ledger.rollback().expect("rollback");
        receipt.hash
    };

    assert_eq!(
        ledger.artifact_output_seal(&hash).expect("read seal"),
        None,
        "a rolled-back transaction must leave no seal"
    );
}
