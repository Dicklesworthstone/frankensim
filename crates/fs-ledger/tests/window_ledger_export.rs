//! Persistence tests for exporting coupled-driver window ledgers into
//! `fs-ledger` (`record_window_ledger` + `metrics_for_op`).

use fs_ledger::{EdgeRole, FiveExplicits, Ledger, LedgerError, OpOutcome};

fn sample_explicits() -> FiveExplicits<'static> {
    FiveExplicits {
        seed: b"window-ledger-seed-56",
        versions: r#"{"kernel":"0.1.0"}"#,
        budget: r#"{"max_steps":1000000}"#,
        capability: r#"{"mode":"deterministic"}"#,
    }
}

#[test]
fn window_ledger_artifact_and_metrics_round_trip_and_reject_invalid_inputs() {
    let ledger = Ledger::open(":memory:").expect("in-memory ledger opens");
    let op = ledger
        .begin_op(None, r#"{"kind":"coupled_run"}"#, &sample_explicits(), 100)
        .expect("begin_op succeeds");

    let canonical_json = concat!(
        r#"{"schema":"frankensim.window_ledger.v1","steps":1000000,"duration_s":1.0e3,"#,
        r#""overall_color":{"color":"verified","payload":{"interval":[-1.0e-12,1.0e-12]}},"#,
        r#""totals":{"energy_residual_j":0.0e0,"dissipated_j":0.0e0,"numerical_loss_j":0.0e0}}"#
    );
    let metrics = vec![
        ("window.steps".to_owned(), 1_000_000.0),
        ("window.energy.residual_j".to_owned(), 0.0),
        ("window.energy.numerical_loss_j".to_owned(), 0.0),
        ("window.entropy.generated_j_per_k".to_owned(), 0.0),
    ];

    let receipt = ledger
        .record_window_ledger(op, 1000, canonical_json, &metrics)
        .expect("record_window_ledger persists artifact + metrics");
    assert!(!receipt.deduped);

    let stored_bytes = ledger
        .get_artifact_bounded(&receipt.hash, 1024 * 1024)
        .expect("get_artifact_bounded succeeds")
        .expect("artifact present");
    assert_eq!(stored_bytes, canonical_json.as_bytes());
    let info = ledger
        .artifact_info(&receipt.hash)
        .expect("artifact_info succeeds")
        .expect("artifact info present");
    assert_eq!(info.kind, Ledger::WINDOW_LEDGER_ARTIFACT_KIND);

    let edges = ledger
        .op_artifact_edges_bounded(op, 16)
        .expect("op_artifact_edges_bounded succeeds");
    assert_eq!(edges.edges.len(), 1);
    assert_eq!(edges.edges[0].role, EdgeRole::Out);
    assert_eq!(edges.edges[0].artifact, receipt.hash);

    let rows = ledger.metrics_for_op(op).expect("metrics_for_op succeeds");
    assert_eq!(rows.len(), 4);
    assert_eq!(
        rows[0],
        (1000, "window.energy.numerical_loss_j".to_owned(), 0.0)
    );
    assert_eq!(rows[1], (1000, "window.energy.residual_j".to_owned(), 0.0));
    assert_eq!(
        rows[2],
        (1000, "window.entropy.generated_j_per_k".to_owned(), 0.0)
    );
    assert_eq!(rows[3], (1000, "window.steps".to_owned(), 1_000_000.0));

    // Malformed JSON and non-finite metrics are refused cleanly with zero partial writes.
    let op2 = ledger
        .begin_op(
            None,
            r#"{"kind":"coupled_run_bad"}"#,
            &sample_explicits(),
            200,
        )
        .expect("begin_op 2 succeeds");
    assert!(matches!(
        ledger.record_window_ledger(op2, 2000, "{not-json", &metrics),
        Err(LedgerError::Invalid { .. })
    ));
    assert!(matches!(
        ledger.record_window_ledger(
            op2,
            2000,
            canonical_json,
            &[("window.bad".to_owned(), f64::NAN)]
        ),
        Err(LedgerError::Invalid { .. })
    ));
    assert!(matches!(
        ledger.record_window_ledger(
            op2,
            2000,
            canonical_json,
            &[
                ("window.dup".to_owned(), 1.0),
                ("window.dup".to_owned(), 2.0),
            ]
        ),
        Err(LedgerError::Invalid { .. })
    ));
    assert_eq!(ledger.metrics_for_op(op2).unwrap(), Vec::new());

    ledger
        .finish_op(op, OpOutcome::Ok, None, 200)
        .expect("finish_op succeeds");
}

#[test]
fn execution_plan_and_causalization_audit_round_trip_and_reject_bad_json() {
    let ledger = Ledger::open(":memory:").expect("in-memory ledger opens");
    let op = ledger
        .begin_op(None, r#"{"kind":"seal_system"}"#, &sample_explicits(), 10)
        .expect("begin_op succeeds");

    let plan_json = concat!(
        r#"{"schema":"frankensim.execution_plan.v1","plan_hash":"0000000000000000000000000000000000000000000000000000000000000001","#,
        r#""block_count":1,"loop_count":1}"#
    );
    let audit_json = r#"{"schema":"frankensim.causalization_audit.v1","schema_version":1,"determination":"well_determined"}"#;

    let (plan_receipt, audit_receipt) = ledger
        .record_execution_plan(op, 10, plan_json, Some(audit_json), 1, 1, 2, 4096)
        .expect("record_execution_plan persists artifacts + metrics");
    let audit_receipt = audit_receipt.expect("audit receipt present");

    assert_eq!(
        ledger
            .artifact_info(&plan_receipt.hash)
            .unwrap()
            .unwrap()
            .kind,
        Ledger::EXECUTION_PLAN_ARTIFACT_KIND
    );
    assert_eq!(
        ledger
            .artifact_info(&audit_receipt.hash)
            .unwrap()
            .unwrap()
            .kind,
        Ledger::CAUSALIZATION_AUDIT_ARTIFACT_KIND
    );
    assert_eq!(
        ledger
            .get_artifact_bounded(&plan_receipt.hash, 1024 * 1024)
            .unwrap()
            .unwrap(),
        plan_json.as_bytes()
    );
    assert_eq!(
        ledger
            .get_artifact_bounded(&audit_receipt.hash, 1024 * 1024)
            .unwrap()
            .unwrap(),
        audit_json.as_bytes()
    );

    let edges = ledger.op_artifact_edges_bounded(op, 16).unwrap();
    assert_eq!(edges.edges.len(), 2);

    let rows = ledger.metrics_for_op(op).unwrap();
    assert_eq!(
        rows,
        vec![
            (10, "plan.blocks".to_owned(), 1.0),
            (10, "plan.kept_unknowns".to_owned(), 2.0),
            (10, "plan.loops".to_owned(), 1.0),
            (10, "plan.workspace_bytes".to_owned(), 4096.0),
        ]
    );

    let op_bad = ledger
        .begin_op(None, r#"{"kind":"seal_bad"}"#, &sample_explicits(), 20)
        .unwrap();
    assert!(matches!(
        ledger.record_execution_plan(op_bad, 20, "{bad", None, 1, 0, 0, 0),
        Err(LedgerError::Invalid { .. })
    ));
    assert!(matches!(
        ledger.record_execution_plan(op_bad, 20, plan_json, Some("{bad"), 1, 0, 0, 0),
        Err(LedgerError::Invalid { .. })
    ));
    assert_eq!(ledger.metrics_for_op(op_bad).unwrap(), Vec::new());

    ledger.finish_op(op, OpOutcome::Ok, None, 30).unwrap();
}

#[test]
fn lowering_receipt_round_trip_and_reject_bad_json() {
    let ledger = Ledger::open(":memory:").expect("in-memory ledger opens");
    let op = ledger
        .begin_op(
            None,
            r#"{"kind":"lower_assembly"}"#,
            &sample_explicits(),
            40,
        )
        .expect("begin_op succeeds");

    let receipt_json = concat!(
        r#"{"schema":"frankensim.lowering_receipt.v1","receipt_hash":"0000000000000000000000000000000000000000000000000000000000000002","#,
        r#""body":{"assembly_name":"slice_a","components":[],"interfaces":[],"joints":[]}}"#
    );

    let receipt = ledger
        .record_lowering_receipt(op, 40, receipt_json, 2, 1, 1, 3)
        .expect("record_lowering_receipt persists artifact + metrics");
    assert!(!receipt.deduped);

    assert_eq!(
        ledger.artifact_info(&receipt.hash).unwrap().unwrap().kind,
        Ledger::LOWERING_RECEIPT_ARTIFACT_KIND
    );
    assert_eq!(
        ledger
            .get_artifact_bounded(&receipt.hash, 1024 * 1024)
            .unwrap()
            .unwrap(),
        receipt_json.as_bytes()
    );

    let edges = ledger.op_artifact_edges_bounded(op, 16).unwrap();
    assert_eq!(edges.edges.len(), 1);
    assert_eq!(edges.edges[0].role, EdgeRole::Out);
    assert_eq!(edges.edges[0].artifact, receipt.hash);

    let rows = ledger.metrics_for_op(op).unwrap();
    assert_eq!(
        rows,
        vec![
            (40, "lowering.components".to_owned(), 2.0),
            (40, "lowering.interfaces".to_owned(), 1.0),
            (40, "lowering.joints".to_owned(), 1.0),
            (40, "lowering.property_receipts".to_owned(), 3.0),
        ]
    );

    let op_bad = ledger
        .begin_op(None, r#"{"kind":"lower_bad"}"#, &sample_explicits(), 50)
        .unwrap();
    assert!(matches!(
        ledger.record_lowering_receipt(op_bad, 50, "{not-valid-json", 1, 0, 0, 0),
        Err(LedgerError::Invalid { .. })
    ));
    assert_eq!(ledger.metrics_for_op(op_bad).unwrap(), Vec::new());

    ledger.finish_op(op, OpOutcome::Ok, None, 60).unwrap();
}
