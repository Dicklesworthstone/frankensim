//! G0/G3/e2e battery for the measured wedge decision audit (bead
//! `frankensim-extreal-program-f85xj.1.5`).

use fs_govern::{
    WEDGE_DECISION_AUDIT_SCHEMA, WedgeDecisionAuditError, WedgeDecisionAuditRequest,
    build_wedge_decision_audit,
};
use fs_wedge::KillVerdict;
use std::process::{Command, Output};

const EVIDENCE: &str = "ledger://runs/reference-cooling/cycle-time-v1";

fn run_binary(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_wedge-audit"))
        .args(arguments)
        .output()
        .expect("wedge-audit binary executes")
}

#[test]
fn complete_artifact_and_derivation_logs_are_byte_deterministic() {
    let request = WedgeDecisionAuditRequest::new(2.0, EVIDENCE).expect("valid request");
    let first = build_wedge_decision_audit(&request).expect("audit renders");
    let second = build_wedge_decision_audit(&request).expect("audit replays");

    assert_eq!(first, second);
    assert_eq!(first.kill_verdict(), KillVerdict::Indeterminate);
    assert!(first.artifact().contains(WEDGE_DECISION_AUDIT_SCHEMA));
    assert!(first.artifact().contains(EVIDENCE));
    assert!(first.artifact().contains("thermal-design-assurance"));
    assert!(
        first
            .artifact()
            .contains("frankensim-vertical-ratification-v1")
    );
    assert!(first.artifact().contains("comparison_report_tsv"));
    assert!(first.artifact().contains("\"verdict\":\"indeterminate\""));
    assert!(first.artifact().contains("commercial-decision-audit-only"));
    assert!(
        first
            .artifact()
            .contains("caller-supplied cycle-time locator")
    );
    assert_eq!(first.logs().len(), 7);
    for step in [
        "request-validated",
        "wedge-self-audit",
        "measured-inputs-and-baseline",
        "scoring-and-sensitivity",
        "ratification-validated",
        "kill-criterion-evaluated",
        "artifact-rendered",
    ] {
        assert!(
            first.logs().iter().any(|log| log.to_json().contains(step)),
            "missing log step {step}"
        );
    }
    let final_log = first.logs().last().expect("artifact log").to_json();
    assert!(final_log.contains(&first.identity().to_string()));
    assert!(final_log.contains("thermal-design-assurance"));
}

#[test]
fn missing_cycle_time_evidence_is_a_typed_actionable_refusal() {
    let refusal = WedgeDecisionAuditRequest::new(2.0, " ");
    assert_eq!(
        refusal,
        Err(WedgeDecisionAuditError::MissingCycleTimeEvidence)
    );
    let error = refusal.expect_err("empty evidence must refuse");
    assert_eq!(error.code(), "missing-cycle-time-evidence");
    assert!(error.fix().contains("--cycle-time-evidence"));
}

#[test]
fn binary_replays_artifact_and_logs_exactly() {
    let arguments = ["--measured-days", "2.0", "--cycle-time-evidence", EVIDENCE];
    let first = run_binary(&arguments);
    let second = run_binary(&arguments);
    assert!(
        first.status.success(),
        "stderr={}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(second.status.success());
    assert_eq!(first.stdout, second.stdout);
    assert_eq!(first.stderr, second.stderr);

    let stdout = String::from_utf8(first.stdout).expect("artifact is UTF-8");
    let stderr = String::from_utf8(first.stderr).expect("logs are UTF-8");
    assert!(stdout.contains(WEDGE_DECISION_AUDIT_SCHEMA));
    assert!(stdout.contains("\"identity\":"));
    assert!(stderr.lines().all(|line| line.starts_with("{\"level\":")));
    assert_eq!(stderr.lines().count(), 7);
    assert!(stderr.contains("\"level\":\"info\""));
    assert!(stderr.contains("scoring-and-sensitivity"));
    assert!(stderr.contains("thermal-design-assurance"));
    assert!(stderr.contains("frankensim-vertical-ratification-v1"));
    assert!(stderr.contains("indeterminate"));
}

#[test]
fn seeded_missing_evidence_fault_fails_with_one_clear_warning() {
    let output = run_binary(&[
        "--measured-days",
        "2.0",
        "--cycle-time-evidence",
        EVIDENCE,
        "--seed-fault",
        "missing-cycle-time-evidence",
    ]);
    assert_eq!(output.status.code(), Some(3));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).expect("warning is UTF-8");
    assert_eq!(stderr.lines().count(), 1);
    assert!(stderr.contains("\"level\":\"warn\""));
    assert!(stderr.contains("\"code\":\"missing-cycle-time-evidence\""));
    assert!(stderr.contains("--cycle-time-evidence"));
    assert!(stderr.contains("retained run or timing receipt locator"));
}

#[test]
fn cli_argument_refusal_is_stable_and_does_not_publish_an_artifact() {
    let output = run_binary(&[
        "--measured-days",
        "2.0",
        "--measured-days",
        "1.0",
        "--cycle-time-evidence",
        EVIDENCE,
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).expect("warning is UTF-8");
    assert!(stderr.contains("\"code\":\"duplicate-flag\""));
    assert!(stderr.contains("usage: wedge-audit"));
}
