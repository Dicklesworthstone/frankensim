//! Integration tests for `frankensim compare` (bead frankensim-extreal-program-f85xj.6.14.1).

use fs_cli::{exit, run};

#[test]
fn compare_text_refuses_instead_of_emitting_fixture_evidence() {
    let output = run(["compare".to_string(), "run_base".to_string(), "run_opt".to_string()]);
    assert_eq!(output.exit_code, exit::UNAVAILABLE);
    assert!(output.stdout.contains("status=unavailable"));
    assert!(output.stderr.contains("cli-stage-unavailable"));
    assert!(!output.stdout.contains("junction_maximum"));
    assert!(!output.stdout.contains("thermal_margin"));
}

#[test]
fn identical_runs_do_not_produce_a_fabricated_change() {
    let output = run([
        "--json".to_string(),
        "compare".to_string(),
        "same_run".to_string(),
        "same_run".to_string(),
    ]);
    assert_eq!(output.exit_code, exit::UNAVAILABLE);
    assert!(output.stdout.contains("\"status\":\"unavailable\""));
    assert!(output.stderr.contains("\"code\":\"cli-stage-unavailable\""));
    assert!(!output.stdout.contains("\"status\":\"changed\""));
    assert!(!output.stdout.contains("junction_maximum"));
}

#[test]
fn nonexistent_ledger_cannot_be_ignored_to_mint_comparison_authority() {
    let output = run([
        "--json".to_string(),
        "compare".to_string(),
        "same_run".to_string(),
        "same_run".to_string(),
        "/definitely/missing/frankensim-ledger.db".to_string(),
    ]);
    assert_eq!(output.exit_code, exit::UNAVAILABLE);
    assert!(output.stdout.contains("\"command\":\"compare\""));
    assert!(!output.stdout.contains("evidence-aware-semantic-run-diff"));
    assert!(!output.stdout.contains("verified"));
}
