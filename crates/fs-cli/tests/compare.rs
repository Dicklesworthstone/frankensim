//! Integration tests for `frankensim compare` (bead frankensim-extreal-program-f85xj.6.14.1).

use fs_cli::{exit, run};

#[test]
fn test_cli_compare_text_output() {
    let output = run(["compare".to_string(), "run_base".to_string(), "run_opt".to_string()]);
    assert_eq!(output.exit_code, exit::SUCCESS);
    assert!(output.stdout.contains("FrankenSim Semantic Run Comparison"));
    assert!(output.stdout.contains("junction_maximum"));
    assert!(output.stdout.contains("thermal_margin"));
}

#[test]
fn test_cli_compare_json_output() {
    let output = run([
        "--json".to_string(),
        "compare".to_string(),
        "run_01".to_string(),
        "run_02".to_string(),
    ]);
    assert_eq!(output.exit_code, exit::SUCCESS);
    assert!(output.stdout.contains("\"schema\": \"frankensim.cli.compare-result.v1\""));
    assert!(output.stdout.contains("\"status\": \"changed\""));
    assert!(output.stdout.contains("\"name\": \"junction_maximum\""));
}

#[test]
fn test_cli_compare_with_ledger_path() {
    let output = run([
        "--json".to_string(),
        "compare".to_string(),
        "run_01".to_string(),
        "run_02".to_string(),
        "/tmp/test_ledger.db".to_string(),
    ]);
    assert_eq!(output.exit_code, exit::SUCCESS);
    assert!(output.stdout.contains("\"command\": \"compare\""));
}
