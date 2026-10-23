//! Reuse the real-binary pilot fixture and its physical cooling producer.
use super::*;

const UNIT_LAYOUT: &str = "2,3,1,0.5,1,6";
fn constant_input(samples: usize, wall: f64, transient: bool) -> String {
    request(samples, wall, transient).replace("\"std_dev\":0.2", "\"std_dev\":0")
        .replace("\"temperature_limit_k\":303", "\"temperature_limit_k\":1000000000")
}
fn decide(base: &Path, input: &Path, minimum: &str, alpha: &str) -> Command {
    let mut c = command(base, input, UNIT_LAYOUT);
    c.args(["--compliance-probability", "0.5", "--confidence-alpha", alpha,
        "--min-decision-samples", minimum]);
    c
}
fn assessed(output: &Output, code: u8) -> J {
    assert_eq!(output.status.code(), Some(i32::from(code)), "{}", String::from_utf8_lossy(&output.stderr));
    let r = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(r.str_field("schema"), Some("frankensim.cooling-network-uq.importance-pilot-compliance.v1"));
    assert_eq!(r.str_field("confidence_scope"), Some("independent-production-only"));
    for field in ["mean_k", "percentiles_p05_p50_p95_k", "failure_probability_standard_error"] {
        assert!(r.get(field).is_none());
    }
    r
}

#[test]
fn learning_does_not_count_toward_the_decision_and_restart_preserves_the_stop() {
    let dir = scratch("decision-recovery");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, constant_input(64, 300.0, false)).unwrap();
    let full_path = dir.join("full.pilot");
    let full = decide(&base, &input, "16", "0.05").arg("--checkpoint").arg(&full_path).output().unwrap();
    let full_report = assessed(&full, fs_cli::exit::SUCCESS);
    assert_eq!(number(&full_report, "pilot_samples_evaluated"), 6.0);
    assert_eq!(number(&full_report, "production_samples_evaluated"), 16.0);
    assert_eq!(number(&full_report, "confidence_samples"), 16.0);
    assert_eq!(number(&full_report, "samples_evaluated"), 22.0);
    assert_eq!(full_report.str_field("decision"), Some("meets-probability-target"));
    let first_path = dir.join("first.pilot");
    let first = decide(&base, &input, "16", "0.05").arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "5"]).output().unwrap();
    let r = assessed(&first, fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("stage"), Some("pilot"));
    assert_eq!(number(&r, "samples_evaluated"), 5.0);
    assert_eq!(number(&r, "confidence_samples"), 0.0);
    assert_eq!(r.get("failure_probability_confidence_sequence"), Some(&J::Null));
    let source = fs::read(&first_path).unwrap();
    let second_path = dir.join("second.pilot");
    let second = decide(&base, &input, "16", "0.05").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "3"]).output().unwrap();
    let r = assessed(&second, fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("stage"), Some("production"));
    assert_eq!(number(&r, "samples_evaluated"), 8.0);
    assert_eq!(number(&r, "confidence_samples"), 2.0);
    let done_path = dir.join("done.pilot");
    let done = decide(&base, &input, "16", "0.05").arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&done_path).output().unwrap();
    assessed(&done, fs_cli::exit::SUCCESS);
    assert_eq!(full.stdout, done.stdout);
    assert_eq!(fs::read(&full_path).unwrap(), fs::read(&done_path).unwrap());
    assert_eq!(fs::read(&first_path).unwrap(), source);
    let short = dir.join("short.json"); fs::write(&short, constant_input(64, 1e-9, false)).unwrap();
    let terminal = decide(&base, &short, "16", "0.05").arg("--resume").arg(&done_path).output().unwrap();
    assessed(&terminal, fs_cli::exit::SUCCESS);
    assert_eq!(terminal.stdout, full.stdout);
    let never = dir.join("changed-alpha.pilot");
    let changed = decide(&base, &input, "16", "0.01").arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&never).output().unwrap();
    assert_eq!(changed.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(!never.exists()); assert!(changed.stdout.is_empty());
    assert!(String::from_utf8_lossy(&changed.stderr).contains("identity differs"));
}

#[test]
fn insufficient_production_is_refused_and_unresolved_production_exits_budget() {
    let dir = scratch("decision-budget");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, constant_input(12, 300.0, false)).unwrap();
    let never = dir.join("unattainable-minimum.pilot");
    let bad = decide(&base, &input, "7", "0.05").arg("--checkpoint").arg(&never).output().unwrap();
    assert_eq!(bad.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(bad.stdout.is_empty()); assert!(!never.exists());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("AFTER every declared pilot"));
    let output = command(&base, &input, UNIT_LAYOUT).args(["--compliance-probability", "0.999",
        "--confidence-alpha", "0.05", "--min-decision-samples", "2"]).output().unwrap();
    let r = assessed(&output, fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("termination"), Some("sample-budget"));
    assert_eq!(r.str_field("decision"), Some("indeterminate"));
    assert_eq!(number(&r, "samples_evaluated"), 12.0);
    assert_eq!(number(&r, "confidence_samples"), 6.0);
    assert_eq!(number(&r, "estimated_failure_probability"), 0.0);
}

#[test]
fn transient_decisions_count_completed_trajectories_and_timeout_has_no_pilot_inference() {
    let dir = scratch("decision-transient");
    let base = example("nonlinear-contact-pulse.json");
    let input = dir.join("uq.json"); fs::write(&input, constant_input(24, 300.0, true)).unwrap();
    let output = decide(&base, &input, "16", "0.05").output().unwrap();
    let r = assessed(&output, fs_cli::exit::SUCCESS);
    assert_eq!(number(&r, "samples_evaluated"), 22.0);
    assert_eq!(number(&r, "confidence_samples"), 16.0);
    assert_eq!(r.get("qoi").unwrap().str_field("observation"), Some("one-completed-trajectory"));
    fs::write(&input, constant_input(24, 1e-9, true)).unwrap();
    let r = assessed(&decide(&base, &input, "16", "0.05").output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("termination"), Some("wall-time-budget"));
    assert_eq!(number(&r, "confidence_samples"), 0.0);
    assert_eq!(r.get("failure_probability_confidence_sequence"), Some(&J::Null));
}
