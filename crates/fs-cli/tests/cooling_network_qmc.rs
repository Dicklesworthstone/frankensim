//! Real-binary QMC consumer checks; no substitute thermal solver.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn command(base: &Path, request: &Path, replicates: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.arg("--json").arg("cooling-network-uq").arg(base).arg(request)
        .args(["--qmc-replicates", replicates]); c
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let p = std::env::temp_dir().join(format!("fs-cooling-qmc-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&p).unwrap(); p
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn number(root: &J, name: &str) -> f64 { root.get(name).unwrap().as_f64().unwrap() }

#[test]
fn actual_cooling_qmc_replays_and_reports_between_replica_error() {
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let first = command(&base, &request, "2").output().unwrap();
    let result = parsed(&first);
    assert_eq!(result.get("schema").unwrap().as_str(), Some("frankensim.cooling-network-uq.qmc.v1"));
    assert_eq!(number(&result, "samples_evaluated"), 8.0);
    assert_eq!(number(&result, "replicates"), 2.0);
    assert_eq!(number(&result, "samples_per_replicate"), 4.0);
    let means = result.get("replicate_means_k").unwrap().as_array().unwrap();
    assert_eq!(means.len(), 2);
    let a = means[0].as_f64().unwrap(); let b = means[1].as_f64().unwrap();
    assert!((number(&result, "mean_k") - 0.5 * (a + b)).abs() < 1e-10);
    assert!((number(&result, "between_replicate_standard_error_k") - 0.5 * (a - b).abs()).abs() < 1e-10);
    assert!(result.get("sampling_standard_error_k").is_none());
    assert!(result.get("confidence_interval").is_none());
    let again = command(&base, &request, "2").output().unwrap();
    parsed(&again);
    assert_eq!(first.stdout, again.stdout);
}

#[test]
fn invalid_net_layout_and_unknown_joint_measure_refuse_without_a_result() {
    let base = example("fan-correlated-hotspot.json");
    for (request, replicates) in [("uq-fan-hotspot.json", "3"), ("uq-unknown-dependence.json", "2")] {
        let output = command(&base, &example(request), replicates).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn incompatible_policies_refuse_before_input_reads_or_file_creation() {
    let dir = scratch("admission");
    let missing = dir.join("does-not-exist.json");
    let checkpoint = dir.join("must-not-be-created.bin");
    let result = command(&missing, &missing, "2").arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(result.stdout.is_empty()); assert!(!checkpoint.exists());
    let result = command(&missing, &missing, "2").args([
        "--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16",
    ]).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(String::from_utf8_lossy(&result.stderr).contains("QMC"));
}

#[test]
fn deadline_and_physical_solver_failure_never_publish_a_partial_distribution() {
    let dir = scratch("failures");
    let base = example("fan-correlated-hotspot.json");
    let original = fs::read_to_string(example("uq-fan-hotspot.json")).unwrap();
    assert!(original.contains("\"wall_seconds\": 300"));
    let short = dir.join("short.json");
    fs::write(&short, original.replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let output = command(&base, &short, "2").output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
    let original_base = fs::read_to_string(&base).unwrap();
    assert!(original_base.contains("\"linear_iterations\": 20000"));
    let invalid = dir.join("invalid-solver.json");
    fs::write(&invalid, original_base.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let output = command(&invalid, &example("uq-fan-hotspot.json"), "2").output().unwrap();
    assert!(!output.status.success());
    assert_ne!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
}
