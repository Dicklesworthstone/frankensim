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
    let result = command(&missing, &missing, "2").arg("--checkpoint").arg(&checkpoint).args([
        "--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16",
    ]).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(String::from_utf8_lossy(&result.stderr).contains("QMC"));
    assert!(result.stdout.is_empty()); assert!(!checkpoint.exists());
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

fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(result.str_field("schema"), Some("frankensim.cooling-network-uq.qmc.progress.v1"));
    for field in ["mean_k", "estimated_probability_of_compliance", "between_replicate_standard_error_k", "replicate_means_k"] {
        assert!(result.get(field).is_none(), "partial fixed-layout report published {field}");
    }
    result
}

#[test]
fn checkpointed_real_cooling_qmc_replays_across_unfinished_nets() {
    let dir = scratch("resume");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let full_path = dir.join("full.qmc");
    let full = command(&base, &request, "2").arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first_path = dir.join("first.qmc");
    let first = command(&base, &request, "2").arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let row = progress(&first);
    assert_eq!(number(&row, "samples_evaluated"), 3.0);
    assert_eq!(number(&row, "completed_replicates"), 0.0);
    assert_eq!(number(&row, "next_replicate_ordinal"), 0.0);
    assert_eq!(number(&row, "next_point_ordinal"), 3.0);
    let retained = fs::read(&first_path).unwrap();
    let second_path = dir.join("second.qmc");
    let second = command(&base, &request, "2").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "2"]).output().unwrap();
    let row = progress(&second);
    assert_eq!(number(&row, "samples_evaluated"), 5.0);
    assert_eq!(number(&row, "samples_evaluated_this_run"), 2.0);
    assert_eq!(number(&row, "completed_replicates"), 1.0);
    assert_eq!(number(&row, "next_replicate_ordinal"), 1.0);
    assert_eq!(number(&row, "next_point_ordinal"), 1.0);
    assert_eq!(fs::read(&first_path).unwrap(), retained);
    let completed_path = dir.join("completed.qmc");
    let completed = command(&base, &request, "2").arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&completed_path).output().unwrap();
    parsed(&completed);
    assert_eq!(completed.stdout, full.stdout);
    assert_eq!(fs::read(&completed_path).unwrap(), fs::read(&full_path).unwrap());
    let fast_request = dir.join("fast.json");
    fs::write(&fast_request, fs::read_to_string(&request).unwrap()
        .replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let terminal = command(&base, &fast_request, "2").arg("--resume").arg(&completed_path).output().unwrap();
    parsed(&terminal);
    assert_eq!(terminal.stdout, full.stdout, "a complete checkpoint must not start another solve");
    // Same total points, different grouping is a different integration rule.
    let changed_path = dir.join("changed-layout.qmc");
    let changed = command(&base, &request, "4").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&changed_path).output().unwrap();
    assert!(!changed.status.success()); assert!(changed.stdout.is_empty());
    assert!(!changed_path.exists());
    assert!(String::from_utf8_lossy(&changed.stderr).contains("identity differs"));
}

#[test]
fn qmc_timeout_retains_progress_and_solver_failure_invalidates_recovery() {
    let dir = scratch("recovery-failure");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let short = dir.join("short.json");
    fs::write(&short, fs::read_to_string(&request).unwrap()
        .replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let saved = dir.join("timeout.qmc");
    let output = command(&base, &short, "2").arg("--checkpoint").arg(&saved).output().unwrap();
    let row = progress(&output);
    assert_eq!(row.str_field("termination"), Some("wall-time-budget"));
    assert_eq!(number(&row, "samples_evaluated"), 0.0);
    assert!(saved.exists());
    let resumed = command(&base, &request, "2").arg("--resume").arg(&saved).output().unwrap();
    assert_eq!(number(&parsed(&resumed), "samples_evaluated"), 8.0);
    let broken = dir.join("broken-solver.json");
    fs::write(&broken, fs::read_to_string(&base).unwrap()
        .replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.qmc");
    let output = command(&broken, &request, "2").arg("--checkpoint").arg(&failed).output().unwrap();
    assert!(!output.status.success()); assert!(output.stdout.is_empty());
    assert!(fs::read_to_string(&failed).unwrap().starts_with("FRANKENSIM-UQ-FAILED"));
    let next = dir.join("must-not-exist.qmc");
    let retry = command(&broken, &request, "2").arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!retry.status.success()); assert!(!next.exists());
}
