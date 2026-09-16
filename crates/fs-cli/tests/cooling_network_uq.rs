use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}

fn uq(base: &Path, request: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.arg("--json").arg("cooling-network-uq").arg(base).arg(request);
    command
}

fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-uq-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}

fn assert_complete(output: &Output) {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    assert!(String::from_utf8_lossy(&output.stdout).contains("\"status\":\"complete\""));
}

fn assert_progress(output: &Output, count: usize) {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("frankensim.cooling-network-uq.progress.v1"));
    assert!(text.contains(&format!("\"samples_evaluated\":{count},")));
    assert!(!text.contains("mean_k"));
    assert!(!text.contains("empirical_probability_of_compliance"));
    assert!(!text.contains("\"status\":\"complete\""));
}

#[test]
fn actual_binary_uq_runs_real_cooling_samples() {
    let output = uq(&example("fan-correlated-hotspot.json"), &example("uq-fan-hotspot.json"))
        .output().expect("run UQ binary");
    assert_complete(&output);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("frankensim.cooling-network-uq.result.v1"));
    assert!(stdout.contains("\"samples_evaluated\":8"));
    assert!(stdout.contains("empirical_probability_of_compliance"));
}

#[test]
fn unknown_multivariate_dependence_refuses_before_sampling() {
    let output = uq(&example("fan-correlated-hotspot.json"), &example("uq-unknown-dependence.json"))
        .output().expect("run UQ negative case");
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("explicit dependence") || stderr.contains("joint probability"));
}

#[test]
fn actual_cooling_chunks_resume_to_identical_output_and_checkpoint_bytes() {
    let dir = scratch("chunks");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let full_path = dir.join("full.bin");
    let full = uq(&base, &request).arg("--checkpoint").arg(&full_path).output().unwrap();
    assert_complete(&full);

    let first_path = dir.join("first.bin");
    let first = uq(&base, &request).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    assert_progress(&first, 3);
    let retained = fs::read(&first_path).unwrap();
    let second_path = dir.join("second.bin");
    let second = uq(&base, &request).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path)
        .args(["--max-new-samples", "2"]).output().unwrap();
    assert_progress(&second, 5);
    assert!(String::from_utf8_lossy(&second.stdout).contains("\"samples_evaluated_this_run\":2"));
    assert_eq!(fs::read(&first_path).unwrap(), retained, "resume input must remain unchanged");

    let final_path = dir.join("final.bin");
    let resumed = uq(&base, &request).arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&final_path).output().unwrap();
    assert_complete(&resumed);
    assert_eq!(resumed.stdout, full.stdout);
    assert_eq!(fs::read(&final_path).unwrap(), fs::read(&full_path).unwrap());
    let terminal = uq(&base, &request).arg("--resume").arg(&final_path).output().unwrap();
    assert_complete(&terminal);
    assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn expired_invocation_keeps_its_completed_prefix_and_original_sample_budget() {
    let dir = scratch("timeout");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let first_path = dir.join("first.bin");
    let first = uq(&base, &request).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    assert_progress(&first, 3);
    let text = fs::read_to_string(&request).unwrap();
    assert!(text.contains("\"wall_seconds\": 300"));
    let short_request = dir.join("short-time.json");
    fs::write(&short_request, text.replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let paused_path = dir.join("paused.bin");
    let paused = uq(&base, &short_request).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&paused_path).output().unwrap();
    assert_progress(&paused, 3);
    assert!(String::from_utf8_lossy(&paused.stdout).contains("\"termination\":\"wall-time-budget\""));
    let final_result = uq(&base, &request).arg("--resume").arg(&paused_path).output().unwrap();
    assert_complete(&final_result);
    let uninterrupted = uq(&base, &request).output().unwrap();
    assert_complete(&uninterrupted);
    assert_eq!(final_result.stdout, uninterrupted.stdout);
}

#[test]
fn resume_binds_the_exact_base_request_and_seed_before_reserving_an_output() {
    let dir = scratch("identity");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let prefix = dir.join("empty.bin");
    let initial = uq(&base, &request).arg("--checkpoint").arg(&prefix)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_progress(&initial, 0);
    let original = fs::read(&prefix).unwrap();
    let changed_base = dir.join("changed-base.json");
    // Exact bytes are deliberately bound, even for a whitespace-only change.
    fs::write(&changed_base, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    let changed_plan = dir.join("changed-plan.json");
    let text = fs::read_to_string(&request).unwrap();
    assert!(text.contains("\"seed\": \"73\""));
    fs::write(&changed_plan, text.replace("\"seed\": \"73\"", "\"seed\": \"74\"")).unwrap();
    for (index, (base, request)) in [(&changed_base, &request), (&base, &changed_plan)].into_iter().enumerate() {
        let destination = dir.join(format!("must-not-exist-{index}.bin"));
        let output = uq(base, request).arg("--resume").arg(&prefix)
            .arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("identity differs"));
        assert!(!destination.exists());
    }
    assert_eq!(fs::read(prefix).unwrap(), original);
}

#[test]
fn existing_checkpoints_and_corrupt_inputs_are_never_silently_replaced() {
    let dir = scratch("files");
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json");
    let path = dir.join("keep.bin");
    fs::write(&path, b"keep-this-user-file").unwrap();
    let output = uq(&base, &request).arg("--checkpoint").arg(&path)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(fs::read(&path).unwrap(), b"keep-this-user-file");
    let output = uq(&base, &request).arg("--resume").arg(&path).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cooling-network-uq-checkpoint"));
}

#[test]
fn a_refused_cooling_sample_invalidates_the_current_resumable_output() {
    let dir = scratch("refusal");
    let source = fs::read_to_string(example("fan-correlated-hotspot.json")).unwrap();
    assert!(source.contains("\"linear_iterations\": 20000"));
    let base = dir.join("bad-budget.json");
    // UQ admits the parameter plan, but the real cooling child must reject
    // this solver budget; it must not be reclassified as a resumable timeout.
    fs::write(&base, source.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let request = example("uq-fan-hotspot.json");
    let path = dir.join("failed.bin");
    let result = uq(&base, &request).arg("--checkpoint").arg(&path).output().unwrap();
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    let bytes = fs::read(&path).unwrap();
    assert!(bytes.starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let result = uq(&base, &request).arg("--resume").arg(&path).output().unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("cooling-network-uq-checkpoint"));
}

fn compliance_request(dir: &Path, name: &str, samples: usize, ceiling: Option<f64>, seconds: f64) -> PathBuf {
    let limit = ceiling.map_or_else(String::new, |value| format!(",\"temperature_limit_k\":{value}"));
    let path = dir.join(name);
    // Degenerate uncertainty makes the real child objective identical across
    // ordinals. The stopping calculation still has to earn its confidence.
    fs::write(&path, format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.v1\",\"seed\":\"73\",\"samples\":{samples},\"wall_seconds\":{seconds}{limit},\"correlation\":{{\"kind\":\"independent\"}},\"parameters\":[{{\"target\":{{\"kind\":\"air-density\"}},\"distribution\":{{\"kind\":\"uniform\",\"lo\":1.2,\"hi\":1.2}}}}]}}"
    )).unwrap();
    path
}

fn sequential(request: &Path, probability: &str, alpha: &str, minimum: &str) -> Command {
    let mut command = uq(&example("fan-correlated-hotspot.json"), request);
    command.args(["--compliance-probability", probability, "--confidence-alpha", alpha, "--min-decision-samples", minimum]);
    command
}

fn assert_inconclusive(output: &Output, count: usize, termination: &str) {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("frankensim.cooling-network-uq.compliance.v1"));
    assert!(text.contains("\"status\":\"inconclusive\""));
    assert!(text.contains("\"decision\":\"indeterminate\""));
    assert!(text.contains(&format!("\"samples_evaluated\":{count},")));
    assert!(text.contains(&format!("\"termination\":\"{termination}\"")));
    assert!(!text.contains("\"status\":\"complete\""));
}

#[test]
fn real_cooling_confidence_decides_both_sides_before_the_sample_cap() {
    let dir = scratch("sequential-decisions");
    for (index, (ceiling, decision)) in [
        (1.0e6, "meets-probability-target"), (1.0, "below-probability-target"),
    ].into_iter().enumerate() {
        let request = compliance_request(&dir, &format!("request-{index}.json"), 64, Some(ceiling), 300.0);
        let output = sequential(&request, "0.5", "0.05", "16").output().unwrap();
        assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
        let text = String::from_utf8_lossy(&output.stdout);
        assert!(text.contains("\"status\":\"decision-reached\""));
        assert!(text.contains(&format!("\"decision\":\"{decision}\"")));
        assert!(text.contains("\"samples_evaluated\":16,"));
        assert!(text.contains("\"samples_planned\":64,"));
        assert!(text.contains("\"scope\":\"probability-of-declared-numerical-model\""));
        assert!(!text.contains("\"status\":\"complete\""));
        assert!(!text.contains("sampling_standard_error_k"));
    }
}

#[test]
fn sequential_chunks_and_terminal_resume_reproduce_the_first_decision_bytes() {
    let dir = scratch("sequential-replay");
    let request = compliance_request(&dir, "request.json", 64, Some(1.0e6), 300.0);
    let full_path = dir.join("full.bin");
    let full = sequential(&request, "0.5", "0.05", "16")
        .arg("--checkpoint").arg(&full_path).output().unwrap();
    assert!(full.status.success(), "{}", String::from_utf8_lossy(&full.stderr));
    let first_path = dir.join("first.bin");
    let first = sequential(&request, "0.5", "0.05", "16")
        .arg("--checkpoint").arg(&first_path).args(["--max-new-samples", "5"]).output().unwrap();
    assert_inconclusive(&first, 5, "sample-chunk");
    let second_path = dir.join("second.bin");
    let second = sequential(&request, "0.5", "0.05", "16")
        .arg("--resume").arg(&first_path).arg("--checkpoint").arg(&second_path)
        .args(["--max-new-samples", "4"]).output().unwrap();
    assert_inconclusive(&second, 9, "sample-chunk");
    let final_path = dir.join("final.bin");
    let resumed = sequential(&request, "0.5", "0.05", "16")
        .arg("--resume").arg(&second_path).arg("--checkpoint").arg(&final_path).output().unwrap();
    assert!(resumed.status.success(), "{}", String::from_utf8_lossy(&resumed.stderr));
    assert_eq!(full.stdout, resumed.stdout);
    assert_eq!(fs::read(&full_path).unwrap(), fs::read(&final_path).unwrap());
    // No time for another real solve: a stopped checkpoint must already be
    // terminal under its unchanged policy, before launching any child.
    let short = compliance_request(&dir, "short.json", 64, Some(1.0e6), 1.0e-9);
    let terminal = sequential(&short, "0.5", "0.05", "16")
        .arg("--resume").arg(&final_path).output().unwrap();
    assert!(terminal.status.success(), "{}", String::from_utf8_lossy(&terminal.stderr));
    assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn sequential_sample_and_time_budgets_cannot_manufacture_a_decision() {
    let dir = scratch("sequential-inconclusive");
    let small = compliance_request(&dir, "small.json", 2, Some(1.0e6), 300.0);
    let output = sequential(&small, "0.99", "0.05", "2").output().unwrap();
    assert_inconclusive(&output, 2, "sample-budget");
    let short = compliance_request(&dir, "short.json", 64, Some(1.0e6), 1.0e-9);
    let output = sequential(&short, "0.5", "0.05", "16").output().unwrap();
    assert_inconclusive(&output, 0, "wall-time-budget");
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("\"probability_confidence_sequence\":null"));
    assert!(text.contains("\"empirical_probability_of_compliance\":null"));
}

#[test]
fn changed_or_removed_sequential_policy_refuses_before_output_reservation() {
    let dir = scratch("sequential-policy-binding");
    let request = compliance_request(&dir, "request.json", 64, Some(1.0e6), 300.0);
    let prefix = dir.join("prefix.bin");
    let output = sequential(&request, "0.5", "0.05", "16")
        .arg("--checkpoint").arg(&prefix).args(["--max-new-samples", "0"]).output().unwrap();
    assert_inconclusive(&output, 0, "sample-chunk");
    let retained = fs::read(&prefix).unwrap();
    let mut commands = vec![
        sequential(&request, "0.6", "0.05", "16"),
        sequential(&request, "0.5", "0.01", "16"),
        sequential(&request, "0.5", "0.05", "17"),
        uq(&example("fan-correlated-hotspot.json"), &request),
    ];
    for (index, command) in commands.iter_mut().enumerate() {
        let destination = dir.join(format!("must-not-exist-{index}.bin"));
        let output = command.arg("--resume").arg(&prefix).arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("identity differs"));
        assert!(!destination.exists());
    }
    assert_eq!(fs::read(prefix).unwrap(), retained);
}

#[test]
fn sequential_admission_needs_a_ceiling_and_an_achievable_minimum_sample_count() {
    let dir = scratch("sequential-admission");
    for (index, (samples, ceiling)) in [(64, None), (8, Some(1.0e6))].into_iter().enumerate() {
        let request = compliance_request(&dir, &format!("request-{index}.json"), samples, ceiling, 300.0);
        let destination = dir.join(format!("must-not-exist-{index}.bin"));
        let output = sequential(&request, "0.5", "0.05", "16")
            .arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!destination.exists());
    }
}

#[test]
fn sequential_model_refusal_does_not_publish_confidence_from_a_filtered_run() {
    let dir = scratch("sequential-model-refusal");
    let source = fs::read_to_string(example("fan-correlated-hotspot.json")).unwrap();
    assert!(source.contains("\"linear_iterations\": 20000"));
    let base = dir.join("bad-budget.json");
    fs::write(&base, source.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let request = compliance_request(&dir, "request.json", 64, Some(1.0e6), 300.0);
    let checkpoint = dir.join("failed.bin");
    let output = uq(&base, &request).args([
        "--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16",
    ]).arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(fs::read(&checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
}
