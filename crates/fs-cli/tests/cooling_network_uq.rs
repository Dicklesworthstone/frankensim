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
