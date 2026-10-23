//! Actual cooling child processes, not a substitute numerical callback.
#[path = "../src/json_read.rs"]
mod json;
#[path = "cooling_parallel/qmc.rs"]
mod qmc;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-parallel-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path, workers: Option<&str>) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request);
    if let Some(n) = workers { c.args(["--workers", n]); }
    c
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let j = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(j.str_field("schema"), Some("frankensim.cooling-network-uq.parallel.progress.v1"));
    for name in ["mean_k", "sampling_standard_error_k", "empirical_probability_of_compliance", "percentiles_p05_p50_p95_k"] {
        assert!(j.get(name).is_none());
    }
    j
}
fn number(j: &J, name: &str) -> f64 { j.get(name).unwrap().as_f64().unwrap() }

#[test]
fn actual_steady_and_transient_cooling_are_identical_across_worker_counts() {
    for (base, request) in [("fan-correlated-hotspot.json", "uq-fan-hotspot.json"),
        ("nonlinear-contact-pulse.json", "uq-nonlinear-contact-pulse.json")]
    {
        let base = example(base); let request = example(request);
        let serial = command(&base, &request, None).output().unwrap(); parsed(&serial);
        for workers in ["1", "2", "3"] {
            let parallel = command(&base, &request, Some(workers)).output().unwrap(); parsed(&parallel);
            assert_eq!(parallel.stdout, serial.stdout, "workers={workers}");
        }
    }
}

#[test]
fn changed_worker_count_resumes_without_refunding_or_repeating_physics() {
    let dir = scratch("resume");
    let base = example("fan-correlated-hotspot.json"); let request = example("uq-fan-hotspot.json");
    let full_path = dir.join("full.batch");
    let full = command(&base, &request, Some("4")).arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first_path = dir.join("first.batch");
    let first = command(&base, &request, Some("3")).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    assert_eq!(number(&progress(&first), "samples_evaluated"), 3.0);
    let original = fs::read(&first_path).unwrap(); assert_eq!(&original[..8], b"FSUQBT01");
    let second_path = dir.join("second.batch");
    let second = command(&base, &request, Some("2")).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "2"]).output().unwrap();
    let p = progress(&second);
    assert_eq!(number(&p, "samples_evaluated"), 5.0);
    assert_eq!(number(&p, "samples_evaluated_this_run"), 2.0);
    assert_eq!(fs::read(&first_path).unwrap(), original);
    let final_path = dir.join("final.batch");
    let final_run = command(&base, &request, Some("1")).arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&final_path).output().unwrap();
    parsed(&final_run);
    assert_eq!(final_run.stdout, full.stdout);
    assert_eq!(fs::read(&final_path).unwrap(), fs::read(&full_path).unwrap());
    let short = dir.join("short.json");
    fs::write(&short, fs::read_to_string(&request).unwrap().replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let complete = command(&base, &short, Some("2")).arg("--resume").arg(&final_path).output().unwrap();
    parsed(&complete); assert_eq!(complete.stdout, full.stdout);
}

#[test]
fn invalid_workers_changed_models_and_existing_files_refuse_before_output() {
    let dir = scratch("admission"); let missing = dir.join("missing.json");
    for n in ["0", "65", "-1", "NaN"] {
        let out = command(&missing, &missing, Some(n)).output().unwrap();
        assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    }
    let base = example("fan-correlated-hotspot.json"); let request = example("uq-fan-hotspot.json");
    let saved = dir.join("empty.batch");
    let out = command(&base, &request, Some("2")).arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(number(&progress(&out), "samples_evaluated"), 0.0);
    let bytes = fs::read(&saved).unwrap();
    let changed = dir.join("changed.json"); fs::write(&changed, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    let unused = dir.join("must-not-exist.batch");
    let out = command(&changed, &request, Some("2")).arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty()); assert!(!unused.exists());
    let out = command(&base, &request, Some("2")).arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&saved).output().unwrap();
    assert!(!out.status.success()); assert_eq!(fs::read(&saved).unwrap(), bytes);
    let out = command(&missing, &missing, Some("2")).args(["--sensitivity", "sobol"])
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::USAGE))); assert!(!unused.exists());
    let help = Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json", "cooling-network-uq", "--help"]).output().unwrap();
    assert!(parsed(&help).str_field("help").unwrap().contains("--workers"));
}

#[test]
fn deadlines_retain_progress_but_physical_failures_cannot_resume_past_the_error() {
    let dir = scratch("failures"); let base = example("fan-correlated-hotspot.json");
    let request = example("uq-fan-hotspot.json"); let short = dir.join("short.json");
    fs::write(&short, fs::read_to_string(&request).unwrap().replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let saved = dir.join("timeout.batch");
    let out = command(&base, &short, Some("2")).arg("--checkpoint").arg(&saved).output().unwrap();
    let p = progress(&out); assert_eq!(p.str_field("termination"), Some("wall-time-budget"));
    let resumed = command(&base, &request, Some("3")).arg("--resume").arg(&saved).output().unwrap();
    assert_eq!(number(&parsed(&resumed), "samples_evaluated"), 8.0);
    let broken = dir.join("broken.json");
    fs::write(&broken, fs::read_to_string(&base).unwrap().replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.batch");
    let out = command(&broken, &request, Some("3")).arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(out.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let next = dir.join("never.batch");
    let out = command(&broken, &request, Some("1")).arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!out.status.success()); assert!(!next.exists());
}
