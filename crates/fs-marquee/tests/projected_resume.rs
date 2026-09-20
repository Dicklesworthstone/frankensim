//! Real-process restart of the original projected study, not a field warm start.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);
fn workspace() -> PathBuf {
    let root = std::env::temp_dir().join(format!("fs-projected-resume-{}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&root).unwrap();
    std::fs::write(root.join("loads.csv"),
        "right,0.375,0.625,0,-1,0.7\nright,0.375,0.625,0.5,0,0.3\n").unwrap();
    root
}
fn start(root: &Path, output: &str, budget: &str, pause: Option<&str>, stress: bool) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"));
    cmd.arg("--projected").arg(root.join(output)).arg(root.join("loads.csv"))
        .args(["3", "2", "0.6", "8", "sum", budget, "--checkpoint"]);
    if let Some(pause) = pause { cmd.args(["--pause-after", pause]); }
    if stress { cmd.args(["--stress-limit", "100000000"]); }
    cmd.output().unwrap()
}
fn resume(checkpoint: &Path, output: &Path, budget: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .args(["--projected", "--resume"]).arg(checkpoint).arg(output)
        .args(["--recovery-solves", budget]).output().unwrap()
}
fn text(path: &Path) -> String { std::fs::read_to_string(path).unwrap() }

#[test]
fn paused_stress_study_resumes_the_exact_uninterrupted_tail() {
    let root = workspace();
    let full = start(&root, "full", "64", None, true);
    assert!(matches!(full.status.code(), Some(0 | 11)), "{}", String::from_utf8_lossy(&full.stderr));
    let paused = start(&root, "paused", "64", Some("1"), true);
    assert_eq!(paused.status.code(), Some(14), "{}", String::from_utf8_lossy(&paused.stderr));
    assert_eq!(text(&root.join("paused/trajectory.jsonl")).lines().count(), 1,
        "test must pause after a real accepted step");
    let source = root.join("paused/checkpoint.fscp");
    let source_bytes = std::fs::read(&source).unwrap();
    let continued = resume(&source, &root.join("resumed"), "4");
    assert_eq!(continued.status.code(), full.status.code(), "{}", String::from_utf8_lossy(&continued.stderr));
    for name in ["baseline-level-set.csv", "level-set.csv", "load-cases.csv", "checkpoint.fscp"] {
        assert_eq!(std::fs::read(root.join("full").join(name)).unwrap(),
            std::fs::read(root.join("resumed").join(name)).unwrap(), "{name}");
    }
    let mut joined = text(&root.join("paused/trajectory.jsonl"));
    joined.push_str(&text(&root.join("resumed/trajectory.jsonl")));
    assert_eq!(joined, text(&root.join("full/trajectory.jsonl")));
    assert_eq!(std::fs::read(&source).unwrap(), source_bytes);
    let summary = text(&root.join("resumed/summary.json"));
    assert!(summary.contains("\"segment_start_iteration\":1"));
    assert!(summary.contains("\"recovery_solves_started\":4"));
    assert!(summary.contains("\"final_stress\":"));
}

#[test]
fn resume_does_not_refund_an_exhausted_original_solve_budget() {
    let root = workspace();
    assert_eq!(start(&root, "paused", "2", Some("0"), false).status.code(), Some(14));
    let source = root.join("paused/checkpoint.fscp");
    let result = resume(&source, &root.join("resumed"), "4");
    assert_eq!(result.status.code(), Some(13), "{}", String::from_utf8_lossy(&result.stderr));
    let summary = text(&root.join("resumed/summary.json"));
    assert!(summary.contains("\"solves_started\":2"));
    assert!(summary.contains("\"accepted_updates\":0"));
    assert!(summary.contains("\"recovery_solves_started\":4"));
    assert_eq!(std::fs::read(&source).unwrap(), std::fs::read(root.join("resumed/checkpoint.fscp")).unwrap());
}

#[test]
fn corrupt_or_underfunded_checkpoint_cannot_create_successful_outputs() {
    let root = workspace();
    assert_eq!(start(&root, "paused", "2", Some("0"), false).status.code(), Some(14));
    let source = root.join("paused/checkpoint.fscp");
    let original = std::fs::read(&source).unwrap();
    let mut corrupted = original.clone();
    *corrupted.last_mut().unwrap() ^= 1;
    std::fs::write(root.join("bad.fscp"), corrupted).unwrap();
    let result = resume(&root.join("bad.fscp"), &root.join("bad-output"), "4");
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("hash mismatch"));
    assert!(result.stdout.is_empty());
    assert!(!root.join("bad-output").exists());
    let result = resume(&source, &root.join("underfunded"), "3");
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("recovery_solves_started=0"));
    assert!(!root.join("underfunded").exists());
    assert_eq!(std::fs::read(&source).unwrap(), original);
}

#[test]
fn resume_refuses_policy_changes_and_existing_destinations() {
    let root = workspace();
    let result = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .args(["--projected", "--resume"]).arg(root.join("missing.fscp"))
        .arg(root.join("out")).args(["--stress-limit", "100"])
        .output().unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("immutable"));
    assert!(!root.join("out").exists());
    std::fs::create_dir(root.join("occupied")).unwrap();
    std::fs::write(root.join("occupied/keep"), b"original").unwrap();
    let result = resume(&root.join("missing.fscp"), &root.join("occupied"), "4");
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("refusing to overwrite"));
    assert_eq!(std::fs::read(root.join("occupied/keep")).unwrap(), b"original");
    assert!(result.stdout.is_empty());
}

#[path = "projected_resume/refinement.rs"]
mod refinement;
