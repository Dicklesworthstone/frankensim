#![cfg(feature = "marquee")]
//! Execute the real stress marquee binary; no stand-in solver or CLI.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!("fs-projected-stress-cli-{}-{}",
            std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir(&path).expect("unique test scratch directory");
        Self(path)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) { let _ = std::fs::remove_dir_all(&self.0); }
}

fn run(path: &Path, limit: &str, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-stress"))
        .arg("--projected").arg(path).arg(limit).args(extra).output().unwrap()
}

fn number(text: &str, name: &str) -> f64 {
    let marker = format!("\"{name}\":");
    let value = text.split_once(&marker).unwrap().1;
    value.split([',', '}']).next().unwrap().parse().unwrap()
}

#[test]
fn real_command_accepts_same_material_step_and_retains_its_geometry_and_stress() {
    let scratch = Scratch::new();
    let output_dir = scratch.0.join("study");
    let output = run(&output_dir, "1e12", &["3", "1", "0.6", "8", "0", "300"]);
    assert!(output.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    let summary = std::fs::read_to_string(output_dir.join("summary.json")).unwrap();
    assert_eq!(String::from_utf8(output.stdout).unwrap(), summary);
    assert!(summary.contains("\"status\":\"iteration_limit\""));
    assert!(summary.contains("\"converged\":false"));
    assert!(summary.contains("\"continuous_max_stress\":false"));
    assert_eq!(number(&summary, "accepted_updates"), 1.0);
    let baseline = summary.split_once("\"baseline\":").unwrap().1;
    let current = summary.split_once("\"current\":").unwrap().1;
    assert!(number(current, "compliance") < number(baseline, "compliance"));
    assert!((number(baseline, "volume") - 0.6).abs() <= 1e-4);
    assert!((number(current, "volume") - 0.6).abs() <= 1e-4);
    assert!(number(current, "sampled_max_von_mises") <= 1e12);
    assert!(number(current, "sample_count") > 0.0);
    let field = std::fs::read(output_dir.join("level-set.csv")).unwrap();
    assert_eq!(field, std::fs::read(output_dir.join("accepted-0001.csv")).unwrap());
    assert_ne!(field, std::fs::read(output_dir.join("baseline-level-set.csv")).unwrap());
    let checks = std::fs::read_to_string(output_dir.join("stress-checks.jsonl")).unwrap();
    assert!(checks.contains("\"refusal\":null"));
    assert_eq!(std::fs::read_to_string(output_dir.join("trajectory.jsonl")).unwrap().lines().count(), 1);
}

#[test]
fn overstressed_baseline_refuses_without_creating_a_study_or_success_result() {
    let scratch = Scratch::new();
    let output_dir = scratch.0.join("refused");
    let output = run(&output_dir, "1e-20", &["3", "1", "0.6"]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("baseline refused"));
    assert!(!output_dir.exists());
}

#[test]
fn malformed_limit_and_budget_refuse_before_side_effects() {
    let scratch = Scratch::new();
    for (i, (limit, extra)) in [
        ("NaN", vec![]),
        ("1e12", vec!["40"]),
        ("1e12", vec!["3", "1", "0.6", "8", "0", "0"]),
        ("1e12", vec!["3", "1", "0.6", "17"]),
    ].into_iter().enumerate() {
        let output_dir = scratch.0.join(format!("refused-{i}"));
        let output = run(&output_dir, limit, &extra);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!output_dir.exists());
    }
}

#[test]
fn existing_output_is_never_overwritten() {
    let scratch = Scratch::new();
    let sentinel = scratch.0.join("summary.json");
    std::fs::write(&sentinel, "owned existing result").unwrap();
    let output = run(&scratch.0, "1e12", &[]);
    assert!(!output.status.success());
    assert_eq!(std::fs::read_to_string(sentinel).unwrap(), "owned existing result");
}
