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

fn resume(checkpoint: &Path, output: &Path, options: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-stress"))
        .args(["--projected", "--resume"]).arg(checkpoint).arg(output)
        .args(options).output().unwrap()
}

fn assert_refused(output: Output, path: &Path, reason: &str) {
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains(reason), "{}",
        String::from_utf8_lossy(&output.stderr));
    assert!(!path.exists(), "refused resume must not create a new study");
}

#[test]
fn disk_resume_reproduces_uninterrupted_geometry_stress_and_candidate_history() {
    let scratch = Scratch::new();
    let full = scratch.0.join("full");
    let first = scratch.0.join("first");
    let second = scratch.0.join("second");
    let args = ["3", "2", "0.6", "8", "0", "300", "--checkpoint"];
    let uninterrupted = run(&full, "1e12", &args);
    assert!(matches!(uninterrupted.status.code(), Some(0 | 11)), "{}",
        String::from_utf8_lossy(&uninterrupted.stderr));
    let mut paused_args = args.to_vec();
    paused_args.extend(["--pause-after", "1"]);
    let paused = run(&first, "1e12", &paused_args);
    assert_eq!(paused.status.code(), Some(6), "{}", String::from_utf8_lossy(&paused.stderr));
    let input = first.join("checkpoint-0001.fscp");
    let original_bytes = std::fs::read(&input).unwrap();
    let resumed = resume(&input, &second, &["--wall-seconds", "300"]);
    assert_eq!(resumed.status.code(), uninterrupted.status.code(), "{}",
        String::from_utf8_lossy(&resumed.stderr));
    assert_eq!(std::fs::read(&input).unwrap(), original_bytes);
    let summary = std::fs::read_to_string(second.join("summary.json")).unwrap();
    assert!(summary.contains("\"resumed\":true"));
    assert!(summary.contains("\"baseline_scope\":\"current_segment\""));
    assert_eq!(number(&summary, "start_iteration"), 1.0);
    let ordinal = number(&summary, "accepted_updates") as usize;
    assert_eq!(std::fs::read(full.join("level-set.csv")).unwrap(),
        std::fs::read(second.join("level-set.csv")).unwrap());
    let checkpoint = format!("checkpoint-{ordinal:04}.fscp");
    // Includes exact node bits, all policies, global ordinal, AL multiplier,
    // independent mechanics and stress -- not just a rounded output comparison.
    assert_eq!(std::fs::read(full.join(&checkpoint)).unwrap(),
        std::fs::read(second.join(&checkpoint)).unwrap());
    for name in ["trajectory.jsonl", "attempts.jsonl", "stress-checks.jsonl"] {
        let mut joined = std::fs::read(first.join(name)).unwrap();
        joined.extend(std::fs::read(second.join(name)).unwrap());
        assert_eq!(std::fs::read(full.join(name)).unwrap(), joined, "{name}");
    }
    // A fully consumed checkpoint may be inspected/resumed without fabricating
    // another accepted iteration or silently extending the original total.
    if ordinal == 2 {
        let completed = scratch.0.join("completed");
        let result = resume(&second.join(checkpoint), &completed, &["--wall-seconds", "300"]);
        assert!(result.status.success());
        assert!(std::fs::read(completed.join("trajectory.jsonl")).unwrap().is_empty());
    }
}

#[test]
fn checkpoint_corruption_counts_schema_build_and_metrics_fail_closed() {
    let scratch = Scratch::new();
    let source = scratch.0.join("source");
    let paused = run(&source, "1e12", &["3", "2", "0.6", "8", "0", "300", "--pause-after", "0"]);
    assert_eq!(paused.status.code(), Some(6), "{}", String::from_utf8_lossy(&paused.stderr));
    let bytes = std::fs::read(source.join("checkpoint-0000.fscp")).unwrap();
    let header = 65 + b"fs-marquee-projected-stress-checkpoint-v1\n".len();
    let words = header + 64;
    let mut cases = Vec::new();
    cases.push((bytes[..64].to_vec(), "envelope"));
    cases.push((bytes[..bytes.len() - 1].to_vec(), "hash mismatch"));
    let mut corrupted = bytes.clone();
    *corrupted.last_mut().unwrap() ^= 1;
    cases.push((corrupted, "hash mismatch"));
    // Recompute the transport digest: these must fail semantic admission, not
    // merely the corruption check. A digest is intentionally not authenticity.
    for (offset, replacement, reason) in [
        (header - 2, vec![b'2'], "schema"),
        (header, vec![b'x'], "original executable"),
        (words + 8 * 8, u64::MAX.to_le_bytes().to_vec(), "count exceeds"),
        (words + 9 * 8, f64::NAN.to_bits().to_le_bytes().to_vec(), "non-finite"),
        (words + 35 * 8, 1u64.to_le_bytes().to_vec(), "fixed traces"),
        (words + 28 * 8, 0.0f64.to_bits().to_le_bytes().to_vec(), "restored mechanics/stress differ"),
    ] {
        let mut changed = bytes.clone();
        changed[offset..offset + replacement.len()].copy_from_slice(&replacement);
        let hash = fs_ledger::hash_bytes(&changed[65..]).to_string();
        changed[..64].copy_from_slice(hash.as_bytes());
        cases.push((changed, reason));
    }
    for (index, (changed, reason)) in cases.into_iter().enumerate() {
        let input = scratch.0.join(format!("bad-{index}.fscp"));
        let output = scratch.0.join(format!("refused-{index}"));
        std::fs::write(&input, changed).unwrap();
        assert_refused(resume(&input, &output, &["--wall-seconds", "300"]), &output, reason);
    }
}

#[test]
fn resume_requires_budget_and_refuses_policy_overrides_and_existing_output() {
    let scratch = Scratch::new();
    let missing = scratch.0.join("not-read.fscp");
    for (i, options) in [
        vec![], vec!["--wall-seconds", "0"],
        vec!["--wall-seconds", "300", "--stress-limit", "1e30"],
        vec!["--wall-seconds", "300", "--pause-after", "201"],
        vec!["--wall-seconds", "300", "--wall-seconds", "400"],
    ].into_iter().enumerate() {
        let output = scratch.0.join(format!("refused-{i}"));
        let result = resume(&missing, &output, &options);
        assert!(!result.status.success());
        assert!(result.stdout.is_empty());
        assert!(!output.exists());
    }
    let sentinel = scratch.0.join("owned.txt");
    std::fs::write(&sentinel, "keep").unwrap();
    let result = resume(&missing, &scratch.0, &["--wall-seconds", "300"]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("already exists"));
    assert_eq!(std::fs::read_to_string(sentinel).unwrap(), "keep");
}
