//! Real-binary checks for the hard-volume elasticity mode, not a shell mock.
#![cfg(feature = "marquee")]

use std::process::Command;

fn output_path(tag: &str) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!(
        "frankensim-projected-{tag}-{}-{}", std::process::id(),
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos(),
    ));
    assert!(!path.exists());
    path
}

#[test]
fn g0_projected_cli_runs_real_same_volume_elasticity_and_exports_accepted_fields() {
    let path = output_path("solve");
    let output = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity"))
        .arg("--projected").arg(&path).args(["3", "1", "0.6", "8"])
        .output().expect("execute the real elasticity binary");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let summary = std::fs::read_to_string(path.join("summary.json")).unwrap();
    assert!(summary.contains("\"accepted_updates\":1"));
    assert!(summary.contains("\"status\":\"iteration_limit\""));
    assert!(summary.contains("\"converged\":false"));
    assert!(summary.contains("\"authority\":\"estimated\""));
    assert_eq!(String::from_utf8(output.stdout).unwrap().trim(), summary.trim());
    let baseline = std::fs::read_to_string(path.join("baseline-level-set.csv")).unwrap();
    let accepted = std::fs::read_to_string(path.join("level-set.csv")).unwrap();
    assert_ne!(baseline, accepted, "the product must export actual changed geometry");
    assert_eq!(baseline.lines().count(), 82);
    assert_eq!(accepted.lines().count(), 82);
    assert_eq!(std::fs::read_to_string(path.join("trajectory.jsonl")).unwrap().lines().count(), 1);
    assert!(!std::fs::read_to_string(path.join("attempts.jsonl")).unwrap().is_empty());
    // Do not remove result artifacts: retained fields support failure diagnosis.
}

#[test]
fn g4_projected_cli_refuses_invalid_budget_without_output_side_effects() {
    let path = output_path("invalid");
    let output = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity"))
        .arg("--projected").arg(&path).args(["3", "0", "0.6", "8"])
        .output().unwrap();
    assert!(!output.status.success());
    assert!(!path.exists());
    assert!(output.stdout.is_empty());
}
