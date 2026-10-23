//! Real cooling processes; no substitute model or proposal-frequency decision.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples/cooling-network/fan-correlated-hotspot.json")
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let p = std::env::temp_dir().join(format!("fs-importance-decision-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&p).unwrap(); p
}
fn request(samples: usize, wall: f64, limit: f64) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{samples},"wall_seconds":{wall},"temperature_limit_k":{limit},"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"inlet-temperature","index":0}},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0}}}}]}}"#)
}
fn command(base: &Path, request: &Path, probability: &str, minimum: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request)
        .args(["--importance", "0.2,0", "--compliance-probability", probability,
            "--confidence-alpha", "0.05", "--min-decision-samples", minimum]);
    c
}
fn parsed(output: &Output, code: u8) -> J {
    assert_eq!(output.status.code(), Some(i32::from(code)), "{}", String::from_utf8_lossy(&output.stderr));
    let r = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(r.str_field("schema"), Some("frankensim.cooling-network-uq.importance-compliance.v1"));
    for name in ["mean_k", "percentiles_p05_p50_p95_k", "empirical_probability_of_compliance", "failure_standard_error"] {
        assert!(r.get(name).is_none(), "stopped proposal responses published as {name}");
    }
    r
}
fn number(r: &J, key: &str) -> f64 { r.get(key).unwrap().as_f64().unwrap() }

#[test]
fn both_probability_decisions_stop_actual_cooling_before_the_lifetime_cap() {
    let dir = scratch("decisions");
    let path = dir.join("uq.json");
    for (limit, expected) in [(1000.0, "meets-probability-target"), (1.0, "below-probability-target")] {
        fs::write(&path, request(64, 300.0, limit)).unwrap();
        let output = command(&example(), &path, "0.5", "16").output().unwrap();
        let r = parsed(&output, fs_cli::exit::SUCCESS);
        assert_eq!(r.str_field("status"), Some("decision-reached"));
        assert_eq!(r.str_field("termination"), Some("probability-target"));
        assert_eq!(r.str_field("decision"), Some(expected));
        assert_eq!(number(&r, "samples_evaluated"), 16.0);
        assert_eq!(number(&r, "confidence_samples"), 16.0);
        let ci = r.get("compliance_probability_confidence_sequence").unwrap().as_array().unwrap();
        if limit > 300.0 { assert!(ci[0].as_f64().unwrap() >= 0.5); }
        else { assert!(ci[1].as_f64().unwrap() < 0.5); }
    }
}

#[test]
fn stopping_checkpoint_replays_without_physics_and_cannot_change_its_policy() {
    let dir = scratch("recovery");
    let path = dir.join("uq.json"); fs::write(&path, request(64, 300.0, 1000.0)).unwrap();
    let full_path = dir.join("full.imp");
    let full = command(&example(), &path, "0.5", "16").arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full, fs_cli::exit::SUCCESS);
    let first_path = dir.join("first.imp");
    let first = command(&example(), &path, "0.5", "16").arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "5"]).output().unwrap();
    let progress = parsed(&first, fs_cli::exit::BUDGET);
    assert_eq!(progress.str_field("decision"), Some("indeterminate"));
    assert_eq!(progress.str_field("termination"), Some("sample-chunk"));
    assert_eq!(number(&progress, "samples_evaluated"), 5.0);
    let original = fs::read(&first_path).unwrap();
    let done_path = dir.join("done.imp");
    let done = command(&example(), &path, "0.5", "16").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&done_path).output().unwrap();
    parsed(&done, fs_cli::exit::SUCCESS);
    assert_eq!(full.stdout, done.stdout);
    assert_eq!(fs::read(&full_path).unwrap(), fs::read(&done_path).unwrap());
    assert_eq!(fs::read(&first_path).unwrap(), original);
    let short = dir.join("short.json"); fs::write(&short, request(64, 1e-9, 1000.0)).unwrap();
    let terminal = command(&example(), &short, "0.5", "16").arg("--resume").arg(&done_path).output().unwrap();
    parsed(&terminal, fs_cli::exit::SUCCESS);
    assert_eq!(terminal.stdout, full.stdout);
    for (probability, minimum) in [("0.6", "16"), ("0.5", "17")] {
        let never = dir.join(format!("changed-{probability}-{minimum}.imp"));
        let changed = command(&example(), &path, probability, minimum).arg("--resume").arg(&first_path)
            .arg("--checkpoint").arg(&never).output().unwrap();
        assert_eq!(changed.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
        assert!(changed.stdout.is_empty()); assert!(!never.exists());
        assert!(String::from_utf8_lossy(&changed.stderr).contains("identity differs"));
    }
    // Removing the policy cannot turn a stopped run into a fixed-count run.
    let plain = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network-uq"]).arg(example()).arg(&path)
        .args(["--importance", "0.2,0", "--resume"]).arg(&done_path).output().unwrap();
    assert_eq!(plain.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(plain.stdout.is_empty());
}

#[test]
fn an_exhausted_experiment_with_zero_failures_is_not_an_approval() {
    let dir = scratch("inconclusive");
    let path = dir.join("uq.json"); fs::write(&path, request(8, 300.0, 1000.0)).unwrap();
    let r = parsed(&command(&example(), &path, "0.999", "2").output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("status"), Some("inconclusive"));
    assert_eq!(r.str_field("termination"), Some("sample-budget"));
    assert_eq!(r.str_field("decision"), Some("indeterminate"));
    assert_eq!(number(&r, "samples_evaluated"), 8.0);
    assert_eq!(number(&r, "estimated_failure_probability"), 0.0);
    let ci = r.get("failure_probability_confidence_sequence").unwrap().as_array().unwrap();
    assert!(ci[1].as_f64().unwrap() > 0.001);
    let too_short = command(&example(), &path, "0.5", "9").output().unwrap();
    assert_eq!(too_short.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(too_short.stdout.is_empty());
}

#[test]
fn timeout_and_failed_physics_never_manufacture_a_decision() {
    let dir = scratch("failures");
    let path = dir.join("uq.json"); fs::write(&path, request(64, 1e-9, 1000.0)).unwrap();
    let r = parsed(&command(&example(), &path, "0.5", "16").output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(r.str_field("termination"), Some("wall-time-budget"));
    assert_eq!(number(&r, "confidence_samples"), 0.0);
    assert_eq!(r.get("failure_probability_confidence_sequence"), Some(&J::Null));
    fs::write(&path, request(64, 300.0, 1000.0)).unwrap();
    let broken = dir.join("broken.json");
    let base = fs::read_to_string(example()).unwrap();
    assert!(base.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, base.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let checkpoint = dir.join("failed.imp");
    let failed = command(&broken, &path, "0.5", "16").arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert_eq!(failed.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(failed.stdout.is_empty());
    assert!(fs::read(&checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
}
