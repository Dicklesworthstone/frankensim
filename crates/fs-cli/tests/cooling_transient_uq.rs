//! Product tests: every observation comes from the actual cooling binary.
//! Direct single runs establish the observable; no alternate heat solver is
//! substituted for the transient, contact or nonlinear material implementation.
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonlinear-contact-pulse.json"));

fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-transient-uq-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}
fn write(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    fs::write(&path, text).unwrap();
    path
}
fn uq(base: &Path, plan: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", "cooling-network-uq"]).arg(base).arg(plan);
    command
}
fn direct(base: &Path) -> J {
    document(&Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(base).output().unwrap())
}
fn document(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn value(document: &J, path: &[&str]) -> f64 {
    document.path(path).and_then(J::as_f64).unwrap()
}
fn close(a: f64, b: f64) { assert!((a - b).abs() < 1.0e-7, "{a} != {b}"); }
fn short_base() -> String {
    assert!(BASE.contains("\"duration_s\":30") && BASE.contains("\"duration_s\":120"));
    BASE.replace("\"duration_s\":30", "\"duration_s\":4")
        .replace("\"duration_s\":120", "\"duration_s\":4")
}
fn plan(target: &str, lo: f64, hi: f64, samples: usize, ceiling: f64) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"79","samples":{samples},"wall_seconds":300,"qoi":{{"kind":"transient-sampled-peak"}},"temperature_limit_k":{ceiling},"correlation":{{"kind":"independent"}},"parameters":[{{"target":{target},"distribution":{{"kind":"uniform","lo":{lo},"hi":{hi}}}}]}}"#)
}
fn progress(output: &Output, count: usize) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)),
        "{}", String::from_utf8_lossy(&output.stderr));
    let document = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(value(&document, &["samples_evaluated"]), count as f64);
    assert!(document.get("mean_k").is_none());
    document
}

#[test]
fn zero_uncertainty_measures_the_peak_not_the_final_temperature() {
    let dir = scratch("peak");
    let base = write(&dir, "base.json", BASE);
    let nominal = direct(&base);
    let peak = value(&nominal, &["transient", "sampled_peak_objective_k"]);
    let final_value = value(&nominal, &["objective", "value_k"]);
    assert!(peak > final_value + 1.0, "fixture must distinguish peak from final temperature");
    let ceiling = 0.5 * (peak + final_value);
    let target = r#"{"kind":"interval-power-scale","interval":0}"#;
    let request = write(&dir, "fixed.json", &plan(target, 1.0, 1.0, 2, ceiling));
    let result = document(&uq(&base, &request).output().unwrap());
    close(value(&result, &["mean_k"]), peak);
    close(value(&result, &["std_dev_k"]), 0.0);
    assert_eq!(value(&result, &["empirical_probability_of_compliance"]), 0.0);
    assert_eq!(value(&result, &["samples_evaluated"]), 2.0);
    assert_eq!(result.path(&["qoi", "kind"]).and_then(J::as_str), Some("transient-sampled-peak"));
    assert_eq!(result.path(&["qoi", "continuous_time_bound"]), Some(&J::Bool(false)));
    // Zero is a real off interval, not an invalid sample or a request to redraw.
    let off = write(&dir, "off.json", &plan(target, 0.0, 0.0, 2, ceiling));
    let result = document(&uq(&base, &off).output().unwrap());
    close(value(&result, &["mean_k"]), 300.0);
    assert_eq!(value(&result, &["empirical_probability_of_compliance"]), 1.0);
}

#[test]
fn fixed_repetition_uses_the_all_cycle_peak_including_initial_history() {
    let dir = scratch("cycles");
    let text = short_base().replace("\"initial_temperature_k\":300", "\"initial_temperature_k\":350")
        .replace("\"power_scale\":1", "\"power_scale\":0")
        .replace("\"max_step_s\":2", "\"repeat\":{\"cycles\":3,\"max_total_steps\":100},\"max_step_s\":2");
    let base = write(&dir, "base.json", &text);
    let nominal = direct(&base);
    let peak = value(&nominal, &["repeated_cycles", "sampled_peak_objective_k"]);
    let last_cycle_peak = value(&nominal, &["transient", "sampled_peak_objective_k"]);
    assert!(peak > last_cycle_peak + 1.0e-5, "fixture must distinguish earlier-cycle heat");
    let request = write(&dir, "fixed.json", &plan(r#"{"kind":"initial-temperature"}"#, 350.0, 350.0, 2, 340.0));
    let result = document(&uq(&base, &request).output().unwrap());
    close(value(&result, &["mean_k"]), peak);
    assert_eq!(value(&result, &["qoi", "cycles"]), 3.0);
}

#[test]
fn uncertain_trajectories_resume_after_chunk_and_timeout_without_rerunning_the_prefix() {
    let dir = scratch("replay");
    let base = write(&dir, "base.json", &short_base());
    let text = plan(r#"{"kind":"interval-power-scale","interval":0}"#, 0.8, 1.2, 4, 305.0);
    let request = write(&dir, "uq.json", &text);
    let full_path = dir.join("full.bin");
    let full = uq(&base, &request).arg("--checkpoint").arg(&full_path).output().unwrap();
    assert!(value(&document(&full), &["std_dev_k"]) > 0.0);
    let prefix = dir.join("prefix.bin");
    let first = uq(&base, &request).arg("--checkpoint").arg(&prefix)
        .args(["--max-new-samples", "2"]).output().unwrap();
    progress(&first, 2);
    let retained = fs::read(&prefix).unwrap();
    let short = write(&dir, "short.json", &text.replace("\"wall_seconds\":300", "\"wall_seconds\":1e-9"));
    let paused_path = dir.join("paused.bin");
    let paused = uq(&base, &short).arg("--resume").arg(&prefix)
        .arg("--checkpoint").arg(&paused_path).output().unwrap();
    let paused_report = progress(&paused, 2);
    assert_eq!(paused_report.str_field("termination"), Some("wall-time-budget"));
    let final_path = dir.join("final.bin");
    let resumed = uq(&base, &request).arg("--resume").arg(&paused_path)
        .arg("--checkpoint").arg(&final_path).output().unwrap();
    document(&resumed);
    assert_eq!(full.stdout, resumed.stdout);
    assert_eq!(fs::read(full_path).unwrap(), fs::read(final_path).unwrap());
    assert_eq!(fs::read(prefix).unwrap(), retained);
}

#[test]
fn sequential_confidence_counts_trajectories_and_stops_before_the_cap() {
    let dir = scratch("sequential");
    let base = write(&dir, "base.json", &short_base());
    let request = write(&dir, "uq.json", &plan(r#"{"kind":"interval-power-scale","interval":0}"#,
        0.0, 0.0, 64, 1000.0));
    let checkpoint = dir.join("decided.bin");
    let options = ["--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "16"];
    let output = uq(&base, &request).args(options).arg("--checkpoint").arg(&checkpoint).output().unwrap();
    let result = document(&output);
    assert_eq!(result.str_field("decision"), Some("meets-probability-target"));
    assert_eq!(value(&result, &["samples_evaluated"]), 16.0);
    assert_eq!(value(&result, &["samples_planned"]), 64.0);
    assert_eq!(result.path(&["qoi", "observation"]).and_then(J::as_str), Some("one-completed-trajectory"));
    let resumed = uq(&base, &request).args(options).arg("--resume").arg(checkpoint).output().unwrap();
    document(&resumed);
    assert_eq!(output.stdout, resumed.stdout);
}

#[test]
fn ambiguous_observables_inactive_controls_and_variable_horizons_refuse_before_output() {
    let dir = scratch("admission");
    let base_text = short_base();
    let target = r#"{"kind":"interval-power-scale","interval":0}"#;
    let valid = plan(target, 1.0, 1.0, 2, 305.0);
    let cases = [
        (base_text.clone(), valid.replace("\"qoi\":{\"kind\":\"transient-sampled-peak\"},", "")),
        (base_text.clone(), plan(r#"{"kind":"fan-speed-ratio"}"#, 1.0, 1.0, 2, 305.0)),
        (base_text.clone(), plan(r#"{"kind":"interval-power-scale","interval":9}"#, 1.0, 1.0, 2, 305.0)),
        (base_text.clone(), plan(r#"{"kind":"volumetric-heat-capacity"}"#, 1e6, 1e6, 2, 305.0)),
        (base_text.replace("\"max_step_s\":2", "\"repeat\":{\"until_periodic\":{}},\"max_step_s\":2"), valid.clone()),
        (base_text.replace("\"max_step_s\":2", "\"power_design\":{},\"max_step_s\":2"), valid),
    ];
    for (i, (base_text, plan_text)) in cases.iter().enumerate() {
        let base = write(&dir, &format!("base-{i}.json"), base_text);
        let request = write(&dir, &format!("uq-{i}.json"), plan_text);
        let destination = dir.join(format!("must-not-exist-{i}.bin"));
        let output = uq(&base, &request).arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(!destination.exists(), "bad admission reserved an output");
    }
}

#[test]
fn child_failure_does_not_turn_a_partial_trajectory_into_a_sample() {
    let dir = scratch("model-failure");
    let base = write(&dir, "base.json", &short_base());
    // The first interval runs normally, but the next interval's enormous load
    // cannot be completed inside the material's bounded temperature support.
    let request = write(&dir, "uq.json", &plan(r#"{"kind":"interval-power-scale","interval":1}"#,
        1e6, 1e6, 2, 305.0));
    let checkpoint = dir.join("failed.bin");
    let output = uq(&base, &request).arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty(), "no peak may be taken from an unfinished trajectory");
    assert!(fs::read(&checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
}
