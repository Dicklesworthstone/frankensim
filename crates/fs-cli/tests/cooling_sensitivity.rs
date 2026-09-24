//! Real binary and physical cooling producer; no surrogate or replacement solver.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-sobol-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", "cooling-network-uq"]).arg(base).arg(request)
        .args(["--sensitivity", "sobol"]); command
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let value = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(value.str_field("schema"), Some("frankensim.cooling-network-uq.sensitivity.v1"));
    value
}
fn number(root: &J, name: &str) -> f64 { root.get(name).unwrap().as_f64().unwrap() }
fn one_input(transient: bool) -> String {
    let target = if transient { r#"{"kind":"initial-temperature"}"# }
        else { r#"{"kind":"inlet-temperature","index":0}"# };
    let qoi = if transient { r#", "qoi":{"kind":"transient-sampled-peak"}"# } else { "" };
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":12,
        "wall_seconds":300,"correlation":{{"kind":"independent"}},
        "parameters":[{{"target":{target},"distribution":{{"kind":"uniform","lo":299,"hi":301}}}}]{qoi}}}"#)
}

#[test]
fn actual_cooling_single_input_attribution_replays_and_has_no_probability_claim() {
    let dir = scratch("steady");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("inlet.json"); fs::write(&request, one_input(false)).unwrap();
    let first = command(&base, &request).output().unwrap();
    let report = parsed(&first);
    assert_eq!(number(&report, "samples_evaluated"), 12.0);
    assert_eq!(number(&report, "base_samples"), 4.0);
    assert_eq!(number(&report, "evaluations_per_row"), 3.0);
    assert!(number(&report, "base_output_std_dev_k") > 0.0);
    let effects = report.get("effects").unwrap().as_array().unwrap();
    assert_eq!(effects.len(), 1);
    assert_eq!(effects[0].str_field("parameter"), Some("inlet[0].temperature_k"));
    // With a single uncertain input every hybrid IS its B base. Agreement
    // requires deterministic complete physical solves, not a surrogate fit.
    assert_eq!(number(&effects[0], "first_order"), 1.0);
    assert!(number(&effects[0], "total_order") > 0.0);
    for field in ["confidence_interval", "probability_of_compliance", "mean_k", "sampling_standard_error_k"] {
        assert!(report.get(field).is_none());
    }
    let again = command(&base, &request).output().unwrap(); parsed(&again);
    assert_eq!(again.stdout, first.stdout);
}

#[test]
fn physical_transient_sensitivity_uses_complete_trajectory_peaks() {
    let dir = scratch("transient");
    let request = dir.join("initial.json"); fs::write(&request, one_input(true)).unwrap();
    let output = command(&example("nonlinear-contact-pulse.json"), &request).output().unwrap();
    let report = parsed(&output);
    let qoi = report.get("qoi").unwrap();
    assert_eq!(qoi.str_field("kind"), Some("transient-sampled-peak"));
    assert_eq!(qoi.str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(qoi.str_field("temporal_scope"), Some("initial-state-and-accepted-endpoints"));
    assert_eq!(number(&report, "samples_evaluated"), 12.0);
    let effects = report.get("effects").unwrap().as_array().unwrap();
    assert_eq!(number(&effects[0], "first_order"), 1.0);
    assert!(number(&report, "base_output_std_dev_k") > 0.0);
}

#[test]
fn multiple_physical_inputs_are_named_and_all_hybrid_work_is_counted() {
    let dir = scratch("multi");
    let request = dir.join("four-inputs.json");
    let text = fs::read_to_string(example("uq-fan-hotspot.json")).unwrap();
    assert!(text.contains("\"samples\": 8"));
    fs::write(&request, text.replace("\"samples\": 8", "\"samples\": 12")).unwrap();
    let report = parsed(&command(&example("fan-correlated-hotspot.json"), &request).output().unwrap());
    assert_eq!(number(&report, "base_samples"), 2.0);
    assert_eq!(number(&report, "samples_evaluated"), 12.0);
    assert_eq!(number(&report, "evaluations_per_row"), 6.0);
    let effects = report.get("effects").unwrap().as_array().unwrap();
    let expected = ["air.density_kg_m3", "inlet[0].temperature_k", "hydraulics.fan.speed_ratio", "component[chip].power_w"];
    assert_eq!(effects.len(), expected.len());
    for (effect, name) in effects.iter().zip(expected) {
        assert_eq!(effect.str_field("parameter"), Some(name));
        assert!(number(effect, "first_order").is_finite());
        assert!(number(effect, "total_order").is_finite());
    }
}

#[test]
fn incompatible_modes_fail_before_reading_inputs_or_creating_outputs() {
    let dir = scratch("modes");
    let missing = dir.join("missing.json");
    let output = dir.join("must-not-exist.bin");
    for tail in [vec!["--qmc-replicates", "2"], vec!["--sensitivity", "sobol"],
        vec!["--compliance-probability", "0.5", "--confidence-alpha", "0.05", "--min-decision-samples", "2"]] {
        let result = command(&missing, &missing).args(tail).output().unwrap();
        assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
        assert!(result.stdout.is_empty());
    }
    let result = command(&missing, &missing).arg("--checkpoint").arg(&output).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(!output.exists());
}

#[test]
fn invalid_layout_dependence_and_physical_failure_never_publish_indices() {
    let dir = scratch("refusals");
    let base = example("fan-correlated-hotspot.json");
    for (name, text) in [
        ("layout", one_input(false).replace("\"samples\":12", "\"samples\":11")),
        ("dependence", one_input(false).replace("\"independent\"", "\"unknown\"")),
    ] {
        let request = dir.join(format!("{name}.json")); fs::write(&request, text).unwrap();
        let result = command(&base, &request).output().unwrap();
        assert!(!result.status.success()); assert!(result.stdout.is_empty());
    }
    let request = dir.join("valid.json"); fs::write(&request, one_input(false)).unwrap();
    let broken = dir.join("broken.json");
    let base_text = fs::read_to_string(&base).unwrap();
    assert!(base_text.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, base_text.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let result = command(&broken, &request).output().unwrap();
    assert!(!result.status.success()); assert!(result.stdout.is_empty());
    assert_ne!(result.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    let short = dir.join("timeout.json");
    fs::write(&short, one_input(false).replace("\"wall_seconds\":300", "\"wall_seconds\":0.000000001")).unwrap();
    let result = command(&base, &short).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(result.stdout.is_empty());
}
