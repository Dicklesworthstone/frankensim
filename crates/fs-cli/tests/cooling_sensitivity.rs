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
fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let value = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(value.str_field("schema"), Some("frankensim.cooling-network-uq.sensitivity.progress.v1"));
    for field in ["effects", "base_output_std_dev_k", "mean_k", "probability_of_compliance", "confidence_interval"] {
        assert!(value.get(field).is_none(), "published {field} from a partial design");
    }
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
    let result = command(&missing, &missing).args(["--qmc-replicates", "2"])
        .arg("--checkpoint").arg(&output).output().unwrap();
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
    let failed_checkpoint = dir.join("failed.sobol");
    let result = command(&broken, &request).arg("--checkpoint").arg(&failed_checkpoint).output().unwrap();
    assert!(!result.status.success()); assert!(result.stdout.is_empty());
    assert_ne!(result.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(fs::read(&failed_checkpoint).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let refused_destination = dir.join("must-not-resume.sobol");
    let resumed = command(&broken, &request).arg("--resume").arg(&failed_checkpoint)
        .arg("--checkpoint").arg(&refused_destination).output().unwrap();
    assert!(!resumed.status.success());
    assert!(!refused_destination.exists());
    let short = dir.join("timeout.json");
    fs::write(&short, one_input(false).replace("\"wall_seconds\":300", "\"wall_seconds\":0.000000001")).unwrap();
    let result = command(&base, &short).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(result.stdout.is_empty());
    let retained_timeout = dir.join("timeout.sobol");
    let result = command(&base, &short).arg("--checkpoint").arg(&retained_timeout).output().unwrap();
    let report = progress(&result);
    assert_eq!(report.str_field("termination"), Some("wall-time-budget"));
    assert_eq!(number(&report, "samples_evaluated"), 0.0);
    let continued = dir.join("continued.sobol");
    let result = command(&base, &request).arg("--resume").arg(&retained_timeout)
        .arg("--checkpoint").arg(&continued).args(["--max-new-samples", "1"]).output().unwrap();
    assert_eq!(number(&progress(&result), "samples_evaluated"), 1.0);
}

#[test]
fn real_cooling_sobol_retains_partial_hybrid_rows_and_replays_exactly() {
    let dir = scratch("recovery");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("two-inputs.json");
    fs::write(&request, r#"{
        "schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":12,
        "wall_seconds":300,"correlation":{"kind":"independent"},
        "parameters":[
            {"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"uniform","lo":299,"hi":301}},
            {"target":{"kind":"component-power","component":"chip"},"distribution":{"kind":"uniform","lo":0.8,"hi":1.2}}
        ]}"#).unwrap();
    let full_path = dir.join("full.sobol");
    let full = command(&base, &request).arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first_path = dir.join("first.sobol");
    let first = command(&base, &request).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let report = progress(&first);
    assert_eq!(number(&report, "samples_evaluated"), 3.0);
    assert_eq!(number(&report, "completed_rows"), 0.0);
    assert_eq!(number(&report, "next_row_slot"), 3.0);
    assert_eq!(report.str_field("next_evaluation_kind"), Some("hybrid"));
    assert_eq!(report.str_field("next_hybrid_parameter"), Some("component[chip].power_w"));
    let retained = fs::read(&first_path).unwrap();
    assert_eq!(&retained[..8], b"FSSOB001");
    let second_path = dir.join("second.sobol");
    let second = command(&base, &request).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "2"]).output().unwrap();
    let report = progress(&second);
    assert_eq!(number(&report, "samples_evaluated"), 5.0);
    assert_eq!(number(&report, "samples_evaluated_this_run"), 2.0);
    assert_eq!(number(&report, "completed_rows"), 1.0);
    assert_eq!(number(&report, "next_row_ordinal"), 1.0);
    assert_eq!(report.str_field("next_evaluation_kind"), Some("base-b"));
    assert_eq!(fs::read(&first_path).unwrap(), retained);
    let completed_path = dir.join("completed.sobol");
    let completed = command(&base, &request).arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&completed_path).output().unwrap();
    parsed(&completed);
    assert_eq!(completed.stdout, full.stdout);
    assert_eq!(fs::read(&completed_path).unwrap(), fs::read(&full_path).unwrap());
    // A completed checkpoint is terminal even with too little time to solve.
    let short = dir.join("short.json");
    fs::write(&short, fs::read_to_string(&request).unwrap().replace("\"wall_seconds\":300", "\"wall_seconds\":0.000000001")).unwrap();
    let terminal = command(&base, &short).arg("--resume").arg(&completed_path).output().unwrap();
    parsed(&terminal);
    assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn changed_identity_corruption_and_existing_destinations_refuse_before_writing() {
    let dir = scratch("identity");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("original.json");
    fs::write(&request, one_input(false)).unwrap();
    let checkpoint = dir.join("empty.sobol");
    let initial = command(&base, &request).arg("--checkpoint").arg(&checkpoint)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(number(&progress(&initial), "samples_evaluated"), 0.0);
    let saved = fs::read(&checkpoint).unwrap();
    let changed_base = dir.join("changed-base.json");
    fs::write(&changed_base, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    let changed_request = dir.join("changed-seed.json");
    fs::write(&changed_request, one_input(false).replace("\"73\"", "\"74\"")).unwrap();
    let corrupt = dir.join("corrupt.sobol");
    let mut corrupt_bytes = saved.clone(); *corrupt_bytes.last_mut().unwrap() ^= 1;
    fs::write(&corrupt, corrupt_bytes).unwrap();
    for (index, (model, plan, source)) in [
        (&changed_base, &request, &checkpoint),
        (&base, &changed_request, &checkpoint),
        (&base, &request, &corrupt),
    ].into_iter().enumerate() {
        let destination = dir.join(format!("refused-{index}.sobol"));
        let output = command(model, plan).arg("--resume").arg(source)
            .arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success()); assert!(output.stdout.is_empty());
        assert!(!destination.exists());
    }
    let overwrite = command(&base, &request).arg("--resume").arg(&checkpoint)
        .arg("--checkpoint").arg(&checkpoint).output().unwrap();
    assert!(!overwrite.status.success());
    assert_eq!(fs::read(&checkpoint).unwrap(), saved);
}
