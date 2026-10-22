//! Actual cooling child solves, not a substitute thermal callback.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use fs_uq::{CorrelationModel, ParameterUncertainty, PceConfig, PceExecution, PropagationMethod, UqPlan, UqStatus};
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-pce-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}
fn command(base: &Path, request: &Path, layout: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request).args(["--pce", layout]);
    c
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn number(root: &J, name: &str) -> f64 { root.get(name).unwrap().as_f64().unwrap() }
fn inlet_request(samples: usize, wall: f64) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{samples},"wall_seconds":{wall},"temperature_limit_k":303,"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"inlet-temperature","index":0}},"distribution":{{"kind":"uniform","lo":299,"hi":301}}}}]}}"#)
}
fn physical(dir: &Path, original: &str, inlet: f64) -> f64 {
    assert_eq!(original.matches("\"temperature_k\": 300").count(), 1);
    let request = original.replace("\"temperature_k\": 300", &format!("\"temperature_k\": {inlet}"));
    let path = dir.join("physical.json");
    fs::write(&path, request).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(path).output().unwrap();
    number(parsed(&output).get("objective").unwrap(), "value_k")
}

#[test]
fn cooling_pce_matches_independent_library_run_on_the_original_physics() {
    let dir = scratch("physical");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("uq.json");
    fs::write(&request, inlet_request(8, 300.0)).unwrap();
    let output = command(&base, &request, "1,6").output().unwrap();
    let result = parsed(&output);
    assert_eq!(result.str_field("schema"), Some("frankensim.cooling-network-uq.pce.v1"));
    assert_eq!(number(&result, "samples_evaluated"), 8.0);
    assert_eq!(number(&result, "training_samples"), 6.0);
    assert_eq!(number(&result, "validation_samples"), 2.0);
    assert_eq!(result.get("germ_basis").unwrap().as_array().unwrap()[0].as_str(), Some("uniform-01-shifted-legendre"));
    assert!(result.get("mean_k").is_none(), "surrogate mean is not a physical MC mean");
    assert!(result.get("confidence_interval").is_none());
    let mut plan = UqPlan::new("temperature", PropagationMethod::PolynomialChaos, 8)
        .with_parameter(ParameterUncertainty::uniform("inlet", 299.0, 301.0, "K"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(303.0);
    plan.seed = 73;
    let original = fs::read_to_string(&base).unwrap();
    let mut reference = PceExecution::new(&plan, PceConfig::new(1, 6)).unwrap();
    let expected = reference.advance(8, || false, |x| Ok::<_, &str>(physical(&dir, &original, x[0])));
    assert_eq!(expected.status, UqStatus::Complete);
    let summary = expected.summary.unwrap();
    assert!((number(&result, "surrogate_mean_k") - summary.surrogate_mean).abs() < 1e-9);
    assert!((number(&result, "surrogate_variance_k2") - 1.0 / 3.0).abs() < 1e-6);
    assert!((number(&result, "validation_rmse_k") - summary.validation_rmse).abs() < 1e-9);
    assert_eq!(result.get("held_out_probability_of_compliance").unwrap().as_f64(), summary.probability_of_compliance);
    let coefficients = result.get("coefficients_k").unwrap().as_array().unwrap();
    for (actual, expected) in coefficients.iter().zip(reference.model().unwrap().coefficients()) {
        assert!((actual.as_f64().unwrap() - expected).abs() < 1e-9);
    }
    let repeat = command(&base, &request, "1,6").output().unwrap();
    parsed(&repeat);
    assert_eq!(output.stdout, repeat.stdout);
}

#[test]
fn complete_transient_trajectories_remain_the_polynomial_observations() {
    let output = command(&example("nonlinear-contact-pulse.json"),
        &example("uq-nonlinear-contact-pulse.json"), "0,6").output().unwrap();
    let result = parsed(&output);
    assert_eq!(number(&result, "samples_evaluated"), 8.0);
    assert_eq!(number(&result, "basis_terms"), 1.0);
    assert!(String::from_utf8_lossy(&output.stdout).contains("transient-sampled-peak"));
    assert!(number(&result, "surrogate_mean_k").is_finite());
    assert!(number(&result, "validation_rmse_k") >= 0.0);
}

#[test]
fn polynomial_admission_does_not_create_recovery_files_or_ignore_other_methods() {
    let dir = scratch("admission");
    let missing = dir.join("missing.json");
    let saved = dir.join("must-not-exist.bin");
    for extra in [vec!["--qmc-replicates", "2"], vec!["--sensitivity", "sobol"],
        vec!["--mean-control", "adjoint"], vec!["--pce", "0,2"]] {
        let output = command(&missing, &missing, "1,6").args(extra).output().unwrap();
        assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
        assert!(output.stdout.is_empty());
    }
    let output = command(&missing, &missing, "1,6").arg("--checkpoint").arg(&saved).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(!saved.exists());
    let request = dir.join("uq.json");
    fs::write(&request, inlet_request(8, 300.0)).unwrap();
    let base = example("fan-correlated-hotspot.json");
    for layout in ["2,5", "1,7"] {
        let output = command(&base, &request, layout).output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn timeout_and_failed_physics_cannot_publish_a_partial_surrogate() {
    let dir = scratch("failures");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("uq.json");
    fs::write(&request, inlet_request(8, 1e-9)).unwrap();
    let output = command(&base, &request, "1,6").output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
    fs::write(&request, inlet_request(8, 300.0)).unwrap();
    let broken = dir.join("broken.json");
    let original = fs::read_to_string(&base).unwrap();
    assert!(original.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, original.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let output = command(&broken, &request, "1,6").output().unwrap();
    assert!(!output.status.success());
    assert_ne!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
}

fn mean_command(base: &Path, request: &Path, layout: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request).args(["--pce-mean", layout]);
    c
}

#[test]
fn controlled_cooling_mean_matches_independent_physics_with_all_work_counted() {
    use fs_uq::{PceMeanConfig, PceMeanExecution};
    let dir = scratch("controlled");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("uq.json");
    fs::write(&request, inlet_request(12, 300.0)).unwrap();
    let output = mean_command(&base, &request, "1,6,8").output().unwrap();
    let result = parsed(&output);
    assert_eq!(result.str_field("schema"), Some("frankensim.cooling-network-uq.pce-mean.v1"));
    assert_eq!(number(&result, "model_calls_evaluated"), 12.0);
    assert_eq!(number(&result, "fit_evaluations"), 8.0);
    assert_eq!(number(&result, "production_samples_planned"), 4.0);
    assert_eq!(number(&result, "production_samples_evaluated"), 4.0);
    let nested = result.get("fit").unwrap();
    assert_eq!(number(nested, "samples_planned"), 8.0);
    assert_eq!(number(nested, "training_samples"), 6.0);
    assert_eq!(number(nested, "validation_samples"), 2.0);
    assert_ne!(result.str_field("fit_seed"), result.str_field("production_seed"));
    let mut plan = UqPlan::new("temperature", PropagationMethod::MonteCarlo, 12)
        .with_parameter(ParameterUncertainty::uniform("inlet", 299.0, 301.0, "K"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(303.0);
    plan.seed = 73;
    let policy = PceMeanConfig { fit_evaluations: 8, regression: PceConfig::new(1, 6), control_weight: 1.0 };
    let mut reference = PceMeanExecution::new(&plan, policy).unwrap();
    let original = fs::read_to_string(&base).unwrap();
    let mut calls = 0;
    let expected = reference.advance(12, || false, |x| {
        calls += 1;
        Ok::<_, &str>(physical(&dir, &original, x[0]))
    });
    assert_eq!(calls, 12);
    assert_eq!(expected.status, UqStatus::Complete);
    let estimate = expected.estimate.unwrap();
    assert!((number(&result, "mean_k") - estimate.mean).abs() < 1e-9);
    assert!((number(&result, "raw_mean_k") - estimate.raw_mean).abs() < 1e-9);
    assert!(number(&result, "sampling_standard_error_k") < number(&result, "raw_sampling_standard_error_k") * 1e-3);
    assert!((number(&result, "mean_k") - physical(&dir, &original, 300.0)).abs() < 1e-6);
    assert_eq!(result.get("raw_production_probability_of_compliance").unwrap().as_f64(),
        reference.production().report().probability_of_compliance);
    assert_eq!(result.str_field("production_seed").unwrap().parse::<u64>().unwrap(), reference.production().plan().seed);
}

#[test]
fn mean_mode_requires_separate_production_and_never_publishes_timed_out_fits() {
    let dir = scratch("controlled-admission");
    let base = example("fan-correlated-hotspot.json");
    let request = dir.join("uq.json");
    fs::write(&request, inlet_request(8, 300.0)).unwrap();
    let output = mean_command(&base, &request, "1,6,8").output().unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("production"));
    fs::write(&request, inlet_request(12, 1e-9)).unwrap();
    let output = mean_command(&base, &request, "1,6,8").output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
    let missing = dir.join("missing.json");
    let saved = dir.join("must-not-exist.bin");
    let output = mean_command(&missing, &missing, "1,6,8").arg("--checkpoint").arg(&saved).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(!saved.exists());
}

#[test]
fn a_constant_control_preserves_complete_transient_production_statistics() {
    let dir = scratch("controlled-transient");
    let base = example("nonlinear-contact-pulse.json");
    let original = fs::read_to_string(example("uq-nonlinear-contact-pulse.json")).unwrap();
    assert!(original.contains("\"samples\": 8"));
    let request = dir.join("uq.json");
    fs::write(&request, original.replace("\"samples\": 8", "\"samples\": 12")).unwrap();
    let output = mean_command(&base, &request, "0,6,8").output().unwrap();
    let result = parsed(&output);
    assert_eq!(number(&result, "model_calls_evaluated"), 12.0);
    assert_eq!(number(&result, "production_samples_evaluated"), 4.0);
    assert_eq!(number(&result, "mean_k"), number(&result, "raw_mean_k"));
    assert_eq!(number(&result, "sampling_standard_error_k"), number(&result, "raw_sampling_standard_error_k"));
    assert!(String::from_utf8_lossy(&output.stdout).contains("transient-sampled-peak"));
}

#[test]
fn supplied_mixed_law_example_runs_both_original_cooling_workflows() {
    let base = example("fan-correlated-hotspot.json");
    let request = example("uq-polynomial-hotspot.json");
    let surrogate = parsed(&command(&base, &request, "1,20").output().unwrap());
    assert_eq!(number(&surrogate, "samples_evaluated"), 32.0);
    assert_eq!(number(&surrogate, "basis_terms"), 3.0);
    let basis = surrogate.get("germ_basis").unwrap().as_array().unwrap();
    assert_eq!(basis[0].as_str(), Some("uniform-01-shifted-legendre"));
    assert_eq!(basis[1].as_str(), Some("standard-normal-hermite"));
    let controlled = parsed(&mean_command(&base, &request, "1,20,24").output().unwrap());
    assert_eq!(number(&controlled, "model_calls_evaluated"), 32.0);
    assert_eq!(number(&controlled, "fit_evaluations"), 24.0);
    assert_eq!(number(&controlled, "production_samples_evaluated"), 8.0);
    assert_eq!(number(controlled.get("fit").unwrap(), "validation_samples"), 4.0);
    assert!(number(&controlled, "mean_k") > 0.0);
    assert!(controlled.get("confidence_interval").is_none());
}
