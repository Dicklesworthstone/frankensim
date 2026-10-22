//! Exercise the actual cooling binary, its child solver and pilot/production split.
#[path = "../src/json_read.rs"]
mod json;
#[path = "cooling_importance_pilot/recovery.rs"]
mod recovery;
use json::JsonValue as J;
use fs_uq::{
    CorrelationModel, ImportanceConfig, ImportanceFitConfig, ParameterUncertainty,
    PilotImportanceConfig, PilotImportanceExecution, PropagationMethod, UqPlan, UqStatus,
};
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

const LAYOUT: &str = "2,3,0.25,0.5,1,6";
fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-importance-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}
fn command(base: &Path, request: &Path, layout: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", "cooling-network-uq"]).arg(base).arg(request)
        .args(["--importance-pilot", layout]);
    command
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn number(value: &J, field: &str) -> f64 { value.get(field).unwrap().as_f64().unwrap() }
fn request(samples: usize, wall: f64, transient: bool) -> String {
    let target = if transient { r#"{"kind":"initial-temperature"}"# }
        else { r#"{"kind":"inlet-temperature","index":0}"# };
    let qoi = if transient { r#", "qoi":{"kind":"transient-sampled-peak"}"# } else { "" };
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{samples},"wall_seconds":{wall},"temperature_limit_k":303,"correlation":{{"kind":"independent"}},"parameters":[{{"target":{target},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0.2}}}}]{qoi}}}"#)
}
fn physical(dir: &Path, original: &str, inlet: f64) -> f64 {
    assert_eq!(original.matches("\"temperature_k\": 300").count(), 1);
    let text = original.replace("\"temperature_k\": 300", &format!("\"temperature_k\": {inlet}"));
    let path = dir.join("physical.json");
    fs::write(&path, text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(path).output().unwrap();
    number(parsed(&output).get("objective").unwrap(), "value_k")
}
fn optional_matches(report: &J, name: &str, expected: Option<f64>) {
    let actual = report.get(name).unwrap().as_f64();
    match (actual, expected) {
        (None, None) => {}
        (Some(a), Some(b)) => assert!((a - b).abs() <= 1e-11 * b.abs().max(1.0), "{name}: {a} vs {b}"),
        other => panic!("{name}: availability mismatch {other:?}"),
    }
}

#[test]
fn physical_pilots_and_production_match_an_independent_library_reconstruction() {
    let dir = scratch("reference");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json");
    fs::write(&input, request(12, 300.0, false)).unwrap();
    let output = command(&base, &input, LAYOUT).output().unwrap();
    let actual = parsed(&output);
    assert_eq!(actual.str_field("schema"), Some("frankensim.cooling-network-uq.importance-pilot.v1"));
    assert_eq!(number(&actual, "samples_evaluated"), 12.0);
    assert_eq!(number(&actual, "pilot_samples_evaluated"), 6.0);
    assert_eq!(number(&actual, "production_samples"), 6.0);
    let mut plan = UqPlan::new("temperature", PropagationMethod::MonteCarlo, 12)
        .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 0.2, "K"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(303.0);
    plan.seed = 73;
    let settings = PilotImportanceConfig {
        pilot_rounds: 2, samples_per_pilot: 3,
        initial_proposal: ImportanceConfig { germ_shift: vec![0.0], defensive_fraction: 0.25 },
        fit: ImportanceFitConfig { elite_fraction: 0.5, relaxation: 1.0, maximum_shift_norm: 6.0 },
    };
    let original = fs::read_to_string(&base).unwrap();
    let mut expected_run = PilotImportanceExecution::new(&plan, settings).unwrap();
    let expected = expected_run.advance(12, || false, |x| Ok::<_, &str>(physical(&dir, &original, x[0])));
    assert_eq!(expected.status, UqStatus::Complete);
    let estimate = expected.estimate.unwrap();
    optional_matches(&actual, "estimated_failure_probability", Some(estimate.failure_probability));
    optional_matches(&actual, "failure_probability_standard_error", estimate.standard_error);
    optional_matches(&actual, "failure_effective_samples", estimate.failure_effective_samples);
    optional_matches(&actual, "mean_likelihood_weight", Some(estimate.mean_likelihood_weight));
    let seed = expected.current_seed.to_string();
    assert_eq!(actual.str_field("production_seed"), Some(seed.as_str()));
    let pilots = actual.get("pilots").unwrap().as_array().unwrap();
    assert_eq!(pilots.len(), 2);
    for (row, reference) in pilots.iter().zip(&expected.completed_pilots) {
        assert_eq!(number(row, "samples"), 3.0);
        assert_eq!(number(row, "selected"), reference.fit.selected as f64);
        assert_eq!(number(row, "observed_failures"), reference.fit.observed_failures as f64);
        optional_matches(row, "cutoff_k", Some(reference.fit.cutoff));
    }
    // No target quantiles or temperature distribution may be inferred from
    // tilted proposal samples or from adaptation observations.
    for field in ["mean_k", "percentiles_p05_p50_p95_k", "empirical_probability_of_compliance", "confidence_interval"] {
        assert!(actual.get(field).is_none(), "unexpected target-distribution claim {field}");
    }
    let repeat = command(&base, &input, LAYOUT).output().unwrap();
    parsed(&repeat);
    assert_eq!(output.stdout, repeat.stdout);
}

#[test]
fn singular_joint_gaussian_inputs_keep_the_original_physical_law() {
    let dir = scratch("joint");
    let input = dir.join("uq.json");
    fs::write(&input, r#"{
        "schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":12,
        "wall_seconds":300,"temperature_limit_k":303,
        "correlation":{"kind":"joint-gaussian","matrix":[[1,1],[1,1]]},
        "parameters":[
            {"target":{"kind":"inlet-temperature","index":0},"distribution":{"kind":"gaussian","mean":300,"std_dev":0.2}},
            {"target":{"kind":"component-power","component":"chip"},"distribution":{"kind":"gaussian","mean":1,"std_dev":0.01}}
        ]}"#).unwrap();
    let actual = parsed(&command(&example("fan-correlated-hotspot.json"), &input, LAYOUT).output().unwrap());
    assert_eq!(actual.str_field("correlation"), Some("joint-gaussian"));
    assert_eq!(number(&actual, "samples_evaluated"), 12.0);
    let proposal = actual.get("production_proposal").unwrap();
    assert_eq!(proposal.get("independent_normal_germ_shift").unwrap().as_array().unwrap().len(), 2);
    let parameters = actual.get("parameters").unwrap().as_array().unwrap();
    assert_eq!(parameters[0].str_field("name"), Some("inlet[0].temperature_k"));
    assert_eq!(parameters[1].str_field("name"), Some("component[chip].power_w"));
}

#[test]
fn complete_transient_peaks_remain_one_observation_per_physical_call() {
    let dir = scratch("transient");
    let input = dir.join("uq.json");
    fs::write(&input, request(8, 300.0, true)).unwrap();
    let output = command(&example("nonlinear-contact-pulse.json"), &input, "1,2,0.25,0.5,1,6").output().unwrap();
    let report = parsed(&output);
    assert_eq!(number(&report, "samples_evaluated"), 8.0);
    assert_eq!(number(&report, "pilot_samples_evaluated"), 2.0);
    assert_eq!(number(&report, "production_samples"), 6.0);
    let qoi = report.get("qoi").unwrap();
    assert_eq!(qoi.str_field("kind"), Some("transient-sampled-peak"));
    assert_eq!(qoi.str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(qoi.str_field("temporal_scope"), Some("initial-state-and-accepted-endpoints"));
}

#[test]
fn incompatible_flags_refuse_before_input_reads_or_checkpoint_reservation() {
    let dir = scratch("modes");
    let missing = dir.join("missing.json");
    for extra in [
        vec!["--qmc-replicates", "2"], vec!["--sensitivity", "sobol"],
        vec!["--sensitivity", "shapley"], vec!["--pce", "0,2"],
        vec!["--pce-mean", "0,2,4"], vec!["--mean-control", "adjoint"],
        vec!["--importance-pilot", LAYOUT], vec!["--importance", "0.2,0"],
        vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2"],
    ] {
        let result = command(&missing, &missing, LAYOUT).args(extra).output().unwrap();
        assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
        assert!(result.stdout.is_empty());
    }
    let destination = dir.join("must-not-exist.bin");
    let result = command(&missing, &missing, LAYOUT).args(["--qmc-replicates", "2"]).arg("--checkpoint").arg(&destination).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(!destination.exists());
    let help = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network-uq", "--help"]).output().unwrap();
    assert!(parsed(&help).str_field("help").unwrap().contains("--importance-pilot"));
}

#[test]
fn missing_event_unsupported_law_and_missing_production_budget_never_publish() {
    let dir = scratch("plan");
    let base = example("fan-correlated-hotspot.json");
    let original = request(12, 300.0, false);
    for text in [
        original.replace(",\"temperature_limit_k\":303", ""),
        request(7, 300.0, false),
        original.replace("\"kind\":\"independent\"", "\"kind\":\"unknown\""),
        original.replace("\"kind\":\"gaussian\",\"mean\":300,\"std_dev\":0.2",
            "\"kind\":\"uniform\",\"lo\":299,\"hi\":301"),
    ] {
        let input = dir.join("uq.json"); fs::write(&input, text).unwrap();
        let result = command(&base, &input, LAYOUT).output().unwrap();
        assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
        assert!(result.stdout.is_empty());
    }
}

#[test]
fn timeout_and_failed_cooling_do_not_become_zero_risk_or_partial_statistics() {
    let dir = scratch("failure");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json");
    fs::write(&input, request(12, 1e-9, false)).unwrap();
    let result = command(&base, &input, LAYOUT).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(result.stdout.is_empty());
    fs::write(&input, request(12, 300.0, false)).unwrap();
    let text = fs::read_to_string(&base).unwrap();
    assert!(text.contains("\"linear_iterations\": 20000"));
    let broken = dir.join("broken.json");
    fs::write(&broken, text.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let result = command(&broken, &input, LAYOUT).output().unwrap();
    assert_eq!(result.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(result.stdout.is_empty());
}
