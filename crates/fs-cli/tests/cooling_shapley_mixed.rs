//! Mixed-law physical attribution through the real cooling binary.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, ShapleyExecution, UqPlan, UqStatus};
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-shapley-mixed-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path, workers: Option<&str>) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", "cooling-network-uq"]).arg(base).arg(request)
        .args(["--sensitivity", "shapley"]);
    if let Some(workers) = workers { command.args(["--workers", workers]); }
    command
}
fn document(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn number(doc: &J, key: &str) -> f64 { doc.get(key).unwrap().as_f64().unwrap() }
fn request(wall: f64) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":12,"wall_seconds":{wall},
        "correlation":{{"kind":"independent"}},"parameters":[
        {{"target":{{"kind":"inlet-temperature","index":0}},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0.2}}}},
        {{"target":{{"kind":"component-power","component":"chip"}},"distribution":{{"kind":"uniform","lo":0.9,"hi":1.1}}}}]}}"#)
}
fn physical(dir: &Path, original: &str, x: &[f64]) -> f64 {
    let mut text = original.to_string();
    for (needle, value) in [("\"temperature_k\": 300", x[0]), ("\"total_w\": 1", x[1]), ("\"watts\": 1", x[1])] {
        assert_eq!(text.matches(needle).count(), 1);
        text = text.replace(needle, &format!("{}: {value}", needle.split(':').next().unwrap()));
    }
    let path = dir.join("physical.json"); fs::write(&path, text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(path).output().unwrap();
    number(document(&output).get("objective").unwrap(), "value_k")
}

#[test]
fn mixed_physical_cooling_matches_independent_reconstruction_and_parallel_output() {
    let dir = scratch("physics"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(300.0)).unwrap();
    let output = command(&base, &input, None).output().unwrap();
    let actual = document(&output);
    assert_eq!(actual.str_field("schema"), Some("frankensim.cooling-network-uq.shapley.v1"));
    assert_eq!(actual.str_field("method"), Some("conditional-product-philox"));
    assert_eq!(actual.str_field("correlation"), Some("independent"));
    assert_eq!(number(&actual, "samples_evaluated"), 12.0);
    let mut plan = UqPlan::new("reference", PropagationMethod::MonteCarlo, 12)
        .with_parameter(ParameterUncertainty::gaussian("inlet[0].temperature_k", 300.0, 0.2, "K"))
        .with_parameter(ParameterUncertainty::uniform("component[chip].power_w", 0.9, 1.1, "W"))
        .with_correlation(CorrelationModel::Independent);
    plan.seed = 73;
    let original = fs::read_to_string(&base).unwrap();
    let mut run = ShapleyExecution::new(&plan).unwrap();
    let report = run.advance(12, || false, |x| Ok::<_, &str>(physical(&dir, &original, x)));
    assert_eq!(report.status, UqStatus::Complete);
    let expected = report.estimate.unwrap();
    for (effect, reference) in actual.get("effects").unwrap().as_array().unwrap().iter().zip(&expected.effects) {
        assert_eq!(effect.str_field("parameter"), Some(reference.parameter.as_str()));
        assert_eq!(effect.str_field("unit"), Some(reference.unit.as_str()));
        assert!((number(effect, "index") - reference.index).abs() < 1e-9);
    }
    for (fraction, reference) in actual.get("conditional_variance_fractions").unwrap()
        .as_array().unwrap().iter().zip(&expected.conditional_variance_fractions)
    { assert!((fraction.as_f64().unwrap() - reference).abs() < 1e-9); }
    let parallel = command(&base, &input, Some("3")).output().unwrap(); document(&parallel);
    assert_eq!(output.stdout, parallel.stdout);
    for field in ["confidence_interval", "probability_of_compliance", "mean_k"] {
        assert!(actual.get(field).is_none());
    }
}

#[test]
fn mixed_serial_and_parallel_checkpoints_resume_the_same_physical_experiment() {
    let dir = scratch("recovery"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(300.0)).unwrap();
    for parallel in [false, true] {
        let workers = parallel.then_some("3"); let suffix = if parallel { "batch" } else { "serial" };
        let full_path = dir.join(format!("full.{suffix}"));
        let full = command(&base, &input, workers).arg("--checkpoint").arg(&full_path).output().unwrap();
        document(&full);
        let saved = dir.join(format!("prefix.{suffix}"));
        let partial = command(&base, &input, workers).arg("--checkpoint").arg(&saved)
            .args(["--max-new-samples", "3"]).output().unwrap();
        assert_eq!(partial.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
        let progress = J::parse(std::str::from_utf8(&partial.stdout).unwrap()).unwrap();
        assert_eq!(number(&progress, "samples_evaluated"), 3.0);
        assert!(progress.get("effects").is_none());
        let before = fs::read(&saved).unwrap();
        let completed = dir.join(format!("completed.{suffix}"));
        let resumed = command(&base, &input, parallel.then_some("1"))
            .arg("--resume").arg(&saved).arg("--checkpoint").arg(&completed).output().unwrap();
        document(&resumed); assert_eq!(resumed.stdout, full.stdout);
        assert_eq!(fs::read(&saved).unwrap(), before);
        assert_eq!(fs::read(&completed).unwrap(), fs::read(&full_path).unwrap());
        let fast = dir.join(format!("fast-{suffix}.json")); fs::write(&fast, request(1e-9)).unwrap();
        let terminal = command(&base, &fast, workers).arg("--resume").arg(&completed).output().unwrap();
        document(&terminal); assert_eq!(terminal.stdout, full.stdout);
    }
}

#[test]
fn uniform_initial_temperature_retains_complete_transient_observations() {
    let dir = scratch("transient"); let input = dir.join("uq.json");
    fs::write(&input, r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,"wall_seconds":300,
        "correlation":{"kind":"independent"},"qoi":{"kind":"transient-sampled-peak"},
        "parameters":[{"target":{"kind":"initial-temperature"},"distribution":{"kind":"uniform","lo":299.8,"hi":300.2}}]}"#).unwrap();
    let base = example("nonlinear-contact-pulse.json");
    let serial = command(&base, &input, None).output().unwrap(); let report = document(&serial);
    assert_eq!(report.str_field("method"), Some("conditional-product-philox"));
    let qoi = report.get("qoi").unwrap();
    assert_eq!(qoi.str_field("kind"), Some("transient-sampled-peak"));
    assert_eq!(qoi.str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(number(&report.get("effects").unwrap().as_array().unwrap()[0], "index"), 1.0);
    let parallel = command(&base, &input, Some("2")).output().unwrap(); document(&parallel);
    assert_eq!(parallel.stdout, serial.stdout);
}

#[test]
fn non_gaussian_dependence_is_not_silently_replaced_with_independence() {
    let dir = scratch("admission"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); let destination = dir.join("must-not-exist.shapley");
    for text in [request(300.0).replace("\"independent\"", "\"unknown\""),
        request(300.0).replace("\"kind\":\"independent\"", "\"kind\":\"joint-gaussian\",\"matrix\":[[1,0],[0,1]]"),
        request(300.0).replace("\"hi\":1.1", "\"hi\":0.9")]
    {
        fs::write(&input, text).unwrap();
        for workers in [None, Some("2")] {
            let output = command(&base, &input, workers).arg("--checkpoint").arg(&destination).output().unwrap();
            assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
            assert!(output.stdout.is_empty()); assert!(!destination.exists());
        }
    }
}
