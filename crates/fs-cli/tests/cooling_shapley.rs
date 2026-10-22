//! Real cooling binary, correlated inputs and restartable conditional pairs.
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
    let path = std::env::temp_dir().join(format!("fs-cooling-shapley-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request).args(["--sensitivity", "shapley"]); c
}
fn document(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn completed(output: &Output) -> J {
    let doc = document(output);
    assert_eq!(doc.str_field("schema"), Some("frankensim.cooling-network-uq.shapley.v1")); doc
}
fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)),
        "{}", String::from_utf8_lossy(&output.stderr));
    let doc = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(doc.str_field("schema"), Some("frankensim.cooling-network-uq.shapley.progress.v1"));
    for key in ["effects", "conditional_variance_fractions", "base_output_std_dev_k", "mean_k", "confidence_interval"] {
        assert!(doc.get(key).is_none(), "partial design published {key}");
    }
    doc
}
fn number(doc: &J, key: &str) -> f64 { doc.get(key).unwrap().as_f64().unwrap() }
fn request(rho: f64, samples: usize, wall: f64) -> String {
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{samples},"wall_seconds":{wall},"temperature_limit_k":303,
        "correlation":{{"kind":"joint-gaussian","matrix":[[1,{rho}],[{rho},1]]}},
        "parameters":[{{"target":{{"kind":"inlet-temperature","index":0}},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0.2}}}},
        {{"target":{{"kind":"component-power","component":"chip"}},"distribution":{{"kind":"gaussian","mean":1,"std_dev":0.02}}}}]}}"#)
}
fn physical(dir: &Path, original: &str, x: &[f64]) -> f64 {
    let mut text = original.to_string();
    for (needle, value) in [("\"temperature_k\": 300", x[0]), ("\"total_w\": 1", x[1]), ("\"watts\": 1", x[1])] {
        assert_eq!(text.matches(needle).count(), 1);
        let key = needle.split(':').next().unwrap();
        text = text.replace(needle, &format!("{key}: {value}"));
    }
    let path = dir.join("physical.json"); fs::write(&path, text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(path).output().unwrap();
    number(document(&output).get("objective").unwrap(), "value_k")
}

#[test]
fn correlated_and_singular_effects_match_independently_executed_physics() {
    let dir = scratch("physical");
    let base = example("fan-correlated-hotspot.json");
    let original = fs::read_to_string(&base).unwrap();
    let input = dir.join("uq.json");
    for rho in [0.6, 1.0] {
        fs::write(&input, request(rho, 12, 300.0)).unwrap();
        let output = command(&base, &input).output().unwrap();
        let actual = completed(&output);
        assert_eq!(number(&actual, "samples_evaluated"), 12.0);
        assert_eq!(number(&actual, "base_samples"), 2.0);
        assert_eq!(number(&actual, "evaluations_per_row"), 6.0);
        assert_eq!(actual.str_field("correlation"), Some("joint-gaussian"));
        let mut plan = UqPlan::new("reference", PropagationMethod::MonteCarlo, 12)
            .with_parameter(ParameterUncertainty::gaussian("inlet[0].temperature_k", 300.0, 0.2, "K"))
            .with_parameter(ParameterUncertainty::gaussian("component[chip].power_w", 1.0, 0.02, "W"))
            .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, rho], vec![rho, 1.0]] });
        plan.seed = 73;
        let mut run = ShapleyExecution::new(&plan).unwrap();
        let expected = run.advance(12, || false, |x| Ok::<_, &str>(physical(&dir, &original, x)));
        assert_eq!(expected.status, UqStatus::Complete);
        let expected = expected.estimate.unwrap();
        let effects = actual.get("effects").unwrap().as_array().unwrap();
        assert_eq!(effects.len(), 2);
        for (effect, reference) in effects.iter().zip(&expected.effects) {
            assert_eq!(effect.str_field("parameter"), Some(reference.parameter.as_str()));
            assert_eq!(effect.str_field("unit"), Some(reference.unit.as_str()));
            assert!((number(effect, "index") - reference.index).abs() < 1e-9);
            assert!(effect.get("first_order").is_none() && effect.get("total_order").is_none());
            if rho == 1.0 { assert!((number(effect, "index") - 0.5).abs() < 1e-12); }
        }
        let fractions = actual.get("conditional_variance_fractions").unwrap().as_array().unwrap();
        assert_eq!(fractions.len(), 4);
        for (got, want) in fractions.iter().zip(&expected.conditional_variance_fractions) {
            assert!((got.as_f64().unwrap() - want).abs() < 1e-9);
        }
        for key in ["probability_of_compliance", "confidence_interval", "mean_k"] { assert!(actual.get(key).is_none()); }
        let again = command(&base, &input).output().unwrap(); completed(&again);
        assert_eq!(output.stdout, again.stdout);
    }
}

#[test]
fn physical_transient_observations_are_complete_trajectory_peaks() {
    let dir = scratch("transient"); let input = dir.join("uq.json");
    fs::write(&input, r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":8,"wall_seconds":300,
        "correlation":{"kind":"independent"},"qoi":{"kind":"transient-sampled-peak"},
        "parameters":[{"target":{"kind":"initial-temperature"},"distribution":{"kind":"gaussian","mean":300,"std_dev":0.2}}]}"#).unwrap();
    let result = completed(&command(&example("nonlinear-contact-pulse.json"), &input).output().unwrap());
    let qoi = result.get("qoi").unwrap();
    assert_eq!(qoi.str_field("kind"), Some("transient-sampled-peak"));
    assert_eq!(qoi.str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(number(&result, "samples_evaluated"), 8.0);
    let effects = result.get("effects").unwrap().as_array().unwrap();
    assert_eq!(number(&effects[0], "index"), 1.0);
}

#[test]
fn recovery_retains_paid_halves_and_matches_uninterrupted_report_and_bytes() {
    let dir = scratch("recovery"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(0.6, 12, 300.0)).unwrap();
    let full_file = dir.join("full.shapley");
    let full = command(&base, &input).arg("--checkpoint").arg(&full_file).output().unwrap(); completed(&full);
    let first_file = dir.join("first.shapley");
    let first = command(&base, &input).arg("--checkpoint").arg(&first_file)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let row = progress(&first);
    assert_eq!(number(&row, "samples_evaluated"), 3.0);
    assert_eq!(number(&row, "next_subset_mask"), 1.0);
    assert_eq!(number(&row, "next_pair_half"), 1.0);
    let saved = fs::read(&first_file).unwrap(); assert_eq!(&saved[..8], b"FSSHP001");
    let second_file = dir.join("second.shapley");
    let second = command(&base, &input).arg("--resume").arg(&first_file)
        .arg("--checkpoint").arg(&second_file).args(["--max-new-samples", "4"]).output().unwrap();
    let row = progress(&second);
    assert_eq!(number(&row, "samples_evaluated_this_run"), 4.0);
    assert_eq!(number(&row, "next_row_ordinal"), 1.0);
    assert_eq!(number(&row, "next_subset_mask"), 0.0);
    assert_eq!(number(&row, "next_pair_half"), 1.0);
    assert_eq!(fs::read(&first_file).unwrap(), saved);
    let final_file = dir.join("final.shapley");
    let final_run = command(&base, &input).arg("--resume").arg(&second_file)
        .arg("--checkpoint").arg(&final_file).output().unwrap(); completed(&final_run);
    assert_eq!(final_run.stdout, full.stdout);
    assert_eq!(fs::read(&final_file).unwrap(), fs::read(&full_file).unwrap());
    let fast = dir.join("fast.json"); fs::write(&fast, request(0.6, 12, 1e-9)).unwrap();
    let terminal = command(&base, &fast).arg("--resume").arg(&final_file).output().unwrap(); completed(&terminal);
    assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn invalid_laws_allocations_and_conflicting_options_do_not_create_outputs() {
    let dir = scratch("admission"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); let destination = dir.join("must-not-exist.shapley");
    for text in [request(0.6, 13, 300.0), request(0.6, 12, 300.0).replace("\"std_dev\":0.2", "\"std_dev\":0"),
        request(0.6, 12, 300.0).replace("joint-gaussian", "unknown"),
        fs::read_to_string(example("uq-fan-hotspot.json")).unwrap()]
    {
        fs::write(&input, text).unwrap();
        let output = command(&base, &input).arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success()); assert!(output.stdout.is_empty()); assert!(!destination.exists());
    }
    let missing = dir.join("missing.json");
    for tail in [vec!["--sensitivity", "sobol"], vec!["--pce", "1,6"], vec!["--pce-mean", "0,2,4"],
        vec!["--mean-control", "adjoint"], vec!["--qmc-replicates", "2"]]
    {
        let output = command(&missing, &missing).args(tail).arg("--checkpoint").arg(&destination).output().unwrap();
        assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
        assert!(output.stdout.is_empty()); assert!(!destination.exists());
    }
}

#[test]
fn timeouts_retain_only_progress_and_failed_physics_invalidates_recovery() {
    let dir = scratch("failures"); let base = example("fan-correlated-hotspot.json");
    let fast = dir.join("fast.json"); fs::write(&fast, request(0.6, 12, 1e-9)).unwrap();
    let output = command(&base, &fast).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET))); assert!(output.stdout.is_empty());
    let prefix = dir.join("timeout.shapley");
    let output = command(&base, &fast).arg("--checkpoint").arg(&prefix).output().unwrap();
    let row = progress(&output); assert_eq!(number(&row, "samples_evaluated"), 0.0);
    assert_eq!(row.str_field("termination"), Some("wall-time-budget"));
    let valid = dir.join("valid.json"); fs::write(&valid, request(0.6, 12, 300.0)).unwrap();
    completed(&command(&base, &valid).arg("--resume").arg(&prefix).output().unwrap());
    let original = fs::read_to_string(&base).unwrap(); assert!(original.contains("\"linear_iterations\": 20000"));
    let broken = dir.join("broken.json");
    fs::write(&broken, original.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.shapley");
    let output = command(&broken, &valid).arg("--checkpoint").arg(&failed).output().unwrap();
    assert!(!output.status.success()); assert!(output.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let next = dir.join("must-not-resume.shapley");
    let output = command(&broken, &valid).arg("--resume").arg(&failed).arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!output.status.success()); assert!(!next.exists());
}

#[test]
fn changed_binding_and_corruption_refuse_before_reserving_a_destination() {
    let dir = scratch("identity"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(0.6, 12, 300.0)).unwrap();
    let prefix = dir.join("empty.shapley");
    progress(&command(&base, &input).arg("--checkpoint").arg(&prefix).args(["--max-new-samples", "0"]).output().unwrap());
    let saved = fs::read(&prefix).unwrap();
    let changed = dir.join("different-law.json"); fs::write(&changed, request(0.5, 12, 300.0)).unwrap();
    let next = dir.join("must-not-exist.shapley");
    let output = command(&base, &changed).arg("--resume").arg(&prefix).arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!output.status.success()); assert!(output.stdout.is_empty()); assert!(!next.exists());
    let corrupt = dir.join("corrupt.shapley"); let mut bytes = saved.clone(); *bytes.last_mut().unwrap() ^= 1;
    fs::write(&corrupt, bytes).unwrap();
    let output = command(&base, &input).arg("--resume").arg(&corrupt).arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!output.status.success()); assert!(!next.exists());
    let output = command(&base, &input).arg("--resume").arg(&prefix).arg("--checkpoint").arg(&prefix).output().unwrap();
    assert!(!output.status.success()); assert_eq!(fs::read(&prefix).unwrap(), saved);
}
