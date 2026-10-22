//! Real cooling child evaluations; reference statistics use the existing UQ owner.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use fs_uq::{CorrelationModel, ImportanceConfig, ImportanceExecution, ParameterUncertainty,
    PropagationMethod, UqPlan, UqStatus};
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-cooling-importance-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path, proposal: &str) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request)
        .args(["--importance", proposal]); c
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn number(root: &J, key: &str) -> f64 { root.get(key).unwrap().as_f64().unwrap() }
fn request(n: usize, wall: f64, threshold: f64, transient: bool) -> String {
    let target = if transient { r#"{"kind":"initial-temperature"}"# }
        else { r#"{"kind":"inlet-temperature","index":0}"# };
    let observable = if transient { r#", "qoi":{"kind":"transient-sampled-peak"}"# } else { "" };
    format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":{n},"wall_seconds":{wall},"temperature_limit_k":{threshold},"correlation":{{"kind":"independent"}},"parameters":[{{"target":{target},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0.25}}}}]{observable}}}"#)
}
fn physical(dir: &Path, original: &str, inlet: f64) -> f64 {
    assert_eq!(original.matches("\"temperature_k\": 300").count(), 1);
    let path = dir.join("physical.json");
    fs::write(&path, original.replace("\"temperature_k\": 300", &format!("\"temperature_k\": {inlet}"))).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).arg(path).output().unwrap();
    number(parsed(&output).get("objective").unwrap(), "value_k")
}

#[test]
fn weighted_probability_matches_independent_run_of_actual_cooling_calls() {
    let dir = scratch("physical");
    let base = example("fan-correlated-hotspot.json");
    let original = fs::read_to_string(&base).unwrap();
    let threshold = physical(&dir, &original, 300.0) + 0.5;
    let path = dir.join("uq.json"); fs::write(&path, request(16, 300.0, threshold, false)).unwrap();
    let output = command(&base, &path, "0.25,2").output().unwrap();
    let result = parsed(&output);
    assert_eq!(result.str_field("schema"), Some("frankensim.cooling-network-uq.importance.v1"));
    assert_eq!(number(&result, "samples_evaluated"), 16.0);
    let mut plan = UqPlan::new("temperature", PropagationMethod::MonteCarlo, 16)
        .with_parameter(ParameterUncertainty::gaussian("inlet", 300.0, 0.25, "K"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(threshold);
    plan.seed = 73;
    let mut reference = ImportanceExecution::new(&plan,
        ImportanceConfig { germ_shift: vec![2.0], defensive_fraction: 0.25 }).unwrap();
    let report = reference.advance(16, || false,
        |x| Ok::<_, &str>(physical(&dir, &original, x[0])));
    assert_eq!(report.status, UqStatus::Complete);
    let expected = report.estimate.unwrap();
    assert!((number(&result, "estimated_failure_probability") - expected.failure_probability).abs() < 1e-13);
    assert!((number(&result, "mean_likelihood_weight") - expected.mean_likelihood_weight).abs() < 1e-13);
    assert_eq!(result.get("failure_standard_error").unwrap().as_f64(), expected.standard_error);
    assert_eq!(number(&result, "proposal_failure_count"), report.proposal_failures as f64);
    for name in ["mean_k", "percentiles_p05_p50_p95_k", "probability_of_compliance", "confidence_interval"] {
        assert!(result.get(name).is_none(), "proposal samples published as {name}");
    }
    let repeat = command(&base, &path, "0.25,2").output().unwrap(); parsed(&repeat);
    assert_eq!(repeat.stdout, output.stdout);
}

#[test]
fn original_target_proposal_agrees_with_raw_mc_failure_frequency() {
    let dir = scratch("mc");
    let base = example("fan-correlated-hotspot.json");
    let path = dir.join("uq.json"); fs::write(&path, request(8, 300.0, 303.0, false)).unwrap();
    let ordinary = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network-uq"]).arg(&base).arg(&path).output().unwrap();
    let mc = parsed(&ordinary);
    for proposal in ["0.2,0", "1,5"] {
        let report = parsed(&command(&base, &path, proposal).output().unwrap());
        assert_eq!(number(&report, "estimated_failure_probability"),
            1.0 - number(&mc, "empirical_probability_of_compliance"));
        assert_eq!(number(&report, "mean_likelihood_weight"), 1.0);
    }
}

#[test]
fn incompatible_laws_modes_timeouts_and_failed_physics_publish_no_probability() {
    let dir = scratch("refusals");
    let missing = dir.join("missing.json");
    let unused = dir.join("unused.checkpoint");
    let output = command(&missing, &missing, "0.2,2").args(["--qmc-replicates", "2"]).arg("--checkpoint").arg(&unused).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::USAGE)));
    assert!(!unused.exists()); assert!(output.stdout.is_empty());
    let base = example("fan-correlated-hotspot.json");
    let path = dir.join("uq.json");
    let valid = request(8, 300.0, 303.0, false);
    for invalid in [
        valid.replace("\"gaussian\",\"mean\":300,\"std_dev\":0.25", "\"uniform\",\"lo\":299,\"hi\":301"),
        valid.replace("\"independent\"", "\"unknown\""),
        valid.replace("\"temperature_limit_k\":303,", ""),
    ] {
        fs::write(&path, invalid).unwrap();
        let output = command(&base, &path, "0.2,2").output().unwrap();
        assert!(!output.status.success()); assert!(output.stdout.is_empty());
    }
    fs::write(&path, request(8, 1e-9, 303.0, false)).unwrap();
    let output = command(&base, &path, "0.2,2").output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET))); assert!(output.stdout.is_empty());
    fs::write(&path, valid).unwrap();
    let wrong_dimension = command(&base, &path, "0.2,1,1").output().unwrap();
    assert!(!wrong_dimension.status.success()); assert!(wrong_dimension.stdout.is_empty());
    let broken = dir.join("broken.json");
    let original = fs::read_to_string(&base).unwrap();
    assert!(original.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, original.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let output = command(&broken, &path, "0.2,2").output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(output.stdout.is_empty());
}

#[test]
fn transient_failure_observations_are_complete_trajectory_peaks() {
    let dir = scratch("transient");
    let path = dir.join("uq.json"); fs::write(&path, request(4, 300.0, 303.0, true)).unwrap();
    let output = command(&example("nonlinear-contact-pulse.json"), &path, "0.25,2").output().unwrap();
    let result = parsed(&output);
    let qoi = result.get("qoi").unwrap();
    assert_eq!(qoi.str_field("kind"), Some("transient-sampled-peak"));
    assert_eq!(qoi.str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(qoi.str_field("temporal_scope"), Some("initial-state-and-accepted-endpoints"));
    assert_eq!(number(&result, "samples_evaluated"), 4.0);
    assert_eq!(result.str_field("failure_event"), Some("qoi > temperature_limit_k"));
}

fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let row = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(row.str_field("schema"), Some("frankensim.cooling-network-uq.importance.progress.v1"));
    for name in ["estimated_failure_probability", "failure_standard_error", "failure_effective_samples", "proposal_failure_count"] {
        assert!(row.get(name).is_none(), "published partial statistic {name}");
    }
    row
}

#[test]
fn durable_importance_keeps_paid_calls_and_final_report_bytes() {
    let dir = scratch("recovery"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(8, 300.0, 303.0, false)).unwrap();
    let full_path = dir.join("full.imp");
    let full = command(&base, &input, "0.25,2").arg("--checkpoint").arg(&full_path).output().unwrap(); parsed(&full);
    let first_path = dir.join("first.imp");
    let first = command(&base, &input, "0.25,2").arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    assert_eq!(number(&progress(&first), "samples_evaluated"), 3.0);
    let retained = fs::read(&first_path).unwrap(); assert_eq!(&retained[..8], b"FSIMP001");
    let second_path = dir.join("second.imp");
    let second = command(&base, &input, "0.25,2").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "2"]).output().unwrap();
    let row = progress(&second); assert_eq!(number(&row, "samples_evaluated"), 5.0);
    assert_eq!(number(&row, "samples_evaluated_this_run"), 2.0);
    assert_eq!(fs::read(&first_path).unwrap(), retained);
    let done_path = dir.join("done.imp");
    let done = command(&base, &input, "0.25,2").arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&done_path).output().unwrap(); parsed(&done);
    assert_eq!(done.stdout, full.stdout); assert_eq!(fs::read(&done_path).unwrap(), fs::read(&full_path).unwrap());
    let plain = command(&base, &input, "0.25,2").output().unwrap(); parsed(&plain);
    assert_eq!(plain.stdout, full.stdout);
    fs::write(&input, request(8, 1e-9, 303.0, false)).unwrap();
    let terminal = command(&base, &input, "0.25,2").arg("--resume").arg(&done_path).output().unwrap(); parsed(&terminal);
    assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn changed_proposals_corruption_and_existing_destinations_cannot_reweight_saved_work() {
    let dir = scratch("binding"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(8, 300.0, 303.0, false)).unwrap();
    let saved = dir.join("original.imp");
    let first = command(&base, &input, "0.25,2").arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "0"]).output().unwrap(); progress(&first);
    let bytes = fs::read(&saved).unwrap();
    for (i, proposal) in ["0.3,2", "0.25,3"].iter().enumerate() {
        let next = dir.join(format!("changed-{i}.imp"));
        let result = command(&base, &input, proposal).arg("--resume").arg(&saved)
            .arg("--checkpoint").arg(&next).output().unwrap();
        assert!(!result.status.success()); assert!(result.stdout.is_empty()); assert!(!next.exists());
    }
    let corrupt = dir.join("corrupt.imp"); let mut damaged = bytes.clone(); damaged[10] ^= 1;
    fs::write(&corrupt, damaged).unwrap(); let next = dir.join("corrupt-next.imp");
    let result = command(&base, &input, "0.25,2").arg("--resume").arg(&corrupt)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!result.status.success()); assert!(!next.exists());
    let result = command(&base, &input, "0.25,2").arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&saved).output().unwrap();
    assert!(!result.status.success()); assert_eq!(fs::read(&saved).unwrap(), bytes);
}

#[test]
fn timeout_retains_recoverable_progress_but_physical_failure_invalidates_it() {
    let dir = scratch("failure-recovery"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(8, 1e-9, 303.0, false)).unwrap();
    let saved = dir.join("timeout.imp");
    let stopped = command(&base, &input, "0.25,2").arg("--checkpoint").arg(&saved).output().unwrap();
    let row = progress(&stopped); assert_eq!(number(&row, "samples_evaluated"), 0.0);
    assert_eq!(row.str_field("termination"), Some("wall-time-budget"));
    fs::write(&input, request(8, 300.0, 303.0, false)).unwrap();
    let resumed = command(&base, &input, "0.25,2").arg("--resume").arg(&saved).output().unwrap();
    assert_eq!(number(&parsed(&resumed), "samples_evaluated"), 8.0);
    let broken = dir.join("broken.json");
    fs::write(&broken, fs::read_to_string(&base).unwrap().replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.imp");
    let output = command(&broken, &input, "0.25,2").arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(output.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let next = dir.join("must-not-exist.imp");
    let retry = command(&broken, &input, "0.25,2").arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!retry.status.success()); assert!(!next.exists());
}

#[test]
fn completed_transient_observations_survive_process_restart() {
    let dir = scratch("transient-recovery"); let base = example("nonlinear-contact-pulse.json");
    let input = dir.join("uq.json"); fs::write(&input, request(4, 300.0, 303.0, true)).unwrap();
    let full = command(&base, &input, "0.25,2").output().unwrap(); parsed(&full);
    let saved = dir.join("transient.imp");
    let partial = command(&base, &input, "0.25,2").arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "1"]).output().unwrap();
    let row = progress(&partial); assert_eq!(number(&row, "samples_evaluated"), 1.0);
    assert_eq!(row.get("qoi").unwrap().str_field("observation"), Some("one-completed-trajectory"));
    let resumed = command(&base, &input, "0.25,2").arg("--resume").arg(&saved).output().unwrap(); parsed(&resumed);
    assert_eq!(resumed.stdout, full.stdout);
}
