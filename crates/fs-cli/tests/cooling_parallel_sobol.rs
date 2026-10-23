//! Real-binary Sobol sensitivity: original physical solves and original reports.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch(name: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-parallel-sobol-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn request(transient: bool, wall: f64) -> String {
    if transient {
        format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":12,"wall_seconds":{wall},"correlation":{{"kind":"independent"}},"qoi":{{"kind":"transient-sampled-peak"}},"parameters":[{{"target":{{"kind":"initial-temperature"}},"distribution":{{"kind":"gaussian","mean":300,"std_dev":0.2}}}}]}}"#)
    } else {
        format!(r#"{{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":16,"wall_seconds":{wall},"correlation":{{"kind":"independent"}},"parameters":[{{"target":{{"kind":"inlet-temperature","index":0}},"distribution":{{"kind":"uniform","lo":299,"hi":301}}}},{{"target":{{"kind":"component-power","component":"chip"}},"distribution":{{"kind":"gaussian","mean":1,"std_dev":0.01}}}}]}}"#)
    }
}
fn command(base: &Path, input: &Path, workers: Option<&str>) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    cmd.args(["--json", "cooling-network-uq"]).arg(base).arg(input).args(["--sensitivity", "sobol"]);
    if let Some(n) = workers { cmd.args(["--workers", n]); }
    cmd
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let value = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(value.str_field("schema"), Some("frankensim.cooling-network-uq.sensitivity.parallel.progress.v1"));
    for name in ["effects", "first_order", "total_order", "base_output_std_dev_k"] {
        assert!(value.get(name).is_none());
    }
    value
}
fn number(value: &J, key: &str) -> f64 { value.get(key).unwrap().as_f64().unwrap() }

#[test]
fn parallel_steady_and_complete_transient_indices_match_serial_result_bytes() {
    let dir = scratch("parity");
    for (transient, name) in [(false, "fan-correlated-hotspot.json"), (true, "nonlinear-contact-pulse.json")] {
        let input = dir.join("uq.json"); fs::write(&input, request(transient, 300.0)).unwrap();
        let base = example(name);
        let serial = command(&base, &input, None).output().unwrap();
        assert_eq!(parsed(&serial).str_field("schema"), Some("frankensim.cooling-network-uq.sensitivity.v1"));
        for workers in ["1", "2", "5"] {
            let parallel = command(&base, &input, Some(workers)).output().unwrap();
            let result = parsed(&parallel);
            assert_eq!(parallel.stdout, serial.stdout, "workers={workers}, transient={transient}");
            if transient {
                assert_eq!(result.get("qoi").unwrap().str_field("observation"), Some("one-completed-trajectory"));
            }
        }
    }
}

#[test]
fn partial_hybrid_rows_resume_with_changed_workers_and_identical_final_checkpoint() {
    let dir = scratch("resume"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(false, 300.0)).unwrap();
    let full_path = dir.join("full.sbatch");
    let full = command(&base, &input, Some("5")).arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first_path = dir.join("first.sbatch");
    let first = command(&base, &input, Some("3")).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let p = progress(&first);
    assert_eq!(number(&p, "samples_evaluated"), 3.0);
    assert_eq!(p.str_field("next_evaluation_kind"), Some("hybrid"));
    assert_eq!(number(&p, "next_row_slot"), 3.0);
    let original = fs::read(&first_path).unwrap(); assert_eq!(&original[..8], b"FSSOBB01");
    let second_path = dir.join("second.sbatch");
    let second = command(&base, &input, Some("2")).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "4"]).output().unwrap();
    let p = progress(&second);
    assert_eq!(number(&p, "samples_evaluated"), 7.0);
    assert_eq!(number(&p, "samples_evaluated_this_run"), 4.0);
    assert_eq!(fs::read(&first_path).unwrap(), original);
    let final_path = dir.join("final.sbatch");
    let done = command(&base, &input, Some("1")).arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&final_path).output().unwrap();
    parsed(&done); assert_eq!(done.stdout, full.stdout);
    assert_eq!(fs::read(&final_path).unwrap(), fs::read(&full_path).unwrap());
    fs::write(&input, request(false, 1e-9)).unwrap();
    let replay = command(&base, &input, Some("4")).arg("--resume").arg(&final_path).output().unwrap();
    parsed(&replay); assert_eq!(replay.stdout, full.stdout);
}

#[test]
fn invalid_laws_layouts_and_changed_identities_refuse_before_output_reservation() {
    let dir = scratch("refusals"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); let unused = dir.join("never.sbatch");
    let original = request(false, 300.0);
    for invalid in [original.replace("\"samples\":16", "\"samples\":15"),
        original.replace("\"independent\"", "\"unknown\""),
        original.replace("\"std_dev\":0.01", "\"std_dev\":0")]
    {
        fs::write(&input, invalid).unwrap();
        let out = command(&base, &input, Some("2")).arg("--checkpoint").arg(&unused).output().unwrap();
        assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
        assert!(out.stdout.is_empty()); assert!(!unused.exists());
    }
    fs::write(&input, &original).unwrap();
    let saved = dir.join("empty.sbatch");
    let out = command(&base, &input, Some("2")).arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(number(&progress(&out), "samples_evaluated"), 0.0);
    let bytes = fs::read(&saved).unwrap();
    fs::write(&input, original.replace("\"seed\":\"73\"", "\"seed\":\"74\"")).unwrap();
    let out = command(&base, &input, Some("3")).arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert!(!out.status.success()); assert!(!unused.exists());
    fs::write(&input, &original).unwrap();
    let out = command(&base, &input, Some("2")).arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&saved).output().unwrap();
    assert!(!out.status.success()); assert_eq!(fs::read(&saved).unwrap(), bytes);
    let out = command(&base, &input, None).arg("--resume").arg(&saved).output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    let missing = dir.join("missing.json");
    let out = command(&missing, &missing, Some("2")).args(["--qmc-replicates", "2"])
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::USAGE))); assert!(!unused.exists());
}

#[test]
fn timeout_retains_progress_but_solver_failure_is_terminal_after_workers_drain() {
    let dir = scratch("failure"); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(false, 1e-9)).unwrap();
    let saved = dir.join("paused.sbatch");
    let paused = command(&base, &input, Some("3")).arg("--checkpoint").arg(&saved).output().unwrap();
    assert_eq!(progress(&paused).str_field("termination"), Some("wall-time-budget"));
    fs::write(&input, request(false, 300.0)).unwrap();
    parsed(&command(&base, &input, Some("2")).arg("--resume").arg(&saved).output().unwrap());
    let original = fs::read_to_string(&base).unwrap();
    assert!(original.contains("\"linear_iterations\": 20000"));
    let broken = dir.join("broken.json");
    fs::write(&broken, original.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.sbatch");
    let out = command(&broken, &input, Some("3")).arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(out.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let next = dir.join("never.sbatch");
    let out = command(&broken, &input, Some("1")).arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!out.status.success()); assert!(!next.exists());
}
