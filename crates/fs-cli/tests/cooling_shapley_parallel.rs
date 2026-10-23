//! Actual cooling evaluations under the same serial/parallel conditional law.
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}
fn scratch() -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-shapley-parallel-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap(); path
}
fn command(base: &Path, request: &Path, workers: Option<&str>) -> Command {
    let mut c = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    c.args(["--json", "cooling-network-uq"]).arg(base).arg(request).args(["--sensitivity", "shapley"]);
    if let Some(n) = workers { c.args(["--workers", n]); } c
}
fn parsed(output: &Output) -> J {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let j = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(j.str_field("schema"), Some("frankensim.cooling-network-uq.shapley.parallel.progress.v1"));
    for field in ["effects", "conditional_variance_fractions", "base_output_std_dev_k", "mean_k"] {
        assert!(j.get(field).is_none());
    }
    j
}
fn number(j: &J, field: &str) -> f64 { j.get(field).unwrap().as_f64().unwrap() }
fn request() -> String {
    let text = fs::read_to_string(example("uq-shapley-hotspot.json")).unwrap();
    assert!(text.contains("\"samples\": 24")); text.replace("\"samples\": 24", "\"samples\": 12")
}

#[test]
fn actual_correlated_singular_and_transient_results_match_serial_bytes() {
    let dir = scratch(); let input = dir.join("uq.json");
    let base = example("fan-correlated-hotspot.json");
    let original = request();
    for text in [original.clone(), original.replace("[[1, 0.6], [0.6, 1]]", "[[1, 1], [1, 1]]")] {
        fs::write(&input, text).unwrap();
        let serial = command(&base, &input, None).output().unwrap();
        assert_eq!(parsed(&serial).str_field("schema"), Some("frankensim.cooling-network-uq.shapley.v1"));
        for n in ["1", "5"] {
            let parallel = command(&base, &input, Some(n)).output().unwrap(); parsed(&parallel);
            assert_eq!(parallel.stdout, serial.stdout, "workers={n}");
        }
    }
    fs::write(&input, r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":4,
        "wall_seconds":300,"qoi":{"kind":"transient-sampled-peak"},"correlation":{"kind":"independent"},
        "parameters":[{"target":{"kind":"initial-temperature"},
        "distribution":{"kind":"gaussian","mean":300,"std_dev":0.2}}]}"#).unwrap();
    let base = example("nonlinear-contact-pulse.json");
    let serial = command(&base, &input, None).output().unwrap(); parsed(&serial);
    let parallel = command(&base, &input, Some("3")).output().unwrap();
    let result = parsed(&parallel);
    assert_eq!(result.get("qoi").unwrap().str_field("observation"), Some("one-completed-trajectory"));
    assert_eq!(parallel.stdout, serial.stdout);
}

#[test]
fn half_pair_and_row_boundary_restarts_keep_the_original_result_and_budget() {
    let dir = scratch(); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request()).unwrap();
    let full_path = dir.join("full.shapley");
    let full = command(&base, &input, Some("5")).arg("--checkpoint").arg(&full_path).output().unwrap(); parsed(&full);
    let first_path = dir.join("first.shapley");
    let first = command(&base, &input, Some("3")).arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let p = progress(&first);
    assert_eq!(number(&p, "samples_evaluated"), 3.0);
    assert_eq!(number(&p, "next_subset_mask"), 1.0);
    assert_eq!(number(&p, "next_pair_half"), 1.0);
    let bytes = fs::read(&first_path).unwrap(); assert_eq!(&bytes[..8], b"FSSHPB01");
    let second_path = dir.join("second.shapley");
    let second = command(&base, &input, Some("2")).arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "4"]).output().unwrap();
    let p = progress(&second);
    assert_eq!(number(&p, "samples_evaluated"), 7.0);
    assert_eq!(number(&p, "samples_evaluated_this_run"), 4.0);
    assert_eq!(number(&p, "next_row_ordinal"), 1.0);
    assert_eq!(number(&p, "next_pair_half"), 1.0);
    assert_eq!(fs::read(&first_path).unwrap(), bytes);
    let final_path = dir.join("final.shapley");
    let final_run = command(&base, &input, Some("1")).arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&final_path).output().unwrap(); parsed(&final_run);
    assert_eq!(final_run.stdout, full.stdout);
    assert_eq!(fs::read(&final_path).unwrap(), fs::read(&full_path).unwrap());
    fs::write(&input, request().replace("\"wall_seconds\": 300", "\"wall_seconds\": 1e-9")).unwrap();
    let terminal = command(&base, &input, Some("3")).arg("--resume").arg(&final_path).output().unwrap();
    parsed(&terminal); assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn changed_laws_models_and_serial_formats_refuse_before_output_reservation() {
    let dir = scratch(); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request()).unwrap();
    let saved = dir.join("empty.shapley");
    let zero = command(&base, &input, Some("2")).arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(number(&progress(&zero), "samples_evaluated"), 0.0);
    let bytes = fs::read(&saved).unwrap(); let unused = dir.join("unused.shapley");
    for text in [request().replace("[[1, 0.6], [0.6, 1]]", "[[1, 1], [1, 1]]"),
        request().replace("\"73\"", "\"74\""), request().replace("\"samples\": 12", "\"samples\": 18")]
    {
        fs::write(&input, text).unwrap();
        let out = command(&base, &input, Some("2")).arg("--resume").arg(&saved)
            .arg("--checkpoint").arg(&unused).output().unwrap();
        assert!(!out.status.success()); assert!(out.stdout.is_empty()); assert!(!unused.exists());
    }
    fs::write(&input, request()).unwrap();
    let changed_base = dir.join("changed.json");
    fs::write(&changed_base, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    for (model, workers, destination) in [(&changed_base, Some("2"), &unused),
        (&base, None, &unused), (&base, Some("2"), &saved)]
    {
        let out = command(model, &input, workers).arg("--resume").arg(&saved)
            .arg("--checkpoint").arg(destination).output().unwrap();
        assert!(!out.status.success()); assert!(out.stdout.is_empty()); assert!(!unused.exists());
        assert_eq!(fs::read(&saved).unwrap(), bytes);
    }
    let missing = dir.join("missing.json");
    let out = command(&missing, &missing, Some("2")).args(["--qmc-replicates", "2"])
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::USAGE))); assert!(!unused.exists());
}

#[test]
fn timeouts_resume_but_physical_failures_never_become_selected_prefixes() {
    let dir = scratch(); let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json");
    fs::write(&input, request().replace("\"wall_seconds\": 300", "\"wall_seconds\": 1e-9")).unwrap();
    let saved = dir.join("timeout.shapley");
    let timeout = command(&base, &input, Some("3")).arg("--checkpoint").arg(&saved).output().unwrap();
    assert_eq!(progress(&timeout).str_field("termination"), Some("wall-time-budget"));
    fs::write(&input, request()).unwrap();
    let resumed = command(&base, &input, Some("2")).arg("--resume").arg(&saved).output().unwrap();
    assert_eq!(number(&parsed(&resumed), "samples_evaluated"), 12.0);
    let broken = dir.join("broken.json"); let original = fs::read_to_string(&base).unwrap();
    assert!(original.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, original.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.shapley");
    let out = command(&broken, &input, Some("3")).arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(out.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let next = dir.join("never.shapley");
    let out = command(&broken, &input, Some("1")).arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&next).output().unwrap();
    assert!(!out.status.success()); assert!(!next.exists());
}
