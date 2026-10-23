//! Original physical QMC results must be independent of the cooling worker count.
use super::*;

fn qmc(base: &Path, request: &Path, workers: Option<&str>, replicates: &str) -> Command {
    let mut c = command(base, request, workers);
    c.args(["--qmc-replicates", replicates]); c
}
fn qmc_progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)), "{}", String::from_utf8_lossy(&output.stderr));
    let report = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(report.str_field("schema"), Some("frankensim.cooling-network-uq.qmc.parallel.progress.v1"));
    for key in ["mean_k", "replicate_means_k", "between_replicate_standard_error_k",
        "estimated_probability_of_compliance", "probability_confidence_sequence"]
    { assert!(report.get(key).is_none(), "partial design published {key}"); }
    report
}

#[test]
fn parallel_qmc_steady_and_complete_transient_outputs_match_serial_qmc() {
    for (base, request) in [("fan-correlated-hotspot.json", "uq-fan-hotspot.json"),
        ("nonlinear-contact-pulse.json", "uq-nonlinear-contact-pulse.json")]
    {
        let base = example(base); let request = example(request);
        let serial = qmc(&base, &request, None, "2").output().unwrap();
        assert_eq!(parsed(&serial).str_field("schema"), Some("frankensim.cooling-network-uq.qmc.v1"));
        for workers in ["1", "2", "3"] {
            let out = qmc(&base, &request, Some(workers), "2").output().unwrap();
            parsed(&out); assert_eq!(out.stdout, serial.stdout);
        }
    }
}

#[test]
fn interrupted_nets_resume_with_different_workers_and_identical_final_bytes() {
    let dir = scratch("qmc-resume");
    let base = example("fan-correlated-hotspot.json"); let request = example("uq-fan-hotspot.json");
    let full_path = dir.join("full.qbatch");
    let full = qmc(&base, &request, Some("4"), "2").arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first_path = dir.join("first.qbatch");
    let out = qmc(&base, &request, Some("3"), "2").arg("--checkpoint").arg(&first_path)
        .args(["--max-new-samples", "3"]).output().unwrap();
    let p = qmc_progress(&out);
    assert_eq!(number(&p, "samples_evaluated"), 3.0);
    assert_eq!(number(&p, "next_replicate_ordinal"), 0.0);
    assert_eq!(number(&p, "next_point_ordinal"), 3.0);
    let source = fs::read(&first_path).unwrap(); assert_eq!(&source[..8], b"FSQMB001");
    let second_path = dir.join("second.qbatch");
    let out = qmc(&base, &request, Some("2"), "2").arg("--resume").arg(&first_path)
        .arg("--checkpoint").arg(&second_path).args(["--max-new-samples", "2"]).output().unwrap();
    let p = qmc_progress(&out);
    assert_eq!(number(&p, "samples_evaluated"), 5.0);
    assert_eq!(number(&p, "samples_evaluated_this_run"), 2.0);
    assert_eq!(number(&p, "next_replicate_ordinal"), 1.0);
    assert_eq!(fs::read(&first_path).unwrap(), source);
    let last_path = dir.join("last.qbatch");
    let out = qmc(&base, &request, Some("1"), "2").arg("--resume").arg(&second_path)
        .arg("--checkpoint").arg(&last_path).output().unwrap();
    parsed(&out); assert_eq!(out.stdout, full.stdout);
    assert_eq!(fs::read(&last_path).unwrap(), fs::read(&full_path).unwrap());
    let short = dir.join("short.json");
    let text = fs::read_to_string(&request).unwrap(); assert!(text.contains("\"wall_seconds\": 300"));
    fs::write(&short, text.replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let out = qmc(&base, &short, Some("2"), "2").arg("--resume").arg(&last_path).output().unwrap();
    parsed(&out); assert_eq!(out.stdout, full.stdout);
}

#[test]
fn qmc_layout_and_checkpoint_kind_stay_bound_and_confidence_options_conflict() {
    let dir = scratch("qmc-admission");
    let base = example("fan-correlated-hotspot.json"); let request = example("uq-fan-hotspot.json");
    let source = dir.join("empty.qbatch");
    let out = qmc(&base, &request, Some("2"), "2").arg("--checkpoint").arg(&source)
        .args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(number(&qmc_progress(&out), "samples_evaluated"), 0.0);
    let bytes = fs::read(&source).unwrap();
    let unused = dir.join("unused.qbatch");
    for (workers, replicates) in [(Some("2"), "4"), (None, "2")] {
        let out = qmc(&base, &request, workers, replicates).arg("--resume").arg(&source)
            .arg("--checkpoint").arg(&unused).output().unwrap();
        assert!(!out.status.success()); assert!(out.stdout.is_empty()); assert!(!unused.exists());
    }
    let out = qmc(&base, &request, Some("2"), "2").arg("--resume").arg(&source)
        .arg("--checkpoint").arg(&source).output().unwrap();
    assert!(!out.status.success()); assert_eq!(fs::read(&source).unwrap(), bytes);
    let out = qmc(&base, &request, Some("2"), "3").arg("--checkpoint").arg(&unused).output().unwrap();
    assert!(!out.status.success()); assert!(!unused.exists());
    let missing = dir.join("missing.json");
    let out = qmc(&missing, &missing, Some("2"), "2").args([
        "--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2",
    ]).arg("--checkpoint").arg(&unused).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::USAGE))); assert!(!unused.exists());
}

#[test]
fn qmc_deadlines_preserve_work_and_real_failures_never_become_missing_points() {
    let dir = scratch("qmc-failures");
    let base = example("fan-correlated-hotspot.json"); let request = example("uq-fan-hotspot.json");
    let short = dir.join("short.json");
    fs::write(&short, fs::read_to_string(&request).unwrap().replace("\"wall_seconds\": 300", "\"wall_seconds\": 0.000000001")).unwrap();
    let saved = dir.join("timeout.qbatch");
    let out = qmc(&base, &short, Some("3"), "2").arg("--checkpoint").arg(&saved).output().unwrap();
    assert_eq!(qmc_progress(&out).str_field("termination"), Some("wall-time-budget"));
    let continued = qmc(&base, &request, Some("1"), "2").arg("--resume").arg(&saved).output().unwrap();
    let serial = qmc(&base, &request, None, "2").output().unwrap();
    parsed(&continued); parsed(&serial); assert_eq!(continued.stdout, serial.stdout);
    let broken = dir.join("broken.json");
    let text = fs::read_to_string(&base).unwrap(); assert!(text.contains("\"linear_iterations\": 20000"));
    fs::write(&broken, text.replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.qbatch");
    let out = qmc(&broken, &request, Some("2"), "2").arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(out.status.code(), Some(i32::from(fs_cli::exit::REFUSED))); assert!(out.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let unused = dir.join("unused.qbatch");
    let out = qmc(&broken, &request, Some("1"), "2").arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&unused).output().unwrap();
    assert!(!out.status.success()); assert!(!unused.exists());
}
