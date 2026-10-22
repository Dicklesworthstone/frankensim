//! Reuse the real-binary pilot fixtures; recovery never substitutes a model.
use super::*;

fn progress(output: &Output) -> J {
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)),
        "{}", String::from_utf8_lossy(&output.stderr));
    let report = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(report.str_field("schema"), Some("frankensim.cooling-network-uq.importance-pilot.progress.v1"));
    for field in ["estimated_failure_probability", "failure_probability_standard_error",
        "mean_k", "mean_likelihood_weight", "production_proposal_failures", "pilots"]
    {
        assert!(report.get(field).is_none(), "partial design published {field}");
    }
    report
}

#[test]
fn pilot_and_production_chunks_preserve_final_report_and_checkpoint_bits() {
    let dir = scratch("pilot-restart");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(12, 300.0, false)).unwrap();
    let full_path = dir.join("full.pilot");
    let full = command(&base, &input, LAYOUT).arg("--checkpoint").arg(&full_path).output().unwrap();
    parsed(&full);
    let first = dir.join("first.pilot");
    let first_report = progress(&command(&base, &input, LAYOUT).arg("--checkpoint").arg(&first)
        .args(["--max-new-samples", "2"]).output().unwrap());
    assert_eq!(first_report.str_field("stage"), Some("pilot"));
    assert_eq!(number(&first_report, "completed_pilots"), 0.0);
    let original = fs::read(&first).unwrap(); assert_eq!(&original[..8], b"FSIPL001");
    let second = dir.join("second.pilot");
    let second_report = progress(&command(&base, &input, LAYOUT).arg("--resume").arg(&first)
        .arg("--checkpoint").arg(&second).args(["--max-new-samples", "2"]).output().unwrap());
    assert_eq!(number(&second_report, "samples_evaluated"), 4.0);
    assert_eq!(number(&second_report, "completed_pilots"), 1.0);
    assert_eq!(number(&second_report, "current_evaluations"), 1.0);
    let third = dir.join("third.pilot");
    let third_report = progress(&command(&base, &input, LAYOUT).arg("--resume").arg(&second)
        .arg("--checkpoint").arg(&third).args(["--max-new-samples", "3"]).output().unwrap());
    assert_eq!(third_report.str_field("stage"), Some("production"));
    assert_eq!(number(&third_report, "completed_pilots"), 2.0);
    assert_eq!(number(&third_report, "samples_evaluated"), 7.0);
    assert_eq!(number(&third_report, "samples_evaluated_this_run"), 3.0);
    let completed = dir.join("completed.pilot");
    let result = command(&base, &input, LAYOUT).arg("--resume").arg(&third)
        .arg("--checkpoint").arg(&completed).output().unwrap();
    parsed(&result); assert_eq!(result.stdout, full.stdout);
    assert_eq!(fs::read(&completed).unwrap(), fs::read(&full_path).unwrap());
    assert_eq!(fs::read(&first).unwrap(), original);
    // A completed restore cannot need another physical solve.
    let short = dir.join("short.json"); fs::write(&short, request(12, 1e-9, false)).unwrap();
    let terminal = command(&base, &short, LAYOUT).arg("--resume").arg(&completed).output().unwrap();
    parsed(&terminal); assert_eq!(terminal.stdout, full.stdout);
}

#[test]
fn all_physics_paid_but_publication_pending_finishes_without_a_new_draw() {
    let dir = scratch("pilot-finalize");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(8, 300.0, false)).unwrap();
    let layout = "1,2,0.25,0.5,1,6";
    let saved = dir.join("full.pilot");
    let full = command(&base, &input, layout).arg("--checkpoint").arg(&saved).output().unwrap();
    parsed(&full);
    // Reproduce the reachable cancellation boundary between the last paid
    // producer call and outer publication, without altering any observations.
    let mut bytes = fs::read(&saved).unwrap(); bytes[40] = 1;
    let end = bytes.len() - 32;
    let sum = fs_blake3::hash_domain("org.frankensim.uq.pilot-importance.checkpoint.v1", &bytes[..end]);
    bytes[end..].copy_from_slice(sum.as_bytes());
    let pending = dir.join("pending.pilot"); fs::write(&pending, &bytes).unwrap();
    let copied = dir.join("zero.pilot");
    let zero = progress(&command(&base, &input, layout).arg("--resume").arg(&pending)
        .arg("--checkpoint").arg(&copied).args(["--max-new-samples", "0"]).output().unwrap());
    assert_eq!(number(&zero, "samples_evaluated_this_run"), 0.0);
    assert_eq!(number(&zero, "remaining_physical_calls"), 0.0);
    assert_eq!(zero.get("phase_work_pending"), Some(&J::Bool(true)));
    assert_eq!(fs::read(&copied).unwrap(), bytes);
    let finished = dir.join("finished.pilot");
    let output = command(&base, &input, layout).arg("--resume").arg(&copied)
        .arg("--checkpoint").arg(&finished).output().unwrap();
    parsed(&output); assert_eq!(output.stdout, full.stdout);
    assert_eq!(fs::read(&finished).unwrap(), fs::read(&saved).unwrap());
}

#[test]
fn changed_policy_model_corruption_and_existing_output_refuse_before_writing() {
    let dir = scratch("pilot-identity");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(12, 300.0, false)).unwrap();
    let saved = dir.join("empty.pilot");
    let zero = progress(&command(&base, &input, LAYOUT).arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "0"]).output().unwrap());
    assert_eq!(number(&zero, "samples_evaluated"), 0.0);
    let original = fs::read(&saved).unwrap();
    for (index, layout) in ["1,3,0.25,0.5,1,6", "2,2,0.25,0.5,1,6", "2,3,0.5,0.5,1,6",
        "2,3,0.25,0.25,1,6", "2,3,0.25,0.5,0.5,6", "2,3,0.25,0.5,1,8"].iter().enumerate()
    {
        let destination = dir.join(format!("changed-{index}.pilot"));
        let output = command(&base, &input, layout).arg("--resume").arg(&saved)
            .arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success()); assert!(output.stdout.is_empty()); assert!(!destination.exists());
    }
    let changed = dir.join("base.json"); fs::write(&changed, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    let corrupt = dir.join("corrupt.pilot");
    let mut bytes = original.clone(); let last = bytes.len() - 1; bytes[last] ^= 1;
    fs::write(&corrupt, bytes).unwrap();
    for (b, cp) in [(&changed, &saved), (&base, &corrupt), (&base, &saved)] {
        let destination = if cp == &saved && b == &base { saved.clone() } else { dir.join("rejected.pilot") };
        let output = command(b, &input, LAYOUT).arg("--resume").arg(cp)
            .arg("--checkpoint").arg(&destination).output().unwrap();
        assert!(!output.status.success()); assert!(output.stdout.is_empty());
        if destination != saved { assert!(!destination.exists()); }
    }
    assert_eq!(fs::read(&saved).unwrap(), original);
}

#[test]
fn timeout_retains_learning_state_but_physical_refusal_invalidates_it() {
    let dir = scratch("pilot-failure");
    let base = example("fan-correlated-hotspot.json");
    let input = dir.join("uq.json"); fs::write(&input, request(12, 1e-9, false)).unwrap();
    let saved = dir.join("timeout.pilot");
    let timed = progress(&command(&base, &input, LAYOUT).arg("--checkpoint").arg(&saved).output().unwrap());
    assert_eq!(timed.str_field("termination"), Some("wall-time-budget"));
    assert_eq!(number(&timed, "samples_evaluated"), 0.0);
    fs::write(&input, request(12, 300.0, false)).unwrap();
    let next = dir.join("continued.pilot");
    let resumed = progress(&command(&base, &input, LAYOUT).arg("--resume").arg(&saved)
        .arg("--checkpoint").arg(&next).args(["--max-new-samples", "2"]).output().unwrap());
    assert_eq!(number(&resumed, "samples_evaluated"), 2.0);
    let broken = dir.join("broken.json");
    fs::write(&broken, fs::read_to_string(&base).unwrap()
        .replace("\"linear_iterations\": 20000", "\"linear_iterations\": 0")).unwrap();
    let failed = dir.join("failed.pilot");
    let output = command(&broken, &input, LAYOUT).arg("--checkpoint").arg(&failed).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::REFUSED)));
    assert!(output.stdout.is_empty());
    assert!(fs::read(&failed).unwrap().starts_with(b"FRANKENSIM-UQ-FAILED\n"));
    let rejected = dir.join("must-not-resume.pilot");
    let output = command(&broken, &input, LAYOUT).arg("--resume").arg(&failed)
        .arg("--checkpoint").arg(&rejected).output().unwrap();
    assert!(!output.status.success()); assert!(!rejected.exists());
}

#[test]
fn complete_transient_observations_survive_a_production_restart() {
    let dir = scratch("pilot-transient-restart");
    let base = example("nonlinear-contact-pulse.json");
    let input = dir.join("uq.json"); fs::write(&input, request(8, 300.0, true)).unwrap();
    let layout = "1,2,0.25,0.5,1,6";
    let full = command(&base, &input, layout).output().unwrap(); parsed(&full);
    let saved = dir.join("partial.pilot");
    let first = progress(&command(&base, &input, layout).arg("--checkpoint").arg(&saved)
        .args(["--max-new-samples", "3"]).output().unwrap());
    assert_eq!(number(&first, "samples_evaluated"), 3.0);
    assert_eq!(first.str_field("stage"), Some("production"));
    assert_eq!(first.get("qoi").unwrap().str_field("observation"), Some("one-completed-trajectory"));
    let resumed = command(&base, &input, layout).arg("--resume").arg(&saved).output().unwrap();
    parsed(&resumed); assert_eq!(resumed.stdout, full.stdout);
}
