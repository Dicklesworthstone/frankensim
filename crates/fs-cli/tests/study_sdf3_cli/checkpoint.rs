use super::*;

fn progress(output: &Output) -> Vec<J> {
    std::str::from_utf8(&output.stderr).unwrap().lines()
        .filter_map(|line| J::parse(line).ok())
        .filter(|doc| doc.str_field("schema") == Some("frankensim.cli.sdf3-progress.v1"))
        .collect()
}

#[test]
fn g1_g5_sdf3_public_stage_budget_resume_and_retained_exports() {
    let dir = scratch("stage-resume");
    let source = FIXTURE.replace(":updates-per-stage 3", ":updates-per-stage 1");
    let (reference_db, reference) = run(&dir, "reference", &source, 0);
    let input = dir.join("split.fsim");
    let ledger = dir.join("split.db");
    fs::write(&input, &source).unwrap();
    let first_output = command("study").arg(&input).arg(&ledger)
        .args(["--budget", "1"]).output().unwrap();
    let first = document(&first_output, 6);
    let first_progress = progress(&first_output);
    assert_eq!(first_progress.len(), 1);
    assert_eq!(first_progress[0].str_field("run_id"), first.str_field("run_id"));
    let partial = document(&command("report").arg(first.str_field("run_id").unwrap())
        .arg(&ledger).output().unwrap(), 0);
    assert_eq!(partial.str_field("study_status"), Some("budget-exhausted"));
    let initial_state = retained(&ledger, &first, "checkpoint");
    let mut current = first.clone();
    for stage in [2, 3] {
        let previous = current.str_field("run_id").unwrap().to_string();
        let output = command("study").args(["--resume", &previous])
            .arg(&ledger).args(["--budget", "1"]).output().unwrap();
        current = document(&output, if stage == 3 { 0 } else { 6 });
        let emitted = progress(&output);
        assert_eq!(emitted.len(), 1, "replayed stages must not be republished");
        assert_eq!(emitted[0].str_field("run_id"), current.str_field("run_id"));
        let receipt = current.get("receipt").unwrap();
        assert_eq!(number(receipt, "stages_completed"), stage as f64);
        assert_eq!(number(receipt, "stages_replayed"), (stage - 1) as f64);
        assert_eq!(number(receipt, "stage_limit_this_invocation"), 1.0);
        assert_eq!(receipt.str_field("predecessor"), previous.strip_prefix("study-"));
        assert_eq!(receipt.str_field("resume_mode"), Some("verified-stage-replay-v1"));
    }
    for key in ["design", "iterations", "checkpoint"] {
        assert_eq!(retained(&reference_db, &reference, key), retained(&ledger, &current, key));
    }
    assert_eq!(initial_state, retained(&ledger, &first, "checkpoint"));
    for key in ["linear_iterations", "linear_solves", "geometry_boxes", "geometry_points"] {
        assert!(number(current.path(&["receipt", "work"]).unwrap(), key)
            > number(reference.path(&["receipt", "work"]).unwrap(), key));
    }
    let id = current.str_field("run_id").unwrap();
    let no_work = command("study").args(["--resume", id]).arg(&ledger).output().unwrap();
    assert_eq!(document(&no_work, 0), current);
    assert!(progress(&no_work).is_empty());
    let exported = document(&command("report").arg(id).arg(&ledger).output().unwrap(), 0);
    for key in ["report_html", "report_json", "design", "iterations"] {
        assert_eq!(fs::read(exported.str_field(key).unwrap()).unwrap(), retained(&ledger, &current, key));
    }
    document(&command("package").arg(id).arg(&ledger).output().unwrap(), 0);
}

#[test]
fn g4_sdf3_progress_receipts_survive_later_stage_failure() {
    let dir = scratch("durable-prefix");
    let input = dir.join("limited.fsim");
    let ledger = dir.join("limited.db");
    fs::write(&input, FIXTURE.replace(":updates-per-stage 3", ":updates-per-stage 1")
        .replace(":maximum-leaves 2048", ":maximum-leaves 8")).unwrap();
    let output = command("study").arg(input).arg(&ledger).output().unwrap();
    let terminal = document(&output, 6);
    let emitted = progress(&output);
    assert_eq!(emitted.len(), 1);
    let durable_id = emitted[0].str_field("run_id").unwrap();
    assert_ne!(Some(durable_id), terminal.str_field("run_id"));
    let report = document(&command("report").arg(durable_id).arg(&ledger).output().unwrap(), 0);
    assert_eq!(report.str_field("study_status"), Some("checkpointed"));
    let summary = J::parse(&fs::read_to_string(report.str_field("report_json").unwrap()).unwrap()).unwrap();
    assert_eq!(number(&summary, "stages_completed"), 1.0);
    assert_eq!(number(&summary, "active_cells"), 8.0);
    let retry = document(&command("study").args(["--resume", durable_id])
        .arg(&ledger).output().unwrap(), 6);
    assert_eq!(retained(&ledger, &terminal, "checkpoint"), retained(&ledger, &retry, "checkpoint"));
    assert!(retry.path(&["receipt", "stop", "refinement"]).and_then(J::as_str).unwrap().contains("LeafBudget"));
}
