//! Native stress-constrained study, durable incumbent, sealed exports and caps.
use super::*;

const STRESS: &str = include_str!("../../../../examples/marquee/bracket-3d-stress.fsim");

#[test]
fn g1_g5_stress_cli_reduces_volume_keeps_independent_fields_and_exports_sealed_results() {
    let dir = scratch("stress");
    let source = STRESS.replace(":max-updates 80", ":max-updates 4");
    let (ledger, result) = run(&dir, "stress", &source, 6);
    assert_eq!(result.str_field("status"), Some("budget-exhausted"));
    let report = data(&ledger, &result, "report_json");
    assert_eq!(report.str_field("objective_unit"), Some("1"));
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
    assert_eq!(
        report.path(&["gradient_check", "passed"]),
        Some(&J::Bool(true))
    );
    assert_eq!(
        report.path(&["stress_measure", "maximum_stress_is_constrained"]),
        Some(&J::Bool(false))
    );
    assert_eq!(
        report.path(&["stop", "reason"]).and_then(J::as_str),
        Some("IterationLimit")
    );
    let iterations = data(&ledger, &result, "iterations");
    let rows = iterations.get("history").and_then(J::as_array).unwrap();
    assert_eq!(rows.len(), 5);
    let selected = report.get("selected").unwrap();
    assert!(number(selected, "volume_fraction") < number(&rows[0], "volume_fraction"));
    assert!(number(selected, "stress_aggregate_pa") <= 8.0 * (1.0 + 2e-6));
    let design = data(&ledger, &result, "design");
    let field = design.get("selected").unwrap();
    let cells = field.get("cells").and_then(J::as_array).unwrap();
    let volume: f64 = cells.iter().map(|c| number(c, "cut_volume_m3")).sum();
    let material: f64 = cells
        .iter()
        .map(|c| number(c, "cut_volume_m3") * number(c, "projected_density"))
        .sum();
    assert!((material / volume - number(selected, "volume_fraction")).abs() < 1e-12);
    let fields = field.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);
    let again = run(&dir, "repeated", &source, 6);
    for key in ["design", "iterations"] {
        assert_eq!(
            retained(&ledger, &result, key),
            retained(&again.0, &again.1, key)
        );
    }
    let id = result.str_field("run_id").unwrap();
    let exported = document(&command("report").arg(id).arg(&ledger).output().unwrap(), 0);
    for key in ["report_html", "report_json", "design", "iterations"] {
        assert_eq!(
            fs::read(exported.str_field(key).unwrap()).unwrap(),
            retained(&ledger, &result, key)
        );
    }
    let packaged = document(
        &command("package").arg(id).arg(&ledger).output().unwrap(),
        0,
    );
    let package_bytes = fs::read(packaged.str_field("package").unwrap()).unwrap();
    assert_eq!(package_bytes, retained(&ledger, &result, "package"));
    let package =
        fs_package::EvidencePackage::from_json(std::str::from_utf8(&package_bytes).unwrap())
            .unwrap();
    assert!(fs_checker::check(&package).passed());
    let claims = package.declared_claims_unverified();
    assert_eq!(claims.len(), 4);
    for claim in claims {
        assert!(matches!(
            claim.declared_color_unverified(),
            fs_evidence::Color::Estimated { .. }
        ));
        let artifact = J::parse(claim.statement()).unwrap();
        let kind = artifact.str_field("kind").unwrap();
        let contents = artifact.str_field("contents").unwrap().as_bytes();
        assert_eq!(contents, retained(&ledger, &result, kind));
        assert_eq!(
            artifact.str_field("hash"),
            Some(fs_blake3::hash_bytes(contents).to_hex().as_str())
        );
    }
    let finished = command("study")
        .arg("--resume")
        .arg(id)
        .arg(&ledger)
        .output()
        .unwrap();
    let finished = document(&finished, 6);
    assert_eq!(finished.str_field("run_id"), Some(id), "an exhausted update target is unchanged");
}

#[test]
fn g5_stress_cli_resumes_in_a_new_process_with_identical_design_and_spent_original_budgets() {
    let dir = scratch("stress-resume");
    let source = STRESS.replace(":max-updates 80", ":max-updates 4");
    let (reference_ledger, reference) = run(&dir, "uninterrupted", &source, 6);
    let path = dir.join("segmented.fsim");
    let ledger = dir.join("segmented.db");
    fs::write(&path, &source).unwrap();
    let first = document(&command("study").arg(&path).arg(&ledger)
        .args(["--budget", "2"]).output().unwrap(), 6);
    let first_id = first.str_field("run_id").unwrap();
    let before = data(&ledger, &first, "checkpoint");
    assert_eq!(first.path(&["receipt", "resume_supported"]), Some(&J::Bool(true)));
    // Resume uses the immutable retained source, even when the working file is
    // later edited into a different design problem.
    fs::write(&path, source.replace(":stress-limit-pa 8.0", ":stress-limit-pa 0.0001")).unwrap();
    let resumed = document(&command("study").arg("--resume").arg(first_id).arg(&ledger)
        .args(["--budget", "2"]).output().unwrap(), 6);
    assert_eq!(resumed.path(&["receipt", "iterations_completed"]).and_then(J::as_f64), Some(4.0));
    for artifact in ["design", "iterations"] {
        assert_eq!(retained(&ledger, &resumed, artifact), retained(&reference_ledger, &reference, artifact));
    }
    assert_eq!(data(&ledger, &first, "checkpoint"), before, "original checkpoint is immutable");
    let report = data(&ledger, &resumed, "report_json");
    let uninterrupted = data(&reference_ledger, &reference, "report_json");
    let costs = report.get("optimizer_work").unwrap();
    let extra = number(costs, "restoration_evaluations");
    assert!((1.0..=2.0).contains(&extra), "only accepted/distinct incumbent endpoints re-solve");
    assert_eq!(number(costs, "evaluations"),
        number(uninterrupted.get("optimizer_work").unwrap(), "evaluations") + extra);
    assert_eq!(report.get("gradient_check"), uninterrupted.get("gradient_check"));
    assert_eq!(report.get("stop"), uninterrupted.get("stop"));
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
    assert!(number(report.get("work").unwrap(), "linear_iterations")
        > number(uninterrupted.get("work").unwrap(), "linear_iterations"));
    assert!(number(costs, "total_evaluations_including_gradient_gate") <= 2000.0);
    let complete_id = resumed.str_field("run_id").unwrap();
    let again = document(&command("study").arg("--resume").arg(complete_id).arg(&ledger)
        .output().unwrap(), 6);
    assert_eq!(again.str_field("run_id"), Some(complete_id));
}

#[test]
fn g4_stress_resume_cannot_renew_an_exhausted_evaluation_allowance() {
    let dir = scratch("stress-resume-budget");
    let (ledger, first) = run(&dir, "limited",
        &STRESS.replace(":max-evaluations 2000", ":max-evaluations 6"), 6);
    assert_eq!(first.path(&["receipt", "iterations_completed"]).and_then(J::as_f64), Some(0.0));
    let old = retained(&ledger, &first, "checkpoint");
    let stopped = command("study").arg("--resume").arg(first.str_field("run_id").unwrap())
        .arg(&ledger).args(["--budget", "80"]).output().unwrap();
    document(&stopped, 6);
    assert!(String::from_utf8_lossy(&stopped.stdout).contains("cli-study-sdf3-resume-budget"));
    assert_eq!(retained(&ledger, &first, "checkpoint"), old);
}

#[test]
fn g4_stress_cli_keeps_a_durable_accepted_prefix_and_never_labels_an_infeasible_stop_completed() {
    let dir = scratch("stress-prefix");
    let source = dir.join("study.fsim");
    let ledger = dir.join("study.db");
    fs::write(&source, STRESS).unwrap();
    let output = command("study")
        .arg(&source)
        .arg(&ledger)
        .args(["--budget", "2"])
        .output()
        .unwrap();
    let result = document(&output, 6);
    assert_eq!(
        result
            .path(&["receipt", "iterations_completed"])
            .and_then(J::as_f64),
        Some(2.0)
    );
    let checkpoints: Vec<_> = std::str::from_utf8(&output.stderr)
        .unwrap()
        .lines()
        .filter_map(|line| J::parse(line).ok())
        .filter(|doc| doc.str_field("schema") == Some("frankensim.cli.sdf3-stress-progress.v1"))
        .collect();
    assert_eq!(
        checkpoints.len(),
        2,
        "initial and first accepted update are durable before later work"
    );
    for checkpoint in checkpoints {
        document(
            &command("report")
                .arg(checkpoint.str_field("run_id").unwrap())
                .arg(&ledger)
                .output()
                .unwrap(),
            0,
        );
    }
    let (bad_ledger, bad) = run(
        &dir,
        "infeasible",
        &STRESS
            .replace(":max-updates 80", ":max-updates 1")
            .replace(":stress-limit-pa 8.0", ":stress-limit-pa 0.0001"),
        6,
    );
    assert_ne!(bad.str_field("status"), Some("completed"));
    let report = data(&bad_ledger, &bad, "report_json");
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(false)));
    assert_eq!(report.str_field("selected_design"), Some("last-accepted"));
    assert!(number(report.get("selected").unwrap(), "constraint_violation") > 0.0);
}

const STRESS_REGIONS: &str =
    include_str!("../../../../examples/marquee/bracket-3d-stress-regions.fsim");

#[test]
fn g5_stress_regions_survive_a_new_process_and_remain_in_sealed_design_exports() {
    let dir = scratch("stress-regions-resume");
    let (reference_ledger, reference) = run(&dir, "uninterrupted", STRESS_REGIONS, 6);
    let path = dir.join("segmented.fsim");
    let ledger = dir.join("segmented.db");
    fs::write(&path, STRESS_REGIONS).unwrap();
    let first = document(&command("study").arg(&path).arg(&ledger)
        .args(["--budget", "2"]).output().unwrap(), 6);
    let first_id = first.str_field("run_id").unwrap();
    let original_checkpoint = retained(&ledger, &first, "checkpoint");
    // Resume must recover the authored region map from the retained source,
    // even when the working file now places the void in a different cell.
    fs::write(&path, STRESS_REGIONS.replace(
        ":void (((0.0 0.5 0.5) (0.5 1.0 1.0)))",
        ":void (((0.5 0.5 0.5) (1.0 1.0 1.0)))",
    )).unwrap();
    let resumed = document(&command("study").arg("--resume").arg(first_id).arg(&ledger)
        .args(["--budget", "2"]).output().unwrap(), 6);
    assert_eq!(resumed.path(&["receipt", "iterations_completed"]).and_then(J::as_f64), Some(4.0));
    assert_eq!(retained(&ledger, &first, "checkpoint"), original_checkpoint);
    for artifact in ["design", "iterations"] {
        assert_eq!(retained(&ledger, &resumed, artifact),
            retained(&reference_ledger, &reference, artifact));
    }
    for result in [&first, &resumed] {
        let report = data(&ledger, result, "report_json");
        assert_eq!(report.path(&["gradient_check", "passed"]), Some(&J::Bool(true)));
        assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
        let design = data(&ledger, result, "design");
        for key in ["selected", "last_accepted"] {
            let field = design.get(key).unwrap();
            let cells = field.get("cells").and_then(J::as_array).unwrap();
            let mut counts = [0; 3];
            for cell in cells {
                let index = cell.get("index").and_then(J::as_array).unwrap();
                let index: Vec<_> = index.iter().map(|v| v.as_f64().unwrap()).collect();
                let expected = if index == [0.0, 0.0, 0.0] { "solid" }
                    else if index == [0.0, 1.0, 1.0] { "void" } else { "design" };
                assert_eq!(cell.str_field("physical_region"), Some(expected));
                match expected {
                    "solid" => { counts[0] += 1; assert_eq!(number(cell, "projected_density"), 1.0); }
                    "void" => { counts[1] += 1; assert_eq!(number(cell, "projected_density"), 0.0); }
                    _ => { counts[2] += 1; assert!(number(cell, "projected_density") > 0.0); }
                }
            }
            assert_eq!(counts, [1, 1, 6]);
        }
    }
    let id = resumed.str_field("run_id").unwrap();
    let exported = document(&command("report").arg(id).arg(&ledger).output().unwrap(), 0);
    assert_eq!(fs::read(exported.str_field("design").unwrap()).unwrap(),
        retained(&ledger, &resumed, "design"));
}
