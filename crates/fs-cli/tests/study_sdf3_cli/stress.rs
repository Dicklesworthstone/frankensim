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
    let refused = command("study")
        .arg("--resume")
        .arg(id)
        .arg(&ledger)
        .output()
        .unwrap();
    document(&refused, 4);
    assert!(String::from_utf8_lossy(&refused.stdout).contains("cli-study-sdf3-resume-unsupported"));
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
