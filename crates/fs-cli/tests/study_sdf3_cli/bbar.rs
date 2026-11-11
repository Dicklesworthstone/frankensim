use super::*;

#[test]
fn bbar_cli_solves_nearly_incompressible_material_updates_density_and_exports_report() {
    let dir = scratch("bbar");
    let legacy = FIXTURE
        .replace(":updates-per-stage 3", ":updates-per-stage 1")
        .replace(
            ":youngs-modulus-pa 1048576.0",
            ":youngs-modulus-pa 1000000.0",
        )
        .replace(":poissons-ratio 0.3", ":poissons-ratio 0.49")
        .replace(
            ":schedule ((1.0 1.0) (2.0 2.0) (3.0 8.0))",
            ":schedule ((1.0 1.0))",
        );

    // The old declaration must still refuse this material before creating a
    // ledger. Near-incompressible admission requires the explicit formulation.
    let legacy_path = dir.join("legacy.fsim");
    let legacy_ledger = dir.join("legacy.db");
    fs::write(&legacy_path, &legacy).unwrap();
    let refused = command("study")
        .arg(&legacy_path)
        .arg(&legacy_ledger)
        .output()
        .unwrap();
    assert_eq!(refused.status.code(), Some(4));
    let diagnostic = String::from_utf8_lossy(&refused.stderr);
    assert!(diagnostic.contains("cli-study-sdf3-input"), "{diagnostic}");
    assert!(diagnostic.contains("Poisson ratio"), "{diagnostic}");
    assert!(!legacy_ledger.exists());

    let source = legacy.replace(
        ":poissons-ratio 0.49",
        ":poissons-ratio 0.49 :formulation bbar",
    );
    let (ledger, result) = run(&dir, "bbar", &source, 0);
    assert_eq!(result.str_field("status"), Some("completed"));
    let iterations = data(&ledger, &result, "iterations");
    let stages = iterations.get("stages").and_then(J::as_array).unwrap();
    assert_eq!(stages.len(), 1);
    assert_eq!(stages[0].get("gradient_check_passed"), Some(&J::Bool(true)));
    let history = stages[0].get("history").and_then(J::as_array).unwrap();
    assert!(history.len() >= 2, "B-bar must perform a real design update");
    let baseline = number(&history[0], "compliance_j");
    let endpoint = history.last().unwrap();
    assert!(baseline > 0.0);
    assert!(number(endpoint, "compliance_j") < baseline);
    assert!(number(endpoint, "volume_fraction") <= 0.50000001);
    assert!(result.path(&["receipt", "work", "linear_iterations"])
        .and_then(J::as_f64).unwrap() > 0.0);

    let design = data(&ledger, &result, "design");
    assert_eq!(design.str_field("state"), Some("solved"));
    assert!(!design.get("cells").and_then(J::as_array).unwrap().is_empty());
    let fields = design.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);

    let exported = document(&command("report")
        .arg(result.str_field("run_id").unwrap())
        .arg(&ledger).output().unwrap(), 0);
    for key in ["report_html", "report_json", "design", "iterations"] {
        assert_eq!(
            fs::read(exported.str_field(key).unwrap()).unwrap(),
            retained(&ledger, &result, key),
        );
    }
}
