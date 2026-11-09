//! The authored affine strength family survives resume and scientific export.
use super::*;
use fs_viz::vtu::{DataAssociation, DataValues, VtuChecker};

const SOURCE: &str = include_str!("../../../../examples/marquee/bracket-3d-affine-strength.fsim");

fn numbers(value: &J) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn g1_g3_g5_affine_strength_cli_restores_its_joint_limits_and_exports_selected_envelopes() {
    let dir = scratch("affine-strength");
    // Explicit headroom keeps this small journey about retained physical
    // results. Native tests separately exercise active and violated caps.
    let source = SOURCE
        .replace(":compliance-limit-j 0.01", ":compliance-limit-j 0.1")
        .replace(":stress-limit-pa 40.0", ":stress-limit-pa 400.0");
    let source_path = dir.join("strength.fsim");
    let ledger = dir.join("strength.db");
    fs::write(&source_path, source).unwrap();
    let first = document(
        &command("study")
            .arg(&source_path)
            .arg(&ledger)
            .args(["--budget", "1"])
            .output()
            .unwrap(),
        6,
    );
    let first_report = data(&ledger, &first, "report_json");
    let original_fields = retained(&ledger, &first, "design");
    assert_eq!(number(&first_report, "iterations_completed"), 1.0);
    assert_eq!(
        first_report.path(&["gradient_check", "passed"]),
        Some(&J::Bool(true))
    );

    // A different process must use the original retained source, even when
    // the authored file is no longer present at its initial path.
    fs::rename(&source_path, dir.join("moved-strength.fsim")).unwrap();
    let resumed_output = command("study")
        .args(["--resume", first.str_field("run_id").unwrap()])
        .arg(&ledger)
        .args(["--budget", "1"])
        .output()
        .unwrap();
    assert!(
        matches!(resumed_output.status.code(), Some(0 | 6)),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&resumed_output.stdout),
        String::from_utf8_lossy(&resumed_output.stderr)
    );
    let result = document(&resumed_output, resumed_output.status.code().unwrap());
    let report = data(&ledger, &result, "report_json");
    assert_eq!(number(&report, "iterations_completed"), 2.0);
    assert_eq!(
        report.get("gradient_check"),
        first_report.get("gradient_check")
    );
    assert!(
        number(
            report.get("optimizer_work").unwrap(),
            "restoration_evaluations"
        ) >= 1.0
    );
    let before_work = first_report.get("work").unwrap();
    let after_work = report.get("work").unwrap();
    assert!(number(after_work, "linear_iterations") > number(before_work, "linear_iterations"));
    assert!(number(after_work, "geometry_boxes") > number(before_work, "geometry_boxes"));
    assert_eq!(retained(&ledger, &first, "design"), original_fields);

    let model = report.get("model").unwrap();
    assert_eq!(number(model, "compliance_limit_j"), 0.1);
    assert_eq!(
        number(model.get("strength_constraint").unwrap(), "stress_limit_pa"),
        400.0
    );
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
    let selected = report.get("selected").unwrap();
    assert!(number(selected, "joint_constraint") <= 0.0);
    assert!(number(selected, "worst_compliance_j") <= 0.1);
    let strength = selected.get("strength").unwrap();
    assert_eq!(strength.get("feasible"), Some(&J::Bool(true)));
    assert_eq!(number(strength, "stress_limit_pa"), 400.0);
    assert!(number(strength, "aggregate_pa") <= 400.0);
    assert!(number(strength, "constraint") <= 0.0);
    assert!(number(strength, "point_count") > 0.0);
    assert_eq!(
        number(strength, "point_fields"),
        3.0 * number(strength, "point_count")
    );
    for (name, envelope) in [
        ("worst_physical", "physical_envelope_pa"),
        ("worst_relaxed", "relaxed_envelope_pa"),
    ] {
        let witness = strength.get(name).unwrap();
        let coefficients = numbers(witness.get("load_coefficients").unwrap());
        assert_eq!(coefficients.len(), 2);
        assert!(coefficients.iter().all(|q| q.is_finite()));
        assert!(coefficients.iter().map(|q| q * q).sum::<f64>() <= 1.0 + 1e-10);
        assert_eq!(numbers(witness.get("position_m").unwrap()).len(), 3);
        assert_eq!(
            witness.str_field("stress_layout"),
            Some("xx, yy, zz, xy, yz, zx")
        );
        let tensor = numbers(witness.get("reference_stress_pa").unwrap());
        assert_eq!(tensor.len(), 6);
        let vm = (0.5
            * ((tensor[0] - tensor[1]).powi(2)
                + (tensor[1] - tensor[2]).powi(2)
                + (tensor[2] - tensor[0]).powi(2))
            + 3.0 * (tensor[3].powi(2) + tensor[4].powi(2) + tensor[5].powi(2)))
        .sqrt();
        let measured = number(witness, "reference_von_mises_pa");
        assert!((vm - measured).abs() <= 1e-9 * vm.max(measured).max(1e-20));
        assert!(number(witness, envelope) <= number(strength, "aggregate_pa"));
        assert!(number(witness, "relative_defect") < 1e-7);
        assert!(number(witness.get("oracle").unwrap(), "relative_duality_gap") < 1e-7);
    }

    let design = data(&ledger, &result, "design");
    assert_eq!(design.get("strength_constrained"), Some(&J::Bool(true)));
    let fields = design.get("selected").unwrap();
    let cells = fields.get("cells").and_then(J::as_array).unwrap();
    let envelopes = fields.get("strength").unwrap();
    let physical = numbers(envelopes.get("cell_physical_envelope_pa").unwrap());
    let relaxed = numbers(envelopes.get("cell_relaxed_envelope_pa").unwrap());
    let counts: Vec<i64> = envelopes
        .get("cell_sample_count")
        .and_then(J::as_array)
        .unwrap()
        .iter()
        .map(|v| v.number_raw().unwrap().parse().unwrap())
        .collect();
    assert_eq!(physical.len(), cells.len());
    assert_eq!(relaxed.len(), cells.len());
    assert_eq!(counts.len(), cells.len());
    assert!(counts.iter().all(|count| *count >= 0));
    assert_eq!(
        counts.iter().sum::<i64>() as f64,
        number(strength, "point_count")
    );
    assert_eq!(
        physical.iter().copied().fold(0.0_f64, f64::max),
        number(strength, "physical_envelope_max_pa")
    );
    assert_eq!(
        relaxed.iter().copied().fold(0.0_f64, f64::max),
        number(strength, "relaxed_envelope_max_pa")
    );

    let exported = document(
        &command("report")
            .arg(result.str_field("run_id").unwrap())
            .arg(&ledger)
            .output()
            .unwrap(),
        0,
    );
    for key in ["report_json", "design", "iterations"] {
        assert_eq!(
            fs::read(exported.str_field(key).unwrap()).unwrap(),
            retained(&ledger, &result, key)
        );
    }
    assert_eq!(exported.str_field("field_vtu_status"), Some("available"));
    assert_eq!(
        exported.get("field_vtu_source_design"),
        result.path(&["receipt", "design"])
    );
    let xml = fs::read_to_string(exported.str_field("field_vtu").unwrap()).unwrap();
    let grid = VtuChecker::parse_ascii(&xml).unwrap();
    for (name, expected) in [
        ("sampled_physical_stress_envelope", physical),
        ("sampled_relaxed_stress_envelope", relaxed),
    ] {
        let array = grid.arrays.iter().find(|a| a.name == name).unwrap();
        assert_eq!(array.association, DataAssociation::CellData);
        assert_eq!(array.unit.as_deref(), Some("Pa"));
        assert_eq!(array.components, 1);
        assert_eq!(array.values, DataValues::Float64(expected));
    }
    let samples = grid
        .arrays
        .iter()
        .find(|a| a.name == "stress_sample_count")
        .unwrap();
    assert_eq!(samples.association, DataAssociation::CellData);
    assert_eq!(samples.values, DataValues::Int64(counts));
}
