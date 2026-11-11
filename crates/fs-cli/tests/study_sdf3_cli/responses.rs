//! Signed displacement, force and moment witnesses survive resume and export.
use super::*;
use fs_viz::vtu::{DataAssociation, DataValues, VtuChecker};

const SOURCE: &str = include_str!("../../../../examples/marquee/bracket-3d-affine-responses.fsim");

fn numbers(value: &J) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .collect()
}

#[test]
fn g1_g3_g5_affine_response_cli_preserves_all_signed_bounds_and_exports_each_adverse_field() {
    let dir = scratch("affine-responses");
    let source_path = dir.join("responses.fsim");
    let ledger = dir.join("responses.db");
    fs::write(&source_path, SOURCE).unwrap();
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
    assert_eq!(number(&first_report, "iterations_completed"), 1.0);
    let first_pointer = first.str_field("run_id").unwrap();
    let first_hash = ContentHash::from_hex(first_pointer.strip_prefix("study-").unwrap()).unwrap();
    let original_receipt = Ledger::open(ledger.to_str().unwrap())
        .unwrap()
        .get_artifact_bounded(&first_hash, 16 * 1024 * 1024)
        .unwrap()
        .unwrap();

    // Both continuation and visualization must use the retained study after
    // its source file has moved, with every original response bound still live.
    fs::rename(&source_path, dir.join("moved-responses.fsim")).unwrap();
    let continued_output = command("study")
        .args(["--resume", first_pointer])
        .arg(&ledger)
        .args(["--budget", "1"])
        .output()
        .unwrap();
    assert!(
        matches!(continued_output.status.code(), Some(0 | 6)),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&continued_output.stdout),
        String::from_utf8_lossy(&continued_output.stderr)
    );
    let result = document(&continued_output, continued_output.status.code().unwrap());
    let report = data(&ledger, &result, "report_json");
    assert_eq!(number(&report, "iterations_completed"), 2.0);
    assert_eq!(report.get("model"), first_report.get("model"));
    assert_eq!(
        report.get("gradient_check"),
        first_report.get("gradient_check")
    );
    assert_eq!(
        result.path(&["receipt", "source"]),
        first.path(&["receipt", "source"])
    );
    assert!(
        number(
            report.get("optimizer_work").unwrap(),
            "restoration_evaluations"
        ) >= 1.0
    );
    let retained_receipt = Ledger::open(ledger.to_str().unwrap())
        .unwrap()
        .get_artifact_bounded(&first_hash, 16 * 1024 * 1024)
        .unwrap()
        .unwrap();
    assert_eq!(retained_receipt, original_receipt);
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
    let selected = report.get("selected").unwrap();
    assert!(number(selected, "joint_constraint") <= 0.0);
    assert_eq!(
        selected.path(&["strength", "feasible"]),
        Some(&J::Bool(true))
    );
    let response = selected.get("responses").unwrap();
    assert_eq!(response.get("feasible"), Some(&J::Bool(true)));
    assert!(number(response, "constraint") <= 0.0);
    assert!(number(response, "worst_residual") <= 0.0);
    let limits = response.get("limits").and_then(J::as_array).unwrap();
    let expected = [
        (
            "tip-z-lower",
            "mean-displacement",
            "m",
            "at-least",
            -0.02 / 1024.0,
        ),
        (
            "tip-z-upper",
            "mean-displacement",
            "m",
            "at-most",
            0.02 / 1024.0,
        ),
        ("support-z-lower", "support-force", "N", "at-least", -2.0),
        ("support-z-upper", "support-force", "N", "at-most", 2.0),
        ("support-y-upper", "support-moment", "N*m", "at-most", 1.0),
    ];
    assert_eq!(limits.len(), expected.len());
    for (index, (limit, (name, kind, unit, sense, bound))) in
        limits.iter().zip(expected).enumerate()
    {
        assert_eq!(number(limit, "index"), index as f64);
        assert_eq!(limit.str_field("name"), Some(name));
        assert_eq!(limit.str_field("kind"), Some(kind));
        assert_eq!(limit.str_field("unit"), Some(unit));
        assert_eq!(limit.str_field("sense"), Some(sense));
        assert_eq!(number(limit, "bound"), bound);
        assert!(number(limit, "scale") > 0.0);
        assert!(number(limit, "residual") <= 0.0);
        assert!(number(limit, "witness_residual") <= 0.0);
        let coefficients = numbers(limit.get("coefficients").unwrap());
        assert_eq!(coefficients.len(), 2);
        assert!(coefficients.iter().all(|q| q.is_finite()));
        assert!(coefficients.iter().map(|q| q * q).sum::<f64>() <= 1.0 + 1e-10);
        for key in ["worst_value", "witness_value"] {
            let value = number(limit, key);
            assert!(value.is_finite());
            assert!(if sense == "at-most" {
                value <= bound
            } else {
                value >= bound
            });
        }
    }
    for (lower, upper) in [(0, 1), (2, 3)] {
        let lower = numbers(limits[lower].get("coefficients").unwrap());
        let upper = numbers(limits[upper].get("coefficients").unwrap());
        assert_eq!(lower, upper.iter().map(|q| -q).collect::<Vec<_>>());
    }

    let design = data(&ledger, &result, "design");
    assert_eq!(design.get("response_constrained"), Some(&J::Bool(true)));
    assert_eq!(design.get("strength_constrained"), Some(&J::Bool(true)));
    let fields = design.get("selected").unwrap();
    let offset = numbers(fields.get("offset_displacement_m").unwrap());
    let basis: Vec<_> = fields
        .get("basis_displacements_m")
        .and_then(J::as_array)
        .unwrap()
        .iter()
        .map(numbers)
        .collect();
    let witnesses = fields
        .path(&["responses", "witnesses"])
        .and_then(J::as_array)
        .unwrap();
    assert_eq!(witnesses.len(), limits.len());
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
    assert_eq!(
        grid.arrays
            .iter()
            .filter(|a| a.name.starts_with("response_displacement_"))
            .count(),
        limits.len()
    );
    for (index, (limit, witness)) in limits.iter().zip(witnesses).enumerate() {
        for key in ["index", "name", "kind", "unit", "sense", "coefficients"] {
            assert_eq!(witness.get(key), limit.get(key));
        }
        let coefficients = numbers(witness.get("coefficients").unwrap());
        let mut expected = offset.clone();
        for (component, value) in expected.iter_mut().enumerate() {
            for (column, coefficient) in basis.iter().zip(&coefficients) {
                *value += coefficient * column[component];
            }
        }
        let name = format!(
            "response_displacement_{index}_{}",
            witness.str_field("name").unwrap()
        );
        let array = grid.arrays.iter().find(|a| a.name == name).unwrap();
        assert_eq!(array.association, DataAssociation::PointData);
        assert_eq!(array.components, 3);
        assert_eq!(array.unit.as_deref(), Some("m"));
        assert_eq!(array.values, DataValues::Float64(expected));
    }
    assert!(
        grid.arrays
            .iter()
            .any(|a| a.name == "sampled_physical_stress_envelope")
    );
    assert!(
        grid.arrays
            .iter()
            .any(|a| a.name == "sampled_relaxed_stress_envelope")
    );
}
