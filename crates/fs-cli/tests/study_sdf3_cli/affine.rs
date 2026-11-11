//! Actual CLI optimization and exports for self-weight plus uncertain loads.
use super::*;
use fs_viz::vtu::{DataValues, VtuChecker};

const SOURCE: &str = include_str!("../../../../examples/marquee/bracket-3d-affine-load-space.fsim");

fn numbers(value: &J) -> Vec<f64> {
    value
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn affine_cli_reduces_feasible_material_and_exports_the_complete_physical_witness() {
    let dir = scratch("affine-load-space");
    // Leave enough physical headroom for a material reduction in two updates.
    // The report must still demonstrate feasibility against this explicit cap.
    let source = SOURCE.replace(
        ":compliance-limit-j 0.000009765625",
        ":compliance-limit-j 0.00009765625",
    );
    let source_path = dir.join("affine.fsim");
    let ledger = dir.join("affine.db");
    fs::write(&source_path, &source).unwrap();
    let result = document(
        &command("study")
            .arg(&source_path)
            .arg(&ledger)
            .args(["--budget", "2"])
            .output()
            .unwrap(),
        6,
    );
    assert_eq!(result.str_field("status"), Some("budget-exhausted"));
    assert_eq!(
        result
            .path(&["receipt", "iteration_limit_this_invocation"])
            .and_then(J::as_f64),
        Some(2.0)
    );
    assert_eq!(
        result.path(&["receipt", "resume_supported"]),
        Some(&J::Bool(true))
    );

    let report = data(&ledger, &result, "report_json");
    assert_eq!(
        report.str_field("driver"),
        Some("affine-load-space-sdf3-elasticity-study-v1")
    );
    assert_eq!(
        report.path(&["gradient_check", "passed"]),
        Some(&J::Bool(true))
    );
    assert_eq!(
        report.str_field("selected_design"),
        Some("least-material-feasible")
    );
    assert_eq!(report.get("selected_feasible"), Some(&J::Bool(true)));
    let selected = report.get("selected").unwrap();
    let cap = number(report.get("model").unwrap(), "compliance_limit_j");
    assert_eq!(cap, 0.1 / 1024.0);
    assert!(number(selected, "worst_compliance_j") > 0.0);
    assert!(number(selected, "worst_compliance_j") <= cap);
    assert!(number(selected, "maximum_ratio") <= 1.0);
    assert!(number(selected, "normalized_offset") > 0.0);
    assert!(number(selected, "witness_relative_residual") < 1e-7);
    let iterations = data(&ledger, &result, "iterations");
    let history = iterations.get("history").and_then(J::as_array).unwrap();
    assert!(
        history.len() >= 2,
        "the CLI must accept a real material update"
    );
    assert!(
        number(selected, "dilated_volume_fraction")
            < number(&history[0], "dilated_volume_fraction")
    );
    assert!(number(&report, "iterations_completed") <= 2.0);

    // Both the common offset and every uncertain column must survive retention.
    // Reconstruct the reported worst physical field from the retained witness.
    let design = data(&ledger, &result, "design");
    let fields = design.get("selected").unwrap();
    let coefficients = numbers(selected.get("worst_coefficients").unwrap());
    let basis: Vec<_> = fields
        .get("basis_displacements_m")
        .and_then(J::as_array)
        .unwrap()
        .iter()
        .map(numbers)
        .collect();
    let offset = numbers(fields.get("offset_displacement_m").unwrap());
    let worst = numbers(fields.get("worst_displacement_m").unwrap());
    assert_eq!(basis.len(), 2);
    assert_eq!(coefficients.len(), basis.len());
    assert_eq!(offset.len(), worst.len());
    assert!(offset.iter().any(|u| *u != 0.0));
    assert!(coefficients.iter().map(|q| q * q).sum::<f64>() <= 1.0 + 1e-10);
    let scale = worst.iter().map(|u| u.abs()).fold(1e-30_f64, f64::max);
    for field in &basis {
        assert_eq!(field.len(), worst.len());
    }
    for i in 0..worst.len() {
        let reconstructed = offset[i]
            + basis
                .iter()
                .zip(&coefficients)
                .map(|(u, q)| u[i] * q)
                .sum::<f64>();
        assert!(
            (worst[i] - reconstructed).abs() <= 1e-10 * scale,
            "affine displacement mismatch at component {i}"
        );
    }

    // A separate process must export from the retained design, including the
    // common offset, even when the authored study file has moved away.
    fs::rename(&source_path, dir.join("relocated.fsim")).unwrap();
    let pointer = result.str_field("run_id").unwrap();
    let exported = document(
        &command("report")
            .arg(pointer)
            .arg(&ledger)
            .output()
            .unwrap(),
        0,
    );
    for key in ["report_html", "report_json", "design", "iterations"] {
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
        ("offset_displacement", &offset),
        ("worst_displacement", &worst),
    ] {
        let array = grid.arrays.iter().find(|array| array.name == name).unwrap();
        let DataValues::Float64(actual) = &array.values else {
            panic!("physical displacement array");
        };
        assert_eq!(actual, expected);
    }
    let package = document(
        &command("package")
            .arg(pointer)
            .arg(&ledger)
            .output()
            .unwrap(),
        0,
    );
    assert_eq!(package.str_field("checker"), Some("pass"));
    assert_eq!(
        fs::read(package.str_field("package").unwrap()).unwrap(),
        retained(&ledger, &result, "package")
    );

    let before = retained(&ledger, &result, "design");
    let resumed_output = command("study")
        .args(["--resume", pointer])
        .arg(&ledger)
        .args(["--budget", "2"])
        .output()
        .unwrap();
    assert!(matches!(resumed_output.status.code(), Some(0 | 6)),
        "{}", String::from_utf8_lossy(&resumed_output.stderr));
    let resumed = document(&resumed_output, resumed_output.status.code().unwrap());
    let work = resumed.path(&["receipt", "optimizer_work"]).unwrap();
    assert!(number(work, "iterations") > number(&report, "iterations_completed"));
    assert!(number(work, "iterations") <= 4.0);
    assert!(number(work, "restoration_evaluations") >= 1.0);
    let whole_ledger = dir.join("whole.db");
    let whole_output = command("study")
        .arg(dir.join("relocated.fsim"))
        .arg(&whole_ledger)
        .args(["--budget", "4"])
        .output()
        .unwrap();
    assert_eq!(whole_output.status.code(), resumed_output.status.code());
    let whole = document(&whole_output, whole_output.status.code().unwrap());
    for key in ["design", "iterations"] {
        assert_eq!(retained(&ledger, &resumed, key), retained(&whole_ledger, &whole, key));
    }
    assert_eq!(resumed.path(&["receipt", "trace_hash"]), whole.path(&["receipt", "trace_hash"]));
    assert_eq!(retained(&ledger, &result, "design"), before);
}

#[test]
fn affine_cli_refuses_missing_self_weight_and_probability_weights_before_ledger_creation() {
    let dir = scratch("affine-admission");
    let prefix = SOURCE.rsplit_once("  (material-acceleration").unwrap().0;
    let no_weight_law = format!("{prefix})\n");
    let probability_weight = SOURCE.replacen(":weight 1.0", ":weight 0.5", 1);
    for (name, source) in [
        ("missing-self-weight", no_weight_law),
        ("nonunit-basis", probability_weight),
    ] {
        let source_path = dir.join(format!("{name}.fsim"));
        let ledger = dir.join(format!("{name}.db"));
        fs::write(&source_path, source).unwrap();
        let refused = command("study")
            .arg(&source_path)
            .arg(&ledger)
            .args(["--budget", "2"])
            .output()
            .unwrap();
        assert_eq!(refused.status.code(), Some(4));
        assert!(String::from_utf8_lossy(&refused.stderr).contains("cli-study-sdf3-input"));
        assert!(
            !ledger.exists(),
            "{name} must refuse before creating a ledger"
        );
    }
}
