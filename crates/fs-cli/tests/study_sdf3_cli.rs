//! G1/G3/G4/G5: actual CLI, implicit quadrature, elasticity, adaptive SIMP and ledger.
#![cfg(feature = "sdf3-study")]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use json::JsonValue as J;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const FIXTURE: &str = include_str!("../../../examples/marquee/bracket-3d-adaptive.fsim");

fn scratch(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "fs-sdf3-cli-{label}-{}-{nonce}",
        std::process::id()
    ));
    fs::create_dir(&path).unwrap();
    path
}
fn command(verb: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", verb]);
    command
}
fn document(output: &Output, exit: i32) -> J {
    assert_eq!(
        output.status.code(),
        Some(exit),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run(dir: &Path, name: &str, text: &str, exit: i32) -> (PathBuf, J) {
    let source = dir.join(format!("{name}.fsim"));
    let ledger = dir.join(format!("{name}.db"));
    fs::write(&source, text).unwrap();
    let result = document(
        &command("study").arg(source).arg(&ledger).output().unwrap(),
        exit,
    );
    (ledger, result)
}
fn retained(ledger: &Path, result: &J, name: &str) -> Vec<u8> {
    let hash = result
        .path(&["receipt", name])
        .and_then(J::as_str)
        .and_then(ContentHash::from_hex)
        .unwrap();
    Ledger::open(ledger.to_str().unwrap())
        .unwrap()
        .get_artifact_bounded(&hash, 16 * 1024 * 1024)
        .unwrap()
        .unwrap()
}
fn data(ledger: &Path, result: &J, name: &str) -> J {
    J::parse(std::str::from_utf8(&retained(ledger, result, name)).unwrap()).unwrap()
}
fn number(value: &J, key: &str) -> f64 {
    value.get(key).and_then(J::as_f64).unwrap()
}

#[test]
fn g1_g5_sdf3_cli_refines_checks_gradients_and_exports_the_solved_density_and_fields() {
    let dir = scratch("adaptive");
    let source = FIXTURE.replace(":updates-per-stage 3", ":updates-per-stage 1");
    let (ledger, result) = run(&dir, "adaptive", &source, 0);
    assert_eq!(result.str_field("status"), Some("completed"));
    let report = data(&ledger, &result, "report_json");
    let iterations = data(&ledger, &result, "iterations");
    let stages = iterations.get("stages").and_then(J::as_array).unwrap();
    assert_eq!(stages.len(), 3);
    for stage in stages {
        assert_eq!(stage.get("gradient_check_passed"), Some(&J::Bool(true)));
        let probes = stage.get("gradient_probes").and_then(J::as_array).unwrap();
        assert_eq!(probes.len(), 2);
        for probe in probes {
            assert!(number(probe, "compliance_relative_error") <= 5e-4);
            assert!(number(probe, "volume_relative_error") <= 5e-4);
        }
        let history = stage.get("history").and_then(J::as_array).unwrap();
        assert!(
            history.len() >= 2,
            "a real update follows every new solved baseline"
        );
        for row in history {
            assert!(number(row, "volume_fraction") <= 0.50000001);
            let cases = row.get("case_compliances_j").and_then(J::as_array).unwrap();
            let independent = 0.3 * cases[0].as_f64().unwrap() + 0.7 * cases[1].as_f64().unwrap();
            assert!((independent - number(row, "compliance_j")).abs() <= independent * 1e-12);
        }
        for rows in history.windows(2) {
            assert!(
                number(&rows[1], "compliance_j")
                    <= number(&rows[0], "compliance_j") * (1.0 + 1e-12)
            );
        }
    }
    let refinements = report.get("refinements").and_then(J::as_array).unwrap();
    assert_eq!(refinements.len(), 2);
    for refinement in refinements {
        assert_eq!(refinement.get("installed"), Some(&J::Bool(true)));
        assert!(
            number(refinement, "target_active_cells") > number(refinement, "source_active_cells")
        );
        for case in refinement.get("cases").and_then(J::as_array).unwrap() {
            assert!(number(case, "identity_relative_defect") <= 1e-8);
            for residual in case.get("field_residuals").and_then(J::as_array).unwrap() {
                assert!(residual.as_f64().unwrap() <= 1e-8);
            }
        }
    }
    let design = data(&ledger, &result, "design");
    let cells = design.get("cells").and_then(J::as_array).unwrap();
    assert!(cells.len() > 8);
    assert_eq!(cells.len() as f64, number(&report, "active_cells"));
    let volume: f64 = cells.iter().map(|cell| number(cell, "cut_volume_m3")).sum();
    let material: f64 = cells
        .iter()
        .map(|cell| number(cell, "cut_volume_m3") * number(cell, "projected_density"))
        .sum();
    assert!((material / volume - number(&report, "final_volume_fraction")).abs() < 1e-12);
    let nodes = design
        .get("physical_nodes_m")
        .and_then(J::as_array)
        .unwrap();
    let fields = design.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(
        fields[0], fields[1],
        "independent load cases have independent displacement fields"
    );
    for field in fields {
        let values = field.as_array().unwrap();
        assert_eq!(values.len(), nodes.len() * 3);
        assert!(
            values
                .iter()
                .any(|value| value.as_f64().unwrap().abs() > 0.0)
        );
        for (index, node) in nodes.iter().enumerate() {
            if node.as_array().unwrap()[0].as_f64() == Some(0.0) {
                for value in &values[3 * index..3 * index + 3] {
                    assert_eq!(value.as_f64(), Some(0.0));
                }
            }
        }
    }
    let again = run(&dir, "repeated", &source, 0);
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
    document(
        &command("package").arg(id).arg(&ledger).output().unwrap(),
        0,
    );
}

#[test]
fn g4_sdf3_leaf_budget_retains_the_exact_solved_coarse_endpoint() {
    let dir = scratch("leaf-budget");
    let source = FIXTURE.replace(":updates-per-stage 3", ":updates-per-stage 1");
    let coarse = source.replace(
        ":schedule ((1.0 1.0) (2.0 2.0) (3.0 8.0))",
        ":schedule ((1.0 1.0))",
    );
    let (coarse_db, baseline) = run(&dir, "coarse", &coarse, 0);
    let (partial_db, partial) = run(
        &dir,
        "partial",
        &source.replace(":maximum-leaves 2048", ":maximum-leaves 8"),
        6,
    );
    assert_eq!(partial.str_field("status"), Some("budget-exhausted"));
    assert_eq!(
        retained(&coarse_db, &baseline, "design"),
        retained(&partial_db, &partial, "design")
    );
    assert_eq!(
        retained(&coarse_db, &baseline, "iterations"),
        retained(&partial_db, &partial, "iterations")
    );
    let summary = data(&partial_db, &partial, "report_json");
    assert_eq!(number(&summary, "active_cells"), 8.0);
    assert!(
        summary
            .path(&["stop", "refinement"])
            .and_then(J::as_str)
            .unwrap()
            .contains("LeafBudget")
    );
    let pointer = partial.str_field("run_id").unwrap();
    let exported = document(
        &command("report")
            .arg(pointer)
            .arg(&partial_db)
            .output()
            .unwrap(),
        0,
    );
    assert_eq!(exported.str_field("study_status"), Some("budget-exhausted"));
    let refused = command("study")
        .arg("--resume")
        .arg(pointer)
        .arg(&partial_db)
        .output()
        .unwrap();
    let resumed = document(&refused, 6);
    for key in ["design", "iterations", "checkpoint"] {
        assert_eq!(retained(&partial_db, &partial, key), retained(&partial_db, &resumed, key));
    }
    assert!(resumed.path(&["receipt", "work", "linear_iterations"]).and_then(J::as_f64).unwrap()
        > partial.path(&["receipt", "work", "linear_iterations"]).and_then(J::as_f64).unwrap());
    assert_eq!(resumed.path(&["receipt", "stages_replayed"]).and_then(J::as_f64), Some(1.0));
    assert_eq!(
        retained(&coarse_db, &baseline, "design"),
        retained(&partial_db, &partial, "design")
    );
}

#[test]
fn g4_sdf3_exhausted_initial_solve_exports_no_unassessed_design() {
    let dir = scratch("linear-budget");
    let (ledger, result) = run(
        &dir,
        "stopped",
        &FIXTURE.replace(":linear-iterations 250000", ":linear-iterations 1"),
        6,
    );
    let design = data(&ledger, &result, "design");
    assert_eq!(design.str_field("state"), Some("no-solved-baseline"));
    assert!(
        design
            .get("cells")
            .and_then(J::as_array)
            .unwrap()
            .is_empty()
    );
    assert!(
        design
            .get("displacements_m")
            .and_then(J::as_array)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        result
            .path(&["receipt", "work", "linear_iterations"])
            .and_then(J::as_f64),
        Some(1.0)
    );
}

#[path = "study_sdf3_cli/checkpoint.rs"]
mod checkpoint;
