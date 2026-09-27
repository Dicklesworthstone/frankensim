//! Real binary and on-disk ledger regressions, not a surrogate optimizer.
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-2d.fsim"));

fn scratch(label: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-native-study-{label}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}
fn source(dir: &Path) -> PathBuf {
    let path = dir.join("study.fsim");
    fs::write(&path, FIXTURE.replace(":mesh-level 4", ":mesh-level 2")
        .replace(":max-iterations 32", ":max-iterations 3")
        .replace(":steps 32", ":steps 3")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":nucleation-period 4", ":nucleation-period 0")).unwrap();
    path
}
fn command(verb: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", verb]);
    command
}
fn document(output: &Output, exit: u8) -> J {
    assert_eq!(output.status.code(), Some(i32::from(exit)), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run_id(result: &J) -> &str { result.str_field("run_id").unwrap() }
/// Admission refusals print no result on stdout; the typed diagnostic is the
/// single JSON line on stderr.
fn refusal(output: &Output, exit: u8) -> J {
    assert_eq!(output.status.code(), Some(i32::from(exit)), "{}", String::from_utf8_lossy(&output.stderr));
    assert!(output.stdout.is_empty(), "a refusal publishes no result: {}", String::from_utf8_lossy(&output.stdout));
    let diagnostic = J::parse(std::str::from_utf8(&output.stderr).unwrap().trim()).unwrap();
    assert_eq!(diagnostic.str_field("schema"), Some("frankensim.cli.diagnostic.v1"));
    assert_eq!(diagnostic.str_field("severity"), Some("error"));
    assert!(diagnostic.str_field("code").is_some_and(|code| code.starts_with("cli-study")));
    diagnostic
}
fn updates(result: &J) -> f64 {
    result.path(&["receipt", "continuation", "updates_this_invocation"]).and_then(J::as_f64).unwrap()
}
fn retained(ledger: &Path, result: &J, key: &str) -> Vec<u8> {
    let hash = result.path(&["receipt", key]).and_then(J::as_str).and_then(ContentHash::from_hex).unwrap();
    Ledger::open(ledger.to_str().unwrap()).unwrap()
        .get_artifact_bounded(&hash, 16 * 1024 * 1024).unwrap().unwrap()
}

#[test]
fn native_study_chunks_keep_identical_geometry_and_rows_then_export_without_solving() {
    let dir = scratch("chunks");
    let source = source(&dir);
    let whole_db = dir.join("whole.db");
    let chunk_db = dir.join("chunk.db");
    let whole = document(&command("study").arg(&source).arg(&whole_db).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(whole.str_field("status"), Some("constraint-unmet"));
    assert_eq!(updates(&whole), 3.0);
    let first = document(&command("study").arg(&source).arg(&chunk_db)
        .args(["--budget", "1"]).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(first.str_field("status"), Some("budget-exhausted"));
    assert_eq!(updates(&first), 1.0);
    let first_design = retained(&chunk_db, &first, "design");
    let resumed = document(&command("study").arg("--resume").arg(run_id(&first)).arg(&chunk_db)
        .args(["--budget", "2"]).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(resumed.str_field("status"), Some("constraint-unmet"));
    assert_eq!(updates(&resumed), 2.0);
    assert_eq!(resumed.path(&["receipt", "continuation", "legacy_prefix_updates_replayed"]).and_then(J::as_f64), Some(0.0));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&whole_db, &whole, key), retained(&chunk_db, &resumed, key));
    }
    assert_eq!(retained(&chunk_db, &first, "design"), first_design);
    assert_eq!(whole.path(&["receipt", "trace_hash"]), resumed.path(&["receipt", "trace_hash"]));
    for verb in ["report", "package"] {
        let exported = document(&command(verb).arg(run_id(&resumed)).arg(&chunk_db)
            .output().unwrap(), fs_cli::exit::SUCCESS);
        assert_eq!(exported.str_field("status"), Some("ok"));
        assert_eq!(exported.str_field("study_status"), Some("constraint-unmet"));
    }
    let again = document(&command("study").arg("--resume").arg(run_id(&resumed)).arg(&chunk_db)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(again, resumed, "a terminal receipt is returned unchanged, including its unmet constraint");
}

#[test]
fn infeasible_material_target_cannot_report_a_successful_study() {
    // Three short updates from two 0.12-radius holes cannot move the area
    // from ~0.91 to 0.45. Both the public status/exit and retained report
    // must distinguish iteration completion from material feasibility.
    let dir = scratch("feasibility");
    let source = source(&dir);
    let database = dir.join("study.db");
    let result = document(&command("study").arg(&source).arg(&database)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(result.str_field("status"), Some("constraint-unmet"));
    let report = J::parse(std::str::from_utf8(&retained(&database, &result, "report_json")).unwrap()).unwrap();
    let area = report.path(&["final_material_area_m2"]).and_then(J::as_f64).unwrap();
    let target = report.path(&["volume_target_m2"]).and_then(J::as_f64).unwrap();
    let violation = report.path(&["volume_violation_m2"]).and_then(J::as_f64).unwrap();
    let tolerance = report.path(&["volume_tolerance_m2"]).and_then(J::as_f64).unwrap();
    assert_eq!(target, 0.45, "unit box times the fixture's volume fraction");
    assert_eq!(violation, area - target);
    assert_eq!(tolerance, 0.01);
    assert!(violation > tolerance, "fixture must end infeasible: area {area}");
    assert_eq!(report.path(&["volume_constraint_satisfied"]), Some(&J::Bool(false)));
    assert_eq!(report.str_field("status"), Some("constraint-unmet"));
    assert_eq!(report.f64_field("iterations_completed"), Some(3.0));
    let html = String::from_utf8(retained(&database, &result, "report_html")).unwrap();
    assert!(html.contains("CONSTRAINT NOT SATISFIED") && html.contains("constraint-unmet"));
}

#[test]
fn public_resume_finalizes_a_committed_final_design_without_repeating_an_update() {
    let dir = scratch("finalize");
    let source = source(&dir);
    let database = dir.join("study.db");
    let complete = document(&command("study").arg(&source).arg(&database)
        .output().unwrap(), fs_cli::exit::BUDGET);
    let predecessor = complete.path(&["receipt", "predecessor"]).and_then(J::as_str).unwrap();
    let pointer = format!("study-{predecessor}");
    let recovered = document(&command("study").arg("--resume").arg(&pointer).arg(&database)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(updates(&recovered), 0.0);
    assert_eq!(recovered.str_field("status"), Some("constraint-unmet"));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&database, &complete, key), retained(&database, &recovered, key));
    }
}

#[path = "study_checkpoint_cli/projected.rs"]
mod projected;

#[path = "study_checkpoint_cli/multi_load.rs"]
mod multi_load;

#[test]
fn elasticity_objective_twins_respond_to_the_load_and_the_material_budget() {
    // q61wp.16 objective-sensitivity twins, run in the projected-volume mode:
    // the plain augmented-Lagrangian mode does not converge its area
    // (feasibility flips with the step count), so its endpoints are not
    // comparable. Spreading the traction band (v1 admits only bands symmetric
    // about y = 0.5) must change the accepted design. A larger material budget
    // must end stiffer at its own area. A study that ignored its load or its
    // volume constraint fails both.
    const PROJECTED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
    let eight = PROJECTED.replace(":steps 2", ":steps 8").replace(":max-iterations 2", ":max-iterations 8");
    assert!(eight.contains(":load-band (0.375 0.625)") && eight.contains(":volume-fraction 0.75"));
    let dir = scratch("objective-twins");
    // A twin may stop early with no-feasible-descent (an honest terminal
    // status) but must have accepted updates at its own area.
    let run = |name: &str, text: String| -> J {
        let path = dir.join(format!("{name}.fsim"));
        fs::write(&path, text).unwrap();
        let output = command("study").arg(&path).arg(dir.join(format!("{name}.db"))).output().unwrap();
        let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
        assert!(
            matches!(result.str_field("status"), Some("completed" | "no-feasible-descent")),
            "{name}: {}", String::from_utf8_lossy(&output.stderr)
        );
        result
    };
    let last = |result: &J| -> (f64, f64) {
        let accepted = result.path(&["receipt", "continuation", "constraints", "accepted"]).unwrap().as_array().unwrap();
        let row = accepted.last().expect("accepted updates");
        (row.f64_field("area_m2").unwrap(), row.f64_field("compliance_j").unwrap())
    };
    let design = |result: &J| result.path(&["receipt", "design"]).and_then(|v| v.as_str()).unwrap().to_string();
    let base = run("base", eight.clone());
    let wide = run("wide-band", eight.replace(":load-band (0.375 0.625)", ":load-band (0.25 0.75)"));
    let richer = run("richer", eight.replace(":volume-fraction 0.75", ":volume-fraction 0.85"));
    let ((base_area, base_j), (wide_area, _), (richer_area, richer_j)) = (last(&base), last(&wide), last(&richer));
    for (area, target) in [(base_area, 0.75), (wide_area, 0.75), (richer_area, 0.85)] {
        assert!((area - target).abs() <= 1e-3, "area {area} vs target {target}");
    }
    assert_ne!(design(&base), design(&wide), "a redistributed traction band must change the accepted design");
    assert!(richer_j < base_j, "more material must end stiffer: {richer_j} vs {base_j}");
    println!("{{\"base_j\":{base_j},\"richer_j\":{richer_j}}}");
}


#[test]
fn final_dwr_assesses_the_retained_design_and_resume_reuses_completed_work() {
    let dir = scratch("final-dwr");
    let path = dir.join("assessed.fsim");
    let source = FIXTURE.replace(":mesh-level 4", ":mesh-level 3")
        .replace(":max-iterations 32", ":max-iterations 1")
        .replace(":steps 32", ":steps 1")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":nucleation-period 4", ":nucleation-period 0");
    fs::write(&path, format!(
        "{}\n  (assessment :type elasticity-dwr :max-solves-per-attempt 2)\n)\n",
        source.trim_end().strip_suffix(')').unwrap(),
    )).unwrap();
    let database = dir.join("assessed.db");
    let result = document(&command("study").arg(&path).arg(&database)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(result.str_field("status"), Some("constraint-unmet"));
    let summary = J::parse(std::str::from_utf8(
        &retained(&database, &result, "report_json")).unwrap()).unwrap();
    let dwr = summary.get("goal_error_assessment").unwrap();
    assert_eq!(dwr.str_field("status"), Some("estimated"));
    assert_eq!(dwr.str_field("snapshot"), summary.str_field("snapshot"));
    assert_eq!(dwr.f64_field("coarse_compliance_j"), summary.f64_field("final_compliance_j"));
    assert_eq!(dwr.f64_field("material_area_m2"), summary.f64_field("final_material_area_m2"));
    assert_eq!(dwr.f64_field("enriched_level"), Some(4.0));
    assert!(dwr.f64_field("enriched_dofs").unwrap() > dwr.f64_field("coarse_dofs").unwrap());
    for key in ["eta_signed_j", "absolute_indicator_sum_j", "enriched_compliance_j"] {
        assert!(dwr.f64_field(key).unwrap().is_finite());
    }
    let terms = dwr.get("residual_terms_j").unwrap();
    let parts: Vec<f64> = ["bulk", "nitsche", "outer_traction", "ghost"]
        .iter().map(|key| terms.f64_field(key).unwrap()).collect();
    let total: f64 = parts.iter().sum();
    let scale: f64 = parts.iter().map(|value| value.abs()).sum();
    assert!((total - dwr.f64_field("eta_signed_j").unwrap()).abs()
        <= 1e-8 * scale.max(1e-8));
    assert_eq!(dwr.path(&["solver", "relative_residual_kind"]).and_then(J::as_str),
        Some("recomputed-euclidean"));
    for key in ["coarse_relative_residual", "enriched_relative_residual"] {
        assert!(dwr.path(&["solver", key]).and_then(J::as_f64).unwrap() < 1e-12);
    }
    let html = String::from_utf8(retained(&database, &result, "report_html")).unwrap();
    assert!(html.contains("Final compliance goal-error estimate")
        && html.contains("not certified continuum-error bounds"));

    // The final optimizer update was already durable before DWR began. Model
    // losing the assessment publication by resuming that exact predecessor:
    // assessment repeats, but no geometry update or trajectory replay occurs.
    let predecessor = result.path(&["receipt", "predecessor"]).and_then(J::as_str).unwrap();
    let resumed = document(&command("study").arg("--resume").arg(format!("study-{predecessor}"))
        .arg(&database).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(updates(&resumed), 0.0);
    for key in ["design", "iterations", "report_json"] {
        assert_eq!(retained(&database, &result, key), retained(&database, &resumed, key));
    }
    let again = document(&command("study").arg("--resume").arg(run_id(&resumed))
        .arg(&database).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(again, resumed, "a completed assessment is returned without another solve");
}
