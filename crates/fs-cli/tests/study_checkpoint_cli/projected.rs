//! Use the same executable/ledger helpers as the existing continuation battery.
use super::*;

const CONSTRAINED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-projected-stress-2d.fsim"));

fn source(dir: &Path, text: &str) -> PathBuf {
    let path = dir.join("projected.fsim");
    fs::write(&path, text).unwrap();
    path
}
fn constraints(result: &J) -> &J {
    result.path(&["receipt", "continuation", "constraints"]).expect("native constraint history")
}

#[test]
fn projected_native_command_retains_constraints_through_disk_resume_and_report_export() {
    let dir = scratch("projected-chunks");
    let input = source(&dir, CONSTRAINED);
    let whole_db = dir.join("whole.db");
    let chunks_db = dir.join("chunks.db");
    let first = document(&command("study").arg(&input).arg(&chunks_db)
        .args(["--budget", "1"]).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(first.str_field("status"), Some("budget-exhausted"));
    assert_eq!(updates(&first), 1.0, "fixture must actually accept a geometry update");
    let state = constraints(&first);
    assert_eq!(state.str_field("mode"), Some("projected-stress-v1"));
    assert_eq!(state.str_field("baseline_scope"), Some("feasible_study_start"));
    let baseline = state.get("baseline").unwrap();
    let accepted = state.get("accepted").and_then(J::as_array).unwrap();
    assert_eq!(accepted.len(), 1);
    assert!(accepted[0].f64_field("compliance_j").unwrap() < baseline.f64_field("compliance_j").unwrap());
    assert!((accepted[0].f64_field("area_m2").unwrap() - 0.75).abs() <= 1e-4);
    assert!(accepted[0].f64_field("sampled_von_mises_pa").unwrap() <= 1e12);
    let original_design = retained(&chunks_db, &first, "design");
    let whole_output = command("study").arg(&input).arg(&whole_db).output().unwrap();
    // A bounded search may stop at the second step, but may not fake the first.
    let code = whole_output.status.code().unwrap();
    assert!(code == i32::from(fs_cli::exit::SUCCESS) || code == i32::from(fs_cli::exit::REFUSED),
        "stdout={} stderr={}", String::from_utf8_lossy(&whole_output.stdout), String::from_utf8_lossy(&whole_output.stderr));
    let whole = document(&whole_output, u8::try_from(code).unwrap());
    assert!(matches!(whole.str_field("status"), Some("completed" | "no-feasible-descent")));
    let resumed = document(&command("study").arg("--resume").arg(run_id(&first)).arg(&chunks_db)
        .output().unwrap(), u8::try_from(code).unwrap());
    assert_eq!(resumed.str_field("status"), whole.str_field("status"));
    assert_eq!(constraints(&resumed), constraints(&whole));
    assert_eq!(resumed.path(&["receipt", "continuation", "legacy_prefix_updates_replayed"])
        .and_then(J::as_f64), Some(0.0));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&whole_db, &whole, key), retained(&chunks_db, &resumed, key));
    }
    assert_eq!(retained(&chunks_db, &first, "design"), original_design);
    for verb in ["report", "package"] {
        let exported = document(&command(verb).arg(run_id(&resumed)).arg(&chunks_db)
            .output().unwrap(), fs_cli::exit::SUCCESS);
        assert_eq!(exported.str_field("study_status"), resumed.str_field("status"));
        assert_eq!(exported.str_field("verification"), Some("sealed-evidence"));
    }
    let summary_path = dir.join(format!("{}.json", run_id(&resumed)));
    let summary_bytes = fs::read(&summary_path).unwrap();
    assert_eq!(summary_bytes, retained(&chunks_db, &resumed, "report_json"));
    let summary = J::parse(std::str::from_utf8(&summary_bytes).unwrap()).unwrap();
    assert_eq!(summary.get("constraints"), Some(constraints(&resumed)));
    let html = fs::read_to_string(dir.join(format!("{}.html", run_id(&resumed)))).unwrap();
    assert!(html.contains("Hard material area"));
    assert!(html.contains("sample-scoped, not a continuous-domain bound"));
}

#[test]
fn projected_public_admission_refuses_missing_controls_before_creating_a_ledger() {
    for (index, bad) in [
        CONSTRAINED.replace("    :max-area-evaluations 64\n", ""),
        CONSTRAINED.replace(":max-candidates 16", ":max-candidates 0"),
        CONSTRAINED.replace(":constraint-mode projected-stress",
            ":constraint-mode projected-stress :constraint-mode projected-stress"),
    ].into_iter().enumerate() {
        let dir = scratch(&format!("projected-refusal-{index}"));
        let input = source(&dir, &bad);
        let db = dir.join("must-not-exist.db");
        refusal(&command("study").arg(&input).arg(&db).output().unwrap(), fs_cli::exit::REFUSED);
        assert!(!db.exists(), "malformed policy must refuse before opening persistent state");
    }
}

#[test]
fn projected_public_stall_is_exportable_but_is_never_a_completed_optimization() {
    let dir = scratch("projected-stall");
    let input = source(&dir, &CONSTRAINED
        .replace(":min-relative-improvement 0.00000001", ":min-relative-improvement 0.999999")
        .replace(":max-candidates 16", ":max-candidates 2"));
    let db = dir.join("stall.db");
    let result = document(&command("study").arg(&input).arg(&db).output().unwrap(), fs_cli::exit::REFUSED);
    assert_eq!(result.str_field("status"), Some("no-feasible-descent"));
    assert_eq!(updates(&result), 0.0);
    let constraints = constraints(&result);
    assert!(constraints.get("accepted").and_then(J::as_array).unwrap().is_empty());
    assert_eq!(constraints.get("terminal_refusals").and_then(J::as_array).unwrap().len(), 2);
    let report = document(&command("report").arg(run_id(&result)).arg(&db).output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(report.str_field("study_status"), Some("no-feasible-descent"));
    let again = document(&command("study").arg("--resume").arg(run_id(&result)).arg(&db)
        .output().unwrap(), fs_cli::exit::REFUSED);
    assert_eq!(again, result, "stalled terminal must not silently enlarge its candidate search");
}

#[test]
fn native_region_example_keeps_exact_prescriptions_after_source_independent_disk_resume() {
    use fs_topols::design_regions::{DesignPhase, DesignRegion, prepare_design_regions};
    use fs_topols::GridSdf;
    const AUTHORED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-protected-regions-2d.fsim"));
    let initial = GridSdf::from_fn(8, &|x, y| {
        [0.35, 0.65].iter().map(|&cx| 0.12 - (x - cx).hypot(y - 0.5)).fold(-1.0, f64::max)
    });
    let boundary: Vec<_> = initial.nodes().iter().copied().enumerate()
        .filter(|(i, _)| i % 9 == 0 || i % 9 == 8).collect();
    let prescribed = prepare_design_regions(&initial, &boundary, &[
        DesignRegion::new(DesignPhase::Material, [0.125, 0.125], [0.25, 0.25], 0.01).unwrap(),
        DesignRegion::new(DesignPhase::Void, [0.3, 0.4], [0.32, 0.42], 0.01).unwrap(),
    ]).unwrap();
    assert!(prescribed.changed_nodes > 0);
    assert_eq!(prescribed.fixed_nodes.len(), 26);
    let check_geometry = |db: &Path, result: &J| {
        let bytes = retained(db, result, "design");
        let design = J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
        assert_eq!(design.f64_field("n"), Some(8.0));
        let bits = design.get("phi_bits").and_then(J::as_array).unwrap();
        assert_eq!(bits.len(), 81);
        for &(i, expected) in &prescribed.fixed_nodes {
            let actual = u64::from_str_radix(bits[i].as_str().unwrap(), 16).unwrap();
            assert_eq!(actual, expected.to_bits(), "protected node {i}");
        }
        assert_eq!(constraints(result).get("design_regions").and_then(J::as_array).unwrap().len(), 2);
    };
    let dir = scratch("native-regions");
    let input = source(&dir, AUTHORED);
    let db = dir.join("regions.db");
    let whole_db = dir.join("whole.db");
    let first = document(&command("study").arg(&input).arg(&db)
        .args(["--budget", "1"]).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(updates(&first), 1.0, "the authored example must accept an actual update");
    check_geometry(&db, &first);
    let original_design = retained(&db, &first, "design");
    let whole_output = command("study").arg(&input).arg(&whole_db).output().unwrap();
    let code = u8::try_from(whole_output.status.code().unwrap()).unwrap();
    assert!(code == fs_cli::exit::SUCCESS || code == fs_cli::exit::REFUSED,
        "{}", String::from_utf8_lossy(&whole_output.stderr));
    let whole = document(&whole_output, code);
    assert!(matches!(whole.str_field("status"), Some("completed" | "no-feasible-descent")));
    // Keep the test input, but prove recovery no longer depends on its path.
    fs::rename(&input, dir.join("original-input.fsim")).unwrap();
    let resumed = document(&command("study").arg("--resume").arg(run_id(&first)).arg(&db)
        .output().unwrap(), code);
    check_geometry(&db, &resumed);
    assert_eq!(constraints(&whole), constraints(&resumed));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&whole_db, &whole, key), retained(&db, &resumed, key));
    }
    assert_eq!(retained(&db, &first, "design"), original_design);
    let report = document(&command("report").arg(run_id(&resumed)).arg(&db)
        .output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(report.str_field("study_status"), resumed.str_field("status"));
    let summary = J::parse(&fs::read_to_string(dir.join(format!("{}.json", run_id(&resumed)))).unwrap()).unwrap();
    assert_eq!(summary.get("constraints"), Some(constraints(&resumed)));
    let html = fs::read_to_string(dir.join(format!("{}.html", run_id(&resumed)))).unwrap();
    assert!(html.contains("protected material/void regions"));
    assert!(html.contains("not a certified physical clearance"));
}
