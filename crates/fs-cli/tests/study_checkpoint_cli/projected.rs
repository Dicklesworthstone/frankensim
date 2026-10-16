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

#[test]
fn a_binding_stress_limit_changes_the_trajectory_exactly_where_the_loose_study_exceeds_it() {
    // The tracked example's 1e12 Pa limit never binds (q61wp.75 item 2). At
    // 1.626 Pa over six updates the loose study's fifth accepted design reaches
    // 1.6269 Pa (measured 2026-09-25). The constrained study must match it
    // until then, refuse that candidate, and keep every accepted design under
    // the limit. At 1.61 Pa no update is admissible and the study says so.
    let dir = scratch("projected-binding");
    let six = CONSTRAINED.replace(":steps 2", ":steps 6").replace(":max-iterations 2", ":max-iterations 6");
    assert_ne!(six, CONSTRAINED);
    let run = |limit: &str, name: &str| -> J {
        let text = six.replace(":sampled-stress-limit-pa 1000000000000.0", &format!(":sampled-stress-limit-pa {limit}"));
        let path = dir.join(format!("{name}.fsim"));
        fs::write(&path, text).unwrap();
        let output = command("study").arg(&path).arg(dir.join(format!("{name}.db"))).output().unwrap();
        J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
    };
    let stress = |result: &J| -> Vec<f64> {
        constraints(result).get("accepted").unwrap().as_array().unwrap().iter()
            .map(|row| row.f64_field("sampled_von_mises_pa").unwrap()).collect()
    };
    let loose = run("1000000000000.0", "loose");
    let bound = run("1.626", "bound");
    assert_eq!(bound.str_field("status"), Some("completed"));
    let (loose_vm, bound_vm) = (stress(&loose), stress(&bound));
    assert_eq!(bound_vm.len(), 6);
    let first_excess = loose_vm.iter().position(|&vm| vm > 1.626).expect("the loose study exceeds the limit");
    assert!(first_excess > 0, "some updates are admissible under the limit");
    assert_eq!(&loose_vm[..first_excess], &bound_vm[..first_excess], "identical prefix before the limit binds");
    assert_ne!(loose_vm[first_excess], bound_vm[first_excess], "the binding limit changes the accepted design");
    assert!(bound_vm.iter().all(|&vm| vm <= 1.626), "{bound_vm:?}");
    let tight = run("1.61", "tight");
    assert_eq!(tight.str_field("status"), Some("no-feasible-descent"));
    assert!(stress(&tight).is_empty());
}


#[test]
fn projected_volume_dwr_assesses_the_feasible_endpoint_and_reuses_terminal_results() {
    // G2/G4: assess both an accepted endpoint and a feasible baseline when the
    // bounded search accepts no update, then reuse that result.
    const PROJECTED: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-projected-volume-2d.fsim"));
    let one = PROJECTED.replace(":steps 2", ":steps 1")
        .replace(":max-iterations 2", ":max-iterations 1");
    let assessed = format!(
        "{}\n  (assessment :type elasticity-dwr :max-solves-per-attempt 2)\n)\n",
        one.trim_end().strip_suffix(')').unwrap());
    for (label, status, exit, count) in [
        ("accepted", "completed", fs_cli::exit::SUCCESS, 1),
        ("stalled", "no-feasible-descent", fs_cli::exit::REFUSED, 0),
    ] {
        let dir = scratch(&format!("projected-dwr-{label}"));
        let text = if count == 0 {
            assessed.replace(":min-relative-improvement 0.00000001",
                ":min-relative-improvement 0.999999")
                .replace(":max-candidates 16", ":max-candidates 2")
        } else { assessed.clone() };
        let input = source(&dir, &text);
        let db = dir.join("assessed.db");
        let result = document(&command("study").arg(&input).arg(&db).output().unwrap(), exit);
        assert_eq!(result.str_field("status"), Some(status));
        assert_eq!(updates(&result), if count == 0 { 0.0 } else { 1.0 });
        let state = constraints(&result);
        let accepted = state.get("accepted").and_then(J::as_array).unwrap();
        assert_eq!(accepted.len(), count);
        let final_state = accepted.last().unwrap_or_else(|| state.get("baseline").unwrap());
        assert_eq!(state.get("area_constraint_satisfied"), Some(&J::Bool(true)));
        assert!((final_state.f64_field("area_m2").unwrap() - state.f64_field("area_target_m2").unwrap()).abs()
            <= state.f64_field("area_tolerance_m2").unwrap());
        if count == 0 {
            assert_eq!(state.get("terminal_refusals").and_then(J::as_array).unwrap().len(), 2);
        }
        let summary = J::parse(std::str::from_utf8(&retained(&db, &result, "report_json")).unwrap()).unwrap();
        let design = J::parse(std::str::from_utf8(&retained(&db, &result, "design")).unwrap()).unwrap();
        let dwr = summary.get("goal_error_assessment").unwrap();
        assert_eq!(dwr.str_field("status"), Some("estimated"));
        assert_eq!(summary.str_field("status"), Some(status));
        let snapshot = dwr.str_field("snapshot").unwrap();
        for actual in [final_state.str_field("snapshot"), summary.str_field("snapshot"), design.str_field("snapshot")] {
            assert_eq!(actual, Some(snapshot));
        }
        for (assessment_key, state_key, summary_key) in [
            ("coarse_compliance_j", "compliance_j", "final_compliance_j"),
            ("material_area_m2", "area_m2", "final_material_area_m2"),
        ] {
            let actual = dwr.f64_field(assessment_key).unwrap();
            assert_eq!(Some(actual), final_state.f64_field(state_key));
            assert_eq!(Some(actual), summary.f64_field(summary_key));
        }
        assert_eq!(result.path(&["receipt", "continuation", "goal_error_assessment"]), Some(dwr));
        assert_eq!(dwr.f64_field("coarse_level"), Some(3.0));
        assert_eq!(dwr.f64_field("enriched_level"), Some(4.0));
        assert!(dwr.f64_field("coarse_dofs").unwrap() > 0.0);
        assert!(dwr.f64_field("enriched_dofs").unwrap() > dwr.f64_field("coarse_dofs").unwrap());
        for key in ["eta_signed_j", "absolute_indicator_sum_j", "enriched_compliance_j"] {
            assert!(dwr.f64_field(key).unwrap().is_finite());
        }
        let solver = dwr.get("solver").unwrap();
        assert_eq!(solver.str_field("relative_residual_kind"), Some("recomputed-euclidean"));
        assert_eq!(solver.f64_field("solves"), Some(2.0));
        for key in ["coarse_iterations", "enriched_iterations"] {
            assert!(solver.f64_field(key).unwrap() > 0.0);
        }
        for key in ["coarse_relative_residual", "enriched_relative_residual"] {
            let residual = solver.f64_field(key).unwrap();
            assert!(residual.is_finite() && (0.0..=1e-12).contains(&residual));
        }
        let exported = document(&command("report").arg(run_id(&result)).arg(&db)
            .output().unwrap(), fs_cli::exit::SUCCESS);
        assert_eq!(exported.str_field("study_status"), Some(status));
        let html = String::from_utf8(retained(&db, &result, "report_html")).unwrap();
        assert!(html.contains("Final compliance goal-error estimate")
            && html.contains("not certified continuum-error bounds"));
        let again = document(&command("study").arg("--resume").arg(run_id(&result)).arg(&db)
            .output().unwrap(), exit);
        assert_eq!(again, result, "terminal assessment is reused without further solves or candidates");
        assert_eq!(retained(&db, &again, "report_json"), retained(&db, &result, "report_json"));
    }
}
