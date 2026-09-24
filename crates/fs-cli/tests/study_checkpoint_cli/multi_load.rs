//! Native executable and on-disk ledger coverage for independent load families.
use super::*;
use fs_topols::{DesignBoxEdge, GridSdf, OptimizeSettings, RobustAggregate, RobustLoadCase};

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-multi-load-2d.fsim"));

fn source(dir: &Path, text: &str) -> PathBuf {
    let path = dir.join("multi-load.fsim");
    fs::write(&path, text).unwrap();
    path
}
fn constraints(value: &J) -> &J {
    value.path(&["receipt", "continuation", "constraints"]).unwrap()
}
fn family(value: &J) -> &J { constraints(value).get("load_family").unwrap() }
fn field(database: &Path, result: &J) -> GridSdf {
    let bytes = retained(database, result, "design");
    let design = J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    let n = design.get("n").and_then(J::number_raw).unwrap().parse::<usize>().unwrap();
    let nodes = design.get("phi_bits").and_then(J::as_array).unwrap();
    let mut phi = GridSdf::from_fn(n, &|_, _| 0.0);
    assert_eq!(nodes.len(), phi.nodes().len());
    for (slot, value) in phi.nodes_mut().iter_mut().zip(nodes) {
        *slot = f64::from_bits(u64::from_str_radix(value.as_str().unwrap(), 16).unwrap());
    }
    phi
}

#[test]
fn noncollinear_native_loads_match_independent_pde_replay_and_retained_exports() {
    let dir = scratch("independent-loads");
    // Exactly the baseline family may run. This exercises the real numerical
    // producer without assuming that an arbitrary second load admits a step.
    let input = source(&dir, &EXAMPLE.replace(":max-solves 512", ":max-solves 2"));
    let db = dir.join("study.db");
    let out = document(&command("study").arg(&input).arg(&db).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(out.str_field("status"), Some("budget-exhausted"));
    assert_eq!(updates(&out), 0.0);
    let measured = family(&out).get("baseline_cases").and_then(J::as_array).unwrap();
    assert_eq!(measured.len(), 2);
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 1.0).unwrap(),
    ];
    let oracle = fs_topols::evaluate_robust_sampled_stress(&field(&db, &out), &cases,
        OptimizeSettings { level: 3, youngs: 2.0, poisson: 0.3, ..OptimizeSettings::default() },
        RobustAggregate::WeightedSum).unwrap();
    for (i, state) in measured.iter().enumerate() {
        for (name, expected) in [("compliance_j", oracle.case_compliances[i]),
            ("sampled_von_mises_pa", oracle.case_sampled_max_von_mises[i]),
            ("area_m2", oracle.volume), ("max_x", oracle.case_max_locations[i][0]),
            ("max_y", oracle.case_max_locations[i][1])]
        { assert_eq!(state.f64_field(name).unwrap().to_bits(), expected.to_bits(), "case {i}: {name}"); }
        assert_eq!(state.get("sample_count").and_then(J::number_raw).unwrap()
            .parse::<usize>().unwrap(), oracle.case_sample_counts[i]);
        assert_eq!(state.str_field("snapshot"), Some(format!("{:#018x}", oracle.snapshot).as_str()));
    }
    let baseline = constraints(&out).get("baseline").unwrap();
    assert_eq!(baseline.f64_field("compliance_j").unwrap().to_bits(), oracle.objective.to_bits());
    assert_eq!(baseline.f64_field("sampled_von_mises_pa").unwrap().to_bits(), oracle.worst_sampled_von_mises.to_bits());
    assert!(oracle.case_compliances.iter().all(|value| *value > 0.0));
    assert_eq!(family(&out).f64_field("solves_started"), Some(2.0));
    for verb in ["report", "package"] {
        let exported = document(&command(verb).arg(run_id(&out)).arg(&db).output().unwrap(), fs_cli::exit::SUCCESS);
        assert_eq!(exported.str_field("study_status"), Some("budget-exhausted"));
    }
    let bytes = fs::read(dir.join(format!("{}.json", run_id(&out)))).unwrap();
    assert_eq!(bytes, retained(&db, &out, "report_json"));
    let report = J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(report.get("constraints"), Some(constraints(&out)));
    let html = fs::read_to_string(dir.join(format!("{}.html", run_id(&out)))).unwrap();
    assert!(html.contains("independent right-edge operating conditions"));
    assert!(html.contains("weights are not probabilities"));
    let repeated = document(&command("study").arg("--resume").arg(run_id(&out)).arg(&db)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(repeated, out, "exhausted study work must not restart on resume");
}

#[test]
fn native_multi_load_split_run_keeps_exact_accepted_cases_geometry_and_study_work() {
    let dir = scratch("multi-load-chunks");
    // Opposite operating conditions remain distinct positive-demand solves.
    // This deliberately informative fixture must accept a real first update.
    let input = source(&dir, &EXAMPLE.replace(":traction-pa (0.5 0.0)", ":traction-pa (0.0 1.0)"));
    let db = dir.join("chunks.db");
    let whole_db = dir.join("whole.db");
    let first = document(&command("study").arg(&input).arg(&db).args(["--budget", "1"])
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(first.str_field("status"), Some("budget-exhausted"));
    assert_eq!(updates(&first), 1.0);
    let before = retained(&db, &first, "design");
    let full_output = command("study").arg(&input).arg(&whole_db).output().unwrap();
    let code = full_output.status.code().unwrap();
    assert!(code == i32::from(fs_cli::exit::SUCCESS) || code == i32::from(fs_cli::exit::REFUSED),
        "stdout={} stderr={}", String::from_utf8_lossy(&full_output.stdout), String::from_utf8_lossy(&full_output.stderr));
    let full = document(&full_output, u8::try_from(code).unwrap());
    assert!(matches!(full.str_field("status"), Some("completed" | "no-feasible-descent")));
    fs::rename(&input, dir.join("original-input-retained.fsim")).unwrap();
    let resumed = document(&command("study").arg("--resume").arg(run_id(&first)).arg(&db)
        .output().unwrap(), u8::try_from(code).unwrap());
    assert_eq!(resumed.str_field("status"), full.str_field("status"));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&db, &resumed, key), retained(&whole_db, &full, key));
    }
    for key in ["baseline", "accepted", "candidate_counts", "terminal_refusals", "relative_reduction"] {
        assert_eq!(constraints(&resumed).get(key), constraints(&full).get(key), "{key}");
    }
    for key in ["cases", "baseline_cases", "accepted_cases", "checkpoint_hex", "solves_started"] {
        assert_eq!(family(&resumed).get(key), family(&full).get(key), "{key}");
    }
    assert_eq!(family(&resumed).f64_field("recovery_solves_used"), Some(4.0));
    assert_eq!(family(&full).f64_field("recovery_solves_used"), Some(0.0));
    assert_eq!(retained(&db, &first, "design"), before);
    assert_eq!(resumed.path(&["receipt", "continuation", "legacy_prefix_updates_replayed"])
        .and_then(J::as_f64), Some(0.0));
}

#[test]
fn invalid_native_load_families_refuse_before_creating_persistent_state() {
    for (index, source_text) in [
        EXAMPLE.replace(":max-solves 512", ":max-solves 1"),
        EXAMPLE.replace("      :aggregate weighted-sum\n", ""),
        EXAMPLE.replace(":weight 1.0", ":weight 1.0 :weight 2.0"),
    ].into_iter().enumerate() {
        let dir = scratch(&format!("load-family-refusal-{index}"));
        let input = source(&dir, &source_text);
        let db = dir.join("must-not-exist.db");
        refusal(&command("study").arg(&input).arg(&db).output().unwrap(), fs_cli::exit::REFUSED);
        assert!(!db.exists());
    }
}

const RESTORATION: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-stress-restoration-2d.fsim"));

#[test]
fn primary_only_restoration_exports_an_infeasible_result_without_inventing_progress() {
    let dir = scratch("restoration-infeasible-export");
    let input = source(&dir, &RESTORATION.replace(":max-solves 256", ":max-solves 1")
        .replace(":sampled-stress-limit-pa 1.0", ":sampled-stress-limit-pa 0.000000000001"));
    let db = dir.join("study.db");
    let out = document(&command("study").arg(&input).arg(&db).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(out.str_field("status"), Some("budget-exhausted"));
    let info = constraints(&out);
    assert_eq!(info.str_field("baseline_scope"), Some("area_feasible_study_start"));
    assert_eq!(info.get("relative_reduction"), Some(&J::Null));
    let phase = info.get("stress_restoration").unwrap();
    assert_eq!(phase.get("stress_feasible"), Some(&J::Bool(false)));
    assert_eq!(phase.get("feasible_baseline"), Some(&J::Null));
    assert_eq!(phase.f64_field("restoration_updates"), Some(0.0));
    assert!(phase.f64_field("remaining_excess_pa").unwrap() > 0.0);
    assert_eq!(family(&out).get("cases").and_then(J::as_array).unwrap().len(), 1);
    assert_eq!(family(&out).f64_field("solves_started"), Some(1.0));
    let exported = document(&command("report").arg(run_id(&out)).arg(&db).output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(exported.str_field("study_status"), Some("budget-exhausted"));
    let report_bytes = fs::read(dir.join(format!("{}.json", run_id(&out)))).unwrap();
    assert_eq!(report_bytes, retained(&db, &out, "report_json"));
    let report = J::parse(std::str::from_utf8(&report_bytes).unwrap()).unwrap();
    assert_eq!(report.get("constraints"), Some(info));
    let html = fs::read_to_string(dir.join(format!("{}.html", run_id(&out)))).unwrap();
    assert!(html.contains("Stress remains infeasible"));
    assert!(html.contains("Area-feasible input compliance"));
    fs::rename(&input, dir.join("original-restoration-input.fsim")).unwrap();
    let again = document(&command("study").arg("--resume").arg(run_id(&out)).arg(&db)
        .output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(again, out, "an infeasible exhausted terminal cannot restart its allowance");
}

#[test]
fn restorative_native_commands_preserve_repair_history_and_checkpoint_across_chunks() {
    let dir = scratch("restoration-chunks");
    // Deliberately impossible engineering limit: this exercises honest partial
    // repair and stopping, not a claimed feasible optimum or successful example.
    let text = RESTORATION.replace(":max-iterations 8", ":max-iterations 2")
        .replace(":steps 8", ":steps 2")
        .replace(":stress-restoration-reduction 0.01", ":stress-restoration-reduction 0.0")
        .replace(":sampled-stress-limit-pa 1.0", ":sampled-stress-limit-pa 0.000000000001");
    let input = source(&dir, &text);
    let chunks = dir.join("chunks.db");
    let whole = dir.join("whole.db");
    let read = |output: Output| {
        let code = output.status.code().unwrap();
        assert!(code == i32::from(fs_cli::exit::BUDGET) || code == i32::from(fs_cli::exit::REFUSED),
            "stdout={} stderr={}", String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
        let result = document(&output, u8::try_from(code).unwrap());
        assert!(matches!(result.str_field("status"), Some("budget-exhausted" | "no-feasible-descent")));
        assert_eq!(constraints(&result).path(&["stress_restoration", "stress_feasible"]), Some(&J::Bool(false)));
        assert_eq!(constraints(&result).get("relative_reduction"), Some(&J::Null));
        result
    };
    let first = read(command("study").arg(&input).arg(&chunks).args(["--budget", "1"]).output().unwrap());
    assert!(family(&first).f64_field("solves_started").unwrap() > 1.0, "must attempt actual restoration physics");
    let previous = retained(&chunks, &first, "design");
    let full = read(command("study").arg(&input).arg(&whole).output().unwrap());
    fs::rename(&input, dir.join("original-input-retained.fsim")).unwrap();
    let resumed = read(command("study").arg("--resume").arg(run_id(&first)).arg(&chunks).output().unwrap());
    assert_eq!(resumed.str_field("status"), full.str_field("status"));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&chunks, &resumed, key), retained(&whole, &full, key));
    }
    for key in ["baseline", "accepted", "candidate_counts", "terminal_refusals", "stress_restoration"] {
        assert_eq!(constraints(&resumed).get(key), constraints(&full).get(key), "{key}");
    }
    for key in ["checkpoint_hex", "solves_started", "accepted_cases"] {
        assert_eq!(family(&resumed).get(key), family(&full).get(key), "{key}");
    }
    let mut previous_stress = constraints(&full).get("baseline").unwrap().f64_field("sampled_von_mises_pa").unwrap();
    for state in constraints(&full).get("accepted").and_then(J::as_array).unwrap() {
        let stress = state.f64_field("sampled_von_mises_pa").unwrap();
        assert!(stress < previous_stress);
        previous_stress = stress;
    }
    let expected_recovery = if first.str_field("status") == Some("no-feasible-descent") { 0.0 } else { 2.0 };
    assert_eq!(family(&resumed).f64_field("recovery_solves_used"), Some(expected_recovery));
    assert_eq!(retained(&chunks, &first, "design"), previous);
}
