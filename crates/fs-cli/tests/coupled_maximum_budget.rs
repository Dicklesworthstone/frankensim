//! Run the actual .fsim import/solve/report/package path; the coupled bound
//! must be attached to the published field, not merely available as an API.
#[path = "../src/json_read.rs"]
mod json_read;

use json_read::JsonValue;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);
fn scratch() -> PathBuf {
    let base = std::env::temp_dir();
    loop {
        let path = base.join(format!("fs-cli-coupled-budget-{}-{}", std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)));
        match std::fs::create_dir(&path) {
            Ok(()) => return path,
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {},
            Err(error) => panic!("scratch directory: {error}"),
        }
    }
}
fn command(args: &[&str]) -> JsonValue {
    let out = fs_cli::run(args.iter().map(|s| s.to_string()).collect::<Vec<_>>());
    assert_eq!(out.exit_code, fs_cli::exit::SUCCESS, "{}\n{}", out.stdout, out.stderr);
    JsonValue::parse(&out.stdout).unwrap()
}
fn execute(dir: &Path) -> (String, Vec<u8>, Vec<u8>, JsonValue) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let project = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let mesh = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let ledger_path = dir.join("run.db");
    command(&["--json", "import", project.to_str().unwrap(), mesh.to_str().unwrap(),
        ledger_path.to_str().unwrap(), "--unit", "m", "--max-hole-edges", "0"]);
    let out = command(&["--json", "run", project.to_str().unwrap(),
        ledger_path.to_str().unwrap(), "--materials", pack.to_str().unwrap()]);
    assert_eq!(out.str_field("status"), Some("completed"));
    let id = out.str_field("run").unwrap().to_string();
    let report_bytes = std::fs::read(dir.join(format!("{id}.report.json"))).unwrap();
    let package_bytes = std::fs::read(dir.join(format!("{id}.fspkg"))).unwrap();
    let report = JsonValue::parse(std::str::from_utf8(&report_bytes).unwrap()).unwrap();
    let stage = report.get("stages").and_then(JsonValue::as_array).unwrap().iter()
        .find(|stage| stage.str_field("stage") == Some("conduction")).unwrap();
    let hash = fs_blake3::ContentHash::from_hex(stage.str_field("receipt_hash").unwrap()).unwrap();
    let ledger = fs_ledger::Ledger::open(ledger_path.to_str().unwrap()).unwrap();
    let bytes = ledger.get_artifact(&hash).unwrap().unwrap();
    let conduction = JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    (id, report_bytes, package_bytes, conduction)
}

#[test]
fn fsim_coupled_budget_uses_published_feedback_and_sealed_reports_replay() {
    let first_dir = scratch();
    let (id, report_bytes, package, conduction) = execute(&first_dir);
    let control = conduction.get("solver_control").unwrap();
    assert_eq!(control.str_field("schema"), Some("frankensim.cli.coupled-maximum-publication.v1"));
    assert_eq!(control.str_field("mode"), Some("physical-goal-correction"));
    assert_eq!(control.get("correction_supported"), Some(&JsonValue::Bool(true)));
    assert!(control.f64_field("primal_iterations").unwrap()
        <= control.f64_field("max_primal_iterations").unwrap());
    assert_published_cooling(&conduction);
    assert!(control.f64_field("air_paths").unwrap() >= 1.0);
    assert!(control.f64_field("ports").unwrap() >= 1.0);
    assert!(control.f64_field("response_iterations").unwrap()
        <= control.f64_field("max_response_iterations").unwrap());
    let report = JsonValue::parse(std::str::from_utf8(&report_bytes).unwrap()).unwrap();
    let terms: Vec<_> = report.get("budget_terms").and_then(JsonValue::as_array).unwrap()
        .iter().filter(|row| row.str_field("kind") == Some("solver-algebraic")).collect();
    assert_eq!(terms.len(), 1, "coupled and tolerance estimates must not both enter the budget");
    match control.f64_field("final_bound_k") {
        Some(bound) => {
            assert!(bound.is_finite() && bound >= 0.0);
            // Propagated half-widths publish as `interval` (QoiTermReceipt::interval).
            assert_eq!(terms[0].str_field("state"), Some("interval"));
            assert_eq!(terms[0].f64_field("value"), Some(bound));
            assert!(terms[0].str_field("reason").unwrap().contains("coupled"));
        }
        None => {
            assert_eq!(control.str_field("status"), Some("coupled-bound-unavailable"));
            assert_eq!(terms[0].str_field("state"), Some("no-data"));
            assert!(terms[0].f64_field("value").is_none());
        }
    }
    assert!(!std::str::from_utf8(&report_bytes).unwrap().contains("tolerance-tightening-resolve"));
    // Two independent ledgers: same field, analysis, report and package.
    let (other_id, other_report, other_package, other_conduction) = execute(&scratch());
    assert_eq!(id, other_id);
    assert_eq!(conduction, other_conduction);
    assert_eq!(report_bytes, other_report);
    assert_eq!(package, other_package);
    // Export remains a projection of the sealed bytes, not a fresh solve.
    command(&["--json", "report", &id, first_dir.join("run.db").to_str().unwrap()]);
    assert_eq!(report_bytes, std::fs::read(first_dir.join(format!("{id}.report.json"))).unwrap());
}

fn branch_rows(value: &JsonValue) -> &[JsonValue] {
    value.get("branches").and_then(JsonValue::as_array)
        .unwrap_or_else(|| std::slice::from_ref(value))
}

/// Inspect the actual retained stage, not the correction API's return value.
/// A changed schema alone must not make an unwired or stale publisher pass.
fn assert_published_cooling(conduction: &JsonValue) {
    let control = conduction.get("solver_control").unwrap();
    assert_eq!(control.get("physical_accepted"), Some(&JsonValue::Bool(true)),
        "the heatsink fixture must reach the physical publication boundary: {control:?}");
    let checks = control.f64_field("physical_checks").unwrap();
    assert!(checks >= 1.0 && checks <= control.f64_field("goal_checks").unwrap() + 1.0);
    assert!(control.f64_field("physical_rejections").unwrap() <= checks);
    assert_eq!(control.get("physical_refusal"), Some(&JsonValue::Null));
    assert_eq!(control.get("candidate_bound_k"), control.get("final_bound_k"));
    let tolerance = control.f64_field("requested_tolerance_k").unwrap();
    let meets_goal = control.f64_field("final_bound_k").is_some_and(|bound| bound <= tolerance);
    assert_eq!(control.get("goal_met"), Some(&JsonValue::Bool(meets_goal)));
    if control.f64_field("final_bound_k").is_some() {
        match control.str_field("inverse_method") {
            Some("state-contraction") => {
                assert!(control.f64_field("feedback_gain_infinity_upper").unwrap() < 1.0);
            }
            Some("port-schur-dominance") => {
                assert!(control.f64_field("schur_inverse_infinity_upper").unwrap() > 0.0);
            }
            other => panic!("a finite coupled bound needs its checked inverse route: {other:?}"),
        }
    }
    let gates = control.get("physical_gates").unwrap();
    assert!(conduction.path(&["energy", "relative_closure"]).unwrap().as_f64().unwrap()
        <= gates.f64_field("energy_relative").unwrap());

    let air = conduction.get("conjugate").unwrap();
    assert_eq!(air.str_field("publication_schema"), Some("frankensim.cli.accepted-cooling.v1"));
    let history = air.get("initial_exchange").unwrap();
    let current = branch_rows(air);
    let previous = branch_rows(history);
    assert!(!current.is_empty());
    assert_eq!(current.len(), previous.len());
    for (branch, old) in current.iter().zip(previous) {
        // Preserve the original transport inputs, but never its result fields.
        for key in ["branch", "path", "inlet_k", "flow_m3_s", "mass_flow_kg_s", "air_properties"] {
            assert!(old.get(key).is_some(), "missing original transport field {key}");
            assert_eq!(branch.get(key), old.get(key), "changed transport field {key}");
        }
        assert_eq!(branch.path(&["acceleration", "method"]).unwrap().as_str(),
            Some("fgmres-physical-goal-polish"));
        assert_eq!(branch.f64_field("iterations"), control.f64_field("primal_iterations"));
        assert_eq!(branch.get("worst_recorded_imbalance_w"), Some(&JsonValue::Null));
        let limit = branch.f64_field("balance_tolerance_w").unwrap();
        assert!(branch.f64_field("enthalpy_imbalance_w").unwrap() <= limit);
        assert!(branch.f64_field("reference_delta_k").unwrap()
            <= gates.f64_field("reference_k").unwrap());
        let rows = branch.get("segments").unwrap().as_array().unwrap();
        let old_rows = old.get("segments").unwrap().as_array().unwrap();
        assert!(!rows.is_empty());
        assert_eq!(rows.len(), old_rows.len());
        let (mut solid_heat, mut air_heat) = (0.0, 0.0);
        for (row, old_row) in rows.iter().zip(old_rows) {
            for key in ["target", "order", "card", "htc_w_m2_k", "wetted_area_m2"] {
                assert!(old_row.get(key).is_some(), "missing segment metadata {key}");
                assert_eq!(row.get(key), old_row.get(key), "changed segment metadata {key}");
            }
            assert!(row.f64_field("wall_temperature_k").unwrap().is_finite());
            assert!(row.f64_field("imbalance_w").unwrap().abs() <= limit);
            assert!((row.f64_field("reference_k").unwrap()
                - row.f64_field("marched_reference_k").unwrap()).abs()
                <= gates.f64_field("reference_k").unwrap());
            solid_heat += row.f64_field("solid_heat_rate_w").unwrap();
            air_heat += row.f64_field("air_heat_rate_w").unwrap();
        }
        assert!((solid_heat - branch.f64_field("solid_total_w").unwrap()).abs() <= limit);
        assert!((air_heat - branch.f64_field("air_total_w").unwrap()).abs() <= limit);
        assert_eq!(rows.last().unwrap().f64_field("air_out_k"), branch.f64_field("outlet_k"));
    }
}
