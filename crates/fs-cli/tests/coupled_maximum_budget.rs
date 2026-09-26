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
    assert_eq!(control.str_field("schema"), Some("frankensim.cli.coupled-maximum-evidence.v1"));
    assert_eq!(control.str_field("mode"), Some("assessment-only"));
    assert_eq!(control.f64_field("primal_iterations"), Some(0.0));
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
            assert_eq!(terms[0].str_field("state"), Some("measured"));
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
