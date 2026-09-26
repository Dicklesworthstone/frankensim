use super::*;
use crate::study::elasticity::parse as parse_study;
use fs_topols::resolution::ResolutionStage;
use fs_topols::evaluated::DesignEvaluationStage;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-projected-stress-2d.fsim"));
const CHECK: &str = "    :mesh-check (refinement\n      :extra-levels 1\n      :absolute-compliance-tolerance-j 0.00000000000001\n      :relative-compliance-tolerance 0.0\n      :area-tolerance-m2 0.005)";
fn source() -> String {
    BASE.replace("    :stress-tolerance-pa 0.0)", &format!("    :stress-tolerance-pa 0.0\n{CHECK}\n  )"))
}
fn spec() -> ElasticitySpec { parse_study(&source()).unwrap() }
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn history(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, ConstraintEvidence) {
    let value = json(out);
    let geometry = document(&linked(ledger, &value, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &value, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &value, &geometry, &rows).unwrap();
    let history = ConstraintEvidence::read(value.path(&["continuation", "constraints"]).unwrap(),
        &report, spec.projected.as_ref().unwrap()).unwrap();
    (phi, history)
}

#[test]
fn stress_mesh_source_is_explicit_and_never_silently_ignores_independent_loads() {
    let example = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-mesh-stress-2d.fsim"));
    assert!(parse_study(example).is_ok());
    let original = parse_study(BASE).unwrap();
    assert!(stress_controls(&original).unwrap().resolution.is_none());
    assert!(!original.canonical.contains("mesh-check"));
    let declared = spec();
    assert_eq!(parse_study(&declared.canonical).unwrap().id, declared.id);
    assert_ne!(declared.id, original.id);
    for (a,b) in [(":extra-levels 1", ":extra-levels 0"),
        (":extra-levels 1", ":extra-levels 1 :extra-levels 1"),
        (":relative-compliance-tolerance 0.0", ""),
        (":area-tolerance-m2 0.005", ":area-tolerance-m2 0.0")] {
        assert!(parse_study(&source().replace(a,b)).is_err());
    }
    let family = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/marquee/bracket-multi-load-2d.fsim"));
    let bad = family.replace("    :load-family", &format!("{CHECK}\n    :load-family"));
    let error = parse_study(&bad).unwrap_err();
    assert!(error.message.contains("single-load"));
}

#[test]
fn stress_mesh_unresolved_baseline_retains_real_fine_stresses_and_never_starts_a_proposal() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive_observed(&spec, &ledger, None, &gate(), None, |stage| {
        if let Stage::Mesh(stage) = stage { assert!(matches!(stage, MeshCheckStage::Baseline(_))); }
    }).unwrap();
    assert_eq!(out.status, "mesh-unresolved");
    let (phi, retained) = history(&ledger, &out, &spec);
    assert!(retained.accepted.is_empty());
    assert_eq!(retained.policy.stress, stress_controls(&spec).unwrap().stress);
    let check = retained.mesh.as_ref().unwrap();
    assert_eq!(check.baseline.rungs.len(), 2);
    assert!(same(&check.baseline.rungs[0].evaluation, retained.current()));
    let fine = fs_topols::refinement::prolongate_level_set(&phi, &[]).unwrap().geometry;
    let expected = fs_topols::evaluate_sampled_stress(&fine, fixture(&spec),
        fs_topols::OptimizeSettings { level: 4, ..settings(&spec, spec.steps) }).unwrap();
    assert_eq!(expected, check.baseline.rungs[1].evaluation);
    let old = load(&ledger, &out.pointer).unwrap();
    let again = drive_observed(&spec, &ledger, None, &gate(), Some(&old), |_| panic!("terminal must not solve")).unwrap();
    assert_eq!(again.receipt, out.receipt);
    assert_eq!(render(out, OutputMode::Json).exit_code, exit::BUDGET);
}

#[test]
fn stress_mesh_fine_sampling_cancellation_retains_the_seed_and_retry_matches_full_execution() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let mut reached = false;
    let cancelled = drive_observed(&spec, &ledger, None, &stop, None, |stage| {
        if matches!(stage, Stage::Mesh(MeshCheckStage::Baseline(ResolutionStage::Evaluate {
            level: 4, stage: DesignEvaluationStage::StressCell(_) }))) { reached = true; stop.request(); }
    }).unwrap();
    assert!(reached);
    assert_eq!(cancelled.status, "cancelled");
    let (_, retained) = history(&ledger, &cancelled, &spec);
    assert!(retained.mesh.is_none());
    assert!(retained.accepted.is_empty());
    let old = load(&ledger, &cancelled.pointer).unwrap();
    let resumed = drive(&spec, &ledger, None, &gate(), Some(&old)).unwrap();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &whole_db, None, &gate(), None).unwrap();
    assert_eq!(resumed.status, whole.status);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations"), ("report_json", "study-report-json")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(), linked(&whole_db, &json(&whole), key, kind).unwrap());
    }
    assert_eq!(load(&ledger, &cancelled.pointer).unwrap().bytes, cancelled.receipt);
}

#[test]
fn stress_mesh_candidate_cancellation_and_resume_preserve_the_real_accepted_search() {
    let loose = source().replace(":absolute-compliance-tolerance-j 0.00000000000001", ":absolute-compliance-tolerance-j 1000000.0")
        .replace(":area-tolerance-m2 0.005", ":area-tolerance-m2 1.0");
    let spec = parse_study(&loose).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let mut reached = false;
    let cancelled = drive_observed(&spec, &ledger, None, &stop, None, |stage| {
        if matches!(stage, Stage::Mesh(MeshCheckStage::Candidate { stage: ResolutionStage::Publish, .. })) {
            reached = true; stop.request();
        }
    }).unwrap();
    assert!(reached, "must complete actual candidate mechanics and stress");
    assert_eq!(cancelled.status, "cancelled");
    let (_, seed) = history(&ledger, &cancelled, &spec);
    assert!(seed.accepted.is_empty());
    assert!(seed.mesh.is_none());
    let old = load(&ledger, &cancelled.pointer).unwrap();
    let resumed = drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &whole_db, Some(1), &gate(), None).unwrap();
    assert_eq!(resumed.status, whole.status);
    let (_, a) = history(&ledger, &resumed, &spec);
    let (_, b) = history(&whole_db, &whole, &spec);
    assert_eq!(a.json(), b.json());
    assert!(a.mesh.as_ref().unwrap().checked > 0);
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(linked(&ledger, &json(&resumed), key, kind).unwrap(), linked(&whole_db, &json(&whole), key, kind).unwrap());
    }
}

#[test]
fn stress_mesh_retained_policy_and_missing_fine_measurements_are_refused() {
    let spec = spec();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = drive(&spec, &ledger, None, &gate(), None).unwrap();
    let (_, h) = history(&ledger, &out, &spec);
    let text = h.json();
    let policy = stress_controls(&spec).unwrap();
    let changed = text.replace("\"extra_levels\":1", "\"extra_levels\":2");
    assert!(read(&document(changed.as_bytes()).unwrap(), policy, &h.baseline, &h.accepted).is_err());
    let changed = text.replace("\"level\":4", "\"level\":5");
    assert!(read(&document(changed.as_bytes()).unwrap(), policy, &h.baseline, &h.accepted).is_err());
    let changed = text.replace("\"outcome\":\"baseline-unresolved\"", "\"outcome\":\"accepted\"");
    assert!(read(&document(changed.as_bytes()).unwrap(), policy, &h.baseline, &h.accepted).is_err());
}

#[test]
fn stress_mesh_public_study_report_package_and_resume_retain_the_measured_levels() {
    use std::{fs, time::{SystemTime, UNIX_EPOCH}};
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-stress-mesh-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let source = dir.join("study.fsim"); let database = dir.join("study.db");
    fs::write(&source, spec().canonical).unwrap();
    let invoke = |args: &[&str]| crate::run(args.iter().map(|s| s.to_string()));
    let out = invoke(&["--json", "study", source.to_str().unwrap(), database.to_str().unwrap()]);
    assert_eq!(out.exit_code, exit::BUDGET, "{}", out.stderr);
    let value = JsonValue::parse(&out.stdout).unwrap();
    assert_eq!(value.str_field("status"), Some("mesh-unresolved"));
    let pointer = value.str_field("run_id").unwrap();
    let report = invoke(&["--json", "report", pointer, database.to_str().unwrap()]);
    assert_eq!(report.exit_code, exit::SUCCESS, "{}", report.stderr);
    let paths = JsonValue::parse(&report.stdout).unwrap();
    let html = fs::read_to_string(paths.str_field("report_html").unwrap()).unwrap();
    assert!(html.contains("Mesh-checked sampled stress"));
    let summary = document(&fs::read(paths.str_field("report_json").unwrap()).unwrap()).unwrap();
    assert_eq!(summary.path(&["constraints", "mesh_resolution", "last_check", "baseline", "rungs"])
        .and_then(JsonValue::as_array).unwrap().len(), 2);
    let package = invoke(&["--json", "package", pointer, database.to_str().unwrap()]);
    assert_eq!(package.exit_code, exit::SUCCESS, "{}", package.stderr);
    let paths = JsonValue::parse(&package.stdout).unwrap();
    let package = EvidencePackage::from_json(&fs::read_to_string(paths.str_field("package").unwrap()).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
    let repeated = invoke(&["--json", "study", "--resume", pointer, database.to_str().unwrap()]);
    assert_eq!(repeated.exit_code, out.exit_code);
    assert_eq!(repeated.stdout, out.stdout);
}
