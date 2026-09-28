//! G0/G3/G4/G5 consumer checks: real native studies and real finer-grid solves.
use super::*;
use crate::study::elasticity::parse as parse_study;
use fs_topols::robust_resolution::{MultiLoadMeshStage, MultiLoadResolutionStage};

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-multi-load-mesh-2d.fsim"));
fn source(strict: bool) -> String {
    EXAMPLE.replace(":absolute-compliance-tolerance-j 0.01",
        if strict { ":absolute-compliance-tolerance-j 1e-30" }
        else { ":absolute-compliance-tolerance-j 1000000.0" })
        .replace(":relative-compliance-tolerance 0.05", ":relative-compliance-tolerance 0.0")
        .replace(":area-tolerance-m2 0.005", ":area-tolerance-m2 1.0")
}
fn gate() -> CancelGate { CancelGate::new_clock_free() }
fn json(out: &Outcome) -> JsonValue { document(out.receipt.as_bytes()).unwrap() }
fn history(ledger: &Ledger, out: &Outcome, spec: &ElasticitySpec) -> (GridSdf, ConstraintEvidence) {
    let value = json(out);
    let geometry = document(&linked(ledger, &value, "design", "study-design").unwrap()).unwrap();
    let rows = document(&linked(ledger, &value, "iterations", "study-iterations").unwrap()).unwrap();
    let (phi, report) = decode(spec, &value, &geometry, &rows).unwrap();
    let retained = ConstraintEvidence::read(value.path(&["continuation", "constraints"]).unwrap(),
        &report, spec.projected.as_ref().unwrap()).unwrap();
    (phi, retained)
}

#[test]
fn native_mesh_family_is_explicit_and_restoration_is_not_silently_ignored() {
    let spec = parse_study(EXAMPLE).unwrap();
    assert_eq!(parse_study(&spec.canonical).unwrap().id, spec.id);
    let policy = stress_controls(&spec).unwrap();
    assert_eq!(policy.resolution.unwrap().extra_levels, 1);
    assert_eq!(policy.family.as_ref().unwrap().additional[0].weight(), 0.0);
    for (from, to) in [(":extra-levels 1", ":extra-levels 3"),
        (":relative-compliance-tolerance 0.05", ":relative-compliance-tolerance 1.0"),
        (":max-recovery-solves 64", ":max-recovery-solves 64 :stress-restoration-reduction 0.01")] {
        assert!(parse_study(&EXAMPLE.replace(from, to)).is_err(), "{to}");
    }
    assert!(parse_study(&EXAMPLE.replace("weighted-sum", "worst-weighted-case")).is_ok());
}

#[test]
fn unresolved_native_baseline_retains_real_independent_fine_stresses_without_search() {
    let spec = parse_study(&source(true)).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let out = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), None,
        |_| panic!("unresolved baseline must not start optimization"),
        |_| panic!("unresolved baseline must not start final DWR"), |stage| {
            assert!(matches!(stage, MultiLoadMeshStage::Baseline(_)));
        }).unwrap();
    assert_eq!(out.status, "mesh-unresolved");
    let (phi, retained) = history(&ledger, &out, &spec);
    assert!(retained.accepted.is_empty());
    let h = retained.family.as_ref().unwrap();
    assert_eq!(h.solves, 4, "two coarse and two actual finer solves");
    let check = h.mesh.as_ref().unwrap();
    assert_eq!(check.baseline.rungs[0].cases, h.baseline);
    let fine = fs_topols::refinement::prolongate_level_set(&phi, &[]).unwrap().geometry;
    let replay = fs_topols::evaluate_robust_sampled_stress(&fine, &h.cases,
        OptimizeSettings { level: 4, ..settings(&spec, spec.steps) }, RobustAggregate::WeightedSum).unwrap();
    assert_eq!(check.baseline.rungs[1].cases, case_states(&replay));
    for rung in &check.baseline.rungs {
        assert!(rung.cases[0].compliance > 0.0 && rung.cases[0].sampled_max_von_mises > 0.0);
        assert!((rung.cases[1].compliance / rung.cases[0].compliance - 16.0).abs() < 1e-7);
        assert!((rung.cases[1].sampled_max_von_mises / rung.cases[0].sampled_max_von_mises - 4.0).abs() < 1e-7);
    }
    // Mutate actual measured evidence to isolate the all-case gate, not to
    // replace a solver. A zero objective weight must not hide fine overstress.
    let mut bad = check.baseline.clone();
    let policy = stress_controls(&spec).unwrap();
    let loose = ResolutionPolicy { absolute_compliance_tolerance: 1e8, ..policy.resolution.unwrap() };
    bad.rungs[1].cases[1].sampled_max_von_mises = 2.0 * policy.stress.admitted_max();
    let reason = bad.refusal(loose, policy.area.target, policy.stress, &h.cases).unwrap().unwrap();
    assert!(reason.contains("case 1") && reason.contains("level 4"));
    let old = load(&ledger, &out.pointer).unwrap();
    let again = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), Some(&old),
        |_| panic!("terminal optimization is inert"), |_| panic!("no assessment"),
        |_| panic!("terminal mesh result is reused")).unwrap();
    assert_eq!(again.receipt, out.receipt);
}

#[test]
fn native_mesh_work_reserves_whole_families_and_exhausted_resume_does_no_recovery() {
    for (cap, expected) in [(2, 2), (3, 2), (5, 4)] {
        let spec = parse_study(&source(false).replace(":max-solves 512", &format!(":max-solves {cap}"))).unwrap();
        let ledger = Ledger::open(":memory:").unwrap();
        let out = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), None,
            |_| panic!("no candidate can fit"), |_| panic!("no final assessment"), |_| {
                assert_eq!(cap, 5, "an unfunded baseline must not start");
            }).unwrap();
        assert_eq!(out.status, "budget-exhausted");
        let (_, retained) = history(&ledger, &out, &spec);
        let h = retained.family.as_ref().unwrap();
        assert_eq!(h.solves, expected);
        assert!(h.accepted.is_empty());
        assert_eq!(h.recovery_solves, 0);
        let old = load(&ledger, &out.pointer).unwrap();
        let again = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), Some(&old),
            |_| panic!("exhausted optimizer is inert"), |_| panic!("no DWR"),
            |_| panic!("exhausted mesh study is inert")).unwrap();
        assert_eq!(again.receipt, out.receipt);
    }
}

#[test]
fn interrupted_fine_work_keeps_geometry_and_retry_reproduces_a_genuine_accepted_update() {
    let spec = parse_study(&source(false)).unwrap();
    let full_db = Ledger::open(":memory:").unwrap();
    let full = driver::drive(&spec, &full_db, Some(1), &gate(), None).unwrap();
    let (_, full_history) = history(&full_db, &full, &spec);
    assert_eq!(full_history.accepted.len(), 1, "exercise actual fine-grid accepted design, not a vacuous replay");
    let expected = full_history.family.as_ref().unwrap();
    assert_eq!(expected.mesh.as_ref().unwrap().outcome, "accepted");
    for stop_at_candidate in [false, true] {
        let ledger = Ledger::open(":memory:").unwrap();
        let stop = gate();
        let mut reached = false;
        let out = driver::drive_mesh_observed(&spec, &ledger, Some(1), &stop, None, |_| {}, |_| {}, |stage| {
            let hit = if stop_at_candidate {
                matches!(stage, MultiLoadMeshStage::Candidate { stage: MultiLoadResolutionStage::Publish, .. })
            } else {
                matches!(stage, MultiLoadMeshStage::Baseline(MultiLoadResolutionStage::CaseIterations { case: 1, iterations, .. }) if iterations > 0)
            };
            if hit { reached = true; stop.request(); }
        }).unwrap();
        assert!(reached);
        assert_eq!(out.status, "cancelled");
        let (phi, retained) = history(&ledger, &out, &spec);
        assert_eq!(snapshot(&phi), retained.baseline.snapshot);
        assert!(retained.accepted.is_empty());
        let h = retained.family.as_ref().unwrap();
        assert!(h.mesh.is_none(), "partial grid/load families are not published");
        let spent = h.solves - h.cases.len();
        assert!(spent >= 2);
        let old = load(&ledger, &out.pointer).unwrap();
        let retry = driver::drive(&spec, &ledger, Some(1), &gate(), Some(&old)).unwrap();
        let (_, retried) = history(&ledger, &retry, &spec);
        let actual = retried.family.as_ref().unwrap();
        assert_eq!(actual.solves, expected.solves + spent);
        assert_eq!(actual.recovery_solves, 4);
        assert_eq!(actual.mesh.as_ref().unwrap().json(), expected.mesh.as_ref().unwrap().json());
        for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
            assert_eq!(linked(&ledger, &json(&retry), key, kind).unwrap(),
                linked(&full_db, &json(&full), key, kind).unwrap());
        }
    }
}

#[test]
fn final_weighted_dwr_follows_mesh_acceptance_and_retries_without_optimizer_recovery() {
    let text = source(false).replace(":steps 2", ":steps 1")
        .replace(":max-iterations 2", ":max-iterations 1")
        .replace(":max-recovery-solves 64", ":max-recovery-solves 0");
    let text = format!("{}\n  (assessment :type elasticity-dwr :max-solves-per-attempt 4)\n)\n",
        text.trim_end().strip_suffix(')').unwrap());
    let spec = parse_study(&text).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let stop = gate();
    let mut reached = false;
    let out = driver::drive_mesh_observed(&spec, &ledger, None, &stop, None, |_| {}, |stage| {
        if matches!(stage, fs_topols::WeightedComplianceDwrStage::BeforeCase { case: 1 }) {
            reached = true; stop.request();
        }
    }, |_| {}).unwrap();
    assert!(reached);
    assert_eq!(out.status, "cancelled");
    let old = load(&ledger, &out.pointer).unwrap();
    let retry = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), Some(&old),
        |_| panic!("no optimizer recovery"), |_| {}, |_| panic!("mesh acceptance already retained")).unwrap();
    assert_eq!(retry.status, "completed");
    let (_, retained) = history(&ledger, &retry, &spec);
    assert_eq!(retained.family.as_ref().unwrap().recovery_solves, 0);
    let result = json(&retry);
    assert_eq!(result.path(&["continuation", "goal_error_assessment", "status"])
        .and_then(JsonValue::as_str), Some("estimated"));
    let old = load(&ledger, &retry.pointer).unwrap();
    let again = driver::drive_mesh_observed(&spec, &ledger, None, &gate(), Some(&old),
        |_| panic!("completed"), |_| panic!("completed"), |_| panic!("completed")).unwrap();
    assert_eq!(again.receipt, retry.receipt);
}

#[test]
fn public_study_export_and_source_free_resume_preserve_all_load_grid_measurements() {
    use std::{fs, time::{SystemTime, UNIX_EPOCH}};
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-multi-load-mesh-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let source_path = dir.join("study.fsim");
    let database = dir.join("study.db");
    fs::write(&source_path, source(true)).unwrap();
    let invoke = |args: &[&str]| crate::run(args.iter().map(|s| s.to_string()));
    let out = invoke(&["--json", "study", source_path.to_str().unwrap(), database.to_str().unwrap()]);
    assert_eq!(out.exit_code, exit::BUDGET, "{}", out.stderr);
    let result = JsonValue::parse(&out.stdout).unwrap();
    assert_eq!(result.str_field("status"), Some("mesh-unresolved"));
    let pointer = result.str_field("run_id").unwrap();
    fs::rename(&source_path, dir.join("relocated.fsim")).unwrap();
    let exported = invoke(&["--json", "report", pointer, database.to_str().unwrap()]);
    assert_eq!(exported.exit_code, exit::SUCCESS, "{}", exported.stderr);
    let paths = JsonValue::parse(&exported.stdout).unwrap();
    let report = document(&fs::read(paths.str_field("report_json").unwrap()).unwrap()).unwrap();
    let rungs = report.path(&["constraints", "load_family", "mesh_resolution", "last_check", "baseline", "rungs"])
        .and_then(JsonValue::as_array).unwrap();
    assert_eq!(rungs.len(), 2);
    assert!(rungs.iter().all(|row| row.get("cases").and_then(JsonValue::as_array).unwrap().len() == 2));
    assert!(fs::read_to_string(paths.str_field("report_html").unwrap()).unwrap().contains("Mesh-checked independent loads"));
    let exported = invoke(&["--json", "package", pointer, database.to_str().unwrap()]);
    assert_eq!(exported.exit_code, exit::SUCCESS, "{}", exported.stderr);
    let paths = JsonValue::parse(&exported.stdout).unwrap();
    let package = EvidencePackage::from_json(&fs::read_to_string(paths.str_field("package").unwrap()).unwrap()).unwrap();
    assert!(fs_checker::check(&package).passed());
    let again = invoke(&["--json", "study", "--resume", pointer, database.to_str().unwrap()]);
    assert_eq!(again.exit_code, out.exit_code);
    assert_eq!(again.stdout, out.stdout);
}
