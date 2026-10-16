//! Native multilevel declarations must reach real equilibrium and recovery.
use super::*;
use fs_topopt::pipeline::LoadCase;

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-multilevel.fsim"));
const LEGACY: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-adaptive.fsim"));
const STRESS: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-stress.fsim"));
const POLICY: &str = "(linear-solver :type multilevel :coarsest-level 0 :max-transfer-terms 500000 :max-matrix-entries 500000 :max-galerkin-products 100000000 :max-diagonal-contributions 100000000)";

fn append(source: &str, section: &str) -> String {
    format!("{}\n{section})", source.trim_end().strip_suffix(')').unwrap())
}
fn funded(source: &str) -> String {
    source.replace(":memory-bytes 536870912", ":memory-bytes 1073741824")
}
fn small_source() -> String { FIXTURE.replace(":initial-level 2", ":initial-level 1") }
fn near(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).map(|v| v.abs()).fold(1e-20_f64, f64::max);
    for (&a, &b) in a.iter().zip(b) {
        assert!((a - b).abs() <= 2e-7 * scale, "{a:e} != {b:e}");
    }
}

#[test]
fn solver_policy_is_explicit_canonical_bounded_and_never_silently_ignored() {
    let spec = spec::parse(FIXTURE).unwrap();
    assert!(matches!(spec.solver, Policy::Multilevel { coarsest_level: 0, .. }));
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    let old = spec::parse(LEGACY).unwrap();
    assert!(matches!(old.solver, Policy::Legacy));
    let enough_memory = funded(LEGACY);
    assert_ne!(spec::parse(&append(&enough_memory, POLICY)).unwrap().id,
        spec::parse(&enough_memory).unwrap().id);
    for bad in [
        FIXTURE.replace(":coarsest-level 0", ":coarsest-level 2"),
        FIXTURE.replace(":type multilevel", ":type unknown"),
        FIXTURE.replace(":max-transfer-terms 500000", ":max-transfer-terms 0"),
        FIXTURE.replace(":max-matrix-entries 500000", ":max-matrix-entries 8000001"),
        FIXTURE.replace(":max-galerkin-products 100000000", ":max-galerkin-products 1000000001"),
        FIXTURE.replace(":max-diagonal-contributions 100000000", ":max-diagonal-contributions 0"),
        FIXTURE.replace(":memory-bytes 1073741824", ":memory-bytes 536870912"),
        append(FIXTURE, POLICY),
        append(FIXTURE, "(design-regions :solid () :void ())"),
    ] { assert!(spec::parse(&bad).is_err(), "must refuse: {bad}"); }
    let with_regions = append(&append(&enough_memory, "(design-regions :solid () :void ())"), POLICY);
    assert!(spec::parse(&with_regions).is_ok());
}

fn problem(spec: &Spec) -> (CutDensityStudy3<AdaptiveSolveSpace3>, Vec<f64>, Vec<f64>, QuadratureWork3) {
    let bounds = HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap();
    let domain = PhysicalDomain { bounds: spec.bounds, height: spec.height, curvature: spec.curvature };
    let material = IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let mut build = |tree: &Octree3, checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
        if checkpoint().is_break() { return Err(EvaluationStop::Cancelled.into()); }
        loading::build_operator(spec, bounds, tree, &domain, &material, &mut q)
            .map_err(GoalRefinementError3::from)
    };
    let corrections = correction_spaces(spec, &mut build, &mut || ControlFlow::Continue(())).unwrap();
    let refs: Vec<_> = corrections.iter().collect();
    let tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves).unwrap();
    let operator = build(&tree, &mut || ControlFlow::Continue(())).unwrap();
    let y = operator.body_load(&|_| spec.loads[0].0, || ControlFlow::Continue(())).unwrap();
    let z = operator.body_load(&|_| spec.loads[1].0, || ControlFlow::Continue(())).unwrap();
    let space = spec.solver.wrap(operator, &refs, || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(space, spec.radius, spec.schedule[0]), y, z, q.work())
}

#[test]
fn current_density_multilevel_matches_jacobi_and_prepares_once_for_independent_loads() {
    let spec = spec::parse(FIXTURE).unwrap();
    let (mut study, y, z, geometry) = problem(&spec);
    let levels = study.operator().level_sizes();
    assert_eq!(levels.len(), 3);
    assert!(levels.windows(2).all(|p| p[0] > p[1]));
    assert!(study.operator().transfer_entries() > 0);
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.38 + 0.06 * (i % 5) as f64).collect();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let actual = study.evaluate(&raw, &[
        LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 },
    ], &mut control).unwrap();
    let work = control.work();
    assert!(work.preconditioner_galerkin_products > 0);
    assert_eq!(work.preconditioner_operator_applications, 0);
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&raw, &[LoadCase { force: &z, weight: 1.0 }], &mut control).unwrap();
    assert_eq!(control.work().preconditioner_galerkin_products, work.preconditioner_galerkin_products,
        "the second independent RHS must reuse the density's preparation");
    let mut legacy = spec.clone();
    legacy.solver = Policy::Legacy;
    let (mut reference, y, z, old_geometry) = problem(&legacy);
    assert!(geometry.boxes > old_geometry.boxes, "coarse geometry must be charged");
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let expected = reference.evaluate(&raw, &[
        LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 },
    ], &mut control).unwrap();
    assert_eq!(actual.projected_rho, expected.projected_rho);
    assert_eq!(actual.volume_fraction, expected.volume_fraction);
    near(&[actual.objective.compliance], &[expected.objective.compliance]);
    near(&actual.objective.gradient, &expected.objective.gradient);
    for (a, b) in actual.objective.displacements.iter().zip(&expected.objective.displacements) { near(a, b); }
    assert_ne!(actual.objective.displacements[0], actual.objective.displacements[1]);
}

fn bytes(ledger: &Ledger, outcome: &Outcome, key: &str, kind: &str) -> Vec<u8> {
    let receipt = JsonValue::parse(&outcome.receipt).unwrap();
    linked(ledger, &receipt, key, kind).unwrap()
}
fn json(ledger: &Ledger, outcome: &Outcome, key: &str, kind: &str) -> JsonValue {
    JsonValue::parse(std::str::from_utf8(&bytes(ledger, outcome, key, kind)).unwrap()).unwrap()
}
fn setup_work(outcome: &Outcome) -> f64 {
    JsonValue::parse(&outcome.receipt).unwrap()
        .path(&["work", "preconditioner_galerkin_products"]).and_then(JsonValue::as_f64).unwrap()
}

#[test]
fn native_multilevel_adaptivity_and_ledger_resume_retain_the_same_accepted_fields() {
    let spec = spec::parse(&small_source()).unwrap();
    let gate = CancelGate::new_clock_free();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = checkpoint::drive(&spec, &whole_db, None, &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(whole.status, "completed");
    assert!(setup_work(&whole) > 0.0);
    let trace = json(&whole_db, &whole, "iterations", "study-iterations");
    let stages = trace.get("stages").and_then(JsonValue::as_array).unwrap();
    assert_eq!(stages.len(), 2);
    for stage in stages {
        assert_eq!(stage.get("gradient_check_passed"), Some(&JsonValue::Bool(true)));
        assert!(stage.get("history").and_then(JsonValue::as_array).unwrap().len() > 1);
    }
    let report = json(&whole_db, &whole, "report_json", "study-report-json");
    assert_eq!(report.get("refinements").and_then(JsonValue::as_array).unwrap()[0].get("installed"),
        Some(&JsonValue::Bool(true)));
    let resumed_db = Ledger::open(":memory:").unwrap();
    let first = checkpoint::drive(&spec, &resumed_db, Some(1), &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(first.status, "budget-exhausted");
    let old = load(&resumed_db, &first.pointer).unwrap();
    let resumed = checkpoint::drive(&spec, &resumed_db, None, &gate, Some(&old), |_, _| Ok(())).unwrap();
    assert_eq!(resumed.status, "completed");
    assert!(setup_work(&resumed) > setup_work(&whole), "prefix replay spends setup again");
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(bytes(&whole_db, &whole, key, kind), bytes(&resumed_db, &resumed, key, kind));
    }
}

#[test]
fn native_multilevel_exhaustion_and_precancellation_do_not_publish_unfunded_fields() {
    let source = small_source().replace(":max-galerkin-products 100000000", ":max-galerkin-products 1");
    let spec = spec::parse(&source).unwrap();
    let run = compute(&spec, &CancelGate::new_clock_free()).unwrap();
    assert_eq!(run.status, "budget-exhausted");
    assert!(run.report.continuation.last.is_none());
    assert_eq!(run.spent.linear.preconditioner_galerkin_products, 1);
    let source = small_source().replace(":max-transfer-terms 500000", ":max-transfer-terms 1");
    let spec = spec::parse(&source).unwrap();
    assert!(compute(&spec, &CancelGate::new_clock_free()).is_err());
    let spec = spec::parse(&small_source()).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    match compute(&spec, &gate) {
        Err(error) => assert_eq!(error.exit, exit::CANCELLED),
        Ok(_) => panic!("precancelled study cannot start correction geometry"),
    }
}

#[test]
fn stress_solver_policy_survives_direct_optimizer_state_restoration() {
    let source = append(&funded(STRESS), POLICY);
    let spec = spec::parse(&source).unwrap();
    let gate = CancelGate::new_clock_free();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = stress::drive(&spec, &whole_db, Some(2), &gate).unwrap();
    assert!(setup_work(&whole) > 0.0);
    let summary = json(&whole_db, &whole, "report_json", "study-report-json");
    assert_eq!(summary.path(&["gradient_check", "passed"]), Some(&JsonValue::Bool(true)));
    assert_eq!(summary.get("iterations_completed").and_then(JsonValue::as_f64), Some(2.0));
    let resumed_db = Ledger::open(":memory:").unwrap();
    let first = stress::drive(&spec, &resumed_db, Some(1), &gate).unwrap();
    let old = load(&resumed_db, &first.pointer).unwrap();
    let resumed = stress::resume(&resumed_db, &old, Some(1), &gate).unwrap();
    assert!(setup_work(&resumed) > setup_work(&whole), "endpoint verification uses current-density factors");
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(bytes(&whole_db, &whole, key, kind), bytes(&resumed_db, &resumed, key, kind));
    }
}
