//! Actual native fixture admission, cut equilibrium and adaptive study recovery.
use super::*;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::continuation::controlled_sdf3_gradient_check;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-embedded.fsim"));

fn fixture() -> Spec { spec::parse(EXAMPLE).unwrap() }
fn build(spec: &Spec, level: u8, control: &mut QuadratureControl3<'_>)
    -> std::result::Result<AdaptiveElasticity3, ElasticityError3> {
    let tree = Octree3::uniform(level, 5, spec.leaves).unwrap();
    // This placeholder is empty. The shared native builder must use the actual
    // authored shape and may not invent a supported height-field substitute.
    let placeholder = Domain { height: -1.0, curvature: 0.0 };
    build_operator(spec, HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap(), &tree,
        &placeholder, &IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap(), control)
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-6 * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn embedded_support_grammar_identity_and_translated_axis_selection_are_explicit() {
    let source = fixture();
    assert_eq!(source.fixed.embedded_penalty(), Some(32.0));
    assert_eq!(spec::parse(&source.canonical).unwrap().id, source.id);
    let changed = spec::parse(&EXAMPLE.replace(":penalty 32.0", ":penalty 64.0")).unwrap();
    assert_ne!(source.id, changed.id, "support method parameters bind study identity");
    for bad in [
        EXAMPLE.replace(":axis x", ":axis diagonal"),
        EXAMPLE.replace(":fraction (0.0 0.5)", ":fraction (0.3 0.5)"),
        EXAMPLE.replace(":fraction (0.0 0.5)", ":fraction (0.5 0.0)"),
        EXAMPLE.replace(":fraction (0.0 0.5)", ":fraction (0.5 0.5)"),
        EXAMPLE.replace(":penalty 32.0", ":penalty 0.0"),
        EXAMPLE.replace(":penalty 32.0", ":penalty 10001.0"),
        EXAMPLE.replace(":penalty 32.0", ":penalty 32.0 :displacement (1.0 0.0 0.0)"),
    ] { assert!(spec::parse(&bad).is_err(), "invalid or unsupported support must refuse"); }
    let bounds = ([1.0, -2.0, 3.0], [3.0, 2.0, 9.0]);
    for axis in 0..3 {
        let support = FixedFace::Embedded { axis, fraction: [0.0, 0.5], beta: 32.0 };
        support.validate(1).unwrap();
        let mut inside = bounds.0;
        inside[axis] += 0.25 * (bounds.1[axis] - bounds.0[axis]);
        let mut outside = inside;
        outside[axis] = bounds.0[axis] + 0.75 * (bounds.1[axis] - bounds.0[axis]);
        assert!(support.embedded_contains(inside, bounds));
        assert!(!support.embedded_contains(outside, bounds));
        assert!(!support.contains(bounds.0, bounds), "embedded support pins no box node");
    }
}

#[test]
fn interior_surface_support_reproduces_affine_extension_without_any_box_clamp() {
    let mut spec = fixture();
    spec.poisson = 0.0;
    for level in [1, 2] {
        let mut poll = |_| ControlFlow::Continue(());
        let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
        let op = build(&spec, level, &mut q).unwrap();
        assert!(op.fixed().iter().all(|pin| !pin));
        near(op.embedded_dirichlet_area().unwrap(), 1.0);
        assert_eq!(op.embedded_dirichlet_penalty(), Some(32.0));
        near(op.volumes().iter().sum(), 0.66);
        let law = spec.surfaces[0].unwrap();
        let force = op.surface_load(&|p, n| law.traction(p, n, spec.bounds),
            || ControlFlow::Continue(())).unwrap();
        near(force.resultant[0], 1.0);
        near(force.resultant[1], 0.0);
        near(force.resultant[2], 0.0);
        let solution = op.solve_controlled(&force.rhs, 1e-10, 10000, 16,
            |_| ControlFlow::Continue(())).unwrap();
        near(solution.compliance(), 0.66);
        let physical = op.physical_displacements(solution.coefficients()).unwrap();
        for (p, u) in op.physical_nodes().iter().zip(physical.chunks_exact(3)) {
            near(u[0], p[0] - 0.17);
            near(u[1], 0.0);
            near(u[2], 0.0);
        }
    }
}

#[test]
fn embedded_nitsche_terms_participate_in_the_complete_density_gradient() {
    let spec = fixture();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = build(&spec, 1, &mut q).unwrap();
    with_laws(&spec, |laws| {
        let forces: Vec<_> = laws.iter().map(|law|
            op.reference_load(law.load, || ControlFlow::Continue(())).unwrap()).collect();
        let loads: Vec<_> = laws.iter().zip(&forces).map(|(law, force)|
            LoadCase { force, weight: law.weight }).collect();
        let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000),
            spec.radius, SimpParams { penal: 3.0, beta: 2.0, ..Default::default() });
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.02 * (i % 8) as f64).collect();
        let before = study.operator().elasticity().scales().to_vec();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let audit = controlled_sdf3_gradient_check(&mut study, &loads, &rho,
            Default::default(), &mut control).unwrap();
        assert!(audit.passed(), "{audit:?}");
        assert_eq!(study.operator().elasticity().scales(), before);
        assert!(control.work().linear_iterations > 0);
    });
}

#[test]
fn absent_or_loaded_support_and_interrupted_geometry_never_publish_an_operator() {
    let spec = fixture();
    let conflict = spec::parse(&EXAMPLE.replace(":fraction (0.0 0.5)", ":fraction (0.0 1.0)")).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    assert!(matches!(build(&conflict, 1, &mut q), Err(ElasticityError3::Invalid(
        "surface load overlaps an embedded support; declare disjoint physical patches"))));
    let empty = spec::parse(&EXAMPLE.replace(":offset-m -0.17", ":offset-m -0.67")
        .replace(":offset-m 0.83", ":offset-m 0.93")).unwrap();
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    assert!(matches!(build(&empty, 1, &mut q), Err(ElasticityError3::Invalid(
        "embedded support has no positive numerical area"))));
    let mut q = QuadratureControl3::new(QuadratureOptions3 { max_points: 0, ..Default::default() },
        &mut poll).unwrap();
    assert!(matches!(build(&spec, 1, &mut q),
        Err(ElasticityError3::Quadrature(QuadratureError3::PointBudget))));
    assert_eq!(q.work().points, 0);
    let mut cancel = |work: QuadratureWork3| if work.points > 0 {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut cancel).unwrap();
    assert!(matches!(build(&spec, 1, &mut q),
        Err(ElasticityError3::Cancelled | ElasticityError3::Quadrature(QuadratureError3::Cancelled))));
}

fn retained(ledger: &Ledger, outcome: &Outcome, key: &str, kind: &str) -> Vec<u8> {
    let receipt = JsonValue::parse(&outcome.receipt).unwrap();
    linked(ledger, &receipt, key, kind).unwrap()
}

#[test]
fn embedded_support_survives_adaptive_stages_and_native_ledger_recovery() {
    let spec = fixture();
    let gate = CancelGate::new_clock_free();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = checkpoint::drive(&spec, &whole_db, None, &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(whole.status, "completed");
    let rows = retained(&whole_db, &whole, "iterations", "study-iterations");
    let rows = JsonValue::parse(std::str::from_utf8(&rows).unwrap()).unwrap();
    let stages = rows.get("stages").and_then(JsonValue::as_array).unwrap();
    assert_eq!(stages.len(), 2);
    for stage in stages {
        assert_eq!(stage.get("gradient_check_passed"), Some(&JsonValue::Bool(true)));
        let history = stage.get("history").and_then(JsonValue::as_array).unwrap();
        assert!(history.len() >= 2, "a real supported design update is required");
        for pair in history.windows(2) {
            let c = |row: &JsonValue| row.get("compliance_j").and_then(JsonValue::as_f64).unwrap();
            assert!(c(&pair[1]) <= c(&pair[0]));
        }
    }
    let segment_db = Ledger::open(":memory:").unwrap();
    let first = checkpoint::drive(&spec, &segment_db, Some(1), &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(first.status, "budget-exhausted");
    let old = load(&segment_db, &first.pointer).unwrap();
    let resumed = checkpoint::drive(&spec, &segment_db, None, &gate, Some(&old), |_, _| Ok(())).unwrap();
    assert_eq!(resumed.status, "completed");
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(retained(&whole_db, &whole, key, kind), retained(&segment_db, &resumed, key, kind));
    }
    let design = retained(&whole_db, &whole, "design", "study-design");
    let design = JsonValue::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    let fields = design.get("displacements_m").and_then(JsonValue::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);
}
