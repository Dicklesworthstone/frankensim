use super::*;

const EXAMPLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-pressure.fsim"));
fn unit_source() -> String {
    EXAMPLE.replace("(0.1 0.05 0.04)", "(1.0 1.0 1.0)")
        .replace(":height-m 0.028", ":height-m 0.7")
        .replace(":curvature-per-m 1.0", ":curvature-per-m 0.1")
        .replace(":filter-radius-m 0.01", ":filter-radius-m 0.15")
        .replace(":youngs-modulus-pa 70000000000.0", ":youngs-modulus-pa 1.0")
        .replace("(0.0 0.0 -26000.0)", "(0.0 0.0 0.0)")
        .replace(":pa 1000.0", ":pa 1.0")
        .replace(":pa (0.0 1000.0 0.0)", ":pa (0.0 1.0 0.0)")
        .replace(":updates-per-stage 3", ":updates-per-stage 1")
        .replace("((1.0 1.0) (2.0 2.0) (3.0 8.0))", "((1.0 1.0) (2.0 2.0))")
}
fn build(spec: &Spec, level: u8, control: &mut QuadratureControl3<'_>)
    -> std::result::Result<AdaptiveElasticity3, ElasticityError3> {
    let domain = PhysicalDomain { bounds: spec.bounds, height: spec.height, curvature: spec.curvature };
    let tree = Octree3::uniform(level, 5, spec.leaves).unwrap();
    build_operator(spec, HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap(), &tree, &domain,
        &IsotropicElastic::new(spec.youngs, spec.poisson, 1.0).unwrap(), control)
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-8 * scale.max(1e-20), "{actual:e} != {expected:e}");
}

#[test]
fn explicit_mixed_loads_preserve_source_identity_and_refuse_missing_or_aliased_physics() {
    let spec = spec::parse(EXAMPLE).unwrap();
    assert_eq!(spec.loads.len(), 2);
    assert!(matches!(spec.surfaces[0].unwrap().force, Force::Pressure(1000.0)));
    assert!(matches!(spec.surfaces[1].unwrap().force, Force::Traction(_)));
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    for bad in [
        EXAMPLE.replace(":pa 1000.0", ":pa 0.0"),
        EXAMPLE.replace(":pa (0.0 1000.0 0.0)", ":pa (0.0 0.0 0.0)"),
        EXAMPLE.replace(":x-fraction (0.5 1.0)", ":x-fraction (0.3 1.0)"),
        EXAMPLE.replace(":x-fraction (0.5 1.0)", ":x-fraction (1.0 0.5)"),
        EXAMPLE.replace(":surface (pressure", ":ignored-surface (pressure"),
        EXAMPLE.replace("(pressure :pa", "(follower-pressure :pa"),
        EXAMPLE.replace(":weight 0.7", ":weight 0.0"),
    ] { assert!(spec::parse(&bad).is_err()); }
    let changed = spec::parse(&EXAMPLE.replace(":pa 1000.0", ":pa -1000.0")).unwrap();
    assert_ne!(changed.id, spec.id);
    let body_only = EXAMPLE.replace("(pressure :pa 1000.0 :x-fraction (0.5 1.0))", "none")
        .replace("(traction :pa (0.0 1000.0 0.0) :x-fraction (0.5 1.0))", "none");
    assert!(spec::parse(&body_only).unwrap().surfaces.iter().all(Option::is_none));
    assert!(spec::parse(&body_only.replace("(0.0 0.0 -26000.0)", "(0.0 0.0 0.0)")).is_err());
}

#[test]
fn curved_patch_pressure_has_correct_force_and_refinement_does_not_change_its_law() {
    let spec = spec::parse(EXAMPLE).unwrap();
    let pressure = spec.surfaces[0].unwrap();
    let (lo, hi) = spec.bounds;
    let lx = hi[0] - lo[0];
    let ly = hi[1] - lo[1];
    let expected = [-1000.0 * ly * spec.curvature * lx * lx / 4.0, 0.0, -1000.0 * ly * lx / 2.0];
    for level in [1, 2] {
        let mut poll = |_| ControlFlow::Continue(());
        let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
        let op = build(&spec, level, &mut q).unwrap();
        let result = op.surface_load(&|p, n| pressure.traction(p, n, spec.bounds),
            || ControlFlow::Continue(())).unwrap();
        for axis in 0..3 { near(result.resultant[axis], expected[axis], expected[2].abs()); }
        assert!(q.work().points > 0);
        // Mixed source loading is a sum inside ONE load case, not an average.
        with_laws(&spec, |laws| {
            let mixed = op.reference_load(laws[0].load, || ControlFlow::Continue(())).unwrap();
            let body = op.body_load(&|_| spec.loads[0].0, || ControlFlow::Continue(())).unwrap();
            for ((a, b), c) in mixed.iter().zip(&body).zip(&result.rhs) {
                near(*a, b + c, expected[2].abs());
            }
        });
    }
}

#[test]
fn flat_pressure_and_traction_agree_including_translated_force_moments() {
    let mut spec = spec::parse(&unit_source()).unwrap();
    spec.bounds = ([0.25, -0.5, 1.0], [1.25, 0.5, 2.0]);
    spec.curvature = 0.0;
    spec.surfaces[1] = Some(SurfaceSpec { force: Force::Traction([0.0, 0.0, -1.0]), x_fraction: [0.5, 1.0] });
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = build(&spec, 1, &mut q).unwrap();
    with_laws(&spec, |laws| {
        let a = op.reference_load(laws[0].load, || ControlFlow::Continue(())).unwrap();
        let b = op.reference_load(laws[1].load, || ControlFlow::Continue(())).unwrap();
        assert_eq!(a, b);
    });
    let pressure = spec.surfaces[0].unwrap();
    let force = op.surface_load(&|p, n| pressure.traction(p, n, spec.bounds),
        || ControlFlow::Continue(())).unwrap();
    near(force.resultant[2], -0.5, 1.0);
    near(force.moment[0], 0.0, 1.0);
    near(force.moment[1], 0.5, 1.0); // x centroid is 1.0 m.
    near(force.moment[2], 0.0, 1.0);
}

#[test]
fn surface_quadrature_obeys_shared_point_budget_and_cancellation() {
    let spec = spec::parse(EXAMPLE).unwrap();
    let mut body = spec.clone();
    body.surfaces.fill(None);
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    build(&body, 1, &mut q).unwrap();
    let bulk_points = q.work().points;
    let mut q = QuadratureControl3::new(QuadratureOptions3 {
        max_points: bulk_points, ..Default::default()
    }, &mut poll).unwrap();
    assert!(matches!(build(&spec, 1, &mut q),
        Err(ElasticityError3::Quadrature(QuadratureError3::PointBudget))));
    assert_eq!(q.work().points, bulk_points);
    let mut cancel = |work: QuadratureWork3| if work.points > bulk_points {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut cancel).unwrap();
    assert!(matches!(build(&spec, 1, &mut q),
        Err(ElasticityError3::Cancelled | ElasticityError3::Quadrature(QuadratureError3::Cancelled))));
}

fn bytes(ledger: &Ledger, outcome: &Outcome, key: &str, kind: &str) -> Vec<u8> {
    let receipt = JsonValue::parse(&outcome.receipt).unwrap();
    linked(ledger, &receipt, key, kind).unwrap()
}

#[test]
fn pressure_and_shear_studies_refine_check_gradients_and_resume_exact_physical_fields() {
    let spec = spec::parse(&unit_source()).unwrap();
    let gate = CancelGate::new_clock_free();
    let whole_db = Ledger::open(":memory:").unwrap();
    let whole = checkpoint::drive(&spec, &whole_db, None, &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(whole.status, "completed");
    let summary = JsonValue::parse(std::str::from_utf8(&bytes(&whole_db, &whole, "report_json", "study-report-json")).unwrap()).unwrap();
    let refinement = summary.get("refinements").and_then(JsonValue::as_array).unwrap();
    assert_eq!(refinement.len(), 1);
    assert_eq!(refinement[0].get("installed"), Some(&JsonValue::Bool(true)));
    let iterations = JsonValue::parse(std::str::from_utf8(&bytes(&whole_db, &whole, "iterations", "study-iterations")).unwrap()).unwrap();
    for stage in iterations.get("stages").and_then(JsonValue::as_array).unwrap() {
        assert_eq!(stage.get("gradient_check_passed"), Some(&JsonValue::Bool(true)));
        for row in stage.get("history").and_then(JsonValue::as_array).unwrap() {
            let cases = row.get("case_compliances_j").and_then(JsonValue::as_array).unwrap();
            let weighted = 0.7 * cases[0].as_f64().unwrap() + 0.3 * cases[1].as_f64().unwrap();
            near(row.get("compliance_j").and_then(JsonValue::as_f64).unwrap(), weighted, weighted);
            assert!(row.get("volume_fraction").and_then(JsonValue::as_f64).unwrap() <= 0.50000001);
        }
    }
    let resumed_db = Ledger::open(":memory:").unwrap();
    let first = checkpoint::drive(&spec, &resumed_db, Some(1), &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(first.status, "budget-exhausted");
    let loaded = load(&resumed_db, &first.pointer).unwrap();
    let resumed = checkpoint::drive(&spec, &resumed_db, None, &gate, Some(&loaded), |_, _| Ok(())).unwrap();
    assert_eq!(resumed.status, "completed");
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(bytes(&whole_db, &whole, key, kind), bytes(&resumed_db, &resumed, key, kind));
    }
    let design = JsonValue::parse(std::str::from_utf8(&bytes(&resumed_db, &resumed, "design", "study-design")).unwrap()).unwrap();
    let fields = design.get("displacements_m").and_then(JsonValue::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);
}
