use super::*;
use fs_cutfem::elastic3::surface::SurfaceForce3;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-preload.fsim"
));
const DEAD: &str = "  (affine-dead-load\n    :body-n-m3 (0.0 0.0 -0.1)\n    :surface (pressure :pa 0.1 :x-fraction (0.5 1.0)))\n";

fn fixed_source() -> String {
    let (source, _) = FIXTURE.split_once("  (affine-grid-assessment").unwrap();
    format!("{source})\n")
        .replace(":maximum-level 3", ":maximum-level 1")
        .replace(":maximum-leaves 512", ":maximum-leaves 8")
}

fn compute(spec: &Spec) -> Computation {
    compute_observed(spec, &CancelGate::new_clock_free(), 0, None, 0.0, |_, _| Ok(()))
        .unwrap()
}

fn same_vector(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    let scale = actual.iter().chain(expected).map(|v| v.abs()).fold(1e-12, f64::max);
    for (a, b) in actual.iter().zip(expected) {
        assert!((a - b).abs() <= 2e-8 * scale, "{a:e} != {b:e}");
    }
}

// Integrate the authored fixed law independently of the CLI adapter and of
// the material body-load map. Neither stiffness nor occupancy enters it.
fn fixed_force(spec: &Spec, op: &AdaptiveElasticity3) -> Vec<f64> {
    let dead = spec.affine_dead_load.unwrap();
    let body = |_: [f64; 3]| dead.body_n_m3;
    let surface = |p, n| dead.surface.unwrap().traction(p, n, spec.bounds);
    op.reference_load(ReferenceLoad3 {
        body: Some(&body),
        surface: dead.surface.map(|_| SurfaceForce3::Traction(&surface)),
    }, || ControlFlow::Continue(())).unwrap()
}

fn weight(spec: &Spec, op: &AdaptiveElasticity3, rho: &[f64]) -> Vec<f64> {
    let law = spec.material_acceleration.unwrap();
    op.prepare_material_body_load(&|_| law.acceleration_m_s2, law.max_terms,
        || ControlFlow::Continue(())).unwrap()
        .load(rho, || ControlFlow::Continue(())).unwrap()
}

fn check_offset(spec: &Spec, op: &AdaptiveElasticity3, rho: &[f64], actual: &[f64]) {
    let fixed = fixed_force(spec, op);
    let body = weight(spec, op, rho);
    let expected: Vec<_> = fixed.iter().zip(body).map(|(f, b)| f + b).collect();
    same_vector(actual, &expected);
}

fn build(spec: &Spec, level: u8, q: &mut QuadratureControl3<'_>)
    -> std::result::Result<AdaptiveElasticity3, ElasticityError3> {
    let tree = Octree3::uniform(level, spec.max_level as u8, spec.leaves).unwrap();
    let domain = PhysicalDomain { bounds: spec.bounds, height: spec.height, curvature: spec.curvature };
    loading::build_operator(spec, HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap(),
        &tree, &domain, &loading::material(spec).unwrap(), q)
}

#[test]
fn g1_affine_preload_is_fixed_while_weight_tracks_eroded_material() {
    let mut spec = spec::parse(&fixed_source()).unwrap();
    let mut offsets = Vec::new();
    let mut weights = Vec::new();
    for density in [0.65, 0.8] {
        spec.density = density;
        let run = compute(&spec);
        let e = &run.state.accepted.evaluation;
        let op = run.study.operator().elasticity();
        assert!(run.state.gradient.passed());
        assert_eq!(e.policy.dimension(), 2);
        assert!(e.responses.is_some() && e.stress.is_some());
        check_offset(&spec, op, &e.physical.eroded_rho, &e.offset_force);
        let body = weight(&spec, op, &e.physical.eroded_rho);
        offsets.push(e.offset_force.iter().zip(&body).map(|(f, b)| f - b).collect::<Vec<_>>());
        weights.push(body);
        loading::with_laws(&spec, |laws| {
            let basis: Vec<_> = laws.iter().map(|law|
                op.reference_load(law.load, || ControlFlow::Continue(())).unwrap()).collect();
            let force: Vec<_> = e.offset_force.iter().enumerate().map(|(d, b)| {
                b + basis.iter().zip(&e.worst_coefficients).map(|(f, q)| f[d] * q).sum::<f64>()
            }).collect();
            same_vector(&e.worst_force, &force);
        });
    }
    same_vector(&offsets[0], &offsets[1]);
    assert!(offsets[0].iter().any(|f| f.abs() > 1e-5));
    assert!(weights[0].iter().zip(&weights[1]).any(|(a, b)| (a - b).abs() > 1e-5));
}

#[test]
fn g1_affine_preload_alone_retains_reference_surface_rules_and_rejects_support_overlap() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    assert!(spec.surfaces.iter().all(Option::is_none));
    spec.fixed = FixedFace::Left;
    let pressure = spec.affine_dead_load.unwrap().surface.unwrap();
    let (lo, hi) = spec.bounds;
    let (lx, ly) = (hi[0] - lo[0], hi[1] - lo[1]);
    let expected = [-0.1 * ly * spec.curvature * lx * lx / 4.0, 0.0, -0.1 * ly * lx / 2.0];
    for level in [1, 2] {
        let mut poll = |_| ControlFlow::Continue(());
        let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
        let op = build(&spec, level, &mut q).unwrap();
        let result = op.surface_load(&|p, n| pressure.traction(p, n, spec.bounds),
            || ControlFlow::Continue(())).unwrap();
        same_vector(&result.resultant, &expected);
        assert!(q.work().points > 0);
        let body = op.body_load(&|_| spec.affine_dead_load.unwrap().body_n_m3,
            || ControlFlow::Continue(())).unwrap();
        let total: Vec<_> = body.iter().zip(&result.rhs).map(|(a, b)| a + b).collect();
        loading::with_affine_laws(&spec, |dead, basis| {
            assert_eq!(basis.len(), 2);
            let actual = op.reference_load(dead, || ControlFlow::Continue(())).unwrap();
            same_vector(&actual, &total);
        });
    }
    let overlap = spec::parse(&FIXTURE.replace(":fraction (0.0 0.5)", ":fraction (0.0 1.0)")).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    assert!(matches!(build(&overlap, 1, &mut q), Err(ElasticityError3::Invalid(
        "surface load overlaps an embedded support; declare disjoint physical patches"))));
}

#[test]
fn g1_g3_affine_preload_is_reintegrated_in_the_fine_assessment_and_refit() {
    let spec = spec::parse(FIXTURE).unwrap();
    let source = compute(&spec);
    let assessment = source.state.assessment.as_ref().unwrap();
    assert_eq!(assessment.status, "checked", "{:?}", assessment.reason);
    let report = assessment.report.as_ref().unwrap();
    assert_eq!(report.goal, Some(spec.affine_load_space.unwrap().goal));
    check_offset(&spec, source.study.operator().elasticity(),
        &report.coarse_eroded_rho, &report.coarse.offset_force);
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let fine = build(&spec, (spec.level + 1) as u8, &mut q).unwrap();
    assert!(report.fine.offset_force.len() > report.coarse.offset_force.len());
    check_offset(&spec, &fine, &report.fine_eroded_rho, &report.fine.offset_force);
    assert!(report.responses.is_some() && report.stress.is_some());

    let gate = CancelGate::new_clock_free();
    let context = grid_refinement::Context::new(&source, hash_bytes(b"preload source")).unwrap();
    let prepared = grid_refinement::prepare(&spec, &source, context, &gate).result.unwrap();
    let candidate_spec = grid_refinement::effective_spec(&spec, Some(&prepared.context)).unwrap();
    let candidate = compute_on_grid(&candidate_spec, &gate, 0, None, 0.0,
        Some(prepared), |_, _| Ok(())).unwrap();
    let e = &candidate.state.accepted.evaluation;
    assert!(candidate.study.cells() > source.study.cells());
    assert!(candidate.state.gradient.passed());
    assert_eq!(e.policy.dimension(), 2);
    check_offset(&candidate_spec, candidate.study.operator().elasticity(),
        &e.physical.eroded_rho, &e.offset_force);
}

#[test]
fn g3_affine_preload_admission_zero_compatibility_and_resume_keep_the_authored_law() {
    let source = fixed_source().replace(":max-updates 80", ":max-updates 2");
    assert!(source.contains(DEAD));
    let omitted = spec::parse(&source.replace(DEAD, "")).unwrap();
    let zero = spec::parse(&source.replace(DEAD,
        "  (affine-dead-load :body-n-m3 (0.0 0.0 0.0) :surface none)\n")).unwrap();
    assert_ne!(omitted.id, zero.id);
    let a = compute(&omitted);
    let b = compute(&zero);
    assert_eq!(a.state.accepted.evaluation.worst_force, b.state.accepted.evaluation.worst_force);
    assert_eq!(a.state.accepted.evaluation.worst_displacement, b.state.accepted.evaluation.worst_displacement);
    assert_eq!(a.state.history, b.state.history);
    loading::with_affine_laws(&zero, |dead, basis| {
        assert!(dead.body.is_none() && dead.surface.is_none());
        assert_eq!(basis.len(), 2);
    });
    for (from, to) in [
        (":body-n-m3 (0.0 0.0 -0.1)", ":body-n-m3 (0.0 0.0 1e13)"),
        ("pressure :pa 0.1", "pressure :pa 0.0"),
        ("pressure :pa 0.1", "follower-pressure :pa 0.1"),
        ("pressure :pa 0.1 :x-fraction (0.5 1.0)", "pressure :pa 0.1 :x-fraction (0.3 1.0)"),
    ] {
        assert!(spec::parse(&source.replace(from, to)).is_err(), "admitted {to}");
    }
    assert!(spec::parse(&source.replace(DEAD, &format!("{DEAD}{DEAD}"))).is_err());
    let late = source.replace(DEAD, "").replace("  (affine-responses", &format!("{DEAD}  (affine-responses"));
    assert!(spec::parse(&late).is_err());
    let traction = spec::parse(&source.replace("pressure :pa 0.1",
        "traction :pa (0.0 0.02 -0.1)")).unwrap();
    assert!(matches!(traction.affine_dead_load.unwrap().surface.unwrap().force,
        loading::Force::Traction([0.0, 0.02, -0.1])));
    let ordinary = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-3d-adaptive.fsim")).trim_end();
    let ordinary = format!("{}\n  (material-acceleration :density-kg-m3 0.2 :acceleration-m-s2 (0.0 0.0 -1.0) :max-stencil-terms 100000)\n",
        ordinary.strip_suffix(')').unwrap());
    assert!(spec::parse(&format!("{ordinary})")).is_ok());
    assert!(spec::parse(&format!("{ordinary}{DEAD})")).is_err());

    let spec = spec::parse(&source).unwrap();
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    let gate = CancelGate::new_clock_free();
    let ledger = Ledger::open(":memory:").unwrap();
    let whole = drive(&spec, &ledger, Some(2), &gate).unwrap();
    let whole = load(&ledger, &whole.pointer).unwrap();
    let first = drive(&spec, &ledger, Some(1), &gate).unwrap();
    let first = load(&ledger, &first.pointer).unwrap();
    let continued = resume(&ledger, &first, Some(1), &gate).unwrap();
    let continued = load(&ledger, &continued.pointer).unwrap();
    assert_eq!(linked(&ledger, &whole.value, "design", "study-design").unwrap(),
        linked(&ledger, &continued.value, "design", "study-design").unwrap());
    let bytes = linked(&ledger, &continued.value, "report_json", "study-report-json").unwrap();
    let report = JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap();
    let model = report.get("model").unwrap();
    assert_eq!(model.str_field("dead_offset"), Some("authored-reference-load"));
    assert_eq!(model.str_field("force_set"), Some("f=f_dead+B*rho_eroded+Fq, ||q||_2<=1"));
    let dead = model.get("affine_dead_load").unwrap();
    assert_eq!(dead.get("surface").unwrap().str_field("kind"), Some("pressure"));
    assert_eq!(dead.get("body_n_m3").unwrap().as_array().unwrap()[2].as_f64(), Some(-0.1));
    let changed = spec::parse(&source.replace(":body-n-m3 (0.0 0.0 -0.1)", ":body-n-m3 (0.0 0.0 -0.2)")).unwrap();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    assert!(resume::recover(&changed, &ledger, &first, producer).is_err());
}
