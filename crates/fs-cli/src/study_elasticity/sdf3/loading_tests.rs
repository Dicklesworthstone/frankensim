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
        &material(spec).unwrap(), control)
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
fn a_refused_study_names_the_cause_of_its_stop() {
    // E = 1 Pa puts the first displacement outside the material's small-strain domain. The stop names its stage, and the report carries the refusal behind it.
    let spec = spec::parse(&unit_source()).unwrap();
    let gate = CancelGate::new_clock_free();
    let db = Ledger::open(":memory:").unwrap();
    let run = checkpoint::drive(&spec, &db, None, &gate, None, |_, _| Ok(())).unwrap();
    assert_eq!(run.status, "numerical-failure");
    let summary = JsonValue::parse(std::str::from_utf8(&bytes(&db, &run, "report_json", "study-report-json")).unwrap()).unwrap();
    let stop = summary.get("stop").unwrap();
    assert_eq!(stop.str_field("evaluation"), Some("Breakdown { stage: \"sdf3-material-domain\" }"));
    let cause = stop.str_field("cause").expect("the refusal behind the stop is reported");
    assert!(cause.contains("affine bulk strain ball exceeds retained material domain"), "{cause}");
}

#[test]
fn pressure_and_shear_studies_refine_check_gradients_and_resume_exact_physical_fields() {
    // E = 2^20 Pa. At the 1 Pa of `unit_source` the unit pressure puts the first displacement outside the material's small-strain domain and the study is refused
    // (`a_refused_study_names_the_cause_of_its_stop` pins that refusal); the checks below are relative, so no other literal moves with the modulus.
    let spec = spec::parse(&unit_source().replace(":youngs-modulus-pa 1.0", ":youngs-modulus-pa 1048576.0")).unwrap();
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

fn with_acceleration(source: &str) -> String {
    format!("{}\n  (material-acceleration :density-kg-m3 2.0 :acceleration-m-s2 (0.0 0.0 -0.5) :max-stencil-terms 100000))",
        source.trim_end().strip_suffix(')').unwrap())
}

#[test]
fn material_acceleration_is_explicit_identity_bound_and_admitted_before_work() {
    let source = with_acceleration(&unit_source());
    let spec = spec::parse(&source).unwrap();
    let law = spec.material_acceleration.unwrap();
    assert_eq!(law.density_kg_m3, 2.0);
    assert_eq!(law.acceleration_m_s2, [0.0, 0.0, -0.5]);
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    assert_ne!(spec::parse(&source.replace(":density-kg-m3 2.0", ":density-kg-m3 3.0")).unwrap().id, spec.id);
    assert_ne!(spec::parse(&source.replace("(0.0 0.0 -0.5)", "(0.0 0.0 -0.6)")).unwrap().id, spec.id);
    for source in [
        source.replace(":density-kg-m3 2.0", ":density-kg-m3 0.0"),
        source.replace("(0.0 0.0 -0.5)", "(0.0 0.0 0.0)"),
        source.replace(":acceleration-m-s2", ":body-force-n-m3"),
        source.replace(":max-stencil-terms 100000", ":max-stencil-terms 0"),
        source.replace(":max-stencil-terms 100000", ":max-stencil-terms 4000001"),
        with_acceleration(&source),
    ] { assert!(spec::parse(&source).is_err()); }

    let gravity_only = source
        .replace("(pressure :pa 1.0 :x-fraction (0.5 1.0))", "none")
        .replace("(traction :pa (0.0 1.0 0.0) :x-fraction (0.5 1.0))", "none");
    let gravity_only = spec::parse(&gravity_only).unwrap();
    assert!(gravity_only.loads.iter().all(|(f, _)| *f == [0.0; 3]));
    assert!(gravity_only.surfaces.iter().all(Option::is_none));

    let load_space = include_str!("../../../../../examples/marquee/bracket-3d-load-space.fsim");
    let error = spec::parse(&with_acceleration(load_space)).unwrap_err();
    assert!(error.message.contains("zero-centered"));
}

#[test]
fn material_acceleration_uses_native_mass_occupancy_and_full_compliance_gradient() {
    use fs_solver::op::LinearOp;
    use fs_topopt::pipeline::LoadCase;

    let source = unit_source()
        .replace(":youngs-modulus-pa 1.0", ":youngs-modulus-pa 1000000.0")
        .replace("(pressure :pa 1.0 :x-fraction (0.5 1.0))", "none")
        .replace("(traction :pa (0.0 1.0 0.0) :x-fraction (0.5 1.0))", "none");
    let source = format!("{}\n (design-regions :solid (((0.0 0.0 0.0) (0.5 0.5 0.5))) :void (((0.0 0.5 0.5) (0.5 1.0 1.0)))))",
        source.trim_end().strip_suffix(')').unwrap());
    let source = with_acceleration(&source);
    let spec = spec::parse(&source).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = build(&spec, 1, &mut q).unwrap();
    for (mass, volume) in op.reference_cell_masses_kg().unwrap().iter().zip(op.volumes()) {
        near(*mass, 2.0 * volume, *mass);
    }
    let native = op.prepare_material_body_load(&|_| [0.0, 0.0, -0.5], 100000,
        || ControlFlow::Continue(())).unwrap();
    let study = regions::bind(CutDensityStudy3::new(op, spec.radius, spec.schedule[0]), &spec).unwrap();
    let mut study = bind_material_acceleration(study, &spec, &CancelGate::new(),
        || ControlFlow::Continue(())).unwrap();
    assert!(study.has_material_body_load());
    assert!(study.can_refine_material_body_load());
    assert_eq!(study.material_body_load_fingerprint(), Some(native.fingerprint()));
    let map = study.material_body_load().unwrap();
    let full = map.load(&vec![1.0; study.cells()], || ControlFlow::Continue(())).unwrap();
    let half = map.load(&vec![0.5; study.cells()], || ControlFlow::Continue(())).unwrap();
    let empty = map.load(&vec![0.0; study.cells()], || ControlFlow::Continue(())).unwrap();
    assert!(full.iter().any(|f| f.abs() > 1e-4));
    for ((full, half), empty) in full.iter().zip(half).zip(empty) {
        assert_eq!(*full * 0.5, half);
        assert_eq!(empty, 0.0);
    }
    let n = study.operator().n();
    let zero = vec![0.0; n];
    let cases = [LoadCase { force: &zero, weight: 1.0 }];
    let raw: Vec<_> = (0..study.cells()).map(|i| 0.45 + 0.02 * i as f64).collect();
    let mut checkpoint = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut checkpoint);
    let value = study.evaluate(&raw, &cases, &mut control).unwrap();
    assert!(value.objective.compliance > 0.0);
    assert!(value.projected_rho.contains(&0.0));
    assert!(value.projected_rho.contains(&1.0));
    let physical_force = study.material_body_load().unwrap().load(&value.projected_rho,
        || ControlFlow::Continue(())).unwrap();
    let defect = study.operator().field_residual(&value.objective.displacements[0], &physical_force,
        || ControlFlow::Continue(())).unwrap();
    assert!(defect < 1e-8, "equilibrium must include the actual projected mass");

    // G3: unchanged stiffness and occupancy, twice the source mass, four times
    // f^T u. Rebuild the declared physical law without an OC sign assumption.
    let mass_spec = spec::parse(&source.replace(":density-kg-m3 2.0", ":density-kg-m3 4.0")).unwrap();
    let heavy_op = build(&mass_spec, 1, &mut q).unwrap();
    let heavy = regions::bind(CutDensityStudy3::new(heavy_op, mass_spec.radius,
        mass_spec.schedule[0]), &mass_spec).unwrap();
    let mut heavy = bind_material_acceleration(heavy, &mass_spec, &CancelGate::new(),
        || ControlFlow::Continue(())).unwrap();
    let doubled = heavy.evaluate(&raw, &cases, &mut control).unwrap();
    assert_eq!(doubled.projected_rho, value.projected_rho);
    near(doubled.objective.compliance, 4.0 * value.objective.compliance,
        4.0 * value.objective.compliance);

    let epsilon = 1e-5;
    let direction: Vec<_> = (0..raw.len()).map(|i| 0.2 + 0.03 * i as f64).collect();
    let shifted = |sign: f64| raw.iter().zip(&direction).map(|(r, d)| r + sign * epsilon * d).collect::<Vec<_>>();
    let high = study.evaluate(&shifted(1.0), &cases, &mut control).unwrap().objective.compliance;
    let low = study.evaluate(&shifted(-1.0), &cases, &mut control).unwrap().objective.compliance;
    let actual = (high - low) / (2.0 * epsilon);
    let analytic: f64 = value.objective.gradient.iter().zip(&direction).map(|(g, d)| g * d).sum();
    assert!((actual - analytic).abs() < 2e-5 * actual.abs().max(1e-8),
        "material force derivative must pass through the shared filter/projection: {actual} vs {analytic}");
}

#[test]
fn material_acceleration_preparation_preserves_cancellation_and_stencil_limits() {
    let mut spec = spec::parse(&with_acceleration(&unit_source())).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = build(&spec, 1, &mut q).unwrap();
    let study = CutDensityStudy3::new(op, spec.radius, spec.schedule[0]);
    let gate = CancelGate::new();
    gate.request();
    let error = match bind_material_acceleration(study, &spec, &gate, || ControlFlow::Break(())) {
        Ok(_) => panic!("cancelled preparation must not publish a study"), Err(e) => e,
    };
    assert_eq!(error.exit, exit::CANCELLED);
    let op = build(&spec, 1, &mut q).unwrap();
    let study = CutDensityStudy3::new(op, spec.radius, spec.schedule[0]);
    spec.material_acceleration.as_mut().unwrap().max_terms = 1;
    let error = match bind_material_acceleration(study, &spec, &CancelGate::new(),
        || ControlFlow::Continue(())) {
        Ok(_) => panic!("bounded stencil preparation must refuse"), Err(e) => e,
    };
    assert_eq!(error.exit, exit::BUDGET);
}

fn with_formulation(source: &str, formulation: &str, poisson: &str) -> String {
    source.replace(":poissons-ratio 0.3", &format!(
        ":poissons-ratio {poisson} :formulation {formulation}"))
}

#[test]
fn bbar_is_explicit_identity_bound_and_nearly_incompressible_admission_is_bounded() {
    use fs_cutfem::elastic3::ElasticityFormulation3;

    let legacy = spec::parse(EXAMPLE).unwrap();
    assert_eq!(legacy.formulation, ElasticityFormulation3::Displacement);
    let explicit = spec::parse(&with_formulation(EXAMPLE, "displacement", "0.3")).unwrap();
    assert_eq!(explicit.formulation, legacy.formulation);
    assert_ne!(explicit.id, legacy.id, "an explicit declaration remains part of the source");
    let source = with_formulation(EXAMPLE, "bbar", "0.4999");
    let bbar = spec::parse(&source).unwrap();
    assert_eq!(bbar.formulation, ElasticityFormulation3::BBar);
    assert_eq!(bbar.poisson, 0.4999);
    assert_eq!(spec::parse(&bbar.canonical).unwrap().id, bbar.id);
    assert_ne!(bbar.id, legacy.id);
    assert_ne!(bbar.id, spec::parse(&with_formulation(EXAMPLE, "bbar", "0.49")).unwrap().id);
    for bad in [
        with_formulation(EXAMPLE, "automatic", "0.3"),
        with_formulation(EXAMPLE, "bbar", "0.5"),
        with_formulation(EXAMPLE, "bbar", "0.49999"),
        with_formulation(EXAMPLE, "bbar", "-0.01"),
        with_formulation(EXAMPLE, "displacement", "0.49"),
        EXAMPLE.replace(":poissons-ratio 0.3", ":poissons-ratio 0.49"),
        source.replace(":formulation bbar", ":formulation bbar :formulation displacement"),
        source.replace(":formulation bbar", ":ignored-formulation bbar"),
    ] { assert!(spec::parse(&bad).is_err(), "unsupported formulations must refuse before physics"); }
}

#[test]
fn bbar_reaches_bulk_surface_embedded_and_mass_backed_builders_on_each_grid() {
    use fs_cutfem::elastic3::ElasticityFormulation3;

    let surface = with_formulation(EXAMPLE, "bbar", "0.4999");
    let body = surface
        .replace("(pressure :pa 1000.0 :x-fraction (0.5 1.0))", "none")
        .replace("(traction :pa (0.0 1000.0 0.0) :x-fraction (0.5 1.0))", "none");
    let embedded = body.replace(":fixed-boundary left",
        ":fixed-boundary (embedded :axis x :fraction (0.0 0.5) :penalty 32.0)");
    for source in [&body, &surface, &embedded] {
        for native_material in [false, true] {
            let source = if native_material { with_acceleration(source) } else { source.clone() };
            let spec = spec::parse(&source).unwrap();
            for level in [1, 2] {
                let mut poll = |_| ControlFlow::Continue(());
                let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
                let op = build(&spec, level, &mut q).unwrap();
                assert_eq!(op.formulation(), ElasticityFormulation3::BBar);
                assert_eq!(op.reference_cell_masses_kg().is_some(), native_material);
                assert!(op.volumes().iter().sum::<f64>() > 0.0);
                if spec.fixed.embedded_penalty().is_some() {
                    assert!(op.embedded_dirichlet_area().unwrap() > 0.0);
                }
            }
        }
    }
}

#[test]
fn bbar_changes_native_volumetric_work_and_explicit_displacement_preserves_legacy_bits() {
    use fs_solver::op::LinearOp;

    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let legacy = build(&spec::parse(EXAMPLE).unwrap(), 1, &mut q).unwrap();
    let explicit = build(&spec::parse(&with_formulation(EXAMPLE, "displacement", "0.3")).unwrap(), 1, &mut q).unwrap();
    let bbar = build(&spec::parse(&with_formulation(EXAMPLE, "bbar", "0.3")).unwrap(), 1, &mut q).unwrap();
    assert_eq!(legacy.nodes(), bbar.nodes());
    let field: Vec<_> = legacy.nodes().iter().flat_map(|&[x, y, z]|
        [x * y, x * z, 0.0]).collect();
    let mut full_action = vec![0.0; legacy.n()];
    let mut explicit_action = full_action.clone();
    let mut bbar_action = full_action.clone();
    legacy.apply(&field, &mut full_action);
    explicit.apply(&field, &mut explicit_action);
    bbar.apply(&field, &mut bbar_action);
    assert_eq!(full_action, explicit_action);
    let work = |action: &[f64]| field.iter().zip(action).map(|(u, f)| u * f).sum::<f64>();
    let full_work = work(&full_action);
    let bbar_work = work(&bbar_action);
    assert!(bbar_work > 0.0 && bbar_work < full_work * (1.0 - 1e-6),
        "cell-averaged dilation must change the actual operator, not only metadata: {bbar_work} vs {full_work}");
}

#[test]
fn bbar_nearly_incompressible_study_keeps_real_compliance_and_volume_gradients() {
    use fs_topopt::pipeline::LoadCase;
    use fs_topopt::sdf3::continuation::controlled_sdf3_gradient_check;

    let source = unit_source().replace(":youngs-modulus-pa 1.0", ":youngs-modulus-pa 1000000.0");
    let spec = spec::parse(&with_formulation(&source, "bbar", "0.4999")).unwrap();
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
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.4 + 0.02 * (i % 8) as f64).collect();
        let before = study.operator().elasticity().scales().to_vec();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget {
            total_iterations: 100_000, per_solve_iterations: 20_000,
        }, &mut poll);
        let check = controlled_sdf3_gradient_check(&mut study, &loads, &rho,
            Default::default(), &mut control).unwrap();
        assert!(check.passed(), "nearly incompressible native derivatives: {check:?}");
        assert!(check.baseline_compliance > 0.0);
        assert!(check.work.linear_iterations > 0);
        assert_eq!(study.operator().elasticity().scales(), before);
    });
}
