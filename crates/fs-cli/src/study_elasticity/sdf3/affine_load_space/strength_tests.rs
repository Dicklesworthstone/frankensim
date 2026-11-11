use super::*;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-strength.fsim"
));

fn fixture() -> Spec {
    spec::parse(FIXTURE).unwrap()
}

fn compute(spec: &Spec, updates: usize) -> Computation {
    compute_observed(
        spec,
        &CancelGate::new_clock_free(),
        updates,
        None,
        0.0,
        |_, _| Ok(()),
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-8 * a.abs().max(b.abs()).max(1e-20),
        "{a:e} != {b:e}"
    );
}

fn report(ledger: &Ledger, receipt: &Loaded) -> JsonValue {
    let bytes = linked(ledger, &receipt.value, "report_json", "study-report-json").unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

#[test]
fn g1_affine_strength_source_binds_the_full_family_and_admits_only_valid_physics() {
    let spec = fixture();
    let config = spec.affine_load_space.unwrap();
    let strength = config.stress.unwrap();
    assert_eq!(strength.options.max_cases, spec.loads.len() + 1);
    assert_eq!(strength.options.max_points, 100_000);
    assert_eq!(config.constraints(None).stress, Some(strength));
    assert!(config.constraints(None).responses.is_none());

    for (from, to) in [
        (":stress-limit-pa 40.0", ":stress-limit-pa 0.0"),
        (":relaxation-power 1.0", ":relaxation-power 0.9"),
        (":relaxation-power 1.0", ":relaxation-power 3.0"),
        (":aggregation-power 8.0", ":aggregation-power 1.0"),
        (":aggregation-power 8.0", ":aggregation-power 65.0"),
        (
            ":max-stress-point-fields 100000",
            ":max-stress-point-fields 0",
        ),
        (
            ":max-stress-point-fields 100000",
            ":max-stress-point-fields 1000001",
        ),
        (":joint-sharpness 40.0", ":joint-sharpness 0.0"),
        (":density-floor 0.05", ":density-floor 0.001"),
        // The point-field allowance is part of admission, even before solving.
        (":memory-bytes 536870912", ":memory-bytes 402653184"),
    ] {
        assert!(
            spec::parse(&FIXTURE.replace(from, to)).is_err(),
            "admitted {to}"
        );
    }

    let start = FIXTURE.find("  (affine-strength").unwrap();
    let without_strength = format!("{})\n", &FIXTURE[..start]);
    assert!(
        spec::parse(&without_strength)
            .unwrap()
            .affine_load_space
            .unwrap()
            .stress
            .is_none()
    );
    let section = &FIXTURE[start..FIXTURE.trim_end().len() - 1];
    let duplicate = format!("{}{}{})\n", &FIXTURE[..start], section, section);
    assert!(spec::parse(&duplicate).is_err());
    let material_start = FIXTURE.find("  (material-acceleration").unwrap();
    let wrong_order = format!(
        "{}{}{})\n",
        &FIXTURE[..material_start],
        section,
        &FIXTURE[material_start..start]
    );
    assert!(spec::parse(&wrong_order).is_err());
}

#[test]
fn g1_g5_affine_strength_cap_controls_the_joint_constraint_and_feasible_archive() {
    let mut spec = fixture();
    let initial = compute(&spec, 0);
    let initial_point = &initial.state.accepted;
    let e = &initial_point.evaluation;
    let stress = e.stress.as_ref().unwrap();
    assert!(e.oracle.gradient_admissible());
    assert!(stress.gradient_admissible());
    assert!(
        stress
            .gradient
            .as_ref()
            .is_some_and(|g| g.len() == initial_point.rho.len())
    );
    assert!(stress.sampled_physical_max > 0.0);
    assert!(stress.aggregate >= stress.sampled_physical_max);
    assert!(stress.aggregate >= stress.sampled_relaxed_max);
    assert_eq!(stress.point_count, stress.witnesses.len());
    assert_eq!(
        stress.point_fields,
        stress.point_count * (spec.loads.len() + 1)
    );
    assert_eq!(stress.adjoints.len(), spec.loads.len() + 1);
    assert!(initial.state.gradient.passed());
    assert_eq!(initial.state.gradient.constraints, e.constraints);
    near(
        initial.state.gradient.baseline_constraint,
        initial_point.constraint(GOAL).unwrap(),
    );

    // Calibrate to actual physical stresses, keeping compliance independently
    // loose. This distinguishes enforcing strength from merely reporting it.
    let config = spec.affine_load_space.as_mut().unwrap();
    config.cap_j = 2.0 * e.worst_compliance;
    config.stress.as_mut().unwrap().limit = 2.0 * stress.aggregate;
    let loose = compute(&spec, 4);
    let best = loose
        .state
        .best
        .as_ref()
        .expect("the initial point is feasible");
    assert!(spec.affine_load_space.unwrap().feasible(best));
    assert!(
        best.evaluation
            .stress
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(best.evaluation.physical.dilated_volume_fraction < e.physical.dilated_volume_fraction);
    assert!(loose.state.gradient.passed());
    assert_eq!(
        loose.state.gradient.constraints,
        spec.affine_load_space.unwrap().constraint_policy()
    );

    let mut strict_spec = spec.clone();
    strict_spec
        .affine_load_space
        .as_mut()
        .unwrap()
        .stress
        .as_mut()
        .unwrap()
        .limit = 0.5 * stress.sampled_physical_max;
    let strict = compute(&strict_spec, 0);
    let point = &strict.state.accepted;
    let strength = point.evaluation.stress.as_ref().unwrap();
    near(strength.aggregate, stress.aggregate);
    near(
        strength.constraint,
        strength.aggregate / strength.limit.limit - 1.0,
    );
    assert!(point.evaluation.worst_compliance < strict_spec.affine_load_space.unwrap().cap_j);
    assert!(strength.sampled_physical_max > strength.limit.limit);
    assert!(!strength.physically_feasible());
    assert!(point.constraint(GOAL).unwrap() >= strength.constraint);
    assert!(point.constraint(GOAL).unwrap() > loose.state.gradient.baseline_constraint);
    assert!(strict.state.gradient.passed());
    assert!(!strict_spec.affine_load_space.unwrap().feasible(point));
    assert!(strict.state.best.is_none());
    assert_ne!(strict.state.status, "completed");
}

#[test]
fn g4_affine_strength_charges_all_point_fields_and_original_evaluation_allowance() {
    let mut spec = fixture();
    let baseline = compute(&spec, 0);
    let stress = baseline.state.accepted.evaluation.stress.as_ref().unwrap();
    let fields = stress.point_fields;
    assert!(fields > 1);

    spec.affine_load_space
        .as_mut()
        .unwrap()
        .stress
        .as_mut()
        .unwrap()
        .options
        .max_points = fields;
    let exact = compute(&spec, 0);
    let exact_stress = exact.state.accepted.evaluation.stress.as_ref().unwrap();
    assert_eq!(exact_stress.point_fields, fields);
    near(exact_stress.aggregate, stress.aggregate);

    spec.affine_load_space
        .as_mut()
        .unwrap()
        .stress
        .as_mut()
        .unwrap()
        .options
        .max_points = fields - 1;
    let refused = compute_observed(
        &spec,
        &CancelGate::new_clock_free(),
        0,
        None,
        0.0,
        |_, _| panic!("a point-field overrun cannot publish a baseline"),
    );
    let error = refused.err().unwrap();
    assert_eq!(error.exit, exit::REFUSED);
    // A point-field overrun is a budget stage of the affine stress physics (fs-topopt `controlled`), not a physics failure.
    assert!(error.message.contains("sdf3-affine-stress-point-budget"));

    let mut spec = fixture();
    spec.affine_load_space
        .as_mut()
        .unwrap()
        .options
        .optimizer
        .max_evaluations = 1;
    let bounded = compute(&spec, 8);
    assert_eq!(bounded.state.optimizer.evaluations, 1);
    assert_eq!(bounded.state.optimizer.iterations, 0);
    assert_eq!(bounded.state.status, "budget-exhausted");
    assert!(bounded.state.accepted.evaluation.stress.is_some());
    assert!(bounded.state.spent.linear.linear_solves > spec.loads.len());
    assert!(bounded.state.spent.linear.linear_iterations <= spec.linear);
}

#[test]
fn g4_affine_strength_cancellation_keeps_the_accepted_native_constraints() {
    let spec = fixture();
    let gate = CancelGate::new_clock_free();
    gate.request();
    let stopped = compute_observed(&spec, &gate, 1, None, 0.0, |_, _| {
        panic!("pre-cancelled strength study cannot publish a baseline")
    });
    assert_eq!(stopped.err().unwrap().exit, exit::CANCELLED);

    let gate = CancelGate::new_clock_free();
    let mut retained = None;
    let stopped = compute_observed(&spec, &gate, 4, None, 0.0, |_, state| {
        let stress = state.accepted.evaluation.stress.as_ref().unwrap();
        retained = Some((state.accepted.rho.clone(), stress.aggregate, stress.limit));
        gate.request();
        Ok(())
    })
    .unwrap();
    let (rho, aggregate, limit) = retained.unwrap();
    let stress = stopped.state.accepted.evaluation.stress.as_ref().unwrap();
    assert_eq!(stopped.state.status, "cancelled");
    assert_eq!(stopped.state.optimizer.iterations, 0);
    assert_eq!(stopped.state.accepted.rho, rho);
    assert_eq!(stress.aggregate, aggregate);
    assert_eq!(stress.limit, limit);
    assert_eq!(
        stopped.state.gradient.constraints,
        spec.affine_load_space.unwrap().constraint_policy()
    );
    assert!(stopped.state.spent.linear.linear_iterations > 0);
}

#[test]
fn g3_g5_affine_strength_resume_retains_policy_and_refuses_another_load_case() {
    let spec = fixture();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let prefix = drive(&spec, &ledger, Some(0), &gate).unwrap();
    let old = load(&ledger, &prefix.pointer).unwrap();
    let before = report(&ledger, &old);
    let resumed = resume(&ledger, &old, Some(1), &gate).unwrap();
    let resumed = load(&ledger, &resumed.pointer).unwrap();
    let after = report(&ledger, &resumed);
    let policy = before
        .get("model")
        .unwrap()
        .get("strength_constraint")
        .unwrap();
    assert_eq!(policy.str_field("kind"), Some("affine-sampled-peak-bound"));
    assert_eq!(
        Some(policy),
        after.get("model").unwrap().get("strength_constraint")
    );
    assert_eq!(before.get("gradient_check"), after.get("gradient_check"));
    assert_eq!(
        after
            .get("gradient_check")
            .unwrap()
            .get("constraint_policy"),
        Some(policy)
    );
    assert_eq!(
        after
            .get("selected")
            .unwrap()
            .get("strength")
            .unwrap()
            .get("stress_limit_pa"),
        policy.get("stress_limit_pa")
    );
    assert_eq!(integer(&resumed.value, "iterations_completed").unwrap(), 1);
    assert_eq!(
        integer(
            resumed.value.get("optimizer_work").unwrap(),
            "restoration_evaluations"
        )
        .unwrap(),
        1
    );
    let before_work = Spent::read(&old.value).unwrap();
    let after_work = Spent::read(&resumed.value).unwrap();
    assert!(after_work.linear.linear_solves > before_work.linear.linear_solves);

    // A new basis vector or strength policy is a different physical problem,
    // even when the caller presents the previous run's valid receipt.
    for source in [
        FIXTURE.replace("(0.0 -1.0 -0.5)", "(0.0 -1.0 -0.6)"),
        FIXTURE.replace(":relaxation-power 1.0", ":relaxation-power 1.5"),
        FIXTURE.replace(":joint-sharpness 40.0", ":joint-sharpness 50.0"),
    ] {
        let changed = spec::parse(&source).unwrap();
        assert_ne!(changed.id, spec.id);
        let error = resume::drive(&changed, &ledger, Some(1), &gate, Some(&old)).unwrap_err();
        assert_eq!(error.code, "cli-study-sdf3-affine-load-space-checkpoint");
        assert_eq!(error.exit, exit::REFUSED);
    }
    assert_eq!(load(&ledger, &prefix.pointer).unwrap().bytes, old.bytes);
}
