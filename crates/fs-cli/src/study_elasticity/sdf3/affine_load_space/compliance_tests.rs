use super::*;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-compliance.fsim"
));

fn fixed_source() -> String {
    let (source, _) = FIXTURE.split_once("  (affine-grid-assessment").unwrap();
    format!("{source})\n")
        .replace(":maximum-level 3", ":maximum-level 1")
        .replace(":maximum-leaves 512", ":maximum-leaves 8")
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

fn field<'a>(value: &'a JsonValue, key: &str) -> &'a JsonValue {
    value.get(key).unwrap_or_else(|| panic!("missing {key}"))
}

fn real(value: &JsonValue, key: &str) -> f64 {
    field(value, key).as_f64().unwrap()
}

fn document(bytes: &[u8]) -> JsonValue {
    JsonValue::parse(std::str::from_utf8(bytes).unwrap()).unwrap()
}

fn report(ledger: &Ledger, receipt: &Loaded) -> JsonValue {
    document(&linked(ledger, &receipt.value, "report_json", "study-report-json").unwrap())
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-8 * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_affine_compliance_minimizes_the_full_load_ball_under_the_actual_volume_cap() {
    let source = fixed_source();
    let mut spec = spec::parse(&source).unwrap();
    let goal = AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit: 0.9 };
    assert_eq!(spec.affine_load_space.unwrap().goal, goal);
    assert_eq!(spec.volume, 0.9);
    assert!(spec::parse(&source.replace(":initial-density 0.75", ":initial-density 0.95")).is_ok());
    for (from, to) in [
        (":volume-fraction 0.9", ":volume-fraction 0.0"),
        (":volume-fraction 0.9", ":volume-fraction 1.1"),
        (":compliance-reference-j 0.01", ":compliance-reference-j 0.0"),
        (":compliance-reference-j", ":compliance-limit-j"),
        (":type compliance :sense minimize :unit \"J\"", ":type compliance :sense minimize :unit \"1\""),
    ] {
        assert!(spec::parse(&source.replace(from, to)).is_err(), "admitted {to}");
    }

    let initial = compute(&spec, 0);
    let initial_point = &initial.state.accepted;
    assert!(spec.affine_load_space.unwrap().feasible(initial_point));
    let reference = 0.5 * initial_point.evaluation.worst_compliance;
    let config = spec.affine_load_space.as_mut().unwrap();
    config.cap_j = reference;
    // Keep two measurable objective steps within the volume headroom. This
    // is an integration check, not a claim of optimizer stationarity.
    config.options.optimizer.max_spectral_step = 0.001;
    let run = compute(&spec, 2);
    let config = spec.affine_load_space.unwrap();
    let best = run.state.best.as_ref().expect("retain the volume-feasible archive");
    let e = &best.evaluation;
    assert_eq!(run.state.gradient.goal, goal);
    assert!(run.state.gradient.passed());
    near(run.state.gradient.baseline_objective, 2.0);
    near(run.state.gradient.baseline_constraint, initial_point.constraint(goal).unwrap());
    assert!(run.state.optimizer.iterations > 0);
    assert!(e.worst_compliance < initial_point.evaluation.worst_compliance);
    assert!(e.maximum > 1.0, "reference scale must not become a compliance cap");
    assert!(e.worst_compliance > reference);
    assert!(e.physical.dilated_volume_fraction <= 0.9);
    assert!(e.responses.as_ref().unwrap().physically_feasible());
    assert!(e.stress.as_ref().unwrap().physically_feasible());
    assert!(best.constraint(goal).unwrap() <= 0.0);
    assert!(config.feasible(best));

    let measurement = document(output::measurement(best, config, &spec).unwrap().as_bytes());
    assert_eq!(field(&measurement, "compliance_limit_j"), &JsonValue::Null);
    assert_eq!(field(&measurement, "compliance_feasible"), &JsonValue::Null);
    assert_eq!(field(&measurement, "volume_feasible"), &JsonValue::Bool(true));
    near(real(&measurement, "compliance_reference_j"), reference);
    near(real(&measurement, "volume_limit_fraction"), 0.9);
    near(real(&measurement, "objective_value"), e.worst_compliance);
    near(real(&measurement, "normalized_objective"), e.maximum);
    assert_eq!(measurement.str_field("objective_unit"), Some("J"));

    // The same physical starting point becomes infeasible when only its
    // dilated-volume cap changes; a loose compliance scale cannot excuse it.
    let volume_limit = 0.5 * initial_point.evaluation.physical.dilated_volume_fraction;
    spec.volume = volume_limit;
    spec.affine_load_space.as_mut().unwrap().goal =
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit };
    let strict = compute(&spec, 0);
    assert!(strict.state.best.is_none());
    assert!(!spec.affine_load_space.unwrap().feasible(&strict.state.accepted));
    assert!(strict.state.accepted.constraint(spec.affine_load_space.unwrap().goal).unwrap() > 0.0);
}

#[test]
fn g3_g5_affine_compliance_resume_preserves_the_objective_archive_and_joint_fields() {
    let source = fixed_source().replace(":max-updates 80", ":max-updates 4");
    let spec = spec::parse(&source).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let full = drive(&spec, &ledger, Some(4), &gate).unwrap();
    let full = load(&ledger, &full.pointer).unwrap();
    let first = drive(&spec, &ledger, Some(1), &gate).unwrap();
    let first = load(&ledger, &first.pointer).unwrap();
    let resumed = resume(&ledger, &first, Some(3), &gate).unwrap();
    let resumed = load(&ledger, &resumed.pointer).unwrap();
    for (key, kind) in [("design", "study-design"), ("iterations", "study-iterations")] {
        assert_eq!(
            linked(&ledger, &resumed.value, key, kind).unwrap(),
            linked(&ledger, &full.value, key, kind).unwrap(),
            "{key}"
        );
    }
    let full_report = report(&ledger, &full);
    let first_report = report(&ledger, &first);
    let resumed_report = report(&ledger, &resumed);
    assert_eq!(resumed_report.str_field("objective"), Some("compliance"));
    assert_eq!(resumed_report.str_field("objective_unit"), Some("J"));
    assert_eq!(resumed_report.str_field("selected_design"), Some("least-compliance-feasible"));
    assert_eq!(field(&resumed_report, "goal").str_field("kind"), Some("worst-compliance"));
    for key in ["selected", "last_accepted", "gradient_check", "stop"] {
        assert_eq!(field(&resumed_report, key), field(&full_report, key), "{key}");
    }
    assert_eq!(field(&resumed_report, "gradient_check"), field(&first_report, "gradient_check"));
    let selected = field(&resumed_report, "selected");
    assert_eq!(field(selected, "volume_feasible"), &JsonValue::Bool(true));
    assert_eq!(field(selected, "compliance_limit_j"), &JsonValue::Null);
    assert_eq!(field(selected, "compliance_feasible"), &JsonValue::Null);
    near(real(selected, "objective_value"), real(selected, "worst_compliance_j"));
    let full_work = field(&full.value, "optimizer_work");
    let restored_work = field(&resumed.value, "optimizer_work");
    let extra = integer(restored_work, "restoration_evaluations").unwrap();
    assert!((1..=2).contains(&extra));
    assert_eq!(integer(restored_work, "evaluations").unwrap(), integer(full_work, "evaluations").unwrap() + extra);
    for key in ["iterations", "rejected_trials", "multiplier_updates"] {
        assert_eq!(field(restored_work, key), field(full_work, key), "{key}");
    }
    let before = Spent::read(&first.value).unwrap();
    let after = Spent::read(&resumed.value).unwrap();
    assert!(after.linear.linear_solves > before.linear.linear_solves);
    assert!(after.geometry.points > before.geometry.points);
    after.validate(&spec).unwrap();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    let changed = spec::parse(&source.replace(":volume-fraction 0.9", ":volume-fraction 0.8")).unwrap();
    assert_ne!(changed.id, spec.id);
    assert!(resume::recover(&changed, &ledger, &first, producer).is_err());
    let finished = resume(&ledger, &resumed, Some(1), &gate).unwrap();
    assert_eq!(finished.receipt, resumed.bytes);
}

#[test]
fn g1_g3_affine_compliance_assessment_and_refit_keep_the_volume_goal_on_each_grid() {
    let mut spec = spec::parse(FIXTURE).unwrap();
    let preliminary = compute(&spec, 0);
    let probe = preliminary.state.assessment.as_ref().unwrap().report.as_ref().unwrap();
    spec.affine_load_space.as_mut().unwrap().cap_j =
        0.25 * probe.coarse.worst_compliance.min(probe.fine.worst_compliance);
    spec.affine_load_space.as_mut().unwrap().options.optimizer.max_spectral_step = 0.001;
    let source = compute(&spec, 0);
    let config = spec.affine_load_space.unwrap();
    let assessment = source.state.assessment.as_ref().unwrap();
    assert_eq!(assessment.status, "checked", "{:?}", assessment.reason);
    let probe = assessment.report.as_ref().unwrap();
    assert_eq!(probe.goal, Some(config.goal));
    assert!(probe.coarse.maximum > 1.0 && probe.fine.maximum > 1.0);
    assert!(probe.responses.as_ref().unwrap().coarse.physically_feasible());
    assert!(probe.responses.as_ref().unwrap().fine.physically_feasible());
    let inventory_change = probe.inventory.as_ref().unwrap().volume_change;
    let objective_change = probe.objective_change().unwrap();
    assert!(objective_change > inventory_change);
    let policy = grid_assessment::Policy {
        compliance_tolerance: 1e6,
        response_tolerance: 1e6,
        stress_tolerance: 1e6,
        inventory_tolerance: 0.5 * (inventory_change + objective_change),
        ..spec.affine_grid_assessment.unwrap()
    };
    // This gate lies between the two measured changes: volume stability passes
    // while compliance change must use its own declared tolerance.
    assert!(assessment.passed(policy, config.goal));
    assert!(!assessment.passed(policy, GOAL));
    let gate = CancelGate::new_clock_free();
    let context = grid_refinement::Context::new(&source, hash_bytes(b"in-memory source")).unwrap();
    let preparation = grid_refinement::prepare(&spec, &source, context, &gate);
    let prepared = preparation.result.unwrap();
    assert!(prepared.study.cells() > source.study.cells());
    let candidate_spec = grid_refinement::effective_spec(&spec, Some(&prepared.context)).unwrap();
    let candidate = compute_on_grid(
        &candidate_spec,
        &gate,
        1,
        None,
        0.0,
        Some(prepared),
        |_, _| Ok(()),
    ).unwrap();
    assert_eq!(candidate.state.gradient.goal, config.goal);
    assert!(candidate.state.gradient.passed());
    near(candidate.state.history[0].objective, candidate.state.gradient.baseline_objective);
    assert!(candidate.state.history[0].objective > 1.0);
    let best = candidate.state.best.as_ref().expect("fine volume-feasible archive");
    assert!(candidate_spec.affine_load_space.unwrap().feasible(best));
    assert!(best.evaluation.physical.dilated_volume_fraction <= spec.volume);
    assert!(best.evaluation.maximum > 1.0);
    assert!(best.evaluation.responses.as_ref().unwrap().physically_feasible());
    assert!(best.evaluation.stress.as_ref().unwrap().physically_feasible());
    assert!(candidate.state.spent.linear.linear_iterations > source.state.spent.linear.linear_iterations);
}
