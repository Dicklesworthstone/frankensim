use super::*;

const FIXTURE: &str = include_str!("../../../../../../examples/marquee/bracket-3d-stress.fsim");

#[test]
fn g0_stress_schema_has_explicit_units_work_and_monotone_material_branch() {
    let spec = spec::parse(FIXTURE).unwrap();
    assert!(spec.stress.is_some());
    assert_eq!(spec::parse(&spec.canonical).unwrap().id, spec.id);
    for source in [
        FIXTURE.replace(":stress-limit-pa 8.0", ":stress-limit-pa 0.0"),
        FIXTURE.replace(":stress-limit-pa 8.0", ":stress-limit-pa NaN"),
        FIXTURE.replace(":density-floor 0.05", ":density-floor 0.0001"),
        FIXTURE.replace(":relaxation-power 1.0", ":relaxation-power 3.0"),
        FIXTURE.replace(":aggregation-power 8.0", ":aggregation-power 1.0"),
        FIXTURE.replace(":maximum-level 1", ":maximum-level 2"),
        FIXTURE.replace(":max-stress-points 50000", ":max-stress-points 1000000000"),
        FIXTURE.replace(":max-evaluations 2000", ":max-evaluations 0"),
        FIXTURE.replace(":unit \"1\"", ":unit \"J\""),
        FIXTURE.replace(":type volume-fraction", ":type compliance"),
        FIXTURE.replace(":beta 0.0", ":beta 16.0"),
    ] {
        assert!(
            spec::parse(&source).is_err(),
            "must refuse unsupported stress declaration: {source}"
        );
    }
    let changed =
        spec::parse(&FIXTURE.replace(":stress-limit-pa 8.0", ":stress-limit-pa 9.0")).unwrap();
    assert_ne!(changed.id, spec.id);
    assert!(ir_for(STRESS3_DRIVER, spec.id, 0).contains("\"objective\":\"1\""));
}

#[test]
fn g1_stress_driver_uses_actual_gradients_and_retains_feasible_material_reduction() {
    let spec = spec::parse(FIXTURE).unwrap();
    let mut checkpoints = Vec::new();
    let run = compute_observed(&spec, &CancelGate::new(), 4, None, 0.0, |_, state| {
        checkpoints.push(state.iterations());
        assert!(state.audit.passed);
        assert!(state.accepted.is_some());
        Ok(())
    })
    .unwrap();
    assert_eq!(checkpoints, vec![0, 1, 2, 3]);
    assert_eq!(run.state.iterations(), 4);
    assert_eq!(run.state.status, "budget-exhausted");
    assert_eq!(run.state.audit.probes.len(), 2);
    let selected = run.state.selected().unwrap();
    assert!(selected.volume_fraction < run.state.history[0].volume_fraction);
    assert!(selected.aggregate <= spec.stress.unwrap().stress_limit * (1.0 + 2e-6));
    assert!(selected.aggregate > 0.0);
    assert!(selected.sampled_relaxed_max >= selected.aggregate);
    assert_ne!(selected.displacements[0], selected.displacements[1]);
    assert!(run.state.spent.linear.linear_iterations > 0);
    assert!(
        run.state.audit.evaluations + run.state.optimizer_work.evaluations
            <= spec.stress.unwrap().optimizer.max_evaluations
    );
}

#[test]
fn g4_stress_evaluation_and_linear_caps_publish_no_unchecked_design() {
    for source in [
        FIXTURE.replace(":max-evaluations 2000", ":max-evaluations 1"),
        FIXTURE.replace(":linear-iterations 250000", ":linear-iterations 1"),
        FIXTURE.replace(":max-stress-points 50000", ":max-stress-points 1"),
    ] {
        let spec = spec::parse(&source).unwrap();
        let run = compute_observed(&spec, &CancelGate::new(), 1, None, 0.0, |_, _| {
            panic!("no accepted design before gradient admission")
        })
        .unwrap();
        assert_eq!(run.state.status, "budget-exhausted");
        assert!(run.state.selected().is_none());
        assert!(!run.state.audit.passed);
        assert!(run.state.spent.linear.linear_iterations <= spec.linear);
        assert!(run.state.audit.evaluations <= spec.stress.unwrap().optimizer.max_evaluations);
    }
}

#[test]
fn g4_infeasible_stress_endpoint_is_never_completed_or_selected_as_feasible() {
    let spec =
        spec::parse(&FIXTURE.replace(":stress-limit-pa 8.0", ":stress-limit-pa 0.0001")).unwrap();
    let run = compute_observed(&spec, &CancelGate::new(), 1, None, 0.0, |_, _| Ok(())).unwrap();
    assert_ne!(run.state.status, "completed");
    assert!(run.state.best.is_none());
    assert!(run.state.accepted.as_ref().unwrap().aggregate > spec.stress.unwrap().stress_limit);
}
