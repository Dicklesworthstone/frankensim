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

const REGIONS_FIXTURE: &str =
    include_str!("../../../../../../examples/marquee/bracket-3d-stress-regions.fsim");

#[test]
fn g0_stress_regions_are_identity_bearing_and_keep_free_density_admission() {
    let spec = spec::parse(REGIONS_FIXTURE).unwrap();
    assert_eq!(spec.regions.len(), 2);
    let canonical = spec::parse(&spec.canonical).unwrap();
    assert_eq!(canonical.id, spec.id);
    assert_eq!(canonical.regions, spec.regions);
    let changed = REGIONS_FIXTURE.replace(
        ":void (((0.0 0.5 0.5) (0.5 1.0 1.0)))",
        ":void (((0.5 0.5 0.5) (1.0 1.0 1.0)))",
    );
    assert_ne!(spec::parse(&changed).unwrap().id, spec.id);
    assert!(spec::parse(&REGIONS_FIXTURE.replace(
        ":density-floor 0.05", ":density-floor 0.0001",
    )).is_err());
}

fn assert_stress_regions(study: &CutDensityStudy3<AdaptiveSolveSpace3>, state: &State) {
    use fs_topopt::sdf3::PhysicalRegion3;
    let labels = study.physical_regions().expect("authored labels must be bound before any solve");
    let leaves = study.operator().elasticity().leaves();
    assert_eq!(labels.len(), leaves.len());
    assert_eq!(labels.iter().filter(|r| **r == PhysicalRegion3::Solid).count(), 1);
    assert_eq!(labels.iter().filter(|r| **r == PhysicalRegion3::Void).count(), 1);
    for field in [state.accepted.as_ref(), state.best.as_ref()].into_iter().flatten() {
        for ((leaf, label), density) in leaves.iter().zip(labels).zip(&field.projected_rho) {
            let expected = match leaf.index() {
                [0, 0, 0] => PhysicalRegion3::Solid,
                [0, 1, 1] => PhysicalRegion3::Void,
                _ => PhysicalRegion3::Design,
            };
            assert_eq!(*label, expected);
            match label {
                PhysicalRegion3::Solid => assert_eq!(*density, 1.0),
                PhysicalRegion3::Void => assert_eq!(*density, 0.0),
                PhysicalRegion3::Design => assert!(*density > 0.0 && *density <= 1.0),
            }
        }
        assert!(field.volume_fraction >= study.prescribed_solid_fraction());
    }
}

#[test]
fn g1_stress_regions_use_the_actual_gradient_gate_and_survive_accepted_material_removal() {
    let spec = spec::parse(REGIONS_FIXTURE).unwrap();
    let mut checkpoints = 0;
    let run = compute_observed(&spec, &CancelGate::new(), 4, None, 0.0, |study, state| {
        checkpoints += 1;
        assert!(state.audit.passed);
        assert_stress_regions(study, state);
        Ok(())
    }).unwrap();
    assert!(checkpoints > 1, "require durable accepted updates");
    assert_eq!(run.state.iterations(), 4);
    assert!(run.state.audit.passed);
    assert_eq!(run.state.audit.probes.len(), 2);
    assert_stress_regions(&run.study, &run.state);
    let best = run.state.best.as_ref().expect("a feasible design must be retained");
    assert!(best.volume_fraction < run.state.history[0].volume_fraction);
    assert!(best.aggregate <= spec.stress.unwrap().stress_limit * (1.0 + 2e-6));
    assert_ne!(best.displacements[0], best.displacements[1]);
}
