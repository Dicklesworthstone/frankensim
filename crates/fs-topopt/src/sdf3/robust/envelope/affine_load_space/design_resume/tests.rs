use super::*;
use crate::sdf3::response::limits::ResponseLimit3;
use crate::sdf3::robust::envelope::affine_load_space::response::AffineResponseOptions3;
use crate::sdf3::robust::envelope::affine_load_space::stress::AffineStressLimit3;
use crate::sdf3::robust::envelope::tests::native_body_operator;
use crate::{SolveBudget, SolveProgress};
use fs_solver::LinearOp;

fn fixture(acceleration: f64, radius: f64) -> (CutDensityStudy3, Vec<f64>, Vec<f64>, Vec<f64>) {
    let op = native_body_operator();
    let first = op
        .body_load(&|p| [0.0, -0.001 * (1.0 + p[1]), 0.0], || {
            ControlFlow::Continue(())
        })
        .unwrap();
    let second = op
        .body_load(
            &|p| [0.0, -0.0006 * (1.0 + p[1]), -0.0008 * (1.0 + p[0])],
            || ControlFlow::Continue(()),
        )
        .unwrap();
    let dead = op
        .body_load(
            &|_| [0.0002, -0.0001, -0.0003],
            || ControlFlow::Continue(()),
        )
        .unwrap();
    let study = CutDensityStudy3::new(
        op,
        radius,
        SimpParams {
            beta: 2.0,
            penal: 3.0,
            e_min: 0.03,
            ..Default::default()
        },
    )
    .with_material_acceleration(
        move |p| [0.0, -acceleration * (1.0 + p[0]), -1.25 * acceleration],
        100_000,
        || ControlFlow::Continue(()),
    )
    .unwrap();
    let mut regions = vec![crate::sdf3::PhysicalRegion3::Design; study.cells()];
    regions[0] = crate::sdf3::PhysicalRegion3::Solid;
    let last = regions.len() - 1;
    regions[last] = crate::sdf3::PhysicalRegion3::Void;
    (
        study.with_physical_regions(regions).unwrap(),
        dead,
        first,
        second,
    )
}

fn options() -> AffineLoadSpaceDesignOptions3 {
    AffineLoadSpaceDesignOptions3 {
        density_floor: 0.3,
        ..Default::default()
    }
}
fn volume_goal() -> AffineLoadSpaceDesignGoal3 {
    AffineLoadSpaceDesignGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::VolumeFraction,
        volume_scale: 1.0,
    }
}
fn cap(
    study: &mut CutDensityStudy3,
    dead: &[f64],
    basis: &[&[f64]],
    rho: &[f64],
    factor: f64,
    control: &mut SolveControl<'_>,
) -> AffineLoadSpace3 {
    let initial = study
        .evaluate_robust_affine_load_space(
            rho,
            dead,
            basis,
            RobustOptions3::default(),
            &AffineLoadSpace3::new(2, 1.0).unwrap(),
            control,
        )
        .unwrap();
    AffineLoadSpace3::new(2, factor * initial.worst_compliance).unwrap()
}

#[test]
fn affine_restore_matches_uninterrupted_volume_mass_and_compliance_searches() {
    for goal in [
        volume_goal(),
        AffineLoadSpaceDesignGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
            volume_scale: 1.0,
        },
        AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit: 0.9 },
    ] {
        let run = |restore: bool| {
            let (mut study, dead, first, second) = fixture(0.0003, 0.15);
            let basis = [&first[..], &second[..]];
            let rho = vec![0.7; study.cells()];
            let mut poll = |_| ControlFlow::Continue(());
            let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
            let policy = cap(&mut study, &dead, &basis, &rho, 1.5, &mut control);
            let mut design = AffineLoadSpaceDesign3::new(
                &mut study,
                &dead,
                &basis,
                &rho,
                policy.clone(),
                goal,
                options(),
                &mut control,
            )
            .unwrap();
            design.run(2).unwrap();
            let mut checkpoint = design.checkpoint();
            let cost = if restore {
                checkpoint.restoration_cost()
            } else {
                0
            };
            let prior = design.work();
            if restore {
                drop(design);
                let fresh = fixture(0.0003, 0.15);
                assert_eq!((&dead, &first, &second), (&fresh.1, &fresh.2, &fresh.3));
                study = fresh.0;
                design = AffineLoadSpaceDesign3::restore(
                    &mut study,
                    &dead,
                    &basis,
                    &mut checkpoint,
                    policy,
                    goal,
                    options(),
                    &mut control,
                )
                .unwrap();
                assert!(design.work().linear_iterations > prior.linear_iterations);
            }
            design.run(4).unwrap();
            (
                design.checkpoint(),
                design.accepted().evaluation.worst_displacement.clone(),
                design
                    .best_feasible()
                    .map(|point| point.evaluation.worst_displacement.clone()),
                cost,
            )
        };
        let original = run(false);
        let restored = run(true);
        assert_eq!(original.0.history, restored.0.history);
        assert_eq!(original.0.optimizer.point, restored.0.optimizer.point);
        assert_eq!(
            original.0.best_feasible_density,
            restored.0.best_feasible_density
        );
        assert_eq!(original.1, restored.1);
        assert_eq!(original.2, restored.2);
        let mut expected = original.0.optimizer.clone();
        expected.work.evaluations += restored.3;
        assert_eq!(expected, restored.0.optimizer);
        assert_eq!(restored.0.restoration_evaluations, restored.3);
        assert_eq!(original.0.binding, restored.0.binding);
        assert_eq!(
            original.0.gradient_check.work,
            restored.0.gradient_check.work
        );
        assert_eq!(
            original.0.gradient_check.objective_relative_error,
            restored.0.gradient_check.objective_relative_error
        );
        assert!(restored.0.history.len() > 2);
    }
}

#[test]
fn affine_restore_refuses_changed_loads_body_filter_policy_and_corrupt_state() {
    let (mut study, dead, first, second) = fixture(0.0003, 0.15);
    let basis = [&first[..], &second[..]];
    let rho = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = cap(&mut study, &dead, &basis, &rho, 1.5, &mut control);
    let checkpoint = {
        let mut design = AffineLoadSpaceDesign3::new(
            &mut study,
            &dead,
            &basis,
            &rho,
            policy.clone(),
            volume_goal(),
            options(),
            &mut control,
        )
        .unwrap();
        design.run(2).unwrap();
        design.checkpoint()
    };
    let scales = study.operator().scales().to_vec();
    let before = control.work();
    // These changed entries are strongly clamped and physically masked out.
    // Their source identity must still bind, rather than relying on one response.
    assert!(study.operator().fixed()[0]);
    let mut changed_dead = dead.clone();
    changed_dead[0] += 0.1;
    let mut changed_first = first.clone();
    changed_first[0] += 0.1;
    for (d, f) in [(&changed_dead, &first), (&dead, &changed_first)] {
        let mut copy = checkpoint.clone();
        let basis = [&f[..], &second[..]];
        assert!(
            AffineLoadSpaceDesign3::restore(
                &mut study,
                d,
                &basis,
                &mut copy,
                policy.clone(),
                volume_goal(),
                options(),
                &mut control
            )
            .is_err()
        );
        assert_eq!(copy.optimizer, checkpoint.optimizer);
    }
    for (acceleration, radius) in [(0.0004, 0.15), (0.0003, 0.16)] {
        let mut fresh = fixture(acceleration, radius).0;
        let mut copy = checkpoint.clone();
        assert!(
            AffineLoadSpaceDesign3::restore(
                &mut fresh,
                &dead,
                &basis,
                &mut copy,
                policy.clone(),
                volume_goal(),
                options(),
                &mut control
            )
            .is_err()
        );
        assert_eq!(copy.optimizer, checkpoint.optimizer);
    }
    let mut changed = options();
    changed.optimizer.max_evaluations += 1;
    assert!(
        AffineLoadSpaceDesign3::restore(
            &mut study,
            &dead,
            &basis,
            &mut checkpoint.clone(),
            policy.clone(),
            volume_goal(),
            changed,
            &mut control
        )
        .is_err()
    );
    let mut history = checkpoint.clone();
    history.history.pop();
    let mut lost = checkpoint.clone();
    assert!(lost.best_feasible_density.take().is_some());
    let mut audit = checkpoint.clone();
    audit.gradient_check.baseline_objective += 0.1;
    for mut copy in [history, lost, audit] {
        let work = copy.optimizer.work;
        assert!(
            AffineLoadSpaceDesign3::restore(
                &mut study,
                &dead,
                &basis,
                &mut copy,
                policy.clone(),
                volume_goal(),
                options(),
                &mut control
            )
            .is_err()
        );
        assert_eq!(copy.optimizer.work, work);
    }
    assert_eq!(
        control.work(),
        before,
        "structural refusals must spend no endpoint physics"
    );
    let mut derivative = checkpoint.clone();
    derivative.optimizer.sample.constraint_gradient[0] += 0.1;
    assert!(
        AffineLoadSpaceDesign3::restore(
            &mut study,
            &dead,
            &basis,
            &mut derivative,
            policy,
            volume_goal(),
            options(),
            &mut control
        )
        .is_err()
    );
    assert!(control.work().linear_iterations > before.linear_iterations);
    assert_eq!(
        derivative.restoration_evaluations,
        checkpoint.restoration_cost()
    );
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn affine_restore_preserves_reserved_work_on_cancellation_and_linear_budget_stops() {
    let (mut study, dead, first, second) = fixture(0.0003, 0.15);
    let basis = [&first[..], &second[..]];
    let rho = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = cap(&mut study, &dead, &basis, &rho, 1.5, &mut control);
    let checkpoint = AffineLoadSpaceDesign3::new(
        &mut study,
        &dead,
        &basis,
        &rho,
        policy.clone(),
        volume_goal(),
        options(),
        &mut control,
    )
    .unwrap()
    .checkpoint();
    let scales = study.operator().scales().to_vec();
    let cost = checkpoint.restoration_cost();
    let mut copy = checkpoint.clone();
    let mut cancel = |progress: SolveProgress| {
        if progress.stage == "sdf3-affine-restore-reserved" {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let mut cancelled = SolveControl::new(SolveBudget::default(), &mut cancel);
    assert!(matches!(
        AffineLoadSpaceDesign3::restore(
            &mut study,
            &dead,
            &basis,
            &mut copy,
            policy.clone(),
            volume_goal(),
            options(),
            &mut cancelled
        ),
        Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))
    ));
    assert_eq!(copy.restoration_evaluations, cost);
    assert_eq!(
        copy.optimizer.work.evaluations,
        checkpoint.optimizer.work.evaluations + cost
    );
    assert_eq!(copy.optimizer.sample, checkpoint.optimizer.sample);
    assert_eq!(copy.history, checkpoint.history);
    assert_eq!(cancelled.work(), SolveWork::default());
    // `control` still holds `poll` and is used below; the one-iteration
    // budget gets its own callback with the same Continue behaviour.
    let mut tiny_poll = |_| ControlFlow::Continue(());
    let mut tiny = SolveControl::new(
        SolveBudget {
            total_iterations: 1,
            ..SolveBudget::default()
        },
        &mut tiny_poll,
    );
    assert!(matches!(
        AffineLoadSpaceDesign3::restore(
            &mut study,
            &dead,
            &basis,
            &mut copy,
            policy.clone(),
            volume_goal(),
            options(),
            &mut tiny
        ),
        Err(ProjectedAlError::Evaluation(
            EvaluationStop::TotalBudget { .. }
        ))
    ));
    assert_eq!(tiny.work().linear_iterations, 1);
    assert_eq!(copy.restoration_evaluations, 2 * cost);
    assert_eq!(study.operator().scales(), scales);

    let mut bounded = options();
    bounded.optimizer.max_evaluations = 1;
    let mut exhausted = AffineLoadSpaceDesign3::new(
        &mut study,
        &dead,
        &basis,
        &rho,
        policy.clone(),
        volume_goal(),
        bounded,
        &mut control,
    )
    .unwrap()
    .checkpoint();
    let original = exhausted.optimizer.clone();
    let before = control.work();
    assert!(
        AffineLoadSpaceDesign3::restore(
            &mut study,
            &dead,
            &basis,
            &mut exhausted,
            policy,
            volume_goal(),
            bounded,
            &mut control
        )
        .is_err()
    );
    assert_eq!(exhausted.optimizer, original);
    assert_eq!(exhausted.restoration_evaluations, 0);
    assert_eq!(control.work(), before);
}

#[test]
fn affine_restore_keeps_joint_response_and_sampled_strength_requirements() {
    let (mut study, dead, first, second) = fixture(0.0003, 0.15);
    let basis = [&first[..], &second[..]];
    let observation: Vec<_> = first.iter().map(|f| -f).collect();
    let limits = [ResponseLimit3 {
        response: BoundedResponse3::Displacement(&observation),
        sense: ResponseLimitSense3::AtLeast,
        bound: -1.0,
        scale: 1.0,
    }];
    let constraints = AffineLoadSpaceConstraints3 {
        responses: Some(AffineResponseLimits3 {
            limits: &limits,
            options: AffineResponseOptions3::default(),
        }),
        stress: Some(AffineStressLimit3::default()),
        ..Default::default()
    };
    let options = AffineLoadSpaceDesignOptions3 {
        density_floor: 0.6,
        ..options()
    };
    let rho = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = cap(&mut study, &dead, &basis, &rho, 2.0, &mut control);
    let goal = AffineLoadSpaceDesignGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::MassKg { scale_kg: 2.0 },
        volume_scale: 1.0,
    };
    let (mut checkpoint, expected, worst) = {
        let mut design = AffineLoadSpaceDesign3::new_with_constraints(
            &mut study,
            &dead,
            &basis,
            &rho,
            policy.clone(),
            goal,
            options,
            constraints,
            &mut control,
        )
        .unwrap();
        design.run(1).unwrap();
        let checkpoint = design.checkpoint();
        design.run(1).unwrap();
        (
            checkpoint,
            design.checkpoint(),
            design.accepted().evaluation.worst_displacement.clone(),
        )
    };
    let mut fresh = fixture(0.0003, 0.15).0;
    let before = control.work();
    assert!(
        AffineLoadSpaceDesign3::restore(
            &mut fresh,
            &dead,
            &basis,
            &mut checkpoint.clone(),
            policy.clone(),
            goal,
            options,
            &mut control
        )
        .is_err()
    );
    let changed_limits = [ResponseLimit3 {
        bound: -2.0,
        ..limits[0]
    }];
    let changed = AffineLoadSpaceConstraints3 {
        responses: Some(AffineResponseLimits3 {
            limits: &changed_limits,
            options: AffineResponseOptions3::default(),
        }),
        ..constraints
    };
    assert!(
        AffineLoadSpaceDesign3::restore_with_constraints(
            &mut fresh,
            &dead,
            &basis,
            &mut checkpoint.clone(),
            policy.clone(),
            goal,
            options,
            changed,
            &mut control
        )
        .is_err()
    );
    assert_eq!(control.work(), before);
    let mut restored = AffineLoadSpaceDesign3::restore_with_constraints(
        &mut fresh,
        &dead,
        &basis,
        &mut checkpoint,
        policy,
        goal,
        options,
        constraints,
        &mut control,
    )
    .unwrap();
    restored.run(1).unwrap();
    assert_eq!(restored.history(), expected.history);
    assert_eq!(restored.accepted().evaluation.worst_displacement, worst);
    let evaluation = &restored.accepted().evaluation;
    assert_eq!(evaluation.constraints, constraints.policy());
    assert!(evaluation.responses.as_ref().unwrap().physically_feasible());
    assert!(evaluation.stress.as_ref().unwrap().physically_feasible());
}
