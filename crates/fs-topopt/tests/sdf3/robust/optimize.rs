use super::*;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlStop};
use fs_topopt::sdf3::robust::{RobustDesignOptions3, RobustDesignStudy3};

fn policy(scale: f64) -> RobustDesignOptions3 {
    RobustDesignOptions3 { compliance_scale: scale, ..Default::default() }
}

#[test]
fn zero_step_reports_eroded_objective_and_dilated_constraint_without_solving_again() {
    let (mut study, y, z) = fixture(SimpParams { beta: 4.0, ..Default::default() });
    let rho = vec![0.5; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    {
        let mut session = RobustDesignStudy3::new(&mut study, &loads, &rho, policy(17.0), &mut control).unwrap();
        let before = session.work(); let optimizer_before = session.optimizer_work();
        let report = session.run(0).unwrap();
        let e = &session.accepted().evaluation;
        assert_eq!(report.objective, e.eroded.objective.compliance / 17.0);
        assert_eq!(report.constraint, e.dilated_volume_fraction / 0.5 - 1.0);
        assert!(e.nominal_volume_fraction < e.dilated_volume_fraction);
        assert!(!session.feasible(), "nominal feasibility must not substitute for dilated feasibility");
        assert!(session.best_feasible().is_none());
        assert_eq!(session.point(), rho); assert_eq!(session.accepted().rho, rho);
        assert_eq!(session.work(), before); assert_eq!(session.optimizer_work(), optimizer_before);
        assert_eq!(session.history().len(), 1);
    }
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn real_eroded_compliance_design_improves_and_retains_a_feasible_complete_incumbent() {
    let (mut study, y, z) = fixture(SimpParams { beta: 3.0, ..Default::default() });
    let rho = vec![0.3; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let baseline = study.evaluate_robust(&rho, &loads, RobustOptions3::default(), &mut control).unwrap();
    let scales = study.operator().scales().to_vec();
    let best = {
        let mut session = RobustDesignStudy3::new(&mut study, &loads, &rho,
            policy(baseline.eroded.objective.compliance), &mut control).unwrap();
        session.run(100).unwrap();
        assert!(session.history().len() > 1);
        let best = session.best_feasible().unwrap().clone();
        assert!(best.evaluation.eroded.objective.compliance < baseline.eroded.objective.compliance);
        assert!(best.evaluation.dilated_volume_fraction <= 0.5*(1.0 + session.options().optimizer.tolerance));
        for row in session.history() {
            assert!(row.eroded_volume_fraction <= row.nominal_volume_fraction);
            assert!(row.nominal_volume_fraction <= row.dilated_volume_fraction);
        }
        best
    };
    assert_eq!(study.operator().scales(), scales);
    let solved = study.evaluate_robust(&best.rho, &loads, RobustOptions3::default(), &mut control).unwrap();
    assert_eq!(solved.eroded.objective.compliance.to_bits(), best.evaluation.eroded.objective.compliance.to_bits());
    assert_eq!(solved.eroded.objective.displacements, best.evaluation.eroded.objective.displacements);
    assert_eq!(solved.dilated_rho, best.evaluation.dilated_rho);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn segmented_runs_preserve_exact_fields_dual_state_and_cumulative_work() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho = vec![0.35; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let run = |study: &mut CutDensityStudy3, segmented: bool| {
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let mut session = RobustDesignStudy3::new(study, &loads, &rho, policy(100.0), &mut control).unwrap();
        let report = if segmented { session.run(4).unwrap(); session.run(6).unwrap() } else { session.run(10).unwrap() };
        (session.accepted().clone(), session.history().to_vec(), session.work(), report)
    };
    let (a, ah, aw, ar) = run(&mut study, false);
    let (b, bh, bw, br) = run(&mut study, true);
    assert_eq!(a.rho, b.rho); assert_eq!(a.evaluation.eroded.objective.displacements, b.evaluation.eroded.objective.displacements);
    assert_eq!(ah, bh); assert_eq!(aw, bw); assert_eq!(ar.work, br.work);
    assert_eq!(ar.multiplier.to_bits(), br.multiplier.to_bits()); assert_eq!(ar.penalty.to_bits(), br.penalty.to_bits());
}

#[test]
fn interrupted_trial_preserves_accepted_fields_and_counts_spent_work() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho = vec![0.35; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let scales = study.operator().scales().to_vec();
    let mut evaluations = 0;
    let mut poll = |p: SolveProgress| {
        if p.stage == "sdf3-robust-evaluation" { evaluations += 1; }
        if evaluations > 1 && p.stage == "sdf3-elasticity" && p.solve_iterations > 0 { ControlFlow::Break(()) }
        else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    {
        let mut session = RobustDesignStudy3::new(&mut study, &loads, &rho, policy(100.0), &mut control).unwrap();
        let before = session.work(); let fields = session.accepted().evaluation.eroded.objective.displacements.clone();
        assert!(session.run(10).is_err());
        assert_eq!(session.point(), rho); assert_eq!(session.accepted().evaluation.eroded.objective.displacements, fields);
        assert_eq!(session.history().len(), 1); assert_eq!(session.optimizer_work().iterations, 0);
        assert_eq!(session.optimizer_work().evaluations, 2);
        assert!(session.work().linear_iterations > before.linear_iterations);
        assert_eq!(session.best_feasible().unwrap().rho, rho);
    }
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn stop_after_acceptance_keeps_matching_fields_and_can_resume_without_restarting() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho = vec![0.35; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut armed = true;
    let mut poll = |p: SolveProgress| {
        if p.stage == "sdf3-robust-design-accepted" && armed { armed = false; ControlFlow::Break(()) }
        else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut session = RobustDesignStudy3::new(&mut study, &loads, &rho, policy(100.0), &mut control).unwrap();
    assert!(matches!(session.run(10), Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(session.optimizer_work().iterations, 1); assert_eq!(session.history().len(), 2);
    assert_eq!(session.accepted().rho, session.point()); assert_ne!(session.point(), rho);
    let before = session.work();
    session.run(0).unwrap(); assert_eq!(session.work(), before);
    session.run(2).unwrap();
    assert!(session.optimizer_work().iterations > 1);
    assert_eq!(session.accepted().rho, session.point());
}

#[test]
fn evaluation_allowance_is_not_reset_and_impossible_protected_volume_costs_no_solve() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho = vec![0.35; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    {
        let mut options = policy(100.0); options.optimizer.max_evaluations = 1;
        let mut session = RobustDesignStudy3::new(&mut study, &loads, &rho, options, &mut control).unwrap();
        let before = session.work();
        assert_eq!(session.run(20).unwrap().stop, ProjectedAlStop::EvaluationLimit);
        assert_eq!(session.run(20).unwrap().stop, ProjectedAlStop::EvaluationLimit);
        assert_eq!(session.optimizer_work().evaluations, 1); assert_eq!(session.work(), before);
    }
    let n = study.cells();
    let mut study = study.with_physical_regions(vec![PhysicalRegion3::Solid; n]).unwrap();
    let before = control.work();
    assert!(matches!(RobustDesignStudy3::new(&mut study, &loads, &rho, policy(100.0), &mut control),
        Err(ProjectedAlError::Invalid(_))));
    assert_eq!(control.work(), before);
}
