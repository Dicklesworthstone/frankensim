use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{SolveBudget, SolveProgress};
fn normalized_policy(study: &mut CutDensityStudy3, raw: &[f64], basis: &[&[f64]], control: &mut SolveControl<'_>) -> LoadSpace3 {
    let reference = study.evaluate_robust_load_space(raw, basis, RobustOptions3::default(),
        &LoadSpace3::new(basis.len(), 1.0, 0.05).unwrap(), control).unwrap().worst_compliance;
    LoadSpace3::new(basis.len(), reference, 0.05).unwrap()
}

#[test]
fn real_design_recovery_improves_load_space_without_rewriting_optimizer_state() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.3; study.cells()];
    let basis = [&y[..], &z[..]]; let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = normalized_policy(&mut study, &raw, &basis, &mut control);
    let options = EnvelopeDesignOptions3::default();
    let mut session = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy.clone(), options, &mut control).unwrap();
    let initial = session.accepted().evaluation.spectrum.smooth_maximum;
    assert!(session.gradient_check().numerical.passed());
    let outcome = session.run(30).unwrap();
    assert!(session.history().len() > 1);
    let before = (session.history().to_vec(), session.accepted().rho.clone(), session.optimizer_work());
    let work = session.work();
    let recovered = session.recover_feasible().unwrap();
    assert_eq!(before, (session.history().to_vec(), session.accepted().rho.clone(), session.optimizer_work()));
    assert!(session.work().linear_solves > work.linear_solves);
    assert!(recovered.evaluation.physical.dilated_volume_fraction <= options.volume_limit);
    assert!(recovered.evaluation.spectrum.smooth_maximum < initial);
    assert_eq!(outcome.objective, session.accepted().evaluation.spectrum.smooth_maximum);
    drop(session);
    let fresh = study.evaluate_robust_load_space(&recovered.rho, &basis, options.robustness, &policy, &mut control).unwrap();
    assert_eq!(fresh.worst_displacement, recovered.evaluation.worst_displacement);
    assert_eq!(fresh.gradient, recovered.evaluation.gradient);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn segmented_optimization_replays_trajectory_and_physical_work() {
    let run = |segments: &[usize]| {
        let (mut study, y, z) = fixture(1.0); let raw = vec![0.3; study.cells()];
        let basis = [&y[..], &z[..]];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let policy = normalized_policy(&mut study, &raw, &basis, &mut control);
        let mut session = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy,
            EnvelopeDesignOptions3::default(), &mut control).unwrap();
        let work = session.work(); session.run(0).unwrap(); assert_eq!(session.work(), work);
        for &steps in segments {
            if session.run(steps).unwrap().stop != ProjectedAlStop::IterationLimit { break; }
        }
        (session.history().to_vec(), session.accepted().rho.clone(), session.work(), session.optimizer_work())
    };
    let one = run(&[12]); let split = run(&[5,7]);
    assert!(one.0.len() > 1); assert_eq!(one, split);
}

#[test]
fn late_cancel_retains_the_accepted_physics_and_budget_refusals_do_not_publish_trials() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.3; study.cells()];
    let basis = [&y[..], &z[..]]; let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-design-accepted" { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = normalized_policy(&mut study, &raw, &basis, &mut control);
    let mut session = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy.clone(),
        EnvelopeDesignOptions3::default(), &mut control).unwrap();
    assert!(matches!(session.run(4), Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(session.history().len(), 2);
    assert_eq!(session.accepted().rho, session.state.point());
    assert_eq!(session.history().last().unwrap().smooth_maximum, session.accepted().evaluation.spectrum.smooth_maximum);
    drop(session); assert_eq!(study.operator().scales(), scales);
    let mut options = EnvelopeDesignOptions3::default(); options.optimizer.max_evaluations = 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut session = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy, options, &mut control).unwrap();
    let work = session.work();
    assert_eq!(session.run(10).unwrap().stop, ProjectedAlStop::EvaluationLimit);
    assert_eq!(session.work(), work); assert_eq!(session.history().len(), 1);
}

#[test]
fn recovery_cannot_override_a_raw_floor_or_renew_exhausted_physics() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.95; study.cells()];
    let basis = [&y[..], &z[..]]; let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let policy = normalized_policy(&mut study, &raw, &basis, &mut control);
    let options = EnvelopeDesignOptions3 { volume_limit: 0.1, density_floor: 0.9, ..Default::default() };
    let mut session = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy.clone(), options, &mut control).unwrap();
    let before = session.accepted().rho.clone();
    assert!(session.recover_feasible().is_err());
    assert_eq!(session.accepted().rho, before); assert!(session.best_feasible().is_none());
    drop(session);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(matches!(LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy, options, &mut control),
        Err(ProjectedAlError::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations,1); assert_eq!(study.operator().scales(),scales);
}
