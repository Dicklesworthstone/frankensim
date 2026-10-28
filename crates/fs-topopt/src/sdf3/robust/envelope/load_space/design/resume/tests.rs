use super::*;
use super::super::limits::{LoadSpaceCap3, LoadSpaceLimitedDesign3, LoadSpaceLimitedOptions3};
use crate::sdf3::robust::envelope::tests::fixture;
use crate::SolveBudget;

#[test]
fn g5_restored_load_space_search_matches_uninterrupted_accepted_designs() {
    let run = |restore: bool| {
        let (mut study, y, z) = fixture(2.0);
        let basis = [&y[..], &z[..]];
        let rho = vec![0.7; study.cells()];
        let options = LoadSpaceLimitedOptions3::default();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let initial = study.evaluate_robust_load_space(&rho, &basis, options.robustness,
            &LoadSpace3::new(2, 1.0, 0.02).unwrap(), &mut control).unwrap();
        let cap = LoadSpaceCap3::new(2, 1.8 * initial.worst_compliance, 0.02).unwrap();
        let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho,
            cap.clone(), options, &mut control).unwrap();
        session.run(2).unwrap();
        let checkpoint = session.checkpoint();
        let cost = if restore { checkpoint.restoration_cost() } else { 0 };
        if restore {
            drop(session);
            session = LoadSpaceLimitedDesign3::restore(&mut study, &basis, checkpoint,
                cap, options, &mut control).unwrap();
        }
        session.run(4).unwrap();
        let checkpoint = session.checkpoint();
        (session.history().to_vec(), session.accepted().rho.clone(),
            session.best_feasible().map(|p| p.rho.clone()),
            session.accepted().evaluation.worst_displacement.clone(), checkpoint, cost)
    };
    let reference = run(false);
    let restored = run(true);
    assert_eq!(reference.0, restored.0);
    assert_eq!(reference.1, restored.1);
    assert_eq!(reference.2, restored.2);
    assert_eq!(reference.3, restored.3);
    assert_eq!(restored.4.restoration_evaluations, restored.5);
    assert_eq!(restored.4.optimizer.work.evaluations,
        reference.4.optimizer.work.evaluations + restored.5);
    assert_eq!(reference.4.optimizer.sample, restored.4.optimizer.sample);
    assert_eq!(reference.4.optimizer.multiplier, restored.4.optimizer.multiplier);
    assert_eq!(reference.4.optimizer.penalty, restored.4.optimizer.penalty);
    assert_eq!(reference.4.optimizer.spectral_step, restored.4.optimizer.spectral_step);
}

#[test]
fn g3_restore_refuses_changed_response_history_and_missing_incumbent() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]];
    let raw = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let options = LoadSpaceLimitedOptions3::default();
    let baseline = study.evaluate_robust_load_space(&raw, &basis, options.robustness,
        &LoadSpace3::new(2, 1.0, 0.02).unwrap(), &mut control).unwrap();
    let cap = LoadSpaceCap3::new(2, 1.8 * baseline.worst_compliance, 0.02).unwrap();
    let checkpoint = {
        let mut run = LoadSpaceLimitedDesign3::new(&mut study, &basis, &raw,
            cap.clone(), options, &mut control).unwrap();
        run.run(2).unwrap(); run.checkpoint()
    };
    let scales = study.operator().scales().to_vec();
    let mut wrong_sample = checkpoint.clone();
    wrong_sample.optimizer.sample.constraint_gradient[0] += 0.1;
    let mut wrong_history = checkpoint.clone();
    wrong_history.history.pop();
    let mut lost_best = checkpoint.clone();
    assert!(lost_best.best_feasible_density.is_some());
    lost_best.best_feasible_density = None;
    for bad in [wrong_sample, wrong_history, lost_best] {
        assert!(LoadSpaceLimitedDesign3::restore(&mut study, &basis, bad,
            cap.clone(), options, &mut control).is_err());
        assert_eq!(study.operator().scales(), scales);
    }
    let changed: Vec<_> = y.iter().map(|v| 1.01 * v).collect();
    let changed_basis = [&changed[..], &z[..]];
    assert!(LoadSpaceLimitedDesign3::restore(&mut study, &changed_basis, checkpoint,
        cap, options, &mut control).is_err());
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn g4_restore_charges_endpoints_and_cannot_renew_exhausted_allowances() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]];
    let raw = vec![0.7; study.cells()];
    let cap = LoadSpaceCap3::new(2, 1e8, 0.02).unwrap();
    let options = LoadSpaceLimitedOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let checkpoint = LoadSpaceLimitedDesign3::new(&mut study, &basis, &raw,
        cap.clone(), options, &mut control).unwrap().checkpoint();
    let mut exhausted = options;
    exhausted.optimizer.max_evaluations = checkpoint.optimizer.work.evaluations;
    let before = control.work();
    assert!(LoadSpaceLimitedDesign3::restore(&mut study, &basis, checkpoint.clone(),
        cap.clone(), exhausted, &mut control).is_err());
    assert_eq!(control.work(), before);
    let scales = study.operator().scales().to_vec();
    let mut cancel = |_| ControlFlow::Break(());
    let mut cancelled = SolveControl::new(SolveBudget::default(), &mut cancel);
    assert!(matches!(LoadSpaceLimitedDesign3::restore(&mut study, &basis, checkpoint.clone(),
        cap.clone(), options, &mut cancelled), Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(study.operator().scales(), scales);
    let mut poll = |_| ControlFlow::Continue(());
    let mut tiny = SolveControl::new(SolveBudget { total_iterations: 1,
        ..SolveBudget::default() }, &mut poll);
    assert!(matches!(LoadSpaceLimitedDesign3::restore(&mut study, &basis, checkpoint,
        cap, options, &mut tiny), Err(ProjectedAlError::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(tiny.work().linear_iterations, 1);
    assert_eq!(study.operator().scales(), scales);
}
