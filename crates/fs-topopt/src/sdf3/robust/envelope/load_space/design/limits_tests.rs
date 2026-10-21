use super::*;
use crate::sdf3::robust::envelope::tests::fixture;
use crate::{SolveBudget, SolveProgress};
use crate::sdf3::PhysicalRegion3;

fn measured_cap<O: Sdf3Elasticity>(study: &mut CutDensityStudy3<O>, basis: &[&[f64]],
    rho: &[f64], factor: f64, control: &mut SolveControl<'_>) -> LoadSpaceCap3 {
    let p = LoadSpace3::new(basis.len(), 1.0, 0.02).unwrap();
    let e = study.evaluate_robust_load_space(rho, basis, RobustOptions3::default(), &p, control).unwrap();
    LoadSpaceCap3::new(basis.len(), factor*e.worst_compliance, 0.02).unwrap()
}

#[test]
fn basis_caps_cannot_hide_a_failed_simultaneous_combination() {
    let (mut study, y, _) = fixture(2.0);
    let rho = vec![0.6; study.cells()];
    let basis = [&y[..], &y[..]];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_robust_load_space(&rho, &basis, RobustOptions3::default(),
        &LoadSpace3::new(2, 1.0, 0.02).unwrap(), &mut control).unwrap();
    let cap = LoadSpaceCap3::new(2, 1.1*e.physical.eroded.objective.case_compliances[0], 0.02).unwrap();
    let mut e = study.evaluate_robust_load_space(&rho, &basis, RobustOptions3::default(), cap.policy(), &mut control).unwrap();
    assert!(e.physical.eroded.objective.case_compliances.iter().all(|c| *c < cap.cap()));
    assert!(e.worst_compliance > cap.cap());
    assert!(!cap.admits(&e, 0.9));
    // Even a purported rounded normalized ratio cannot waive physical excess.
    e.spectrum.maximum = 1.0; e.spectrum.smooth_maximum = 1.0;
    e.worst_compliance = f64::from_bits(cap.cap().to_bits()+1);
    assert!(!cap.admits(&e, 1e-6));
    for (n, c, t) in [(0, 1.0, 0.02), (17, 1.0, 0.02), (2, f64::NAN, 0.02), (2, 1.0, 2.0)] {
        assert!(LoadSpaceCap3::new(n, c, t).is_err());
    }
}

#[test]
fn material_descent_rechecks_the_complete_load_set_not_a_frozen_direction() {
    let (mut study, y, z) = fixture(2.0);
    let second: Vec<_> = y.iter().zip(&z).map(|(y,z)| 0.4*y+z).collect();
    let basis = [&y[..], &second[..]];
    let rho = vec![0.75; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = measured_cap(&mut study, &basis, &rho, 1.8, &mut control);
    let options = LoadSpaceLimitedOptions3::default();
    let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap.clone(), options, &mut control).unwrap();
    assert!(session.gradient_check().numerical.passed());
    let first = session.accepted().evaluation.physical.dilated_volume_fraction;
    let terminal = session.run(12).unwrap();
    assert!(session.history().len() > 1);
    let best = session.best_feasible().unwrap().clone();
    assert!(cap.admits(&best.evaluation, options.optimizer.tolerance));
    assert!(best.evaluation.physical.dilated_volume_fraction < first);
    assert_eq!(terminal.objective, session.accepted().evaluation.physical.dilated_volume_fraction);
    assert_eq!(terminal.constraint, session.accepted().evaluation.spectrum.smooth_maximum-1.0);
    let before = (session.history().to_vec(), session.optimizer_work(), session.accepted().rho.clone());
    let rechecked = session.recheck_best().unwrap().unwrap();
    assert_eq!(session.history(), before.0);
    assert_eq!(session.optimizer_work(), before.1);
    assert_eq!(session.accepted().rho, before.2);
    assert_eq!(rechecked.evaluation.worst_displacement, best.evaluation.worst_displacement);
    drop(session);
    let direct = study.evaluate_robust(&best.rho, &[LoadCase { force: &best.evaluation.worst_force, weight: 1.0 }],
        options.robustness, &mut control).unwrap();
    assert!((direct.eroded.objective.compliance-best.evaluation.worst_compliance).abs()
        < 1e-8*best.evaluation.worst_compliance);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn segmented_material_runs_share_the_same_optimizer_and_physical_work() {
    let run = |segments: &[usize]| {
        let (mut study, y, z) = fixture(2.0);
        let basis = [&y[..], &z[..]]; let rho = vec![0.7; study.cells()];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let cap = measured_cap(&mut study, &basis, &rho, 1.8, &mut control);
        let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap,
            LoadSpaceLimitedOptions3::default(), &mut control).unwrap();
        let work = session.work(); session.run(0).unwrap(); assert_eq!(work, session.work());
        for &steps in segments {
            if session.run(steps).unwrap().stop != ProjectedAlStop::IterationLimit { break; }
        }
        (session.history().to_vec(), session.accepted().rho.clone(), session.best_feasible().unwrap().rho.clone(),
            session.work(), session.optimizer_work())
    };
    assert_eq!(run(&[10]), run(&[4, 6]));
}

#[test]
fn late_cancellation_keeps_an_accepted_cap_bound_prefix() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]]; let rho = vec![0.7; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = measured_cap(&mut study, &basis, &rho, 1.8, &mut control);
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-design-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap.clone(),
        LoadSpaceLimitedOptions3::default(), &mut control).unwrap();
    assert!(session.run(8).is_err());
    assert_eq!(session.history().len(), 2);
    assert_eq!(session.optimizer_work().iterations, 1);
    assert!(cap.admits(&session.best_feasible().unwrap().evaluation, 1e-6));
    assert_eq!(session.accepted().evaluation.physical.dilated_volume_fraction, session.history()[1].dilated_volume);
    drop(session); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn infeasible_starts_and_limits_do_not_publish_a_feasible_design_or_renew_work() {
    let (mut study, y, z) = fixture(2.0);
    let basis = [&y[..], &z[..]]; let rho = vec![0.6; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = measured_cap(&mut study, &basis, &rho, 0.7, &mut control);
    let mut options = LoadSpaceLimitedOptions3::default(); options.optimizer.max_evaluations = 1;
    let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap.clone(), options, &mut control).unwrap();
    assert!(session.best_feasible().is_none());
    assert_eq!(session.run(5).unwrap().stop, ProjectedAlStop::EvaluationLimit);
    let work = session.work();
    assert_eq!(session.run(5).unwrap().stop, ProjectedAlStop::EvaluationLimit);
    assert_eq!(session.work(), work);
    assert!(session.recheck_best().unwrap().is_none());
    drop(session);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap, options, &mut control).is_err());
    assert_eq!(control.work().linear_iterations, 1);
}

#[test]
fn prescribed_regions_remain_physical_in_the_minimum_material_problem() {
    let (study, y, z) = fixture(2.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let basis = [&y[..], &z[..]]; let rho = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let cap = measured_cap(&mut study, &basis, &rho, 1.8, &mut control);
    let mut session = LoadSpaceLimitedDesign3::new(&mut study, &basis, &rho, cap,
        LoadSpaceLimitedOptions3::default(), &mut control).unwrap();
    session.run(4).unwrap();
    let e = &session.best_feasible().unwrap().evaluation.physical;
    for field in [&e.eroded.projected_rho, &e.nominal_rho, &e.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
}
