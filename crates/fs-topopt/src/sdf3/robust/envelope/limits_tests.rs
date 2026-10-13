use super::*;
use super::super::tests::fixture;
use crate::pipeline::LoadCase;
use crate::{SolveBudget, SolveProgress};

fn caps(study: &mut CutDensityStudy3, rho: &[f64], y: &[f64], z: &[f64],
    control: &mut SolveControl<'_>) -> ComplianceCaps3 {
    let e = study.evaluate_robust(rho, &[LoadCase { force: y, weight: 0.5 },
        LoadCase { force: z, weight: 0.5 }], RobustOptions3::default(), control).unwrap();
    ComplianceCaps3::new(e.eroded.objective.case_compliances.iter().map(|c| 1.5*c).collect(), 0.03, 2).unwrap()
}

#[test]
fn minimum_material_descent_respects_each_eroded_cap_and_replays_the_fields() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.65; study.cells()];
    let keys = study.operator().cell_keys(); let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &raw, &y, &z, &mut control);
    let forces = [&y[..], &z[..]];
    let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits.clone(),
        ComplianceLimitedOptions3::default(), &mut control).unwrap();
    let initial = session.accepted().evaluation.physical.dilated_volume_fraction;
    assert!(session.gradient_check().numerical.passed());
    let outcome = session.run(40).unwrap();
    assert!(session.history().len() > 1);
    assert_eq!(outcome.objective, session.accepted().evaluation.physical.dilated_volume_fraction);
    assert_eq!(outcome.constraint, session.accepted().evaluation.envelope.value - 1.0);
    assert_eq!(session.state.sample().gradient, session.accepted().evaluation.physical.dilated_volume_gradient);
    assert_eq!(session.state.sample().constraint_gradient, session.accepted().evaluation.gradient);
    let best = session.best_feasible().unwrap().clone();
    assert!(best.evaluation.physical.dilated_volume_fraction < 0.99*initial);
    assert!(limits.satisfied_by(&best.evaluation.physical.eroded.objective.case_compliances));
    let work = session.work(); let optimizer_work = session.optimizer_work();
    let rechecked = session.recheck_best().unwrap().unwrap();
    assert_eq!(rechecked.rho, best.rho); assert_eq!(rechecked.evaluation.envelope, best.evaluation.envelope);
    assert!(session.work().linear_solves > work.linear_solves);
    assert_eq!(session.optimizer_work(), optimizer_work);
    drop(session);
    let fresh = study.evaluate_robust_envelope(&best.rho, &forces, RobustOptions3::default(),
        limits.envelope(), &mut control).unwrap();
    assert_eq!(fresh.envelope, best.evaluation.envelope);
    assert_eq!(fresh.physical.eroded.objective.displacements, best.evaluation.physical.eroded.objective.displacements);
    assert_eq!(fresh.physical.dilated_rho, best.evaluation.physical.dilated_rho);
    assert_eq!(study.operator().cell_keys(), keys); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn a_good_average_and_small_optimizer_violation_cannot_waive_physical_caps() {
    let limits = ComplianceCaps3::new(vec![1.0, 1.0], 0.03, 2).unwrap();
    assert!(!limits.satisfied_by(&[1.2, 0.1])); // Average 0.65 is irrelevant.
    assert!(!limits.satisfied_by(&[1.0 + 1e-12, 0.0]));
    assert!(!limits.satisfied_by(&[0.0]));
    assert!(!limits.satisfied_by(&[0.0, f64::NAN]));
    assert!(limits.satisfied_by(&[1.0, 0.0]));
    let (mut study, _, z) = fixture(1.0);
    let raw = vec![0.5; study.cells()];
    let options = ComplianceLimitedOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let e = study.evaluate_robust(&raw, &[LoadCase { force: &z, weight: 1.0 }], options.robustness, &mut control).unwrap();
    let limits = ComplianceCaps3::new(vec![e.eroded.objective.compliance/(1.0 + 0.25*options.optimizer.tolerance)], 0.03, 1).unwrap();
    let forces = [&z[..]];
    let session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits, options, &mut control).unwrap();
    assert!(session.history()[0].smooth_violation <= options.optimizer.tolerance);
    assert!(!session.history()[0].feasible);
    assert!(session.best_feasible().is_none());
}

#[test]
fn segmented_minimum_material_runs_do_not_restart_multipliers_or_budgets() {
    let run = |chunks: &[usize]| {
        let (mut study, y, z) = fixture(1.0); let raw = vec![0.65; study.cells()];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let limits = caps(&mut study, &raw, &y, &z, &mut control); let forces = [&y[..], &z[..]];
        let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits,
            ComplianceLimitedOptions3::default(), &mut control).unwrap();
        let before = (session.work(), session.optimizer_work()); session.run(0).unwrap();
        assert_eq!(before, (session.work(), session.optimizer_work()));
        for &n in chunks { if session.run(n).unwrap().stop != ProjectedAlStop::IterationLimit { break; } }
        (session.history().to_vec(), session.accepted().rho.clone(),
            session.best_feasible().unwrap().rho.clone(), session.work(), session.optimizer_work())
    };
    let whole = run(&[12]); assert!(whole.0.len() > 1);
    assert_eq!(whole, run(&[5, 7]));
}

#[test]
fn late_cancellation_preserves_accepted_fields_and_the_feasible_incumbent() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.65; study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-compliance-limited-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &raw, &y, &z, &mut control); let forces = [&y[..], &z[..]];
    let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits.clone(),
        ComplianceLimitedOptions3::default(), &mut control).unwrap();
    assert!(matches!(session.run(20), Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(session.history().len(), 2); assert_eq!(session.accepted().rho, session.state.point());
    assert_eq!(session.history()[1].dilated_volume, session.accepted().evaluation.physical.dilated_volume_fraction);
    assert!(limits.satisfied_by(&session.best_feasible().unwrap().evaluation.physical.eroded.objective.case_compliances));
    drop(session); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn a_cancelled_trial_is_charged_but_cannot_replace_the_initial_point() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.65; study.cells()];
    let armed = std::cell::Cell::new(false);
    let mut poll = |p: SolveProgress| if armed.get() && p.stage == "sdf3-envelope-adjoint" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &raw, &y, &z, &mut control); let forces = [&y[..], &z[..]];
    let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits,
        ComplianceLimitedOptions3::default(), &mut control).unwrap();
    let work = session.work(); armed.set(true);
    assert!(session.run(20).is_err());
    assert_eq!(session.history().len(), 1); assert_eq!(session.best_feasible().unwrap().rho, raw);
    assert!(session.work().linear_solves > work.linear_solves);
    assert_eq!(session.optimizer_work().evaluations, 2);
}

#[test]
fn protected_physical_regions_survive_material_minimization() {
    use crate::sdf3::PhysicalRegion3;
    let (study, y, z) = fixture(1.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let raw = vec![0.65; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &raw, &y, &z, &mut control); let forces = [&y[..], &z[..]];
    let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits,
        ComplianceLimitedOptions3::default(), &mut control).unwrap();
    session.run(8).unwrap();
    let e = &session.best_feasible().unwrap().evaluation.physical;
    for field in [&e.eroded.projected_rho, &e.nominal_rho, &e.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
}

#[test]
fn exhausted_evaluation_and_physics_budgets_are_not_renewed() {
    let (mut study, y, z) = fixture(1.0); let raw = vec![0.65; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let limits = caps(&mut study, &raw, &y, &z, &mut control); let forces = [&y[..], &z[..]];
    let mut options = ComplianceLimitedOptions3::default(); options.optimizer.max_evaluations = 1;
    let mut session = ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits.clone(), options, &mut control).unwrap();
    let before = session.work();
    assert_eq!(session.run(20).unwrap().stop, ProjectedAlStop::EvaluationLimit);
    assert_eq!(session.work(), before); assert_eq!(session.history().len(), 1);
    drop(session);
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(matches!(ComplianceLimitedDesign3::new(&mut study, &forces, &raw, limits, options, &mut control),
        Err(ProjectedAlError::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations, 1); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn invalid_caps_and_force_shapes_refuse_before_any_solve() {
    assert!(ComplianceCaps3::new(vec![1.0; 2], 2.0, 2).is_err());
    for cap in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(ComplianceCaps3::new(vec![cap], 0.03, 1).is_err());
    }
    let (mut study, _, z) = fixture(1.0); let raw = vec![0.65; study.cells()];
    let limits = ComplianceCaps3::new(vec![1000.0], 0.03, 1).unwrap();
    let bad_forces = [&z[..z.len()-1]];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(ComplianceLimitedDesign3::new(&mut study, &bad_forces, &raw, limits,
        ComplianceLimitedOptions3::default(), &mut control), Err(ProjectedAlError::Invalid(_))));
    assert_eq!(control.work().linear_solves, 0);
}
