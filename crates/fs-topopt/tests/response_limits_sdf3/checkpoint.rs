//! Reuse the real minimum-material fixtures, solver and constraint evaluator.
use super::*;

#[test]
fn g5_reconstruction_preserves_search_history_and_dual_schedule_without_replay() {
    let (mut sa, mut f, _) = fixture(false); let (mut sb, _, _) = fixture(false); f.fill(0.0);
    let limits = lower(0.002); let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &limits }];
    let rho = vec![0.7; sa.cells()]; let options = MinimumMaterialOptions3::default();
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut a = MinimumMaterialStudy3::new(&mut sa, &cases, &rho, options, &mut ca).unwrap();
    a.run(4).unwrap();
    let mut b = MinimumMaterialStudy3::new(&mut sb, &cases, &rho, options, &mut cb).unwrap();
    b.run(2).unwrap(); let before = b.work(); let mut saved = b.checkpoint(); let cost = saved.restoration_cost(); drop(b);
    let mut b = MinimumMaterialStudy3::restore(&mut sb, &cases, &mut saved, options, &mut cb).unwrap();
    assert!(b.work().linear_iterations > before.linear_iterations);
    assert_eq!(b.optimizer_work().iterations, 2); b.run(2).unwrap();
    same(a.accepted(), b.accepted()); same(a.best_feasible().unwrap(), b.best_feasible().unwrap());
    assert_eq!(a.history(), b.history());
    let ac = a.checkpoint(); let mut bc = b.checkpoint();
    assert_eq!(bc.restoration_evaluations, cost);
    bc.optimizer.work.evaluations -= cost;
    assert_eq!(ac.optimizer, bc.optimizer, "spectral and multiplier states must survive reconstruction");
    let spent = b.work(); b.run(0).unwrap(); assert_eq!(b.work(), spent);
}

#[test]
fn g0_inactive_original_limits_and_policy_changes_refuse_before_physics() {
    let (mut s, mut f, q) = fixture(false); f.fill(0.0); let options = MinimumMaterialOptions3::default();
    let limits = [lower(0.002)[0], ResponseLimit3 { response: BoundedResponse3::Displacement(&q),
        sense: ResponseLimitSense3::AtMost, bound: 1e6, scale: 1.0 }];
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &limits }];
    let rho = vec![0.7; s.cells()];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
    let saved = session.checkpoint(); drop(session); let work = c.work();
    let mut changed = limits; changed[1].bound *= 2.0;
    let wrong = [ResponseLimitCase3 { limits: &changed, ..cases[0] }];
    assert!(matches!(MinimumMaterialStudy3::restore(&mut s, &wrong, &mut saved.clone(), options, &mut c), Err(ProjectedAlError::Invalid(_))));
    let mut altered = options; altered.limits.sharpness *= 2.0;
    assert!(matches!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved.clone(), altered, &mut c), Err(ProjectedAlError::Invalid(_))));
    let mut altered = options; altered.optimizer.tolerance *= 2.0;
    assert!(matches!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved.clone(), altered, &mut c), Err(ProjectedAlError::Invalid(_))));
    assert_eq!(c.work(), work);
}

#[test]
fn g4_cancellation_and_changed_endpoint_leave_material_intact_and_charge_retries() {
    for stage in ["minimum-material-restore-reserved", "response-limits-adjoint", "minimum-material-restore-publish"] {
        let (mut s, mut f, _) = fixture(false); f.fill(0.0);
        let rows = lower(0.002); let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let rho = vec![0.7; s.cells()]; let options = MinimumMaterialOptions3::default(); let enabled = Cell::new(false);
        let mut p = |p: SolveProgress| if enabled.get() && p.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
        let expected = session.accepted().clone(); let mut saved = session.checkpoint(); let evaluations = saved.optimizer.work.evaluations;
        let cost = saved.restoration_cost(); drop(session);
        let before = s.operator().elasticity().scales().to_vec(); enabled.set(true);
        assert!(matches!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved, options, &mut c),
            Err(ProjectedAlError::Evaluation(ResponseError3::Evaluation(EvaluationStop::Cancelled)))));
        assert_eq!(s.operator().elasticity().scales(), before);
        assert_eq!(saved.optimizer.work.evaluations, evaluations + cost);
        assert_eq!(saved.restoration_evaluations, cost);
        enabled.set(false);
        let restored = MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved, options, &mut c).unwrap();
        same(restored.accepted(), &expected);
        assert_eq!(restored.checkpoint().restoration_evaluations, 2*cost);
    }
}

#[test]
fn g0_stale_physics_and_missing_incumbent_are_not_repaired_from_cached_values() {
    let (mut s, mut f, _) = fixture(false); f.fill(0.0); let rows = lower(0.002);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.7; s.cells()]; let options = MinimumMaterialOptions3::default();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
    let saved = session.checkpoint(); drop(session); let before = s.operator().elasticity().scales().to_vec();
    let mut absent = saved.clone(); absent.best_feasible = None; let work = c.work();
    assert!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut absent, options, &mut c).is_err());
    assert_eq!(c.work(), work);
    // Deliberately corrupt one current derivative. It has valid shape but must
    // be rejected by the physical re-evaluation, not fed into the resumed AL.
    let mut stale = saved.clone(); stale.optimizer.sample.constraint_gradient[0] += 1.0;
    assert!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut stale, options, &mut c).is_err());
    assert_eq!(stale.restoration_evaluations, saved.restoration_cost());
    assert_eq!(s.operator().elasticity().scales(), before);
    // An observation witness is checked even when aggregate physics is intact.
    let mut stale = saved.clone(); stale.accepted.measurements[0][0] += 0.001;
    if let Some(best) = &mut stale.best_feasible { best.measurements[0][0] += 0.001; }
    assert!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut stale, options, &mut c).is_err());
    assert_eq!(s.operator().elasticity().scales(), before);
}

#[test]
fn g4_restoration_reservation_cannot_exceed_the_original_evaluation_budget() {
    let (mut s, mut f, _) = fixture(false); f.fill(0.0); let rows = lower(0.002);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.7; s.cells()]; let mut options = MinimumMaterialOptions3::default();
    options.optimizer.max_evaluations = 2;
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
    let mut saved = session.checkpoint(); drop(session);
    let restored = MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved, options, &mut c).unwrap();
    assert_eq!(restored.optimizer_work().evaluations, 2); let mut saved = restored.checkpoint(); drop(restored);
    let before = saved.optimizer.clone(); let spent = c.work();
    assert!(matches!(MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved, options, &mut c), Err(ProjectedAlError::Invalid(_))));
    assert_eq!(c.work(), spent); assert_eq!(saved.optimizer, before);
}

#[test]
fn g0_a_genuinely_infeasible_checkpoint_restores_without_inventing_an_incumbent() {
    let (mut s, mut f, _) = fixture(true); f.fill(0.0); let rows = lower(1.0);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.5; s.cells()]; let options = MinimumMaterialOptions3::default();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
    assert!(session.best_feasible().is_none()); let mut saved = session.checkpoint();
    let expected = session.accepted().clone(); drop(session);
    let restored = MinimumMaterialStudy3::restore(&mut s, &cases, &mut saved, options, &mut c).unwrap();
    same(restored.accepted(), &expected); assert!(restored.best_feasible().is_none());
    assert!(restored.into_best_feasible().is_err());
}
