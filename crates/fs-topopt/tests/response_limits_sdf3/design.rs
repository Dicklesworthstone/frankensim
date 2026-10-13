use super::*;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlStop};
use fs_topopt::sdf3::response::limits::design::*;

fn lower(bound: f64) -> [ResponseLimit3<'static>; 1] {
    [ResponseLimit3 { response: BoundedResponse3::Reaction(&right),
        sense: ResponseLimitSense3::AtLeast, bound, scale: 0.002 }]
}
fn same(a: &ResponseLimitEvaluation3, b: &ResponseLimitEvaluation3) {
    assert_eq!(a.rho, b.rho); assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.constraint.to_bits(), b.constraint.to_bits());
    assert_eq!(a.constraint_gradient, b.constraint_gradient);
    assert_eq!(a.measurements, b.measurements); assert_eq!(a.residuals, b.residuals);
    assert_eq!(a.displacements, b.displacements); assert_eq!(a.adjoints, b.adjoints);
}

#[test]
fn g1_minimum_material_reduces_volume_without_losing_the_original_force_limit() {
    let (mut s, mut f, _) = fixture(false); f.fill(0.0); let r = lower(0.002);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &r }];
    let rho = vec![0.7; s.cells()]; let options = MinimumMaterialOptions3::default();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, options, &mut c).unwrap();
    let initial = session.accepted().volume_fraction;
    session.run(12).unwrap();
    let best = session.best_feasible().unwrap();
    assert!(best.feasible(options.optimizer.tolerance)); assert!(best.volume_fraction < initial-1e-4);
    assert!(best.measurements[0][0] >= r[0].bound-options.optimizer.tolerance*r[0].scale);
    let expected = best.clone(); let spent = session.work();
    let installed = session.into_best_feasible().unwrap(); same(&installed, &expected); assert_eq!(installed.work, spent);
    let scales = s.operator().elasticity().scales().to_vec();
    let again = s.evaluate_response_limits(&installed.rho, &cases, options.limits, &mut c).unwrap();
    same(&installed, &again); assert_eq!(s.operator().elasticity().scales(), scales);
}

#[test]
fn g5_split_runs_keep_the_same_optimizer_state_history_and_work() {
    let (mut a, mut f, _) = fixture(false); let (mut b, _, _) = fixture(false); f.fill(0.0);
    let rows = lower(0.002); let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.7; a.cells()];
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut a = MinimumMaterialStudy3::new(&mut a, &cases, &rho, Default::default(), &mut ca).unwrap();
    let mut b = MinimumMaterialStudy3::new(&mut b, &cases, &rho, Default::default(), &mut cb).unwrap();
    a.run(4).unwrap(); b.run(2).unwrap(); b.run(2).unwrap();
    same(a.accepted(), b.accepted()); same(a.best_feasible().unwrap(), b.best_feasible().unwrap());
    assert_eq!(a.history(), b.history()); assert_eq!(a.work(), b.work()); assert_eq!(a.optimizer_work(), b.optimizer_work());
    let before = b.work(); b.run(0).unwrap(); assert_eq!(b.work(), before);
}

#[test]
fn g4_cancelled_probes_and_accepted_transitions_leave_matching_material_and_fields() {
    for stage in ["response-limits-adjoint", "response-limit-design-accepted"] {
        let (mut s, mut f, _) = fixture(false); f.fill(0.0); let rows = lower(0.002);
        let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
        let rho = vec![0.7; s.cells()]; let cancel = Cell::new(false);
        let mut p = |s: SolveProgress| if cancel.get() && s.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let mut session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, Default::default(), &mut c).unwrap();
        let initial = session.accepted().clone(); cancel.set(true);
        assert!(session.run(2).is_err()); cancel.set(false);
        let retained = session.accepted().clone();
        if stage == "response-limits-adjoint" { same(&initial, &retained); }
        else { assert_eq!(session.history().len(), 2); assert_ne!(initial.rho, retained.rho); }
        let scales = session.study().operator().elasticity().scales().to_vec(); drop(session);
        let again = s.evaluate_response_limits(&retained.rho, &cases, Default::default(), &mut c).unwrap();
        same(&again, &retained); assert_eq!(s.operator().elasticity().scales(), scales);
    }
}

#[test]
fn g0_an_infeasible_stop_cannot_produce_a_feasible_material_design() {
    let (mut s, mut f, _) = fixture(false); f.fill(0.0); let rows = lower(1.0);
    let cases = [ResponseLimitCase3 { force: &f, prescribed: Some(&motion), limits: &rows }];
    let rho = vec![0.5; s.cells()];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut session = MinimumMaterialStudy3::new(&mut s, &cases, &rho, Default::default(), &mut c).unwrap();
    let report = session.run(0).unwrap();
    assert_ne!(report.stop, ProjectedAlStop::Converged);
    assert!(!session.accepted().feasible(1e-6)); assert!(session.best_feasible().is_none());
    assert!(session.into_best_feasible().is_err());
    let mut bad = MinimumMaterialOptions3::default(); bad.volume_scale = f64::NAN;
    assert!(matches!(MinimumMaterialStudy3::new(&mut s, &cases, &rho, bad, &mut c), Err(ProjectedAlError::Invalid(_))));
}
