//! The real projected optimizer, not an independent surrogate search.
use super::*;
use super::super::design::{ResponseStressDesignOptions3, ResponseStressDesignStudy3};
use fs_ascent::projected_al::ProjectedAlError;
fn options() -> ResponseStressDesignOptions3 {
    ResponseStressDesignOptions3 { density_floor:0.2,
        limits:ResponseStressOptions3{stress_limit:0.2,..Default::default()},..Default::default() }
}

#[test]
fn g1_joint_design_retains_less_material_with_original_limits_and_fields() {
    let (mut s,_,q)=fixture(false);let force=vec![0.0;s.operator.n()];let rho=vec![0.65;s.cells()];
    let a=limits(&q);let cases=[ResponseLimitCase3{force:&force,prescribed:Some(&motion),limits:&a}];
    let o=options();let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let mut session=ResponseStressDesignStudy3::new(&mut s,&cases,&rho,o,&mut c).unwrap();
    assert!(session.feasible());let initial=session.accepted().response.volume_fraction;
    session.run(8).unwrap();let best=session.best_feasible().unwrap();
    assert!(best.feasible(o.optimizer.tolerance));assert!(best.response.volume_fraction<initial-1e-4);
    assert_eq!(session.study().operator.scales(),session.accepted().stress.scales);
    let selected=session.into_best_feasible().unwrap();assert_eq!(s.operator.scales(),selected.stress.scales);
    let again=s.evaluate_response_stress_limits(&selected.response.rho,&cases,o.limits,&mut c).unwrap();same(&selected,&again);
}

#[test]
fn g5_split_runs_keep_the_same_joint_search_and_exact_work() {
    let (mut a,_,q)=fixture(false);let (mut b,_,_)=fixture(false);let f=vec![0.0;a.operator.n()];
    let limits=limits(&q);let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&limits}];
    let rho=vec![0.65;a.cells()];let o=options();
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let mut a=ResponseStressDesignStudy3::new(&mut a,&cases,&rho,o,&mut ca).unwrap();
    let mut b=ResponseStressDesignStudy3::new(&mut b,&cases,&rho,o,&mut cb).unwrap();
    a.run(4).unwrap();b.run(2).unwrap();b.run(2).unwrap();
    same(a.accepted(),b.accepted());same(a.best_feasible().unwrap(),b.best_feasible().unwrap());
    assert_eq!(a.history(),b.history());assert_eq!(a.optimizer_work(),b.optimizer_work());assert_eq!(a.work(),b.work());
    let before=b.work();b.run(0).unwrap();assert_eq!(b.work(),before);
}

#[test]
fn g4_cancellation_after_acceptance_preserves_the_complete_new_physical_state() {
    let (mut s,_,q)=fixture(false);let f=vec![0.0;s.operator.n()];let a=limits(&q);let rho=vec![0.65;s.cells()];
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];let enabled=Cell::new(false);
    let mut p=|s:SolveProgress| if enabled.get() && s.stage=="response-stress-design-accepted" {ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let mut session=ResponseStressDesignStudy3::new(&mut s,&cases,&rho,options(),&mut c).unwrap();
    enabled.set(true);
    assert!(matches!(session.run(2),Err(ProjectedAlError::Evaluation(ResponseStressError3::Evaluation(EvaluationStop::Cancelled)))));
    assert_eq!(session.optimizer_work().iterations,1);assert_eq!(session.history().len(),2);
    assert_eq!(session.study().operator.scales(),session.accepted().stress.scales);
    let endpoint=session.accepted().clone();drop(session);enabled.set(false);
    let again=s.evaluate_response_stress_limits(&endpoint.response.rho,&cases,options().limits,&mut c).unwrap();same(&endpoint,&again);
}

#[test]
fn g0_infeasible_requirements_and_foldback_admission_are_not_successes() {
    let (mut s,_,q)=fixture(false);let f=vec![0.0;s.operator.n()];let rho=vec![0.65;s.cells()];let calls=Cell::new(0);
    let law=|p,n| {calls.set(calls.get()+1);motion(p,n)};let a=limits(&q);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&law),limits:&a}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(ResponseStressDesignStudy3::new(&mut s,&cases,&rho,
        ResponseStressDesignOptions3{density_floor:1e-6,..options()},&mut c).is_err());
    assert_eq!(calls.get(),0);assert_eq!(c.work().linear_solves,0);
    let zero=vec![0.0;f.len()];let bad=[ResponseLimit3{response:BoundedResponse3::Displacement(&zero),
        sense:ResponseLimitSense3::AtLeast,bound:1.0,scale:1.0}];
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&bad}];
    let mut session=ResponseStressDesignStudy3::new(&mut s,&cases,&rho,options(),&mut c).unwrap();
    session.run(0).unwrap();assert!(!session.feasible());assert!(session.best_feasible().is_none());
    assert!(session.into_best_feasible().is_err());
}
