use super::*;
fn scan_options() -> BinaryStressScanOptions3 {
    BinaryStressScanOptions3 { stress: BinaryResponseStressOptions3 {stress_limit:100.0,..options()},
        ..Default::default() }
}
#[test]
fn g3_joint_scan_matches_explicit_endpoints_and_filters_only_once() {
    let (mut s,f,q)=fixture(true);let a=limits(&q);let rho=raw(&s);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];
    let thresholds=[0.1,0.5,0.9,0.5];let o=BinaryStressScanOptions3{volume_cap:0.75,..scan_options()};
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let before=s.operator.scales().to_vec();let report=s.scan_binary_response_stress_limits(&rho,&cases,&thresholds,o,&mut c).unwrap();
    assert!(report.outcome.is_ok());assert_eq!(report.candidates.len(),thresholds.len());
    assert_eq!(c.work().linear_solves,1+thresholds.len()*cases.len());
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let mut expected:Option<(usize,BinaryResponseStressEvaluation3)>=None;
    for (i,&t) in thresholds.iter().enumerate() {
        let mut policy=o.stress;policy.endpoint.threshold=t;
        let e=s.evaluate_binary_response_stress_limits(&rho,&cases,policy,&mut c).unwrap();
        let r=&report.candidates[i];assert_eq!(r.constraint,e.constraint);assert_eq!(r.stress_aggregate,e.stress_aggregate);
        assert_eq!(r.volume_fraction,e.response.volume_fraction);assert_eq!(r.joint_feasible,e.feasible(o.feasibility_tolerance));
        if e.feasible(o.feasibility_tolerance) && e.response.volume_fraction<=o.volume_cap
            && expected.as_ref().is_none_or(|(_,best)|e.response.volume_fraction<best.response.volume_fraction) {expected=Some((i,e));}
    }
    let (index,endpoint)=expected.unwrap();assert_eq!(report.best_index,Some(index));same(report.best.as_ref().unwrap(),&endpoint);
    assert_eq!(report.work.linear_solves+thresholds.len()-1,c.work().linear_solves);assert_eq!(s.operator.scales(),before);
}
#[test]
fn g0_response_feasible_candidates_are_not_selected_when_every_stress_check_fails() {
    let (mut s,f,q)=fixture(false);let a=limits(&q);let rho=raw(&s);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];let thresholds=[0.1,0.5,0.9];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let seed=s.scan_binary_response_stress_limits(&rho,&cases,&thresholds,scan_options(),&mut c).unwrap();
    assert!(seed.outcome.is_ok());let minimum=seed.candidates.iter().map(|r|r.stress_aggregate).reduce(f64::min).unwrap();
    assert!(minimum>0.0);
    let o=BinaryStressScanOptions3{stress:BinaryResponseStressOptions3{stress_limit:0.5*minimum,..scan_options().stress},..scan_options()};
    let report=s.scan_binary_response_stress_limits(&rho,&cases,&thresholds,o,&mut c).unwrap();
    assert!(report.outcome.is_ok());assert!(report.best.is_none());assert!(report.best_index.is_none());
    assert!(report.candidates.iter().all(|r|r.response_feasible&&!r.joint_feasible&&r.stress_residual>0.99));
}
#[test]
fn g4_scan_cancellation_keeps_only_complete_candidates_and_their_incumbent() {
    let (mut s,f,q)=fixture(false);let a=limits(&q);let rho=raw(&s);
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&a}];let before=s.operator.scales().to_vec();
    for stage in ["response-binary-stress-scan-retain","response-binary-stress-scan-accepted"] {
        let mut p=|p:SolveProgress|if p.stage==stage {ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let r=s.scan_binary_response_stress_limits(&rho,&cases,&[0.5,0.9],scan_options(),&mut c).unwrap();
        assert!(matches!(r.outcome,Err(ResponseStressError3::Evaluation(EvaluationStop::Cancelled))));
        assert!(r.work.linear_iterations>0);assert_eq!(s.operator.scales(),before);
        if stage.ends_with("accepted") {assert_eq!(r.candidates.len(),1);assert_eq!(r.best_index,Some(0));assert!(r.best.unwrap().feasible(0.0));}
        else {assert!(r.candidates.is_empty());assert!(r.best.is_none());}
    }
}
#[test]
fn g5_equal_volume_ties_keep_the_first_and_bad_schedules_spend_no_physics() {
    let (mut s,f,q)=fixture(false);let a=limits(&q);let rho=raw(&s);let calls=Cell::new(0);
    let g=|p,n|{calls.set(calls.get()+1);motion(p,n)};
    let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&g),limits:&a}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    for thresholds in [&[][..],&[0.5,f64::NAN][..],&[0.0][..]] {
        assert!(s.scan_binary_response_stress_limits(&rho,&cases,thresholds,scan_options(),&mut c).is_err());
    }
    assert!(s.scan_binary_response_stress_limits(&rho,&cases,&[0.5,0.9],BinaryStressScanOptions3{max_candidates:1,..scan_options()},&mut c).is_err());
    assert_eq!(calls.get(),0);assert_eq!(c.work().linear_solves,0);
    let a=s.scan_binary_response_stress_limits(&rho,&cases,&[0.5,0.5],scan_options(),&mut c).unwrap();
    assert!(a.outcome.is_ok());assert_eq!(a.best_index,Some(0));assert_eq!(a.candidates[0],a.candidates[1]);
    let b=s.scan_binary_response_stress_limits(&rho,&cases,&[0.5,0.5],scan_options(),&mut c).unwrap();
    assert_eq!(a.candidates,b.candidates);same(a.best.as_ref().unwrap(),b.best.as_ref().unwrap());
}
