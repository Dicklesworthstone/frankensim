use super::*;

fn loose(q: &[f64]) -> [ResponseLimit3<'_>; 2] {
    [ResponseLimit3 { response: BoundedResponse3::Reaction(&right), sense: ResponseLimitSense3::AtMost, bound: 1e6, scale: 1.0 },
     ResponseLimit3 { response: BoundedResponse3::Displacement(q), sense: ResponseLimitSense3::AtMost, bound: 1e6, scale: 1.0 }]
}
#[test]
fn g1_scan_matches_independent_trials_and_keeps_the_smallest_feasible_endpoint() {
    let (mut s,f,q)=fixture(false); let rho:Vec<_>=(0..s.cells()).map(|i|0.25+0.055*(i%8) as f64).collect();
    let rows=loose(&q); let negative:Vec<_>=f.iter().map(|v|-v).collect();
    let cases=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&rows },
        ResponseLimitCase3 { force:&negative,prescribed:None,limits:&rows }];
    let thresholds=[0.45,0.25,0.6,0.6];
    let options=BinaryScanOptions3 { void_stiffness:0.01,..Default::default() };
    let before=s.operator().elasticity().scales().to_vec();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let scan=s.scan_binary_response_limits(&rho,&cases,&thresholds,options,&mut c).unwrap();
    assert!(scan.outcome.is_ok());assert_eq!(scan.candidates.len(),thresholds.len());
    assert_eq!(scan.work.linear_solves,1+2*thresholds.len());
    assert_eq!(s.operator().elasticity().scales(),before);
    let best=scan.best.as_ref().unwrap(); assert_eq!(scan.best_index,Some(2));
    for (i,&threshold) in thresholds.iter().enumerate() {
        let e=s.evaluate_binary_response_limits(&rho,&cases,BinaryEndpointOptions3 { threshold,
            void_stiffness:options.void_stiffness,limits:options.limits },&mut c).unwrap();
        assert_eq!(scan.candidates[i].constraint,e.constraint);
        assert_eq!(scan.candidates[i].volume_fraction,e.volume_fraction);
        assert!(e.volume_fraction>=best.volume_fraction);
        if i==2 {assert_eq!(best.displacements,e.displacements);assert_eq!(best.solid,e.solid);}
    }
    // Material rejection is independent of otherwise satisfied response bounds.
    let capped=s.scan_binary_response_limits(&rho,&cases,&[0.25],BinaryScanOptions3 { volume_cap:0.0,..options },&mut c).unwrap();
    assert!(capped.outcome.is_ok() && capped.best.is_none());
    assert!(capped.candidates[0].response_feasible && !capped.candidates[0].volume_feasible);
}
#[test]
fn g4_interruption_keeps_completed_feasible_fields_and_accounts_failed_trial_work() {
    for stage in ["response-binary-scan-accepted","response-limits-primal"] {
        let (mut s,f,q)=fixture(false);let rho=vec![0.5;s.cells()];let rows=loose(&q);
        let cases=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&rows }];
        let before=s.operator().elasticity().scales().to_vec();let finished=Cell::new(0);
        let mut p=|s:SolveProgress| {
            if s.stage=="response-binary-scan-accepted" {finished.set(finished.get()+1);}
            let interrupt=if stage=="response-binary-scan-accepted" {s.stage==stage}
                else {finished.get()==1 && s.stage==stage && s.solve_iterations>0};
            if interrupt {ControlFlow::Break(())}else{ControlFlow::Continue(())}
        };
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let r=s.scan_binary_response_limits(&rho,&cases,&[0.25,0.75],Default::default(),&mut c).unwrap();
        assert!(matches!(r.outcome,Err(ResponseError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(r.candidates.len(),1);assert_eq!(r.best_index,Some(0));
        assert!(r.best.as_ref().unwrap().feasible(0.0));
        assert_eq!(s.operator().elasticity().scales(),before);
        if stage=="response-limits-primal" {assert!(r.work.linear_iterations>r.best.as_ref().unwrap().work.linear_iterations);}
    }
}
#[test]
fn g0_invalid_late_threshold_or_exhausted_schedule_refuses_before_solving() {
    let (mut s,f,q)=fixture(false);let rho=vec![0.5;s.cells()];let rows=loose(&q);
    let cases=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&rows }];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    for thresholds in [&[][..],&[0.2,f64::NAN],&[0.2,0.0]] {
        assert!(s.scan_binary_response_limits(&rho,&cases,thresholds,Default::default(),&mut c).is_err());
    }
    assert!(s.scan_binary_response_limits(&rho,&cases,&[0.2,0.3],BinaryScanOptions3 { max_candidates:1,..Default::default() },&mut c).is_err());
    assert_eq!(c.work().linear_solves,0);
    let impossible=[ResponseLimit3 { bound:-1e6,..rows[0] }];
    let impossible=[ResponseLimitCase3 { force:&f,prescribed:Some(&motion),limits:&impossible }];
    let r=s.scan_binary_response_limits(&rho,&impossible,&[0.2,0.8],Default::default(),&mut c).unwrap();
    assert!(r.outcome.is_ok() && r.best.is_none() && r.best_index.is_none());
    assert!(r.candidates.iter().all(|r|!r.response_feasible));
}
