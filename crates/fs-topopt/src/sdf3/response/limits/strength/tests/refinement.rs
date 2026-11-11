use super::*;
use crate::sdf3::response::limits::refinement::{ReferenceBoundedResponse3,
    ReferenceResponseLimit3,ReferenceResponseLimitCase3};
use crate::sdf3::response::limits::strength::refinement::{
    ResponseStressRefinementError3,ResponseStressRefinementOptions3};
use fs_cutfem::elastic3::surface::ReferenceLoad3;

fn body(p:[f64;3])->[f64;3]{[0.001*(1.0+p[1]),0.0,-0.0007]}
fn observe(p:[f64;3])->[f64;3]{[1.0+p[1],0.0,0.0]}

#[test]
fn g3_joint_enrichment_reuses_response_primals_and_replays_coarse_stress() {
    let (mut coarse,f,q)=fixture(false);let (fine,_,_)=fixture(true);let rho=raw(&coarse);
    let limits=limits(&q);let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&limits}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=coarse.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    let reference=[ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Reaction(&right),
        sense:limits[0].sense,bound:limits[0].bound,scale:limits[0].scale},
        ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
        sense:limits[1].sense,bound:limits[1].bound,scale:limits[1].scale}];
    let refs=[ReferenceResponseLimitCase3{load:ReferenceLoad3::body(&body),prescribed:Some(&motion),limits:&reference}];
    let before=coarse.operator.scales().to_vec();
    // One fine primal solve, read from the work counter: the stage's first event reports linear_solves before its solve is counted, its last after it.
    // The number of events is not the number of solves: SolveControl::solve_preconditioned (control.rs) polls at the start, at every loop turn and
    // after each batch of at most 32 iterations, so one solve reports at least three events.
    let mut goal_primal_solves:Option<(usize,usize)>=None;
    drop(c);let mut p=|w:SolveProgress|{if w.stage=="response-limit-goal-primal"{
        let n=w.work.linear_solves;goal_primal_solves=Some(goal_primal_solves.map_or((n,n),|(first,_)|(first,n)));}
        ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let options=ResponseStressRefinementOptions3::default();let tol=options.numerical.identity_tolerance;
    let report=coarse.estimate_response_stress_enrichment(fine.operator,&accepted,&refs,options,&mut c).unwrap();
    let (first,last)=goal_primal_solves.expect("the fine primal stage reports progress");
    assert_eq!(last-first,1);assert_eq!(report.response.cases.len(),1);
    // The retained coarse stress is replayed within the options' identity tolerance, the check the code itself makes (strength/refinement.rs `near`), not bit for bit:
    // the replay differs from the accepted aggregate by three units in the last place, 1.0e-17 (0.029471428621175858 against 0.02947142862117587), about a billionth of the 1e-8 tolerance.
    assert!((report.coarse_stress.aggregate-accepted.stress.aggregate).abs()<=tol*accepted.stress.aggregate.abs().max(1.0));
    assert!(report.fine_stress.aggregate>=report.fine_stress.sampled_relaxed_max);
    assert!(report.fine_stress.aggregate>=report.fine_stress.sampled_physical_max);
    assert!(report.fine_stress.point_count>0);assert_eq!(coarse.operator.scales(),before);
    assert!(report.max_normalized_change().is_finite());
}

#[test]
fn g0_changed_stress_policy_is_refused_before_enriched_physics() {
    let (mut coarse,f,q)=fixture(false);let (fine,_,_)=fixture(true);let rho=raw(&coarse);
    let limits=limits(&q);let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&limits}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=coarse.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    let reference=[ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Reaction(&right),
        sense:limits[0].sense,bound:limits[0].bound,scale:limits[0].scale},
        ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
        sense:limits[1].sense,bound:limits[1].bound,scale:limits[1].scale}];
    let refs=[ReferenceResponseLimitCase3{load:ReferenceLoad3::body(&body),prescribed:Some(&motion),limits:&reference}];
    let mut options=ResponseStressRefinementOptions3::default();options.joint.stress.aggregation_power=16.0;
    let before=c.work();
    assert!(matches!(coarse.estimate_response_stress_enrichment(fine.operator,&accepted,&refs,options,&mut c),
        Err(ResponseStressRefinementError3::Invalid(_))));
    assert_eq!(c.work(),before);
}

#[test]
fn g4_late_stress_sampling_cancellation_restores_source_material() {
    let (mut coarse,f,q)=fixture(false);let (fine,_,_)=fixture(true);let rho=raw(&coarse);
    let limits=limits(&q);let cases=[ResponseLimitCase3{force:&f,prescribed:Some(&motion),limits:&limits}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=coarse.evaluate_response_stress_limits(&rho,&cases,Default::default(),&mut c).unwrap();
    let reference=[ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Reaction(&right),
        sense:limits[0].sense,bound:limits[0].bound,scale:limits[0].scale},
        ReferenceResponseLimit3{response:ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
        sense:limits[1].sense,bound:limits[1].bound,scale:limits[1].scale}];
    let refs=[ReferenceResponseLimitCase3{load:ReferenceLoad3::body(&body),prescribed:Some(&motion),limits:&reference}];
    let incoming=coarse.operator.scales().to_vec();drop(c);
    let mut p=|w:SolveProgress|if w.stage=="response-stress-refinement-stress"{ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(coarse.estimate_response_stress_enrichment(fine.operator,&accepted,&refs,
        ResponseStressRefinementOptions3::default(),&mut c).is_err());
    assert_eq!(coarse.operator.scales(),incoming);
}
