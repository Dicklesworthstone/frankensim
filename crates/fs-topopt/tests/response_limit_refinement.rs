//! Real CutFEM limit observations, affine reactions, and original-grid marking.
use std::{cell::Cell, ops::ControlFlow};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::limits::*;
use fs_topopt::sdf3::response::limits::refinement::*;
use fs_topopt::sdf3_goal::GoalRefinementError3;

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self, lo:[f64;3], hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]); (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo:[f64;3], hi:[f64;3], a:HeightAxis)->Interval {
        if a==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)} else {Interval::new(0.0,0.0)}
    }
}
fn motion(p:[f64;3], _:[f64;3])->[f64;3] {[0.02*(p[0]-0.17),0.0,0.0]}
fn right(_:[f64;3], n:[f64;3])->[f64;3] {if n[0]>0.0 {[1.0,0.0,0.0]} else {[0.0;3]}}
fn body(_:[f64;3])->[f64;3] {[0.002,0.0,-0.001]}
fn observe(p:[f64;3])->[f64;3] {[1.0+p[1],0.0,0.0]}
fn tree(mixed:bool)->Octree3 {
    let t=Octree3::uniform(1,4,4096).unwrap();
    if mixed {t.refined(&[*t.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap()} else {t}
}
fn build(t:&Octree3)->AdaptiveSolveSpace3 {
    let mut p=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),t,&Slab,
        &IsotropicElastic::new(1.0,0.0,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op,100_000_000)
}
fn fine()->AdaptiveSolveSpace3 {build(&Octree3::uniform(2,4,4096).unwrap())}
fn study(mixed:bool)->Study {CutDensityStudy3::new(build(&tree(mixed)),0.15,SimpParams {penal:3.0,beta:2.0,..Default::default()})}
fn raw(s:&Study)->Vec<f64> {(0..s.cells()).map(|i|0.35+0.025*(i%7) as f64).collect()}
fn laws()->[ReferenceResponseLimit3<'static>;3] {
    [ReferenceResponseLimit3 {response:ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),sense:ResponseLimitSense3::AtMost,bound:0.006,scale:0.02},
     ReferenceResponseLimit3 {response:ReferenceBoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtLeast,bound:0.001,scale:0.002},
     ReferenceResponseLimit3 {response:ReferenceBoundedResponse3::Reaction(&right),sense:ResponseLimitSense3::AtMost,bound:0.004,scale:0.002}]
}
fn cases<'a>(limits:&'a [ReferenceResponseLimit3<'a>])->[ReferenceResponseLimitCase3<'a>;1] {
    [ReferenceResponseLimitCase3 {load:ReferenceLoad3::body(&body),prescribed:Some(&motion),limits}]
}
fn evaluate(s:&mut Study, laws:&[ReferenceResponseLimit3<'_>], c:&mut SolveControl<'_>)->ResponseLimitEvaluation3 {
    let op=s.operator().elasticity();
    let f=op.reference_load(ReferenceLoad3::body(&body),||ControlFlow::Continue(())).unwrap();
    let q=op.reference_load(ReferenceLoad3::body(&observe),||ControlFlow::Continue(())).unwrap();
    let limits:Vec<_>=laws.iter().map(|l|ResponseLimit3 {sense:l.sense,bound:l.bound,scale:l.scale,response:match l.response {
        ReferenceBoundedResponse3::Displacement(_)=>BoundedResponse3::Displacement(&q),
        ReferenceBoundedResponse3::Reaction(mode)=>BoundedResponse3::Reaction(mode),
    }}).collect();
    s.evaluate_response_limits(&raw(s),&[ResponseLimitCase3 {force:&f,prescribed:Some(&motion),limits:&limits}],Default::default(),c).unwrap()
}
fn close(a:f64,b:f64,scale:f64) {assert!((a-b).abs()<=1e-8*scale.max(1e-12),"{a} != {b}");}

#[test]
fn individual_limits_match_reintegrated_goals_and_full_affine_identities() {
    for mixed in [false,true] {
        let mut s=study(mixed);let laws=laws();let scales=s.operator().elasticity().scales().to_vec();
        let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let accepted=evaluate(&mut s,&laws,&mut c);
        let e=s.estimate_response_limit_enrichment(fine(),&accepted,&cases(&laws),Default::default(),&mut c).unwrap();
        let mut op=fine();
        let parent=AdaptiveTransfer3::new(s.operator().elasticity(),op.elasticity(),2_000_000,||ControlFlow::Continue(())).unwrap().parents().to_vec();
        let par=s.params();let values:Vec<_>=parent.iter().map(|&i|par.e_min+(1.0-par.e_min)*fs_math::det::pow(accepted.projected_rho[i],par.penal)).collect();
        op.set_scales(&values).unwrap();let u=&e.cases[0].fine_displacement;
        let rhs=op.elasticity().reference_load_with_motion(ReferenceLoad3::body(&body),Some(&motion),||ControlFlow::Continue(())).unwrap();
        assert!(op.elasticity().field_residual(u,&rhs,||ControlFlow::Continue(())).unwrap()<1e-8);
        let q=op.elasticity().reference_load(ReferenceLoad3::body(&observe),||ControlFlow::Continue(())).unwrap();
        let d:f64=q.iter().zip(u).map(|(q,u)|q*u).sum();
        let r=op.elasticity().embedded_reaction(u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value;
        for (j,change) in e.cases[0].limits.iter().enumerate() {
            close(change.goal.fine_value,if j==0 {d}else{r},laws[j].scale);
            close(change.fine_residual-change.coarse_residual,change.correction(),1.0);
        }
        close(e.cases[0].limits[1].correction(),-e.cases[0].limits[2].correction(),1.0);
        assert!(e.marking_mass.values().sum::<f64>()>0.0);
        assert!(e.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.iter().all(|l|s.operator().elasticity().leaves().contains(l)));
        assert_eq!(s.operator().elasticity().scales(),scales);
    }
}

#[test]
fn a_coarse_feasible_limit_can_fail_the_enriched_physics() {
    let mut s=study(false);let mut laws=[laws()[1]];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=evaluate(&mut s,&laws,&mut c);
    let probe=s.estimate_response_limit_enrichment(fine(),&a,&cases(&laws),Default::default(),&mut c).unwrap();
    let goal=&probe.cases[0].limits[0].goal;assert!((goal.fine_value-goal.coarse_value).abs()>1e-12);
    laws[0].bound=f64::midpoint(goal.coarse_value,goal.fine_value);
    laws[0].sense=if goal.fine_value>goal.coarse_value {ResponseLimitSense3::AtMost}else{ResponseLimitSense3::AtLeast};
    let a=evaluate(&mut s,&laws,&mut c);assert!(a.feasible(0.0));
    let e=s.estimate_response_limit_enrichment(fine(),&a,&cases(&laws),Default::default(),&mut c).unwrap();
    assert!(!e.fine_feasible(0.0));assert!(e.fine_worst_residual>0.0);assert!(e.max_normalized_change()>0.0);
}

#[test]
fn complete_family_must_match_before_any_preconditioner_is_prepared() {
    let mut s=study(false);let base=laws();let setups=Cell::new(0);
    let mut p=|p:SolveProgress|{if p.stage=="sdf3-preconditioner-start"{setups.set(setups.get()+1);}ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);let a=evaluate(&mut s,&base,&mut c);
    for field in 0..3 {
        let mut wrong=base;
        match field {0=>wrong[2].bound+=0.001,1=>wrong[2].sense=ResponseLimitSense3::AtLeast,_=>wrong[2].scale*=2.0}
        setups.set(0);
        assert!(matches!(s.estimate_response_limit_enrichment(fine(),&a,&cases(&wrong),Default::default(),&mut c),Err(GoalRefinementError3::Invalid(_))));
        assert_eq!(setups.get(),0);
    }
    let mut stale=a.clone();stale.displacements[0][0]+=1.0;setups.set(0);
    assert!(s.estimate_response_limit_enrichment(fine(),&stale,&cases(&base),Default::default(),&mut c).is_err());
    assert_eq!(setups.get(),0);
    let mut malformed=a.clone();malformed.measurements[0].pop();
    let work=c.work();assert!(s.estimate_response_limit_enrichment(fine(),&malformed,&cases(&base),Default::default(),&mut c).is_err());
    assert_eq!(c.work(),work);
}

#[test]
fn cancellation_and_work_exhaustion_leave_no_partial_estimate_or_changed_scales() {
    let mut s=study(false);let laws=laws();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=evaluate(&mut s,&laws,&mut c);let scales=s.operator().elasticity().scales().to_vec();
    for stage in ["response-limit-goal-primal","response-limit-goal-fine-adjoint","response-limit-goal-publish"] {
        let mut p=|p:SolveProgress|if p.stage==stage&&(stage.ends_with("publish")||p.solve_iterations>0){ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(s.estimate_response_limit_enrichment(fine(),&a,&cases(&laws),Default::default(),&mut c),Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
        assert!(c.work().linear_iterations>0);assert_eq!(s.operator().elasticity().scales(),scales);
    }
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:1,..Default::default()},&mut p);
    assert!(matches!(s.estimate_response_limit_enrichment(fine(),&a,&cases(&laws),Default::default(),&mut c),Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget{..}))));
    assert_eq!(c.work().linear_iterations,1);assert_eq!(s.operator().elasticity().scales(),scales);
}

#[test]
fn bound_direction_changes_do_not_cancel_marking_and_identical_runs_replay() {
    let mut s=study(false);let laws=laws();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);let a=evaluate(&mut s,&laws,&mut c);
    let run=|s:&mut Study| {let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        s.estimate_response_limit_enrichment(fine(),&a,&cases(&laws),Default::default(),&mut c).unwrap()};
    let x=run(&mut s);let y=run(&mut s);
    assert_eq!(x.work,y.work);assert_eq!(x.marking_mass,y.marking_mass);assert_eq!(x.fine_constraint.to_bits(),y.fine_constraint.to_bits());
    assert_eq!(x.cases[0].fine_displacement,y.cases[0].fine_displacement);
    assert_eq!(x.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked,y.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked);
    let mass:f64=x.marking_mass.values().sum();assert!(mass>0.0);
    assert!(!x.fine_feasible(f64::NAN));
}
