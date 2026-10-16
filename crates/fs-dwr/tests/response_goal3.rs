//! G0/G1/G4/G5: actual embedded-motion elasticity and quadratic goal refinement.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3, ElasticityError3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_dwr::elasticity3::{GoalFields3, GoalError3};
use fs_dwr::elasticity3::response::{EquilibriumLoad3, ResponseObservation3, estimate_response_goal3, linearize_responses3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        if a==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}else{Interval::new(0.0,0.0)}
    }
}
fn build(tree:&Octree3)->AdaptiveElasticity3 {
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),tree,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,n|n[0]<0.0,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn pair(mixed:bool)->(AdaptiveElasticity3,AdaptiveElasticity3) {
    let tree=Octree3::uniform(1,4,4096).unwrap();
    let marks=if mixed {vec![*tree.leaves().iter().next().unwrap()]}else{tree.leaves().iter().copied().collect()};
    let mut c=build(&tree);let mut f=build(&tree.refined(&marks,||ControlFlow::Continue(())).unwrap());
    c.set_scales(&(0..c.cells()).map(|i|0.3+0.07*(i%7) as f64).collect::<Vec<_>>()).unwrap();
    let scales=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap().inherited_scales();f.set_scales(&scales).unwrap();(c,f)
}
fn shear(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.0,0.0,-0.2]}else{[0.0;3]}}
fn observe(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.0,0.0,1.0]}else{[0.0;3]}}
fn motion(p:[f64;3],_:[f64;3])->[f64;3] {[0.01,0.015*p[2],0.04*p[1]*p[2]]}
fn solve(op:&AdaptiveElasticity3,b:&[f64])->Vec<f64> {
    op.solve_controlled(b,1e-11,20_000,32,|_|ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
fn law()->EquilibriumLoad3<'static> {EquilibriumLoad3 {external:ReferenceLoad3::traction(&shear),prescribed:Some(&motion)}}
fn target(value:f64)->ResponseObservation3<'static> {ResponseObservation3 {functional:ReferenceLoad3::traction(&observe),target:value,scale:0.7,weight:1.3}}
fn states(c:&AdaptiveElasticity3,f:&AdaptiveElasticity3,obs:&[ResponseObservation3<'_>])->[Vec<f64>;4] {
    let uc=solve(c,&law().assemble(c,||ControlFlow::Continue(())).unwrap());let uf=solve(f,&law().assemble(f,||ControlFlow::Continue(())).unwrap());
    let zc=solve(c,&linearize_responses3(c,obs,&uc,8,||ControlFlow::Continue(())).unwrap().adjoint_rhs);
    let zf=solve(f,&linearize_responses3(f,obs,&uf,8,||ControlFlow::Continue(())).unwrap().adjoint_rhs);[uc,uf,zc,zf]
}
fn fields(s:&[Vec<f64>;4])->GoalFields3<'_> {GoalFields3 {coarse_primal:&s[0],fine_primal:&s[1],coarse_adjoint:&s[2],fine_adjoint:&s[3]}}
#[test]
fn g1_nonlinear_loss_identity_includes_motion_and_quadratic_remainder() {
    for mixed in [false,true] {
        let(c,f)=pair(mixed);let obs=[target(-0.05)];let s=states(&c,&f,&obs);
        let t=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap();
        let report=estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())).unwrap();
        let delta=report.fine_value-report.coarse_value;
        assert!((report.correction()-delta).abs()<1e-8*report.coarse_value.max(report.fine_value));
        assert!(report.quadratic_remainder<0.0);assert!(report.identity_relative_defect<1e-8);
        assert!((report.correction()-report.quadratic_remainder-delta).abs()>1e-12);
        let replay=estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())).unwrap();
        assert_eq!(report.correction().to_bits(),replay.correction().to_bits());
        assert_eq!(report.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked,replay.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked);
        assert!(estimate_response_goal3(&t,EquilibriumLoad3 {prescribed:None,..law()},&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())).is_err());
    }
}
#[test]
fn g0_exact_fine_fit_keeps_the_nonzero_curvature_marking_signal() {
    let(c,f)=pair(false);let preliminary=[target(0.0)];let initial=states(&c,&f,&preliminary);
    let value=linearize_responses3(&f,&preliminary,&initial[1],8,||ControlFlow::Continue(())).unwrap().responses[0];
    let obs=[target(value)];let s=states(&c,&f,&obs);assert!(s[3].iter().all(|v|*v==0.0));
    let t=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap();
    let r=estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())).unwrap();
    assert_eq!(r.fine_value,0.0);assert!(r.coarse_value>0.0&&r.quadratic_remainder<0.0);
    assert!(r.cells.values().any(|c|c.quadratic_remainder!=0.0));assert!(!r.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
}
#[test]
fn g0_stale_adjoint_wrong_targets_and_changed_physical_scales_are_refused() {
    let(c,mut f)=pair(true);let obs=[target(0.0)];let mut s=states(&c,&f,&obs);
    let t=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap();
    assert!(matches!(estimate_response_goal3(&t,law(),&[target(0.2)],fields(&s),Default::default(),8,||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{..})));
    s[2][0]+=0.1;
    assert!(matches!(estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{..})));
    drop(t);let scales=vec![0.9;f.cells()];f.set_scales(&scales).unwrap();
    let t=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap();
    assert!(matches!(estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||ControlFlow::Continue(())),Err(GoalError3::Invalid(_))));
}
#[test]
fn g4_cancellation_and_observation_caps_never_publish_partial_estimates() {
    let(c,f)=pair(true);let obs=[target(0.0)];let s=states(&c,&f,&obs);
    let t=AdaptiveTransfer3::new(&c,&f,100_000,||ControlFlow::Continue(())).unwrap();
    for limit in [0,20,500] {let mut calls=0;let r=estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),8,||{
        calls+=1;if calls>limit{ControlFlow::Break(())}else{ControlFlow::Continue(())}});
        assert!(matches!(r,Err(GoalError3::Physics(ElasticityError3::Cancelled))));}
    assert!(matches!(estimate_response_goal3(&t,law(),&obs,fields(&s),Default::default(),0,||ControlFlow::Continue(())),Err(GoalError3::Invalid(_))));
    let mut bad=obs;bad[0].scale=f64::NAN;assert!(linearize_responses3(&c,&bad,&s[0],8,||ControlFlow::Continue(())).is_err());
}
