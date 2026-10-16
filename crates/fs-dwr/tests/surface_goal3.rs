//! Real mixed reference loads and independent surface observation adjoints.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{ElasticityError3,ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::{ReferenceLoad3,SurfaceForce3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_dwr::elasticity3::{GoalFields3,GoalError3,estimate_goal3,estimate_reference_goal3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval{let d=if a==HeightAxis::Z{1.0}else{0.0};Interval::new(d,d)}
}
fn build(level:u8,surface:bool)->AdaptiveElasticity3{
    let tree=Octree3::uniform(level,4,4096).unwrap();let mut poll=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    let d=HexCell::try_new([0.0;3],[1.0;3]).unwrap();let m=IsotropicElastic::new(1.0,0.3,1.0).unwrap();
    if surface{AdaptiveElasticity3::build_with_surface(d,&tree,&Slab,&m,&|p|p[0]==0.0,ElasticityOptions3::default(),Default::default(),&mut q).unwrap()}
    else{AdaptiveElasticity3::build(d,&tree,&Slab,&m,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()}
}
fn pair()->(AdaptiveElasticity3,AdaptiveElasticity3){
    let mut c=build(1,true);let scales:Vec<_>=(0..c.cells()).map(|i|0.3+0.1*(i%5) as f64).collect();c.set_scales(&scales).unwrap();
    let mut f=build(2,true);let scales=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap().inherited_scales();f.set_scales(&scales).unwrap();(c,f)
}
fn body(p:[f64;3])->[f64;3]{[0.0,-0.2*(1.0+p[1]),0.0]}
fn pressure(p:[f64;3])->f64{1.0+0.2*p[0]+0.3*p[1]*p[1]}
fn observation(p:[f64;3],_:[f64;3])->[f64;3]{[0.0,1.0+0.4*p[0],0.2*p[1]]}
fn solve(op:&AdaptiveElasticity3,law:ReferenceLoad3<'_>)->Vec<f64>{
    let rhs=op.reference_load(law,||ControlFlow::Continue(())).unwrap();
    op.solve_controlled(&rhs,1e-10,20000,32,|_|ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
#[test]
fn g1_surface_and_mixed_compliance_reconstruct_the_actual_two_grid_change(){
    let (c,f)=pair();let t=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap();
    for body in [None,Some(&body as &dyn Fn([f64;3])->[f64;3])]{
        let law=ReferenceLoad3{body,surface:Some(SurfaceForce3::Pressure(&pressure))};let u=solve(&c,law);let v=solve(&f,law);
        let g=estimate_reference_goal3(&t,law,law,GoalFields3::compliance(&u,&v),Default::default(),||ControlFlow::Continue(())).unwrap();
        assert!(g.fine_value>0.0);assert!(g.identity_relative_defect<1e-8);
        assert!((g.correction()-(g.fine_value-g.coarse_value)).abs()<1e-8*g.fine_value);
        assert!(!g.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
        assert!(matches!(estimate_reference_goal3(&t,ReferenceLoad3::default(),law,GoalFields3::compliance(&u,&v),Default::default(),||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{field:"coarse-primal",..})));
    }
}
#[test]
fn g3_independent_surface_observation_requires_its_own_adjoint(){
    let (c,f)=pair();let load=ReferenceLoad3{body:Some(&body),surface:Some(SurfaceForce3::Pressure(&pressure))};let goal=ReferenceLoad3::traction(&observation);
    let uc=solve(&c,load);let uf=solve(&f,load);let zc=solve(&c,goal);let zf=solve(&f,goal);
    let t=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap();
    let fields=GoalFields3{coarse_primal:&uc,fine_primal:&uf,coarse_adjoint:&zc,fine_adjoint:&zf};
    let g=estimate_reference_goal3(&t,load,goal,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    assert!(g.identity_relative_defect<1e-8);assert!(g.field_residuals.iter().all(|r|*r<1e-8));
    assert!(matches!(estimate_reference_goal3(&t,load,goal,GoalFields3::compliance(&uc,&uf),Default::default(),||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{field:"coarse-adjoint",..})));
}
#[test]
fn g5_body_adapter_and_general_entry_point_replay_exactly(){
    let (c,f)=pair();let law=ReferenceLoad3::body(&body);let u=solve(&c,law);let v=solve(&f,law);let fields=GoalFields3::compliance(&u,&v);
    let t=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap();
    let a=estimate_goal3(&t,&body,&body,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    let b=estimate_reference_goal3(&t,law,law,fields,Default::default(),||ControlFlow::Continue(())).unwrap();
    assert_eq!(a.correction().to_bits(),b.correction().to_bits());assert_eq!(a.field_residuals,b.field_residuals);
    assert_eq!(a.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked,b.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked);
}
#[test]
fn g4_missing_surface_geometry_changed_load_and_cancellation_do_not_mint_evidence(){
    let c=build(1,false);let f=build(2,false);let u=solve(&c,ReferenceLoad3::body(&body));let v=solve(&f,ReferenceLoad3::body(&body));
    let t=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap();let law=ReferenceLoad3::pressure(&pressure);
    assert!(matches!(estimate_reference_goal3(&t,law,law,GoalFields3::compliance(&u,&v),Default::default(),||ControlFlow::Continue(())),Err(GoalError3::Physics(ElasticityError3::Invalid(_)))));
    let (c,f)=pair();let u=solve(&c,law);let v=solve(&f,law);let t=AdaptiveTransfer3::new(&c,&f,1_000_000,||ControlFlow::Continue(())).unwrap();
    let different=|p|2.0*pressure(p);let wrong=ReferenceLoad3::pressure(&different);
    assert!(matches!(estimate_reference_goal3(&t,wrong,wrong,GoalFields3::compliance(&u,&v),Default::default(),||ControlFlow::Continue(())),Err(GoalError3::FieldResidual{..})));
    let mut polls=0;let result=estimate_reference_goal3(&t,law,law,GoalFields3::compliance(&u,&v),Default::default(),||{polls+=1;if polls>40{ControlFlow::Break(())}else{ControlFlow::Continue(())}});
    assert!(matches!(result,Err(GoalError3::Physics(ElasticityError3::Cancelled))));
}
