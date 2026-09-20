//! G0/G3/G4/G5: real CutFEM response adjoints, not differentiated solver traces.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;
use fs_topopt::{SimpParams, SolveBudget, SolveControl, SolveProgress, EvaluationStop};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use fs_topopt::sdf3::response::{ResponseCase3, ResponseTarget3, ResponseOptions3, ResponseError3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]); (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        if a==HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(0.17,0.17)-Interval::new(0.83,0.83) }
        else { Interval::new(0.0,0.0) }
    }
}
fn build(mixed:bool)->AdaptiveElasticity3 {
    let t=Octree3::uniform(1,4,4096).unwrap();
    let t=if mixed { t.refined(&[*t.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap() } else {t};
    let mut p=|_|ControlFlow::Continue(()); let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&t,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn motion(p:[f64;3],n:[f64;3])->[f64;3] {
    if n[0]>0.0 {[0.02,0.0,0.0]} else {[0.0,0.0,0.007*p[1]]}
}
fn other_motion(p:[f64;3],n:[f64;3])->[f64;3] {
    if n[0]>0.0 {[0.0,0.01,0.005*p[2]]} else {[-0.01,0.0,0.0]}
}
fn vectors(op:&AdaptiveElasticity3)->(Vec<f64>,Vec<f64>,Vec<f64>) {
    let load=|f:&dyn Fn([f64;3])->[f64;3]|op.body_load(f,||ControlFlow::Continue(())).unwrap();
    (load(&|_|[0.002,0.0,-0.003]),load(&|p|[0.0,0.0,1.0+p[0]]),load(&|p|[1.0+p[1],0.0,0.0]))
}
fn study(op:AdaptiveElasticity3)->CutDensityStudy3<AdaptiveSolveSpace3> {
    CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default())
}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn close(a:&[f64],b:&[f64],tol:f64) {
    assert_eq!(a.len(),b.len()); let scale=b.iter().map(|v|v.abs()).fold(1e-12_f64,f64::max);
    assert!(a.iter().zip(b).all(|(a,b)|(a-b).abs()<=tol*scale),"{a:?} vs {b:?}");
}
#[test]
fn g3_full_motion_response_gradient_matches_every_density_resolve() {
    let op=build(false);let (f,q,r)=vectors(&op);let negative:Vec<_>=f.iter().map(|v|-0.5*v).collect();let mut s=study(op);
    let targets=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:0.7},ResponseTarget3{q:&r,target:-0.004,scale:0.02,weight:0.3}];
    let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&targets},ResponseCase3{force:&negative,prescribed:Some(&other_motion),targets:&targets}];
    let rho:Vec<_>=(0..s.cells()).map(|i|0.35+0.03*i as f64).collect();let original=s.operator().scales().to_vec();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let options=ResponseOptions3{volume_weight:0.03,..Default::default()};let a=s.evaluate_responses(&rho,&cases,options,&mut c).unwrap();
    for i in 0..rho.len() {
        let mut x=rho.clone();x[i]+=1e-4;let plus=s.evaluate_responses(&x,&cases,options,&mut c).unwrap();
        x[i]-=2e-4;let minus=s.evaluate_responses(&x,&cases,options,&mut c).unwrap();let fd=(plus.objective-minus.objective)/2e-4;
        assert!((fd-a.gradient[i]).abs()<2e-4*fd.abs().max(1e-7),"cell {i}: {fd} vs {}",a.gradient[i]);
    }
    assert_eq!(s.operator().scales(),original);assert_eq!(a.responses.len(),2);assert_ne!(a.displacements[0],a.displacements[1]);
    let b=s.evaluate_responses(&rho,&cases,options,&mut c).unwrap();assert_eq!(a.gradient,b.gradient);assert_eq!(a.displacements,b.displacements);
}
#[test]
fn g1_rigid_prescribed_translation_cancels_load_and_stiffness_derivatives() {
    let op=build(true);let (_,q,_)=vectors(&op);let zero=vec![0.0;op.n()];let mut s=study(op);
    let g=|_:[f64;3],_:[f64;3]|[0.03,-0.01,0.02];
    let targets=[ResponseTarget3{q:&q,target:0.0,scale:0.02,weight:1.0}];
    let cases=[ResponseCase3{force:&zero,prescribed:Some(&g),targets:&targets}];let rho=vec![0.5;s.cells()];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=s.evaluate_responses(&rho,&cases,Default::default(),&mut c).unwrap();
    for u in a.displacements[0].chunks_exact(3) {close(u,&[0.03,-0.01,0.02],1e-7);}
    assert!(a.gradient.iter().all(|g|g.abs()<1e-6));assert_eq!(a.external_work,[0.0]);
    let omitted=s.operator().elasticity().prescribed_displacement_scale_work(&g,&a.adjoints[0],||ControlFlow::Continue(())).unwrap();
    assert!(omitted.iter().any(|v|v.abs()>1e-5),"zero gradient must require the lifting derivative");
}
#[test]
fn g0_mixed_bilinear_contractions_match_operator_perturbations_with_unequal_fields() {
    let mut op=build(true);let scales=vec![0.6;op.cells()];op.set_scales(&scales).unwrap();
    let u:Vec<_>=(0..op.n()).map(|i|(i%11)as f64/11.0).collect();let z:Vec<_>=(0..op.n()).map(|i|1e-8*((i%7)as f64-3.0)).collect();
    let exact=op.scale_bilinear_forms(&z,&u,||ControlFlow::Continue(())).unwrap();
    let symmetric=op.scale_bilinear_forms(&u,&z,||ControlFlow::Continue(())).unwrap();close(&exact,&symmetric,1e-12);
    for i in 0..op.cells() {
        let mut s=scales.clone();s[i]+=1e-4;op.set_scales(&s).unwrap();let mut ax=vec![0.0;op.n()];op.apply(&u,&mut ax);let a=dot(&z,&ax);
        s[i]-=2e-4;op.set_scales(&s).unwrap();op.apply(&u,&mut ax);let b=dot(&z,&ax);
        assert!(((a-b)/2e-4-exact[i]).abs()<1e-8*exact[i].abs().max(1e-12));
    }
}
#[test]
fn g3_homogeneous_response_matches_existing_compliance_chain() {
    let op=build(false);let (f,_,_)=vectors(&op);let mut s=study(op);let rho=vec![0.5;s.cells()];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let base=s.evaluate(&rho,&[LoadCase{force:&f,weight:1.0}],&mut c).unwrap();
    let targets=[ResponseTarget3{q:&f,target:0.0,scale:1.0,weight:1.0}];
    let a=s.evaluate_responses(&rho,&[ResponseCase3{force:&f,prescribed:None,targets:&targets}],Default::default(),&mut c).unwrap();
    close(&[a.objective],&[0.5*base.objective.compliance.powi(2)],1e-9);
    close(&a.gradient,&base.objective.gradient.iter().map(|g|g*base.objective.compliance).collect::<Vec<_>>(),1e-7);
}
#[test]
fn g4_response_stops_never_replace_incoming_material_or_publish_partial_cases() {
    for stage in ["sdf3-response-case","sdf3-response-adjoint","sdf3-response-publish"] {
        let op=build(false);let(f,q,_)=vectors(&op);let mut s=study(op);let before=s.operator().scales().to_vec();
        let targets=[ResponseTarget3{q:&q,target:0.01,scale:0.01,weight:1.0}];let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&targets}];
        let mut p=|s:SolveProgress|if s.stage==stage&&(stage!="sdf3-response-adjoint"||s.solve_iterations>0){ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(s.evaluate_responses(&vec![0.5;s.cells()],&cases,Default::default(),&mut c),Err(ResponseError3::Evaluation(EvaluationStop::Cancelled))));
        assert_eq!(s.operator().scales(),before);
    }
}
#[test]
fn g0_late_invalid_observation_refuses_before_filter_or_elasticity_work() {
    let op=build(false);let(f,q,_)=vectors(&op);let mut s=study(op);
    let targets=[ResponseTarget3{q:&q,target:0.0,scale:1.0,weight:1.0}];let bad=[ResponseTarget3{scale:f64::NAN,..targets[0]}];
    let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&targets},ResponseCase3{force:&f,prescribed:None,targets:&bad}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(matches!(s.evaluate_responses(&vec![0.5;s.cells()],&cases,Default::default(),&mut c),Err(ResponseError3::Invalid(_))));
    assert_eq!(c.work().linear_solves,0);
}
