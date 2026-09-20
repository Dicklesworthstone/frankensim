//! Real SIMP, independent elastic loads, sparse hierarchy and DWR integration.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{AdaptiveMultilevelOptions3,AdaptiveSolveSpace3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop,MultiLoadOcOptions,SimpParams,SolveBudget,SolveControl,SolveProgress};
use fs_topopt::sdf3::{CutDensityStudy3,Sdf3Elasticity,controlled_sdf3_optimality_criteria};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3_goal::{GoalBodyLoad3,GoalPreconditioner3,GoalRefinementOptions3,GoalRefinementError3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],axis:HeightAxis)->Interval{
        let d=if axis==HeightAxis::Z{1.0}else{0.0};Interval::new(d,d)
    }
}
fn build(level:u8)->AdaptiveElasticity3{
    let tree=Octree3::uniform(level,5,8192).unwrap();let mut poll=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&tree,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()
}
fn space(options:AdaptiveMultilevelOptions3)->AdaptiveSolveSpace3{
    AdaptiveSolveSpace3::multilevel(build(2),&[&build(1),&build(0)],options,||ControlFlow::Continue(())).unwrap()
}
fn y(_:[f64;3])->[f64;3]{[0.0,-1.0,0.0]}
fn z(_:[f64;3])->[f64;3]{[0.0,0.0,-1.0]}
fn forces(op:&AdaptiveElasticity3)->(Vec<f64>,Vec<f64>){
    (op.body_load(&y,||ControlFlow::Continue(())).unwrap(),op.body_load(&z,||ControlFlow::Continue(())).unwrap())
}
fn loads(f:&(Vec<f64>,Vec<f64>))->[LoadCase<'_>;2]{[LoadCase{force:&f.0,weight:0.3},LoadCase{force:&f.1,weight:0.7}]}
fn params()->SimpParams{SimpParams{penal:3.0,beta:2.0,..Default::default()}}
fn close(a:&[f64],b:&[f64],tol:f64){
    assert_eq!(a.len(),b.len());let scale=b.iter().map(|v|v.abs()).fold(1e-30_f64,f64::max);
    assert!(a.iter().zip(b).all(|(a,b)|(a-b).abs()<tol*scale));
}
#[test]
fn g3_recursive_density_evaluation_preserves_objective_fields_and_full_chain_gradient(){
    let op=space(Default::default());let f=forces(op.elasticity());let mut study=CutDensityStudy3::new(op,0.15,params());
    let mut plain=CutDensityStudy3::new(build(2),0.15,params());let rho:Vec<f64>=(0..study.cells()).map(|i|0.35+0.02*(i%9) as f64).collect();
    let mut starts=0;let mut poll=|p:SolveProgress|{if p.stage=="sdf3-preconditioner-start"{starts+=1;}ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);let a=study.evaluate(&rho,&loads(&f),&mut c).unwrap();
    assert!(c.work().preconditioner_galerkin_products>0);assert_eq!(c.work().preconditioner_operator_applications,0);drop(c);assert_eq!(starts,1);
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let b=plain.evaluate(&rho,&loads(&f),&mut c).unwrap();close(&[a.objective.compliance],&[b.objective.compliance],1e-8);
    close(&a.objective.gradient,&b.objective.gradient,1e-7);assert_eq!(a.projected_rho,b.projected_rho);
    for (a,b) in a.objective.displacements.iter().zip(&b.objective.displacements){close(a,b,1e-7);}
    for i in [0,7,19,33]{
        let mut plus=rho.clone();let mut minus=rho.clone();plus[i]+=1e-4;minus[i]-=1e-4;
        let cp=study.evaluate(&plus,&loads(&f),&mut c).unwrap().objective.compliance;
        let cm=study.evaluate(&minus,&loads(&f),&mut c).unwrap().objective.compliance;
        let fd=(cp-cm)/2e-4;assert!((fd-a.objective.gradient[i]).abs()<2e-4*fd.abs());
    }
}
#[test]
fn g1_real_optimizer_descends_and_replays_without_rebuilding_geometry(){
    let op=space(Default::default());let f=forces(op.elasticity());let nodes=op.elasticity().nodes().to_vec();let sizes=op.level_sizes();
    let mut study=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.5;study.cells()];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,&loads(&f),&rho,MultiLoadOcOptions{max_iterations:3,..Default::default()},&mut c);
    assert!(report.history.len()>1,"{report:?}");assert!(report.history.last().unwrap().compliance<report.history[0].compliance);
    assert!(report.history.iter().all(|r|r.volume_fraction<=0.50000001));
    for p in report.history.windows(2){assert!(p[1].compliance<=p[0].compliance);}
    assert!(report.work.preconditioner_galerkin_products>0);assert_eq!(report.work.preconditioner_operator_applications,0);
    let again=study.evaluate(&report.rho,&loads(&f),&mut c).unwrap();assert_eq!(again.objective.displacements,report.displacements);
    assert_eq!(study.operator().elasticity().nodes(),nodes);assert_eq!(study.operator().level_sizes(),sizes);
}
#[test]
fn g4_sparse_setup_caps_and_cancellation_preserve_scales_and_charge_discarded_work(){
    let mut options=AdaptiveMultilevelOptions3::default();options.hierarchy.max_galerkin_products=3;
    let op=space(options);let f=forces(op.elasticity());let mut study=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.5;study.cells()];
    let before=study.operator().scales().to_vec();let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    assert!(matches!(study.evaluate(&rho,&loads(&f),&mut c),Err(EvaluationStop::TotalBudget{stage:"sdf3-preconditioner-setup"})));
    assert_eq!(c.work().preconditioner_galerkin_products,3);assert_eq!(study.operator().scales(),before);
    let mut study=CutDensityStudy3::new(space(Default::default()),0.15,params());
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let accepted=study.evaluate(&rho,&loads(&f),&mut c).unwrap();let before=study.operator().scales().to_vec();
    let mut poll=|p:SolveProgress|if p.work.preconditioner_galerkin_products>0{ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    assert!(matches!(study.evaluate(&vec![0.4;study.cells()],&loads(&f),&mut c),Err(EvaluationStop::Cancelled)));
    assert!(c.work().preconditioner_galerkin_products>0);assert_eq!(study.operator().scales(),before);
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    assert_eq!(study.evaluate(&rho,&loads(&f),&mut c).unwrap().objective.displacements,accepted.objective.displacements);
}
#[test]
fn g3_recursive_goal_refinement_reuses_physics_and_refuses_unrequested_ladders(){
    let op=build(1);let f=forces(&op);let mut study=CutDensityStudy3::new(op,0.15,params());let rho=vec![0.5;study.cells()];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let accepted=study.evaluate(&rho,&loads(&f),&mut c).unwrap();let body=[GoalBodyLoad3{density:&y,weight:0.3},GoalBodyLoad3{density:&z,weight:0.7}];
    let plain=study.estimate_compliance_enrichment(build(2),&body,&accepted.objective.displacements,Default::default(),&mut c).unwrap();
    let before=study.operator().scales().to_vec();let options=GoalRefinementOptions3{
        preconditioner:GoalPreconditioner3::Multilevel{options:Default::default()},..Default::default()};
    let goal=study.estimate_compliance_enrichment_with_coarse_levels(build(2),&[&build(0)],&body,&accepted.objective.displacements,options,&mut c).unwrap();
    close(&[goal.fine_value,goal.coarse_value,goal.correction],&[plain.fine_value,plain.coarse_value,plain.correction],1e-7);
    assert!(goal.work.preconditioner_galerkin_products>0);assert_eq!(study.operator().scales(),before);
    assert!(!goal.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
    let spent=c.work();
    assert!(matches!(study.estimate_compliance_enrichment_with_coarse_levels(build(2),&[&build(0)],&body,&accepted.objective.displacements,Default::default(),&mut c),Err(GoalRefinementError3::Invalid(_))));
    assert_eq!(c.work(),spent);
}
