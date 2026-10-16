//! Real embedded supports through the existing density/preconditioning/DWR stack.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{MultiLoadOcOptions,MultiLoadOcTermination,SimpParams,SolveBudget,SolveControl,SolveProgress};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,controlled_sdf3_optimality_criteria};
use fs_topopt::sdf3_goal::{GoalPreconditioner3,GoalReferenceLoad3,GoalRefinementOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval {
        if a==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(0.17,0.17)-Interval::new(0.83,0.83)}else{Interval::new(0.0,0.0)}
    }
}
fn build(level:u8)->AdaptiveElasticity3 {
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&Octree3::uniform(level,4,4096).unwrap(),&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,n|n[0]<0.0,ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn x(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[-1.0,0.0,0.0]}else{[0.0;3]}}
fn z(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.0,0.0,-1.0]}else{[0.0;3]}}
fn forces(op:&AdaptiveElasticity3)->(Vec<f64>,Vec<f64>) {
    (op.surface_load(&x,||ControlFlow::Continue(())).unwrap().rhs,op.surface_load(&z,||ControlFlow::Continue(())).unwrap().rhs)
}
fn loads(f:&(Vec<f64>,Vec<f64>))->[LoadCase<'_>;2] {[LoadCase{force:&f.0,weight:0.5},LoadCase{force:&f.1,weight:0.5}]}
fn close(a:&[f64],b:&[f64],tol:f64){let scale=b.iter().map(|v|v.abs()).fold(1e-30_f64,f64::max);assert_eq!(a.len(),b.len());assert!(a.iter().zip(b).all(|(a,b)|(a-b).abs()<tol*scale));}
#[test]
fn g3_embedded_boundaries_enter_all_preconditioners_and_the_full_density_pullback() {
    let op=build(1);let f=forces(&op);let rho:Vec<_>=(0..op.cells()).map(|i|0.35+0.03*i as f64).collect();
    let mut plain=CutDensityStudy3::new(op,0.15,SimpParams::default());let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let expected=plain.evaluate(&rho,&loads(&f),&mut c).unwrap();
    for mode in 0..3 {
        let op=build(1);let coarse=build(0);
        let space=match mode {
            0=>AdaptiveSolveSpace3::jacobi(op,100_000_000),
            1=>AdaptiveSolveSpace3::two_level(op,&coarse,Default::default(),||ControlFlow::Continue(())).unwrap(),
            _=>AdaptiveSolveSpace3::multilevel(op,&[&coarse],Default::default(),||ControlFlow::Continue(())).unwrap(),
        };
        let mut study=CutDensityStudy3::new(space,0.15,SimpParams::default());
        let a=study.evaluate(&rho,&loads(&f),&mut c).unwrap();close(&[a.objective.compliance],&[expected.objective.compliance],1e-8);close(&a.objective.gradient,&expected.objective.gradient,1e-7);
        for (u,v) in a.objective.displacements.iter().zip(&expected.objective.displacements){close(u,v,1e-7);}
        for i in [0,4,7] {
            let mut p=rho.clone();let mut m=rho.clone();p[i]+=1e-4;m[i]-=1e-4;
            let cp=study.evaluate(&p,&loads(&f),&mut c).unwrap().objective.compliance;let cm=study.evaluate(&m,&loads(&f),&mut c).unwrap().objective.compliance;
            let fd=(cp-cm)/2e-4;assert!((fd-a.objective.gradient[i]).abs()<5e-4*fd.abs());
        }
    }
}
#[test]
fn g1_optimize_with_no_box_clamp_then_estimate_surface_goal_on_enriched_space() {
    let op=build(1);assert!(!op.fixed().iter().any(|x|*x));let f=forces(&op);
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());let rho=vec![0.5;study.cells()];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=controlled_sdf3_optimality_criteria(&mut study,&loads(&f),&rho,MultiLoadOcOptions{max_iterations:3,..Default::default()},&mut c);
    assert!(report.history.len()>1,"{report:?}");assert!(report.history.last().unwrap().compliance<report.history[0].compliance);
    assert!(report.history.iter().all(|r|r.volume_fraction<=0.50000001));
    let goal=study.estimate_reference_compliance_enrichment(build(2),&[
        GoalReferenceLoad3{load:ReferenceLoad3::traction(&x),weight:0.5},GoalReferenceLoad3{load:ReferenceLoad3::traction(&z),weight:0.5}],&report.displacements,
        GoalRefinementOptions3{preconditioner:GoalPreconditioner3::Jacobi{max_contributions:100_000_000},..Default::default()},&mut c).unwrap();
    close(&[goal.coarse_value],&[report.history.last().unwrap().compliance],1e-8);
    assert_eq!(goal.cases.len(),2);assert!(goal.cases.iter().all(|g|g.identity_relative_defect<1e-8));
    assert!(!goal.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
}
#[test]
fn g4_interrupted_accepted_step_preserves_the_embedded_model_and_fields() {
    let op=build(1);let f=forces(&op);let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());let rho=vec![0.5;study.cells()];
    let mut poll=|p:SolveProgress|if p.stage=="sdf3-accept"{ControlFlow::Break(())}else{ControlFlow::Continue(())};let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let r=controlled_sdf3_optimality_criteria(&mut study,&loads(&f),&rho,MultiLoadOcOptions{max_iterations:3,..Default::default()},&mut c);
    assert_eq!(r.termination,MultiLoadOcTermination::Cancelled);assert_eq!(r.history.len(),1);assert_eq!(r.rho,rho);
    let before=study.operator().elasticity().scales().to_vec();let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let again=study.evaluate(&r.rho,&loads(&f),&mut c).unwrap();assert_eq!(again.objective.displacements,r.displacements);assert_eq!(study.operator().elasticity().scales(),before);
    assert_eq!(study.operator().elasticity().embedded_dirichlet_penalty(),Some(32.0));
}
#[test]
fn g5_repeated_embedded_density_studies_replay_exactly() {
    let op=build(1);let f=forces(&op);let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());let rho=vec![0.5;study.cells()];
    let run=|s:&mut CutDensityStudy3<AdaptiveSolveSpace3>|{let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
        controlled_sdf3_optimality_criteria(s,&loads(&f),&rho,MultiLoadOcOptions{max_iterations:2,..Default::default()},&mut c)};
    let a=run(&mut study);let b=run(&mut study);assert!(a.history.len()>1);assert_eq!(a.rho,b.rho);assert_eq!(a.projected_rho,b.projected_rho);assert_eq!(a.displacements,b.displacements);assert_eq!(a.work,b.work);
}
