//! Real mixed-load optimization, independent enrichment and selective refinement.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::{ReferenceLoad3,SurfaceForce3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_dwr::elasticity3::GoalError3;
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop,MultiLoadOcOptions,SimpParams,SolveBudget,SolveControl,SolveProgress};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3,Sdf3Elasticity,controlled_sdf3_optimality_criteria,inherit_raw_densities3};
use fs_topopt::sdf3_goal::{GoalBodyLoad3,GoalReferenceLoad3,GoalRefinementOptions3,GoalPreconditioner3,GoalRefinementError3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval{let d=if a==HeightAxis::Z{1.0}else{0.0};Interval::new(d,d)}
}
fn tree(l:u8)->Octree3{Octree3::uniform(l,4,4096).unwrap()}
fn build(t:&Octree3,surface:bool)->AdaptiveElasticity3{
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut p).unwrap();
    let d=HexCell::try_new([0.0;3],[1.0;3]).unwrap();let m=IsotropicElastic::new(1.0,0.3,1.0).unwrap();
    if surface{AdaptiveElasticity3::build_with_surface(d,t,&Slab,&m,&|p|p[0]==0.0,ElasticityOptions3::default(),Default::default(),&mut q).unwrap()}
    else{AdaptiveElasticity3::build(d,t,&Slab,&m,&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()}
}
fn body(_:[f64;3])->[f64;3]{[0.0,-0.15,0.0]}
fn pressure(p:[f64;3])->f64{1.0+0.2*p[0]}
fn negative(p:[f64;3])->f64{-pressure(p)}
fn shear(p:[f64;3],_:[f64;3])->[f64;3]{[0.0,1.0+0.2*p[1],0.0]}
fn laws()->[GoalReferenceLoad3<'static>;3]{[
    GoalReferenceLoad3{load:ReferenceLoad3::pressure(&pressure),weight:0.2},
    GoalReferenceLoad3{load:ReferenceLoad3::pressure(&negative),weight:0.3},
    GoalReferenceLoad3{load:ReferenceLoad3{body:Some(&body),surface:Some(SurfaceForce3::Traction(&shear))},weight:0.5},
]}
fn forces(op:&AdaptiveElasticity3,ls:&[GoalReferenceLoad3<'_>])->Vec<Vec<f64>>{
    ls.iter().map(|l|op.reference_load(l.load,||ControlFlow::Continue(())).unwrap()).collect()
}
fn loads<'a>(f:&'a [Vec<f64>],ls:&[GoalReferenceLoad3<'_>])->Vec<LoadCase<'a>>{
    f.iter().zip(ls).map(|(f,l)|LoadCase{force:f,weight:l.weight}).collect()
}
fn study(t:&Octree3)->CutDensityStudy3<AdaptiveSolveSpace3>{
    CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(build(t,true),100_000_000),0.15,SimpParams::default())
}
#[test]
fn g3_reference_families_use_all_solver_policies_without_load_cancellation(){
    let mut s=study(&tree(1));let ls=laws();let f=forces(s.operator().elasticity(),&ls);let rho=vec![0.5;s.cells()];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=s.evaluate(&rho,&loads(&f,&ls),&mut c).unwrap();let scales=s.operator().scales().to_vec();
    let coarse=build(&tree(0),false);let mut reference=None;
    for policy in [GoalPreconditioner3::Identity,GoalPreconditioner3::Jacobi{max_contributions:100_000_000},
        GoalPreconditioner3::TwoLevel{budget:Default::default(),max_diagonal_contributions:100_000_000},GoalPreconditioner3::Multilevel{options:Default::default()}]{
        let extra=if matches!(policy,GoalPreconditioner3::Multilevel{..}){vec![&coarse]}else{Vec::new()};
        let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        let report=s.estimate_reference_compliance_enrichment_with_coarse_levels(build(&tree(2),true),&extra,&ls,&accepted.objective.displacements,
            GoalRefinementOptions3{preconditioner:policy,..Default::default()},&mut c).unwrap();
        assert_eq!(report.work.linear_solves,3);assert_eq!(report.cases.len(),3);assert!(report.fine_value>0.0);
        assert!((report.cases[0].fine_value-report.cases[1].fine_value).abs()<1e-8*report.cases[0].fine_value);
        assert!((report.coarse_value-accepted.objective.compliance).abs()<1e-8*report.coarse_value);
        if let Some(value)=reference{assert!((report.fine_value-value).abs()<1e-7*report.fine_value);}reference=Some(report.fine_value);
        assert_eq!(s.operator().scales(),scales);assert!(!report.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
    }
}
#[test]
fn g4_stale_second_field_is_refused_before_setup_or_a_fine_solve(){
    let mut s=study(&tree(1));let ls=laws();let f=forces(s.operator().elasticity(),&ls);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let mut a=s.evaluate(&vec![0.5;s.cells()],&loads(&f,&ls),&mut c).unwrap();a.objective.displacements[1].fill(0.0);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);let before=c.work();
    let result=s.estimate_reference_compliance_enrichment(build(&tree(2),true),&ls,&a.objective.displacements,
        GoalRefinementOptions3{preconditioner:GoalPreconditioner3::TwoLevel{budget:Default::default(),max_diagonal_contributions:100_000_000},..Default::default()},&mut c);
    assert!(matches!(result,Err(GoalRefinementError3::Estimate(GoalError3::FieldResidual{field:"coarse-primal",..}))));assert_eq!(c.work(),before);
}
#[test]
fn g4_interrupted_second_load_and_exhausted_budget_return_no_estimate_family(){
    let mut s=study(&tree(1));let ls=laws();let f=forces(s.operator().elasticity(),&ls);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=s.evaluate(&vec![0.5;s.cells()],&loads(&f,&ls),&mut c).unwrap();let scales=s.operator().scales().to_vec();
    let mut p=|s:SolveProgress|if s.stage=="sdf3-goal-elasticity"&&s.work.linear_solves>=2&&s.solve_iterations>0{ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(matches!(s.estimate_reference_compliance_enrichment(build(&tree(2),true),&ls,&a.objective.displacements,Default::default(),&mut c),Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(c.work().linear_solves,2);assert!(c.work().linear_iterations>0);assert_eq!(s.operator().scales(),scales);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:1,..Default::default()},&mut p);
    assert!(matches!(s.estimate_reference_compliance_enrichment(build(&tree(2),true),&ls,&a.objective.displacements,Default::default(),&mut c),Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget{..}))));
    assert_eq!(c.work().linear_iterations,1);assert_eq!(s.operator().scales(),scales);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(matches!(s.estimate_reference_compliance_enrichment(build(&tree(2),false),&ls,&a.objective.displacements,Default::default(),&mut c),Err(GoalRefinementError3::Physics(_))));
}
#[test]
fn g1_pressure_optimization_marks_original_cells_then_restores_refined_volume(){
    let t=tree(1);let mut s=study(&t);let ls=laws();let f=forces(s.operator().elasticity(),&ls);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let rho=vec![0.5;s.cells()];let report=controlled_sdf3_optimality_criteria(&mut s,&loads(&f,&ls),&rho,MultiLoadOcOptions{max_iterations:2,..Default::default()},&mut c);
    assert!(report.history.len()>1,"{report:?}");
    for pair in report.history.windows(2){assert!(pair[1].compliance<=pair[0].compliance);assert!(pair[1].volume_fraction<=0.50000001);}
    let goal=s.estimate_reference_compliance_enrichment(build(&tree(2),true),&ls,&report.displacements,Default::default(),&mut c).unwrap();
    let mark=goal.mark(0.5,2,||ControlFlow::Continue(())).unwrap();assert!(!mark.marked.is_empty());assert!(mark.marked.iter().all(|leaf|t.leaves().contains(leaf)));
    let next=t.refined(&mark.marked,||ControlFlow::Continue(())).unwrap();assert_eq!(t.leaves().len(),8);let mut refined=study(&next);
    let raw=inherit_raw_densities3(s.operator().elasticity().leaves(),&report.rho,refined.operator().elasticity().leaves(),&mut c).unwrap();
    let rho=refined.feasible_start(&raw,0.5,1e-8,&mut c).unwrap();let f=forces(refined.operator().elasticity(),&ls);
    let a=refined.evaluate(&rho,&loads(&f,&ls),&mut c).unwrap();assert!(a.volume_fraction<=0.50000001);assert_eq!(a.objective.displacements.len(),3);
}
#[test]
fn g5_old_body_api_matches_reference_load_family_and_work(){
    let mut s=study(&tree(1));let ls=[GoalReferenceLoad3{load:ReferenceLoad3::body(&body),weight:1.0}];let f=forces(s.operator().elasticity(),&ls);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=s.evaluate(&vec![0.5;s.cells()],&loads(&f,&ls),&mut c).unwrap();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let old=s.estimate_compliance_enrichment(build(&tree(2),true),&[GoalBodyLoad3{density:&body,weight:1.0}],&a.objective.displacements,Default::default(),&mut c).unwrap();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let new=s.estimate_reference_compliance_enrichment(build(&tree(2),true),&ls,&a.objective.displacements,Default::default(),&mut c).unwrap();
    assert_eq!(old.fine_value.to_bits(),new.fine_value.to_bits());assert_eq!(old.marking_mass,new.marking_mass);assert_eq!(old.work,new.work);
}
