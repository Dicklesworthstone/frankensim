//! G1/G3/G4/G5: response fitting's ACTUAL objective, including imposed motion.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{ResponseCase3, ResponseOptions3, ResponseTarget3};
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3, ReferenceResponseTarget3, ResponseRefinementOptions3};
use fs_topopt::sdf3_goal::GoalRefinementError3;

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
fn build(level:u8)->AdaptiveElasticity3 {
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&Octree3::uniform(level,4,4096).unwrap(),&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|_|false,&|_,n|n[0]<0.0,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn force(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[-0.03,0.0,-0.02]}else{[0.0;3]}}
fn observation(p:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.0,0.0,1.0+p[1]]}else{[0.0;3]}}
fn motion(p:[f64;3],_:[f64;3])->[f64;3] {[0.004*p[1],0.007*p[2],0.003]}
fn target(value:f64)->ReferenceResponseTarget3<'static> {
    ReferenceResponseTarget3{observation:ReferenceLoad3::traction(&observation),target:value,scale:1.0,weight:1.0}
}

#[test]
fn g1_exact_coarse_fit_still_has_refinement_signal_and_a_complete_quadratic_identity() {
    let op=build(1);let force_vector=op.reference_load(ReferenceLoad3::traction(&force),||ControlFlow::Continue(())).unwrap();
    let q=op.reference_load(ReferenceLoad3::traction(&observation),||ControlFlow::Continue(())).unwrap();
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
    let rho=vec![0.5;study.cells()];let mut poll=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let initial=[ResponseTarget3{q:&q,target:0.0,scale:1.0,weight:1.0}];
    let first=study.evaluate_responses(&rho,&[ResponseCase3{force:&force_vector,prescribed:Some(&motion),targets:&initial}],Default::default(),&mut control).unwrap();
    let exact=[ResponseTarget3{target:first.responses[0][0],..initial[0]}];
    let accepted=study.evaluate_responses(&rho,&[ResponseCase3{force:&force_vector,prescribed:Some(&motion),targets:&exact}],Default::default(),&mut control).unwrap();
    assert_eq!(accepted.objective,0.0);assert!(accepted.adjoints[0].iter().all(|v|*v==0.0));
    let targets=[target(exact[0].target)];let laws=[ReferenceResponseCase3{load:ReferenceLoad3::traction(&force),prescribed:Some(&motion),targets:&targets}];
    let before=study.operator().elasticity().scales().to_vec();
    let result=study.estimate_response_enrichment(AdaptiveSolveSpace3::jacobi(build(2),100_000_000),&accepted,&laws,Default::default(),&mut control).unwrap();
    assert_eq!(result.coarse_objective,0.0);assert!(result.fine_objective>1e-12,"{result:?}");
    assert!(result.cases[0].secant_weights.iter().any(|w|*w!=0.0));
    assert!(result.identity_relative_defect<1e-8);
    assert!((result.correction()-result.fine_objective).abs()<1e-8*result.fine_objective);
    assert!(!result.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
    assert_eq!(study.operator().elasticity().scales(),before);
}

#[test]
fn g3_all_existing_fine_preparations_agree_without_mutating_the_coarse_design() {
    let op=build(1);let f=op.reference_load(ReferenceLoad3::traction(&force),||ControlFlow::Continue(())).unwrap();
    let q=op.reference_load(ReferenceLoad3::traction(&observation),||ControlFlow::Continue(())).unwrap();
    let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());
    let targets=[ResponseTarget3{q:&q,target:0.1,scale:0.8,weight:0.7}];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let options=ResponseOptions3{volume_weight:0.05,..Default::default()};
    let accepted=study.evaluate_responses(&vec![0.5;study.cells()],&[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&targets}],options,&mut c).unwrap();
    let t=[ReferenceResponseTarget3{target:0.1,scale:0.8,weight:0.7,..target(0.1)}];
    let laws=[ReferenceResponseCase3{load:ReferenceLoad3::traction(&force),prescribed:Some(&motion),targets:&t}];
    let settings=ResponseRefinementOptions3{response:options,..Default::default()};
    let reference=study.estimate_response_enrichment(build(2),&accepted,&laws,settings,&mut c).unwrap();
    for mode in 0..3 {
        let fine=build(2);let coarse=build(1);
        let prepared=match mode {
            0=>AdaptiveSolveSpace3::jacobi(fine,100_000_000),
            1=>AdaptiveSolveSpace3::two_level(fine,&coarse,Default::default(),||ControlFlow::Continue(())).unwrap(),
            _=>AdaptiveSolveSpace3::multilevel(fine,&[&coarse],Default::default(),||ControlFlow::Continue(())).unwrap(),
        };
        let actual=study.estimate_response_enrichment(prepared,&accepted,&laws,settings,&mut c).unwrap();
        assert!((actual.fine_objective-reference.fine_objective).abs()<1e-8*reference.fine_objective.abs());
        assert!(actual.identity_relative_defect<1e-8);
    }
}

#[test]
fn g4_second_case_failure_is_rejected_before_setup_and_cancellation_restores_state() {
    let op=build(1);let f=op.reference_load(ReferenceLoad3::traction(&force),||ControlFlow::Continue(())).unwrap();
    let q=op.reference_load(ReferenceLoad3::traction(&observation),||ControlFlow::Continue(())).unwrap();
    let mut study=CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op,100_000_000),0.15,SimpParams::default());
    let t=[ResponseTarget3{q:&q,target:0.1,scale:1.0,weight:1.0}];
    let case=ResponseCase3{force:&f,prescribed:Some(&motion),targets:&t};
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=study.evaluate_responses(&vec![0.5;study.cells()],&[case,case],Default::default(),&mut c).unwrap();
    let targets=[target(0.1)];let law=ReferenceResponseCase3{load:ReferenceLoad3::traction(&force),prescribed:Some(&motion),targets:&targets};
    let before=study.operator().elasticity().scales().to_vec();let mut bad=accepted.clone();bad.displacements[1][0]+=0.02;
    let mut preparations=0;let mut p=|s:SolveProgress|{if s.stage=="sdf3-preconditioner-start"{preparations+=1;}ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(study.estimate_response_enrichment(build(2),&bad,&[law,law],Default::default(),&mut c).is_err());
    drop(c);assert_eq!(preparations,0);assert_eq!(study.operator().elasticity().scales(),before);
    let mut p=|s:SolveProgress|if s.stage=="response-goal-fine-adjoint"&&s.solve_iterations>0{ControlFlow::Break(())}else{ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    assert!(matches!(study.estimate_response_enrichment(build(2),&accepted,&[law,law],Default::default(),&mut c),
        Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert!(c.work().linear_iterations>0);assert_eq!(study.operator().elasticity().scales(),before);
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:1,..Default::default()},&mut p);
    assert!(study.estimate_response_enrichment(build(2),&accepted,&[law,law],Default::default(),&mut c).is_err());
    assert_eq!(c.work().linear_iterations,1);assert_eq!(study.operator().elasticity().scales(),before);
}

#[test]
fn g5_completed_response_goal_replays_while_corrupted_density_evidence_refuses() {
    let op=build(1);let f=op.reference_load(ReferenceLoad3::traction(&force),||ControlFlow::Continue(())).unwrap();
    let q=op.reference_load(ReferenceLoad3::traction(&observation),||ControlFlow::Continue(())).unwrap();
    let mut study=CutDensityStudy3::new(op,0.15,SimpParams::default());let t=[ResponseTarget3{q:&q,target:0.0,scale:1.0,weight:1.0}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=study.evaluate_responses(&vec![0.5;study.cells()],&[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&t}],Default::default(),&mut c).unwrap();
    let targets=[target(0.0)];let laws=[ReferenceResponseCase3{load:ReferenceLoad3::traction(&force),prescribed:Some(&motion),targets:&targets}];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let a=study.estimate_response_enrichment(build(2),&accepted,&laws,Default::default(),&mut c).unwrap();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let b=study.estimate_response_enrichment(build(2),&accepted,&laws,Default::default(),&mut c).unwrap();
    assert_eq!(a.fine_objective.to_bits(),b.fine_objective.to_bits());assert_eq!(a.work,b.work);
    assert_eq!(a.cases[0].fine_displacement,b.cases[0].fine_displacement);
    let mut bad=accepted;bad.projected_rho[0]+=0.1;
    assert!(study.estimate_response_enrichment(build(2),&bad,&laws,Default::default(),&mut c).is_err());
}
