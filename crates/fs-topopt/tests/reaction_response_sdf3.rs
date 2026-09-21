//! Real mixed observation losses: direct reaction terms plus one adjoint/case.
use std::{cell::Cell,ops::ControlFlow};
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3,ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams,SolveBudget,SolveControl,SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{ReactionTarget3,ResponseCase3,ResponseTarget3,ProjectedResponseStudy3,ProjectedResponseOptions3};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64{(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval{
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],a:HeightAxis)->Interval{
        if a==HeightAxis::X{Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}else{Interval::new(0.0,0.0)}
    }
}
fn motion(p:[f64;3],_:[f64;3])->[f64;3]{[0.02*(p[0]-0.17),0.0,0.0]}
fn right(_:[f64;3],n:[f64;3])->[f64;3]{if n[0]>0.0{[1.0,0.0,0.0]}else{[0.0;3]}}
fn backend()->(AdaptiveSolveSpace3,Vec<f64>,Vec<f64>){
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let tree=Octree3::uniform(1,4,4096).unwrap();
    let tree=tree.refined(&[*tree.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap();
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&tree,&Slab,
        &IsotropicElastic::new(1.0,0.0,1.0).unwrap(),&|_|false,&|_,_|true,ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    let force=op.body_load(&|_|[0.0;3],||ControlFlow::Continue(())).unwrap();
    let observation=op.body_load(&|p|[1.0+p[1],0.0,0.0],||ControlFlow::Continue(())).unwrap();
    (AdaptiveSolveSpace3::jacobi(op,100_000_000),force,observation)
}
fn fixture()->(CutDensityStudy3<AdaptiveSolveSpace3>,Vec<f64>,Vec<f64>){
    let (op,f,q)=backend();
    (CutDensityStudy3::new(op,0.15,SimpParams{beta:2.0,..Default::default()}),f,q)
}
#[test]
fn g3_mixed_targets_have_full_chain_coordinate_derivatives_and_one_preparation(){
    let (mut study,f,q)=fixture();let d=[ResponseTarget3{q:&q,target:0.004,scale:0.02,weight:0.3}];
    let r=[ReactionTarget3{mode:&right,target:0.001,scale:0.02,weight:0.7}];let rr:[&[ReactionTarget3<'_>];2]=[&r,&r];
    let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&d};2];
    let rho:Vec<f64>=(0..study.cells()).map(|i|0.35+0.01*i as f64).collect();let incoming=study.operator().elasticity().scales().to_vec();
    let mut preparations=0;let mut poll=|p:SolveProgress|{if p.stage=="sdf3-preconditioner-start"{preparations+=1;}ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    let base=study.evaluate_responses_with_reactions(&rho,&cases,&rr,Default::default(),&mut c).unwrap();drop(c);assert_eq!(preparations,1);
    assert_eq!(base.responses.len(),2);assert_eq!(base.reaction_responses.len(),2);assert!(base.reaction_responses[0][0]>0.0);
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    for i in 0..rho.len(){
        let h=1e-4;let mut a=rho.clone();let mut b=rho.clone();a[i]+=h;b[i]-=h;
        let plus=study.evaluate_responses_with_reactions(&a,&cases,&rr,Default::default(),&mut c).unwrap();
        let minus=study.evaluate_responses_with_reactions(&b,&cases,&rr,Default::default(),&mut c).unwrap();let fd=(plus.objective-minus.objective)/(2.0*h);
        assert!((fd-base.gradient[i]).abs()<2e-7+2e-4*fd.abs(),"{i}: {fd} vs {}",base.gradient[i]);
    }
    assert_eq!(study.operator().elasticity().scales(),incoming);
}
#[test]
fn g5_empty_reaction_rows_preserve_displacement_only_results_and_work_exactly(){
    let (mut study,f,q)=fixture();let d=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:1.0}];
    let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&d}];let rho=vec![0.4;study.cells()];
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let a=study.evaluate_responses(&rho,&cases,Default::default(),&mut ca).unwrap();
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let b=study.evaluate_responses_with_reactions(&rho,&cases,&[&[]],Default::default(),&mut cb).unwrap();
    assert_eq!(a.objective,b.objective);assert_eq!(a.gradient,b.gradient);assert_eq!(a.displacements,b.displacements);assert_eq!(a.adjoints,b.adjoints);
    assert_eq!(a.work,b.work);assert_eq!(a.reaction_responses,b.reaction_responses);
}
#[test]
fn g1_reaction_only_fit_installs_complete_fields_and_preserves_split_run_state(){
    let (mut study,f,_)=fixture();let (mut split,_,_)=fixture();let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&[]}];
    let r=[ReactionTarget3{mode:&right,target:0.001,scale:0.002,weight:1.0}];let rr:[&[ReactionTarget3<'_>];1]=[&r];let rho=vec![0.5;study.cells()];
    let options=ProjectedResponseOptions3{volume_cap:0.6,..Default::default()};
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let mut a=ProjectedResponseStudy3::new_with_reactions(&mut study,&cases,&rr,&rho,options,&mut ca).unwrap();
    let mut b=ProjectedResponseStudy3::new_with_reactions(&mut split,&cases,&rr,&rho,options,&mut cb).unwrap();
    let initial=a.accepted().objective;let ra=a.run(8).unwrap();b.run(3).unwrap();let rb=b.run(5).unwrap();
    assert!(a.accepted().objective<initial);assert!(a.constraint_violation()<=options.optimizer.tolerance);
    assert_eq!(ra.stop,rb.stop);assert_eq!(a.point(),b.point());assert_eq!(a.accepted().reaction_responses,b.accepted().reaction_responses);
    assert_eq!(a.work(),b.work());assert_eq!(a.optimizer_work(),b.optimizer_work());
    let accepted=a.accepted().clone();drop(a);
    let again=study.evaluate_responses_with_reactions(&accepted.rho,&cases,&rr,options.response,&mut ca).unwrap();
    assert_eq!(accepted.displacements,again.displacements);assert_eq!(accepted.gradient,again.gradient);assert_eq!(accepted.reaction_responses,again.reaction_responses);
}
#[test]
fn g4_reaction_and_adjoint_cancellation_restore_accepted_scales_and_resume(){
    for stage in ["sdf3-response-reaction","sdf3-response-adjoint"] {
        let (mut study,f,_)=fixture();let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&[]}];
        let targets=[ReactionTarget3{mode:&right,target:0.001,scale:0.002,weight:1.0}];let rr:[&[ReactionTarget3<'_>];1]=[&targets];
        let enabled=Cell::new(false);let mut poll=|p:SolveProgress|if enabled.get()&&p.stage==stage{ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut poll);let rho=vec![0.5;study.cells()];
        let mut session=ProjectedResponseStudy3::new_with_reactions(&mut study,&cases,&rr,&rho,Default::default(),&mut c).unwrap();
        let before=session.study().operator().elasticity().scales().to_vec();let work=session.work();enabled.set(true);assert!(session.run(1).is_err());enabled.set(false);
        assert_eq!(session.point(),rho);assert_eq!(session.study().operator().elasticity().scales(),before);assert!(session.work().linear_iterations>work.linear_iterations);
        session.run(1).unwrap();assert_eq!(session.accepted().reaction_responses[0].len(),1);
    }
}
#[test]
fn g0_admission_and_unsupported_reaction_refinement_fail_before_physics(){
    let (mut study,f,q)=fixture();let d=[ResponseTarget3{q:&q,target:0.004,scale:0.02,weight:0.3}];let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&d}];
    let r=[ReactionTarget3{mode:&right,target:0.001,scale:0.02,weight:0.7}];let rr:[&[ReactionTarget3<'_>];1]=[&r];let rho=vec![0.4;study.cells()];
    let mut poll=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut poll);
    assert!(study.evaluate_responses_with_reactions(&rho,&cases,&[],Default::default(),&mut c).is_err());assert_eq!(c.work().linear_solves,0);
    let accepted=study.evaluate_responses_with_reactions(&rho,&cases,&rr,Default::default(),&mut c).unwrap();let spent=c.work();let (fine,_,_)=backend();
    let result=study.estimate_response_enrichment(fine,&accepted,&[],Default::default(),&mut c);
    assert!(result.is_err());assert_eq!(c.work(),spent);
}

#[path = "reaction_response_sdf3/refinement.rs"]
mod refinement;
