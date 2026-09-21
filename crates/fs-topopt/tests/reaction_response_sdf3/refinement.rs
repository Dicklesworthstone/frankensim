//! The same reaction-fitting fixtures, now through real enriched goal solves.
use super::*;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3,ReferenceResponseTarget3};
use fs_topopt::sdf3_goal::GoalRefinementError3;
use fs_topopt::EvaluationStop;
fn zero(_:[f64;3])->[f64;3]{[0.0;3]}
fn observe(p:[f64;3])->[f64;3]{[1.0+p[1],0.0,0.0]}
fn enriched()->AdaptiveSolveSpace3 {
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),
        &Octree3::uniform(2,4,4096).unwrap(),&Slab,&IsotropicElastic::new(1.0,0.0,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op,100_000_000)
}
#[test]
fn mixed_loss_reconstructs_reintegrated_reactions_and_uses_one_preparation_per_grid() {
    let(mut study,f,q)=fixture();let d=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:0.3}];
    let r=[ReactionTarget3{mode:&right,target:0.001,scale:0.002,weight:0.7}];let rr:[&[ReactionTarget3<'_>];2]=[&r,&r];
    let nodal=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&d};2];
    let rho:Vec<_>=(0..study.cells()).map(|i|0.25+0.025*i as f64).collect();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=study.evaluate_responses_with_reactions(&rho,&nodal,&rr,Default::default(),&mut c).unwrap();
    let target=[ReferenceResponseTarget3{observation:ReferenceLoad3::body(&observe),target:d[0].target,scale:d[0].scale,weight:d[0].weight}];
    let cases=[ReferenceResponseCase3{load:ReferenceLoad3::body(&zero),prescribed:Some(&motion),targets:&target};2];
    let incoming=study.operator().elasticity().scales().to_vec();let mut preparations=0;
    let mut p=|s:SolveProgress|{if s.stage=="sdf3-preconditioner-start"{preparations+=1;}ControlFlow::Continue(())};
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let result=study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&cases,&rr,Default::default(),&mut c).unwrap();
    assert_eq!(result.cases.len(),2);assert!((result.coarse_objective-accepted.objective).abs()<1e-8*accepted.objective);
    assert!(result.identity_relative_defect<1e-8);
    assert!((result.correction()-(result.fine_objective-result.coarse_objective)).abs()<1e-8*result.coarse_objective.max(result.fine_objective));
    assert!(result.cases.iter().all(|c|c.coarse_reactions.len()==1&&c.fine_reactions.len()==1));
    assert!(result.cases.iter().any(|c|c.reaction_offset_correction.abs()>1e-6));
    assert!(!result.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
    assert_eq!(study.operator().elasticity().scales(),incoming);drop(c);assert_eq!(preparations,2);
}
#[test]
fn exact_coarse_reaction_fit_does_not_erase_enriched_misfit_or_refinement_signal() {
    let(mut study,f,_)=fixture();let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&[]}];
    let rho:Vec<_>=(0..study.cells()).map(|i|0.25+0.025*i as f64).collect();
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let provisional=[ReactionTarget3{mode:&right,target:0.0,scale:0.002,weight:1.0}];
    let baseline=study.evaluate_responses_with_reactions(&rho,&cases,&[&provisional],Default::default(),&mut c).unwrap();
    let target=[ReactionTarget3{target:baseline.reaction_responses[0][0],..provisional[0]}];let rr:[&[ReactionTarget3<'_>];1]=[&target];
    let exact=study.evaluate_responses_with_reactions(&rho,&cases,&rr,Default::default(),&mut c).unwrap();assert_eq!(exact.objective,0.0);
    let reference=[ReferenceResponseCase3{load:ReferenceLoad3::body(&zero),prescribed:Some(&motion),targets:&[]}];
    let report=study.estimate_response_enrichment_with_reactions(enriched(),&exact,&reference,&rr,Default::default(),&mut c).unwrap();
    let delta=(report.cases[0].fine_reactions[0]-target[0].target)/target[0].scale;
    assert!(report.fine_objective>1e-16);assert_eq!(report.coarse_objective,0.0);
    assert!((report.fine_objective-0.5*delta*delta).abs()<1e-12*report.fine_objective);
    assert!(report.cases[0].reaction_secant_weights[0]!=0.0);
    assert!(!report.mark(0.5,2,||ControlFlow::Continue(())).unwrap().marked.is_empty());
}
#[test]
fn missing_changed_or_interrupted_reaction_families_never_replace_source_fields() {
    let(mut study,f,_)=fixture();let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&[]}];
    let r=[ReactionTarget3{mode:&right,target:0.001,scale:0.002,weight:1.0}];let rr:[&[ReactionTarget3<'_>];1]=[&r];
    let rho:Vec<_>=(0..study.cells()).map(|i|0.3+0.02*i as f64).collect();let mut p=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=study.evaluate_responses_with_reactions(&rho,&cases,&rr,Default::default(),&mut c).unwrap();
    let reference=[ReferenceResponseCase3{load:ReferenceLoad3::body(&zero),prescribed:Some(&motion),targets:&[]}];
    let before=study.operator().elasticity().scales().to_vec();let spent=c.work();
    assert!(study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&reference,&[],Default::default(),&mut c).is_err());assert_eq!(c.work(),spent);
    let changed=[ReactionTarget3{target:0.2,..r[0]}];
    assert!(matches!(study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&reference,&[&changed],Default::default(),&mut c),Err(GoalRefinementError3::Invalid(_))));
    assert_eq!(study.operator().elasticity().scales(),before);
    for stage in ["response-goal-fine-adjoint","response-goal-publish"] {
        let mut p=|s:SolveProgress|if s.stage==stage{ControlFlow::Break(())}else{ControlFlow::Continue(())};
        let mut c=SolveControl::new(SolveBudget::default(),&mut p);
        assert!(matches!(study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&reference,&rr,Default::default(),&mut c),Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
        assert!(c.work().linear_iterations>0);assert_eq!(study.operator().elasticity().scales(),before);
    }
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget{total_iterations:1,..Default::default()},&mut p);
    assert!(study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&reference,&rr,Default::default(),&mut c).is_err());
    assert_eq!(c.work().linear_iterations,1);assert_eq!(study.operator().elasticity().scales(),before);
}
#[test]
fn empty_reaction_rows_preserve_displacement_estimation_values_and_work() {
    let(mut study,f,q)=fixture();let target=[ResponseTarget3{q:&q,target:0.003,scale:0.02,weight:1.0}];
    let cases=[ResponseCase3{force:&f,prescribed:Some(&motion),targets:&target}];let rho=vec![0.4;study.cells()];
    let mut p=|_|ControlFlow::Continue(());let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let accepted=study.evaluate_responses(&rho,&cases,Default::default(),&mut c).unwrap();
    let target=[ReferenceResponseTarget3{observation:ReferenceLoad3::body(&observe),target:0.003,scale:0.02,weight:1.0}];
    let reference=[ReferenceResponseCase3{load:ReferenceLoad3::body(&zero),prescribed:Some(&motion),targets:&target}];
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let a=study.estimate_response_enrichment(enriched(),&accepted,&reference,Default::default(),&mut ca).unwrap();
    let b=study.estimate_response_enrichment_with_reactions(enriched(),&accepted,&reference,&[&[]],Default::default(),&mut cb).unwrap();
    assert_eq!(a.correction().to_bits(),b.correction().to_bits());assert_eq!(a.fine_objective,b.fine_objective);assert_eq!(a.work,b.work);
    assert_eq!(b.cases[0].reaction_offset_correction,0.0);assert!(b.cases[0].fine_reactions.is_empty());
}
