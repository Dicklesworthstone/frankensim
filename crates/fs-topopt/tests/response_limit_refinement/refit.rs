use super::*;
use fs_topopt::sdf3::response::limits::design::*;
use fs_topopt::sdf3::PhysicalRegion3;

#[test]
fn reference_minimum_material_uses_the_same_nodal_optimizer_and_fields() {
    let mut a=study(false);let mut b=study(false);let laws=laws();let rho=raw(&a);
    let mut pa=|_|ControlFlow::Continue(());let mut ca=SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb=|_|ControlFlow::Continue(());let mut cb=SolveControl::new(SolveBudget::default(),&mut pb);
    let actual=a.fit_reference_minimum_material(&cases(&laws),&rho,Default::default(),3,&mut ca).unwrap();
    let f=b.operator().elasticity().reference_load(ReferenceLoad3::body(&body),||ControlFlow::Continue(())).unwrap();
    let q=b.operator().elasticity().reference_load(ReferenceLoad3::body(&observe),||ControlFlow::Continue(())).unwrap();
    let rows:Vec<_>=laws.iter().map(|l|ResponseLimit3 {sense:l.sense,bound:l.bound,scale:l.scale,response:match l.response {
        ReferenceBoundedResponse3::Displacement(_)=>BoundedResponse3::Displacement(&q),
        ReferenceBoundedResponse3::Reaction(mode)=>BoundedResponse3::Reaction(mode),
    }}).collect();
    let nodal=[ResponseLimitCase3 {force:&f,prescribed:Some(&motion),limits:&rows}];
    let mut expected=MinimumMaterialStudy3::new(&mut b,&nodal,&rho,Default::default(),&mut cb).unwrap();
    let result=expected.run(3);assert_eq!(actual.outcome.is_ok(),result.is_ok());
    assert_eq!(actual.history,expected.history());assert_eq!(actual.accepted.rho,expected.accepted().rho);
    assert_eq!(actual.accepted.measurements,expected.accepted().measurements);
    assert_eq!(actual.accepted.constraint_gradient,expected.accepted().constraint_gradient);
    assert_eq!(actual.accepted.displacements,expected.accepted().displacements);assert_eq!(actual.work,expected.work());
}

#[test]
fn marked_refit_transfers_only_raw_controls_and_solves_new_constraints() {
    let mut s=study(false);let laws=laws();let mut p=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);let source=evaluate(&mut s,&laws,&mut c);
    let before=s.operator().elasticity().scales().to_vec();
    let e=s.estimate_response_limit_enrichment(fine(),&source,&cases(&laws),Default::default(),&mut c).unwrap();
    let marked=e.mark(0.5,2,||ControlFlow::Continue(())).unwrap();assert!(!marked.marked.is_empty());
    let t=tree(false).refined(&marked.marked,||ControlFlow::Continue(())).unwrap();
    let candidate=CutDensityStudy3::new(build(&t),0.15,s.params());
    let parents=AdaptiveTransfer3::new(s.operator().elasticity(),candidate.operator().elasticity(),2_000_000,||ControlFlow::Continue(())).unwrap().parents().to_vec();
    let mut result=s.refit_reference_minimum_material(candidate,&source,&cases(&laws),Default::default(),0,2_000_000,&mut c).unwrap();
    assert_eq!(result.inherited_rho,parents.iter().map(|&i|source.rho[i]).collect::<Vec<_>>());
    assert_eq!(result.fit.initial.rho,result.inherited_rho);assert_eq!(result.fit.history[0].iteration,0);
    assert_eq!(result.fit.history.len(),1);assert!(result.study.cells()>s.cells());
    assert!(result.fit.accepted.displacements[0].len()>source.displacements[0].len());
    let repeated=result.study.fit_reference_minimum_material(&cases(&laws),&result.inherited_rho,Default::default(),0,&mut c).unwrap();
    assert_eq!(result.fit.initial.measurements,repeated.initial.measurements);
    assert_eq!(result.fit.initial.displacements,repeated.initial.displacements);
    assert_eq!(s.operator().elasticity().scales(),before);
    assert_eq!(result.fit.best_feasible.is_some(),result.fit.initial.feasible(MinimumMaterialOptions3::default().optimizer.tolerance));
}

#[test]
fn interrupted_candidate_retains_its_baseline_and_never_replaces_source() {
    let mut s=study(false);let laws=laws();let mut p=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);let source=evaluate(&mut s,&laws,&mut c);
    let before=s.operator().elasticity().scales().to_vec();let initialized=Cell::new(false);
    let mut p=|p:SolveProgress| {
        if p.stage=="response-limit-design-initialize" {initialized.set(true);}
        if initialized.get()&&p.stage=="response-limit-design-control" {ControlFlow::Break(())}else{ControlFlow::Continue(())}
    };
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);
    let result=s.refit_reference_minimum_material(CutDensityStudy3::new(fine(),0.15,s.params()),&source,&cases(&laws),Default::default(),3,2_000_000,&mut c).unwrap();
    assert!(result.fit.outcome.is_err());assert_eq!(result.fit.accepted.rho,result.fit.initial.rho);
    assert_eq!(result.fit.accepted.displacements,result.fit.initial.displacements);assert_eq!(result.fit.history.len(),1);
    assert!(result.fit.work.linear_iterations>0);assert_eq!(s.operator().elasticity().scales(),before);
}

#[test]
fn refit_inherits_protected_regions_and_rejects_stale_designs() {
    let s=study(false);let mut labels=vec![PhysicalRegion3::Design;s.cells()];labels[0]=PhysicalRegion3::Solid;labels[1]=PhysicalRegion3::Void;
    let mut s=s.with_physical_regions(labels.clone()).unwrap();let laws=laws();let mut p=|_|ControlFlow::Continue(());
    let mut c=SolveControl::new(SolveBudget::default(),&mut p);let source=evaluate(&mut s,&laws,&mut c);
    let op=fine();let parents=AdaptiveTransfer3::new(s.operator().elasticity(),op.elasticity(),2_000_000,||ControlFlow::Continue(())).unwrap().parents().to_vec();
    let result=s.refit_reference_minimum_material(CutDensityStudy3::new(op,0.15,s.params()),&source,&cases(&laws),Default::default(),0,2_000_000,&mut c).unwrap();
    assert_eq!(result.study.physical_regions().unwrap(),parents.iter().map(|&i|labels[i]).collect::<Vec<_>>());
    for (i,&p) in parents.iter().enumerate() {match labels[p] {
        PhysicalRegion3::Solid=>assert_eq!(result.fit.initial.projected_rho[i],1.0),
        PhysicalRegion3::Void=>assert_eq!(result.fit.initial.projected_rho[i],0.0),PhysicalRegion3::Design=>{},
    }}
    let mut stale=source.clone();stale.rho[2]=0.8;
    assert!(s.refit_reference_minimum_material(CutDensityStudy3::new(fine(),0.15,s.params()),&stale,&cases(&laws),Default::default(),0,2_000_000,&mut c).is_err());
}
