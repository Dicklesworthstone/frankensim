use super::*;
use super::super::tests::{build, trees, params, laws, options, close};
use crate::{SolveBudget, SolveProgress};
use crate::sdf3::{PhysicalRegion3, checked_solve};

#[test]
fn fine_binary_witness_replays_exact_parent_stiffness_and_protected_masks() {
    let (tree,fine_tree)=trees();
    let study=CutDensityStudy3::new(build(&tree),0.15,SimpParams::default());
    let mut regions=vec![PhysicalRegion3::Design;study.cells()];
    regions[0]=PhysicalRegion3::Solid;regions[1]=PhysicalRegion3::Void;
    let mut study=study.with_physical_regions(regions).unwrap();
    let original=study.params();let scales=study.operator().scales().to_vec();
    let raw:Vec<_>=(0..study.cells()).map(|i|0.35+0.05*i as f64).collect();
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let binary=BinaryComplianceOptions3 { threshold:0.45,void_stiffness:0.005,..Default::default() };
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=study.estimate_binary_load_space_cap_enrichment(build(&fine_tree),&laws(),&raw,params(),
        &cap,binary,options(),&mut control).unwrap();
    assert_eq!(report.work.linear_solves,5); // one source filter, two coarse and two fine states
    for mask in [&report.physical.eroded,&report.physical.nominal,&report.physical.dilated] {
        assert!(mask.solid[0]);assert!(!mask.solid[1]);
    }
    assert_eq!(study.operator().scales(),scales);
    assert_eq!(study.params().beta,original.beta);assert_eq!(study.params().eta,original.eta);
    assert_eq!(report.physical.params.beta,params().beta);
    let mut coarse=build(&tree);coarse.set_scales(&report.physical.eroded_stiffness_scales()).unwrap();
    let mut fine=build(&fine_tree);
    let transfer=AdaptiveTransfer3::new(&coarse,&fine,options().max_transfer_terms,
        ||ControlFlow::Continue(())).unwrap();
    let inherited=transfer.inherited_scales();drop(transfer);
    assert!(inherited.iter().all(|v|*v==1.0 || *v==binary.void_stiffness));
    fine.set_scales(&inherited).unwrap();
    let u=checked_solve(&fine,&report.fine.worst_force,1e-10,"binary-fine-replay",&mut control).unwrap();
    close(u.iter().zip(&report.fine.worst_force).map(|(u,f)|u*f).sum(),report.fine.worst_compliance);
    let error:f64=u.iter().zip(&report.fine.worst_displacement).map(|(a,b)|(a-b)*(a-b)).sum();
    assert!(error<1e-15*u.iter().map(|u|u*u).sum::<f64>());
    assert!(report.stable_and_feasible(1e-6,1.0));
}

#[test]
fn binary_enrichment_checks_combinations_instead_of_passing_basis_cases() {
    let (tree,fine)=trees();let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.8;study.cells()];let law=laws()[0];let laws=[law,law];
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let endpoint=BinaryComplianceOptions3::default();
    let first=study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),
        &LoadSpaceCap3::new(2,1e9,0.02).unwrap(),endpoint,options(),&mut control).unwrap();
    let individual=first.basis_assessment.cases[0].fine_value;
    let cap=LoadSpaceCap3::new(2,1.25*individual,0.02).unwrap();
    let report=study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),
        &cap,endpoint,options(),&mut control).unwrap();
    assert!(report.basis_assessment.cases.iter().all(|e|e.fine_value<cap.cap()));
    close(report.fine.worst_compliance,2.0*individual);
    assert!(!report.stable_and_feasible(0.9,f64::MAX));
}

#[test]
fn binary_refinement_restores_state_after_fine_cancellation_and_work_exhaustion() {
    let (tree,fine)=trees();let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.6;study.cells()];let scales=study.operator().scales().to_vec();
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll=|p:SolveProgress| if p.stage=="sdf3-goal-elasticity" && p.solve_iterations>0 {
        ControlFlow::Break(())
    } else {ControlFlow::Continue(())};
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let result=study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws(),&raw,params(),
        &cap,BinaryComplianceOptions3::default(),options(),&mut control);
    assert!(matches!(result,Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(study.operator().scales(),scales);assert!(control.work().linear_iterations>0);
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget {total_iterations:1,..Default::default()},&mut poll);
    let result=study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws(),&raw,params(),
        &cap,BinaryComplianceOptions3::default(),options(),&mut control);
    assert!(matches!(result,Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget {..}))));
    assert_eq!(control.work().linear_iterations,1);assert_eq!(study.operator().scales(),scales);
}

#[test]
fn invalid_binary_options_refuse_before_reference_load_callbacks() {
    let (tree,fine)=trees();let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.6;study.cells()];let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let forbidden=|_:[f64;3]|->[f64;3]{panic!("invalid binary request invoked a load")};
    let laws=[ReferenceLoad3::body(&forbidden),ReferenceLoad3::body(&forbidden)];
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    for binary in [BinaryComplianceOptions3 {threshold:f64::NAN,..Default::default()},
        BinaryComplianceOptions3 {void_stiffness:0.0,..Default::default()},
        BinaryComplianceOptions3 {max_load_cases:1,..Default::default()}] {
        assert!(study.estimate_binary_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),
            &cap,binary,options(),&mut control).is_err());
    }
    assert!(study.estimate_binary_load_space_cap_enrichment(build(&tree),&laws,&raw,params(),
        &cap,BinaryComplianceOptions3::default(),options(),&mut control).is_err());
    assert_eq!(control.work().linear_solves,0);
}
