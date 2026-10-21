use super::*;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::ElasticityOptions3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SolveBudget, SolveProgress};
use crate::sdf3::{PhysicalRegion3, checked_solve};
use crate::sdf3_goal::GoalPreconditioner3;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { p[2]-0.73 }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)
    }
    fn derivative_enclose(&self, _: [f64;3], _: [f64;3], a: HeightAxis) -> Interval {
        let d=if a==HeightAxis::Z {1.0} else {0.0}; Interval::new(d,d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(), &|p| p[0]==0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn params() -> SimpParams { SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() } }
fn y(p: [f64;3]) -> [f64;3] { [0.0,-1.0-p[1],0.0] }
fn yz(p: [f64;3]) -> [f64;3] { [0.0,-0.4*(1.0+p[1]),-1.0-p[0]] }
fn laws() -> [ReferenceLoad3<'static>;2] { [ReferenceLoad3::body(&y), ReferenceLoad3::body(&yz)] }
fn options() -> GoalRefinementOptions3 {
    GoalRefinementOptions3 { preconditioner: GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 }, ..Default::default() }
}
fn trees() -> (Octree3, Octree3) {
    let tree=Octree3::uniform(1,4,4096).unwrap();
    let fine=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(), || ControlFlow::Continue(())).unwrap();
    (tree,fine)
}
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 1e-8*a.abs().max(b.abs()).max(1e-12), "{a} != {b}"); }

#[test]
fn fine_worst_combination_replays_inherited_eroded_material_without_raw_refiltering() {
    let (tree,fine_tree)=trees();
    let mut study=CutDensityStudy3::new(build(&tree),0.15,SimpParams::default());
    let raw: Vec<_>=(0..study.cells()).map(|i|0.35+0.04*i as f64).collect();
    let original=study.params(); let scales=study.operator().scales().to_vec();
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let robust=RobustOptions3::default();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=study.estimate_load_space_cap_enrichment(build(&fine_tree),&laws(),&raw,params(),robust,
        &cap,options(),&mut control).unwrap();
    assert_eq!(report.basis_assessment.cases.len(),2);
    assert_eq!(report.coarse.spectrum.gram.len(),4);
    assert!(report.fine.worst_displacement.len()>report.coarse.worst_displacement.len());
    assert!(report.stable_and_feasible(1e-6,1.0));
    assert_eq!(study.operator().scales(),scales);
    assert_eq!(study.params().eta,original.eta); assert_eq!(study.params().beta,original.beta);
    // Reconstruct the SAME physical fine operator independently. No fine raw
    // control/filter evaluation occurs here; every child inherits its parent E.
    let shifted=SimpParams { eta: params().eta+robust.threshold_offset, ..params() };
    let mut reference=CutDensityStudy3::new(build(&tree),0.15,shifted);
    let fy=reference.operator().reference_load(laws()[0],||ControlFlow::Continue(())).unwrap();
    let fz=reference.operator().reference_load(laws()[1],||ControlFlow::Continue(())).unwrap();
    reference.evaluate(&raw,&[LoadCase { force:&fy,weight:0.5 },LoadCase { force:&fz,weight:0.5 }],&mut control).unwrap();
    let mut fine=build(&fine_tree);
    let inherited=AdaptiveTransfer3::new(reference.operator(),&fine,options().max_transfer_terms,
        ||ControlFlow::Continue(())).unwrap().inherited_scales();
    fine.set_scales(&inherited).unwrap();
    let direct=checked_solve(&fine,&report.fine.worst_force,1e-10,"fine-worst-replay",&mut control).unwrap();
    let work: f64=direct.iter().zip(&report.fine.worst_force).map(|(u,f)|u*f).sum();
    close(work,report.fine.worst_compliance);
    let error: f64=direct.iter().zip(&report.fine.worst_displacement).map(|(a,b)|(a-b)*(a-b)).sum();
    assert!(error <= 1e-15*direct.iter().map(|x|x*x).sum::<f64>());
    for (i,case) in report.basis_assessment.cases.iter().enumerate() {
        close(case.coarse_value/cap.cap(),report.coarse.spectrum.gram[i*2+i]);
        close(case.fine_value/cap.cap(),report.fine.spectrum.gram[i*2+i]);
    }
}

#[test]
fn individually_passing_refined_basis_loads_do_not_pass_the_combined_cap() {
    let (tree,fine)=trees();
    let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.55;study.cells()];
    let laws=[ReferenceLoad3::body(&y),ReferenceLoad3::body(&y)];
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let first=study.estimate_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control).unwrap();
    let individual=first.basis_assessment.cases[0].fine_value;
    let cap=LoadSpaceCap3::new(2,1.25*individual,0.02).unwrap();
    let report=study.estimate_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control).unwrap();
    assert!(report.basis_assessment.cases.iter().all(|case|case.fine_value<cap.cap()));
    close(report.fine.worst_compliance,2.0*individual);
    assert!(!report.stable_and_feasible(0.9,f64::MAX));
}

#[test]
fn equal_maxima_cannot_hide_rotating_adverse_directions_or_off_diagonal_change() {
    let a=[0.8,0.0,0.0,0.1]; let b=[0.1,0.0,0.0,0.8];
    let cap=LoadSpaceCap3::new(2,1.0,0.02).unwrap();
    close(spectrum(a.to_vec(),cap.policy()).unwrap().maximum,spectrum(b.to_vec(),cap.policy()).unwrap().maximum);
    close(gram_change(&a,&b).unwrap(),0.7*2.0_f64.sqrt());
    close(gram_change(&a,&[0.8,0.2,0.2,0.1]).unwrap(),0.2*2.0_f64.sqrt());
    assert!(gram_change(&a,&[f64::NAN;4]).is_err());
    assert!(gram_change(&a,&[]).is_err());
}

#[test]
fn invalid_or_unenriched_requests_refuse_before_any_load_callback() {
    let (tree,fine)=trees();
    let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.5;study.cells()];
    let forbidden=|_: [f64;3]| -> [f64;3] { panic!("refused request invoked load law") };
    let laws=[ReferenceLoad3::body(&forbidden),ReferenceLoad3::body(&forbidden)];
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    assert!(study.estimate_load_space_cap_enrichment(build(&tree),&laws,&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control).is_err());
    let bad=GoalRefinementOptions3 { max_load_cases:1,..options() };
    assert!(study.estimate_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),RobustOptions3::default(),
        &cap,bad,&mut control).is_err());
    let bad=GoalRefinementOptions3 { max_transfer_terms:1,..options() };
    assert!(study.estimate_load_space_cap_enrichment(build(&fine),&laws,&raw,params(),RobustOptions3::default(),
        &cap,bad,&mut control).is_err());
    assert_eq!(control.work().linear_solves,0);
}

#[test]
fn late_cancellation_and_shared_budget_restore_the_incoming_nonunit_state() {
    let (tree,fine)=trees();
    let mut study=CutDensityStudy3::new(build(&tree),0.15,params());
    let raw=vec![0.5;study.cells()];
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let force=study.operator().reference_load(laws()[0],||ControlFlow::Continue(())).unwrap();
    study.evaluate(&raw,&[LoadCase { force:&force,weight:1.0 }],&mut control).unwrap();
    let scales=study.operator().scales().to_vec(); let original=study.params();
    let mut poll=|p:SolveProgress| if p.stage=="sdf3-load-space-grid-publish" {
        ControlFlow::Break(())
    } else {ControlFlow::Continue(())};
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let result=study.estimate_load_space_cap_enrichment(build(&fine),&laws(),&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control);
    assert!(matches!(result,Err(GoalRefinementError3::Evaluation(EvaluationStop::Cancelled))));
    assert!(control.work().linear_solves>=5);
    assert_eq!(study.operator().scales(),scales);
    assert_eq!(study.params().eta,original.eta);
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget { total_iterations:1,..Default::default() },&mut poll);
    let result=study.estimate_load_space_cap_enrichment(build(&fine),&laws(),&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control);
    assert!(matches!(result,Err(GoalRefinementError3::Evaluation(EvaluationStop::TotalBudget { .. }))));
    assert_eq!(control.work().linear_iterations,1); assert_eq!(study.operator().scales(),scales);
}

#[test]
fn protected_regions_and_original_reference_estimator_match_without_extra_state_solves() {
    let (tree,fine)=trees();
    let study=CutDensityStudy3::new(build(&tree),0.15,params());
    let mut regions=vec![PhysicalRegion3::Design;study.cells()];
    regions[0]=PhysicalRegion3::Solid;regions[1]=PhysicalRegion3::Void;
    let mut study=study.with_physical_regions(regions).unwrap();
    let raw=vec![0.55;study.cells()];let scales=study.operator().scales().to_vec();
    let cap=LoadSpaceCap3::new(2,1e9,0.02).unwrap();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let report=study.estimate_load_space_cap_enrichment(build(&fine),&laws(),&raw,params(),RobustOptions3::default(),
        &cap,options(),&mut control).unwrap();
    assert_eq!(report.eroded.projected_rho[0],1.0);assert_eq!(report.eroded.projected_rho[1],0.0);
    let loads: Vec<_>=laws().iter().map(|&load|GoalReferenceLoad3 {load,weight:0.5}).collect();
    let mut poll=|_|ControlFlow::Continue(());
    let mut control=SolveControl::new(SolveBudget::default(),&mut poll);
    let old=study.estimate_robust_reference_compliance_enrichment(build(&fine),&loads,&raw,params(),
        RobustOptions3::default(),options(),&mut control).unwrap();
    assert_eq!(old.estimate.fine_value,report.basis_assessment.fine_value);
    assert_eq!(old.estimate.marking_mass,report.basis_assessment.marking_mass);
    assert_eq!(old.eroded.objective.displacements,report.eroded.objective.displacements);
    assert_eq!(old.estimate.work,report.work);
    assert_eq!(study.operator().scales(),scales);
}
