use super::*;
use crate::sdf3::robust::envelope::load_space::LoadSpace3;
use crate::sdf3::PhysicalRegion3;
use crate::{SolveBudget, SolveProgress};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { p[2]-0.73 }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)
    }
    fn derivative_enclose(&self, _: [f64;3], _: [f64;3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d,d)
    }
}
fn build(tree: &Octree3) -> AdaptiveElasticity3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth:1, ..Default::default() }, &mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q).unwrap()
}
fn y(p: [f64;3]) -> [f64;3] { [0.0,-1.0-p[1],0.0] }
fn z(p: [f64;3]) -> [f64;3] { [0.0,-0.4*(1.0+p[1]),-1.0-p[0]] }
fn laws() -> [ReferenceLoad3<'static>;2] { [ReferenceLoad3::body(&y),ReferenceLoad3::body(&z)] }
fn params() -> SimpParams { SimpParams { penal:2.0,beta:2.0,e_min:1e-3,..Default::default() } }
fn policy(steps: usize) -> LoadSpaceCapRefitOptions3 {
    let mut o = LoadSpaceCapRefitOptions3 { accepted_steps:steps,..Default::default() };
    o.start.max_bisections = 6; o
}
fn fixture() -> (Octree3,Octree3,CutDensityStudy3<AdaptiveElasticity3>,Vec<f64>) {
    let tree = Octree3::uniform(1,4,4096).unwrap();
    let op = build(&tree);
    let fine = tree.refined(&[op.leaves()[0]],||ControlFlow::Continue(())).unwrap();
    let study = CutDensityStudy3::new(op,0.15,params());
    let rho = vec![0.75;study.cells()];
    (tree,fine,study,rho)
}
fn cap_for(tree: &Octree3, control: &mut SolveControl<'_>) -> LoadSpaceCap3 {
    let mut reference = CutDensityStudy3::new(build(tree),0.15,params());
    let forces: Vec<_> = laws().iter().map(|law|
        reference.operator().reference_load(*law,||ControlFlow::Continue(())).unwrap()).collect();
    let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let raw = vec![0.75;reference.cells()];
    let e = reference.evaluate_robust_load_space(&raw,&basis,policy(0).design.robustness,
        &LoadSpace3::new(2,1.0,0.02).unwrap(),control).unwrap();
    LoadSpaceCap3::new(2,1.8*e.worst_compliance,0.02).unwrap()
}

#[test]
fn selective_refit_descends_and_replays_the_full_combined_load_problem() {
    let (_,fine,source,rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let cap = cap_for(&fine,&mut control);
    let result = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),&cap,policy(8),&mut control).unwrap();
    assert!(result.completed_and_rechecked(), "{:?} {:?}",result.outcome,result.recheck);
    assert!(result.study.cells() > source.cells());
    assert_eq!(result.cap,cap);
    assert_eq!(result.study.filter_radius,source.filter_radius);
    assert!(result.history.len() > 1);
    let best = result.best.as_ref().unwrap();
    assert!(cap.admits(&best.evaluation,policy(8).design.optimizer.tolerance));
    assert!(best.evaluation.physical.dilated_volume_fraction < result.history[0].dilated_volume);
    assert_eq!(best.evaluation.worst_displacement,
        result.recheck.as_ref().unwrap().as_ref().unwrap().evaluation.worst_displacement);
    // An independent fine owner verifies loads, cross terms and paired fields.
    let mut direct = CutDensityStudy3::new(build(&fine),0.15,params());
    let forces: Vec<_> = laws().iter().map(|law|
        direct.operator().reference_load(*law,||ControlFlow::Continue(())).unwrap()).collect();
    let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let fresh = direct.evaluate_robust_load_space(&best.rho,&basis,policy(8).design.robustness,cap.policy(),&mut control).unwrap();
    assert_eq!(fresh.spectrum.gram,best.evaluation.spectrum.gram);
    assert_eq!(fresh.worst_displacement,best.evaluation.worst_displacement);
    assert_eq!(fresh.physical.dilated_rho,best.evaluation.physical.dilated_rho);
    assert_eq!(source.operator().scales(),scales);
    assert_eq!(source.params().beta,params().beta);
}

#[test]
fn raw_parent_map_preserves_tiny_controls_and_authored_solid_void_regions() {
    let (_,fine,source,mut rho) = fixture();
    rho[0] = 2e-5;
    let mut regions = vec![PhysicalRegion3::Design;source.cells()];
    regions[0] = PhysicalRegion3::Solid; regions[1] = PhysicalRegion3::Void;
    let source = source.with_physical_regions(regions.clone()).unwrap();
    let operator = build(&fine);
    let parents = AdaptiveTransfer3::new(source.operator(),&operator,2_000_000,
        ||ControlFlow::Continue(())).unwrap().parents().to_vec();
    let mut options = policy(0); options.design.density_floor = 1e-5;
    let cap = LoadSpaceCap3::new(2,1e12,0.02).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let result = source.refit_load_space_cap(operator,&rho,params(),&laws(),&cap,options,&mut control).unwrap();
    assert!(result.completed_and_rechecked());
    let point = result.best.as_ref().unwrap();
    for (i,&parent) in parents.iter().enumerate() {
        assert_eq!(result.inherited_rho[i],rho[parent]);
        assert_eq!(result.study.physical_regions().unwrap()[i],regions[parent]);
        for field in [&point.evaluation.physical.eroded.projected_rho,
            &point.evaluation.physical.nominal_rho,&point.evaluation.physical.dilated_rho] {
            match regions[parent] {
                PhysicalRegion3::Solid => assert_eq!(field[i],1.0),
                PhysicalRegion3::Void => assert_eq!(field[i],0.0),
                PhysicalRegion3::Design => {}
            }
        }
    }
}

#[test]
fn initial_and_refined_owners_share_the_same_numerical_path() {
    let (_,fine,source,rho) = fixture();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let cap = cap_for(&fine,&mut control);
    let refit = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),&cap,policy(4),&mut control).unwrap();
    let initial = CutDensityStudy3::new(build(&fine),0.15,params())
        .solve_reference_load_space_cap(&refit.inherited_rho,&laws(),&cap,policy(4),&mut control).unwrap();
    assert!(initial.completed_and_rechecked());
    assert_eq!(initial.history,refit.history);
    assert_eq!(initial.best.as_ref().unwrap().rho,refit.best.as_ref().unwrap().rho);
    assert_eq!(initial.best.unwrap().evaluation.spectrum.gram,refit.best.unwrap().evaluation.spectrum.gram);
}

#[test]
fn unchanged_space_and_invalid_policy_refuse_before_load_callbacks() {
    let (tree,fine,source,rho) = fixture();
    let forbidden = |_: [f64;3]| -> [f64;3] { panic!("invalid refit invoked a load law") };
    let laws = [ReferenceLoad3::body(&forbidden);2];
    let cap = LoadSpaceCap3::new(2,1e12,0.02).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    assert!(source.refit_load_space_cap(build(&tree),&rho,params(),&laws,&cap,policy(0),&mut control).is_err());
    let mut options = policy(0); options.max_transfer_terms = 1;
    assert!(source.refit_load_space_cap(build(&fine),&rho,params(),&laws,&cap,options,&mut control).is_err());
    options = policy(0); options.max_load_cases = 1;
    assert!(source.refit_load_space_cap(build(&fine),&rho,params(),&laws,&cap,options,&mut control).is_err());
    assert!(source.refit_load_space_cap(build(&fine),&rho,params(),&laws[..1],&cap,policy(0),&mut control).is_err());
    assert_eq!(control.work().linear_solves,0);
}

#[test]
fn accepted_and_recheck_cancellation_preserve_distinct_solved_prefixes() {
    for stage in ["sdf3-load-space-design-accepted","sdf3-load-space-cap-recheck"] {
        let (_,fine,source,rho) = fixture();
        let scales = source.operator().scales().to_vec();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        let cap = cap_for(&fine,&mut control);
        let mut poll = |p: SolveProgress| if p.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        let result = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),&cap,policy(4),&mut control).unwrap();
        assert!(!result.completed_and_rechecked());
        assert!(result.restoration.point.is_some() && result.best.is_some());
        if stage == "sdf3-load-space-design-accepted" {
            assert_eq!(result.history.len(),2);
            assert!(result.outcome.as_ref().unwrap().is_err()); assert!(result.recheck.is_none());
        } else {
            assert!(result.outcome.as_ref().unwrap().is_ok());
            assert!(matches!(result.recheck,Some(Err(EvaluationStop::Cancelled))));
        }
        assert_eq!(source.operator().scales(),scales);
    }
}

#[test]
fn exhaustion_and_rejected_full_density_do_not_claim_a_completed_refit() {
    let (_,fine,source,rho) = fixture();
    let scales = source.operator().scales().to_vec();
    let cap = LoadSpaceCap3::new(2,1e12,0.02).unwrap();
    let mut options = policy(4); options.design.optimizer.max_evaluations = 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let limited = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),&cap,options,&mut control).unwrap();
    assert!(!limited.completed_and_rechecked());
    assert_eq!(limited.outcome.as_ref().unwrap().as_ref().unwrap().stop,ProjectedAlStop::EvaluationLimit);
    assert!(limited.best.is_some() && limited.recheck.is_none());
    let refused = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),
        &LoadSpaceCap3::new(2,1e-20,0.02).unwrap(),policy(0),&mut control).unwrap();
    assert_eq!(refused.restoration.stop,LoadSpaceStartStop3::FullDensityRejected);
    assert!(refused.best.is_none() && refused.outcome.is_none());
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations:1,..Default::default() },&mut poll);
    let stopped = source.refit_load_space_cap(build(&fine),&rho,params(),&laws(),&cap,policy(4),&mut control).unwrap();
    assert_eq!(stopped.work.linear_iterations,1);
    assert!(stopped.restoration.point.is_none() && stopped.best.is_none());
    assert!(!stopped.completed_and_rechecked());
    assert_eq!(source.operator().scales(),scales);
}
