use super::*;
use std::cell::Cell;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityError3, ElasticityOptions3, adaptive::{AdaptiveElasticity3,
    enrichment::precondition::AdaptiveSolveSpace3}};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use crate::{SimpParams, SolveBudget, SolveProgress};
use crate::sdf3::CutDensityStudy3;
use crate::sdf3_goal::GoalPreconditioner3;
use super::super::{LoadSpaceCap3, refit::LoadSpaceCapRefitOptions3};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
fn build(tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>) -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
    let mut poll = |_| cp();
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll)
        .map_err(ElasticityError3::from)?;
    let op = AdaptiveElasticity3::build(HexCell::try_new([0.0;3], [1.0;3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q)?;
    Ok(AdaptiveSolveSpace3::jacobi(op, 100_000_000))
}
fn y(p: [f64;3]) -> [f64;3] { [0.0, -1.0-p[1], 0.0] }
fn yz(p: [f64;3]) -> [f64;3] { [0.0, -0.6*(1.0+p[1]), -1.0-p[0]] }
fn laws() -> [ReferenceLoad3<'static>;2] { [ReferenceLoad3::body(&y), ReferenceLoad3::body(&yz)] }
fn initial(control: &mut SolveControl<'_>) -> (Octree3, LoadSpaceCapRefit3<AdaptiveSolveSpace3>) {
    let tree = Octree3::uniform(0,4,4096).unwrap();
    let params = SimpParams { penal:2.0, beta:2.0, e_min:1e-3, ..Default::default() };
    let study = CutDensityStudy3::new(build(&tree, &mut || ControlFlow::Continue(())).unwrap(),0.15,params);
    let raw = vec![0.3; study.cells()];
    let options = LoadSpaceCapRefitOptions3 { accepted_steps:0, ..Default::default() };
    let result = study.solve_reference_load_space_cap(&raw,&laws(),
        &LoadSpaceCap3::new(2,1e12,0.03).unwrap(),options,control).unwrap();
    assert!(result.completed_and_rechecked());
    (tree,result)
}
fn options(max_refits: usize) -> LoadSpaceCapCampaignOptions3 {
    LoadSpaceCapCampaignOptions3 { max_refits, change_tolerance:0.0,
        enrichment:GoalRefinementOptions3 { preconditioner:GoalPreconditioner3::Jacobi {
            max_contributions:100_000_000 }, ..Default::default() }, ..Default::default() }
}

#[test]
fn segmented_refits_match_uninterrupted_grids_fields_and_cumulative_work() {
    let run = |segments: &[usize]| {
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        let (tree,current) = initial(&mut control); let laws = laws(); let builds = Cell::new(0);
        let builder = |t: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>| {
            builds.set(builds.get()+1); build(t,cp)
        };
        let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,options(2),&mut control,builder).unwrap();
        let before = campaign.work();
        assert_eq!(campaign.run(0).unwrap(),LoadSpaceCapCampaignStop3::Paused);
        assert_eq!(campaign.work(),before); assert_eq!(builds.get(),0);
        for (index,&steps) in segments.iter().enumerate() {
            let stop = campaign.run(steps).unwrap();
            assert_eq!(stop,if index+1 == segments.len() {LoadSpaceCapCampaignStop3::RefinementLimit}
                else {LoadSpaceCapCampaignStop3::Paused});
        }
        assert_eq!(campaign.grids().len(),3);
        assert_eq!(campaign.probes().len(),3); // includes the final installed grid
        assert_eq!(campaign.candidate_attempts(),2); assert_eq!(builds.get(),5);
        let point = campaign.current().recheck.as_ref().unwrap().as_ref().unwrap();
        assert!(campaign.current().cap.admits(&point.evaluation,1e-6));
        let saved = (campaign.grids().to_vec(),campaign.probes().to_vec(),point.rho.clone(),
            point.evaluation.worst_displacement.clone(),point.evaluation.spectrum.gram.clone(),campaign.work());
        assert_eq!(campaign.run(10).unwrap(),LoadSpaceCapCampaignStop3::RefinementLimit);
        assert_eq!(campaign.work(),saved.5); assert_eq!(builds.get(),5);
        saved
    };
    assert_eq!(run(&[2]),run(&[1,1]));
}

#[test]
fn cancelled_promotion_keeps_the_prepared_candidate_for_retry_without_reoptimization() {
    let stopped = Cell::new(false);
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-campaign-promote" && !stopped.replace(true) {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let (tree,current) = initial(&mut control); let laws = laws(); let builds = Cell::new(0);
    let builder = |t: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>| {
        builds.set(builds.get()+1); build(t,cp)
    };
    let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,options(1),&mut control,builder).unwrap();
    assert!(campaign.run(1).is_err());
    assert_eq!(campaign.grids().len(),1); assert_eq!(builds.get(),2);
    let retained = campaign.pending().unwrap().1.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::RefinementLimit);
    assert_eq!(campaign.candidate_attempts(),1);
    assert_eq!(builds.get(),3); // only the new current-grid probe was built
    assert_eq!(campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho,retained);
    assert!(campaign.pending().is_none());
}

#[test]
fn cancellation_after_promotion_keeps_the_new_grid_and_retries_only_its_probe() {
    let stopped = Cell::new(false);
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-load-space-campaign-promoted" && !stopped.replace(true) {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let (tree,current) = initial(&mut control); let laws = laws();
    let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,options(1),&mut control,build).unwrap();
    assert!(campaign.run(1).is_err());
    assert_eq!(campaign.grids().len(),2); assert_eq!(campaign.probes().len(),1);
    assert!(campaign.pending().is_none()); assert!(campaign.current_probe().is_none());
    let raw = campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    assert_eq!(raw.len(),campaign.current().study.cells());
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::RefinementLimit);
    assert_eq!(campaign.candidate_attempts(),1); assert_eq!(campaign.probes().len(),2);
    assert_eq!(campaign.probes()[1].source_grid,1);
    assert_eq!(campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho,raw);
}

#[test]
fn failed_candidate_build_is_charged_and_does_not_repeat_the_completed_probe() {
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let (tree,current) = initial(&mut control); let laws = laws(); let builds = Cell::new(0);
    let builder = |t: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>| {
        builds.set(builds.get()+1);
        if builds.get() == 2 { return Err(GoalRefinementError3::Invalid("injected candidate build failure")); }
        build(t,cp)
    };
    let mut policy = options(2); policy.max_candidate_attempts = 1;
    let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,policy,&mut control,builder).unwrap();
    assert!(campaign.run(2).is_err());
    assert_eq!(campaign.grids().len(),1); assert!(campaign.current_probe().is_some());
    let work = campaign.work();
    assert_eq!(campaign.run(2).unwrap(),LoadSpaceCapCampaignStop3::CandidateAttemptLimit);
    assert_eq!(campaign.candidate_attempts(),1); assert_eq!(campaign.work(),work); assert_eq!(builds.get(),2);
}

#[test]
fn physical_work_exhaustion_cannot_be_renewed_by_a_new_segment_or_grid() {
    let budget = {
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        let (tree,current) = initial(&mut control); let laws = laws();
        let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,options(0),&mut control,build).unwrap();
        assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::RefinementLimit);
        campaign.work().linear_iterations+1
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {total_iterations:budget,..Default::default()},&mut poll);
    let (tree,current) = initial(&mut control); let laws = laws();
    let mut policy = options(2); policy.max_candidate_attempts = 1;
    let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,policy,&mut control,build).unwrap();
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::CandidateStopped);
    assert_eq!(campaign.work().linear_iterations,budget);
    assert_eq!(campaign.grids().len(),1); assert!(campaign.rejected().is_some());
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::CandidateAttemptLimit);
    assert_eq!(campaign.work().linear_iterations,budget);
}

#[test]
fn current_grid_pass_is_checked_without_refits_and_bad_policy_refuses_before_builds() {
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
    let (tree,current) = initial(&mut control); let laws = laws();
    let mut policy = options(0); policy.change_tolerance = 1.0;
    let mut campaign = LoadSpaceCapCampaign3::new(current,tree,&laws,policy,&mut control,build).unwrap();
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::StableAndFeasible);
    assert_eq!(campaign.grids().len(),1); assert_eq!(campaign.probes().len(),1);
    assert_eq!(campaign.candidate_attempts(),0);
    let work = campaign.work();
    assert_eq!(campaign.run(1).unwrap(),LoadSpaceCapCampaignStop3::StableAndFeasible);
    assert_eq!(campaign.work(),work);
    drop(campaign);
    let (tree,current) = initial(&mut control);
    policy.change_tolerance = f64::NAN;
    let forbidden = |_: &Octree3, _: &mut dyn FnMut() -> ControlFlow<()>| -> Result<AdaptiveSolveSpace3,GoalRefinementError3> {
        panic!("invalid policy reached geometry builder")
    };
    assert!(LoadSpaceCapCampaign3::new(current,tree,&laws,policy,&mut control,forbidden).is_err());
}
