use super::*;
use crate::sdf3::robust::envelope::limits::campaign::{ComplianceCapCampaign3,
    ComplianceCapCampaignOptions3, ComplianceCapCampaignStop3 as Stop};
use crate::sdf3::robust::envelope::limits::refit::ComplianceCapRefitOptions3;
use std::cell::Cell;

fn policies(refits: usize) -> (ComplianceCapRefitOptions3, ComplianceCapCampaignOptions3) {
    let design = ComplianceCapRefitOptions3 { accepted_steps: 1, ..Default::default() };
    let campaign = ComplianceCapCampaignOptions3 { max_refits: refits,
        max_marks: 1, change_tolerance: 0.0, ..Default::default() };
    (design, campaign)
}
fn builder(tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>)
    -> Result<AdaptiveElasticity3, GoalRefinementError3> {
    let mut poll = |_| cp();
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll)
        .map_err(fs_cutfem::elastic3::ElasticityError3::from)?;
    Ok(AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
        ElasticityOptions3::default(), &mut q)?)
}
fn requirements() -> ComplianceCaps3 { ComplianceCaps3::new(vec![1e6, 2e6], 0.02, 2).unwrap() }

#[test]
fn split_refits_match_uninterrupted_grids_fields_probes_and_work() {
    let run = |chunks: &[usize]| {
        let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
        let (design, options) = policies(2);
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
        let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, builder).unwrap();
        let before = campaign.work();
        assert_eq!(campaign.run(0).unwrap(), Stop::Paused);
        assert_eq!(campaign.work(), before); assert!(campaign.probes().is_empty());
        let mut stop = Stop::Paused;
        for &n in chunks { stop = campaign.run(n).unwrap(); }
        assert_eq!(stop, Stop::RefinementLimit);
        assert_eq!(campaign.grids().len(), 3); assert_eq!(campaign.probes().len(), 3);
        assert_eq!(campaign.current().caps, caps);
        assert_eq!(campaign.probes().last().unwrap().source_grid, 2);
        let point = campaign.current().recheck.as_ref().unwrap().as_ref().unwrap();
        assert!(caps.satisfied_by(&point.evaluation.physical.eroded.objective.case_compliances));
        assert_eq!(point.rho.len(), campaign.current().study.cells());
        (campaign.grids().to_vec(), campaign.probes().to_vec(), point.rho.clone(),
            point.evaluation.physical.eroded.objective.displacements.clone(), campaign.work(), campaign.candidate_attempts())
    };
    assert_eq!(run(&[2]), run(&[1, 1]));
}

#[test]
fn publication_retry_reuses_the_complete_candidate_without_rebuilding() {
    let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
    let (design, options) = policies(1);
    let cancel = Cell::new(true); let builds = Cell::new(0);
    let mut poll = |p: SolveProgress| if cancel.get() && p.stage == "sdf3-cap-campaign-promote" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
    let original = initial.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    let build = |t: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>| {
        builds.set(builds.get()+1); builder(t, cp)
    };
    let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, build).unwrap();
    assert!(campaign.run(1).is_err());
    assert_eq!(campaign.grids().len(), 1); assert_eq!(builds.get(), 2);
    assert_eq!(campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho, original);
    assert!(campaign.pending().unwrap().1.completed_and_rechecked());
    let pending = campaign.pending().unwrap().1.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    cancel.set(false);
    assert_eq!(campaign.run(1).unwrap(), Stop::RefinementLimit);
    assert_eq!(builds.get(), 3, "retry should build only the new design's diagnostic grid");
    assert_eq!(campaign.candidate_attempts(), 1);
    assert!(campaign.pending().is_none());
    assert_eq!(campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho, pending);
}

#[test]
fn cancellation_after_promotion_keeps_and_then_assesses_the_new_design() {
    let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
    let (design, options) = policies(1); let cancel = Cell::new(true);
    let mut poll = |p: SolveProgress| if cancel.get() && p.stage == "sdf3-cap-campaign-promoted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
    let initial_cells = initial.study.cells();
    let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, builder).unwrap();
    assert!(campaign.run(1).is_err());
    assert_eq!(campaign.grids().len(), 2); assert!(campaign.current_probe().is_none());
    assert!(campaign.current().study.cells() > initial_cells);
    let latest = campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    cancel.set(false);
    assert_eq!(campaign.run(1).unwrap(), Stop::RefinementLimit);
    assert_eq!(campaign.probes().last().unwrap().source_grid, 1);
    assert_eq!(campaign.current().recheck.as_ref().unwrap().as_ref().unwrap().rho, latest);
}

#[test]
fn failed_candidates_consume_attempts_while_the_source_probe_stays_cached() {
    let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
    let (design, mut options) = policies(2); options.max_candidate_attempts = 2;
    let builds = Cell::new(0);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
    let build = |t: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>| {
        builds.set(builds.get()+1);
        if builds.get() > 1 { return Err(GoalRefinementError3::Invalid("candidate geometry refused")); }
        builder(t, cp)
    };
    let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, build).unwrap();
    assert!(campaign.run(1).is_err()); let work = campaign.work();
    assert!(campaign.run(1).is_err());
    assert_eq!(campaign.run(1).unwrap(), Stop::CandidateAttemptLimit);
    assert_eq!(campaign.candidate_attempts(), 2); assert_eq!(builds.get(), 3);
    assert_eq!(campaign.probes().len(), 1); assert_eq!(campaign.grids().len(), 1);
    assert_eq!(campaign.work(), work);
}

#[test]
fn only_a_current_two_grid_check_passes_and_terminal_calls_spend_no_more_work() {
    for tolerance in [0.0, 1.0] {
        let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
        let (design, mut options) = policies(0); options.change_tolerance = tolerance;
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
        let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, builder).unwrap();
        let expected = if tolerance == 0.0 { Stop::RefinementLimit } else { Stop::StableAndFeasible };
        assert_eq!(campaign.run(1).unwrap(), expected);
        assert_eq!(campaign.probes().len(), 1); assert_eq!(campaign.grids().len(), 1);
        let work = campaign.work();
        assert_eq!(campaign.run(100).unwrap(), expected); assert_eq!(campaign.work(), work);
    }
}

#[test]
fn campaign_refuses_renewed_physical_work_ledger_before_building() {
    let (tree, _, study, rho) = fixture(); let laws = laws(); let caps = requirements();
    let (design, options) = policies(1);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let initial = study.solve_reference_compliance_caps(&rho, &laws, &caps, design, &mut control).unwrap();
    drop(control);
    let mut replacement = SolveControl::new(SolveBudget::default(), &mut poll);
    let forbidden = |_: &Octree3, _: &mut dyn FnMut() -> ControlFlow<()>| -> Result<AdaptiveElasticity3, GoalRefinementError3> {
        panic!("renewed work allowance reached builder")
    };
    assert!(ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut replacement, forbidden).is_err());
    assert_eq!(replacement.work().linear_solves, 0);
}
