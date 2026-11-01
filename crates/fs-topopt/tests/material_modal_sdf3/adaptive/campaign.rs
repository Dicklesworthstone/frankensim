//! Adaptive campaign boundaries retain complete rechecked native designs.
use super::*;
use super::refit::{assert_physical, point, seed};
use std::cell::{Cell, RefCell};
use std::rc::Rc;
use fs_topopt::sdf3::modal::campaign::{
    Sdf3EigenfrequencyBuild3, Sdf3EigenfrequencyCampaign3,
    Sdf3EigenfrequencyCampaignOptions3, Sdf3EigenfrequencyCampaignStop3,
    Sdf3EigenfrequencyMarking3,
};
use fs_topopt::sdf3::modal::refinement::Sdf3EigenfrequencyAdaptiveError3;
use fs_topopt::SolveControlCheckpoint;

fn forced_refit() -> Sdf3EigenfrequencyCampaignOptions3 {
    Sdf3EigenfrequencyCampaignOptions3 {
        max_refits: 1,
        max_candidate_attempts: 2,
        // An actual change in the independently solved spectrum forces a
        // refit, even when all other diagnostic gates admit the first probe.
        eigenvalue_tolerance: 0.0,
        objective_tolerance: 1e6,
        volume_tolerance: 1.0,
        mass_relative_tolerance: 1.0,
        subspace_tolerance: 1.0,
        marking: Sdf3EigenfrequencyMarking3::UniformActive,
        ..Default::default()
    }
}

#[test]
fn campaigns_probe_refit_and_reprobe_both_physical_goals() {
    for minimum_mass in [false, true] {
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let (tree, source) = seed(minimum_mass, &mut control);
        let initial_n = source.study.operator().n();
        let builds = Rc::new(RefCell::new(Vec::new()));
        let recorded = Rc::clone(&builds);
        let mut options = forced_refit();
        if minimum_mass {
            options.marking = Sdf3EigenfrequencyMarking3::ModalResidual { theta: 0.5, max_marks: 4 };
        }
        let mut campaign = Sdf3EigenfrequencyCampaign3::new(
            source, tree, options, &mut control,
            move |tree, kind, checkpoint| {
                recorded.borrow_mut().push(kind);
                build(tree, 3.0, checkpoint)
            },
        ).unwrap();
        assert_eq!(campaign.run(1).unwrap(), Sdf3EigenfrequencyCampaignStop3::RefinementLimit);
        assert_eq!(campaign.grids().len(), 2);
        assert_eq!(campaign.probes().len(), 2);
        assert_eq!(campaign.candidate_attempts(), 1);
        assert!(campaign.pending().is_none());
        assert!(campaign.rejected().is_none());
        assert!(campaign.current_probe().is_some());
        assert!(campaign.current().study.operator().n() > initial_n);
        assert_physical(campaign.current());
        assert_eq!(&*builds.borrow(), &[
            Sdf3EigenfrequencyBuild3::Probe,
            Sdf3EigenfrequencyBuild3::Refit,
            Sdf3EigenfrequencyBuild3::Probe,
        ]);
    }
}

#[test]
fn cancellation_before_promotion_retains_the_rechecked_candidate_for_resume() {
    let armed = Cell::new(false);
    let cancelled = Cell::new(false);
    let mut callback = |progress: SolveProgress| {
        if armed.get() && !cancelled.get() && progress.stage == "sdf3-modal-campaign-promote" {
            cancelled.set(true);
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (tree, source) = seed(false, &mut control);
    let initial_n = source.study.operator().n();
    let initial_rho = point(&source).rho.clone();
    let builds = Rc::new(RefCell::new(Vec::new()));
    let recorded = Rc::clone(&builds);
    let mut campaign = Sdf3EigenfrequencyCampaign3::new(
        source, tree, forced_refit(), &mut control,
        move |tree, kind, checkpoint| {
            recorded.borrow_mut().push(kind);
            build(tree, 3.0, checkpoint)
        },
    ).unwrap();
    armed.set(true);
    assert!(matches!(
        campaign.run(1),
        Err(Sdf3EigenfrequencyAdaptiveError3::Refinement(
            GoalRefinementError3::Evaluation(EvaluationStop::Cancelled),
        )),
    ));
    assert!(cancelled.get());
    assert_eq!(campaign.grids().len(), 1);
    assert_eq!(campaign.current().study.operator().n(), initial_n);
    assert_eq!(point(campaign.current()).rho, initial_rho);
    assert!(campaign.pending().is_some());
    assert_eq!(campaign.candidate_attempts(), 1);

    armed.set(false);
    assert_eq!(campaign.run(1).unwrap(), Sdf3EigenfrequencyCampaignStop3::RefinementLimit);
    assert_eq!(campaign.grids().len(), 2);
    assert_eq!(campaign.probes().len(), 2);
    assert_eq!(campaign.candidate_attempts(), 1);
    assert!(campaign.pending().is_none());
    assert_physical(campaign.current());
    assert_eq!(&*builds.borrow(), &[
        Sdf3EigenfrequencyBuild3::Probe,
        Sdf3EigenfrequencyBuild3::Refit,
        Sdf3EigenfrequencyBuild3::Probe,
    ]);
}


#[test]
fn detached_pending_candidate_resumes_without_rebuild_or_budget_reset() {
    let armed = Cell::new(false);
    let cancelled = Cell::new(false);
    let mut first_callback = |progress: SolveProgress| {
        if armed.get() && !cancelled.get() && progress.stage == "sdf3-modal-campaign-promote" {
            cancelled.set(true);
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut first_control = SolveControl::new(SolveBudget::default(), &mut first_callback);
    let (tree, source) = seed(false, &mut first_control);
    let builds = Rc::new(RefCell::new(Vec::new()));
    let first_builds = Rc::clone(&builds);
    let mut campaign = Sdf3EigenfrequencyCampaign3::new(
        source, tree, forced_refit(), &mut first_control,
        move |tree, kind, checkpoint| {
            first_builds.borrow_mut().push(kind);
            build(tree, 3.0, checkpoint)
        },
    ).unwrap();
    armed.set(true);
    assert!(campaign.run(1).is_err());
    assert!(campaign.pending().is_some());
    let spent = campaign.work();
    let state = campaign.detach();
    let receipt: SolveControlCheckpoint = state.control_checkpoint();
    assert_eq!(receipt.work(), spent);
    assert!(state.pending().is_some());
    assert_eq!(state.candidate_attempts(), 1);
    drop(first_control);

    armed.set(false);
    let mut second_callback = |_| poll();
    let mut second_control = SolveControl::resume(receipt, &mut second_callback);
    let second_builds = Rc::clone(&builds);
    let mut resumed = Sdf3EigenfrequencyCampaign3::resume(
        state, &mut second_control,
        move |tree, kind, checkpoint| {
            second_builds.borrow_mut().push(kind);
            build(tree, 3.0, checkpoint)
        },
    ).unwrap();
    assert_eq!(resumed.work(), spent);
    assert_eq!(resumed.run(1).unwrap(), Sdf3EigenfrequencyCampaignStop3::RefinementLimit);
    assert_eq!(resumed.candidate_attempts(), 1);
    assert_eq!(resumed.grids().len(), 2);
    assert_eq!(resumed.probes().len(), 2);
    assert!(resumed.pending().is_none());
    assert_physical(resumed.current());
    assert_eq!(&*builds.borrow(), &[
        Sdf3EigenfrequencyBuild3::Probe,
        Sdf3EigenfrequencyBuild3::Refit,
        Sdf3EigenfrequencyBuild3::Probe,
    ]);
}
