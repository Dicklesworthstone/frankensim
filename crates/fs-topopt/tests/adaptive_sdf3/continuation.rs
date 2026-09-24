use super::*;
use fs_cutfem::elastic3::{ElasticityError3, surface::ReferenceLoad3};
use fs_topopt::sdf3::adaptive_continuation::{
    AdaptiveContinuationError3, AdaptiveContinuationOptions3, controlled_adaptive_sdf3_continuation,
};
use fs_topopt::sdf3_goal::{GoalReferenceLoad3, GoalRefinementError3};
use fs_topopt::{ContinuationTermination, GradientCheckOptions};

fn build(
    t: &Octree3,
    checkpoint: &mut dyn FnMut() -> ControlFlow<()>,
) -> Result<AdaptiveElasticity3, GoalRefinementError3> {
    let mut poll = |_| checkpoint();
    let mut qc = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 1,
            ..Default::default()
        },
        &mut poll,
    )
    .map_err(ElasticityError3::from)?;
    Ok(AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        t,
        &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3::default(),
        &mut qc,
    )?)
}

fn policy() -> AdaptiveContinuationOptions3 {
    AdaptiveContinuationOptions3 {
        optimization: MultiLoadOcOptions {
            max_iterations: 1,
            change_tolerance: 0.0,
            ..Default::default()
        },
        gradient: GradientCheckOptions::default(),
        max_marks: 1,
        ..Default::default()
    }
}

fn schedule() -> [SimpParams; 3] {
    [(1.0, 1.0), (2.0, 2.0), (3.0, 4.0)].map(|(penal, beta)| SimpParams {
        penal,
        beta,
        ..Default::default()
    })
}

#[test]
fn g1_dwr_continuation_performs_two_refinements_and_checks_every_new_design_map() {
    let mut tree = Octree3::uniform(1, 4, 1000).unwrap();
    let (mut study, force) = fixture(&tree, SimpParams::default());
    let original_cells = study.cells();
    let raw = vec![0.5; study.cells()];
    let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [GoalReferenceLoad3 {
        load: ReferenceLoad3::body(&body),
        weight: 1.0,
    }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_adaptive_sdf3_continuation(
        &mut study,
        &mut tree,
        &loads,
        &raw,
        &schedule(),
        policy(),
        &mut control,
        build,
    );
    assert_eq!(
        report.continuation.termination,
        ContinuationTermination::ScheduleComplete,
        "{report:?}"
    );
    assert_eq!(report.continuation.stages.len(), 3);
    assert_eq!(report.refinements.len(), 2);
    assert!(report.refinement_error.is_none());
    for refinement in &report.refinements {
        assert!(refinement.installed);
        assert!(refinement.target_active_cells > refinement.source_active_cells);
        assert!(refinement.target_background_cells > refinement.source_background_cells);
        assert_eq!(refinement.marking.marked.len(), 1);
        let marked = refinement.marking.marked[0];
        let mass = refinement.estimate.marking_mass[&marked];
        assert!(mass > 0.0);
        assert!(
            refinement
                .estimate
                .marking_mass
                .values()
                .all(|m| *m <= mass)
        );
        assert!(
            (refinement.estimate.fine_value
                - refinement.estimate.coarse_value
                - refinement.estimate.correction)
                .abs()
                < 1e-8
        );
    }
    for (index, stage) in report.continuation.stages.iter().enumerate() {
        assert_eq!(stage.stage, index);
        assert!(stage.gradient_check.as_ref().unwrap().passed());
        assert_eq!(stage.history.len(), 2);
        assert!(stage.history[1].compliance < stage.history[0].compliance);
        assert!(
            stage
                .history
                .iter()
                .all(|row| row.volume_fraction <= 0.5 + 1e-8)
        );
    }
    assert!(study.cells() > original_cells);
    assert_eq!(report.continuation.work, control.work());
    let accepted = report.continuation.last.as_ref().unwrap();
    let assembled = study
        .operator()
        .reference_load(loads[0].load, || ControlFlow::Continue(()))
        .unwrap();
    assert_ne!(assembled.len(), force.len());
    let fresh = study
        .evaluate(
            &accepted.rho,
            &[LoadCase {
                force: &assembled,
                weight: 1.0,
            }],
            &mut control,
        )
        .unwrap();
    assert_eq!(accepted.displacements, fresh.objective.displacements);
    assert_eq!(accepted.projected_rho, fresh.projected_rho);
}

#[test]
fn g4_cancelled_second_proposal_preserves_first_refined_stage_and_spent_work() {
    let execute = |cancel: bool| {
        let mut tree = Octree3::uniform(1, 4, 1000).unwrap();
        let (mut study, _) = fixture(&tree, SimpParams::default());
        let raw = vec![0.5; study.cells()];
        let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
        let loads = [GoalReferenceLoad3 {
            load: ReferenceLoad3::body(&body),
            weight: 1.0,
        }];
        let mut stage = 0;
        let mut poll = |p: SolveProgress| {
            if p.stage == "sdf3-adaptive-stage" {
                stage += 1;
            }
            if cancel && stage == 3 && p.stage == "gradient-check-publish" {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let stages = schedule();
        let stages = if cancel { &stages[..] } else { &stages[..2] };
        let report = controlled_adaptive_sdf3_continuation(
            &mut study,
            &mut tree,
            &loads,
            &raw,
            stages,
            policy(),
            &mut control,
            build,
        );
        (tree, study, report)
    };
    let (tree, study, accepted) = execute(false);
    let (stopped_tree, stopped_study, stopped) = execute(true);
    assert_eq!(
        accepted.continuation.termination,
        ContinuationTermination::ScheduleComplete
    );
    assert_eq!(
        stopped.continuation.evaluation_stop,
        Some(EvaluationStop::Cancelled)
    );
    assert_eq!(stopped.continuation.stopped_stage, Some(2));
    assert_eq!(stopped.continuation.stages.len(), 2);
    assert_eq!(stopped.refinements.len(), 2);
    assert!(stopped.refinements[0].installed);
    assert!(!stopped.refinements[1].installed);
    assert_eq!(tree.leaves(), stopped_tree.leaves());
    assert_eq!(study.operator().scales(), stopped_study.operator().scales());
    assert_eq!(study.params().beta, stopped_study.params().beta);
    let good = accepted.continuation.last.unwrap();
    let retained = stopped.continuation.last.unwrap();
    assert_eq!(good.rho, retained.rho);
    assert_eq!(good.projected_rho, retained.projected_rho);
    assert_eq!(good.displacements, retained.displacements);
    assert!(
        stopped.continuation.work.linear_iterations > accepted.continuation.work.linear_iterations
    );
}

#[test]
fn g4_enrichment_leaf_budget_retains_the_solved_coarse_result() {
    let mut tree = Octree3::uniform(1, 4, 8).unwrap();
    let (mut study, _) = fixture(&tree, SimpParams::default());
    let raw = vec![0.5; study.cells()];
    let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [GoalReferenceLoad3 {
        load: ReferenceLoad3::body(&body),
        weight: 1.0,
    }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_adaptive_sdf3_continuation(
        &mut study,
        &mut tree,
        &loads,
        &raw,
        &schedule(),
        policy(),
        &mut control,
        build,
    );
    assert_eq!(
        report.refinement_error,
        Some(AdaptiveContinuationError3::Background(
            fs_cutfem::octree3::OctreeError3::LeafBudget
        ))
    );
    assert_eq!(report.continuation.stopped_stage, Some(1));
    assert_eq!(report.continuation.stages.len(), 1);
    assert!(report.continuation.last.is_some());
    assert_eq!(tree.leaves().len(), 8);
    assert_eq!(report.continuation.work, control.work());
}

#[test]
fn g4_enrichment_spends_the_shared_linear_budget_without_replacing_the_coarse_design() {
    let execute = |total_iterations, stages: usize| {
        let mut tree = Octree3::uniform(1, 4, 1000).unwrap();
        let (mut study, _) = fixture(&tree, SimpParams::default());
        let raw = vec![0.5; study.cells()];
        let body = |_: [f64; 3]| [0.0, 0.0, -1.0];
        let loads = [GoalReferenceLoad3 {
            load: ReferenceLoad3::body(&body),
            weight: 1.0,
        }];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(
            SolveBudget {
                total_iterations,
                ..Default::default()
            },
            &mut poll,
        );
        let report = controlled_adaptive_sdf3_continuation(
            &mut study,
            &mut tree,
            &loads,
            &raw,
            &schedule()[..stages],
            policy(),
            &mut control,
            build,
        );
        (study, report)
    };
    let (coarse, reference) = execute(usize::MAX, 1);
    assert_eq!(
        reference.continuation.termination,
        ContinuationTermination::ScheduleComplete
    );
    let cap = reference.continuation.work.linear_iterations + 1;
    let (stopped_study, stopped) = execute(cap, 3);
    assert_eq!(stopped.continuation.stopped_stage, Some(1));
    assert!(matches!(
        stopped.continuation.evaluation_stop,
        Some(EvaluationStop::TotalBudget {
            stage: "sdf3-goal-elasticity"
        })
    ));
    assert_eq!(stopped.continuation.work.linear_iterations, cap);
    assert_eq!(
        coarse.operator().leaves(),
        stopped_study.operator().leaves()
    );
    assert_eq!(
        coarse.operator().scales(),
        stopped_study.operator().scales()
    );
    let old = reference.continuation.last.unwrap();
    let retained = stopped.continuation.last.unwrap();
    assert_eq!(old.rho, retained.rho);
    assert_eq!(old.displacements, retained.displacements);
}
