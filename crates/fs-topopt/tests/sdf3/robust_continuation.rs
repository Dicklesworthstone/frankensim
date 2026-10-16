use super::*;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlStop};
use fs_topopt::sdf3::robust::{RobustDesignOptions3, RobustDesignStudy3};
use fs_topopt::sdf3::robust::continuation::{RobustContinuationOptions3,
    RobustContinuationStop3, controlled_sdf3_robust_continuation,
    controlled_sdf3_robust_feasible_start};

fn policy() -> RobustContinuationOptions3 {
    RobustContinuationOptions3 { steps_per_stage: 0,
        design: RobustDesignOptions3 { volume_limit: 0.4, compliance_scale: 1000.0,
            ..Default::default() }, ..Default::default() }
}
fn models() -> [SimpParams; 2] {
    [0.0, 8.0].map(|beta| SimpParams { beta, penal: 2.0, eta: 0.35,
        ..Default::default() })
}
fn same_model(a: SimpParams, b: SimpParams) {
    assert_eq!([a.penal, a.beta, a.eta, a.e_min], [b.penal, b.beta, b.eta, b.e_min]);
}

#[test]
fn sharper_stages_restore_dilated_not_nominal_volume_and_recheck_real_fields() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let keys = study.operator().cell_keys();
    let nodes = study.operator().nodes().to_vec();
    let scales = study.operator().scales().to_vec();
    let original = study.params();
    let rho = vec![0.4; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study, &loads, &rho, &models(), policy(), &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::ScheduleComplete, "{report:?}");
    assert_eq!(report.stages.len(), 2);
    assert!(report.stages[1].incoming_dilated_volume_fraction > 0.8);
    assert!(report.stages[1].restoration_scale < 1.0);
    for stage in &report.stages {
        assert!(stage.gradient_check.numerical.passed());
        same_model(stage.params, stage.gradient_check.numerical.params);
        assert_eq!(stage.history.len(), 1);
        assert_eq!(stage.selected_iteration, 0);
        assert!(stage.history[0].dilated_volume_fraction <= 0.4 + 1e-12);
        assert!(stage.history[0].feasible);
        assert_eq!(stage.outcome.as_ref().unwrap().work.evaluations, 1);
    }
    let last = report.last.as_ref().unwrap();
    assert!(last.evaluation.nominal_volume_fraction < last.evaluation.dilated_volume_fraction - 0.1);
    let (mut reference, _, _) = fixture(report.params);
    let checked = reference.evaluate_robust(&last.rho, &loads, policy().design.robustness, &mut control).unwrap();
    assert_eq!(last.evaluation.eroded.objective.displacements, checked.eroded.objective.displacements);
    assert_eq!(last.evaluation.dilated_rho, checked.dilated_rho);
    assert_eq!(last.evaluation.eroded.objective.compliance.to_bits(), checked.eroded.objective.compliance.to_bits());
    assert_eq!(study.operator().cell_keys(), keys);
    assert_eq!(study.operator().nodes(), nodes);
    assert_eq!(study.operator().scales(), scales);
    same_model(study.params(), original);
}

#[test]
fn single_stage_runs_the_existing_optimizer_and_selects_its_incumbent() {
    let params = models()[0];
    let (mut study, _, z) = fixture(params);
    let rho = vec![0.4; study.cells()];
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let options = RobustContinuationOptions3 { steps_per_stage: 8, ..policy() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study, &loads, &rho, &[params], options, &mut control).unwrap();
    assert_eq!(report.stages.len(), 1, "{report:?}");
    let mut direct = RobustDesignStudy3::new(&mut study, &loads,
        &report.stages[0].starting_rho, options.design, &mut control).unwrap();
    let direct_outcome = direct.run(options.steps_per_stage).unwrap();
    assert!(direct.history().len() > 1);
    assert_eq!(report.stages[0].history, direct.history());
    assert_eq!(report.last.as_ref().unwrap().rho, direct.best_feasible().unwrap().rho);
    assert_eq!(report.stages[0].outcome.as_ref().unwrap().stop, direct_outcome.stop);
    assert!(report.stages[0].history[report.stages[0].selected_iteration].feasible);
}

#[test]
fn new_model_probe_cancellation_keeps_previous_feasible_model_and_nominal_state() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let original = study.params();
    let scales = study.operator().scales().to_vec();
    let rho = vec![0.4; study.cells()];
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let mut stages = 0;
    let mut poll = |p: SolveProgress| {
        if p.stage == "sdf3-robust-continuation-stage" { stages += 1; }
        if stages == 2 && p.stage == "gradient-probe" { ControlFlow::Break(()) }
        else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study, &loads, &rho, &models(), policy(), &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::EvaluationStopped);
    assert_eq!(report.stopped_stage, Some(1));
    assert_eq!(report.stages.len(), 1);
    same_model(report.params, models()[0]);
    assert_eq!(report.last.as_ref().unwrap().rho, report.stages[0].starting_rho);
    assert!(report.work.linear_solves > report.last.as_ref().unwrap().evaluation.work.linear_solves);
    assert_eq!(study.operator().scales(), scales);
    same_model(study.params(), original);
}

#[test]
fn interruption_after_an_accepted_step_retains_an_aligned_stage_prefix() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.4; study.cells()];
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-robust-design-accepted" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study,
        &[LoadCase { force: &z, weight: 1.0 }], &rho, &models(),
        RobustContinuationOptions3 { steps_per_stage: 5, ..policy() }, &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::EvaluationStopped, "{report:?}");
    assert_eq!(report.stages.len(), 1);
    assert_eq!(report.stages[0].history.len(), 2);
    assert!(report.stages[0].outcome.is_err());
    assert!(report.error.is_some());
    let selected = &report.stages[0].history[report.stages[0].selected_iteration];
    assert!(selected.feasible);
    assert_eq!(report.last.as_ref().unwrap().evaluation.eroded.objective.compliance, selected.eroded_compliance);
}

#[test]
fn cumulative_physics_budget_is_not_renewed_for_a_new_material_model() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.4; study.cells()];
    let loads = [LoadCase { force: &z, weight: 1.0 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let first = controlled_sdf3_robust_continuation(&mut study, &loads, &rho, &models()[..1], policy(), &mut control).unwrap();
    let budget = first.work.linear_iterations + 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study, &loads, &rho, &models(), policy(), &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::EvaluationStopped);
    assert_eq!(report.stages.len(), 1);
    assert_eq!(report.work.linear_iterations, budget);
    assert_eq!(report.last.as_ref().unwrap().rho, first.last.as_ref().unwrap().rho);
    assert!(matches!(report.error, Some(ProjectedAlError::Evaluation(EvaluationStop::TotalBudget { .. }))));
}

#[test]
fn floor_and_physical_solid_costs_refuse_without_publishing_a_fake_design() {
    let (study, _, z) = fixture(SimpParams::default());
    let regions = vec![fs_topopt::sdf3::PhysicalRegion3::Solid; study.cells()];
    let mut study = study.with_physical_regions(regions).unwrap();
    let scales = study.operator().scales().to_vec();
    let rho = vec![0.7; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study,
        &[LoadCase { force: &z, weight: 1.0 }], &rho, &models(), policy(), &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::RestorationFailed);
    assert!(report.last.is_none() && report.stages.is_empty());
    assert_eq!(study.operator().scales(), scales);
    assert!(report.work.linear_solves > 0);
}

#[test]
fn feasible_restoration_preserves_authored_regions_and_is_nonmutating_on_cancel() {
    use fs_topopt::sdf3::PhysicalRegion3;
    let (study, _, z) = fixture(models()[1]);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid;
    regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let scales = study.operator().scales().to_vec();
    let original = study.params();
    let rho = vec![0.8; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let restored = controlled_sdf3_robust_feasible_start(&mut study, &rho, policy().design, &mut control).unwrap();
    let evaluated = study.evaluate_robust(&restored, &[LoadCase { force: &z, weight: 1.0 }],
        policy().design.robustness, &mut control).unwrap();
    assert!(evaluated.dilated_volume_fraction <= 0.4 + 1e-12);
    for field in [&evaluated.eroded.projected_rho, &evaluated.nominal_rho, &evaluated.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
    let mut poll = |p: SolveProgress| if p.stage == "sdf3-volume-restoration" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(controlled_sdf3_robust_feasible_start(&mut study, &rho, policy().design, &mut control).is_err());
    same_model(study.params(), original);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn optimizer_limits_are_not_reported_as_completed_continuation() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.4; study.cells()];
    let mut options = policy();
    options.steps_per_stage = 5;
    options.design.optimizer.max_evaluations = 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_robust_continuation(&mut study,
        &[LoadCase { force: &z, weight: 1.0 }], &rho, &models(), options, &mut control).unwrap();
    assert_eq!(report.stop, RobustContinuationStop3::OptimizerStopped(ProjectedAlStop::EvaluationLimit));
    assert_eq!(report.stopped_stage, Some(0));
    assert_eq!(report.stages.len(), 1);
    assert_eq!(report.last.as_ref().unwrap().rho, report.stages[0].starting_rho);
}
