//! Real-elasticity continuation, feasibility, rollback and replay regressions.
use std::ops::ControlFlow;

use fs_topopt::continuation::{
    ContinuationTermination, controlled_multi_load_continuation, multi_load_continuation,
};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop, MultiLoadOcOptions,
    SimpParams, SolveBudget, SolveControl, SolveProgress,
};

fn fixture() -> (DesignPipeline, DensityElasticity, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (complex, positions) = fs_feec::kuhn_cube(2);
    let elasticity = DensityElasticity::new(&complex, &positions, 1.0, 0.3, &|p| p[0] < 1e-12);
    let mut force = vec![0.0; elasticity.n()];
    for (v, p) in positions.iter().enumerate() {
        if p[0] > 1.0 - 1e-12 {
            force[3 * v + 2] = -1.0;
        }
    }
    let rho = vec![0.5; elasticity.cells()];
    let volumes = fs_feec::element_geometry(&complex, &positions)
        .vol_signed.iter().map(|v| v.abs()).collect();
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
    };
    (pipeline, elasticity, rho, force, volumes)
}

fn options() -> MultiLoadOcOptions {
    MultiLoadOcOptions { max_iterations: 2, change_tolerance: 0.0, ..Default::default() }
}

fn schedule() -> [SimpParams; 3] {
    [(1.0, 1.0), (2.0, 2.0), (3.0, 4.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..SimpParams::default() })
}

#[test]
fn g3_continuation_preserves_independent_loads_and_stagewise_descent() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let opposite: Vec<f64> = force.iter().map(|f| -f).collect();
    let loads = [LoadCase { force: &force, weight: 0.25 },
        LoadCase { force: &opposite, weight: 0.75 }];
    let report = multi_load_continuation(&mut pipeline, &mut elasticity, &loads,
        &rho, &volumes, &schedule(), options(), || ControlFlow::Continue(()));
    assert_eq!(report.termination, ContinuationTermination::ScheduleComplete);
    assert_eq!(report.stages.len(), 3);
    assert!(report.stages[0].history.len() > 1);
    for stage in &report.stages {
        assert!(stage.history[0].compliance > 0.0);
        for row in &stage.history {
            assert!(row.volume_fraction <= options().volume_fraction + options().volume_tolerance);
            assert_eq!(row.case_compliances.len(), 2);
            assert!((row.case_compliances[0] - row.case_compliances[1]).abs()
                <= 1e-8 * row.case_compliances[0]);
        }
        for pair in stage.history.windows(2) {
            assert!(pair[1].compliance <= pair[0].compliance);
        }
    }
    let last = report.last.as_ref().unwrap();
    let (_, projected, moduli) = pipeline.forward(&last.rho);
    assert_eq!(last.projected_rho, projected);
    assert_eq!(elasticity.moduli, moduli);
    let solved = pipeline.multi_load_compliance_and_gradient(&mut elasticity, &last.rho, &loads);
    assert_eq!(last.displacements, solved.displacements);
    assert_eq!(last.history.last().unwrap().compliance.to_bits(), solved.compliance.to_bits());
    assert_eq!(pipeline.params.beta.to_bits(), schedule()[2].beta.to_bits());
}

#[test]
fn g0_sharper_projection_restores_actual_volume_before_new_baseline() {
    let (mut pipeline, mut elasticity, mut rho, force, volumes) = fixture();
    rho.fill(0.4);
    let stages = [0.0, 8.0].map(|beta| SimpParams {
        beta, eta: 0.25, ..SimpParams::default()
    });
    let report = multi_load_continuation(&mut pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes, &stages,
        MultiLoadOcOptions { max_iterations: 0, volume_fraction: 0.4, ..options() },
        || ControlFlow::Continue(()));
    assert_eq!(report.termination, ContinuationTermination::ScheduleComplete);
    assert_eq!(report.stages.len(), 2);
    assert!(report.stages[1].incoming_volume_fraction > 0.4);
    assert!(report.stages[1].restoration_scale < 1.0);
    assert!(report.stages[1].history[0].volume_fraction <= 0.4 + options().volume_tolerance);
    assert!(report.last.as_ref().unwrap().rho.iter().all(|r| *r < 0.4));
    assert_eq!(elasticity.moduli, pipeline.forward(&report.last.unwrap().rho).2);
}

#[test]
fn g4_cancel_between_stages_restores_matching_parameters_and_fields() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let mut starts = 0;
    let mut callback = |progress: SolveProgress| {
        if progress.stage == "continuation-stage" { starts += 1; }
        if starts == 2 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = controlled_multi_load_continuation(&mut pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        &schedule(), options(), &mut control);
    assert_eq!(report.termination, ContinuationTermination::EvaluationStopped);
    assert_eq!(report.evaluation_stop, Some(EvaluationStop::Cancelled));
    assert_eq!(report.stopped_stage, Some(1));
    assert_eq!(report.stages.len(), 1);
    assert_eq!(pipeline.params.penal.to_bits(), schedule()[0].penal.to_bits());
    assert_eq!(elasticity.moduli, pipeline.forward(&report.last.unwrap().rho).2);
}

#[test]
fn g4_cancellation_inside_transition_cannot_publish_a_new_model() {
    let (mut pipeline, mut elasticity, mut rho, force, volumes) = fixture();
    rho.fill(0.4);
    let stages = [0.0, 8.0].map(|beta| SimpParams {
        beta, eta: 0.25, ..SimpParams::default()
    });
    let mut callback = |progress: SolveProgress| {
        if progress.stage == "continuation-volume-search" {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = controlled_multi_load_continuation(&mut pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes, &stages,
        MultiLoadOcOptions { max_iterations: 0, volume_fraction: 0.4, ..options() },
        &mut control);
    assert_eq!(report.evaluation_stop, Some(EvaluationStop::Cancelled));
    assert_eq!(report.stages.len(), 1);
    assert_eq!(pipeline.params.beta.to_bits(), 0.0_f64.to_bits());
    let last = report.last.unwrap();
    assert_eq!(last.rho, rho);
    assert!(report.work.linear_iterations > last.work.linear_iterations);
    assert_eq!(elasticity.moduli, pipeline.forward(&rho).2);
}

#[test]
fn g4_cumulative_linear_budget_is_not_reset_by_a_new_stage() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let loads = [LoadCase { force: &force, weight: 1.0 }];
    let opts = MultiLoadOcOptions { max_iterations: 0, ..options() };
    let baseline = multi_load_continuation(&mut pipeline, &mut elasticity, &loads,
        &rho, &volumes, &schedule()[..1], opts, || ControlFlow::Continue(()));
    let limit = baseline.work.linear_iterations + 1;
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {
        total_iterations: limit, ..SolveBudget::default()
    }, &mut callback);
    let report = controlled_multi_load_continuation(&mut pipeline, &mut elasticity,
        &loads, &rho, &volumes, &schedule(), opts, &mut control);
    assert_eq!(report.termination, ContinuationTermination::EvaluationStopped);
    assert_eq!(report.stopped_stage, Some(1));
    assert!(matches!(report.evaluation_stop, Some(EvaluationStop::TotalBudget { .. })));
    assert_eq!(report.work.linear_iterations, limit);
    assert_eq!(report.last.unwrap().rho, baseline.last.unwrap().rho);
    assert_eq!(pipeline.params.beta.to_bits(), schedule()[0].beta.to_bits());
}

#[test]
fn g0_failed_restoration_does_not_claim_infeasibility_or_change_operator() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let original_params = pipeline.params;
    let report = multi_load_continuation(&mut pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        &[SimpParams { beta: 0.0, ..SimpParams::default() }],
        MultiLoadOcOptions { volume_fraction: 1e-6, ..options() },
        || ControlFlow::Continue(()));
    assert_eq!(report.termination, ContinuationTermination::FeasibilityRestorationFailed);
    assert!(report.last.is_none() && report.stages.is_empty());
    assert_eq!(elasticity.moduli, original);
    assert_eq!(pipeline.params.beta.to_bits(), original_params.beta.to_bits());
}

#[test]
fn g5_schedule_replays_designs_fields_work_and_restoration_bitwise() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let loads = [LoadCase { force: &force, weight: 1.0 }];
    let a = multi_load_continuation(&mut pipeline, &mut elasticity, &loads,
        &rho, &volumes, &schedule(), options(), || ControlFlow::Continue(()));
    let b = multi_load_continuation(&mut pipeline, &mut elasticity, &loads,
        &rho, &volumes, &schedule(), options(), || ControlFlow::Continue(()));
    assert_eq!(a.termination, b.termination);
    assert_eq!(a.work, b.work);
    for (a, b) in a.stages.iter().zip(&b.stages) {
        assert_eq!(a.restoration_scale.to_bits(), b.restoration_scale.to_bits());
        assert_eq!(a.incoming_volume_fraction.to_bits(), b.incoming_volume_fraction.to_bits());
        assert_eq!(a.history.len(), b.history.len());
        for (a, b) in a.history.iter().zip(&b.history) {
            assert_eq!(a.compliance.to_bits(), b.compliance.to_bits());
            assert_eq!(a.volume_fraction.to_bits(), b.volume_fraction.to_bits());
        }
    }
    let (a, b) = (a.last.unwrap(), b.last.unwrap());
    assert_eq!(a.rho, b.rho);
    assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.displacements, b.displacements);
}

#[test]
fn g0_invalid_later_stage_is_rejected_before_any_model_mutation() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let original_params = pipeline.params;
    let mut stages = schedule();
    stages[2].beta = f64::NAN;
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        multi_load_continuation(&mut pipeline, &mut elasticity,
            &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
            &stages, options(), || ControlFlow::Continue(()))
    }));
    assert!(result.is_err());
    assert_eq!(elasticity.moduli, original);
    assert_eq!(pipeline.params.beta.to_bits(), original_params.beta.to_bits());
}
