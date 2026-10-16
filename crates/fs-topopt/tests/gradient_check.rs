//! Public-consumer numerical gradient gates, rollback, and shared-work tests.
use std::ops::ControlFlow;

use fs_topopt::continuation::{
    ContinuationTermination, controlled_gradient_checked_multi_load_continuation,
};
use fs_topopt::gradient_check::{
    GradientCheckOptions, GradientDirection, controlled_multi_load_gradient_check,
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
        if p[0] > 1.0 - 1e-12 { force[3 * v + 2] = -1.0; }
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

#[test]
fn g3_real_independent_load_gradients_pass_at_multiple_continuation_models() {
    let (mut pipeline, mut elasticity, mut rho, force, volumes) = fixture();
    // A nonuniform design prevents a constant-preserving filter from making
    // the whole experiment a uniform stiffness-scaling special case.
    for (index, r) in rho.iter_mut().enumerate() {
        *r = 0.35 + 0.03 * (index % 9) as f64;
    }
    let opposite: Vec<f64> = force.iter().map(|f| -f).collect();
    let loads = [LoadCase { force: &force, weight: 0.25 },
        LoadCase { force: &opposite, weight: 0.75 }];
    let original = elasticity.moduli.clone();
    for (penal, beta) in [(1.0, 1.0), (3.0, 2.0), (3.0, 8.0)] {
        pipeline.params.penal = penal;
        pipeline.params.beta = beta;
        let mut callback = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let audit = controlled_multi_load_gradient_check(&pipeline, &mut elasticity,
            &loads, &rho, &volumes, GradientCheckOptions::default(), &mut control).unwrap();
        assert!(audit.passed(), "p={penal} beta={beta}: {audit:?}");
        assert_eq!(audit.probes.len(), 2);
        assert!(audit.baseline_compliance > 0.0);
        assert!(audit.work.linear_iterations > 0);
        assert_eq!(elasticity.moduli, original);
        for probe in &audit.probes {
            assert_eq!(probe.active_densities, rho.len());
            assert!(probe.volume_analytic.abs() > 0.0);
        }
    }
}

#[test]
fn g0_upper_bound_uses_only_the_inward_probe_without_counting_an_empty_mask() {
    let (pipeline, mut elasticity, mut rho, force, volumes) = fixture();
    rho.fill(1.0);
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let audit = controlled_multi_load_gradient_check(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        GradientCheckOptions::default(), &mut control).unwrap();
    assert_eq!(audit.probes.len(), 1);
    assert_eq!(audit.probes[0].direction, GradientDirection::Decreasing);
    assert_eq!(audit.probes[0].active_densities, rho.len());
    assert!(audit.passed(), "{audit:?}");
}

#[test]
fn g4_cancel_inside_a_real_probe_restores_the_original_operator() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let mut evaluations = 0;
    let mut callback = |progress: SolveProgress| {
        if progress.stage == "evaluation" { evaluations += 1; }
        if evaluations == 2 && progress.stage == "elasticity" && progress.solve_iterations > 0 {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let outcome = controlled_multi_load_gradient_check(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        GradientCheckOptions::default(), &mut control);
    assert!(matches!(outcome, Err(EvaluationStop::Cancelled)));
    assert_eq!(elasticity.moduli, original);
    assert!(control.work().linear_iterations > 0);
}

#[test]
fn g4_unfunded_gradient_audit_cannot_return_partial_evidence() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {
        total_iterations: 1, ..SolveBudget::default()
    }, &mut callback);
    let outcome = controlled_multi_load_gradient_check(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        GradientCheckOptions::default(), &mut control);
    assert!(matches!(outcome, Err(EvaluationStop::TotalBudget { .. })));
    assert_eq!(control.work().linear_iterations, 1);
    assert_eq!(elasticity.moduli, original);
}

#[test]
fn g0_unrepresentable_probe_step_is_refused_not_reported_as_zero_error() {
    let (pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let outcome = controlled_multi_load_gradient_check(&pipeline, &mut elasticity,
        &[LoadCase { force: &force, weight: 1.0 }], &rho, &volumes,
        GradientCheckOptions { step: f64::MIN_POSITIVE, ..GradientCheckOptions::default() },
        &mut control);
    assert!(matches!(outcome, Err(EvaluationStop::Breakdown { stage: "gradient-check-stencil" })));
    assert_eq!(elasticity.moduli, original);
}

#[test]
fn g3_checked_continuation_records_an_audit_of_every_new_stage_baseline() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let stages = [(1.0, 1.0), (3.0, 4.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..SimpParams::default() });
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = controlled_gradient_checked_multi_load_continuation(
        &mut pipeline, &mut elasticity, &[LoadCase { force: &force, weight: 1.0 }],
        &rho, &volumes, &stages,
        MultiLoadOcOptions { max_iterations: 1, ..MultiLoadOcOptions::default() },
        GradientCheckOptions::default(), &mut control,
    );
    assert_eq!(report.termination, ContinuationTermination::ScheduleComplete, "{report:?}");
    assert_eq!(report.stages.len(), 2);
    for stage in &report.stages {
        let audit = stage.gradient_check.as_ref().unwrap();
        assert!(audit.passed());
        assert_eq!(audit.params.beta.to_bits(), stage.params.beta.to_bits());
        assert_eq!(audit.baseline_compliance.to_bits(), stage.history[0].compliance.to_bits());
        assert_eq!(audit.baseline_volume_fraction.to_bits(), stage.history[0].volume_fraction.to_bits());
    }
    assert!(report.rejected_gradient_check.is_none());
    let last = report.last.unwrap();
    assert_eq!(elasticity.moduli, pipeline.forward(&last.rho).2);
}

#[test]
fn g0_failed_numerical_gate_refuses_stage_before_publishing_new_design() {
    let (mut pipeline, mut elasticity, rho, force, volumes) = fixture();
    let original = elasticity.moduli.clone();
    let original_params = pipeline.params;
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = controlled_gradient_checked_multi_load_continuation(
        &mut pipeline, &mut elasticity, &[LoadCase { force: &force, weight: 1.0 }],
        &rho, &volumes, &[SimpParams { penal: 1.0, beta: 0.0, ..SimpParams::default() }],
        MultiLoadOcOptions::default(),
        // Reciprocal stiffness has nonzero cubic derivatives: this deliberately
        // coarse stencil cannot meet a near-machine-precision comparison gate.
        GradientCheckOptions { step: 0.01, relative_tolerance: 1e-12 }, &mut control,
    );
    assert_eq!(report.termination, ContinuationTermination::GradientCheckFailed);
    assert_eq!(report.stopped_stage, Some(0));
    assert!(!report.rejected_gradient_check.unwrap().passed());
    assert!(report.last.is_none() && report.stages.is_empty());
    assert!(report.evaluation_stop.is_none());
    assert_eq!(elasticity.moduli, original);
    assert_eq!(pipeline.params.beta.to_bits(), original_params.beta.to_bits());
}
