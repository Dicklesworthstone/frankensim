use super::*;
use fs_topopt::sdf3::continuation::{
    controlled_gradient_checked_sdf3_continuation, controlled_sdf3_continuation,
    controlled_sdf3_gradient_check,
};
use fs_topopt::{ContinuationTermination, GradientCheckOptions};

fn schedule() -> [SimpParams; 3] {
    [(1.0, 1.0), (2.0, 2.0), (3.0, 4.0)].map(|(penal, beta)| SimpParams {
        penal,
        beta,
        ..Default::default()
    })
}

#[test]
fn cut_cell_continuation_gates_each_model_and_retains_real_fields() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let keys = study.operator().cell_keys();
    let nodes = study.operator().nodes().to_vec();
    let rho = vec![0.5; study.cells()];
    let loads = [
        LoadCase {
            force: &y,
            weight: 0.3,
        },
        LoadCase {
            force: &z,
            weight: 0.7,
        },
    ];
    let run = |study: &mut CutDensityStudy3| {
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        controlled_gradient_checked_sdf3_continuation(
            study,
            &loads,
            &rho,
            &schedule(),
            MultiLoadOcOptions {
                max_iterations: 2,
                ..options()
            },
            GradientCheckOptions::default(),
            &mut control,
        )
    };
    let report = run(&mut study);
    assert_eq!(
        report.termination,
        ContinuationTermination::ScheduleComplete,
        "{report:?}"
    );
    assert_eq!(report.stages.len(), 3);
    for stage in &report.stages {
        let check = stage.gradient_check.as_ref().unwrap();
        assert!(check.passed());
        assert_eq!(check.probes.len(), 2);
        assert!(stage.history.len() > 1);
        assert!(stage.history.last().unwrap().compliance < stage.history[0].compliance);
        for row in &stage.history {
            assert!(row.volume_fraction <= 0.5 + options().volume_tolerance);
            assert_eq!(row.case_compliances.len(), 2);
            assert!(
                (row.compliance - (0.3 * row.case_compliances[0] + 0.7 * row.case_compliances[1]))
                    .abs()
                    <= 1e-12 * row.compliance
            );
        }
        for pair in stage.history.windows(2) {
            assert!(pair[1].compliance <= pair[0].compliance);
        }
    }
    assert_eq!(study.operator().cell_keys(), keys);
    assert_eq!(study.operator().nodes(), nodes);
    let last = report.last.as_ref().unwrap();
    let accepted_scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let solved = study.evaluate(&last.rho, &loads, &mut control).unwrap();
    assert_eq!(last.displacements, solved.objective.displacements);
    assert_eq!(last.projected_rho, solved.projected_rho);
    assert_eq!(
        last.history.last().unwrap().compliance.to_bits(),
        solved.objective.compliance.to_bits()
    );
    assert_eq!(study.operator().scales(), accepted_scales);
    assert_eq!(study.params().beta, 4.0);
    let repeated = run(&mut study);
    assert_eq!(report.work, repeated.work);
    assert_eq!(last.rho, repeated.last.as_ref().unwrap().rho);
    assert_eq!(
        last.displacements,
        repeated.last.as_ref().unwrap().displacements
    );
}

#[test]
fn sharper_projection_restores_new_model_volume_before_baseline() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.4; study.cells()];
    let stages = [0.0, 8.0].map(|beta| SimpParams {
        beta,
        eta: 0.25,
        ..Default::default()
    });
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let report = controlled_sdf3_continuation(
        &mut study,
        &[LoadCase {
            force: &z,
            weight: 1.0,
        }],
        &rho,
        &stages,
        MultiLoadOcOptions {
            volume_fraction: 0.4,
            max_iterations: 0,
            ..options()
        },
        &mut control,
    );
    assert_eq!(
        report.termination,
        ContinuationTermination::ScheduleComplete
    );
    assert_eq!(report.stages.len(), 2);
    assert!(report.stages[1].incoming_volume_fraction > 0.4);
    assert!(report.stages[1].restoration_scale < 1.0);
    assert!(report.stages[1].history[0].volume_fraction <= 0.4 + options().volume_tolerance);
    assert!(report.last.as_ref().unwrap().rho.iter().all(|r| *r < 0.4));
}

#[test]
fn interruption_inside_new_model_gradient_probe_restores_old_model_and_fields() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.5; study.cells()];
    let mut stage_count = 0;
    let mut poll = |p: SolveProgress| {
        if p.stage == "sdf3-continuation-stage" {
            stage_count += 1;
        }
        if stage_count == 2 && p.stage == "gradient-probe" {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let loads = [LoadCase {
        force: &z,
        weight: 1.0,
    }];
    let report = controlled_gradient_checked_sdf3_continuation(
        &mut study,
        &loads,
        &rho,
        &schedule(),
        MultiLoadOcOptions {
            max_iterations: 1,
            ..options()
        },
        GradientCheckOptions::default(),
        &mut control,
    );
    assert_eq!(
        report.termination,
        ContinuationTermination::EvaluationStopped
    );
    assert_eq!(report.evaluation_stop, Some(EvaluationStop::Cancelled));
    assert_eq!(report.stopped_stage, Some(1));
    assert_eq!(report.stages.len(), 1);
    assert_eq!(study.params().penal, 1.0);
    let last = report.last.unwrap();
    assert!(report.work.linear_iterations > last.work.linear_iterations);
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let evaluated = study.evaluate(&last.rho, &loads, &mut control).unwrap();
    assert_eq!(last.displacements, evaluated.objective.displacements);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn failed_gradient_gate_and_exhausted_budget_never_install_an_unchecked_model() {
    let (mut study, _, z) = fixture(SimpParams::default());
    let rho = vec![0.4; study.cells()];
    let scales = study.operator().scales().to_vec();
    let params = study.params();
    let loads = [LoadCase {
        force: &z,
        weight: 1.0,
    }];
    for budget in [0, usize::MAX] {
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(
            SolveBudget {
                total_iterations: budget,
                ..Default::default()
            },
            &mut poll,
        );
        let report = controlled_gradient_checked_sdf3_continuation(
            &mut study,
            &loads,
            &rho,
            &schedule(),
            options(),
            GradientCheckOptions {
                step: 0.1,
                relative_tolerance: 1e-8,
            },
            &mut control,
        );
        assert!(report.last.is_none());
        assert_eq!(report.stopped_stage, Some(0));
        if budget == 0 {
            assert_eq!(
                report.termination,
                ContinuationTermination::EvaluationStopped
            );
            assert!(report.rejected_gradient_check.is_none());
        } else {
            assert_eq!(
                report.termination,
                ContinuationTermination::GradientCheckFailed
            );
            assert!(!report.rejected_gradient_check.unwrap().passed());
        }
        assert_eq!(study.operator().scales(), scales);
        assert_eq!(study.params().penal, params.penal);
        assert_eq!(study.params().beta, params.beta);
    }
}

#[test]
fn standalone_gradient_audit_preserves_incoming_scales() {
    let (mut study, y, z) = fixture(SimpParams::default());
    let rho: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.03 * i as f64).collect();
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let check = controlled_sdf3_gradient_check(
        &mut study,
        &[
            LoadCase {
                force: &y,
                weight: 0.3,
            },
            LoadCase {
                force: &z,
                weight: 0.7,
            },
        ],
        &rho,
        GradientCheckOptions::default(),
        &mut control,
    )
    .unwrap();
    assert!(check.passed(), "{check:?}");
    assert_eq!(study.operator().scales(), scales);
}
