//! Full nuisance covariance and signed-Jacobian integration regressions.
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, JacobianRecord, SharedNoiseModel};
use fs_assimilate::{AssimError, Belief, Observation, assimilate_all, point_sensor};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_budget<R>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 6416, kernel_id: 16, tile: 0, iteration: 0 },
            budget, ExecMode::Deterministic,
        ).with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2.0e-10 * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
}

fn record(model: &SharedNoiseModel, row: &str, jacobian: Vec<f64>) -> JacobianRecord {
    JacobianRecord::new(
        point_sensor(0, 1, 1.0, 0.25, "scan-instrument").unwrap(),
        RowId::try_new(row).unwrap(), model, jacobian,
    ).unwrap()
}

#[test]
fn full_covariance_propagates_parameter_correlations_and_signed_sensitivities() {
    with_budget(Budget::INFINITE, |cx| {
        let c = vec![vec![4.0, 1.0], vec![1.0, 2.0]];
        let model = SharedNoiseModel::new("registration-fit", vec!["translation".into(), "rotation".into()], c.clone(), cx).unwrap();
        let j = [vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, -2.0]];
        let rows = (0..3).map(|i| record(&model, &format!("row-{i}"), j[i].clone())).collect();
        let batch = JacobianBatch::new(model, rows, cx).unwrap();
        for i in 0..3 {
            for q in 0..3 {
                let mut expected = if i == q { 0.25 } else { 0.0 };
                for a in 0..2 { for b in 0..2 { expected += j[i][a] * c[a][b] * j[q][b]; } }
                close(batch.covariance()[i][q], expected);
                assert_eq!(batch.covariance()[i][q].to_bits(), batch.covariance()[q][i].to_bits());
            }
        }
        close(batch.covariance()[1][2], -3.0);
        let prior = Belief::scalar(0.0, 10.0).unwrap();
        let first = batch.assimilate(&prior, cx).unwrap();
        let replay = batch.assimilate(&prior, cx).unwrap();
        assert_eq!(first, replay);
        assert!(first.robust().posterior().is_some());
        assert!(first.robust().factor_diagnostics().is_some());
        assert_eq!(first.robust().audit().effective_dof(), 3);
    });
}

#[test]
fn opposite_signed_common_errors_cancel_without_discarding_covariance() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("offset-fit", vec!["offset".into()], vec![vec![4.0]], cx).unwrap();
        let rows = vec![record(&model, "left", vec![1.0]), record(&model, "right", vec![-1.0])];
        let batch = JacobianBatch::new(model, rows, cx).unwrap();
        close(batch.covariance()[0][1], -4.0);
        let prior = Belief::scalar(0.0, 10.0).unwrap();
        let result = batch.assimilate(&prior, &cx).unwrap();
        let posterior = result.robust().posterior().unwrap();
        // The average rejects the nuisance; only instrument noise remains.
        let expected_variance = 1.0 / (0.1 + 2.0 / 0.25);
        close(posterior.variance(0).unwrap(), expected_variance);
        close(posterior.mean()[0], expected_variance * 2.0 / 0.25);
        let discarded_correlation = 1.0 / (0.1 + 2.0 / 4.25);
        assert!(posterior.variance(0).unwrap() < discarded_correlation / 10.0);
    });
}

#[test]
fn common_mode_repetitions_respect_the_shared_uncertainty_floor() {
    with_budget(Budget::INFINITE, |cx| {
        let prior = Belief::scalar(0.0, 10.0).unwrap();
        for count in [2, 8, 32] {
            let model = SharedNoiseModel::new("offset-fit", vec!["offset".into()], vec![vec![4.0]], cx).unwrap();
            let rows = (0..count).map(|i| record(&model, &format!("sample-{i}"), vec![1.0])).collect();
            let batch = JacobianBatch::new(model, rows, cx).unwrap();
            let result = batch.assimilate(&prior, cx).unwrap();
            let information = count as f64 / (0.25 + 4.0 * count as f64);
            let variance = result.robust().posterior().unwrap().variance(0).unwrap();
            close(variance, 1.0 / (0.1 + information));
            assert!(variance > 1.0 / (0.1 + 1.0 / 4.0));
        }
    });
}

#[test]
fn singular_nuisance_covariance_preserves_its_exact_null_direction() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("coupled-offsets", vec!["a".into(), "b".into()], vec![vec![0.25, 0.25], vec![0.25, 0.25]], cx).unwrap();
        let rows = vec![record(&model, "difference", vec![1.0, -1.0]), record(&model, "sum", vec![1.0, 1.0])];
        let batch = JacobianBatch::new(model, rows, cx).unwrap();
        assert_eq!(batch.covariance(), &[vec![0.25, 0.0], vec![0.0, 1.25]]);
        let result = batch.assimilate(&Belief::scalar(0.0, 10.0).unwrap(), cx).unwrap();
        close(result.robust().posterior().unwrap().variance(0).unwrap(), 1.0 / (0.1 + 1.0 / 0.25 + 1.0 / 1.25));
    });
}

#[test]
fn zero_shared_covariance_reduces_to_independent_sensor_fusion() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("zero-offset", vec!["offset".into()], vec![vec![-0.0]], cx).unwrap();
        let rows = vec![record(&model, "row-a", vec![2.0]), record(&model, "row-b", vec![-3.0])];
        let independent: Vec<Observation> = rows.iter().map(|r| r.observation().clone()).collect();
        let batch = JacobianBatch::new(model, rows, cx).unwrap();
        assert_eq!(batch.covariance(), &[vec![0.25, 0.0], vec![0.0, 0.25]]);
        let prior = Belief::scalar(0.0, 10.0).unwrap();
        let reference = assimilate_all(&prior, &independent, cx).unwrap();
        let result = batch.assimilate(&prior, cx).unwrap();
        assert_eq!(result.robust().posterior(), Some(&reference));
        assert!(result.robust().factor_diagnostics().is_none());
    });
}

#[test]
fn duplicate_rows_mixed_models_and_invalid_jacobians_fail_closed() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("fit", vec!["x".into(), "y".into()], vec![vec![1.0, 0.0], vec![0.0, 1.0]], cx).unwrap();
        let reordered = SharedNoiseModel::new("fit", vec!["y".into(), "x".into()], vec![vec![1.0, 0.0], vec![0.0, 1.0]], cx).unwrap();
        assert_ne!(model.content_id(), reordered.content_id());
        let row = record(&model, "same-row", vec![1.0, 0.0]);
        assert!(matches!(JacobianBatch::new(model.clone(), vec![row.clone(), row.clone()], cx), Err(AssimError::DuplicateDatasetRow { .. })));
        assert!(matches!(JacobianBatch::new(reordered, vec![row], cx), Err(AssimError::SharedSourceDeclaration { .. })));
        for jacobian in [vec![1.0], vec![1.0, 0.0, 0.0], vec![f64::NAN, 0.0], vec![0.0, f64::INFINITY]] {
            assert!(JacobianRecord::new(point_sensor(0, 1, 0.0, 0.25, "instrument").unwrap(), RowId::try_new("row").unwrap(), &model, jacobian).is_err());
        }
        assert!(SharedNoiseModel::new("fit", vec!["x".into(), "x".into()], vec![vec![1.0, 0.0], vec![0.0, 1.0]], cx).is_err());
        assert!(SharedNoiseModel::new("fit", vec!["x".into(), "y".into()], vec![vec![1.0, 2.0], vec![2.0, 1.0]], cx).is_err());
    });
}

#[test]
fn model_and_row_jacobians_are_retained_even_when_the_covariance_is_identical() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("fit", vec!["offset".into()], vec![vec![1.0]], cx).unwrap();
        let renamed = SharedNoiseModel::new("other-fit", vec!["offset".into()], vec![vec![1.0]], cx).unwrap();
        assert_ne!(model.content_id(), renamed.content_id());
        let a = JacobianBatch::new(model.clone(), vec![record(&model, "a", vec![1.0])], cx).unwrap();
        let b = JacobianBatch::new(model.clone(), vec![record(&model, "a", vec![-1.0])], cx).unwrap();
        assert_eq!(a.covariance(), b.covariance());
        assert_ne!(a.content_id(), b.content_id());
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        let ar = a.assimilate(&prior, cx).unwrap();
        let br = b.assimilate(&prior, cx).unwrap();
        assert_eq!(ar.robust().posterior(), br.robust().posterior());
        assert_ne!(ar.content_id(), br.content_id());
        let minus_zero = SharedNoiseModel::new("zero", vec!["offset".into()], vec![vec![-0.0]], cx).unwrap();
        let plus_zero = SharedNoiseModel::new("zero", vec!["offset".into()], vec![vec![0.0]], cx).unwrap();
        assert_eq!(minus_zero.content_id(), plus_zero.content_id());
    });
}

#[test]
fn bounded_projection_refuses_exhausted_budgets_and_nonfinite_products() {
    let (model, rows) = with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("fit", vec!["offset".into()], vec![vec![1.0]], cx).unwrap();
        let rows = (0..16).map(|i| record(&model, &format!("row-{i}"), vec![1.0])).collect::<Vec<_>>();
        (model, rows)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            assert!(matches!(JacobianBatch::new(model.clone(), rows.clone(), cx), Err(AssimError::BudgetRefused(_)) | Err(AssimError::Cancelled { .. })));
        });
    }
    with_budget(Budget::INFINITE, |cx| {
        let row = record(&model, "overflow", vec![f64::MAX]);
        assert!(matches!(JacobianBatch::new(model, vec![row], cx), Err(AssimError::NonFiniteComputation { .. })));
    });
}
