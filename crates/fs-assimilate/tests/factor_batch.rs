//! Full-batch regression tests use the dense Joseph and direct block formulas
//! as independent references, not a second copy of the UD update.
use fs_assimilate::factor::{FactorBelief, assimilate_scalar, batch::assimilate_sequence};
use fs_assimilate::{
    AssimError, Belief, Observation, ObservationBatch, RobustObservation, assimilate_all,
    assimilate_observation_batch, misfit, point_sensor,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    let stream = StreamKey { seed: 6416, kernel_id: 64, tile: 0, iteration: 0 };
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, stream, budget, ExecMode::Deterministic)
            .with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2.0e-11 * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn retained_sequence_matches_scalar_and_dense_oracles() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::new(vec![0.5, -1.0], vec![vec![2.0, 0.4], vec![0.4, 1.0]], cx).unwrap();
        let observations: Vec<_> = (0..12).map(|i| {
            Observation::new(vec![1.0, 0.1 * i as f64 - 0.5], i as f64 / 4.0, 0.4, format!("sensor-{i}")).unwrap()
        }).collect();
        let factor = FactorBelief::from_belief(&prior, cx).unwrap();
        let saved = factor.clone();
        let result = assimilate_sequence(&factor, &observations, cx).unwrap();
        let mut sequential = factor.clone();
        for (observation, receipt) in observations.iter().zip(result.diagnostics().receipts()) {
            let scalar = assimilate_scalar(&sequential, observation, cx).unwrap();
            assert_eq!(receipt, scalar.receipt());
            sequential = scalar.belief().clone();
        }
        assert_eq!(result.belief(), &sequential);
        assert_eq!(factor, saved);
        assert_eq!(result, assimilate_sequence(&factor, &observations, cx).unwrap());
        let dense = result.belief().to_belief(cx).unwrap();
        let oracle = assimilate_all(&prior, &observations, cx).unwrap();
        for i in 0..2 {
            close(dense.mean()[i], oracle.mean()[i]);
            for j in 0..2 { close(dense.covariance()[i][j], oracle.covariance()[i][j]); }
        }
        close(result.diagnostics().joint_misfit_before(), misfit(&prior, &observations, cx).unwrap());
        close(result.diagnostics().joint_misfit_after(), misfit(&dense, &observations, cx).unwrap());
    });
}

fn quadratic(v: [f64; 2], a: [[f64; 2]; 2], w: [f64; 2]) -> f64 {
    v[0] * (a[0][0] * w[0] + a[0][1] * w[1]) + v[1] * (a[1][0] * w[0] + a[1][1] * w[1])
}

#[test]
fn correlated_update_and_joint_misfit_match_direct_block_conditioning() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        for correlation in [-0.7, -0.1, 0.2, 0.7] {
            let p = [[2.0, 0.4], [0.4, 1.0]];
            let r = [[1.0, correlation], [correlation, 0.8]];
            let h = [[1.0, 0.2], [-0.3, 1.0]];
            let mean = [0.5, -1.0];
            let y = [1.0, 2.0];
            let prior = Belief::new(mean.to_vec(), p.iter().map(|row| row.to_vec()).collect(), cx).unwrap();
            let records = (0..2).map(|i| RobustObservation::available(
                Observation::new(h[i].to_vec(), y[i], r[i][i], format!("probe-{i}")).unwrap()
            )).collect();
            let batch = ObservationBatch::new(records, r.iter().map(|row| row.to_vec()).collect(), cx).unwrap();
            let result = assimilate_observation_batch(&prior, &batch, cx).unwrap();
            let posterior = result.posterior().unwrap();
            let s00 = quadratic(h[0], p, h[0]) + r[0][0];
            let s01 = quadratic(h[0], p, h[1]) + r[0][1];
            let s11 = quadratic(h[1], p, h[1]) + r[1][1];
            let det = s00 * s11 - s01 * s01;
            let inverse = [[s11 / det, -s01 / det], [-s01 / det, s00 / det]];
            let residual = [y[0] - h[0][0] * mean[0] - h[0][1] * mean[1], y[1] - h[1][0] * mean[0] - h[1][1] * mean[1]];
            let ph = [
                [p[0][0] * h[0][0] + p[0][1] * h[0][1], p[0][0] * h[1][0] + p[0][1] * h[1][1]],
                [p[1][0] * h[0][0] + p[1][1] * h[0][1], p[1][0] * h[1][0] + p[1][1] * h[1][1]],
            ];
            for i in 0..2 {
                close(posterior.mean()[i], mean[i] + quadratic(ph[i], inverse, residual));
                for j in 0..2 { close(posterior.covariance()[i][j], p[i][j] - quadratic(ph[i], inverse, ph[j])); }
            }
            let det_r = r[0][0] * r[1][1] - correlation * correlation;
            let inverse_r = [[r[1][1] / det_r, -correlation / det_r], [-correlation / det_r, r[0][0] / det_r]];
            let post_residual = [y[0] - h[0][0] * posterior.mean()[0] - h[0][1] * posterior.mean()[1], y[1] - h[1][0] * posterior.mean()[0] - h[1][1] * posterior.mean()[1]];
            let diagnostics = result.factor_diagnostics().unwrap();
            assert_eq!(diagnostics.receipts().len(), 2);
            close(diagnostics.joint_misfit_before(), quadratic(residual, inverse_r, residual));
            close(diagnostics.joint_misfit_after(), quadratic(post_residual, inverse_r, post_residual));
        }
    });
}

#[test]
fn correlated_whitening_rejects_both_short_and_long_operator_rows() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::diagonal(vec![0.0; 2], &[1.0; 2], cx).unwrap();
        for bad_dim in [1, 3] {
            let records = vec![
                RobustObservation::available(point_sensor(0, 2, 1.0, 1.0, "good").unwrap()),
                RobustObservation::available(point_sensor(0, bad_dim, 2.0, 1.0, "bad").unwrap()),
            ];
            let batch = ObservationBatch::new(records, vec![vec![1.0, 0.4], vec![0.4, 1.0]], cx).unwrap();
            assert_eq!(assimilate_observation_batch(&prior, &batch, cx), Err(AssimError::DimMismatch { state: 2, operator: bad_dim }));
        }
    });
}

#[test]
fn constrained_state_stays_exact_through_correlated_updates() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::diagonal(vec![0.0, 3.0], &[1.0, 0.0], cx).unwrap();
        let records = vec![
            RobustObservation::available(Observation::new(vec![1.0, 1.0], 4.0, 1.0, "a").unwrap()),
            RobustObservation::available(Observation::new(vec![0.5, -1.0], -2.0, 1.0, "b").unwrap()),
        ];
        let batch = ObservationBatch::new(records, vec![vec![1.0, 0.2], vec![0.2, 1.0]], cx).unwrap();
        let result = assimilate_observation_batch(&prior, &batch, cx).unwrap();
        let posterior = result.posterior().unwrap();
        assert_eq!(posterior.mean()[1], 3.0);
        assert_eq!(posterior.variance(1), Ok(0.0));
        assert_eq!(posterior.covariance()[0][1], 0.0);
    });
}

#[test]
fn factor_sequence_handles_batches_beyond_the_dense_cubic_envelope() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let prior = Belief::diagonal(vec![0.0; 128], &[1.0; 128], cx).unwrap();
        let observations: Vec<_> = (0..40).map(|i| point_sensor(i, 128, 1.0, 1.0, format!("sensor-{i}")).unwrap()).collect();
        assert!(matches!(assimilate_all(&prior, &observations, cx), Err(AssimError::AssimilationWorkLimit { .. })));
        let factor = FactorBelief::from_belief(&prior, cx).unwrap();
        let result = assimilate_sequence(&factor, &observations, cx).unwrap();
        for i in 0..128 {
            assert_eq!(result.belief().mean()[i], if i < 40 { 0.5 } else { 0.0 });
            assert_eq!(result.belief().variance(i), Some(if i < 40 { 0.5 } else { 1.0 }));
        }
    });
}

#[test]
fn sequence_does_not_reissue_scalar_poll_or_cost_budgets() {
    let gate = CancelGate::new();
    let (prior, observations) = with_cx(&gate, Budget::INFINITE, |cx| {
        let factor = FactorBelief::diagonal(vec![0.0], &[1.0], cx).unwrap();
        let rows = (0..16).map(|i| point_sensor(0, 1, 1.0, 1.0, format!("probe-{i}")).unwrap()).collect::<Vec<_>>();
        (factor, rows)
    });
    let saved = prior.clone();
    for budget in [Budget::INFINITE.with_poll_quota(4), Budget::INFINITE.with_cost_quota(1_000)] {
        with_cx(&gate, budget, |cx| {
            assimilate_scalar(&prior, &observations[0], cx).expect("one row fits");
            assert!(matches!(assimilate_sequence(&prior, &observations, cx), Err(AssimError::BudgetRefused(_)) | Err(AssimError::Cancelled { .. })));
        });
    }
    assert_eq!(prior, saved);
}

#[test]
fn state_independent_middle_row_is_kept_in_whitening_and_joint_diagnostics() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        // L = [[1,0,0],[1,1,0],[0,1,1]], so whitening H=[1,1,1]
        // gives [1,0,1], and whitening y=[2,5,7] gives [2,3,4].
        let records = [2.0, 5.0, 7.0].into_iter().enumerate().map(|(i, value)| {
            RobustObservation::available(point_sensor(0, 1, value, if i == 0 { 1.0 } else { 2.0 }, format!("row-{i}")).unwrap())
        }).collect();
        let batch = ObservationBatch::new(records, vec![vec![1.0, 1.0, 0.0], vec![1.0, 2.0, 1.0], vec![0.0, 1.0, 2.0]], cx).unwrap();
        let result = assimilate_observation_batch(&Belief::scalar(0.0, 1.0).unwrap(), &batch, cx).unwrap();
        let posterior = result.posterior().unwrap();
        close(posterior.mean()[0], 2.0);
        close(posterior.variance(0).unwrap(), 1.0 / 3.0);
        let diagnostics = result.factor_diagnostics().unwrap();
        assert_eq!(diagnostics.receipts().len(), 2);
        assert_eq!(diagnostics.updated_row_indices(), &[0, 2]);
        assert_eq!(diagnostics.state_independent_row_indices(), &[1]);
        assert_eq!(diagnostics.state_independent_misfit(), 9.0);
        close(diagnostics.joint_misfit_before(), 29.0);
        close(diagnostics.joint_misfit_after(), 13.0);
        // The audit counts actual readings, not informative whitened rows.
        assert_eq!(result.audit().effective_dof(), 3);
    });
}

#[test]
fn tiny_nonzero_whitened_sensitivity_is_not_thresholded_away() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let epsilon = 2.0_f64.powi(-40);
        let records = vec![
            RobustObservation::available(Observation::new(vec![1.0], 2.0, 1.0, "row-a").unwrap()),
            RobustObservation::available(Observation::new(vec![1.0 + epsilon], 5.0, 2.0, "row-b").unwrap()),
        ];
        let batch = ObservationBatch::new(records, vec![vec![1.0, 1.0], vec![1.0, 2.0]], cx).unwrap();
        let result = assimilate_observation_batch(&Belief::scalar(0.0, 1.0).unwrap(), &batch, cx).unwrap();
        let diagnostics = result.factor_diagnostics().unwrap();
        assert_eq!(diagnostics.receipts().len(), 2);
        assert_eq!(diagnostics.updated_row_indices(), &[0, 1]);
        assert!(diagnostics.state_independent_row_indices().is_empty());
        close(result.posterior().unwrap().mean()[0], (2.0 + 3.0 * epsilon) / (2.0 + epsilon * epsilon));
    });
}

#[test]
fn whitening_that_erases_all_coefficients_by_underflow_still_refuses() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let records = vec![
            RobustObservation::available(Observation::new(vec![f64::from_bits(1)], 0.0, 1.0e308, "tiny-row").unwrap()),
            RobustObservation::available(Observation::new(vec![1.0], 0.0, 1.0, "normal-row").unwrap()),
        ];
        let batch = ObservationBatch::new(records, vec![vec![1.0e308, 0.25], vec![0.25, 1.0]], cx).unwrap();
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        assert_eq!(assimilate_observation_batch(&prior, &batch, cx), Err(AssimError::ZeroObservationOperator));
    });
}

#[test]
fn overflow_in_a_state_independent_residual_never_publishes_a_posterior() {
    with_cx(&CancelGate::new(), Budget::INFINITE, |cx| {
        let records = vec![
            RobustObservation::available(point_sensor(0, 1, 0.0, 1.0, "row-a").unwrap()),
            RobustObservation::available(point_sensor(0, 1, f64::MAX, 2.0, "row-b").unwrap()),
        ];
        let batch = ObservationBatch::new(records, vec![vec![1.0, 1.0], vec![1.0, 2.0]], cx).unwrap();
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        assert!(matches!(assimilate_observation_batch(&prior, &batch, cx), Err(AssimError::NonFiniteComputation { .. })));
    });
}
