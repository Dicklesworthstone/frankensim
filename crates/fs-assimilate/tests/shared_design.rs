//! Structural-factor and multi-measurement design regressions.
use fs_assimilate::factor::{FactorBelief, assimilate_scalar};
use fs_assimilate::groups::RowId;
use fs_assimilate::groups::jacobian::{JacobianBatch, JacobianRecord, SharedNoiseModel};
use fs_assimilate::groups::stream::SharedNoiseStream;
use fs_assimilate::{Belief, Observation};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};

fn with_budget<R>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 6420, kernel_id: 20, tile: 0, iteration: 0 },
            budget, ExecMode::Deterministic).with_time_source(&clock);
        f(&cx)
    });
    assert!(pool.stats().quiescent());
    result
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2.0e-10 * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn independent_product_copies_correlated_blocks_and_exact_zeros() {
    with_budget(Budget::INFINITE, |cx| {
        let a = Belief::new(vec![1.0, -2.0], vec![vec![2.0, 0.5], vec![0.5, 1.0]], cx).unwrap();
        let b = Belief::new(vec![3.0, 4.0], vec![vec![0.25, 0.25], vec![0.25, 0.25]], cx).unwrap();
        let fa = FactorBelief::from_belief(&a, cx).unwrap();
        let fb = FactorBelief::from_belief(&b, cx).unwrap();
        let joint = fa.independent_product(&fb, cx).unwrap();
        assert_eq!(joint.mean(), &[1.0, -2.0, 3.0, 4.0]);
        for (start, block) in [(0, &fa), (2, &fb)] {
            for i in 0..2 {
                assert_eq!(joint.diag(start + i), block.diag(i));
                for j in i + 1..2 { assert_eq!(joint.upper_entry(start + i, start + j), block.upper_entry(i, j)); }
            }
        }
        let dense = joint.to_dense_covariance();
        for i in 0..4 { for j in 0..4 {
            let expected = if i < 2 && j < 2 { a.covariance()[i][j] }
                else if i >= 2 && j >= 2 { b.covariance()[i - 2][j - 2] } else { 0.0 };
            assert_eq!(dense[i][j], expected);
        } }
        assert_eq!(joint.prefix_marginal(2, cx).unwrap(), a);
    });
}

#[test]
fn singular_shared_stream_publishes_physical_marginal_before_and_after_updates() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("singular-fit", vec!["a".into(), "b".into()],
            vec![vec![0.25, 0.25], vec![0.25, 0.25]], cx).unwrap();
        let prior = Belief::scalar(0.0, 2.0).unwrap();
        let mut stream = SharedNoiseStream::new(&prior, &model, 2, cx).unwrap();
        assert_eq!(stream.physical_belief(cx).unwrap(), prior);
        for i in 0..2 {
            let row = JacobianRecord::new(Observation::new(vec![1.0], 1.0, 0.25, "sensor").unwrap(),
                RowId::try_new(format!("row-{i}")).unwrap(), &model, vec![1.0, -1.0]).unwrap();
            stream.assimilate(&JacobianBatch::new(model.clone(), vec![row], cx).unwrap(), cx).unwrap();
            close(stream.physical_belief(cx).unwrap().variance(0).unwrap(), 1.0 / (0.5 + 4.0 * (i + 1) as f64));
            assert_eq!(stream.joint().variance(1), Some(0.25));
            assert_eq!(stream.joint().variance(2), Some(0.25));
        }
    });
}

#[test]
fn marginal_uses_tail_factor_columns_and_operations_refuse_exhausted_budgets() {
    let joint = with_budget(Budget::INFINITE, |cx| {
        let a = FactorBelief::diagonal(vec![0.0], &[2.0], cx).unwrap();
        let b = FactorBelief::diagonal(vec![0.0], &[1.0], cx).unwrap();
        let joint = a.independent_product(&b, cx).unwrap();
        let observation = Observation::new(vec![1.0, 1.0], 3.0, 0.25, "sensor").unwrap();
        let updated = assimilate_scalar(&joint, &observation, cx).unwrap();
        let joint = updated.belief().clone();
        let marginal = joint.prefix_marginal(1, cx).unwrap();
        close(marginal.mean()[0], 6.0 / 3.25);
        close(marginal.variance(0).unwrap(), 2.0 - 4.0 / 3.25);
        assert!(marginal.variance(0).unwrap() > joint.diag(0).unwrap());
        assert!(joint.prefix_marginal(0, cx).is_err());
        assert!(joint.prefix_marginal(3, cx).is_err());
        joint
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            let before = joint.clone();
            assert!(joint.independent_product(&joint, cx).is_err());
            assert!(joint.prefix_marginal(1, cx).is_err());
            assert_eq!(joint, before);
        });
    }
}

#[path = "shared_design/batch.rs"]
mod batch;
