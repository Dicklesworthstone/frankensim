use super::*;
use fs_assimilate::groups::stream::{MAX_STREAM_ROWS, SharedNoiseStream};

fn offset(cx: &Cx<'_>, name: &str, variance: f64) -> SharedNoiseModel {
    SharedNoiseModel::new(name, vec!["offset".into()], vec![vec![variance]], cx).unwrap()
}

fn batch(model: &SharedNoiseModel, rows: Vec<JacobianRecord>, cx: &Cx<'_>) -> JacobianBatch {
    JacobianBatch::new(model.clone(), rows, cx).unwrap()
}

#[test]
fn shared_error_floor_survives_thirty_two_separate_arrivals() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "one-fixed-calibration", 4.0);
        let prior = Belief::scalar(0.0, 10.0).unwrap();
        let mut stream = SharedNoiseStream::new(&prior, &model, 32, cx).unwrap();
        let mut incorrectly_reset = prior.clone();
        for i in 0..32 {
            let arriving = batch(&model, vec![record(&model, &format!("row-{i}"), vec![1.0])], cx);
            let update = stream.assimilate(&arriving, cx).unwrap();
            assert_eq!(update.accepted_indices(), &[0]);
            assert_eq!(update.diagnostics().unwrap().receipts().len(), 1);
            incorrectly_reset = arriving.assimilate(&incorrectly_reset, cx).unwrap().robust().posterior().unwrap().clone();
            let count = (i + 1) as f64;
            let expected = 1.0 / (0.1 + count / (0.25 + 4.0 * count));
            close(stream.joint().variance(0).unwrap(), expected);
            close(stream.joint().mean()[0], expected * count / (0.25 + 4.0 * count));
        }
        assert_eq!(stream.accepted_rows(), 32);
        assert!(stream.joint().variance(0).unwrap() > 1.0 / (0.1 + 0.25));
        assert!(incorrectly_reset.variance(0).unwrap() < stream.joint().variance(0).unwrap() / 10.0);
        close(stream.physical_belief(cx).unwrap().variance(0).unwrap(), stream.joint().variance(0).unwrap());
    });
}

#[test]
fn split_and_join_match_full_covariance_conditioning_and_joint_factors() {
    with_budget(Budget::INFINITE, |cx| {
        let prior = Belief::new(vec![0.2, -0.3], vec![vec![2.0, 0.4], vec![0.4, 1.0]], cx).unwrap();
        let model = SharedNoiseModel::new("offset-tilt", vec!["offset".into(), "tilt".into()], vec![vec![0.8, 0.2], vec![0.2, 0.5]], cx).unwrap();
        let rows: Vec<_> = (0..6).map(|i| JacobianRecord::new(
            Observation::new(vec![1.0, (i as f64 - 2.0) * 0.25], i as f64 * 0.3, 0.2, "instrument").unwrap(),
            RowId::try_new(format!("row-{i}")).unwrap(), &model, vec![1.0, i as f64 - 3.0],
        ).unwrap()).collect();
        let full = batch(&model, rows.clone(), cx);
        let reference = full.assimilate(&prior, cx).unwrap();
        let mut joined = SharedNoiseStream::new(&prior, &model, 8, cx).unwrap();
        joined.assimilate(&full, cx).unwrap();
        for chunk_size in [1, 2, 3] {
            let mut split = SharedNoiseStream::new(&prior, &model, 8, cx).unwrap();
            for chunk in rows.chunks(chunk_size) {
                split.assimilate(&batch(&model, chunk.to_vec(), cx), cx).unwrap();
            }
            assert_eq!(split.joint(), joined.joint());
            assert_eq!(split.history_id(), joined.history_id());
            assert_eq!(split, joined);
        }
        let marginal = joined.physical_belief(cx).unwrap();
        let expected = reference.robust().posterior().unwrap();
        for i in 0..2 {
            close(marginal.mean()[i], expected.mean()[i]);
            for j in 0..2 { close(marginal.covariance()[i][j], expected.covariance()[i][j]); }
        }
    });
}

#[test]
fn exact_replays_are_noops_and_mixed_arrivals_only_update_new_rows() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let a = record(&model, "row-a", vec![1.0]);
        let b = record(&model, "row-b", vec![-1.0]);
        let first = batch(&model, vec![a.clone()], cx);
        let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 2, cx).unwrap();
        stream.assimilate(&first, cx).unwrap();
        let before = stream.clone();
        let replay = stream.assimilate(&first, cx).unwrap();
        assert_eq!(replay.replayed_indices(), &[0]);
        assert!(replay.accepted_indices().is_empty());
        assert!(replay.diagnostics().is_none());
        assert_eq!(stream, before);
        let update = stream.assimilate(&batch(&model, vec![a, b.clone()], cx), cx).unwrap();
        assert_eq!(update.accepted_indices(), &[1]);
        assert_eq!(update.replayed_indices(), &[0]);
        let mut reference = before;
        reference.assimilate(&batch(&model, vec![b], cx), cx).unwrap();
        assert_eq!(stream, reference);
        // Replaying does not need another row-cap slot, even at capacity.
        stream.assimilate(&first, cx).unwrap();
        assert_eq!(stream, reference);
    });
}

#[test]
fn changed_reading_noise_operator_jacobian_or_instrument_refuses_atomically() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 4, cx).unwrap();
        stream.assimilate(&batch(&model, vec![record(&model, "old", vec![1.0])], cx), cx).unwrap();
        let original = stream.clone();
        for (h, value, noise, j, instrument) in [
            (1.0, 2.0, 0.25, 1.0, "scan-instrument"),
            (1.0, 1.0, 0.5, 1.0, "scan-instrument"),
            (2.0, 1.0, 0.25, 1.0, "scan-instrument"),
            (1.0, 1.0, 0.25, -1.0, "scan-instrument"),
            (1.0, 1.0, 0.25, 1.0, "different-instrument"),
        ] {
            let changed = JacobianRecord::new(Observation::new(vec![h], value, noise, instrument).unwrap(), RowId::try_new("old").unwrap(), &model, vec![j]).unwrap();
            let incoming = batch(&model, vec![record(&model, "new", vec![1.0]), changed], cx);
            assert!(matches!(stream.assimilate(&incoming, cx), Err(AssimError::DuplicateDatasetRow { .. })));
            assert_eq!(stream, original);
        }
    });
}

#[test]
fn model_dimension_and_capacity_refusals_preserve_the_stream() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let prior = Belief::scalar(0.0, 1.0).unwrap();
        assert!(SharedNoiseStream::new(&prior, &model, 0, cx).is_err());
        assert!(SharedNoiseStream::new(&prior, &model, MAX_STREAM_ROWS + 1, cx).is_err());
        let mut stream = SharedNoiseStream::new(&prior, &model, 1, cx).unwrap();
        let original = stream.clone();
        let other = offset(cx, "different-fit", 1.0);
        assert!(matches!(stream.assimilate(&batch(&other, vec![record(&other, "a", vec![1.0])], cx), cx), Err(AssimError::SharedSourceDeclaration { .. })));
        let wrong_dim = JacobianRecord::new(point_sensor(0, 2, 0.0, 0.25, "sensor").unwrap(), RowId::try_new("b").unwrap(), &model, vec![1.0]).unwrap();
        assert!(matches!(stream.assimilate(&batch(&model, vec![wrong_dim], cx), cx), Err(AssimError::DimMismatch { .. })));
        let over_cap = batch(&model, vec![record(&model, "a", vec![1.0]), record(&model, "b", vec![1.0])], cx);
        assert!(matches!(stream.assimilate(&over_cap, cx), Err(AssimError::ObservationCountLimit { .. })));
        assert_eq!(stream, original);
    });
}

#[test]
fn failed_factor_update_and_exhausted_quotas_never_consume_rows() {
    let (stream, incoming) = with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 4, cx).unwrap();
        let incoming = batch(&model, vec![record(&model, "a", vec![1.0]), record(&model, "b", vec![1.0])], cx);
        let huge = JacobianRecord::new(Observation::new(vec![1e200], 1.0, 0.25, "large-operator").unwrap(), RowId::try_new("huge").unwrap(), &model, vec![1.0]).unwrap();
        let mut candidate = stream.clone();
        assert!(candidate.assimilate(&batch(&model, vec![record(&model, "good", vec![1.0]), huge], cx), cx).is_err());
        assert_eq!(candidate, stream);
        (stream, incoming)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(1), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            let mut candidate = stream.clone();
            assert!(candidate.assimilate(&incoming, cx).is_err());
            assert_eq!(candidate, stream);
        });
    }
    with_budget(Budget::INFINITE, |cx| {
        let mut retry = stream;
        retry.assimilate(&incoming, cx).unwrap();
        assert_eq!(retry.accepted_rows(), 2);
    });
}

#[test]
fn zero_nuisance_and_constrained_physical_components_remain_exact() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "zero-error", 0.0);
        let prior = Belief::diagonal(vec![3.0, 0.0], &[0.0, 2.0], cx).unwrap();
        let mut stream = SharedNoiseStream::new(&prior, &model, 2, cx).unwrap();
        for i in 0..2 {
            let row = JacobianRecord::new(Observation::new(vec![1.0, 1.0], 5.0, 0.25, "sensor").unwrap(), RowId::try_new(format!("row-{i}")).unwrap(), &model, vec![100.0]).unwrap();
            stream.assimilate(&batch(&model, vec![row], cx), cx).unwrap();
        }
        assert_eq!(stream.joint().mean()[0], 3.0);
        assert_eq!(stream.joint().variance(0), Some(0.0));
        assert_eq!(stream.joint().variance(2), Some(0.0));
        close(stream.joint().variance(1).unwrap(), 1.0 / (0.5 + 2.0 / 0.25));
    });
}

#[test]
fn explicit_cancellation_leaves_the_retained_state_untouched() {
    let (mut stream, incoming) = with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 2, cx).unwrap();
        let incoming = batch(&model, vec![record(&model, "row", vec![1.0])], cx);
        (stream, incoming)
    });
    let before = stream.clone();
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let clock = VirtualClock::new();
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 6416, kernel_id: 16, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        assert!(stream.assimilate(&incoming, &cx).is_err());
    });
    assert!(pool.stats().quiescent());
    assert_eq!(stream, before);
}

#[test]
fn admitted_singular_nuisance_keeps_its_null_direction() {
    with_budget(Budget::INFINITE, |cx| {
        let model = SharedNoiseModel::new("singular-fit", vec!["a".into(), "b".into()], vec![vec![0.25, 0.25], vec![0.25, 0.25]], cx).unwrap();
        let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 2, cx).unwrap();
        for i in 0..2 {
            stream.assimilate(&batch(&model, vec![record(&model, &format!("row-{i}"), vec![1.0, -1.0])], cx), cx).unwrap();
        }
        close(stream.joint().variance(0).unwrap(), 1.0 / (0.5 + 2.0 / 0.25));
        assert_eq!(stream.joint().mean()[1], 0.0);
        assert_eq!(stream.joint().mean()[2], 0.0);
        assert_eq!(stream.joint().variance(1), Some(0.25));
        assert_eq!(stream.joint().variance(2), Some(0.25));
    });
}
