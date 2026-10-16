//! G0/G4/G5: actual-reading recovery, bounded untrusted bytes, and numerical replay.
use fs_blake3::DomainHasher;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::measurement::{
    MeasurementAlternative, MeasurementDecision, MeasurementError, MeasurementStudy,
};
use fs_oed_e2e::plan::measurement::checkpoint::MAX_MEASUREMENT_CHECKPOINT_BYTES;
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::{Dims, QUANTITY_SPEC_ENCODED_LEN, QtyAny};
use fs_qty::semantic::{QuantityKind, SemanticQty, SemanticType, ValueForm};

fn with_cx<R>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey { seed: 920, kernel_id: 20, tile: 0, iteration: 0 },
            budget, ExecMode::Deterministic).with_time_source(&clock))
    });
    assert!(pool.stats().quiescent());
    result
}
fn value(x: f64) -> ObjectiveValue { ObjectiveValue::dimensionless(x).unwrap() }
fn menu() -> Vec<MeasurementAlternative> {
    [("A", 0.0, 4.0), ("B", 0.1, 9.0)].into_iter().map(|(name, mean, variance)| {
        MeasurementAlternative::new(name, value(mean), QtyAny::dimensionless(variance),
            QtyAny::dimensionless(1.0), AcquisitionCost::relative(1.0).unwrap()).unwrap()
    }).collect()
}
fn initialized(cx: &Cx<'_>) -> MeasurementStudy {
    MeasurementStudy::new(menu(), value(0.0), 8, cx).unwrap()
}
fn advance(study: &mut MeasurementStudy, count: usize, cx: &Cx<'_>) {
    for _ in 0..count {
        let MeasurementDecision::Acquire(request) = study.recommend(cx).unwrap() else { panic!("expected request") };
        let reading = request.reading(value(0.2), "actual-meter").unwrap();
        study.assimilate(&request, &reading, cx).unwrap();
    }
}
fn same_science(a: &MeasurementStudy, b: &MeasurementStudy) {
    assert_eq!(a.history_id(), b.history_id());
    assert_eq!(a.best_alternative(), b.best_alternative());
    assert_eq!(a.evpi(), b.evpi());
    assert_eq!(a.updates().len(), b.updates().len());
    for name in ["A", "B"] {
        assert_eq!(a.posterior_mean(name), b.posterior_mean(name));
        assert_eq!(a.posterior_variance(name), b.posterior_variance(name));
    }
    for (left, right) in a.updates().iter().zip(b.updates()) {
        assert_eq!(left.request(), right.request());
        assert_eq!(left.outcome(), right.outcome());
        assert_eq!(left.posterior_mean(), right.posterior_mean());
        assert_eq!(left.posterior_variance(), right.posterior_variance());
    }
    // Execution colors may change with the recovery context's budget/seed.
    // They are recomputed evidence, not serialized claims about the old run.
}
fn reseal(bytes: &mut [u8]) {
    let split = bytes.len() - 32;
    let mut hash = DomainHasher::new("org.frankensim.sensorforge.measurement-checkpoint.v1");
    hash.update(&bytes[..split]);
    bytes[split..].copy_from_slice(hash.finalize().as_bytes());
}
fn first_record() -> usize {
    24 + 2 * QUANTITY_SPEC_ENCODED_LEN + 1 + 8 + 2 * (2 + 1 + 32)
}

#[test]
fn initial_and_measured_studies_round_trip_and_continue_identically() {
    with_cx(Budget::INFINITE, |cx| {
        for count in [0, 1, 4] {
            let mut original = initialized(cx);
            advance(&mut original, count, cx);
            let bytes = original.to_checkpoint(cx).unwrap();
            let mut restored = MeasurementStudy::from_checkpoint(&bytes, cx).unwrap();
            same_science(&original, &restored);
            assert_eq!(bytes, restored.to_checkpoint(cx).unwrap());
            assert_eq!(original.recommend(cx), restored.recommend(cx));
            advance(&mut original, 1, cx);
            advance(&mut restored, 1, cx);
            same_science(&original, &restored);
        }
    });
}

#[test]
fn canonical_menu_order_and_cost_roles_survive_storage() {
    with_cx(Budget::INFINITE, |cx| {
        let original = initialized(cx);
        let mut reversed = menu();
        reversed.reverse();
        let reversed = MeasurementStudy::new(reversed, value(0.0), 8, cx).unwrap();
        assert_eq!(original.to_checkpoint(cx), reversed.to_checkpoint(cx));
        let physical_cost = menu().into_iter().enumerate().map(|(i, a)| {
            MeasurementAlternative::new(a.name(), value(i as f64 * 0.1), QtyAny::dimensionless(4.0),
                QtyAny::dimensionless(1.0), AcquisitionCost::dimensional(QtyAny::dimensionless(1.0)).unwrap()).unwrap()
        }).collect();
        let physical = MeasurementStudy::new(physical_cost, value(0.0), 8, cx).unwrap();
        let restored = MeasurementStudy::from_checkpoint(&physical.to_checkpoint(cx).unwrap(), cx).unwrap();
        assert_ne!(restored.history_id(), original.history_id());
        assert_eq!(restored.recommend(cx), physical.recommend(cx));
    });
}

#[test]
fn semantic_objectives_and_physical_costs_are_revalidated() {
    with_cx(Budget::INFINITE, |cx| {
        let absolute = SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static);
        let make = |x| ObjectiveValue::semantic(SemanticQty::new(QtyAny::new(x, absolute.expected_dims()), absolute).unwrap());
        let initial = make(300.0);
        let variance = initial.spec().variance_dims().unwrap();
        let cost_dims = Dims([0, 0, 1, 0, 0, 0]);
        let rows = [("A", 300.0), ("B", 300.1)].into_iter().map(|(name, mean)| {
            MeasurementAlternative::new(name, make(mean), QtyAny::new(4.0, variance),
                QtyAny::new(1.0, variance), AcquisitionCost::dimensional(QtyAny::new(2.0, cost_dims)).unwrap()).unwrap()
        }).collect();
        let mut study = MeasurementStudy::new(rows, initial.spec().decision_value(0.0).unwrap(), 2, cx).unwrap();
        let MeasurementDecision::Acquire(request) = study.recommend(cx).unwrap() else { panic!("expected request") };
        study.assimilate(&request, &request.reading(make(300.2), "thermometer").unwrap(), cx).unwrap();
        let restored = MeasurementStudy::from_checkpoint(&study.to_checkpoint(cx).unwrap(), cx).unwrap();
        same_science(&study, &restored);
        assert_eq!(restored.evpi().spec(), initial.spec().decision_spec());
        assert_eq!(restored.posterior_mean("A").unwrap().spec(), initial.spec());
    });
}

#[test]
fn truncation_version_checksum_trailing_bytes_and_count_bombs_refuse() {
    with_cx(Budget::INFINITE, |cx| {
        let mut study = initialized(cx);
        advance(&mut study, 1, cx);
        let bytes = study.to_checkpoint(cx).unwrap();
        for end in [0, 8, 24, bytes.len() - 1] {
            assert!(MeasurementStudy::from_checkpoint(&bytes[..end], cx).is_err());
        }
        for offset in [0, 8, 24, bytes.len() - 1] {
            let mut changed = bytes.clone(); changed[offset] ^= 1;
            assert!(MeasurementStudy::from_checkpoint(&changed, cx).is_err());
        }
        let mut trailing = bytes.clone(); trailing.insert(trailing.len() - 32, 0); reseal(&mut trailing);
        assert!(MeasurementStudy::from_checkpoint(&trailing, cx).is_err());
        for offset in [12, 16, 20] {
            let mut changed = bytes.clone(); changed[offset..offset + 4].copy_from_slice(&u32::MAX.to_le_bytes());
            reseal(&mut changed);
            assert!(MeasurementStudy::from_checkpoint(&changed, cx).is_err());
        }
        assert!(MeasurementStudy::from_checkpoint(&vec![0; MAX_MEASUREMENT_CHECKPOINT_BYTES + 1], cx).is_err());
    });
}

#[test]
fn a_new_checksum_does_not_make_altered_readings_or_requests_replayable() {
    with_cx(Budget::INFINITE, |cx| {
        let mut study = initialized(cx); advance(&mut study, 2, cx);
        let bytes = study.to_checkpoint(cx).unwrap();
        let record = first_record();
        for offset in [record, bytes.len() - 64] {
            let mut changed = bytes.clone(); changed[offset] ^= 1; reseal(&mut changed);
            assert!(MeasurementStudy::from_checkpoint(&changed, cx).is_err());
        }
        for number in [0.5, f64::NAN, f64::INFINITY, -0.0] {
            let mut changed = bytes.clone();
            changed[record + 32..record + 40].copy_from_slice(&number.to_le_bytes());
            reseal(&mut changed);
            assert!(MeasurementStudy::from_checkpoint(&changed, cx).is_err());
        }
        assert_eq!(bytes, study.to_checkpoint(cx).unwrap());
        same_science(&study, &MeasurementStudy::from_checkpoint(&bytes, cx).unwrap());
    });
}

#[test]
fn recovery_shares_poll_and_work_budgets_across_the_entire_history() {
    let (one, four) = with_cx(Budget::INFINITE, |cx| {
        let mut study = initialized(cx); advance(&mut study, 1, cx);
        let one = study.to_checkpoint(cx).unwrap(); advance(&mut study, 3, cx);
        (one, study.to_checkpoint(cx).unwrap())
    });
    let mut quota = 1;
    while !with_cx(Budget::INFINITE.with_poll_quota(quota), |cx| MeasurementStudy::from_checkpoint(&one, cx).is_ok()) {
        quota *= 2; assert!(quota <= 16384);
    }
    assert!(with_cx(Budget::INFINITE.with_poll_quota(quota), |cx| MeasurementStudy::from_checkpoint(&four, cx)).is_err());
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(1)] {
        assert!(with_cx(budget, |cx| MeasurementStudy::from_checkpoint(&one, cx)).is_err());
    }
    assert!(with_cx(Budget::INFINITE, |cx| MeasurementStudy::from_checkpoint(&four, cx)).is_ok());
}

#[test]
fn cancellation_preserves_the_source_and_saved_data_for_retry() {
    let study = with_cx(Budget::INFINITE, |cx| {
        let mut study = initialized(cx); advance(&mut study, 2, cx); study
    });
    let bytes = with_cx(Budget::INFINITE, |cx| study.to_checkpoint(cx).unwrap());
    let before = study.clone();
    let gate = CancelGate::new(); gate.request();
    let clock = VirtualClock::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 920, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        assert!(study.to_checkpoint(&cx).is_err());
        assert!(MeasurementStudy::from_checkpoint(&bytes, &cx).is_err());
    });
    assert!(pool.stats().quiescent()); assert_eq!(study, before);
    with_cx(Budget::INFINITE, |cx| same_science(&study, &MeasurementStudy::from_checkpoint(&bytes, cx).unwrap()));
}

use fs_oed_e2e::plan::acquisition::{AcquiredMeasurement, AcquisitionError, MeasurementProvider};
use fs_oed_e2e::plan::measurement::{MeasurementOutcome, MeasurementRequest};

struct Meter { calls: usize }
impl MeasurementProvider for Meter {
    type Error = MeasurementError;
    fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
        self.calls += 1;
        request.reading(value(0.2), "actual-meter")
    }
}
fn request(study: &MeasurementStudy, cx: &Cx<'_>) -> MeasurementRequest {
    let MeasurementDecision::Acquire(request) = study.recommend(cx).unwrap() else { panic!("expected request") };
    request
}

#[test]
fn acquired_but_unapplied_readings_survive_restart_without_another_provider_call() {
    with_cx(Budget::INFINITE, |cx| {
        let mut meter = Meter { calls: 0 };
        let mut original = initialized(cx);
        advance(&mut original, 1, cx);
        let source = original.to_checkpoint(cx).unwrap();
        let acquired = original.acquire(&request(&original, cx), &mut meter, cx).unwrap();
        let pending = acquired.to_checkpoint(cx).unwrap();
        let expected = acquired.apply(&mut original, cx).unwrap().clone();
        drop(original);
        let mut restored = MeasurementStudy::from_checkpoint(&source, cx).unwrap();
        let recovered = AcquiredMeasurement::from_checkpoint(&restored, &pending, cx).unwrap();
        assert_eq!(recovered, acquired);
        let applied = recovered.apply(&mut restored, cx).unwrap();
        assert_eq!(applied.request(), expected.request());
        assert_eq!(applied.outcome(), expected.outcome());
        assert_eq!(applied.posterior_variance(), expected.posterior_variance());
        assert_eq!(meter.calls, 1);
        assert!(matches!(AcquiredMeasurement::from_checkpoint(&restored, &pending, cx), Err(MeasurementError::StaleRequest)));
        assert_eq!(restored.updates().len(), 2);
    });
}

#[test]
fn post_io_cancellation_can_be_saved_and_recovered_before_numerical_application() {
    struct CancellingMeter<'a> { gate: &'a CancelGate, calls: usize }
    impl MeasurementProvider for CancellingMeter<'_> {
        type Error = MeasurementError;
        fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
            self.calls += 1;
            let reading = request.reading(value(0.3), "cancel-after-read")?;
            self.gate.request();
            Ok(reading)
        }
    }
    let (study, planned) = with_cx(Budget::INFINITE, |cx| { let study = initialized(cx); let planned = request(&study, cx); (study, planned) });
    let gate = CancelGate::new();
    let mut meter = CancellingMeter { gate: &gate, calls: 0 };
    let clock = VirtualClock::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let retained = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 920, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        match study.acquire(&planned, &mut meter, &cx) {
            Err(AcquisitionError::After { acquired, .. }) => acquired,
            other => panic!("expected acquired reading after cancellation: {other:?}"),
        }
    });
    assert!(pool.stats().quiescent()); assert_eq!(meter.calls, 1);
    assert!(study.updates().is_empty());
    with_cx(Budget::INFINITE, |cx| {
        let source = study.to_checkpoint(cx).unwrap();
        let pending = retained.to_checkpoint(cx).unwrap();
        let mut restored = MeasurementStudy::from_checkpoint(&source, cx).unwrap();
        AcquiredMeasurement::from_checkpoint(&restored, &pending, cx).unwrap().apply(&mut restored, cx).unwrap();
        assert_eq!(restored.updates().len(), 1);
        assert_eq!(restored.updates()[0].outcome().instrument(), "cancel-after-read");
    });
    assert_eq!(meter.calls, 1);
}

#[test]
fn pending_transport_refuses_foreign_history_corruption_and_cross_type_frames() {
    with_cx(Budget::INFINITE, |cx| {
        let study = initialized(cx);
        let mut meter = Meter { calls: 0 };
        let pending = study.acquire(&request(&study, cx), &mut meter, cx).unwrap().to_checkpoint(cx).unwrap();
        let foreign = MeasurementStudy::new(menu(), value(0.0), 7, cx).unwrap();
        assert!(matches!(AcquiredMeasurement::from_checkpoint(&foreign, &pending, cx), Err(MeasurementError::StaleRequest)));
        assert!(MeasurementStudy::from_checkpoint(&pending, cx).is_err());
        assert!(AcquiredMeasurement::from_checkpoint(&study, &study.to_checkpoint(cx).unwrap(), cx).is_err());
        for end in [0, 12, pending.len() - 1] {
            assert!(AcquiredMeasurement::from_checkpoint(&study, &pending[..end], cx).is_err());
        }
        for offset in [0, 8, 12, pending.len() - 1] {
            let mut changed = pending.clone(); changed[offset] ^= 1;
            assert!(AcquiredMeasurement::from_checkpoint(&study, &changed, cx).is_err());
        }
        assert_eq!(meter.calls, 1); assert!(study.updates().is_empty());
    });
}

#[test]
fn pending_data_rechecks_schema_and_current_request_even_with_a_replaced_checksum() {
    with_cx(Budget::INFINITE, |cx| {
        let study = initialized(cx);
        let pending = study.acquire(&request(&study, cx), &mut Meter { calls: 0 }, cx).unwrap().to_checkpoint(cx).unwrap();
        for offset in [44, 76 + 1] { // request identity, then one dimension exponent
            let mut changed = pending.clone(); changed[offset] ^= 1;
            let split = changed.len() - 32;
            let mut hash = DomainHasher::new("org.frankensim.sensorforge.acquired-measurement.v1");
            hash.update(&changed[..split]); changed[split..].copy_from_slice(hash.finalize().as_bytes());
            assert!(AcquiredMeasurement::from_checkpoint(&study, &changed, cx).is_err());
        }
        assert!(study.updates().is_empty());
    });
}

#[test]
fn pending_recovery_budget_failure_keeps_reading_available_for_retry() {
    let (study, pending) = with_cx(Budget::INFINITE, |cx| {
        let study = initialized(cx);
        let pending = study.acquire(&request(&study, cx), &mut Meter { calls: 0 }, cx).unwrap().to_checkpoint(cx).unwrap();
        (study, pending)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(1), Budget::INFINITE.with_cost_quota(1)] {
        assert!(with_cx(budget, |cx| AcquiredMeasurement::from_checkpoint(&study, &pending, cx)).is_err());
    }
    let mut restored = study.clone();
    with_cx(Budget::INFINITE, |cx| {
        AcquiredMeasurement::from_checkpoint(&restored, &pending, cx).unwrap().apply(&mut restored, cx).unwrap();
    });
    assert!(study.updates().is_empty()); assert_eq!(restored.updates().len(), 1);
}
