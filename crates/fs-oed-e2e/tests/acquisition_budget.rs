//! G0/G4/G5: actual acquisition caps and durable pre-I/O commitments.
use fs_blake3::DomainHasher;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::acquisition::{AcquiredMeasurement, AcquisitionError, MeasurementProvider};
use fs_oed_e2e::plan::acquisition::budget::{AcquisitionAccount, AcquisitionLimitError};
use fs_oed_e2e::plan::measurement::{MeasurementAlternative, MeasurementDecision, MeasurementError,
    MeasurementOutcome, MeasurementRequest, MeasurementStudy};
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::{Dims, QtyAny};

fn with_cx<R>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new(); let clock = VirtualClock::new();
    let result = pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 942, kernel_id: 20, tile: 0, iteration: 0 },
        budget, ExecMode::Deterministic).with_time_source(&clock)));
    assert!(pool.stats().quiescent()); result
}
fn value(x: f64) -> ObjectiveValue { ObjectiveValue::dimensionless(x).unwrap() }
fn cost(x: f64) -> AcquisitionCost { AcquisitionCost::relative(x).unwrap() }
fn study(cost: AcquisitionCost, cx: &Cx<'_>) -> MeasurementStudy {
    let alternatives = [("A", 0.0), ("B", 0.1)].into_iter().map(|(name, mean)|
        MeasurementAlternative::new(name, value(mean), QtyAny::dimensionless(4.0),
            QtyAny::dimensionless(1.0), cost).unwrap()).collect();
    MeasurementStudy::new(alternatives, value(0.0), 8, cx).unwrap()
}
fn request(study: &MeasurementStudy, cx: &Cx<'_>) -> MeasurementRequest {
    let MeasurementDecision::Acquire(request) = study.recommend(cx).unwrap() else { panic!("expected acquisition") }; request
}
#[derive(Default)]
struct Meter { calls: usize, fail: bool }
impl MeasurementProvider for Meter {
    type Error = MeasurementError;
    fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
        self.calls += 1;
        if self.fail { Err(MeasurementError::Invalid("instrument failed")) }
        else { request.reading(value(0.2), "actual-meter") }
    }
}

#[test]
fn cumulative_cost_stops_io_without_changing_the_scientific_study() {
    with_cx(Budget::INFINITE, |cx| {
        let mut study = study(cost(2.0), cx);
        let mut account = AcquisitionAccount::new(cost(3.0), 8).unwrap();
        let mut meter = Meter::default();
        let permit = account.reserve(&study, &request(&study, cx), "attempt-0", cx).unwrap();
        assert_eq!(account.committed_cost_upper(), 2.0);
        assert_eq!(meter.calls, 0); // Reservation alone never acquires or assimilates.
        let acquired = permit.acquire(&study, &mut meter, cx).unwrap();
        acquired.apply(&mut study, cx).unwrap();
        let before = study.clone();
        assert!(matches!(account.reserve(&study, &request(&study, cx), "attempt-1", cx), Err(AcquisitionLimitError::CostLimit { .. })));
        assert_eq!(meter.calls, 1); assert_eq!(study, before);
        assert_eq!(account.charges().len(), 1);
    });
}

#[test]
fn failed_attempts_are_not_refunded_and_retries_need_new_ids_and_allowance() {
    with_cx(Budget::INFINITE, |cx| {
        let study = study(cost(2.0), cx); let planned = request(&study, cx);
        let mut account = AcquisitionAccount::new(cost(4.0), 2).unwrap();
        let mut meter = Meter { calls: 0, fail: true };
        let permit = account.reserve(&study, &planned, "attempt-0", cx).unwrap();
        assert!(matches!(permit.acquire(&study, &mut meter, cx), Err(AcquisitionError::Provider { .. })));
        assert!(matches!(account.reserve(&study, &planned, "attempt-0", cx), Err(AcquisitionLimitError::DuplicateAttempt { .. })));
        assert_eq!(account.committed_cost_upper(), 2.0); assert_eq!(meter.calls, 1);
        account.reserve(&study, &planned, "explicit-retry", cx).unwrap().acquire(&study, &mut meter, cx).unwrap_err();
        assert_eq!(meter.calls, 2); assert_eq!(account.committed_cost_upper(), 4.0);
        assert!(matches!(account.reserve(&study, &planned, "attempt-2", cx), Err(AcquisitionLimitError::AttemptLimit { maximum: 2 })));
        assert!(study.updates().is_empty());
    });
}

#[test]
fn account_and_pending_reading_restart_do_not_reacquire_or_recharge() {
    with_cx(Budget::INFINITE, |cx| {
        let mut study = study(cost(2.0), cx); let planned = request(&study, cx);
        let mut account = AcquisitionAccount::new(cost(4.0), 2).unwrap();
        let source = study.to_checkpoint(cx).unwrap();
        let permit = account.reserve(&study, &planned, "sample-0", cx).unwrap();
        let charged = account.to_checkpoint(cx).unwrap(); // Persist BEFORE I/O.
        let mut meter = Meter::default();
        let acquired = permit.acquire(&study, &mut meter, cx).unwrap();
        let pending = acquired.to_checkpoint(cx).unwrap();
        acquired.apply(&mut study, cx).unwrap();
        let mut restored = MeasurementStudy::from_checkpoint(&source, cx).unwrap();
        let mut restored_account = AcquisitionAccount::from_checkpoint(cost(4.0), 2, &charged, cx).unwrap();
        assert_eq!(restored_account, account);
        assert!(matches!(restored_account.reserve(&restored, &planned, "sample-0", cx), Err(AcquisitionLimitError::DuplicateAttempt { .. })));
        AcquiredMeasurement::from_checkpoint(&restored, &pending, cx).unwrap().apply(&mut restored, cx).unwrap();
        assert_eq!(restored.history_id(), study.history_id());
        assert_eq!(restored_account.committed_cost_upper(), 2.0); assert_eq!(meter.calls, 1);
        assert_eq!(restored_account.to_checkpoint(cx).unwrap(), charged);
    });
}

#[test]
fn cancellation_before_reservation_is_free_but_abandoned_permits_remain_charged() {
    let (mut account, study, planned) = with_cx(Budget::INFINITE, |cx| {
        let study = study(cost(1.0), cx); let planned = request(&study, cx);
        (AcquisitionAccount::new(cost(2.0), 2).unwrap(), study, planned)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(0)] {
        assert!(with_cx(budget, |cx| account.reserve(&study, &planned, "unused", cx)).is_err());
        assert!(account.charges().is_empty()); assert_eq!(account.committed_cost_upper(), 0.0);
    }
    with_cx(Budget::INFINITE, |cx| { drop(account.reserve(&study, &planned, "abandoned", cx).unwrap()); });
    assert_eq!(account.committed_cost_upper(), 1.0);
    let gate = CancelGate::new(); gate.request();
    let clock = VirtualClock::new(); let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let mut meter = Meter::default();
    let permit = with_cx(Budget::INFINITE, |cx| account.reserve(&study, &planned, "cancel-before-io", cx).unwrap());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 942, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        assert!(matches!(permit.acquire(&study, &mut meter, &cx), Err(AcquisitionError::Before(_))));
    });
    assert!(pool.stats().quiescent()); assert_eq!(meter.calls, 0);
    assert_eq!(account.committed_cost_upper(), 2.0);
}

#[test]
fn post_io_cancellation_retains_both_the_charge_and_actual_reading() {
    struct CancellingMeter<'a> { gate: &'a CancelGate, calls: usize }
    impl MeasurementProvider for CancellingMeter<'_> {
        type Error = MeasurementError;
        fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
            self.calls += 1; let outcome = request.reading(value(0.2), "cancel-after-reading")?;
            self.gate.request(); Ok(outcome)
        }
    }
    let (mut study, account, permit) = with_cx(Budget::INFINITE, |cx| {
        let study = study(cost(1.0), cx); let mut account = AcquisitionAccount::new(cost(1.0), 1).unwrap();
        let permit = account.reserve(&study, &request(&study, cx), "sample", cx).unwrap(); (study, account, permit)
    });
    let gate = CancelGate::new(); let clock = VirtualClock::new();
    let mut meter = CancellingMeter { gate: &gate, calls: 0 };
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let acquired = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 942, kernel_id: 20, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        match permit.acquire(&study, &mut meter, &cx) {
            Err(AcquisitionError::After { acquired, .. }) => acquired,
            other => panic!("expected retained acquisition: {other:?}"),
        }
    });
    assert!(pool.stats().quiescent()); assert_eq!(meter.calls, 1);
    with_cx(Budget::INFINITE, |cx| { acquired.apply(&mut study, cx).unwrap(); });
    assert_eq!(account.committed_cost_upper(), 1.0); assert_eq!(account.charges().len(), 1);
    assert_eq!(meter.calls, 1); assert_eq!(study.updates().len(), 1);
}

#[test]
fn schema_identity_staleness_and_zero_attempts_refuse_before_commitment() {
    with_cx(Budget::INFINITE, |cx| {
        let mut study = study(cost(1.0), cx); let planned = request(&study, cx);
        let physical = AcquisitionCost::dimensional(QtyAny::dimensionless(2.0)).unwrap();
        let mut account = AcquisitionAccount::new(physical, 2).unwrap();
        assert!(matches!(account.reserve(&study, &planned, "sample", cx), Err(AcquisitionLimitError::Schema { .. })));
        assert!(account.charges().is_empty());
        let mut disabled = AcquisitionAccount::new(cost(1.0), 0).unwrap();
        assert!(matches!(disabled.reserve(&study, &planned, "sample", cx), Err(AcquisitionLimitError::AttemptLimit { maximum: 0 })));
        let mut account = AcquisitionAccount::new(cost(1.0), 1).unwrap();
        for id in ["".to_owned(), "bad id".to_owned(), "a".repeat(257)] { assert!(account.reserve(&study, &planned, id, cx).is_err()); }
        study.assimilate(&planned, &planned.reading(value(0.2), "meter").unwrap(), cx).unwrap();
        assert!(matches!(account.reserve(&study, &planned, "stale", cx), Err(AcquisitionLimitError::Measurement(MeasurementError::StaleRequest))));
        assert!(account.charges().is_empty());
        assert!(AcquisitionAccount::new(cost(1.0), 4097).is_err());
    });
}

#[test]
fn tiny_positive_costs_cannot_disappear_in_a_large_committed_balance() {
    with_cx(Budget::INFINITE, |cx| {
        let big = study(cost(1e16), cx); let small = study(cost(1.0), cx);
        let mut account = AcquisitionAccount::new(cost(1e16), 2).unwrap();
        drop(account.reserve(&big, &request(&big, cx), "large", cx).unwrap());
        assert!(matches!(account.reserve(&small, &request(&small, cx), "tiny", cx), Err(AcquisitionLimitError::CostLimit { .. })));
        assert_eq!(account.charges().len(), 1);
        let mut account = AcquisitionAccount::new(cost(1e16 + 2.0), 2).unwrap();
        drop(account.reserve(&big, &request(&big, cx), "large", cx).unwrap());
        drop(account.reserve(&small, &request(&small, cx), "tiny", cx).unwrap());
        assert_eq!(account.committed_cost_upper(), 1e16 + 2.0);
    });
}

#[test]
fn durable_allowance_revalidates_configuration_costs_and_complete_history() {
    with_cx(Budget::INFINITE, |cx| {
        let seconds = AcquisitionCost::dimensional(QtyAny::new(2.0, Dims([0, 0, 1, 0, 0, 0]))).unwrap();
        let limit = AcquisitionCost::dimensional(QtyAny::new(4.0, seconds.spec().dims())).unwrap();
        let study = study(seconds, cx); let planned = request(&study, cx);
        let mut account = AcquisitionAccount::new(limit, 2).unwrap();
        drop(account.reserve(&study, &planned, "a", cx).unwrap());
        drop(account.reserve(&study, &planned, "b", cx).unwrap());
        let bytes = account.to_checkpoint(cx).unwrap();
        assert_eq!(account, AcquisitionAccount::from_checkpoint(limit, 2, &bytes, cx).unwrap());
        assert!(AcquisitionAccount::from_checkpoint(cost(4.0), 2, &bytes, cx).is_err());
        assert!(AcquisitionAccount::from_checkpoint(limit, 3, &bytes, cx).is_err());
        for end in [0, 8, bytes.len() - 1] { assert!(AcquisitionAccount::from_checkpoint(limit, 2, &bytes[..end], cx).is_err()); }
        let reseal = |bytes: &mut [u8]| { let end = bytes.len() - 32;
            let mut hash = DomainHasher::new("org.frankensim.sensorforge.acquisition-account.v1");
            hash.update(&bytes[..end]); bytes[end..].copy_from_slice(hash.finalize().as_bytes()); };
        let mut invalid_cost = bytes.clone();
        invalid_cost[76..84].copy_from_slice(&f64::INFINITY.to_le_bytes()); reseal(&mut invalid_cost);
        assert!(AcquisitionAccount::from_checkpoint(limit, 2, &invalid_cost, cx).is_err());
        let mut over_budget = bytes.clone(); over_budget[76..84].copy_from_slice(&3.0_f64.to_le_bytes()); reseal(&mut over_budget);
        assert!(AcquisitionAccount::from_checkpoint(limit, 2, &over_budget, cx).is_err());
        let mut duplicate = bytes.clone(); duplicate[86] = b'a'; reseal(&mut duplicate);
        assert!(matches!(AcquisitionAccount::from_checkpoint(limit, 2, &duplicate, cx), Err(AcquisitionLimitError::DuplicateAttempt { .. })));
        let mut trailing = bytes.clone(); trailing.insert(trailing.len() - 32, 0); reseal(&mut trailing);
        assert!(AcquisitionAccount::from_checkpoint(limit, 2, &trailing, cx).is_err());
        with_cx(Budget::INFINITE.with_poll_quota(0), |limited| { assert!(AcquisitionAccount::from_checkpoint(limit, 2, &bytes, limited).is_err()); });
    });
}
