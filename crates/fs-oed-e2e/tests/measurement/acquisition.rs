use super::*;
use fs_oed_e2e::plan::acquisition::{AcquisitionError, MeasurementProvider};
use std::io;

struct Meter { calls: usize, fail: bool, value: f64 }
impl MeasurementProvider for Meter {
    type Error = io::Error;
    fn read(&mut self, request: &MeasurementRequest, _cx: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
        self.calls += 1;
        if self.fail { return Err(io::Error::other("instrument unavailable")); }
        request.reading(value(self.value), "external-meter").map_err(io::Error::other)
    }
}

#[test]
fn external_read_and_numerical_application_are_separate() {
    context(Budget::INFINITE, |cx| {
        let mut study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx);
        let original = study.clone();
        let mut meter = Meter { calls: 0, fail: false, value: 10.1 };
        let acquired = study.acquire(&selected, &mut meter, cx).unwrap();
        assert_eq!(study, original);
        assert_eq!(meter.calls, 1);
        assert_eq!(acquired.outcome().value(), value(10.1));
        acquired.apply(&mut study, cx).unwrap();
        assert_eq!(meter.calls, 1);
        let committed = study.clone();
        assert!(acquired.apply(&mut study, cx).is_err());
        assert_eq!(study, committed);
    });
}

#[test]
fn stale_requests_and_pre_io_cancellation_never_call_the_provider() {
    let (study, selected) = context(Budget::INFINITE, |cx| {
        let mut study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx);
        let mut meter = Meter { calls: 0, fail: false, value: 10.1 };
        let acquired = study.acquire(&selected, &mut meter, cx).unwrap();
        acquired.apply(&mut study, cx).unwrap();
        assert!(matches!(study.acquire(&selected, &mut meter, cx), Err(AcquisitionError::Before(MeasurementError::StaleRequest))));
        assert_eq!(meter.calls, 1);
        let selected = request(&study, cx);
        (study, selected)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(1)] {
        context(budget, |cx| {
            let mut meter = Meter { calls: 0, fail: false, value: 10.1 };
            assert!(matches!(study.acquire(&selected, &mut meter, cx), Err(AcquisitionError::Before(_))));
            assert_eq!(meter.calls, 0);
        });
    }
}

#[test]
fn provider_error_keeps_accepted_state_and_names_the_unfulfilled_request() {
    context(Budget::INFINITE, |cx| {
        let study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx);
        let original = study.clone();
        let mut meter = Meter { calls: 0, fail: true, value: 0.0 };
        match study.acquire(&selected, &mut meter, cx) {
            Err(AcquisitionError::Provider { request, source }) => {
                assert_eq!(request, selected);
                assert_eq!(source.to_string(), "instrument unavailable");
            }
            other => panic!("expected provider error, got {other:?}"),
        }
        assert_eq!(study, original);
        assert_eq!(meter.calls, 1);
    });
}

#[test]
fn post_io_cancellation_returns_the_real_reading_for_retry_without_resampling() {
    struct CancellingMeter<'a>(&'a CancelGate, usize);
    impl MeasurementProvider for CancellingMeter<'_> {
        type Error = MeasurementError;
        fn read(&mut self, request: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> {
            self.1 += 1;
            let outcome = request.reading(value(10.1), "completed-acquisition")?;
            self.0.request();
            Ok(outcome)
        }
    }
    let (mut study, selected) = context(Budget::INFINITE, |cx| {
        let study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx); (study, selected)
    });
    let original = study.clone();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new(); let clock = VirtualClock::new();
    let mut meter = CancellingMeter(&gate, 0);
    let acquired = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 20, kernel_id: 4, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic).with_time_source(&clock);
        match study.acquire(&selected, &mut meter, &cx) {
            Err(AcquisitionError::After { acquired, .. }) => acquired,
            other => panic!("expected retained reading, got {other:?}"),
        }
    });
    assert!(pool.stats().quiescent());
    assert_eq!(study, original);
    assert_eq!(meter.1, 1);
    context(Budget::INFINITE, |cx| { acquired.apply(&mut study, cx).unwrap(); });
    assert_eq!(meter.1, 1);
    assert_eq!(study.updates()[0].outcome().instrument(), "completed-acquisition");
}

#[test]
fn wrong_provider_request_is_retained_but_never_silently_rebound() {
    struct Wrong(MeasurementOutcome);
    impl MeasurementProvider for Wrong {
        type Error = io::Error;
        fn read(&mut self, _: &MeasurementRequest, _: &Cx<'_>) -> Result<MeasurementOutcome, Self::Error> { Ok(self.0.clone()) }
    }
    context(Budget::INFINITE, |cx| {
        let mut study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let other = MeasurementStudy::new(alternatives(), value(0.001), 2, cx).unwrap();
        let selected = request(&study, cx); let foreign = request(&other, cx);
        let original = study.clone();
        let mut provider = Wrong(foreign.reading(value(10.1), "other-meter").unwrap());
        let acquired = match study.acquire(&selected, &mut provider, cx) {
            Err(AcquisitionError::After { acquired, source: MeasurementError::RequestMismatch }) => acquired,
            other => panic!("expected wrong-request refusal, got {other:?}"),
        };
        assert_eq!(acquired.outcome().request_identity(), foreign.identity());
        assert!(matches!(acquired.apply(&mut study, cx), Err(MeasurementError::RequestMismatch)));
        assert_eq!(study, original);
    });
}
