//! G0/G3/G4: stochastic execution through the actual acquisition/update APIs.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::acquisition::simulation::{GaussianMeasurementSimulator, SimulationTruth};
use fs_oed_e2e::plan::acquisition::simulation::study::SimulatedStudy;
use fs_oed_e2e::plan::measurement::{MeasurementAlternative, MeasurementDecision, MeasurementStop};
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::{Dims, QtyAny};

fn context<T>(budget: Budget, cancelled: bool, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    if cancelled { gate.request(); }
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 20, kernel_id: 4, tile: 0, iteration: 0 }, budget,
        ExecMode::Deterministic).with_time_source(&clock)));
    assert!(pool.stats().quiescent());
    result
}
fn value(x: f64) -> ObjectiveValue { ObjectiveValue::dimensionless(x).unwrap() }
fn menu() -> Vec<MeasurementAlternative> {
    [("A", 10.0), ("B", 10.2)].into_iter().map(|(name, mean)| {
        MeasurementAlternative::new(name, value(mean), QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(0.25), AcquisitionCost::relative(1.0).unwrap()).unwrap()
    }).collect()
}
fn source(seed: u64, a: f64, b: f64, cx: &Cx<'_>) -> GaussianMeasurementSimulator {
    GaussianMeasurementSimulator::new(seed, vec![
        SimulationTruth::new("A", value(a)).unwrap(),
        SimulationTruth::new("B", value(b)).unwrap(),
    ], cx).unwrap()
}
fn session(cap: usize, cx: &Cx<'_>) -> SimulatedStudy {
    SimulatedStudy::new(menu(), value(0.0), cap, source(42, 10.1, 10.2, cx), cx).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-11 * (1.0 + a.abs().max(b.abs())), "{a} != {b}");
}

#[test]
fn segments_and_canonical_menus_replay_the_actual_provider_loop() {
    context(Budget::INFINITE, false, |cx| {
        let mut whole = session(8, cx);
        let mut split = whole.clone();
        let mut manual = whole.study().clone();
        let mut provider = source(42, 10.1, 10.2, cx);
        whole.advance(4, cx).unwrap();
        for _ in 0..4 {
            split.advance(1, cx).unwrap();
            match manual.recommend(cx).unwrap() {
                MeasurementDecision::Stop(_) => break,
                MeasurementDecision::Acquire(request) => {
                    manual.acquire(&request, &mut provider, cx).unwrap().apply(&mut manual, cx).unwrap();
                }
            }
        }
        assert_eq!(whole, split);
        assert_eq!(whole.study(), &manual);
        assert!(!whole.study().updates().is_empty());
        let mut reversed = menu(); reversed.reverse();
        let mut permuted = SimulatedStudy::new(reversed, value(0.0), 8, provider, cx).unwrap();
        permuted.advance(4, cx).unwrap();
        assert_eq!(permuted, whole);
        for update in whole.study().updates() {
            close(update.posterior_variance().value, update.request().predicted_posterior_variance().value);
            assert!(update.outcome().instrument().starts_with("simulated-gaussian-v1-"));
        }
    });
}

#[test]
fn truth_is_not_a_planning_input_but_real_readings_change_the_posterior() {
    context(Budget::INFINITE, false, |cx| {
        let mut ordinary = session(4, cx);
        let mut other = SimulatedStudy::new(menu(), value(0.0), 4,
            source(42, 20.0, -20.0, cx), cx).unwrap();
        assert_eq!(ordinary.study(), other.study());
        assert_eq!(ordinary.study().recommend(cx), other.study().recommend(cx));
        assert_ne!(ordinary.source_identity(), other.source_identity());
        ordinary.advance(1, cx).unwrap(); other.advance(1, cx).unwrap();
        assert_ne!(ordinary.study().history_id(), other.study().history_id());
        let name = ordinary.study().updates()[0].request().alternative();
        assert_ne!(ordinary.study().posterior_mean(name), other.study().posterior_mean(name));
    });
}

#[test]
fn segment_reading_and_threshold_stops_are_not_conflated() {
    context(Budget::INFINITE, false, |cx| {
        let mut live = session(1, cx);
        let before = live.clone();
        let empty = live.advance(0, cx).unwrap();
        assert_eq!(empty.applied(), 0); assert_eq!(empty.stop(), None);
        assert_eq!(live, before);
        let first = live.advance(1, cx).unwrap();
        assert_eq!(first.applied(), 1); assert_eq!(first.stop(), None);
        let stop = live.advance(1, cx).unwrap();
        assert_eq!(stop.applied(), 0); assert_eq!(stop.stop(), Some(MeasurementStop::ReadingLimit));
        assert!(live.study().evpi().value() > 0.0);
        let mut satisfied = SimulatedStudy::new(menu(), value(100.0), 4,
            source(42, 10.1, 10.2, cx), cx).unwrap();
        assert_eq!(satisfied.advance(1, cx).unwrap().stop(), Some(MeasurementStop::ThresholdReached));
        assert!(satisfied.study().updates().is_empty());
    });
}

#[test]
fn cancelled_or_refused_resume_keeps_the_prefix_and_retries_identically() {
    let mut live = context(Budget::INFINITE, false, |cx| {
        let mut live = session(8, cx); live.advance(1, cx).unwrap(); live
    });
    let mut expected = live.clone();
    context(Budget::INFINITE, true, |cx| assert!(live.advance(2, cx).is_err()));
    assert_eq!(live, expected);
    context(Budget::INFINITE.with_cost_quota(1), false, |cx| assert!(live.advance(2, cx).is_err()));
    assert_eq!(live, expected);
    context(Budget::INFINITE, false, |cx| {
        live.advance(2, cx).unwrap(); expected.advance(2, cx).unwrap();
    });
    assert_eq!(live, expected);
}

#[test]
fn complete_same_schema_truth_menu_is_required_before_execution() {
    context(Budget::INFINITE, false, |cx| {
        for names in [vec!["A"], vec!["A", "C"]] {
            let truths = names.into_iter().map(|name| SimulationTruth::new(name, value(1.0)).unwrap()).collect();
            let provider = GaussianMeasurementSimulator::new(42, truths, cx).unwrap();
            assert!(SimulatedStudy::new(menu(), value(0.0), 4, provider, cx).is_err());
        }
        let wrong = ObjectiveValue::dimensional(QtyAny::new(1.0, Dims([1, 0, 0, 0, 0, 0]))).unwrap();
        let truths = ["A", "B"].into_iter().map(|name| SimulationTruth::new(name, wrong).unwrap()).collect();
        let provider = GaussianMeasurementSimulator::new(42, truths, cx).unwrap();
        assert!(SimulatedStudy::new(menu(), value(0.0), 4, provider, cx).is_err());
        let mut live = session(4, cx);
        assert!(live.advance(fs_oed_e2e::MAX_CAMPAIGN_SENSORS + 1, cx).is_err());
        assert!(live.study().updates().is_empty());
    });
}

#[test]
fn realized_loss_only_reads_the_selected_synthetic_diagnostic() {
    context(Budget::INFINITE, false, |cx| {
        let mut live = SimulatedStudy::new(menu(), value(0.0), 4,
            source(42, 12.0, 10.5, cx), cx).unwrap();
        assert_eq!(live.study().best_alternative(), "A");
        let before = live.clone(); close(live.realized_loss().unwrap().value(), 1.5);
        assert_eq!(live, before);
        live.advance(3, cx).unwrap();
        let loss = if live.study().best_alternative() == "A" { 1.5 } else { 0.0 };
        close(live.realized_loss().unwrap().value(), loss);
    });
}
