use fs_evidence::Color;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::measurement::*;
use fs_oed_e2e::{AcquisitionCost, Candidate, ObjectiveValue, run_campaign};
use fs_qty::QtyAny;

fn context<T>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| f(&Cx::new(&gate, arena,
        StreamKey { seed: 20, kernel_id: 4, tile: 0, iteration: 0 }, budget,
        ExecMode::Deterministic).with_time_source(&clock)));
    assert!(pool.stats().quiescent());
    result
}
fn value(x: f64) -> ObjectiveValue { ObjectiveValue::dimensionless(x).unwrap() }
fn alternative(name: &str, mean: f64, variance: f64, noise: f64) -> MeasurementAlternative {
    MeasurementAlternative::new(name, value(mean), QtyAny::dimensionless(variance),
        QtyAny::dimensionless(noise), AcquisitionCost::relative(1.0).unwrap()).unwrap()
}
fn alternatives() -> Vec<MeasurementAlternative> {
    vec![alternative("A", 10.0, 1.0, 0.01), alternative("B", 10.2, 1.0, 0.01)]
}
fn request(study: &MeasurementStudy, cx: &Cx<'_>) -> MeasurementRequest {
    match study.recommend(cx).unwrap() {
        MeasurementDecision::Acquire(request) => request,
        stop => panic!("expected measurement, got {stop:?}"),
    }
}
fn close(a: f64, b: f64) { assert!((a-b).abs() <= 2e-11 * (1.0 + a.abs().max(b.abs())), "{a} != {b}"); }

#[test]
fn production_decisions_match_existing_planner_and_real_kalman_updates() {
    context(Budget::INFINITE, |cx| {
        let legacy = vec![
            Candidate::new("A", value(10.1), value(10.0), QtyAny::dimensionless(1.0), QtyAny::dimensionless(0.01), QtyAny::dimensionless(1.0)).unwrap(),
            Candidate::new("B", value(10.1), value(10.2), QtyAny::dimensionless(1.0), QtyAny::dimensionless(0.01), QtyAny::dimensionless(1.0)).unwrap(),
        ];
        let reference = run_campaign(&legacy, value(0.0), 2, cx).unwrap();
        let mut study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        for name in reference.placements() {
            let selected = request(&study, cx);
            assert_eq!(selected.alternative(), name);
            let acquired = selected.reading(value(10.1), "actual-meter").unwrap();
            let result = study.assimilate(&selected, &acquired, cx).unwrap();
            close(result.posterior_variance().value, selected.predicted_posterior_variance().value);
            assert_eq!(result.outcome().value(), value(10.1));
        }
        close(study.evpi().value(), reference.final_evpi().value());
        assert_eq!(study.updates().len(), reference.placements().len());
        assert!(matches!(study.evpi_color(), Color::Estimated { dispersion, .. } if dispersion.is_infinite()));
    });
}

#[test]
fn noise_changes_the_actual_decision_and_readings_change_the_choice() {
    context(Budget::INFINITE, |cx| {
        let mut study = MeasurementStudy::new(vec![alternative("precise", 0.0, 1.0, 0.01), alternative("noisy", 0.0, 1.0, 100.0)], value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx);
        assert_eq!(selected.alternative(), "precise");
        let original = study.clone();
        study.assimilate(&selected, &selected.reading(value(-3.0), "meter").unwrap(), cx).unwrap();
        assert_eq!(study.best_alternative(), "precise");
        let mut other = original.clone();
        other.assimilate(&selected, &selected.reading(value(3.0), "meter").unwrap(), cx).unwrap();
        assert_eq!(other.best_alternative(), "noisy");
        close(study.posterior_variance("precise").unwrap().value, other.posterior_variance("precise").unwrap().value);
        assert_ne!(study.history_id(), other.history_id());
        assert_eq!(request(&original, cx), selected); // planning read no future outcome
    });
}

#[test]
fn menu_permutations_and_exact_resume_preserve_requests_and_state() {
    context(Budget::INFINITE, |cx| {
        let menu = alternatives();
        let mut reverse = menu.clone(); reverse.reverse();
        let mut a = MeasurementStudy::new(menu, value(0.0), 3, cx).unwrap();
        let mut b = MeasurementStudy::new(reverse, value(0.0), 3, cx).unwrap();
        assert_eq!(a, b);
        for _ in 0..2 {
            let ra = request(&a, cx); let rb = request(&b, cx); assert_eq!(ra, rb);
            let outcome = ra.reading(value(10.1), "meter").unwrap();
            a.assimilate(&ra, &outcome, cx).unwrap();
            b = b.clone(); // cloning may discard reserved spare Vec capacity
            b.assimilate(&rb, &outcome, cx).unwrap();
            assert_eq!(a, b);
        }
    });
}

#[test]
fn stale_foreign_and_wrong_request_readings_never_mutate_the_study() {
    context(Budget::INFINITE, |cx| {
        let mut a = MeasurementStudy::new(alternatives(), value(0.0), 3, cx).unwrap();
        let b = MeasurementStudy::new(alternatives(), value(0.001), 3, cx).unwrap();
        let ra = request(&a, cx); let rb = request(&b, cx);
        let original = a.clone();
        let wrong = rb.reading(value(10.1), "meter").unwrap();
        assert!(matches!(a.assimilate(&ra, &wrong, cx), Err(MeasurementError::RequestMismatch)));
        assert!(matches!(a.assimilate(&rb, &wrong, cx), Err(MeasurementError::StaleRequest)));
        assert_eq!(a, original);
        let reading = ra.reading(value(10.1), "meter").unwrap();
        a.assimilate(&ra, &reading, cx).unwrap();
        let accepted = a.clone();
        assert!(matches!(a.assimilate(&ra, &reading, cx), Err(MeasurementError::StaleRequest)));
        assert_eq!(a, accepted);
    });
}

#[test]
fn cancellation_and_budget_exhaustion_preserve_readings_for_retry() {
    let (study, selected, reading) = context(Budget::INFINITE, |cx| {
        let study = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&study, cx);
        let reading = selected.reading(value(10.1), "meter").unwrap();
        (study, selected, reading)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(2), Budget::INFINITE.with_cost_quota(1)] {
        context(budget, |cx| {
            let mut candidate = study.clone();
            assert!(candidate.assimilate(&selected, &reading, cx).is_err());
            assert_eq!(candidate, study);
        });
    }
    context(Budget::INFINITE, |cx| {
        let mut retry = study;
        retry.assimilate(&selected, &reading, cx).unwrap();
        assert_eq!(retry.updates().len(), 1);
    });
}

#[test]
fn whole_menu_schema_and_noise_admission_precede_any_acquisition() {
    context(Budget::INFINITE, |cx| {
        for noise in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(MeasurementAlternative::new("x", value(0.0), QtyAny::dimensionless(1.0), QtyAny::dimensionless(noise), AcquisitionCost::relative(1.0).unwrap()).is_err());
        }
        assert!(MeasurementStudy::new(vec![alternative("x", 0.0, 1.0, 1.0); 2], value(0.0), 1, cx).is_err());
        let dimensional_cost = MeasurementAlternative::new("y", value(0.0), QtyAny::dimensionless(1.0), QtyAny::dimensionless(1.0), AcquisitionCost::dimensional(QtyAny::dimensionless(1.0)).unwrap()).unwrap();
        assert!(MeasurementStudy::new(vec![alternative("x", 0.0, 1.0, 1.0), dimensional_cost], value(0.0), 1, cx).is_err());
        let study = MeasurementStudy::new(alternatives(), value(0.0), 1, cx).unwrap();
        let selected = request(&study, cx);
        assert!(selected.reading(ObjectiveValue::dimensional(fs_qty::parse::parse_qty("2 m").unwrap()).unwrap(), "meter").is_err());
        assert!(selected.reading(value(1.0), "").is_err());
    });
}

#[test]
fn zero_reading_cap_is_not_a_robustness_claim() {
    context(Budget::INFINITE, |cx| {
        let ambiguous = MeasurementStudy::new(alternatives(), value(0.0), 0, cx).unwrap();
        assert!(ambiguous.evpi().value() > 0.0);
        assert_eq!(ambiguous.recommend(cx).unwrap(), MeasurementDecision::Stop(MeasurementStop::ReadingLimit));
        let exact = MeasurementStudy::new(vec![alternative("a", 1.0, 0.0, 1.0), alternative("b", 2.0, 0.0, 1.0)], value(0.0), 0, cx).unwrap();
        assert_eq!(exact.recommend(cx).unwrap(), MeasurementDecision::Stop(MeasurementStop::ThresholdReached));
    });
}

#[test]
fn actual_instrument_and_noise_declarations_are_bound_to_history() {
    context(Budget::INFINITE, |cx| {
        let original = MeasurementStudy::new(alternatives(), value(0.0), 2, cx).unwrap();
        let selected = request(&original, cx);
        let mut a = original.clone(); let mut b = original;
        a.assimilate(&selected, &selected.reading(value(10.1), "instrument-one").unwrap(), cx).unwrap();
        b.assimilate(&selected, &selected.reading(value(10.1), "instrument-two").unwrap(), cx).unwrap();
        assert_eq!(a.posterior_mean(selected.alternative()), b.posterior_mean(selected.alternative()));
        assert_ne!(a.history_id(), b.history_id());
        assert_ne!(a.updates()[0].assimilation_color(), b.updates()[0].assimilation_color());
    });
}

#[path = "measurement/acquisition.rs"]
mod acquisition;
