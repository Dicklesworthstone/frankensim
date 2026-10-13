use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_oed_e2e::plan::acquisition::simulation::{GaussianMeasurementSimulator, SimulationTruth};
use fs_oed_e2e::plan::measurement::*;
use fs_oed_e2e::{AcquisitionCost, ObjectiveValue};
use fs_qty::{Dims, QtyAny};
use fs_qty::semantic::{QuantityKind, SemanticQty, SemanticType, ValueForm};

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
fn study(cx: &Cx<'_>) -> MeasurementStudy {
    let menu = [("A", 10.0), ("B", 10.2)].into_iter().map(|(name, mean)| {
        MeasurementAlternative::new(name, value(mean), QtyAny::dimensionless(1.0),
            QtyAny::dimensionless(0.25), AcquisitionCost::relative(1.0).unwrap()).unwrap()
    }).collect();
    MeasurementStudy::new(menu, value(0.0), 8, cx).unwrap()
}
fn request(study: &MeasurementStudy, cx: &Cx<'_>) -> MeasurementRequest {
    match study.recommend(cx).unwrap() {
        MeasurementDecision::Acquire(r) => r, other => panic!("expected request: {other:?}"),
    }
}
fn simulator(seed: u64, latent: f64, cx: &Cx<'_>) -> GaussianMeasurementSimulator {
    GaussianMeasurementSimulator::new(seed, ["A", "B"].into_iter()
        .map(|name| SimulationTruth::new(name, value(latent)).unwrap()).collect(), cx).unwrap()
}
fn close(a: f64, b: f64) { assert!((a-b).abs() < 2e-11 * (1.0 + a.abs().max(b.abs())), "{a} != {b}"); }

#[test]
fn replay_uses_the_native_counter_based_normal_and_not_provider_call_order() {
    context(Budget::INFINITE, |cx| {
        let mut study = study(cx);
        let mut source = simulator(123, 10.1, cx);
        let first = request(&study, cx);
        let a = source.sample(&first, cx).unwrap();
        assert_eq!(a.standard_normal(), a.stream_key().stream().next_normal());
        let acquired = study.acquire(&first, &mut source, cx).unwrap();
        assert_eq!(acquired.outcome(), a.outcome());
        acquired.apply(&mut study, cx).unwrap();
        let second = request(&study, cx);
        let b = source.sample(&second, cx).unwrap();
        assert_ne!(a.stream_key(), b.stream_key());
        assert_eq!(a, source.sample(&first, cx).unwrap());
        assert_eq!(b, source.clone().sample(&second, cx).unwrap());
        assert_ne!(a, simulator(124, 10.1, cx).sample(&first, cx).unwrap());
        assert!(a.outcome().instrument().starts_with("simulated-gaussian-v1-"));
    });
}

#[test]
fn noise_moments_match_the_declared_instrument_variance_not_exact_truth() {
    context(Budget::INFINITE, |cx| {
        let r = request(&study(cx), cx);
        let (mut mean, mut second) = (0.0, 0.0);
        for seed in 0..2048 {
            let sample = simulator(seed, 10.1, cx).sample(&r, cx).unwrap();
            let z = (sample.outcome().value().value() - 10.1) / r.noise_variance().value.sqrt();
            close(z, sample.standard_normal());
            mean += z; second += z * z;
        }
        mean /= 2048.0; second /= 2048.0;
        assert!(mean.abs() < 0.08, "mean={mean}");
        assert!((second - 1.0).abs() < 0.15, "second moment={second}");
    });
}

#[test]
fn hidden_truth_changes_outcomes_but_never_the_request_or_noise_address() {
    context(Budget::INFINITE, |cx| {
        let study = study(cx);
        let original = study.clone();
        let r = request(&study, cx);
        let a = simulator(8, 9.0, cx).sample(&r, cx).unwrap();
        let b = simulator(8, 11.0, cx).sample(&r, cx).unwrap();
        assert_eq!(a.stream_key(), b.stream_key());
        assert_eq!(a.standard_normal(), b.standard_normal());
        assert_ne!(a.simulator_identity(), b.simulator_identity());
        assert_ne!(a.outcome(), b.outcome());
        assert_eq!(study, original);
        assert_eq!(r, request(&study, cx));
        let mut left = study.clone(); let mut right = study;
        left.assimilate(&r, a.outcome(), cx).unwrap();
        right.assimilate(&r, b.outcome(), cx).unwrap();
        assert_ne!(left.history_id(), right.history_id());
        close(left.posterior_variance(r.alternative()).unwrap().value,
            right.posterior_variance(r.alternative()).unwrap().value);
    });
}

#[test]
fn seeded_closed_loop_replays_after_cloning_without_reusing_noise() {
    context(Budget::INFINITE, |cx| {
        let mut first = study(cx); let mut replay = first.clone();
        let mut source = simulator(32, 10.1, cx);
        for _ in 0..4 {
            let r = request(&first, cx);
            let rr = request(&replay, cx);
            assert_eq!(r, rr);
            let prior = first.posterior_mean(r.alternative()).unwrap().value();
            let variance = first.posterior_variance(r.alternative()).unwrap().value;
            let acquired = first.acquire(&r, &mut source, cx).unwrap();
            let expected = prior + variance / (variance + r.noise_variance().value)
                * (acquired.outcome().value().value() - prior);
            acquired.apply(&mut first, cx).unwrap();
            let mut restored = source.clone();
            replay.acquire(&rr, &mut restored, cx).unwrap().apply(&mut replay, cx).unwrap();
            assert_eq!(first, replay);
            close(first.posterior_mean(r.alternative()).unwrap().value(), expected);
        }
    });
}

#[test]
fn input_schema_identity_and_numerical_range_are_checked_without_resampling() {
    context(Budget::INFINITE, |cx| {
        let r = request(&study(cx), cx);
        let row = SimulationTruth::new("A", value(1.0)).unwrap();
        assert!(GaussianMeasurementSimulator::new(1, vec![], cx).is_err());
        assert!(GaussianMeasurementSimulator::new(1, vec![row.clone(), row], cx).is_err());
        let unknown = GaussianMeasurementSimulator::new(1,
            vec![SimulationTruth::new("absent", value(1.0)).unwrap()], cx).unwrap();
        assert!(unknown.sample(&r, cx).is_err());
        let wrong = ObjectiveValue::dimensional(QtyAny::new(1.0, Dims([1, 0, 0, 0, 0, 0]))).unwrap();
        let source = GaussianMeasurementSimulator::new(1,
            vec![SimulationTruth::new(r.alternative(), wrong).unwrap()], cx).unwrap();
        assert!(source.sample(&r, cx).is_err());
        assert!(simulator(1, f64::MAX, cx).sample(&r, cx).is_err());
        let a = simulator(5, 1.0, cx);
        let b = GaussianMeasurementSimulator::new(5,
            vec![SimulationTruth::new("B", value(1.0)).unwrap(), SimulationTruth::new("A", value(1.0)).unwrap()], cx).unwrap();
        assert_eq!(a, b);
    });
}

#[test]
fn semantic_absolute_temperature_gets_delta_noise_and_a_checked_absolute_reading() {
    context(Budget::INFINITE, |cx| {
        let kind = SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static);
        let temp = ObjectiveValue::semantic(SemanticQty::new(QtyAny::new(295.0, Dims([0, 0, 0, 1, 0, 0])), kind).unwrap());
        let dims = temp.spec().variance_dims().unwrap();
        let menu = ["A", "B"].into_iter().map(|name| MeasurementAlternative::new(name, temp,
            QtyAny::new(1.0, dims), QtyAny::new(0.25, dims), AcquisitionCost::relative(1.0).unwrap()).unwrap()).collect();
        let study = MeasurementStudy::new(menu, temp.spec().decision_value(0.0).unwrap(), 2, cx).unwrap();
        let r = request(&study, cx);
        let source = GaussianMeasurementSimulator::new(5,
            ["A", "B"].into_iter().map(|name| SimulationTruth::new(name, temp).unwrap()).collect(), cx).unwrap();
        let result = source.sample(&r, cx).unwrap();
        assert_eq!(result.outcome().value().spec(), temp.spec());
        assert_eq!(result.perturbation().spec(), temp.spec().decision_spec());
        close(result.outcome().value().value(), 295.0 + result.perturbation().value());
    });
}

#[test]
fn exhausted_simulation_budgets_preserve_requests_for_identical_retry() {
    let (source, r, expected) = context(Budget::INFINITE, |cx| {
        let r = request(&study(cx), cx); let source = simulator(19, 10.1, cx);
        let expected = source.sample(&r, cx).unwrap(); (source, r, expected)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_cost_quota(1)] {
        context(budget, |cx| assert!(source.sample(&r, cx).is_err()));
    }
    context(Budget::INFINITE, |cx| assert_eq!(source.sample(&r, cx).unwrap(), expected));
}
