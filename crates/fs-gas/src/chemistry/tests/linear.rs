use super::*;
use crate::chemistry::isobaric::IsobaricReactor;
use crate::chemistry::linear::{FirstOrderLimits, PreparedFirstOrder};
use fs_thermochem::kinetics::{PressureDependenceV1, ThirdBodyEfficienciesV1};

fn exact_controls(step: f64) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::ExactFirstOrder {
            limits: FirstOrderLimits::default(),
        },
        maximum_step_s: step,
        initial_step_s: step,
        maximum_steps: 64,
        ..controls(1e-6)
    }
}

fn half_reaction(
    name: &str,
    parent: &str,
    products: &[(&str, u32)],
    half: f64,
) -> ElementaryReactionV1 {
    ElementaryReactionV1::from_half_life(
        ReactionId::new(name).unwrap(),
        id(parent),
        products
            .iter()
            .map(|(name, order)| (id(name), *order))
            .collect(),
        Time::new(half),
        Temperature::new(500.0),
    )
    .unwrap()
}

fn half_isomer(heat: f64, entropy: f64) -> IdealGasKineticsV1 {
    mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -heat, entropy),
        ],
        vec![1, 1],
        vec![half_reaction("A-B", "A", &[("B", 1)], 1.0)],
    )
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_half_life_uses_one_exponential_step_and_original_isochoric_calorics() {
    let m = half_isomer(2000.0, 5.0);
    let mut control = exact_controls(1.0);
    control.maximum_steps = 1;
    let result = run(&m, 500.0, 1.0, control);
    relative(result.amounts()[0].value(), 0.5, 2e-13);
    relative(result.amounts()[1].value(), 0.5, 2e-13);
    relative(
        result.temperature().value(),
        500.0 + 1000.0 / (3.0 * R),
        2e-12,
    );
    assert_eq!(result.balance().accepted_steps, 1);
    assert_eq!(result.balance().nonlinear_iterations, 0);
    assert_eq!(result.balance().rejected_steps, 0);
    assert_eq!(result.elapsed_time(), Time::new(1.0));
    assert!(result.balance().entropy_production_j_k > 0.0);
    assert!(result.balance().energy_residual_j.abs() <= control.energy_tolerance_j);
    assert!(result.balance().maximum_relative_element_residual < 1e-13);
    assert_eq!(result.mechanism(), &m);
    assert_eq!(result.sources().len(), 2);
}

#[test]
fn g1_repeated_rate_ingrowth_retains_independent_channel_events_and_heat() {
    let m = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -2000.0, 5.0),
            source("C", 0.01, 4.0, -4000.0, 10.0),
        ],
        vec![1, 1, 1],
        vec![
            half_reaction("A-B", "A", &[("B", 1)], 1.0),
            half_reaction("B-C", "B", &[("C", 1)], 1.0),
        ],
    );
    let expected_a = 0.25;
    let expected_b = 0.5 * core::f64::consts::LN_2;
    let expected_c = 1.0 - expected_a - expected_b;
    let result = run(&m, 500.0, 2.0, exact_controls(2.0));
    for (amount, expected) in result
        .amounts()
        .iter()
        .zip([expected_a, expected_b, expected_c])
    {
        relative(amount.value(), expected, 3e-13);
    }
    relative(
        result.temperature().value(),
        500.0 + (2000.0 * expected_b + 4000.0 * expected_c) / (3.0 * R),
        2e-12,
    );
    let mut prepared =
        PreparedFirstOrder::new(&m, FirstOrderLimits::default(), &mut || false).unwrap();
    let mut events = [0.0; 2];
    let _ = prepared
        .advance(&m, &[1.0, 0.0, 0.0], 2.0, &mut events, &mut || false)
        .unwrap();
    relative(events[0], 0.75, 3e-13);
    relative(events[1], expected_c, 3e-13);
}

#[test]
fn g1_competing_and_cyclic_events_are_not_net_population_differences() {
    let m = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 5.0),
            source("C", 0.01, 4.0, 0.0, 5.0),
        ],
        vec![1, 1, 1],
        vec![
            reaction("A-B", &[("A", 1)], &[("B", 1)], 1000.0).irreversible(),
            reaction("A-C", &[("A", 1)], &[("C", 1)], 1e-20).irreversible(),
        ],
    );
    let mut prepared =
        PreparedFirstOrder::new(&m, FirstOrderLimits::default(), &mut || false).unwrap();
    let mut events = [0.0; 2];
    let answer = prepared
        .advance(&m, &[1.0, 0.0, 0.0], 1e-4, &mut events, &mut || false)
        .unwrap();
    let loss = -(-0.1_f64).exp_m1();
    relative(events[0], loss, 3e-13);
    relative(events[1], loss * 1e-23, 3e-13);
    relative(answer[2], events[1], 3e-13);

    let cycle = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 0.0),
        ],
        vec![1, 1],
        vec![
            reaction("A-B", &[("A", 1)], &[("B", 1)], 1.0).irreversible(),
            reaction("B-A", &[("B", 1)], &[("A", 1)], 1.0).irreversible(),
        ],
    );
    let mut prepared =
        PreparedFirstOrder::new(&cycle, FirstOrderLimits::default(), &mut || false).unwrap();
    let answer = prepared
        .advance(&cycle, &[1.0, 0.0], 2.0, &mut events, &mut || false)
        .unwrap();
    let net = 0.5 * (1.0 - (-4.0_f64).exp());
    relative(answer[1], net, 5e-13);
    relative(events[0], 1.0 + 0.5 * net, 5e-13);
    relative(events[1], 1.0 - 0.5 * net, 5e-13);
    assert!(events[0] > 1.0);
    assert!(events[1] > 0.5);
}

#[test]
fn g3_direct_populations_preserve_trace_survivors_and_tiny_loss_events() {
    let m = half_isomer(0.0, 200.0);
    let elapsed = 100.0 / core::f64::consts::LN_2;
    let result = run(&m, 500.0, elapsed, exact_controls(elapsed));
    relative(result.amounts()[0].value(), (-100.0_f64).exp(), 2e-12);
    assert!(result.amounts()[0].value() > 0.0);
    assert_eq!(result.amounts()[1].value(), 1.0);

    let late_chain = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -2000.0, 100.0),
            source("C", 0.01, 4.0, -4000.0, 200.0),
        ],
        vec![1, 1, 1],
        vec![
            half_reaction("A-B", "A", &[("B", 1)], 1.0),
            half_reaction("B-C", "B", &[("C", 1)], 1.0),
        ],
    );
    let late = run(&late_chain, 500.0, 100.0, exact_controls(100.0));
    let expected_a = 2.0_f64.powi(-100);
    let expected_b = 100.0 * core::f64::consts::LN_2 * expected_a;
    relative(late.amounts()[0].value(), expected_a, 2e-12);
    relative(late.amounts()[1].value(), expected_b, 2e-12);
    assert!(late.amounts()[1].value() > 0.0);
    assert!(late.balance().energy_residual_j.abs() < 1e-8);

    let slow = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, 0.0, 5.0),
        ],
        vec![1, 1],
        vec![reaction("A-B", &[("A", 1)], &[("B", 1)], 1e-30).irreversible()],
    );
    let mut prepared =
        PreparedFirstOrder::new(&slow, FirstOrderLimits::default(), &mut || false).unwrap();
    let mut events = [0.0];
    let answer = prepared
        .advance(&slow, &[1e100, 0.0], 1e-20, &mut events, &mut || false)
        .unwrap();
    assert_eq!(answer[0], 1e100);
    relative(answer[1], 1e50, 2e-13);
    relative(events[0], 1e50, 2e-13);
}

#[test]
fn g3_restart_and_isobaric_closure_reuse_original_energy_and_volume() {
    let m = half_isomer(2000.0, 5.0);
    let control = exact_controls(3.0);
    let whole = run(&m, 500.0, 3.0, control);
    let first = run(&m, 500.0, 1.0, control);
    let mut reactor = IsochoricReactor::new(m.clone(), control, || false).unwrap();
    let restarted = reactor
        .advance_adiabatic_energy(
            first.amounts(),
            first.internal_energy(),
            Volume::new(0.01),
            first.temperature(),
            Time::new(2.0),
            || false,
        )
        .unwrap();
    for (a, b) in restarted.amounts().iter().zip(whole.amounts()) {
        relative(a.value(), b.value(), 2e-12);
    }
    relative(
        restarted.temperature().value(),
        whole.temperature().value(),
        2e-12,
    );

    let mut reactor = IsobaricReactor::new(m, exact_controls(1.0), || false).unwrap();
    let result = reactor
        .advance_adiabatic(
            &[Amount::new(1.0), Amount::new(0.0)],
            Pressure::new(100000.0),
            Temperature::new(500.0),
            Time::new(1.0),
            || false,
        )
        .unwrap();
    let expected_temperature = 500.0 + 1000.0 / (4.0 * R);
    relative(result.amounts()[0].value(), 0.5, 3e-13);
    relative(result.temperature().value(), expected_temperature, 2e-12);
    relative(
        result.volume().value(),
        R * expected_temperature / 100000.0,
        2e-12,
    );
    relative(
        result.boundary_work().value(),
        R * (expected_temperature - 500.0),
        3e-11,
    );
    assert!(result.balance().energy_residual_j.abs() <= control.energy_tolerance_j);
    assert_eq!(result.balance().nonlinear_iterations, 0);
}

#[test]
fn g1_isobaric_dissociation_preserves_multiplicity_and_actual_diluting_volume() {
    let m = mechanism(
        vec![
            source("A", 0.02, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -2000.0, 5.0),
        ],
        vec![2, 1],
        vec![half_reaction("A-2B", "A", &[("B", 2)], 1.0)],
    );
    let mut reactor = IsobaricReactor::new(m, exact_controls(1.0), || false).unwrap();
    let result = reactor
        .advance_adiabatic(
            &[Amount::new(1.0), Amount::new(0.0)],
            Pressure::new(100000.0),
            Temperature::new(500.0),
            Time::new(1.0),
            || false,
        )
        .unwrap();
    let temperature = (4.0 * R * 500.0 + 2000.0) / (6.0 * R);
    relative(result.amounts()[0].value(), 0.5, 3e-13);
    relative(result.amounts()[1].value(), 1.0, 3e-13);
    relative(result.temperature().value(), temperature, 2e-12);
    relative(
        result.volume().value(),
        1.5 * R * temperature / 100000.0,
        2e-12,
    );
    relative(
        result.boundary_work().value(),
        R * (1.5 * temperature - 500.0),
        3e-11,
    );
    assert!(result.balance().mass_residual_kg.abs() < 1e-15);
    assert!(result.balance().energy_residual_j.abs() < 1e-8);
}

#[test]
fn g0_ineligible_mechanisms_heat_exchange_and_numerical_budgets_refuse() {
    let control = exact_controls(1.0);
    assert!(matches!(
        IsochoricReactor::new(isomers(0.0, 1.0), control, || false),
        Err(IsochoricKineticsError::Invalid(_))
    ));
    let with_law = |rate: ArrheniusRateV1, pressure| {
        mechanism(
            vec![
                source("A", 0.01, 4.0, 0.0, 0.0),
                source("B", 0.01, 4.0, 0.0, 5.0),
            ],
            vec![1, 1],
            vec![
                ElementaryReactionV1::new(
                    ReactionId::new("A-B").unwrap(),
                    vec![(id("A"), 1)],
                    vec![(id("B"), 1)],
                    rate,
                )
                .unwrap()
                .irreversible()
                .with_pressure_dependence(pressure),
            ],
        )
    };
    let law = |exponent, energy| {
        ArrheniusRateV1::new(
            Temperature::new(500.0),
            ReactionProgressRateQuantityV1::new(1.0),
            exponent,
            MolarEnergyQuantityV1::new(energy),
        )
        .unwrap()
    };
    for m in [
        with_law(law(1.0, 0.0), PressureDependenceV1::Independent),
        with_law(law(0.0, 1000.0), PressureDependenceV1::Independent),
        with_law(
            law(0.0, 0.0),
            PressureDependenceV1::ThirdBody(ThirdBodyEfficienciesV1::default()),
        ),
    ] {
        assert!(matches!(
            IsochoricReactor::new(m, control, || false),
            Err(IsochoricKineticsError::Invalid(_))
        ));
    }
    let mut reactor = IsochoricReactor::new(half_isomer(2000.0, 5.0), control, || false).unwrap();
    let feed = [Amount::new(1.0), Amount::new(0.0)];
    assert!(matches!(
        reactor.advance_with_heat_bath(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(1.0),
            HeatBath::new(ThermalConductance::new(1.0), Temperature::new(300.0)).unwrap(),
            || false,
        ),
        Err(IsochoricKineticsError::Invalid(_))
    ));
    let mut limited = control;
    limited.method = IsochoricKineticsMethod::ExactFirstOrder {
        limits: FirstOrderLimits {
            maximum_series_terms: 1,
            ..FirstOrderLimits::default()
        },
    };
    let mut reactor = IsochoricReactor::new(half_isomer(2000.0, 5.0), limited, || false).unwrap();
    assert!(matches!(
        reactor.advance_adiabatic(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(1.0),
            || false,
        ),
        Err(IsochoricKineticsError::FirstOrder(
            fs_la::decay::DecayError::SeriesNotConverged { .. }
        ))
    ));
    assert_eq!(feed, [Amount::new(1.0), Amount::new(0.0)]);
}

#[test]
fn g0_large_cyclic_counters_must_resolve_the_requested_caloric_energy() {
    let m = mechanism(
        vec![
            source("A", 0.01, 4.0, 0.0, 0.0),
            source("B", 0.01, 4.0, -1.0, 0.0),
        ],
        vec![1, 1],
        vec![
            reaction("A-B", &[("A", 1)], &[("B", 1)], 1.0).irreversible(),
            reaction("B-A", &[("B", 1)], &[("A", 1)], 1.0).irreversible(),
        ],
    );
    let mut reactor = IsochoricReactor::new(m, exact_controls(2.0), || false).unwrap();
    let small = reactor
        .advance_adiabatic(
            &[Amount::new(1.0), Amount::new(0.0)],
            Volume::new(1.0),
            Temperature::new(500.0),
            Time::new(2.0),
            || false,
        )
        .unwrap();
    assert!(small.balance().first_order_counter_roundoff_j > 0.0);
    assert!(
        small.balance().energy_residual_j.abs() + small.balance().first_order_counter_roundoff_j
            < 1e-8
    );
    // Same concentrations and generator, but macroscopic event counts near
    // 1e16 mol cannot resolve a 1e-8 J energy account at ±1 J/mol per event.
    let feed = [Amount::new(1e16), Amount::new(0.0)];
    assert!(matches!(
        reactor.advance_adiabatic(
            &feed,
            Volume::new(1e16),
            Temperature::new(500.0),
            Time::new(2.0),
            || false,
        ),
        Err(IsochoricKineticsError::StepFailure {
            reason: "first-order counter energy precision exhausts the energy allowance",
            ..
        })
    ));
    assert_eq!(feed, [Amount::new(1e16), Amount::new(0.0)]);
}

#[test]
fn g4_cancellation_and_refusal_leave_reusable_reactor_and_untouched_feed() {
    let m = half_isomer(2000.0, 5.0);
    let control = exact_controls(1.0);
    let feed = [Amount::new(1.0), Amount::new(0.0)];
    let mut reactor = IsochoricReactor::new(m, control, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        reactor.advance_adiabatic(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(1.0),
            || {
                polls += 1;
                polls == 40
            },
        ),
        Err(IsochoricKineticsError::Cancelled)
    ));
    let expected = reactor
        .advance_adiabatic(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(1.0),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    let replay = reactor
        .advance_adiabatic(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(1.0),
            || false,
        )
        .unwrap()
        .to_owned(|| false)
        .unwrap();
    assert_eq!(expected, replay);
    assert_eq!(feed, [Amount::new(1.0), Amount::new(0.0)]);
    let zero = reactor
        .advance_adiabatic(
            &feed,
            Volume::new(0.01),
            Temperature::new(500.0),
            Time::new(0.0),
            || false,
        )
        .unwrap();
    assert_eq!(zero.amounts(), &feed);
    assert_eq!(zero.balance().accepted_steps, 0);
}
