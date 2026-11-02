use super::*;
use fs_fft::harmonic_window::{
    HarmonicProductInput, HarmonicProductLaw, HarmonicProductSpectrumLimits,
};

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn limits() -> LumpedInterferenceLimits {
    LumpedInterferenceLimits {
        maximum_nodes: MAX_LUMPED_NODES,
        maximum_spectral_terms: MAX_LUMPED_INTERFERENCE_TERMS,
        maximum_workspace_bytes: 2 * 1024 * 1024,
    }
}

fn node(name: &str, capacity: f64, conductance: f64) -> LumpedNode {
    LumpedNode::new(name, capacity, conductance, 0.001, 100.0, 1.0).unwrap()
}

fn spectrum(products: Vec<HarmonicProductInput>) -> PreparedHarmonicProductSpectrum {
    let count = products.len();
    PreparedHarmonicProductSpectrum::new(
        products,
        HarmonicProductSpectrumLimits {
            maximum_products: count,
            maximum_generated_terms: 2 * count,
            maximum_workspace_bytes: 64 * 1024,
        },
    )
    .unwrap()
}

fn constant(value: f64) -> PreparedHarmonicProductSpectrum {
    spectrum(vec![HarmonicProductInput {
        left_frequency_hz: 0.0,
        right_frequency_hz: 0.0,
        left_amplitude: [Dd::from_f64(value), Dd::ZERO],
        right_amplitude: [Dd::ONE, Dd::ZERO],
        law: HarmonicProductLaw::Real,
    }])
}

fn square(
    carriers: &[(f64, [f64; 2])],
    law: HarmonicProductLaw,
) -> PreparedHarmonicProductSpectrum {
    let mut products = Vec::new();
    for &(left_frequency_hz, left) in carriers {
        for &(right_frequency_hz, right) in carriers {
            products.push(HarmonicProductInput {
                left_frequency_hz,
                right_frequency_hz,
                left_amplitude: left.map(Dd::from_f64),
                right_amplitude: right.map(Dd::from_f64),
                law,
            });
        }
    }
    spectrum(products)
}

fn binding<'a>(
    name: &'a str,
    spectrum: &'a PreparedHarmonicProductSpectrum,
    scale: f64,
) -> LumpedPowerSpectrumBinding<'a> {
    LumpedPowerSpectrumBinding {
        node_name: name,
        spectrum,
        input_unit: "W/m^3",
        watts_per_spectrum_unit: scale,
        source_clock: "laboratory trigger",
    }
}

fn prepared<'a>(
    network: &'a LumpedNetwork,
    bindings: &[LumpedPowerSpectrumBinding<'a>],
) -> PreparedLumpedInterference<'a> {
    PreparedLumpedInterference::new(network, BiotGate::corpus_default(), bindings, limits())
        .unwrap()
}

// An independent sampled ODE reference, used only to check the analytic solver.
// The state carries temperature, integrated source power and ambient exchange.
fn rk4(
    power: impl Fn(f64) -> f64,
    start: f64,
    duration: f64,
    initial: f64,
    ambient: f64,
    capacity: f64,
    conductance: f64,
) -> [f64; 3] {
    let steps = 8_000;
    let h = duration / steps as f64;
    let rhs = |t: f64, state: [f64; 3]| {
        let input = power(t);
        let exchange = conductance * (state[0] - ambient);
        [(input - exchange) / capacity, input, exchange]
    };
    let add =
        |a: [f64; 3], b: [f64; 3], factor: f64| core::array::from_fn(|i| a[i] + factor * b[i]);
    let mut state = [initial, 0.0, 0.0];
    for step in 0..steps {
        let t = start + step as f64 * h;
        let k1 = rhs(t, state);
        let k2 = rhs(t + 0.5 * h, add(state, k1, 0.5 * h));
        let k3 = rhs(t + 0.5 * h, add(state, k2, 0.5 * h));
        let k4 = rhs(t + h, add(state, k3, h));
        for i in 0..3 {
            state[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    state
}

#[test]
fn constant_power_matches_native_network_and_keeps_original_identity_and_units() {
    let network =
        LumpedNetwork::new(vec![node("z", 9.0, 0.3), node("a", 2.0, 0.8)], 298.0).unwrap();
    let density_a = constant(2.5e6);
    let density_z = constant(1.5e6);
    // Bindings deliberately reverse the original deterministic node order.
    let owner = prepared(
        &network,
        &[
            binding("z", &density_z, 2.0e-6),
            binding("a", &density_a, 2.0e-6),
        ],
    );
    let initial = [310.0, 289.0];
    let response = owner.response_at(1234.5, 3.2, &initial).unwrap();
    let native = network.response_at(&[5.0, 3.0], &initial, 3.2).unwrap();
    assert!(core::ptr::eq(owner.network(), &network));
    assert!(core::ptr::eq(response.network(), &network));
    assert_eq!(response.source_clock(), "laboratory trigger");
    assert_eq!(response.start_s(), 1234.5);
    assert_eq!(response.duration_s(), 3.2);
    assert_eq!(response.gate(), BiotGate::corpus_default());
    assert_eq!(
        response.verdict(),
        &BiotGate::corpus_default().adjudicate(&network)
    );
    for (index, result) in response.nodes().iter().enumerate() {
        assert!(core::ptr::eq(result.node(), &network.nodes()[index]));
        assert_eq!(result.input_unit(), "W/m^3");
        assert_eq!(result.watts_per_spectrum_unit(), 2.0e-6);
        assert_eq!(result.initial_temperature_k_dd().to_f64(), initial[index]);
        near(
            result.temperature_k_dd().to_f64(),
            native[index],
            1e-15,
            1e-13,
        );
        near(
            result.input_energy_j_dd().to_f64(),
            [5.0, 3.0][index] * 3.2,
            2e-15,
            0.0,
        );
        near(
            result.energy_balance_residual_j_dd().to_f64(),
            0.0,
            0.0,
            1e-28,
        );
    }
}

#[test]
fn real_carriers_drive_full_beat_and_sum_frequency_thermal_response() {
    let network = LumpedNetwork::new(vec![node("body", 2.0, 0.8)], 300.0).unwrap();
    let source = square(
        &[(2.5, [1.0, 0.0]), (2.8, [0.0, 0.7])],
        HarmonicProductLaw::Real,
    );
    let owner = prepared(&network, &[binding("body", &source, 3.2)]);
    let start = 0.31;
    let duration = 1.7;
    let actual = owner.response_at(start, duration, &[303.0]).unwrap();
    let tau = core::f64::consts::TAU;
    let expected = rk4(
        |t| {
            let field = fs_math::det::cos(tau * 2.5 * t) + 0.7 * fs_math::det::sin(tau * 2.8 * t);
            3.2 * field * field
        },
        start,
        duration,
        303.0,
        300.0,
        2.0,
        0.8,
    );
    let result = &actual.nodes()[0];
    near(result.temperature_k_dd().to_f64(), expected[0], 0.0, 2e-10);
    near(result.input_energy_j_dd().to_f64(), expected[1], 0.0, 2e-10);
    near(
        result.ambient_exchange_j_dd().to_f64(),
        expected[2],
        0.0,
        2e-10,
    );
    near(
        result.stored_energy_change_j_dd().to_f64(),
        2.0 * (expected[0] - 303.0),
        0.0,
        4e-10,
    );
    let reset = owner.response_at(0.0, duration, &[303.0]).unwrap();
    assert!(
        (result.temperature_k_dd().to_f64() - reset.nodes()[0].temperature_k_dd().to_f64()).abs()
            > 1e-3
    );
}

#[test]
fn coherent_proper_optical_carriers_keep_thermal_beats_without_fast_sampling() {
    let network = LumpedNetwork::new(vec![node("absorber", 1.2, 0.6)], 300.0).unwrap();
    let source = square(
        &[(1.0e12, [1.0, 0.0]), (1.0e12 + 0.25, [1.0, 0.0])],
        HarmonicProductLaw::ProperComplex,
    );
    let dc = constant(1.0);
    let owner = prepared(&network, &[binding("absorber", &source, 4.0)]);
    let dc_owner = prepared(&network, &[binding("absorber", &dc, 4.0)]);
    assert_eq!(source.tones().len(), 1);
    assert_eq!(owner.spectral_terms(), 2);
    let start = 0.25;
    let duration = 0.75;
    let actual = owner.response_at(start, duration, &[300.0]).unwrap();
    let expected = rk4(
        |t| 4.0 * (1.0 + fs_math::det::cos(core::f64::consts::TAU * 0.25 * t)),
        start,
        duration,
        300.0,
        300.0,
        1.2,
        0.6,
    );
    near(
        actual.nodes()[0].temperature_k_dd().to_f64(),
        expected[0],
        0.0,
        2e-10,
    );
    near(
        actual.nodes()[0].input_energy_j_dd().to_f64(),
        expected[1],
        0.0,
        2e-10,
    );
    let averaged = dc_owner.response_at(start, duration, &[300.0]).unwrap();
    assert!(
        actual.nodes()[0].temperature_k_dd().to_f64()
            > averaged.nodes()[0].temperature_k_dd().to_f64() + 0.5
    );
}

#[test]
fn successive_source_clock_intervals_compose_without_double_debiting_energy() {
    let network = LumpedNetwork::new(vec![node("body", 3.0, 0.3)], 300.0).unwrap();
    let source = square(
        &[(5.0, [1.0, 0.0]), (5.5, [0.6, 0.8])],
        HarmonicProductLaw::ProperComplex,
    );
    let owner = prepared(&network, &[binding("body", &source, 2.0)]);
    let whole = owner.response_at(0.125, 1.0, &[302.0]).unwrap();
    let first = owner.response_at(0.125, 0.375, &[302.0]).unwrap();
    let second = owner
        .response_at_dd(0.5, 0.625, &[first.nodes()[0].temperature_k_dd()])
        .unwrap();
    let a = &first.nodes()[0];
    let b = &second.nodes()[0];
    let w = &whole.nodes()[0];
    near(
        (b.temperature_k_dd() - w.temperature_k_dd()).to_f64(),
        0.0,
        0.0,
        2e-14,
    );
    for difference in [
        a.input_energy_j_dd() + b.input_energy_j_dd() - w.input_energy_j_dd(),
        a.stored_energy_change_j_dd() + b.stored_energy_change_j_dd()
            - w.stored_energy_change_j_dd(),
        a.ambient_exchange_j_dd() + b.ambient_exchange_j_dd() - w.ambient_exchange_j_dd(),
    ] {
        near(difference.to_f64(), 0.0, 0.0, 2e-14);
    }
    // A repeated observation is pure; it does not consume energy or advance state.
    let repeated = owner.response_at(0.125, 1.0, &[302.0]).unwrap();
    assert_eq!(repeated.nodes()[0].temperature_k_dd(), w.temperature_k_dd());
    assert_eq!(
        repeated.nodes()[0].input_energy_j_dd(),
        w.input_energy_j_dd()
    );
}

#[test]
fn sub_ulp_heating_accumulates_through_dd_initial_temperatures() {
    let network = LumpedNetwork::new(vec![node("body", 1.0e6, 1.0)], 300.0).unwrap();
    let source = constant(1.0);
    let owner = prepared(&network, &[binding("body", &source, 1.0)]);
    let duration = 1.0e-10;
    let first = owner.response_at(0.0, duration, &[300.0]).unwrap();
    let first_node = &first.nodes()[0];
    assert_eq!(first_node.temperature_k_dd().to_f64(), 300.0);
    near(
        first_node.temperature_change_k_dd().to_f64(),
        1.0e-16,
        1e-14,
        0.0,
    );
    assert!(first_node.temperature_k_dd().lo > 0.0);
    let mut temperature = Dd::from_f64(300.0);
    for step in 0..1_000 {
        temperature = owner
            .response_at_dd(step as f64 * duration, duration, &[temperature])
            .unwrap()
            .nodes()[0]
            .temperature_k_dd();
    }
    let whole = owner
        .response_at(0.0, 1_000.0 * duration, &[300.0])
        .unwrap();
    assert!(temperature.to_f64() > 300.0);
    near(
        (temperature - whole.nodes()[0].temperature_k_dd()).to_f64(),
        0.0,
        0.0,
        1e-25,
    );
}

#[test]
fn tiny_relaxation_keeps_energy_when_rounded_exponential_equals_one() {
    let network = LumpedNetwork::new(vec![node("body", 1.0e20, 1.0)], 300.0).unwrap();
    let source = constant(0.0);
    let owner = prepared(&network, &[binding("body", &source, 1.0)]);
    let actual = owner.response_at(0.0, 1.0e-6, &[310.0]).unwrap();
    let result = &actual.nodes()[0];
    assert_eq!(result.temperature_k_dd().to_f64(), 310.0);
    near(
        result.temperature_change_k_dd().to_f64(),
        -1e-25,
        2e-15,
        0.0,
    );
    near(
        result.stored_energy_change_j_dd().to_f64(),
        -1e-5,
        2e-15,
        0.0,
    );
    near(result.ambient_exchange_j_dd().to_f64(), 1e-5, 2e-15, 0.0);
}

#[test]
fn signed_extraction_and_inward_ambient_heat_keep_physical_signs() {
    let network = LumpedNetwork::new(vec![node("body", 4.0, 0.5)], 300.0).unwrap();
    let extraction = constant(-2.0);
    let owner = prepared(&network, &[binding("body", &extraction, 1.0)]);
    let response = owner.response_at(2.0, 3.0, &[290.0]).unwrap();
    let result = &response.nodes()[0];
    let expected = 296.0 + (290.0 - 296.0) * fs_math::det::exp(-3.0 * 0.125);
    near(result.temperature_k_dd().to_f64(), expected, 0.0, 1e-13);
    near(result.input_energy_j_dd().to_f64(), -6.0, 0.0, 1e-14);
    assert!(result.stored_energy_change_j_dd().to_f64() > 0.0);
    assert!(result.ambient_exchange_j_dd().to_f64() < -6.0);
    near(
        result.energy_balance_residual_j_dd().to_f64(),
        0.0,
        0.0,
        1e-28,
    );
}

#[test]
fn zero_duration_is_an_exact_identity_with_signed_spectral_power() {
    let network = LumpedNetwork::new(vec![node("body", 2.0, 0.5)], 300.0).unwrap();
    let source = constant(-1.0e9);
    let owner = prepared(&network, &[binding("body", &source, 1.0)]);
    let initial = Dd::from_f64(300.0) + Dd::from_f64(1.0e-20);
    let result = owner.response_at_dd(1.0e16, 0.0, &[initial]).unwrap();
    let node = &result.nodes()[0];
    assert_eq!(node.temperature_k_dd(), initial);
    assert_eq!(node.temperature_change_k_dd(), Dd::ZERO);
    assert_eq!(node.input_energy_j_dd(), Dd::ZERO);
    assert_eq!(node.stored_energy_change_j_dd(), Dd::ZERO);
    assert_eq!(node.ambient_exchange_j_dd(), Dd::ZERO);
}

#[test]
fn original_biot_and_finite_decay_gate_cannot_be_bypassed_by_a_spectrum() {
    let source = constant(1.0);
    let excessive_biot = LumpedNetwork::new(
        vec![LumpedNode::new("body", 1.0, 2.0, 1.0, 1.0, 1.0).unwrap()],
        300.0,
    )
    .unwrap();
    assert!(
        PreparedLumpedInterference::new(
            &excessive_biot,
            BiotGate::corpus_default(),
            &[binding("body", &source, 1.0)],
            limits(),
        )
        .is_err()
    );
    let explicit_gate = BiotGate::at(2.0).unwrap();
    let explicit = PreparedLumpedInterference::new(
        &excessive_biot,
        explicit_gate,
        &[binding("body", &source, 1.0)],
        limits(),
    )
    .unwrap();
    assert_eq!(explicit.gate().ceiling(), 2.0);
    assert_eq!(
        explicit.verdict(),
        &ValidityVerdict::Admitted { worst_biot: 2.0 }
    );
    for original in [
        LumpedNode::new("body", 1.0, 1.0e308, 1.0, 1.0, 1.0e-300).unwrap(),
        node("body", 1.0e300, 1.0e-300),
        node("body", 1.0e-300, 1.0e300),
    ] {
        let network = LumpedNetwork::new(vec![original], 300.0).unwrap();
        assert!(
            PreparedLumpedInterference::new(
                &network,
                BiotGate::corpus_default(),
                &[binding("body", &source, 1.0)],
                limits(),
            )
            .is_err()
        );
    }
}

#[test]
fn exact_original_node_coverage_units_clocks_and_scales_are_required() {
    let network =
        LumpedNetwork::new(vec![node("a", 2.0, 0.5), node("b", 3.0, 0.5)], 300.0).unwrap();
    let source = constant(1.0);
    let a = binding("a", &source, 1.0);
    let b = binding("b", &source, 1.0);
    let invalid = [
        [a, a],
        [a, binding("foreign", &source, 1.0)],
        [
            a,
            LumpedPowerSpectrumBinding {
                input_unit: " ",
                ..b
            },
        ],
        [
            a,
            LumpedPowerSpectrumBinding {
                source_clock: "other trigger",
                ..b
            },
        ],
        [
            a,
            LumpedPowerSpectrumBinding {
                source_clock: "",
                ..b
            },
        ],
    ];
    for bindings in invalid {
        assert!(
            PreparedLumpedInterference::new(
                &network,
                BiotGate::corpus_default(),
                &bindings,
                limits(),
            )
            .is_err()
        );
    }
    assert!(
        PreparedLumpedInterference::new(&network, BiotGate::corpus_default(), &[a], limits(),)
            .is_err()
    );
    for scale in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedLumpedInterference::new(
                &network,
                BiotGate::corpus_default(),
                &[a, binding("b", &source, scale)],
                limits(),
            )
            .is_err()
        );
    }
}

#[test]
fn node_tone_and_combined_workspace_budgets_refuse_before_results() {
    let network = LumpedNetwork::new(vec![node("body", 2.0, 0.5)], 300.0).unwrap();
    let source = square(&[(1.0, [1.0, 0.0])], HarmonicProductLaw::Real);
    for limit in [
        LumpedInterferenceLimits {
            maximum_nodes: 0,
            ..limits()
        },
        LumpedInterferenceLimits {
            maximum_spectral_terms: 1,
            ..limits()
        },
        LumpedInterferenceLimits {
            maximum_workspace_bytes: 0,
            ..limits()
        },
    ] {
        assert!(
            PreparedLumpedInterference::new(
                &network,
                BiotGate::corpus_default(),
                &[binding("body", &source, 1.0)],
                limit,
            )
            .is_err()
        );
    }
    let required = core::mem::size_of::<PreparedNodePower<'_>>()
        + core::mem::size_of::<LumpedInterferenceNodeResponse<'_>>();
    let exact = PreparedLumpedInterference::new(
        &network,
        BiotGate::corpus_default(),
        &[binding("body", &source, 1.0)],
        LumpedInterferenceLimits {
            maximum_workspace_bytes: required,
            ..limits()
        },
    )
    .unwrap();
    assert!(exact.response_at(0.0, 1.0, &[300.0]).is_ok());
}

#[test]
fn invalid_time_and_absolute_temperatures_fail_without_clipping() {
    let network = LumpedNetwork::new(vec![node("body", 2.0, 0.5)], 300.0).unwrap();
    let cooling = constant(-1.0e6);
    let owner = prepared(&network, &[binding("body", &cooling, 1.0)]);
    for (start, duration, initial) in [
        (f64::NAN, 1.0, 300.0),
        (0.0, f64::INFINITY, 300.0),
        (0.0, -1.0, 300.0),
        (0.0, 1.0, f64::NAN),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, -3.0),
        (0.0, 1.0, 300.0), // Endpoint below absolute zero.
    ] {
        assert!(owner.response_at(start, duration, &[initial]).is_err());
    }
    assert!(owner.response_at(0.0, 0.0, &[]).is_err());
    assert!(
        owner
            .response_at_dd(
                0.0,
                0.0,
                &[Dd {
                    hi: 300.0,
                    lo: f64::NAN
                }]
            )
            .is_err()
    );
    for ambient in [0.0, -1.0] {
        let invalid = LumpedNetwork::new(vec![node("body", 2.0, 0.5)], ambient).unwrap();
        assert!(
            PreparedLumpedInterference::new(
                &invalid,
                BiotGate::corpus_default(),
                &[binding("body", &cooling, 1.0)],
                limits(),
            )
            .is_err()
        );
    }
}

#[test]
fn cancellation_is_atomic_in_preparation_and_inside_spectral_response() {
    let network =
        LumpedNetwork::new(vec![node("a", 2.0, 0.5), node("b", 3.0, 0.5)], 300.0).unwrap();
    let original = network.clone();
    let source = square(
        &[(2.0, [1.0, 0.0]), (3.0, [1.0, 0.0])],
        HarmonicProductLaw::Real,
    );
    let bindings = [binding("b", &source, 1.0), binding("a", &source, 1.0)];
    assert!(matches!(
        PreparedLumpedInterference::new_cancellable(
            &network,
            BiotGate::corpus_default(),
            &bindings,
            limits(),
            || true,
        ),
        Err(ConductionError::Cancelled { .. })
    ));
    let mut preparation_polls = 0;
    assert!(matches!(
        PreparedLumpedInterference::new_cancellable(
            &network,
            BiotGate::corpus_default(),
            &bindings,
            limits(),
            || {
                preparation_polls += 1;
                preparation_polls == 5
            },
        ),
        Err(ConductionError::Cancelled { .. })
    ));
    let owner = prepared(&network, &bindings);
    let previous = owner.response_at(0.0, 0.1, &[300.0, 300.0]).unwrap();
    let retained = previous.nodes()[0].temperature_k_dd();
    let mut polls = 0;
    assert!(matches!(
        owner.response_at_cancellable(0.0, 0.1, &[300.0, 300.0], || {
            polls += 1;
            polls >= 6
        }),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(network, original);
    assert_eq!(previous.nodes()[0].temperature_k_dd(), retained);
    assert!(owner.response_at(0.0, 0.1, &[300.0, 300.0]).is_ok());
}
