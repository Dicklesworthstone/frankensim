//! Independent calorimetry, NTU and manufactured-Q exchanger regressions.

mod phase;
mod response;

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityQuery;
use fs_qty::{Pressure, QuantitySpec, Temperature};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 719,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new_clock_free(), f)
}

fn source(knots: &[(f64, f64)], datum: f64, pressure: f64) -> Arc<IsobaricHeatCapacityLaw> {
    let lower = knots[0].0;
    let upper = knots[knots.len() - 1].0;
    let key = PropertyKey::new("Cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    lower,
                    upper,
                )
                .with_quantity(
                    "p",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    pressure,
                    pressure,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "independent finite-stream calorimetry fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured exchanger fluid".into(),
            phase: "single-phase".into(),
            process: "fixed source isobar".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |temperature| {
        QueryPoint::new()
            .with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                temperature,
            )
            .unwrap()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), pressure)
            .unwrap()
    };
    Arc::new(
        IsobaricHeatCapacityLaw::try_new(
            &card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &key,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &point(lower),
                upper: &point(upper),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            datum,
        )
        .unwrap(),
    )
}

fn affine(intercept: f64, slope: f64, datum: f64, pressure: f64) -> Arc<IsobaricHeatCapacityLaw> {
    source(
        &[250.0, 350.0, 450.0, 550.0, 700.0].map(|t| (t, intercept + slope * t)),
        datum,
        pressure,
    )
}

fn stream(initial: f64, mass: f64, law: Arc<IsobaricHeatCapacityLaw>) -> ExchangerStream {
    ExchangerStream::new(initial, mass, law).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}

fn audit(state: ExchangerState, exchanger: &CaloricHeatExchanger) {
    let config = ExchangerConfig::default();
    let heat = state.heat_rate_first_to_second_w;
    close(
        exchanger.first.mass_flow_kg_per_s * state.first_specific_enthalpy_change_j_kg,
        -heat,
        8.0 * f64::EPSILON * heat.abs(),
    );
    close(
        exchanger.second.mass_flow_kg_per_s * state.second_specific_enthalpy_change_j_kg,
        heat,
        8.0 * f64::EPSILON * heat.abs(),
    );
    assert!(
        state.first_endpoint_energy_residual_w.abs()
            <= config.heat_tolerance(heat).unwrap()
                + state.first_endpoint_energy_roundoff_allowance_w
    );
    assert!(
        state.second_endpoint_energy_residual_w.abs()
            <= config.heat_tolerance(heat).unwrap()
                + state.second_endpoint_energy_roundoff_allowance_w
    );
    assert!(state.entropy_generation_w_per_k >= -state.entropy_roundoff_allowance_w_per_k);
    assert!(state.heat_error_bound_w <= config.heat_tolerance(heat).unwrap());
    assert!(
        state.conductance_residual_w_per_k.abs() + state.conductance_error_allowance_w_per_k
            <= config.ua_tolerance(exchanger.conductance_w_per_k).unwrap()
    );
}

#[test]
fn constant_sources_match_independent_ntu_and_reverse_reciprocally() {
    with_cx(|cx| {
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            let first = stream(600.0, 1.0, affine(1000.0, 0.0, 400.0, 100_000.0));
            let second = stream(300.0, 1.5, affine(2000.0, 0.0, 650.0, 250_000.0));
            let exchanger =
                CaloricHeatExchanger::new(first.clone(), second.clone(), arrangement, 750.0)
                    .unwrap();
            let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
            let ntu: f64 = 0.75;
            let ratio: f64 = 1.0 / 3.0;
            let effectiveness = match arrangement {
                FlowArrangement::CoCurrent => -(-ntu * (1.0 + ratio)).exp_m1() / (1.0 + ratio),
                FlowArrangement::Counterflow => {
                    let e = (-ntu * (1.0 - ratio)).exp();
                    (1.0 - e) / (1.0 - ratio * e)
                }
            };
            let expected = effectiveness * 300_000.0;
            close(state.heat_rate_first_to_second_w, expected, 2e-8);
            close(
                state.first_outlet_temperature_k,
                600.0 - expected / 1000.0,
                2e-11,
            );
            close(
                state.second_outlet_temperature_k,
                300.0 + expected / 3000.0,
                2e-11,
            );
            close(
                state.entropy_generation_w_per_k,
                1000.0 * (state.first_outlet_temperature_k / 600.0).ln()
                    + 3000.0 * (state.second_outlet_temperature_k / 300.0).ln(),
                2e-10,
            );
            assert_eq!(state.iterations, 0);
            audit(state, &exchanger);
            let reversed = CaloricHeatExchanger::new(second, first, arrangement, 750.0).unwrap();
            let reverse = reversed.solve(cx, ExchangerConfig::default()).unwrap();
            assert_eq!(
                reverse.heat_rate_first_to_second_w,
                -state.heat_rate_first_to_second_w
            );
            assert_eq!(
                reverse.first_outlet_temperature_k,
                state.second_outlet_temperature_k
            );
            assert_eq!(
                reverse.second_outlet_temperature_k,
                state.first_outlet_temperature_k
            );
            close(
                reverse.entropy_generation_w_per_k,
                state.entropy_generation_w_per_k,
                1e-12,
            );
            audit(reverse, &reversed);
        }
    });
}

#[test]
fn equal_and_nearly_equal_capacity_counterflow_keep_the_correct_limit() {
    with_cx(|cx| {
        let law = affine(1000.0, 0.0, 300.0, 100_000.0);
        for ratio_change in [0.0, 1e-12, 1e-7] {
            let exchanger = CaloricHeatExchanger::new(
                stream(600.0, 1.0, law.clone()),
                stream(300.0, 1.0 + ratio_change, law.clone()),
                FlowArrangement::Counterflow,
                1700.0,
            )
            .unwrap();
            let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
            let a = 1.0 / 1000.0;
            let b = 1.0 / (1000.0 * (1.0 + ratio_change));
            // Independent constant-Cp boundary-value solution, expressed in
            // reciprocal capacity rates rather than epsilon/capacity ratio.
            let exponent = (1700.0_f64 * (a - b)).exp_m1();
            let expected = if ratio_change == 0.0 {
                300.0 / (1.0 / 1700.0 + a)
            } else {
                300.0 * exponent / ((a - b) + a * exponent)
            };
            close(state.heat_rate_first_to_second_w, expected, 2e-7);
            audit(state, &exchanger);
        }
    });
}

#[test]
fn weak_exchange_retains_heat_and_entropy_when_outlet_temperatures_round() {
    with_cx(|cx| {
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            for ua in [1e-20, 1e-13] {
                let exchanger = CaloricHeatExchanger::new(
                    stream(500.0, 1.0, affine(1000.0, 0.0, 700.0, 100_000.0)),
                    stream(300.0, 1.0, affine(1000.0, 0.0, 250.0, 400_000.0)),
                    arrangement,
                    ua,
                )
                .unwrap();
                let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
                let expected_heat = ua * 200.0;
                close(
                    state.heat_rate_first_to_second_w,
                    expected_heat,
                    expected_heat * 2e-14,
                );
                assert_eq!(state.first_outlet_temperature_k, 500.0);
                assert_eq!(state.second_outlet_temperature_k, 300.0);
                let entropy = expected_heat * (1.0 / 300.0 - 1.0 / 500.0);
                close(state.entropy_generation_w_per_k, entropy, entropy * 2e-14);
                assert!(state.first_endpoint_energy_residual_w > 0.0);
                audit(state, &exchanger);
            }
        }
    });
}

// Independent quadratic primitive for Cp(T)=a+bT. It never queries the native
// material inverse, heat-exchanger residual or quadrature implementation.
fn after_heat(initial: f64, heat: f64, a: f64, b: f64) -> f64 {
    let cp = a + b * initial;
    initial + 2.0 * heat / (cp + (cp * cp + 2.0 * b * heat).sqrt())
}

fn manufactured_ua(heat: f64, arrangement: FlowArrangement, panels: usize) -> f64 {
    let integrand = |q: f64| {
        let hot = after_heat(620.0, -q / 0.8, 600.0, 2.0);
        let progress = if arrangement == FlowArrangement::CoCurrent {
            q
        } else {
            heat - q
        };
        let cold = after_heat(280.0, progress / 1.1, 1800.0, -1.0);
        1.0 / (hot - cold)
    };
    let step = heat / panels as f64;
    let mut sum = integrand(0.0) + integrand(heat);
    for i in 1..panels {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(step * i as f64);
    }
    step * sum / 3.0
}

#[test]
fn variable_original_cp_matches_manufactured_heat_and_independent_entropy() {
    with_cx(|cx| {
        let expected_heat = 180_000.0;
        for arrangement in [FlowArrangement::CoCurrent, FlowArrangement::Counterflow] {
            let ua = manufactured_ua(expected_heat, arrangement, 16_384);
            close(ua, manufactured_ua(expected_heat, arrangement, 8192), 2e-9);
            let mut reference: Option<ExchangerState> = None;
            for (hot_datum, cold_datum) in [(250.0, 700.0), (700.0, 250.0)] {
                let exchanger = CaloricHeatExchanger::new(
                    stream(620.0, 0.8, affine(600.0, 2.0, hot_datum, 100_000.0)),
                    stream(280.0, 1.1, affine(1800.0, -1.0, cold_datum, 300_000.0)),
                    arrangement,
                    ua,
                )
                .unwrap();
                let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
                close(state.heat_rate_first_to_second_w, expected_heat, 4e-5);
                let expected_hot = after_heat(620.0, -expected_heat / 0.8, 600.0, 2.0);
                let expected_cold = after_heat(280.0, expected_heat / 1.1, 1800.0, -1.0);
                close(state.first_outlet_temperature_k, expected_hot, 4e-8);
                close(state.second_outlet_temperature_k, expected_cold, 4e-8);
                let entropy = 0.8
                    * (600.0 * (expected_hot / 620.0).ln() + 2.0 * (expected_hot - 620.0))
                    + 1.1 * (1800.0 * (expected_cold / 280.0).ln() - (expected_cold - 280.0));
                close(state.entropy_generation_w_per_k, entropy, 1e-7);
                assert!(state.iterations > 0);
                audit(state, &exchanger);
                if let Some(previous) = reference {
                    close(
                        state.heat_rate_first_to_second_w,
                        previous.heat_rate_first_to_second_w,
                        4e-5,
                    );
                    close(
                        state.entropy_generation_w_per_k,
                        previous.entropy_generation_w_per_k,
                        1e-7,
                    );
                }
                reference = Some(state);
            }
        }
    });
}

#[test]
fn actual_cp_kinks_match_an_independently_split_manufactured_integral() {
    with_cx(|cx| {
        let target = 320_000.0;
        let hot_temperature = |q: f64| {
            if q <= 280_000.0 {
                650.0 - 2.0 * q / (1800.0 + (1800.0 * 1800.0 - 8.0 * q).sqrt())
            } else {
                450.0 - (q - 280_000.0) / 1000.0
            }
        };
        let cold_temperature = |q: f64| {
            let heat = q / 1.3;
            if heat <= 100_100.0 {
                280.0 + 2.0 * heat / (1360.0 + (1360.0 * 1360.0 + 4.0 * heat).sqrt())
            } else {
                let extra = heat - 100_100.0;
                350.0 + 2.0 * extra / (1500.0 + (1500.0 * 1500.0 - 2.0 * extra).sqrt())
            }
        };
        let integrand = |q: f64| 1.0 / (hot_temperature(q) - cold_temperature(target - q));
        // These heat breakpoints follow independent trapezoidal calorimetry
        // on the declared Cp pieces, not the production source-knot index.
        let breaks = [0.0, target - 1.3 * 100_100.0, 280_000.0, target];
        let mut ua = 0.0;
        for interval in breaks.windows(2) {
            let step = (interval[1] - interval[0]) / 8192.0;
            let mut sum = integrand(interval[0]) + integrand(interval[1]);
            for i in 1..8192 {
                sum +=
                    (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(interval[0] + step * i as f64);
            }
            ua += step * sum / 3.0;
        }
        let exchanger = CaloricHeatExchanger::new(
            stream(
                650.0,
                1.0,
                source(
                    &[(250.0, 1000.0), (450.0, 1000.0), (700.0, 2000.0)],
                    450.0,
                    100_000.0,
                ),
            ),
            stream(
                280.0,
                1.3,
                source(
                    &[(250.0, 1300.0), (350.0, 1500.0), (700.0, 1150.0)],
                    350.0,
                    200_000.0,
                ),
            ),
            FlowArrangement::Counterflow,
            ua,
        )
        .unwrap();
        let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
        close(state.heat_rate_first_to_second_w, target, 6e-5);
        close(
            state.first_outlet_temperature_k,
            hot_temperature(target),
            7e-8,
        );
        close(
            state.second_outlet_temperature_k,
            cold_temperature(target),
            7e-8,
        );
        audit(state, &exchanger);
    });
}

#[test]
fn counterflow_rejects_an_interior_pinch_with_positive_endpoint_gaps() {
    with_cx(|cx| {
        let hot = stream(
            600.0,
            1.0,
            source(
                &[(250.0, 360.0), (400.0, 210.0), (600.0, 10.0)],
                500.0,
                100_000.0,
            ),
        );
        let cold = stream(
            300.0,
            1.0,
            source(
                &[(300.0, 10.0), (450.0, 160.0), (650.0, 360.0)],
                450.0,
                200_000.0,
            ),
        );
        let config = ExchangerConfig::default();
        let pair = Pair {
            hot: SourcePath::new(cx, &hot, 300.0, -1.0, config).unwrap(),
            cold: SourcePath::new(cx, &cold, 600.0, 1.0, config).unwrap(),
            arrangement: FlowArrangement::Counterflow,
            driving: 300.0,
        };
        let q = 40_000.0;
        let (ha, ca) = pair.point(q, 0.0).unwrap();
        let (hb, cb) = pair.point(q, q).unwrap();
        assert!(ha.temperature > ca.temperature && hb.temperature > cb.temperature);
        assert_eq!(
            pair.measure(cx, q, 1e-6, config, None).unwrap_err(),
            ExchangerError::Pinch
        );
        // Independent symmetric inverse gives the first interior pinch at
        // Q=2*(10*150+150^2/2)=25500 W, despite the 48000 W endpoint bound.
        let exchanger =
            CaloricHeatExchanger::new(hot, cold, FlowArrangement::Counterflow, 100.0).unwrap();
        let state = exchanger.solve(cx, config).unwrap();
        assert!(state.heat_rate_first_to_second_w < 25_500.0);
        audit(state, &exchanger);
    });
}

#[test]
fn each_actual_outlet_requires_its_own_source_support() {
    with_cx(|cx| {
        let hot_limited = stream(
            400.0,
            1.0,
            source(&[(390.0, 1780.0), (400.0, 1800.0)], 395.0, 100_000.0),
        );
        let cold = stream(300.0, 1.0, affine(1200.0, 1.0, 400.0, 200_000.0));
        let valid = CaloricHeatExchanger::new(
            hot_limited.clone(),
            cold.clone(),
            FlowArrangement::Counterflow,
            20.0,
        )
        .unwrap();
        audit(valid.solve(cx, ExchangerConfig::default()).unwrap(), &valid);
        let outside =
            CaloricHeatExchanger::new(hot_limited, cold, FlowArrangement::Counterflow, 1000.0)
                .unwrap();
        assert_eq!(
            outside.solve(cx, ExchangerConfig::default()).unwrap_err(),
            ExchangerError::OutsideSourceSupport
        );
        let cold_limited = stream(
            300.0,
            1.0,
            source(&[(300.0, 1500.0), (305.0, 1505.0)], 302.0, 300_000.0),
        );
        let outside = CaloricHeatExchanger::new(
            stream(400.0, 1.0, affine(1000.0, 0.0, 250.0, 100_000.0)),
            cold_limited,
            FlowArrangement::CoCurrent,
            1000.0,
        )
        .unwrap();
        assert_eq!(
            outside.solve(cx, ExchangerConfig::default()).unwrap_err(),
            ExchangerError::OutsideSourceSupport
        );
    });
}

#[test]
fn zero_transfer_invalid_inputs_and_resolved_saturation_are_explicit() {
    with_cx(|cx| {
        let law = affine(1000.0, 0.0, 400.0, 100_000.0);
        for (second_inlet, ua) in [(300.0, 0.0), (500.0, 500.0)] {
            let exchanger = CaloricHeatExchanger::new(
                stream(500.0, 1.0, law.clone()),
                stream(second_inlet, 1.0, law.clone()),
                FlowArrangement::CoCurrent,
                ua,
            )
            .unwrap();
            let state = exchanger.solve(cx, ExchangerConfig::default()).unwrap();
            assert_eq!(state.heat_rate_first_to_second_w, 0.0);
            assert_eq!(state.entropy_generation_w_per_k, 0.0);
            audit(state, &exchanger);
        }
        assert!(ExchangerStream::new(500.0, 0.0, law.clone()).is_err());
        assert!(ExchangerStream::new(800.0, 1.0, law.clone()).is_err());
        assert!(
            CaloricHeatExchanger::new(
                stream(500.0, 1.0, law.clone()),
                stream(300.0, 1.0, law.clone()),
                FlowArrangement::Counterflow,
                f64::NAN
            )
            .is_err()
        );
        let saturated = CaloricHeatExchanger::new(
            stream(500.0, 1.0, law.clone()),
            stream(300.0, 1.0, law),
            FlowArrangement::CoCurrent,
            1e9,
        )
        .unwrap();
        assert_eq!(
            saturated.solve(cx, ExchangerConfig::default()).unwrap_err(),
            ExchangerError::Pinch
        );
    });
}

#[test]
fn bounded_work_nonconvergence_and_cancellation_publish_no_state() {
    let exchanger = CaloricHeatExchanger::new(
        stream(620.0, 0.8, affine(600.0, 2.0, 250.0, 100_000.0)),
        stream(280.0, 1.1, affine(1800.0, -1.0, 700.0, 300_000.0)),
        FlowArrangement::CoCurrent,
        350.0,
    )
    .unwrap();
    with_cx(|cx| {
        assert_eq!(
            exchanger
                .solve(
                    cx,
                    ExchangerConfig {
                        max_quadrature_panels: 2,
                        ..ExchangerConfig::default()
                    }
                )
                .unwrap_err(),
            ExchangerError::WorkLimit
        );
        assert_eq!(
            exchanger
                .solve(
                    cx,
                    ExchangerConfig {
                        max_iterations: 1,
                        ..ExchangerConfig::default()
                    }
                )
                .unwrap_err(),
            ExchangerError::NonConvergence { iterations: 1 }
        );
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert_eq!(
            exchanger.solve(cx, ExchangerConfig::default()).unwrap_err(),
            ExchangerError::Interrupted
        );
    });
    with_cx(|cx| {
        audit(
            exchanger.solve(cx, ExchangerConfig::default()).unwrap(),
            &exchanger,
        )
    });
}
