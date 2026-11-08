use super::super::co_current_water::*;
use super::super::phase::{WaterPhaseExchangerConfig, WaterPhaseHeatExchanger};
use super::*;
use fs_material::phase::water::{IsobaricWater, WaterChartControls, WaterPhase, WaterState};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn water_window(minimum: f64, maximum: f64) -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(1e6),
        WaterChartControls {
            minimum_temperature_k: minimum,
            maximum_temperature_k: maximum,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )
    .unwrap()
}
fn water() -> IsobaricWater {
    water_window(273.16, 800.0)
}
fn at(chart: &IsobaricWater, t: f64) -> WaterState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn cp_stream(t: f64, mass: f64, cp: f64, datum: f64) -> ExchangerStream {
    stream(t, mass, source(&[(250.0, cp), (1000.0, cp)], datum, 2e5))
}
fn controls(tolerance: f64) -> CoCurrentWaterExchangerConfig {
    CoCurrentWaterExchangerConfig {
        enthalpy_absolute_tolerance_j_kg: tolerance,
        temperature_absolute_tolerance_k: None,
        relative_tolerance: 0.0,
        ..CoCurrentWaterExchangerConfig::default()
    }
}
fn audit(exchanger: &CoCurrentWaterHeatExchanger, state: &CoCurrentWaterExchangerState) {
    let q = state.heat_rate_water_to_sensible_w;
    let tolerance = WaterPhaseExchangerConfig::default().relative_heat_tolerance * q.abs() + 1e-5;
    for residual in state
        .water_source_energy_residuals_w
        .into_iter()
        .chain(state.water_endpoint_energy_residuals_w)
    {
        assert!(residual.abs() <= tolerance);
    }
    assert!(
        state.sensible_endpoint_energy_residual_w.abs()
            <= tolerance + state.sensible_endpoint_roundoff_allowance_w
    );
    close(
        state.water_internal_energy_change_rate_w + state.water_pressure_flow_work_change_rate_w,
        -q,
        tolerance,
    );
    close(
        exchanger.water_mass_flow_kg_per_s() * state.water_specific_enthalpy_change_j_kg,
        -q,
        8.0 * f64::EPSILON * q.abs(),
    );
    close(
        exchanger.sensible_stream().mass_flow_kg_per_s()
            * state.sensible_specific_enthalpy_change_j_kg,
        q,
        8.0 * f64::EPSILON * q.abs(),
    );
    assert!(state.thermal_entropy_generation_w_per_k >= -state.entropy_numerical_allowance_w_per_k);
    assert!(state.source_entropy_generation_w_per_k >= -state.entropy_numerical_allowance_w_per_k);
    close(
        state.source_entropy_generation_w_per_k
            + state.water_thermal_entropy_correction_rate_w_per_k,
        state.thermal_entropy_generation_w_per_k,
        2e-11,
    );
    if let Some(report) = &state.integration {
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.work.energy.incoming_energy_j, 0.0);
        assert_eq!(report.work.energy.external_heat_j, 0.0);
        assert!(report.work.energy.residual_j.abs() <= 1e-5);
    }
}

// Independent enthalpy-coordinate Simpson rating. The production calculation
// instead evolves the original two-body energy equations with native SSPRK2.
// Split exactly at both ORIGINAL latent endpoints before fixed-grid quadrature.
fn required_ua(
    initial: WaterState,
    final_state: WaterState,
    water_mass: f64,
    tin: f64,
    capacity_rate: f64,
) -> (f64, f64) {
    required_ua_with_temperature(
        initial,
        final_state,
        water_mass,
        |q| tin + q / capacity_rate,
        &[],
    )
}

fn required_ua_with_temperature(
    initial: WaterState,
    final_state: WaterState,
    water_mass: f64,
    sensible_temperature: impl Fn(f64) -> f64,
    extra_breaks: &[f64],
) -> (f64, f64) {
    let h0 = initial.specific_enthalpy_j_kg();
    let heat = water_mass * (h0 - final_state.specific_enthalpy_j_kg());
    let sign = heat.signum();
    let magnitude = heat.abs();
    let mut breaks = vec![0.0, magnitude];
    breaks.extend(
        extra_breaks
            .iter()
            .copied()
            .filter(|q| *q > 0.0 && *q < magnitude),
    );
    for quality in [0.0, 1.0] {
        let h = initial
            .law()
            .state_at_quality(quality, || false)
            .unwrap()
            .specific_enthalpy_j_kg();
        let q = sign * water_mass * (h0 - h);
        if q > 0.0 && q < magnitude {
            breaks.push(q);
        }
    }
    breaks.sort_by(f64::total_cmp);
    let integrand = |q: f64| {
        let water = initial
            .law()
            .state_at_specific_enthalpy(h0 - sign * q / water_mass, || false)
            .unwrap()
            .into_state();
        let gap = sign * (water.temperature_k() - sensible_temperature(sign * q));
        assert!(gap > 0.0);
        1.0 / gap
    };
    let mut ua = 0.0;
    for range in breaks.windows(2) {
        let count = 2048;
        let step = (range[1] - range[0]) / count as f64;
        let mut sum = integrand(range[0]) + integrand(range[1]);
        for i in 1..count {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(range[0] + i as f64 * step);
        }
        ua += step * sum / 3.0;
    }
    (heat, ua)
}

#[test]
fn original_nonlinear_cp_crosses_knots_and_requires_only_its_actual_temperature_path() {
    let chart = water();
    let initial = at(&chart, 700.0);
    let target = at(&chart, 350.0);
    let knots = [
        (250.0, 3800.0),
        (310.0, 4400.0),
        (330.0, 3600.0),
        (400.0, 5000.0),
    ];
    // Independent trapezoidal calorimetry and stable quadratic inverses of
    // the three authored affine-Cp intervals. No production caloric inverse.
    let segments = [
        (300.0, 310.0, 4300.0, 10.0),
        (310.0, 330.0, 4400.0, -40.0),
        (330.0, 400.0, 3600.0, 20.0),
    ];
    let temperature = |heat: f64| {
        let mut remaining = heat / 20.0;
        for (a, b, cp, slope) in segments {
            let width = b - a;
            let amount = (cp + 0.5 * slope * width) * width;
            if remaining <= amount {
                return a + 2.0 * remaining / (cp + (cp * cp + 2.0 * slope * remaining).sqrt());
            }
            remaining -= amount;
        }
        panic!("reference exhausted authored sensible support");
    };
    let breaks = [20.0 * 43_500.0, 20.0 * (43_500.0 + 80_000.0)];
    let (heat, ua) = required_ua_with_temperature(initial, target, 1.0, temperature, &breaks);
    let exchanger = CoCurrentWaterHeatExchanger::new(
        stream(300.0, 20.0, source(&knots, 350.0, 2e5)),
        initial,
        1.0,
        ua,
    )
    .unwrap();
    let state = with_cx(|cx| exchanger.solve(cx, controls(1e-4)).unwrap());
    assert!(
        initial.temperature_k()
            > exchanger
                .sensible_stream()
                .caloric_law()
                .temperature_span_k()[1]
    );
    assert!(state.sensible_outlet_temperature_k > 330.0);
    close(state.heat_rate_water_to_sensible_w, heat, 0.5);
    close(state.sensible_outlet_temperature_k, temperature(heat), 1e-5);
    close(state.water_outlet.temperature_k(), 350.0, 2e-4);
    audit(&exchanger, &state);
}

#[test]
fn continuous_plateau_matches_the_original_exact_condenser_and_boiler() {
    let chart = water();
    for (tin, cp, mass, quality, water_mass, ua) in [
        (300.0, 4200.0, 1.5, 0.8, 0.4, 4000.0),
        (550.0, 1100.0, 2.0, 0.1, 0.25, 1800.0),
    ] {
        let inlet = chart.state_at_quality(quality, || false).unwrap();
        let sensible = cp_stream(tin, mass, cp, 350.0);
        let exact =
            WaterPhaseHeatExchanger::new(sensible.clone(), inlet.into(), water_mass, ua).unwrap();
        let exchanger = CoCurrentWaterHeatExchanger::new(sensible, inlet, water_mass, ua).unwrap();
        let expected = with_cx(|cx| {
            exact
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap()
        });
        let state = with_cx(|cx| exchanger.solve(cx, controls(1e-4)).unwrap());
        close(
            state.heat_rate_water_to_sensible_w,
            expected.heat_rate_water_to_sensible_w,
            0.1,
        );
        close(
            state.water_outlet.vapor_quality().unwrap(),
            expected.water_outlet.vapor_quality().unwrap(),
            5e-7,
        );
        assert_eq!(state.water_outlet.temperature_k(), inlet.temperature_k());
        audit(&exchanger, &state);
    }
}

#[test]
fn full_desuperheating_condensation_and_subcooling_match_independent_ua_with_refinement() {
    let chart = water();
    let initial = at(&chart, 700.0);
    let target = at(&chart, 350.0);
    let (heat, ua) = required_ua(initial, target, 1.0, 300.0, 20.0 * 4200.0);
    let exchanger =
        CoCurrentWaterHeatExchanger::new(cp_stream(300.0, 20.0, 4200.0, 350.0), initial, 1.0, ua)
            .unwrap();
    let mut errors = vec![];
    for accuracy in [0.5, 1e-4] {
        let state = with_cx(|cx| exchanger.solve(cx, controls(accuracy)).unwrap());
        errors.push((state.heat_rate_water_to_sensible_w - heat).abs());
        assert_eq!(state.water_outlet.phase(), WaterPhase::CompressedLiquid);
        assert!(state.water_pressure_flow_work_change_rate_w < -100_000.0);
        assert!(state.sensible_outlet_temperature_k < state.water_outlet.temperature_k());
        audit(&exchanger, &state);
    }
    assert!(errors[1] < errors[0] / 8.0, "heat refinement {errors:?}");
    assert!(errors[1] < 0.5, "fine heat error {errors:?}");
}

#[test]
fn feedwater_boiling_and_superheating_match_independent_ua_and_preserve_mass_rate_scaling() {
    let chart = water();
    let initial = at(&chart, 300.0);
    let target = at(&chart, 600.0);
    let (heat, ua) = required_ua(initial, target, 1.0, 900.0, 20.0 * 1500.0);
    let mut outcomes = vec![];
    for (scale, datum) in [(1.0, 300.0), (2.0, 700.0)] {
        let exchanger = CoCurrentWaterHeatExchanger::new(
            cp_stream(900.0, 20.0 * scale, 1500.0, datum),
            initial,
            scale,
            ua * scale,
        )
        .unwrap();
        let mut config = controls(1e-4);
        config.endpoint.absolute_heat_tolerance_w *= scale;
        let state = with_cx(|cx| exchanger.solve(cx, config).unwrap());
        close(state.heat_rate_water_to_sensible_w / scale, heat, 0.5);
        close(state.water_outlet.temperature_k(), 600.0, 3e-4);
        assert_eq!(state.water_outlet.phase(), WaterPhase::SuperheatedVapor);
        assert!(state.water_pressure_flow_work_change_rate_w > 100_000.0 * scale);
        outcomes.push(state);
    }
    close(
        outcomes[0].water_outlet.specific_enthalpy_j_kg(),
        outcomes[1].water_outlet.specific_enthalpy_j_kg(),
        1e-3,
    );
    close(
        outcomes[0].sensible_outlet_temperature_k,
        outcomes[1].sensible_outlet_temperature_k,
        1e-7,
    );
}

#[test]
fn zero_exchange_source_exhaustion_bounded_work_and_cancellation_publish_no_partial_state() {
    let chart = water();
    let initial = at(&chart, 700.0);
    for (tin, ua) in [(300.0, 0.0), (700.0, 1000.0)] {
        let exchanger =
            CoCurrentWaterHeatExchanger::new(cp_stream(tin, 1.0, 2000.0, 300.0), initial, 1.0, ua)
                .unwrap();
        let state = with_cx(|cx| exchanger.solve(cx, controls(1e-3)).unwrap());
        assert_eq!(state.water_outlet, initial);
        assert_eq!(state.sensible_outlet_temperature_k, tin);
        assert_eq!(state.heat_rate_water_to_sensible_w, 0.0);
        assert!(state.integration.is_none());
        audit(&exchanger, &state);
    }
    let exchanger = CoCurrentWaterHeatExchanger::new(
        cp_stream(300.0, 20.0, 4200.0, 350.0),
        initial,
        1.0,
        10_000.0,
    )
    .unwrap();
    with_cx(|cx| {
        assert!(matches!(
            exchanger.solve(
                cx,
                CoCurrentWaterExchangerConfig {
                    maximum_stage_evaluations: 1,
                    ..controls(1e-3)
                }
            ),
            Err(CoCurrentWaterExchangerError::Thermal(_))
        ))
    });
    let limited = CoCurrentWaterHeatExchanger::new(
        stream(
            300.0,
            1.0,
            source(&[(290.0, 2000.0), (310.0, 2000.0)], 300.0, 1e5),
        ),
        initial,
        1.0,
        1000.0,
    )
    .unwrap();
    with_cx(|cx| {
        assert!(matches!(
            limited.solve(cx, controls(1e-3)),
            Err(CoCurrentWaterExchangerError::Thermal(_))
        ))
    });
    assert_eq!(limited.water_inlet(), initial);
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert_eq!(
            exchanger.solve(cx, controls(1e-3)).unwrap_err(),
            CoCurrentWaterExchangerError::Interrupted
        )
    });
}
