use super::super::countercurrent_water::*;
use super::super::phase::{WaterPhaseExchangerConfig, WaterPhaseHeatExchanger};
use super::*;
use fs_material::phase::water::{IsobaricWater, WaterChartControls, WaterPhase, WaterState};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn water() -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(1e6),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 800.0,
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
fn at(chart: &IsobaricWater, t: f64) -> WaterState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn cp_stream(t: f64, mass: f64, cp: f64, datum: f64) -> ExchangerStream {
    stream(t, mass, source(&[(250.0, cp), (1000.0, cp)], datum, 2e5))
}
fn audit(owner: &CountercurrentWaterHeatExchanger, state: &CountercurrentWaterExchangerState) {
    let balance = &state.balance;
    let heat = balance.heat_rate_water_to_sensible_w;
    let tolerance = 1e-5 + 1e-10 * heat.abs();
    for residual in balance
        .water_source_energy_residuals_w
        .into_iter()
        .chain(balance.water_endpoint_energy_residuals_w)
    {
        assert!(residual.abs() <= tolerance, "source energy {residual:e}");
    }
    assert!(
        balance.sensible_endpoint_energy_residual_w.abs()
            <= tolerance + balance.sensible_endpoint_roundoff_allowance_w
    );
    close(
        balance.water_internal_energy_change_rate_w
            + balance.water_pressure_flow_work_change_rate_w,
        -heat,
        tolerance,
    );
    close(
        owner.water_mass_flow_kg_per_s() * balance.water_specific_enthalpy_change_j_kg,
        -heat,
        8.0 * f64::EPSILON * heat.abs(),
    );
    close(
        owner.sensible_stream().mass_flow_kg_per_s()
            * balance.sensible_specific_enthalpy_change_j_kg,
        heat,
        8.0 * f64::EPSILON * heat.abs(),
    );
    assert!(
        balance.source_entropy_generation_w_per_k >= -balance.entropy_numerical_allowance_w_per_k
    );
    assert!(
        balance.thermal_entropy_generation_w_per_k >= -balance.entropy_numerical_allowance_w_per_k
    );
    close(
        balance.source_entropy_generation_w_per_k
            + balance.water_thermal_entropy_correction_rate_w_per_k,
        balance.thermal_entropy_generation_w_per_k,
        2e-11,
    );
    assert!(balance.integration.is_none());
    assert!(state.minimum_temperature_gap_lower_bound_k > 0.0 || heat == 0.0);
    assert!(state.quadrature_error_estimate_w_per_k >= 0.0);
    assert!(state.source_inverse_allowance_w_per_k >= 0.0);
}

// Independent FIXED-grid Simpson rating, split at original phase endpoints.
// Opposing sensible temperature is supplied independently of SourcePath and
// the production adaptive-panel implementation; its argument is signed heat
// accumulated from that stream's inlet, which is at the WATER outlet.
fn required_ua(
    initial: WaterState,
    target: WaterState,
    water_mass: f64,
    sensible_temperature: impl Fn(f64) -> f64,
    extra_progress_breaks: &[f64],
) -> (f64, f64) {
    let h0 = initial.specific_enthalpy_j_kg();
    let heat = water_mass * (h0 - target.specific_enthalpy_j_kg());
    let sign = heat.signum();
    let total = heat.abs();
    let mut breaks = vec![0.0, total];
    for quality in [0.0, 1.0] {
        let h = initial
            .law()
            .state_at_quality(quality, || false)
            .unwrap()
            .specific_enthalpy_j_kg();
        let q = sign * water_mass * (h0 - h);
        if q > 0.0 && q < total {
            breaks.push(q);
        }
    }
    breaks.extend(
        extra_progress_breaks
            .iter()
            .map(|p| total - p)
            .filter(|q| *q > 0.0 && *q < total),
    );
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    let value = |q: f64| {
        let water = initial
            .law()
            .state_at_specific_enthalpy(h0 - sign * q / water_mass, || false)
            .unwrap()
            .into_state();
        let gap = sign * (water.temperature_k() - sensible_temperature(sign * (total - q)));
        assert!(gap > 0.0, "reference interior gap {gap:e}");
        1.0 / gap
    };
    let mut ua = 0.0;
    for range in breaks.windows(2) {
        let count = 2048;
        let step = (range[1] - range[0]) / count as f64;
        let mut sum = value(range[0]) + value(range[1]);
        for i in 1..count {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * value(range[0] + i as f64 * step);
        }
        ua += step * sum / 3.0;
    }
    (heat, ua)
}

#[test]
fn full_condensation_and_subcooling_allow_countercurrent_outlet_crossover() {
    let chart = water();
    let inlet = at(&chart, 700.0);
    let target = at(&chart, 350.0);
    let (heat, ua) = required_ua(inlet, target, 1.0, |q| 300.0 + q / (5.0 * 4200.0), &[]);
    let owner =
        CountercurrentWaterHeatExchanger::new(cp_stream(300.0, 5.0, 4200.0, 350.0), inlet, 1.0, ua)
            .unwrap();
    let config = CountercurrentWaterExchangerConfig::default();
    let rated = with_cx(|cx| owner.rate_duty(cx, heat, config).unwrap());
    close(rated.required_conductance_w_per_k, ua, 2e-5);
    let state = with_cx(|cx| owner.solve(cx, config).unwrap());
    close(state.balance.heat_rate_water_to_sensible_w, heat, 0.05);
    close(state.balance.water_outlet.temperature_k(), 350.0, 2e-5);
    assert_eq!(
        state.balance.water_outlet.phase(),
        WaterPhase::CompressedLiquid
    );
    assert!(
        state.balance.sensible_outlet_temperature_k > state.balance.water_outlet.temperature_k()
    );
    assert!(state.balance.water_pressure_flow_work_change_rate_w < -100_000.0);
    assert!(state.heat_bracket_width_w <= config.numerics.heat_tolerance(heat).unwrap());
    audit(&owner, &state);
    audit(&owner, &rated);
}

#[test]
fn feedwater_boiling_superheat_preserve_mass_rate_scaling_and_independent_datums() {
    let chart = water();
    let inlet = at(&chart, 300.0);
    let target = at(&chart, 600.0);
    let (heat, ua) = required_ua(inlet, target, 1.0, |q| 900.0 + q / (4.0 * 1500.0), &[]);
    let mut outlets = vec![];
    for (scale, datum) in [(1.0, 300.0), (2.0, 700.0)] {
        let owner = CountercurrentWaterHeatExchanger::new(
            cp_stream(900.0, 4.0 * scale, 1500.0, datum),
            inlet,
            scale,
            ua * scale,
        )
        .unwrap();
        let state = with_cx(|cx| {
            owner
                .solve(cx, CountercurrentWaterExchangerConfig::default())
                .unwrap()
        });
        close(
            state.balance.heat_rate_water_to_sensible_w / scale,
            heat,
            0.05,
        );
        close(state.balance.water_outlet.temperature_k(), 600.0, 3e-5);
        assert_eq!(
            state.balance.water_outlet.phase(),
            WaterPhase::SuperheatedVapor
        );
        assert!(
            state.balance.sensible_outlet_temperature_k
                < state.balance.water_outlet.temperature_k()
        );
        assert!(state.balance.water_pressure_flow_work_change_rate_w > 100_000.0 * scale);
        audit(&owner, &state);
        outlets.push(state.balance);
    }
    close(
        outlets[0].water_outlet.specific_enthalpy_j_kg(),
        outlets[1].water_outlet.specific_enthalpy_j_kg(),
        1e-3,
    );
    close(
        outlets[0].sensible_outlet_temperature_k,
        outlets[1].sensible_outlet_temperature_k,
        1e-7,
    );
}

#[test]
fn pure_plateau_countercurrent_matches_the_original_exact_condenser_and_boiler() {
    let chart = water();
    for (tin, cp, mass, quality, water_mass, ua) in [
        (300.0, 4200.0, 1.5, 0.8, 0.4, 4000.0),
        (550.0, 1100.0, 2.0, 0.1, 0.25, 1800.0),
    ] {
        let inlet = chart.state_at_quality(quality, || false).unwrap();
        let sensible = cp_stream(tin, mass, cp, 350.0);
        let exact =
            WaterPhaseHeatExchanger::new(sensible.clone(), inlet.into(), water_mass, ua).unwrap();
        let expected = with_cx(|cx| {
            exact
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap()
        });
        let owner = CountercurrentWaterHeatExchanger::new(sensible, inlet, water_mass, ua).unwrap();
        let state = with_cx(|cx| {
            owner
                .solve(cx, CountercurrentWaterExchangerConfig::default())
                .unwrap()
        });
        close(
            state.balance.heat_rate_water_to_sensible_w,
            expected.heat_rate_water_to_sensible_w,
            0.01,
        );
        close(
            state.balance.water_outlet.vapor_quality().unwrap(),
            expected.water_outlet.vapor_quality().unwrap(),
            2e-8,
        );
        assert_eq!(
            state.balance.water_outlet.temperature_k(),
            inlet.temperature_k()
        );
        audit(&owner, &state);
    }
}

#[test]
fn nonlinear_original_cp_crosses_knots_with_independent_source_temperature_domains() {
    let chart = water();
    let inlet = at(&chart, 700.0);
    let target = at(&chart, 350.0);
    let knots = [
        (250.0, 3800.0),
        (310.0, 4400.0),
        (330.0, 3600.0),
        (400.0, 5000.0),
    ];
    let segments = [
        (300.0, 310.0, 4300.0, 10.0),
        (310.0, 330.0, 4400.0, -40.0),
        (330.0, 400.0, 3600.0, 20.0),
    ];
    // Independent analytic integral and stable quadratic inverse on authored
    // affine Cp segments. The 700 K WATER inlet is outside Cp support.
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
        panic!("reference source exhausted");
    };
    let (heat, ua) = required_ua(
        inlet,
        target,
        1.0,
        temperature,
        &[20.0 * 43_500.0, 20.0 * 123_500.0],
    );
    let owner = CountercurrentWaterHeatExchanger::new(
        stream(300.0, 20.0, source(&knots, 350.0, 2e5)),
        inlet,
        1.0,
        ua,
    )
    .unwrap();
    let state = with_cx(|cx| {
        owner
            .solve(cx, CountercurrentWaterExchangerConfig::default())
            .unwrap()
    });
    close(state.balance.heat_rate_water_to_sensible_w, heat, 0.05);
    close(
        state.balance.sensible_outlet_temperature_k,
        temperature(heat),
        1e-6,
    );
    close(state.balance.water_outlet.temperature_k(), 350.0, 2e-5);
    assert!(state.balance.sensible_outlet_temperature_k > 330.0);
    assert!(inlet.temperature_k() > owner.sensible_stream().caloric_law().temperature_span_k()[1]);
    audit(&owner, &state);
}

#[test]
fn original_interior_pinch_refuses_even_when_both_endpoint_gaps_are_positive() {
    let chart = water();
    let inlet = at(&chart, 700.0);
    let target = at(&chart, 350.0);
    let heat = inlet.specific_enthalpy_j_kg() - target.specific_enthalpy_j_kg();
    let law = source(
        &[
            (250.0, 9000.0),
            (350.0, 11000.0),
            (450.0, 9500.0),
            (650.0, 10500.0),
            (800.0, 10000.0),
        ],
        350.0,
        2e5,
    );
    let sensible_out = law.temperature_after_specific_heat(300.0, heat).unwrap();
    assert!(inlet.temperature_k() > sensible_out);
    assert!(target.temperature_k() > 300.0);
    let owner =
        CountercurrentWaterHeatExchanger::new(stream(300.0, 1.0, law), inlet, 1.0, 10_000.0)
            .unwrap();
    let config = CountercurrentWaterExchangerConfig::default();
    with_cx(|cx| {
        assert_eq!(
            owner.rate_duty(cx, heat, config).unwrap_err(),
            CountercurrentWaterExchangerError::Caloric(ExchangerError::Pinch)
        )
    });
    // A smaller realizable duty is still found: the rejected interior pinch
    // is a bracket restriction, never an accepted clipped phase inventory.
    let state = with_cx(|cx| owner.solve(cx, config).unwrap());
    assert!(state.balance.heat_rate_water_to_sensible_w < heat);
    audit(&owner, &state);
}

#[test]
fn source_bounds_work_limits_zero_exchange_and_cancellation_are_transactional() {
    let chart = water();
    let inlet = at(&chart, 700.0);
    let config = CountercurrentWaterExchangerConfig::default();
    for (tin, ua) in [(300.0, 0.0), (700.0, 1000.0)] {
        let owner = CountercurrentWaterHeatExchanger::new(
            cp_stream(tin, 1.0, 2000.0, 300.0),
            inlet,
            1.0,
            ua,
        )
        .unwrap();
        let state = with_cx(|cx| owner.solve(cx, config).unwrap());
        assert_eq!(state.balance.water_outlet, inlet);
        assert_eq!(state.balance.sensible_outlet_temperature_k, tin);
        assert_eq!(state.balance.heat_rate_water_to_sensible_w, 0.0);
        assert_eq!(state.water_source_queries, 0);
        audit(&owner, &state);
    }
    let limited = CountercurrentWaterHeatExchanger::new(
        stream(
            300.0,
            1.0,
            source(&[(290.0, 2000.0), (310.0, 2000.0)], 300.0, 1e5),
        ),
        inlet,
        1.0,
        1000.0,
    )
    .unwrap();
    with_cx(|cx| {
        assert_eq!(
            limited.solve(cx, config).unwrap_err(),
            CountercurrentWaterExchangerError::Caloric(ExchangerError::OutsideSourceSupport)
        )
    });
    let owner = CountercurrentWaterHeatExchanger::new(
        cp_stream(300.0, 5.0, 4200.0, 350.0),
        inlet,
        1.0,
        10_000.0,
    )
    .unwrap();
    with_cx(|cx| {
        assert_eq!(
            owner
                .solve(
                    cx,
                    CountercurrentWaterExchangerConfig {
                        maximum_water_source_queries: 4,
                        ..config
                    }
                )
                .unwrap_err(),
            CountercurrentWaterExchangerError::Caloric(ExchangerError::WorkLimit)
        )
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert_eq!(
            owner.solve(cx, config).unwrap_err(),
            CountercurrentWaterExchangerError::Interrupted
        )
    });
    assert_eq!(owner.water_inlet(), inlet);
    assert_eq!(limited.water_inlet(), inlet);
}
