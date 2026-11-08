use super::super::phase::*;
use super::*;
use fs_material::phase::water::flow::across::{WaterFlowChart, WaterFlowState};
use fs_material::phase::water::pressure::{DenseIsobaricControls, DenseIsobaricWater};
use fs_material::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_material::phase::water::{IsobaricWater, WaterChartControls};
use fs_thermochem::water::if97::If97SaturationPair;
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn water() -> WaterFlowChart {
    IsobaricWater::new(
        Pressure::new(1e6),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 623.15,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )
    .unwrap()
    .into()
}

fn phase_exchanger(
    chart: WaterFlowChart,
    quality: f64,
    sensible: ExchangerStream,
    water_mass: f64,
    ua: f64,
) -> WaterPhaseHeatExchanger {
    WaterPhaseHeatExchanger::new(
        sensible,
        chart.state_at_quality(quality, || false).unwrap(),
        water_mass,
        ua,
    )
    .unwrap()
}

fn phase_audit(state: WaterPhaseExchangerState, exchanger: &WaterPhaseHeatExchanger) {
    let heat = state.heat_rate_water_to_sensible_w;
    let config = WaterPhaseExchangerConfig::default();
    let tolerance = config.absolute_heat_tolerance_w + config.relative_heat_tolerance * heat.abs();
    close(
        exchanger.sensible_stream().mass_flow_kg_per_s()
            * state.sensible_specific_enthalpy_change_j_kg,
        heat,
        4.0 * f64::EPSILON * heat.abs(),
    );
    close(
        exchanger.water_mass_flow_kg_per_s() * state.water_specific_enthalpy_change_j_kg,
        -heat,
        4.0 * f64::EPSILON * heat.abs(),
    );
    assert!(
        state.water_endpoint_energy_residual_w.abs()
            <= tolerance + state.water_endpoint_roundoff_allowance_w
    );
    assert!(
        state.water_first_law_residual_w.abs()
            <= tolerance
                + state.water_endpoint_roundoff_allowance_w
                + state.water_pressure_closure_allowance_w
    );
    assert!(
        state.sensible_endpoint_energy_residual_w.abs()
            <= tolerance + state.sensible_endpoint_roundoff_allowance_w
    );
    assert!(state.thermal_entropy_generation_w_per_k >= -state.entropy_roundoff_allowance_w_per_k);
    close(
        state.source_entropy_generation_w_per_k
            + state.water_thermal_entropy_correction_rate_w_per_k,
        state.thermal_entropy_generation_w_per_k,
        2e-12,
    );
    assert_eq!(
        state.water_outlet.temperature_k(),
        exchanger.water_inlet().temperature_k()
    );
    assert!(state.water_outlet.vapor_quality().is_some());
}

#[test]
fn condensation_and_boiling_match_independent_ntu_and_original_latent_heat() {
    let chart = water();
    let pair = If97SaturationPair::at_pressure(Pressure::new(1e6)).unwrap();
    let ts = pair.temperature().value();
    let latent =
        pair.vapor().specific_enthalpy().value() - pair.liquid().specific_enthalpy().value();
    for (tin, cp, mass, quality, water_mass, ua) in [
        (300.0, 4200.0, 1.5, 0.8, 0.4, 4000.0),
        (550.0, 1100.0, 2.0, 0.1, 0.25, 1800.0),
    ] {
        let exchanger = phase_exchanger(
            chart,
            quality,
            stream(tin, mass, source(&[(280.0, cp), (600.0, cp)], 350.0, 2e5)),
            water_mass,
            ua,
        );
        let state = with_cx(|cx| {
            exchanger
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap()
        });
        let fraction = -(-ua / (mass * cp)).exp_m1();
        let expected_heat = mass * cp * (ts - tin) * fraction;
        let expected_quality = quality - expected_heat / (water_mass * latent);
        close(state.heat_rate_water_to_sensible_w, expected_heat, 2e-8);
        close(
            state.sensible_outlet_temperature_k,
            tin + (ts - tin) * fraction,
            2e-11,
        );
        close(
            state.water_outlet.vapor_quality().unwrap(),
            expected_quality,
            2e-14,
        );
        let expected_internal = water_mass
            * (expected_quality - quality)
            * (pair.vapor().specific_internal_energy().value()
                - pair.liquid().specific_internal_energy().value());
        let expected_work = water_mass
            * 1e6
            * (expected_quality - quality)
            * (pair.vapor().specific_volume().value() - pair.liquid().specific_volume().value());
        close(
            state.water_internal_energy_change_rate_w,
            expected_internal,
            2e-8,
        );
        close(
            state.water_pressure_flow_work_change_rate_w,
            expected_work,
            2e-8,
        );
        assert!(state.water_pressure_flow_work_change_rate_w.abs() > 1000.0);
        close(expected_internal + expected_work, -expected_heat, 2e-8);
        phase_audit(state, &exchanger);
    }
}

#[test]
fn nonlinear_cp_crosses_original_knots_without_freezing_capacity_or_shifting_datum() {
    let chart = water();
    let ts = chart
        .state_at_quality(0.0, || false)
        .unwrap()
        .temperature_k();
    let knots = [
        (280.0, 1000.0),
        (330.0, 1500.0),
        (390.0, 900.0),
        (470.0, 2100.0),
        (550.0, 1800.0),
    ];
    let (tin, tout, mass) = (300.0, 420.0, 2.0);
    // Independent analytic integrals of each original affine segment, for a
    // manufactured outlet that crosses two genuine Cp slope changes.
    let mut exposure = 0.0;
    let mut heat_per_mass = 0.0;
    for segment in knots.windows(2) {
        let [(left, cp_left), (right, cp_right)] = [segment[0], segment[1]];
        let low = f64::max(left, tin);
        let high = f64::min(right, tout);
        if high <= low {
            continue;
        }
        let slope = (cp_right - cp_left) / (right - left);
        let intercept = cp_left - slope * left;
        exposure +=
            (intercept + slope * ts) * ((ts - low) / (ts - high)).ln() - slope * (high - low);
        heat_per_mass += (intercept + 0.5 * slope * (low + high)) * (high - low);
    }
    for datum in [300.0, 430.0] {
        let exchanger = phase_exchanger(
            chart,
            0.75,
            stream(tin, mass, source(&knots, datum, 2e5)),
            0.4,
            mass * exposure,
        );
        let state = with_cx(|cx| {
            exchanger
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap()
        });
        close(state.sensible_outlet_temperature_k, tout, 2e-11);
        close(
            state.heat_rate_water_to_sensible_w,
            mass * heat_per_mass,
            2e-8,
        );
        phase_audit(state, &exchanger);
    }
}

#[test]
fn original_dense_and_regional_water_keep_their_own_plateau_and_pressure_closure() {
    let pressure_controls = DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    };
    let dense = DenseIsobaricWater::new(
        Pressure::new(20e6),
        DenseIsobaricControls {
            minimum_temperature_k: 625.0,
            maximum_temperature_k: 645.0,
            enthalpy_tolerance_j_kg: 1e-5,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            pressure: pressure_controls,
        },
        || false,
    )
    .unwrap();
    let regional = RegionalWater::new(
        Pressure::new(20e6),
        RegionalWaterControls {
            minimum_temperature_k: 300.0,
            maximum_temperature_k: 873.15,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: pressure_controls,
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap();
    for chart in [WaterFlowChart::from(dense), WaterFlowChart::from(regional)] {
        let exchanger = phase_exchanger(
            chart,
            0.6,
            stream(
                500.0,
                1.0,
                source(&[(450.0, 2000.0), (700.0, 2000.0)], 500.0, 1e5),
            ),
            2.0,
            1200.0,
        );
        let state = with_cx(|cx| {
            exchanger
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap()
        });
        assert_eq!(state.water_outlet.law(), chart);
        assert!(state.water_outlet.vapor_quality().unwrap() < 0.6);
        assert_eq!(state.enthalpy_iterations, 0);
        assert_eq!(state.dense_source_evaluations, 0);
        assert!(
            state
                .water_pressure_closure_defects_j_kg
                .into_iter()
                .all(|d| d.abs() <= 1e-4)
        );
        assert!(
            state
                .water_pressure_closure_defects_j_kg
                .into_iter()
                .any(|d| d != 0.0)
        );
        with_cx(|cx| {
            assert!(matches!(
                exchanger.solve(
                    cx,
                    WaterPhaseExchangerConfig {
                        maximum_pressure_closure_defect_j_kg: 0.0,
                        ..WaterPhaseExchangerConfig::default()
                    }
                ),
                Err(WaterPhaseExchangerError::PressureClosure { .. })
            ));
        });
        phase_audit(state, &exchanger);
    }
}

#[test]
fn finite_latent_throughput_refuses_exhaustion_in_both_directions_even_below_h_resolution() {
    let chart = water();
    for (quality, tin, water_mass, ua) in [
        (0.01, 300.0, 0.1, 5000.0),
        (0.99, 600.0, 0.1, 5000.0),
        (0.0, 300.0, 1e8, 1e-9),
        (1.0, 600.0, 1e8, 1e-9),
    ] {
        let exchanger = phase_exchanger(
            chart,
            quality,
            stream(
                tin,
                1.0,
                source(&[(280.0, 2000.0), (650.0, 2000.0)], 300.0, 1e5),
            ),
            water_mass,
            ua,
        );
        with_cx(|cx| {
            assert_eq!(
                exchanger
                    .solve(cx, WaterPhaseExchangerConfig::default())
                    .unwrap_err(),
                WaterPhaseExchangerError::LatentThroughputExceeded
            )
        });
    }
    let WaterFlowChart::Gibbs(source_chart) = chart else {
        unreachable!()
    };
    let compressed: WaterFlowState = source_chart
        .state_at_temperature(Temperature::new(300.0), || false)
        .unwrap()
        .into();
    assert_eq!(
        WaterPhaseHeatExchanger::new(
            stream(300.0, 1.0, affine(2000.0, 0.0, 300.0, 1e5)),
            compressed,
            1.0,
            100.0
        )
        .unwrap_err(),
        WaterPhaseExchangerError::OutsideWaterPlateau
    );
}

#[test]
fn raw_source_entropy_defect_is_reported_and_requires_its_own_explicit_allowance() {
    let chart = water();
    let liquid = chart.state_at_quality(0.0, || false).unwrap();
    let vapor = chart.state_at_quality(1.0, || false).unwrap();
    let ts = liquid.temperature_k();
    let latent = vapor.specific_enthalpy_j_kg() - liquid.specific_enthalpy_j_kg();
    let seam = chart.saturation_entropy_defect_j_kg_k();
    // Choose the physical heat direction and a small temperature gap so the
    // signed source-fit mismatch dominates positive irreversible production.
    let gap = 0.02 * seam.abs() * ts * ts / latent;
    let tin = ts + gap.copysign(seam);
    let exchanger = phase_exchanger(
        chart,
        0.5,
        stream(tin, 1.0, affine(2000.0, 0.0, 300.0, 1e5)),
        1.0,
        1000.0,
    );
    with_cx(|cx| {
        assert!(matches!(
            exchanger.solve(cx, WaterPhaseExchangerConfig::default()),
            Err(WaterPhaseExchangerError::SourceEntropyDecrease { .. })
        ));
        let state = exchanger
            .solve(
                cx,
                WaterPhaseExchangerConfig {
                    source_entropy_tolerance_w_per_k: 1e-4,
                    ..WaterPhaseExchangerConfig::default()
                },
            )
            .unwrap();
        assert!(state.source_entropy_generation_w_per_k < 0.0);
        assert!(
            state.thermal_entropy_generation_w_per_k >= -state.entropy_roundoff_allowance_w_per_k
        );
        close(
            state.water_thermal_entropy_correction_rate_w_per_k,
            -state.heat_rate_water_to_sensible_w * seam / latent,
            2e-15,
        );
        phase_audit(state, &exchanger);
    });
}

#[test]
fn source_support_zero_transfer_invalid_controls_and_cancellation_are_transactional() {
    let chart = water();
    let limited = phase_exchanger(
        chart,
        0.9,
        stream(
            300.0,
            1.0,
            source(&[(280.0, 2000.0), (320.0, 2000.0)], 300.0, 1e5),
        ),
        10.0,
        1000.0,
    );
    with_cx(|cx| {
        assert!(matches!(
            limited.solve(cx, WaterPhaseExchangerConfig::default()),
            Err(WaterPhaseExchangerError::Caloric(ExchangerError::Source(_)))
        ))
    });
    let exchanger = phase_exchanger(
        chart,
        0.7,
        stream(300.0, 1.0, affine(2000.0, 0.0, 300.0, 1e5)),
        1.0,
        0.0,
    );
    assert!(
        WaterPhaseHeatExchanger::new(
            exchanger.sensible_stream().clone(),
            exchanger.water_inlet(),
            0.0,
            0.0,
        )
        .is_err()
    );
    with_cx(|cx| {
        let state = exchanger
            .solve(cx, WaterPhaseExchangerConfig::default())
            .unwrap();
        assert_eq!(state.water_outlet, exchanger.water_inlet());
        assert_eq!(state.sensible_outlet_temperature_k, 300.0);
        assert_eq!(state.heat_rate_water_to_sensible_w, 0.0);
        assert_eq!(state.source_entropy_generation_w_per_k, 0.0);
        phase_audit(state, &exchanger);
        assert!(
            exchanger
                .solve(
                    cx,
                    WaterPhaseExchangerConfig {
                        source_entropy_tolerance_w_per_k: -1.0,
                        ..WaterPhaseExchangerConfig::default()
                    }
                )
                .is_err()
        );
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert_eq!(
            exchanger
                .solve(cx, WaterPhaseExchangerConfig::default())
                .unwrap_err(),
            WaterPhaseExchangerError::Interrupted
        )
    });
    assert_eq!(
        with_cx(|cx| exchanger
            .solve(cx, WaterPhaseExchangerConfig::default())
            .unwrap())
        .water_outlet,
        exchanger.water_inlet()
    );
}
