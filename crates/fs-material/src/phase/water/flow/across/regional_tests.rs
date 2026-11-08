use super::*;
use crate::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls, RegionalWaterState,
};
use crate::phase::water::{WaterChartControls, WaterEntropySolve};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::regions::If97FluidRegion;

fn controls(minimum: f64, maximum: f64) -> RegionalWaterControls {
    RegionalWaterControls {
        minimum_temperature_k: minimum,
        maximum_temperature_k: maximum,
        enthalpy_tolerance_j_kg: 1e-6,
        maximum_iterations: 96,
        maximum_saturation_entropy_defect_j_kg_k: 0.1,
        pressure: DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
        junctions: RegionalJunctionControls {
            maximum_enthalpy_jump_j_kg: 200.0,
            maximum_entropy_jump_j_kg_k: 0.25,
            maximum_volume_jump_m3_kg: 1e-5,
            maximum_temperature_gap_k: 1e-8,
            maximum_iterations: 96,
        },
    }
}
fn chart(p: f64) -> RegionalWater {
    RegionalWater::new(Pressure::new(p), controls(273.16, 873.15), || false).unwrap()
}
fn shaft(efficiency: f64) -> WaterShaftSolve {
    WaterShaftSolve {
        isentropic_efficiency: efficiency,
        entropy: WaterEntropySolve {
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_iterations: 96,
        },
        balance: WaterFlowSolve {
            energy_tolerance_j_kg: 2e-5,
            entropy_tolerance_j_kg_k: 2e-8,
        },
    }
}
fn regional(state: WaterFlowState) -> RegionalWaterState {
    match state {
        WaterFlowState::Regional(state) => state,
        _ => panic!("regional pressure device changed source owner"),
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, tolerance {tolerance:.4e}"
    );
}

#[test]
fn uniform_regional_ingress_preserves_original_plateau_and_inverse() {
    let source = chart(20e6);
    let uniform = WaterFlowChart::from(source);
    let native = source.state_at_quality(0.375, || false).unwrap();
    assert_eq!(
        uniform.state_at_quality(0.375, || false).unwrap(),
        native.into()
    );
    let h = native.specific_enthalpy_j_kg();
    let solved = source.state_at_specific_enthalpy(h, || false).unwrap();
    let work = (solved.iterations(), solved.dense_source_evaluations());
    assert_eq!(
        uniform.state_at_specific_enthalpy(h, || false).unwrap(),
        (solved.into_state().into(), work.0, work.1)
    );
    assert!(
        WaterFlowChart::from(chart(25e6))
            .state_at_quality(0.5, || false)
            .is_err()
    );
}

#[test]
fn complete_high_pressure_cycles_include_cold_feed_and_all_source_region_crossings() {
    // G3: these are complete four-device cycles, unlike the earlier hot-inlet
    // steam train. Every heat operation uses one chart across all source pieces.
    for (p, hot_t, maximum) in [
        (20e6, 823.15, 873.15),
        (25e6, 823.15, 873.15),
        (20e6, 1500.0, 2273.15),
    ] {
        let low = chart(10_000.0);
        let high =
            RegionalWater::new(Pressure::new(p), controls(273.16, maximum), || false).unwrap();
        let liquid = low.state_at_quality(0.0, || false).unwrap();
        let pump = WaterFlowState::from(liquid)
            .adiabatic_shaft_to(
                &high.into(),
                WaterShaftDevice::Compressor,
                shaft(0.8),
                || false,
            )
            .unwrap();
        let feed = regional(pump.outlet);
        assert!(feed.temperature_k() < 330.0);
        assert_eq!(
            feed.single_phase_source().unwrap().region(),
            If97FluidRegion::Region1
        );
        for (t, region) in [
            (600.0, If97FluidRegion::Region1),
            (630.0, If97FluidRegion::Region3),
            (823.15, If97FluidRegion::Region2),
        ] {
            let source = high
                .state_at_temperature(Temperature::new(t), || false)
                .unwrap()
                .single_phase_source()
                .unwrap();
            assert_eq!(source.region(), region);
        }
        let target = high
            .state_at_temperature(Temperature::new(hot_t), || false)
            .unwrap();
        if hot_t > 1073.15 {
            assert_eq!(
                target.single_phase_source().unwrap().region(),
                If97FluidRegion::Region5
            );
            assert_eq!(high.junctions().count(), 3);
        }
        let q_in = target.specific_enthalpy_j_kg() - feed.specific_enthalpy_j_kg();
        let boiler = high.advance_heat(feed, 1.0, q_in, 2e-5, || false).unwrap();
        near(boiler.state.temperature_k(), hot_t, 2e-8);
        near(
            boiler.sensible_enthalpy_change_j
                + boiler.latent_enthalpy_change_j
                + boiler.regional_enthalpy_matching_change_j,
            q_in,
            2e-5,
        );
        assert!(boiler.regional_enthalpy_matching_change_j.abs() > 1e-3);
        assert!(
            boiler.regional_enthalpy_matching_change_j.abs()
                < 200.0 * high.junctions().count() as f64
        );
        if p < 22.064e6 {
            assert!(high.saturation().is_some());
            assert!(boiler.latent_enthalpy_change_j > 500_000.0);
        } else {
            assert!(high.saturation().is_none());
            assert_eq!(boiler.latent_enthalpy_change_j, 0.0);
            assert_eq!(boiler.state.vapor_mass_fraction(), None);
            assert_eq!(boiler.state.vapor_quality(), None);
        }
        let turbine = WaterFlowState::from(boiler.state)
            .adiabatic_shaft_to(&low.into(), WaterShaftDevice::Expander, shaft(0.85), || {
                false
            })
            .unwrap();
        let exhaust = regional(turbine.outlet);
        if hot_t > 1073.15 {
            assert_eq!(exhaust.vapor_mass_fraction(), Some(1.0));
            assert_eq!(exhaust.vapor_quality(), None);
        } else {
            assert!(exhaust.vapor_quality().is_some_and(|x| x > 0.0 && x < 1.0));
        }
        let q_out = liquid.specific_enthalpy_j_kg() - exhaust.specific_enthalpy_j_kg();
        let condenser = low
            .advance_heat(exhaust, 1.0, q_out, 2e-5, || false)
            .unwrap();
        assert_eq!(condenser.state, liquid);
        let net_work = pump.balance.shaft_work_j_kg + turbine.balance.shaft_work_j_kg;
        assert!(pump.balance.shaft_work_j_kg < 0.0);
        assert!(turbine.balance.shaft_work_j_kg > 0.0);
        assert!(net_work > 0.0 && net_work < q_in);
        assert!(q_out < 0.0);
        near(q_in + q_out, net_work, 8e-5);
        assert!(pump.balance.source_entropy_change_j_kg_k > 0.0);
        assert!(turbine.balance.source_entropy_change_j_kg_k > 0.0);

        // Independent original endpoint accounting includes p_out*v_out -
        // p_in*v_in exactly once, with no regional defect added as heat again.
        let mut totals = [0.0; 4];
        for (a, b, heat, work) in [
            (liquid, feed, 0.0, pump.balance.shaft_work_j_kg),
            (feed, boiler.state, q_in, 0.0),
            (boiler.state, exhaust, 0.0, turbine.balance.shaft_work_j_kg),
            (exhaust, condenser.state, q_out, 0.0),
        ] {
            let dh = b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg();
            let du = b.specific_internal_energy_j_kg() - a.specific_internal_energy_j_kg();
            let dpv = b.pressure_pa() * b.specific_volume_m3_kg()
                - a.pressure_pa() * a.specific_volume_m3_kg();
            let ds = b.specific_entropy_j_kg_k() - a.specific_entropy_j_kg_k();
            near(dh + work, heat, 2e-5);
            near(du + dpv + work, heat, 2e-5);
            for (total, change) in totals.iter_mut().zip([dh, du, dpv, ds]) {
                *total += change;
            }
        }
        for total in totals {
            near(total, 0.0, 8e-5);
        }
    }
}

#[test]
fn high_temperature_throttle_and_turbine_recover_independent_gibbs_balances() {
    let high = RegionalWater::new(Pressure::new(30e6), controls(900.0, 2273.15), || false).unwrap();
    let low = RegionalWater::new(Pressure::new(1e6), controls(600.0, 2273.15), || false).unwrap();
    let inlet = WaterFlowState::from(
        high.state_at_temperature(Temperature::new(1500.0), || false)
            .unwrap(),
    );
    // Independent 65-digit evaluation of original IAPWS R7-97(2012) Gibbs
    // equations, with scalar h/s roots. The inlet is Table 42's 30 MPa state.
    near(inlet.specific_enthalpy_j_kg(), 5_167_235.140089517, 1e-6);
    near(inlet.specific_entropy_j_kg_k(), 7729.701326182763, 1e-8);
    let check_first_law = |outlet: WaterFlowState, work| {
        near(
            outlet.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg() + work,
            0.0,
            2e-5,
        );
        near(
            outlet.specific_internal_energy_j_kg() - inlet.specific_internal_energy_j_kg()
                + outlet.pressure_pa() * outlet.specific_volume_m3_kg()
                - inlet.pressure_pa() * inlet.specific_volume_m3_kg()
                + work,
            0.0,
            2e-5,
        );
        assert!(outlet.specific_entropy_j_kg_k() - inlet.specific_entropy_j_kg_k() >= -2e-8);
    };
    let throttle = inlet
        .throttle_to(&low.into(), shaft(1.0).balance, || false)
        .unwrap();
    near(throttle.outlet.temperature_k(), 1480.238263631873, 2e-8);
    near(
        throttle.outlet.specific_entropy_j_kg_k(),
        9298.941428651517,
        2e-8,
    );
    assert_eq!(
        regional(throttle.outlet)
            .single_phase_source()
            .unwrap()
            .region(),
        If97FluidRegion::Region5
    );
    assert!(throttle.balance.source_entropy_change_j_kg_k > 0.0);
    check_first_law(throttle.outlet, 0.0);
    for (eta, expected_t, expected_work, expected_s) in [
        (
            1.0,
            761.0112142691587,
            1_714_514.0425254162,
            7729.701326182763,
        ),
        (
            0.85,
            878.2454305113812,
            1_457_336.9361466037,
            8043.898627729243,
        ),
    ] {
        let turbine = inlet
            .adiabatic_shaft_to(&low.into(), WaterShaftDevice::Expander, shaft(eta), || {
                false
            })
            .unwrap();
        near(turbine.outlet.temperature_k(), expected_t, 2e-8);
        near(turbine.outlet.specific_entropy_j_kg_k(), expected_s, 2e-8);
        near(turbine.balance.shaft_work_j_kg, expected_work, 2e-5);
        near(
            turbine.isentropic_outlet.temperature_k(),
            761.0112142691587,
            2e-8,
        );
        assert_eq!(
            regional(turbine.outlet)
                .single_phase_source()
                .unwrap()
                .region(),
            If97FluidRegion::Region2
        );
        assert_eq!(turbine.dense_source_evaluations, 0);
        check_first_law(turbine.outlet, turbine.balance.shaft_work_j_kg);
        if eta < 1.0 {
            assert!(turbine.balance.source_entropy_change_j_kg_k > 0.0);
        }
    }
}

#[test]
fn regional_devices_recover_the_existing_four_mpa_rankine_reference_without_datum_changes() {
    let old_chart = |p| {
        IsobaricWater::new(
            Pressure::new(p),
            WaterChartControls {
                minimum_temperature_k: 273.16,
                maximum_temperature_k: 873.15,
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
    };
    let (old_low, old_high) = (old_chart(10_000.0), old_chart(4e6));
    let (low, high) = (chart(10_000.0), chart(4e6));
    let old_liquid = old_low
        .state_at_specific_enthalpy(
            old_low.saturation().liquid().specific_enthalpy().value(),
            || false,
        )
        .unwrap()
        .into_state();
    let liquid = low.state_at_quality(0.0, || false).unwrap();
    assert_eq!(
        old_liquid.specific_enthalpy_j_kg(),
        liquid.specific_enthalpy_j_kg()
    );
    assert_eq!(
        old_liquid.specific_entropy_j_kg_k(),
        liquid.specific_entropy_j_kg_k()
    );
    let old_pump = old_liquid
        .adiabatic_shaft_to(&old_high, WaterShaftDevice::Compressor, shaft(0.8), || {
            false
        })
        .unwrap();
    let pump = WaterFlowState::from(liquid)
        .adiabatic_shaft_to(
            &high.into(),
            WaterShaftDevice::Compressor,
            shaft(0.8),
            || false,
        )
        .unwrap();
    let old_steam = old_high
        .state_at_temperature(Temperature::new(773.15), || false)
        .unwrap();
    let steam = high
        .state_at_temperature(Temperature::new(773.15), || false)
        .unwrap();
    assert_eq!(
        steam.specific_enthalpy_j_kg(),
        old_steam.specific_enthalpy_j_kg()
    );
    assert_eq!(
        steam.specific_entropy_j_kg_k(),
        old_steam.specific_entropy_j_kg_k()
    );
    let old_turbine = old_steam
        .adiabatic_shaft_to(&old_low, WaterShaftDevice::Expander, shaft(0.85), || false)
        .unwrap();
    let turbine = WaterFlowState::from(steam)
        .adiabatic_shaft_to(&low.into(), WaterShaftDevice::Expander, shaft(0.85), || {
            false
        })
        .unwrap();
    near(
        pump.balance.shaft_work_j_kg,
        old_pump.balance.shaft_work_j_kg,
        2e-6,
    );
    near(
        turbine.balance.shaft_work_j_kg,
        old_turbine.balance.shaft_work_j_kg,
        2e-6,
    );
    near(
        turbine.outlet.vapor_quality().unwrap(),
        old_turbine.outlet.vapor_quality().unwrap(),
        2e-12,
    );
    near(
        pump.outlet.temperature_k(),
        old_pump.outlet.temperature_k(),
        2e-8,
    );
}

#[test]
fn regional_endpoint_domain_and_cancellation_refusals_leave_the_feed_state_intact() {
    let low = chart(10_000.0);
    let liquid = low.state_at_quality(0.0, || false).unwrap();
    let inlet = WaterFlowState::from(liquid);
    let high = chart(20e6);
    assert!(matches!(
        inlet.adiabatic_shaft_to(
            &high.into(),
            WaterShaftDevice::Compressor,
            shaft(0.8),
            || true
        ),
        Err(AcrossWaterFlowError::Flow(WaterFlowError::Cancelled))
    ));
    let hot_only =
        RegionalWater::new(Pressure::new(20e6), controls(600.0, 873.15), || false).unwrap();
    assert!(matches!(
        inlet.adiabatic_shaft_to(
            &hot_only.into(),
            WaterShaftDevice::Compressor,
            shaft(0.8),
            || false
        ),
        Err(AcrossWaterFlowError::Regional(_))
    ));
    assert_eq!(regional(inlet), liquid);
    assert_eq!(inlet.pressure_pa(), 10_000.0);
    assert_eq!(inlet.vapor_quality(), Some(0.0));
}
