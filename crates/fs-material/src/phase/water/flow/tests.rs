use super::*;
use crate::phase::water::{WaterChartControls, WaterPhase};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn chart(p: f64) -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(p),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 1073.15,
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
fn entropy() -> WaterEntropySolve {
    WaterEntropySolve {
        entropy_tolerance_j_kg_k: 1e-9,
        maximum_iterations: 96,
    }
}
fn accuracy() -> WaterFlowSolve {
    WaterFlowSolve {
        energy_tolerance_j_kg: 2e-7,
        entropy_tolerance_j_kg_k: 2e-9,
    }
}
fn shaft(efficiency: f64) -> WaterShaftSolve {
    WaterShaftSolve {
        isentropic_efficiency: efficiency,
        entropy: entropy(),
        balance: accuracy(),
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e} (budget {tolerance:.5e})"
    );
}
fn saturated_liquid(law: &IsobaricWater) -> WaterState {
    law.state_at_specific_enthalpy(
        law.saturation().liquid().specific_enthalpy().value(),
        || false,
    )
    .unwrap()
    .into_state()
}

#[test]
fn pressure_throttling_reaches_wet_and_superheated_original_source_outlets() {
    let hot = chart(1e6);
    let low = chart(101_325.0);
    let liquid = saturated_liquid(&hot);
    let wet = liquid.throttle_to(&low, accuracy(), || false).unwrap();
    assert_eq!(wet.inlet, liquid);
    assert_eq!(wet.outlet.phase(), WaterPhase::LiquidVapor);
    assert!(wet.outlet.vapor_mass_fraction() > 0.1);
    assert!(wet.outlet.vapor_mass_fraction() < 0.2);
    assert_eq!(wet.balance.shaft_work_j_kg, 0.0);
    near(
        wet.outlet.specific_enthalpy_j_kg(),
        liquid.specific_enthalpy_j_kg(),
        accuracy().energy_tolerance_j_kg,
    );
    assert!(wet.balance.source_entropy_change_j_kg_k > 0.0);
    assert!(wet.balance.flow_work_change_j_kg > 0.0);
    assert!(wet.balance.internal_energy_change_j_kg < 0.0);
    near(
        wet.balance.first_law_residual_j_kg,
        0.0,
        accuracy().energy_tolerance_j_kg,
    );

    let steam = hot
        .state_at_temperature(Temperature::new(800.0), || false)
        .unwrap();
    let expanded = steam.throttle_to(&low, accuracy(), || false).unwrap();
    assert_eq!(expanded.outlet.phase(), WaterPhase::SuperheatedVapor);
    assert_eq!(expanded.outlet.vapor_quality(), None);
    assert!(expanded.outlet.temperature_k() > low.saturation().temperature().value());
    assert!(expanded.balance.source_entropy_change_j_kg_k > 0.0);
    near(
        expanded.balance.enthalpy_residual_j_kg,
        0.0,
        accuracy().energy_tolerance_j_kg,
    );
    // Using p_out*Delta v would omit the actual inlet flow work.
    let wrong_work = expanded.outlet.pressure_pa()
        * (expanded.outlet.specific_volume_m3_kg() - steam.specific_volume_m3_kg());
    assert!((wrong_work - expanded.balance.flow_work_change_j_kg).abs() > 100_000.0);
    near(
        expanded.balance.source_entropy_change_j_kg_k,
        expanded.outlet.specific_entropy_j_kg_k() - steam.specific_entropy_j_kg_k(),
        0.0,
    );
}

#[test]
fn isentropic_liquid_pump_work_agrees_with_independent_pressure_volume_integral() {
    let inlet = chart(100_000.0)
        .state_at_temperature(Temperature::new(300.0), || false)
        .unwrap();
    let outlet = chart(4e6);
    let step = inlet
        .adiabatic_shaft_to(&outlet, WaterShaftDevice::Compressor, shaft(1.0), || false)
        .unwrap();
    assert_eq!(step.outlet.phase(), WaterPhase::CompressedLiquid);
    assert_eq!(step.outlet, step.isentropic_outlet);
    assert_eq!(step.enthalpy_iterations, 0);
    assert!(step.balance.shaft_work_j_kg < 0.0);
    assert!(step.outlet.temperature_k() > inlet.temperature_k());
    near(
        step.balance.source_entropy_change_j_kg_k,
        0.0,
        entropy().entropy_tolerance_j_kg_k,
    );
    // Maxwell relation on a liquid isentrope: dh = v dp. Composite Simpson
    // integration is independent of the device's endpoint enthalpy subtraction.
    let count = 16;
    let dp = (4e6 - 100_000.0) / f64::from(count);
    let mut integral = 0.0;
    for i in 0..=count {
        let pressure = 100_000.0 + f64::from(i) * dp;
        let state = chart(pressure)
            .state_at_specific_entropy(inlet.specific_entropy_j_kg_k(), entropy(), || false)
            .unwrap()
            .into_state();
        let weight = if i == 0 || i == count {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        integral += weight * state.specific_volume_m3_kg();
    }
    integral *= dp / 3.0;
    near(-step.balance.shaft_work_j_kg, integral, 2e-5);
    let lossy = inlet
        .adiabatic_shaft_to(&outlet, WaterShaftDevice::Compressor, shaft(0.8), || false)
        .unwrap();
    near(
        lossy.balance.shaft_work_j_kg,
        step.balance.shaft_work_j_kg / 0.8,
        1e-10,
    );
    assert!(lossy.balance.source_entropy_change_j_kg_k > 0.0);
    assert!(lossy.outlet.temperature_k() > step.outlet.temperature_k());
    near(
        lossy.balance.first_law_residual_j_kg,
        0.0,
        accuracy().energy_tolerance_j_kg,
    );
}

#[test]
fn superheated_steam_expansion_and_compression_obey_signed_efficiency_and_entropy() {
    let high = chart(4e6);
    let low = chart(10_000.0);
    let inlet = high
        .state_at_temperature(Temperature::new(773.15), || false)
        .unwrap();
    let ideal = inlet
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(1.0), || false)
        .unwrap();
    let turbine = inlet
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.85), || false)
        .unwrap();
    assert_eq!(ideal.outlet.phase(), WaterPhase::LiquidVapor);
    assert_eq!(turbine.outlet.phase(), WaterPhase::LiquidVapor);
    assert!(turbine.outlet.vapor_mass_fraction() > ideal.outlet.vapor_mass_fraction());
    assert!(turbine.balance.shaft_work_j_kg > 0.0);
    near(
        turbine.balance.shaft_work_j_kg,
        0.85 * ideal.balance.shaft_work_j_kg,
        1e-8,
    );
    assert!(turbine.balance.source_entropy_change_j_kg_k > 0.0);
    near(
        ideal.isentropic_entropy_residual_j_kg_k,
        0.0,
        entropy().entropy_tolerance_j_kg_k,
    );
    near(
        turbine.balance.first_law_residual_j_kg,
        0.0,
        accuracy().energy_tolerance_j_kg,
    );

    let steam = chart(100_000.0)
        .state_at_temperature(Temperature::new(500.0), || false)
        .unwrap();
    let compressor = steam
        .adiabatic_shaft_to(
            &chart(1e6),
            WaterShaftDevice::Compressor,
            shaft(0.8),
            || false,
        )
        .unwrap();
    assert_eq!(compressor.outlet.phase(), WaterPhase::SuperheatedVapor);
    assert!(compressor.outlet.temperature_k() > steam.temperature_k());
    assert!(compressor.balance.shaft_work_j_kg < 0.0);
    assert!(compressor.balance.source_entropy_change_j_kg_k > 0.0);
    near(
        compressor.balance.enthalpy_residual_j_kg,
        0.0,
        accuracy().energy_tolerance_j_kg,
    );
}

#[test]
fn complete_rankine_cycle_balances_heat_and_shaft_work_without_double_counting_flow_work() {
    let low = chart(10_000.0);
    let high = chart(4e6);
    let liquid = saturated_liquid(&low);
    let pump = liquid
        .adiabatic_shaft_to(&high, WaterShaftDevice::Compressor, shaft(0.8), || false)
        .unwrap();
    let steam = high
        .state_at_temperature(Temperature::new(773.15), || false)
        .unwrap();
    let turbine = steam
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.85), || false)
        .unwrap();
    let heat_in = steam.specific_enthalpy_j_kg() - pump.outlet.specific_enthalpy_j_kg();
    let heat_out = liquid.specific_enthalpy_j_kg() - turbine.outlet.specific_enthalpy_j_kg();
    let net_work = turbine.balance.shaft_work_j_kg + pump.balance.shaft_work_j_kg;
    assert!(heat_in > 0.0 && heat_out < 0.0 && net_work > 0.0);
    near(
        heat_in + heat_out,
        net_work,
        2.0 * accuracy().energy_tolerance_j_kg,
    );
    assert!((0.2..0.5).contains(&(net_work / heat_in)));
    let pv = |s: WaterState| s.pressure_pa() * s.specific_volume_m3_kg();
    let boiler_pv = pv(steam) - pv(pump.outlet);
    let condenser_pv = pv(liquid) - pv(turbine.outlet);
    near(
        pump.balance.flow_work_change_j_kg
            + boiler_pv
            + turbine.balance.flow_work_change_j_kg
            + condenser_pv,
        0.0,
        1e-9,
    );
    let boiler_u =
        steam.specific_internal_energy_j_kg() - pump.outlet.specific_internal_energy_j_kg();
    let condenser_u =
        liquid.specific_internal_energy_j_kg() - turbine.outlet.specific_internal_energy_j_kg();
    near(
        boiler_u + boiler_pv,
        heat_in,
        accuracy().energy_tolerance_j_kg,
    );
    near(
        condenser_u + condenser_pv,
        heat_out,
        accuracy().energy_tolerance_j_kg,
    );
    // Boiler/condenser p*Delta v is already flow work within h; adding it to
    // turbine/pump shaft output would give a materially wrong cycle efficiency.
    assert!((boiler_pv + condenser_pv).abs() > 10_000.0);
}

#[test]
fn devices_refuse_invalid_direction_efficiency_source_accuracy_and_late_cancellation() {
    let high = chart(4e6);
    let low = chart(10_000.0);
    let state = high
        .state_at_temperature(Temperature::new(773.15), || false)
        .unwrap();
    let saved = state;
    assert_eq!(
        state.throttle_to(&high, accuracy(), || false),
        Err(WaterFlowError::PressureDirection)
    );
    assert_eq!(
        state.adiabatic_shaft_to(&low, WaterShaftDevice::Compressor, shaft(0.8), || false),
        Err(WaterFlowError::PressureDirection)
    );
    for efficiency in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            state.adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(efficiency), || {
                false
            }),
            Err(WaterFlowError::Invalid(_))
        ));
    }
    assert!(matches!(
        state.adiabatic_shaft_to(
            &low,
            WaterShaftDevice::Expander,
            WaterShaftSolve {
                balance: WaterFlowSolve {
                    energy_tolerance_j_kg: 0.0,
                    ..accuracy()
                },
                ..shaft(0.83)
            },
            || false
        ),
        Err(WaterFlowError::EnergyBalance { .. })
    ));
    let mut count = 0;
    let accepted = state
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.85), || {
            count += 1;
            false
        })
        .unwrap();
    assert!(count > 2);
    let mut observed = 0;
    assert_eq!(
        state.adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.85), || {
            observed += 1;
            observed == count
        }),
        Err(WaterFlowError::Cancelled)
    );
    assert_eq!(state, saved);
    assert_eq!(accepted.inlet, saved);
    let cold = low
        .state_at_temperature(Temperature::new(280.0), || false)
        .unwrap();
    // A very inefficient pump asks for energy beyond the caller's source chart.
    assert!(matches!(
        cold.adiabatic_shaft_to(&high, WaterShaftDevice::Compressor, shaft(1e-6), || false),
        Err(WaterFlowError::Source(WaterError::OutsideDomain))
    ));
}
