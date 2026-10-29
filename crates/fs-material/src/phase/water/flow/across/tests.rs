use super::*;
use crate::phase::water::pressure::DenseIsobaricControls;
use crate::phase::water::{WaterChartControls, WaterPhase};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;

fn gibbs(p: f64) -> IsobaricWater {
    gibbs_window(p, 273.16, 1073.15)
}

fn gibbs_window(p: f64, minimum: f64, maximum: f64) -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(p),
        WaterChartControls {
            minimum_temperature_k: minimum,
            maximum_temperature_k: maximum,
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
}

fn dense_controls(minimum: f64, maximum: f64) -> DenseIsobaricControls {
    DenseIsobaricControls {
        minimum_temperature_k: minimum,
        maximum_temperature_k: maximum,
        enthalpy_tolerance_j_kg: 1e-5,
        maximum_iterations: 96,
        maximum_entropy_seam_j_kg_k: 0.1,
        pressure: DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-8,
            maximum_iterations: 128,
        },
    }
}

fn dense(p: f64, minimum: f64, maximum: f64) -> DenseIsobaricWater {
    DenseIsobaricWater::new(Pressure::new(p), dense_controls(minimum, maximum), || false).unwrap()
}

fn balance() -> WaterFlowSolve {
    WaterFlowSolve {
        energy_tolerance_j_kg: 2e-5,
        entropy_tolerance_j_kg_k: 2e-8,
    }
}

fn shaft(efficiency: f64) -> WaterShaftSolve {
    WaterShaftSolve {
        isentropic_efficiency: efficiency,
        entropy: WaterEntropySolve {
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_iterations: 96,
        },
        balance: balance(),
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; allowance {tolerance:.3e}"
    );
}

#[test]
fn dense_liquid_throttles_into_a_gibbs_two_phase_outlet_with_flow_work_once() {
    let inlet = WaterFlowState::from(
        dense(20e6, 625.0, 645.0)
            .state_at_temperature(Temperature::new(630.0), || false)
            .unwrap(),
    );
    let low = gibbs(1e6);
    let step = inlet.throttle_to(&low.into(), balance(), || false).unwrap();
    let WaterFlowState::Gibbs(outlet) = step.outlet else {
        panic!("low-pressure outlet must retain its Gibbs source");
    };
    assert_eq!(outlet.phase(), WaterPhase::LiquidVapor);
    assert!(outlet.vapor_mass_fraction() > 0.3 && outlet.vapor_mass_fraction() < 0.6);
    assert_eq!(step.balance.shaft_work_j_kg, 0.0);
    near(
        step.outlet.specific_enthalpy_j_kg(),
        inlet.specific_enthalpy_j_kg(),
        balance().energy_tolerance_j_kg,
    );
    near(
        step.balance.first_law_residual_j_kg,
        0.0,
        balance().energy_tolerance_j_kg,
    );
    assert!(step.balance.source_entropy_change_j_kg_k > 0.0);
    for defect in step.pressure_closure_defects_j_kg {
        assert!(defect.abs() <= balance().energy_tolerance_j_kg);
    }
    let incorrect_flow_work =
        outlet.pressure_pa() * (outlet.specific_volume_m3_kg() - inlet.specific_volume_m3_kg());
    assert!((incorrect_flow_work - step.balance.flow_work_change_j_kg).abs() > 10_000.0);
    near(
        step.balance.internal_energy_change_j_kg + step.balance.flow_work_change_j_kg,
        0.0,
        balance().energy_tolerance_j_kg,
    );
}

#[test]
fn supercritical_expansion_uses_raw_entropy_and_retains_efficiency_loss() {
    let inlet = WaterFlowState::from(
        dense(25e6, 630.0, 675.0)
            .state_at_temperature(Temperature::new(670.0), || false)
            .unwrap(),
    );
    assert_eq!(inlet.vapor_mass_fraction(), None);
    let low = WaterFlowChart::from(gibbs(10_000.0));
    let ideal = inlet
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(1.0), || false)
        .unwrap();
    let actual = inlet
        .adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.85), || false)
        .unwrap();
    assert_eq!(ideal.outlet, ideal.isentropic_outlet);
    assert_eq!(ideal.enthalpy_iterations, 0);
    near(ideal.isentropic_entropy_residual_j_kg_k, 0.0, 1e-8);
    assert!(actual.balance.shaft_work_j_kg > 0.0);
    near(
        actual.balance.shaft_work_j_kg,
        0.85 * ideal.balance.shaft_work_j_kg,
        1e-8,
    );
    assert!(actual.balance.source_entropy_change_j_kg_k > 0.0);
    assert!(actual.outlet.vapor_quality().unwrap() > ideal.outlet.vapor_quality().unwrap());
    near(
        actual.balance.first_law_residual_j_kg,
        0.0,
        balance().energy_tolerance_j_kg,
    );
    let WaterFlowState::Gibbs(reference) = ideal.isentropic_outlet else {
        panic!("wet reference should use original Gibbs chart");
    };
    // The thermal primitive carries a pressure-specific seam correction. Using
    // it for an isentropic device would change the actual original entropy.
    assert!(reference.thermal_entropy_correction_j_kg_k().abs() > 1e-8);
    near(
        reference.specific_entropy_j_kg_k(),
        inlet.specific_entropy_j_kg_k(),
        1e-8,
    );
}

#[test]
fn hot_liquid_compression_crosses_into_dense_water_and_reverses_at_unit_efficiency() {
    let low = gibbs_window(16e6, 615.0, 626.0);
    let inlet = WaterFlowState::from(
        low.state_at_temperature(Temperature::new(620.0), || false)
            .unwrap(),
    );
    let high = WaterFlowChart::from(dense(50e6, 623.15, 700.0));
    let ideal = inlet
        .adiabatic_shaft_to(&high, WaterShaftDevice::Compressor, shaft(1.0), || false)
        .unwrap();
    assert!(matches!(ideal.outlet, WaterFlowState::Dense(_)));
    assert!(ideal.outlet.temperature_k() > 623.15);
    assert!(ideal.balance.shaft_work_j_kg < 0.0);
    assert!(ideal.dense_source_evaluations > 0);
    near(ideal.isentropic_entropy_residual_j_kg_k, 0.0, 1e-8);
    near(
        ideal.balance.first_law_residual_j_kg,
        0.0,
        balance().energy_tolerance_j_kg,
    );
    let lossy = inlet
        .adiabatic_shaft_to(&high, WaterShaftDevice::Compressor, shaft(0.8), || false)
        .unwrap();
    near(
        lossy.balance.shaft_work_j_kg,
        ideal.balance.shaft_work_j_kg / 0.8,
        1e-8,
    );
    assert!(lossy.balance.source_entropy_change_j_kg_k > 0.0);
    assert!(lossy.outlet.temperature_k() > ideal.outlet.temperature_k());
    let returned = ideal
        .outlet
        .adiabatic_shaft_to(&low.into(), WaterShaftDevice::Expander, shaft(1.0), || {
            false
        })
        .unwrap();
    near(returned.outlet.temperature_k(), inlet.temperature_k(), 2e-8);
    near(
        ideal.balance.shaft_work_j_kg + returned.balance.shaft_work_j_kg,
        0.0,
        balance().energy_tolerance_j_kg,
    );
}

#[test]
fn cross_chart_devices_refuse_bad_direction_budgets_domains_and_cancellation() {
    let high = dense(20e6, 625.0, 645.0);
    let inlet = WaterFlowState::from(
        high.state_at_temperature(Temperature::new(630.0), || false)
            .unwrap(),
    );
    let saved = inlet;
    let high_chart = WaterFlowChart::from(high);
    let low = WaterFlowChart::from(gibbs(1e6));
    assert!(matches!(
        inlet.throttle_to(&high_chart, balance(), || false),
        Err(AcrossWaterFlowError::Flow(
            WaterFlowError::PressureDirection
        ))
    ));
    assert!(matches!(
        inlet.adiabatic_shaft_to(&low, WaterShaftDevice::Expander, shaft(0.0), || false),
        Err(AcrossWaterFlowError::Flow(WaterFlowError::Invalid(_)))
    ));
    assert!(matches!(
        inlet.throttle_to(&low, balance(), || true),
        Err(AcrossWaterFlowError::Flow(WaterFlowError::Cancelled))
    ));
    let mut polls = 0;
    assert!(matches!(
        inlet.throttle_to(&low, balance(), || {
            polls += 1;
            polls >= 3
        }),
        Err(AcrossWaterFlowError::Flow(WaterFlowError::Cancelled))
    ));
    assert_eq!(inlet, saved);
    let mut impossible = balance();
    impossible.energy_tolerance_j_kg = 0.0;
    assert!(inlet.throttle_to(&low, impossible, || false).is_err());

    let cold = WaterFlowState::from(
        gibbs(100_000.0)
            .state_at_temperature(Temperature::new(300.0), || false)
            .unwrap(),
    );
    assert!(matches!(
        cold.adiabatic_shaft_to(
            &high_chart,
            WaterShaftDevice::Compressor,
            shaft(0.8),
            || false
        ),
        Err(AcrossWaterFlowError::Dense(_))
    ));

    let hot = WaterFlowState::from(
        gibbs_window(16e6, 615.0, 626.0)
            .state_at_temperature(Temperature::new(620.0), || false)
            .unwrap(),
    );
    let mut one_iteration = dense_controls(623.15, 700.0);
    one_iteration.maximum_iterations = 1;
    let limited = WaterFlowChart::from(
        DenseIsobaricWater::new(Pressure::new(50e6), one_iteration, || false).unwrap(),
    );
    assert!(matches!(
        hot.adiabatic_shaft_to(&limited, WaterShaftDevice::Compressor, shaft(0.8), || false),
        Err(AcrossWaterFlowError::Dense(_))
    ));
}
