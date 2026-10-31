use super::*;
use crate::phase::water::WaterChartControls;
use fs_qty::Pressure;
use fs_thermal::evaporation::{
    GasFilm, LiquidSurfaceState, SolverConfig, advance_surface, surface_flux,
};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

fn chart(pressure: f64) -> IsobaricWater {
    IsobaricWater::new(
        Pressure::new(pressure),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 450.0,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 65535,
            },
        },
        || false,
    )
    .unwrap()
}
fn water() -> WaterEvaporation {
    WaterEvaporation::try_new(chart(101_325.0), [280.0, 350.0]).unwrap()
}
fn film(pv: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 101_325.0,
        vapor_partial_pressure_pa: pv,
        temperature_k: 300.0,
        mass_transfer_coefficient_m_s: 0.005,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn g1_water_rate_uses_original_if97_liquid_and_vapor_calorics() {
    let water = water();
    for t in [280.0, 300.0, 350.0] {
        let p = water.properties(t).unwrap();
        let original_liquid = water
            .liquid_chart()
            .state_at_temperature(Temperature::new(t), || false)
            .unwrap();
        let original_pair = If97SaturationPair::at_temperature(Temperature::new(t)).unwrap();
        assert_eq!(p.saturation_pressure_pa, original_pair.pressure().value());
        assert_eq!(p.liquid_density_kg_m3, original_liquid.bulk_density_kg_m3());
        assert_eq!(
            p.latent_heat_j_kg,
            original_pair.vapor().specific_enthalpy().value()
                - original_liquid.specific_enthalpy_j_kg()
        );
        assert_eq!(
            water.vapor_specific_enthalpy_j_kg(t).unwrap(),
            original_pair.vapor().specific_enthalpy().value()
        );
        let h = original_liquid.specific_enthalpy_j_kg();
        near(
            water.temperature_at_specific_enthalpy_k(h).unwrap(),
            t,
            3e-8,
        );
    }
    // Independent published-scale checks: IAPWS/NIST 300 K saturation is
    // 3.5368 kPa with approximately 2437 kJ/kg vaporization enthalpy.
    let p = water.properties(300.0).unwrap();
    near(p.saturation_pressure_pa, 3536.8, 0.5);
    near(p.latent_heat_j_kg, 2_437_000.0, 1000.0);
    let flux = surface_flux(&water, 300.0, film(1000.0)).unwrap();
    let expected = 0.005 * 101_325.0 * 0.018_015_28 / (8.314_462_618_153_24 * 300.0)
        * ((101_325.0_f64 - 1000.0) / (101_325.0 - p.saturation_pressure_pa)).ln();
    near(flux.mass_flux_kg_m2_s, expected, 2e-18);
}

#[test]
fn g1_actual_subboiling_water_cools_and_condenses_on_its_original_thermal_chart() {
    let water = water();
    let initial = LiquidSurfaceState {
        mass_kg: 0.02,
        specific_enthalpy_j_kg: water.liquid_state(300.0).unwrap().specific_enthalpy_j_kg(),
        area_m2: 0.01,
    };
    let controls = SolverConfig {
        specific_enthalpy_bracket_j_kg: water.enthalpy_range_j_kg(),
        relative_tolerance: 1e-10,
        absolute_energy_tolerance_j: 1e-7,
        max_iterations: 128,
    };
    for vapor_pressure in [500.0, 5000.0] {
        let step = advance_surface(
            &water,
            film(vapor_pressure),
            initial,
            5.0,
            controls,
            |h| water.temperature_at_specific_enthalpy_k(h),
            || Ok(()),
        )
        .unwrap();
        assert!(step.energy_residual_j.abs() < 2e-7);
        near(
            step.liquid.mass_kg * step.liquid.specific_enthalpy_j_kg + step.vapor_energy_j,
            initial.mass_kg * initial.specific_enthalpy_j_kg,
            2e-7,
        );
        assert_eq!(
            initial.mass_kg - step.liquid.mass_kg,
            step.transferred_mass_kg
        );
        near(
            step.vapor_energy_j,
            step.transferred_mass_kg
                * water
                    .vapor_specific_enthalpy_j_kg(step.surface_temperature_k)
                    .unwrap(),
            1e-10,
        );
        if vapor_pressure < 3000.0 {
            assert!(step.surface_temperature_k < 300.0);
        } else {
            assert!(step.surface_temperature_k > 300.0);
        }
        near(
            ThermalLaw::temperature_k(&water, step.liquid.specific_enthalpy_j_kg).unwrap(),
            step.surface_temperature_k,
            1e-12,
        );
        assert!(
            ThermalLaw::specific_entropy_change_j_kg_k(
                &water,
                initial.specific_enthalpy_j_kg,
                step.liquid.specific_enthalpy_j_kg
            )
            .unwrap()
            .is_some()
        );
    }
}

#[test]
fn g1_water_analytic_saturation_tangent_matches_original_pressure_variation() {
    let water = water();
    let h = 0.001;
    for t in [285.0, 300.0, 340.0] {
        let left = surface_flux(&water, t - h, film(1000.0)).unwrap();
        let right = surface_flux(&water, t + h, film(1000.0)).unwrap();
        let fd = (right.mass_flux_kg_m2_s - left.mass_flux_kg_m2_s) / (2.0 * h);
        let slope = surface_flux(&water, t, film(1000.0))
            .unwrap()
            .mass_flux_derivative_kg_m2_s_k;
        near(slope, fd, slope.abs() * 3e-8);
    }
}

#[test]
fn g4_water_adapter_refuses_boiling_high_pressure_ice_and_caloric_extrapolation() {
    let liquid = chart(101_325.0);
    for range in [
        [273.0, 300.0],
        [280.0, 373.15],
        [300.0, 300.0],
        [f64::NAN, 330.0],
    ] {
        assert!(WaterEvaporation::try_new(liquid, range).is_err());
    }
    assert!(WaterEvaporation::try_new(chart(300_000.0), [280.0, 350.0]).is_err());
    let water = water();
    assert!(water.properties(279.0).is_err());
    assert!(
        water
            .temperature_at_specific_enthalpy_k(water.enthalpy_range_j_kg()[0] - 1.0)
            .is_err()
    );
    assert!(
        ThermalLaw::temperature_k(
            &water,
            liquid.saturation().vapor().specific_enthalpy().value()
        )
        .is_err()
    );
    let low_pressure = GasFilm {
        total_pressure_pa: 3000.0,
        ..film(1000.0)
    };
    assert!(matches!(
        surface_flux(&water, 300.0, low_pressure),
        Err(EvaporationError::Boiling { .. })
    ));
}
