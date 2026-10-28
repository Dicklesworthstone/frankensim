//! Routing regressions, not new experimental-accuracy claims. Both humidity
//! entry points ultimately use the shared gas-mixture constitutive kernel.

use fs_material::gas::{
    GasSpec, GasState, HumidityReference, saturation_pressure_water_pa,
    saturation_pressure_water_pa_with_reference,
};
use fs_qty::semantic::{Composition, CompositionBasis};

#[test]
fn humid_front_end_and_general_composition_produce_identical_states() {
    for (temperature, pressure, humidity) in [
        (253.15, 101_325.0, 0.2),
        (293.15, 101_325.0, 0.5),
        (310.15, 90_000.0, 0.8),
        (323.15, 200_000.0, 1.0),
    ] {
        let vapor_fraction = humidity * saturation_pressure_water_pa(temperature).unwrap() / pressure;
        let composition = Composition::new(
            CompositionBasis::MoleFraction,
            vec![1.0 - vapor_fraction, vapor_fraction],
        ).unwrap();
        let general = GasState::try_new_mixture(
            &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
            &composition,
            temperature,
            pressure,
            Some(1),
        ).unwrap();
        let humid = GasState::try_new_moist_air(temperature, pressure, humidity).unwrap();
        assert_eq!(general, humid, "mixture routing at T={temperature}, p={pressure}, RH={humidity}");
    }
}

#[test]
fn zero_humidity_keeps_the_original_pure_gas_path_and_domain() {
    for (temperature, pressure) in [
        (50.0, 1000.0),
        (293.15, 101_325.0),
        (2000.0, 1.0e7),
    ] {
        let humid = GasState::try_new_moist_air(temperature, pressure, 0.0).unwrap();
        let pure = GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap();
        assert_eq!(humid, pure);
    }
}

#[test]
fn humid_front_end_retains_its_phase_and_fit_refusals() {
    assert!(GasState::try_new_moist_air(293.15, 101_325.0, 1.1).is_err());
    assert!(GasState::try_new_moist_air(240.0, 101_325.0, 0.5).is_err());
    assert!(GasState::try_new_moist_air(293.15, 100.0, 1.0).is_err());
}

#[test]
fn explicit_liquid_reference_preserves_existing_results_and_dry_domain() {
    for (temperature, pressure, humidity) in [
        (253.15, 101_325.0, 0.2),
        (293.15, 101_325.0, 0.5),
        (310.15, 90_000.0, 0.8),
        (323.15, 200_000.0, 1.0),
    ] {
        let original = GasState::try_new_moist_air(temperature, pressure, humidity).unwrap();
        let explicit = GasState::try_new_moist_air_with_reference(
            temperature,
            pressure,
            humidity,
            HumidityReference::LiquidWater,
        )
        .unwrap();
        assert_eq!(original, explicit);
    }
    for reference in [HumidityReference::LiquidWater, HumidityReference::Ice] {
        for temperature in [50.0, 230.0, 293.15, 2000.0] {
            let humid =
                GasState::try_new_moist_air_with_reference(temperature, 101_325.0, 0.0, reference)
                    .unwrap();
            let dry =
                GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap();
            assert_eq!(humid, dry, "zero vapor needs no condensed-phase branch");
        }
    }
}

#[test]
fn ice_relative_humidity_uses_iapws_and_the_original_mixture_owner() {
    let (temperature, pressure, humidity) = (230.0, 101_325.0, 0.8);
    let saturation =
        saturation_pressure_water_pa_with_reference(temperature, HumidityReference::Ice).unwrap();
    // IAPWS 2011 numerical check, rounded to the source's reported digits.
    assert!((saturation - 8.94735).abs() < 0.000005);
    let vapor_fraction = humidity * saturation / pressure;
    let composition = Composition::new(
        CompositionBasis::MoleFraction,
        vec![1.0 - vapor_fraction, vapor_fraction],
    )
    .unwrap();
    let general = GasState::try_new_mixture(
        &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()],
        &composition,
        temperature,
        pressure,
        Some(1),
    )
    .unwrap();
    let frost_relative = GasState::try_new_moist_air_with_reference(
        temperature,
        pressure,
        humidity,
        HumidityReference::Ice,
    )
    .unwrap();
    assert_eq!(general, frost_relative);
    assert!((frost_relative.water_mole_fraction * pressure - humidity * saturation).abs() < 1e-13);
    assert!(frost_relative.sound_speed.is_finite());
    assert!(frost_relative.dynamic_viscosity > 0.0);
    assert!(frost_relative.thermal_conductivity > 0.0);
    assert!(GasState::try_new_moist_air(temperature, pressure, humidity).is_err());
}

#[test]
fn humidity_basis_changes_reference_without_silently_switching_phase() {
    let temperature = 263.15;
    let pressure = 101_325.0;
    let liquid_pressure = saturation_pressure_water_pa(temperature).unwrap();
    let ice_pressure =
        saturation_pressure_water_pa_with_reference(temperature, HumidityReference::Ice).unwrap();
    assert!(ice_pressure < liquid_pressure);
    let liquid = GasState::try_new_moist_air(temperature, pressure, 0.7).unwrap();
    let ice = GasState::try_new_moist_air_with_reference(
        temperature,
        pressure,
        0.7,
        HumidityReference::Ice,
    )
    .unwrap();
    assert!(ice.water_mole_fraction < liquid.water_mole_fraction);
    // The SAME vapor pressure expressed on the other RH basis gives the same
    // physical mixture, up to the rounding of the RH conversion.
    let ice_humidity = 0.7 * liquid_pressure / ice_pressure;
    assert!(ice_humidity < 1.0);
    let same_vapor = GasState::try_new_moist_air_with_reference(
        temperature,
        pressure,
        ice_humidity,
        HumidityReference::Ice,
    )
    .unwrap();
    assert!((same_vapor.water_mole_fraction - liquid.water_mole_fraction).abs() < 1e-17);
    assert!((same_vapor.density - liquid.density).abs() < 1e-14);
    assert!((same_vapor.sound_speed - liquid.sound_speed).abs() < 1e-12);
}

#[test]
fn ice_reference_retains_domain_dilution_and_humidity_refusals() {
    for temperature in [49.99, 273.1601, 300.0, f64::NAN] {
        assert!(
            saturation_pressure_water_pa_with_reference(temperature, HumidityReference::Ice,)
                .is_err()
        );
        assert!(
            GasState::try_new_moist_air_with_reference(
                temperature,
                101_325.0,
                0.5,
                HumidityReference::Ice,
            )
            .is_err()
        );
    }
    for humidity in [-0.1, 1.01, f64::INFINITY, f64::NAN] {
        assert!(
            GasState::try_new_moist_air_with_reference(
                230.0,
                101_325.0,
                humidity,
                HumidityReference::Ice,
            )
            .is_err()
        );
    }
    assert!(
        GasState::try_new_moist_air_with_reference(230.0, 10.0, 1.0, HumidityReference::Ice,)
            .is_err(),
        "the existing 0.15 mole-fraction ceiling still applies"
    );
}
