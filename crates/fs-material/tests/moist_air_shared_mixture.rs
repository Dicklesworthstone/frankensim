//! Routing regressions, not new experimental-accuracy claims. Both humidity
//! entry points ultimately use the shared gas-mixture constitutive kernel.

use fs_material::gas::{GasSpec, GasState, saturation_pressure_water_pa};
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
