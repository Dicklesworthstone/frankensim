use super::*;

fn relative_eq(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual} != {expected}",
    );
}

#[test]
fn g1_iapws_official_verification_point_and_triple_point() {
    let ice = WaterIceSublimation::try_new(917.0).unwrap();
    let point = ice.properties(230.0).unwrap();
    // Published verification value is rounded to six significant figures.
    relative_eq(point.saturation_pressure_pa, 8.947_35, 5.0e-7);
    relative_eq(
        ice.properties(273.16).unwrap().saturation_pressure_pa,
        611.657,
        1.0e-14,
    );
    relative_eq(
        ice.properties(50.0).unwrap().saturation_pressure_pa,
        1.934_958_486_808_867_2e-40,
        2.0e-13,
    );
    assert_eq!(point.solid_density_kg_m3, 917.0);
    assert_eq!(point.molar_mass_kg_mol, 0.018_015_28);
    assert_eq!(
        WaterIceSublimation::saturation_pressure_pa(230.0).unwrap(),
        point.saturation_pressure_pa,
    );
    assert!(WaterIceSublimation::saturation_pressure_pa(280.0).is_err());
}

fn check_derivatives(law: &impl SublimationLaw, temperature: f64) {
    let delta = 1.0e-4;
    let center = law.properties(temperature).unwrap();
    let lower = law.properties(temperature - delta).unwrap();
    let upper = law.properties(temperature + delta).unwrap();
    relative_eq(
        center.saturation_pressure_derivative_pa_k,
        (upper.saturation_pressure_pa - lower.saturation_pressure_pa) / (2.0 * delta),
        2.0e-8,
    );
    let heat_derivative = (upper.latent_heat_j_kg - lower.latent_heat_j_kg) / (2.0 * delta);
    assert!((center.latent_heat_derivative_j_kg_k - heat_derivative).abs() < 2.0e-5,);
    // Independent ideal-vapor Clausius-Clapeyron identity.
    relative_eq(
        center.latent_heat_j_kg,
        GAS_CONSTANT_J_MOL_K / center.molar_mass_kg_mol
            * temperature
            * temperature
            * center.saturation_pressure_derivative_pa_k
            / center.saturation_pressure_pa,
        1.0e-14,
    );
}

#[test]
fn g1_pressure_and_latent_heat_derivatives_use_the_same_curve() {
    let ice = WaterIceSublimation::try_new(917.0).unwrap();
    for temperature in [100.0, 200.0, 230.0, 260.0] {
        check_derivatives(&ice, temperature);
    }
    let co2 = CarbonDioxideSublimation::try_new(1560.0).unwrap();
    for temperature in [160.0, 180.0, 195.0] {
        check_derivatives(&co2, temperature);
    }
}

#[test]
fn g1_co2_antoine_units_and_latent_heat_trend() {
    let co2 = CarbonDioxideSublimation::try_new(1560.0).unwrap();
    let at_normal_sublimation = co2.properties(194.67).unwrap();
    // The NIST fit has pressure in bar. This approximate atmospheric anchor
    // detects an omitted bar-to-Pa conversion without reusing the formula.
    assert!((at_normal_sublimation.saturation_pressure_pa - 101_325.0).abs() < 1000.0,);
    let cold = co2.properties(154.26).unwrap();
    let warm = co2.properties(195.89).unwrap();
    assert!(cold.saturation_pressure_pa < warm.saturation_pressure_pa);
    assert!(cold.latent_heat_j_kg > warm.latent_heat_j_kg);
    assert!(warm.latent_heat_j_kg > 500_000.0);
    assert!(cold.latent_heat_j_kg < 650_000.0);
    assert!(warm.latent_heat_derivative_j_kg_k < 0.0);
}

#[test]
fn g1_ice_latent_heat_is_positive_and_not_artificially_constant() {
    let ice = WaterIceSublimation::try_new(917.0).unwrap();
    let cold = ice.properties(100.0).unwrap();
    let middle = ice.properties(230.0).unwrap();
    let warm = ice.properties(273.16).unwrap();
    assert!(cold.latent_heat_j_kg < middle.latent_heat_j_kg);
    assert!(warm.latent_heat_j_kg < middle.latent_heat_j_kg);
    assert!(cold.latent_heat_derivative_j_kg_k > 0.0);
    assert!(warm.latent_heat_derivative_j_kg_k < 0.0);
    for temperature in [50.0, 100.0, 150.0, 200.0, 250.0, 273.16] {
        let point = ice.properties(temperature).unwrap();
        assert!(point.saturation_pressure_derivative_pa_k > 0.0);
        assert!(point.latent_heat_j_kg > 2.6e6);
        assert!(point.latent_heat_j_kg < 2.9e6);
    }
}

#[test]
fn g1_generic_reference_and_constant_latent_heat() {
    let law = ClausiusClapeyronSublimation::try_new(
        200.0,
        100.0,
        500_000.0,
        0.04,
        1500.0,
        [150.0, 250.0],
    )
    .unwrap();
    assert_eq!(law.properties(200.0).unwrap().saturation_pressure_pa, 100.0,);
    let point = law.properties(225.0).unwrap();
    let expected =
        100.0 * ((500_000.0 * 0.04 / GAS_CONSTANT_J_MOL_K) * (1.0 / 200.0 - 1.0 / 225.0)).exp();
    relative_eq(point.saturation_pressure_pa, expected, 2.0e-14);
    assert_eq!(point.latent_heat_j_kg, 500_000.0);
    assert_eq!(point.latent_heat_derivative_j_kg_k, 0.0);
    check_derivatives(&law, 225.0);
}

#[test]
fn g0_domains_and_nonfinite_parameters_refuse() {
    let ice = WaterIceSublimation::try_new(917.0).unwrap();
    let co2 = CarbonDioxideSublimation::try_new(1560.0).unwrap();
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, 49.99, 273.160_001] {
        assert!(ice.properties(temperature).is_err());
    }
    for temperature in [154.25, 195.90, f64::NAN] {
        assert!(co2.properties(temperature).is_err());
    }
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(WaterIceSublimation::try_new(density).is_err());
        assert!(CarbonDioxideSublimation::try_new(density).is_err());
    }
    for (reference_t, reference_p, heat, mass, density, range) in [
        (200.0, 100.0, 500_000.0, 0.04, 1500.0, [250.0, 150.0]),
        (200.0, 100.0, 500_000.0, 0.04, 1500.0, [150.0, 150.0]),
        (200.0, 100.0, 500_000.0, 0.04, 1500.0, [0.0, 250.0]),
        (
            200.0,
            100.0,
            500_000.0,
            0.04,
            1500.0,
            [150.0, f64::INFINITY],
        ),
        (300.0, 100.0, 500_000.0, 0.04, 1500.0, [150.0, 250.0]),
        (f64::NAN, 100.0, 500_000.0, 0.04, 1500.0, [150.0, 250.0]),
        (200.0, 0.0, 500_000.0, 0.04, 1500.0, [150.0, 250.0]),
        (200.0, 100.0, -1.0, 0.04, 1500.0, [150.0, 250.0]),
        (200.0, 100.0, 500_000.0, 0.0, 1500.0, [150.0, 250.0]),
        (200.0, 100.0, 500_000.0, 0.04, f64::NAN, [150.0, 250.0]),
        // Underflow cannot become a calibrated zero-vapor-pressure material.
        (200.0, 100.0, 1.0e10, 0.04, 1500.0, [150.0, 250.0]),
    ] {
        assert!(
            ClausiusClapeyronSublimation::try_new(
                reference_t,
                reference_p,
                heat,
                mass,
                density,
                range,
            )
            .is_err(),
        );
    }
}
