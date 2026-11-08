use super::*;

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn complete_source_retains_exact_original_knots_and_condition() {
    assert_eq!(SOURCE_TEMPERATURE_K, 298.15);
    assert_eq!(
        HaleQuerryWaterAbsorption::wavelength_bounds_nm(),
        (200.0, 200_000.0)
    );
    let knots = HaleQuerryWaterAbsorption::knots();
    assert_eq!(knots.len(), 169);
    assert!(knots.windows(2).all(|w| w[0].0 < w[1].0));
    for &(wavelength, k) in knots {
        assert!(k.is_finite() && k > 0.0);
        assert_eq!(
            HaleQuerryWaterAbsorption::extinction_index(wavelength).unwrap(),
            k
        );
    }
    // Independent entries across the original Table I, including resonances.
    for (wavelength, k) in [
        (200.0, 1.10e-7),
        (475.0, 9.35e-10),
        (700.0, 3.35e-8),
        (975.0, 3.48e-6),
        (2800.0, 0.115),
        (3000.0, 0.272),
        (6000.0, 0.107),
        (17_000.0, 0.429),
        (50_000.0, 0.514),
        (130_000.0, 0.514),
        (200_000.0, 0.504),
    ] {
        assert_eq!(
            HaleQuerryWaterAbsorption::extinction_index(wavelength).unwrap(),
            k
        );
    }
    let identity = HaleQuerryWaterAbsorption::source_identity();
    assert_ne!(identity, ContentHash([0; 32]));
    assert_eq!(identity, HaleQuerryWaterAbsorption::source_identity());
}

#[test]
fn absorption_uses_intensity_units_and_interpolates_k_before_conversion() {
    // Independently evaluated 4*pi*k/lambda_m controls, not a table lookup.
    for (wavelength, expected) in [
        (400.0, 0.058_433_623_356_770_16),
        (500.0, 0.025_132_741_228_718_346),
        (975.0, 44.852_276_654_328_13),
        (3000.0, 1_139_350.935_701_898_4),
        (200_000.0, 31_667.253_948_185_115),
    ] {
        near(
            HaleQuerryWaterAbsorption::extinction_per_m(wavelength).unwrap(),
            expected,
            3.0e-15,
        );
    }
    let midpoint = 512.5;
    let k = HaleQuerryWaterAbsorption::extinction_index(midpoint).unwrap();
    near(k, 1.16e-9, 2.0e-16);
    let alpha = HaleQuerryWaterAbsorption::extinction_per_m(midpoint).unwrap();
    near(
        alpha,
        4.0 * core::f64::consts::PI * 1.16e-9 / 512.5e-9,
        3.0e-16,
    );
    let interpolated_alpha = 0.5
        * (HaleQuerryWaterAbsorption::extinction_per_m(500.0).unwrap()
            + HaleQuerryWaterAbsorption::extinction_per_m(525.0).unwrap());
    assert!((alpha - interpolated_alpha).abs() > 5.0e-5);
}

#[test]
fn source_refuses_invalid_queries_without_clamping_or_an_absorption_floor() {
    for wavelength in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        -1.0,
        f64::MIN_POSITIVE,
        199.999,
        200_000.001,
        f64::MAX,
    ] {
        assert_eq!(
            HaleQuerryWaterAbsorption::extinction_index(wavelength),
            Err(WaterAbsorptionError::WavelengthOutsideRange)
        );
        assert_eq!(
            HaleQuerryWaterAbsorption::extinction_per_m(wavelength),
            Err(WaterAbsorptionError::WavelengthOutsideRange)
        );
    }
    for wavelength in [200.0, 200_000.0] {
        assert!(
            HaleQuerryWaterAbsorption::extinction_per_m(wavelength)
                .unwrap()
                .is_finite()
        );
    }
    // The measured visible minimum remains small and positive.
    assert!(HaleQuerryWaterAbsorption::extinction_per_m(475.0).unwrap() < 0.025);
}
