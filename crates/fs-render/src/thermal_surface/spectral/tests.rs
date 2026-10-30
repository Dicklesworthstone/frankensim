//! Independent SI source, Kirchhoff equilibrium and declared-band admission.

use super::*;

fn source(temperature: f64, curve: SpectralEmissivity) -> SpectralDiffuseThermalSurface {
    SpectralDiffuseThermalSurface::try_new(temperature, curve, ContentHash([47; 32])).unwrap()
}

fn selective() -> SpectralEmissivity {
    SpectralEmissivity::table(&[(380.0, 0.1), (500.0, 0.7), (600.0, 0.2), (780.0, 1.0)]).unwrap()
}

// Independent dimensional expression, evaluated directly in metres. The
// production kernel uses a logarithmic wavelength-in-nm implementation.
fn si_blackbody(temperature: f64, wavelength_nm: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength_m = wavelength_nm * 1.0e-9;
    let exponent = h * c / (wavelength_m * k * temperature);
    2.0 * h * c * c / (wavelength_m.powi(5) * exponent.exp_m1()) * 1.0e-9
}

#[test]
fn g0_linear_knots_preserve_endpoints_support_and_exact_constant_laws() {
    let curve = selective();
    assert_eq!(
        curve.knots(),
        &[(380.0, 0.1), (500.0, 0.7), (600.0, 0.2), (780.0, 1.0)]
    );
    for &(wavelength, emissivity) in curve.knots() {
        assert_eq!(
            curve.at(wavelength).unwrap().to_bits(),
            emissivity.to_bits()
        );
    }
    for (wavelength, expected) in [(440.0, 0.4), (470.0, 0.55), (550.0, 0.45), (690.0, 0.6)] {
        assert!((curve.at(wavelength).unwrap() - expected).abs() < 2.0e-16);
    }
    for emissivity in [0.0, 0.17, 0.5, 1.0] {
        let curve = SpectralEmissivity::constant(emissivity).unwrap();
        for index in 0..=100 {
            assert_eq!(
                curve.at(380.0 + 4.0 * f64::from(index)).unwrap(),
                emissivity
            );
        }
    }
    let narrow = SpectralEmissivity::table(&[
        (380.0, 0.0),
        (553.0 - 1.0e-9, 0.0),
        (553.0, 1.0),
        (553.0 + 1.0e-9, 0.0),
        (780.0, 0.0),
    ])
    .unwrap();
    assert_eq!(narrow.at(552.0), Ok(0.0));
    assert_eq!(narrow.at(553.0), Ok(1.0));
    assert!(narrow.at(553.0 + 0.5e-9).unwrap() > 0.0);
}

#[test]
fn g1_spectral_emission_and_reciprocal_reflection_match_independent_si_equilibrium() {
    let curve = selective();
    for temperature in [1200.0, 1800.0, 6000.0] {
        let surface = source(temperature, curve);
        for (wavelength, epsilon) in [
            (380.0, 0.1),
            (440.0, 0.4),
            (500.0, 0.7),
            (550.0, 0.45),
            (600.0, 0.2),
            (690.0, 0.6),
            (780.0, 1.0),
        ] {
            let blackbody = si_blackbody(temperature, wavelength);
            let emitted = surface.radiance_w_m2_sr_nm(wavelength).unwrap();
            let reflected_fraction = surface.reflectance(wavelength).unwrap();
            assert!(
                (emitted / (epsilon * blackbody) - 1.0).abs() < 3.0e-13,
                "T={temperature}, lambda={wavelength}: {emitted} vs {}",
                epsilon * blackbody
            );
            assert!((reflected_fraction - (1.0 - epsilon)).abs() < 2.0e-16);
            // A surrounding blackbody at the same temperature remains a
            // blackbody: absorbed/reemitted and reflected fractions balance
            // separately at each wavelength, including the selective slopes.
            let outgoing = emitted + reflected_fraction * blackbody;
            assert!((outgoing / blackbody - 1.0).abs() < 3.0e-13);
        }
    }
}

#[test]
fn g3_constant_spectral_sources_and_gray_views_retain_exact_gray_radiance() {
    for temperature in [300.0, 1800.0, 6000.0] {
        for emissivity in [0.0, 0.17, 0.5, 1.0] {
            let gray =
                GrayDiffuseThermalSurface::try_new(temperature, emissivity, ContentHash([47; 32]))
                    .unwrap();
            let spectral = source(
                temperature,
                SpectralEmissivity::constant(emissivity).unwrap(),
            );
            let gray_view = DiffuseThermalSurface::from(gray);
            let spectral_view = DiffuseThermalSurface::from(&spectral);
            assert_eq!(gray_view.identity(), gray.identity());
            assert_eq!(spectral_view.identity(), spectral.identity());
            for wavelength in [380.0, 417.5, 500.0, 610.0, 779.5, 780.0] {
                let expected = gray.radiance_w_m2_sr_nm(wavelength).unwrap().to_bits();
                assert_eq!(
                    spectral.radiance_w_m2_sr_nm(wavelength).unwrap().to_bits(),
                    expected
                );
                assert_eq!(
                    gray_view.radiance_w_m2_sr_nm(wavelength).unwrap().to_bits(),
                    expected
                );
                assert_eq!(
                    spectral_view
                        .radiance_w_m2_sr_nm(wavelength)
                        .unwrap()
                        .to_bits(),
                    expected
                );
                assert_eq!(gray_view.reflectance(wavelength), Ok(gray.reflectance()));
                assert_eq!(
                    spectral_view.reflectance(wavelength),
                    Ok(gray.reflectance())
                );
            }
            // The compact gray view retains the pre-existing infrared domain.
            assert_eq!(
                gray_view.radiance_w_m2_sr_nm(10_000.0),
                gray.radiance_w_m2_sr_nm(10_000.0)
            );
        }
    }
}

#[test]
fn g0_tables_refuse_missing_duplicate_unsorted_excessive_or_nonphysical_knots() {
    for invalid in [
        vec![],
        vec![(380.0, 0.5)],
        vec![(381.0, 0.5), (780.0, 0.5)],
        vec![(380.0, 0.5), (779.0, 0.5)],
        vec![(380.0, 0.5), (500.0, 0.5), (500.0, 0.5), (780.0, 0.5)],
        vec![(380.0, 0.5), (600.0, 0.5), (500.0, 0.5), (780.0, 0.5)],
        vec![(380.0, 0.5), (f64::NAN, 0.5), (780.0, 0.5)],
        vec![(380.0, 0.5), (f64::INFINITY, 0.5), (780.0, 0.5)],
    ] {
        assert_eq!(
            SpectralEmissivity::table(&invalid),
            Err(ThermalSurfaceError::SpectralEmissivity)
        );
    }
    let maximum: Vec<_> = (0..16)
        .map(|index| (380.0 + f64::from(index) * (400.0 / 15.0), 0.5))
        .collect();
    // Construct the exact declared endpoint independently of rounding in the
    // evenly spaced interior, since silently clipping endpoints is forbidden.
    let mut maximum = maximum;
    maximum[15].0 = 780.0;
    assert_eq!(
        SpectralEmissivity::table(&maximum).unwrap().knots().len(),
        16
    );
    let excess: Vec<_> = (0..17)
        .map(|index| (380.0 + f64::from(index) * 25.0, 0.5))
        .collect();
    assert_eq!(
        SpectralEmissivity::table(&excess),
        Err(ThermalSurfaceError::SpectralEmissivity)
    );
    for invalid in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert_eq!(
            SpectralEmissivity::constant(invalid),
            Err(ThermalSurfaceError::Emissivity)
        );
    }
}

#[test]
fn g0_visible_law_never_extrapolates_even_when_it_is_constant_or_nonemitting() {
    for emissivity in [0.0, 0.5] {
        let curve = SpectralEmissivity::constant(emissivity).unwrap();
        let surface = source(1800.0, curve);
        let view = DiffuseThermalSurface::from(&surface);
        for wavelength in [
            0.0,
            -1.0,
            379.999,
            780.001,
            10_000.0,
            f64::NAN,
            f64::INFINITY,
        ] {
            assert_eq!(curve.at(wavelength), Err(ThermalSurfaceError::Wavelength));
            assert_eq!(
                surface.reflectance(wavelength),
                Err(ThermalSurfaceError::Wavelength)
            );
            assert_eq!(
                view.radiance_w_m2_sr_nm(wavelength),
                Err(ThermalSurfaceError::Wavelength)
            );
        }
    }
}

#[test]
fn g0_visible_source_representability_does_not_invent_broadband_cooling() {
    let curve = SpectralEmissivity::constant(0.5).unwrap();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            SpectralDiffuseThermalSurface::try_new(invalid, curve, ContentHash([47; 32])),
            Err(ThermalSurfaceError::Temperature)
        );
    }
    assert_eq!(
        SpectralDiffuseThermalSurface::try_new(1800.0, curve, ContentHash([0; 32])),
        Err(ThermalSurfaceError::SourceIdentity)
    );
    assert_eq!(
        SpectralDiffuseThermalSurface::try_new(f64::MAX, curve, ContentHash([47; 32])),
        Err(ThermalSurfaceError::Unrepresentable)
    );

    // The visible source is finite even where an all-wavelength T^4 cooling
    // power would overflow. No false bolometric limit enters this constructor.
    let high = source(1.0e100, curve);
    assert!(GrayDiffuseThermalSurface::try_new(1.0e100, 0.5, ContentHash([47; 32])).is_err());
    for wavelength in [380.0, 550.0, 780.0] {
        let got = high.radiance_w_m2_sr_nm(wavelength).unwrap();
        assert!(got > 0.0 && got.is_finite());
    }
    let dark = source(f64::MAX, SpectralEmissivity::constant(0.0).unwrap());
    assert!(!dark.is_emitting());
    assert_eq!(dark.radiance_w_m2_sr_nm(380.0), Ok(0.0));
    assert_eq!(dark.reflectance(780.0), Ok(1.0));
    let cold = source(1.0, curve);
    assert!(cold.is_emitting());
    assert_eq!(cold.radiance_w_m2_sr_nm(550.0), Ok(0.0));
}

#[test]
fn g3_identities_bind_source_temperature_and_exact_spectral_declaration() {
    let curve = selective();
    let original = source(1800.0, curve);
    assert_eq!(original.identity(), source(1800.0, curve).identity());
    assert_ne!(original.identity(), source(1900.0, curve).identity());
    assert_ne!(
        original.identity(),
        SpectralDiffuseThermalSurface::try_new(1800.0, curve, ContentHash([48; 32]))
            .unwrap()
            .identity()
    );
    let changed =
        SpectralEmissivity::table(&[(380.0, 0.1), (510.0, 0.7), (600.0, 0.2), (780.0, 1.0)])
            .unwrap();
    assert_ne!(curve.identity(), changed.identity());
    assert_ne!(original.identity(), source(1800.0, changed).identity());
    assert_eq!(
        SpectralEmissivity::constant(0.0).unwrap().identity(),
        SpectralEmissivity::constant(-0.0).unwrap().identity()
    );
    let constant = SpectralEmissivity::constant(0.5).unwrap();
    let explicitly_segmented =
        SpectralEmissivity::table(&[(380.0, 0.5), (500.0, 0.5), (780.0, 0.5)]).unwrap();
    assert_ne!(constant.identity(), explicitly_segmented.identity());
}
