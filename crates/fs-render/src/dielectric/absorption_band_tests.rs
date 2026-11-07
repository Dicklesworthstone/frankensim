//! G1/G3: original wavelength support, SI attenuation and shared bulk laws.

use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 16.0 * f64::EPSILON * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

fn infrared_ior() -> TabulatedIor {
    let samples = std::array::from_fn(|index| {
        DielectricIorSample::try_new_in_band(8000.0 + 750.0 * index as f64, 1.5).unwrap()
    });
    TabulatedIor::try_new_in_band(samples, 750.0).unwrap()
}

#[test]
fn infrared_bulk_retains_original_knots_and_physical_path_scaling() {
    let knots = [
        (8000.0, 0.0),
        (9990.0, 0.0),
        (10_000.0, 24.0),
        (14_000.0, 8.0),
    ];
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&knots).unwrap();
    let BeerLambertParameters::Spectral { spectrum } = absorption.parameters() else {
        panic!("the original absorption spectrum must be retained");
    };
    assert_eq!(spectrum.knots(), knots);
    assert_eq!(spectrum.wavelength_bounds_nm(), (8000.0, 14_000.0));
    assert_eq!(absorption.wavelength_bounds_nm(), Some((8000.0, 14_000.0)));
    for &(wavelength, value) in &knots {
        assert_eq!(absorption.extinction_per_m(wavelength), Ok(value));
    }
    assert_eq!(absorption.extinction_per_m(9995.0), Ok(12.0));
    assert_eq!(absorption.extinction_per_m(12_000.0), Ok(16.0));
    close(
        absorption.transmittance(12_000.0, 0.125).unwrap(),
        (-2.0_f64).exp(),
    );
    close(
        absorption.absorbed_fraction(12_000.0, 0.125).unwrap(),
        -(-2.0_f64).exp_m1(),
    );
    let first = absorption.transmittance(9995.0, 0.01).unwrap();
    let second = absorption.transmittance(9995.0, 0.025).unwrap();
    close(
        absorption.transmittance(9995.0, 0.035).unwrap(),
        first * second,
    );
    for wavelength in [0.0, -1.0, 780.0, 7999.0, 14_001.0, f64::NAN, f64::INFINITY] {
        // Even a zero-length segment must not hide unsupported material data.
        assert_eq!(
            absorption.transmittance(wavelength, 0.0),
            Err(DielectricError::InvalidWavelength)
        );
    }
}

#[test]
fn infrared_extinction_index_uses_vacuum_wavelength_and_intensity_decay() {
    let knots = [
        (8000.0, 1.0e-6),
        (10_000.0, 2.0e-6),
        (14_000.0, 1.0e-6),
    ];
    let absorption =
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&knots).unwrap();
    let BeerLambertParameters::ExtinctionIndex { spectrum } = absorption.parameters() else {
        panic!("the original imaginary-index source must be retained");
    };
    assert_eq!(spectrum.knots(), knots);
    close(spectrum.extinction_index(9000.0).unwrap(), 1.5e-6);
    // E(z) decays as exp(-2*pi*k*z/lambda_vac); intensity decays twice as fast.
    close(absorption.extinction_per_m(10_000.0).unwrap(), 0.8 * PI);
    close(
        absorption.transmittance(10_000.0, 0.25).unwrap(),
        (-0.2 * PI).exp(),
    );
    // The original k(9000 nm) is 1.5e-6. Retabulating alpha would incorrectly
    // return 0.65*pi instead of the source law's (2/3)*pi.
    close(absorption.extinction_per_m(9000.0).unwrap(), (2.0 / 3.0) * PI);
    assert_eq!(absorption.wavelength_bounds_nm(), Some((8000.0, 14_000.0)));
}

#[test]
fn constant_imaginary_index_retains_inverse_wavelength_opacity_and_path_composition() {
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        (400.0, 1.0e-6),
        (800.0, 1.0e-6),
    ])
    .unwrap();
    for wavelength in [400.0, 450.0, 600.0, 799.0, 800.0] {
        let alpha = (4.0 * PI * 1.0e3) / wavelength;
        close(absorption.extinction_per_m(wavelength).unwrap(), alpha);
        close(
            absorption.transmittance(wavelength, 0.025).unwrap(),
            (-alpha * 0.025).exp(),
        );
        close(
            absorption.absorbed_fraction(wavelength, 0.025).unwrap(),
            -(-alpha * 0.025).exp_m1(),
        );
        close(
            absorption.transmittance(wavelength, 0.025).unwrap(),
            absorption.transmittance(wavelength, 0.01).unwrap()
                * absorption.transmittance(wavelength, 0.015).unwrap(),
        );
    }
    // At 600 nm an interpolated-alpha substitute overstates opacity by 12.5%.
    let alpha_table = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (400.0, 10.0 * PI),
        (800.0, 5.0 * PI),
    ])
    .unwrap();
    close(
        alpha_table.extinction_per_m(600.0).unwrap()
            / absorption.extinction_per_m(600.0).unwrap(),
        1.125,
    );
    let material = |law| crate::tracer::Material::Dielectric {
        glass: DielectricGlass::new(
            CauchyIor::try_constant(1.5).unwrap(),
            law,
            GlassProvenance::Custom,
        ),
        surface: DielectricSurface::SMOOTH,
    };
    assert_ne!(
        material(absorption).content_identity(),
        material(alpha_table).content_identity(),
    );
}

#[test]
fn imaginary_index_source_refinement_keeps_the_same_constitutive_law() {
    // k = 1e-6 + lambda_nm * 1e-9, so alpha = 4*pi*(1000/lambda_nm + 1).
    let endpoints = [(400.0, 1.4e-6), (800.0, 1.8e-6)];
    let coarse =
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&endpoints).unwrap();
    let refined = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        endpoints[0],
        (550.0, 1.55e-6),
        (600.0, 1.6e-6),
        endpoints[1],
    ])
    .unwrap();
    for wavelength in [400.0, 479.0, 550.0, 600.0, 623.0, 800.0] {
        let expected = 4.0 * PI * (1000.0 / wavelength + 1.0);
        close(coarse.extinction_per_m(wavelength).unwrap(), expected);
        close(refined.extinction_per_m(wavelength).unwrap(), expected);
        close(
            coarse.transmittance(wavelength, 0.02).unwrap(),
            refined.transmittance(wavelength, 0.02).unwrap(),
        );
    }
}

#[test]
fn extinction_index_clear_and_thin_limits_retain_the_declared_source_band() {
    let clear = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        (400.0, -0.0),
        (800.0, 0.0),
    ])
    .unwrap();
    assert_ne!(clear.parameters(), BeerLambertParameters::Clear);
    let weak = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        (400.0, 1.0e-28),
        (800.0, 1.0e-28),
    ])
    .unwrap();
    for wavelength in [400.0, 600.0, 800.0] {
        assert_eq!(clear.extinction_per_m(wavelength).unwrap().to_bits(), 0);
        assert_eq!(clear.transmittance(wavelength, 1.0), Ok(1.0));
        assert_eq!(clear.absorbed_fraction(wavelength, 1.0), Ok(0.0));
        assert_eq!(weak.transmittance(wavelength, 1.0), Ok(1.0));
        close(
            weak.absorbed_fraction(wavelength, 1.0).unwrap(),
            (4.0 * PI * 1.0e-19) / wavelength,
        );
    }
    for absorption in [clear, weak] {
        for wavelength in [0.0, 399.0, 801.0, f64::INFINITY, f64::NAN] {
            assert_eq!(
                absorption.transmittance(wavelength, 0.0),
                Err(DielectricError::InvalidWavelength),
            );
        }
    }
}

#[test]
fn positive_extinction_index_cannot_underflow_into_a_clear_material() {
    let wavelengths = [1.0e100, 2.0e100];
    let positive = wavelengths.map(|wavelength| (wavelength, f64::from_bits(1)));
    assert_eq!(
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&positive),
        Err(DielectricError::InvalidAbsorption),
    );
    // An explicitly zero source remains exactly transparent on the same band.
    let clear = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(
        &wavelengths.map(|wavelength| (wavelength, 0.0)),
    )
    .unwrap();
    assert_eq!(clear.extinction_per_m(1.5e100), Ok(0.0));
    assert_eq!(clear.absorbed_fraction(1.5e100, 1.0), Ok(0.0));
}

#[test]
fn explicit_achromatic_laws_do_not_extrapolate_rgb_or_sampled_sources() {
    let constant = BeerLambertAbsorption::try_constant(2.0).unwrap();
    assert_eq!(constant.wavelength_bounds_nm(), None);
    assert_eq!(BeerLambertAbsorption::CLEAR.wavelength_bounds_nm(), None);
    for wavelength in [200.0, 550.0, 10_000.0, 1.0e6] {
        assert_eq!(
            BeerLambertAbsorption::CLEAR.extinction_per_m(wavelength),
            Ok(0.0)
        );
        assert_eq!(constant.extinction_per_m(wavelength), Ok(2.0));
        close(
            constant.transmittance(wavelength, 0.5).unwrap(),
            (-1.0_f64).exp(),
        );
    }
    let rgb = BeerLambertAbsorption::try_from_rgb_transmittance([0.8; 3], 1.0).unwrap();
    let table =
        BeerLambertAbsorption::try_from_spectral_extinction(&[(380.0, 2.0), (780.0, 2.0)]).unwrap();
    for absorption in [rgb, table] {
        assert_eq!(absorption.wavelength_bounds_nm(), Some((380.0, 780.0)));
        for wavelength in [200.0, 10_000.0] {
            assert_eq!(
                absorption.extinction_per_m(wavelength),
                Err(DielectricError::InvalidWavelength)
            );
        }
    }
    for absorption in [BeerLambertAbsorption::CLEAR, constant, rgb, table] {
        for wavelength in [-1.0, -0.0, 0.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                absorption.extinction_per_m(wavelength),
                Err(DielectricError::InvalidWavelength)
            );
        }
    }
}

#[test]
fn dielectric_common_support_intersects_both_original_laws() {
    let ior = infrared_ior();
    let make = |absorption| DielectricGlass::new(ior, absorption, GlassProvenance::Custom);
    for absorption in [
        BeerLambertAbsorption::CLEAR,
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
    ] {
        assert_eq!(
            make(absorption).wavelength_bounds_nm(),
            Ok((8000.0, 14_000.0))
        );
    }
    let narrowed = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (9000.0, 0.1),
        (12_000.0, 0.5),
    ])
    .unwrap();
    assert_eq!(
        make(narrowed).wavelength_bounds_nm(),
        Ok((9000.0, 12_000.0))
    );
    let partial = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (6000.0, 0.1),
        (10_000.0, 0.5),
    ])
    .unwrap();
    assert_eq!(make(partial).wavelength_bounds_nm(), Ok((8000.0, 10_000.0)));
    let endpoint = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (14_000.0, 0.1),
        (16_000.0, 0.5),
    ])
    .unwrap();
    assert_eq!(
        make(endpoint).wavelength_bounds_nm(),
        Ok((14_000.0, 14_000.0))
    );
    let disjoint =
        BeerLambertAbsorption::try_from_spectral_extinction(&[(380.0, 0.1), (780.0, 0.5)]).unwrap();
    assert_eq!(
        make(disjoint).wavelength_bounds_nm(),
        Err(DielectricError::InvalidWavelength)
    );
    let visible = DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        narrowed,
        GlassProvenance::Custom,
    );
    assert_eq!(
        visible.wavelength_bounds_nm(),
        Err(DielectricError::InvalidWavelength)
    );
}

#[test]
fn visible_constructors_preserve_canonical_spectral_parameters() {
    let knots = [(380.0, -0.0), (501.0, 2.0), (780.0, 4.0)];
    let legacy = BeerLambertAbsorption::try_from_spectral_extinction(&knots).unwrap();
    let explicit = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&knots).unwrap();
    assert_eq!(legacy, explicit);
    assert_eq!(legacy.parameters(), explicit.parameters());
    assert_eq!(
        SpectralAbsorption::try_new(&knots),
        SpectralAbsorption::try_new_in_band(&knots)
    );
    assert_eq!(legacy.extinction_per_m(380.0).unwrap().to_bits(), 0);
    let indices = [(380.0, 1.0e-6), (780.0, 2.0e-6)];
    assert_eq!(
        BeerLambertAbsorption::try_from_spectral_extinction_index(&indices),
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&indices)
    );
    assert_eq!(
        BeerLambertAbsorption::try_from_rgb_transmittance([1.0; 3], 1.0).unwrap(),
        BeerLambertAbsorption::CLEAR
    );
    assert_eq!(
        BeerLambertAbsorption::try_from_spectral_extinction(&[(8000.0, 0.1), (14_000.0, 0.2)]),
        Err(DielectricError::InvalidAbsorption)
    );
}

#[test]
fn declared_absorption_band_rejects_invalid_sources_and_preserves_clear_limits() {
    for knots in [
        vec![],
        vec![(8000.0, 0.0)],
        vec![(0.0, 0.0), (14_000.0, 1.0)],
        vec![(8000.0, 1.0), (8000.0, 0.0)],
        vec![(14_000.0, 0.0), (8000.0, 1.0)],
        vec![(8000.0, -1.0), (14_000.0, 1.0)],
        vec![(8000.0, f64::INFINITY), (14_000.0, 1.0)],
        vec![(8000.0, 0.0), (f64::NAN, 1.0)],
        (0..17)
            .map(|i| (8000.0 + 100.0 * f64::from(i), 1.0))
            .collect(),
    ] {
        assert_eq!(
            SpectralAbsorption::try_new_in_band(&knots),
            Err(DielectricError::InvalidAbsorption)
        );
    }
    let weak = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (8000.0, 1.0e-20),
        (14_000.0, 1.0e-20),
    ])
    .unwrap();
    assert_eq!(weak.transmittance(10_000.0, 1.0), Ok(1.0));
    close(weak.absorbed_fraction(10_000.0, 1.0).unwrap(), 1.0e-20);
    let tiny = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        (1.0e-310, 0.0),
        (2.0e-310, 1.0e-320),
    ])
    .unwrap();
    assert_eq!(tiny.extinction_per_m(1.0e-310), Ok(0.0));
    close(
        tiny.extinction_per_m(2.0e-310).unwrap(),
        (1.0e-320 / 2.0e-310) * (4.0 * PI * 1.0e9),
    );
    assert_eq!(
        BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
            (8000.0, f64::MAX),
            (14_000.0, f64::MAX)
        ]),
        Err(DielectricError::InvalidAbsorption)
    );
}
