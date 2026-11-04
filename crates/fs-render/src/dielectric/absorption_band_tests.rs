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
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_index_in_band(&[
        (8000.0, 1.0e-6),
        (10_000.0, 2.0e-6),
        (14_000.0, 1.0e-6),
    ])
    .unwrap();
    // E(z) decays as exp(-2*pi*k*z/lambda_vac); intensity decays twice as fast.
    close(absorption.extinction_per_m(10_000.0).unwrap(), 0.8 * PI);
    close(
        absorption.transmittance(10_000.0, 0.25).unwrap(),
        (-0.2 * PI).exp(),
    );
    // The retained law interpolates alpha, not k. This distinction matters
    // between knots; no continuous complex-index fit has been introduced.
    close(absorption.extinction_per_m(9000.0).unwrap(), 0.65 * PI);
    assert_eq!(absorption.wavelength_bounds_nm(), Some((8000.0, 14_000.0)));
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
