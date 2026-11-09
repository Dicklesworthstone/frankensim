use super::*;

const B: [f64; 3] = [1.03961212, 0.231792344, 1.01046945];
const C: [f64; 3] = [0.00600069867, 0.0200179144, 103.560653];

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance {tolerance:e}"
    );
}
fn air(t: f64, p: f64, water: f64, co2: f64) -> CiddorAirState {
    CiddorAirState::try_new(CiddorAir::try_new(water, co2).unwrap(), t, p).unwrap()
}
fn law(reference: CiddorAirState, coordinate: CatalogWavelength) -> AirRelativeSellmeier {
    AirRelativeSellmeier::try_new(B, C, (320.0, 1800.0), reference, coordinate).unwrap()
}

// Independent oscillator equations, in the source's declared coordinate.
fn relative_reference(lambda_nm: f64) -> (f64, f64) {
    let lambda = lambda_nm * 0.001;
    let s = lambda * lambda;
    let mut n2 = 1.0;
    let mut derivative = 0.0;
    for (b, c) in B.into_iter().zip(C) {
        n2 += b * s / (s - c);
        derivative -= 2.0 * b * lambda * c / (s - c).powi(2);
    }
    let n = n2.sqrt();
    (n, n - lambda * derivative / (2.0 * n))
}

#[test]
fn g1_absolute_catalog_jet_keeps_independent_atmospheres_and_ordinary_curvature() {
    let reference = air(293.15, 101_330.0, 0.012, 650.0);
    let coordinate_air = air(288.15, 101_325.0, 0.0, 450.0);
    for coordinate in [
        CatalogWavelength::Vacuum,
        CatalogWavelength::Air(reference),
        CatalogWavelength::Air(coordinate_air),
    ] {
        let law = law(reference, coordinate);
        let phase = |lambda| {
            let a = reference.sample(lambda).unwrap().phase_index;
            let b = match coordinate {
                CatalogWavelength::Vacuum => 1.0,
                CatalogWavelength::Air(air) => air.sample(lambda).unwrap().phase_index,
            };
            a * relative_reference(lambda / b).0
        };
        for lambda in [350.0, 550.0, 1000.0, 1600.0] {
            let sample = law.refractive_index_jet(lambda).unwrap();
            let a = reference.sample(lambda).unwrap();
            let (b, gb) = match coordinate {
                CatalogWavelength::Vacuum => (1.0, 1.0),
                CatalogWavelength::Air(air) => {
                    let b = air.sample(lambda).unwrap();
                    (b.phase_index, b.group_index)
                }
            };
            let (r, gr) = relative_reference(lambda / b);
            near(sample.phase_index, a.phase_index * r, 1e-14);
            near(
                sample.group_index,
                a.group_index * r + (a.phase_index * gb / b) * (gr - r),
                1e-14,
            );
            assert_eq!(sample.phase_index, law.eval(lambda).unwrap());
            assert_eq!(sample.group_index, law.group_index(lambda).unwrap());
            let h = 0.2;
            let first = (phase(lambda - 2.0 * h) - 8.0 * phase(lambda - h)
                + 8.0 * phase(lambda + h)
                - phase(lambda + 2.0 * h))
                / (12.0 * h);
            let second = (-phase(lambda - 2.0 * h) + 16.0 * phase(lambda - h)
                - 30.0 * phase(lambda)
                + 16.0 * phase(lambda + h)
                - phase(lambda + 2.0 * h))
                / (12.0 * h * h);
            near(
                sample.scaled_wavelength_derivatives[0],
                lambda * first,
                5e-11,
            );
            near(
                sample.scaled_wavelength_derivatives[1],
                lambda * lambda * second,
                6e-8,
            );
            near(
                sample.group_index,
                sample.phase_index - sample.scaled_wavelength_derivatives[0],
                8e-16,
            );
        }
    }
}

#[test]
fn g0_catalog_admits_the_full_source_band_and_maps_coordinate_edges_inward() {
    let reference = air(293.15, 101_330.0, 0.012, 650.0);
    let coordinate_air = air(288.15, 101_325.0, 0.0, 450.0);
    let coordinate = CatalogWavelength::Air(coordinate_air);
    let original =
        AirRelativeSellmeier::try_new(B, C, (380.0, 780.0), reference, coordinate).unwrap();
    let (lo, hi) = original.wavelength_bounds_nm();
    let mapped = |lambda| lambda / coordinate_air.sample(lambda).unwrap().phase_index;
    assert!(mapped(lo) >= 380.0 && mapped(hi) <= 780.0);
    assert!(mapped(lo.next_down()) < 380.0);
    assert!(mapped(hi.next_up()) > 780.0);
    for edge in [lo, hi] {
        assert!(original.refractive_index_jet(edge).is_ok());
    }
    for outside in [lo.next_down(), hi.next_up(), f64::NAN] {
        assert_eq!(
            original.refractive_index_jet(outside),
            Err(CatalogRefractionError::InvalidWavelength)
        );
    }
    // An active 2000 nm pole is outside Ciddor's intersection, but inside the
    // ORIGINAL catalog fit, so intersecting first must never hide that pole.
    assert_eq!(
        AirRelativeSellmeier::try_new(
            [1.25, 0.1, 0.0],
            [0.0, 4.0, 0.0],
            (320.0, 2200.0),
            reference,
            coordinate
        ),
        Err(CatalogRefractionError::InvalidIndex)
    );
    assert_eq!(
        AirRelativeSellmeier::try_new(B, C, (1900.0, 2100.0), reference, coordinate),
        Err(CatalogRefractionError::InvalidWavelength)
    );
    assert_eq!(
        law(reference, CatalogWavelength::Vacuum).wavelength_bounds_nm(),
        (320.0, 1700.0)
    );
}

#[test]
fn g3_constant_relative_glass_keeps_air_curvature_and_complete_source_identity() {
    let reference = air(302.15, 95_000.0, 0.012, 650.0);
    let other = air(288.15, 101_325.0, 0.0, 450.0);
    for coordinate in [
        CatalogWavelength::Vacuum,
        CatalogWavelength::Air(reference),
        CatalogWavelength::Air(other),
    ] {
        let law = AirRelativeSellmeier::try_new(
            [1.25, 0.0, -0.0],
            [0.0, 1.0, 2.0],
            (200.0, 2000.0),
            reference,
            coordinate,
        )
        .unwrap();
        assert_eq!(law.coefficients(), ([1.25, 0.0, 0.0], [0.0; 3]));
        for lambda in [400.0, 550.0, 1550.0] {
            let sample = law.refractive_index_jet(lambda).unwrap();
            let a = reference.sample(lambda).unwrap();
            near(sample.group_index, 1.5 * a.group_index, 5e-16);
            near(
                sample.scaled_wavelength_derivatives[0],
                1.5 * lambda * a.wavelength_derivative_per_nm,
                2e-19,
            );
            near(
                sample.scaled_wavelength_derivatives[1],
                1.5 * lambda * lambda * a.wavelength_second_derivative_per_nm2,
                2e-19,
            );
            assert!(sample.scaled_wavelength_derivatives[1] > 0.0);
        }
    }
    let base = law(reference, CatalogWavelength::Air(other));
    assert_eq!(base, law(reference, CatalogWavelength::Air(other)));
    assert_ne!(
        base.content_identity(),
        law(other, CatalogWavelength::Air(other)).content_identity()
    );
    assert_ne!(
        base.content_identity(),
        law(reference, CatalogWavelength::Vacuum).content_identity()
    );
    assert_ne!(
        base.content_identity(),
        law(reference, CatalogWavelength::Air(reference)).content_identity()
    );
}
