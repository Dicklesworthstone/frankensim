use super::*;
use crate::conductor::{CoatingIorSample, ConductorIorSample};

fn source() -> ConductorSource {
    ConductorSource::try_new(
        fs_blake3::hash_domain(
            "org.frankensim.test.multilayer",
            b"declared-synthetic-optics",
        ),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

fn substrate(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        source(),
    )
    .unwrap()
}

fn layer(eta: f64, k: f64, thickness_nm: f64) -> CauchyCoatingLayer {
    CauchyCoatingLayer::try_new(
        CauchyIor::try_constant(eta).unwrap(),
        k,
        thickness_nm * 1.0e-9,
    )
    .unwrap()
}

fn absorbing(eta: f64, k: f64, thickness_nm: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        thickness_nm * 1.0e-9,
        source(),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn check_mueller(actual: Mueller, expected: [f64; 4], tolerance: f64) {
    for (value, reference) in [
        actual.0[0][0],
        actual.0[0][1],
        actual.0[2][2],
        actual.0[3][2],
    ]
    .into_iter()
    .zip(expected)
    {
        close(value, reference, tolerance);
    }
    assert_eq!(actual.0[0][1], actual.0[1][0]);
    assert_eq!(actual.0[2][3], -actual.0[3][2]);
}

#[test]
fn g1_arbitrary_order_matches_independent_complex_airy_and_signed_polarization() {
    // Independent nested-interface Airy oracle in Python cmath. This sums
    // r=(r_interface+r_below*exp(2i delta))/(1+r_interface*r_below*exp(2i delta)),
    // rather than propagating the production boundary-field vectors.
    let layers = [
        layer(2.1, 0.05, 45.0),
        layer(1.45, 0.0, 90.0),
        layer(2.5, 0.3, 20.0),
    ];
    let bare = substrate(0.17, 3.1);
    let coated =
        bare.with_periodic_coating(CoherentPeriodicCoating::try_new(&layers, 1, source()).unwrap());
    for (wavelength, cosine, expected) in [
        (
            380.0,
            1.0,
            [
                0.214_456_698_238_881_36,
                0.0,
                -0.214_456_698_238_881_36,
                0.0,
            ],
        ),
        (
            550.0,
            0.83,
            [
                0.723_204_065_298_868_3,
                0.050_531_840_327_248_945,
                -0.712_076_933_324_895_7,
                0.115_832_181_213_008_1,
            ],
        ),
        (
            780.0,
            0.17,
            [
                0.808_287_355_120_614_9,
                0.091_656_187_820_391_98,
                0.380_763_460_836_430_87,
                -0.707_069_146_953_801,
            ],
        ),
    ] {
        close(
            coated.fresnel(wavelength, 1.0, cosine).unwrap(),
            expected[0],
            2.0e-12,
        );
        check_mueller(
            coated.fresnel_mueller(wavelength, 1.0, cosine).unwrap(),
            expected,
            2.0e-12,
        );
    }
    let reversed = bare.with_periodic_coating(
        CoherentPeriodicCoating::try_new(&[layers[2], layers[1], layers[0]], 1, source()).unwrap(),
    );
    close(
        reversed.fresnel(550.0, 1.0, 0.83).unwrap(),
        0.746_782_604_347_273_6,
        2.0e-12,
    );
    assert_ne!(coated.content_identity(), reversed.content_identity());
}

#[test]
fn g1_bragg_quarter_wave_stack_has_closed_form_admittance() {
    let bare = substrate(0.25, 2.7);
    let (high, low, wavelength) = (2.4_f64, 1.4_f64, 560.0);
    for repeats in [1, 3, 8, 32] {
        for reverse in [false, true] {
            let (first, second) = if reverse { (low, high) } else { (high, low) };
            let cell = [
                layer(first, 0.0, wavelength / (4.0 * first)),
                layer(second, 0.0, wavelength / (4.0 * second)),
            ];
            let coated = bare.with_periodic_coating(
                CoherentPeriodicCoating::try_new(&cell, repeats, source()).unwrap(),
            );
            // Two quarter-wave layers transform normal-incidence admittance
            // by (n_first/n_second)^2, independently of their common phase.
            let mut gain = 1.0;
            for _ in 0..repeats {
                gain *= (first / second).powi(2);
            }
            let admittance = Complex::new(0.25, 2.7) * gain;
            let amplitude =
                (Complex::new(1.0, 0.0) - admittance) / (Complex::new(1.0, 0.0) + admittance);
            close(
                coated.fresnel(wavelength, 1.0, 1.0).unwrap(),
                amplitude.norm_squared(),
                3.0e-13,
            );
            let matrix = coated.fresnel_mueller(wavelength, 1.0, 1.0).unwrap();
            close(matrix.0[0][0], amplitude.norm_squared(), 3.0e-13);
            assert_eq!(matrix.0[0][1], 0.0);
            assert_eq!(matrix.0[2][2], -matrix.0[0][0]);
        }
    }
}

#[test]
fn g1_borrowed_stack_keeps_full_eta_k_spectra_and_layer_sources() {
    let films = [
        (1.65, 0.001, 0.04, 0.0008, 62.0),
        (2.05, -0.0003, 0.18, -0.0002, 28.0),
    ]
    .map(|(eta, slope, k, k_slope, thickness)| {
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                let offset = 50.0 * i as f64;
                CoatingIorSample::try_new(
                    380.0 + offset,
                    eta + slope * offset,
                    k + k_slope * offset,
                )
                .unwrap()
            }),
            thickness * 1.0e-9,
            source(),
        )
        .unwrap()
    });
    let layers = films.each_ref().map(CoatingLayer::Absorbing);
    let stack = CoherentCoatingStack::try_new(&layers).unwrap();
    let bare = substrate(0.25, 3.7);
    let expected = [
        0.265_739_642_769_449,
        -0.064_202_235_846_048_05,
        -0.254_295_690_241_335_7,
        -0.042_770_697_636_179_2,
    ];
    close(
        stack.fresnel(&bare, 505.0, 1.33, 0.63).unwrap(),
        expected[0],
        2.0e-12,
    );
    check_mueller(
        stack.fresnel_mueller(&bare, 505.0, 1.33, 0.63).unwrap(),
        expected,
        2.0e-12,
    );
    let reordered = [layers[1], layers[0]];
    let reversed = CoherentCoatingStack::try_new(&reordered).unwrap();
    assert_ne!(stack.content_identity(), reversed.content_identity());
    assert!(
        (stack.fresnel(&bare, 505.0, 1.33, 0.63).unwrap()
            - reversed.fresnel(&bare, 505.0, 1.33, 0.63).unwrap())
        .abs()
            > 0.01
    );
}

#[test]
fn g3_repeated_dispersive_cell_equals_explicit_arbitrary_stack() {
    let laws = [
        CauchyIor::try_new(1.4, 0.012, 0.0002).unwrap(),
        CauchyIor::try_new(2.0, 0.021, 0.0004).unwrap(),
    ];
    let cell = [
        CauchyCoatingLayer::try_new(laws[0], 0.0, 90.0e-9).unwrap(),
        CauchyCoatingLayer::try_new(laws[1], 0.0, 50.0e-9).unwrap(),
    ];
    let films = [
        CoherentDielectricCoating::try_new(laws[0], 90.0e-9, source()).unwrap(),
        CoherentDielectricCoating::try_new(laws[1], 50.0e-9, source()).unwrap(),
    ];
    let bare = substrate(0.3, 3.0);
    let compact =
        bare.with_periodic_coating(CoherentPeriodicCoating::try_new(&cell, 32, source()).unwrap());
    let expanded: [CoatingLayer<'_>; 64] =
        core::array::from_fn(|i| CoatingLayer::Dielectric(&films[i % 2]));
    let general = CoherentCoatingStack::try_new(&expanded).unwrap();
    for wavelength in [380.0, 425.0, 550.0, 780.0] {
        for cosine in [0.0, 0.2, 0.8, 1.0] {
            assert_eq!(
                compact.fresnel(wavelength, 1.5, cosine).unwrap().to_bits(),
                general
                    .fresnel(&bare, wavelength, 1.5, cosine)
                    .unwrap()
                    .to_bits()
            );
            assert_eq!(
                compact.fresnel_mueller(wavelength, 1.5, cosine).unwrap(),
                general
                    .fresnel_mueller(&bare, wavelength, 1.5, cosine)
                    .unwrap()
            );
        }
    }
}

#[test]
fn g1_critical_evanescent_and_weakly_absorbing_layers_are_continuous() {
    let bare = substrate(0.25, 3.7);
    let cosine = (1.0_f64 - (1.0_f64 / 1.5).powi(2)).sqrt();
    // Independent critical limit: replace the first layer by the analytic
    // admittance transform Y/(1-i*k0*d*Y); obtain Y below it from a single
    // Airy interface pair. No 0/0 interface ratios are evaluated at q=0.
    let expected = [
        0.732_359_807_014_451_3,
        0.128_992_986_469_564_53,
        -0.382_675_552_272_025,
        0.610_959_178_722_447_4,
    ];
    for loss in [0.0, 1.0e-14] {
        let coated = bare.with_periodic_coating(
            CoherentPeriodicCoating::try_new(
                &[layer(1.0, loss, 85.0), layer(1.9, 0.13, 37.0)],
                1,
                source(),
            )
            .unwrap(),
        );
        close(
            coated.fresnel(550.0, 1.5, cosine).unwrap(),
            expected[0],
            2.0e-12,
        );
        check_mueller(
            coated.fresnel_mueller(550.0, 1.5, cosine).unwrap(),
            expected,
            2.0e-12,
        );
        for adjacent in [cosine - 1.0e-9, cosine + 1.0e-9] {
            close(
                coated.fresnel(550.0, 1.5, adjacent).unwrap(),
                expected[0],
                2.0e-8,
            );
        }
    }
}

#[test]
fn g3_zero_layers_and_matched_boundaries_keep_bare_bits() {
    let bare = substrate(1.8, 0.6);
    let zero = absorbing(2.2, 0.3, 0.0);
    let zero_layers = [CoatingLayer::Absorbing(&zero); 64];
    for layers in [&[][..], &zero_layers[..]] {
        let stack = CoherentCoatingStack::try_new(layers).unwrap();
        for cosine in [0.0, 0.2, 1.0] {
            assert_eq!(
                stack.fresnel(&bare, 550.0, 1.33, cosine).unwrap().to_bits(),
                bare.fresnel(550.0, 1.33, cosine).unwrap().to_bits()
            );
            assert_eq!(
                stack.fresnel_mueller(&bare, 550.0, 1.33, cosine).unwrap(),
                bare.fresnel_mueller(550.0, 1.33, cosine).unwrap()
            );
        }
    }
    let matched = bare.with_periodic_coating(
        CoherentPeriodicCoating::try_new(
            &[
                layer(1.33, 0.0, 170.0),
                layer(2.2, 0.4, 0.0),
                layer(1.8, 0.6, 80.0),
            ],
            1,
            source(),
        )
        .unwrap(),
    );
    for cosine in [0.0, 0.35, 1.0] {
        assert_eq!(
            matched.fresnel(550.0, 1.33, cosine),
            bare.fresnel(550.0, 1.33, cosine)
        );
        assert_eq!(
            matched.fresnel_mueller(550.0, 1.33, cosine),
            bare.fresnel_mueller(550.0, 1.33, cosine)
        );
    }
}

#[test]
fn g3_single_film_and_split_layer_limits_preserve_existing_response() {
    let bare = substrate(0.2, 3.0);
    for (eta, k) in [(1.7, 0.0), (1.8, 0.4), (1.0, 0.0)] {
        let original = bare.with_absorbing_coating(absorbing(eta, k, 130.0));
        let one = bare.with_periodic_coating(
            CoherentPeriodicCoating::try_new(&[layer(eta, k, 130.0)], 1, source()).unwrap(),
        );
        let split = bare.with_periodic_coating(
            CoherentPeriodicCoating::try_new(
                &[
                    layer(eta, k, 35.0),
                    layer(2.7, 0.5, 0.0),
                    layer(eta, k, 95.0),
                ],
                1,
                source(),
            )
            .unwrap(),
        );
        for wavelength in [380.0, 550.0, 780.0] {
            for cosine in [0.0, 0.2, 0.75, 1.0] {
                let expected = original.fresnel(wavelength, 1.5, cosine).unwrap();
                close(
                    one.fresnel(wavelength, 1.5, cosine).unwrap(),
                    expected,
                    2.0e-12,
                );
                close(
                    split.fresnel(wavelength, 1.5, cosine).unwrap(),
                    expected,
                    2.0e-12,
                );
                let reference = original.fresnel_mueller(wavelength, 1.5, cosine).unwrap();
                for row in 0..4 {
                    for col in 0..4 {
                        close(
                            split.fresnel_mueller(wavelength, 1.5, cosine).unwrap().0[row][col],
                            reference.0[row][col],
                            2.0e-12,
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g1_opaque_front_layer_screens_arbitrary_deep_layers_without_overflow() {
    let bare = substrate(0.25, 3.7);
    let front_halfspace = substrate(1.8, 0.6);
    let coating = CoherentPeriodicCoating::try_new(
        &[
            CauchyCoatingLayer::try_new(
                CauchyIor::try_constant(1.8).unwrap(),
                0.6,
                MAX_COHERENT_COATING_THICKNESS_M,
            )
            .unwrap(),
            layer(1.0, 0.0, 80.0),
            layer(4.0, 0.0, 31.0),
            layer(1.5, 100.0, 99_999.0),
        ],
        16,
        source(),
    )
    .unwrap();
    let coated = bare.with_periodic_coating(coating);
    for wavelength in [380.0, 550.0, 780.0] {
        for cosine in [1.0e-5, 0.2, 0.8, 1.0] {
            close(
                coated.fresnel(wavelength, 1.33, cosine).unwrap(),
                front_halfspace.fresnel(wavelength, 1.33, cosine).unwrap(),
                2.0e-12,
            );
            let actual = coated.fresnel_mueller(wavelength, 1.33, cosine).unwrap();
            let expected = front_halfspace
                .fresnel_mueller(wavelength, 1.33, cosine)
                .unwrap();
            for row in 0..4 {
                for col in 0..4 {
                    close(actual.0[row][col], expected.0[row][col], 2.0e-12);
                }
            }
        }
    }
}

#[test]
fn g0_bounded_compact_storage_and_identity_bind_authored_optics() {
    assert!(
        core::mem::size_of::<CoherentPeriodicCoating>()
            <= core::mem::size_of::<CoherentAbsorbingCoating>()
    );
    let first = layer(2.1, 0.0, 60.0);
    let second = layer(1.4, 0.1, 90.0);
    assert_eq!(
        CoherentPeriodicCoating::try_new(&[], 1, source()),
        Err(ConductorError::InvalidCoatingStack)
    );
    assert_eq!(
        CoherentPeriodicCoating::try_new(&[first; 5], 1, source()),
        Err(ConductorError::InvalidCoatingStack)
    );
    for repeats in [0, 33, usize::MAX] {
        assert_eq!(
            CoherentPeriodicCoating::try_new(&[first, second], repeats, source()),
            Err(ConductorError::InvalidCoatingStack)
        );
    }
    let film = absorbing(1.4, 0.1, 90.0);
    let over_budget = [CoatingLayer::Absorbing(&film); 65];
    assert_eq!(
        CoherentCoatingStack::try_new(&over_budget),
        Err(ConductorError::InvalidCoatingStack)
    );
    for invalid in [f64::NAN, -0.01, 101.0] {
        assert!(CauchyCoatingLayer::try_new(first.ior(), invalid, 1.0e-7).is_err());
    }
    let bare = substrate(0.2, 3.0);
    let identity = |cell: &[CauchyCoatingLayer], repeats, source| {
        bare.with_periodic_coating(CoherentPeriodicCoating::try_new(cell, repeats, source).unwrap())
            .content_identity()
    };
    let original = identity(&[first, second], 3, source());
    for changed in [
        identity(&[second, first], 3, source()),
        identity(&[first, second], 4, source()),
        identity(&[layer(2.1, 0.01, 60.0), second], 3, source()),
        identity(&[layer(2.1, 0.0, 61.0), second], 3, source()),
        identity(
            &[
                CauchyCoatingLayer::try_new(
                    CauchyIor::try_new(2.1, 0.001, 0.0).unwrap(),
                    0.0,
                    60.0e-9,
                )
                .unwrap(),
                second,
            ],
            3,
            source(),
        ),
        identity(
            &[first, second],
            3,
            ConductorSource::try_new(
                source().identity(),
                ConductorDataStatus::CallerAssertedMeasured,
            )
            .unwrap(),
        ),
    ] {
        assert_ne!(original, changed);
    }
}
