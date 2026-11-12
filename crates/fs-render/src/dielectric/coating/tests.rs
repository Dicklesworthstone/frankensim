use super::*;
use crate::conductor::{CauchyCoatingLayer, ConductorDataStatus, ConductorSource};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricIorSample,
    DielectricSurface, GlassProvenance, SellmeierIor, TabulatedIor,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use fs_geom::Point3;

fn source() -> ConductorSource {
    ConductorSource::try_new(
        fs_blake3::hash_domain(
            "org.frankensim.test.coated-dielectric",
            b"authored-transparent-films",
        ),
        ConductorDataStatus::Representative,
    )
    .unwrap()
}

fn layer(index: f64, thickness_nm: f64) -> CauchyCoatingLayer {
    CauchyCoatingLayer::try_new(
        CauchyIor::try_constant(index).unwrap(),
        0.0,
        thickness_nm * 1.0e-9,
    )
    .unwrap()
}

fn surface(layers: &[CauchyCoatingLayer], repeats: usize) -> DielectricSurface {
    DielectricSurface::SMOOTH.with_coating(
        DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(layers, repeats, source()).unwrap(),
        )
        .unwrap(),
    )
}

fn glass(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn entries(matrix: Mueller) -> [f64; 4] {
    [
        matrix.0[0][0],
        matrix.0[0][1],
        matrix.0[2][2],
        matrix.0[3][2],
    ]
}

#[test]
fn g3_bare_interfaces_keep_native_smooth_and_rough_transport_outside_visible_light() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.3, 0.0, (1.0_f64 - 0.09).sqrt());
    let bare = DielectricSurface::try_rough(0.2).unwrap();
    for reverse in [false, true] {
        let (eta_i, eta_t) = if reverse { (1.5, 1.0) } else { (1.0, 1.5) };
        for wavelength in [250.0, 550.0, 1550.0, 10_000.0] {
            let interface = bare.interface(wavelength, reverse).unwrap();
            for branch in [0.0, 0.9] {
                assert_eq!(
                    interface
                        .sample_smooth(normal, wo, eta_i, eta_t, branch)
                        .unwrap(),
                    crate::dielectric::sample_smooth_dielectric(normal, wo, eta_i, eta_t, branch,)
                        .unwrap(),
                );
                let native = crate::dielectric::sample_rough_dielectric(
                    normal, wo, eta_i, eta_t, 0.2, 0.17, 0.63, branch,
                )
                .unwrap();
                let sample = interface
                    .sample_rough(normal, wo, eta_i, eta_t, 0.2, 0.17, 0.63, branch)
                    .unwrap();
                assert_eq!(sample, native);
                if let Some(sample) = sample {
                    assert_eq!(
                        interface
                            .evaluate_rough(normal, wo, sample.direction, eta_i, eta_t, 0.2)
                            .unwrap(),
                        crate::dielectric::evaluate_rough_dielectric(
                            normal,
                            wo,
                            sample.direction,
                            eta_i,
                            eta_t,
                            0.2,
                        )
                        .unwrap(),
                    );
                }
            }
        }
    }
    for invalid in [0.0, -0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            bare.interface(invalid, false),
            Err(DielectricError::InvalidWavelength)
        ));
    }
}

#[test]
fn g0_nonvisible_bare_admission_does_not_extend_original_film_support() {
    // A zero-thickness compact film still owns its original Cauchy source.
    for thickness_nm in [0.0, 120.0] {
        let coated = surface(&[layer(1.4, thickness_nm)], 1);
        for reverse in [false, true] {
            assert!(coated.interface(550.0, reverse).is_ok());
            assert!(matches!(
                coated.interface(1550.0, reverse),
                Err(DielectricError::InvalidWavelength)
            ));
        }
    }
    let empty = DielectricCoatingStack::try_new_in_band(&[], (1200.0, 1800.0)).unwrap();
    let indices = core::array::from_fn(|i| {
        DielectricIorSample::try_new_in_band(1000.0 + 125.0 * i as f64, 1.4).unwrap()
    });
    let layer = CoherentDielectricCoating::try_new(
        TabulatedIor::try_new_in_band(indices, 125.0).unwrap(),
        120e-9,
        source(),
    )
    .unwrap();
    let layers = [layer];
    let infrared = DielectricCoatingStack::try_new_in_band(&layers, (1200.0, 1800.0)).unwrap();
    for stack in [empty, infrared] {
        for reverse in [false, true] {
            assert!(
                stack
                    .interface(1550.0, reverse)
                    .unwrap()
                    .fresnel(0.8, 1.0, 1.5)
                    .is_ok()
            );
            for outside in [550.0, 1199.0, 1801.0] {
                assert!(matches!(
                    stack.interface(outside, reverse),
                    Err(DielectricError::InvalidWavelength)
                ));
            }
        }
    }
}

#[test]
fn g1_transmitted_and_reflected_signed_fields_match_independent_airy() {
    // Independent Python cmath nested-interface scattering recursion, retaining
    // t01*t_below*exp(i delta)/(1+r01*r_below*exp(2i delta)). Transmission is
    // normalized by Re(Y_exit)/Y_incident separately for s and p.
    let coated = surface(&[layer(1.3, 72.0), layer(2.1, 43.0), layer(1.45, 110.0)], 1);
    for (wavelength, cosine, eta_i, eta_t, reflection, transmission) in [
        (
            550.0,
            0.83,
            1.0,
            1.7,
            [
                0.200_151_050_961_076_55,
                0.053_077_043_112_156_04,
                -0.191_495_167_719_472_8,
                0.023_934_732_824_595_795,
            ],
            [
                0.799_848_949_038_923_8,
                -0.053_077_043_112_156_22,
                0.797_898_309_929_365_7,
                -0.017_304_790_810_729_61,
            ],
        ),
        (
            415.0,
            0.31,
            1.0,
            1.7,
            [
                0.327_484_691_202_295_3,
                0.207_056_733_824_791_28,
                0.057_502_950_401_120_606,
                0.247_117_669_632_987_76,
            ],
            [
                0.672_515_308_797_704_9,
                -0.207_056_733_824_791_48,
                0.607_690_561_203_008_1,
                -0.200_291_116_552_519_1,
            ],
        ),
        (
            780.0,
            0.87,
            1.5,
            1.0,
            [
                0.131_691_385_773_215_35,
                0.091_292_489_068_166_34,
                -0.087_543_753_823_200_01,
                -0.036_665_974_597_049_91,
            ],
            [
                0.868_308_614_226_784_6,
                -0.091_292_489_068_165_95,
                0.863_496_085_833_234_3,
                0.000_201_818_883_152_349_3,
            ],
        ),
    ] {
        let response = coated
            .interface(wavelength, false)
            .unwrap()
            .fresnel_polarized(cosine, eta_i, eta_t)
            .unwrap();
        for (actual, expected) in entries(response.reflection)
            .into_iter()
            .zip(reflection)
            .chain(entries(response.transmission).into_iter().zip(transmission))
        {
            close(actual, expected, 2.0e-12);
        }
        close(response.scalar.reflectance, reflection[0], 2.0e-12);
        close(response.scalar.transmittance, transmission[0], 2.0e-12);
    }
}

#[test]
fn g1_quarter_wave_antireflection_changes_native_event_probabilities() {
    let coated = surface(&[layer(1.5, 550.0 / 6.0)], 1);
    let interface = coated.interface(550.0, false).unwrap();
    let result = interface.fresnel(1.0, 1.0, 2.25).unwrap();
    assert!(result.reflectance < 1.0e-28);
    close(result.transmittance, 1.0, 2.0e-15);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let bare =
        crate::dielectric::sample_smooth_dielectric(normal, normal, 1.0, 2.25, 0.05).unwrap();
    let actual = interface
        .sample_smooth(normal, normal, 1.0, 2.25, 0.05)
        .unwrap();
    assert_eq!(bare.event, DielectricEvent::Reflection);
    assert_eq!(actual.event, DielectricEvent::Transmission);
    assert_eq!(actual.direction, normal.scale(-1.0));
    close(actual.probability, 1.0, 2.0e-15);
    close(actual.radiance_weight, (1.0_f64 / 2.25).powi(2), 2.0e-15);
    let exit = coated
        .interface(550.0, true)
        .unwrap()
        .sample_smooth(normal, normal, 2.25, 1.0, 0.5)
        .unwrap();
    close(actual.radiance_weight * exit.radiance_weight, 1.0, 2.0e-15);
}

#[test]
fn g1_equal_bulk_indices_keep_coherent_etalon_reflection_and_exact_nulls() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let etalon = surface(&[layer(1.5, 550.0 / 6.0)], 1);
    let interface = etalon.interface(550.0, false).unwrap();
    let f = interface.fresnel(1.0, 1.0, 1.0).unwrap();
    let expected = ((2.25_f64 - 1.0) / (2.25 + 1.0)).powi(2);
    close(f.reflectance, expected, 2.0e-14);
    close(f.transmittance, 1.0 - expected, 2.0e-14);
    assert_eq!(
        interface
            .sample_smooth(normal, normal, 1.0, 1.0, 0.01)
            .unwrap()
            .event,
        DielectricEvent::Reflection
    );
    assert_eq!(
        interface
            .sample_smooth(normal, normal, 1.0, 1.0, 0.9)
            .unwrap()
            .direction,
        normal.scale(-1.0)
    );
    assert_eq!(
        interface.sample_rough(normal, normal, 1.0, 1.0, 0.2, 0.2, 0.3, 0.9),
        Err(DielectricError::CoatedNullRoughInterface)
    );
    for trivial in [
        surface(&[layer(1.5, 0.0)], 1),
        surface(&[layer(1.0, 200.0)], 1),
    ] {
        let null = trivial
            .interface(550.0, false)
            .unwrap()
            .sample_rough(normal, normal, 1.0, 1.0, 0.2, 0.2, 0.3, 0.9)
            .unwrap()
            .unwrap();
        assert!(null.delta);
        assert_eq!(null.direction, normal.scale(-1.0));
        assert_eq!(null.radiance_weight, 1.0);
    }
    let half_wave = surface(&[layer(1.5, 550.0 / 3.0)], 1);
    assert!(
        half_wave
            .interface(550.0, false)
            .unwrap()
            .fresnel(1.0, 1.0, 1.0)
            .unwrap()
            .reflectance
            < 1.0e-28
    );
}

#[test]
fn g1_internal_critical_and_evanescent_films_transmit_correct_power() {
    let film = surface(&[layer(1.0, 95.0)], 1);
    let response = film
        .interface(550.0, false)
        .unwrap()
        .fresnel_polarized(0.5, 1.5, 1.5)
        .unwrap();
    assert!(!response.scalar.total_internal_reflection);
    close(
        response.scalar.reflectance,
        0.601_431_462_981_405_3,
        2.0e-12,
    );
    close(
        response.scalar.transmittance,
        0.398_568_537_018_594_8,
        2.0e-12,
    );
    close(
        response.transmission.0[3][2],
        0.209_894_669_844_168_02,
        2.0e-12,
    );
    // At the film's exact critical angle, its characteristic layer has the
    // analytic form [1,-i*a;0,1]. Equal external media give T=4/(4+(a*Y)^2).
    let critical = (1.0_f64 - 1.0 / 2.25).sqrt();
    let f = film
        .interface(550.0, false)
        .unwrap()
        .fresnel(critical, 1.5, 1.5)
        .unwrap();
    let a = core::f64::consts::TAU * 95.0 / 550.0;
    let ys = 1.5 * critical;
    let yp = critical / 1.5;
    let expected = 0.5 * (4.0 / (4.0 + (a * ys).powi(2)) + 4.0 / (4.0 + (a * yp).powi(2)));
    close(f.transmittance, expected, 2.0e-13);
    for adjacent in [critical - 1.0e-9, critical + 1.0e-9] {
        close(
            film.interface(550.0, false)
                .unwrap()
                .fresnel(adjacent, 1.5, 1.5)
                .unwrap()
                .transmittance,
            expected,
            2.0e-8,
        );
    }
}

#[test]
fn g1_external_tir_retains_stack_phase_and_zero_transmitted_flux() {
    let coated = surface(&[layer(1.3, 72.0), layer(2.1, 43.0), layer(1.45, 110.0)], 1);
    let f = coated
        .interface(550.0, false)
        .unwrap()
        .fresnel_polarized(0.5, 1.5, 1.0)
        .unwrap();
    assert!(f.scalar.total_internal_reflection);
    close(f.scalar.reflectance, 1.0, 2.0e-15);
    assert_eq!(f.scalar.transmittance, 0.0);
    assert_eq!(f.transmission, Mueller([[0.0; 4]; 4]));
    close(f.reflection.0[2][2], -0.513_595_207_316_458_6, 2.0e-12);
    close(f.reflection.0[3][2], 0.858_032_611_863_654_4, 2.0e-12);
}

#[test]
fn g3_lossless_flux_and_front_back_reciprocity_hold_for_all_stokes_inputs() {
    let coated = surface(&[layer(1.2, 48.0), layer(2.4, 75.0), layer(1.6, 120.0)], 7);
    for wavelength in [380.0, 550.0, 780.0] {
        for cosine in [0.12, 0.5, 0.91, 1.0] {
            let forward = coated
                .interface(wavelength, false)
                .unwrap()
                .fresnel_polarized(cosine, 1.0, 1.7)
                .unwrap();
            let reverse = coated
                .interface(wavelength, true)
                .unwrap()
                .fresnel_polarized(forward.scalar.transmitted_cosine, 1.7, 1.0)
                .unwrap();
            close(
                forward.scalar.reflectance,
                reverse.scalar.reflectance,
                2.0e-12,
            );
            close(
                forward.scalar.transmittance,
                reverse.scalar.transmittance,
                2.0e-12,
            );
            for state in [
                Stokes::unpolarized(1.0),
                Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap(),
                Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap(),
                Stokes::try_new(1.0, 0.0, 0.6, 0.8).unwrap(),
            ] {
                let r = forward.reflection.apply(state);
                let t = forward.transmission.apply(state);
                close(r.i + t.i, 1.0, 2.0e-14);
                assert!(Stokes::try_new(r.i, r.q, r.u, r.v).is_some());
                assert!(Stokes::try_new(t.i, t.q, t.u, t.v).is_some());
            }
        }
    }
}

#[test]
fn g1_deep_bragg_stack_retains_transmission_below_one_minus_r_resolution() {
    let coated = surface(&[layer(3.0, 550.0 / 12.0), layer(1.0, 550.0 / 4.0)], 32);
    let f = coated
        .interface(550.0, false)
        .unwrap()
        .fresnel(1.0, 1.5, 1.5)
        .unwrap();
    let ratio = 3.0_f64.powi(64);
    let expected = 4.0 * ratio / (1.0 + ratio).powi(2);
    assert_eq!(f.reflectance, 1.0);
    assert!(f.transmittance > 0.0);
    close(f.transmittance / expected, 1.0, 2.0e-12);
}

#[test]
fn g1_weak_transmission_keeps_smooth_and_rough_sampling_support() {
    let coated = surface(&[layer(3.0, 550.0 / 12.0), layer(1.0, 550.0 / 4.0)], 32);
    let interface = coated.interface(550.0, false).unwrap();
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let f = interface.fresnel(1.0, 1.0, 1.5).unwrap();
    assert_eq!(f.reflectance, 1.0);
    assert!(f.transmittance > 0.0 && f.transmittance < f64::EPSILON);
    let probability = f.transmittance / (f.reflectance + f.transmittance);
    for event in [0.0, probability * 0.5] {
        let smooth = interface
            .sample_smooth(normal, normal, 1.0, 1.5, event)
            .unwrap();
        assert_eq!(smooth.event, DielectricEvent::Transmission);
        assert_eq!(smooth.probability, probability);
        close(
            smooth.probability * smooth.radiance_weight / f.transmittance,
            (1.0_f64 / 1.5).powi(2),
            2.0e-14,
        );
        let rough = interface
            .sample_rough(normal, normal, 1.0, 1.5, 0.2, 0.0, 0.0, event)
            .unwrap()
            .unwrap();
        assert_eq!(rough.event, DielectricEvent::Transmission);
        assert!(rough.pdf > 0.0 && rough.value > 0.0);
        let evaluation = interface
            .evaluate_rough(normal, normal, rough.direction, 1.0, 1.5, 0.2)
            .unwrap();
        assert_eq!(rough.pdf, evaluation.pdf);
        close(
            rough.radiance_weight,
            evaluation.value / evaluation.pdf,
            2.0e-14,
        );
        let polarized = interface
            .sample_smooth_polarized(normal, normal, 1.0, 1.5, event)
            .unwrap();
        close(
            polarized.radiance_weight.0[0][0],
            smooth.radiance_weight,
            2.0e-14,
        );
    }
    // The half-open threshold belongs to reflection, with no missing interval
    // or zero-probability selected event. Both native samplers make the same choice.
    assert_eq!(
        interface
            .sample_smooth(normal, normal, 1.0, 1.5, probability)
            .unwrap()
            .event,
        DielectricEvent::Reflection
    );
    assert_eq!(
        interface
            .sample_rough(normal, normal, 1.0, 1.5, 0.2, 0.0, 0.0, probability)
            .unwrap()
            .unwrap()
            .event,
        DielectricEvent::Reflection
    );
}

#[test]
fn g1_exit_critical_roundoff_uses_native_snell_classification() {
    // The ratio-based sin² calculation rounds to one at the original cosine,
    // but the exact binary64 inputs give positive transmitted normal squared
    // 3.7434066644227618394460881557844072642e-16. Retain this physical
    // propagating branch and the adjacent binary64 rays on either side.
    let (eta_i, eta_t, cosine) = (
        2.048_747_980_327_861,
        1.846_015_323_252_327_5,
        0.433_724_397_570_750_8_f64,
    );
    let coated = surface(&[layer(2.4, 100.0)], 1);
    let interface = coated.interface(550.0, false).unwrap();
    let evanescent_cosine = f64::from_bits(cosine.to_bits() - 2);
    let bare = crate::dielectric::fresnel_dielectric(evanescent_cosine, eta_i, eta_t).unwrap();
    assert!(bare.total_internal_reflection);
    let evanescent = interface
        .fresnel_polarized(evanescent_cosine, eta_i, eta_t)
        .unwrap();
    assert!(evanescent.scalar.total_internal_reflection);
    assert_eq!(
        evanescent.scalar.transmitted_cosine,
        bare.transmitted_cosine
    );
    assert_eq!(evanescent.scalar.transmittance, 0.0);
    assert_eq!(evanescent.transmission, Mueller([[0.0; 4]; 4]));
    close(evanescent.scalar.reflectance, 1.0, 2.0e-15);
    // Independent 90-digit decimal evaluation from exact binary64 inputs.
    // The first entry is adjacent to the evanescent cosine above.
    for (propagating_cosine, transmitted_cosine) in [
        (
            f64::from_bits(cosine.to_bits() - 1),
            7.109_058_727_775_48e-9,
        ),
        (cosine, 1.048_089_096_473_242_3e-8),
        (
            f64::from_bits(cosine.to_bits() + 3),
            1.696_408_422_737_746_4e-8,
        ),
    ] {
        let bare = crate::dielectric::fresnel_dielectric(propagating_cosine, eta_i, eta_t).unwrap();
        assert!(!bare.total_internal_reflection);
        close(bare.transmitted_cosine, transmitted_cosine, 4.0e-23);
        let propagating = interface
            .fresnel_polarized(propagating_cosine, eta_i, eta_t)
            .unwrap();
        assert!(!propagating.scalar.total_internal_reflection);
        assert_eq!(
            propagating.scalar.transmitted_cosine,
            bare.transmitted_cosine
        );
        assert!(propagating.scalar.transmittance > 0.0);
        close(
            propagating.scalar.reflectance + propagating.scalar.transmittance,
            1.0,
            2.0e-14,
        );
    }
}

#[test]
fn g3_borrowed_full_spectrum_laws_and_zero_layers_keep_native_limits() {
    let laws: [crate::dielectric::DielectricIor; 3] = [
        CauchyIor::try_new(1.3, 0.01, 0.0001).unwrap().into(),
        SellmeierIor::try_new(
            [0.696_166_3, 0.407_942_6, 0.897_479_4],
            [
                0.068_404_3_f64.powi(2),
                0.116_241_4_f64.powi(2),
                9.896_161_f64.powi(2),
            ],
        )
        .unwrap()
        .into(),
        TabulatedIor::try_new(core::array::from_fn(|i| {
            DielectricIorSample::try_new(380.0 + 50.0 * i as f64, 1.8 - 0.02 * i as f64).unwrap()
        }))
        .unwrap()
        .into(),
    ];
    let films = core::array::from_fn::<_, 3, _>(|i| {
        CoherentDielectricCoating::try_new(laws[i], (60.0 + 15.0 * i as f64) * 1.0e-9, source())
            .unwrap()
    });
    let stack = DielectricCoatingStack::try_new(&films).unwrap();
    for wavelength in [380.0, 505.0, 780.0] {
        // At this one wavelength, each full law must equal an explicitly
        // authored constant phase index, while retaining its original data.
        let equivalent = core::array::from_fn::<_, 3, _>(|i| {
            layer(laws[i].eval(wavelength).unwrap(), 60.0 + 15.0 * i as f64)
        });
        let compact = surface(&equivalent, 1);
        let expected = compact
            .interface(wavelength, false)
            .unwrap()
            .fresnel(0.8, 1.0, 1.6)
            .unwrap();
        let actual = stack.fresnel(wavelength, false, 0.8, 1.0, 1.6).unwrap();
        close(actual.reflectance, expected.reflectance, 2.0e-14);
        close(actual.transmittance, expected.transmittance, 2.0e-14);
    }
    let zero = [CoherentDielectricCoating::try_new(laws[2], 0.0, source()).unwrap(); 64];
    for layers in [&[][..], &zero[..]] {
        for cosine in [0.0, 0.3, 1.0] {
            let actual = DielectricCoatingStack::try_new(layers)
                .unwrap()
                .fresnel_polarized(550.0, false, cosine, 1.5, 1.0)
                .unwrap();
            assert_eq!(
                actual,
                fresnel_dielectric_polarized(cosine, 1.5, 1.0).unwrap()
            );
        }
    }
}

#[test]
fn g3_rough_coating_reuses_only_identical_angles_and_preserves_exact_pdf() {
    use crate::dielectric::{
        SampledFresnel, evaluate_rough_with_fresnel, normalize_admitted, sample_ggx_visible_normal,
    };

    let absorbing =
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.0).unwrap(), 0.01, 50.0e-9).unwrap();
    let coatings = [
        surface(&[layer(1.3, 90.0), layer(2.0, 50.0)], 32),
        surface(&[layer(1.3, 90.0), absorbing], 32),
    ];
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let mut accepted = 0_u32;
    let mut reused = 0_u32;
    let mut reflections = 0_u32;
    let mut transmissions = 0_u32;
    let mut tir = 0_u32;
    for coated in coatings {
        for (eta_i, eta_t, reverse) in [(1.0, 1.5, false), (1.5, 1.0, true)] {
            let interface = coated.interface(550.0, reverse).unwrap();
            for cosine in [1.0_f64, 0.8, 0.01] {
                let wo = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, cosine);
                let normalized_wo = normalize_admitted(wo);
                for alpha in [0.03, 0.5] {
                    for (u, v) in [(0.0, 0.0), (0.17, 0.63), (0.37, 0.61), (0.99, 0.1)] {
                        let facet =
                            sample_ggx_visible_normal(normal, normalized_wo, alpha, u, v).unwrap();
                        let sampled_cosine = normalized_wo.dot(facet).clamp(0.0, 1.0);
                        let response = interface.fresnel(sampled_cosine, eta_i, eta_t).unwrap();
                        for event_sample in [0.0, 0.999] {
                            let Some(sample) = interface
                                .sample_rough(normal, wo, eta_i, eta_t, alpha, u, v, event_sample)
                                .unwrap()
                            else {
                                continue;
                            };
                            let full = interface
                                .evaluate_rough(
                                    normal,
                                    normalized_wo,
                                    sample.direction,
                                    eta_i,
                                    eta_t,
                                    alpha,
                                )
                                .unwrap();
                            let cached = SampledFresnel {
                                cosine: sampled_cosine,
                                response,
                            };
                            let mut recomputed = false;
                            let replay = evaluate_rough_with_fresnel(
                                normal,
                                normalized_wo,
                                sample.direction,
                                eta_i,
                                eta_t,
                                alpha,
                                Some(interface),
                                |angle| {
                                    cached.resolve(angle, || {
                                        recomputed = true;
                                        interface.fresnel(angle, eta_i, eta_t)
                                    })
                                },
                            )
                            .unwrap();
                            assert_eq!(replay.value.to_bits(), full.value.to_bits());
                            assert_eq!(replay.pdf.to_bits(), full.pdf.to_bits());
                            assert_eq!(sample.value.to_bits(), full.value.to_bits());
                            assert_eq!(sample.pdf.to_bits(), full.pdf.to_bits());
                            assert_eq!(
                                sample.radiance_weight.to_bits(),
                                (full.value * normal.dot(sample.direction).abs() / full.pdf)
                                    .to_bits()
                            );
                            accepted += 1;
                            reused += u32::from(!recomputed);
                            reflections += u32::from(sample.event == DielectricEvent::Reflection);
                            transmissions +=
                                u32::from(sample.event == DielectricEvent::Transmission);
                            tir += u32::from(response.total_internal_reflection);
                        }
                    }
                }
            }
        }
    }
    assert!(reflections > 0 && transmissions > 0 && tir > 0);
    assert!(reused > 0 && reused < accepted);
    eprintln!(
        "rough 64-layer coatings: {accepted} accepted draws, {reused} exact response reuses; optical solves {} -> {}",
        2 * accepted,
        2 * accepted - reused
    );
}

#[test]
fn g0_rough_coating_response_recomputes_for_one_ulp_angle_change() {
    let coated = surface(&[layer(3.0, 550.0 / 12.0), layer(1.0, 550.0 / 4.0)], 32);
    let interface = coated.interface(550.0, false).unwrap();
    let cosine = 0.8_f64;
    let response = interface.fresnel(cosine, 1.0, 1.5).unwrap();
    let cache = crate::dielectric::SampledFresnel { cosine, response };
    assert_eq!(
        cache.resolve(cosine, || panic!(
            "identical angle repeated the coating solve"
        )),
        Ok(response)
    );
    let adjacent = f64::from_bits(cosine.to_bits() + 1);
    let mut calls = 0;
    let cache = crate::dielectric::SampledFresnel { cosine, response };
    let actual = cache
        .resolve(adjacent, || {
            calls += 1;
            interface.fresnel(adjacent, 1.0, 1.5)
        })
        .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(actual, interface.fresnel(adjacent, 1.0, 1.5).unwrap());
}

#[test]
fn g3_rough_samples_replay_pdf_and_reciprocal_transmitted_bsdf() {
    let coated = surface(&[layer(1.3, 90.0), layer(2.0, 50.0)], 2);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    let interface = coated.interface(550.0, false).unwrap();
    for event_sample in [0.0, 0.999] {
        let sample = interface
            .sample_rough(normal, wo, 1.0, 1.6, 0.17, 0.21, 0.67, event_sample)
            .unwrap()
            .unwrap();
        let evaluation = interface
            .evaluate_rough(normal, wo, sample.direction, 1.0, 1.6, 0.17)
            .unwrap();
        assert_eq!(sample.pdf, evaluation.pdf);
        assert_eq!(sample.value, evaluation.value);
        close(
            sample.radiance_weight,
            evaluation.value * normal.dot(sample.direction).abs() / sample.pdf,
            2.0e-14,
        );
        if sample.event == DielectricEvent::Transmission {
            let reverse = coated
                .interface(550.0, true)
                .unwrap()
                .evaluate_rough(normal.scale(-1.0), sample.direction, wo, 1.6, 1.0, 0.17)
                .unwrap();
            close(evaluation.value * 1.6_f64.powi(2), reverse.value, 2.0e-10);
        }
    }
}

#[test]
fn g1_native_polarized_beam_uses_coated_transmission_and_correct_frame() {
    let coated = surface(&[layer(1.3, 72.0), layer(2.1, 43.0), layer(1.45, 110.0)], 1);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let cosine = 0.83_f64;
    let direction = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, -cosine);
    let input = Stokes::try_new(1.0, 0.0, 0.6, 0.8).unwrap();
    let beam = PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(direction, Vec3::new(0.0, -1.0, 0.0)).unwrap(),
        550.0,
        glass(1.0),
        input,
    )
    .unwrap();
    let expected = coated
        .interface(550.0, false)
        .unwrap()
        .fresnel_polarized(cosine, 1.0, 1.7)
        .unwrap();
    for (event_sample, event, matrix, probability) in [
        (
            0.0,
            DielectricEvent::Reflection,
            expected.reflection,
            expected.scalar.reflectance,
        ),
        (
            0.999,
            DielectricEvent::Transmission,
            expected.transmission,
            expected.scalar.transmittance,
        ),
    ] {
        let sampled = beam
            .sample_interface(normal, glass(1.7), coated, [0.2, 0.7, event_sample])
            .unwrap()
            .unwrap();
        assert_eq!(sampled.event, event);
        close(sampled.pdf, probability, 2.0e-14);
        let projected = matrix.apply(input);
        let actual = sampled.beam.power();
        for (actual, expected) in [actual.i, actual.q, actual.u, actual.v].into_iter().zip([
            projected.i,
            projected.q,
            projected.u,
            projected.v,
        ]) {
            close(actual * sampled.pdf, expected, 2.0e-12);
        }
    }
}

#[test]
fn g0_absorbing_dielectric_films_retain_loss_and_work_is_bounded() {
    let absorbing =
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.5).unwrap(), 0.01, 80.0e-9).unwrap();
    let coating = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(&[absorbing], 1, source()).unwrap(),
    )
    .unwrap();
    let response = DielectricSurface::SMOOTH
        .with_coating(coating)
        .interface(550.0, false)
        .unwrap()
        .fresnel(0.8, 1.0, 1.5)
        .unwrap();
    assert!(response.absorption_residual() > 0.0);
    let film = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.5).unwrap(),
        80.0e-9,
        source(),
    )
    .unwrap();
    assert!(DielectricCoatingStack::try_new(&[film; 65]).is_err());
    let coated = surface(&[layer(1.4, 80.0)], 1);
    assert!(
        coated
            .interface(550.0, false)
            .unwrap()
            .fresnel(0.0, 1.0, 1.5)
            .is_err()
    );
    assert_eq!(
        core::mem::size_of::<DielectricCoating>(),
        core::mem::size_of::<CoherentPeriodicCoating>()
    );
}

#[test]
fn g1_lossless_closure_allowance_is_derived_for_the_solve_and_refuses_planted_excess() {
    // g3's lossless 1.3/2.0 x 32 stack between 1.0 and 1.5 at 550 nm. A
    // lossless solve is checked against its own derived first-order allowance
    // (multilayer.rs ClosureBound), not the flat 512 eps (N + 1).
    let coating = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(&[layer(1.3, 90.0), layer(2.0, 50.0)], 32, source())
            .unwrap(),
    )
    .unwrap();
    let layers = coating.0.layer_count();
    let flat = 512.0 * f64::EPSILON * (layers + 1) as f64;
    let solve = |cosine: f64| {
        let bare = fresnel_dielectric_polarized(cosine, 1.0, 1.5).unwrap();
        let exit = exit_normal_wave(bare.scalar, cosine, 1.0, 1.5);
        (coating.0.transmission(550.0, 1.0, 1.5, cosine, false, exit).unwrap(), bare)
    };
    // Rescale pol 0's transmitted amplitude so that its closure is 1 + error.
    let planted = |cosine: f64, error: f64| {
        let (mut response, bare) = solve(cosine);
        let [r, t] = [response.reflected[0], response.transmitted[0]];
        let reflected = r[0] * r[0] + r[1] * r[1];
        let scale = ((1.0 + error - reflected) / (t[0] * t[0] + t[1] * t[1])).sqrt();
        response.transmitted[0] = t.map(|component| component * scale);
        (response, bare)
    };
    // The facet cosine of g3's refused call. Its closure error, 7.99e-12, was
    // above the flat allowance; the derived allowance there is larger and
    // admits it. Planted errors of twice the allowance are refused on either
    // side; half of it is admitted.
    let grazing = 4.39970846435952696e-3;
    let (response, bare) = solve(grazing);
    assert!(response.lossless && !response.bare_equivalent);
    assert!(response.closure_allowance.iter().all(Option::is_some));
    let allowance = response.closure_allowance[0].unwrap();
    assert!(allowance > flat, "{allowance:e}");
    assert!(polarized_response(response, bare, layers).is_ok());
    for (factor, admitted) in [(2.0, false), (-2.0, false), (0.5, true), (-0.5, true)] {
        let (response, bare) = planted(grazing, factor * allowance);
        assert_eq!(
            polarized_response(response, bare, layers).is_ok(),
            admitted,
            "{factor} x {allowance:e}"
        );
    }
    // Normal incidence: the derived allowance is below 1e-12 and the flat one
    // above it. A planted 1e-12 is refused, and admitted on the flat path.
    let (response, _) = solve(1.0);
    let allowance = response.closure_allowance[0].unwrap();
    assert!(allowance < 1.0e-12 && 1.0e-12 < flat, "{allowance:e}");
    let (response, bare) = planted(1.0, 1.0e-12);
    assert!(matches!(
        polarized_response(response, bare, layers),
        Err(DielectricError::InvalidInterface)
    ));
    let (mut response, bare) = planted(1.0, 1.0e-12);
    response.closure_allowance = [None; 2];
    assert!(polarized_response(response, bare, layers).is_ok());
    // A deep reflector leaves the first-order regime. The quarter-wave
    // 3.0/1.0 x 32 stack between 1.5 and 1.5 has T ~ 1.2e-30 and reaches
    // u kappa ~ 1e14, where the stored J has no resolution left. Its
    // allowance is None, and the flat check, the one applied before the
    // derived allowance existed, admits it.
    let deep = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[layer(3.0, 550.0 / 12.0), layer(1.0, 550.0 / 4.0)],
            32,
            source(),
        )
        .unwrap(),
    )
    .unwrap();
    let bare = fresnel_dielectric_polarized(1.0, 1.5, 1.5).unwrap();
    let exit = exit_normal_wave(bare.scalar, 1.0, 1.5, 1.5);
    let response = deep.0.transmission(550.0, 1.5, 1.5, 1.0, false, exit).unwrap();
    assert!(response.lossless && response.closure_allowance == [None; 2]);
    assert!(polarized_response(response, bare, deep.0.layer_count()).is_ok());
}
