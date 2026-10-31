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
    // Native (n_i/n_t)^2*sin² rounds to one here while n_t²-n_i²*sin²
    // rounds positive. The latter must not invent transmitted flux for a
    // native TIR ray; the signed response uses the same critical zero root.
    let (eta_i, eta_t, cosine) = (
        2.048_747_980_327_861,
        1.846_015_323_252_327_5,
        0.433_724_397_570_750_8_f64,
    );
    let coated = surface(&[layer(2.4, 100.0)], 1);
    let interface = coated.interface(550.0, false).unwrap();
    let bare = crate::dielectric::fresnel_dielectric(cosine, eta_i, eta_t).unwrap();
    assert!(bare.total_internal_reflection);
    assert!(eta_t * eta_t - eta_i * eta_i * (1.0 - cosine * cosine) > 0.0);
    let critical = interface.fresnel_polarized(cosine, eta_i, eta_t).unwrap();
    assert!(critical.scalar.total_internal_reflection);
    assert_eq!(critical.scalar.transmitted_cosine, bare.transmitted_cosine);
    assert_eq!(critical.scalar.transmittance, 0.0);
    assert_eq!(critical.transmission, Mueller([[0.0; 4]; 4]));
    close(critical.scalar.reflectance, 1.0, 2.0e-15);
    // Three adjacent representable cosines reach the first native propagating
    // result in this fixture. Its small positive flux must also be retained.
    let propagating_cosine = f64::from_bits(cosine.to_bits() + 3);
    let bare = crate::dielectric::fresnel_dielectric(propagating_cosine, eta_i, eta_t).unwrap();
    assert!(!bare.total_internal_reflection);
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
fn g0_absorbing_dielectric_films_are_refused_and_work_is_bounded() {
    let absorbing =
        CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.5).unwrap(), 0.01, 80.0e-9).unwrap();
    assert_eq!(
        DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(&[absorbing], 1, source()).unwrap()
        ),
        Err(DielectricError::InvalidCoating)
    );
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
