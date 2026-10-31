use super::*;
use crate::conductor::{ConductorDataStatus, ConductorSource};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIorSample, GlassProvenance,
    TabulatedIor,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{BeamMedium, PolarizationFrame, PolarizedBeam};
use fs_blake3::ContentHash;
use fs_geom::Point3;
use fs_math::det;

fn film(index: f64, thickness: f64) -> CoherentDielectricCoating {
    CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(index).unwrap(),
        thickness,
        ConductorSource::try_new(ContentHash([104; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}
fn glass(index: f64, extinction: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(extinction).unwrap(),
        GlassProvenance::Custom,
    )
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12, "{a:.17e} != {b:.17e}");
}
fn terms(m: Mueller) -> [f64; 4] {
    [m.0[0][0], m.0[0][1], m.0[2][2], m.0[3][2]]
}
fn same(a: Mueller, b: Mueller) {
    for (a, b) in a.0.into_iter().flatten().zip(b.0.into_iter().flatten()) {
        near(a, b);
    }
}

#[test]
fn g1_complex_coating_matrices_match_independent_airy_answers() {
    // Independent complex Airy multiple-reflection calculation. Both output
    // matrices retain the signed cross phase, not merely scalar reflectance.
    for (ni, nf, nt, c, d, r, t) in [
        (
            1.0,
            2.0,
            1.5,
            0.6,
            80e-9,
            [
                0.22490029475325646,
                0.17065949664811386,
                -0.14628927640299488,
                0.007411234269822986,
            ],
            [
                0.7750997052467434,
                -0.17065949664811397,
                0.7559969556627637,
                -0.011112709212171737,
            ],
        ),
        (
            1.5,
            1.0,
            1.33,
            0.7,
            180e-9,
            [
                0.522439790490337,
                0.08161248024603973,
                0.43363568221384163,
                0.27971920360714086,
            ],
            [
                0.47756020950966305,
                -0.08161248024603987,
                0.4232654225862848,
                0.20554692607250713,
            ],
        ),
        (
            1.5,
            1.2,
            1.0,
            0.5,
            100e-9,
            [1.0, 0.0, 0.8554977493746024, 0.5178065283626596],
            [0.0; 4],
        ),
        (
            1.0,
            1.7,
            1.0,
            1.0,
            550e-9 / (4.0 * 1.7),
            [0.23606108867903322, 0.0, -0.23606108867903322, 0.0],
            [0.7639389113209665, 0.0, 0.7639389113209665, 0.0],
        ),
    ] {
        let f = fresnel_coated_dielectric_polarized(c, ni, nt, 550.0, &film(nf, d)).unwrap();
        for (a, b) in terms(f.reflection).into_iter().zip(r) {
            near(a, b);
        }
        for (a, b) in terms(f.transmission).into_iter().zip(t) {
            near(a, b);
        }
        for (q, u, v) in [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (-1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.3, -0.4, 0.2),
        ] {
            let input = Stokes::try_new(1.0, q, u, v).unwrap();
            let reflected = f.reflection.apply(input);
            let transmitted = f.transmission.apply(input);
            near(reflected.i + transmitted.i, 1.0);
            for s in [reflected, transmitted] {
                assert!(s.i >= -1e-14 && s.i <= 1.0 + 1e-14);
                assert!(s.q.hypot(s.u).hypot(s.v) <= s.i + 1e-14);
            }
        }
    }
}

#[test]
fn g3_zero_matched_antireflection_and_equal_surrounding_limits_are_distinct() {
    for (ni, nt) in [(1.0, 1.5), (1.5, 1.0), (1.33, 1.33)] {
        for c in [0.0, 0.3, 0.7, 1.0] {
            let bare = fresnel_dielectric_polarized(c, ni, nt).unwrap();
            for coating in [film(1.7, 0.0), film(ni, 100e-9), film(nt, 100e-9)] {
                assert_eq!(
                    fresnel_coated_dielectric_polarized(c, ni, nt, 550.0, &coating).unwrap(),
                    bare
                );
            }
        }
    }
    let n = det::sqrt(1.5);
    let ar = film(n, 550e-9 / (4.0 * n));
    let f = fresnel_coated_dielectric_polarized(1.0, 1.0, 1.5, 550.0, &ar).unwrap();
    assert!(f.scalar.reflectance < 1e-28);
    same(f.transmission, Mueller::IDENTITY);
    let half = film(1.7, 550e-9 / (2.0 * 1.7));
    let f = fresnel_coated_dielectric_polarized(1.0, 1.0, 1.0, 550.0, &half).unwrap();
    assert!(f.scalar.reflectance < 1e-28);
    same(f.transmission, Mueller::IDENTITY);
    let quarter = film(1.7, 550e-9 / (4.0 * 1.7));
    let f = fresnel_coated_dielectric_polarized(1.0, 1.0, 1.0, 550.0, &quarter).unwrap();
    assert!(f.scalar.reflectance > 0.23);
    assert!(!f.scalar.total_internal_reflection);
    let n = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.6, 0.0, 0.8);
    for u in [0.0, 0.99] {
        assert_eq!(
            sample_smooth_coated_dielectric_polarized(n, wo, 1.0, 1.5, 550.0, &film(1.7, 0.0), u)
                .unwrap(),
            sample_smooth_dielectric_polarized(n, wo, 1.0, 1.5, u).unwrap()
        );
    }
}

#[test]
fn g1_critical_and_evanescent_film_keeps_transmission_and_tir_phase() {
    let coating = film(1.0, 200e-9);
    let critical = fresnel_coated_dielectric_polarized(0.6, 1.25, 1.5, 550.0, &coating).unwrap();
    // Independent q_f=0 characteristic matrix: diagonal=1,
    // upper=-i*(2*pi*d/lambda), lower=0, for n_f=1.
    for (a, b) in terms(critical.reflection).into_iter().zip([
        0.3845567348348546,
        0.14699752361403012,
        0.35233505364479833,
        0.046212772382798695,
    ]) {
        near(a, b);
    }
    for (a, b) in terms(critical.transmission).into_iter().zip([
        0.6154432651651457,
        -0.14699752361403007,
        0.5728131539054875,
        0.17043248341007883,
    ]) {
        near(a, b);
    }
    for dc in [-1e-9, 1e-9] {
        let f = fresnel_coated_dielectric_polarized(0.6 + dc, 1.25, 1.5, 550.0, &coating).unwrap();
        for (a, b) in f
            .transmission
            .0
            .iter()
            .flatten()
            .zip(critical.transmission.0.iter().flatten())
        {
            assert!((a - b).abs() < 1e-7);
        }
    }
    let tunnel =
        fresnel_coated_dielectric_polarized(0.7, 1.5, 1.33, 550.0, &film(1.0, 180e-9)).unwrap();
    assert!(tunnel.transmission.0[0][0] > 0.47);
    let opaque =
        fresnel_coated_dielectric_polarized(0.7, 1.5, 1.33, 550.0, &film(1.0, 1e-4)).unwrap();
    same(
        opaque.reflection,
        fresnel_dielectric_polarized(0.7, 1.5, 1.0)
            .unwrap()
            .reflection,
    );
    assert_eq!(opaque.transmission, Mueller([[0.0; 4]; 4]));
    let tir =
        fresnel_coated_dielectric_polarized(0.5, 1.5, 1.0, 550.0, &film(1.2, 100e-9)).unwrap();
    assert!(tir.scalar.total_internal_reflection);
    near(tir.reflection.0[0][0], 1.0);
    assert!(tir.reflection.0[3][2] > 0.5);
    assert_eq!(tir.transmission, Mueller([[0.0; 4]; 4]));
}

#[test]
fn g3_original_film_is_reciprocal_when_incident_and_exit_media_reverse() {
    for (ni, nf, nt, c) in [
        (1.0, 2.0, 1.5, 0.6),
        (1.5, 1.0, 1.33, 0.7),
        (1.0, 1.7, 1.0, 0.2),
    ] {
        let coating = film(nf, 180e-9);
        let forward = fresnel_coated_dielectric_polarized(c, ni, nt, 550.0, &coating).unwrap();
        let reverse = fresnel_coated_dielectric_polarized(
            forward.scalar.transmitted_cosine,
            nt,
            ni,
            550.0,
            &coating,
        )
        .unwrap();
        same(forward.transmission, reverse.transmission);
        near(forward.reflection.0[0][0], reverse.reflection.0[0][0]);
        near(forward.reflection.0[0][1], reverse.reflection.0[0][1]);
    }
}

#[test]
fn g1_beam_samples_close_polarized_branch_work_and_change_only_transmitted_medium() {
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let coating = film(2.0, 80e-9);
    for (ni, nt) in [(1.0, 1.5), (1.5, 1.0), (1.33, 1.33)] {
        let incident = glass(ni, 0.0);
        let target = glass(nt, 0.25);
        let input = Stokes::try_new(1.0, 0.3, 0.4, 0.2).unwrap();
        let beam = PolarizedBeam::try_new(
            Point3::new(0.0, 0.0, 0.0),
            PolarizationFrame::try_new(Vec3::new(0.6, 0.0, -0.8), Vec3::new(0.0, -1.0, 0.0))
                .unwrap(),
            550.0,
            incident,
            input,
        )
        .unwrap();
        let expected = fresnel_coated_dielectric_polarized(0.8, ni, nt, 550.0, &coating).unwrap();
        let mut power = 0.0;
        let mut reflected = false;
        let mut transmitted = false;
        for event_sample in [0.0, 0.999999] {
            let out = beam
                .sample_coated_interface(normal, target, &coating, event_sample)
                .unwrap();
            let (matrix, medium) = match out.event {
                DielectricEvent::Reflection => {
                    reflected = true;
                    (expected.reflection, BeamMedium::Homogeneous(incident))
                }
                DielectricEvent::Transmission => {
                    transmitted = true;
                    (expected.transmission, BeamMedium::Homogeneous(target))
                }
            };
            let reference = matrix.apply(input);
            let actual = out.beam.power();
            for (a, b) in [actual.i, actual.q, actual.u, actual.v].into_iter().zip([
                reference.i,
                reference.q,
                reference.u,
                reference.v,
            ]) {
                near(a * out.pdf, b);
            }
            assert_eq!(out.beam.medium(), medium);
            assert_eq!(out.beam.ray().origin, beam.ray().origin);
            assert!(out.delta);
            near(out.beam.ray().dir.norm(), 1.0);
            power += actual.i * out.pdf;
        }
        assert!(reflected && transmitted);
        near(power, input.i);
    }
}

#[test]
fn g1_positive_weak_branches_survive_a_rounded_unit_complement() {
    let n = Vec3::new(0.0, 0.0, 1.0);
    let wo = Vec3::new(0.8, 0.0, 0.6);
    let opaque = film(1.0, 1e-5);
    let f = fresnel_coated_dielectric_polarized(0.6, 1.5, 1.5, 550.0, &opaque).unwrap();
    assert_eq!(1.0 - f.scalar.reflectance, 0.0);
    assert!(f.transmission.0[0][0] > 0.0 && f.transmission.0[0][0] < f64::EPSILON);
    let t =
        sample_smooth_coated_dielectric_polarized(n, wo, 1.5, 1.5, 550.0, &opaque, 0.0).unwrap();
    assert_eq!(t.event, DielectricEvent::Transmission);
    assert!(t.probability > 0.0 && t.probability < f64::EPSILON);
    near(t.radiance_weight.apply(Stokes::unpolarized(1.0)).i, 1.0);
    let index = det::sqrt(1.5);
    let ar = film(index, 550e-9 / (4.0 * index));
    let r = sample_smooth_coated_dielectric_polarized(n, n, 1.0, 1.5, 550.0, &ar, 0.0).unwrap();
    assert_eq!(r.event, DielectricEvent::Reflection);
    assert!(r.probability > 0.0 && r.probability < f64::EPSILON);
    near(r.radiance_weight.apply(Stokes::unpolarized(1.0)).i, 1.0);
}

#[test]
fn g3_sampled_original_film_law_drives_interference_without_a_dispersion_fit() {
    let law = TabulatedIor::try_new(core::array::from_fn(|i| {
        DielectricIorSample::try_new(380.0 + 50.0 * i as f64, 1.2 + 0.05 * i as f64).unwrap()
    }))
    .unwrap();
    let source = film(1.5, 100e-9).source();
    let coating = CoherentDielectricCoating::try_new(law, 100e-9, source).unwrap();
    assert_eq!(coating.source(), source);
    for i in 0..9 {
        let w = 380.0 + 50.0 * i as f64;
        assert_eq!(
            fresnel_coated_dielectric_polarized(0.8, 1.0, 1.5, w, &coating).unwrap(),
            fresnel_coated_dielectric_polarized(
                0.8,
                1.0,
                1.5,
                w,
                &film(1.2 + 0.05 * i as f64, 100e-9)
            )
            .unwrap()
        );
    }
    assert!(fresnel_coated_dielectric_polarized(0.8, 1.0, 1.5, 379.0, &coating).is_err());
    assert!(fresnel_coated_dielectric_polarized(0.0, 1.0, 1.5, 550.0, &coating).is_err());
    assert!(fresnel_coated_dielectric_polarized(f64::NAN, 1.0, 1.5, 550.0, &coating).is_err());
    assert!(fresnel_coated_dielectric_polarized(0.8, 0.5, 1.5, 550.0, &coating).is_err());
}

#[test]
fn g3_original_film_adapter_shares_explicit_transmission_and_native_critical_root() {
    // Independent quarter-wave antireflection limit. Before the shared-owner
    // repair scalar T incorrectly retained the bare 0.96 while its matrix had 1.
    let index = det::sqrt(1.5);
    let ar = film(index, 550e-9 / (4.0 * index));
    let response = fresnel_coated_dielectric_polarized(1.0, 1.0, 1.5, 550.0, &ar).unwrap();
    assert!(response.scalar.reflectance < 1e-28);
    near(response.scalar.transmittance, 1.0);
    let coating = film(2.4, 100e-9);
    let cosine = 0.433_724_397_570_750_8_f64;
    let (eta_i, eta_t) = (2.048_747_980_327_861, 1.846_015_323_252_327_5);
    for (cosine, tir) in [
        (cosine, true),
        (f64::from_bits(cosine.to_bits() + 3), false),
    ] {
        let actual =
            fresnel_coated_dielectric_polarized(cosine, eta_i, eta_t, 550.0, &coating).unwrap();
        let expected = DielectricCoatingStack::try_new(core::slice::from_ref(&coating))
            .unwrap()
            .fresnel_polarized(550.0, false, cosine, eta_i, eta_t)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(actual.scalar.total_internal_reflection, tir);
        assert_eq!(actual.scalar.transmittance, actual.transmission.0[0][0]);
        if actual.scalar.total_internal_reflection {
            assert_eq!(actual.scalar.transmittance, 0.0);
        } else {
            assert!(actual.scalar.transmittance > 0.0);
        }
    }
    let barrier = film(1.0, 1e-5);
    let weak = fresnel_coated_dielectric_polarized(0.6, 1.5, 1.5, 550.0, &barrier).unwrap();
    assert_eq!(weak.scalar.reflectance, 1.0);
    assert!(weak.scalar.transmittance > 0.0 && weak.scalar.transmittance < f64::EPSILON);
    assert_eq!(weak.scalar.transmittance, weak.transmission.0[0][0]);
}
