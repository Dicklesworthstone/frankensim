//! Analytic ghost power/delay, reciprocal film order, and proposal measure.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensCoating, LensCoatingLayer, LensSurface, PhysicalCamera,
};
use crate::charts::Ray;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, CoherentDielectricCoating,
    ConductorDataStatus, ConductorSource,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricCoatingStack, DielectricGlass, GlassProvenance,
    SellmeierIor,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 61,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx);
    });
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn glass(ior: CauchyIor, sigma: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(sigma).unwrap(),
        GlassProvenance::Custom,
    )
}

fn slab(material: DielectricGlass) -> SequentialLens {
    SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.02, 0.0, 0.2, Some(material)).unwrap(),
            LensSurface::try_new(0.05, 0.0, 0.2, None).unwrap(),
        ],
    )
    .unwrap()
}

fn ray(sine: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(sine, 0.0, (1.0 - sine * sine).sqrt()),
    }
}

#[test]
fn slab_ghost_retains_polarized_fresnel_and_triple_beer_and_group_path() {
    let pair = LensGhostPath::try_new(1, 0).unwrap();
    with_cx(|cx| {
        for ior in [
            CauchyIor::try_constant(1.5).unwrap(),
            CauchyIor::try_new(1.5, 0.01, 0.0).unwrap(),
        ] {
            for sigma in [0.0, 7.0] {
                let lens = slab(glass(ior, sigma));
                for wavelength in [450.0, 550.0, 650.0] {
                    let n = ior.eval(wavelength).unwrap();
                    let group = ior.group_index(wavelength).unwrap();
                    for sine in [0.0, 0.6] {
                        let incoming = ray(sine);
                        let ci = incoming.dir.z;
                        let st = sine / n;
                        let ct = (1.0 - st * st).sqrt();
                        let rs = ((ci - n * ct) / (ci + n * ct)).powi(2);
                        let rp = ((n * ci - ct) / (n * ci + ct)).powi(2);
                        let interior = 3.0 * 0.03 / ct;
                        let sensor = 0.02 / ci;
                        let expected = 0.5
                            * ((1.0 - rs).powi(2) * rs.powi(2) + (1.0 - rp).powi(2) * rp.powi(2))
                            * (-sigma * interior).exp();
                        let actual = lens
                            .trace_ghost_ray(cx, incoming, wavelength, pair)
                            .unwrap()
                            .unwrap();
                        close(actual.transmittance, expected, 2e-14);
                        close(actual.distance_m, sensor + interior, 3e-15);
                        close(actual.optical_path_m, sensor + n * interior, 4e-15);
                        close(
                            actual.ray.origin.x,
                            0.02 * sine / ci + 3.0 * 0.03 * st / ct,
                            3e-15,
                        );
                        close(actual.ray.origin.z, 0.05, 3e-15);
                        close(actual.ray.dir.x, sine, 3e-15);
                        close(actual.ray.dir.z, ci, 3e-15);
                        let timed = lens
                            .trace_seeded_path(
                                cx,
                                incoming,
                                wavelength,
                                None,
                                Stokes::unpolarized(1.0),
                                true,
                                Some(pair),
                            )
                            .unwrap()
                            .unwrap();
                        close(
                            timed.group_delay_s.unwrap().to_f64(),
                            (sensor + group * interior) / 299_792_458.0,
                            3e-23,
                        );
                        let direct = lens.trace_ray(cx, incoming, wavelength).unwrap().unwrap();
                        assert!(direct.transmittance + actual.transmittance <= 1.0);
                    }
                }
            }
        }
        let clear = slab(glass(CauchyIor::try_constant(1.5).unwrap(), 0.0));
        close(
            clear
                .trace_ghost_ray(cx, ray(0.0), 550.0, pair)
                .unwrap()
                .unwrap()
                .transmittance,
            0.96_f64.powi(2) * 0.04_f64.powi(2),
            2e-16,
        );
    });
}

#[test]
fn ghost_returns_through_original_spectral_films_in_reverse_order() {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([11; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let transparent = CoherentDielectricCoating::try_new(
        SellmeierIor::try_new_in_band([0.9, 0.05, 0.0], [0.02, 0.09, 0.0], (800.0, 1600.0))
            .unwrap(),
        170e-9,
        source,
    )
    .unwrap();
    let absorbing = CoherentAbsorbingCoating::try_new_in_band(
        std::array::from_fn(|i| {
            CoatingIorSample::try_new_in_band(
                800.0 + 100.0 * i as f64,
                2.1 - 0.03 * i as f64,
                0.04 + 0.01 * i as f64,
            )
            .unwrap()
        }),
        130e-9,
        100.0,
        source,
    )
    .unwrap();
    let layers: [LensCoatingLayer; 2] = [transparent.into(), absorbing.into()];
    let coating = LensCoating::try_new(&layers).unwrap();
    let borrowed = [
        CoatingLayer::Dielectric(&transparent),
        CoatingLayer::Absorbing(&absorbing),
    ];
    let native = DielectricCoatingStack::try_mixed_in_band(&borrowed, (800.0, 1600.0)).unwrap();
    let material = DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (800.0, 1600.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.4).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        1064.0,
        &[
            LensSurface::try_new(0.02, 0.0, 0.2, Some(material)).unwrap(),
            LensSurface::try_new(0.05, 0.0, 0.2, None).unwrap(),
        ],
    )
    .unwrap()
    .with_spectral_coating(0, coating.clone())
    .unwrap()
    .with_spectral_coating(1, coating)
    .unwrap();
    let pair = LensGhostPath::try_new(1, 0).unwrap();
    let incoming = ray(0.6);
    let n = material.ior().eval(1064.0).unwrap();
    let ct = (1.0 - (0.6 / n).powi(2)).sqrt();
    let enter = native
        .fresnel_polarized(1064.0, false, incoming.dir.z, 1.0, n)
        .unwrap();
    let front = native.fresnel_polarized(1064.0, false, ct, n, 1.0).unwrap();
    let rear = native.fresnel_polarized(1064.0, true, ct, n, 1.0).unwrap();
    let seed = Stokes {
        i: 0.6,
        q: 0.1,
        u: 0.2,
        v: 0.25,
    };
    let first_return = front.reflection.apply(enter.transmission.apply(seed));
    let expected = front
        .transmission
        .apply(rear.reflection.apply(first_return));
    let wrong_order = front
        .transmission
        .apply(front.reflection.apply(first_return));
    assert!((expected.i - wrong_order.i).abs() > 1e-5);
    let attenuation = (-0.4 * 3.0 * 0.03 / ct).exp();
    with_cx(|cx| {
        let actual = lens
            .trace_seeded_path(
                cx,
                incoming,
                1064.0,
                Some(Vec3::new(0.0, 1.0, 0.0)),
                seed,
                false,
                Some(pair),
            )
            .unwrap()
            .unwrap();
        for (actual, expected) in [
            actual.power.i,
            actual.power.q,
            actual.power.u,
            actual.power.v,
        ]
        .into_iter()
        .zip([expected.i, expected.q, expected.u, expected.v])
        {
            close(actual, attenuation * expected, 8e-13);
        }
        assert!(matches!(
            lens.trace_seeded_path(
                cx,
                incoming,
                1064.0,
                None,
                Stokes::unpolarized(1.0),
                true,
                Some(pair)
            ),
            Err(CameraError::SequentialLensOutsideModel)
        ));
    });
}

fn camera(profile: usize) -> PhysicalCamera {
    let material = glass(CauchyIor::try_new(1.5, 0.01, 0.0).unwrap(), 0.7);
    let rear = if profile == 0 {
        LensSurface::try_new(0.08, 20.0, 0.012, Some(material))
    } else {
        let mut coefficients = [0.0; 7];
        if profile == 2 {
            coefficients[0] = 2e3;
        }
        LensSurface::try_new_aspheric(0.08, 20.0, -0.5, coefficients, 0.012, Some(material))
    }
    .unwrap();
    let lens = SequentialLens::try_new(
        550.0,
        &[
            rear,
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.024).unwrap(),
        1.0,
        Aperture::try_circular(0.001).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

#[test]
fn camera_proposal_mean_is_the_explicit_direct_plus_ghost_sum() {
    let pair = LensGhostPath::try_new(1, 0).unwrap();
    with_cx(|cx| {
        for profile in 0..3 {
            let direct_camera = camera(profile);
            let lens = direct_camera.sequential_lens().unwrap();
            let sampled_camera = direct_camera
                .clone()
                .with_sequential_lens(lens.clone().with_two_reflection_ghosts(0.5).unwrap())
                .unwrap();
            for (u, v) in [(0.5, 0.5), (0.35, 0.4), (0.65, 0.6)] {
                let pupil = LensSample::try_new(u, v).unwrap();
                let [x, y] = direct_camera.aperture.sample(pupil);
                let rear = &lens.surfaces[0];
                let target = Point3::new(x, y, rear.vertex_z_m + rear.sag(x * x + y * y));
                let local_ray = Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: super::super::super::scaled_unit(
                        target.delta_from(Point3::new(0.0, 0.0, 0.0)),
                    )
                    .unwrap(),
                };
                for wavelength in [450.0, 650.0] {
                    let direct = direct_camera
                        .generate_sequential_spectral_beams_from_tangent_offsets(
                            cx,
                            0.0,
                            0.0,
                            pupil,
                            &[wavelength],
                        )
                        .unwrap()[0]
                        .unwrap();
                    let ghost = lens
                        .trace_ghost_ray(cx, local_ray, wavelength, pair)
                        .unwrap();
                    if u == 0.5 {
                        assert!(ghost.is_some());
                    }
                    let expected = direct.geometry_weight()
                        * (direct.power().i + ghost.map_or(0.0, |path| path.transmittance));
                    let mut sampled_sum = 0.0;
                    for joint_u in [0.5 * u, 0.5 + 0.5 * u] {
                        let sample = LensSample::try_new(joint_u, v).unwrap();
                        let beam = sampled_camera
                            .generate_sequential_spectral_beams_from_tangent_offsets(
                                cx,
                                0.0,
                                0.0,
                                sample,
                                &[wavelength],
                            )
                            .unwrap()[0];
                        sampled_sum +=
                            beam.map_or(0.0, |beam| beam.geometry_weight() * beam.power().i);
                    }
                    close(0.5 * sampled_sum, expected, 3e-12);
                }
            }
            // Both conditional midpoint samples hit the same on-axis pupil.
            // A normal ghost crosses the interior three times even on an asphere.
            let direct = sampled_camera
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::try_new(0.25, 0.5).unwrap(),
                    &[550.0],
                )
                .unwrap()[0]
                .unwrap();
            let ghost = sampled_camera
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::try_new(0.75, 0.5).unwrap(),
                    &[550.0],
                )
                .unwrap()[0]
                .unwrap();
            let group = lens.surfaces[0]
                .medium_after()
                .unwrap()
                .ior()
                .group_index(550.0)
                .unwrap();
            close(
                (ghost.group_delay_s().unwrap() - direct.group_delay_s().unwrap()).to_f64(),
                group * 0.02 / 299_792_458.0,
                3e-23,
            );
        }
    });
}

#[test]
fn proposal_preserves_each_pupil_interval_and_has_bounded_physical_steps() {
    for surfaces in [2, 3, 32] {
        let proposal = GhostSampling::try_new(surfaces, 0.4).unwrap();
        let count = surfaces * (surfaces - 1) / 2;
        let mut pair_index = 0;
        let mut total_probability = 0.4;
        for first in 1..surfaces {
            for second in 0..first {
                let expected = LensGhostPath::try_new(first, second).unwrap();
                for pupil_u in [0.25, 0.5, 0.75] {
                    let input_u = 0.4 + 0.6 * (pair_index as f64 + pupil_u) / count as f64;
                    let (path, pupil, probability) = proposal
                        .sample(LensSample::try_new(input_u, 0.37).unwrap())
                        .unwrap();
                    assert_eq!(path, Some(expected));
                    close(pupil.u(), pupil_u, 2e-12);
                    assert_eq!(pupil.v(), 0.37);
                    close(probability, 0.6 / count as f64, 0.0);
                }
                let steps: Vec<_> = PathSteps::new(surfaces, Some(expected)).collect();
                assert_eq!(steps.len(), surfaces + 2 * (first - second));
                assert!(steps.len() <= 3 * surfaces - 2);
                let reflections: Vec<_> = steps
                    .iter()
                    .filter(|step| step.reflect)
                    .map(|step| (step.surface, step.reverse))
                    .collect();
                assert_eq!(reflections, [(first, false), (second, true)]);
                for step in steps {
                    assert!(
                        (if step.reverse {
                            step.surface + 1
                        } else {
                            step.surface
                        }) < surfaces
                    );
                }
                total_probability += 0.6 / count as f64;
                pair_index += 1;
            }
        }
        assert_eq!(pair_index, count);
        close(total_probability, 1.0, 2e-13);
        for pupil_u in [0.0, 0.25, 0.5, 0.75] {
            let (path, pupil, probability) = proposal
                .sample(LensSample::try_new(0.4 * pupil_u, 0.37).unwrap())
                .unwrap();
            assert_eq!(path, None);
            close(pupil.u(), pupil_u, 2e-16);
            assert_eq!(probability, 0.4);
        }
        let (_, pupil, _) = proposal
            .sample(LensSample::try_new(1.0_f64.next_down(), 0.37).unwrap())
            .unwrap();
        assert!(pupil.u() < 1.0);
    }
}

#[test]
fn exact_zero_reflection_is_a_miss_and_tir_does_not_invent_group_delay() {
    // Native Jones reconstruction retains a weak channel whose mean +/-
    // difference rounded away. Its later underflow must not become a miss.
    let weak = Mueller([
        [0.05, 0.05, 0.0, 0.0],
        [0.05, 0.05, 0.0, 0.0],
        [0.0, 0.0, 1e-12, 0.0],
        [0.0, 0.0, 0.0, 1e-12],
    ]);
    let weak_input = Stokes {
        i: f64::MIN_POSITIVE,
        q: -f64::MIN_POSITIVE,
        u: 0.0,
        v: 0.0,
    };
    assert_eq!(weak.apply_diagonal_jones(weak_input).unwrap().i, 0.0);
    assert!(!reflection_has_zero_support(weak, weak_input));
    let pair = LensGhostPath::try_new(1, 0).unwrap();
    let matched = slab(glass(CauchyIor::try_constant(1.0).unwrap(), 0.0));
    let tir = SequentialLens::try_new_with_sensor_medium(
        550.0,
        Some(glass(CauchyIor::try_constant(1.6).unwrap(), 0.0)),
        &[
            LensSurface::try_new(
                0.02,
                0.0,
                0.2,
                Some(glass(CauchyIor::try_constant(2.0).unwrap(), 0.0)),
            )
            .unwrap(),
            LensSurface::try_new(0.05, 0.0, 0.2, None).unwrap(),
        ],
    )
    .unwrap();
    with_cx(|cx| {
        assert_eq!(
            matched.trace_ghost_ray(cx, ray(0.0), 550.0, pair).unwrap(),
            None
        );
        assert!(matches!(
            tir.trace_seeded_path(
                cx,
                ray(0.8),
                550.0,
                None,
                Stokes::unpolarized(1.0),
                true,
                Some(pair)
            ),
            Err(CameraError::SequentialLensOutsideModel)
        ));
        assert_eq!(
            tir.trace_ghost_ray(cx, ray(0.8), 550.0, pair).unwrap(),
            None
        );
    });
}
