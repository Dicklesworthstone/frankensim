use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
    fresnel_dielectric,
};
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, GradientIndexOptions, GradientIndexSample,
};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::SceneMedium;
use crate::polarization::dielectric::scene::{SceneError, SurfaceBinding};
use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 13,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, z: f64, center_x: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("curved-camera-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center_x - radius, -radius, z],
                [center_x + radius, -radius, z],
                [center_x, radius, z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 4.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}
fn ray(z: f64, direction: f64) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, z),
        dir: Vec3::new(0.0, 0.0, direction),
    }
}
fn emitter(id: u64, radiance: f64, sides: EmissionSides) -> SurfaceEmitter {
    SurfaceEmitter::new(id, RadianceSpectrum::try_constant(radiance).unwrap(), sides)
}
fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn g0_spectral_radiance_keeps_original_values_and_rejects_incomplete_or_invalid_sources() {
    let spectrum =
        RadianceSpectrum::try_new(vec![(380.0, 0.0), (500.0, 4.0), (780.0, 2.0)]).unwrap();
    assert_eq!(spectrum.evaluate(380.0).unwrap(), 0.0);
    assert_eq!(spectrum.evaluate(440.0).unwrap(), 2.0);
    assert_eq!(spectrum.evaluate(640.0).unwrap(), 3.0);
    assert_eq!(spectrum.evaluate(780.0).unwrap(), 2.0);
    for knots in [
        vec![],
        vec![(380.0, 1.0)],
        vec![(400.0, 1.0), (780.0, 1.0)],
        vec![(380.0, 1.0), (500.0, 1.0), (500.0, 2.0), (780.0, 1.0)],
        vec![(380.0, -1.0), (780.0, 1.0)],
        vec![(380.0, f64::NAN), (780.0, 1.0)],
        vec![(380.0, 1.0), (780.0, f64::INFINITY)],
    ] {
        assert!(matches!(
            RadianceSpectrum::try_new(knots),
            Err(CameraRadianceError::InvalidSpectrum)
        ));
    }
    for wavelength in [379.0, 781.0, f64::NAN] {
        assert!(spectrum.evaluate(wavelength).is_err());
    }
}

#[test]
fn g1_native_homogeneous_rays_preserve_beer_radiance_and_per_side_emission() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 1.0)]).unwrap();
        let material = glass(1.3, 0.4);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        for sides in [EmissionSides::Front, EmissionSides::Both] {
            let emitters = [emitter(9, 3.0, sides)];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            for (z, direction, visible) in
                [(0.0, 1.0, sides == EmissionSides::Both), (2.0, -1.0, true)]
            {
                let sample = camera
                    .sample_ray(cx, ray(z, direction), 550.0, 0, settings(), 31, 1)
                    .unwrap();
                near(
                    sample.radiance_w_m2_sr_nm(),
                    if visible { 3.0 * (-0.4_f64).exp() } else { 0.0 },
                    2e-6,
                );
                near(sample.transport_weight(), (-0.4_f64).exp(), 1e-6);
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
                assert_eq!(sample.path().termination(), PathTermination::Opaque);
            }
        }
    });
}

#[test]
fn g1_branch_average_uses_one_fresnel_probability_one_eta_ratio_and_physical_beer_distance() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, 1.0, 0.0, 1.0)])
                .unwrap();
        let incident = glass(1.0, 0.1);
        let transmitted = glass(1.5, 0.4);
        let media = [
            SceneMedium::Homogeneous(&incident),
            SceneMedium::Homogeneous(&transmitted),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9, 3.5, EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let reflected = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [0.5, 0.5, 0.0])
            .unwrap();
        let transmitted = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                [0.5, 0.5, 0.99]
            })
            .unwrap();
        assert_eq!(reflected.radiance_w_m2_sr_nm(), 0.0);
        let beer = (-0.1_f64 * 0.4 - 0.4 * 0.6).exp();
        let conditional = 3.5 * beer / (1.5 * 1.5);
        near(transmitted.radiance_w_m2_sr_nm(), conditional, 3e-6);
        near(transmitted.path().beam().power().i, beer, 2e-6);
        near(transmitted.transport_weight(), beer / (1.5 * 1.5), 2e-6);
        // The sampler already divided by its branch proposal. Multiplying
        // source radiance by another Fresnel transmission would bias this mean.
        let r = 0.04;
        near(
            r * reflected.radiance_w_m2_sr_nm() + (1.0 - r) * transmitted.radiance_w_m2_sr_nm(),
            (1.0 - r) * conditional,
            3e-6,
        );
        let mut limited = settings();
        limited.max_surface_events = 1;
        assert!(matches!(
            camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, limited, |_| [
                0.5, 0.5, 0.99
            ]),
            Err(CameraRadianceError::Unresolved(
                PathTermination::SurfaceLimit
            ))
        ));
    });
}

struct AffineField;
impl GradientIndexField for AffineField {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        GradientIndexSample::try_new(
            1.2 + 0.2 * p.x,
            Vec3::new(0.2, 0.0, 0.0),
            0.3,
            Some(1.2 + 0.2 * p.x),
        )
    }
}

#[test]
fn g1_grin_camera_sees_a_small_emitter_missed_by_a_straight_ray_and_uses_the_endpoint_index() {
    with_cx(|cx, _| {
        // Analytic n(x)=n0+g*x, launch along z: n(z)=n0*cosh(g*z/n0),
        // s(z)=n0/g*sinh(g*z/n0), x(z)=(n(z)-n0)/g.
        let n0 = 1.2_f64;
        let g = 0.2_f64;
        let z = 1.0_f64;
        let n_end = n0 * (g * z / n0).cosh();
        let distance = n0 / g * (g * z / n0).sinh();
        let x = (n_end - n0) / g;
        let geometry = InstanceScene::try_new(vec![plane(9, z, x, 0.01)]).unwrap();
        let field = AffineField;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([19; 32]),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = [emitter(9, 5.0, EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let sample = camera
            .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 1)
            .unwrap();
        let expected = 5.0 * (n0 / n_end).powi(2) * (-0.3 * distance).exp();
        near(sample.radiance_w_m2_sr_nm(), expected, 3e-6);
        near(sample.path().beam().ray().origin.x, x, 2e-6);
        near(sample.path().traveled_m(), distance, 2e-6);
        assert_eq!(
            sample.path().medium_snapshot_identity(),
            Some(ContentHash([19; 32]))
        );
        assert!(x > 0.01, "straight camera ray misses the finite emitter");
        let mut limited = settings();
        limited.max_distance_m = 0.5;
        assert!(matches!(
            camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, limited, 31, 1),
            Err(CameraRadianceError::Unresolved(
                PathTermination::DistanceLimit
            ))
        ));
    });
}

#[test]
fn g1_grin_slab_environment_uses_refracted_direction_and_one_total_phase_index_ratio() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 10.0)]).unwrap();
        let field = AffineField;
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::GradientIndex {
                field: &field,
                identity: ContentHash([23; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 1,
                interaction: Interaction::Dielectric {
                    inside: 0,
                    outside: 1,
                    surface: DielectricSurface::SMOOTH,
                },
            }],
            100_000,
        )
        .unwrap();
        let pixels = (0..32).map(|i| [0.5 + 0.5 * (i % 8) as f32; 3]).collect();
        let environment = EnvironmentMap::try_from_linear_srgb(8, 4, pixels, 0.0).unwrap();
        let calibrated = EnvironmentRadiance::try_new(&environment, 2.5).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], Some(calibrated)).unwrap();
        let sample = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                [0.5, 0.5, 0.99]
            })
            .unwrap();
        let distance = 6.0 * (1.0_f64 / 6.0).sinh();
        let x = 0.2 * distance;
        let direction = Vec3::new(x, 0.0, (1.0 - x * x).sqrt());
        let source = environment.evaluate(direction).unwrap().emission;
        let straight = environment
            .evaluate(Vec3::new(0.0, 0.0, 1.0))
            .unwrap()
            .emission;
        let emitted = source.0.eval(550.0) * source.1 * calibrated.scale_w_m2_sr_nm();
        assert_ne!(source.1, straight.1);
        near(
            sample.radiance_w_m2_sr_nm(),
            emitted * 1.2 * 1.2 * (-0.3 * distance).exp(),
            1e-5,
        );
        near(sample.path().beam().ray().dir.x, direction.x, 2e-6);
        assert_eq!(sample.endpoint(), RadianceEndpoint::Environment);
    });
}

#[test]
fn g0_opaque_bindings_unsupported_scattering_and_infinite_absorption_are_explicit() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 10.0, 1.0)]).unwrap();
        let environment = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let calibrated = EnvironmentRadiance::try_new(&environment, 1.0).unwrap();
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                EnvironmentRadiance::try_new(&environment, invalid),
                Err(CameraRadianceError::InvalidEnvironmentScale)
            ));
        }
        for absorption in [0.0, 0.2] {
            let medium = glass(1.0, absorption);
            let media = [SceneMedium::Homogeneous(&medium)];
            let scene =
                CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000)
                    .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], Some(calibrated)).unwrap();
            let sample = camera
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 1, 1)
                .unwrap();
            if absorption > 0.0 {
                assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
                assert_eq!(sample.transport_weight(), 0.0);
            } else {
                assert!(sample.radiance_w_m2_sr_nm() > 0.0);
                near(sample.transport_weight(), 1.0, 1e-14);
            }
            let wrong = [emitter(77, 1.0, EmissionSides::Both)];
            assert!(matches!(
                CurvedCameraScene::try_new(cx, &scene, &wrong, None),
                Err(CameraRadianceError::InvalidEmitter { object_id: 77 })
            ));
            let duplicate = [
                emitter(9, 1.0, EmissionSides::Both),
                emitter(9, 2.0, EmissionSides::Both),
            ];
            assert!(CurvedCameraScene::try_new(cx, &scene, &duplicate, None).is_err());
        }
        let medium = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&medium)];
        let rough = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 0,
                    surface: DielectricSurface::POLISHED,
                },
            }],
            100_000,
        )
        .unwrap();
        assert!(matches!(
            CurvedCameraScene::try_new(cx, &rough, &[], None),
            Err(CameraRadianceError::UnsupportedSurface { object_id: 9 })
        ));
        let smooth = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 0,
                    surface: DielectricSurface::SMOOTH,
                },
            }],
            100_000,
        )
        .unwrap();
        let source = [emitter(9, 1.0, EmissionSides::Both)];
        assert!(matches!(
            CurvedCameraScene::try_new(cx, &smooth, &source, None),
            Err(CameraRadianceError::InvalidEmitter { object_id: 9 })
        ));
        let camera = CurvedCameraScene::try_new(cx, &smooth, &[], None).unwrap();
        gate.request();
        assert!(matches!(
            camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 1, 1),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}

#[test]
fn g0_finite_source_underflow_is_unresolved_but_infinite_absorption_is_exact_zero() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 1.0)]).unwrap();
        let material = glass(1.0, 800.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        // exp(-800) underflows before this finite source can rescale it into
        // a representable radiance. The consumer must not call that black.
        let emitters = [emitter(9, 1e308, EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        assert!(matches!(
            camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 1, 1),
            Err(CameraRadianceError::Unresolved(
                PathTermination::Extinguished
            ))
        ));

        let missed_geometry = InstanceScene::try_new(vec![plane(9, 1.0, 10.0, 1.0)]).unwrap();
        let escaped =
            CurvedPolarizedScene::try_new(cx, &missed_geometry, &media, vec![opaque(9)], 100_000)
                .unwrap();
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1e30; 3]], 0.0).unwrap();
        let calibrated = EnvironmentRadiance::try_new(&map, f64::MAX).unwrap();
        // The unattenuated calibrated source would overflow, but propagation
        // through an infinite absorbing medium is exactly zero before lookup.
        let camera = CurvedCameraScene::try_new(cx, &escaped, &[], Some(calibrated)).unwrap();
        let sample = camera
            .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 1, 1)
            .unwrap();
        assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(sample.transport_weight(), 0.0);
    });
}

#[test]
fn g0_cancelled_camera_never_queries_the_source_index_field() {
    struct UnqueriedField;
    impl GradientIndexField for UnqueriedField {
        fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            panic!("cancelled camera queried its field");
        }
    }
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 1.0)]).unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &UnqueriedField,
            identity: ContentHash([31; 32]),
        }];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        gate.request();
        for result in [
            camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 1, 1),
            camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [0.5; 3]),
        ] {
            assert!(matches!(
                result,
                Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                    SceneError::Cancelled
                )))
            ));
        }
    });
}

#[test]
fn g3_reciprocal_non_coplanar_interfaces_agree_for_unpolarized_endpoints() {
    let first = glass(1.0, 0.0);
    let second = glass(1.5, 0.0);
    let last = glass(1.2, 0.0);
    let normal1 = Vec3::new(0.4, 0.1, -(1.0_f64 - 0.4 * 0.4 - 0.1 * 0.1).sqrt());
    let normal2 = Vec3::new(0.1, 0.35, -(1.0_f64 - 0.1 * 0.1 - 0.35 * 0.35).sqrt());
    let frame =
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap();
    let beam = PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        frame,
        550.0,
        first,
        Stokes::unpolarized(1.0),
    )
    .unwrap();
    let a = beam.sample_smooth_interface(normal1, second, 0.99).unwrap();
    let b = a.beam.sample_smooth_interface(normal2, last, 0.99).unwrap();
    let reverse_frame =
        PolarizationFrame::try_new(b.beam.ray().dir.scale(-1.0), b.beam.frame().horizontal())
            .unwrap();
    let reverse = PolarizedBeam::try_new(
        beam.ray().origin,
        reverse_frame,
        550.0,
        last,
        Stokes::unpolarized(1.0),
    )
    .unwrap();
    let c = reverse
        .sample_smooth_interface(normal2, second, 0.99)
        .unwrap();
    let d = c
        .beam
        .sample_smooth_interface(normal1, first, 0.99)
        .unwrap();
    near(
        b.beam.power().i * a.pdf * b.pdf,
        d.beam.power().i * c.pdf * d.pdf,
        3e-14,
    );
    near(d.beam.ray().dir.x, 0.0, 3e-14);
    near(d.beam.ray().dir.y, 0.0, 3e-14);
    near(d.beam.ray().dir.z, -1.0, 3e-14);
    // This is reciprocal POWER; radiance itself differs by the endpoint n²
    // ratio. The camera consumer supplies that ratio once after its native walk.
    let reflectance = fresnel_dielectric(1.0, 1.0, 1.5).unwrap().reflectance;
    near(reflectance, 0.04, 1e-14);
}
