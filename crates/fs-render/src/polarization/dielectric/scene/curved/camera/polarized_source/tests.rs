//! G0/G1/G3: source fractions, physical frames and native reciprocal transport.

use super::super::{
    CameraAnalyzer, CurvedCameraScene, EmissionSides, RadianceEndpoint, RadianceSpectrum,
    SurfaceEmitter,
};
use super::SourcePolarization;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 89,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
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
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn patch(
    id: u64,
    center: Point3,
    normal: Vec3,
    radius: f64,
    placement: RigidTransform,
) -> GeometryInstance {
    let reference = if normal.x.abs() < 0.9 { X } else { Y };
    let u = cross(reference, normal);
    let u = u.scale(1.0 / u.norm());
    let v = cross(normal, u);
    let vertices = [(-1.0, -1.0), (1.0, -1.0), (0.0, 2.0)]
        .map(|(a, b)| {
            [
                center.x + radius * (a * u.x + b * v.x),
                center.y + radius * (a * u.y + b * v.y),
                center.z + radius * (a * u.z + b * v.z),
            ]
        })
        .to_vec();
    GeometryInstance::try_new(
        id,
        hash_domain("native-polarized-emitter", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2]])),
        placement,
    )
    .unwrap()
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}

fn source(state: Stokes, axis: Vec3, radiance: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        9,
        RadianceSpectrum::try_constant(radiance).unwrap(),
        EmissionSides::Both,
    )
    .with_polarization(SourcePolarization::try_constant(axis, state).unwrap())
    .unwrap()
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Z,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g0_source_admission_and_spectral_interpolation_preserve_the_physical_cone() {
    let linear = Stokes {
        i: 1.0,
        q: 1.0,
        u: 0.0,
        v: 0.0,
    };
    let circular = Stokes {
        i: 1.0,
        q: 0.0,
        u: 0.0,
        v: 1.0,
    };
    let law =
        SourcePolarization::try_spectral(X, vec![(400.0, linear), (700.0, circular)]).unwrap();
    assert_eq!(law.reference_axis(), X);
    assert_eq!(law.normalized_stokes(400.0).unwrap(), linear);
    assert_eq!(law.normalized_stokes(700.0).unwrap(), circular);
    assert_eq!(
        law.normalized_stokes(550.0).unwrap(),
        Stokes {
            i: 1.0,
            q: 0.5,
            u: 0.0,
            v: 0.5,
        }
    );
    for wavelength in [400.0, 450.0, 550.0, 625.0, 700.0] {
        assert!(law.normalized_stokes(wavelength).unwrap().degree() <= 1.0);
    }
    for wavelength in [0.0, 399.0, 701.0, f64::NAN, f64::INFINITY] {
        assert!(law.normalized_stokes(wavelength).is_err());
    }
    for state in [
        Stokes::unpolarized(0.0),
        Stokes::unpolarized(2.0),
        Stokes { q: 1.01, ..linear },
        Stokes {
            v: f64::NAN,
            ..linear
        },
    ] {
        assert!(SourcePolarization::try_constant(X, state).is_err());
    }
    for axis in [Vec3::new(0.0, 0.0, 0.0), Vec3::new(f64::INFINITY, 0.0, 0.0)] {
        assert!(SourcePolarization::try_constant(axis, linear).is_err());
    }
    assert!(SourcePolarization::try_spectral(X, vec![(550.0, linear)]).is_err());
    assert!(
        SourcePolarization::try_spectral(X, vec![(550.0, linear), (550.0, circular)],).is_err()
    );
    let calibrated = SurfaceEmitter::new(
        9,
        RadianceSpectrum::try_new_in_band(vec![(400.0, 2.0), (700.0, 3.0)]).unwrap(),
        EmissionSides::Both,
    );
    let bound = calibrated.clone().with_polarization(law.clone()).unwrap();
    assert_eq!(bound.polarization(), Some(&law));
    assert!(
        calibrated
            .with_polarization(
                SourcePolarization::try_spectral(X, vec![(450.0, linear), (650.0, circular)],)
                    .unwrap()
            )
            .is_err()
    );
}

#[test]
fn g1_direct_sources_resolve_full_analyzer_overlap_and_exact_extinction() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(0.0, 0.0, 1.0),
            Z,
            0.5,
            RigidTransform::identity(),
        )])
        .unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let analyzers = [
            CameraAnalyzer::try_new(X, Stokes::unpolarized(1.0)).unwrap(),
            CameraAnalyzer::try_linear(X, 0.0).unwrap(),
            CameraAnalyzer::try_new(
                X,
                Stokes {
                    i: 0.5,
                    q: -0.5,
                    u: 0.0,
                    v: 0.0,
                },
            )
            .unwrap(),
            CameraAnalyzer::try_linear(X, core::f64::consts::FRAC_PI_4).unwrap(),
            CameraAnalyzer::try_circular(X, true).unwrap(),
            CameraAnalyzer::try_circular(X, false).unwrap(),
            CameraAnalyzer::try_new(
                X,
                Stokes {
                    i: 0.6,
                    q: 0.2,
                    u: 0.1,
                    v: 0.1,
                },
            )
            .unwrap(),
        ];
        for state in [
            Stokes::unpolarized(1.0),
            Stokes {
                i: 1.0,
                q: 1.0,
                u: 0.0,
                v: 0.0,
            },
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 1.0,
                v: 0.0,
            },
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: 1.0,
            },
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: -1.0,
            },
            Stokes {
                i: 1.0,
                q: 0.2,
                u: -0.4,
                v: 0.3,
            },
        ] {
            let emitters = [source(state, X, 5.0)];
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let scalar = camera
                .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| {
                    panic!("opaque source has no surface draw")
                })
                .unwrap();
            near(scalar.radiance_w_m2_sr_nm(), 5.0, 2e-14);
            for analyzer in analyzers {
                let sample = camera
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("opaque source has no surface draw"),
                    )
                    .unwrap();
                let row = analyzer.row();
                let expected = 5.0 * (row.i + row.q * state.q + row.u * state.u + row.v * state.v);
                near(sample.radiance_w_m2_sr_nm(), expected, 3e-14);
                if expected == 0.0 {
                    assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
                }
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
                assert_eq!(sample.path().surface_events(), 1);
            }
        }
        let legacy = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Both,
        )];
        let explicit = [source(Stokes::unpolarized(1.0), X, 5.0)];
        let before = CurvedCameraScene::try_new(cx, &native, &legacy, None).unwrap();
        let after = CurvedCameraScene::try_new(cx, &native, &explicit, None).unwrap();
        assert_eq!(
            before
                .sample_ray(cx, ray(), 550.0, 0, settings(), 37, 1)
                .unwrap(),
            after
                .sample_ray(cx, ray(), 550.0, 0, settings(), 37, 1)
                .unwrap(),
        );
    });
}

#[test]
fn g3_source_axes_follow_body_placement_and_parallel_polarized_axes_refuse() {
    with_cx(|cx| {
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let state = Stokes {
            i: 1.0,
            q: 1.0,
            u: 0.0,
            v: 0.0,
        };
        let analyzer = CameraAnalyzer::try_linear(X, 0.0).unwrap();
        for angle in [
            0.0,
            core::f64::consts::FRAC_PI_4,
            core::f64::consts::FRAC_PI_2,
        ] {
            let placement = RigidTransform::try_new(
                [0.0, 0.0, (0.5 * angle).sin(), (0.5 * angle).cos()],
                [0.0, 0.0, 1.0],
            )
            .unwrap();
            let geometry = InstanceScene::try_new(vec![patch(
                9,
                Point3::new(0.0, 0.0, 0.0),
                Z,
                0.5,
                placement,
            )])
            .unwrap();
            let native =
                CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000)
                    .unwrap();
            let emitters = [source(state, X, 4.0)];
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let sample = camera
                .sample_ray_analyzed_with_samples(cx, ray(), 550.0, 0, analyzer, settings(), |_| {
                    panic!("opaque source")
                })
                .unwrap();
            near(
                sample.radiance_w_m2_sr_nm(),
                4.0 * angle.cos().powi(2),
                2e-13,
            );

            let singular = [source(state, Z, 4.0)];
            let camera = CurvedCameraScene::try_new(cx, &native, &singular, None).unwrap();
            assert!(
                camera
                    .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| panic!(
                        "opaque source"
                    ),)
                    .is_err()
            );
        }
    });
}

#[test]
fn g1_native_camera_queries_source_polarization_at_the_actual_wavelength() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(0.0, 0.0, 1.0),
            Z,
            0.5,
            RigidTransform::identity(),
        )])
        .unwrap();
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let law = SourcePolarization::try_spectral(
            X,
            vec![
                (
                    400.0,
                    Stokes {
                        i: 1.0,
                        q: 1.0,
                        u: 0.0,
                        v: 0.0,
                    },
                ),
                (
                    700.0,
                    Stokes {
                        i: 1.0,
                        q: 0.0,
                        u: 0.0,
                        v: 1.0,
                    },
                ),
            ],
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new_in_band(vec![(400.0, 8.0), (700.0, 8.0)]).unwrap(),
            EmissionSides::Both,
        )
        .with_polarization(law)
        .unwrap()];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        for (wavelength, q, v) in [(400.0, 1.0, 0.0), (550.0, 0.5, 0.5), (700.0, 0.0, 1.0)] {
            for (analyzer, expected) in [
                (CameraAnalyzer::try_linear(X, 0.0).unwrap(), 4.0 * (1.0 + q)),
                (
                    CameraAnalyzer::try_circular(X, true).unwrap(),
                    4.0 * (1.0 + v),
                ),
            ] {
                let observed = camera
                    .sample_ray_analyzed_with_samples(
                        cx,
                        ray(),
                        wavelength,
                        0,
                        analyzer,
                        settings(),
                        |_| panic!("opaque source"),
                    )
                    .unwrap();
                near(observed.radiance_w_m2_sr_nm(), expected, 3e-14);
            }
        }
    });
}

#[test]
fn g1_native_tir_preserves_emitted_linear_and_circular_phase() {
    with_cx(|cx| {
        // Forward propagation is +z -> +x at 45-degree incidence, n=3/2 -> 1.
        // In the common +y s-axis rs=4/5-3i/5 and rp=7/25-24i/25, hence
        // rs*conj(rp)=4/5+3i/5. Thus U'=4U/5-3V/5 and V'=3U/5+4V/5.
        // TIR has unit branch probability: there is no unpolarized Fresnel
        // factor or extra one-half source/analyzer normalization to divide out.
        let normal = Vec3::new(1.0, 0.0, -1.0).scale(core::f64::consts::FRAC_1_SQRT_2);
        let geometry = InstanceScene::try_new(vec![
            patch(
                1,
                Point3::new(0.0, 0.0, 0.0),
                normal,
                0.08,
                RigidTransform::identity(),
            ),
            patch(
                9,
                Point3::new(0.0, 0.0, -0.4),
                Z,
                0.08,
                RigidTransform::identity(),
            ),
        ])
        .unwrap();
        let bulk = glass(1.5, 0.0);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&bulk),
            SceneMedium::Homogeneous(&air),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let camera_ray = Ray {
            origin: Point3::new(0.4, 0.0, 0.0),
            dir: X.scale(-1.0),
        };
        for state in [
            Stokes {
                i: 1.0,
                q: 0.2,
                u: 0.4,
                v: 0.6,
            },
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: 1.0,
            },
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 1.0,
                v: 0.0,
            },
        ] {
            let emitters = [source(state, Y, 5.0)];
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let out_u = 0.8 * state.u - 0.6 * state.v;
            let out_v = 0.6 * state.u + 0.8 * state.v;
            for (analyzer, expected) in [
                (
                    CameraAnalyzer::try_linear(Y, core::f64::consts::FRAC_PI_4).unwrap(),
                    2.5 * (1.0 + out_u),
                ),
                (
                    CameraAnalyzer::try_circular(Y, true).unwrap(),
                    2.5 * (1.0 + out_v),
                ),
                (
                    CameraAnalyzer::try_circular(Y, false).unwrap(),
                    2.5 * (1.0 - out_v),
                ),
            ] {
                let sample = camera
                    .sample_ray_analyzed_with_samples(
                        cx,
                        camera_ray,
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| [0.5, 0.5, 0.0],
                    )
                    .unwrap();
                near(sample.radiance_w_m2_sr_nm(), expected, 2e-10);
                near(sample.path().traveled_m(), 0.8, 2e-7);
                assert_eq!(sample.path().surface_events(), 2);
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
            }
        }
    });
}

#[test]
fn g1_polarized_source_histories_keep_native_group_delay_and_shutter_mean() {
    use super::super::transient::exposure::ExposureWindow;
    use super::super::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
    use super::super::transient::{SurfaceEmissionHistoryBinding, TransientCurvedCameraScene};
    use fs_math::dd::Dd;

    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(0.0, 0.0, 1.0),
            Z,
            0.5,
            RigidTransform::identity(),
        )])
        .unwrap();
        let material = glass(1.2, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let emitters = [source(
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: 1.0,
            },
            X,
            5.0,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let knots = [(0.0, 0.0), (2e-8, 2.0)];
        let history = PiecewiseLinearHistory::try_new(cx, Dd::from_f64(0.0), &knots).unwrap();
        let histories = [SurfaceEmissionHistoryBinding {
            object_id: 9,
            source_identity: hash_domain("polarized-source-history", b"linear-ramp"),
            history: &history,
        }];
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &histories)
            .unwrap()
            .with_analyzer(CameraAnalyzer::try_circular(X, true).unwrap());
        let instant = timed
            .sample_ray_with_samples_at_time(
                cx,
                ray(),
                550.0,
                0,
                Dd::from_f64(1e-8),
                settings(),
                |_| panic!("opaque source"),
            )
            .unwrap();
        let window = ExposureWindow::try_new(Dd::from_f64(1e-8), Dd::from_f64(1.4e-8)).unwrap();
        let exposure = timed
            .sample_ray_with_samples_during_exposure(
                cx,
                ray(),
                550.0,
                0,
                window,
                settings(),
                |_| panic!("opaque source"),
            )
            .unwrap();
        let delay = 1.2 / 299_792_458.0;
        near(
            instant.radiance_sample().radiance_w_m2_sr_nm(),
            5.0 * (1e-8 - delay) * 1e8,
            2e-7,
        );
        let expected_mean = 5.0 * (1.2e-8 - delay) * 1e8;
        near(
            exposure.mean_radiance_sample().radiance_w_m2_sr_nm(),
            expected_mean,
            2e-7,
        );
        near(
            exposure.integrated_radiance_j_m2_sr_nm(),
            expected_mean * 4e-9,
            1e-15,
        );
        assert_eq!(
            instant.radiance_sample().path(),
            exposure.mean_radiance_sample().path(),
        );
    });
}

#[cfg(feature = "volumes")]
#[test]
fn g1_molecular_camera_contracts_actual_source_polarization_after_scattering() {
    use super::super::participating::CurvedRayleighCameraScene;
    use crate::polarization::dielectric::scene::curved::participating::{
        CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
    };
    use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};

    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![patch(
            9,
            Point3::new(1.0, 0.0, 0.5),
            X,
            0.5,
            RigidTransform::identity(),
        )])
        .unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let delta = 0.03;
        let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        for q in [-0.6, 0.0, 0.6] {
            let emitters = [source(
                Stokes {
                    i: 1.0,
                    q,
                    u: 0.0,
                    v: 0.4,
                },
                Y,
                4.0,
            )];
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let sample = molecular
                .sample_ray_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    CurvedRayleighTraceSettings {
                        surfaces: settings(),
                        max_volume_events: 8,
                    },
                    |_| panic!("opaque source"),
                    |flight| {
                        if flight == 0 {
                            let isotropic = 3.0 * delta / (2.0 + delta);
                            [
                                1.0 - (-0.35_f64).exp(),
                                isotropic + 0.5 * (1.0 - isotropic) + f64::EPSILON,
                                0.0,
                            ]
                        } else {
                            [1.0 - 1e-12, 0.5, 0.0]
                        }
                    },
                )
                .unwrap();
            // Right-angle phase M01/M00=(1-delta)/(1+delta); the native
            // scalar proposal already divides by M00. Only the 1.5 m true
            // absorption path remains, not a second scattering survival.
            let expected = 4.0 * (-0.3_f64).exp() * (1.0 + (1.0 - delta) / (1.0 + delta) * q);
            near(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                expected,
                3e-6,
            );
            near(sample.radiance_sample().path().traveled_m(), 1.5, 2e-7);
            assert_eq!(sample.volume_events(), 1);
            assert_eq!(sample.free_flight_draws(), 2);
        }
    });
}
