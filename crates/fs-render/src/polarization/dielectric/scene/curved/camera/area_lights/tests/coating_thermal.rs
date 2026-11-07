//! Directional Kirchhoff power through original area proposals and native MIS.

use super::*;
use crate::conductor::{
    CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus,
    ConductorIorSample, ConductorOptics, ConductorSource,
};
use crate::dielectric::{DielectricEvent, DielectricSurface};
use crate::polarization::dielectric::scene::curved::camera::CoatingThermalEmitter;
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::{
    AreaLightBinding, AreaLightCameraSettings, AreaLightCameraWorkspace, CurvedAreaLightCameraScene,
};
use crate::spectral::SpectralReflectance;
use crate::thermal_surface::CoatingThermalEmission;
use fs_blake3::ContentHash;

const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);

fn metal() -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.3, 3.0).unwrap()
        }),
        ConductorSource::try_new(ContentHash([117; 32]), ConductorDataStatus::Representative)
            .unwrap(),
    )
    .unwrap()
}

fn thermal(
    optics: ConductorOptics,
    side: BoundarySide,
    temperature: f64,
    identity: u8,
) -> CoatingThermalEmitter {
    CoatingThermalEmitter::try_new(
        2,
        side,
        CoatingThermalEmission::try_opaque(
            &[],
            &[],
            optics,
            temperature,
            ContentHash([identity; 32]),
        )
        .unwrap(),
    )
    .unwrap()
}

fn source_binding() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 2,
        interaction: Interaction::SmoothConductor {
            optics: 0,
            medium: 0,
        },
    }
}

fn placed_source(placement: RigidTransform) -> GeometryInstance {
    let original = source(0.0);
    GeometryInstance::try_new(
        2,
        original.geometry_identity(),
        original.geometry().clone(),
        placement,
    )
    .unwrap()
}

#[test]
fn coating_area_emission_uses_converged_angle_and_counts_only_needed_channel_replays() {
    with_cx(|cx| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let bulk = metal();
        let film = CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.6, 0.35).unwrap()
            }),
            90e-9,
            ConductorSource::try_new(ContentHash([122; 32]), ConductorDataStatus::Representative)
                .unwrap(),
        )
        .unwrap();
        let metals = [bulk.with_absorbing_coating(film)];
        let emitter = CoatingThermalEmitter::try_new(
            2,
            BoundarySide::Positive,
            CoatingThermalEmission::try_opaque(
                &[CoatingLayer::Absorbing(&film)],
                &[1400.0],
                bulk,
                2400.0,
                ContentHash([118; 32]),
            )
            .unwrap(),
        )
        .unwrap();
        for cosine in [1.0_f64, 0.6] {
            let angle = cosine.acos();
            let point = Point3::new(2.0 * angle.tan(), 0.0, 2.0);
            let geometry = InstanceScene::try_new(vec![source(0.0), detector(point)]).unwrap();
            let scene = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                vec![source_binding(), detector_binding()],
                100_000,
            )
            .unwrap();
            let light =
                CurvedAreaEmitter::try_new_coating_thermal(cx, &scene, emitter.clone(), 1).unwrap();
            let sample = light.sample_area([4.0 / 9.0, 0.5]).unwrap();
            let seeds = [seed(angle - 0.02), seed(angle + 0.02)];
            let target =
                BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
            let image_options = BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            };
            let mut workspace = BundleImageWorkspace::try_new(1, 2).unwrap();
            let (launch, source, uniform) = light
                .launch(sample, BoundarySide::Positive, 0, Y, 550.0, &seeds)
                .unwrap();
            assert_eq!(uniform, None);
            assert_eq!(launch.normalized_source_stokes, Stokes::unpolarized(1.0));
            let catalog = scene
                .discover_bundle_images_from_surface_with_workspace(
                    cx,
                    launch,
                    source,
                    target,
                    options(),
                    &seeds,
                    image_options,
                    &mut workspace,
                )
                .unwrap();
            assert_eq!(catalog.images().len(), 1);
            let geometric_traces = catalog.trace_evaluations();
            let observed = light
                .incident_irradiance(
                    cx,
                    sample,
                    BoundarySide::Positive,
                    0,
                    Y,
                    550.0,
                    target,
                    options(),
                    &seeds,
                    image_options,
                    &mut workspace,
                )
                .unwrap();
            let radiance = emitter.source().evaluate(550.0, cosine).unwrap();
            let expected = radiance.stokes().i * 0.27 * cosine.powi(4) / 4.0;
            near(observed.spectral_irradiance_w_m2_nm(), expected, 5e-6);
            assert_eq!(observed.image_count(), 1);
            let channels = radiance.polarized_radiance_w_m2_sr_nm();
            let replays = if channels[0] == channels[1] { 0 } else { 2 };
            assert_eq!(observed.trace_evaluations(), geometric_traces + replays);
            if replays == 2 {
                let limited = BundleImageOptions {
                    max_trace_evaluations: geometric_traces + 1,
                    ..image_options
                };
                assert!(matches!(
                    light.incident_irradiance(
                        cx,
                        sample,
                        BoundarySide::Positive,
                        0,
                        Y,
                        550.0,
                        target,
                        options(),
                        &seeds,
                        limited,
                        &mut workspace
                    ),
                    Err(AreaEmitterError::ReplayTraceLimit {
                        required: 2,
                        remaining: 1
                    })
                ));
            }
            assert!(matches!(
                light.launch(sample, BoundarySide::Negative, 0, Y, 550.0, &seeds),
                Err(AreaEmitterError::UnsupportedEmissionSide)
            ));
        }
    });
}

#[test]
fn thermal_s_p_replay_preserves_weak_brewster_power_and_rotated_emitting_normal() {
    with_cx(|cx| {
        // Raw face normal is exactly (2,0,-1): vacuum -> n=2 is Brewster at
        // cos(theta)=1/sqrt(5). Ideal Rs=9/25, Rp=0; reflected light reaches
        // (4/5,0,8/5). The represented cosine retains a tiny positive p residue.
        let point = Point3::new(0.8, 0.0, 1.6);
        let air = glass(1.0, 0.0);
        let dielectric = glass(2.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&dielectric),
        ];
        let metals = [metal()];
        let emitter = thermal(metals[0], BoundarySide::Positive, 2400.0, 119);
        let cosine = 0.6_f64;
        let (s, c) = (0.5 * cosine.acos()).sin_cos();
        let events = [DielectricEvent::Reflection];
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&events),
        }];
        for rotated in [false, true] {
            let quaternion = if rotated {
                [-s, 0.0, 0.0, c]
            } else {
                [0.0, s, 0.0, c]
            };
            let geometry = InstanceScene::try_new(vec![
                mesh(
                    1,
                    vec![[-0.25, -0.25, 0.5], [-0.25, 0.75, 0.5], [0.75, -0.25, 2.5]],
                    vec![[0, 1, 2]],
                    RigidTransform::identity(),
                ),
                placed_source(RigidTransform::try_new(quaternion, [0.0; 3]).unwrap()),
                detector(point),
            ])
            .unwrap();
            let scene = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::Dielectric {
                            outside: 0,
                            inside: 1,
                            surface: DielectricSurface::SMOOTH,
                        },
                    },
                    source_binding(),
                    detector_binding(),
                ],
                100_000,
            )
            .unwrap();
            let light =
                CurvedAreaEmitter::try_new_coating_thermal(cx, &scene, emitter.clone(), 1).unwrap();
            let sample = light.sample_area([4.0 / 9.0, 0.5]).unwrap();
            let target =
                BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
            let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
            let image_options = BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            };
            let (launch, source, uniform) = light
                .launch(sample, BoundarySide::Positive, 0, Y, 550.0, &seeds)
                .unwrap();
            let catalog = scene
                .discover_bundle_images_from_surface_with_workspace(
                    cx,
                    launch,
                    source,
                    target,
                    options(),
                    &seeds,
                    image_options,
                    &mut workspace,
                )
                .unwrap();
            let original = catalog.images()[0].connection();
            let emitted = light
                .image_emission(
                    cx,
                    sample,
                    BoundarySide::Positive,
                    launch,
                    original,
                    target,
                    options(),
                    seeds[0].policy,
                    uniform,
                    2,
                )
                .unwrap();
            let weak = if rotated { 0 } else { 1 };
            let weak_log_power = emitted.channels[weak]
                .unwrap()
                .0
                .path()
                .bundle()
                .interface_log_power();
            assert!(weak_log_power.is_finite());
            assert!(weak_log_power.exp() > 0.0 && weak_log_power.exp() < 1e-28);
            assert_eq!(emitted.traces, 2);
            let radiance = emitter
                .source()
                .evaluate(550.0, cosine)
                .unwrap()
                .polarized_radiance_w_m2_sr_nm();
            let expected =
                radiance[usize::from(rotated)] * 0.27 * cosine / 4.0 * 0.6 * (9.0 / 25.0);
            for chart_axis in [Y, Vec3::new(1.0, 0.0, 0.0)] {
                let observed = light
                    .incident_irradiance(
                        cx,
                        sample,
                        BoundarySide::Positive,
                        0,
                        chart_axis,
                        550.0,
                        target,
                        options(),
                        &seeds,
                        image_options,
                        &mut workspace,
                    )
                    .unwrap();
                near(observed.spectral_irradiance_w_m2_nm(), expected, 5e-6);
                assert_eq!(
                    (observed.image_count(), observed.trace_evaluations()),
                    (1, 3)
                );
            }
        }
    });
}

#[test]
fn coating_nee_and_camera_hits_share_mis_and_the_exact_selected_thermal_snapshot() {
    with_cx(|cx| {
        let (cosine, sine) = (0.6_f64, 0.8_f64);
        let geometry = InstanceScene::try_new(vec![
            placed_source(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-2.0 * sine / cosine, 0.0, -1.0])
                    .unwrap(),
            ),
            detector(Point3::new(0.0, 0.0, 1.0)),
        ])
        .unwrap();
        let host = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal()];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &metals,
            &[],
            &reflectances,
            vec![
                source_binding(),
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
            ],
            100_000,
        )
        .unwrap();
        let emitter = thermal(metals[0], BoundarySide::Positive, 2400.0, 120);
        let back = thermal(metals[0], BoundarySide::Negative, 1400.0, 121);
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
            .unwrap()
            .with_coating_thermal_emitters(cx, &[emitter.clone(), back.clone()])
            .unwrap();
        assert!(CurvedAreaEmitter::from_camera(cx, &camera, 2, 1).is_err());
        let light = CurvedAreaEmitter::from_camera_side(
            cx,
            &camera,
            2,
            BoundarySide::Positive,
            AreaEmitterLimits::default(),
        )
        .unwrap();
        assert_eq!(
            light
                .source_coating_thermal_emitter()
                .unwrap()
                .source()
                .identity(),
            emitter.source().identity()
        );
        let negative = CurvedAreaEmitter::from_camera_side(
            cx,
            &camera,
            2,
            BoundarySide::Negative,
            AreaEmitterLimits::default(),
        )
        .unwrap();
        assert_eq!(
            negative
                .source_coating_thermal_emitter()
                .unwrap()
                .source()
                .identity(),
            back.source().identity()
        );
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(sine, 0.0, cosine),
            policy: BundleBranchPolicy::Sequence(&[]),
        }];
        let binding = AreaLightBinding {
            source: &light,
            receiver_object_id: 9,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Y,
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 2e-7,
        };
        let bindings = [binding];
        let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let settings = AreaLightCameraSettings {
            surfaces: options().scene.trace,
            connection: options(),
            max_area_samples: 1,
            max_trace_evaluations: 128,
        };
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let continuation = |event| {
            assert_eq!(event, 0);
            [sine * sine, 0.75, 0.0]
        };
        let baseline = camera
            .sample_ray_with_samples(cx, ray, 550.0, 0, settings.surfaces, continuation)
            .unwrap();
        let observed = owner
            .sample_ray_with_samples(
                cx,
                ray,
                550.0,
                0,
                settings,
                &mut workspace,
                continuation,
                |event, source| {
                    assert_eq!((event, source), (0, 0));
                    [4.0 / 9.0, 0.5]
                },
            )
            .unwrap();
        let le = emitter.source().evaluate(550.0, cosine).unwrap().stokes().i;
        let distance = 2.0 / cosine;
        let raw_hit = 0.6 * le * (-0.2 * (1.0 + distance)).exp();
        let p_light = distance * distance / (0.27 * cosine);
        let p_camera = cosine / core::f64::consts::PI;
        let denominator = p_light * p_light + p_camera * p_camera;
        let expected_nee = raw_hit * 0.27 * cosine * cosine
            / (distance * distance * core::f64::consts::PI)
            * p_light
            * p_light
            / denominator;
        let expected_hit = raw_hit * p_camera * p_camera / denominator;
        near(baseline.radiance_w_m2_sr_nm(), raw_hit, 5e-6);
        near(observed.area_radiance_w_m2_sr_nm(), expected_nee, 5e-6);
        near(
            observed.radiance_sample().radiance_w_m2_sr_nm(),
            expected_nee + expected_hit,
            5e-6,
        );
        assert_eq!(observed.radiance_sample().path(), baseline.path());
        assert_eq!((observed.area_samples(), observed.images()), (1, 1));
        assert!(observed.trace_evaluations() >= 5);

        let changed = CurvedAreaEmitter::try_new_coating_thermal(
            cx,
            &scene,
            thermal(metals[0], BoundarySide::Positive, 2600.0, 122),
            1,
        )
        .unwrap();
        let stale = [AreaLightBinding {
            source: &changed,
            ..binding
        }];
        assert!(CurvedAreaLightCameraScene::try_new(cx, &camera, &stale).is_err());
        let wrong_side = [AreaLightBinding {
            source_side: BoundarySide::Negative,
            ..binding
        }];
        assert!(CurvedAreaLightCameraScene::try_new(cx, &camera, &wrong_side).is_err());

        #[cfg(feature = "volumes")]
        {
            use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
                CurvedRayleighAreaLightCameraScene, RayleighAreaLightCameraSettings, RayleighAreaLightCameraWorkspace,
            };
            let volumes = [None];
            let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let participating =
                CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &owner).unwrap();
            let mut work = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let replay = participating
                .sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    RayleighAreaLightCameraSettings {
                        area_lights: settings,
                        max_volume_events: 1,
                        max_source_free_flight_draws: 0,
                        max_source_candidates: 0,
                    },
                    &mut work,
                    continuation,
                    |_| panic!("no molecular flight"),
                    |_, _| [4.0 / 9.0, 0.5],
                    |_, _| panic!("no molecular vertex"),
                    |_, _, _| panic!("no scattering source draw"),
                )
                .unwrap();
            assert_eq!(replay.radiance_sample(), observed.radiance_sample());
            assert_eq!(
                replay.area_radiance_w_m2_sr_nm(),
                observed.area_radiance_w_m2_sr_nm()
            );
            assert_eq!(replay.source_free_flight_draws(), 0);
        }
    });
}
