//! Held physical states use the camera's time for the entire optical path.

use super::*;
use fs_exec::RunId;
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraKeyframe, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::tracer::{
    MediumExposure, MediumExposureFilm, RenderExecutionConfig, RenderExecutionError,
    ScatteringFilm, render_cinematic, render_medium_exposure, render_medium_exposure_range,
    render_medium_exposure_range_with_execution, render_medium_exposure_with_execution,
    render_motion, render_scattering_range,
};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::timeline::{
    MediumFrame, MediumTimeline, MediumTimelineBudget, PhysicalMedium, TimelineError,
};

fn shutter(open: f64, close: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        open,
        close - open,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 16 },
        ShotTimeBounds::try_new(-1.0, 2.0).unwrap(),
    )
    .unwrap()
}

fn slab<'field>(
    cx: &Cx<'_>,
    density: &'field [f64],
    absorption: f64,
    temperature: Option<&'field [f64]>,
) -> AbsorbingGrid<'field> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1; 3], density, [-2.0, -2.0, 0.0], [4.0, 4.0, 1.0]),
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        temperature,
        budget(),
    )
    .unwrap()
}

fn timeline<'view, 'field>(
    cx: &Cx<'_>,
    frames: &'view [MediumFrame<'view, 'field>],
) -> MediumTimeline<'view, 'field> {
    MediumTimeline::try_new(cx, frames, MediumTimelineBudget { max_frames: 3 }).unwrap()
}

fn physical_camera(camera: &Camera, shift_x: f64) -> PhysicalCamera {
    PhysicalCamera::try_legacy_compatible(
        Point3::new(camera.eye.x + shift_x, camera.eye.y, camera.eye.z),
        camera.forward,
        camera.up,
        camera.half_tan,
        1.5,
        Aperture::try_circular(0.03).unwrap(),
    )
    .unwrap()
}

fn cinematic(camera: &Camera, shift_x: f64) -> AnimatedCamera {
    AnimatedCamera::try_static(77, 0.0, 1.0, physical_camera(camera, shift_x)).unwrap()
}

fn execution(width: u32, height: u32, workers: usize) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(width, height, workers, 8 << 20, RunId(773)).unwrap()
}

#[test]
fn g1_exposure_integrates_transmittance_and_planck_source_before_time_averaging() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperatures = [[1800.0], [2600.0]];
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 16_384;
        for thermal in [false, true] {
            let first = slab(cx, &density, 0.1, thermal.then_some(&temperatures[0]));
            let second = slab(cx, &density, 2.0, thermal.then_some(&temperatures[1]));
            let frames = [
                MediumFrame {
                    start_s: 0.0,
                    end_s: 0.25,
                    medium: PhysicalMedium::Absorbing(&first),
                },
                MediumFrame {
                    start_s: 0.25,
                    end_s: 1.0,
                    medium: PhysicalMedium::Absorbing(&second),
                },
            ];
            let history = timeline(cx, &frames);
            let exposure = MediumExposure::try_new(cx, &history, shutter(0.0, 1.0)).unwrap();
            let (environment, (emission, scale)) = environment_scene();
            let scene = if thermal { empty_scene() } else { environment };
            let expected = xyz_of_spectrum(|lambda| {
                if thermal {
                    0.25 * planck(1800.0, lambda) * -(-0.1_f64).exp_m1()
                        + 0.75 * planck(2600.0, lambda) * -(-2.0_f64).exp_m1()
                } else {
                    emission.eval(lambda)
                        * scale
                        * (0.25 * (-0.1_f64).exp() + 0.75 * (-2.0_f64).exp())
                }
            });
            let film = render_medium_exposure(&scene, &exposure, cx, &settings).unwrap();
            assert_mean(&film, expected, 0.025);
        }
        // Averaging density/temperature before solving transport is physically
        // different. These gray predictions differ by a factor above 1.4.
        let mean_transmittance = 0.25 * (-0.1_f64).exp() + 0.75 * (-2.0_f64).exp();
        assert!(mean_transmittance > 1.4 * (-(0.25 * 0.1 + 0.75 * 2.0_f64)).exp());
    });
}

#[test]
fn g1_every_camera_scatter_and_shadow_ray_uses_its_sampled_snapshot() {
    with_cx(|cx, _| {
        let density = [1.0];
        let vacuum = slab(cx, &density, 0.0, None);
        let absorbing = slab(cx, &density, 0.7, None);
        let scattering = ScatteringGrid::try_new(
            cx,
            slab(cx, &density, 0.2, None),
            MassScatteringSpectrum::constant(0.5).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 0.25,
                medium: PhysicalMedium::Absorbing(&vacuum),
            },
            MediumFrame {
                start_s: 0.25,
                end_s: 0.75,
                medium: PhysicalMedium::Absorbing(&absorbing),
            },
            MediumFrame {
                start_s: 0.75,
                end_s: 1.0,
                medium: PhysicalMedium::Scattering(&scattering),
            },
        ];
        let history = timeline(cx, &frames);
        let shutter = shutter(0.0, 1.0);
        let exposure = MediumExposure::try_new(cx, &history, shutter).unwrap();
        let (scene, _, _) = lit_receiver();
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.spp = 64;
        settings.max_depth = 4;
        let mut clear = AbsorbingFilm::try_new(&vacuum, &settings).unwrap();
        let mut absorbed = AbsorbingFilm::try_new(&absorbing, &settings).unwrap();
        let mut scattered = ScatteringFilm::try_new(&scattering, &settings).unwrap();
        let mut actual = MediumExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut expected = [0.0; 3];
        let mut counts = [0; 3];
        for sample in 0..settings.spp {
            let before = [
                clear.film().xyz[0],
                absorbed.film().xyz[0],
                scattered.film().xyz[0],
            ];
            render_absorbing_range(
                &scene,
                &vacuum,
                cx,
                &settings,
                &mut clear,
                sample,
                sample + 1,
            )
            .unwrap();
            render_absorbing_range(
                &scene,
                &absorbing,
                cx,
                &settings,
                &mut absorbed,
                sample,
                sample + 1,
            )
            .unwrap();
            render_scattering_range(
                &scene,
                &scattering,
                cx,
                &settings,
                &mut scattered,
                sample,
                sample + 1,
            )
            .unwrap();
            let after = [
                clear.film().xyz[0],
                absorbed.film().xyz[0],
                scattered.film().xyz[0],
            ];
            let time =
                shutter.time_at(shutter.sample_for_stream(settings.seed, 0, u64::from(sample)));
            let state = if time < 0.25 {
                0
            } else if time < 0.75 {
                1
            } else {
                2
            };
            counts[state] += 1;
            render_medium_exposure_range(
                &scene,
                &exposure,
                cx,
                &settings,
                &mut actual,
                sample,
                sample + 1,
            )
            .unwrap();
            // Independent static renders provide the complete selected path,
            // including surface NEE through the same medium as the camera ray.
            // Differencing their accumulated sums entails roundoff, not sampling error.
            for channel in 0..3 {
                expected[channel] += after[state][channel] - before[state][channel];
                assert!(
                    (actual.film().xyz[0][channel] - expected[channel]).abs()
                        <= 2.0e-12 * expected[channel].abs().max(1.0),
                    "sample {sample}, state {state}, channel {channel}"
                );
            }
        }
        assert_eq!(counts, [16, 32, 16]);
        assert!(expected[1] > 0.0);
    });
}

#[test]
fn g3_constant_vacuum_preserves_motion_and_cinematic_pixels() {
    with_cx(|cx, _| {
        let density = [1.0];
        let vacuum = slab(cx, &density, 0.0, None);
        let frames = [MediumFrame {
            start_s: 0.0,
            end_s: 1.0,
            medium: PhysicalMedium::Absorbing(&vacuum),
        }];
        let history = timeline(cx, &frames);
        let shutter = shutter(0.0, 1.0);
        let exposure = MediumExposure::try_new(cx, &history, shutter).unwrap();
        let (mut scene, _, _) = lit_receiver();
        scene.camera.half_tan = 0.3;
        let camera = AnimatedCamera::try_new(vec![
            CameraShot::try_new(
                77,
                0.0,
                1.0,
                vec![
                    CameraKeyframe::try_new(0.0, physical_camera(&scene.camera, -0.75)).unwrap(),
                    CameraKeyframe::try_new(1.0, physical_camera(&scene.camera, 0.75)).unwrap(),
                ],
            )
            .unwrap(),
        ])
        .unwrap();
        let fixed_camera = cinematic(&scene.camera, 0.0);
        let cinematic_exposure = exposure.with_camera(cx, &camera, CutSide::After).unwrap();
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let mut settings = settings(DirectStrategy::Mis);
            settings.spp = 32;
            settings.width = 3;
            settings.height = 2;
            settings.sampler = sampler;
            assert_eq!(
                render_medium_exposure(&scene, &exposure, cx, &settings).unwrap(),
                render_motion(&scene, cx, &settings, shutter).unwrap()
            );
            let moving =
                render_medium_exposure(&scene, &cinematic_exposure, cx, &settings).unwrap();
            assert_eq!(
                moving,
                render_cinematic(&scene, &camera, CutSide::After, cx, &settings, shutter).unwrap()
            );
            assert_ne!(
                moving.xyz,
                render_cinematic(
                    &scene,
                    &fixed_camera,
                    CutSide::After,
                    cx,
                    &settings,
                    shutter
                )
                .unwrap()
                .xyz,
                "the fixture must distinguish sampled camera motion from a frozen camera"
            );
        }
    });
}

#[test]
fn g3_physical_exposure_is_identical_across_tiles_workers_and_append_boundaries() {
    with_cx(|cx, _| {
        let density = [1.0];
        let absorbing = slab(cx, &density, 0.2, None);
        let scattering = ScatteringGrid::try_new(
            cx,
            slab(cx, &density, 0.1, None),
            MassScatteringSpectrum::constant(0.7).unwrap(),
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 0.25,
                medium: PhysicalMedium::Absorbing(&absorbing),
            },
            MediumFrame {
                start_s: 0.25,
                end_s: 1.0,
                medium: PhysicalMedium::Scattering(&scattering),
            },
        ];
        let history = timeline(cx, &frames);
        let exposure = MediumExposure::try_new(cx, &history, shutter(0.0, 1.0)).unwrap();
        let (mut scene, _, _) = lit_receiver();
        scene.camera.half_tan = 0.3;
        let camera = cinematic(&scene.camera, 0.0);
        let exposure = exposure.with_camera(cx, &camera, CutSide::After).unwrap();
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.width = 3;
        settings.height = 2;
        settings.spp = 32;
        settings.max_depth = 5;
        let serial = render_medium_exposure(&scene, &exposure, cx, &settings).unwrap();
        for config in [execution(1, 1, 1), execution(2, 2, 3)] {
            let output =
                render_medium_exposure_with_execution(&scene, &exposure, cx, &settings, &config)
                    .unwrap();
            assert_eq!(output.film, serial);
            assert_eq!(output.report.retained_film_bytes, 0);
            assert_eq!(output.report.memory.used_bytes, 0);
        }
        let mut progressive = MediumExposureFilm::try_new(&exposure, &settings).unwrap();
        render_medium_exposure_range(&scene, &exposure, cx, &settings, &mut progressive, 0, 7)
            .unwrap();
        render_medium_exposure_range_with_execution(
            &scene,
            &exposure,
            cx,
            &settings,
            &mut progressive,
            7,
            19,
            &execution(2, 1, 3),
        )
        .unwrap();
        render_medium_exposure_range(&scene, &exposure, cx, &settings, &mut progressive, 19, 32)
            .unwrap();
        assert_eq!(*progressive.film(), serial);
    });
}

#[test]
fn g4_changed_physical_history_shutter_same_shot_camera_and_cancel_never_append() {
    with_cx(|cx, gate| {
        let density = [1.0];
        let absorbing = slab(cx, &density, 0.2, None);
        let changed = slab(cx, &density, 0.8, None);
        let frames = [MediumFrame {
            start_s: 0.0,
            end_s: 1.0,
            medium: PhysicalMedium::Absorbing(&absorbing),
        }];
        let changed_frames = [MediumFrame {
            medium: PhysicalMedium::Absorbing(&changed),
            ..frames[0]
        }];
        let history = timeline(cx, &frames);
        let changed_history = timeline(cx, &changed_frames);
        let (scene, _, _) = lit_receiver();
        let camera = cinematic(&scene.camera, 0.0);
        let changed_camera = cinematic(&scene.camera, 0.125); // Same shot ID 77.
        let base = MediumExposure::try_new(cx, &history, shutter(0.0, 1.0)).unwrap();
        let exposure = base.with_camera(cx, &camera, CutSide::After).unwrap();
        let wrong_camera = base
            .with_camera(cx, &changed_camera, CutSide::After)
            .unwrap();
        let wrong_history = MediumExposure::try_new(cx, &changed_history, shutter(0.0, 1.0))
            .unwrap()
            .with_camera(cx, &camera, CutSide::After)
            .unwrap();
        let wrong_shutter = MediumExposure::try_new(cx, &history, shutter(0.0, 0.75))
            .unwrap()
            .with_camera(cx, &camera, CutSide::After)
            .unwrap();
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 24;
        let mut film = MediumExposureFilm::try_new(&exposure, &settings).unwrap();
        render_medium_exposure_range(&scene, &exposure, cx, &settings, &mut film, 0, 7).unwrap();
        let before = film.clone();
        for invalid in [wrong_camera, wrong_history, wrong_shutter] {
            assert_eq!(
                render_medium_exposure_range(&scene, &invalid, cx, &settings, &mut film, 7, 24),
                Err(TracerError::InvalidInput)
            );
            assert_eq!(film, before);
        }
        let changed_settings = Settings {
            seed: settings.seed + 1,
            ..settings
        };
        assert_eq!(
            render_medium_exposure_range(
                &scene,
                &exposure,
                cx,
                &changed_settings,
                &mut film,
                7,
                24,
            ),
            Err(TracerError::InvalidInput)
        );
        assert_eq!(film, before);
        let config = execution(1, 1, 2);
        assert!(matches!(
            render_medium_exposure_range_with_execution(
                &scene,
                &wrong_camera,
                cx,
                &settings,
                &mut film,
                7,
                24,
                &config,
            ),
            Err(RenderExecutionError::Tracer(TracerError::InvalidInput))
        ));
        assert_eq!(film, before);
        gate.request();
        assert_eq!(
            render_medium_exposure_range(&scene, &exposure, cx, &settings, &mut film, 7, 24),
            Err(TracerError::Cancelled)
        );
        assert!(matches!(
            render_medium_exposure_range_with_execution(
                &scene, &exposure, cx, &settings, &mut film, 7, 24, &config,
            ),
            Err(RenderExecutionError::Tracer(TracerError::Cancelled))
        ));
        assert_eq!(film, before);
    });
}

#[test]
fn g4_exposure_owns_exact_discontinuity_endpoints_and_requires_full_coverage() {
    with_cx(|cx, _| {
        let density = [1.0];
        let first = slab(cx, &density, 0.1, None);
        let second = slab(cx, &density, 2.0, None);
        let frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 0.5,
                medium: PhysicalMedium::Absorbing(&first),
            },
            MediumFrame {
                start_s: 0.5,
                end_s: 1.0,
                medium: PhysicalMedium::Absorbing(&second),
            },
        ];
        let history = timeline(cx, &frames);
        let (scene, _) = environment_scene();
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 32;
        for (time, state) in [(0.0, &first), (0.5, &second), (1.0, &second)] {
            let exposure = MediumExposure::try_new(cx, &history, shutter(time, time)).unwrap();
            let actual = render_medium_exposure(&scene, &exposure, cx, &settings).unwrap();
            let expected = render_absorbing(&scene, state, cx, &settings).unwrap();
            assert_eq!(actual.xyz, expected.xyz, "instantaneous exposure at {time}");
        }
        for interval in [shutter(-0.1, 0.5), shutter(0.5, 1.1)] {
            assert!(matches!(
                MediumExposure::try_new(cx, &history, interval),
                Err(TracerError::MediumTimeline(
                    TimelineError::ShutterOutsideHistory
                ))
            ));
        }
    });
}
