//! Accepted scene changes share one clock with lighting, camera and physical media.

use super::*;
use fs_exec::RunId;
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraKeyframe, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::tracer::{
    Film, RenderExecutionConfig, RenderExecutionError, SceneExposure, SceneExposureBudget,
    SceneExposureFilm, SceneFrame, render_cinematic, render_motion, render_motion_range,
    render_scene_exposure, render_scene_exposure_range, render_scene_exposure_range_with_execution,
    render_scene_exposure_with_execution,
};

fn shutter(open: f64, close: f64) -> ShutterInterval {
    ShutterInterval::resolve(
        open,
        close - open,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::StratifiedCounterV1 { strata: 16 },
        ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
    )
    .unwrap()
}

fn frames<'a>(first: &'a Scene, second: &'a Scene, split: f64) -> [SceneFrame<'a>; 2] {
    [
        SceneFrame {
            start_s: 0.0,
            end_s: split,
            scene: first,
        },
        SceneFrame {
            start_s: split,
            end_s: 1.0,
            scene: second,
        },
    ]
}

fn exposure<'a>(
    cx: &Cx<'_>,
    frames: &'a [SceneFrame<'a>],
    shutter: ShutterInterval,
) -> SceneExposure<'a> {
    SceneExposure::try_new(
        cx,
        frames,
        shutter,
        SceneExposureBudget {
            max_frames: 4,
            max_total_primitives: 16,
            max_total_light_triangles: 32,
        },
    )
    .unwrap()
}

fn quad(x: f64, z: f64, half_width: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [x - half_width, -half_width, z],
                [x + half_width, -half_width, z],
                [x + half_width, half_width, z],
                [x - half_width, half_width, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material,
        emission: None,
    }
}

fn illuminated(temperature: f64, emissivity: f64) -> Scene {
    let mut s = scene(Material::Lambertian {
        reflectance: lift_rgb([0.4; 3]),
    });
    s.camera.eye.z = 0.5;
    s.primitives.push(quad(
        0.0,
        2.0,
        1.0,
        Material::GrayDiffuseThermal {
            surface: surface(temperature, emissivity),
        },
    ));
    s
}

fn execution(width: u32, height: u32, workers: usize) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(width, height, workers, 8 << 20, RunId(540)).unwrap()
}

fn assert_sum(actual: [f64; 3], expected: [f64; 3], sample: u32) {
    for channel in 0..3 {
        assert!(
            (actual[channel] - expected[channel]).abs()
                <= 2.0e-12 * expected[channel].abs().max(1.0),
            "sample {sample}, channel {channel}: {} vs {}",
            actual[channel],
            expected[channel]
        );
    }
}

#[test]
fn g1_scene_exposure_averages_physical_emission_and_reflection_not_material_parameters() {
    with_cx(|cx| {
        let make = |temperature, emissivity| {
            let mut s = scene(Material::GrayDiffuseThermal {
                surface: surface(temperature, emissivity),
            });
            s.environment =
                Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[20.0; 3]], 0.0).unwrap());
            s
        };
        let cold = make(900.0, 0.25);
        let hot = make(2500.0, 0.8);
        let states = frames(&cold, &hot, 0.25);
        let exposure = exposure(cx, &states, shutter(0.0, 1.0));
        let exterior = cold
            .environment
            .as_ref()
            .unwrap()
            .evaluate(Vec3::new(0.0, 0.0, 1.0))
            .unwrap()
            .emission;
        // Independent SI Planck formula, including metres-to-nanometres conversion.
        let planck = |temperature: f64, wavelength_nm: f64| {
            let wavelength_m = wavelength_nm * 1.0e-9;
            let h = 6.626_070_15e-34;
            let c = 299_792_458.0;
            let k = 1.380_649e-23;
            2.0 * h * c * c
                / wavelength_m.powi(5)
                / (h * c / (wavelength_m * k * temperature)).exp_m1()
                * 1.0e-9
        };
        let expected = xyz_of_spectrum(|lambda| {
            let reflected = exterior.0.eval(lambda) * exterior.1;
            0.25 * (0.25 * planck(900.0, lambda) + 0.75 * reflected)
                + 0.75 * (0.8 * planck(2500.0, lambda) + 0.2 * reflected)
        });
        let settings = Settings {
            spp: 8192,
            ..settings(DirectStrategy::BsdfOnly)
        };
        let image = render_scene_exposure(&exposure, cx, &settings).unwrap();
        for (actual, expected) in image.xyz[0].iter().zip(expected) {
            let mean = actual / f64::from(image.spp_done);
            assert!(
                (mean / expected - 1.0).abs() < 0.025,
                "{mean} vs {expected}"
            );
        }
    });
}

#[test]
fn g1_every_thermal_light_and_secondary_ray_uses_the_selected_scene() {
    with_cx(|cx| {
        let first = illuminated(2200.0, 0.5);
        let mut second = illuminated(2800.0, 0.8);
        second.primitives.push(quad(
            3.0,
            2.5,
            0.6,
            Material::GrayDiffuseThermal {
                surface: surface(3100.0, 0.4),
            },
        ));
        let states = frames(&first, &second, 0.25);
        let shutter = shutter(0.0, 1.0);
        let exposure = exposure(cx, &states, shutter);
        let settings = Settings {
            spp: 64,
            max_depth: 3,
            sampler: Sampler::OwenSobolFullPath,
            ..settings(DirectStrategy::PowerMis)
        };
        let mut frozen = [Film::new(1, 1), Film::new(1, 1)];
        let mut actual = SceneExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut expected = [0.0; 3];
        let mut counts = [0; 2];
        for sample in 0..settings.spp {
            let before = [frozen[0].xyz[0], frozen[1].xyz[0]];
            for (state, s) in [&first, &second].into_iter().enumerate() {
                render_motion_range(
                    s,
                    cx,
                    &settings,
                    &mut frozen[state],
                    sample,
                    sample + 1,
                    shutter,
                )
                .unwrap();
            }
            let time =
                shutter.time_at(shutter.sample_for_stream(settings.seed, 0, u64::from(sample)));
            let state = usize::from(time >= 0.25);
            counts[state] += 1;
            for (channel, value) in expected.iter_mut().enumerate() {
                *value += frozen[state].xyz[0][channel] - before[state][channel];
            }
            render_scene_exposure_range(&exposure, cx, &settings, &mut actual, sample, sample + 1)
                .unwrap();
            // Frozen scene renders supply the whole path and its actual light
            // proposal. Only accumulated-sum differencing introduces roundoff.
            assert_sum(actual.film().xyz[0], expected, sample);
        }
        assert_eq!(counts, [16, 48]);
        assert!(expected[1] > 0.0);
        assert_ne!(frozen[0].xyz, frozen[1].xyz);
    });
}

#[test]
fn g3_discontinuous_topology_changes_visibility_at_the_exact_accepted_boundary() {
    with_cx(|cx| {
        let clear = illuminated(2500.0, 1.0);
        let mut blocked = illuminated(2500.0, 1.0);
        blocked.primitives.push(quad(
            0.0,
            1.0,
            2.0,
            Material::Lambertian {
                reflectance: lift_rgb([0.0; 3]),
            },
        ));
        let states = frames(&clear, &blocked, 0.5);
        let settings = Settings {
            spp: 32,
            max_depth: 1,
            ..settings(DirectStrategy::NeeOnly)
        };
        for (time, selected) in [(0.0, &clear), (0.5, &blocked), (1.0, &blocked)] {
            let shutter = shutter(time, time);
            let exposure = exposure(cx, &states, shutter);
            let image = render_scene_exposure(&exposure, cx, &settings).unwrap();
            assert_eq!(
                image,
                render_motion(selected, cx, &settings, shutter).unwrap()
            );
            if time == 0.0 {
                assert!(image.xyz[0][1] > 0.0);
            } else {
                assert_eq!(image.xyz, vec![[0.0; 3]]);
            }
        }
    });
}

#[cfg(feature = "volumes")]
#[test]
fn g1_scene_state_moving_source_and_medium_history_share_the_same_sample_time() {
    use fs_render::animated_instances::{
        AnimatedGeometryInstance, RigidTransformTrajectory, TransformKeyframe,
    };
    use fs_render::instances::{RigidTransform, SharedGeometry};
    use fs_render::tracer::{MediumExposure, MediumExposureFilm, render_medium_exposure_range};
    use fs_render::volumes::VolumeGrid;
    use fs_render::volumes::transport::timeline::{
        MediumFrame, MediumTimeline, MediumTimelineBudget, PhysicalMedium,
    };
    use fs_render::volumes::transport::{
        AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum,
    };

    with_cx(|cx| {
        let moving = |temperature| {
            let mut s = illuminated(temperature, 1.0);
            let Shape::Mesh(mesh) = quad(0.0, 0.0, 1.0, s.primitives[1].material).shape else {
                unreachable!();
            };
            let trajectory = RigidTransformTrajectory::try_new(vec![
                TransformKeyframe::try_new(
                    0.0,
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-0.6, 0.0, 2.0]).unwrap(),
                    [1.2, 0.0, 0.4],
                )
                .unwrap(),
                TransformKeyframe::try_new(
                    1.0,
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.6, 0.0, 2.4]).unwrap(),
                    [1.2, 0.0, 0.4],
                )
                .unwrap(),
            ])
            .unwrap();
            s.primitives[1].shape = Shape::AnimatedInstance(
                AnimatedGeometryInstance::try_new(
                    541,
                    ContentHash([54; 32]),
                    SharedGeometry::mesh(mesh),
                    trajectory,
                )
                .unwrap(),
            );
            s
        };
        let first = moving(2200.0);
        let second = moving(2800.0);
        let states = frames(&first, &second, 0.25);
        let density = [1.0];
        let make_grid = |sigma| {
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-6.0, -4.0, 0.0], [12.0, 8.0, 4.0]),
                MassAbsorptionSpectrum::constant(sigma).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 1,
                    max_steps_per_segment: 8,
                },
            )
            .unwrap()
        };
        let thin = make_grid(0.1);
        let dense = make_grid(0.8);
        let medium_frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 0.5,
                medium: PhysicalMedium::Absorbing(&thin),
            },
            MediumFrame {
                start_s: 0.5,
                end_s: 1.0,
                medium: PhysicalMedium::Absorbing(&dense),
            },
        ];
        let history =
            MediumTimeline::try_new(cx, &medium_frames, MediumTimelineBudget { max_frames: 2 })
                .unwrap();
        let shutter = shutter(0.0, 1.0);
        let exposure = exposure(cx, &states, shutter)
            .with_medium_history(cx, &history)
            .unwrap();
        let medium_exposure = MediumExposure::try_new(cx, &history, shutter).unwrap();
        let settings = Settings {
            spp: 64,
            max_depth: 2,
            sampler: Sampler::OwenSobolFullPath,
            ..settings(DirectStrategy::PowerMis)
        };
        let mut frozen = [
            MediumExposureFilm::try_new(&medium_exposure, &settings).unwrap(),
            MediumExposureFilm::try_new(&medium_exposure, &settings).unwrap(),
        ];
        let mut actual = SceneExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut expected = [0.0; 3];
        for sample in 0..settings.spp {
            let before = [frozen[0].film().xyz[0], frozen[1].film().xyz[0]];
            for (state, s) in [&first, &second].into_iter().enumerate() {
                render_medium_exposure_range(
                    s,
                    &medium_exposure,
                    cx,
                    &settings,
                    &mut frozen[state],
                    sample,
                    sample + 1,
                )
                .unwrap();
            }
            let time =
                shutter.time_at(shutter.sample_for_stream(settings.seed, 0, u64::from(sample)));
            let state = usize::from(time >= 0.25);
            for (channel, value) in expected.iter_mut().enumerate() {
                *value += frozen[state].film().xyz[0][channel] - before[state][channel];
            }
            render_scene_exposure_range(&exposure, cx, &settings, &mut actual, sample, sample + 1)
                .unwrap();
            assert_sum(actual.film().xyz[0], expected, sample);
        }
        assert!(expected[1] > 0.0);
    });
}

#[test]
fn g5_cinematic_scene_exposure_replays_across_tiles_workers_and_progressive_ranges() {
    with_cx(|cx| {
        let mut first = illuminated(2200.0, 0.5);
        let mut second = illuminated(2800.0, 0.8);
        first.camera.half_tan = 0.2;
        second.camera.half_tan = 0.2;
        let physical = |x| {
            PhysicalCamera::try_legacy_compatible(
                Point3::new(x, 0.0, 0.5),
                first.camera.forward,
                first.camera.up,
                first.camera.half_tan,
                1.5,
                Aperture::try_circular(0.03).unwrap(),
            )
            .unwrap()
        };
        let camera = AnimatedCamera::try_new(vec![
            CameraShot::try_new(
                540,
                0.0,
                1.0,
                vec![
                    CameraKeyframe::try_new(0.0, physical(-0.4)).unwrap(),
                    CameraKeyframe::try_new(1.0, physical(0.4)).unwrap(),
                ],
            )
            .unwrap(),
        ])
        .unwrap();
        let states = frames(&first, &second, 0.25);
        let constant = [SceneFrame {
            start_s: 0.0,
            end_s: 1.0,
            scene: &first,
        }];
        let shutter = shutter(0.0, 1.0);
        let varying = exposure(cx, &states, shutter)
            .with_camera(cx, &camera, CutSide::After)
            .unwrap();
        let fixed_scene = exposure(cx, &constant, shutter)
            .with_camera(cx, &camera, CutSide::After)
            .unwrap();
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let settings = Settings {
                width: 3,
                height: 2,
                spp: 32,
                sampler,
                ..settings(DirectStrategy::PowerMis)
            };
            assert_eq!(
                render_scene_exposure(&fixed_scene, cx, &settings).unwrap(),
                render_cinematic(&first, &camera, CutSide::After, cx, &settings, shutter).unwrap()
            );
            let serial = render_scene_exposure(&varying, cx, &settings).unwrap();
            for config in [execution(1, 1, 1), execution(2, 1, 3)] {
                assert_eq!(
                    render_scene_exposure_with_execution(&varying, cx, &settings, &config)
                        .unwrap()
                        .film,
                    serial
                );
            }
            let mut resumed = SceneExposureFilm::try_new(&varying, &settings).unwrap();
            render_scene_exposure_range(&varying, cx, &settings, &mut resumed, 0, 7).unwrap();
            render_scene_exposure_range_with_execution(
                &varying,
                cx,
                &settings,
                &mut resumed,
                7,
                19,
                &execution(2, 1, 3),
            )
            .unwrap();
            render_scene_exposure_range(&varying, cx, &settings, &mut resumed, 19, settings.spp)
                .unwrap();
            assert_eq!(*resumed.film(), serial);
        }
    });
}

#[test]
fn g4_changed_scene_exposure_settings_and_cancellation_never_publish_a_partial_append() {
    with_cx(|cx| {
        let first = illuminated(2200.0, 0.5);
        let second = illuminated(2800.0, 0.8);
        let states = frames(&first, &second, 0.25);
        let changed_states = frames(&first, &second, 0.5);
        let accepted = exposure(cx, &states, shutter(0.0, 1.0));
        let changed_history = exposure(cx, &changed_states, shutter(0.0, 1.0));
        let changed_shutter = exposure(cx, &states, shutter(0.0, 0.75));
        let settings = Settings {
            spp: 24,
            ..settings(DirectStrategy::PowerMis)
        };
        let mut film = SceneExposureFilm::try_new(&accepted, &settings).unwrap();
        render_scene_exposure_range(&accepted, cx, &settings, &mut film, 0, 7).unwrap();
        let before = film.clone();
        for invalid in [&changed_history, &changed_shutter] {
            assert!(render_scene_exposure_range(invalid, cx, &settings, &mut film, 7, 24).is_err());
            assert_eq!(film, before);
            assert!(
                render_scene_exposure_range_with_execution(
                    invalid,
                    cx,
                    &settings,
                    &mut film,
                    7,
                    24,
                    &execution(1, 1, 2)
                )
                .is_err()
            );
            assert_eq!(film, before);
        }
        let changed_settings = Settings {
            seed: settings.seed + 1,
            ..settings
        };
        assert!(
            render_scene_exposure_range(&accepted, cx, &changed_settings, &mut film, 7, 24)
                .is_err()
        );
        assert_eq!(film, before);
        let gate = CancelGate::new_clock_free();
        let cancelled = Cx::new(
            &gate,
            cx.arena(),
            cx.stream_key(),
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        gate.request();
        assert_eq!(
            render_scene_exposure_range(&accepted, &cancelled, &settings, &mut film, 7, 24),
            Err(TracerError::Cancelled)
        );
        assert!(matches!(
            render_scene_exposure_range_with_execution(
                &accepted,
                &cancelled,
                &settings,
                &mut film,
                7,
                24,
                &execution(1, 1, 2)
            ),
            Err(RenderExecutionError::Tracer(TracerError::Cancelled))
        ));
        assert_eq!(film, before);
    });
}
