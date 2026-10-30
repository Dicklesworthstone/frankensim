//! Accepted geometry and participating-field owners share the physical shutter.

use super::*;
use fs_render::camera::{
    AnimatedCamera, Aperture, CameraKeyframe, CameraShot, CutSide, PhysicalCamera,
};
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::tracer::{
    SceneExposure, SceneExposureBudget, SceneExposureFilm, SceneFrame, render_scene_exposure,
    render_scene_exposure_range, render_scene_exposure_range_with_execution,
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

fn state_scenes() -> [Scene; 2] {
    let make_body = |index: f64, thickness: f64| {
        slab(
            0.0,
            -thickness,
            4.0,
            glass(
                CauchyIor::try_constant(index).unwrap(),
                80.0,
                DielectricSurface::SMOOTH,
            ),
        )
    };
    [
        empty_scene(vec![make_body(1.3, 0.2)]),
        empty_scene(vec![
            // A genuine primitive-list/topology change moves the glass owner
            // from index zero to one. Reusing the earlier binding is invalid.
            rectangle(
                Point3::new(20.0, -1.0, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 2.0, 0.0),
            ),
            make_body(1.5, 0.6),
        ]),
    ]
}

fn state_grids(cx: &Cx<'_>) -> [AbsorbingGrid<'static>; 2] {
    [
        field(
            cx,
            &[0.5],
            [-4.0, -4.0, -0.2],
            [8.0, 8.0, 0.2],
            MassAbsorptionSpectrum::constant(0.8).unwrap(),
            Some(&[1700.0]),
        ),
        field(
            cx,
            &[1.5],
            [-4.0, -4.0, -0.6],
            [8.0, 8.0, 0.6],
            MassAbsorptionSpectrum::constant(0.6).unwrap(),
            Some(&[2300.0]),
        ),
    ]
}

fn frames(scenes: &[Scene; 2]) -> [SceneFrame<'_>; 2] {
    [
        SceneFrame {
            start_s: 0.0,
            end_s: 0.25,
            scene: &scenes[0],
        },
        SceneFrame {
            start_s: 0.25,
            end_s: 1.0,
            scene: &scenes[1],
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
            max_frames: 2,
            max_total_primitives: 8,
            max_total_light_triangles: 16,
        },
    )
    .unwrap()
}

#[test]
fn g1_exposure_uses_the_selected_geometry_and_field_for_every_path_sample() {
    with_cx(|cx, _| {
        let scenes = state_scenes();
        let [first_grid, second_grid] = state_grids(cx);
        let second_grid = ScatteringGrid::try_new(
            cx,
            second_grid,
            MassScatteringSpectrum::constant(0.5).unwrap(),
            PhaseFunction::HenyeyGreenstein { g: 0.2 },
        )
        .unwrap();
        let ambient = [0.1, 0.4].map(|coefficient| {
            field(
                cx,
                &[1.0],
                [-4.0, -4.0, -1.0],
                [8.0, 8.0, 2.0],
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                None,
            )
        });
        let bindings = [
            [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(&first_grid),
            }],
            [InteriorMedium {
                boundary_primitive: 1,
                medium: PhysicalMedium::Scattering(&second_grid),
            }],
        ];
        let media = [
            SceneMedia::try_new(
                cx,
                &scenes[0],
                Some(PhysicalMedium::Absorbing(&ambient[0])),
                &bindings[0],
                media_budget(),
            )
            .unwrap(),
            SceneMedia::try_new(
                cx,
                &scenes[1],
                Some(PhysicalMedium::Absorbing(&ambient[1])),
                &bindings[1],
                media_budget(),
            )
            .unwrap(),
        ];
        let owners = [&media[0], &media[1]];
        let frames = frames(&scenes);
        let shutter = shutter(0.0, 1.0);
        let exposure = exposure(cx, &frames, shutter)
            .with_scene_media(cx, &owners)
            .unwrap();
        let settings = settings(DirectStrategy::PowerMis, 64, 6);
        let mut frozen = [
            SceneMediaFilm::try_new(&media[0], &settings).unwrap(),
            SceneMediaFilm::try_new(&media[1], &settings).unwrap(),
        ];
        let mut actual = SceneExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut expected = [0.0; 3];
        let mut counts = [0; 2];
        for sample in 0..settings.spp {
            let before = [frozen[0].film().xyz[0], frozen[1].film().xyz[0]];
            for state in 0..2 {
                render_scene_media_range(
                    &media[state],
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
            counts[state] += 1;
            for channel in 0..3 {
                expected[channel] += frozen[state].film().xyz[0][channel] - before[state][channel];
            }
            render_scene_exposure_range(&exposure, cx, &settings, &mut actual, sample, sample + 1)
                .unwrap();
            for channel in 0..3 {
                let got = actual.film().xyz[0][channel];
                assert!(
                    (got - expected[channel]).abs() <= 2.0e-12 * expected[channel].abs().max(1.0),
                    "sample {sample}, state {state}, channel {channel}: {got} vs {}",
                    expected[channel]
                );
            }
        }
        assert!(counts.iter().all(|count| *count > 0));
        assert!(expected[1] > 0.0);
        assert_ne!(frozen[0].film().xyz, frozen[1].film().xyz);
    });
}

#[test]
fn g1_zero_duration_exposure_selects_later_and_final_medium_owners() {
    with_cx(|cx, _| {
        let scenes = state_scenes();
        let grids = state_grids(cx);
        let bindings = [
            [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(&grids[0]),
            }],
            [InteriorMedium {
                boundary_primitive: 1,
                medium: PhysicalMedium::Absorbing(&grids[1]),
            }],
        ];
        let media = [
            SceneMedia::try_new(cx, &scenes[0], None, &bindings[0], media_budget()).unwrap(),
            SceneMedia::try_new(cx, &scenes[1], None, &bindings[1], media_budget()).unwrap(),
        ];
        let owners = [&media[0], &media[1]];
        let frames = frames(&scenes);
        let settings = settings(DirectStrategy::BsdfOnly, 32, 3);
        let expected = render_scene_media(&media[1], cx, &settings).unwrap();
        for endpoint in [0.25, 1.0] {
            let exposure = exposure(cx, &frames, shutter(endpoint, endpoint))
                .with_scene_media(cx, &owners)
                .unwrap();
            let actual = render_scene_exposure(&exposure, cx, &settings).unwrap();
            assert_eq!(actual.xyz, expected.xyz, "accepted endpoint {endpoint}");
            assert_eq!(actual.spp_done, expected.spp_done);
        }
    });
}

#[test]
fn g5_cinematic_exposure_with_interior_fields_replays_across_tiles_and_ranges() {
    with_cx(|cx, _| {
        let scenes = state_scenes();
        let [first_grid, second_grid] = state_grids(cx);
        let second_grid = ScatteringGrid::try_new(
            cx,
            second_grid,
            MassScatteringSpectrum::constant(0.5).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let bindings = [
            [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(&first_grid),
            }],
            [InteriorMedium {
                boundary_primitive: 1,
                medium: PhysicalMedium::Scattering(&second_grid),
            }],
        ];
        let media = [
            SceneMedia::try_new(cx, &scenes[0], None, &bindings[0], media_budget()).unwrap(),
            SceneMedia::try_new(cx, &scenes[1], None, &bindings[1], media_budget()).unwrap(),
        ];
        let owners = [&media[0], &media[1]];
        let frames = frames(&scenes);
        let physical = |x| {
            PhysicalCamera::try_legacy_compatible(
                Point3::new(x, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -1.0),
                Vec3::new(0.0, 1.0, 0.0),
                0.2,
                1.5,
                Aperture::try_circular(0.02).unwrap(),
            )
            .unwrap()
        };
        let camera = AnimatedCamera::try_new(vec![
            CameraShot::try_new(
                713,
                0.0,
                1.0,
                vec![
                    CameraKeyframe::try_new(0.0, physical(-0.2)).unwrap(),
                    CameraKeyframe::try_new(1.0, physical(0.2)).unwrap(),
                ],
            )
            .unwrap(),
        ])
        .unwrap();
        let exposure = exposure(cx, &frames, shutter(0.0, 1.0))
            .with_scene_media(cx, &owners)
            .unwrap()
            .with_camera(cx, &camera, CutSide::After)
            .unwrap();
        let execution = RenderExecutionConfig::try_new(2, 1, 3, 8 << 20, RunId(714)).unwrap();
        let mut settings = settings(DirectStrategy::PowerMis, 24, 6);
        settings.width = 3;
        settings.height = 2;
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            settings.sampler = sampler;
            let serial = render_scene_exposure(&exposure, cx, &settings).unwrap();
            assert_eq!(
                render_scene_exposure_with_execution(&exposure, cx, &settings, &execution)
                    .unwrap()
                    .film,
                serial
            );
            let mut progressive = SceneExposureFilm::try_new(&exposure, &settings).unwrap();
            render_scene_exposure_range(&exposure, cx, &settings, &mut progressive, 0, 7).unwrap();
            render_scene_exposure_range_with_execution(
                &exposure,
                cx,
                &settings,
                &mut progressive,
                7,
                19,
                &execution,
            )
            .unwrap();
            render_scene_exposure_range(&exposure, cx, &settings, &mut progressive, 19, 24)
                .unwrap();
            assert_eq!(progressive.film(), &serial);
        }
    });
}

#[test]
fn g4_changed_exposure_field_and_cancellation_preserve_the_completed_prefix() {
    with_cx(|cx, gate| {
        let scenes = state_scenes();
        let grids = state_grids(cx);
        let changed_grid = field(
            cx,
            &[3.0],
            [-4.0, -4.0, -0.6],
            [8.0, 8.0, 0.6],
            MassAbsorptionSpectrum::constant(0.6).unwrap(),
            Some(&[2300.0]),
        );
        let bindings = [
            [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(&grids[0]),
            }],
            [InteriorMedium {
                boundary_primitive: 1,
                medium: PhysicalMedium::Absorbing(&grids[1]),
            }],
            [InteriorMedium {
                boundary_primitive: 1,
                medium: PhysicalMedium::Absorbing(&changed_grid),
            }],
        ];
        let media = [
            SceneMedia::try_new(cx, &scenes[0], None, &bindings[0], media_budget()).unwrap(),
            SceneMedia::try_new(cx, &scenes[1], None, &bindings[1], media_budget()).unwrap(),
            SceneMedia::try_new(cx, &scenes[1], None, &bindings[2], media_budget()).unwrap(),
        ];
        let owners = [&media[0], &media[1]];
        let changed_owners = [&media[0], &media[2]];
        let frames = frames(&scenes);
        let accepted = exposure(cx, &frames, shutter(0.0, 1.0))
            .with_scene_media(cx, &owners)
            .unwrap();
        let changed = exposure(cx, &frames, shutter(0.0, 1.0))
            .with_scene_media(cx, &changed_owners)
            .unwrap();
        let settings = settings(DirectStrategy::PowerMis, 24, 4);
        let mut film = SceneExposureFilm::try_new(&accepted, &settings).unwrap();
        render_scene_exposure_range(&accepted, cx, &settings, &mut film, 0, 7).unwrap();
        let before = film.clone();
        assert!(render_scene_exposure_range(&changed, cx, &settings, &mut film, 7, 24).is_err());
        assert_eq!(film, before);
        let execution = RenderExecutionConfig::try_new(1, 1, 2, 8 << 20, RunId(715)).unwrap();
        assert!(
            render_scene_exposure_range_with_execution(
                &changed, cx, &settings, &mut film, 7, 24, &execution
            )
            .is_err()
        );
        assert_eq!(film, before);
        gate.request();
        assert!(matches!(
            render_scene_exposure_range(&accepted, cx, &settings, &mut film, 7, 24),
            Err(TracerError::Cancelled)
        ));
        assert_eq!(film, before);
    });
}
