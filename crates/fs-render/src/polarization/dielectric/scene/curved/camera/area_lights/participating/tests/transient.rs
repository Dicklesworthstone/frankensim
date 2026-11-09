//! Retarded histories on actual molecular vertices and refracted area images.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::gradient_index::bundle::{IndexHessian, IndexJet, SecondOrderIndexField};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::TransientAreaLightCameraError;
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::participating::TransientCurvedRayleighAreaLightCameraScene;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use fs_math::dd::Dd;
use std::sync::Mutex;

const C: f64 = 299_792_458.0;
const AREA_DRAW: [f64; 2] = [4.0 / 9.0, 0.5];

struct History<F> {
    evaluate: F,
    queries: Mutex<Vec<SurfaceEmissionQuery>>,
}

impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            queries: Mutex::new(Vec::new()),
        }
    }
}
impl<F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync> SurfaceEmissionHistory
    for History<F>
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.evaluate)(query)
    }
}
fn history_binding(history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: 9,
        source_identity: ContentHash([93; 32]),
        history,
    }
}

fn physical_camera() -> PhysicalCamera {
    let lens_glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(
        SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(lens_glass)).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
                LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_molecular_histories_retain_null_candidates_and_independent_source_survival() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(0)]).unwrap();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(LE).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let seed = seeds(Vec3::new(-1.0, 0.0, 0.0));
        let bindings = [binding(&source, &seed)];
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new_heterogeneous(
                &spectrum,
                &HalfDensity,
                1.0,
                ContentHash([83; 32]),
                &phase,
            )
            .unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &bindings,
        )
        .unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(20e-9);
        let history = History::new(|_| Ok(2.0));
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer =
            TransientCurvedRayleighAreaLightCameraScene::try_new(cx, &timed, &owner).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let flight = |draw| match draw {
            0 => [1.0 - (-0.5_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - (-0.25_f64).exp(), 0.1, 0.4, 0.9],
            2 => [0.99, 0.5, 0.0, 0.9],
            _ => panic!("unexpected camera flight"),
        };
        for candidate in [false, true] {
            let source_draw = |area, image, draw| {
                assert_eq!((area, image), (0, 0));
                if candidate && draw == 0 {
                    1.0 - (-0.25_f64).exp()
                } else {
                    0.99
                }
            };
            let baseline = owner
                .sample_ray_with_tracking_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("opaque emitter"),
                    flight,
                    |_, _| panic!("no surface receiver"),
                    |_, _| AREA_DRAW,
                    source_draw,
                )
                .unwrap();
            let mut flights = Vec::new();
            let mut area_draws = Vec::new();
            let mut source_draws = Vec::new();
            let result = observer
                .sample_ray_with_tracking_samples_at_time(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    clock,
                    settings(),
                    &mut workspace,
                    |_| panic!("opaque emitter"),
                    |draw| {
                        flights.push(draw);
                        flight(draw)
                    },
                    |_, _| panic!("no surface receiver"),
                    |event, binding| {
                        area_draws.push((event, binding));
                        AREA_DRAW
                    },
                    |area, image, draw| {
                        source_draws.push((area, image, draw));
                        source_draw(area, image, draw)
                    },
                )
                .unwrap();
            let area = result.area_sample();
            near(
                area.area_radiance_w_m2_sr_nm(),
                2.0 * baseline.area_radiance_w_m2_sr_nm(),
            );
            near(
                area.radiance_sample().radiance_w_m2_sr_nm(),
                2.0 * baseline.radiance_sample().radiance_w_m2_sr_nm(),
            );
            assert_eq!(
                area.radiance_sample().path(),
                baseline.radiance_sample().path()
            );
            assert_eq!((area.volume_events(), area.free_flight_draws()), (1, 3));
            assert_eq!(area.trace_evaluations(), baseline.trace_evaluations());
            assert_eq!(area.field_evaluations(), baseline.field_evaluations());
            assert_eq!(area.query_work(), baseline.query_work());
            assert_eq!(
                area.source_free_flight_draws(),
                baseline.source_free_flight_draws()
            );
            assert_eq!(area.source_candidates(), u32::from(candidate));
            assert_eq!(flights, [0, 1, 2]);
            assert_eq!(area_draws, [(0, 0)]);
            assert_eq!(
                source_draws.len(),
                baseline.source_free_flight_draws() as usize
            );
            let queries = history.queries.lock().unwrap();
            let pair = &queries[queries.len() - 2..];
            for query in pair {
                assert_eq!(query.observation_time_s, clock);
                near(query.group_delay_s.to_f64(), 1.4 * 1.5 / C);
                near(
                    (query.emission_time_s - epoch).to_f64(),
                    20e-9 - 1.4 * 1.5 / C,
                );
            }
            assert_eq!(
                pair[0].position_m,
                source.sample_area(AREA_DRAW).unwrap().position_m()
            );
        }
        assert_eq!(
            history.queries.lock().unwrap().len(),
            4,
            "null candidates do not query a history"
        );
    });
}

#[test]
fn g1_timed_lens_molecular_images_keep_distinct_clocks_and_integrate_the_source_shutter() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![panel(0)]).unwrap();
        let material = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(LE).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let seed = seeds(Vec3::new(-1.0, 0.0, 0.0));
        let bindings = [binding(&source, &seed)];
        let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
            cx, &molecular, &surface, &bindings,
        )
        .unwrap();
        let physical = physical_camera();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(20e-9);
        assert_eq!(clock.to_f64(), epoch.to_f64());
        let history = History::new(|_| Ok(1.0));
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer =
            TransientCurvedRayleighAreaLightCameraScene::try_new(cx, &timed, &owner).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let flight = |draw| match draw {
            0 => [1.0 - (-0.35_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("opaque emitter precedes another collision"),
        };
        for wavelength in [450.0, 650.0] {
            let beam = physical
                .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[wavelength],
                    analyzer.reference_axis(),
                    analyzer.row(),
                )
                .unwrap()[0]
                .unwrap();
            assert!(beam.power().q != 0.0 || beam.power().u != 0.0);
            let baseline = owner
                .sample_sequential_ray_with_tracking_samples(
                    cx,
                    &beam,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("opaque emitter"),
                    flight,
                    |_, _| panic!("no surface receiver"),
                    |_, _| AREA_DRAW,
                    |_, _, _| panic!("uniform survival is analytic"),
                )
                .unwrap();
            let result = observer
                .sample_sequential_ray_with_tracking_samples_at_time(
                    cx,
                    &beam,
                    0,
                    clock,
                    settings(),
                    &mut workspace,
                    |_| panic!("opaque emitter"),
                    flight,
                    |_, _| panic!("no surface receiver"),
                    |_, _| AREA_DRAW,
                    |_, _, _| panic!("uniform survival is analytic"),
                )
                .unwrap();
            assert_eq!(result.area_sample(), &baseline);
            assert!(baseline.area_radiance_w_m2_sr_nm() > 0.0);
            assert!(
                baseline.radiance_sample().radiance_w_m2_sr_nm()
                    > baseline.area_radiance_w_m2_sr_nm()
            );
            let ng = 1.48 + 0.018 / (wavelength * 1e-3).powi(2);
            let lens_delay = (0.09 + 0.01 * ng) / C;
            let queries = history.queries.lock().unwrap();
            let pair = &queries[queries.len() - 2..];
            // The first collision is at z=.6, the area sample at (1,0,.5),
            // and the sampled continuation hits (1,0,.6). Their clocks differ.
            near(
                pair[0].group_delay_s.to_f64(),
                lens_delay + (0.5 + 1.01_f64.sqrt()) / C,
            );
            near(pair[1].group_delay_s.to_f64(), lens_delay + 1.5 / C);
            assert_ne!(pair[0].group_delay_s, pair[1].group_delay_s);
            for query in pair {
                assert_eq!(query.observation_time_s, clock);
            }
        }

        let beam = physical
            .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
                analyzer.reference_axis(),
                analyzer.row(),
            )
            .unwrap()[0]
            .unwrap();
        let baseline = owner
            .sample_sequential_ray_with_tracking_samples(
                cx,
                &beam,
                0,
                settings(),
                &mut workspace,
                |_| panic!("opaque emitter"),
                flight,
                |_, _| panic!("no surface receiver"),
                |_, _| AREA_DRAW,
                |_, _, _| panic!("uniform survival is analytic"),
            )
            .unwrap();
        let pulse = PiecewiseLinearHistory::try_new(
            cx,
            epoch,
            &[
                (-1e-6, 0.0),
                (0.0, 0.0),
                (1e-9, 2.0),
                (2e-9, 0.0),
                (1e-6, 0.0),
            ],
        )
        .unwrap();
        let pulse_camera =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&pulse)]).unwrap();
        let pulse_observer =
            TransientCurvedRayleighAreaLightCameraScene::try_new(cx, &pulse_camera, &owner)
                .unwrap();
        let arrival =
            epoch + Dd::from_f64((0.09 + 0.01 * (1.48 + 0.018 / 0.55_f64.powi(2)) + 1.5) / C);
        let window = ExposureWindow::try_new(
            arrival - Dd::from_f64(0.5e-9),
            arrival + Dd::from_f64(2.5e-9),
        )
        .unwrap();
        let exposed = pulse_observer
            .sample_sequential_ray_with_tracking_samples_during_exposure(
                cx,
                &beam,
                0,
                window,
                settings(),
                &mut workspace,
                |_| panic!("opaque emitter"),
                flight,
                |_, _| panic!("no surface receiver"),
                |_, _| AREA_DRAW,
                |_, _, _| panic!("uniform survival is analytic"),
            )
            .unwrap();
        near(
            exposed.mean_area_sample().area_radiance_w_m2_sr_nm(),
            baseline.area_radiance_w_m2_sr_nm() * 2.0 / 3.0,
        );
        near(
            exposed.mean_radiance_sample().radiance_w_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * 2.0 / 3.0,
        );
        near(
            exposed.integrated_radiance_j_m2_sr_nm(),
            baseline.radiance_sample().radiance_w_m2_sr_nm() * 2e-9,
        );
        assert_eq!(
            exposed.mean_area_sample().trace_evaluations(),
            baseline.trace_evaluations()
        );
        assert_eq!(
            exposed.mean_area_sample().free_flight_draws(),
            baseline.free_flight_draws()
        );
        assert_eq!(
            exposed.mean_radiance_sample().path(),
            baseline.radiance_sample().path()
        );

        let untimed = physical
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert!(
            owner
                .sample_sequential_ray_with_tracking_samples(
                    cx,
                    &untimed,
                    0,
                    settings(),
                    &mut workspace,
                    |_| panic!("opaque emitter"),
                    flight,
                    |_, _| panic!("no surface receiver"),
                    |_, _| AREA_DRAW,
                    |_, _, _| panic!("uniform source"),
                )
                .unwrap()
                .radiance_sample()
                .radiance_w_m2_sr_nm()
                > 0.0
        );
        assert_eq!(
            observer.sample_sequential_ray_with_tracking_samples_at_time(
                cx,
                &untimed,
                0,
                clock,
                settings(),
                &mut workspace,
                |_| panic!("untimed beam"),
                |_| panic!("untimed beam"),
                |_, _| panic!("untimed beam"),
                |_, _| panic!("untimed beam"),
                |_, _, _| panic!("untimed beam"),
            ),
            Err(TransientAreaLightCameraError::Source(
                TransientCameraError::MissingLensGroupDelay
            ))
        );
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &[])
            .unwrap()
            .with_analyzer(analyzer);
        let analyzed_observer =
            TransientCurvedRayleighAreaLightCameraScene::try_new(cx, &analyzed, &owner).unwrap();
        assert_eq!(
            analyzed_observer.sample_sequential_ray_with_tracking_samples_at_time(
                cx,
                &beam,
                0,
                clock,
                settings(),
                &mut workspace,
                |_| panic!("sensor analyzer"),
                |_| panic!("sensor analyzer"),
                |_, _| panic!("sensor analyzer"),
                |_, _| panic!("sensor analyzer"),
                |_, _, _| panic!("sensor analyzer"),
            ),
            Err(TransientAreaLightCameraError::Source(
                TransientCameraError::AnalyzerMustPrecedeLens
            ))
        );
    });
}

#[test]
fn g0_molecular_source_connections_refuse_missing_group_data_and_original_history_failures() {
    struct Field {
        missing_group: bool,
    }
    impl GradientIndexField for Field {
        fn sample(
            &self,
            position: Point3,
            _: f64,
        ) -> Result<GradientIndexSample, GradientIndexError> {
            GradientIndexSample::try_new(
                1.4,
                Vec3::new(0.0, 0.0, 0.0),
                0.2,
                if self.missing_group && position.x > 0.1 {
                    None
                } else {
                    Some(1.4)
                },
            )
        }
    }
    impl SecondOrderIndexField for Field {
        fn sample_second_order(
            &self,
            position: Point3,
            wavelength: f64,
        ) -> Result<IndexJet, GradientIndexError> {
            Ok(IndexJet::new(
                self.sample(position, wavelength)?,
                IndexHessian::ZERO,
            ))
        }
    }
    with_cx(|cx| {
        for missing_group in [false, true] {
            let field = Field { missing_group };
            let geometry = InstanceScene::try_new(vec![panel(0)]).unwrap();
            let media = [SceneMedium::Differentiable {
                field: &field,
                identity: ContentHash([95; 32]),
            }];
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                }],
                100_000,
            )
            .unwrap();
            let emitters = [SurfaceEmitter::new(
                9,
                RadianceSpectrum::try_constant(LE).unwrap(),
                EmissionSides::Front,
            )];
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            let surface = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
            let seed = seeds(Vec3::new(-1.0, 0.0, 0.0));
            let bindings = [binding(&source, &seed)];
            let spectrum = MassScatteringSpectrum::constant(0.7).unwrap();
            let phase = PhaseFunction::Rayleigh;
            let laws = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let owner = CurvedRayleighAreaLightCameraScene::try_new_with_volume_bindings(
                cx, &molecular, &surface, &bindings,
            )
            .unwrap();
            let history = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
            let timed =
                TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)])
                    .unwrap();
            let observer =
                TransientCurvedRayleighAreaLightCameraScene::try_new(cx, &timed, &owner).unwrap();
            let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let flight = |draw| {
                if draw == 0 {
                    [1.0 - (-0.35_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0]
                } else {
                    [1.0 - 1e-12, 0.5, 0.0, 0.0]
                }
            };
            let result = observer.sample_ray_with_tracking_samples_at_time(
                cx,
                camera_ray(),
                550.0,
                0,
                Dd::ZERO,
                settings(),
                &mut workspace,
                |_| panic!("opaque emitter"),
                flight,
                |_, _| panic!("no surface receiver"),
                |_, _| AREA_DRAW,
                |_, _, _| panic!("uniform source survival"),
            );
            let expected = if missing_group {
                TransientCameraError::MissingGroupDelay { object_id: 9 }
            } else {
                TransientCameraError::History {
                    object_id: 9,
                    error: SurfaceHistoryError::OutsideTimeDomain,
                }
            };
            assert_eq!(result, Err(TransientAreaLightCameraError::Source(expected)));
            assert_eq!(
                history.queries.lock().unwrap().len(),
                usize::from(!missing_group)
            );
        }
    });
}
