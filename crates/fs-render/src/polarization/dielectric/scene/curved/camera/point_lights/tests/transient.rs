//! Exact source identity and native group clocks for calibrated point fixtures.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, piecewise_linear::PiecewiseLinearHistory,
};
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::*;
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use fs_math::dd::Dd;
use std::sync::Mutex;

const C: f64 = 299_792_458.0;
const ID: ContentHash = ContentHash([123; 32]);

struct History<'a> {
    source: &'a PhotometricSource,
    epoch: Dd,
    offset: f64,
    slope: f64,
    calls: Mutex<Vec<(f64, Vec3, Dd, Dd, Dd)>>,
}
impl PointEmissionHistory for History<'_> {
    fn sample(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        assert!(core::ptr::eq(query.source, self.source));
        assert_eq!(query.source_identity, ID);
        assert_eq!(query.position_m, self.source.position_m());
        self.calls.lock().unwrap().push((
            query.wavelength_nm,
            query.outgoing_direction,
            query.observation_time_s,
            query.emission_time_s,
            query.group_delay_s,
        ));
        Ok(self.offset + self.slope * (query.emission_time_s - self.epoch).to_f64())
    }
}
fn history_binding<'a>(
    source: &'a PhotometricSource,
    history: &'a dyn PointEmissionHistory,
) -> PointEmissionHistoryBinding<'a> {
    PointEmissionHistoryBinding {
        source,
        source_identity: ID,
        history,
    }
}
struct SurfaceHistory {
    epoch: Dd,
}
impl SurfaceEmissionHistory for SurfaceHistory {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        Ok(1.0 + 2e7 * (query.emission_time_s - self.epoch).to_f64())
    }
}

#[test]
fn g1_point_and_surface_histories_keep_distinct_original_source_clocks_and_work() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 2.0), plane(9, -1.0, 0.0, 2.0)])
                .unwrap();
        let material = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.006, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9, 2.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let lamp = source(Point3::new(0.6, 0.0, 0.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
        let lights = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = History {
            source: &lamp,
            epoch,
            offset: 1.0,
            slope: 1e7,
            calls: Mutex::new(Vec::new()),
        };
        let surface = SurfaceHistory { epoch };
        let timed = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[SurfaceEmissionHistoryBinding {
                object_id: 9,
                source_identity: ContentHash([124; 32]),
                history: &surface,
            }],
        )
        .unwrap();
        let histories = [history_binding(&lamp, &history)];
        let observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &histories).unwrap();
        for wavelength in [450.0_f64, 650.0] {
            let baseline = lights
                .sample_ray_with_samples(cx, ray(0.0), wavelength, 0, settings(), |_| [0.0; 3])
                .unwrap();
            let mut draws = Vec::new();
            let result = observer
                .sample_ray_with_samples_at_time(
                    cx,
                    ray(0.0),
                    wavelength,
                    0,
                    clock,
                    settings(),
                    |event| {
                        draws.push(event);
                        [0.0; 3]
                    },
                )
                .unwrap();
            let group_index = 1.4 + 0.018 / (wavelength * 1e-3).powi(2);
            let distance = 1.36_f64.sqrt();
            let point_delay = group_index * (1.0 + distance) / C;
            let surface_delay = group_index * 3.0 / C;
            let point_scale = 1.0 + 1e7 * (30e-9 - point_delay);
            let surface_scale = 1.0 + 2e7 * (30e-9 - surface_delay);
            near(
                result.point_light_sample().fixture_radiance_w_m2_sr_nm(),
                baseline.fixture_radiance_w_m2_sr_nm() * point_scale,
            );
            near(
                result.radiance_sample().radiance_w_m2_sr_nm(),
                baseline.fixture_radiance_w_m2_sr_nm() * point_scale
                    + (baseline.radiance_sample().radiance_w_m2_sr_nm()
                        - baseline.fixture_radiance_w_m2_sr_nm())
                        * surface_scale,
            );
            assert_eq!(draws, [0]);
            assert_eq!(
                result.radiance_sample().path(),
                baseline.radiance_sample().path()
            );
            assert_eq!(
                result.point_light_sample().connections(),
                baseline.connections()
            );
            assert_eq!(
                result.point_light_sample().trace_evaluations(),
                baseline.trace_evaluations()
            );
            assert_eq!(
                result.point_light_sample().field_evaluations(),
                baseline.field_evaluations()
            );
            assert_eq!(
                result.point_light_sample().query_work(),
                baseline.query_work()
            );
            let queries = history.calls.lock().unwrap();
            let query = queries.last().unwrap();
            assert_eq!(query.0, wavelength);
            near(query.1.x, -0.6 / distance);
            near(query.1.y, 0.0);
            near(query.1.z, 1.0 / distance);
            assert_eq!(query.2, clock);
            assert_eq!(query.3, clock - query.4);
            near(query.4.to_f64(), point_delay);
        }
        let clone = lamp.clone();
        let cloned_binding = [history_binding(&clone, &history)];
        assert!(matches!(
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &cloned_binding),
            Err(TransientPointLightCameraError::InvalidHistoryBinding { index: 0 })
        ));
        let duplicate = [histories[0], histories[0]];
        assert!(matches!(
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &duplicate),
            Err(TransientPointLightCameraError::InvalidHistoryBinding { index: 1 })
        ));
        let no_identity = [PointEmissionHistoryBinding {
            source_identity: ContentHash([0; 32]),
            ..histories[0]
        }];
        assert!(matches!(
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &no_identity),
            Err(TransientPointLightCameraError::MissingHistoryIdentity { index: 0 })
        ));
    });
}

fn physical_camera() -> PhysicalCamera {
    let lens_glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(lens_glass)).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
            LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

#[test]
fn g1_real_lens_point_history_keeps_sensor_prefix_and_integrates_narrow_source_pulse() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 2.0), plane(9, -1.0, 0.0, 2.0)])
                .unwrap();
        let material = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let lamp = source(Point3::new(0.6, 0.0, 0.0));
        let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
        let lights = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let history = History {
            source: &lamp,
            epoch,
            offset: 1.0,
            slope: 0.0,
            calls: Mutex::new(Vec::new()),
        };
        let histories = [history_binding(&lamp, &history)];
        let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
        let observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &histories).unwrap();
        let physical = physical_camera();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let beams = physical
            .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[450.0, 650.0],
                analyzer.reference_axis(),
                analyzer.row(),
            )
            .unwrap();
        let knots = [
            (-10e-9, 0.0),
            (0.0, 0.0),
            (1e-9, 2.0),
            (2e-9, 0.0),
            (10e-9, 0.0),
        ];
        let pulse = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
        let pulses = [history_binding(&lamp, &pulse)];
        let pulse_observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &pulses).unwrap();
        for optical in beams.into_iter().flatten() {
            assert!(optical.power().q != 0.0 || optical.power().u != 0.0);
            let baseline = lights
                .sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| [0.0; 3])
                .unwrap();
            let clock = epoch + Dd::from_f64(30e-9);
            let result = observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &optical,
                    0,
                    clock,
                    settings(),
                    |_| [0.0; 3],
                )
                .unwrap();
            assert_eq!(result.point_light_sample(), &baseline);
            let delay = optical.group_delay_s().unwrap()
                + Dd::from_f64(0.9 + 1.36_f64.sqrt()) / Dd::from_f64(C);
            let queries = history.calls.lock().unwrap();
            let query = queries.last().unwrap();
            near(query.4.to_f64(), delay.to_f64());
            assert_eq!(query.2, clock);
            assert_eq!(query.3, clock - query.4);
            drop(queries);
            let dark = pulse_observer
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &optical,
                    0,
                    epoch + delay - Dd::from_f64(2e-9),
                    settings(),
                    |_| [0.0; 3],
                )
                .unwrap();
            assert_eq!(dark.radiance_sample().radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(
                dark.point_light_sample().trace_evaluations(),
                baseline.trace_evaluations()
            );
            let window = ExposureWindow::try_new(
                epoch + delay - Dd::from_f64(0.5e-9),
                epoch + delay + Dd::from_f64(2.5e-9),
            )
            .unwrap();
            let exposed = pulse_observer
                .sample_sequential_ray_with_samples_during_exposure(
                    cx,
                    &optical,
                    0,
                    window,
                    settings(),
                    |_| [0.0; 3],
                )
                .unwrap();
            near(
                exposed.mean_radiance_sample().radiance_w_m2_sr_nm(),
                baseline.radiance_sample().radiance_w_m2_sr_nm() * (2.0 / 3.0),
            );
            near(
                exposed.integrated_radiance_j_m2_sr_nm(),
                baseline.radiance_sample().radiance_w_m2_sr_nm() * 2e-9,
            );
            assert_eq!(
                exposed.mean_point_light_sample().trace_evaluations(),
                baseline.trace_evaluations()
            );
            assert_eq!(
                exposed.mean_radiance_sample().path(),
                baseline.radiance_sample().path()
            );
        }
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
            lights
                .sample_sequential_ray_with_samples(cx, &untimed, 0, settings(), |_| [0.0; 3])
                .unwrap()
                .fixture_radiance_w_m2_sr_nm()
                > 0.0
        );
        assert!(matches!(
            observer.sample_sequential_ray_with_samples_at_time(
                cx,
                &untimed,
                0,
                epoch,
                settings(),
                |_| panic!("unqualified timing must refuse before camera draws")
            ),
            Err(TransientPointLightCameraError::Source(
                TransientCameraError::MissingLensGroupDelay
            ))
        ));
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &[])
            .unwrap()
            .with_analyzer(analyzer);
        let analyzed_observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &analyzed, &lights, &histories)
                .unwrap();
        assert!(matches!(
            analyzed_observer.sample_sequential_ray_with_samples_at_time(
                cx,
                &beams[0].unwrap(),
                0,
                epoch,
                settings(),
                |_| panic!("analyzer must precede lens")
            ),
            Err(TransientPointLightCameraError::Source(
                TransientCameraError::AnalyzerMustPrecedeLens
            ))
        ));
    });
}

#[test]
fn g0_point_connection_group_data_and_original_history_errors_refuse_without_phase_fallback() {
    use crate::gradient_index::bundle::{IndexJet, SecondOrderIndexField};
    use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
    use std::sync::atomic::{AtomicUsize, Ordering};
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
    struct Failure(AtomicUsize);
    impl PointEmissionHistory for Failure {
        fn sample(
            &self,
            cx: &Cx<'_>,
            _: PointEmissionQuery<'_>,
        ) -> Result<f64, SurfaceHistoryError> {
            cx.checkpoint()?;
            self.0.fetch_add(1, Ordering::Relaxed);
            Err(SurfaceHistoryError::OutsideTimeDomain)
        }
    }
    with_cx(|cx, _| {
        for missing_group in [false, true] {
            let geometry =
                InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 2.0), plane(9, -1.0, 0.0, 2.0)])
                    .unwrap();
            let field = Field { missing_group };
            let media = [SceneMedium::Differentiable {
                field: &field,
                identity: ContentHash([125; 32]),
            }];
            let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
            let native = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![diffuse(1, 0), opaque(9)],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
            let lamp = source(Point3::new(0.6, 0.0, 0.0));
            let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
            let lights = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
            let baseline = lights
                .sample_ray_with_samples(cx, ray(0.0), 550.0, 0, settings(), |_| [0.0; 3])
                .unwrap();
            let timed = TransientCurvedCameraScene::try_new(cx, &camera, &[]).unwrap();
            let unbound =
                TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &[]).unwrap();
            assert_eq!(
                unbound
                    .sample_ray_with_samples_at_time(
                        cx,
                        ray(0.0),
                        550.0,
                        0,
                        Dd::ZERO,
                        settings(),
                        |_| [0.0; 3]
                    )
                    .unwrap()
                    .point_light_sample(),
                &baseline
            );
            let failure = Failure(AtomicUsize::new(0));
            let histories = [history_binding(&lamp, &failure)];
            let observer =
                TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &histories)
                    .unwrap();
            let result = observer.sample_ray_with_samples_at_time(
                cx,
                ray(0.0),
                550.0,
                0,
                Dd::ZERO,
                settings(),
                |_| [0.0; 3],
            );
            let expected = if missing_group {
                TransientPointLightCameraError::MissingGroupDelay {
                    source_identity: ID,
                }
            } else {
                TransientPointLightCameraError::History {
                    source_identity: ID,
                    error: SurfaceHistoryError::OutsideTimeDomain,
                }
            };
            assert_eq!(result, Err(expected));
            assert_eq!(
                failure.0.load(Ordering::Relaxed),
                usize::from(!missing_group)
            );
            if !missing_group {
                let window = ExposureWindow::try_new(Dd::ZERO, Dd::from_f64(1e-9)).unwrap();
                assert!(matches!(
                    observer.sample_ray_with_samples_during_exposure(
                        cx,
                        ray(0.0),
                        550.0,
                        0,
                        window,
                        settings(),
                        |_| [0.0; 3]
                    ),
                    Err(TransientPointLightCameraError::History {
                        source_identity: ID,
                        error: SurfaceHistoryError::ExposureUnsupported,
                    })
                ));
                assert_eq!(
                    failure.0.load(Ordering::Relaxed),
                    1,
                    "shutter must not call point samples"
                );
            }
        }
    });
}
