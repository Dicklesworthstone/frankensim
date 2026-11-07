use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct AverageOnly<'a> {
    history: PiecewiseLinearHistory<'a>,
    queries: Mutex<Vec<SurfaceExposureQuery>>,
}

impl SurfaceEmissionHistory for AverageOnly<'_> {
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("source-owned interval integration must not draw times")
    }

    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        self.queries.lock().unwrap().push(query);
        self.history.average(cx, query)
    }
}

#[test]
fn real_lens_exposure_retards_the_original_source_interval_after_two_molecular_turns() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let physical = lens(glass);
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [(0.0, 0.0), (1e-6, 10.0)];
        let history = AverageOnly {
            history: PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap(),
            queries: Mutex::new(Vec::new()),
        };
        let bindings = [SurfaceEmissionHistoryBinding {
            object_id: 9,
            source_identity: ContentHash([96; 32]),
            history: &history,
        }];
        let transient = TransientCurvedCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(5e-8), epoch + Dd::from_f64(9e-8))
                .unwrap();
        for wavelength in [450.0, 650.0] {
            let optical = physical
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[wavelength],
                )
                .unwrap()[0]
                .unwrap();
            let base = steady
                .sample_sequential_ray_with_samples(
                    cx,
                    &optical,
                    0,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            let mut draws = Vec::new();
            let sample = scene
                .sample_sequential_ray_with_samples_during_exposure(
                    cx,
                    &optical,
                    0,
                    window,
                    settings(),
                    |_| panic!("opaque"),
                    |ordinal| {
                        draws.push(ordinal);
                        two_turns(ordinal)
                    },
                )
                .unwrap();
            assert_eq!(draws, [0, 1, 2]);
            assert_eq!(sample.volume_events(), 2);
            assert_eq!(sample.free_flight_draws(), 3);
            assert_eq!(
                sample.mean_radiance_sample().path(),
                base.radiance_sample().path()
            );
            let delay = optical.group_delay_s().unwrap()
                + Dd::from_f64(
                    base.radiance_sample()
                        .path()
                        .group_optical_path_length_m()
                        .unwrap(),
                ) / Dd::from_f64(LIGHT_M_S);
            let mean = (((window.open_s() - epoch) - delay).to_f64()
                + ((window.close_s() - epoch) - delay).to_f64())
                * 0.5e7;
            let expected = base.radiance_sample().radiance_w_m2_sr_nm() * mean;
            near(
                sample.mean_radiance_sample().radiance_w_m2_sr_nm(),
                expected,
                1e-11,
            );
            near(
                sample.integrated_radiance_j_m2_sr_nm(),
                expected * window.duration_s().to_f64(),
                1e-18,
            );
            let query = *history.queries.lock().unwrap().last().unwrap();
            assert_eq!(query.opening.observation_time_s, window.open_s());
            assert_eq!(query.observation_close_s, window.close_s());
            assert_eq!(query.opening.group_delay_s, delay);
            assert_eq!(query.opening.emission_time_s, window.open_s() - delay);
            assert_eq!(query.emission_close_s, window.close_s() - delay);
            let tracking = scene
                .sample_sequential_ray_with_tracking_samples_during_exposure(
                    cx,
                    &optical,
                    0,
                    window,
                    settings(),
                    |_| panic!("opaque"),
                    |ordinal| {
                        let [a, b, c] = two_turns(ordinal);
                        [a, b, c, 0.75]
                    },
                )
                .unwrap();
            assert_eq!(tracking, sample);
        }
        assert_eq!(history.queries.lock().unwrap().len(), 4);
        let optical = physical
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        let ordinary = physical
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(
            scene
                .sample_sequential_ray_during_exposure(cx, &ordinary, 0, window, settings(), 31, 42)
                .unwrap_err(),
            TransientCameraError::MissingLensGroupDelay,
        );
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let analyzed = TransientCurvedCameraScene::try_new(cx, &camera, &bindings)
            .unwrap()
            .with_analyzer(analyzer);
        let analyzed_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &analyzed, &scattering).unwrap();
        assert_eq!(
            analyzed_scene
                .sample_sequential_ray_during_exposure(cx, &optical, 0, window, settings(), 31, 42)
                .unwrap_err(),
            TransientCameraError::AnalyzerMustPrecedeLens,
        );
        gate.request();
        assert!(matches!(
            scene.sample_sequential_ray_during_exposure(
                cx,
                &optical,
                0,
                window,
                settings(),
                31,
                42
            ),
            Err(TransientCameraError::Camera(CameraRadianceError::Scene(
                CurvedSceneError::Scene(SceneError::Cancelled)
            )))
        ));
    });
}
