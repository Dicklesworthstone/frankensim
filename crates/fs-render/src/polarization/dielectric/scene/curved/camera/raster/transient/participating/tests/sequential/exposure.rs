use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct AverageOnly<'a> {
    history: PiecewiseLinearHistory<'a>,
    queries: Mutex<Vec<SurfaceExposureQuery>>,
}

impl AverageOnly<'_> {
    fn binding(&self) -> SurfaceEmissionHistoryBinding<'_> {
        SurfaceEmissionHistoryBinding {
            object_id: 7,
            source_identity: ContentHash([97; 32]),
            history: self,
        }
    }
}

impl SurfaceEmissionHistory for AverageOnly<'_> {
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("shutter integration must not take temporal samples")
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
fn analyzed_molecular_exposure_keeps_clipped_packet_measure_and_original_work() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let vacuum = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.12).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let surface = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let wide = lens_camera(visible_glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &surface, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let mut radii = packet(&wide, cx, full.wavelengths_nm()).map(|beam| {
            let ray = beam.unwrap().ray();
            let p = ray
                .origin
                .offset(ray.dir.scale((0.12 - ray.origin.z) / ray.dir.z));
            p.x.hypot(p.y)
        });
        radii.sort_by(f64::total_cmp);
        assert!(radii[2] - radii[1] > 1e-8);
        let physical = lens_camera(visible_glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [(0.0, 1.0), (1.0, 1.0)];
        let history = AverageOnly {
            history: PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap(),
            queries: Mutex::new(Vec::new()),
        };
        let transient = TransientCurvedCameraScene::try_new(cx, &surface, &[history.binding()])
            .unwrap()
            .with_analyzer(analyzer);
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(1e-6), epoch + Dd::from_f64(2e-6))
                .unwrap();
        let phase = PhaseFunction::Rayleigh;
        for coefficient in [0.0, 0.5] {
            let sigma = MassScatteringSpectrum::constant(coefficient).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let steady = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
            let scene =
                TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
            let expected = raster
                .sample_pixel_sequential_rayleigh_analyzed(
                    cx,
                    &steady,
                    0,
                    analyzer,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let before = history.queries.lock().unwrap().len();
            let actual = raster
                .sample_pixel_sequential_rayleigh_during_exposure(
                    cx,
                    &scene,
                    0,
                    window,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(actual.window(), window);
            assert_eq!(actual.mean_pixel_sample(), expected.pixel_sample());
            assert_eq!(actual.volume_events(), expected.volume_events());
            assert_eq!(actual.free_flight_draws(), expected.free_flight_draws());
            assert_eq!(
                actual
                    .mean_pixel_sample()
                    .camera_rays()
                    .iter()
                    .flatten()
                    .count(),
                2
            );
            let queries = history.queries.lock().unwrap()[before..].to_vec();
            assert_eq!(queries.len(), 2);
            let beams = packet(&physical, cx, actual.mean_pixel_sample().wavelengths_nm());
            let mut source = 0;
            for lane in 0..CAMERA_WAVELENGTHS {
                if let Some(beam) = beams[lane] {
                    let path = actual.mean_pixel_sample().spectral_samples()[lane].unwrap();
                    let delay = beam.group_delay_s().unwrap()
                        + Dd::from_f64(path.path().group_optical_path_length_m().unwrap())
                            / Dd::from_f64(C0);
                    assert_eq!(queries[source].opening.observation_time_s, window.open_s());
                    assert_eq!(queries[source].observation_close_s, window.close_s());
                    assert_eq!(queries[source].opening.group_delay_s, delay);
                    assert_eq!(
                        queries[source].opening.emission_time_s,
                        window.open_s() - delay
                    );
                    assert_eq!(queries[source].emission_close_s, window.close_s() - delay);
                    if coefficient > 0.0 {
                        assert!(actual.volume_events()[lane] > 0);
                    }
                    source += 1;
                } else {
                    assert_eq!(actual.volume_events()[lane], 0);
                    assert_eq!(actual.free_flight_draws()[lane], 0);
                }
            }
        }
    });
}

#[test]
fn infrared_exposure_shares_refracted_molecular_paths_and_integrates_qe_channels_once() {
    with_cx(|cx, _| {
        let ior = |base: f64, slope: f64| {
            TabulatedIor::try_new_in_band(
                core::array::from_fn(|i| {
                    DielectricIorSample::try_new_in_band(
                        8000.0 + 750.0 * i as f64,
                        base + slope * i as f64,
                    )
                    .unwrap()
                }),
                750.0,
            )
            .unwrap()
        };
        let material = DielectricGlass::new(
            ior(1.5625, -1.0 / 128.0),
            BeerLambertAbsorption::try_constant(2.0).unwrap(),
            GlassProvenance::Custom,
        );
        let vacuum = DielectricGlass::new(
            ior(1.0, 0.0),
            BeerLambertAbsorption::try_constant(0.12).unwrap(),
            GlassProvenance::Custom,
        );
        let geometry = geometry();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_new_in_band(vec![(8000.0, 3.0), (14000.0, 3.0)]).unwrap(),
            EmissionSides::Both,
        )];
        let surface = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::source_table(
            &[(8000.0, 0.5), (14000.0, 0.5)],
            ContentHash([92; 32]),
        )
        .unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let physical = lens_camera(material, 11000.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [(0.0, 0.0), (1.0, 1e6)];
        let history = AverageOnly {
            history: PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap(),
            queries: Mutex::new(Vec::new()),
        };
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[history.binding()]).unwrap();
        let scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(1e-6), epoch + Dd::from_f64(2e-6))
                .unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let base = raster
            .sample_pixel_sensor_sequential_rayleigh(
                cx,
                &steady,
                &sensor,
                0,
                trace_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sensor_sequential_rayleigh_during_exposure(
                cx,
                &scene,
                &sensor,
                0,
                window,
                trace_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(actual.window(), window);
        assert_eq!(actual.volume_events(), base.volume_events());
        assert_eq!(actual.free_flight_draws(), base.free_flight_draws());
        assert!(actual.volume_events().iter().all(|&n| n > 0));
        let queries = history.queries.lock().unwrap();
        assert_eq!(queries.len(), 4);
        let output = actual.mean_pixel_sample();
        let wavelengths = output
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(&physical, cx, &wavelengths);
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for lane in 0..CAMERA_WAVELENGTHS {
            let beam = beams[lane].unwrap();
            let reference = base.pixel_sample().spectral_samples()[lane].unwrap();
            assert_eq!(
                output.spectral_samples()[lane].unwrap().path(),
                reference.path()
            );
            let delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(reference.path().group_optical_path_length_m().unwrap())
                    / Dd::from_f64(C0);
            assert_eq!(queries[lane].opening.group_delay_s, delay);
            assert_eq!(queries[lane].opening.observation_time_s, window.open_s());
            assert_eq!(queries[lane].observation_close_s, window.close_s());
            let mean = (((window.open_s() - epoch) - delay).to_f64()
                + ((window.close_s() - epoch) - delay).to_f64())
                * 0.5e6;
            let radiance = reference.radiance_w_m2_sr_nm() * mean;
            near(
                output.spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                radiance,
            );
            let nm = wavelengths[lane];
            let qe = [(nm - 8000.0) / 6000.0, (14000.0 - nm) / 6000.0, 0.0];
            for c in 0..3 {
                let contribution = radiance * beam.geometry_weight() * qe[c] * 6000.0 / 4.0;
                energy[c] += contribution;
                electrons[c] += contribution * nm * 1e-9 / (6.62607015e-34 * C0);
            }
        }
        for c in 0..3 {
            near(output.qe_weighted_radiance_w_m2_sr()[c], energy[c]);
            near(output.photoelectron_rate_per_m2_sr()[c], electrons[c]);
            near(
                actual.qe_weighted_energy_j_m2_sr()[c],
                energy[c] * window.duration_s().to_f64(),
            );
            near(
                actual.photoelectrons_per_m2_sr()[c],
                electrons[c] * window.duration_s().to_f64(),
            );
        }
        assert_eq!(actual.photoelectrons_per_m2_sr()[2], 0.0);
        drop(queries);
        let refusal = History::new(|_, _| panic!("no temporal point samples"));
        let failed =
            TransientCurvedCameraScene::try_new(cx, &surface, &[refusal.binding()]).unwrap();
        let failed_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &failed, &scattering).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_rayleigh_during_exposure(
                    cx,
                    &failed_scene,
                    &sensor,
                    0,
                    window,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::ExposureUnsupported,
            }),
        );
        assert_eq!(refusal.calls.load(Ordering::Relaxed), 0);
    });
}
