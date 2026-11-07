//! G1/G3/G4: sensor shutter means retain optics, clipping and original QE measure.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery, piecewise_linear::PiecewiseLinearHistory,
};

struct Recording<'a> {
    source: &'a dyn SurfaceEmissionHistory,
    queries: Mutex<Vec<SurfaceExposureQuery>>,
}
impl<'a> Recording<'a> {
    fn new(source: &'a dyn SurfaceEmissionHistory) -> Self {
        Self {
            source,
            queries: Mutex::new(Vec::new()),
        }
    }
}
impl SurfaceEmissionHistory for Recording<'_> {
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("shutter integration must never point-sample the source history")
    }
    fn average(&self, cx: &Cx<'_>, q: SurfaceExposureQuery) -> Result<f64, SurfaceHistoryError> {
        self.queries.lock().unwrap().push(q);
        self.source.average(cx, q)
    }
}

#[test]
fn finite_sensor_shutter_keeps_analyzed_clipped_packet_and_original_clock() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let vacuum = vacuum();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [(0.0, 1.0), (1.0, 1.0)];
        let constant = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
        let history = Recording::new(&constant);
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let transient = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)])
            .unwrap()
            .with_analyzer(analyzer);
        let wide = camera(glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let mut radii = packet(&wide, cx, full.wavelengths_nm(), None).map(|beam| {
            let ray = beam.unwrap().ray();
            let p = ray
                .origin
                .offset(ray.dir.scale((0.12 - ray.origin.z) / ray.dir.z));
            p.x.hypot(p.y)
        });
        radii.sort_by(f64::total_cmp);
        assert!(radii[2] - radii[1] > 1e-8);
        let physical = camera(glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(1e-8), epoch + Dd::from_f64(2e-8))
                .unwrap();
        let actual = raster
            .sample_pixel_sequential_during_exposure(
                cx,
                &transient,
                0,
                window,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let steady = raster
            .sample_pixel_sequential_analyzed(cx, &scene, 0, analyzer, settings(), SEED, 0, SAMPLE)
            .unwrap();
        assert_eq!(*actual.mean_pixel_sample(), steady);
        assert_eq!(actual.window(), window);
        assert_eq!(steady.camera_rays().iter().flatten().count(), 2);
        let queries = history.queries.lock().unwrap().clone();
        assert_eq!(queries.len(), 2);
        let beams = packet(&physical, cx, steady.wavelengths_nm(), Some(analyzer));
        for (lane, beam) in beams.iter().enumerate() {
            let Some(beam) = beam else {
                continue;
            };
            let q = queries
                .iter()
                .find(|q| q.opening.wavelength_nm == beam.wavelength_nm())
                .unwrap();
            let group = steady.spectral_samples()[lane]
                .unwrap()
                .path()
                .group_optical_path_length_m()
                .unwrap();
            let delay = beam.group_delay_s().unwrap() + Dd::from_f64(group) / Dd::from_f64(C);
            assert_eq!(q.opening.observation_time_s, window.open_s());
            assert_eq!(q.observation_close_s, window.close_s());
            assert_eq!(q.duration_s, window.duration_s());
            assert_eq!(q.opening.group_delay_s, delay);
            assert_eq!(q.opening.emission_time_s, window.open_s() - delay);
            assert_eq!(q.emission_close_s, window.close_s() - delay);
        }
        assert_eq!(
            actual,
            raster
                .sample_pixel_sequential_during_exposure(
                    cx,
                    &transient,
                    0,
                    window,
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap()
        );
        let closed_camera = camera(glass(), 550.0, Some(radii[0] * 0.5));
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let before = history.queries.lock().unwrap().len();
        let dark = closed
            .sample_pixel_sequential_during_exposure(
                cx,
                &transient,
                0,
                window,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.mean_pixel_sample().xyz(), [0.0; 3]);
        assert_eq!(dark.mean_pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(dark.window(), window);
        assert_eq!(history.queries.lock().unwrap().len(), before);
        gate.request();
        assert!(
            closed
                .sample_pixel_sequential_during_exposure(
                    cx,
                    &transient,
                    0,
                    window,
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .is_err()
        );
    });
}

#[test]
fn infrared_pulse_exposure_integrates_lens_delay_and_detector_energy_electrons() {
    with_cx(|cx, _| {
        let table = |glass: bool| {
            TabulatedIor::try_new_in_band(
                core::array::from_fn(|i| {
                    DielectricIorSample::try_new_in_band(
                        8000.0 + 750.0 * i as f64,
                        if glass {
                            1.5625 - i as f64 / 128.0
                        } else {
                            1.0
                        },
                    )
                    .unwrap()
                }),
                750.0,
            )
            .unwrap()
        };
        let material = DielectricGlass::new(
            table(true),
            BeerLambertAbsorption::try_constant(2.0).unwrap(),
            GlassProvenance::Custom,
        );
        let vacuum = DielectricGlass::new(
            table(false),
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
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [
            (0.0, 0.0),
            (1e-9, 0.0),
            (2e-9, 4.0),
            (3e-9, 0.0),
            (1e-7, 0.0),
        ];
        let pulse = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
        let history = Recording::new(&pulse);
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)]).unwrap();
        let physical = camera(material, 11000.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(8.2e-9), epoch + Dd::from_f64(1e-8))
                .unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let mut trace = settings();
        trace.event.surface_tolerance_m = 1e-9;
        let exposed = raster
            .sample_pixel_sensor_sequential_during_exposure(
                cx, &transient, &sensor, 0, window, trace, SEED, 0, SAMPLE,
            )
            .unwrap();
        let output = exposed.mean_pixel_sample();
        let wavelengths = output
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(&physical, cx, &wavelengths, None);
        // Independent primitive of the unit-nanosecond triangular pulse.
        let primitive = |t: f64| {
            if t <= 1.0 {
                0.0
            } else if t <= 2.0 {
                2.0 * (t - 1.0).powi(2)
            } else if t <= 3.0 {
                4.0 - 2.0 * (3.0 - t).powi(2)
            } else {
                4.0
            }
        };
        let mut energy_rate = [0.0; 3];
        let mut electron_rate = [0.0; 3];
        for (i, beam) in beams.iter().enumerate() {
            let beam = beam.unwrap();
            let ray = beam.ray();
            let length = (2.0 - ray.origin.z) / ray.dir.z;
            let delay = beam.group_delay_s().unwrap() + Dd::from_f64(length) / Dd::from_f64(C);
            let a = ((window.open_s() - epoch) - delay).to_f64() * 1e9;
            let b = ((window.close_s() - epoch) - delay).to_f64() * 1e9;
            let mean = (primitive(b) - primitive(a)) / (b - a);
            assert!(mean > 0.0);
            let radiance = 3.0 * beam.power().i * (-0.12 * length).exp() * mean;
            near(
                output.spectral_samples()[i].unwrap().radiance_w_m2_sr_nm(),
                radiance,
                2e-7,
            );
            let nm = beam.wavelength_nm();
            let qe = [(nm - 8000.0) / 6000.0, (14000.0 - nm) / 6000.0, 0.0];
            for channel in 0..3 {
                let contribution = radiance * beam.geometry_weight() * qe[channel] * 6000.0 / 4.0;
                energy_rate[channel] += contribution;
                electron_rate[channel] += contribution * nm * 1e-9 / (6.62607015e-34 * C);
            }
        }
        assert_eq!(exposed.window(), window);
        for channel in 0..3 {
            near(
                output.qe_weighted_radiance_w_m2_sr()[channel],
                energy_rate[channel],
                2e-7 * energy_rate[channel].max(1.0),
            );
            near(
                output.photoelectron_rate_per_m2_sr()[channel],
                electron_rate[channel],
                2e-7 * electron_rate[channel].max(1.0),
            );
            let energy = energy_rate[channel] * window.duration_s().to_f64();
            let electrons = electron_rate[channel] * window.duration_s().to_f64();
            near(
                exposed.qe_weighted_energy_j_m2_sr()[channel],
                energy,
                2e-7 * energy.abs().max(1e-20),
            );
            near(
                exposed.photoelectrons_per_m2_sr()[channel],
                electrons,
                2e-7 * electrons.abs().max(1e-20),
            );
        }
        assert_eq!(exposed.photoelectrons_per_m2_sr()[2], 0.0);
        assert_eq!(history.queries.lock().unwrap().len(), 4);
        let outside = ExposureWindow::try_new(epoch + Dd::ONE, epoch + Dd::from_f64(2.0)).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_during_exposure(
                    cx, &transient, &sensor, 0, outside, trace, SEED, 0, SAMPLE,
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::OutsideTimeDomain,
            })
        );
        assert_eq!(history.queries.lock().unwrap().len(), 5);
    });
}
