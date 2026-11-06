use super::*;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::analyzer::CameraAnalyzer;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};
use fs_geom::Vec3;

struct CountedPulse<'a> {
    source: PiecewiseLinearHistory<'a>,
    points: AtomicUsize,
    intervals: AtomicUsize,
}
impl SurfaceEmissionHistory for CountedPulse<'_> {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        self.points.fetch_add(1, Ordering::Relaxed);
        self.source.sample(cx, query)
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        self.intervals.fetch_add(1, Ordering::Relaxed);
        self.source.average(cx, query)
    }
}

#[test]
fn g1_analyzed_chromatic_pulses_reach_native_xyz_and_qe_shutters_in_four_paths() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let optical_camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 3, 2).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.0, 1.0, 0.0]),
            (680.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let one = QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [
            (-10e-9, 0.0),
            (0.0, 0.0),
            (0.5e-9, 1.0),
            (1e-9, 0.0),
            (10e-9, 0.0),
        ];
        let pulse = CountedPulse {
            source: PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap(),
            points: AtomicUsize::new(0),
            intervals: AtomicUsize::new(0),
        };
        let scalar = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&pulse)]).unwrap();
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(10.6e-9), epoch + Dd::from_f64(11.2e-9))
                .unwrap();
        let time = epoch + Dd::from_f64(10.9e-9);
        let original_xyz = raster
            .sample_pixel_during_exposure(cx, &scalar, 0, window, settings(), 79, 4, 6)
            .unwrap();
        let original_sensor = raster
            .sample_pixel_sensor_during_exposure(
                cx,
                &scalar,
                &sensor,
                0,
                window,
                settings(),
                79,
                4,
                6,
            )
            .unwrap();
        let original_point = raster
            .sample_pixel_sensor_transient_at_time(
                cx,
                &scalar,
                &sensor,
                0,
                time,
                settings(),
                79,
                4,
                6,
            )
            .unwrap();
        let mut xyz = [0.0; 3];
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for sign in [-1.0, 1.0] {
            let analyzer = CameraAnalyzer::try_new(
                Vec3::new(1.0, 0.0, 0.0),
                Stokes {
                    i: 0.5,
                    q: sign * 0.2,
                    u: sign * 0.15,
                    v: sign * 0.1,
                },
            )
            .unwrap();
            let selected = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&pulse)])
                .unwrap()
                .with_analyzer(analyzer);
            let before = pulse.points.load(Ordering::Relaxed);
            let point = raster
                .sample_pixel_sensor_transient_at_time(
                    cx,
                    &selected,
                    &sensor,
                    0,
                    time,
                    settings(),
                    79,
                    4,
                    6,
                )
                .unwrap();
            assert_eq!(pulse.points.load(Ordering::Relaxed) - before, 4);
            let point_calls = pulse.points.load(Ordering::Relaxed);
            let before = pulse.intervals.load(Ordering::Relaxed);
            let color = raster
                .sample_pixel_during_exposure(cx, &selected, 0, window, settings(), 79, 4, 6)
                .unwrap();
            assert_eq!(pulse.intervals.load(Ordering::Relaxed) - before, 4);
            let exposed = raster
                .sample_pixel_sensor_during_exposure(
                    cx,
                    &selected,
                    &sensor,
                    0,
                    window,
                    settings(),
                    79,
                    4,
                    6,
                )
                .unwrap();
            assert_eq!(pulse.intervals.load(Ordering::Relaxed) - before, 8);
            let single = raster
                .sample_pixel_sensor_during_exposure(
                    cx,
                    &selected,
                    &one,
                    0,
                    window,
                    settings(),
                    79,
                    4,
                    6,
                )
                .unwrap();
            assert_eq!(pulse.intervals.load(Ordering::Relaxed) - before, 12);
            assert_eq!(pulse.points.load(Ordering::Relaxed), point_calls);
            assert_eq!(
                single.mean_pixel_sample().spectral_samples(),
                exposed.mean_pixel_sample().spectral_samples()
            );
            assert_eq!(
                exposed.mean_pixel_sample().camera_rays(),
                original_sensor.mean_pixel_sample().camera_rays()
            );
            assert_eq!(
                exposed.mean_pixel_sample().detector_wavelengths(),
                original_sensor.mean_pixel_sample().detector_wavelengths()
            );
            for channel in 0..3 {
                near(
                    point.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                    0.5 * original_point.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                    original_point.pixel_sample().photoelectron_rate_per_m2_sr()[channel]
                        .abs()
                        .max(1.0)
                        * 1e-13,
                );
                xyz[channel] += color.mean_pixel_sample().xyz()[channel];
                energy[channel] += exposed.qe_weighted_energy_j_m2_sr()[channel];
                electrons[channel] += exposed.photoelectrons_per_m2_sr()[channel];
            }
            assert_eq!(exposed.photoelectrons_per_m2_sr()[2], 0.0);
        }
        assert!(energy[0] + energy[1] > 0.0);
        for channel in 0..3 {
            near(
                xyz[channel],
                original_xyz.mean_pixel_sample().xyz()[channel],
                1e-12,
            );
            near(
                energy[channel],
                original_sensor.qe_weighted_energy_j_m2_sr()[channel],
                1e-20,
            );
            near(
                electrons[channel],
                original_sensor.photoelectrons_per_m2_sr()[channel],
                1e-2,
            );
        }
        let unsupported = History::new(|_, _| Ok(1.0));
        let selected = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&unsupported)])
            .unwrap()
            .with_analyzer(CameraAnalyzer::try_circular(Vec3::new(1.0, 0.0, 0.0), true).unwrap());
        assert_eq!(
            raster.sample_pixel_sensor_during_exposure(
                cx,
                &selected,
                &sensor,
                0,
                window,
                settings(),
                79,
                4,
                6
            ),
            Err(TransientRasterError::Source(
                TransientCameraError::History {
                    object_id: 7,
                    error: SurfaceHistoryError::ExposureUnsupported,
                }
            )),
        );
    });
}
