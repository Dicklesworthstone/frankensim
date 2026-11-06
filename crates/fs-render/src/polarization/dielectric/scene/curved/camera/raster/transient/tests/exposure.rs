//! G1/G3: original pupil/QE paths, shutter energy and chromatic pulse area.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct Integral<F> {
    evaluate: F,
    calls: AtomicUsize,
}
impl<F> Integral<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            calls: AtomicUsize::new(0),
        }
    }
}
impl<F: Fn(SurfaceExposureQuery, usize) -> Result<f64, SurfaceHistoryError> + Sync>
    SurfaceEmissionHistory for Integral<F>
{
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("shutter must use the source integral")
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        (self.evaluate)(query, self.calls.fetch_add(1, Ordering::Relaxed))
    }
}

#[test]
fn g3_unit_shutter_preserves_native_color_and_sensor_with_four_shared_paths() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let unit = Integral::new(|_, _| Ok(1.0));
        let transient = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&unit)]).unwrap();
        let optical_camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 3, 2).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let window = ExposureWindow::try_new(epoch, epoch + Dd::from_f64(1e-8)).unwrap();
        let original = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 9)
            .unwrap();
        let exposed = raster
            .sample_pixel_during_exposure(cx, &transient, 0, window, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(*exposed.mean_pixel_sample(), original);
        assert_eq!(unit.calls.load(Ordering::Relaxed), 4);
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.2, 0.3, 0.0]),
            (680.0, [0.2, 0.3, 0.0]),
        ])
        .unwrap();
        let static_sensor = raster
            .sample_pixel_sensor(cx, &scene, &sensor, 0, settings(), 73, 4, 9)
            .unwrap();
        let exposed = raster
            .sample_pixel_sensor_during_exposure(
                cx,
                &transient,
                &sensor,
                0,
                window,
                settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(*exposed.mean_pixel_sample(), static_sensor);
        assert_eq!(unit.calls.load(Ordering::Relaxed), 8);
        for channel in 0..3 {
            near(
                exposed.qe_weighted_energy_j_m2_sr()[channel],
                static_sensor.qe_weighted_radiance_w_m2_sr()[channel] * 1e-8,
                1e-20,
            );
            near(
                exposed.photoelectrons_per_m2_sr()[channel],
                static_sensor.photoelectron_rate_per_m2_sr()[channel] * 1e-8,
                0.01,
            );
        }
        let one = QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let first = raster
            .sample_pixel_sensor_during_exposure(
                cx,
                &transient,
                &one,
                0,
                window,
                settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(unit.calls.load(Ordering::Relaxed), 12);
        assert_eq!(
            first.mean_pixel_sample().spectral_samples(),
            exposed.mean_pixel_sample().spectral_samples()
        );
        assert_eq!(exposed.photoelectrons_per_m2_sr()[2], 0.0);
        assert_eq!(exposed.qe_weighted_energy_j_m2_sr()[2], 0.0);
    });
}

// Independent continuous triangle primitive, with peak one at half_width.
fn triangle_area(time: f64, half_width: f64) -> f64 {
    if time <= 0.0 {
        0.0
    } else if time <= half_width {
        time * time / (2.0 * half_width)
    } else if time < 2.0 * half_width {
        half_width - (2.0 * half_width - time).powi(2) / (2.0 * half_width)
    } else {
        half_width
    }
}

#[test]
fn g1_dispersion_shifts_pulse_area_before_native_qe_and_photon_conversion() {
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
        let knots = [
            (-10e-9, 0.0),
            (0.0, 0.0),
            (0.5e-9, 1.0),
            (1e-9, 0.0),
            (10e-9, 0.0),
        ];
        let open = 10.6e-9;
        let close = 11.2e-9;
        let mut settings = settings();
        settings.event.surface_tolerance_m = 1e-9;
        for epoch in [Dd::ZERO, Dd::from_f64(2.0_f64.powi(30))] {
            let pulse = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
            let transient =
                TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&pulse)]).unwrap();
            let window =
                ExposureWindow::try_new(epoch + Dd::from_f64(open), epoch + Dd::from_f64(close))
                    .unwrap();
            let observed = raster
                .sample_pixel_sensor_during_exposure(
                    cx, &transient, &sensor, 0, window, settings, 79, 4, 6,
                )
                .unwrap();
            let pixel = observed.mean_pixel_sample();
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            let mut phase_energy = 0.0;
            for lane in 0..4 {
                let wavelength = pixel.detector_wavelengths()[lane].wavelength_nm();
                let ray = pixel.camera_rays()[lane];
                let length = (2.0 - ray.origin.z) / ray.dir.z;
                let group = 1.4 + 0.06 / (wavelength * 1e-3).powi(2);
                let phase = 1.4 + 0.02 / (wavelength * 1e-3).powi(2);
                let area_for = |index: f64| {
                    let delay = index * length / 299_792_458.0;
                    triangle_area(close - delay, 0.5e-9) - triangle_area(open - delay, 0.5e-9)
                };
                let fluence = 3.0 * (-0.12 * length).exp() * area_for(group);
                near(
                    pixel.spectral_samples()[lane].radiance_w_m2_sr_nm(),
                    fluence / (close - open),
                    2e-7,
                );
                let qe = [
                    (wavelength - 420.0) / 260.0,
                    (680.0 - wavelength) / 260.0,
                    0.0,
                ];
                for channel in 0..3 {
                    let weighted = fluence * qe[channel] * 260.0 / 4.0;
                    energy[channel] += weighted;
                    electrons[channel] +=
                        weighted * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
                }
                phase_energy += 3.0 * (-0.12 * length).exp() * area_for(phase) * 260.0 / 4.0;
            }
            for channel in 0..3 {
                near(
                    observed.qe_weighted_energy_j_m2_sr()[channel],
                    energy[channel],
                    3e-14,
                );
                near(
                    observed.photoelectrons_per_m2_sr()[channel],
                    electrons[channel],
                    2e-6 * electrons[channel].max(1.0),
                );
            }
            assert!((phase_energy - energy[0] - energy[1]).abs() > 1e-12);
        }
    });
}

#[test]
fn g4_a_failed_interval_refuses_the_complete_detector_without_time_retry() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = Integral::new(|_, call| {
            if call == 1 {
                Err(SurfaceHistoryError::OutsideTimeDomain)
            } else {
                Ok(1.0)
            }
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)]).unwrap();
        let optical_camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(420.0, [0.5, 0.5]), (680.0, [0.5, 0.5])])
                .unwrap();
        let window = ExposureWindow::try_new(Dd::ZERO, Dd::from_f64(1e-8)).unwrap();
        assert!(matches!(
            raster.sample_pixel_sensor_during_exposure(
                cx,
                &transient,
                &sensor,
                0,
                window,
                settings(),
                17,
                0,
                0
            ),
            Err(TransientRasterError::Source(
                TransientCameraError::History {
                    object_id: 7,
                    error: SurfaceHistoryError::OutsideTimeDomain
                }
            ))
        ));
        assert_eq!(history.calls.load(Ordering::Relaxed), 2);
    });
}
