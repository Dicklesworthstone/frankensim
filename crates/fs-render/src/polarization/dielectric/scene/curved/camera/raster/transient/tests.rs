//! G1/G3: actual chromatic source times through native pixel and QE integrals.

mod exposure;

use super::super::tests::{bindings, camera, geometry, settings, with_cx};
use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use fs_blake3::ContentHash;
use std::sync::atomic::{AtomicUsize, Ordering};

struct History<F> {
    evaluate: F,
    calls: AtomicUsize,
}

impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            calls: AtomicUsize::new(0),
        }
    }
}

impl<F> SurfaceEmissionHistory for History<F>
where
    F: Fn(SurfaceEmissionQuery, usize) -> Result<f64, SurfaceHistoryError> + Sync,
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        (self.evaluate)(query, self.calls.fetch_add(1, Ordering::Relaxed))
    }
}

fn binding(history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: 7,
        source_identity: ContentHash([87; 32]),
        history,
    }
}

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.12).unwrap(),
        GlassProvenance::Custom,
    )
}

fn emitter() -> SurfaceEmitter {
    SurfaceEmitter::new(
        7,
        RadianceSpectrum::try_constant(3.0).unwrap(),
        EmissionSides::Both,
    )
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn g3_steady_histories_preserve_pixels_channels_and_absolute_sampling_identity() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let unit = History::new(|_, _| Ok(1.0));
        let transient = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&unit)]).unwrap();
        let optical_camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 3, 2).unwrap();
        let clock = Dd::from_f64(2.0_f64.powi(30)) + Dd::from_f64(1e-8);
        let original = raster
            .sample_pixel(cx, &scene, 0, settings(), 73, 4, 9)
            .unwrap();
        let timed = raster
            .sample_pixel_transient_at_time(cx, &transient, 0, clock, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(*timed.pixel_sample(), original);
        assert_eq!(timed.observation_time_s(), clock);
        assert_eq!(unit.calls.load(Ordering::Relaxed), 4);
        raster
            .sample_pixel_transient_at_time(cx, &transient, 0, clock, settings(), 73, 1, 3)
            .unwrap();
        assert_eq!(
            raster
                .sample_pixel_transient_at_time(cx, &transient, 0, clock, settings(), 73, 4, 9)
                .unwrap(),
            timed
        );

        let one = QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let three = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.2, 0.3, 0.0]),
            (680.0, [0.2, 0.3, 0.0]),
        ])
        .unwrap();
        let before = unit.calls.load(Ordering::Relaxed);
        let first = raster
            .sample_pixel_sensor_transient_at_time(
                cx,
                &transient,
                &one,
                0,
                clock,
                settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(unit.calls.load(Ordering::Relaxed) - before, 4);
        let before = unit.calls.load(Ordering::Relaxed);
        let second = raster
            .sample_pixel_sensor_transient_at_time(
                cx,
                &transient,
                &three,
                0,
                clock,
                settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(unit.calls.load(Ordering::Relaxed) - before, 4);
        let static_sensor = raster
            .sample_pixel_sensor(cx, &scene, &three, 0, settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(*second.pixel_sample(), static_sensor);
        assert_eq!(
            first.pixel_sample().camera_rays(),
            second.pixel_sample().camera_rays()
        );
        assert_eq!(
            first.pixel_sample().spectral_samples(),
            second.pixel_sample().spectral_samples()
        );
        near(
            first.pixel_sample().qe_weighted_radiance_w_m2_sr()[0],
            second.pixel_sample().qe_weighted_radiance_w_m2_sr()[0]
                + second.pixel_sample().qe_weighted_radiance_w_m2_sr()[1],
            1e-12,
        );
        assert_eq!(second.pixel_sample().photoelectron_rate_per_m2_sr()[2], 0.0);
    });
}

#[test]
fn g1_chromatic_arrival_drives_native_qe_energy_and_photon_rates_at_one_clock() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let pulse = History::new(|query: SurfaceEmissionQuery, _| {
            Ok((-(query.emission_time_s.to_f64() / 0.6e-9).powi(2)).exp())
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&pulse)]).unwrap();
        let optical_camera = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 3, 2).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.0, 1.0, 0.0]),
            (680.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let clock = 11.0e-9;
        let mut settings = settings();
        settings.event.surface_tolerance_m = 1e-9;
        let sampled = raster
            .sample_pixel_sensor_transient_at_time(
                cx,
                &transient,
                &sensor,
                0,
                Dd::from_f64(clock),
                settings,
                79,
                4,
                6,
            )
            .unwrap();
        let pixel = sampled.pixel_sample();
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        let mut incorrect_phase_energy = 0.0;
        for lane in 0..4 {
            let wavelength = pixel.detector_wavelengths()[lane].wavelength_nm();
            let ray = pixel.camera_rays()[lane];
            let length = (2.0 - ray.origin.z) / ray.dir.z;
            let phase = 1.4 + 0.02 / (wavelength * 1e-3).powi(2);
            let group = 1.4 + 0.06 / (wavelength * 1e-3).powi(2);
            let source_time = clock - group * length / 299_792_458.0;
            let radiance = 3.0 * (-0.12 * length - (source_time / 0.6e-9).powi(2)).exp();
            near(
                pixel.spectral_samples()[lane].radiance_w_m2_sr_nm(),
                radiance,
                2e-7,
            );
            let eta = [
                (wavelength - 420.0) / 260.0,
                (680.0 - wavelength) / 260.0,
                0.0,
            ];
            // Sum of the two exact linear QEs is one: proposal 1/260nm.
            for channel in 0..3 {
                let weighted = radiance * eta[channel] * 260.0 / 4.0;
                energy[channel] += weighted;
                electrons[channel] +=
                    weighted * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
            let wrong_time = clock - phase * length / 299_792_458.0;
            incorrect_phase_energy +=
                3.0 * (-0.12 * length - (wrong_time / 0.6e-9).powi(2)).exp() * 260.0 / 4.0;
        }
        for channel in 0..3 {
            near(
                pixel.qe_weighted_radiance_w_m2_sr()[channel],
                energy[channel],
                5e-5,
            );
            near(
                pixel.photoelectron_rate_per_m2_sr()[channel],
                electrons[channel],
                1e-7 * electrons[channel].max(1.0),
            );
        }
        assert_eq!(pulse.calls.load(Ordering::Relaxed), 4);
        assert!((incorrect_phase_energy - energy[0] - energy[1]).abs() > 1.0);
        assert_eq!(pixel.qe_weighted_radiance_w_m2_sr()[2], 0.0);
        assert_eq!(pixel.photoelectron_rate_per_m2_sr()[2], 0.0);
    });
}

#[test]
fn g4_a_failed_history_lane_aborts_the_original_complete_pixel_without_retry() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let material = glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [emitter()];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let interrupted = History::new(|_, call| {
            if call == 1 {
                Err(SurfaceHistoryError::OutsideTimeDomain)
            } else {
                Ok(1.0)
            }
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&interrupted)]).unwrap();
        let optical_camera = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&optical_camera, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(420.0, [1.0]), (680.0, [1.0])]).unwrap();
        let result = raster.sample_pixel_sensor_transient_at_time(
            cx,
            &transient,
            &sensor,
            0,
            Dd::ZERO,
            settings(),
            79,
            0,
            6,
        );
        assert!(matches!(
            result,
            Err(TransientRasterError::Source(
                TransientCameraError::History {
                    object_id: 7,
                    error: SurfaceHistoryError::OutsideTimeDomain,
                }
            ))
        ));
        assert_eq!(interrupted.calls.load(Ordering::Relaxed), 2);
        assert!(matches!(
            raster.sample_pixel_transient_at_time(
                cx,
                &transient,
                0,
                Dd::ZERO,
                settings(),
                79,
                1,
                6
            ),
            Err(TransientRasterError::Raster(
                CurvedRasterError::InvalidRaster
            ))
        ));
        assert_eq!(interrupted.calls.load(Ordering::Relaxed), 2);
        let cancel = History::new(|_, _| {
            gate.request();
            Ok(1.0)
        });
        let cancelled =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&cancel)]).unwrap();
        assert!(
            raster
                .sample_pixel_transient_at_time(cx, &cancelled, 0, Dd::ZERO, settings(), 79, 0, 6)
                .is_err()
        );
        assert_eq!(cancel.calls.load(Ordering::Relaxed), 1);
    });
}
