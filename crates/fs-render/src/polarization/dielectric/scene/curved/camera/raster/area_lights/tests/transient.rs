//! Retarded source connections through the original four-lane detectors.

use super::*;

mod sequential;
use crate::polarization::dielectric::scene::curved::camera::raster::transient::area_lights::TransientAreaLightRasterError;
use crate::polarization::dielectric::scene::curved::camera::transient::area_lights::{
    TransientAreaLightCameraError, TransientCurvedAreaLightCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_blake3::ContentHash;
use fs_math::dd::Dd;
use std::sync::atomic::{AtomicUsize, Ordering};

struct LinearHistory {
    epoch: Dd,
    offset: f64,
    slope: f64,
    calls: AtomicUsize,
    fail_at: AtomicUsize,
}

impl LinearHistory {
    fn evaluate(&self, cx: &Cx<'_>, time: Dd) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        let count = self.calls.fetch_add(1, Ordering::Relaxed) + 1;
        if count == self.fail_at.load(Ordering::Relaxed) {
            return Err(SurfaceHistoryError::OutsideTimeDomain);
        }
        Ok(self.offset + self.slope * (time - self.epoch).to_f64())
    }
}

impl SurfaceEmissionHistory for LinearHistory {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        self.evaluate(cx, query.emission_time_s)
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        // Exact mean of this source's affine law over its own retarded interval.
        self.evaluate(
            cx,
            query.opening.emission_time_s + query.duration_s * Dd::from_f64(0.5),
        )
    }
}

fn history_binding(history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: 9,
        source_identity: ContentHash([95; 32]),
        history,
    }
}

#[test]
fn g3_retarded_area_pixels_keep_native_cie_qe_paths_and_channel_independent_work() {
    fixtures::with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let glass = material();
        let media = [SceneMedium::Homogeneous(&glass)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
                SurfaceBinding {
                    object_id: 10,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 2e-7,
        }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let history = LinearHistory {
            epoch,
            offset: 1.0,
            slope: 0.0,
            calls: AtomicUsize::new(0),
            fail_at: AtomicUsize::new(0),
        };
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &area).unwrap();
        let physical = fixtures::camera(0.01);
        let settings = || scene_settings(&geometry, &physical);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let clock = epoch + Dd::from_f64(30e-9);
        let ordinary = raster
            .sample_pixel_area_lights(
                cx,
                &area,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let result = raster
            .sample_pixel_area_lights_at_time(
                cx,
                &observer,
                0,
                clock,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(result.area_sample(), &ordinary);
        assert_eq!(result.observation_time_s(), clock);
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = raster
            .sample_pixel_area_lights_during_exposure(
                cx,
                &observer,
                0,
                window,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(
            exposed.pixel_sample().mean_pixel_sample(),
            ordinary.pixel_sample()
        );
        assert_eq!(
            exposed.area_radiance_w_m2_sr_nm(),
            ordinary.area_radiance_w_m2_sr_nm()
        );
        assert_eq!(exposed.trace_evaluations(), ordinary.trace_evaluations());
        let one = QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [1.0]), (LAMBDA_MAX, [1.0])])
            .unwrap();
        let three = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let static_three = raster
            .sample_pixel_sensor_area_lights(
                cx,
                &area,
                &three,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let timed_three = raster
            .sample_pixel_sensor_area_lights_at_time(
                cx,
                &observer,
                &three,
                0,
                clock,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(timed_three.area_sample(), &static_three);
        history.calls.store(0, Ordering::Relaxed);
        let a = raster
            .sample_pixel_sensor_area_lights_during_exposure(
                cx,
                &observer,
                &one,
                0,
                window,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let calls = history.calls.swap(0, Ordering::Relaxed);
        let b = raster
            .sample_pixel_sensor_area_lights_during_exposure(
                cx,
                &observer,
                &three,
                0,
                window,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert!(calls >= CAMERA_WAVELENGTHS);
        assert_eq!(history.calls.load(Ordering::Relaxed), calls);
        assert_eq!(a.images(), b.images());
        assert_eq!(a.trace_evaluations(), b.trace_evaluations());
        assert_eq!(a.query_work(), b.query_work());
        assert_eq!(
            a.pixel_sample().mean_pixel_sample().spectral_samples(),
            b.pixel_sample().mean_pixel_sample().spectral_samples()
        );
        near(
            b.pixel_sample().photoelectrons_per_m2_sr()[0],
            a.pixel_sample().photoelectrons_per_m2_sr()[0] * 0.25,
        );
        near(
            b.pixel_sample().photoelectrons_per_m2_sr()[1],
            a.pixel_sample().photoelectrons_per_m2_sr()[0] * 0.75,
        );
        assert_eq!(b.pixel_sample().photoelectrons_per_m2_sr()[2], 0.0);
        assert_eq!(b.pixel_sample().qe_weighted_energy_j_m2_sr()[2], 0.0);

        // At most two source evaluations occur in a lane of this plane fixture.
        // Refuse after at least one complete lane and expose no partial pixel.
        history.calls.store(0, Ordering::Relaxed);
        history.fail_at.store(3, Ordering::Relaxed);
        assert!(
            matches!(raster.sample_pixel_sensor_area_lights_at_time(cx, &observer, &three, 0,
            clock, settings(), &mut workspace, seed, pixel, sample),
            Err(TransientAreaLightRasterError::Transport { lane, error:
                TransientAreaLightCameraError::Source(TransientCameraError::History {
                    object_id: 9, error: SurfaceHistoryError::OutsideTimeDomain }) }) if lane > 0)
        );
        history.fail_at.store(0, Ordering::Relaxed);
        let replay = raster
            .sample_pixel_sensor_area_lights_at_time(
                cx,
                &observer,
                &three,
                0,
                clock,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(replay, timed_three);
        gate.request();
        assert!(
            raster
                .sample_pixel_area_lights_during_exposure(
                    cx,
                    &observer,
                    0,
                    window,
                    settings(),
                    &mut workspace,
                    seed,
                    pixel,
                    sample
                )
                .is_err()
        );
    });
}

#[test]
fn g1_area_detector_shutter_uses_each_chromatic_group_delay_and_native_energy_conversion() {
    fixtures::with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.2, 0.01, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&glass)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
                SurfaceBinding {
                    object_id: 10,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new(vec![(LAMBDA_MIN, 1.0), (LAMBDA_MAX, 5.0)]).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 2e-7,
        }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let history = LinearHistory {
            epoch,
            offset: 2.0,
            slope: 1.0 / 30e-9,
            calls: AtomicUsize::new(0),
            fail_at: AtomicUsize::new(0),
        };
        let timed =
            TransientCurvedCameraScene::try_new(cx, &camera, &[history_binding(&history)]).unwrap();
        let observer = TransientCurvedAreaLightCameraScene::try_new(cx, &timed, &area).unwrap();
        let physical = fixtures::camera(0.01);
        let settings = || scene_settings(&geometry, &physical);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let midpoint = epoch + Dd::from_f64(30e-9);
        let window =
            ExposureWindow::try_new(midpoint - Dd::from_f64(2e-9), midpoint + Dd::from_f64(2e-9))
                .unwrap();
        let baseline = raster
            .sample_pixel_sensor_area_lights(
                cx,
                &area,
                &sensor,
                0,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let exposed = raster
            .sample_pixel_sensor_area_lights_during_exposure(
                cx,
                &observer,
                &sensor,
                0,
                window,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let instant = raster
            .sample_pixel_sensor_area_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                midpoint,
                settings(),
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let mean = exposed.pixel_sample().mean_pixel_sample();
        let mut independent_power = [0.0; 3];
        let mut independent_electrons = [0.0; 3];
        let key = [seed as u32, (seed >> 32) as u32];
        for lane in 0..CAMERA_WAVELENGTHS {
            let wavelength = mean.detector_wavelengths()[lane].wavelength_nm();
            let ray = mean.camera_rays()[lane];
            let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
            let receiver = Point3::new(
                ray.origin.x + camera_length * ray.dir.x,
                ray.origin.y + camera_length * ray.dir.y,
                1.0,
            );
            let block = philox4x32_10([pixel, sample, 0, 0x6372_6100 + lane as u32 * 64], key);
            let point = source
                .sample_area([unit(block[0]), unit(block[1])])
                .unwrap()
                .position_m();
            let delta = receiver.delta_from(point);
            let source_length = delta.dot(delta).sqrt();
            let group_index = 1.2 + 3.0 * 0.01 / (wavelength * 1e-3).powi(2);
            let group_delay = (camera_length + source_length) * group_index / 299_792_458.0;
            let scale = 2.0 + (30e-9 - group_delay) / 30e-9;
            near(
                exposed.area_radiance_w_m2_sr_nm()[lane],
                baseline.area_radiance_w_m2_sr_nm()[lane] * scale,
            );
            near(
                mean.spectral_samples()[lane].radiance_w_m2_sr_nm(),
                instant.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
            );
            let radiance = mean.spectral_samples()[lane].radiance_w_m2_sr_nm();
            for (channel, qe) in [0.25, 0.75, 0.0].into_iter().enumerate() {
                // Sum-QE is one: the original uniform proposal contributes
                // its 400nm width divided across the same four paths.
                let energy = radiance * qe * 100.0 * window.duration_s().to_f64();
                independent_power[channel] += energy;
                independent_electrons[channel] +=
                    energy * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        assert_eq!(exposed.trace_evaluations(), baseline.trace_evaluations());
        assert_eq!(exposed.images(), baseline.images());
        for channel in 0..3 {
            near(
                exposed.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                independent_power[channel],
            );
            near(
                exposed.pixel_sample().photoelectrons_per_m2_sr()[channel],
                independent_electrons[channel],
            );
        }
    });
}
