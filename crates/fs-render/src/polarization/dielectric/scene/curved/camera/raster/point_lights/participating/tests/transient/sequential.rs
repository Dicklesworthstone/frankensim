//! Real-lens point/image pixels retain native camera and source null tracking.

use super::*;
use crate::gradient_index::AffineIndexField;
use crate::polarization::dielectric::scene::curved::camera::raster::tests::enclosed_geometry;
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use crate::polarization::dielectric::scene::curved::camera::raster::point_lights::tests::transient::sequential::{
    SAMPLE, SEED, packet, physical_camera,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighBinding, CurvedRayleighField, CurvedRayleighSample,
};
use crate::volumes::transport::VolumeTransportError;
use fs_blake3::ContentHash;
use fs_geom::Aabb;

const FIELD_ID: ContentHash = ContentHash([128; 32]);

// An empty cell inside a positive authored majorant still consumes native
// camera/null and source-survival draws. Its exact physical attenuation is one.
struct EmptyCell;
impl CurvedRayleighField for EmptyCell {
    fn source_field_identity(&self) -> ContentHash {
        FIELD_ID
    }
    fn scattering_majorant_per_m(&self, _: f64) -> Result<f64, VolumeTransportError> {
        Ok(1.0)
    }
    fn sample_scattering(
        &self,
        _: Point3,
        _: f64,
    ) -> Result<CurvedRayleighSample, VolumeTransportError> {
        Ok(CurvedRayleighSample {
            scattering_per_m: 0.0,
            depolarization_ratio: 0.0,
        })
    }
}

#[test]
fn g1_molecular_lens_point_images_replay_all_native_domains_and_share_qe_shutter_work() {
    with_cx(|cx, _| {
        // In the open plane scene the camera's diffuse continuation from the receiver runs
        // down through the empty cell's positive majorant, is never certified to escape, and
        // ends at the distance limit: Camera(Unresolved(DistanceLimit)). The black box of
        // enclosed_geometry() closes the scene below; the continuation ends on it with no
        // radiance, so the expected radiance is unchanged. The 10 m trace budget covers
        // every camera path: about 2 m to the receiver, then at most the box diagonal,
        // sqrt(41) = 6.4 m.
        let geometry = enclosed_geometry();
        let material = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let field = AffineIndexField::try_new(
            &material,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Aabb {
                min: Point3::new(-50.0, -50.0, -50.0),
                max: Point3::new(50.0, 50.0, 50.0),
            },
        )
        .unwrap();
        let media = [SceneMedium::Differentiable {
            field: &field,
            identity: FIELD_ID,
        }];
        let reflectances = [SpectralReflectance::constant(0.37).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 7,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                SurfaceBinding {
                    object_id: 8,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let source = fixture_source();
        let branches = [PointLightBranch {
            source: &source,
            receiver_object_id: 7,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::Sequence(&[]),
            position_tolerance_m: 1e-6,
        }];
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
        let empty_fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
        let cell = EmptyCell;
        let laws = [Some(CurvedRayleighBinding::Field(&cell))];
        let scattering = CurvedRayleighScene::try_new_with_fields(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let lit = CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let unbound =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &empty_fixtures).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&[]),
        }; 2];
        let groups = [PointLightImageGroup {
            source: &source,
            receiver: PointLightImageReceiver::SurfaceObject(7),
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 2,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 1e-6,
        }];
        let images = CurvedRayleighImageLightCameraScene::try_new(cx, &unbound, &groups).unwrap();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(2, 2).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = LinearHistory::new(&source, epoch, 1e7);
        let histories = [history_binding(&source, &history)];
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let timed = TransientCurvedCameraScene::try_new(cx, &surface, &[])
            .unwrap()
            .with_analyzer(analyzer);
        let observer =
            TransientCurvedRayleighPointLightCameraScene::try_new(cx, &timed, &lit, &histories)
                .unwrap();
        let image_observer =
            TransientCurvedRayleighImageLightCameraScene::try_new(cx, &timed, &images, &histories)
                .unwrap();
        let raw = TransientCurvedCameraScene::try_new(cx, &surface, &[]).unwrap();
        let raw_observer =
            TransientCurvedRayleighImageLightCameraScene::try_new(cx, &raw, &images, &histories)
                .unwrap();
        let physical = physical_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let result = raster
            .sample_pixel_sequential_rayleigh_point_lights_at_time(
                cx,
                &observer,
                0,
                clock,
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let discovered = raster
            .sample_pixel_sequential_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                0,
                clock,
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(result.molecular_sample().volume_events(), &[0; 4]);
        assert_eq!(
            result.molecular_sample().volume_events(),
            discovered.molecular_sample().volume_events()
        );
        assert_eq!(
            result.molecular_sample().free_flight_draws(),
            discovered.molecular_sample().free_flight_draws()
        );
        assert_eq!(
            result.molecular_sample().source_free_flight_draws(),
            discovered.molecular_sample().source_free_flight_draws()
        );
        assert_eq!(
            result.molecular_sample().source_candidates(),
            discovered.molecular_sample().source_candidates()
        );
        assert!(
            result
                .molecular_sample()
                .source_free_flight_draws()
                .iter()
                .all(|draws| *draws > 0)
        );
        let beams = packet(
            &physical,
            cx,
            result.pixel_sample().wavelengths_nm(),
            analyzer,
        );
        let render_cx = cx.with_stream_seed(SEED);
        let mut xyz = [0.0; 3];
        for lane in 0..4 {
            let beam = beams[lane].unwrap();
            let ray = beam.ray();
            let camera_length = (2.0 - ray.origin.z) / ray.dir.z;
            let receiver = ray.origin.offset(ray.dir.scale(camera_length));
            let delta = receiver.delta_from(source.position_m());
            let distance_squared = delta.dot(delta);
            let distance = distance_squared.sqrt();
            let direction = delta.scale(1.0 / distance);
            let delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(camera_length + distance) / Dd::from_f64(LIGHT_SPEED);
            let expected = beam.power().i * 0.37 / core::f64::consts::PI
                * source
                    .spectral_intensity_w_sr_nm(direction, beam.wavelength_nm())
                    .unwrap()
                * direction.z
                / distance_squared
                * (-0.2 * (camera_length + distance)).exp()
                * (1.0 + 1e7 * ((clock - epoch) - delay).to_f64());
            near(
                result.pixel_sample().spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                expected,
            );
            near(
                discovered.pixel_sample().spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                expected,
            );
            let manual = raw_observer
                .sample_sequential_ray_with_tracking_samples_at_time(
                    &render_cx,
                    &beam,
                    0,
                    clock,
                    point_settings(),
                    &mut workspace,
                    |event| {
                        let block = philox4x32_10(
                            [0, SAMPLE, event, 0x6372_7340 + lane as u32],
                            [SEED as u32, 0],
                        );
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10(
                            [0, SAMPLE, flight, 0x6372_7350 + lane as u32],
                            [SEED as u32, 0],
                        )
                        .map(unit)
                    },
                    |connection, flight| {
                        unit(
                            philox4x32_10(
                                [0, SAMPLE, connection, flight],
                                [SEED as u32, 0x6372_7360 + lane as u32],
                            )[0],
                        )
                    },
                )
                .unwrap();
            let work = manual.point_light_sample();
            let pixel_work = discovered.molecular_sample();
            assert_eq!(
                discovered.pixel_sample().spectral_samples()[lane],
                Some(*manual.radiance_sample())
            );
            assert_eq!(
                pixel_work.point_light_sample().connections()[lane],
                u64::from(work.connections())
            );
            assert_eq!(
                pixel_work.point_light_sample().trace_evaluations()[lane],
                u64::from(work.trace_evaluations())
            );
            assert_eq!(
                pixel_work.point_light_sample().field_evaluations()[lane],
                work.field_evaluations()
            );
            assert_eq!(
                pixel_work.point_light_sample().query_work()[lane],
                work.query_work()
            );
            assert_eq!(
                pixel_work.free_flight_draws()[lane],
                work.free_flight_draws()
            );
            assert_eq!(
                pixel_work.source_free_flight_draws()[lane],
                work.source_free_flight_draws()
            );
            assert_eq!(
                pixel_work.source_candidates()[lane],
                work.source_candidates()
            );
            for (sum, response) in xyz.iter_mut().zip([
                cie_x(beam.wavelength_nm()),
                cie_y(beam.wavelength_nm()),
                cie_z(beam.wavelength_nm()),
            ]) {
                *sum += expected * beam.geometry_weight() * response * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0
                    / y_integral();
            }
        }
        for channel in 0..3 {
            near(result.pixel_sample().xyz()[channel], xyz[channel]);
            near(discovered.pixel_sample().xyz()[channel], xyz[channel]);
        }
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposure = raster
            .sample_pixel_sequential_rayleigh_point_lights_during_exposure(
                cx,
                &observer,
                0,
                window,
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let image_exposure = raster
            .sample_pixel_sequential_rayleigh_image_lights_during_exposure(
                cx,
                &image_observer,
                0,
                window,
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let middle = raster
            .sample_pixel_sequential_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                0,
                clock + Dd::from_f64(2e-9),
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        for channel in 0..3 {
            near(
                exposure.pixel_sample().mean_pixel_sample().xyz()[channel],
                middle.pixel_sample().xyz()[channel],
            );
            near(
                image_exposure.pixel_sample().mean_pixel_sample().xyz()[channel],
                middle.pixel_sample().xyz()[channel],
            );
        }
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let detected = raster
            .sample_pixel_sensor_sequential_rayleigh_point_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let detected_images = raster
            .sample_pixel_sensor_sequential_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                &sensor,
                0,
                clock,
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let energy = raster
            .sample_pixel_sensor_sequential_rayleigh_point_lights_during_exposure(
                cx,
                &observer,
                &sensor,
                0,
                window,
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let image_energy = raster
            .sample_pixel_sensor_sequential_rayleigh_image_lights_during_exposure(
                cx,
                &image_observer,
                &sensor,
                0,
                window,
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let middle = raster
            .sample_pixel_sensor_sequential_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                &sensor,
                0,
                clock + Dd::from_f64(2e-9),
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        for channel in 0..3 {
            near(
                detected.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                detected_images
                    .pixel_sample()
                    .qe_weighted_radiance_w_m2_sr()[channel],
            );
            near(
                energy.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                middle.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel] * 4e-9,
            );
            near(
                image_energy.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                middle.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel] * 4e-9,
            );
            near(
                image_energy.pixel_sample().photoelectrons_per_m2_sr()[channel],
                middle.pixel_sample().photoelectron_rate_per_m2_sr()[channel] * 4e-9,
            );
        }
        assert_eq!(
            energy.source_free_flight_draws(),
            detected.molecular_sample().source_free_flight_draws()
        );
        assert_eq!(
            image_energy.source_candidates(),
            detected_images.molecular_sample().source_candidates()
        );
        let closed_camera = physical_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1)
            .unwrap()
            .sample_pixel_sensor_sequential_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                &sensor,
                0,
                clock,
                point_settings(),
                &mut workspace,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(closed.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            closed.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(closed.molecular_sample().free_flight_draws(), &[0; 4]);
        assert_eq!(
            closed.molecular_sample().source_free_flight_draws(),
            &[0; 4]
        );
        assert_eq!(closed.molecular_sample().source_candidates(), &[0; 4]);
        assert_eq!(
            closed.molecular_sample().point_light_sample().connections(),
            &[0; 4]
        );
        assert_eq!(
            closed
                .molecular_sample()
                .point_light_sample()
                .trace_evaluations(),
            &[0; 4]
        );
    });
}
