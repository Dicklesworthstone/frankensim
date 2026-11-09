//! Original molecular vertices, source-image survival and all three counter domains.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::raster::{LAMBDA_MAX, LAMBDA_MIN};
use crate::polarization::dielectric::scene::curved::camera::raster::point_lights::tests::transient::{
    LIGHT_SPEED, LinearHistory, fixture_source, history_binding, near,
};
use crate::polarization::dielectric::scene::curved::camera::raster::transient::point_lights::TransientPointLightRasterError;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::participating::TransientCurvedRayleighPointLightCameraScene;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::participating::images::TransientCurvedRayleighImageLightCameraScene;
use crate::polarization::dielectric::scene::curved::camera::transient::TransientCurvedCameraScene;
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_math::dd::Dd;

mod sequential;

#[test]
fn g1_molecular_point_and_image_pixels_keep_retarded_single_scattering_and_shutter_energy() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                9,
                hash_domain("molecular-point-pixel", b"black-enclosure"),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.6,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let material = material();
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(
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
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
        let spectrum = MassScatteringSpectrum::constant(0.2).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let source = fixture_source();
        let branches = [MolecularPointLightBranch {
            source: &source,
            receiver_medium: 0,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::Sequence(&[]),
            position_tolerance_m: 1e-6,
        }];
        let lit = CurvedRayleighPointLightCameraScene::try_new_with_volume_branches(
            cx, &molecular, &fixtures, &branches,
        )
        .unwrap();
        let unbound =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&[]),
        }; 2];
        let groups = [PointLightImageGroup {
            source: &source,
            receiver: PointLightImageReceiver::MolecularMedium(0),
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
        let timed = TransientCurvedCameraScene::try_new(cx, &surface, &[]).unwrap();
        let observer =
            TransientCurvedRayleighPointLightCameraScene::try_new(cx, &timed, &lit, &histories)
                .unwrap();
        let image_observer =
            TransientCurvedRayleighImageLightCameraScene::try_new(cx, &timed, &images, &histories)
                .unwrap();
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        // Original established fixture: one real collision and two flights per lane.
        let (seed, pixel_id, sample_id) = (120_u64, 4, 9);
        let ordinary = raster
            .sample_pixel_rayleigh_point_lights(
                cx,
                &lit,
                0,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let original_images = raster
            .sample_pixel_rayleigh_image_lights(
                cx,
                &images,
                0,
                point_settings(),
                &mut workspace,
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let result = raster
            .sample_pixel_rayleigh_point_lights_at_time(
                cx,
                &observer,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let discovered = raster
            .sample_pixel_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                0,
                clock,
                point_settings(),
                &mut workspace,
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        assert_eq!(result.observation_time_s(), clock);
        assert_eq!(discovered.observation_time_s(), clock);
        assert_eq!(result.molecular_sample().volume_events(), &[1; 4]);
        assert_eq!(result.molecular_sample().free_flight_draws(), &[2; 4]);
        assert_eq!(
            result.molecular_sample().source_free_flight_draws(),
            &[0; 4]
        );
        assert_eq!(result.molecular_sample().source_candidates(), &[0; 4]);
        assert_eq!(
            result
                .molecular_sample()
                .point_light_sample()
                .trace_evaluations(),
            ordinary.point_light_sample().trace_evaluations()
        );
        assert_eq!(
            discovered
                .molecular_sample()
                .point_light_sample()
                .trace_evaluations(),
            original_images.point_light_sample().trace_evaluations()
        );
        assert_eq!(
            result.molecular_sample().volume_events(),
            discovered.molecular_sample().volume_events()
        );
        assert_eq!(
            result.molecular_sample().free_flight_draws(),
            discovered.molecular_sample().free_flight_draws()
        );
        let ray = result.pixel_sample().camera_ray();
        let expected = |lane: usize, wavelength: f64, observation: Dd| {
            let block = philox4x32_10(
                [pixel_id, sample_id, 0, 0x6372_7350 + lane as u32],
                [seed as u32, 0],
            );
            let distance = -(1.0 - unit(block[0])).ln() / 0.2;
            assert!((0.05..1.4).contains(&distance));
            let vertex = ray.origin.offset(ray.dir.scale(distance));
            let source_delta = vertex.delta_from(source.position_m());
            let source_length = source_delta.dot(source_delta).sqrt();
            let direction = source_delta.scale(1.0 / source_length);
            let cosine = -direction.dot(ray.dir);
            let phase_per_sr = 3.0 * (1.0 + cosine * cosine) / (16.0 * core::f64::consts::PI);
            let delay = 1.25 * (distance + source_length) / LIGHT_SPEED;
            let scale = 1.0 + 1e7 * ((observation - epoch).to_f64() - delay);
            // Camera survival is sampled, while the independent source leg
            // receives exp(-beta_s*s) once. No phase-PDF divisor or delta MIS.
            phase_per_sr
                * source
                    .spectral_intensity_w_sr_nm(direction, wavelength)
                    .unwrap()
                / source_length.powi(2)
                * (-0.1 * distance - 0.3 * source_length).exp()
                * scale
        };
        let mut xyz = [0.0; 3];
        let render_cx = cx.with_stream_seed(seed);
        for (lane, &wavelength) in result.pixel_sample().wavelengths_nm().iter().enumerate() {
            let radiance = expected(lane, wavelength, clock);
            near(
                result.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                radiance,
            );
            near(
                discovered.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                radiance,
            );
            assert_eq!(
                result.pixel_sample().spectral_samples()[lane].path(),
                ordinary.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(
                discovered.pixel_sample().spectral_samples()[lane].path(),
                original_images.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(
                discovered
                    .molecular_sample()
                    .point_light_sample()
                    .connections()[lane],
                1
            );
            assert!(
                discovered
                    .molecular_sample()
                    .point_light_sample()
                    .trace_evaluations()[lane]
                    > result
                        .molecular_sample()
                        .point_light_sample()
                        .trace_evaluations()[lane]
            );
            let manual = image_observer
                .sample_ray_with_tracking_samples_at_time(
                    &render_cx,
                    ray,
                    wavelength,
                    0,
                    clock,
                    point_settings(),
                    &mut workspace,
                    |event| {
                        let block = philox4x32_10(
                            [pixel_id, sample_id, event, 0x6372_7340 + lane as u32],
                            [seed as u32, 0],
                        );
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |flight| {
                        philox4x32_10(
                            [pixel_id, sample_id, flight, 0x6372_7350 + lane as u32],
                            [seed as u32, 0],
                        )
                        .map(unit)
                    },
                    |connection, flight| {
                        unit(
                            philox4x32_10(
                                [pixel_id, sample_id, connection, flight],
                                [seed as u32, 0x6372_7360 + lane as u32],
                            )[0],
                        )
                    },
                )
                .unwrap();
            assert_eq!(
                discovered.pixel_sample().spectral_samples()[lane],
                *manual.radiance_sample()
            );
            assert_eq!(
                discovered
                    .molecular_sample()
                    .point_light_sample()
                    .query_work()[lane],
                manual.point_light_sample().query_work()
            );
            for (sum, response) in
                xyz.iter_mut()
                    .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
            {
                *sum += radiance * response * (LAMBDA_MAX - LAMBDA_MIN) / 4.0 / y_integral();
            }
        }
        for channel in 0..3 {
            near(result.pixel_sample().xyz()[channel], xyz[channel]);
            near(discovered.pixel_sample().xyz()[channel], xyz[channel]);
        }
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let integrated = raster
            .sample_pixel_rayleigh_point_lights_during_exposure(
                cx,
                &observer,
                0,
                window,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let image_integrated = raster
            .sample_pixel_rayleigh_image_lights_during_exposure(
                cx,
                &image_observer,
                0,
                window,
                point_settings(),
                &mut workspace,
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        for lane in 0..4 {
            let wavelength = result.pixel_sample().wavelengths_nm()[lane];
            let radiance = expected(lane, wavelength, clock + Dd::from_f64(2e-9));
            near(
                integrated
                    .pixel_sample()
                    .mean_pixel_sample()
                    .spectral_samples()[lane]
                    .radiance_w_m2_sr_nm(),
                radiance,
            );
            near(
                image_integrated
                    .pixel_sample()
                    .mean_pixel_sample()
                    .spectral_samples()[lane]
                    .radiance_w_m2_sr_nm(),
                radiance,
            );
        }
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let detected = raster
            .sample_pixel_sensor_rayleigh_point_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let detected_images = raster
            .sample_pixel_sensor_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                &sensor,
                0,
                clock,
                point_settings(),
                &mut workspace,
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let energy = raster
            .sample_pixel_sensor_rayleigh_point_lights_during_exposure(
                cx,
                &observer,
                &sensor,
                0,
                window,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let image_energy = raster
            .sample_pixel_sensor_rayleigh_image_lights_during_exposure(
                cx,
                &image_observer,
                &sensor,
                0,
                window,
                point_settings(),
                &mut workspace,
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let mut joules = [0.0; 3];
        let mut electrons = [0.0; 3];
        for (lane, draw) in detected
            .pixel_sample()
            .detector_wavelengths()
            .iter()
            .enumerate()
        {
            let wavelength = draw.wavelength_nm();
            near(
                detected.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                expected(lane, wavelength, clock),
            );
            near(
                detected_images.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                expected(lane, wavelength, clock),
            );
            let mean = expected(lane, wavelength, clock + Dd::from_f64(2e-9));
            for channel in 0..3 {
                let contribution =
                    mean * [0.25, 0.75, 0.0][channel] * (LAMBDA_MAX - LAMBDA_MIN) / 4.0 * 4e-9;
                joules[channel] += contribution;
                electrons[channel] +=
                    contribution * wavelength * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
            }
        }
        for channel in 0..3 {
            near(
                energy.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                joules[channel],
            );
            near(
                image_energy.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                joules[channel],
            );
            near(
                energy.pixel_sample().photoelectrons_per_m2_sr()[channel],
                electrons[channel],
            );
            near(
                image_energy.pixel_sample().photoelectrons_per_m2_sr()[channel],
                electrons[channel],
            );
        }
        assert_eq!(
            energy.volume_events(),
            detected.molecular_sample().volume_events()
        );
        assert_eq!(
            image_energy.point_light_sample().trace_evaluations(),
            detected_images
                .molecular_sample()
                .point_light_sample()
                .trace_evaluations()
        );
        let mut too_small = CurvedTransmittanceImageWorkspace::try_new(1, 1).unwrap();
        assert!(matches!(
            raster.sample_pixel_rayleigh_image_lights_at_time(
                cx,
                &image_observer,
                0,
                clock,
                point_settings(),
                &mut too_small,
                seed,
                pixel_id,
                sample_id
            ),
            Err(TransientPointLightRasterError::Transport { lane: 0, .. })
        ));
    });
}
