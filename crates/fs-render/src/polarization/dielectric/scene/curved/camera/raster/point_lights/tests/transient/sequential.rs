//! The sensor clock and analyzer precede the full physical lens once.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;

pub(in super::super::super) const SEED: u64 = 79;
pub(in super::super::super) const SAMPLE: u32 = 6;

pub(in super::super::super) fn physical_camera(stop: f64) -> PhysicalCamera {
    let material = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(material)).unwrap(),
            LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
            LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            LensSurface::try_new(0.12, 0.0, stop, None).unwrap(),
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
pub(in super::super::super) fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; 4],
    analyzer: CameraAnalyzer,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    camera
        .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
            cx,
            (2.0 * unit(raster[0]) - 1.0) * half,
            (1.0 - 2.0 * unit(raster[1])) * half,
            LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap(),
            wavelengths,
            analyzer.reference_axis(),
            analyzer.row(),
        )
        .unwrap()
}

#[test]
fn g1_real_lens_point_pixels_keep_dd_clock_analyzer_native_work_and_clipped_qe_mass() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let material = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.37).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(),
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
            policy: BundleBranchPolicy::TransmitOrReflectTir,
            position_tolerance_m: 1e-6,
        }];
        let lights = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let history = LinearHistory::new(&source, epoch, 1e7);
        let histories = [history_binding(&source, &history)];
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let timed = TransientCurvedCameraScene::try_new(cx, &surface, &[])
            .unwrap()
            .with_analyzer(analyzer);
        let observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &histories).unwrap();
        let raw = TransientCurvedCameraScene::try_new(cx, &surface, &[]).unwrap();
        let raw_observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &raw, &lights, &histories).unwrap();
        let physical = physical_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let result = raster
            .sample_pixel_sequential_point_lights_at_time(
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
        let pixel = result.pixel_sample();
        let beams = packet(&physical, cx, pixel.wavelengths_nm(), analyzer);
        let render_cx = cx.with_stream_seed(SEED);
        let mut xyz = [0.0; 3];
        for lane in 0..4 {
            let beam = beams[lane].unwrap();
            assert!(beam.power().q != 0.0 || beam.power().u != 0.0);
            assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
            assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
            let ray = beam.ray();
            let camera_length = (2.0 - ray.origin.z) / ray.dir.z;
            let receiver = ray.origin.offset(ray.dir.scale(camera_length));
            let delta = receiver.delta_from(source.position_m());
            let distance_squared = delta.dot(delta);
            let distance = distance_squared.sqrt();
            let direction = delta.scale(1.0 / distance);
            let delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(camera_length + distance) / Dd::from_f64(LIGHT_SPEED);
            let scale = 1.0 + 1e7 * ((clock - epoch) - delay).to_f64();
            let expected = beam.power().i * 0.37 / core::f64::consts::PI
                * source
                    .spectral_intensity_w_sr_nm(direction, beam.wavelength_nm())
                    .unwrap()
                * direction.z
                / distance_squared
                * (-0.2 * (camera_length + distance)).exp()
                * scale;
            near(
                pixel.spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                expected,
            );
            let manual = raw_observer
                .sample_sequential_ray_with_samples_at_time(
                    &render_cx,
                    &beam,
                    0,
                    clock,
                    point_settings(),
                    |event| {
                        let block = philox4x32_10(
                            [0, SAMPLE, event, 0x6372_7340 + lane as u32],
                            [SEED as u32, 0],
                        );
                        core::array::from_fn(|i| unit(block[i]))
                    },
                )
                .unwrap();
            assert_eq!(
                pixel.spectral_samples()[lane],
                Some(*manual.radiance_sample())
            );
            assert_eq!(
                result.point_light_sample().connections()[lane],
                u64::from(manual.point_light_sample().connections())
            );
            assert_eq!(
                result.point_light_sample().trace_evaluations()[lane],
                u64::from(manual.point_light_sample().trace_evaluations())
            );
            assert_eq!(
                result.point_light_sample().field_evaluations()[lane],
                manual.point_light_sample().field_evaluations()
            );
            assert_eq!(
                result.point_light_sample().query_work()[lane],
                manual.point_light_sample().query_work()
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
            near(pixel.xyz()[channel], xyz[channel]);
        }
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = raster
            .sample_pixel_sequential_point_lights_during_exposure(
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
        let middle = raster
            .sample_pixel_sequential_point_lights_at_time(
                cx,
                &observer,
                0,
                clock + Dd::from_f64(2e-9),
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        for channel in 0..3 {
            near(
                exposed.pixel_sample().mean_pixel_sample().xyz()[channel],
                middle.pixel_sample().xyz()[channel],
            );
        }
        assert_eq!(
            exposed.trace_evaluations(),
            result.point_light_sample().trace_evaluations()
        );

        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let full = raster
            .sample_pixel_sensor_sequential_point_lights_at_time(
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
        let wavelengths = full
            .pixel_sample()
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let sensor_beams = packet(&physical, cx, &wavelengths, analyzer);
        let radii = sensor_beams.map(|beam| {
            let ray = beam.unwrap().ray();
            let d = (0.12 - ray.origin.z) / ray.dir.z;
            (ray.origin.x + d * ray.dir.x).hypot(ray.origin.y + d * ray.dir.y)
        });
        let lo = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = radii.iter().copied().fold(0.0, f64::max);
        assert!(hi - lo > 1e-10);
        let stop = (lo + hi) * 0.5;
        let partial_camera = physical_camera(stop);
        let partial = CurvedCameraRaster::try_new(&partial_camera, 1, 1)
            .unwrap()
            .sample_pixel_sensor_sequential_point_lights_at_time(
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
        let mut survivors = 0;
        let mut energy = [0.0; 3];
        let mut photons = [0.0; 3];
        for lane in 0..4 {
            let present = radii[lane] < stop;
            assert_eq!(
                partial.pixel_sample().spectral_samples()[lane].is_some(),
                present
            );
            if !present {
                assert_eq!(partial.point_light_sample().connections()[lane], 0);
                assert_eq!(partial.point_light_sample().trace_evaluations()[lane], 0);
                assert_eq!(partial.point_light_sample().field_evaluations()[lane], 0);
                assert_eq!(partial.point_light_sample().query_work()[lane], 0);
                continue;
            }
            survivors += 1;
            let original = full.pixel_sample().spectral_samples()[lane].unwrap();
            assert_eq!(
                partial.pixel_sample().spectral_samples()[lane],
                Some(original)
            );
            assert_eq!(
                partial.point_light_sample().trace_evaluations()[lane],
                full.point_light_sample().trace_evaluations()[lane]
            );
            for channel in 0..3 {
                let contribution = original.radiance_w_m2_sr_nm()
                    * sensor_beams[lane].unwrap().geometry_weight()
                    * [0.25, 0.75, 0.0][channel]
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0;
                energy[channel] += contribution;
                photons[channel] +=
                    contribution * wavelengths[lane] * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
            }
        }
        assert!(survivors > 0 && survivors < 4);
        for channel in 0..3 {
            near(
                partial.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                energy[channel],
            );
            near(
                partial.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                photons[channel],
            );
        }
        let integrated = raster
            .sample_pixel_sensor_sequential_point_lights_during_exposure(
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
        let middle = raster
            .sample_pixel_sensor_sequential_point_lights_at_time(
                cx,
                &observer,
                &sensor,
                0,
                clock + Dd::from_f64(2e-9),
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(
            integrated.trace_evaluations(),
            full.point_light_sample().trace_evaluations()
        );
        for channel in 0..3 {
            near(
                integrated.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                middle.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel] * 4e-9,
            );
            near(
                integrated.pixel_sample().photoelectrons_per_m2_sr()[channel],
                middle.pixel_sample().photoelectron_rate_per_m2_sr()[channel] * 4e-9,
            );
        }

        let closed_camera = physical_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let before = history.calls.load(Ordering::Relaxed);
        let dark = closed
            .sample_pixel_sensor_sequential_point_lights_at_time(
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
        assert_eq!(dark.observation_time_s(), clock);
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(dark.point_light_sample().connections(), &[0; 4]);
        assert_eq!(dark.point_light_sample().trace_evaluations(), &[0; 4]);
        assert_eq!(history.calls.load(Ordering::Relaxed), before);
        assert!(matches!(
            closed.sample_pixel_sequential_point_lights_at_time(
                cx,
                &observer,
                0,
                Dd::from_f64(f64::NAN),
                point_settings(),
                SEED,
                0,
                SAMPLE
            ),
            Err(TransientPointLightRasterError::Source(
                TransientCameraError::InvalidObservationTime
            ))
        ));
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_point_lights_at_time(
                    cx,
                    &observer,
                    &sensor,
                    0,
                    clock,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}
