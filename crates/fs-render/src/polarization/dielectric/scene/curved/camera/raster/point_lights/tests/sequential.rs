//! Active coated lenses reach steady point CIE/QE without a fictitious clock.

use super::transient::sequential::{SAMPLE, SEED};
use super::transient::{LIGHT_SPEED, fixture_source, near};
use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::DielectricCoating;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::analyzer::AnalyzerResponse;
use crate::polarization::dielectric::scene::curved::camera::raster::{LAMBDA_MAX, LAMBDA_MIN};
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::{
    TransientCurvedPointLightCameraScene, TransientPointLightCameraError,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    TransientCameraError, TransientCurvedCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, SpectralCameraAnalyzer,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_blake3::hash_domain;
use fs_geom::Point3;
use fs_math::dd::Dd;

fn coated_camera(stop_radius: f64) -> PhysicalCamera {
    let film = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.02, 110e-9)
                    .unwrap(),
            ],
            1,
            ConductorSource::try_new(
                hash_domain(
                    "org.frankensim.test.steady-point-pixel",
                    b"authored absorbing coating",
                ),
                ConductorDataStatus::Representative,
            )
            .unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(glass))
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.09, -20.0, 0.012, None)
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.12, 0.0, stop_radius, None).unwrap(),
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

fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; 4],
    response: Option<AnalyzerResponse<'_>>,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    let x = (2.0 * unit(raster[0]) - 1.0) * half;
    let y = (1.0 - 2.0 * unit(raster[1])) * half;
    let lens = LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap();
    if let Some(response) = response {
        let rows = wavelengths.map(|wavelength| response.at(wavelength).unwrap().row());
        let axis = response.at(wavelengths[0]).unwrap().reference_axis();
        camera
            .generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets(
                cx,
                x,
                y,
                lens,
                wavelengths,
                axis,
                &rows,
            )
            .unwrap()
    } else {
        camera
            .generate_sequential_spectral_beams_from_tangent_offsets(cx, x, y, lens, wavelengths)
            .unwrap()
    }
}

fn expected_point(source: &PhotometricSource, beam: SequentialCameraBeam, scattering: f64) -> f64 {
    let ray = beam.ray();
    let camera_length = (2.0 - ray.origin.z) / ray.dir.z;
    let receiver = ray.origin.offset(ray.dir.scale(camera_length));
    let delta = receiver.delta_from(source.position_m());
    let distance_squared = delta.dot(delta);
    let distance = distance_squared.sqrt();
    let direction = delta.scale(1.0 / distance);
    // The lens row already contains coating polarization and absorption.
    // Conditional camera survival and point delta mass introduce no PDF factor.
    beam.power().i * 0.37 / core::f64::consts::PI
        * source
            .spectral_intensity_w_sr_nm(direction, beam.wavelength_nm())
            .unwrap()
        * direction.z
        / distance_squared
        * (-0.2 * (camera_length + distance)).exp()
        * (-scattering * distance).exp()
}

#[test]
fn g1_steady_active_coating_keeps_point_radiometry_detector_order_and_native_work() {
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
            policy: BundleBranchPolicy::Sequence(&[]),
            position_tolerance_m: 1e-6,
        }];
        let lights = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
        let physical = coated_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let analyzer = CameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            Stokes {
                i: 0.55,
                q: 0.2,
                u: -0.1,
                v: 0.05,
            },
        )
        .unwrap();
        let spectral = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (
                    LAMBDA_MIN,
                    Stokes {
                        i: 0.45,
                        q: 0.1,
                        u: 0.0,
                        v: 0.05,
                    },
                ),
                (
                    LAMBDA_MAX,
                    Stokes {
                        i: 0.65,
                        q: -0.15,
                        u: 0.08,
                        v: -0.03,
                    },
                ),
            ],
        )
        .unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let render_cx = cx.with_stream_seed(SEED);
        let mut cie_work = None;
        let mut qe_work = None;
        for mode in 0..3 {
            let response = match mode {
                0 => None,
                1 => Some(AnalyzerResponse::Constant(analyzer)),
                _ => Some(AnalyzerResponse::Spectral(&spectral)),
            };
            let visible = match mode {
                0 => raster.sample_pixel_sequential_point_lights(
                    cx,
                    &lights,
                    0,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
                1 => raster.sample_pixel_sequential_point_lights_analyzed(
                    cx,
                    &lights,
                    0,
                    analyzer,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
                _ => raster.sample_pixel_sequential_point_lights_spectrally_analyzed(
                    cx,
                    &lights,
                    0,
                    &spectral,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
            }
            .unwrap();
            let pixel = visible.pixel_sample();
            let beams = packet(&physical, cx, pixel.wavelengths_nm(), response);
            let mut xyz = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                assert_eq!(beam.group_delay_s(), None);
                assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
                assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
                let expected = expected_point(&source, beam, 0.0);
                assert!(expected > 0.0);
                near(
                    pixel.spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                let manual = lights
                    .sample_sequential_ray_with_samples(
                        &render_cx,
                        &beam,
                        0,
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
                assert_eq!(visible.connections()[lane], u64::from(manual.connections()));
                assert_eq!(
                    visible.trace_evaluations()[lane],
                    u64::from(manual.trace_evaluations())
                );
                assert_eq!(
                    visible.field_evaluations()[lane],
                    manual.field_evaluations()
                );
                assert_eq!(visible.query_work()[lane], manual.query_work());
                for (sum, observer) in xyz.iter_mut().zip([
                    cie_x(beam.wavelength_nm()),
                    cie_y(beam.wavelength_nm()),
                    cie_z(beam.wavelength_nm()),
                ]) {
                    *sum +=
                        expected * beam.geometry_weight() * observer * (LAMBDA_MAX - LAMBDA_MIN)
                            / 4.0
                            / y_integral();
                }
            }
            for channel in 0..3 {
                near(pixel.xyz()[channel], xyz[channel]);
            }
            let work = (
                *visible.connections(),
                *visible.trace_evaluations(),
                *visible.query_work(),
            );
            if let Some(original) = cie_work {
                assert_eq!(work, original);
            }
            cie_work = Some(work);
            let detected = match mode {
                0 => raster.sample_pixel_sensor_sequential_point_lights(
                    cx,
                    &lights,
                    &sensor,
                    0,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
                1 => raster.sample_pixel_sensor_sequential_point_lights_analyzed(
                    cx,
                    &lights,
                    &sensor,
                    0,
                    analyzer,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
                _ => raster.sample_pixel_sensor_sequential_point_lights_spectrally_analyzed(
                    cx,
                    &lights,
                    &sensor,
                    0,
                    &spectral,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE,
                ),
            }
            .unwrap();
            let pixel = detected.pixel_sample();
            let wavelengths = pixel
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(&physical, cx, &wavelengths, response);
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                assert_eq!(beam.group_delay_s(), None);
                let expected = expected_point(&source, beam, 0.0);
                near(
                    pixel.spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                for channel in 0..3 {
                    let contribution = expected
                        * beam.geometry_weight()
                        * [0.25, 0.75, 0.0][channel]
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0;
                    energy[channel] += contribution;
                    electrons[channel] +=
                        contribution * wavelengths[lane] * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
                }
            }
            for channel in 0..3 {
                near(
                    pixel.qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    pixel.photoelectron_rate_per_m2_sr()[channel],
                    electrons[channel],
                );
            }
            let work = (
                *detected.connections(),
                *detected.trace_evaluations(),
                *detected.query_work(),
            );
            if let Some(original) = qe_work {
                assert_eq!(work, original);
            }
            qe_work = Some(work);
        }

        // A steady coated beam is still unqualified for source arrival time.
        let timed = TransientCurvedCameraScene::try_new(cx, &surface, &[]).unwrap();
        let observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &[]).unwrap();
        let beam = packet(&physical, cx, &[450.0, 550.0, 650.0, 750.0], None)[0].unwrap();
        assert!(matches!(
            observer.sample_sequential_ray_with_samples_at_time(
                cx,
                &beam,
                0,
                Dd::ZERO,
                point_settings(),
                |_| panic!("missing group model must precede sampling")
            ),
            Err(TransientPointLightCameraError::Source(
                TransientCameraError::MissingLensGroupDelay
            ))
        ));
        assert!(
            raster
                .sample_pixel_sequential_point_lights_at_time(
                    cx,
                    &observer,
                    0,
                    Dd::ZERO,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );

        #[cfg(feature = "volumes")]
        {
            use crate::polarization::dielectric::scene::curved::bundle::connection::images::{BundleImageOptions, BundleImageSeed};
            use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::{CurvedRayleighPointLightCameraScene, RayleighPointLightCameraSettings};
            use crate::polarization::dielectric::scene::curved::camera::point_lights::participating::images::{CurvedRayleighImageLightCameraScene, PointLightImageGroup, PointLightImageReceiver};
            use crate::polarization::dielectric::scene::curved::participating::{CurvedRayleighMedium, CurvedRayleighScene};
            use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
            use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
            use crate::polarization::dielectric::scene::curved::camera::raster::tests::{
                distance_budget_m, enclosed_geometry,
            };

            // In the open plane scene the camera's diffuse continuation from the receiver
            // runs down into unbounded scattering medium, is never certified to escape, and
            // ends at the distance limit: Camera(Unresolved(DistanceLimit)). The scattering
            // block closes the scene below with the black box of enclosed_geometry(). The
            // continuation ends on the box with no radiance, so the expected point radiance
            // is unchanged. The camera path's trace budget is the scene's own
            // (distance_budget_m, as in the area-light raster tests): 4 m does not reach from
            // the camera past the receiver to the box's far corners.
            let geometry = enclosed_geometry();
            let mut enclosed_bindings = bindings();
            enclosed_bindings.push(SurfaceBinding {
                object_id: 8,
                interaction: Interaction::Opaque,
            });
            let native = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                enclosed_bindings,
                100_000,
            )
            .unwrap();
            let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
            let lights = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
            let mut point_lights = point_settings();
            point_lights.surfaces.max_distance_m =
                distance_budget_m(&geometry, physical.eye(), point_lights.surfaces);

            // Choose a positive homogeneous rate whose first two original
            // camera flights exceed twice the box diagonal in every lane. The
            // camera path's two segments (camera to receiver, receiver to box)
            // lie in the box, so neither collides. Only the independent source
            // segment then has explicit survival.
            let diagonal = (4.0_f64 * 4.0 + 4.0 * 4.0 + 3.0 * 3.0).sqrt();
            let mut rate = 0.01_f64;
            for lane in 0..4 {
                for flight in 0..2 {
                    let draw =
                        philox4x32_10([0, SAMPLE, flight, 0x6372_7350 + lane], [SEED as u32, 0]);
                    rate = rate.min(-(1.0 - unit(draw[0])).ln() / (2.0 * diagonal));
                }
            }
            assert!(rate > 0.0);
            let spectrum = MassScatteringSpectrum::constant(rate).unwrap();
            let phase = PhaseFunction::Rayleigh;
            let laws = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
            let lit =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &lights).unwrap();
            let empty = CurvedPointLightCameraScene::try_new(cx, &surface, &[]).unwrap();
            let unbound =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &empty).unwrap();
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
            let images =
                CurvedRayleighImageLightCameraScene::try_new(cx, &unbound, &groups).unwrap();
            let mut workspace = CurvedTransmittanceImageWorkspace::try_new(2, 2).unwrap();
            let settings = RayleighPointLightCameraSettings {
                point_lights,
                max_volume_events: 8,
                max_source_free_flight_draws: 256,
                max_source_candidates: 128,
            };
            let visible = raster
                .sample_pixel_sequential_rayleigh_point_lights_spectrally_analyzed(
                    cx, &lit, 0, &spectral, settings, SEED, 0, SAMPLE,
                )
                .unwrap();
            let image_visible = raster
                .sample_pixel_sequential_rayleigh_image_lights_spectrally_analyzed(
                    cx,
                    &images,
                    0,
                    &spectral,
                    settings,
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(visible.volume_events(), &[0; 4]);
            assert_eq!(
                visible.free_flight_draws(),
                image_visible.free_flight_draws()
            );
            let beams = packet(
                &physical,
                cx,
                visible.pixel_sample().wavelengths_nm(),
                Some(AnalyzerResponse::Spectral(&spectral)),
            );
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                let expected = expected_point(&source, beam, rate);
                near(
                    visible.pixel_sample().spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                near(
                    image_visible.pixel_sample().spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                let manual = images
                    .sample_sequential_ray_with_tracking_samples(
                        &render_cx,
                        &beam,
                        0,
                        settings,
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
                assert_eq!(
                    image_visible.pixel_sample().spectral_samples()[lane],
                    Some(*manual.radiance_sample())
                );
                assert_eq!(
                    image_visible.point_light_sample().trace_evaluations()[lane],
                    u64::from(manual.trace_evaluations())
                );
                assert_eq!(
                    image_visible.point_light_sample().query_work()[lane],
                    manual.query_work()
                );
                assert_eq!(
                    image_visible.free_flight_draws()[lane],
                    manual.free_flight_draws()
                );
                assert_eq!(
                    image_visible.source_free_flight_draws()[lane],
                    manual.source_free_flight_draws()
                );
                assert_eq!(
                    image_visible.source_candidates()[lane],
                    manual.source_candidates()
                );
            }
            let detected = raster
                .sample_pixel_sensor_sequential_rayleigh_point_lights_spectrally_analyzed(
                    cx, &lit, &sensor, 0, &spectral, settings, SEED, 0, SAMPLE,
                )
                .unwrap();
            let image_detected = raster
                .sample_pixel_sensor_sequential_rayleigh_image_lights_spectrally_analyzed(
                    cx,
                    &images,
                    &sensor,
                    0,
                    &spectral,
                    settings,
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let wavelengths = detected
                .pixel_sample()
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(
                &physical,
                cx,
                &wavelengths,
                Some(AnalyzerResponse::Spectral(&spectral)),
            );
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            for lane in 0..4 {
                let expected = expected_point(&source, beams[lane].unwrap(), rate);
                near(
                    detected.pixel_sample().spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                near(
                    image_detected.pixel_sample().spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    expected,
                );
                for channel in 0..3 {
                    let contribution = expected
                        * beams[lane].unwrap().geometry_weight()
                        * [0.25, 0.75, 0.0][channel]
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0;
                    energy[channel] += contribution;
                    electrons[channel] +=
                        contribution * wavelengths[lane] * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
                }
            }
            for channel in 0..3 {
                near(
                    detected.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    image_detected.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    image_detected.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                    electrons[channel],
                );
            }
        }

        let closed_camera = coated_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let dark = closed
            .sample_pixel_sensor_sequential_point_lights(
                cx,
                &lights,
                &sensor,
                0,
                point_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(dark.connections(), &[0; 4]);
        assert_eq!(dark.trace_evaluations(), &[0; 4]);
        let incomplete = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (500.0, Stokes::unpolarized(0.5)),
                (600.0, Stokes::unpolarized(0.5)),
            ],
        )
        .unwrap();
        assert!(
            closed
                .sample_pixel_sensor_sequential_point_lights_spectrally_analyzed(
                    cx,
                    &lights,
                    &sensor,
                    0,
                    &incomplete,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_point_lights(
                    cx,
                    &lights,
                    &sensor,
                    0,
                    point_settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}
