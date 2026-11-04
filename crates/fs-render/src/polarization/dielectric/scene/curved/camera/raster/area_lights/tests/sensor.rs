//! Source-band detector integration through the original area-light walk.

use super::*;
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::spectral::detector::QuantumEfficiencySensor;

fn infrared_medium() -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.0).unwrap()
            }),
            750.0,
        )
        .unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

#[test]
fn g1_infrared_area_detector_keeps_nee_mis_and_physical_channel_weights() {
    fixtures::with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let glass = infrared_medium();
        let media = [SceneMedium::Homogeneous(&glass)];
        let reflectances = [SpectralReflectance::constant_in_band(0.6, (8000.0, 14000.0)).unwrap()];
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
            RadianceSpectrum::try_new_in_band(vec![(8000.0, 1.0), (14000.0, 5.0)]).unwrap(),
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
        let physical = fixtures::camera(0.01);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let (seed, pixel, sample) = (0x1234_5678_0000_0049_u64, 4, 9);
        let mut controls = settings();
        controls.surfaces.max_distance_m = 12.0;
        controls.connection.scene.trace = controls.surfaces;
        let result = raster
            .sample_pixel_sensor_area_lights(
                cx,
                &area,
                &sensor,
                0,
                controls,
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let observed = result.pixel_sample();
        let key = [seed as u32, (seed >> 32) as u32];
        let render_cx = cx.with_stream_seed(seed);
        let mut power = [0.0; 3];
        let mut electrons = [0.0; 3];
        let mut source_distances = [0.0; CAMERA_WAVELENGTHS];
        for lane in 0..CAMERA_WAVELENGTHS {
            let wavelength = observed.detector_wavelengths()[lane].wavelength_nm();
            assert!((8000.0..=14000.0).contains(&wavelength));
            let ray = observed.camera_rays()[lane];
            let camera_length = (1.0 - ray.origin.z) / ray.dir.z;
            let receiver = Point3::new(
                ray.origin.x + camera_length * ray.dir.x,
                ray.origin.y + camera_length * ray.dir.y,
                1.0,
            );
            let draw = philox4x32_10([pixel, sample, 0, 0x6372_6100 + lane as u32 * 64], key);
            let root = unit(draw[0]).sqrt();
            let v = unit(draw[1]);
            let point = Point3::new(
                -0.3 + 0.6 * root * (1.0 - v) + 0.3 * root * v,
                -0.3 + 0.9 * root * v,
                -1.0,
            );
            let delta = receiver.delta_from(point);
            let r2 = delta.dot(delta);
            let distance = r2.sqrt();
            source_distances[lane] = distance;
            let cosine = 2.0 / distance;
            let le = 1.0 + 4.0 * (wavelength - 8000.0) / 6000.0;
            let p_light = r2 / (0.27 * cosine);
            let p_bsdf = cosine / core::f64::consts::PI;
            let mis = p_light.powi(2) / (p_light.powi(2) + p_bsdf.powi(2));
            let expected = le * 0.27 * cosine.powi(2) / r2 * 0.6 / core::f64::consts::PI
                * (-0.2 * (camera_length + distance)).exp()
                * mis;
            near(result.area_radiance_w_m2_sr_nm()[lane], expected);
            assert_eq!(result.area_samples()[lane], 1);
            assert_eq!(result.images()[lane], 1);
            let manual = area
                .sample_ray_with_samples(
                    &render_cx,
                    ray,
                    wavelength,
                    0,
                    controls,
                    &mut workspace,
                    |event| {
                        let block =
                            philox4x32_10([pixel, sample, event, 0x6372_7340 + lane as u32], key);
                        core::array::from_fn(|i| unit(block[i]))
                    },
                    |event, binding| {
                        let draw = philox4x32_10(
                            [
                                pixel,
                                sample,
                                event,
                                0x6372_6100 + lane as u32 * 64 + binding as u32,
                            ],
                            key,
                        );
                        [unit(draw[0]), unit(draw[1])]
                    },
                )
                .unwrap();
            assert_eq!(&observed.spectral_samples()[lane], manual.radiance_sample());
            let eta = [
                (wavelength - 8000.0) / 6000.0,
                (14000.0 - wavelength) / 6000.0,
                0.0,
            ];
            for c in 0..3 {
                // Sum-QE is one: p(lambda)=1/6000nm. Use the full native
                // radiance, including the complementary emitter-hit term.
                let energy = manual.radiance_sample().radiance_w_m2_sr_nm() * eta[c] * 1500.0;
                power[c] += energy;
                electrons[c] += energy * wavelength * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for c in 0..3 {
            near(observed.qe_weighted_radiance_w_m2_sr()[c], power[c]);
            near(observed.photoelectron_rate_per_m2_sr()[c], electrons[c]);
        }
        assert_eq!(observed.photoelectron_rate_per_m2_sr()[2], 0.0);
        assert_eq!(
            result,
            raster
                .sample_pixel_sensor_area_lights(
                    cx,
                    &area,
                    &sensor,
                    0,
                    controls,
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .unwrap()
        );
        let empty = CurvedAreaLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let ordinary = raster
            .sample_pixel_sensor(
                cx,
                &camera,
                &sensor,
                0,
                controls.surfaces,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let unconnected = raster
            .sample_pixel_sensor_area_lights(
                cx,
                &empty,
                &sensor,
                0,
                controls,
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        assert_eq!(unconnected.pixel_sample(), &ordinary);
        assert_eq!(unconnected.area_samples(), &[0; CAMERA_WAVELENGTHS]);

        #[cfg(feature = "volumes")]
        {
            use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
                CurvedRayleighAreaLightCameraScene, RayleighAreaLightCameraSettings,
                RayleighAreaLightCameraWorkspace,
            };
            use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::participating::{
                CurvedRayleighMedium, CurvedRayleighScene,
            };
            use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
            let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            let phase = PhaseFunction::Rayleigh;
            let settings = RayleighAreaLightCameraSettings {
                area_lights: controls,
                max_volume_events: 8,
                max_source_free_flight_draws: 16,
                max_source_candidates: 16,
            };
            for coefficient in [0.0, 0.01] {
                let spectrum = MassScatteringSpectrum::source_table(
                    &[(8000.0, coefficient), (14000.0, 2.0 * coefficient)],
                    hash_domain("area-detector-scattering", b"infrared"),
                )
                .unwrap();
                let volumes = [Some(
                    CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
                )];
                let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
                let molecular =
                    CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
                let scene =
                    CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &area).unwrap();
                let scattered = raster
                    .sample_pixel_sensor_rayleigh_area_lights(
                        cx,
                        &scene,
                        &sensor,
                        0,
                        settings,
                        &mut workspace,
                        seed,
                        pixel,
                        sample,
                    )
                    .unwrap();
                assert_eq!(scattered.volume_events(), &[0; CAMERA_WAVELENGTHS]);
                for lane in 0..CAMERA_WAVELENGTHS {
                    let lambda = observed.detector_wavelengths()[lane].wavelength_nm();
                    let extinction = coefficient * (1.0 + (lambda - 8000.0) / 6000.0);
                    near(
                        scattered.area_pixel_sample().area_radiance_w_m2_sr_nm()[lane],
                        result.area_radiance_w_m2_sr_nm()[lane]
                            * (-extinction * source_distances[lane]).exp(),
                    );
                }
                if coefficient == 0.0 {
                    assert_eq!(
                        scattered.pixel_sample().camera_rays(),
                        observed.camera_rays()
                    );
                    assert_eq!(
                        scattered.pixel_sample().detector_wavelengths(),
                        observed.detector_wavelengths()
                    );
                    for c in 0..3 {
                        near(
                            scattered.pixel_sample().qe_weighted_radiance_w_m2_sr()[c],
                            power[c],
                        );
                        near(
                            scattered.pixel_sample().photoelectron_rate_per_m2_sr()[c],
                            electrons[c],
                        );
                    }
                }
                assert_eq!(
                    scattered,
                    raster
                        .sample_pixel_sensor_rayleigh_area_lights(
                            cx,
                            &scene,
                            &sensor,
                            0,
                            settings,
                            &mut workspace,
                            seed,
                            pixel,
                            sample,
                        )
                        .unwrap()
                );
            }
        }

        let mut limited = controls;
        limited.max_trace_evaluations = 1;
        assert!(matches!(
            raster.sample_pixel_sensor_area_lights(
                cx,
                &area,
                &sensor,
                0,
                limited,
                &mut workspace,
                seed,
                pixel,
                sample,
            ),
            Err(AreaLightRasterError::Transport { lane: 0, .. })
        ));
        // Reusing the worker after a refused sample retains absolute identity.
        assert_eq!(
            result,
            raster
                .sample_pixel_sensor_area_lights(
                    cx,
                    &area,
                    &sensor,
                    0,
                    controls,
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .unwrap()
        );
        gate.request();
        assert!(
            raster
                .sample_pixel_sensor_area_lights(
                    cx,
                    &area,
                    &sensor,
                    0,
                    controls,
                    &mut workspace,
                    seed,
                    pixel,
                    sample,
                )
                .is_err()
        );
    });
}
