//! Native thermal NEE and emitter hits use the same accepted spectral law.

use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::polarization::dielectric::scene::curved::camera::DiffuseThermalEmitter;
use crate::polarization::dielectric::scene::curved::camera::raster::{
    CurvedCameraRaster, sensor::CurvedSensorPixelSample,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::thermal_surface::{SpectralDiffuseThermalSurface, SpectralEmissivity};
use fs_blake3::ContentHash;

fn infrared_glass() -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(7000.0 + 1000.0 * i as f64, 1.0).unwrap()
            }),
            1000.0,
        )
        .unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn thermal(temperature: f64, identity: u8) -> DiffuseThermalEmitter {
    DiffuseThermalEmitter::new(
        9,
        SpectralDiffuseThermalSurface::try_new(
            temperature,
            SpectralEmissivity::table_in_band(&[(8000.0, 0.25), (11000.0, 0.6), (14000.0, 0.9)])
                .unwrap(),
            ContentHash([identity; 32]),
        )
        .unwrap(),
    )
}

fn planck(temperature: f64, wavelength: f64) -> f64 {
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    let lambda = wavelength * 1e-9;
    2.0 * h * c * c * 1e-9 / (lambda.powi(5) * (h * c / (lambda * k * temperature)).exp_m1())
}

fn continuation(event: u32) -> [f64; 3] {
    // Hit the source at normal incidence, then reflect outside the finite
    // receiver. Emission is added before that source's native reflection.
    if event == 0 {
        [0.0; 3]
    } else {
        [0.96, 0.5, 0.0]
    }
}

#[test]
fn g1_infrared_planck_nee_and_actual_emitter_hits_share_complementary_native_mis() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let material = infrared_glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(1e-4).unwrap(),
            1.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [0.25, 0.75]), (14000.0, [0.75, 0.25])])
                .unwrap();
        let mut cold_pixel: Option<CurvedSensorPixelSample<2>> = None;
        for temperature in [300.0, 450.0] {
            let emitter = thermal(temperature, 51);
            let reflectances = [
                SpectralReflectance::table_in_band(&[(7000.0, RHO), (15000.0, RHO)]).unwrap(),
                emitter.reflectance_spectrum().unwrap(),
            ];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![
                    diffuse(1, 0),
                    SurfaceBinding {
                        object_id: 9,
                        interaction: Interaction::DepolarizingDiffuse {
                            medium: 0,
                            reflectance: 1,
                        },
                    },
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
                .unwrap()
                .with_diffuse_thermal_emitters(cx, &[emitter])
                .unwrap();
            let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
            assert!(light.source_emitter().is_none());
            assert_eq!(light.source_thermal_emitter(), Some(&emitter));
            let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
            let bindings = [AreaLightBinding {
                source: &light,
                receiver_object_id: 1,
                source_side: BoundarySide::Positive,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
            let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
            for wavelength in [8000.0, 9750.0, 12250.0, 14000.0] {
                let baseline = camera
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        wavelength,
                        0,
                        settings().surfaces,
                        continuation,
                    )
                    .unwrap();
                let result = owner
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        wavelength,
                        0,
                        settings(),
                        &mut workspace,
                        continuation,
                        |event, binding| {
                            assert_eq!((event, binding), (0, 0));
                            [4.0 / 9.0, 0.5]
                        },
                    )
                    .unwrap();
                let epsilon = if wavelength <= 11000.0 {
                    0.25 + (wavelength - 8000.0) * 0.35 / 3000.0
                } else {
                    0.6 + (wavelength - 11000.0) * 0.3 / 3000.0
                };
                let le = epsilon * planck(temperature, wavelength);
                let p_light: f64 = 4.0 / AREA;
                let p_bsdf = 1.0 / core::f64::consts::PI;
                let denominator = p_light * p_light + p_bsdf * p_bsdf;
                let raw_hit = RHO * le * (-0.6_f64).exp();
                let nee = raw_hit * AREA / (4.0 * core::f64::consts::PI) * p_light * p_light
                    / denominator;
                let hit = raw_hit * p_bsdf * p_bsdf / denominator;
                near(baseline.radiance_w_m2_sr_nm(), raw_hit, 3e-6);
                near(result.area_radiance_w_m2_sr_nm(), nee, 4e-6);
                near(
                    result.radiance_sample().radiance_w_m2_sr_nm(),
                    nee + hit,
                    4e-6,
                );
                assert_eq!(result.radiance_sample().path(), baseline.path());
                assert_eq!((result.area_samples(), result.images()), (1, 1));
                assert!(result.trace_evaluations() >= 3);

                #[cfg(feature = "volumes")]
                {
                    use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
                    use crate::polarization::dielectric::scene::curved::camera::area_lights::participating::{
                        CurvedRayleighAreaLightCameraScene, RayleighAreaLightCameraSettings,
                        RayleighAreaLightCameraWorkspace,
                    };
                    use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighScene;
                    let volumes = [None];
                    let scattering = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
                    let molecular =
                        CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
                    let participating =
                        CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &owner)
                            .unwrap();
                    let mut work = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
                    let replay = participating
                        .sample_ray_with_tracking_samples(
                            cx,
                            camera_ray(),
                            wavelength,
                            0,
                            RayleighAreaLightCameraSettings {
                                area_lights: settings(),
                                max_volume_events: 1,
                                max_source_free_flight_draws: 0,
                                max_source_candidates: 0,
                            },
                            &mut work,
                            continuation,
                            |_| panic!("no volume flight"),
                            |_, _| [4.0 / 9.0, 0.5],
                            |_, _| panic!("no volume vertex"),
                            |_, _, _| panic!("no source scattering"),
                        )
                        .unwrap();
                    assert_eq!(replay.radiance_sample(), result.radiance_sample());
                    assert_eq!(
                        replay.area_radiance_w_m2_sr_nm(),
                        result.area_radiance_w_m2_sr_nm()
                    );
                    assert_eq!(
                        (replay.volume_events(), replay.source_free_flight_draws()),
                        (0, 0)
                    );
                }
            }
            let detector = raster
                .sample_pixel_sensor_area_lights(
                    cx,
                    &owner,
                    &sensor,
                    0,
                    settings(),
                    &mut workspace,
                    73,
                    0,
                    4,
                )
                .unwrap();
            assert!(
                detector
                    .area_radiance_w_m2_sr_nm()
                    .iter()
                    .all(|&radiance| radiance > 0.0)
            );
            let pixel = detector.pixel_sample();
            if let Some(cold) = cold_pixel {
                assert_eq!(cold.detector_wavelengths(), pixel.detector_wavelengths());
                let mut energy = [0.0; 2];
                let mut electrons = [0.0; 2];
                for ((draw, cold_path), warm_path) in cold
                    .detector_wavelengths()
                    .iter()
                    .zip(cold.spectral_samples())
                    .zip(pixel.spectral_samples())
                {
                    let wavelength = draw.wavelength_nm();
                    let expected = cold_path.radiance_w_m2_sr_nm() * planck(450.0, wavelength)
                        / planck(300.0, wavelength);
                    near(warm_path.radiance_w_m2_sr_nm(), expected, 3e-11);
                    assert_eq!(cold_path.path(), warm_path.path());
                    for channel in 0..2 {
                        energy[channel] += expected * draw.radiance_weights_nm()[channel] / 4.0;
                        electrons[channel] +=
                            expected * draw.electron_weights_nm_per_joule()[channel] / 4.0;
                    }
                }
                for channel in 0..2 {
                    near(
                        pixel.qe_weighted_radiance_w_m2_sr()[channel],
                        energy[channel],
                        3e-11,
                    );
                    near(
                        pixel.photoelectron_rate_per_m2_sr()[channel],
                        electrons[channel],
                        3e-11,
                    );
                }
            } else {
                cold_pixel = Some(*pixel);
            }
            assert!(
                owner
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        7999.0,
                        0,
                        settings(),
                        &mut workspace,
                        continuation,
                        |_, _| [4.0 / 9.0, 0.5]
                    )
                    .is_err()
            );
            let exhausted = AreaLightCameraSettings {
                max_area_samples: 0,
                ..settings()
            };
            assert!(
                owner
                    .sample_ray_with_samples(
                        cx,
                        camera_ray(),
                        11000.0,
                        0,
                        exhausted,
                        &mut workspace,
                        continuation,
                        |_, _| [4.0 / 9.0, 0.5]
                    )
                    .is_err()
            );
            assert_eq!(emitter.surface().temperature_k(), temperature);
        }
    });
}

#[test]
fn g0_area_camera_refuses_stale_temperature_identity_or_calibrated_surrogate() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let material = infrared_glass();
        let media = [SceneMedium::Homogeneous(&material)];
        let emitter = thermal(300.0, 51);
        let reflectances = [
            SpectralReflectance::table_in_band(&[(7000.0, RHO), (15000.0, RHO)]).unwrap(),
            emitter.reflectance_spectrum().unwrap(),
        ];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                diffuse(1, 0),
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 1,
                    },
                },
            ],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[emitter])
            .unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        for candidate in [thermal(301.0, 51), thermal(300.0, 52)] {
            let source = CurvedAreaEmitter::try_new_thermal(cx, &scene, candidate, 1).unwrap();
            let bindings = [AreaLightBinding {
                source: &source,
                receiver_object_id: 1,
                source_side: BoundarySide::Positive,
                initial_medium: 0,
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                seeds: &seeds,
                image_options: image_options(1),
                position_tolerance_m: 2e-7,
            }];
            assert!(matches!(
                CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings),
                Err(AreaLightCameraError::InvalidBinding { index: 0 })
            ));
        }
        let fake = SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_new_in_band(vec![(8000.0, 1.0), (14000.0, 1.0)]).unwrap(),
            EmissionSides::Both,
        );
        let source = CurvedAreaEmitter::try_new(cx, &scene, &fake, 1).unwrap();
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        assert!(matches!(
            CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings),
            Err(AreaLightCameraError::InvalidBinding { index: 0 })
        ));
    });
}
