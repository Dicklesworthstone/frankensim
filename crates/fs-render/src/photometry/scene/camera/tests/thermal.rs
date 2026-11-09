//! G1/G3: original Kirchhoff surfaces observed through native scene transport.

use super::*;
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::gradient_index::event::BoundarySide;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::DiffuseThermalEmitter;
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::{
    AreaLightBinding, AreaLightCameraSettings, AreaLightCameraWorkspace, CurvedAreaLightCameraScene,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use crate::thermal_surface::{
    DiffuseThermalSurface, GrayDiffuseThermalSurface, SpectralDiffuseThermalSurface,
    SpectralEmissivity, ThermalSurfaceError,
};

fn gray(temperature: f64, epsilon: f64) -> GrayDiffuseThermalSurface {
    GrayDiffuseThermalSurface::try_new(
        temperature,
        epsilon,
        hash_domain("native-original-gray-source", &temperature.to_le_bytes()),
    )
    .unwrap()
}

fn spectral(temperature: f64, epsilon: f64) -> SpectralDiffuseThermalSurface {
    SpectralDiffuseThermalSurface::try_new(
        temperature,
        SpectralEmissivity::constant_in_band(epsilon, (8_000.0, 14_000.0)).unwrap(),
        hash_domain(
            "native-original-spectral-source",
            &temperature.to_le_bytes(),
        ),
    )
    .unwrap()
}

fn thermal_panel(z: f64, width: f64, material: Material) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-width, -width, z],
                [width, -width, z],
                [0.0, 2.0 * width, z],
            ],
            vec![[0, 1, 2]],
        )),
        material,
        emission: None,
    }
}

fn infrared() -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(7_000.0 + 1_000.0 * i as f64, 1.0).unwrap()
            }),
            1_000.0,
        )
        .unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

fn blackbody(temperature: f64, wavelength_nm: f64) -> f64 {
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    let lambda = wavelength_nm * 1e-9;
    2.0 * h * c * c * 1e-9 / (lambda.powi(5) * (h * c / (lambda * k * temperature)).exp_m1())
}

#[test]
fn g1_original_gray_source_refracts_with_actual_beer_and_unscaled_si_radiance() {
    with_cx(|cx, _| {
        let source = gray(1_800.0, 1.0);
        let owner = admit(
            cx,
            scene(vec![
                Primitive {
                    shape: Shape::Chart(Box::new(SphereChart {
                        center: Point3::new(0.0, 0.0, 0.0),
                        radius: 1.0,
                    })),
                    material: Material::Dielectric {
                        glass: glass(1.5, 0.2),
                        surface: DielectricSurface::SMOOTH,
                    },
                    emission: None,
                },
                thermal_panel(2.0, 4.0, Material::GrayDiffuseThermal { surface: source }),
            ]),
            glass(1.0, 0.1),
        );
        let inside = owner
            .medium_index(NativeOpticalMedium::InteriorOf(0))
            .unwrap();
        for wavelength in [401.25, 550.0, 699.125] {
            for (origin, medium, depth, index) in [
                (Point3::new(0.0, 0.0, -3.0), 0, 0.7_f64, 1.0_f64),
                (Point3::new(0.0, 0.0, 0.0), inside, 0.3, 1.5),
            ] {
                let mut prior = None;
                for authored_scale in [0.0, 1e9] {
                    let sample = owner
                        .with_curved_camera(
                            cx,
                            100_000,
                            NativeRadianceCalibration::try_new(authored_scale, authored_scale)
                                .unwrap(),
                            |native, camera| {
                                let area =
                                    CurvedAreaEmitter::from_camera(cx, native, 2, 1).unwrap();
                                let thermal = area.source_thermal_emitter().unwrap();
                                assert!(area.source_emitter().is_none());
                                assert!(matches!(
                                    thermal.surface(),
                                    DiffuseThermalSurface::Gray(_)
                                ));
                                assert_eq!(thermal.surface().identity(), source.identity());
                                assert!(core::ptr::eq(
                                    area.source_instance(),
                                    &owner.geometry().instances()[1]
                                ));
                                native
                                    .sample_ray_with_samples(
                                        cx,
                                        Ray {
                                            origin,
                                            dir: camera.forward(),
                                        },
                                        wavelength,
                                        medium,
                                        settings(),
                                        |_| [0.0, 0.0, 0.99],
                                    )
                                    .unwrap()
                            },
                        )
                        .unwrap();
                    // Conditional smooth transmission retains its native
                    // Fresnel/PDF cancellation; actual Beer lengths and the
                    // camera/source index ratio still belong to transport.
                    close(
                        sample.radiance_w_m2_sr_nm(),
                        blackbody(1_800.0, wavelength) * (-depth).exp() * index * index,
                        2e-7,
                    );
                    close(sample.path().optical_depth(), depth, 2e-7);
                    if let Some(original) = prior {
                        assert_eq!(
                            sample, original,
                            "authored light calibration cannot scale heat"
                        );
                    } else {
                        prior = Some(sample);
                    }
                }
            }
        }
    });
}

#[test]
fn g1_original_gray_and_spectral_sources_share_kirchhoff_and_complementary_area_mis() {
    with_cx(|cx, _| {
        let epsilon = 0.37;
        let front = gray(300.0, epsilon);
        let back = spectral(300.0, 1.0);
        let owner = admit(
            cx,
            scene(vec![
                thermal_panel(0.0, 4.0, Material::GrayDiffuseThermal { surface: front }),
                thermal_panel(
                    -4.0,
                    0.3,
                    Material::SpectralDiffuseThermal { surface: back },
                ),
            ]),
            infrared(),
        );
        owner
            .with_curved_camera(
                cx,
                100_000,
                NativeRadianceCalibration::try_new(0.0, 0.0).unwrap(),
                |native, camera| {
                    let light = CurvedAreaEmitter::from_camera(cx, native, 2, 1).unwrap();
                    assert_eq!(
                        light.source_thermal_emitter().unwrap().surface().identity(),
                        back.identity()
                    );
                    let front_area = CurvedAreaEmitter::from_camera(cx, native, 1, 1).unwrap();
                    let gray_source = front_area.source_thermal_emitter().unwrap();
                    assert_eq!(gray_source.surface().identity(), front.identity());
                    assert_eq!(
                        gray_source
                            .reflectance_spectrum()
                            .unwrap()
                            .wavelength_bounds_nm(),
                        infrared().ior().wavelength_bounds_nm()
                    );
                    let seeds = [BundleImageSeed {
                        initial_direction: Vec3::new(0.0, 0.0, 1.0),
                        policy: BundleBranchPolicy::Sequence(&[]),
                    }];
                    let bindings = [AreaLightBinding {
                        source: &light,
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
                    let lit = CurvedAreaLightCameraScene::try_new(cx, native, &bindings).unwrap();
                    let mut work = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
                    let mut connection = BundleConnectionOptions::default();
                    connection.scene.trace = settings();
                    let limits = AreaLightCameraSettings {
                        surfaces: settings(),
                        connection,
                        max_area_samples: 1,
                        max_trace_evaluations: 4096,
                    };
                    let ray = Ray {
                        origin: camera.eye(),
                        dir: camera.forward(),
                    };
                    for wavelength in [8_000.0, 9_750.0, 12_250.0, 14_000.0] {
                        let baseline = native
                            .sample_ray_with_samples(cx, ray, wavelength, 0, settings(), |_| {
                                [0.0; 3]
                            })
                            .unwrap();
                        let b = blackbody(300.0, wavelength);
                        // epsilon*B at the front plus (1-epsilon)*B from the
                        // other original body closes the equilibrium radiance.
                        close(baseline.radiance_w_m2_sr_nm(), b, 2e-11);
                        let sample = lit
                            .sample_ray_with_samples(
                                cx,
                                ray,
                                wavelength,
                                0,
                                limits,
                                &mut work,
                                |_| [0.0; 3],
                                |event, binding| {
                                    assert_eq!((event, binding), (0, 0));
                                    [4.0 / 9.0, 0.5]
                                },
                            )
                            .unwrap();
                        let area = light.total_area_m2();
                        let p_light = 16.0 / area;
                        let p_bsdf = 1.0 / core::f64::consts::PI;
                        let denominator = p_light * p_light + p_bsdf * p_bsdf;
                        let reflected = (1.0 - epsilon) * b;
                        let nee =
                            reflected * area / (16.0 * core::f64::consts::PI) * p_light * p_light
                                / denominator;
                        close(sample.area_radiance_w_m2_sr_nm(), nee, 4e-6);
                        close(
                            sample.radiance_sample().radiance_w_m2_sr_nm(),
                            epsilon * b + reflected * p_bsdf * p_bsdf / denominator + nee,
                            4e-6,
                        );
                        assert_eq!(sample.radiance_sample().path(), baseline.path());
                        assert_eq!((sample.area_samples(), sample.images()), (1, 1));
                    }
                    // The original medium supports 7–15 um. The original spectral
                    // source's narrower 8–14 um law must still refuse at 7.5 um.
                    assert!(matches!(
                        native
                            .sample_ray_with_samples(cx, ray, 7_500.0, 0, settings(), |_| [0.0; 3]),
                        Err(CameraRadianceError::Thermal(
                            ThermalSurfaceError::Wavelength
                        ))
                    ));
                },
            )
            .unwrap();
    });
}

#[test]
fn g1_original_gray_infrared_sensor_keeps_photon_weights_two_sides_and_temperature() {
    with_cx(|cx, _| {
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(8_000.0, [0.0, 1.0]), (14_000.0, [1.0, 0.0])])
                .unwrap();
        let mut cold: Option<[f64; 2]> = None;
        for temperature in [300.0, 450.0] {
            let source = gray(temperature, 0.42);
            let owner = admit(
                cx,
                scene(vec![thermal_panel(
                    0.0,
                    4.0,
                    Material::GrayDiffuseThermal { surface: source },
                )]),
                infrared(),
            );
            let pixel = owner
                .with_curved_camera(
                    cx,
                    100_000,
                    NativeRadianceCalibration::try_new(0.0, 0.0).unwrap(),
                    |native, camera| {
                        for side in [-1.0, 1.0] {
                            let sample = native
                                .sample_ray_with_samples(
                                    cx,
                                    Ray {
                                        origin: Point3::new(0.0, 0.0, side),
                                        dir: Vec3::new(0.0, 0.0, -side),
                                    },
                                    11_000.0,
                                    0,
                                    settings(),
                                    |_| [0.0; 3],
                                )
                                .unwrap();
                            close(
                                sample.radiance_w_m2_sr_nm(),
                                0.42 * blackbody(temperature, 11_000.0),
                                2e-11,
                            );
                        }
                        let raster = CurvedCameraRaster::try_new(camera, 1, 1).unwrap();
                        let pixel = raster
                            .sample_pixel_sensor(cx, native, &sensor, 0, settings(), 44, 0, 3)
                            .unwrap();
                        assert_eq!(
                            pixel,
                            raster
                                .sample_pixel_sensor(cx, native, &sensor, 0, settings(), 44, 0, 3)
                                .unwrap()
                        );
                        let (mut power, mut electrons) = ([0.0; 2], [0.0; 2]);
                        for (draw, lane) in pixel
                            .detector_wavelengths()
                            .iter()
                            .zip(pixel.spectral_samples())
                        {
                            let wavelength = draw.wavelength_nm();
                            let radiance = 0.42 * blackbody(temperature, wavelength);
                            close(lane.radiance_w_m2_sr_nm(), radiance, 2e-11);
                            let fraction = (wavelength - 8_000.0) / 6_000.0;
                            for (channel, qe) in [fraction, 1.0 - fraction].into_iter().enumerate()
                            {
                                let weighted = radiance * qe * 6_000.0 / 4.0;
                                power[channel] += weighted;
                                electrons[channel] += weighted * wavelength * 1e-9
                                    / (6.626_070_15e-34 * 299_792_458.0);
                            }
                        }
                        for channel in 0..2 {
                            close(
                                pixel.qe_weighted_radiance_w_m2_sr()[channel],
                                power[channel],
                                2e-11,
                            );
                            close(
                                pixel.photoelectron_rate_per_m2_sr()[channel],
                                electrons[channel],
                                2e-11,
                            );
                        }
                        *pixel.photoelectron_rate_per_m2_sr()
                    },
                )
                .unwrap();
            if let Some(previous) = cold {
                for channel in 0..2 {
                    assert!(pixel[channel] > previous[channel]);
                }
            } else {
                cold = Some(pixel);
            }
        }
    });
}

#[test]
fn g0_original_thermal_sources_refuse_duplicate_emission_wrong_medium_and_material() {
    let source = gray(1_800.0, 0.4);
    for bounds in [
        (0.0, 780.0),
        (780.0, 380.0),
        (550.0, 550.0),
        (380.0, f64::NAN),
    ] {
        assert!(DiffuseThermalEmitter::try_gray_in_band(1, source, bounds).is_err());
    }
    with_cx(|cx, gate| {
        for material in [
            Material::GrayDiffuseThermal { surface: source },
            Material::SpectralDiffuseThermal {
                surface: spectral(300.0, 0.4),
            },
        ] {
            let mut duplicate = thermal_panel(0.0, 4.0, material);
            duplicate.emission = Some((SPECTRUM, 0.0));
            assert!(matches!(
                NativeOpticalScene::try_from_scene(
                    cx,
                    scene(vec![duplicate]),
                    glass(1.0, 0.0),
                    &[binding(1)],
                    1
                ),
                Err(NativeOpticalSceneError::ConflictingThermalEmission { primitive: 0 })
            ));
        }
        let owner = admit(
            cx,
            scene(vec![thermal_panel(
                0.0,
                4.0,
                Material::GrayDiffuseThermal { surface: source },
            )]),
            glass(1.3, 0.0),
        );
        let calibration = NativeRadianceCalibration::try_new(0.0, 0.0).unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, |native, camera| {
                assert!(matches!(
                    native.sample_ray_with_samples(
                        cx,
                        Ray {
                            origin: camera.eye(),
                            dir: camera.forward()
                        },
                        550.0,
                        0,
                        settings(),
                        |_| [0.0; 3]
                    ),
                    Err(CameraRadianceError::UnsupportedSurface { object_id: 1 })
                ));
                let wrong =
                    DiffuseThermalEmitter::try_gray_in_band(1, gray(1_800.0, 0.5), (380.0, 780.0))
                        .unwrap();
                let unbound =
                    CurvedCameraScene::try_new(cx, native.source_scene(), &[], None).unwrap();
                assert!(matches!(
                    unbound.with_diffuse_thermal_emitters(cx, &[wrong]),
                    Err(CameraRadianceError::InvalidEmitter { object_id: 1 })
                ));
            })
            .unwrap();
        assert!(
            owner
                .with_curved_camera(cx, 0, calibration, |_, _| panic!("no partial view"))
                .is_err()
        );
        gate.request();
        assert!(matches!(
            owner.with_curved_camera(cx, 100_000, calibration, |_, _| ()),
            Err(NativeCameraError::OpticalScene(
                NativeOpticalSceneError::Cancelled
            ))
        ));
    });
}
