use super::*;
use crate::charts::Ray;
use crate::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::camera::CurvedCameraScene;
use crate::polarization::dielectric::scene::curved::camera::point_lights::{
    CurvedPointLightCameraScene, PointLightBranch, PointLightCameraSettings,
};
use crate::spectral::SpectralReflectance;
use core::f64::consts::PI;

fn opaque(object_id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Opaque,
    }
}

fn lamp(z: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, z]).unwrap(),
        100.0,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
    )
    .unwrap()
}

fn intensity(lamp: &PhotometricSource) -> f64 {
    lamp.flux_lm() * lamp.spectrum().bins_w_per_lumen_nm()[0] / (4.0 * PI)
}

fn light_settings() -> PointLightCameraSettings {
    PointLightCameraSettings {
        surfaces: curved_settings(),
        connection: BundleConnectionOptions {
            scene: CurvedBundleSettings {
                trace: curved_settings(),
                ..CurvedBundleSettings::default()
            },
            ..BundleConnectionOptions::default()
        },
        max_connections: 16,
        max_trace_evaluations: 256,
    }
}

#[test]
fn g1_calibrated_rough_stack_receivers_keep_original_spectra_and_both_media() {
    with_cx(|cx| {
        let (transparent, absorbing) = films();
        let layers = [
            CoatingLayer::Dielectric(&transparent),
            CoatingLayer::Absorbing(&absorbing),
        ];
        let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, true),
            plane(2, 2.0, true),
            plane(3, -2.0, true),
        ])
        .unwrap();
        let indices = [1.0_f64, 1.6];
        let absorption = [0.1_f64, 0.3];
        let materials = [
            glass(indices[0], absorption[0]),
            glass(indices[1], absorption[1]),
        ];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let alpha = 0.2_f64;
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::try_rough(alpha).unwrap(),
                    },
                },
                opaque(2),
                opaque(3),
            ],
            100_000,
        )
        .unwrap()
        .with_dielectric_coating_stacks(
            cx,
            &[DielectricStackBinding {
                object_id: 1,
                stack,
            }],
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        for (wavelength, r_out, r_in, transmission) in POWERS {
            let mut reciprocal = [0.0; 2];
            for camera_medium in 0..2 {
                let camera_side = if camera_medium == 0 { 1.0 } else { -1.0 };
                let ray = Ray {
                    origin: Point3::new(0.0, 0.0, 0.5 * camera_side),
                    dir: Vec3::new(0.0, 0.0, -camera_side),
                };
                let ordinary = camera
                    .sample_ray_with_samples(
                        cx,
                        ray,
                        wavelength,
                        camera_medium,
                        curved_settings(),
                        |_| [0.0, 0.0, 0.99],
                    )
                    .unwrap();
                for source_medium in 0..2 {
                    let source_side = if source_medium == 0 { 1.0 } else { -1.0 };
                    let source = lamp(source_side);
                    let branches = [PointLightBranch {
                        source: &source,
                        receiver_object_id: 1,
                        initial_medium: source_medium,
                        initial_direction: Vec3::new(0.0, 0.0, -source_side),
                        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                        policy: BundleBranchPolicy::Sequence(&[]),
                        position_tolerance_m: 2e-6,
                    }];
                    let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
                    let sample = lit
                        .sample_ray_with_samples(
                            cx,
                            ray,
                            wavelength,
                            camera_medium,
                            light_settings(),
                            |_| [0.0, 0.0, 0.99],
                        )
                        .unwrap();
                    // Independent normal-incidence GGX radiance law with the
                    // original lossy Airy powers, never 1-R for transmission.
                    let response = if source_medium == camera_medium {
                        [r_out, r_in][camera_medium] / (4.0 * PI * alpha * alpha)
                    } else {
                        transmission
                            / (PI
                                * alpha
                                * alpha
                                * (1.0 - indices[source_medium] / indices[camera_medium]).powi(2))
                    };
                    let arriving = intensity(&source)
                        * (-absorption[source_medium] - 0.5 * absorption[camera_medium]).exp();
                    close(
                        sample.fixture_radiance_w_m2_sr_nm(),
                        arriving * response,
                        5e-6,
                    );
                    close(
                        sample.radiance_sample().radiance_w_m2_sr_nm(),
                        arriving * response,
                        5e-6,
                    );
                    assert_eq!(sample.radiance_sample().path(), ordinary.path());
                    assert_eq!(sample.connections(), 1);
                    if camera_medium != source_medium {
                        reciprocal[camera_medium] = sample.fixture_radiance_w_m2_sr_nm()
                            / arriving
                            / indices[camera_medium].powi(2);
                    }
                }
            }
            close(reciprocal[0], reciprocal[1], 5e-6);
        }
    });
}

#[test]
fn g1_original_stack_source_connection_retains_snell_area_and_physical_transmission() {
    with_cx(|cx| {
        let (transparent, absorbing) = films();
        let layers = [
            CoatingLayer::Dielectric(&transparent),
            CoatingLayer::Absorbing(&absorbing),
        ];
        let stack = DielectricCoatingStack::try_mixed(&layers).unwrap();
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.25, false),
            plane(2, 1.0, true),
            plane(3, -0.5, true),
        ])
        .unwrap();
        let materials = [glass(1.0, 0.1), glass(1.6, 0.3)];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let reflectances = [SpectralReflectance::constant(0.5).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                boundary(1, 0, 1),
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 1,
                        reflectance: 0,
                    },
                },
                opaque(3),
            ],
            100_000,
        )
        .unwrap()
        .with_dielectric_coating_stacks(
            cx,
            &[DielectricStackBinding {
                object_id: 1,
                stack,
            }],
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let source = lamp(0.0);
        let events = [DielectricEvent::Transmission];
        let branches = [PointLightBranch {
            source: &source,
            receiver_object_id: 2,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::Sequence(&events),
            position_tolerance_m: 2e-6,
        }];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        for (wavelength, _, _, transmission) in POWERS {
            let sample = lit
                .sample_ray_with_samples(
                    cx,
                    Ray {
                        origin: Point3::new(0.0, 0.0, 0.75),
                        dir: Vec3::new(0.0, 0.0, 1.0),
                    },
                    wavelength,
                    1,
                    light_settings(),
                    |_| [0.0, 0.0, 0.99],
                )
                .unwrap();
            // Launch-angle derivatives use the two bulk indices. Coating
            // optical thickness changes T, not this geometric source footprint.
            let area = (0.25_f64 + 0.75 / 1.6).powi(2);
            let expected = intensity(&source) * transmission / area * 0.5 / PI
                * (-0.1_f64 * 0.25 - 0.3 * (0.75 + 0.25)).exp();
            close(sample.fixture_radiance_w_m2_sr_nm(), expected, 8e-6);
            close(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                expected,
                8e-6,
            );
            assert_eq!(sample.connections(), 1);
        }
    });
}
