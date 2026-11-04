use super::*;
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::thermal_surface::{
    SpectralDiffuseThermalSurface, SpectralEmissivity, ThermalSurfaceError,
};
use fs_blake3::ContentHash;

fn infrared_medium(index: f64) -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(7000.0 + 1000.0 * i as f64, index).unwrap()
            }),
            1000.0,
        )
        .unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn thermal(epsilon: f64) -> DiffuseThermalEmitter {
    thermal_at(epsilon, 320.0)
}

fn thermal_at(epsilon: f64, temperature: f64) -> DiffuseThermalEmitter {
    DiffuseThermalEmitter::new(
        2,
        SpectralDiffuseThermalSurface::try_new(
            temperature,
            SpectralEmissivity::constant_in_band(epsilon, (8000.0, 14000.0)).unwrap(),
            ContentHash([61; 32]),
        )
        .unwrap(),
    )
}

fn planck(wavelength: f64) -> f64 {
    let (h, c, k) = (6.626_070_15e-34, 299_792_458.0, 1.380_649e-23);
    let lambda = wavelength * 1e-9;
    2.0 * h * c * c * 1e-9 / (lambda.powi(5) * (h * c / (lambda * k * 320.0)).exp_m1())
}

#[test]
fn g1_native_thermal_area_irradiance_keeps_per_face_planck_units_and_local_vacuum() {
    with_cx(|cx| {
        for direction in [1.0, -1.0] {
            let point = Point3::new(0.0, 0.0, 4.0 * direction);
            let geometry = InstanceScene::try_new(vec![source(0.0), detector(point)]).unwrap();
            for index in [1.0, 1.3] {
                let material = infrared_medium(index);
                let media = [SceneMedium::Homogeneous(&material)];
                let scene = CurvedPolarizedScene::try_new(
                    cx,
                    &geometry,
                    &media,
                    vec![
                        SurfaceBinding {
                            object_id: 2,
                            interaction: Interaction::Opaque,
                        },
                        detector_binding(),
                    ],
                    100_000,
                )
                .unwrap();
                let emitter = thermal(1.0);
                let light = CurvedAreaEmitter::try_new_thermal(cx, &scene, emitter, 1).unwrap();
                let sample = light.sample_area([4.0 / 9.0, 0.5]).unwrap();
                let side = if direction > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::Negative
                };
                let seeds = [BundleImageSeed {
                    initial_direction: Vec3::new(0.0, 0.0, direction),
                    policy: BundleBranchPolicy::TransmitOrReflectTir,
                }];
                let target = BundleConnectionTarget::try_new(
                    9,
                    point,
                    Vec3::new(0.0, 0.0, -direction),
                    2e-7,
                )
                .unwrap();
                let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
                let result = light.incident_irradiance(
                    cx,
                    sample,
                    side,
                    0,
                    Vec3::new(0.0, 1.0, 0.0),
                    9750.0,
                    target,
                    options(),
                    &seeds,
                    BundleImageOptions {
                        max_images: 1,
                        ..BundleImageOptions::default()
                    },
                    &mut workspace,
                );
                if index == 1.0 {
                    let result = result.unwrap();
                    near(
                        result.spectral_irradiance_w_m2_nm(),
                        planck(9750.0) * 0.27 / 16.0 * (-0.8_f64).exp(),
                        3e-6,
                    );
                    assert_eq!(result.image_count(), 1);
                    assert!(matches!(
                        light.launch(sample, side, 0, Vec3::new(0.0, 1.0, 0.0), 7500.0, &seeds),
                        Err(AreaEmitterError::Camera(CameraRadianceError::Thermal(
                            ThermalSurfaceError::Wavelength
                        )))
                    ));
                } else {
                    assert!(matches!(
                        result,
                        Err(AreaEmitterError::Camera(
                            CameraRadianceError::UnsupportedSurface { object_id: 2 }
                        ))
                    ));
                }
                assert!(matches!(
                    CurvedAreaEmitter::try_new_thermal(cx, &scene, thermal(0.5), 1),
                    Err(AreaEmitterError::Camera(
                        CameraRadianceError::InvalidEmitter { object_id: 2 }
                    ))
                ));
            }
        }
    });
}

#[test]
fn g1_thermal_chart_and_face_proposals_keep_original_area_and_exact_sample_ownership() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![
            GeometryInstance::try_new(
                2,
                ContentHash([62; 32]),
                SharedGeometry::chart(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 0.5,
                }),
                RigidTransform::identity(),
            )
            .unwrap(),
        ])
        .unwrap();
        let material = infrared_medium(1.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 2,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let emitter = thermal(1.0);
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[emitter])
            .unwrap();
        let limits = AreaEmitterLimits {
            max_triangles: 0,
            max_preparation_work: 1,
            max_lookup_work: 1,
        };
        let mut light =
            CurvedAreaEmitter::try_new_thermal_with_limits(cx, &scene, emitter, limits).unwrap();
        let other = CurvedAreaEmitter::from_camera(cx, &camera, 2, 0).unwrap();
        assert!(light.source_emitter().is_none());
        assert_eq!(
            light.source_thermal_emitter().unwrap().surface().identity(),
            emitter.surface().identity()
        );
        assert!(core::ptr::eq(
            light.source_instance(),
            &geometry.instances()[0]
        ));
        let sample = light.sample_area_with_context(cx, [0.25, 0.125]).unwrap();
        assert_eq!(sample.triangle_index(), None);
        near(sample.pdf_area_per_m2(), 1.0 / core::f64::consts::PI, 1e-14);
        near(
            sample
                .position_m()
                .delta_from(Point3::new(0.0, 0.0, 0.0))
                .dot(sample.normal()),
            0.5,
            1e-14,
        );
        let seeds = [BundleImageSeed {
            initial_direction: sample.normal(),
            policy: BundleBranchPolicy::Sequence(&[]),
        }];
        let (_, _, radiance) = light
            .launch(
                sample,
                BoundarySide::Positive,
                0,
                Vec3::new(0.0, 1.0, 0.0),
                9750.0,
                &seeds,
            )
            .unwrap();
        near(radiance, planck(9750.0), 2e-12);
        assert!(matches!(
            other.launch(
                sample,
                BoundarySide::Positive,
                0,
                Vec3::new(0.0, 1.0, 0.0),
                9750.0,
                &seeds
            ),
            Err(AreaEmitterError::SampleOwnerMismatch)
        ));
        assert!(matches!(
            CurvedAreaEmitter::try_new_thermal_faces(cx, &scene, emitter, &[0], limits),
            Err(AreaEmitterError::UnsupportedFaceSubset)
        ));
        let original_sample = (
            sample.position_m(),
            sample.normal(),
            sample.pdf_area_per_m2(),
        );
        let hot = thermal_at(1.0, 480.0);
        light.refresh_thermal_emitter(cx, hot).unwrap();
        let refreshed_sample = light.sample_area_with_context(cx, [0.25, 0.125]).unwrap();
        assert_eq!(
            original_sample,
            (
                refreshed_sample.position_m(),
                refreshed_sample.normal(),
                refreshed_sample.pdf_area_per_m2()
            )
        );
        assert_eq!(light.source_thermal_emitter(), Some(&hot));
        assert!(core::ptr::eq(
            light.source_instance(),
            &geometry.instances()[0]
        ));
        assert!(light.refresh_thermal_emitter(cx, thermal(0.5)).is_err());
        assert_eq!(light.source_thermal_emitter(), Some(&hot));

        let mesh = InstanceScene::try_new(vec![source(0.0)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &mesh,
            &media,
            vec![SurfaceBinding {
                object_id: 2,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
            .unwrap()
            .with_diffuse_thermal_emitters(cx, &[emitter])
            .unwrap();
        let mut face = CurvedAreaEmitter::from_camera_faces(
            cx,
            &camera,
            2,
            &[0],
            AreaEmitterLimits {
                max_triangles: 1,
                ..AreaEmitterLimits::default()
            },
        )
        .unwrap();
        assert_eq!(
            face.sample_area([4.0 / 9.0, 0.5]).unwrap().triangle_index(),
            Some(0)
        );
        near(face.total_area_m2(), 0.27, 1e-14);
        let cdf = match &face.support {
            SourceSupport::Mesh(triangles) => triangles.as_ptr(),
            SourceSupport::Chart => panic!("original mesh source"),
        };
        face.refresh_thermal_emitter(cx, hot).unwrap();
        assert_eq!(face.source_thermal_emitter(), Some(&hot));
        assert!(
            matches!(&face.support, SourceSupport::Mesh(triangles) if triangles.as_ptr() == cdf)
        );
        near(face.total_area_m2(), 0.27, 1e-14);
    });
}
