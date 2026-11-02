use super::super::super::SceneMedium;
use super::super::tests::{
    ball_reference, bindings, close, close_vec, glass, plane, reference_columns, settings, sphere,
    unit, with_cx,
};
use super::*;
use crate::instances::InstanceScene;
use crate::polarization::dielectric::scene::{Interaction, SceneError, SurfaceBinding};
use crate::spectral::SpectralReflectance;

fn options() -> BundleConnectionOptions {
    BundleConnectionOptions {
        scene: settings(),
        ..BundleConnectionOptions::default()
    }
}

fn launch(ray: Ray, initial_medium: usize) -> BundleConnectionLaunch {
    BundleConnectionLaunch {
        ray,
        transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        wavelength_nm: 550.0,
        normalized_source_stokes: Stokes::unpolarized(1.0),
        initial_medium,
    }
}

fn diffuse(object_id: u64, medium: usize) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::DepolarizingDiffuse {
            medium,
            reflectance: 0,
        },
    }
}

fn target(object_id: u64, point: Point3, tolerance: f64) -> BundleConnectionTarget {
    BundleConnectionTarget::try_diffuse_receiver(
        object_id,
        point,
        Vec3::new(0.0, 0.0, -1.0),
        tolerance,
    )
    .unwrap()
}

#[test]
fn g1_diffuse_endpoint_retains_incident_bundle_and_detector_irradiance_bits() {
    with_cx(|cx, _| {
        let origin = Point3::new(-0.3, 0.1, -1.0);
        let point = Point3::new(0.4, -0.2, 2.0);
        let displacement = point.delta_from(origin);
        let distance = displacement.dot(displacement).sqrt();
        let direction = unit(displacement);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 2.0)]).unwrap();
        let material = glass(1.3, 0.12);
        let media = [SceneMedium::Homogeneous(&material)];
        let detector = SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Detector {
                reference_axis: Vec3::new(0.0, 1.0, 0.0),
            },
        };
        let detector_scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector], 100_000).unwrap();
        let initial = launch(
            Ray {
                origin,
                dir: unit(add(displacement, Vec3::new(0.15, -0.1, 0.0))),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        assert!(receiver.is_diffuse_receiver());
        let detector_target = BundleConnectionTarget::try_new(
            9,
            point,
            receiver.normal(),
            receiver.position_tolerance_m(),
        )
        .unwrap();
        assert!(!detector_target.is_diffuse_receiver());
        let detected = detector_scene
            .connect_bundle_transmitted(cx, initial, detector_target, options())
            .unwrap();
        assert!(!detected.path().is_diffuse_receiver());
        assert_eq!(detected.path().termination(), PathTermination::Detected);
        for rho in [0.0, 0.6, 1.0] {
            let reflectances = [SpectralReflectance::constant(rho).unwrap()];
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![diffuse(9, 0)],
                100_000,
            )
            .unwrap();
            let received = scene
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap();
            let path = received.path();
            assert!(path.is_diffuse_receiver());
            assert_eq!(path.termination(), PathTermination::Opaque);
            assert_eq!(path.detector_reference_axis(), None);
            assert_eq!(path.surface_events(), 1);
            assert_eq!(path.bundle(), detected.path().bundle());
            assert_eq!(received.residual_m(), detected.residual_m());
            assert_eq!(received.iterations(), detected.iterations());
            assert_eq!(received.trace_evaluations(), detected.trace_evaluations());
            assert_ne!(path.branch_signature(), detected.path().branch_signature());
            let geometric = path.bundle().propagation();
            close_vec(geometric.initial_ray().dir, direction, 1e-7);
            close(geometric.propagation().distance_m, distance, 3e-7);
            let area = distance * distance / direction.z;
            close(
                geometric
                    .project_at_endpoint(receiver.normal(), options().projection)
                    .unwrap()
                    .area_per_solid_angle_m2(),
                area,
                5e-6,
            );
            // This is incident irradiance, independent even of a black
            // receiver's reflectance; receiving does not execute a BSDF.
            close(
                path.bundle()
                    .point_source_irradiance(receiver.normal(), options().projection, 2.0)
                    .unwrap(),
                2.0 * (-0.12 * distance).exp() / area,
                1e-7,
            );
        }
        // Both original detector connection APIs keep their exact result.
        let prescribed = detector_scene
            .connect_bundle_with_policy(
                cx,
                initial,
                detector_target,
                options(),
                BundleBranchPolicy::Sequence(&[]),
            )
            .unwrap();
        assert_eq!(prescribed, detected);
    });
}

#[test]
fn g1_ball_lens_diffuse_connection_preserves_fresnel_beer_jacobian_and_exact_plan() {
    with_cx(|cx, _| {
        let expected_ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let (endpoint, outside, inside, transfer) = ball_reference(expected_ray, 1.5);
        let columns = reference_columns(expected_ray, |ray| ball_reference(ray, 1.5).0.origin);
        let area = cross(columns[0], columns[1]).z.abs();
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, endpoint.origin, 0.08)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let reflectances = [SpectralReflectance::constant(0.4).unwrap()];
        let mut surfaces = bindings(Vec3::new(0.0, 1.0, 0.0));
        surfaces[1] = diffuse(9, 1);
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            surfaces,
            100_000,
        )
        .unwrap();
        let initial = launch(
            Ray {
                dir: unit(add(expected_ray.dir, Vec3::new(0.012, -0.009, 0.0))),
                ..expected_ray
            },
            1,
        );
        let receiver = target(9, endpoint.origin, 2e-7);
        let plan = [DielectricEvent::Transmission; 2];
        let received = scene
            .connect_bundle_with_policy(
                cx,
                initial,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&plan),
            )
            .unwrap();
        let path = received.path();
        assert!(path.is_diffuse_receiver());
        assert_eq!(path.surface_events(), 3);
        assert_eq!(path.bundle().propagation().interface_events(), 2);
        close_vec(
            path.bundle().propagation().initial_ray().dir,
            expected_ray.dir,
            2e-6,
        );
        close(
            path.bundle().propagation().propagation().distance_m,
            outside + inside,
            3e-6,
        );
        close(path.bundle().interface_log_power().exp(), transfer, 2e-6);
        let expected = transfer * (-0.2 * inside).exp() / area;
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), options().projection, 1.0)
                .unwrap(),
            expected,
            7e-5 * expected,
        );
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    initial,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[DielectricEvent::Transmission; 3])
                )
                .unwrap_err(),
            BundleConnectionError::BranchPlanNotConsumed {
                planned: 3,
                consumed: 2
            }
        );
        assert!(matches!(
            scene.connect_bundle_with_policy(
                cx,
                initial,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[DielectricEvent::Transmission])
            ),
            Err(BundleConnectionError::Scene(
                CurvedBundleSceneError::BranchPlanExhausted {
                    dielectric_event: 1,
                }
            ))
        ));
    });
}

#[test]
fn g0_only_the_selected_original_diffuse_surface_can_receive_the_bundle() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![
            plane(5, Point3::new(0.0, 0.0, 1.0), 1.0),
            plane(9, point, 1.0),
        ])
        .unwrap();
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&air),
        ];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        for (interaction, expected) in [
            (
                Interaction::Opaque,
                BundleConnectionError::ReceiverNotReached(PathTermination::Opaque),
            ),
            (
                diffuse(5, 0).interaction,
                BundleConnectionError::Scene(CurvedBundleSceneError::UnsupportedInteraction {
                    object_id: 5,
                }),
            ),
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![
                    SurfaceBinding {
                        object_id: 5,
                        interaction,
                    },
                    diffuse(9, 0),
                ],
                100_000,
            )
            .unwrap();
            assert_eq!(
                scene
                    .connect_bundle_transmitted(cx, initial, receiver, options())
                    .unwrap_err(),
                expected
            );
        }
        let geometry = InstanceScene::try_new(vec![plane(9, point, 1.0)]).unwrap();
        let opaque = CurvedPolarizedScene::try_new(
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
        assert_eq!(
            opaque
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap_err(),
            BundleConnectionError::ReceiverNotReached(PathTermination::Opaque)
        );
        let mismatch = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(9, 1)],
            100_000,
        )
        .unwrap();
        assert_eq!(
            mismatch
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap_err(),
            BundleConnectionError::Scene(
                SceneError::MediumMismatch {
                    object_id: 9,
                    expected: 1,
                    actual: 0,
                }
                .into()
            )
        );
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(9, 0)],
            100_000,
        )
        .unwrap();
        let detector = BundleConnectionTarget::try_new(9, point, receiver.normal(), 2e-7).unwrap();
        assert!(matches!(
            scene.connect_bundle_transmitted(cx, initial, detector, options()),
            Err(BundleConnectionError::Scene(
                CurvedBundleSceneError::UnsupportedInteraction { object_id: 9 }
            ))
        ));
        let missed = BundleConnectionLaunch {
            ray: Ray {
                dir: unit(Vec3::new(2.0, 0.0, 1.0)),
                ..initial.ray
            },
            ..initial
        };
        assert!(matches!(
            scene.connect_bundle_transmitted(cx, missed, receiver, options()),
            Err(BundleConnectionError::ReceiverNotReached(_))
        ));
    });
}

#[test]
fn g0_diffuse_receivers_keep_planarity_orientation_resolution_work_and_cancellation_checks() {
    with_cx(|cx, gate| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0)],
            100_000,
        )
        .unwrap();
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, -3.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        assert_eq!(
            scene
                .connect_bundle_transmitted(
                    cx,
                    initial,
                    target(1, Point3::new(0.0, 0.0, -1.0), 2e-7),
                    options()
                )
                .unwrap_err(),
            BundleConnectionError::NonPlanarReceiver
        );
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 1.0)]).unwrap();
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(9, 0)],
            100_000,
        )
        .unwrap();
        let initial = BundleConnectionLaunch {
            ray: Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                ..initial.ray
            },
            ..initial
        };
        for (normal, expected) in [
            (
                Vec3::new(0.1, 0.0, -1.0),
                BundleConnectionError::ReceiverPlaneMismatch,
            ),
            (
                Vec3::new(0.0, 0.0, 1.0),
                BundleConnectionError::BackFacingReceiver,
            ),
        ] {
            let receiver =
                BundleConnectionTarget::try_diffuse_receiver(9, point, normal, 2e-7).unwrap();
            assert_eq!(
                scene
                    .connect_bundle_transmitted(cx, initial, receiver, options())
                    .unwrap_err(),
                expected
            );
        }
        let mut coarse = options();
        coarse.scene.trace.event.surface_tolerance_m = 1e-4;
        let received = scene
            .connect_bundle_transmitted(cx, initial, target(9, point, 1e-3), coarse)
            .unwrap();
        assert!(
            received.residual_m() > 1e-12,
            "do not snap the native endpoint onto the plane"
        );
        assert!(matches!(
            scene.connect_bundle_transmitted(
                cx,
                initial,
                target(9, point, received.residual_m() * 0.5),
                coarse
            ),
            Err(BundleConnectionError::ContactResolution { .. })
        ));
        let mut limited = options();
        limited.max_trace_evaluations = 1;
        let perturbed = BundleConnectionLaunch {
            ray: Ray {
                dir: unit(Vec3::new(0.1, 0.0, 1.0)),
                ..initial.ray
            },
            ..initial
        };
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, perturbed, target(9, point, 2e-7), limited)
                .unwrap_err(),
            BundleConnectionError::TraceBudgetExceeded
        );
        gate.request();
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, target(9, point, 2e-7), options())
                .unwrap_err(),
            BundleConnectionError::Cancelled
        );
    });
}
