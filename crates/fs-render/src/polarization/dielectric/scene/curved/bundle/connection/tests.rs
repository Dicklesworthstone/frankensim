use super::super::super::SceneMedium;
use super::super::tests::{
    ball_reference, bindings, close, close_vec, glass, grin_reference, launch_basis, plane,
    reference_columns, settings, sphere, unit, with_cx,
};
use super::*;
use crate::dielectric::DielectricSurface;
use crate::gradient_index::AffineIndexField;
use crate::instances::InstanceScene;
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use fs_blake3::ContentHash;
use fs_geom::Aabb;

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

fn target(object_id: u64, point: Point3, tolerance: f64) -> BundleConnectionTarget {
    BundleConnectionTarget::try_new(object_id, point, Vec3::new(0.0, 0.0, -1.0), tolerance).unwrap()
}

fn detector(object_id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Detector {
            reference_axis: Vec3::new(0.0, 1.0, 0.0),
        },
    }
}

#[test]
fn g1_homogeneous_connection_solves_point_with_oblique_area_and_beer_attenuation() {
    with_cx(|cx, _| {
        let origin = Point3::new(-0.4, 0.2, -2.0);
        let point = Point3::new(0.7, -0.3, 2.5);
        let displacement = point.delta_from(origin);
        let distance = displacement.dot(displacement).sqrt();
        let direction = unit(displacement);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 2.0)]).unwrap();
        let material = glass(1.3, 0.12);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let initial = launch(
            Ray {
                origin,
                dir: unit(add(displacement, Vec3::new(0.15, -0.1, 0.0))),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let connection = scene
            .connect_bundle_transmitted(cx, initial, receiver, options())
            .unwrap();
        assert!(connection.iterations() > 0);
        assert!(connection.trace_evaluations() > 1);
        assert!(connection.residual_m() <= receiver.position_tolerance_m());
        let path = connection.path();
        assert_eq!(path.termination(), PathTermination::Detected);
        assert_eq!(path.surface_events(), 1);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, direction, 1e-7);
        close(geometric.propagation().distance_m, distance, 3e-7);
        close(
            geometric.propagation().optical_path_length_m,
            1.3 * distance,
            4e-7,
        );
        close(geometric.propagation().optical_depth, 0.12 * distance, 5e-8);
        let projection = geometric
            .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
            .unwrap();
        let area = distance * distance / direction.z;
        close(projection.area_per_solid_angle_m2(), area, 5e-6);
        let irradiance = 2.0 * (-0.12 * distance).exp() / area;
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 2.0)
                .unwrap(),
            irradiance,
            1e-7,
        );
        assert!(connection.query_work() > path.query_work());
        assert!(connection.field_evaluations() > path.field_evaluations());
    });
}

#[test]
fn g1_ball_lens_connection_recovers_independently_chosen_off_axis_launch() {
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
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        let initial_ray = Ray {
            dir: unit(add(expected_ray.dir, Vec3::new(0.012, -0.009, 0.0))),
            ..expected_ray
        };
        let initial = BundleConnectionLaunch {
            transverse_axis: launch_basis(expected_ray.dir)[0],
            ..launch(initial_ray, 1)
        };
        let receiver = target(9, endpoint.origin, 2e-7);
        let connection = scene
            .connect_bundle_transmitted(cx, initial, receiver, options())
            .unwrap();
        assert!(connection.iterations() > 0);
        assert!(connection.residual_m() <= receiver.position_tolerance_m());
        let path = connection.path();
        assert_eq!(path.surface_events(), 3);
        assert_eq!(path.bundle().propagation().interface_events(), 2);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, expected_ray.dir, 2e-6);
        close_vec(geometric.propagation().ray.dir, endpoint.dir, 2e-6);
        close(geometric.propagation().distance_m, outside + inside, 3e-6);
        close(path.bundle().interface_log_power().exp(), transfer, 2e-6);
        let projection = geometric
            .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
            .unwrap();
        close(projection.area_per_solid_angle_m2(), area, 5e-6);
        let expected_irradiance = transfer * (-0.2 * inside).exp() / area;
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 1.0)
                .unwrap(),
            expected_irradiance,
            7e-5 * expected_irradiance,
        );
    });
}

#[test]
fn g1_grin_connection_recovers_analytic_native_sphere_branch() {
    with_cx(|cx, _| {
        let expected_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: unit(Vec3::new(0.025, -0.01, 1.0)),
        };
        let (endpoint, inside, outside, optical, transfer) = grin_reference(expected_ray, true);
        let columns = reference_columns(expected_ray, |ray| grin_reference(ray, true).0.origin);
        let area = cross(columns[0], columns[1]).z.abs();
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, endpoint.origin, 0.2)]).unwrap();
        let material = glass(1.4, 0.3);
        let field = AffineIndexField::try_new(
            &material,
            expected_ray.origin,
            Vec3::new(0.2, 0.0, 0.0),
            Aabb {
                min: Point3::new(-1.5, -1.5, -2.0),
                max: Point3::new(1.5, 1.5, 3.0),
            },
        )
        .unwrap();
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Differentiable {
                field: &field,
                identity: ContentHash([56; 32]),
            },
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            bindings(Vec3::new(0.0, 1.0, 0.0)),
            100_000,
        )
        .unwrap();
        let initial = launch(
            Ray {
                dir: unit(add(expected_ray.dir, Vec3::new(0.02, 0.015, 0.0))),
                ..expected_ray
            },
            0,
        );
        let receiver = target(9, endpoint.origin, 2e-7);
        let connection = scene
            .connect_bundle_transmitted(cx, initial, receiver, options())
            .unwrap();
        assert!(connection.iterations() > 0);
        assert!(connection.residual_m() <= receiver.position_tolerance_m());
        let path = connection.path();
        assert_eq!(path.surface_events(), 2);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, expected_ray.dir, 2e-6);
        close(geometric.propagation().distance_m, inside + outside, 4e-6);
        close(geometric.propagation().optical_path_length_m, optical, 5e-6);
        let projection = geometric
            .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
            .unwrap();
        close(projection.area_per_solid_angle_m2(), area, 3e-5);
        let expected_irradiance = transfer * (-0.3 * inside).exp() / area;
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 1.0)
                .unwrap(),
            expected_irradiance,
            2e-5 * expected_irradiance,
        );
    });
}

#[test]
fn g0_connection_refuses_occlusion_wrong_detector_curvature_and_unresolved_contact() {
    with_cx(|cx, _| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let initial = launch(ray, 0);
        let receiver = target(9, Point3::new(0.0, 0.0, 2.0), 2e-7);
        let geometry = InstanceScene::try_new(vec![
            plane(5, Point3::new(0.0, 0.0, 1.0), 1.0),
            plane(9, receiver.position_m(), 1.0),
        ])
        .unwrap();
        for (interaction, error) in [
            (
                Interaction::Opaque,
                BundleConnectionError::ReceiverNotReached(PathTermination::Opaque),
            ),
            (
                detector(5).interaction,
                BundleConnectionError::WrongDetector,
            ),
        ] {
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![
                    SurfaceBinding {
                        object_id: 5,
                        interaction,
                    },
                    detector(9),
                ],
                100_000,
            )
            .unwrap();
            assert_eq!(
                scene
                    .connect_bundle_transmitted(cx, initial, receiver, options())
                    .unwrap_err(),
                error
            );
        }
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(1)], 100_000)
                .unwrap();
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, -3.0),
                ..ray
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

        let geometry = InstanceScene::try_new(vec![plane(9, receiver.position_m(), 1.0)]).unwrap();
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let initial = launch(ray, 0);
        let wrong_normal = BundleConnectionTarget::try_new(
            9,
            receiver.position_m(),
            Vec3::new(0.1, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, wrong_normal, options())
                .unwrap_err(),
            BundleConnectionError::ReceiverPlaneMismatch
        );
        let back_face = BundleConnectionTarget::try_new(
            9,
            receiver.position_m(),
            Vec3::new(0.0, 0.0, 1.0),
            2e-7,
        )
        .unwrap();
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, back_face, options())
                .unwrap_err(),
            BundleConnectionError::BackFacingReceiver
        );
        let mut coarse = options();
        coarse.scene.trace.event.surface_tolerance_m = 1e-4;
        let admissible = target(9, receiver.position_m(), 1e-3);
        let actual = scene
            .connect_bundle_transmitted(cx, initial, admissible, coarse)
            .unwrap();
        assert!(
            actual.residual_m() > 1e-12,
            "native contact must retain its nonzero unsnapped endpoint gap"
        );
        let strict = target(9, receiver.position_m(), actual.residual_m() * 0.5);
        assert!(matches!(
            scene.connect_bundle_transmitted(cx, initial, strict, coarse),
            Err(BundleConnectionError::ContactResolution { .. })
        ));
    });
}

#[test]
fn g4_connection_enforces_trace_iteration_and_path_limits_and_cancellation() {
    with_cx(|cx, gate| {
        let point = Point3::new(0.5, 0.25, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 2.0)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let mut limited = options();
        limited.max_trace_evaluations = 1;
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, receiver, limited)
                .unwrap_err(),
            BundleConnectionError::TraceBudgetExceeded
        );
        limited = options();
        limited.max_iterations = 1;
        limited.max_angular_step_rad = 1e-4;
        assert!(matches!(
            scene.connect_bundle_transmitted(cx, initial, receiver, limited),
            Err(BundleConnectionError::IterationLimit { .. })
        ));
        limited = options();
        limited.scene.trace.max_distance_m = 0.5;
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, receiver, limited)
                .unwrap_err(),
            BundleConnectionError::ReceiverNotReached(PathTermination::DistanceLimit)
        );
        gate.request();
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap_err(),
            BundleConnectionError::Cancelled
        );
    });
}

#[test]
fn g1_prescribed_reflection_connects_the_dielectric_virtual_source_image() {
    with_cx(|cx, _| {
        let origin = Point3::new(-0.5, 0.0, 1.0);
        let point = Point3::new(0.6, 0.0, 1.0);
        let virtual_displacement = Vec3::new(point.x - origin.x, 0.0, -2.0);
        let expected_direction = unit(virtual_displacement);
        let distance = virtual_displacement.dot(virtual_displacement).sqrt();
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(0.0, 0.0, 0.0), 3.0),
            plane(9, point, 0.25),
        ])
        .unwrap();
        let air = glass(1.0, 0.04);
        let substrate = glass(1.5, 0.3);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&substrate),
        ];
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
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let input = launch(
            Ray {
                origin,
                dir: unit(add(expected_direction, Vec3::new(0.025, -0.015, 0.0))),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let reflected = scene
            .connect_bundle_with_policy(
                cx,
                input,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[DielectricEvent::Reflection]),
            )
            .unwrap();
        assert!(reflected.iterations() > 0);
        assert!(reflected.trace_evaluations() > 1);
        let path = reflected.path();
        assert_eq!(path.medium_index(), 0);
        assert_eq!(path.surface_events(), 2);
        close_vec(
            path.bundle().propagation().initial_ray().dir,
            expected_direction,
            2e-7,
        );
        let ci = -expected_direction.z;
        let ct = (1.0 - (1.0 - ci * ci) / 2.25).sqrt();
        let rs = ((ci - 1.5 * ct) / (ci + 1.5 * ct)).powi(2);
        let rp = ((1.5 * ci - ct) / (1.5 * ci + ct)).powi(2);
        let transfer = 0.5 * (rs + rp);
        let area = distance * distance / ci;
        close(path.bundle().interface_log_power().exp(), transfer, 2e-8);
        close(
            path.bundle().propagation().propagation().distance_m,
            distance,
            5e-7,
        );
        close(
            path.bundle().propagation().propagation().optical_depth,
            0.04 * distance,
            5e-8,
        );
        close(
            path.bundle()
                .propagation()
                .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2(),
            area,
            5e-6,
        );
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 1.0)
                .unwrap(),
            transfer * (-0.04 * distance).exp() / area,
            2e-8,
        );

        // Exact plans never invent missing choices or ignore a requested event.
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[])
                )
                .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::BranchPlanExhausted {
                dielectric_event: 0,
            })
        );
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[
                        DielectricEvent::Reflection,
                        DielectricEvent::Transmission
                    ])
                )
                .unwrap_err(),
            BundleConnectionError::BranchPlanNotConsumed {
                planned: 2,
                consumed: 1
            }
        );
    });
}

#[test]
fn g1_prescribed_slab_ghost_retains_two_internal_reflections_and_polarized_power() {
    with_cx(|cx, _| {
        let expected = Ray {
            origin: Point3::new(0.0, 0.0, -1.0),
            dir: unit(Vec3::new(0.2, -0.04, 1.0)),
        };
        // Independent unfolded Snell construction: two outside axial metres,
        // three inside axial metres for the T-R-R-T slab image.
        let endpoint = |ray: Ray| {
            let inner_z = (1.0 - (ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y) / 2.25).sqrt();
            Point3::new(
                ray.origin.x + 2.0 * ray.dir.x / ray.dir.z + 3.0 * ray.dir.x / (1.5 * inner_z),
                ray.origin.y + 2.0 * ray.dir.y / ray.dir.z + 3.0 * ray.dir.y / (1.5 * inner_z),
                2.0,
            )
        };
        let point = endpoint(expected);
        let columns = reference_columns(expected, endpoint);
        let expected_area = cross(columns[0], columns[1]).z.abs();
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(0.0, 0.0, 0.0), 4.0),
            plane(2, Point3::new(0.0, 0.0, 1.0), 4.0),
            plane(9, point, 0.5),
        ])
        .unwrap();
        let air = glass(1.0, 0.03);
        let slab = glass(1.5, 0.2);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&slab),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                // Both native triangles face +Z. At the front plane the positive
                // side is the slab; at the rear plane it is the outside medium.
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let input = launch(
            Ray {
                dir: unit(add(expected.dir, Vec3::new(0.015, 0.012, 0.0))),
                ..expected
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let plan = [
            DielectricEvent::Transmission,
            DielectricEvent::Reflection,
            DielectricEvent::Reflection,
            DielectricEvent::Transmission,
        ];
        let connected = scene
            .connect_bundle_with_policy(
                cx,
                input,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&plan),
            )
            .unwrap();
        assert!(connected.iterations() > 0);
        let path = connected.path();
        assert_eq!(path.medium_index(), 0);
        assert_eq!(path.surface_events(), 5);
        assert_eq!(path.bundle().propagation().interface_events(), 4);
        close_vec(
            path.bundle().propagation().initial_ray().dir,
            expected.dir,
            2e-7,
        );
        close_vec(
            path.bundle().propagation().propagation().ray.dir,
            expected.dir,
            2e-7,
        );
        let ci = expected.dir.z;
        let ct = (1.0 - (1.0 - ci * ci) / 2.25).sqrt();
        let rs = ((ci - 1.5 * ct) / (ci + 1.5 * ct)).powi(2);
        let rp = ((1.5 * ci - ct) / (1.5 * ci + ct)).powi(2);
        let transfer = 0.5 * ((1.0 - rs).powi(2) * rs * rs + (1.0 - rp).powi(2) * rp * rp);
        let outside = 2.0 / ci;
        let inside = 3.0 / ct;
        let geometric = path.bundle().propagation();
        close(geometric.propagation().distance_m, outside + inside, 1e-6);
        close(
            geometric.propagation().optical_path_length_m,
            outside + 1.5 * inside,
            1e-6,
        );
        close(
            geometric.propagation().optical_depth,
            0.03 * outside + 0.2 * inside,
            2e-7,
        );
        close(path.bundle().interface_log_power().exp(), transfer, 1e-9);
        close(
            geometric
                .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2(),
            expected_area,
            2e-5,
        );
        close(
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 1.0)
                .unwrap(),
            transfer * (-0.03 * outside - 0.2 * inside).exp() / expected_area,
            2e-9,
        );

        // A shorter plan may not silently finish the remaining internal image.
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&plan[..3])
                )
                .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::BranchPlanExhausted {
                dielectric_event: 3,
            })
        );
    });
}

#[test]
fn g0_prescribed_transmission_does_not_become_reflection_at_tir() {
    use crate::gradient_index::bundle::interface::BundleInterfaceError;
    use crate::polarization::dielectric::transport::bundle::PolarizedBundleError;
    with_cx(|cx, _| {
        let origin = Point3::new(0.0, 0.0, 1.0);
        let direction = Vec3::new(0.8, 0.0, -0.6);
        let point = Point3::new(8.0 / 3.0, 0.0, 1.0);
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(0.0, 0.0, 0.0), 6.0),
            plane(9, point, 0.3),
        ])
        .unwrap();
        let inside = glass(1.5, 0.0);
        let outside = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&inside),
            SceneMedium::Homogeneous(&outside),
        ];
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
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let input = launch(
            Ray {
                origin,
                dir: direction,
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let reflected = scene
            .connect_bundle_with_policy(
                cx,
                input,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[DielectricEvent::Reflection]),
            )
            .unwrap();
        let automatic = scene
            .connect_bundle_transmitted(cx, input, receiver, options())
            .unwrap();
        assert_eq!(reflected.path(), automatic.path());
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[DielectricEvent::Transmission])
                )
                .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::Bundle(
                PolarizedBundleError::Interface(BundleInterfaceError::TotalInternalReflection)
            ))
        );
    });
}

#[test]
fn g1_dark_brewster_connection_keeps_zero_power_through_later_interfaces() {
    with_cx(|cx, _| {
        let cosine = 1.0 / 3.25_f64.sqrt();
        let point = Point3::new(31.0 / 12.0, 0.0, 1.0);
        // Plane 1 is centred at (2, 0, 0): vertices (-1, -3, 0), (5, -3, 0), (2, 3, 0).
        // Centred at (1.5, 0, 0), as first written, its right edge ran through (3, 0, 0),
        // exactly where the second Brewster reflection lands, and that contact refuses
        // (g1_brewster_second_landing_on_a_triangle_edge_refuses_as_ambiguous_contact).
        // Here both landings are interior. At y = 0 the edges are at x = 0.5 and
        // x = 3.5, so the landing at (1.5, 0, 0) is 2/sqrt(5) m from an edge and
        // the one at (3, 0, 0) is 1/sqrt(5) m, more than 1e3 contact tolerances (1e-8 m).
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(2.0, 0.0, 0.0), 3.0),
            // This finite transmitting patch lies on the reflected leg only.
            plane(2, Point3::new(2.25, 0.0, 0.5), 0.5),
            plane(9, point, 0.3),
        ])
        .unwrap();
        // Error model of the p amplitude at this Brewster angle, as in
        // transport::bundle g1_p_polarized_brewster_branch_is_dark_without_invalid_normalized_stokes:
        // |rp| <= (1.25 / 2.25) * 2 |dc| / (36/13), with |dc| <= 8 eps.
        // Each reflection's power is then at most rp_bound^2.
        let rp_bound = (1.25 / 2.25) * (16.0 * f64::EPSILON) * (1.0 + f64::EPSILON) / (36.0 / 13.0);
        let air = glass(1.0, 0.02);
        let body = glass(1.5, 0.1);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&body),
        ];
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
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let input = BundleConnectionLaunch {
            transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            normalized_source_stokes: Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap(),
            ..launch(
                Ray {
                    origin: Point3::new(0.0, 0.0, 1.0),
                    dir: Vec3::new(1.5 * cosine, 0.0, -cosine),
                },
                0,
            )
        };
        let forward = scene
            .trace_bundle_with_branches(
                cx,
                scene
                    .launch_bundle(
                        input.ray,
                        input.transverse_axis,
                        input.wavelength_nm,
                        input.normalized_source_stokes,
                        input.initial_medium,
                    )
                    .unwrap(),
                input.initial_medium,
                options().scene,
                |_, _| DielectricEvent::Reflection,
            )
            .unwrap();
        // An exactly dark first reflection stops the trace there (Extinguished).
        // A weak one goes on: the reflection at plane 2's centre and the one at
        // (3, 0, 0) on plane 1 are at Brewster's angle too. Nothing lies on the
        // upward leg from (3, 0, 0): it crosses z = 0.5 at x = 3.75, where plane 2
        // spans x in [2.0, 2.5] at y = 0, and z = 1 at x = 4.5, where plane 9 spans
        // [2.43, 2.73]. So the trace escapes after three reflections, each bounded.
        let forward_log_power = forward.bundle().interface_log_power();
        if forward_log_power == f64::NEG_INFINITY {
            assert_eq!(forward.termination(), PathTermination::Extinguished);
            assert_eq!(forward.surface_events(), 1);
        } else {
            assert_eq!(forward.termination(), PathTermination::Escaped);
            assert_eq!(forward.surface_events(), 3);
            assert!(
                forward_log_power <= 3.0 * 2.0 * rp_bound.ln(),
                "ln power after three Brewster reflections = {forward_log_power}"
            );
        }

        let receiver = target(9, point, 2e-7);
        let plan = [DielectricEvent::Reflection, DielectricEvent::Transmission];
        let connect = |normalized_source_stokes: Stokes| {
            scene
                .connect_bundle_with_policy(
                    cx,
                    BundleConnectionLaunch {
                        normalized_source_stokes,
                        ..input
                    },
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&plan),
                )
                .unwrap()
        };
        let connected = connect(input.normalized_source_stokes);
        // Sibling rule: the s-polarized beam on the same rays has the same geometry
        // and Beer factors. At this angle rs^2 = 25/169 at plane 1, and plane 2 passes
        // ts = 4 (2/6.5)(4.5/6.5) = 144/169, while the p beam passes tp = 1 there.
        // So the p power over the s power is at most rp_bound^2 / (25/169 * 144/169).
        let bright = connect(Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap());
        let sibling_ratio = rp_bound * rp_bound / (25.0 / 169.0 * (144.0 / 169.0)) * (1.0 + 1e-12);
        let path = connected.path();
        assert_eq!(path.termination(), PathTermination::Detected);
        assert_eq!(path.surface_events(), 3);
        assert_eq!(path.medium_index(), 1);
        assert_eq!(path.bundle().propagation().interface_events(), 2);
        assert_eq!(bright.path().termination(), PathTermination::Detected);
        let log_power = path.bundle().interface_log_power();
        assert!(log_power <= 2.0 * rp_bound.ln(), "ln power = {log_power}");
        let irradiance = |path: &CurvedBundlePathResult| {
            path.bundle()
                .point_source_irradiance(receiver.normal(), ProjectionOptions::default(), 1e200)
                .unwrap()
        };
        let (dark_irradiance, bright_irradiance) = (irradiance(path), irradiance(bright.path()));
        assert!(bright_irradiance.is_finite() && bright_irradiance > 0.0);
        assert!(
            dark_irradiance >= 0.0 && dark_irradiance <= bright_irradiance * sibling_ratio,
            "dark {dark_irradiance:.17e} bright {bright_irradiance:.17e}"
        );
        let state = path.bundle().interface_state();
        assert!(Stokes::try_new(state.i, state.q, state.u, state.v).is_some());

        // Geometrical validation still applies to a physically dark branch.
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    target(10, point, 2e-7),
                    options(),
                    BundleBranchPolicy::Sequence(&plan)
                )
                .unwrap_err(),
            BundleConnectionError::WrongDetector
        );
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    input,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&plan[..1])
                )
                .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::BranchPlanExhausted {
                dielectric_event: 1,
            })
        );

        #[cfg(feature = "photometry")]
        {
            use crate::instances::RigidTransform;
            use crate::photometry::source::bundle::PlanarIrradianceReceiver;
            use crate::photometry::source::{
                AngularDistribution, PhotometricSource, SourceSpectrum,
            };
            let source = PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0]).unwrap(),
                1200.0,
                1.0,
                AngularDistribution::isotropic(),
                SourceSpectrum::equal_energy(),
            )
            .unwrap();
            let receiver = PlanarIrradianceReceiver::try_new(
                point,
                receiver.normal(),
                receiver.position_tolerance_m(),
            )
            .unwrap();
            // The same sibling rule bounds the photometric irradiance.
            let photometric = |path: &CurvedBundlePathResult| {
                source
                    .spectral_irradiance_from_refractive_bundle(
                        path.bundle(),
                        receiver,
                        ProjectionOptions::default(),
                    )
                    .unwrap()
            };
            let (dark_flux, bright_flux) = (photometric(path), photometric(bright.path()));
            assert!(bright_flux.is_finite() && bright_flux > 0.0);
            assert!(
                dark_flux >= 0.0 && dark_flux <= bright_flux * sibling_ratio,
                "dark {dark_flux:.17e} bright {bright_flux:.17e}"
            );
        }
    });
}

#[test]
fn g1_brewster_second_landing_on_a_triangle_edge_refuses_as_ambiguous_contact() {
    // The geometry first written for g1_dark_brewster_connection_keeps_zero_power_through_later_interfaces,
    // with plane 1 centred at (1.5, 0, 0) and radius 3. Its right edge
    // x = 4.5 - (y + 3) / 2 runs through (3, 0, 0). A ray from (0, 0, 1) along
    // (1.5 c, 0, -c) reflects at (1.5, 0, 0) and at plane 2's centre (2.25, 0, 0.5),
    // then lands exactly on that edge.
    // The proximity query refuses a contact whose departure radius is at most 4 contact
    // tolerances: a normal is not a unique surface event at an unresolved edge
    // (instances/proximity.rs:1144). The bundle trace passes event.surface_tolerance_m
    // as that tolerance (bundle.rs:563-567).
    // Unpolarized light keeps the s power rs^2 = 25/169 at each reflection, so the
    // trace reaches the edge whatever the rounding of the Brewster cosine.
    with_cx(|cx, _| {
        let cosine = 1.0 / 3.25_f64.sqrt();
        let point = Point3::new(31.0 / 12.0, 0.0, 1.0);
        let geometry = InstanceScene::try_new(vec![
            plane(1, Point3::new(1.5, 0.0, 0.0), 3.0),
            plane(2, Point3::new(2.25, 0.0, 0.5), 0.5),
            plane(9, point, 0.3),
        ])
        .unwrap();
        let air = glass(1.0, 0.02);
        let body = glass(1.5, 0.1);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&body),
        ];
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
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 1.0),
            dir: Vec3::new(1.5 * cosine, 0.0, -cosine),
        };
        let result = scene.trace_bundle_with_branches(
            cx,
            scene
                .launch_bundle(
                    ray,
                    Vec3::new(0.0, 1.0, 0.0),
                    550.0,
                    Stokes::unpolarized(1.0),
                    0,
                )
                .unwrap(),
            0,
            options().scene,
            |_, _| DielectricEvent::Reflection,
        );
        match result {
            Err(CurvedBundleSceneError::Scene(
                crate::polarization::dielectric::scene::curved::CurvedSceneError::Geometry(
                    crate::instances::proximity::ProximityError::AmbiguousContact,
                ),
            )) => {}
            Err(other) => panic!("expected AmbiguousContact at the edge, got {other:?}"),
            Ok(path) => panic!(
                "expected AmbiguousContact at the edge, got {:?} after {} events",
                path.termination(),
                path.surface_events()
            ),
        }
    });
}

mod surface_sources;
