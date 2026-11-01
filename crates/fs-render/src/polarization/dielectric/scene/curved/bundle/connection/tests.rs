use super::super::super::SceneMedium;
use super::super::tests::{
    ball_reference, bindings, close, close_vec, glass, grin_reference, launch_basis, plane,
    reference_columns, settings, sphere, unit, with_cx,
};
use super::*;
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
