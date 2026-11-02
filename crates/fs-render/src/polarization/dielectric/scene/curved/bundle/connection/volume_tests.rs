use super::super::super::SceneMedium;
use super::super::tests::{
    ball_reference, close, close_vec, glass, grin_reference, plane, reference_columns, settings,
    sphere, unit, with_cx,
};
use super::*;
use crate::dielectric::DielectricSurface;
use crate::gradient_index::AffineIndexField;
use crate::instances::InstanceScene;
use crate::polarization::dielectric::scene::{Interaction, SceneError, SurfaceBinding};
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

fn target(medium: usize, point: Point3) -> BundleConnectionTarget {
    BundleConnectionTarget::try_volume_receiver(medium, point, Vec3::new(0.0, 0.0, -1.0), 2e-7)
        .unwrap()
}

fn sphere_binding() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::Dielectric {
            outside: 1,
            inside: 0,
            surface: DielectricSurface::SMOOTH,
        },
    }
}

#[test]
fn g1_volume_receiver_has_transverse_inverse_square_power_without_scene_surface() {
    with_cx(|cx, _| {
        let origin = Point3::new(-0.4, 0.2, -2.0);
        let point = Point3::new(0.7, -0.3, 2.5);
        let displacement = point.delta_from(origin);
        let distance = displacement.dot(displacement).sqrt();
        let direction = unit(displacement);
        let geometry = InstanceScene::try_new(vec![]).unwrap();
        let material = glass(1.3, 0.12);
        let media = [SceneMedium::Homogeneous(&material)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![], 100_000).unwrap();
        let initial = launch(
            Ray {
                origin,
                dir: unit(add(displacement, Vec3::new(0.15, -0.1, 0.0))),
            },
            0,
        );
        let receiver = target(0, point);
        let connection = scene
            .connect_bundle_transmitted(cx, initial, receiver, options())
            .unwrap();
        assert!(connection.iterations() > 0);
        assert!(connection.residual_m() <= receiver.position_tolerance_m());
        assert!(receiver.is_volume_receiver());
        assert_eq!(receiver.volume_medium_index(), Some(0));
        let path = connection.path();
        assert!(path.is_volume_receiver());
        assert!(!path.is_diffuse_receiver());
        assert_eq!(path.surface_events(), 0);
        assert_eq!(path.last_contact(), None);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, direction, 1e-7);
        close(geometric.propagation().distance_m, distance, 3e-7);
        close(geometric.propagation().optical_depth, 0.12 * distance, 5e-8);
        let transverse = geometric.propagation().ray.dir.scale(-1.0);
        let projection = geometric
            .project_at_endpoint(transverse, ProjectionOptions::default())
            .unwrap();
        close(
            projection.area_per_solid_angle_m2(),
            distance * distance,
            5e-6,
        );
        let section = geometric
            .project_at_endpoint(receiver.normal(), ProjectionOptions::default())
            .unwrap();
        close(
            section.area_per_solid_angle_m2(),
            distance * distance / direction.z,
            5e-6,
        );
        close(
            path.bundle()
                .point_source_irradiance(transverse, ProjectionOptions::default(), 2.0)
                .unwrap(),
            2.0 * (-0.12 * distance).exp() / (distance * distance),
            1e-7,
        );
    });
}

#[test]
fn g1_volume_receiver_behind_ball_lens_retains_native_refraction_and_transverse_jacobian() {
    with_cx(|cx, _| {
        let ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let (endpoint, outside, inside, transfer) = ball_reference(ray, 1.5);
        let columns = reference_columns(ray, |ray| ball_reference(ray, 1.5).0.origin);
        let transverse_area = cross(columns[0], columns[1]).dot(endpoint.dir).abs();
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![sphere_binding()], 100_000)
                .unwrap();
        let initial = launch(
            Ray {
                dir: unit(add(ray.dir, Vec3::new(0.012, -0.009, 0.0))),
                ..ray
            },
            1,
        );
        let connection = scene
            .connect_bundle_transmitted(cx, initial, target(1, endpoint.origin), options())
            .unwrap();
        let path = connection.path();
        assert!(path.is_volume_receiver());
        assert_eq!(path.surface_events(), 2);
        assert_eq!(path.last_contact().unwrap().object_id, 1);
        assert_eq!(path.medium_index(), 1);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, ray.dir, 1e-6);
        close_vec(geometric.propagation().ray.dir, endpoint.dir, 2e-6);
        close(geometric.propagation().distance_m, outside + inside, 3e-6);
        let transverse = geometric.propagation().ray.dir.scale(-1.0);
        close(
            geometric
                .project_at_endpoint(transverse, ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2(),
            transverse_area,
            4e-5,
        );
        let expected = transfer * (-0.2 * inside).exp() / transverse_area;
        close(
            path.bundle()
                .point_source_irradiance(transverse, ProjectionOptions::default(), 1.0)
                .unwrap(),
            expected,
            7e-5 * expected,
        );
    });
}

#[test]
fn g1_volume_receiver_after_grin_field_preserves_original_curvature_and_medium() {
    with_cx(|cx, _| {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: unit(Vec3::new(0.025, -0.01, 1.0)),
        };
        let (endpoint, inside, outside, optical, transfer) = grin_reference(ray, true);
        let columns = reference_columns(ray, |r| grin_reference(r, true).0.origin);
        let transverse_area = cross(columns[0], columns[1]).dot(endpoint.dir).abs();
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let material = glass(1.4, 0.3);
        let field = AffineIndexField::try_new(
            &material,
            ray.origin,
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
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![sphere_binding()], 100_000)
                .unwrap();
        let initial = launch(
            Ray {
                dir: unit(add(ray.dir, Vec3::new(0.02, 0.015, 0.0))),
                ..ray
            },
            0,
        );
        let connection = scene
            .connect_bundle_transmitted(cx, initial, target(1, endpoint.origin), options())
            .unwrap();
        let path = connection.path();
        assert!(path.is_volume_receiver());
        assert_eq!(path.surface_events(), 1);
        assert_eq!(path.medium_index(), 1);
        let geometric = path.bundle().propagation();
        close_vec(geometric.initial_ray().dir, ray.dir, 2e-6);
        close(geometric.propagation().distance_m, inside + outside, 4e-6);
        close(geometric.propagation().optical_path_length_m, optical, 5e-6);
        let transverse = geometric.propagation().ray.dir.scale(-1.0);
        close(
            geometric
                .project_at_endpoint(transverse, ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2(),
            transverse_area,
            3e-5,
        );
        let expected = transfer * (-0.3 * inside).exp() / transverse_area;
        close(
            path.bundle()
                .point_source_irradiance(transverse, ProjectionOptions::default(), 1.0)
                .unwrap(),
            expected,
            2e-5 * expected,
        );
    });
}

#[test]
fn g0_volume_receiver_does_not_hide_geometry_medium_plan_or_work_refusals() {
    with_cx(|cx, gate| {
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&air),
        ];
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        let receiver = target(0, Point3::new(0.0, 0.0, 2.0));
        let geometry =
            InstanceScene::try_new(vec![plane(5, Point3::new(0.0, 0.0, 1.0), 2.0)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 5,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap_err(),
            BundleConnectionError::ReceiverNotReached(PathTermination::Opaque)
        );

        let geometry = InstanceScene::try_new(vec![]).unwrap();
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![], 100_000).unwrap();
        assert!(matches!(
            scene.connect_bundle_transmitted(
                cx,
                initial,
                target(1, receiver.position_m()),
                options(),
            ),
            Err(BundleConnectionError::ReceiverNotReached(_))
        ));
        assert!(matches!(
            scene.connect_bundle_transmitted(
                cx,
                initial,
                target(2, receiver.position_m()),
                options(),
            ),
            Err(BundleConnectionError::Scene(CurvedBundleSceneError::Scene(
                CurvedSceneError::Scene(SceneError::InvalidInitialMedium)
            )))
        ));
        assert_eq!(
            scene
                .connect_bundle_with_policy(
                    cx,
                    initial,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[DielectricEvent::Transmission]),
                )
                .unwrap_err(),
            BundleConnectionError::BranchPlanNotConsumed {
                planned: 1,
                consumed: 0
            }
        );
        let mut capped = options();
        capped.max_trace_evaluations = 1;
        let perturbed = BundleConnectionLaunch {
            ray: Ray {
                dir: unit(Vec3::new(0.1, 0.0, 1.0)),
                ..initial.ray
            },
            ..initial
        };
        assert_eq!(
            scene
                .connect_bundle_transmitted(cx, perturbed, receiver, capped)
                .unwrap_err(),
            BundleConnectionError::TraceBudgetExceeded
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
