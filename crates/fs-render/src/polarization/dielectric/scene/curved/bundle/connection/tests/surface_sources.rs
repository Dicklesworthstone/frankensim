//! Native emission departures retain the source point and all later geometry.

use super::*;
use crate::charts::TriMesh;
use crate::gradient_index::event::BoundarySide;
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::SceneError;
use crate::spectral::SpectralReflectance;
use fs_blake3::hash_domain;

#[test]
fn g1_surface_sources_keep_native_origin_area_and_absorption_without_a_source_bsdf() {
    with_cx(|cx, _| {
        let material = glass(1.3, 0.12);
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.03).unwrap()];
        let origin = Point3::new(0.0, 0.0, 0.0);
        for interaction in [
            Interaction::Opaque,
            Interaction::DepolarizingDiffuse {
                medium: 0,
                reflectance: 0,
            },
        ] {
            for sign in [-1.0, 1.0] {
                let point = Point3::new(0.35, -0.2, 2.0 * sign);
                let displacement = point.delta_from(origin);
                let distance = displacement.dot(displacement).sqrt();
                let direction = unit(displacement);
                let geometry = InstanceScene::try_new(vec![
                    plane(2, origin, 0.5),
                    plane(9, point, 1.0),
                ])
                .unwrap();
                let scene = CurvedPolarizedScene::try_new_with_materials(
                    cx,
                    &geometry,
                    &media,
                    &[],
                    &[],
                    &reflectances,
                    vec![
                        SurfaceBinding {
                            object_id: 2,
                            interaction,
                        },
                        detector(9),
                    ],
                    100_000,
                )
                .unwrap();
                let side = if sign > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::Negative
                };
                let source = BundleSurfaceSource::try_new(2, side).unwrap();
                let receiver = BundleConnectionTarget::try_new(
                    9,
                    point,
                    Vec3::new(0.0, 0.0, -sign),
                    2e-7,
                )
                .unwrap();
                let initial = launch(
                    Ray {
                        origin,
                        dir: unit(Vec3::new(
                            displacement.x + 0.04,
                            displacement.y - 0.03,
                            displacement.z,
                        )),
                    },
                    0,
                );
                let connection = scene
                    .connect_bundle_from_surface_with_policy(
                        cx,
                        initial,
                        source,
                        receiver,
                        options(),
                        BundleBranchPolicy::Sequence(&[]),
                    )
                    .unwrap();
                let path = connection.path();
                let geometric = path.bundle().propagation();
                assert_eq!(geometric.initial_ray().origin, origin);
                assert_eq!(path.termination(), PathTermination::Detected);
                assert_eq!(path.surface_events(), 1, "emission is not a surface response");
                assert_eq!(geometric.interface_events(), 0);
                assert_eq!(path.bundle().interface_log_power(), 0.0);
                assert_eq!(path.medium_index(), 0);
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
                let area = distance * distance / direction.z.abs();
                close(projection.area_per_solid_angle_m2(), area, 5e-6);
                close(
                    path.bundle()
                        .point_source_irradiance(
                            receiver.normal(),
                            ProjectionOptions::default(),
                            2.0,
                        )
                        .unwrap(),
                    2.0 * (-0.12 * distance).exp() / area,
                    1e-7,
                );
                assert!(connection.iterations() > 0);
                assert!(connection.residual_m() <= receiver.position_tolerance_m());
            }
        }
    });
}

#[test]
fn g1_surface_source_ball_lens_matches_independent_snell_area_and_power() {
    with_cx(|cx, _| {
        let expected_ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let (endpoint, outside, inside, transfer) = ball_reference(expected_ray, 1.5);
        let columns = reference_columns(expected_ray, |ray| ball_reference(ray, 1.5).0.origin);
        let area = cross(columns[0], columns[1]).z.abs();
        let geometry = InstanceScene::try_new(vec![
            sphere(),
            plane(2, expected_ray.origin, 0.5),
            plane(9, endpoint.origin, 0.08),
        ])
        .unwrap();
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
            vec![
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
                    interaction: Interaction::Opaque,
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let initial = BundleConnectionLaunch {
            transverse_axis: launch_basis(expected_ray.dir)[0],
            ..launch(
                Ray {
                    dir: unit(Vec3::new(
                        expected_ray.dir.x + 0.012,
                        expected_ray.dir.y - 0.009,
                        expected_ray.dir.z,
                    )),
                    ..expected_ray
                },
                1,
            )
        };
        let receiver = target(9, endpoint.origin, 2e-7);
        let connection = scene
            .connect_bundle_from_surface_with_policy(
                cx,
                initial,
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[
                    DielectricEvent::Transmission,
                    DielectricEvent::Transmission,
                ]),
            )
            .unwrap();
        let path = connection.path();
        let geometric = path.bundle().propagation();
        assert_eq!(geometric.initial_ray().origin, expected_ray.origin);
        assert_eq!(path.surface_events(), 3);
        assert_eq!(geometric.interface_events(), 2);
        assert_eq!(path.medium_index(), 1);
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
fn g0_surface_source_rejects_wrong_object_point_medium_and_departure_side() {
    assert_eq!(
        BundleSurfaceSource::try_new(2, BoundarySide::OnSurface).unwrap_err(),
        CurvedBundleSceneError::InvalidSurfaceSource { object_id: 2 }
    );
    with_cx(|cx, _| {
        let material = glass(1.2, 0.0);
        let other = glass(1.4, 0.0);
        let media = [
            SceneMedium::Homogeneous(&material),
            SceneMedium::Homogeneous(&other),
        ];
        let reflectances = [SpectralReflectance::constant(0.4).unwrap()];
        let origin = Point3::new(0.0, 0.0, 0.0);
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![
            plane(2, origin, 0.5),
            plane(9, point, 1.0),
        ])
        .unwrap();
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                detector(9),
            ],
            100_000,
        )
        .unwrap();
        let initial = launch(
            Ray {
                origin,
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        let receiver = target(9, point, 2e-7);
        let connect = |initial, source| {
            scene.connect_bundle_from_surface_with_policy(
                cx,
                initial,
                source,
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[]),
            )
        };
        for (object, side, ray) in [
            (42, BoundarySide::Positive, initial.ray),
            (9, BoundarySide::Positive, initial.ray),
            (
                2,
                BoundarySide::Positive,
                Ray {
                    origin: Point3::new(0.0, 0.0, 1e-3),
                    ..initial.ray
                },
            ),
        ] {
            assert_eq!(
                connect(
                    BundleConnectionLaunch { ray, ..initial },
                    BundleSurfaceSource::try_new(object, side).unwrap(),
                )
                .unwrap_err(),
                BundleConnectionError::Scene(CurvedBundleSceneError::InvalidSurfaceSource {
                    object_id: object,
                })
            );
        }
        let excluded = connect(
            initial,
            BundleSurfaceSource::try_new(2, BoundarySide::Negative).unwrap(),
        ).unwrap_err();
        let BundleConnectionError::Scene(CurvedBundleSceneError::SourceSideExcluded {
            object_id: 2, field_evaluations, query_work,
        }) = excluded else { panic!("wrong source side must retain its explicit exclusion"); };
        assert!(field_evaluations > 0 && query_work > 0);
        assert_eq!(
            connect(
                BundleConnectionLaunch {
                    initial_medium: 1,
                    ..initial
                },
                BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
            )
            .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::from(
                SceneError::MediumMismatch {
                    object_id: 2,
                    expected: 0,
                    actual: 1,
                }
            ))
        );

        // A real detector contact is not an emitting departure, even if its
        // position and outward side have otherwise been described correctly.
        let detector_scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![detector(2), detector(9)],
            100_000,
        )
        .unwrap();
        assert_eq!(
            detector_scene
                .connect_bundle_from_surface_with_policy(
                    cx,
                    initial,
                    BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                )
                .unwrap_err(),
            BundleConnectionError::Scene(CurvedBundleSceneError::InvalidSurfaceSource {
                object_id: 2,
            })
        );
    });
}

#[test]
fn g0_surface_departure_retains_other_blockers_and_later_source_object_patches() {
    with_cx(|cx, _| {
        let material = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&material)];
        let origin = Point3::new(0.0, 0.0, 0.0);
        let point = Point3::new(0.0, 0.0, 2.0);
        for same_object in [false, true] {
            let (instances, bindings) = if same_object {
                let source = GeometryInstance::try_new(
                    2,
                    hash_domain("surface-source-departure-test", b"two-original-source-patches"),
                    SharedGeometry::mesh(TriMesh::new(
                        vec![
                            [-0.5, -0.5, 0.0],
                            [0.5, -0.5, 0.0],
                            [0.0, 0.5, 0.0],
                            [-0.5, -0.5, 1.0],
                            [0.5, -0.5, 1.0],
                            [0.0, 0.5, 1.0],
                        ],
                        vec![[0, 1, 2], [3, 4, 5]],
                    )),
                    RigidTransform::identity(),
                )
                .unwrap();
                (
                    vec![source, plane(9, point, 1.0)],
                    vec![
                        SurfaceBinding {
                            object_id: 2,
                            interaction: Interaction::Opaque,
                        },
                        detector(9),
                    ],
                )
            } else {
                (
                    vec![
                        plane(2, origin, 0.5),
                        plane(5, Point3::new(0.0, 0.0, 1.0), 0.5),
                        plane(9, point, 1.0),
                    ],
                    vec![
                        SurfaceBinding {
                            object_id: 2,
                            interaction: Interaction::Opaque,
                        },
                        SurfaceBinding {
                            object_id: 5,
                            interaction: Interaction::Opaque,
                        },
                        detector(9),
                    ],
                )
            };
            let geometry = InstanceScene::try_new(instances).unwrap();
            let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings, 100_000)
                .unwrap();
            assert_eq!(
                scene
                    .connect_bundle_from_surface_with_policy(
                        cx,
                        launch(
                            Ray {
                                origin,
                                dir: Vec3::new(0.0, 0.0, 1.0),
                            },
                            0,
                        ),
                        BundleSurfaceSource::try_new(2, BoundarySide::Positive).unwrap(),
                        target(9, point, 2e-7),
                        options(),
                        BundleBranchPolicy::Sequence(&[]),
                    )
                    .unwrap_err(),
                BundleConnectionError::ReceiverNotReached(PathTermination::Opaque),
                "a departure only admits the initial source contact, same_object={same_object}",
            );
        }
    });
}
