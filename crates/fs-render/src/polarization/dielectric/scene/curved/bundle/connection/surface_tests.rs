use super::super::super::SceneMedium;
use super::super::tests::{
    ball_reference, bindings, close, close_vec, glass, plane, reference_columns, settings, sphere,
    unit, with_cx,
};
use super::*;
use super::images::{BundleImageOptions, BundleImageSeed};
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::DielectricSurface;
use crate::instances::InstanceScene;
use crate::polarization::dielectric::scene::{Interaction, SceneError, SurfaceBinding};
use crate::spectral::SpectralReflectance;
use crate::surface_detail::{BrushingPattern, SurfaceDetail};

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

fn rough_glass(object_id: u64, outside: usize, inside: usize, alpha: f64) -> SurfaceBinding {
    SurfaceBinding {
        object_id,
        interaction: Interaction::Dielectric {
            outside,
            inside,
            surface: DielectricSurface::try_rough(alpha).unwrap(),
        },
    }
}

fn target(object_id: u64, point: Point3, section_normal: Vec3) -> BundleConnectionTarget {
    BundleConnectionTarget::try_surface_receiver(object_id, point, section_normal, 2e-7).unwrap()
}

#[test]
fn g1_rough_glass_receives_both_actual_media_before_any_bsdf_or_branch_choice() {
    with_cx(|cx, _| {
        let point = Point3::new(0.1, -0.05, 0.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 2.0)]).unwrap();
        let solid = glass(1.5, 0.12);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&solid),
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&air),
        ];
        // An oblique section is a numerical coordinate choice, even when its
        // orientation points away from one image's source.
        let receiver = target(9, point, Vec3::new(0.3, 0.1, 1.0));
        assert!(receiver.is_surface_receiver());
        for (side, medium) in [(1.0, 1), (-1.0, 0)] {
            let origin = Point3::new(0.2, 0.2, side * 2.0);
            let displacement = point.delta_from(origin);
            let distance = norm(displacement);
            let direction = unit(displacement);
            let initial = launch(
                Ray {
                    origin,
                    dir: unit(add(displacement, Vec3::new(0.06, -0.04, 0.0))),
                },
                medium,
            );
            let mut previous: Option<BundleConnection> = None;
            for alpha in [0.08, 0.65] {
                let scene = CurvedPolarizedScene::try_new(
                    cx,
                    &geometry,
                    &media,
                    vec![rough_glass(9, 1, 0, alpha)],
                    100_000,
                )
                .unwrap();
                let received = scene
                    .connect_bundle_with_policy(
                        cx,
                        initial,
                        receiver,
                        options(),
                        BundleBranchPolicy::Sequence(&[]),
                    )
                    .unwrap();
                let path = received.path();
                assert!(path.is_surface_receiver());
                assert!(!path.is_diffuse_receiver());
                assert_eq!(path.termination(), PathTermination::Opaque);
                assert_eq!(path.medium_index(), medium);
                assert_eq!(path.surface_events(), 1);
                assert_eq!(path.bundle().propagation().interface_events(), 0);
                assert_eq!(path.bundle().interface_log_power(), 0.0);
                assert_eq!(path.detector_reference_axis(), None);
                close_vec(path.bundle().propagation().initial_ray().dir, direction, 1e-7);
                let projection = received
                    .receiver_projection(receiver, options().projection)
                    .unwrap();
                close_vec(projection.normal(), Vec3::new(0.0, 0.0, side), 1e-12);
                close(
                    projection.area_per_solid_angle_m2(),
                    distance * distance / direction.z.abs(),
                    3e-6,
                );
                if let Some(before) = previous {
                    assert_eq!(path.bundle(), before.path().bundle());
                    assert_eq!(path.branch_signature(), before.path().branch_signature());
                }
                previous = Some(received);
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
                        consumed: 0,
                    }
                );
                // Equal IOR does not identify an original scene medium.
                assert_eq!(
                    scene
                        .connect_bundle_transmitted(
                            cx,
                            BundleConnectionLaunch {
                                initial_medium: 2,
                                ..initial
                            },
                            receiver,
                            options(),
                        )
                        .unwrap_err(),
                    BundleConnectionError::Scene(
                        SceneError::MediumMismatch {
                            object_id: 9,
                            expected: medium,
                            actual: 2,
                        }
                        .into()
                    )
                );
            }
        }
    });
}

// Independent quadratic intersection, including departure from the inside.
fn sphere_endpoint(ray: Ray) -> Point3 {
    let origin = Vec3::new(ray.origin.x, ray.origin.y, ray.origin.z);
    let b = origin.dot(ray.dir);
    let root = (b * b - origin.dot(origin) + 1.0).sqrt();
    let near = -b - root;
    ray.at(if near > 0.0 { near } else { -b + root })
}

#[test]
fn g1_native_curved_rough_receiver_uses_actual_tangent_jacobian_and_image_resolution() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![sphere()]).unwrap();
        let solid = glass(1.5, 0.0);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&solid),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![rough_glass(1, 1, 0, 0.2)],
            100_000,
        )
        .unwrap();
        for (origin, direction, medium) in [
            (
                Point3::new(0.25, 0.07, -3.0),
                Vec3::new(0.11, -0.035, 1.0),
                1,
            ),
            (
                Point3::new(0.15, -0.1, 0.1),
                Vec3::new(0.22, 0.13, 1.0),
                0,
            ),
        ] {
            let ray = Ray {
                origin,
                dir: unit(direction),
            };
            let point = sphere_endpoint(ray);
            let normal = unit(Vec3::new(point.x, point.y, point.z));
            let columns = reference_columns(ray, sphere_endpoint);
            let distance = norm(point.delta_from(origin));
            let initial = launch(
                Ray {
                    dir: unit(add(ray.dir, Vec3::new(0.01, -0.008, 0.0))),
                    ..ray
                },
                medium,
            );
            let receiver = target(1, point, normal);
            let received = scene
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap();
            assert!(received.iterations() > 0);
            assert!(received.path().is_surface_receiver());
            assert!(
                received.path().last_contact().unwrap().witness.departure_radius_m.is_none()
            );
            assert_eq!(received.path().medium_index(), medium);
            let projection = received
                .receiver_projection(receiver, options().projection)
                .unwrap();
            for (actual, expected) in projection.angular_derivatives().iter().zip(columns) {
                close_vec(actual.position_m_per_rad, expected, 5e-6);
                close(actual.position_m_per_rad.dot(projection.normal()), 0.0, 1e-12);
            }
            close(
                projection.area_per_solid_angle_m2(),
                distance * distance / normal.dot(ray.dir).abs(),
                1e-5,
            );
            // An off-surface target with the same tangential residual is still
            // refused by the full unsnapped three-dimensional endpoint check.
            assert!(matches!(
                received.receiver_projection(
                    target(1, point.offset(normal.scale(0.01)), normal),
                    options().projection,
                ),
                Err(BundleConnectionError::ContactResolution { .. })
            ));
            // Changing the numerical section orientation cannot change native
            // image parity or turn a curved receiver into a projected plane.
            for section in [normal, normal.scale(-1.0)] {
                let catalog = scene
                    .discover_bundle_images(
                        cx,
                        launch(ray, medium),
                        target(1, point, section),
                        options(),
                        &[BundleImageSeed {
                            initial_direction: ray.dir,
                            policy: BundleBranchPolicy::Sequence(&[]),
                        }],
                        BundleImageOptions::default(),
                    )
                    .unwrap();
                assert_eq!(catalog.images().len(), 1);
                let image = &catalog.images()[0];
                assert_eq!(image.projection_parity(), -1);
                // A straight point-source tangent map has sigma_min = distance.
                close(
                    image.angular_resolution_rad(),
                    (receiver.position_tolerance_m()
                        + image.connection().path().max_contact_resolution_m())
                        / distance,
                    1e-12,
                );
            }
        }
    });
}

#[test]
fn g1_smooth_ball_lens_connection_retains_two_interfaces_before_rough_receiver() {
    with_cx(|cx, _| {
        let ray = Ray {
            origin: Point3::new(0.12, 0.07, -3.0),
            dir: unit(Vec3::new(0.025, -0.015, 1.0)),
        };
        let (endpoint, outside, inside, transfer) = ball_reference(ray, 1.5);
        let geometry =
            InstanceScene::try_new(vec![sphere(), plane(9, endpoint.origin, 0.08)]).unwrap();
        let lens = glass(1.5, 0.2);
        let air = glass(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&lens),
            SceneMedium::Homogeneous(&air),
        ];
        let mut surfaces = bindings(Vec3::new(0.0, 1.0, 0.0));
        surfaces[1] = rough_glass(9, 0, 1, 0.3);
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, surfaces, 100_000).unwrap();
        let receiver = target(9, endpoint.origin, Vec3::new(0.0, 0.0, 1.0));
        let received = scene
            .connect_bundle_with_policy(
                cx,
                launch(
                    Ray {
                        dir: unit(add(ray.dir, Vec3::new(0.012, -0.009, 0.0))),
                        ..ray
                    },
                    1,
                ),
                receiver,
                options(),
                BundleBranchPolicy::Sequence(&[DielectricEvent::Transmission; 2]),
            )
            .unwrap();
        let path = received.path();
        assert!(path.is_surface_receiver());
        assert_eq!(path.medium_index(), 1);
        assert_eq!(path.surface_events(), 3);
        assert_eq!(path.bundle().propagation().interface_events(), 2);
        close(
            path.bundle().propagation().propagation().distance_m,
            outside + inside,
            3e-6,
        );
        close(path.bundle().interface_log_power().exp(), transfer, 2e-6);
        let columns = reference_columns(ray, |r| ball_reference(r, 1.5).0.origin);
        let area = cross(columns[0], columns[1]).z.abs();
        let projection = received
            .receiver_projection(receiver, options().projection)
            .unwrap();
        close(projection.area_per_solid_angle_m2(), area, 5e-5);
        close(
            path.bundle()
                .point_source_irradiance(projection.normal(), options().projection, 1.0)
                .unwrap(),
            transfer * (-0.2 * inside).exp() / area,
            5e-6,
        );
    });
}

#[test]
fn g0_surface_receiver_admits_only_original_continuous_materials_and_preserves_blockers() {
    with_cx(|cx, _| {
        let point = Point3::new(0.0, 0.0, 2.0);
        let geometry = InstanceScene::try_new(vec![plane(9, point, 1.0)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([91; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let conductors = [ConductorOptics::try_new(
            core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
            }),
            source,
        )
        .unwrap()];
        let reflectances = [SpectralReflectance::constant(0.0).unwrap()];
        let initial = launch(
            Ray {
                origin: Point3::new(0.0, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            0,
        );
        let receiver = target(9, point, Vec3::new(0.0, 0.0, 1.0));
        for interaction in [
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::try_rough(0.2).unwrap(),
            },
            Interaction::BrushedConductor {
                optics: 0,
                medium: 0,
                detail: SurfaceDetail::try_new(0.2, 3.0, BrushingPattern::Radial, 0.2, None).unwrap(),
            },
            Interaction::DepolarizingDiffuse { medium: 0, reflectance: 0 },
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &conductors,
                &[],
                &reflectances,
                vec![SurfaceBinding { object_id: 9, interaction }],
                100_000,
            )
            .unwrap();
            let received = scene
                .connect_bundle_transmitted(cx, initial, receiver, options())
                .unwrap();
            assert!(received.path().is_surface_receiver());
            assert_eq!(received.path().bundle().interface_log_power(), 0.0);
            assert_eq!(received.path().bundle().propagation().interface_events(), 0);
        }
        for interaction in [
            Interaction::Dielectric {
                outside: 0,
                inside: 0,
                surface: DielectricSurface::SMOOTH,
            },
            Interaction::SmoothConductor { optics: 0, medium: 0 },
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::SMOOTH,
            },
            Interaction::Opaque,
            Interaction::Detector { reference_axis: Vec3::new(1.0, 0.0, 0.0) },
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &conductors,
                &[],
                &reflectances,
                vec![SurfaceBinding { object_id: 9, interaction }],
                100_000,
            )
            .unwrap();
            assert_eq!(
                scene
                    .connect_bundle_transmitted(cx, initial, receiver, options())
                    .unwrap_err(),
                BundleConnectionError::Scene(CurvedBundleSceneError::UnsupportedInteraction {
                    object_id: 9,
                }),
            );
        }
        let blocked = InstanceScene::try_new(vec![
            plane(5, Point3::new(0.0, 0.0, 1.0), 1.0),
            plane(9, point, 1.0),
        ])
        .unwrap();
        for (blocker, expected) in [
            (
                Interaction::Opaque,
                BundleConnectionError::ReceiverNotReached(PathTermination::Opaque),
            ),
            (
                rough_glass(5, 0, 0, 0.2).interaction,
                BundleConnectionError::Scene(CurvedBundleSceneError::UnsupportedInteraction {
                    object_id: 5,
                }),
            ),
        ] {
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &blocked,
                &media,
                vec![
                    SurfaceBinding { object_id: 5, interaction: blocker },
                    rough_glass(9, 0, 0, 0.2),
                ],
                100_000,
            )
            .unwrap();
            assert_eq!(
                scene
                    .connect_bundle_transmitted(cx, initial, receiver, options())
                    .unwrap_err(),
                expected,
            );
        }
    });
}
