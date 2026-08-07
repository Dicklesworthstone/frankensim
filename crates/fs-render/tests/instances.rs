//! G0/G3 rigid-instance, transformed-backend, and scene-equivalence tests.
#![cfg(feature = "chart-backends")]

use asupersync::types::Budget;
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, Vec3};
use fs_render::charts::{Ray, TraceTermination, TriMesh};
use fs_render::instances::{
    GeometryInstance, InstanceBackendAudit, InstanceError, InstanceScene, RigidTransform,
    SharedGeometry,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn content(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.render-instance", label.as_bytes())
}

fn triangle() -> TriMesh {
    TriMesh::new(
        vec![[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
}

fn z_rotation(angle: f64, translation: [f64; 3]) -> RigidTransform {
    let half = angle * 0.5;
    RigidTransform::try_new([0.0, 0.0, half.sin(), half.cos()], translation).unwrap()
}

fn x_rotation(angle: f64, translation: [f64; 3]) -> RigidTransform {
    let half = angle * 0.5;
    RigidTransform::try_new([half.sin(), 0.0, 0.0, half.cos()], translation).unwrap()
}

fn close(left: f64, right: f64) {
    assert!((left - right).abs() <= 2.0e-12 * left.abs().max(right.abs()).max(1.0));
}

fn close_point(left: Point3, right: Point3) {
    close(left.x, right.x);
    close(left.y, right.y);
    close(left.z, right.z);
}

#[test]
fn rigid_transform_rejects_invalid_input_and_canonicalizes_quaternion_sign() {
    assert_eq!(
        RigidTransform::try_new([0.0; 4], [0.0; 3]),
        Err(InstanceError::InvalidTransform)
    );
    assert_eq!(
        RigidTransform::try_new([0.0, 0.0, 0.0, 2.0], [0.0; 3]),
        Err(InstanceError::InvalidTransform)
    );
    assert_eq!(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [f64::NAN, 0.0, 0.0]),
        Err(InstanceError::InvalidTransform)
    );
    let positive = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [1.0, 2.0, 3.0]).unwrap();
    let negative = RigidTransform::try_new([0.0, 0.0, 0.0, -1.0], [1.0, 2.0, 3.0]).unwrap();
    assert_eq!(positive, negative);
    assert_eq!(positive.content_identity(), negative.content_identity());
}

#[test]
fn transform_and_frame_identities_match_the_legacy_concatenated_preimages() {
    let transform = z_rotation(core::f64::consts::FRAC_PI_4, [1.25, -2.5, 3.75]);
    let mut transform_preimage = Vec::with_capacity(56);
    for value in transform
        .rotation_xyzw()
        .into_iter()
        .chain(transform.translation_m())
    {
        transform_preimage.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    let expected_transform = hash_domain(
        "org.frankensim.render.rigid-transform.v1",
        &transform_preimage,
    );
    assert_eq!(transform.content_identity(), expected_transform);

    let geometry_identity = content("legacy-frame-preimage");
    let instance = GeometryInstance::try_new(
        0x0102_0304_0506_0708,
        geometry_identity,
        SharedGeometry::mesh(triangle()),
        transform,
    )
    .unwrap();
    let mut frame_preimage = Vec::with_capacity(72);
    frame_preimage.extend_from_slice(&instance.object_id().to_le_bytes());
    frame_preimage.extend_from_slice(geometry_identity.as_bytes());
    frame_preimage.extend_from_slice(expected_transform.as_bytes());
    let expected_frame = hash_domain(
        "org.frankensim.render.geometry-instance.v1",
        &frame_preimage,
    );
    assert_eq!(instance.frame_identity(), expected_frame);
}

#[test]
fn inverse_and_composition_round_trip_points_vectors_and_identity() {
    let transform = z_rotation(core::f64::consts::FRAC_PI_2, [1.0, -2.0, 0.5]);
    let inverse = transform.inverse().unwrap();
    let point = Point3::new(0.25, -0.5, 2.0);
    let vector = Vec3::new(-0.5, 0.75, 1.0);
    close_point(
        inverse.transform_point(transform.transform_point(point)),
        point,
    );
    let round_trip_vector = inverse.transform_vector(transform.transform_vector(vector));
    close(round_trip_vector.x, vector.x);
    close(round_trip_vector.y, vector.y);
    close(round_trip_vector.z, vector.z);
    let composed = transform.compose(inverse).unwrap();
    close_point(composed.transform_point(point), point);
}

#[test]
fn identity_mesh_instance_matches_direct_hit_and_transforms_differentials() {
    with_cx(|cx| {
        let mesh = triangle();
        let bvh_fingerprint = mesh.bvh_fingerprint();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 2.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let direct = mesh.intersect(&ray).unwrap();
        let instance = GeometryInstance::try_new(
            1,
            content("triangle"),
            SharedGeometry::mesh(mesh),
            RigidTransform::identity(),
        )
        .unwrap();
        let placed = instance.intersect(cx, &ray, 4.0, 1.0e-9).unwrap().unwrap();
        assert_eq!(direct, placed.hit);
        assert_eq!(placed.object_id, 1);
        assert_eq!(
            placed.backend_audit,
            InstanceBackendAudit::Mesh { bvh_fingerprint }
        );
        assert!(placed.hit.tangent_u.is_some() && placed.hit.dp_du.is_some());
    });
}

#[test]
fn translated_rotated_mesh_preserves_t_and_world_surface_frame() {
    with_cx(|cx| {
        let transform = z_rotation(core::f64::consts::FRAC_PI_2, [2.0, 3.0, 1.0]);
        let instance = GeometryInstance::try_new(
            8,
            content("placed-triangle"),
            SharedGeometry::mesh(triangle()),
            transform,
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(2.0, 3.0, 3.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let hit = instance
            .intersect(cx, &ray, 5.0, 1.0e-9)
            .unwrap()
            .unwrap()
            .hit;
        close(hit.t, 2.0);
        close_point(hit.point, Point3::new(2.0, 3.0, 1.0));
        let normal = hit.normal.unwrap();
        close(normal.x, 0.0);
        close(normal.y, 0.0);
        close(normal.z, 1.0);
        let dp_du = hit.dp_du.unwrap();
        close(dp_du.x, 0.0);
        close(dp_du.y, 2.0);
        close(dp_du.z, 0.0);
    });
}

#[test]
fn transforming_geometry_and_camera_ray_is_rigidly_equivariant() {
    with_cx(|cx| {
        let mesh = triangle();
        let transform = x_rotation(core::f64::consts::FRAC_PI_4, [-3.0, 2.0, 0.75]);
        let local_ray = Ray {
            origin: Point3::new(0.125, -0.25, 2.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let direct = mesh.intersect(&local_ray).unwrap();
        let world_ray = Ray {
            origin: transform.transform_point(local_ray.origin),
            dir: transform.transform_vector(local_ray.dir),
        };
        let instance = GeometryInstance::try_new(
            11,
            content("equivariant-triangle"),
            SharedGeometry::mesh(mesh),
            transform,
        )
        .unwrap();
        let placed = instance
            .intersect(cx, &world_ray, 4.0, 1.0e-9)
            .unwrap()
            .unwrap()
            .hit;
        close(placed.t, direct.t);
        close_point(placed.point, transform.transform_point(direct.point));
        let expected_normal = transform.transform_vector(direct.normal.unwrap());
        let actual_normal = placed.normal.unwrap();
        close(actual_normal.x, expected_normal.x);
        close(actual_normal.y, expected_normal.y);
        close(actual_normal.z, expected_normal.z);
    });
}

#[test]
fn tangent_chart_ray_propagates_the_bounded_backend_refusal() {
    with_cx(|cx| {
        let instance = GeometryInstance::try_new(
            12,
            content("tangent-sphere"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
            RigidTransform::identity(),
        )
        .unwrap();
        assert_eq!(
            instance.intersect(
                cx,
                &Ray {
                    origin: Point3::new(-2.0, 1.0, 0.0),
                    dir: Vec3::new(1.0, 0.0, 0.0),
                },
                4.0,
                1.0e-7,
            ),
            Err(InstanceError::BackendFailure(TraceTermination::StepLimit))
        );
    });
}

#[test]
fn transformed_chart_retains_certified_audit_without_tessellation() {
    with_cx(|cx| {
        let instance = GeometryInstance::try_new(
            3,
            content("analytic-sphere"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
            z_rotation(core::f64::consts::FRAC_PI_4, [4.0, -1.0, 0.5]),
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(4.0, -1.0, 4.5),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let hit = instance.intersect(cx, &ray, 8.0, 1.0e-9).unwrap().unwrap();
        close(hit.hit.t, 3.0);
        assert!(matches!(
            hit.backend_audit,
            InstanceBackendAudit::Chart(audit)
                if audit.certified
                    && audit.certifies_hit()
                    && audit.termination == TraceTermination::Hit
        ));
    });
}

#[test]
fn shared_geometry_object_collisions_and_exact_ties_are_deterministic() {
    with_cx(|cx| {
        let shared = SharedGeometry::mesh(triangle());
        let clone = shared.clone();
        assert!(shared.ptr_eq(&clone));
        let high = GeometryInstance::try_new(
            9,
            content("same-triangle"),
            shared,
            RigidTransform::identity(),
        )
        .unwrap();
        let low = GeometryInstance::try_new(
            2,
            content("same-triangle"),
            clone,
            RigidTransform::identity(),
        )
        .unwrap();
        let duplicate = GeometryInstance::try_new(
            2,
            content("other-triangle"),
            SharedGeometry::mesh(triangle()),
            RigidTransform::identity(),
        )
        .unwrap();
        assert!(matches!(
            InstanceScene::try_new(vec![low.clone(), duplicate]),
            Err(InstanceError::DuplicateObjectId)
        ));
        let scene = InstanceScene::try_new(vec![high, low]).unwrap();
        assert_eq!(scene.instances()[0].object_id(), 2);
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        assert_eq!(
            scene
                .intersect(cx, &ray, 2.0, 1.0e-9)
                .unwrap()
                .unwrap()
                .object_id,
            2
        );
    });
}

#[test]
fn missing_normals_chart_refusals_and_cancellation_propagate() {
    with_cx(|cx| {
        let skinny = TriMesh::new(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0e-13, 0.0]],
            vec![[0, 1, 2]],
        );
        let missing = GeometryInstance::try_new(
            4,
            content("missing-normal-skinny-triangle"),
            SharedGeometry::mesh(skinny),
            RigidTransform::identity(),
        )
        .unwrap();
        let ray = Ray {
            origin: Point3::new(0.25, 2.5e-14, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        assert_eq!(
            missing.intersect(cx, &ray, 2.0, 1.0e-9),
            Err(InstanceError::MissingNormal)
        );

        let no_claim = GeometryInstance::try_new(
            6,
            content("no-claim-chart"),
            SharedGeometry::chart(NoClaimChart),
            RigidTransform::identity(),
        )
        .unwrap();
        assert_eq!(
            no_claim.intersect(cx, &ray, 2.0, 1.0e-9),
            Err(InstanceError::UncertifiedTrace)
        );
    });

    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let instance = GeometryInstance::try_new(
            5,
            content("cancelled-mesh"),
            SharedGeometry::mesh(triangle()),
            RigidTransform::identity(),
        )
        .unwrap();
        gate.request();
        assert_eq!(
            instance.intersect(
                &cx,
                &Ray {
                    origin: Point3::new(0.0, 0.0, 1.0),
                    dir: Vec3::new(0.0, 0.0, -1.0),
                },
                2.0,
                1.0e-9,
            ),
            Err(InstanceError::Cancelled)
        );
    });
}

struct NoClaimChart;

impl Chart for NoClaimChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        ChartSample {
            signed_distance: 0.0,
            gradient: None,
            lipschitz: None,
            error: NumericalCertificate::no_claim(),
        }
    }

    fn support(&self) -> Aabb {
        Aabb::WHOLE_SPACE
    }

    fn name(&self) -> &'static str {
        "instance-no-claim"
    }
}

#[test]
fn geometry_and_transform_changes_move_frame_identity_but_not_geometry_identity() {
    let shared = SharedGeometry::mesh(triangle());
    let first = GeometryInstance::try_new(
        1,
        content("immutable-mesh"),
        shared.clone(),
        RigidTransform::identity(),
    )
    .unwrap();
    let moved = GeometryInstance::try_new(
        1,
        content("immutable-mesh"),
        shared,
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.1, 0.0, 0.0]).unwrap(),
    )
    .unwrap();
    assert_eq!(first.geometry_identity(), moved.geometry_identity());
    assert_ne!(first.frame_identity(), moved.frame_identity());
}

#[test]
fn zero_object_and_geometry_identities_refuse() {
    assert!(matches!(
        GeometryInstance::try_new(
            0,
            content("mesh"),
            SharedGeometry::mesh(triangle()),
            RigidTransform::identity(),
        ),
        Err(InstanceError::InvalidObjectId)
    ));
    assert!(matches!(
        GeometryInstance::try_new(
            1,
            ContentHash([0; 32]),
            SharedGeometry::mesh(triangle()),
            RigidTransform::identity(),
        ),
        Err(InstanceError::InvalidGeometryIdentity)
    ));
}

#[test]
fn invalid_ray_and_intersection_limits_refuse_consistently() {
    with_cx(|cx| {
        let instance = GeometryInstance::try_new(
            13,
            content("input-validation-mesh"),
            SharedGeometry::mesh(triangle()),
            RigidTransform::identity(),
        )
        .unwrap();
        let valid_ray = Ray {
            origin: Point3::new(0.0, 0.0, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        for (ray, t_max, eps) in [
            (
                Ray {
                    origin: valid_ray.origin,
                    dir: Vec3::new(0.0, 0.0, 0.0),
                },
                2.0,
                1.0e-9,
            ),
            (
                Ray {
                    origin: Point3::new(f64::NAN, 0.0, 1.0),
                    dir: valid_ray.dir,
                },
                2.0,
                1.0e-9,
            ),
            (valid_ray, f64::INFINITY, 1.0e-9),
            (valid_ray, 2.0, 0.0),
        ] {
            assert_eq!(
                instance.intersect(cx, &ray, t_max, eps),
                Err(InstanceError::InvalidIntersectionInput)
            );
        }
    });
}

#[cfg(feature = "tracer")]
#[test]
fn tracer_e2e_renders_one_pose_and_observes_a_moved_pose() {
    use fs_render::spectral::lift_rgb;
    use fs_render::tracer::{
        Camera, DirectStrategy, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
        render,
    };

    fn scene(shared: SharedGeometry, translation: [f64; 3]) -> Scene {
        let transform = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], translation).unwrap();
        let instance =
            GeometryInstance::try_new(41, content("emissive-quad"), shared, transform).unwrap();
        let white = lift_rgb([1.0, 1.0, 1.0]);
        let emission = (white, 4.0);
        Scene {
            primitives: vec![Primitive {
                shape: Shape::Instance(instance),
                material: Material::Lambertian { reflectance: white },
                emission: Some(emission),
            }],
            lights: vec![RectLight {
                corner: Point3::new(translation[0] - 1.0, translation[1] - 1.0, translation[2]),
                edge_u: Vec3::new(2.0, 0.0, 0.0),
                edge_v: Vec3::new(0.0, 2.0, 0.0),
                prim: 0,
                emission,
            }],
            environment: None,
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 2.0),
                forward: Vec3::new(0.0, 0.0, -1.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.1,
            },
        }
    }

    let quad = SharedGeometry::mesh(TriMesh::new(
        vec![
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    ));
    let settings = Settings {
        width: 2,
        height: 2,
        spp: 2,
        max_depth: 2,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x1a57,
    };
    with_cx(|cx| {
        let visible = render(&scene(quad.clone(), [0.0, 0.0, 0.0]), cx, &settings).unwrap();
        let moved = render(&scene(quad, [4.0, 0.0, 0.0]), cx, &settings).unwrap();
        let visible_energy: f64 = visible.xyz.iter().flatten().sum();
        let moved_energy: f64 = moved.xyz.iter().flatten().sum();
        assert!(visible_energy.is_finite() && visible_energy > 0.0);
        assert_eq!(moved_energy.to_bits(), 0.0f64.to_bits());
        assert_ne!(visible.xyz, moved.xyz);
    });
}

#[cfg(feature = "tracer")]
#[test]
fn tracer_rejects_duplicate_ids_and_makes_exact_ties_order_independent() {
    use fs_render::spectral::lift_rgb;
    use fs_render::tracer::{
        Camera, DirectStrategy, Material, Primitive, RectLight, Sampler, Scene, Settings, Shape,
        TracerError, render,
    };

    fn tied_scene(shared: &SharedGeometry, ids: [u64; 2]) -> Scene {
        let white = lift_rgb([1.0, 1.0, 1.0]);
        let low_emission = (lift_rgb([0.1, 0.8, 0.2]), 3.0);
        let high_emission = (lift_rgb([0.8, 0.1, 0.2]), 3.0);
        let primitives = ids
            .into_iter()
            .map(|object_id| Primitive {
                shape: Shape::Instance(
                    GeometryInstance::try_new(
                        object_id,
                        content("tied-emitter"),
                        shared.clone(),
                        RigidTransform::identity(),
                    )
                    .unwrap(),
                ),
                material: Material::Lambertian { reflectance: white },
                emission: Some(if object_id == 2 {
                    low_emission
                } else {
                    high_emission
                }),
            })
            .collect::<Vec<_>>();
        let low_index = ids.iter().position(|id| *id == 2).unwrap();
        Scene {
            primitives,
            lights: vec![RectLight {
                corner: Point3::new(-1.0, -1.0, 0.0),
                edge_u: Vec3::new(2.0, 0.0, 0.0),
                edge_v: Vec3::new(0.0, 2.0, 0.0),
                prim: low_index,
                emission: low_emission,
            }],
            environment: None,
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 2.0),
                forward: Vec3::new(0.0, 0.0, -1.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.1,
            },
        }
    }

    let shared = SharedGeometry::mesh(TriMesh::new(
        vec![
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    ));
    let settings = Settings {
        width: 1,
        height: 1,
        spp: 1,
        max_depth: 1,
        sampler: Sampler::Iid,
        strategy: DirectStrategy::Mis,
        seed: 0x71e,
    };
    with_cx(|cx| {
        let high_first = render(&tied_scene(&shared, [9, 2]), cx, &settings).unwrap();
        let low_first = render(&tied_scene(&shared, [2, 9]), cx, &settings).unwrap();
        assert_eq!(high_first, low_first);
        assert!(high_first.xyz.iter().flatten().sum::<f64>() > 0.0);

        let mut duplicate = tied_scene(&shared, [2, 9]);
        duplicate.primitives[1].shape = Shape::Instance(
            GeometryInstance::try_new(
                2,
                content("duplicate-emitter"),
                shared,
                RigidTransform::identity(),
            )
            .unwrap(),
        );
        assert_eq!(
            render(&duplicate, cx, &settings),
            Err(TracerError::InvalidInstance)
        );
    });
}
