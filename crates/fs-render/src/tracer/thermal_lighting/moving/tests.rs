//! G1/G3/G4: source density, shared rigid time, replay, and explicit refusals.

use std::sync::Arc;

use super::*;
use crate::animated_instances::{RigidTransformTrajectory, TransformKeyframe};
use crate::charts::TriMesh;
use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
use crate::thermal_surface::GrayDiffuseThermalSurface;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        operation(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 417,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn source(id: u64, size: f64, end: RigidTransform, temperature: f64) -> Primitive {
    let trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(0.0, RigidTransform::identity(), end.translation_m()).unwrap(),
        TransformKeyframe::try_new(1.0, end, end.translation_m()).unwrap(),
    ])
    .unwrap();
    Primitive {
        shape: Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                id,
                ContentHash([41; 32]),
                SharedGeometry::Mesh(Arc::new(TriMesh::new(
                    vec![[0.0; 3], [size, 0.0, 0.0], [0.0, size, 0.0]],
                    vec![[0, 1, 2]],
                ))),
                trajectory,
            )
            .unwrap(),
        ),
        material: Material::GrayDiffuseThermal {
            surface: GrayDiffuseThermalSurface::try_new(temperature, 0.7, ContentHash([42; 32]))
                .unwrap(),
        },
        emission: None,
    }
}

fn scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 3.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn time(t: f64) -> PathTime {
    PathTime {
        interval: ShutterInterval::resolve(
            0.0,
            1.0,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap(),
        normalized: NormalizedShutterTime::try_new(t).unwrap(),
        cached_animated: core::array::from_fn(|_| None),
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_rigid_source_density_matches_area_jacobian_at_the_actual_shutter_time() {
    with_cx(|cx, _| {
        let end = RigidTransform::try_new([0.0, 0.5, 0.0, 3.0_f64.sqrt() / 2.0], [2.0, -1.0, 0.5])
            .unwrap();
        let scene = scene(vec![source(1, 1.0, end, 2200.0)]);
        let lights = MovingThermalLights::try_new(&scene, cx).unwrap();
        let local_origin = Point3::new(0.15, 0.25, 3.0);
        for t in [0.0, 0.37, 1.0] {
            let mut time = time(t);
            let transform = lights.sources[0]
                .instance
                .trajectory()
                .evaluate(t)
                .unwrap()
                .transform;
            let origin = transform.transform_point(local_origin);
            let sample = lights
                .sample(cx, Some(&time), origin, [0.25, 0.36, 0.7])
                .unwrap()
                .unwrap();
            let displacement = sample.point.delta_from(origin);
            let distance = displacement.norm();
            // A unit right triangle has area 1/2. Proper rigid motion leaves
            // dA/domega = r² / |n.w| unchanged at every evaluated instant.
            let expected = 2.0 * distance * distance
                / sample.normal.dot(displacement.scale(1.0 / distance)).abs();
            near(sample.pdf_solid_angle, expected);
            near(
                lights
                    .hit_pdf(cx, Some(&time), 0, origin, sample.point)
                    .unwrap(),
                expected,
            );
            // Geometry tracing's cached pose is the identical transform used
            // when this source is queried outside the four-instance cache.
            time.cached_animated[0] = Some(CachedAnimatedInstance {
                primitive_index: 0,
                instance: lights.sources[0].instance.instance_at(cx, t).unwrap(),
            });
            let cached = lights
                .sample(cx, Some(&time), origin, [0.25, 0.36, 0.7])
                .unwrap()
                .unwrap();
            assert_eq!(cached.point, sample.point);
            assert_eq!(cached.normal, sample.normal);
            assert_eq!(cached.pdf_solid_angle, sample.pdf_solid_angle);
        }
    });
}

#[test]
fn g3_body_power_mixture_matches_reverse_density_and_primitive_reordering() {
    with_cx(|cx, _| {
        let make = || {
            vec![
                source(8, 1.0, RigidTransform::identity(), 2200.0),
                source(
                    2,
                    2.0,
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [3.0, 0.0, 0.0]).unwrap(),
                    2200.0,
                ),
            ]
        };
        let original = scene(make());
        let mut reversed_sources = make();
        reversed_sources.reverse();
        let reversed = scene(reversed_sources);
        let a = MovingThermalLights::try_new(&original, cx).unwrap();
        let b = MovingThermalLights::try_new(&reversed, cx).unwrap();
        let time = time(1.0);
        let origin = Point3::new(0.15, 0.25, 3.0);
        // Sources are ordered by object ID. The larger body's area is four
        // times the smaller's, before the explicit uniform-support mixture.
        near(a.probability(0), 1.0 / 32.0 + (15.0 / 16.0) * (4.0 / 5.0));
        let mut count_large = 0;
        for k in 0..4096 {
            let u = [(f64::from(k) + 0.5) / 4096.0, 0.36, 0.7];
            let x = a.sample(cx, Some(&time), origin, u).unwrap().unwrap();
            let y = b.sample(cx, Some(&time), origin, u).unwrap().unwrap();
            count_large += usize::from(x.primitive_index == 1);
            assert_eq!(x.point, y.point);
            assert_eq!(x.pdf_solid_angle, y.pdf_solid_angle);
            assert_eq!(x.primitive_index, 1 - y.primitive_index);
            near(
                a.hit_pdf(cx, Some(&time), x.primitive_index, origin, x.point)
                    .unwrap(),
                x.pdf_solid_angle,
            );
        }
        assert_eq!(count_large, 3200);
        assert_eq!(a.hit_pdf(cx, None, 91, origin, origin).unwrap(), 0.0);
    });
}

#[test]
fn g3_world_roundoff_retains_forward_reverse_support_for_rotated_distant_meshes() {
    with_cx(|cx, _| {
        let inverse_norm = 1.0 / 30.0_f64.sqrt();
        let end = RigidTransform::try_new(
            [
                inverse_norm,
                2.0 * inverse_norm,
                3.0 * inverse_norm,
                4.0 * inverse_norm,
            ],
            [1.0e9, -3.0e9, 2.0e9],
        )
        .unwrap();
        let mut lost_without_world_uncertainty = 0;
        // Exercise both the exact small-source loop and the spatial hierarchy.
        for faces in [1, 12] {
            let mut primitive = source(1, 1.0, end, 2200.0);
            let Shape::AnimatedInstance(instance) = &primitive.shape else {
                unreachable!()
            };
            let mut vertices = Vec::new();
            let mut triangles = Vec::new();
            for face in 0..faces {
                let x = 2.0 * f64::from(face);
                vertices.extend([[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]]);
                let first = 3 * face as u32;
                triangles.push([first, first + 1, first + 2]);
            }
            primitive.shape = Shape::AnimatedInstance(
                AnimatedGeometryInstance::try_new(
                    1,
                    ContentHash([41; 32]),
                    SharedGeometry::Mesh(Arc::new(TriMesh::new(vertices, triangles))),
                    instance.trajectory().clone(),
                )
                .unwrap(),
            );
            let scene = scene(vec![primitive]);
            let lights = MovingThermalLights::try_new(&scene, cx).unwrap();
            for t in [0.2, 0.65, 1.0] {
                let time = time(t);
                let transform = lights.sources[0].transform(Some(&time)).unwrap();
                let inverse = transform.inverse().unwrap();
                let origin = transform.transform_point(Point3::new(0.15, 0.25, 3.0));
                for u in [[0.23, 0.0, 0.0], [0.23, 0.36, 0.7], [0.23, 0.81, 0.999999]] {
                    let sample = lights.sample(cx, Some(&time), origin, u).unwrap().unwrap();
                    assert!(sample.pdf_solid_angle > 0.0);
                    assert_eq!(
                        lights
                            .hit_pdf(cx, Some(&time), 0, origin, sample.point)
                            .unwrap(),
                        sample.pdf_solid_angle,
                    );
                    let strict = lights.sources[0].local.thermal_mixture_pdf(
                        0,
                        inverse.transform_point(origin),
                        inverse.transform_point(sample.point),
                    );
                    lost_without_world_uncertainty += usize::from(strict == 0.0);
                }
                // The envelope accounts for arithmetic error, not a physical
                // expansion of the mesh to unrelated local positions.
                let outside = transform.transform_point(Point3::new(-1.0, -1.0, 0.0));
                assert_eq!(
                    lights.hit_pdf(cx, Some(&time), 0, origin, outside).unwrap(),
                    0.0
                );
            }
        }
        assert!(lost_without_world_uncertainty > 0);
    });
}

#[test]
fn g4_empty_cold_invalid_time_and_cancel_have_explicit_outcomes() {
    with_cx(|cx, gate| {
        let cold_scene = scene(vec![source(1, 1.0, RigidTransform::identity(), 1.0)]);
        let cold = MovingThermalLights::try_new(&cold_scene, cx).unwrap();
        let origin = Point3::new(0.0, 0.0, 3.0);
        assert!(cold.is_empty());
        assert!(cold.sample(cx, None, origin, [0.5; 3]).unwrap().is_none());
        let hot_scene = scene(vec![source(1, 1.0, RigidTransform::identity(), 2200.0)]);
        let hot = MovingThermalLights::try_new(&hot_scene, cx).unwrap();
        assert!(matches!(
            hot.sample(cx, None, origin, [0.5; 3]),
            Err(TracerError::MissingRayTime)
        ));
        assert!(matches!(
            hot.sample(cx, Some(&time(0.5)), origin, [1.0, 0.5, 0.5]),
            Err(TracerError::InvalidInput)
        ));
        gate.request();
        assert!(matches!(
            hot.sample(cx, Some(&time(0.5)), origin, [0.5; 3]),
            Err(TracerError::Cancelled)
        ));
        assert!(matches!(
            hot.hit_pdf(cx, Some(&time(0.5)), 0, origin, origin),
            Err(TracerError::Cancelled)
        ));
        assert!(matches!(
            MovingThermalLights::try_new(&hot_scene, cx),
            Err(TracerError::Cancelled)
        ));
    });
}
