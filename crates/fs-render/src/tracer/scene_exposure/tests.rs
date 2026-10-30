//! Accepted-state selection, preparation bounds and numeric continuity.

use super::*;
use crate::animated_instances::{RigidTransformTrajectory, TransformKeyframe};
use crate::instances::RigidTransform;
use crate::motion::{ShutterConvention, ShutterDistribution};
use crate::spectral::lift_rgb;
use crate::thermal_surface::GrayDiffuseThermalSurface;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 43,
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

fn shutter(open: f64, close: f64) -> ShutterInterval {
    ShutterInterval::try_from_canonical_parts(
        open,
        close,
        ShutterConvention::FrontLoaded,
        ShutterDistribution::UniformCounterV1,
    )
    .unwrap()
}

fn budget() -> SceneExposureBudget {
    SceneExposureBudget {
        max_frames: 8,
        max_total_primitives: 16,
        max_total_light_triangles: 32,
    }
}

fn empty_scene() -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn mesh(height: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [-1.0, -1.0, height],
            [1.0, -1.0, height],
            [1.0, 1.0, height],
            [-1.0, 1.0, height],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn thermal(temperature: f64, emissivity: f64) -> Material {
    Material::GrayDiffuseThermal {
        surface: GrayDiffuseThermalSurface::try_new(temperature, emissivity, ContentHash([43; 32]))
            .unwrap(),
    }
}

fn thermal_scene(temperature: f64, emissivity: f64, height: f64) -> Scene {
    let mut scene = empty_scene();
    scene.primitives.push(Primitive {
        shape: Shape::Mesh(mesh(height)),
        material: thermal(temperature, emissivity),
        emission: None,
    });
    scene
}

fn moving_scene(start: f64, end: f64) -> Scene {
    let trajectory = RigidTransformTrajectory::try_new(vec![
        TransformKeyframe::try_new(
            start,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0; 3]).unwrap(),
            [0.0; 3],
        )
        .unwrap(),
        TransformKeyframe::try_new(
            end,
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.1, 0.0, 0.0]).unwrap(),
            [0.0; 3],
        )
        .unwrap(),
    ])
    .unwrap();
    let mut scene = empty_scene();
    scene.primitives.push(Primitive {
        shape: Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                43,
                ContentHash([44; 32]),
                SharedGeometry::mesh(mesh(0.0)),
                trajectory,
            )
            .unwrap(),
        ),
        material: thermal(1800.0, 0.4),
        emission: None,
    });
    scene
}

#[test]
fn g0_exact_state_boundaries_select_later_scenes_and_include_final_time() {
    with_cx(|cx, _| {
        let scenes = [empty_scene(), empty_scene(), empty_scene()];
        let frames: [_; 3] = core::array::from_fn(|index| SceneFrame {
            start_s: index as f64,
            end_s: index as f64 + 1.0,
            scene: &scenes[index],
        });
        // Closing exactly on a change still prepares that later state's
        // proposals, despite its zero-duration intersection with the shutter.
        let exposure = SceneExposure::try_new(cx, &frames, shutter(0.5, 2.0), budget()).unwrap();
        assert_eq!(exposure.frame_count(), 3);
        assert!(exposure.lighting.iter().all(Option::is_some));
        assert!(core::ptr::eq(exposure.initial_scene(), &scenes[0]));
        for (time, index) in [(0.5, 0), (1.0, 1), (1.75, 1), (2.0, 2)] {
            assert!(core::ptr::eq(
                exposure.at_time(time).unwrap().0,
                &scenes[index]
            ));
        }
        for time in [0.0, 2.5, f64::NAN, f64::INFINITY] {
            assert_eq!(
                exposure.at_time(time).err(),
                Some(SceneExposureError::TimeOutsideExposure.into())
            );
        }
        let final_point = SceneExposure::try_new(cx, &frames, shutter(3.0, 3.0), budget()).unwrap();
        assert!(final_point.lighting[0].is_none());
        assert!(final_point.lighting[1].is_none());
        assert!(final_point.lighting[2].is_some());
        assert!(core::ptr::eq(
            final_point.at_time(3.0).unwrap().0,
            &scenes[2]
        ));
    });
}

#[test]
fn g0_ambiguous_or_uncovered_physical_time_is_refused() {
    with_cx(|cx, _| {
        let scene = empty_scene();
        let frame = |start_s, end_s| SceneFrame {
            start_s,
            end_s,
            scene: &scene,
        };
        let reject = |frames: &[SceneFrame<'_>]| {
            SceneExposure::try_new(cx, frames, shutter(0.0, 1.0), budget()).err()
        };
        assert_eq!(reject(&[]), Some(SceneExposureError::Empty.into()));
        for (start, end) in [
            (0.0, 0.0),
            (1.0, 0.0),
            (f64::NAN, 1.0),
            (0.0, f64::INFINITY),
        ] {
            assert_eq!(
                reject(&[frame(start, end)]),
                Some(SceneExposureError::InvalidInterval { index: 0 }.into())
            );
        }
        assert_eq!(
            reject(&[frame(0.0, 0.4), frame(0.5, 1.0)]),
            Some(SceneExposureError::Gap { index: 1 }.into())
        );
        assert_eq!(
            reject(&[frame(0.0, 0.6), frame(0.5, 1.0)]),
            Some(SceneExposureError::Overlap { index: 1 }.into())
        );
        assert_eq!(
            reject(&[frame(0.1, 1.0)]),
            Some(SceneExposureError::ShutterOutsideHistory.into())
        );
        assert_eq!(
            reject(&[frame(0.0, 0.9)]),
            Some(SceneExposureError::ShutterOutsideHistory.into())
        );
    });
}

#[test]
fn g0_budget_counts_static_and_moving_thermal_faces_and_rectangle_metadata() {
    with_cx(|cx, _| {
        let mut scene = moving_scene(0.0, 1.0);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(mesh(0.0)),
            material: thermal(1800.0, 0.4),
            emission: None,
        });
        let emission = (lift_rgb([1.0; 3]), 2.0);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(mesh(0.0)),
            material: Material::Lambertian {
                reflectance: lift_rgb([0.0; 3]),
            },
            emission: Some(emission),
        });
        scene.lights.push(RectLight {
            corner: Point3::new(-1.0, -1.0, 0.0),
            edge_u: Vec3::new(2.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 2.0, 0.0),
            prim: 2,
            emission,
        });
        let frames = [SceneFrame {
            start_s: 0.0,
            end_s: 1.0,
            scene: &scene,
        }];
        let exact = SceneExposureBudget {
            max_frames: 1,
            max_total_primitives: 3,
            max_total_light_triangles: 6,
        };
        SceneExposure::try_new(cx, &frames, shutter(0.0, 1.0), exact).unwrap();
        for insufficient in [
            SceneExposureBudget {
                max_frames: 0,
                ..exact
            },
            SceneExposureBudget {
                max_total_primitives: 2,
                ..exact
            },
            SceneExposureBudget {
                max_total_light_triangles: 5,
                ..exact
            },
        ] {
            assert_eq!(
                SceneExposure::try_new(cx, &frames, shutter(0.0, 1.0), insufficient).err(),
                Some(SceneExposureError::Budget.into())
            );
        }

        // A zero-emissivity accepted material allocates no thermal proposal.
        let nonemitting = thermal_scene(1800.0, 0.0, 0.0);
        let frames = [SceneFrame {
            start_s: 0.0,
            end_s: 1.0,
            scene: &nonemitting,
        }];
        SceneExposure::try_new(
            cx,
            &frames,
            shutter(0.0, 1.0),
            SceneExposureBudget {
                max_total_light_triangles: 0,
                ..exact
            },
        )
        .unwrap();
    });
}

#[test]
fn g0_each_trajectory_is_admitted_only_over_its_selected_frame_intersection() {
    with_cx(|cx, _| {
        let scenes = [
            moving_scene(0.0, 1.0),
            moving_scene(1.0, 2.0),
            moving_scene(9.0, 10.0),
        ];
        let frames: [_; 3] = core::array::from_fn(|index| SceneFrame {
            start_s: index as f64,
            end_s: index as f64 + 1.0,
            scene: &scenes[index],
        });
        // Neither selected trajectory covers the entire exposure. Each covers
        // exactly its own accepted support, and the unused scene is not prepared.
        let exposure = SceneExposure::try_new(cx, &frames, shutter(0.25, 1.75), budget()).unwrap();
        assert!(exposure.lighting[0].is_some());
        assert!(exposure.lighting[1].is_some());
        assert!(exposure.lighting[2].is_none());
        assert!(core::ptr::eq(exposure.at_time(1.0).unwrap().0, &scenes[1]));
        let at_change = SceneExposure::try_new(cx, &frames, shutter(1.0, 1.0), budget()).unwrap();
        assert!(at_change.lighting[0].is_none());
        assert!(at_change.lighting[1].is_some());
        assert!(at_change.lighting[2].is_none());
    });
}

#[test]
fn g3_history_guard_binds_numeric_material_geometry_camera_and_interval_state() {
    with_cx(|cx, _| {
        let identity = |scene: &Scene, start_s: f64| {
            let frames = [SceneFrame {
                start_s,
                end_s: 1.0,
                scene,
            }];
            SceneExposure::try_new(cx, &frames, shutter(0.25, 0.75), budget())
                .unwrap()
                .history_identity()
        };
        let original = thermal_scene(1800.0, 0.4, 0.0);
        let same_values = thermal_scene(1800.0, 0.4, 0.0);
        let expected = identity(&original, 0.0);
        assert_eq!(expected, identity(&same_values, -0.0));
        assert_ne!(expected, identity(&original, 0.125));
        let mut other_camera = thermal_scene(1800.0, 0.4, 0.0);
        other_camera.camera.half_tan = 0.25;
        for changed in [
            thermal_scene(1900.0, 0.4, 0.0),
            thermal_scene(1800.0, 0.5, 0.0),
            thermal_scene(1800.0, 0.4, 0.1),
            other_camera,
        ] {
            assert_ne!(expected, identity(&changed, 0.0));
        }
    });
}

#[test]
fn g4_cancelled_admission_publishes_no_prepared_exposure() {
    with_cx(|cx, gate| {
        let scene = thermal_scene(1800.0, 0.4, 0.0);
        let frames = [SceneFrame {
            start_s: 0.0,
            end_s: 1.0,
            scene: &scene,
        }];
        gate.request();
        assert_eq!(
            SceneExposure::try_new(cx, &frames, shutter(0.0, 1.0), budget()).err(),
            Some(TracerError::Cancelled)
        );
    });
}
