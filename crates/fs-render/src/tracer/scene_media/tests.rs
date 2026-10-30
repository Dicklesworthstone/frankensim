//! Ownership, bounded admission and immutable scene/field continuity.

use super::*;
use crate::animated_instances::{RigidTransformTrajectory, TransformKeyframe};
use crate::instances::RigidTransform;
use crate::spectral::lift_rgb;
use crate::thermal_surface::GrayDiffuseThermalSurface;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
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
                    seed: 47,
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

fn budget() -> SceneMediaBudget {
    SceneMediaBudget {
        max_interior_bindings: 4,
        max_primitives: 4,
        max_light_triangles: 32,
    }
}

fn field<'a>(cx: &Cx<'_>, density: &'a [f64]) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([1; 3], density, [-1.0; 3], [2.0; 3]),
        MassAbsorptionSpectrum::constant(0.4).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: 1,
            max_steps_per_segment: 16,
        },
    )
    .unwrap()
}

fn cube() -> TriMesh {
    TriMesh::new(
        vec![
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    )
}

fn scene() -> Scene {
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(cube()),
            material: Material::Dielectric {
                glass: DielectricGlass::representative_crown(),
                surface: DielectricSurface::SMOOTH,
            },
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

#[test]
fn g0_bound_vacuum_retains_replacement_ownership_and_borrows_the_original_field() {
    with_cx(|cx, _| {
        let density = [0.0];
        let vacuum = field(cx, &density);
        let scene = scene();
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&vacuum),
        }];
        let media = SceneMedia::try_new(
            cx,
            &scene,
            Some(PhysicalMedium::Absorbing(&vacuum)),
            &bindings,
            budget(),
        )
        .unwrap();
        assert!(media.ambient_path().is_none());
        assert!(
            matches!(media.interior_path(0), Some(PathMedium::Absorbing(grid)) if core::ptr::eq(grid, &vacuum))
        );
        assert!(media.interior_path(1).is_none());
        assert!(core::ptr::eq(media.interiors().as_ptr(), bindings.as_ptr()));
        assert!(core::ptr::eq(media.scene(), &scene));
        assert!(!media.has_thermal_emission());
        let unbound = SceneMedia::try_new(cx, &scene, None, &[], budget()).unwrap();
        assert_ne!(media.content_identity(), unbound.content_identity());
    });
}

#[test]
fn g0_ambiguous_or_missing_medium_owners_are_refused() {
    with_cx(|cx, _| {
        let density = [1.0];
        let field = field(cx, &density);
        let scene = scene();
        let binding = InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&field),
        };
        assert_eq!(
            SceneMedia::try_new(cx, &scene, None, &[binding, binding], budget()).err(),
            Some(SceneMediaError::UnsortedBindings { index: 1 }.into())
        );
        assert_eq!(
            SceneMedia::try_new(
                cx,
                &scene,
                None,
                &[InteriorMedium {
                    boundary_primitive: 1,
                    ..binding
                }],
                budget()
            )
            .err(),
            Some(
                SceneMediaError::MissingBoundary {
                    boundary_primitive: 1
                }
                .into()
            )
        );
        let mut opaque = scene;
        opaque.primitives[0].material = Material::Lambertian {
            reflectance: lift_rgb([0.4; 3]),
        };
        assert_eq!(
            SceneMedia::try_new(cx, &opaque, None, &[binding], budget()).err(),
            Some(
                SceneMediaError::NotDielectric {
                    boundary_primitive: 0
                }
                .into()
            )
        );
    });
}

#[test]
fn g0_metadata_budgets_and_static_world_field_support_are_admitted_before_render() {
    with_cx(|cx, _| {
        let density = [1.0];
        let field = field(cx, &density);
        let mut scene = scene();
        let bindings = [InteriorMedium {
            boundary_primitive: 0,
            medium: PhysicalMedium::Absorbing(&field),
        }];
        for budget in [
            SceneMediaBudget {
                max_interior_bindings: 0,
                ..budget()
            },
            SceneMediaBudget {
                max_primitives: 0,
                ..budget()
            },
        ] {
            assert_eq!(
                SceneMedia::try_new(cx, &scene, None, &bindings, budget).err(),
                Some(SceneMediaError::Budget.into())
            );
        }
        scene.primitives[0].material = Material::GrayDiffuseThermal {
            surface: GrayDiffuseThermalSurface::try_new(1800.0, 0.4, ContentHash([47; 32]))
                .unwrap(),
        };
        assert_eq!(
            SceneMedia::try_new(
                cx,
                &scene,
                None,
                &[],
                SceneMediaBudget {
                    max_light_triangles: 11,
                    ..budget()
                }
            )
            .err(),
            Some(SceneMediaError::Budget.into())
        );
        let pose = RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0; 3]).unwrap();
        let trajectory = RigidTransformTrajectory::try_new(vec![
            TransformKeyframe::try_new(0.0, pose, [0.0; 3]).unwrap(),
            TransformKeyframe::try_new(1.0, pose, [0.0; 3]).unwrap(),
        ])
        .unwrap();
        scene.primitives[0].shape = Shape::AnimatedInstance(
            AnimatedGeometryInstance::try_new(
                47,
                ContentHash([48; 32]),
                SharedGeometry::mesh(cube()),
                trajectory,
            )
            .unwrap(),
        );
        assert_eq!(
            SceneMedia::try_new(cx, &scene, None, &[], budget()).err(),
            Some(SceneMediaError::AnimatedGeometry { primitive: 0 }.into())
        );
    });
}

#[test]
fn g3_numeric_identity_binds_actual_scene_and_field_values_without_addresses() {
    with_cx(|cx, _| {
        let density = [1.0];
        let same_density = [1.0];
        let changed_density = [1.1];
        let field = field(cx, &density);
        let same_field = self::field(cx, &same_density);
        let changed_field = self::field(cx, &changed_density);
        let scene = scene();
        let mut changed_scene = self::scene();
        changed_scene.camera.half_tan = 0.2;
        let identity = |scene: &Scene, field: &AbsorbingGrid<'_>| {
            let bindings = [InteriorMedium {
                boundary_primitive: 0,
                medium: PhysicalMedium::Absorbing(field),
            }];
            SceneMedia::try_new(cx, scene, None, &bindings, budget())
                .unwrap()
                .content_identity()
        };
        let original = identity(&scene, &field);
        assert_eq!(original, identity(&self::scene(), &same_field));
        assert_ne!(original, identity(&scene, &changed_field));
        assert_ne!(original, identity(&changed_scene, &field));
    });
}

#[test]
fn g4_cancelled_admission_publishes_no_media_binding() {
    with_cx(|cx, gate| {
        let scene = scene();
        gate.request();
        assert_eq!(
            SceneMedia::try_new(cx, &scene, None, &[], budget()).err(),
            Some(TracerError::Cancelled)
        );
    });
}
