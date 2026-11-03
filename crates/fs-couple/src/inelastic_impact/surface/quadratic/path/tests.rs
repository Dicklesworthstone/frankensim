use super::*;
use crate::inelastic_impact::surface::FiniteSolidSurfaceControls;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::frictional_multicontact::NormalVelocityLaw;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_mbd::{Pose, UnitQuaternion, Vec3};
use fs_query::triangle::{
    quadratic::QuadraticQueryControls,
    sweep::{LinearSweepControls, linear_sphere_triangle},
};
use fs_solid::linear3::{TetAssemblyBudget, tet10::promote_p1_to_p2};

fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 19,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn fixture(cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, QuadraticSolidContactSurface) {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let (nodes, elements) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, cx).unwrap();
    let surface = QuadraticSolidContactSurface::from_reference(
        &nodes,
        &elements,
        TetAssemblyBudget::standard(),
        4,
        cx,
    )
    .unwrap();
    // x=X-8tYZ, y=Y, z=Z: a globally invertible quadratic shear with det F=1
    // for every time. All corner coordinates remain unchanged.
    let end: Vec<_> = nodes
        .iter()
        .map(|p| [p[0] - 8.0 * p[1] * p[2], p[1], p[2]])
        .collect();
    assert_eq!(&nodes[..4], &end[..4]);
    (nodes, end, surface)
}

fn rigid(center: [f64; 3]) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            Vec3::new(center[0], center[1], center[2]),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap()
}

fn sphere(id: u64, body: usize) -> FiniteSolidSphere {
    FiniteSolidSphere {
        id,
        rigid_body: body,
        center_body_m: Vec3::ZERO,
        radius_m: 1.0 / 128.0,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn controls() -> QuadraticSolidPathControls {
    QuadraticSolidPathControls {
        initial: QuadraticSolidSurfaceControls {
            surface: FiniteSolidSurfaceControls {
                activation_tolerance_m: 1e-4,
                maximum_contacts: 32,
                ..Default::default()
            },
            query: QuadraticQueryControls {
                distance_tolerance: 1e-5,
                witness_tolerance: 1e-4,
                maximum_cells: 1 << 18,
                maximum_active_cells: 1 << 14,
            },
            maximum_witness_span_m: 5e-4,
        },
        sweep: QuadraticSweepControls {
            time: LinearSweepControls {
                maximum_windows: 2048,
                fraction_tolerance: 1e-3,
            },
            distance: QuadraticQueryControls {
                distance_tolerance: 1e-4,
                witness_tolerance: 1.0,
                maximum_cells: 32768,
                maximum_active_cells: 8192,
            },
            maximum_distance_cells: 1 << 20,
        },
    }
}

#[test]
fn fixed_corners_cannot_hide_an_interior_curved_encounter_or_its_earliest_sphere() {
    run(|_, cx| {
        let (start, end, surface) = fixture(cx);
        let centers_start = [[-1.0 / 32.0, 0.25, 0.25], [-1.0 / 16.0, 0.25, 0.25]];
        let centers_end = [[-1.0 / 32.0, -0.25, 0.25], [-1.0 / 16.0, -0.25, 0.25]];
        let rigid_start = centers_start.map(rigid);
        let rigid_end = centers_end.map(rigid);
        // Canonical sphere 17 is examined first but meets the face later.
        // Sphere 99 must win even though input and temporal orders differ.
        let spheres = [sphere(99, 0), sphere(17, 1)];
        let endpoint = surface
            .sphere_separations_at_positions(&end, &rigid_end, &spheres, controls().initial, cx)
            .unwrap();
        assert!(
            endpoint
                .closest
                .iter()
                .all(|witness| witness.gap_lower_m > 0.0)
        );
        for (index, sphere) in spheres.iter().enumerate() {
            for face in surface.faces() {
                let corners = [face[0], face[1], face[2]].map(|node| start[node as usize]);
                let linear = linear_sphere_triangle(
                    centers_start[index],
                    centers_end[index],
                    sphere.radius_m,
                    corners,
                    corners,
                    controls().sweep.time,
                    cx,
                )
                .unwrap();
                assert!(matches!(linear.verdict, LinearSweepVerdict::Clear { .. }));
            }
        }
        let report = surface
            .affine_sphere_paths_at_positions(
                &start,
                &end,
                &rigid_start,
                &rigid_end,
                &spheres,
                controls(),
                cx,
            )
            .unwrap();
        match report.verdict {
            QuadraticSolidPathVerdict::PossibleContact {
                sphere_id,
                face,
                fraction_start,
                fraction_end,
            } => {
                assert_eq!(sphere_id, 99);
                assert!(fraction_start > 0.0 && fraction_end < 0.125);
                assert!(fraction_end - fraction_start <= controls().sweep.time.fraction_tolerance);
                let nodes = surface.faces()[face];
                let mut corners = [nodes[0], nodes[1], nodes[2]];
                corners.sort_unstable();
                assert_eq!(corners, [0, 2, 3]);
            }
            other => panic!("missed isochoric midside crossing: {other:?}"),
        }
        assert!(report.tested_paths >= 2);
        assert!(report.examined_windows >= report.tested_paths);
        assert!(report.examined_distance_cells <= controls().sweep.maximum_distance_cells);
    });
}

#[test]
fn swept_culling_reports_finite_authored_clearance_after_initial_admission() {
    run(|_, cx| {
        let (start, end, surface) = fixture(cx);
        let report = surface
            .affine_sphere_paths_at_positions(
                &start,
                &end,
                &[rigid([-10.0, 0.25, 0.25])],
                &[rigid([-10.0, -0.25, 0.25])],
                &[sphere(7, 0)],
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            report.verdict,
            QuadraticSolidPathVerdict::Clear {
                minimum_clearance_m: controls().initial.surface.activation_tolerance_m,
            }
        );
        assert_eq!(report.tested_paths, 0);
        assert_eq!(report.examined_windows, 0);
        assert!(report.examined_distance_cells > 0);
    });
}

#[test]
fn exact_initial_plane_touch_can_recede_without_an_epsilon_allowance() {
    run(|_, cx| {
        let (nodes, _, surface) = fixture(cx);
        let mut descriptor = sphere(7, 0);
        descriptor.radius_m = 0.125;
        let end = [rigid([-0.25, 0.25, 0.25])];
        let report = surface
            .affine_sphere_paths_at_positions(
                &nodes,
                &nodes,
                &[rigid([-0.125, 0.25, 0.25])],
                &end,
                &[descriptor],
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(
            report.verdict,
            QuadraticSolidPathVerdict::Clear {
                minimum_clearance_m: 0.0,
            }
        );
        // Both the actual touching face and the opposite face's bounding hull
        // are swept candidates. Each is clear in one window, but their work
        // cannot each receive a fresh copy of the same one-window budget.
        let mut limits = controls();
        limits.sweep.time.maximum_windows = 1;
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(
                &nodes,
                &nodes,
                &[rigid([-0.125, 0.25, 0.25])],
                &end,
                &[descriptor],
                limits,
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Query(
                QuadraticQueryError::Geometry(TriangleQueryError::WorkBudgetExhausted { .. })
            ))
        ));
        let penetrating = surface.affine_sphere_paths_at_positions(
            &nodes,
            &nodes,
            &[rigid([-0.125_f64.next_down(), 0.25, 0.25])],
            &end,
            &[descriptor],
            controls(),
            cx,
        );
        assert!(matches!(
            penetrating,
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::Penetration { sphere_id: 7, .. }
            )) | Ok(QuadraticSolidPathReport {
                verdict: QuadraticSolidPathVerdict::PossibleContact {
                    fraction_start: 0.0,
                    ..
                },
                ..
            })
        ));
    });
}

#[test]
fn initially_contained_or_penetrating_spheres_cannot_be_admitted_by_surface_distance() {
    run(|_, cx| {
        let (start, _, surface) = fixture(cx);
        let contained = rigid([0.125, 0.125, 0.125]);
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(
                &start,
                &start,
                &[contained],
                &[contained],
                &[sphere(7, 0)],
                controls(),
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::UnresolvedExterior { sphere_id: 7 }
            ))
        ));
        let penetrating = rigid([-1.0 / 256.0, 0.25, 0.25]);
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(
                &start,
                &start,
                &[penetrating],
                &[penetrating],
                &[sphere(7, 0)],
                controls(),
                cx,
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::Penetration { sphere_id: 7, .. }
            ))
        ));
    });
}

#[test]
fn aggregate_work_and_cancellation_refuse_without_partial_clearance() {
    run(|gate, cx| {
        let (start, end, surface) = fixture(cx);
        let a = [rigid([-1.0 / 32.0, 0.25, 0.25])];
        let b = [rigid([-1.0 / 32.0, -0.25, 0.25])];
        let spheres = [sphere(7, 0)];
        let initial = surface
            .sphere_separations_at_positions(&start, &a, &spheres, controls().initial, cx)
            .unwrap();
        let mut limits = controls();
        limits.sweep.time.maximum_windows = 1;
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(&start, &end, &a, &b, &spheres, limits, cx),
            Err(QuadraticSolidSurfaceError::Query(
                QuadraticQueryError::Geometry(TriangleQueryError::WorkBudgetExhausted { .. })
            ))
        ));
        limits = controls();
        limits.initial.surface.maximum_query_visits =
            initial.visited_nodes + initial.tested_triangles;
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(&start, &end, &a, &b, &spheres, limits, cx),
            Err(QuadraticSolidSurfaceError::Query(
                QuadraticQueryError::Geometry(TriangleQueryError::WorkBudgetExhausted { .. })
            ))
        ));
        limits = controls();
        limits.sweep.maximum_distance_cells = initial.examined_cells;
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(&start, &end, &a, &b, &spheres, limits, cx),
            Err(QuadraticSolidSurfaceError::Query(
                QuadraticQueryError::CellBudgetExhausted { .. }
            )) | Err(QuadraticSolidSurfaceError::Query(
                QuadraticQueryError::ActiveCellBudgetExhausted { .. }
            ))
        ));
        gate.request();
        assert!(matches!(
            surface.affine_sphere_paths_at_positions(
                &start,
                &end,
                &a,
                &b,
                &spheres,
                controls(),
                cx
            ),
            Err(QuadraticSolidSurfaceError::Surface(
                FiniteSolidSurfaceError::Cancelled
            ))
        ));
    });
}

#[test]
fn empty_sphere_sets_still_validate_authored_controls_and_all_end_nodes() {
    run(|_, cx| {
        let (start, mut end, surface) = fixture(cx);
        let mut limits = controls();
        limits.sweep.time.maximum_windows = 0;
        assert!(
            surface
                .affine_sphere_paths_at_positions(&start, &end, &[], &[], &[], limits, cx)
                .is_err()
        );
        limits = controls();
        limits.initial.query.maximum_active_cells = limits.initial.query.maximum_cells + 1;
        assert!(
            surface
                .affine_sphere_paths_at_positions(&start, &end, &[], &[], &[], limits, cx)
                .is_err()
        );
        end[4][0] = f64::NAN;
        assert!(
            surface
                .affine_sphere_paths_at_positions(&start, &end, &[], &[], &[], controls(), cx)
                .is_err()
        );
    });
}
