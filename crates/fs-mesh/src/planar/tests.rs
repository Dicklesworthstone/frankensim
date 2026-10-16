use super::*;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 55,
                kernel_id: 55,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}
fn options() -> PlanarOptions<'static> {
    PlanarOptions {
        maximum_edge_m: 2.0,
        minimum_angle_deg: 20.0,
        size_field: None,
        budget: PlanarBudget {
            max_vertices: 10000,
            max_triangles: 20000,
            max_refinement_steps: 10000,
        },
    }
}
const SQUARE: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

#[test]
fn g4_cancelled_outline_refuses_before_work() {
    let gate = fs_exec::CancelGate::new_clock_free();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 55,
                kernel_id: 55,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        assert!(matches!(
            mesh_outline(
                &PlanarOutline {
                    outer: Loop::Polyline(&SQUARE),
                    holes: vec![]
                },
                &options(),
                &cx
            ),
            Err(PlanarError::Cancelled)
        ));
    });
}

#[test]
fn g0_uniform_h_refinement_and_caps() {
    with_cx(|cx| {
        let original = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&SQUARE),
                holes: vec![],
            },
            &options(),
            cx,
        )
        .unwrap();
        let budget = PlanarBudget {
            max_vertices: 9,
            max_triangles: 8,
            max_refinement_steps: 5,
        };
        let fine = original.refine_uniform(budget, cx).unwrap();
        assert_eq!(fine.vertices.len(), 9);
        assert_eq!(fine.triangles.len(), 8);
        assert_eq!(fine.receipt.boundary_loops[0].len(), 8);
        assert!((fine.receipt.mesh_area_m2 - 1.0).abs() < 1e-14);
        assert!(matches!(
            original.refine_uniform(
                PlanarBudget {
                    max_vertices: 8,
                    ..budget
                },
                cx
            ),
            Err(PlanarError::Budget {
                what: "vertices",
                ..
            })
        ));
        assert!(matches!(
            original.refine_uniform(
                PlanarBudget {
                    max_triangles: 7,
                    ..budget
                },
                cx
            ),
            Err(PlanarError::Budget {
                what: "triangles",
                ..
            })
        ));
        assert!(matches!(
            original.refine_uniform(
                PlanarBudget {
                    max_refinement_steps: 4,
                    ..budget
                },
                cx
            ),
            Err(PlanarError::Budget {
                what: "refinement steps",
                ..
            })
        ));
    });
}

#[test]
fn g0_exact_predicates_on_collinear_and_cocircular_input() {
    with_cx(|cx| {
        for ring in [
            SQUARE.to_vec(),
            vec![[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        ] {
            let outline = PlanarOutline {
                outer: Loop::Polyline(&ring),
                holes: vec![],
            };
            let a = mesh_outline(&outline, &options(), cx).unwrap();
            let b = mesh_outline(&outline, &options(), cx).unwrap();
            assert_eq!(a.vertices, b.vertices);
            assert_eq!(a.triangles, b.triangles);
            assert_eq!(a.vertices.len(), ring.len());
            assert!((a.receipt.mesh_area_m2 - 1.0).abs() < 1e-14);
        }
    });
}
#[test]
fn g0_hole_seed_outside_refuses() {
    with_cx(|cx| {
        let hole = [[2.0, 2.0], [3.0, 2.0], [2.0, 3.0]];
        assert!(matches!(
            mesh_outline(
                &PlanarOutline {
                    outer: Loop::Polyline(&SQUARE),
                    holes: vec![Loop::Polyline(&hole)]
                },
                &options(),
                cx
            ),
            Err(PlanarError::HoleOutside { .. })
        ));
    });
}
#[test]
fn g0_overlapping_holes_refuse() {
    with_cx(|cx| {
        let a = [[0.1, 0.1], [0.6, 0.1], [0.6, 0.6], [0.1, 0.6]];
        let b = [[0.4, 0.4], [0.9, 0.4], [0.9, 0.9], [0.4, 0.9]];
        assert!(matches!(
            mesh_outline(
                &PlanarOutline {
                    outer: Loop::Polyline(&SQUARE),
                    holes: vec![Loop::Polyline(&a), Loop::Polyline(&b)]
                },
                &options(),
                cx
            ),
            Err(PlanarError::OverlappingHoles { .. })
        ));
    });
}
#[test]
fn g0_self_intersecting_outline_refuses() {
    with_cx(|cx| {
        let bow = [[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]];
        assert!(matches!(
            mesh_outline(
                &PlanarOutline {
                    outer: Loop::Polyline(&bow),
                    holes: vec![]
                },
                &options(),
                cx
            ),
            Err(PlanarError::SelfIntersection { .. })
        ));
        let invalid = [[0.0, 0.0], [1.0, f64::NAN], [0.0, 1.0]];
        assert!(
            mesh_outline(
                &PlanarOutline {
                    outer: Loop::Polyline(&invalid),
                    holes: vec![]
                },
                &options(),
                cx
            )
            .is_err()
        );
    });
}
#[test]
fn g0_budget_cap_and_cap_plus_one() {
    with_cx(|cx| {
        let outline = PlanarOutline {
            outer: Loop::Polyline(&SQUARE),
            holes: vec![],
        };
        let mut opts = options();
        opts.budget = PlanarBudget {
            max_vertices: 4,
            max_triangles: 2,
            max_refinement_steps: 0,
        };
        assert!(mesh_outline(&outline, &opts, cx).is_ok());
        opts.budget.max_vertices = 3;
        assert!(matches!(
            mesh_outline(&outline, &opts, cx),
            Err(PlanarError::Budget { .. })
        ));
        opts.budget.max_vertices = 4;
        opts.budget.max_triangles = 1;
        assert!(matches!(
            mesh_outline(&outline, &opts, cx),
            Err(PlanarError::Budget {
                what: "triangles",
                ..
            })
        ));
        opts = options();
        opts.maximum_edge_m = 0.25;
        opts.budget.max_refinement_steps = 0;
        assert!(matches!(
            mesh_outline(&outline, &opts, cx),
            Err(PlanarError::Budget {
                what: "refinement steps",
                ..
            })
        ));
    });
}
#[test]
fn g3_input_rotation_and_translation_invariance() {
    with_cx(|cx| {
        let transformed: Vec<_> = SQUARE.iter().map(|p| [4.0 - p[1], p[0] - 2.0]).collect();
        let a = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&SQUARE),
                holes: vec![],
            },
            &options(),
            cx,
        )
        .unwrap();
        let b = mesh_outline(
            &PlanarOutline {
                outer: Loop::Polyline(&transformed),
                holes: vec![],
            },
            &options(),
            cx,
        )
        .unwrap();
        assert_eq!(a.triangles.len(), b.triangles.len());
        assert!((a.receipt.mesh_area_m2 - b.receipt.mesh_area_m2).abs() < 1e-12);
        assert!((a.receipt.minimum_angle_deg - b.receipt.minimum_angle_deg).abs() < 1e-12);
    });
}
