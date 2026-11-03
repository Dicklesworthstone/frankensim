use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena, StreamKey { seed: 7, kernel_id: 941, tile: 0, iteration: 0 }, Budget::INFINITE, ExecMode::Deterministic))
    })
}

fn graph(f: impl Fn(f64, f64) -> f64) -> [[f64; 3]; 6] {
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5]]
        .map(|[u, v]| [u, v, f(u, v)])
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e} within {tolerance:e}");
}

#[test]
fn signed_quadratic_interpolation_and_actual_normal_match_analytic_graph() {
    let nodes = graph(|u, v| 0.2 * u * v);
    let locations = [CORNERS[0], CORNERS[1], CORNERS[2], [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]];
    for (node, location) in nodes.into_iter().zip(locations) {
        let evaluated = evaluate_quadratic_triangle(nodes, location).unwrap();
        for axis in 0..3 { near(evaluated.point[axis], node[axis], 2.0e-15); }
    }
    let evaluated = evaluate_quadratic_triangle(nodes, [0.5, 0.25, 0.25]).unwrap();
    near(evaluated.point[2], 0.0125, 2.0e-15);
    let length = 1.005_f64.sqrt();
    for (actual, expected) in evaluated.unit_normal.into_iter().zip([-0.05 / length, -0.05 / length, 1.0 / length]) {
        near(actual, expected, 2.0e-14);
    }
}

#[test]
fn affine_interior_edge_and_vertex_are_globally_localized() {
    with_cx(|_, cx| {
        let nodes = graph(|_, _| 0.0);
        for (query, expected, feature) in [
            ([0.3, 0.2, 1.0], [0.3, 0.2, 0.0], TriangleFeature::Face),
            ([0.5, -0.3, 0.2], [0.5, 0.0, 0.0], TriangleFeature::Edge(0, 1)),
            ([-1.0, -2.0, 1.0], [0.0; 3], TriangleFeature::Vertex(0)),
        ] {
            let result = closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls { witness_tolerance: 1.0e-8, ..Default::default() }, cx).unwrap();
            let distance = norm(core::array::from_fn(|i| query[i] - expected[i])).unwrap();
            assert!(result.distance_lower <= distance + 2.0e-15);
            assert!(result.distance_upper >= distance - 2.0e-15);
            assert!(result.distance_upper - result.distance_lower <= 1.0e-9);
            assert!(result.witness_error_bound <= 1.0e-8);
            assert_eq!(result.feature, feature);
            assert!(result.oriented_side_bounds[0] > 0.0);
            for axis in 0..3 {
                near(result.point[axis], expected[axis], 1.0e-8);
                assert!(result.minimizer_world_bounds[axis][0] <= expected[axis] + 1.0e-14);
                assert!(result.minimizer_world_bounds[axis][1] >= expected[axis] - 1.0e-14);
            }
            assert!(result.examined_cells < 4_096, "affine KKT localization should not tessellate a whole face");
        }
    });
}

#[test]
fn curved_interior_gap_uses_all_six_nodes_and_original_material_coordinates() {
    with_cx(|_, cx| {
        let nodes = graph(|u, v| 0.2 * u * v);
        let expected = [0.25, 0.25, 0.0125];
        let length = 1.005_f64.sqrt();
        let normal = [-0.05 / length, -0.05 / length, 1.0 / length];
        let distance = 0.050_000_1;
        let query = core::array::from_fn(|i| expected[i] + distance * normal[i]);
        let result = closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls { witness_tolerance: 1.0e-5, ..Default::default() }, cx).unwrap();
        for axis in 0..3 { near(result.point[axis], expected[axis], result.witness_error_bound + 1.0e-12); }
        near(result.distance, distance, 2.0e-12);
        assert!(result.distance_lower <= distance + 1.0e-14);
        assert!(result.distance_upper >= distance - 1.0e-14);
        assert!(result.oriented_side_bounds[0] > 0.0);
        assert_eq!(result.feature, TriangleFeature::Face);
        let flat = closest_triangle_point(query, [nodes[0], nodes[1], nodes[2]]).unwrap();
        assert!(flat.distance - result.distance > 0.01, "the corner plane misses the actual bulge");
    });
}

#[test]
fn bowed_boundary_edge_has_its_own_quadratic_closest_point() {
    with_cx(|_, cx| {
        let nodes = graph(|u, v| 0.4 * u * (1.0 - u - v));
        let query = [0.5, -0.2, 0.2];
        let result = closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls::default(), cx).unwrap();
        assert_eq!(result.feature, TriangleFeature::Edge(0, 1));
        for (actual, expected) in result.point.into_iter().zip([0.5, 0.0, 0.1]) { near(actual, expected, 1.0e-8); }
        near(result.distance, 0.05_f64.sqrt(), 2.0e-13);
        assert!(result.unit_normal[1] > 0.19, "normal must include the actual edge slope");
    });
}

#[test]
fn separated_equal_minima_cannot_be_mislabeled_a_localized_witness() {
    with_cx(|_, cx| {
        let nodes = graph(|u, _| 4.0 * (u - 0.5).powi(2));
        let query = [0.5, 0.0, 1.0];
        let distance_only = closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls {
            distance_tolerance: 1.0e-7, witness_tolerance: 10.0, ..Default::default()
        }, cx).unwrap();
        let exact_distance = 0.234_375_f64.sqrt();
        assert!(distance_only.distance_lower <= exact_distance);
        assert!(distance_only.distance_upper >= exact_distance);
        assert!(distance_only.minimizer_bounds[1][0] < 0.04);
        assert!(distance_only.minimizer_bounds[1][1] > 0.96);
        assert!(distance_only.witness_error_bound > 0.9);
        assert!(matches!(closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls {
            distance_tolerance: 1.0e-7, witness_tolerance: 1.0e-4,
            maximum_cells: 4_096, maximum_active_cells: 1_024,
        }, cx), Err(QuadraticQueryError::CellBudgetExhausted { .. }
            | QuadraticQueryError::ActiveCellBudgetExhausted { .. }
            | QuadraticQueryError::ParameterResolution)));
    });
}

#[test]
fn bvh_refit_retains_patch_ids_complete_candidates_and_global_budgets() {
    with_cx(|gate, cx| {
        let mut positions = Vec::new();
        let mut triangles = Vec::new();
        for i in 0..24 {
            let z = if i < 2 { 0.1 * i as f64 } else { 10.0 * i as f64 };
            let start = positions.len() as u32;
            positions.extend(graph(|_, _| z));
            triangles.push(core::array::from_fn(|node| start + node as u32));
        }
        let prepared = PreparedQuadraticTriangleSurface::new(&positions, &triangles, cx).unwrap();
        let current = prepared.refit(&positions, cx).unwrap();
        let query = [0.2, 0.2, 0.4];
        let result = current.point_query(query, 0.5, 256, 4, QuadraticQueryControls::default(), cx).unwrap();
        assert_eq!(result.nearest_triangle, 1);
        assert_eq!(result.nearby.iter().map(|(i, _)| *i).collect::<Vec<_>>(), vec![0, 1]);
        assert!(result.tested_triangles < 12);
        assert!(result.nearest_oriented_side_bounds[0] > 0.0);
        assert!(matches!(current.point_query(query, 0.5, 256, 1, QuadraticQueryControls::default(), cx),
            Err(QuadraticQueryError::Geometry(TriangleQueryError::CandidateBudgetExhausted { .. }))));
        assert!(matches!(current.point_query(query, 0.5, 0, 4, QuadraticQueryControls::default(), cx),
            Err(QuadraticQueryError::Geometry(TriangleQueryError::WorkBudgetExhausted { .. }))));
        for node in &mut positions[..6] { node[2] = 0.39; }
        let refitted = prepared.refit(&positions, cx).unwrap();
        assert_eq!(refitted.point_query(query, 0.5, 256, 4, QuadraticQueryControls::default(), cx).unwrap().nearest_triangle, 0);
        gate.request();
        assert!(matches!(prepared.refit(&positions, cx), Err(QuadraticQueryError::Geometry(TriangleQueryError::Cancelled))));
    });
}

#[test]
fn near_tied_opposite_windings_are_both_present_in_global_side_enclosure() {
    with_cx(|_, cx| {
        let mut positions = graph(|_, _| 0.0).to_vec();
        positions.extend(graph(|_, _| 0.2));
        // Both normals point +z: the lower patch faces the query; the equally
        // near upper patch faces away. One selected upper witness cannot prove
        // that every possible globally nearest patch faces the query.
        let triangles = [[0, 1, 2, 3, 4, 5], [6, 7, 8, 9, 10, 11]];
        let prepared = PreparedQuadraticTriangleSurface::new(&positions, &triangles, cx).unwrap();
        let current = prepared.refit(&positions, cx).unwrap();
        let result = current.point_query([0.25, 0.25, 0.1], 0.01, 32, 2, QuadraticQueryControls::default(), cx).unwrap();
        assert!(result.nearby.is_empty());
        assert!(result.nearest_oriented_side_bounds[0] < 0.0);
        assert!(result.nearest_oriented_side_bounds[1] > 0.0);
    });
}

#[test]
fn swept_bounds_include_midside_motion_and_queries_refuse_exhausted_cells() {
    with_cx(|_, cx| {
        let start = graph(|_, _| 0.0);
        let mut end = start;
        end[3][2] = 1.0;
        let prepared = PreparedQuadraticTriangleSurface::new(&start, &[[0, 1, 2, 3, 4, 5]], cx).unwrap();
        let swept = prepared.refit_swept(&start, &end, cx).unwrap();
        assert_eq!(swept.sphere_candidates([0.5, 0.0, 0.75], [0.5, 0.0, 0.75], 0.1, 10, 1, cx).unwrap().triangles, vec![0]);
        assert!(matches!(closest_quadratic_triangle_point([0.3, 0.2, 1.0], start, QuadraticQueryControls {
            maximum_cells: 1, maximum_active_cells: 1, ..Default::default()
        }, cx), Err(QuadraticQueryError::CellBudgetExhausted { maximum_cells: 1 })));
    });
}

#[test]
fn large_translation_preserves_resolved_small_geometry() {
    with_cx(|_, cx| {
        let shift = [1.0e9, -1.0e9, 1.0e9];
        let nodes = graph(|_, _| 0.0).map(|p| core::array::from_fn(|i| p[i] + shift[i]));
        let query = [shift[0] + 0.25, shift[1] + 0.25, shift[2] + 0.5];
        let result = closest_quadratic_triangle_point(query, nodes, QuadraticQueryControls { witness_tolerance: 1.0e-5, ..Default::default() }, cx).unwrap();
        assert!(result.distance_lower <= 0.5 && result.distance_upper >= 0.5);
        near(result.distance, 0.5, 1.0e-13);
        near(result.point[2], shift[2], 1.0e-5);
    });
}
