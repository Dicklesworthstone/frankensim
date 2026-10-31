//! Analytic triangle witnesses, cached refit and complete swept candidates.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_query::{
    PreparedTriangleSurface, TriangleFeature, TriangleQueryError, closest_triangle_point,
};

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 41,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-13 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn plane_edge_and_vertex_witnesses_match_independent_geometry() {
    let triangle = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let cases = [
        (
            [0.5, 0.25, 3.0],
            [0.5, 0.25, 0.0],
            [0.625, 0.25, 0.125],
            3.0,
            TriangleFeature::Face,
        ),
        (
            [1.0, -3.0, 2.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            13.0_f64.sqrt(),
            TriangleFeature::Edge(0, 1),
        ),
        (
            [-1.0, -2.0, 1.0],
            [0.0; 3],
            [1.0, 0.0, 0.0],
            6.0_f64.sqrt(),
            TriangleFeature::Vertex(0),
        ),
        (
            [2.0, 2.0, 1.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.5, 0.5],
            3.0_f64.sqrt(),
            TriangleFeature::Edge(1, 2),
        ),
    ];
    for (query, position, barycentric, distance, feature) in cases {
        let p = closest_triangle_point(query, triangle).unwrap();
        for axis in 0..3 {
            near(p.point[axis], position[axis]);
            near(p.barycentric[axis], barycentric[axis]);
        }
        near(p.distance, distance);
        assert_eq!(p.feature, feature);
        assert_eq!(p.unit_normal, [0.0, 0.0, 1.0]);
    }
    let reversed =
        closest_triangle_point([0.5, 0.25, 3.0], [triangle[0], triangle[2], triangle[1]]).unwrap();
    assert_eq!(reversed.unit_normal, [0.0, 0.0, -1.0]);
    assert_eq!(reversed.barycentric, [0.625, 0.125, 0.25]);
}

#[test]
fn triangle_distance_avoids_squaring_physical_scales() {
    for scale in [1e-120, 1.0, 1e120] {
        let p = closest_triangle_point(
            [0.25 * scale, 0.25 * scale, 2.0 * scale],
            [[0.0; 3], [scale, 0.0, 0.0], [0.0, scale, 0.0]],
        )
        .unwrap();
        near(p.distance / scale, 2.0);
        for (actual, expected) in p.barycentric.into_iter().zip([0.5, 0.25, 0.25]) {
            near(actual, expected);
        }
    }
    assert!(matches!(
        closest_triangle_point([0.0; 3], [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
        Err(TriangleQueryError::DegenerateTriangle)
    ));
    assert!(closest_triangle_point([f64::NAN, 0.0, 0.0], [[0.0; 3]; 3]).is_err());
}

fn separated_triangles() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut positions = Vec::new();
    let mut triangles = Vec::new();
    for i in 0..32u32 {
        let x = 8.0 * f64::from(i);
        positions.extend([[x, 0.0, 0.0], [x + 2.0, 0.0, 0.0], [x, 2.0, 0.0]]);
        triangles.push([3 * i, 3 * i + 1, 3 * i + 2]);
    }
    (positions, triangles)
}

#[test]
fn cached_partition_refits_large_deformation_without_losing_nearest_or_nearby_faces() {
    with_cx(|_, cx| {
        let (positions, triangles) = separated_triangles();
        let prepared = PreparedTriangleSurface::new(&positions, &triangles, cx).unwrap();
        let refit = prepared.refit(&positions, cx).unwrap();
        let first = refit
            .point_query([136.5, 0.5, 3.0], 3.0, 1000, 32, cx)
            .unwrap();
        assert_eq!(first.nearest_triangle, 17);
        assert_eq!(
            first.nearby.iter().map(|(i, _)| *i).collect::<Vec<_>>(),
            vec![17]
        );
        near(first.nearest.distance, 3.0);
        assert!(first.tested_triangles < triangles.len());

        // Swap two separated material triangles after the spatial partition was
        // prepared. Triangle IDs and barycentrics stay material, not spatial.
        let mut moved = positions.clone();
        for local in 0..3 {
            moved.swap(local, 3 * 17 + local);
        }
        let refit = prepared.refit(&moved, cx).unwrap();
        let after = refit
            .point_query([136.5, 0.5, 3.0], 3.0, 1000, 32, cx)
            .unwrap();
        assert_eq!(after.nearest_triangle, 0);
        assert_eq!(
            after.nearby.iter().map(|(i, _)| *i).collect::<Vec<_>>(),
            vec![0]
        );
        near(after.nearest.distance, 3.0);
        assert_eq!(after.nearest.barycentric, [0.5, 0.25, 0.25]);
        assert_eq!(
            after.nearby,
            refit
                .point_query([136.5, 0.5, 3.0], 3.0, 1000, 32, cx)
                .unwrap()
                .nearby
        );
    });
}

#[test]
fn swept_bvh_retains_fast_through_shot_and_moving_triangles() {
    with_cx(|_, cx| {
        let (start, triangles) = separated_triangles();
        let mut end = start.clone();
        // Triangle 10 sweeps across the query line even though neither endpoint
        // box overlaps its centre. Triangle 17 is stationary on that line.
        for p in &mut end[30..33] {
            p[0] += 112.0;
        }
        let prepared = PreparedTriangleSurface::new(&start, &triangles, cx).unwrap();
        let swept = prepared.refit_swept(&start, &end, cx).unwrap();
        let candidates = swept
            .sphere_candidates([136.5, 0.5, 10.0], [136.5, 0.5, -10.0], 0.125, 1000, 32, cx)
            .unwrap();
        assert_eq!(candidates.triangles, vec![10, 17]);
        assert!(candidates.tested_triangles < triangles.len());
        let empty = swept
            .sphere_candidates(
                [136.5, 20.0, 10.0],
                [136.5, 20.0, -10.0],
                0.125,
                1000,
                0,
                cx,
            )
            .unwrap();
        assert!(empty.triangles.is_empty());
    });
}

#[test]
fn bounded_query_refuses_incomplete_output_or_cancelled_refit() {
    with_cx(|gate, cx| {
        let (positions, triangles) = separated_triangles();
        let prepared = PreparedTriangleSurface::new(&positions, &triangles, cx).unwrap();
        let refit = prepared.refit(&positions, cx).unwrap();
        assert!(matches!(
            refit.point_query([0.5, 0.5, 1.0], 1.0, 0, 32, cx),
            Err(TriangleQueryError::WorkBudgetExhausted { .. })
        ));
        assert!(matches!(
            refit.point_query([0.5, 0.5, 1.0], 1.0, 1000, 0, cx),
            Err(TriangleQueryError::CandidateBudgetExhausted { .. })
        ));
        let swept = prepared.refit_swept(&positions, &positions, cx).unwrap();
        assert!(matches!(
            swept.sphere_candidates([0.5, 0.5, 1.0], [0.5, 0.5, -1.0], 0.125, 1000, 0, cx),
            Err(TriangleQueryError::CandidateBudgetExhausted { .. })
        ));
        gate.request();
        assert!(matches!(
            prepared.refit(&positions, cx),
            Err(TriangleQueryError::Cancelled)
        ));
    });
}
