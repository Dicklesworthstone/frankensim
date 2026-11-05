use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 17,
                    kernel_id: 93,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

const CORNERS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn mapped(
    complex: &TetComplex,
    original: &[[f64; 3]],
    map: impl Fn([f64; 3]) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let positions = original.iter().copied().map(&map).collect();
    let midpoints = complex
        .edges
        .iter()
        .map(|&[a, b]| {
            map(core::array::from_fn(|i| {
                f64::midpoint(original[a as usize][i], original[b as usize][i])
            }))
        })
        .collect();
    (positions, midpoints)
}

fn close(got: f64, expected: f64) {
    assert!(
        (got - expected).abs() <= 3.0e-13 * expected.abs().max(1.0),
        "got {got:e}, expected {expected:e}"
    );
}

fn vector_close(got: [f64; 3], expected: [f64; 3]) {
    for i in 0..3 {
        close(got[i], expected[i]);
    }
}

#[test]
fn isochoric_quadratic_shear_retains_exact_volume_gradient_and_curved_normal() {
    run(|_, cx| {
        // This globally specified map is exactly quadratic and det J = 1.
        // The lambda3=0 face is a curved graph, despite unchanged corners.
        let a = 0.75;
        let map = |[x, y, z]: [f64; 3]| [x, y, z + a * x * y];
        let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, midpoints) = mapped(&complex, &CORNERS, map);
        assert_eq!(positions, CORNERS);
        let geometry = QuadraticGeometry::new(
            &complex,
            &positions,
            &midpoints,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        for lambda in [
            [0.125, 0.25, 0.125, 0.5],
            [0.25, 0.25, 0.25, 0.25],
            [0.25, 0.5, 0.25, 0.0],
            [1.0, 0.0, 0.0, 0.0],
        ] {
            let [_, x, y, z] = lambda;
            let sample = geometry.point(&complex, &positions, 0, lambda).unwrap();
            vector_close(sample.position_m, map([x, y, z]));
            close(sample.determinant_m3, 1.0);
            let gradients = [
                [-1.0 + a * y, -1.0 + a * x, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-a * y, -a * x, 1.0],
            ];
            for i in 0..4 {
                vector_close(sample.barycentric_gradients[i], gradients[i]);
            }
        }
        // All six nodal values survive the canonical-edge to local-edge map.
        for (i, &[u, v]) in EDGES.iter().enumerate() {
            let mut lambda = [0.0; 4];
            lambda[u] = 0.5;
            lambda[v] = 0.5;
            let expected = map(core::array::from_fn(|j| {
                0.5 * (CORNERS[u][j] + CORNERS[v][j])
            }));
            vector_close(
                geometry
                    .point(&complex, &positions, 0, lambda)
                    .unwrap()
                    .position_m,
                expected,
            );
            assert_eq!(
                geometry.elements[0].edges[i],
                complex
                    .edges
                    .binary_search(&{
                        let mut key = [u as u32, v as u32];
                        key.sort_unstable();
                        key
                    })
                    .unwrap()
            );
        }
        let lambda = [0.25, 0.5, 0.25, 0.0];
        let face = geometry
            .face_point(&complex, &positions, 0, 3, lambda)
            .unwrap();
        let outward = [a * lambda[2], a * lambda[1], -1.0];
        let length = det::sqrt(dot(outward, outward));
        close(face.surface_jacobian_m2, length);
        vector_close(face.outward_normal, outward.map(|v| v / length));
        assert!(!geometry.face_is_affine(&complex, &positions, 0, 3).unwrap());
        assert!(!geometry.face_is_affine(&complex, &positions, 0, 0).unwrap());
        assert!(geometry.face_is_affine(&complex, &positions, 0, 1).unwrap());
        assert!(geometry.face_is_affine(&complex, &positions, 0, 2).unwrap());
    });
}

#[test]
fn affine_midsides_recover_original_metrics_for_all_stored_corner_orientations() {
    run(|_, cx| {
        let map = |[x, y, z]: [f64; 3]| {
            [
                1024.0 + 2.0 * x + 0.25 * y,
                -2048.0 + 0.125 * x + y + 0.5 * z,
                4096.0 - 0.25 * y + 1.5 * z,
            ]
        };
        let weights = [0.125, 0.25, 0.25, 0.375];
        let expected_point = map([weights[1], weights[2], weights[3]]);
        let mut count = 0;
        for a in 0..4_u32 {
            for b in 0..4_u32 {
                for c in 0..4_u32 {
                    for d in 0..4_u32 {
                        let tet = [a, b, c, d];
                        let mut sorted = tet;
                        sorted.sort_unstable();
                        if sorted != [0, 1, 2, 3] {
                            continue;
                        }
                        count += 1;
                        let complex = TetComplex::from_tets(4, vec![tet]);
                        let (positions, midpoints) = mapped(&complex, &CORNERS, map);
                        let geometry = QuadraticGeometry::new(
                            &complex,
                            &positions,
                            &midpoints,
                            QuadraticGeometryControls::default(),
                            cx,
                        )
                        .unwrap();
                        let metrics =
                            crate::whitney::try_element_geometry(&complex, &positions).unwrap();
                        let lambda = tet.map(|v| weights[v as usize]);
                        let original =
                            affine_point(&complex, &positions, 0, lambda, &metrics).unwrap();
                        assert_eq!(original.barycentric_gradients, metrics.grads[0]);
                        let curved = geometry.point(&complex, &positions, 0, lambda).unwrap();
                        vector_close(curved.position_m, expected_point);
                        close(curved.determinant_m3, 6.0 * metrics.vol_signed[0].abs());
                        for i in 0..4 {
                            vector_close(curved.barycentric_gradients[i], metrics.grads[0][i]);
                            assert!(geometry.face_is_affine(&complex, &positions, 0, i).unwrap());
                        }
                    }
                }
            }
        }
        assert_eq!(count, 24);
    });
}

#[test]
fn a_shared_curved_face_has_identical_positions_measures_and_opposite_normals() {
    run(|_, cx| {
        let original = [
            CORNERS[0],
            CORNERS[1],
            CORNERS[2],
            CORNERS[3],
            [0.0, 0.0, -1.0],
        ];
        let complex = TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]);
        let (positions, midpoints) =
            mapped(&complex, &original, |[x, y, z]| [x, y, z + 0.75 * x * y]);
        let geometry = QuadraticGeometry::new(
            &complex,
            &positions,
            &midpoints,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        for [x, y] in [[0.125, 0.25], [0.25, 0.5], [0.5, 0.0], [0.0, 0.0]] {
            let a = geometry
                .face_point(&complex, &positions, 0, 3, [1.0 - x - y, x, y, 0.0])
                .unwrap();
            let b = geometry
                .face_point(&complex, &positions, 1, 3, [1.0 - x - y, y, x, 0.0])
                .unwrap();
            vector_close(a.point.position_m, b.point.position_m);
            vector_close(a.outward_normal, b.outward_normal.map(|v| -v));
            close(a.surface_jacobian_m2, b.surface_jacobian_m2);
        }
    });
}

#[test]
fn whole_element_proof_refuses_folds_exhaustion_invalid_data_and_cancellation() {
    run(|gate, cx| {
        let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, folded) = mapped(&complex, &CORNERS, |[x, y, z]| {
            [x + 1.5 * x * (1.0 - x), y, z]
        });
        assert_eq!(positions, CORNERS);
        // The center has positive determinant 1.75 and the corner tetrahedron
        // is unchanged, but det J = 2.5 - 3x is negative next to x=1.
        assert!(matches!(
            QuadraticGeometry::new(
                &complex,
                &positions,
                &folded,
                QuadraticGeometryControls::default(),
                cx
            ),
            Err(AcousticFeError::CurvedGeometry {
                tet: 0,
                reason: QuadraticGeometryRefusal::Unresolved(_)
            })
        ));

        // det J = (x-1/2)^2 + e^2 + e*y is positive on the entire simplex.
        // Initial interval/Bernstein bounds straddle zero, so a one-cell
        // budget must refuse instead of promoting sampled positivity.
        let e = 0.125;
        let (positions, midpoints) = mapped(&complex, &CORNERS, |[x, y, z]| {
            [0.5 * x * x - 0.5 * x - e * y, e * x + (x - 0.5) * y, z]
        });
        let controls = QuadraticGeometryControls {
            maximum_cells: 1,
            ..Default::default()
        };
        assert!(matches!(
            QuadraticGeometry::new(&complex, &positions, &midpoints, controls, cx),
            Err(AcousticFeError::CurvedGeometry {
                tet: 0,
                reason: QuadraticGeometryRefusal::Query(
                    QuadraticTetPathError::CellBudgetExhausted { maximum_cells: 1 }
                )
            })
        ));
        assert!(
            QuadraticGeometry::new(
                &complex,
                &positions,
                &midpoints,
                QuadraticGeometryControls::default(),
                cx
            )
            .is_ok()
        );

        assert!(
            QuadraticGeometry::new(
                &complex,
                &positions,
                &midpoints[..5],
                QuadraticGeometryControls::default(),
                cx
            )
            .is_err()
        );
        let mut nonfinite = midpoints.clone();
        nonfinite[0][0] = f64::NAN;
        assert!(
            QuadraticGeometry::new(
                &complex,
                &positions,
                &nonfinite,
                QuadraticGeometryControls::default(),
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            QuadraticGeometry::new(
                &complex,
                &positions,
                &midpoints,
                QuadraticGeometryControls::default(),
                cx
            ),
            Err(AcousticFeError::Cancelled)
        ));
    });
}
