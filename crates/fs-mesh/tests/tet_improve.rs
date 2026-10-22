//! IRP M1.3-W4: the boundary-frozen tet quality pass on hand-built
//! configurations whose answer is known: a valid and a non-convex 2-3 flip,
//! a 3-2 flip, an inverted interior vertex, and the refusals.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::tet_improve::{ImproveError, ImproveOptions, improve_tets};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x1_3B,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn volume(p: &[[f64; 3]], t: [u32; 4]) -> f64 {
    let [a, b, c, d] = t.map(|v| p[v as usize]);
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0]))
        / 6.0
}

fn options() -> ImproveOptions {
    ImproveOptions {
        target_min_dihedral_deg: 15.0,
        max_passes: 8,
    }
}

/// Faces of a tet set that only one tet uses (sorted triples).
fn hull_faces(tets: &[[u32; 4]]) -> Vec<[u32; 3]> {
    let mut count = std::collections::BTreeMap::new();
    for t in tets {
        for f in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut key = f.map(|i| t[i]);
            key.sort_unstable();
            *count.entry(key).or_insert(0) += 1;
        }
    }
    count
        .into_iter()
        .filter(|(_, c)| *c == 1)
        .map(|(f, _)| f)
        .collect()
}

/// Two tets on triangle (0,1,2): apex 3 barely above it (a sliver), apex 4
/// well below. All five points are frozen, so only a flip can help.
fn bipyramid(below: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.3, 0.3, 0.02],
        below,
    ];
    let mut tets = vec![[0, 1, 2, 3], [0, 2, 1, 4]];
    for t in &mut tets {
        if volume(&positions, *t) < 0.0 {
            t.swap(2, 3);
        }
    }
    assert!(tets.iter().all(|t| volume(&positions, *t) > 0.0));
    (positions, tets)
}

#[test]
fn g0_flip23_valid_and_nonconvex_refused() {
    with_cx(|cx| {
        let (mut positions, mut tets) = bipyramid([0.3, 0.3, -1.0]);
        let before: f64 = tets.iter().map(|t| volume(&positions, *t)).sum();
        let boundary = hull_faces(&tets);
        let stats = improve_tets(
            &mut positions,
            &mut tets,
            &[true; 5],
            &boundary,
            1.0,
            options(),
            cx,
        )
        .expect("valid flip");
        assert_eq!((stats.flips_23, stats.smoothed), (1, 0));
        assert_eq!(tets.len(), 3);
        assert!(stats.min_dihedral_after_deg > stats.min_dihedral_before_deg);
        let after: f64 = tets.iter().map(|t| volume(&positions, *t)).sum();
        assert!(((after - before) / before).abs() <= 1e-14);
        assert!(tets.iter().all(|t| volume(&positions, *t) > 0.0));
        assert_eq!(hull_faces(&tets), boundary, "the hull is unchanged");

        // Apex 4 far outside the prism over the triangle: segment 3-4 misses
        // the shared face, the union is not convex there, and every
        // candidate set would overlap. Nothing changes.
        let (mut positions, mut tets) = bipyramid([3.0, 3.0, -0.2]);
        let original = tets.clone();
        let boundary = hull_faces(&tets);
        let stats = improve_tets(
            &mut positions,
            &mut tets,
            &[true; 5],
            &boundary,
            1.0,
            options(),
            cx,
        )
        .expect("refusing to flip is not an error");
        assert_eq!(stats.flips_23, 0);
        assert_eq!(tets, original);
    });
}

#[test]
fn g0_flip32_removes_a_degree_three_edge() {
    with_cx(|cx| {
        // Edge (3,4) pierces triangle (0,1,2) with its ends at the regular
        // apex height: the three ring tets reach 35.3 degrees, the two tets
        // that replace them 70.5 (computed independently).
        let mut positions = vec![
            [1.0, 0.0, 0.0],
            [-0.5, 0.866, 0.0],
            [-0.5, -0.866, 0.0],
            [0.0, 0.0, 1.414],
            [0.0, 0.0, -1.414],
        ];
        let mut tets: Vec<[u32; 4]> = vec![[3, 4, 0, 1], [3, 4, 1, 2], [3, 4, 2, 0]];
        for t in &mut tets {
            if volume(&positions, *t) < 0.0 {
                t.swap(2, 3);
            }
        }
        let before: f64 = tets.iter().map(|t| volume(&positions, *t)).sum();
        let boundary = hull_faces(&tets);
        let stats = improve_tets(
            &mut positions,
            &mut tets,
            &[true; 5],
            &boundary,
            1.0,
            ImproveOptions {
                target_min_dihedral_deg: 45.0,
                max_passes: 8,
            },
            cx,
        )
        .expect("3-2 flip");
        assert!(stats.min_dihedral_before_deg < 36.0 && stats.min_dihedral_after_deg > 70.0);
        assert_eq!(stats.flips_32, 1);
        assert_eq!(tets.len(), 2);
        let after: f64 = tets.iter().map(|t| volume(&positions, *t)).sum();
        assert!(((after - before) / before).abs() <= 1e-14);
        assert_eq!(hull_faces(&tets), boundary);
    });
}

#[test]
fn g0_inverted_interior_vertex_is_untangled_and_boundary_frozen() {
    with_cx(|cx| {
        // Octahedron split into eight tets around its centre (vertex 6).
        let mut positions = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 0.0],
        ];
        let mut tets = Vec::new();
        for &x in &[0u32, 1] {
            for &y in &[2u32, 3] {
                for &z in &[4u32, 5] {
                    let mut t = [6, x, y, z];
                    if volume(&positions, t) < 0.0 {
                        t.swap(2, 3);
                    }
                    tets.push(t);
                }
            }
        }
        let exact = 4.0 / 3.0;
        // Push the centre outside the octahedron's +x face region: several
        // tets invert.
        positions[6] = [1.4, 0.1, 0.05];
        let frozen = [true, true, true, true, true, true, false];
        let boundary = hull_faces(&tets);
        let hull_before: Vec<[f64; 3]> = positions[..6].to_vec();
        let stats = improve_tets(
            &mut positions,
            &mut tets,
            &frozen,
            &boundary,
            1.0,
            options(),
            cx,
        )
        .expect("untangled");
        assert!(stats.inverted_before > 0);
        assert!(stats.smoothed > 0);
        assert!(tets.iter().all(|t| volume(&positions, *t) > 0.0));
        let total: f64 = tets.iter().map(|t| volume(&positions, *t)).sum();
        assert!(((total - exact) / exact).abs() <= 1e-14);
        assert_eq!(
            &positions[..6],
            &hull_before[..],
            "frozen vertices never move"
        );
    });
}

#[test]
fn g0_refusals() {
    with_cx(|cx| {
        let (mut positions, mut tets) = bipyramid([0.3, 0.3, -1.0]);
        let boundary = hull_faces(&tets);
        for (sign, target, frozen_len) in [(0.0, 15.0, 5), (1.0, 61.0, 5), (1.0, 15.0, 4)] {
            let frozen = vec![true; frozen_len];
            let result = improve_tets(
                &mut positions,
                &mut tets,
                &frozen,
                &boundary,
                sign,
                ImproveOptions {
                    target_min_dihedral_deg: target,
                    max_passes: 1,
                },
                cx,
            );
            assert!(matches!(result, Err(ImproveError::InvalidInput(_))));
        }
        // An inverted tet with every vertex frozen cannot be repaired.
        let mut tets = vec![[0, 2, 1, 3]];
        assert!(volume(&positions, tets[0]) < 0.0);
        let result = improve_tets(
            &mut positions,
            &mut tets,
            &[true; 5],
            &[],
            1.0,
            options(),
            cx,
        );
        assert!(matches!(result, Err(ImproveError::Inverted { tet: 0 })));
    });
}
