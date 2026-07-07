//! Deterministic tet-mesh fixtures: the Kuhn/Freudenthal 6-tet cube
//! subdivision on structured grids (conforming, covers the unit cube,
//! refinement ladder for G1 convergence studies) plus the minimal
//! single-tet and two-tet complexes. Fixture generation is pure
//! combinatorics — no RNG, no floating-point decisions — so meshes are
//! identical across runs and ISAs.

use fs_rep_mesh::TetComplex;

/// One reference tetrahedron (vertices 0..4).
#[must_use]
pub fn single_tet() -> (TetComplex, Vec<[f64; 3]>) {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    (TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), positions)
}

/// Two tets sharing the face {1, 2, 3}.
#[must_use]
pub fn two_tets() -> (TetComplex, Vec<[f64; 3]>) {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    (
        TetComplex::from_tets(5, vec![[0, 1, 2, 3], [1, 2, 3, 4]]),
        positions,
    )
}

/// Kuhn/Freudenthal subdivision of the unit cube into an n×n×n grid of
/// cells, each split into 6 tets along the sorted-coordinate paths
/// from corner (0,0,0) to corner (1,1,1). Conforming across cell faces
/// (neighbouring cells induce the same diagonal on shared faces), and
/// every tet has POSITIVE volume in stored order.
///
/// Vertex (i, j, k) has index `i·(n+1)² + j·(n+1) + k` and position
/// (i/n, j/n, k/n).
///
/// # Panics
/// If `n == 0`.
#[must_use]
pub fn kuhn_cube(n: usize) -> (TetComplex, Vec<[f64; 3]>) {
    assert!(n > 0, "kuhn_cube needs at least one cell");
    let np = n + 1;
    let idx = |i: usize, j: usize, k: usize| -> u32 {
        u32::try_from(i * np * np + j * np + k).expect("grid fits u32")
    };
    let h = 1.0 / n as f64;
    let mut positions = Vec::with_capacity(np * np * np);
    for i in 0..np {
        for j in 0..np {
            for k in 0..np {
                positions.push([i as f64 * h, j as f64 * h, k as f64 * h]);
            }
        }
    }
    // The 6 permutations of unit steps (x, y, z): each ordering of the
    // three axes gives one tet 0 → e_{p0} → e_{p0}+e_{p1} → 1.
    const PERMS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut tets = Vec::with_capacity(6 * n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let base = [i, j, k];
                for perm in &PERMS {
                    let mut corners = [[0usize; 3]; 4];
                    corners[0] = base;
                    for (step, &axis) in perm.iter().enumerate() {
                        corners[step + 1] = corners[step];
                        corners[step + 1][axis] += 1;
                    }
                    let v: Vec<u32> = corners
                        .iter()
                        .map(|c| idx(c[0], c[1], c[2]))
                        .collect();
                    // Positive orientation in stored order: the path
                    // tets alternate parity with the permutation sign;
                    // swap the middle pair on odd permutations.
                    let odd = matches!(perm, [0, 2, 1] | [1, 0, 2] | [2, 1, 0]);
                    let tet = if odd {
                        [v[0], v[2], v[1], v[3]]
                    } else {
                        [v[0], v[1], v[2], v[3]]
                    };
                    tets.push(tet);
                }
            }
        }
    }
    (TetComplex::from_tets(np * np * np, tets), positions)
}

/// True when the vertex at `p` lies on the boundary of the unit cube
/// (fixture helper for Dirichlet pinning in tests).
#[must_use]
pub fn on_unit_cube_boundary(p: [f64; 3]) -> bool {
    p.iter().any(|&c| c == 0.0 || c == 1.0)
}
