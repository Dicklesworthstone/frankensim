//! Deterministic tet-mesh fixtures: the Kuhn/Freudenthal 6-tet cube
//! subdivision on structured grids (conforming, covers the unit cube,
//! refinement ladder for G1 convergence studies) plus the minimal
//! single-tet and two-tet complexes. Fixture generation is pure
//! combinatorics — no RNG, no floating-point decisions — so meshes are
//! identical across runs and ISAs.

use fs_rep_mesh::TetComplex;

/// The 6 permutations of unit steps (x, y, z): each axis ordering
/// gives one Kuhn tet 0 → e_{p0} → e_{p0}+e_{p1} → 1.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

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
                    let v: Vec<u32> = corners.iter().map(|c| idx(c[0], c[1], c[2])).collect();
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
    // A face coordinate is 0 or 1, but `kuhn_cube(n)` reconstructs the far face
    // as `n * (1.0/n)`, which for many n is 1 ULP below 1.0 (first at n = 49),
    // so an EXACT bit compare misclassifies those face vertices as INTERIOR
    // (and, under Dirichlet pinning, yields a near-singular reduced system).
    // Interior vertices sit ≥ 1/n from a face — far outside this tolerance for
    // any realistic mesh — so a small absolute band is unambiguous.
    p.iter()
        .any(|&c| c.abs() < 1e-9 || (c - 1.0).abs() < 1e-9)
}
/// A masked cube grid: `nx × ny × nz` cells, keeping only cells where
/// `keep(i, j, k)` — the MULTIPLY-CONNECTED fixture builder (rings,
/// slabs with holes, exterior-like domains). Unused vertices are
/// compacted away so Betti numbers reflect the kept region only.
/// Positions are in cell units (vertex (i,j,k) at (i, j, k)).
///
/// # Panics
/// If no cell is kept.
#[must_use]
pub fn masked_cube_grid(
    nx: usize,
    ny: usize,
    nz: usize,
    keep: &dyn Fn(usize, usize, usize) -> bool,
) -> (TetComplex, Vec<[f64; 3]>) {
    let (npx, npy, npz) = (nx + 1, ny + 1, nz + 1);
    let raw = |i: usize, j: usize, k: usize| i * npy * npz + j * npz + k;
    let mut used = vec![false; npx * npy * npz];
    let mut cells = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if keep(i, j, k) {
                    cells.push([i, j, k]);
                    for di in 0..2 {
                        for dj in 0..2 {
                            for dk in 0..2 {
                                used[raw(i + di, j + dj, k + dk)] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(!cells.is_empty(), "masked_cube_grid kept no cells");
    let mut remap = vec![u32::MAX; used.len()];
    let mut positions = Vec::new();
    let mut next = 0u32;
    for i in 0..npx {
        for j in 0..npy {
            for k in 0..npz {
                let r = raw(i, j, k);
                if used[r] {
                    remap[r] = next;
                    next += 1;
                    positions.push([i as f64, j as f64, k as f64]);
                }
            }
        }
    }
    let mut tets = Vec::with_capacity(6 * cells.len());
    for base in cells {
        for perm in &PERMS {
            let mut corners = [[0usize; 3]; 4];
            corners[0] = base;
            for (step, &axis) in perm.iter().enumerate() {
                corners[step + 1] = corners[step];
                corners[step + 1][axis] += 1;
            }
            let v: Vec<u32> = corners
                .iter()
                .map(|c| remap[raw(c[0], c[1], c[2])])
                .collect();
            let odd = matches!(perm, [0, 2, 1] | [1, 0, 2] | [2, 1, 0]);
            let tet = if odd {
                [v[0], v[2], v[1], v[3]]
            } else {
                [v[0], v[1], v[2], v[3]]
            };
            tets.push(tet);
        }
    }
    (TetComplex::from_tets(positions.len(), tets), positions)
}
