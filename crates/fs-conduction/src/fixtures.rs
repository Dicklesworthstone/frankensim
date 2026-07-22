//! Deterministic mesh fixtures for the conduction batteries.
//!
//! Everything here is built from ONE primitive: a structured `nx × ny ×
//! nz` index grid, each cell split into the 6 Kuhn/Freudenthal
//! tetrahedra, with vertex POSITIONS supplied by a caller mapping from
//! the index grid. Because the diagonals are chosen by index order (not
//! geometry), the subdivision stays conforming under any injective,
//! orientation-preserving mapping — which is how the same routine
//! produces the unit cube, a rectangular fin, and a curved annular
//! sector.
//!
//! Generation is pure combinatorics plus the caller's mapping: no RNG,
//! no floating-point branching, identical across runs. Trigonometry in
//! [`annulus_sector`] routes through `fs_math::det`, so even the curved
//! fixture is bit-identical across ISAs.

use fs_rep_mesh::TetComplex;

/// The 6 axis orderings; each gives one Kuhn tet
/// `0 → e_{p0} → e_{p0}+e_{p1} → 1`.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// Kuhn-subdivide an `nx × ny × nz` index grid, positioning vertex
/// `(i, j, k)` at `place(i, j, k)`.
///
/// Vertex `(i, j, k)` gets index `i·(ny+1)(nz+1) + j·(nz+1) + k`.
///
/// # Panics
/// If any count is zero.
#[must_use]
pub fn structured_tets(
    counts: [usize; 3],
    place: &dyn Fn(usize, usize, usize) -> [f64; 3],
) -> (TetComplex, Vec<[f64; 3]>) {
    let [nx, ny, nz] = counts;
    assert!(
        nx > 0 && ny > 0 && nz > 0,
        "structured_tets needs at least one cell per axis"
    );
    let (px, py, pz) = (nx + 1, ny + 1, nz + 1);
    let idx = |i: usize, j: usize, k: usize| -> u32 {
        u32::try_from(i * py * pz + j * pz + k).expect("grid fits u32")
    };
    let mut positions = Vec::with_capacity(px * py * pz);
    for i in 0..px {
        for j in 0..py {
            for k in 0..pz {
                positions.push(place(i, j, k));
            }
        }
    }
    let mut tets = Vec::with_capacity(6 * nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                for perm in &PERMS {
                    let mut corners = [[0usize; 3]; 4];
                    corners[0] = [i, j, k];
                    for (step, &axis) in perm.iter().enumerate() {
                        corners[step + 1] = corners[step];
                        corners[step + 1][axis] += 1;
                    }
                    let v: Vec<u32> = corners.iter().map(|c| idx(c[0], c[1], c[2])).collect();
                    let odd = matches!(perm, [0, 2, 1] | [1, 0, 2] | [2, 1, 0]);
                    tets.push(if odd {
                        [v[0], v[2], v[1], v[3]]
                    } else {
                        [v[0], v[1], v[2], v[3]]
                    });
                }
            }
        }
    }
    (TetComplex::from_tets(positions.len(), tets), positions)
}

/// An axis-aligned box `[0, ex] × [0, ey] × [0, ez]` on an `nx × ny × nz`
/// Kuhn grid.
///
/// # Panics
/// If any count is zero.
#[must_use]
pub fn box_grid(counts: [usize; 3], extent: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
    let [nx, ny, nz] = counts;
    let hx = extent[0] / nx as f64;
    let hy = extent[1] / ny as f64;
    let hz = extent[2] / nz as f64;
    structured_tets(counts, &|i, j, k| {
        [i as f64 * hx, j as f64 * hy, k as f64 * hz]
    })
}

/// The unit cube on an `n × n × n` Kuhn grid — the G1 ladder's base
/// fixture. Positions are reproduced as `i/n` so a face coordinate is
/// exactly `0.0` or (for most `n`) within one ULP of `1.0`; see
/// [`on_box_face`] for the classification band this requires.
///
/// # Panics
/// If `n` is zero.
#[must_use]
pub fn unit_cube(n: usize) -> (TetComplex, Vec<[f64; 3]>) {
    box_grid([n, n, n], [1.0, 1.0, 1.0])
}

/// Face classification for box fixtures. A face coordinate reconstructed
/// as `n · (extent/n)` can sit one ULP off `extent`, so an exact compare
/// misclassifies boundary vertices as interior; interior nodes are at
/// least one cell away, so an absolute band of `1e-9` is unambiguous.
#[must_use]
pub fn on_box_face(coordinate: f64, target: f64) -> bool {
    (coordinate - target).abs() < 1e-9
}

/// An annular sector: `r ∈ [r_inner, r_outer]`, `θ ∈ [0, sweep]`,
/// `z ∈ [0, height]`, on an `nr × ntheta × nz` Kuhn grid.
///
/// This is the curved-geometry fixture for the analytic cylindrical
/// case. The radial faces (`r = r_inner`, `r = r_outer`) carry the
/// Dirichlet data; the `θ` and `z` end faces are adiabatic, which is
/// EXACT for a radially symmetric solution because `∇T · n = 0` there.
///
/// # Panics
/// If any count is zero, or the radii/sweep are not positive and
/// ordered.
#[must_use]
pub fn annulus_sector(
    counts: [usize; 3],
    r_inner: f64,
    r_outer: f64,
    sweep: f64,
    height: f64,
) -> (TetComplex, Vec<[f64; 3]>) {
    assert!(
        r_inner > 0.0 && r_outer > r_inner && sweep > 0.0 && height > 0.0,
        "annulus_sector needs 0 < r_inner < r_outer and positive sweep/height"
    );
    let [nr, nt, nz] = counts;
    let dr = (r_outer - r_inner) / nr as f64;
    let dt = sweep / nt as f64;
    let dz = height / nz as f64;
    structured_tets(counts, &|i, j, k| {
        let r = r_inner + i as f64 * dr;
        let theta = j as f64 * dt;
        [
            r * fs_math::det::cos(theta),
            r * fs_math::det::sin(theta),
            k as f64 * dz,
        ]
    })
}

/// Radius of a point in the `xy` plane — the classifier the annular
/// fixture's boundary tagging uses.
#[must_use]
pub fn cylindrical_radius(p: [f64; 3]) -> f64 {
    fs_math::det::sqrt(p[0].mul_add(p[0], p[1] * p[1]))
}
