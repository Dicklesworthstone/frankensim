//! BDM1 element machinery: the minimal H(div)-conforming pair with a
//! FULL-RANK gradient. Per cell, BDM1 = all of P1² (6 dofs: two
//! normal-component moments per edge against q ∈ {1, s}, with s the
//! signed arclength from the GLOBAL edge's lower vertex — orientation
//! consistency for free). div(BDM1) = P0 exactly, so discrete
//! velocities are EXACTLY divergence-free against P0 pressures — the
//! de Rham exactness that buys pressure-robustness. (RT0 was the first
//! candidate and is rank-deficient for the viscous form: its per-cell
//! gradient is c·I — measured and documented in the bead trail.)
//!
//! The per-cell basis is built numerically: a 6×6 dof-Vandermonde over
//! the monomial basis {(1,0),(x,0),(y,0),(0,1),(0,x),(0,y)}, inverted
//! once per cell (fs-la LU) — no hand-derived formulas to get wrong.

use crate::trimesh::TriMesh;
use fs_la::factor::lu;

/// Per-cell BDM1 basis: coefficients of each of the 6 basis functions
/// in the monomial basis, plus cached gradients (constant per cell).
pub struct CellBasis {
    /// `coef[i]` = 6 monomial coefficients of basis function i.
    pub coef: [[f64; 6]; 6],
    /// Constant gradient of each basis function: [[du_x/dx, du_x/dy],
    /// [du_y/dx, du_y/dy]].
    pub grad: [[[f64; 2]; 2]; 6],
    /// Constant divergence of each basis function.
    pub div: [f64; 6],
}

/// Value of a monomial-coefficient field at a point.
fn eval_mono(c: &[f64; 6], x: [f64; 2]) -> [f64; 2] {
    [
        c[0] + c[1] * x[0] + c[2] * x[1],
        c[3] + c[4] * x[0] + c[5] * x[1],
    ]
}

/// Two-point Gauss on an edge (exact for the P1 moments).
pub fn edge_gauss_pub(a: [f64; 2], b: [f64; 2]) -> [([f64; 2], f64); 2] {
    edge_gauss(a, b)
}

fn edge_gauss(a: [f64; 2], b: [f64; 2]) -> [([f64; 2], f64); 2] {
    let g = 0.5 / 3.0f64.sqrt();
    let mid = [f64::midpoint(a[0], b[0]), f64::midpoint(a[1], b[1])];
    let half = [0.5 * (b[0] - a[0]), 0.5 * (b[1] - a[1])];
    let len = (b[0] - a[0]).hypot(b[1] - a[1]);
    [
        (
            [mid[0] - 2.0 * g * half[0], mid[1] - 2.0 * g * half[1]],
            0.5 * len,
        ),
        (
            [mid[0] + 2.0 * g * half[0], mid[1] + 2.0 * g * half[1]],
            0.5 * len,
        ),
    ]
}

/// Build the BDM1 basis for triangle `t`: dof (2k) is the mean normal
/// moment on local edge k, dof (2k+1) the linear (signed-s) moment —
/// both against the GLOBAL edge normal and orientation.
///
/// # Panics
/// On degenerate cells (singular Vandermonde).
#[must_use]
pub fn cell_basis(mesh: &TriMesh, t: usize) -> CellBasis {
    let mut vand = [[0.0f64; 6]; 6];
    for k in 0..3 {
        let (e, _) = mesh.tri_edges[t][k];
        let edge = &mesh.edges[e];
        let (va, vb) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
        let n = edge.normal;
        for (m, mono) in MONOS.iter().enumerate() {
            let mut m0 = 0.0;
            let mut m1 = 0.0;
            for (gx, w) in edge_gauss(va, vb) {
                let v = mono(gx);
                let un = v[0] * n[0] + v[1] * n[1];
                // Signed arclength coordinate from the lower vertex,
                // centered: s ∈ [−1/2, 1/2].
                let sl = ((gx[0] - va[0]) * (vb[0] - va[0]) + (gx[1] - va[1]) * (vb[1] - va[1]))
                    / (edge.len * edge.len)
                    - 0.5;
                m0 += w * un / edge.len;
                m1 += w * un * sl / edge.len;
            }
            vand[2 * k][m] = m0;
            vand[2 * k + 1][m] = m1;
        }
    }
    // Invert: coef = V⁻¹ (columns = basis functions in monomials).
    let flat: Vec<f64> = vand.iter().flatten().copied().collect();
    let f = lu(&flat, 6).expect("BDM1 Vandermonde nonsingular");
    let mut coef = [[0.0f64; 6]; 6];
    for i in 0..6 {
        let mut rhs = [0.0f64; 6];
        rhs[i] = 1.0;
        let mut x = rhs.to_vec();
        f.solve(&mut x);
        for m in 0..6 {
            coef[i][m] = x[m];
        }
    }
    let mut grad = [[[0.0f64; 2]; 2]; 6];
    let mut div = [0.0f64; 6];
    for i in 0..6 {
        let c = coef[i];
        grad[i] = [[c[1], c[2]], [c[4], c[5]]];
        div[i] = c[1] + c[5];
    }
    CellBasis { coef, grad, div }
}

/// Monomial basis (component-wise P1).
type Mono = fn([f64; 2]) -> [f64; 2];
const MONOS: [Mono; 6] = [
    |_x| [1.0, 0.0],
    |x| [x[0], 0.0],
    |x| [x[1], 0.0],
    |_x| [0.0, 1.0],
    |x| [0.0, x[0]],
    |x| [0.0, x[1]],
];

/// Evaluate basis function `i` of a built cell basis at `x`.
#[must_use]
pub fn eval_basis(basis: &CellBasis, i: usize, x: [f64; 2]) -> [f64; 2] {
    eval_mono(&basis.coef[i], x)
}

/// Degree-4 six-point triangle quadrature (barycentric), exact for
/// quartics — enough for every product this crate assembles.
#[must_use]
pub fn tri_quad(p: [[f64; 2]; 3], area: f64) -> [([f64; 2], f64); 6] {
    const A1: f64 = 0.445_948_490_915_965;
    const A2: f64 = 0.091_576_213_509_771;
    const W1: f64 = 0.223_381_589_678_011;
    const W2: f64 = 0.109_951_743_655_322;
    let bary = [
        [1.0 - 2.0 * A1, A1, A1, W1],
        [A1, 1.0 - 2.0 * A1, A1, W1],
        [A1, A1, 1.0 - 2.0 * A1, W1],
        [1.0 - 2.0 * A2, A2, A2, W2],
        [A2, 1.0 - 2.0 * A2, A2, W2],
        [A2, A2, 1.0 - 2.0 * A2, W2],
    ];
    core::array::from_fn(|q| {
        let b = bary[q];
        (
            [
                b[0] * p[0][0] + b[1] * p[1][0] + b[2] * p[2][0],
                b[0] * p[0][1] + b[1] * p[1][1] + b[2] * p[2][1],
            ],
            b[3] * area,
        )
    })
}
