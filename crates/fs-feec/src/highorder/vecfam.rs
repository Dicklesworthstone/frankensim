//! Simplicial high-order VECTOR families (bead dcng): hierarchical-
//! by-entity first-kind Nédélec H(curl) and Raviart–Thomas H(div) on
//! tets at r = 1..4, plus discontinuous L² P_{r−1} — completing the
//! simplicial P_rΛᵏ tree whose H¹ member lives in `simplex.rs`.
//!
//! NO hand-derived shape functions (the classic FEEC bug farm): per
//! element the space is spanned by Koszul spanning sets in centered/
//! scaled Cartesian monomials — RT_r = (P_{r−1})³ ⊕ x·H̃_{r−1} (an
//! EXACT direct sum: square dof-Vandermonde, LU) and
//! N_r = (P_{r−1})³ + {x × x^α eᵢ} (overspans by r(r−1)/2: the basis
//! is the least-norm solution c = Aᵀ(AAᵀ)⁻¹e via Cholesky) — and the
//! basis is whatever makes the classical moment dofs Kronecker.
//!
//! Orientation is the sorted-GLOBAL convention extended to FRAMES:
//! edge dofs integrate against Legendre P_k(s) with s signed toward
//! the larger global index (deram1's direction); face dofs use the
//! frame (p_b−p_a, p_c−p_a) and the circulation normal of the SORTED
//! global triple (deram2's normal). Two elements sharing an entity
//! therefore see IDENTICAL functionals — conformity by construction,
//! checked by the battery rather than trusted.

use crate::highorder::quad1d::{gauss_legendre, legendre};
use crate::highorder::simplex::{SimplexSpace, duffy_quadrature};
use fs_la::factor::{cholesky, lu};
use fs_rep_mesh::TetComplex;
use fs_sparse::{Coo, Csr};

/// Highest supported order (dense per-element solves stay small).
pub const MAX_R: usize = 4;

// ------------------------------------------------------------------
// Monomial vector polynomials in centered/scaled local coordinates.
// ------------------------------------------------------------------

/// Monomial exponent list for total degree ≤ d (deterministic order).
#[must_use]
pub fn monomials(d: usize) -> Vec<[usize; 3]> {
    let mut out = Vec::new();
    for total in 0..=d {
        for a in (0..=total).rev() {
            for b in (0..=(total - a)).rev() {
                out.push([a, b, total - a - b]);
            }
        }
    }
    out
}

/// A vector polynomial: 3 components, each a coefficient vector over
/// `monomials(deg)` in LOCAL coordinates ξ = (x − x₀)/h.
#[derive(Clone)]
pub struct VecPoly {
    /// Max total degree of the coefficient tables.
    pub deg: usize,
    /// Component coefficient vectors (length = monomials(deg).len()).
    pub comp: [Vec<f64>; 3],
}

impl VecPoly {
    /// Zero polynomial of degree `deg`.
    #[must_use]
    pub fn zero(deg: usize) -> VecPoly {
        let n = monomials(deg).len();
        VecPoly {
            deg,
            comp: [vec![0.0; n], vec![0.0; n], vec![0.0; n]],
        }
    }

    /// Evaluate at LOCAL point ξ.
    #[must_use]
    pub fn eval_local(&self, monos: &[[usize; 3]], xi: [f64; 3]) -> [f64; 3] {
        let mut v = [0.0f64; 3];
        for (m, &[a, b, c]) in monos.iter().enumerate() {
            let p = xi[0].powi(i32::try_from(a).expect("small"))
                * xi[1].powi(i32::try_from(b).expect("small"))
                * xi[2].powi(i32::try_from(c).expect("small"));
            for k in 0..3 {
                v[k] = self.comp[k][m].mul_add(p, v[k]);
            }
        }
        v
    }

    /// Exact curl in LOCAL coordinates (physical curl = local/h).
    #[must_use]
    pub fn curl_local(&self, monos: &[[usize; 3]]) -> VecPoly {
        let mut out = VecPoly::zero(self.deg);
        let idx = |e: [usize; 3]| -> Option<usize> { monos.iter().position(|&m| m == e) };
        // (curl u)_0 = ∂1 u2 − ∂2 u1, cyclic.
        for (m, &[a, b, c]) in monos.iter().enumerate() {
            let e = [a, b, c];
            let terms: [(usize, usize, usize); 3] = [(0, 1, 2), (1, 2, 0), (2, 0, 1)];
            for &(k, d1, d2) in &terms {
                // +∂_{d1} comp[d2]
                if e[d1] >= 1 {
                    let mut e2 = e;
                    e2[d1] -= 1;
                    if let Some(t) = idx(e2) {
                        out.comp[k][t] += self.comp[d2][m] * e[d1] as f64;
                    }
                }
                // −∂_{d2} comp[d1]
                if e[d2] >= 1 {
                    let mut e2 = e;
                    e2[d2] -= 1;
                    if let Some(t) = idx(e2) {
                        out.comp[k][t] -= self.comp[d1][m] * e[d2] as f64;
                    }
                }
            }
        }
        out
    }

    /// Exact divergence in LOCAL coordinates → scalar coefficients
    /// (physical div = local/h).
    #[must_use]
    pub fn div_local(&self, monos: &[[usize; 3]]) -> Vec<f64> {
        let mut out = vec![0.0f64; monos.len()];
        let idx = |e: [usize; 3]| -> Option<usize> { monos.iter().position(|&m| m == e) };
        for (m, &[a, b, c]) in monos.iter().enumerate() {
            let e = [a, b, c];
            for k in 0..3 {
                if e[k] >= 1 {
                    let mut e2 = e;
                    e2[k] -= 1;
                    if let Some(t) = idx(e2) {
                        out[t] += self.comp[k][m] * e[k] as f64;
                    }
                }
            }
        }
        out
    }
}

// ------------------------------------------------------------------
// Element frames and quadrature.
// ------------------------------------------------------------------

/// Per-element chart: centroid and scale for local coordinates.
#[derive(Clone, Copy)]
pub struct Chart {
    /// Centroid.
    pub x0: [f64; 3],
    /// Isotropic scale (max corner distance from centroid).
    pub h: f64,
}

impl Chart {
    fn of(corners: &[[f64; 3]; 4]) -> Chart {
        let mut x0 = [0.0f64; 3];
        for p in corners {
            for k in 0..3 {
                x0[k] += p[k] / 4.0;
            }
        }
        let mut h = 0.0f64;
        for p in corners {
            let d = fs_math::det::sqrt(
                (p[0] - x0[0]).powi(2) + (p[1] - x0[1]).powi(2) + (p[2] - x0[2]).powi(2),
            );
            h = h.max(d);
        }
        Chart { x0, h }
    }

    /// Map a physical point into local coordinates.
    #[must_use]
    pub fn local(&self, p: [f64; 3]) -> [f64; 3] {
        [
            (p[0] - self.x0[0]) / self.h,
            (p[1] - self.x0[1]) / self.h,
            (p[2] - self.x0[2]) / self.h,
        ]
    }
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -(a[2] * b[1])),
        a[2].mul_add(b[0], -(a[0] * b[2])),
        a[0].mul_add(b[1], -(a[1] * b[0])),
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

/// Collapsed Gauss quadrature on the triangle with vertices `p`:
/// physical points + weights summing to the triangle area.
#[must_use]
pub fn tri_quad3d(p: [[f64; 3]; 3], n: usize) -> Vec<([f64; 3], f64)> {
    let (qx, qw) = gauss_legendre(n);
    let map = |x: f64| f64::midpoint(1.0, x);
    let area = {
        let e1 = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
        let e2 = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
        let c = cross3(e1, e2);
        0.5 * fs_math::det::sqrt(dot3(c, c))
    };
    let mut out = Vec::with_capacity(n * n);
    for (&xu, &wu) in qx.iter().zip(&qw) {
        let u = map(xu);
        for (&xv, &wv) in qx.iter().zip(&qw) {
            let v = map(xv);
            // Barycentrics of the collapse: l1 = u, l2 = v(1−u).
            let (l1, l2) = (u, v * (1.0 - u));
            let l0 = 1.0 - l1 - l2;
            let pt = [
                l0.mul_add(p[0][0], l1.mul_add(p[1][0], l2 * p[2][0])),
                l0.mul_add(p[0][1], l1.mul_add(p[1][1], l2 * p[2][1])),
                l0.mul_add(p[0][2], l1.mul_add(p[1][2], l2 * p[2][2])),
            ];
            // Reference-triangle jacobian (1−u)/4, ref area 1/2.
            out.push((pt, wu * wv * (1.0 - u) / 4.0 * 2.0 * area));
        }
    }
    out
}

// ------------------------------------------------------------------
// Dof functionals (sorted-global frames).
// ------------------------------------------------------------------

/// One moment functional, applied to arbitrary evaluable fields by
/// quadrature (exact for the polynomial degrees in play).
enum Dof {
    /// ∫ u·(p_b−p_a) P_k(2τ−1) dτ along the sorted-global edge.
    Edge { a: [f64; 3], b: [f64; 3], k: usize },
    /// ∫_f (u·dir) q dA / area, q a face monomial in the sorted frame.
    FaceMoment {
        tri: [[f64; 3]; 3],
        dir: [f64; 3],
        qi: usize,
        qj: usize,
    },
    /// ∫_T (u·dir) q dV / vol over the tet, q a cell monomial.
    CellMoment {
        corners: [[f64; 3]; 4],
        dir: [f64; 3],
        qm: [usize; 3],
        chart: Chart,
    },
}

impl Dof {
    /// Apply to a field by quadrature; `n` Gauss points per direction.
    fn apply<F: Fn([f64; 3]) -> [f64; 3]>(&self, u: &F, n: usize) -> f64 {
        match self {
            Dof::Edge { a, b, k } => {
                let (qx, qw) = gauss_legendre(n);
                let t = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let mut acc = 0.0f64;
                for (&xi, &wi) in qx.iter().zip(&qw) {
                    let tau = f64::midpoint(1.0, xi);
                    let p = [
                        tau.mul_add(t[0], a[0]),
                        tau.mul_add(t[1], a[1]),
                        tau.mul_add(t[2], a[2]),
                    ];
                    let (pk, _) = legendre(*k, xi);
                    acc = (wi / 2.0 * pk).mul_add(dot3(u(p), t), acc);
                }
                acc
            }
            Dof::FaceMoment { tri, dir, qi, qj } => {
                let quad = tri_quad3d(*tri, n);
                let area: f64 = quad.iter().map(|&(_, w)| w).sum();
                let mut acc = 0.0f64;
                for (p, w) in quad {
                    // Face barycentrics via the sorted frame: solve for
                    // (s, t) in p = p0 + s(p1−p0) + t(p2−p0).
                    let (s, t) = face_params(*tri, p);
                    let q = (s - 0.5).powi(i32::try_from(*qi).expect("small"))
                        * (t - 0.5).powi(i32::try_from(*qj).expect("small"));
                    acc = (w * q).mul_add(dot3(u(p), *dir), acc);
                }
                acc / area
            }
            Dof::CellMoment {
                corners,
                dir,
                qm,
                chart,
            } => {
                let quad = duffy_quadrature(n);
                let e = [
                    [
                        corners[1][0] - corners[0][0],
                        corners[1][1] - corners[0][1],
                        corners[1][2] - corners[0][2],
                    ],
                    [
                        corners[2][0] - corners[0][0],
                        corners[2][1] - corners[0][1],
                        corners[2][2] - corners[0][2],
                    ],
                    [
                        corners[3][0] - corners[0][0],
                        corners[3][1] - corners[0][1],
                        corners[3][2] - corners[0][2],
                    ],
                ];
                let vol6 = dot3(e[0], cross3(e[1], e[2])).abs();
                let mut acc = 0.0f64;
                let mut wsum = 0.0f64;
                for &(lam, w) in &quad {
                    let mut p = [0.0f64; 3];
                    for (a, corner) in corners.iter().enumerate() {
                        for c in 0..3 {
                            p[c] = lam[a].mul_add(corner[c], p[c]);
                        }
                    }
                    let xi = chart.local(p);
                    let q = xi[0].powi(i32::try_from(qm[0]).expect("small"))
                        * xi[1].powi(i32::try_from(qm[1]).expect("small"))
                        * xi[2].powi(i32::try_from(qm[2]).expect("small"));
                    acc = (w * vol6 * q).mul_add(dot3(u(p), *dir), acc);
                    wsum += w * vol6;
                }
                acc / wsum
            }
        }
    }
}

/// Solve p = p0 + s(p1−p0) + t(p2−p0) in-plane (least squares via the
/// 2×2 Gram system of the frame).
fn face_params(tri: [[f64; 3]; 3], p: [f64; 3]) -> (f64, f64) {
    let e1 = [
        tri[1][0] - tri[0][0],
        tri[1][1] - tri[0][1],
        tri[1][2] - tri[0][2],
    ];
    let e2 = [
        tri[2][0] - tri[0][0],
        tri[2][1] - tri[0][1],
        tri[2][2] - tri[0][2],
    ];
    let d = [p[0] - tri[0][0], p[1] - tri[0][1], p[2] - tri[0][2]];
    let (g11, g12, g22) = (dot3(e1, e1), dot3(e1, e2), dot3(e2, e2));
    let (b1, b2) = (dot3(d, e1), dot3(d, e2));
    let det = g11.mul_add(g22, -(g12 * g12));
    (
        g22.mul_add(b1, -(g12 * b2)) / det,
        g11.mul_add(b2, -(g12 * b1)) / det,
    )
}

// ------------------------------------------------------------------
// Per-entity dof counts.
// ------------------------------------------------------------------

/// Nédélec (first kind) entity dof counts at order r: (edge, face,
/// interior). Global dim = 6r + 4r(r−1) + r(r−1)(r−2)/2 per tet-mesh
/// entity table = r(r+2)(r+3)/2 on one tet.
#[must_use]
pub fn nedelec_entity_dofs(r: usize) -> (usize, usize, usize) {
    (
        r,
        r * r.saturating_sub(1),
        r * r.saturating_sub(1) * r.saturating_sub(2) / 2,
    )
}

/// RT entity dof counts at order r: (face, interior). One-tet dim
/// r(r+1)(r+3)/2.
#[must_use]
pub fn rt_entity_dofs(r: usize) -> (usize, usize) {
    (r * (r + 1) / 2, r * (r + 1) * r.saturating_sub(1) / 2)
}

/// L² (P_{r−1}) dofs per cell.
#[must_use]
pub fn dg_cell_dofs(r: usize) -> usize {
    r * (r + 1) * (r + 2) / 6
}

// ------------------------------------------------------------------
// Element construction.
// ------------------------------------------------------------------

/// A built element basis: each basis function a `VecPoly` in the
/// element chart, ordered to match the space's `element_dofs`.
pub struct ElementBasis {
    /// The chart (local coordinates).
    pub chart: Chart,
    /// Basis functions (dof-dual: functional_i(basis_j) = δᵢⱼ).
    pub funcs: Vec<VecPoly>,
    /// Monomial table for `funcs`.
    pub monos: Vec<[usize; 3]>,
}

/// Family selector.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Family {
    /// First-kind Nédélec H(curl).
    Nedelec,
    /// Raviart–Thomas H(div).
    Rt,
}

/// Build the dof functional list for one tet in the space's canonical
/// local order (edges by sorted local pair / faces by omitted vertex
/// with sorted triple / interior), all frames sorted-global.
fn element_dof_functionals(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    t: usize,
    r: usize,
    family: Family,
    chart: Chart,
) -> Vec<Dof> {
    let tet = complex.tets[t];
    let corners: [[f64; 3]; 4] = core::array::from_fn(|k| positions[tet[k] as usize]);
    let mut dofs = Vec::new();
    if family == Family::Nedelec {
        // Edge moments, sorted-global direction.
        for p in 0..4 {
            for q in (p + 1)..4 {
                let (ga, gb) = if tet[p] < tet[q] {
                    (tet[p], tet[q])
                } else {
                    (tet[q], tet[p])
                };
                for k in 0..r {
                    dofs.push(Dof::Edge {
                        a: positions[ga as usize],
                        b: positions[gb as usize],
                        k,
                    });
                }
            }
        }
    }
    // Face moments.
    let (nface_dirs, qdeg) = match family {
        Family::Nedelec => (2usize, r.saturating_sub(2)),
        Family::Rt => (1usize, r - 1),
    };
    let per_face = match family {
        Family::Nedelec => nedelec_entity_dofs(r).1,
        Family::Rt => rt_entity_dofs(r).0,
    };
    if per_face > 0 {
        for omit in 0..4 {
            let mut tri_g: Vec<u32> = (0..4).filter(|&l| l != omit).map(|l| tet[l]).collect();
            tri_g.sort_unstable();
            let tri: [[f64; 3]; 3] = core::array::from_fn(|k| positions[tri_g[k] as usize]);
            let e1 = [
                tri[1][0] - tri[0][0],
                tri[1][1] - tri[0][1],
                tri[1][2] - tri[0][2],
            ];
            let e2 = [
                tri[2][0] - tri[0][0],
                tri[2][1] - tri[0][1],
                tri[2][2] - tri[0][2],
            ];
            let dirs: Vec<[f64; 3]> = match family {
                Family::Nedelec => vec![e1, e2],
                Family::Rt => {
                    let n = cross3(e1, e2);
                    let nn = fs_math::det::sqrt(dot3(n, n));
                    vec![[n[0] / nn, n[1] / nn, n[2] / nn]]
                }
            };
            let mut count = 0usize;
            for qi in 0..=qdeg {
                for qj in 0..=(qdeg - qi) {
                    for dir in dirs.iter().take(nface_dirs) {
                        dofs.push(Dof::FaceMoment {
                            tri,
                            dir: *dir,
                            qi,
                            qj,
                        });
                        count += 1;
                    }
                }
            }
            debug_assert_eq!(count, per_face, "face dof count");
        }
    }
    // Interior moments.
    let per_cell = match family {
        Family::Nedelec => nedelec_entity_dofs(r).2,
        Family::Rt => rt_entity_dofs(r).1,
    };
    if per_cell > 0 {
        let qdeg = match family {
            Family::Nedelec => r - 3,
            Family::Rt => r - 2,
        };
        let mut count = 0usize;
        for qm in monomials(qdeg) {
            for dir in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
                dofs.push(Dof::CellMoment {
                    corners,
                    dir,
                    qm,
                    chart,
                });
                count += 1;
            }
        }
        debug_assert_eq!(count, per_cell, "interior dof count");
    }
    dofs
}

/// Build the spanning set of the trimmed space in LOCAL coordinates.
fn spanning_set(r: usize, family: Family) -> (Vec<VecPoly>, Vec<[usize; 3]>) {
    let monos = monomials(r);
    let n = monos.len();
    let mut span = Vec::new();
    // (P_{r−1})³.
    for (m, &[a, b, c]) in monos.iter().enumerate() {
        if a + b + c <= r - 1 {
            for k in 0..3 {
                let mut vp = VecPoly::zero(r);
                vp.comp[k][m] = 1.0;
                span.push(vp);
            }
        }
    }
    // Koszul part: homogeneous degree r−1 monomials ξ^α.
    for &[a, b, c] in &monos {
        if a + b + c != r - 1 {
            continue;
        }
        match family {
            Family::Rt => {
                // ξ·ξ^α: component k gets ξ^α·ξ_k.
                let mut vp = VecPoly::zero(r);
                for k in 0..3 {
                    let mut e = [a, b, c];
                    e[k] += 1;
                    let idx = monos.iter().position(|&mm| mm == e).expect("in table");
                    vp.comp[k][idx] = 1.0;
                }
                span.push(vp);
            }
            Family::Nedelec => {
                // ξ × (ξ^α eᵢ) over the INDEPENDENT subset: the Koszul
                // kernel is {W = ξh}, and {α₀ = 0 when i = 0} ∪ {all α
                // for i = 1, 2} is an exact complement (dimension
                // r + r(r+1) = r(r+2) = dim S_r) — the element
                // Vandermonde is SQUARE and solved by pivoted LU, not
                // normal equations (which cost √ε in the traces:
                // 1.6e-7 tangential jumps measured before this fix).
                for i in 0..3 {
                    if i == 0 && a > 0 {
                        continue;
                    }
                    let mut vp = VecPoly::zero(r);
                    // (ξ × w)_k = ε_{k,l,i} ξ_l w with w = ξ^α along eᵢ.
                    for (k, l, sign) in [
                        (0usize, 1usize, 1.0f64),
                        (0, 2, -1.0),
                        (1, 2, 1.0),
                        (1, 0, -1.0),
                        (2, 0, 1.0),
                        (2, 1, -1.0),
                    ] {
                        // ε_{klm}: (ξ × w)_k = Σ_l ε_{k l m} ξ_l w_m; here
                        // w_m ≠ 0 only for m = i.
                        let m_ = 3 - k - l; // the remaining index
                        if m_ != i {
                            continue;
                        }
                        let mut e = [a, b, c];
                        e[l] += 1;
                        let idx = monos.iter().position(|&mm| mm == e).expect("in table");
                        vp.comp[k][idx] += sign;
                    }
                    span.push(vp);
                }
            }
        }
    }
    let _ = n;
    (span, monos)
}

/// Build one element's dof-dual basis.
///
/// # Panics
/// On degenerate geometry or if r is outside 1..=[`MAX_R`].
#[must_use]
pub fn build_element(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    t: usize,
    r: usize,
    family: Family,
) -> ElementBasis {
    assert!(r >= 1 && r <= MAX_R, "vecfam supports r = 1..={MAX_R}");
    let tet = complex.tets[t];
    let corners: [[f64; 3]; 4] = core::array::from_fn(|k| positions[tet[k] as usize]);
    let chart = Chart::of(&corners);
    let dofs = element_dof_functionals(complex, positions, t, r, family, chart);
    let (span, monos) = spanning_set(r, family);
    let ndof = dofs.len();
    let nspan = span.len();
    let nq = r + 3;
    // A[dof][span] = functional(span function).
    let mut a = vec![0.0f64; ndof * nspan];
    for (j, sp) in span.iter().enumerate() {
        let field = |p: [f64; 3]| sp.eval_local(&monos, chart.local(p));
        for (i, d) in dofs.iter().enumerate() {
            a[i * nspan + j] = d.apply(&field, nq);
        }
    }
    // Basis coefficients in span coordinates.
    let coefs: Vec<Vec<f64>> = if nspan == ndof {
        let f = lu(&a, ndof).expect("vector-family Vandermonde nonsingular");
        (0..ndof)
            .map(|i| {
                let mut rhs = vec![0.0f64; ndof];
                rhs[i] = 1.0;
                f.solve(&mut rhs);
                rhs
            })
            .collect()
    } else {
        // Least-norm: c = Aᵀ(AAᵀ)⁻¹ e (A has full row rank by
        // unisolvence of the classical dofs).
        let mut aat = vec![0.0f64; ndof * ndof];
        for i in 0..ndof {
            for j in 0..ndof {
                let mut s = 0.0f64;
                for k in 0..nspan {
                    s = a[i * nspan + k].mul_add(a[j * nspan + k], s);
                }
                aat[i * ndof + j] = s;
            }
        }
        let f = cholesky(&aat, ndof).expect("dof Gram SPD (unisolvence)");
        (0..ndof)
            .map(|i| {
                let mut y = vec![0.0f64; ndof];
                y[i] = 1.0;
                f.solve(&mut y);
                let mut c = vec![0.0f64; nspan];
                for (row, &yr) in y.iter().enumerate() {
                    if yr != 0.0 {
                        for (k, ck) in c.iter_mut().enumerate() {
                            *ck = a[row * nspan + k].mul_add(yr, *ck);
                        }
                    }
                }
                c
            })
            .collect()
    };
    // Expand span coordinates into monomial VecPolys.
    let funcs: Vec<VecPoly> = coefs
        .iter()
        .map(|c| {
            let mut vp = VecPoly::zero(r);
            for (k, &ck) in c.iter().enumerate() {
                if ck != 0.0 {
                    for comp in 0..3 {
                        for (m, coef) in span[k].comp[comp].iter().enumerate() {
                            vp.comp[comp][m] = coef.mul_add(ck, vp.comp[comp][m]);
                        }
                    }
                }
            }
            vp
        })
        .collect();
    ElementBasis {
        chart,
        funcs,
        monos,
    }
}

// ------------------------------------------------------------------
// Global spaces.
// ------------------------------------------------------------------

/// A global H(curl) or H(div) space on a tet complex.
pub struct VecSpace<'c> {
    /// The complex.
    pub complex: &'c TetComplex,
    /// Order.
    pub r: usize,
    /// Family.
    pub family: Family,
    /// Per-edge dofs (0 for RT).
    pub per_edge: usize,
    /// Per-face dofs.
    pub per_face: usize,
    /// Per-cell dofs.
    pub per_cell: usize,
    /// Block offsets: edges at 0, then faces, then cells.
    pub face_off: usize,
    /// Cell block offset.
    pub cell_off: usize,
    /// Total dofs.
    pub ndof: usize,
    /// Built element bases (deterministic order).
    pub elements: Vec<ElementBasis>,
}

impl<'c> VecSpace<'c> {
    /// Build the global space (all element bases up front).
    ///
    /// # Panics
    /// On unsupported r or degenerate elements.
    #[must_use]
    pub fn new(
        complex: &'c TetComplex,
        positions: &[[f64; 3]],
        r: usize,
        family: Family,
    ) -> VecSpace<'c> {
        let (pe, pf, pc) = match family {
            Family::Nedelec => nedelec_entity_dofs(r),
            Family::Rt => {
                let (f, c) = rt_entity_dofs(r);
                (0, f, c)
            }
        };
        let face_off = complex.edges.len() * pe;
        let cell_off = face_off + complex.faces.len() * pf;
        let ndof = cell_off + complex.tets.len() * pc;
        let elements = (0..complex.tets.len())
            .map(|t| build_element(complex, positions, t, r, family))
            .collect();
        VecSpace {
            complex,
            r,
            family,
            per_edge: pe,
            per_face: pf,
            per_cell: pc,
            face_off,
            cell_off,
            ndof,
            elements,
        }
    }

    /// Global dof ids of element `t`, in the element's local dof order.
    #[must_use]
    pub fn element_dofs(&self, t: usize) -> Vec<usize> {
        let tet = self.complex.tets[t];
        let mut out = Vec::new();
        if self.per_edge > 0 {
            for p in 0..4 {
                for q in (p + 1)..4 {
                    let key = if tet[p] < tet[q] {
                        [tet[p], tet[q]]
                    } else {
                        [tet[q], tet[p]]
                    };
                    let e = self
                        .complex
                        .edges
                        .binary_search(&key)
                        .expect("edge in table");
                    for k in 0..self.per_edge {
                        out.push(e * self.per_edge + k);
                    }
                }
            }
        }
        if self.per_face > 0 {
            for omit in 0..4 {
                let mut tri = [0u32; 3];
                let mut c = 0;
                for (i, &v) in tet.iter().enumerate() {
                    if i != omit {
                        tri[c] = v;
                        c += 1;
                    }
                }
                tri.sort_unstable();
                let f = self
                    .complex
                    .faces
                    .binary_search(&tri)
                    .expect("face in table");
                for k in 0..self.per_face {
                    out.push(self.face_off + f * self.per_face + k);
                }
            }
        }
        for k in 0..self.per_cell {
            out.push(self.cell_off + t * self.per_cell + k);
        }
        out
    }

    /// Evaluate the global field with dof vector `u` inside element `t`.
    #[must_use]
    pub fn eval_in(&self, t: usize, u: &[f64], p: [f64; 3]) -> [f64; 3] {
        let el = &self.elements[t];
        let xi = el.chart.local(p);
        let mut v = [0.0f64; 3];
        for (l, &g) in self.element_dofs(t).iter().enumerate() {
            if u[g] != 0.0 {
                let b = el.funcs[l].eval_local(&el.monos, xi);
                for k in 0..3 {
                    v[k] = u[g].mul_add(b[k], v[k]);
                }
            }
        }
        v
    }

    /// Assemble the global mass matrix ∫ uᵢ·uⱼ dV.
    #[must_use]
    pub fn mass(&self, positions: &[[f64; 3]]) -> Csr {
        let quad = duffy_quadrature(self.r + 3);
        let mut coo = Coo::new(self.ndof, self.ndof);
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let el = &self.elements[t];
            let corners: [[f64; 3]; 4] = core::array::from_fn(|k| positions[tet[k] as usize]);
            let e1 = [
                corners[1][0] - corners[0][0],
                corners[1][1] - corners[0][1],
                corners[1][2] - corners[0][2],
            ];
            let e2 = [
                corners[2][0] - corners[0][0],
                corners[2][1] - corners[0][1],
                corners[2][2] - corners[0][2],
            ];
            let e3 = [
                corners[3][0] - corners[0][0],
                corners[3][1] - corners[0][1],
                corners[3][2] - corners[0][2],
            ];
            let vol6 = dot3(e1, cross3(e2, e3)).abs();
            let dofs = self.element_dofs(t);
            let nl = dofs.len();
            let mut me = vec![0.0f64; nl * nl];
            for &(lam, w) in &quad {
                let mut p = [0.0f64; 3];
                for (a, corner) in corners.iter().enumerate() {
                    for c in 0..3 {
                        p[c] = lam[a].mul_add(corner[c], p[c]);
                    }
                }
                let xi = el.chart.local(p);
                let vals: Vec<[f64; 3]> = el
                    .funcs
                    .iter()
                    .map(|f| f.eval_local(&el.monos, xi))
                    .collect();
                for i in 0..nl {
                    for j in 0..nl {
                        me[i * nl + j] =
                            (w * vol6).mul_add(dot3(vals[i], vals[j]), me[i * nl + j]);
                    }
                }
            }
            for (i, &gi) in dofs.iter().enumerate() {
                for (j, &gj) in dofs.iter().enumerate() {
                    coo.push(gi, gj, me[i * nl + j]);
                }
            }
        }
        coo.assemble()
    }

    /// Canonical interpolation of an analytic field: apply every
    /// GLOBAL dof functional (shared-entity functionals are identical
    /// from all incident elements — sorted-global frames — so each is
    /// computed once from its first incident element).
    #[must_use]
    pub fn interpolate<F: Fn([f64; 3]) -> [f64; 3]>(
        &self,
        positions: &[[f64; 3]],
        u: &F,
    ) -> Vec<f64> {
        let mut out = vec![0.0f64; self.ndof];
        let mut done = vec![false; self.ndof];
        let nq = self.r + 3;
        for t in 0..self.complex.tets.len() {
            let el = &self.elements[t];
            let fns = element_dof_functionals(
                self.complex,
                positions,
                t,
                self.r,
                self.family,
                el.chart,
            );
            for (l, &g) in self.element_dofs(t).iter().enumerate() {
                if !done[g] {
                    out[g] = fns[l].apply(u, nq);
                    done[g] = true;
                }
            }
        }
        out
    }

    /// L2 error of a dof vector against an analytic field.
    #[must_use]
    pub fn l2_error<F: Fn([f64; 3]) -> [f64; 3]>(
        &self,
        positions: &[[f64; 3]],
        u: &[f64],
        exact: &F,
    ) -> f64 {
        let quad = duffy_quadrature(self.r + 4);
        let mut total = 0.0f64;
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let corners: [[f64; 3]; 4] = core::array::from_fn(|k| positions[tet[k] as usize]);
            let e1 = [
                corners[1][0] - corners[0][0],
                corners[1][1] - corners[0][1],
                corners[1][2] - corners[0][2],
            ];
            let e2 = [
                corners[2][0] - corners[0][0],
                corners[2][1] - corners[0][1],
                corners[2][2] - corners[0][2],
            ];
            let e3 = [
                corners[3][0] - corners[0][0],
                corners[3][1] - corners[0][1],
                corners[3][2] - corners[0][2],
            ];
            let vol6 = dot3(e1, cross3(e2, e3)).abs();
            for &(lam, w) in &quad {
                let mut p = [0.0f64; 3];
                for (a, corner) in corners.iter().enumerate() {
                    for c in 0..3 {
                        p[c] = lam[a].mul_add(corner[c], p[c]);
                    }
                }
                let uh = self.eval_in(t, u, p);
                let ue = exact(p);
                let d = [uh[0] - ue[0], uh[1] - ue[1], uh[2] - ue[2]];
                total += w * vol6 * dot3(d, d);
            }
        }
        fs_math::det::sqrt(total)
    }
}

/// The discontinuous L² space P_{r−1} per cell (chart monomials).
pub struct DgSpace<'c> {
    /// The complex.
    pub complex: &'c TetComplex,
    /// Order r (polynomials of degree ≤ r−1).
    pub r: usize,
    /// Dofs per cell.
    pub per_cell: usize,
    /// Total dofs.
    pub ndof: usize,
}

impl<'c> DgSpace<'c> {
    /// Build.
    #[must_use]
    pub fn new(complex: &'c TetComplex, r: usize) -> DgSpace<'c> {
        let per_cell = dg_cell_dofs(r);
        DgSpace {
            complex,
            r,
            per_cell,
            ndof: complex.tets.len() * per_cell,
        }
    }
}

// ------------------------------------------------------------------
// Chain maps.
// ------------------------------------------------------------------

/// grad: H¹ (`SimplexSpace`, order r) → Nédélec order r, as the
/// interpolation matrix D[i][j] = dofᵢ(grad φⱼ) (grad P_r ⊂ N_r, so
/// the interpolant IS the gradient and curl∘grad = 0 to roundoff).
#[must_use]
pub fn grad_matrix(h1: &SimplexSpace<'_>, ned: &VecSpace<'_>, positions: &[[f64; 3]]) -> Csr {
    assert!(matches!(ned.family, Family::Nedelec), "grad lands in N_r");
    assert_eq!(h1.r, ned.r, "order mismatch");
    let geo = crate::whitney::element_geometry(ned.complex, positions);
    let mut coo = Coo::new(ned.ndof, h1.ndof);
    let mut done = vec![false; ned.ndof];
    let nq = ned.r + 3;
    for t in 0..ned.complex.tets.len() {
        let el = &ned.elements[t];
        let fns =
            element_dof_functionals(ned.complex, positions, t, ned.r, Family::Nedelec, el.chart);
        let h1_dofs = h1.element_dofs(t);
        let vdofs = ned.element_dofs(t);
        let corners: [[f64; 3]; 4] =
            core::array::from_fn(|k| positions[ned.complex.tets[t][k] as usize]);
        // Barycentric coordinates of a physical point (affine solve).
        let bary = make_bary(&corners);
        for (l, &g) in vdofs.iter().enumerate() {
            if done[g] {
                continue;
            }
            done[g] = true;
            for &(gj, lf) in &h1_dofs {
                let grad_phi = |p: [f64; 3]| -> [f64; 3] {
                    let lam = bary(p);
                    let dl = lf.d_lambda(lam, h1.r);
                    let mut gvec = [0.0f64; 3];
                    for (a, da) in dl.iter().enumerate() {
                        for c in 0..3 {
                            gvec[c] = da.mul_add(geo.grads[t][a][c], gvec[c]);
                        }
                    }
                    gvec
                };
                let v = fns[l].apply(&grad_phi, nq);
                if v.abs() > 1e-13 {
                    coo.push(g, gj, v);
                }
            }
        }
    }
    coo.assemble()
}

/// curl: Nédélec order r → RT order r (curl N_r ⊂ RT_r).
#[must_use]
pub fn curl_matrix(ned: &VecSpace<'_>, rt: &VecSpace<'_>, positions: &[[f64; 3]]) -> Csr {
    assert!(matches!(ned.family, Family::Nedelec), "domain is N_r");
    assert!(matches!(rt.family, Family::Rt), "codomain is RT_r");
    assert_eq!(ned.r, rt.r, "order mismatch");
    let mut coo = Coo::new(rt.ndof, ned.ndof);
    let mut done = vec![false; rt.ndof];
    let nq = rt.r + 3;
    for t in 0..rt.complex.tets.len() {
        let rt_el = &rt.elements[t];
        let fns =
            element_dof_functionals(rt.complex, positions, t, rt.r, Family::Rt, rt_el.chart);
        let rt_dofs = rt.element_dofs(t);
        let ned_dofs = ned.element_dofs(t);
        let ned_el = &ned.elements[t];
        // Exact curls of the Nédélec basis (physical = local / h).
        let curls: Vec<VecPoly> = ned_el
            .funcs
            .iter()
            .map(|f| f.curl_local(&ned_el.monos))
            .collect();
        for (l, &g) in rt_dofs.iter().enumerate() {
            if done[g] {
                continue;
            }
            done[g] = true;
            for (j, &gj) in ned_dofs.iter().enumerate() {
                let curl_j = |p: [f64; 3]| -> [f64; 3] {
                    let v = curls[j].eval_local(&ned_el.monos, ned_el.chart.local(p));
                    [
                        v[0] / ned_el.chart.h,
                        v[1] / ned_el.chart.h,
                        v[2] / ned_el.chart.h,
                    ]
                };
                let v = fns[l].apply(&curl_j, nq);
                if v.abs() > 1e-13 {
                    coo.push(g, gj, v);
                }
            }
        }
    }
    coo.assemble()
}

/// div: RT order r → DG P_{r−1} (L² moments in the cell chart).
#[must_use]
pub fn div_matrix(rt: &VecSpace<'_>, dg: &DgSpace<'_>, positions: &[[f64; 3]]) -> Csr {
    assert!(matches!(rt.family, Family::Rt), "domain is RT_r");
    assert_eq!(rt.r, dg.r, "order mismatch");
    let quad = duffy_quadrature(rt.r + 3);
    let dg_monos = monomials(rt.r - 1);
    let mut coo = Coo::new(dg.ndof, rt.ndof);
    for (t, tet) in rt.complex.tets.iter().enumerate() {
        let el = &rt.elements[t];
        let corners: [[f64; 3]; 4] = core::array::from_fn(|k| positions[tet[k] as usize]);
        let e1 = [
            corners[1][0] - corners[0][0],
            corners[1][1] - corners[0][1],
            corners[1][2] - corners[0][2],
        ];
        let e2 = [
            corners[2][0] - corners[0][0],
            corners[2][1] - corners[0][1],
            corners[2][2] - corners[0][2],
        ];
        let e3 = [
            corners[3][0] - corners[0][0],
            corners[3][1] - corners[0][1],
            corners[3][2] - corners[0][2],
        ];
        let vol6 = dot3(e1, cross3(e2, e3)).abs();
        let rt_dofs = rt.element_dofs(t);
        let divs: Vec<Vec<f64>> = el.funcs.iter().map(|f| f.div_local(&el.monos)).collect();
        // DG dof (t, m) = (1/vol)∫ q_m · div u — a moment; the DG "mass"
        // in these dofs is the monomial Gram, but for dd = 0 only the
        // LINEAR map matters.
        for (m, &[qa, qb, qc]) in dg_monos.iter().enumerate() {
            let row = t * dg.per_cell + m;
            for (j, &gj) in rt_dofs.iter().enumerate() {
                let mut acc = 0.0f64;
                let mut wsum = 0.0f64;
                for &(lam, w) in &quad {
                    let mut p = [0.0f64; 3];
                    for (a, corner) in corners.iter().enumerate() {
                        for c in 0..3 {
                            p[c] = lam[a].mul_add(corner[c], p[c]);
                        }
                    }
                    let xi = el.chart.local(p);
                    let q = xi[0].powi(i32::try_from(qa).expect("small"))
                        * xi[1].powi(i32::try_from(qb).expect("small"))
                        * xi[2].powi(i32::try_from(qc).expect("small"));
                    // div in physical coordinates = local / h.
                    let mut dv = 0.0f64;
                    for (mm, &[a2, b2, c2]) in el.monos.iter().enumerate() {
                        if divs[j][mm] != 0.0 {
                            let pm = xi[0].powi(i32::try_from(a2).expect("small"))
                                * xi[1].powi(i32::try_from(b2).expect("small"))
                                * xi[2].powi(i32::try_from(c2).expect("small"));
                            dv = divs[j][mm].mul_add(pm, dv);
                        }
                    }
                    dv /= el.chart.h;
                    acc = (w * vol6 * q).mul_add(dv, acc);
                    wsum += w * vol6;
                }
                let v = acc / wsum;
                if v.abs() > 1e-13 {
                    coo.push(row, gj, v);
                }
            }
        }
    }
    coo.assemble()
}

/// Barycentric-coordinate evaluator for a tet (affine inverse).
fn make_bary(corners: &[[f64; 3]; 4]) -> impl Fn([f64; 3]) -> [f64; 4] + '_ {
    // λ_a solves the 3×3 system with columns (p1−p0, p2−p0, p3−p0).
    let e = [
        [
            corners[1][0] - corners[0][0],
            corners[2][0] - corners[0][0],
            corners[3][0] - corners[0][0],
        ],
        [
            corners[1][1] - corners[0][1],
            corners[2][1] - corners[0][1],
            corners[3][1] - corners[0][1],
        ],
        [
            corners[1][2] - corners[0][2],
            corners[2][2] - corners[0][2],
            corners[3][2] - corners[0][2],
        ],
    ];
    let det = e[0][0].mul_add(
        e[1][1].mul_add(e[2][2], -(e[1][2] * e[2][1])),
        e[0][1].mul_add(
            -(e[1][0].mul_add(e[2][2], -(e[1][2] * e[2][0]))),
            e[0][2] * e[1][0].mul_add(e[2][1], -(e[1][1] * e[2][0])),
        ),
    );
    let inv = [
        [
            e[1][1].mul_add(e[2][2], -(e[1][2] * e[2][1])) / det,
            e[0][2].mul_add(e[2][1], -(e[0][1] * e[2][2])) / det,
            e[0][1].mul_add(e[1][2], -(e[0][2] * e[1][1])) / det,
        ],
        [
            e[1][2].mul_add(e[2][0], -(e[1][0] * e[2][2])) / det,
            e[0][0].mul_add(e[2][2], -(e[0][2] * e[2][0])) / det,
            e[0][2].mul_add(e[1][0], -(e[0][0] * e[1][2])) / det,
        ],
        [
            e[1][0].mul_add(e[2][1], -(e[1][1] * e[2][0])) / det,
            e[0][1].mul_add(e[2][0], -(e[0][0] * e[2][1])) / det,
            e[0][0].mul_add(e[1][1], -(e[0][1] * e[1][0])) / det,
        ],
    ];
    move |p: [f64; 3]| {
        let d = [
            p[0] - corners[0][0],
            p[1] - corners[0][1],
            p[2] - corners[0][2],
        ];
        let l1 = inv[0][0].mul_add(d[0], inv[0][1].mul_add(d[1], inv[0][2] * d[2]));
        let l2 = inv[1][0].mul_add(d[0], inv[1][1].mul_add(d[1], inv[1][2] * d[2]));
        let l3 = inv[2][0].mul_add(d[0], inv[2][1].mul_add(d[1], inv[2][2] * d[2]));
        [1.0 - l1 - l2 - l3, l1, l2, l3]
    }
}
