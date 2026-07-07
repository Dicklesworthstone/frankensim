//! Whitney forms (lowest order, P₁Λᵏ) on tet complexes: element
//! geometry batched through fs-la's small-dense kernels (Jacobian
//! determinants and inverses per element — the layout those kernels
//! exist for), closed-form element mass matrices from the constant
//! barycentric gradients, and de-Rham maps (field → cochain by exact
//! low-order quadrature).
//!
//! Orientation bookkeeping follows fs-rep-mesh's sorted conventions
//! everywhere: edges [u, v] u<v oriented u→v, faces [a, b, c] sorted
//! with the a→b→c circulation, cells signed by the stored-order
//! volume. The commutation tests (R∘d = d∘R) in the battery are what
//! pin every sign.

use fs_la::batched::{BatchMat, batch_det, batch_inv};
use fs_rep_mesh::TetComplex;
use fs_sparse::{Coo, Csr};

/// Per-element geometry: signed volumes, barycentric gradients, and
/// the gradient Gram matrix — computed through batched kernels.
pub struct ElementGeometry {
    /// Signed volume of each tet in STORED vertex order (det J / 6).
    pub vol_signed: Vec<f64>,
    /// Barycentric gradients ∇λ_a (constant per tet), local a = 0..4.
    pub grads: Vec<[[f64; 3]; 4]>,
    /// Gram matrix g[a][b] = ∇λ_a · ∇λ_b per tet.
    pub gram: Vec<[[f64; 4]; 4]>,
}

/// Build element geometry for every tet.
///
/// # Panics
/// If a tet is degenerate (exactly singular Jacobian) — a mesh bug,
/// not a runtime condition.
#[must_use]
pub fn element_geometry(complex: &TetComplex, positions: &[[f64; 3]]) -> ElementGeometry {
    let nt = complex.tets.len();
    // Jacobian columns are edge vectors from vertex 0.
    let jac = BatchMat::from_fn(3, nt, |m, i, j| {
        let t = complex.tets[m];
        positions[t[j + 1] as usize][i] - positions[t[0] as usize][i]
    });
    let dets = batch_det(&jac);
    let mut jinv = BatchMat::zeros(3, nt);
    let flags = batch_inv(&jac, &mut jinv);
    assert!(
        flags.is_empty(),
        "degenerate tets in element_geometry: {flags:?}"
    );
    let mut vol_signed = Vec::with_capacity(nt);
    let mut grads = Vec::with_capacity(nt);
    let mut gram = Vec::with_capacity(nt);
    for m in 0..nt {
        vol_signed.push(dets[m] / 6.0);
        // ∇λ_a for a = 1..3 is row a−1 of J⁻¹; ∇λ_0 = −Σ.
        let mut g = [[0.0f64; 3]; 4];
        for a in 1..4 {
            for c in 0..3 {
                g[a][c] = jinv.get(m, a - 1, c);
            }
        }
        for c in 0..3 {
            g[0][c] = -(g[1][c] + g[2][c] + g[3][c]);
        }
        let mut gr = [[0.0f64; 4]; 4];
        for a in 0..4 {
            for b in 0..4 {
                gr[a][b] = g[a][0].mul_add(g[b][0], g[a][1].mul_add(g[b][1], g[a][2] * g[b][2]));
            }
        }
        grads.push(g);
        gram.push(gr);
    }
    ElementGeometry { vol_signed, grads, gram }
}

/// Local index (0..4) of a global vertex within a tet.
fn local_of(tet: [u32; 4], v: u32) -> usize {
    tet.iter().position(|&x| x == v).expect("vertex in tet")
}

fn edge_id(complex: &TetComplex, u: u32, v: u32) -> usize {
    let key = if u < v { [u, v] } else { [v, u] };
    complex.edges.binary_search(&key).expect("edge in table")
}

fn face_id(complex: &TetComplex, mut f: [u32; 3]) -> usize {
    f.sort_unstable();
    complex.faces.binary_search(&f).expect("face in table")
}

const fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

/// Assemble the P₁Λᵏ Whitney mass matrix (∫ wᵢ · wⱼ dV summed over
/// elements) into a deterministic CSR. Degrees: 0 = vertex hats,
/// 1 = edge Whitney forms, 2 = face Whitney forms, 3 = per-cell
/// constant densities (diagonal 1/|V|).
///
/// # Panics
/// If `degree > 3` or geometry/complex sizes mismatch.
#[must_use]
pub fn mass_matrix(complex: &TetComplex, geo: &ElementGeometry, degree: u8) -> Csr {
    assert_eq!(geo.vol_signed.len(), complex.tets.len(), "geometry mismatch");
    let n = crate::cochain::cell_count(complex, degree);
    let mut coo = Coo::new(n, n);
    // ∫ λ_p λ_q dV = V/20·(1 + δ_pq) — the only scalar integral needed.
    for (m, &tet) in complex.tets.iter().enumerate() {
        let vol = geo.vol_signed[m].abs();
        let s = |p: usize, q: usize| -> f64 {
            if p == q { vol / 10.0 } else { vol / 20.0 }
        };
        match degree {
            0 => {
                for a in 0..4 {
                    for b in 0..4 {
                        coo.push(tet[a] as usize, tet[b] as usize, s(a, b));
                    }
                }
            }
            1 => {
                // Global-sorted edge (u, v): w = λ_u ∇λ_v − λ_v ∇λ_u.
                let mut locals = [(0usize, 0usize, 0usize); 6];
                let mut c = 0;
                for p in 0..4 {
                    for q in (p + 1)..4 {
                        let (gu, gv) = if tet[p] < tet[q] { (p, q) } else { (q, p) };
                        locals[c] = (edge_id(complex, tet[p], tet[q]), gu, gv);
                        c += 1;
                    }
                }
                let gr = &geo.gram[m];
                for &(e, a, b) in &locals {
                    for &(f, cc, d) in &locals {
                        let val = s(a, cc).mul_add(
                            gr[b][d],
                            s(b, d).mul_add(
                                gr[a][cc],
                                -s(a, d) * gr[b][cc] - s(b, cc) * gr[a][d],
                            ),
                        );
                        coo.push(e, f, val);
                    }
                }
            }
            2 => {
                // Face [a, b, c] sorted: w = 2(λ_a u_a + λ_b u_b + λ_c u_c),
                // u_a = ∇λ_b × ∇λ_c (cyclic in the sorted order).
                let g = &geo.grads[m];
                let mut faces = Vec::with_capacity(4);
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
                    let la = local_of(tet, tri[0]);
                    let lb = local_of(tet, tri[1]);
                    let lc = local_of(tet, tri[2]);
                    let u = [
                        cross(g[lb], g[lc]),
                        cross(g[lc], g[la]),
                        cross(g[la], g[lb]),
                    ];
                    faces.push((face_id(complex, tri), [la, lb, lc], u));
                }
                for (fi, li, ui) in &faces {
                    for (fj, lj, uj) in &faces {
                        let mut val = 0.0f64;
                        for p in 0..3 {
                            for q in 0..3 {
                                val = (4.0 * s(li[p], lj[q])).mul_add(dot(ui[p], uj[q]), val);
                            }
                        }
                        coo.push(*fi, *fj, val);
                    }
                }
            }
            3 => {
                coo.push(m, m, 1.0 / vol);
            }
            _ => panic!("mass_matrix degree must be 0..=3"),
        }
    }
    coo.assemble()
}

/// De-Rham map, degree 0: vertex point values.
#[must_use]
pub fn deram0<F: Fn([f64; 3]) -> f64>(positions: &[[f64; 3]], f: &F) -> Vec<f64> {
    positions.iter().map(|&p| f(p)).collect()
}

/// De-Rham map, degree 1: edge line integrals ∫ A·dl along u→v (u<v),
/// Simpson quadrature (exact for cubic integrands, hence exact for
/// quadratic vector fields).
#[must_use]
pub fn deram1<F: Fn([f64; 3]) -> [f64; 3]>(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    a: &F,
) -> Vec<f64> {
    complex
        .edges
        .iter()
        .map(|&[u, v]| {
            let (pu, pv) = (positions[u as usize], positions[v as usize]);
            let t = [pv[0] - pu[0], pv[1] - pu[1], pv[2] - pu[2]];
            let mid = [
                f64::midpoint(pu[0], pv[0]),
                f64::midpoint(pu[1], pv[1]),
                f64::midpoint(pu[2], pv[2]),
            ];
            let (a0, am, a1) = (a(pu), a(mid), a(pv));
            let simpson = [
                (4.0 * am[0] + a0[0] + a1[0]) / 6.0,
                (4.0 * am[1] + a0[1] + a1[1]) / 6.0,
                (4.0 * am[2] + a0[2] + a1[2]) / 6.0,
            ];
            dot(simpson, t)
        })
        .collect()
}

/// De-Rham map, degree 2: face fluxes ∫ B·n dA with the vector area
/// ((p_b−p_a)×(p_c−p_a))/2 of the SORTED triple (a→b→c circulation),
/// edge-midpoint quadrature (exact for quadratic fields).
#[must_use]
pub fn deram2<F: Fn([f64; 3]) -> [f64; 3]>(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    b: &F,
) -> Vec<f64> {
    complex
        .faces
        .iter()
        .map(|&[x, y, z]| {
            let (pa, pb, pc) = (
                positions[x as usize],
                positions[y as usize],
                positions[z as usize],
            );
            let e1 = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
            let e2 = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
            let ndA = cross(e1, e2).map(|c| 0.5 * c);
            let mids = [
                [
                    f64::midpoint(pa[0], pb[0]),
                    f64::midpoint(pa[1], pb[1]),
                    f64::midpoint(pa[2], pb[2]),
                ],
                [
                    f64::midpoint(pb[0], pc[0]),
                    f64::midpoint(pb[1], pc[1]),
                    f64::midpoint(pb[2], pc[2]),
                ],
                [
                    f64::midpoint(pc[0], pa[0]),
                    f64::midpoint(pc[1], pa[1]),
                    f64::midpoint(pc[2], pa[2]),
                ],
            ];
            let mut avg = [0.0f64; 3];
            for m in &mids {
                let bm = b(*m);
                for c in 0..3 {
                    avg[c] += bm[c] / 3.0;
                }
            }
            dot(avg, ndA)
        })
        .collect()
}

/// Parity of the permutation taking the stored tet order to sorted
/// order: +1 even, −1 odd (inversion count).
#[must_use]
pub fn sort_parity(tet: [u32; 4]) -> f64 {
    let mut inversions = 0usize;
    for i in 0..4 {
        for j in (i + 1)..4 {
            if tet[i] > tet[j] {
                inversions += 1;
            }
        }
    }
    if inversions % 2 == 0 { 1.0 } else { -1.0 }
}

/// De-Rham map, degree 3: signed cell integrals — centroid quadrature
/// (exact for affine integrands) times the SORTED-order orientation
/// sign (d2's convention is built on sorted tets, so the stored-order
/// signed volume is corrected by the stored→sorted parity).
#[must_use]
pub fn deram3<F: Fn([f64; 3]) -> f64>(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    f: &F,
) -> Vec<f64> {
    complex
        .tets
        .iter()
        .enumerate()
        .map(|(m, tet)| {
            let mut c = [0.0f64; 3];
            for &v in tet {
                for k in 0..3 {
                    c[k] += positions[v as usize][k] / 4.0;
                }
            }
            f(c) * geo.vol_signed[m] * sort_parity(*tet)
        })
        .collect()
}
