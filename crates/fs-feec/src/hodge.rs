//! Discrete Hodge stars. Two variants land here, with their tradeoffs
//! stated rather than implied:
//!
//! - GALERKIN star = the Whitney mass matrix M_k (SPD, consistent,
//!   couples neighbours — the accuracy-first choice; it is literally
//!   `whitney::mass_matrix`, re-exported under its Hodge name so call
//!   sites say what they mean).
//! - DIAGONAL BARYCENTRIC star: per-cell positive ratios built from
//!   uniform volume shares (each tet distributes its volume equally
//!   among its k-cells) over primal measures. Always positive on any
//!   mesh (the monotonicity-first choice), only low-order accurate.
//!
//! The circumcentric diagonal star is deliberately NOT here: on
//! non-well-centered meshes (including the Kuhn fixtures this crate
//! tests with) it produces negative dual measures, so landing it
//! without the well-centeredness machinery would be a certificate
//! without evidence. Recorded as follow-up scope in the CONTRACT.

use crate::whitney::ElementGeometry;
use fs_rep_mesh::TetComplex;
use fs_sparse::Csr;

/// Galerkin Hodge star: the P₁Λᵏ Whitney mass matrix.
#[must_use]
pub fn galerkin_star(complex: &TetComplex, geo: &ElementGeometry, degree: u8) -> Csr {
    crate::whitney::mass_matrix(complex, geo, degree)
}

fn norm(v: [f64; 3]) -> f64 {
    fs_math::det::sqrt(v[0].mul_add(v[0], v[1].mul_add(v[1], v[2] * v[2])))
}

/// Diagonal barycentric star: entry per k-cell = (uniform volume
/// share of incident tets) / (primal k-measure). Positive on every
/// valid mesh; k = 0 uses primal measure 1, k = 3 is 1/|V|.
///
/// # Panics
/// If `degree > 3`.
#[must_use]
pub fn hodge_diagonal_barycentric(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    degree: u8,
) -> Vec<f64> {
    let n = crate::cochain::cell_count(complex, degree);
    let mut dual = vec![0.0f64; n];
    match degree {
        0 => {
            for (m, tet) in complex.tets.iter().enumerate() {
                let share = geo.vol_signed[m].abs() / 4.0;
                for &v in tet {
                    dual[v as usize] += share;
                }
            }
            dual
        }
        1 => {
            for (m, tet) in complex.tets.iter().enumerate() {
                let share = geo.vol_signed[m].abs() / 6.0;
                for p in 0..4 {
                    for q in (p + 1)..4 {
                        let key = if tet[p] < tet[q] {
                            [tet[p], tet[q]]
                        } else {
                            [tet[q], tet[p]]
                        };
                        let e = complex.edges.binary_search(&key).expect("edge");
                        dual[e] += share;
                    }
                }
            }
            for (e, d) in dual.iter_mut().enumerate() {
                let [u, v] = complex.edges[e];
                let (pu, pv) = (positions[u as usize], positions[v as usize]);
                let len = norm([pv[0] - pu[0], pv[1] - pu[1], pv[2] - pu[2]]);
                *d /= len;
            }
            dual
        }
        2 => {
            for (m, tet) in complex.tets.iter().enumerate() {
                let share = geo.vol_signed[m].abs() / 4.0;
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
                    let f = complex.faces.binary_search(&tri).expect("face");
                    dual[f] += share;
                }
            }
            for (f, d) in dual.iter_mut().enumerate() {
                let [a, b, c] = complex.faces[f];
                let (pa, pb, pc) = (
                    positions[a as usize],
                    positions[b as usize],
                    positions[c as usize],
                );
                let e1 = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                let e2 = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
                let cx = [
                    e1[1].mul_add(e2[2], -(e1[2] * e2[1])),
                    e1[2].mul_add(e2[0], -(e1[0] * e2[2])),
                    e1[0].mul_add(e2[1], -(e1[1] * e2[0])),
                ];
                let area = 0.5 * norm(cx);
                *d /= area;
            }
            dual
        }
        3 => geo.vol_signed.iter().map(|v| 1.0 / v.abs()).collect(),
        _ => panic!("hodge degree must be 0..=3"),
    }
}
