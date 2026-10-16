//! COHOMOLOGY (plan §8.1, bead tfz.7): harmonic cochains and the
//! discrete Hodge decomposition — the CORRECT treatment of
//! multiply-connected domains, where naive FEM silently produces the
//! wrong answer or an underdetermined system.
//!
//! In the diagonal-star inner products `⟨u, v⟩_k = uᵀ M_k v`:
//! - EXACT part: the M-orthogonal projection onto `im d_{k−1}`
//!   (normal equations `d ᵀ M d a = d ᵀ M x`, matrix-free CG);
//! - COEXACT part: the projection onto `im δ_{k+1}` via the SPD system
//!   `(d_k M_k⁻¹ d_kᵀ) w = d_k x`, coexact `= M_k⁻¹ d_kᵀ w`;
//! - HARMONIC part: the remainder — orthogonal to both by projection,
//!   with dimension equal to the Betti number `b_k` (cross-checked
//!   against the integer-rank Betti computation: geometry and physics
//!   agreeing is the internal consistency check the plan calls
//!   beautiful).
//!
//! The physical payoff is wired first-class: the CIRCULATION functional
//! (a cycle pairing) extracts `Γ` from the harmonic component of a flow
//! cochain — Kutta–Joukowski `L' = ρ V Γ` — and harmonic bases deflate
//! saddle systems so they are well-posed on domains with handles.

use crate::cochain::cell_count;
use crate::hodge::hodge_diagonal_barycentric;
use crate::whitney::ElementGeometry;
use fs_rep_mesh::{Incidence, TetComplex};

/// CG relative tolerance for the projection solves.
const CG_TOL: f64 = 1e-12;

/// Acceptance threshold for a harmonic candidate (relative M-norm).
const HARMONIC_FLOOR: f64 = 1e-8;

/// A sparse ±1 operator in CSR-like row form (from the exact integer
/// incidences) with its transpose.
struct Op {
    rows: Vec<Vec<(usize, f64)>>,
    cols: Vec<Vec<(usize, f64)>>,
    ncols: usize,
}

impl Op {
    fn of(inc: &Incidence, ncols: usize) -> Op {
        let rows: Vec<Vec<(usize, f64)>> = inc
            .rows
            .iter()
            .map(|r| r.iter().map(|&(c, s)| (c, f64::from(s))).collect())
            .collect();
        let mut cols = vec![Vec::new(); ncols];
        for (r, row) in rows.iter().enumerate() {
            for &(c, s) in row {
                cols[c].push((r, s));
            }
        }
        Op { rows, cols, ncols }
    }

    fn apply(&self, x: &[f64]) -> Vec<f64> {
        self.rows
            .iter()
            .map(|row| row.iter().map(|&(c, s)| s * x[c]).sum())
            .collect()
    }

    fn apply_t(&self, y: &[f64]) -> Vec<f64> {
        (0..self.ncols)
            .map(|c| self.cols[c].iter().map(|&(r, s)| s * y[r]).sum())
            .collect()
    }
}

/// Deterministic matrix-free CG (SPD / consistent-semidefinite).
fn cg(apply: &dyn Fn(&[f64]) -> Vec<f64>, b: &[f64], max_iters: usize) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let b_norm = rr.max(f64::MIN_POSITIVE);
    for _ in 0..max_iters {
        if rr <= CG_TOL * CG_TOL * b_norm {
            break;
        }
        let ap = apply(&p);
        let pap: f64 = p.iter().zip(&ap).map(|(a, c)| a * c).sum();
        if pap.abs() < f64::MIN_POSITIVE {
            break;
        }
        let alpha = rr / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
    }
    x
}

/// The operators + stars around degree k.
struct Setting {
    /// `d_{k−1}` (absent for k = 0).
    d_down: Option<Op>,
    /// `d_k` (absent for k = 3).
    d_up: Option<Op>,
    /// Diagonal star `M_k`.
    m: Vec<f64>,
    /// Diagonal star `M_{k−1}` (unused for k = 0).
    m_down: Vec<f64>,
}

fn setting(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    degree: u8,
) -> Setting {
    let counts = |d: u8| -> usize { cell_count(complex, d) };
    let inc = |d: u8| match d {
        0 => complex.d0(),
        1 => complex.d1(),
        2 => complex.d2(),
        _ => unreachable!("d only up to degree 2"),
    };
    let d_down = (degree > 0).then(|| Op::of(&inc(degree - 1), counts(degree - 1)));
    let d_up = (degree < 3).then(|| Op::of(&inc(degree), counts(degree)));
    let m = hodge_diagonal_barycentric(complex, positions, geo, degree);
    let m_down = if degree > 0 {
        hodge_diagonal_barycentric(complex, positions, geo, degree - 1)
    } else {
        Vec::new()
    };
    Setting {
        d_down,
        d_up,
        m,
        m_down,
    }
}

/// The Hodge decomposition of one degree-k cochain value vector.
#[derive(Debug, Clone, PartialEq)]
pub struct HodgeParts {
    /// The exact component `d a`.
    pub exact: Vec<f64>,
    /// The coexact component `δ b`.
    pub coexact: Vec<f64>,
    /// The harmonic remainder.
    pub harmonic: Vec<f64>,
    /// M-weighted orthogonality residuals (exact·coexact, exact·harmonic,
    /// coexact·harmonic), relative to ‖x‖²_M.
    pub ortho_residuals: [f64; 3],
}

fn m_dot(m: &[f64], a: &[f64], b: &[f64]) -> f64 {
    m.iter().zip(a).zip(b).map(|((w, x), y)| w * x * y).sum()
}

/// Hodge-decompose `x` at `degree` on the complex (diagonal stars).
///
/// # Panics
/// On a length mismatch between `x` and the k-cell count.
#[must_use]
pub fn hodge_decompose(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    degree: u8,
    x: &[f64],
) -> HodgeParts {
    let st = setting(complex, positions, geo, degree);
    assert_eq!(x.len(), st.m.len(), "cochain length");
    let n = x.len();
    let iters = 6 * n + 60;
    // Exact: solve dᵀ M d a = dᵀ M x.
    let exact = if let Some(d) = &st.d_down {
        let mx: Vec<f64> = st.m.iter().zip(x).map(|(w, v)| w * v).collect();
        let rhs = d.apply_t(&mx);
        let apply = |a: &[f64]| {
            let da = d.apply(a);
            let mda: Vec<f64> = st.m.iter().zip(&da).map(|(w, v)| w * v).collect();
            d.apply_t(&mda)
        };
        let a = cg(&apply, &rhs, iters);
        d.apply(&a)
    } else {
        vec![0.0; n]
    };
    let r1: Vec<f64> = x.iter().zip(&exact).map(|(a, b)| a - b).collect();
    // Coexact: solve (d M⁻¹ dᵀ) w = d r, coexact = M⁻¹ dᵀ w.
    let coexact = if let Some(d) = &st.d_up {
        let rhs = d.apply(&r1);
        let apply = |w: &[f64]| {
            let dtw = d.apply_t(w);
            let mi: Vec<f64> = st.m.iter().zip(&dtw).map(|(m, v)| v / m).collect();
            d.apply(&mi)
        };
        let w = cg(&apply, &rhs, iters);
        let dtw = d.apply_t(&w);
        st.m.iter().zip(&dtw).map(|(m, v)| v / m).collect()
    } else {
        vec![0.0; n]
    };
    let harmonic: Vec<f64> = r1.iter().zip(&coexact).map(|(a, b)| a - b).collect();
    let total = m_dot(&st.m, x, x).max(f64::MIN_POSITIVE);
    let ortho_residuals = [
        m_dot(&st.m, &exact, &coexact).abs() / total,
        m_dot(&st.m, &exact, &harmonic).abs() / total,
        m_dot(&st.m, &coexact, &harmonic).abs() / total,
    ];
    let _ = &st.m_down;
    HodgeParts {
        exact,
        coexact,
        harmonic,
        ortho_residuals,
    }
}

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

/// An M-orthonormal basis of the degree-k harmonic space, found by
/// Hodge-projecting deterministic pseudo-random seeds and Gram–Schmidt
/// filtering. `hint` bounds the search (Betti number + margin is
/// plenty). The returned dimension IS the computed `b_k` — cross-check
/// it against [`crate::betti_numbers`].
#[must_use]
pub fn harmonic_basis(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    degree: u8,
    hint: usize,
) -> Vec<Vec<f64>> {
    let m = hodge_diagonal_barycentric(complex, positions, geo, degree);
    let n = m.len();
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut state = 0x00c0_ffee_u64 + u64::from(degree);
    for _ in 0..(hint + 3) {
        let seed: Vec<f64> = (0..n).map(|_| lcg(&mut state)).collect();
        let parts = hodge_decompose(complex, positions, geo, degree, &seed);
        let mut h = parts.harmonic;
        // Modified Gram–Schmidt against the accepted basis.
        for b in &basis {
            let proj = m_dot(&m, &h, b);
            for (hi, bi) in h.iter_mut().zip(b) {
                *hi -= proj * bi;
            }
        }
        let norm = m_dot(&m, &h, &h).sqrt();
        let seed_norm = m_dot(&m, &seed, &seed).sqrt().max(f64::MIN_POSITIVE);
        if norm > HARMONIC_FLOOR * seed_norm {
            for hi in &mut h {
                *hi /= norm;
            }
            basis.push(h);
        }
    }
    basis
}

/// The circulation functional: the pairing of a 1-cochain with an edge
/// CYCLE (edge index, sign). Depends only on the cohomology class for
/// closed cochains — the wing payoff: `Γ` from the harmonic component,
/// lift `L' = ρ V Γ`.
#[must_use]
pub fn circulation(values: &[f64], cycle: &[(usize, f64)]) -> f64 {
    cycle.iter().map(|&(e, s)| s * values[e]).sum()
}

/// Deflate a cochain against a harmonic basis (M-orthogonal projection
/// removal) — the cohomology-aware solver setup: appending these
/// orthogonality constraints (rows `hᵀM`) makes saddle systems on
/// handled domains well-posed.
#[must_use]
pub fn deflate_harmonics(
    complex: &TetComplex,
    positions: &[[f64; 3]],
    geo: &ElementGeometry,
    degree: u8,
    basis: &[Vec<f64>],
    x: &[f64],
) -> Vec<f64> {
    let m = hodge_diagonal_barycentric(complex, positions, geo, degree);
    let mut out = x.to_vec();
    for b in basis {
        let proj = m_dot(&m, &out, b);
        for (oi, bi) in out.iter_mut().zip(b) {
            *oi -= proj * bi;
        }
    }
    out
}
