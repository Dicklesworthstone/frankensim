//! Impedance-form passivity: certification and convex residue repair.
//!
//! A rational impedance is passive iff it is positive-real; for a
//! STABLE model that reduces to `Re H(i*omega) >= 0` for all real
//! omega, plus non-negative asymptotic terms. Two detection layers:
//!
//! 1. Dense grid: log-spaced band sweep plus refinement around each
//!    pole's resonant frequency — fast, but a grid can miss a narrow
//!    dip.
//! 2. Hamiltonian eigenvalue test (exact): the purely imaginary
//!    eigenvalues of the positive-real Hamiltonian matrix are exactly
//!    the frequencies where `Re H` crosses zero. No grid to fool.
//!
//! DESCRIPTOR-FORM STATEMENT (improper models, `e != 0`): on the
//! imaginary axis `Re(i*omega*e) = 0`, so the `s*e` term is LOSSLESS
//! (an ideal inductance) and contributes NOTHING to `Re H`. The
//! Hamiltonian test therefore applies to the PROPER part `(A, B, C,
//! d)` unchanged, with the descriptor conditions reducing to `e >= 0`
//! (checked separately). This is the impedance-form specialization of
//! the descriptor-Hamiltonian test; the general improper case needs
//! the even-matrix-pencil machinery, recorded as a no-claim.
//!
//! The Hamiltonian arm needs `R = 2 d > 0`; when `d <= tol` the
//! certificate is honestly downgraded to `GridOnly` (named weaker
//! class) rather than silently claimed exact.
//!
//! Repair: minimize the weighted L2 residue perturbation subject to
//! `Re H(i*omega_v) >= margin` at the violation frequencies — a convex
//! QP solved by an active-set of KKT equality solves, iterated with
//! re-certification until the Hamiltonian test is green.

use crate::model::{PoleTerm, RationalModel};
use fs_la::eigen_complex::{EigFailure, eig};
use fs_la::factor::lu;
use fs_math::c64::C64;

/// How strong the passivity certificate is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertificateClass {
    /// Grid AND Hamiltonian crossing test both clean.
    HamiltonianExact,
    /// Grid clean but `d <= tol` makes the Hamiltonian arm unavailable
    /// — a named weaker certificate, not a silent claim.
    GridOnly,
}

/// Outcome of a passivity check.
#[derive(Debug, Clone)]
pub struct PassivityReport {
    /// True iff no violation was found by any available arm.
    pub passive: bool,
    /// Certificate strength when `passive`.
    pub class: CertificateClass,
    /// Worst (most negative) `Re H` seen on the grid and its frequency.
    pub worst: (f64, f64),
    /// Zero-crossing frequencies from the Hamiltonian test (rad/s,
    /// positive axis; empty when the arm is unavailable or clean).
    pub crossings: Vec<f64>,
    /// Frequencies at which violations were sampled (repair
    /// constraints are built here).
    pub violation_freqs: Vec<f64>,
}

/// Typed passivity-machinery failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PassivityError {
    /// Model has a pole on or right of the imaginary axis — passivity
    /// is undefined; fix stability first.
    Unstable,
    /// Hamiltonian eigensolve failed.
    Eig(EigFailure),
    /// The repair QP could not reach feasibility.
    RepairFailed {
        /// Outer certification rounds attempted.
        rounds: usize,
    },
    /// Asymptotic terms negative (`d < 0` or `e < 0`) — the LS solve
    /// should have prevented this; refuse rather than repair silently.
    NegativeAsymptote,
}

impl core::fmt::Display for PassivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PassivityError::Unstable => write!(f, "model is not stable"),
            PassivityError::Eig(e) => write!(f, "hamiltonian eigensolve failed: {e}"),
            PassivityError::RepairFailed { rounds } => {
                write!(f, "passivity repair failed after {rounds} rounds")
            }
            PassivityError::NegativeAsymptote => write!(f, "negative d or e asymptote"),
        }
    }
}

impl std::error::Error for PassivityError {}

const D_TOL: f64 = 1.0e-12;

/// Certify passivity of a stable impedance-form model.
///
/// # Errors
/// [`PassivityError::Unstable`] / [`PassivityError::NegativeAsymptote`]
/// on precondition failure, [`PassivityError::Eig`] if the Hamiltonian
/// eigensolve stalls.
pub fn check_passivity(
    model: &RationalModel,
    band: (f64, f64),
) -> Result<PassivityReport, PassivityError> {
    if !model.is_stable() {
        return Err(PassivityError::Unstable);
    }
    if model.d < 0.0 || model.e < 0.0 {
        return Err(PassivityError::NegativeAsymptote);
    }
    // Grid arm: log sweep over ~[band.0/10, band.1*10] plus per-pole
    // refinement (resonances are where Re H moves fastest).
    let grid = build_grid(model, band);
    let mut worst = (f64::INFINITY, 0.0f64);
    let mut violation_freqs: Vec<f64> = Vec::new();
    // Grid violations are compressed to ONE representative per
    // contiguous negative band (its argmin) plus the band edges: a
    // wide violation would otherwise hand the repair QP hundreds of
    // near-duplicate constraints and push the KKT solves cubic (the
    // executed hang that motivated this shape).
    let mut band_min: Option<(f64, f64)> = None; // (re, w) of current band
    let mut band_edges: (f64, f64) = (0.0, 0.0);
    for &w in &grid {
        let re = model.eval_iw(w).re;
        if re < worst.0 {
            worst = (re, w);
        }
        if re < 0.0 {
            match &mut band_min {
                None => {
                    band_min = Some((re, w));
                    band_edges = (w, w);
                }
                Some(m) => {
                    if re < m.0 {
                        *m = (re, w);
                    }
                    band_edges.1 = w;
                }
            }
        } else if let Some((_, wmin)) = band_min.take() {
            push_band(&mut violation_freqs, band_edges, wmin);
        }
    }
    if let Some((_, wmin)) = band_min.take() {
        push_band(&mut violation_freqs, band_edges, wmin);
    }
    // Hamiltonian arm (proper part; see descriptor-form statement).
    let mut crossings = Vec::new();
    let class = if model.d > D_TOL {
        crossings = hamiltonian_crossings(model).map_err(PassivityError::Eig)?;
        // Between/around crossings, sample midpoints to classify which
        // side is negative (crossings alone don't say).
        for pair in crossings.windows(2) {
            let mid = fs_math::det::sqrt(pair[0] * pair[1]);
            if model.eval_iw(mid).re < 0.0 && !violation_freqs.contains(&mid) {
                violation_freqs.push(mid);
            }
        }
        if let Some(&first) = crossings.first() {
            let below = first * 0.5;
            if model.eval_iw(below).re < 0.0 {
                violation_freqs.push(below);
            }
        }
        if let Some(&last) = crossings.last() {
            let above = last * 2.0;
            if model.eval_iw(above).re < 0.0 {
                violation_freqs.push(above);
            }
        }
        CertificateClass::HamiltonianExact
    } else {
        CertificateClass::GridOnly
    };
    violation_freqs.sort_by(f64::total_cmp);
    violation_freqs.dedup();
    let passive = violation_freqs.is_empty()
        && worst.0 >= 0.0
        && (class == CertificateClass::GridOnly || crossings.is_empty());
    Ok(PassivityReport {
        passive,
        class,
        worst,
        crossings,
        violation_freqs,
    })
}

/// Emit a violation band's representatives: both edges plus the
/// argmin when it is interior. The comparisons are IDENTITY checks on
/// values copied from the same grid, so strict equality is exact.
#[allow(clippy::float_cmp)]
fn push_band(out: &mut Vec<f64>, edges: (f64, f64), wmin: f64) {
    out.push(edges.0);
    if wmin != edges.0 && wmin != edges.1 {
        out.push(wmin);
    }
    out.push(edges.1);
}

fn build_grid(model: &RationalModel, band: (f64, f64)) -> Vec<f64> {
    let lo = (band.0 / 10.0).max(1.0e-3);
    let hi = band.1 * 10.0;
    let n = 2000usize;
    let mut grid = Vec::with_capacity(n + model.terms.len() * 32);
    let lr = fs_math::det::ln(hi / lo);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        grid.push(lo * fs_math::det::exp(t * lr));
    }
    for t in &model.terms {
        if let PoleTerm::Pair { pole, .. } = t {
            let w0 = pole.abs();
            for k in 0u32..32 {
                let f = 0.9 + 0.2 * (f64::from(k) / 31.0);
                grid.push(w0 * f);
            }
        }
    }
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    grid
}

/// Positive imaginary-axis crossings of `Re H` via the positive-real
/// Hamiltonian matrix of the proper part (requires `d > 0`):
/// `M = [[A - B C / (2d), -B B^T / (2d)], [C^T C / (2d), -(A - B C /
/// (2d))^T]]`; eigenvalues within the imaginary-axis tolerance report
/// their `|Im|` as crossing frequencies.
fn hamiltonian_crossings(model: &RationalModel) -> Result<Vec<f64>, EigFailure> {
    let ss = model.state_space();
    let n = ss.n;
    let r_inv = 1.0 / (2.0 * ss.d);
    let dim = 2 * n;
    let mut m = vec![C64::ZERO; dim * dim];
    // F = A - B C * r_inv
    let mut f = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            f[i * n + j] = ss.a[i * n + j] - ss.b[i] * ss.c[j] * r_inv;
        }
    }
    for i in 0..n {
        for j in 0..n {
            m[i * dim + j] = C64::from_re(f[i * n + j]);
            m[i * dim + n + j] = C64::from_re(-ss.b[i] * ss.b[j] * r_inv);
            m[(n + i) * dim + j] = C64::from_re(ss.c[i] * ss.c[j] * r_inv);
            m[(n + i) * dim + n + j] = C64::from_re(-f[j * n + i]);
        }
    }
    let eigs = eig(&m, dim)?;
    let mut crossings: Vec<f64> = Vec::new();
    for lam in eigs {
        let tol = 1.0e-7 * (1.0 + lam.abs());
        if lam.re.abs() <= tol && lam.im > 0.0 {
            crossings.push(lam.im);
        }
    }
    crossings.sort_by(f64::total_cmp);
    // Merge near-duplicates (conjugate symmetry already halved).
    crossings.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-6 * (*a + *b));
    Ok(crossings)
}

/// Repair outcome.
#[derive(Debug, Clone)]
pub struct RepairReport {
    /// Certification rounds (QP solve + re-check) used.
    pub rounds: usize,
    /// L2 norm of the residue perturbation, relative to the residue
    /// vector norm.
    pub relative_perturbation: f64,
    /// Worst KKT stationarity residual of the final QP solve
    /// (diagnostic; near machine-zero for a converged active set).
    pub kkt_residual: f64,
    /// The post-repair passivity report (green by construction on
    /// success).
    pub certificate: PassivityReport,
}

/// Convex residue-perturbation repair: returns the repaired model and
/// report. The poles, `d`, and `e` are held FIXED; only residues move
/// (the perturbation that cannot destabilize the model).
///
/// # Errors
/// [`PassivityError::RepairFailed`] when feasibility is not reached
/// within the round budget; precondition errors as in
/// [`check_passivity`].
pub fn repair_passivity(
    model: &RationalModel,
    band: (f64, f64),
) -> Result<(RationalModel, RepairReport), PassivityError> {
    const MAX_ROUNDS: usize = 12;
    let mut current = model.clone();
    let mut constraint_freqs: Vec<f64> = Vec::new();
    let mut last_kkt = 0.0f64;
    let base_norm = residue_norm(model);
    for round in 0..MAX_ROUNDS {
        let report = check_passivity(&current, band)?;
        if report.passive {
            let pert = perturbation_norm(model, &current);
            return Ok((
                current,
                RepairReport {
                    rounds: round,
                    relative_perturbation: pert / base_norm.max(f64::MIN_POSITIVE),
                    kkt_residual: last_kkt,
                    certificate: report,
                },
            ));
        }
        // Accumulate constraints: every violation frequency seen so
        // far, plus the grid-worst point, with a margin that grows a
        // little each round to overcome crossing regeneration.
        for &w in &report.violation_freqs {
            if !constraint_freqs
                .iter()
                .any(|&c| (c - w).abs() <= 1.0e-9 * w)
            {
                constraint_freqs.push(w);
            }
        }
        if report.worst.0 < 0.0 {
            let w = report.worst.1;
            if !constraint_freqs
                .iter()
                .any(|&c| (c - w).abs() <= 1.0e-9 * w)
            {
                constraint_freqs.push(w);
            }
        }
        let margin = report.worst.0.abs() * 0.05 * (round as f64 + 1.0);
        let (repaired, kkt) = qp_step(model, &constraint_freqs, margin);
        last_kkt = kkt;
        current = repaired;
    }
    Err(PassivityError::RepairFailed { rounds: MAX_ROUNDS })
}

fn residue_norm(model: &RationalModel) -> f64 {
    let mut sq = 0.0;
    for t in &model.terms {
        match t {
            PoleTerm::Real { residue, .. } => sq += residue * residue,
            PoleTerm::Pair { residue, .. } => sq += residue.norm_sq(),
        }
    }
    fs_math::det::sqrt(sq)
}

fn perturbation_norm(a: &RationalModel, b: &RationalModel) -> f64 {
    let mut sq = 0.0;
    for (ta, tb) in a.terms.iter().zip(&b.terms) {
        match (ta, tb) {
            (PoleTerm::Real { residue: ra, .. }, PoleTerm::Real { residue: rb, .. }) => {
                sq += (ra - rb) * (ra - rb);
            }
            (PoleTerm::Pair { residue: ra, .. }, PoleTerm::Pair { residue: rb, .. }) => {
                sq += (*ra - *rb).norm_sq();
            }
            _ => {}
        }
    }
    fs_math::det::sqrt(sq)
}

/// Real coordinates of the residue vector (same layout as vf).
fn residue_coords(model: &RationalModel) -> Vec<f64> {
    let mut v = Vec::new();
    for t in &model.terms {
        match *t {
            PoleTerm::Real { residue, .. } => v.push(residue),
            PoleTerm::Pair { residue, .. } => {
                v.push(residue.re);
                v.push(residue.im);
            }
        }
    }
    v
}

fn with_residue_coords(model: &RationalModel, coords: &[f64]) -> RationalModel {
    let mut terms = Vec::with_capacity(model.terms.len());
    let mut at = 0usize;
    for t in &model.terms {
        match *t {
            PoleTerm::Real { pole, .. } => {
                terms.push(PoleTerm::Real {
                    pole,
                    residue: coords[at],
                });
                at += 1;
            }
            PoleTerm::Pair { pole, .. } => {
                terms.push(PoleTerm::Pair {
                    pole,
                    residue: C64::new(coords[at], coords[at + 1]),
                });
                at += 2;
            }
        }
    }
    RationalModel {
        terms,
        d: model.d,
        e: model.e,
    }
}

/// Row of `d Re H(i*w) / d residue-coords` — linear because `H` is
/// linear in the residues.
fn re_gradient_row(model: &RationalModel, w: f64) -> Vec<f64> {
    let s = C64::new(0.0, w);
    let mut row = Vec::new();
    for t in &model.terms {
        match *t {
            PoleTerm::Real { pole, .. } => {
                row.push((s - C64::from_re(pole)).recip().re);
            }
            PoleTerm::Pair { pole, .. } => {
                let a = (s - pole).recip();
                let b = (s - pole.conj()).recip();
                // d/d(rho): Re(a + b); d/d(sigma): Re(i a - i b)
                row.push((a + b).re);
                let diff = a - b;
                row.push(-diff.im);
            }
        }
    }
    row
}

/// One QP solve: minimize ||delta||^2 (identity Hessian in residue
/// coordinates) subject to `Re H_base(w_v) + g(w_v) . delta >= margin`
/// for all accumulated constraint frequencies, treating the ORIGINAL
/// model as the anchor. Active-set over the inequality constraints:
/// all-violated start, KKT equality solves via fs-la LU, multipliers
/// pruned until sign-feasible.
fn qp_step(anchor: &RationalModel, freqs: &[f64], margin: f64) -> (RationalModel, f64) {
    let base = residue_coords(anchor);
    let nv = base.len();
    let rows: Vec<Vec<f64>> = freqs.iter().map(|&w| re_gradient_row(anchor, w)).collect();
    let vals: Vec<f64> = freqs.iter().map(|&w| anchor.eval_iw(w).re).collect();
    // Constraint: vals[j] + rows[j] . delta >= margin  =>
    //   rows[j] . delta >= margin - vals[j] = rhs_j
    let rhs: Vec<f64> = vals.iter().map(|&v| margin - v).collect();
    // Start with every constraint whose unperturbed slack is negative.
    let mut active: Vec<usize> = (0..freqs.len()).filter(|&j| rhs[j] > 0.0).collect();
    let mut delta = vec![0.0f64; nv];
    let mut kkt_res = 0.0f64;
    for _ in 0..50 {
        if active.is_empty() {
            break;
        }
        let na = active.len();
        // KKT: [ 2I  A^T ] [delta ]   [ 0 ]
        //      [ A   0   ] [lambda] = [rhs_active]
        let dim = nv + na;
        let mut k = vec![0.0f64; dim * dim];
        let mut b = vec![0.0f64; dim];
        for i in 0..nv {
            k[i * dim + i] = 2.0;
        }
        for (aj, &j) in active.iter().enumerate() {
            for i in 0..nv {
                k[i * dim + nv + aj] = rows[j][i];
                k[(nv + aj) * dim + i] = rows[j][i];
            }
            b[nv + aj] = rhs[j];
        }
        let Ok(fact) = lu(&k, dim) else {
            // Degenerate active set (duplicate constraints): drop the
            // newest and retry.
            active.pop();
            continue;
        };
        let mut sol = b.clone();
        fact.solve(&mut sol);
        delta = sol[..nv].to_vec();
        let lambda = &sol[nv..];
        // KKT stationarity residual: 2*delta + A^T lambda = 0.
        let mut worst = 0.0f64;
        for i in 0..nv {
            let mut acc = 2.0 * delta[i];
            for (aj, &j) in active.iter().enumerate() {
                acc += rows[j][i] * lambda[aj];
            }
            worst = worst.max(acc.abs());
        }
        kkt_res = worst;
        // Multiplier sign: for `>=` constraints in a min problem the
        // multipliers must be <= 0 in this sign convention exactly when
        // the constraint pushes AGAINST the objective; prune the most
        // wrong-signed one and re-solve.
        let mut drop_idx: Option<usize> = None;
        let mut worst_lam = 0.0f64;
        for (aj, &lam) in lambda.iter().enumerate() {
            if lam > worst_lam {
                worst_lam = lam;
                drop_idx = Some(aj);
            }
        }
        if let Some(aj) = drop_idx {
            active.remove(aj);
            continue;
        }
        // Add most-violated inactive constraint, if any.
        let mut add_idx: Option<usize> = None;
        let mut worst_gap = 1.0e-12;
        for j in 0..freqs.len() {
            if active.contains(&j) {
                continue;
            }
            let mut lhs = 0.0;
            for i in 0..nv {
                lhs += rows[j][i] * delta[i];
            }
            let gap = rhs[j] - lhs;
            if gap > worst_gap {
                worst_gap = gap;
                add_idx = Some(j);
            }
        }
        match add_idx {
            Some(j) => active.push(j),
            None => break,
        }
    }
    let coords: Vec<f64> = base.iter().zip(&delta).map(|(&b0, &d0)| b0 + d0).collect();
    (with_residue_coords(anchor, &coords), kkt_res)
}
