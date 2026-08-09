//! Loewner-framework identification: the data-driven matrix-pencil
//! front end used as an INDEPENDENT cross-check against vector fitting.
//!
//! Data on the imaginary axis is conjugate-augmented (so the underlying
//! system is real), split into interleaved left/right partitions, and
//! assembled into the Loewner pencil `(L, Ls)`. The standard real
//! transformation (per conjugate pair, rows/columns `[x; conj x]` map
//! to `[sqrt2 Re x; sqrt2 Im x]`) makes both matrices real without
//! changing the pencil's eigenvalues. Rank is revealed by SVD of the
//! stacked `[L; Ls]`; the projected pencil's generalized eigenvalues
//! are the identified poles. Residues then come from the SAME final
//! least-squares residue pass vector fitting uses, so the two front
//! ends share nothing upstream of the pole estimates.

use crate::model::RationalModel;
use crate::vf::{FitOptions, FitOutcome, VfError};
use fs_la::eigen_complex::{eig, lu_complex};
use fs_la::factor::svd_jacobi;
use fs_math::c64::C64;

/// Typed Loewner failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LoewnerError {
    /// Fewer than 2 samples per partition after augmentation.
    TooFewSamples {
        /// Samples provided.
        samples: usize,
    },
    /// The projected `E` pencil factor was singular (rank decision
    /// admitted too many states for the data).
    SingularPencil,
    /// Eigensolve failure inside the pencil reduction.
    Eig(fs_la::eigen_complex::EigFailure),
    /// The residue pass rejected the identified poles.
    Vf(VfError),
}

impl core::fmt::Display for LoewnerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LoewnerError::TooFewSamples { samples } => {
                write!(f, "loewner needs >= 4 samples, got {samples}")
            }
            LoewnerError::SingularPencil => write!(f, "projected loewner pencil is singular"),
            LoewnerError::Eig(e) => write!(f, "pencil eigensolve failed: {e}"),
            LoewnerError::Vf(e) => write!(f, "residue pass failed: {e}"),
        }
    }
}

impl std::error::Error for LoewnerError {}

/// Loewner identification at a caller-chosen order cap.
///
/// Returns the fitted model (poles from the pencil, residues from the
/// shared LS pass) plus the singular-value ratios used for the rank
/// decision (diagnostic: a clean system shows a hard drop at the true
/// order).
///
/// # Errors
/// Typed [`LoewnerError`].
pub fn loewner_fit(
    omega: &[f64],
    h: &[C64],
    max_order: usize,
    rank_tol: f64,
    opts: &FitOptions,
) -> Result<(FitOutcome, Vec<f64>), LoewnerError> {
    if omega.len() < 4 {
        return Err(LoewnerError::TooFewSamples {
            samples: omega.len(),
        });
    }
    // DIRECT-TERM STRIPPING, iterated: the Loewner pencil assumes
    // strictly proper data — a constant `d` shifts `Ls` by a rank-one
    // block and an improper `e*s` term corrupts every entry, biasing
    // the identified poles (executed: ~6e-4 relative pole bias,
    // Q-amplified to percent-level response error). So: crude d/e
    // estimate from the top of the band, strip, pencil, then let the
    // shared residue pass refit d/e at the identified poles and strip
    // again — converging ~8x per round on clean data.
    let mut d_est = if opts.fit_d { h[h.len() - 1].re } else { 0.0 };
    let mut e_est = if opts.fit_e && omega.len() >= 2 {
        let n = omega.len();
        ((h[n - 1].im - h[n - 2].im) / (omega[n - 1] - omega[n - 2])).max(0.0)
    } else {
        0.0
    };
    let mut best: Option<(FitOutcome, Vec<f64>)> = None;
    for _round in 0..6 {
        let (outcome, ratios) = pencil_round(omega, h, max_order, rank_tol, opts, d_est, e_est)?;
        d_est = outcome.model.d;
        e_est = outcome.model.e;
        let better = match &best {
            None => true,
            Some((b, _)) => outcome.report.weighted_rms < b.report.weighted_rms,
        };
        let converged = best.as_ref().is_some_and(|(b, _)| {
            (outcome.report.weighted_rms - b.report.weighted_rms).abs()
                <= 1.0e-3 * b.report.weighted_rms.max(f64::MIN_POSITIVE)
        });
        if better {
            best = Some((outcome, ratios));
        }
        if converged {
            break;
        }
    }
    Ok(best.expect("at least one round ran"))
}

/// One strip-then-pencil round: subsample the band for the pencil
/// (Loewner needs interpolation points, not the whole grid — the SVDs
/// are cubic), identify poles, then residue-fit on the FULL data.
fn pencil_round(
    omega: &[f64],
    h: &[C64],
    max_order: usize,
    rank_tol: f64,
    opts: &FitOptions,
    d_est: f64,
    e_est: f64,
) -> Result<(FitOutcome, Vec<f64>), LoewnerError> {
    // Subsample ~48 points evenly over the (assumed sorted) grid.
    let target = 48usize.min(omega.len());
    let stride = omega.len().div_ceil(target);
    let picks: Vec<usize> = (0..omega.len()).step_by(stride).collect();
    // Conjugate augmentation: (i*w, H) and (-i*w, conj H), kept
    // adjacent so the real transformation acts on 2x2 blocks.
    // Interleave sample pairs into left/right partitions.
    let mut left_pts: Vec<(C64, C64)> = Vec::new();
    let mut right_pts: Vec<(C64, C64)> = Vec::new();
    for (k, &i) in picks.iter().enumerate() {
        let w = omega[i];
        let hv = h[i] - C64::from_re(d_est) - C64::new(0.0, w * e_est);
        let bucket = if k % 2 == 0 {
            &mut left_pts
        } else {
            &mut right_pts
        };
        bucket.push((C64::new(0.0, w), hv));
        bucket.push((C64::new(0.0, -w), hv.conj()));
    }
    let (nl, nr) = (left_pts.len(), right_pts.len());
    // Complex Loewner and shifted Loewner.
    let mut lw = vec![C64::ZERO; nl * nr];
    let mut ls = vec![C64::ZERO; nl * nr];
    for (j, &(mu, v)) in left_pts.iter().enumerate() {
        for (i, &(lam, w)) in right_pts.iter().enumerate() {
            let denom = (mu - lam).recip();
            lw[j * nr + i] = (v - w) * denom;
            ls[j * nr + i] = (mu * v - lam * w) * denom;
        }
    }
    // Real transformation: per conjugate pair of LEFT rows,
    // [x; conj x] -> [sqrt2 Re x; sqrt2 Im x]; same on RIGHT columns.
    // After both, entries are real up to roundoff (asserted by taking
    // the real part; the imaginary residue is a diagnostic bound).
    let lw_r = realify(&lw, nl, nr);
    let ls_r = realify(&ls, nl, nr);
    // Rank decision on stacked [L; Ls] (2*nl x nr).
    let mut stacked = vec![0.0f64; 2 * nl * nr];
    stacked[..nl * nr].copy_from_slice(&lw_r);
    stacked[nl * nr..].copy_from_slice(&ls_r);
    let svd = svd_jacobi(&stacked, 2 * nl, nr);
    let s0 = svd.sigma[0].max(f64::MIN_POSITIVE);
    let ratios: Vec<f64> = svd.sigma.iter().map(|&s| s / s0).collect();
    let mut r = ratios.iter().filter(|&&t| t > rank_tol).count();
    r = r.min(max_order).min(nl).min(nr);
    if r == 0 {
        return Err(LoewnerError::SingularPencil);
    }
    // Left projector Y: nl x r from the leading left singular
    // directions of [L, Ls] side-by-side (nl x 2nr), obtained as right
    // singular vectors of its transpose.
    let mut side = vec![0.0f64; nl * 2 * nr];
    for j in 0..nl {
        side[j * 2 * nr..j * 2 * nr + nr].copy_from_slice(&lw_r[j * nr..(j + 1) * nr]);
        side[j * 2 * nr + nr..(j + 1) * 2 * nr].copy_from_slice(&ls_r[j * nr..(j + 1) * nr]);
    }
    let svd_side = svd_jacobi(&transpose(&side, nl, 2 * nr), 2 * nr, nl);
    // Columns of svd_side.v are right singular vectors of side^T, i.e.
    // LEFT singular vectors of side (nl-dimensional).
    let y = leading_cols(&svd_side.v, nl, r);
    let x = leading_cols(&svd.v, nr, r);
    // Projected pencil: E = -Y^T L X, A = -Y^T Ls X (signs cancel in
    // the generalized eigenproblem; poles = eig(E^{-1} A)).
    let e_p = project(&lw_r, nl, nr, &y, &x, r);
    let a_p = project(&ls_r, nl, nr, &y, &x, r);
    // Solve E^{-1} A via complex LU (E real, embedded).
    let e_c: Vec<C64> = e_p.iter().map(|&v| C64::from_re(v)).collect();
    let lu = lu_complex(&e_c, r).map_err(|_| LoewnerError::SingularPencil)?;
    let mut m = vec![C64::ZERO; r * r];
    for col in 0..r {
        let mut b: Vec<C64> = (0..r).map(|row| C64::from_re(a_p[row * r + col])).collect();
        lu.solve(&mut b);
        for row in 0..r {
            m[row * r + col] = b[row];
        }
    }
    let mut poles = eig(&m, r).map_err(LoewnerError::Eig)?;
    for p in &mut poles {
        if p.re > 0.0 {
            p.re = -p.re;
        }
        if p.re == 0.0 {
            p.re = -1.0e-6 * (1.0 + p.im.abs());
        }
    }
    let terms = crate::vf::terms_from_poles(&poles);
    // Residue pass on the FULL (unstripped) data: d and e are refit
    // there, which is what feeds the next stripping round.
    let outcome =
        crate::vf::residue_fit_at_poles(omega, h, &terms, opts).map_err(LoewnerError::Vf)?;
    Ok((outcome, ratios))
}

fn realify(m: &[C64], nl: usize, nr: usize) -> Vec<f64> {
    // Left rows come in adjacent conjugate pairs (2k, 2k+1); right
    // columns likewise. Apply the pair map on rows, then on columns,
    // and take real parts.
    let sqrt2 = fs_math::det::sqrt(2.0);
    let mut rows = vec![C64::ZERO; nl * nr];
    for k in 0..nl / 2 {
        for c in 0..nr {
            let a = m[(2 * k) * nr + c];
            let b = m[(2 * k + 1) * nr + c];
            // b = conj-partner row.
            rows[(2 * k) * nr + c] = (a + b).scale(1.0 / sqrt2);
            let diff = a - b;
            rows[(2 * k + 1) * nr + c] = C64::new(diff.im, -diff.re).scale(1.0 / sqrt2);
        }
    }
    let mut out = vec![0.0f64; nl * nr];
    for rrow in 0..nl {
        for k in 0..nr / 2 {
            let a = rows[rrow * nr + 2 * k];
            let b = rows[rrow * nr + 2 * k + 1];
            out[rrow * nr + 2 * k] = ((a + b).scale(1.0 / sqrt2)).re;
            let diff = a - b;
            out[rrow * nr + 2 * k + 1] = (C64::new(diff.im, -diff.re).scale(1.0 / sqrt2)).re;
        }
    }
    out
}

fn transpose(m: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; rows * cols];
    for i in 0..rows {
        for j in 0..cols {
            out[j * rows + i] = m[i * cols + j];
        }
    }
    out
}

/// First `r` columns of a row-major `n x n` matrix, as row-major `n x r`.
fn leading_cols(v: &[f64], n: usize, r: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * r];
    for i in 0..n {
        for j in 0..r {
            out[i * r + j] = v[i * n + j];
        }
    }
    out
}

/// `Y^T M X` with `Y: nl x r`, `M: nl x nr`, `X: nr x r` (row-major).
fn project(m: &[f64], nl: usize, nr: usize, y: &[f64], x: &[f64], r: usize) -> Vec<f64> {
    let mut mx = vec![0.0f64; nl * r];
    for i in 0..nl {
        for j in 0..r {
            let mut acc = 0.0;
            for k in 0..nr {
                acc += m[i * nr + k] * x[k * r + j];
            }
            mx[i * r + j] = acc;
        }
    }
    let mut out = vec![0.0f64; r * r];
    for i in 0..r {
        for j in 0..r {
            let mut acc = 0.0;
            for k in 0..nl {
                acc += y[k * r + i] * mx[k * r + j];
            }
            out[i * r + j] = acc;
        }
    }
    out
}

/// Agreement report between the two identification front ends.
#[derive(Debug, Clone)]
pub struct CrossCheck {
    /// Worst pole distance (each VF pole to its nearest Loewner pole),
    /// relative to the pole magnitude.
    pub worst_pole_mismatch: f64,
    /// Max relative response deviation |H_vf - H_loew| / |H_data| over
    /// the samples.
    pub worst_response_mismatch: f64,
}

/// Compare two fitted models against each other on the data grid.
#[must_use]
pub fn cross_check(omega: &[f64], h: &[C64], a: &RationalModel, b: &RationalModel) -> CrossCheck {
    let pa = a.poles_expanded();
    let pb = b.poles_expanded();
    let mut worst_pole = 0.0f64;
    for p in &pa {
        let mut best = f64::INFINITY;
        for q in &pb {
            best = best.min((*p - *q).abs());
        }
        worst_pole = worst_pole.max(best / p.abs().max(1.0));
    }
    let mut worst_resp = 0.0f64;
    for (&w, &hv) in omega.iter().zip(h) {
        let dev = (a.eval_iw(w) - b.eval_iw(w)).abs() / hv.abs().max(f64::MIN_POSITIVE);
        worst_resp = worst_resp.max(dev);
    }
    CrossCheck {
        worst_pole_mismatch: worst_pole,
        worst_response_mismatch: worst_resp,
    }
}

// Re-exported helper hooks implemented in vf.rs (shared residue pass).
pub use crate::vf::{residue_fit_at_poles, terms_from_poles};
