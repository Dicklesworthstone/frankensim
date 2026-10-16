//! Mixed-precision solves with iterative refinement AS A POLICY
//! (plan §6.1): factor cheap, refine to target, escalate only where the
//! evidence demands it. The ladder is
//!
//!   f32-factor / f64-refine  →  f64-direct  →  f64-factor / dd-refine
//!
//! and every solve returns a [`RefineReport`] recording the ladder chosen,
//! why, the residual trajectory, and what was achieved — a precision
//! decision is EVIDENCE, never a silent downgrade.
//!
//! The dd rung computes residuals in double-double (`fs_math::dd`, the
//! relocated single implementation) — extended-precision residuals push
//! forward accuracy toward eps even at condition numbers where plain f64
//! refinement stalls.
//!
//! Determinism: fixed thresholds, fixed iteration order, deterministic
//! stall rules — same input, same ladder, same bits (tested).

use crate::factor::{FactorError, lu};
use fs_math::dd::Dd;

/// Which rung of the precision ladder produced the returned solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ladder {
    /// f32 factorization, f64 residuals + corrections.
    F32Refine,
    /// Plain f64 factorization and solve.
    F64Direct,
    /// f64 factorization, double-double residuals.
    F64DdRefine,
}

/// Accuracy request. `backward` bounds
/// ‖b − A·x‖∞ / (‖A‖∞·‖x‖∞ + ‖b‖∞); `forward` (optional) bounds the
/// RELATIVE FORWARD error ‖x − x*‖∞/‖x*‖∞. The distinction matters:
/// any stable f64 solve is backward-accurate to ~eps, but forward error
/// scales with κ — beating κ·eps forward REQUIRES extended-precision
/// residuals (the dd rung), which is exactly when it engages.
#[derive(Debug, Clone, Copy)]
pub struct ResidualTarget {
    /// Backward-error bound to hit (e.g. 1e-14 for f64-grade answers).
    pub backward: f64,
    /// Optional forward-error bound; None = backward-only.
    pub forward: Option<f64>,
}

/// The evidence object: what was decided, what happened, what was achieved.
#[derive(Debug, Clone)]
pub struct RefineReport {
    /// The rung that produced the answer.
    pub ladder: Ladder,
    /// Refinement steps taken on the final rung (0 = direct solve only).
    pub steps: usize,
    /// Achieved backward error.
    pub achieved: f64,
    /// Whether the target was met.
    pub converged: bool,
    /// True if a cheaper rung was tried and abandoned (stall/failure).
    pub escalated: bool,
    /// 1-norm condition estimate used for the policy decision (from the
    /// f32 factorization when available; f64 otherwise).
    pub condition_estimate: f64,
    /// Forward-error estimate for the returned x: on the dd rung this is
    /// the last correction ratio ‖d‖∞/‖x‖∞ (the Demmel-style estimate);
    /// on cheaper rungs it is the κ·backward bound. None when no forward
    /// target was requested.
    pub forward_estimate: Option<f64>,
    /// Backward error after each refinement step (the ledger trajectory).
    pub trajectory: Vec<f64>,
}

/// f32 unit roundoff.
const EPS32: f64 = 5.960_464_477_539_063e-8; // 2⁻²⁴
/// f64 unit roundoff.
const EPS64: f64 = 1.110_223_024_625_156_5e-16; // 2⁻⁵³
/// Refinement contraction is roughly κ·eps per step; demand headroom of
/// 1/16 before trusting a rung (deterministic policy constant).
const HEADROOM: f64 = 16.0;
/// Max refinement steps before declaring a stall.
const MAX_STEPS: usize = 30;

// ---------------------------------------------------------------------------
// Compact f32 LU (private: the cheap rung's factorization)
// ---------------------------------------------------------------------------

struct LuF32 {
    n: usize,
    data: Vec<f32>,
    perm: Vec<usize>,
}

fn lu_f32(a: &[f64], n: usize) -> Result<LuF32, FactorError> {
    let mut m: Vec<f32> = a.iter().map(|&v| v as f32).collect();
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        // Partial pivot, lowest index on ties (strict > scan).
        let (mut piv_row, mut piv_val) = (k, m[k * n + k].abs());
        for r in k + 1..n {
            let v = m[r * n + k].abs();
            if v > piv_val {
                piv_val = v;
                piv_row = r;
            }
        }
        if piv_val == 0.0 {
            return Err(FactorError::Singular { index: k });
        }
        if piv_row != k {
            for c in 0..n {
                m.swap(k * n + c, piv_row * n + c);
            }
            perm.swap(k, piv_row);
        }
        let piv = m[k * n + k];
        for r in k + 1..n {
            let mult = m[r * n + k] / piv;
            m[r * n + k] = mult;
            for c in k + 1..n {
                m[r * n + c] = (-mult).mul_add(m[k * n + c], m[r * n + c]);
            }
        }
    }
    Ok(LuF32 { n, data: m, perm })
}

impl LuF32 {
    /// Solve in f32 (input/output f64; the cast IS the precision policy).
    fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut y: Vec<f32> = self.perm.iter().map(|&p| b[p] as f32).collect();
        for i in 0..n {
            let mut v = y[i];
            for (k, &yk) in y.iter().enumerate().take(i) {
                v = (-self.data[i * n + k]).mul_add(yk, v);
            }
            y[i] = v;
        }
        for i in (0..n).rev() {
            let mut v = y[i];
            for (k, &yk) in y.iter().enumerate().take(n).skip(i + 1) {
                v = (-self.data[i * n + k]).mul_add(yk, v);
            }
            y[i] = v / self.data[i * n + i];
        }
        y.iter().map(|&v| f64::from(v)).collect()
    }

    /// 1-norm condition estimate from the f32 factors (cheap policy input;
    /// widened to f64 arithmetic outside the solves).
    fn condition_1(&self, a: &[f64]) -> f64 {
        let n = self.n;
        let norm_a = (0..n)
            .map(|j| (0..n).map(|i| a[i * n + j].abs()).sum::<f64>())
            .fold(0.0f64, f64::max);
        let mut x = vec![1.0 / n as f64; n];
        let mut est = 0.0f64;
        for _ in 0..4 {
            x = self.solve(&x);
            let new_est: f64 = x.iter().map(|v| v.abs()).sum();
            if new_est <= est {
                break;
            }
            est = new_est;
            // Steepest 1-norm direction via the sign vector (Hager).
            let z: Vec<f64> = x
                .iter()
                .map(|&v| if v >= 0.0 { 1.0 } else { -1.0 })
                .collect();
            let (mut best, mut best_i) = (0.0f64, 0usize);
            for (i, &v) in self.solve(&z).iter().enumerate() {
                if v.abs() > best {
                    best = v.abs();
                    best_i = i;
                }
            }
            x = vec![0.0; n];
            x[best_i] = 1.0;
        }
        est * norm_a
    }
}

// ---------------------------------------------------------------------------
// Residuals
// ---------------------------------------------------------------------------

/// r = b − A·x in f64 (fused, fixed order); returns (r, ‖r‖∞).
fn residual_f64(a: &[f64], n: usize, x: &[f64], b: &[f64]) -> (Vec<f64>, f64) {
    let mut r = vec![0.0f64; n];
    let mut norm = 0.0f64;
    for i in 0..n {
        let mut acc = b[i];
        for (j, &xj) in x.iter().enumerate() {
            acc = (-a[i * n + j]).mul_add(xj, acc);
        }
        r[i] = acc;
        norm = norm.max(acc.abs());
    }
    (r, norm)
}

/// r = b − A·x with the accumulation carried in DOUBLE-DOUBLE — the
/// extended-precision residual that lets refinement push forward error
/// toward eps even for badly conditioned systems.
fn residual_dd(a: &[f64], n: usize, x: &[f64], b: &[f64]) -> (Vec<f64>, f64) {
    let mut r = vec![0.0f64; n];
    let mut norm = 0.0f64;
    for i in 0..n {
        let mut acc = Dd::from_f64(b[i]);
        for (j, &xj) in x.iter().enumerate() {
            acc = acc - Dd::from_f64(a[i * n + j]) * Dd::from_f64(xj);
        }
        r[i] = acc.to_f64();
        norm = norm.max(r[i].abs());
    }
    (r, norm)
}

fn inf_norm_mat(a: &[f64], n: usize) -> f64 {
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j].abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
}

fn inf_norm_vec(v: &[f64]) -> f64 {
    v.iter().fold(0.0f64, |m, &x| m.max(x.abs()))
}

/// Backward error ‖r‖∞ / (‖A‖∞·‖x‖∞ + ‖b‖∞).
fn backward_error(r_norm: f64, a_norm: f64, x: &[f64], b_norm: f64) -> f64 {
    r_norm
        / a_norm
            .mul_add(inf_norm_vec(x), b_norm)
            .max(f64::MIN_POSITIVE)
}

// ---------------------------------------------------------------------------
// The adaptive driver
// ---------------------------------------------------------------------------

/// Solve A·x = b to the requested backward error, choosing the cheapest
/// precision ladder the CONDITION EVIDENCE supports and escalating on
/// stall. Returns the solution and the decision/evidence report.
///
/// # Errors
/// [`FactorError`] if the matrix is singular at every attempted precision.
pub fn solve_adaptive(
    a: &[f64],
    n: usize,
    b: &[f64],
    target: ResidualTarget,
) -> Result<(Vec<f64>, RefineReport), FactorError> {
    assert_eq!(a.len(), n * n, "a must be n*n = {}", n * n);
    assert_eq!(b.len(), n, "b must have length {n}");
    let a_norm = inf_norm_mat(a, n);
    let b_norm = inf_norm_vec(b);
    let mut escalated = false;

    // A rung refining with WORKING-precision residuals tops out at forward
    // error ~ κ·eps64; it can only serve a forward target with headroom.
    let forward_ok_without_dd =
        |cond: f64| target.forward.is_none_or(|f| cond * EPS64 * HEADROOM < f);

    // Rung 1: f32 factor + f64 refinement — if the f32 factorization holds
    // together and the condition estimate leaves headroom.
    if let Ok(f32_fact) = lu_f32(a, n) {
        let cond = f32_fact.condition_1(a);
        if cond * EPS32 * HEADROOM < 1.0 && forward_ok_without_dd(cond) {
            let mut x = f32_fact.solve(b);
            let mut trajectory = Vec::new();
            let mut prev = f64::INFINITY;
            let mut stalls = 0;
            for step in 0..MAX_STEPS {
                let (r, r_norm) = residual_f64(a, n, &x, b);
                let be = backward_error(r_norm, a_norm, &x, b_norm);
                trajectory.push(be);
                if be <= target.backward {
                    let forward_estimate = target.forward.map(|_| cond * be);
                    return Ok((
                        x,
                        RefineReport {
                            ladder: Ladder::F32Refine,
                            steps: step,
                            achieved: be,
                            converged: true,
                            escalated,
                            condition_estimate: cond,
                            trajectory,
                            forward_estimate,
                        },
                    ));
                }
                // Stall rule: two consecutive steps without halving.
                if be > 0.5 * prev {
                    stalls += 1;
                    if stalls >= 2 {
                        break;
                    }
                } else {
                    stalls = 0;
                }
                prev = be;
                let d = f32_fact.solve(&r);
                for (xi, di) in x.iter_mut().zip(&d) {
                    *xi += di;
                }
            }
            escalated = true; // tried the cheap rung, it wasn't enough
        } else {
            escalated = true; // condition evidence vetoed the cheap rung
        }
    } else {
        escalated = true; // singular in f32 (but maybe fine in f64)
    }

    // Rung 2/3: f64 factorization; direct answer, then dd refinement if
    // the target demands more than plain f64 delivers.
    let f = lu(a, n)?;
    let cond64 = f.condition_1(a);
    let mut x = b.to_vec();
    f.solve(&mut x);
    let (_, r_norm) = residual_f64(a, n, &x, b);
    let be_direct = backward_error(r_norm, a_norm, &x, b_norm);
    if be_direct <= target.backward && forward_ok_without_dd(cond64) {
        let forward_estimate = target.forward.map(|_| cond64 * be_direct);
        return Ok((
            x,
            RefineReport {
                ladder: Ladder::F64Direct,
                steps: 0,
                achieved: be_direct,
                converged: true,
                escalated,
                condition_estimate: cond64,
                trajectory: vec![be_direct],
                forward_estimate,
            },
        ));
    }
    // dd-residual refinement (top rung) — extracted driver.
    Ok(dd_refine(
        a,
        n,
        b,
        &f,
        target,
        DdCtx {
            a_norm,
            b_norm,
            cond64,
            be_direct,
            x,
        },
    ))
}

/// Context handed from the direct rung to the dd-refinement driver.
struct DdCtx {
    a_norm: f64,
    b_norm: f64,
    cond64: f64,
    be_direct: f64,
    x: Vec<f64>,
}

/// The top rung: dd-residual refinement with the correction-ratio forward
/// estimate (Demmel-style: when the correction stops mattering, x is as
/// good as f64 storage allows).
fn dd_refine(
    a: &[f64],
    n: usize,
    b: &[f64],
    f: &crate::factor::Lu,
    target: ResidualTarget,
    ctx: DdCtx,
) -> (Vec<f64>, RefineReport) {
    let DdCtx {
        a_norm,
        b_norm,
        cond64,
        be_direct,
        mut x,
    } = ctx;
    let mut trajectory = vec![be_direct];
    let mut prev = be_direct;
    let mut stalls = 0;
    let mut best = be_direct;
    let mut fwd_est = f64::INFINITY;
    for step in 1..=MAX_STEPS {
        let (r, r_norm) = residual_dd(a, n, &x, b);
        let mut d = r.clone();
        f.solve(&mut d);
        for (xi, di) in x.iter_mut().zip(&d) {
            *xi += di;
        }
        fwd_est = inf_norm_vec(&d) / inf_norm_vec(&x).max(f64::MIN_POSITIVE);
        let be = backward_error(r_norm, a_norm, &x, b_norm);
        trajectory.push(be);
        best = best.min(be);
        let fwd_met = target.forward.is_none_or(|ft| fwd_est <= ft);
        if be <= target.backward && fwd_met {
            return (
                x,
                RefineReport {
                    ladder: Ladder::F64DdRefine,
                    steps: step,
                    achieved: be,
                    converged: true,
                    escalated: true,
                    condition_estimate: cond64,
                    trajectory,
                    forward_estimate: Some(fwd_est),
                },
            );
        }
        if be > 0.5 * prev {
            stalls += 1;
            if stalls >= 2 {
                break;
            }
        } else {
            stalls = 0;
        }
        prev = be;
    }
    // Ran dry: report honestly (converged = false, best achieved recorded).
    (
        x,
        RefineReport {
            ladder: Ladder::F64DdRefine,
            steps: trajectory.len() - 1,
            achieved: best,
            converged: false,
            escalated: true,
            condition_estimate: cond64,
            trajectory,
            forward_estimate: Some(fwd_est),
        },
    )
}
