//! Tape-differentiated acquisition gradients (feature `tape-acq`,
//! a2g2 lane a): the fixed-bank q-EI surface is differentiable by
//! construction — μ(X), the posterior covariance, ITS CHOLESKY, and
//! the reparameterized samples are all scalar arithmetic, so fs-ad's
//! FrankenTorch scalar tape differentiates the whole chain with no
//! matrix-level Cholesky backward needed (one reverse pass per
//! gradient, O(cost) instead of O(q·d·cost) FD).
//!
//! DETERMINISM CLASS (inherited from the bridge, declared): tape
//! primals/gradients use ft's elementary functions, NOT fs-math's
//! det kernels — verified against FD to TOLERANCES, never bitwise,
//! and excluded from cross-ISA goldens. The production f64 q-EI
//! (`acq::q_expected_improvement`) stays on the det kernels.
//!
//! The q-EI surface is piecewise smooth (per-sample min over the
//! batch, hinge at zero improvement): at ties the tape takes the
//! branch the primal takes — a valid subgradient, and FD gates are
//! evaluated away from ties.

use crate::gp::Gp;
use fs_ad::Real;
use fs_ad::bridge::{TapeReal, reverse_gradient};

/// Matérn-5⁄2 over a generic `Real` (the taped twin of
/// `Kernel::eval`; the f64 path keeps fs-math's det kernels).
fn matern52<R: Real>(signal: f64, ell: &[f64], a: &[R], b: &[R]) -> R {
    let mut acc = R::zero();
    for ((ai, bi), l) in a.iter().zip(b).zip(ell) {
        let d = (*ai - *bi) * R::from_f64(1.0 / l);
        acc = d.mul_add(d, acc);
    }
    // r = 0 (the Σ diagonal, coincident candidates): √'s backward is
    // infinite at 0 and poisoned the whole gradient with NaN
    // (measured — and the FD gates passed VACUOUSLY on flat regions
    // before this was caught). Matérn-5⁄2 is C², so the true gradient
    // of k at r = 0 is exactly zero: the constant is CORRECT, not a
    // subgradient choice.
    if acc.value() < 1e-280 {
        return R::from_f64(signal);
    }
    let r = acc.sqrt();
    let s5 = R::from_f64(fs_math::det::sqrt(5.0));
    let a5 = s5 * r;
    let poly = R::one() + a5 + a5 * a5 * R::from_f64(1.0 / 3.0);
    poly * (-a5).exp() * R::from_f64(signal)
}

/// Fixed-bank q-EI over a generic `Real`: candidates are the
/// variables; the trained GP (data, α, training Cholesky) and the
/// z-bank are constants. Mirrors `q_expected_improvement` exactly.
#[allow(clippy::too_many_lines)]
fn qei_generic<R: Real>(
    gp: &Gp,
    train_x: &[Vec<f64>],
    l_train: &[f64],
    alpha: &[f64],
    cands: &[R],
    q: usize,
    dim: usize,
    f_best: f64,
    bank: &[f64],
) -> R {
    let n = train_x.len();
    let signal = gp.kernel.signal;
    let ell = &gp.kernel.lengthscales;
    // k*(x_b) against the training set; μ_b and V columns.
    let mut mu = Vec::with_capacity(q);
    let mut vcols: Vec<Vec<R>> = Vec::with_capacity(q);
    for b in 0..q {
        let xb = &cands[b * dim..(b + 1) * dim];
        let mut kstar: Vec<R> = Vec::with_capacity(n);
        for xi in train_x {
            let xi_r: Vec<R> = xi.iter().map(|&v| R::from_f64(v)).collect();
            kstar.push(matern52(signal, ell, &xi_r, xb));
        }
        let mut m = R::zero();
        for (ks, &al) in kstar.iter().zip(alpha) {
            m = ks.mul_add(R::from_f64(al), m);
        }
        mu.push(m);
        // v = L_train⁻¹ k* (forward substitution, constant L).
        let mut v = kstar;
        for i in 0..n {
            let mut acc = v[i];
            for j in 0..i {
                acc = R::from_f64(-l_train[i * n + j]).mul_add(v[j], acc);
            }
            v[i] = acc * R::from_f64(1.0 / l_train[i * n + i]);
        }
        vcols.push(v);
    }
    // Posterior covariance Σ (jittered) and its Cholesky in R.
    let mut sigma = vec![R::zero(); q * q];
    for a in 0..q {
        for b in 0..=a {
            let xa = &cands[a * dim..(a + 1) * dim];
            let xb = &cands[b * dim..(b + 1) * dim];
            let prior = matern52(signal, ell, xa, xb);
            let mut dot = R::zero();
            for (p, r) in vcols[a].iter().zip(&vcols[b]) {
                dot = p.mul_add(*r, dot);
            }
            let val = prior - dot;
            sigma[a * q + b] = val;
            sigma[b * q + a] = val;
        }
        sigma[a * q + a] = sigma[a * q + a] + R::from_f64(1e-8);
    }
    // Generic Cholesky (the tape differentiates straight through it).
    let mut lq = vec![R::zero(); q * q];
    for i in 0..q {
        for j in 0..=i {
            let mut acc = sigma[i * q + j];
            for k in 0..j {
                acc = acc - lq[i * q + k] * lq[j * q + k];
            }
            if i == j {
                lq[i * q + i] = acc.sqrt();
            } else {
                lq[i * q + j] = acc * lq[j * q + j].recip();
            }
        }
    }
    // Reparameterized samples over the fixed bank; per-sample min via
    // primal-branch selection (a valid subgradient at ties).
    let samples = bank.len() / q;
    // Anchor the root to the variables: when NO sample improves, a
    // bare zero constant has no grad node and backward fails
    // (RootDoesNotRequireGrad — measured); a 0·x term keeps the graph
    // rooted with exactly zero gradient contribution.
    let mut acc = cands[0] * R::from_f64(0.0);
    for s in 0..samples {
        let z = &bank[s * q..(s + 1) * q];
        let mut best: Option<R> = None;
        for i in 0..q {
            let mut f = mu[i];
            for j in 0..=i {
                f = lq[i * q + j].mul_add(R::from_f64(z[j]), f);
            }
            best = Some(match best {
                None => f,
                Some(cur) => {
                    if f.value() < cur.value() {
                        f
                    } else {
                        cur
                    }
                }
            });
        }
        let imp = R::from_f64(f_best) - best.expect("q >= 1");
        if imp.value() > 0.0 {
            acc = acc + imp;
        }
    }
    acc * R::from_f64(1.0 / samples as f64)
}

/// Exact q-EI gradient w.r.t. the flattened candidate block (q·d),
/// by ONE reverse pass through the taped chain (kernels → posterior
/// → Cholesky → reparameterized hinge). Returns (q-EI, gradient).
#[must_use]
pub fn qei_gradient(gp: &Gp, cands: &[Vec<f64>], f_best: f64, bank: &[f64]) -> (f64, Vec<f64>) {
    let q = cands.len();
    let dim = cands[0].len();
    let (train_x, l_train, alpha) = gp.training_view();
    let flat: Vec<f64> = cands.iter().flatten().copied().collect();
    reverse_gradient(&flat, |vars: &[TapeReal]| {
        qei_generic(gp, train_x, &l_train, alpha, vars, q, dim, f_best, bank)
    })
}

/// L-BFGS ASCENT on the fixed-bank q-EI surface — the tape-powered
/// argmax (replaces CMA-ES where dimension makes gradients
/// worthwhile): maximize by minimizing −q-EI with fs-ascent's L-BFGS,
/// the box handled by clamping inside the objective (polish-phase
/// optima are interior in practice; boundary-clamp gradient
/// inconsistency at the walls is accepted and documented). Returns
/// the best block, its q-EI, and the reverse-pass count.
#[must_use]
pub fn qei_ascent(
    gp: &Gp,
    start: &[Vec<f64>],
    f_best: f64,
    bank: &[f64],
    bounds: (f64, f64),
    iters: usize,
) -> (Vec<Vec<f64>>, f64, usize) {
    let (lo, hi) = bounds;
    let q = start.len();
    let dim = start[0].len();
    let flat0: Vec<f64> = start.iter().flatten().copied().collect();
    let mut evals = 0usize;
    let mut fg = |flat: &[f64]| -> (f64, Vec<f64>) {
        let cands: Vec<Vec<f64>> = (0..q)
            .map(|b| {
                flat[b * dim..(b + 1) * dim]
                    .iter()
                    .map(|v| v.clamp(lo, hi))
                    .collect()
            })
            .collect();
        evals += 1;
        let (val, grad) = qei_gradient(gp, &cands, f_best, bank);
        (-val, grad.iter().map(|g| -g).collect())
    };
    let mut st = fs_ascent::LbfgsState::new(&flat0, 8, &mut fg);
    st.run(&mut fg, &fs_ascent::StopRule::GradNorm(1e-10), iters);
    let x: Vec<Vec<f64>> = (0..q)
        .map(|b| {
            st.x[b * dim..(b + 1) * dim]
                .iter()
                .map(|v| v.clamp(lo, hi))
                .collect()
        })
        .collect();
    let fx = -st.f;
    (x, fx, evals)
}
