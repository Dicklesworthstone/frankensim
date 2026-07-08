//! Symmetric eigensolvers (plan §6.1): matrix-free Lanczos and LOBPCG for
//! extremal eigenpairs, plus a dense cyclic-Jacobi symmetric eigensolver
//! (the small-problem workhorse and the Rayleigh–Ritz kernel).
//!
//! MATRIX-FREE: operators are closures `apply(x, y)` writing y = A·x —
//! fs-sparse SpMV, stencils, and composed operators all plug in without
//! coupling fs-la to any matrix format.
//!
//! RESUMABLE: both iterative solvers expose plain-data state objects;
//! `checkpoint = clone`, and resuming N+M steps is BITWISE identical to
//! running N+M straight (tested) — the P2 property extended to
//! eigeniteration.
//!
//! Honesty: every returned [`EigenPair`] carries its TRUE residual norm
//! ‖A·v − λ·v‖₂ (recomputed through the operator, not the internal
//! estimate) — callers certify against that, not against our optimism.
//!
//! v1 notes: Lanczos uses FULL reorthogonalization (the conservative
//! special case of selective reorth; the ω-recurrence optimization is a
//! recorded refinement). LOBPCG uses the X−X_prev conjugate direction and
//! identity preconditioning by default. Eigenvector adjoints are the
//! recorded follow-up; dλ/dp = vᵀ(∂A/∂p)v composes caller-side with
//! fs-ad duals.

use crate::factor::qr;

/// One converged (or best-effort) eigenpair with its certified residual.
#[derive(Debug, Clone)]
pub struct EigenPair {
    /// Ritz value.
    pub value: f64,
    /// Ritz vector (unit 2-norm).
    pub vector: Vec<f64>,
    /// TRUE residual ‖A·v − λ·v‖₂, recomputed through the operator.
    pub residual: f64,
}

// ---------------------------------------------------------------------------
// Dense cyclic Jacobi (symmetric)
// ---------------------------------------------------------------------------

/// Eigendecomposition of a dense SYMMETRIC matrix (row-major n×n) by
/// cyclic Jacobi: returns (eigenvalues ascending, eigenvectors as columns
/// of a row-major n×n matrix). Deterministic sweep order; ascending sort
/// with lowest-index tie-break (P2).
#[must_use]
pub fn jacobi_eigh(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    assert_eq!(a.len(), n * n, "a must be n*n = {}", n * n);
    let mut m = a.to_vec();
    let mut v = vec![0.0f64; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..60 {
        let mut rotated = false;
        for p in 0..n {
            for q in p + 1..n {
                let apq = m[p * n + q];
                let scale = (m[p * n + p].abs() + m[q * n + q].abs()).max(f64::MIN_POSITIVE);
                if apq.abs() <= 1e-16 * scale {
                    continue;
                }
                rotated = true;
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (1.0 + theta * theta).sqrt())
                } else {
                    -1.0 / (-theta + (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                // Update rows/cols p and q of the symmetric matrix.
                for k in 0..n {
                    let (mkp, mkq) = (m[k * n + p], m[k * n + q]);
                    m[k * n + p] = c.mul_add(mkp, -(s * mkq));
                    m[k * n + q] = s.mul_add(mkp, c * mkq);
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p * n + k], m[q * n + k]);
                    m[p * n + k] = c.mul_add(mpk, -(s * mqk));
                    m[q * n + k] = s.mul_add(mpk, c * mqk);
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c.mul_add(vkp, -(s * vkq));
                    v[k * n + q] = s.mul_add(vkp, c * vkq);
                }
            }
        }
        if !rotated {
            break;
        }
    }
    // Sort ascending, lowest original index first on ties.
    let mut order: Vec<usize> = (0..n).collect();
    let vals: Vec<f64> = (0..n).map(|i| m[i * n + i]).collect();
    order.sort_by(|&x, &y| vals[x].partial_cmp(&vals[y]).unwrap().then(x.cmp(&y)));
    let mut evals = vec![0.0f64; n];
    let mut evecs = vec![0.0f64; n * n];
    for (slot, &src) in order.iter().enumerate() {
        evals[slot] = vals[src];
        for i in 0..n {
            evecs[i * n + slot] = v[i * n + src];
        }
    }
    (evals, evecs)
}

// ---------------------------------------------------------------------------
// Lanczos (full reorthogonalization)
// ---------------------------------------------------------------------------

/// Resumable Lanczos state: plain data; `clone()` IS a checkpoint.
#[derive(Debug, Clone)]
pub struct LanczosState {
    n: usize,
    basis: Vec<Vec<f64>>,
    alphas: Vec<f64>,
    betas: Vec<f64>,
    /// Residual direction for the next expansion step.
    next: Vec<f64>,
}

impl LanczosState {
    /// Fresh state with the DETERMINISTIC start vector (fixed pattern,
    /// normalized — no RNG, replay-stable by construction).
    #[must_use]
    pub fn new(n: usize) -> LanczosState {
        // fs-math strict sin, NOT std sin: the start vector must be
        // cross-ISA bit-identical (platform libm here cost us a golden-hash
        // divergence on x86 — caught by the trj verification pipeline).
        let mut v0: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.5 * fs_math::det::sin((i as f64) * 0.7))
            .collect();
        let nrm = norm2(&v0);
        for x in &mut v0 {
            *x /= nrm;
        }
        LanczosState {
            n,
            basis: Vec::new(),
            alphas: Vec::new(),
            betas: Vec::new(),
            next: v0,
        }
    }

    /// Krylov dimension built so far.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.basis.len()
    }
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        acc = x.mul_add(*y, acc);
    }
    acc
}

/// Advance Lanczos by `steps` expansions (full reorthogonalization —
/// deterministic, conservative), then extract the `k_want` extremal Ritz
/// pairs (`largest` picks which end). Residuals in the returned pairs are
/// TRUE operator residuals.
pub fn lanczos_run<Op>(
    op: &Op,
    state: &mut LanczosState,
    steps: usize,
    k_want: usize,
    largest: bool,
) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
{
    let n = state.n;
    let mut w = vec![0.0f64; n];
    for _ in 0..steps {
        let v = state.next.clone();
        op(&v, &mut w);
        let alpha = dot(&v, &w);
        for (wi, &vi) in w.iter_mut().zip(&v) {
            *wi = (-alpha).mul_add(vi, *wi);
        }
        if let Some(prev) = state.basis.last() {
            let beta_prev = *state.betas.last().unwrap_or(&0.0);
            for (wi, &pi) in w.iter_mut().zip(prev) {
                *wi = (-beta_prev).mul_add(pi, *wi);
            }
        }
        state.basis.push(v);
        state.alphas.push(alpha);
        // FULL reorthogonalization (two passes — "twice is enough").
        for _ in 0..2 {
            for b in &state.basis {
                let c = dot(&w, b);
                for (wi, &bi) in w.iter_mut().zip(b) {
                    *wi = (-c).mul_add(bi, *wi);
                }
            }
        }
        let beta = norm2(&w);
        if beta < 1e-300 {
            // Invariant subspace found: restart direction is irrelevant;
            // stop expanding (dim() stops growing; callers see it).
            state.betas.push(0.0);
            state.next = vec![0.0; n];
            break;
        }
        state.betas.push(beta);
        state.next = w.iter().map(|&x| x / beta).collect();
    }
    ritz_pairs(op, state, k_want, largest)
}

/// Rayleigh–Ritz on the current Krylov basis.
fn ritz_pairs<Op>(op: &Op, state: &LanczosState, k_want: usize, largest: bool) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
{
    let m = state.basis.len();
    if m == 0 {
        return Vec::new();
    }
    // Dense tridiagonal T.
    let mut t = vec![0.0f64; m * m];
    for i in 0..m {
        t[i * m + i] = state.alphas[i];
        if i + 1 < m {
            t[i * m + i + 1] = state.betas[i];
            t[(i + 1) * m + i] = state.betas[i];
        }
    }
    let (evals, evecs) = jacobi_eigh(&t, m);
    let n = state.n;
    let take = k_want.min(m);
    let idx: Vec<usize> = if largest {
        (m - take..m).rev().collect()
    } else {
        (0..take).collect()
    };
    let mut out = Vec::with_capacity(take);
    let mut av = vec![0.0f64; n];
    for &j in &idx {
        let mut y = vec![0.0f64; n];
        for (i, b) in state.basis.iter().enumerate() {
            let c = evecs[i * m + j];
            for (yk, &bk) in y.iter_mut().zip(b) {
                *yk = c.mul_add(bk, *yk);
            }
        }
        let nrm = norm2(&y);
        for v in &mut y {
            *v /= nrm;
        }
        op(&y, &mut av);
        let lam = evals[j];
        let mut rsq = 0.0f64;
        for (a, &v) in av.iter().zip(&y) {
            let r = (-lam).mul_add(v, *a);
            rsq = r.mul_add(r, rsq);
        }
        out.push(EigenPair {
            value: lam,
            vector: y,
            residual: rsq.sqrt(),
        });
    }
    out
}

// ---------------------------------------------------------------------------
// LOBPCG
// ---------------------------------------------------------------------------

/// Resumable LOBPCG state (plain data; `clone()` IS a checkpoint).
#[derive(Debug, Clone)]
pub struct LobpcgState {
    n: usize,
    b: usize,
    /// Current block X, row-major n×b, orthonormal columns.
    x: Vec<f64>,
    /// Previous block (for the conjugate direction); empty on iteration 0.
    x_prev: Vec<f64>,
    /// Iterations completed.
    pub iters: usize,
}

impl LobpcgState {
    /// Deterministic orthonormal start block (fixed pattern, QR-cleaned).
    #[must_use]
    pub fn new(n: usize, b: usize) -> LobpcgState {
        // The search block Z = [X | W | P] has up to 3b columns, and its QR
        // requires n >= columns; enforce that here rather than panicking deep
        // inside `lobpcg_run`'s second iteration.
        assert!(
            b >= 1 && 3 * b <= n,
            "block size {b} needs 3b <= n for the [X|W|P] search space; got n={n}"
        );
        let mut x0 = vec![0.0f64; n * b];
        for j in 0..b {
            for i in 0..n {
                // Strict sin (see LanczosState::new — same libm hazard).
                x0[i * b + j] = 1.0
                    + fs_math::det::sin((i as f64) * 0.31 + (j as f64) * 1.7)
                    + if i % (j + 2) == 0 { 0.5 } else { 0.0 };
            }
        }
        LobpcgState {
            n,
            b,
            x: orthonormalize(&x0, n, b),
            x_prev: Vec::new(),
            iters: 0,
        }
    }
}

/// Orthonormalize the columns of a row-major n×k block via QR (returns
/// Q's first k columns explicitly).
fn orthonormalize(block: &[f64], n: usize, k: usize) -> Vec<f64> {
    let f = qr(block, n, k);
    let mut q = vec![0.0f64; n * k];
    for j in 0..k {
        let mut e = vec![0.0f64; n];
        if j < k {
            e[j] = 1.0;
        }
        f.apply_q(&mut e);
        for i in 0..n {
            q[i * k + j] = e[i];
        }
    }
    q
}

/// Run LOBPCG for `steps` iterations toward the `largest` or smallest
/// end of the spectrum, with an optional preconditioner `prec` (maps a
/// residual block column to a preconditioned direction; identity =
/// `|r, out| out.copy_from_slice(r)`). Returns the current Ritz pairs
/// with TRUE residuals.
pub fn lobpcg_run<Op, Pr>(
    op: &Op,
    state: &mut LobpcgState,
    steps: usize,
    largest: bool,
    prec: &Pr,
) -> Vec<EigenPair>
where
    Op: Fn(&[f64], &mut [f64]),
    Pr: Fn(&[f64], &mut [f64]),
{
    let (n, b) = (state.n, state.b);
    for _ in 0..steps {
        // Ritz values for the current X.
        let ax = apply_block(op, &state.x, n, b);
        let lam = block_rayleigh(&state.x, &ax, n, b);
        // Residual block W = AX − X·diag(λ), preconditioned column-wise.
        let mut w = vec![0.0f64; n * b];
        for j in 0..b {
            let mut col = vec![0.0f64; n];
            for i in 0..n {
                col[i] = (-lam[j]).mul_add(state.x[i * b + j], ax[i * b + j]);
            }
            let mut pcol = vec![0.0f64; n];
            prec(&col, &mut pcol);
            for i in 0..n {
                w[i * b + j] = pcol[i];
            }
        }
        // Search space Z = [X | W | P], QR-orthonormalized.
        let p_cols = if state.x_prev.is_empty() { 0 } else { b };
        let zc = 2 * b + p_cols;
        let mut z = vec![0.0f64; n * zc];
        for i in 0..n {
            for j in 0..b {
                z[i * zc + j] = state.x[i * b + j];
                z[i * zc + b + j] = w[i * b + j];
            }
            for j in 0..p_cols {
                // P = X − X_prev (the classic conjugate direction).
                z[i * zc + 2 * b + j] = state.x[i * b + j] - state.x_prev[i * b + j];
            }
        }
        let zq = orthonormalize(&z, n, zc);
        // Rayleigh–Ritz on the reduced space.
        let az = apply_block(op, &zq, n, zc);
        let mut s = vec![0.0f64; zc * zc];
        for p in 0..zc {
            for q in 0..zc {
                let mut acc = 0.0f64;
                for i in 0..n {
                    acc = zq[i * zc + p].mul_add(az[i * zc + q], acc);
                }
                s[p * zc + q] = acc;
            }
        }
        // Exact symmetrization (Jacobi assumes it; roundoff breaks it).
        for p in 0..zc {
            for q in p + 1..zc {
                let avg = f64::midpoint(s[p * zc + q], s[q * zc + p]);
                s[p * zc + q] = avg;
                s[q * zc + p] = avg;
            }
        }
        let (evals, evecs) = jacobi_eigh(&s, zc);
        // Take b extremal Ritz vectors.
        let cols: Vec<usize> = if largest {
            (zc - b..zc).rev().collect()
        } else {
            (0..b).collect()
        };
        let mut x_new = vec![0.0f64; n * b];
        for (slot, &c) in cols.iter().enumerate() {
            for i in 0..n {
                let mut acc = 0.0f64;
                for p in 0..zc {
                    acc = zq[i * zc + p].mul_add(evecs[p * zc + c], acc);
                }
                x_new[i * b + slot] = acc;
            }
        }
        let _ = evals;
        state.x_prev = std::mem::take(&mut state.x);
        state.x = orthonormalize(&x_new, n, b);
        state.iters += 1;
    }
    // Final certified pairs.
    let ax = apply_block(op, &state.x, n, b);
    let lam = block_rayleigh(&state.x, &ax, n, b);
    let mut out = Vec::with_capacity(b);
    for j in 0..b {
        let mut v = vec![0.0f64; n];
        let mut rsq = 0.0f64;
        for i in 0..n {
            v[i] = state.x[i * b + j];
            let r = (-lam[j]).mul_add(v[i], ax[i * b + j]);
            rsq = r.mul_add(r, rsq);
        }
        out.push(EigenPair {
            value: lam[j],
            vector: v,
            residual: rsq.sqrt(),
        });
    }
    // Deterministic presentation: ascending (or descending for largest).
    out.sort_by(|a, c| {
        if largest {
            c.value.partial_cmp(&a.value).unwrap()
        } else {
            a.value.partial_cmp(&c.value).unwrap()
        }
    });
    out
}

fn apply_block<Op>(op: &Op, x: &[f64], n: usize, k: usize) -> Vec<f64>
where
    Op: Fn(&[f64], &mut [f64]),
{
    let mut out = vec![0.0f64; n * k];
    let mut col = vec![0.0f64; n];
    let mut res = vec![0.0f64; n];
    for j in 0..k {
        for i in 0..n {
            col[i] = x[i * k + j];
        }
        op(&col, &mut res);
        for i in 0..n {
            out[i * k + j] = res[i];
        }
    }
    out
}

fn block_rayleigh(x: &[f64], ax: &[f64], n: usize, b: usize) -> Vec<f64> {
    let mut lam = vec![0.0f64; b];
    for (j, l) in lam.iter_mut().enumerate() {
        let mut acc = 0.0f64;
        for i in 0..n {
            acc = x[i * b + j].mul_add(ax[i * b + j], acc);
        }
        *l = acc;
    }
    lam
}
