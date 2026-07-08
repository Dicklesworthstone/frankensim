//! fs-wasm — a thin browser surface over FrankenSim's pure numerical leaves.
//!
//! Every function here runs the *real* kernel code (`fs-sparse`, `fs-cheb`,
//! `fs-rand`, `fs-math`) — the same code the native workspace compiles — just
//! targeted at `wasm32-unknown-unknown`. No mocks, no re-implementations: the
//! browser demos on the website are driven by these functions.
//!
//! The plain `pub fn`s below compile natively (rlib) and to wasm (cdylib). The
//! `#[wasm_bindgen]` layer at the bottom is compiled only for wasm32 and exposes
//! them to JavaScript as `Float64Array`-returning functions.

use fs_cheb::{orr_sommerfeld, Cheb1};
use fs_sparse::{Coo, Csr};

/* ----------------------------------------------------------------------- */
/*  L1 · BEDROCK — sparse linear algebra: a real 2D Poisson solve           */
/* ----------------------------------------------------------------------- */

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Assemble the 5-point Laplacian (SPD, `-Δ` up to the `1/h²` scale) on an
/// `n×n` interior grid of the unit square with zero Dirichlet boundaries.
fn laplacian_5pt(n: usize) -> Csr {
    let m = n * n;
    let mut coo = Coo::new(m, m);
    let idx = |i: usize, j: usize| i * n + j;
    for i in 0..n {
        for j in 0..n {
            let k = idx(i, j);
            coo.push(k, k, 4.0);
            if i > 0 {
                coo.push(k, idx(i - 1, j), -1.0);
            }
            if i + 1 < n {
                coo.push(k, idx(i + 1, j), -1.0);
            }
            if j > 0 {
                coo.push(k, idx(i, j - 1), -1.0);
            }
            if j + 1 < n {
                coo.push(k, idx(i, j + 1), -1.0);
            }
        }
    }
    coo.assemble()
}

/// Conjugate gradients against a `Csr` operator (matrix-free via `spmv`).
fn cg(a: &Csr, b: &[f64], maxit: usize, tol: f64) -> Vec<f64> {
    let m = b.len();
    let mut x = vec![0.0f64; m];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut ap = vec![0.0f64; m];
    let mut rs = dot(&r, &r);
    let bnorm = dot(b, b).sqrt().max(1e-30);
    for _ in 0..maxit {
        a.spmv(&p, &mut ap);
        let denom = dot(&p, &ap);
        if denom.abs() < 1e-300 {
            break;
        }
        let alpha = rs / denom;
        for i in 0..m {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new = dot(&r, &r);
        if rs_new.sqrt() / bnorm < tol {
            break;
        }
        let beta = rs_new / rs;
        for i in 0..m {
            p[i] = r[i] + beta * p[i];
        }
        rs = rs_new;
    }
    x
}

/// Solve `-Δu = f` on an `n×n` grid (two Gaussian sources) with conjugate
/// gradients on the assembled Laplacian. Returns the field row-major (`n*n`).
/// Real `fs-sparse` assembly + SpMV + a matrix-free CG.
pub fn poisson2d(n_in: usize) -> Vec<f64> {
    let n = n_in.clamp(3, 110);
    let m = n * n;
    let a = laplacian_5pt(n);
    let h = 1.0 / (n as f64 + 1.0);
    let mut b = vec![0.0f64; m];
    for i in 0..n {
        for j in 0..n {
            let x = (i as f64 + 1.0) * h;
            let y = (j as f64 + 1.0) * h;
            let g = |cx: f64, cy: f64, s: f64| {
                (-(((x - cx).powi(2) + (y - cy).powi(2)) / (2.0 * s * s))).exp()
            };
            b[i * n + j] = (g(0.32, 0.34, 0.10) - 0.85 * g(0.68, 0.66, 0.12)) * h * h * 60.0;
        }
    }
    cg(&a, &b, 600, 1e-9)
}

/// Explicit heat diffusion of an initial two-spot field, stepped with the real
/// Laplacian SpMV. Returns `frames` snapshots concatenated (`frames * n*n`) so
/// the browser can animate the diffusion in time.
pub fn heat_frames(n_in: usize, frames_in: usize, steps_per_frame_in: usize) -> Vec<f64> {
    let n = n_in.clamp(3, 96);
    let frames = frames_in.clamp(1, 240);
    let spf = steps_per_frame_in.clamp(1, 40);
    let m = n * n;
    let a = laplacian_5pt(n);
    // Initial condition: a hot and a cold blob.
    let mut u = vec![0.0f64; m];
    let h = 1.0 / (n as f64 + 1.0);
    for i in 0..n {
        for j in 0..n {
            let x = (i as f64 + 1.0) * h;
            let y = (j as f64 + 1.0) * h;
            let g = |cx: f64, cy: f64, s: f64| {
                (-(((x - cx).powi(2) + (y - cy).powi(2)) / (2.0 * s * s))).exp()
            };
            u[i * n + j] = g(0.3, 0.3, 0.07) - g(0.7, 0.68, 0.08);
        }
    }
    let dt = 0.20; // A has ~4 on the diagonal → 8 spectral bound → dt < 0.25.
    let mut au = vec![0.0f64; m];
    let mut out = Vec::with_capacity(frames * m);
    for _ in 0..frames {
        out.extend_from_slice(&u);
        for _ in 0..spf {
            a.spmv(&u, &mut au); // Au approximates (-Δu)·something ≥ 0
            for i in 0..m {
                u[i] -= dt * au[i];
            }
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  L1 · BEDROCK — Chebyshev spectral: hydrodynamic stability               */
/* ----------------------------------------------------------------------- */

/// Maximum temporal growth rate of plane Poiseuille flow at `(Re, α)` via a
/// Chebyshev-collocation Orr–Sommerfeld solve. Positive ⇒ unstable. Real
/// `fs-cheb` spectral eigensolve — the physics behind "the spout that never
/// dribbles."
pub fn orr_sommerfeld_max_growth(re: f64, alpha: f64, n_in: usize) -> f64 {
    let n = n_in.clamp(16, 120);
    orr_sommerfeld::max_growth(re, alpha, n).unwrap_or(f64::NAN)
}

/// A growth-rate curve `max_growth(Re)` at fixed `α`, sampled over
/// `[re_min, re_max]` in `steps` points — for a neutral-stability plot.
pub fn orr_sommerfeld_curve(
    alpha: f64,
    n_in: usize,
    re_min: f64,
    re_max: f64,
    steps_in: usize,
) -> Vec<f64> {
    let n = n_in.clamp(16, 100);
    let steps = steps_in.clamp(2, 400);
    (0..steps)
        .map(|k| {
            let t = k as f64 / ((steps - 1) as f64);
            let re = re_min + (re_max - re_min) * t;
            orr_sommerfeld::max_growth(re, alpha, n).unwrap_or(f64::NAN)
        })
        .collect()
}

/// Adaptive Chebyshev approximation of a chosen test function on `[-1,1]`,
/// sampled at `samples` points, returning `[f, p, p']` concatenated so the
/// browser can show the spectral fit and its exact derivative. `kind`:
/// 0 = runge 1/(1+25x²), 1 = |x|·sin(6x), 2 = tanh(8x). Real `fs-cheb`.
pub fn chebyshev_fit(kind: u32, max_degree_in: usize, samples_in: usize) -> Vec<f64> {
    let max_degree = max_degree_in.clamp(2, 128);
    let samples = samples_in.clamp(8, 1024);
    let f = move |x: f64| match kind {
        0 => 1.0 / (1.0 + 25.0 * x * x),
        1 => x.abs() * (6.0 * x).sin(),
        _ => (8.0 * x).tanh(),
    };
    let cheb = Cheb1::build(&f, -1.0, 1.0, max_degree);
    let dcheb = cheb.differentiate();
    let mut out = Vec::with_capacity(samples * 3);
    for k in 0..samples {
        let x = -1.0 + 2.0 * (k as f64) / ((samples - 1) as f64);
        out.push(f(x));
        out.push(cheb.eval(x));
        out.push(dcheb.eval(x));
    }
    out
}

/// The magnitude of the Chebyshev coefficients (spectral decay) for the same
/// test function — the classic "spectral accuracy" fingerprint.
pub fn chebyshev_spectrum(kind: u32, max_degree_in: usize) -> Vec<f64> {
    let max_degree = max_degree_in.clamp(2, 128);
    let f = move |x: f64| match kind {
        0 => 1.0 / (1.0 + 25.0 * x * x),
        1 => x.abs() * (6.0 * x).sin(),
        _ => (8.0 * x).tanh(),
    };
    Cheb1::build(&f, -1.0, 1.0, max_degree)
        .coeffs()
        .iter()
        .map(|c| c.abs())
        .collect()
}

/* ----------------------------------------------------------------------- */
/*  The JavaScript boundary (wasm32 only)                                   */
/* ----------------------------------------------------------------------- */

#[cfg(target_arch = "wasm32")]
mod wasm {
    use wasm_bindgen::prelude::*;

    #[wasm_bindgen]
    pub fn poisson2d(n: usize) -> Vec<f64> {
        super::poisson2d(n)
    }

    #[wasm_bindgen]
    pub fn heat_frames(n: usize, frames: usize, steps_per_frame: usize) -> Vec<f64> {
        super::heat_frames(n, frames, steps_per_frame)
    }

    #[wasm_bindgen]
    pub fn orr_sommerfeld_max_growth(re: f64, alpha: f64, n: usize) -> f64 {
        super::orr_sommerfeld_max_growth(re, alpha, n)
    }

    #[wasm_bindgen]
    pub fn orr_sommerfeld_curve(
        alpha: f64,
        n: usize,
        re_min: f64,
        re_max: f64,
        steps: usize,
    ) -> Vec<f64> {
        super::orr_sommerfeld_curve(alpha, n, re_min, re_max, steps)
    }

    #[wasm_bindgen]
    pub fn chebyshev_fit(kind: u32, max_degree: usize, samples: usize) -> Vec<f64> {
        super::chebyshev_fit(kind, max_degree, samples)
    }

    #[wasm_bindgen]
    pub fn chebyshev_spectrum(kind: u32, max_degree: usize) -> Vec<f64> {
        super::chebyshev_spectrum(kind, max_degree)
    }

    /// A build stamp so the page can prove it's running the real engine.
    #[wasm_bindgen]
    pub fn engine() -> String {
        "fs-wasm · FrankenSim numerical kernels (fs-sparse · fs-cheb · fs-rand)".into()
    }
}
