//! fs-wasm — a thin browser surface over FrankenSim's pure numerical leaves.
//!
//! Every function here runs the *real* kernel code — the same code the native
//! workspace compiles — just targeted at `wasm32-unknown-unknown`. No mocks,
//! no re-implementations: the in-browser showcase on the website is driven end
//! to end by these functions.
//!
//! Kernels surfaced: `fs-sparse` (sparse assembly + SpMV + CG), `fs-cheb`
//! (Chebyshev spectral + Orr–Sommerfeld), `fs-rand` (Philox streams + Sobol'
//! QMC), `fs-math` (strict elementary functions + error-free transforms),
//! `fs-ivl` (certified interval / Taylor-model arithmetic + exact geometric
//! predicates), `fs-ad` (forward-mode automatic differentiation), `fs-fft`
//! (real FFT), and `fs-la` (randomized SVD + symmetric eigensolve).
//!
//! The plain `pub fn`s below compile natively (rlib) and to wasm (cdylib). The
//! `#[wasm_bindgen]` layer at the bottom is compiled only for wasm32 and
//! exposes them to JavaScript as `Float64Array`-returning functions. Every
//! input is clamped to a safe range and every fallible kernel result is
//! folded to `NaN` / an empty vector — nothing here can trap at runtime.

use fs_ad::{second_directional, Real};
use fs_cheb::{orr_sommerfeld, Cheb1};
use fs_fft::RealFft;
use fs_ivl::{orient2d, Interval, Sign, TaylorModel1};
use fs_la::{eigen::jacobi_eigh, rand_nla::rsvd};
use fs_math::{det, eft::two_sum};
use fs_rand::{qmc::Sobol, StreamKey};
use fs_sparse::{Coo, Csr};

// ---------------------------------------------------------------------------
// Tier-2 demos: bespoke 3D/WebGL-oriented computations, each still running the
// real compiled kernel math. Split into topical modules for readability; every
// public entry point is re-exported at the crate root and re-surfaced through
// the `#[wasm_bindgen]` layer at the bottom of this file.
// ---------------------------------------------------------------------------
pub mod certified;
pub mod deep;
pub mod dynamics;
pub mod geom;
pub mod pde;

pub use certified::mandelbrot_certified;
pub use deep::{
    betti_shapes, cmaes_trace, cutfem_quadtree, cyclic_symmetry, ffd_deform, gp_regression,
    hodge_decomposition, krylov_convergence, navier_stokes_cavity, optimal_transport,
};
pub use dynamics::{ga_motor_orbit, lorenz_points, symplectic_vs_euler};
pub use geom::{marching_cubes, sdf_volume};
pub use pde::{fluid_frames, gray_scott_frames, topopt_frames, wave2d_frames};

/* ----------------------------------------------------------------------- */
/*  L1 · BEDROCK — sparse linear algebra: a real 2D Poisson solve           */
/* ----------------------------------------------------------------------- */

pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Assemble the 5-point Laplacian (SPD, `-Δ` up to the `1/h²` scale) on an
/// `n×n` interior grid of the unit square with zero Dirichlet boundaries.
pub(crate) fn laplacian_5pt(n: usize) -> Csr {
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
pub(crate) fn cg(a: &Csr, b: &[f64], maxit: usize, tol: f64) -> Vec<f64> {
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

/// The chosen smooth test function on `[-1, 1]`. `kind`:
/// 0 = Runge `1/(1+25x²)`, 1 = `sin(6x)`, 2 = `tanh(6x)`, 3 = `exp(-4x²)`.
/// All four resolve to machine precision well under 300 Chebyshev modes, so
/// the adaptive builder never exceeds its degree cap (no runtime trap).
fn cheb_testfn(kind: u32) -> impl Fn(f64) -> f64 {
    move |x: f64| match kind {
        0 => 1.0 / (1.0 + 25.0 * x * x),
        1 => (6.0 * x).sin(),
        2 => (6.0 * x).tanh(),
        _ => (-4.0 * x * x).exp(),
    }
}

/// Adaptive Chebyshev approximation of a chosen smooth test function on
/// `[-1,1]`, sampled at `samples` points, returning `[f, p, p']` interleaved
/// (`3 * samples`) so the browser can overlay the spectral fit and its exact
/// derivative on the truth. Real `fs-cheb` (DCT-based analysis + Clenshaw).
pub fn chebyshev_fit(kind: u32, samples_in: usize) -> Vec<f64> {
    // 256 modes resolve every listed test function; cap at 300 so the adaptive
    // builder's `n <= max_degree` guard never fires (it doubles 16→…→256).
    let max_degree = 300usize;
    let samples = samples_in.clamp(8, 1024);
    let f = cheb_testfn(kind);
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

/// Magnitudes of the Chebyshev coefficients (spectral decay) of the same
/// smooth test function — the classic "spectral accuracy" fingerprint, one
/// value per retained mode (length = resolved degree + 1). Real `fs-cheb`.
pub fn chebyshev_spectrum(kind: u32) -> Vec<f64> {
    let max_degree = 300usize;
    let f = cheb_testfn(kind);
    Cheb1::build(&f, -1.0, 1.0, max_degree)
        .coeffs()
        .iter()
        .map(|c| c.abs())
        .collect()
}

/* ----------------------------------------------------------------------- */
/*  fs-ivl — a PROVEN enclosure via a Taylor model                          */
/* ----------------------------------------------------------------------- */

/// Certified enclosure of the nonlinear function `f(x) = exp(sin(x))` over a
/// domain `[center-radius, center+radius]`, built with a rigorous
/// [`TaylorModel1`] (polynomial part + outward-rounded Lagrange remainder).
///
/// Output layout (length `K + 2`, where `K = order`):
/// - `[0 .. K]`  — rigorous remainder WIDTH at expansion orders `1..=K`. This
///   is a *proven* upper bound on the truncation+rounding error, and it
///   shrinks like `O(w^{n+1})` as the order climbs (the convergence curve).
/// - `[K]`       — certified lower bound of `f` over the whole domain.
/// - `[K+1]`     — certified upper bound of `f` over the whole domain.
///
/// The last two bracket the TRUE range with a machine-checkable guarantee —
/// not a float estimate. Real `fs-ivl` interval + Taylor-model arithmetic.
pub fn taylor_bound(center: f64, radius: f64, order_in: usize) -> Vec<f64> {
    let c = center.clamp(-3.0, 3.0);
    let r = radius.clamp(1.0e-3, 2.0);
    let k = order_in.clamp(1, 16);
    let dom = Interval::new(c - r, c + r);
    let mut out = Vec::with_capacity(k + 2);
    let mut lo = f64::NAN;
    let mut hi = f64::NAN;
    for order in 1..=k {
        let x = TaylorModel1::variable(dom, order);
        let f = x.sin().exp(); // exp∘sin — genuinely nonlinear, non-monotone
        out.push(f.remainder().width());
        if order == k {
            let b = f.bound();
            lo = b.lo();
            hi = b.hi();
        }
    }
    out.push(lo);
    out.push(hi);
    out
}

/* ----------------------------------------------------------------------- */
/*  fs-ad — exact automatic differentiation vs finite differences           */
/* ----------------------------------------------------------------------- */

/// The differentiation test function `f(x) = sin(3x)·exp(-x²/4)`, written once
/// generic over `Real` so it runs on `f64` (primal) and on nested duals
/// (exact derivatives) from the SAME source — the fs-ad contract.
fn ad_testfn<T: Real>(x: T) -> T {
    let three = T::from_f64(3.0);
    let quarter = T::from_f64(0.25);
    (three * x).sin() * (-(x * x * quarter)).exp()
}

/// Exact value + first + second derivative of `f(x) = sin(3x)·exp(-x²/4)` at
/// `samples` points across `[xmin, xmax]`, via forward-mode AD through NESTED
/// dual numbers (fs-ad `second_directional`). No finite differences — these
/// derivatives are exact to machine precision.
///
/// Output layout: `4 * samples` values, interleaved per point as
/// `[x, f(x), f'(x), f''(x)]`.
pub fn autodiff_derivatives(xmin: f64, xmax: f64, samples_in: usize) -> Vec<f64> {
    let a = xmin.clamp(-20.0, 20.0);
    let b = xmax.clamp(-20.0, 20.0);
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    let samples = samples_in.clamp(2, 2048);
    let mut out = Vec::with_capacity(samples * 4);
    for i in 0..samples {
        let t = i as f64 / ((samples - 1) as f64);
        let x = a + (b - a) * t;
        let (f, d1, d2) = second_directional([x], [1.0], |v| ad_testfn(v[0]));
        out.push(x);
        out.push(f);
        out.push(d1);
        out.push(d2);
    }
    out
}

/// Finite-difference error of `f'(x0)` across a log-spaced sweep of step sizes
/// `h`, contrasted with automatic differentiation (which is exact). At each
/// `h`, the central difference `(f(x0+h) − f(x0−h)) / 2h` is compared against
/// the AD-exact derivative; the error curve is U-shaped — truncation error
/// falls as `h²`, then rounding error blows up as `1/h`, bottoming near
/// `√ε ≈ 1e-8`. AD has none of that: it is exact for every `h`.
///
/// Output layout: `2 * steps` values, interleaved as `[h, |fd_error|]`,
/// with `h` running log-spaced from `1e-1` down to `1e-13`. (AD's error is
/// ~machine epsilon and flat — the demo draws it as a reference line.)
pub fn finite_difference_error(x0: f64, steps_in: usize) -> Vec<f64> {
    let x = x0.clamp(-10.0, 10.0);
    let steps = steps_in.clamp(4, 400);
    let (_, dexact, _) = second_directional([x], [1.0], |v| ad_testfn(v[0]));
    let mut out = Vec::with_capacity(steps * 2);
    for i in 0..steps {
        let t = i as f64 / ((steps - 1) as f64);
        let h = 10.0f64.powf(-1.0 - 12.0 * t); // 1e-1 → 1e-13
        let fd = (ad_testfn(x + h) - ad_testfn(x - h)) / (2.0 * h);
        let err = (fd - dexact).abs();
        out.push(h);
        out.push(err.max(1.0e-18)); // floor keeps log-plots finite
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  fs-la — randomized SVD of a smooth (numerically low-rank) kernel matrix  */
/* ----------------------------------------------------------------------- */

/// Randomized SVD (fs-la `rsvd`) of an `n×n` symmetric RBF kernel matrix
/// `A[i,j] = exp(-((i-j)/bw)²) + small symmetric noise` — a smooth,
/// numerically low-rank operator with fast-decaying spectrum. Returns the
/// TRUE top-`k` singular values (from a dense symmetric eigensolve,
/// `jacobi_eigh`) alongside the RANDOMIZED estimates, plus the rank-`k`
/// reconstruction error — randomized numerical linear algebra, in the browser.
///
/// Output layout (length `2*k + 1`):
/// - `[0 .. k]`     — true singular values, descending (|eigenvalues| of A).
/// - `[k .. 2k]`    — randomized singular values from `rsvd`, descending.
/// - `[2k]`         — relative Frobenius reconstruction error
///   `‖A − U diag(σ) Vᵀ‖_F / ‖A‖_F` at rank `k`.
pub fn randomized_svd(n_in: usize, rank_in: usize, seed_in: u32) -> Vec<f64> {
    let n = n_in.clamp(6, 64);
    let k = rank_in.clamp(1, n.min(16));
    let seed = seed_in as u64;
    // Build the symmetric RBF kernel + tiny symmetric noise (low-rank + noise).
    let bw = (n as f64 / 6.0).max(1.0);
    let noise_amp = 1.0e-3;
    let mut nz = StreamKey {
        seed,
        kernel: 0xF5_01,
        tile: 0,
    }
    .stream();
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in i..n {
            let d = (i as f64 - j as f64) / bw;
            let val = det::exp(-(d * d)) + noise_amp * nz.next_normal();
            a[i * n + j] = val;
            a[j * n + i] = val;
        }
    }
    // TRUE spectrum: symmetric eigenvalues → singular values are |λ|.
    let (evals, _) = jacobi_eigh(&a, n);
    let mut sv: Vec<f64> = evals.iter().map(|v| v.abs()).collect();
    sv.sort_by(|x, y| y.total_cmp(x)); // descending
    // RANDOMIZED SVD.
    let (u, sigma, v, _report) = rsvd(&a, n, n, k, 5, 2, seed ^ 0x9E37_79B9);
    let keep = sigma.len();
    // Rank-k reconstruction error (Frobenius, relative).
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let mut ak = 0.0f64;
            for l in 0..keep {
                ak += u[i * keep + l] * sigma[l] * v[j * keep + l];
            }
            let d = a[i * n + j] - ak;
            num += d * d;
            den += a[i * n + j] * a[i * n + j];
        }
    }
    let relerr = if den > 0.0 { (num / den).sqrt() } else { f64::NAN };
    let mut out = Vec::with_capacity(2 * k + 1);
    for i in 0..k {
        out.push(sv.get(i).copied().unwrap_or(0.0));
    }
    for i in 0..k {
        out.push(sigma.get(i).copied().unwrap_or(0.0));
    }
    out.push(relerr);
    out
}

/* ----------------------------------------------------------------------- */
/*  fs-fft — power spectrum of a synthetic multi-tone signal                 */
/* ----------------------------------------------------------------------- */

/// Power spectrum of a synthetic multi-tone signal plus Gaussian noise,
/// length `n` rounded to a power of two, computed with the real packed
/// real-input FFT (fs-fft `RealFft`). Three pure tones are injected; the
/// spectrum shows clean peaks at exactly their bins above the noise floor.
///
/// Output layout: `n/2 + 1` values — the per-bin power `|X[k]|²` for
/// `k = 0 ..= n/2` (bin `k` is frequency `k` cycles across the window).
pub fn fft_power_spectrum(n_in: usize, seed_in: u32) -> Vec<f64> {
    let mut n = n_in.clamp(64, 4096).next_power_of_two();
    if n > 4096 {
        n = 4096;
    }
    let seed = seed_in as u64;
    // Injected tones (bins strictly below n/2), with amplitudes.
    let k1 = (n / 32).max(2);
    let k2 = (n / 12).max(5);
    let k3 = (n / 6).max(11);
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut nz = StreamKey {
        seed,
        kernel: 0xF5_02,
        tile: 0,
    }
    .stream();
    let mut sig = vec![0.0f64; n];
    for (j, s) in sig.iter_mut().enumerate() {
        let jf = j as f64;
        let nf = n as f64;
        *s = 1.0 * det::cos(two_pi * k1 as f64 * jf / nf)
            + 0.6 * det::cos(two_pi * k2 as f64 * jf / nf)
            + 0.8 * det::sin(two_pi * k3 as f64 * jf / nf)
            + 0.15 * nz.next_normal();
    }
    let spectrum = RealFft::new(n).forward(&sig);
    spectrum.iter().map(|c| c.norm_sq()).collect()
}

/* ----------------------------------------------------------------------- */
/*  fs-la — vibration modes of the 1D Laplacian (symmetric eigensolve)       */
/* ----------------------------------------------------------------------- */

/// The lowest `k` eigenpairs of the 1D Dirichlet Laplacian (the symmetric
/// tridiagonal `[-1, 2, -1]` on `n` interior nodes), from the dense
/// cyclic-Jacobi symmetric eigensolver (fs-la `jacobi_eigh`). The
/// eigenvectors are the exact discrete standing-wave / vibration modes — a
/// plucked string's harmonics — ready to animate.
///
/// Output layout (length `k + k*n`):
/// - `[0 .. k]` — the `k` smallest eigenvalues, ascending.
/// - then `k` blocks of `n` values each: eigenvector `m` occupies
///   `[k + m*n .. k + (m+1)*n]` (mode shape sampled at the `n` nodes).
pub fn laplacian_modes(n_in: usize, k_in: usize) -> Vec<f64> {
    let n = n_in.clamp(4, 80);
    let k = k_in.clamp(1, n.min(8));
    // Dense symmetric tridiagonal 1D Laplacian.
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        a[i * n + i] = 2.0;
        if i + 1 < n {
            a[i * n + i + 1] = -1.0;
            a[(i + 1) * n + i] = -1.0;
        }
    }
    let (evals, evecs) = jacobi_eigh(&a, n);
    let mut out = Vec::with_capacity(k + k * n);
    for m in 0..k {
        out.push(evals.get(m).copied().unwrap_or(f64::NAN));
    }
    // eigenvectors are columns of the row-major evecs matrix.
    for m in 0..k {
        for i in 0..n {
            out.push(evecs[i * n + m]);
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  fs-rand — QMC vs Monte-Carlo convergence (π via the unit disk)           */
/* ----------------------------------------------------------------------- */

/// Convergence of a π estimate (fraction of points inside the unit quarter
/// disk, ×4) as the sample count grows, comparing plain pseudorandom Monte
/// Carlo (fs-rand Philox stream) against scrambled Sobol' quasi-Monte-Carlo
/// (fs-rand `Sobol::scrambled`, true Owen scrambling). QMC converges faster
/// than the `1/√N` MC rate — visibly, in the browser.
///
/// Output layout: `3 * steps` values, interleaved as `[N, mc_error, qmc_error]`
/// for `N = 2^6, 2^7, …, 2^max_log2` (so `steps = max_log2 − 5`). Errors are
/// `|estimate − π|`.
pub fn qmc_vs_mc(max_log2_in: usize, seed_in: u32) -> Vec<f64> {
    let lg = max_log2_in.clamp(8, 18);
    let seed = seed_in as u64;
    let nmax = 1u32 << lg;
    let sob = Sobol::scrambled(2, seed);
    let mut mc = StreamKey {
        seed,
        kernel: 0xF5_03,
        tile: 0,
    }
    .stream();
    let mut buf = [0.0f64; 2];
    let mut qmc_in = 0u64;
    let mut mc_in = 0u64;
    let pi = std::f64::consts::PI;
    let mut out = Vec::new();
    for i in 0..nmax {
        // QMC point.
        sob.point(i, &mut buf);
        if buf[0] * buf[0] + buf[1] * buf[1] < 1.0 {
            qmc_in += 1;
        }
        // MC point (two fresh uniforms).
        let mx = mc.next_f64();
        let my = mc.next_f64();
        if mx * mx + my * my < 1.0 {
            mc_in += 1;
        }
        // Record at every power-of-two checkpoint from 2^6 up.
        let count = i + 1;
        if count >= 64 && count.is_power_of_two() {
            let nf = count as f64;
            let mc_est = 4.0 * mc_in as f64 / nf;
            let qmc_est = 4.0 * qmc_in as f64 / nf;
            out.push(nf);
            out.push((mc_est - pi).abs());
            out.push((qmc_est - pi).abs());
        }
    }
    out
}

/* ----------------------------------------------------------------------- */
/*  fs-ivl — provably robust convex hull via EXACT geometric predicates      */
/* ----------------------------------------------------------------------- */

/// The exact convex hull of every integer lattice point inside a disk of the
/// given radius (a maximally-degenerate input: hull edges carry many exactly
/// collinear points that break naive floating-point hulls). Turn tests use
/// the adaptive-precision, provably-exact `orient2d` predicate (fs-ivl), so
/// the monotone-chain hull is correct with zero risk of the misordered /
/// self-crossing polygons a float cross-product produces on collinear input.
///
/// Output layout (a single flat array the viz slices):
/// - `[0]`                         — `P`, the number of input points.
/// - `[1 .. 1+2P]`                 — input points, interleaved `x,y`.
/// - `[1+2P]`                      — `H`, the number of hull vertices.
/// - `[2+2P .. 2+2P+2H]`           — hull vertices in counter-clockwise order,
///   interleaved `x,y` (strict hull: no collinear points on the edges).
pub fn robust_hull(radius_in: usize) -> Vec<f64> {
    let r = radius_in.clamp(2, 12) as i64;
    let r2 = r * r;
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for x in -r..=r {
        for y in -r..=r {
            if x * x + y * y <= r2 {
                pts.push((x as f64, y as f64));
            }
        }
    }
    let hull = convex_hull_exact(&pts);
    let mut out = Vec::with_capacity(1 + 2 * pts.len() + 1 + 2 * hull.len());
    out.push(pts.len() as f64);
    for &(x, y) in &pts {
        out.push(x);
        out.push(y);
    }
    out.push(hull.len() as f64);
    for &(x, y) in &hull {
        out.push(x);
        out.push(y);
    }
    out
}

/// Andrew's monotone-chain convex hull using the EXACT `orient2d` predicate
/// for every turn decision (strict hull — collinear boundary points dropped).
fn convex_hull_exact(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut p = points.to_vec();
    p.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    p.dedup();
    if p.len() <= 2 {
        return p;
    }
    // A left turn (strictly CCW) keeps the middle vertex; anything else pops.
    let keep = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| -> bool {
        orient2d([o.0, o.1], [a.0, a.1], [b.0, b.1]) == Sign::Positive
    };
    let mut lower: Vec<(f64, f64)> = Vec::new();
    for &pt in &p {
        while lower.len() >= 2 {
            let o = lower[lower.len() - 2];
            let a = lower[lower.len() - 1];
            if keep(o, a, pt) {
                break;
            }
            lower.pop();
        }
        lower.push(pt);
    }
    let mut upper: Vec<(f64, f64)> = Vec::new();
    for &pt in p.iter().rev() {
        while upper.len() >= 2 {
            let o = upper[upper.len() - 2];
            let a = upper[upper.len() - 1];
            if keep(o, a, pt) {
                break;
            }
            upper.pop();
        }
        upper.push(pt);
    }
    // Drop each chain's last point (it is the first of the other).
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/* ----------------------------------------------------------------------- */
/*  fs-math — error-free transforms: exact vs naive summation                */
/* ----------------------------------------------------------------------- */

/// Deterministic compensated summation vs naive accumulation on a
/// catastrophically-cancelling series `[B, 1, 1, …, 1 (×count), −B]`, whose
/// exact sum is `count`. Naive left-to-right addition loses every `1` under
/// the huge `B` and returns `0`; the compensated sum — built on fs-math's
/// `two_sum` error-free transform (Kahan–Babuška) — recovers `count` exactly.
///
/// Output layout (length 5):
/// `[naive_sum, compensated_sum, true_value, |naive_error|, |compensated_error|]`.
pub fn compensated_sum(count_in: usize, log10_big_in: i32) -> Vec<f64> {
    let count = count_in.clamp(1, 4_000_000);
    let log10_big = log10_big_in.clamp(6, 300);
    let big = 10.0f64.powi(log10_big);
    let truth = count as f64;

    let mut xs = Vec::with_capacity(count + 2);
    xs.push(big);
    for _ in 0..count {
        xs.push(1.0);
    }
    xs.push(-big);

    // Naive left-to-right accumulation.
    let mut naive = 0.0f64;
    for &x in &xs {
        naive += x;
    }
    // Compensated: accumulate the exact rounding error via two_sum.
    let mut s = 0.0f64;
    let mut c = 0.0f64;
    for &x in &xs {
        let (t, e) = two_sum(s, x);
        c += e;
        s = t;
    }
    let comp = s + c;

    vec![
        naive,
        comp,
        truth,
        (naive - truth).abs(),
        (comp - truth).abs(),
    ]
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
    pub fn chebyshev_fit(kind: u32, samples: usize) -> Vec<f64> {
        super::chebyshev_fit(kind, samples)
    }

    #[wasm_bindgen]
    pub fn chebyshev_spectrum(kind: u32) -> Vec<f64> {
        super::chebyshev_spectrum(kind)
    }

    #[wasm_bindgen]
    pub fn taylor_bound(center: f64, radius: f64, order: usize) -> Vec<f64> {
        super::taylor_bound(center, radius, order)
    }

    #[wasm_bindgen]
    pub fn autodiff_derivatives(xmin: f64, xmax: f64, samples: usize) -> Vec<f64> {
        super::autodiff_derivatives(xmin, xmax, samples)
    }

    #[wasm_bindgen]
    pub fn finite_difference_error(x0: f64, steps: usize) -> Vec<f64> {
        super::finite_difference_error(x0, steps)
    }

    #[wasm_bindgen]
    pub fn randomized_svd(n: usize, rank: usize, seed: u32) -> Vec<f64> {
        super::randomized_svd(n, rank, seed)
    }

    #[wasm_bindgen]
    pub fn fft_power_spectrum(n: usize, seed: u32) -> Vec<f64> {
        super::fft_power_spectrum(n, seed)
    }

    #[wasm_bindgen]
    pub fn laplacian_modes(n: usize, k: usize) -> Vec<f64> {
        super::laplacian_modes(n, k)
    }

    #[wasm_bindgen]
    pub fn qmc_vs_mc(max_log2: usize, seed: u32) -> Vec<f64> {
        super::qmc_vs_mc(max_log2, seed)
    }

    #[wasm_bindgen]
    pub fn robust_hull(radius: usize) -> Vec<f64> {
        super::robust_hull(radius)
    }

    #[wasm_bindgen]
    pub fn compensated_sum(count: usize, log10_big: i32) -> Vec<f64> {
        super::compensated_sum(count, log10_big)
    }

    /* ------------------------------------------------------------------- */
    /*  Tier-2 demos (3D / WebGL): fields, meshes, trajectories             */
    /* ------------------------------------------------------------------- */

    #[wasm_bindgen]
    pub fn topopt_frames(nx: usize, ny: usize, iters: usize, volfrac: f64) -> Vec<f64> {
        super::topopt_frames(nx, ny, iters, volfrac)
    }

    #[wasm_bindgen]
    pub fn marching_cubes(res: usize, kind: u32, iso: f64) -> Vec<f64> {
        super::marching_cubes(res, kind, iso)
    }

    #[wasm_bindgen]
    pub fn sdf_volume(res: usize, kind: u32, t: f64) -> Vec<f64> {
        super::sdf_volume(res, kind, t)
    }

    #[wasm_bindgen]
    pub fn ga_motor_orbit(n_points: usize, steps: usize) -> Vec<f64> {
        super::ga_motor_orbit(n_points, steps)
    }

    #[wasm_bindgen]
    pub fn symplectic_vs_euler(steps: usize, dt: f64) -> Vec<f64> {
        super::symplectic_vs_euler(steps, dt)
    }

    #[wasm_bindgen]
    pub fn lorenz_points(steps: usize, dt: f64, rho: f64) -> Vec<f64> {
        super::lorenz_points(steps, dt, rho)
    }

    #[wasm_bindgen]
    pub fn wave2d_frames(n: usize, frames: usize, steps_per_frame: usize) -> Vec<f64> {
        super::wave2d_frames(n, frames, steps_per_frame)
    }

    #[wasm_bindgen]
    pub fn gray_scott_frames(n: usize, frames: usize, feed: f64, kill: f64) -> Vec<f64> {
        super::gray_scott_frames(n, frames, feed, kill)
    }

    #[wasm_bindgen]
    pub fn mandelbrot_certified(
        w: usize,
        h: usize,
        cx: f64,
        cy: f64,
        scale: f64,
        maxiter: usize,
    ) -> Vec<f64> {
        super::mandelbrot_certified(w, h, cx, cy, scale, maxiter)
    }

    #[wasm_bindgen]
    pub fn fluid_frames(n: usize, frames: usize) -> Vec<f64> {
        super::fluid_frames(n, frames)
    }

    /// A build stamp so the page can prove it's running the real engine.
    #[wasm_bindgen]
    pub fn engine() -> String {
        "fs-wasm · FrankenSim numerical kernels (fs-sparse · fs-cheb · fs-rand · fs-ivl · \
         fs-ad · fs-fft · fs-la · fs-ga · fs-math)"
            .into()
    }
}
