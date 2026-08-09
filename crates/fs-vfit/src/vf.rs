//! Relaxed vector fitting (Gustavsen-Semlyen) over fs-la least squares.
//!
//! One iteration: with current conjugate-closed poles, solve the linear
//! weighted LS for the model residues (plus `d`, `e`) AND the residues
//! of the scalar sigma function `sigma(s) = d_tilde + sum c_tilde_k
//! phi_k(s)` in `sigma(s)*H(s) ~= N(s)`; the relocated poles are the
//! zeros of `sigma`, computed as eigenvalues of `A - B c_tilde /
//! d_tilde` on sigma's real block realization. The RELAXED
//! non-triviality row (weighted real-part sum of `sigma` pinned to the
//! sample count) replaces the classic `sigma(inf) = 1` normalization.
//!
//! Real arithmetic throughout: conjugate-pair poles use the two real
//! basis functions `phi = 1/(s-p) + 1/(s-conj p)` and
//! `phi' = i/(s-p) - i/(s-conj p)`, so LS unknowns are real and the
//! fitted model is real by construction.
//!
//! Sign constraints `d >= 0`, `e >= 0` (impedance-form asymptotic
//! passivity) are enforced INSIDE the solve as an exact active-set on
//! the two bound constraints: if the unconstrained optimum violates a
//! bound, the bound is active at the optimum of this convex problem, so
//! the column is removed and the system re-solved.

use crate::model::{PoleTerm, RationalModel};
use fs_la::eigen_complex::{EigFailure, eig};
use fs_la::factor::qr;
use fs_math::c64::C64;

/// Weighting presets for the LS rows (per-sample weight `w_i`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightPreset {
    /// `w_i = 1`.
    Uniform,
    /// `w_i = 1/|H_i|` — relative error; what makes impedance fits
    /// honest at ANTIRESONANCES (which set waveguide reflection nulls).
    InverseMagnitude,
    /// `w_i = 1/sqrt(|H_i|)` — a compromise emphasizing both peaks and
    /// valleys on log-spaced bands.
    LogBand,
}

impl WeightPreset {
    /// Human-readable label recorded in fit reports.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            WeightPreset::Uniform => "uniform",
            WeightPreset::InverseMagnitude => "inverse-magnitude",
            WeightPreset::LogBand => "log-band",
        }
    }

    fn weight(self, h: C64) -> f64 {
        let m = h.abs();
        match self {
            WeightPreset::Uniform => 1.0,
            WeightPreset::InverseMagnitude => {
                if m > 0.0 {
                    1.0 / m
                } else {
                    1.0
                }
            }
            WeightPreset::LogBand => {
                if m > 0.0 {
                    1.0 / fs_math::det::sqrt(m)
                } else {
                    1.0
                }
            }
        }
    }
}

/// Structural options for one vector-fitting run.
#[derive(Debug, Clone, Copy)]
pub struct FitOptions {
    /// Pole count (pairs count 2); the starting pole set has exactly
    /// this order.
    pub order: usize,
    /// Pole-relocation iterations (each is one LS + one eigensolve).
    pub iterations: usize,
    /// Row weighting preset (recorded in the report).
    pub weights: WeightPreset,
    /// Fit the improper `s*e` term (asymptotically linear responses);
    /// `false` pins `e = 0`.
    pub fit_e: bool,
    /// Fit the direct `d` term; `false` pins `d = 0`.
    pub fit_d: bool,
}

impl FitOptions {
    /// Defaults: 20 relocation iterations, inverse-magnitude weights,
    /// both direct terms fitted.
    #[must_use]
    pub fn new(order: usize) -> Self {
        FitOptions {
            order,
            iterations: 20,
            weights: WeightPreset::InverseMagnitude,
            fit_e: true,
            fit_d: true,
        }
    }
}

/// Everything a fit run reports besides the model itself.
#[derive(Debug, Clone)]
pub struct FitReport {
    /// Weighted root-mean-square misfit of the FINAL residue pass.
    pub weighted_rms: f64,
    /// Unweighted maximum absolute misfit over the samples.
    pub max_abs_error: f64,
    /// Weight preset label actually used.
    pub weights: &'static str,
    /// Pole-relocation iterations actually run.
    pub iterations_run: usize,
    /// Largest pole movement (absolute, rad/s) in the LAST relocation.
    pub final_pole_movement: f64,
}

/// A fitted model plus its report.
#[derive(Debug, Clone)]
pub struct FitOutcome {
    /// The conjugate-closed, stable fitted model.
    pub model: RationalModel,
    /// Fit diagnostics.
    pub report: FitReport,
}

/// Typed vector-fitting failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VfError {
    /// Fewer samples than unknowns (order + 2 direct terms, doubled
    /// for the sigma system).
    TooFewSamples {
        /// Samples provided.
        samples: usize,
        /// Minimum required.
        needed: usize,
    },
    /// Zero or non-finite frequency/response data.
    BadSample {
        /// Index of the offending sample.
        index: usize,
    },
    /// `order == 0` or odd order with too few real poles to place.
    BadOrder {
        /// Requested order.
        order: usize,
    },
    /// The sigma-zero eigensolve failed to converge.
    Eig(EigFailure),
    /// The relaxed d_tilde collapsed to ~0 (degenerate sigma).
    SigmaCollapse,
}

impl core::fmt::Display for VfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            VfError::TooFewSamples { samples, needed } => {
                write!(f, "vector fitting needs >= {needed} samples, got {samples}")
            }
            VfError::BadSample { index } => write!(
                f,
                "non-finite or non-positive-frequency sample at index {index}"
            ),
            VfError::BadOrder { order } => write!(f, "unusable model order {order}"),
            VfError::Eig(e) => write!(f, "pole-relocation eigensolve failed: {e}"),
            VfError::SigmaCollapse => write!(f, "relaxed sigma normalization collapsed"),
        }
    }
}

impl std::error::Error for VfError {}

/// Deterministic starting poles: weakly damped conjugate pairs with
/// imaginary parts log-spaced over the sampled band (the classic VF
/// initialization; odd orders add one real pole at the low edge).
#[must_use]
pub fn initial_poles(omega: &[f64], order: usize) -> Vec<PoleTerm> {
    let (lo, hi) = band_edges(omega);
    let mut terms = Vec::new();
    let pairs = order / 2;
    for k in 0..pairs {
        let t = if pairs == 1 {
            0.5
        } else {
            k as f64 / (pairs - 1) as f64
        };
        let w = lo * fs_math::det::exp(t * fs_math::det::ln(hi / lo));
        terms.push(PoleTerm::Pair {
            pole: C64::new(-w / 100.0, w),
            residue: C64::ZERO,
        });
    }
    if order % 2 == 1 {
        terms.push(PoleTerm::Real {
            pole: -lo,
            residue: 0.0,
        });
    }
    terms
}

fn band_edges(omega: &[f64]) -> (f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = 0.0f64;
    for &w in omega {
        if w > 0.0 {
            lo = lo.min(w);
            hi = hi.max(w);
        }
    }
    if !lo.is_finite() || hi <= 0.0 {
        (1.0, 10.0)
    } else if lo == hi {
        (lo * 0.5, hi * 2.0)
    } else {
        (lo, hi)
    }
}

/// Per-pole real basis values at `s = i*omega`: for each term, one
/// column (real pole) or two columns (pair), matching the coordinate
/// layout used everywhere in this crate.
fn basis_columns(terms: &[PoleTerm], s: C64) -> Vec<C64> {
    let mut cols = Vec::new();
    for t in terms {
        match *t {
            PoleTerm::Real { pole, .. } => {
                cols.push((s - C64::from_re(pole)).recip());
            }
            PoleTerm::Pair { pole, .. } => {
                let a = (s - pole).recip();
                let b = (s - pole.conj()).recip();
                cols.push(a + b);
                // i/(s-p) - i/(s-conj p)
                let diff = a - b;
                cols.push(C64::new(-diff.im, diff.re));
            }
        }
    }
    cols
}

/// Rebuild terms from real residue coordinates (same layout as
/// [`basis_columns`]).
fn terms_with_residues(terms: &[PoleTerm], coords: &[f64]) -> Vec<PoleTerm> {
    let mut out = Vec::with_capacity(terms.len());
    let mut at = 0usize;
    for t in terms {
        match *t {
            PoleTerm::Real { pole, .. } => {
                out.push(PoleTerm::Real {
                    pole,
                    residue: coords[at],
                });
                at += 1;
            }
            PoleTerm::Pair { pole, .. } => {
                out.push(PoleTerm::Pair {
                    pole,
                    residue: C64::new(coords[at], coords[at + 1]),
                });
                at += 2;
            }
        }
    }
    out
}

/// One full vector-fitting run: relocation iterations then a final
/// residue-only pass with `sigma == 1`.
///
/// # Errors
/// Typed [`VfError`] on degenerate input, order, eigensolve failure, or
/// sigma collapse. Never panics on finite input.
pub fn vector_fit(omega: &[f64], h: &[C64], opts: &FitOptions) -> Result<FitOutcome, VfError> {
    validate(omega, h, opts)?;
    let mut terms = initial_poles(omega, opts.order);
    let weights: Vec<f64> = h.iter().map(|&v| opts.weights.weight(v)).collect();
    let mut iterations_run = 0usize;
    let mut final_move = 0.0f64;
    for _ in 0..opts.iterations {
        let (new_terms, movement) = relocate_once(omega, h, &weights, &terms, opts)?;
        terms = new_terms;
        iterations_run += 1;
        final_move = movement;
        let (_, hi) = band_edges(omega);
        if movement < 1.0e-9 * hi {
            break;
        }
    }
    let (model, weighted_rms, max_abs) = residue_pass(omega, h, &weights, &terms, opts)?;
    Ok(FitOutcome {
        model,
        report: FitReport {
            weighted_rms,
            max_abs_error: max_abs,
            weights: opts.weights.label(),
            iterations_run,
            final_pole_movement: final_move,
        },
    })
}

fn validate(omega: &[f64], h: &[C64], opts: &FitOptions) -> Result<(), VfError> {
    if opts.order == 0 {
        return Err(VfError::BadOrder { order: opts.order });
    }
    for (i, (&w, v)) in omega.iter().zip(h).enumerate() {
        if !w.is_finite() || w <= 0.0 || !v.re.is_finite() || !v.im.is_finite() {
            return Err(VfError::BadSample { index: i });
        }
    }
    // Sigma system unknowns: order (model residues) + 2 (d, e) + order
    // (sigma residues) + 1 (relaxed d_tilde); each sample gives 2 real
    // rows and the relaxation adds 1.
    let needed = opts.order + 1 + (opts.order + 2).div_ceil(2);
    if omega.len() < needed || omega.len() != h.len() {
        return Err(VfError::TooFewSamples {
            samples: omega.len(),
            needed,
        });
    }
    Ok(())
}

/// One relocation: LS solve of the relaxed sigma system, then the new
/// poles as zeros of sigma. Returns the conjugate-closed, stability-
/// flipped new terms and the largest pole movement.
fn relocate_once(
    omega: &[f64],
    h: &[C64],
    weights: &[f64],
    terms: &[PoleTerm],
    opts: &FitOptions,
) -> Result<(Vec<PoleTerm>, f64), VfError> {
    let order: usize = terms.iter().map(PoleTerm::state_dim).sum();
    // Residue coords + optional d + optional e; a pinned term gets NO
    // column (an all-zero column would make R exactly singular).
    let n_model = order + usize::from(opts.fit_d) + usize::from(opts.fit_e);
    let d_col = order;
    let e_col = order + usize::from(opts.fit_d);
    let n_sigma = order + 1; // sigma residue coords + relaxed d_tilde
    let ncols = n_model + n_sigma;
    let nrows = 2 * omega.len() + 1;
    let mut a = vec![0.0f64; nrows * ncols];
    let mut rhs = vec![0.0f64; nrows];
    for (i, ((&w, &hv), &wt)) in omega.iter().zip(h).zip(weights).enumerate() {
        let s = C64::new(0.0, w);
        let cols = basis_columns(terms, s);
        let (r_re, r_im) = (2 * i, 2 * i + 1);
        for (k, &phi) in cols.iter().enumerate() {
            a[r_re * ncols + k] = wt * phi.re;
            a[r_im * ncols + k] = wt * phi.im;
            // Sigma columns multiply -H(s).
            let m = -(hv * phi);
            a[r_re * ncols + n_model + k] = wt * m.re;
            a[r_im * ncols + n_model + k] = wt * m.im;
        }
        if opts.fit_d {
            a[r_re * ncols + d_col] = wt;
        }
        if opts.fit_e {
            // s*e at s = i*w is purely imaginary.
            a[r_im * ncols + e_col] = wt * w;
        }
        // Relaxed d_tilde column multiplies -H(s).
        a[r_re * ncols + ncols - 1] = wt * -hv.re;
        a[r_im * ncols + ncols - 1] = wt * -hv.im;
        rhs[r_re] = wt * hv.re;
        rhs[r_im] = wt * hv.im;
    }
    // Relaxed non-triviality row: sum_i Re(sigma(s_i)) = len (weighted
    // by the mean weight so the row's scale matches the block).
    let mean_wt = weights.iter().sum::<f64>() / weights.len() as f64;
    let last = nrows - 1;
    for (i, &w) in omega.iter().enumerate() {
        let s = C64::new(0.0, w);
        let cols = basis_columns(terms, s);
        for (k, &phi) in cols.iter().enumerate() {
            a[last * ncols + n_model + k] += mean_wt * phi.re;
        }
        a[last * ncols + ncols - 1] += mean_wt;
        let _ = i;
    }
    rhs[last] = mean_wt * omega.len() as f64;
    let x = qr(&a, nrows, ncols).solve_ls(&rhs);
    let d_tilde = x[ncols - 1];
    if d_tilde.abs() < 1.0e-12 {
        return Err(VfError::SigmaCollapse);
    }
    let sigma_coords = &x[n_model..n_model + order];
    // Zeros of sigma = eig(A - B * c_tilde / d_tilde) on sigma's real
    // block realization (same block layout as RationalModel).
    let sigma_model = RationalModel {
        terms: terms_with_residues(terms, sigma_coords),
        d: d_tilde,
        e: 0.0,
    };
    let ss = sigma_model.state_space();
    let n = ss.n;
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = C64::from_re(ss.a[i * n + j] - ss.b[i] * ss.c[j] / d_tilde);
        }
    }
    let mut zeros = eig(&m, n).map_err(VfError::Eig)?;
    // Stability flip: reflect any right-half-plane zero.
    for z in &mut zeros {
        if z.re > 0.0 {
            z.re = -z.re;
        }
        if z.re == 0.0 {
            // A marginal pole would make the next LS singular at a
            // sample frequency; nudge into the left half-plane.
            z.re = -1.0e-6 * (1.0 + z.im.abs());
        }
    }
    let new_terms = conjugate_close(&zeros);
    let movement = pole_movement(terms, &new_terms);
    Ok((new_terms, movement))
}

/// Group a conjugate-symmetric eigenvalue list into canonical terms:
/// near-real eigenvalues become real poles; the rest pair up by sorted
/// order (the solver's canonical (re, im) ordering makes the pairing
/// deterministic).
fn conjugate_close(zeros: &[C64]) -> Vec<PoleTerm> {
    let mut real_poles: Vec<f64> = Vec::new();
    let mut upper: Vec<C64> = Vec::new();
    for &z in zeros {
        let imag_scale = 1.0e-9 * (1.0 + z.re.abs());
        if z.im.abs() <= imag_scale {
            real_poles.push(z.re);
        } else if z.im > 0.0 {
            upper.push(z);
        }
        // Lower-half members are the implied conjugates: dropped.
    }
    // An odd leftover (numerical asymmetry) is folded to real.
    let lower_count = zeros.len() - real_poles.len() - upper.len();
    if lower_count != upper.len() {
        // Asymmetric pairing: rebuild conservatively by magnitude.
        // (Deterministic fallback; exercised only under numerical
        // degeneracy.)
        for &z in zeros {
            let imag_scale = 1.0e-9 * (1.0 + z.re.abs());
            if z.im < -imag_scale
                && !upper
                    .iter()
                    .any(|u| (u.conj() - z).abs() <= 1.0e-6 * z.abs())
            {
                upper.push(z.conj());
            }
        }
    }
    real_poles.sort_by(f64::total_cmp);
    upper.sort_by(|x, y| x.im.total_cmp(&y.im).then(x.re.total_cmp(&y.re)));
    let mut terms: Vec<PoleTerm> = real_poles
        .into_iter()
        .map(|p| PoleTerm::Real {
            pole: p,
            residue: 0.0,
        })
        .collect();
    terms.extend(upper.into_iter().map(|p| PoleTerm::Pair {
        pole: p,
        residue: C64::ZERO,
    }));
    terms
}

fn pole_movement(old: &[PoleTerm], new: &[PoleTerm]) -> f64 {
    // Greedy nearest-match movement metric (diagnostic only).
    let expand = |ts: &[PoleTerm]| -> Vec<C64> {
        let mut v = Vec::new();
        for t in ts {
            match *t {
                PoleTerm::Real { pole, .. } => v.push(C64::from_re(pole)),
                PoleTerm::Pair { pole, .. } => v.push(pole),
            }
        }
        v
    };
    let a = expand(old);
    let b = expand(new);
    let mut worst = 0.0f64;
    for p in &a {
        let mut best = f64::INFINITY;
        for q in &b {
            best = best.min((*p - *q).abs());
        }
        if best.is_finite() {
            worst = worst.max(best);
        }
    }
    worst
}

/// Final residue-only pass (`sigma == 1`) with the `d >= 0`, `e >= 0`
/// bound constraints resolved by exact active-set on the two bounds.
fn residue_pass(
    omega: &[f64],
    h: &[C64],
    weights: &[f64],
    terms: &[PoleTerm],
    opts: &FitOptions,
) -> Result<(RationalModel, f64, f64), VfError> {
    let order: usize = terms.iter().map(PoleTerm::state_dim).sum();
    // Try with both direct columns (when enabled); clamp any that come
    // out negative and re-solve without them.
    let mut use_d = opts.fit_d;
    let mut use_e = opts.fit_e;
    loop {
        let (coords, d, e) = residue_solve(omega, h, weights, terms, order, use_d, use_e);
        if d < 0.0 {
            use_d = false;
            continue;
        }
        if e < 0.0 {
            use_e = false;
            continue;
        }
        let model = RationalModel {
            terms: terms_with_residues(terms, &coords),
            d,
            e,
        };
        let (wrms, maxabs) = fit_errors(omega, h, weights, &model);
        return Ok((model, wrms, maxabs));
    }
}

fn residue_solve(
    omega: &[f64],
    h: &[C64],
    weights: &[f64],
    terms: &[PoleTerm],
    order: usize,
    use_d: bool,
    use_e: bool,
) -> (Vec<f64>, f64, f64) {
    let ncols = order + usize::from(use_d) + usize::from(use_e);
    let nrows = 2 * omega.len();
    let mut a = vec![0.0f64; nrows * ncols];
    let mut rhs = vec![0.0f64; nrows];
    for (i, ((&w, &hv), &wt)) in omega.iter().zip(h).zip(weights).enumerate() {
        let s = C64::new(0.0, w);
        let cols = basis_columns(terms, s);
        let (r_re, r_im) = (2 * i, 2 * i + 1);
        for (k, &phi) in cols.iter().enumerate() {
            a[r_re * ncols + k] = wt * phi.re;
            a[r_im * ncols + k] = wt * phi.im;
        }
        let mut at = order;
        if use_d {
            a[r_re * ncols + at] = wt;
            at += 1;
        }
        if use_e {
            a[r_im * ncols + at] = wt * w;
        }
        rhs[r_re] = wt * hv.re;
        rhs[r_im] = wt * hv.im;
    }
    let x = qr(&a, nrows, ncols).solve_ls(&rhs);
    let coords = x[..order].to_vec();
    let mut at = order;
    let d = if use_d {
        let v = x[at];
        at += 1;
        v
    } else {
        0.0
    };
    let e = if use_e { x[at] } else { 0.0 };
    (coords, d, e)
}

fn fit_errors(omega: &[f64], h: &[C64], weights: &[f64], model: &RationalModel) -> (f64, f64) {
    let mut sq = 0.0f64;
    let mut maxabs = 0.0f64;
    for ((&w, &hv), &wt) in omega.iter().zip(h).zip(weights) {
        let err = (model.eval_iw(w) - hv).abs();
        maxabs = maxabs.max(err);
        sq += (wt * err) * (wt * err);
    }
    (fs_math::det::sqrt(sq / (2.0 * omega.len() as f64)), maxabs)
}

/// Canonical conjugate-closed terms from an expanded pole list (public
/// hook shared with the Loewner front end).
#[must_use]
pub fn terms_from_poles(poles: &[C64]) -> Vec<PoleTerm> {
    conjugate_close(poles)
}

/// Residue-only fit at FIXED poles (the shared final pass both front
/// ends use), with the same `d >= 0` / `e >= 0` active-set handling.
///
/// # Errors
/// Typed [`VfError`] on degenerate input.
pub fn residue_fit_at_poles(
    omega: &[f64],
    h: &[C64],
    terms: &[PoleTerm],
    opts: &FitOptions,
) -> Result<FitOutcome, VfError> {
    for (i, (&w, v)) in omega.iter().zip(h).enumerate() {
        if !w.is_finite() || w <= 0.0 || !v.re.is_finite() || !v.im.is_finite() {
            return Err(VfError::BadSample { index: i });
        }
    }
    let order: usize = terms.iter().map(PoleTerm::state_dim).sum();
    if order == 0 {
        return Err(VfError::BadOrder { order });
    }
    if 2 * omega.len() < order + 2 || omega.len() != h.len() {
        return Err(VfError::TooFewSamples {
            samples: omega.len(),
            needed: (order + 2).div_ceil(2),
        });
    }
    let weights: Vec<f64> = h.iter().map(|&v| opts.weights.weight(v)).collect();
    let (model, weighted_rms, max_abs) = residue_pass(omega, h, &weights, terms, opts)?;
    Ok(FitOutcome {
        model,
        report: FitReport {
            weighted_rms,
            max_abs_error: max_abs,
            weights: opts.weights.label(),
            iterations_run: 0,
            final_pole_movement: 0.0,
        },
    })
}

/// Ascending-order automatic model selection: fit at each order in
/// `orders`, stop when the weighted error stops improving by at least
/// `plateau_ratio` (or dips below `noise_floor`, the overfit refusal),
/// and return the selected outcome plus the full order-vs-error curve.
///
/// # Errors
/// First fit error encountered; at least one order must succeed.
pub fn fit_auto_order(
    omega: &[f64],
    h: &[C64],
    orders: &[usize],
    base: &FitOptions,
    plateau_ratio: f64,
    noise_floor: f64,
) -> Result<(FitOutcome, Vec<(usize, f64)>), VfError> {
    assert!(!orders.is_empty(), "orders list must be non-empty");
    let mut curve: Vec<(usize, f64)> = Vec::new();
    let mut best: Option<FitOutcome> = None;
    for &order in orders {
        let opts = FitOptions { order, ..*base };
        let outcome = vector_fit(omega, h, &opts)?;
        let err = outcome.report.weighted_rms;
        curve.push((order, err));
        let improved = match &best {
            None => true,
            Some(b) => err < b.report.weighted_rms * (1.0 - plateau_ratio),
        };
        let hit_floor = err <= noise_floor;
        if improved {
            best = Some(outcome);
        }
        if hit_floor || !improved {
            break;
        }
    }
    Ok((best.expect("at least one order fitted"), curve))
}
