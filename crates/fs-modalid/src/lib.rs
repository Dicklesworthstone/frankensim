//! # fs-modalid — experimental modal identification
//!
//! Fits measured frequency-response functions (impact-hammer / laser
//! FRFs of real structures) into modal parameters `(f_k, eta_k,
//! phi_k)` — the measurement side of the simulate-vs-measure loop.
//!
//! Three INDEPENDENT identifiers over one shared pole/residue core
//! (fs-vfit's `RationalModel` — explicitly the same code, tested
//! once): relaxed vector fitting, the Loewner matrix pencil (both
//! re-exported from fs-vfit), and the classical rational-fraction-
//! polynomial method on a Forsythe-orthogonal basis ([`rfp_fit`]).
//! Disagreement between identifiers is a DIAGNOSTIC, not something to
//! average away.
//!
//! Mode-quality machinery: stabilization-diagram automation with
//! machine-readable accept/reject per pole ([`stabilization`]), MAC
//! matrices and pairing ([`mac`]), exponential-window damping-bias
//! correction ([`correct_exponential_window`]), split-sample
//! confidence intervals, and an SNR gate that REFUSES to identify
//! below the noise floor (named refusal, never a fabricated mode).
//!
//! Convention: measured FRF data enters in the engineering
//! `e^{-i*omega*t}` convention typical of acoustics instrumentation;
//! [`FrfData::new`] CONJUGATES it onto the Laplace axis fs-vfit uses
//! (the executed-lesson doctrine from the vector-fitting bead). Data
//! already on the Laplace axis uses [`FrfData::new_laplace`].

use fs_math::c64::C64;
use fs_vfit::vf::{FitOptions, residue_fit_at_poles, terms_from_poles};
use fs_vfit::{RationalModel, vector_fit};

/// One channel of measured FRF samples (shared frequency grid).
#[derive(Debug, Clone)]
pub struct FrfChannel {
    /// Complex FRF samples on the Laplace axis (already conjugated).
    pub h: Vec<C64>,
    /// Optional per-sample coherence in `[0, 1]` (feeds the SNR
    /// estimate when present).
    pub coherence: Option<Vec<f64>>,
}

/// A measured multi-channel FRF data set on a shared frequency grid.
#[derive(Debug, Clone)]
pub struct FrfData {
    /// Angular frequencies [rad/s], ascending.
    pub omega: Vec<f64>,
    /// Channels (e.g. response points of a roving-hammer test).
    pub channels: Vec<FrfChannel>,
}

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum ModalIdError {
    /// Grid/channel shape mismatch or empty data.
    Shape {
        /// What disagreed.
        what: &'static str,
    },
    /// Estimated SNR below the caller's floor: identification REFUSES
    /// rather than fabricating modes.
    SnrTooLow {
        /// Estimated SNR (linear power ratio).
        snr: f64,
        /// The floor that was not met.
        floor: f64,
    },
    /// The underlying fit failed.
    Fit(String),
    /// CSV parse failure at a line.
    Csv {
        /// 1-based line number.
        line: usize,
    },
    /// The stabilization diagram accepted no poles — the data carries
    /// no identifiable stable modes under the rule.
    NoStablePoles,
}

impl core::fmt::Display for ModalIdError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ModalIdError::Shape { what } => write!(f, "shape mismatch: {what}"),
            ModalIdError::SnrTooLow { snr, floor } => {
                write!(
                    f,
                    "SNR {snr:.2e} below floor {floor:.2e}: refusing to identify"
                )
            }
            ModalIdError::Fit(e) => write!(f, "identification fit failed: {e}"),
            ModalIdError::Csv { line } => write!(f, "CSV parse failure at line {line}"),
            ModalIdError::NoStablePoles => {
                write!(f, "no poles survived the stabilization rule")
            }
        }
    }
}

impl std::error::Error for ModalIdError {}

impl FrfData {
    /// Build from engineering-convention (`e^{-i*omega*t}`) samples:
    /// CONJUGATES onto the Laplace axis (the load-bearing convention
    /// bridge from the vector-fitting bead).
    ///
    /// # Errors
    /// [`ModalIdError::Shape`] on mismatched lengths or empty data.
    pub fn new(
        omega: Vec<f64>,
        channels_engineering: Vec<(Vec<C64>, Option<Vec<f64>>)>,
    ) -> Result<Self, ModalIdError> {
        let conv: Vec<(Vec<C64>, Option<Vec<f64>>)> = channels_engineering
            .into_iter()
            .map(|(h, c)| (h.iter().map(|v| v.conj()).collect(), c))
            .collect();
        Self::new_laplace(omega, conv)
    }

    /// Build from samples already on the Laplace axis.
    ///
    /// # Errors
    /// [`ModalIdError::Shape`].
    pub fn new_laplace(
        omega: Vec<f64>,
        channels: Vec<(Vec<C64>, Option<Vec<f64>>)>,
    ) -> Result<Self, ModalIdError> {
        if omega.is_empty() || channels.is_empty() {
            return Err(ModalIdError::Shape { what: "empty data" });
        }
        let n = omega.len();
        let mut out = Vec::with_capacity(channels.len());
        for (h, coherence) in channels {
            if h.len() != n {
                return Err(ModalIdError::Shape {
                    what: "channel length vs grid",
                });
            }
            if coherence.as_ref().is_some_and(|c| c.len() != n) {
                return Err(ModalIdError::Shape {
                    what: "coherence length vs grid",
                });
            }
            out.push(FrfChannel { h, coherence });
        }
        Ok(FrfData {
            omega,
            channels: out,
        })
    }

    /// Parse the v1 CSV FRF table: header-free lines
    /// `omega,re_1,im_1[,coh_1][,re_2,im_2[,coh_2],...]` with a fixed
    /// channel count and `with_coherence` flag. Engineering
    /// convention (conjugated on ingest).
    ///
    /// # Errors
    /// [`ModalIdError::Csv`] with the offending line.
    pub fn parse_csv(
        text: &str,
        n_channels: usize,
        with_coherence: bool,
    ) -> Result<Self, ModalIdError> {
        let per = if with_coherence { 3 } else { 2 };
        let mut omega = Vec::new();
        let mut chans: Vec<(Vec<C64>, Vec<f64>)> = vec![(Vec::new(), Vec::new()); n_channels];
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split(',').map(str::trim).collect();
            if fields.len() != 1 + per * n_channels {
                return Err(ModalIdError::Csv { line: i + 1 });
            }
            let parse = |s: &str| {
                s.parse::<f64>()
                    .map_err(|_| ModalIdError::Csv { line: i + 1 })
            };
            omega.push(parse(fields[0])?);
            for (c, chan) in chans.iter_mut().enumerate() {
                let base = 1 + per * c;
                let re = parse(fields[base])?;
                let im = parse(fields[base + 1])?;
                chan.0.push(C64::new(re, im));
                if with_coherence {
                    chan.1.push(parse(fields[base + 2])?);
                }
            }
        }
        let channels = chans
            .into_iter()
            .map(|(h, coh)| (h, if with_coherence { Some(coh) } else { None }))
            .collect();
        Self::new(omega, channels)
    }
}

/// Per-channel SNR estimate (median over channels): with coherence,
/// `gamma^2/(1-gamma^2)` averaged; without, a SECOND-DIFFERENCE
/// roughness estimator — white noise passes through `h[i+1] - 2 h[i]
/// + h[i-1]` amplified by 6 in power while a smooth FRF contributes
/// only curvature, so `median |d2|^2 / 6` estimates the noise power
/// and `median |h|^2` the signal power.
///
/// The naive median-over-lowest-decile magnitude ratio is WRONG for
/// FRFs: antiresonance valleys are genuine signal, not noise
/// (executed failure: clean data read as SNR 7).
#[must_use]
pub fn estimate_snr(data: &FrfData) -> f64 {
    let mut snrs = Vec::new();
    for ch in &data.channels {
        if let Some(coh) = &ch.coherence {
            let mut acc = 0.0;
            let mut n = 0usize;
            for &g in coh {
                let g = g.clamp(0.0, 1.0 - 1.0e-12);
                acc += g / (1.0 - g);
                n += 1;
            }
            snrs.push(acc / n as f64);
        } else {
            let n = ch.h.len();
            if n < 3 {
                snrs.push(f64::INFINITY);
                continue;
            }
            let mut d2: Vec<f64> = (1..n - 1)
                .map(|i| {
                    let v = ch.h[i + 1] - ch.h[i].scale(2.0) + ch.h[i - 1];
                    v.norm_sq()
                })
                .collect();
            d2.sort_by(f64::total_cmp);
            let noise = (d2[d2.len() / 2] / 6.0).max(f64::MIN_POSITIVE);
            let mut mags: Vec<f64> = ch.h.iter().map(|v| v.norm_sq()).collect();
            mags.sort_by(f64::total_cmp);
            snrs.push(mags[mags.len() / 2] / noise);
        }
    }
    snrs.sort_by(f64::total_cmp);
    snrs[snrs.len() / 2]
}

// ---------------------------------------------------------------------
// RFP: rational fraction polynomial on a Forsythe-orthogonal basis
// ---------------------------------------------------------------------

/// Forsythe-orthogonal polynomial basis evaluated on the (scaled)
/// imaginary axis: three-term recurrence orthonormal under the
/// sample-weighted inner product. Returns `deg + 1` rows of length
/// `points.len()` plus the recurrence coefficients needed to convert
/// back to monomials.
fn forsythe_basis(points: &[f64], weights: &[f64], deg: usize) -> (Vec<Vec<C64>>, Vec<(f64, f64)>) {
    // Basis in the variable s = i*w_scaled; polynomials of i*w with
    // REAL recurrence coefficients (even/odd structure preserved by
    // the imaginary axis).
    let n = points.len();
    let mut rows: Vec<Vec<C64>> = Vec::with_capacity(deg + 1);
    let mut coeffs: Vec<(f64, f64)> = Vec::with_capacity(deg + 1);
    let mut prev: Vec<C64> = vec![C64::ZERO; n];
    let mut cur: Vec<C64> = vec![C64::ONE; n];
    // Normalize p0.
    let norm0 = fs_math::det::sqrt(weights.iter().sum::<f64>());
    for v in &mut cur {
        *v = v.scale(1.0 / norm0);
    }
    rows.push(cur.clone());
    coeffs.push((norm0, 0.0));
    for _d in 0..deg {
        // next = s*cur - beta*prev; beta = <s*cur, prev> (real by the
        // even/odd structure); then normalize.
        let mut next: Vec<C64> = (0..n).map(|i| C64::new(0.0, points[i]) * cur[i]).collect();
        let mut beta = 0.0;
        for i in 0..n {
            // <s cur, prev> under sum w * conj(prev) * (s cur); the
            // real part carries the coefficient.
            let ip = prev[i].conj() * next[i];
            beta += weights[i] * ip.re;
        }
        for i in 0..n {
            next[i] = next[i] - prev[i].scale(beta);
        }
        let mut norm_sq = 0.0;
        for i in 0..n {
            norm_sq += weights[i] * next[i].norm_sq();
        }
        let norm = fs_math::det::sqrt(norm_sq).max(f64::MIN_POSITIVE);
        for v in &mut next {
            *v = v.scale(1.0 / norm);
        }
        prev = core::mem::take(&mut cur);
        cur = next;
        rows.push(cur.clone());
        coeffs.push((norm, beta));
    }
    (rows, coeffs)
}

/// Convert Forsythe coefficients to monomial coefficients (in the
/// scaled variable `s`): rebuild each basis polynomial's monomial
/// expansion through the recurrence and accumulate.
fn forsythe_to_monomial(coeffs: &[(f64, f64)], x: &[f64]) -> Vec<f64> {
    let deg = coeffs.len() - 1;
    // poly[d] = monomial coefficients (real, in powers of s where odd
    // powers carry i — but the recurrence s*p keeps REAL monomial
    // coefficient vectors in the variable s directly).
    let mut polys: Vec<Vec<f64>> = Vec::with_capacity(deg + 1);
    let p0 = vec![1.0 / coeffs[0].0];
    polys.push(p0);
    for d in 1..=deg {
        let (norm, beta) = coeffs[d];
        let mut next = vec![0.0; d + 1];
        // s * polys[d-1]
        for (k, &c) in polys[d - 1].iter().enumerate() {
            next[k + 1] += c;
        }
        // - beta * polys[d-2]
        if d >= 2 {
            for (k, &c) in polys[d - 2].iter().enumerate() {
                next[k] -= beta * c;
            }
        } else {
            // d == 1: prev is the zero polynomial; nothing to subtract.
        }
        for v in &mut next {
            *v /= norm;
        }
        polys.push(next);
    }
    let mut out = vec![0.0; deg + 1];
    for (d, poly) in polys.iter().enumerate() {
        for (k, &c) in poly.iter().enumerate() {
            out[k] += x[d] * c;
        }
    }
    out
}

/// Classical rational-fraction-polynomial identification on the
/// Forsythe-orthogonal basis: common-denominator LS over one channel,
/// denominator roots (in the SCALED variable) via the companion
/// matrix, stable-flipped, then residues through the SHARED fs-vfit
/// residue pass. The orthogonal basis is what keeps the normal
/// equations well-conditioned at useful orders (the monomial-basis
/// comparison is pinned in the battery).
///
/// # Errors
/// [`ModalIdError::Fit`] on degenerate data or eigensolve failure.
pub fn rfp_fit(
    omega: &[f64],
    h: &[C64],
    order: usize,
    opts: &FitOptions,
) -> Result<fs_vfit::FitOutcome, ModalIdError> {
    if omega.len() < 2 * (order + 1) {
        return Err(ModalIdError::Fit("too few samples for RFP".to_string()));
    }
    // Scale the axis to O(1) for the polynomial basis.
    let w_max = omega.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    if w_max <= 0.0 {
        return Err(ModalIdError::Fit("degenerate frequency grid".to_string()));
    }
    let ws: Vec<f64> = omega.iter().map(|&w| w / w_max).collect();
    let weights: Vec<f64> = h
        .iter()
        .map(|v| {
            let m = v.abs();
            if m > 0.0 { 1.0 / m } else { 1.0 }
        })
        .collect();
    // Numerator degree = order, denominator degree = order (monic-ish
    // via the normalization row below).
    let (num_basis, _) = forsythe_basis(&ws, &weights, order);
    let (den_basis, den_coeffs) = forsythe_basis(&ws, &weights, order);
    let n = ws.len();
    let n_num = order + 1;
    let n_den = order + 1;
    let ncols = n_num + n_den;
    // Rows: Re/Im of num(s_i) - H_i * den(s_i) = 0, plus the
    // non-triviality row sum Re(den) = n (the relaxed normalization,
    // same trick as vector fitting).
    let nrows = 2 * n + 1;
    let mut a = vec![0.0f64; nrows * ncols];
    let mut rhs = vec![0.0f64; nrows];
    for i in 0..n {
        let wt = weights[i];
        for k in 0..n_num {
            let v = num_basis[k][i];
            a[(2 * i) * ncols + k] = wt * v.re;
            a[(2 * i + 1) * ncols + k] = wt * v.im;
        }
        for k in 0..n_den {
            let v = -(h[i] * den_basis[k][i]);
            a[(2 * i) * ncols + n_num + k] = wt * v.re;
            a[(2 * i + 1) * ncols + n_num + k] = wt * v.im;
        }
    }
    let mean_wt = weights.iter().sum::<f64>() / n as f64;
    for k in 0..n_den {
        let acc: f64 = den_basis[k].iter().map(|v| v.re).sum();
        a[(nrows - 1) * ncols + n_num + k] = mean_wt * acc;
    }
    rhs[nrows - 1] = mean_wt * n as f64;
    let x = fs_la::factor::qr(&a, nrows, ncols).solve_ls(&rhs);
    // Denominator monomial coefficients in the scaled variable.
    let den_mono = forsythe_to_monomial(&den_coeffs, &x[n_num..]);
    // Companion-matrix roots of den(s) (degree = order).
    let lead = den_mono[order];
    if lead.abs() < 1.0e-14 * den_mono.iter().fold(0.0f64, |acc, &v| acc.max(v.abs())) {
        return Err(ModalIdError::Fit(
            "denominator degree collapsed".to_string(),
        ));
    }
    let mut comp = vec![C64::ZERO; order * order];
    for i in 1..order {
        comp[i * order + (i - 1)] = C64::ONE;
    }
    for i in 0..order {
        comp[i * order + (order - 1)] = C64::from_re(-den_mono[i] / lead);
    }
    let mut roots =
        fs_la::eigen_complex::eig(&comp, order).map_err(|e| ModalIdError::Fit(e.to_string()))?;
    // Unscale, stable-flip.
    for r in &mut roots {
        *r = r.scale(w_max);
        if r.re > 0.0 {
            r.re = -r.re;
        }
        if r.re == 0.0 {
            r.re = -1.0e-6 * (1.0 + r.im.abs());
        }
    }
    let terms = terms_from_poles(&roots);
    residue_fit_at_poles(omega, h, &terms, opts).map_err(|e| ModalIdError::Fit(e.to_string()))
}

// ---------------------------------------------------------------------
// Modal parameters, stabilization, MAC
// ---------------------------------------------------------------------

/// One identified mode.
#[derive(Debug, Clone)]
pub struct Mode {
    /// Natural frequency [Hz].
    pub frequency_hz: f64,
    /// Viscous damping ratio `zeta` (loss factor `eta = 2 zeta` for
    /// light damping).
    pub damping_ratio: f64,
    /// Complex residue per channel — the (unscaled) mode shape.
    pub shape: Vec<C64>,
    /// Split-sample half-width confidence interval on the frequency
    /// [Hz] (see [`identify`]).
    pub frequency_ci_hz: f64,
    /// Split-sample half-width confidence interval on `zeta`.
    pub damping_ci: f64,
}

/// Machine-readable verdict for one stabilization-diagram pole.
#[derive(Debug, Clone)]
pub struct PoleVerdict {
    /// Frequency [Hz].
    pub frequency_hz: f64,
    /// Damping ratio.
    pub damping_ratio: f64,
    /// Number of consecutive orders the pole was stable over.
    pub stable_orders: usize,
    /// Accepted into the mode table?
    pub accepted: bool,
    /// Why (accepted or rejected) — machine-readable code.
    pub reason: &'static str,
}

/// Stabilization tolerances: a pole at order `k+1` matches one at
/// order `k` if frequency agrees within `freq_tol` (relative) and
/// damping within `damp_tol` (relative).
#[derive(Debug, Clone, Copy)]
pub struct StabilizationRule {
    /// Relative frequency tolerance (classic 0.5% default).
    pub freq_tol: f64,
    /// Relative damping tolerance (classic 5% default).
    pub damp_tol: f64,
    /// Consecutive stable orders required for acceptance.
    pub required_runs: usize,
}

impl Default for StabilizationRule {
    fn default() -> Self {
        StabilizationRule {
            freq_tol: 0.005,
            damp_tol: 0.05,
            required_runs: 3,
        }
    }
}

/// Pole tables per model order -> stabilization verdicts. Input: for
/// each order (ascending), the identified poles as `(omega_n, zeta)`.
#[must_use]
pub fn stabilization(orders: &[Vec<(f64, f64)>], rule: &StabilizationRule) -> Vec<PoleVerdict> {
    // Track runs: for each pole in the FINAL order, count how many
    // consecutive previous orders contain a matching pole.
    let Some(last) = orders.last() else {
        return Vec::new();
    };
    let mut verdicts = Vec::new();
    for &(wn, zeta) in last {
        let mut runs = 0usize;
        for prev in orders.iter().rev().skip(1) {
            let matched = prev.iter().any(|&(w2, z2)| {
                (w2 - wn).abs() <= rule.freq_tol * wn
                    && (z2 - zeta).abs() <= rule.damp_tol * zeta.max(1.0e-6)
            });
            if matched {
                runs += 1;
            } else {
                break;
            }
        }
        let accepted = runs + 1 >= rule.required_runs;
        verdicts.push(PoleVerdict {
            frequency_hz: wn / (2.0 * core::f64::consts::PI),
            damping_ratio: zeta,
            stable_orders: runs + 1,
            accepted,
            reason: if accepted {
                "stable-run"
            } else {
                "unstable-across-orders"
            },
        });
    }
    verdicts
}

/// Modal assurance criterion between two complex shapes:
/// `|a^H b|^2 / (a^H a)(b^H b)`.
#[must_use]
pub fn mac(a: &[C64], b: &[C64]) -> f64 {
    let mut num = C64::ZERO;
    let mut na = 0.0;
    let mut nb = 0.0;
    for (x, y) in a.iter().zip(b) {
        num = num + x.conj() * *y;
        na += x.norm_sq();
        nb += y.norm_sq();
    }
    num.norm_sq() / (na * nb).max(f64::MIN_POSITIVE)
}

/// Full MAC matrix (rows = `a` modes, cols = `b` modes).
#[must_use]
pub fn mac_matrix(a: &[Vec<C64>], b: &[Vec<C64>]) -> Vec<Vec<f64>> {
    a.iter()
        .map(|sa| b.iter().map(|sb| mac(sa, sb)).collect())
        .collect()
}

/// Greedy MAC pairing: each `a` mode takes its best `b` match above
/// `floor`, without reuse; unmatched modes are reported with index
/// `None` — the calibration diff primitive.
#[must_use]
pub fn mac_pairing(a: &[Vec<C64>], b: &[Vec<C64>], floor: f64) -> Vec<(usize, Option<usize>, f64)> {
    let m = mac_matrix(a, b);
    let mut used = vec![false; b.len()];
    let mut out = Vec::with_capacity(a.len());
    for (i, row) in m.iter().enumerate() {
        let mut best: Option<(usize, f64)> = None;
        for (j, &v) in row.iter().enumerate() {
            if !used[j] && best.is_none_or(|(_, bv)| v > bv) {
                best = Some((j, v));
            }
        }
        match best {
            Some((j, v)) if v >= floor => {
                used[j] = true;
                out.push((i, Some(j), v));
            }
            Some((_, v)) => out.push((i, None, v)),
            None => out.push((i, None, 0.0)),
        }
    }
    out
}

/// Exponential-window damping-bias correction: an impact-hammer
/// window `e^{-t/tau}` ADDS `1/tau` to every modal decay rate, so
/// `zeta_measured = zeta_true + 1/(tau * omega_n)`. Returns the
/// corrected ratio and the (logged) delta; refuses nothing — a
/// negative corrected value is CLAMPED to zero and flagged by the
/// returned delta exceeding the measured value.
#[must_use]
pub fn correct_exponential_window(zeta_measured: f64, omega_n: f64, tau: f64) -> (f64, f64) {
    let delta = 1.0 / (tau * omega_n);
    ((zeta_measured - delta).max(0.0), delta)
}

/// Identification options.
#[derive(Debug, Clone, Copy)]
pub struct IdentifyOptions {
    /// Model orders for the stabilization ladder (ascending).
    pub min_order: usize,
    /// Final (largest) order — also the identification order.
    pub max_order: usize,
    /// Order step.
    pub order_step: usize,
    /// SNR floor (linear power ratio); below it identification
    /// REFUSES.
    pub snr_floor: f64,
    /// Stabilization rule.
    pub rule: StabilizationRule,
    /// Exponential-window time constant `tau` [s] (`None` = no window
    /// correction).
    pub window_tau: Option<f64>,
}

impl Default for IdentifyOptions {
    fn default() -> Self {
        IdentifyOptions {
            // Step 2 with headroom past the expected order: the
            // required_runs=3 acceptance needs >= 3 consecutive
            // orders that can HOLD every true mode (executed: a
            // step-4 ladder ending exactly at the modal count
            // accepted 1 of 10 modes).
            min_order: 4,
            max_order: 28,
            order_step: 2,
            snr_floor: 10.0,
            rule: StabilizationRule::default(),
            window_tau: None,
        }
    }
}

/// Identification outcome: the mode table plus the full
/// machine-readable evidence trail.
#[derive(Debug, Clone)]
pub struct Identification {
    /// Accepted modes, ascending in frequency.
    pub modes: Vec<Mode>,
    /// Stabilization verdicts for every final-order pole.
    pub verdicts: Vec<PoleVerdict>,
    /// Estimated SNR that passed the gate.
    pub snr: f64,
    /// Per-order pole tables `(omega_n, zeta)` (the stabilization
    /// diagram, machine-readable).
    pub diagram: Vec<Vec<(f64, f64)>>,
    /// Raw (pre-window-correction) damping ratios per accepted mode
    /// (equal to the corrected ones when no window was declared).
    pub damping_raw: Vec<f64>,
}

/// Poles of a fitted model as `(omega_n, zeta)` pairs (conjugate
/// pairs only — real poles are not vibration modes).
fn model_poles(model: &RationalModel) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for t in &model.terms {
        if let fs_vfit::PoleTerm::Pair { pole, .. } = t {
            let wn = pole.abs();
            let zeta = -pole.re / wn.max(f64::MIN_POSITIVE);
            out.push((wn, zeta));
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

/// Full identification pipeline on the REFERENCE channel (channel 0):
/// SNR gate, per-order vector fits (stabilization ladder), verdicts,
/// then per-channel residues at the accepted poles (shapes), window
/// correction, and split-sample confidence intervals (even/odd
/// half-grids re-identified at fixed poles; the CI is the half-spread
/// of the refit pole parameters).
///
/// # Errors
/// [`ModalIdError::SnrTooLow`] below the floor; fit failures.
#[allow(clippy::too_many_lines)] // one linear pipeline: gate -> ladder -> verdicts -> shapes -> CI
pub fn identify(data: &FrfData, opts: &IdentifyOptions) -> Result<Identification, ModalIdError> {
    let snr = estimate_snr(data);
    if snr < opts.snr_floor {
        return Err(ModalIdError::SnrTooLow {
            snr,
            floor: opts.snr_floor,
        });
    }
    let omega = &data.omega;
    let href = &data.channels[0].h;
    let fit_opts = FitOptions {
        fit_e: false,
        ..FitOptions::new(opts.max_order)
    };
    // Stabilization ladder.
    let mut diagram: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut order = opts.min_order;
    while order <= opts.max_order {
        let fo = FitOptions { order, ..fit_opts };
        let outcome = vector_fit(omega, href, &fo).map_err(|e| ModalIdError::Fit(e.to_string()))?;
        diagram.push(model_poles(&outcome.model));
        order += opts.order_step;
    }
    let verdicts = stabilization(&diagram, &opts.rule);
    // Accepted poles -> full pole set for the residue passes.
    let accepted: Vec<(f64, f64)> = verdicts
        .iter()
        .filter(|v| v.accepted)
        .map(|v| {
            (
                v.frequency_hz * 2.0 * core::f64::consts::PI,
                v.damping_ratio,
            )
        })
        .collect();
    if accepted.is_empty() {
        return Err(ModalIdError::NoStablePoles);
    }
    // Merge near-duplicate accepted poles (a spare pole parked on an
    // already-identified resonance stabilizes too — executed: 201.000
    // and 201.002 Hz both accepted). The merge tolerance is FIVE
    // TIMES TIGHTER than the cross-order match tolerance: the two
    // roles differ, and merging at the match tolerance swallowed a
    // genuine 0.5%-separated close pair (executed).
    let merge_tol = opts.rule.freq_tol / 5.0;
    let mut merged: Vec<(f64, f64)> = Vec::with_capacity(accepted.len());
    for &(wn, zeta) in &accepted {
        if merged
            .last()
            .is_none_or(|&(wp, _)| (wn - wp).abs() > merge_tol * wn)
        {
            merged.push((wn, zeta));
        }
    }
    let accepted = merged;
    // terms_from_poles expects the EXPANDED conjugate-complete list
    // (its count reconciliation folds pairs to real otherwise —
    // executed: upper-half-only input silently produced zero shapes).
    let poles: Vec<C64> = accepted
        .iter()
        .flat_map(|&(wn, zeta)| {
            let re = -zeta * wn;
            let im = wn * fs_math::det::sqrt((1.0 - zeta * zeta).max(0.0));
            [C64::new(re, im), C64::new(re, -im)]
        })
        .collect();
    let terms = terms_from_poles(&poles);
    // Shapes: per-channel residues at the SHARED accepted poles.
    let mut shapes: Vec<Vec<C64>> = vec![Vec::new(); accepted.len()];
    for ch in &data.channels {
        let fo = FitOptions {
            order: 2 * accepted.len(),
            ..fit_opts
        };
        let outcome = residue_fit_at_poles(omega, &ch.h, &terms, &fo)
            .map_err(|e| ModalIdError::Fit(e.to_string()))?;
        // Match model pair terms back to accepted poles by frequency.
        for (k, &(wn, _)) in accepted.iter().enumerate() {
            let mut best: Option<(f64, C64)> = None;
            for t in &outcome.model.terms {
                if let fs_vfit::PoleTerm::Pair { pole, residue } = t {
                    let d = (pole.abs() - wn).abs();
                    if best.is_none_or(|(bd, _)| d < bd) {
                        best = Some((d, *residue));
                    }
                }
            }
            shapes[k].push(best.map_or(C64::ZERO, |(_, r)| r));
        }
    }
    // Residue-significance gate: a stable spare pole that models the
    // noise floor carries residues orders below the physical modes
    // (executed: a band-edge spare survived stabilization with a
    // ~1e-5-relative shape and poisoned MAC pairing). Drop modes whose
    // largest channel residue is below 1e-4 of the global largest.
    let mut peak = f64::MIN_POSITIVE;
    for sh in &shapes {
        for r in sh {
            peak = peak.max(r.abs());
        }
    }
    let significant: Vec<bool> = shapes
        .iter()
        .map(|sh| sh.iter().any(|r| r.abs() >= 1.0e-4 * peak))
        .collect();
    let accepted: Vec<(f64, f64)> = accepted
        .iter()
        .zip(&significant)
        .filter_map(|(&p, &sig)| if sig { Some(p) } else { None })
        .collect();
    let shapes: Vec<Vec<C64>> = shapes
        .into_iter()
        .zip(&significant)
        .filter_map(|(sh, &sig)| if sig { Some(sh) } else { None })
        .collect();
    // Split-sample CI: refit poles on even/odd half-grids at the final
    // order and report the half-spread of the matched parameters.
    let half = |parity: usize| -> Result<Vec<(f64, f64)>, ModalIdError> {
        let om: Vec<f64> = omega
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 2 == parity)
            .map(|(_, &w)| w)
            .collect();
        let hh: Vec<C64> = href
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 2 == parity)
            .map(|(_, &v)| v)
            .collect();
        let outcome =
            vector_fit(&om, &hh, &fit_opts).map_err(|e| ModalIdError::Fit(e.to_string()))?;
        Ok(model_poles(&outcome.model))
    };
    let (even, odd) = (half(0)?, half(1)?);
    let nearest = |set: &[(f64, f64)], wn: f64| -> Option<(f64, f64)> {
        set.iter()
            .min_by(|a, b| (a.0 - wn).abs().total_cmp(&(b.0 - wn).abs()))
            .copied()
    };
    let mut modes = Vec::with_capacity(accepted.len());
    let mut damping_raw = Vec::with_capacity(accepted.len());
    for (k, &(wn, zeta)) in accepted.iter().enumerate() {
        let (fe, ze) = nearest(&even, wn).unwrap_or((wn, zeta));
        let (fo_, zo) = nearest(&odd, wn).unwrap_or((wn, zeta));
        let f_ci = 0.5 * (fe - fo_).abs() / (2.0 * core::f64::consts::PI);
        let z_ci = 0.5 * (ze - zo).abs();
        damping_raw.push(zeta);
        let (zeta_corr, _delta) = match opts.window_tau {
            Some(tau) => correct_exponential_window(zeta, wn, tau),
            None => (zeta, 0.0),
        };
        modes.push(Mode {
            frequency_hz: wn / (2.0 * core::f64::consts::PI),
            damping_ratio: zeta_corr,
            shape: shapes[k].clone(),
            frequency_ci_hz: f_ci,
            damping_ci: z_ci,
        });
    }
    Ok(Identification {
        modes,
        verdicts,
        snr,
        diagram,
        damping_raw,
    })
}
