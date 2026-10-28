//! Multiport impedance models with one pole set for every entry (IRP
//! M2.5-W3b).
//!
//! `Z(s) = SUM_k [R_k / (s - p_k) (+ conjugate term)] + D + s E` for an
//! `n x n` port matrix, on this crate's `s = i omega` axis. Conjugate
//! `e^{-i omega t}` acoustic data before fitting. Because one pole set
//! serves every entry, one state-space realization carries the whole matrix.
//!
//! - Fitting is relaxed vector fitting with a single sigma shared by every
//!   fitted entry. A reciprocal fit uses the upper triangle of the data's
//!   symmetric part. Its report still measures the misfit against the
//!   original data, asymmetry included.
//! - Each frequency is weighted by the inverse Frobenius norm of its data
//!   matrix, so the controlled quantity is the matrix relative error.
//! - Passivity is impedance-form: `Z(i w) + Z(i w)^H >= 0` at every real
//!   `w`. The Hamiltonian arm needs `D + D^T` positive definite; otherwise
//!   the certificate is the named weaker `CertificateClass::GridOnly`.
//! - Repair moves residues only; poles, `D` and `E` stay fixed, so
//!   stability cannot change. Each `(frequency, eigenvector)` pair gives a
//!   linear constraint. The constraints accumulate over rounds against the
//!   original model, and every round is re-certified.
use crate::model::{PoleTerm, RationalModel};
use crate::passivity::{CertificateClass, PassivityError};
use crate::vf::{
    FitOptions, FitReport, VfError, WeightPreset, basis_columns, conjugate_close, initial_poles,
    magnitude_scale, pole_movement, solve_ls_equilibrated, terms_with_residues,
};
use fs_la::eigen::jacobi_eigh;
use fs_la::eigen_complex::{EigFailure, eig, lu_complex};
use fs_la::factor::lu;
use fs_math::{c64::C64, det};

/// Shared bounded admission and original exact ZOH runtime for causal loads.
pub mod load;

/// A real multiport rational impedance with shared poles.
#[derive(Debug, Clone, PartialEq)]
pub struct MimoModel {
    /// Port count `n`.
    pub ports: usize,
    /// Shared pole terms. Their `residue` fields are zero; residues live in
    /// [`Self::residues`].
    pub poles: Vec<PoleTerm>,
    /// Per pole term, the row-major `n x n` residue matrix of the stored pole
    /// (real terms use only the real parts).
    pub residues: Vec<Vec<C64>>,
    /// Row-major direct term `D`.
    pub d: Vec<f64>,
    /// Row-major improper term `E` (`s E`).
    pub e: Vec<f64>,
    /// True when every fitted entry is mirrored (`Z = Z^T` by construction).
    pub reciprocal: bool,
}

/// Real state-space realization `x' = A x + B u`, `y = C x + D u + E u'`.
#[derive(Debug, Clone, PartialEq)]
pub struct MimoStateSpace {
    /// State dimension (model order times ports).
    pub n: usize,
    /// Port count.
    pub ports: usize,
    /// Row-major `n x n`.
    pub a: Vec<f64>,
    /// Row-major `n x ports`.
    pub b: Vec<f64>,
    /// Row-major `ports x n`.
    pub c: Vec<f64>,
    /// Row-major `ports x ports`.
    pub d: Vec<f64>,
    /// Row-major `ports x ports`.
    pub e: Vec<f64>,
}

/// Fit diagnostics in matrix-relative terms.
#[derive(Debug, Clone)]
pub struct MimoFitReport {
    /// Largest `||Z_fit - Z_data||_F / ||Z_data||_F` over the samples.
    pub max_relative_error: f64,
    /// Root-mean-square of the same per-sample ratio.
    pub rms_relative_error: f64,
    /// Largest `||Z_data - Z_data^T||_F / ||Z_data||_F` over the samples.
    pub data_reciprocity_defect: f64,
    /// Pole-relocation iterations actually run.
    pub iterations_run: usize,
    /// Largest pole movement in the last relocation [rad/s].
    pub final_pole_movement: f64,
}

impl MimoModel {
    /// Total pole count (pairs count 2).
    #[must_use]
    pub fn order(&self) -> usize {
        self.poles.iter().map(PoleTerm::state_dim).sum()
    }

    /// True iff every shared pole is in the open left half-plane.
    #[must_use]
    pub fn is_stable(&self) -> bool {
        self.poles.iter().all(PoleTerm::is_stable)
    }

    /// Row-major `Z(s)`.
    #[must_use]
    pub fn eval(&self, s: C64) -> Vec<C64> {
        let n = self.ports;
        let mut z: Vec<C64> = (0..n * n)
            .map(|i| C64::from_re(self.d[i]) + s.scale(self.e[i]))
            .collect();
        for (term, r) in self.poles.iter().zip(&self.residues) {
            match *term {
                PoleTerm::Real { pole, .. } => {
                    let phi = (s - C64::from_re(pole)).recip();
                    for (zi, ri) in z.iter_mut().zip(r) {
                        *zi = *zi + phi.scale(ri.re);
                    }
                }
                PoleTerm::Pair { pole, .. } => {
                    let (a, b) = ((s - pole).recip(), (s - pole.conj()).recip());
                    for (zi, ri) in z.iter_mut().zip(r) {
                        *zi = *zi + *ri * a + ri.conj() * b;
                    }
                }
            }
        }
        z
    }

    /// `Z(i omega)`.
    #[must_use]
    pub fn eval_iw(&self, omega: f64) -> Vec<C64> {
        self.eval(C64::new(0.0, omega))
    }

    /// Real block realization: each pole term contributes its scalar block
    /// (as `RationalModel::state_space`) times the `ports x ports` identity.
    #[must_use]
    pub fn state_space(&self) -> MimoStateSpace {
        let p = self.ports;
        let n = self.order() * p;
        let (mut a, mut b, mut c) = (vec![0.0; n * n], vec![0.0; n * p], vec![0.0; p * n]);
        let mut at = 0;
        for (term, r) in self.poles.iter().zip(&self.residues) {
            match *term {
                PoleTerm::Real { pole, .. } => {
                    for q in 0..p {
                        a[(at + q) * n + at + q] = pole;
                        b[(at + q) * p + q] = 1.0;
                        for row in 0..p {
                            c[row * n + at + q] = r[row * p + q].re;
                        }
                    }
                    at += p;
                }
                PoleTerm::Pair { pole, .. } => {
                    for q in 0..p {
                        let (i1, i2) = (at + q, at + p + q);
                        a[i1 * n + i1] = pole.re;
                        a[i1 * n + i2] = pole.im;
                        a[i2 * n + i1] = -pole.im;
                        a[i2 * n + i2] = pole.re;
                        b[i1 * p + q] = 1.0;
                        for row in 0..p {
                            c[row * n + i1] = 2.0 * r[row * p + q].re;
                            c[row * n + i2] = 2.0 * r[row * p + q].im;
                        }
                    }
                    at += 2 * p;
                }
            }
        }
        MimoStateSpace {
            n,
            ports: p,
            a,
            b,
            c,
            d: self.d.clone(),
            e: self.e.clone(),
        }
    }

    fn entries(&self) -> Vec<(usize, usize)> {
        entries(self.ports, self.reciprocal)
    }
}

impl MimoStateSpace {
    /// `C (sI - A)^{-1} B + D + s E` by dense complex solves: an independent
    /// route to [`MimoModel::eval`].
    ///
    /// # Errors
    /// The fs-la singularity signal when `s` is an eigenvalue of `A`.
    pub fn eval(&self, s: C64) -> Result<Vec<C64>, EigFailure> {
        let (n, p) = (self.n, self.ports);
        let mut m = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n {
                m[i * n + j] = C64::from_re(-self.a[i * n + j]);
            }
            m[i * n + i] = m[i * n + i] + s;
        }
        let factor = lu_complex(&m, n)?;
        let mut z: Vec<C64> = (0..p * p)
            .map(|i| C64::from_re(self.d[i]) + s.scale(self.e[i]))
            .collect();
        for q in 0..p {
            let mut x: Vec<C64> = (0..n).map(|i| C64::from_re(self.b[i * p + q])).collect();
            factor.solve(&mut x);
            for row in 0..p {
                let mut acc = C64::ZERO;
                for (i, xi) in x.iter().enumerate() {
                    acc = acc + xi.scale(self.c[row * n + i]);
                }
                z[row * p + q] = z[row * p + q] + acc;
            }
        }
        Ok(z)
    }
}

fn entries(ports: usize, reciprocal: bool) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for a in 0..ports {
        for b in 0..ports {
            if !reciprocal || a <= b {
                out.push((a, b));
            }
        }
    }
    out
}

fn frobenius(z: &[C64]) -> f64 {
    det::sqrt(z.iter().map(|v| v.norm_sq()).sum())
}

/// Joint common-pole fit of multiple FRF channels on one angular-frequency grid.
///
/// Each entry of responses is one channel, in caller order. A single relaxed
/// denominator is fitted to ALL channels simultaneously through the same
/// relocation kernel as [mimo_vector_fit]. A mode absent at one sensor can
/// therefore be identified from another; responses are never averaged together.
///
/// The weighting preset applies to the Euclidean norm of the channel vector
/// at each frequency. Channels must have compatible calibrated units; this is
/// not automatic sensor-noise whitening. A jointly zero sample uses unit weight.
/// The returned models have identical poles and channel-specific residues.
/// Direct and improper terms are unconstrained, as for the general MIMO fit;
/// no impedance passivity is inferred from a displacement or transfer FRF.
///
/// The report's weighted_rms is the weighted channel-vector error on normalized
/// data (relative vector error with InverseMagnitude); max_abs_error is the
/// largest physical channel-vector error norm. Coincident physical modes with
/// the same pole appear as a combined residue, not separately resolved shapes.
///
/// # Errors
/// [VfError] for invalid order, empty/ragged/nonfinite or entirely zero data,
/// unrepresentable normalization, sigma collapse, or eigensolve failure.
pub fn common_pole_vector_fit(
    omega: &[f64],
    responses: &[&[C64]],
    opts: &FitOptions,
) -> Result<(Vec<RationalModel>, FitReport), VfError> {
    fit_common_responses(omega, responses, opts, None)
}

/// Fit channel-specific residues at one declared shared pole set.
///
/// This uses the same vector weighting, normalization and unconstrained
/// residue owner as [common_pole_vector_fit], without relocating the poles.
/// Signed direct terms remain admissible for opposite-polarity transfer or
/// accelerance sensors; this is not an impedance passivity certificate.
///
/// # Errors
/// [VfError] for invalid poles/data or unrepresentable numerical results.
pub fn common_pole_residue_fit(
    omega: &[f64],
    responses: &[&[C64]],
    poles: &[PoleTerm],
    opts: &FitOptions,
) -> Result<(Vec<RationalModel>, FitReport), VfError> {
    for (index, pole) in poles.iter().enumerate() {
        let valid = match pole {
            PoleTerm::Real { pole, .. } => pole.is_finite(),
            PoleTerm::Pair { pole, .. } => {
                pole.re.is_finite() && pole.im.is_finite() && pole.im > 0.0
            }
        };
        if !valid {
            return Err(VfError::BadSample { index });
        }
    }
    let order = poles.iter().map(PoleTerm::state_dim).sum();
    let options = FitOptions { order, ..*opts };
    fit_common_responses(omega, responses, &options, Some(poles))
}

fn fit_common_responses(
    omega: &[f64],
    responses: &[&[C64]],
    opts: &FitOptions,
    fixed_poles: Option<&[PoleTerm]>,
) -> Result<(Vec<RationalModel>, FitReport), VfError> {
    if opts.order == 0 {
        return Err(VfError::BadOrder { order: opts.order });
    }
    if omega.len() < opts.order.saturating_add(2) {
        return Err(VfError::TooFewSamples {
            samples: omega.len(),
            needed: opts.order.saturating_add(2),
        });
    }
    if responses.is_empty() || responses.iter().any(|h| h.len() != omega.len()) {
        return Err(VfError::BadSample { index: 0 });
    }
    for (i, &w) in omega.iter().enumerate() {
        if !w.is_finite()
            || w <= 0.0
            || responses
                .iter()
                .any(|h| !h[i].re.is_finite() || !h[i].im.is_finite())
        {
            return Err(VfError::BadSample { index: i });
        }
    }
    let all: Vec<C64> = responses.iter().flat_map(|h| h.iter().copied()).collect();
    if !all.iter().any(|v| v.re != 0.0 || v.im != 0.0) {
        return Err(VfError::BadSample { index: 0 });
    }
    let scale = magnitude_scale(&all);
    if !scale.is_finite() || scale <= 0.0 {
        return Err(VfError::BadSample { index: 0 });
    }
    let h: Vec<Vec<C64>> = responses
        .iter()
        .map(|channel| {
            channel
                .iter()
                .map(|v| C64::new(v.re / scale, v.im / scale))
                .collect()
        })
        .collect();
    let mut weights = Vec::with_capacity(omega.len());
    for i in 0..omega.len() {
        let norm = h
            .iter()
            .fold(0.0, |norm, channel| det::hypot(norm, channel[i].abs()));
        if !norm.is_finite() {
            return Err(VfError::BadSample { index: i });
        }
        let weight = match opts.weights {
            WeightPreset::Uniform => 1.0,
            WeightPreset::InverseMagnitude if norm > 0.0 => 1.0 / norm,
            WeightPreset::LogBand if norm > 0.0 => 1.0 / det::sqrt(norm),
            _ => 1.0,
        };
        if !weight.is_finite() {
            return Err(VfError::BadSample { index: i });
        }
        weights.push(weight);
    }
    let per_entry = vec![*opts; responses.len()];
    let mut terms =
        fixed_poles.map_or_else(|| initial_poles(omega, opts.order), <[PoleTerm]>::to_vec);
    let (mut iterations_run, mut final_pole_movement) = (0, 0.0);
    let hi = omega.iter().copied().fold(0.0, f64::max);
    let iterations = if fixed_poles.is_some() {
        0
    } else {
        opts.iterations
    };
    for _ in 0..iterations {
        let (next, movement) = relocate_shared(omega, &h, &weights, &terms, &per_entry)?;
        terms = next;
        iterations_run += 1;
        final_pole_movement = movement;
        if movement < 1.0e-9 * hi {
            break;
        }
    }
    let normalized: Vec<_> = h
        .iter()
        .map(|channel| {
            let (coords, d, e) = residue_solve(omega, channel, &weights, &terms, opts);
            RationalModel {
                terms: terms_with_residues(&terms, &coords),
                d,
                e,
            }
        })
        .collect();
    let (mut weighted_sq, mut max_abs_error) = (0.0_f64, 0.0_f64);
    for (i, &w) in omega.iter().enumerate() {
        let error = normalized
            .iter()
            .zip(&h)
            .fold(0.0, |norm, (model, channel)| {
                det::hypot(norm, (model.eval_iw(w) - channel[i]).abs())
            });
        if !error.is_finite() || !(error * scale).is_finite() {
            return Err(VfError::BadSample { index: i });
        }
        weighted_sq += (weights[i] * error).powi(2);
        max_abs_error = max_abs_error.max(error * scale);
    }
    let models = normalized
        .into_iter()
        .map(|model| RationalModel {
            terms: model
                .terms
                .into_iter()
                .map(|term| match term {
                    PoleTerm::Real { pole, residue } => PoleTerm::Real {
                        pole,
                        residue: residue * scale,
                    },
                    PoleTerm::Pair { pole, residue } => PoleTerm::Pair {
                        pole,
                        residue: residue.scale(scale),
                    },
                })
                .collect(),
            d: model.d * scale,
            e: model.e * scale,
        })
        .collect();
    Ok((
        models,
        FitReport {
            weighted_rms: det::sqrt(weighted_sq / omega.len() as f64),
            max_abs_error,
            weights: opts.weights.label(),
            iterations_run,
            final_pole_movement,
        },
    ))
}

/// Shared-pole relaxed vector fit of `data[i]` (row-major `ports x ports`)
/// at `omega[i]`, on the `s = i omega` axis. `opts.weights` is ignored: rows
/// are weighted by the inverse Frobenius norm of each data matrix. `fit_d`
/// and `fit_e` apply to every entry, without sign constraints; the passivity
/// check owns the asymptotic conditions.
///
/// # Errors
/// [`VfError`] on empty, nonfinite or mis-shaped data, a zero order, too
/// few samples, sigma collapse or an eigensolve failure.
pub fn mimo_vector_fit(
    omega: &[f64],
    data: &[Vec<C64>],
    ports: usize,
    opts: &FitOptions,
    reciprocal: bool,
) -> Result<(MimoModel, MimoFitReport), VfError> {
    fit_shared(omega, data, ports, opts, reciprocal, &|_, _| true)
}

/// [`mimo_vector_fit`] with the improper `s E` term restricted to the
/// entries `improper(a, b)` admits.
fn fit_shared(
    omega: &[f64],
    data: &[Vec<C64>],
    ports: usize,
    opts: &FitOptions,
    reciprocal: bool,
    improper: &dyn Fn(usize, usize) -> bool,
) -> Result<(MimoModel, MimoFitReport), VfError> {
    if opts.order == 0 || ports == 0 {
        return Err(VfError::BadOrder { order: opts.order });
    }
    if omega.len() != data.len() || omega.len() < opts.order + 2 {
        return Err(VfError::TooFewSamples {
            samples: omega.len(),
            needed: opts.order + 2,
        });
    }
    for (i, (&w, z)) in omega.iter().zip(data).enumerate() {
        if !w.is_finite()
            || w <= 0.0
            || z.len() != ports * ports
            || z.iter().any(|v| !v.re.is_finite() || !v.im.is_finite())
            || !(frobenius(z) > 0.0)
        {
            return Err(VfError::BadSample { index: i });
        }
    }
    let list = entries(ports, reciprocal);
    let raw: Vec<Vec<C64>> = list
        .iter()
        .map(|&(a, b)| {
            data.iter()
                .map(|z| {
                    if reciprocal {
                        (z[a * ports + b] + z[b * ports + a]).scale(0.5)
                    } else {
                        z[a * ports + b]
                    }
                })
                .collect()
        })
        .collect();
    let all: Vec<C64> = raw.iter().flatten().copied().collect();
    let scale = magnitude_scale(&all);
    let h: Vec<Vec<C64>> = raw
        .iter()
        .map(|r| r.iter().map(|v| v.scale(1.0 / scale)).collect())
        .collect();
    let weights: Vec<f64> = data.iter().map(|z| scale / frobenius(z)).collect();
    let per_entry: Vec<FitOptions> = list
        .iter()
        .map(|&(a, b)| FitOptions {
            fit_e: opts.fit_e && improper(a, b),
            ..*opts
        })
        .collect();
    let mut terms = initial_poles(omega, opts.order);
    let (mut iterations_run, mut final_pole_movement) = (0, 0.0);
    let hi = omega.iter().copied().fold(0.0, f64::max);
    for _ in 0..opts.iterations {
        let (next, movement) = relocate_shared(omega, &h, &weights, &terms, &per_entry)?;
        terms = next;
        iterations_run += 1;
        final_pole_movement = movement;
        if movement < 1.0e-9 * hi {
            break;
        }
    }
    let n = ports;
    let mut residues = vec![vec![C64::ZERO; n * n]; terms.len()];
    let (mut d, mut e) = (vec![0.0; n * n], vec![0.0; n * n]);
    for ((&(a, b), hj), entry_opts) in list.iter().zip(&h).zip(&per_entry) {
        let (coords, dj, ej) = residue_solve(omega, hj, &weights, &terms, entry_opts);
        let fitted = terms_with_residues(&terms, &coords);
        let slots = if reciprocal && a != b {
            vec![a * n + b, b * n + a]
        } else {
            vec![a * n + b]
        };
        for slot in slots {
            for (k, t) in fitted.iter().enumerate() {
                residues[k][slot] = match *t {
                    PoleTerm::Real { residue, .. } => C64::from_re(residue * scale),
                    PoleTerm::Pair { residue, .. } => residue.scale(scale),
                };
            }
            d[slot] = dj * scale;
            e[slot] = ej * scale;
        }
    }
    let poles = terms
        .iter()
        .map(|t| match *t {
            PoleTerm::Real { pole, .. } => PoleTerm::Real { pole, residue: 0.0 },
            PoleTerm::Pair { pole, .. } => PoleTerm::Pair {
                pole,
                residue: C64::ZERO,
            },
        })
        .collect();
    let model = MimoModel {
        ports,
        poles,
        residues,
        d,
        e,
        reciprocal,
    };
    let (mut worst, mut sq, mut defect) = (0.0f64, 0.0, 0.0f64);
    for (&w, z) in omega.iter().zip(data) {
        let fit = model.eval_iw(w);
        let gap: Vec<C64> = fit.iter().zip(z).map(|(&f, &v)| f - v).collect();
        let rel = frobenius(&gap) / frobenius(z);
        worst = worst.max(rel);
        sq += rel * rel;
        let asym: Vec<C64> = (0..n * n).map(|i| z[i] - z[(i % n) * n + i / n]).collect();
        defect = defect.max(frobenius(&asym) / frobenius(z));
    }
    Ok((
        model,
        MimoFitReport {
            max_relative_error: worst,
            rms_relative_error: det::sqrt(sq / omega.len() as f64),
            data_reciprocity_defect: defect,
            iterations_run,
            final_pole_movement,
        },
    ))
}

/// One relocation with a sigma shared by every response: the stacked
/// version of the SISO relaxed system (same column layout per response,
/// with each response's own direct and improper columns).
fn relocate_shared(
    omega: &[f64],
    h: &[Vec<C64>],
    weights: &[f64],
    terms: &[PoleTerm],
    per_entry: &[FitOptions],
) -> Result<(Vec<PoleTerm>, f64), VfError> {
    let order: usize = terms.iter().map(PoleTerm::state_dim).sum();
    let widths: Vec<usize> = per_entry
        .iter()
        .map(|o| order + usize::from(o.fit_d) + usize::from(o.fit_e))
        .collect();
    let offsets: Vec<usize> = widths
        .iter()
        .scan(0, |at, &w| {
            *at += w;
            Some(*at - w)
        })
        .collect();
    let responses = h.len();
    let sigma0: usize = widths.iter().sum();
    let ncols = sigma0 + order + 1;
    let nrows = 2 * omega.len() * responses + 1;
    let mut a = vec![0.0f64; nrows * ncols];
    let mut rhs = vec![0.0f64; nrows];
    let mean_wt = weights.iter().sum::<f64>() / weights.len() as f64;
    for (i, (&w, &wt)) in omega.iter().zip(weights).enumerate() {
        let cols = basis_columns(terms, C64::new(0.0, w));
        for (j, hj) in h.iter().enumerate() {
            let hv = hj[i];
            let (r_re, r_im) = (2 * (j * omega.len() + i), 2 * (j * omega.len() + i) + 1);
            let (model0, opts) = (offsets[j], &per_entry[j]);
            for (k, &phi) in cols.iter().enumerate() {
                a[r_re * ncols + model0 + k] = wt * phi.re;
                a[r_im * ncols + model0 + k] = wt * phi.im;
                let m = -(hv * phi);
                a[r_re * ncols + sigma0 + k] = wt * m.re;
                a[r_im * ncols + sigma0 + k] = wt * m.im;
            }
            if opts.fit_d {
                a[r_re * ncols + model0 + order] = wt;
            }
            if opts.fit_e {
                a[r_im * ncols + model0 + order + usize::from(opts.fit_d)] = wt * w;
            }
            a[r_re * ncols + ncols - 1] = wt * -hv.re;
            a[r_im * ncols + ncols - 1] = wt * -hv.im;
        }
        let last = nrows - 1;
        for (k, &phi) in cols.iter().enumerate() {
            a[last * ncols + sigma0 + k] += mean_wt * phi.re;
        }
        a[last * ncols + ncols - 1] += mean_wt;
    }
    rhs[nrows - 1] = mean_wt * omega.len() as f64;
    let x = solve_ls_equilibrated(&mut a, nrows, ncols, &rhs);
    let d_tilde = x[ncols - 1];
    if d_tilde.abs() < 1.0e-12 {
        return Err(VfError::SigmaCollapse);
    }
    let sigma = crate::model::RationalModel {
        terms: terms_with_residues(terms, &x[sigma0..sigma0 + order]),
        d: d_tilde,
        e: 0.0,
    };
    let ss = sigma.state_space();
    let n = ss.n;
    let mut m = vec![C64::ZERO; n * n];
    for i in 0..n {
        for j in 0..n {
            m[i * n + j] = C64::from_re(ss.a[i * n + j] - ss.b[i] * ss.c[j] / d_tilde);
        }
    }
    let mut zeros = eig(&m, n).map_err(VfError::Eig)?;
    for z in &mut zeros {
        if z.re > 0.0 {
            z.re = -z.re;
        }
        if z.re == 0.0 {
            z.re = -1.0e-6 * (1.0 + z.im.abs());
        }
    }
    let next = conjugate_close(&zeros);
    let movement = pole_movement(terms, &next);
    Ok((next, movement))
}

/// Unconstrained residue pass for one response at fixed poles.
fn residue_solve(
    omega: &[f64],
    h: &[C64],
    weights: &[f64],
    terms: &[PoleTerm],
    opts: &FitOptions,
) -> (Vec<f64>, f64, f64) {
    let order: usize = terms.iter().map(PoleTerm::state_dim).sum();
    let ncols = order + usize::from(opts.fit_d) + usize::from(opts.fit_e);
    let nrows = 2 * omega.len();
    let mut a = vec![0.0f64; nrows * ncols];
    let mut rhs = vec![0.0f64; nrows];
    for (i, ((&w, &hv), &wt)) in omega.iter().zip(h).zip(weights).enumerate() {
        let cols = basis_columns(terms, C64::new(0.0, w));
        for (k, &phi) in cols.iter().enumerate() {
            a[2 * i * ncols + k] = wt * phi.re;
            a[(2 * i + 1) * ncols + k] = wt * phi.im;
        }
        if opts.fit_d {
            a[2 * i * ncols + order] = wt;
        }
        if opts.fit_e {
            a[(2 * i + 1) * ncols + order + usize::from(opts.fit_d)] = wt * w;
        }
        rhs[2 * i] = wt * hv.re;
        rhs[2 * i + 1] = wt * hv.im;
    }
    let x = solve_ls_equilibrated(&mut a, nrows, ncols, &rhs);
    let d = if opts.fit_d { x[order] } else { 0.0 };
    let e = if opts.fit_e {
        x[order + usize::from(opts.fit_d)]
    } else {
        0.0
    };
    (x[..order].to_vec(), d, e)
}

/// Outcome of a multiport passivity check.
#[derive(Debug, Clone)]
pub struct MimoPassivityReport {
    /// True iff no arm found a violation.
    pub passive: bool,
    /// Certificate strength.
    pub class: CertificateClass,
    /// Smallest grid eigenvalue of `Z + Z^H` and its frequency [rad/s].
    pub worst: (f64, f64),
    /// Hamiltonian imaginary-axis crossings [rad/s].
    pub crossings: Vec<f64>,
    /// Frequencies at which violations were sampled.
    pub violation_freqs: Vec<f64>,
}

/// Eigenvalues (ascending) and eigenvectors of the Hermitian part
/// `Z + Z^H` (not halved), via the real symmetric embedding. Each value of
/// the embedding appears twice; this returns one copy per pair.
fn hermitian_part(z: &[C64], n: usize) -> (Vec<f64>, Vec<Vec<C64>>) {
    let mut m = vec![0.0; 4 * n * n];
    for a in 0..n {
        for b in 0..n {
            let h = z[a * n + b] + z[b * n + a].conj();
            m[a * 2 * n + b] = h.re;
            m[(a + n) * 2 * n + b + n] = h.re;
            m[a * 2 * n + b + n] = -h.im;
            m[(a + n) * 2 * n + b] = h.im;
        }
    }
    let (values, vectors) = jacobi_eigh(&m, 2 * n);
    let (mut vals, mut vecs) = (Vec::with_capacity(n), Vec::with_capacity(n));
    let mut taken: Vec<Vec<C64>> = Vec::new();
    for (index, &value) in values.iter().enumerate() {
        if vals.len() == n {
            break;
        }
        let v: Vec<C64> = (0..n)
            .map(|i| C64::new(vectors[i * 2 * n + index], vectors[(i + n) * 2 * n + index]))
            .collect();
        // The partner [-y; x] is i*v: the same complex direction. Skip it.
        let duplicate = taken.iter().any(|u| {
            let overlap = u
                .iter()
                .zip(&v)
                .fold(C64::ZERO, |s, (&a, &b)| s + a.conj() * b);
            overlap.abs() > 0.5
        });
        if duplicate {
            continue;
        }
        let norm = frobenius(&v);
        let v: Vec<C64> = v.iter().map(|x| x.scale(1.0 / norm)).collect();
        taken.push(v.clone());
        vals.push(value);
        vecs.push(v);
    }
    (vals, vecs)
}

/// Smallest grid eigenvalue with its frequency, and, per contiguous negative
/// band, its first, smallest and last frequencies. A wide band would otherwise
/// hand the repair hundreds of near-duplicate constraints.
fn grid_scan(grid: &[f64], lambda_min: impl Fn(f64) -> f64) -> ((f64, f64), Vec<f64>) {
    let mut worst = (f64::INFINITY, 0.0);
    let mut violation_freqs = Vec::new();
    let mut run: Option<(f64, f64, f64)> = None; // (first, argmin, min)
    let mut last = 0.0;
    for &w in grid {
        let lam = lambda_min(w);
        if lam < worst.0 {
            worst = (lam, w);
        }
        if lam < 0.0 {
            run = Some(match run {
                None => (w, w, lam),
                Some((first, _, min)) if lam < min => (first, w, lam),
                Some(r) => r,
            });
        } else if let Some((first, arg, _)) = run.take() {
            violation_freqs.extend([first, arg, last]);
        }
        last = w;
    }
    if let Some((first, arg, _)) = run {
        violation_freqs.extend([first, arg, last]);
    }
    violation_freqs.dedup();
    (worst, violation_freqs)
}

fn build_grid(model: &MimoModel, band: (f64, f64)) -> Vec<f64> {
    let lo = (band.0 / 10.0).max(1.0e-3);
    let hi = band.1 * 10.0;
    let count = 2000;
    let span = det::ln(hi / lo);
    let mut grid: Vec<f64> = (0..count)
        .map(|i| lo * det::exp(span * i as f64 / (count - 1) as f64))
        .collect();
    for t in &model.poles {
        if let PoleTerm::Pair { pole, .. } = t {
            let w0 = pole.abs();
            grid.extend((0u32..32).map(|k| w0 * (0.9 + 0.2 * f64::from(k) / 31.0)));
        }
    }
    grid.sort_by(f64::total_cmp);
    grid.dedup();
    grid
}

/// Certify impedance-form passivity of a stable multiport model.
///
/// # Errors
/// [`PassivityError::Unstable`] for a pole on or right of the imaginary
/// axis; [`PassivityError::NegativeAsymptote`] when `E` is not symmetric
/// positive semidefinite or `D + D^T` has a negative eigenvalue;
/// [`PassivityError::Eig`] when the Hamiltonian eigensolve fails.
pub fn check_mimo_passivity(
    model: &MimoModel,
    band: (f64, f64),
) -> Result<MimoPassivityReport, PassivityError> {
    if !model.is_stable() {
        return Err(PassivityError::Unstable);
    }
    let n = model.ports;
    let scale = model
        .d
        .iter()
        .chain(&model.e)
        .fold(0.0f64, |m, v| m.max(v.abs()));
    let tol = 1.0e-12 * scale.max(f64::MIN_POSITIVE);
    let sym =
        |m: &[f64]| -> Vec<f64> { (0..n * n).map(|i| m[i] + m[(i % n) * n + i / n]).collect() };
    let e_skew = (0..n * n).any(|i| (model.e[i] - model.e[(i % n) * n + i / n]).abs() > tol);
    let (e_values, _) = jacobi_eigh(&sym(&model.e), n);
    let (r_values, r_vectors) = jacobi_eigh(&sym(&model.d), n);
    if e_skew || e_values[0] < -tol || r_values[0] < -tol {
        return Err(PassivityError::NegativeAsymptote);
    }
    let lambda_min = |w: f64| hermitian_part(&model.eval_iw(w), n).0[0];
    let (worst, mut violation_freqs) = grid_scan(&build_grid(model, band), lambda_min);
    let mut crossings = Vec::new();
    let class = if r_values[0] > 1.0e-9 * r_values[n - 1].abs().max(f64::MIN_POSITIVE) {
        crossings =
            hamiltonian_crossings(model, &r_values, &r_vectors).map_err(PassivityError::Eig)?;
        let mut probes: Vec<f64> = crossings
            .windows(2)
            .map(|p| det::sqrt(p[0] * p[1]))
            .collect();
        if let Some(&first) = crossings.first() {
            probes.push(0.5 * first);
        }
        if let Some(&final_crossing) = crossings.last() {
            probes.push(2.0 * final_crossing);
        }
        violation_freqs.extend(probes.into_iter().filter(|&w| lambda_min(w) < 0.0));
        CertificateClass::HamiltonianExact
    } else {
        CertificateClass::GridOnly
    };
    violation_freqs.sort_by(f64::total_cmp);
    violation_freqs.dedup();
    // As in the SISO check, an exact certificate also needs no crossing: an
    // eigenvalue touching zero between grid points is not a clean margin.
    let passive = violation_freqs.is_empty()
        && worst.0 >= 0.0
        && (class == CertificateClass::GridOnly || crossings.is_empty());
    Ok(MimoPassivityReport {
        passive,
        class,
        worst,
        crossings,
        violation_freqs,
    })
}

/// Crossings of an eigenvalue of `Z + Z^H` through zero, from the
/// eigenvalues of `M = [[F, -G], [Q, -F^T]]` on the imaginary axis, with
/// `R0 = D + D^T`, `F = A - B R0^-1 C`, `G = B R0^-1 B^T` and
/// `Q = C^T R0^-1 C` (the SISO matrix with `2d` replaced by `R0`).
fn hamiltonian_crossings(
    model: &MimoModel,
    r_values: &[f64],
    r_vectors: &[f64],
) -> Result<Vec<f64>, EigFailure> {
    let ss = model.state_space();
    let (n, p) = (ss.n, ss.ports);
    let mut r_inv = vec![0.0; p * p];
    for i in 0..p {
        for j in 0..p {
            r_inv[i * p + j] = (0..p)
                .map(|k| r_vectors[i * p + k] * r_vectors[j * p + k] / r_values[k])
                .sum();
        }
    }
    // B R0^-1 (n x p) and R0^-1 C (p x n).
    let mut b_r = vec![0.0; n * p];
    let mut r_c = vec![0.0; p * n];
    for i in 0..n {
        for j in 0..p {
            b_r[i * p + j] = (0..p).map(|k| ss.b[i * p + k] * r_inv[k * p + j]).sum();
        }
    }
    for i in 0..p {
        for j in 0..n {
            r_c[i * n + j] = (0..p).map(|k| r_inv[i * p + k] * ss.c[k * n + j]).sum();
        }
    }
    let dim = 2 * n;
    let mut m = vec![C64::ZERO; dim * dim];
    for i in 0..n {
        for j in 0..n {
            let mut f = ss.a[i * n + j];
            let (mut g, mut q) = (0.0, 0.0);
            for k in 0..p {
                f -= ss.b[i * p + k] * r_c[k * n + j];
                g += b_r[i * p + k] * ss.b[j * p + k];
                q += ss.c[k * n + i] * r_c[k * n + j];
            }
            m[i * dim + j] = C64::from_re(f);
            m[(n + j) * dim + n + i] = C64::from_re(-f);
            m[i * dim + n + j] = C64::from_re(-g);
            m[(n + i) * dim + j] = C64::from_re(q);
        }
    }
    let mut crossings: Vec<f64> = eig(&m, dim)?
        .into_iter()
        .filter(|lam| lam.re.abs() <= 1.0e-7 * (1.0 + lam.abs()) && lam.im > 0.0)
        .map(|lam| lam.im)
        .collect();
    crossings.sort_by(f64::total_cmp);
    crossings.dedup_by(|a, b| (*a - *b).abs() <= 1.0e-6 * (*a + *b));
    Ok(crossings)
}

/// Repair outcome.
#[derive(Debug, Clone)]
pub struct MimoRepairReport {
    /// Certification rounds used.
    pub rounds: usize,
    /// `||delta residue coords|| / ||residue coords||`.
    pub relative_perturbation: f64,
    /// Largest `||Z_repaired - Z_original||_F / ||Z_original||_F` over the
    /// requested band's log grid: the perturbation in operator terms.
    pub max_relative_change: f64,
    /// The post-repair certificate (passive on success).
    pub certificate: MimoPassivityReport,
}

fn coords(model: &MimoModel) -> Vec<f64> {
    let mut v = Vec::new();
    for (term, r) in model.poles.iter().zip(&model.residues) {
        for (a, b) in model.entries() {
            let z = r[a * model.ports + b];
            match term {
                PoleTerm::Real { .. } => v.push(z.re),
                PoleTerm::Pair { .. } => v.extend([z.re, z.im]),
            }
        }
    }
    v
}

fn with_coords(model: &MimoModel, x: &[f64]) -> MimoModel {
    let mut out = model.clone();
    let n = model.ports;
    let list = model.entries();
    let mut at = 0;
    for (term, r) in model.poles.iter().zip(&mut out.residues) {
        for &(a, b) in &list {
            let z = match term {
                PoleTerm::Real { .. } => {
                    at += 1;
                    C64::from_re(x[at - 1])
                }
                PoleTerm::Pair { .. } => {
                    at += 2;
                    C64::new(x[at - 2], x[at - 1])
                }
            };
            r[a * n + b] = z;
            if model.reciprocal {
                r[b * n + a] = z;
            }
        }
    }
    out
}

/// Gradient of `v^H (Z + Z^H) v` with respect to the residue coordinates at
/// `omega`. It is exact, because `Z` is linear in the residues.
fn rayleigh_gradient(model: &MimoModel, delays: Option<&[f64]>, omega: f64, v: &[C64]) -> Vec<f64> {
    let s = C64::new(0.0, omega);
    let n = model.ports;
    let list = model.entries();
    let lag = |a: usize, b: usize| {
        delays.map_or(C64::ONE, |t| {
            let (s, c) = det::sin_cos(omega * t[a * n + b]);
            C64::new(c, -s)
        })
    };
    let mut g = Vec::new();
    for term in &model.poles {
        let basis: Vec<C64> = match *term {
            PoleTerm::Real { pole, .. } => vec![(s - C64::from_re(pole)).recip()],
            PoleTerm::Pair { pole, .. } => {
                let (a, b) = ((s - pole).recip(), (s - pole.conj()).recip());
                let diff = a - b;
                vec![a + b, C64::new(-diff.im, diff.re)]
            }
        };
        for &(a, b) in &list {
            let mut weight = v[a].conj() * v[b] * lag(a, b);
            if model.reciprocal && a != b {
                weight = weight + v[b].conj() * v[a] * lag(b, a);
            }
            g.extend(basis.iter().map(|&phi| 2.0 * (weight * phi).re));
        }
    }
    g
}

/// Minimize `||delta||^2` subject to `rows[j] . delta >= rhs[j]` by an
/// active set of KKT equality solves (the SISO repair's QP shape).
fn qp_min_norm(rows: &[Vec<f64>], rhs: &[f64], nv: usize) -> Vec<f64> {
    let mut active: Vec<usize> = (0..rows.len()).filter(|&j| rhs[j] > 0.0).collect();
    let mut delta = vec![0.0; nv];
    for _ in 0..(4 * rows.len() + 8) {
        delta = vec![0.0; nv];
        if !active.is_empty() {
            let na = active.len();
            let mut gram = vec![0.0; na * na];
            for (i, &ri) in active.iter().enumerate() {
                for (j, &rj) in active.iter().enumerate() {
                    gram[i * na + j] = rows[ri].iter().zip(&rows[rj]).map(|(x, y)| x * y).sum();
                }
                gram[i * na + i] += 1.0e-14 * gram[i * na + i].abs().max(f64::MIN_POSITIVE);
            }
            let mut mu: Vec<f64> = active.iter().map(|&j| rhs[j]).collect();
            match lu(&gram, na) {
                Ok(factor) => factor.solve(&mut mu),
                // Degenerate active set (repeated constraints): drop the
                // newest and retry, as the SISO repair does.
                Err(_) => {
                    active.pop();
                    continue;
                }
            }
            if let Some(worst) = (0..na)
                .filter(|&i| mu[i] < 0.0)
                .min_by(|&a, &b| mu[a].total_cmp(&mu[b]))
            {
                active.remove(worst);
                continue;
            }
            for (i, &j) in active.iter().enumerate() {
                for (d, r) in delta.iter_mut().zip(&rows[j]) {
                    *d += mu[i] * r;
                }
            }
        }
        let violated = (0..rows.len())
            .filter(|j| !active.contains(j))
            .find(|&j| rows[j].iter().zip(&delta).map(|(x, y)| x * y).sum::<f64>() < rhs[j]);
        match violated {
            Some(j) => active.push(j),
            None => break,
        }
    }
    delta
}

/// Residue-perturbation passivity repair with poles, `D` and `E` fixed.
///
/// # Errors
/// Precondition errors as for [`check_mimo_passivity`], or
/// [`PassivityError::RepairFailed`] when 12 rounds do not certify.
pub fn repair_mimo_passivity(
    model: &MimoModel,
    band: (f64, f64),
) -> Result<(MimoModel, MimoRepairReport), PassivityError> {
    repair_with(model, None, band)
}

/// [`repair_mimo_passivity`] for a delayed model: the same residue QP with
/// every basis term carrying its entry's `exp(-i w tau)`, certified by the
/// grid-only [`check_delayed_mimo_passivity`].
///
/// # Errors
/// [`PassivityError::Unstable`] or [`PassivityError::RepairFailed`].
pub fn repair_delayed_mimo_passivity(
    delayed: &DelayedMimoModel,
    band: (f64, f64),
) -> Result<(DelayedMimoModel, MimoRepairReport), PassivityError> {
    let (model, report) = repair_with(&delayed.model, Some(&delayed.delays), band)?;
    Ok((
        DelayedMimoModel {
            delays: delayed.delays.clone(),
            model,
        },
        report,
    ))
}

fn repair_with(
    model: &MimoModel,
    delays: Option<&[f64]>,
    band: (f64, f64),
) -> Result<(MimoModel, MimoRepairReport), PassivityError> {
    const MAX_ROUNDS: usize = 12;
    let restore = |m: &MimoModel| DelayedMimoModel {
        delays: delays.map_or_else(|| vec![0.0; m.ports * m.ports], <[f64]>::to_vec),
        model: m.clone(),
    };
    let eval = |m: &MimoModel, w: f64| match delays {
        Some(_) => restore(m).eval_iw(w),
        None => m.eval_iw(w),
    };
    let check = |m: &MimoModel| match delays {
        Some(_) => check_delayed_mimo_passivity(&restore(m), band),
        None => check_mimo_passivity(m, band),
    };
    let base = coords(model);
    let base_norm = det::sqrt(base.iter().map(|v| v * v).sum());
    let mut current = model.clone();
    let (mut rows, mut values): (Vec<Vec<f64>>, Vec<f64>) = (Vec::new(), Vec::new());
    for round in 0..MAX_ROUNDS {
        let report = check(&current)?;
        if report.passive {
            let x = coords(&current);
            let delta = det::sqrt(x.iter().zip(&base).map(|(a, b)| (a - b) * (a - b)).sum());
            let span = det::ln(band.1 / band.0);
            let change = (0..200)
                .map(|i| {
                    let w = band.0 * det::exp(span * f64::from(i) / 199.0);
                    let (z0, z1) = (eval(model, w), eval(&current, w));
                    let gap: Vec<C64> = z1.iter().zip(&z0).map(|(&a, &b)| a - b).collect();
                    frobenius(&gap) / frobenius(&z0)
                })
                .fold(0.0, f64::max);
            return Ok((
                current,
                MimoRepairReport {
                    rounds: round,
                    relative_perturbation: delta / base_norm.max(f64::MIN_POSITIVE),
                    max_relative_change: change,
                    certificate: report,
                },
            ));
        }
        let mut freqs = report.violation_freqs.clone();
        freqs.extend(&report.crossings);
        if report.worst.0 < 0.0 {
            freqs.push(report.worst.1);
        }
        let floor = 1.0e-9 * frobenius(&eval(model, band.1));
        let dip = report.worst.0.min(0.0).abs();
        let margin = dip.max(floor) * 0.05 * (round as f64 + 1.0);
        for &w in &freqs {
            // Constrain every direction below this round's margin at w,
            // anchored on the original model (its Rayleigh quotient along v
            // is linear in the residues). A crossing frequency sits at zero.
            let (vals, vecs) = hermitian_part(&eval(&current, w), model.ports);
            for (lam, v) in vals.iter().zip(&vecs) {
                if *lam < margin {
                    let z0 = eval(model, w);
                    let n = model.ports;
                    let mut quotient = C64::ZERO;
                    for a in 0..n {
                        for b in 0..n {
                            quotient = quotient
                                + v[a].conj() * (z0[a * n + b] + z0[b * n + a].conj()) * v[b];
                        }
                    }
                    rows.push(rayleigh_gradient(model, delays, w, v));
                    values.push(quotient.re);
                }
            }
        }
        let rhs: Vec<f64> = values.iter().map(|&v| margin - v).collect();
        let delta = qp_min_norm(&rows, &rhs, base.len());
        let x: Vec<f64> = base.iter().zip(&delta).map(|(a, d)| a + d).collect();
        current = with_coords(model, &x);
    }
    Err(PassivityError::RepairFailed { rounds: MAX_ROUNDS })
}

/// A shared-pole model of the delay-extracted entries:
/// `Z_ab(s) = exp(-s tau_ab) Zhat_ab(s)`.
#[derive(Debug, Clone, PartialEq)]
pub struct DelayedMimoModel {
    /// Row-major per-entry propagation delays [s].
    pub delays: Vec<f64>,
    /// The rational model of the delay-free entries.
    pub model: MimoModel,
}

impl DelayedMimoModel {
    /// Row-major `Z(i omega)`, with the extracted delays restored.
    #[must_use]
    pub fn eval_iw(&self, omega: f64) -> Vec<C64> {
        self.model
            .eval_iw(omega)
            .into_iter()
            .zip(&self.delays)
            .map(|(z, &tau)| {
                let (s, c) = det::sin_cos(omega * tau);
                z * C64::new(c, -s)
            })
            .collect()
    }
}

/// Fit the delay-extracted entries `Z_ab(i w) exp(+i w tau_ab)` with shared
/// poles. The report measures the restored model against the original data.
/// Delays are the caller's declared propagation times, for example the
/// distance between aperture reference points over the sound speed.
///
/// # Errors
/// As for [`mimo_vector_fit`], plus `BadSample` for delays that are not
/// `ports x ports`, finite and nonnegative, or asymmetric in a reciprocal fit.
pub fn delayed_mimo_vector_fit(
    omega: &[f64],
    data: &[Vec<C64>],
    ports: usize,
    delays: &[f64],
    opts: &FitOptions,
    reciprocal: bool,
) -> Result<(DelayedMimoModel, MimoFitReport), VfError> {
    let n = ports;
    if delays.len() != n * n
        || delays.iter().any(|t| !t.is_finite() || *t < 0.0)
        || (reciprocal
            && (0..n * n).any(|i| delays[i].to_bits() != delays[(i % n) * n + i / n].to_bits()))
    {
        return Err(VfError::BadSample { index: 0 });
    }
    let shifted: Vec<Vec<C64>> = omega
        .iter()
        .zip(data)
        .map(|(&w, z)| {
            z.iter()
                .zip(delays)
                .map(|(&v, &tau)| {
                    let (s, c) = det::sin_cos(w * tau);
                    v * C64::new(c, s)
                })
                .collect()
        })
        .collect();
    // A delayed entry with an improper term would put `2 w E_ab sin(w tau)`
    // into the Hermitian part, which no residue can bound. Only undelayed
    // (self) entries carry `s E`.
    let (model, mut report) = fit_shared(omega, &shifted, ports, opts, reciprocal, &|a, b| {
        !(delays[a * n + b] > 0.0)
    })?;
    let delayed = DelayedMimoModel {
        delays: delays.to_vec(),
        model,
    };
    let (mut worst, mut sq) = (0.0f64, 0.0);
    for (&w, z) in omega.iter().zip(data) {
        let gap: Vec<C64> = delayed
            .eval_iw(w)
            .iter()
            .zip(z)
            .map(|(&f, &v)| f - v)
            .collect();
        let rel = frobenius(&gap) / frobenius(z);
        worst = worst.max(rel);
        sq += rel * rel;
    }
    report.max_relative_error = worst;
    report.rms_relative_error = det::sqrt(sq / omega.len() as f64);
    Ok((delayed, report))
}

/// Grid passivity of a delayed model: `Z + Z^H >= 0` on the dense log grid
/// of [`check_mimo_passivity`]. Delays put the model outside the rational
/// Hamiltonian test, so this is always the named `GridOnly` class.
///
/// # Errors
/// [`PassivityError::Unstable`] for an unstable delay-free part.
pub fn check_delayed_mimo_passivity(
    delayed: &DelayedMimoModel,
    band: (f64, f64),
) -> Result<MimoPassivityReport, PassivityError> {
    if !delayed.model.is_stable() {
        return Err(PassivityError::Unstable);
    }
    let n = delayed.model.ports;
    let (worst, violation_freqs) = grid_scan(&build_grid(&delayed.model, band), |w| {
        hermitian_part(&delayed.eval_iw(w), n).0[0]
    });
    Ok(MimoPassivityReport {
        passive: violation_freqs.is_empty(),
        class: CertificateClass::GridOnly,
        worst,
        crossings: Vec::new(),
        violation_freqs,
    })
}

/// Exact zero-order-hold stepping of a [`MimoModel`]: one complex (pair) or
/// real state per pole term per input port, advanced by
/// `x <- exp(p T) x + (exp(p T) - 1)/p u` (the exponential of the state
/// matrix, never a bilinear map). Output `y_n = SUM 2 Re(R x_n) (pairs),
/// R x_n (real) + D u_n + E (u_n - u_(n-1)) / T`. The improper term uses the
/// same first-order differentiator as the SISO discretization.
#[derive(Debug, Clone)]
pub struct MimoZoh {
    ports: usize,
    time_step_s: f64,
    /// Per term: `exp(p T)`, `(exp(p T) - 1)/p`, residue matrix, pair flag.
    terms: Vec<(C64, C64, Vec<C64>, bool)>,
    /// phi1(p T) maps each original pole state to its exact step mean.
    mean_state_gains: Vec<C64>,
    /// D + SUM R T phi2(p T), including conjugate pairs.
    mean_d: Vec<f64>,
    d: Vec<f64>,
    e: Vec<f64>,
    state: Vec<C64>,
    candidate_state: Vec<C64>,
    last_input: Vec<f64>,
}

// Entire phi functions. Direct subtraction loses the mean input term near zero.
// The fixed near-zero series retains complex pairs and deterministic traversal.
fn zoh_mean_phi(z: C64, exponential: C64) -> (C64, C64) {
    if z.abs() < 0.5 {
        let (mut first, mut second) = (C64::ONE, C64::from_re(0.5));
        let (mut first_term, mut second_term) = (first, second);
        for k in 1..=18 {
            first_term = (first_term * z).scale(1.0 / f64::from(k + 1));
            second_term = (second_term * z).scale(1.0 / f64::from(k + 2));
            first = first + first_term;
            second = second + second_term;
        }
        (first, second)
    } else {
        let first = (exponential - C64::ONE) * z.recip();
        (first, (first - C64::ONE) * z.recip())
    }
}

impl MimoModel {
    /// Exact ZOH stepping at `time_step_s`.
    ///
    /// # Errors
    /// [`PassivityError::Unstable`] for a pole on or right of the imaginary
    /// axis or a non-positive, non-finite step.
    pub fn zoh(&self, time_step_s: f64) -> Result<MimoZoh, PassivityError> {
        if !self.is_stable() || !(time_step_s.is_finite() && time_step_s > 0.0) {
            return Err(PassivityError::Unstable);
        }
        let mut terms = Vec::with_capacity(self.poles.len());
        let mut mean_state_gains = Vec::with_capacity(self.poles.len());
        let mut mean_d = self.d.clone();
        for (term, r) in self.poles.iter().zip(&self.residues) {
            let (pole, pair) = match *term {
                PoleTerm::Real { pole, .. } => (C64::from_re(pole), false),
                PoleTerm::Pair { pole, .. } => (pole, true),
            };
            let x = pole.scale(time_step_s);
            let (s, c) = det::sin_cos(x.im);
            let decay = C64::new(c, s).scale(det::exp(x.re));
            let (phi1, phi2) = zoh_mean_phi(x, decay);
            let mean_input = phi2.scale(time_step_s);
            for (d, residue) in mean_d.iter_mut().zip(r) {
                *d += if pair {
                    2.0 * (*residue * mean_input).re
                } else {
                    residue.re * mean_input.re
                };
            }
            mean_state_gains.push(phi1);
            terms.push((decay, phi1.scale(time_step_s), r.clone(), pair));
        }
        let n = self.ports;
        Ok(MimoZoh {
            ports: n,
            time_step_s,
            state: vec![C64::ZERO; terms.len() * n],
            candidate_state: vec![C64::ZERO; terms.len() * n],
            terms,
            mean_state_gains,
            mean_d,
            d: self.d.clone(),
            e: self.e.clone(),
            last_input: vec![0.0; n],
        })
    }
}

impl MimoZoh {
    /// Actual sampled transfer matrix at angular frequency `omega` [rad/s].
    /// This includes the held-input state update and start-of-interval output;
    /// it is not the continuous rational response relabeled as a digital filter.
    #[must_use]
    pub fn eval_discrete_iw(&self, omega: f64) -> Vec<C64> {
        self.eval_output_discrete_iw(omega, false)
    }

    /// Actual sampled transfer for the exact interval-mean proper output.
    /// The state transition is identical to `eval_discrete_iw`; only the output
    /// map changes. An improper E term retains the existing differentiated-input
    /// convention; the physical terminal adapter separately refuses that term.
    #[must_use]
    pub fn eval_mean_discrete_iw(&self, omega: f64) -> Vec<C64> {
        self.eval_output_discrete_iw(omega, true)
    }

    fn eval_output_discrete_iw(&self, omega: f64, mean: bool) -> Vec<C64> {
        let (s, c) = det::sin_cos(omega * self.time_step_s);
        let z = C64::new(c, s);
        let direct = if mean { &self.mean_d } else { &self.d };
        let mut y: Vec<C64> = direct
            .iter()
            .zip(&self.e)
            .map(|(&d, &e)| C64::from_re(d) + (C64::ONE - z.conj()).scale(e / self.time_step_s))
            .collect();
        for (k, (decay, gain, r, pair)) in self.terms.iter().enumerate() {
            let weight = if mean {
                self.mean_state_gains[k]
            } else {
                C64::ONE
            };
            let first = weight * *gain * (z - *decay).recip();
            let second = weight.conj() * gain.conj() * (z - decay.conj()).recip();
            for (y, r) in y.iter_mut().zip(r) {
                *y = *y
                    + if *pair {
                        *r * first + r.conj() * second
                    } else {
                        first.scale(r.re)
                    };
            }
        }
        y
    }

    /// Check the actual sampled impedance on an explicit uniform DC-to-Nyquist
    /// grid. A continuous Hamiltonian certificate does not certify a sampled
    /// interconnection. This diagnostic is always `GridOnly`, even if every
    /// point passes; no gain scaling or residue repair is performed here.
    ///
    /// # Panics
    /// Fewer than two grid points.
    #[must_use]
    pub fn sampled_passivity(&self, points: usize) -> MimoPassivityReport {
        self.output_sampled_passivity(points, false)
    }

    /// DC-to-Nyquist diagnostic of the exact interval-mean output that exchanges
    /// work with held input. This is still `GridOnly`; changing output quadrature
    /// does not manufacture a discrete storage certificate.
    ///
    /// # Panics
    /// Fewer than two grid points.
    #[must_use]
    pub fn mean_sampled_passivity(&self, points: usize) -> MimoPassivityReport {
        self.output_sampled_passivity(points, true)
    }

    fn output_sampled_passivity(&self, points: usize, mean: bool) -> MimoPassivityReport {
        assert!(points >= 2, "passivity grid needs DC and Nyquist");
        let mut worst = (f64::INFINITY, 0.0);
        let mut violation_freqs = Vec::new();
        for i in 0..points {
            let w = core::f64::consts::PI / self.time_step_s * (i as f64 / (points - 1) as f64);
            let response = self.eval_output_discrete_iw(w, mean);
            let lambda = if response
                .iter()
                .all(|z| z.re.is_finite() && z.im.is_finite())
            {
                hermitian_part(&response, self.ports).0[0]
            } else {
                f64::NAN
            };
            if !lambda.is_finite() || lambda < 0.0 {
                violation_freqs.push(w);
            }
            if lambda < worst.0 {
                worst = (lambda, w);
            }
        }
        MimoPassivityReport {
            passive: violation_freqs.is_empty(),
            class: CertificateClass::GridOnly,
            worst,
            crossings: Vec::new(),
            violation_freqs,
        }
    }

    pub(crate) fn mean_history_output_into(&self, y: &mut [f64]) {
        y.fill(0.0);
        let n = self.ports;
        for (k, (_, _, r, pair)) in self.terms.iter().enumerate() {
            let x = &self.state[k * n..(k + 1) * n];
            for (a, ya) in y.iter_mut().enumerate() {
                let mut acc = C64::ZERO;
                for (b, xb) in x.iter().enumerate() {
                    acc = acc + r[a * n + b] * *xb;
                }
                acc = acc * self.mean_state_gains[k];
                *ya += if *pair { 2.0 * acc.re } else { acc.re };
            }
        }
    }

    pub(crate) fn mean_direct_term(&self) -> &[f64] {
        &self.mean_d
    }

    /// Advance the same pole states and return the interval-mean output for a
    /// held input. For proper models, `dot(mean_output, input) * dt` is the exact
    /// continuous port work, up to floating-point evaluation. This is neither
    /// a trapezoid of endpoint samples nor a second propagation runtime.
    /// The improper E term follows `step`'s original differentiated-input rule;
    /// the coupled terminal owner does not admit improper models.
    ///
    /// # Panics
    /// If `u.len()` differs from the port count.
    pub fn step_mean(&mut self, u: &[f64]) -> Vec<f64> {
        assert_eq!(u.len(), self.ports, "one held input per port");
        let n = self.ports;
        let mut y = vec![0.0; n];
        self.mean_history_output_into(&mut y);
        for (a, ya) in y.iter_mut().enumerate() {
            *ya += (0..n)
                .map(|b| {
                    self.mean_d[a * n + b] * u[b]
                        + self.e[a * n + b] * (u[b] - self.last_input[b]) / self.time_step_s
                })
                .sum::<f64>();
        }
        self.preview_state(u);
        self.accept_state(u);
        y
    }

    // The network previews the SAME ZOH state update before publishing waves.
    // No heap allocation or physical mutation occurs here.
    pub(crate) fn preview_state(&mut self, u: &[f64]) -> bool {
        let n = self.ports;
        for (k, (decay, gain, _, _)) in self.terms.iter().enumerate() {
            for (b, &u) in u.iter().enumerate() {
                self.candidate_state[k * n + b] = *decay * self.state[k * n + b] + gain.scale(u);
            }
        }
        self.candidate_state
            .iter()
            .all(|z| z.re.is_finite() && z.im.is_finite())
    }
    pub(crate) fn accept_state(&mut self, u: &[f64]) {
        self.state.copy_from_slice(&self.candidate_state);
        self.last_input.copy_from_slice(u);
    }
    /// Accepted original pole states, in pole-term then input-port order.
    #[must_use]
    pub fn pole_states(&self) -> &[C64] {
        &self.state
    }

    /// Advance one sample with the held input `u` (one value per port) and
    /// return the output at the start of the interval.
    ///
    /// # Panics
    /// If `u.len()` differs from the port count.
    pub fn step(&mut self, u: &[f64]) -> Vec<f64> {
        let n = self.ports;
        assert_eq!(u.len(), n, "one held input per port");
        let mut y: Vec<f64> = (0..n)
            .map(|a| {
                (0..n)
                    .map(|b| {
                        self.d[a * n + b] * u[b]
                            + self.e[a * n + b] * (u[b] - self.last_input[b]) / self.time_step_s
                    })
                    .sum()
            })
            .collect();
        for (k, (_, _, r, pair)) in self.terms.iter().enumerate() {
            let x = &self.state[k * n..(k + 1) * n];
            for (a, ya) in y.iter_mut().enumerate() {
                let mut acc = C64::ZERO;
                for (b, xb) in x.iter().enumerate() {
                    acc = acc + r[a * n + b] * *xb;
                }
                *ya += if *pair { 2.0 * acc.re } else { acc.re };
            }
        }
        self.preview_state(u);
        self.accept_state(u);
        y
    }
}

#[cfg(test)]
#[path = "mimo_tests.rs"]
mod tests;
