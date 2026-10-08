//! Koopman / dynamic-mode-decomposition reduced models (plan §9.7: "Koopman/
//! DMD for unsteady flows"), wrapped in conformal forecast bands so they obey
//! the crate's certify-or-escalate discipline.
//!
//! - [`dmd`] — EXACT DMD (Tu et al., J. Comput. Dyn. 2014): from a snapshot
//!   sequence `x₀ … x_m` sampled every `dt`, the rank-`r` linear operator
//!   `Ã = Uᵀ Y V Σ⁻¹` on the POD subspace of `X = [x₀ … x_{m−1}]` that best
//!   maps `X` to `Y = [x₁ … x_m]`. Forecasts iterate `x ↦ U Ã Uᵀ x`.
//! - [`edmd`] — EXTENDED DMD (Williams–Kevrekidis–Rowley 2015) with a
//!   polynomial dictionary: a finite-section Koopman matrix `K` with
//!   `ψ(x_{k+1}) ≈ K ψ(x_k)`, ridge-regularized least squares. The state is
//!   read back from the dictionary's linear monomials; multi-step forecasts
//!   re-lift every step.
//! - Spectra: discrete eigenvalues `λ` (complex, via Hessenberg + Francis QR)
//!   and continuous-time `(growth, frequency) = (ln|λ|, arg λ)/dt`.
//! - [`forecast_bands`] — per-horizon split-conformal bands on the forecast
//!   error norm over calibration trajectories (exchangeable initial
//!   conditions), optionally Bonferroni-simultaneous over the horizon.
//!
//! No-claim: DMD/EDMD are LEARNED models. Their spectra describe the fitted
//! operator, not the physics, unless the dynamics are (finite-section)
//! linear in the dictionary; the conformal band is the only accuracy
//! statement, it is marginal over the calibration distribution, and it says
//! nothing outside that distribution — consumers escalate when the query
//! leaves the calibrated regime (see [`crate::certify_or_escalate`]).

use crate::linalg;
use crate::{ConformalBand, SurrogateError, conformal_band, jacobi_eig};

/// A model that forecasts a state trajectory.
pub trait Forecaster {
    /// The predicted states `x̂₁ … x̂_steps` from the initial state `x0`.
    fn forecast(&self, x0: &[f64], steps: usize) -> Vec<Vec<f64>>;
}

/// How many POD modes the DMD keeps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DmdRank {
    /// Exactly this many (capped by the numerical rank).
    Fixed(usize),
    /// The fewest modes capturing this fraction of `X`'s energy, in `(0, 1]`.
    Energy(f64),
}

/// An exact-DMD model.
#[derive(Debug, Clone, PartialEq)]
pub struct Dmd {
    /// Orthonormal POD modes of `X` (each of state dimension).
    basis: Vec<Vec<f64>>,
    /// Reduced operator `Ã` (`r × r`, row-major).
    a_tilde: Vec<f64>,
    dt: f64,
    eigenvalues: Vec<(f64, f64)>,
    fit_residual: f64,
}

fn validate(snapshots: &[Vec<f64>]) -> Result<usize, SurrogateError> {
    let Some(first) = snapshots.first() else {
        return Err(SurrogateError::NoSnapshots);
    };
    let n = first.len();
    for s in snapshots {
        if s.len() != n {
            return Err(SurrogateError::DimMismatch {
                expected: n,
                found: s.len(),
            });
        }
    }
    Ok(n)
}

/// Fit exact DMD to a uniformly sampled snapshot sequence.
///
/// # Errors
/// [`SurrogateError`] for fewer than two snapshots, ragged snapshots, a bad
/// rank request or step, or a numerically zero snapshot matrix.
pub fn dmd(snapshots: &[Vec<f64>], dt: f64, rank: DmdRank) -> Result<Dmd, SurrogateError> {
    let n = validate(snapshots)?;
    if snapshots.len() < 2 {
        return Err(SurrogateError::NoSnapshots);
    }
    if !(dt.is_finite() && dt > 0.0) {
        return Err(SurrogateError::BadThreshold);
    }
    match rank {
        DmdRank::Fixed(0) => return Err(SurrogateError::BadThreshold),
        DmdRank::Energy(e) if !(e > 0.0 && e <= 1.0) => return Err(SurrogateError::BadThreshold),
        _ => {}
    }
    let m = snapshots.len() - 1;
    let x = &snapshots[..m];
    let y = &snapshots[1..];
    // Method of snapshots: XᵀX = V Σ² Vᵀ.
    let mut c = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in i..m {
            let d: f64 = x[i].iter().zip(&x[j]).map(|(a, b)| a * b).sum();
            c[i][j] = d;
            c[j][i] = d;
        }
    }
    let (lambda, vecs) = jacobi_eig(c);
    let lmax = lambda.first().copied().unwrap_or(0.0);
    if !(lmax > 0.0) {
        return Err(SurrogateError::NoSnapshots);
    }
    let numerical: Vec<usize> = (0..m).filter(|&k| lambda[k] > 1e-20 * lmax).collect();
    let total: f64 = numerical.iter().map(|&k| lambda[k]).sum();
    let r = match rank {
        DmdRank::Fixed(r) => r.min(numerical.len()),
        DmdRank::Energy(e) => {
            let mut cum = 0.0;
            let mut r = 0;
            for &k in &numerical {
                cum += lambda[k];
                r += 1;
                if cum / total >= e {
                    break;
                }
            }
            r
        }
    };
    // U = X V Σ⁻¹ (n × r), stored as r modes.
    let sig: Vec<f64> = (0..r).map(|k| lambda[k].sqrt()).collect();
    let basis: Vec<Vec<f64>> = (0..r)
        .map(|k| {
            (0..n)
                .map(|row| (0..m).map(|i| x[i][row] * vecs[k][i]).sum::<f64>() / sig[k])
                .collect()
        })
        .collect();
    // Ã = Uᵀ Y V Σ⁻¹: entry (a, b) = Σ_i (U_aᵀ y_i) V[i][b] / σ_b.
    let uty: Vec<Vec<f64>> = basis
        .iter()
        .map(|u| {
            y.iter()
                .map(|yi| u.iter().zip(yi).map(|(p, q)| p * q).sum())
                .collect()
        })
        .collect();
    let mut a_tilde = vec![0.0; r * r];
    for a in 0..r {
        for b in 0..r {
            a_tilde[a * r + b] = (0..m).map(|i| uty[a][i] * vecs[b][i]).sum::<f64>() / sig[b];
        }
    }
    let eigenvalues = linalg::eigenvalues_real(&a_tilde, r).unwrap_or_default();
    let mut model = Dmd {
        basis,
        a_tilde,
        dt,
        eigenvalues,
        fit_residual: 0.0,
    };
    // Relative one-step residual over the training pairs.
    let (mut num, mut den) = (0.0, 0.0);
    for (xi, yi) in x.iter().zip(y) {
        let p = model.step(xi);
        num += p
            .iter()
            .zip(yi)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>();
        den += yi.iter().map(|b| b * b).sum::<f64>();
    }
    model.fit_residual = if den > 0.0 { (num / den).sqrt() } else { 0.0 };
    Ok(model)
}

fn continuous(ev: &[(f64, f64)], dt: f64) -> Vec<(f64, f64)> {
    ev.iter()
        .map(|&(re, im)| {
            let modulus = re.hypot(im);
            (modulus.ln() / dt, im.atan2(re) / dt)
        })
        .collect()
}

impl Dmd {
    /// Retained rank.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.basis.len()
    }

    /// Sampling interval.
    #[must_use]
    pub fn dt(&self) -> f64 {
        self.dt
    }

    /// Discrete-time eigenvalues `(re, im)` of `Ã`, descending modulus.
    #[must_use]
    pub fn eigenvalues(&self) -> &[(f64, f64)] {
        &self.eigenvalues
    }

    /// Continuous-time `(growth rate, angular frequency)` per eigenvalue.
    #[must_use]
    pub fn continuous_spectrum(&self) -> Vec<(f64, f64)> {
        continuous(&self.eigenvalues, self.dt)
    }

    /// Largest eigenvalue modulus (`< 1` ⇔ the FITTED operator is stable).
    #[must_use]
    pub fn spectral_radius(&self) -> f64 {
        self.eigenvalues
            .iter()
            .fold(0.0f64, |m, (re, im)| m.max(re.hypot(*im)))
    }

    /// Relative one-step residual on the training pairs.
    #[must_use]
    pub fn fit_residual(&self) -> f64 {
        self.fit_residual
    }

    /// One step `x ↦ U Ã Uᵀ x`.
    #[must_use]
    pub fn step(&self, x: &[f64]) -> Vec<f64> {
        let r = self.rank();
        let z: Vec<f64> = self
            .basis
            .iter()
            .map(|u| u.iter().zip(x).map(|(a, b)| a * b).sum())
            .collect();
        let zn: Vec<f64> = (0..r)
            .map(|a| (0..r).map(|b| self.a_tilde[a * r + b] * z[b]).sum())
            .collect();
        let n = x.len();
        let mut out = vec![0.0; n];
        for (c, u) in zn.iter().zip(&self.basis) {
            for (o, ui) in out.iter_mut().zip(u) {
                *o += c * ui;
            }
        }
        out
    }
}

impl Forecaster for Dmd {
    fn forecast(&self, x0: &[f64], steps: usize) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(steps);
        let mut x = x0.to_vec();
        for _ in 0..steps {
            x = self.step(&x);
            out.push(x.clone());
        }
        out
    }
}

/// An EDMD model over a polynomial dictionary.
#[derive(Debug, Clone, PartialEq)]
pub struct Edmd {
    nvars: usize,
    dictionary: Vec<Vec<u32>>,
    /// `N × N` Koopman matrix, row-major: `ψ(x⁺) ≈ K ψ(x)`.
    k: Vec<f64>,
    linear_index: Vec<usize>,
    dt: f64,
    eigenvalues: Vec<(f64, f64)>,
}

fn dictionary(nvars: usize, degree: u32) -> Vec<Vec<u32>> {
    fn rec(k: usize, left: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if k == cur.len() {
            out.push(cur.clone());
            return;
        }
        for e in 0..=left {
            cur[k] = e;
            rec(k + 1, left - e, cur, out);
        }
        cur[k] = 0;
    }
    let mut out = Vec::new();
    let mut cur = vec![0; nvars];
    rec(0, degree, &mut cur, &mut out);
    // Graded order: degree, then reverse-lex (1, x₀, x₁, …, x₀², …).
    out.sort_by(|a, b| {
        let (da, db): (u32, u32) = (a.iter().sum(), b.iter().sum());
        da.cmp(&db).then_with(|| b.cmp(a))
    });
    out
}

/// Fit EDMD with all monomials of total degree `≤ degree` from snapshot
/// pairs `(x_k, x_{k+1})`; `ridge` is relative Tikhonov regularization of the
/// Gram matrix (escalated automatically if it is singular).
///
/// # Errors
/// [`SurrogateError`] for empty/ragged data, `degree == 0`, or a bad step.
pub fn edmd(
    pairs: &[(Vec<f64>, Vec<f64>)],
    degree: u32,
    ridge: f64,
    dt: f64,
) -> Result<Edmd, SurrogateError> {
    let Some((x0, _)) = pairs.first() else {
        return Err(SurrogateError::NoSnapshots);
    };
    let nvars = x0.len();
    for (a, b) in pairs {
        for v in [a, b] {
            if v.len() != nvars {
                return Err(SurrogateError::DimMismatch {
                    expected: nvars,
                    found: v.len(),
                });
            }
        }
    }
    if degree == 0 || nvars == 0 || !(dt.is_finite() && dt > 0.0) || !(ridge >= 0.0) {
        return Err(SurrogateError::BadThreshold);
    }
    let dict = dictionary(nvars, degree);
    let nd = dict.len();
    let lift = |x: &[f64]| -> Vec<f64> {
        dict.iter()
            .map(|e| e.iter().zip(x).map(|(&p, &xi)| xi.powi(p as i32)).product())
            .collect()
    };
    // G = Σ ψ(x)ψ(x)ᵀ, A = Σ ψ(y)ψ(x)ᵀ; K = A G⁻¹ ⇔ G Kᵀ = Aᵀ.
    let mut g = vec![0.0; nd * nd];
    let mut a = vec![0.0; nd * nd];
    for (x, y) in pairs {
        let px = lift(x);
        let py = lift(y);
        for i in 0..nd {
            for j in 0..nd {
                g[i * nd + j] += px[i] * px[j];
                a[i * nd + j] += py[i] * px[j];
            }
        }
    }
    let mut k = vec![0.0; nd * nd];
    for row in 0..nd {
        // Row `row` of K solves G k_row = A[row, :]ᵀ (G symmetric).
        let rhs: Vec<f64> = (0..nd).map(|j| a[row * nd + j]).collect();
        let sol = linalg::ridge_solve(&g, nd, &rhs, ridge).ok_or(SurrogateError::NoSnapshots)?;
        k[row * nd..(row + 1) * nd].copy_from_slice(&sol);
    }
    let linear_index: Vec<usize> = (0..nvars)
        .map(|v| {
            dict.iter()
                .position(|e| e.iter().sum::<u32>() == 1 && e[v] == 1)
                .unwrap_or_else(|| unreachable!("degree >= 1 contains every linear monomial"))
        })
        .collect();
    let eigenvalues = linalg::eigenvalues_real(&k, nd).unwrap_or_default();
    Ok(Edmd {
        nvars,
        dictionary: dict,
        k,
        linear_index,
        dt,
        eigenvalues,
    })
}

impl Edmd {
    /// Dictionary size `N`.
    #[must_use]
    pub fn dictionary_size(&self) -> usize {
        self.dictionary.len()
    }

    /// The dictionary's exponent vectors (graded order).
    #[must_use]
    pub fn dictionary(&self) -> &[Vec<u32>] {
        &self.dictionary
    }

    /// The Koopman matrix (`N × N`, row-major).
    #[must_use]
    pub fn koopman_matrix(&self) -> &[f64] {
        &self.k
    }

    /// Discrete eigenvalues of `K`, descending modulus.
    #[must_use]
    pub fn eigenvalues(&self) -> &[(f64, f64)] {
        &self.eigenvalues
    }

    /// Continuous-time `(growth rate, angular frequency)` per eigenvalue.
    #[must_use]
    pub fn continuous_spectrum(&self) -> Vec<(f64, f64)> {
        continuous(&self.eigenvalues, self.dt)
    }

    /// Lift a state into the dictionary.
    #[must_use]
    pub fn lift(&self, x: &[f64]) -> Vec<f64> {
        self.dictionary
            .iter()
            .map(|e| e.iter().zip(x).map(|(&p, &xi)| xi.powi(p as i32)).product())
            .collect()
    }

    /// One step: lift, apply `K`, read the linear observables.
    #[must_use]
    pub fn step(&self, x: &[f64]) -> Vec<f64> {
        let psi = self.lift(x);
        let nd = self.dictionary.len();
        self.linear_index
            .iter()
            .map(|&row| (0..nd).map(|j| self.k[row * nd + j] * psi[j]).sum())
            .collect()
    }

    /// Number of state variables.
    #[must_use]
    pub fn nvars(&self) -> usize {
        self.nvars
    }
}

impl Forecaster for Edmd {
    fn forecast(&self, x0: &[f64], steps: usize) -> Vec<Vec<f64>> {
        let mut out = Vec::with_capacity(steps);
        let mut x = x0.to_vec();
        for _ in 0..steps {
            x = self.step(&x);
            out.push(x.clone());
        }
        out
    }
}

/// Per-horizon conformal bands on the forecast error norm `‖x̂_h − x_h‖₂`.
#[derive(Debug, Clone, PartialEq)]
pub struct ForecastBands {
    /// `bands[h−1]` covers horizon `h`.
    pub bands: Vec<ConformalBand>,
    /// Whether the bands were Bonferroni-adjusted to hold SIMULTANEOUSLY
    /// over all horizons at the requested level.
    pub simultaneous: bool,
}

impl ForecastBands {
    /// Does a forecast trajectory lie inside the bands around the truth?
    #[must_use]
    pub fn covers(&self, forecast: &[Vec<f64>], truth: &[Vec<f64>]) -> bool {
        self.bands
            .iter()
            .zip(forecast.iter().zip(truth))
            .all(|(b, (f, t))| {
                let e = f
                    .iter()
                    .zip(t)
                    .map(|(a, c)| (a - c) * (a - c))
                    .sum::<f64>()
                    .sqrt();
                e <= b.half_width
            })
    }

    /// The widest band (the decision-relevant one for whole-horizon use).
    #[must_use]
    pub fn max_half_width(&self) -> f64 {
        self.bands.iter().fold(0.0f64, |m, b| m.max(b.half_width))
    }
}

/// Calibrate per-horizon split-conformal bands for `model` on independent
/// calibration trajectories (each `traj[0]` is an initial state, `traj[h]` the
/// truth at horizon `h`). With `simultaneous`, each horizon uses `α/H` so the
/// whole forecast is covered with probability `≥ 1 − α`.
///
/// # Panics
/// If no trajectory reaches `horizon` or `alpha ∉ (0, 1)`.
#[must_use]
pub fn forecast_bands<F: Forecaster>(
    model: &F,
    calibration: &[Vec<Vec<f64>>],
    horizon: usize,
    alpha: f64,
    simultaneous: bool,
) -> ForecastBands {
    assert!(horizon >= 1, "horizon must be at least 1");
    let usable: Vec<&Vec<Vec<f64>>> = calibration.iter().filter(|t| t.len() > horizon).collect();
    assert!(
        !usable.is_empty(),
        "no calibration trajectory reaches the horizon"
    );
    let a = if simultaneous {
        alpha / horizon as f64
    } else {
        alpha
    };
    let forecasts: Vec<Vec<Vec<f64>>> = usable
        .iter()
        .map(|t| model.forecast(&t[0], horizon))
        .collect();
    let bands = (1..=horizon)
        .map(|h| {
            let residuals: Vec<f64> = usable
                .iter()
                .zip(&forecasts)
                .map(|(t, f)| {
                    f[h - 1]
                        .iter()
                        .zip(&t[h])
                        .map(|(p, q)| (p - q) * (p - q))
                        .sum::<f64>()
                        .sqrt()
                })
                .collect();
            conformal_band(&residuals, a)
        })
        .collect();
    ForecastBands {
        bands,
        simultaneous,
    }
}
