//! Finite-strain Prony branches (plan M0.3-W1b): a hyperelastic equilibrium
//! law plus Maxwell branches, each a multiplicative split `F = F_e F_v` with
//! isotropic linear evolution (Reese & Govindjee 1998, Int J Solids Struct
//! 35:3455). A branch relaxes its isochoric elastic stretches; its stress
//! at the end of a step is the gradient of an incremental potential, so the
//! algorithmic tangent is its (symmetric) Hessian.
//!
//! Integration is backward Euler in the principal logarithmic elastic
//! stretches `eps` (the exponential map), with `det F_v = 1` held by
//! construction. On the volume-preserving plane the local problem is
//!
//! `eps* = argmin psi(eps) + (G tau / h) |eps - eps_trial|^2`,
//!
//! `G` the branch's small-strain shear modulus, `tau` its relaxation time
//! and `h` the step. It is strictly convex, so a damped Newton iteration
//! converges globally; a stall or an unattained tolerance is a typed
//! refusal, never a best iterate. The value function's gradient and Hessian
//! in the trial stretches feed the same second-order spectral lift as the
//! Ogden law (Simo & Taylor 1991): no eigenvector derivative, repeated
//! stretches included. The divided difference of the gradient at close or
//! equal roots is the integral of a Hessian entry along the swap path
//! (five-point Gauss rule), so it has no cancellation.
//!
//! Energy ledger, booked as the exact frozen-state split of each step:
//! loading at frozen internal state (`work`), then relaxation at the new
//! `F` (`dissipated`, a sum of non-negative terms computed directly, not as
//! a difference of energies). `stored` is the sum of the equilibrium and
//! branch energies; `work - Delta stored - dissipated` therefore measures
//! the local solve's consistency, not a tautology.

use crate::hyper::{Hyperelastic, MAX_OGDEN_TERMS, OgdenTerm};
use crate::{MaterialAdmissibility, MaterialError};
use fs_ad::{Dual, Dual64, Real};
use fs_la::eigen::jacobi_eigh;
use fs_math::det;

/// Maximum number of Maxwell branches (as `ThermoMaxwell`).
pub const MAX_FINITE_STRAIN_BRANCHES: usize = 64;
/// Default relative residual of the local Newton iteration.
pub const DEFAULT_LOCAL_NEWTON_TOLERANCE: f64 = 1e-13;
/// Default cap on the local Newton iterations.
pub const DEFAULT_LOCAL_NEWTON_ITERATIONS: usize = 64;

// Relative root gap below which the divided difference of the principal
// gradient is evaluated as a path integral instead of a difference.
const COALESCENCE_GAP: f64 = 0.1;
const S2: f64 = core::f64::consts::FRAC_1_SQRT_2;
const S6: f64 = 0.408_248_290_463_863_1;
// Five-point Gauss-Legendre rule on [-1, 1].
const GL5_XI: [f64; 5] = [
    -0.906_179_845_938_664,
    -0.538_469_310_105_683_1,
    0.0,
    0.538_469_310_105_683_1,
    0.906_179_845_938_664,
];
const GL5_W: [f64; 5] = [
    0.236_926_885_056_189_1,
    0.478_628_670_499_366_5,
    0.568_888_888_888_888_9,
    0.478_628_670_499_366_5,
    0.236_926_885_056_189_1,
];

/// Isochoric non-equilibrium energy of one branch, per reference volume.
#[derive(Debug, Clone, PartialEq)]
pub enum BranchEnergy {
    /// `mu/2 (tr b_e_bar - 3)`.
    NeoHookean {
        /// Shear modulus (Pa).
        mu_pa: f64,
    },
    /// Decoupled Ogden terms in the `mu/alpha` convention of the
    /// equilibrium law.
    Ogden {
        /// The terms (at most [`MAX_OGDEN_TERMS`]).
        terms: Vec<OgdenTerm>,
    },
}

/// One Maxwell branch: an isochoric energy and a relaxation time.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteStrainBranch {
    /// The non-equilibrium energy.
    pub energy: BranchEnergy,
    /// Relaxation time (s), resolved at the state point.
    pub tau_s: f64,
}

/// Internal variables and energy ledger of a [`ViscoHyperelastic`] history.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoState {
    /// `C_v^{-1}` per branch, `[xx, yy, zz, xy, xz, yz]`, determinant 1.
    pub inverse_viscous_cauchy_green: Vec<[f64; 6]>,
    /// The deformation gradient the state was last advanced to.
    pub deformation_gradient: [f64; 9],
    /// Accumulated loading work at frozen internal state (J/m^3).
    pub work: f64,
    /// Current stored energy, equilibrium plus branches (J/m^3).
    pub stored: f64,
    /// Accumulated relaxation dissipation (J/m^3), non-negative and
    /// non-decreasing.
    pub dissipated: f64,
}

/// A hyperelastic equilibrium law with finite-strain Maxwell branches.
#[derive(Debug, Clone, PartialEq)]
pub struct ViscoHyperelastic {
    /// The equilibrium law (decoupled: isochoric energy plus `U(J)`).
    pub equilibrium: Hyperelastic,
    /// The Maxwell branches (at most [`MAX_FINITE_STRAIN_BRANCHES`]).
    pub branches: Vec<FiniteStrainBranch>,
    newton_tolerance: f64,
    newton_max_iterations: usize,
}

fn parameters(what: String) -> MaterialError {
    MaterialError::Parameters { what }
}

fn state_error(what: String) -> MaterialError {
    MaterialError::State { what }
}

impl BranchEnergy {
    /// Small-strain shear modulus `G` (Pa): `mu` or `sum mu_p alpha_p / 2`.
    #[must_use]
    pub fn shear_modulus_pa(&self) -> f64 {
        match self {
            Self::NeoHookean { mu_pa } => *mu_pa,
            Self::Ogden { terms } => 0.5 * terms.iter().map(|t| t.mu_pa * t.alpha).sum::<f64>(),
        }
    }

    fn validate(&self) -> Result<(), MaterialError> {
        let ok = match self {
            Self::NeoHookean { mu_pa } => mu_pa.is_finite() && *mu_pa > 0.0,
            Self::Ogden { terms } => {
                !terms.is_empty()
                    && terms.len() <= MAX_OGDEN_TERMS
                    && terms.iter().all(|t| {
                        t.mu_pa.is_finite()
                            && t.alpha.is_finite()
                            && (t.mu_pa / t.alpha).is_finite()
                            && (t.mu_pa * t.alpha).is_finite()
                            && t.mu_pa * t.alpha > 0.0
                    })
                    && self.shear_modulus_pa().is_finite()
                    && self.shear_modulus_pa() > 0.0
            }
        };
        if ok {
            Ok(())
        } else {
            Err(parameters(format!(
                "{self:?}: branch moduli must be finite and positive (Ogden: 1..={MAX_OGDEN_TERMS} terms, mu*alpha > 0)"
            )))
        }
    }

    // Energy and the diagonal of its Hessian in the three principal
    // logarithmic stretches x (full, unprojected).
    fn terms(&self, x: &[f64; 3]) -> (f64, [f64; 3]) {
        let mut psi = 0.0;
        let mut hess = [0.0; 3];
        match self {
            Self::NeoHookean { mu_pa } => {
                for a in 0..3 {
                    psi += 0.5 * mu_pa * det::expm1(2.0 * x[a]);
                    hess[a] = 2.0 * mu_pa * det::exp(2.0 * x[a]);
                }
            }
            Self::Ogden { terms } => {
                for t in terms {
                    for a in 0..3 {
                        psi += t.mu_pa / t.alpha * det::expm1(t.alpha * x[a]);
                        hess[a] += t.mu_pa * t.alpha * det::exp(t.alpha * x[a]);
                    }
                }
            }
        }
        (psi, hess)
    }

    // The gradient of the energy projected on the volume-preserving plane,
    // E^T grad psi, written without cancellation: each exponential
    // difference is an expm1 of a stretch difference, so the result keeps
    // its relative accuracy as it tends to zero. A step that is huge against
    // tau drives the relaxed stretches to ~1e-17, where a difference of
    // exponentials would carry an absolute error of mu * eps and, times the
    // step factor c, an unreachable Newton residual.
    fn reduced_gradient(&self, x: &[f64; 3]) -> [f64; 2] {
        let mut g = [0.0; 2];
        let mut add = |coefficient: f64, rate: f64| {
            g[0] += S2 * coefficient * det::exp(rate * x[1]) * det::expm1(rate * (x[0] - x[1]));
            g[1] += S6
                * coefficient
                * det::exp(rate * x[2])
                * (det::expm1(rate * (x[0] - x[2])) + det::expm1(rate * (x[1] - x[2])));
        };
        match self {
            Self::NeoHookean { mu_pa } => add(*mu_pa, 2.0),
            Self::Ogden { terms } => {
                for t in terms {
                    add(t.mu_pa, t.alpha);
                }
            }
        }
        g
    }

    // psi(b + delta) - psi(b) - grad psi(b) . delta, term by term, without
    // subtracting nearly equal energies. Non-negative for convex psi.
    fn bregman(&self, b: &[f64; 3], delta: &[f64; 3]) -> f64 {
        let mut sum = 0.0;
        match self {
            Self::NeoHookean { mu_pa } => {
                for a in 0..3 {
                    sum += 0.5 * mu_pa * det::exp(2.0 * b[a]) * expm1_minus_x(2.0 * delta[a]);
                }
            }
            Self::Ogden { terms } => {
                for t in terms {
                    for a in 0..3 {
                        sum += t.mu_pa / t.alpha
                            * det::exp(t.alpha * b[a])
                            * expm1_minus_x(t.alpha * delta[a]);
                    }
                }
            }
        }
        sum
    }
}

// exp(x) - 1 - x, accurate for small |x| (series) and large |x| (expm1).
fn expm1_minus_x(x: f64) -> f64 {
    if x.abs() < 0.25 {
        let mut term = 0.5 * x * x;
        let mut sum = term;
        for n in 3..40u32 {
            term *= x / f64::from(n);
            sum += term;
            if term.abs() <= 1e-18 * sum.abs() {
                break;
            }
        }
        sum
    } else {
        det::expm1(x) - x
    }
}

// Orthonormal basis of the volume-preserving plane: x = E eta.
fn expand(eta: [f64; 2]) -> [f64; 3] {
    [
        S2 * eta[0] + S6 * eta[1],
        -S2 * eta[0] + S6 * eta[1],
        -2.0 * S6 * eta[1],
    ]
}

// E^T x (also discards the mean).
fn reduce(x: [f64; 3]) -> [f64; 2] {
    [S2 * (x[0] - x[1]), S6 * (x[0] + x[1] - 2.0 * x[2])]
}

// E^T diag(h) E.
fn reduced_hessian(h: &[f64; 3]) -> [[f64; 2]; 2] {
    let off = S2 * S6 * (h[0] - h[1]);
    [
        [0.5 * (h[0] + h[1]), off],
        [off, (h[0] + h[1] + 4.0 * h[2]) / 6.0],
    ]
}

// E m E^T for a symmetric 2x2.
fn embed_symmetric(m: &[[f64; 2]; 2]) -> [[f64; 3]; 3] {
    let e = [[S2, S6], [-S2, S6], [0.0, -2.0 * S6]];
    let mut out = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let mut s = 0.0;
            for p in 0..2 {
                for q in 0..2 {
                    s += e[a][p] * m[p][q] * e[b][q];
                }
            }
            out[a][b] = s;
        }
    }
    out
}

struct Local {
    eta: [f64; 2],
    psi_tr: f64,
    psi_star: f64,
    // The step's dissipation, a sum of non-negative terms.
    dissipation: f64,
    // Principal (volume-preserving) Kirchhoff stress of the branch.
    tau: [f64; 3],
    // d tau / d eps_trial, symmetric, supported on the plane.
    k: [[f64; 3]; 3],
    // Value of the incremental potential.
    pi: f64,
}

fn norm_inf(v: [f64; 2]) -> f64 {
    v[0].abs().max(v[1].abs())
}

fn local_solve(
    energy: &BranchEnergy,
    tau_s: f64,
    dt: f64,
    eta_tr: [f64; 2],
    tolerance: f64,
    max_iterations: usize,
) -> Result<Local, MaterialError> {
    let x_tr = expand(eta_tr);
    let (psi_tr, hess_tr) = energy.terms(&x_tr);
    let g_tr = energy.reduced_gradient(&x_tr);
    if !(psi_tr.is_finite() && g_tr.iter().all(|g| g.is_finite())) {
        return Err(state_error(
            "branch energy overflowed at the trial stretches".to_string(),
        ));
    }
    if dt == 0.0 {
        let tau = expand(g_tr);
        return Ok(Local {
            eta: eta_tr,
            psi_tr,
            psi_star: psi_tr,
            dissipation: 0.0,
            tau,
            k: embed_symmetric(&reduced_hessian(&hess_tr)),
            pi: psi_tr,
        });
    }
    let g_mod = energy.shear_modulus_pa();
    let c = dt / (2.0 * g_mod * tau_s);
    if !(c.is_finite() && c > 0.0) {
        return Err(state_error(format!(
            "step {dt} s against relaxation time {tau_s} s is outside the local solve's range"
        )));
    }
    // eval: energy terms and residual r = eta - eta_tr + c g(eta).
    let eval = |eta: [f64; 2]| -> Result<(f64, [f64; 2], [f64; 3], [f64; 2]), MaterialError> {
        let x = expand(eta);
        let (psi, hess) = energy.terms(&x);
        let g2 = energy.reduced_gradient(&x);
        let r = [
            eta[0] - eta_tr[0] + c * g2[0],
            eta[1] - eta_tr[1] + c * g2[1],
        ];
        if psi.is_finite() && r.iter().all(|v| v.is_finite()) {
            Ok((psi, g2, hess, r))
        } else {
            Err(state_error(
                "branch energy overflowed in the local solve".to_string(),
            ))
        }
    };
    let newton_step = |h: &[f64; 3], r: [f64; 2]| -> Result<[f64; 2], MaterialError> {
        let hr = reduced_hessian(h);
        let n = [
            [1.0 + c * hr[0][0], c * hr[0][1]],
            [c * hr[1][0], 1.0 + c * hr[1][1]],
        ];
        let d = n[0][0] * n[1][1] - n[0][1] * n[1][0];
        if !(d.is_finite() && d > 0.0) {
            return Err(state_error(
                "local Newton matrix is not positive definite".to_string(),
            ));
        }
        Ok([
            -(n[1][1] * r[0] - n[0][1] * r[1]) / d,
            -(-n[1][0] * r[0] + n[0][0] * r[1]) / d,
        ])
    };

    let mut eta = eta_tr;
    let (mut psi, mut g2, mut hess, mut r) = eval(eta)?;
    let mut polished = false;
    let mut iteration = 0usize;
    loop {
        let scale = norm_inf(eta_tr).max(norm_inf(eta));
        let rnorm = norm_inf(r);
        if rnorm <= tolerance * scale || rnorm == 0.0 {
            if polished {
                break;
            }
            // One more full step: quadratic convergence takes the residual
            // to roundoff, which keeps the finite-difference tangent clean.
            polished = true;
            let delta = newton_step(&hess, r)?;
            let cand = [eta[0] + delta[0], eta[1] + delta[1]];
            if let Ok((p2, g2c, h2, r2)) = eval(cand) {
                if norm_inf(r2) <= rnorm {
                    eta = cand;
                    psi = p2;
                    g2 = g2c;
                    hess = h2;
                    r = r2;
                }
            }
            continue;
        }
        if iteration >= max_iterations {
            return Err(state_error(format!(
                "local Newton did not reach relative residual {tolerance:e} in {max_iterations} iterations (residual {rnorm:e})"
            )));
        }
        iteration += 1;
        let delta = newton_step(&hess, r)?;
        let mut t = 1.0;
        let mut accepted = false;
        for _ in 0..50 {
            let cand = [eta[0] + t * delta[0], eta[1] + t * delta[1]];
            if let Ok((p2, g2c, h2, r2)) = eval(cand) {
                if norm_inf(r2) < rnorm {
                    eta = cand;
                    psi = p2;
                    g2 = g2c;
                    hess = h2;
                    r = r2;
                    accepted = true;
                    break;
                }
            }
            t *= 0.5;
        }
        if !accepted {
            return Err(state_error(format!(
                "local Newton stalled at residual {rnorm:e} before reaching {tolerance:e}"
            )));
        }
    }

    let tau = expand(g2);
    // eps_trial - eps* = c * tau exactly at the stationary point.
    let delta = [c * tau[0], c * tau[1], c * tau[2]];
    let x_star = expand(eta);
    let dissipation = c * (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2])
        + energy.bregman(&x_star, &delta);
    // K = H (I + c H)^-1 on the plane.
    let hr = reduced_hessian(&hess);
    let n = [
        [1.0 + c * hr[0][0], c * hr[0][1]],
        [c * hr[1][0], 1.0 + c * hr[1][1]],
    ];
    let d = n[0][0] * n[1][1] - n[0][1] * n[1][0];
    let ninv = [[n[1][1] / d, -n[0][1] / d], [-n[1][0] / d, n[0][0] / d]];
    let mut k2 = [[0.0; 2]; 2];
    for p in 0..2 {
        for q in 0..2 {
            k2[p][q] = hr[p][0] * ninv[0][q] + hr[p][1] * ninv[1][q];
        }
    }
    let off = 0.5 * (k2[0][1] + k2[1][0]);
    k2[0][1] = off;
    k2[1][0] = off;
    let pi = psi
        + ((eta[0] - eta_tr[0]).powi(2) + (eta[1] - eta_tr[1]).powi(2)) / (2.0 * c);
    Ok(Local {
        eta,
        psi_tr,
        psi_star: psi,
        dissipation,
        tau,
        k: embed_symmetric(&k2),
        pi,
    })
}

fn det3<T: Real>(f: &[T; 9]) -> T {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn mat_mul(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        let (i, j) = (ij / 3, ij % 3);
        (0..3).map(|k| a[3 * i + k] * b[3 * k + j]).sum()
    })
}

fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| a[3 * (ij % 3) + ij / 3])
}

fn inverse3(a: &[f64; 9]) -> Option<[f64; 9]> {
    let d = det3(a);
    if !(d.is_finite() && d != 0.0) {
        return None;
    }
    let cof = |r0: usize, c0: usize, r1: usize, c1: usize| {
        a[3 * r0 + c0] * a[3 * r1 + c1] - a[3 * r0 + c1] * a[3 * r1 + c0]
    };
    Some([
        cof(1, 1, 2, 2) / d,
        -cof(0, 1, 2, 2) / d,
        cof(0, 1, 1, 2) / d,
        -cof(1, 0, 2, 2) / d,
        cof(0, 0, 2, 2) / d,
        -cof(0, 0, 1, 2) / d,
        cof(1, 0, 2, 1) / d,
        -cof(0, 0, 2, 1) / d,
        cof(0, 0, 1, 1) / d,
    ])
}

fn sym_expand(s: &[f64; 6]) -> [f64; 9] {
    [s[0], s[3], s[4], s[3], s[1], s[5], s[4], s[5], s[2]]
}

fn sym_pack(m: &[f64; 9]) -> [f64; 6] {
    [
        m[0],
        m[4],
        m[8],
        0.5 * (m[1] + m[3]),
        0.5 * (m[2] + m[6]),
        0.5 * (m[5] + m[7]),
    ]
}

// J^{-2/3} F C_v^{-1} F^T, generic so the lift can differentiate it.
fn trial_stretch_tensor<T: Real>(f: &[T; 9], cv: &[f64; 9]) -> [T; 9] {
    let j = det3(f);
    let jm23 = (T::from_f64(-2.0 / 3.0) * j.ln()).exp();
    let m: [T; 9] = core::array::from_fn(|ij| {
        let (r, c) = (ij / 3, ij % 3);
        (0..3).fold(T::zero(), |s, k| s + f[3 * r + k] * T::from_f64(cv[3 * k + c]))
    });
    core::array::from_fn(|ij| {
        let (r, c) = (ij / 3, ij % 3);
        jm23 * (0..3).fold(T::zero(), |s, k| s + m[3 * r + k] * f[3 * c + k])
    })
}

// Second-order spectral lift of one branch's incremental potential.
struct Lift {
    v: [f64; 9],
    ct: [f64; 9],
    pi0: f64,
    omega: [f64; 3],
    w: [[f64; 3]; 3],
    dd: [[f64; 3]; 3],
}

impl Lift {
    fn evaluate<T: Real>(&self, ct: &[T; 9]) -> T {
        let delta: [T; 9] = core::array::from_fn(|ij| ct[ij] - T::from_f64(self.ct[ij]));
        let h: [T; 9] = core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            let mut value = T::zero();
            for p in 0..3 {
                for q in 0..3 {
                    value = value
                        + T::from_f64(self.v[3 * p + i] * self.v[3 * q + j]) * delta[3 * p + q];
                }
            }
            value
        });
        let mut first = T::zero();
        let mut second = T::zero();
        for i in 0..3 {
            first = first + T::from_f64(self.omega[i]) * h[4 * i];
            for j in 0..3 {
                second = second + T::from_f64(self.w[i][j]) * h[4 * i] * h[4 * j];
                if i != j {
                    second = second + T::from_f64(self.dd[i][j]) * h[3 * i + j] * h[3 * j + i];
                }
            }
        }
        T::from_f64(self.pi0) + first + T::from_f64(0.5) * second
    }
}

// Gradient and Hessian entries of the potential in the three eigenvalues d
// of the trial tensor, from a local solution at those eigenvalues. The
// potential is extended off det = 1 through the volume-preserving
// projection of (1/2) ln d, which no admissible direction can sample.
fn principal_derivatives(local: &Local, d: &[f64; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let omega = core::array::from_fn(|i| local.tau[i] / (2.0 * d[i]));
    let w = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let mut v = local.k[i][j] / (4.0 * d[i] * d[j]);
            if i == j {
                v -= local.tau[i] / (2.0 * d[i] * d[i]);
            }
            v
        })
    });
    (omega, w)
}

struct BranchPrimal {
    lift: Lift,
    cv: [f64; 9],
    cv_new: [f64; 6],
    psi_tr: f64,
    psi_star: f64,
    dissipation: f64,
}

struct Prepared {
    branches: Vec<BranchPrimal>,
}

impl Prepared {
    fn potential<T: Real>(&self, f: &[T; 9]) -> T {
        self.branches.iter().fold(T::zero(), |s, b| {
            s + b.lift.evaluate(&trial_stretch_tensor(f, &b.cv))
        })
    }
}

impl ViscoHyperelastic {
    /// Construct with parameter checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on more than
    /// [`MAX_FINITE_STRAIN_BRANCHES`] branches, a non-finite or non-positive
    /// relaxation time, or a non-positive branch modulus.
    pub fn new(
        equilibrium: Hyperelastic,
        branches: Vec<FiniteStrainBranch>,
    ) -> Result<Self, MaterialError> {
        if branches.len() > MAX_FINITE_STRAIN_BRANCHES {
            return Err(parameters(format!(
                "{} branches exceed the cap of {MAX_FINITE_STRAIN_BRANCHES}",
                branches.len()
            )));
        }
        for b in &branches {
            if !(b.tau_s.is_finite() && b.tau_s > 0.0) {
                return Err(parameters(format!(
                    "relaxation time {} s must be finite and positive",
                    b.tau_s
                )));
            }
            b.energy.validate()?;
        }
        Ok(Self {
            equilibrium,
            branches,
            newton_tolerance: DEFAULT_LOCAL_NEWTON_TOLERANCE,
            newton_max_iterations: DEFAULT_LOCAL_NEWTON_ITERATIONS,
        })
    }

    /// Replace the local Newton tolerance (relative residual) and cap.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on a non-finite or non-positive
    /// tolerance or a zero cap.
    pub fn with_local_newton(
        mut self,
        tolerance: f64,
        max_iterations: usize,
    ) -> Result<Self, MaterialError> {
        if !(tolerance.is_finite() && tolerance > 0.0) || max_iterations == 0 {
            return Err(parameters(format!(
                "local Newton tolerance {tolerance} and cap {max_iterations} must be positive"
            )));
        }
        self.newton_tolerance = tolerance;
        self.newton_max_iterations = max_iterations;
        Ok(self)
    }

    /// The local Newton relative-residual tolerance.
    #[must_use]
    pub fn local_newton_tolerance(&self) -> f64 {
        self.newton_tolerance
    }

    /// The local Newton iteration cap.
    #[must_use]
    pub fn local_newton_max_iterations(&self) -> usize {
        self.newton_max_iterations
    }

    /// The unstrained, relaxed state: `C_v^{-1} = I`, `F = I`, empty ledger.
    #[must_use]
    pub fn initial_state(&self) -> ViscoState {
        ViscoState {
            inverse_viscous_cauchy_green: vec![[1.0, 1.0, 1.0, 0.0, 0.0, 0.0]; self.branches.len()],
            deformation_gradient: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            work: 0.0,
            stored: 0.0,
            dissipated: 0.0,
        }
    }

    fn check_inputs(&self, s: &ViscoState, f: &[f64; 9], dt: f64) -> Result<(), MaterialError> {
        if !(dt.is_finite() && dt >= 0.0) {
            return Err(state_error(format!(
                "step {dt} s must be finite and non-negative"
            )));
        }
        let j = det3(f);
        if !f.iter().all(|x| x.is_finite()) || !(j.is_finite() && j > 0.0) {
            return Err(state_error(format!(
                "det F = {j} (orientation lost or degenerate)"
            )));
        }
        if s.inverse_viscous_cauchy_green.len() != self.branches.len() {
            return Err(state_error(format!(
                "state holds {} branches, the law has {}",
                s.inverse_viscous_cauchy_green.len(),
                self.branches.len()
            )));
        }
        for cv in &s.inverse_viscous_cauchy_green {
            let m = sym_expand(cv);
            let dm = det3(&m);
            if !cv.iter().all(|x| x.is_finite()) || !(dm.is_finite() && (dm - 1.0).abs() < 1e-6) {
                return Err(state_error(format!(
                    "inverse viscous Cauchy-Green tensor has determinant {dm}, not 1"
                )));
            }
        }
        Ok(())
    }

    // Eigenvalues and eigenvectors of the trial tensor of one branch.
    fn trial_spectrum(
        cv: &[f64; 9],
        f: &[f64; 9],
    ) -> Result<([f64; 9], [f64; 3], [f64; 9]), MaterialError> {
        let ct = trial_stretch_tensor(f, cv);
        let sym: [f64; 9] =
            core::array::from_fn(|ij| 0.5 * (ct[ij] + ct[3 * (ij % 3) + ij / 3]));
        let (d, v) = jacobi_eigh(&sym, 3);
        if d.len() != 3 || d.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(state_error(
                "non-positive or non-finite trial stretch eigenvalue".to_string(),
            ));
        }
        let vv: [f64; 9] = core::array::from_fn(|i| v[i]);
        Ok((sym, [d[0], d[1], d[2]], vv))
    }

    fn local_at(
        &self,
        branch: &FiniteStrainBranch,
        dt: f64,
        d: &[f64; 3],
    ) -> Result<Local, MaterialError> {
        let eta_tr = reduce([
            0.5 * det::ln(d[0]),
            0.5 * det::ln(d[1]),
            0.5 * det::ln(d[2]),
        ]);
        local_solve(
            &branch.energy,
            branch.tau_s,
            dt,
            eta_tr,
            self.newton_tolerance,
            self.newton_max_iterations,
        )
    }

    fn prepare_branch(
        &self,
        branch: &FiniteStrainBranch,
        cv6: &[f64; 6],
        f: &[f64; 9],
        dt: f64,
    ) -> Result<BranchPrimal, MaterialError> {
        let cv = sym_expand(cv6);
        let (sym, d, v) = Self::trial_spectrum(&cv, f)?;
        let local = self.local_at(branch, dt, &d)?;
        let (omega, w) = principal_derivatives(&local, &d);
        let mut dd = [[0.0; 3]; 3];
        for (i, j) in [(0usize, 1usize), (0, 2), (1, 2)] {
            let u = d[i] - d[j];
            let big = d[i].max(d[j]);
            let value = if u.abs() > COALESCENCE_GAP * big {
                (omega[i] - omega[j]) / u
            } else if u == 0.0 {
                w[i][i] - w[i][j]
            } else {
                // (omega_i - omega_j) / (d_i - d_j) is the integral over the
                // swap path of (omega_ii - omega_ij); a smooth integrand.
                let mut sum = 0.0;
                for (xi, wt) in GL5_XI.iter().zip(GL5_W.iter()) {
                    let s = 0.5 * (1.0 + xi);
                    let mut point = d;
                    point[i] = d[j] + s * u;
                    point[j] = d[i] - s * u;
                    let at = self.local_at(branch, dt, &point)?;
                    let (_, wp) = principal_derivatives(&at, &point);
                    sum += 0.5 * wt * (wp[i][i] - wp[i][j]);
                }
                sum
            };
            dd[i][j] = value;
            dd[j][i] = value;
        }
        // State update: relaxed principal stretches in the trial eigenbasis.
        let cv_new = if dt == 0.0 {
            *cv6
        } else {
            let x_star = expand(local.eta);
            let j = det3(f);
            let jm23 = det::exp(-2.0 / 3.0 * det::ln(j));
            let mut bt = [0.0; 9];
            for a in 0..3 {
                let e = det::exp(2.0 * x_star[a]);
                for r in 0..3 {
                    for c in 0..3 {
                        bt[3 * r + c] += e * v[3 * r + a] * v[3 * c + a];
                    }
                }
            }
            let be: [f64; 9] = core::array::from_fn(|ij| bt[ij] / jm23);
            let finv = inverse3(f).ok_or_else(|| state_error("singular F".to_string()))?;
            let cv_full = mat_mul(&mat_mul(&finv, &be), &transpose(&finv));
            let packed = sym_pack(&cv_full);
            let dm = det3(&sym_expand(&packed));
            if !(dm.is_finite() && dm > 0.0) {
                return Err(state_error(
                    "relaxed viscous tensor lost positive definiteness".to_string(),
                ));
            }
            // det F_v = 1 exactly, up to the last bit of this rescaling.
            let scale = det::exp(-det::ln(dm) / 3.0);
            core::array::from_fn(|i| packed[i] * scale)
        };
        Ok(BranchPrimal {
            lift: Lift {
                v,
                ct: sym,
                pi0: local.pi,
                omega,
                w,
                dd,
            },
            cv,
            cv_new,
            psi_tr: local.psi_tr,
            psi_star: local.psi_star,
            dissipation: local.dissipation,
        })
    }

    fn prepare(
        &self,
        s: &ViscoState,
        f: &[f64; 9],
        dt: f64,
    ) -> Result<Prepared, MaterialError> {
        self.check_inputs(s, f, dt)?;
        let mut branches = Vec::with_capacity(self.branches.len());
        for (branch, cv6) in self.branches.iter().zip(&s.inverse_viscous_cauchy_green) {
            branches.push(self.prepare_branch(branch, cv6, f, dt)?);
        }
        Ok(Prepared { branches })
    }

    /// Advance one step to `f_new` over `dt` seconds. Returns the new state
    /// and the first Piola-Kirchhoff stress at `f_new` (row-major).
    /// `dt = 0` is an elastic jump (equilibrium plus the instantaneous
    /// response of every branch) with no viscous loss.
    ///
    /// # Errors
    /// [`MaterialError::State`] on a negative or non-finite `dt`, det F <= 0,
    /// a state that does not match the law, an equilibrium-law domain
    /// refusal, or a local Newton that does not converge.
    pub fn step(
        &self,
        s: &ViscoState,
        f_new: &[f64; 9],
        dt: f64,
    ) -> Result<(ViscoState, [f64; 9]), MaterialError> {
        let prep = self.prepare(s, f_new, dt)?;
        let p = self.stress_from(&prep, f_new)?;
        let psi_eq = self.equilibrium.energy(f_new);
        let psi_tr: f64 = prep.branches.iter().map(|b| b.psi_tr).sum();
        let psi_star: f64 = prep.branches.iter().map(|b| b.psi_star).sum();
        let dissipation: f64 = prep.branches.iter().map(|b| b.dissipation).sum();
        let next = ViscoState {
            inverse_viscous_cauchy_green: prep.branches.iter().map(|b| b.cv_new).collect(),
            deformation_gradient: *f_new,
            work: s.work + (psi_eq + psi_tr - s.stored),
            stored: psi_eq + psi_star,
            dissipated: s.dissipated + dissipation,
        };
        Ok((next, p))
    }

    fn stress_from(&self, prep: &Prepared, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        let mut p = self.equilibrium.piola(f)?;
        let (_, grad) = fs_ad::gradient::<9>(*f, |x| prep.potential(&x));
        for (a, g) in p.iter_mut().zip(grad) {
            *a += g;
        }
        if p.iter().all(|x| x.is_finite()) {
            Ok(p)
        } else {
            Err(state_error("nonfinite Piola stress".to_string()))
        }
    }

    /// The algorithmic tangent `dP/dF` of [`Self::step`] at fixed state and
    /// step (9x9, symmetric): the exact Hessian of the incremental potential.
    ///
    /// # Errors
    /// As [`Self::step`].
    pub fn algorithmic_tangent(
        &self,
        s: &ViscoState,
        f_new: &[f64; 9],
        dt: f64,
    ) -> Result<[[f64; 9]; 9], MaterialError> {
        let prep = self.prepare(s, f_new, dt)?;
        let mut a = self.equilibrium.tangent(f_new)?;
        for col in 0..9 {
            let mut vars = [Dual::<Dual64<9>, 1>::constant(Dual64::constant(0.0)); 9];
            for (i, v) in vars.iter_mut().enumerate() {
                *v = Dual {
                    re: Dual64::variable(f_new[i], i),
                    eps: [if i == col {
                        Dual64::constant(1.0)
                    } else {
                        Dual64::constant(0.0)
                    }],
                };
            }
            let w = prep.potential(&vars);
            for (row, out) in a.iter_mut().enumerate() {
                out[col] += w.eps[0].eps[row];
            }
        }
        if a.iter().flatten().all(|x| x.is_finite()) {
            Ok(a)
        } else {
            Err(state_error("nonfinite algorithmic tangent".to_string()))
        }
    }

    /// Stored energy density at `f` with the internal state frozen:
    /// equilibrium plus branch energies (J/m^3).
    ///
    /// # Errors
    /// As [`Self::step`] (with `dt = 0`).
    pub fn free_energy(&self, s: &ViscoState, f: &[f64; 9]) -> Result<f64, MaterialError> {
        self.check_inputs(s, f, 0.0)?;
        let mut total = self.equilibrium.energy(f);
        for (branch, cv6) in self.branches.iter().zip(&s.inverse_viscous_cauchy_green) {
            let cv = sym_expand(cv6);
            let (_, d, _) = Self::trial_spectrum(&cv, f)?;
            let x = expand(reduce([
                0.5 * det::ln(d[0]),
                0.5 * det::ln(d[1]),
                0.5 * det::ln(d[2]),
            ]));
            total += branch.energy.terms(&x).0;
        }
        if total.is_finite() {
            Ok(total)
        } else {
            Err(state_error("nonfinite stored energy".to_string()))
        }
    }

    /// Thermodynamic declarations.
    #[must_use]
    pub fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "validity ends at the equilibrium law's stretch_limit; branch energies carry no separate limit",
        }
    }
}
