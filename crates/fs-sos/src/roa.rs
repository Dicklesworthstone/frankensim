//! Certified regions of attraction (plan §9.8's show-stealer, §15.1 step 4):
//! "stable" becomes a PROVEN SET, not an eigenvalue vibe.
//!
//! For a polynomial vector field `ẋ = f(x)` with `f(0) = 0` and a quadratic
//! Lyapunov candidate `V(x) = xᵀPx`, the sublevel set `Ω_c = {V ≤ c}` lies in
//! the region of attraction of the origin if
//!
//! ```text
//! −V̇(x) − ε‖x‖² − s(x)·(c − V(x)) = σ₀(x),   s, σ₀ SOS,  ε > 0,  P ≻ 0
//! ```
//!
//! because then `V̇ ≤ −ε‖x‖² < 0` on `Ω_c \ {0}` (S-procedure: `s ≥ 0` and
//! `c − V ≥ 0` there). `Ω_c` is a compact ellipsoid, so it is positively
//! invariant and every trajectory starting in it converges to the origin
//! (Lyapunov/LaSalle).
//!
//! `V̇ = ∇V·f` is enclosed with interval arithmetic from the exact `f64`
//! coefficients of `f` and `P`, and each candidate level `c` is proved by
//! [`crate::verify::verify`] — the bisection only ever REPORTS a level whose
//! identity was verified. `P` defaults to the solution of the Lyapunov
//! equation `AᵀP + PA = −Q` of the linearization `A = Df(0)`, and its positive
//! definiteness is itself proved by interval Cholesky.
//!
//! No-claim: the certified set is an INNER estimate; failing levels above the
//! returned one say nothing about the true region (the multiplier degree, `ε`
//! and the fixed quadratic `V` all limit tightness). Model-form validity of
//! `f` (that the polynomial model is the physics) is outside this theorem.

use crate::dense;
use crate::mpoly::{IPoly, MPoly, Monomial, lie_derivative, monomials_in_degree_range};
use crate::program::{Identity, SosProgram};
use crate::sdp::SdpSettings;
use crate::verify::{Certificate, VerifyError, interval_cholesky_pd, verify};
use fs_ivl::Interval;

/// Options for [`certify_roa`].
#[derive(Debug, Clone, PartialEq)]
pub struct RoaOptions {
    /// Lyapunov-equation weight `Q` (`n×n` row-major, symmetric PD);
    /// `None` = identity.
    pub q: Option<Vec<f64>>,
    /// Decay margin `ε`; `None` = `10⁻³·λ_min(Q)`.
    pub decay: Option<f64>,
    /// Largest level tried.
    pub level_cap: f64,
    /// Maximum bisection steps after a feasible level is bracketed.
    pub bisection_steps: usize,
    /// Relative bracket width at which bisection stops.
    pub relative_tolerance: f64,
    /// SDP controls for each level's feasibility solve.
    pub sdp: SdpSettings,
}

impl Default for RoaOptions {
    fn default() -> Self {
        RoaOptions {
            q: None,
            decay: None,
            level_cap: 1e6,
            bisection_steps: 40,
            relative_tolerance: 1e-3,
            sdp: SdpSettings::default(),
        }
    }
}

/// Why no region of attraction was certified.
#[derive(Debug, Clone, PartialEq)]
pub enum RoaError {
    /// Inconsistent dimensions or non-finite data.
    BadInput {
        /// Description.
        what: String,
    },
    /// `f(0) ≠ 0`: the origin is not an equilibrium.
    NotEquilibrium {
        /// First component with a nonzero constant term.
        component: usize,
    },
    /// The Lyapunov equation is singular (an eigenvalue pair of `A` sums to 0).
    LyapunovSolveFailed,
    /// `P` is not provably positive definite (the linearization is not
    /// Hurwitz, or the supplied `P` is not PD).
    NotPositiveDefinite,
    /// No level could be proved (the last refusal is attached).
    NoCertifiedLevel {
        /// The last verification refusal, if any solve got that far.
        last: Option<VerifyError>,
    },
}

impl core::fmt::Display for RoaError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            RoaError::BadInput { what } => write!(f, "bad ROA input: {what}"),
            RoaError::NotEquilibrium { component } => write!(
                f,
                "f({component})(0) != 0: shift coordinates so the trim state is the origin"
            ),
            RoaError::LyapunovSolveFailed => {
                write!(f, "Lyapunov equation AᵀP + PA = −Q is singular")
            }
            RoaError::NotPositiveDefinite => write!(
                f,
                "Lyapunov matrix P is not provably positive definite \
                 (linearization not Hurwitz?)"
            ),
            RoaError::NoCertifiedLevel { last } => {
                write!(f, "no sublevel set could be certified")?;
                if let Some(e) = last {
                    write!(f, " (last refusal: {e})")?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for RoaError {}

/// A PROVED inner estimate `{x : xᵀPx ≤ level}` of the region of attraction.
#[derive(Debug, Clone, PartialEq)]
pub struct RoaCertificate {
    /// Lyapunov matrix `P` (row-major).
    pub lyapunov: Vec<f64>,
    /// The proved level `c`.
    pub level: f64,
    /// The decay margin `ε` in `V̇ ≤ −ε‖x‖²` on the set.
    pub decay: f64,
    /// Semi-axes `√(c/λᵢ(P))`, ascending eigenvalue order (so descending
    /// length). Floating-point derived quantities, not enclosures.
    pub semi_axes: Vec<f64>,
    /// Volume of the ellipsoid (floating-point derived).
    pub volume: f64,
    /// Smallest level at which a proof was attempted and failed (`None` when
    /// the cap itself was proved). Says nothing about the true region.
    pub unproved_level: Option<f64>,
    /// Number of SOS feasibility solves performed.
    pub sos_solves: usize,
    /// The verified certificate at `level`.
    pub certificate: Certificate,
}

/// Solve `AᵀP + PA = −Q` for symmetric `P` (`n×n`, row-major).
#[must_use]
pub fn solve_lyapunov(a: &[f64], q: &[f64], n: usize) -> Option<Vec<f64>> {
    let nn = n * n;
    let mut k = vec![0.0; nn * nn];
    for i in 0..n {
        for j in 0..n {
            let row = i * n + j;
            for l in 0..n {
                // (AᵀP)_ij = Σ_l A_li P_lj
                k[row * nn + l * n + j] += a[l * n + i];
                // (PA)_ij = Σ_l P_il A_lj
                k[row * nn + i * n + l] += a[l * n + j];
            }
        }
    }
    let rhs: Vec<f64> = q.iter().map(|v| -v).collect();
    let mut p = dense::lu_solve(&k, nn, &rhs)?;
    dense::symmetrize(&mut p, n);
    Some(p)
}

fn quadratic_form(p: &[f64], n: usize) -> MPoly {
    let mut v = MPoly::zero(n);
    for i in 0..n {
        for j in i..n {
            let m = Monomial::var(n, i).mul(&Monomial::var(n, j));
            // 2·P_ij is exact in binary floating point.
            let c = if i == j {
                p[i * n + i]
            } else {
                2.0 * p[i * n + j]
            };
            v.add_term(m, c);
        }
    }
    v
}

fn unit_ball_volume(n: usize) -> f64 {
    // V_n = π^{n/2} / Γ(n/2 + 1), via V_n = 2π/n · V_{n−2}.
    let mut v = if n.is_multiple_of(2) { 1.0 } else { 2.0 };
    let mut k = if n.is_multiple_of(2) { 2 } else { 3 };
    while k <= n {
        v *= 2.0 * std::f64::consts::PI / k as f64;
        k += 2;
    }
    v
}

/// Certify a region of attraction of the origin with the default quadratic
/// Lyapunov function from the linearization.
///
/// # Errors
/// [`RoaError`] naming why no set was proved.
pub fn certify_roa(f: &[MPoly], opts: &RoaOptions) -> Result<RoaCertificate, RoaError> {
    let n = f.len();
    if n == 0 || f.iter().any(|fi| fi.nvars() != n) {
        return Err(RoaError::BadInput {
            what: format!("vector field must have n components in n variables (n = {n})"),
        });
    }
    // Linearization A = Df(0).
    let mut a = vec![0.0; n * n];
    for (i, fi) in f.iter().enumerate() {
        for k in 0..n {
            a[i * n + k] = fi.coeff(&Monomial::var(n, k));
        }
    }
    let q = opts.q.clone().unwrap_or_else(|| {
        let mut id = vec![0.0; n * n];
        for i in 0..n {
            id[i * n + i] = 1.0;
        }
        id
    });
    if q.len() != n * n || q.iter().any(|v| !v.is_finite()) {
        return Err(RoaError::BadInput {
            what: "Q must be a finite n×n matrix".into(),
        });
    }
    let p = solve_lyapunov(&a, &q, n).ok_or(RoaError::LyapunovSolveFailed)?;
    let decay = match opts.decay {
        Some(e) => e,
        None => 1e-3 * dense::sym_eigenvalues(&q, n)[0],
    };
    certify_roa_with(f, &p, decay, opts)
}

/// Certify a region of attraction with a caller-supplied Lyapunov matrix `P`
/// and decay margin `ε`.
///
/// # Errors
/// [`RoaError`] naming why no set was proved.
#[allow(clippy::too_many_lines)]
pub fn certify_roa_with(
    f: &[MPoly],
    p: &[f64],
    decay: f64,
    opts: &RoaOptions,
) -> Result<RoaCertificate, RoaError> {
    let n = f.len();
    if n == 0 || f.iter().any(|fi| fi.nvars() != n) || p.len() != n * n {
        return Err(RoaError::BadInput {
            what: "dimension mismatch between f and P".into(),
        });
    }
    if p.iter().any(|v| !v.is_finite()) || !(decay.is_finite() && decay > 0.0) {
        return Err(RoaError::BadInput {
            what: "P must be finite and the decay margin finite and positive".into(),
        });
    }
    if !(opts.level_cap.is_finite() && opts.level_cap > 0.0) {
        return Err(RoaError::BadInput {
            what: "level cap must be finite and positive".into(),
        });
    }
    for (i, fi) in f.iter().enumerate() {
        if fi.terms().any(|(_, c)| !c.is_finite()) {
            return Err(RoaError::BadInput {
                what: format!("f[{i}] has a non-finite coefficient"),
            });
        }
        if fi.coeff(&Monomial::one(n)) != 0.0 {
            return Err(RoaError::NotEquilibrium { component: i });
        }
    }
    // The symmetric matrix the certificate is ABOUT: upper triangle of P.
    let mut psym = p.to_vec();
    for i in 0..n {
        for j in (i + 1)..n {
            psym[j * n + i] = psym[i * n + j];
        }
    }
    let pmat: Vec<Interval> = psym.iter().map(|&v| Interval::point(v)).collect();
    if interval_cholesky_pd(&pmat, n).is_err() {
        return Err(RoaError::NotPositiveDefinite);
    }
    let v = quadratic_form(&psym, n);
    let vdot = lie_derivative(&v, f);
    // constant = −V̇ − ε‖x‖² (enclosed).
    let mut constant = IPoly::zero(n).sub(&vdot);
    for i in 0..n {
        constant.add_term(
            Monomial::var(n, i).mul(&Monomial::var(n, i)),
            Interval::point(-decay),
        );
    }
    let deg_vdot = vdot
        .terms()
        .map(|(m, _)| m.degree())
        .max()
        .unwrap_or(2)
        .max(2);
    let d0 = deg_vdot.div_ceil(2);
    let ds = d0.saturating_sub(1);
    let basis0 = monomials_in_degree_range(n, 1, d0);
    let basis_s = if ds >= 1 {
        monomials_in_degree_range(n, 1, ds)
    } else {
        Vec::new()
    };
    let mut solves = 0usize;
    let mut last_err: Option<VerifyError> = None;
    let mut check = |c: f64| -> Option<Certificate> {
        solves += 1;
        let mut prog = SosProgram::new(n);
        let sigma0 = prog.sos("sigma0", basis0.clone()).ok()?;
        let mut ident = Identity::new(constant.clone());
        if !basis_s.is_empty() {
            let s = prog.sos("s", basis_s.clone()).ok()?;
            // −s·(c − V) = (V − c)·s; V's coefficients and c are exact.
            ident = ident.term(v.sub(&MPoly::constant(n, c)), s);
        }
        ident = ident.term(MPoly::constant(n, -1.0), sigma0);
        prog.add_identity(ident).ok()?;
        let sol = prog.solve_centered(&[], &opts.sdp).ok()?;
        if sol.margin.is_none_or(|t| t <= 0.0) {
            return None;
        }
        match verify(&prog, &sol.values) {
            Ok(cert) => Some(cert),
            Err(e) => {
                last_err = Some(e);
                None
            }
        }
    };
    let cap = opts.level_cap;
    let (mut lo, mut lo_cert, mut hi): (f64, Certificate, Option<f64>);
    if let Some(cert) = check(cap) {
        lo = cap;
        lo_cert = cert;
        hi = None;
    } else {
        let mut c = cap;
        let mut found = None;
        while c > cap * 1e-14 {
            let failed = c;
            c *= 0.25;
            if let Some(cert) = check(c) {
                found = Some((c, cert, failed));
                break;
            }
        }
        let Some((c_ok, cert, failed)) = found else {
            return Err(RoaError::NoCertifiedLevel { last: last_err });
        };
        lo = c_ok;
        lo_cert = cert;
        hi = Some(failed);
        for _ in 0..opts.bisection_steps {
            let h = hi.unwrap_or(lo);
            if h - lo <= opts.relative_tolerance * lo {
                break;
            }
            let mid = 0.5 * (lo + h);
            if let Some(cert) = check(mid) {
                lo = mid;
                lo_cert = cert;
            } else {
                hi = Some(mid);
            }
        }
    }
    let lambdas = dense::sym_eigenvalues(&psym, n);
    let semi_axes: Vec<f64> = lambdas.iter().map(|l| (lo / l).sqrt()).collect();
    let det: f64 = lambdas.iter().product();
    let volume = unit_ball_volume(n) * lo.powf(n as f64 / 2.0) / det.sqrt();
    Ok(RoaCertificate {
        lyapunov: psym,
        level: lo,
        decay,
        semi_axes,
        volume,
        unproved_level: hi,
        sos_solves: solves,
        certificate: lo_cert,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lyapunov_equation_residual() {
        let a = vec![0.0, 1.0, -2.0, -0.5];
        let q = vec![1.0, 0.0, 0.0, 1.0];
        let p = solve_lyapunov(&a, &q, 2).unwrap();
        // AᵀP + PA + Q == 0
        for i in 0..2 {
            for j in 0..2 {
                let mut s = q[i * 2 + j];
                for l in 0..2 {
                    s += a[l * 2 + i] * p[l * 2 + j] + p[i * 2 + l] * a[l * 2 + j];
                }
                assert!(s.abs() < 1e-12, "residual {s}");
            }
        }
    }

    #[test]
    fn unit_ball_volumes() {
        assert!((unit_ball_volume(1) - 2.0).abs() < 1e-15);
        assert!((unit_ball_volume(2) - std::f64::consts::PI).abs() < 1e-14);
        assert!((unit_ball_volume(3) - 4.0 / 3.0 * std::f64::consts::PI).abs() < 1e-14);
    }
}
