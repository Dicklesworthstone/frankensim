//! Isolated eigenvalues of a dense nonlinear matrix function: find `s` and `v != 0` with `T(s) v = 0`, where `T` is
//! complex symmetric (`T^T = T`) and analytic near the wanted eigenvalue.
//!
//! Method: Newton's method on the Rayleigh functional combined with nonlinear inverse iteration (Ruhe 1973, SIAM J.
//! Numer. Anal. 10:674-689, section 3). From `(s_k, v_k)`:
//! `s_{k+1} = s_k - (v_k^T T(s_k) v_k) / (v_k^T T'(s_k) v_k)`, then `v_{k+1} = T(s_{k+1})^{-1} T'(s_{k+1}) v_k`,
//! normalized. The unconjugated bilinear form is the Rayleigh functional of a complex-symmetric `T`. For a simple
//! eigenvalue the iteration converges quadratically.
//!
//! The stopping test is the relative step `|s_{k+1} - s_k| <= tol |s_{k+1}|`. The returned residual
//! `||T(s) v|| / (||T(s)||_F ||v||)` is computed at the returned pair, never assumed. A pivot that vanishes exactly means
//! `s` is an eigenvalue to working precision; the iteration then keeps its last direction. Refusals: inconsistent sizes,
//! a non-finite start, a vanishing functional derivative, a caller function that refuses, or no convergence within
//! the iteration cap. Deterministic: fixed operation order, no randomness.

use crate::eigen_complex::lu_complex;
use fs_math::c64::C64;

/// Typed refusals of [`rayleigh_newton_pole`].
#[derive(Debug, Clone, PartialEq)]
pub enum NonlinearEigError {
    /// Sizes, tolerance, iteration cap or start are inconsistent or non-finite.
    BadInput {
        /// What is wrong.
        what: &'static str,
    },
    /// The caller's matrix function refused at this argument.
    FunctionRefused {
        /// The argument it refused.
        at: C64,
    },
    /// `v^T T'(s) v` vanished or was not finite, so the Newton step is undefined.
    DegenerateFunctional {
        /// Iteration at which it happened.
        iteration: usize,
    },
    /// The step did not fall below the tolerance within the cap.
    NotConverged {
        /// Iterations spent.
        iterations: usize,
        /// The last relative step.
        last_relative_step: f64,
    },
}

impl core::fmt::Display for NonlinearEigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LA-NONLINEAR-EIG: {self:?}")
    }
}

impl std::error::Error for NonlinearEigError {}

/// A converged eigenpair of `T(s) v = 0`.
#[derive(Debug, Clone)]
pub struct NonlinearPole {
    /// The eigenvalue `s`.
    pub value: C64,
    /// The eigenvector, unit 2-norm.
    pub vector: Vec<C64>,
    /// `||T(s) v|| / (||T(s)||_F ||v||)` at the returned pair.
    pub residual: f64,
    /// Newton iterations spent.
    pub iterations: usize,
}

fn bilinear(v: &[C64], a: &[C64], w: &[C64]) -> C64 {
    let n = v.len();
    let mut acc = C64::ZERO;
    for i in 0..n {
        let mut row = C64::ZERO;
        for j in 0..n {
            row = row + a[i * n + j] * w[j];
        }
        acc = acc + v[i] * row;
    }
    acc
}

fn matvec(a: &[C64], x: &[C64]) -> Vec<C64> {
    let n = x.len();
    (0..n)
        .map(|i| {
            let mut s = C64::ZERO;
            for j in 0..n {
                s = s + a[i * n + j] * x[j];
            }
            s
        })
        .collect()
}

fn norm2(x: &[C64]) -> f64 {
    x.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt()
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Newton on the Rayleigh functional with nonlinear inverse iteration, started from `(s0, v0)`.
///
/// `t(s)` returns `(T(s), T'(s))`, each row-major `n x n`, or `None` to refuse that argument.
///
/// # Errors
/// [`NonlinearEigError`] as listed in the module documentation.
pub fn rayleigh_newton_pole<F>(
    n: usize,
    s0: C64,
    v0: &[C64],
    tol: f64,
    max_iter: usize,
    mut t: F,
) -> Result<NonlinearPole, NonlinearEigError>
where
    F: FnMut(C64) -> Option<(Vec<C64>, Vec<C64>)>,
{
    if n == 0 || v0.len() != n {
        return Err(NonlinearEigError::BadInput {
            what: "start vector length must equal n > 0",
        });
    }
    if !(tol > 0.0 && tol.is_finite()) || max_iter == 0 {
        return Err(NonlinearEigError::BadInput {
            what: "tolerance must be positive and finite, iteration cap positive",
        });
    }
    let v0_norm = norm2(v0);
    if !finite(s0) || !(v0_norm > 0.0 && v0_norm.is_finite()) {
        return Err(NonlinearEigError::BadInput {
            what: "start value and vector must be finite and nonzero",
        });
    }
    let mut eval = |s: C64| -> Result<(Vec<C64>, Vec<C64>), NonlinearEigError> {
        let (a, d) = t(s).ok_or(NonlinearEigError::FunctionRefused { at: s })?;
        if a.len() != n * n || d.len() != n * n {
            return Err(NonlinearEigError::BadInput {
                what: "T(s) and T'(s) must be n x n",
            });
        }
        Ok((a, d))
    };
    let mut v: Vec<C64> = v0.iter().map(|z| z.scale(1.0 / v0_norm)).collect();
    let mut s = s0;
    let (mut tm, mut dt) = eval(s)?;
    let mut last = f64::INFINITY;
    for iteration in 1..=max_iter {
        let num = bilinear(&v, &tm, &v);
        let den = bilinear(&v, &dt, &v);
        if !(finite(den) && den.abs() > 0.0 && finite(num)) {
            return Err(NonlinearEigError::DegenerateFunctional { iteration });
        }
        let s_new = s - num / den;
        if !finite(s_new) {
            return Err(NonlinearEigError::DegenerateFunctional { iteration });
        }
        let (tm_new, dt_new) = eval(s_new)?;
        let mut u = matvec(&dt_new, &v);
        if let Ok(lu) = lu_complex(&tm_new, n) {
            lu.solve(&mut u);
            let nu = norm2(&u);
            if nu > 0.0 && nu.is_finite() {
                v = u.iter().map(|z| z.scale(1.0 / nu)).collect();
            }
        }
        last = (s_new - s).abs() / s_new.abs().max(f64::MIN_POSITIVE);
        s = s_new;
        tm = tm_new;
        dt = dt_new;
        if last <= tol {
            let r = norm2(&matvec(&tm, &v));
            let fro = tm.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt();
            let residual = if fro > 0.0 { r / fro } else { r };
            return Ok(NonlinearPole {
                value: s,
                vector: v,
                residual,
                iterations: iteration,
            });
        }
    }
    Err(NonlinearEigError::NotConverged {
        iterations: max_iter,
        last_relative_step: last,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigen_complex::eig;

    fn re(x: f64) -> C64 {
        C64::from_re(x)
    }

    /// A damped two-degree-of-freedom system `T(s) = K + s C + s^2 M` (complex symmetric): every eigenvalue found by
    /// the Newton iteration must be one of the companion-matrix eigenvalues from the dense QR solver.
    #[test]
    fn quadratic_pencil_matches_companion_eigenvalues() {
        let k = [5.0, -1.0, -1.0, 3.0];
        let c = [0.2, 0.05, 0.05, 0.1];
        let m = [2.0, 0.0, 0.0, 1.0];
        let t = |s: C64| {
            let a: Vec<C64> = (0..4)
                .map(|i| re(k[i]) + s * re(c[i]) + s * s * re(m[i]))
                .collect();
            let d: Vec<C64> = (0..4)
                .map(|i| re(c[i]) + (s * re(m[i])).scale(2.0))
                .collect();
            Some((a, d))
        };
        // Companion form of M s^2 + C s + K: [[0, I], [-M^-1 K, -M^-1 C]].
        let minv = [0.5, 0.0, 0.0, 1.0];
        let mut comp = vec![C64::ZERO; 16];
        comp[2] = C64::ONE;
        comp[4 + 3] = C64::ONE;
        for i in 0..2 {
            for j in 0..2 {
                let mut mk = 0.0;
                let mut mc = 0.0;
                for l in 0..2 {
                    mk += minv[i * 2 + l] * k[l * 2 + j];
                    mc += minv[i * 2 + l] * c[l * 2 + j];
                }
                comp[(2 + i) * 4 + j] = re(-mk);
                comp[(2 + i) * 4 + 2 + j] = re(-mc);
            }
        }
        let exact = eig(&comp, 4).expect("companion eigenvalues");
        for start in [C64::new(-0.05, 1.2), C64::new(-0.05, 2.0)] {
            let pole =
                rayleigh_newton_pole(2, start, &[C64::ONE, C64::new(0.3, 0.0)], 1e-14, 40, t)
                    .expect("converges");
            let nearest = exact
                .iter()
                .map(|z| (*z - pole.value).abs())
                .fold(f64::INFINITY, f64::min);
            assert!(
                nearest < 1e-12,
                "pole {:?} misses the companion spectrum by {nearest:e}",
                pole.value
            );
            assert!(pole.residual < 1e-13, "residual {:e}", pole.residual);
        }
    }

    /// A scalar with a square-root term, the shape of a boundary-layer law: `T(s) = s^2 - 4 + b sqrt(s)`.
    #[test]
    fn scalar_with_square_root_term_converges_to_a_root() {
        let b = C64::new(0.05, 0.05);
        let t = |s: C64| {
            let r = s.sqrt();
            Some((
                vec![s * s - re(4.0) + b * r],
                vec![s.scale(2.0) + b * r.recip().scale(0.5)],
            ))
        };
        let pole = rayleigh_newton_pole(1, re(2.0), &[C64::ONE], 1e-15, 30, t).expect("converges");
        let s = pole.value;
        let f = s * s - re(4.0) + b * s.sqrt();
        assert!(f.abs() < 1e-14 * 4.0, "|T(s)| = {:e}", f.abs());
    }

    #[test]
    fn refusals_are_typed() {
        let ok = |_s: C64| Some((vec![C64::ONE], vec![C64::ONE]));
        assert!(matches!(
            rayleigh_newton_pole(1, C64::ONE, &[], 1e-12, 5, ok),
            Err(NonlinearEigError::BadInput { .. })
        ));
        assert!(matches!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ONE], 0.0, 5, ok),
            Err(NonlinearEigError::BadInput { .. })
        ));
        assert!(matches!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ZERO], 1e-12, 5, ok),
            Err(NonlinearEigError::BadInput { .. })
        ));
        let refuse = |_s: C64| None;
        assert_eq!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ONE], 1e-12, 5, refuse).err(),
            Some(NonlinearEigError::FunctionRefused { at: C64::ONE })
        );
        let flat = |_s: C64| Some((vec![C64::ONE], vec![C64::ZERO]));
        assert_eq!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ONE], 1e-12, 5, flat).err(),
            Some(NonlinearEigError::DegenerateFunctional { iteration: 1 })
        );
        // s^2 - 2 from s = 1: one iteration cannot reach 1e-15.
        let quad = |s: C64| Some((vec![s * s - re(2.0)], vec![s.scale(2.0)]));
        assert!(matches!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ONE], 1e-15, 1, quad),
            Err(NonlinearEigError::NotConverged { iterations: 1, .. })
        ));
        let wrong = |_s: C64| Some((vec![C64::ONE; 2], vec![C64::ONE]));
        assert!(matches!(
            rayleigh_newton_pole(1, C64::ONE, &[C64::ONE], 1e-12, 5, wrong),
            Err(NonlinearEigError::BadInput { .. })
        ));
    }

    #[test]
    fn repeat_is_bitwise() {
        let t = |s: C64| {
            Some((
                vec![re(3.0) - s * s, re(0.5), re(0.5), re(1.0) - s * s],
                vec![s.scale(-2.0), C64::ZERO, C64::ZERO, s.scale(-2.0)],
            ))
        };
        let a = rayleigh_newton_pole(2, re(1.6), &[C64::ONE, C64::ZERO], 1e-14, 40, t).expect("a");
        let b = rayleigh_newton_pole(2, re(1.6), &[C64::ONE, C64::ZERO], 1e-14, 40, t).expect("b");
        assert_eq!(a.value.re.to_bits(), b.value.re.to_bits());
        assert_eq!(a.value.im.to_bits(), b.value.im.to_bits());
    }
}
