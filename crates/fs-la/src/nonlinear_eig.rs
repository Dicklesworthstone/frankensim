//! Isolated finite poles of an analytic dense matrix function T(s).
//!
//! The single numerical owner uses Newton on the vector equation and one
//! normalization constraint, with Jacobian [T, T'v; e_p^T, 0]. This includes
//! nonsymmetric functions and regular descriptor pencils without inverting
//! their singular mass matrix. The pivot p is the largest vector component.
//! See Guettel and Tisseur (2017), Acta Numerica 26, section 4.2,
//! doi:10.1017/S0962492917000034.
//!
//! Both the original-equation residual and the scaled pole correction must
//! pass. A small Rayleigh quotient or a zero pole correction alone is not an
//! eigenpair. All matrices, scaling, budgets and cancellation are explicit.
//! This is local refinement of one simple finite pole from a supplied seed,
//! not a spectrum count, eigenvalue enclosure or guarantee of convergence
//! from an arbitrary start. Multiple/defective poles may refuse.

use crate::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, det};

/// Hard ceiling for this dense local method, before the normalization row.
pub const MAX_NONLINEAR_POLE_DIMENSION: usize = 1024;

/// Typed refusals of the isolated-pole owner.
#[derive(Debug, Clone, PartialEq)]
pub enum NonlinearEigError {
    /// Invalid dimensions, finite data, normalization, tolerances or capacity.
    BadInput { what: &'static str },
    /// The caller cannot evaluate its analytic operator at this argument.
    FunctionRefused { at: C64 },
    /// T'(s)v vanishes, so the local pole correction is undefined.
    DegenerateFunctional { iteration: usize },
    /// The normalized Newton Jacobian is singular; no isolated pair is claimed.
    BorderedSingular { iteration: usize, index: usize },
    /// The requested refinement did not satisfy BOTH convergence conditions.
    NotConverged {
        iterations: usize,
        last_relative_step: f64,
    },
    /// Cooperative cancellation; no partial eigenpair is published.
    Cancelled,
}
impl core::fmt::Display for NonlinearEigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LA-NONLINEAR-EIG: {self:?}")
    }
}
impl std::error::Error for NonlinearEigError {}

/// Explicit accuracy and resources for one isolated-pole refinement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PoleOptions {
    /// Positive reference scale in the same units as s; also resolves zero poles.
    pub spectral_scale: f64,
    /// Bound on |delta s| / max(spectral_scale, |s|).
    pub step_tolerance: f64,
    /// Bound on ||T(s)v||_2 / (evaluation.residual_scale * ||v||_2).
    pub residual_tolerance: f64,
    /// Maximum number of Newton evaluations/factorizations.
    pub max_iterations: usize,
    /// Allowed original dimension, at most MAX_NONLINEAR_POLE_DIMENSION.
    pub max_dimension: usize,
    /// Conservative internal numerical payload, excluding callback-owned caches.
    pub max_dense_bytes: usize,
    /// Bound max_iterations * (n+1)^2 * (n+17), before calling the operator.
    pub max_work: u64,
}

/// Original operator, analytic derivative, and its physical normalization.
#[derive(Debug, Clone)]
pub struct PoleEvaluation {
    /// Row-major T(s), n by n.
    pub matrix: Vec<C64>,
    /// Row-major analytic T'(s), n by n.
    pub derivative: Vec<C64>,
    /// Positive finite operator magnitude in the units of T, supplied by the
    /// producer. For a polynomial, use a coefficient envelope such as
    /// ||K|| + |s|||C|| + |s|^2||M||, with a declared nonzero reference scale
    /// when that envelope vanishes. Using ||T(s)|| alone is unsuitable for a
    /// scalar root, where numerator and denominator vanish together.
    pub residual_scale: f64,
}

/// An accepted eigenpair in the ORIGINAL coordinates.
#[derive(Debug, Clone)]
pub struct NonlinearPole {
    pub value: C64,
    /// Unit Euclidean norm, with largest component real and positive.
    pub vector: Vec<C64>,
    /// Original equation residual relative to the producer's explicit scale.
    pub residual: f64,
    /// Scaled Newton correction at the accepted pair.
    pub relative_step: f64,
    /// Operator scale used for the reported residual.
    pub residual_scale: f64,
    pub iterations: usize,
}

fn bad(what: &'static str) -> NonlinearEigError {
    NonlinearEigError::BadInput { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), NonlinearEigError> {
    if cancelled() {
        Err(NonlinearEigError::Cancelled)
    } else {
        Ok(())
    }
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn norm2(x: &[C64]) -> f64 {
    let scale = x.iter().map(|z| z.abs()).fold(0.0, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    if !scale.is_finite() {
        return f64::INFINITY;
    }
    scale * det::sqrt(x.iter().map(|&z| divide(z, scale).norm_sq()).sum())
}

// Set v_p=1 without squaring an unscaled coefficient or losing its phase.
// Repeating this after each Newton step keeps every component bounded by one.
fn normalize(v: &mut [C64]) -> Result<usize, NonlinearEigError> {
    if v.is_empty() || v.iter().any(|&z| !finite(z)) {
        return Err(bad("eigenvector must be finite and nonempty"));
    }
    let mut pivot = 0;
    for i in 1..v.len() {
        if v[i].abs() > v[pivot].abs() {
            pivot = i;
        }
    }
    let scale = v[pivot].abs();
    if scale == 0.0 || !scale.is_finite() {
        return Err(bad("eigenvector must have a representable nonzero scale"));
    }
    for value in v.iter_mut() {
        *value = divide(*value, scale);
    }
    let phase = v[pivot];
    for value in v.iter_mut() {
        *value = *value / phase;
    }
    v[pivot] = C64::ONE;
    Ok(pivot)
}

fn admit(n: usize, options: &PoleOptions) -> Result<usize, NonlinearEigError> {
    let d = n
        .checked_add(1)
        .ok_or_else(|| bad("bordered dimension overflow"))?;
    let square = d
        .checked_mul(d)
        .ok_or_else(|| bad("bordered matrix extent overflow"))?;
    let bytes = square
        .checked_mul(128)
        .and_then(|v| d.checked_mul(256).and_then(|b| v.checked_add(b)))
        .ok_or_else(|| bad("pole payload bound overflow"))?;
    let work = (square as u64)
        .checked_mul(d as u64 + 16)
        .and_then(|v| v.checked_mul(options.max_iterations as u64))
        .ok_or_else(|| bad("pole work bound overflow"))?;
    if n == 0
        || n > MAX_NONLINEAR_POLE_DIMENSION
        || n > options.max_dimension
        || options.max_iterations == 0
        || bytes > options.max_dense_bytes
        || work > options.max_work
        || !options.spectral_scale.is_finite()
        || options.spectral_scale <= 0.0
        || !options.step_tolerance.is_finite()
        || options.step_tolerance <= 0.0
        || options.step_tolerance >= 1.0
        || !options.residual_tolerance.is_finite()
        || options.residual_tolerance <= 0.0
        || options.residual_tolerance >= 1.0
    {
        return Err(bad(
            "invalid pole accuracy or insufficient dimension, byte or work budget",
        ));
    }
    Ok(d)
}

fn lu_error(error: LuComplexError, iteration: usize) -> NonlinearEigError {
    match error {
        LuComplexError::Cancelled => NonlinearEigError::Cancelled,
        LuComplexError::Singular { index } => {
            NonlinearEigError::BorderedSingular { iteration, index }
        }
    }
}

/// Refine one isolated finite pole of a general analytic matrix function.
///
/// The callback fills the ORIGINAL T and analytic T', plus a declared physical
/// coefficient scale; it can poll the supplied cancellation function while
/// assembling. No complex symmetry, Hermiticity, positive mass or invertible
/// descriptor matrix is assumed. The callback owns its analytic domain and
/// coefficient validity. It must not silently drop states or change branches.
///
/// Fixed-order row/column equilibration affects only the Newton solve. The
/// acceptance residual is evaluated before equilibration on the original
/// equations. Exact zero-frequency poles use the explicit spectral scale.
/// A singular normalized Jacobian refuses, including unresolved multiple poles.
pub fn newton_pole_with_cancel<F, C>(
    n: usize,
    s0: C64,
    v0: &[C64],
    options: &PoleOptions,
    mut evaluate: F,
    mut cancelled: C,
) -> Result<NonlinearPole, NonlinearEigError>
where
    F: FnMut(C64, &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, NonlinearEigError>,
    C: FnMut() -> bool,
{
    poll(&mut cancelled)?;
    let d = admit(n, options)?;
    if !finite(s0) || v0.len() != n {
        return Err(bad(
            "finite start value and exactly n start coordinates required",
        ));
    }
    let mut v = v0.to_vec();
    let mut pivot = normalize(&mut v)?;
    let mut s = s0;
    let mut last = f64::INFINITY;
    for iteration in 1..=options.max_iterations {
        poll(&mut cancelled)?;
        let evaluation = evaluate(s, &mut cancelled)?;
        poll(&mut cancelled)?;
        let PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        } = evaluation;
        if matrix.len() != n * n
            || derivative.len() != n * n
            || matrix.iter().chain(&derivative).any(|&z| !finite(z))
            || !residual_scale.is_finite()
            || residual_scale <= 0.0
        {
            return Err(bad(
                "T, T' and the physical residual scale must be finite and correctly sized",
            ));
        }
        let mut border = vec![C64::ZERO; d * d];
        let mut rhs = vec![C64::ZERO; d];
        let mut derivative_norm = 0.0_f64;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut r = C64::ZERO;
            let mut dv = C64::ZERO;
            for j in 0..n {
                let a = divide(matrix[i * n + j], residual_scale);
                let b = divide(derivative[i * n + j], residual_scale).scale(options.spectral_scale);
                border[i * d + j] = a;
                r = r + a * v[j];
                dv = dv + b * v[j];
            }
            rhs[i] = -r;
            border[i * d + n] = dv;
            derivative_norm = derivative_norm.max(dv.abs());
        }
        border[n * d + pivot] = C64::ONE;
        if border.iter().chain(&rhs).any(|&z| !finite(z)) {
            return Err(bad("normalized pole equations are not representable"));
        }
        if derivative_norm == 0.0 {
            return Err(NonlinearEigError::DegenerateFunctional { iteration });
        }
        let relative = norm2(&rhs[..n]) / norm2(&v);
        if !relative.is_finite() {
            return Err(bad("original pole residual is not representable"));
        }

        for i in 0..d {
            poll(&mut cancelled)?;
            let scale = border[i * d..(i + 1) * d]
                .iter()
                .map(|z| z.abs())
                .fold(0.0, f64::max);
            if scale == 0.0 || !scale.is_finite() {
                return Err(NonlinearEigError::BorderedSingular {
                    iteration,
                    index: i,
                });
            }
            for a in &mut border[i * d..(i + 1) * d] {
                *a = divide(*a, scale);
            }
            rhs[i] = divide(rhs[i], scale);
        }
        let mut columns = vec![0.0; d];
        for j in 0..d {
            poll(&mut cancelled)?;
            let scale = (0..d).map(|i| border[i * d + j].abs()).fold(0.0, f64::max);
            if scale == 0.0 || !scale.is_finite() {
                return Err(NonlinearEigError::BorderedSingular {
                    iteration,
                    index: j,
                });
            }
            columns[j] = scale;
            for i in 0..d {
                border[i * d + j] = divide(border[i * d + j], scale);
            }
        }
        let lu = lu_complex_with_cancel(&border, d, &mut cancelled)
            .map_err(|error| lu_error(error, iteration))?;
        lu.solve_with_cancel(&mut rhs, &mut cancelled)
            .map_err(|error| lu_error(error, iteration))?;
        for (delta, scale) in rhs.iter_mut().zip(columns) {
            *delta = divide(*delta, scale);
        }
        if rhs.iter().any(|&z| !finite(z)) {
            return Err(bad("Newton correction is not representable"));
        }
        let step = rhs[n].scale(options.spectral_scale);
        if !finite(step) {
            return Err(bad("pole correction overflows"));
        }
        last = step.abs() / options.spectral_scale.max(s.abs());
        if relative <= options.residual_tolerance && last <= options.step_tolerance {
            let length = norm2(&v);
            for value in &mut v {
                *value = divide(*value, length);
            }
            poll(&mut cancelled)?;
            return Ok(NonlinearPole {
                value: s,
                vector: v,
                residual: relative,
                relative_step: last,
                residual_scale,
                iterations: iteration,
            });
        }
        s = s + step;
        if !finite(s) {
            return Err(bad("updated pole is not representable"));
        }
        for (value, &delta) in v.iter_mut().zip(&rhs[..n]) {
            *value = *value + delta;
        }
        pivot = normalize(&mut v)?;
    }
    Err(NonlinearEigError::NotConverged {
        iterations: options.max_iterations,
        last_relative_step: last,
    })
}

/// Compatibility entry point for the original complex-symmetric local solver.
///
/// This now uses the same residual-gated bordered Newton owner as the general
/// API. The unchanged tuple callback supplies T and T'. Its normalization is
/// fixed from the initial coefficient norms, max(||T(s0)||_F,
/// max(1,|s0|)||T'(s0)||_F), so scalar and zero roots have a meaningful residual.
/// New physical callers should supply their own coefficient envelope and
/// resources through newton_pole_with_cancel. The dense dimension cap applies.
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
    let spectral_scale = s0.abs().max(1.0);
    let options = PoleOptions {
        spectral_scale,
        step_tolerance: tol,
        residual_tolerance: tol,
        max_iterations: max_iter,
        max_dimension: n,
        max_dense_bytes: usize::MAX,
        max_work: u64::MAX,
    };
    let mut fixed_scale = None;
    newton_pole_with_cancel(
        n,
        s0,
        v0,
        &options,
        |s, _| {
            let (matrix, derivative) = t(s).ok_or(NonlinearEigError::FunctionRefused { at: s })?;
            let scale = *fixed_scale
                .get_or_insert_with(|| norm2(&matrix).max(spectral_scale * norm2(&derivative)));
            if scale == 0.0 {
                return Err(NonlinearEigError::DegenerateFunctional { iteration: 1 });
            }
            Ok(PoleEvaluation {
                matrix,
                derivative,
                residual_scale: scale,
            })
        },
        || false,
    )
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

#[cfg(test)]
mod isolated_pole_regressions {
    use super::*;
    use std::cell::Cell;

    // Keep the regression executable compact while exercising the same public owner.
    fn solve<F, C>(
        n: usize,
        s: C64,
        v: &[C64],
        options: &PoleOptions,
        mut evaluate: F,
        mut cancelled: C,
    ) -> Result<NonlinearPole, NonlinearEigError>
    where
        F: FnMut(C64, &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, NonlinearEigError>,
        C: FnMut() -> bool,
    {
        let evaluate: &mut dyn FnMut(
            C64,
            &mut dyn FnMut() -> bool,
        ) -> Result<PoleEvaluation, NonlinearEigError> = &mut evaluate;
        let cancelled: &mut dyn FnMut() -> bool = &mut cancelled;
        newton_pole_with_cancel(n, s, v, options, evaluate, cancelled)
    }

    fn options() -> PoleOptions {
        PoleOptions {
            spectral_scale: 2.0,
            step_tolerance: 1e-12,
            residual_tolerance: 1e-12,
            max_iterations: 24,
            max_dimension: 8,
            max_dense_bytes: 1 << 20,
            max_work: 1_000_000,
        }
    }

    // E=diag(1,0), A=[[-2,1],[3,-4]]. The algebraic row gives
    // v1=3*v0/4 and the only finite pencil pole is -5/4.
    fn descriptor(
        s: C64,
        _: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        Ok(PoleEvaluation {
            matrix: vec![
                s + C64::from_re(2.0),
                C64::from_re(-1.0),
                C64::from_re(-3.0),
                C64::from_re(4.0),
            ],
            derivative: vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO],
            residual_scale: 7.0 + s.abs(),
        })
    }

    #[test]
    fn stationary_rayleigh_value_is_refined_instead_of_published() {
        // At s=0, v=(1,1), v^T T(s)v vanishes but T(s)v does not.
        // The former Rayleigh-step acceptance returned this non-eigenpair.
        let pole = rayleigh_newton_pole(2, C64::ZERO, &[C64::ONE, C64::ONE], 1e-12, 12, |s| {
            Some((
                vec![s - C64::ONE, C64::ZERO, C64::ZERO, s + C64::ONE],
                vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
            ))
        })
        .expect("the bordered equation resolves the stationary Rayleigh value");
        let r0 = (pole.value - C64::ONE) * pole.vector[0];
        let r1 = (pole.value + C64::ONE) * pole.vector[1];
        assert!(r0.abs().hypot(r1.abs()) < 1e-11, "{pole:?}");
        assert!(
            (pole.value - C64::ONE)
                .abs()
                .min((pole.value + C64::ONE).abs())
                < 1e-11
        );
        assert!(pole.residual <= 1e-12);
    }

    #[test]
    fn simple_zero_pole_uses_a_nonzero_reference_scale() {
        let pole = rayleigh_newton_pole(
            1,
            C64::new(0.1, 0.2),
            &[C64::new(3.0, 4.0)],
            1e-12,
            12,
            |s| Some((vec![s], vec![C64::ONE])),
        )
        .expect("a simple zero pole is a finite eigenvalue");
        assert!(pole.value.abs() < 1e-12);
        assert!(pole.residual.is_finite() && pole.residual <= 1e-12);
        assert!(pole.relative_step <= 1e-12);
        assert!((pole.vector[0].abs() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn scalar_pole_preserves_extreme_coefficient_and_vector_scaling() {
        let exact = C64::new(-0.125, 3.0);
        for scale in [1e-250, 1.0, 1e250] {
            let pole = rayleigh_newton_pole(
                1,
                C64::new(-0.2, 3.3),
                &[C64::new(0.6 * scale, 0.8 * scale)],
                1e-12,
                12,
                |s| Some((vec![(s - exact).scale(scale)], vec![C64::from_re(scale)])),
            )
            .expect("finite changes of coefficient and eigenvector scale preserve the pole");
            assert!(
                (pole.value - exact).abs() < 1e-11,
                "scale={scale:e}: {pole:?}"
            );
            assert!(pole.residual.is_finite() && pole.residual <= 1e-12);
            assert!((pole.vector[0].abs() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn nonsymmetric_singular_mass_pencil_retains_its_algebraic_row() {
        let pole = solve(
            2,
            C64::from_re(-1.1),
            &[C64::ONE, C64::from_re(0.5)],
            &options(),
            descriptor,
            || false,
        )
        .expect("a regular descriptor pencil has an isolated finite pole");
        assert!((pole.value - C64::from_re(-1.25)).abs() < 1e-12, "{pole:?}");
        let dynamic = (pole.value + C64::from_re(2.0)) * pole.vector[0] - pole.vector[1];
        let algebraic = pole.vector[0].scale(-3.0) + pole.vector[1].scale(4.0);
        assert!(dynamic.abs() < 1e-12 && algebraic.abs() < 1e-12);
        assert!(pole.residual <= options().residual_tolerance);
        assert!(pole.relative_step <= options().step_tolerance);
    }

    #[test]
    fn cancellation_never_returns_a_partial_or_cached_pair() {
        let polls = Cell::new(0);
        let reference = solve(
            2,
            C64::from_re(-1.1),
            &[C64::ONE, C64::from_re(0.5)],
            &options(),
            descriptor,
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
        let total = polls.get();
        assert!(total > 12, "fixture must enter the factorization and solve");
        for cut in [1, 2, 4, 12, total / 2, total] {
            let calls = Cell::new(0);
            let result = solve(
                2,
                C64::from_re(-1.1),
                &[C64::ONE, C64::from_re(0.5)],
                &options(),
                descriptor,
                || {
                    calls.set(calls.get() + 1);
                    calls.get() >= cut
                },
            );
            assert_eq!(
                result.unwrap_err(),
                NonlinearEigError::Cancelled,
                "poll {cut}"
            );
        }
        let replay = solve(
            2,
            C64::from_re(-1.1),
            &[C64::ONE, C64::from_re(0.5)],
            &options(),
            descriptor,
            || false,
        )
        .unwrap();
        assert_eq!(reference.value, replay.value);
        assert_eq!(reference.vector, replay.vector);
        assert_eq!(reference.residual.to_bits(), replay.residual.to_bits());
    }

    #[test]
    fn operator_callback_can_poll_the_original_cancellation_source() {
        let requested = Cell::new(false);
        let result = solve(
            2,
            C64::from_re(-1.1),
            &[C64::ONE, C64::from_re(0.5)],
            &options(),
            |_s, cancelled| {
                requested.set(true);
                assert!(cancelled());
                Err(NonlinearEigError::Cancelled)
            },
            || requested.get(),
        );
        assert_eq!(result.unwrap_err(), NonlinearEigError::Cancelled);
    }

    #[test]
    fn work_memory_and_dimension_refuse_before_operator_evaluation() {
        let mut limits = [options(); 4];
        limits[0].max_dense_bytes = 0;
        limits[1].max_work = 0;
        limits[2].max_dimension = 1;
        limits[3].max_iterations = 0;
        for settings in limits {
            let result = solve(
                2,
                C64::from_re(-1.1),
                &[C64::ONE, C64::from_re(0.5)],
                &settings,
                |_, _| panic!("inadmissible work must not call the producer"),
                || false,
            );
            assert!(matches!(result, Err(NonlinearEigError::BadInput { .. })));
        }
        let result = solve(
            usize::MAX,
            C64::ZERO,
            &[],
            &options(),
            |_, _| panic!("overflowing shape must not call the producer"),
            || false,
        );
        assert!(matches!(result, Err(NonlinearEigError::BadInput { .. })));
    }

    #[test]
    fn malformed_nonfinite_and_refused_evaluations_are_typed() {
        for defect in 0..6 {
            let result = solve(
                1,
                C64::ZERO,
                &[C64::ONE],
                &options(),
                |s, _| {
                    let mut e = PoleEvaluation {
                        matrix: vec![s + C64::ONE],
                        derivative: vec![C64::ONE],
                        residual_scale: 2.0,
                    };
                    match defect {
                        0 => e.matrix.clear(),
                        1 => e.derivative.clear(),
                        2 => e.residual_scale = 0.0,
                        3 => e.residual_scale = f64::NAN,
                        4 => e.matrix[0] = C64::new(f64::NAN, 0.0),
                        _ => e.derivative[0] = C64::new(0.0, f64::INFINITY),
                    }
                    Ok(e)
                },
                || false,
            );
            assert!(
                matches!(result, Err(NonlinearEigError::BadInput { .. })),
                "defect {defect}"
            );
        }
        let at = C64::new(-0.5, 0.25);
        let result = solve(
            1,
            at,
            &[C64::ONE],
            &options(),
            |s, _| Err(NonlinearEigError::FunctionRefused { at: s }),
            || false,
        );
        assert_eq!(
            result.unwrap_err(),
            NonlinearEigError::FunctionRefused { at }
        );
    }

    #[test]
    fn vanishing_derivative_and_unresolved_multiple_poles_refuse() {
        let derivative_zero = solve(
            1,
            C64::ONE,
            &[C64::ONE],
            &options(),
            |s, _| {
                let x = s - C64::ONE;
                Ok(PoleEvaluation {
                    matrix: vec![x * x],
                    derivative: vec![x.scale(2.0)],
                    residual_scale: 1.0,
                })
            },
            || false,
        );
        assert!(matches!(
            derivative_zero,
            Err(NonlinearEigError::DegenerateFunctional { .. })
        ));
        let repeated = solve(
            2,
            C64::from_re(-1.0),
            &[C64::ONE, C64::ZERO],
            &options(),
            |s, _| {
                Ok(PoleEvaluation {
                    matrix: vec![s + C64::ONE, C64::ZERO, C64::ZERO, s + C64::ONE],
                    derivative: vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ONE],
                    residual_scale: 2.0,
                })
            },
            || false,
        );
        assert!(matches!(
            repeated,
            Err(NonlinearEigError::BorderedSingular { .. })
        ));
    }
}
