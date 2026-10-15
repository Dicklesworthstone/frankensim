//! Fallible residuals on the existing Newton--Krylov interface.
//!
//! Source-domain errors are retained, never replaced with a zero residual.
//! Directional differences are numerical Jacobian actions, not analytic
//! tangents or adjoints. The caller owns the coordinate scale, domain and work
//! budget; the existing nonlinear owner still owns iteration/globalization.

use crate::NonlinearProblem;
use std::cell::{Cell, RefCell};

/// First source or stencil refusal. No Clone bound is imposed on source errors.
#[derive(Debug)]
pub enum FiniteDifferenceError<E> {
    /// The physical residual owner refused an evaluation.
    Source(E),
    /// Invalid dimensions, nonfinite values or an unrepresentable stencil.
    Invalid(&'static str),
}
impl<E: core::fmt::Display> core::fmt::Display for FiniteDifferenceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => write!(f, "numerical residual source: {e}"),
            Self::Invalid(s) => write!(f, "numerical residual: {s}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FiniteDifferenceError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Source(e) => Some(e), Self::Invalid(_) => None }
    }
}

/// A repeatable physical residual with bounded-direction numerical differences.
///
/// `difference_step` is the maximum-component input perturbation in the caller's
/// coordinates, not a dimensionless tolerance. Each Jacobian action makes at
/// most two residual calls. Callbacks must not advance accepted physical history.
/// Use `run_controlled` with `has_failed()` in its stop predicate, then inspect
/// `take_failure()` BEFORE interpreting any solver report or accepting a state.
/// A failure latches permanently, even after its original error has been taken.
pub struct FallibleFiniteDifference<F, E> {
    dimension: usize,
    difference_step: f64,
    evaluate: F,
    failed: Cell<bool>,
    failure: RefCell<Option<FiniteDifferenceError<E>>>,
}
impl<F, E> FallibleFiniteDifference<F, E>
where F: Fn(&[f64]) -> Result<Vec<f64>, E> {
    /// Admit dimensions and a caller-selected positive finite difference scale.
    pub fn new(dimension: usize, difference_step: f64, evaluate: F)
        -> Result<Self, FiniteDifferenceError<E>> {
        if dimension == 0 || !difference_step.is_finite() || difference_step <= 0.0 {
            return Err(FiniteDifferenceError::Invalid("positive dimension and finite positive difference scale required"));
        }
        Ok(Self { dimension, difference_step, evaluate,
            failed: Cell::new(false), failure: RefCell::new(None) })
    }
    /// Whether any source or numerical callback has refused this solve.
    pub fn has_failed(&self) -> bool { self.failed.get() }
    /// Move the first concrete refusal out without permitting more evaluations.
    pub fn take_failure(&self) -> Option<FiniteDifferenceError<E>> {
        self.failure.borrow_mut().take()
    }
    fn fail(&self, error: FiniteDifferenceError<E>) {
        if !self.failed.replace(true) { *self.failure.borrow_mut() = Some(error); }
    }
    fn checked(&self, x: &[f64]) -> Option<Vec<f64>> {
        if self.has_failed() { return None; }
        if x.len() != self.dimension || x.iter().any(|v| !v.is_finite()) {
            self.fail(FiniteDifferenceError::Invalid("finite complete residual coordinates required"));
            return None;
        }
        match (self.evaluate)(x) {
            Err(e) => { self.fail(FiniteDifferenceError::Source(e)); None },
            Ok(y) if y.len() != self.dimension || y.iter().any(|v| !v.is_finite()) => {
                self.fail(FiniteDifferenceError::Invalid("residual owner returned incomplete or nonfinite values")); None
            },
            Ok(y) => Some(y),
        }
    }
}
impl<F, E> NonlinearProblem for FallibleFiniteDifference<F, E>
where F: Fn(&[f64]) -> Result<Vec<f64>, E> {
    fn dimension(&self) -> usize { self.dimension }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if output.len() != self.dimension {
            self.fail(FiniteDifferenceError::Invalid("residual output dimension mismatch")); return;
        }
        if let Some(y) = self.checked(x) { output.copy_from_slice(&y); }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.has_failed() { return; }
        if x.len() != self.dimension || direction.len() != self.dimension || output.len() != self.dimension
            || x.iter().chain(direction).any(|v| !v.is_finite()) {
            self.fail(FiniteDifferenceError::Invalid("finite complete Jacobian coordinates required")); return;
        }
        let scale = direction.iter().copied().map(f64::abs).fold(0.0, f64::max);
        if scale == 0.0 { output.fill(0.0); return; }
        let delta = self.difference_step / scale;
        if !delta.is_finite() || delta <= 0.0 {
            self.fail(FiniteDifferenceError::Invalid("unrepresentable Jacobian difference scale")); return;
        }
        let shifted: Vec<_> = x.iter().zip(direction).map(|(&a, &d)| a + delta*d).collect();
        if shifted.iter().any(|v| !v.is_finite()) || shifted.as_slice() == x {
            self.fail(FiniteDifferenceError::Invalid("Jacobian stencil is nonfinite or does not change coordinates")); return;
        }
        let Some(base) = self.checked(x) else { return; };
        let Some(other) = self.checked(&shifted) else { return; };
        for ((value, a), b) in output.iter_mut().zip(base).zip(other) { *value = (b-a)/delta; }
        if output.iter().any(|v| !v.is_finite()) {
            self.fail(FiniteDifferenceError::Invalid("Jacobian action is unrepresentable")); output.fill(f64::NAN);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::Infallible;
    #[test]
    fn affine_directional_action_matches_the_original_operator() {
        let p = FallibleFiniteDifference::new(2, 1e-5,
            |x: &[f64]| Ok::<_, Infallible>(vec![2.0*x[0]-3.0*x[1], x[0]+4.0*x[1]])).unwrap();
        let mut y = [0.0; 2]; p.jacobian_apply(&[2.0, -1.0], &[0.5, 2.0], &mut y);
        assert!((y[0]+5.0).abs() < 1e-8 && (y[1]-8.5).abs() < 1e-8);
        p.jacobian_apply(&[2.0, -1.0], &[0.0; 2], &mut y); assert_eq!(y, [0.0; 2]);
        assert!(!p.has_failed());
    }
    #[test]
    fn nonclone_source_refusal_is_retained_and_latched() {
        #[derive(Debug, PartialEq)] struct Refusal(u32);
        let calls = Cell::new(0);
        let p = FallibleFiniteDifference::new(1, 0.1, |x: &[f64]| {
            calls.set(calls.get()+1);
            if x[0]>1.0 { Err(Refusal(7)) } else { Ok(vec![x[0]]) }
        }).unwrap();
        let mut y = [0.0]; p.jacobian_apply(&[1.0], &[1.0], &mut y);
        assert!(y[0].is_nan() && p.has_failed());
        assert!(matches!(p.take_failure(), Some(FiniteDifferenceError::Source(Refusal(7)))));
        let before = calls.get(); p.residual(&[0.0], &mut y);
        assert_eq!(calls.get(), before); assert!(y[0].is_nan()); assert!(p.has_failed());
    }
    #[test]
    fn bad_residuals_and_stagnant_stencils_do_not_become_convergence() {
        for invalid in [vec![], vec![f64::NAN]] {
            let p = FallibleFiniteDifference::new(1, 0.1,
                |_: &[f64]| Ok::<_, Infallible>(invalid.clone())).unwrap();
            let mut y = [0.0]; p.residual(&[0.0], &mut y);
            assert!(p.has_failed() && y[0].is_nan());
        }
        let p = FallibleFiniteDifference::new(1, f64::MIN_POSITIVE,
            |x: &[f64]| Ok::<_, Infallible>(x.to_vec())).unwrap();
        let mut y = [0.0]; p.jacobian_apply(&[1.0], &[1.0], &mut y);
        assert!(p.has_failed() && y[0].is_nan());
    }
}
