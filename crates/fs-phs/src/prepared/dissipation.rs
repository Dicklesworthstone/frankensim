//! An additional passive force in the SAME implicit Gonzalez equation.
//! No explicit lag, operator split, changed storage, or second time integrator.
use super::{StepWorkspace, PreparedStepError, PreparedStepRecord, dimensions, norm};
use crate::PortHamiltonian;

type Force<'a> = &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool;
type Tangent<'a> = &'a dyn Fn(&[f64], &[f64], &[f64], &[f64], &mut [f64]) -> bool;

#[derive(Clone, Copy)]
pub(super) struct Dissipation<'a> {
    pub force: Force<'a>,
    pub tangent: Option<Tangent<'a>>,
}

/// Conjugate power, not the norm of a force. Refuse a materially active port.
/// A last-bit negative dot product within its summation roundoff is numerical
/// zero; the caller's independent Hamiltonian/work gate remains unchanged.
pub(super) fn power(effort: &[f64], force: &[f64]) -> Result<f64, crate::PhsError> {
    let mut value = 0.0;
    let mut absolute = 0.0;
    for (&e, &f) in effort.iter().zip(force) {
        let product = e*f;
        value += product;
        absolute += product.abs();
    }
    norm(&[value, absolute])?;
    let tolerance = 8.0*f64::EPSILON*effort.len() as f64*absolute;
    if !tolerance.is_finite() || value < -tolerance {
        return Err(dimensions("nonlinear dissipative port supplies energy"));
    }
    Ok(value.max(0.0))
}

impl Dissipation<'_> {
    pub(super) fn evaluate(self, state: &[f64], effort: &[f64], out: &mut [f64])
        -> Result<(), crate::PhsError>
    {
        out.fill(0.0);
        if !(self.force)(state, effort, out) {
            return Err(dimensions("nonlinear dissipative port refused its trial"));
        }
        norm(out)?;
        power(effort, out)?;
        Ok(())
    }
    pub(super) fn directional(self, state: &[f64], effort: &[f64], ds: &[f64], de: &[f64],
        out: &mut [f64]) -> Result<(), crate::PhsError>
    {
        let tangent = self.tangent.ok_or_else(|| dimensions("missing nonlinear dissipative tangent"))?;
        out.fill(0.0);
        if !tangent(state, effort, ds, de, out) {
            return Err(dimensions("nonlinear dissipative tangent refused its trial"));
        }
        norm(out)?;
        Ok(())
    }
}

impl StepWorkspace {
    /// Solve `dx = dt*((J-R)*dg - D(midpoint,dg) + G*u)` in one Newton system.
    ///
    /// `dissipation` fills the positive resisting state-flow D. It must be a
    /// passive, read-only law: `dg dot D >= 0`. Its work uses the SAME discrete
    /// effort as the state equation, not an endpoint velocity or previous force.
    /// The returned loss includes both R and D; D is NOT external supplied work.
    /// No constitutive history is committed here. Returning false refuses.
    ///
    /// Passing `hessian` selects analytic Newton and requires `tangent`, which
    /// fills `D_state*state_direction + D_effort*effort_direction`. It receives
    /// independent directions, including the full Gonzalez energy correction.
    /// With no Hessian the canonical finite-difference residual probes include
    /// D automatically; an unused tangent is not evaluated. There is no silent
    /// analytic-to-FD fallback. Both paths retain the existing solver limits.
    ///
    /// Preparation allocates; calls add no solver allocation. This guarantee
    /// requires allocation-free physical callbacks. It is not a real-time claim.
    ///
    /// # Errors
    /// Original dimension, finite-state, Newton and cancellation refusals, plus
    /// a refused/nonfinite/active dissipative law or missing analytic tangent.
    /// Neither caller output buffer changes on any failure.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_into_dissipative_controlled<F: FnMut() -> bool>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64],
        dissipation: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
        hessian: Option<&dyn Fn(&[f64], &[f64], &mut [f64]) -> bool>,
        tangent: Option<&dyn Fn(&[f64], &[f64], &[f64], &[f64], &mut [f64]) -> bool>,
        cancelled: F,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        self.step_with_hessian(sys, x0, u, dt, x_next, y, hessian,
            Some(Dissipation { force: dissipation, tangent }), cancelled)
    }
}

#[cfg(test)]
#[path = "dissipation_tests.rs"]
mod tests;
