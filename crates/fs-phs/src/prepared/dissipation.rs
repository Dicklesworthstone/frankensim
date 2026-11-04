//! Passive or externally supplied flow in the SAME implicit Gonzalez equation.
//! No explicit lag, operator split, changed storage, or second time integrator.
use super::{StepWorkspace, PreparedStepError, PreparedStepRecord, dimensions, norm};
use crate::PortHamiltonian;

type Force<'a> = &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> bool;
type DrivenForce<'a> = &'a dyn Fn(&[f64], &[f64], &mut [f64]) -> Option<NonlinearPortPower>;
type Tangent<'a> = &'a dyn Fn(&[f64], &[f64], &[f64], &[f64], &mut [f64]) -> bool;

#[derive(Clone, Copy)]
pub(super) enum Flow<'a> {
    Passive(Force<'a>),
    Driven(DrivenForce<'a>),
}

#[derive(Clone, Copy)]
pub(super) struct Dissipation<'a> {
    pub force: Flow<'a>,
    pub tangent: Option<Tangent<'a>>,
}

/// Physical powers from the SAME nonlinear flow evaluation and joint effort.
/// The explicit loss preserves small positive dissipation beside much larger
/// opposing storage/supply powers. The owner must satisfy
/// `effort dot resisting_flow + supplied_power = dissipated_power`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct NonlinearPortPower {
    /// Signed external power into the system, in watts.
    pub supplied_power_w: f64,
    /// Nonnegative physical dissipated power, in watts.
    pub dissipated_power_w: f64,
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

// Check the supplied owner's loss without reconstructing a small heat rate by
// subtracting two large powers. Passive callers retain their original path.
fn admit_supplied_power(effort: &[f64], force: &[f64], power: NonlinearPortPower)
    -> Result<(), crate::PhsError>
{
    let NonlinearPortPower { supplied_power_w: supplied, dissipated_power_w: loss } = power;
    norm(&[supplied, loss])?;
    if loss < 0.0 { return Err(dimensions("nonlinear supplied port declares negative loss")); }
    let mut value = 0.0;
    let mut absolute = 0.0;
    for (&e, &f) in effort.iter().zip(force) {
        let product = e*f;
        value += product;
        absolute += product.abs();
    }
    let defect = (value + supplied) - loss;
    absolute += supplied.abs() + loss;
    norm(&[value, absolute, defect])?;
    let tolerance = 8.0*f64::EPSILON*(effort.len()+2) as f64*absolute;
    if !tolerance.is_finite() || defect.abs() > tolerance {
        return Err(dimensions("nonlinear supplied port violates its power identity"));
    }
    Ok(())
}

impl Dissipation<'_> {
    pub(super) fn evaluate(self, state: &[f64], effort: &[f64], out: &mut [f64])
        -> Result<NonlinearPortPower, crate::PhsError>
    {
        out.fill(0.0);
        let powers = match self.force {
            Flow::Passive(force) => {
                if !force(state, effort, out) {
                    return Err(dimensions("nonlinear dissipative port refused its trial"));
                }
                norm(out)?;
                NonlinearPortPower { supplied_power_w: 0.0, dissipated_power_w: power(effort,out)? }
            }
            Flow::Driven(force) => {
                let powers = force(state, effort, out)
                    .ok_or_else(|| dimensions("nonlinear supplied port refused its trial"))?;
                norm(out)?;
                admit_supplied_power(effort,out,powers)?;
                powers
            }
        };
        Ok(powers)
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
    /// With no Hessian, an available source-owned effort Jacobian and supplied
    /// dissipative tangent enable exact Newton; otherwise the finite-difference
    /// residual probes include D automatically. Explicit analytic calls never
    /// silently fall back to finite differences. Both paths retain the existing
    /// solver limits and physical equation.
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
            Some(Dissipation { force: Flow::Passive(dissipation), tangent }), cancelled)
    }

    /// Advance a nonlinear externally supplied port in the ORIGINAL implicit
    /// equation `dx = dt*((J-R)*dg - D(midpoint,dg) + G*u)`.
    ///
    /// `port` fills resisting flow D and returns its signed external supplied
    /// power and nonnegative physical dissipation, in watts. `None` refuses.
    /// Every trial checks `dg dot D + supplied_power = dissipated_power` within
    /// arithmetic roundoff; undeclared active power still refuses. Returning
    /// both powers preserves a small loss beside large opposing storage/supply
    /// powers without subtractive reconstruction. An independently pressurized fluid
    /// port can therefore admit reservoir work while retaining its positive
    /// viscous loss. The ordinary linear `G*u` input remains separate.
    ///
    /// Both rates come from the SAME callback and discrete effort used for the
    /// accepted equation. The receipt adds `dt*supplied_power` to supplied work
    /// and `dt*dissipated_power` to dissipation. The callback must
    /// include every nonlinear source exactly once and hold source parameters
    /// fixed during the step. It owns physical source validity and work units;
    /// this accounting does not certify a caller-authored constitutive law.
    ///
    /// The tangent differentiates the COMPLETE D, including state-dependent
    /// source coupling, at independent state and discrete-effort directions.
    /// No derivative of the scalar supply is needed for the state equation.
    /// Analytic/source-owned/finite-difference Newton, condensation, bounded
    /// backtracking and cancellation use the same existing solver. Callbacks
    /// are immutable trial queries; no history or caller output publishes on
    /// failure. Prepared stepping adds no solver allocation or extra port solve.
    ///
    /// # Errors
    /// The existing implicit-step refusals plus refused/nonfinite supplied
    /// power, negative physical loss or an inconsistent conjugate power identity.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_into_supplied_controlled<F: FnMut() -> bool>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64],
        port: &dyn Fn(&[f64], &[f64], &mut [f64]) -> Option<NonlinearPortPower>,
        hessian: Option<&dyn Fn(&[f64], &[f64], &mut [f64]) -> bool>,
        tangent: Option<&dyn Fn(&[f64], &[f64], &[f64], &[f64], &mut [f64]) -> bool>,
        cancelled: F,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        self.step_with_hessian(sys, x0, u, dt, x_next, y, hessian,
            Some(Dissipation { force: Flow::Driven(port), tangent }), cancelled)
    }
}

#[cfg(test)]
#[path = "dissipation_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "supplied_tests.rs"]
mod supplied_tests;
