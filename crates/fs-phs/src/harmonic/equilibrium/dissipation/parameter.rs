//! Complete fixed-state parameter chain of the original nonlinear loss law.

use super::{LinearizationError, NonlinearDissipation, invalid, poll};
use crate::interconnection::accumulation::Sum;

/// Explicit partial `D_p(x,e)` at FIXED original state and FIXED effort.
/// Storage and pump effort changes are supplied separately and enter through
/// the original `D_e` action exactly once. This signed derivative need not be
/// dissipative. Every coordinate, including material memory, must be written.
/// The callback is smooth, deterministic and read-only; refuse corners or
/// unsupported source states rather than choosing a one-sided surrogate.
pub type DissipationParameterAction<'a> = dyn Fn(&[f64], &[f64], &mut [f64]) -> bool + 'a;

/// A successful original-source sample, borrowed until the next mutable query.
/// The separate flow terms let the physical equation owner combine all J/R/G
/// and resisting contributions before its final rounding. A rounded sum of
/// the two derivatives can erase weak motion when other large terms oppose it.
pub struct DissipationParameterEvaluation<'a> {
    pub flow: &'a [f64],
    /// Original D_e * e_p, with the independent state direction held at zero.
    pub effort_chain: &'a [f64],
    /// Explicit original D_p. Exactly zero when the constitutive law is fixed.
    pub direct_flow: &'a [f64],
    pub power: f64,
    /// e_p dot D + e dot (D_e e_p + D_p), independently contracted with all
    /// products retained. It is signed and never clamped to nonnegative heat.
    pub power_derivative: f64,
}

/// Reusable original nonlinear-loss parameter direction. Preparation stores
/// four n-vectors; warm evaluation has no owner allocations and makes one
/// complete original tangent action. Callback resources remain caller-owned.
/// Native harmonic dimensions are bounded independently of any solver budget.
pub struct PreparedDissipationParameterDirection<'a, 'p> {
    law: NonlinearDissipation<'a>,
    direct: Option<&'p DissipationParameterAction<'p>>,
    zero: Vec<f64>,
    flow: Vec<f64>,
    effort_chain: Vec<f64>,
    direct_flow: Vec<f64>,
}

impl<'a> NonlinearDissipation<'a> {
    /// Bind one signed constitutive parameter direction to THIS loss law.
    /// None explicitly holds the law fixed at fixed state/effort; it does not
    /// suppress D_e e_p from changes to the original storage or pump.
    pub fn prepare_parameter_direction<'p>(
        &self,
        dimension: usize,
        direct: Option<&'p DissipationParameterAction<'p>>,
    ) -> Result<PreparedDissipationParameterDirection<'a, 'p>, LinearizationError> {
        if dimension == 0
            || dimension > crate::harmonic::MAX_HARMONIC_DIMENSION
            || self.max_tangent_actions == 0
        {
            return Err(invalid(
                "bounded complete loss parameter dimension and one tangent action required",
            ));
        }
        Ok(PreparedDissipationParameterDirection {
            law: *self,
            direct,
            zero: vec![0.0; dimension],
            flow: vec![0.0; dimension],
            effort_chain: vec![0.0; dimension],
            direct_flow: vec![0.0; dimension],
        })
    }
}

impl PreparedDissipationParameterDirection<'_, '_> {
    /// Evaluate at fixed state. The complete effort derivative may contain
    /// storage parameters, pump coefficients and explicit memory couplings.
    /// Unknown waveform changes belong to the separate original state tangent.
    /// Original e dot D admission precedes the derivative; passivity of the
    /// original sample does not imply a nonnegative parameter or state slope.
    /// Cancellation/refusal returns no sample; subsequent queries are reusable.
    pub fn evaluate(
        &mut self,
        state: &[f64],
        effort: &[f64],
        effort_derivative: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipationParameterEvaluation<'_>, LinearizationError> {
        poll(&mut cancelled)?;
        let n = self.zero.len();
        if state.len() != n
            || effort.len() != n
            || effort_derivative.len() != n
            || state
                .iter()
                .chain(effort)
                .chain(effort_derivative)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite fixed-state loss parameter inputs required",
            ));
        }
        let power = self
            .law
            .evaluate_into(state, effort, &mut self.flow, &mut cancelled)?;
        self.effort_chain.fill(f64::NAN);
        if !(self.law.tangent)(
            state,
            effort,
            &self.zero,
            effort_derivative,
            &mut self.effort_chain,
        ) || self.effort_chain.iter().any(|v| !v.is_finite())
        {
            return Err(LinearizationError::DissipationParameterRefused {
                what: "complete original effort-direction loss action required",
            });
        }
        poll(&mut cancelled)?;
        if let Some(direct) = self.direct {
            self.direct_flow.fill(f64::NAN);
            if !direct(state, effort, &mut self.direct_flow)
                || self.direct_flow.iter().any(|v| !v.is_finite())
            {
                return Err(LinearizationError::DissipationParameterRefused {
                    what: "complete original fixed-state/fixed-effort loss partial required",
                });
            }
        } else {
            self.direct_flow.fill(0.0);
        }
        poll(&mut cancelled)?;
        let mut derivative = Sum::default();
        for i in 0..n {
            derivative
                .product(effort_derivative[i], self.flow[i])
                .and_then(|_| derivative.product(effort[i], self.effort_chain[i]))
                .and_then(|_| derivative.product(effort[i], self.direct_flow[i]))
                .ok_or_else(|| invalid("nonlinear loss parameter work is unrepresentable"))?;
        }
        let power_derivative = derivative
            .value()
            .ok_or_else(|| invalid("nonlinear loss parameter work is unrepresentable"))?;
        poll(&mut cancelled)?;
        Ok(DissipationParameterEvaluation {
            flow: &self.flow,
            effort_chain: &self.effort_chain,
            direct_flow: &self.direct_flow,
            power,
            power_derivative,
        })
    }
}
