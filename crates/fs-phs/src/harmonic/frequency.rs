//! Exact frequency differentiation of the already factored physical pencil.
//!
//! For A(omega) = -i omega E - (J-R)W and A x = G u, the tangent satisfies
//! A x' = G u' + i E x. Only differential coordinates enter E x; algebraic
//! efforts and every hereditary state remain in the shared native solve.

use super::{
    C64, HarmonicError, HarmonicFactor, HarmonicResponse, Sum, finite, invalid, poll,
    resource_bounds,
};

/// Original response and its derivative with respect to angular frequency.
///
/// All derivative fields are per rad/s, in the original coordinate/port units.
/// The physical model and its equilibrium remain fixed; an optional supplied
/// input derivative can describe a frequency-dependent excitation. No finite
/// difference, pole fit, modal reduction or differentiation of LU iterations is
/// used. The residual is a backward-error test, not a near-pole forward bound.
#[derive(Clone, Debug, PartialEq)]
pub struct HarmonicFrequencyResponse {
    /// Identical physical response to [`HarmonicFactor::response`].
    pub response: HarmonicResponse,
    /// Differential coordinates and algebraic reactions per rad/s.
    pub state_derivative: Vec<C64>,
    /// `[Q_diff x'_diff, x'_alg]` per rad/s.
    pub effort_derivative: Vec<C64>,
    /// Original power-conjugate output derivatives, in original port order.
    pub output_derivative: Vec<C64>,
    /// Derivative of cycle-mean stored energy [J per rad/s].
    pub mean_stored_energy_derivative: f64,
    /// Per-port supplied mean-power derivatives [W per rad/s].
    pub port_power_derivative: Vec<f64>,
    /// Derivative of aggregate original port exchange [W per rad/s], with
    /// coupling remainders retained through the effort contraction.
    pub supplied_power_derivative: f64,
    /// Independent derivative of `Re(e* R e + e* T x)/2` [W per rad/s].
    /// Nonlinear equilibrium images retain signed incremental work; this is
    /// not a derivative of the biased source's total heat.
    pub dissipated_power_derivative: f64,
    /// Dissipation derivative minus supply derivative [W per rad/s].
    pub power_residual_derivative: f64,
    /// Componentwise backward residual of `A x' = G u' + i E x`.
    pub relative_derivative_residual: f64,
}

impl HarmonicFrequencyResponse {
    /// `y'/y` at one original port, without an amplitude floor. The real part
    /// is the derivative of log amplitude; the imaginary part is the phase
    /// slope in seconds. With `exp(-i omega t)`, a pure delay has POSITIVE
    /// phase slope. This includes any supplied input-frequency dependence;
    /// it is a transfer delay only when the input phasors are held fixed.
    /// Exactly zero output has no defined phase/log amplitude and returns None.
    pub fn output_logarithmic_derivative(&self, port: usize) -> Result<Option<C64>, HarmonicError> {
        let output = *self
            .response
            .output
            .get(port)
            .ok_or_else(|| invalid("frequency-slope port is out of range"))?;
        let derivative = *self
            .output_derivative
            .get(port)
            .ok_or_else(|| invalid("frequency-slope derivative port is out of range"))?;
        if !finite(output) || !finite(derivative) {
            return Err(invalid("finite output and frequency derivative required"));
        }
        if output == C64::ZERO {
            return Ok(None);
        }
        let slope = derivative / output;
        if !finite(slope) {
            return Err(invalid("output frequency slope is not representable"));
        }
        Ok(Some(slope))
    }
}

impl HarmonicFactor<'_> {
    /// Evaluate the original driven response and its continuous frequency tangent.
    ///
    /// `input_derivative` gives `du/domega` for every original source port;
    /// None holds every input phasor fixed. Full coupled storage, dissipation,
    /// internal memory and algebraic reactions remain in the same original
    /// equations. The derivative reuses this factor's LU and row/column scales,
    /// including the native bounded residual refinement. No new factor is made.
    ///
    /// This operation admits twice the ordinary response work score and the
    /// ordinary payload bound plus `128*n + 64*ports` bytes for retained
    /// response/tangent arrays before allocation. Cancellation, invalid data,
    /// arithmetic or residual refusal publishes neither a partial response nor
    /// a derivative and leaves the reusable factor unchanged.
    pub fn response_with_frequency_derivative(
        &self,
        input: &[C64],
        input_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFrequencyResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = self.model;
        let (n, m) = (model.n, model.m);
        if input_derivative
            .is_some_and(|values| values.len() != m || values.iter().any(|&value| !finite(value)))
        {
            return Err(invalid(
                "one finite input derivative per source port required",
            ));
        }
        let (base_bytes, base_work) = resource_bounds(n, m)?;
        let bytes = n
            .checked_mul(128)
            .and_then(|v| m.checked_mul(64).and_then(|p| v.checked_add(p)))
            .and_then(|v| base_bytes.checked_add(v))
            .ok_or_else(|| invalid("frequency derivative payload bound overflow"))?;
        let work = base_work
            .checked_mul(2)
            .ok_or_else(|| invalid("frequency derivative work bound overflow"))?;
        if bytes > model.budget.max_dense_bytes || work > model.budget.max_dense_work {
            return Err(invalid(
                "insufficient frequency derivative byte or work budget",
            ));
        }
        let response = self.response(input, &mut cancelled)?;
        let mut rhs = vec![C64::ZERO; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            let initial = if i < model.n_diff {
                C64::new(-response.state[i].im, response.state[i].re)
            } else { C64::ZERO };
            let terms = super::input_tangent_sums(
                &model.g[i*m..(i+1)*m], None, input, input_derivative, |_| true,
            )?;
            rhs[i] = super::input_value(terms, initial)?;
        }
        let (state_derivative, relative_derivative_residual) =
            self.solve_rhs(&rhs, &mut cancelled)?;
        let mut effort_derivative = state_derivative.clone();
        let mut energy_derivative = Sum::default();
        for i in 0..model.n_diff {
            poll(&mut cancelled)?;
            let original = super::storage_row(model, i, &response.state)?;
            let derivative = super::storage_row(model, i, &state_derivative)?;
            effort_derivative[i] = super::input_value(derivative, C64::ZERO)?;
            super::add_stored_work(&mut energy_derivative, original, state_derivative[i])?;
            super::add_stored_work(&mut energy_derivative, derivative, response.state[i])?;
        }
        let energy_derivative = energy_derivative.value()
            .ok_or_else(|| invalid("stored-energy tangent is not representable"))?;
        let mut loss_derivative = Sum::default();
        let mut output_derivative = vec![C64::ZERO; m];
        for i in 0..n {
            poll(&mut cancelled)?;
            let original = super::resisting_row(model, i, &response.state, &response.effort)?;
            let derivative = super::resisting_tangent_row(model, i, &response.state,
                &response.effort, &state_derivative, &effort_derivative,
                super::HarmonicParameterDirection::default())?;
            super::add_dissipative_work(&mut loss_derivative, original, effort_derivative[i])?;
            super::add_dissipative_work(&mut loss_derivative, derivative, response.effort[i])?;
        }
        let loss_derivative = loss_derivative.value()
            .ok_or_else(|| invalid("dissipative-power tangent is not representable"))?;
        for (p, value) in output_derivative.iter_mut().enumerate() {
            poll(&mut cancelled)?;
            *value = super::port_output(&model.g, m, p, &effort_derivative)?;
        }
        let mut port_power_derivative = Vec::with_capacity(m);
        for p in 0..m {
            poll(&mut cancelled)?;
            let derivative = super::peak_power_tangent(input[p], response.output[p],
                input_derivative.map(|values| values[p]), output_derivative[p])?;
            port_power_derivative.push(derivative);
        }
        let supply_derivative = super::physical_power_tangent(model, input, input_derivative,
            &response.effort, &effort_derivative, None, &mut cancelled)?;
        let balance_derivative = loss_derivative - supply_derivative;
        if effort_derivative
            .iter()
            .chain(&output_derivative)
            .any(|&v| !finite(v))
            || port_power_derivative.iter().any(|v| !v.is_finite())
            || [
                energy_derivative,
                loss_derivative,
                supply_derivative,
                balance_derivative,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid("frequency derivative observation overflow"));
        }
        poll(&mut cancelled)?;
        Ok(HarmonicFrequencyResponse {
            response,
            state_derivative,
            effort_derivative,
            output_derivative,
            mean_stored_energy_derivative: energy_derivative,
            port_power_derivative,
            supplied_power_derivative: supply_derivative,
            dissipated_power_derivative: loss_derivative,
            power_residual_derivative: balance_derivative,
            relative_derivative_residual,
        })
    }
}
