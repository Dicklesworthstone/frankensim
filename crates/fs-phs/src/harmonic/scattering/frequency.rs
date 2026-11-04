//! Frequency tangents of actual loaded physical inputs and outgoing waves.
//! The reference closes the native descriptor before either solve. It remains
//! boundary exchange, never extra material dissipation in the observations.

use super::{
    C64, HarmonicError, HarmonicPortPower, ScatteringFactor, ScatteringModel, ScatteringResponse,
    finite, invalid, poll,
};
use crate::harmonic::HarmonicFrequencyResponse;

/// Derivatives at one selected wave port, per rad/s at fixed reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringWaveDerivative {
    /// Original input/output slot; order matches `ScatteringResponse::waves`.
    pub original_port: usize,
    /// Derivative of the incident output-unit amplitude a.
    pub incident: C64,
    /// Derivative of b=y-a, including the direct incident-wave contribution.
    pub outgoing: C64,
    /// Derivative of |a|^2/(2Z) [W per rad/s].
    pub incident_power_w: f64,
    /// Independent derivative of |b|^2/(2Z) [W per rad/s].
    pub outgoing_power_w: f64,
    /// Derivative of original Re(conj(u)*y)/2 [W per rad/s].
    pub physical_power_w: f64,
    /// Incident minus outgoing minus physical power derivatives [W per rad/s].
    pub power_residual_w: f64,
}

/// Loaded response and its continuous frequency tangent in original coordinates.
///
/// All derivatives are per rad/s. Material coefficients, equilibrium and wave
/// references remain fixed. Ordinary source and incident-wave derivatives may
/// be supplied explicitly. Physical work uses actual `u=(2a-y)/Z`; the effective
/// `2a/Z` drive is never reported as a physical boundary source.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringFrequencyResponse {
    /// Same physical response as [`ScatteringFactor::response`].
    pub response: ScatteringResponse,
    pub state_derivative: Vec<C64>,
    pub effort_derivative: Vec<C64>,
    pub output_derivative: Vec<C64>,
    /// Actual original input derivatives, including the reciprocal reaction
    /// `(2a'-y')/Z` at selected ports. Ordinary slots retain authored derivatives.
    pub physical_input_derivative: Vec<C64>,
    /// Selected waves in the same authored order as the original response.
    pub wave_derivatives: Vec<ScatteringWaveDerivative>,
    pub mean_stored_energy_derivative: f64,
    pub port_power_derivative: Vec<f64>,
    pub supplied_power_derivative: f64,
    /// Original internal dissipation/incremental work derivative only. The
    /// reference-boundary loading is excluded, just as in the base response.
    pub dissipated_power_derivative: f64,
    pub power_residual_derivative: f64,
    /// Internal loss minus incident/outgoing and ordinary source derivatives.
    pub wave_power_residual_derivative: f64,
    /// Maximum of the closed tangent solve and independently checked original
    /// equations with actual input derivatives. No forward-error claim near poles.
    pub relative_derivative_residual: f64,
}

impl ScatteringFactor<'_> {
    /// Differentiate the complete closed physical response, reusing its native LU.
    ///
    /// Original input slots contain incident output waves on selected ports and
    /// ordinary inputs elsewhere. `input_derivative` follows the same convention;
    /// None means fixed phasors. Full descriptor constraints, mutual storage and
    /// material memory remain in both original-coordinate solves. References are
    /// held fixed; frequency-dependent physical loads require explicit retained
    /// states in the supplied system, not a changing numerical wave reference.
    ///
    /// In addition to the shared harmonic derivative admission, this operation
    /// admits twice the constructor's additional scattering payload/work bound
    /// against its original `ScatteringBudget`. No partial result escapes a
    /// cancellation, input, budget, representability or physical-residual refusal.
    pub fn response_with_frequency_derivative(
        &self,
        input: &[C64],
        input_derivative: Option<&[C64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringFrequencyResponse, HarmonicError> {
        poll(&mut cancelled)?;
        self.model.admit_frequency_observation()?;
        let effective = self.effective_input(input, &mut cancelled)?;
        let effective_derivative = input_derivative
            .map(|values| self.effective_input(values, &mut cancelled))
            .transpose()?;
        let raw = self.factor.response_with_frequency_derivative(
            &effective,
            effective_derivative.as_deref(),
            &mut cancelled,
        )?;
        let mut physical_input_derivative =
            input_derivative.map_or_else(|| vec![C64::ZERO; input.len()], <[C64]>::to_vec);
        for &(p, z) in &self.model.ports {
            poll(&mut cancelled)?;
            physical_input_derivative[p] =
                (physical_input_derivative[p].scale(2.0) - raw.output_derivative[p]).scale(1.0 / z);
        }
        if physical_input_derivative.iter().any(|&u| !finite(u)) {
            return Err(invalid("actual wave input derivative is not representable"));
        }
        let relative_derivative_residual =
            raw.relative_derivative_residual
                .max(self.model.original_frequency_residual(
                    &raw,
                    &physical_input_derivative,
                    &mut cancelled,
                )?);
        let response = self.observe_response(input, raw.response, &mut cancelled)?;
        let mut port_power_derivative = Vec::with_capacity(input.len());
        for p in 0..input.len() {
            poll(&mut cancelled)?;
            port_power_derivative.push(
                half_pair(physical_input_derivative[p], response.physical.output[p])?
                    + half_pair(response.physical_inputs[p], raw.output_derivative[p])?,
            );
        }
        let mut wave_derivatives = Vec::with_capacity(self.model.ports.len());
        let mut wave_supply_derivative = 0.0;
        for wave in &response.waves {
            poll(&mut cancelled)?;
            let p = wave.original_port;
            let incident = input_derivative.map_or(C64::ZERO, |values| values[p]);
            let outgoing = raw.output_derivative[p] - incident;
            let incident_power_w =
                2.0 * half_pair(wave.incident.scale(1.0 / wave.reference), incident)?;
            let outgoing_power_w =
                2.0 * half_pair(wave.outgoing.scale(1.0 / wave.reference), outgoing)?;
            let physical_power_w = port_power_derivative[p];
            let power_residual_w = incident_power_w - outgoing_power_w - physical_power_w;
            if !finite(outgoing)
                || [
                    incident_power_w,
                    outgoing_power_w,
                    physical_power_w,
                    power_residual_w,
                ]
                .iter()
                .any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "wave amplitude or power derivative is not representable",
                ));
            }
            wave_supply_derivative += incident_power_w - outgoing_power_w;
            wave_derivatives.push(ScatteringWaveDerivative {
                original_port: p,
                incident,
                outgoing,
                incident_power_w,
                outgoing_power_w,
                physical_power_w,
                power_residual_w,
            });
        }
        for (p, &power) in port_power_derivative.iter().enumerate() {
            poll(&mut cancelled)?;
            if !self.model.selected[p] {
                wave_supply_derivative += power;
            }
        }
        let supplied_power_derivative: f64 = port_power_derivative.iter().sum();
        let dissipated_power_derivative = raw.dissipated_power_derivative;
        let power_residual_derivative = dissipated_power_derivative - supplied_power_derivative;
        let wave_power_residual_derivative = dissipated_power_derivative - wave_supply_derivative;
        if port_power_derivative.iter().any(|x| !x.is_finite())
            || [
                supplied_power_derivative,
                power_residual_derivative,
                wave_power_residual_derivative,
            ]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "physical scattering power derivative is not representable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ScatteringFrequencyResponse {
            response,
            state_derivative: raw.state_derivative,
            effort_derivative: raw.effort_derivative,
            output_derivative: raw.output_derivative,
            physical_input_derivative,
            wave_derivatives,
            mean_stored_energy_derivative: raw.mean_stored_energy_derivative,
            port_power_derivative,
            supplied_power_derivative,
            dissipated_power_derivative,
            power_residual_derivative,
            wave_power_residual_derivative,
            relative_derivative_residual,
        })
    }
}

impl ScatteringModel {
    fn admit_frequency_observation(&self) -> Result<(), HarmonicError> {
        let (n, m, k) = (
            self.closed.n as u128,
            self.closed.m as u128,
            self.ports.len() as u128,
        );
        let bytes = 2 * (8 * n * k + 256 * k + 64 * m);
        let work = 2 * (2 * n * n * (k + 1) + n * m + k * k + 32 * m);
        if bytes > self.budget.maximum_payload_bytes as u128
            || work > self.budget.maximum_work as u128
        {
            return Err(invalid(
                "insufficient loaded frequency-derivative observation resources",
            ));
        }
        Ok(())
    }

    // Independent derivative of the ORIGINAL unclosed physical equations:
    // -i omega E dx - (J-R) de + T dx - G du - i E x = 0.
    // The closed generator is deliberately excluded from this sign/work check.
    fn original_frequency_residual(
        &self,
        value: &HarmonicFrequencyResponse,
        input_derivative: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, HarmonicError> {
        let original = &self.closed;
        let mut relative = 0.0_f64;
        let laplace = C64::new(0.0, -value.response.omega_rad_s);
        for i in 0..original.n {
            poll(cancelled)?;
            let mut residual = C64::ZERO;
            let mut scale = 0.0;
            if i < original.n_diff {
                residual = laplace * value.state_derivative[i];
                let forcing = C64::new(-value.response.state[i].im, value.response.state[i].re);
                scale = residual.abs() + forcing.abs();
                residual = residual - forcing;
            }
            for (j, &e) in value.effort_derivative.iter().enumerate() {
                let term = e.scale(original.structure[i * original.n + j]);
                residual = residual - term;
                scale += term.abs();
            }
            if original.resisting_tangent.is_some() {
                let term = original.nonlinear_loss_row(i, &value.state_derivative);
                residual = residual + term;
                scale += term.abs();
            }
            for (p, &u) in input_derivative.iter().enumerate() {
                let term = u.scale(original.g[i * original.m + p]);
                residual = residual - term;
                scale += term.abs();
            }
            if !finite(residual) || !scale.is_finite() {
                return Err(invalid(
                    "original scattering tangent residual is not representable",
                ));
            }
            relative = relative.max(if scale == 0.0 {
                0.0
            } else {
                residual.abs() / scale
            });
        }
        if relative > original.budget.relative_residual_tolerance {
            return Err(HarmonicError::Residual {
                relative,
                tolerance: original.budget.relative_residual_tolerance,
            });
        }
        Ok(relative)
    }
}

fn half_pair(left: C64, right: C64) -> Result<f64, HarmonicError> {
    Ok(HarmonicPortPower::from_peak(left, right)?.active_w)
}
