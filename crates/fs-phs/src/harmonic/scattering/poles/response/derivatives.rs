//! Exact original parameter actions of a fully closed outgoing mode.

use super::{HarmonicPoleResponseError, PreparedScatteringPole, invalid};
use crate::harmonic::scattering::parameter::{ParameterResidualState, add_parameter_source};
use crate::harmonic::{HarmonicParameterDirection, finite, parameter, poll};
use fs_math::c64::C64;

/// One local physical direction per caller-declared real parameter. All original
/// material/operating-point contributions belong in `physical`; references are
/// differentiated in authored wave-port order. Topology, state chart and
/// descriptor rank stay fixed. Source amplitudes cannot move a free linear pole.
#[derive(Clone, Copy, Debug, Default)]
pub struct ScatteringPoleParameterDirection<'a> {
    pub physical: HarmonicParameterDirection<'a>,
    pub reference_derivative: Option<&'a [f64]>,
}

/// Derivative of one outgoing boundary's exponential-mode account, in the
/// prepared fixed pivot gauge. Powers have no driven cycle-average factor 1/2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecayWaveParameterDerivative {
    pub original_port: usize,
    pub reference_derivative: f64,
    pub outgoing: C64,
    pub physical_input: C64,
    pub outgoing_power_w: f64,
    pub physical_power_w: f64,
    pub power_residual_w: f64,
}

/// Local derivatives per parameter unit of the loaded free pole and original
/// physical mode in v[pivot]=1 gauge. Stored-energy/work slopes depend on that
/// gauge; they are not derivatives of a driven harmonic power or unit-energy
/// normalization. Boundary work is kept separate from original R/T material work.
#[derive(Clone, Debug)]
pub struct ScatteringPoleDerivative {
    pub laplace_pole_per_s: C64,
    pub angular_frequency_rad_s: C64,
    pub state: Vec<C64>,
    pub effort: Vec<C64>,
    pub output: Vec<C64>,
    pub physical_inputs: Vec<C64>,
    pub waves: Vec<DecayWaveParameterDerivative>,
    pub stored_energy_j: f64,
    pub internal_dissipated_power_w: f64,
    pub outgoing_power_w: f64,
    pub total_decay_power_w: f64,
    pub decay_power_residual_w: f64,
    pub physical_power_residual_w: f64,
    pub relative_residual: f64,
    pub linear_residual: f64,
}

impl PreparedScatteringPole<'_> {
    /// Differentiate original dJ/dR/dQ/dG/dT and selected references through the
    /// complete free outgoing mode, including original reactions and memory.
    ///
    /// With u=-y/Z, the closed generator action is
    /// dL*v+dG*u+G*(-dG^T*e-G^T*dQ*v-u*dZ)/Z at selected ports. The existing
    /// shared source-action owner supplies both changing-G terms. The prepared
    /// bordered LU solves that action in its fixed pivot gauge; no boundary
    /// exchange is inserted into physical R/T or passed off as internal heat.
    ///
    /// Actual du, outgoing amplitudes, dH, original material work and outward
    /// power are observed independently. Every tangent rechecks the unclosed
    /// physical rows, including ds*E*v. Finite signed parameter directions are
    /// allowed; they do not admit a finite change, a changing state chart, a
    /// moving bias or a repeated-pole continuation. No finite difference, mode
    /// matching, extra eigenproblem or factorization is used. Batches are bounded
    /// and transactional; a refused/cancelled call leaves this factor reusable.
    pub fn parameter_derivatives(
        &self,
        directions: &[ScatteringPoleParameterDirection<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<ScatteringPoleDerivative>, HarmonicPoleResponseError> {
        poll(&mut cancelled)?;
        if directions.is_empty() || directions.len() > self.maximum_vectors {
            return Err(invalid(
                "complete nonempty loaded parameter batch within its ceiling required",
            ));
        }
        let model = &self.model.closed;
        let mut results = Vec::with_capacity(directions.len());
        for direction in directions {
            poll(&mut cancelled)?;
            direction.physical.validate(model, &mut cancelled)?;
            if direction.reference_derivative.is_some_and(|v| {
                v.len() != self.model.ports.len() || v.iter().any(|z| !z.is_finite())
            }) {
                return Err(invalid(
                    "complete finite selected-reference derivative required",
                ));
            }
            let mut action = parameter::generator_action(
                model,
                &direction.physical,
                self.right_state(),
                self.mode_effort(),
                &mut cancelled,
            )?;
            add_parameter_source(
                self.model,
                direction.physical,
                self.mode_effort(),
                &self.observation.physical_inputs,
                &action.storage,
                None,
                direction.reference_derivative,
                &mut action.generator,
                Some(&mut action.envelope),
                &mut cancelled,
            )?;
            let physical = self.physical.parameter_derivative_from_action(
                direction.physical,
                action,
                &mut cancelled,
            )?;
            let mut physical_inputs = vec![C64::ZERO; model.m];
            let mut waves = Vec::with_capacity(self.model.ports.len());
            let (mut outgoing_power, mut physical_power) = (0.0, 0.0);
            for (index, wave) in self.observation.waves.iter().enumerate() {
                poll(&mut cancelled)?;
                let p = wave.original_port;
                let z = wave.reference;
                let dz = direction.reference_derivative.map_or(0.0, |v| v[index]);
                let output = physical.output[p];
                let input = -output.scale(1.0 / z) - wave.physical_input.scale(dz / z);
                physical_inputs[p] = input;
                let outward = 2.0 * (wave.outgoing.scale(1.0 / z).conj() * output).re
                    - wave.outgoing_power_w * (dz / z);
                let physical_work =
                    (input.conj() * wave.outgoing + wave.physical_input.conj() * output).re;
                let residual = outward + physical_work;
                if !finite(input)
                    || [outward, physical_work, residual]
                        .iter()
                        .any(|v| !v.is_finite())
                {
                    return Err(invalid(
                        "loaded pole boundary derivative is not representable",
                    ));
                }
                outgoing_power += outward;
                physical_power += physical_work;
                waves.push(DecayWaveParameterDerivative {
                    original_port: p,
                    reference_derivative: dz,
                    outgoing: output,
                    physical_input: input,
                    outgoing_power_w: outward,
                    physical_power_w: physical_work,
                    power_residual_w: residual,
                });
            }
            let relative = physical
                .relative_residual
                .max(self.model.original_parameter_residual(
                    ParameterResidualState {
                        laplace: self.laplace_pole_per_s(),
                        laplace_derivative: physical.laplace_pole_per_s,
                        state: self.right_state(),
                        effort: self.mode_effort(),
                        input: &self.observation.physical_inputs,
                        state_derivative: &physical.state,
                        effort_derivative: &physical.effort,
                        input_derivative: &physical_inputs,
                    },
                    direction.physical,
                    &mut cancelled,
                )?);
            let energy_rate =
                2.0 * physical.laplace_pole_per_s.re * self.observation.stored_energy_j
                    + 2.0 * self.laplace_pole_per_s().re * physical.stored_energy_j;
            let total_decay = physical.dissipated_power_w + outgoing_power;
            let decay_balance = energy_rate + total_decay;
            let physical_balance = energy_rate + physical.dissipated_power_w - physical_power;
            if [
                outgoing_power,
                physical_power,
                total_decay,
                decay_balance,
                physical_balance,
            ]
            .iter()
            .any(|v| !v.is_finite())
            {
                return Err(invalid("loaded pole decay derivative is not representable"));
            }
            results.push(ScatteringPoleDerivative {
                laplace_pole_per_s: physical.laplace_pole_per_s,
                angular_frequency_rad_s: physical.angular_frequency_rad_s,
                state: physical.state,
                effort: physical.effort,
                output: physical.output,
                physical_inputs,
                waves,
                stored_energy_j: physical.stored_energy_j,
                internal_dissipated_power_w: physical.dissipated_power_w,
                outgoing_power_w: outgoing_power,
                total_decay_power_w: total_decay,
                decay_power_residual_w: decay_balance,
                physical_power_residual_w: physical_balance,
                relative_residual: relative,
                linear_residual: physical.linear_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(results)
    }
}
