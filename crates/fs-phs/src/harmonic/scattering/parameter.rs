//! Physical parameter tangents with the entire wave boundary closed.
//!
//! Write dy=Y dx+d, d=dG^T e+G^T dQ x. Differentiating u=(2a-y)/Z
//! gives du=(2da-dy-u*dZ)/Z. The SAME closed factor therefore solves
//! A_closed dx=dL x+dG*u+G*(2da-d-u*dZ)/Z at selected ports, with
//! ordinary source tangents elsewhere. Both appearances of G and the change
//! of the physical observation map remain; boundary loading is not R/heat.

use super::{
    C64, HarmonicError, HarmonicPortPower, ScatteringFactor, ScatteringModel, ScatteringResponse,
    finite, invalid, poll,
};
use crate::harmonic::{HarmonicParameterDirection, parameter};
use parameter::{ParameterStateDerivative, accumulate, at};

/// One selected boundary differentiated per unit of the caller's real parameter.
/// The reference derivative follows the same authored order as wave_ports().
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringWaveParameterDerivative {
    pub original_port: usize,
    pub reference_derivative: f64,
    pub incident: C64,
    pub outgoing: C64,
    /// d(|a|²/(2Z)), including -P_incident*dZ/Z.
    pub incident_power_w: f64,
    /// d(|b|²/(2Z)), independently including -P_outgoing*dZ/Z.
    pub outgoing_power_w: f64,
    pub physical_power_w: f64,
    pub power_residual_w: f64,
}

/// Complete loaded response and a local physical parameter derivative.
/// All slopes are per caller-declared parameter unit, at fixed angular frequency,
/// state chart, topology and descriptor rank. Signed heat slopes are not losses.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringParameterResponse {
    pub response: ScatteringResponse,
    pub state_derivative: Vec<C64>,
    pub effort_derivative: Vec<C64>,
    pub output_derivative: Vec<C64>,
    /// Actual original du, including feedback, changing G/Q, incident wave and Z.
    pub physical_input_derivative: Vec<C64>,
    pub wave_derivatives: Vec<ScatteringWaveParameterDerivative>,
    pub mean_stored_energy_derivative: f64,
    pub port_power_derivative: Vec<f64>,
    pub supplied_power_derivative: f64,
    /// Independent original R/T material work. Boundary feedback is excluded.
    pub dissipated_power_derivative: f64,
    pub power_residual_derivative: f64,
    pub wave_power_residual_derivative: f64,
    /// Largest closed tangent and independently audited original-row residual.
    /// This is a backward check, not a near-resonance forward-error certificate.
    pub relative_derivative_residual: f64,
}

impl ScatteringFactor<'_> {
    /// Differentiate physical coefficients, source phasors and wave references
    /// through the complete loaded response, sharing its existing LU.
    ///
    /// `direction` describes the ORIGINAL dJ/dR/dQ/dG/dT, including every material
    /// and operating-point term. `input_derivative` has one original source slot:
    /// incident output-wave derivative on selected ports, ordinary du elsewhere.
    /// `reference_derivative` is dZ in wave_ports() order, NOT source-port order.
    /// None holds the respective values fixed. Finite signed dZ is a local
    /// direction at an already admitted positive reference, not a finite update.
    ///
    /// Full constraints, mutual storage and memory remain. Changing G affects
    /// forcing, feedback AND observation. Changing Z also changes incident and
    /// outgoing watt accounts at fixed output-wave amplitude. This remains
    /// valid at a bare pole when the closed operator is regular. No impedance
    /// inverse, perturbed model, derivative matrix or finite difference is used.
    ///
    /// The original physical parameter-work/payload admission applies. Added
    /// scattering payload is twice its constructor bound plus 96n+128m+128k
    /// bytes, and work is twice its constructor score plus 16n(n+m)+64(m+k),
    /// where k is the selected port count. All are checked before allocation.
    /// No bias solve, moving chart/rank or finite parameter interval is inferred.
    pub fn response_with_parameter_derivative(
        &self,
        input: &[C64],
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringParameterResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let model = &self.model.closed;
        let (n, m, k) = (model.n, model.m, self.model.ports.len());
        if input_derivative.is_some_and(|v| v.len() != m || v.iter().any(|&z| !finite(z)))
            || reference_derivative
                .is_some_and(|v| v.len() != k || v.iter().any(|z| !z.is_finite()))
        {
            return Err(invalid(
                "complete finite source and selected-reference derivatives required",
            ));
        }
        parameter::admit_response(model)?;
        let (nn, mm, kk) = (n as u128, m as u128, k as u128);
        let bytes = 2 * (8 * nn * kk + 256 * kk + 64 * mm) + 96 * nn + 128 * mm + 128 * kk;
        let work = 2 * (2 * nn * nn * (kk + 1) + nn * mm + kk * kk + 32 * mm)
            + 16 * nn * (nn + mm)
            + 64 * (mm + kk);
        if bytes > self.model.budget.maximum_payload_bytes as u128
            || work > self.model.budget.maximum_work as u128
        {
            return Err(invalid(
                "insufficient loaded physical parameter observation resources",
            ));
        }
        direction.validate(model, &mut cancelled)?;
        let response = self.response(input, &mut cancelled)?;
        let action = parameter::generator_action(
            model,
            &direction,
            &response.physical.state,
            &response.physical.effort,
            &mut cancelled,
        )?;
        let mut rhs = action.generator;
        add_parameter_source(
            self.model,
            direction,
            &response.physical.effort,
            &response.physical_inputs,
            &action.storage,
            input_derivative,
            reference_derivative,
            &mut rhs,
            &mut cancelled,
        )?;
        let (state, relative) = self.factor.solve_rhs(&rhs, &mut cancelled)?;
        let mut observed = parameter::observe_state(
            model,
            &response.physical,
            direction,
            state,
            &action.storage,
            relative,
            &mut cancelled,
        )?;
        let mut physical_input_derivative =
            input_derivative.map_or_else(|| vec![C64::ZERO; m], <[C64]>::to_vec);
        for (index, &(p, z)) in self.model.ports.iter().enumerate() {
            let dz = reference_derivative.map_or(0.0, |v| v[index]);
            physical_input_derivative[p] =
                (physical_input_derivative[p].scale(2.0) - observed.output[p]).scale(1.0 / z)
                    - response.physical_inputs[p].scale(dz / z);
        }
        if physical_input_derivative.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "actual loaded parameter input is not representable",
            ));
        }
        observed.relative_residual =
            observed
                .relative_residual
                .max(self.original_parameter_residual(
                    &response,
                    &observed,
                    &physical_input_derivative,
                    direction,
                    &mut cancelled,
                )?);
        let ScatteringResponse {
            physical,
            physical_inputs,
            waves,
            wave_power_residual_w,
        } = response;
        let physical = parameter::finish_response(
            physical,
            observed,
            &physical_inputs,
            Some(&physical_input_derivative),
            &mut cancelled,
        )?;
        let mut wave_derivatives = Vec::with_capacity(k);
        let mut wave_supply = 0.0;
        for (index, wave) in waves.iter().enumerate() {
            poll(&mut cancelled)?;
            let p = wave.original_port;
            let dz = reference_derivative.map_or(0.0, |v| v[index]);
            let incident = input_derivative.map_or(C64::ZERO, |v| v[p]);
            let outgoing = physical.output_derivative[p] - incident;
            let incident_power_w = 2.0
                * half_pair(wave.incident.scale(1.0 / wave.reference), incident)?
                - wave.incident_power_w * (dz / wave.reference);
            let outgoing_power_w = 2.0
                * half_pair(wave.outgoing.scale(1.0 / wave.reference), outgoing)?
                - wave.outgoing_power_w * (dz / wave.reference);
            let physical_power_w = physical.port_power_derivative[p];
            let power_residual_w = incident_power_w - outgoing_power_w - physical_power_w;
            if !finite(outgoing)
                || [
                    incident_power_w,
                    outgoing_power_w,
                    physical_power_w,
                    power_residual_w,
                ]
                .iter()
                .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "parameter wave or independent boundary power overflows",
                ));
            }
            wave_supply += incident_power_w - outgoing_power_w;
            wave_derivatives.push(ScatteringWaveParameterDerivative {
                original_port: p,
                reference_derivative: dz,
                incident,
                outgoing,
                incident_power_w,
                outgoing_power_w,
                physical_power_w,
                power_residual_w,
            });
        }
        for (p, &work) in physical.port_power_derivative.iter().enumerate() {
            if !self.model.selected[p] {
                wave_supply += work;
            }
        }
        let wave_power_residual_derivative = physical.dissipated_power_derivative - wave_supply;
        if !wave_supply.is_finite() || !wave_power_residual_derivative.is_finite() {
            return Err(invalid("loaded parameter wave balance overflows"));
        }
        poll(&mut cancelled)?;
        Ok(ScatteringParameterResponse {
            response: ScatteringResponse {
                physical: physical.response,
                physical_inputs,
                waves,
                wave_power_residual_w,
            },
            state_derivative: physical.state_derivative,
            effort_derivative: physical.effort_derivative,
            output_derivative: physical.output_derivative,
            physical_input_derivative,
            wave_derivatives,
            mean_stored_energy_derivative: physical.mean_stored_energy_derivative,
            port_power_derivative: physical.port_power_derivative,
            supplied_power_derivative: physical.supplied_power_derivative,
            dissipated_power_derivative: physical.dissipated_power_derivative,
            power_residual_derivative: physical.power_residual_derivative,
            wave_power_residual_derivative,
            relative_derivative_residual: physical.relative_derivative_residual,
        })
    }

    // Independently audit the differentiated ORIGINAL physical equations:
    // sE dx-(J-R)de-(dJ-dR)e+T dx+dT x-G du-dG u=0. No closed
    // generator or imposed wave-power balance can validate its own sign here.
    fn original_parameter_residual(
        &self,
        response: &ScatteringResponse,
        observed: &ParameterStateDerivative,
        input_derivative: &[C64],
        direction: HarmonicParameterDirection<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, HarmonicError> {
        let model = &self.model.closed;
        let (n, m) = (model.n, model.m);
        let s = C64::new(0.0, -response.physical.omega_rad_s);
        let mut relative = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut residual = if i < model.n_diff {
                s * observed.state[i]
            } else {
                C64::ZERO
            };
            let mut scale = residual.abs();
            for j in 0..n {
                let terms = [
                    -observed.effort[j].scale(model.structure[i * n + j]),
                    -response.physical.effort[j].scale(
                        at(direction.interconnection, i * n + j)
                            - at(direction.dissipation, i * n + j),
                    ),
                    observed.state[j].scale(
                        model
                            .resisting_tangent
                            .as_ref()
                            .map_or(0.0, |t| t[i * n + j]),
                    ),
                    response.physical.state[j].scale(at(direction.resisting_flow, i * n + j)),
                ];
                for term in terms {
                    residual = residual + term;
                    scale += term.abs();
                }
            }
            for p in 0..m {
                let force = input_derivative[p].scale(model.g[i * m + p]);
                let map = response.physical_inputs[p].scale(at(direction.ports, i * m + p));
                residual = residual - force - map;
                scale += force.abs() + map.abs();
            }
            if !finite(residual) || !scale.is_finite() {
                return Err(invalid(
                    "original loaded parameter row is not representable",
                ));
            }
            relative = relative.max(if scale == 0.0 {
                0.0
            } else {
                residual.abs() / scale
            });
        }
        if relative > model.budget.relative_residual_tolerance {
            return Err(HarmonicError::Residual {
                relative,
                tolerance: model.budget.relative_residual_tolerance,
            });
        }
        Ok(relative)
    }
}

// Add dG*u+G*du_explicit to an admitted full-state physical action. The explicit
// changing output map excludes Y*dx, already closed on the left. Reusing this
// same action on frequency tangents supplies the mixed frequency/parameter RHS.
#[allow(clippy::too_many_arguments)]
pub(super) fn add_parameter_source(
    owner: &ScatteringModel,
    direction: HarmonicParameterDirection<'_>,
    effort: &[C64],
    physical_input: &[C64],
    storage_action: &[C64],
    input_derivative: Option<&[C64]>,
    reference_derivative: Option<&[f64]>,
    rhs: &mut [C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HarmonicError> {
    let model = &owner.closed;
    let (n, m) = (model.n, model.m);
    let mut direct_output = vec![C64::ZERO; m];
    for p in 0..m {
        poll(cancelled)?;
        let mut correction = C64::ZERO;
        for i in 0..n {
            accumulate(
                &mut direct_output[p],
                &mut correction,
                effort[i].scale(at(direction.ports, i * m + p))
                    + storage_action[i].scale(model.g[i * m + p]),
            );
        }
    }
    let mut source = input_derivative.map_or_else(|| vec![C64::ZERO; m], <[C64]>::to_vec);
    for (index, &(p, z)) in owner.ports.iter().enumerate() {
        poll(cancelled)?;
        let dz = reference_derivative.map_or(0.0, |v| v[index]);
        source[p] = (source[p].scale(2.0) - direct_output[p]).scale(1.0 / z)
            - physical_input[p].scale(dz / z);
    }
    if direct_output.iter().chain(&source).any(|&z| !finite(z)) {
        return Err(invalid(
            "loaded parameter source action is not representable",
        ));
    }
    for i in 0..n {
        poll(cancelled)?;
        let mut correction = C64::ZERO;
        for p in 0..m {
            accumulate(
                &mut rhs[i],
                &mut correction,
                physical_input[p].scale(at(direction.ports, i * m + p))
                    + source[p].scale(model.g[i * m + p]),
            );
        }
    }
    Ok(())
}

fn half_pair(left: C64, right: C64) -> Result<f64, HarmonicError> {
    Ok(HarmonicPortPower::from_peak(left, right)?.active_w)
}
