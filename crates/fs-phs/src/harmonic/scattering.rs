//! Continuous wave response of the complete original physical system.
//!
//! At each selected power port, let y be the original output, u the original
//! input, and Z>0 an output/input reference. With peak exp(-i omega t) waves,
//! u=(2a-y)/Z and b=y-a. Close ALL selected ports in the original descriptor
//! pencil before factorization: L_closed=L-sum G_p G_p^T W/Z_p, then
//! A_closed=-i omega E-L_closed. This remains well defined at a bare component
//! resonance; no isolated impedance is inverted.
//!
//! The original R and all material/constraint coordinates remain unchanged.
//! The reference feedback is a boundary exchange, not internal material heat.
//! Frequency response has no time-step warping or interval-average substitution.

use super::{
    C64, HarmonicError, HarmonicFactor, HarmonicModel, HarmonicPortPower, HarmonicResponse, finite,
    invalid, poll,
};

mod frequency;
mod curvature;
mod mixed;
mod parameter;
mod poles;
mod resonance;
pub use frequency::{ScatteringFrequencyResponse, ScatteringWaveDerivative};
pub use curvature::{
    ScatteringComplexFrequencyJet, ScatteringFrequencyCurvatureResponse,
    ScatteringPowerFrequencyJet, ScatteringWaveFrequencyCurvature,
};
pub use mixed::{
    ScatteringComplexJet, ScatteringFrequencyParameterResponse, ScatteringPowerJet,
    ScatteringWaveFrequencyParameterJet,
};
pub use parameter::{ScatteringParameterResponse, ScatteringWaveParameterDerivative};
pub use poles::response::{
    DecayWaveParameterDerivative, PreparedScatteringPole, ScatteringPoleDerivative,
    ScatteringPoleObservation, ScatteringPoleParameterDirection, ScatteringPortResidue,
};
pub use poles::{DecayWavePort, ScatteringPole, ScatteringSpectrum};
pub use resonance::ScatteringObservation;

/// Additional resources for wave closure, beyond the consumed harmonic model's
/// existing budget. Payload excludes that model, allocator overhead and RSS.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringBudget {
    /// Every selected physical port is retained; no truncation to this cap.
    pub maximum_ports: usize,
    /// Peak added output rows, port records, response input and wave accounts.
    pub maximum_payload_bytes: usize,
    /// Bounds closure plus original-row observation: 2*n²*(selected_ports+1)
    /// + n*source_ports + selected_ports² + 32*source_ports.
    pub maximum_work: u64,
}

/// One independently evaluated boundary account; amplitudes use y's units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicWavePort {
    pub original_port: usize,
    /// Positive y/u reference. This is impedance when y=pressure and u=flow;
    /// it is admittance when y=velocity and u=force.
    pub reference: f64,
    pub incident: C64,
    pub outgoing: C64,
    /// |a|²/(2Z), in watts for power-conjugate source units.
    pub incident_power_w: f64,
    /// |b|²/(2Z), independently computed from the outgoing complex amplitude.
    pub outgoing_power_w: f64,
    /// Re(conj(u)*y)/2, positive INTO the original physical system.
    pub physical_power_w: f64,
    /// incident - outgoing - physical power; never used to manufacture heat.
    pub power_residual_w: f64,
}

/// Original state and physical losses, with actual boundary inputs restored.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringResponse {
    /// Every original state, output and material loss. Its source powers use
    /// actual u, including wave feedback, rather than the effective 2a/Z drive.
    /// Equilibrium linearizations retain their original incremental semantics.
    pub physical: HarmonicResponse,
    pub physical_inputs: Vec<C64>,
    /// Selected ports in their authored order; unselected ports are ordinary inputs.
    pub waves: Vec<HarmonicWavePort>,
    /// Original loss minus (incident-outgoing wave power + ordinary source work).
    /// This is independent of the original conjugate input/output power account.
    pub wave_power_residual_w: f64,
}

/// One immutable full-state closure, reusable across frequencies and drives.
pub struct ScatteringModel {
    closed: HarmonicModel,
    ports: Vec<(usize, f64)>,
    selected: Vec<bool>,
    budget: ScatteringBudget,
}

/// One shared LU of the fully closed physical problem, reusable across drives.
pub struct ScatteringFactor<'a> {
    model: &'a ScatteringModel,
    factor: HarmonicFactor<'a>,
}

impl ScatteringModel {
    /// Consume the complete admitted linear image and close the selected ports.
    /// All other source ports retain their original ordinary-input semantics.
    /// Supports descriptor constraints without inverting the descriptor mass.
    ///
    /// Z is an authored wave reference, not a newly inferred material or
    /// radiation law. Zero incident waves impose the corresponding matched
    /// boundary. Geometry, source bands and model validity remain the caller's
    /// source contract. No accepted time history is read or advanced.
    pub fn new(
        mut model: HarmonicModel,
        ports: &[(usize, f64)],
        budget: ScatteringBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HarmonicError> {
        poll(&mut cancelled)?;
        let (n, m, k) = (model.n, model.m, ports.len());
        let bytes = 8_u128 * n as u128 * k as u128 + 256_u128 * k as u128 + 64_u128 * m as u128;
        let work = 2_u128 * n as u128 * n as u128 * (k as u128 + 1)
            + n as u128 * m as u128
            + k as u128 * k as u128
            + 32_u128 * m as u128;
        if k == 0
            || k > budget.maximum_ports
            || bytes > budget.maximum_payload_bytes as u128
            || work > budget.maximum_work as u128
        {
            return Err(invalid(
                "complete scattering closure exceeds its explicit resources",
            ));
        }
        for (at, &(p, z)) in ports.iter().enumerate() {
            poll(&mut cancelled)?;
            if p >= m
                || !z.is_finite()
                || z <= 0.0
                || !(2.0 / z).is_finite()
                || 1.0 / z == 0.0
                || ports[..at].iter().any(|&(q, _)| q == p)
            {
                return Err(invalid(
                    "distinct original ports and positive representable references required",
                ));
            }
        }
        // Fixed original output rows G^T W. Algebraic multipliers are efforts,
        // so their W block is identity, never an artificial energy matrix.
        let mut output = vec![0.0; k * n];
        for (at, &(p, _)) in ports.iter().enumerate() {
            for j in 0..n {
                poll(&mut cancelled)?;
                output[at * n + j] = if j < model.n_diff {
                    (0..model.n_diff)
                        .map(|i| model.g[i * m + p] * model.q[i * n + j])
                        .sum()
                } else {
                    model.g[j * m + p]
                };
            }
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "original scattering output row is not representable",
            ));
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                for (at, &(p, z)) in ports.iter().enumerate() {
                    let g = model.g[i * m + p];
                    let y = output[at * n + j];
                    let feedback = g * (y / z);
                    if !feedback.is_finite() || (g != 0.0 && y != 0.0 && feedback == 0.0) {
                        return Err(invalid("scattering feedback is not representable"));
                    }
                    model.generator[i * n + j] -= feedback;
                }
                if !model.generator[i * n + j].is_finite() {
                    return Err(invalid("closed scattering operator is not representable"));
                }
            }
        }
        let mut selected = vec![false; m];
        for &(p, _) in ports {
            selected[p] = true;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            closed: model,
            ports: ports.to_vec(),
            selected,
            budget,
        })
    }

    pub fn dimension(&self) -> usize {
        self.closed.n
    }
    /// All original source slots, including ordinary inputs.
    pub fn port_count(&self) -> usize {
        self.closed.m
    }
    /// Original source slot and positive output/input reference, in wave order.
    pub fn wave_ports(&self) -> &[(usize, f64)] {
        &self.ports
    }

    pub fn factor(
        &self,
        omega_rad_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringFactor<'_>, HarmonicError> {
        let factor = self.closed.factor(omega_rad_s, cancelled)?;
        Ok(ScatteringFactor {
            model: self,
            factor,
        })
    }

    // Recheck the original physical equations with actual boundary inputs.
    // Both forced response and free decay use this same independent account;
    // the closed generator is deliberately not used to validate its own sign.
    fn original_residual(
        &self,
        laplace: C64,
        state: &[C64],
        effort: &[C64],
        inputs: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, HarmonicError> {
        let original = &self.closed;
        let mut relative = 0.0_f64;
        for i in 0..original.n {
            poll(cancelled)?;
            let mut residual = if i < original.n_diff {
                laplace * state[i]
            } else {
                C64::ZERO
            };
            let mut scale = residual.abs();
            for (j, &e) in effort.iter().enumerate() {
                let term = e.scale(original.structure[i * original.n + j]);
                residual = residual - term;
                scale += term.abs();
            }
            if original.resisting_tangent.is_some() {
                let term = original.nonlinear_loss_row(i, state);
                residual = residual + term;
                scale += term.abs();
            }
            for (p, u) in inputs.iter().enumerate() {
                let term = u.scale(original.g[i * original.m + p]);
                residual = residual - term;
                scale += term.abs();
            }
            if !finite(residual) || !scale.is_finite() {
                return Err(invalid(
                    "original scattering equation residual is not representable",
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

impl ScatteringFactor<'_> {
    /// One finite peak amplitude per ORIGINAL port: incident output wave a at
    /// selected ports, ordinary input u elsewhere. Multiple drives share this
    /// factorization. All response storage is private until acceptance.
    pub fn response(
        &self,
        input: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScatteringResponse, HarmonicError> {
        let effective = self.effective_input(input, &mut cancelled)?;
        let physical = self.factor.response(&effective, &mut cancelled)?;
        self.observe_response(input, physical, &mut cancelled)
    }

    // The same source convention applies to amplitudes and their derivatives:
    // incident waves on selected slots, ordinary physical inputs elsewhere.
    fn effective_input(
        &self,
        input: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HarmonicError> {
        poll(cancelled)?;
        if input.len() != self.model.closed.m || input.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "one finite wave or ordinary input per original port required",
            ));
        }
        let mut effective = input.to_vec();
        for &(p, z) in &self.model.ports {
            effective[p] = input[p].scale(2.0 / z);
        }
        Ok(effective)
    }

    // Restore original boundary inputs and physical work for either an ordinary
    // response or the base of its frequency tangent. No second solve is needed.
    fn observe_response(
        &self,
        input: &[C64],
        mut physical: HarmonicResponse,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ScatteringResponse, HarmonicError> {
        let mut physical_inputs = input.to_vec();
        let mut waves = Vec::with_capacity(self.model.ports.len());
        for &(p, z) in &self.model.ports {
            poll(cancelled)?;
            let a = input[p];
            let b = physical.output[p] - a;
            let u = (a.scale(2.0) - physical.output[p]).scale(1.0 / z);
            physical_inputs[p] = u;
            let incident = 0.5 * a.norm_sq() / z;
            let outgoing = 0.5 * b.norm_sq() / z;
            let power = HarmonicPortPower::from_peak(u, physical.output[p])?.active_w;
            let residual = incident - outgoing - power;
            if !finite(b)
                || !finite(u)
                || [incident, outgoing, power, residual]
                    .iter()
                    .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "physical scattering wave or power is not representable",
                ));
            }
            waves.push(HarmonicWavePort {
                original_port: p,
                reference: z,
                incident: a,
                outgoing: b,
                incident_power_w: incident,
                outgoing_power_w: outgoing,
                physical_power_w: power,
                power_residual_w: residual,
            });
        }
        for (p, &u) in physical_inputs.iter().enumerate() {
            poll(cancelled)?;
            let power = HarmonicPortPower::from_peak(u, physical.output[p])?;
            physical.port_power_w[p] = power.active_w;
            physical.port_power_components[p] = power;
        }
        physical.supplied_power_w = super::physical_power(
            &self.model.closed, &physical_inputs, &physical.effort, cancelled,
        )?;
        physical.power_residual_w = physical.dissipated_power_w - physical.supplied_power_w;
        let mut wave_supply = super::physical_power_sum(&self.model.closed,
            &physical_inputs, &physical.effort, |p| !self.model.selected[p], cancelled)?;
        for wave in &waves {
            poll(cancelled)?;
            wave_supply.product(wave.incident_power_w, 1.0)
                .and_then(|()| wave_supply.product(-wave.outgoing_power_w, 1.0))
                .ok_or_else(|| invalid("physical wave supply is not representable"))?;
        }
        let wave_power_residual_w = physical.dissipated_power_w
            - wave_supply.value().ok_or_else(|| invalid("physical wave supply is not representable"))?;
        if physical.port_power_w.iter().any(|v| !v.is_finite())
            || [
                physical.supplied_power_w,
                physical.power_residual_w,
                wave_power_residual_w,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid("physical scattering balance is not representable"));
        }
        physical.relative_residual = physical.relative_residual.max(self.model.original_residual(
            C64::new(0.0, -physical.omega_rad_s),
            &physical.state,
            &physical.effort,
            &physical_inputs,
            cancelled,
        )?);
        poll(cancelled)?;
        Ok(ScatteringResponse {
            physical,
            physical_inputs,
            waves,
            wave_power_residual_w,
        })
    }
}
