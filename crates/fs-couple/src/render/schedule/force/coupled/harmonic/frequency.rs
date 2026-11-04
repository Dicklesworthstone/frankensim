//! Analytic frequency tangents of the original mechanical/material/air pencil.
//!
//! Differentiate D(omega) q = F as D q' = F' - D' q, then reuse the SAME
//! condensed factor. Maxwell histories remain causal and air remains reciprocal;
//! neither is replaced with a frozen complex stiffness or an isolated impedance.

use super::*;

/// The original AC response and its derivative with respect to angular frequency.
///
/// Every derivative is per rad/s, in the original field's physical units and
/// coordinate/channel order. The physical coefficients and DC equilibrium stay
/// fixed. An optional force derivative describes a frequency-dependent drive.
/// Signed derivatives of heat/storage are slopes, not negative physical losses.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalHarmonicFrequencyResponse {
    pub response: ModalHarmonicResponse,
    pub displacement_derivative: Vec<C64>,
    pub velocity_derivative: Vec<C64>,
    pub material_memory_derivative: Vec<C64>,
    pub radiation_state_derivative: Vec<C64>,
    pub connection_forces_derivative: Vec<C64>,
    pub connection_spring_forces_derivative: Vec<C64>,
    pub connection_damping_forces_derivative: Vec<C64>,
    pub generalized_connection_forces_derivative: Vec<C64>,
    pub mean_stored_energy_derivative: f64,
    pub external_power_derivative: f64,
    pub component_dissipation_derivative: f64,
    pub connection_dissipation_derivative: f64,
    pub material_dissipation_derivative: f64,
    pub material_branch_dissipation_derivative: Vec<f64>,
    pub material_branch_mean_storage_derivative: Vec<f64>,
    pub radiation_dissipation_derivative: f64,
    /// Derivative of total mean loss minus external mean power [W per rad/s].
    pub power_residual_derivative: f64,
    /// Infinity-norm backward residual of D q' = F' - D' q. This is not a
    /// forward-error certificate for an ill-conditioned near-pole response.
    pub relative_derivative_residual: f64,
}

impl ModalHarmonicFrequencyResponse {
    pub fn total_dissipation_derivative(&self) -> f64 {
        self.component_dissipation_derivative
            + self.connection_dissipation_derivative
            + self.material_dissipation_derivative
            + self.radiation_dissipation_derivative
    }

    /// Derivative of log amplitude (real part) and phase (imaginary part) of a
    /// FIXED real projection of the original displacements. With exp(-i omega t)
    /// a pure delay has positive phase slope in seconds. A zero observation has
    /// no phase/log amplitude and returns None; there is no amplitude floor.
    /// Includes any supplied drive slope, so it is a transfer delay only for a
    /// frequency-independent force phasor.
    pub fn displacement_logarithmic_derivative(
        &self,
        observation: &[f64],
    ) -> Result<Option<C64>, ModalCouplingError> {
        if observation.len() != self.response.displacement_m_sqrt_kg.len()
            || observation.len() != self.displacement_derivative.len()
            || observation.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "frequency observation must match every original finite mode",
            ));
        }
        let value = projection(observation, &self.response.displacement_m_sqrt_kg);
        let derivative = projection(observation, &self.displacement_derivative);
        if !complex_finite(value) || !complex_finite(derivative) {
            return Err(invalid("frequency observation is not representable"));
        }
        if value == C64::ZERO {
            return Ok(None);
        }
        let slope = derivative / value;
        if !complex_finite(slope) {
            return Err(invalid("frequency observation slope is not representable"));
        }
        Ok(Some(slope))
    }
}

fn pole_derivative(omega: f64, damping: f64) -> C64 {
    C64::new(-2.0 * omega, -damping)
}

fn velocity_derivative(omega: f64, displacement: C64, derivative: C64) -> C64 {
    C64::new(0.0, -omega) * derivative + C64::new(0.0, -1.0) * displacement
}

// g=s/(s+r), h=r/(s+r)=1-g. Evaluate g'=-i r/(s+r)^2 as
// two bounded ratios, avoiding a squared denominator and r*omega products.
fn relaxation_gain_derivative(omega: f64, rate: f64) -> C64 {
    let denominator = C64::new(rate, -omega);
    (C64::new(0.0, -1.0) / denominator) * (C64::from_re(rate) / denominator)
}

// The derivative of |x|^2/2; all mean-energy and mean-power factors below
// follow from this bilinear form, including each original Maxwell channel.
fn half_norm_derivative(value: C64, derivative: C64) -> f64 {
    (value.conj() * derivative).re
}

fn subtract_rank_action(rhs: &mut [C64], column: &[f64], state: &[C64], gain_derivative: C64) {
    let effort = gain_derivative * projection(column, state);
    for (value, &coefficient) in rhs.iter_mut().zip(column) {
        *value = *value - effort.scale(coefficient);
    }
}

fn subtract_reciprocal_air_action(
    rhs: &mut [C64],
    state: &[C64],
    mechanical: usize,
    acoustic: usize,
    coupling: f64,
) {
    let derivative = C64::new(0.0, -coupling);
    rhs[mechanical] = rhs[mechanical] - derivative * state[acoustic];
    rhs[acoustic] = rhs[acoustic] + derivative * state[mechanical];
}

impl CoupledModalSystem {
    /// Original coupled AC solution and its analytic frequency derivative.
    ///
    /// Uses the same physical scope and phasor convention as `harmonic_response`:
    /// all original mechanical coordinates, connection directions, Maxwell
    /// histories and two-way air poles; fixed coefficients and equilibrium.
    /// Nonlinear contacts and external DG terminals have their own response APIs
    /// and are not part of this particular linear response. Material/source fit
    /// frequency validity remains the coefficient producer's responsibility.
    ///
    /// `force_derivative` is dF/domega in the original generalized-force order;
    /// None holds the force fixed. A second solve reuses the original LU without
    /// another factorization or a dense derivative matrix. Maxwell gain slopes
    /// and BOTH air feedback directions are differentiated before that solve.
    /// Every physical observation includes its explicit frequency dependence.
    ///
    /// For d=mechanical_modes+air_poles, C links, G anchored-loss links and M
    /// Maxwell branches, `maximum_work` bounds d²*(d+G+2*(C+M+6)). Complete
    /// response/tangent storage is bounded by MAX_HARMONIC_BATCH_RESPONSE_BYTES
    /// before assembly. Each solve must meet the existing residual tolerance.
    /// No accepted state, material history, sample clock or ledger is changed.
    pub fn harmonic_response_with_frequency_derivative(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ModalHarmonicFrequencyResponse, ModalCouplingError> {
        poll(Some(gate))?;
        let n = self.mode_count();
        let poles = self
            .radiation
            .as_ref()
            .map_or(0, |air| air.model.poles.len());
        let dimension = n
            .checked_add(poles)
            .ok_or_else(|| invalid("harmonic tangent dimension overflow"))?;
        if force_derivative.is_some_and(|values| {
            values.len() != n || values.iter().any(|&value| !complex_finite(value))
        }) {
            return Err(invalid(
                "harmonic force derivative must match every original finite mode",
            ));
        }
        let work = harmonic_work(
            dimension,
            self.connections.len(),
            self.connections
                .iter()
                .filter(|law| law.constant_loss_reference_omega_rad_s.is_some())
                .count(),
            self.relaxation_branches.len(),
            6,
            2,
        )?;
        let bytes = harmonic_response_bytes(
            n,
            self.connections.len(),
            self.relaxation_branches.len(),
            poles,
            2,
        )?
        .checked_add(core::mem::size_of::<ModalHarmonicFrequencyResponse>())
        .ok_or_else(|| invalid("harmonic tangent storage overflow"))?;
        if dimension > MAX_HARMONIC_UNKNOWNS
            || work > maximum_work
            || bytes > MAX_HARMONIC_BATCH_RESPONSE_BYTES
        {
            return Err(invalid(
                "complete harmonic tangent exceeds its work or retained-response storage bound",
            ));
        }
        let operator = self.harmonic_operator(omega_rad_s, force, 0, maximum_work, gate)?;
        let factor = HarmonicFactor::new(&operator.values, operator.dimension, gate)?;
        let (state, relative) = factor.solve(force, self.config.solve_relative_tolerance, gate)?;
        let response = self.observe_harmonic_state(omega_rad_s, force, &state, relative, gate)?;
        let rhs = self.harmonic_frequency_rhs(omega_rad_s, &state, force_derivative, gate)?;
        let (derivative, relative_derivative) =
            factor.solve(&rhs, self.config.solve_relative_tolerance, gate)?;
        self.observe_harmonic_frequency_derivative(
            response,
            force,
            force_derivative,
            &state,
            &derivative,
            relative_derivative,
            gate,
        )
    }

    pub(super) fn harmonic_frequency_rhs(
        &self,
        omega: f64,
        state: &[C64],
        force_derivative: Option<&[C64]>,
        gate: &CancelGate,
    ) -> Result<Vec<C64>, ModalCouplingError> {
        let n = self.mode_count();
        let mut rhs = vec![C64::ZERO; state.len()];
        if let Some(force) = force_derivative {
            rhs[..n].copy_from_slice(force);
        }
        for (index, &damping) in self
            .models
            .iter()
            .flat_map(|model| model.damping_rates_per_s())
            .enumerate()
        {
            poll(Some(gate))?;
            rhs[index] = rhs[index]
                - pole_derivative(omega, damping) * state[index];
        }
        for (index, law) in self.connections.iter().enumerate() {
            poll(Some(gate))?;
            subtract_rank_action(
                &mut rhs,
                self.harmonic_damping_column(index)?,
                &state[..n],
                C64::new(0.0, -law.damping_n_s_m),
            );
        }
        let mut column = vec![0.0; n];
        for arm in &self.relaxation_branches {
            poll(Some(gate))?;
            let root = det::sqrt(arm.stiffness);
            for (value, &p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) {
                *value = root * p;
            }
            subtract_rank_action(
                &mut rhs,
                &column,
                &state[..n],
                relaxation_gain_derivative(omega, arm.relaxation_time_s.recip()),
            );
        }
        if let Some(air) = &self.radiation {
            for (j, pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let index = n + j;
                rhs[index] = rhs[index]
                    - pole_derivative(omega, 2.0 * pole.zeta * pole.omega) * state[index];
                for (&i, &coupling) in air.source_modes.iter().zip(&pole.coupling) {
                    subtract_reciprocal_air_action(&mut rhs, state, i, index, coupling);
                }
            }
        }
        if rhs.iter().any(|&value| !complex_finite(value)) {
            return Err(invalid(
                "harmonic frequency tangent forcing is not representable",
            ));
        }
        poll(Some(gate))?;
        Ok(rhs)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn observe_harmonic_frequency_derivative(
        &self,
        response: ModalHarmonicResponse,
        force: &[C64],
        force_derivative: Option<&[C64]>,
        state: &[C64],
        derivative: &[C64],
        relative_derivative: f64,
        gate: &CancelGate,
    ) -> Result<ModalHarmonicFrequencyResponse, ModalCouplingError> {
        let n = self.mode_count();
        let omega = response.omega_rad_s;
        let s = C64::new(0.0, -omega);
        let mut result = ModalHarmonicFrequencyResponse {
            displacement_derivative: derivative[..n].to_vec(),
            velocity_derivative: state[..n]
                .iter()
                .zip(&derivative[..n])
                .map(|(&q, &dq)| velocity_derivative(omega, q, dq))
                .collect(),
            material_memory_derivative: Vec::with_capacity(response.material_memory.len()),
            radiation_state_derivative: Vec::with_capacity(response.radiation_state.len()),
            connection_forces_derivative: Vec::with_capacity(self.connections.len()),
            connection_spring_forces_derivative: Vec::with_capacity(self.connections.len()),
            connection_damping_forces_derivative: Vec::with_capacity(self.connections.len()),
            generalized_connection_forces_derivative: vec![C64::ZERO; n],
            mean_stored_energy_derivative: 0.0,
            external_power_derivative: 0.0,
            component_dissipation_derivative: 0.0,
            connection_dissipation_derivative: 0.0,
            material_dissipation_derivative: 0.0,
            material_branch_dissipation_derivative: Vec::with_capacity(
                self.relaxation_branches.len(),
            ),
            material_branch_mean_storage_derivative: Vec::with_capacity(
                self.relaxation_branches.len(),
            ),
            radiation_dissipation_derivative: 0.0,
            power_residual_derivative: 0.0,
            relative_derivative_residual: relative_derivative,
            response,
        };
        for (index, (mode, &damping)) in self
            .models
            .iter()
            .flat_map(|model| model.modes().iter().zip(model.damping_rates_per_s()))
            .enumerate()
        {
            poll(Some(gate))?;
            let q = state[index];
            let dq = derivative[index];
            let v = result.response.velocity_m_sqrt_kg_per_s[index];
            let dv = result.velocity_derivative[index];
            let w = mode.angular_frequency_rad_s;
            result.mean_stored_energy_derivative +=
                0.5 * (half_norm_derivative(q.scale(w), dq.scale(w)) + half_norm_derivative(v, dv));
            result.component_dissipation_derivative +=
                damping * half_norm_derivative(v, dv);
            let df = force_derivative.map_or(C64::ZERO, |values| values[index]);
            result.external_power_derivative += 0.5 * (df * v.conj() + force[index] * dv.conj()).re;
        }
        for (index, (column, law)) in self.columns.iter().zip(&self.connections).enumerate() {
            poll(Some(gate))?;
            let damping_column = self.harmonic_damping_column(index)?;
            let extension = projection(column, &state[..n]);
            let d_extension = projection(column, &derivative[..n]);
            let damping_extension = projection(damping_column, &state[..n]);
            let velocity = s * damping_extension;
            let d_velocity = velocity_derivative(
                omega,
                damping_extension,
                projection(damping_column, &derivative[..n]),
            );
            let spring = -d_extension.scale(law.stiffness_n_m);
            let damping = -d_velocity.scale(law.damping_n_s_m);
            result.connection_forces_derivative.push(spring + damping);
            result.connection_spring_forces_derivative.push(spring);
            result.connection_damping_forces_derivative.push(damping);
            for ((value, &b), &d) in result
                .generalized_connection_forces_derivative
                .iter_mut()
                .zip(column)
                .zip(damping_column)
            {
                *value = *value + spring.scale(b) + damping.scale(d);
            }
            result.mean_stored_energy_derivative +=
                0.5 * law.stiffness_n_m * half_norm_derivative(extension, d_extension);
            result.connection_dissipation_derivative +=
                law.damping_n_s_m * half_norm_derivative(velocity, d_velocity);
        }
        let mut column = vec![0.0; n];
        for arm in &self.relaxation_branches {
            poll(Some(gate))?;
            let root = det::sqrt(arm.stiffness);
            let rate = arm.relaxation_time_s.recip();
            for (value, &p) in column.iter_mut().zip(arm.projection.iter().step_by(2)) {
                *value = root * p;
            }
            let strain = projection(&column, &state[..n]);
            let d_strain = projection(&column, &derivative[..n]);
            let denominator = s + C64::from_re(rate);
            let gain = s / denominator;
            let gain_derivative = relaxation_gain_derivative(omega, rate);
            let memory_derivative =
                (C64::from_re(rate) / denominator) * d_strain - gain_derivative * strain;
            let mismatch = gain * strain;
            let mismatch_derivative = gain * d_strain + gain_derivative * strain;
            let stored = 0.5 * half_norm_derivative(mismatch, mismatch_derivative);
            let heat = rate * half_norm_derivative(mismatch, mismatch_derivative);
            result.material_memory_derivative.push(memory_derivative);
            result.mean_stored_energy_derivative += stored;
            result.material_dissipation_derivative += heat;
            result.material_branch_mean_storage_derivative.push(stored);
            result.material_branch_dissipation_derivative.push(heat);
        }
        if let Some(air) = &self.radiation {
            for (j, pole) in air.model.poles.iter().enumerate() {
                poll(Some(gate))?;
                let q = state[n + j];
                let dq = derivative[n + j];
                let p = s * q;
                let dp = velocity_derivative(omega, q, dq);
                result.radiation_state_derivative.extend([dq, dp]);
                result.mean_stored_energy_derivative += 0.5
                    * (half_norm_derivative(q.scale(pole.omega), dq.scale(pole.omega))
                        + half_norm_derivative(p, dp));
                result.radiation_dissipation_derivative +=
                    2.0 * pole.zeta * pole.omega * half_norm_derivative(p, dp);
            }
        }
        result.power_residual_derivative =
            result.total_dissipation_derivative() - result.external_power_derivative;
        if result
            .displacement_derivative
            .iter()
            .chain(&result.velocity_derivative)
            .chain(&result.material_memory_derivative)
            .chain(&result.radiation_state_derivative)
            .chain(&result.connection_forces_derivative)
            .chain(&result.connection_spring_forces_derivative)
            .chain(&result.connection_damping_forces_derivative)
            .chain(&result.generalized_connection_forces_derivative)
            .any(|&value| !complex_finite(value))
            || result
                .material_branch_dissipation_derivative
                .iter()
                .chain(&result.material_branch_mean_storage_derivative)
                .any(|value| !value.is_finite())
            || [
                result.mean_stored_energy_derivative,
                result.external_power_derivative,
                result.component_dissipation_derivative,
                result.connection_dissipation_derivative,
                result.material_dissipation_derivative,
                result.radiation_dissipation_derivative,
                result.total_dissipation_derivative(),
                result.power_residual_derivative,
            ]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "harmonic frequency state, power or storage derivative overflow",
            ));
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod integration_tests;
