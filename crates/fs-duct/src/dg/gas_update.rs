//! Uniform temporal material changes in the original acoustic coordinates.
//!
//! The fixed geometry stores compression chi=p/B and axial momentum rho*q/A.
//! Keeping those quantities through a zero-duration material intervention gives
//! p_new=(B_new/B_old)*p_old and q_new=(rho_old/rho_new)*q_old. The resulting
//! compressive/kinetic storage change is parameter work, never heat or damping.
use super::*;

mod viscothermal;

/// Measured storage change from one accepted uniform gas intervention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgGasMaterialUpdate {
    /// Complete acoustic and wall storage before the intervention [J].
    pub energy_before_j: f64,
    /// Complete acoustic and unchanged wall storage afterwards [J].
    pub energy_after_j: f64,
    /// Original acoustic storage, including any admitted diffusion memory [J].
    pub acoustic_energy_before_j: f64,
    /// Updated acoustic storage, including any admitted diffusion memory [J].
    pub acoustic_energy_after_j: f64,
    /// Signed external work INTO the represented acoustic system [J].
    /// The owner changing the background gas owes this transfer separately.
    pub parameter_work_j: f64,
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}
fn scaled(value: f64, factor: f64) -> Result<f64, DuctError> {
    let next = value * factor;
    if !next.is_finite() || (value != 0.0 && next == 0.0) {
        return Err(invalid(
            "uniform gas update lost a finite physical coefficient or coordinate",
        ));
    }
    Ok(next)
}
fn energy(values: &[f64]) -> f64 {
    0.5 * values.iter().map(|x| x * x).sum::<f64>()
}

impl DgDuct {
    /// Refresh a live spatially uniform numerical gas, at the SAME geometry.
    ///
    /// Preserve acoustic compression p/B and axial momentum rho*q/A; in energy
    /// coordinates these scale by sqrt(B_new/B_old) and sqrt(rho_old/rho_new).
    /// Wall displacement, velocity and every local Maxwell history are retained
    /// unchanged. Original geometry, node ordering, boundary laws and side-port
    /// declarations survive. Pressure effort maps and matched impedances follow
    /// the same new gas, with no extra forcing or instrument-specific retuning.
    ///
    /// Actual storage change enters cumulative_parameter_work_j ONCE, without
    /// resetting initial energy, prior dissipation or the clock. The finite
    /// nonnegative work budget is an admission, not supplied thermal energy.
    /// Rebuild the native propagator before continuing. All checks, including
    /// cancellation, precede publishing the operator and accepted state.
    ///
    /// This is a prescribed temporal material interface, not a heat, species
    /// transport or mean-flow solve. A change of background pressure does not
    /// supply a static wall traction, heating pulse or mass-injection source.
    /// Nonuniform gas and viscothermal gas memory refuse: neither is averaged,
    /// reset or assigned a guessed internal-variable transfer. Elastic/KV and
    /// regional Maxwell WALL memory remain supported and dissipate normally.
    /// Source-backed ducts require `rebind_uniform_gas_field` so their material
    /// catalog, query and phase policy cannot be silently erased.
    pub fn rebind_uniform_gas(
        &mut self,
        state: &mut DgState,
        gas: GasState,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        if self.gas_field.is_some() {
            return Err(invalid(
                "source-backed gas updates require a complete replacement gas field",
            ));
        }
        self.rebind_uniform_gas_data(
            state,
            gas,
            None,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    /// Source-retaining counterpart of `rebind_uniform_gas`. Both the original
    /// and replacement gas must be spatially uniform on the complete original
    /// length. Temperature, pressure, explicit composition and source models
    /// come from the shared field owner. Its complete catalog remains attached
    /// even when an explicitly zero species has no active local fraction.
    pub fn rebind_uniform_gas_field(
        &mut self,
        state: &mut DgState,
        field: &AxialGasField,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        poll(cancelled)?;
        if field.length_m().to_bits() != self.total_length.to_bits() {
            return Err(invalid(
                "live uniform gas field must cover the exact original duct length",
            ));
        }
        let gas = *field
            .uniform_state()
            .ok_or_else(|| {
                invalid(
                    "live gas field must be uniform; spatial history transfer is not represented",
                )
            })?
            .state();
        self.rebind_uniform_gas_data(
            state,
            gas,
            Some(field.clone()),
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    fn rebind_uniform_gas_data(
        &mut self,
        state: &mut DgState,
        gas: GasState,
        field: Option<AxialGasField>,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        poll(cancelled)?;
        crate::validate_gas(&gas)?;
        let m = self.field_dofs;
        let n = self.total_dofs();
        if self.viscothermal.is_some()
            || self
                .gas_field
                .as_ref()
                .is_some_and(|field| field.uniform_state().is_none())
        {
            return Err(invalid(
                "live gas update requires uniform gas without viscothermal gas histories",
            ));
        }
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != n
            || state.p.len() != m
            || state.q.len() != m
            || state.tilde.iter().any(|v| !v.is_finite())
            || !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
            || !state.cumulative_parameter_work_j.is_finite()
        {
            return Err(invalid(
                "live gas update needs its original finite state and a finite nonnegative work budget",
            ));
        }
        let energy_before_j = state.stored_energy_j();
        let acoustic_energy_before_j = energy(&state.tilde[..2 * m]);
        if !energy_before_j.is_finite() {
            return Err(invalid("original acoustic storage is not finite"));
        }
        let old_bulk = self.gas.density * self.gas.sound_speed * self.gas.sound_speed;
        let new_bulk = gas.density * gas.sound_speed * gas.sound_speed;
        let pressure_scale = det::sqrt(new_bulk) / det::sqrt(old_bulk);
        let momentum_scale = det::sqrt(self.gas.density) / det::sqrt(gas.density);
        let speed_scale = gas.sound_speed / self.gas.sound_speed;
        if [
            old_bulk,
            new_bulk,
            pressure_scale,
            momentum_scale,
            speed_scale,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "updated gas bulk modulus or coordinate scales are unrepresentable",
            ));
        }
        let mut candidate = self.clone();
        // Uniform rho and B commute with each ORIGINAL geometric mass matrix.
        // These are the exact material factors of the native weak operator,
        // including matched boundaries; no new spatial basis is introduced.
        for element in &mut candidate.elements {
            poll(cancelled)?;
            for (values, factor) in [
                (&mut element.mp_sqrt, 1.0 / pressure_scale),
                (&mut element.mp_inv_sqrt, pressure_scale),
                (&mut element.mq_sqrt, 1.0 / momentum_scale),
                (&mut element.mq_inv_sqrt, momentum_scale),
            ] {
                for value in values {
                    *value = scaled(*value, factor)?;
                }
            }
        }
        for i in 0..2 * m {
            poll(cancelled)?;
            for j in 0..2 * m {
                candidate.l_tilde[i * n + j] = scaled(self.l_tilde[i * n + j], speed_scale)?;
            }
        }
        // The reciprocal lateral wall connection touches gas pressure only.
        // Every internal wall spring, viscosity and history block is unchanged.
        for i in 0..m {
            poll(cancelled)?;
            for j in 2 * m..n {
                candidate.l_tilde[i * n + j] = scaled(self.l_tilde[i * n + j], pressure_scale)?;
                candidate.l_tilde[j * n + i] = scaled(self.l_tilde[j * n + i], pressure_scale)?;
            }
        }
        if let Some(g) = &mut candidate.g_block {
            for value in g {
                *value = scaled(*value, speed_scale)?;
            }
        }
        for (_, column) in &mut candidate.ports {
            poll(cancelled)?;
            if column[m..].iter().any(|value| *value != 0.0) {
                return Err(invalid(
                    "live gas update requires the original pressure/volume-flow port convention",
                ));
            }
            for value in &mut column[..m] {
                *value = scaled(*value, pressure_scale)?;
            }
        }
        if let Some(wall) = &mut candidate.compliant_wall {
            wall.rescale_uniform_gas_observation(pressure_scale)?;
        }
        candidate.gas = gas;
        candidate.outlet_gas = gas;
        candidate.gas_field = field;
        candidate.operator_identity = hash_operator(
            candidate.mesh,
            candidate.inlet,
            candidate.outlet,
            &candidate.elements,
            &candidate.l_tilde,
            &candidate.ports,
        );
        let mut next = state.clone();
        let mut parameter_work_j = 0.0;
        for (start, factor) in [(0, pressure_scale), (m, momentum_scale)] {
            poll(cancelled)?;
            for i in start..start + m {
                let before = state.tilde[i];
                let after = scaled(before, factor)?;
                next.tilde[i] = after;
                // Avoid subtracting a potentially huge unchanged wall storage.
                parameter_work_j += 0.5 * (after - before) * (after + before);
            }
        }
        next.operator_identity = candidate.operator_identity;
        (next.p, next.q) = candidate.from_energy_coords(&next.tilde);
        next.nodal_dirty = false;
        next.modal_anchor_p.clear();
        next.modal_anchor_q.clear();
        next.modal_anchor_dt_bits = None;
        next.unforced_steps = 0;
        let energy_after_j = next.stored_energy_j();
        let acoustic_energy_after_j = energy(&next.tilde[..2 * m]);
        let accumulated = state.cumulative_parameter_work_j + parameter_work_j;
        if !energy_after_j.is_finite()
            || next.p.iter().chain(&next.q).any(|v| !v.is_finite())
            || !parameter_work_j.is_finite()
            || parameter_work_j.abs() > maximum_absolute_parameter_work_j
            || !accumulated.is_finite()
            || (parameter_work_j != 0.0 && accumulated == state.cumulative_parameter_work_j)
            || (state.cumulative_parameter_work_j != 0.0 && accumulated == parameter_work_j)
        {
            return Err(invalid(
                "live gas parameter work exceeded its budget or finite physical resolution",
            ));
        }
        if candidate.compliant_wall.is_some() {
            candidate.compliant_wall_state(&next)?;
        }
        next.cumulative_parameter_work_j = accumulated;
        poll(cancelled)?;
        let _old_duct = core::mem::replace(self, candidate);
        let _old_state = core::mem::replace(state, next);
        Ok(DgGasMaterialUpdate {
            energy_before_j,
            energy_after_j,
            acoustic_energy_before_j,
            acoustic_energy_after_j,
            parameter_work_j,
        })
    }
}

#[cfg(test)]
mod tests;
