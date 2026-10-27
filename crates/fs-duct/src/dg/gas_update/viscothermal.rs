//! Live transfer of the native uniform circular Bessel diffusion histories.
//!
//! The original dimensionless Jacobi sections, radius and pole count identify
//! each arm. This is not a match between arbitrary state-space realizations.
//! Series arm inertance is (rho/A)*(gain/rate); shunt arm compliance is
//! (A/B)*(gain/rate). Preserve each arm's momentum or compression volume at a
//! zero-duration material intervention, and account its actual storage change.
use super::*;
use super::super::viscothermal::ViscothermalMemory;

impl DgDuct {
    /// Change a live uniform gas without resetting its circular wall diffusion.
    ///
    /// Radius, geometry, mesh, exact angular band and dimensionless Bessel
    /// sections are retained. Each series arm preserves inertive momentum and
    /// each shunt arm preserves compression volume. The main fields retain p/B
    /// and rho*q/A; normalized wall coordinates and Maxwell histories survive.
    /// Actual wave AND diffusion-memory storage change enters parameter work
    /// once. Viscosity/conductivity-only changes therefore change future losses
    /// without instantaneous work when density, bulk modulus and gamma agree.
    ///
    /// The new gas must independently pass the original band's complex-response
    /// AND real-loss admission under maximum_relative_error. Changing pole count,
    /// radius, thermal-arm support, nonuniform gas or spatial memory refuses.
    /// This fixed-section physical transfer does not apply to generic fitted
    /// realizations. Rebuild the original propagator before continuing.
    ///
    /// This prescribed background intervention supplies no heat/species solve,
    /// mass injection or mean-pressure traction. Its external owner owes the
    /// signed parameter work separately. All refusals leave both owners intact.
    pub fn rebind_uniform_viscothermal_gas(
        &mut self,
        state: &mut DgState,
        gas: GasState,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        if self.gas_field.is_some() {
            return Err(invalid(
                "source-backed viscothermal updates require the complete gas field",
            ));
        }
        self.rebind_uniform_viscothermal_data(
            state,
            gas,
            None,
            maximum_relative_error,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    /// Source-retaining counterpart with the complete original-length field.
    /// The material catalog, query, explicit composition and phase policy remain
    /// attached, including explicitly zero species. Both fields must be uniform.
    pub fn rebind_uniform_viscothermal_gas_field(
        &mut self,
        state: &mut DgState,
        field: &AxialGasField,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        poll(cancelled)?;
        if field.length_m().to_bits() != self.total_length.to_bits() {
            return Err(invalid(
                "live viscothermal gas field must cover the exact original length",
            ));
        }
        let gas = *field.uniform_state().ok_or_else(|| {
            invalid("live viscothermal transfer requires a uniform replacement gas field")
        })?.state();
        self.rebind_uniform_viscothermal_data(
            state,
            gas,
            Some(field.clone()),
            maximum_relative_error,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn rebind_uniform_viscothermal_data(
        &mut self,
        state: &mut DgState,
        gas: GasState,
        field: Option<AxialGasField>,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, DuctError> {
        poll(cancelled)?;
        crate::validate_gas(&gas)?;
        let Some(ViscothermalMemory::Uniform(previous)) = &self.viscothermal else {
            return Err(invalid(
                "live diffusion transfer requires the native uniform circular Bessel owner",
            ));
        };
        if self.gas_field.as_ref().is_some_and(|f| f.uniform_state().is_none()) {
            return Err(invalid("live diffusion transfer cannot average a spatial gas field"));
        }
        let radius = self.uniform_radius_m.ok_or_else(|| {
            invalid("live diffusion transfer requires the original circular radius")
        })?;
        let m = self.field_dofs;
        let n = self.total_dofs();
        let series_count = previous.series_terms().len();
        let shunt_count = previous.shunt_terms().len();
        let prefix = m * (2 + series_count + shunt_count);
        if prefix > n
            || state.operator_identity != self.operator_identity
            || state.tilde.len() != n
            || state.p.len() != m
            || state.q.len() != m
            || state.tilde.iter().any(|value| !value.is_finite())
            || !state.cumulative_parameter_work_j.is_finite()
            || !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
        {
            return Err(invalid(
                "live diffusion transfer needs its original finite state and finite work budget",
            ));
        }
        if (self.gas.gamma == 1.0) != (gas.gamma == 1.0) {
            return Err(invalid(
                "live diffusion transfer cannot create or discard thermal arm support",
            ));
        }
        let energy_before_j = state.stored_energy_j();
        let acoustic_energy_before_j = energy(&state.tilde[..prefix]);
        let old_bulk = self.gas.density * self.gas.sound_speed * self.gas.sound_speed;
        let new_bulk = gas.density * gas.sound_speed * gas.sound_speed;
        let pressure_scale = det::sqrt(new_bulk) / det::sqrt(old_bulk);
        let momentum_scale = det::sqrt(self.gas.density) / det::sqrt(gas.density);
        let thermal_scale = if shunt_count == 0 {
            pressure_scale
        } else {
            pressure_scale * (det::sqrt(self.gas.gamma - 1.0) / det::sqrt(gas.gamma - 1.0))
        };
        if !energy_before_j.is_finite()
            || [old_bulk, new_bulk, pressure_scale, momentum_scale, thermal_scale]
                .iter().any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid("live diffusion transfer has unrepresentable storage scales"));
        }
        // Rebuild the ORIGINAL dimensionless sections at exactly the same count.
        // Their ordering follows the unchanged Jacobi problem, not a match of
        // gas-dependent relaxation rates. The old angular band is not rounded
        // through a Hz conversion. Admission uses the NEW gas.
        let model = fs_phs::ZkDiffusive::new(
            radius,
            gas.density,
            gas.dynamic_viscosity,
            gas.gamma,
            gas.prandtl,
            series_count,
            previous.band_rad_s(),
            maximum_relative_error,
        ).map_err(|_| invalid(
            "replacement gas does not resolve the original Bessel diffusion band",
        ))?;
        if model.series_terms().len() != series_count || model.shunt_terms().len() != shunt_count {
            return Err(invalid("replacement diffusion changed native arm support"));
        }
        poll(cancelled)?;
        // The native cold assembly owns the new mass matrices, reciprocal wave
        // operator, matched impedances, physical ports and complete diffusion
        // block. No second spatial or stepping solver is introduced.
        let mut candidate = DgDuct::new(
            &self.geometry, &gas, self.mesh, self.inlet, self.outlet,
        )?.with_side_ports(&self.side_ports)?.with_uniform_viscothermal_model(model)?;
        if candidate.total_dofs() != prefix || candidate.elements.len() != self.elements.len() {
            return Err(invalid("replacement diffusion changed the original spatial basis"));
        }
        for (old, new) in self.elements.iter().zip(&candidate.elements) {
            poll(cancelled)?;
            if old.x_left.to_bits() != new.x_left.to_bits()
                || old.length.to_bits() != new.length.to_bits()
                || old.nodes_x != new.nodes_x
                || old.areas != new.areas
            {
                return Err(invalid("replacement diffusion changed the original geometry"));
            }
        }
        candidate.gas_field = field;
        if let Some(old_wall) = &self.compliant_wall {
            // Retain every structural coordinate, internal spring/dashpot block
            // and regional Maxwell arm. Only reciprocal pressure loading changes.
            let mut operator = self.l_tilde.clone();
            for i in 0..prefix {
                poll(cancelled)?;
                operator[i * n..i * n + prefix]
                    .copy_from_slice(&candidate.l_tilde[i * prefix..(i + 1) * prefix]);
            }
            for i in 0..m {
                for j in prefix..n {
                    operator[i * n + j] = scaled(self.l_tilde[i * n + j], pressure_scale)?;
                    operator[j * n + i] = scaled(self.l_tilde[j * n + i], pressure_scale)?;
                }
            }
            let mut wall = old_wall.clone();
            wall.rescale_uniform_gas_observation(pressure_scale)?;
            candidate.compliant_wall = Some(wall);
            candidate.memory_dofs = self.memory_dofs;
            candidate.l_tilde = operator;
            for (_, column) in &mut candidate.ports {
                column.resize(n, 0.0);
            }
        } else if prefix != n {
            return Err(invalid("unidentified extra coordinates in live diffusion state"));
        }
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
        for (start, end, factor) in [
            (0, m, pressure_scale),
            (m, (2 + series_count) * m, momentum_scale),
            ((2 + series_count) * m, prefix, thermal_scale),
        ] {
            poll(cancelled)?;
            for i in start..end {
                let before = state.tilde[i];
                let after = scaled(before, factor)?;
                next.tilde[i] = after;
                // Include all actual diffusion storage, without subtracting a
                // potentially much larger unchanged structural wall energy.
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
        let acoustic_energy_after_j = energy(&next.tilde[..prefix]);
        let accumulated = state.cumulative_parameter_work_j + parameter_work_j;
        if !energy_after_j.is_finite()
            || next.p.iter().chain(&next.q).any(|value| !value.is_finite())
            || !parameter_work_j.is_finite()
            || parameter_work_j.abs() > maximum_absolute_parameter_work_j
            || !accumulated.is_finite()
            || (parameter_work_j != 0.0 && accumulated == state.cumulative_parameter_work_j)
            || (state.cumulative_parameter_work_j != 0.0 && accumulated == parameter_work_j)
        {
            return Err(invalid(
                "live diffusion parameter work exceeds its budget or finite resolution",
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
