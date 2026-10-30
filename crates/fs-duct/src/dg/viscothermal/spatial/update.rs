//! Live material transfer on the original spatial quadrature and DG basis.
//!
//! The conserved weak acoustic variables are Mp*p and Mq*q, not pointwise
//! values sampled at the interpolation nodes. With x=M^(1/2)*u their transfer
//! is x_new=M_new^(-1/2)*M_old^(1/2)*x_old. These matrices need not commute.
//! Each retained quadrature arm instead has a scalar physical inertia or
//! compliance, and preserves its original momentum or compression volume.
use super::*;
use crate::dg::gas_update::{DgGasMaterialUpdate, conserved_map, energy, poll, scaled};

impl SpatialViscothermal {
    /// The molecular owner replaces thermal channels but retains these exact
    /// native viscosity coordinates. Match physical sections, never pole rates.
    pub(in crate::dg) fn transfer_molecular_viscous_momenta(
        &self,
        replacement: &Self,
        field_dofs: usize,
        before: &[f64],
        after: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        let old: Vec<_> = self
            .channels
            .iter()
            .filter(|channel| channel.field_start >= field_dofs)
            .collect();
        let new: Vec<_> = replacement
            .channels
            .iter()
            .filter(|channel| channel.field_start >= field_dofs)
            .collect();
        if old.len() != new.len() {
            return Err(invalid("molecular refresh changed viscous channel support"));
        }
        for (old, new) in old.into_iter().zip(new) {
            poll(cancelled)?;
            if old.point_m.to_bits() != new.point_m.to_bits()
                || old.radius_m.to_bits() != new.radius_m.to_bits()
                || old.field_start != new.field_start
                || old.memory_start != new.memory_start
                || old.terms.len() != new.terms.len()
                || old.memory_start + old.terms.len() > before.len()
                || new.memory_start + new.terms.len() > after.len()
            {
                return Err(invalid(
                    "molecular refresh changed a physical viscous section",
                ));
            }
            let factor = old.storage_root / new.storage_root;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid(
                    "molecular viscous momentum scale is unrepresentable",
                ));
            }
            for arm in 0..old.terms.len() {
                let index = old.memory_start + arm;
                after[index] = scaled(before[index], factor)?;
            }
        }
        Ok(())
    }
}

impl DgDuct {
    /// Refresh the actual regional gas field without resetting spatial diffusion.
    ///
    /// The original geometry, element boundaries, polynomial basis, quadrature
    /// points, radius at each point, angular band and native Bessel pole support
    /// must be identical. New source data are evaluated by the original gas-mass
    /// assembly and local wall producer; no inlet average or new damping fit is
    /// substituted. Every new point must satisfy the explicit complex-response
    /// AND real-loss error ceiling on the retained band.
    ///
    /// Preserve each original weak compression Mp*p and momentum Mq*q. This is
    /// a finite-basis conservation statement, not exact pointwise p/B or rho*q
    /// conservation for arbitrary spatial variation within an element. Each
    /// fixed Gauss-point viscous arm preserves its weighted physical momentum;
    /// each thermal arm preserves compression volume. Their dimensionless
    /// Jacobi sections remain fixed, not matched by gas-dependent pole values.
    ///
    /// Native structural displacement/velocity and all regional Maxwell wall
    /// histories are retained. Reciprocal pressure loading uses the NEW mass
    /// root, while structural observations retain their physical coordinate
    /// map. Actual main and diffusion storage changes enter parameter work once,
    /// without advancing time or resetting previous loss, work or initial energy.
    ///
    /// The external material owner owes this signed work. No background heat,
    /// species transport, mean-pressure traction or mass injection is inferred.
    /// Changing mesh/quadrature/pole topology or creating/dropping thermal arms
    /// refuses. This API requires the spatial memory owner; uniform exact-mass
    /// memories use their existing explicit transfer instead. Callers rebuild
    /// the original propagator after acceptance. Cancellation is checked around
    /// the existing bounded cold assemblies and through transfer work; every
    /// refusal leaves the accepted duct and state unchanged.
    pub fn rebind_spatial_viscothermal_gas_field(
        &mut self,
        state: &mut DgState,
        field: &AxialGasField,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, GasFieldDuctError> {
        poll(cancelled)?;
        if self.molecular.is_some() {
            return Err(
                invalid("live gas updates require a molecular history-transfer law").into(),
            );
        }
        let Some(ViscothermalMemory::Spatial(old_memory)) = &self.viscothermal else {
            return Err(
                invalid("spatial gas refresh requires its native spatial diffusion owner").into(),
            );
        };
        let m = self.field_dofs;
        let n = self.total_dofs();
        let prefix = 2 * m + old_memory.admission.memory_dofs;
        if field.length_m().to_bits() != self.total_length.to_bits()
            || prefix > n
            || state.operator_identity != self.operator_identity
            || state.tilde.len() != n
            || state.p.len() != m
            || state.q.len() != m
            || state.tilde.iter().any(|value| !value.is_finite())
            || !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
            || !state.cumulative_parameter_work_j.is_finite()
        {
            return Err(invalid(
                "spatial gas refresh needs its original finite state, length and work budget",
            )
            .into());
        }
        let energy_before_j = state.stored_energy_j();
        let acoustic_energy_before_j = energy(&state.tilde[..prefix]);
        if !energy_before_j.is_finite() {
            return Err(invalid("original spatial acoustic storage is not finite").into());
        }
        let resolution = DgWallResolution {
            maximum_relative_error,
            ..old_memory.admission.resolution
        };
        // Reuse both original physical owners, including field/geometry span
        // subdivision, exact mass quadrature, matched endpoints and side ports.
        let mut candidate =
            Self::new_with_gas_field(&self.geometry, field, self.mesh, self.inlet, self.outlet)?
                .with_side_ports(&self.side_ports)?
                .with_spatial_viscothermal(resolution)?;
        poll(cancelled)?;
        if candidate.total_dofs() != prefix || candidate.elements.len() != self.elements.len() {
            return Err(invalid(
                "spatial gas refresh cannot change the original basis or memory support",
            )
            .into());
        }
        let local = self.mesh.order + 1;
        let mut next = state.clone();
        let mut pressure_maps = Vec::with_capacity(self.elements.len());
        for (element, (old, new)) in self.elements.iter().zip(&candidate.elements).enumerate() {
            poll(cancelled)?;
            if old.x_left.to_bits() != new.x_left.to_bits()
                || old.length.to_bits() != new.length.to_bits()
                || old.nodes_x != new.nodes_x
                || old.areas != new.areas
            {
                return Err(invalid("spatial gas refresh cannot remesh a continuing state").into());
            }
            let pressure = conserved_map(&old.mp_sqrt, &new.mp_sqrt, &new.mp_inv_sqrt, local)?;
            let flow = conserved_map(&old.mq_sqrt, &new.mq_sqrt, &new.mq_inv_sqrt, local)?;
            for (offset, map) in [(0, &pressure), (m, &flow)] {
                let start = offset + element * local;
                for i in 0..local {
                    next.tilde[start + i] = (0..local)
                        .map(|j| map[i * local + j] * state.tilde[start + j])
                        .sum();
                }
            }
            pressure_maps.push(pressure);
        }
        let Some(ViscothermalMemory::Spatial(new_memory)) = &candidate.viscothermal else {
            return Err(invalid("replacement lost its native spatial diffusion owner").into());
        };
        if old_memory.channels.len() != new_memory.channels.len() {
            return Err(invalid("replacement changed spatial diffusion channel support").into());
        }
        for (old, new) in old_memory.channels.iter().zip(&new_memory.channels) {
            poll(cancelled)?;
            if old.point_m.to_bits() != new.point_m.to_bits()
                || old.radius_m.to_bits() != new.radius_m.to_bits()
                || old.field_start != new.field_start
                || old.memory_start != new.memory_start
                || old.terms.len() != new.terms.len()
            {
                return Err(
                    invalid("replacement changed physical diffusion section identity").into(),
                );
            }
            // Fixed Jacobi weights cancel. This preserves scalar arm momentum
            // or compression and gives exactly one for transport-only changes.
            let factor = old.storage_root / new.storage_root;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(
                    invalid("replacement diffusion storage scale is unrepresentable").into(),
                );
            }
            for arm in 0..old.terms.len() {
                let index = old.memory_start + arm;
                next.tilde[index] = scaled(state.tilde[index], factor)?;
            }
        }
        if let Some(old_wall) = &self.compliant_wall {
            let mut operator = self.l_tilde.clone();
            for i in 0..prefix {
                poll(cancelled)?;
                operator[i * n..i * n + prefix]
                    .copy_from_slice(&candidate.l_tilde[i * prefix..(i + 1) * prefix]);
            }
            // C_new = Mp_new^(-1/2) Mp_old^(1/2) C_old. Pressure
            // receives -C_new; the structural velocity receives +C_new^T.
            for (element, map) in pressure_maps.iter().enumerate() {
                poll(cancelled)?;
                let base = element * local;
                for i in 0..local {
                    for j in prefix..n {
                        let value: f64 = (0..local)
                            .map(|k| map[i * local + k] * self.l_tilde[(base + k) * n + j])
                            .sum();
                        if !value.is_finite() {
                            return Err(invalid(
                                "replacement pressure/wall coupling is unrepresentable",
                            )
                            .into());
                        }
                        operator[(base + i) * n + j] = value;
                        operator[j * n + base + i] = -value;
                    }
                }
            }
            let mut wall = old_wall.clone();
            wall.retain_physical_observation_maps(&self.elements);
            candidate.compliant_wall = Some(wall);
            candidate.memory_dofs = self.memory_dofs;
            candidate.l_tilde = operator;
            for (_, port) in &mut candidate.ports {
                port.resize(n, 0.0);
            }
        } else if prefix != n {
            return Err(invalid("unidentified coordinates after spatial gas memory").into());
        }
        candidate.operator_identity = hash_operator(
            candidate.mesh,
            candidate.inlet,
            candidate.outlet,
            &candidate.elements,
            &candidate.l_tilde,
            &candidate.ports,
        );
        let mut parameter_work_j = 0.0;
        for i in 0..prefix {
            if i % local == 0 {
                poll(cancelled)?;
            }
            let before = state.tilde[i];
            let after = next.tilde[i];
            parameter_work_j += 0.5 * (after - before) * (after + before);
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
                "spatial gas parameter work exceeds its budget or finite resolution",
            )
            .into());
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
