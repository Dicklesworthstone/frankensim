//! Spatial counterpart of the uniform entropy-preserving material interface.
//!
//! E=M_inf^(1/2) and F=M_inf^(-1/2) are the ORIGINAL element mass maps.
//! The weak compression is E (x_p - sum a_n P_n^T y_n). Thus
//! x_p,new = F_new E_old (x_p,old - sum a_old P_old^T y_old)
//!           + sum a_new P_new^T y_new.
//! Matrix order matters: spatial heat capacity and density do not commute
//! with the polynomial mass. Each y and z retains its local entropy density;
//! the original positive wall rule fixes each physical volume and radial mode.

use super::*;
use crate::dg::gas_update::conserved_map;

fn paired_materials(
    old: &AxialGasField,
    new: &AxialGasField,
    law: &RelaxationAbsorption,
    point: f64,
) -> Result<(RelaxationState, RelaxationState), GasFieldDuctError> {
    let a = RelaxingMedium::from_thermochemical(&old.state_at(point)?)?;
    let b = RelaxingMedium::from_thermochemical(&new.state_at(point)?)?;
    if a.components() != b.components() {
        return Err(
            invalid("spatial molecular refresh cannot transport or relabel species").into(),
        );
    }
    let a = law.evaluate(&a)?;
    let b = law.evaluate(&b)?;
    if a.processes().len() != b.processes().len()
        || a.processes().iter().zip(b.processes()).any(|(a, b)| {
            a.species != b.species
                || a.index != b.index
                || (a.heat_capacity_j_kg_k > 0.0) != (b.heat_capacity_j_kg_k > 0.0)
        })
    {
        return Err(
            invalid("spatial molecular refresh cannot create or discard reservoirs").into(),
        );
    }
    Ok((a, b))
}

impl DgDuct {
    /// Continue a nonuniform molecular gas on its original spatial wall rule.
    ///
    /// Primitive temperature and common pressure may change. Original station
    /// positions and every station's species fractions remain fixed, including
    /// a spatial composition gradient. No inlet average, pole matching, history
    /// reset, or pointwise mass lumping is used. The same material declaration
    /// reevaluates heat capacities, relaxation rates and Bessel wall diffusion.
    /// Every point must still pass the original angular band's complex AND
    /// real-loss admission. Uniform/nonuniform transitions are supported while
    /// retaining the same spatial quadrature; exact-mass circular memory uses
    /// `rebind_uniform_molecular_gas_field` instead.
    ///
    /// Preserve weak total compression (including thermal histories), weak
    /// axial momentum, and every local thermal/molecular entropy perturbation
    /// and viscous-arm momentum. Structural coordinates and histories remain
    /// unchanged and pressure loading follows the new frozen pressure mass.
    /// The exact quadratic storage change is one signed external parameter
    /// work entry; clock, initial energy and previous loss/work survive.
    ///
    /// This is an instantaneous prescribed background intervention. It does
    /// not solve heat/species transport, mean flow, condensation, background
    /// inventory or mean-pressure wall traction. Remeshing and reservoir
    /// support changes refuse. A changed operator invalidates its propagator.
    /// All refusals, including cancellation and work budget, leave both owners
    /// unchanged. Between interventions the original linear frozen-field
    /// approximation still applies.
    pub fn rebind_spatial_molecular_gas_field(
        &mut self,
        state: &mut DgState,
        field: &AxialGasField,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, GasFieldDuctError> {
        poll(cancelled)?;
        let m = self.field_dofs;
        let n = self.total_dofs();
        if field.length_m().to_bits() != self.total_length.to_bits()
            || !maximum_relative_error.is_finite()
            || !(0.0 < maximum_relative_error && maximum_relative_error <= 0.5)
            || !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
            || state.operator_identity != self.operator_identity
            || state.tilde.len() != n
            || state.p.len() != m
            || state.q.len() != m
            || state.tilde.iter().any(|v| !v.is_finite())
            || !state.cumulative_parameter_work_j.is_finite()
        {
            return Err(invalid("spatial molecular refresh needs its original finite state, length and accuracy/work budgets").into());
        }
        let original_field = self.gas_field.as_ref().ok_or_else(|| {
            invalid("spatial molecular refresh requires its original complete source field")
        })?;
        if original_field.stations_m() != field.stations_m()
            || original_field
                .primitive_states()
                .iter()
                .zip(field.primitive_states())
                .any(|(a, b)| a.composition != b.composition)
        {
            return Err(invalid("spatial molecular refresh preserves station positions and complete species fractions").into());
        }
        let previous = self.molecular.as_ref().ok_or_else(|| {
            invalid("spatial molecular refresh requires original molecular histories")
        })?;
        let law = previous
            .field_law
            .as_ref()
            .ok_or_else(|| {
                invalid("spatial molecular refresh requires its original constitutive declaration")
            })?
            .clone();
        let Some(ViscothermalMemory::Spatial(old_wall)) = &self.viscothermal else {
            return Err(
                invalid("spatial molecular refresh requires native spatial diffusion").into(),
            );
        };
        for &x in field.stations_m() {
            poll(cancelled)?;
            paired_materials(original_field, field, &law, x)?;
        }
        let resolution = DgWallResolution {
            maximum_relative_error,
            ..self
                .spatial_viscothermal_admission()
                .expect("spatial wall")
                .resolution
        };
        let mut candidate =
            Self::new_with_gas_field(&self.geometry, field, self.mesh, self.inlet, self.outlet)?
                .with_side_ports(&self.side_ports)?
                .with_spatial_viscothermal(resolution)?
                .with_molecular_relaxation_field_controlled(law, cancelled)?;
        let renewed = candidate
            .molecular
            .as_ref()
            .expect("native molecular assembly");
        let prefix = candidate.total_dofs();
        if prefix > n
            || renewed.modes.len() != previous.modes.len()
            || self.elements.len() != candidate.elements.len()
        {
            return Err(
                invalid("spatial molecular refresh changed original coordinate support").into(),
            );
        }
        // An empty law has no extra thermodynamic coordinates. Preserve the
        // original nonmolecular transfer and arithmetic, plus the source law.
        if previous.modes.is_empty() {
            let molecular = candidate.molecular;
            let mut owner = self.clone();
            owner.molecular = None;
            let mut next = state.clone();
            let report = owner.rebind_spatial_viscothermal_gas_field(
                &mut next,
                field,
                maximum_relative_error,
                maximum_absolute_parameter_work_j,
                cancelled,
            )?;
            owner.molecular = molecular;
            poll(cancelled)?;
            *self = owner;
            *state = next;
            return Ok(report);
        }

        let local = self.mesh.order + 1;
        let mut next = state.clone();
        let mut pressure_maps = Vec::with_capacity(self.elements.len());
        for (e, (old, new)) in self.elements.iter().zip(&candidate.elements).enumerate() {
            poll(cancelled)?;
            if old.x_left.to_bits() != new.x_left.to_bits()
                || old.length.to_bits() != new.length.to_bits()
                || old.nodes_x != new.nodes_x
                || old.areas != new.areas
            {
                return Err(
                    invalid("spatial molecular refresh cannot remesh a continuing state").into(),
                );
            }
            let pressure = conserved_map(&old.mp_sqrt, &new.mp_sqrt, &new.mp_inv_sqrt, local)?;
            let flow = conserved_map(&old.mq_sqrt, &new.mq_sqrt, &new.mq_inv_sqrt, local)?;
            for (offset, map) in [(0, &pressure), (m, &flow)] {
                let start = offset + e * local;
                for i in 0..local {
                    next.tilde[start + i] = (0..local)
                        .map(|j| map[i * local + j] * state.tilde[start + j])
                        .sum();
                }
            }
            pressure_maps.push(pressure);
        }
        let Some(ViscothermalMemory::Spatial(new_wall)) = &candidate.viscothermal else {
            unreachable!("native spatial molecular assembly")
        };
        old_wall.transfer_molecular_viscous_momenta(
            new_wall,
            m,
            &state.tilde,
            &mut next.tilde,
            cancelled,
        )?;
        let old_points = old_wall.thermal_quadrature(m);
        let new_points = new_wall.thermal_quadrature(m);
        if old_points.len() != new_points.len() {
            return Err(
                invalid("spatial molecular refresh changed the original wall quadrature").into(),
            );
        }
        let mut at = 0;
        for (old_point, new_point) in old_points.iter().zip(&new_points) {
            poll(cancelled)?;
            if old_point.point_m.to_bits() != new_point.point_m.to_bits()
                || old_point.field_start != new_point.field_start
                || old_point.terms.len() != new_point.terms.len()
                || old_point.projection.len() != local
                || new_point.projection.len() != local
                || old_point.field_start % local != 0
                || old_point.field_start + local > m
            {
                return Err(invalid("spatial molecular refresh changed a thermal section").into());
            }
            let (a, b) = paired_materials(
                original_field,
                field,
                previous.field_law.as_ref().expect("source law"),
                old_point.point_m,
            )?;
            let thermal = entropy_scale(
                a.gas(),
                a.frozen_cv_j_kg_k() + a.gas().specific_gas_constant,
                b.gas(),
                b.frozen_cv_j_kg_k() + b.gas().specific_gas_constant,
            )?;
            let reservoir: Vec<_> = a
                .processes()
                .iter()
                .zip(b.processes())
                .filter(|(a, _)| a.heat_capacity_j_kg_k > 0.0)
                .map(|(a_mode, b_mode)| {
                    entropy_scale(
                        a.gas(),
                        a_mode.heat_capacity_j_kg_k,
                        b.gas(),
                        b_mode.heat_capacity_j_kg_k,
                    )
                })
                .collect::<Result<_, _>>()?;
            let first = old_point.field_start;
            let map = &pressure_maps[first / local];
            for _ in &old_point.terms {
                let old = previous
                    .modes
                    .get(at)
                    .ok_or_else(|| invalid("original molecular mode order is incomplete"))?;
                let new = renewed
                    .modes
                    .get(at)
                    .ok_or_else(|| invalid("replacement molecular mode order is incomplete"))?;
                at += 1;
                if old.start != new.start
                    || old.rates.len() != reservoir.len()
                    || new.rates.len() != reservoir.len()
                    || old.projection.width(m) != 1
                    || new.projection.width(m) != 1
                    || !same_projection_support(&old.projection, &new.projection)
                {
                    return Err(invalid(
                        "spatial molecular refresh changed native reservoir support",
                    )
                    .into());
                }
                next.tilde[new.start] = scaled(state.tilde[old.start], thermal)?;
                for (j, &factor) in reservoir.iter().enumerate() {
                    next.tilde[new.start + j + 1] = scaled(state.tilde[old.start + j + 1], factor)?;
                }
                let old_entries = old.projection.entries(0, old.pressure_scale);
                let new_entries = new.projection.entries(0, new.pressure_scale);
                if old_entries.len() != local
                    || new_entries.len() != local
                    || old_entries
                        .iter()
                        .zip(&new_entries)
                        .enumerate()
                        .any(|(i, ((a, _), (b, _)))| *a != first + i || *b != first + i)
                {
                    return Err(invalid(
                        "spatial molecular refresh lost the native pressure functional",
                    )
                    .into());
                }
                for i in 0..local {
                    let mapped: f64 = (0..local)
                        .map(|j| map[i * local + j] * old_entries[j].1)
                        .sum();
                    // Difference form preserves an exact identity intervention.
                    next.tilde[first + i] +=
                        new_entries[i].1 * next.tilde[new.start] - mapped * state.tilde[old.start];
                }
            }
        }
        if at != previous.modes.len()
            || previous
                .modes
                .last()
                .is_none_or(|mode| mode.start + mode.rates.len() + 1 != prefix)
            || next.tilde[..prefix].iter().any(|v| !v.is_finite())
        {
            return Err(
                invalid("spatial molecular physical transfer is incomplete or nonfinite").into(),
            );
        }

        if let Some(old_wall) = &self.compliant_wall {
            let mut operator = self.l_tilde.clone();
            for i in 0..prefix {
                poll(cancelled)?;
                operator[i * n..i * n + prefix]
                    .copy_from_slice(&candidate.l_tilde[i * prefix..(i + 1) * prefix]);
            }
            for (e, map) in pressure_maps.iter().enumerate() {
                poll(cancelled)?;
                let first = e * local;
                for i in 0..local {
                    for j in prefix..n {
                        let value: f64 = (0..local)
                            .map(|k| map[i * local + k] * self.l_tilde[(first + k) * n + j])
                            .sum();
                        if !value.is_finite() {
                            return Err(invalid(
                                "spatial molecular pressure/wall coupling is unrepresentable",
                            )
                            .into());
                        }
                        operator[(first + i) * n + j] = value;
                        operator[j * n + first + i] = -value;
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
            return Err(invalid("spatial molecular refresh found unowned coordinates").into());
        }
        candidate.operator_identity = hash_operator(
            candidate.mesh,
            candidate.inlet,
            candidate.outlet,
            &candidate.elements,
            &candidate.l_tilde,
            &candidate.ports,
        );
        let energy_before_j = state.stored_energy_j();
        let acoustic_energy_before_j = energy(&state.tilde[..prefix]);
        let parameter_work_j: f64 = state.tilde[..prefix]
            .iter()
            .zip(&next.tilde[..prefix])
            .map(|(a, b)| 0.5 * (b - a) * (b + a))
            .sum();
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
        if !energy_before_j.is_finite()
            || !energy_after_j.is_finite()
            || next.p.iter().chain(&next.q).any(|v| !v.is_finite())
            || !parameter_work_j.is_finite()
            || parameter_work_j.abs() > maximum_absolute_parameter_work_j
            || !accumulated.is_finite()
            || (parameter_work_j != 0.0 && accumulated == state.cumulative_parameter_work_j)
            || (state.cumulative_parameter_work_j != 0.0 && accumulated == parameter_work_j)
        {
            return Err(invalid(
                "spatial molecular parameter work exceeds its allowance or finite resolution",
            )
            .into());
        }
        if candidate.compliant_wall.is_some() {
            candidate.compliant_wall_state(&next)?;
        }
        next.cumulative_parameter_work_j = accumulated;
        poll(cancelled)?;
        *self = candidate;
        *state = next;
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
