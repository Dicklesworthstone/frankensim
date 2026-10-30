//! A prescribed uniform background intervention on the original gas memories.
//!
//! The fixed-volume physical coordinates are compression, momentum and entropy
//! perturbations. For each native radial mode, t_n=y_n+a_n*x_p is temperature
//! in energy coordinates; a_n=sqrt(w_n*R/cv_inf). Compression is proportional
//! to (x_p-sum a_n*P_n^T*y_n)/sqrt(B_inf), NOT x_p/sqrt(B_inf); P_n is the
//! original exact-mass or positive spatial-quadrature functional. Translational entropy
//! density is sqrt(rho*cp_inf/T)*y_n and reservoir entropy density is
//! sqrt(rho*C_j/T)*z_nj. Preserving these quantities gives the affine pressure
//! correction below. Simply scaling pressure would inject compression when
//! heat capacities change. Every actual quadratic storage change is external
//! parameter work, distinct from molecular loss and gas-to-wall heat.

use super::super::gas_update::{energy, scaled};
use super::*;
use fs_material::relaxation::RelaxingMedium;

impl DgDuct {
    /// Change a continuing uniform molecular gas without discarding its thermal
    /// or internal-temperature histories.
    ///
    /// The complete original molecular declaration is evaluated on the new
    /// source-backed field. Species identities AND fractions, reservoir support,
    /// geometry, exact-mass basis, native circular radial sections and their
    /// ordering must survive. Each thermal and molecular mode preserves its
    /// physical entropy perturbation; total weak compression and every axial
    /// or viscous-arm momentum are preserved. Pressure consequently includes
    /// the heat-capacity-dependent correction from ALL thermal histories.
    /// Structural motion, viscosity and Maxwell histories remain unchanged.
    ///
    /// Supports the original exact-mass circular diffusion, positive spatial
    /// wall quadrature on variable geometry, and adiabatic molecular gas.
    /// Spatial sections retain their exact original point, radius and radial
    /// support. The BACKGROUND GAS must remain uniform; a nonuniform field
    /// refuses without averaging. Every new local diffusion must pass the
    /// original angular band's complex AND loss admission.
    ///
    /// This is a zero-duration, prescribed background intervention with no
    /// injected acoustic entropy or compression. It does not solve background
    /// heat/species/inventory transport or supply mean pressure wall
    /// traction. Its external owner owes the reported signed parameter work.
    /// Original sample count and cumulative work/loss survive. Rebuild the
    /// propagator before continuing. Every refusal leaves both owners intact.
    pub fn rebind_uniform_molecular_gas_field(
        &mut self,
        state: &mut DgState,
        field: &AxialGasField,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, GasFieldDuctError> {
        poll(cancelled)?;
        if field.length_m().to_bits() != self.total_length.to_bits()
            || !maximum_relative_error.is_finite()
            || !(0.0 < maximum_relative_error && maximum_relative_error <= 0.5)
            || !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
        {
            return Err(invalid("molecular refresh requires the original length and finite accuracy/work allowances").into());
        }
        let original_field = self
            .gas_field
            .as_ref()
            .ok_or_else(|| invalid("molecular refresh requires the complete original gas field"))?;
        let old_resolved = original_field
            .uniform_state()
            .ok_or_else(|| invalid("molecular refresh cannot average a nonuniform original gas"))?;
        let new_resolved = field
            .uniform_state()
            .ok_or_else(|| invalid("molecular refresh requires a uniform replacement gas field"))?;
        let old_medium = RelaxingMedium::from_thermochemical(old_resolved)?;
        let new_medium = RelaxingMedium::from_thermochemical(new_resolved)?;
        if old_medium.components() != new_medium.components() {
            return Err(invalid("molecular refresh preserves the complete species fractions; acoustic species transport is not represented").into());
        }
        let previous = self
            .molecular
            .as_ref()
            .ok_or_else(|| invalid("molecular refresh requires original molecular histories"))?;
        let old_material = previous
            .state
            .as_ref()
            .ok_or_else(|| invalid("molecular refresh requires one uniform evaluated material"))?;
        let law = previous
            .field_law
            .as_ref()
            .ok_or_else(|| {
                invalid("molecular refresh needs the retained original constitutive declaration")
            })?
            .clone();
        if law.evaluate(&old_medium)? != *old_material || old_material.gas() != &self.gas {
            return Err(
                invalid("molecular refresh lost its exact original material binding").into(),
            );
        }
        let replacement = law.evaluate(&new_medium)?;
        let same_support = old_material.processes().len() == replacement.processes().len()
            && old_material
                .processes()
                .iter()
                .zip(replacement.processes())
                .all(|(old, new)| {
                    old.species == new.species
                        && old.index == new.index
                        && (old.heat_capacity_j_kg_k > 0.0) == (new.heat_capacity_j_kg_k > 0.0)
                });
        if !same_support {
            return Err(
                invalid("molecular refresh cannot create, discard or relabel a reservoir").into(),
            );
        }
        let m = self.field_dofs;
        let n = self.total_dofs();
        if state.operator_identity != self.operator_identity
            || state.tilde.len() != n
            || state.p.len() != m
            || state.q.len() != m
            || state.tilde.iter().any(|v| !v.is_finite())
            || !state.cumulative_parameter_work_j.is_finite()
        {
            return Err(
                invalid("molecular refresh needs its original finite acoustic state").into(),
            );
        }

        // Zero-capacity declarations have no additional coordinates. Preserve
        // the existing nonmolecular transfer and its arithmetic exactly, while
        // retaining the source law and every explicitly inactive reservoir.
        if previous.modes.is_empty() {
            let mut candidate = self.clone();
            candidate.molecular = None;
            let mut next = state.clone();
            let report = match &candidate.viscothermal {
                Some(ViscothermalMemory::Uniform(_)) => candidate
                    .rebind_uniform_viscothermal_gas_field(
                        &mut next,
                        field,
                        maximum_relative_error,
                        maximum_absolute_parameter_work_j,
                        cancelled,
                    )?,
                None => candidate.rebind_uniform_gas_field(
                    &mut next,
                    field,
                    maximum_absolute_parameter_work_j,
                    cancelled,
                )?,
                Some(ViscothermalMemory::Spatial(_)) => candidate
                    .rebind_spatial_viscothermal_gas_field(
                        &mut next,
                        field,
                        maximum_relative_error,
                        maximum_absolute_parameter_work_j,
                        cancelled,
                    )?,
            };
            candidate.molecular = Some(MolecularMemory {
                state: Some(replacement),
                field_law: Some(law),
                bulk: candidate.pressure_bulk_modulus(),
                modes: Vec::new(),
            });
            poll(cancelled)?;
            *self = candidate;
            *state = next;
            return Ok(report);
        }

        let old_gas = self.gas;
        let new_gas = *replacement.gas();
        let cp_old = old_material.frozen_cv_j_kg_k() + old_gas.specific_gas_constant;
        let cp_new = replacement.frozen_cv_j_kg_k() + new_gas.specific_gas_constant;
        let speed_new = replacement.frozen_sound_speed_m_s();
        let bulk_new = new_gas.density * speed_new * speed_new;
        let pressure_scale = det::sqrt(bulk_new) / det::sqrt(previous.bulk);
        let momentum_scale = det::sqrt(old_gas.density) / det::sqrt(new_gas.density);
        let thermal_entropy_scale = entropy_scale(&old_gas, cp_old, &new_gas, cp_new)?;
        let reservoir_scales: Vec<_> = old_material
            .processes()
            .iter()
            .zip(replacement.processes())
            .filter(|(old, _)| old.heat_capacity_j_kg_k > 0.0)
            .map(|(old, new)| {
                entropy_scale(
                    &old_gas,
                    old.heat_capacity_j_kg_k,
                    &new_gas,
                    new.heat_capacity_j_kg_k,
                )
            })
            .collect::<Result<_, _>>()?;
        if [pressure_scale, momentum_scale]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(
                invalid("molecular refresh has unrepresentable physical storage scales").into(),
            );
        }

        // Reuse the original geometry/mass/wall/material assembly. Radial
        // section identity is the unchanged dimensionless Jacobi problem,
        // not a nearest-rate match between fitted state-space coordinates.
        let mut candidate =
            DgDuct::new_with_gas_field(&self.geometry, field, self.mesh, self.inlet, self.outlet)?
                .with_side_ports(&self.side_ports)?;
        match &self.viscothermal {
            Some(ViscothermalMemory::Uniform(model)) => {
                let radius = self.uniform_radius_m.ok_or_else(|| {
                    invalid("uniform molecular diffusion lost its original circular radius")
                })?;
                let count = model.series_terms().len();
                let rebound = fs_phs::ZkDiffusive::new(
                    radius, new_gas.density, new_gas.dynamic_viscosity, new_gas.gamma,
                    new_gas.prandtl, count, model.band_rad_s(), maximum_relative_error,
                ).map_err(|_| invalid("replacement molecular gas does not resolve the original Bessel diffusion band"))?;
                if rebound.series_terms().len() != count
                    || rebound.shunt_terms().len() != model.shunt_terms().len()
                {
                    return Err(invalid("molecular refresh changed native radial support").into());
                }
                candidate = candidate.with_uniform_viscothermal_model(rebound)?;
            }
            None => {}
            Some(ViscothermalMemory::Spatial(_)) => {
                let resolution = DgWallResolution {
                    maximum_relative_error,
                    ..self
                        .spatial_viscothermal_admission()
                        .expect("original spatial wall")
                        .resolution
                };
                candidate = candidate.with_spatial_viscothermal(resolution)?;
                let (
                    Some(ViscothermalMemory::Spatial(old)),
                    Some(ViscothermalMemory::Spatial(new)),
                ) = (&self.viscothermal, &candidate.viscothermal)
                else {
                    unreachable!("native spatial assembly")
                };
                let old_points = old.thermal_quadrature(m);
                let new_points = new.thermal_quadrature(m);
                if old_points.len() != new_points.len()
                    || old_points.iter().zip(&new_points).any(|(old, new)| {
                        old.point_m.to_bits() != new.point_m.to_bits()
                            || old.field_start != new.field_start
                            || old.terms.len() != new.terms.len()
                            || old.projection.len() != new.projection.len()
                    })
                {
                    return Err(invalid(
                        "molecular refresh changed an original spatial thermal section",
                    )
                    .into());
                }
            }
        }
        candidate = candidate.with_molecular_relaxation_controlled(replacement, cancelled)?;
        candidate
            .molecular
            .as_mut()
            .expect("original molecular assembly")
            .field_law = Some(law);
        let renewed = candidate
            .molecular
            .as_ref()
            .expect("original molecular assembly");
        let prefix = candidate.total_dofs();
        let viscous_end = previous.modes[0].start;
        if prefix > n
            || renewed.modes.len() != previous.modes.len()
            || previous.modes.iter().zip(&renewed.modes).any(|(old, new)| {
                old.start != new.start
                    || old.rates.len() != new.rates.len()
                    || old.rates.len() != reservoir_scales.len()
                    || !same_projection_support(&old.projection, &new.projection)
            })
            || viscous_end < 2 * m
            || previous.modes.last().is_none_or(|mode| {
                mode.start + mode.projection.width(m) * (mode.rates.len() + 1) != prefix
            })
            || candidate.elements.len() != self.elements.len()
            || self
                .elements
                .iter()
                .zip(&candidate.elements)
                .any(|(old, new)| {
                    old.x_left.to_bits() != new.x_left.to_bits()
                        || old.length.to_bits() != new.length.to_bits()
                        || old.nodes_x != new.nodes_x
                        || old.areas != new.areas
                })
        {
            return Err(invalid(
                "molecular refresh changed the original physical coordinate support",
            )
            .into());
        }

        let mut next = state.clone();
        for i in m..viscous_end {
            next.tilde[i] = scaled(state.tilde[i], momentum_scale)?;
        }
        for (old, new) in previous.modes.iter().zip(&renewed.modes) {
            poll(cancelled)?;
            let width = old.projection.width(m);
            for i in 0..width {
                next.tilde[new.start + i] =
                    scaled(state.tilde[old.start + i], thermal_entropy_scale)?;
                for (j, &factor) in reservoir_scales.iter().enumerate() {
                    let at = old.start + (j + 1) * width + i;
                    next.tilde[at] = scaled(state.tilde[at], factor)?;
                }
            }
        }
        for i in 0..m {
            next.tilde[i] = scaled(state.tilde[i], pressure_scale)?;
        }
        for (old, new) in previous.modes.iter().zip(&renewed.modes) {
            poll(cancelled)?;
            for i in 0..old.projection.width(m) {
                let old_entries = old.projection.entries(i, old.pressure_scale);
                let new_entries = new.projection.entries(i, new.pressure_scale);
                for ((row, a_old), (_, a_new)) in old_entries.into_iter().zip(new_entries) {
                    // Difference form retains exact no-op/transport-only histories.
                    next.tilde[row] += a_new * next.tilde[new.start + i]
                        - pressure_scale * a_old * state.tilde[old.start + i];
                }
            }
        }
        if next.tilde[..m].iter().any(|v| !v.is_finite()) {
            return Err(invalid("molecular compression transfer is not finite").into());
        }

        if let Some(old_wall) = &self.compliant_wall {
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
            for (_, port) in &mut candidate.ports {
                port.resize(n, 0.0);
            }
        } else if prefix != n {
            return Err(invalid("molecular refresh found unowned acoustic coordinates").into());
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
            .map(|(before, after)| 0.5 * (after - before) * (after + before))
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
            || next.p.iter().chain(&next.q).any(|value| !value.is_finite())
            || !parameter_work_j.is_finite()
            || parameter_work_j.abs() > maximum_absolute_parameter_work_j
            || !accumulated.is_finite()
            || (parameter_work_j != 0.0 && accumulated == state.cumulative_parameter_work_j)
            || (state.cumulative_parameter_work_j != 0.0 && accumulated == parameter_work_j)
        {
            return Err(invalid(
                "molecular parameter work exceeded its allowance or finite physical resolution",
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

fn same_projection_support(old: &PressureProjection, new: &PressureProjection) -> bool {
    match (old, new) {
        (PressureProjection::Uniform, PressureProjection::Uniform) => true,
        (
            PressureProjection::Quadrature {
                first: a,
                weights: x,
            },
            PressureProjection::Quadrature {
                first: b,
                weights: y,
            },
        ) => a == b && x.len() == y.len(),
        _ => false,
    }
}

fn entropy_scale(
    old: &GasState,
    capacity_old: f64,
    new: &GasState,
    capacity_new: f64,
) -> Result<f64, DuctError> {
    // Square roots of separate positive factors avoid forming a large product
    // before taking its square root. No relaxation time enters storage.
    let scale = (det::sqrt(old.density) / det::sqrt(new.density))
        * (det::sqrt(capacity_old) / det::sqrt(capacity_new))
        * (det::sqrt(new.temperature) / det::sqrt(old.temperature));
    if !scale.is_finite() || scale <= 0.0 {
        Err(invalid(
            "molecular entropy transfer has an unrepresentable storage scale",
        ))
    } else {
        Ok(scale)
    }
}

#[cfg(test)]
mod tests;
