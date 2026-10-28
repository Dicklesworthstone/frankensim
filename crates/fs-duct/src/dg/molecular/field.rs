//! Pointwise molecular thermodynamics on the original primitive gas field.
//!
//! Let E=Mp_eq^(1/2), F=Mp_inf^(-1/2), T=F E for each original element.
//! The retained generator and flow ports use the congruence T L T^T and T B,
//! not a similarity or element-average sound speed. A local pressure channel
//! becomes sqrt(C_inf(x)/C_eq(x)) T v_eq = sqrt(w C_inf(x)) F phi(x).
//! Its heat and molecular reservoirs use that point's material cp and rates.
//! The exact mass and the original positive wall quadrature remain separate
//! refinement axes. At finite quadrature the effective compliance is the
//! exact frozen mass plus the quadrature of the local dynamic correction.

use super::*;
use fs_material::relaxation::RelaxingMedium;

fn evaluate(
    field: &AxialGasField,
    law: &RelaxationAbsorption,
    x: f64,
) -> Result<RelaxationState, GasFieldDuctError> {
    Ok(law.evaluate(&RelaxingMedium::from_thermochemical(&field.state_at(x)?)?)?)
}

fn active(state: &RelaxationState) -> usize {
    state
        .processes()
        .iter()
        .filter(|p| p.heat_capacity_j_kg_k > 0.0)
        .count()
}

fn add_states(total: usize, modes: usize, reservoirs: usize) -> Result<usize, DuctError> {
    modes
        .checked_mul(reservoirs)
        .and_then(|extra| total.checked_add(extra))
        .filter(|n| *n <= MAX_DG_PROPAGATOR_DOFS)
        .ok_or_else(|| invalid("molecular plus thermal gas exceeds the DG state budget"))
}

impl DgDuct {
    /// Total cold state count with the original field and declared molecular
    /// law, before allocating the augmented dense generator. Count only active
    /// local reservoirs at the actual original thermal quadrature points.
    /// Requires the original source-backed field; a varying field additionally
    /// requires spatial viscothermal memory. Uniform fields use the existing
    /// exact-mass or spatial route. Repeated or post-structural binding refuses.
    /// Material evaluation checks every field station and wall quadrature point;
    /// construction additionally evaluates the original adaptive mass rule.
    pub fn molecular_relaxation_field_dofs(
        &self,
        law: &RelaxationAbsorption,
    ) -> Result<usize, GasFieldDuctError> {
        self.molecular_relaxation_field_dofs_controlled(law, &mut || false)
    }

    /// Cancellable count-only counterpart; no accepted operator is changed.
    pub fn molecular_relaxation_field_dofs_controlled(
        &self,
        law: &RelaxationAbsorption,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, GasFieldDuctError> {
        poll(cancelled)?;
        if self.molecular.is_some() || self.compliant_wall.is_some() {
            return Err(
                invalid("attach molecular gas once, before compliant wall mechanics").into(),
            );
        }
        let field = self.gas_field.as_ref().ok_or_else(|| {
            invalid("pointwise molecular law requires the original source-backed gas field")
        })?;
        for &x in field.stations_m() {
            poll(cancelled)?;
            evaluate(field, law, x)?;
        }
        if field.uniform_state().is_some() {
            let count = active(&evaluate(field, law, 0.0)?);
            if count == 0 {
                return Ok(self.total_dofs());
            }
            return Ok(match &self.viscothermal {
                Some(ViscothermalMemory::Uniform(model)) => add_states(
                    self.total_dofs(),
                    self.field_dofs * model.shunt_terms().len(),
                    count,
                )?,
                Some(ViscothermalMemory::Spatial(wall)) => {
                    let modes = wall
                        .thermal_quadrature(self.field_dofs)
                        .iter()
                        .map(|p| p.terms.len())
                        .sum();
                    add_states(self.total_dofs(), modes, count)?
                }
                None => add_states(self.total_dofs(), self.field_dofs, count + 1)?,
            });
        }
        let Some(ViscothermalMemory::Spatial(wall)) = &self.viscothermal else {
            return Err(invalid(
                "nonuniform molecular gas requires the original spatial thermal quadrature",
            )
            .into());
        };
        let mut total = self.total_dofs();
        for point in wall.thermal_quadrature(self.field_dofs) {
            poll(cancelled)?;
            total = add_states(
                total,
                point.terms.len(),
                active(&evaluate(field, law, point.point_m)?),
            )?;
        }
        Ok(total)
    }

    /// Attach one declared molecular law to the complete original gas field.
    /// Each point resolves its actual primitive temperature, composition and
    /// material sources. No inlet state, average cp, or average rate is used.
    /// This is cold, frozen-background construction; species/heat transport and
    /// live reservoir transfer remain separate, unsupported operations.
    pub fn with_molecular_relaxation_field(
        self,
        law: RelaxationAbsorption,
    ) -> Result<Self, GasFieldDuctError> {
        self.with_molecular_relaxation_field_controlled(law, &mut || false)
    }

    /// Cancellable counterpart of [`Self::with_molecular_relaxation_field`].
    /// The original local-radius Bessel measure, viscosity histories, element
    /// basis, physical ports and endpoint boundary impedances are retained.
    /// Pressure mass is integrated to the existing adaptive mass tolerance;
    /// relaxation corrections use the declared wall quadrature resolution.
    /// An all-zero law preserves the original numerical operator bitwise.
    pub fn with_molecular_relaxation_field_controlled(
        mut self,
        law: RelaxationAbsorption,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasFieldDuctError> {
        let total = self.molecular_relaxation_field_dofs_controlled(&law, cancelled)?;
        let field = self
            .gas_field
            .as_ref()
            .expect("field count validated original field")
            .clone();
        if field.uniform_state().is_some() {
            return Ok(
                self.with_molecular_relaxation_controlled(evaluate(&field, &law, 0.0)?, cancelled)?
            );
        }
        let inlet = evaluate(&field, &law, 0.0)?;
        // Use the original mass integrator and element boundaries. The temporary
        // sampling record supplies only rho and c_inf to element_masses; it is
        // never exposed as a new equilibrium GasState or material receipt.
        let local = self.mesh.order + 1;
        let m = self.field_dofs;
        let mut maps = Vec::with_capacity(self.elements.len());
        let mut roots = Vec::with_capacity(self.elements.len());
        let mut segment = 0;
        let mut start = 0.0;
        let mut sampled_active = false;
        for &x in field.stations_m() {
            poll(cancelled)?;
            sampled_active |= active(&evaluate(&field, &law, x)?) > 0;
        }
        for element in &self.elements {
            poll(cancelled)?;
            while segment + 1 < self.geometry.segments.len()
                && element.x_left + 0.5 * element.length
                    >= start + self.geometry.segments[segment].axial_length()
            {
                start += self.geometry.segments[segment].axial_length();
                segment += 1;
            }
            let span = AxialSpan {
                x_start: start,
                length: self.geometry.segments[segment].axial_length(),
                geometry_start: 0.0,
                segment: self.geometry.segments[segment].clone(),
            };
            let (mass, _) = element_masses(
                &span,
                element.x_left - start,
                element.length,
                &self.lgl,
                &mut |x| -> Result<GasState, GasFieldDuctError> {
                    poll(cancelled)?;
                    let state = evaluate(&field, &law, x)?;
                    sampled_active |= active(&state) > 0;
                    let mut gas = *state.gas();
                    gas.sound_speed = state.frozen_sound_speed_m_s();
                    Ok(gas)
                },
            )?;
            let (root, inverse) = spd_sqrt_and_inv_sqrt(&mass, local)?;
            let mut map = vec![0.0; local * local];
            for i in 0..local {
                for j in 0..local {
                    map[i * local + j] = (0..local)
                        .map(|k| inverse[i * local + k] * element.mp_sqrt[k * local + j])
                        .sum();
                }
            }
            if map.iter().any(|v| !v.is_finite()) {
                return Err(
                    invalid("pointwise frozen-pressure mass map is unrepresentable").into(),
                );
            }
            maps.push(map);
            roots.push((root, inverse));
        }
        if total == self.total_dofs() {
            // A coarse wall rule cannot silently erase relaxation seen by the
            // original stations or independently refined pressure mass rule.
            if sampled_active {
                return Err(invalid(
                    "thermal quadrature misses active molecular support; refine its axial points",
                )
                .into());
            }
            self.molecular = Some(MolecularMemory {
                state: None,
                field_law: Some(law),
                bulk: self.gas.density * self.gas.sound_speed.powi(2),
                modes: Vec::new(),
            });
            poll(cancelled)?;
            return Ok(self);
        }
        let old_total = self.total_dofs();
        let Some(ViscothermalMemory::Spatial(wall)) = &mut self.viscothermal else {
            unreachable!()
        };
        let (kept, channels) = wall.split_for_molecular(m);
        let prefix = kept.len();
        let mut cursor = prefix;
        let mut modes = Vec::new();
        for point in channels {
            poll(cancelled)?;
            let state = evaluate(&field, &law, point.point_m)?;
            let gas = state.gas();
            let cv = state.frozen_cv_j_kg_k();
            let cp = cv + gas.specific_gas_constant;
            let ratio = gas.sound_speed / state.frozen_sound_speed_m_s();
            let map = &maps[point.field_start / local];
            let projection: Vec<f64> = (0..local)
                .map(|i| {
                    ratio
                        * (0..local)
                            .map(|j| map[i * local + j] * point.projection[j])
                            .sum::<f64>()
                })
                .collect();
            let reservoirs: Vec<_> = state
                .processes()
                .iter()
                .filter(|p| p.heat_capacity_j_kg_k > 0.0)
                .collect();
            let scales: Vec<_> = reservoirs
                .iter()
                .map(|p| det::sqrt(p.heat_capacity_j_kg_k / cp))
                .collect();
            let rates: Vec<_> = reservoirs
                .iter()
                .map(|p| 1.0 / p.relaxation_time_s)
                .collect();
            let weight_sum: f64 = point
                .terms
                .iter()
                .map(|t| t.gain_per_s / (t.rate_per_s * (gas.gamma - 1.0)))
                .sum();
            if point.terms.is_empty()
                || (weight_sum - 1.0).abs() > 4e-12
                || projection.iter().any(|v| !v.is_finite())
                || scales
                    .iter()
                    .chain(&rates)
                    .any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(
                    invalid("pointwise molecular thermal measure is unrepresentable").into(),
                );
            }
            for term in point.terms {
                let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                let pressure_scale = det::sqrt(weight * gas.specific_gas_constant / cv);
                let heat_rate = term.rate_per_s * gas.specific_heat_cp / cp;
                if [pressure_scale, heat_rate]
                    .iter()
                    .any(|v| !v.is_finite() || *v <= 0.0)
                {
                    return Err(invalid(
                        "pointwise molecular thermal coefficients are unrepresentable",
                    )
                    .into());
                }
                modes.push(ThermalMode {
                    start: cursor,
                    pressure_scale,
                    heat_rate,
                    projection: PressureProjection::Quadrature {
                        first: point.field_start,
                        weights: projection.clone(),
                    },
                    reservoir_scales: scales.clone(),
                    rates: rates.clone(),
                });
                cursor = add_states(cursor, 1, reservoirs.len() + 1)?;
            }
        }
        if cursor != total {
            return Err(
                invalid("pointwise molecular state count changed during preparation").into(),
            );
        }
        let mut operator = vec![0.0; total * total];
        for i in 0..prefix {
            poll(cancelled)?;
            for j in 0..prefix {
                operator[i * total + j] = self.l_tilde[kept[i] * old_total + kept[j]];
            }
        }
        restore_pressure_boundaries(&self, &mut operator, total);
        // Left multiplication then right multiplication by the TRANSPOSE of
        // the same local map. One scratch vector avoids another dense matrix.
        let mut scratch = vec![0.0; local];
        for (e, map) in maps.iter().enumerate() {
            poll(cancelled)?;
            let first = e * local;
            for j in 0..prefix {
                for i in 0..local {
                    scratch[i] = (0..local)
                        .map(|k| map[i * local + k] * operator[(first + k) * total + j])
                        .sum();
                }
                for i in 0..local {
                    operator[(first + i) * total + j] = scratch[i];
                }
            }
        }
        for (e, map) in maps.iter().enumerate() {
            poll(cancelled)?;
            let first = e * local;
            for i in 0..prefix {
                for j in 0..local {
                    scratch[j] = (0..local)
                        .map(|k| operator[i * total + first + k] * map[j * local + k])
                        .sum();
                }
                for j in 0..local {
                    operator[i * total + first + j] = scratch[j];
                }
            }
        }
        for mode in &modes {
            poll(cancelled)?;
            stamp_mode(&mut operator, total, mode, m);
        }
        for (_, port) in &mut self.ports {
            poll(cancelled)?;
            port.resize(total, 0.0);
            for (e, map) in maps.iter().enumerate() {
                let first = e * local;
                for i in 0..local {
                    scratch[i] = (0..local)
                        .map(|j| map[i * local + j] * port[first + j])
                        .sum();
                }
                port[first..first + local].copy_from_slice(&scratch);
            }
            if port.iter().any(|v| !v.is_finite()) {
                return Err(invalid("pointwise molecular pressure port is unrepresentable").into());
            }
        }
        if operator.iter().any(|v| !v.is_finite()) {
            return Err(invalid("pointwise molecular generator is unrepresentable").into());
        }
        for (element, (root, inverse)) in self.elements.iter_mut().zip(roots) {
            element.mp_sqrt = root;
            element.mp_inv_sqrt = inverse;
        }
        self.l_tilde = operator;
        self.memory_dofs = total - 2 * m;
        self.g_block = None;
        self.molecular = Some(MolecularMemory {
            state: None,
            field_law: Some(law),
            bulk: inlet.gas().density * inlet.frozen_sound_speed_m_s().powi(2),
            modes,
        });
        self.operator_identity = hash_operator(
            self.mesh,
            self.inlet,
            self.outlet,
            &self.elements,
            &self.l_tilde,
            &self.ports,
        );
        poll(cancelled)?;
        Ok(self)
    }

    /// Original pointwise law for a nonuniform material field. Uniform fields
    /// report their evaluated state through [`Self::molecular_relaxation`].
    #[must_use]
    pub fn molecular_relaxation_field(&self) -> Option<&RelaxationAbsorption> {
        self.molecular.as_ref().and_then(|m| m.field_law.as_ref())
    }
}
