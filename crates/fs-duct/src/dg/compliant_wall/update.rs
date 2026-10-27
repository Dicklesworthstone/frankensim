//! Live material changes on the original reciprocal gas/wall state.
//!
//! Fixed geometry makes each wall energy coordinate proportional to sqrt(K) w
//! or sqrt(sigma) v, with the SAME geometric Gram root. Preserve displacement
//! and conjugate surface momentum sigma*v, not normalized amplitudes. The
//! resulting storage change belongs to the external parameter-work account.

use super::*;

/// A zero-duration material intervention on an already vibrating wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgWallMaterialUpdate {
    /// Total acoustic, viscothermal and wall energy before the update [J].
    pub energy_before_j: f64,
    /// Total energy after preserving wall displacement and momentum [J].
    pub energy_after_j: f64,
    /// Original wall kinetic and spring storage [J].
    pub wall_energy_before_j: f64,
    /// Updated wall kinetic and spring storage [J].
    pub wall_energy_after_j: f64,
    /// Signed work INTO the acoustic/wall system [J]. The parameter-changing
    /// owner owes this energy; it is neither dissipated heat nor port work.
    pub parameter_work_j: f64,
}

impl DgDuct {
    /// Change elastic/Kelvin--Voigt material coefficients on a live DG wall.
    ///
    /// One pin is required per ORIGINAL element, with geometry, gas, boundary
    /// ports and the viscothermal discretization unchanged. All gas wave and
    /// diffusion-memory coordinates are copied exactly. Radial displacement w
    /// and conjugate surface momentum sigma*v are preserved at every node:
    /// spring coordinates scale by sqrt(K_new/K_old), velocity coordinates by
    /// sqrt(sigma_old/sigma_new). No amplitude is normalized to keep a sound
    /// level or to hide material work. Viscosity changes only future loss.
    ///
    /// The measured spring/kinetic storage change is added ONCE to the state's
    /// separate cumulative parameter-work ledger, never to inferred dissipation
    /// or initial energy. `maximum_absolute_parameter_work_j` is a finite,
    /// nonnegative caller budget for this intervention, not a supplied heat.
    /// The external material/thermal owner must settle that signed work itself.
    /// No mass-transfer, thermal expansion, finite-temperature heat or elapsed
    /// physical time is inferred by this fixed-geometry state transfer.
    ///
    /// Existing propagators become stale when the operator changes and must be
    /// rebuilt. This method publishes only after the complete operator/state
    /// proposal, work budget and cancellation checks pass; callers needing an
    /// atomic propagator rebuild can prepare this update on private clones.
    /// Hereditary Maxwell walls refuse: no branch history is dropped or matched
    /// by array position. Mesh/geometry/gas changes require their own transfer.
    pub fn rebind_compliant_wall_regions(
        &mut self,
        state: &mut DgState,
        pins: &[WallPin],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgWallMaterialUpdate, DuctError> {
        check_cancelled(cancelled)?;
        if !maximum_absolute_parameter_work_j.is_finite()
            || maximum_absolute_parameter_work_j < 0.0
            || pins.len() != self.mesh.elements
            || state.operator_identity != self.operator_identity
            || state.tilde.len() != self.total_dofs()
            || state.tilde.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "live wall update needs the original finite state, one pin per element and a finite nonnegative work budget",
            ));
        }
        let wall = self.compliant_wall.as_ref().ok_or_else(|| {
            invalid("live material update requires an already attached compliant wall")
        })?;
        if !wall.relaxation.is_empty() {
            return Err(invalid(
                "live Maxwell material changes require a source-identified branch-history transfer",
            ));
        }
        let energy_before_j = state.stored_energy_j();
        let wall_energy_before_j = self.compliant_wall_state(state)?.stored_energy_j;
        if !energy_before_j.is_finite() || !state.cumulative_parameter_work_j.is_finite() {
            return Err(invalid("live wall energy account is not finite"));
        }
        for pin in pins {
            check_cancelled(cancelled)?;
            fs_phs::wall_specific_impedance(pin, 1.0)
                .map_err(|_| invalid("live wall needs passive finite local coefficients"))?;
        }
        if pins
            .iter()
            .enumerate()
            .all(|(element, pin)| *pin == wall.pin_at(element))
        {
            check_cancelled(cancelled)?;
            return Ok(DgWallMaterialUpdate {
                energy_before_j,
                energy_after_j: energy_before_j,
                wall_energy_before_j,
                wall_energy_after_j: wall_energy_before_j,
                parameter_work_j: 0.0,
            });
        }

        // The original gas and diffusion block precedes the attached mechanics.
        // Reuse it verbatim, including all original port vectors and source maps.
        let m = self.field_dofs;
        let prefix = wall.start;
        let n = self.total_dofs();
        let mut candidate = self.clone();
        let mut operator = vec![0.0; prefix * prefix];
        for row in 0..prefix {
            check_cancelled(cancelled)?;
            operator[row * prefix..(row + 1) * prefix]
                .copy_from_slice(&self.l_tilde[row * n..row * n + prefix]);
        }
        candidate.l_tilde = operator;
        candidate.compliant_wall = None;
        candidate.memory_dofs = prefix - 2 * m;
        for (_, port) in &mut candidate.ports {
            port.truncate(prefix);
        }
        let candidate = candidate.with_compliant_wall_regions_controlled(pins, cancelled)?;
        let mut next = state.clone();
        let local = self.mesh.order + 1;
        let mut parameter_work_j = 0.0;
        for (element, pin) in pins.iter().enumerate() {
            check_cancelled(cancelled)?;
            let old = wall.pin_at(element);
            // Ratios of square roots avoid overflowing the squared ratio first.
            let spring_scale =
                det::sqrt(pin.stiffness_per_area) / det::sqrt(old.stiffness_per_area);
            let momentum_scale = det::sqrt(old.surface_density) / det::sqrt(pin.surface_density);
            for (start, scale) in [(prefix, spring_scale), (prefix + m, momentum_scale)] {
                if !scale.is_finite() || scale <= 0.0 {
                    return Err(invalid("live wall coordinate scale is unrepresentable"));
                }
                for i in 0..local {
                    let index = start + element * local + i;
                    let before = state.tilde[index];
                    let after = before * scale;
                    if !after.is_finite() || (before != 0.0 && after == 0.0) {
                        return Err(invalid("live wall displacement or momentum was lost"));
                    }
                    next.tilde[index] = after;
                    // Difference the changed stores before adding the potentially
                    // much larger unchanged acoustic/diffusive background energy.
                    parameter_work_j += 0.5 * (after - before) * (after + before);
                }
            }
        }
        next.operator_identity = candidate.operator_identity;
        next.modal_anchor_p.clear();
        next.modal_anchor_q.clear();
        next.modal_anchor_dt_bits = None;
        next.unforced_steps = 0;
        let energy_after_j = next.stored_energy_j();
        let wall_energy_after_j = candidate.compliant_wall_state(&next)?.stored_energy_j;
        let accumulated = state.cumulative_parameter_work_j + parameter_work_j;
        if !energy_after_j.is_finite()
            || !parameter_work_j.is_finite()
            || parameter_work_j.abs() > maximum_absolute_parameter_work_j
            || !accumulated.is_finite()
            || (parameter_work_j != 0.0 && accumulated == state.cumulative_parameter_work_j)
            || (state.cumulative_parameter_work_j != 0.0 && accumulated == parameter_work_j)
        {
            return Err(invalid(
                "live wall parameter work exceeded its budget or resolution",
            ));
        }
        next.cumulative_parameter_work_j = accumulated;
        check_cancelled(cancelled)?;
        let _old_duct = core::mem::replace(self, candidate);
        let _old_state = core::mem::replace(state, next);
        Ok(DgWallMaterialUpdate {
            energy_before_j,
            energy_after_j,
            wall_energy_before_j,
            wall_energy_after_j,
            parameter_work_j,
        })
    }
}
