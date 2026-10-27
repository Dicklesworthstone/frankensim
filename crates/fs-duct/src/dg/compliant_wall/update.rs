//! Live material changes on the original reciprocal gas/wall state.
//!
//! Fixed geometry makes each wall energy coordinate proportional to sqrt(K) w
//! or sqrt(sigma) v, with the SAME geometric Gram root. Preserve displacement
//! and conjugate surface momentum sigma*v, not normalized amplitudes. The
//! Maxwell arm z=K_j*(w-alpha_j) retains its elastic extension z/K_j and
//! viscous offset alpha_j. The complete storage change belongs to external
//! parameter work, with no instantaneous dashpot loss.

use super::*;

/// A zero-duration material intervention on an already vibrating wall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgWallMaterialUpdate {
    /// Total acoustic, viscothermal and wall energy before the update [J].
    pub energy_before_j: f64,
    /// Total energy after preserving wall displacement and momentum [J].
    pub energy_after_j: f64,
    /// Original wall kinetic, equilibrium and Maxwell spring storage [J].
    pub wall_energy_before_j: f64,
    /// Updated wall kinetic, equilibrium and Maxwell spring storage [J].
    pub wall_energy_after_j: f64,
    /// Signed work INTO the acoustic/wall system [J]. The parameter-changing
    /// owner owes this energy; it is neither dissipated heat nor port work.
    pub parameter_work_j: f64,
}

/// Explicit physical correspondence of one local Maxwell history.
/// The source/material owner must establish arm identity; this is not inferred
/// from sorted relaxation times, equal coefficients or storage-vector positions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DgWallRelaxationTransfer {
    /// Original element whose complete nodal history is transferred.
    pub element: usize,
    /// Arm index in that element's OLD supplied GeneralizedMaxwell.terms,
    /// including any exact zero arms that have no stored coordinates.
    pub previous_arm: usize,
    /// Corresponding index in the replacement element's source spectrum.
    pub replacement_arm: usize,
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
        self.rebind_wall_data(
            state,
            pins,
            None,
            &[],
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    /// Update a live Maxwell wall with a complete, explicit local arm mapping.
    ///
    /// The original source owner declares which arm remains the SAME physical
    /// internal variable on each fixed element. Every positive old and new arm
    /// must be covered exactly once on every node of its support; zero/positive
    /// transitions, missing/duplicate histories and moved supports refuse.
    /// Uniform full-field and regional storage layouts may differ: the transfer
    /// only reorders existing element blocks, without projection or averaging.
    ///
    /// Preserve w, surface momentum and every arm's elastic extension z_j/K_j.
    /// Thus z_new=(K_new/K_old)*z_old and alpha_j=w-z_j/K_j is unchanged.
    /// Each arm's normalized state scales by sqrt(K_new/K_old); its actual
    /// spring-energy change joins the existing signed parameter-work account.
    /// A relaxation-time-only change owes zero instantaneous work and changes
    /// future dissipation through the original native exponential. No history
    /// is reset and no second constitutive time step is run.
    ///
    /// Gas/geometry, clock, budget, atomicity and propagator-rebuild rules are
    /// those of `rebind_compliant_wall_regions`. This raw coefficient API does
    /// not invent source identities: source-backed consumers must verify their
    /// immutable material/claim lineage before supplying `correspondence`.
    #[allow(clippy::too_many_arguments)]
    pub fn rebind_compliant_wall_relaxation_regions(
        &mut self,
        state: &mut DgState,
        pins: &[WallPin],
        spectra: &[Option<&fs_material::visco::GeneralizedMaxwell>],
        correspondence: &[DgWallRelaxationTransfer],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgWallMaterialUpdate, DuctError> {
        self.rebind_wall_data(
            state,
            pins,
            Some(spectra),
            correspondence,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn rebind_wall_data(
        &mut self,
        state: &mut DgState,
        pins: &[WallPin],
        spectra: Option<&[Option<&fs_material::visco::GeneralizedMaxwell>]>,
        correspondence: &[DgWallRelaxationTransfer],
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
        if spectra.is_none() && !wall.relaxation.is_empty() {
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
        if spectra.is_none()
            && pins
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
        let candidate = match spectra {
            Some(spectra) => candidate
                .with_compliant_wall_relaxation_regions_controlled(pins, spectra, cancelled)?,
            None => candidate.with_compliant_wall_regions_controlled(pins, cancelled)?,
        };
        let mut next = state.clone();
        next.tilde.resize(candidate.total_dofs(), 0.0);
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
        if spectra.is_some() {
            parameter_work_j += transfer_relaxation(
                wall,
                candidate
                    .compliant_wall
                    .as_ref()
                    .expect("replacement wall attached"),
                state,
                &mut next,
                self.mesh,
                correspondence,
                cancelled,
            )?;
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

// Each match covers the original local polynomial basis exactly. This handles
// uniform <-> regional storage changes without matching by flattened position.
#[allow(clippy::too_many_arguments)]
fn transfer_relaxation(
    old: &CompliantWall,
    new: &CompliantWall,
    before: &DgState,
    after: &mut DgState,
    mesh: DgMesh,
    correspondence: &[DgWallRelaxationTransfer],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, DuctError> {
    let local = mesh.order + 1;
    let memory_start = old.start + 2 * mesh.elements * local;
    if correspondence.len() > before.tilde.len() / local {
        return Err(invalid(
            "Maxwell correspondence exceeds the original history extent",
        ));
    }
    let mut old_seen = vec![false; before.tilde.len() - memory_start];
    let mut new_seen = vec![false; after.tilde.len() - memory_start];
    let find = |wall: &CompliantWall, element: usize, arm: usize| {
        wall.relaxation.iter().position(|branch| {
            branch.source_arm == arm
                && element * local >= branch.field_start
                && element * local + local <= branch.field_start + branch.dofs
        })
    };
    let mut work = 0.0;
    for entry in correspondence {
        check_cancelled(cancelled)?;
        if entry.element >= mesh.elements {
            return Err(invalid(
                "Maxwell history correspondence is outside the original elements",
            ));
        }
        let a = find(old, entry.element, entry.previous_arm).ok_or_else(|| {
            invalid("Maxwell correspondence names an absent original positive arm")
        })?;
        let b = find(new, entry.element, entry.replacement_arm).ok_or_else(|| {
            invalid("Maxwell correspondence names an absent replacement positive arm")
        })?;
        let (a, b) = (&old.relaxation[a], &new.relaxation[b]);
        let scale = det::sqrt(b.stiffness_per_area) / det::sqrt(a.stiffness_per_area);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(invalid(
                "Maxwell extension transfer scale is unrepresentable",
            ));
        }
        let old_start = a.start + entry.element * local - a.field_start;
        let new_start = b.start + entry.element * local - b.field_start;
        for i in 0..local {
            let (old_index, new_index) = (old_start + i, new_start + i);
            if old_seen[old_index - memory_start] || new_seen[new_index - memory_start] {
                return Err(invalid(
                    "Maxwell history correspondence duplicates an internal variable",
                ));
            }
            old_seen[old_index - memory_start] = true;
            new_seen[new_index - memory_start] = true;
            let old_value = before.tilde[old_index];
            let new_value = old_value * scale;
            if !new_value.is_finite() || (old_value != 0.0 && new_value == 0.0) {
                return Err(invalid(
                    "Maxwell elastic extension was lost during state transfer",
                ));
            }
            after.tilde[new_index] = new_value;
            work += 0.5 * (new_value - old_value) * (new_value + old_value);
        }
    }
    if old_seen.iter().chain(&new_seen).any(|seen| !seen) {
        return Err(invalid(
            "Maxwell update must preserve every original and replacement history",
        ));
    }
    Ok(work)
}
