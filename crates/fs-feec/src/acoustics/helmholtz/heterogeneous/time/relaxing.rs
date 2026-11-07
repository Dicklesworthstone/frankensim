//! Full spatial pressure waves with causal, spatially varying wall admittance.
//!
//! The original forms implement
//! `M p_dot + K psi + C0 p + sum Cj qj = b`, `psi_dot = p`,
//! and `tau_j qj_dot + qj = trace(p)`. Each memory retains the COMPLETE
//! original geometric trace, including coordinates shared with adjacent
//! pressure-release faces and coordinates with zero residue coefficient.
//! No arm matrix is inverted, lumped, diagonalized or assumed definite.
//!
//! One midpoint elimination adds `h/(2 tau_j+h) Cj` to the original sparse
//! pressure factor and subtracts `2 tau_j/(2 tau_j+h) Cj qj0` from its load.
//! The effective Robin loss is only an algebraic check. Reported heat uses
//! the physical C0 and independently reconstructed memories; each memory
//! stores `tau_j qj^T Cj qj/2` and dissipates `qj^T Cj qj`.
//!
//! This is second-order midpoint, not an exact relaxation propagator. Large
//! h/tau can alternate endpoint memories even though the scheme is passive.
//! Resolve temporal relaxation and wave phase independently of the original
//! surface/material quadrature and spatial pressure order.

use super::super::boundary_admittance::relaxation::midpoint_weights;
use super::super::{BoundaryVelocitySource, SpatialRelaxingBoundary, VolumeVelocitySource};
use super::mechanical::{balance, paired_work};
use super::{
    BoundaryAdmittance, C64, Coo, Csr, Dd, HeterogeneousError, HeterogeneousHelmholtzFe,
    HeterogeneousTimeEnergy, HeterogeneousTimeLimits, HeterogeneousTimeState,
    PreparedHeterogeneousTime, RealVelocityTrace, ResistiveBoundary, SpatialResistiveBoundary,
    energy_change, invalid, poll, positive_form, relative,
};

/// Original sparse acoustic limits plus a cap on complete boundary memories.
#[derive(Clone, Copy, Debug)]
pub struct RelaxingTimeLimits {
    pub acoustic: HeterogeneousTimeLimits,
    /// Sum of original geometric trace DOFs over every patch/relaxation arm.
    /// Sparse residue forms and the existing acoustic factor are additional.
    pub max_memory_dofs: usize,
}

impl Default for RelaxingTimeLimits {
    fn default() -> Self {
        Self {
            acoustic: HeterogeneousTimeLimits::default(),
            max_memory_dofs: 16384,
        }
    }
}

struct PreparedArm {
    patch: usize,
    term: usize,
    dofs: Vec<usize>,
    conductance: Csr,
    tau: f64,
    history_weight: f64,
    pressure_weight: f64,
}

/// One cached sparse acoustic factor and every original boundary-memory form.
pub struct PreparedRelaxingTime<'m, 'f, 'c> {
    condensed: PreparedHeterogeneousTime<'m, 'f, 'c>,
    instantaneous: Csr,
    arms: Vec<PreparedArm>,
}

/// Physical fluid and boundary memories bound to this exact preparation.
#[derive(Clone)]
pub struct RelaxingTimeState<'p, 'm, 'f, 'c> {
    owner: &'p PreparedRelaxingTime<'m, 'f, 'c>,
    acoustic: HeterogeneousTimeState<'m, 'f, 'c>,
    memories: Vec<Vec<f64>>,
}

impl core::fmt::Debug for RelaxingTimeState<'_, '_, '_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("RelaxingTimeState")
            .field("acoustic", &self.acoustic)
            .field("relaxation_terms", &self.memories.len())
            .finish_non_exhaustive()
    }
}

impl RelaxingTimeState<'_, '_, '_, '_> {
    /// Original acoustic state for existing point and surface probes.
    pub fn acoustic_state(&self) -> &HeterogeneousTimeState<'_, '_, '_> {
        &self.acoustic
    }
    pub fn pressure_pa(&self) -> &[f64] {
        self.acoustic.pressure_pa()
    }
    pub fn pressure_potential_pa_s(&self) -> &[f64] {
        self.acoustic.pressure_potential_pa_s()
    }
    /// Memory pressure q_j [Pa] in the prepared arm's original trace order.
    pub fn memory_pressure_pa(&self, arm: usize) -> Option<&[f64]> {
        self.memories.get(arm).map(Vec::as_slice)
    }
}

/// Actual physical acoustic and boundary storage [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxingTimeEnergy {
    pub acoustic: HeterogeneousTimeEnergy,
    pub boundary_memory_j: f64,
}

impl RelaxingTimeEnergy {
    pub fn total_j(self) -> f64 {
        self.acoustic.total_j() + self.boundary_memory_j
    }
}

/// Independent storage, heat and source-work accounts over one time window.
#[derive(Clone, Copy, Debug)]
pub struct RelaxingTimeWork {
    pub initial_energy_j: f64,
    pub final_energy_j: f64,
    pub delta_energy_j: f64,
    pub acoustic_delta_energy_j: f64,
    pub memory_delta_energy_j: f64,
    pub input_work_j: f64,
    /// Subtotal already included in input_work_j.
    pub boundary_source_work_j: f64,
    /// Subtotal already included in input_work_j.
    pub volume_source_work_j: f64,
    /// Heat h p_mid^T C0 p_mid of the original instantaneous laws.
    pub instantaneous_dissipation_j: f64,
    /// Heat h sum q_mid^T Cj q_mid of the original relaxation branches.
    pub relaxation_dissipation_j: f64,
    /// Signed transfer h sum p_mid^T Cj q_mid FROM fluid TO memories.
    /// May be negative when a previously excited boundary returns energy.
    pub memory_work_j: f64,
    pub backward_error: f64,
}

impl RelaxingTimeWork {
    pub fn balance_residual_j(self) -> f64 {
        self.delta_energy_j + self.instantaneous_dissipation_j + self.relaxation_dissipation_j
            - self.input_work_j
    }
    pub fn acoustic_balance_residual_j(self) -> f64 {
        self.acoustic_delta_energy_j + self.instantaneous_dissipation_j + self.memory_work_j
            - self.input_work_j
    }
    pub fn memory_balance_residual_j(self) -> f64 {
        self.memory_delta_energy_j + self.relaxation_dissipation_j - self.memory_work_j
    }
}

#[derive(Clone, Debug)]
pub struct RelaxingTimeStep<'p, 'm, 'f, 'c> {
    pub state: RelaxingTimeState<'p, 'm, 'f, 'c>,
    pub work: RelaxingTimeWork,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare causal spatial wall absorption with ONE original sparse factor.
    ///
    /// Uniform, memoryless spatial, relaxing and pressure-release patches must
    /// have disjoint faces. Internal components of one relaxing patch share its
    /// full pressure trace. Boundary and volume sources retain their original
    /// projected profiles and may drive any non-release face. Material fields
    /// are sampled before this call; no warm factorization or sampling recurs.
    /// These boundary memories are not combined with mechanical-wall states by
    /// this owner, and the midpoint condensate is not an exact-time material law.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_relaxing_boundaries<'m>(
        &'m self,
        duration_s: f64,
        uniform: &[ResistiveBoundary],
        spatial: &[&SpatialResistiveBoundary<'_, '_, '_>],
        relaxing: &[&SpatialRelaxingBoundary<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: RelaxingTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedRelaxingTime<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let mut effective = Vec::with_capacity(relaxing.len());
        let mut arms = Vec::new();
        let mut memory_dofs = 0usize;
        for (patch, law) in relaxing.iter().enumerate() {
            poll(cancelled)?;
            law.instantaneous.validate(self)?;
            // Geometric support is independent of the residue. Zero/semidefinite
            // coefficients do not authorize dropping original memory coordinates.
            let (trace, _) = self.admittance_matrices(
                &[BoundaryAdmittance {
                    faces: law.faces().to_vec(),
                    admittance_m_pa_s: C64::ONE,
                }],
                cancelled,
            )?;
            let dofs = (0..self.ndof())
                .filter(|&row| {
                    let (cols, values) = trace.row(row);
                    cols.binary_search(&row).is_ok_and(|i| values[i] > 0.0)
                })
                .collect::<Vec<_>>();
            memory_dofs = dofs
                .len()
                .checked_mul(law.terms.len())
                .and_then(|count| memory_dofs.checked_add(count))
                .ok_or_else(|| invalid("spatial admittance memory extent overflow"))?;
            if memory_dofs > limits.max_memory_dofs {
                return Err(HeterogeneousError::Limit {
                    what: "spatial admittance memory degrees of freedom",
                    requested: memory_dofs,
                    maximum: limits.max_memory_dofs,
                });
            }
            for (term, original) in law.terms.iter().enumerate() {
                poll(cancelled)?;
                let (history_weight, pressure_weight) =
                    midpoint_weights(duration_s, original.relaxation_time_s)?;
                let mut local = Coo::new(dofs.len(), dofs.len());
                for (row, &global_row) in dofs.iter().enumerate() {
                    poll(cancelled)?;
                    let (cols, entries) = original.conductance.row(global_row);
                    for (&col, &entry) in cols.iter().zip(entries) {
                        if entry != 0.0 {
                            let col = dofs.binary_search(&col).map_err(|_| {
                                invalid("incomplete original admittance trace support")
                            })?;
                            local.push(row, col, entry);
                        }
                    }
                }
                arms.push(PreparedArm {
                    patch,
                    term,
                    dofs: dofs.clone(),
                    conductance: local.assemble(),
                    tau: original.relaxation_time_s,
                    history_weight,
                    pressure_weight,
                });
            }
            effective.push(law.midpoint_condensate(duration_s, cancelled)?);
        }
        let effective_fields = spatial
            .iter()
            .copied()
            .chain(effective.iter())
            .collect::<Vec<_>>();
        let condensed = self.prepare_time_with_spatial_boundaries(
            duration_s,
            uniform,
            &effective_fields,
            pressure_release_faces,
            limits.acoustic,
            cancelled,
        )?;
        let physical_fields = spatial
            .iter()
            .copied()
            .chain(relaxing.iter().map(|law| &law.instantaneous))
            .collect::<Vec<_>>();
        let instantaneous = self.resistive_matrix(
            uniform,
            &physical_fields,
            &condensed.pressure_release,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(PreparedRelaxingTime {
            condensed,
            instantaneous,
            arms,
        })
    }
}

impl<'m, 'f, 'c> PreparedRelaxingTime<'m, 'f, 'c> {
    pub fn duration_s(&self) -> f64 {
        self.condensed.duration_s
    }
    pub fn relaxation_terms(&self) -> usize {
        self.arms.len()
    }
    /// Flattened arm order is patch order, then original authored term order.
    pub fn memory_patch_and_term(&self, arm: usize) -> Option<(usize, usize)> {
        self.arms.get(arm).map(|a| (a.patch, a.term))
    }
    /// Original geometric trace coordinates, including release-shared DOFs.
    pub fn memory_dofs(&self, arm: usize) -> Option<&[usize]> {
        self.arms.get(arm).map(|a| a.dofs.as_slice())
    }
    pub fn memory_relaxation_time_s(&self, arm: usize) -> Option<f64> {
        self.arms.get(arm).map(|a| a.tau)
    }

    /// Import original fluid fields with initially unexcited boundary memories.
    /// Nonzero pressure does not silently initialize a relaxed boundary state.
    pub fn state<'p>(
        &'p self,
        potential_pa_s: &[f64],
        pressure_pa: &[f64],
    ) -> Result<RelaxingTimeState<'p, 'm, 'f, 'c>, HeterogeneousError> {
        let state = RelaxingTimeState {
            owner: self,
            acoustic: self.condensed.state(potential_pa_s, pressure_pa)?,
            memories: self.arms.iter().map(|a| vec![0.0; a.dofs.len()]).collect(),
        };
        self.energy(&state)?;
        Ok(state)
    }

    /// Import every original memory explicitly in `memory_dofs` order.
    /// Independent memories at release-shared coordinates are retained.
    pub fn state_with_memories<'p>(
        &'p self,
        potential_pa_s: &[f64],
        pressure_pa: &[f64],
        memories: &[&[f64]],
    ) -> Result<RelaxingTimeState<'p, 'm, 'f, 'c>, HeterogeneousError> {
        if memories.len() != self.arms.len()
            || memories
                .iter()
                .zip(&self.arms)
                .any(|(q, a)| q.len() != a.dofs.len() || q.iter().any(|x| !x.is_finite()))
        {
            return Err(invalid("finite complete original admittance memories required").into());
        }
        let state = RelaxingTimeState {
            owner: self,
            acoustic: self.condensed.state(potential_pa_s, pressure_pa)?,
            memories: memories.iter().map(|q| q.to_vec()).collect(),
        };
        self.energy(&state)?;
        Ok(state)
    }

    fn validate_state(
        &self,
        state: &RelaxingTimeState<'_, '_, '_, '_>,
    ) -> Result<(), HeterogeneousError> {
        if !core::ptr::eq(self, state.owner) {
            return Err(invalid(
                "admittance memory state belongs to another prepared boundary law",
            )
            .into());
        }
        self.condensed.validate_state(&state.acoustic)
    }

    pub fn energy(
        &self,
        state: &RelaxingTimeState<'_, '_, '_, '_>,
    ) -> Result<RelaxingTimeEnergy, HeterogeneousError> {
        self.validate_state(state)?;
        self.energy_with_roundoff(state, &mut || false).map(|v| v.0)
    }

    fn energy_with_roundoff(
        &self,
        state: &RelaxingTimeState<'_, '_, '_, '_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(RelaxingTimeEnergy, f64, f64), HeterogeneousError> {
        let (acoustic, fluid_roundoff) = self
            .condensed
            .energy_with_roundoff(&state.acoustic, cancelled)?;
        let mut memory = Dd::from_f64(0.0);
        let mut memory_roundoff = 0.0;
        for (arm, q) in self.arms.iter().zip(&state.memories) {
            let (stored, roundoff) = positive_form(&arm.conductance, q, cancelled)?;
            memory = memory + Dd::from_f64(0.5) * Dd::from_f64(arm.tau) * Dd::from_f64(stored);
            memory_roundoff += 0.5 * (arm.tau * roundoff);
        }
        let energy = RelaxingTimeEnergy {
            acoustic,
            boundary_memory_j: memory.to_f64(),
        };
        if !energy.total_j().is_finite() || !memory_roundoff.is_finite() {
            return Err(invalid("spatial admittance stored energy overflow").into());
        }
        Ok((energy, fluid_roundoff, memory_roundoff))
    }

    pub fn step<'p>(
        &'p self,
        initial: &RelaxingTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RelaxingTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_boundary_sources(initial, inward_velocity, &[], &[], cancelled)
    }

    pub fn step_with_volume_sources<'p>(
        &'p self,
        initial: &RelaxingTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RelaxingTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        self.step_with_boundary_sources(initial, inward_velocity, &[], sources, cancelled)
    }

    /// Advance all original fluid and wall-memory coordinates with midpoint
    /// linear/quadratic boundary velocities and volumetric source amplitudes.
    /// Accepted inputs remain immutable on cancellation or numerical refusal.
    pub fn step_with_boundary_sources<'p>(
        &'p self,
        initial: &RelaxingTimeState<'_, '_, '_, '_>,
        inward_velocity: &[RealVelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, f64)],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RelaxingTimeStep<'p, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        self.validate_state(initial)?;
        let core = &self.condensed;
        let source =
            core.projected_source_load(inward_velocity, boundary_sources, sources, cancelled)?;
        let mut load = source.clone();
        for (arm, q) in self.arms.iter().zip(&initial.memories) {
            poll(cancelled)?;
            let mut projected = vec![0.0; arm.dofs.len()];
            arm.conductance.spmv(q, &mut projected);
            for (&dof, value) in arm.dofs.iter().zip(projected) {
                load[dof] = (Dd::from_f64(load[dof])
                    - Dd::from_f64(arm.history_weight) * Dd::from_f64(value))
                .to_f64();
            }
        }
        let acoustic = core.step_projected(&initial.acoustic, &load, cancelled)?;
        let midpoint_p = initial
            .pressure_pa()
            .iter()
            .zip(acoustic.state.pressure_pa())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let midpoint_psi = initial
            .pressure_potential_pa_s()
            .iter()
            .zip(acoustic.state.pressure_potential_pa_s())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        let mut next = RelaxingTimeState {
            owner: self,
            acoustic: acoustic.state,
            memories: Vec::with_capacity(self.arms.len()),
        };
        let h = self.duration_s();
        let mut outward = vec![0.0; source.len()];
        let mut memory_change = Dd::from_f64(0.0);
        let mut memory_loss = Dd::from_f64(0.0);
        let mut memory_work = Dd::from_f64(0.0);
        let mut loss_roundoff = 0.0;
        let mut work_roundoff = 0.0;
        let mut backward_error: f64 = 0.0;
        for (arm, initial_q) in self.arms.iter().zip(&initial.memories) {
            poll(cancelled)?;
            let mut endpoint = Vec::with_capacity(arm.dofs.len());
            let mut middle = Vec::with_capacity(arm.dofs.len());
            let mut pressure = Vec::with_capacity(arm.dofs.len());
            for (&dof, &q0) in arm.dofs.iter().zip(initial_q) {
                let p = midpoint_p[dof];
                let qm = (Dd::from_f64(arm.history_weight) * Dd::from_f64(q0)
                    + Dd::from_f64(arm.pressure_weight) * Dd::from_f64(p))
                .to_f64();
                let q1 = (Dd::from_f64(2.0) * Dd::from_f64(qm) - Dd::from_f64(q0)).to_f64();
                let defect = (Dd::from_f64(arm.tau) * (Dd::from_f64(q1) - Dd::from_f64(q0))
                    + Dd::from_f64(h) * (Dd::from_f64(qm) - Dd::from_f64(p)))
                .to_f64();
                let scale = arm.tau * (q0.abs() + q1.abs()) + h * (qm.abs() + p.abs());
                if ![qm, q1, defect, scale].into_iter().all(f64::is_finite) {
                    return Err(invalid("spatial admittance memory candidate overflow").into());
                }
                backward_error = backward_error.max(relative(defect.abs(), scale));
                endpoint.push(q1);
                middle.push(qm);
                pressure.push(p);
            }
            let mut projected = vec![0.0; arm.dofs.len()];
            arm.conductance.spmv(&middle, &mut projected);
            for (&dof, &value) in arm.dofs.iter().zip(&projected) {
                outward[dof] = (Dd::from_f64(outward[dof]) + Dd::from_f64(value)).to_f64();
            }
            let (work, roundoff) = paired_work(h, &pressure, &projected);
            memory_work = memory_work + Dd::from_f64(work);
            work_roundoff += roundoff;
            let (loss, roundoff) = positive_form(&arm.conductance, &middle, cancelled)?;
            memory_loss = memory_loss + Dd::from_f64(h) * Dd::from_f64(loss);
            loss_roundoff += h * roundoff;
            memory_change = memory_change
                + Dd::from_f64(arm.tau)
                    * Dd::from_f64(energy_change(
                        &arm.conductance,
                        initial_q,
                        &endpoint,
                        cancelled,
                    )?);
            next.memories.push(endpoint);
        }
        let (_, fluid_error) = core.physical_residual(
            &initial.acoustic,
            &next.acoustic,
            &midpoint_p,
            &midpoint_psi,
            &self.instantaneous,
            &source,
            Some(&outward),
            cancelled,
        )?;
        backward_error = backward_error.max(fluid_error);
        if backward_error > core.limits.residual_tolerance {
            return Err(HeterogeneousError::Residual {
                observed: backward_error,
                limit: core.limits.residual_tolerance,
            });
        }
        let (initial_energy, initial_fluid_error, initial_memory_error) =
            self.energy_with_roundoff(initial, cancelled)?;
        let (final_energy, final_fluid_error, final_memory_error) =
            self.energy_with_roundoff(&next, cancelled)?;
        let (instant_loss, instant_roundoff) =
            positive_form(&self.instantaneous, &midpoint_p, cancelled)?;
        let (input, input_roundoff) = paired_work(h, &midpoint_p, &source);
        let mut boundary_work = Dd::from_f64(0.0);
        for &(boundary, amplitude) in boundary_sources {
            poll(cancelled)?;
            let p = boundary.midpoint_pressure(initial.pressure_pa(), next.pressure_pa());
            boundary_work =
                boundary_work + Dd::from_f64(h) * Dd::from_f64(amplitude) * Dd::from_f64(p);
        }
        let mut volume_work = Dd::from_f64(0.0);
        for &(volume, amplitude) in sources {
            poll(cancelled)?;
            let p = volume.midpoint_pressure(initial.pressure_pa(), next.pressure_pa());
            volume_work = volume_work + Dd::from_f64(h) * Dd::from_f64(amplitude) * Dd::from_f64(p);
        }
        let work = RelaxingTimeWork {
            initial_energy_j: initial_energy.total_j(),
            final_energy_j: final_energy.total_j(),
            delta_energy_j: (Dd::from_f64(acoustic.work.delta_energy_j) + memory_change).to_f64(),
            acoustic_delta_energy_j: acoustic.work.delta_energy_j,
            memory_delta_energy_j: memory_change.to_f64(),
            input_work_j: input,
            boundary_source_work_j: boundary_work.to_f64(),
            volume_source_work_j: volume_work.to_f64(),
            instantaneous_dissipation_j: h * instant_loss,
            relaxation_dissipation_j: memory_loss.to_f64(),
            memory_work_j: memory_work.to_f64(),
            backward_error,
        };
        if !work.boundary_source_work_j.is_finite() || !work.volume_source_work_j.is_finite() {
            return Err(invalid("spatial admittance source work overflow").into());
        }
        let fluid_roundoff =
            initial_fluid_error + final_fluid_error + h * instant_roundoff + input_roundoff;
        let memory_roundoff = initial_memory_error + final_memory_error + loss_roundoff;
        let tolerance = core.limits.energy_balance_tolerance;
        balance(
            work.acoustic_balance_residual_j(),
            work.acoustic_delta_energy_j.abs()
                + work.instantaneous_dissipation_j
                + work.memory_work_j.abs()
                + input.abs(),
            fluid_roundoff + work_roundoff,
            tolerance,
        )?;
        balance(
            work.memory_balance_residual_j(),
            work.memory_delta_energy_j.abs()
                + work.relaxation_dissipation_j
                + work.memory_work_j.abs(),
            memory_roundoff + work_roundoff,
            tolerance,
        )?;
        balance(
            work.balance_residual_j(),
            work.delta_energy_j.abs()
                + work.instantaneous_dissipation_j
                + work.relaxation_dissipation_j
                + input.abs(),
            fluid_roundoff + memory_roundoff,
            tolerance,
        )?;
        poll(cancelled)?;
        Ok(RelaxingTimeStep { state: next, work })
    }
}

#[cfg(test)]
mod tests;
