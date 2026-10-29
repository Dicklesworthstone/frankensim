//! Original native mechanics continues after cohesive failure with sourced
//! residual friction. Enabling the handoff changes no body, clock or history.
//! Actual initial and endpoint slip velocities preserve kinetic occupancy even
//! when the audited path has zero net slip. Intrastep stopping and reversal are
//! not event-resolved by this held-force transaction.

use super::*;
use crate::cohesive_contact::reclosure::{
    ReclosureActivation, ReclosureConfig, ReclosureError, ReclosureLoadLine, ReclosureRegime,
    ReclosureState, ResidualCoulombSource,
};

pub(super) struct ReclosureBinding<'a> {
    load_line: ReclosureLoadLine<'a>,
    state: ReclosureState<'a>,
}

impl<'a> CohesiveModalSystem<'a> {
    /// Explicitly select residual contact on the SAME fully failed interface.
    /// Both material assignments, original interface card and fixed query must
    /// agree. Initial residual occupancy is explicit; a failed material point
    /// contains no velocity or prescribed-motion history from which to infer it.
    /// Later opening clears occupancy, and reclosure enters a fresh Coulomb
    /// selection. No fracture, body state, sample or observation is reset.
    pub fn enable_frictional_reclosure(
        &mut self,
        friction: &'a ResidualCoulombSource,
        receiving_materials: [&MaterialStateId; 2],
        activation: ReclosureActivation,
        mut config: ReclosureConfig,
        gate: &CancelGate,
    ) -> Result<(), CohesiveModalError> {
        poll(Some(gate))?;
        if self.residual.is_some() {
            return Err(invalid("an existing residual interface cannot be replaced").into());
        }
        if !self.point.failed() {
            return Err(ReclosureError::NotFullyFailed.into());
        }
        let source = self
            .source
            .ok_or_else(|| invalid("residual handoff requires the original cohesive source"))?;
        if config.max_evaluations < 3 {
            return Err(
                invalid("residual native budget must reserve its actual-motion audit").into(),
            );
        }
        config.max_evaluations -= 1;
        let load_line = ReclosureLoadLine::new(
            source,
            friction,
            self.area_m2(),
            self.compliance_m_per_n(),
            config,
        )?;
        let state = load_line.adopt(self.point, receiving_materials, activation)?;
        poll(Some(gate))?;
        self.residual = Some(ReclosureBinding { load_line, state });
        Ok(())
    }

    #[must_use]
    pub fn residual_friction_source(&self) -> Option<&'a ResidualCoulombSource> {
        self.residual
            .as_ref()
            .map(|binding| binding.load_line.friction_source())
    }
    #[must_use]
    pub fn residual_contact_occupancy(&self) -> Option<ReclosureRegime> {
        self.residual
            .as_ref()
            .map(|binding| binding.state.occupancy())
    }

    pub(super) fn step_reclosure(
        &mut self,
        external: &[f64],
        prescribed_rate: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<&CohesiveModalFrame, CohesiveModalError> {
        poll(gate)?;
        let before = self.total_energy_j()?;
        let old = self.point.jump_m();
        let area = self.area_m2();
        let dt = self.network.dt;
        let binding = self
            .residual
            .as_ref()
            .ok_or_else(|| invalid("residual handoff is absent"))?;
        if binding.state.point() != self.point {
            return Err(ReclosureError::StateMismatch.into());
        }
        let mut prescribed_delta = [0.0; 3];
        for i in 0..3 {
            prescribed_delta[i] = finite(prescribed_rate[i] * dt)?;
        }
        let previous = binding
            .state
            .with_endpoint_tangential_velocity(tangential_velocity(
                &self.columns,
                self.network.components(),
                prescribed_rate,
            )?)?;
        self.network.prepare_prediction(external, gate)?;
        let mut free = [0.0; 3];
        for i in 0..3 {
            free[i] = finite(
                old[i] + prescribed_delta[i] + dot(&self.columns[i], &self.network.free_delta)?,
            )?;
        }
        let solved = binding.load_line.trial(previous, free, dt, || {
            !gate.is_some_and(CancelGate::is_requested)
        })?;
        self.forces.copy_from_slice(external);
        for i in 0..3 {
            limit(
                "original interface force budget",
                solved.mean_force_n[i].abs(),
                self.config.contact.maximum_force_n,
            )?;
            for (force, weight) in self.forces.iter_mut().zip(&self.columns[i]) {
                *force = finite(*force - weight * solved.mean_force_n[i])?;
            }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let mut actual_jump = [0.0; 3];
        let mut compatibility = [0.0; 3];
        let solved_jump = solved.audit.state.point().jump_m();
        for i in 0..3 {
            actual_jump[i] = finite(
                old[i] + prescribed_delta[i] + dot(&self.columns[i], &self.network.free_delta)?,
            )?;
            compatibility[i] = finite(actual_jump[i] - solved_jump[i])?;
            if compatibility[i].abs() > solved.compatibility_tolerance_m {
                return Err(CohesiveModalError::JumpMismatch {
                    component: i,
                    residual_m: compatibility[i],
                    tolerance_m: solved.compatibility_tolerance_m,
                });
            }
            limit(
                "original cohesive jump budget",
                (actual_jump[i] - old[i]).abs(),
                self.config.contact.maximum_jump_increment_m,
            )?;
        }
        poll(gate)?;
        let actual = binding.load_line.audit(
            previous,
            actual_jump,
            dt,
            solved.mean_force_n,
            solved.audit.regime,
        )?;
        // Held-force motion may reverse velocity while its net slip is zero.
        // Preserve that actual endpoint motion in the next-step occupancy;
        // the frame still reports the separately audited path regime.
        let next_state = actual
            .state
            .with_endpoint_tangential_velocity(tangential_velocity(
                &self.columns,
                self.network.staged_components(),
                prescribed_rate,
            )?)?;
        let mut force_tolerance = [0.0; 3];
        for i in 0..3 {
            limit(
                "original constitutive force budget",
                actual.expected_force_n[i].abs(),
                self.config.contact.maximum_force_n,
            )?;
            let native_tolerance = finite(
                self.config.force_absolute_tolerance_n
                    + self.config.force_relative_tolerance
                        * actual.expected_force_n[i]
                            .abs()
                            .max(solved.mean_force_n[i].abs()),
            )?;
            force_tolerance[i] = native_tolerance.min(actual.force_tolerance_n[i]);
            if actual.force_residual_n[i].abs() > force_tolerance[i] {
                return Err(CohesiveModalError::ForceMismatch {
                    component: i,
                    residual_n: actual.force_residual_n[i],
                    tolerance_n: force_tolerance[i],
                });
            }
        }
        let frame = &self.network.staged_frame;
        let network_energy = finite(frame.total_energy_j())?;
        let cohesive_energy = finite(area * actual.state.point().stored_energy_j_m2())?;
        let network_loss = finite(frame.total_dissipation_j())?;
        let fracture = finite(area * actual.material.equilibrium.energy.dissipated)?;
        let viscous = finite(area * actual.material.viscous_dissipated_j_m2)?;
        let friction = actual.friction_dissipation_j;
        let external_work = dot(external, &self.network.free_delta)?;
        let prescribed_work = dot(&solved.mean_force_n, &prescribed_delta)?;
        let after = finite(network_energy + cohesive_energy)?;
        limit(
            "total energy including failed interface",
            after,
            self.network.config.maximum_total_energy_j,
        )?;
        let loss = finite(network_loss + fracture + viscous + friction)?;
        let residual = finite((after - before) + loss - external_work - prescribed_work)?;
        let scale = before
            .max(after)
            .max(external_work.abs())
            .max(prescribed_work.abs())
            .max(loss.abs());
        let tolerance = finite(
            self.network.config.energy_absolute_tolerance_j
                + self.network.config.energy_relative_tolerance * scale,
        )?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            }
            .into());
        }
        let accepted = CohesiveModalFrame {
            sample: frame.sample,
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: network_energy,
            cohesive_energy_j: cohesive_energy,
            external_work_j: external_work,
            prescribed_work_j: prescribed_work,
            network_dissipation_j: network_loss,
            fracture_dissipation_j: fracture,
            viscous_dissipation_j: viscous,
            residual_friction_dissipation_j: friction,
            residual_friction_regime: Some(actual.regime),
            residual_friction_work_defect_j: actual.work_residual_j,
            cohesive_work_j: finite(area * actual.material.total_work_j_m2())?,
            energy_residual_j: residual,
            energy_tolerance_j: tolerance,
            mean_force_n: solved.mean_force_n,
            endpoint_traction_pa: actual.material.traction_pa,
            jump_m: actual_jump,
            damage: actual.state.point().damage(),
            compatibility_residual_m: compatibility,
            compatibility_tolerance_m: solved.compatibility_tolerance_m,
            constitutive_residual_n: actual.force_residual_n,
            force_tolerance_n: force_tolerance,
            iterations: solved.iterations,
            material_evaluations: solved.evaluations + 1,
        };
        self.network.validate_staged_material_points(gate)?;
        poll(gate)?;
        self.network.publish_staged();
        self.point = next_state.point();
        self.residual
            .as_mut()
            .expect("retained residual binding")
            .state = next_state;
        Ok(self.last.insert(accepted))
    }
}

fn tangential_velocity(
    columns: &[Vec<f64>; 3],
    components: &[ModalAcousticTimeModel],
    prescribed_rate: [f64; 3],
) -> Result<[f64; 2], ModalCouplingError> {
    let mut result = [prescribed_rate[1], prescribed_rate[2]];
    for (index, mode) in components
        .iter()
        .flat_map(|model| model.states())
        .enumerate()
    {
        for tangent in 0..2 {
            result[tangent] = finite(
                result[tangent] + columns[tangent + 1][index] * mode.velocity_m_sqrt_kg_per_s,
            )?;
        }
    }
    Ok(result)
}
