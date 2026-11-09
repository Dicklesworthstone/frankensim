//! Second-order Strang transport and energy-consistent implicit midpoint.
//!
//! At fixed transported areas, solve Y_mid-Y_start=(h/2)G(Y_mid) with the
//! original joint nonlinear owner, then reflect conserved momentum and energy.
//! The midpoint temperature is the original inverse of midpoint caloric energy.
//! The identity Delta K = h <u_mid,Gq_mid> closes against h*P_mid, with no
//! backward-Euler damping. Nonlinear/roundoff defects remain reported residuals.
//!
//! The dissipative method is A-stable, NOT L-stable. A very large step can reverse
//! an unresolved viscous mode or exchange hot/cold cells while conserving total
//! energy. Temporal accuracy and all reflected endpoint source checks remain
//! necessary; adaptive comparison is provided by the shared original driver.

use super::*;

/// The original input of a complete split; source preparation remains shared.
pub(super) struct Checkpoint {
    fields: Fields,
    energy: Vec<f64>,
}
impl Checkpoint {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        Ok(Self {
            fields: Fields::zeros(n, gate)?,
            energy: zeros(n, gate)?,
        })
    }
    fn save(
        &mut self,
        jet: &ThermalCarreauJet<'_>,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        for i in 0..self.energy.len() {
            poll_at(gate, i)?;
            self.fields.a[i] = jet.flow.candidate.a[i];
            self.fields.q[i] = jet.flow.candidate.q[i];
            self.energy[i] = jet.candidate[i];
        }
        Ok(())
    }
    fn restore(
        &self,
        jet: &mut ThermalCarreauJet<'_>,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        for i in 0..self.energy.len() {
            poll_at(gate, i)?;
            jet.flow.candidate.a[i] = self.fields.a[i];
            jet.flow.candidate.q[i] = self.fields.q[i];
            jet.candidate[i] = self.energy[i];
        }
        Ok(())
    }
}

impl ThermalCarreauJet<'_> {
    /// Candidate restoration keeps the original fixed and adaptive retry contract.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn try_midpoint_step(
        &mut self,
        dt: f64,
        base: Evaluation,
        cfl: f64,
        controls: ThermalJetImplicitControls,
        heat: TrialHeat,
        work: &mut Workspace,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Option<TrialStep>, ThermalJetImplicitError> {
        let mut checkpoint = match work.midpoint.take() {
            Some(checkpoint) => checkpoint,
            None => Checkpoint::new(self.accepted.len(), gate)?,
        };
        let result = (|| {
            checkpoint.save(self, gate)?;
            self.midpoint_step_inner(dt, base, cfl, controls, heat, work, evaluations, gate)
        })();
        let restore = matches!(
            &result,
            Ok(None)
                | Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
                    NonlinearJetError::InitialDomain(_)
                )))
        );
        let restored = if restore {
            checkpoint.restore(self, gate)
        } else {
            Ok(())
        };
        work.midpoint = Some(checkpoint);
        restored?;
        result
    }

    /// One SSPRK2 conservative transport step, leaving its endpoint in stage.
    /// No physical viscous heat is integrated during the transport subflow.
    fn transport_rk2(
        &mut self,
        dt: f64,
        base: Evaluation,
        cfl: f64,
        gate: &CancelGate,
    ) -> Result<Option<(Evaluation, f64)>, ThermalJetImplicitError> {
        if dt > stable_step(base.rate_per_s, cfl) {
            return Ok(None);
        }
        for i in 0..self.candidate.len() {
            poll_at(gate, i)?;
            self.flow.stage.a[i] = self.flow.candidate.a[i] + dt * self.flow.rhs.a[i];
            self.flow.stage.q[i] = self.flow.candidate.q[i] + dt * self.flow.rhs.q[i];
            self.stage[i] = finite(self.candidate[i] + dt * self.rhs[i])?;
        }
        let first = self.evaluate_transport(true, gate)?;
        if dt > stable_step(first.rate_per_s, cfl) {
            return Ok(None);
        }
        self.blend(dt, 0.5, gate)?;
        let end = self.evaluate_transport(true, gate)?;
        if dt > stable_step(end.rate_per_s, cfl) {
            return Ok(None);
        }
        let upwind = finite(0.5 * dt * base.upwind_power_w + 0.5 * dt * first.upwind_power_w)?;
        Ok(Some((end, upwind)))
    }

    #[allow(clippy::too_many_arguments)]
    fn midpoint_step_inner(
        &mut self,
        dt: f64,
        base: Evaluation,
        cfl: f64,
        controls: ThermalJetImplicitControls,
        heat: TrialHeat,
        work: &mut Workspace,
        evaluations: &mut usize,
        gate: &CancelGate,
    ) -> Result<Option<TrialStep>, ThermalJetImplicitError> {
        let half = positive(0.5 * dt)?;
        let Some((_, first_upwind)) = self.transport_rk2(half, base, cfl, gate)? else {
            return Ok(None);
        };
        core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
        core::mem::swap(&mut self.candidate, &mut self.stage);

        // A half-duration joint BE solve is exactly the implicit midpoint
        // equation. Its source remains the transported q/E, not a frozen T.
        let midpoint_controls = ThermalJetImplicitControls {
            velocity_tolerance_m_s: 0.5 * controls.velocity_tolerance_m_s,
            cell_energy_tolerance_j: 0.5 * controls.cell_energy_tolerance_j,
            ..controls
        };
        let mut solved = solve_stage(self, work, half, midpoint_controls, evaluations, gate)?;
        let physical_heat = finite(dt * solved.evaluation.viscous_power_w)?;
        let mut velocity_residual = 0.0_f64;
        let mut energy_residual = 0.0_f64;
        for i in 0..self.candidate.len() {
            poll_at(gate, i)?;
            // Centered differences preserve a large common velocity or caloric
            // reference. Check the ACTUAL reflected f64 state, not just 2*x.
            let dq = finite(self.flow.stage.q[i] - self.flow.candidate.q[i])?;
            let de = finite(self.stage[i] - self.candidate[i])?;
            self.flow.stage.q[i] = finite(self.flow.candidate.q[i] + 2.0 * dq)?;
            self.stage[i] = finite(self.candidate[i] + 2.0 * de)?;
            let momentum_defect = finite(
                ((self.flow.stage.q[i] - self.flow.candidate.q[i]) - dt * self.flow.rhs.q[i])
                    / work.cells[i].face_mass,
            )?;
            let caloric_defect = finite((self.stage[i] - self.candidate[i]) - dt * self.rhs[i])?;
            velocity_residual = velocity_residual.max(momentum_defect.abs());
            energy_residual = energy_residual.max(caloric_defect.abs());
        }
        if velocity_residual > controls.velocity_tolerance_m_s
            || energy_residual > controls.cell_energy_tolerance_j
        {
            return Err(ThermalJetImplicitError::ComponentResidual {
                velocity_m_s: velocity_residual,
                energy_j: energy_residual,
            });
        }
        solved.velocity_residual = solved.velocity_residual.max(velocity_residual);
        solved.energy_residual = solved.energy_residual.max(energy_residual);
        // Reflection does not change the retained spatial midpoint power.
        // Capture it BEFORE source admission overwrites that scratch.
        self.record_trial_heat(heat, dt, gate)?;
        let reflected = self.evaluate_transport(true, gate)?;
        core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
        core::mem::swap(&mut self.candidate, &mut self.stage);

        let Some((end, second_upwind)) = self.transport_rk2(half, reflected, cfl, gate)? else {
            return Ok(None);
        };
        solved.evaluation = end;
        Ok(Some(TrialStep {
            solved,
            heat: physical_heat,
            upwind: finite(first_upwind + second_upwind)?,
            damping: 0.0,
        }))
    }
}

#[cfg(test)]
mod tests;
