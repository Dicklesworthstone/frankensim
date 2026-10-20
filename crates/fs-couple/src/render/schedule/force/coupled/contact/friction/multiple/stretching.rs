//! The existing DG reaction iteration around the original simultaneous solve.
use super::*;
use super::super::nonlinear::{NonlinearContactWork, NonlinearTrial, NonlinearTrialWork};

struct MultipleTrial<'a> {
    inner: &'a mut MultiFrictionModalSystem,
    drives: &'a [FrictionDrive],
    obstacle_offsets: Option<&'a [f64]>,
    initial_contact_energy_j: Option<f64>,
    frame: Option<MultiFrictionFrame>,
}
impl NonlinearTrial for MultipleTrial<'_> {
    fn core(&self) -> &FrictionModalSystem { &self.inner.core }
    fn stage(&mut self, forces: &[f64], gate: Option<&CancelGate>)
        -> Result<NonlinearTrialWork, ModalFrictionError>
    {
        let frame = self.inner.stage_multiple(forces, self.drives, self.obstacle_offsets, gate)?;
        let mut absolute_drive_work_j = 0.0;
        for port in &self.inner.staged_ports {
            absolute_drive_work_j = finite(absolute_drive_work_j + port.drive_work_j.abs())?;
        }
        let work = NonlinearTrialWork {
            observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: frame.network_energy_j,
            drive_work_j: frame.drive_work_j,
            absolute_drive_work_j,
            network_dissipation_j: frame.network_dissipation_j,
            friction_dissipation_j: frame.friction_dissipation_j,
            contacts: match (&self.inner.obstacles, self.initial_contact_energy_j) {
                (Some(obstacles), Some(before_energy_j)) => Some(NonlinearContactWork {
                    before_energy_j,
                    after_energy_j: frame.obstacle_energy_j,
                    dissipation_j: frame.obstacle_dissipation_j,
                    work_j: frame.obstacle_work_j,
                    absolute_work_j: obstacles.staged_absolute_work_j()?,
                }),
                (None, None) => None,
                _ => return Err(invalid("nonlinear trial lost its original contact storage").into()),
            },
        };
        self.frame = Some(frame);
        Ok(work)
    }
}
impl MultiFrictionModalSystem {
    pub(super) fn stage_with_nonlinear(&mut self, external: &[f64],
        drives: &[FrictionDrive], obstacle_offsets: Option<&[f64]>, gate: Option<&CancelGate>)
        -> Result<MultiFrictionFrame, ModalFrictionError>
    {
        // Accepted support offsets, never the pending commanded endpoints.
        let initial_contact_energy_j = self.obstacles.as_ref()
            .map(|obstacles| obstacles.energy(self.components())).transpose()?;
        let state = self.core.nonlinear.take()
            .ok_or_else(|| invalid("nonlinear trial requires retained constitutive storage"))?;
        let result = (|| {
            let mut trial = MultipleTrial { inner: &mut *self, drives,
                obstacle_offsets, initial_contact_energy_j, frame: None };
            let correction = state.stage_trial(&mut trial, external, gate)?;
            let mut frame = trial.frame
                .ok_or_else(|| invalid("nonlinear solve accepted no simultaneous mechanical trial"))?;
            frame.observer_pressure_pa = correction.observer_pressure_pa;
            frame.network_energy_j = correction.network_energy_j;
            frame.external_work_j = correction.external_work_j;
            frame.energy_residual_j = correction.energy_residual_j;
            frame.energy_tolerance_j = correction.energy_tolerance_j;
            frame.nonlinear = Some(correction.diagnostics);
            Ok(frame)
        })();
        // Restore even on domain, budget, convergence or cancellation failure.
        // Neither the accepted state nor its report was changed by the trials.
        self.core.nonlinear = Some(state);
        result
    }
}
