//! The existing DG reaction iteration around the original simultaneous solve.
use super::*;
use super::super::nonlinear::{NonlinearTrial, NonlinearTrialWork};

struct MultipleTrial<'a> {
    inner: &'a mut MultiFrictionModalSystem,
    drives: &'a [FrictionDrive],
    frame: Option<MultiFrictionFrame>,
}
impl NonlinearTrial for MultipleTrial<'_> {
    fn core(&self) -> &FrictionModalSystem { &self.inner.core }
    fn stage(&mut self, forces: &[f64], gate: Option<&CancelGate>)
        -> Result<NonlinearTrialWork, ModalFrictionError>
    {
        let frame = self.inner.stage_multiple(forces, self.drives, None, gate)?;
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
        };
        self.frame = Some(frame);
        Ok(work)
    }
}
impl MultiFrictionModalSystem {
    pub(super) fn stage_with_nonlinear(&mut self, external: &[f64],
        drives: &[FrictionDrive], gate: Option<&CancelGate>)
        -> Result<MultiFrictionFrame, ModalFrictionError>
    {
        // Those compositions need a joint normal-contact/constitutive residual.
        // Admission already rejects them; never turn a future bypass into a
        // dropped stop, normal work port or material potential.
        if self.obstacles.is_some() || self.contact_loads.iter().any(Option::is_some) {
            return Err(invalid("nonlinear storage currently requires prescribed bow loads without obstacles").into());
        }
        let state = self.core.nonlinear.take()
            .ok_or_else(|| invalid("nonlinear trial requires retained constitutive storage"))?;
        let result = (|| {
            let mut trial = MultipleTrial { inner: &mut *self, drives, frame: None };
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
