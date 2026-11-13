//! Complete original physical power and per-cycle work gradients.

use super::{
    ForcedMeanObjectiveControls, ForcedSensitivityError, PhsForcedSensitivityError,
    PhsPowerObjective, PhsStorageHessianAction, PreparedPhsForcedSensitivity,
};
use crate::phs_orbit::sensitivity::parameter::{ode_clock_column, with_parameter_samples};
use crate::phs_orbit::{PeriodicPortInput, PhsForcedParameterDirection};
use fs_orbit::sensitivity::{
    ForcedObjectiveGradient, ForcedSensitivityDirection, PreparedForcedAdjoint,
};
use fs_phs::PortHamiltonian;

/// One original physical power/work objective and its reusable forced adjoint.
/// The original drive, source, complete state covector and fixed-phase clock
/// chart are bound together; queries supply only the original physical
/// parameter direction. No manually assembled direct objective partial is needed.
pub struct PreparedPhsForcedPowerAdjoint<'a, 's> {
    system: &'s PortHamiltonian,
    input: &'a PeriodicPortInput,
    harmonics: usize,
    objective: PhsPowerObjective,
    value: f64,
    value_sampling_defect: f64,
    gradient_sampling_defect: f64,
    prepared: PreparedForcedAdjoint<'a>,
}

impl<'s> PreparedPhsForcedSensitivity<'s> {
    /// Sample one original physical power or per-cycle work objective and
    /// prepare its transpose solve using the existing complete coupled factor.
    /// Objective units/tolerances apply to the selected power or energy, not to
    /// an unscaled intermediate mean. All state and material memories remain.
    pub fn prepare_power_adjoint(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsForcedPowerAdjoint<'_, 's>, PhsForcedSensitivityError> {
        self.prepare_physical_power_adjoint(objective, controls, None, cancelled)
    }

    /// The same objective with an explicit COMPLETE original storage Hessian.
    /// This action is used only to prepare the objective's state covector;
    /// subsequent physical parameter queries reuse its cached adjoint.
    pub fn prepare_power_adjoint_analytic(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsForcedPowerAdjoint<'_, 's>, PhsForcedSensitivityError> {
        self.prepare_physical_power_adjoint(objective, controls, Some(hessian), cancelled)
    }

    fn prepare_physical_power_adjoint(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsForcedPowerAdjoint<'_, 's>, PhsForcedSensitivityError> {
        let observation = self.mean_physical_power(objective, controls, hessian, &mut cancelled)?;
        let harmonics = self.prepared.harmonics();
        let prepared = self.prepared.prepare_adjoint(
            &observation.gradient,
            controls.objective_scale,
            &mut cancelled,
        )?;
        Ok(PreparedPhsForcedPowerAdjoint {
            system: self.system,
            input: &self.input,
            harmonics,
            objective,
            value: observation.value,
            value_sampling_defect: observation.value_sampling_defect,
            gradient_sampling_defect: observation.gradient_sampling_defect,
            prepared,
        })
    }
}

impl PreparedPhsForcedPowerAdjoint<'_, '_> {
    pub fn objective(&self) -> PhsPowerObjective {
        self.objective
    }

    /// Accepted fine-grid original objective, in power or energy units.
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Coarse/fine value difference normalized by the requested objective unit.
    pub fn value_sampling_defect(&self) -> f64 {
        self.value_sampling_defect
    }

    /// Complete normalized Fourier-covector sampling defect at preparation.
    pub fn gradient_sampling_defect(&self) -> f64 {
        self.gradient_sampling_defect
    }

    /// Original complete dual weights and independent dual-equation diagnostics.
    pub fn prepared(&self) -> &PreparedForcedAdjoint<'_> {
        &self.prepared
    }

    /// Complete gradient under simultaneous original J/R/G, storage, Fourier
    /// input and prescribed-clock changes. The pHS owner supplies both the
    /// direct power partial and physical residual column. Each cached grid
    /// contracts its OWN direct mean with the unchanged original adjoint.
    ///
    /// For per-cycle work, the explicit partial is
    /// `T * (power_p - Omega_p/Omega * power)` at fixed phase; the residual
    /// includes `i*n*Omega_p*X_n` exactly once. Mean power has no extra drive
    /// clock term at fixed Fourier coefficients. No original orbit, factor or
    /// objective covector is recomputed, and source refusal leaves this reusable.
    pub fn gradient_parameter(
        &mut self,
        chart: ForcedSensitivityDirection,
        direction: PhsForcedParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient, PhsForcedSensitivityError> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let n = self.system.state_dim();
        let objective = self.objective;
        with_parameter_samples(
            self.system,
            self.input,
            self.harmonics,
            direction,
            |samples| {
                self.prepared.gradient_mean_split_direction(
                    chart,
                    |time, state, out| {
                        samples.borrow_mut().objective_partial(
                            time,
                            state,
                            objective,
                            chart.angular_frequency_derivative,
                            out,
                        )
                    },
                    |time, state, out| samples.borrow_mut().vector_field(time, state, out),
                    |point, out| {
                        ode_clock_column(n, chart.angular_frequency_derivative, point, out)
                    },
                    &mut cancelled,
                )
            },
        )
    }
}
