//! Complete mean-power and period-dependent work gradients of physical cycles.

use super::{
    ForcedMeanObjectiveControls, ForcedSensitivityError, PhsFeedbackStateJacobian,
    PhsPowerObjective, PhsStorageHessianAction, PreparedPhsAutonomousSensitivity, Result,
};
use crate::phs_orbit::autonomous::sensitivity::parameter::with_parameter_samples;
use crate::phs_orbit::autonomous::{AutonomousPortFeedback, PhsAutonomousParameterDirection};
use fs_orbit::sensitivity::{AutonomousObjectiveGradient, PreparedAutonomousAdjoint};
use fs_phs::PortHamiltonian;

/// One original physical power/work objective with its unchanged phase-bordered
/// adjoint. The full storage, feedback, memory and solved-frequency chain are
/// retained. Queries need no hand-assembled objective parameter partial.
pub struct PreparedPhsAutonomousPowerAdjoint<'a, 's> {
    system: &'s PortHamiltonian,
    feedback: &'s dyn AutonomousPortFeedback,
    omega: f64,
    objective: PhsPowerObjective,
    value: f64,
    angular_frequency_partial: f64,
    value_sampling_defect: f64,
    gradient_sampling_defect: f64,
    prepared: PreparedAutonomousAdjoint<'a>,
}

impl<'s> PreparedPhsAutonomousSensitivity<'s> {
    /// Sample original mean power or energy per ACTUAL cycle and reuse the
    /// original bordered factor for one objective transpose solve. Supply uses
    /// the same complete total feedback state Jacobian as orbit preparation;
    /// resistance power does not evaluate that argument.
    ///
    /// For cycle energy the sampled covector includes `q_Omega=-q/Omega` at
    /// fixed coefficients. Both its state and clock components are accepted in
    /// the selected energy units, preserving the solved period chain once.
    pub fn prepare_power_adjoint(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsAutonomousPowerAdjoint<'_, 's>> {
        self.prepare_physical_power_adjoint(objective, controls, feedback_jacobian, None, cancelled)
    }

    /// The same original objective with an explicitly authored COMPLETE storage
    /// Hessian. No quadratic approximation or numerical perturbation is used.
    pub fn prepare_power_adjoint_analytic(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsAutonomousPowerAdjoint<'_, 's>> {
        self.prepare_physical_power_adjoint(
            objective,
            controls,
            feedback_jacobian,
            Some(hessian),
            cancelled,
        )
    }

    fn prepare_physical_power_adjoint(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsAutonomousPowerAdjoint<'_, 's>> {
        let observation = self.mean_physical_power(
            objective,
            controls,
            hessian,
            objective.supplied().then_some(feedback_jacobian),
            &mut cancelled,
        )?;
        let omega = self.prepared.angular_frequency();
        let prepared = self.prepared.prepare_adjoint(
            &observation.gradient,
            observation.angular_frequency_partial,
            controls.objective_scale,
            &mut cancelled,
        )?;
        Ok(PreparedPhsAutonomousPowerAdjoint {
            system: self.system,
            feedback: self.feedback,
            omega,
            objective,
            value: observation.value,
            angular_frequency_partial: observation.angular_frequency_partial,
            value_sampling_defect: observation.value_sampling_defect,
            gradient_sampling_defect: observation.gradient_sampling_defect,
            prepared,
        })
    }
}

impl PreparedPhsAutonomousPowerAdjoint<'_, '_> {
    pub fn objective(&self) -> PhsPowerObjective {
        self.objective
    }

    /// Accepted original objective, in physical power or energy units.
    pub fn value(&self) -> f64 {
        self.value
    }

    /// Explicit physical-frequency partial retained in the bordered dual.
    pub fn angular_frequency_partial(&self) -> f64 {
        self.angular_frequency_partial
    }

    /// Coarse/fine value difference normalized by the requested objective unit.
    pub fn value_sampling_defect(&self) -> f64 {
        self.value_sampling_defect
    }

    /// Complete normalized coefficient AND clock covector sampling defect.
    pub fn gradient_sampling_defect(&self) -> f64 {
        self.gradient_sampling_defect
    }

    /// Original bordered weights, phase weight and dual-equation diagnostics.
    pub fn prepared(&self) -> &PreparedAutonomousAdjoint<'_> {
        &self.prepared
    }

    /// Differentiate simultaneous original J/R/G/storage/feedback changes.
    /// Complete `G_p^T e + G^T h_p` enters the authored feedback parameter law
    /// once, and the original pHS owner contracts both direct power and residual
    /// derivatives. The source is never replaced by a neighboring/fitted model.
    ///
    /// Cycle energy's explicit fixed-frequency partial is `T*power_p`; its
    /// unknown period change belongs to the already prepared bordered adjoint.
    /// Each grid's direct mean and residual column are contracted together and
    /// their TOTAL gradient must agree. No new state solve, factor or objective
    /// Hessian/feedback Jacobian evaluation occurs for a parameter query.
    pub fn gradient_parameter(
        &mut self,
        parameter_scale: f64,
        direction: PhsAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let objective = self.objective;
        let omega = self.omega;
        with_parameter_samples(self.system, self.feedback, direction, |samples| {
            self.prepared.gradient_mean_island_direction(
                parameter_scale,
                |_time, state, out| {
                    samples
                        .borrow_mut()
                        .objective_partial(state, objective, omega, out)
                },
                |_time, state, out| samples.borrow_mut().vector_field(state, out),
                &mut cancelled,
            )
        })
    }
}
