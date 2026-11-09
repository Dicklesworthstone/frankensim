//! Original physical parameter gradients of complete autonomous pHS objectives.

use super::{
    AutonomousPortFeedback, PhsAutonomousParameterDirection, PreparedPhsAutonomousSensitivity,
    Result, with_parameter_direction,
};
use fs_orbit::sensitivity::{
    AutonomousObjectiveGradient, ForcedSensitivityError, PreparedAutonomousAdjoint,
};
use fs_phs::PortHamiltonian;

/// Source-bound scalar objective derivative through every original harmonic,
/// material memory and the solved physical frequency. The original reference
/// phase remains fixed; no phase equation is discarded or regularized.
pub struct PreparedPhsAutonomousAdjoint<'a, 's> {
    system: &'s PortHamiltonian,
    feedback: &'s dyn AutonomousPortFeedback,
    prepared: PreparedAutonomousAdjoint<'a>,
}

impl<'s> PreparedPhsAutonomousSensitivity<'s> {
    /// Prepare one original objective's bordered transpose solve. The covector
    /// includes every ordinary real Fourier coordinate, and the separate
    /// frequency partial holds coefficients and phase fixed. A mean_objective
    /// result supplies both directly. The solved frequency chain must not be
    /// counted again in a parameter's explicit objective partial.
    pub fn prepare_adjoint(
        &mut self,
        objective_gradient: &[f64],
        objective_frequency_derivative: f64,
        objective_scale: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsAutonomousAdjoint<'_, 's>> {
        Ok(PreparedPhsAutonomousAdjoint {
            system: self.system,
            feedback: self.feedback,
            prepared: self.prepared.prepare_adjoint(
                objective_gradient,
                objective_frequency_derivative,
                objective_scale,
                cancelled,
            )?,
        })
    }
}

impl PreparedPhsAutonomousAdjoint<'_, '_> {
    /// Shared complete dual diagnostics, phase weight and original residual
    /// direction entry. All weights retain physical objective/equation units.
    pub fn prepared(&self) -> &PreparedAutonomousAdjoint<'_> {
        &self.prepared
    }

    /// Differentiate an objective along one original J/R/G/storage/feedback
    /// direction. The same exact source assembly as solve_parameter supplies
    /// `(J-R)h_p + (J_p-R_p)e + G_p u + G u_p`, including the full flow chain
    /// in authored feedback. No state-direction solve or material refit occurs.
    ///
    /// explicit_derivative is q_p at fixed Fourier coefficients AND frequency.
    /// The bordered adjoint includes the unknown frequency change exactly once.
    /// Source refusal/cancellation publishes no partial result and leaves the
    /// original factor and dual reusable.
    pub fn gradient_parameter(
        &mut self,
        parameter_scale: f64,
        explicit_derivative: f64,
        direction: PhsAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        with_parameter_direction(self.system, self.feedback, direction, |action| {
            self.prepared.gradient_island_direction(
                parameter_scale,
                explicit_derivative,
                action,
                &mut cancelled,
            )
        })
    }
}
