//! Complete supply and resistive-loss objectives of the original feedback loop.

use super::{
    AutonomousMeanObjective, ForcedMeanObjectiveControls, ForcedSensitivityError,
    PhsFeedbackStateJacobian, PhsStorageHessianAction, PreparedPhsAutonomousSensitivity, Result,
};
use crate::phs_orbit::sensitivity::power::{PhsPowerObjective, PowerObservation};
use fs_math::dd::Dd;

mod adjoint;
pub use adjoint::PreparedPhsAutonomousPowerAdjoint;

impl PreparedPhsAutonomousSensitivity<'_> {
    /// Mean original `grad(H)^T R grad(H)` and complete coefficient covector.
    /// Every nonlinear/material-memory state participates through the exact
    /// original Hamiltonian Hessian. Its explicit fixed-state clock partial is
    /// zero; the bordered adjoint still includes the solved frequency chain.
    /// This observes signed original R and does not certify a heat destination.
    pub fn mean_resistance_power(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.mean_physical_power(
            PhsPowerObjective::MeanResistancePower,
            controls,
            None,
            None,
            cancelled,
        )
    }

    /// Original resistance power with an explicit complete Hamiltonian Hessian.
    pub fn mean_resistance_power_analytic(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.mean_physical_power(
            PhsPowerObjective::MeanResistancePower,
            controls,
            Some(hessian),
            None,
            cancelled,
        )
    }

    /// Mean supplied power of the ORIGINAL instantaneous physical feedback.
    /// Uses its actual input and original conjugate flows. The supplied complete
    /// feedback Jacobian must be the same total state derivative as in orbit
    /// preparation, including its flow/Hessian chain exactly once.
    ///
    /// The complete gradient is `Hxx G u + (du_total/dx)^T y`. It includes
    /// explicit controller/constitutive memories; treating feedback as an
    /// externally prescribed drive would omit its second term. Fixed-state
    /// clock dependence is zero for this autonomous source. Direct parameter
    /// changes of H/G/feedback remain separate explicit adjoint partials.
    pub fn mean_supplied_power(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.mean_physical_power(
            PhsPowerObjective::MeanSuppliedPower,
            controls,
            None,
            Some(feedback_jacobian),
            cancelled,
        )
    }

    /// Same original supply observation with an explicit complete Hessian.
    pub fn mean_supplied_power_analytic(
        &mut self,
        controls: ForcedMeanObjectiveControls,
        feedback_jacobian: &PhsFeedbackStateJacobian<'_>,
        hessian: &PhsStorageHessianAction<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        self.mean_physical_power(
            PhsPowerObjective::MeanSuppliedPower,
            controls,
            Some(hessian),
            Some(feedback_jacobian),
            cancelled,
        )
    }

    fn mean_physical_power(
        &mut self,
        objective: PhsPowerObjective,
        controls: ForcedMeanObjectiveControls,
        hessian: Option<&PhsStorageHessianAction<'_>>,
        feedback_jacobian: Option<&PhsFeedbackStateJacobian<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (n, m) = (self.system.state_dim(), self.system.port_dim());
        let mut observation = PowerObservation::new(self.system, hessian)?;
        let mut flows = vec![f64::NAN; m];
        let mut input = vec![f64::NAN; m];
        let mut jacobian = if feedback_jacobian.is_some() {
            vec![f64::NAN; m * n]
        } else {
            Vec::new()
        };
        let feedback = self.feedback;
        let omega = self.prepared.angular_frequency();
        self.prepared
            .mean_objective(
                controls,
                |_time, state, value, gradient, frequency| {
                    *frequency = 0.0;
                    if !observation.load(state) {
                        return false;
                    }
                    if !objective.supplied() {
                        return observation.resistance(state, value, gradient)
                            && objective.normalize(omega, value, gradient, Some(frequency));
                    }
                    let Some(derivative) = feedback_jacobian else {
                        return false;
                    };
                    input.fill(f64::NAN);
                    if !observation.flows(&mut flows)
                        || !feedback.input(state, &flows, &mut input)
                        || input.iter().any(|v| !v.is_finite())
                        || !observation.supplied(state, &input, value, gradient)
                    {
                        return false;
                    }
                    jacobian.fill(f64::NAN);
                    if !derivative(state, &flows, &mut jacobian)
                        || jacobian.iter().any(|v| !v.is_finite())
                    {
                        return false;
                    }
                    for (column, value) in gradient.iter_mut().enumerate() {
                        let mut sum = Dd::from_f64(*value);
                        for (port, &flow) in flows.iter().enumerate() {
                            let derivative = jacobian[port * n + column];
                            if derivative != 0.0 && flow != 0.0 {
                                sum = sum + Dd::from_f64(derivative) * Dd::from_f64(flow);
                            }
                        }
                        *value = sum.to_f64();
                    }
                    objective.normalize(omega, value, gradient, Some(frequency))
                },
                cancelled,
            )
            .map_err(Into::into)
    }
}
