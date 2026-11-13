//! Independent physical supply, native resistance and nonlinear internal heat.

use super::{
    AutonomousMeanObjective, AutonomousObjectiveGradient, DissipativeAutonomousHamiltonian,
    DissipativeAutonomousParameterDirection, DissipativeError, ForcedMeanObjectiveControls,
    ForcedSensitivityError, PreparedAutonomousAdjoint, PreparedDissipativeAutonomousSensitivity,
    Result, parameter,
};
use crate::phs_orbit::sensitivity::power::PhsPowerObjective;
use fs_math::dd::Dd;
use fs_phs::harmonic::equilibrium::LinearizationError;

/// Separate original accounts, either cycle-mean power or work over the actual
/// solved period. No account is inferred by balancing the other two.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DissipativeAutonomousPowerObjective {
    MeanSuppliedPower,
    MeanResistancePower,
    MeanNonlinearPower,
    SuppliedEnergyPerCycle,
    ResistanceEnergyPerCycle,
    NonlinearEnergyPerCycle,
}
#[derive(Clone, Copy)]
pub(super) enum Kind {
    Supplied,
    Resistance,
    Nonlinear,
}
impl DissipativeAutonomousPowerObjective {
    pub(super) fn kind(self) -> Kind {
        match self {
            Self::MeanSuppliedPower | Self::SuppliedEnergyPerCycle => Kind::Supplied,
            Self::MeanResistancePower | Self::ResistanceEnergyPerCycle => Kind::Resistance,
            Self::MeanNonlinearPower | Self::NonlinearEnergyPerCycle => Kind::Nonlinear,
        }
    }
    // Unit normalization only. Original R and D observations remain independent.
    pub(super) fn normalization(self) -> PhsPowerObjective {
        match self {
            Self::MeanSuppliedPower => PhsPowerObjective::MeanSuppliedPower,
            Self::MeanResistancePower | Self::MeanNonlinearPower => {
                PhsPowerObjective::MeanResistancePower
            }
            Self::SuppliedEnergyPerCycle => PhsPowerObjective::SuppliedEnergyPerCycle,
            Self::ResistanceEnergyPerCycle | Self::NonlinearEnergyPerCycle => {
                PhsPowerObjective::ResistanceEnergyPerCycle
            }
        }
    }
}

impl DissipativeAutonomousHamiltonian<'_> {
    fn observe_power(
        &self,
        objective: DissipativeAutonomousPowerObjective,
        omega: f64,
        state: &[f64],
        value: &mut f64,
        gradient: &mut [f64],
        frequency: &mut f64,
    ) -> Result<()> {
        let n = self.system.state_dim();
        if gradient.len() != n {
            return Err(DissipativeError::SourceRefused.into());
        }
        let mut scratch = self.scratch.borrow_mut();
        let nonlinear_power = self.evaluate(state, &mut scratch)?;
        let super::super::Scratch {
            effort,
            resistance,
            state_direction,
            effort_direction,
            resistance_direction,
            field,
            flows,
            input,
            feedback_derivative,
            ..
        } = &mut *scratch;
        match objective.kind() {
            Kind::Supplied => {
                *value = self.system.supplied_power_from_effort(effort, input)?;
                self.system.input_flow_into(input, field)?;
                if !self.tangent.apply_hessian(state, field, gradient) {
                    return Err(DissipativeError::SourceRefused.into());
                }
                feedback_derivative.fill(f64::NAN);
                if !(self.feedback_jacobian)(state, flows, feedback_derivative)
                    || feedback_derivative.iter().any(|v| !v.is_finite())
                {
                    return Err(DissipativeError::SourceRefused.into());
                }
                // du_total already contains its flow/Hessian state chain.
                for (column, derivative) in gradient.iter_mut().enumerate() {
                    let mut sum = Dd::from_f64(*derivative);
                    for (port, &flow) in flows.iter().enumerate() {
                        sum = sum
                            + Dd::from_f64(flow)
                                * Dd::from_f64(feedback_derivative[port * n + column]);
                    }
                    *derivative = sum.to_f64();
                }
            }
            Kind::Resistance => {
                *value = self.system.resistance_power_from_effort(effort)?;
                self.system
                    .resistance_power_effort_gradient_into(effort, field)?;
                if !self.tangent.apply_hessian(state, field, gradient) {
                    return Err(DissipativeError::SourceRefused.into());
                }
            }
            Kind::Nonlinear => {
                *value = nonlinear_power;
                state_direction.fill(0.0);
                for column in 0..n {
                    state_direction[column] = 1.0;
                    if !self
                        .tangent
                        .apply_hessian(state, state_direction, effort_direction)
                    {
                        return Err(DissipativeError::SourceRefused.into());
                    }
                    resistance_direction.fill(f64::NAN);
                    if !(self.dissipation.tangent)(
                        state,
                        effort,
                        state_direction,
                        effort_direction,
                        resistance_direction,
                    ) || resistance_direction.iter().any(|v| !v.is_finite())
                    {
                        return Err(LinearizationError::DissipationRefused {
                            column: Some(column),
                        }
                        .into());
                    }
                    let mut sum = Dd::ZERO;
                    for i in 0..n {
                        sum = sum + Dd::from_f64(effort_direction[i]) * Dd::from_f64(resistance[i]);
                        sum = sum + Dd::from_f64(effort[i]) * Dd::from_f64(resistance_direction[i]);
                    }
                    gradient[column] = sum.to_f64();
                    state_direction[column] = 0.0;
                }
            }
        }
        // Work per cycle includes q_omega=-q/omega at fixed coefficients.
        // The bordered dual supplies the unknown omega_p, exactly once.
        if !objective
            .normalization()
            .normalize(omega, value, gradient, Some(frequency))
        {
            return Err(DissipativeError::SourceRefused.into());
        }
        Ok(())
    }
}

impl<'s, 'a> PreparedDissipativeAutonomousSensitivity<'s, 'a> {
    /// Sample the original account and its complete coefficient/frequency
    /// covector on both native AFT grids. Nonlinear heat uses
    /// `de dot D + e dot (D_x dx + D_e de)`, de=Hxx dx. Every memory and
    /// prestress cross term remains. Pointwise positive heat does not impose
    /// a positive directional slope or certify a thermal destination.
    pub fn power_objective(
        &mut self,
        objective: DissipativeAutonomousPowerObjective,
        controls: ForcedMeanObjectiveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousMeanObjective> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let source = self.source;
        let omega = self.prepared.angular_frequency();
        let mut error = None;
        let result = self.prepared.mean_objective(
            controls,
            |_time, state, value, gradient, frequency| match source
                .observe_power(objective, omega, state, value, gradient, frequency)
            {
                Ok(()) => true,
                Err(cause) => {
                    if error.is_none() {
                        error = Some(cause);
                    }
                    false
                }
            },
            &mut cancelled,
        );
        if matches!(result, Err(ForcedSensitivityError::Cancelled)) {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        if let Some(error) = error {
            return Err(error);
        }
        result.map_err(Into::into)
    }

    /// Reuse the original bordered factor for one complete physical account.
    /// Source storage/feedback/loss tangents cannot be replaced at this query.
    /// Parameter queries reuse both cached trajectories and the same dual.
    pub fn prepare_power_adjoint(
        &mut self,
        objective: DissipativeAutonomousPowerObjective,
        controls: ForcedMeanObjectiveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDissipativeAutonomousPowerAdjoint<'_, 's, 'a>> {
        let observation = self.power_objective(objective, controls, &mut cancelled)?;
        let omega = self.prepared.angular_frequency();
        let prepared = self.prepared.prepare_adjoint(
            &observation.gradient,
            observation.angular_frequency_partial,
            controls.objective_scale,
            &mut cancelled,
        )?;
        Ok(PreparedDissipativeAutonomousPowerAdjoint {
            source: self.source,
            objective,
            omega,
            value: observation.value,
            angular_frequency_partial: observation.angular_frequency_partial,
            value_sampling_defect: observation.value_sampling_defect,
            gradient_sampling_defect: observation.gradient_sampling_defect,
            prepared,
        })
    }
}

/// One original power/work observation and its complete physical bordered dual.
pub struct PreparedDissipativeAutonomousPowerAdjoint<'p, 's, 'a> {
    source: &'s DissipativeAutonomousHamiltonian<'a>,
    objective: DissipativeAutonomousPowerObjective,
    omega: f64,
    value: f64,
    angular_frequency_partial: f64,
    value_sampling_defect: f64,
    gradient_sampling_defect: f64,
    prepared: PreparedAutonomousAdjoint<'p>,
}
impl PreparedDissipativeAutonomousPowerAdjoint<'_, '_, '_> {
    pub fn objective(&self) -> DissipativeAutonomousPowerObjective {
        self.objective
    }
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn angular_frequency_partial(&self) -> f64 {
        self.angular_frequency_partial
    }
    pub fn value_sampling_defect(&self) -> f64 {
        self.value_sampling_defect
    }
    pub fn gradient_sampling_defect(&self) -> f64 {
        self.gradient_sampling_defect
    }
    pub fn prepared(&self) -> &PreparedAutonomousAdjoint<'_> {
        &self.prepared
    }

    /// Complete physical q_p and residual R_p at fixed state/frequency.
    /// For D heat, q_p=e_p dot D+e dot(D_e e_p+D_p); the original law supplies
    /// all terms. Each grid's direct mean participates in its own adjoint
    /// contraction. The unknown waveform AND period response stay in the dual.
    pub fn gradient_parameter(
        &mut self,
        parameter_scale: f64,
        direction: DissipativeAutonomousParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AutonomousObjectiveGradient> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let (objective, omega) = (self.objective, self.omega);
        parameter::with_samples(self.source, direction, |samples| {
            self.prepared.gradient_mean_island_direction(
                parameter_scale,
                |_time, state, out| {
                    samples
                        .borrow_mut()
                        .objective_partial(state, objective, omega, out)
                },
                |_time, state, out| samples.borrow_mut().field(state, out),
                &mut cancelled,
            )
        })
    }
}
