//! Independent original pump, port and resistance power/work objectives.

use super::super::{Scratch, add_matrix_action, quadratic};
use super::{
    ParametricError, ParametricHamiltonian, ParametricParameterDirection,
    PreparedParametricSensitivity, Result, direction, ode_clock_column,
};
use crate::phs_orbit::sensitivity::power::PhsPowerObjective;
use fs_math::dd::Dd;
use fs_orbit::sensitivity::{
    ForcedMeanObjective, ForcedMeanObjectiveControls, ForcedObjectiveGradient,
    ForcedSensitivityDirection, ForcedSensitivityError, PreparedForcedAdjoint,
};
use fs_phs::PreparedPhsParameterDirection;

/// Independent original physical accounts. Pump and port power enter the
/// Hamiltonian system; resistance power leaves it. A finite-order balance
/// defect stays observable; no account is reconstructed from the other two.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ParametricPowerObjective {
    MeanPumpPower,
    MeanPortPower,
    MeanResistancePower,
    PumpEnergyPerCycle,
    PortEnergyPerCycle,
    ResistanceEnergyPerCycle,
}

#[derive(Clone, Copy)]
enum PowerKind {
    Pump,
    Port,
    Resistance,
}

impl ParametricPowerObjective {
    fn kind(self) -> PowerKind {
        match self {
            Self::MeanPumpPower | Self::PumpEnergyPerCycle => PowerKind::Pump,
            Self::MeanPortPower | Self::PortEnergyPerCycle => PowerKind::Port,
            Self::MeanResistancePower | Self::ResistanceEnergyPerCycle => PowerKind::Resistance,
        }
    }
    // Pump and port work share the existing physical supply normalization. The
    // actual laws below remain separate; only power/period units are shared.
    fn normalization(self) -> PhsPowerObjective {
        match self {
            Self::MeanPumpPower | Self::MeanPortPower => PhsPowerObjective::MeanSuppliedPower,
            Self::MeanResistancePower => PhsPowerObjective::MeanResistancePower,
            Self::PumpEnergyPerCycle | Self::PortEnergyPerCycle => {
                PhsPowerObjective::SuppliedEnergyPerCycle
            }
            Self::ResistanceEnergyPerCycle => PhsPowerObjective::ResistanceEnergyPerCycle,
        }
    }

    fn direct_partial(
        self,
        source: &ParametricHamiltonian<'_>,
        chart: ForcedSensitivityDirection,
        state: &[f64],
        physical: &PreparedPhsParameterDirection<'_, '_>,
        base: &Scratch,
        change: &Scratch,
    ) -> Result<f64> {
        let omega = source.angular_frequency();
        let clock = chart.angular_frequency_derivative;
        if matches!(self.kind(), PowerKind::Pump) {
            // Qdot contains the physical clock even at fixed phase. For work
            // per cycle its explicit clock derivative cancels T_p analytically;
            // the separate unknown waveform/clock chain is still in the dual.
            let mut derivative = quadratic(&change.q_rate, state) * Dd::from_f64(0.5);
            if self.normalization().per_cycle() {
                derivative = derivative * Dd::from_f64(core::f64::consts::TAU / omega);
            } else if clock != 0.0 {
                derivative = derivative
                    + quadratic(&base.q_rate, state) * Dd::from_f64(0.5) * Dd::from_f64(clock)
                        / Dd::from_f64(omega);
            }
            let value = derivative.to_f64();
            return if value.is_finite() {
                Ok(value)
            } else {
                Err(ParametricError::SourceRefused)
            };
        }
        let direct = match self.kind() {
            PowerKind::Port => physical.supplied_power_from_effort(
                &base.effort,
                &change.effort,
                &base.ports,
                &change.ports,
            )?,
            PowerKind::Resistance => {
                physical.resistance_power_from_effort(&base.effort, &change.effort)?
            }
            PowerKind::Pump => unreachable!(),
        };
        let original = if self.normalization().per_cycle() && clock != 0.0 {
            match self.kind() {
                PowerKind::Port => source
                    .system
                    .supplied_power_from_effort(&base.effort, &base.ports)?,
                PowerKind::Resistance => {
                    source.system.resistance_power_from_effort(&base.effort)?
                }
                PowerKind::Pump => unreachable!(),
            }
        } else {
            0.0
        };
        self.normalization()
            .direct_partial(omega, clock, direct, original)
            .ok_or(ParametricError::SourceRefused)
    }
}

impl ParametricHamiltonian<'_> {
    fn observe_power(
        &self,
        objective: ParametricPowerObjective,
        time: f64,
        state: &[f64],
        value: &mut f64,
        gradient: &mut [f64],
    ) -> Result<()> {
        let mut scratch = self.scratch.borrow_mut();
        self.effort(
            time,
            state,
            &mut scratch,
            matches!(objective.kind(), PowerKind::Pump),
        )?;
        let Scratch {
            q,
            q_rate,
            effort,
            ports,
            flow,
            ..
        } = &mut *scratch;
        match objective.kind() {
            PowerKind::Pump => {
                *value = (quadratic(q_rate, state) * Dd::from_f64(0.5)).to_f64();
                gradient.fill(0.0);
                add_matrix_action(q_rate, state, gradient)?;
            }
            kind => {
                match kind {
                    PowerKind::Port => {
                        *value = self.system.supplied_power_from_effort(effort, ports)?;
                        self.system.input_flow_into(ports, flow)?;
                    }
                    PowerKind::Resistance => {
                        *value = self.system.resistance_power_from_effort(effort)?;
                        self.system
                            .resistance_power_effort_gradient_into(effort, flow)?;
                    }
                    PowerKind::Pump => unreachable!(),
                }
                if !self.tangent.apply_hessian(state, flow, gradient) {
                    return Err(ParametricError::SourceRefused);
                }
                add_matrix_action(q, flow, gradient)?;
            }
        }
        if !objective
            .normalization()
            .normalize(self.angular_frequency(), value, gradient, None)
        {
            return Err(ParametricError::SourceRefused);
        }
        Ok(())
    }
}

impl<'s, 'a> PreparedParametricSensitivity<'s, 'a> {
    /// Original physical value and complete real Fourier covector, sampled on
    /// the same independent AFT grids as the prepared response. Pump uses
    /// `xᵀQdot*x/2` with state gradient `Qdot*x`; port and loss use the COMPLETE
    /// effort `grad H0+Q*x` and storage Hessian `H0_xx+Q`, including memories.
    ///
    /// Requested power/energy units condition both value and full covector
    /// sampling checks. Physical parameter and prescribed-clock partials are
    /// assembled automatically by `prepare_power_adjoint`'s parameter queries.
    pub fn power_objective(
        &mut self,
        objective: ParametricPowerObjective,
        controls: ForcedMeanObjectiveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedMeanObjective> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let source = self.source;
        let mut error = None;
        let result = self.prepared.mean_objective(
            controls,
            |time, state, value, gradient| match source
                .observe_power(objective, time, state, value, gradient)
            {
                Ok(()) => true,
                Err(cause) => {
                    error = Some(cause);
                    false
                }
            },
            &mut cancelled,
        );
        if let Some(error) = error {
            return Err(error);
        }
        result.map_err(Into::into)
    }

    /// Prepare one complete original pump/port/loss objective transpose solve.
    /// Uses the H0 tangent already bound by this source; no independently
    /// replaceable drive, pump or tangent can be supplied at this query.
    pub fn prepare_power_adjoint(
        &mut self,
        objective: ParametricPowerObjective,
        controls: ForcedMeanObjectiveControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedParametricPowerAdjoint<'_, 's, 'a>> {
        let observation = self.power_objective(objective, controls, &mut cancelled)?;
        let harmonics = self.prepared.harmonics();
        let prepared = self.prepared.prepare_adjoint(
            &observation.gradient,
            controls.objective_scale,
            &mut cancelled,
        )?;
        Ok(PreparedParametricPowerAdjoint {
            source: self.source,
            harmonics,
            objective,
            value: observation.value,
            value_sampling_defect: observation.value_sampling_defect,
            gradient_sampling_defect: observation.gradient_sampling_defect,
            prepared,
        })
    }
}

/// One original pumped physical objective, reused for complete physical
/// parameter directions without new state solves, factorization or Hessian
/// evaluations. Original constitutive/source callbacks still sample the cached
/// physical states; pump, port and resistance accounts retain their own laws.
pub struct PreparedParametricPowerAdjoint<'p, 's, 'a> {
    source: &'s ParametricHamiltonian<'a>,
    harmonics: usize,
    objective: ParametricPowerObjective,
    value: f64,
    value_sampling_defect: f64,
    gradient_sampling_defect: f64,
    prepared: PreparedForcedAdjoint<'p>,
}

impl PreparedParametricPowerAdjoint<'_, '_, '_> {
    pub fn objective(&self) -> ParametricPowerObjective {
        self.objective
    }
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn value_sampling_defect(&self) -> f64 {
        self.value_sampling_defect
    }
    pub fn gradient_sampling_defect(&self) -> f64 {
        self.gradient_sampling_defect
    }
    pub fn prepared(&self) -> &PreparedForcedAdjoint<'_> {
        &self.prepared
    }

    /// Automatic TOTAL derivative along original J/R/G/H0, Q Fourier, input
    /// Fourier and clock changes. Each AFT grid contracts its own direct scalar
    /// mean and complete physical residual direction with the unchanged dual.
    /// Explicit pump power contains Qdot's clock derivative; pump work per
    /// cycle cancels that explicit clock term against the period derivative.
    /// The original sI residual clock column remains in both cases, exactly once.
    pub fn gradient_parameter(
        &mut self,
        chart: ForcedSensitivityDirection,
        direction: ParametricParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedObjectiveGradient> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let source = self.source;
        let objective = self.objective;
        let n = source.system.state_dim();
        direction::with_samples(source, self.harmonics, direction, |samples| {
            self.prepared.gradient_mean_split_direction(
                chart,
                |time, state, out| {
                    samples.borrow_mut().scalar(
                        time,
                        state,
                        matches!(objective.kind(), PowerKind::Pump),
                        |physical, base, change| {
                            objective.direct_partial(source, chart, state, physical, base, change)
                        },
                        out,
                    )
                },
                |time, state, out| samples.borrow_mut().field(time, state, out),
                |point, out| ode_clock_column(n, chart.angular_frequency_derivative, point, out),
                &mut cancelled,
            )
        })
    }
}
