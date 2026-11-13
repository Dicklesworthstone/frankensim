//! One original fixed-state source chain for waveform and objective queries.

use super::{
    DissipativeAutonomousHamiltonian, DissipativeAutonomousParameterDirection,
    DissipativeAutonomousPowerObjective, DissipativeAutonomousSensitivityError, DissipativeError,
    ForcedSensitivityError, Result,
};
use core::cell::RefCell;
use fs_phs::PreparedPhsParameterDirection;
use fs_phs::harmonic::equilibrium::PreparedDissipationParameterDirection;

pub(super) fn with_samples<'s, 'a, 'd, T>(
    source: &'s DissipativeAutonomousHamiltonian<'a>,
    direction: DissipativeAutonomousParameterDirection<'d>,
    run: impl FnOnce(&RefCell<Samples<'s, 'a, 'd>>) -> core::result::Result<T, ForcedSensitivityError>,
) -> Result<T> {
    let (n, m) = (source.system.state_dim(), source.system.port_dim());
    let physical = source
        .system
        .prepare_parameter_direction(direction.physical.structure)?;
    let loss = source
        .dissipation
        .prepare_parameter_direction(n, direction.dissipation)?;
    let samples = RefCell::new(Samples {
        source,
        direction,
        physical,
        loss,
        effort: vec![f64::NAN; n],
        partial_effort: vec![0.0; n],
        flows: vec![f64::NAN; m],
        partial_flows: vec![f64::NAN; m],
        input: vec![f64::NAN; m],
        partial_input: vec![f64::NAN; m],
        error: None,
    });
    let result = run(&samples);
    if matches!(result, Err(ForcedSensitivityError::Cancelled)) {
        return Err(ForcedSensitivityError::Cancelled.into());
    }
    if let Some(error) = samples.into_inner().error {
        return Err(error);
    }
    result.map_err(Into::into)
}

pub(super) struct Samples<'s, 'a, 'd> {
    source: &'s DissipativeAutonomousHamiltonian<'a>,
    direction: DissipativeAutonomousParameterDirection<'d>,
    physical: PreparedPhsParameterDirection<'a, 'd>,
    loss: PreparedDissipationParameterDirection<'a, 'd>,
    effort: Vec<f64>,
    partial_effort: Vec<f64>,
    flows: Vec<f64>,
    partial_flows: Vec<f64>,
    input: Vec<f64>,
    partial_input: Vec<f64>,
    error: Option<DissipativeAutonomousSensitivityError>,
}

impl Samples<'_, '_, '_> {
    fn load(&mut self, state: &[f64]) -> Result<()> {
        if state.len() != self.effort.len() || state.iter().any(|v| !v.is_finite()) {
            return Err(DissipativeError::SourceRefused.into());
        }
        self.effort.fill(f64::NAN);
        self.source
            .system
            .storage()
            .gradient(state, &mut self.effort);
        self.partial_effort.fill(0.0);
        if let Some(storage) = self.direction.physical.storage {
            self.partial_effort.fill(f64::NAN);
            if !storage(state, &mut self.partial_effort) {
                return Err(DissipativeError::SourceRefused.into());
            }
        }
        if self
            .effort
            .iter()
            .chain(&self.partial_effort)
            .any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused.into());
        }
        self.source
            .system
            .output_from_effort_into(&self.effort, &mut self.flows)?;
        self.input.fill(f64::NAN);
        if !self
            .source
            .feedback
            .input(state, &self.flows, &mut self.input)
            || self.input.iter().any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused.into());
        }
        self.physical.output_from_effort_into(
            &self.effort,
            &self.partial_effort,
            &mut self.partial_flows,
        )?;
        self.partial_input.fill(f64::NAN);
        if !(self.direction.physical.feedback)(
            state,
            &self.flows,
            &self.partial_flows,
            &mut self.partial_input,
        ) || self.partial_input.iter().any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused.into());
        }
        Ok(())
    }

    fn record(&mut self, result: Result<()>) -> bool {
        if let Err(error) = result {
            if self.error.is_none() {
                self.error = Some(error);
            }
            false
        } else {
            true
        }
    }

    pub(super) fn field(&mut self, state: &[f64], out: &mut [f64]) -> bool {
        let result = (|| {
            self.load(state)?;
            let loss = self
                .loss
                .evaluate(state, &self.effort, &self.partial_effort, || false)?;
            // Keep the independent D_e*e_p and D_p contributions separate until
            // the same physical row owner rounds ALL terms exactly once.
            self.physical
                .vector_field_with_dissipation_from_effort_into(
                    &self.effort,
                    &self.partial_effort,
                    &self.input,
                    &self.partial_input,
                    loss.effort_chain,
                    loss.direct_flow,
                    out,
                )?;
            Ok(())
        })();
        self.record(result)
    }

    pub(super) fn objective_partial(
        &mut self,
        state: &[f64],
        objective: DissipativeAutonomousPowerObjective,
        omega: f64,
        out: &mut f64,
    ) -> bool {
        let result = (|| {
            self.load(state)?;
            let direct = match objective.kind() {
                super::power::Kind::Supplied => self.physical.supplied_power_from_effort(
                    &self.effort,
                    &self.partial_effort,
                    &self.input,
                    &self.partial_input,
                )?,
                super::power::Kind::Resistance => self
                    .physical
                    .resistance_power_from_effort(&self.effort, &self.partial_effort)?,
                super::power::Kind::Nonlinear => {
                    self.loss
                        .evaluate(state, &self.effort, &self.partial_effort, || false)?
                        .power_derivative
                }
            };
            // This is q_p at fixed coefficients AND frequency. The bordered
            // dual owns the unknown period change through q_omega exactly once.
            *out = objective
                .normalization()
                .direct_partial(omega, 0.0, direct, 0.0)
                .ok_or(DissipativeError::SourceRefused)?;
            Ok(())
        })();
        self.record(result)
    }
}
