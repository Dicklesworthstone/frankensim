//! Shared original physical sample assembly for waveform and power derivatives.

use super::{
    AutonomousPortFeedback, ForcedSensitivityError, PhsAutonomousParameterDirection, PhsOrbitError,
    Result,
};
use crate::phs_orbit::sensitivity::power::PhsPowerObjective;
use fs_phs::{PhsError, PortHamiltonian, PreparedPhsParameterDirection};
use std::cell::RefCell;

pub(super) fn with_parameter_samples<'s, 'd, T>(
    system: &'s PortHamiltonian,
    feedback: &'s dyn AutonomousPortFeedback,
    direction: PhsAutonomousParameterDirection<'d>,
    run: impl FnOnce(
        &RefCell<ParameterSamples<'s, 'd>>,
    ) -> core::result::Result<T, ForcedSensitivityError>,
) -> Result<T> {
    let physical = system.prepare_parameter_direction(direction.structure)?;
    let (n, m) = (system.state_dim(), system.port_dim());
    let samples = RefCell::new(ParameterSamples {
        system,
        feedback,
        direction,
        physical,
        effort: vec![f64::NAN; n],
        partial_effort: vec![0.0; n],
        flows: vec![f64::NAN; m],
        partial_flows: vec![f64::NAN; m],
        input: vec![f64::NAN; m],
        partial_input: vec![f64::NAN; m],
        source_refused: false,
        physical_error: None,
    });
    let result = run(&samples);
    let samples = samples.into_inner();
    if samples.source_refused {
        return Err(PhsOrbitError::SourceRefused.into());
    }
    if let Some(error) = samples.physical_error {
        return Err(error.into());
    }
    result.map_err(Into::into)
}

// The scalar and vector callbacks are sequential consumers of this one scratch
// owner. In particular, BOTH use the same complete Gp/Hp -> yp -> feedback chain.
pub(super) struct ParameterSamples<'s, 'd> {
    system: &'s PortHamiltonian,
    feedback: &'s dyn AutonomousPortFeedback,
    direction: PhsAutonomousParameterDirection<'d>,
    physical: PreparedPhsParameterDirection<'s, 'd>,
    effort: Vec<f64>,
    partial_effort: Vec<f64>,
    flows: Vec<f64>,
    partial_flows: Vec<f64>,
    input: Vec<f64>,
    partial_input: Vec<f64>,
    source_refused: bool,
    physical_error: Option<PhsError>,
}

impl ParameterSamples<'_, '_> {
    fn load(&mut self, state: &[f64], needs_input: bool) -> bool {
        self.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut self.effort);
        if self.effort.iter().any(|v| !v.is_finite()) {
            self.source_refused = true;
            return false;
        }
        if let Some(storage) = self.direction.storage {
            self.partial_effort.fill(f64::NAN);
            if !storage(state, &mut self.partial_effort)
                || self.partial_effort.iter().any(|v| !v.is_finite())
            {
                return false;
            }
        }
        if !needs_input {
            return true;
        }
        if self
            .system
            .output_from_effort_into(&self.effort, &mut self.flows)
            .is_err()
        {
            self.source_refused = true;
            return false;
        }
        self.input.fill(f64::NAN);
        if !self.feedback.input(state, &self.flows, &mut self.input)
            || self.input.iter().any(|v| !v.is_finite())
        {
            self.source_refused = true;
            return false;
        }
        if let Err(error) = self.physical.output_from_effort_into(
            &self.effort,
            &self.partial_effort,
            &mut self.partial_flows,
        ) {
            self.physical_error = Some(error);
            return false;
        }
        self.partial_input.fill(f64::NAN);
        if !(self.direction.feedback)(
            state,
            &self.flows,
            &self.partial_flows,
            &mut self.partial_input,
        ) || self.partial_input.iter().any(|v| !v.is_finite())
        {
            return false;
        }
        true
    }

    pub(super) fn vector_field(&mut self, state: &[f64], out: &mut [f64]) -> bool {
        if !self.load(state, true) {
            return false;
        }
        match self.physical.vector_field_from_effort_into(
            &self.effort,
            &self.partial_effort,
            &self.input,
            &self.partial_input,
            out,
        ) {
            Ok(()) => true,
            Err(error) => {
                self.physical_error = Some(error);
                false
            }
        }
    }

    pub(super) fn objective_partial(
        &mut self,
        state: &[f64],
        objective: PhsPowerObjective,
        omega: f64,
        out: &mut f64,
    ) -> bool {
        if !self.load(state, objective.supplied()) {
            return false;
        }
        let power = if objective.supplied() {
            self.physical.supplied_power_from_effort(
                &self.effort,
                &self.partial_effort,
                &self.input,
                &self.partial_input,
            )
        } else {
            self.physical
                .resistance_power_from_effort(&self.effort, &self.partial_effort)
        };
        let direct = match power {
            Ok(value) => value,
            Err(error) => {
                self.physical_error = Some(error);
                return false;
            }
        };
        // The original bordered objective owns the UNKNOWN frequency chain.
        // No solved Omega_p belongs in this fixed-frequency direct partial.
        let Some(value) = objective.direct_partial(omega, 0.0, direct, 0.0) else {
            return false;
        };
        *out = value;
        true
    }
}
