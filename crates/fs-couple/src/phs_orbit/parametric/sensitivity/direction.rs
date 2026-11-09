//! One original sample chain for pumped residual and objective directions.

use super::super::{Scratch, add_matrix_action};
use super::{ParametricError, ParametricHamiltonian, ParametricParameterDirection, Result};
use core::cell::RefCell;
use fs_orbit::sensitivity::ForcedSensitivityError;
use fs_phs::PreparedPhsParameterDirection;
use fs_phs::harmonic::equilibrium::{
    DissipationParameterAction, DissipationParameterEvaluation,
    PreparedDissipationParameterDirection,
};

pub(super) fn with_samples<'s, 'a, 'd, T>(
    source: &'s ParametricHamiltonian<'a>,
    harmonics: usize,
    direction: ParametricParameterDirection<'d>,
    dissipation: Option<&'d DissipationParameterAction<'d>>,
    run: impl FnOnce(&RefCell<Samples<'s, 'a, 'd>>) -> core::result::Result<T, ForcedSensitivityError>,
) -> Result<T> {
    let (n, m) = (source.system.state_dim(), source.system.port_dim());
    if let Some(pump) = direction.modulation {
        if pump.dimension != n
            || pump.omega != source.modulation.omega
            || pump.harmonics > harmonics / source.periods
        {
            return Err(ParametricError::Invalid(
                "pump direction requires original dimension, clock and retained bandwidth",
            ));
        }
    }
    if let Some(input) = direction.physical.input {
        if input.port_count() != m
            || input.omega_rad_s() != source.modulation.omega
            || input.highest_nonzero_harmonic() > harmonics / source.periods
        {
            return Err(ParametricError::Invalid(
                "input direction requires original ports, clock and retained bandwidth",
            ));
        }
    }
    let physical = source
        .system
        .prepare_parameter_direction(direction.physical.structure)?;
    let loss = match source.dissipation {
        Some(law) => Some(
            law.prepare_parameter_direction(n, dissipation)
                .map_err(ParametricError::Dissipation)?,
        ),
        None if dissipation.is_some() => {
            return Err(ParametricError::Invalid(
                "a nonlinear loss direction requires the original attached loss law",
            ));
        }
        None => None,
    };
    let pump_order = direction.modulation.map_or(0, |pump| pump.harmonics);
    let samples = RefCell::new(Samples {
        source,
        physical,
        loss,
        direction,
        scratch: Scratch::new(n, m, pump_order),
        error: None,
    });
    let result = run(&samples);
    if let Some(error) = samples.into_inner().error {
        return Err(error);
    }
    result.map_err(Into::into)
}

// Both consumers run sequentially through one direction scratch, while the
// existing source retains its own original Q/effort/input scratch and sampler.
pub(super) struct Samples<'s, 'a, 'd> {
    source: &'s ParametricHamiltonian<'a>,
    physical: PreparedPhsParameterDirection<'a, 'd>,
    loss: Option<PreparedDissipationParameterDirection<'a, 'd>>,
    direction: ParametricParameterDirection<'d>,
    scratch: Scratch,
    error: Option<ParametricError>,
}

impl Samples<'_, '_, '_> {
    fn sample<T>(
        &mut self,
        time: f64,
        state: &[f64],
        rate: bool,
        evaluate_loss: bool,
        observe: impl FnOnce(
            &PreparedPhsParameterDirection<'_, '_>,
            &Scratch,
            &Scratch,
            Option<&DissipationParameterEvaluation<'_>>,
        ) -> Result<T>,
    ) -> Result<T> {
        let mut original = self.source.scratch.borrow_mut();
        self.source.effort(time, state, &mut original, rate)?;
        if let Some(pump) = self.direction.modulation {
            pump.sample(time, &mut self.scratch, rate)?;
        }
        self.scratch.effort.fill(0.0);
        if let Some(storage) = self.direction.physical.storage {
            self.scratch.effort.fill(f64::NAN);
            if !storage(state, &mut self.scratch.effort)
                || self.scratch.effort.iter().any(|v| !v.is_finite())
            {
                return Err(ParametricError::SourceRefused);
            }
        }
        add_matrix_action(&self.scratch.q, state, &mut self.scratch.effort)?;
        self.scratch.ports.fill(0.0);
        if let Some(input) = self.direction.physical.input {
            if !input.sample_into(time, &mut self.scratch.ports) {
                return Err(ParametricError::SourceRefused);
            }
        }
        let loss = self
            .loss
            .as_mut()
            .filter(|_| evaluate_loss)
            .map(|loss| loss.evaluate(state, &original.effort, &self.scratch.effort, || false))
            .transpose()
            .map_err(ParametricError::Dissipation)?;
        observe(&self.physical, &original, &self.scratch, loss.as_ref())
    }

    pub(super) fn field(&mut self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        match self.sample(time, state, false, true, |physical, base, change, loss| {
            if let Some(loss) = loss {
                physical
                    .vector_field_with_dissipation_from_effort_into(
                        &base.effort,
                        &change.effort,
                        &base.ports,
                        &change.ports,
                        loss.effort_chain,
                        loss.direct_flow,
                        out,
                    )
                    .map_err(Into::into)
            } else {
                physical
                    .vector_field_from_effort_into(
                        &base.effort,
                        &change.effort,
                        &base.ports,
                        &change.ports,
                        out,
                    )
                    .map_err(Into::into)
            }
        }) {
            Ok(()) => true,
            Err(error) => {
                self.error = Some(error);
                false
            }
        }
    }

    pub(super) fn scalar(
        &mut self,
        time: f64,
        state: &[f64],
        rate: bool,
        evaluate_loss: bool,
        observe: impl FnOnce(
            &PreparedPhsParameterDirection<'_, '_>,
            &Scratch,
            &Scratch,
            Option<&DissipationParameterEvaluation<'_>>,
        ) -> Result<f64>,
        out: &mut f64,
    ) -> bool {
        match self.sample(time, state, rate, evaluate_loss, observe) {
            Ok(value) if value.is_finite() => {
                *out = value;
                true
            }
            Ok(_) => {
                self.error = Some(ParametricError::SourceRefused);
                false
            }
            Err(error) => {
                self.error = Some(error);
                false
            }
        }
    }
}
