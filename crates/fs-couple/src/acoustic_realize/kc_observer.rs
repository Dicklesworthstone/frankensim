//! Endpoint radiation from the complete fixed-end string storage law.
//!
//! The accepted model has constant J, R and G over a held-input interval:
//! x' = (J-R) grad H + G u, and x'' = (J-R) D(grad H) x'. Its canonical,
//! mass-normalized momentum rows are modal acceleration and jerk. Asking the
//! original storage for its tangent retains stretching, elastic contact and
//! every bending-memory coordinate without reconstructing their constitutive
//! laws here. This is an exact derivative at the numerical endpoint, not an
//! exact nonlinear trajectory or a derivative of the discrete step map.

use super::{AcousticRealizeError, StringObserver, render_phs_error};
use fs_phs::PortHamiltonian;

/// Reused scratch for the original flow and its directional derivative.
/// The admitted KC, contact and linear bending-memory actions are borrowed and
/// allocate nothing during observation; all buffers are prepared once.
pub(super) struct KcObservation {
    effort: Vec<f64>,
    state_rate: Vec<f64>,
    effort_rate: Vec<f64>,
    state_second: Vec<f64>,
    zero_input: Vec<f64>,
}

impl KcObservation {
    pub(super) fn new(system: &PortHamiltonian) -> Result<Self, AcousticRealizeError> {
        if system.storage().continuous_tangent().is_none() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "KC endpoint radiation requires the complete storage tangent",
            });
        }
        let n = system.state_dim();
        Ok(Self {
            effort: vec![0.0; n],
            state_rate: vec![0.0; n],
            effort_rate: vec![0.0; n],
            state_second: vec![0.0; n],
            zero_input: vec![0.0; system.port_dim()],
        })
    }

    /// Smooth limit at the end of this sample's held force. Bow, tangential
    /// friction and contact loss are already prescribed ZOH input here; this
    /// does not differentiate an unimplemented continuously varying drive law.
    /// Ideal input jumps need a resolved transition; a contact corner may have
    /// no classical endpoint derivative. A declared tangent refusal never
    /// falls back to differences.
    pub(super) fn observe(
        &mut self,
        observer: &StringObserver,
        system: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
    ) -> Result<f64, AcousticRealizeError> {
        let n = self.effort.len();
        if system.state_dim() != n
            || state.len() != n
            || system.port_dim() != observer.weights.len()
            || input.len() != self.zero_input.len()
            || system.port_dim() != input.len()
            || observer
                .weights
                .len()
                .checked_mul(2)
                .is_none_or(|prefix| prefix > n)
            || state.iter().chain(input).any(|value| !value.is_finite())
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "KC endpoint radiation needs complete finite state and held input",
            });
        }
        let storage = system.storage();
        let tangent =
            storage
                .continuous_tangent()
                .ok_or(AcousticRealizeError::InvalidDescription {
                    what: "KC endpoint radiation requires the complete storage tangent",
                })?;
        self.effort.fill(f64::NAN);
        storage.gradient(state, &mut self.effort);
        system
            .vector_field_from_effort_into(&self.effort, input, &mut self.state_rate)
            .map_err(render_phs_error)?;
        self.effort_rate.fill(f64::NAN);
        if !tangent.apply(state, &self.state_rate, &mut self.effort_rate)
            || self.effort_rate.iter().any(|value| !value.is_finite())
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "KC endpoint radiation requires a finite smooth storage tangent",
            });
        }
        system
            .vector_field_from_effort_into(
                &self.effort_rate,
                &self.zero_input,
                &mut self.state_second,
            )
            .map_err(render_phs_error)?;
        let pressure = observer
            .weights
            .iter()
            .enumerate()
            .map(|(mode, &(area, moment))| {
                let row = 2 * mode + 1;
                area * self.state_rate[row] + moment * self.state_second[row]
            })
            .sum();
        StringObserver::check_pressure(pressure)
    }
}

#[cfg(test)]
mod tests;
