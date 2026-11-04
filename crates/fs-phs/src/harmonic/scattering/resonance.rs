//! Driven peaks of the physically closed system, including interference
//! between an incident wave and the resonator's outgoing wave.

use super::{ScatteringFrequencyResponse, ScatteringModel, ScatteringResponse, finite, poll};
use crate::harmonic::resonance::{
    DrivenResonance, DrivenResonanceError, ResonanceControls, project, search_resonance,
};
use fs_math::c64::C64;

/// Fixed real observation weights in the original physical units. Coherent
/// complex amplitudes are combined before taking a magnitude; mutual phase,
/// direct incident-wave interference and ordinary-source contributions survive.
#[derive(Clone, Copy, Debug)]
pub enum ScatteringObservation<'a> {
    /// All original differential states and algebraic multipliers.
    State(&'a [f64]),
    /// All original storage efforts and algebraic constraint reactions.
    Effort(&'a [f64]),
    /// Actual conjugate inputs, including `(2a-y)/Z` at selected wave ports.
    /// Weights follow original source order, including ordinary input ports.
    Input(&'a [f64]),
    /// Actual conjugate outputs in original source order.
    Output(&'a [f64]),
    /// Outgoing waves `b=y-a`, in [`ScatteringModel::wave_ports`] order.
    /// A reflected-wave peak can differ from an internal-motion peak because
    /// its incident feedthrough remains in this complex observation.
    Outgoing(&'a [f64]),
}

impl ScatteringModel {
    /// Resolve one driven peak of the complete wave-loaded physical model.
    ///
    /// `input` has one fixed incident output-wave amplitude at each selected
    /// original source port, and one fixed ordinary input at every other port.
    /// The references stay fixed. All boundary feedback is closed before every
    /// factorization, retaining bare resonances, material memory and constraints.
    /// Physical frequency tangents check the original equations with actual
    /// boundary inputs, rather than treating the effective `2a/Z` drive as work.
    ///
    /// The caller must isolate the desired maximum with positive amplitude
    /// slope at the lower endpoint and negative slope at the upper endpoint.
    /// If requested, both -3 dB crossings must also lie in that interval.
    /// The shared [`HarmonicModel::driven_resonance`](crate::harmonic::HarmonicModel::driven_resonance)
    /// search owns the bounded refinement and all numerical/no-uniqueness limits.
    /// Observed bandwidth/Q describes this selected amplitude, not material heat
    /// or a free-pole decay rate. Reflection minima are not silently called peaks.
    pub fn driven_resonance(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, DrivenResonanceError> {
        poll(&mut cancelled)?;
        if input.len() != self.port_count() || input.iter().any(|&v| !finite(v)) {
            return Err(DrivenResonanceError::Invalid(
                "one finite incident wave or ordinary input per original port required",
            ));
        }
        let (weights, dimension) = match observation {
            ScatteringObservation::State(w) | ScatteringObservation::Effort(w) => {
                (w, self.dimension())
            }
            ScatteringObservation::Input(w) | ScatteringObservation::Output(w) => {
                (w, self.port_count())
            }
            ScatteringObservation::Outgoing(w) => (w, self.wave_ports().len()),
        };
        if weights.len() != dimension
            || weights.iter().any(|w| !w.is_finite())
            || !weights.iter().any(|w| *w != 0.0)
        {
            return Err(DrivenResonanceError::Invalid(
                "finite nonzero weights in the selected physical observation order required",
            ));
        }
        search_resonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, derivative, cancelled| {
                let factor = self.factor(omega, &mut *cancelled)?;
                let (response, slope, residual) = if derivative {
                    let result =
                        factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
                    let slope = observation.derivative(&result, cancelled)?;
                    (result.response, slope, result.relative_derivative_residual)
                } else {
                    (factor.response(input, &mut *cancelled)?, C64::ZERO, 0.0)
                };
                let value = observation.value(&response, cancelled)?;
                Ok((
                    value,
                    slope,
                    residual.max(response.physical.relative_residual),
                ))
            },
        )
        .map_err(DrivenResonanceError::from)
    }
}

impl ScatteringObservation<'_> {
    fn value(
        self,
        result: &ScatteringResponse,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, DrivenResonanceError> {
        let (weights, values) = match self {
            Self::State(w) => (w, &result.physical.state),
            Self::Effort(w) => (w, &result.physical.effort),
            Self::Input(w) => (w, &result.physical_inputs),
            Self::Output(w) => (w, &result.physical.output),
            Self::Outgoing(w) => {
                return project(w, result.waves.iter().map(|wave| wave.outgoing), cancelled);
            }
        };
        project(weights, values.iter().copied(), cancelled)
    }

    fn derivative(
        self,
        result: &ScatteringFrequencyResponse,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, DrivenResonanceError> {
        let (weights, values) = match self {
            Self::State(w) => (w, &result.state_derivative),
            Self::Effort(w) => (w, &result.effort_derivative),
            Self::Input(w) => (w, &result.physical_input_derivative),
            Self::Output(w) => (w, &result.output_derivative),
            Self::Outgoing(w) => {
                return project(
                    w,
                    result.wave_derivatives.iter().map(|wave| wave.outgoing),
                    cancelled,
                );
            }
        };
        project(weights, values.iter().copied(), cancelled)
    }
}
