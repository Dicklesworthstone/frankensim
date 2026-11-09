//! Coherent minima of the complete wave-boundary response.

use super::{ScatteringModel, ScatteringObservation, finite, poll};
use crate::harmonic::resonance::{
    AntiresonanceControls, DrivenAntiresonance, DrivenResonanceError, search_antiresonance,
};
use fs_math::c64::C64;

impl ScatteringModel {
    /// Resolve one coherent notch in the complete physically closed wave model.
    ///
    /// `input` supplies a fixed incident output-wave amplitude at each selected
    /// source port and a fixed ordinary input at every other original port.
    /// Fixed wave references close all reciprocal load feedback before each
    /// factorization. The original response and tangent retain constraints,
    /// material memory, actual conjugate inputs `(2a-y)/Z`, and outgoing waves
    /// `y-a`; phase interference is never replaced by an intensity sum.
    ///
    /// The caller isolates a strictly decreasing-to-increasing minimum of the
    /// selected squared amplitude. The shared [`search_antiresonance`] owner
    /// admits exact reflected-wave cancellation and nonzero dissipative minima
    /// without dividing by the observation. A silent reflected wave can coexist
    /// with finite internal motion, stored energy, input work and material heat.
    /// Recovery levels use the selected observation's original squared units;
    /// they are not automatically watts or an internal material quality factor.
    /// Source refusals, query caps, cancellation and original backward-residual
    /// checks remain active, and no partial minimum or band is published.
    pub fn driven_antiresonance(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, DrivenResonanceError> {
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
        search_antiresonance(
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
