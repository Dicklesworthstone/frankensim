//! Native physical observations for the shared bounded band survey.

use super::{ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError, survey_resonances};
use crate::harmonic::resonance::{
    DrivenResonanceError, HarmonicObservation, ResonanceSearchError, project,
};
use crate::harmonic::scattering::{ScatteringModel, ScatteringObservation};
use crate::harmonic::{HarmonicModel, finite, poll};
use fs_math::c64::C64;

impl HarmonicModel {
    /// Survey driven peaks of one physical observation under fixed input phasors.
    ///
    /// The authored positive, increasing frequencies establish the coverage;
    /// the shared survey refines detected increasing-to-decreasing brackets.
    /// Peaks narrower than the sampling coverage can remain undetected. No
    /// completeness certificate, linewidth or free-pole interpretation is added.
    /// Use `driven_resonance` on an isolated wider interval to measure bandwidth.
    ///
    /// Every query retains the original descriptor constraints, storage and
    /// physical losses. Its analytic tangent reuses that frequency's factor,
    /// and original response/tangent residuals remain enforced. Observation
    /// weights retain their physical units and coherent interference. Exact
    /// antiresonances and a silent drive are valid survey samples. The total
    /// evaluation and peak caps apply to the complete band; any refusal leaves
    /// the immutable source unchanged and publishes no partial survey.
    pub fn driven_resonances(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        frequencies_rad_s: &[f64],
        controls: ResonanceSurveyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<DrivenResonanceError>> {
        poll(&mut cancelled)
            .map_err(|e| ResonanceSearchError::Source(DrivenResonanceError::from(e)))?;
        let (weights, dimension) = match observation {
            HarmonicObservation::State(w) | HarmonicObservation::Effort(w) => (w, self.n),
            HarmonicObservation::Output(w) => (w, self.m),
        };
        validate(input, self.m, weights, dimension)?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, derivative, cancelled| {
                let factor = self.factor(omega, &mut *cancelled)?;
                let (response, slope, residual) = if derivative {
                    let result =
                        factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
                    let values = match observation {
                        HarmonicObservation::State(_) => &result.state_derivative,
                        HarmonicObservation::Effort(_) => &result.effort_derivative,
                        HarmonicObservation::Output(_) => &result.output_derivative,
                    };
                    let slope = project(weights, values.iter().copied(), cancelled)?;
                    (result.response, slope, result.relative_derivative_residual)
                } else {
                    (factor.response(input, &mut *cancelled)?, C64::ZERO, 0.0)
                };
                let values = match observation {
                    HarmonicObservation::State(_) => &response.state,
                    HarmonicObservation::Effort(_) => &response.effort,
                    HarmonicObservation::Output(_) => &response.output,
                };
                let value = project(weights, values.iter().copied(), cancelled)?;
                Ok((value, slope, residual.max(response.relative_residual)))
            },
        )
    }
}

impl ScatteringModel {
    /// Survey peaks of the complete physically wave-loaded response.
    ///
    /// Each selected source receives its fixed incident output-wave amplitude;
    /// other sources retain their fixed ordinary inputs. Reference impedances
    /// remain fixed and all boundary feedback is closed before factorization.
    /// State, effort and conjugate input/output observations keep original
    /// source order. Outgoing observations use `wave_ports` order and retain
    /// coherent incident-wave feedthrough, so their peaks can differ from
    /// internal motion peaks. Reflection minima are not reported as resonances.
    ///
    /// The shared bounded survey owns detection and refinement on the authored
    /// positive increasing frequency grid. Its coverage is not a completeness
    /// guarantee and peaks have no inferred linewidth. Every physical response
    /// and analytic tangent retains native admission and original-equation
    /// residual checks. A lossless body's loaded resonance does not imply
    /// material dissipation. All caps and cancellation apply to the whole band.
    pub fn driven_resonances(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        frequencies_rad_s: &[f64],
        controls: ResonanceSurveyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<DrivenResonanceError>> {
        poll(&mut cancelled)
            .map_err(|e| ResonanceSearchError::Source(DrivenResonanceError::from(e)))?;
        let (weights, dimension) = match observation {
            ScatteringObservation::State(w) | ScatteringObservation::Effort(w) => {
                (w, self.dimension())
            }
            ScatteringObservation::Input(w) | ScatteringObservation::Output(w) => {
                (w, self.port_count())
            }
            ScatteringObservation::Outgoing(w) => (w, self.wave_ports().len()),
        };
        validate(input, self.port_count(), weights, dimension)?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, derivative, cancelled| {
                let factor = self.factor(omega, &mut *cancelled)?;
                let (response, slope, residual) = if derivative {
                    let result =
                        factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
                    let slope = match observation {
                        ScatteringObservation::Outgoing(_) => project(
                            weights,
                            result.wave_derivatives.iter().map(|wave| wave.outgoing),
                            cancelled,
                        )?,
                        _ => {
                            let values = match observation {
                                ScatteringObservation::State(_) => &result.state_derivative,
                                ScatteringObservation::Effort(_) => &result.effort_derivative,
                                ScatteringObservation::Input(_) => {
                                    &result.physical_input_derivative
                                }
                                _ => &result.output_derivative,
                            };
                            project(weights, values.iter().copied(), cancelled)?
                        }
                    };
                    (result.response, slope, result.relative_derivative_residual)
                } else {
                    (factor.response(input, &mut *cancelled)?, C64::ZERO, 0.0)
                };
                let value = match observation {
                    ScatteringObservation::Outgoing(_) => project(
                        weights,
                        response.waves.iter().map(|wave| wave.outgoing),
                        cancelled,
                    )?,
                    _ => {
                        let values = match observation {
                            ScatteringObservation::State(_) => &response.physical.state,
                            ScatteringObservation::Effort(_) => &response.physical.effort,
                            ScatteringObservation::Input(_) => &response.physical_inputs,
                            _ => &response.physical.output,
                        };
                        project(weights, values.iter().copied(), cancelled)?
                    }
                };
                Ok((
                    value,
                    slope,
                    residual.max(response.physical.relative_residual),
                ))
            },
        )
    }
}

fn validate(
    input: &[C64],
    ports: usize,
    weights: &[f64],
    dimension: usize,
) -> Result<(), ResonanceSearchError<DrivenResonanceError>> {
    if input.len() != ports || input.iter().any(|&z| !finite(z)) {
        return Err(ResonanceSearchError::Invalid(
            "one finite excitation per original source port required",
        ));
    }
    if weights.len() != dimension
        || weights.iter().any(|w| !w.is_finite())
        || !weights.iter().any(|w| *w != 0.0)
    {
        return Err(ResonanceSearchError::Invalid(
            "finite nonzero weights in the selected physical observation order required",
        ));
    }
    Ok(())
}
