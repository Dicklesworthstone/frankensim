//! Physical tuning of coherent extrema with the entire wave boundary closed.
//!
//! The implicit extremum/crossing equations belong to the shared resonance
//! owner. This module supplies their actual loaded observation derivatives,
//! retaining changing source coupling, storage, references and feedthrough.

use super::{ScatteringModel, ScatteringObservation};
use crate::harmonic::resonance::antiresonance_parameter;
use crate::harmonic::resonance::sensitivity::{
    ParameterObservationJet, admit_parameter_response, parameter_search_controls,
    resolve_parameter_response,
};
use crate::harmonic::resonance::{
    AntiresonanceControls, DrivenAntiresonanceParameterResponse, DrivenResonanceError,
    DrivenResonanceParameterResponse, ResonanceControls, project,
};
use crate::harmonic::scattering::{
    ScatteringComplexJet, ScatteringFrequencyParameterResponse, ScatteringParameterResponse,
};
use crate::harmonic::{HarmonicParameterDirection, finite, poll};
use fs_math::c64::C64;

impl ScatteringModel {
    /// Resolve and differentiate one actual wave-loaded observed peak and its
    /// optional squared-amplitude half-height crossings, width and response Q.
    /// The result includes the TOTAL amplitude and phase derivative along the
    /// moving peak, not just a fixed-frequency tangent or free-pole motion.
    ///
    /// `direction` supplies the complete original physical dJ/dR/dQ/dG/dT.
    /// `input_derivative` has one slot per original port: incident-wave slopes
    /// at selected ports and ordinary source slopes elsewhere. The optional
    /// `reference_derivative` supplies dZ in `wave_ports()` order. All are local
    /// derivatives per one real parameter, at fixed state chart, topology and
    /// descriptor rank. Include every material and operating-point chain-rule
    /// term; no bias solve or parameter-dependent observation weight is inferred.
    ///
    /// Sources, source parameter slopes and references are constant in frequency.
    /// Changing G affects forcing, boundary feedback AND the original output.
    /// Actual boundary-input and outgoing-wave observations retain incident
    /// feedthrough; their peaks need not follow the internal motion peak. The
    /// original material loss and every memory/constraint state remain present.
    ///
    /// Uses the existing bounded search and shared implicit peak/width equations.
    /// One peak query, or three including both edges, is reserved before search.
    /// A strictly negative finite peak curvature and two nonzero crossing slopes
    /// are required. The caller must isolate a smoothly continuing branch; peak
    /// switching, flat maxima and reflection minima are not differentiated.
    ///
    /// Each query reuses one closed LU for all original response tangents. In
    /// addition to the ordinary resonance-parameter admission, simultaneous
    /// loaded observations use the complete public wave-jet admission: six
    /// constructor payload bounds plus 512(n+m+k) bytes, and six constructor
    /// work scores plus 48n(n+m)+192(m+k), where k counts selected ports.
    /// All caps are checked
    /// before search. Factors are released before another frequency is queried.
    /// No perturbed model, finite difference or extra dense matrix is used.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_resonance_with_parameter_derivative(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonanceParameterResponse, DrivenResonanceError> {
        self.validate_resonance_parameter_inputs(
            direction,
            input_derivative,
            reference_derivative,
            &mut cancelled,
        )?;
        let controls = parameter_search_controls(controls)?;
        self.admit_loaded_resonance_parameter_response()?;
        let resonance =
            self.driven_resonance(input, observation, interval_rad_s, controls, &mut cancelled)?;
        resolve_parameter_response(
            resonance,
            |omega, second_order, cancelled| {
                self.resonance_parameter_observation(
                    input,
                    observation,
                    direction,
                    input_derivative,
                    reference_derivative,
                    omega,
                    second_order,
                    cancelled,
                )
            },
            &mut cancelled,
        )
    }

    /// Resolve and differentiate one coherent notch of the fully wave-loaded
    /// physical model, including an exact simple reflected-wave zero.
    ///
    /// The returned frequency derivative follows the observed minimum, and the
    /// complex observation derivative is TOTAL along that moving minimum.
    /// Exact cancellation does not erase internal motion, material heat, stored
    /// energy or the original source and boundary-reference parameter slopes.
    /// No logarithmic amplitude, phase derivative or notch Q is formed at zero.
    ///
    /// The complete original direction, source slopes, selected-reference dZ
    /// and fixed real observation have the same meaning as
    /// [Self::driven_resonance_with_parameter_derivative]. Source/parameter laws
    /// and references remain constant in frequency. Changes in G affect the
    /// drive, reciprocal boundary feedback and original output simultaneously;
    /// actual input and outgoing-wave feedthrough remain in every derivative.
    /// All material memory, storage and algebraic reactions remain present.
    ///
    /// The caller isolates a smoothly continuing minimum with finite positive
    /// squared-amplitude curvature. Optional recovery levels are FIXED absolute
    /// squared-observation values: both moving edges need nonzero slopes away
    /// from the minimum. Flat minima, branch switching, missing edges and native
    /// refusals publish no partial sensitivity or band.
    ///
    /// One minimum factor and optional two recovery-edge factors are reserved
    /// from the original total query cap before search. The same simultaneous
    /// loaded tangent payload/work admission as the peak API applies. Each
    /// scalar jet reuses one closed LU and releases it before the next query;
    /// no perturbed model, finite difference or extra dense operator is made.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_antiresonance_with_parameter_derivative(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        interval_rad_s: (f64, f64),
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonanceParameterResponse, DrivenResonanceError> {
        self.validate_resonance_parameter_inputs(
            direction,
            input_derivative,
            reference_derivative,
            &mut cancelled,
        )?;
        let controls = antiresonance_parameter::parameter_controls(controls)?;
        self.admit_loaded_resonance_parameter_response()?;
        let antiresonance = self.driven_antiresonance(
            input,
            observation,
            interval_rad_s,
            controls,
            &mut cancelled,
        )?;
        antiresonance_parameter::resolve_parameter_response(
            antiresonance,
            |omega, second_order, cancelled| {
                self.resonance_parameter_observation(
                    input,
                    observation,
                    direction,
                    input_derivative,
                    reference_derivative,
                    omega,
                    second_order,
                    cancelled,
                )
            },
            &mut cancelled,
        )
    }

    fn validate_resonance_parameter_inputs(
        &self,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DrivenResonanceError> {
        poll(cancelled)?;
        let model = &self.closed;
        direction.validate(model, cancelled)?;
        if input_derivative.is_some_and(|v| v.len() != model.m || v.iter().any(|&z| !finite(z)))
            || reference_derivative
                .is_some_and(|v| v.len() != self.ports.len() || v.iter().any(|z| !z.is_finite()))
        {
            return Err(DrivenResonanceError::Invalid(
                "complete finite source and selected-reference derivatives required",
            ));
        }
        Ok(())
    }

    fn admit_loaded_resonance_parameter_response(&self) -> Result<(), DrivenResonanceError> {
        admit_parameter_response(&self.closed)?;
        self.admit_mixed_observation().map_err(|_| {
            DrivenResonanceError::Invalid(
                "insufficient simultaneous loaded resonance parameter payload or work budget",
            )
        })
    }

    // One original loaded frequency/parameter/second/mixed observation owner
    // serves both extrema. Keep G, dG, source, storage and dZ feedback together.
    #[allow(clippy::too_many_arguments)]
    fn resonance_parameter_observation(
        &self,
        input: &[C64],
        observation: ScatteringObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        reference_derivative: Option<&[f64]>,
        omega: f64,
        second_order: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ParameterObservationJet, DrivenResonanceError> {
        let factor = self.factor(omega, &mut *cancelled)?;
        if second_order {
            let full = factor.response_with_frequency_parameter_derivative(
                input,
                direction,
                None,
                None,
                input_derivative,
                None,
                reference_derivative,
                &mut *cancelled,
            )?;
            return observation.full_parameter_jet(&full, cancelled);
        }
        let frequency = factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
        let parameter = factor.response_with_parameter_derivative(
            input,
            direction,
            input_derivative,
            reference_derivative,
            &mut *cancelled,
        )?;
        let jet = ParameterObservationJet {
            value: observation.value(&frequency.response, cancelled)?,
            frequency: observation.derivative(&frequency, cancelled)?,
            parameter: observation.parameter_derivative(&parameter, cancelled)?,
            second_frequency: C64::ZERO,
            mixed: C64::ZERO,
            relative_residual: frequency
                .response
                .physical
                .relative_residual
                .max(frequency.relative_derivative_residual)
                .max(parameter.response.physical.relative_residual)
                .max(parameter.relative_derivative_residual),
        };
        Ok(jet)
    }
}

impl ScatteringObservation<'_> {
    fn parameter_derivative(
        self,
        result: &ScatteringParameterResponse,
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

impl ScatteringObservation<'_> {
    fn full_parameter_jet(
        self,
        result: &ScatteringFrequencyParameterResponse,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ParameterObservationJet, DrivenResonanceError> {
        let physical = &result.physical;
        let curvature = &physical.frequency_response;
        let frequency = &curvature.frequency_response;
        let parameter = &physical.parameter_response;
        let relative_residual = frequency
            .response
            .relative_residual
            .max(frequency.relative_derivative_residual)
            .max(curvature.relative_second_derivative_residual)
            .max(parameter.relative_derivative_residual)
            .max(physical.relative_mixed_derivative_residual);
        let fields: [&[C64]; 5];
        let weights = match self {
            Self::State(w) => {
                fields = [
                    &frequency.response.state,
                    &frequency.state_derivative,
                    &curvature.state_second_derivative,
                    &parameter.state_derivative,
                    &physical.state_mixed_derivative,
                ];
                w
            }
            Self::Effort(w) => {
                fields = [
                    &frequency.response.effort,
                    &frequency.effort_derivative,
                    &curvature.effort_second_derivative,
                    &parameter.effort_derivative,
                    &physical.effort_mixed_derivative,
                ];
                w
            }
            Self::Output(w) => {
                fields = [
                    &frequency.response.output,
                    &frequency.output_derivative,
                    &curvature.output_second_derivative,
                    &parameter.output_derivative,
                    &physical.output_mixed_derivative,
                ];
                w
            }
            Self::Input(w) => {
                return project_full_jet(
                    w,
                    result.physical_inputs.iter().copied(),
                    relative_residual,
                    cancelled,
                );
            }
            Self::Outgoing(w) => {
                return project_full_jet(
                    w,
                    result.waves.iter().map(|wave| wave.outgoing),
                    relative_residual,
                    cancelled,
                );
            }
        };
        Ok(ParameterObservationJet {
            value: project(weights, fields[0].iter().copied(), cancelled)?,
            frequency: project(weights, fields[1].iter().copied(), cancelled)?,
            second_frequency: project(weights, fields[2].iter().copied(), cancelled)?,
            parameter: project(weights, fields[3].iter().copied(), cancelled)?,
            mixed: project(weights, fields[4].iter().copied(), cancelled)?,
            relative_residual,
        })
    }
}

fn project_full_jet(
    weights: &[f64],
    values: impl Iterator<Item = ScatteringComplexJet> + Clone,
    relative_residual: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ParameterObservationJet, DrivenResonanceError> {
    Ok(ParameterObservationJet {
        value: project(weights, values.clone().map(|v| v.value), cancelled)?,
        frequency: project(weights, values.clone().map(|v| v.frequency), cancelled)?,
        second_frequency: project(
            weights,
            values.clone().map(|v| v.second_frequency),
            cancelled,
        )?,
        parameter: project(weights, values.clone().map(|v| v.parameter), cancelled)?,
        mixed: project(weights, values.map(|v| v.mixed), cancelled)?,
        relative_residual,
    })
}
