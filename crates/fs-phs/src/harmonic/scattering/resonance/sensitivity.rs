//! Physical tuning of coherent extrema with the entire wave boundary closed.
//!
//! The implicit extremum/crossing equations belong to the shared resonance
//! owner. This module supplies their actual loaded observation derivatives,
//! retaining changing source coupling, storage, references and feedthrough.

use super::{ScatteringModel, ScatteringObservation};
use crate::harmonic::resonance::antiresonance_parameter;
use crate::harmonic::resonance::sensitivity::{
    ParameterObservationJet, admit_parameter_response, parameter_search_controls,
    resolve_parameter_response, tangent_effort, tangent_output,
};
use crate::harmonic::resonance::{
    AntiresonanceControls, DrivenAntiresonanceParameterResponse, DrivenResonanceError,
    DrivenResonanceParameterResponse, ResonanceControls, project,
};
use crate::harmonic::scattering::{ScatteringParameterResponse, parameter::add_parameter_source};
use crate::harmonic::{HarmonicParameterDirection, finite, parameter, poll};
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
    /// loaded observations require four constructor payload bounds plus
    /// 256n+512m+256k bytes, and four constructor work scores plus
    /// 32n(n+m)+128(m+k), where k counts selected ports. All caps are checked
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
            direction, input_derivative, reference_derivative, &mut cancelled,
        )?;
        let controls = parameter_search_controls(controls)?;
        self.admit_loaded_resonance_parameter_response()?;
        let resonance =
            self.driven_resonance(input, observation, interval_rad_s, controls, &mut cancelled)?;
        resolve_parameter_response(
            resonance,
            |omega, second_order, cancelled| {
                self.resonance_parameter_observation(
                    input, observation, direction, input_derivative, reference_derivative,
                    omega, second_order, cancelled,
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
            direction, input_derivative, reference_derivative, &mut cancelled,
        )?;
        let controls = antiresonance_parameter::parameter_controls(controls)?;
        self.admit_loaded_resonance_parameter_response()?;
        let antiresonance =
            self.driven_antiresonance(input, observation, interval_rad_s, controls, &mut cancelled)?;
        antiresonance_parameter::resolve_parameter_response(
            antiresonance,
            |omega, second_order, cancelled| {
                self.resonance_parameter_observation(
                    input, observation, direction, input_derivative, reference_derivative,
                    omega, second_order, cancelled,
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
        let model = &self.closed;
        admit_parameter_response(model)?;
        let (n, m, k) = (model.n as u128, model.m as u128, self.ports.len() as u128);
        let bytes = 4 * (8 * n * k + 256 * k + 64 * m) + 256 * n + 512 * m + 256 * k;
        let work =
            4 * (2 * n * n * (k + 1) + n * m + k * k + 32 * m) + 32 * n * (n + m) + 128 * (m + k);
        if bytes > self.budget.maximum_payload_bytes as u128
            || work > self.budget.maximum_work as u128
        {
            return Err(DrivenResonanceError::Invalid(
                "insufficient simultaneous loaded resonance parameter payload or work budget",
            ));
        }
        Ok(())
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
        let model = &self.closed;
        let factor = self.factor(omega, &mut *cancelled)?;
        let frequency =
            factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
        let parameter = factor.response_with_parameter_derivative(
            input,
            direction,
            input_derivative,
            reference_derivative,
            &mut *cancelled,
        )?;
        let mut jet = ParameterObservationJet {
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
        if !second_order {
            return Ok(jet);
        }

        // A_closed x_omegaomega = 2 i E x_omega. The authored source
        // law and boundary reference have zero frequency derivatives.
        let mut rhs = vec![C64::ZERO; model.n];
        for (target, &dx) in rhs[..model.n_diff]
            .iter_mut()
            .zip(&frequency.state_derivative)
        {
            *target = C64::new(-2.0 * dx.im, 2.0 * dx.re);
        }
        let (second, second_residual) = factor.factor.solve_rhs(&rhs, &mut *cancelled)?;
        jet.second_frequency = self.observe_higher_tangent(
            observation,
            &second,
            None,
            None,
            &frequency.effort_derivative,
            &frequency.physical_input_derivative,
            None,
            cancelled,
        )?;

        // dL*x_omega + dG*u_omega + G*(-dY*x_omega-u_omega*dZ)/Z
        // is the parameter action of the complete closed operator.
        // The same source-action owner supplies both G appearances.
        let action = parameter::generator_action(
            model,
            &direction,
            &frequency.state_derivative,
            &frequency.effort_derivative,
            cancelled,
        )?;
        rhs = action.generator;
        add_parameter_source(
            self,
            direction,
            &frequency.effort_derivative,
            &frequency.physical_input_derivative,
            &action.storage,
            None,
            reference_derivative,
            &mut rhs,
            None,
            cancelled,
        )?;
        for (target, &dx) in rhs[..model.n_diff]
            .iter_mut()
            .zip(&parameter.state_derivative)
        {
            *target = *target + C64::new(-dx.im, dx.re);
        }
        let (mixed, mixed_residual) = factor.factor.solve_rhs(&rhs, &mut *cancelled)?;
        jet.mixed = self.observe_higher_tangent(
            observation,
            &mixed,
            Some(&action.storage),
            direction.ports,
            &frequency.effort_derivative,
            &frequency.physical_input_derivative,
            reference_derivative,
            cancelled,
        )?;
        jet.relative_residual = jet
            .relative_residual
            .max(second_residual)
            .max(mixed_residual);
        Ok(jet)
    }

    // Selected original-coordinate second/mixed observation. The shared effort
    // and output maps preserve full storage, algebraic reactions and dG^T e.
    #[allow(clippy::too_many_arguments)]
    fn observe_higher_tangent(
        &self,
        observation: ScatteringObservation<'_>,
        state: &[C64],
        additive_effort: Option<&[C64]>,
        port_direction: Option<&[f64]>,
        frequency_effort: &[C64],
        frequency_input: &[C64],
        reference_derivative: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, DrivenResonanceError> {
        if let ScatteringObservation::State(weights) = observation {
            return project(weights, state.iter().copied(), cancelled);
        }
        let effort = tangent_effort(&self.closed, state, additive_effort, cancelled)?;
        if let ScatteringObservation::Effort(weights) = observation {
            return project(weights, effort, cancelled);
        }
        let output = tangent_output(
            &self.closed,
            &effort,
            port_direction,
            frequency_effort,
            cancelled,
        )?;
        match observation {
            ScatteringObservation::Output(weights) => project(weights, output, cancelled),
            ScatteringObservation::Outgoing(weights) => {
                // b=y-a, and the fixed source/parameter law has no mixed or
                // second frequency derivative. Direct feedthrough is already
                // retained in the primal and first parameter observation.
                project(
                    weights,
                    self.ports.iter().map(|&(p, _)| output[p]),
                    cancelled,
                )
            }
            ScatteringObservation::Input(weights) => {
                let mut input = vec![C64::ZERO; self.closed.m];
                for (index, &(p, z)) in self.ports.iter().enumerate() {
                    poll(cancelled)?;
                    let dz = reference_derivative.map_or(0.0, |v| v[index]);
                    input[p] = -output[p].scale(1.0 / z) - frequency_input[p].scale(dz / z);
                }
                project(weights, input, cancelled)
            }
            ScatteringObservation::State(_) | ScatteringObservation::Effort(_) => unreachable!(),
        }
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
