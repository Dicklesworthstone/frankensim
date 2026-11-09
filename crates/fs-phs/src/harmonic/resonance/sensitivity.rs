//! Physical parameter derivatives of an observed peak and its linewidth.
//!
//! At a strict peak, Re(y_omega/y)=0. Implicit differentiation uses original
//! response tangents, A x_omegaomega=2 i E x_omega and
//! A x_omega_theta=dL x_omega+i E x_theta. These linear systems share the
//! same physical factor; no perturbed model, mode matching or finite difference.

use super::{
    DrivenResonance, DrivenResonanceError, HarmonicObservation, ResonanceControls, project,
};
use crate::harmonic::{
    HarmonicModel, finite, poll, resource_bounds, response::HarmonicParameterDirection,
};
use fs_math::c64::C64;

/// Local derivatives per unit of the caller's real physical parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfPowerBandwidthDerivative {
    pub lower_frequency_derivative_rad_s: f64,
    pub upper_frequency_derivative_rad_s: f64,
    pub bandwidth_derivative_rad_s: f64,
    pub response_quality_factor_derivative: f64,
}

/// Response of the selected, smoothly continuing observed peak to one parameter.
/// These are local numerical derivatives, not enclosures or finite-change fits.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenResonanceParameterResponse {
    /// Original bracketed peak and optional squared-amplitude half-height width.
    pub resonance: DrivenResonance,
    /// d omega_peak / d theta [rad/s per parameter unit].
    pub peak_frequency_derivative_rad_s: f64,
    /// TOTAL complex observation derivative along the moving peak, including
    /// y_omega * d omega_peak/d theta. This retains source and observation phase.
    pub observation_derivative: C64,
    /// Total derivative of log |y_peak| per parameter unit.
    pub log_amplitude_derivative: f64,
    /// Total phase derivative [rad per parameter unit] in exp(-i omega t).
    pub phase_derivative: f64,
    /// Second frequency derivative of log amplitude at the evaluated
    /// peak [s^2]. Its actual signed value is returned: a strict maximum is <0.
    pub log_amplitude_curvature: f64,
    /// Re(y_omega/y) [s] at the numerically resolved peak. Its departure from
    /// zero exposes the finite peak bracket resolution; no exact-root claim.
    pub stationarity_residual: f64,
    pub half_power_derivative: Option<HalfPowerBandwidthDerivative>,
    /// Search queries plus one peak factor and, if requested, two edge factors.
    /// Never exceeds the ORIGINAL controls.maximum_evaluations.
    pub evaluations: usize,
    /// Largest native response, parameter, frequency and mixed tangent residual.
    /// This remains a backward-error diagnostic, not a near-pole forward bound.
    pub maximum_relative_residual: f64,
}

/// One original observation and its local partial derivatives. Second/mixed
/// entries are requested only at the peak; edge queries need first derivatives.
/// The source evaluator releases its factor before returning this scalar record.
pub(in crate::harmonic) struct ParameterObservationJet {
    pub value: C64,
    pub frequency: C64,
    pub parameter: C64,
    pub second_frequency: C64,
    pub mixed: C64,
    pub relative_residual: f64,
}

pub(in crate::harmonic) fn parameter_search_controls(
    mut controls: ResonanceControls,
) -> Result<ResonanceControls, DrivenResonanceError> {
    let additional = if controls.measure_half_power_bandwidth {
        3
    } else {
        1
    };
    if !(3 + additional..=4096).contains(&controls.maximum_evaluations) {
        return Err(DrivenResonanceError::Invalid(
            "evaluation capacity must include the peak search and parameter queries",
        ));
    }
    controls.maximum_evaluations -= additional;
    Ok(controls)
}

pub(in crate::harmonic) fn admit_parameter_response(
    model: &HarmonicModel,
) -> Result<(), DrivenResonanceError> {
    let (base_bytes, base_work) = resource_bounds(model.n, model.m)?;
    let bytes = base_bytes as u128 + 512 * model.n as u128 + 256 * model.m as u128;
    let work = 6 * base_work as u128;
    if bytes > model.budget.max_dense_bytes as u128 || work > model.budget.max_dense_work as u128 {
        return Err(DrivenResonanceError::Invalid(
            "insufficient simultaneous resonance parameter payload or work budget",
        ));
    }
    Ok(())
}

/// One implicit differentiation owner for ordinary and physically loaded peaks.
/// Each evaluation owns at most one factor and drops it before the next call.
pub(in crate::harmonic) fn resolve_parameter_response<C, F>(
    resonance: DrivenResonance,
    mut evaluate: F,
    cancelled: &mut C,
) -> Result<DrivenResonanceParameterResponse, DrivenResonanceError>
where
    C: FnMut() -> bool,
    F: FnMut(f64, bool, &mut C) -> Result<ParameterObservationJet, DrivenResonanceError>,
{
    poll(cancelled)?;
    let omega = resonance.peak.angular_frequency_rad_s;
    let peak = evaluate(omega, true, cancelled)?;
    let value = resonance.observation;
    let frequency_slope = resonance.observation_derivative;
    let r = ratio(frequency_slope, value)?;
    let t = ratio(peak.parameter, value)?;
    let curvature = (ratio(peak.second_frequency, value)? - r * r).re;
    let mixed_log_slope = (ratio(peak.mixed, value)? - r * t).re;
    if !curvature.is_finite() || curvature >= 0.0 || !mixed_log_slope.is_finite() {
        return Err(DrivenResonanceError::Invalid(
            "a finite strictly negative local peak curvature is required",
        ));
    }
    let peak_derivative = -mixed_log_slope / curvature;
    let observation_derivative = peak.parameter + frequency_slope.scale(peak_derivative);
    let total_log = t + r.scale(peak_derivative);
    if !peak_derivative.is_finite() || !finite(observation_derivative) || !finite(total_log) {
        return Err(DrivenResonanceError::Invalid(
            "local peak parameter derivative is not representable",
        ));
    }
    let mut residual = resonance
        .maximum_relative_residual
        .max(peak.relative_residual);
    let half_power_derivative = if let Some(band) = resonance.half_power {
        let mut crossing = |frequency: f64, lower: bool| {
            let sample = evaluate(frequency, false, &mut *cancelled)?;
            let slope = ratio(sample.frequency, sample.value)?.re;
            if (lower && slope <= 0.0) || (!lower && slope >= 0.0) {
                return Err(DrivenResonanceError::NotBracketed(
                    "half-height crossing must have a nonzero slope toward the selected peak",
                ));
            }
            let result = (total_log.re - ratio(sample.parameter, sample.value)?.re) / slope;
            if !result.is_finite() {
                return Err(DrivenResonanceError::Invalid(
                    "crossing parameter derivative is not representable",
                ));
            }
            residual = residual.max(sample.relative_residual);
            Ok(result)
        };
        let lower = crossing(band.lower.angular_frequency_rad_s, true)?;
        let upper = crossing(band.upper.angular_frequency_rad_s, false)?;
        let width = upper - lower;
        let quality =
            band.response_quality_factor * (peak_derivative / omega - width / band.bandwidth_rad_s);
        if !width.is_finite() || !quality.is_finite() {
            return Err(DrivenResonanceError::Invalid(
                "bandwidth parameter derivative is not representable",
            ));
        }
        Some(HalfPowerBandwidthDerivative {
            lower_frequency_derivative_rad_s: lower,
            upper_frequency_derivative_rad_s: upper,
            bandwidth_derivative_rad_s: width,
            response_quality_factor_derivative: quality,
        })
    } else {
        None
    };
    let evaluations = resonance.evaluations + if resonance.half_power.is_some() { 3 } else { 1 };
    poll(cancelled)?;
    Ok(DrivenResonanceParameterResponse {
        resonance,
        peak_frequency_derivative_rad_s: peak_derivative,
        observation_derivative,
        log_amplitude_derivative: total_log.re,
        phase_derivative: total_log.im,
        log_amplitude_curvature: curvature,
        stationarity_residual: r.re,
        half_power_derivative,
        evaluations,
        maximum_relative_residual: residual,
    })
}

impl HarmonicModel {
    /// Resolve an observed resonance and differentiate its frequency, complex
    /// amplitude, phase, and optional -3 dB bandwidth/Q with respect to a real
    /// physical parameter. The observation weights and input frequency law stay
    /// fixed. `input_derivative` supplies du/dtheta, constant in frequency.
    ///
    /// The direction uses the same complete original dJ/dR/dQ/dG/dT chain rule
    /// as the forced parameter response. It must include every constitutive or
    /// operating-point contribution. State chart, topology and descriptor rank
    /// are fixed. Neither a parameter-dependent observation weight nor a bias
    /// solve is inferred. This query concerns the actual observed driven peak;
    /// it does not substitute the motion of a free pole or natural frequency.
    ///
    /// A strictly negative, representable log-amplitude curvature is required.
    /// Each requested crossing must have nonzero slope toward this peak. Flat
    /// peaks, unresolved roots, exact antiresonances and native numerical/budget
    /// refusals publish no partial result. The caller must isolate one smooth
    /// peak/width branch; nearby mode or peak switching is not differentiated.
    ///
    /// Reserves one additional frequency query, or three with bandwidth, from
    /// controls.maximum_evaluations BEFORE the existing bounded search. The
    /// model additionally admits six ordinary work scores and its ordinary
    /// payload plus 512*n+256*ports bytes for simultaneous tangent arrays. One
    /// factor per evaluated frequency is reused; no extra dense matrix is made.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_resonance_with_parameter_derivative(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonanceParameterResponse, DrivenResonanceError> {
        poll(&mut cancelled)?;
        direction.validate(self, &mut cancelled)?;
        if input_derivative.is_some_and(|u| u.len() != self.m || u.iter().any(|&z| !finite(z))) {
            return Err(DrivenResonanceError::Invalid(
                "one finite parameter input derivative per original port required",
            ));
        }
        let controls = parameter_search_controls(controls)?;
        admit_parameter_response(self)?;
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
                    omega,
                    second_order,
                    cancelled,
                )
            },
            &mut cancelled,
        )
    }

    // Same original second/mixed physical observation for peaks and notches.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::harmonic) fn resonance_parameter_observation(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        omega: f64,
        second_order: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ParameterObservationJet, DrivenResonanceError> {
        let factor = self.factor(omega, &mut *cancelled)?;
        if second_order {
            let jet = factor.response_with_frequency_parameter_derivative(
                input,
                direction,
                None,
                None,
                input_derivative,
                None,
                &mut *cancelled,
            )?;
            let curvature = &jet.frequency_response;
            let frequency = &curvature.frequency_response;
            let parameter = &jet.parameter_response;
            return Ok(ParameterObservationJet {
                value: selected(
                    observation,
                    &frequency.response.state,
                    &frequency.response.effort,
                    &frequency.response.output,
                    cancelled,
                )?,
                frequency: selected(
                    observation,
                    &frequency.state_derivative,
                    &frequency.effort_derivative,
                    &frequency.output_derivative,
                    cancelled,
                )?,
                parameter: selected(
                    observation,
                    &parameter.state_derivative,
                    &parameter.effort_derivative,
                    &parameter.output_derivative,
                    cancelled,
                )?,
                second_frequency: selected(
                    observation,
                    &curvature.state_second_derivative,
                    &curvature.effort_second_derivative,
                    &curvature.output_second_derivative,
                    cancelled,
                )?,
                mixed: selected(
                    observation,
                    &jet.state_mixed_derivative,
                    &jet.effort_mixed_derivative,
                    &jet.output_mixed_derivative,
                    cancelled,
                )?,
                relative_residual: frequency
                    .response
                    .relative_residual
                    .max(frequency.relative_derivative_residual)
                    .max(parameter.relative_derivative_residual)
                    .max(curvature.relative_second_derivative_residual)
                    .max(jet.relative_mixed_derivative_residual),
            });
        }
        // Width edges need only first tangents. Reuse their original
        // response rather than another primal solve at the same factor.
        let frequency = factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
        let parameter = factor.response_parameter_from_primal(
            input,
            direction,
            input_derivative,
            frequency.response.clone(),
            &mut *cancelled,
        )?;
        Ok(ParameterObservationJet {
            value: selected(
                observation,
                &frequency.response.state,
                &frequency.response.effort,
                &frequency.response.output,
                cancelled,
            )?,
            frequency: selected(
                observation,
                &frequency.state_derivative,
                &frequency.effort_derivative,
                &frequency.output_derivative,
                cancelled,
            )?,
            parameter: selected(
                observation,
                &parameter.state_derivative,
                &parameter.effort_derivative,
                &parameter.output_derivative,
                cancelled,
            )?,
            second_frequency: C64::ZERO,
            mixed: C64::ZERO,
            relative_residual: frequency
                .response
                .relative_residual
                .max(frequency.relative_derivative_residual)
                .max(parameter.relative_derivative_residual),
        })
    }
}

fn selected(
    observation: HarmonicObservation<'_>,
    state: &[C64],
    effort: &[C64],
    output: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<C64, DrivenResonanceError> {
    let (weights, values) = match observation {
        HarmonicObservation::State(w) => (w, state),
        HarmonicObservation::Effort(w) => (w, effort),
        HarmonicObservation::Output(w) => (w, output),
    };
    project(weights, values.iter().copied(), cancelled)
}

fn ratio(numerator: C64, denominator: C64) -> Result<C64, DrivenResonanceError> {
    if denominator == C64::ZERO {
        return Err(DrivenResonanceError::NotBracketed(
            "a silent observation has no local logarithmic derivative",
        ));
    }
    let result = numerator / denominator;
    if finite(result) {
        Ok(result)
    } else {
        Err(DrivenResonanceError::Invalid(
            "observed logarithmic derivative is not representable",
        ))
    }
}
