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
    HarmonicModel, finite, parameter, poll, resource_bounds, response::HarmonicParameterDirection,
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
                let factor = self.factor(omega, &mut *cancelled)?;
                let frequency =
                    factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
                let parameter = factor.response_with_parameter_derivative(
                    input,
                    direction,
                    input_derivative,
                    &mut *cancelled,
                )?;
                let mut jet = ParameterObservationJet {
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
                        .max(parameter.response.relative_residual)
                        .max(parameter.relative_derivative_residual),
                };
                if !second_order {
                    return Ok(jet);
                }
                let mut rhs = vec![C64::ZERO; self.n];
                for (target, &dx) in rhs[..self.n_diff]
                    .iter_mut()
                    .zip(&frequency.state_derivative)
                {
                    *target = C64::new(-2.0 * dx.im, 2.0 * dx.re);
                }
                let (second, second_residual) = factor.solve_rhs(&rhs, &mut *cancelled)?;
                let second_observation = self.observe_parameter_tangent(
                    observation,
                    &second,
                    None,
                    None,
                    &frequency.effort_derivative,
                    cancelled,
                )?;
                let action = parameter::generator_action(
                    self,
                    &direction,
                    &frequency.state_derivative,
                    &frequency.effort_derivative,
                    cancelled,
                )?;
                rhs = action.generator;
                for (target, &dx) in rhs[..self.n_diff]
                    .iter_mut()
                    .zip(&parameter.state_derivative)
                {
                    *target = *target + C64::new(-dx.im, dx.re);
                }
                let (mixed, mixed_residual) = factor.solve_rhs(&rhs, &mut *cancelled)?;
                let mixed_observation = self.observe_parameter_tangent(
                    observation,
                    &mixed,
                    Some(&action.storage),
                    direction.ports,
                    &frequency.effort_derivative,
                    cancelled,
                )?;
                jet.second_frequency = second_observation;
                jet.mixed = mixed_observation;
                jet.relative_residual = jet
                    .relative_residual
                    .max(second_residual)
                    .max(mixed_residual);
                Ok(jet)
            },
            &mut cancelled,
        )
    }

    // Observe W dx + dQ x and G^T de + dG^T e in original coordinates.
    // The additive effort/action is already computed by the shared parameter
    // owner, so there is no second implementation of its physical chain rule.
    fn observe_parameter_tangent(
        &self,
        observation: HarmonicObservation<'_>,
        state: &[C64],
        additive_effort: Option<&[C64]>,
        port_direction: Option<&[f64]>,
        original_effort: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<C64, DrivenResonanceError> {
        if let HarmonicObservation::State(weights) = observation {
            return project(weights, state.iter().copied(), cancelled);
        }
        let effort = tangent_effort(self, state, additive_effort, cancelled)?;
        match observation {
            HarmonicObservation::Effort(weights) => project(weights, effort, cancelled),
            HarmonicObservation::Output(weights) => {
                let output =
                    tangent_output(self, &effort, port_direction, original_effort, cancelled)?;
                project(weights, output, cancelled)
            }
            HarmonicObservation::State(_) => unreachable!(),
        }
    }
}

pub(in crate::harmonic) fn tangent_effort(
    model: &HarmonicModel,
    state: &[C64],
    additive_effort: Option<&[C64]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, DrivenResonanceError> {
    let mut effort = state.to_vec();
    for i in 0..model.n_diff {
        poll(cancelled)?;
        effort[i] = project(
            &model.q[i * model.n..i * model.n + model.n_diff],
            state[..model.n_diff].iter().copied(),
            cancelled,
        )?;
        if let Some(extra) = additive_effort {
            effort[i] = effort[i] + extra[i];
        }
    }
    Ok(effort)
}

pub(in crate::harmonic) fn tangent_output(
    model: &HarmonicModel,
    effort: &[C64],
    port_direction: Option<&[f64]>,
    original_effort: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, DrivenResonanceError> {
    let mut output = vec![C64::ZERO; model.m];
    for i in 0..model.n {
        poll(cancelled)?;
        for p in 0..model.m {
            output[p] = output[p] + effort[i].scale(model.g[i * model.m + p]);
            if let Some(dg) = port_direction {
                output[p] = output[p] + original_effort[i].scale(dg[i * model.m + p]);
            }
        }
    }
    Ok(output)
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
