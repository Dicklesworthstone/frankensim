//! Exact physical tuning of a coherent notch, including a simple response zero.
//! Stationarity is Re(conj(H) H_omega)=0. Its parameter derivative never
//! divides by H; both nonzero minima and simple destructive zeros are smooth.

use super::{
    AntiresonanceControls, DrivenAntiresonance, DrivenResonanceError, HarmonicObservation,
};
use super::sensitivity::{ParameterObservationJet, admit_parameter_response};
use crate::harmonic::{HarmonicModel, HarmonicParameterDirection, finite, poll};
use fs_math::{c64::C64, dd::Dd};

/// Recovery-edge motion at a FIXED absolute squared-observation level.
/// Units are rad/s per unit of the caller's original real parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AntiresonanceRecoveryDerivative {
    pub lower_frequency_derivative_rad_s: f64,
    pub upper_frequency_derivative_rad_s: f64,
    pub bandwidth_derivative_rad_s: f64,
}

/// Local parameter motion of one smoothly continuing original coherent notch.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenAntiresonanceParameterResponse {
    pub antiresonance: DrivenAntiresonance,
    pub minimum_frequency_derivative_rad_s: f64,
    /// TOTAL complex derivative along the moving minimum: H_p + H_omega omega_p.
    /// This is meaningful at H=0, where amplitude and phase derivatives need
    /// not exist. No phase, logarithmic amplitude or notch Q is manufactured.
    pub observation_derivative: C64,
    /// Positive numerical scale a applied uniformly to response/frequency jets.
    /// It carries the observation's units; its choice is only conditioning.
    pub curvature_amplitude_scale: f64,
    /// 2*(|H_omega/a|² + Re(conj(H/a) H_omegaomega/a)) [s²].
    /// Multiplying by a² recovers physical squared-amplitude curvature when
    /// representable. Scaling avoids squaring extreme source amplitudes.
    pub scaled_power_curvature_s2: f64,
    /// Re(conj(H/a) H_omega/a) [s] at the actual numerical minimum.
    /// This exposes finite root resolution and is not an enclosure.
    pub scaled_stationarity_residual_s: f64,
    pub recovery_derivative: Option<AntiresonanceRecoveryDerivative>,
    /// Search queries plus one minimum factor and optional two edge factors.
    /// The original caller's total maximum_evaluations is never exceeded.
    pub evaluations: usize,
    pub maximum_relative_residual: f64,
}

pub(in crate::harmonic) fn parameter_controls(
    mut controls: AntiresonanceControls,
) -> Result<AntiresonanceControls, DrivenResonanceError> {
    let additional = if controls.recovery_power.is_some() { 3 } else { 1 };
    if !(3 + additional..=4096).contains(&controls.maximum_evaluations) {
        return Err(DrivenResonanceError::Invalid(
            "evaluation capacity must include notch search and parameter queries",
        ));
    }
    controls.maximum_evaluations -= additional;
    Ok(controls)
}

pub(in crate::harmonic) fn resolve_parameter_response<C, F>(
    antiresonance: DrivenAntiresonance,
    mut evaluate: F,
    cancelled: &mut C,
) -> Result<DrivenAntiresonanceParameterResponse, DrivenResonanceError>
where
    C: FnMut() -> bool,
    F: FnMut(f64, bool, &mut C) -> Result<ParameterObservationJet, DrivenResonanceError>,
{
    poll(cancelled)?;
    let jet = evaluate(antiresonance.minimum.angular_frequency_rad_s, true, cancelled)?;
    let h = antiresonance.observation;
    let hw = antiresonance.observation_derivative;
    if [h, hw, jet.parameter, jet.second_frequency, jet.mixed]
        .iter()
        .any(|&value| !finite(value))
    {
        return Err(DrivenResonanceError::Invalid("finite original notch jets required"));
    }
    // At an exact cancellation these jets have identically zero weights.
    // Including them in a numerical scale can erase a small, representable
    // frequency shift even though they contribute nothing to it.
    let (second_frequency, mixed) = if h == C64::ZERO {
        (C64::ZERO, C64::ZERO)
    } else {
        (jet.second_frequency, jet.mixed)
    };
    let scale = conditioning_scale(maximum(h).max(maximum(hw)).max(maximum(second_frequency)));
    if !scale.is_finite() || scale == 0.0 {
        return Err(DrivenResonanceError::Invalid("resolved nonflat notch jets required"));
    }
    let h = divided(h, scale)?;
    let hw = divided(hw, scale)?;
    let hww = divided(second_frequency, scale)?;
    let curvature = dot(hw, hw)? + dot(h, hww)?;
    let scaled_power_curvature_s2 = 2.0 * checked(curvature)?;
    if !scaled_power_curvature_s2.is_finite() || scaled_power_curvature_s2 <= 0.0 {
        return Err(DrivenResonanceError::Invalid(
            "strictly positive resolved notch power curvature required",
        ));
    }
    let parameter_scale = conditioning_scale(maximum(jet.parameter).max(maximum(mixed)));
    if !parameter_scale.is_finite() {
        return Err(DrivenResonanceError::Invalid("finite original parameter jets required"));
    }
    let minimum_derivative = if parameter_scale == 0.0 {
        0.0
    } else {
        let hp = divided(jet.parameter, parameter_scale)?;
        let hwp = divided(mixed, parameter_scale)?;
        // g_p = Re(conj(H_p) H_omega + conj(H) H_omega_p).
        // g_omega = |H_omega|² + Re(conj(H) H_omegaomega).
        scaled_ratio(-(dot(hp, hw)? + dot(h, hwp)?), curvature, parameter_scale, scale)?
    };
    let observation_derivative = C64::new(
        checked(Dd::from_f64(jet.parameter.re)
            + Dd::from_f64(antiresonance.observation_derivative.re)
                * Dd::from_f64(minimum_derivative))?,
        checked(Dd::from_f64(jet.parameter.im)
            + Dd::from_f64(antiresonance.observation_derivative.im)
                * Dd::from_f64(minimum_derivative))?,
    );
    let mut residual = antiresonance.maximum_relative_residual.max(jet.relative_residual);
    let recovery_derivative = if let Some(band) = antiresonance.recovery {
        let mut edge = |omega: f64, lower: bool| {
            let jet = evaluate(omega, false, &mut *cancelled)?;
            if [jet.value, jet.frequency, jet.parameter]
                .iter()
                .any(|&value| !finite(value))
            {
                return Err(DrivenResonanceError::Invalid("finite original recovery jets required"));
            }
            let a = conditioning_scale(maximum(jet.value));
            let b = conditioning_scale(maximum(jet.frequency));
            if a == 0.0 || b == 0.0 {
                return Err(DrivenResonanceError::NotBracketed(
                    "recovery edge requires nonzero coherent response and frequency slope",
                ));
            }
            let value = divided(jet.value, a)?;
            let slope = dot(value, divided(jet.frequency, b)?)?;
            let slope_sign = checked(slope)?;
            if (lower && slope_sign >= 0.0) || (!lower && slope_sign <= 0.0) {
                return Err(DrivenResonanceError::NotBracketed(
                    "recovery edge must slope away from the selected minimum",
                ));
            }
            let parameter_scale = conditioning_scale(maximum(jet.parameter));
            let derivative = if parameter_scale == 0.0 {
                0.0
            } else {
                let parameter = divided(jet.parameter, parameter_scale)?;
                // The authored ABSOLUTE recovery level is fixed. Unlike peak
                // half-height sensitivity, no moving-height term belongs here.
                scaled_ratio(-dot(value, parameter)?, slope, parameter_scale, b)?
            };
            residual = residual.max(jet.relative_residual);
            Ok(derivative)
        };
        let lower = edge(band.lower.angular_frequency_rad_s, true)?;
        let upper = edge(band.upper.angular_frequency_rad_s, false)?;
        let width = checked(Dd::from_f64(upper) - Dd::from_f64(lower))?;
        Some(AntiresonanceRecoveryDerivative {
            lower_frequency_derivative_rad_s: lower,
            upper_frequency_derivative_rad_s: upper,
            bandwidth_derivative_rad_s: width,
        })
    } else {
        None
    };
    let evaluations = antiresonance.evaluations
        + if antiresonance.recovery.is_some() { 3 } else { 1 };
    poll(cancelled)?;
    Ok(DrivenAntiresonanceParameterResponse {
        antiresonance,
        minimum_frequency_derivative_rad_s: minimum_derivative,
        observation_derivative,
        curvature_amplitude_scale: scale,
        scaled_power_curvature_s2,
        scaled_stationarity_residual_s: checked(dot(h, hw)?)?,
        recovery_derivative,
        evaluations,
        maximum_relative_residual: residual,
    })
}

impl HarmonicModel {
    /// Differentiate an original coherent notch and both optional recovery
    /// edges with the complete physical dJ/dR/dQ/dG/dT direction.
    ///
    /// The source and its optional parameter derivative are fixed in frequency;
    /// original receiver weights, topology, state chart and descriptor rank
    /// stay fixed. Every material, source and operating-point contribution
    /// belongs in the supplied direction. No bias solve or moving receiver is
    /// inferred. Simple exact cancellations are admitted without dividing by H.
    ///
    /// A strictly positive local squared-amplitude curvature is required.
    /// Flat or higher-order zeros, zero-slope recovery edges and unresolved
    /// numerical jet products refuse instead of publishing a false zero. The
    /// caller must isolate one smoothly continuing minimum; mode/notch branch
    /// switching is not differentiated. Recovery_power stays an ABSOLUTE fixed
    /// squared-observation level through the parameter change.
    ///
    /// Reserves one parameter query, or three with recovery edges, before the
    /// bounded search. Original second/mixed tangent solves share one factor
    /// per frequency with the existing peak sensitivity owner. Its six-work
    /// score and simultaneous payload admission apply here unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn driven_antiresonance_with_parameter_derivative(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        direction: HarmonicParameterDirection<'_>,
        input_derivative: Option<&[C64]>,
        interval_rad_s: (f64, f64),
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonanceParameterResponse, DrivenResonanceError> {
        poll(&mut cancelled)?;
        direction.validate(self, &mut cancelled)?;
        if input_derivative.is_some_and(|u| u.len() != self.m || u.iter().any(|&z| !finite(z))) {
            return Err(DrivenResonanceError::Invalid(
                "one finite parameter input derivative per original port required",
            ));
        }
        let controls = parameter_controls(controls)?;
        admit_parameter_response(self)?;
        let antiresonance =
            self.driven_antiresonance(input, observation, interval_rad_s, controls, &mut cancelled)?;
        resolve_parameter_response(
            antiresonance,
            |omega, second_order, cancelled| {
                self.resonance_parameter_observation(
                    input, observation, direction, input_derivative,
                    omega, second_order, cancelled,
                )
            },
            &mut cancelled,
        )
    }
}

fn maximum(value: C64) -> f64 {
    value.re.abs().max(value.im.abs())
}
// A binary scale preserves every admitted normal component exactly. An
// arbitrary maximum would round components before their compensated products.
fn conditioning_scale(value: f64) -> f64 {
    if !value.is_finite() || value <= 0.0 {
        return value;
    }
    let bits = value.to_bits();
    let exponent = bits & (0x7ff_u64 << 52);
    if exponent == 0 {
        f64::from_bits(1_u64 << (63 - bits.leading_zeros()))
    } else {
        f64::from_bits(exponent)
    }
}

fn divided(value: C64, scale: f64) -> Result<C64, DrivenResonanceError> {
    if !finite(value) || !scale.is_finite() || scale <= 0.0 {
        return Err(DrivenResonanceError::Invalid("finite nonzero jet scale required"));
    }
    let result = C64::new(value.re / scale, value.im / scale);
    if finite(result)
        && (value.re == 0.0 || result.re.is_normal())
        && (value.im == 0.0 || result.im.is_normal())
    {
        Ok(result)
    } else {
        Err(DrivenResonanceError::Invalid("scaled notch jet is not representable"))
    }
}
fn dot(a: C64, b: C64) -> Result<Dd, DrivenResonanceError> {
    // DD preserves cancellation of representable products, but does not extend
    // their exponent range. Reserve all 106 product bits above the smallest
    // subnormal: losing even the low product word can erase a physical shift
    // when opposite terms cancel and the final parameter scale restores it.
    // Exact cancellation between resolved products remains a valid zero.
    const MIN_PRODUCT_MAGNITUDE: f64 = f64::MIN_POSITIVE / (0.5 * f64::EPSILON);
    let product = |x: f64, y: f64| {
        let value = Dd::from_f64(x) * Dd::from_f64(y);
        if x != 0.0 && y != 0.0
            && (!value.hi.is_finite() || value.hi.abs() < MIN_PRODUCT_MAGNITUDE)
        {
            Err(DrivenResonanceError::Invalid(
                "scaled notch jet product is not representable",
            ))
        } else {
            Ok(value)
        }
    };
    Ok(product(a.re, b.re)? + product(a.im, b.im)?)
}
fn checked(value: Dd) -> Result<f64, DrivenResonanceError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(DrivenResonanceError::Invalid("notch parameter arithmetic is not representable"))
    }
}
fn scaled_ratio(
    numerator: Dd,
    denominator: Dd,
    numerator_scale: f64,
    denominator_scale: f64,
) -> Result<f64, DrivenResonanceError> {
    if numerator == Dd::ZERO {
        return Ok(0.0);
    }
    let a = Dd::from_f64(numerator_scale);
    let b = Dd::from_f64(denominator_scale);
    // Retain representable answers when a premature scale quotient overflows
    // or underflows. No exponent, frequency, damping or derivative is clipped.
    for order in 0..4 {
        let first = match order {
            0 | 1 => numerator / denominator,
            2 => numerator * a,
            _ => numerator / b,
        };
        if !first.hi.is_normal() || !first.lo.is_finite() {
            continue;
        }
        let second = match order {
            0 | 3 => first * a,
            _ => first / b,
        };
        if !second.hi.is_normal() || !second.lo.is_finite() {
            continue;
        }
        let value = match order {
            0 => second / b,
            1 => second * a,
            _ => second / denominator,
        };
        if let Ok(value) = checked(value) {
            if value != 0.0 {
                return Ok(value);
            }
        }
    }
    Err(DrivenResonanceError::Invalid(
        "nonzero notch parameter derivative is unrepresentable",
    ))
}

#[cfg(test)]
mod tests;
