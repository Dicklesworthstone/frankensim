//! Analytic parameter sensitivity of an original observed peak or notch.
//!
//! Implicit differentiation uses g = Re(conj(H) H_w), not H_w/H.
//! A simple exact transfer zero therefore retains its finite notch-frequency
//! derivative. The caller supplies the complete source/receiver/material jet;
//! no frequency or parameter perturbation, pole fit, or mode matching occurs.

use super::query::refine_half_power_bandwidth_samples_with_cancel;
use super::{
    HalfPowerBandwidth, HarmonicTransferSample, ResonanceControls, ResonanceError,
    TransferExtremum, TransferExtremumKind, refine_transfer_extremum_with_cancel,
};
use fs_math::{c64::C64, dd::Dd};

/// One original observation and its partial derivatives at fixed parameter/frequency.
/// All entries use the same source, receiver, coordinates, and exp(-i*w*t) convention.
/// Frequency-dependent source/receiver/load laws must include every chain-rule term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicTransferParameterJet {
    /// H and its complete partial derivative H_w, per rad/s.
    pub sample: HarmonicTransferSample,
    /// H_ww, per (rad/s)^2.
    pub angular_second_derivative: C64,
    /// H_theta at fixed frequency, per unit real physical parameter.
    pub parameter_derivative: C64,
    /// H_wtheta, per rad/s per unit parameter.
    pub angular_parameter_derivative: C64,
}

/// Local response of the same smoothly continuing observed extremum.
/// These are numerical derivatives, not finite-change predictions or error enclosures.
#[derive(Clone, Copy, Debug)]
pub struct TransferExtremumParameterResponse {
    /// Original bracketed physical peak or notch.
    pub extremum: TransferExtremum,
    /// d w_extremum/d theta, in rad/s per parameter unit.
    pub angular_frequency_derivative_rad_s: f64,
    /// TOTAL complex derivative H_theta + H_w * d w_extremum/d theta.
    pub amplitude_derivative: C64,
    /// Total derivative of log |H|; absent at an exact transfer zero.
    pub log_amplitude_derivative: Option<f64>,
    /// Total phase derivative in radians; absent at an exact transfer zero.
    pub phase_derivative: Option<f64>,
    /// Signed g_w divided by the sum of absolute real product terms forming it.
    /// Dimensionless, in [-1,1]; its sign follows the strict peak/notch curvature
    /// when representable. This diagnostic can round to zero for subnormal
    /// curvature/envelope ratios; orientation validation retains the underlying
    /// nonzero signed curvature without a diagnostic floor.
    /// A small magnitude exposes cancellation in the local curvature, not a
    /// forward-error certificate. No curvature floor is introduced.
    pub normalized_curvature: f64,
    /// g divided by the sum of absolute real product terms forming it, or zero
    /// when all terms vanish. This exposes finite search resolution without
    /// squaring the original response. It is not an exact stationarity claim.
    pub normalized_stationarity_residual: f64,
    /// Original first-derivative search calls plus the one final jet call.
    pub evaluations: usize,
}

/// Actual half-height width and its local parameter derivatives.
/// This differentiates the selected observable's width, not a free-pole decay Q.
#[derive(Clone, Copy, Debug)]
pub struct HalfPowerBandwidthParameterResponse {
    /// Actual original-observation half-power crossings and width.
    pub bandwidth: HalfPowerBandwidth,
    /// d w_lower/d theta, in rad/s per parameter unit.
    pub lower_angular_frequency_derivative_rad_s: f64,
    /// d w_upper/d theta, in rad/s per parameter unit.
    pub upper_angular_frequency_derivative_rad_s: f64,
    /// d (w_upper - w_lower)/d theta, in rad/s per parameter unit.
    pub angular_bandwidth_derivative_rad_s: f64,
    /// d (w_peak / bandwidth)/d theta, per parameter unit.
    pub response_quality_factor_derivative: f64,
    /// Flank-search calls plus two edge-jet calls; excludes the prior peak query.
    pub evaluations: usize,
}

/// Refine and differentiate one observed resonance peak OR antiresonance notch.
///
/// The first callback supplies H,H_w during the existing bounded sign search.
/// Exactly one full jet is requested at the returned frequency. It must reproduce
/// that same sample to 32 machine epsilons in each complex observation's largest
/// quadrature, without a unity floor. This admits ordinary summation roundoff;
/// an exact zero versus a nonzero observation still refuses. The accepted full
/// jet becomes the returned sample. Changing the physical model, drive or
/// receiver beyond that allowance refuses. Both callbacks retain their native work,
/// domain, solve-residual and internal cancellation obligations.
///
/// One call is reserved from controls.maximum_evaluations BEFORE searching.
/// Returned evaluations include it. No full parameter/curvature solve is done
/// during frequency refinement. Cancellation, source errors, mismatched samples,
/// zero/wrong-signed curvature or unrepresentable derivatives publish no result.
///
/// With g=Re(conj(H)H_w), g_w=|H_w|^2+Re(conj(H)H_ww) and
/// g_theta=Re(conj(H_theta)H_w+conj(H)H_wtheta), the frequency derivative is
/// -g_theta/g_w. This remains valid at a simple exact notch. Log-amplitude and
/// phase derivatives are then absent, while the complex derivative is retained.
/// The caller isolates one smooth branch; peak switching and topology changes
/// are not differentiable through this query.
pub fn refine_transfer_extremum_with_parameter_derivative_with_cancel<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    mut evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<TransferExtremumParameterResponse, ResonanceError<E>> {
    poll(&mut cancelled)?;
    let search_controls = reserve(controls, 1)?;
    let extremum = refine_transfer_extremum_with_cancel(
        bracket_rad_s,
        search_controls,
        &mut cancelled,
        evaluate,
    )
    .map_err(|error| original_budget(error, controls.maximum_evaluations))?;
    poll(&mut cancelled)?;
    let jet = evaluate_jet(extremum.angular_frequency_rad_s).map_err(ResonanceError::Evaluation)?;
    poll(&mut cancelled)?;
    let response = resolve(extremum, jet)?;
    poll(&mut cancelled)?;
    Ok(response)
}

/// Differentiate the actual half-power crossings of an already differentiated peak.
///
/// Both callbacks describe the SAME immutable transfer/parameter as peak.
/// The lower and upper brackets independently isolate the rising/falling
/// half-height crossings; the existing search supplies their actual line shape.
/// Two edge-jet calls are reserved from this query's evaluation cap. The peak's
/// earlier search/evaluation count is not charged again. Only H,H_w,H_theta
/// are consumed at an edge; unused second/mixed entries may be zero.
/// Each edge jet must reproduce its original crossing's H and H_w under the
/// same relative consistency check used by the extremum query. Retaining those
/// two samples requires neither extra source evaluations nor a search history.
///
/// The crossing derivative includes motion of the peak height. A source gain
/// change therefore cannot masquerade as a linewidth change. Zero or incorrectly
/// oriented crossing slopes refuse. This is a local smooth-branch derivative,
/// not a Lorentzian assumption, damping fit, or certificate of crossing uniqueness.
#[allow(clippy::too_many_arguments)]
pub fn refine_half_power_bandwidth_with_parameter_derivative_with_cancel<E>(
    peak: &TransferExtremumParameterResponse,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    mut evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<E>> {
    poll(&mut cancelled)?;
    if !finite(peak.amplitude_derivative) || !peak.angular_frequency_derivative_rad_s.is_finite() {
        return Err(ResonanceError::Invalid(
            "finite original peak parameter response required",
        ));
    }
    let search_controls = reserve(controls, 2)?;
    let (bandwidth, samples) = refine_half_power_bandwidth_samples_with_cancel(
        &peak.extremum,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        search_controls,
        &mut cancelled,
        evaluate,
    )
    .map_err(|error| original_budget(error, controls.maximum_evaluations))?;
    let mut derivatives = [0.0; 2];
    for (index, frequency) in [
        bandwidth.lower_angular_frequency_rad_s,
        bandwidth.upper_angular_frequency_rad_s,
    ]
    .into_iter()
    .enumerate()
    {
        poll(&mut cancelled)?;
        let jet = evaluate_jet(frequency).map_err(ResonanceError::Evaluation)?;
        poll(&mut cancelled)?;
        if !finite(jet.sample.amplitude)
            || !finite(jet.sample.angular_derivative)
            || !finite(jet.parameter_derivative)
        {
            return Err(ResonanceError::Invalid(
                "finite original half-power parameter jet required",
            ));
        }
        if !same_observation(jet.sample.amplitude, samples[index].amplitude)
            || !same_observation(
                jet.sample.angular_derivative,
                samples[index].angular_derivative,
            )
        {
            return Err(ResonanceError::Invalid(
                "half-power parameter jet must reproduce the same original crossing transfer and frequency derivative",
            ));
        }
        let h = jet.sample.amplitude;
        let w = jet.sample.angular_derivative;
        let t = jet.parameter_derivative;
        let p = peak.extremum.sample.amplitude;
        let pt = peak.amplitude_derivative;
        let slope = Products::sum(&[(h.re, w.re, 0), (h.im, w.im, 0)]);
        if (index == 0 && slope.fraction.hi <= 0.0) || (index == 1 && slope.fraction.hi >= 0.0) {
            return Err(ResonanceError::InvalidBracket(
                "nonzero half-power slope toward the same peak required",
            ));
        }
        let numerator = Products::sum(&[
            (p.re, pt.re, -1),
            (p.im, pt.im, -1),
            (-h.re, t.re, 0),
            (-h.im, t.im, 0),
        ]);
        derivatives[index] = numerator.ratio(slope)?;
    }
    let width_derivative =
        Products::sum(&[(derivatives[1], 1.0, 0), (-derivatives[0], 1.0, 0)]).value()?;
    let width = bandwidth.angular_bandwidth_rad_s;
    let quality_derivative = Products::sum(&[
        (peak.angular_frequency_derivative_rad_s, width, 0),
        (-peak.extremum.angular_frequency_rad_s, width_derivative, 0),
    ])
    .ratio(Products::sum(&[(width, width, 0)]))?;
    poll(&mut cancelled)?;
    Ok(HalfPowerBandwidthParameterResponse {
        bandwidth,
        lower_angular_frequency_derivative_rad_s: derivatives[0],
        upper_angular_frequency_derivative_rad_s: derivatives[1],
        angular_bandwidth_derivative_rad_s: width_derivative,
        response_quality_factor_derivative: quality_derivative,
        evaluations: bandwidth.evaluations + 2,
    })
}

fn resolve<E>(
    mut extremum: TransferExtremum,
    jet: HarmonicTransferParameterJet,
) -> Result<TransferExtremumParameterResponse, ResonanceError<E>> {
    if [
        jet.sample.amplitude,
        jet.sample.angular_derivative,
        jet.angular_second_derivative,
        jet.parameter_derivative,
        jet.angular_parameter_derivative,
    ]
    .into_iter()
    .any(|z| !finite(z))
    {
        return Err(ResonanceError::Invalid(
            "finite complete original parameter jet required",
        ));
    }
    if !same_observation(jet.sample.amplitude, extremum.sample.amplitude)
        || !same_observation(
            jet.sample.angular_derivative,
            extremum.sample.angular_derivative,
        )
    {
        return Err(ResonanceError::Invalid(
            "parameter jet must reproduce the same original transfer and frequency derivative",
        ));
    }
    extremum.sample = jet.sample;
    let h = jet.sample.amplitude;
    let w = jet.sample.angular_derivative;
    let ww = jet.angular_second_derivative;
    let t = jet.parameter_derivative;
    let wt = jet.angular_parameter_derivative;
    let curvature = Products::sum(&[
        (w.re, w.re, 0),
        (w.im, w.im, 0),
        (h.re, ww.re, 0),
        (h.im, ww.im, 0),
    ]);
    let sign_matches = match extremum.kind {
        TransferExtremumKind::Maximum => curvature.fraction.hi < 0.0,
        TransferExtremumKind::Minimum => curvature.fraction.hi > 0.0,
    };
    if !sign_matches {
        return Err(ResonanceError::Invalid(
            "nonzero local curvature with the original extremum orientation required",
        ));
    }
    let parameter_slope = Products::sum(&[
        (t.re, w.re, 0),
        (t.im, w.im, 0),
        (h.re, wt.re, 0),
        (h.im, wt.im, 0),
    ]);
    let frequency_derivative = -parameter_slope.ratio(curvature)?;
    let amplitude_derivative = C64::new(
        Products::sum(&[(t.re, 1.0, 0), (w.re, frequency_derivative, 0)]).value()?,
        Products::sum(&[(t.im, 1.0, 0), (w.im, frequency_derivative, 0)]).value()?,
    );
    let (log_amplitude_derivative, phase_derivative) = if h == C64::ZERO {
        (None, None)
    } else {
        let ratio = amplitude_derivative / h;
        if !finite(ratio) {
            return Err(ResonanceError::Invalid(
                "log-amplitude or phase parameter derivative is not representable",
            ));
        }
        (Some(ratio.re), Some(ratio.im))
    };
    Ok(TransferExtremumParameterResponse {
        extremum,
        angular_frequency_derivative_rad_s: frequency_derivative,
        amplitude_derivative,
        log_amplitude_derivative,
        phase_derivative,
        normalized_curvature: curvature.normalized(),
        normalized_stationarity_residual: Products::sum(&[(h.re, w.re, 0), (h.im, w.im, 0)])
            .normalized(),
        evaluations: extremum.evaluations + 1,
    })
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn same_observation(a: C64, b: C64) -> bool {
    let scale = a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs());
    scale == 0.0
        || ((a.re / scale - b.re / scale).abs() <= 32.0 * f64::EPSILON
            && (a.im / scale - b.im / scale).abs() <= 32.0 * f64::EPSILON)
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResonanceError<E>> {
    if cancelled() {
        Err(ResonanceError::Cancelled)
    } else {
        Ok(())
    }
}

fn reserve<E>(
    mut controls: ResonanceControls,
    count: usize,
) -> Result<ResonanceControls, ResonanceError<E>> {
    if controls.maximum_evaluations < 3 + count {
        return Err(ResonanceError::Invalid(
            "evaluation capacity must include the search and original parameter jets",
        ));
    }
    controls.maximum_evaluations -= count;
    Ok(controls)
}

fn original_budget<E>(error: ResonanceError<E>, maximum_evaluations: usize) -> ResonanceError<E> {
    match error {
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations,
        },
        other => other,
    }
}

// Short signed contractions are accumulated as binary-normalized double-double
// products before their ratio is taken. Thus common response/parameter/time-unit
// changes do not overflow or erase the implicit equation by squaring H or H_w.
// Constant scratch, no allocation. Cancellation beyond the retained precision
// remains a conditioning limitation, not a relative-accuracy certificate.
#[derive(Clone, Copy)]
pub(super) struct Products {
    fraction: Dd,
    envelope: Dd,
    exponent: i32,
}

impl Products {
    pub(super) fn sum(terms: &[(f64, f64, i32)]) -> Self {
        let exponent = terms
            .iter()
            .filter(|(a, b, _)| *a != 0.0 && *b != 0.0)
            .map(|&(a, b, shift)| parts(a).1 + parts(b).1 + shift)
            .max()
            .unwrap_or(0);
        let mut fraction = Dd::ZERO;
        let mut envelope = Dd::ZERO;
        for &(a, b, shift) in terms {
            if a == 0.0 || b == 0.0 {
                continue;
            }
            let (a, ae) = parts(a);
            let (b, be) = parts(b);
            let product = Dd::from_f64(a) * Dd::from_f64(b);
            let shift = ae + be + shift - exponent;
            let value = Dd::from_pair(
                binary_scale(product.hi, shift),
                binary_scale(product.lo, shift),
            );
            fraction = fraction + value;
            envelope = envelope + value.abs();
        }
        Self {
            fraction,
            envelope,
            exponent,
        }
    }

    fn normalized(self) -> f64 {
        if self.envelope.hi == 0.0 {
            0.0
        } else {
            (self.fraction / self.envelope).to_f64()
        }
    }

    pub(super) fn ratio<E>(self, denominator: Self) -> Result<f64, ResonanceError<E>> {
        if denominator.fraction.hi == 0.0 {
            return Err(ResonanceError::Invalid(
                "zero implicit parameter derivative denominator",
            ));
        }
        if self.fraction.hi == 0.0 {
            return Ok(0.0);
        }
        // Cancellation can leave either contraction's fraction subnormal even
        // though its separately stored power of two makes the physical ratio
        // ordinary. Normalize both DD fractions BEFORE dividing; otherwise an
        // intermediate overflow/underflow would discard that finite ratio.
        let numerator_shift = parts(self.fraction.hi).1;
        let denominator_shift = parts(denominator.fraction.hi).1;
        let numerator = Dd::from_pair(
            binary_scale(self.fraction.hi, -numerator_shift),
            binary_scale(self.fraction.lo, -numerator_shift),
        );
        let denominator_fraction = Dd::from_pair(
            binary_scale(denominator.fraction.hi, -denominator_shift),
            binary_scale(denominator.fraction.lo, -denominator_shift),
        );
        let fraction = (numerator / denominator_fraction).to_f64();
        representable(
            fraction,
            self.exponent - denominator.exponent + numerator_shift - denominator_shift,
        )
    }

    pub(super) fn value<E>(self) -> Result<f64, ResonanceError<E>> {
        representable(self.fraction.to_f64(), self.exponent)
    }
}

fn representable<E>(fraction: f64, exponent: i32) -> Result<f64, ResonanceError<E>> {
    let value = binary_scale(fraction, exponent);
    if !value.is_finite() || (fraction != 0.0 && value == 0.0) {
        Err(ResonanceError::Invalid(
            "physical parameter derivative is not representable",
        ))
    } else {
        Ok(value)
    }
}

// Exact significand/exponent separation, including subnormal source entries.
fn parts(value: f64) -> (f64, i32) {
    if value == 0.0 {
        return (value, 0);
    }
    let bits = value.to_bits();
    let encoded = ((bits >> 52) & 0x7ff) as i32;
    if encoded == 0 {
        let (fraction, exponent) = parts(value * 18_014_398_509_481_984.0); // 2^54
        return (fraction, exponent - 54);
    }
    let fraction = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (fraction, encoded - 1023)
}

fn binary_scale(value: f64, shift: i32) -> f64 {
    if value == 0.0 || !value.is_finite() {
        return value;
    }
    let (fraction, exponent) = parts(value);
    let exponent = exponent + shift;
    if exponent > 1023 {
        return f64::INFINITY.copysign(value);
    }
    if exponent < -1075 {
        return 0.0_f64.copysign(value);
    }
    if exponent >= -1022 {
        fraction * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else {
        // First product is normal and exact; only the final subnormal rounds.
        (fraction * f64::MIN_POSITIVE) * f64::from_bits(((exponent + 1022 + 1023) as u64) << 52)
    }
}
