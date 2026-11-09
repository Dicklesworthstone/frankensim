//! Actual positive physical-power half-height crossings and design derivatives.

use super::super::query::{Query, QuerySample, valid_interval};
use super::super::sensitivity::Products;
use super::super::{
    HalfPowerBandwidth, HalfPowerBandwidthParameterResponse, ResonanceControls, ResonanceError,
};
use super::{
    HarmonicPowerParameterJet, HarmonicPowerSample, PowerExtremum, PowerExtremumKind,
    PowerExtremumParameterResponse, poll, same,
};
use fs_exec::Cx;

/// Resolve both actual half-height crossings of one positive physical-power peak.
///
/// The same immutable physical work channel, model and source must be used for
/// the peak and both flanks. Each flank must separately isolate a continuous
/// crossing of P_peak/2, rising below the peak and falling above it. The shared
/// numerical query owner refines the signed level P-P_peak/2 directly, with
/// common scaling before subtraction. No amplitude, squared-power proxy or
/// assumed symmetric/Lorentzian line shape is introduced.
///
/// A positive finite maximum at positive frequency is required. Its original
/// curvature must be finite and negative. Missing, stationary or wrongly
/// oriented crossings and touching final crossing brackets refuse. The returned
/// Q is peak frequency divided by the measured width; it is not an inferred
/// material damping ratio. Brackets retain the shared numerical, non-certified
/// interpretation and the caller's responsibility to isolate each feature.
///
/// One evaluation cap covers both flanks. A shared peak endpoint is sampled
/// once, exact crossing samples retain the query owner's constant-size cache,
/// and no final full jet is needed by this value-only entry. Original evaluator
/// failures and cancellation publish no partial bandwidth.
pub fn refine_power_bandwidth_with_cancel<E>(
    peak: &PowerExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
) -> Result<HalfPowerBandwidth, ResonanceError<E>> {
    refine_samples(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        cancelled,
        evaluate,
    )
    .map(|(bandwidth, _)| bandwidth)
}

/// Differentiate both physical-power half-height crossings, measured width and Q.
///
/// Exactly two final edge callbacks are reserved BEFORE refinement, within the
/// same total evaluation ceiling. Each must reproduce its original P/P_omega
/// sample within 32 machine epsilons without a unity floor. Only the sample and
/// P_theta entries are needed from each supplied parameter jet; curvature and
/// mixed derivatives do not enter a simple level-crossing derivative.
///
/// At an edge, omega_theta=(P_peak_total_theta/2-P_theta)/P_omega. The TOTAL
/// moving-peak power derivative therefore moves the threshold, and a common
/// source gain does not invent a linewidth change. Both edge slopes must remain
/// nonzero and oriented toward the same peak. Original parameter, source and
/// material units are retained through the shared normalized DD contractions,
/// including the width and Q ratios. No finite differences or extra perturbed
/// physical solves are introduced.
///
/// The prior peak query has its own evaluation count. This query is local to
/// the same smoothly continuing peak and isolated crossings; it does not resolve
/// feature switching, a flat crossing, or finite parameter steps.
pub fn refine_power_bandwidth_with_parameter_derivative_with_cancel<E>(
    peak: &PowerExtremumParameterResponse,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    mut evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<E>> {
    poll(&mut cancelled)?;
    if !peak.average_power_derivative.is_finite()
        || !peak.angular_frequency_derivative_rad_s.is_finite()
    {
        return Err(ResonanceError::Invalid(
            "finite original physical-power peak parameter response required",
        ));
    }
    if controls.maximum_evaluations < 5 {
        return Err(ResonanceError::Invalid(
            "evaluation capacity must include the power-flank search and two edge jets",
        ));
    }
    let mut search_controls = controls;
    search_controls.maximum_evaluations -= 2;
    let (bandwidth, samples) = refine_samples(
        &peak.extremum,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        search_controls,
        &mut cancelled,
        evaluate,
    )
    .map_err(|error| match error {
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations: controls.maximum_evaluations,
        },
        other => other,
    })?;

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
        if !jet.sample.finite() || !jet.parameter_derivative.is_finite() {
            return Err(ResonanceError::Invalid(
                "finite original physical-power crossing parameter data required",
            ));
        }
        if !same(jet.sample.average_power_w, samples[index].average_power_w)
            || !same(
                jet.sample.angular_derivative,
                samples[index].angular_derivative,
            )
        {
            return Err(ResonanceError::Invalid(
                "power edge jet must reproduce the same original crossing power and frequency derivative",
            ));
        }
        let slope = jet.sample.angular_derivative;
        if (index == 0 && slope <= 0.0) || (index == 1 && slope >= 0.0) {
            return Err(ResonanceError::InvalidBracket(
                "nonzero physical-power crossing slopes toward the same peak required",
            ));
        }
        derivatives[index] = Products::sum(&[
            (peak.average_power_derivative, 1.0, -1),
            (-jet.parameter_derivative, 1.0, 0),
        ])
        .ratio(Products::sum(&[(slope, 1.0, 0)]))?;
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

/// Executor-context entry for actual positive physical-power bandwidth and Q.
pub fn refine_power_bandwidth<E>(
    peak: &PowerExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
) -> Result<HalfPowerBandwidth, ResonanceError<E>> {
    refine_power_bandwidth_with_cancel(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
    )
}

/// Executor-context entry for both power crossings, bandwidth and Q sensitivities.
pub fn refine_power_bandwidth_with_parameter_derivative<E>(
    peak: &PowerExtremumParameterResponse,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<HalfPowerBandwidthParameterResponse, ResonanceError<E>> {
    refine_power_bandwidth_with_parameter_derivative_with_cancel(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}

/// Retain the two original crossing samples for the final edge-jet checks.
fn refine_samples<E>(
    peak: &PowerExtremum,
    lower_bracket: (f64, f64),
    upper_bracket: (f64, f64),
    controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
) -> Result<(HalfPowerBandwidth, [HarmonicPowerSample; 2]), ResonanceError<E>> {
    let mut query = Query::new(controls, cancelled, evaluate)?;
    valid_interval(lower_bracket)?;
    valid_interval(upper_bracket)?;
    if peak.kind != PowerExtremumKind::Maximum
        || !peak.angular_frequency_rad_s.is_finite()
        || peak.angular_frequency_rad_s <= 0.0
        || !peak.sample.finite()
        || peak.sample.average_power_w <= 0.0
        || !peak.angular_second_derivative.is_finite()
        || peak.angular_second_derivative >= 0.0
        || lower_bracket.1 > peak.angular_frequency_rad_s
        || upper_bracket.0 < peak.angular_frequency_rad_s
    {
        return Err(ResonanceError::Invalid(
            "positive finite physical-power maximum and separate flanks required",
        ));
    }
    let level = |sample: HarmonicPowerSample| {
        let scale = sample
            .average_power_w
            .abs()
            .max(peak.sample.average_power_w);
        sample.average_power_w / scale - 0.5 * (peak.sample.average_power_w / scale)
    };
    let l0 = query.sample(lower_bracket.0, level)?;
    let l1 = query.sample(lower_bracket.1, level)?;
    // The two isolated flanks often share the same peak endpoint. Retain its
    // original sample once without another model evaluation or a heap cache.
    let u0 = if upper_bracket.0 == lower_bracket.1 {
        l1
    } else {
        query.sample(upper_bracket.0, level)?
    };
    let u1 = query.sample(upper_bracket.1, level)?;
    if !(l0.level < 0.0 && l1.level > 0.0 && u0.level > 0.0 && u1.level < 0.0) {
        return Err(ResonanceError::InvalidBracket(
            "each physical-power flank must strictly bracket its half-height crossing",
        ));
    }
    let (lower, lower_bracket) = query.refine(l0, l1, level, false)?;
    let (upper, upper_bracket) = query.refine(u0, u1, level, false)?;
    if lower_bracket.1 >= upper_bracket.0 {
        return Err(ResonanceError::InvalidBracket(
            "power crossing brackets do not resolve a positive linewidth",
        ));
    }
    if lower.sample.angular_derivative <= 0.0 || upper.sample.angular_derivative >= 0.0 {
        return Err(ResonanceError::InvalidBracket(
            "nonzero physical-power crossing slopes toward the same peak required",
        ));
    }
    let width = upper.frequency - lower.frequency;
    let quality = peak.angular_frequency_rad_s / width;
    if width <= 0.0 || !quality.is_finite() || quality <= 0.0 {
        return Err(ResonanceError::Invalid(
            "physical-power linewidth or quality factor is not representable",
        ));
    }
    query.poll()?;
    Ok((
        HalfPowerBandwidth {
            lower_angular_frequency_rad_s: lower.frequency,
            upper_angular_frequency_rad_s: upper.frequency,
            angular_bandwidth_rad_s: width,
            response_quality_factor: quality,
            lower_bracket_rad_s: lower_bracket,
            upper_bracket_rad_s: upper_bracket,
            evaluations: query.evaluations,
        },
        [lower.sample, upper.sample],
    ))
}
