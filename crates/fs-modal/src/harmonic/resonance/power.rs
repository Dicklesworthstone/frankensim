//! Stationary extrema of actual signed cycle-mean physical power.
//!
//! These queries root P_omega, never the derivative of P squared or of a
//! fabricated amplitude. They reuse the original numerical sign-search owner.
//! Source, load and constitutive domains remain the evaluator's responsibility.

use super::super::{
    HarmonicFrequencyCurvatureResponse, HarmonicFrequencyParameterResponse,
    HarmonicFrequencyResponse, HarmonicResponseError, invalid, parameter::PowerSum,
};
use super::query::{Query, QuerySample, valid_interval};
use super::{ResonanceControls, ResonanceError};
use fs_exec::Cx;
use fs_math::c64::C64;

mod bandwidth;
pub use bandwidth::{
    refine_power_bandwidth, refine_power_bandwidth_with_cancel,
    refine_power_bandwidth_with_parameter_derivative,
    refine_power_bandwidth_with_parameter_derivative_with_cancel,
};

/// Actual signed cycle-mean power and its complete fixed-parameter slope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPowerSample {
    /// Signed physical power [W], such as input, viscosity or net boundary work.
    pub average_power_w: f64,
    /// Partial derivative of the same power per rad/s.
    pub angular_derivative: f64,
}

impl QuerySample for HarmonicPowerSample {
    const NONFINITE_MESSAGE: &'static str =
        "finite signed physical power and its frequency derivative required";

    fn finite(self) -> bool {
        self.average_power_w.is_finite() && self.angular_derivative.is_finite()
    }
}

/// One original signed power sample with its complete frequency curvature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPowerCurvatureSample {
    pub sample: HarmonicPowerSample,
    /// Partial second derivative of power per (rad/s)^2.
    pub angular_second_derivative: f64,
}

/// Complete scalar physical-power jet at fixed frequency and parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicPowerParameterJet {
    pub sample: HarmonicPowerSample,
    /// P_omegaomega, per (rad/s)^2.
    pub angular_second_derivative: f64,
    /// P_theta at fixed physical frequency, per caller-declared parameter unit.
    pub parameter_derivative: f64,
    /// P_omegatheta, per rad/s per parameter unit.
    pub angular_parameter_derivative: f64,
}

/// Orientation of an extremum of the signed physical power itself.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PowerExtremumKind {
    Maximum,
    Minimum,
}

/// One locally resolved stationary physical-power extremum.
///
/// The bracket is numerical sign evidence, not an interval enclosure of model
/// or evaluator error. Negative and zero powers retain their original meaning.
#[derive(Clone, Copy, Debug)]
pub struct PowerExtremum {
    pub kind: PowerExtremumKind,
    pub angular_frequency_rad_s: f64,
    /// The original sample reproduced by the accepted final curvature/jet call.
    pub sample: HarmonicPowerSample,
    /// Nonzero original curvature with the required maximum/minimum sign.
    pub angular_second_derivative: f64,
    pub bracket_rad_s: (f64, f64),
    /// All original evaluator calls, including the one final complete sample.
    pub evaluations: usize,
}

/// Local continuation of the same stationary physical-power extremum.
#[derive(Clone, Copy, Debug)]
pub struct PowerExtremumParameterResponse {
    pub extremum: PowerExtremum,
    /// -P_omegatheta/P_omegaomega, in rad/s per parameter unit.
    pub angular_frequency_derivative_rad_s: f64,
    /// TOTAL P_theta + P_omega * d(omega_extremum)/dtheta, in W per parameter.
    /// The finite search sample's P_omega is retained; it is not forced to zero.
    pub average_power_derivative: f64,
}

#[derive(Clone, Copy)]
struct Search {
    kind: PowerExtremumKind,
    frequency: f64,
    sample: HarmonicPowerSample,
    bracket: (f64, f64),
    evaluations: usize,
}

/// Refine ONE isolated extremum of actual signed physical power.
///
/// Only P and P_omega are evaluated during the shared bounded sign search.
/// One call is reserved for complete P/P_omega/P_omegaomega at the returned
/// frequency. It must reproduce the search sample within 32 machine epsilons,
/// with no unity floor or allowance for zero versus nonzero power/slope.
///
/// Both endpoints must have strictly opposite slopes. A maximum requires
/// negative final curvature; a minimum requires positive final curvature.
/// Thus a divergent passive P~1/(omega-pole)^2 slope reversal is not reported
/// as a finite maximum. A sampled pole propagates its original evaluator error.
/// Flat/degenerate extrema and inconsistent final samples refuse explicitly.
///
/// The caller supplies one continuous physical branch and isolates the extremum;
/// an unresolved pole or additional extremum inside a coarse bracket is not
/// ruled out globally. No absolute power, artificial loss, or heat/storage
/// interpretation is introduced. Cancellation or any refusal publishes no result.
pub fn refine_power_extremum_with_cancel<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    mut evaluate_curvature: impl FnMut(f64) -> Result<HarmonicPowerCurvatureSample, E>,
) -> Result<PowerExtremum, ResonanceError<E>> {
    poll(&mut cancelled)?;
    let search = search(bracket_rad_s, controls, &mut cancelled, evaluate)?;
    poll(&mut cancelled)?;
    let sample = evaluate_curvature(search.frequency).map_err(ResonanceError::Evaluation)?;
    poll(&mut cancelled)?;
    let result = resolve(search, sample)?;
    poll(&mut cancelled)?;
    Ok(result)
}

/// Refine and differentiate the same signed physical-power extremum.
///
/// The first callback supplies only P/P_omega; exactly one reserved final call
/// supplies P/P_omega/P_omegaomega/P_theta/P_omegatheta from the same physical
/// source, model, work channel and parameter. All constitutive, attachment and
/// source derivatives belong in that jet. A peak of net absorbed boundary work
/// need not coincide with a displacement, transmitted-effort or pressure peak.
///
/// At a simple stationary extremum d(omega)/dtheta=-P_omegatheta/P_omegaomega.
/// The returned power tangent includes the shift term at the actual numerical
/// sample. No P denominator is used, so signed and exactly zero extrema retain
/// ordinary derivatives. Singular/incorrectly oriented curvature and
/// unrepresentable derivatives refuse. This local derivative does not track
/// peak switching, topology changes or a finite parameter step.
pub fn refine_power_extremum_with_parameter_derivative_with_cancel<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    mut evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<PowerExtremumParameterResponse, ResonanceError<E>> {
    poll(&mut cancelled)?;
    let search = search(bracket_rad_s, controls, &mut cancelled, evaluate)?;
    poll(&mut cancelled)?;
    let jet = evaluate_jet(search.frequency).map_err(ResonanceError::Evaluation)?;
    poll(&mut cancelled)?;
    if !jet.parameter_derivative.is_finite() || !jet.angular_parameter_derivative.is_finite() {
        return Err(ResonanceError::Invalid(
            "finite complete physical power parameter jet required",
        ));
    }
    let extremum = resolve(
        search,
        HarmonicPowerCurvatureSample {
            sample: jet.sample,
            angular_second_derivative: jet.angular_second_derivative,
        },
    )?;
    let frequency_derivative = -jet.angular_parameter_derivative / jet.angular_second_derivative;
    if !frequency_derivative.is_finite()
        || (jet.angular_parameter_derivative != 0.0 && frequency_derivative == 0.0)
    {
        return Err(ResonanceError::Invalid(
            "physical power extremum frequency derivative is not representable",
        ));
    }
    // Reuse the physical signed-work owner: the product may overflow before
    // cancellation with P_theta even when the total derivative is finite.
    let mut total = PowerSum::default();
    let power_derivative = (|| {
        total.add_scalar(jet.parameter_derivative)?;
        total.half_inner(
            C64::from_re(jet.sample.angular_derivative),
            C64::from_re(frequency_derivative),
            2.0,
        )?;
        total.value()
    })()
    .map_err(|_| {
        ResonanceError::Invalid("total physical power extremum derivative is not representable")
    })?;
    poll(&mut cancelled)?;
    Ok(PowerExtremumParameterResponse {
        extremum,
        angular_frequency_derivative_rad_s: frequency_derivative,
        average_power_derivative: power_derivative,
    })
}

fn search<E>(
    bracket: (f64, f64),
    mut controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
) -> Result<Search, ResonanceError<E>> {
    let maximum_evaluations = controls.maximum_evaluations;
    if maximum_evaluations < 4 {
        return Err(ResonanceError::Invalid(
            "evaluation capacity must include the search and final power curvature/parameter jet",
        ));
    }
    controls.maximum_evaluations -= 1;
    let run = || {
        let mut query = Query::new(controls, cancelled, evaluate)?;
        valid_interval(bracket)?;
        let slope = |sample: HarmonicPowerSample| sample.angular_derivative;
        let left = query.sample(bracket.0, slope)?;
        let right = query.sample(bracket.1, slope)?;
        let kind = match (left.level > 0.0, right.level > 0.0) {
            (true, false) if right.level < 0.0 => PowerExtremumKind::Maximum,
            (false, true) if left.level < 0.0 => PowerExtremumKind::Minimum,
            _ => {
                return Err(ResonanceError::InvalidBracket(
                    "strict opposite physical power slopes required",
                ));
            }
        };
        // Actual power is signed and may have a finite extremum at zero.
        // No log-amplitude guard is applicable. The required original final
        // curvature distinguishes a simple extremum from a pole/flat reversal.
        let (point, bracket) = query.refine(left, right, slope, false)?;
        Ok(Search {
            kind,
            frequency: point.frequency,
            sample: point.sample,
            bracket,
            evaluations: query.evaluations,
        })
    };
    run().map_err(|error| match error {
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations,
        },
        other => other,
    })
}

fn resolve<E>(
    search: Search,
    final_sample: HarmonicPowerCurvatureSample,
) -> Result<PowerExtremum, ResonanceError<E>> {
    let sample = final_sample.sample;
    let curvature = final_sample.angular_second_derivative;
    if !sample.finite() || !curvature.is_finite() {
        return Err(ResonanceError::Invalid(
            "finite complete physical power frequency curvature required",
        ));
    }
    if !same(sample.average_power_w, search.sample.average_power_w)
        || !same(sample.angular_derivative, search.sample.angular_derivative)
    {
        return Err(ResonanceError::Invalid(
            "final power jet must reproduce the same original signed power and frequency derivative",
        ));
    }
    let oriented = match search.kind {
        PowerExtremumKind::Maximum => curvature < 0.0,
        PowerExtremumKind::Minimum => curvature > 0.0,
    };
    if !oriented {
        return Err(ResonanceError::Invalid(
            "nonzero physical power curvature with the original extremum orientation required",
        ));
    }
    Ok(PowerExtremum {
        kind: search.kind,
        angular_frequency_rad_s: search.frequency,
        sample,
        angular_second_derivative: curvature,
        bracket_rad_s: search.bracket,
        evaluations: search.evaluations + 1,
    })
}

fn same(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs());
    scale == 0.0 || (a / scale - b / scale).abs() <= 32.0 * f64::EPSILON
}

fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ResonanceError<E>> {
    if cancelled() {
        Err(ResonanceError::Cancelled)
    } else {
        Ok(())
    }
}

/// Executor-context entry sharing the same signed-power query owner.
pub fn refine_power_extremum<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_curvature: impl FnMut(f64) -> Result<HarmonicPowerCurvatureSample, E>,
) -> Result<PowerExtremum, ResonanceError<E>> {
    refine_power_extremum_with_cancel(
        bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_curvature,
    )
}

/// Executor-context entry sharing the same complete signed-power sensitivity query.
pub fn refine_power_extremum_with_parameter_derivative<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<PowerExtremumParameterResponse, ResonanceError<E>> {
    refine_power_extremum_with_parameter_derivative_with_cancel(
        bracket_rad_s,
        controls,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}

impl HarmonicFrequencyResponse {
    /// Observe signed cycle-mean source work and its complete frequency slope.
    pub fn input_power_sample(&self) -> HarmonicPowerSample {
        HarmonicPowerSample {
            average_power_w: self.response.average_input_power_w,
            angular_derivative: self.average_input_power_derivative,
        }
    }

    /// Observe structural viscous work, with gyroscopic direct work excluded.
    pub fn viscous_power_sample(&self) -> HarmonicPowerSample {
        HarmonicPowerSample {
            average_power_w: self.response.average_viscous_power_w,
            angular_derivative: self.average_viscous_power_derivative,
        }
    }

    /// Observe signed net load-port work, without inferring its heat/storage split.
    pub fn port_power_sample(&self) -> Result<HarmonicPowerSample, HarmonicResponseError> {
        let (Some(original), Some(tangent)) = (&self.response.port_load, &self.port_load) else {
            return Err(invalid(
                "complete original load-port power and frequency slope required",
            ));
        };
        Ok(HarmonicPowerSample {
            average_power_w: original.average_absorbed_power_w,
            angular_derivative: tangent.average_absorbed_power_derivative,
        })
    }
}

impl HarmonicFrequencyCurvatureResponse {
    pub fn input_power_curvature_sample(&self) -> HarmonicPowerCurvatureSample {
        HarmonicPowerCurvatureSample {
            sample: self.frequency_response.input_power_sample(),
            angular_second_derivative: self.average_input_power_second_derivative,
        }
    }

    pub fn viscous_power_curvature_sample(&self) -> HarmonicPowerCurvatureSample {
        HarmonicPowerCurvatureSample {
            sample: self.frequency_response.viscous_power_sample(),
            angular_second_derivative: self.average_viscous_power_second_derivative,
        }
    }

    pub fn port_power_curvature_sample(
        &self,
    ) -> Result<HarmonicPowerCurvatureSample, HarmonicResponseError> {
        let curvature = self
            .port_load
            .as_ref()
            .ok_or_else(|| invalid("complete load-port power frequency curvature required"))?;
        Ok(HarmonicPowerCurvatureSample {
            sample: self.frequency_response.port_power_sample()?,
            angular_second_derivative: curvature.average_absorbed_power_second_derivative,
        })
    }
}

impl HarmonicFrequencyParameterResponse {
    pub fn input_power_parameter_jet(&self) -> HarmonicPowerParameterJet {
        HarmonicPowerParameterJet {
            sample: self
                .frequency_response
                .frequency_response
                .input_power_sample(),
            angular_second_derivative: self
                .frequency_response
                .average_input_power_second_derivative,
            parameter_derivative: self.parameter_response.average_input_power_derivative,
            angular_parameter_derivative: self.average_input_power_mixed_derivative,
        }
    }

    pub fn viscous_power_parameter_jet(&self) -> HarmonicPowerParameterJet {
        HarmonicPowerParameterJet {
            sample: self
                .frequency_response
                .frequency_response
                .viscous_power_sample(),
            angular_second_derivative: self
                .frequency_response
                .average_viscous_power_second_derivative,
            parameter_derivative: self.parameter_response.average_viscous_power_derivative,
            angular_parameter_derivative: self.average_viscous_power_mixed_derivative,
        }
    }

    pub fn port_power_parameter_jet(
        &self,
    ) -> Result<HarmonicPowerParameterJet, HarmonicResponseError> {
        let curvature = self.frequency_response.port_power_curvature_sample()?;
        let (Some(parameter), Some(mixed)) = (&self.parameter_response.port_load, &self.port_load)
        else {
            return Err(invalid(
                "complete load-port power parameter and mixed derivatives required",
            ));
        };
        Ok(HarmonicPowerParameterJet {
            sample: curvature.sample,
            angular_second_derivative: curvature.angular_second_derivative,
            parameter_derivative: parameter.average_absorbed_power_derivative,
            angular_parameter_derivative: mixed.average_absorbed_power_mixed_derivative,
        })
    }
}
