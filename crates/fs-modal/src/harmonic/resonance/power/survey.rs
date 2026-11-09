//! Grid-resolved extrema of original signed physical power.
//!
//! The complete grid and all refinements share the original query owner.
//! A final curvature or mixed physical jet is reserved for each discovered
//! candidate before any refinement; no partial collection is published.

use super::super::query::Query;
use super::{
    HarmonicPowerCurvatureSample, HarmonicPowerParameterJet, HarmonicPowerSample, PowerExtremum,
    PowerExtremumKind, PowerExtremumParameterResponse, ResonanceControls, ResonanceError, Search,
    resolve, resolve_parameter,
};
use fs_exec::Cx;

/// Simple signed physical-power extrema resolved by an authored frequency grid.
#[derive(Clone, Debug)]
pub struct PowerSurvey {
    /// Ascending frequency, retaining both signed maxima and minima.
    /// Each extremum's evaluation count excludes the shared cached grid.
    pub extrema: Vec<PowerExtremum>,
    /// First and last supplied angular frequency [rad/s].
    pub interval_rad_s: (f64, f64),
    /// Number of original grid frequencies, all evaluated on success.
    pub grid_points: usize,
    /// All first-sample and final-curvature calls across the complete survey.
    pub evaluations: usize,
}

/// Local physical-parameter continuation of every grid-resolved power extremum.
///
/// Each result concerns its own isolated simple extremum of the same original
/// power channel. This does not assign modes or track switches between extrema.
#[derive(Clone, Debug)]
pub struct PowerParameterSurvey {
    /// Ascending frequency, including each extremum's frequency and TOTAL
    /// observed power derivative at the actual returned numerical sample.
    pub extrema: Vec<PowerExtremumParameterResponse>,
    pub interval_rad_s: (f64, f64),
    pub grid_points: usize,
    /// All first-sample and final-complete-jet calls, including the shared grid.
    pub evaluations: usize,
}

/// Find the signed power maxima and minima resolved by a supplied frequency grid.
///
/// The grid is finite, strictly increasing and nonnegative, with at least two
/// entries. Every first callback supplies the same original P and complete
/// P_omega. Strict slope reversals identify candidates; stationary grid samples
/// remain between their strictly oriented neighbours and are not counted twice.
/// Exactly zero and negative powers retain their physical meaning.
///
/// The input grid is cached once. Each candidate uses the existing bounded
/// sign-search owner and one final original curvature call, checked against
/// the same sampled P/P_omega and required nonzero maximum/minimum curvature.
/// All final calls are reserved before refinement. The evaluation cap covers
/// grid, refinement and final calls together; the positive extremum cap refuses
/// truncation before any refinement or final call. Cancellation, a source
/// refusal, an unresolved flat interval or invalid curvature returns no list.
///
/// Storage is linear in grid and discovered extrema; each query keeps only
/// constant refinement scratch. A grid must separate the features of interest:
/// equal-sign samples can hide a pair, and a sign reversal alone does not prove
/// uniqueness or absence of a pole. Band endpoints are not two-sided extrema.
/// Constant power gives an empty survey. No absolute-value transform, modal
/// attribution, inferred linewidth or global completeness claim is introduced.
pub fn survey_power_extrema_with_cancel<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_curvature: impl FnMut(f64) -> Result<HarmonicPowerCurvatureSample, E>,
) -> Result<PowerSurvey, ResonanceError<E>> {
    let (extrema, evaluations) = survey(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        cancelled,
        evaluate,
        evaluate_curvature,
        resolve,
    )?;
    Ok(PowerSurvey {
        extrema,
        interval_rad_s: (
            frequencies_rad_s[0],
            frequencies_rad_s[frequencies_rad_s.len() - 1],
        ),
        grid_points: frequencies_rad_s.len(),
        evaluations,
    })
}

/// Find and differentiate every simple signed power extremum resolved by the grid.
///
/// First callbacks supply only P/P_omega. Exactly one reserved complete
/// P/P_omega/P_omegaomega/P_theta/P_omegatheta callback is made at each resolved
/// extremum, using the same physical source, work channel and parameter.
/// The existing local power owner validates these jets and computes
/// -P_omegatheta/P_omegaomega and the total P_theta+P_omega*omega_theta.
///
/// All candidates share the same cache, evaluation cap and atomic refusal as
/// survey_power_extrema_with_cancel. Full mixed jets are never requested at
/// every refinement trial. Local derivatives do not continue through a
/// degenerate extremum, feature switching or topology change.
pub fn survey_power_extrema_with_parameter_derivative_with_cancel<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<PowerParameterSurvey, ResonanceError<E>> {
    let (extrema, evaluations) = survey(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        cancelled,
        evaluate,
        evaluate_jet,
        resolve_parameter,
    )?;
    Ok(PowerParameterSurvey {
        extrema,
        interval_rad_s: (
            frequencies_rad_s[0],
            frequencies_rad_s[frequencies_rad_s.len() - 1],
        ),
        grid_points: frequencies_rad_s.len(),
        evaluations,
    })
}

/// Executor-context entry sharing the original signed-power survey owner.
pub fn survey_power_extrema<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_curvature: impl FnMut(f64) -> Result<HarmonicPowerCurvatureSample, E>,
) -> Result<PowerSurvey, ResonanceError<E>> {
    survey_power_extrema_with_cancel(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_curvature,
    )
}

/// Executor-context entry sharing the complete signed-power sensitivity survey.
pub fn survey_power_extrema_with_parameter_derivative<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicPowerParameterJet, E>,
) -> Result<PowerParameterSurvey, ResonanceError<E>> {
    survey_power_extrema_with_parameter_derivative_with_cancel(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}

#[allow(clippy::too_many_arguments)]
fn survey<E, J, R>(
    frequencies: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicPowerSample, E>,
    mut evaluate_final: impl FnMut(f64) -> Result<J, E>,
    mut resolve_final: impl FnMut(Search, J) -> Result<R, ResonanceError<E>>,
) -> Result<(Vec<R>, usize), ResonanceError<E>> {
    let run = || {
        let mut query = Query::new(controls, cancelled, evaluate)?;
        let slope = |sample: HarmonicPowerSample| sample.angular_derivative;
        let brackets = query.survey_brackets(frequencies, maximum_extrema, slope)?;

        let final_calls = brackets.len();
        query.reserve_evaluations(final_calls)?;
        let mut extrema = Vec::with_capacity(final_calls);
        for (left, right) in brackets {
            query.poll()?;
            let left = query.grid[left];
            let right = query.grid[right];
            let kind = if left.level > 0.0 {
                PowerExtremumKind::Maximum
            } else {
                PowerExtremumKind::Minimum
            };
            let before = query.evaluations;
            let (point, bracket) = query.refine(left, right, slope, false)?;
            let search = Search {
                kind,
                frequency: point.frequency,
                sample: point.sample,
                bracket,
                evaluations: query.evaluations - before,
            };
            query.poll()?;
            let final_sample =
                evaluate_final(point.frequency).map_err(ResonanceError::Evaluation)?;
            query.poll()?;
            let extremum = resolve_final(search, final_sample)?;
            query.poll()?;
            extrema.push(extremum);
        }
        query.poll()?;
        Ok((extrema, query.evaluations + final_calls))
    };
    run().map_err(|error| match error {
        ResonanceError::EvaluationBudget { .. } => ResonanceError::EvaluationBudget {
            maximum_evaluations: controls.maximum_evaluations,
        },
        other => other,
    })
}
