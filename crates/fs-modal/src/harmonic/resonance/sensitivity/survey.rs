//! Parameter derivatives of every transfer extremum resolved by an authored grid.

use super::super::query::Query;
use super::super::{
    HarmonicTransferSample, ResonanceControls, ResonanceError, TransferExtremum,
    TransferExtremumKind,
};
use super::{
    HarmonicTransferParameterJet, TransferExtremumParameterResponse, original_budget, resolve,
};
use fs_exec::Cx;

/// Original observed peaks/notches and their local physical design derivatives.
///
/// Entries are ordered by frequency. The supplied grid must resolve every
/// feature the caller needs; equal-sign endpoints may hide a peak/notch pair.
/// This is not mode matching or a free-spectrum/global completeness claim.
#[derive(Clone, Debug)]
pub struct TransferParameterSurvey {
    pub extrema: Vec<TransferExtremumParameterResponse>,
    pub interval_rad_s: (f64, f64),
    pub grid_points: usize,
    /// All unique grid/refinement first-callback calls plus one final jet per
    /// returned extremum. Shared input samples are counted only once.
    pub evaluations: usize,
}

/// Resolve and differentiate all actual transfer extrema separated by a grid.
///
/// First callbacks supply the original coherent H/H_omega. The shared survey
/// owner caches every grid sample, discovers strictly oriented candidates and
/// enforces the extremum cap before refinement. One complete analytic parameter
/// jet is then reserved for EVERY candidate within the same global evaluation
/// limit; full jets are never evaluated at every refinement trial.
///
/// The same complete physical model, source, receiver and parameter must be
/// represented by both callbacks. Each final jet reproduces the original
/// H/H_omega sample and passes the existing signed-curvature/representability
/// checks. A simple exact cancellation zero retains its notch-frequency and
/// total complex derivatives, with absent log-amplitude and phase derivatives.
///
/// Original grid samples and exact stationary samples are reused. Each entry's
/// evaluation count includes its new refinement calls and one final jet; shared
/// grid calls belong to the survey total. Memory is linear in grid points and
/// resolved extrema, with no dense jet cache, perturbed solves or new root owner.
///
/// Band endpoints and silent/flat observations yield no invented features.
/// Unresolved pairs between grid points remain the caller's resolution limit.
/// Sampled poles retain the evaluator error and maxima retain the original
/// finite-response guard. A failed jet, budget or cancellation returns no
/// partial survey. Derivatives continue each isolated simple feature locally;
/// finite parameter changes, peak switching and topology changes are not inferred.
pub fn survey_transfer_extrema_with_parameter_derivative_with_cancel<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    mut evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<TransferParameterSurvey, ResonanceError<E>> {
    let run = || {
        let mut query = Query::new(controls, cancelled, evaluate)?;
        let brackets = query.survey_brackets(
            frequencies_rad_s,
            maximum_extrema,
            HarmonicTransferSample::slope,
        )?;
        let count = brackets.len();
        query.reserve_evaluations(count)?;
        let mut extrema = Vec::with_capacity(count);
        for (left, right) in brackets {
            query.poll()?;
            let left = query.grid[left];
            let right = query.grid[right];
            let kind = if left.level > 0.0 {
                TransferExtremumKind::Maximum
            } else {
                TransferExtremumKind::Minimum
            };
            let before = query.evaluations;
            let (point, bracket_rad_s) = query.refine(
                left,
                right,
                HarmonicTransferSample::slope,
                kind == TransferExtremumKind::Maximum,
            )?;
            let extremum = TransferExtremum {
                kind,
                angular_frequency_rad_s: point.frequency,
                sample: point.sample,
                bracket_rad_s,
                evaluations: query.evaluations - before,
            };
            query.poll()?;
            let jet = evaluate_jet(point.frequency).map_err(ResonanceError::Evaluation)?;
            query.poll()?;
            extrema.push(resolve(extremum, jet)?);
        }
        query.poll()?;
        Ok(TransferParameterSurvey {
            extrema,
            interval_rad_s: (
                frequencies_rad_s[0],
                frequencies_rad_s[frequencies_rad_s.len() - 1],
            ),
            grid_points: frequencies_rad_s.len(),
            evaluations: query.evaluations + count,
        })
    };
    run().map_err(|error| original_budget(error, controls.maximum_evaluations))
}

/// Executor-context entry sharing the same complete transfer-parameter survey.
pub fn survey_transfer_extrema_with_parameter_derivative<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
    evaluate_jet: impl FnMut(f64) -> Result<HarmonicTransferParameterJet, E>,
) -> Result<TransferParameterSurvey, ResonanceError<E>> {
    survey_transfer_extrema_with_parameter_derivative_with_cancel(
        frequencies_rad_s,
        controls,
        maximum_extrema,
        || cx.checkpoint().is_err(),
        evaluate,
        evaluate_jet,
    )
}
