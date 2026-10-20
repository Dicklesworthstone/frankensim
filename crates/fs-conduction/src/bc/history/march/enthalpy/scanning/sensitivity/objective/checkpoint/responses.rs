//! Objective and all state constraints share each replayed short tape.
//! Point/peak adjoint seeds coexist under the original gradient-entry cap.
//! Peak sample values, not full nodal states, are retained under max_peak_samples.
//! No extra forward solve is performed for an additional point or peak response.

use super::{ScanCheckpointBudget, ScanCheckpointWork, ScanObjectiveBudget,
    ScanObjectiveEvaluation, ScanObservation, ScanObjectiveSampling,
    ScanSensitivityError, ScanningEnthalpyRun, EnthalpyStepConfig, LinearConfig,
    ReverseState, ReplayedScan, admit_checkpoints, admit_layout, observed,
    same_scalar, finite, invalid, poll, zeros};
use super::super::add;
use super::super::responses::{ScanResponse, ScanResponseBudget, ScanResponseEvaluation};
use super::super::responses::peaks::{ScanPeakResponse, ScanPeakEvaluation,
    PeakSummary, admit_rows, summarize_samples, weight};
use fs_exec::Cx;
use super::super::spatial::ScanObservationSpec;

/// Original objective, physical constraint values/gradients and actual replay work.
#[derive(Debug)]
pub struct CheckpointedScanResponses {
    pub objective: ScanObjectiveEvaluation,
    pub responses: Vec<ScanResponseEvaluation>,
    pub peaks: Vec<ScanPeakEvaluation>,
    pub work: ScanCheckpointWork,
}

/// Checkpointed equivalent of evaluate_scan_objective_with_peaks.
///
/// Each segment is replayed ONCE, then all applicable adjoints traverse its
/// tangents before they are discarded. Gradient accumulation for each output
/// still proceeds in original reverse chronological order. The same physical
/// samples, peak reduction and endpoint transpose owners are reused. Original
/// and replayed selected response samples must match bitwise; no partial
/// objective, gradient or constraint report escapes a mismatch or refusal.
///
/// In addition to the short tape/checkpoints, this retains selected peak values
/// (max_peak_samples) and one running gradient per response (max_gradient_entries).
/// These are shape bounds, not an allocation peak certificate. Constraints still
/// cover only specified vertices and endpoints, not an unsampled maximum.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_objective_checkpointed_with_peaks(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    responses: &[ScanResponse], peaks: &[ScanPeakResponse], response_budget: ScanResponseBudget,
    max_peak_samples: usize, checkpoints: ScanCheckpointBudget,
) -> Result<CheckpointedScanResponses, ScanSensitivityError> {
    evaluate_observations_checkpointed(cx, run, endpoints, observations, sampling,
        step, transpose, budget, responses, peaks, response_budget, max_peak_samples, checkpoints)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate_observations_checkpointed<O: ScanObservationSpec>(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    responses: &[ScanResponse], peaks: &[ScanPeakResponse], response_budget: ScanResponseBudget,
    max_peak_samples: usize, checkpoints: ScanCheckpointBudget,
) -> Result<CheckpointedScanResponses, ScanSensitivityError> {
    let segment = admit_checkpoints(cx, &run, endpoints.len(), checkpoints)?;
    let (at_endpoint, per_step) = admit_layout(cx, &run, endpoints, observations, budget, segment)?;
    let sample_counts = admit_rows(cx, &run, endpoints.len(), responses, peaks, response_budget, max_peak_samples)?;
    let vertices = run.conduction().temperature().len();
    let mut predictions = zeros(observations.len())?;
    let mut residuals = zeros(observations.len())?;
    let mut point_values = zeros(responses.len())?;
    let mut peak_values: Vec<Vec<f64>> = Vec::with_capacity(peaks.len());
    for &count in &sample_counts { peak_values.push(zeros(count)?); }
    let mut value = 0.0;
    let mut replay = ReplayedScan::forward(cx, run, endpoints, sampling, step, segment,
        |index, report| {
            for &i in &at_endpoint[index] {
                poll(cx, i)?;
                let o = &observations[i];
                let prediction = o.prediction(cx, report)?;
                let residual = finite((prediction-o.target())/o.scale())?;
                let (term, _) = o.loss().value_derivative(residual)?;
                value = finite(value + term)?;
                predictions[i] = prediction; residuals[i] = residual;
            }
            for (i, response) in responses.iter().enumerate() {
                poll(cx, i)?;
                if index == response.endpoint {
                    point_values[i] = finite(observed(report, response.observable, response.vertex))?;
                }
            }
            for (peak, values) in peaks.iter().zip(&mut peak_values) {
                poll(cx, index)?;
                if index < peak.first_endpoint || index > peak.last_endpoint { continue; }
                let offset = (index-peak.first_endpoint)*peak.vertices.len();
                for (j, &vertex) in peak.vertices.iter().enumerate() {
                    if j % 256 == 0 { poll(cx, j)?; }
                    values[offset+j] = finite(observed(report, peak.observable, vertex))?;
                }
            }
            Ok(())
        })?;
    let mut summaries = Vec::with_capacity(peaks.len());
    for (peak, values) in peaks.iter().zip(&peak_values) {
        summaries.push(summarize_samples(cx, values.len(), peak.smoothing_scale, |i| values[i])?);
    }
    let mut objective = ReverseState::new(vertices)?;
    let mut point_reverse = Vec::with_capacity(responses.len());
    let mut peak_reverse = Vec::with_capacity(peaks.len());
    for _ in responses { poll(cx, point_reverse.len())?; point_reverse.push(ReverseState::new(vertices)?); }
    for _ in peaks { poll(cx, peak_reverse.len())?; peak_reverse.push(ReverseState::new(vertices)?); }
    for segment_index in (0..replay.segments()).rev() {
        let (begin, tape) = replay.segment(cx, segment_index, endpoints, sampling, step, per_step)?;
        for local in (0..tape.len()).rev() {
            let index = begin+local;
            let binding = &tape[local];
            for &i in &at_endpoint[index] {
                poll(cx, i)?;
                let o = &observations[i];
                same_scalar(o.prediction(cx, binding.report())?, predictions[i])?;
            }
            objective.inject_observations(cx, binding, observations, &at_endpoint[index], &residuals)?;
            objective.pull(cx, &replay.run, binding, sampling, transpose, budget.max_schedule_controls)?;
            for (i, (response, reverse)) in responses.iter().zip(&mut point_reverse).enumerate() {
                poll(cx, i)?;
                if index > response.endpoint { continue; }
                if index == response.endpoint {
                    same_scalar(observed(binding.report(), response.observable, response.vertex), point_values[i])?;
                    reverse.seed[response.vertex] = 1.0;
                    reverse.seed = response.observable.seed_in_enthalpy(cx, binding,
                        std::mem::take(&mut reverse.seed))?;
                }
                reverse.pull(cx, &replay.run, binding, sampling, transpose, budget.max_schedule_controls)?;
            }
            for (i, ((peak, reverse), summary)) in peaks.iter().zip(&mut peak_reverse).zip(&summaries).enumerate() {
                poll(cx, i)?;
                if index > peak.last_endpoint { continue; }
                if index >= peak.first_endpoint {
                    let offset = (index-peak.first_endpoint)*peak.vertices.len();
                    let mut seed = zeros(vertices)?;
                    for (j, &vertex) in peak.vertices.iter().enumerate() {
                        if j % 256 == 0 { poll(cx, j)?; }
                        let sample = observed(binding.report(), peak.observable, vertex);
                        same_scalar(sample, peak_values[i][offset+j])?;
                        seed[vertex] = weight(sample, summary.maximum, peak.smoothing_scale)/summary.sum;
                    }
                    seed = peak.observable.seed_in_enthalpy(cx, binding, seed)?;
                    // The carried downstream seed is already in h coordinates.
                    add(cx, &mut reverse.seed, &seed)?;
                }
                reverse.pull(cx, &replay.run, binding, sampling, transpose, budget.max_schedule_controls)?;
            }
        }
    }
    let mut point_results = Vec::with_capacity(responses.len());
    for ((&response, value), reverse) in responses.iter().zip(point_values).zip(point_reverse) {
        poll(cx, point_results.len())?;
        point_results.push(ScanResponseEvaluation { response, value,
            sources: reverse.sources.ok_or_else(|| invalid("missing checkpointed point response sweep"))?,
            initial_specific_enthalpy: reverse.seed, transpose_iterations: reverse.iterations,
            max_transpose_relative_residual: reverse.residual });
    }
    let mut peak_results = Vec::with_capacity(peaks.len());
    for (((peak, reverse), summary), &sample_count) in peaks.iter().zip(peak_reverse).zip(summaries).zip(&sample_counts) {
        poll(cx, peak_results.len())?;
        let PeakSummary { maximum, gap, value, bound, .. } = summary;
        peak_results.push(ScanPeakEvaluation { peak: peak.clone(), value,
            sampled_maximum: maximum, smoothing_gap: gap, maximum_smoothing_gap: bound, sample_count,
            sources: reverse.sources.ok_or_else(|| invalid("missing checkpointed peak response sweep"))?,
            initial_specific_enthalpy: reverse.seed, transpose_iterations: reverse.iterations,
            max_transpose_relative_residual: reverse.residual });
    }
    poll(cx, endpoints.len())?;
    Ok(CheckpointedScanResponses {
        objective: ScanObjectiveEvaluation { value, predictions, standardized_residuals: residuals,
            sources: objective.sources.ok_or_else(|| invalid("missing checkpointed objective sweep"))?,
            initial_specific_enthalpy: objective.seed, final_step: replay.final_step,
            endpoint_count: endpoints.len(), transpose_iterations: objective.iterations,
            max_transpose_relative_residual: objective.residual },
        responses: point_results, peaks: peak_results, work: replay.work,
    })
}
