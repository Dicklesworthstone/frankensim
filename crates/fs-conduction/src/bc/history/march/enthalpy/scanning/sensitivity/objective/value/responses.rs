//! Point and sampled-peak values without assuming a differentiable trajectory.
//!
//! The fitting loss and every requested response share one primal pass. Only
//! the selected scalar peak samples are retained; no endpoint tangent, full
//! field history, transpose solve or replay is needed. Peak reduction is the
//! SAME ordered log-sum-exp used by the derivative path, not a new estimate.

use super::{ScanValueBudget, ScanValueEvaluation, ScanObjectiveSampling,
    ScanSensitivityError, ScanningEnthalpyRun, EnthalpyStepConfig,
    evaluate_with_observer, finite, invalid, poll, zeros};
use super::super::{observed, spatial::ScanObservationSpec};
use super::super::responses::ScanResponse;
use super::super::responses::peaks::{ScanPeakResponse, summarize_samples};
use fs_exec::Cx;

/// Value-only shape allowances, independent of any derivative-storage budget.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanValueResponseBudget {
    /// Combined number of point and peak responses.
    pub max_responses: usize,
    /// Total selected vertex/time pairs retained for all peak responses.
    pub max_peak_samples: usize,
}

/// A physical response, with no fabricated gradient at a constitutive corner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanResponseValue {
    pub response: ScanResponse,
    /// Kelvin, J/kg or dimensionless liquid mass fraction, as declared.
    pub value: f64,
}

/// Same sampled maximum and smoothing convention as the adjoint peak report.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanPeakValue {
    pub peak: ScanPeakResponse,
    pub value: f64,
    pub sampled_maximum: f64,
    pub smoothing_gap: f64,
    pub maximum_smoothing_gap: f64,
    pub sample_count: usize,
}

/// Evaluate fitting and all point/peak requirements on ONE physical trajectory.
///
/// Invalid selections or shape allowances refuse before the first endpoint.
/// Valid phase-corner values remain available even when the corresponding
/// derivative does not exist. Nonfinite data, physical failures and cancellation
/// return no partial result. Robust fitting losses never modify response values.
///
/// A peak covers only its explicit vertices and endpoint window. Its smoothing
/// margin is not an interval certificate or a bound on unsampled physical peaks.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_values_with_responses<O: ScanObservationSpec>(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[O], sampling: ScanObjectiveSampling, step: EnthalpyStepConfig,
    budget: ScanValueBudget, responses: &[ScanResponse], peaks: &[ScanPeakResponse],
    response_budget: ScanValueResponseBudget,
) -> Result<(ScanValueEvaluation, Vec<ScanResponseValue>, Vec<ScanPeakValue>), ScanSensitivityError> {
    poll(cx, 0)?;
    if responses.len().checked_add(peaks.len())
        .is_none_or(|n| n > response_budget.max_responses)
    { return Err(invalid("value-only point/peak response allowance exceeded")); }
    // Response admission precedes peak allocation. The shared forward owner
    // still admits the original objective and run before solving any endpoint.
    let vertices = run.conduction().temperature().len();
    for (i, response) in responses.iter().enumerate() {
        poll(cx, i)?;
        if response.endpoint >= endpoints.len() || response.vertex >= vertices {
            return Err(invalid("value-only response must name an existing endpoint and vertex"));
        }
    }
    let mut sample_counts = Vec::new();
    sample_counts.try_reserve_exact(peaks.len())
        .map_err(|_| invalid("value-only peak count allocation refused"))?;
    let mut samples = 0usize;
    for peak in peaks {
        let count = peak.validate(cx, endpoints.len(), vertices,
            response_budget.max_peak_samples - samples)?;
        samples += count;
        sample_counts.push(count);
    }
    let mut point_values = zeros(responses.len())?;
    let mut peak_values = Vec::new();
    peak_values.try_reserve_exact(peaks.len())
        .map_err(|_| invalid("value-only peak sample allocation refused"))?;
    for &count in &sample_counts { poll(cx, count)?; peak_values.push(zeros(count)?); }
    let objective = evaluate_with_observer(cx, run, endpoints, observations, sampling,
        step, budget, |index, report| {
            for (i, response) in responses.iter().enumerate() {
                poll(cx, i)?;
                if index == response.endpoint {
                    point_values[i] = finite(observed(report, response.observable, response.vertex))?;
                }
            }
            for (peak, values) in peaks.iter().zip(&mut peak_values) {
                poll(cx, index)?;
                if index < peak.first_endpoint || index > peak.last_endpoint { continue; }
                let offset = (index - peak.first_endpoint) * peak.vertices.len();
                for (j, &vertex) in peak.vertices.iter().enumerate() {
                    if j % 256 == 0 { poll(cx, j)?; }
                    values[offset + j] = finite(observed(report, peak.observable, vertex))?;
                }
            }
            Ok(())
        })?;
    let mut points = Vec::new();
    points.try_reserve_exact(responses.len())
        .map_err(|_| invalid("value-only response report allocation refused"))?;
    for (&response, value) in responses.iter().zip(point_values) {
        poll(cx, points.len())?;
        points.push(ScanResponseValue { response, value });
    }
    let mut output = Vec::new();
    output.try_reserve_exact(peaks.len())
        .map_err(|_| invalid("value-only peak report allocation refused"))?;
    for (peak, values) in peaks.iter().zip(peak_values) {
        let summary = summarize_samples(cx, values.len(), peak.smoothing_scale, |i| values[i])?;
        output.push(ScanPeakValue { peak: peak.clone(), value: summary.value,
            sampled_maximum: summary.maximum, smoothing_gap: summary.gap,
            maximum_smoothing_gap: summary.bound, sample_count: values.len() });
    }
    poll(cx, endpoints.len())?;
    Ok((objective, points, output))
}
