//! Smooth upper envelopes of explicitly sampled thermal hot spots.
//!
//! For sampled values y, p = max(y) + tau*ln(sum(exp((y-max(y))/tau))).
//! In exact arithmetic max(y) <= p <= max(y)+tau*ln(N). Unlike a log-MEAN-
//! exponential, this unnormalized sum cannot hide the largest selected sample.
//! Its softmax derivatives are injected at EVERY selected endpoint in one
//! reverse sweep. Neither a winner selection nor a smoothed forward field is
//! differentiated. The original nonlinear/radiative trajectory is unchanged.
//!
//! The envelope concerns supplied nodal samples, not unsampled time/space,
//! continuum error or machine safety. Floating-point evaluation is not an
//! outward-rounded certificate. Tau is fixed, positive, and in response units.

use super::{ScanResponse, ScanResponseBudget, ScanResponseEvaluation, pull_response};
use super::super::{ScanObservable, ScanObjectiveBudget, ScanObjectiveEvaluation, ScanObservation,
    ScanObjectiveSampling, ScanSensitivityError, ScanStepLinearization, ScanningEnthalpyRun,
    MovingHeatScheduleGradient, EnthalpyStepConfig, LinearConfig, add, evaluate_with_tape,
    finite, invalid, poll, zeros};
use fs_exec::Cx;

/// A maximum over a fixed vertex set and an inclusive range of retained times.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanPeakResponse {
    pub first_endpoint: usize,
    pub last_endpoint: usize,
    /// Nonempty, strictly increasing original vertex IDs. Explicitly supply all
    /// IDs to cover the full mesh; an empty selection never means "everywhere".
    pub vertices: Vec<usize>,
    pub observable: ScanObservable,
    /// Kelvin, J/kg or dimensionless mass fraction. Smaller values reduce the
    /// envelope bias but can increase optimization curvature. Never inferred.
    pub smoothing_scale: f64,
}

impl ScanPeakResponse {
    /// Admit sample count before traversing the vertex list or doing physics.
    /// The return is the number of selected vertex/time pairs, not a byte bound.
    pub fn validate(&self, cx: &Cx<'_>, endpoints: usize, vertices: usize,
        max_samples: usize) -> Result<usize, ScanSensitivityError>
    {
        poll(cx, 0)?;
        if self.first_endpoint > self.last_endpoint || self.last_endpoint >= endpoints
            || self.vertices.is_empty() || !self.smoothing_scale.is_finite() || self.smoothing_scale <= 0.0
        { return Err(invalid("scan peak needs a nonempty endpoint/vertex selection and positive finite physical smoothing scale")); }
        let count = (self.last_endpoint - self.first_endpoint + 1).checked_mul(self.vertices.len())
            .filter(|&n| n <= max_samples)
            .ok_or_else(|| invalid("scan peak samples exceed their explicit allowance"))?;
        let mut previous = None;
        for (i, &vertex) in self.vertices.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if vertex >= vertices || previous.is_some_and(|old| old >= vertex) {
                return Err(invalid("scan peak vertices must be distinct increasing IDs inside the original mesh"));
            }
            previous = Some(vertex);
        }
        poll(cx, count)?;
        Ok(count)
    }
}

/// Envelope value and gradients through one original physical trajectory.
#[derive(Debug)]
pub struct ScanPeakEvaluation {
    pub peak: ScanPeakResponse,
    /// Smooth upper envelope, in the selected response's original units.
    pub value: f64,
    /// Actual maximum among the selected samples, before adding smoothing bias.
    pub sampled_maximum: f64,
    /// Computed tau*ln(sum(exp(...))), not an estimate of physical solve error.
    pub smoothing_gap: f64,
    /// tau*ln(N): exact-arithmetic bias bound, evaluated in ordinary binary64.
    pub maximum_smoothing_gap: f64,
    pub sample_count: usize,
    pub sources: Vec<MovingHeatScheduleGradient>,
    pub initial_specific_enthalpy: Vec<f64>,
    pub transpose_iterations: usize,
    pub max_transpose_relative_residual: f64,
}

/// Evaluate fitting, point limits and sampled-peak limits on ONE forward tape.
/// Each peak needs last_endpoint+1 transpose calls, not one sweep per sample.
/// Point and peak rows share the response/gradient/transpose allowance; the
/// independent max_peak_samples bounds the sum of all selected sample counts.
/// No partial result escapes failed sampling, a derivative gate or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_scan_objective_with_peaks(
    cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>, endpoints: &[f64],
    observations: &[ScanObservation], sampling: ScanObjectiveSampling,
    step: EnthalpyStepConfig, transpose: LinearConfig, budget: ScanObjectiveBudget,
    responses: &[ScanResponse], peaks: &[ScanPeakResponse], response_budget: ScanResponseBudget,
    max_peak_samples: usize,
) -> Result<(ScanObjectiveEvaluation, Vec<ScanResponseEvaluation>, Vec<ScanPeakEvaluation>), ScanSensitivityError> {
    admit_rows(cx, &run, endpoints.len(), responses, peaks, response_budget, max_peak_samples)?;
    let (objective, (points, peaks)) = evaluate_with_tape(cx, run, endpoints, observations,
        sampling, step, transpose, budget, |run, tape| {
            let mut points = Vec::with_capacity(responses.len());
            for &response in responses {
                poll(cx, points.len())?;
                points.push(pull_response(cx, run, tape, response, sampling, transpose, budget.max_schedule_controls)?);
            }
            let mut output = Vec::with_capacity(peaks.len());
            for peak in peaks {
                poll(cx, output.len())?;
                output.push(pull_peak(cx, run, tape, peak, sampling, transpose, budget.max_schedule_controls)?);
            }
            Ok((points, output))
        })?;
    Ok((objective, points, peaks))
}

// Same physical response/work admission for full and checkpointed tapes.
pub(crate) fn admit_rows(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    endpoint_count: usize, responses: &[ScanResponse], peaks: &[ScanPeakResponse],
    response_budget: ScanResponseBudget, max_peak_samples: usize,
) -> Result<Vec<usize>, ScanSensitivityError> {
    poll(cx, 0)?;
    let count = responses.len().checked_add(peaks.len())
        .filter(|&n| n <= response_budget.max_responses)
        .ok_or_else(|| invalid("combined scan point/peak response allowance exceeded"))?;
    let vertices = run.conduction().temperature().len();
    let mut entries = vertices;
    for (i, source) in run.sources().iter().enumerate() {
        poll(cx, i)?;
        entries = entries.checked_add(source.schedule_control_count()?)
            .ok_or_else(|| invalid("scan peak gradient shape overflow"))?;
    }
    if entries.checked_mul(count).is_none_or(|n| n > response_budget.max_gradient_entries) {
        return Err(invalid("combined scan point/peak gradients exceed their retained scalar allowance"));
    }
    let mut solves = endpoint_count;
    for (i, response) in responses.iter().enumerate() {
        poll(cx, i)?;
        if response.endpoint >= endpoint_count || response.vertex >= vertices {
            return Err(invalid("scan point response must name an existing endpoint and vertex"));
        }
        solves = solves.checked_add(response.endpoint + 1)
            .ok_or_else(|| invalid("scan point response work overflow"))?;
    }
    let mut samples = 0usize;
    let mut sample_counts = Vec::with_capacity(peaks.len());
    for peak in peaks {
        let count = peak.validate(cx, endpoint_count, vertices, max_peak_samples - samples)?;
        samples += count;
        sample_counts.push(count);
        solves = solves.checked_add(peak.last_endpoint + 1)
            .ok_or_else(|| invalid("scan peak response work overflow"))?;
    }
    if solves > response_budget.max_transpose_solves {
        return Err(invalid("combined objective/point/peak sweeps exceed their transpose allowance"));
    }
    Ok(sample_counts)
}

fn sample(binding: &ScanStepLinearization<'_>, observable: ScanObservable, vertex: usize) -> f64 {
    let solution = binding.report().solution();
    match observable {
        ScanObservable::Temperature => solution.temperature[vertex],
        ScanObservable::SpecificEnthalpy => solution.specific_enthalpy_j_kg[vertex],
        ScanObservable::LiquidMassFraction => solution.liquid_mass_fraction[vertex],
    }
}

// A finite ordered difference may overflow toward -infinity. Its exponential
// is zero; do not square, exponentiate positive magnitudes or invent clipping.
pub(crate) fn weight(value: f64, maximum: f64, scale: f64) -> f64 {
    let argument = (value - maximum) / scale;
    if argument == f64::NEG_INFINITY { 0.0 } else { fs_math::det::exp(argument) }
}

pub(crate) struct PeakSummary {
    pub maximum: f64,
    pub sum: f64,
    pub gap: f64,
    pub value: f64,
    pub bound: f64,
}

// Shared reduction: checkpoint samples and full-tape samples use the same
// summation order, normalization and binary64 envelope, not two approximations.
pub(crate) fn summarize_samples(cx: &Cx<'_>, count: usize, scale: f64,
    sample: impl Fn(usize) -> f64,
) -> Result<PeakSummary, ScanSensitivityError> {
    if count == 0 || !scale.is_finite() || scale <= 0.0 {
        return Err(invalid("scan peak reduction needs samples and a positive finite scale"));
    }
    let mut maximum = f64::NEG_INFINITY;
    for i in 0..count {
        if i % 256 == 0 { poll(cx, i)?; }
        maximum = maximum.max(finite(sample(i))?);
    }
    let mut sum = 0.0;
    for i in 0..count {
        if i % 256 == 0 { poll(cx, i)?; }
        sum += weight(sample(i), maximum, scale);
    }
    if !sum.is_finite() || sum < 1.0 { return Err(invalid("scan peak normalization is not representable")); }
    let gap = finite(scale * fs_math::det::ln(sum))?;
    let value = finite(maximum + gap)?;
    let bound = finite(scale * fs_math::det::ln(count as f64))?;
    poll(cx, count)?;
    Ok(PeakSummary { maximum, sum, gap, value, bound })
}

fn pull_peak(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    tape: &[ScanStepLinearization<'_>], peak: &ScanPeakResponse,
    sampling: ScanObjectiveSampling, transpose: LinearConfig, max_controls: usize,
) -> Result<ScanPeakEvaluation, ScanSensitivityError> {
    let sample_count = (peak.last_endpoint - peak.first_endpoint + 1) * peak.vertices.len();
    let PeakSummary { maximum, sum, gap, value, bound } = summarize_samples(cx,
        sample_count, peak.smoothing_scale, |i| sample(
            &tape[peak.first_endpoint+i/peak.vertices.len()], peak.observable,
            peak.vertices[i%peak.vertices.len()]))?;
    let mut seed = zeros(run.conduction().temperature().len())?;
    let mut accumulated: Option<Vec<MovingHeatScheduleGradient>> = None;
    let mut iterations = 0usize;
    let mut residual = 0.0_f64;
    for index in (0..=peak.last_endpoint).rev() {
        poll(cx, index)?;
        let binding = &tape[index];
        if index >= peak.first_endpoint {
            let mut local = zeros(seed.len())?;
            for (i, &vertex) in peak.vertices.iter().enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                local[vertex] = weight(sample(binding, peak.observable, vertex), maximum, peak.smoothing_scale) / sum;
            }
            local = peak.observable.seed_in_enthalpy(cx, binding, local)?;
            // Earlier downstream derivatives are already in enthalpy coordinates.
            // Convert only the NEW local observable seed, never the whole seed.
            add(cx, &mut seed, &local)?;
        }
        let gradient = binding.pullback(cx, &seed, transpose)?;
        iterations = iterations.checked_add(gradient.storage.transport.iterations)
            .ok_or_else(|| invalid("scan peak transpose iteration count overflow"))?;
        residual = residual.max(gradient.storage.transport.relative_residual);
        if gradient.sources.len() != run.sources().len() { return Err(invalid("scan peak source layout changed")); }
        let mut schedules = Vec::with_capacity(gradient.sources.len());
        for (i, (source, gradient)) in run.sources().iter().zip(gradient.sources).enumerate() {
            poll(cx, i)?;
            schedules.push(match sampling {
                ScanObjectiveSampling::RightEndpoint => source.pullback_schedule_at(cx, binding.end_time_s(), gradient, max_controls)?,
                ScanObjectiveSampling::EventLeft => source.pullback_schedule_before(cx, binding.end_time_s(), gradient, max_controls)?,
            });
        }
        match &mut accumulated {
            None => accumulated = Some(schedules),
            Some(total) => {
                if total.len() != schedules.len() { return Err(invalid("scan peak source identity changed")); }
                for (a, b) in total.iter_mut().zip(&schedules) {
                    for axis in 0..3 { add(cx, &mut a.center_values_m[axis], &b.center_values_m[axis])?; }
                    add(cx, &mut a.half_width_m, &b.half_width_m)?;
                    add(cx, &mut a.absorbed_power_values_w, &b.absorbed_power_values_w)?;
                }
            }
        }
        seed = gradient.storage.transport.previous_specific_enthalpy;
    }
    poll(cx, sample_count)?;
    Ok(ScanPeakEvaluation { peak: peak.clone(), value, sampled_maximum: maximum,
        smoothing_gap: gap, maximum_smoothing_gap: bound, sample_count,
        sources: accumulated.ok_or_else(|| invalid("missing scan peak sweep"))?,
        initial_specific_enthalpy: seed, transpose_iterations: iterations,
        max_transpose_relative_residual: residual })
}
