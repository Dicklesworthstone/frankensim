//! Check real design feasibility without requiring physical derivatives.
//!
//! This is an evaluator, not a penalty objective or a constrained optimizer.
//! A low measurement loss never substitutes for satisfying the returned rows.
//! The original SQP input admission/energy weights and physical value collector
//! remain the owners; no second thermal model or constraint convention is used.

use super::{ScanPowerStudy, ScanObservation, ScanObservationSpec, ScanSensitivityError,
    ScanControlBound, ScanEnergyCap, ScanStateBound, ScanPeakBound, ScanResponseBudget,
    ScanConstraintEvaluation, ScanningEnthalpyRun, MovingHeatSource,
    Limits, LimitAdmission, finite, invalid, poll};
use super::super::value::{ScanPowerValueEvaluation, ScanValueResponseBudget};
use fs_exec::Cx;

/// Explicit retained scalar/selection limits, with no irrelevant KKT allowance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanValueConstraintBudget {
    /// All inequality rows, including automatic nonnegative-domain rows.
    pub max_rows: usize,
    /// Energy windows times selected power values. Fixed energy is one scalar
    /// per window; no row-by-coordinate derivative matrix is allocated.
    pub max_energy_weights: usize,
    pub responses: ScanValueResponseBudget,
}

/// Physical values and signed inequalities, intentionally without Jacobians.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanConstraintValues {
    pub nonnegative_power_values_w: Vec<f64>,
    pub bound_values: Vec<f64>,
    pub absorbed_energy_j: Vec<f64>,
    pub state_values: Vec<f64>,
    pub peak_values: Vec<f64>,
    pub sampled_peak_values: Vec<f64>,
    pub peak_smoothing_gaps: Vec<f64>,
    /// Same order/scales as SQP: nonnegative powers, control lower/upper,
    /// energy caps, state lower/upper, then peak caps. Each must be <= 0.
    pub residuals: Vec<f64>,
}

impl ScanConstraintValues {
    fn from_inputs(input: ScanConstraintEvaluation) -> Self {
        debug_assert!(input.jacobian.is_empty());
        Self { nonnegative_power_values_w: input.nonnegative_power_values_w,
            bound_values: input.bound_values, absorbed_energy_j: input.absorbed_energy_j,
            state_values: input.state_values, peak_values: input.peak_values,
            sampled_peak_values: input.sampled_peak_values,
            peak_smoothing_gaps: input.peak_smoothing_gaps, residuals: input.residuals }
    }

    /// Largest positive scaled inequality. This is not fitting loss or a
    /// physical error estimate. Zero means every represented row is satisfied.
    pub fn max_scaled_violation(&self) -> f64 {
        self.residuals.iter().map(|v| v.max(0.0)).fold(0.0, f64::max)
    }

    /// Check all declared rows with an explicit nonnegative scaled tolerance.
    /// Tolerance permits tolerance*row.scale physical violation; it does not
    /// include solver error or provide coverage of unselected times/vertices.
    pub fn is_feasible(&self, tolerance: f64) -> Result<bool, ScanSensitivityError> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(invalid("scan feasibility tolerance must be finite and nonnegative"));
        }
        for &value in &self.residuals { finite(value)?; }
        Ok(self.max_scaled_violation() <= tolerance)
    }
}

/// Immutable physical constraints around the original source-controlled model.
/// Unlike an SQP constructor, admission performs no physical solve or adjoint,
/// so initial states on a phase corner are admitted. Every later evaluation
/// shares the study's lifetime reservation, including failed physical attempts.
pub struct ValueConstrainedScanStudy<F, O = ScanObservation> {
    study: ScanPowerStudy<F, O>,
    limits: Limits,
}

impl<F, O: ScanObservationSpec> ValueConstrainedScanStudy<F, O> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(cx: &Cx<'_>, study: ScanPowerStudy<F, O>, bounds: Vec<ScanControlBound>,
        energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>, peaks: Vec<ScanPeakBound>,
        budget: ScanValueConstraintBudget) -> Result<Self, ScanSensitivityError>
    {
        let limits = Limits::compile(cx, &study, bounds, energy, states, peaks,
            ScanResponseBudget { max_responses: budget.responses.max_responses,
                max_gradient_entries: 0, max_transpose_solves: 0 },
            budget.responses.max_peak_samples, LimitAdmission::Values {
                max_rows: budget.max_rows, max_energy_weights: budget.max_energy_weights })?;
        Ok(Self { study, limits })
    }

    pub fn initial_coordinates(&self) -> &[f64] { self.study.initial_coordinates() }
    pub fn evaluations(&self) -> usize { self.study.evaluations() }
    pub fn evaluation_limit(&self) -> usize { self.study.evaluation_limit() }
    pub fn state_bounds(&self) -> &[ScanStateBound] { &self.limits.states }
    pub fn peak_bounds(&self) -> &[ScanPeakBound] { &self.limits.peaks }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), ScanSensitivityError> {
        self.study.extend_evaluation_budget(total)
    }

    /// Export the original candidate declarations without asserting feasibility.
    pub fn sources_at(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<Vec<MovingHeatSource>, ScanSensitivityError>
    { self.study.sources_at(cx, point) }

    /// Cheap input-only check. Never silently omit attached state/peak limits.
    pub fn constraints_at(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanConstraintValues, ScanSensitivityError>
    {
        poll(cx)?;
        if !self.limits.states.is_empty() || !self.limits.peaks.is_empty() {
            return Err(invalid("state constraints require evaluate_at on the physical trajectory"));
        }
        Ok(ScanConstraintValues::from_inputs(self.limits.at_impl(cx, &self.study, point, false)?))
    }

    /// Evaluate objective and every constraint on one original trajectory.
    /// No dense Jacobian, transpose solve, replay or phase derivative is needed.
    /// A successful but infeasible simulation returns its signed violations;
    /// failed physics remains an error and never becomes a fake penalty value.
    pub fn evaluate_at<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, point: &[f64])
        -> Result<(ScanPowerValueEvaluation, ScanConstraintValues), ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        let mut out = ScanConstraintValues::from_inputs(self.limits.at_impl(cx, &self.study, point, false)?);
        let requests: Vec<_> = self.limits.states.iter().map(|state| state.response).collect();
        let peaks: Vec<_> = self.limits.peaks.iter().map(|peak| peak.response.clone()).collect();
        let (physical, responses, peak_responses) = self.study.evaluate_value_with_responses(cx, point,
            &requests, &peaks, ScanValueResponseBudget {
                max_responses: self.limits.response_budget.max_responses,
                max_peak_samples: self.limits.max_peak_samples })?;
        if responses.len() != self.limits.states.len() || peak_responses.len() != self.limits.peaks.len() {
            return Err(invalid("value-only state response layout changed"));
        }
        for (bound, response) in self.limits.states.iter().zip(responses) {
            poll(cx)?;
            if response.response != bound.response { return Err(invalid("value-only response identity changed")); }
            out.state_values.push(response.value);
            for (limit, sign) in [(bound.lower, -1.0), (bound.upper, 1.0)] {
                if let Some(limit) = limit {
                    out.residuals.push(finite(sign * finite(response.value - limit)? / bound.scale)?);
                }
            }
        }
        for (bound, response) in self.limits.peaks.iter().zip(peak_responses) {
            poll(cx)?;
            if response.peak != bound.response { return Err(invalid("value-only peak identity changed")); }
            out.peak_values.push(response.value);
            out.sampled_peak_values.push(response.sampled_maximum);
            out.peak_smoothing_gaps.push(response.smoothing_gap);
            out.residuals.push(finite(finite(response.value - bound.maximum)? / bound.scale)?);
        }
        if out.residuals.len() != self.limits.rows { return Err(invalid("value-only constraint row count changed")); }
        poll(cx)?;
        Ok((physical, out))
    }
}
