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
use super::super::value::{ScanPowerValueEvaluation, ScanValueResponseBudget,
    ScanValueBatch, ScanResponseValue, ScanPeakValue};
use fs_exec::Cx;

/// Shared authored designs evaluated under multiple physical operating cases.
pub mod scenarios;

/// Borrow the same physical study for constrained adjoint-based local refinement.
pub mod polish;

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
    /// Signed W/s in source-declaration/segment order, including fixed ramps.
    pub power_slew_values_w_s: Vec<f64>,
    pub state_values: Vec<f64>,
    pub peak_values: Vec<f64>,
    pub sampled_peak_values: Vec<f64>,
    pub peak_smoothing_gaps: Vec<f64>,
    /// Same order/scales as SQP: nonnegative powers, control lower/upper,
    /// energy caps, slew rise/fall, state lower/upper, then peak caps. Each must be <= 0.
    pub residuals: Vec<f64>,
}

impl ScanConstraintValues {
    fn from_inputs(input: ScanConstraintEvaluation) -> Self {
        debug_assert!(input.jacobian.is_empty());
        Self { nonnegative_power_values_w: input.nonnegative_power_values_w,
            bound_values: input.bound_values, absorbed_energy_j: input.absorbed_energy_j,
            power_slew_values_w_s: input.power_slew_values_w_s,
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

/// Borrow one immutable objective/constraint declaration for concurrent value
/// evaluations. Workers receive their own Cx and independent physical runs.
/// The underlying batch reserves and settles the original lifetime budget;
/// no worker can replace source controls, constraints or feasibility scales.
pub struct ScanConstraintValueBatch<'a, F, O = ScanObservation> {
    values: ScanValueBatch<'a, F, O>,
    limits: &'a Limits,
}

impl<F, O: ScanObservationSpec> ScanConstraintValueBatch<'_, F, O> {
    pub fn started_evaluations(&self) -> usize { self.values.started_evaluations() }
    pub fn evaluations(&self) -> usize { self.values.evaluations() }
    pub fn evaluation_allowance(&self) -> usize { self.values.evaluation_allowance() }

    /// Evaluate every original inequality, without a derivative or penalty.
    /// F must be callable through an immutable reference; sharing the batch
    /// between threads additionally requires F and its observation layout Sync.
    pub fn evaluate_at<'h, 'm, 'c>(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<(ScanPowerValueEvaluation, ScanConstraintValues), ScanSensitivityError>
    where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        let out = ScanConstraintValues::from_inputs(
            self.limits.at_impl(cx, self.values.study(), point, false)?);
        let requests: Vec<_> = self.limits.states.iter().map(|state| state.response).collect();
        let peaks: Vec<_> = self.limits.peaks.iter().map(|peak| peak.response.clone()).collect();
        let (physical, responses, peak_responses) = self.values.evaluate_value_with_responses(cx, point,
            &requests, &peaks, ScanValueResponseBudget {
                max_responses: self.limits.response_budget.max_responses,
                max_peak_samples: self.limits.max_peak_samples })?;
        Ok((physical, self.limits.finish_values(cx, out, &responses, &peak_responses)?))
    }
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

    /// Add rise/fall bounds on every segment of each named continuous power
    /// schedule. Fixed and unselected values participate, including shutdown.
    /// Held tables refuse: timestamp differences do not bound their jumps.
    /// No physical solve is run during admission. All later serial/batch/scenario
    /// evaluations return these rows through the same original value collector.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_power_slew_bounds(cx: &Cx<'_>, study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>,
        peaks: Vec<ScanPeakBound>, slew_bounds: Vec<super::ScanPowerSlewBound>,
        max_slew_segments: usize, budget: ScanValueConstraintBudget)
        -> Result<Self, ScanSensitivityError>
    {
        let limits = Limits::compile_with_slew(cx, &study, bounds, energy, states, peaks,
            ScanResponseBudget { max_responses: budget.responses.max_responses,
                max_gradient_entries: 0, max_transpose_solves: 0 },
            budget.responses.max_peak_samples, slew_bounds, max_slew_segments,
            LimitAdmission::Values { max_rows: budget.max_rows,
                max_energy_weights: budget.max_energy_weights })?;
        Ok(Self { study, limits })
    }

    /// Exact immutable inequality count, including automatic domain rows.
    /// Available before running physics for bounded population admission.
    pub fn constraint_count(&self) -> usize { self.limits.rows }

    pub fn initial_coordinates(&self) -> &[f64] { self.study.initial_coordinates() }
    pub fn evaluations(&self) -> usize { self.study.evaluations() }
    pub fn evaluation_limit(&self) -> usize { self.study.evaluation_limit() }
    pub fn state_bounds(&self) -> &[ScanStateBound] { &self.limits.states }
    pub fn peak_bounds(&self) -> &[ScanPeakBound] { &self.limits.peaks }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), ScanSensitivityError> {
        self.study.extend_evaluation_budget(total)
    }

    /// Import spent work when rebuilding a fresh study from an application
    /// checkpoint. The caller must first bind the original model and limits.
    /// This charges work only: it cannot reset an evaluated study, increase an
    /// allowance or undo failed/replayed evaluations. It runs no physics.
    pub fn restore_evaluation_count(&mut self, cx: &Cx<'_>, spent: usize)
        -> Result<(), ScanSensitivityError>
    {
        poll(cx)?;
        if self.study.evaluations != 0 || spent > self.study.evaluation_limit {
            return Err(invalid("restore spent physical work only into an unevaluated study within its original allowance"));
        }
        self.study.evaluations = spent;
        Ok(())
    }

    /// Admit bounded shared evaluation without cloning factories or constraints.
    /// Drop after draining workers restores only unused reservation slots;
    /// started failed/cancelled work remains charged to this same study.
    pub fn value_batch(&mut self, cx: &Cx<'_>, allowance: usize)
        -> Result<ScanConstraintValueBatch<'_, F, O>, ScanSensitivityError>
    {
        let Self { study, limits } = self;
        Ok(ScanConstraintValueBatch { values: study.value_batch(cx, allowance)?, limits })
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
        let out = ScanConstraintValues::from_inputs(self.limits.at_impl(cx, &self.study, point, false)?);
        let requests: Vec<_> = self.limits.states.iter().map(|state| state.response).collect();
        let peaks: Vec<_> = self.limits.peaks.iter().map(|peak| peak.response.clone()).collect();
        let (physical, responses, peak_responses) = self.study.evaluate_value_with_responses(cx, point,
            &requests, &peaks, ScanValueResponseBudget {
                max_responses: self.limits.response_budget.max_responses,
                max_peak_samples: self.limits.max_peak_samples })?;
        Ok((physical, self.limits.finish_values(cx, out, &responses, &peak_responses)?))
    }
}

impl Limits {
    // Keep row ordering, physical scales and smoothing reports identical for
    // ordinary FnMut serial studies and immutable concurrently shared batches.
    fn finish_values(&self, cx: &Cx<'_>, mut out: ScanConstraintValues,
        responses: &[ScanResponseValue], peak_responses: &[ScanPeakValue])
        -> Result<ScanConstraintValues, ScanSensitivityError>
    {
        if responses.len() != self.states.len() || peak_responses.len() != self.peaks.len() {
            return Err(invalid("value-only state response layout changed"));
        }
        for (bound, response) in self.states.iter().zip(responses) {
            poll(cx)?;
            if response.response != bound.response { return Err(invalid("value-only response identity changed")); }
            out.state_values.push(response.value);
            for (limit, sign) in [(bound.lower, -1.0), (bound.upper, 1.0)] {
                if let Some(limit) = limit {
                    out.residuals.push(finite(sign * finite(response.value - limit)? / bound.scale)?);
                }
            }
        }
        for (bound, response) in self.peaks.iter().zip(peak_responses) {
            poll(cx)?;
            if response.peak != bound.response { return Err(invalid("value-only peak identity changed")); }
            out.peak_values.push(response.value);
            out.sampled_peak_values.push(response.sampled_maximum);
            out.peak_smoothing_gaps.push(response.smoothing_gap);
            out.residuals.push(finite(finite(response.value - bound.maximum)? / bound.scale)?);
        }
        if out.residuals.len() != self.rows { return Err(invalid("value-only constraint row count changed")); }
        poll(cx)?;
        Ok(out)
    }
}
