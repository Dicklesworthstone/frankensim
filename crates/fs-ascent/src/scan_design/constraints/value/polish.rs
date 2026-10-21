//! Switch a value-only physical search to the original constrained adjoints.
//!
//! The session borrows, rather than clones or rebuilds, its physical study.
//! All source controls, measurements, constraints and spent work stay attached.
//! A derivative refusal leaves the value-only owner available after the borrow
//! ends. No one-sided derivative or simulation penalty is substituted.

use super::{ValueConstrainedScanStudy, ScanConstraintValues, ScanObservation,
    ScanObservationSpec, ScanSensitivityError, ScanResponseBudget, Limits,
    MovingHeatSource, ScanningEnthalpyRun, ScanPowerValueEvaluation, poll};
use crate::{SqpError, SqpRunReport, SqpSample, SqpState};
use fs_exec::Cx;

/// Resumable local refinement borrowing the original physical evaluator.
/// The caller must provision the original objective's adjoint-storage budgets
/// before starting the study. Value-only population evaluations do not spend
/// those storage allowances or build a derivative tape.
pub struct ScanSqpPolish<'a, F, O = ScanObservation> {
    owner: &'a mut ValueConstrainedScanStudy<F, O>,
    limits: Limits,
    boxes: Vec<[f64; 2]>,
    state: SqpState,
}

/// Value-only replay of the accepted SQP point, independently of the adjoint.
/// KKT status belongs to the separate SqpRunReport, not to a successful replay.
#[derive(Debug)]
pub struct ScanSqpReplay {
    pub coordinates: Vec<f64>,
    pub physical: ScanPowerValueEvaluation,
    pub constraints: ScanConstraintValues,
    pub sources: Vec<MovingHeatSource>,
    /// Additional lower/upper inequalities in optimizer-coordinate units.
    /// The original physical rows remain unchanged and separately reported.
    pub box_residuals: Vec<f64>,
    pub physical_evaluations: usize,
}

impl<F, O: ScanObservationSpec> ValueConstrainedScanStudy<F, O> {
    /// Start native SQP at an already discovered design. This costs one actual
    /// trajectory plus its admitted adjoints; one physical slot remains reserved
    /// for replay. Failure retains all spent work in this same owner.
    ///
    /// Nonempty boxes must cover every optimizer coordinate. They remain hard
    /// evaluation-domain bounds and explicit KKT inequalities, so a population
    /// search box cannot disappear during polishing. Empty means no extra box.
    /// The original nonnegative-power and physical-domain checks still apply.
    ///
    /// Constitutive corners can refuse this derivative path even when the
    /// original value-only evaluation succeeds. No accepted population state
    /// needs to be discarded: end this borrow and continue that search.
    pub fn begin_sqp_at<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, point: &[f64],
        coordinate_bounds: &[[f64; 2]], response_budget: ScanResponseBudget,
        max_kkt_dimension: usize)
        -> Result<ScanSqpPolish<'_, F, O>, SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        poll(cx).map_err(SqpError::Evaluation)?;
        if self.evaluation_limit().saturating_sub(self.evaluations()) < 2 {
            return Err(SqpError::Invalid("SQP handoff needs an initial physical evaluation and one replay"));
        }
        let n = self.initial_coordinates().len();
        admit_boxes(point, coordinate_bounds, n)?;
        let extra_rows = coordinate_bounds.len().checked_mul(2)
            .ok_or(SqpError::Invalid("SQP box row count overflow"))?;
        let original_cap = max_kkt_dimension.checked_sub(extra_rows)
            .ok_or(SqpError::Invalid("SQP box rows exceed the dense dimension allowance"))?;
        // Reuse original bounds/energy admission and exact row order. Only the
        // derivative-work allowances change; this does not alter owner.limits.
        let limits = Limits::new(cx, &self.study, self.limits.bounds.clone(),
            self.limits.energy.iter().map(|row| row.cap).collect(),
            self.limits.states.clone(), self.limits.peaks.clone(), response_budget,
            self.limits.max_peak_samples, original_cap).map_err(SqpError::Evaluation)?;
        let boxes = coordinate_bounds.to_vec();
        let state = SqpState::try_new(point, max_kkt_dimension, &mut |point| {
            sample(cx, self, &limits, &boxes, point)
        }, Some(cx))?;
        Ok(ScanSqpPolish { owner: self, limits, boxes, state })
    }
}

impl<F, O: ScanObservationSpec> ScanSqpPolish<'_, F, O> {
    pub fn point(&self) -> &[f64] { self.state.point() }
    pub fn sample(&self) -> &SqpSample { self.state.sample() }
    pub fn evaluations(&self) -> usize { self.state.evaluations() }
    pub fn physical_evaluations(&self) -> usize { self.owner.evaluations() }
    pub fn physical_evaluation_limit(&self) -> usize { self.owner.evaluation_limit() }
    pub fn original_constraint_count(&self) -> usize { self.limits.rows }

    /// Extending physical work never refunds earlier population or SQP probes.
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), ScanSensitivityError> {
        self.owner.extend_evaluation_budget(total)
    }

    /// Continue the same accepted SQP/BFGS state under an additional callback
    /// ceiling AND the original remaining physical allowance, excluding replay.
    /// Out-of-box and negative-power trials may spend callback attempts without
    /// running physics. Original solve/derivative failures always propagate.
    pub fn advance<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, tolerance: f64,
        additional_iterations: usize, additional_evaluations: usize)
        -> Result<SqpRunReport, SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        let remaining = self.owner.evaluation_limit().saturating_sub(self.owner.evaluations())
            .saturating_sub(1).min(additional_evaluations);
        let ceiling = self.state.evaluations().checked_add(remaining)
            .ok_or(SqpError::Invalid("SQP polish callback ceiling overflow"))?;
        let Self { owner, limits, boxes, state } = self;
        state.try_run(&mut |point| sample(cx, owner, limits, boxes, point),
            tolerance, additional_iterations, ceiling, Some(cx))
    }

    /// Run the accepted design through the original value-only owner. Require
    /// bit-identical objective and EVERY inequality, not merely the largest
    /// violation. A failed check publishes no replay or source declaration.
    pub fn replay<'h, 'm, 'c>(&mut self, cx: &Cx<'_>)
        -> Result<ScanSqpReplay, SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        poll(cx).map_err(SqpError::Evaluation)?;
        let coordinates = self.state.point().to_vec();
        let (physical, constraints) = self.owner.evaluate_at(cx, &coordinates)
            .map_err(SqpError::Evaluation)?;
        let box_residuals = box_rows(&coordinates, &self.boxes);
        let sample = self.state.sample();
        let rows = constraints.residuals.iter().chain(&box_residuals);
        if physical.physical.value.to_bits() != sample.f.to_bits()
            || constraints.residuals.len() + box_residuals.len() != sample.ci.len()
            || rows.zip(&sample.ci).any(|(a, b)| a.to_bits() != b.to_bits())
        { return Err(SqpError::Invalid("SQP physical replay changed objective or constraint values")); }
        let sources = self.owner.sources_at(cx, &coordinates).map_err(SqpError::Evaluation)?;
        poll(cx).map_err(SqpError::Evaluation)?;
        Ok(ScanSqpReplay { coordinates, physical, constraints, sources, box_residuals,
            physical_evaluations: self.owner.evaluations() })
    }
}

fn admit_boxes<E>(point: &[f64], bounds: &[[f64; 2]], dimension: usize) -> Result<(), SqpError<E>> {
    if point.len() != dimension || point.iter().any(|x| !x.is_finite())
        || (!bounds.is_empty() && bounds.len() != dimension)
    { return Err(SqpError::Invalid("SQP handoff point and box must match the original coordinates")); }
    for (&x, &[lo, hi]) in point.iter().zip(bounds) {
        if !lo.is_finite() || !hi.is_finite() || lo >= hi || !(hi-lo).is_finite() || x < lo || x > hi {
            return Err(SqpError::Invalid("SQP handoff needs finite ordered boxes containing its start"));
        }
    }
    Ok(())
}

fn box_rows(point: &[f64], boxes: &[[f64; 2]]) -> Vec<f64> {
    point.iter().zip(boxes).flat_map(|(&x, &[lo, hi])| [lo-x, x-hi]).collect()
}

fn sample<'h, 'm, 'c, F, O: ScanObservationSpec>(cx: &Cx<'_>,
    owner: &mut ValueConstrainedScanStudy<F, O>, limits: &Limits,
    boxes: &[[f64; 2]], point: &[f64]) -> Result<Option<SqpSample>, ScanSensitivityError>
where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>)
    -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
{
    poll(cx)?;
    // Only this declared analytic domain is unavailable. Never swallow an
    // original physical failure or a failed derivative as an out-of-box trial.
    if point.iter().zip(boxes).any(|(&x, &[lo, hi])| !x.is_finite() || x < lo || x > hi) {
        return Ok(None);
    }
    let Some(mut out) = limits.sample(cx, &mut owner.study, point)? else { return Ok(None); };
    let n = point.len();
    for (i, &[lo, hi]) in boxes.iter().enumerate() {
        poll(cx)?;
        for (residual, derivative) in [(lo-point[i], -1.0), (point[i]-hi, 1.0)] {
            out.ci.push(residual);
            let offset = out.ji.len();
            out.ji.resize(offset+n, 0.0);
            out.ji[offset+i] = derivative;
        }
    }
    Ok(Some(out))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn box_admission_preserves_coordinate_domain_and_signed_rows() {
        assert!(admit_boxes::<()>( &[-1.0, 2.0], &[[-1.0, 0.0], [1.0, 2.0]], 2).is_ok());
        assert_eq!(box_rows(&[-1.0, 2.0], &[[-1.0, 0.0], [1.0, 2.0]]), vec![0.0, -1.0, -1.0, 0.0]);
        assert!(admit_boxes::<()>(&[0.0], &[], 1).is_ok());
        for bounds in [vec![[0.0, 0.0]], vec![[1.0, 2.0]], vec![[0.0, f64::INFINITY]],
            vec![[0.0, 1.0], [0.0, 1.0]]] {
            assert!(admit_boxes::<()>(&[0.0], &bounds, 1).is_err());
        }
    }
}
