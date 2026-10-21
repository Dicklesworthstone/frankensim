//! Physical scan scores without requiring a classical endpoint derivative.
use super::{ScanPowerStudy, ScanObservationSpec, ScanSensitivityError,
    ScanObservation, ScanningEnthalpyRun, MovingHeatSource, invalid, poll};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::{
    ScanValueBudget, ScanValueEvaluation};
pub use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::responses::{
    ScanValueResponseBudget, ScanResponseValue, ScanPeakValue};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::responses::
    evaluate_scan_values_with_responses;
use super::{ScanResponse, ScanPeakResponse};
use fs_exec::Cx;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Same model and source controls, but intentionally no derivative/KKT fields.
#[derive(Debug)]
pub struct ScanPowerValueEvaluation {
    pub physical: ScanValueEvaluation,
    /// Actual selected absorbed powers, in watts, without snapping or clipping.
    pub powers_w: Vec<f64>,
    /// Actual selected path coordinates/half widths in metres.
    pub geometry_m: Vec<f64>,
}

/// Scoped, thread-shareable access to one fixed study's value-only producer.
///
/// This borrows the study exclusively, so declarations and lifetime allowances
/// cannot change while workers are evaluating. It shares only immutable inputs:
/// every invocation builds a fresh physical run using its caller's own Cx.
/// Sharing across workers requires a `Fn + Sync` factory and `Sync` readouts;
/// an arbitrary `FnMut` factory is never cloned, locked or shared implicitly.
///
/// The full allowance is reserved up front. Drop returns only unstarted slots,
/// including after errors or unwinding. Forgetting the guard conservatively
/// spends the full reservation, never refunds work. Callback-owned physical
/// memory and deterministic factory behavior remain the caller's obligations.
pub struct ScanValueBatch<'a, F, O = ScanObservation> {
    study: &'a mut ScanPowerStudy<F, O>,
    started: AtomicUsize,
    allowance: usize,
    previously_spent: usize,
}

impl<F, O> Drop for ScanValueBatch<'_, F, O> {
    fn drop(&mut self) {
        // The exclusive borrow cannot end while a worker still borrows this
        // guard. Every successful reservation is <= allowance by construction.
        self.study.evaluations -= self.allowance - self.started.load(Ordering::Relaxed);
    }
}

impl<F, O: ScanObservationSpec> ScanValueBatch<'_, F, O> {
    /// Started physical reservations, including failed rebuilds and solves.
    pub fn started_evaluations(&self) -> usize { self.started.load(Ordering::Relaxed) }
    pub const fn evaluation_allowance(&self) -> usize { self.allowance }
    /// Total study work including reservations actually started in this batch.
    pub fn evaluations(&self) -> usize { self.previously_spent + self.started_evaluations() }

    pub(crate) fn study(&self) -> &ScanPowerStudy<F, O> { self.study }

    pub fn evaluate_value<'h, 'm, 'c>(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanPowerValueEvaluation, ScanSensitivityError>
    where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        self.evaluate_value_with_responses(cx, point, &[], &[],
            ScanValueResponseBudget { max_responses: 0, max_peak_samples: 0 })
            .map(|(objective, _, _)| objective)
    }

    /// Same producer, response reductions and phase-corner semantics as serial
    /// evaluation. Each started reservation is charged once even on failure.
    pub fn evaluate_value_with_responses<'h, 'm, 'c>(&self, cx: &Cx<'_>, point: &[f64],
        responses: &[ScanResponse], peaks: &[ScanPeakResponse], budget: ScanValueResponseBudget)
        -> Result<(ScanPowerValueEvaluation, Vec<ScanResponseValue>, Vec<ScanPeakValue>), ScanSensitivityError>
    where F: Fn(&Cx<'_>, Vec<MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        poll(cx)?;
        let powers_w = self.study.powers_at(point)?;
        let geometry_m = self.study.geometry_at(point)?;
        self.started.fetch_update(Ordering::Relaxed, Ordering::Relaxed,
            |used| (used < self.allowance).then(|| used + 1))
            .map_err(|_| invalid("scan value batch evaluation allowance exhausted"))?;
        let sources = self.study.rebuild_sources(cx, &powers_w, &geometry_m)?;
        let run = (self.study.factory)(cx, sources.clone())?;
        self.study.finish_value_with_responses(cx, run, &sources, powers_w, geometry_m,
            responses, peaks, budget)
    }
}

impl<F, O: ScanObservationSpec> ScanPowerStudy<F, O> {
    /// Borrow this immutable physical declaration for a bounded evaluation
    /// batch. No model is run. Zero allowance is valid for draining cached work.
    /// All started work returns to this study's existing cumulative counter;
    /// unused reservations return when the batch is dropped after worker drain.
    pub fn value_batch(&mut self, cx: &Cx<'_>, allowance: usize)
        -> Result<ScanValueBatch<'_, F, O>, ScanSensitivityError>
    {
        poll(cx)?;
        let reserved = self.evaluations.checked_add(allowance)
            .filter(|&total| total <= self.evaluation_limit)
            .ok_or_else(|| invalid("scan value batch exceeds the remaining physical allowance"))?;
        let previously_spent = self.evaluations;
        self.evaluations = reserved;
        Ok(ScanValueBatch { study: self, started: AtomicUsize::new(0), allowance, previously_spent })
    }

    /// Evaluate the original physical fitting loss once without any adjoint.
    ///
    /// Valid primal states at a phase/material corner remain evaluable. Failed
    /// physical solves still propagate; no NaN or penalty replaces their error.
    /// The same cumulative evaluation counter is shared with gradient probes,
    /// including failed probes. A later derivative probe retains all its original
    /// smoothness checks. Checkpointing is unnecessary here: no tape is retained.
    ///
    /// This evaluates only the fitting objective, not a ConstrainedScanStudy's
    /// inequalities. A scalar score must not be interpreted as feasibility.
    pub fn evaluate_value<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanPowerValueEvaluation, ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<super::MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        self.evaluate_value_with_responses(cx, point, &[], &[],
            ScanValueResponseBudget { max_responses: 0, max_peak_samples: 0 })
            .map(|(objective, _, _)| objective)
    }

    /// Sample the objective and physical point/peak requirements with one
    /// value-only trajectory and one cumulative evaluation reservation.
    /// Response values are never softened by a robust fitting loss. No adjoint
    /// tape or smoothness assumption is introduced, including at phase corners.
    pub fn evaluate_value_with_responses<'h, 'm, 'c>(&mut self, cx: &Cx<'_>, point: &[f64],
        responses: &[ScanResponse], peaks: &[ScanPeakResponse], budget: ScanValueResponseBudget)
        -> Result<(ScanPowerValueEvaluation, Vec<ScanResponseValue>, Vec<ScanPeakValue>), ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<super::MovingHeatSource>)
        -> Result<ScanningEnthalpyRun<'h, 'm, 'c>, ScanSensitivityError>
    {
        poll(cx)?;
        let powers_w = self.powers_at(point)?;
        let geometry_m = self.geometry_at(point)?;
        if self.evaluations >= self.evaluation_limit {
            return Err(invalid("scan inverse evaluation budget exhausted"));
        }
        self.evaluations += 1;
        let sources = self.rebuild_sources(cx, &powers_w, &geometry_m)?;
        let run = (self.factory)(cx, sources.clone())?;
        self.finish_value_with_responses(cx, run, &sources, powers_w, geometry_m,
            responses, peaks, budget)
    }

    // Serial and shared batches use the same source-identity checks and original
    // forward/observation owner. Only booking and calling the factory differ.
    #[allow(clippy::too_many_arguments)]
    fn finish_value_with_responses(&self, cx: &Cx<'_>, run: ScanningEnthalpyRun<'_, '_, '_>,
        sources: &[MovingHeatSource], powers_w: Vec<f64>, geometry_m: Vec<f64>,
        responses: &[ScanResponse], peaks: &[ScanPeakResponse], budget: ScanValueResponseBudget)
        -> Result<(ScanPowerValueEvaluation, Vec<ScanResponseValue>, Vec<ScanPeakValue>), ScanSensitivityError>
    {
        if run.conduction().accepted_steps() != 0 {
            return Err(invalid("scan factory must return an unadvanced physical run"));
        }
        if run.sources() != sources {
            return Err(invalid("scan factory replaced the supplied controlled sources"));
        }
        let (physical, responses, peaks) = evaluate_scan_values_with_responses(cx, run, &self.objective.endpoints_s,
            &self.objective.observations, self.objective.sampling, self.objective.step,
            ScanValueBudget { max_steps: self.objective.budget.max_steps,
                max_observations: self.objective.budget.max_observations,
                max_schedule_controls: self.objective.budget.max_schedule_controls },
            responses, peaks, budget)?;
        poll(cx)?;
        Ok((ScanPowerValueEvaluation { physical, powers_w, geometry_m }, responses, peaks))
    }
}
