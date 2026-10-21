//! Physical scan scores without requiring a classical endpoint derivative.
use super::{ScanPowerStudy, ScanObservationSpec, ScanSensitivityError,
    ScanningEnthalpyRun, invalid, poll};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::{
    ScanValueBudget, ScanValueEvaluation};
pub use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::responses::{
    ScanValueResponseBudget, ScanResponseValue, ScanPeakValue};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::value::responses::
    evaluate_scan_values_with_responses;
use super::{ScanResponse, ScanPeakResponse};
use fs_exec::Cx;

/// Same model and source controls, but intentionally no derivative/KKT fields.
#[derive(Debug)]
pub struct ScanPowerValueEvaluation {
    pub physical: ScanValueEvaluation,
    /// Actual selected absorbed powers, in watts, without snapping or clipping.
    pub powers_w: Vec<f64>,
    /// Actual selected path coordinates/half widths in metres.
    pub geometry_m: Vec<f64>,
}

impl<F, O: ScanObservationSpec> ScanPowerStudy<F, O> {
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
        if run.conduction().accepted_steps() != 0 {
            return Err(invalid("scan factory must return an unadvanced physical run"));
        }
        if run.sources() != sources.as_slice() {
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
