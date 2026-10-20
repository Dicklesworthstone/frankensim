//! Joint absorbed-power, path and footprint fitting through native scan adjoints.
//!
//! The factory builds the same physical model from supplied source declarations
//! for every candidate. Its returned sources are checked against those exact
//! declarations; the factory must hold mesh, materials, history and boundaries
//! fixed. No approximate forward model, finite-difference gradient, clipping,
//! or altered primal acceptance policy is used by the optimizer adapter.

use fs_conduction::bc::history::march::enthalpy::scanning::ScanningEnthalpyRun;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::{ScanSensitivityError,
    objective::{evaluate_scan_objective, ScanObjectiveEvaluation, ScanObjectiveBudget,
        ScanObjectiveSampling, ScanObservation}};
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::responses::{
    ScanResponse, ScanResponseBudget, ScanResponseEvaluation,
    peaks::{evaluate_scan_objective_with_peaks, ScanPeakResponse, ScanPeakEvaluation}};
use fs_conduction::bc::history::march::moving::sensitivity::schedule::MovingHeatScheduleGradient;
use fs_conduction::bc::history::march::moving::MovingHeatSource;
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint::{
    ScanCheckpointBudget, ScanCheckpointWork, evaluate_scan_objective_checkpointed_with_peaks};
use fs_conduction::transient::enthalpy::EnthalpyStepConfig;
use fs_conduction::{ConductionError, LinearConfig};
use fs_exec::Cx;
use crate::{LbfgsState, LbfgsReport, LbfgsError, StopRule};
use std::collections::BTreeSet;

mod geometry;
pub use geometry::ScanGeometryControl;
/// Physical bounds and absorbed-energy constraints through the existing SQP engine.
pub mod constraints;

/// One selected authored power value; unselected values (including off states)
/// remain exactly as declared. Coordinates obey P = reference_w * exp(x).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanPowerControl {
    pub source: usize,
    pub value: usize,
    pub reference_w: f64,
}

/// Original fixed-grid objective and independent physical/adjoint allowances.
pub struct ScanPowerObjective {
    pub endpoints_s: Vec<f64>,
    pub observations: Vec<ScanObservation>,
    pub sampling: ScanObjectiveSampling,
    pub step: EnthalpyStepConfig,
    pub transpose: LinearConfig,
    pub budget: ScanObjectiveBudget,
    /// Instant used only to validate new source declarations, not as a time shift.
    pub source_reference_time_s: f64,
}

/// One native objective with the complete optimizer-coordinate chain rule.
#[derive(Debug)]
pub struct ScanPowerEvaluation {
    pub physical: ScanObjectiveEvaluation,
    /// Log-power derivatives followed by geometry derivatives in declaration order.
    pub gradient: Vec<f64>,
    /// Selected powers in control declaration order, W.
    pub powers_w: Vec<f64>,
    /// Selected center values and half widths in geometry-control order, metres.
    pub geometry_m: Vec<f64>,
    /// Actual original/recomputed endpoint work when short-tape mode is selected.
    /// None means the existing full-storage path was used.
    pub checkpoint_work: Option<ScanCheckpointWork>,
}

/// Bounded inverse-design evaluator over immutable template sources.
///
/// F is a caller-owned constructor of the physical run, not an objective or
/// gradient callback. Each accepted/rejected optimizer probe starts from its
/// original state. The study's evaluation counter is independent of optimizer
/// checkpoints; cloning or restarting L-BFGS cannot refund physical work.
/// No uniqueness, global optimum, or experimentally calibrated result is claimed.
pub struct ScanPowerStudy<F> {
    sources: Vec<MovingHeatSource>,
    controls: Vec<ScanPowerControl>,
    geometry_controls: Vec<ScanGeometryControl>,
    objective: ScanPowerObjective,
    factory: F,
    initial: Vec<f64>,
    evaluations: usize,
    evaluation_limit: usize,
    checkpointing: Option<ScanCheckpointBudget>,
}

impl<F> ScanPowerStudy<F> {
    /// Bind controls without a primal solve. Selected initial values and their
    /// reference scales must be positive; leave exact zero/off entries unselected.
    pub fn new(cx: &Cx<'_>, sources: Vec<MovingHeatSource>, controls: Vec<ScanPowerControl>,
        objective: ScanPowerObjective, factory: F, evaluation_limit: usize)
        -> Result<Self, ScanSensitivityError>
    {
        Self::new_joint(cx, sources, controls, Vec::new(), objective, factory, evaluation_limit)
    }

    /// Optimize selected powers, authored path values and footprint widths in
    /// one physical study. Either control group may be empty, but not both.
    /// Coordinates are log powers first, followed by scaled signed centers or
    /// log widths in geometry declaration order. No event time is a variable.
    /// Original source support, motion, material and derivative refusals remain
    /// authoritative; an invalid design is not silently projected or penalized.
    #[allow(clippy::too_many_arguments)]
    pub fn new_joint(cx: &Cx<'_>, sources: Vec<MovingHeatSource>, controls: Vec<ScanPowerControl>,
        geometry_controls: Vec<ScanGeometryControl>, objective: ScanPowerObjective,
        factory: F, evaluation_limit: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx)?;
        let count = controls.len().checked_add(geometry_controls.len())
            .filter(|&n| n > 0 && n <= objective.budget.max_schedule_controls)
            .ok_or_else(|| invalid("scan study needs a nonempty bounded set of power or geometry controls"))?;
        if sources.is_empty()
            || !objective.source_reference_time_s.is_finite()
        { return Err(invalid("power study needs sources, bounded nonempty controls and an explicit finite reference time")); }
        let mut total_controls = 0usize;
        for source in &sources {
            poll(cx)?;
            total_controls = total_controls.checked_add(source.schedule_control_count()?)
                .filter(|&n| n <= objective.budget.max_schedule_controls)
                .ok_or_else(|| invalid("power study sources exceed their combined authored-control allowance"))?;
        }
        if objective.endpoints_s.is_empty() || objective.endpoints_s.len() > objective.budget.max_steps
            || objective.observations.is_empty() || objective.observations.len() > objective.budget.max_observations
        { return Err(invalid("power study needs bounded nonempty endpoints and observations")); }
        let mut seen = BTreeSet::new();
        let mut initial = Vec::with_capacity(count);
        for control in &controls {
            poll(cx)?;
            let source = sources.get(control.source).ok_or_else(|| invalid("power control names an absent source"))?;
            let values = source.absorbed_power_values(objective.budget.max_schedule_controls)?;
            let value = *values.get(control.value).ok_or_else(|| invalid("power control names an absent authored value"))?;
            if !seen.insert((control.source,control.value)) || !value.is_finite() || value <= 0.0
                || !control.reference_w.is_finite() || control.reference_w <= 0.0
            { return Err(invalid("power controls must be distinct with positive initial and reference watts")); }
            // Separate logarithms avoid overflow/underflow in value/reference.
            initial.push(finite(fs_math::det::ln(value)-fs_math::det::ln(control.reference_w))?);
        }
        let mut geometry_seen = BTreeSet::new();
        for &control in &geometry_controls {
            poll(cx)?;
            let source = sources.get(control.source()).ok_or_else(|| invalid("geometry control names an absent source"))?;
            if !geometry_seen.insert((control.source(),control.key())) {
                return Err(invalid("geometry controls must name distinct authored coordinates"));
            }
            initial.push(control.initial(source.geometry_value(control.key())?)?);
        }
        Ok(Self { sources, controls, geometry_controls, objective, factory, initial, evaluations:0,
            evaluation_limit, checkpointing: None })
    }

    pub fn initial_coordinates(&self) -> &[f64] { &self.initial }
    pub const fn evaluations(&self) -> usize { self.evaluations }
    pub const fn evaluation_limit(&self) -> usize { self.evaluation_limit }

    /// Select short-tape differentiation before the first physical evaluation.
    /// Applies to unconstrained fitting AND point/peak-constrained SQP without
    /// discarding any state requirements. The objective's retained-tangent caps
    /// bound one segment; its endpoint cap still bounds the whole trajectory.
    /// Checkpoint arrays and the original-plus-replayed solve work have separate
    /// explicit allowances. Mesh-dependent shapes are checked by the evaluator.
    ///
    /// The policy cannot change after a successful or failed physical probe, so
    /// an optimizer checkpoint cannot silently resume under a different policy.
    pub fn set_checkpointing(&mut self, cx: &Cx<'_>, policy: ScanCheckpointBudget)
        -> Result<(), ScanSensitivityError>
    {
        poll(cx)?;
        if self.evaluations != 0 || policy.segment_steps == 0 || policy.max_checkpoint_vertices == 0
            || self.objective.endpoints_s.len().checked_mul(2)
                .is_none_or(|n| n > policy.max_forward_solves)
        { return Err(invalid("select positive checkpoint allowances before any physical probe, including both original and replayed endpoints")); }
        self.checkpointing = Some(policy);
        Ok(())
    }

    pub const fn checkpointing(&self) -> Option<ScanCheckpointBudget> { self.checkpointing }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), ScanSensitivityError> {
        if total <= self.evaluation_limit { return Err(invalid("scan evaluation allowance must increase")); }
        self.evaluation_limit = total;
        Ok(())
    }
    /// Decode selected positive watts without running or mutating the model.
    pub fn powers_at(&self, point: &[f64]) -> Result<Vec<f64>, ScanSensitivityError> {
        self.check_point(point)?;
        point.iter().zip(&self.controls).map(|(&x,c)| positive_coordinate(x,c.reference_w)).collect()
    }

    /// Decode selected path/width values without a physical solve. The complete
    /// point uses power coordinates first, then geometry coordinates.
    pub fn geometry_at(&self, point: &[f64]) -> Result<Vec<f64>, ScanSensitivityError> {
        self.check_point(point)?;
        point[self.controls.len()..].iter().zip(&self.geometry_controls)
            .map(|(&x,&c)| c.decode(x)).collect()
    }

    fn check_point(&self, point: &[f64]) -> Result<(), ScanSensitivityError> {
        if point.len() != self.initial.len() { return Err(invalid("scan design point has the wrong shape")); }
        for &value in point { finite(value)?; }
        Ok(())
    }

    /// Export reusable physical source declarations at an optimizer point.
    /// No trajectory is solved and no physical-evaluation slot is consumed.
    /// All unselected geometry/power values, knots and interpolation are
    /// preserved by the same builder used for optimizer probes.
    pub fn sources_at(&self, cx: &Cx<'_>, point: &[f64])
        -> Result<Vec<MovingHeatSource>, ScanSensitivityError>
    {
        poll(cx)?;
        self.rebuild_sources(cx, &self.powers_at(point)?, &self.geometry_at(point)?)
    }

    fn rebuild_sources(&self, cx: &Cx<'_>, powers: &[f64], geometry: &[f64])
        -> Result<Vec<MovingHeatSource>, ScanSensitivityError>
    {
        let limit = self.objective.budget.max_schedule_controls;
        let mut values: Vec<Option<Vec<f64>>> = vec![None; self.sources.len()];
        for (control, &power) in self.controls.iter().zip(powers) {
            poll(cx)?;
            if values[control.source].is_none() {
                values[control.source] = Some(self.sources[control.source].absorbed_power_values(limit)?);
            }
            values[control.source].as_mut().expect("power values prepared")[control.value] = power;
        }
        let mut edits = vec![Vec::new(); self.sources.len()];
        for (&control, &value) in self.geometry_controls.iter().zip(geometry) {
            poll(cx)?;
            edits[control.source()].push((control.key(), value));
        }
        let mut sources = Vec::with_capacity(self.sources.len());
        for (index, source) in self.sources.iter().enumerate() {
            poll(cx)?;
            let powered = match &values[index] {
                Some(values) => source.with_absorbed_power_values(values, self.objective.source_reference_time_s, limit)?,
                None => source.clone(),
            };
            sources.push(if edits[index].is_empty() { powered } else {
                powered.with_geometry_values(cx, &edits[index], self.objective.source_reference_time_s, limit)?
            });
        }
        poll(cx)?;
        Ok(sources)
    }

    /// Evaluate actual physics and the full history adjoint once. Failed probes
    /// consume a reserved evaluation slot; they do not change an accepted
    /// optimizer point. The factory must not vary non-source model data.
    pub fn evaluate<'h,'m,'c>(&mut self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanPowerEvaluation, ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>, ScanSensitivityError>
    {
        self.evaluate_with_responses(cx, point, None).map(|(objective, _, _)| objective)
    }

    // Constraint rows and the fitting objective share the exact supplied-source
    // check, one original trajectory, and one cumulative evaluation reservation.
    fn evaluate_with_responses<'h,'m,'c>(&mut self, cx: &Cx<'_>, point: &[f64],
        responses: Option<(&[ScanResponse], &[ScanPeakResponse], ScanResponseBudget, usize)>)
        -> Result<(ScanPowerEvaluation, Vec<ScanResponseEvaluation>, Vec<ScanPeakEvaluation>), ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>, ScanSensitivityError>
    {
        poll(cx)?;
        let powers_w = self.powers_at(point)?;
        let geometry_m = self.geometry_at(point)?;
        if self.evaluations >= self.evaluation_limit { return Err(invalid("scan inverse evaluation budget exhausted")); }
        self.evaluations += 1;
        let sources = self.rebuild_sources(cx, &powers_w, &geometry_m)?;
        let run = (self.factory)(cx,sources.clone())?;
        if run.conduction().accepted_steps() != 0 {
            return Err(invalid("scan factory must return an unadvanced physical run"));
        }
        if run.sources() != sources.as_slice() { return Err(invalid("scan factory replaced the supplied controlled sources")); }
        let mut checkpoint_work = None;
        let (physical, rows, peaks) = if let Some(checkpoints) = self.checkpointing {
            let (rows, peaks, response_budget, max_samples) = responses.unwrap_or((
                &[], &[], ScanResponseBudget { max_responses: 0, max_gradient_entries: 0,
                    max_transpose_solves: self.objective.endpoints_s.len() }, 0));
            let result = evaluate_scan_objective_checkpointed_with_peaks(cx, run,
                &self.objective.endpoints_s, &self.objective.observations, self.objective.sampling,
                self.objective.step, self.objective.transpose, self.objective.budget,
                rows, peaks, response_budget, max_samples, checkpoints)?;
            checkpoint_work = Some(result.work);
            (result.objective, result.responses, result.peaks)
        } else { match responses {
            None => (evaluate_scan_objective(cx,run,&self.objective.endpoints_s,&self.objective.observations,
                self.objective.sampling,self.objective.step,self.objective.transpose,self.objective.budget)?, Vec::new(), Vec::new()),
            Some((responses, peaks, response_budget, max_samples)) => evaluate_scan_objective_with_peaks(cx,run,
                &self.objective.endpoints_s,&self.objective.observations,self.objective.sampling,
                self.objective.step,self.objective.transpose,self.objective.budget,responses,peaks,response_budget,max_samples)?,
        }};
        let gradient = self.chain_gradient(cx, &powers_w, &geometry_m, &physical.sources)?;
        poll(cx)?;
        Ok((ScanPowerEvaluation { physical,gradient,powers_w,geometry_m,checkpoint_work }, rows, peaks))
    }

    // Both objective and unsquared state-response derivatives use this same
    // positive-power, signed-position and positive-width coordinate chain.
    fn chain_gradient(&self, cx: &Cx<'_>, powers_w: &[f64], geometry_m: &[f64],
        sources: &[MovingHeatScheduleGradient]) -> Result<Vec<f64>, ScanSensitivityError>
    {
        let mut gradient = Vec::with_capacity(self.initial.len());
        for (control,&power) in self.controls.iter().zip(powers_w) {
            poll(cx)?;
            let g = sources.get(control.source)
                .and_then(|s|s.absorbed_power_values_w.get(control.value))
                .ok_or_else(||invalid("native authored-power gradient changed its control layout"))?;
            gradient.push(finite(power * *g)?);
        }
        for (&control, &value) in self.geometry_controls.iter().zip(geometry_m) {
            poll(cx)?;
            let source_gradient = sources.get(control.source())
                .ok_or_else(||invalid("native geometry gradient changed its source layout"))?;
            gradient.push(control.chain(value, source_gradient)?);
        }
        poll(cx)?;
        Ok(gradient)
    }

    /// Initialize the existing fallible L-BFGS engine on these exact coordinates.
    pub fn begin_optimizer<'h,'m,'c>(&mut self, cx: &Cx<'_>, memory: usize)
        -> Result<LbfgsState,LbfgsError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>, ScanSensitivityError>
    {
        let initial = self.initial.clone();
        LbfgsState::try_new(&initial,memory,&mut |point| {
            self.evaluate(cx,point).map(|r|(r.physical.value,r.gradient))
        })
    }

    /// Continue this study's L-BFGS checkpoint under both native and cumulative
    /// physical-evaluation limits. Typed source/phase/solver errors propagate;
    /// they are never replaced by an invented penalty or a widened physical gate.
    pub fn advance_optimizer<'h,'m,'c>(&mut self, cx: &Cx<'_>, state: &mut LbfgsState,
        rule: &StopRule, max_iterations: usize) -> Result<LbfgsReport,LbfgsError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>, ScanSensitivityError>
    {
        poll(cx).map_err(LbfgsError::Evaluation)?;
        let remaining = self.evaluation_limit.saturating_sub(self.evaluations);
        let ceiling = state.evals.checked_add(remaining)
            .ok_or(LbfgsError::InvalidInput("scan optimizer evaluation ceiling overflow"))?;
        state.try_run(&mut |point|self.evaluate(cx,point).map(|r|(r.physical.value,r.gradient)),
            rule,max_iterations,ceiling)
    }
}

fn positive_coordinate(x: f64, reference: f64) -> Result<f64, ScanSensitivityError> {
    finite(x)?;
    let value = finite(fs_math::det::exp(finite(x+fs_math::det::ln(reference))?))?;
    if value <= 0.0 { return Err(invalid("positive scan coordinate underflowed to zero")); }
    Ok(value)
}
fn finite(value: f64) -> Result<f64, ScanSensitivityError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("scan inverse coordinate or gradient is not representable")) }
}
fn invalid(what: &str) -> ScanSensitivityError {
    ConductionError::Config { parameter:"scan inverse design",what:what.into() }.into()
}
fn poll(cx: &Cx<'_>) -> Result<(),ScanSensitivityError> {
    cx.checkpoint().map_err(|_|ConductionError::Cancelled {stage:"scan-inverse-design",at:0}.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn positive_coordinates_preserve_extreme_reference_scales_and_the_chain_rule() {
        for reference in [1e-200_f64,1.,1e200] {
            for x in [-2.,0.,2.] {
                let p=positive_coordinate(x,reference).unwrap();
                assert!((p/(reference*x.exp())-1.).abs()<2e-13);
                let eps=1e-4;
                let fd=(positive_coordinate(x+eps,reference).unwrap()-positive_coordinate(x-eps,reference).unwrap())/(2.*eps);
                assert!((fd/p-1.).abs()<1e-7);
            }
        }
        assert!(positive_coordinate(f64::NAN,1.).is_err());
        assert!(positive_coordinate(1000.,1.).is_err());
        assert!(positive_coordinate(-1000.,1.).is_err());
    }
}
