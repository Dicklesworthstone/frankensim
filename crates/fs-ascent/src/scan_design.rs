//! Positive absorbed-power schedule fitting through native scan adjoints.
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
use fs_conduction::bc::history::march::moving::MovingHeatSource;
use fs_conduction::transient::enthalpy::EnthalpyStepConfig;
use fs_conduction::{ConductionError, LinearConfig};
use fs_exec::Cx;
use crate::{LbfgsState, LbfgsReport, LbfgsError, StopRule};
use std::collections::BTreeSet;

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

/// One native objective with the exact positive-coordinate chain rule.
#[derive(Debug)]
pub struct ScanPowerEvaluation {
    pub physical: ScanObjectiveEvaluation,
    pub gradient: Vec<f64>,
    /// Selected powers in control declaration order, W.
    pub powers_w: Vec<f64>,
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
    objective: ScanPowerObjective,
    factory: F,
    initial: Vec<f64>,
    evaluations: usize,
    evaluation_limit: usize,
}

impl<F> ScanPowerStudy<F> {
    /// Bind controls without a primal solve. Selected initial values and their
    /// reference scales must be positive; leave exact zero/off entries unselected.
    pub fn new(cx: &Cx<'_>, sources: Vec<MovingHeatSource>, controls: Vec<ScanPowerControl>,
        objective: ScanPowerObjective, factory: F, evaluation_limit: usize)
        -> Result<Self, ScanSensitivityError>
    {
        poll(cx)?;
        if sources.is_empty() || controls.is_empty() || controls.len() > objective.budget.max_schedule_controls
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
        let mut initial = Vec::with_capacity(controls.len());
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
        Ok(Self { sources, controls, objective, factory, initial, evaluations:0, evaluation_limit })
    }

    pub fn initial_coordinates(&self) -> &[f64] { &self.initial }
    pub const fn evaluations(&self) -> usize { self.evaluations }
    pub const fn evaluation_limit(&self) -> usize { self.evaluation_limit }
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), ScanSensitivityError> {
        if total <= self.evaluation_limit { return Err(invalid("scan evaluation allowance must increase")); }
        self.evaluation_limit = total;
        Ok(())
    }
    /// Decode selected positive watts without running or mutating the model.
    pub fn powers_at(&self, point: &[f64]) -> Result<Vec<f64>, ScanSensitivityError> {
        if point.len() != self.controls.len() { return Err(invalid("log-power point has the wrong shape")); }
        point.iter().zip(&self.controls).map(|(&x,c)| positive_power(x,c.reference_w)).collect()
    }

    /// Evaluate actual physics and the full history adjoint once. Failed probes
    /// consume a reserved evaluation slot; they do not change an accepted
    /// optimizer point. The factory must not vary non-source model data.
    pub fn evaluate<'h,'m,'c>(&mut self, cx: &Cx<'_>, point: &[f64])
        -> Result<ScanPowerEvaluation, ScanSensitivityError>
    where F: FnMut(&Cx<'_>, Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>, ScanSensitivityError>
    {
        poll(cx)?;
        let powers_w = self.powers_at(point)?;
        if self.evaluations >= self.evaluation_limit { return Err(invalid("scan inverse evaluation budget exhausted")); }
        self.evaluations += 1;
        let limit = self.objective.budget.max_schedule_controls;
        let mut values = Vec::with_capacity(self.sources.len());
        for source in &self.sources { poll(cx)?; values.push(source.absorbed_power_values(limit)?); }
        for (control,&power) in self.controls.iter().zip(&powers_w) { values[control.source][control.value] = power; }
        let mut sources = Vec::with_capacity(self.sources.len());
        for (source,values) in self.sources.iter().zip(&values) {
            poll(cx)?;
            sources.push(source.with_absorbed_power_values(values,self.objective.source_reference_time_s,limit)?);
        }
        let run = (self.factory)(cx,sources.clone())?;
        if run.conduction().accepted_steps() != 0 {
            return Err(invalid("scan factory must return an unadvanced physical run"));
        }
        if run.sources() != sources.as_slice() { return Err(invalid("scan factory replaced the supplied controlled sources")); }
        let physical = evaluate_scan_objective(cx,run,&self.objective.endpoints_s,&self.objective.observations,
            self.objective.sampling,self.objective.step,self.objective.transpose,self.objective.budget)?;
        let mut gradient = Vec::with_capacity(self.controls.len());
        for (control,&power) in self.controls.iter().zip(&powers_w) {
            poll(cx)?;
            let g = physical.sources.get(control.source)
                .and_then(|s|s.absorbed_power_values_w.get(control.value))
                .ok_or_else(||invalid("native authored-power gradient changed its control layout"))?;
            gradient.push(finite(power * *g)?);
        }
        poll(cx)?;
        Ok(ScanPowerEvaluation { physical,gradient,powers_w })
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

fn positive_power(x: f64, reference: f64) -> Result<f64, ScanSensitivityError> {
    finite(x)?;
    let value = finite(fs_math::det::exp(finite(x+fs_math::det::ln(reference))?))?;
    if value <= 0.0 { return Err(invalid("log-power coordinate underflowed to zero")); }
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
                let p=positive_power(x,reference).unwrap();
                assert!((p/(reference*x.exp())-1.).abs()<2e-13);
                let eps=1e-4;
                let fd=(positive_power(x+eps,reference).unwrap()-positive_power(x-eps,reference).unwrap())/(2.*eps);
                assert!((fd/p-1.).abs()<1e-7);
            }
        }
        assert!(positive_power(f64::NAN,1.).is_err());
        assert!(positive_power(1000.,1.).is_err());
        assert!(positive_power(-1000.,1.).is_err());
    }
}
