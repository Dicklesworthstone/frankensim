//! Finite-scenario worst-case design over the existing equilibrium and SQP owners.
//!
//! Minimize t subject to J_s(x + delta_s) <= t and EVERY physical response
//! constraint in EVERY supplied scenario. The extra epigraph coordinate avoids
//! inventing a smooth gradient for max(J_s) at a tie. No probability, interval
//! coverage, unseen-scenario guarantee, new contact law or optimizer is implied.
use super::{sample as physical_sample, DesignControl, DesignError,
    DesignEvaluation, DesignWork, EquilibriumDesign, SqpError, SqpRunReport, SqpSample,
    SqpState, SqpStop};
use fs_exec::CancelGate;

/// One fixed realization of additive parameter tolerances. A shared variable's
/// offset follows ALL its bound fields, exactly like its nominal value.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumScenario {
    /// Unique nonempty label, at most 128 bytes.
    pub name: String,
    /// One offset per design variable, in that variable's original physical
    /// units. Zero denotes no perturbation. No weights/probabilities are inferred.
    pub physical_offsets: Vec<f64>,
}

/// Preserve the original case error AND the scenario that caused it.
#[derive(Debug)]
pub enum ScenarioError {
    /// Malformed scenario data or nonrepresentable arithmetic.
    Invalid(&'static str),
    /// A source design refusal. None names the unperturbed nominal parameter
    /// decode; Some(i) names scenario i in construction order.
    Design { scenario: Option<usize>, source: DesignError },
    /// Cancellation before returning a complete scenario family.
    Cancelled,
}
impl core::fmt::Display for ScenarioError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "equilibrium scenarios: {what}"),
            Self::Design { scenario, source } => write!(f, "equilibrium scenario {scenario:?}: {source}"),
            Self::Cancelled => write!(f, "equilibrium scenario evaluation cancelled"),
        }
    }
}
impl std::error::Error for ScenarioError {}
impl ScenarioError {
    fn domain_rejection(&self) -> bool {
        matches!(self, Self::Design { source: DesignError::OutsideBounds { .. }, .. })
    }
}
fn poll(gate: &CancelGate) -> Result<(), ScenarioError> {
    if gate.is_requested() { Err(ScenarioError::Cancelled) } else { Ok(()) }
}
fn design_error(scenario: Option<usize>, source: DesignError) -> ScenarioError {
    if matches!(&source, DesignError::Cancelled) { ScenarioError::Cancelled }
    else { ScenarioError::Design { scenario, source } }
}

/// A complete scenario-family evaluation. The nominal parameter vector is NOT
/// an extra physical solve unless the caller supplies a zero-offset scenario.
#[derive(Clone, Debug, PartialEq)]
pub struct ScenarioEvaluation {
    /// Unperturbed physical design values in variable declaration order.
    pub nominal_parameters: Vec<f64>,
    /// Largest actual scenario objective, not the optimizer's epigraph estimate.
    pub worst_objective: f64,
    /// Complete original physical results, in scenario declaration order.
    pub scenarios: Vec<DesignEvaluation>,
}

/// Immutable finite uncertainty set over one immutable physical problem.
/// Domain intersection keeps both nominal and realized parameters within the
/// original bounds. There is no clipping, silent scenario removal or resampling.
#[derive(Clone)]
pub struct ScenarioProblem<'a> {
    problem: &'a EquilibriumDesign,
    scenarios: Vec<EquilibriumScenario>,
    shifts: Vec<Vec<f64>>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    maximum_kkt_dimension: usize,
}
impl<'a> ScenarioProblem<'a> {
    /// Admit 1..=32 declared scenarios and a bounded dense SQP system before
    /// evaluating physics. The KKT screen is 3*n + 1 + s*(1+c), including the
    /// epigraph, common box faces, s objective rows and all c physical rows.
    /// Zero-width or empty common parameter domains refuse explicitly.
    pub fn new(problem: &'a EquilibriumDesign, scenarios: Vec<EquilibriumScenario>,
        maximum_scenarios: usize, maximum_kkt_dimension: usize, gate: &CancelGate)
        -> Result<Self, ScenarioError>
    {
        poll(gate)?;
        let n = problem.variables().len();
        if n == 0 || !(1..=32).contains(&maximum_scenarios) || scenarios.is_empty()
            || scenarios.len() > maximum_scenarios || !(1..=1024).contains(&maximum_kkt_dimension) {
            return Err(ScenarioError::Invalid("nonempty bounded scenario and decision families are required"));
        }
        let rows = scenarios.len().checked_mul(problem.constraints().len() + 1)
            .and_then(|v| n.checked_mul(3).and_then(|w| w.checked_add(1)).and_then(|w| v.checked_add(w)))
            .ok_or(ScenarioError::Invalid("scenario KKT size overflow"))?;
        if rows > maximum_kkt_dimension {
            return Err(ScenarioError::Invalid("all scenario, physical and box rows must fit the dense KKT cap"));
        }
        let mut lower = Vec::with_capacity(n);
        let mut upper = Vec::with_capacity(n);
        for variable in problem.variables() {
            let lo = (variable.minimum - variable.reference) / variable.scale;
            let hi = (variable.maximum - variable.reference) / variable.scale;
            if !lo.is_finite() || !hi.is_finite() || !((hi-lo).is_finite() && lo < hi) {
                return Err(ScenarioError::Invalid("parameter domain is not representable in decision coordinates"));
            }
            lower.push(lo); upper.push(hi);
        }
        let mut shifts = Vec::with_capacity(scenarios.len());
        for (i, scenario) in scenarios.iter().enumerate() {
            poll(gate)?;
            if scenario.name.trim().is_empty() || scenario.name.len() > 128
                || scenarios[..i].iter().any(|s| s.name == scenario.name)
                || scenario.physical_offsets.len() != n {
                return Err(ScenarioError::Invalid("each scenario needs a unique label and one offset per variable"));
            }
            let mut shift = Vec::with_capacity(n);
            for (j, (offset, variable)) in scenario.physical_offsets.iter().zip(problem.variables()).enumerate() {
                let delta = offset / variable.scale;
                let lo = (variable.minimum - variable.reference) / variable.scale - delta;
                let hi = (variable.maximum - variable.reference) / variable.scale - delta;
                if !offset.is_finite() || !delta.is_finite() || (*offset != 0.0 && delta == 0.0)
                    || !lo.is_finite() || !hi.is_finite() {
                    return Err(ScenarioError::Invalid("physical tolerance offset is not representable"));
                }
                lower[j] = lower[j].max(lo); upper[j] = upper[j].min(hi);
                shift.push(delta);
            }
            shifts.push(shift);
        }
        if lower.iter().zip(&upper).any(|(lo, hi)| lo >= hi) {
            return Err(ScenarioError::Invalid("scenarios have no common positive-width nominal parameter domain"));
        }
        Ok(Self { problem, scenarios, shifts, lower, upper, maximum_kkt_dimension })
    }
    /// Original physical model, response constraints and declared variable units.
    #[must_use]
    pub fn problem(&self) -> &EquilibriumDesign { self.problem }
    /// Exact supplied scenario names and physical offsets.
    #[must_use]
    pub fn scenarios(&self) -> &[EquilibriumScenario] { &self.scenarios }
    /// Common nominal bounds in the source problem's scaled decision units.
    #[must_use]
    pub fn decision_bounds(&self) -> (&[f64], &[f64]) { (&self.lower, &self.upper) }

    /// Solve all independent load cases at every scenario. Every original
    /// evaluate call consumes the caller's DesignControl; failed work remains
    /// charged. All shifted domains are checked BEFORE the first physical solve.
    /// Insufficient mid-family budget returns no partial family, without refunds.
    pub fn evaluate(&self, nominal: &[f64], control: &mut DesignControl, gate: &CancelGate)
        -> Result<ScenarioEvaluation, ScenarioError>
    {
        poll(gate)?;
        let nominal_parameters = self.problem.physical_parameters(nominal)
            .map_err(|e| design_error(None, e))?;
        let mut points = Vec::with_capacity(self.scenarios.len());
        for (i, shift) in self.shifts.iter().enumerate() {
            poll(gate)?;
            let point: Vec<f64> = nominal.iter().zip(shift).map(|(x, d)| if *d == 0.0 { *x } else { x+d }).collect();
            self.problem.physical_parameters(&point).map_err(|e| design_error(Some(i), e))?;
            points.push(point);
        }
        let mut scenarios = Vec::with_capacity(points.len());
        let mut worst_objective = f64::NEG_INFINITY;
        for (i, point) in points.iter().enumerate() {
            poll(gate)?;
            let value = self.problem.evaluate(point, control, gate).map_err(|e| design_error(Some(i), e))?;
            worst_objective = worst_objective.max(value.value);
            scenarios.push(value);
        }
        poll(gate)?;
        Ok(ScenarioEvaluation { nominal_parameters, worst_objective, scenarios })
    }

    fn sample(&self, point: &[f64], evaluation: &ScenarioEvaluation) -> SqpSample {
        let n = self.lower.len();
        let width = n+1;
        let mut gradient = vec![0.0; width]; gradient[n] = 1.0;
        let mut ci = Vec::new(); let mut ji = Vec::new();
        let mut ce = Vec::new(); let mut je = Vec::new();
        for i in 0..n {
            ci.push(self.lower[i]-point[i]); ci.push(point[i]-self.upper[i]);
            let mut row = vec![0.0; width]; row[i] = -1.0; ji.extend_from_slice(&row);
            row[i] = 1.0; ji.extend_from_slice(&row);
        }
        // Within each inequality scenario: epigraph row, then original physical
        // inequalities. Equalities retain scenario order and source order. The
        // existing adapter owns physical-row signs, normalizations and Jacobians.
        for physical in &evaluation.scenarios {
            ci.push(physical.value-point[n]);
            ji.extend_from_slice(&physical.gradient); ji.push(-1.0);
            let original = physical_sample(self.problem, physical);
            ci.extend_from_slice(&original.ci[2*n..]);
            for row in original.ji[2*n*n..].chunks_exact(n) { ji.extend_from_slice(row); ji.push(0.0); }
            ce.extend_from_slice(&original.ce);
            for row in original.je.chunks_exact(n) { je.extend_from_slice(row); je.push(0.0); }
        }
        SqpSample { f: point[n], gradient, ce, ci, je, ji }
    }
}

/// Original SQP and source scenario/case failures, including spent work.
pub type ScenarioStudyError = SqpError<ScenarioError>;
fn normalize(error: ScenarioStudyError) -> ScenarioStudyError {
    match error { SqpError::Evaluation(ScenarioError::Cancelled) => SqpError::Cancelled, other => other }
}

/// Resumable minimax study. A single accepted checkpoint binds all scenarios;
/// the epigraph t is the last optimizer coordinate, never a physical parameter.
/// This searches a local finite-scenario KKT point, not a continuum guarantee.
pub struct ScenarioEquilibriumStudy<'problem, 'work> {
    ensemble: ScenarioProblem<'problem>,
    control: &'work mut DesignControl,
    state: SqpState,
    accepted: ScenarioEvaluation,
    last_domain_rejection: Option<ScenarioError>,
}
impl<'problem, 'work> ScenarioEquilibriumStudy<'problem, 'work> {
    /// Solve the initial complete scenario family once, setting t to its actual
    /// maximum objective. The one initial SQP sample reuses that result; it does
    /// not spend a second family of primal/adjoint solves.
    pub fn new(ensemble: ScenarioProblem<'problem>, nominal: &[f64], control: &'work mut DesignControl,
        gate: &CancelGate) -> Result<Self, ScenarioStudyError>
    {
        let accepted = ensemble.evaluate(nominal, control, gate).map_err(SqpError::Evaluation).map_err(normalize)?;
        let mut initial = nominal.to_vec(); initial.push(accepted.worst_objective);
        let state = SqpState::try_new(&initial, ensemble.maximum_kkt_dimension,
            &mut |_| Ok::<_, ScenarioError>(Some(ensemble.sample(&initial, &accepted))), None)?;
        poll(gate).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(Self { ensemble, control, state, accepted, last_domain_rejection: None })
    }
    /// Accepted SQP state, with the epigraph in its final coordinate.
    #[must_use]
    pub fn optimizer(&self) -> &SqpState { &self.state }
    /// Full physical evidence for exactly the accepted nominal point.
    #[must_use]
    pub fn accepted(&self) -> &ScenarioEvaluation { &self.accepted }
    /// Immutable finite scenario set; cannot be replaced on continuation.
    #[must_use]
    pub fn ensemble(&self) -> &ScenarioProblem<'problem> { &self.ensemble }
    /// Actual physical evaluations and case solves; not ensemble callback count.
    #[must_use]
    pub fn work(&self) -> DesignWork { self.control.work() }
    /// Most recent unavailable realized parameter point, not a physical failure.
    #[must_use]
    pub fn last_domain_rejection(&self) -> Option<&ScenarioError> { self.last_domain_rejection.as_ref() }
    /// Extend physical limits without losing or refunding any spent work.
    pub fn extend_physics_budget(&mut self, evaluations: usize, case_solves: usize) -> Result<(), DesignError> {
        self.control.extend(evaluations, case_solves)
    }
    /// Complete callback count is bounded separately from physical scenario and
    /// load-case work. Only original out-of-domain errors become unavailable
    /// trials. Force, contact-margin, budget and adjoint failures all propagate.
    /// Interrupted searches restart from accepted state with their work charged.
    pub fn run(&mut self, tolerance: f64, additional_iterations: usize, maximum_evaluations: usize,
        gate: &CancelGate) -> Result<SqpRunReport, ScenarioStudyError>
    {
        poll(gate).map_err(SqpError::Evaluation).map_err(normalize)?;
        let mut report = self.advance(tolerance, 0, maximum_evaluations, gate)?;
        for _ in 0..additional_iterations {
            if report.stop != SqpStop::IterationLimit { break; }
            poll(gate).map_err(SqpError::Evaluation).map_err(normalize)?;
            report = self.advance(tolerance, 1, maximum_evaluations, gate)?;
        }
        poll(gate).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(report)
    }
    fn advance(&mut self, tolerance: f64, steps: usize, maximum: usize, gate: &CancelGate)
        -> Result<SqpRunReport, ScenarioStudyError>
    {
        let ensemble = &self.ensemble; let control = &mut *self.control;
        let rejected = &mut self.last_domain_rejection;
        let n = ensemble.lower.len(); let mut candidate = None;
        let result = self.state.try_run(&mut |point| {
            match ensemble.evaluate(&point[..n], control, gate) {
                Ok(value) => {
                    let sample = ensemble.sample(point, &value);
                    candidate = Some((point.to_vec(), value)); Ok(Some(sample))
                }
                Err(error) if error.domain_rejection() => { *rejected = Some(error); Ok(None) }
                Err(error) => Err(error),
            }
        }, tolerance, steps, maximum, None);
        if let Some((point, value)) = candidate {
            if point.as_slice() == self.state.point() { self.accepted = value; }
        }
        let report = result.map_err(normalize)?;
        poll(gate).map_err(SqpError::Evaluation).map_err(normalize)?;
        Ok(report)
    }
    /// Re-solve the accepted nominal point under ALL scenarios. Charged to the
    /// original physical ledger; no cached gradient stands in for this work.
    pub fn recheck(&mut self, gate: &CancelGate) -> Result<ScenarioEvaluation, ScenarioError> {
        let n = self.ensemble.lower.len();
        let evaluation = self.ensemble.evaluate(&self.state.point()[..n], self.control, gate)?;
        if evaluation != self.accepted { return Err(ScenarioError::Invalid("scenario re-solve differs from accepted evidence")); }
        Ok(evaluation)
    }
}

#[cfg(test)]
mod tests;
