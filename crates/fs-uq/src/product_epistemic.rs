//! Worst-case propagation of epistemic parameter boxes, not probability sampling.
//!
//! The callback must enclose ONE fixed real-valued model on every supplied closed
//! box, including singleton boxes. Ordinary point samples or corner evaluations
//! are not enclosures. Bounds below are conditional on that callback contract;
//! this producer does not mint Verified evidence or validate physical models.
//!
//! Bead: `frankensim-extreal-program-f85xj.6.7`; plan §§9.8 and 9.9.
use core::fmt::Display;
use crate::{CorrelationModel, PropagationMethod, UncertaintyKind, UqPlan};

/// Explicit accuracy and retained subdivision limits. Tolerance has QoI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EpistemicConfig {
    pub absolute_tolerance: f64,
    /// At most 4096 live leaves; two prospective children may coexist with them.
    pub maximum_cells: usize,
}

/// Completion means BOTH extrema are enclosed to the requested absolute width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpistemicStatus {
    Paused,
    Complete,
    /// A requested upper-limit decision is established, not extrema convergence.
    /// Ordinary `advance` may resume the same paid work to refine the extrema.
    DecisionReached,
    EvaluationLimit,
    CellLimit,
    ParameterResolution,
    Cancelled,
    Refused,
}
impl EpistemicStatus {
    fn terminal(self) -> bool {
        !matches!(self, Self::Paused | Self::Cancelled | Self::DecisionReached)
    }
}

/// Upper-limit decisions conditional on the model's whole-box enclosures.
/// These refer to the declared Cartesian domain, not a probability or an
/// unrepresented jointly feasible set. Equality satisfies the upper limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpistemicCompliance {
    NotRequested,
    Undetermined,
    AllCompliant,
    AllViolating,
    WitnessedViolation,
}
impl EpistemicCompliance {
    fn decided(self) -> bool {
        matches!(self, Self::AllCompliant | Self::AllViolating | Self::WitnessedViolation)
    }
}

/// An actual parameter vector and its enclosed response, not a midpoint estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct EpistemicWitness {
    /// Original parameter order and physical units, within the Cartesian box.
    /// Additional joint physical feasibility constraints are not represented.
    pub parameters: Vec<f64>,
    pub value: [f64; 2],
}

/// Conditional bounds on the model over the full declared Cartesian product.
/// None means root coverage is not yet established, or the evaluator refused.
/// No probability, quantile, confidence interval, or evidence promotion is made.
#[derive(Debug, Clone, PartialEq)]
pub struct EpistemicReport {
    pub qoi_name: String,
    pub status: EpistemicStatus,
    pub evaluations_planned: usize,
    pub evaluations_attempted: usize,
    pub evaluations_accepted: usize,
    /// All covered leaves, including retired interior domains.
    pub covered_cells: usize,
    pub live_cells: usize,
    /// Cells whose entire ranges lie between accepted extremum witnesses.
    pub retired_cells: usize,
    pub pending_split: bool,
    pub output_range: Option<[f64; 2]>,
    /// Encloses the infimum; no continuity or attainment is assumed.
    pub minimum: Option<[f64; 2]>,
    /// Encloses the supremum; this is NOT an interval of observed samples.
    pub maximum: Option<[f64; 2]>,
    pub minimizing_witness: Option<EpistemicWitness>,
    pub maximizing_witness: Option<EpistemicWitness>,
    pub compliance: EpistemicCompliance,
    pub rejection_reason: Option<String>,
}

#[derive(Debug, Clone)]
struct Cell {
    domain: Vec<[f64; 2]>,
    depth: Vec<u32>,
    range: [f64; 2],
}
#[derive(Debug, Clone)]
struct Draft {
    domain: Vec<[f64; 2]>,
    depth: Vec<u32>,
    range: Option<[f64; 2]>,
}
#[derive(Debug, Clone)]
struct Batch {
    // The complete parent remains in `cells` until BOTH children are complete.
    parent: Option<usize>,
    children: Vec<Draft>,
    next: usize,
}

/// Bounded branch-and-bound execution for 1..=64 finite epistemic intervals.
///
/// The domain is explicitly their full Cartesian product. Unknown dependence
/// does not become probabilistic independence; any additional feasible-set
/// restriction is outside this API. Correlation matrices and probability laws
/// refuse. The seed is retained in the plan but no random draws are performed.
///
/// The plan's lifetime budget counts box AND singleton enclosure calls, at most
/// 65536. A completed child call survives cancellation/clone resumption. An
/// unfinished callback retries the identical box. A failed callback is terminal:
/// it cannot be silently skipped. Reports retain prior complete coverage during
/// interruption; no numerical bound is published after a terminal refusal.
///
/// Selection scans at most 4096 cells and 64 axes, polling between cells. Report
/// reduction is bounded synchronous work. Long callbacks own their cancellation.
/// Storage is bounded by live cells plus two drafts and two witnesses; allocator
/// overhead and infallible Clone allocation are not a caller-owned byte budget.
#[derive(Debug, Clone)]
pub struct EpistemicExecution {
    plan: UqPlan,
    config: EpistemicConfig,
    cells: Vec<Cell>,
    retired_range: Option<[f64; 2]>,
    retired: usize,
    pending: Option<Batch>,
    smallest: Option<EpistemicWitness>,
    largest: Option<EpistemicWitness>,
    attempted: usize,
    accepted: usize,
    status: EpistemicStatus,
    failure: Option<String>,
}

impl EpistemicExecution {
    pub fn new(plan: &UqPlan, config: EpistemicConfig) -> Result<Self, &'static str> {
        if plan.method != PropagationMethod::EpistemicBounding {
            return Err("epistemic execution requires an EpistemicBounding plan");
        }
        if plan.target_qoi.trim().is_empty() || !(1..=64).contains(&plan.parameters.len()) {
            return Err("epistemic execution needs a named QoI and 1..=64 parameters");
        }
        if !(1..=65_536).contains(&plan.budget_max_samples)
            || !(1..=4096).contains(&config.maximum_cells)
            || !config.absolute_tolerance.is_finite() || config.absolute_tolerance <= 0.0
            || plan.compliance_threshold.is_some_and(|x| !x.is_finite())
        {
            return Err("invalid epistemic tolerance, cell limit, evaluation budget or threshold");
        }
        if !matches!(plan.correlation, CorrelationModel::Unknown | CorrelationModel::Independent) {
            return Err("a correlation matrix does not define an epistemic feasible set");
        }
        let mut domain = Vec::with_capacity(plan.parameters.len());
        for (i, parameter) in plan.parameters.iter().enumerate() {
            if parameter.name.trim().is_empty() || parameter.unit.trim().is_empty()
                || plan.parameters[..i].iter().any(|p| p.name == parameter.name)
            {
                return Err("epistemic parameters need unique names and explicit units");
            }
            let UncertaintyKind::EpistemicInterval { lo, hi } = parameter.kind else {
                return Err("epistemic bounding requires intervals, not sampled probability laws");
            };
            validate([lo, hi])?;
            domain.push([lo, hi]);
        }
        let mut cells = Vec::new();
        cells.try_reserve_exact(config.maximum_cells)
            .map_err(|_| "epistemic frontier allocation refused")?;
        let root = Draft { depth: vec![0; domain.len()], domain, range: None };
        Ok(Self {
            plan: plan.clone(), config, cells, retired_range: None, retired: 0,
            pending: Some(Batch { parent: None, children: vec![root], next: 0 }),
            smallest: None, largest: None, attempted: 0, accepted: 0,
            status: EpistemicStatus::Paused, failure: None,
        })
    }

    #[must_use]
    pub fn plan(&self) -> &UqPlan { &self.plan }

    #[must_use]
    pub fn report(&self) -> EpistemicReport {
        let extrema = if self.failure.is_none() { self.extrema().ok().flatten() } else { None };
        EpistemicReport {
            qoi_name: self.plan.target_qoi.clone(), status: self.status,
            evaluations_planned: self.plan.budget_max_samples,
            evaluations_attempted: self.attempted, evaluations_accepted: self.accepted,
            covered_cells: (self.cells.len() + self.retired).max(usize::from(extrema.is_some())),
            live_cells: self.cells.len(), retired_cells: self.retired,
            pending_split: self.pending.as_ref().is_some_and(|p| p.parent.is_some()),
            output_range: extrema.map(|(minimum, maximum)| [minimum[0], maximum[1]]),
            minimum: extrema.map(|e| e.0), maximum: extrema.map(|e| e.1),
            minimizing_witness: if self.failure.is_none() { self.smallest.clone() } else { None },
            maximizing_witness: if self.failure.is_none() { self.largest.clone() } else { None },
            compliance: self.compliance(),
            rejection_reason: self.failure.clone(),
        }
    }

    pub fn advance<F, E, C>(&mut self, allowance: usize, cancelled: C, mut enclose: F)
        -> EpistemicReport
    where F: FnMut(&[[f64; 2]]) -> Result<[f64; 2], E>, E: Display, C: FnMut() -> bool {
        self.advance_interruptible(allowance, cancelled, |domain| enclose(domain).map(Some))
    }

    /// Every callback result encloses the same model on the ENTIRE closed box.
    /// `Ok(None)` pauses without consuming an evaluation; already paid results
    /// remain. Errors, nonfinite/reversed bounds and detected inconsistency
    /// refuse permanently. Panics propagate; side effects must tolerate retries.
    pub fn advance_interruptible<F, E, C>(&mut self, allowance: usize,
        cancelled: C, enclose: F) -> EpistemicReport
    where F: FnMut(&[[f64; 2]]) -> Result<Option<[f64; 2]>, E>,
        E: Display, C: FnMut() -> bool {
        self.advance_for(allowance, cancelled, enclose, false)
    }

    /// Stop as soon as the upper compliance limit is established for the whole
    /// domain or contradicted by a point witness. This may need only the root
    /// enclosure, without solving either extremum to the requested tolerance.
    /// Missing threshold returns an error without mutating the execution.
    pub fn advance_until_compliance<F, E, C>(&mut self, allowance: usize,
        cancelled: C, mut enclose: F) -> Result<EpistemicReport, &'static str>
    where F: FnMut(&[[f64; 2]]) -> Result<[f64; 2], E>, E: Display, C: FnMut() -> bool {
        self.advance_until_compliance_interruptible(allowance, cancelled,
            |domain| enclose(domain).map(Some))
    }

    /// Interruptible decision-driven variant; all paid calls remain resumable.
    pub fn advance_until_compliance_interruptible<F, E, C>(&mut self, allowance: usize,
        cancelled: C, enclose: F) -> Result<EpistemicReport, &'static str>
    where F: FnMut(&[[f64; 2]]) -> Result<Option<[f64; 2]>, E>,
        E: Display, C: FnMut() -> bool {
        if self.plan.compliance_threshold.is_none() {
            return Err("decision-driven epistemic execution needs a compliance threshold");
        }
        Ok(self.advance_for(allowance, cancelled, enclose, true))
    }

    fn advance_for<F, E, C>(&mut self, allowance: usize,
        mut cancelled: C, mut enclose: F, stop_on_decision: bool) -> EpistemicReport
    where F: FnMut(&[[f64; 2]]) -> Result<Option<[f64; 2]>, E>,
        E: Display, C: FnMut() -> bool {
        if allowance == 0 || self.status.terminal() { return self.report(); }
        let mut remaining = allowance.min(self.plan.budget_max_samples - self.attempted);
        self.status = EpistemicStatus::Paused;
        loop {
            if cancelled() { self.status = EpistemicStatus::Cancelled; break; }
            if stop_on_decision && self.compliance().decided() {
                self.status = EpistemicStatus::DecisionReached; break;
            }
            if self.pending.is_none() {
                let extrema = match self.extrema() {
                    Ok(Some(value)) => value,
                    _ => return self.fail("epistemic coverage invariant failed".into()),
                };
                if width(extrema.0) <= self.config.absolute_tolerance
                    && width(extrema.1) <= self.config.absolute_tolerance {
                    self.status = EpistemicStatus::Complete; break;
                }
            }
            if self.attempted == self.plan.budget_max_samples {
                self.status = EpistemicStatus::EvaluationLimit; break;
            }
            if remaining == 0 { break; }
            if self.pending.is_none() {
                match self.start_split(&mut cancelled) {
                    Ok(()) => {},
                    Err(EpistemicStatus::Refused) =>
                        return self.fail("epistemic coverage invariant failed".into()),
                    Err(status) => { self.status = status; break; }
                }
            }
            let batch = self.pending.as_ref().expect("root or admitted split");
            let draft = &batch.children[batch.next];
            let request: Vec<_> = if draft.range.is_some() {
                draft.domain.iter().map(|p| { let m = p[0].midpoint(p[1]); [m, m] }).collect()
            } else { draft.domain.clone() };
            let result = enclose(&request);
            if matches!(&result, Ok(None)) { self.status = EpistemicStatus::Cancelled; break; }
            self.attempted += 1;
            let value = match result {
                Ok(Some(value)) => value,
                Err(error) => return self.fail(format!("epistemic model refused: {error}")),
                Ok(None) => unreachable!("unfinished evaluation returned before accounting"),
            };
            if let Err(reason) = self.accept(value) { return self.fail(reason.into()); }
            self.accepted += 1;
            remaining -= 1;
        }
        self.report()
    }

    fn fail(&mut self, reason: String) -> EpistemicReport {
        self.status = EpistemicStatus::Refused;
        self.failure = Some(reason);
        self.report()
    }

    fn extrema(&self) -> Result<Option<([f64; 2], [f64; 2])>, &'static str> {
        let mut hull = self.retired_range;
        for cell in &self.cells { include(&mut hull, cell.range); }
        if self.cells.is_empty() {
            if let Some(batch) = &self.pending {
                if let Some(range) = batch.children[0].range { include(&mut hull, range); }
            }
        }
        let Some(range) = hull else { return Ok(None); };
        let minimum = [range[0], self.smallest.as_ref().map_or(range[1], |w| w.value[1].min(range[1]))];
        let maximum = [self.largest.as_ref().map_or(range[0], |w| w.value[0].max(range[0])), range[1]];
        validate(minimum).and_then(|_| validate(maximum))
            .map_err(|_| "model enclosures contradict an accepted point witness")?;
        Ok(Some((minimum, maximum)))
    }

    fn start_split(&mut self, cancelled: &mut impl FnMut() -> bool)
        -> Result<(), EpistemicStatus> {
        self.retire_interior(cancelled)?;
        let (minimum, maximum) = self.extrema().ok().flatten().ok_or(EpistemicStatus::Refused)?;
        let mut chosen: Option<(usize, usize, f64)> = None;
        for (index, cell) in self.cells.iter().enumerate() {
            if cancelled() { return Err(EpistemicStatus::Cancelled); }
            // Priority arithmetic is heuristic only; no bound depends on it.
            let score = (minimum[1] - cell.range[0]).max(cell.range[1] - maximum[0]);
            if score <= 0.0 { continue; }
            let axis = cell.domain.iter().enumerate()
                .filter(|(_, p)| { let m = p[0].midpoint(p[1]); p[0] < m && m < p[1] })
                .min_by_key(|(axis, _)| (cell.depth[*axis], *axis)).map(|(axis, _)| axis);
            if let Some(axis) = axis {
                if chosen.is_none_or(|(_, _, old)| score > old) { chosen = Some((index, axis, score)); }
            }
        }
        let Some((parent, axis, _)) = chosen else { return Err(EpistemicStatus::ParameterResolution); };
        if self.cells.len() == self.config.maximum_cells { return Err(EpistemicStatus::CellLimit); }
        let cell = &self.cells[parent];
        let midpoint = cell.domain[axis][0].midpoint(cell.domain[axis][1]);
        let mut left = Draft { domain: cell.domain.clone(), depth: cell.depth.clone(), range: None };
        let mut right = Draft { domain: cell.domain.clone(), depth: cell.depth.clone(), range: None };
        left.domain[axis][1] = midpoint;
        right.domain[axis][0] = midpoint;
        left.depth[axis] += 1;
        right.depth[axis] += 1;
        self.pending = Some(Batch { parent: Some(parent), children: vec![left, right], next: 0 });
        Ok(())
    }

    fn compliance(&self) -> EpistemicCompliance {
        let Some(threshold) = self.plan.compliance_threshold else {
            return EpistemicCompliance::NotRequested;
        };
        if self.failure.is_some() { return EpistemicCompliance::Undetermined; }
        let Ok(Some((minimum, maximum))) = self.extrema() else {
            return EpistemicCompliance::Undetermined;
        };
        if maximum[1] <= threshold { EpistemicCompliance::AllCompliant }
        else if minimum[0] > threshold { EpistemicCompliance::AllViolating }
        else if self.largest.as_ref().is_some_and(|w| w.value[0] > threshold) {
            EpistemicCompliance::WitnessedViolation
        } else { EpistemicCompliance::Undetermined }
    }

    fn retire_interior(&mut self, cancelled: &mut impl FnMut() -> bool)
        -> Result<(), EpistemicStatus> {
        let (Some(smallest), Some(largest)) = (&self.smallest, &self.largest) else {
            return Ok(());
        };
        let floor = smallest.value[1];
        let ceiling = largest.value[0];
        let mut index = 0;
        while index < self.cells.len() {
            if cancelled() { return Err(EpistemicStatus::Cancelled); }
            let range = self.cells[index].range;
            if range[0] >= floor && range[1] <= ceiling {
                // Preserve an aggregate inclusion of every retired subdomain.
                // Its endpoints cannot bind either extremum beyond the current
                // witnesses; subsequent better witnesses only strengthen that.
                include(&mut self.retired_range, range);
                self.retired += 1;
                self.cells.swap_remove(index);
            } else { index += 1; }
        }
        Ok(())
    }

    fn accept(&mut self, value: [f64; 2]) -> Result<(), &'static str> {
        validate(value)?;
        let batch = self.pending.as_mut().expect("one admitted pending evaluation");
        let draft = &mut batch.children[batch.next];
        let point_range = if let Some(range) = draft.range {
            intersect(value, range)?
        } else {
            let range = match batch.parent {
                Some(parent) => intersect(value, self.cells[parent].range)?,
                None => value,
            };
            draft.range = Some(range);
            if draft.domain.iter().any(|p| p[0] != p[1]) { return Ok(()); }
            range // A singleton box is already its own point evaluation.
        };
        let witness = EpistemicWitness {
            parameters: draft.domain.iter().map(|p| p[0].midpoint(p[1])).collect(),
            value: point_range,
        };
        if self.smallest.as_ref().is_none_or(|old| point_range[1] < old.value[1]) {
            self.smallest = Some(witness.clone());
        }
        if self.largest.as_ref().is_none_or(|old| point_range[0] > old.value[0]) {
            self.largest = Some(witness);
        }
        batch.next += 1;
        if batch.next == batch.children.len() {
            let batch = self.pending.take().expect("completed batch");
            for (i, draft) in batch.children.into_iter().enumerate() {
                let cell = Cell { domain: draft.domain, depth: draft.depth,
                    range: draft.range.expect("box evaluated before its point") };
                if i == 0 {
                    if let Some(parent) = batch.parent { self.cells[parent] = cell; }
                    else { self.cells.push(cell); }
                } else { self.cells.push(cell); }
            }
        }
        let _ = self.extrema()?;
        Ok(())
    }
}

fn validate(range: [f64; 2]) -> Result<(), &'static str> {
    if range.iter().any(|x| !x.is_finite()) || range[0] > range[1] {
        Err("epistemic enclosure must have finite ordered endpoints")
    } else { Ok(()) }
}
fn intersect(a: [f64; 2], b: [f64; 2]) -> Result<[f64; 2], &'static str> {
    let range = [a[0].max(b[0]), a[1].min(b[1])];
    validate(range).map_err(|_| "model returned disjoint nested enclosures")?;
    Ok(range)
}
fn include(hull: &mut Option<[f64; 2]>, range: [f64; 2]) {
    *hull = Some(hull.map_or(range, |h| [h[0].min(range[0]), h[1].max(range[1])]));
}
fn width(range: [f64; 2]) -> f64 {
    if range[0] == range[1] { 0.0 } else { (range[1] - range[0]).next_up() }
}
