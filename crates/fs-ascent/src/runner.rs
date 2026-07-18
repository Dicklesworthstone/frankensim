//! The Problem-IR STUDY RUNNER (bead ijil): drive an `fs_opt::Problem`
//! end-to-end through the fs-ascent engines — the seam where "problems
//! are data" meets "optimizers are engines".
//!
//! - VARIABLE PACKING: all declared variables concatenate into one
//!   point vector across the manifold product, with a distinct packed
//!   optimizer-parameter vector. Per-block operations delegate to fs-opt's
//!   authoritative manifold runtime so orientations stay orientations.
//! - GRADIENTS: central differences through `fs_opt::eval` in the
//!   retraction-parameter representation (the live IR is evaluation-only — no
//!   reverse mode; documented, deterministic fixed h).
//! - CONSTRAINTS: `EqZero`/`LeZero` roots route to the constrained
//!   engines (AL always; IP/SQP by option) through raw point-coordinate
//!   adapters for Euclidean packings; non-Euclidean packings fail closed.
//! - BUDGET: `Problem.budget.limit` threads into the stop algebra
//!   (`StopRule::Budget`) alongside the caller's rules.
//! - RESUMABLE: [`Study`] is the checkpoint (clone = checkpoint); a
//!   split run is BITWISE equal to a straight run (house pattern,
//!   gated).
//! - CANCELLABLE: `try_run_cancellable` observes `fs_exec::Cx` only at
//!   complete iteration boundaries and returns a replay-complete typed pause.
//! - BOUND + FORKABLE: a checkpoint retains its admitted semantic problem
//!   identity, so cached objective state cannot cross problem meanings;
//!   reweight steering is an explicit immutable-parent world fork.

use crate::stop::{StopObservation, StopReason, StopRule};
use fs_exec::Cx;
use fs_opt::{
    AdmissionReport, ConstraintKind, EvalLimit, Manifold, OptError, Problem, ProblemSemanticId,
    Sense, Variable, eval,
};

/// Version of the resume-safe Study cancellation boundary and pause receipt.
///
/// Bump when polling placement or receipt semantics change; retained G4
/// evidence must never silently cross a different pause protocol.
pub const STUDY_CANCELLATION_BOUNDARY_VERSION: u16 = 1;

/// One packed variable block.
#[derive(Debug, Clone)]
struct Block {
    manifold: Manifold,
    /// Offset in the packed point vector.
    point_start: usize,
    /// Point storage length.
    point_len: usize,
    /// Offset in the packed optimizer-parameter vector.
    parameter_start: usize,
    /// Optimizer-parameter storage length.
    parameter_len: usize,
}

/// The packed view of a problem's variables.
#[derive(Debug, Clone)]
pub struct Packing {
    blocks: Vec<Block>,
    /// Total packed point-storage length.
    pub dim: usize,
    /// Total packed optimizer-parameter length.
    pub param_dim: usize,
}

impl Packing {
    /// Total packed point-storage length.
    ///
    /// This is the explicit spelling of the retained public `dim` field.
    #[must_use]
    pub const fn point_dim(&self) -> usize {
        self.dim
    }

    /// Build the packing for a problem's variable list.
    #[must_use]
    pub fn new(problem: &Problem) -> Packing {
        let mut blocks = Vec::new();
        let mut point_start = 0usize;
        let mut parameter_start = 0usize;
        for v in problem.vars() {
            let layout = v
                .manifold
                .layout()
                .expect("sealed problems carry validated manifold layouts");
            let point_len = usize::try_from(layout.point_dim().get())
                .expect("point storage dimension fits usize");
            let parameter_len = usize::try_from(layout.param_dim().get())
                .expect("parameter storage dimension fits usize");
            let point_end = point_start
                .checked_add(point_len)
                .expect("sealed problem admission bounds total packed storage for this target");
            let parameter_end = parameter_start.checked_add(parameter_len).expect(
                "sealed problem admission bounds total packed parameter storage for this target",
            );
            blocks.push(Block {
                manifold: v.manifold,
                point_start,
                point_len,
                parameter_start,
                parameter_len,
            });
            point_start = point_end;
            parameter_start = parameter_end;
        }
        Packing {
            blocks,
            dim: point_start,
            param_dim: parameter_start,
        }
    }

    /// Split a packed point into per-variable bindings for `eval`.
    #[must_use]
    pub fn unpack(&self, x: &[f64]) -> Vec<Vec<f64>> {
        assert_eq!(
            x.len(),
            self.dim,
            "packed point length must match Packing::dim"
        );
        self.blocks
            .iter()
            .map(|block| x[block.point_start..block.point_start + block.point_len].to_vec())
            .collect()
    }

    /// Pull a packed ambient point-storage gradient back to the packed
    /// optimizer-parameter representation block by block.
    #[must_use]
    pub fn project(&self, x: &[f64], g: &[f64]) -> Vec<f64> {
        assert_eq!(
            x.len(),
            self.dim,
            "packed point length must match Packing::dim"
        );
        assert_eq!(
            g.len(),
            self.dim,
            "packed ambient gradient length must match Packing::dim"
        );
        let mut out = vec![0.0f64; self.param_dim];
        for block in &self.blocks {
            let point = &x[block.point_start..block.point_start + block.point_len];
            let ambient = &g[block.point_start..block.point_start + block.point_len];
            let parameter = block
                .manifold
                .parameter_gradient(point, ambient)
                .unwrap_or_else(|error| {
                    panic!("packed gradient failed fs-opt manifold authority: {error}")
                });
            assert_eq!(parameter.len(), block.parameter_len);
            out[block.parameter_start..block.parameter_start + block.parameter_len]
                .copy_from_slice(&parameter);
        }
        out
    }

    /// Retract a packed optimizer-parameter step into packed point storage.
    #[must_use]
    pub fn retract(&self, x: &[f64], step: &[f64]) -> Vec<f64> {
        assert_eq!(
            x.len(),
            self.dim,
            "packed point length must match Packing::dim"
        );
        assert_eq!(
            step.len(),
            self.param_dim,
            "packed step length must match Packing::param_dim"
        );
        let mut out = vec![0.0f64; self.dim];
        for block in &self.blocks {
            let point = &x[block.point_start..block.point_start + block.point_len];
            let parameter =
                &step[block.parameter_start..block.parameter_start + block.parameter_len];
            let landed = block
                .manifold
                .retract(point, parameter)
                .unwrap_or_else(|error| {
                    panic!("packed retraction failed fs-opt manifold authority: {error}")
                });
            assert_eq!(landed.len(), block.point_len);
            out[block.point_start..block.point_start + block.point_len].copy_from_slice(&landed);
        }
        out
    }

    fn authoritative_point(&self, point: &[f64]) -> Result<Vec<f64>, (usize, OptError)> {
        assert_eq!(
            point.len(),
            self.dim,
            "packed point length must match Packing::dim"
        );
        let mut authoritative = point.to_vec();
        for (variable, block) in self.blocks.iter().enumerate() {
            let stored = &point[block.point_start..block.point_start + block.point_len];
            let zero_parameter = vec![0.0; block.parameter_len];
            let validated = block
                .manifold
                .retract(stored, &zero_parameter)
                .map_err(|source| (variable, source))?;
            assert_eq!(validated.len(), block.point_len);
            if matches!(block.manifold, Manifold::So3) {
                authoritative[block.point_start..block.point_start + block.point_len]
                    .copy_from_slice(&validated);
            }
        }
        Ok(authoritative)
    }

    fn matches_problem_schema(&self, problem: &Problem) -> bool {
        self.blocks.len() == problem.vars().len()
            && self
                .blocks
                .iter()
                .zip(problem.vars())
                .all(|(block, variable)| block.manifold == variable.manifold)
    }
}

/// A resumable study over an IR problem: state clone = checkpoint.
#[derive(Debug, Clone)]
pub struct Study {
    packing: Packing,
    problem_id: ProblemSemanticId,
    variables: Vec<Variable>,
    /// Current packed iterate.
    pub x: Vec<f64>,
    /// Objective history (most recent last).
    pub history: Vec<f64>,
    /// Evaluations spent (budget accounting).
    pub evals: usize,
    /// Steps taken.
    pub steps: usize,
    current_f: Option<f64>,
    current_grad_norm: Option<f64>,
    fd_h: f64,
    lr: f64,
}

/// Typed refusal from constructing, continuing, or steering a [`Study`].
#[derive(Debug, Clone, PartialEq)]
pub enum StudyError {
    /// The supplied problem failed the common `fs-opt` admission gate.
    ProblemRejected(AdmissionReport),
    /// The packed point cannot bind the problem's declared variables.
    PackedPointLength {
        /// Required packed storage length.
        expected: usize,
        /// Supplied packed storage length.
        actual: usize,
    },
    /// A packed variable point violates its authoritative manifold contract.
    ManifoldPointInvalid {
        /// Variable-list index whose point was refused.
        variable: usize,
        /// Typed fs-opt manifold refusal.
        source: OptError,
    },
    /// An objective root cannot be evaluated at the branch point.
    ObjectiveUnevaluable {
        /// Objective-list index.
        index: usize,
        /// Evaluation refusal.
        source: OptError,
    },
    /// A constraint root cannot be evaluated at the branch point.
    ConstraintUnevaluable {
        /// Constraint-list index.
        index: usize,
        /// Evaluation refusal.
        source: OptError,
    },
    /// A checkpoint was presented to a semantically different problem.
    ProblemMismatch {
        /// Problem identity stored in the checkpoint.
        bound: ProblemSemanticId,
        /// Problem identity supplied for continuation.
        supplied: ProblemSemanticId,
    },
    /// A fork was requested without changing the semantic problem.
    ForkProblemUnchanged {
        /// Identity shared by the parent and proposed child.
        problem_id: ProblemSemanticId,
    },
    /// Steering attempted to change the variable schema, not only the study objective.
    ForkVariableSchemaMismatch {
        /// First unequal variable index, or the first missing index when counts differ.
        first_mismatch: usize,
        /// Number of variables in the parent study schema.
        parent_variables: usize,
        /// Number of variables in the proposed child problem.
        child_variables: usize,
    },
}

impl core::fmt::Display for StudyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            StudyError::ProblemRejected(report) => write!(f, "{report}"),
            StudyError::PackedPointLength { expected, actual } => write!(
                f,
                "packed study point has length {actual}, but the problem requires {expected}"
            ),
            StudyError::ManifoldPointInvalid { variable, source } => write!(
                f,
                "packed study variable {variable} is not an authoritative manifold point: {source}"
            ),
            StudyError::ObjectiveUnevaluable { index, source } => {
                write!(
                    f,
                    "objective {index} is unevaluable at the study point: {source}"
                )
            }
            StudyError::ConstraintUnevaluable { index, source } => {
                write!(
                    f,
                    "constraint {index} is unevaluable at the study point: {source}"
                )
            }
            StudyError::ProblemMismatch { bound, supplied } => write!(
                f,
                "study checkpoint is bound to problem {bound}, not supplied problem {supplied}; use Study::fork_for to steer"
            ),
            StudyError::ForkProblemUnchanged { problem_id } => write!(
                f,
                "study fork keeps semantic problem {problem_id} unchanged; resume the checkpoint instead of resetting branch-local accounting"
            ),
            StudyError::ForkVariableSchemaMismatch {
                first_mismatch,
                parent_variables,
                child_variables,
            } => write!(
                f,
                "study fork changes the variable schema at index {first_mismatch} (parent variables {parent_variables}, child variables {child_variables})"
            ),
        }
    }
}

impl std::error::Error for StudyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            StudyError::ProblemRejected(report) => Some(report),
            StudyError::ManifoldPointInvalid { source, .. }
            | StudyError::ObjectiveUnevaluable { source, .. }
            | StudyError::ConstraintUnevaluable { source, .. } => Some(source),
            StudyError::PackedPointLength { .. }
            | StudyError::ProblemMismatch { .. }
            | StudyError::ForkProblemUnchanged { .. }
            | StudyError::ForkVariableSchemaMismatch { .. } => None,
        }
    }
}

/// Replay-complete description of a world-fork steering operation.
///
/// The child starts at `branch_point_bits` with the same numerical driver
/// configuration but fresh branch-local history, step, evaluation, and cache
/// state. The parent is borrowed immutably and remains untouched.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StudyForkReceipt {
    /// Semantic problem identity bound to the parent checkpoint.
    pub parent_problem_id: ProblemSemanticId,
    /// Semantic problem identity bound to the child branch.
    pub child_problem_id: ProblemSemanticId,
    /// Parent steps already landed when the branch was created.
    pub parent_steps: usize,
    /// Parent evaluations already spent when the branch was created.
    pub parent_evals: usize,
    /// Parent objective-history length at the branch point.
    pub parent_history_len: usize,
    /// Exact packed branch point as IEEE-754 bit patterns.
    pub branch_point_bits: Vec<u64>,
    /// Finite-difference step as an IEEE-754 bit pattern.
    pub fd_h_bits: u64,
    /// Learning rate as an IEEE-754 bit pattern.
    pub learning_rate_bits: u64,
}

/// Replay-complete state observed at a resume-safe cancellation boundary.
///
/// The receipt contains every mutable numerical field that is not already
/// fixed by `problem_id`. Cancellation observation itself does not mutate the
/// [`Study`]; callers may retain this receipt, replace the cancelled context,
/// and resume the same in-memory checkpoint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StudyPauseReceipt {
    /// Boundary/receipt semantics version.
    pub boundary_version: u16,
    /// Semantic problem identity bound to the checkpoint.
    pub problem_id: ProblemSemanticId,
    /// Completed retraction-parameter-gradient steps.
    pub steps: usize,
    /// Objective evaluations spent so far.
    pub evals: usize,
    /// Exact packed point bits.
    pub point_bits: Vec<u64>,
    /// Exact objective-history bits.
    pub history_bits: Vec<u64>,
    /// Cached current objective, when populated.
    pub current_objective_bits: Option<u64>,
    /// Cached retraction-parameter-gradient norm, when populated.
    pub current_gradient_norm_bits: Option<u64>,
    /// Finite-difference step bits.
    pub fd_h_bits: u64,
    /// Learning-rate bits.
    pub learning_rate_bits: u64,
}

/// Typed outcome of a cancellation-aware study segment.
#[derive(Debug, Clone)]
#[must_use]
pub enum StudyRunProgress {
    /// The caller's stop rule, problem budget, or segment cap stopped the run.
    Stopped(StudyReport),
    /// The context requested cancellation at a resume-safe boundary.
    Paused(StudyPauseReceipt),
}

/// Outcome of a study segment.
#[derive(Debug, Clone)]
pub struct StudyReport {
    /// Final objective.
    pub f: f64,
    /// Final retraction-parameter gradient infinity norm.
    pub grad_norm: f64,
    /// Why the segment stopped.
    pub reason: StopReason,
    /// Evaluations spent so far (cumulative).
    pub evals: usize,
}

/// The weighted, sense-corrected total objective of a problem at
/// packed point `x` (all objectives, not just the first).
fn objective(problem: &Problem, packing: &Packing, x: &[f64], evals: &mut usize) -> f64 {
    let bindings = packing.unpack(x);
    let mut total = 0.0f64;
    for o in problem.objectives() {
        let v = eval(problem, o.node, &bindings)
            .expect("objective evaluable (checked at study start)")
            .scalar()
            .expect("objective roots are scalar");
        let sign = match o.sense {
            Sense::Minimize => 1.0,
            Sense::Maximize => -1.0,
        };
        total += sign * o.weight * v;
    }
    *evals += 1;
    total
}

impl Study {
    /// Start a study at `x0` (packed). Verifies evaluability of every
    /// objective and constraint root up front (fail loud, fail early).
    ///
    /// # Panics
    /// If the problem is refused by admission, `x0` has the wrong packed
    /// length or violates a declared manifold, or a root is unevaluable. Use
    /// [`Study::try_new`] for a typed refusal.
    #[must_use]
    pub fn new(problem: &Problem, x0: &[f64], fd_h: f64, lr: f64) -> Study {
        Study::try_new(problem, x0, fd_h, lr)
            .unwrap_or_else(|error| panic!("study construction refused: {error}"))
    }

    /// Start a study with typed admission, manifold, binding, and evaluation
    /// refusals. SO(3) starts are stored in the authority's deterministic
    /// antipodal representative before any objective or constraint is read.
    #[must_use]
    pub fn try_new(problem: &Problem, x0: &[f64], fd_h: f64, lr: f64) -> Result<Study, StudyError> {
        let admission = problem.admit().map_err(StudyError::ProblemRejected)?;
        Study::try_new_admitted(problem, admission.semantic_id(), x0, fd_h, lr)
    }

    fn try_new_admitted(
        problem: &Problem,
        problem_id: ProblemSemanticId,
        x0: &[f64],
        fd_h: f64,
        lr: f64,
    ) -> Result<Study, StudyError> {
        let packing = Packing::new(problem);
        if x0.len() != packing.dim {
            return Err(StudyError::PackedPointLength {
                expected: packing.dim,
                actual: x0.len(),
            });
        }
        let point = packing
            .authoritative_point(x0)
            .map_err(|(variable, source)| StudyError::ManifoldPointInvalid { variable, source })?;
        let bindings = packing.unpack(&point);
        for (index, objective) in problem.objectives().iter().enumerate() {
            let _ = eval(problem, objective.node, &bindings)
                .map_err(|source| StudyError::ObjectiveUnevaluable { index, source })?;
        }
        for (index, constraint) in problem.constraints().iter().enumerate() {
            let _ = eval(problem, constraint.node, &bindings)
                .map_err(|source| StudyError::ConstraintUnevaluable { index, source })?;
        }
        Ok(Study {
            packing,
            problem_id,
            variables: problem.vars().to_vec(),
            x: point,
            history: Vec::new(),
            evals: 0,
            steps: 0,
            current_f: None,
            current_grad_norm: None,
            fd_h,
            lr,
        })
    }

    /// Semantic identity of the problem this checkpoint may continue.
    #[must_use]
    pub const fn problem_id(&self) -> ProblemSemanticId {
        self.problem_id
    }

    /// Fork this checkpoint into a fresh branch for a reweighted problem.
    ///
    /// Steering may change objectives, constraints, tags, and budgets, but it
    /// may not reinterpret the checkpoint's packed coordinates: variable names,
    /// manifolds, dimensions, order, and count must remain exactly equal. The
    /// child inherits only the current point and numerical driver configuration;
    /// problem-dependent caches and branch-local counters start empty.
    #[must_use]
    pub fn fork_for(&self, problem: &Problem) -> Result<(Study, StudyForkReceipt), StudyError> {
        let admission = problem.admit().map_err(StudyError::ProblemRejected)?;
        let child_problem_id = admission.semantic_id();
        if child_problem_id == self.problem_id {
            return Err(StudyError::ForkProblemUnchanged {
                problem_id: child_problem_id,
            });
        }
        if self.variables.as_slice() != problem.vars() {
            let first_mismatch = self
                .variables
                .iter()
                .zip(problem.vars())
                .position(|(parent, child)| parent != child)
                .unwrap_or_else(|| self.variables.len().min(problem.vars().len()));
            return Err(StudyError::ForkVariableSchemaMismatch {
                first_mismatch,
                parent_variables: self.variables.len(),
                child_variables: problem.vars().len(),
            });
        }
        let child =
            Study::try_new_admitted(problem, child_problem_id, &self.x, self.fd_h, self.lr)?;
        let receipt = StudyForkReceipt {
            parent_problem_id: self.problem_id,
            child_problem_id,
            parent_steps: self.steps,
            parent_evals: self.evals,
            parent_history_len: self.history.len(),
            branch_point_bits: self.x.iter().map(|value| value.to_bits()).collect(),
            fd_h_bits: self.fd_h.to_bits(),
            learning_rate_bits: self.lr.to_bits(),
        };
        Ok((child, receipt))
    }

    fn pause_receipt(&self) -> StudyPauseReceipt {
        StudyPauseReceipt {
            boundary_version: STUDY_CANCELLATION_BOUNDARY_VERSION,
            problem_id: self.problem_id,
            steps: self.steps,
            evals: self.evals,
            point_bits: self.x.iter().map(|value| value.to_bits()).collect(),
            history_bits: self.history.iter().map(|value| value.to_bits()).collect(),
            current_objective_bits: self.current_f.map(f64::to_bits),
            current_gradient_norm_bits: self.current_grad_norm.map(f64::to_bits),
            fd_h_bits: self.fd_h.to_bits(),
            learning_rate_bits: self.lr.to_bits(),
        }
    }

    /// Central-difference retraction-parameter gradient of the total objective.
    fn retraction_parameter_gradient(&mut self, problem: &Problem) -> Vec<f64> {
        let n = self.packing.param_dim;
        let mut g = vec![0.0f64; n];
        for i in 0..n {
            let mut t = vec![0.0f64; n];
            t[i] = self.fd_h;
            let xp = self.packing.retract(&self.x, &t);
            t[i] = -self.fd_h;
            let xm = self.packing.retract(&self.x, &t);
            let fp = objective(problem, &self.packing, &xp, &mut self.evals);
            let fm = objective(problem, &self.packing, &xm, &mut self.evals);
            g[i] = (fp - fm) / (2.0 * self.fd_h);
        }
        g
    }

    /// Run a segment, failing loud if the supplied problem does not match the
    /// semantic identity bound to this checkpoint.
    ///
    /// # Panics
    /// If the problem is refused by admission or differs from the checkpoint's
    /// bound problem. Use [`Study::try_run`] for a typed refusal.
    pub fn run(&mut self, problem: &Problem, rule: &StopRule, max_steps: usize) -> StudyReport {
        self.try_run(problem, rule, max_steps)
            .unwrap_or_else(|error| panic!("study continuation refused: {error}"))
    }

    /// Run a segment with typed problem-admission and semantic-binding refusal.
    ///
    /// Retraction-parameter-gradient steps retract until the stop rule fires, the bound
    /// problem budget runs out, or `max_steps` is reached. The problem's own
    /// `budget.limit` is always added to the caller's rules (P4).
    pub fn try_run(
        &mut self,
        problem: &Problem,
        rule: &StopRule,
        max_steps: usize,
    ) -> Result<StudyReport, StudyError> {
        self.ensure_problem_matches(problem)?;
        match self.run_bound_progress(problem, rule, max_steps, None) {
            StudyRunProgress::Stopped(report) => Ok(report),
            StudyRunProgress::Paused(_) => {
                unreachable!("a run without a cancellation context cannot pause")
            }
        }
    }

    /// Run a segment under an explicit execution context.
    ///
    /// Cancellation is observed only at resume-safe boundaries: before any
    /// segment work, after the initial objective cache is complete, before an
    /// iteration, and after a complete retract/evaluate iteration. A request
    /// arriving inside the finite-difference gradient therefore drains the
    /// current iteration before [`StudyRunProgress::Paused`] is returned. The
    /// optimizer's [`StopReason`] remains reserved for optimization semantics;
    /// cancellation is a distinct typed outcome.
    ///
    /// # Errors
    ///
    /// Returns [`StudyError`] if the problem fails admission or its semantic
    /// identity differs from the checkpoint's bound problem.
    pub fn try_run_cancellable(
        &mut self,
        problem: &Problem,
        rule: &StopRule,
        max_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<StudyRunProgress, StudyError> {
        self.ensure_problem_matches(problem)?;
        Ok(self.run_bound_progress(problem, rule, max_steps, Some(cx)))
    }

    fn ensure_problem_matches(&self, problem: &Problem) -> Result<(), StudyError> {
        let supplied = problem
            .admit()
            .map_err(StudyError::ProblemRejected)?
            .semantic_id();
        if supplied != self.problem_id {
            return Err(StudyError::ProblemMismatch {
                bound: self.problem_id,
                supplied,
            });
        }
        Ok(())
    }

    fn run_bound_progress(
        &mut self,
        problem: &Problem,
        rule: &StopRule,
        max_steps: usize,
        cx: Option<&Cx<'_>>,
    ) -> StudyRunProgress {
        if cx.is_some_and(Cx::is_cancel_requested) {
            return StudyRunProgress::Paused(self.pause_receipt());
        }
        let mut rules = vec![rule.clone()];
        if let EvalLimit::Limited(maximum) = problem.budget().limit {
            rules.push(StopRule::Budget(
                usize::try_from(maximum.get()).unwrap_or(usize::MAX),
            ));
        }
        let combined = StopRule::Any(rules);
        let mut reason = StopReason::IterationCap;
        let mut f = if let Some(f) = self.current_f {
            f
        } else {
            let f = objective(problem, &self.packing, &self.x, &mut self.evals);
            self.current_f = Some(f);
            f
        };
        let mut gnorm = self.current_grad_norm.unwrap_or(f64::INFINITY);
        if cx.is_some_and(Cx::is_cancel_requested) {
            return StudyRunProgress::Paused(self.pause_receipt());
        }
        let obs = StopObservation {
            grad_norm: gnorm,
            objective: f,
            evals: self.evals,
            history: &self.history,
        };
        if let Some(r) = combined.check(&obs) {
            return StudyRunProgress::Stopped(StudyReport {
                f,
                grad_norm: gnorm,
                reason: r,
                evals: self.evals,
            });
        }
        for _ in 0..max_steps {
            if cx.is_some_and(Cx::is_cancel_requested) {
                return StudyRunProgress::Paused(self.pause_receipt());
            }
            let g = self.retraction_parameter_gradient(problem);
            gnorm = g.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
            self.current_grad_norm = Some(gnorm);
            self.history.push(f);
            let obs = StopObservation {
                grad_norm: gnorm,
                objective: f,
                evals: self.evals,
                history: &self.history,
            };
            if let Some(r) = combined.check(&obs) {
                reason = r;
                break;
            }
            let step: Vec<f64> = g.iter().map(|v| -self.lr * v).collect();
            self.x = self.packing.retract(&self.x, &step);
            self.steps += 1;
            f = objective(problem, &self.packing, &self.x, &mut self.evals);
            self.current_f = Some(f);
            self.current_grad_norm = None;
            if cx.is_some_and(Cx::is_cancel_requested) {
                return StudyRunProgress::Paused(self.pause_receipt());
            }
        }
        StudyRunProgress::Stopped(StudyReport {
            f,
            grad_norm: gnorm,
            reason,
            evals: self.evals,
        })
    }

    /// Constrained adapters: expose the problem's `EqZero`/`LeZero`
    /// roots as packed callbacks for [`crate::augmented_lagrangian`],
    /// [`crate::interior_point`] or [`crate::sqp`] — with FD
    /// Jacobian-transpose actions in raw packed point coordinates. These
    /// adapters deliberately do not reinterpret constrained-manifold
    /// derivatives as optimizer parameters; constrained manifold solving
    /// remains outside this compatibility path's claim. Returns (ce, ce_jt,
    /// ci, ci_jt) closures over the problem.
    ///
    /// # Panics
    /// If the problem or `packing` contains a non-Euclidean block, or if the
    /// packing's variable schema differs from the problem's.
    #[allow(clippy::type_complexity)]
    pub fn constraint_adapters<'p>(
        problem: &'p Problem,
        packing: &'p Packing,
        fd_h: f64,
    ) -> (
        impl Fn(&[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64], &[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64]) -> Vec<f64> + 'p,
        impl Fn(&[f64], &[f64]) -> Vec<f64> + 'p,
    ) {
        assert!(
            problem
                .vars()
                .iter()
                .all(|variable| matches!(variable.manifold, Manifold::Rn { .. })),
            "raw point-coordinate constraint adapters support Rn problem variables only"
        );
        assert!(
            packing.matches_problem_schema(problem),
            "constraint adapter packing must match the problem variable schema"
        );
        assert!(
            packing
                .blocks
                .iter()
                .all(|block| matches!(block.manifold, Manifold::Rn { .. })),
            "raw point-coordinate constraint adapters support Rn blocks only"
        );
        let eval_kind = move |x: &[f64], kind: ConstraintKind| -> Vec<f64> {
            let bindings = packing.unpack(x);
            problem
                .constraints()
                .iter()
                .filter(|c| c.kind == kind)
                .map(|c| {
                    eval(problem, c.node, &bindings)
                        .expect("constraint evaluable")
                        .scalar()
                        .expect("constraint roots are scalar")
                })
                .collect()
        };
        let jt_kind = move |x: &[f64], w: &[f64], kind: ConstraintKind| -> Vec<f64> {
            // FD of wᵀc(x) — one directional pass per packed dim.
            let n = x.len();
            let mut out = vec![0.0f64; n];
            for i in 0..n {
                let mut xp = x.to_vec();
                xp[i] += fd_h;
                let mut xm = x.to_vec();
                xm[i] -= fd_h;
                let cp = eval_kind(&xp, kind);
                let cm = eval_kind(&xm, kind);
                out[i] = cp
                    .iter()
                    .zip(&cm)
                    .zip(w)
                    .map(|((p, m2), wi)| wi * (p - m2) / (2.0 * fd_h))
                    .sum();
            }
            out
        };
        (
            move |x: &[f64]| eval_kind(x, ConstraintKind::EqZero),
            move |x: &[f64], w: &[f64]| jt_kind(x, w, ConstraintKind::EqZero),
            move |x: &[f64]| eval_kind(x, ConstraintKind::LeZero),
            move |x: &[f64], w: &[f64]| jt_kind(x, w, ConstraintKind::LeZero),
        )
    }
}
