//! The solve driver: nonlinear iteration with a DECLARED globalization
//! and a DECLARED stop rule, linear solves through `fs-solver`, and a
//! resumable state that travels inside `fs-exec`'s versioned snapshot
//! envelope.
//!
//! # The stop rule (declared, not tuned in secret)
//!
//! An iterate is accepted when
//!
//! ```text
//!   ‖R(T)‖₂ ≤ residual_rtol · ‖b_f(T)‖₂ + residual_atol
//! ```
//!
//! or when the accepted step satisfies `‖ΔT‖∞ ≤ step_atol` — the second
//! clause exists so a problem whose residual scale is degenerate still
//! terminates on a reason, never on the iteration budget alone.
//! Exhausting `max_iterations` is [`crate::ConductionError::NotConverged`],
//! never a quiet "here is your answer".
//!
//! # The globalization (declared per nonlinearity)
//!
//! - [`Nonlinearity::FixedPoint`]: Picard with a FIXED relaxation `ω`,
//!   `T ← T + ω(y − T)` where `A(T)y = b(T)`. The only backtracking is
//!   the ADMISSIBILITY guard: a trial iterate that leaves the material's
//!   sampled temperature span halves the step rather than extrapolating.
//! - [`Nonlinearity::Newton`]: Armijo backtracking on `‖R‖₂` with the
//!   declared constant `c`, shrink factor, and backtrack budget. The
//!   same admissibility guard applies, and an inadmissible trial counts
//!   as a rejected step.
//!
//! # Which Krylov method, and why
//!
//! The Picard operator `A(T)` is symmetric positive definite (the
//! conduction block is `V gᵀ K g` with SPD `K`; the Robin block is
//! `h ∫λλ` with `h > 0`), so the fixed-point path runs preconditioned CG.
//! The Newton Jacobian carries the `(V/4) gᵀ_a K'(T̄) ∇T_h` term, which is
//! rank-one per element and NOT symmetric, so the Newton path runs
//! FGMRES. `fs-solver` reports CG's residual as a
//! [`ResidualClaim::RecursiveEstimate`]; this crate therefore RECOMPUTES
//! `‖b − Ax‖₂/‖b‖₂` after every solve and gates on that recomputed
//! number ([`LinearSolveEvidence::true_relative_residual`]).

use fs_exec::Cx;
use fs_exec::solver::{SnapshotError, SolverState, codec};
use fs_solver::krylov::{CgState, ResidualClaim, StallDiagnosis};
use fs_solver::nonlinear::FgmresState;
use fs_solver::{CsrOp, norm2};
use fs_sparse::Csr;
use fs_sparse::precond::{Ilu0, Precond, ilu0};

use crate::ConductionError;
use crate::assemble::{
    AssembledSystem, DofMap, assemble_jacobian, assemble_operator, full_residual, reduce,
    reduce_matrix_and_lift, residual,
};
use crate::bc::{ThermalBc, ThermalBoundary};
use crate::field::ScalarField;
use crate::material::{ConductivityModel, ProvenanceClass};
use crate::mesh::ConductionMesh;

/// The pieces a conduction solve needs, borrowed together.
#[derive(Debug, Clone, Copy)]
pub struct ConductionProblem<'m> {
    /// The prepared mesh.
    pub mesh: &'m ConductionMesh,
    /// The boundary partition.
    pub boundary: &'m ThermalBoundary,
    /// The conductivity model.
    pub material: &'m ConductivityModel,
    /// The volumetric source, W/m³.
    pub source: &'m ScalarField,
}

/// Armijo backtracking parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineSearch {
    /// Sufficient-decrease constant `c ∈ (0, 1)`: accept when
    /// `‖R(T + αδ)‖ ≤ (1 − c α)‖R(T)‖`.
    pub armijo_c: f64,
    /// Step shrink factor per backtrack, in `(0, 1)`.
    pub shrink: f64,
    /// Maximum backtracks before the step is refused.
    pub max_backtracks: usize,
}

impl Default for LineSearch {
    fn default() -> Self {
        LineSearch {
            armijo_c: 1e-4,
            shrink: 0.5,
            max_backtracks: 24,
        }
    }
}

/// How the `k(T)` nonlinearity is driven, with its globalization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Nonlinearity {
    /// Damped Picard: `T ← T + ω(y − T)`, `A(T)y = b(T)`, SPD, CG.
    FixedPoint {
        /// Relaxation `ω ∈ (0, 1]`.
        relaxation: f64,
        /// Admissibility backtracks allowed when a trial iterate leaves
        /// the material's sampled span.
        max_backtracks: usize,
    },
    /// Newton with Armijo backtracking, nonsymmetric Jacobian, FGMRES.
    Newton {
        /// The line-search parameters.
        line_search: LineSearch,
    },
}

impl Default for Nonlinearity {
    /// Newton: the only path that converges quadratically on a strongly
    /// temperature-dependent `k(T)`. A linear material makes damped
    /// Picard with `ω = 1` the cheaper equivalent, and callers who know
    /// that say so.
    fn default() -> Self {
        Nonlinearity::Newton {
            line_search: LineSearch::default(),
        }
    }
}

/// The declared stop rule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StopRule {
    /// Relative residual tolerance against `‖b_f‖₂`.
    pub residual_rtol: f64,
    /// Absolute residual floor (W, in the assembled load's units).
    pub residual_atol: f64,
    /// Accepted-step infinity-norm tolerance, K.
    pub step_atol: f64,
    /// Iteration budget; exhausting it is a refusal.
    pub max_iterations: usize,
}

impl Default for StopRule {
    fn default() -> Self {
        StopRule {
            residual_rtol: 1e-10,
            residual_atol: 1e-30,
            step_atol: 1e-12,
            max_iterations: 50,
        }
    }
}

/// Linear-solve budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearConfig {
    /// Relative-residual tolerance handed to the Krylov method AND
    /// gated on this crate's recomputed Euclidean residual.
    pub tolerance: f64,
    /// Krylov iteration (CG) or restart-cycle (FGMRES) budget.
    pub max_iterations: usize,
    /// FGMRES restart length.
    pub restart: usize,
}

impl Default for LinearConfig {
    fn default() -> Self {
        LinearConfig {
            tolerance: 1e-12,
            max_iterations: 20_000,
            restart: 60,
        }
    }
}

/// Where the nonlinear iteration starts.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum InitialGuess {
    /// A uniform temperature, K.
    Uniform(f64),
    /// The arithmetic mean of the prescribed Dirichlet values; falls
    /// back to the midpoint of the material's sampled span, then to
    /// `0 K` when neither exists. The choice is recorded, never hidden.
    #[default]
    DirichletMean,
    /// An explicit free-dof vector.
    Free(Vec<f64>),
}

/// The full solve configuration.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SolveConfig {
    /// Nonlinearity driver and globalization.
    pub nonlinearity: Nonlinearity,
    /// Stop rule.
    pub stop: StopRule,
    /// Linear-solve budget.
    pub linear: LinearConfig,
    /// Starting iterate.
    pub initial: InitialGuess,
}

impl SolveConfig {
    fn validate(&self) -> Result<(), ConductionError> {
        match self.nonlinearity {
            Nonlinearity::FixedPoint { relaxation, .. } => {
                if !(relaxation.is_finite() && relaxation > 0.0 && relaxation <= 1.0) {
                    return Err(ConductionError::Config {
                        parameter: "relaxation",
                        what: format!("omega = {relaxation} must lie in (0, 1]"),
                    });
                }
            }
            Nonlinearity::Newton { line_search } => {
                if !(line_search.armijo_c.is_finite()
                    && line_search.armijo_c > 0.0
                    && line_search.armijo_c < 1.0)
                {
                    return Err(ConductionError::Config {
                        parameter: "armijo_c",
                        what: format!("c = {} must lie in (0, 1)", line_search.armijo_c),
                    });
                }
                if !(line_search.shrink.is_finite()
                    && line_search.shrink > 0.0
                    && line_search.shrink < 1.0)
                {
                    return Err(ConductionError::Config {
                        parameter: "shrink",
                        what: format!("shrink = {} must lie in (0, 1)", line_search.shrink),
                    });
                }
            }
        }
        if !(self.stop.residual_rtol.is_finite() && self.stop.residual_rtol > 0.0) {
            return Err(ConductionError::Config {
                parameter: "residual_rtol",
                what: "must be finite and positive".to_string(),
            });
        }
        if self.stop.max_iterations == 0 {
            return Err(ConductionError::Config {
                parameter: "max_iterations",
                what: "must be at least one".to_string(),
            });
        }
        if !(self.linear.tolerance.is_finite() && self.linear.tolerance > 0.0) {
            return Err(ConductionError::Config {
                parameter: "linear tolerance",
                what: "must be finite and positive".to_string(),
            });
        }
        if self.linear.restart == 0 {
            return Err(ConductionError::Config {
                parameter: "restart",
                what: "must be at least one".to_string(),
            });
        }
        Ok(())
    }
}

/// Why the nonlinear iteration stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StopReason {
    /// `‖R‖₂` met the residual clause of the stop rule.
    ResidualTolerance,
    /// The accepted step met the step clause of the stop rule.
    StepTolerance,
}

impl StopReason {
    /// A stable tag for structured logs.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            StopReason::ResidualTolerance => "residual-tolerance",
            StopReason::StepTolerance => "step-tolerance",
        }
    }
}

/// What one linear solve actually established.
///
/// `reported` is `fs-solver`'s own typed claim — for CG a
/// [`ResidualClaim::RecursiveEstimate`], which drifts in the dangerous
/// direction. `true_relative_residual` is this crate's RECOMPUTED
/// `‖b − Ax‖₂/‖b‖₂`, and `converged_true` gates on that number alone.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearSolveEvidence {
    /// Which nonlinear iteration produced it.
    pub nonlinear_iteration: usize,
    /// `"pcg"` or `"fgmres"`.
    pub method: &'static str,
    /// Krylov iterations performed.
    pub iterations: usize,
    /// The producing solver's typed claim, carried verbatim.
    pub reported: ResidualClaim,
    /// Recomputed `‖b − Ax‖₂/‖b‖₂`.
    pub true_relative_residual: f64,
    /// `true_relative_residual < tolerance`.
    pub converged_true: bool,
    /// The producing solver's stall diagnosis, when it did not converge.
    pub stall: Option<StallDiagnosis>,
}

/// The steady energy balance over the assembled operator.
///
/// The identity being checked: summing the discrete equations against
/// the constant test function `1 = Σ_i φ_i` gives
/// `Σ_i r_i = Q_robin_out + Q_neumann_out − Q_source`, and free rows have
/// `r_i = 0` only up to the ALGEBRAIC residual. So
/// `closure_w = −Σ_{free} r_i` exactly: this number measures how well
/// the linear solves closed the balance, not whether the boundary data
/// itself is right.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyBalance {
    /// `∫_Ω f dV`, W.
    pub source_w: f64,
    /// `∫_{Γ_N} q_n dA`, W (positive = leaving).
    pub neumann_out_w: f64,
    /// `∫_{Γ_R} h (T_h − T_ref) dA`, W (positive = leaving).
    pub robin_out_w: f64,
    /// `Σ_{v ∈ Γ_D} r_v`, W — net heat entering through Dirichlet rows.
    pub dirichlet_in_w: f64,
    /// `source + dirichlet_in − neumann_out − robin_out`, W.
    pub closure_w: f64,
    /// The largest term, used to normalize `closure_w`.
    pub scale_w: f64,
}

impl EnergyBalance {
    /// `|closure_w| / scale_w`.
    #[must_use]
    pub fn relative_closure(&self) -> f64 {
        self.closure_w.abs() / self.scale_w
    }
}

/// The report a solve returns alongside the field.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionReport {
    /// Nonlinear iterations performed.
    pub iterations: usize,
    /// `‖R‖₂` at the start of every iteration, in order.
    pub residual_history: Vec<f64>,
    /// `‖R‖₂` at the accepted solution.
    pub final_residual: f64,
    /// The stop rule's residual threshold at the accepted solution.
    pub residual_threshold: f64,
    /// Why it stopped.
    pub stop_reason: StopReason,
    /// One record per linear solve.
    pub linear: Vec<LinearSolveEvidence>,
    /// The energy-balance closure.
    pub energy: EnergyBalance,
    /// Whether every conductivity number carried an `fs-matdb` receipt.
    pub material_provenance: ProvenanceClass,
    /// How many `fs-matdb` receipts travel with this solve.
    pub material_receipts: usize,
    /// Free degrees of freedom.
    pub free_dofs: usize,
    /// Element count.
    pub elements: usize,
}

/// A solved field with its report.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionSolution {
    /// Nodal temperature over ALL vertices, K (Dirichlet values in
    /// place).
    pub temperature: Vec<f64>,
    /// The report.
    pub report: ConductionReport,
}

/// The resumable nonlinear-iteration state.
///
/// This is the whole thing: `clone()` is a checkpoint and
/// [`SolverState::seal`] is a migratable snapshot inside `fs-exec`'s
/// versioned envelope. A run resumed from a snapshot reproduces the
/// uninterrupted trajectory bitwise, because every linear solve starts
/// from `x₀ = 0` and is a pure function of the assembled system, which
/// is a pure function of this state.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionState {
    /// Iterations completed.
    pub iteration: usize,
    /// The free-dof temperature iterate, K.
    pub free_temperature: Vec<f64>,
    /// `‖R‖₂` observed at the start of every completed iteration.
    pub residual_history: Vec<f64>,
    /// Infinity norm of the last accepted step, K (`f64::INFINITY`
    /// before the first step, so the step clause cannot fire early).
    pub last_step_norm: f64,
}

impl SolverState for ConductionState {
    const TYPE_ID: u64 = 0xf5c0_4d75_c710_0001;
    const SCHEMA_VERSION: u32 = 1;

    fn encode(&self, enc: &mut codec::Enc) {
        enc.put_u64(self.iteration as u64);
        enc.put_f64_slice(&self.free_temperature);
        enc.put_f64_slice(&self.residual_history);
        enc.put_f64(self.last_step_norm);
    }

    fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
        let iteration = dec.get_u64()? as usize;
        let free_temperature = dec.get_f64_vec()?;
        let residual_history = dec.get_f64_vec()?;
        let last_step_norm = dec.get_f64()?;
        Ok(ConductionState {
            iteration,
            free_temperature,
            residual_history,
            last_step_norm,
        })
    }
}

/// Diagonal (Jacobi) preconditioner. Deterministic, allocation-free at
/// apply time, and never a divide by zero: a zero diagonal falls back to
/// `1.0`, which is the identity on that row rather than an infinity.
#[derive(Debug, Clone)]
struct Jacobi {
    inv_diag: Vec<f64>,
}

impl Jacobi {
    fn new(a: &Csr) -> Jacobi {
        let inv_diag = (0..a.nrows())
            .map(|i| {
                let d = a.get(i, i);
                if d.is_finite() && d != 0.0 {
                    1.0 / d
                } else {
                    1.0
                }
            })
            .collect();
        Jacobi { inv_diag }
    }
}

impl Precond for Jacobi {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((zi, ri), di) in z.iter_mut().zip(r).zip(&self.inv_diag) {
            *zi = ri * di;
        }
    }
}

enum ChosenPrecond {
    Incomplete(Box<Ilu0>),
    Diagonal(Jacobi),
}

impl Precond for ChosenPrecond {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        match self {
            ChosenPrecond::Incomplete(p) => p.apply(r, z),
            ChosenPrecond::Diagonal(p) => p.apply(r, z),
        }
    }
}

/// ILU(0) where the factorization exists, Jacobi where it breaks down.
/// Both are deterministic and the choice is a pure function of the
/// matrix, so a replayed solve picks the same one.
pub(crate) fn spd_preconditioner(a: &Csr) -> impl Precond + use<> {
    match ilu0(a) {
        Ok(p) => ChosenPrecond::Incomplete(Box::new(p)),
        Err(_) => ChosenPrecond::Diagonal(Jacobi::new(a)),
    }
}

fn true_relative_residual(a: &Csr, b: &[f64], x: &[f64]) -> f64 {
    let mut ax = vec![0.0f64; b.len()];
    a.spmv(x, &mut ax);
    let r: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
    let bnorm = norm2(b);
    if bnorm > 0.0 {
        norm2(&r) / bnorm
    } else {
        norm2(&r)
    }
}

/// One assembled evaluation of the nonlinear residual at the current
/// iterate: the system, its Dirichlet-reduced form, the free-dof
/// residual, and the stop rule's threshold at this point.
struct Evaluated {
    system: AssembledSystem,
    a_ff: Csr,
    b_f: Vec<f64>,
    r: Vec<f64>,
    rnorm: f64,
    threshold: f64,
    full: Vec<f64>,
}

/// The resumable steady-conduction solver.
pub struct ConductionSolver<'m> {
    problem: ConductionProblem<'m>,
    dofs: DofMap,
    config: SolveConfig,
    state: ConductionState,
    linear: Vec<LinearSolveEvidence>,
    last_system: Option<AssembledSystem>,
    stop_reason: Option<StopReason>,
    final_threshold: f64,
}

/// What one [`ConductionSolver::step`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepReport {
    /// `‖R‖₂` at the START of the step.
    pub residual: f64,
    /// The stop rule's threshold at that point.
    pub threshold: f64,
    /// True when the stop rule fired and NO update was applied.
    pub converged: bool,
    /// Accepted step length `α` (1 when no backtracking happened).
    pub step_length: f64,
    /// Backtracks used.
    pub backtracks: usize,
    /// `‖ΔT‖∞` of the accepted step, K.
    pub step_norm: f64,
}

impl<'m> ConductionSolver<'m> {
    /// Build a solver over a problem.
    ///
    /// # Errors
    /// [`ConductionError::Config`] for an inadmissible configuration;
    /// [`ConductionError::NoFreeDofs`] when everything is pinned;
    /// [`ConductionError::SingularPureNeumann`] when neither a Dirichlet
    /// nor a Robin row exists (the operator is singular up to a
    /// constant, and a Krylov method would return a plausible wrong
    /// answer rather than refuse);
    /// [`ConductionError::FieldLength`] for a mis-sized initial guess.
    pub fn new(
        problem: ConductionProblem<'m>,
        config: SolveConfig,
    ) -> Result<ConductionSolver<'m>, ConductionError> {
        config.validate()?;
        problem
            .source
            .validate("volumetric source", problem.mesh.vertex_count())?;
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count())?;
        if dofs.fixed().is_empty() && !problem.boundary.has_robin() {
            return Err(ConductionError::SingularPureNeumann);
        }
        let free_temperature = initial_free(&problem, &dofs, &config.initial)?;
        Ok(ConductionSolver {
            problem,
            dofs,
            config,
            state: ConductionState {
                iteration: 0,
                free_temperature,
                residual_history: Vec::new(),
                last_step_norm: f64::INFINITY,
            },
            linear: Vec::new(),
            last_system: None,
            stop_reason: None,
            final_threshold: f64::NAN,
        })
    }

    /// The degree-of-freedom map.
    #[must_use]
    pub const fn dofs(&self) -> &DofMap {
        &self.dofs
    }

    /// The current state (a checkpoint by value).
    #[must_use]
    pub const fn state(&self) -> &ConductionState {
        &self.state
    }

    /// The full nodal temperature of the current iterate.
    #[must_use]
    pub fn temperature(&self) -> Vec<f64> {
        self.dofs.scatter(&self.state.free_temperature)
    }

    /// Seal the current state into `fs-exec`'s versioned envelope.
    #[must_use]
    pub fn snapshot(&self, provenance: u64) -> Vec<u8> {
        self.state.seal(provenance)
    }

    /// Restore a sealed state, replacing the current iterate.
    ///
    /// # Errors
    /// [`ConductionError::Snapshot`] for any envelope or payload
    /// refusal, and [`ConductionError::FieldLength`] when the restored
    /// state does not match this problem's free-dof count.
    pub fn restore(&mut self, bytes: &[u8]) -> Result<u64, ConductionError> {
        let (state, provenance) = ConductionState::unseal(bytes).map_err(|e: SnapshotError| {
            ConductionError::Snapshot {
                upstream: e.to_string(),
            }
        })?;
        if state.free_temperature.len() != self.dofs.n() {
            return Err(ConductionError::FieldLength {
                field: "restored free temperature",
                expected: self.dofs.n(),
                found: state.free_temperature.len(),
            });
        }
        self.state = state;
        self.last_system = None;
        self.stop_reason = None;
        Ok(provenance)
    }

    fn merit(&self, cx: &Cx<'_>, free: &[f64]) -> Result<Option<f64>, ConductionError> {
        let full = self.dofs.scatter(free);
        match assemble_operator(
            cx,
            self.problem.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            &full,
        ) {
            Ok(system) => Ok(Some(norm2(&residual(&system, &self.dofs, &full)))),
            Err(ConductionError::OutsideTemperatureSpan { .. }) => Ok(None),
            Err(other) => Err(other),
        }
    }

    fn evaluate(&self, cx: &Cx<'_>) -> Result<Evaluated, ConductionError> {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "nonlinear-iteration",
            at: self.state.iteration,
        })?;
        let full = self.dofs.scatter(&self.state.free_temperature);
        let system = assemble_operator(
            cx,
            self.problem.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            &full,
        )?;
        let (a_ff, b_f) = reduce(&system, &self.dofs);
        let r = residual(&system, &self.dofs, &full);
        let rnorm = norm2(&r);
        let threshold = self
            .config
            .stop
            .residual_rtol
            .mul_add(norm2(&b_f), self.config.stop.residual_atol);
        Ok(Evaluated {
            system,
            a_ff,
            b_f,
            r,
            rnorm,
            threshold,
            full,
        })
    }

    /// One nonlinear iteration. Reports convergence WITHOUT updating
    /// when the stop rule already holds.
    ///
    /// # Errors
    /// Cancellation, material, line-search, and linear-solve refusals.
    pub fn step(&mut self, cx: &Cx<'_>) -> Result<StepReport, ConductionError> {
        let ev = self.evaluate(cx)?;
        let rnorm = ev.rnorm;
        let threshold = ev.threshold;
        self.final_threshold = threshold;

        if rnorm <= threshold {
            self.stop_reason = Some(StopReason::ResidualTolerance);
            self.last_system = Some(ev.system);
            return Ok(StepReport {
                residual: rnorm,
                threshold,
                converged: true,
                step_length: 0.0,
                backtracks: 0,
                step_norm: 0.0,
            });
        }
        if self.state.last_step_norm <= self.config.stop.step_atol {
            self.stop_reason = Some(StopReason::StepTolerance);
            self.last_system = Some(ev.system);
            return Ok(StepReport {
                residual: rnorm,
                threshold,
                converged: true,
                step_length: 0.0,
                backtracks: 0,
                step_norm: self.state.last_step_norm,
            });
        }

        let direction = self.direction(cx, &ev)?;
        let (accepted, alpha, backtracks) = self.globalize(cx, &direction, rnorm)?;

        let step_norm = accepted
            .iter()
            .zip(&self.state.free_temperature)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        self.state.free_temperature = accepted;
        self.state.residual_history.push(rnorm);
        self.state.last_step_norm = step_norm;
        self.state.iteration += 1;
        self.last_system = None;
        Ok(StepReport {
            residual: rnorm,
            threshold,
            converged: false,
            step_length: alpha,
            backtracks,
            step_norm,
        })
    }

    /// The search direction for one iteration: the damped Picard update
    /// or the Newton step, each through its own Krylov method.
    fn direction(&mut self, cx: &Cx<'_>, ev: &Evaluated) -> Result<Vec<f64>, ConductionError> {
        match self.config.nonlinearity {
            Nonlinearity::FixedPoint { relaxation, .. } => {
                let y = self.solve_spd(&ev.a_ff, &ev.b_f)?;
                Ok(y.iter()
                    .zip(&self.state.free_temperature)
                    .map(|(yi, xi)| relaxation * (yi - xi))
                    .collect())
            }
            Nonlinearity::Newton { .. } => {
                let jacobian = assemble_jacobian(
                    cx,
                    self.problem.mesh,
                    self.problem.boundary,
                    self.problem.material,
                    &ev.full,
                )?;
                let (j_ff, _) = reduce_matrix_and_lift(&jacobian, &self.dofs);
                let rhs: Vec<f64> = ev.r.iter().map(|v| -v).collect();
                self.solve_general(&j_ff, &rhs)
            }
        }
    }

    /// Apply the DECLARED globalization: Armijo sufficient decrease for
    /// Newton (`c > 0`), plain acceptance for damped Picard (`c = 0`).
    /// Both share the ADMISSIBILITY guard — a trial iterate outside the
    /// material's sampled span has no merit value at all, so it is
    /// rejected and the step is shrunk rather than extrapolated.
    ///
    /// Returns `(accepted iterate, step length, backtracks used)`.
    fn globalize(
        &mut self,
        cx: &Cx<'_>,
        direction: &[f64],
        base_residual: f64,
    ) -> Result<(Vec<f64>, f64, usize), ConductionError> {
        let (armijo_c, shrink, budget) = match self.config.nonlinearity {
            Nonlinearity::FixedPoint { max_backtracks, .. } => (0.0, 0.5, max_backtracks),
            Nonlinearity::Newton { line_search } => (
                line_search.armijo_c,
                line_search.shrink,
                line_search.max_backtracks,
            ),
        };
        let mut alpha = 1.0f64;
        let mut backtracks = 0usize;
        loop {
            let candidate: Vec<f64> = self
                .state
                .free_temperature
                .iter()
                .zip(direction)
                .map(|(x, d)| alpha.mul_add(*d, *x))
                .collect();
            let admissible = self.merit(cx, &candidate)?;
            if let Some(trial) = admissible
                && trial <= armijo_c.mul_add(-alpha, 1.0) * base_residual
            {
                return Ok((candidate, alpha, backtracks));
            }
            if backtracks >= budget {
                return Err(ConductionError::LineSearchFailed {
                    iteration: self.state.iteration,
                    backtracks,
                    smallest_step: alpha,
                });
            }
            backtracks += 1;
            alpha *= shrink;
        }
    }

    fn solve_spd(&mut self, a: &Csr, b: &[f64]) -> Result<Vec<f64>, ConductionError> {
        let precond = spd_preconditioner(a);
        let op = CsrOp::symmetric(a.clone());
        let mut cg = CgState::new(&op, &precond, b);
        let report = cg.run(
            &op,
            &precond,
            self.config.linear.tolerance,
            self.config.linear.max_iterations,
        );
        let truth = true_relative_residual(a, b, &cg.x);
        let converged_true = truth < self.config.linear.tolerance;
        self.linear.push(LinearSolveEvidence {
            nonlinear_iteration: self.state.iteration,
            method: "pcg",
            iterations: report.iters,
            reported: report.residual_claim(),
            true_relative_residual: truth,
            converged_true,
            stall: report.diagnosis,
        });
        if converged_true {
            Ok(cg.x)
        } else {
            Err(ConductionError::LinearSolveFailed {
                iteration: self.state.iteration,
                krylov_iterations: report.iters,
                true_relative_residual: truth,
                tolerance: self.config.linear.tolerance,
            })
        }
    }

    fn solve_general(&mut self, a: &Csr, b: &[f64]) -> Result<Vec<f64>, ConductionError> {
        let precond = spd_preconditioner(a);
        let op = CsrOp::general(a.clone());
        let mut gm = FgmresState::new(b, self.config.linear.restart);
        let cycles = self
            .config
            .linear
            .max_iterations
            .div_ceil(self.config.linear.restart.max(1));
        let report = gm.run(
            &op,
            &precond,
            b,
            self.config.linear.tolerance,
            cycles.max(1),
        );
        let truth = true_relative_residual(a, b, &gm.x);
        let converged_true = truth < self.config.linear.tolerance;
        self.linear.push(LinearSolveEvidence {
            nonlinear_iteration: self.state.iteration,
            method: "fgmres",
            iterations: report.iters,
            reported: report.residual_claim(),
            true_relative_residual: truth,
            converged_true,
            stall: report.diagnosis,
        });
        if converged_true {
            Ok(gm.x)
        } else {
            Err(ConductionError::LinearSolveFailed {
                iteration: self.state.iteration,
                krylov_iterations: report.iters,
                true_relative_residual: truth,
                tolerance: self.config.linear.tolerance,
            })
        }
    }

    /// Iterate until the stop rule fires.
    ///
    /// # Errors
    /// [`ConductionError::NotConverged`] when the iteration budget runs
    /// out; every refusal [`ConductionSolver::step`] can produce.
    pub fn run(&mut self, cx: &Cx<'_>) -> Result<ConductionSolution, ConductionError> {
        loop {
            let report = self.step(cx)?;
            if report.converged {
                break;
            }
            if self.state.iteration >= self.config.stop.max_iterations {
                // Budget exhausted: evaluate the final iterate WITHOUT
                // updating it, so the refusal reports the residual of
                // the field the caller would have received.
                let ev = self.evaluate(cx)?;
                self.final_threshold = ev.threshold;
                if ev.rnorm <= ev.threshold {
                    self.stop_reason = Some(StopReason::ResidualTolerance);
                    self.last_system = Some(ev.system);
                    break;
                }
                return Err(ConductionError::NotConverged {
                    iterations: self.state.iteration,
                    residual: ev.rnorm,
                    threshold: ev.threshold,
                });
            }
        }
        Ok(self.finish())
    }

    fn finish(&mut self) -> ConductionSolution {
        let system = self
            .last_system
            .as_ref()
            .expect("a converged step caches its assembled system");
        let temperature = self.dofs.scatter(&self.state.free_temperature);
        let energy = energy_balance(
            self.problem.mesh,
            self.problem.boundary,
            self.problem.source,
            system,
            &self.dofs,
            &temperature,
        );
        let final_residual = norm2(&residual(system, &self.dofs, &temperature));
        ConductionSolution {
            temperature,
            report: ConductionReport {
                iterations: self.state.iteration,
                residual_history: self.state.residual_history.clone(),
                final_residual,
                residual_threshold: self.final_threshold,
                stop_reason: self
                    .stop_reason
                    .expect("a converged run records its stop reason"),
                linear: self.linear.clone(),
                energy,
                material_provenance: self.problem.material.provenance(),
                material_receipts: self.problem.material.receipts().len(),
                free_dofs: self.dofs.n(),
                elements: self.problem.mesh.element_count(),
            },
        }
    }
}

fn initial_free(
    problem: &ConductionProblem<'_>,
    dofs: &DofMap,
    guess: &InitialGuess,
) -> Result<Vec<f64>, ConductionError> {
    let n = dofs.n();
    match guess {
        InitialGuess::Uniform(t) => {
            crate::require_finite("initial temperature", *t)?;
            Ok(vec![*t; n])
        }
        InitialGuess::Free(values) => {
            if values.len() != n {
                return Err(ConductionError::FieldLength {
                    field: "initial free temperature",
                    expected: n,
                    found: values.len(),
                });
            }
            for &v in values {
                crate::require_finite("initial temperature", v)?;
            }
            Ok(values.clone())
        }
        InitialGuess::DirichletMean => {
            let fixed = dofs.fixed();
            if !fixed.is_empty() {
                let sum: f64 = fixed.iter().map(|&v| dofs.prescribed()[v]).sum();
                return Ok(vec![sum / fixed.len() as f64; n]);
            }
            let t0 = match problem.material.temperature_span() {
                crate::material::TemperatureSpan::Sampled { low, high } => f64::midpoint(low, high),
                crate::material::TemperatureSpan::Unbounded => 0.0,
            };
            Ok(vec![t0; n])
        }
    }
}

fn energy_balance(
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    source: &ScalarField,
    system: &AssembledSystem,
    dofs: &DofMap,
    temperature: &[f64],
) -> EnergyBalance {
    let mut source_w = 0.0f64;
    for (e, tet) in mesh.complex().tets.iter().enumerate() {
        let volume = mesh.element_volume(e);
        let mut acc = 0.0f64;
        for &v in tet {
            acc += source.at(v as usize);
        }
        source_w += volume * acc / 4.0;
    }

    let mut neumann_out_w = 0.0f64;
    let mut robin_out_w = 0.0f64;
    for (slot, face) in mesh.boundary().iter().enumerate() {
        let Some(condition) = boundary.condition_for(slot) else {
            continue;
        };
        let verts = [
            face.vertices[0] as usize,
            face.vertices[1] as usize,
            face.vertices[2] as usize,
        ];
        match condition {
            ThermalBc::Dirichlet { .. } => {}
            ThermalBc::Neumann { outward_flux } => {
                let mean = verts.iter().map(|&v| outward_flux.at(v)).sum::<f64>() / 3.0;
                neumann_out_w = face.area.mul_add(mean, neumann_out_w);
            }
            ThermalBc::Robin { htc, t_ref } => {
                let h_bar = verts.iter().map(|&v| htc.at(v)).sum::<f64>() / 3.0;
                let delta = verts
                    .iter()
                    .map(|&v| temperature[v] - t_ref.at(v))
                    .sum::<f64>()
                    / 3.0;
                robin_out_w = (h_bar * face.area).mul_add(delta, robin_out_w);
            }
        }
    }

    let r = full_residual(system, temperature);
    let dirichlet_in_w: f64 = dofs.fixed().iter().map(|&v| r[v]).sum();
    let closure_w = source_w + dirichlet_in_w - neumann_out_w - robin_out_w;
    let scale_w = source_w
        .abs()
        .max(neumann_out_w.abs())
        .max(robin_out_w.abs())
        .max(dirichlet_in_w.abs())
        .max(f64::MIN_POSITIVE);
    EnergyBalance {
        source_w,
        neumann_out_w,
        robin_out_w,
        dirichlet_in_w,
        closure_w,
        scale_w,
    }
}

/// Recover the per-element heat-flux vector `q_e = −K(T̄_e) ∇T_h`, W/m².
///
/// P₁ gradients are element-wise constant, so this is the natural
/// (piecewise-constant, one-order-lower) flux. It is NOT a
/// flux-conserving recovery and this crate makes no superconvergence
/// claim about it.
///
/// # Errors
/// [`ConductionError::OutsideTemperatureSpan`] when an element mean
/// leaves the material's sampled span.
pub fn element_heat_flux(
    mesh: &ConductionMesh,
    material: &ConductivityModel,
    temperature: &[f64],
) -> Result<Vec<[f64; 3]>, ConductionError> {
    let mut out = Vec::with_capacity(mesh.element_count());
    for e in 0..mesh.element_count() {
        let tet = mesh.complex().tets[e];
        let t_e = crate::assemble::element_temperature(mesh, e, temperature);
        let k = material.tensor_at(t_e)?;
        let g = &mesh.geometry().grads[e];
        let mut grad = [0.0f64; 3];
        for (b, gb) in g.iter().enumerate() {
            let tb = temperature[tet[b] as usize];
            for i in 0..3 {
                grad[i] = tb.mul_add(gb[i], grad[i]);
            }
        }
        let mut q = [0.0f64; 3];
        for (i, qi) in q.iter_mut().enumerate() {
            *qi = -k[i][0].mul_add(grad[0], k[i][1].mul_add(grad[1], k[i][2] * grad[2]));
        }
        out.push(q);
    }
    Ok(out)
}

/// Solve a steady conduction problem to the configured stop rule.
///
/// # Errors
/// Every refusal [`ConductionSolver::new`] and
/// [`ConductionSolver::run`] can produce.
pub fn solve(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    config: SolveConfig,
) -> Result<ConductionSolution, ConductionError> {
    ConductionSolver::new(problem, config)?.run(cx)
}
