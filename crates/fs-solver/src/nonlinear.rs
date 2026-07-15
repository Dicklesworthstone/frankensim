//! Resumable Newton--Krylov, flexible GMRES, and solver admission.
//!
//! The nonlinear driver owns globalization mechanics but no physics. Problem
//! crates provide residual and Jacobian actions; callers provide independent
//! linear-system verification before selecting a named Krylov method.

use crate::{LinearOp, SolveReport, StallDiagnosis, dot, norm2};
use fs_sparse::precond::Precond;
use std::fmt;

/// Symmetry conclusion from an injected linear-system verifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetryEvidence {
    /// Symmetry was established within the verifier's declared scope.
    Symmetric,
    /// The operator is nonsymmetric.
    Nonsymmetric,
    /// No symmetry conclusion is available.
    Unknown,
}

/// Definiteness conclusion from an injected verifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DefinitenessEvidence {
    /// Symmetric positive definite on the admitted working space.
    PositiveDefinite,
    /// Symmetric indefinite on the admitted working space.
    Indefinite,
    /// No definiteness conclusion is available.
    Unknown,
}

/// Nullspace handling recorded by the verifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullspaceEvidence {
    /// The admitted working operator has trivial nullspace.
    Trivial,
    /// A known nullspace of this dimension was projected from operator and
    /// source.
    Projected {
        /// Projected nullity.
        dimension: usize,
    },
    /// A possible nullspace remains unresolved.
    Unresolved,
}

/// Compatibility of the right-hand side with constraints/nullspaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceCompatibility {
    /// Source compatibility was checked by the injected verifier.
    Compatible,
    /// The source violates a known compatibility condition.
    Incompatible,
    /// No source-compatibility conclusion is available.
    Unknown,
}

/// Preconditioner behavior relevant to Krylov admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreconditionerClass {
    /// No preconditioner.
    None,
    /// One fixed SPD preconditioner.
    FixedSpd,
    /// One fixed general preconditioner.
    FixedGeneral,
    /// The preconditioner may vary by logical iteration.
    Variable,
}

/// Verifier-produced structural finding. The wrapper minted by
/// [`verify_linear_system`] retains verifier identity and validates internally
/// contradictory combinations before admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LinearSystemFinding {
    /// Symmetry conclusion.
    pub symmetry: SymmetryEvidence,
    /// Definiteness conclusion.
    pub definiteness: DefinitenessEvidence,
    /// Nullspace treatment.
    pub nullspace: NullspaceEvidence,
    /// Source compatibility.
    pub source: SourceCompatibility,
    /// Preconditioner behavior.
    pub preconditioner: PreconditionerClass,
}

/// Injected verification capability for one exact operator/source pair.
///
/// This crate validates the returned grammar and applies admission rules. It
/// does not authenticate implementations of this trait; external authority is
/// a caller/ledger responsibility.
pub trait LinearSystemVerifier {
    /// Stable machine-readable verifier identity.
    fn verifier_id(&self) -> &str;
    /// Verify an exact operator/source pair.
    fn verify(
        &self,
        operator: &dyn LinearOp,
        rhs: &[f64],
    ) -> Result<LinearSystemFinding, LinearVerificationError>;
}

/// Refusal while producing the local verification receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LinearVerificationError {
    /// Operator and source dimensions differ.
    Dimension {
        /// Operator dimension.
        operator: usize,
        /// Source length.
        rhs: usize,
    },
    /// A source entry is NaN or infinite.
    NonFiniteSource {
        /// Entry index.
        index: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// Verifier identity is empty.
    EmptyVerifierId,
    /// The verifier refused the pair.
    Rejected(&'static str),
    /// A positive-definite finding omitted required symmetry/trivial-nullspace
    /// evidence.
    ContradictoryPositiveDefiniteFinding,
    /// An indefinite finding omitted required symmetry evidence.
    ContradictoryIndefiniteFinding,
    /// A projected nullspace had zero dimension.
    EmptyProjectedNullspace,
    /// A projected nullspace exceeded the operator dimension.
    ProjectedNullspaceTooLarge {
        /// Reported nullity.
        dimension: usize,
        /// Operator dimension.
        operator: usize,
    },
}

impl fmt::Display for LinearVerificationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dimension { operator, rhs } => {
                write!(
                    f,
                    "linear operator dimension {operator} differs from rhs {rhs}"
                )
            }
            Self::NonFiniteSource { index, bits } => write!(
                f,
                "linear source entry {index} is non-finite (bits 0x{bits:016x})"
            ),
            Self::EmptyVerifierId => f.write_str("linear-system verifier identity is empty"),
            Self::Rejected(reason) => write!(f, "linear-system verifier refused: {reason}"),
            Self::ContradictoryPositiveDefiniteFinding => f.write_str(
                "positive-definite finding requires symmetric evidence and a trivial nullspace",
            ),
            Self::ContradictoryIndefiniteFinding => {
                f.write_str("indefinite finding requires symmetric evidence")
            }
            Self::EmptyProjectedNullspace => {
                f.write_str("projected nullspace dimension must be positive")
            }
            Self::ProjectedNullspaceTooLarge {
                dimension,
                operator,
            } => write!(
                f,
                "projected nullspace dimension {dimension} exceeds operator dimension {operator}"
            ),
        }
    }
}

impl core::error::Error for LinearVerificationError {}

/// Local verification receipt consumed by solver admission.
///
/// This value retains the dimension, verifier identity, and finding, but not a
/// content identity for an opaque operator or source. It is therefore an
/// admission decision input rather than an unforgeable execution capability.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedLinearSystem {
    dimension: usize,
    verifier_id: String,
    finding: LinearSystemFinding,
}

impl VerifiedLinearSystem {
    /// Verified operator/source dimension.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Injected verifier identity.
    #[must_use]
    pub fn verifier_id(&self) -> &str {
        &self.verifier_id
    }

    /// Retained structural finding.
    #[must_use]
    pub const fn finding(&self) -> LinearSystemFinding {
        self.finding
    }
}

/// Verify shape/finiteness, invoke a verifier, and retain its coherent finding.
#[must_use]
pub fn verify_linear_system<V: LinearSystemVerifier>(
    operator: &dyn LinearOp,
    rhs: &[f64],
    verifier: &V,
) -> Result<VerifiedLinearSystem, LinearVerificationError> {
    if operator.n() != rhs.len() {
        return Err(LinearVerificationError::Dimension {
            operator: operator.n(),
            rhs: rhs.len(),
        });
    }
    for (index, value) in rhs.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(LinearVerificationError::NonFiniteSource {
                index,
                bits: value.to_bits(),
            });
        }
    }
    let verifier_id = verifier.verifier_id().to_owned();
    if verifier_id.is_empty() {
        return Err(LinearVerificationError::EmptyVerifierId);
    }
    let finding = verifier.verify(operator, rhs)?;
    if finding.definiteness == DefinitenessEvidence::PositiveDefinite
        && (finding.symmetry != SymmetryEvidence::Symmetric
            || finding.nullspace != NullspaceEvidence::Trivial)
    {
        return Err(LinearVerificationError::ContradictoryPositiveDefiniteFinding);
    }
    if finding.definiteness == DefinitenessEvidence::Indefinite
        && finding.symmetry != SymmetryEvidence::Symmetric
    {
        return Err(LinearVerificationError::ContradictoryIndefiniteFinding);
    }
    if matches!(
        finding.nullspace,
        NullspaceEvidence::Projected { dimension: 0 }
    ) {
        return Err(LinearVerificationError::EmptyProjectedNullspace);
    }
    if let NullspaceEvidence::Projected { dimension } = finding.nullspace
        && dimension > operator.n()
    {
        return Err(LinearVerificationError::ProjectedNullspaceTooLarge {
            dimension,
            operator: operator.n(),
        });
    }
    Ok(VerifiedLinearSystem {
        dimension: operator.n(),
        verifier_id,
        finding,
    })
}

/// Named Krylov algorithm subject to structural admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearSolverKind {
    /// Conjugate gradients.
    Cg,
    /// MINRES.
    Minres,
    /// Restarted GMRES with a fixed preconditioner.
    Gmres,
    /// Flexible GMRES with a logical-iteration-keyed preconditioner.
    Fgmres,
}

/// Typed refusal from solver selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverAdmissionError {
    /// Source compatibility was absent or negative.
    SourceNotCompatible,
    /// A possible nullspace remains unresolved.
    NullspaceUnresolved,
    /// The requested method requires verified symmetry.
    SymmetryRequired,
    /// CG requires verified positive definiteness and trivial nullspace.
    PositiveDefiniteRequired,
    /// MINRES is reserved for verified symmetric-indefinite systems.
    IndefiniteRequired,
    /// The preconditioner class is incompatible with the method.
    PreconditionerIncompatible {
        /// Requested solver.
        solver: LinearSolverKind,
        /// Observed preconditioner class.
        preconditioner: PreconditionerClass,
    },
}

impl fmt::Display for SolverAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceNotCompatible => f.write_str("linear source compatibility is not verified"),
            Self::NullspaceUnresolved => {
                f.write_str("linear operator nullspace remains unresolved")
            }
            Self::SymmetryRequired => f.write_str("requested solver requires verified symmetry"),
            Self::PositiveDefiniteRequired => {
                f.write_str("CG requires verified positive definiteness and trivial nullspace")
            }
            Self::IndefiniteRequired => {
                f.write_str("MINRES requires a verified symmetric-indefinite operator")
            }
            Self::PreconditionerIncompatible {
                solver,
                preconditioner,
            } => write!(
                f,
                "{solver:?} does not admit {preconditioner:?} preconditioner behavior"
            ),
        }
    }
}

impl core::error::Error for SolverAdmissionError {}

/// Solver-choice decision receipt for one local verification finding.
///
/// Raw Krylov state constructors remain available for low-level numerical
/// kernels. Callers that require authorization must bind this receipt to the
/// exact operator/source identity in their ledger before execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedLinearSolver {
    kind: LinearSolverKind,
    verified: VerifiedLinearSystem,
}

impl AdmittedLinearSolver {
    /// Admitted algorithm.
    #[must_use]
    pub const fn kind(&self) -> LinearSolverKind {
        self.kind
    }

    /// Verification wrapper used for admission.
    #[must_use]
    pub const fn verified(&self) -> &VerifiedLinearSystem {
        &self.verified
    }
}

/// Apply algorithm-specific admission rules to a retained finding.
#[must_use]
pub fn admit_linear_solver(
    kind: LinearSolverKind,
    verified: VerifiedLinearSystem,
) -> Result<AdmittedLinearSolver, SolverAdmissionError> {
    let finding = verified.finding;
    if finding.source != SourceCompatibility::Compatible {
        return Err(SolverAdmissionError::SourceNotCompatible);
    }
    if finding.nullspace == NullspaceEvidence::Unresolved {
        return Err(SolverAdmissionError::NullspaceUnresolved);
    }
    match kind {
        LinearSolverKind::Cg => {
            if finding.symmetry != SymmetryEvidence::Symmetric {
                return Err(SolverAdmissionError::SymmetryRequired);
            }
            if finding.definiteness != DefinitenessEvidence::PositiveDefinite
                || finding.nullspace != NullspaceEvidence::Trivial
            {
                return Err(SolverAdmissionError::PositiveDefiniteRequired);
            }
            if !matches!(
                finding.preconditioner,
                PreconditionerClass::None | PreconditionerClass::FixedSpd
            ) {
                return Err(SolverAdmissionError::PreconditionerIncompatible {
                    solver: kind,
                    preconditioner: finding.preconditioner,
                });
            }
        }
        LinearSolverKind::Minres => {
            if finding.symmetry != SymmetryEvidence::Symmetric {
                return Err(SolverAdmissionError::SymmetryRequired);
            }
            if finding.definiteness != DefinitenessEvidence::Indefinite {
                return Err(SolverAdmissionError::IndefiniteRequired);
            }
            if !matches!(
                finding.preconditioner,
                PreconditionerClass::None | PreconditionerClass::FixedSpd
            ) {
                return Err(SolverAdmissionError::PreconditionerIncompatible {
                    solver: kind,
                    preconditioner: finding.preconditioner,
                });
            }
        }
        LinearSolverKind::Gmres => {
            if finding.preconditioner == PreconditionerClass::Variable {
                return Err(SolverAdmissionError::PreconditionerIncompatible {
                    solver: kind,
                    preconditioner: finding.preconditioner,
                });
            }
        }
        LinearSolverKind::Fgmres => {}
    }
    Ok(AdmittedLinearSolver { kind, verified })
}

/// A deterministic flexible preconditioner keyed by logical inner iteration.
/// Implementations must be pure functions of `(logical_iteration, residual)`
/// for checkpoint/replay equality.
pub trait FlexiblePreconditioner {
    /// Apply the selected inverse approximation.
    fn apply(&self, logical_iteration: usize, residual: &[f64], output: &mut [f64]);
}

impl<T> FlexiblePreconditioner for T
where
    T: Precond + ?Sized,
{
    fn apply(&self, _logical_iteration: usize, residual: &[f64], output: &mut [f64]) {
        Precond::apply(self, residual, output);
    }
}

/// Resumable flexible GMRES. Checkpoints occur at restart-cycle boundaries;
/// `clone()` is the complete state.
///
/// Caller invariant: every resume must use the exact operator, right-hand
/// side, and logical-iteration preconditioner policy used to construct the
/// state. The opaque operator/RHS identity is not retained or authenticated.
#[derive(Debug, Clone)]
pub struct FgmresState {
    /// Current iterate.
    pub x: Vec<f64>,
    /// Restart length.
    pub restart: usize,
    bnorm: f64,
    rel: f64,
    /// Inner iterations completed across resumes.
    pub iters: usize,
    /// True relative residual after each completed cycle.
    pub history: Vec<f64>,
}

impl FgmresState {
    /// Start from the zero vector.
    #[must_use]
    pub fn new(rhs: &[f64], restart: usize) -> Self {
        assert!(restart >= 1, "FGMRES restart length must be positive");
        Self {
            x: vec![0.0; rhs.len()],
            restart,
            bnorm: norm2(rhs).max(f64::MIN_POSITIVE),
            rel: 1.0,
            iters: 0,
            history: Vec::new(),
        }
    }

    /// Current true relative residual from the last completed cycle.
    #[must_use]
    pub fn rel_residual(&self) -> f64 {
        self.rel
    }

    /// Run up to `max_cycles` additional restart cycles.
    #[allow(clippy::too_many_lines)] // The Arnoldi/Givens cycle is one invariant.
    pub fn run<A: LinearOp, P: FlexiblePreconditioner>(
        &mut self,
        operator: &A,
        preconditioner: &P,
        rhs: &[f64],
        tolerance: f64,
        max_cycles: usize,
    ) -> SolveReport {
        let n = operator.n();
        assert_eq!(rhs.len(), n, "FGMRES rhs length mismatch");
        assert_eq!(self.x.len(), n, "FGMRES checkpoint dimension mismatch");
        let mut broken = false;
        let mut applied = vec![0.0; n];
        for _ in 0..max_cycles {
            operator.apply(&self.x, &mut applied);
            let mut residual: Vec<f64> = rhs
                .iter()
                .zip(&applied)
                .map(|(right, left)| right - left)
                .collect();
            let beta = norm2(&residual);
            self.rel = beta / self.bnorm;
            if self.rel < tolerance {
                break;
            }
            if !beta.is_finite() || beta == 0.0 {
                broken = true;
                break;
            }
            for value in &mut residual {
                *value /= beta;
            }

            let m = self.restart;
            let mut basis = vec![residual];
            let mut preconditioned_basis: Vec<Vec<f64>> = Vec::with_capacity(m);
            let mut h = vec![0.0; (m + 1) * m];
            let mut cs = vec![0.0; m];
            let mut sn = vec![0.0; m];
            let mut g = vec![0.0; m + 1];
            g[0] = beta;
            let mut columns = 0usize;
            for column in 0..m {
                let mut z = vec![0.0; n];
                preconditioner.apply(self.iters, &basis[column], &mut z);
                if z.iter().any(|value| !value.is_finite()) {
                    broken = true;
                    break;
                }
                let mut w = vec![0.0; n];
                operator.apply(&z, &mut w);
                for (row, vector) in basis.iter().enumerate() {
                    let coefficient = dot(vector, &w);
                    h[row * m + column] = coefficient;
                    for (value, basis_value) in w.iter_mut().zip(vector) {
                        *value = coefficient.mul_add(-basis_value, *value);
                    }
                }
                let next_norm = norm2(&w);
                h[(column + 1) * m + column] = next_norm;
                for row in 0..column {
                    let upper = h[row * m + column];
                    let lower = h[(row + 1) * m + column];
                    h[row * m + column] = cs[row].mul_add(upper, sn[row] * lower);
                    h[(row + 1) * m + column] = (-sn[row]).mul_add(upper, cs[row] * lower);
                }
                let diagonal = h[column * m + column];
                let denominator = diagonal.hypot(next_norm);
                if !denominator.is_finite() || denominator == 0.0 {
                    broken = true;
                    break;
                }
                cs[column] = diagonal / denominator;
                sn[column] = next_norm / denominator;
                h[column * m + column] = denominator;
                h[(column + 1) * m + column] = 0.0;
                g[column + 1] = -sn[column] * g[column];
                g[column] *= cs[column];
                preconditioned_basis.push(z);
                columns = column + 1;
                self.iters += 1;
                if next_norm == 0.0 || g[column + 1].abs() / self.bnorm < tolerance {
                    break;
                }
                for value in &mut w {
                    *value /= next_norm;
                }
                basis.push(w);
            }
            if broken || columns == 0 {
                break;
            }
            let mut coefficients = vec![0.0; columns];
            for row in (0..columns).rev() {
                let mut value = g[row];
                for column in (row + 1)..columns {
                    value = h[row * m + column].mul_add(-coefficients[column], value);
                }
                let diagonal = h[row * m + row];
                if !diagonal.is_finite() || diagonal == 0.0 {
                    broken = true;
                    break;
                }
                coefficients[row] = value / diagonal;
            }
            if broken {
                break;
            }
            for (coefficient, vector) in coefficients.iter().zip(&preconditioned_basis) {
                for (value, direction) in self.x.iter_mut().zip(vector) {
                    *value = coefficient.mul_add(*direction, *value);
                }
            }
            operator.apply(&self.x, &mut applied);
            let true_residual: Vec<f64> = rhs
                .iter()
                .zip(&applied)
                .map(|(right, left)| right - left)
                .collect();
            self.rel = norm2(&true_residual) / self.bnorm;
            self.history.push(self.rel);
            if self.rel < tolerance {
                break;
            }
        }
        SolveReport {
            iters: self.iters,
            rel_residual: self.rel,
            converged: self.rel < tolerance,
            history: self.history.clone(),
            diagnosis: if self.rel < tolerance {
                None
            } else if broken || !self.rel.is_finite() {
                Some(StallDiagnosis::Breakdown)
            } else {
                Some(StallDiagnosis::BudgetExhausted)
            },
        }
    }
}

/// Nonlinear residual and Jacobian actions at an exact state.
pub trait NonlinearProblem {
    /// State/residual dimension.
    fn dimension(&self) -> usize;
    /// Overwrite `residual` with `F(x)`.
    fn residual(&self, x: &[f64], residual: &mut [f64]);
    /// Overwrite `output` with `J(x) direction`.
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]);
}

/// Backtracking line-search controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineSearchConfig {
    /// Armijo sufficient-decrease coefficient in `(0, 1)`.
    pub armijo: f64,
    /// Backtracking multiplier in `(0, 1)`.
    pub contraction: f64,
    /// Smallest admitted step length in `(0, 1]`.
    pub minimum_step: f64,
}

impl Default for LineSearchConfig {
    fn default() -> Self {
        Self {
            armijo: 1.0e-4,
            contraction: 0.5,
            minimum_step: 1.0 / 1024.0,
        }
    }
}

/// Trust-region controls.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrustRegionConfig {
    /// Initial radius.
    pub initial_radius: f64,
    /// Minimum radius before terminal refusal.
    pub minimum_radius: f64,
    /// Maximum radius.
    pub maximum_radius: f64,
    /// Minimum actual/predicted reduction ratio for acceptance.
    pub acceptance_ratio: f64,
    /// Ratio above which a boundary step expands the radius.
    pub expansion_ratio: f64,
    /// Rejection shrink factor in `(0, 1)`.
    pub shrink: f64,
    /// Acceptance expansion factor greater than one.
    pub expansion: f64,
}

impl Default for TrustRegionConfig {
    fn default() -> Self {
        Self {
            initial_radius: 1.0,
            minimum_radius: 1.0e-12,
            maximum_radius: 1.0e6,
            acceptance_ratio: 0.1,
            expansion_ratio: 0.75,
            shrink: 0.25,
            expansion: 2.0,
        }
    }
}

/// Nonlinear globalization policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Globalization {
    /// Armijo backtracking.
    LineSearch(LineSearchConfig),
    /// Radius-limited Newton step with actual/predicted reduction checks.
    TrustRegion(TrustRegionConfig),
}

/// Newton--Krylov configuration retained in every checkpoint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NewtonKrylovConfig {
    /// Absolute residual tolerance.
    pub absolute_tolerance: f64,
    /// Relative residual tolerance against the initial norm.
    pub relative_tolerance: f64,
    /// FGMRES restart length.
    pub linear_restart: usize,
    /// Maximum FGMRES cycles per outer iteration.
    pub max_linear_cycles: usize,
    /// Minimum Eisenstat--Walker forcing term.
    pub forcing_minimum: f64,
    /// Maximum Eisenstat--Walker forcing term.
    pub forcing_maximum: f64,
    /// Eisenstat--Walker multiplier.
    pub forcing_gamma: f64,
    /// Eisenstat--Walker exponent.
    pub forcing_exponent: f64,
    /// Globalization policy.
    pub globalization: Globalization,
}

impl Default for NewtonKrylovConfig {
    fn default() -> Self {
        Self {
            absolute_tolerance: 1.0e-12,
            relative_tolerance: 1.0e-10,
            linear_restart: 24,
            max_linear_cycles: 8,
            forcing_minimum: 1.0e-10,
            forcing_maximum: 0.5,
            forcing_gamma: 0.9,
            forcing_exponent: 1.5,
            globalization: Globalization::LineSearch(LineSearchConfig::default()),
        }
    }
}

/// One retained outer-iteration decision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GlobalizationDecision {
    /// Line search accepted a step.
    LineSearchAccepted,
    /// Line search exhausted its minimum-step budget.
    LineSearchRejected,
    /// Trust region accepted a step.
    TrustRegionAccepted,
    /// Trust region rejected a step and shrank the radius.
    TrustRegionRejected,
    /// A trial/model quantity was non-finite, so no state update occurred.
    NonFiniteRejected,
}

/// Complete deterministic telemetry for one Newton attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct NewtonIteration {
    /// Zero-based outer iteration.
    pub iteration: usize,
    /// Residual norm before the attempt.
    pub residual_before: f64,
    /// Residual norm after the attempt (unchanged on rejection).
    pub residual_after: f64,
    /// Eisenstat--Walker forcing term.
    pub forcing: f64,
    /// Inner true relative residual.
    pub linear_relative_residual: f64,
    /// Inner iterations.
    pub linear_iterations: usize,
    /// Norm of the unscaled Newton direction.
    pub newton_step_norm: f64,
    /// Accepted/scaled step length relative to the Newton direction.
    pub step_length: f64,
    /// Trust actual/predicted ratio, or `None` for line search/an invalid
    /// nonpositive trust-model prediction.
    pub reduction_ratio: Option<f64>,
    /// Globalization decision.
    pub decision: GlobalizationDecision,
}

/// Typed nonlinear stall/refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonStallDiagnosis {
    /// Caller-provided iteration budget ended while progress remained possible.
    BudgetExhausted,
    /// Accepted residuals remained within five percent over eight accepted steps.
    Plateau,
    /// Inner FGMRES broke down or exhausted its forcing budget.
    LinearSolveFailed(StallDiagnosis),
    /// Backtracking crossed the minimum step without sufficient decrease.
    LineSearchRejected,
    /// Repeated trust rejection exhausted the radius budget.
    TrustRegionRadiusExhausted,
    /// A residual or Jacobian-derived quantity became non-finite.
    NonFinite,
}

/// Construction error for a Newton checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NewtonError {
    /// State and problem dimensions differ.
    Dimension {
        /// Problem dimension.
        problem: usize,
        /// State length.
        state: usize,
    },
    /// Configuration is inconsistent or non-finite.
    InvalidConfig(&'static str),
    /// A state entry is NaN or infinite.
    NonFiniteState {
        /// Outer iteration.
        iteration: usize,
        /// Entry index.
        index: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// A residual entry is NaN or infinite.
    NonFiniteResidual {
        /// Outer iteration.
        iteration: usize,
        /// Entry index.
        index: usize,
        /// Exact refused bits.
        bits: u64,
    },
    /// Finite residual entries had a norm that was not representable.
    NonFiniteResidualNorm {
        /// Outer iteration.
        iteration: usize,
    },
}

impl fmt::Display for NewtonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Dimension { problem, state } => write!(
                f,
                "nonlinear problem dimension {problem} differs from state {state}"
            ),
            Self::InvalidConfig(reason) => write!(f, "invalid Newton--Krylov config: {reason}"),
            Self::NonFiniteState {
                iteration,
                index,
                bits,
            } => write!(
                f,
                "Newton state iteration {iteration} entry {index} is non-finite \
                 (bits 0x{bits:016x})"
            ),
            Self::NonFiniteResidual {
                iteration,
                index,
                bits,
            } => write!(
                f,
                "Newton residual iteration {iteration} entry {index} is non-finite \
                 (bits 0x{bits:016x})"
            ),
            Self::NonFiniteResidualNorm { iteration } => write!(
                f,
                "Newton residual norm at iteration {iteration} is not finite"
            ),
        }
    }
}

impl core::error::Error for NewtonError {}

/// Resumable Newton--Krylov state. Every field needed for the next outer
/// attempt is plain cloneable data.
///
/// Caller invariant: every resume must use the exact problem implementation
/// and parameterization used at construction. Same-dimension problem identity
/// is not retained or authenticated by this lower-layer numerical state.
#[derive(Debug, Clone)]
pub struct NewtonKrylovState {
    /// Current state vector.
    pub x: Vec<f64>,
    residual: Vec<f64>,
    residual_norm: f64,
    initial_residual_norm: f64,
    previous_residual_norm: f64,
    trust_radius: f64,
    config: NewtonKrylovConfig,
    terminal: Option<NewtonStallDiagnosis>,
    /// Outer attempts completed across resumes.
    pub iterations: usize,
    /// Complete outer telemetry.
    pub history: Vec<NewtonIteration>,
}

/// Snapshot report from a Newton state.
#[derive(Debug, Clone, PartialEq)]
pub struct NewtonReport {
    /// Outer attempts completed.
    pub iterations: usize,
    /// Current residual norm.
    pub residual_norm: f64,
    /// Absolute/relative target met.
    pub converged: bool,
    /// Full retained iteration telemetry.
    pub history: Vec<NewtonIteration>,
    /// Present iff not converged.
    pub diagnosis: Option<NewtonStallDiagnosis>,
}

struct JacobianAt<'a, P> {
    problem: &'a P,
    x: &'a [f64],
}

impl<P: NonlinearProblem> LinearOp for JacobianAt<'_, P> {
    fn n(&self) -> usize {
        self.problem.dimension()
    }

    fn apply(&self, direction: &[f64], output: &mut [f64]) {
        self.problem.jacobian_apply(self.x, direction, output);
    }
}

#[derive(Debug, Clone, Copy)]
struct IdentityFlexible;

impl FlexiblePreconditioner for IdentityFlexible {
    fn apply(&self, _logical_iteration: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

impl NewtonKrylovState {
    /// Evaluate and retain the initial residual.
    #[must_use]
    pub fn new<P: NonlinearProblem>(
        problem: &P,
        x: Vec<f64>,
        config: NewtonKrylovConfig,
    ) -> Result<Self, NewtonError> {
        validate_newton_config(config)?;
        if problem.dimension() != x.len() {
            return Err(NewtonError::Dimension {
                problem: problem.dimension(),
                state: x.len(),
            });
        }
        let (residual, residual_norm) = evaluate_residual(problem, &x, 0)?;
        let trust_radius = match config.globalization {
            Globalization::LineSearch(_) => 0.0,
            Globalization::TrustRegion(trust) => trust.initial_radius,
        };
        Ok(Self {
            x,
            residual,
            residual_norm,
            initial_residual_norm: residual_norm.max(f64::MIN_POSITIVE),
            previous_residual_norm: residual_norm,
            trust_radius,
            config,
            terminal: None,
            iterations: 0,
            history: Vec::new(),
        })
    }

    /// Current residual norm.
    #[must_use]
    pub fn residual_norm(&self) -> f64 {
        self.residual_norm
    }

    /// Retained configuration.
    #[must_use]
    pub const fn config(&self) -> NewtonKrylovConfig {
        self.config
    }

    /// Run up to `max_iterations` additional outer attempts.
    pub fn run<P: NonlinearProblem>(&mut self, problem: &P, max_iterations: usize) -> NewtonReport {
        assert_eq!(
            problem.dimension(),
            self.x.len(),
            "Newton checkpoint/problem dimension mismatch"
        );
        for _ in 0..max_iterations {
            if self.is_converged() || self.terminal.is_some() {
                break;
            }
            self.step(problem);
        }
        self.report()
    }

    fn target(&self) -> f64 {
        self.config
            .absolute_tolerance
            .max(self.config.relative_tolerance * self.initial_residual_norm)
    }

    fn is_converged(&self) -> bool {
        self.residual_norm() <= self.target()
    }

    fn record_nonfinite_attempt(
        &mut self,
        iteration: usize,
        residual_before: f64,
        forcing: f64,
        linear_report: &SolveReport,
        direction_norm: f64,
        step_length: f64,
    ) {
        self.history.push(NewtonIteration {
            iteration,
            residual_before,
            residual_after: residual_before,
            forcing,
            linear_relative_residual: linear_report.rel_residual,
            linear_iterations: linear_report.iters,
            newton_step_norm: direction_norm,
            step_length,
            reduction_ratio: None,
            decision: GlobalizationDecision::NonFiniteRejected,
        });
        self.iterations += 1;
        self.terminal = Some(NewtonStallDiagnosis::NonFinite);
    }

    #[allow(clippy::too_many_lines)] // One atomic outer attempt and its receipt.
    fn step<P: NonlinearProblem>(&mut self, problem: &P) {
        let residual_before = self.residual_norm();
        if !residual_before.is_finite() {
            self.terminal = Some(NewtonStallDiagnosis::NonFinite);
            return;
        }
        let forcing = if self.iterations == 0 {
            self.config.forcing_maximum
        } else {
            let ratio = residual_before / self.previous_residual_norm.max(f64::MIN_POSITIVE);
            (self.config.forcing_gamma * fs_math::det::pow(ratio, self.config.forcing_exponent))
                .clamp(self.config.forcing_minimum, self.config.forcing_maximum)
        };
        let rhs: Vec<f64> = self.residual.iter().map(|value| -*value).collect();
        let jacobian = JacobianAt {
            problem,
            x: &self.x,
        };
        let mut linear = FgmresState::new(&rhs, self.config.linear_restart);
        let linear_report = linear.run(
            &jacobian,
            &IdentityFlexible,
            &rhs,
            forcing,
            self.config.max_linear_cycles,
        );
        if !linear_report.converged {
            self.terminal = Some(NewtonStallDiagnosis::LinearSolveFailed(
                linear_report
                    .diagnosis
                    .unwrap_or(StallDiagnosis::BudgetExhausted),
            ));
            return;
        }
        let direction = linear.x;
        let Some(direction_norm) = finite_norm(&direction) else {
            self.terminal = Some(NewtonStallDiagnosis::NonFinite);
            return;
        };

        let iteration = self.iterations;
        let outcome = match self.config.globalization {
            Globalization::LineSearch(line) => {
                let mut step_length = 1.0;
                let mut last_step_length = step_length;
                let mut accepted = None;
                let mut saw_finite_trial = false;
                while step_length >= line.minimum_step {
                    last_step_length = step_length;
                    let trial_x = stepped(&self.x, &direction, step_length);
                    match evaluate_residual(problem, &trial_x, iteration + 1) {
                        Ok((trial_residual, trial_norm)) => {
                            saw_finite_trial = true;
                            if trial_norm <= (1.0 - line.armijo * step_length) * residual_before {
                                accepted = Some((trial_x, trial_residual, trial_norm));
                                break;
                            }
                        }
                        Err(_) => {}
                    }
                    step_length *= line.contraction;
                }
                if let Some((trial_x, trial_residual, trial_norm)) = accepted {
                    Some((
                        trial_x,
                        trial_residual,
                        trial_norm,
                        step_length,
                        None,
                        GlobalizationDecision::LineSearchAccepted,
                    ))
                } else {
                    self.history.push(NewtonIteration {
                        iteration,
                        residual_before,
                        residual_after: residual_before,
                        forcing,
                        linear_relative_residual: linear_report.rel_residual,
                        linear_iterations: linear_report.iters,
                        newton_step_norm: direction_norm,
                        step_length: last_step_length,
                        reduction_ratio: None,
                        decision: if saw_finite_trial {
                            GlobalizationDecision::LineSearchRejected
                        } else {
                            GlobalizationDecision::NonFiniteRejected
                        },
                    });
                    self.iterations += 1;
                    self.terminal = Some(if saw_finite_trial {
                        NewtonStallDiagnosis::LineSearchRejected
                    } else {
                        NewtonStallDiagnosis::NonFinite
                    });
                    None
                }
            }
            Globalization::TrustRegion(trust) => {
                let step_length = if direction_norm > self.trust_radius {
                    self.trust_radius / direction_norm
                } else {
                    1.0
                };
                let trial_x = stepped(&self.x, &direction, step_length);
                let (trial_residual, trial_norm) =
                    match evaluate_residual(problem, &trial_x, iteration + 1) {
                        Ok(evaluated) => evaluated,
                        Err(_) => {
                            self.record_nonfinite_attempt(
                                iteration,
                                residual_before,
                                forcing,
                                &linear_report,
                                direction_norm,
                                step_length,
                            );
                            return;
                        }
                    };
                let mut jacobian_step = vec![0.0; direction.len()];
                let scaled_direction: Vec<f64> =
                    direction.iter().map(|value| step_length * value).collect();
                problem.jacobian_apply(&self.x, &scaled_direction, &mut jacobian_step);
                if finite_norm(&jacobian_step).is_none() {
                    self.record_nonfinite_attempt(
                        iteration,
                        residual_before,
                        forcing,
                        &linear_report,
                        direction_norm,
                        step_length,
                    );
                    return;
                }
                let linearized: Vec<f64> = self
                    .residual
                    .iter()
                    .zip(jacobian_step)
                    .map(|(residual, change)| residual + change)
                    .collect();
                let Some(linearized_norm) = finite_norm(&linearized) else {
                    self.record_nonfinite_attempt(
                        iteration,
                        residual_before,
                        forcing,
                        &linear_report,
                        direction_norm,
                        step_length,
                    );
                    return;
                };
                let predicted = residual_before - linearized_norm;
                let actual = residual_before - trial_norm;
                if !predicted.is_finite() || !actual.is_finite() {
                    self.record_nonfinite_attempt(
                        iteration,
                        residual_before,
                        forcing,
                        &linear_report,
                        direction_norm,
                        step_length,
                    );
                    return;
                }
                let ratio = if predicted > 0.0 {
                    let ratio = actual / predicted;
                    if !ratio.is_finite() {
                        self.record_nonfinite_attempt(
                            iteration,
                            residual_before,
                            forcing,
                            &linear_report,
                            direction_norm,
                            step_length,
                        );
                        return;
                    }
                    Some(ratio)
                } else {
                    None
                };
                if actual > 0.0 && ratio.is_some_and(|value| value >= trust.acceptance_ratio) {
                    if ratio.is_some_and(|value| value >= trust.expansion_ratio)
                        && step_length < 1.0
                    {
                        self.trust_radius =
                            (self.trust_radius * trust.expansion).min(trust.maximum_radius);
                    }
                    Some((
                        trial_x,
                        trial_residual,
                        trial_norm,
                        step_length,
                        ratio,
                        GlobalizationDecision::TrustRegionAccepted,
                    ))
                } else {
                    self.trust_radius *= trust.shrink;
                    self.history.push(NewtonIteration {
                        iteration,
                        residual_before,
                        residual_after: residual_before,
                        forcing,
                        linear_relative_residual: linear_report.rel_residual,
                        linear_iterations: linear_report.iters,
                        newton_step_norm: direction_norm,
                        step_length,
                        reduction_ratio: ratio,
                        decision: GlobalizationDecision::TrustRegionRejected,
                    });
                    self.iterations += 1;
                    if self.trust_radius < trust.minimum_radius {
                        self.terminal = Some(NewtonStallDiagnosis::TrustRegionRadiusExhausted);
                    }
                    None
                }
            }
        };

        if let Some((trial_x, trial_residual, trial_norm, step_length, ratio, decision)) = outcome {
            self.previous_residual_norm = residual_before;
            self.x = trial_x;
            self.residual = trial_residual;
            self.residual_norm = trial_norm;
            self.history.push(NewtonIteration {
                iteration,
                residual_before,
                residual_after: trial_norm,
                forcing,
                linear_relative_residual: linear_report.rel_residual,
                linear_iterations: linear_report.iters,
                newton_step_norm: direction_norm,
                step_length,
                reduction_ratio: ratio,
                decision,
            });
            self.iterations += 1;
        }
    }

    fn report(&self) -> NewtonReport {
        let residual_norm = self.residual_norm();
        let converged = self.is_converged();
        let diagnosis = if converged {
            None
        } else if let Some(terminal) = self.terminal {
            Some(terminal)
        } else if nonlinear_plateau(&self.history) {
            Some(NewtonStallDiagnosis::Plateau)
        } else {
            Some(NewtonStallDiagnosis::BudgetExhausted)
        };
        NewtonReport {
            iterations: self.iterations,
            residual_norm,
            converged,
            history: self.history.clone(),
            diagnosis,
        }
    }
}

fn stepped(state: &[f64], direction: &[f64], length: f64) -> Vec<f64> {
    state
        .iter()
        .zip(direction)
        .map(|(value, step)| length.mul_add(*step, *value))
        .collect()
}

fn evaluate_residual<P: NonlinearProblem>(
    problem: &P,
    state: &[f64],
    iteration: usize,
) -> Result<(Vec<f64>, f64), NewtonError> {
    for (index, value) in state.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(NewtonError::NonFiniteState {
                iteration,
                index,
                bits: value.to_bits(),
            });
        }
    }
    let mut residual = vec![0.0; problem.dimension()];
    problem.residual(state, &mut residual);
    for (index, value) in residual.iter().copied().enumerate() {
        if !value.is_finite() {
            return Err(NewtonError::NonFiniteResidual {
                iteration,
                index,
                bits: value.to_bits(),
            });
        }
    }
    let Some(residual_norm) = finite_norm(&residual) else {
        return Err(NewtonError::NonFiniteResidualNorm { iteration });
    };
    Ok((residual, residual_norm))
}

fn finite_norm(values: &[f64]) -> Option<f64> {
    if values.iter().any(|value| !value.is_finite()) {
        return None;
    }
    let norm = norm2(values);
    norm.is_finite().then_some(norm)
}

fn nonlinear_plateau(history: &[NewtonIteration]) -> bool {
    let accepted: Vec<&NewtonIteration> = history
        .iter()
        .filter(|entry| {
            matches!(
                entry.decision,
                GlobalizationDecision::LineSearchAccepted
                    | GlobalizationDecision::TrustRegionAccepted
            )
        })
        .rev()
        .take(8)
        .collect();
    if accepted.len() < 8 {
        return false;
    }
    let first = accepted[7].residual_before;
    let last = accepted[0].residual_after;
    first.is_finite() && last >= 0.95 * first && last <= 1.05 * first
}

fn validate_newton_config(config: NewtonKrylovConfig) -> Result<(), NewtonError> {
    if !config.absolute_tolerance.is_finite() || config.absolute_tolerance < 0.0 {
        return Err(NewtonError::InvalidConfig(
            "absolute_tolerance must be finite and nonnegative",
        ));
    }
    if !config.relative_tolerance.is_finite()
        || config.relative_tolerance <= 0.0
        || config.relative_tolerance > 1.0
    {
        return Err(NewtonError::InvalidConfig(
            "relative_tolerance must be finite and in (0,1]",
        ));
    }
    if config.linear_restart == 0 || config.max_linear_cycles == 0 {
        return Err(NewtonError::InvalidConfig(
            "linear restart and cycle budgets must be positive",
        ));
    }
    if !config.forcing_minimum.is_finite()
        || !config.forcing_maximum.is_finite()
        || config.forcing_minimum <= 0.0
        || config.forcing_minimum > config.forcing_maximum
        || config.forcing_maximum >= 1.0
    {
        return Err(NewtonError::InvalidConfig(
            "forcing terms require 0 < minimum <= maximum < 1",
        ));
    }
    if !config.forcing_gamma.is_finite()
        || config.forcing_gamma <= 0.0
        || !config.forcing_exponent.is_finite()
        || config.forcing_exponent <= 0.0
    {
        return Err(NewtonError::InvalidConfig(
            "forcing gamma and exponent must be finite and positive",
        ));
    }
    match config.globalization {
        Globalization::LineSearch(line) => {
            if !line.armijo.is_finite()
                || line.armijo <= 0.0
                || line.armijo >= 1.0
                || !line.contraction.is_finite()
                || line.contraction <= 0.0
                || line.contraction >= 1.0
                || !line.minimum_step.is_finite()
                || line.minimum_step <= 0.0
                || line.minimum_step > 1.0
            {
                return Err(NewtonError::InvalidConfig(
                    "line search requires Armijo/contraction in (0,1) and minimum step in (0,1]",
                ));
            }
        }
        Globalization::TrustRegion(trust) => {
            if !trust.initial_radius.is_finite()
                || !trust.minimum_radius.is_finite()
                || !trust.maximum_radius.is_finite()
                || trust.minimum_radius <= 0.0
                || trust.initial_radius < trust.minimum_radius
                || trust.initial_radius > trust.maximum_radius
                || !trust.acceptance_ratio.is_finite()
                || trust.acceptance_ratio <= 0.0
                || trust.acceptance_ratio >= 1.0
                || !trust.expansion_ratio.is_finite()
                || trust.expansion_ratio <= trust.acceptance_ratio
                || trust.expansion_ratio >= 1.0
                || !trust.shrink.is_finite()
                || trust.shrink <= 0.0
                || trust.shrink >= 1.0
                || !trust.expansion.is_finite()
                || trust.expansion <= 1.0
            {
                return Err(NewtonError::InvalidConfig(
                    "trust radii/ratios/shrink/expansion are inconsistent",
                ));
            }
        }
    }
    Ok(())
}
