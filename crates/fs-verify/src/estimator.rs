//! The VERIFIER: equilibrated-flux a-posteriori bounds (Prager–Synge,
//! 1D elliptic class), interval-evaluated to VERIFIED color.
//!
//! The rigor structure: ANY σ with `σ′ = −f` yields the guaranteed
//! bound `‖(u − u_h)′‖ ≤ ‖σ − u_h′‖` — so the free constant in
//! `σ = c − F` is optimized in plain f64 for TIGHTNESS while the bound
//! itself is evaluated with outward-rounded intervals over exact Gauss
//! quadrature (polynomial data ⇒ the quadrature identity is exact;
//! only rounding needs enclosing). Malformed inputs and unusable
//! enclosures FAIL CLOSED as structured refusals: no color, ever.

use crate::fem1d::{
    Fem1dError, MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS, MmsProblem, gauss5,
    require_converged, true_energy_error, try_zeroed, validate_candidate, validate_problem,
};
use crate::interval::Iv;
use fs_evidence::Color;
use std::fmt::Write as _;

/// Largest mesh admitted by the synchronous v0 verifier.
pub const MAX_VERIFIER_MESH_NODES: usize = MAX_FEM1D_MESH_NODES;
/// Exactness envelope for the manufactured solution: degree at most five.
pub const MAX_VERIFIER_POLY_COEFFICIENTS: usize = MAX_FEM1D_POLY_COEFFICIENTS;
/// Semantic version of the verifier's complete bounded-work accounting.
pub const VERIFIER_WORK_PLAN_VERSION: u32 = 1;
/// Semantic version of the verifier's callback/checkpoint schedule.
pub const VERIFIER_POLL_POLICY_VERSION: u32 = 1;
/// Maximum completed logical work between verifier work-boundary callbacks.
pub const VERIFIER_POLL_STRIDE_WORK_UNITS: u128 = 256;

/// One phase of the bounded equilibrated-flux verification workflow.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VerifierPhase {
    /// Input and derived-polynomial validation.
    Validation,
    /// Rounded optimizer for the free equilibrated-flux constant.
    Tightness,
    /// Outward-rounded equilibrated-bound construction.
    Equilibrated,
    /// Deterministic reconstructed-flux identity.
    Hash,
    /// Final report construction and publication.
    Finalization,
}

/// Why a verifier progress callback is being invoked.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VerifierCheckpointKind {
    /// Mandatory callback before a phase begins.
    PhaseEntry,
    /// Invocation-global multiple of [`VERIFIER_POLL_STRIDE_WORK_UNITS`].
    WorkBoundary,
    /// Mandatory callback after a structured refusal has inspected all work it reports.
    RefusalFlush,
    /// Mandatory final callback after the complete report is ready to publish.
    Publication,
}

/// Immutable invocation-global verifier progress passed to a callback by value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VerifierProgress {
    /// Callback reason.
    pub kind: VerifierCheckpointKind,
    /// Phase active at this callback.
    pub phase: VerifierPhase,
    /// Exact logical work completed across all phases.
    pub completed_work_units: u128,
    /// Complete constant-time preflighted work for this invocation.
    pub planned_work_units: u128,
}

/// Complete checked logical-work shape for one verifier invocation.
///
/// Counts are exact credited logical progress, not instruction counters. The
/// fixed-cap polynomial helpers process at most six coefficients atomically and
/// are credited only after success; this can conservatively lag physical work
/// by one such micro-tile but never overstates completed work.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VerifierWorkPlan {
    validation: u128,
    tightness: u128,
    equilibrated: u128,
    hash: u128,
    finalization: u128,
    total: u128,
}

impl VerifierWorkPlan {
    /// Preflight the complete work shape without running a callback.
    ///
    /// Shape refusals therefore have exact zero-work semantics. Content
    /// validation happens later through [`verify_with_checkpoint`].
    ///
    /// # Errors
    /// Returns a structured shape refusal for an inadmissible mesh, candidate,
    /// polynomial envelope, or checked work-plan overflow.
    pub fn for_inputs(problem: &MmsProblem, candidate: &[f64]) -> Result<Self, VerifierRefusal> {
        let mesh_nodes = problem.mesh().len();
        if !(2..=MAX_VERIFIER_MESH_NODES).contains(&mesh_nodes) {
            return Err(VerifierRefusal::MeshNodeCount);
        }
        if candidate.len() != mesh_nodes {
            return Err(VerifierRefusal::CandidateLength);
        }
        for (role, polynomial) in [
            (VerifierPolynomial::ExactSolution, problem.exact_solution()),
            (VerifierPolynomial::Forcing, problem.forcing()),
            (
                VerifierPolynomial::ForcingAntiderivative,
                problem.rounded_forcing_antiderivative(),
            ),
        ] {
            if !(1..=MAX_VERIFIER_POLY_COEFFICIENTS).contains(&polynomial.coefficients().len()) {
                return Err(VerifierRefusal::PolynomialCoefficientCount { polynomial: role });
            }
        }

        let mesh_nodes =
            u128::try_from(mesh_nodes).map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let cells = mesh_nodes
            .checked_sub(1)
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let exact_coefficients = u128::try_from(problem.exact_solution().coefficients().len())
            .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let forcing_coefficients = u128::try_from(problem.forcing().coefficients().len())
            .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;
        let antiderivative_coefficients = u128::try_from(
            problem
                .rounded_forcing_antiderivative()
                .coefficients()
                .len(),
        )
        .map_err(|_| VerifierRefusal::WorkPlanOverflow)?;

        let validation_work_units = 3_u128
            .checked_add(
                mesh_nodes
                    .checked_mul(2)
                    .ok_or(VerifierRefusal::WorkPlanOverflow)?,
            )
            .and_then(|work| work.checked_add(cells))
            .and_then(|work| work.checked_add(exact_coefficients.checked_mul(3)?))
            .and_then(|work| work.checked_add(forcing_coefficients.checked_mul(3)?))
            .and_then(|work| work.checked_add(antiderivative_coefficients.checked_mul(2)?))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let tightness_work_units = cells;
        let equilibrated_work_units = cells;
        let hash_work_units = 3_u128
            .checked_add(forcing_coefficients)
            .and_then(|work| work.checked_add(antiderivative_coefficients))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        let finalization_work_units = 1;
        let planned_work_units = validation_work_units
            .checked_add(tightness_work_units)
            .and_then(|work| work.checked_add(equilibrated_work_units))
            .and_then(|work| work.checked_add(hash_work_units))
            .and_then(|work| work.checked_add(finalization_work_units))
            .ok_or(VerifierRefusal::WorkPlanOverflow)?;
        Ok(Self {
            validation: validation_work_units,
            tightness: tightness_work_units,
            equilibrated: equilibrated_work_units,
            hash: hash_work_units,
            finalization: finalization_work_units,
            total: planned_work_units,
        })
    }

    /// Validation work in the plan.
    #[must_use]
    pub const fn validation_work_units(self) -> u128 {
        self.validation
    }

    /// Tightness-optimizer work in the plan.
    #[must_use]
    pub const fn tightness_work_units(self) -> u128 {
        self.tightness
    }

    /// Equilibrated-bound work in the plan.
    #[must_use]
    pub const fn equilibrated_work_units(self) -> u128 {
        self.equilibrated
    }

    /// Flux-identity work in the plan.
    #[must_use]
    pub const fn hash_work_units(self) -> u128 {
        self.hash
    }

    /// Final report/publication work in the plan.
    #[must_use]
    pub const fn finalization_work_units(self) -> u128 {
        self.finalization
    }

    /// Total work in the plan.
    #[must_use]
    pub const fn planned_work_units(self) -> u128 {
        self.total
    }

    /// Stable phase counts used by downstream evidence identities.
    #[must_use]
    pub const fn identity_fields(self) -> [u128; 6] {
        [
            self.validation,
            self.tightness,
            self.equilibrated,
            self.hash,
            self.finalization,
            self.total,
        ]
    }
}

/// Estimator families (Proposal D's independence escalation needs at
/// least two registered per class).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimatorFamily {
    /// Equilibrated flux (guaranteed, constant-free — the verifier).
    EquilibratedFlux,
    /// Hierarchical (refined-mesh comparison — independent, NOT
    /// guaranteed; the falsifier's cross-check).
    Hierarchical,
}

/// Polynomial role carried by a structured verifier refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierPolynomial {
    /// Manufactured exact-solution metadata (`u`).
    ExactSolution,
    /// Canonical forcing (`f = -u''`).
    Forcing,
    /// Canonical zero-constant antiderivative of the forcing (`big_f`).
    ForcingAntiderivative,
}

/// Stable, structured reason why no verifier authority was issued.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifierRefusal {
    /// Mesh length is outside the bounded `2..=MAX_VERIFIER_MESH_NODES` class.
    MeshNodeCount,
    /// One polynomial is empty or exceeds the degree-five exactness class.
    PolynomialCoefficientCount {
        /// Polynomial whose resource envelope was violated.
        polynomial: VerifierPolynomial,
    },
    /// Candidate and mesh lengths differ.
    CandidateLength,
    /// Tolerance is non-finite or non-positive.
    InvalidTolerance,
    /// Mesh endpoints are not canonical `+0.0` and `1.0`.
    MeshDomain,
    /// A mesh coordinate is non-finite or the mesh is not strictly increasing.
    MeshCoordinates,
    /// A candidate value is non-finite.
    CandidateNonFinite,
    /// Candidate endpoints are not canonical homogeneous `+0.0` values.
    CandidateBoundary,
    /// One polynomial contains a non-finite coefficient.
    PolynomialNonFinite {
        /// Polynomial containing the non-finite coefficient.
        polynomial: VerifierPolynomial,
    },
    /// The exact-solution polynomial does not vanish canonically at both ends.
    ExactSolutionBoundary,
    /// A public derived polynomial differs from the canonical value recomputed from `u`.
    DerivedPolynomialMismatch {
        /// Public derived polynomial that did not match its canonical value.
        polynomial: VerifierPolynomial,
    },
    /// The optional tightness constant could not be computed finitely.
    NonFiniteTightness,
    /// Interval construction produced a non-finite, reversed, or unusable enclosure.
    InvalidEnclosure,
    /// Complete verifier work-shape arithmetic overflowed.
    WorkPlanOverflow,
    /// Executed logical work did not match the complete preflighted plan.
    WorkPlanMismatch,
}

impl VerifierRefusal {
    /// Stable identifier for diagnostics and ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            Self::MeshNodeCount => "mesh-node-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-coefficient-count",
            Self::PolynomialCoefficientCount {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-coefficient-count",
            Self::CandidateLength => "candidate-length",
            Self::InvalidTolerance => "invalid-tolerance",
            Self::MeshDomain => "mesh-domain",
            Self::MeshCoordinates => "mesh-coordinates",
            Self::CandidateNonFinite => "candidate-non-finite",
            Self::CandidateBoundary => "candidate-boundary",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "u-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::Forcing,
            } => "f-non-finite",
            Self::PolynomialNonFinite {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "big-f-non-finite",
            Self::ExactSolutionBoundary => "exact-solution-boundary",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::Forcing,
            } => "derived-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            } => "derived-big-f-mismatch",
            Self::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ExactSolution,
            } => "derived-u-mismatch",
            Self::NonFiniteTightness => "non-finite-tightness",
            Self::InvalidEnclosure => "invalid-enclosure",
            Self::WorkPlanOverflow => "work-plan-overflow",
            Self::WorkPlanMismatch => "work-plan-mismatch",
        }
    }
}

impl EstimatorFamily {
    /// Stable id for ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            EstimatorFamily::EquilibratedFlux => "equilibrated-flux-1d",
            EstimatorFamily::Hierarchical => "hierarchical-h2",
        }
    }
}

/// The verifier's verdict on one candidate.
#[derive(Debug, Clone)]
pub struct VerifierReport {
    /// The certified error-bound enclosure (energy norm).
    pub bound: Iv,
    /// Accept ⟺ `bound.hi ≤ tolerance` for an admitted finite report.
    pub accept: bool,
    /// The verified color carried by an ACCEPT (`None` on reject or refusal).
    pub color: Option<Color>,
    /// The tolerance tested against (feeds the planner).
    pub tolerance: f64,
    /// Estimator family id.
    pub family: &'static str,
    /// FNV hash of the reconstructed flux (ledger identity).
    pub flux_hash: u64,
    /// Structured refusal (`None` only when a finite bound was produced).
    pub refusal: Option<VerifierRefusal>,
}

impl VerifierReport {
    /// The review-round-3 ledger row (structured, never stdout).
    #[must_use]
    pub fn to_row(&self, problem: &str, oracle_error: f64) -> String {
        let problem = json_escape(problem);
        let family = json_escape(self.family);
        let bound_lo = finite_scientific(self.bound.lo);
        let bound_hi = finite_scientific(self.bound.hi);
        let oracle = finite_scientific(oracle_error);
        let tolerance = finite_scientific(self.tolerance);
        let effectivity = if self.refusal.is_none()
            && oracle_error.is_finite()
            && oracle_error > 0.0
            && self.bound.hi.is_finite()
        {
            finite_fixed(self.bound.hi / oracle_error)
        } else if self.refusal.is_none() && oracle_error == 0.0 {
            "1.0000".to_string()
        } else {
            "null".to_string()
        };
        let refusal = self.refusal.map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", reason.id()),
        );
        let verdict = if self.refusal.is_some() {
            "refused"
        } else if self.accept {
            "accept"
        } else {
            "reject"
        };
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"problem\":\"{problem}\",\"estimator_family_id\":\"{}\",\
             \"flux_hash\":\"{:016X}\",\"bound_lo\":{bound_lo},\"bound_hi\":{bound_hi},\
             \"oracle_true_error\":{oracle},\"effectivity\":{effectivity},\
             \"verdict\":\"{verdict}\",\"tolerance\":{tolerance},\"refusal\":{refusal}}}",
            family, self.flux_hash,
        );
        s
    }
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

fn finite_scientific(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        "null".to_string()
    }
}

fn finite_fixed(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.4}")
    } else {
        "null".to_string()
    }
}

fn fnv_extend(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

enum VerifierRunError<E> {
    Callback(E),
    Refusal(VerifierRefusal),
}

struct VerifierDriver<F> {
    callback: F,
    plan: VerifierWorkPlan,
    completed_work_units: u128,
    phase: VerifierPhase,
}

impl<F> VerifierDriver<F> {
    fn new(plan: VerifierWorkPlan, callback: F) -> Self {
        Self {
            callback,
            plan,
            completed_work_units: 0,
            phase: VerifierPhase::Validation,
        }
    }

    fn emit<E>(&mut self, kind: VerifierCheckpointKind) -> Result<(), E>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        (self.callback)(VerifierProgress {
            kind,
            phase: self.phase,
            completed_work_units: self.completed_work_units,
            planned_work_units: self.plan.total,
        })
    }

    fn enter<E>(&mut self, phase: VerifierPhase) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.phase = phase;
        self.emit(VerifierCheckpointKind::PhaseEntry)
            .map_err(VerifierRunError::Callback)
    }

    fn complete_one<E>(&mut self) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.completed_work_units = self
            .completed_work_units
            .checked_add(1)
            .filter(|completed| *completed <= self.plan.total)
            .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
        if self
            .completed_work_units
            .is_multiple_of(VERIFIER_POLL_STRIDE_WORK_UNITS)
        {
            self.emit(VerifierCheckpointKind::WorkBoundary)
                .map_err(VerifierRunError::Callback)?;
        }
        Ok(())
    }

    fn refusal_flush<E>(&mut self) -> Result<(), E>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.emit(VerifierCheckpointKind::RefusalFlush)
    }

    fn require_completed<E>(&self, expected: u128) -> Result<(), VerifierRunError<E>> {
        if self.completed_work_units == expected {
            Ok(())
        } else {
            Err(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))
        }
    }

    fn publication<E>(&mut self) -> Result<(), VerifierRunError<E>>
    where
        F: FnMut(VerifierProgress) -> Result<(), E>,
    {
        self.emit(VerifierCheckpointKind::Publication)
            .map_err(VerifierRunError::Callback)
    }
}

#[allow(clippy::too_many_lines)] // One auditable refusal order with exact partial-work semantics.
fn validate_inputs_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<(crate::fem1d::Poly, crate::fem1d::Poly), VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        let valid_count =
            (1..=MAX_VERIFIER_POLY_COEFFICIENTS).contains(&polynomial.coefficients().len());
        driver.complete_one()?;
        if !valid_count {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::PolynomialCoefficientCount { polynomial: role },
            ));
        }
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(VerifierRunError::Refusal(VerifierRefusal::InvalidTolerance));
    }
    if problem.mesh().first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || problem.mesh().last().map(|value| value.to_bits()) != Some(1.0_f64.to_bits())
    {
        return Err(VerifierRunError::Refusal(VerifierRefusal::MeshDomain));
    }
    for value in problem.mesh() {
        let finite = value.is_finite();
        driver.complete_one()?;
        if !finite {
            return Err(VerifierRunError::Refusal(VerifierRefusal::MeshCoordinates));
        }
    }
    for pair in problem.mesh().windows(2) {
        let increasing = pair[0] < pair[1];
        driver.complete_one()?;
        if !increasing {
            return Err(VerifierRunError::Refusal(VerifierRefusal::MeshCoordinates));
        }
    }
    for value in candidate {
        let finite = value.is_finite();
        driver.complete_one()?;
        if !finite {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::CandidateNonFinite,
            ));
        }
    }
    if candidate.first().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
        || candidate.last().map(|value| value.to_bits()) != Some(0.0_f64.to_bits())
    {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::CandidateBoundary,
        ));
    }
    for (role, polynomial) in [
        (VerifierPolynomial::ExactSolution, problem.exact_solution()),
        (VerifierPolynomial::Forcing, problem.forcing()),
        (
            VerifierPolynomial::ForcingAntiderivative,
            problem.rounded_forcing_antiderivative(),
        ),
    ] {
        for value in polynomial.coefficients() {
            let finite = value.is_finite();
            driver.complete_one()?;
            if !finite {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::PolynomialNonFinite { polynomial: role },
                ));
            }
        }
    }
    // The fixed 34-limb boundary certifier is an atomic U<=6 micro-tile. Credit
    // its coefficient units only after the result exists, so cancellation can
    // conservatively under-report by at most one tile but never over-report.
    let first_is_zero = problem
        .exact_solution()
        .coefficients()
        .first()
        .map(|value| value.to_bits())
        == Some(0.0_f64.to_bits());
    let zero_at_one = problem.exact_solution().is_exactly_zero_at_one();
    let exact_solution_has_boundary = first_is_zero && zero_at_one;
    for _ in problem.exact_solution().coefficients() {
        driver.complete_one()?;
    }
    if !exact_solution_has_boundary {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::ExactSolutionBoundary,
        ));
    }

    let expected_f = problem
        .exact_solution()
        .derive()
        .and_then(|derivative| derivative.derive())
        .map(crate::fem1d::Poly::neg)
        .map_err(|_| {
            VerifierRunError::Refusal(VerifierRefusal::PolynomialNonFinite {
                polynomial: VerifierPolynomial::Forcing,
            })
        })?;
    for _ in problem.exact_solution().coefficients() {
        driver.complete_one()?;
    }
    if problem.forcing().coefficients().len() != expected_f.coefficients().len() {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::Forcing,
            },
        ));
    }
    for (declared, expected) in problem
        .forcing()
        .coefficients()
        .iter()
        .zip(expected_f.coefficients())
    {
        let equal = declared.to_bits() == expected.to_bits();
        driver.complete_one()?;
        if !equal {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::DerivedPolynomialMismatch {
                    polynomial: VerifierPolynomial::Forcing,
                },
            ));
        }
    }
    let expected_big_f = expected_f.antiderive().map_err(|_| {
        VerifierRunError::Refusal(VerifierRefusal::PolynomialNonFinite {
            polynomial: VerifierPolynomial::ForcingAntiderivative,
        })
    })?;
    for _ in problem.forcing().coefficients() {
        driver.complete_one()?;
    }
    if problem
        .rounded_forcing_antiderivative()
        .coefficients()
        .len()
        != expected_big_f.coefficients().len()
    {
        return Err(VerifierRunError::Refusal(
            VerifierRefusal::DerivedPolynomialMismatch {
                polynomial: VerifierPolynomial::ForcingAntiderivative,
            },
        ));
    }
    for (declared, expected) in problem
        .rounded_forcing_antiderivative()
        .coefficients()
        .iter()
        .zip(expected_big_f.coefficients())
    {
        let equal = declared.to_bits() == expected.to_bits();
        driver.complete_one()?;
        if !equal {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::DerivedPolynomialMismatch {
                    polynomial: VerifierPolynomial::ForcingAntiderivative,
                },
            ));
        }
    }
    Ok((expected_f, expected_big_f))
}

fn tightness_constant_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    big_f: &crate::fem1d::Poly,
    driver: &mut VerifierDriver<F>,
) -> Result<f64, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut mean = 0.0;
    for element in 0..problem.mesh().len() - 1 {
        let (x0, x1) = (problem.mesh()[element], problem.mesh()[element + 1]);
        let h = x1 - x0;
        let slope = (candidate[element + 1] - candidate[element]) / h;
        if !h.is_finite() || h <= 0.0 || !slope.is_finite() {
            return Err(VerifierRunError::Refusal(
                VerifierRefusal::NonFiniteTightness,
            ));
        }
        for (point, weight) in gauss5(x0, x1) {
            let value = big_f.eval(point) + slope;
            let contribution = weight * value;
            if !point.is_finite()
                || !weight.is_finite()
                || !value.is_finite()
                || !contribution.is_finite()
            {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::NonFiniteTightness,
                ));
            }
            mean += contribution;
            if !mean.is_finite() {
                return Err(VerifierRunError::Refusal(
                    VerifierRefusal::NonFiniteTightness,
                ));
            }
        }
        driver.complete_one()?;
    }
    Ok(mean)
}

fn finite_interval(interval: Iv) -> Result<Iv, VerifierRefusal> {
    if interval.lo.is_finite() && interval.hi.is_finite() && interval.lo <= interval.hi {
        Ok(interval)
    } else {
        Err(VerifierRefusal::InvalidEnclosure)
    }
}

fn interval_element_geometry(x0: f64, x1: f64) -> Result<(Iv, Iv, Iv), VerifierRefusal> {
    let x0 = Iv::point(x0);
    let x1 = Iv::point(x1);
    let h = finite_interval(x1.sub(x0))?;
    if h.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    let midpoint = finite_interval(x0.add(x1).scale_pos(0.5))?;
    let half = finite_interval(h.scale_pos(0.5))?;
    if half.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((h, midpoint, half))
}

fn interval_candidate_slope(first: f64, second: f64, h: Iv) -> Result<Iv, VerifierRefusal> {
    let difference = finite_interval(Iv::point(second).sub(Iv::point(first)))?;
    finite_interval(difference.div_pos(h))
}

fn interval_quadrature_geometry(
    midpoint: Iv,
    half: Iv,
    node_constant: f64,
    weight_constant: f64,
) -> Result<(Iv, Iv), VerifierRefusal> {
    let node = finite_interval(midpoint.add(half.mul(iv_c(node_constant))))?;
    let weight = finite_interval(half.mul(iv_c(weight_constant)))?;
    if weight.lo <= 0.0 {
        return Err(VerifierRefusal::InvalidEnclosure);
    }
    Ok((node, weight))
}

fn interval_antiderivative_coefficient(
    coefficient: f64,
    exponent: usize,
) -> Result<Iv, VerifierRefusal> {
    finite_interval(Iv::point(coefficient).div_pos(Iv::point(exponent as f64)))
}

fn interval_forcing_antiderivative(
    forcing: &crate::fem1d::Poly,
    x: Iv,
) -> Result<Iv, VerifierRefusal> {
    // F(x) = x * Horner(f_k / (k + 1)). Coefficient division is itself
    // intervalized: the rounded coefficients in `big_f` are replay metadata,
    // not point enclosures of the exact antiderivative of the authoritative f.
    let mut accumulated = Iv::zero();
    for (degree, coefficient) in forcing.coefficients().iter().copied().enumerate().rev() {
        let antiderivative_coefficient =
            interval_antiderivative_coefficient(coefficient, degree + 1)?;
        accumulated = finite_interval(accumulated.mul(x).add(antiderivative_coefficient))?;
    }
    finite_interval(x.mul(accumulated))
}

fn equilibrated_bound_with_checkpoint<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    forcing: &crate::fem1d::Poly,
    c_star: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<Iv, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut eta_sq = Iv::zero();
    for element in 0..problem.mesh().len() - 1 {
        let (h, midpoint, half) =
            interval_element_geometry(problem.mesh()[element], problem.mesh()[element + 1])
                .map_err(VerifierRunError::Refusal)?;
        let slope = interval_candidate_slope(candidate[element], candidate[element + 1], h)
            .map_err(VerifierRunError::Refusal)?;
        for (node_constant, weight_constant) in GAUSS5_REF {
            let (node, weight) =
                interval_quadrature_geometry(midpoint, half, node_constant, weight_constant)
                    .map_err(VerifierRunError::Refusal)?;
            let antiderivative = interval_forcing_antiderivative(forcing, node)
                .map_err(VerifierRunError::Refusal)?;
            let residual = finite_interval(Iv::point(c_star).sub(antiderivative).sub(slope))
                .map_err(VerifierRunError::Refusal)?;
            let contribution =
                finite_interval(weight.mul(residual.sq())).map_err(VerifierRunError::Refusal)?;
            eta_sq =
                finite_interval(eta_sq.add(contribution)).map_err(VerifierRunError::Refusal)?;
        }
        driver.complete_one()?;
    }
    let bound = finite_interval(eta_sq.sqrt()).map_err(VerifierRunError::Refusal)?;
    if bound.lo < 0.0 {
        Err(VerifierRunError::Refusal(VerifierRefusal::InvalidEnclosure))
    } else {
        Ok(bound)
    }
}

fn flux_hash_with_checkpoint<F, E>(
    c_star: f64,
    forcing: &crate::fem1d::Poly,
    antiderivative: &crate::fem1d::Poly,
    driver: &mut VerifierDriver<F>,
) -> Result<u64, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let mut hash = fnv_extend(0xcbf2_9ce4_8422_2325, &c_star.to_bits().to_le_bytes());
    driver.complete_one()?;
    for polynomial in [forcing, antiderivative] {
        let length = u64::try_from(polynomial.coefficients().len())
            .map_err(|_| VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
        hash = fnv_extend(hash, &length.to_le_bytes());
        driver.complete_one()?;
        for coefficient in polynomial.coefficients() {
            hash = fnv_extend(hash, &coefficient.to_bits().to_le_bytes());
            driver.complete_one()?;
        }
    }
    Ok(hash)
}

fn refused(tolerance: f64, reason: VerifierRefusal) -> VerifierReport {
    VerifierReport {
        bound: Iv {
            lo: f64::NEG_INFINITY,
            hi: f64::INFINITY,
        },
        accept: false,
        color: None,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: 0,
        refusal: Some(reason),
    }
}

fn run_verifier<F, E>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    driver: &mut VerifierDriver<F>,
) -> Result<VerifierReport, VerifierRunError<E>>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    driver.enter(VerifierPhase::Validation)?;
    let (canonical_f, canonical_big_f) =
        validate_inputs_with_checkpoint(problem, candidate, tolerance, driver)?;
    driver.require_completed(driver.plan.validation)?;

    driver.enter(VerifierPhase::Tightness)?;
    // Any finite c is sound. This rounded optimizer affects tightness only.
    let c_star = tightness_constant_with_checkpoint(problem, candidate, &canonical_big_f, driver)?;
    let after_tightness = driver
        .plan
        .validation
        .checked_add(driver.plan.tightness)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_tightness)?;

    driver.enter(VerifierPhase::Equilibrated)?;
    let bound =
        equilibrated_bound_with_checkpoint(problem, candidate, &canonical_f, c_star, driver)?;
    let after_equilibrated = after_tightness
        .checked_add(driver.plan.equilibrated)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_equilibrated)?;

    driver.enter(VerifierPhase::Hash)?;
    let flux_hash = flux_hash_with_checkpoint(c_star, &canonical_f, &canonical_big_f, driver)?;
    let after_hash = after_equilibrated
        .checked_add(driver.plan.hash)
        .ok_or(VerifierRunError::Refusal(VerifierRefusal::WorkPlanMismatch))?;
    driver.require_completed(after_hash)?;

    let accept = bound.hi <= tolerance;
    let color = if accept {
        Some(Color::Verified {
            lo: 0.0,
            hi: bound.hi,
        })
    } else {
        None
    };
    driver.enter(VerifierPhase::Finalization)?;
    let report = VerifierReport {
        bound,
        accept,
        color,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash,
        refusal: None,
    };
    driver.complete_one()?;
    driver.require_completed(driver.plan.total)?;
    driver.publication()?;
    Ok(report)
}

/// The equilibrated-flux VERIFIER with an explicit sparse progress callback.
///
/// The callback runs at every phase entry, each invocation-global multiple of
/// [`VERIFIER_POLL_STRIDE_WORK_UNITS`], every structured-refusal flush, and the
/// final publication gate. Callback failure wins over any pending scientific
/// refusal and no report is returned. Shape refusals happen during constant-time
/// preflight and invoke no callback.
///
/// # Errors
/// Returns the callback's error unchanged. Scientific input or arithmetic
/// failures remain fail-closed [`VerifierReport`] values with a structured
/// [`VerifierReport::refusal`].
pub fn verify_with_checkpoint<E, F>(
    problem: &MmsProblem,
    candidate: &[f64],
    tolerance: f64,
    callback: F,
) -> Result<VerifierReport, E>
where
    F: FnMut(VerifierProgress) -> Result<(), E>,
{
    let plan = match VerifierWorkPlan::for_inputs(problem, candidate) {
        Ok(plan) => plan,
        Err(reason) => return Ok(refused(tolerance, reason)),
    };
    let mut driver = VerifierDriver::new(plan, callback);
    match run_verifier(problem, candidate, tolerance, &mut driver) {
        Ok(report) => Ok(report),
        Err(VerifierRunError::Callback(error)) => Err(error),
        Err(VerifierRunError::Refusal(reason)) => {
            driver.refusal_flush()?;
            Ok(refused(tolerance, reason))
        }
    }
}

/// The equilibrated-flux VERIFIER: certify (or reject) a candidate's
/// nodal values against `tolerance`. The returned bound is a TRUE
/// upper bound on `‖(u − u_h)′‖` whenever the candidate satisfies the
/// boundary conditions; the enclosure is rigorous by outward rounding.
///
/// This convenience wrapper is bitwise equivalent to
/// [`verify_with_checkpoint`] with an infallible no-op callback.
#[must_use]
pub fn verify(problem: &MmsProblem, candidate: &[f64], tolerance: f64) -> VerifierReport {
    match verify_with_checkpoint(problem, candidate, tolerance, |_| {
        Ok::<(), core::convert::Infallible>(())
    }) {
        Ok(report) => report,
        Err(never) => match never {},
    }
}

const GAUSS5_REF: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_08),
    (-0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.0, 0.568_888_888_888_888_9),
    (0.538_469_310_105_683_1, 0.478_628_670_499_366_47),
    (0.906_179_845_938_664, 0.236_926_885_056_189_08),
];

/// One-ulp-widened constant (the tabulated Gauss data carries ~1 ulp
/// of transcription error; widening keeps enclosures honest).
fn iv_c(v: f64) -> Iv {
    Iv {
        lo: crate::interval::down(v),
        hi: crate::interval::up(v),
    }
}

/// The INDEPENDENT second family: hierarchical estimate from a
/// uniformly refined solve (`h/2`). Not guaranteed — the falsifier's
/// cross-check, never a color source.
///
/// # Errors
/// Returns [`Fem1dError`] for malformed inputs, refinement overflow/resource
/// excess, allocation failure, or a non-finite estimate.
pub fn hierarchical_estimate(problem: &MmsProblem, candidate: &[f64]) -> Result<f64, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let fine_nodes = problem
        .mesh()
        .len()
        .checked_mul(2)
        .and_then(|nodes| nodes.checked_sub(1))
        .ok_or(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: usize::MAX,
            limit: MAX_FEM1D_MESH_NODES,
        })?;
    if fine_nodes > MAX_FEM1D_MESH_NODES {
        return Err(Fem1dError::ResourceLimit {
            resource: "hierarchical mesh nodes",
            requested: fine_nodes,
            limit: MAX_FEM1D_MESH_NODES,
        });
    }
    let mut fine_mesh = Vec::new();
    fine_mesh
        .try_reserve_exact(fine_nodes)
        .map_err(|_| Fem1dError::AllocationFailed {
            stage: "hierarchical mesh",
            requested: fine_nodes,
        })?;
    for w in problem.mesh().windows(2) {
        fine_mesh.push(w[0]);
        fine_mesh.push(f64::midpoint(w[0], w[1]));
    }
    fine_mesh.push(problem.mesh()[problem.mesh().len() - 1]);
    let fine = problem.with_mesh(fine_mesh)?;
    let fine_u = crate::fem1d::solve_p1(&fine)?;
    // ‖u_{h/2}′ − u_h′‖ over the fine mesh.
    let mut acc = 0.0;
    for e in 0..fine.mesh().len() - 1 {
        let (x0, x1) = (fine.mesh()[e], fine.mesh()[e + 1]);
        let h = x1 - x0;
        let fine_slope = (fine_u[e + 1] - fine_u[e]) / h;
        // The coarse element containing this fine element.
        let coarse_e = e / 2;
        let ch = problem.mesh()[coarse_e + 1] - problem.mesh()[coarse_e];
        let coarse_slope = (candidate[coarse_e + 1] - candidate[coarse_e]) / ch;
        let d = fine_slope - coarse_slope;
        let updated = (h * d).mul_add(d, acc);
        if !fine_slope.is_finite()
            || !coarse_slope.is_finite()
            || !d.is_finite()
            || !updated.is_finite()
        {
            return Err(Fem1dError::NonFiniteIntermediate {
                stage: "hierarchical estimate",
                index: Some(e),
            });
        }
        acc = updated;
    }
    let estimate = acc.sqrt();
    if estimate.is_finite() {
        Ok(estimate)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "hierarchical estimate",
            index: None,
        })
    }
}

/// The nonlinear WARM-START fallback: the candidate is accepted only
/// as a starting point; the measured value is iteration savings and
/// the color is ESTIMATED, never verified (the honest R1 boundary).
#[derive(Debug, Clone)]
pub struct WarmStartReport {
    /// Newton iterations from a cold start (zero).
    pub cold_iterations: u32,
    /// Newton iterations from the candidate.
    pub warm_iterations: u32,
    /// The color of the claim (always `Estimated`).
    pub color: Color,
}

/// Measure warm-start savings on the toy nonlinear class.
///
/// # Errors
/// Returns [`Fem1dError`] when either run is malformed, unusable, or does not
/// converge within the admitted budget. Nonconvergence never becomes savings.
pub fn warm_start(
    problem: &MmsProblem,
    candidate: &[f64],
    max_iter: u32,
) -> Result<WarmStartReport, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let zero = try_zeroed("cold nonlinear start", problem.mesh().len())?;
    let cold = crate::fem1d::solve_nonlinear(problem, &zero, max_iter)?;
    require_converged(&cold, "cold nonlinear solve")?;
    let warm = crate::fem1d::solve_nonlinear(problem, candidate, max_iter)?;
    require_converged(&warm, "warm nonlinear solve")?;
    Ok(WarmStartReport {
        cold_iterations: cold.iterations,
        warm_iterations: warm.iterations,
        color: Color::Estimated {
            estimator: "warm-start-iteration-savings".to_string(),
            dispersion: f64::INFINITY,
        },
    })
}

/// Convenience for the batteries: effectivity of a report against the
/// oracle.
///
/// # Errors
/// Returns [`Fem1dError`] when the independent oracle or report bound is not a
/// usable finite value. Oracle failure is never mapped to effectivity `1.0`.
pub fn effectivity(
    problem: &MmsProblem,
    candidate: &[f64],
    report: &VerifierReport,
) -> Result<f64, Fem1dError> {
    if report.refusal.is_some() {
        return Err(Fem1dError::InvalidScalar {
            field: "verifier report",
            reason: "refused reports have no defined effectivity",
        });
    }
    let truth = true_energy_error(problem, candidate)?;
    if !report.bound.hi.is_finite() || report.bound.hi < 0.0 {
        return Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity report bound",
            index: None,
        });
    }
    if truth == 0.0 {
        return Err(Fem1dError::InvalidScalar {
            field: "oracle true error",
            reason: "effectivity is undefined for a zero denominator",
        });
    }
    let value = report.bound.hi / truth;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "effectivity ratio",
            index: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::dd::Dd;

    fn contains_dd(interval: Iv, exact: Dd) -> bool {
        !exact.lt(Dd::from_f64(interval.lo)) && !Dd::from_f64(interval.hi).lt(exact)
    }

    #[test]
    fn intervalized_element_inputs_cover_independent_dd_oracle() {
        // These decimal-looking f64 inputs make every legacy point computation
        // below round away a nonzero residual. The double-double oracle therefore
        // detects removal of intervalized mesh, slope, node, or weight arithmetic.
        let (x0, x1) = (0.1, 0.4);
        let (candidate0, candidate1) = (0.1, 0.2);
        let (node_constant, weight_constant) = GAUSS5_REF[0];
        let (dx0, dx1) = (Dd::from_f64(x0), Dd::from_f64(x1));
        let half_constant = Dd::from_f64(0.5);

        let (h, midpoint, half) = interval_element_geometry(x0, x1).unwrap();
        let exact_h = dx1 - dx0;
        let exact_midpoint = (dx0 + dx1) * half_constant;
        let exact_half = exact_h * half_constant;
        assert_ne!(exact_h, Dd::from_f64(x1 - x0));
        assert_ne!(exact_midpoint, Dd::from_f64(f64::midpoint(x0, x1)));
        assert_ne!(exact_half, Dd::from_f64((x1 - x0) * 0.5));
        assert!(contains_dd(h, exact_h));
        assert!(contains_dd(midpoint, exact_midpoint));
        assert!(contains_dd(half, exact_half));

        let slope = interval_candidate_slope(candidate0, candidate1, h).unwrap();
        let exact_difference = Dd::from_f64(candidate1) - Dd::from_f64(candidate0);
        let rounded_slope = (candidate1 - candidate0) / (x1 - x0);
        assert_ne!(Dd::from_f64(rounded_slope) * exact_h, exact_difference);
        assert!(!(exact_difference).lt(Dd::from_f64(slope.lo) * exact_h));
        assert!(!(Dd::from_f64(slope.hi) * exact_h).lt(exact_difference));

        let (node, weight) =
            interval_quadrature_geometry(midpoint, half, node_constant, weight_constant).unwrap();
        let exact_node = exact_midpoint + exact_half * Dd::from_f64(node_constant);
        let exact_weight = exact_half * Dd::from_f64(weight_constant);
        let rounded_node = f64::midpoint(x0, x1) + (x1 - x0) * 0.5 * node_constant;
        let rounded_weight = (x1 - x0) * 0.5 * weight_constant;
        assert_ne!(exact_node, Dd::from_f64(rounded_node));
        assert_ne!(exact_weight, Dd::from_f64(rounded_weight));
        assert!(contains_dd(node, exact_node));
        assert!(contains_dd(weight, exact_weight));

        // `1/3` is not representable. The coefficient interval must reach the
        // side of the rounded quotient selected by the exact FMA residual;
        // treating the rounded antiderivative coefficient as a point fails it.
        let coefficient = interval_antiderivative_coefficient(1.0, 3).unwrap();
        let rounded = 1.0_f64 / 3.0;
        let residual = rounded.mul_add(3.0, -1.0);
        if residual > 0.0 {
            assert!(coefficient.lo <= crate::interval::down(rounded));
        } else if residual < 0.0 {
            assert!(coefficient.hi >= crate::interval::up(rounded));
        } else {
            assert!(coefficient.lo <= rounded && rounded <= coefficient.hi);
        }
    }

    #[test]
    fn gauss_constants_enclose_independent_truth_brackets() {
        // Each bit pair is the adjacent-f64 bracket around the corresponding
        // high-precision Gauss-Legendre constant, derived independently from
        // the decimal reference values. Fifteen-digit literals miss some
        // weights by up to eight ulps, so this locks the certified quadrature
        // inputs rather than merely checking that `iv_c` widens its input.
        let positive_constants = [
            (
                GAUSS5_REF[4].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[3].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
            (
                GAUSS5_REF[0].1,
                0x3fce_539e_c36e_038c,
                0x3fce_539e_c36e_038d,
            ),
            (
                GAUSS5_REF[1].1,
                0x3fde_a1da_25ae_415a,
                0x3fde_a1da_25ae_415b,
            ),
            (
                GAUSS5_REF[2].1,
                0x3fe2_3456_789a_bcdf,
                0x3fe2_3456_789a_bce0,
            ),
        ];
        for (constant, lower_bits, upper_bits) in positive_constants {
            let interval = iv_c(constant);
            assert!(interval.lo <= f64::from_bits(lower_bits));
            assert!(interval.hi >= f64::from_bits(upper_bits));
        }

        for (constant, positive_lower_bits, positive_upper_bits) in [
            (
                GAUSS5_REF[0].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                GAUSS5_REF[1].0,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
        ] {
            let interval = iv_c(constant);
            assert!(interval.lo <= -f64::from_bits(positive_upper_bits));
            assert!(interval.hi >= -f64::from_bits(positive_lower_bits));
        }
    }
}
