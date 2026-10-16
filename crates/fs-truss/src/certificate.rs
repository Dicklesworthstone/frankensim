//! Outward-verified optimum bounds for the truss layout LP.
//!
//! PDHG iterates are diagnostics: a small residual does not make their
//! objective a primal upper bound.  This module repairs the signed member
//! forces through a deterministically selected square equilibrium basis and
//! proves that repair with an interval Neumann bound.  It independently
//! shrinks and rechecks the dual iterate.  A finite certificate is published
//! only after both witnesses pass and their canonical inputs are hashed.

use crate::lp::{LayoutLp, MAX_PDHG_ITERS, PdhgReport, PdhgSettings};
use fs_blake3::{Blake3, ContentHash, hash_domain};
use fs_exec::Cx;
use fs_ivl::Interval;
use fs_sparse::Csr;

const CERTIFICATE_POLL_STRIDE: usize = 1_024;
const CERTIFICATE_ID_DOMAIN: &str = "frankensim.fs-truss.layout-certificate.v1";
const PROBLEM_ID_DOMAIN: &str = "frankensim.fs-truss.layout-problem.v1";
const INPUT_ID_DOMAIN: &str = "frankensim.fs-truss.layout-certificate-input.v1";
const SOLVER_STATE_ID_DOMAIN: &str = "frankensim.fs-truss.pdhg-state.v1";

#[cfg(test)]
std::thread_local! {
    static IDENTITY_HASH_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Hard ceiling for active equilibrium rows in one certificate attempt.
pub const HARD_MAX_CERTIFICATE_ACTIVE_ROWS: usize = 4_096;
/// Hard ceiling for physical members in one certificate attempt.
pub const HARD_MAX_CERTIFICATE_MEMBERS: usize = 1_000_000;
/// Hard ceiling for retained dense scalar entries in the cold proof path.
pub const HARD_MAX_CERTIFICATE_DENSE_ENTRIES: usize = 16_777_216;
/// Hard ceiling for arithmetic work in one certificate attempt.
pub const HARD_MAX_CERTIFICATE_OPERATIONS: usize = 1_000_000_000;

/// Caller-selected limits beneath the certificate hard ceilings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LayoutCertificateLimits {
    active_rows: usize,
    members: usize,
    dense_entries: usize,
    operations: usize,
}

impl LayoutCertificateLimits {
    /// Admit a bounded certificate workload.
    ///
    /// # Errors
    /// Returns a structured error when any limit is zero or exceeds its hard
    /// ceiling.
    pub fn try_new(
        max_active_rows: usize,
        max_members: usize,
        max_dense_entries: usize,
        max_operations: usize,
    ) -> Result<Self, LayoutCertificateError> {
        validate_limit(
            "max_active_rows",
            max_active_rows,
            HARD_MAX_CERTIFICATE_ACTIVE_ROWS,
        )?;
        validate_limit("max_members", max_members, HARD_MAX_CERTIFICATE_MEMBERS)?;
        validate_limit(
            "max_dense_entries",
            max_dense_entries,
            HARD_MAX_CERTIFICATE_DENSE_ENTRIES,
        )?;
        validate_limit(
            "max_operations",
            max_operations,
            HARD_MAX_CERTIFICATE_OPERATIONS,
        )?;
        Ok(Self {
            active_rows: max_active_rows,
            members: max_members,
            dense_entries: max_dense_entries,
            operations: max_operations,
        })
    }

    /// Maximum nonzero equilibrium rows admitted by this attempt.
    #[must_use]
    pub const fn max_active_rows(self) -> usize {
        self.active_rows
    }

    /// Maximum physical members admitted by this attempt.
    #[must_use]
    pub const fn max_members(self) -> usize {
        self.members
    }

    /// Maximum dense scalar entries admitted by this attempt.
    #[must_use]
    pub const fn max_dense_entries(self) -> usize {
        self.dense_entries
    }

    /// Maximum arithmetic operations admitted by this attempt.
    #[must_use]
    pub const fn max_operations(self) -> usize {
        self.operations
    }
}

impl Default for LayoutCertificateLimits {
    fn default() -> Self {
        Self {
            active_rows: 128,
            members: 32_768,
            dense_entries: 2_097_152,
            operations: 50_000_000,
        }
    }
}

/// Borrowed authoritative LP arrays for consumers that already own canonical
/// sparse assembly (notably the browser campaign surface).
///
/// Certification uses `A`, `b`, and `c` only; a cached transpose is neither
/// accepted nor trusted. Construction checks dimensions, while the certificate
/// attempt performs the full paired-column and finite-domain validation.
#[derive(Debug, Clone, Copy)]
pub struct LayoutCertificateProblem<'a> {
    a: &'a Csr,
    c: &'a [f64],
    b: &'a [f64],
}

impl<'a> LayoutCertificateProblem<'a> {
    /// Borrow canonical equilibrium, objective, and load arrays.
    ///
    /// # Errors
    /// Rejects empty/odd split-variable layouts or mismatched dimensions.
    pub fn try_new(a: &'a Csr, c: &'a [f64], b: &'a [f64]) -> Result<Self, LayoutCertificateError> {
        if c.is_empty() || !c.len().is_multiple_of(2) || a.ncols() != c.len() {
            return Err(LayoutCertificateError::InvalidProblem {
                requirement: "must have a nonempty even [B | -B] variable layout",
            });
        }
        if a.nrows() != b.len() {
            return Err(LayoutCertificateError::InvalidProblem {
                requirement: "must have one load for every equilibrium row",
            });
        }
        Ok(Self { a, c, b })
    }

    /// Authoritative equilibrium matrix.
    #[must_use]
    pub const fn a(&self) -> &'a Csr {
        self.a
    }

    /// Symmetric split-variable costs.
    #[must_use]
    pub const fn c(&self) -> &'a [f64] {
        self.c
    }

    /// Equilibrium right-hand side.
    #[must_use]
    pub const fn b(&self) -> &'a [f64] {
        self.b
    }
}

trait CertificateProblemView {
    fn a(&self) -> &Csr;
    fn c(&self) -> &[f64];
    fn b(&self) -> &[f64];
}

impl CertificateProblemView for LayoutLp {
    fn a(&self) -> &Csr {
        LayoutLp::a(self)
    }

    fn c(&self) -> &[f64] {
        LayoutLp::c(self)
    }

    fn b(&self) -> &[f64] {
        LayoutLp::b(self)
    }
}

impl CertificateProblemView for LayoutCertificateProblem<'_> {
    fn a(&self) -> &Csr {
        self.a
    }

    fn c(&self) -> &[f64] {
        self.c
    }

    fn b(&self) -> &[f64] {
        self.b
    }
}

/// Stable primal-repair algorithm recorded in the certificate identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimalCorrectionMethod {
    /// Signed-force basis correction with a Neumann inverse enclosure.
    SignedForceBasisNeumannV1,
}

/// Stable arithmetic enclosure recorded in the certificate identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArithmeticEnclosureMethod {
    /// `fs-ivl::Interval` outward arithmetic.
    FsIvlIntervalV1,
}

/// A collision-resistant identity for a problem, input, or proof receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LayoutCertificateIdentity(ContentHash);

impl LayoutCertificateIdentity {
    /// Raw 32-byte identity.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal identity.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

/// A finite, outward-verified optimum interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedObjectiveBounds {
    lower: f64,
    upper: f64,
}

impl CertifiedObjectiveBounds {
    /// Outward lower endpoint, obtained from a verified dual witness.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Outward upper endpoint, obtained from an exactly feasible primal
    /// witness.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }
}

/// Structured hard errors: malformed caller state, invalid limits,
/// allocation failure, or cancellation.  Numerical inability to prove a
/// bound is represented by [`LayoutCertificateRefusal`] instead.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LayoutCertificateError {
    /// A limit is outside its admitted domain.
    InvalidLimit {
        /// Stable field name.
        field: &'static str,
        /// Stable requirement.
        requirement: &'static str,
    },
    /// An input vector has the wrong length.
    VectorLength {
        /// `x` or `y`.
        vector: &'static str,
        /// Required length.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// An input vector entry is outside its numerical domain.
    InvalidVector {
        /// `x` or `y`.
        vector: &'static str,
        /// Offending entry.
        index: usize,
        /// Stable requirement.
        requirement: &'static str,
    },
    /// The immutable LP violates the split-force structural contract.
    InvalidProblem {
        /// Stable requirement.
        requirement: &'static str,
    },
    /// A mutable report was not minted by this LP/settings/output state.
    ReportMismatch {
        /// Stable requirement.
        requirement: &'static str,
    },
    /// A fallible reservation failed.
    AllocationFailed {
        /// Stable resource name.
        resource: &'static str,
        /// Requested element count.
        requested: usize,
    },
    /// The enclosing execution scope was cancelled.
    Cancelled {
        /// Stable proof stage.
        stage: &'static str,
    },
}

impl core::fmt::Display for LayoutCertificateError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidLimit { field, requirement } => {
                write!(formatter, "certificate limit {field} {requirement}")
            }
            Self::VectorLength {
                vector,
                expected,
                actual,
            } => write!(
                formatter,
                "certificate vector {vector} length {actual}; expected {expected}"
            ),
            Self::InvalidVector {
                vector,
                index,
                requirement,
            } => write!(
                formatter,
                "certificate vector {vector}[{index}] {requirement}"
            ),
            Self::InvalidProblem { requirement } => {
                write!(formatter, "layout certificate problem {requirement}")
            }
            Self::ReportMismatch { requirement } => {
                write!(formatter, "layout certificate report {requirement}")
            }
            Self::AllocationFailed {
                resource,
                requested,
            } => write!(
                formatter,
                "layout certificate could not reserve {requested} elements for {resource}"
            ),
            Self::Cancelled { stage } => {
                write!(formatter, "layout certificate cancelled during {stage}")
            }
        }
    }
}

impl std::error::Error for LayoutCertificateError {}

/// Honest reasons why a well-formed attempt did not prove a finite enclosure.
#[derive(Debug, Clone, PartialEq)]
pub enum LayoutCertificateRefusal {
    /// A configured size cap was exceeded.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Required amount.
        required: usize,
        /// Admitted amount.
        limit: usize,
    },
    /// A structurally zero row has a nonzero load and is infeasible.
    InconsistentZeroRow {
        /// Original equilibrium row.
        row: usize,
    },
    /// The active equilibrium matrix has no full-row-rank member basis.
    RankDeficient {
        /// Number of nonzero equilibrium rows.
        active_rows: usize,
        /// Rank reached by deterministic elimination.
        rank: usize,
    },
    /// The approximate inverse could not be certified by `rho < 1`.
    IllConditioned {
        /// Outward upper bound on `||I - H M||_inf`.
        contraction_bound: f64,
    },
    /// Interval or floating arithmetic became non-finite.
    NonFiniteArithmetic {
        /// Stable proof stage.
        stage: &'static str,
    },
    /// The repaired witness did not survive the independent residual check.
    ResidualVerificationFailed {
        /// Original equilibrium row.
        row: usize,
    },
    /// Verified endpoints were not finite and ordered.
    InvalidObjectiveBounds {
        /// Candidate lower endpoint.
        lower: f64,
        /// Candidate upper endpoint.
        upper: f64,
    },
}

/// Result of a well-formed certificate attempt.
#[derive(Debug, Clone)]
// Boxing the successful certificate would add a non-reportable allocation to
// a path whose allocation failures are otherwise explicit structured errors.
#[allow(clippy::large_enum_variant)]
pub enum LayoutCertificateStatus {
    /// Both independently checked witnesses produced finite ordered bounds.
    Certified(LayoutOptimalityCertificate),
    /// No finite optimum interval was minted.
    Unavailable(LayoutCertificateRefusal),
}

/// Private-by-construction primal/dual witness and content-addressed receipt.
#[derive(Debug, Clone)]
pub struct LayoutOptimalityCertificate {
    bounds: CertifiedObjectiveBounds,
    repaired_member_forces: Vec<Interval>,
    repaired_split_forces: Vec<Interval>,
    equilibrium_residuals: Vec<Interval>,
    scaled_dual: Vec<f64>,
    dual_slacks: Vec<Interval>,
    dual_scale: f64,
    selected_members: Vec<usize>,
    contraction_bound: f64,
    problem_identity: LayoutCertificateIdentity,
    input_identity: LayoutCertificateIdentity,
    certificate_identity: LayoutCertificateIdentity,
    limits: LayoutCertificateLimits,
    correction_method: PrimalCorrectionMethod,
    enclosure_method: ArithmeticEnclosureMethod,
}

impl LayoutOptimalityCertificate {
    fn verifies_identity_checked(
        &self,
        cx: Option<&Cx<'_>>,
    ) -> Result<bool, LayoutCertificateError> {
        Ok(certificate_is_structurally_valid_checked(self, cx)?
            && certificate_identity(self, cx)? == self.certificate_identity)
    }

    fn verifies_for_problem_view_checked(
        &self,
        problem: &(impl CertificateProblemView + ?Sized),
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        cx: Option<&Cx<'_>>,
    ) -> Result<bool, LayoutCertificateError> {
        if let Some(cx) = cx {
            poll(cx, "certificate verification admission")?;
        }
        let Ok(members) = validate_input_shape(problem, x, y, settings) else {
            return Ok(false);
        };
        let Some(expected_variables) = self.repaired_member_forces.len().checked_mul(2) else {
            return Ok(false);
        };
        let expected_rows = self.equilibrium_residuals.len();
        if members != self.repaired_member_forces.len()
            || problem.a().nrows() != expected_rows
            || problem.a().ncols() != expected_variables
            || problem.b().len() != expected_rows
            || problem.c().len() != expected_variables
            || x.len() != expected_variables
            || y.len() != expected_rows
            || self.repaired_split_forces.len() != expected_variables
            || self.dual_slacks.len() != expected_variables
            || self.scaled_dual.len() != expected_rows
            || members > self.limits.max_members()
            || self.selected_members.len() > self.limits.max_active_rows()
        {
            return Ok(false);
        }
        let required_work = checked_work(
            self.selected_members.len(),
            members,
            problem.a().nrows(),
            problem.a().ncols(),
            problem.a().nnz(),
        )
        .unwrap_or(usize::MAX);
        if required_work > self.limits.max_operations() || !self.verifies_identity_checked(cx)? {
            return Ok(false);
        }
        let problem_identity = problem_identity(problem, cx)?;
        let input_identity = input_identity(
            self.problem_identity,
            x,
            y,
            settings,
            self.correction_method,
            self.enclosure_method,
            self.limits,
            cx,
        )?;
        Ok(self.problem_identity == problem_identity && self.input_identity == input_identity)
    }

    /// Verified optimum interval.
    #[must_use]
    pub const fn bounds(&self) -> CertifiedObjectiveBounds {
        self.bounds
    }

    /// Outward enclosures of the repaired signed member forces.
    #[must_use]
    pub fn repaired_member_forces(&self) -> &[Interval] {
        &self.repaired_member_forces
    }

    /// Outward nonnegative enclosures in `[q+ | q-]` order.
    #[must_use]
    pub fn repaired_split_forces(&self) -> &[Interval] {
        &self.repaired_split_forces
    }

    /// Independent outward residual sanity check for every original row.
    ///
    /// `0` containment alone is not an existence proof because these boxes
    /// forget cross-component correlation. Exact feasibility comes from the
    /// recorded Neumann basis method and its correlated signed-force repair.
    #[must_use]
    pub fn equilibrium_residuals(&self) -> &[Interval] {
        &self.equilibrium_residuals
    }

    /// Representable, uniformly scaled dual witness.
    #[must_use]
    pub fn scaled_dual(&self) -> &[f64] {
        &self.scaled_dual
    }

    /// Outward `c + A^T y` enclosures in split-variable order.
    #[must_use]
    pub fn dual_slacks(&self) -> &[Interval] {
        &self.dual_slacks
    }

    /// Uniform scale applied to the solver dual iterate.
    #[must_use]
    pub const fn dual_scale(&self) -> f64 {
        self.dual_scale
    }

    /// Deterministically selected physical member columns.
    #[must_use]
    pub fn selected_members(&self) -> &[usize] {
        &self.selected_members
    }

    /// Outward contraction bound proving the selected basis invertible.
    #[must_use]
    pub const fn contraction_bound(&self) -> f64 {
        self.contraction_bound
    }

    /// Identity of matrix, load, and costs.
    #[must_use]
    pub const fn problem_identity(&self) -> LayoutCertificateIdentity {
        self.problem_identity
    }

    /// Identity additionally binding settings and input iterates.
    #[must_use]
    pub const fn input_identity(&self) -> LayoutCertificateIdentity {
        self.input_identity
    }

    /// Identity of the complete retained witness.
    #[must_use]
    pub const fn certificate_identity(&self) -> LayoutCertificateIdentity {
        self.certificate_identity
    }

    /// Work and retained-state budget under which the proof was minted.
    #[must_use]
    pub const fn limits(&self) -> LayoutCertificateLimits {
        self.limits
    }

    /// Recorded correction algorithm.
    #[must_use]
    pub const fn correction_method(&self) -> PrimalCorrectionMethod {
        self.correction_method
    }

    /// Recorded outward arithmetic implementation.
    #[must_use]
    pub const fn enclosure_method(&self) -> ArithmeticEnclosureMethod {
        self.enclosure_method
    }

    /// Recompute the receipt hash and structural predicates with bounded
    /// cancellation polling.
    ///
    /// # Errors
    /// Returns a structured cancellation or allocation error while checking
    /// the retained receipt. A clean `false` means the receipt is malformed or
    /// its hash no longer matches.
    pub fn verifies_identity(&self, cx: &Cx<'_>) -> Result<bool, LayoutCertificateError> {
        self.verifies_identity_checked(Some(cx))
    }

    /// Confirm that this private proof object is bound to the supplied LP,
    /// settings, and input iterates as well as to its retained witness.
    ///
    /// # Errors
    /// Returns a structured cancellation or allocation error while hashing or
    /// checking the retained receipt. A clean `false` means it is not bound to
    /// this LP, settings, and state.
    pub fn verifies_for(
        &self,
        lp: &LayoutLp,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        cx: &Cx<'_>,
    ) -> Result<bool, LayoutCertificateError> {
        self.verifies_for_problem_view_checked(lp, x, y, settings, Some(cx))
    }

    /// Confirm that this private proof object is bound to borrowed canonical
    /// LP arrays, settings, and input iterates as well as to its retained
    /// witness.
    ///
    /// This is the transcribed-consumer sibling of [`Self::verifies_for`]. It
    /// is a cheap identity/structure preflight, not a replacement for the
    /// authoritative proof replay available on [`LayoutLp`].
    ///
    /// # Errors
    /// Returns a structured cancellation or allocation error while hashing or
    /// checking the retained receipt. A clean `false` means it is not bound to
    /// these arrays, settings, and state.
    pub fn verifies_for_problem(
        &self,
        problem: &LayoutCertificateProblem<'_>,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        cx: &Cx<'_>,
    ) -> Result<bool, LayoutCertificateError> {
        self.verifies_for_problem_view_checked(problem, x, y, settings, Some(cx))
    }
}

fn validate_limit(
    field: &'static str,
    requested: usize,
    hard_ceiling: usize,
) -> Result<(), LayoutCertificateError> {
    if requested == 0 || requested > hard_ceiling {
        return Err(LayoutCertificateError::InvalidLimit {
            field,
            requirement: "must be positive and no greater than its hard ceiling",
        });
    }
    Ok(())
}

fn poll(cx: &Cx<'_>, stage: &'static str) -> Result<(), LayoutCertificateError> {
    cx.checkpoint()
        .map_err(|_| LayoutCertificateError::Cancelled { stage })
}

fn empty_with_capacity<T>(
    requested: usize,
    resource: &'static str,
) -> Result<Vec<T>, LayoutCertificateError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(requested)
        .map_err(|_| LayoutCertificateError::AllocationFailed {
            resource,
            requested,
        })?;
    Ok(values)
}

fn zeroed_f64(
    requested: usize,
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<f64>, LayoutCertificateError> {
    let mut values = empty_with_capacity(requested, resource)?;
    for index in 0..requested {
        if index.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, stage)?;
        }
        values.push(0.0);
    }
    Ok(values)
}

fn point_intervals(
    requested: usize,
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<Interval>, LayoutCertificateError> {
    let mut values = empty_with_capacity(requested, resource)?;
    for index in 0..requested {
        if index.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, stage)?;
        }
        values.push(Interval::point(0.0));
    }
    Ok(values)
}

fn copied_f64(
    source: &[f64],
    resource: &'static str,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Vec<f64>, LayoutCertificateError> {
    let mut values = empty_with_capacity(source.len(), resource)?;
    for chunk in source.chunks(CERTIFICATE_POLL_STRIDE) {
        poll(cx, stage)?;
        values.extend_from_slice(chunk);
    }
    Ok(values)
}

fn checked_work(
    active_rows: usize,
    members: usize,
    total_rows: usize,
    variables: usize,
    nonzeros: usize,
) -> Option<usize> {
    let dense = active_rows.checked_mul(members)?;
    let square = active_rows.checked_mul(active_rows)?;
    let selection = square.checked_mul(members)?;
    let cubic = square.checked_mul(active_rows)?;
    // Sparse validation includes paired-column binary searches. Dual repair
    // performs at most 65 complete A scans (candidate plus 64 halvings), and
    // identities/residuals scan the same retained state again. These weights
    // deliberately overestimate those bounded loops.
    let sparse_and_hash = nonzeros.checked_mul(160)?;
    let vector_scans = total_rows
        .checked_mul(96)?
        .checked_add(variables.checked_mul(192)?)?;
    dense
        .checked_mul(8)?
        .checked_add(selection.checked_mul(2)?)?
        .checked_add(cubic.checked_mul(12)?)?
        .checked_add(sparse_and_hash)?
        .checked_add(vector_scans)
}

struct PollCounter {
    steps: usize,
}

impl PollCounter {
    fn new() -> Self {
        Self { steps: 0 }
    }

    fn step(&mut self, cx: &Cx<'_>, stage: &'static str) -> Result<(), LayoutCertificateError> {
        if self.steps.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, stage)?;
        }
        self.steps = self.steps.saturating_add(1);
        Ok(())
    }
}

struct IdentityWriter {
    hasher: Blake3,
    domain: &'static str,
}

impl IdentityWriter {
    fn new(domain: &'static str) -> Self {
        Self {
            hasher: Blake3::new(),
            domain,
        }
    }

    fn usize(&mut self, value: usize) {
        self.hasher.update(&(value as u64).to_le_bytes());
    }

    fn f64(&mut self, value: f64) {
        self.hasher.update(&value.to_bits().to_le_bytes());
    }

    fn bytes(&mut self, value: &[u8]) {
        self.usize(value.len());
        self.hasher.update(value);
    }

    fn identity(&mut self, value: LayoutCertificateIdentity) {
        self.hasher.update(value.as_bytes());
    }

    fn interval(&mut self, value: Interval) {
        self.f64(value.lo());
        self.f64(value.hi());
    }

    fn finish(self) -> LayoutCertificateIdentity {
        let stream_hash = self.hasher.finalize();
        LayoutCertificateIdentity(hash_domain(self.domain, stream_hash.as_bytes()))
    }
}

fn identity_poll(
    cx: Option<&Cx<'_>>,
    steps: usize,
    stage: &'static str,
) -> Result<(), LayoutCertificateError> {
    if steps.is_multiple_of(CERTIFICATE_POLL_STRIDE)
        && let Some(cx) = cx
    {
        poll(cx, stage)?;
    }
    Ok(())
}

fn problem_identity<P: CertificateProblemView + ?Sized>(
    lp: &P,
    cx: Option<&Cx<'_>>,
) -> Result<LayoutCertificateIdentity, LayoutCertificateError> {
    #[cfg(test)]
    IDENTITY_HASH_CALLS.with(|calls| calls.set(calls.get().saturating_add(1)));
    let mut writer = IdentityWriter::new(PROBLEM_ID_DOMAIN);
    let mut steps = 0usize;
    writer.usize(lp.a().nrows());
    writer.usize(lp.a().ncols());
    writer.usize(lp.a().nnz());
    for row in 0..lp.a().nrows() {
        identity_poll(cx, steps, "certificate problem-identity hash")?;
        steps = steps.saturating_add(1);
        let (columns, values) = lp.a().row(row);
        writer.usize(columns.len());
        for (&column, &value) in columns.iter().zip(values) {
            identity_poll(cx, steps, "certificate problem-identity hash")?;
            steps = steps.saturating_add(1);
            writer.usize(column);
            writer.f64(value);
        }
    }
    writer.usize(lp.b().len());
    for &value in lp.b() {
        identity_poll(cx, steps, "certificate problem-identity hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    writer.usize(lp.c().len());
    for &value in lp.c() {
        identity_poll(cx, steps, "certificate problem-identity hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    Ok(writer.finish())
}

#[allow(clippy::too_many_arguments)] // Every explicit field is identity-bearing.
fn input_identity(
    problem: LayoutCertificateIdentity,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    correction_method: PrimalCorrectionMethod,
    enclosure_method: ArithmeticEnclosureMethod,
    limits: LayoutCertificateLimits,
    cx: Option<&Cx<'_>>,
) -> Result<LayoutCertificateIdentity, LayoutCertificateError> {
    #[cfg(test)]
    IDENTITY_HASH_CALLS.with(|calls| calls.set(calls.get().saturating_add(1)));
    let mut writer = IdentityWriter::new(INPUT_ID_DOMAIN);
    let mut steps = 0usize;
    writer.identity(problem);
    writer.usize(settings.max_iters);
    writer.f64(settings.gap_tol);
    writer.usize(settings.check_every);
    writer.usize(limits.max_active_rows());
    writer.usize(limits.max_members());
    writer.usize(limits.max_dense_entries());
    writer.usize(limits.max_operations());
    writer.bytes(match correction_method {
        PrimalCorrectionMethod::SignedForceBasisNeumannV1 => b"signed-force-basis-neumann-v1",
    });
    writer.bytes(match enclosure_method {
        ArithmeticEnclosureMethod::FsIvlIntervalV1 => b"fs-ivl-interval-v1",
    });
    writer.bytes(fs_ivl::VERSION.as_bytes());
    writer.bytes(env!("CARGO_PKG_VERSION").as_bytes());
    writer.usize(x.len());
    for &value in x {
        identity_poll(cx, steps, "certificate input-identity hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    writer.usize(y.len());
    for &value in y {
        identity_poll(cx, steps, "certificate input-identity hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    Ok(writer.finish())
}

pub(crate) fn solver_state_identity(
    lp: &LayoutLp,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
) -> [u8; 32] {
    solver_state_identity_checked(lp, x, y, settings, None)
        .expect("infallible identity without checkpoint")
}

fn solver_state_identity_checked(
    lp: &LayoutLp,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    cx: Option<&Cx<'_>>,
) -> Result<[u8; 32], LayoutCertificateError> {
    let mut writer = IdentityWriter::new(SOLVER_STATE_ID_DOMAIN);
    let problem = problem_identity(lp, cx)?;
    let mut steps = 0usize;
    writer.identity(problem);
    writer.usize(settings.max_iters);
    writer.f64(settings.gap_tol);
    writer.usize(settings.check_every);
    writer.usize(x.len());
    for &value in x {
        identity_poll(cx, steps, "certificate solver-state hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    writer.usize(y.len());
    for &value in y {
        identity_poll(cx, steps, "certificate solver-state hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    Ok(*writer.finish().as_bytes())
}

fn certificate_identity(
    certificate: &LayoutOptimalityCertificate,
    cx: Option<&Cx<'_>>,
) -> Result<LayoutCertificateIdentity, LayoutCertificateError> {
    #[cfg(test)]
    IDENTITY_HASH_CALLS.with(|calls| calls.set(calls.get().saturating_add(1)));
    let mut writer = IdentityWriter::new(CERTIFICATE_ID_DOMAIN);
    let mut steps = 0usize;
    writer.identity(certificate.problem_identity);
    writer.identity(certificate.input_identity);
    writer.usize(certificate.limits.max_active_rows());
    writer.usize(certificate.limits.max_members());
    writer.usize(certificate.limits.max_dense_entries());
    writer.usize(certificate.limits.max_operations());
    writer.f64(certificate.bounds.lower);
    writer.f64(certificate.bounds.upper);
    writer.usize(certificate.repaired_member_forces.len());
    for &value in &certificate.repaired_member_forces {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.interval(value);
    }
    writer.usize(certificate.repaired_split_forces.len());
    for &value in &certificate.repaired_split_forces {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.interval(value);
    }
    writer.usize(certificate.equilibrium_residuals.len());
    for &value in &certificate.equilibrium_residuals {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.interval(value);
    }
    writer.usize(certificate.scaled_dual.len());
    for &value in &certificate.scaled_dual {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.f64(value);
    }
    writer.usize(certificate.dual_slacks.len());
    for &value in &certificate.dual_slacks {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.interval(value);
    }
    writer.f64(certificate.dual_scale);
    writer.usize(certificate.selected_members.len());
    for &member in &certificate.selected_members {
        identity_poll(cx, steps, "certificate witness-identity hash")?;
        steps = steps.saturating_add(1);
        writer.usize(member);
    }
    writer.f64(certificate.contraction_bound);
    writer.bytes(match certificate.correction_method {
        PrimalCorrectionMethod::SignedForceBasisNeumannV1 => b"signed-force-basis-neumann-v1",
    });
    writer.bytes(match certificate.enclosure_method {
        ArithmeticEnclosureMethod::FsIvlIntervalV1 => b"fs-ivl-interval-v1",
    });
    Ok(writer.finish())
}

fn certificate_is_structurally_valid_checked(
    certificate: &LayoutOptimalityCertificate,
    cx: Option<&Cx<'_>>,
) -> Result<bool, LayoutCertificateError> {
    let bounds = certificate.bounds;
    if !(bounds.lower.is_finite()
        && bounds.upper.is_finite()
        && bounds.lower <= bounds.upper
        && certificate.dual_scale.is_finite()
        && (0.0..=1.0).contains(&certificate.dual_scale)
        && certificate.contraction_bound.is_finite()
        && (0.0..1.0).contains(&certificate.contraction_bound)
        && certificate.repaired_split_forces.len() == 2 * certificate.repaired_member_forces.len()
        && certificate.dual_slacks.len() == certificate.repaired_split_forces.len()
        && certificate.scaled_dual.len() == certificate.equilibrium_residuals.len())
    {
        return Ok(false);
    }
    let mut steps = 0usize;
    for value in &certificate.repaired_member_forces {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if !value.lo().is_finite() || !value.hi().is_finite() {
            return Ok(false);
        }
    }
    for value in &certificate.repaired_split_forces {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if !value.lo().is_finite() || !value.hi().is_finite() || value.lo() < 0.0 {
            return Ok(false);
        }
    }
    for value in &certificate.equilibrium_residuals {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if !value.lo().is_finite() || !value.hi().is_finite() || !value.contains_zero() {
            return Ok(false);
        }
    }
    for value in &certificate.scaled_dual {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if !value.is_finite() {
            return Ok(false);
        }
    }
    for value in &certificate.dual_slacks {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if !value.lo().is_finite() || !value.hi().is_finite() || value.lo() < 0.0 {
            return Ok(false);
        }
    }
    for (index, &member) in certificate.selected_members.iter().enumerate() {
        identity_poll(cx, steps, "certificate structural verification")?;
        steps = steps.saturating_add(1);
        if member >= certificate.repaired_member_forces.len() {
            return Ok(false);
        }
        for &previous in &certificate.selected_members[..index] {
            identity_poll(cx, steps, "certificate structural verification")?;
            steps = steps.saturating_add(1);
            if previous == member {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn unavailable(reason: LayoutCertificateRefusal) -> LayoutCertificateStatus {
    LayoutCertificateStatus::Unavailable(reason)
}

enum ProofFailure {
    Error(LayoutCertificateError),
    Refusal(LayoutCertificateRefusal),
}

impl From<LayoutCertificateError> for ProofFailure {
    fn from(error: LayoutCertificateError) -> Self {
        Self::Error(error)
    }
}

type ProofResult<T> = Result<T, ProofFailure>;

fn refuse<T>(reason: LayoutCertificateRefusal) -> ProofResult<T> {
    Err(ProofFailure::Refusal(reason))
}

fn validate_input_shape<P: CertificateProblemView + ?Sized>(
    lp: &P,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
) -> Result<usize, LayoutCertificateError> {
    let variables = lp.c().len();
    if variables == 0 || !variables.is_multiple_of(2) || lp.a().ncols() != variables {
        return Err(LayoutCertificateError::InvalidProblem {
            requirement: "must have a nonempty even [B | -B] variable layout",
        });
    }
    if lp.a().nrows() != lp.b().len() {
        return Err(LayoutCertificateError::InvalidProblem {
            requirement: "must have one finite load for every equilibrium row",
        });
    }
    if x.len() != variables {
        return Err(LayoutCertificateError::VectorLength {
            vector: "x",
            expected: variables,
            actual: x.len(),
        });
    }
    if y.len() != lp.b().len() {
        return Err(LayoutCertificateError::VectorLength {
            vector: "y",
            expected: lp.b().len(),
            actual: y.len(),
        });
    }
    if settings.max_iters == 0 || settings.max_iters > MAX_PDHG_ITERS || settings.check_every == 0 {
        return Err(LayoutCertificateError::InvalidProblem {
            requirement: "must be bound to admitted nonzero solver controls",
        });
    }
    if !settings.gap_tol.is_finite() || !(0.0..=1.0).contains(&settings.gap_tol) {
        return Err(LayoutCertificateError::InvalidProblem {
            requirement: "must be bound to a finite gap tolerance in 0..=1",
        });
    }
    Ok(variables / 2)
}

fn validate_inputs<P: CertificateProblemView + ?Sized>(
    lp: &P,
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    cx: &Cx<'_>,
) -> Result<usize, LayoutCertificateError> {
    let members = validate_input_shape(lp, x, y, settings)?;
    for (index, value) in x.iter().enumerate() {
        if index.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate primal-input validation")?;
        }
        if !value.is_finite() || *value < 0.0 {
            return Err(LayoutCertificateError::InvalidVector {
                vector: "x",
                index,
                requirement: "must be finite and non-negative",
            });
        }
    }
    for (index, value) in y.iter().enumerate() {
        if index.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate dual-input validation")?;
        }
        if !value.is_finite() {
            return Err(LayoutCertificateError::InvalidVector {
                vector: "y",
                index,
                requirement: "must be finite",
            });
        }
    }
    for member in 0..members {
        if member.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate objective validation")?;
        }
        let positive = lp.c()[member];
        let negative = lp.c()[members + member];
        if !positive.is_finite() || positive <= 0.0 || positive.to_bits() != negative.to_bits() {
            return Err(LayoutCertificateError::InvalidProblem {
                requirement: "must have finite positive symmetric split-variable costs",
            });
        }
    }
    for (row, value) in lp.b().iter().enumerate() {
        if row.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate load validation")?;
        }
        if !value.is_finite() {
            return Err(LayoutCertificateError::InvalidProblem {
                requirement: "must have finite equilibrium loads",
            });
        }
    }
    Ok(members)
}

#[allow(clippy::float_cmp)] // Exact paired-column and structural-zero invariants.
fn collect_active_rows(
    lp: &(impl CertificateProblemView + ?Sized),
    members: usize,
    max_active_rows: usize,
    cx: &Cx<'_>,
) -> ProofResult<Vec<usize>> {
    let capacity = lp.a().nrows().min(max_active_rows.saturating_add(1));
    let mut active_rows = empty_with_capacity(capacity, "active equilibrium rows")?;
    let mut visits = 0usize;
    for row in 0..lp.a().nrows() {
        if row.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate split-matrix validation")?;
        }
        let (columns, values) = lp.a().row(row);
        let mut nonzero = false;
        for (&column, &value) in columns.iter().zip(values) {
            if visits.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
                poll(cx, "certificate split-matrix validation")?;
            }
            visits = visits.saturating_add(1);
            if column >= 2 * members || !value.is_finite() {
                return Err(LayoutCertificateError::InvalidProblem {
                    requirement: "must contain only finite in-range sparse entries",
                }
                .into());
            }
            let paired_column = if column < members {
                nonzero |= value != 0.0;
                column + members
            } else {
                column - members
            };
            let Ok(paired_index) = columns.binary_search(&paired_column) else {
                return Err(LayoutCertificateError::InvalidProblem {
                    requirement: "must retain exactly paired B and -B sparse columns",
                }
                .into());
            };
            if values[paired_index] != -value {
                return Err(LayoutCertificateError::InvalidProblem {
                    requirement: "must retain exact opposite B and -B coefficients",
                }
                .into());
            }
        }
        if nonzero {
            active_rows.push(row);
            if active_rows.len() > max_active_rows {
                return refuse(LayoutCertificateRefusal::ResourceLimit {
                    resource: "active equilibrium rows",
                    required: active_rows.len(),
                    limit: max_active_rows,
                });
            }
        } else if lp.b()[row] != 0.0 {
            return refuse(LayoutCertificateRefusal::InconsistentZeroRow { row });
        }
    }
    Ok(active_rows)
}

fn dense_first_block(
    lp: &(impl CertificateProblemView + ?Sized),
    active_rows: &[usize],
    members: usize,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, LayoutCertificateError> {
    let entries =
        active_rows
            .len()
            .checked_mul(members)
            .ok_or(LayoutCertificateError::InvalidProblem {
                requirement: "dense proof dimensions must not overflow",
            })?;
    let mut dense = zeroed_f64(
        entries,
        "dense equilibrium block",
        cx,
        "certificate dense-block allocation",
    )?;
    let mut visits = 0usize;
    for (dense_row, &sparse_row) in active_rows.iter().enumerate() {
        let (columns, values) = lp.a().row(sparse_row);
        for (&column, &value) in columns.iter().zip(values) {
            if visits.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
                poll(cx, "certificate dense-block fill")?;
            }
            visits = visits.saturating_add(1);
            if column < members {
                dense[dense_row * members + column] = value;
            }
        }
    }
    Ok(dense)
}

#[allow(clippy::float_cmp, clippy::too_many_lines)] // Complete-pivot rank pipeline; exact zero is its sole cutoff.
fn select_basis(
    dense: &[f64],
    rows: usize,
    members: usize,
    cx: &Cx<'_>,
) -> ProofResult<Vec<usize>> {
    if rows > members {
        return refuse(LayoutCertificateRefusal::RankDeficient {
            active_rows: rows,
            rank: members,
        });
    }
    let mut work = copied_f64(
        dense,
        "basis-selection matrix",
        cx,
        "certificate basis copy",
    )?;
    let mut column_order = empty_with_capacity(members, "basis column permutation")?;
    for column in 0..members {
        if column.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate column-permutation initialization")?;
        }
        column_order.push(column);
    }
    let mut row_order = empty_with_capacity(rows, "basis row permutation")?;
    for row in 0..rows {
        if row.is_multiple_of(CERTIFICATE_POLL_STRIDE) {
            poll(cx, "certificate row-permutation initialization")?;
        }
        row_order.push(row);
    }
    let mut counter = PollCounter::new();
    for pivot in 0..rows {
        let mut best_row = pivot;
        let mut best_column = pivot;
        let mut best_magnitude = 0.0f64;
        let mut best_key = (usize::MAX, usize::MAX);
        for row in pivot..rows {
            for column in pivot..members {
                counter.step(cx, "certificate basis search")?;
                let magnitude = work[row * members + column].abs();
                if !magnitude.is_finite() {
                    return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                        stage: "basis search",
                    });
                }
                let key = (column_order[column], row_order[row]);
                if magnitude > best_magnitude || (magnitude == best_magnitude && key < best_key) {
                    best_magnitude = magnitude;
                    best_row = row;
                    best_column = column;
                    best_key = key;
                }
            }
        }
        if best_magnitude == 0.0 {
            return refuse(LayoutCertificateRefusal::RankDeficient {
                active_rows: rows,
                rank: pivot,
            });
        }
        if best_row != pivot {
            for column in 0..members {
                counter.step(cx, "certificate basis row swap")?;
                work.swap(pivot * members + column, best_row * members + column);
            }
            row_order.swap(pivot, best_row);
        }
        if best_column != pivot {
            for row in 0..rows {
                counter.step(cx, "certificate basis column swap")?;
                work.swap(row * members + pivot, row * members + best_column);
            }
            column_order.swap(pivot, best_column);
        }
        let pivot_value = work[pivot * members + pivot];
        for row in (pivot + 1)..rows {
            counter.step(cx, "certificate basis elimination")?;
            let factor = work[row * members + pivot] / pivot_value;
            if !factor.is_finite() {
                return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                    stage: "basis elimination",
                });
            }
            work[row * members + pivot] = 0.0;
            for column in (pivot + 1)..members {
                counter.step(cx, "certificate basis elimination")?;
                let updated =
                    work[row * members + column] - factor * work[pivot * members + column];
                if !updated.is_finite() {
                    return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                        stage: "basis elimination",
                    });
                }
                work[row * members + column] = updated;
            }
        }
    }
    let mut selected = empty_with_capacity(rows, "selected basis members")?;
    for &member in &column_order[..rows] {
        counter.step(cx, "certificate basis publication")?;
        selected.push(member);
    }
    Ok(selected)
}

fn selected_matrix(
    dense: &[f64],
    rows: usize,
    members: usize,
    selected: &[usize],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, LayoutCertificateError> {
    let square = rows
        .checked_mul(rows)
        .ok_or(LayoutCertificateError::InvalidProblem {
            requirement: "selected basis dimensions must not overflow",
        })?;
    let mut matrix = zeroed_f64(
        square,
        "selected basis matrix",
        cx,
        "certificate basis-matrix allocation",
    )?;
    let mut counter = PollCounter::new();
    for row in 0..rows {
        for (column, &member) in selected.iter().enumerate() {
            counter.step(cx, "certificate basis-matrix fill")?;
            matrix[row * rows + column] = dense[row * members + member];
        }
    }
    Ok(matrix)
}

#[allow(clippy::float_cmp)] // Exact-zero pivot is fail-closed; intervals prove authority.
fn approximate_inverse(matrix: &[f64], size: usize, cx: &Cx<'_>) -> ProofResult<Vec<f64>> {
    let mut work = copied_f64(
        matrix,
        "inverse work matrix",
        cx,
        "certificate inverse copy",
    )?;
    let mut inverse = zeroed_f64(
        size * size,
        "approximate inverse",
        cx,
        "certificate inverse allocation",
    )?;
    let mut counter = PollCounter::new();
    for diagonal in 0..size {
        counter.step(cx, "certificate inverse initialization")?;
        inverse[diagonal * size + diagonal] = 1.0;
    }
    for pivot in 0..size {
        let mut best_row = pivot;
        let mut best = work[pivot * size + pivot].abs();
        for row in (pivot + 1)..size {
            counter.step(cx, "certificate inverse pivot search")?;
            let candidate = work[row * size + pivot].abs();
            if candidate > best {
                best = candidate;
                best_row = row;
            }
        }
        if best == 0.0 || !best.is_finite() {
            return refuse(LayoutCertificateRefusal::RankDeficient {
                active_rows: size,
                rank: pivot,
            });
        }
        if best_row != pivot {
            for column in 0..size {
                counter.step(cx, "certificate inverse row swap")?;
                work.swap(pivot * size + column, best_row * size + column);
                inverse.swap(pivot * size + column, best_row * size + column);
            }
        }
        let pivot_value = work[pivot * size + pivot];
        for column in 0..size {
            counter.step(cx, "certificate inverse normalization")?;
            work[pivot * size + column] /= pivot_value;
            inverse[pivot * size + column] /= pivot_value;
            if !work[pivot * size + column].is_finite()
                || !inverse[pivot * size + column].is_finite()
            {
                return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                    stage: "basis inversion",
                });
            }
        }
        for row in 0..size {
            if row == pivot {
                continue;
            }
            let factor = work[row * size + pivot];
            for column in 0..size {
                counter.step(cx, "certificate inverse elimination")?;
                work[row * size + column] -= factor * work[pivot * size + column];
                inverse[row * size + column] -= factor * inverse[pivot * size + column];
                if !work[row * size + column].is_finite()
                    || !inverse[row * size + column].is_finite()
                {
                    return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                        stage: "basis inversion",
                    });
                }
            }
        }
    }
    Ok(inverse)
}

struct PrimalWitness {
    member_forces: Vec<Interval>,
    split_forces: Vec<Interval>,
    residuals: Vec<Interval>,
    selected_members: Vec<usize>,
    contraction_bound: f64,
    upper_bound: f64,
}

#[allow(clippy::too_many_arguments)] // One checked matrix entry plus polling state.
fn interval_matrix_product_entry(
    left: &[f64],
    right: &[f64],
    size: usize,
    row: usize,
    column: usize,
    counter: &mut PollCounter,
    cx: &Cx<'_>,
    stage: &'static str,
) -> Result<Interval, LayoutCertificateError> {
    let mut sum = Interval::point(0.0);
    for inner in 0..size {
        counter.step(cx, stage)?;
        sum = sum
            + Interval::point(left[row * size + inner])
                * Interval::point(right[inner * size + column]);
    }
    Ok(sum)
}

#[allow(clippy::too_many_lines)] // One auditable Neumann existence proof pipeline.
fn prove_primal(
    lp: &(impl CertificateProblemView + ?Sized),
    x: &[f64],
    members: usize,
    active_rows: &[usize],
    dense: &[f64],
    cx: &Cx<'_>,
) -> ProofResult<PrimalWitness> {
    let rows = active_rows.len();
    let selected = select_basis(dense, rows, members, cx)?;
    let basis = selected_matrix(dense, rows, members, &selected, cx)?;
    let inverse = approximate_inverse(&basis, rows, cx)?;

    let mut counter = PollCounter::new();
    let mut contraction_bound = 0.0f64;
    for row in 0..rows {
        let mut row_sum = Interval::point(0.0);
        for column in 0..rows {
            let product = interval_matrix_product_entry(
                &inverse,
                &basis,
                rows,
                row,
                column,
                &mut counter,
                cx,
                "certificate contraction product",
            )?;
            let identity = if row == column { 1.0 } else { 0.0 };
            let entry = Interval::point(identity) - product;
            row_sum = row_sum + entry.abs();
        }
        contraction_bound = contraction_bound.max(row_sum.hi());
    }
    if !contraction_bound.is_finite() {
        return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
            stage: "basis contraction bound",
        });
    }
    if contraction_bound >= 1.0 {
        return refuse(LayoutCertificateRefusal::IllConditioned { contraction_bound });
    }

    let mut signed_center = zeroed_f64(
        members,
        "signed iterate",
        cx,
        "certificate signed-iterate allocation",
    )?;
    for member in 0..members {
        counter.step(cx, "certificate signed-iterate conversion")?;
        let value = x[member] - x[members + member];
        if !value.is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "signed iterate conversion",
            });
        }
        signed_center[member] = value;
    }

    let mut rhs_residual = point_intervals(
        rows,
        "basis right-hand residual",
        cx,
        "certificate residual allocation",
    )?;
    for (dense_row, &original_row) in active_rows.iter().enumerate() {
        let mut residual = Interval::point(lp.b()[original_row]);
        for member in 0..members {
            counter.step(cx, "certificate right-hand residual")?;
            residual = residual
                - Interval::point(dense[dense_row * members + member])
                    * Interval::point(signed_center[member]);
        }
        if !residual.lo().is_finite() || !residual.hi().is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "right-hand residual",
            });
        }
        rhs_residual[dense_row] = residual;
    }

    // Any finite center is admissible.  Solving the midpoint reduces the
    // interval remainder before the Neumann existence bound is applied.
    let mut correction_center = zeroed_f64(
        rows,
        "basis correction center",
        cx,
        "certificate correction-center allocation",
    )?;
    for row in 0..rows {
        let mut value = 0.0f64;
        for column in 0..rows {
            counter.step(cx, "certificate correction-center solve")?;
            value = inverse[row * rows + column].mul_add(rhs_residual[column].midpoint(), value);
        }
        if !value.is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "correction center",
            });
        }
        correction_center[row] = value;
    }

    let mut centered_residual = point_intervals(
        rows,
        "centered basis residual",
        cx,
        "certificate centered-residual allocation",
    )?;
    for row in 0..rows {
        let mut value = rhs_residual[row];
        for column in 0..rows {
            counter.step(cx, "certificate centered residual")?;
            value = value
                - Interval::point(basis[row * rows + column])
                    * Interval::point(correction_center[column]);
        }
        centered_residual[row] = value;
    }

    let mut residual_norm = 0.0f64;
    for row in 0..rows {
        let mut value = Interval::point(0.0);
        for column in 0..rows {
            counter.step(cx, "certificate preconditioned residual")?;
            value =
                value + Interval::point(inverse[row * rows + column]) * centered_residual[column];
        }
        if !value.lo().is_finite() || !value.hi().is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "preconditioned residual",
            });
        }
        residual_norm = residual_norm.max(value.abs().hi());
    }
    let denominator = Interval::point(1.0) - Interval::point(contraction_bound);
    if denominator.lo() <= 0.0 || !residual_norm.is_finite() {
        return refuse(LayoutCertificateRefusal::IllConditioned { contraction_bound });
    }
    let radius_interval = Interval::point(residual_norm) / denominator;
    let radius = radius_interval.hi();
    if !radius.is_finite() || radius < 0.0 {
        return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
            stage: "correction radius",
        });
    }
    let correction_box = Interval::new(-radius, radius);

    let mut member_forces = point_intervals(
        members,
        "repaired member-force enclosures",
        cx,
        "certificate member-force allocation",
    )?;
    for member in 0..members {
        counter.step(cx, "certificate member-force initialization")?;
        member_forces[member] = Interval::point(signed_center[member]);
    }
    for (basis_column, &member) in selected.iter().enumerate() {
        counter.step(cx, "certificate repaired-force enclosure")?;
        member_forces[member] = member_forces[member]
            + Interval::point(correction_center[basis_column])
            + correction_box;
        if !member_forces[member].lo().is_finite() || !member_forces[member].hi().is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "repaired member-force enclosure",
            });
        }
    }

    let mut split_forces = point_intervals(
        2 * members,
        "repaired split-force enclosures",
        cx,
        "certificate split-force allocation",
    )?;
    for (member, &force) in member_forces.iter().enumerate() {
        counter.step(cx, "certificate nonnegative split enclosure")?;
        split_forces[member] = Interval::new(force.lo().max(0.0), force.hi().max(0.0));
        split_forces[members + member] =
            Interval::new((-force.hi()).max(0.0), (-force.lo()).max(0.0));
    }

    let mut residuals = point_intervals(
        lp.a().nrows(),
        "equilibrium residual enclosures",
        cx,
        "certificate residual-check allocation",
    )?;
    for (row, residual_slot) in residuals.iter_mut().enumerate() {
        counter.step(cx, "certificate independent equilibrium row")?;
        let (columns, values) = lp.a().row(row);
        let mut residual = -Interval::point(lp.b()[row]);
        for (&column, &value) in columns.iter().zip(values) {
            counter.step(cx, "certificate independent equilibrium check")?;
            residual = residual + Interval::point(value) * split_forces[column];
        }
        if !residual.lo().is_finite() || !residual.hi().is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "equilibrium residual verification",
            });
        }
        if !residual.contains_zero() {
            return refuse(LayoutCertificateRefusal::ResidualVerificationFailed { row });
        }
        *residual_slot = residual;
    }

    let mut objective = Interval::point(0.0);
    for (member, &force) in member_forces.iter().enumerate() {
        counter.step(cx, "certificate primal objective enclosure")?;
        objective = objective + Interval::point(lp.c()[member]) * force.abs();
    }
    let upper_bound = objective.hi();
    if !upper_bound.is_finite() {
        return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
            stage: "primal objective enclosure",
        });
    }

    Ok(PrimalWitness {
        member_forces,
        split_forces,
        residuals,
        selected_members: selected,
        contraction_bound,
        upper_bound,
    })
}

struct DualWitness {
    scaled_dual: Vec<f64>,
    slacks: Vec<Interval>,
    scale: f64,
    lower_bound: f64,
}

fn dual_slacks(
    lp: &(impl CertificateProblemView + ?Sized),
    scaled_dual: &[f64],
    cx: &Cx<'_>,
) -> Result<Vec<Interval>, LayoutCertificateError> {
    if scaled_dual.len() != lp.a().nrows() {
        return Err(LayoutCertificateError::VectorLength {
            vector: "y",
            expected: lp.a().nrows(),
            actual: scaled_dual.len(),
        });
    }
    let mut slacks = point_intervals(
        lp.c().len(),
        "dual slack enclosures",
        cx,
        "certificate dual-slack allocation",
    )?;
    let mut counter = PollCounter::new();
    for (slack, &cost) in slacks.iter_mut().zip(lp.c()) {
        counter.step(cx, "certificate dual-slack initialization")?;
        *slack = Interval::point(cost);
    }
    // Derive A^T y from the authoritative matrix instead of trusting its
    // cached transpose. The proof identity binds these exact A entries.
    for (row, &dual) in scaled_dual.iter().enumerate() {
        counter.step(cx, "certificate dual-slack row")?;
        let (columns, values) = lp.a().row(row);
        for (&variable, &value) in columns.iter().zip(values) {
            counter.step(cx, "certificate dual-slack evaluation")?;
            slacks[variable] = slacks[variable] + Interval::point(value) * Interval::point(dual);
        }
    }
    Ok(slacks)
}

fn prove_dual(
    lp: &(impl CertificateProblemView + ?Sized),
    y: &[f64],
    cx: &Cx<'_>,
) -> ProofResult<DualWitness> {
    let original_slacks = dual_slacks(lp, y, cx)?;
    let mut scale = 1.0f64;
    let mut counter = PollCounter::new();
    for (variable, slack) in original_slacks.iter().enumerate() {
        counter.step(cx, "certificate dual-scale derivation")?;
        // `slack - c` encloses A^T y.  Bounding its absolute value gives one
        // conservative scale valid for both members of each paired column.
        let activity = *slack - Interval::point(lp.c()[variable]);
        let magnitude = activity.abs().hi();
        if !magnitude.is_finite() {
            return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
                stage: "dual scale derivation",
            });
        }
        if magnitude > 0.0 {
            let ratio = Interval::point(lp.c()[variable]) / Interval::point(magnitude);
            scale = scale.min(ratio.lo().max(0.0));
        }
    }
    if scale < 1.0 {
        scale *= 0.5;
    }
    let mut scaled_dual = zeroed_f64(
        y.len(),
        "scaled dual witness",
        cx,
        "certificate scaled-dual allocation",
    )?;
    let mut accepted_slacks = None;
    for _ in 0..64 {
        let mut finite = true;
        for (destination, &source) in scaled_dual.iter_mut().zip(y) {
            counter.step(cx, "certificate scaled-dual construction")?;
            *destination = scale * source;
            finite &= destination.is_finite();
        }
        if !finite {
            scale *= 0.5;
            continue;
        }
        let candidate = dual_slacks(lp, &scaled_dual, cx)?;
        let mut feasible = true;
        for slack in &candidate {
            counter.step(cx, "certificate dual-feasibility publication check")?;
            feasible &= slack.lo().is_finite() && slack.hi().is_finite() && slack.lo() >= 0.0;
        }
        if feasible {
            accepted_slacks = Some(candidate);
            break;
        }
        scale *= 0.5;
    }
    let slacks = if let Some(slacks) = accepted_slacks {
        slacks
    } else {
        // The exact zero dual is always feasible because costs are positive.
        scale = 0.0;
        for value in &mut scaled_dual {
            counter.step(cx, "certificate zero-dual construction")?;
            *value = 0.0;
        }
        let mut slacks = point_intervals(
            lp.c().len(),
            "zero-dual slack enclosures",
            cx,
            "certificate zero-dual allocation",
        )?;
        for (slack, &cost) in slacks.iter_mut().zip(lp.c()) {
            counter.step(cx, "certificate zero-dual slack construction")?;
            *slack = Interval::point(cost);
        }
        slacks
    };

    let mut objective = Interval::point(0.0);
    for (&load, &dual) in lp.b().iter().zip(&scaled_dual) {
        counter.step(cx, "certificate dual objective enclosure")?;
        objective = objective - Interval::point(load) * Interval::point(dual);
    }
    let lower_bound = objective.lo().max(0.0);
    if !lower_bound.is_finite() {
        return refuse(LayoutCertificateRefusal::NonFiniteArithmetic {
            stage: "dual objective enclosure",
        });
    }
    Ok(DualWitness {
        scaled_dual,
        slacks,
        scale,
        lower_bound,
    })
}

#[allow(clippy::too_many_lines)] // Admission, both witnesses, identity, and atomic publication.
fn build_certificate(
    lp: &(impl CertificateProblemView + ?Sized),
    x: &[f64],
    y: &[f64],
    settings: PdhgSettings,
    limits: LayoutCertificateLimits,
    cx: &Cx<'_>,
) -> ProofResult<LayoutOptimalityCertificate> {
    poll(cx, "certificate admission")?;
    let members = validate_input_shape(lp, x, y, settings)?;
    if members > limits.max_members() {
        return refuse(LayoutCertificateRefusal::ResourceLimit {
            resource: "physical members",
            required: members,
            limit: limits.max_members(),
        });
    }
    let admission_work = checked_work(0, members, lp.a().nrows(), lp.a().ncols(), lp.a().nnz())
        .unwrap_or(usize::MAX);
    if admission_work > limits.max_operations() {
        return refuse(LayoutCertificateRefusal::ResourceLimit {
            resource: "certificate arithmetic operations",
            required: admission_work,
            limit: limits.max_operations(),
        });
    }
    let _ = validate_inputs(lp, x, y, settings, cx)?;
    let active_rows = collect_active_rows(lp, members, limits.max_active_rows(), cx)?;
    let rows = active_rows.len();
    let dense = rows.saturating_mul(members);
    let square = rows.saturating_mul(rows);
    let retained_dense = dense
        .checked_mul(2)
        .and_then(|value| value.checked_add(square.checked_mul(6)?))
        .and_then(|value| value.checked_add(members.checked_mul(16)?))
        .and_then(|value| value.checked_add(lp.a().nrows().checked_mul(4)?))
        .unwrap_or(usize::MAX);
    if retained_dense > limits.max_dense_entries() {
        return refuse(LayoutCertificateRefusal::ResourceLimit {
            resource: "dense proof entries",
            required: retained_dense,
            limit: limits.max_dense_entries(),
        });
    }
    let required_work = checked_work(rows, members, lp.a().nrows(), lp.a().ncols(), lp.a().nnz())
        .unwrap_or(usize::MAX);
    if required_work > limits.max_operations() {
        return refuse(LayoutCertificateRefusal::ResourceLimit {
            resource: "certificate arithmetic operations",
            required: required_work,
            limit: limits.max_operations(),
        });
    }
    let dense_block = dense_first_block(lp, &active_rows, members, cx)?;
    let primal = prove_primal(lp, x, members, &active_rows, &dense_block, cx)?;
    let dual = prove_dual(lp, y, cx)?;
    let bounds = CertifiedObjectiveBounds {
        lower: dual.lower_bound,
        upper: primal.upper_bound,
    };
    if !bounds.lower.is_finite() || !bounds.upper.is_finite() || bounds.lower > bounds.upper {
        return refuse(LayoutCertificateRefusal::InvalidObjectiveBounds {
            lower: bounds.lower,
            upper: bounds.upper,
        });
    }

    poll(cx, "certificate identity preflight")?;
    let correction_method = PrimalCorrectionMethod::SignedForceBasisNeumannV1;
    let enclosure_method = ArithmeticEnclosureMethod::FsIvlIntervalV1;
    let problem_identity = problem_identity(lp, Some(cx))?;
    let input_identity = input_identity(
        problem_identity,
        x,
        y,
        settings,
        correction_method,
        enclosure_method,
        limits,
        Some(cx),
    )?;
    let mut certificate = LayoutOptimalityCertificate {
        bounds,
        repaired_member_forces: primal.member_forces,
        repaired_split_forces: primal.split_forces,
        equilibrium_residuals: primal.residuals,
        scaled_dual: dual.scaled_dual,
        dual_slacks: dual.slacks,
        dual_scale: dual.scale,
        selected_members: primal.selected_members,
        contraction_bound: primal.contraction_bound,
        problem_identity,
        input_identity,
        certificate_identity: LayoutCertificateIdentity(ContentHash([0; 32])),
        limits,
        correction_method,
        enclosure_method,
    };
    certificate.certificate_identity = certificate_identity(&certificate, Some(cx))?;
    if !certificate_is_structurally_valid_checked(&certificate, Some(cx))? {
        return refuse(LayoutCertificateRefusal::InvalidObjectiveBounds {
            lower: bounds.lower,
            upper: bounds.upper,
        });
    }
    poll(cx, "certificate publication")?;
    Ok(certificate)
}

impl LayoutCertificateProblem<'_> {
    /// Repair and outward-verify witnesses over these borrowed canonical LP
    /// arrays. This is the same proof kernel used by [`LayoutLp`].
    ///
    /// # Errors
    /// Returns a hard error for malformed state, invalid proof limits,
    /// allocation failure, or cancellation. Sound numerical inability is an
    /// `Unavailable` status.
    pub fn certify_optimum(
        &self,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        limits: LayoutCertificateLimits,
        cx: &Cx<'_>,
    ) -> Result<LayoutCertificateStatus, LayoutCertificateError> {
        match build_certificate(self, x, y, settings, limits, cx) {
            Ok(certificate) => Ok(LayoutCertificateStatus::Certified(certificate)),
            Err(ProofFailure::Refusal(reason)) => Ok(unavailable(reason)),
            Err(ProofFailure::Error(error)) => Err(error),
        }
    }
}

impl LayoutLp {
    /// Repair and outward-verify primal and dual witnesses for this layout LP.
    ///
    /// A small PDHG residual alone is never promoted.
    ///
    /// # Errors
    /// Returns a hard error for malformed state, invalid proof limits,
    /// allocation failure, or observed cancellation.  Numerical proof failure
    /// is a typed [`LayoutCertificateStatus::Unavailable`] result.
    pub fn certify_optimum(
        &self,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        limits: LayoutCertificateLimits,
        cx: &Cx<'_>,
    ) -> Result<LayoutCertificateStatus, LayoutCertificateError> {
        match build_certificate(self, x, y, settings, limits, cx) {
            Ok(certificate) => Ok(LayoutCertificateStatus::Certified(certificate)),
            Err(ProofFailure::Refusal(reason)) => Ok(unavailable(reason)),
            Err(ProofFailure::Error(error)) => Err(error),
        }
    }

    /// Certify one solver output and retain its verified bounds in the report
    /// that was minted with exactly the same LP, settings, and returned state.
    ///
    /// The report is updated only after a complete certificate passes; an
    /// unavailable, mismatched, or cancelled attempt clears older bounds.
    ///
    /// # Errors
    /// In addition to [`LayoutLp::certify_optimum`] errors, rejects a report
    /// whose private solve identity or public final diagnostics do not match
    /// `self`, `x`, `y`, and `settings`.
    pub fn certify_optimum_for_report(
        &self,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        report: &mut PdhgReport,
        limits: LayoutCertificateLimits,
        cx: &Cx<'_>,
    ) -> Result<LayoutCertificateStatus, LayoutCertificateError> {
        report.clear_certified_bounds();
        let status = self.certify_optimum(x, y, settings, limits, cx)?;
        if let LayoutCertificateStatus::Certified(certificate) = &status {
            let expected_identity = solver_state_identity_checked(self, x, y, settings, Some(cx))?;
            if !report.matches_solver_state_checked(expected_identity, || {
                poll(cx, "certificate report-snapshot hash")
            })? {
                return Err(LayoutCertificateError::ReportMismatch {
                    requirement: "must match this LP, settings, returned state, and final diagnostics",
                });
            }
            let bounds = certificate.bounds();
            poll(cx, "certificate report publication")?;
            report.retain_certified_bounds(
                bounds.lower(),
                bounds.upper(),
                certificate.dual_scale(),
                *certificate.certificate_identity().as_bytes(),
            );
        }
        Ok(status)
    }

    /// Re-run the mathematical proof and compare the complete receipt.
    ///
    /// This is the authoritative verifier for retained in-memory
    /// certificates; [`LayoutOptimalityCertificate::verifies_for`] is only a
    /// cheap identity/structure preflight.
    ///
    /// # Errors
    /// Returns the same hard errors as [`LayoutLp::certify_optimum`].  A clean
    /// `false` means the proof now refuses or its recomputed receipt differs.
    pub fn verify_optimum_certificate(
        &self,
        certificate: &LayoutOptimalityCertificate,
        x: &[f64],
        y: &[f64],
        settings: PdhgSettings,
        limits: LayoutCertificateLimits,
        cx: &Cx<'_>,
    ) -> Result<bool, LayoutCertificateError> {
        match build_certificate(self, x, y, settings, limits, cx) {
            Ok(recomputed) => {
                let problem_identity = problem_identity(self, Some(cx))?;
                let input_identity = input_identity(
                    problem_identity,
                    x,
                    y,
                    settings,
                    certificate.correction_method,
                    certificate.enclosure_method,
                    certificate.limits,
                    Some(cx),
                )?;
                Ok(
                    recomputed.certificate_identity == certificate.certificate_identity
                        && certificate.problem_identity == problem_identity
                        && certificate.input_identity == input_identity
                        && certificate_is_structurally_valid_checked(certificate, Some(cx))?
                        && certificate_identity(certificate, Some(cx))?
                            == certificate.certificate_identity,
                )
            }
            Err(ProofFailure::Refusal(_)) => Ok(false),
            Err(ProofFailure::Error(error)) => Err(error),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_fixture() -> LayoutOptimalityCertificate {
        let zero = LayoutCertificateIdentity(ContentHash([0; 32]));
        let mut certificate = LayoutOptimalityCertificate {
            bounds: CertifiedObjectiveBounds {
                lower: 0.0,
                upper: 1.0,
            },
            repaired_member_forces: vec![Interval::point(0.0)],
            repaired_split_forces: vec![Interval::point(0.0), Interval::point(0.0)],
            equilibrium_residuals: vec![Interval::point(0.0)],
            scaled_dual: vec![0.0],
            dual_slacks: vec![Interval::point(1.0), Interval::point(1.0)],
            dual_scale: 0.0,
            selected_members: vec![0],
            contraction_bound: 0.5,
            problem_identity: zero,
            input_identity: zero,
            certificate_identity: zero,
            limits: LayoutCertificateLimits::default(),
            correction_method: PrimalCorrectionMethod::SignedForceBasisNeumannV1,
            enclosure_method: ArithmeticEnclosureMethod::FsIvlIntervalV1,
        };
        certificate.certificate_identity =
            certificate_identity(&certificate, None).expect("infallible local identity");
        certificate
    }

    #[test]
    fn private_receipt_detects_endpoint_and_witness_tamper() {
        let certificate = identity_fixture();
        assert!(certificate.verifies_identity_checked(None).unwrap());

        let mut endpoint_tamper = certificate.clone();
        endpoint_tamper.bounds.upper = f64::from_bits(endpoint_tamper.bounds.upper.to_bits() + 1);
        assert!(!endpoint_tamper.verifies_identity_checked(None).unwrap());

        let mut witness_tamper = certificate;
        witness_tamper.scaled_dual[0] = f64::EPSILON;
        assert!(!witness_tamper.verifies_identity_checked(None).unwrap());
    }

    #[test]
    fn verification_work_excess_is_rejected_before_identity_hashing() {
        let matrix = Csr::from_parts(1, 2, vec![0, 2], vec![0, 1], vec![1.0, -1.0]);
        let costs = [1.0, 1.0];
        let loads = [1.0];
        let problem = LayoutCertificateProblem::try_new(&matrix, &costs, &loads)
            .expect("dimensionally valid split problem");
        let mut certificate = identity_fixture();
        certificate.limits =
            LayoutCertificateLimits::try_new(1, 1, 1, 1).expect("positive retained limits");
        certificate.certificate_identity =
            certificate_identity(&certificate, None).expect("fixture identity");
        IDENTITY_HASH_CALLS.with(|calls| calls.set(0));

        assert!(
            !certificate
                .verifies_for_problem_view_checked(
                    &problem,
                    &[0.0, 0.0],
                    &[0.0],
                    PdhgSettings::default(),
                    None,
                )
                .expect("work excess is a clean verification miss")
        );
        IDENTITY_HASH_CALLS.with(|calls| assert_eq!(calls.get(), 0));
    }
}
