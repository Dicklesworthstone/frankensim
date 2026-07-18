//! Typed verification journeys and scoped receipt algebra (V.4.1).
//!
//! This module binds user intent before execution: one journey identity,
//! one claim stratum, one manifest-resolved profile, the Five Explicits,
//! an artifact sandbox, a public-surface identity, workload scale, and
//! inspectable defaults.  The phase graph is explicit and fail-closed.
//!
//! Operation, attempt, job, and campaign receipts are distinct types.
//! Each owns its own execution and requested-predicate axes; a receipt may
//! content-reference another scope but cannot mutate or reinterpret it.
//! Scientific truth remains orthogonal to execution, evidence integrity,
//! evidence completeness, domain applicability, operational support, and
//! promotion.  Process codes project only an [`OperationReceipt`].
//!
//! No-claims: these records specify and content-address intent and results.
//! They do not execute a journey, authenticate evidence, or prove a claim.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_domain};

use crate::FiveExplicits;
use crate::v1::{
    ClaimRelationReceipt, ClaimRevision, ClaimRevisionId, JourneyId, NormalizedGraph, V1Error,
    admit_graph,
};
use crate::v1_selection::{ProfileId, Stratum};

const JOURNEY_DOMAIN: &str = "org.frankensim.fs-vmanifest.journey.v1";
const CLAIM_RECORD_DOMAIN: &str = "org.frankensim.fs-vmanifest.claim-record.v1";
const OPERATION_RECEIPT_DOMAIN: &str = "org.frankensim.fs-vmanifest.operation-receipt.v1";
const ATTEMPT_RECEIPT_DOMAIN: &str = "org.frankensim.fs-vmanifest.attempt-receipt.v1";
const JOB_RECEIPT_DOMAIN: &str = "org.frankensim.fs-vmanifest.job-receipt.v1";
const CAMPAIGN_RECEIPT_DOMAIN: &str = "org.frankensim.fs-vmanifest.campaign-receipt.v1";
/// Domain for the sealed v0-to-v1 scientific-axis migration receipt.
pub const LEGACY_ASSESSMENT_MIGRATION_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-vmanifest.legacy-assessment-migration.v1";
/// Schema version of the sealed v0-to-v1 scientific-axis migration receipt.
pub const LEGACY_ASSESSMENT_MIGRATION_SCHEMA_VERSION: u32 = 1;
const MAX_ID_BYTES: usize = 128;
const MAX_TEXT_BYTES: usize = 4096;
const MAX_REFERENCES: usize = 65_536;

/// Stable refusal from the Journey DSL or receipt algebra.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JourneyError {
    rule: &'static str,
    detail: String,
}

impl JourneyError {
    fn new(rule: &'static str, detail: impl Into<String>) -> Self {
        Self {
            rule,
            detail: detail.into(),
        }
    }

    /// Stable rule slug suitable for agent handling.
    #[must_use]
    pub const fn rule(&self) -> &'static str {
        self.rule
    }

    /// Human-readable refusal detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for JourneyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule, self.detail)
    }
}

impl std::error::Error for JourneyError {}

impl From<V1Error> for JourneyError {
    fn from(value: V1Error) -> Self {
        Self::new(value.rule(), value.to_string())
    }
}

fn validate_id(kind: &'static str, value: &str) -> Result<(), JourneyError> {
    if value.is_empty() || value.len() > MAX_ID_BYTES {
        return Err(JourneyError::new(
            "journey-id-bounds",
            format!(
                "{kind} id length {} outside 1..={MAX_ID_BYTES}",
                value.len()
            ),
        ));
    }
    if !value
        .bytes()
        .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b"-_.:/".contains(&b))
    {
        return Err(JourneyError::new(
            "journey-id-bounds",
            format!("{kind} id {value:?} outside [a-z0-9-_.:/]"),
        ));
    }
    Ok(())
}

fn validate_text(field: &'static str, value: &str) -> Result<(), JourneyError> {
    if value.trim().is_empty() || value.len() > MAX_TEXT_BYTES {
        return Err(JourneyError::new(
            "journey-field-bounds",
            format!(
                "{field} length {} outside the nonblank 1..={MAX_TEXT_BYTES} range",
                value.len()
            ),
        ));
    }
    Ok(())
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u32).to_be_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_hash(bytes: &mut Vec<u8>, value: ContentHash) {
    bytes.extend_from_slice(value.as_bytes());
}

fn validate_profile(profile: &ProfileId) -> Result<(), JourneyError> {
    if let ProfileId::Composite(composite) = profile {
        validate_id("composite profile", &composite.id)?;
        if composite.version == 0 {
            return Err(JourneyError::new(
                "journey-profile-version",
                "composite profile version must be nonzero",
            ));
        }
        composite.validate()?;
    }
    Ok(())
}

/// Production phases admitted by a verification journey.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum JourneyPhase {
    /// Discover available public surfaces and journey templates.
    Discover = 0,
    /// Check identities, capabilities, budgets, and local prerequisites.
    Preflight = 1,
    /// Author a new journey manifest.
    Author = 2,
    /// Import an existing typed manifest.
    Import = 3,
    /// Validate schema and semantic admission constraints.
    Validate = 4,
    /// Estimate cost, error, and resource envelopes.
    Estimate = 5,
    /// Freeze the execution and evidence plan.
    Plan = 6,
    /// Submit the planned request.
    Submit = 7,
    /// Apply capability, budget, and policy admission.
    Admit = 8,
    /// Await an execution allocation.
    Queue = 9,
    /// Execute admitted work.
    Execute = 10,
    /// Observe live bounded telemetry and progress.
    Observe = 11,
    /// Persist a resumable checkpoint.
    Checkpoint = 12,
    /// Pause after draining to a declared boundary.
    Pause = 13,
    /// Move resumable state to another admitted execution context.
    Migrate = 14,
    /// Request cancellation and drain/finalize its scope.
    Cancel = 15,
    /// Resume from an admitted checkpoint.
    Resume = 16,
    /// Fork a new lineage from an admitted checkpoint.
    Fork = 17,
    /// Inspect retained receipts and artifacts.
    Inspect = 18,
    /// Verify requested predicates and scientific evidence.
    Verify = 19,
    /// Produce a bounded result report.
    Report = 20,
    /// Apply publication policy and share admitted artifacts.
    Share = 21,
    /// Replay retained intent and provenance.
    Replay = 22,
}

impl JourneyPhase {
    /// Complete stable phase inventory.
    pub const ALL: [Self; 23] = [
        Self::Discover,
        Self::Preflight,
        Self::Author,
        Self::Import,
        Self::Validate,
        Self::Estimate,
        Self::Plan,
        Self::Submit,
        Self::Admit,
        Self::Queue,
        Self::Execute,
        Self::Observe,
        Self::Checkpoint,
        Self::Pause,
        Self::Migrate,
        Self::Cancel,
        Self::Resume,
        Self::Fork,
        Self::Inspect,
        Self::Verify,
        Self::Report,
        Self::Share,
        Self::Replay,
    ];

    /// Whether the production graph admits this direct transition.
    #[must_use]
    pub const fn allows(self, next: Self) -> bool {
        matches!(
            (self, next),
            (Self::Discover, Self::Preflight)
                | (Self::Preflight, Self::Author | Self::Import)
                | (Self::Author | Self::Import, Self::Validate)
                | (Self::Validate, Self::Estimate | Self::Plan)
                | (Self::Estimate, Self::Plan)
                | (Self::Plan, Self::Submit)
                | (Self::Submit, Self::Admit)
                | (Self::Admit, Self::Queue)
                | (Self::Queue, Self::Execute | Self::Cancel)
                | (
                    Self::Execute,
                    Self::Observe | Self::Checkpoint | Self::Pause | Self::Cancel
                )
                | (
                    Self::Observe,
                    Self::Execute
                        | Self::Checkpoint
                        | Self::Pause
                        | Self::Cancel
                        | Self::Inspect
                        | Self::Verify
                )
                | (
                    Self::Checkpoint,
                    Self::Execute
                        | Self::Pause
                        | Self::Migrate
                        | Self::Cancel
                        | Self::Resume
                        | Self::Fork
                )
                | (
                    Self::Pause,
                    Self::Resume | Self::Migrate | Self::Fork | Self::Cancel
                )
                | (Self::Migrate, Self::Resume | Self::Cancel)
                | (Self::Cancel, Self::Inspect | Self::Report | Self::Replay)
                | (Self::Resume, Self::Execute | Self::Observe)
                | (Self::Fork, Self::Plan | Self::Submit)
                | (Self::Inspect, Self::Verify | Self::Report | Self::Replay)
                | (Self::Verify, Self::Report | Self::Share | Self::Replay)
                | (Self::Report, Self::Share | Self::Replay)
                | (Self::Share, Self::Replay)
                | (Self::Replay, Self::Inspect | Self::Verify | Self::Report)
        )
    }
}

/// Stateful phase cursor that retains the exact transition history.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JourneyCursor {
    journey: JourneyId,
    phase: JourneyPhase,
    history: Vec<JourneyPhase>,
}

impl JourneyCursor {
    /// Start a validated manifest in `discover`.
    pub fn new(manifest: &JourneyManifest) -> Result<Self, JourneyError> {
        manifest.validate()?;
        Ok(Self {
            journey: manifest.journey.clone(),
            phase: JourneyPhase::Discover,
            history: vec![JourneyPhase::Discover],
        })
    }

    /// Advance exactly one legal edge; illegal jumps leave the cursor unchanged.
    pub fn transition(&mut self, next: JourneyPhase) -> Result<(), JourneyError> {
        if !self.phase.allows(next) {
            return Err(JourneyError::new(
                "journey-illegal-transition",
                format!("{:?} cannot transition directly to {next:?}", self.phase),
            ));
        }
        self.phase = next;
        self.history.push(next);
        Ok(())
    }

    #[must_use]
    /// Current phase.
    pub const fn phase(&self) -> JourneyPhase {
        self.phase
    }

    #[must_use]
    /// Exact admitted phase history, including the initial discovery phase.
    pub fn history(&self) -> &[JourneyPhase] {
        &self.history
    }

    #[must_use]
    /// Governing stable journey identity.
    pub fn journey(&self) -> &JourneyId {
        &self.journey
    }
}

/// Bounded relative output namespace for one journey.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArtifactSandbox {
    /// Normalized path relative to the admitted artifact root.
    pub relative_root: String,
    /// Maximum artifact objects admitted for the journey.
    pub max_artifacts: u32,
    /// Maximum cumulative artifact bytes admitted for the journey.
    pub max_bytes: u64,
    /// Explicit retention rule bound into manifest identity.
    pub retention_policy: String,
}

impl ArtifactSandbox {
    /// Validate relative-path normalization and nonzero bounds.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_text("artifact-sandbox root", &self.relative_root)?;
        if self.relative_root.starts_with('/')
            || self.relative_root.starts_with('\\')
            || self.relative_root.contains('\\')
            || self.relative_root.contains(':')
            || self
                .relative_root
                .split('/')
                .any(|part| part.is_empty() || part == "." || part == "..")
        {
            return Err(JourneyError::new(
                "journey-artifact-sandbox",
                "artifact root must be a normalized relative path without parent traversal",
            ));
        }
        if self.max_artifacts == 0 || self.max_bytes == 0 {
            return Err(JourneyError::new(
                "journey-artifact-sandbox",
                "artifact count and byte budgets must both be nonzero",
            ));
        }
        validate_text("artifact retention policy", &self.retention_policy)
    }
}

/// Content identity of the public command/catalog surface used by the journey.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublicSurfaceIdentity {
    /// Stable public catalog/grammar id.
    pub catalog: String,
    /// Version of the public schema.
    pub schema_version: u32,
    /// Content digest of the exact public surface.
    pub content_digest: ContentHash,
}

impl PublicSurfaceIdentity {
    /// Validate the catalog id and nonzero schema version.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_id("public-surface catalog", &self.catalog)?;
        if self.schema_version == 0 {
            return Err(JourneyError::new(
                "journey-public-surface",
                "public-surface schema version must be nonzero",
            ));
        }
        Ok(())
    }
}

/// Workload magnitude, deliberately orthogonal to campaign profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkloadScale {
    /// Logical cases selected before shard placement.
    pub logical_cases: u64,
    /// Exact number of deterministic shards.
    pub shards: u32,
    /// Maximum concurrent workers requested by this scale.
    pub max_concurrency: u32,
}

impl WorkloadScale {
    /// Refuse zero-valued scale axes.
    pub fn validate(self) -> Result<(), JourneyError> {
        if self.logical_cases == 0 || self.shards == 0 || self.max_concurrency == 0 {
            return Err(JourneyError::new(
                "journey-workload-scale",
                "logical cases, shards, and max concurrency must all be nonzero",
            ));
        }
        Ok(())
    }
}

/// Inspectable defaults; every duration has an explicit unit in its name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JourneyDefaults {
    /// Default per-operation timeout in nanoseconds.
    pub operation_timeout_ns: u64,
    /// Default cancellation-drain timeout in nanoseconds.
    pub drain_timeout_ns: u64,
    /// Default maximum attempt count.
    pub max_attempts: u32,
}

impl JourneyDefaults {
    /// Refuse defaults whose bounded quantities are zero.
    pub fn validate(self) -> Result<(), JourneyError> {
        if self.operation_timeout_ns == 0 || self.drain_timeout_ns == 0 || self.max_attempts == 0 {
            return Err(JourneyError::new(
                "journey-defaults",
                "operation timeout, drain timeout, and attempt count must be nonzero",
            ));
        }
        Ok(())
    }
}

/// Frozen user intent for one verification journey.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JourneyManifest {
    /// Stable identity of the verification journey.
    pub journey: JourneyId,
    /// Journey DSL schema version.
    pub schema_version: u32,
    /// Mandatory units, seeds, budgets, versions, and capabilities.
    pub explicits: FiveExplicits,
    /// Bounded artifact namespace.
    pub artifact_sandbox: ArtifactSandbox,
    /// Exact public command/catalog surface.
    pub public_surface: PublicSurfaceIdentity,
    /// Scientific claim surface, orthogonal to profile intensity.
    pub stratum: Stratum,
    /// Exactly one atomic or manifest-defined composite profile.
    pub profile: ProfileId,
    /// Digest of the deterministic selection expansion and semantic diff.
    pub selection_digest: ContentHash,
    /// Workload magnitude, separate from profile identity.
    pub scale: WorkloadScale,
    /// Inspectable default operation policy.
    pub defaults: JourneyDefaults,
}

impl JourneyManifest {
    /// Validate every manifest authority and bound.
    pub fn validate(&self) -> Result<(), JourneyError> {
        if self.schema_version == 0 {
            return Err(JourneyError::new(
                "journey-schema-version",
                "journey schema version must be nonzero",
            ));
        }
        for (field, value) in [
            ("units", self.explicits.units),
            ("seeds", self.explicits.seeds),
            ("budgets", self.explicits.budgets),
            ("versions", self.explicits.versions),
            ("capabilities", self.explicits.capabilities),
        ] {
            validate_text(field, value).map_err(|_| {
                JourneyError::new(
                    "journey-five-explicits",
                    format!("Five Explicits field {field} is missing or out of bounds"),
                )
            })?;
        }
        self.artifact_sandbox.validate()?;
        self.public_surface.validate()?;
        self.scale.validate()?;
        self.defaults.validate()?;
        validate_profile(&self.profile)
    }

    /// Domain-separated identity of every semantic manifest field.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        self.validate()?;
        let mut bytes = Vec::new();
        push_text(&mut bytes, self.journey.as_str());
        bytes.extend_from_slice(&self.schema_version.to_be_bytes());
        for value in [
            self.explicits.units,
            self.explicits.seeds,
            self.explicits.budgets,
            self.explicits.versions,
            self.explicits.capabilities,
        ] {
            push_text(&mut bytes, value);
        }
        push_text(&mut bytes, &self.artifact_sandbox.relative_root);
        bytes.extend_from_slice(&self.artifact_sandbox.max_artifacts.to_be_bytes());
        bytes.extend_from_slice(&self.artifact_sandbox.max_bytes.to_be_bytes());
        push_text(&mut bytes, &self.artifact_sandbox.retention_policy);
        push_text(&mut bytes, &self.public_surface.catalog);
        bytes.extend_from_slice(&self.public_surface.schema_version.to_be_bytes());
        push_hash(&mut bytes, self.public_surface.content_digest);
        push_text(&mut bytes, self.stratum.name());
        push_text(&mut bytes, &self.profile.render());
        if let ProfileId::Composite(composite) = &self.profile {
            push_hash(&mut bytes, composite.digest());
        }
        push_hash(&mut bytes, self.selection_digest);
        bytes.extend_from_slice(&self.scale.logical_cases.to_be_bytes());
        bytes.extend_from_slice(&self.scale.shards.to_be_bytes());
        bytes.extend_from_slice(&self.scale.max_concurrency.to_be_bytes());
        bytes.extend_from_slice(&self.defaults.operation_timeout_ns.to_be_bytes());
        bytes.extend_from_slice(&self.defaults.drain_timeout_ns.to_be_bytes());
        bytes.extend_from_slice(&self.defaults.max_attempts.to_be_bytes());
        Ok(hash_domain(JOURNEY_DOMAIN, &bytes))
    }
}

/// One exact claim revision plus the exact graph needed to interpret it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimRecord {
    /// Exact revision being adjudicated.
    pub subject: ClaimRevisionId,
    /// Revisions needed to interpret the subject's relation graph.
    pub revisions: Vec<ClaimRevision>,
    /// Typed, checker-receipted graph edges between exact revisions.
    pub relations: Vec<ClaimRelationReceipt>,
}

impl ClaimRecord {
    /// Validate and normalize the exact relation graph.
    pub fn normalized_graph(&self) -> Result<NormalizedGraph, JourneyError> {
        let graph = admit_graph(&self.revisions, &self.relations)?;
        if !graph.revisions().contains(&self.subject) {
            return Err(JourneyError::new(
                "journey-claim-subject",
                "claim subject is not one of the exact admitted revisions",
            ));
        }
        Ok(graph)
    }

    /// Claim identity includes the exact subject and normalized relation graph.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        let graph = self.normalized_graph()?;
        let mut bytes = Vec::new();
        push_hash(&mut bytes, self.subject);
        push_hash(&mut bytes, graph.digest());
        Ok(hash_domain(CLAIM_RECORD_DOMAIN, &bytes))
    }

    /// Return the exact subject revision when present.
    #[must_use]
    pub fn subject_revision(&self) -> Option<&ClaimRevision> {
        self.revisions.iter().find(|revision| {
            revision
                .revision_id()
                .is_ok_and(|revision_id| revision_id == self.subject)
        })
    }
}

/// Execution state at exactly one receipt scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ExecutionDisposition {
    /// Scope exists but has not entered a queue.
    Pending = 0,
    /// Scope is admitted and queued.
    Queued = 1,
    /// Scope is actively executing.
    Running = 2,
    /// Scope drained to a resumable paused boundary.
    Paused = 3,
    /// Scope state migrated and awaits/undergoes resume.
    Migrated = 4,
    /// Scope reached its normal completion boundary.
    Completed = 5,
    /// Scope completed with a non-infrastructure operational failure.
    Failed = 6,
    /// Cancellation won and all scoped work drained/finalized.
    CancelledAndDrained = 7,
    /// Independent timeout won and finalization completed.
    TimedOutFinalized = 8,
    /// Supervisor, transport, authentication, or publication failed.
    InfrastructureFailed = 9,
    /// Declared total resource budget was exhausted and finalized.
    BudgetExhaustedFinalized = 10,
}

/// Truth of the predicate requested by exactly one receipt scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum RequestedPredicateOutcome {
    /// The exact requested operation predicate is true.
    Satisfied = 0,
    /// The exact requested operation predicate is false.
    Unsatisfied = 1,
    /// Schema validation or admission rejected the request.
    InvalidSchemaOrAdmission = 2,
    /// Truth is unknown or required operation evidence is incomplete.
    IndeterminateOrIncomplete = 3,
    /// Domain applicability or admitted capabilities do not support it.
    Unsupported = 4,
    /// Cancellation won and its scope drained.
    CancelledAndDrained = 5,
    /// Timeout won and the scope finalized.
    TimeoutFinalized = 6,
    /// Infrastructure prevented a trustworthy operation result.
    InfrastructureError = 7,
    /// Integrity or security checks failed.
    IntegrityOrSecurityFailure = 8,
    /// The operation exhausted its declared total budget and finalized.
    BudgetExhaustedFinalized = 9,
}

/// The two orthogonal operational axes common to each distinct receipt type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReceiptOutcome {
    /// Execution disposition at this receipt's scope only.
    pub execution: ExecutionDisposition,
    /// Requested predicate at this receipt's scope only.
    pub predicate: RequestedPredicateOutcome,
}

impl ReceiptOutcome {
    /// Refuse favorable contradictions between terminal execution and predicate.
    pub fn validate(self) -> Result<(), JourneyError> {
        let required = match self.execution {
            ExecutionDisposition::Pending
            | ExecutionDisposition::Queued
            | ExecutionDisposition::Running
            | ExecutionDisposition::Paused
            | ExecutionDisposition::Migrated => {
                Some(RequestedPredicateOutcome::IndeterminateOrIncomplete)
            }
            ExecutionDisposition::CancelledAndDrained => {
                Some(RequestedPredicateOutcome::CancelledAndDrained)
            }
            ExecutionDisposition::TimedOutFinalized => {
                Some(RequestedPredicateOutcome::TimeoutFinalized)
            }
            ExecutionDisposition::InfrastructureFailed => {
                Some(RequestedPredicateOutcome::InfrastructureError)
            }
            ExecutionDisposition::BudgetExhaustedFinalized => {
                Some(RequestedPredicateOutcome::BudgetExhaustedFinalized)
            }
            ExecutionDisposition::Completed | ExecutionDisposition::Failed => None,
        };
        if required.is_some_and(|required| required != self.predicate) {
            return Err(JourneyError::new(
                "journey-outcome-inconsistent",
                format!(
                    "execution {:?} requires predicate {:?}, not {:?}",
                    self.execution, required, self.predicate
                ),
            ));
        }
        Ok(())
    }
}

/// Scientific adjudication; deliberately has no `Unsupported` or `Partial` state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ClaimAdjudication {
    /// No scientific verdict has been requested or completed.
    Pending = 0,
    /// Admitted evidence supports the exact claim revision.
    Supported = 1,
    /// The frozen acceptance predicate failed without a matched refutation.
    Failed = 2,
    /// An admitted strength-matched counterexample refutes the exact revision.
    Refuted = 3,
    /// Admitted information does not determine the claim.
    Unknown = 4,
}

/// Method family that produced admitted evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum EvidenceMethod {
    /// Focused local/unit checks.
    Unit = 0,
    /// Generative algebraic/property checks.
    Property = 1,
    /// Manufactured-solution convergence evidence.
    ManufacturedSolution = 2,
    /// Canonical benchmark evidence.
    Benchmark = 3,
    /// Metamorphic relation evidence.
    Metamorphic = 4,
    /// Fault, cancellation, and chaos evidence.
    Chaos = 5,
    /// Cross-run/thread/ISA determinism audit.
    DeterminismAudit = 6,
    /// Evidence from an independent admitted oracle.
    IndependentOracle = 7,
    /// Machine-checkable formal proof evidence.
    FormalProof = 8,
}

/// Canonically ordered set of admitted evidence methods.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct EvidenceMethodSet(BTreeSet<EvidenceMethod>);

impl EvidenceMethodSet {
    /// Canonicalize an arbitrary method iterator into a set.
    #[must_use]
    pub fn new(methods: impl IntoIterator<Item = EvidenceMethod>) -> Self {
        Self(methods.into_iter().collect())
    }

    /// Iterate in stable enum order.
    #[must_use]
    pub fn iter(&self) -> impl Iterator<Item = EvidenceMethod> + '_ {
        self.0.iter().copied()
    }

    /// Whether no evidence method was admitted.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Epistemic strength independent of claim truth and execution success.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum EpistemicGrade {
    /// No epistemic authority.
    None = 0,
    /// Exploratory evidence only.
    Exploratory = 1,
    /// Multiple admitted methods corroborate the claim.
    Corroborated = 2,
    /// Frozen certification obligations are discharged.
    Certified = 3,
}

/// Whether the exact frozen claim domain applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum DomainApplicability {
    /// Input and conditions are inside the exact claim domain.
    Admitted = 0,
    /// Well-formed input lies outside the exact claim domain.
    OutsideDomain = 1,
    /// Applicability cannot be determined.
    Indeterminate = 2,
}

/// Whether admitted capabilities implement the requested operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum OperationalSupport {
    /// Admitted capabilities implement the requested operation.
    Supported = 0,
    /// At least one required admitted capability is unavailable.
    MissingCapability = 1,
}

/// Presence of evidence required by the frozen acceptance card.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum EvidenceCompleteness {
    /// Every required evidence component is present.
    Complete = 0,
    /// Some, but not all, required evidence is present.
    Partial = 1,
    /// No required evidence is present.
    None = 2,
}

/// Integrity of evidence bytes, custody, lineage, and checker execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum EvidenceIntegrity {
    /// Bytes, custody, lineage, and checker execution verified.
    Verified = 0,
    /// At least one applicable integrity predicate failed.
    Failed = 1,
    /// Integrity could not be determined.
    Unknown = 2,
}

/// Exact authority transition caused by a scoped receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum PromotionEffect {
    /// Admit promotion to the targeted authority tier.
    Promote = 0,
    /// Preserve current authority without promotion.
    Hold = 1,
    /// Block promotion while retaining current authority.
    Block = 2,
    /// Revoke authority previously granted to this exact revision/scope.
    Revoke = 3,
}

/// Exact scientific axes carried by job and campaign receipts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScientificAssessment {
    /// Exact subject revision and its typed relation graph.
    pub claim: ClaimRecord,
    /// Scientific adjudication, independent of execution/support.
    pub adjudication: ClaimAdjudication,
    /// Exact admitted evidence methods.
    pub methods: EvidenceMethodSet,
    /// Epistemic strength of the admitted evidence set.
    pub grade: EpistemicGrade,
    /// Applicability to the frozen claim domain.
    pub domain: DomainApplicability,
    /// Capability support for the requested operation.
    pub support: OperationalSupport,
    /// Presence of every required evidence component.
    pub completeness: EvidenceCompleteness,
    /// Integrity of evidence, custody, lineage, and checker execution.
    pub integrity: EvidenceIntegrity,
    /// Exact effect on promotion authority.
    pub promotion: PromotionEffect,
}

impl ScientificAssessment {
    /// Validate the claim graph and fail-closed promotion conjunction.
    pub fn validate(&self) -> Result<(), JourneyError> {
        self.claim.normalized_graph()?;
        if self.promotion == PromotionEffect::Promote
            && !(self.adjudication == ClaimAdjudication::Supported
                && self.domain == DomainApplicability::Admitted
                && self.support == OperationalSupport::Supported
                && self.completeness == EvidenceCompleteness::Complete
                && self.integrity == EvidenceIntegrity::Verified
                && !self.methods.is_empty()
                && matches!(
                    self.grade,
                    EpistemicGrade::Corroborated | EpistemicGrade::Certified
                ))
        {
            return Err(JourneyError::new(
                "journey-promotion-laundering",
                "promotion requires supported, applicable, operationally supported, complete, verified, nonempty corroborated evidence",
            ));
        }
        Ok(())
    }

    fn encode(&self, bytes: &mut Vec<u8>) -> Result<(), JourneyError> {
        push_hash(bytes, self.claim.digest()?);
        bytes.push(self.adjudication as u8);
        bytes.extend_from_slice(&(self.methods.0.len() as u32).to_be_bytes());
        for method in self.methods.iter() {
            bytes.push(method as u8);
        }
        bytes.push(self.grade as u8);
        bytes.push(self.domain as u8);
        bytes.push(self.support as u8);
        bytes.push(self.completeness as u8);
        bytes.push(self.integrity as u8);
        bytes.push(self.promotion as u8);
        Ok(())
    }
}

/// Collapsed scientific verdict stored by the legacy v0 journey schema.
///
/// `Partial` and `Unsupported` deliberately have no direct counterparts in
/// [`ClaimAdjudication`]. The migration adapter retains them as unresolved
/// loss instead of manufacturing a weaker claim or guessing whether support
/// failed because of domain applicability or capability availability.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum LegacyClaimStateV0 {
    /// No scientific adjudication had completed.
    Pending = 0,
    /// The legacy aggregate reported support.
    Supported = 1,
    /// The frozen acceptance predicate failed without a refutation claim.
    Failed = 2,
    /// The legacy aggregate reported a strength-matched refutation.
    Refuted = 3,
    /// Available information did not determine the claim.
    Unknown = 4,
    /// The old aggregate collapsed a weaker claim and/or partial evidence.
    Partial = 5,
    /// The old aggregate collapsed outside-domain and missing-capability causes.
    Unsupported = 6,
}

/// Combined evidence state stored by the legacy v0 journey schema.
///
/// Separate completeness and integrity axes can be recovered only when the
/// legacy variant explicitly carried that fact. Missing axes remain absent in
/// the migration receipt and are listed in its loss set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum LegacyCombinedEvidenceV0 {
    /// No required evidence was present; integrity was not recorded.
    Absent = 0,
    /// Some required evidence was present; integrity was not recorded.
    Partial = 1,
    /// All required evidence was present; integrity was not recorded.
    Complete = 2,
    /// Integrity was verified; completeness was not recorded.
    Verified = 3,
    /// Completeness and integrity were both explicitly affirmative.
    CompleteAndVerified = 4,
    /// Integrity failed; completeness was not recorded.
    IntegrityFailed = 5,
}

/// Exact semantic information unavailable in a legacy migration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum LegacyAssessmentMigrationLoss {
    /// `Partial` did not distinguish a weaker claim from incomplete evidence.
    PartialClaimMeaningCollapsed = 0,
    /// `Unsupported` did not distinguish domain from capability support.
    UnsupportedCauseCollapsed = 1,
    /// The legacy evidence state did not determine completeness.
    EvidenceCompletenessMissing = 2,
    /// The legacy evidence state did not determine integrity.
    EvidenceIntegrityMissing = 3,
    /// The legacy aggregate did not identify evidence methods.
    EvidenceMethodsMissing = 4,
    /// The legacy aggregate did not retain an epistemic grade.
    EpistemicGradeMissing = 5,
    /// The legacy aggregate did not retain domain applicability separately.
    DomainApplicabilityMissing = 6,
    /// The legacy aggregate did not retain operational support separately.
    OperationalSupportMissing = 7,
    /// The legacy aggregate did not retain a promotion effect.
    PromotionEffectMissing = 8,
}

/// Sealed, content-addressed report of a legacy assessment migration.
///
/// Optional v1 axes contain only facts exactly recoverable from v0. Every
/// absent or collapsed fact is named in `losses`; consequently this receipt
/// cannot itself be converted into [`ScientificAssessment`] or promotion
/// authority. A higher layer must supply and receipt independent resolution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyAssessmentMigration {
    legacy_claim_state: LegacyClaimStateV0,
    legacy_evidence: LegacyCombinedEvidenceV0,
    adjudication: Option<ClaimAdjudication>,
    methods: Option<EvidenceMethodSet>,
    grade: Option<EpistemicGrade>,
    domain: Option<DomainApplicability>,
    support: Option<OperationalSupport>,
    completeness: Option<EvidenceCompleteness>,
    integrity: Option<EvidenceIntegrity>,
    promotion: Option<PromotionEffect>,
    losses: BTreeSet<LegacyAssessmentMigrationLoss>,
}

impl LegacyAssessmentMigration {
    /// Exact legacy claim state retained by the migration receipt.
    #[must_use]
    pub const fn legacy_claim_state(&self) -> LegacyClaimStateV0 {
        self.legacy_claim_state
    }

    /// Exact legacy combined evidence state retained by the receipt.
    #[must_use]
    pub const fn legacy_evidence(&self) -> LegacyCombinedEvidenceV0 {
        self.legacy_evidence
    }

    /// Exactly recoverable v1 adjudication, if the old state was not collapsed.
    #[must_use]
    pub const fn adjudication(&self) -> Option<ClaimAdjudication> {
        self.adjudication
    }

    /// Exactly recoverable evidence methods; v0 did not carry them.
    #[must_use]
    pub fn methods(&self) -> Option<&EvidenceMethodSet> {
        self.methods.as_ref()
    }

    /// Exactly recoverable epistemic grade; v0 did not carry it.
    #[must_use]
    pub const fn grade(&self) -> Option<EpistemicGrade> {
        self.grade
    }

    /// Exactly recoverable domain applicability; v0 did not separate it.
    #[must_use]
    pub const fn domain(&self) -> Option<DomainApplicability> {
        self.domain
    }

    /// Exactly recoverable operational support; v0 did not separate it.
    #[must_use]
    pub const fn support(&self) -> Option<OperationalSupport> {
        self.support
    }

    /// Exactly recoverable evidence completeness.
    #[must_use]
    pub const fn completeness(&self) -> Option<EvidenceCompleteness> {
        self.completeness
    }

    /// Exactly recoverable evidence integrity.
    #[must_use]
    pub const fn integrity(&self) -> Option<EvidenceIntegrity> {
        self.integrity
    }

    /// Exactly recoverable promotion effect; v0 did not carry it.
    #[must_use]
    pub const fn promotion(&self) -> Option<PromotionEffect> {
        self.promotion
    }

    /// Canonically ordered missing/collapsed semantic facts.
    #[must_use]
    pub fn losses(&self) -> impl Iterator<Item = LegacyAssessmentMigrationLoss> + '_ {
        self.losses.iter().copied()
    }

    /// Whether every axis required for a v1 scientific assessment was retained.
    #[must_use]
    pub fn is_authority_complete(&self) -> bool {
        self.adjudication.is_some()
            && self
                .methods
                .as_ref()
                .is_some_and(|methods| !methods.is_empty())
            && self.grade.is_some()
            && self.domain.is_some()
            && self.support.is_some()
            && self.completeness.is_some()
            && self.integrity.is_some()
            && self.promotion.is_some()
            && self.losses.is_empty()
    }

    /// Domain-separated identity of the exact old states, recovered axes, and
    /// complete loss set.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&LEGACY_ASSESSMENT_MIGRATION_SCHEMA_VERSION.to_be_bytes());
        bytes.push(self.legacy_claim_state as u8);
        bytes.push(self.legacy_evidence as u8);
        encode_optional_tag(&mut bytes, self.adjudication.map(|value| value as u8));
        match &self.methods {
            Some(methods) => {
                bytes.push(1);
                bytes.extend_from_slice(&(methods.0.len() as u32).to_be_bytes());
                for method in methods.iter() {
                    bytes.push(method as u8);
                }
            }
            None => bytes.push(0),
        }
        encode_optional_tag(&mut bytes, self.grade.map(|value| value as u8));
        encode_optional_tag(&mut bytes, self.domain.map(|value| value as u8));
        encode_optional_tag(&mut bytes, self.support.map(|value| value as u8));
        encode_optional_tag(&mut bytes, self.completeness.map(|value| value as u8));
        encode_optional_tag(&mut bytes, self.integrity.map(|value| value as u8));
        encode_optional_tag(&mut bytes, self.promotion.map(|value| value as u8));
        bytes.extend_from_slice(&(self.losses.len() as u32).to_be_bytes());
        for loss in &self.losses {
            bytes.push(*loss as u8);
        }
        hash_domain(LEGACY_ASSESSMENT_MIGRATION_IDENTITY_DOMAIN, &bytes)
    }
}

/// Migrate the collapsed legacy v0 claim/evidence pair without guessing.
///
/// The returned sealed receipt preserves all exactly recoverable axes and
/// lists every unavailable or ambiguous axis. In particular `Partial` never
/// becomes a favorable adjudication and `Unsupported` never chooses between
/// outside-domain and missing-capability.
#[must_use]
pub fn migrate_legacy_assessment_v0(
    claim_state: LegacyClaimStateV0,
    evidence: LegacyCombinedEvidenceV0,
) -> LegacyAssessmentMigration {
    let mut losses = BTreeSet::from([
        LegacyAssessmentMigrationLoss::EvidenceMethodsMissing,
        LegacyAssessmentMigrationLoss::EpistemicGradeMissing,
        LegacyAssessmentMigrationLoss::DomainApplicabilityMissing,
        LegacyAssessmentMigrationLoss::OperationalSupportMissing,
        LegacyAssessmentMigrationLoss::PromotionEffectMissing,
    ]);
    let adjudication = match claim_state {
        LegacyClaimStateV0::Pending => Some(ClaimAdjudication::Pending),
        LegacyClaimStateV0::Supported => Some(ClaimAdjudication::Supported),
        LegacyClaimStateV0::Failed => Some(ClaimAdjudication::Failed),
        LegacyClaimStateV0::Refuted => Some(ClaimAdjudication::Refuted),
        LegacyClaimStateV0::Unknown => Some(ClaimAdjudication::Unknown),
        LegacyClaimStateV0::Partial => {
            losses.insert(LegacyAssessmentMigrationLoss::PartialClaimMeaningCollapsed);
            None
        }
        LegacyClaimStateV0::Unsupported => {
            losses.insert(LegacyAssessmentMigrationLoss::UnsupportedCauseCollapsed);
            None
        }
    };
    let (completeness, integrity) = match evidence {
        LegacyCombinedEvidenceV0::Absent => (Some(EvidenceCompleteness::None), None),
        LegacyCombinedEvidenceV0::Partial => (Some(EvidenceCompleteness::Partial), None),
        LegacyCombinedEvidenceV0::Complete => (Some(EvidenceCompleteness::Complete), None),
        LegacyCombinedEvidenceV0::Verified => (None, Some(EvidenceIntegrity::Verified)),
        LegacyCombinedEvidenceV0::CompleteAndVerified => (
            Some(EvidenceCompleteness::Complete),
            Some(EvidenceIntegrity::Verified),
        ),
        LegacyCombinedEvidenceV0::IntegrityFailed => (None, Some(EvidenceIntegrity::Failed)),
    };
    if completeness.is_none() {
        losses.insert(LegacyAssessmentMigrationLoss::EvidenceCompletenessMissing);
    }
    if integrity.is_none() {
        losses.insert(LegacyAssessmentMigrationLoss::EvidenceIntegrityMissing);
    }
    LegacyAssessmentMigration {
        legacy_claim_state: claim_state,
        legacy_evidence: evidence,
        adjudication,
        methods: None,
        grade: None,
        domain: None,
        support: None,
        completeness,
        integrity,
        promotion: None,
        losses,
    }
}

/// Retained conditional skip; never represented by a missing row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TypedSkip {
    /// Stable identity of the retained skip.
    pub id: String,
    /// Exact conditional predicate that triggered the skip.
    pub predicate: String,
    /// Human/agent-readable reason.
    pub reason: String,
    /// Authority responsible for discharging or reviewing it.
    pub owner: String,
    /// Explicit effect on promotion while the skip remains live.
    pub promotion_effect: PromotionEffect,
}

impl TypedSkip {
    fn validate(&self) -> Result<(), JourneyError> {
        validate_id("skip", &self.id)?;
        validate_text("skip predicate", &self.predicate)?;
        validate_text("skip reason", &self.reason)?;
        validate_text("skip owner", &self.owner)
    }

    fn encode(&self, bytes: &mut Vec<u8>) {
        push_text(bytes, &self.id);
        push_text(bytes, &self.predicate);
        push_text(bytes, &self.reason);
        push_text(bytes, &self.owner);
        bytes.push(self.promotion_effect as u8);
    }
}

/// User-requested predicate family for an operation receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum OperationVerb {
    /// Query retained state without re-adjudicating it.
    Status = 0,
    /// Request cancellation; acceptance precedes terminal drain evidence.
    Cancel = 1,
    /// Request evidence that supports the exact claim.
    Prove = 2,
    /// Request a determinate adjudication, including valid refutation.
    Adjudicate = 3,
    /// Execute admitted planned work.
    Execute = 4,
    /// Verify an exact requested predicate.
    Verify = 5,
    /// Apply publication policy and share artifacts.
    Share = 6,
    /// Replay retained intent and provenance.
    Replay = 7,
}

/// Stable process projection of the current operation only.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ProcessCode {
    /// Requested operation predicate satisfied.
    Satisfied = 0,
    /// Requested operation predicate unsatisfied.
    Unsatisfied = 10,
    /// Invalid schema or admission.
    InvalidSchemaOrAdmission = 11,
    /// Indeterminate result or incomplete required operation evidence.
    IndeterminateOrIncomplete = 12,
    /// Requested operation unsupported.
    Unsupported = 13,
    /// Requested operation cancelled and drained.
    CancelledAndDrained = 14,
    /// Requested operation timed out and finalized.
    TimeoutFinalized = 15,
    /// Infrastructure failure prevented a trustworthy result.
    InfrastructureError = 16,
    /// Integrity or security failure prevented a trustworthy result.
    IntegrityOrSecurityFailure = 17,
    /// Requested operation exhausted its budget and finalized.
    BudgetExhaustedFinalized = 18,
}

impl ProcessCode {
    /// Numeric process status used by command surfaces.
    #[must_use]
    pub const fn value(self) -> u8 {
        self as u8
    }
}

/// Receipt for exactly one user-requested operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperationReceipt {
    /// Stable id of this user-requested operation.
    pub operation_id: String,
    /// Journey whose intent governs the operation.
    pub journey: JourneyId,
    /// Exact requested operation predicate family.
    pub verb: OperationVerb,
    /// Phase in which the operation was receipted.
    pub phase: JourneyPhase,
    /// Operation-scoped execution and predicate axes.
    pub outcome: ReceiptOutcome,
    /// Optional immutable reference to another scoped receipt.
    pub referenced_receipt: Option<ContentHash>,
    /// Every conditional skip retained by this operation.
    pub skips: Vec<TypedSkip>,
}

impl OperationReceipt {
    /// Validate id, scoped outcome, caps, and retained skips.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_id("operation", &self.operation_id)?;
        self.outcome.validate()?;
        if self.skips.len() > MAX_REFERENCES {
            return Err(JourneyError::new(
                "journey-receipt-cap",
                "operation skip count exceeds the receipt cap",
            ));
        }
        let mut ids = BTreeSet::new();
        for skip in &self.skips {
            skip.validate()?;
            if !ids.insert(skip.id.as_str()) {
                return Err(JourneyError::new(
                    "journey-duplicate-skip",
                    format!("duplicate typed skip {:?}", skip.id),
                ));
            }
        }
        Ok(())
    }

    /// Create a scoped query/adjudication operation referencing, not embedding,
    /// a job receipt. Status and cancellation acceptance describe the operation
    /// itself and therefore remain satisfied even when the job is not green.
    pub fn for_job(
        operation_id: impl Into<String>,
        verb: OperationVerb,
        phase: JourneyPhase,
        job: &JobReceipt,
    ) -> Result<Self, JourneyError> {
        job.validate()?;
        let predicate = match verb {
            OperationVerb::Status | OperationVerb::Cancel => RequestedPredicateOutcome::Satisfied,
            OperationVerb::Prove | OperationVerb::Adjudicate => predicate_for_science(verb, job),
            OperationVerb::Execute
            | OperationVerb::Verify
            | OperationVerb::Share
            | OperationVerb::Replay => job.outcome.predicate,
        };
        let receipt = Self {
            operation_id: operation_id.into(),
            journey: job.journey.clone(),
            verb,
            phase,
            outcome: ReceiptOutcome {
                execution: ExecutionDisposition::Completed,
                predicate,
            },
            referenced_receipt: Some(job.digest()?),
            skips: Vec::new(),
        };
        receipt.validate()?;
        Ok(receipt)
    }

    /// Domain-separated content identity of the operation receipt.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        self.validate()?;
        let mut bytes = Vec::new();
        push_text(&mut bytes, &self.operation_id);
        push_text(&mut bytes, self.journey.as_str());
        bytes.push(self.verb as u8);
        bytes.push(self.phase as u8);
        bytes.push(self.outcome.execution as u8);
        bytes.push(self.outcome.predicate as u8);
        match self.referenced_receipt {
            Some(reference) => {
                bytes.push(1);
                push_hash(&mut bytes, reference);
            }
            None => bytes.push(0),
        }
        bytes.extend_from_slice(&(self.skips.len() as u32).to_be_bytes());
        for skip in &self.skips {
            skip.encode(&mut bytes);
        }
        Ok(hash_domain(OPERATION_RECEIPT_DOMAIN, &bytes))
    }

    /// Project only this operation's requested predicate. Referenced job or
    /// campaign receipts are deliberately not traversed.
    #[must_use]
    pub const fn process_code(&self) -> ProcessCode {
        match self.outcome.predicate {
            RequestedPredicateOutcome::Satisfied => ProcessCode::Satisfied,
            RequestedPredicateOutcome::Unsatisfied => ProcessCode::Unsatisfied,
            RequestedPredicateOutcome::InvalidSchemaOrAdmission => {
                ProcessCode::InvalidSchemaOrAdmission
            }
            RequestedPredicateOutcome::IndeterminateOrIncomplete => {
                ProcessCode::IndeterminateOrIncomplete
            }
            RequestedPredicateOutcome::Unsupported => ProcessCode::Unsupported,
            RequestedPredicateOutcome::CancelledAndDrained => ProcessCode::CancelledAndDrained,
            RequestedPredicateOutcome::TimeoutFinalized => ProcessCode::TimeoutFinalized,
            RequestedPredicateOutcome::InfrastructureError => ProcessCode::InfrastructureError,
            RequestedPredicateOutcome::IntegrityOrSecurityFailure => {
                ProcessCode::IntegrityOrSecurityFailure
            }
            RequestedPredicateOutcome::BudgetExhaustedFinalized => {
                ProcessCode::BudgetExhaustedFinalized
            }
        }
    }

    /// Human projection bound to the same receipt digest as the JSONL surface.
    pub fn render_pretty(&self) -> Result<String, JourneyError> {
        Ok(format!(
            "operation={} verb={:?} phase={:?} code={} digest={}",
            self.operation_id,
            self.verb,
            self.phase,
            self.process_code().value(),
            self.digest()?.to_hex()
        ))
    }

    /// Bounded machine projection. Restricted ids make its string fields
    /// JSON-safe without a second escaping grammar.
    pub fn render_json_line(&self) -> Result<String, JourneyError> {
        Ok(format!(
            "{{\"schema\":\"journey-operation-v1\",\"operation\":\"{}\",\"verb\":\"{:?}\",\"phase\":\"{:?}\",\"process_code\":{},\"digest\":\"{}\"}}",
            self.operation_id,
            self.verb,
            self.phase,
            self.process_code().value(),
            self.digest()?.to_hex()
        ))
    }
}

fn predicate_for_science(verb: OperationVerb, job: &JobReceipt) -> RequestedPredicateOutcome {
    if job.science.integrity == EvidenceIntegrity::Failed {
        return RequestedPredicateOutcome::IntegrityOrSecurityFailure;
    }
    match job.outcome.execution {
        ExecutionDisposition::CancelledAndDrained => {
            return RequestedPredicateOutcome::CancelledAndDrained;
        }
        ExecutionDisposition::TimedOutFinalized => {
            return RequestedPredicateOutcome::TimeoutFinalized;
        }
        ExecutionDisposition::InfrastructureFailed => {
            return RequestedPredicateOutcome::InfrastructureError;
        }
        ExecutionDisposition::BudgetExhaustedFinalized => {
            return RequestedPredicateOutcome::BudgetExhaustedFinalized;
        }
        _ => {}
    }
    if job.science.domain == DomainApplicability::OutsideDomain
        || job.science.support == OperationalSupport::MissingCapability
    {
        return RequestedPredicateOutcome::Unsupported;
    }
    if job.science.domain == DomainApplicability::Indeterminate
        || job.science.integrity == EvidenceIntegrity::Unknown
        || job.science.completeness != EvidenceCompleteness::Complete
    {
        return RequestedPredicateOutcome::IndeterminateOrIncomplete;
    }
    match (verb, job.science.adjudication) {
        (OperationVerb::Prove, ClaimAdjudication::Supported) => {
            RequestedPredicateOutcome::Satisfied
        }
        (OperationVerb::Prove, ClaimAdjudication::Failed | ClaimAdjudication::Refuted) => {
            RequestedPredicateOutcome::Unsatisfied
        }
        (
            OperationVerb::Adjudicate,
            ClaimAdjudication::Supported | ClaimAdjudication::Failed | ClaimAdjudication::Refuted,
        ) => RequestedPredicateOutcome::Satisfied,
        _ => RequestedPredicateOutcome::IndeterminateOrIncomplete,
    }
}

/// Receipt for one concrete execution attempt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttemptReceipt {
    /// Stable id of the concrete attempt.
    pub attempt_id: String,
    /// Stable parent job id.
    pub job_id: String,
    /// Governing journey id.
    pub journey: JourneyId,
    /// Exact attempt-level requested predicate id.
    pub requested_predicate: String,
    /// Phase at which this attempt receipt was finalized or observed.
    pub phase: JourneyPhase,
    /// Attempt-scoped execution and predicate axes.
    pub outcome: ReceiptOutcome,
    /// Optional immutable predecessor attempt reference.
    pub parent_attempt: Option<ContentHash>,
    /// Ordered content identities of produced artifacts.
    pub artifacts: Vec<ContentHash>,
}

impl AttemptReceipt {
    /// Validate ids, scoped outcome, bounds, and duplicate references.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_id("attempt", &self.attempt_id)?;
        validate_id("job", &self.job_id)?;
        validate_id("attempt requested predicate", &self.requested_predicate)?;
        self.outcome.validate()?;
        validate_hash_references("attempt artifacts", &self.artifacts)
    }

    /// Domain-separated content identity of the attempt receipt.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        self.validate()?;
        let mut bytes = Vec::new();
        push_text(&mut bytes, &self.attempt_id);
        push_text(&mut bytes, &self.job_id);
        push_text(&mut bytes, self.journey.as_str());
        push_text(&mut bytes, &self.requested_predicate);
        bytes.push(self.phase as u8);
        encode_outcome(&mut bytes, self.outcome);
        encode_optional_hash(&mut bytes, self.parent_attempt);
        encode_hashes(&mut bytes, &self.artifacts);
        Ok(hash_domain(ATTEMPT_RECEIPT_DOMAIN, &bytes))
    }
}

/// Receipt aggregating attempts for exactly one job and one exact claim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JobReceipt {
    /// Stable id of the aggregate job.
    pub job_id: String,
    /// Governing journey id.
    pub journey: JourneyId,
    /// Exact job-level requested predicate id.
    pub requested_predicate: String,
    /// Phase at which this job receipt was finalized or observed.
    pub phase: JourneyPhase,
    /// Job-scoped execution and predicate axes.
    pub outcome: ReceiptOutcome,
    /// Ordered immutable references to attempt receipts.
    pub attempts: Vec<ContentHash>,
    /// Exact scientific axes for this job scope.
    pub science: ScientificAssessment,
    /// Every conditional skip retained by this job.
    pub skips: Vec<TypedSkip>,
}

impl JobReceipt {
    /// Validate ids, scope axes, references, claim science, and skips.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_id("job", &self.job_id)?;
        validate_id("job requested predicate", &self.requested_predicate)?;
        self.outcome.validate()?;
        validate_hash_references("job attempts", &self.attempts)?;
        self.science.validate()?;
        validate_skips(&self.skips)
    }

    /// Domain-separated content identity of the job receipt.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        self.validate()?;
        let mut bytes = Vec::new();
        push_text(&mut bytes, &self.job_id);
        push_text(&mut bytes, self.journey.as_str());
        push_text(&mut bytes, &self.requested_predicate);
        bytes.push(self.phase as u8);
        encode_outcome(&mut bytes, self.outcome);
        encode_hashes(&mut bytes, &self.attempts);
        self.science.encode(&mut bytes)?;
        encode_skips(&mut bytes, &self.skips);
        Ok(hash_domain(JOB_RECEIPT_DOMAIN, &bytes))
    }
}

/// Receipt aggregating jobs for one manifest-resolved campaign run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CampaignReceipt {
    /// Stable id of this one-profile campaign run.
    pub campaign_id: String,
    /// Governing journey id.
    pub journey: JourneyId,
    /// Exact manifest authority used for selection and execution.
    pub manifest_digest: ContentHash,
    /// Core or Max claim surface.
    pub stratum: Stratum,
    /// Exactly one manifest-resolved profile.
    pub profile: ProfileId,
    /// Exact deterministic selection expansion.
    pub selection_digest: ContentHash,
    /// Exact campaign-level aggregate predicate id.
    pub requested_predicate: String,
    /// Phase at which this campaign receipt was finalized or observed.
    pub phase: JourneyPhase,
    /// Campaign-scoped execution and predicate axes.
    pub outcome: ReceiptOutcome,
    /// Ordered immutable references to job receipts.
    pub jobs: Vec<ContentHash>,
    /// Exact scientific axes for this campaign scope.
    pub science: ScientificAssessment,
    /// Every conditional skip retained by this campaign.
    pub skips: Vec<TypedSkip>,
}

impl CampaignReceipt {
    /// Validate ids, selection/profile, scope axes, science, and skips.
    pub fn validate(&self) -> Result<(), JourneyError> {
        validate_id("campaign", &self.campaign_id)?;
        validate_id("campaign requested predicate", &self.requested_predicate)?;
        self.outcome.validate()?;
        validate_hash_references("campaign jobs", &self.jobs)?;
        validate_profile(&self.profile)?;
        self.science.validate()?;
        validate_skips(&self.skips)
    }

    /// Domain-separated content identity of the campaign receipt.
    pub fn digest(&self) -> Result<ContentHash, JourneyError> {
        self.validate()?;
        let mut bytes = Vec::new();
        push_text(&mut bytes, &self.campaign_id);
        push_text(&mut bytes, self.journey.as_str());
        push_hash(&mut bytes, self.manifest_digest);
        push_text(&mut bytes, self.stratum.name());
        push_text(&mut bytes, &self.profile.render());
        if let ProfileId::Composite(composite) = &self.profile {
            push_hash(&mut bytes, composite.digest());
        }
        push_hash(&mut bytes, self.selection_digest);
        push_text(&mut bytes, &self.requested_predicate);
        bytes.push(self.phase as u8);
        encode_outcome(&mut bytes, self.outcome);
        encode_hashes(&mut bytes, &self.jobs);
        self.science.encode(&mut bytes)?;
        encode_skips(&mut bytes, &self.skips);
        Ok(hash_domain(CAMPAIGN_RECEIPT_DOMAIN, &bytes))
    }
}

fn validate_hash_references(
    what: &'static str,
    references: &[ContentHash],
) -> Result<(), JourneyError> {
    if references.len() > MAX_REFERENCES {
        return Err(JourneyError::new(
            "journey-receipt-cap",
            format!("{what} exceeds {MAX_REFERENCES} entries"),
        ));
    }
    let mut seen = BTreeSet::new();
    for reference in references {
        if !seen.insert(*reference) {
            return Err(JourneyError::new(
                "journey-duplicate-reference",
                format!("{what} contains duplicate digest {}", reference.to_hex()),
            ));
        }
    }
    Ok(())
}

fn validate_skips(skips: &[TypedSkip]) -> Result<(), JourneyError> {
    if skips.len() > MAX_REFERENCES {
        return Err(JourneyError::new(
            "journey-receipt-cap",
            "typed skip count exceeds the receipt cap",
        ));
    }
    let mut seen = BTreeSet::new();
    for skip in skips {
        skip.validate()?;
        if !seen.insert(skip.id.as_str()) {
            return Err(JourneyError::new(
                "journey-duplicate-skip",
                format!("duplicate typed skip {:?}", skip.id),
            ));
        }
    }
    Ok(())
}

fn encode_outcome(bytes: &mut Vec<u8>, outcome: ReceiptOutcome) {
    bytes.push(outcome.execution as u8);
    bytes.push(outcome.predicate as u8);
}

fn encode_optional_tag(bytes: &mut Vec<u8>, value: Option<u8>) {
    match value {
        Some(value) => {
            bytes.push(1);
            bytes.push(value);
        }
        None => bytes.push(0),
    }
}

fn encode_optional_hash(bytes: &mut Vec<u8>, value: Option<ContentHash>) {
    match value {
        Some(value) => {
            bytes.push(1);
            push_hash(bytes, value);
        }
        None => bytes.push(0),
    }
}

fn encode_hashes(bytes: &mut Vec<u8>, values: &[ContentHash]) {
    bytes.extend_from_slice(&(values.len() as u32).to_be_bytes());
    for value in values {
        push_hash(bytes, *value);
    }
}

fn encode_skips(bytes: &mut Vec<u8>, skips: &[TypedSkip]) {
    bytes.extend_from_slice(&(skips.len() as u32).to_be_bytes());
    for skip in skips {
        skip.encode(bytes);
    }
}
