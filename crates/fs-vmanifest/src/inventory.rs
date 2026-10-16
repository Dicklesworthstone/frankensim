//! Pure V.1.2 verification-inventory compiler.
//!
//! Source-specific adapters normalize exact, content-addressed inputs into
//! this module's draft types. Compilation validates the complete input before
//! sealing one canonical inventory. It records disagreements; it never treats
//! tracker state, generated code, or an observed artifact as scientific proof.

use core::fmt::{self, Write as _};
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, DomainHasher};

use crate::v1::{
    CaseId, ClaimId, ClaimRelationReceipt, ClaimRevision, ClaimRevisionId, FieldSpec, JourneyId,
    MANIFEST_RECORD_FIELDS, MANIFEST_V1_SCHEMA_VERSION, NormalizedGraph, RelationKind,
    SourceAuthority, SourcePin, V1Error, admit_graph,
};

/// Schema version bound into every compiled inventory and projection row.
pub const INVENTORY_SCHEMA_VERSION: u32 = 1;

/// Implementation version of this normalized-input compiler.
pub const INVENTORY_COMPILER_VERSION: u32 = 1;

/// Version of the field-authority policy implemented by this module.
pub const INVENTORY_AUTHORITY_POLICY_VERSION: u32 = 1;

/// Version of the alias/rename/split/merge reconciliation policy.
pub const INVENTORY_RECONCILIATION_POLICY_VERSION: u32 = 1;

/// Maximum UTF-8 bytes in an inventory locator, adapter id, rationale, or
/// normalized field value.
pub const MAX_INVENTORY_TEXT_BYTES: usize = 4096;

const SOURCE_ID_DOMAIN: &str = "org.frankensim.fs-vmanifest.inventory-source.v1";
const SOURCE_SET_DOMAIN: &str = "org.frankensim.fs-vmanifest.inventory-source-set.v1";
const INVENTORY_ID_DOMAIN: &str = "org.frankensim.fs-vmanifest.inventory.v1";
const SEMANTIC_PROJECTION_DOMAIN: &str =
    "org.frankensim.fs-vmanifest.inventory-projection.semantic.v1";
const HUMAN_PROJECTION_DOMAIN: &str = "org.frankensim.fs-vmanifest.inventory-projection.human.v1";
const JSON_LINES_PROJECTION_DOMAIN: &str =
    "org.frankensim.fs-vmanifest.inventory-projection.json-lines.v1";
const LEDGER_PROJECTION_DOMAIN: &str = "org.frankensim.fs-vmanifest.inventory-projection.ledger.v1";
const DERIVED_SOURCE_RESOLUTION_BYTES: u64 = 128;
const DERIVED_REVISION_RESOLUTION_BYTES: u64 = 64;

/// Fixed protocol ceilings for one V.1.2 compilation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InventoryLimits {
    /// Maximum normalized source snapshots.
    pub max_sources: u32,
    /// Maximum exact claim revisions.
    pub max_revisions: u32,
    /// Maximum typed v1 relation receipts.
    pub max_relations: u32,
    /// Maximum source-bound field assertions.
    pub max_facts: u32,
    /// Maximum source-bound evidence observations.
    pub max_observations: u32,
    /// Maximum alias/rename/split/merge receipts.
    pub max_reconciliations: u32,
    /// Maximum lineage or exact-revision endpoints across reconciliation
    /// receipts. This separately bounds the vectors carried by split/merge.
    pub max_reconciliation_endpoints: u32,
    /// Maximum values in one complete field assertion.
    pub max_values_per_fact: u32,
    /// Maximum aggregate normalized UTF-8 and conservative derived-resolution
    /// logical bytes outside `ClaimRevision` and relation payloads (those are
    /// bounded by v1 graph admission).
    pub max_semantic_bytes: u64,
    /// Maximum rows in any complete semantic projection.
    pub max_projection_rows: u32,
}

impl InventoryLimits {
    /// Fixed V.1.2 protocol envelope. Callers may tighten but not loosen it.
    pub const DEFAULT: Self = Self {
        max_sources: 512,
        max_revisions: 4096,
        max_relations: 8192,
        max_facts: 65_536,
        max_observations: 65_536,
        max_reconciliations: 8192,
        max_reconciliation_endpoints: 32_768,
        max_values_per_fact: 256,
        max_semantic_bytes: 64 * 1024 * 1024,
        max_projection_rows: 262_144,
    };
}

impl Default for InventoryLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Machine-readable resource refusal evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InventoryResourceRefusal {
    /// Quantity that exceeded its envelope.
    pub quantity: &'static str,
    /// Exact required amount.
    pub required: u64,
    /// Admitted amount.
    pub admitted: u64,
    /// Stable unit label.
    pub unit: &'static str,
}

/// Structural or resource refusal. Semantic disagreements are retained in a
/// successfully compiled inventory's conflict set instead.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryRefusal {
    rule: &'static str,
    detail: String,
    fixes: Vec<String>,
    resource: Option<InventoryResourceRefusal>,
}

impl InventoryRefusal {
    fn new(rule: &'static str, detail: impl Into<String>) -> Self {
        Self {
            rule,
            detail: detail.into(),
            fixes: Vec::new(),
            resource: None,
        }
    }

    fn with_fix(mut self, fix: impl Into<String>) -> Self {
        self.fixes.push(fix.into());
        self
    }

    fn with_resource(mut self, resource: InventoryResourceRefusal) -> Self {
        self.resource = Some(resource);
        self
    }

    /// Stable rule id.
    #[must_use]
    pub const fn rule(&self) -> &'static str {
        self.rule
    }

    /// Human-readable refusal detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }

    /// Ranked repairs, most likely first.
    #[must_use]
    pub fn ranked_fixes(&self) -> &[String] {
        &self.fixes
    }

    /// Exact required/admitted resource evidence when applicable.
    #[must_use]
    pub const fn resource_refusal(&self) -> Option<InventoryResourceRefusal> {
        self.resource
    }
}

impl fmt::Display for InventoryRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule, self.detail)?;
        for (index, fix) in self.fixes.iter().enumerate() {
            write!(f, "\n  fix[{index}]: {fix}")?;
        }
        Ok(())
    }
}

impl std::error::Error for InventoryRefusal {}

/// Heterogeneous source class. The class fixes the source authority and
/// semantic role; callers cannot promote a generated artifact by relabeling
/// its public `SourcePin`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum InventorySourceKind {
    /// Normalized Beads obligation records.
    Beads,
    /// Crate or subsystem contract declarations.
    Contract,
    /// Typed IR or catalog registry declarations.
    TypedRegistry,
    /// Executable code registrations.
    CodeRegistration,
    /// Test, fixture, oracle, or checker registrations.
    TestRegistration,
    /// V&V artifact observations.
    VvArtifact,
    /// Benchmark registry observations.
    BenchmarkRegistry,
    /// Ledger receipt observations.
    LedgerReceipt,
    /// A previously sealed inventory used only for exact historical replay.
    FrozenInventory,
}

impl InventorySourceKind {
    /// Canonical source-kind label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Beads => "beads",
            Self::Contract => "contract",
            Self::TypedRegistry => "typed-registry",
            Self::CodeRegistration => "code-registration",
            Self::TestRegistration => "test-registration",
            Self::VvArtifact => "vv-artifact",
            Self::BenchmarkRegistry => "benchmark-registry",
            Self::LedgerReceipt => "ledger-receipt",
            Self::FrozenInventory => "frozen-inventory",
        }
    }

    /// Authority fixed by the source class.
    #[must_use]
    pub const fn authority(self) -> SourceAuthority {
        match self {
            Self::Beads => SourceAuthority::BeadObligation,
            Self::Contract => SourceAuthority::Contract,
            Self::TypedRegistry | Self::CodeRegistration | Self::LedgerReceipt => {
                SourceAuthority::GeneratedArtifact
            }
            Self::TestRegistration | Self::VvArtifact | Self::BenchmarkRegistry => {
                SourceAuthority::TestSource
            }
            Self::FrozenInventory => SourceAuthority::FrozenSnapshot,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Beads => 1,
            Self::Contract => 2,
            Self::TypedRegistry => 3,
            Self::CodeRegistration => 4,
            Self::TestRegistration => 5,
            Self::VvArtifact => 6,
            Self::BenchmarkRegistry => 7,
            Self::LedgerReceipt => 8,
            Self::FrozenInventory => 9,
        }
    }
}

/// Semantic role of a normalized source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum InventoryRole {
    /// Requirement and ownership obligations.
    Obligation,
    /// Contract-declared semantics and no-claim boundaries.
    DeclaredSemantics,
    /// Generated/registered executable surfaces.
    ExecutableRegistration,
    /// Tests, fixtures, oracles, and checker context.
    ValidationContext,
    /// Observed V&V, benchmark, or ledger evidence.
    ObservedEvidence,
    /// Exact replay of a previously compiled inventory.
    FrozenReplay,
}

impl InventoryRole {
    /// Canonical role label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Obligation => "obligation",
            Self::DeclaredSemantics => "declared-semantics",
            Self::ExecutableRegistration => "executable-registration",
            Self::ValidationContext => "validation-context",
            Self::ObservedEvidence => "observed-evidence",
            Self::FrozenReplay => "frozen-replay",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Obligation => 1,
            Self::DeclaredSemantics => 2,
            Self::ExecutableRegistration => 3,
            Self::ValidationContext => 4,
            Self::ObservedEvidence => 5,
            Self::FrozenReplay => 6,
        }
    }
}

/// Caller-normalized immutable source description. Fields are private so raw
/// callers cannot independently choose a class, role, and authority triple.
/// Non-frozen constructors remain an adapter trust boundary: the compiler
/// proves canonicalization and consistency, not that a caller read the named
/// filesystem source honestly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventorySourceDraft {
    kind: InventorySourceKind,
    role: InventoryRole,
    pin: SourcePin,
    adapter_version: String,
    adapter_policy_version: u32,
}

impl InventorySourceDraft {
    fn normalized(
        kind: InventorySourceKind,
        role: InventoryRole,
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self {
            kind,
            role,
            pin: SourcePin {
                source: source.into(),
                authority: kind.authority(),
                snapshot,
            },
            adapter_version: adapter_version.into(),
            adapter_policy_version,
        }
    }

    /// Normalize a Beads obligation snapshot.
    #[must_use]
    pub fn beads(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::Beads,
            InventoryRole::Obligation,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a crate/subsystem contract snapshot.
    #[must_use]
    pub fn contract(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::Contract,
            InventoryRole::DeclaredSemantics,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a typed IR/catalog registry snapshot.
    #[must_use]
    pub fn typed_registry(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::TypedRegistry,
            InventoryRole::ExecutableRegistration,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize an executable code-registration snapshot.
    #[must_use]
    pub fn code_registration(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::CodeRegistration,
            InventoryRole::ExecutableRegistration,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a test/fixture/oracle/checker registration snapshot.
    #[must_use]
    pub fn test_registration(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::TestRegistration,
            InventoryRole::ValidationContext,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a V&V artifact snapshot.
    #[must_use]
    pub fn vv_artifact(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::VvArtifact,
            InventoryRole::ObservedEvidence,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a benchmark-registry snapshot.
    #[must_use]
    pub fn benchmark_registry(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::BenchmarkRegistry,
            InventoryRole::ObservedEvidence,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Normalize a ledger-receipt snapshot.
    #[must_use]
    pub fn ledger_receipt(
        source: impl Into<String>,
        snapshot: ContentHash,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> Self {
        Self::normalized(
            InventorySourceKind::LedgerReceipt,
            InventoryRole::ObservedEvidence,
            source,
            snapshot,
            adapter_version,
            adapter_policy_version,
        )
    }

    /// Closed source class.
    #[must_use]
    pub const fn kind(&self) -> InventorySourceKind {
        self.kind
    }

    /// Fixed semantic role.
    #[must_use]
    pub const fn role(&self) -> InventoryRole {
        self.role
    }

    /// Exact class-derived authority and source snapshot.
    #[must_use]
    pub const fn pin(&self) -> &SourcePin {
        &self.pin
    }

    /// Adapter implementation/version id.
    #[must_use]
    pub fn adapter_version(&self) -> &str {
        &self.adapter_version
    }

    /// Source-specific normalization policy version.
    #[must_use]
    pub const fn adapter_policy_version(&self) -> u32 {
        self.adapter_policy_version
    }
}

/// Exact identity of one admitted normalized source description.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InventorySourceId(ContentHash);

impl InventorySourceId {
    /// Underlying domain-separated content hash.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.0
    }

    /// Lowercase hexadecimal encoding.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

impl fmt::Display for InventorySourceId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Digest of the complete canonical input source set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InventorySourceSetDigest(ContentHash);

impl InventorySourceSetDigest {
    /// Underlying domain-separated content hash.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.0
    }
}

impl fmt::Display for InventorySourceSetDigest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// Digest of one complete canonical inventory, including its conflicts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InventoryDigest(ContentHash);

impl InventoryDigest {
    /// Underlying domain-separated content hash.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.0
    }
}

impl fmt::Display for InventoryDigest {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// The exact 22-field VerificationManifest-v1 record registry as a typed key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum InventoryField {
    /// Exact source snapshot identities (compiler-derived).
    SourceSnapshots,
    /// Owning Bead obligation.
    BeadObligation,
    /// Exact claim revision id (compiler-derived).
    ClaimRevision,
    /// Core or max claim stratum.
    Stratum,
    /// Admitted campaign profiles.
    CampaignProfiles,
    /// S/F/M ambition tag.
    Ambition,
    /// Public code/contract surface.
    PublicSurface,
    /// Stable conformance case ids.
    CaseIds,
    /// Stable journey ids.
    JourneyIds,
    /// Explicit verification owner.
    Ownership,
    /// Fixture ids.
    FixtureIds,
    /// Oracle ids.
    OracleIds,
    /// Checker ids.
    CheckerIds,
    /// TCB overlap declaration.
    TcbOverlap,
    /// Tolerance derivation.
    ToleranceDerivation,
    /// Explicit budgets.
    Budgets,
    /// Required capabilities.
    Capabilities,
    /// Event kinds.
    EventKinds,
    /// Retention policy.
    Retention,
    /// Exact replay command.
    ReplayCommand,
    /// DSR lane.
    DsrLane,
    /// Scoped receipt expectations.
    ReceiptExpectations,
}

impl InventoryField {
    /// Complete registry order. Tests pin this one-to-one to
    /// `MANIFEST_RECORD_FIELDS`.
    pub const ALL: [Self; 22] = [
        Self::SourceSnapshots,
        Self::BeadObligation,
        Self::ClaimRevision,
        Self::Stratum,
        Self::CampaignProfiles,
        Self::Ambition,
        Self::PublicSurface,
        Self::CaseIds,
        Self::JourneyIds,
        Self::Ownership,
        Self::FixtureIds,
        Self::OracleIds,
        Self::CheckerIds,
        Self::TcbOverlap,
        Self::ToleranceDerivation,
        Self::Budgets,
        Self::Capabilities,
        Self::EventKinds,
        Self::Retention,
        Self::ReplayCommand,
        Self::DsrLane,
        Self::ReceiptExpectations,
    ];

    /// Canonical field name from the v1 registry.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::SourceSnapshots => "source-snapshots",
            Self::BeadObligation => "bead-obligation",
            Self::ClaimRevision => "claim-revision",
            Self::Stratum => "stratum",
            Self::CampaignProfiles => "campaign-profiles",
            Self::Ambition => "ambition",
            Self::PublicSurface => "public-surface",
            Self::CaseIds => "case-ids",
            Self::JourneyIds => "journey-ids",
            Self::Ownership => "ownership",
            Self::FixtureIds => "fixture-ids",
            Self::OracleIds => "oracle-ids",
            Self::CheckerIds => "checker-ids",
            Self::TcbOverlap => "tcb-overlap",
            Self::ToleranceDerivation => "tolerance-derivation",
            Self::Budgets => "budgets",
            Self::Capabilities => "capabilities",
            Self::EventKinds => "event-kinds",
            Self::Retention => "retention",
            Self::ReplayCommand => "replay-command",
            Self::DsrLane => "dsr-lane",
            Self::ReceiptExpectations => "receipt-expectations",
        }
    }

    /// Authoritative registry metadata for this field.
    #[must_use]
    pub fn spec(self) -> &'static FieldSpec {
        MANIFEST_RECORD_FIELDS
            .iter()
            .find(|spec| spec.name == self.as_str())
            .expect("InventoryField is pinned one-to-one to MANIFEST_RECORD_FIELDS")
    }

    /// Whether the compiler, rather than a source adapter, derives the field.
    #[must_use]
    pub const fn is_derived(self) -> bool {
        matches!(self, Self::SourceSnapshots | Self::ClaimRevision)
    }
}

/// One source's complete normalized assertion for one manifest-record field.
/// Values are canonicalized as a sorted set; adapters must represent order-
/// sensitive semantics through an already content-identified value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryFactDraft {
    /// Exact claim revision receiving the field.
    pub revision: ClaimRevisionId,
    /// Typed v1 record field.
    pub field: InventoryField,
    /// Complete asserted value set for the field.
    pub values: Vec<String>,
    /// Index into `InventoryDraft::sources` before source canonicalization.
    pub source_index: usize,
}

/// Terminal execution disposition retained independently from scientific
/// adjudication.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ObservationExecution {
    /// Work has not produced a terminal execution result.
    NotRun,
    /// The admitted execution completed.
    Completed,
    /// Execution failed operationally.
    Failed,
    /// Execution cancelled and drained.
    Cancelled,
    /// Execution finalized after its time budget.
    TimedOut,
    /// Execution exhausted an explicit resource budget.
    BudgetExhausted,
    /// Infrastructure prevented a trustworthy result.
    InfrastructureFailed,
}

impl ObservationExecution {
    const fn as_str(self) -> &'static str {
        match self {
            Self::NotRun => "not-run",
            Self::Completed => "completed",
            Self::Failed => "failed",
            Self::Cancelled => "cancelled",
            Self::TimedOut => "timed-out",
            Self::BudgetExhausted => "budget-exhausted",
            Self::InfrastructureFailed => "infrastructure-failed",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::NotRun => 1,
            Self::Completed => 2,
            Self::Failed => 3,
            Self::Cancelled => 4,
            Self::TimedOut => 5,
            Self::BudgetExhausted => 6,
            Self::InfrastructureFailed => 7,
        }
    }
}

/// Completeness of the exact retained evidence payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ObservationCompleteness {
    /// Every required evidence component is retained.
    Complete,
    /// A durable partial payload exists.
    Partial,
    /// No scientific evidence payload is available.
    NoEvidence,
}

impl ObservationCompleteness {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Partial => "partial",
            Self::NoEvidence => "no-evidence",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Complete => 1,
            Self::Partial => 2,
            Self::NoEvidence => 3,
        }
    }
}

/// Integrity of the retained evidence and its custody/checker envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ObservationIntegrity {
    /// Integrity checks completed successfully.
    Verified,
    /// Integrity, custody, or checker validation failed.
    Failed,
}

impl ObservationIntegrity {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Failed => "failed",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Verified => 1,
            Self::Failed => 2,
        }
    }
}

/// Scientific adjudication retained separately from execution and evidence
/// health. Compilation records this claim; it does not prove it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ObservationAdjudication {
    /// No terminal scientific adjudication.
    Pending,
    /// Evidence reports support for the exact revision.
    Supported,
    /// The frozen predicate reports failure without refutation.
    Failed,
    /// A reported counterexample refutes the exact revision.
    Refuted,
    /// Evidence does not determine the exact revision.
    Unknown,
}

impl ObservationAdjudication {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Pending => "pending",
            Self::Supported => "supported",
            Self::Failed => "failed",
            Self::Refuted => "refuted",
            Self::Unknown => "unknown",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Pending => 1,
            Self::Supported => 2,
            Self::Failed => 3,
            Self::Refuted => 4,
            Self::Unknown => 5,
        }
    }

    const fn requires_complete_evidence(self) -> bool {
        matches!(self, Self::Supported | Self::Failed | Self::Refuted)
    }
}

/// One caller-normalized evidence observation. The exact artifact digest and
/// source snapshot are retained; current metadata can never substitute for
/// them during replay.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryObservationDraft {
    /// Exact claim revision observed.
    pub revision: ClaimRevisionId,
    /// Stable observation identity within the source.
    pub observation_id: String,
    /// Exact evidence/failure-bundle/receipt artifact digest.
    pub artifact_digest: ContentHash,
    /// Orthogonal execution disposition.
    pub execution: ObservationExecution,
    /// Orthogonal evidence completeness.
    pub completeness: ObservationCompleteness,
    /// Orthogonal evidence integrity.
    pub integrity: ObservationIntegrity,
    /// Orthogonal scientific adjudication.
    pub adjudication: ObservationAdjudication,
    /// Index into `InventoryDraft::sources` before canonicalization.
    pub source_index: usize,
}

/// Alias/rename/split/merge topology is presentation/lineage metadata, not a
/// scientific relation and not a substitute for `ClaimRelationReceipt`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ReconciliationKind {
    /// Alternate name for one continuing claim lineage.
    Alias,
    /// Explicit presentation rename.
    Rename,
    /// One predecessor lineage split into two or more successors.
    Split,
    /// Two or more predecessor lineages merged into one successor.
    Merge,
}

impl ReconciliationKind {
    /// Canonical kind label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Alias => "alias",
            Self::Rename => "rename",
            Self::Split => "split",
            Self::Merge => "merge",
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Alias => 1,
            Self::Rename => 2,
            Self::Split => 3,
            Self::Merge => 4,
        }
    }
}

/// One exact reconciliation assertion. Presentation-only alias/rename edges
/// bind continuing lineages, while history-changing split/merge edges bind
/// immutable revision ids. The type prevents callers from confusing the two.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReconciliationDraft {
    /// An alternate lineage id names one canonical continuing lineage.
    Alias {
        /// Alternate lineage id.
        alias: ClaimId,
        /// Canonical lineage id.
        canonical: ClaimId,
        /// Source index before canonicalization.
        source_index: usize,
        /// Nonzero reconciliation policy version.
        policy_version: u32,
        /// Bounded rationale retained in identity and projections.
        rationale: String,
    },
    /// A presentation id was explicitly renamed without rewriting history.
    Rename {
        /// Former lineage id.
        previous: ClaimId,
        /// Current lineage id.
        current: ClaimId,
        /// Source index before canonicalization.
        source_index: usize,
        /// Nonzero reconciliation policy version.
        policy_version: u32,
        /// Bounded rationale retained in identity and projections.
        rationale: String,
    },
    /// One exact predecessor revision split into at least two exact successor
    /// revisions. Each successor must supersede the predecessor.
    Split {
        /// Exact predecessor revision.
        predecessor: ClaimRevisionId,
        /// Exact successor revisions.
        successors: Vec<ClaimRevisionId>,
        /// Source index before canonicalization.
        source_index: usize,
        /// Nonzero reconciliation policy version.
        policy_version: u32,
        /// Bounded rationale retained in identity and projections.
        rationale: String,
    },
    /// At least two exact predecessor revisions merged into one exact
    /// successor. The successor's single v1 `supersedes` anchor must name one
    /// of the complete predecessor set retained here.
    Merge {
        /// Exact predecessor revisions.
        predecessors: Vec<ClaimRevisionId>,
        /// Exact successor revision.
        successor: ClaimRevisionId,
        /// Source index before canonicalization.
        source_index: usize,
        /// Nonzero reconciliation policy version.
        policy_version: u32,
        /// Bounded rationale retained in identity and projections.
        rationale: String,
    },
}

impl ReconciliationDraft {
    const fn kind(&self) -> ReconciliationKind {
        match self {
            Self::Alias { .. } => ReconciliationKind::Alias,
            Self::Rename { .. } => ReconciliationKind::Rename,
            Self::Split { .. } => ReconciliationKind::Split,
            Self::Merge { .. } => ReconciliationKind::Merge,
        }
    }

    const fn source_index(&self) -> usize {
        match self {
            Self::Alias { source_index, .. }
            | Self::Rename { source_index, .. }
            | Self::Split { source_index, .. }
            | Self::Merge { source_index, .. } => *source_index,
        }
    }

    const fn policy_version(&self) -> u32 {
        match self {
            Self::Alias { policy_version, .. }
            | Self::Rename { policy_version, .. }
            | Self::Split { policy_version, .. }
            | Self::Merge { policy_version, .. } => *policy_version,
        }
    }

    fn rationale(&self) -> &str {
        match self {
            Self::Alias { rationale, .. }
            | Self::Rename { rationale, .. }
            | Self::Split { rationale, .. }
            | Self::Merge { rationale, .. } => rationale,
        }
    }

    fn endpoint_count(&self) -> Result<u64, InventoryRefusal> {
        match self {
            Self::Alias { .. } | Self::Rename { .. } => Ok(2),
            Self::Split { successors, .. } => checked_collection_with_anchor(successors.len()),
            Self::Merge { predecessors, .. } => checked_collection_with_anchor(predecessors.len()),
        }
    }
}

/// Complete pure compiler input. Concrete source adapters live elsewhere.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryDraft {
    /// Exact heterogeneous normalized source set.
    pub sources: Vec<InventorySourceDraft>,
    /// Every live and historical immutable claim revision.
    pub revisions: Vec<ClaimRevision>,
    /// Typed scientific relation receipts admitted by V.1.1.
    pub relations: Vec<ClaimRelationReceipt>,
    /// Source-bound complete manifest field assertions.
    pub facts: Vec<InventoryFactDraft>,
    /// Source-bound V&V, benchmark, or ledger observations.
    pub observations: Vec<InventoryObservationDraft>,
    /// Explicit alias/rename/split/merge assertions.
    pub reconciliations: Vec<ReconciliationDraft>,
    /// Authority/reconciliation versions expected by the caller.
    pub authority_policy_version: u32,
    /// Reconciliation policy version expected by the caller.
    pub reconciliation_policy_version: u32,
}

/// Canonically retained source.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventorySource {
    id: InventorySourceId,
    kind: InventorySourceKind,
    role: InventoryRole,
    pin: SourcePin,
    adapter_version: String,
    adapter_policy_version: u32,
}

impl InventorySource {
    /// Exact source identity.
    #[must_use]
    pub const fn id(&self) -> InventorySourceId {
        self.id
    }

    /// Closed source class.
    #[must_use]
    pub const fn kind(&self) -> InventorySourceKind {
        self.kind
    }

    /// Source semantic role.
    #[must_use]
    pub const fn role(&self) -> InventoryRole {
        self.role
    }

    /// Exact v1 source pin.
    #[must_use]
    pub const fn pin(&self) -> &SourcePin {
        &self.pin
    }

    /// Adapter implementation/version id.
    #[must_use]
    pub fn adapter_version(&self) -> &str {
        &self.adapter_version
    }

    /// Adapter policy version.
    #[must_use]
    pub const fn adapter_policy_version(&self) -> u32 {
        self.adapter_policy_version
    }
}

/// Canonically retained source-bound fact.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct InventoryFact {
    /// Exact claim revision.
    pub revision: ClaimRevisionId,
    /// Typed field.
    pub field: InventoryField,
    /// Canonical sorted, duplicate-free value set.
    pub values: Vec<String>,
    /// Exact normalized source identity.
    pub source: InventorySourceId,
    /// Source role retained for conflict diagnostics.
    pub role: InventoryRole,
    /// Source authority retained for field-policy checks.
    pub authority: SourceAuthority,
}

/// Canonically retained source-bound evidence observation. Construction is
/// sealed behind inventory admission so terminal adjudications with
/// contradictory evidence-health axes cannot be fabricated as admitted
/// observations.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct InventoryObservation {
    revision: ClaimRevisionId,
    observation_id: String,
    artifact_digest: ContentHash,
    execution: ObservationExecution,
    completeness: ObservationCompleteness,
    integrity: ObservationIntegrity,
    adjudication: ObservationAdjudication,
    source: InventorySourceId,
}

impl InventoryObservation {
    /// Exact observed claim revision.
    #[must_use]
    pub const fn revision(&self) -> ClaimRevisionId {
        self.revision
    }

    /// Stable source-local observation id.
    #[must_use]
    pub fn observation_id(&self) -> &str {
        &self.observation_id
    }

    /// Exact evidence artifact digest.
    #[must_use]
    pub const fn artifact_digest(&self) -> ContentHash {
        self.artifact_digest
    }

    /// Execution disposition.
    #[must_use]
    pub const fn execution(&self) -> ObservationExecution {
        self.execution
    }

    /// Evidence completeness.
    #[must_use]
    pub const fn completeness(&self) -> ObservationCompleteness {
        self.completeness
    }

    /// Evidence integrity.
    #[must_use]
    pub const fn integrity(&self) -> ObservationIntegrity {
        self.integrity
    }

    /// Scientific adjudication as reported by the source.
    #[must_use]
    pub const fn adjudication(&self) -> ObservationAdjudication {
        self.adjudication
    }

    /// Exact normalized source identity.
    #[must_use]
    pub const fn source(&self) -> InventorySourceId {
        self.source
    }
}

/// One resolved or explicitly unresolved v1 record field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryFieldResolution {
    /// Exact claim revision.
    pub revision: ClaimRevisionId,
    /// Typed field.
    pub field: InventoryField,
    /// Canonical value set, or `None` when a blocking conflict prevents a
    /// truthful choice.
    pub values: Option<Vec<String>>,
    /// Sources corroborating the selected value.
    pub sources: Vec<InventorySourceId>,
}

/// Conflict classification. Conflicts are identity-forming and never
/// last-writer-wins side effects.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum InventoryConflictKind {
    /// Required field has no assertion from its declared role/authority.
    MissingRequiredField,
    /// Equal-authority eligible sources asserted different values.
    EqualAuthorityDisagreement,
    /// A non-owning semantic role reported a value different from the
    /// authoritative resolved field value (for example stale code versus its
    /// contract declaration).
    CrossRoleDisagreement,
    /// Distinct complete verified observations report incompatible terminal
    /// scientific adjudications for the same exact revision.
    ObservationAdjudicationConflict,
    /// One logical source locator was supplied at incompatible snapshots.
    SourceSnapshotConflict,
    /// A higher-authority pin explicitly resolved an older pin for the same
    /// logical source; both remain in provenance.
    SourceSnapshotResolved,
    /// Competing complete alias/rename targets or split/merge endpoint sets.
    ReconciliationConflict,
}

impl InventoryConflictKind {
    /// Canonical conflict label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::MissingRequiredField => "missing-required-field",
            Self::EqualAuthorityDisagreement => "equal-authority-disagreement",
            Self::CrossRoleDisagreement => "cross-role-disagreement",
            Self::ObservationAdjudicationConflict => "observation-adjudication-conflict",
            Self::SourceSnapshotConflict => "source-snapshot-conflict",
            Self::SourceSnapshotResolved => "source-snapshot-resolved",
            Self::ReconciliationConflict => "reconciliation-conflict",
        }
    }

    /// Whether the conflict prevents a claim record from being complete.
    #[must_use]
    pub const fn blocking(self) -> bool {
        !matches!(self, Self::SourceSnapshotResolved)
    }

    const fn tag(self) -> u8 {
        match self {
            Self::MissingRequiredField => 1,
            Self::EqualAuthorityDisagreement => 2,
            Self::CrossRoleDisagreement => 3,
            Self::ObservationAdjudicationConflict => 4,
            Self::SourceSnapshotConflict => 5,
            Self::SourceSnapshotResolved => 6,
            Self::ReconciliationConflict => 7,
        }
    }
}

/// One deterministic retained conflict.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct InventoryConflict {
    /// Conflict class.
    pub kind: InventoryConflictKind,
    /// Affected revision, if record-local.
    pub revision: Option<ClaimRevisionId>,
    /// Affected v1 field, if field-local.
    pub field: Option<InventoryField>,
    /// Every involved exact source identity in canonical order.
    pub sources: Vec<InventorySourceId>,
    /// Stable actionable detail bound into inventory identity.
    pub detail: String,
}

/// Canonically retained typed reconciliation receipt.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReconciliationReceipt {
    topology: ReconciliationTopology,
    source: InventorySourceId,
    policy_version: u32,
    rationale: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum ReconciliationTopology {
    /// Presentation alias between continuing lineages.
    Alias {
        /// Alternate lineage id.
        alias: ClaimId,
        /// Canonical lineage id.
        canonical: ClaimId,
    },
    /// Presentation rename between continuing lineages.
    Rename {
        /// Former lineage id.
        previous: ClaimId,
        /// Current lineage id.
        current: ClaimId,
    },
    /// Exact one-to-many revision transition.
    Split {
        /// Exact predecessor revision.
        predecessor: ClaimRevisionId,
        /// Sorted, duplicate-free exact successor revisions.
        successors: Vec<ClaimRevisionId>,
    },
    /// Exact many-to-one revision transition.
    Merge {
        /// Sorted, duplicate-free exact predecessor revisions.
        predecessors: Vec<ClaimRevisionId>,
        /// Exact successor revision.
        successor: ClaimRevisionId,
    },
}

impl ReconciliationReceipt {
    /// Reconciliation topology.
    #[must_use]
    pub const fn kind(&self) -> ReconciliationKind {
        match &self.topology {
            ReconciliationTopology::Alias { .. } => ReconciliationKind::Alias,
            ReconciliationTopology::Rename { .. } => ReconciliationKind::Rename,
            ReconciliationTopology::Split { .. } => ReconciliationKind::Split,
            ReconciliationTopology::Merge { .. } => ReconciliationKind::Merge,
        }
    }

    /// Exact source that asserted the receipt.
    #[must_use]
    pub const fn source(&self) -> InventorySourceId {
        self.source
    }

    /// Reconciliation policy version.
    #[must_use]
    pub const fn policy_version(&self) -> u32 {
        self.policy_version
    }

    /// Human rationale retained in identity.
    #[must_use]
    pub fn rationale(&self) -> &str {
        &self.rationale
    }

    /// Alias endpoints, when this is an alias receipt.
    #[must_use]
    pub const fn alias(&self) -> Option<(&ClaimId, &ClaimId)> {
        match &self.topology {
            ReconciliationTopology::Alias { alias, canonical } => Some((alias, canonical)),
            _ => None,
        }
    }

    /// Rename endpoints, when this is a rename receipt.
    #[must_use]
    pub const fn rename(&self) -> Option<(&ClaimId, &ClaimId)> {
        match &self.topology {
            ReconciliationTopology::Rename { previous, current } => Some((previous, current)),
            _ => None,
        }
    }

    /// Split predecessor and complete successor set, when this is a split.
    #[must_use]
    pub fn split(&self) -> Option<(ClaimRevisionId, &[ClaimRevisionId])> {
        match &self.topology {
            ReconciliationTopology::Split {
                predecessor,
                successors,
            } => Some((*predecessor, successors)),
            _ => None,
        }
    }

    /// Complete predecessor set and successor, when this is a merge.
    #[must_use]
    pub fn merge(&self) -> Option<(&[ClaimRevisionId], ClaimRevisionId)> {
        match &self.topology {
            ReconciliationTopology::Merge {
                predecessors,
                successor,
            } => Some((predecessors, *successor)),
            _ => None,
        }
    }
}

/// Semantic projection row shared by all physical views.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventorySemanticRow {
    /// Global canonical ordinal.
    pub ordinal: u32,
    /// Stable row kind.
    pub kind: String,
    /// Stable subject identity/locator.
    pub subject: String,
    /// Field or relation label.
    pub field: String,
    /// Canonical semantic value.
    pub value: String,
}

/// Domain-separated hashes of the canonical semantic projection and each
/// complete physical view. Physical digests are intentionally distinct; they
/// prove exact bytes, not cross-format byte equality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InventoryProjectionDigests {
    /// Canonical framed semantic-row sequence.
    pub semantic: ContentHash,
    /// Exact human projection bytes.
    pub human: ContentHash,
    /// Exact strict JSON-lines projection bytes.
    pub json_lines: ContentHash,
    /// Exact ledger projection bytes.
    pub ledger: ContentHash,
}

/// Reviewable semantic inventory change classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum InventoryDiffKind {
    /// A semantic row key exists only in the newer inventory.
    Added,
    /// A semantic row key exists only in the older inventory.
    Removed,
    /// The same semantic row key has different canonical values.
    Changed,
}

impl InventoryDiffKind {
    /// Stable projection label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Added => "added",
            Self::Removed => "removed",
            Self::Changed => "changed",
        }
    }
}

/// One sorted semantic inventory diff entry. Values are vectors because
/// corroborating or conflicting sources may share one logical row key.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryDiffEntry {
    kind: InventoryDiffKind,
    row_kind: String,
    subject: String,
    field: String,
    before: Vec<String>,
    after: Vec<String>,
}

impl InventoryDiffEntry {
    /// Change class.
    #[must_use]
    pub const fn kind(&self) -> InventoryDiffKind {
        self.kind
    }

    /// Semantic row kind.
    #[must_use]
    pub fn row_kind(&self) -> &str {
        &self.row_kind
    }

    /// Stable semantic subject.
    #[must_use]
    pub fn subject(&self) -> &str {
        &self.subject
    }

    /// Stable semantic field.
    #[must_use]
    pub fn field(&self) -> &str {
        &self.field
    }

    /// Canonical values in the older inventory.
    #[must_use]
    pub fn before(&self) -> &[String] {
        &self.before
    }

    /// Canonical values in the newer inventory.
    #[must_use]
    pub fn after(&self) -> &[String] {
        &self.after
    }
}

/// Sorted semantic diff between two immutable inventories.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InventoryDiff {
    from: InventoryDigest,
    to: InventoryDigest,
    entries: Vec<InventoryDiffEntry>,
}

impl InventoryDiff {
    /// Older inventory digest.
    #[must_use]
    pub const fn from(&self) -> InventoryDigest {
        self.from
    }

    /// Newer inventory digest.
    #[must_use]
    pub const fn to(&self) -> InventoryDigest {
        self.to
    }

    /// Sorted reviewable entries.
    #[must_use]
    pub fn entries(&self) -> &[InventoryDiffEntry] {
        &self.entries
    }

    /// Whether the two inventories have identical semantic projections.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Deterministic escaped human rendering.
    #[must_use]
    pub fn render_human(&self) -> String {
        let mut out = String::new();
        for entry in &self.entries {
            writeln!(
                out,
                "from={} to={} change={} kind={} subject={} field={} before={} after={}",
                self.from,
                self.to,
                entry.kind.as_str(),
                escape_record_text(&entry.row_kind),
                escape_record_text(&entry.subject),
                escape_record_text(&entry.field),
                escape_record_text(&encode_values(&entry.before)),
                escape_record_text(&encode_values(&entry.after)),
            )
            .expect("writing to String is infallible");
        }
        out
    }
}

/// Sealed canonical inventory. All projections derive from this one value.
#[derive(Debug, Clone)]
pub struct FrozenInventory {
    limits: InventoryLimits,
    sources: Vec<InventorySource>,
    source_set_digest: InventorySourceSetDigest,
    revisions: Vec<(ClaimRevisionId, ClaimRevision)>,
    graph: NormalizedGraph,
    facts: Vec<InventoryFact>,
    observations: Vec<InventoryObservation>,
    resolutions: Vec<InventoryFieldResolution>,
    reconciliations: Vec<ReconciliationReceipt>,
    conflicts: Vec<InventoryConflict>,
    digest: InventoryDigest,
    authority_policy_version: u32,
    reconciliation_policy_version: u32,
}

impl PartialEq for FrozenInventory {
    fn eq(&self, other: &Self) -> bool {
        self.digest == other.digest
    }
}

impl Eq for FrozenInventory {}

impl FrozenInventory {
    /// Complete canonical inventory digest.
    #[must_use]
    pub const fn digest(&self) -> InventoryDigest {
        self.digest
    }

    /// Exact canonical input-source-set digest.
    #[must_use]
    pub const fn source_set_digest(&self) -> InventorySourceSetDigest {
        self.source_set_digest
    }

    /// Admitted source snapshots in canonical order.
    #[must_use]
    pub fn sources(&self) -> &[InventorySource] {
        &self.sources
    }

    /// Every full immutable claim revision in revision-id order.
    #[must_use]
    pub fn revisions(&self) -> &[(ClaimRevisionId, ClaimRevision)] {
        &self.revisions
    }

    /// V.1.1 normalized scientific relation graph.
    #[must_use]
    pub const fn graph(&self) -> &NormalizedGraph {
        &self.graph
    }

    /// Every retained source fact, including losing/conflicting assertions.
    #[must_use]
    pub fn facts(&self) -> &[InventoryFact] {
        &self.facts
    }

    /// Every admitted evidence observation, including negative, incomplete,
    /// unavailable, and refuted history.
    #[must_use]
    pub fn observations(&self) -> &[InventoryObservation] {
        &self.observations
    }

    /// Complete 22-field resolutions for every revision.
    #[must_use]
    pub fn resolutions(&self) -> &[InventoryFieldResolution] {
        &self.resolutions
    }

    /// Explicit alias/rename/split/merge receipts.
    #[must_use]
    pub fn reconciliations(&self) -> &[ReconciliationReceipt] {
        &self.reconciliations
    }

    /// Canonical retained conflict set.
    #[must_use]
    pub fn conflicts(&self) -> &[InventoryConflict] {
        &self.conflicts
    }

    /// Whether every record is complete and disagreement-free.
    #[must_use]
    pub fn conflict_free(&self) -> bool {
        self.conflicts.is_empty()
    }

    /// Whether any conflict prevents a complete authoritative record.
    #[must_use]
    pub fn has_blocking_conflicts(&self) -> bool {
        self.conflicts
            .iter()
            .any(|conflict| conflict.kind.blocking())
    }

    /// Authority policy version retained by the artifact.
    #[must_use]
    pub const fn authority_policy_version(&self) -> u32 {
        self.authority_policy_version
    }

    /// Reconciliation policy version retained by the artifact.
    #[must_use]
    pub const fn reconciliation_policy_version(&self) -> u32 {
        self.reconciliation_policy_version
    }

    /// Effective compiler limits retained for diagnostics. Like the v1 graph
    /// envelope, tightened ceilings are admission metadata and do not change
    /// canonical identity after the same semantics successfully admit.
    #[must_use]
    pub const fn limits(&self) -> InventoryLimits {
        self.limits
    }

    /// Create a metadata-only frozen-replay source draft. Its canonical
    /// locator and snapshot both bind this sealed inventory. Frozen replay
    /// sources cannot author facts, observations, or reconciliation; callers
    /// replay the sealed inventory itself rather than attaching new privileged
    /// assertions to its identity.
    #[must_use]
    pub fn replay_source_draft(
        &self,
        adapter_version: impl Into<String>,
        adapter_policy_version: u32,
    ) -> InventorySourceDraft {
        InventorySourceDraft::normalized(
            InventorySourceKind::FrozenInventory,
            InventoryRole::FrozenReplay,
            format!("inventory:{}", self.digest),
            self.digest.content_hash(),
            adapter_version,
            adapter_policy_version,
        )
    }
}

/// Compile caller-normalized heterogeneous sources into one sealed inventory.
///
/// Structural invalidity and resource exhaustion return `Err` without a
/// partial artifact. Semantically valid disagreement returns `Ok` with an
/// explicit, digest-bound conflict set.
pub fn compile_inventory(
    draft: &InventoryDraft,
    limits: InventoryLimits,
) -> Result<FrozenInventory, InventoryRefusal> {
    validate_limits(limits)?;
    preflight_counts(draft, limits)?;
    preflight_semantic_bytes(draft, limits)?;
    preflight_projection_floor(draft, limits)?;
    if draft.authority_policy_version != INVENTORY_AUTHORITY_POLICY_VERSION {
        return Err(InventoryRefusal::new(
            "inventory-authority-policy-version",
            format!(
                "authority policy version {} is unsupported; compiler implements {}",
                draft.authority_policy_version, INVENTORY_AUTHORITY_POLICY_VERSION
            ),
        )
        .with_fix("recompile the normalized inputs under the implemented authority policy"));
    }
    if draft.reconciliation_policy_version != INVENTORY_RECONCILIATION_POLICY_VERSION {
        return Err(InventoryRefusal::new(
            "inventory-reconciliation-policy-version",
            format!(
                "reconciliation policy version {} is unsupported; compiler implements {}",
                draft.reconciliation_policy_version, INVENTORY_RECONCILIATION_POLICY_VERSION
            ),
        )
        .with_fix("re-admit reconciliation receipts under the implemented policy"));
    }
    if draft.sources.is_empty() {
        return Err(InventoryRefusal::new(
            "inventory-empty-sources",
            "a verification inventory must bind at least one exact source snapshot",
        )
        .with_fix("supply caller-normalized, content-addressed source records"));
    }
    if draft.revisions.is_empty() {
        return Err(InventoryRefusal::new(
            "inventory-empty-revisions",
            "an empty claim inventory is non-green and cannot be sealed",
        )
        .with_fix("supply every live and historical ClaimRevision in scope"));
    }

    let mut semantic_bytes = derived_resolution_bytes(draft.sources.len(), draft.revisions.len())?;
    let (sources_by_input, mut sources) = admit_sources(&draft.sources, &mut semantic_bytes)?;
    enforce_semantic_budget(semantic_bytes, limits)?;
    let source_set_digest = inventory_source_set_digest(&sources)?;

    let mut revisions = admit_revisions(&draft.revisions)?;
    let graph = admit_graph(&draft.revisions, &draft.relations).map_err(inventory_graph_refusal)?;

    let mut facts = admit_facts(
        &draft.facts,
        &sources_by_input,
        &revisions,
        limits,
        &mut semantic_bytes,
    )?;
    let mut observations = admit_observations(
        &draft.observations,
        &sources_by_input,
        &revisions,
        &mut semantic_bytes,
    )?;
    let mut reconciliations = admit_reconciliations(
        &draft.reconciliations,
        &sources_by_input,
        &revisions,
        draft.reconciliation_policy_version,
        &mut semantic_bytes,
    )?;
    enforce_semantic_budget(semantic_bytes, limits)?;

    sources.sort_by_key(InventorySource::id);
    revisions.sort_by_key(|(id, _)| *id);
    facts.sort();
    observations.sort();
    reconciliations.sort();

    let mut conflicts = source_snapshot_conflicts(&sources);
    conflicts.extend(observation_conflicts(&observations));
    conflicts.extend(reconciliation_conflicts(&reconciliations));
    let resolutions = resolve_fields(&sources, &revisions, &facts, &mut conflicts)?;
    conflicts.sort();
    conflicts.dedup();

    let conflict_bytes = conflicts.iter().try_fold(0u64, |total, conflict| {
        total
            .checked_add(checked_usize_to_u64(
                conflict.detail.len(),
                "conflict semantic bytes",
                "bytes",
            )?)
            .ok_or_else(|| inventory_overflow("conflict semantic bytes", "bytes"))
    })?;
    semantic_bytes = semantic_bytes
        .checked_add(conflict_bytes)
        .ok_or_else(|| inventory_overflow("aggregate inventory semantic bytes", "bytes"))?;
    enforce_semantic_budget(semantic_bytes, limits)?;

    let projection_rows = projection_row_count(
        sources.len(),
        revisions.len(),
        graph.edges().len(),
        graph.representatives().len(),
        facts.len(),
        observations.len(),
        resolutions.len(),
        reconciliations.len(),
        conflicts.len(),
    )?;
    if projection_rows > u64::from(limits.max_projection_rows) {
        return Err(inventory_resource_refusal(
            "projection row count",
            projection_rows,
            u64::from(limits.max_projection_rows),
            "rows",
        ));
    }

    let digest = inventory_digest(
        draft.authority_policy_version,
        draft.reconciliation_policy_version,
        source_set_digest,
        &revisions,
        &graph,
        &facts,
        &observations,
        &resolutions,
        &reconciliations,
        &conflicts,
    )?;

    Ok(FrozenInventory {
        limits,
        sources,
        source_set_digest,
        revisions,
        graph,
        facts,
        observations,
        resolutions,
        reconciliations,
        conflicts,
        digest,
        authority_policy_version: draft.authority_policy_version,
        reconciliation_policy_version: draft.reconciliation_policy_version,
    })
}

fn inventory_graph_refusal(error: V1Error) -> InventoryRefusal {
    let mut refusal = InventoryRefusal::new(
        "inventory-claim-graph",
        format!(
            "v1 graph admission refused at {}: {}",
            error.rule(),
            error.detail()
        ),
    );
    refusal.fixes.extend(error.ranked_fixes().iter().cloned());
    if let Some(resource) = error.resource_refusal() {
        refusal.resource = Some(InventoryResourceRefusal {
            quantity: resource.quantity,
            required: resource.required,
            admitted: resource.admitted,
            unit: resource.unit,
        });
    }
    refusal
}

fn validate_limits(limits: InventoryLimits) -> Result<(), InventoryRefusal> {
    for (quantity, requested, maximum, unit) in [
        (
            "configured source ceiling",
            u64::from(limits.max_sources),
            u64::from(InventoryLimits::DEFAULT.max_sources),
            "sources",
        ),
        (
            "configured revision ceiling",
            u64::from(limits.max_revisions),
            u64::from(InventoryLimits::DEFAULT.max_revisions),
            "revisions",
        ),
        (
            "configured relation ceiling",
            u64::from(limits.max_relations),
            u64::from(InventoryLimits::DEFAULT.max_relations),
            "relations",
        ),
        (
            "configured fact ceiling",
            u64::from(limits.max_facts),
            u64::from(InventoryLimits::DEFAULT.max_facts),
            "facts",
        ),
        (
            "configured observation ceiling",
            u64::from(limits.max_observations),
            u64::from(InventoryLimits::DEFAULT.max_observations),
            "observations",
        ),
        (
            "configured reconciliation ceiling",
            u64::from(limits.max_reconciliations),
            u64::from(InventoryLimits::DEFAULT.max_reconciliations),
            "receipts",
        ),
        (
            "configured reconciliation-endpoint ceiling",
            u64::from(limits.max_reconciliation_endpoints),
            u64::from(InventoryLimits::DEFAULT.max_reconciliation_endpoints),
            "endpoints",
        ),
        (
            "configured fact-value ceiling",
            u64::from(limits.max_values_per_fact),
            u64::from(InventoryLimits::DEFAULT.max_values_per_fact),
            "values",
        ),
        (
            "configured semantic-byte ceiling",
            limits.max_semantic_bytes,
            InventoryLimits::DEFAULT.max_semantic_bytes,
            "bytes",
        ),
        (
            "configured projection-row ceiling",
            u64::from(limits.max_projection_rows),
            u64::from(InventoryLimits::DEFAULT.max_projection_rows),
            "rows",
        ),
    ] {
        if requested > maximum {
            return Err(inventory_resource_refusal(
                quantity, requested, maximum, unit,
            ));
        }
    }
    Ok(())
}

fn preflight_counts(
    draft: &InventoryDraft,
    limits: InventoryLimits,
) -> Result<(), InventoryRefusal> {
    let reconciliation_endpoints =
        draft
            .reconciliations
            .iter()
            .try_fold(0u64, |total, reconciliation| {
                total
                    .checked_add(reconciliation.endpoint_count()?)
                    .ok_or_else(|| inventory_overflow("reconciliation endpoint count", "endpoints"))
            })?;
    for (quantity, required, admitted, unit) in [
        (
            "source count",
            checked_usize_to_u64(draft.sources.len(), "source count", "sources")?,
            u64::from(limits.max_sources),
            "sources",
        ),
        (
            "revision count",
            checked_usize_to_u64(draft.revisions.len(), "revision count", "revisions")?,
            u64::from(limits.max_revisions),
            "revisions",
        ),
        (
            "relation count",
            checked_usize_to_u64(draft.relations.len(), "relation count", "relations")?,
            u64::from(limits.max_relations),
            "relations",
        ),
        (
            "fact count",
            checked_usize_to_u64(draft.facts.len(), "fact count", "facts")?,
            u64::from(limits.max_facts),
            "facts",
        ),
        (
            "observation count",
            checked_usize_to_u64(
                draft.observations.len(),
                "observation count",
                "observations",
            )?,
            u64::from(limits.max_observations),
            "observations",
        ),
        (
            "reconciliation count",
            checked_usize_to_u64(
                draft.reconciliations.len(),
                "reconciliation count",
                "receipts",
            )?,
            u64::from(limits.max_reconciliations),
            "receipts",
        ),
        (
            "reconciliation endpoint count",
            reconciliation_endpoints,
            u64::from(limits.max_reconciliation_endpoints),
            "endpoints",
        ),
    ] {
        if required > admitted {
            return Err(inventory_resource_refusal(
                quantity, required, admitted, unit,
            ));
        }
    }
    for fact in &draft.facts {
        let required =
            checked_usize_to_u64(fact.values.len(), "values in one field assertion", "values")?;
        if required > u64::from(limits.max_values_per_fact) {
            return Err(inventory_resource_refusal(
                "values in one field assertion",
                required,
                u64::from(limits.max_values_per_fact),
                "values",
            ));
        }
    }
    Ok(())
}

fn preflight_semantic_bytes(
    draft: &InventoryDraft,
    limits: InventoryLimits,
) -> Result<(), InventoryRefusal> {
    let mut required = derived_resolution_bytes(draft.sources.len(), draft.revisions.len())?;
    let mut charge = |value: &str| -> Result<(), InventoryRefusal> {
        required = required
            .checked_add(checked_usize_to_u64(
                value.len(),
                "aggregate inventory semantic bytes",
                "bytes",
            )?)
            .ok_or_else(|| inventory_overflow("aggregate inventory semantic bytes", "bytes"))?;
        Ok(())
    };
    for source in &draft.sources {
        charge(&source.pin.source)?;
        charge(&source.adapter_version)?;
    }
    for fact in &draft.facts {
        for value in &fact.values {
            charge(value)?;
        }
    }
    for observation in &draft.observations {
        charge(&observation.observation_id)?;
    }
    for reconciliation in &draft.reconciliations {
        charge(reconciliation.rationale())?;
        match reconciliation {
            ReconciliationDraft::Alias {
                alias, canonical, ..
            } => {
                charge(alias.as_str())?;
                charge(canonical.as_str())?;
            }
            ReconciliationDraft::Rename {
                previous, current, ..
            } => {
                charge(previous.as_str())?;
                charge(current.as_str())?;
            }
            ReconciliationDraft::Split { .. } | ReconciliationDraft::Merge { .. } => {}
        }
    }
    if required > limits.max_semantic_bytes {
        return Err(inventory_resource_refusal(
            "aggregate inventory semantic bytes",
            required,
            limits.max_semantic_bytes,
            "bytes",
        ));
    }
    Ok(())
}

fn preflight_projection_floor(
    draft: &InventoryDraft,
    limits: InventoryLimits,
) -> Result<(), InventoryRefusal> {
    let resolutions = draft
        .revisions
        .len()
        .checked_mul(InventoryField::ALL.len())
        .ok_or_else(|| inventory_overflow("projection row count", "rows"))?;
    let floor = [
        1usize,
        draft.sources.len(),
        draft.revisions.len(),
        draft.relations.len(),
        draft.revisions.len(),
        draft.facts.len(),
        draft.observations.len(),
        resolutions,
        draft.reconciliations.len(),
    ]
    .into_iter()
    .try_fold(0u64, |total, count| {
        total
            .checked_add(checked_usize_to_u64(count, "projection row count", "rows")?)
            .ok_or_else(|| inventory_overflow("projection row count", "rows"))
    })?;
    if floor > u64::from(limits.max_projection_rows) {
        return Err(inventory_resource_refusal(
            "projection row count",
            floor,
            u64::from(limits.max_projection_rows),
            "rows",
        ));
    }
    Ok(())
}

fn derived_resolution_bytes(sources: usize, revisions: usize) -> Result<u64, InventoryRefusal> {
    let sources = checked_usize_to_u64(sources, "derived source resolution", "sources")?;
    let revisions = checked_usize_to_u64(revisions, "derived source resolution", "revisions")?;
    let source_entries = sources
        .checked_mul(revisions)
        .and_then(|entries| entries.checked_mul(DERIVED_SOURCE_RESOLUTION_BYTES))
        .ok_or_else(|| inventory_overflow("derived source resolution bytes", "bytes"))?;
    let revision_entries = revisions
        .checked_mul(DERIVED_REVISION_RESOLUTION_BYTES)
        .ok_or_else(|| inventory_overflow("derived revision resolution bytes", "bytes"))?;
    source_entries
        .checked_add(revision_entries)
        .ok_or_else(|| inventory_overflow("derived resolution bytes", "bytes"))
}

fn inventory_resource_refusal(
    quantity: &'static str,
    required: u64,
    admitted: u64,
    unit: &'static str,
) -> InventoryRefusal {
    InventoryRefusal::new(
        "inventory-resource-limit",
        format!("{quantity} requires {required} {unit}, admitted {admitted} {unit}"),
    )
    .with_resource(InventoryResourceRefusal {
        quantity,
        required,
        admitted,
        unit,
    })
    .with_fix(format!("reduce {quantity} to at most {admitted} {unit}"))
}

fn inventory_overflow(quantity: &'static str, unit: &'static str) -> InventoryRefusal {
    InventoryRefusal::new(
        "inventory-resource-overflow",
        format!("{quantity} cannot be represented in u64 {unit}"),
    )
    .with_fix("split the inventory; wrapped resource accounting is forbidden")
}

fn checked_usize_to_u64(
    value: usize,
    quantity: &'static str,
    unit: &'static str,
) -> Result<u64, InventoryRefusal> {
    u64::try_from(value).map_err(|_| inventory_overflow(quantity, unit))
}

fn checked_collection_with_anchor(length: usize) -> Result<u64, InventoryRefusal> {
    checked_usize_to_u64(length, "reconciliation endpoint count", "endpoints")?
        .checked_add(1)
        .ok_or_else(|| inventory_overflow("reconciliation endpoint count", "endpoints"))
}

fn inventory_allocation_refusal(payload: &'static str, count: usize) -> InventoryRefusal {
    InventoryRefusal::new(
        "inventory-allocation-refused",
        format!("allocator refused {payload} reservation for {count} entries"),
    )
    .with_fix("release memory pressure and retry the same bounded input")
}

fn enforce_semantic_budget(required: u64, limits: InventoryLimits) -> Result<(), InventoryRefusal> {
    if required > limits.max_semantic_bytes {
        Err(inventory_resource_refusal(
            "aggregate inventory semantic bytes",
            required,
            limits.max_semantic_bytes,
            "bytes",
        ))
    } else {
        Ok(())
    }
}

fn charge_text(total: &mut u64, field: &'static str, value: &str) -> Result<(), InventoryRefusal> {
    if value.is_empty() || value.len() > MAX_INVENTORY_TEXT_BYTES {
        return Err(InventoryRefusal::new(
            "inventory-text-bounds",
            format!(
                "{field} length {} outside 1..={MAX_INVENTORY_TEXT_BYTES}",
                value.len()
            ),
        )
        .with_fix(format!("supply a non-empty bounded {field}")));
    }
    *total = total
        .checked_add(checked_usize_to_u64(
            value.len(),
            "aggregate inventory semantic bytes",
            "bytes",
        )?)
        .ok_or_else(|| inventory_overflow("aggregate inventory semantic bytes", "bytes"))?;
    Ok(())
}

fn allowed_role(kind: InventorySourceKind, role: InventoryRole) -> bool {
    matches!(
        (kind, role),
        (InventorySourceKind::Beads, InventoryRole::Obligation)
            | (
                InventorySourceKind::Contract,
                InventoryRole::DeclaredSemantics
            )
            | (
                InventorySourceKind::TypedRegistry | InventorySourceKind::CodeRegistration,
                InventoryRole::ExecutableRegistration
            )
            | (
                InventorySourceKind::TestRegistration,
                InventoryRole::ValidationContext
            )
            | (
                InventorySourceKind::VvArtifact
                    | InventorySourceKind::BenchmarkRegistry
                    | InventorySourceKind::LedgerReceipt,
                InventoryRole::ObservedEvidence
            )
            | (
                InventorySourceKind::FrozenInventory,
                InventoryRole::FrozenReplay
            )
    )
}

fn admit_sources(
    drafts: &[InventorySourceDraft],
    semantic_bytes: &mut u64,
) -> Result<(Vec<InventorySource>, Vec<InventorySource>), InventoryRefusal> {
    let mut by_input = Vec::new();
    by_input
        .try_reserve_exact(drafts.len())
        .map_err(|_| inventory_allocation_refusal("source input map", drafts.len()))?;
    let mut seen = BTreeSet::new();
    for source in drafts {
        if source.pin.authority != source.kind.authority() {
            return Err(InventoryRefusal::new(
                "inventory-source-authority",
                format!(
                    "source {:?} presents authority {}, but kind {} fixes authority {}",
                    source.pin.source,
                    source_authority_name(source.pin.authority),
                    source.kind.as_str(),
                    source_authority_name(source.kind.authority())
                ),
            )
            .with_fix("retain the source's real class; content addressing is not authority"));
        }
        if !allowed_role(source.kind, source.role) {
            return Err(InventoryRefusal::new(
                "inventory-source-role",
                format!(
                    "source kind {} cannot occupy semantic role {}",
                    source.kind.as_str(),
                    source.role.as_str()
                ),
            )
            .with_fix("use the closed kind/role matrix; do not relabel observed evidence"));
        }
        charge_text(semantic_bytes, "source locator", &source.pin.source)?;
        charge_text(semantic_bytes, "adapter version", &source.adapter_version)?;
        if source.adapter_policy_version == 0 {
            return Err(InventoryRefusal::new(
                "inventory-adapter-policy-version",
                format!(
                    "source {:?} has zero adapter policy version",
                    source.pin.source
                ),
            )
            .with_fix("pin a nonzero normalization-policy version"));
        }
        if source.pin.snapshot == ContentHash([0; 32]) {
            return Err(InventoryRefusal::new(
                "inventory-source-snapshot",
                format!(
                    "source {:?} uses the all-zero missing hash",
                    source.pin.source
                ),
            )
            .with_fix("hash the exact admitted source bytes before compilation"));
        }
        let id = inventory_source_id(source)?;
        if !seen.insert(id) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-source",
                format!("normalized source {id} was supplied more than once"),
            )
            .with_fix("deduplicate inside the source adapter before compilation"));
        }
        by_input.push(InventorySource {
            id,
            kind: source.kind,
            role: source.role,
            pin: source.pin.clone(),
            adapter_version: source.adapter_version.clone(),
            adapter_policy_version: source.adapter_policy_version,
        });
    }
    Ok((by_input.clone(), by_input))
}

fn inventory_source_id(
    source: &InventorySourceDraft,
) -> Result<InventorySourceId, InventoryRefusal> {
    let mut hasher = DomainHasher::new(SOURCE_ID_DOMAIN);
    hasher.update(&INVENTORY_SCHEMA_VERSION.to_be_bytes());
    hasher.update(&[source.kind.tag(), source.role.tag()]);
    hash_text(&mut hasher, &source.pin.source)?;
    hasher.update(&[authority_tag(source.pin.authority)]);
    hasher.update(source.pin.snapshot.as_bytes());
    hash_text(&mut hasher, &source.adapter_version)?;
    hasher.update(&source.adapter_policy_version.to_be_bytes());
    Ok(InventorySourceId(hasher.finalize()))
}

fn inventory_source_set_digest(
    sources: &[InventorySource],
) -> Result<InventorySourceSetDigest, InventoryRefusal> {
    let mut ids: Vec<InventorySourceId> = sources.iter().map(InventorySource::id).collect();
    ids.sort_unstable();
    let mut hasher = DomainHasher::new(SOURCE_SET_DOMAIN);
    hasher.update(&INVENTORY_SCHEMA_VERSION.to_be_bytes());
    hash_count(&mut hasher, ids.len())?;
    for id in ids {
        hasher.update(id.0.as_bytes());
    }
    Ok(InventorySourceSetDigest(hasher.finalize()))
}

fn admit_revisions(
    drafts: &[ClaimRevision],
) -> Result<Vec<(ClaimRevisionId, ClaimRevision)>, InventoryRefusal> {
    let mut revisions = Vec::new();
    revisions
        .try_reserve_exact(drafts.len())
        .map_err(|_| inventory_allocation_refusal("claim revisions", drafts.len()))?;
    let mut ids = BTreeSet::new();
    for revision in drafts {
        let id = revision.revision_id().map_err(inventory_graph_refusal)?;
        if !ids.insert(id) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-revision",
                format!("claim revision {} was supplied more than once", id.to_hex()),
            )
            .with_fix("deduplicate normalized revision registrations"));
        }
        revisions.push((id, revision.clone()));
    }
    let by_id: BTreeMap<ClaimRevisionId, &ClaimRevision> = revisions
        .iter()
        .map(|(id, revision)| (*id, revision))
        .collect();
    for (id, revision) in &revisions {
        if let Some(predecessor) = revision.supersedes {
            if predecessor == *id {
                return Err(InventoryRefusal::new(
                    "inventory-self-supersession",
                    format!("revision {} supersedes itself", id.to_hex()),
                ));
            }
            if !by_id.contains_key(&predecessor) {
                return Err(InventoryRefusal::new(
                    "inventory-dangling-supersession",
                    format!(
                        "revision {} supersedes unavailable revision {}",
                        id.to_hex(),
                        predecessor.to_hex()
                    ),
                )
                .with_fix(
                    "include the exact historical predecessor; never bind to current metadata",
                ));
            }
        }
    }
    Ok(revisions)
}

fn admit_facts(
    drafts: &[InventoryFactDraft],
    sources: &[InventorySource],
    revisions: &[(ClaimRevisionId, ClaimRevision)],
    limits: InventoryLimits,
    semantic_bytes: &mut u64,
) -> Result<Vec<InventoryFact>, InventoryRefusal> {
    let known_revisions: BTreeSet<ClaimRevisionId> = revisions.iter().map(|(id, _)| *id).collect();
    let mut facts = Vec::new();
    facts
        .try_reserve_exact(drafts.len())
        .map_err(|_| inventory_allocation_refusal("inventory facts", drafts.len()))?;
    let mut seen_slots = BTreeSet::new();
    for fact in drafts {
        if fact.field.is_derived() {
            return Err(InventoryRefusal::new(
                "inventory-derived-field",
                format!(
                    "field {} is compiler-derived and cannot be source-authored",
                    fact.field.as_str()
                ),
            )
            .with_fix(
                "remove the assertion; the compiler binds the exact source/revision identities",
            ));
        }
        if !known_revisions.contains(&fact.revision) {
            return Err(InventoryRefusal::new(
                "inventory-orphan-fact",
                format!(
                    "field {} targets unavailable revision {}",
                    fact.field.as_str(),
                    fact.revision.to_hex()
                ),
            )
            .with_fix("include the exact revision or remove the dangling assertion"));
        }
        let source = sources.get(fact.source_index).ok_or_else(|| {
            InventoryRefusal::new(
                "inventory-source-index",
                format!(
                    "field {} source index {} outside 0..{}",
                    fact.field.as_str(),
                    fact.source_index,
                    sources.len()
                ),
            )
        })?;
        if source.role == InventoryRole::FrozenReplay {
            return Err(InventoryRefusal::new(
                "inventory-frozen-replay-payload",
                format!(
                    "frozen inventory source {} cannot author field {}",
                    source.id,
                    fact.field.as_str()
                ),
            )
            .with_fix("use the sealed inventory directly; never attach new facts to its digest"));
        }
        if !seen_slots.insert((fact.revision, fact.field, source.id)) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-fact",
                format!(
                    "source {} repeats field {} for revision {}",
                    source.id,
                    fact.field.as_str(),
                    fact.revision.to_hex()
                ),
            )
            .with_fix("emit one complete assertion per source, revision, and field"));
        }
        let value_count =
            checked_usize_to_u64(fact.values.len(), "values in one field assertion", "values")?;
        if value_count > u64::from(limits.max_values_per_fact) {
            return Err(inventory_resource_refusal(
                "values in one field assertion",
                value_count,
                u64::from(limits.max_values_per_fact),
                "values",
            ));
        }
        validate_cardinality(fact.field, fact.values.len())?;
        let mut values = fact.values.clone();
        for value in &values {
            charge_text(semantic_bytes, "field value", value)?;
            validate_field_value(fact.field, value)?;
        }
        values.sort();
        if values.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-field-value",
                format!(
                    "source {} repeats a {} value for revision {}",
                    source.id,
                    fact.field.as_str(),
                    fact.revision.to_hex()
                ),
            )
            .with_fix(
                "deduplicate within the source adapter; corroboration requires distinct sources",
            ));
        }
        facts.push(InventoryFact {
            revision: fact.revision,
            field: fact.field,
            values,
            source: source.id,
            role: source.role,
            authority: source.pin.authority,
        });
    }
    Ok(facts)
}

fn validate_cardinality(field: InventoryField, count: usize) -> Result<(), InventoryRefusal> {
    let valid = match field.spec().cardinality {
        "1" => count == 1,
        "0..1" => count <= 1,
        "1..n" => count >= 1,
        "0..n" => true,
        other => {
            return Err(InventoryRefusal::new(
                "inventory-registry-cardinality",
                format!(
                    "field {} has unsupported cardinality {other}",
                    field.as_str()
                ),
            ));
        }
    };
    if valid {
        Ok(())
    } else {
        Err(InventoryRefusal::new(
            "inventory-field-cardinality",
            format!(
                "field {} requires cardinality {}, got {count}",
                field.as_str(),
                field.spec().cardinality
            ),
        ))
    }
}

fn admit_observations(
    drafts: &[InventoryObservationDraft],
    sources: &[InventorySource],
    revisions: &[(ClaimRevisionId, ClaimRevision)],
    semantic_bytes: &mut u64,
) -> Result<Vec<InventoryObservation>, InventoryRefusal> {
    let known_revisions: BTreeSet<ClaimRevisionId> = revisions.iter().map(|(id, _)| *id).collect();
    let mut observations = Vec::new();
    observations
        .try_reserve_exact(drafts.len())
        .map_err(|_| inventory_allocation_refusal("inventory observations", drafts.len()))?;
    let mut seen_slots = BTreeSet::new();
    for draft in drafts {
        if !known_revisions.contains(&draft.revision) {
            return Err(InventoryRefusal::new(
                "inventory-orphan-observation",
                format!(
                    "observation {:?} targets unavailable revision {}",
                    draft.observation_id,
                    draft.revision.to_hex()
                ),
            )
            .with_fix("retain the exact historical revision named by the observation"));
        }
        let source = sources.get(draft.source_index).ok_or_else(|| {
            InventoryRefusal::new(
                "inventory-source-index",
                format!(
                    "observation {:?} source index {} outside 0..{}",
                    draft.observation_id,
                    draft.source_index,
                    sources.len()
                ),
            )
        })?;
        if !matches!(source.role, InventoryRole::ObservedEvidence) {
            return Err(InventoryRefusal::new(
                "inventory-observation-source-role",
                format!(
                    "observation {:?} source {} has role {}, not observed-evidence",
                    draft.observation_id,
                    source.id,
                    source.role.as_str()
                ),
            )
            .with_fix("bind results to a V&V, benchmark, or ledger observation source"));
        }
        charge_text(semantic_bytes, "observation id", &draft.observation_id)?;
        validate_inventory_id("observation", &draft.observation_id)?;
        if draft.artifact_digest == ContentHash([0; 32]) {
            return Err(InventoryRefusal::new(
                "inventory-observation-artifact",
                format!(
                    "observation {:?} uses the all-zero missing artifact digest",
                    draft.observation_id
                ),
            )
            .with_fix("retain the exact evidence or failure-bundle artifact digest"));
        }
        if draft.adjudication.requires_complete_evidence()
            && (draft.execution != ObservationExecution::Completed
                || draft.completeness != ObservationCompleteness::Complete
                || draft.integrity != ObservationIntegrity::Verified)
        {
            return Err(InventoryRefusal::new(
                "inventory-observation-axes",
                format!(
                    "observation {:?} reports adjudication {} with execution {}, completeness {}, integrity {}",
                    draft.observation_id,
                    draft.adjudication.as_str(),
                    draft.execution.as_str(),
                    draft.completeness.as_str(),
                    draft.integrity.as_str()
                ),
            )
            .with_fix("use Unknown or Pending unless execution completed with complete verified evidence"));
        }
        if !seen_slots.insert((draft.observation_id.clone(), source.id)) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-observation",
                format!(
                    "source {} repeats observation {:?} for revision {}",
                    source.id,
                    draft.observation_id,
                    draft.revision.to_hex()
                ),
            )
            .with_fix("emit one immutable observation per source-local id and revision"));
        }
        observations.push(InventoryObservation {
            revision: draft.revision,
            observation_id: draft.observation_id.clone(),
            artifact_digest: draft.artifact_digest,
            execution: draft.execution,
            completeness: draft.completeness,
            integrity: draft.integrity,
            adjudication: draft.adjudication,
            source: source.id,
        });
    }
    Ok(observations)
}

fn validate_field_value(field: InventoryField, value: &str) -> Result<(), InventoryRefusal> {
    match field {
        InventoryField::Stratum if !matches!(value, "core" | "max") => {
            return Err(InventoryRefusal::new(
                "inventory-stratum",
                format!("stratum {value:?} is outside core|max"),
            ));
        }
        InventoryField::Ambition if !matches!(value, "S" | "F" | "M") => {
            return Err(InventoryRefusal::new(
                "inventory-ambition",
                format!("ambition {value:?} is outside S|F|M"),
            ));
        }
        InventoryField::CaseIds => {
            CaseId::new(value).map_err(inventory_graph_refusal)?;
        }
        InventoryField::JourneyIds => {
            JourneyId::new(value).map_err(inventory_graph_refusal)?;
        }
        InventoryField::BeadObligation
        | InventoryField::CampaignProfiles
        | InventoryField::FixtureIds
        | InventoryField::OracleIds
        | InventoryField::CheckerIds
        | InventoryField::Capabilities
        | InventoryField::EventKinds => validate_generic_id(field, value)?,
        InventoryField::SourceSnapshots | InventoryField::ClaimRevision => {
            return Err(InventoryRefusal::new(
                "inventory-derived-field",
                format!("field {} is compiler-derived", field.as_str()),
            ));
        }
        InventoryField::Stratum
        | InventoryField::Ambition
        | InventoryField::PublicSurface
        | InventoryField::Ownership
        | InventoryField::TcbOverlap
        | InventoryField::ToleranceDerivation
        | InventoryField::Budgets
        | InventoryField::Retention
        | InventoryField::ReplayCommand
        | InventoryField::DsrLane
        | InventoryField::ReceiptExpectations => {}
    }
    Ok(())
}

fn validate_generic_id(field: InventoryField, value: &str) -> Result<(), InventoryRefusal> {
    validate_inventory_id(field.as_str(), value)
}

fn validate_inventory_id(kind: &str, value: &str) -> Result<(), InventoryRefusal> {
    if value.len() > 128
        || value.is_empty()
        || !value.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || b"-_.:/".contains(&byte)
        })
    {
        return Err(InventoryRefusal::new(
            "inventory-identity-bounds",
            format!(
                "{kind} value {value:?} is outside the nonempty canonical [a-z0-9-_.:/] identity grammar"
            ),
        ));
    }
    Ok(())
}

fn admit_reconciliations(
    drafts: &[ReconciliationDraft],
    sources: &[InventorySource],
    revisions: &[(ClaimRevisionId, ClaimRevision)],
    expected_policy_version: u32,
    semantic_bytes: &mut u64,
) -> Result<Vec<ReconciliationReceipt>, InventoryRefusal> {
    let known_claims: BTreeSet<&str> = revisions
        .iter()
        .map(|(_, revision)| revision.claim.as_str())
        .collect();
    let revisions_by_id: BTreeMap<ClaimRevisionId, &ClaimRevision> = revisions
        .iter()
        .map(|(id, revision)| (*id, revision))
        .collect();
    let mut receipts = Vec::new();
    receipts
        .try_reserve_exact(drafts.len())
        .map_err(|_| inventory_allocation_refusal("reconciliation receipts", drafts.len()))?;
    for draft in drafts {
        let kind = draft.kind();
        let source = sources.get(draft.source_index()).ok_or_else(|| {
            InventoryRefusal::new(
                "inventory-source-index",
                format!(
                    "{} reconciliation source index {} outside 0..{}",
                    kind.as_str(),
                    draft.source_index(),
                    sources.len()
                ),
            )
        })?;
        if !matches!(source.role, InventoryRole::Obligation) {
            return Err(InventoryRefusal::new(
                "inventory-reconciliation-source-role",
                format!(
                    "{} reconciliation source {} has role {}; only obligation sources may declare lineage topology",
                    kind.as_str(),
                    source.id,
                    source.role.as_str()
                ),
            )
            .with_fix(
                "move the declaration to its owning obligation; code, observations, and prior inventory pointers cannot rewrite lineage",
            ));
        }
        if draft.policy_version() == 0 || draft.policy_version() != expected_policy_version {
            return Err(InventoryRefusal::new(
                "inventory-reconciliation-policy-version",
                format!(
                    "{} receipt policy {} does not match compiled policy {}",
                    kind.as_str(),
                    draft.policy_version(),
                    expected_policy_version
                ),
            ));
        }
        charge_text(
            semantic_bytes,
            "reconciliation rationale",
            draft.rationale(),
        )?;
        let topology = match draft {
            ReconciliationDraft::Alias {
                alias, canonical, ..
            } => {
                charge_text(semantic_bytes, "alias lineage endpoint", alias.as_str())?;
                charge_text(semantic_bytes, "alias lineage endpoint", canonical.as_str())?;
                validate_lineage_reconciliation(
                    ReconciliationKind::Alias,
                    alias,
                    canonical,
                    &known_claims,
                )?;
                ReconciliationTopology::Alias {
                    alias: alias.clone(),
                    canonical: canonical.clone(),
                }
            }
            ReconciliationDraft::Rename {
                previous, current, ..
            } => {
                charge_text(semantic_bytes, "rename lineage endpoint", previous.as_str())?;
                charge_text(semantic_bytes, "rename lineage endpoint", current.as_str())?;
                validate_lineage_reconciliation(
                    ReconciliationKind::Rename,
                    previous,
                    current,
                    &known_claims,
                )?;
                ReconciliationTopology::Rename {
                    previous: previous.clone(),
                    current: current.clone(),
                }
            }
            ReconciliationDraft::Split {
                predecessor,
                successors,
                ..
            } => {
                let successors = validate_split(*predecessor, successors, &revisions_by_id)?;
                ReconciliationTopology::Split {
                    predecessor: *predecessor,
                    successors,
                }
            }
            ReconciliationDraft::Merge {
                predecessors,
                successor,
                ..
            } => {
                let predecessors = validate_merge(predecessors, *successor, &revisions_by_id)?;
                ReconciliationTopology::Merge {
                    predecessors,
                    successor: *successor,
                }
            }
        };
        receipts.push(ReconciliationReceipt {
            topology,
            source: source.id,
            policy_version: draft.policy_version(),
            rationale: draft.rationale().to_owned(),
        });
    }
    receipts.sort();
    let mut seen_topology = BTreeSet::new();
    for receipt in &receipts {
        if !seen_topology.insert((
            receipt.topology.clone(),
            receipt.source,
            receipt.policy_version,
        )) {
            return Err(InventoryRefusal::new(
                "inventory-duplicate-reconciliation-topology",
                format!(
                    "source {} repeats {} topology under policy {} with a second rationale",
                    receipt.source,
                    receipt.kind().as_str(),
                    receipt.policy_version
                ),
            )
            .with_fix("emit one complete topology assertion per source and policy"));
        }
    }
    validate_reconciliation_topology(&receipts)?;
    Ok(receipts)
}

fn validate_lineage_reconciliation(
    kind: ReconciliationKind,
    from: &ClaimId,
    to: &ClaimId,
    known_claims: &BTreeSet<&str>,
) -> Result<(), InventoryRefusal> {
    if from == to {
        return Err(InventoryRefusal::new(
            "inventory-reconciliation-self-edge",
            format!(
                "{} reconciliation maps {} to itself",
                kind.as_str(),
                from.as_str()
            ),
        ));
    }
    for endpoint in [from, to] {
        if !known_claims.contains(endpoint.as_str()) {
            return Err(InventoryRefusal::new(
                "inventory-reconciliation-orphan",
                format!(
                    "{} lineage endpoint {:?} has no retained ClaimRevision",
                    kind.as_str(),
                    endpoint.as_str()
                ),
            )
            .with_fix("retain every historical lineage named by reconciliation"));
        }
    }
    Ok(())
}

fn validate_split(
    predecessor: ClaimRevisionId,
    successors: &[ClaimRevisionId],
    revisions: &BTreeMap<ClaimRevisionId, &ClaimRevision>,
) -> Result<Vec<ClaimRevisionId>, InventoryRefusal> {
    require_revision(predecessor, "split predecessor", revisions)?;
    let successors = canonical_revision_endpoints("split successors", successors, 2)?;
    let mut lineages = BTreeSet::new();
    for successor in &successors {
        if *successor == predecessor {
            return Err(InventoryRefusal::new(
                "inventory-reconciliation-self-edge",
                format!(
                    "split successor {} is its own predecessor",
                    successor.to_hex()
                ),
            ));
        }
        let revision = require_revision(*successor, "split successor", revisions)?;
        if revision.supersedes != Some(predecessor) {
            return Err(InventoryRefusal::new(
                "inventory-split-history",
                format!(
                    "split successor {} does not supersede exact predecessor {}",
                    successor.to_hex(),
                    predecessor.to_hex()
                ),
            )
            .with_fix(
                "set every split successor's immutable supersedes anchor to the predecessor",
            ));
        }
        lineages.insert(revision.claim.as_str());
    }
    if lineages.len() != successors.len() {
        return Err(InventoryRefusal::new(
            "inventory-split-lineages",
            "split successors do not occupy distinct claim lineages",
        )
        .with_fix("use ordinary supersession for multiple revisions of one lineage"));
    }
    Ok(successors)
}

fn validate_merge(
    predecessors: &[ClaimRevisionId],
    successor: ClaimRevisionId,
    revisions: &BTreeMap<ClaimRevisionId, &ClaimRevision>,
) -> Result<Vec<ClaimRevisionId>, InventoryRefusal> {
    let predecessors = canonical_revision_endpoints("merge predecessors", predecessors, 2)?;
    if predecessors.binary_search(&successor).is_ok() {
        return Err(InventoryRefusal::new(
            "inventory-reconciliation-self-edge",
            format!(
                "merge successor {} is also a predecessor",
                successor.to_hex()
            ),
        ));
    }
    let successor_revision = require_revision(successor, "merge successor", revisions)?;
    let supersedes = successor_revision.supersedes.ok_or_else(|| {
        InventoryRefusal::new(
            "inventory-merge-history",
            format!(
                "merge successor {} has no supersedes anchor",
                successor.to_hex()
            ),
        )
        .with_fix("anchor the successor to one exact predecessor; the receipt retains the full set")
    })?;
    if predecessors.binary_search(&supersedes).is_err() {
        return Err(InventoryRefusal::new(
            "inventory-merge-history",
            format!(
                "merge successor {} supersedes {}, outside the declared predecessor set",
                successor.to_hex(),
                supersedes.to_hex()
            ),
        )
        .with_fix("include the immutable supersedes anchor in the merge predecessor set"));
    }
    let mut lineages = BTreeSet::new();
    for predecessor in &predecessors {
        let revision = require_revision(*predecessor, "merge predecessor", revisions)?;
        lineages.insert(revision.claim.as_str());
    }
    if lineages.len() != predecessors.len() {
        return Err(InventoryRefusal::new(
            "inventory-merge-lineages",
            "merge predecessors do not occupy distinct claim lineages",
        )
        .with_fix("use ordinary supersession for multiple revisions of one lineage"));
    }
    Ok(predecessors)
}

fn canonical_revision_endpoints(
    label: &'static str,
    endpoints: &[ClaimRevisionId],
    minimum: usize,
) -> Result<Vec<ClaimRevisionId>, InventoryRefusal> {
    if endpoints.len() < minimum {
        return Err(InventoryRefusal::new(
            "inventory-reconciliation-arity",
            format!(
                "{label} has {}, requires at least {minimum}",
                endpoints.len()
            ),
        ));
    }
    let mut canonical = endpoints.to_vec();
    canonical.sort_unstable();
    if canonical.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(InventoryRefusal::new(
            "inventory-reconciliation-duplicate-endpoint",
            format!("{label} contains an exact endpoint more than once"),
        )
        .with_fix("deduplicate the source record before compilation"));
    }
    Ok(canonical)
}

fn require_revision<'a>(
    revision: ClaimRevisionId,
    label: &'static str,
    revisions: &'a BTreeMap<ClaimRevisionId, &ClaimRevision>,
) -> Result<&'a ClaimRevision, InventoryRefusal> {
    revisions.get(&revision).copied().ok_or_else(|| {
        InventoryRefusal::new(
            "inventory-reconciliation-orphan",
            format!("{label} {} is not retained", revision.to_hex()),
        )
        .with_fix("retain the exact historical revision named by reconciliation")
    })
}

fn validate_reconciliation_topology(
    receipts: &[ReconciliationReceipt],
) -> Result<(), InventoryRefusal> {
    let mut presentation_graph: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for receipt in receipts {
        match &receipt.topology {
            ReconciliationTopology::Alias {
                alias, canonical, ..
            } => {
                presentation_graph
                    .entry(alias.as_str())
                    .or_default()
                    .insert(canonical.as_str());
            }
            ReconciliationTopology::Rename {
                previous, current, ..
            } => {
                presentation_graph
                    .entry(previous.as_str())
                    .or_default()
                    .insert(current.as_str());
            }
            ReconciliationTopology::Split { .. } | ReconciliationTopology::Merge { .. } => {}
        }
    }
    ensure_acyclic_presentation_graph(&presentation_graph)
}

fn ensure_acyclic_presentation_graph(
    graph: &BTreeMap<&str, BTreeSet<&str>>,
) -> Result<(), InventoryRefusal> {
    let mut colors: BTreeMap<&str, u8> = BTreeMap::new();
    let mut nodes = BTreeSet::new();
    for (from, targets) in graph {
        nodes.insert(*from);
        nodes.extend(targets.iter().copied());
    }
    for root in nodes {
        if colors.get(root).copied().unwrap_or(0) != 0 {
            continue;
        }
        let mut stack = vec![(root, false)];
        while let Some((node, exiting)) = stack.pop() {
            if exiting {
                colors.insert(node, 2);
                continue;
            }
            match colors.get(node).copied().unwrap_or(0) {
                1 => {
                    return Err(InventoryRefusal::new(
                        "inventory-reconciliation-cycle",
                        format!("alias/rename cycle reaches {node:?}"),
                    )
                    .with_fix(
                        "choose one canonical presentation id and direct aliases toward it",
                    ));
                }
                2 => continue,
                _ => {}
            }
            colors.insert(node, 1);
            stack.push((node, true));
            if let Some(targets) = graph.get(node) {
                for target in targets.iter().rev() {
                    if colors.get(target).copied().unwrap_or(0) == 1 {
                        return Err(InventoryRefusal::new(
                            "inventory-reconciliation-cycle",
                            format!("alias/rename cycle reaches {target:?}"),
                        )
                        .with_fix(
                            "choose one canonical presentation id and direct aliases toward it",
                        ));
                    }
                    stack.push((target, false));
                }
            }
        }
    }
    Ok(())
}

fn source_snapshot_conflicts(sources: &[InventorySource]) -> Vec<InventoryConflict> {
    let mut by_locator: BTreeMap<&str, Vec<&InventorySource>> = BTreeMap::new();
    for source in sources {
        by_locator
            .entry(&source.pin.source)
            .or_default()
            .push(source);
    }
    let mut conflicts = Vec::new();
    for (locator, group) in by_locator {
        let snapshots: BTreeSet<ContentHash> =
            group.iter().map(|source| source.pin.snapshot).collect();
        if snapshots.len() <= 1 {
            continue;
        }
        let maximum = group
            .iter()
            .map(|source| source.pin.authority)
            .max()
            .expect("nonempty source locator group");
        let winning_snapshots: BTreeSet<ContentHash> = group
            .iter()
            .filter(|source| source.pin.authority == maximum)
            .map(|source| source.pin.snapshot)
            .collect();
        let mut source_ids: Vec<InventorySourceId> = group.iter().map(|source| source.id).collect();
        source_ids.sort_unstable();
        let (kind, detail) = if winning_snapshots.len() == 1 {
            let winner = winning_snapshots
                .iter()
                .next()
                .expect("one winning snapshot");
            (
                InventoryConflictKind::SourceSnapshotResolved,
                format!(
                    "source {locator:?} has multiple snapshots; authority {} explicitly selects {} while retaining all pins",
                    source_authority_name(maximum),
                    winner.to_hex()
                ),
            )
        } else {
            (
                InventoryConflictKind::SourceSnapshotConflict,
                format!(
                    "source {locator:?} has {} distinct snapshots at equal maximum authority {}",
                    winning_snapshots.len(),
                    source_authority_name(maximum)
                ),
            )
        };
        conflicts.push(InventoryConflict {
            kind,
            revision: None,
            field: Some(InventoryField::SourceSnapshots),
            sources: source_ids,
            detail,
        });
    }
    conflicts
}

fn observation_conflicts(observations: &[InventoryObservation]) -> Vec<InventoryConflict> {
    let mut by_revision: BTreeMap<
        ClaimRevisionId,
        BTreeMap<ObservationAdjudication, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    for observation in observations
        .iter()
        .filter(|observation| observation.adjudication.requires_complete_evidence())
    {
        by_revision
            .entry(observation.revision)
            .or_default()
            .entry(observation.adjudication)
            .or_default()
            .push(observation.source);
    }
    let mut conflicts = Vec::new();
    for (revision, adjudications) in by_revision {
        if adjudications.len() <= 1 {
            continue;
        }
        let mut sources: Vec<InventorySourceId> =
            adjudications.values().flatten().copied().collect();
        sources.sort_unstable();
        sources.dedup();
        let states: Vec<String> = adjudications
            .keys()
            .map(|adjudication| adjudication.as_str().to_owned())
            .collect();
        conflicts.push(InventoryConflict {
            kind: InventoryConflictKind::ObservationAdjudicationConflict,
            revision: Some(revision),
            field: None,
            sources,
            detail: format!(
                "revision {} has incompatible complete verified adjudications {}",
                revision.to_hex(),
                encode_values(&states)
            ),
        });
    }
    conflicts
}

fn reconciliation_conflicts(receipts: &[ReconciliationReceipt]) -> Vec<InventoryConflict> {
    let mut presentation: BTreeMap<
        &str,
        BTreeMap<(ReconciliationKind, &str), Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut splits: BTreeMap<
        ClaimRevisionId,
        BTreeMap<Vec<ClaimRevisionId>, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut splits_by_successors: BTreeMap<
        Vec<ClaimRevisionId>,
        BTreeMap<ClaimRevisionId, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut merges: BTreeMap<
        ClaimRevisionId,
        BTreeMap<Vec<ClaimRevisionId>, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut merges_by_predecessors: BTreeMap<
        Vec<ClaimRevisionId>,
        BTreeMap<ClaimRevisionId, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut merges_by_predecessor: BTreeMap<
        ClaimRevisionId,
        BTreeMap<usize, Vec<InventorySourceId>>,
    > = BTreeMap::new();
    let mut merge_topology_ordinals: BTreeMap<ReconciliationTopology, usize> = BTreeMap::new();
    for receipt in receipts {
        match &receipt.topology {
            ReconciliationTopology::Alias { alias, canonical } => {
                presentation
                    .entry(alias.as_str())
                    .or_default()
                    .entry((ReconciliationKind::Alias, canonical.as_str()))
                    .or_default()
                    .push(receipt.source);
            }
            ReconciliationTopology::Rename { previous, current } => {
                presentation
                    .entry(previous.as_str())
                    .or_default()
                    .entry((ReconciliationKind::Rename, current.as_str()))
                    .or_default()
                    .push(receipt.source);
            }
            ReconciliationTopology::Split {
                predecessor,
                successors,
            } => {
                splits
                    .entry(*predecessor)
                    .or_default()
                    .entry(successors.clone())
                    .or_default()
                    .push(receipt.source);
                splits_by_successors
                    .entry(successors.clone())
                    .or_default()
                    .entry(*predecessor)
                    .or_default()
                    .push(receipt.source);
            }
            ReconciliationTopology::Merge {
                predecessors,
                successor,
            } => {
                let next_ordinal = merge_topology_ordinals.len();
                let topology_ordinal = match merge_topology_ordinals.get(&receipt.topology) {
                    Some(ordinal) => *ordinal,
                    None => {
                        merge_topology_ordinals.insert(receipt.topology.clone(), next_ordinal);
                        next_ordinal
                    }
                };
                merges
                    .entry(*successor)
                    .or_default()
                    .entry(predecessors.clone())
                    .or_default()
                    .push(receipt.source);
                merges_by_predecessors
                    .entry(predecessors.clone())
                    .or_default()
                    .entry(*successor)
                    .or_default()
                    .push(receipt.source);
                for predecessor in predecessors {
                    merges_by_predecessor
                        .entry(*predecessor)
                        .or_default()
                        .entry(topology_ordinal)
                        .or_default()
                        .push(receipt.source);
                }
            }
        }
    }
    let mut conflicts = Vec::new();
    for (from, targets) in presentation {
        if targets.len() <= 1 {
            continue;
        }
        let mut sources: Vec<InventorySourceId> = targets.values().flatten().copied().collect();
        sources.sort_unstable();
        sources.dedup();
        let mappings: Vec<String> = targets
            .keys()
            .map(|(kind, target)| format!("{}:{target}", kind.as_str()))
            .collect();
        conflicts.push(InventoryConflict {
            kind: InventoryConflictKind::ReconciliationConflict,
            revision: None,
            field: None,
            sources,
            detail: format!(
                "presentation source {} maps through multiple kind/target pairs {}",
                from,
                encode_values(&mappings)
            ),
        });
    }
    for (predecessor, successor_sets) in splits {
        if successor_sets.len() <= 1 {
            continue;
        }
        conflicts.push(reconciliation_set_conflict(
            format!(
                "split predecessor {} maps to multiple exact successor sets {}",
                predecessor.to_hex(),
                encode_revision_sets(successor_sets.keys())
            ),
            successor_sets.values(),
        ));
    }
    for (successors, predecessors) in splits_by_successors {
        if predecessors.len() <= 1 {
            continue;
        }
        let predecessor_ids: Vec<String> = predecessors.keys().map(ContentHash::to_hex).collect();
        conflicts.push(reconciliation_set_conflict(
            format!(
                "split successor set {} maps from multiple predecessors {}",
                encode_values(
                    &successors
                        .iter()
                        .map(ContentHash::to_hex)
                        .collect::<Vec<_>>()
                ),
                encode_values(&predecessor_ids)
            ),
            predecessors.values(),
        ));
    }
    for (successor, predecessor_sets) in merges {
        if predecessor_sets.len() <= 1 {
            continue;
        }
        conflicts.push(reconciliation_set_conflict(
            format!(
                "merge successor {} maps from multiple exact predecessor sets {}",
                successor.to_hex(),
                encode_revision_sets(predecessor_sets.keys())
            ),
            predecessor_sets.values(),
        ));
    }
    for (predecessors, successors) in merges_by_predecessors {
        if successors.len() <= 1 {
            continue;
        }
        let successor_ids: Vec<String> = successors.keys().map(ContentHash::to_hex).collect();
        conflicts.push(reconciliation_set_conflict(
            format!(
                "merge predecessor set {} maps into multiple successors {}",
                encode_values(
                    &predecessors
                        .iter()
                        .map(ContentHash::to_hex)
                        .collect::<Vec<_>>()
                ),
                encode_values(&successor_ids)
            ),
            successors.values(),
        ));
    }
    for (predecessor, topologies) in merges_by_predecessor {
        if topologies.len() <= 1 {
            continue;
        }
        let ordinals: Vec<String> = topologies.keys().map(ToString::to_string).collect();
        conflicts.push(reconciliation_set_conflict(
            format!(
                "merge predecessor {} participates in {} exact merge topologies at canonical merge ordinals {}",
                predecessor.to_hex(),
                topologies.len(),
                encode_values(&ordinals)
            ),
            topologies.values(),
        ));
    }
    conflicts
}

fn reconciliation_set_conflict<'a>(
    detail: String,
    source_sets: impl Iterator<Item = &'a Vec<InventorySourceId>>,
) -> InventoryConflict {
    let mut sources: Vec<InventorySourceId> = source_sets.flatten().copied().collect();
    sources.sort_unstable();
    sources.dedup();
    InventoryConflict {
        kind: InventoryConflictKind::ReconciliationConflict,
        revision: None,
        field: None,
        sources,
        detail,
    }
}

fn encode_revision_sets<'a>(sets: impl Iterator<Item = &'a Vec<ClaimRevisionId>>) -> String {
    let encoded: Vec<String> = sets
        .map(|set| {
            let revisions: Vec<String> = set.iter().map(ContentHash::to_hex).collect();
            encode_values(&revisions)
        })
        .collect();
    encode_values(&encoded)
}

fn resolve_fields(
    sources: &[InventorySource],
    revisions: &[(ClaimRevisionId, ClaimRevision)],
    facts: &[InventoryFact],
    conflicts: &mut Vec<InventoryConflict>,
) -> Result<Vec<InventoryFieldResolution>, InventoryRefusal> {
    let source_ids: Vec<String> = sources.iter().map(|source| source.id.to_hex()).collect();
    let mut facts_by_slot: BTreeMap<(ClaimRevisionId, InventoryField), Vec<&InventoryFact>> =
        BTreeMap::new();
    for fact in facts {
        facts_by_slot
            .entry((fact.revision, fact.field))
            .or_default()
            .push(fact);
    }
    let resolution_capacity = revisions
        .len()
        .checked_mul(InventoryField::ALL.len())
        .ok_or_else(|| inventory_overflow("field resolution count", "resolutions"))?;
    let mut resolutions = Vec::new();
    resolutions
        .try_reserve_exact(resolution_capacity)
        .map_err(|_| inventory_allocation_refusal("field resolutions", resolution_capacity))?;
    for (revision_id, _) in revisions {
        for field in InventoryField::ALL {
            if field == InventoryField::SourceSnapshots {
                resolutions.push(InventoryFieldResolution {
                    revision: *revision_id,
                    field,
                    values: Some(source_ids.clone()),
                    sources: sources.iter().map(InventorySource::id).collect(),
                });
                continue;
            }
            if field == InventoryField::ClaimRevision {
                resolutions.push(InventoryFieldResolution {
                    revision: *revision_id,
                    field,
                    values: Some(vec![revision_id.to_hex()]),
                    sources: Vec::new(),
                });
                continue;
            }
            let candidates = facts_by_slot
                .get(&(*revision_id, field))
                .map(Vec::as_slice)
                .unwrap_or(&[]);
            let expected_authority = field.spec().authority;
            let expected_role = role_for_authority(expected_authority);
            let mut eligible = Vec::new();
            let mut cross_role = Vec::new();
            for candidate in candidates {
                let declared_source =
                    candidate.authority == expected_authority && candidate.role == expected_role;
                if declared_source {
                    eligible.push(*candidate);
                } else {
                    cross_role.push(*candidate);
                }
            }
            if eligible.is_empty() {
                if field_required(field) {
                    let mut provenance: Vec<InventorySourceId> =
                        candidates.iter().map(|fact| fact.source).collect();
                    provenance.sort_unstable();
                    provenance.dedup();
                    conflicts.push(InventoryConflict {
                        kind: InventoryConflictKind::MissingRequiredField,
                        revision: Some(*revision_id),
                        field: Some(field),
                        sources: provenance,
                        detail: format!(
                            "revision {} has no authoritative {} assertion",
                            revision_id.to_hex(),
                            field.as_str()
                        ),
                    });
                    resolutions.push(InventoryFieldResolution {
                        revision: *revision_id,
                        field,
                        values: None,
                        sources: Vec::new(),
                    });
                } else {
                    let mut disagreements: Vec<InventorySourceId> = cross_role
                        .iter()
                        .filter(|candidate| !candidate.values.is_empty())
                        .map(|candidate| candidate.source)
                        .collect();
                    disagreements.sort_unstable();
                    disagreements.dedup();
                    if !disagreements.is_empty() {
                        conflicts.push(InventoryConflict {
                            kind: InventoryConflictKind::CrossRoleDisagreement,
                            revision: Some(*revision_id),
                            field: Some(field),
                            sources: disagreements,
                            detail: format!(
                                "revision {} optional field {} is absent from authoritative role {}, but one or more non-owning semantic roles report nonempty values",
                                revision_id.to_hex(),
                                field.as_str(),
                                expected_role.as_str()
                            ),
                        });
                    }
                    resolutions.push(InventoryFieldResolution {
                        revision: *revision_id,
                        field,
                        values: Some(Vec::new()),
                        sources: Vec::new(),
                    });
                }
                continue;
            }
            let mut by_value: BTreeMap<Vec<String>, Vec<InventorySourceId>> = BTreeMap::new();
            for candidate in eligible {
                by_value
                    .entry(candidate.values.clone())
                    .or_default()
                    .push(candidate.source);
            }
            if by_value.len() == 1 {
                let (values, mut provenance) = by_value
                    .into_iter()
                    .next()
                    .expect("one authoritative field value");
                provenance.sort_unstable();
                provenance.dedup();
                let mut disagreements: Vec<InventorySourceId> = cross_role
                    .iter()
                    .filter(|candidate| candidate.values.as_slice() != values.as_slice())
                    .map(|candidate| candidate.source)
                    .collect();
                if !disagreements.is_empty() {
                    disagreements.extend(provenance.iter().copied());
                    disagreements.sort_unstable();
                    disagreements.dedup();
                    conflicts.push(InventoryConflict {
                        kind: InventoryConflictKind::CrossRoleDisagreement,
                        revision: Some(*revision_id),
                        field: Some(field),
                        sources: disagreements,
                        detail: format!(
                            "revision {} field {} differs between authoritative role {} and one or more non-owning semantic roles",
                            revision_id.to_hex(),
                            field.as_str(),
                            expected_role.as_str()
                        ),
                    });
                }
                resolutions.push(InventoryFieldResolution {
                    revision: *revision_id,
                    field,
                    values: Some(values),
                    sources: provenance,
                });
            } else {
                let mut provenance: Vec<InventorySourceId> =
                    by_value.values().flatten().copied().collect();
                provenance.sort_unstable();
                provenance.dedup();
                conflicts.push(InventoryConflict {
                    kind: InventoryConflictKind::EqualAuthorityDisagreement,
                    revision: Some(*revision_id),
                    field: Some(field),
                    sources: provenance.clone(),
                    detail: format!(
                        "revision {} has {} distinct {} assertions at declared authority {}",
                        revision_id.to_hex(),
                        by_value.len(),
                        field.as_str(),
                        source_authority_name(expected_authority)
                    ),
                });
                resolutions.push(InventoryFieldResolution {
                    revision: *revision_id,
                    field,
                    values: None,
                    sources: provenance,
                });
            }
        }
    }
    Ok(resolutions)
}

fn field_required(field: InventoryField) -> bool {
    matches!(field.spec().cardinality, "1" | "1..n")
}

fn role_for_authority(authority: SourceAuthority) -> InventoryRole {
    match authority {
        SourceAuthority::GeneratedArtifact => InventoryRole::ExecutableRegistration,
        SourceAuthority::TestSource => InventoryRole::ValidationContext,
        SourceAuthority::Contract => InventoryRole::DeclaredSemantics,
        SourceAuthority::BeadObligation => InventoryRole::Obligation,
        SourceAuthority::FrozenSnapshot => InventoryRole::FrozenReplay,
    }
}

fn projection_row_count(
    sources: usize,
    revisions: usize,
    relations: usize,
    representatives: usize,
    facts: usize,
    observations: usize,
    resolutions: usize,
    reconciliations: usize,
    conflicts: usize,
) -> Result<u64, InventoryRefusal> {
    [
        1usize,
        sources,
        revisions,
        relations,
        representatives,
        facts,
        observations,
        resolutions,
        reconciliations,
        conflicts,
    ]
    .into_iter()
    .try_fold(0u64, |total, count| {
        total
            .checked_add(checked_usize_to_u64(count, "projection row count", "rows")?)
            .ok_or_else(|| inventory_overflow("projection row count", "rows"))
    })
}

fn inventory_digest(
    authority_policy_version: u32,
    reconciliation_policy_version: u32,
    source_set_digest: InventorySourceSetDigest,
    revisions: &[(ClaimRevisionId, ClaimRevision)],
    graph: &NormalizedGraph,
    facts: &[InventoryFact],
    observations: &[InventoryObservation],
    resolutions: &[InventoryFieldResolution],
    reconciliations: &[ReconciliationReceipt],
    conflicts: &[InventoryConflict],
) -> Result<InventoryDigest, InventoryRefusal> {
    let mut hasher = DomainHasher::new(INVENTORY_ID_DOMAIN);
    hasher.update(&INVENTORY_SCHEMA_VERSION.to_be_bytes());
    hasher.update(&INVENTORY_COMPILER_VERSION.to_be_bytes());
    hasher.update(&authority_policy_version.to_be_bytes());
    hasher.update(&reconciliation_policy_version.to_be_bytes());
    hasher.update(source_set_digest.0.as_bytes());
    hasher.update(graph.digest().as_bytes());

    hash_count(&mut hasher, revisions.len())?;
    for (id, revision) in revisions {
        hasher.update(id.as_bytes());
        hash_text(&mut hasher, revision.claim.as_str())?;
        match revision.supersedes {
            Some(predecessor) => {
                hasher.update(&[1]);
                hasher.update(predecessor.as_bytes());
            }
            None => hasher.update(&[0]),
        }
    }

    hash_count(&mut hasher, facts.len())?;
    for fact in facts {
        hasher.update(fact.revision.as_bytes());
        hasher.update(&[
            field_tag(fact.field),
            fact.role.tag(),
            authority_tag(fact.authority),
        ]);
        hasher.update(fact.source.0.as_bytes());
        hash_values(&mut hasher, &fact.values)?;
    }

    hash_count(&mut hasher, observations.len())?;
    for observation in observations {
        hasher.update(observation.revision.as_bytes());
        hash_text(&mut hasher, &observation.observation_id)?;
        hasher.update(observation.artifact_digest.as_bytes());
        hasher.update(&[
            observation.execution.tag(),
            observation.completeness.tag(),
            observation.integrity.tag(),
            observation.adjudication.tag(),
        ]);
        hasher.update(observation.source.0.as_bytes());
    }

    hash_count(&mut hasher, resolutions.len())?;
    for resolution in resolutions {
        hasher.update(resolution.revision.as_bytes());
        hasher.update(&[field_tag(resolution.field)]);
        match &resolution.values {
            Some(values) => {
                hasher.update(&[1]);
                hash_values(&mut hasher, values)?;
            }
            None => hasher.update(&[0]),
        }
        hash_count(&mut hasher, resolution.sources.len())?;
        for source in &resolution.sources {
            hasher.update(source.0.as_bytes());
        }
    }

    hash_count(&mut hasher, reconciliations.len())?;
    for receipt in reconciliations {
        hasher.update(&[receipt.kind().tag()]);
        match &receipt.topology {
            ReconciliationTopology::Alias {
                alias, canonical, ..
            } => {
                hash_text(&mut hasher, alias.as_str())?;
                hash_text(&mut hasher, canonical.as_str())?;
            }
            ReconciliationTopology::Rename {
                previous, current, ..
            } => {
                hash_text(&mut hasher, previous.as_str())?;
                hash_text(&mut hasher, current.as_str())?;
            }
            ReconciliationTopology::Split {
                predecessor,
                successors,
                ..
            } => {
                hasher.update(predecessor.as_bytes());
                hash_count(&mut hasher, successors.len())?;
                for successor in successors {
                    hasher.update(successor.as_bytes());
                }
            }
            ReconciliationTopology::Merge {
                predecessors,
                successor,
                ..
            } => {
                hash_count(&mut hasher, predecessors.len())?;
                for predecessor in predecessors {
                    hasher.update(predecessor.as_bytes());
                }
                hasher.update(successor.as_bytes());
            }
        }
        hasher.update(receipt.source().0.as_bytes());
        hasher.update(&receipt.policy_version().to_be_bytes());
        hash_text(&mut hasher, receipt.rationale())?;
    }

    hash_count(&mut hasher, conflicts.len())?;
    for conflict in conflicts {
        hasher.update(&[conflict.kind.tag()]);
        match conflict.revision {
            Some(revision) => {
                hasher.update(&[1]);
                hasher.update(revision.as_bytes());
            }
            None => hasher.update(&[0]),
        }
        match conflict.field {
            Some(field) => hasher.update(&[1, field_tag(field)]),
            None => hasher.update(&[0]),
        }
        hash_count(&mut hasher, conflict.sources.len())?;
        for source in &conflict.sources {
            hasher.update(source.0.as_bytes());
        }
        hash_text(&mut hasher, &conflict.detail)?;
    }
    Ok(InventoryDigest(hasher.finalize()))
}

fn hash_count(hasher: &mut DomainHasher, count: usize) -> Result<(), InventoryRefusal> {
    let count = u64::try_from(count)
        .map_err(|_| inventory_overflow("canonical collection count", "entries"))?;
    hasher.update(&count.to_be_bytes());
    Ok(())
}

fn hash_text(hasher: &mut DomainHasher, value: &str) -> Result<(), InventoryRefusal> {
    let length = u64::try_from(value.len())
        .map_err(|_| inventory_overflow("canonical text length", "bytes"))?;
    hasher.update(&length.to_be_bytes());
    hasher.update(value.as_bytes());
    Ok(())
}

fn hash_values(hasher: &mut DomainHasher, values: &[String]) -> Result<(), InventoryRefusal> {
    hash_count(hasher, values.len())?;
    for value in values {
        hash_text(hasher, value)?;
    }
    Ok(())
}

const fn authority_tag(authority: SourceAuthority) -> u8 {
    match authority {
        SourceAuthority::GeneratedArtifact => 1,
        SourceAuthority::TestSource => 2,
        SourceAuthority::Contract => 3,
        SourceAuthority::BeadObligation => 4,
        SourceAuthority::FrozenSnapshot => 5,
    }
}

const fn source_authority_name(authority: SourceAuthority) -> &'static str {
    match authority {
        SourceAuthority::GeneratedArtifact => "generated-artifact",
        SourceAuthority::TestSource => "test-source",
        SourceAuthority::Contract => "contract",
        SourceAuthority::BeadObligation => "bead-obligation",
        SourceAuthority::FrozenSnapshot => "frozen-snapshot",
    }
}

const fn field_tag(field: InventoryField) -> u8 {
    match field {
        InventoryField::SourceSnapshots => 1,
        InventoryField::BeadObligation => 2,
        InventoryField::ClaimRevision => 3,
        InventoryField::Stratum => 4,
        InventoryField::CampaignProfiles => 5,
        InventoryField::Ambition => 6,
        InventoryField::PublicSurface => 7,
        InventoryField::CaseIds => 8,
        InventoryField::JourneyIds => 9,
        InventoryField::Ownership => 10,
        InventoryField::FixtureIds => 11,
        InventoryField::OracleIds => 12,
        InventoryField::CheckerIds => 13,
        InventoryField::TcbOverlap => 14,
        InventoryField::ToleranceDerivation => 15,
        InventoryField::Budgets => 16,
        InventoryField::Capabilities => 17,
        InventoryField::EventKinds => 18,
        InventoryField::Retention => 19,
        InventoryField::ReplayCommand => 20,
        InventoryField::DsrLane => 21,
        InventoryField::ReceiptExpectations => 22,
    }
}

fn claim_kind_name(kind: crate::v1::ClaimKind) -> &'static str {
    match kind {
        crate::v1::ClaimKind::Behavioral => "behavioral",
        crate::v1::ClaimKind::QuantitativeBound => "quantitative-bound",
        crate::v1::ClaimKind::Determinism => "determinism",
        crate::v1::ClaimKind::Refusal => "refusal",
        crate::v1::ClaimKind::Theorem => "theorem",
    }
}

fn relation_kind_name(kind: RelationKind) -> &'static str {
    match kind {
        RelationKind::Implication => "implication",
        RelationKind::Refinement => "refinement",
        RelationKind::Restriction => "restriction",
        RelationKind::Counterexample => "counterexample",
        RelationKind::CertifiedEquivalence => "certified-equivalence",
    }
}

fn variance_name(variance: crate::v1::QuantifierVariance) -> &'static str {
    match variance {
        crate::v1::QuantifierVariance::Preserved => "preserved",
        crate::v1::QuantifierVariance::Weakened => "weakened",
        crate::v1::QuantifierVariance::Strengthened => "strengthened",
    }
}

fn encode_values(values: &[String]) -> String {
    let mut out = format!("{}#", values.len());
    for value in values {
        write!(out, "{}:{value}", value.len()).expect("writing to String is infallible");
    }
    out
}

fn encode_optional_values(values: Option<&[String]>) -> String {
    match values {
        Some(values) => format!("resolved:{}", encode_values(values)),
        None => "unresolved".to_owned(),
    }
}

fn push_semantic_row(
    rows: &mut Vec<InventorySemanticRow>,
    kind: &str,
    subject: impl Into<String>,
    field: impl Into<String>,
    value: impl Into<String>,
) -> Result<(), InventoryRefusal> {
    let ordinal = u32::try_from(rows.len())
        .map_err(|_| inventory_overflow("semantic row ordinal", "u32 rows"))?;
    rows.push(InventorySemanticRow {
        ordinal,
        kind: kind.to_owned(),
        subject: subject.into(),
        field: field.into(),
        value: value.into(),
    });
    Ok(())
}

impl FrozenInventory {
    /// Canonical semantic rows shared by human, JSON-lines, and ledger views.
    /// The physical renderings are deliberately not byte-identical.
    pub fn semantic_rows(&self) -> Result<Vec<InventorySemanticRow>, InventoryRefusal> {
        let count = projection_row_count(
            self.sources.len(),
            self.revisions.len(),
            self.graph.edges().len(),
            self.graph.representatives().len(),
            self.facts.len(),
            self.observations.len(),
            self.resolutions.len(),
            self.reconciliations.len(),
            self.conflicts.len(),
        )?;
        let count = usize::try_from(count)
            .map_err(|_| inventory_overflow("semantic row count", "usize rows"))?;
        let mut rows = Vec::new();
        rows.try_reserve_exact(count)
            .map_err(|_| inventory_allocation_refusal("semantic projection rows", count))?;
        push_semantic_row(
            &mut rows,
            "header",
            self.digest.to_string(),
            "inventory",
            format!(
                "inventory_schema={};manifest_schema={};compiler={};authority_policy={};reconciliation_policy={};source_set={};graph={};observations={};conflicts={}",
                INVENTORY_SCHEMA_VERSION,
                MANIFEST_V1_SCHEMA_VERSION,
                INVENTORY_COMPILER_VERSION,
                self.authority_policy_version,
                self.reconciliation_policy_version,
                self.source_set_digest,
                self.graph.digest(),
                self.observations.len(),
                self.conflicts.len()
            ),
        )?;
        for source in &self.sources {
            push_semantic_row(
                &mut rows,
                "source",
                source.id.to_string(),
                source.kind.as_str(),
                encode_values(&[
                    source.role.as_str().to_owned(),
                    source.pin.source.clone(),
                    source_authority_name(source.pin.authority).to_owned(),
                    source.pin.snapshot.to_hex(),
                    source.adapter_version.clone(),
                    source.adapter_policy_version.to_string(),
                ]),
            )?;
        }
        for (revision_id, revision) in &self.revisions {
            push_semantic_row(
                &mut rows,
                "revision",
                revision_id.to_hex(),
                revision.claim.as_str(),
                encode_values(&[
                    claim_kind_name(revision.kind).to_owned(),
                    revision.statement.clone(),
                    revision.quantifiers.clone(),
                    revision.units_conventions.clone(),
                    revision.hypotheses.clone(),
                    revision.domain.clone(),
                    revision.surface.clone(),
                    revision.no_claim.clone(),
                    revision
                        .supersedes
                        .map_or_else(String::new, |id| id.to_hex()),
                ]),
            )?;
        }
        for edge in self.graph.edges() {
            push_semantic_row(
                &mut rows,
                "relation",
                edge.from.to_hex(),
                relation_kind_name(edge.kind),
                encode_values(&[
                    edge.to.to_hex(),
                    edge.checker.clone(),
                    edge.tcb.clone(),
                    variance_name(edge.variance).to_owned(),
                    edge.domain_note.clone(),
                    edge.policy_version.to_string(),
                ]),
            )?;
        }
        for (member, representative) in self.graph.representatives() {
            push_semantic_row(
                &mut rows,
                "equivalence-representative",
                member.to_hex(),
                "representative",
                representative.to_hex(),
            )?;
        }
        for fact in &self.facts {
            push_semantic_row(
                &mut rows,
                "fact",
                fact.revision.to_hex(),
                fact.field.as_str(),
                encode_values(&[
                    fact.source.to_string(),
                    fact.role.as_str().to_owned(),
                    source_authority_name(fact.authority).to_owned(),
                    encode_values(&fact.values),
                ]),
            )?;
        }
        for observation in &self.observations {
            push_semantic_row(
                &mut rows,
                "observation",
                observation.revision.to_hex(),
                &observation.observation_id,
                encode_values(&[
                    observation.artifact_digest.to_hex(),
                    observation.execution.as_str().to_owned(),
                    observation.completeness.as_str().to_owned(),
                    observation.integrity.as_str().to_owned(),
                    observation.adjudication.as_str().to_owned(),
                    observation.source.to_string(),
                ]),
            )?;
        }
        for resolution in &self.resolutions {
            push_semantic_row(
                &mut rows,
                "resolution",
                resolution.revision.to_hex(),
                resolution.field.as_str(),
                encode_values(&[
                    encode_optional_values(resolution.values.as_deref()),
                    encode_values(
                        &resolution
                            .sources
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                    ),
                ]),
            )?;
        }
        for receipt in &self.reconciliations {
            let (subject, endpoints) = match &receipt.topology {
                ReconciliationTopology::Alias {
                    alias, canonical, ..
                } => (alias.as_str().to_owned(), canonical.as_str().to_owned()),
                ReconciliationTopology::Rename {
                    previous, current, ..
                } => (previous.as_str().to_owned(), current.as_str().to_owned()),
                ReconciliationTopology::Split {
                    predecessor,
                    successors,
                    ..
                } => (
                    predecessor.to_hex(),
                    encode_values(
                        &successors
                            .iter()
                            .map(ContentHash::to_hex)
                            .collect::<Vec<_>>(),
                    ),
                ),
                ReconciliationTopology::Merge {
                    predecessors,
                    successor,
                    ..
                } => (
                    successor.to_hex(),
                    encode_values(
                        &predecessors
                            .iter()
                            .map(ContentHash::to_hex)
                            .collect::<Vec<_>>(),
                    ),
                ),
            };
            push_semantic_row(
                &mut rows,
                "reconciliation",
                subject,
                receipt.kind().as_str(),
                encode_values(&[
                    endpoints,
                    receipt.source().to_string(),
                    receipt.policy_version().to_string(),
                    receipt.rationale().to_owned(),
                ]),
            )?;
        }
        for conflict in &self.conflicts {
            push_semantic_row(
                &mut rows,
                "conflict",
                conflict
                    .revision
                    .map_or_else(|| "<inventory>".to_owned(), |id| id.to_hex()),
                conflict.kind.as_str(),
                encode_values(&[
                    conflict
                        .field
                        .map_or_else(String::new, |field| field.as_str().to_owned()),
                    conflict.kind.blocking().to_string(),
                    encode_values(
                        &conflict
                            .sources
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                    ),
                    conflict.detail.clone(),
                ]),
            )?;
        }
        debug_assert_eq!(rows.len(), count);
        Ok(rows)
    }

    /// Deterministic human projection of the canonical semantic rows.
    pub fn render_human(&self) -> Result<String, InventoryRefusal> {
        render_rows(self, ProjectionFormat::Human)
    }

    /// Strict JSON-lines projection of the canonical semantic rows.
    pub fn render_json_lines(&self) -> Result<String, InventoryRefusal> {
        render_rows(self, ProjectionFormat::JsonLines)
    }

    /// Deterministic ledger projection. Rows are compilation-operation
    /// metadata and never a favorable scientific Job result.
    pub fn render_ledger_rows(&self) -> Result<String, InventoryRefusal> {
        render_rows(self, ProjectionFormat::Ledger)
    }

    /// Hash the canonical semantic rows and exact bytes of all three physical
    /// views under distinct domains. All four hashes derive from the same
    /// sealed inventory snapshot.
    pub fn projection_digests(&self) -> Result<InventoryProjectionDigests, InventoryRefusal> {
        let rows = self.semantic_rows()?;
        let semantic = semantic_projection_digest(&rows)?;
        let human = projection_bytes_digest(
            HUMAN_PROJECTION_DOMAIN,
            &render_semantic_rows(self, &rows, ProjectionFormat::Human)?,
        )?;
        let json_lines = projection_bytes_digest(
            JSON_LINES_PROJECTION_DOMAIN,
            &render_semantic_rows(self, &rows, ProjectionFormat::JsonLines)?,
        )?;
        let ledger = projection_bytes_digest(
            LEDGER_PROJECTION_DOMAIN,
            &render_semantic_rows(self, &rows, ProjectionFormat::Ledger)?,
        )?;
        Ok(InventoryProjectionDigests {
            semantic,
            human,
            json_lines,
            ledger,
        })
    }

    /// Compare two immutable inventories by semantic row key. No rename is
    /// inferred: absent/present keys remain remove/add unless an explicit
    /// reconciliation row exists in the inventories themselves.
    pub fn diff(&self, newer: &Self) -> Result<InventoryDiff, InventoryRefusal> {
        let before = semantic_row_map(self.semantic_rows()?);
        let after = semantic_row_map(newer.semantic_rows()?);
        let mut keys = BTreeSet::new();
        keys.extend(before.keys().cloned());
        keys.extend(after.keys().cloned());
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(keys.len())
            .map_err(|_| inventory_allocation_refusal("inventory diff entries", keys.len()))?;
        for (row_kind, subject, field) in keys {
            let before_values = before.get(&(row_kind.clone(), subject.clone(), field.clone()));
            let after_values = after.get(&(row_kind.clone(), subject.clone(), field.clone()));
            let kind = match (before_values, after_values) {
                (Some(before), Some(after)) if before == after => continue,
                (Some(_), Some(_)) => InventoryDiffKind::Changed,
                (Some(_), None) => InventoryDiffKind::Removed,
                (None, Some(_)) => InventoryDiffKind::Added,
                (None, None) => continue,
            };
            entries.push(InventoryDiffEntry {
                kind,
                row_kind,
                subject,
                field,
                before: before_values.cloned().unwrap_or_default(),
                after: after_values.cloned().unwrap_or_default(),
            });
        }
        Ok(InventoryDiff {
            from: self.digest,
            to: newer.digest,
            entries,
        })
    }
}

type SemanticRowKey = (String, String, String);

fn semantic_row_map(rows: Vec<InventorySemanticRow>) -> BTreeMap<SemanticRowKey, Vec<String>> {
    let mut map = BTreeMap::new();
    for row in rows {
        if row.kind == "header" {
            map.entry((row.kind, "<inventory>".to_owned(), row.field))
                .or_insert_with(Vec::new)
                .push(format!("digest={};{}", row.subject, row.value));
        } else {
            map.entry((row.kind, row.subject, row.field))
                .or_insert_with(Vec::new)
                .push(row.value);
        }
    }
    for values in map.values_mut() {
        values.sort();
    }
    map
}

#[derive(Debug, Clone, Copy)]
enum ProjectionFormat {
    Human,
    JsonLines,
    Ledger,
}

fn render_rows(
    inventory: &FrozenInventory,
    format: ProjectionFormat,
) -> Result<String, InventoryRefusal> {
    let rows = inventory.semantic_rows()?;
    render_semantic_rows(inventory, &rows, format)
}

fn render_semantic_rows(
    inventory: &FrozenInventory,
    rows: &[InventorySemanticRow],
    format: ProjectionFormat,
) -> Result<String, InventoryRefusal> {
    let required = rows.iter().try_fold(0usize, |total, row| {
        let semantic = row
            .kind
            .len()
            .checked_add(row.subject.len())
            .and_then(|value| value.checked_add(row.field.len()))
            .and_then(|value| value.checked_add(row.value.len()))
            .and_then(|value| value.checked_mul(6))
            .and_then(|value| value.checked_add(256))
            .ok_or_else(|| inventory_overflow("projection render bytes", "bytes"))?;
        total
            .checked_add(semantic)
            .ok_or_else(|| inventory_overflow("projection render bytes", "bytes"))
    })?;
    let mut out = String::new();
    out.try_reserve_exact(required)
        .map_err(|_| inventory_allocation_refusal("projection bytes", required))?;
    for row in rows {
        match format {
            ProjectionFormat::Human => {
                writeln!(
                    out,
                    "inventory={} ordinal={} kind={} subject={} field={} value={}",
                    inventory.digest,
                    row.ordinal,
                    escape_record_text(&row.kind),
                    escape_record_text(&row.subject),
                    escape_record_text(&row.field),
                    escape_record_text(&row.value),
                )
                .expect("pre-reserved String write is infallible");
            }
            ProjectionFormat::JsonLines => {
                writeln!(
                    out,
                    "{{\"schema_version\":{},\"inventory_digest\":\"{}\",\"ordinal\":{},\"kind\":\"{}\",\"subject\":\"{}\",\"field\":\"{}\",\"value\":\"{}\"}}",
                    INVENTORY_SCHEMA_VERSION,
                    inventory.digest,
                    row.ordinal,
                    escape_json(&row.kind),
                    escape_json(&row.subject),
                    escape_json(&row.field),
                    escape_json(&row.value),
                )
                .expect("pre-reserved String write is infallible");
            }
            ProjectionFormat::Ledger => {
                writeln!(
                    out,
                    "scope=operation outcome=inventory-metadata inventory={} ordinal={} kind={} subject={} field={} value={}",
                    inventory.digest,
                    row.ordinal,
                    escape_ledger_token(&row.kind),
                    escape_ledger_token(&row.subject),
                    escape_ledger_token(&row.field),
                    escape_ledger_token(&row.value),
                )
                .expect("pre-reserved String write is infallible");
            }
        }
    }
    Ok(out)
}

fn semantic_projection_digest(
    rows: &[InventorySemanticRow],
) -> Result<ContentHash, InventoryRefusal> {
    let mut hasher = DomainHasher::new(SEMANTIC_PROJECTION_DOMAIN);
    hasher.update(&INVENTORY_SCHEMA_VERSION.to_be_bytes());
    hash_count(&mut hasher, rows.len())?;
    for row in rows {
        hasher.update(&row.ordinal.to_be_bytes());
        hash_text(&mut hasher, &row.kind)?;
        hash_text(&mut hasher, &row.subject)?;
        hash_text(&mut hasher, &row.field)?;
        hash_text(&mut hasher, &row.value)?;
    }
    Ok(hasher.finalize())
}

fn projection_bytes_digest(
    domain: &'static str,
    bytes: &str,
) -> Result<ContentHash, InventoryRefusal> {
    let mut hasher = DomainHasher::new(domain);
    hasher.update(&INVENTORY_SCHEMA_VERSION.to_be_bytes());
    hash_text(&mut hasher, bytes)?;
    Ok(hasher.finalize())
}

fn escape_record_text(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '\\' => escaped.push_str("\\\\"),
            ' ' => escaped.push_str("\\x20"),
            '=' => escaped.push_str("\\x3d"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{0000}'..='\u{001f}' | '\u{007f}'..='\u{009f}' | '\u{2028}' | '\u{2029}' => {
                write!(escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to String is infallible");
            }
            _ if ch.is_whitespace() => {
                write!(escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to String is infallible");
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn escape_ledger_token(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || b"-_.:/".contains(&byte) {
            escaped.push(char::from(byte));
        } else {
            write!(escaped, "%{byte:02X}").expect("writing to String is infallible");
        }
    }
    escaped
}

fn escape_json(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{0008}' => escaped.push_str("\\b"),
            '\u{000c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{0000}'..='\u{001f}' | '\u{007f}'..='\u{009f}' | '\u{2028}' | '\u{2029}' => {
                write!(escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to String is infallible");
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

/// Replay receipt for one sealed inventory. Matching is exact; callers must
/// not substitute current sources for an unavailable historical source set.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InventoryReceipt {
    /// Schema version.
    pub schema_version: u32,
    /// Compiler version.
    pub compiler_version: u32,
    /// Authority policy version.
    pub authority_policy_version: u32,
    /// Reconciliation policy version.
    pub reconciliation_policy_version: u32,
    /// Exact input source-set digest.
    pub source_set_digest: InventorySourceSetDigest,
    /// Exact V.1.1 graph digest.
    pub graph_digest: ContentHash,
    /// Complete inventory digest.
    pub inventory_digest: InventoryDigest,
}

impl FrozenInventory {
    /// Mint an exact replay receipt for this immutable inventory.
    #[must_use]
    pub fn receipt(&self) -> InventoryReceipt {
        InventoryReceipt {
            schema_version: INVENTORY_SCHEMA_VERSION,
            compiler_version: INVENTORY_COMPILER_VERSION,
            authority_policy_version: self.authority_policy_version,
            reconciliation_policy_version: self.reconciliation_policy_version,
            source_set_digest: self.source_set_digest,
            graph_digest: self.graph.digest(),
            inventory_digest: self.digest,
        }
    }

    /// Require byte-exact receipt identity; no current-source fallback or
    /// alias/equivalence rebinding is attempted.
    pub fn verify_replay_receipt(&self, receipt: InventoryReceipt) -> Result<(), InventoryRefusal> {
        let expected = self.receipt();
        if receipt == expected {
            Ok(())
        } else {
            Err(InventoryRefusal::new(
                "inventory-replay-mismatch",
                format!(
                    "receipt inventory {} / source set {} does not match retained inventory {} / source set {}",
                    receipt.inventory_digest,
                    receipt.source_set_digest,
                    expected.inventory_digest,
                    expected.source_set_digest
                ),
            )
            .with_fix("load the exact archived source set and compiled inventory named by the receipt"))
        }
    }

    /// Verify that every raw source pin needed for exact replay is available.
    /// A current snapshot at the same locator is reported but never
    /// substituted for the historical snapshot.
    pub fn verify_replay_source_availability(
        &self,
        available: &[SourcePin],
    ) -> Result<(), InventoryRefusal> {
        for required in &self.sources {
            if available.iter().any(|candidate| candidate == &required.pin) {
                continue;
            }
            let mut current: Vec<String> = available
                .iter()
                .filter(|candidate| candidate.source == required.pin.source)
                .map(|candidate| {
                    format!(
                        "{}:{}",
                        source_authority_name(candidate.authority),
                        candidate.snapshot.to_hex()
                    )
                })
                .collect();
            current.sort();
            return Err(InventoryRefusal::new(
                "inventory-replay-source-unavailable",
                format!(
                    "replay requires source {:?} at authority {} snapshot {}; available same-locator pins are {}",
                    required.pin.source,
                    source_authority_name(required.pin.authority),
                    required.pin.snapshot.to_hex(),
                    encode_values(&current)
                ),
            )
            .with_fix("restore the exact archived source bytes and authority pin; current metadata is not replay input"));
        }
        Ok(())
    }
}
