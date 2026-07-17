//! fs-obs — structured observability: the ONE event schema for kernels,
//! solvers, test suites, and (once it lands) the ledger `events` table.
//!
//! Roughly forty beads specify per-suite logging ("structured JSON records
//! sufficient to reproduce any failure from logs alone"). Without a single
//! owned schema, forty suites invent forty dialects and that promise dies.
//! This crate makes "diagnosable from logs alone" a CHECKABLE property:
//! every emitter produces [`Event`]s, every event serializes to one JSON-line
//! dialect, and the [`validate_line`] / [`lint_failure_record`] gates run in
//! CI.
//!
//! Determinism split (Decalogue P2): an event's CONTENT (kind, payload,
//! scope, seq) is deterministic in deterministic mode and is what
//! [`Event::content_hash`] covers; wall-clock time lives in the envelope
//! only, excluded from the hash, so logs from two runs of the same seed diff
//! cleanly.
//!
//! Serialization is in-house (Decalogue P1: std + constellation only — serde
//! is not on that list). The wire format is JSON-lines with CANONICAL field
//! order; the strict validator treats deviation as corruption, not dialect.
//!
//! The [`ident`] module owns the CANONICAL REPLAY IDENTITY encoding
//! (bead gp3.14): versioned, typed, length-prefixed — the shared
//! replacement for ad hoc delimiter-concatenation identities.

pub mod ident;

use core::fmt;
use std::fmt::Write as _;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Schema version stamped into every line; bump on any non-additive change.
pub const SCHEMA_VERSION: u32 = 1;

/// Exact typed event-content identity semantics. Version 2 replaces the
/// legacy display-JSON hash, which collapsed distinct NaN payload bits.
pub const EVENT_CONTENT_IDENTITY_VERSION: u32 = 4;

/// Domain-separated artifact kind framed into the typed event identity.
pub const EVENT_CONTENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-obs.event-content.v4";

/// Owner-local declaration consumed by `xtask check-identities`.
pub const EVENT_CONTENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-obs:event-content",
    "version_const=EVENT_CONTENT_IDENTITY_VERSION",
    "version=4",
    "domain=org.frankensim.fs-obs.event-content.v4",
    "domain_const=EVENT_CONTENT_IDENTITY_DOMAIN",
    "encoder=Event::content_identity",
    "encoder_helpers=Event::content_identity_with_versions,Event::content_identity_with_schema,Severity::name,EventKind::kind_name,ReceiptScope::name,ExecutionDisposition::name,PredicateOutcome::name,EpistemicGrade::name,DomainApplicability::name,OperationalSupport::name,EvidenceCompleteness::name,EvidenceIntegrity::name,PromotionEffect::name",
    "schema_constants=EVENT_CONTENT_IDENTITY_VERSION,EVENT_CONTENT_IDENTITY_DOMAIN,SCHEMA_VERSION",
    "schema_functions=check_event_content_identity_version,Event::content_identity_receipt,Event::admit_content_identity,fnv1a64",
    "schema_dependencies=fs-obs:replay-identity-frame",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=Event,EventIdentityReceipt",
    "source_fields=Event.session:semantic,Event.scope:semantic,Event.seq:semantic,Event.severity:semantic,Event.kind:semantic,Event.wall_ns:nonsemantic:wall-clock-envelope-only,EventIdentityReceipt.declared_identity_version:semantic,EventIdentityReceipt.canonical_bytes:semantic,EventIdentityReceipt.root:derived:validated-fnv-root-of-retained-canonical-bytes",
    "source_bindings=Event.session>session,Event.scope>scope,Event.seq>seq,Event.severity>severity,Event.kind>kind+solver-residual-solver+solver-residual-iter+solver-residual-residual+tile-complete-tile+tile-complete-kernel+cancellation-reason+budget-delta-resource+budget-delta-spent+budget-delta-remaining+gradient-check-op+gradient-check-max-rel-err+gradient-check-pass+conformance-case-suite+conformance-case-case+conformance-case-pass+conformance-case-detail+conformance-case-seed+benchmark-result-kernel+benchmark-result-metric+benchmark-result-value+benchmark-result-machine+storm-assertion-name+storm-assertion-pass+storm-assertion-seed+manifest-selection-manifest+manifest-selection-source-snapshot+stratum-expansion-stratum+stratum-expansion-profile+stratum-expansion-cases+dsr-run-run+dsr-run-outcome-scope+dsr-run-outcome-receipt+dsr-run-outcome-disposition+dsr-run-outcome-predicate+dsr-run-outcome-evidence-methods+dsr-run-outcome-grade+dsr-run-outcome-applicability+dsr-run-outcome-support+dsr-run-outcome-completeness+dsr-run-outcome-integrity+dsr-run-outcome-promotion+dsr-run-outcome-detail+race-record-resource+race-record-schedule+race-record-pass+race-record-seed+degradation-event-resource+degradation-event-limit+degradation-event-observed+degradation-event-action+import-receipt-format+import-receipt-artifact+import-receipt-accepted+import-receipt-detail+certificate-verdict-certificate+certificate-verdict-pass+certificate-verdict-bound+certificate-verdict-detail+custom-name+custom-json-exact-opaque-utf8,EventIdentityReceipt.declared_identity_version>retained-producer-version,EventIdentityReceipt.canonical_bytes>retained-canonical-bytes",
    "external_semantic_fields=artifact-domain,identity-version,wire-schema",
    "semantic_fields=artifact-domain,identity-version,wire-schema,session,scope,seq,severity,kind,solver-residual-solver,solver-residual-iter,solver-residual-residual,tile-complete-tile,tile-complete-kernel,cancellation-reason,budget-delta-resource,budget-delta-spent,budget-delta-remaining,gradient-check-op,gradient-check-max-rel-err,gradient-check-pass,conformance-case-suite,conformance-case-case,conformance-case-pass,conformance-case-detail,conformance-case-seed,benchmark-result-kernel,benchmark-result-metric,benchmark-result-value,benchmark-result-machine,storm-assertion-name,storm-assertion-pass,storm-assertion-seed,manifest-selection-manifest,manifest-selection-source-snapshot,stratum-expansion-stratum,stratum-expansion-profile,stratum-expansion-cases,dsr-run-run,dsr-run-outcome-scope,dsr-run-outcome-receipt,dsr-run-outcome-disposition,dsr-run-outcome-predicate,dsr-run-outcome-evidence-methods,dsr-run-outcome-grade,dsr-run-outcome-applicability,dsr-run-outcome-support,dsr-run-outcome-completeness,dsr-run-outcome-integrity,dsr-run-outcome-promotion,dsr-run-outcome-detail,race-record-resource,race-record-schedule,race-record-pass,race-record-seed,degradation-event-resource,degradation-event-limit,degradation-event-observed,degradation-event-action,import-receipt-format,import-receipt-artifact,import-receipt-accepted,import-receipt-detail,certificate-verdict-certificate,certificate-verdict-pass,certificate-verdict-bound,certificate-verdict-detail,custom-name,custom-json-exact-opaque-utf8,retained-producer-version,retained-canonical-bytes",
    "excluded_fields=to-jsonl:display-transport-only",
    "consumers=Event::content_hash,EventIdentityReceipt,Event::admit_content_identity,ledger-event-sinks,replay-comparison",
    "mutations=artifact-domain:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,identity-version:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,wire-schema:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,session:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,scope:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,seq:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,severity:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,kind:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,solver-residual-solver:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,solver-residual-iter:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,solver-residual-residual:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tile-complete-tile:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tile-complete-kernel:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,cancellation-reason:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-spent:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-remaining:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-op:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-max-rel-err:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-suite:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-case:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-kernel:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-metric:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-value:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-machine:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-name:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,manifest-selection-manifest:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,manifest-selection-source-snapshot:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-stratum:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-profile:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-cases:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-run:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-scope:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-receipt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-disposition:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-predicate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-evidence-methods:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-grade:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-applicability:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-support:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-completeness:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-integrity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-promotion:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-schedule:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-limit:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-observed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-action:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-format:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-artifact:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-accepted:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-certificate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-bound:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,custom-name:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,custom-json-exact-opaque-utf8:crates/fs-obs/src/lib.rs#custom_payload_identity_is_exact_opaque_utf8,retained-producer-version:crates/fs-obs/src/lib.rs#retained_event_identity_receipts_admit_exactly_or_fail_closed,retained-canonical-bytes:crates/fs-obs/src/lib.rs#retained_event_identity_receipts_admit_exactly_or_fail_closed",
    "nonsemantic_mutations=Event.wall_ns:crates/fs-obs/src/lib.rs#wall_clock_is_envelope_only,to-jsonl:crates/fs-obs/src/lib.rs#content_identity_preserves_bits_that_display_json_collapses",
    "field_guard=classify_event_identity_fields",
    "transport_guard=Event::admit_content_identity",
    "version_guard=crates/fs-obs/src/lib.rs#event_content_identity_versions_fail_closed",
    "coupling_surface=fs-obs:event-content-identity",
];

/// Structured refusal for event identities produced under unknown semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EventIdentityVersionError {
    /// Version declared by retained evidence.
    pub declared: u32,
    /// Exact version supported by this build.
    pub supported: u32,
}

impl fmt::Display for EventIdentityVersionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "event content identity v{} is unsupported; this build accepts exactly v{}",
            self.declared, self.supported
        )
    }
}

impl core::error::Error for EventIdentityVersionError {}

/// Retained proof of one event's exact content identity.
///
/// A receipt is deliberately more than a naked root: it carries the declared
/// event-identity semantics and the complete canonical preimage. Consumers
/// must call [`Event::admit_content_identity`] before trusting retained data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EventIdentityReceipt {
    declared_identity_version: u32,
    canonical_bytes: Vec<u8>,
    root: u64,
}

impl EventIdentityReceipt {
    /// Reconstruct a receipt loaded from retained evidence.
    ///
    /// This constructor intentionally does not bless the parts. Call
    /// [`Event::admit_content_identity`] to verify the declared version, root,
    /// exact bytes, and event content together.
    #[must_use]
    pub fn from_retained_parts(
        declared_identity_version: u32,
        canonical_bytes: Vec<u8>,
        root: u64,
    ) -> Self {
        Self {
            declared_identity_version,
            canonical_bytes,
            root,
        }
    }

    /// Event-content identity semantics declared by the retained producer.
    #[must_use]
    pub fn declared_identity_version(&self) -> u32 {
        self.declared_identity_version
    }

    /// Complete canonical typed preimage retained by the producer.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Retained FNV-1a root over [`Self::canonical_bytes`].
    #[must_use]
    pub fn root(&self) -> u64 {
        self.root
    }
}

/// Fail-closed refusal when retained event-identity evidence is not exactly
/// the identity of the event being admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventIdentityAdmissionError {
    /// The producer declared identity semantics unknown to this build.
    UnsupportedVersion(EventIdentityVersionError),
    /// The retained root is not derived from the retained canonical bytes.
    RootMismatch {
        /// Root recorded in the receipt.
        declared: u64,
        /// Root recomputed from the receipt's canonical bytes.
        computed: u64,
    },
    /// The self-consistent receipt names different exact event content.
    CanonicalBytesMismatch {
        /// Root recorded in the receipt.
        declared_root: u64,
        /// Root computed from the event supplied to admission.
        expected_root: u64,
    },
}

impl fmt::Display for EventIdentityAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedVersion(error) => fmt::Display::fmt(error, f),
            Self::RootMismatch { declared, computed } => write!(
                f,
                "retained event identity root {declared:016x} does not match canonical bytes root {computed:016x}"
            ),
            Self::CanonicalBytesMismatch {
                declared_root,
                expected_root,
            } => write!(
                f,
                "retained event identity {declared_root:016x} does not bind the admitted event identity {expected_root:016x}"
            ),
        }
    }
}

impl core::error::Error for EventIdentityAdmissionError {}

/// Refuse retained event identities whose exact typed semantics are unknown.
///
/// # Errors
/// [`EventIdentityVersionError`] for any version other than the current one.
pub fn check_event_content_identity_version(
    declared: u32,
) -> Result<(), EventIdentityVersionError> {
    if declared == EVENT_CONTENT_IDENTITY_VERSION {
        Ok(())
    } else {
        Err(EventIdentityVersionError {
            declared,
            supported: EVENT_CONTENT_IDENTITY_VERSION,
        })
    }
}

/// Receipt scope for one verification outcome (i94v.7.3.1): a successful
/// status query, a failed worker attempt, and a refuted scientific job can
/// coexist, so every outcome names the exact scope its receipt covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReceiptScope {
    /// One scoped operation.
    Operation,
    /// One worker attempt.
    Attempt,
    /// One job across attempts.
    Job,
    /// One whole campaign.
    Campaign,
}

impl ReceiptScope {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Operation => "operation",
            Self::Attempt => "attempt",
            Self::Job => "job",
            Self::Campaign => "campaign",
        }
    }
}

/// How the scoped execution ended, independent of what it proved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionDisposition {
    /// Ran to completion.
    Completed,
    /// Ran and failed.
    Failed,
    /// Cancelled before completion.
    Cancelled,
    /// Exceeded its budget/deadline.
    TimedOut,
    /// Refused before execution.
    Refused,
    /// Ended in a state the harness cannot classify; closure incomplete.
    Indeterminate,
}

impl ExecutionDisposition {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Completed => "completed",
            Self::Failed => "failed",
            Self::Cancelled => "cancelled",
            Self::TimedOut => "timed_out",
            Self::Refused => "refused",
            Self::Indeterminate => "indeterminate",
        }
    }
}

/// What the requested predicate concluded, independent of how it ran.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PredicateOutcome {
    /// The predicate held.
    Satisfied,
    /// The predicate was refuted.
    Refuted,
    /// Ran but could not decide.
    Indeterminate,
    /// This outcome carries no predicate.
    NotApplicable,
}

impl PredicateOutcome {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Satisfied => "satisfied",
            Self::Refuted => "refuted",
            Self::Indeterminate => "indeterminate",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Epistemic grade of the produced evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpistemicGrade {
    /// Machine-checked enclosure/certificate.
    Verified,
    /// Validated against admitted external data.
    Validated,
    /// Estimated by an uncertified method.
    Estimated,
    /// Reported without evidence machinery.
    Reported,
    /// No evidence claim in this outcome.
    NotApplicable,
}

impl EpistemicGrade {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Validated => "validated",
            Self::Estimated => "estimated",
            Self::Reported => "reported",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the run stayed inside its declared applicability domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DomainApplicability {
    /// In-domain.
    InDomain,
    /// Out of the declared domain.
    OutOfDomain,
    /// Applicability was not evaluated.
    Unevaluated,
    /// No domain claim applies.
    NotApplicable,
}

impl DomainApplicability {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::InDomain => "in_domain",
            Self::OutOfDomain => "out_of_domain",
            Self::Unevaluated => "unevaluated",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Operational support level under which the outcome was produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperationalSupport {
    /// Fully supported configuration.
    Supported,
    /// Degraded configuration; outcome stands with caveats.
    Degraded,
    /// Unsupported configuration.
    Unsupported,
    /// Support classification does not apply.
    NotApplicable,
}

impl OperationalSupport {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Supported => "supported",
            Self::Degraded => "degraded",
            Self::Unsupported => "unsupported",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the evidence set is complete (SEPARATE from integrity).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceCompleteness {
    /// Every required artifact/receipt is present.
    Complete,
    /// Required evidence is missing; closure incomplete.
    Incomplete,
    /// Completeness does not apply.
    NotApplicable,
}

impl EvidenceCompleteness {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Incomplete => "incomplete",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the evidence set is intact (SEPARATE from completeness).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceIntegrity {
    /// All present evidence authenticated.
    Intact,
    /// Tamper/corruption detected.
    Compromised,
    /// Integrity not yet checked.
    Unchecked,
    /// Integrity does not apply.
    NotApplicable,
}

impl EvidenceIntegrity {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Intact => "intact",
            Self::Compromised => "compromised",
            Self::Unchecked => "unchecked",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Effect of this outcome on promotion state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PromotionEffect {
    /// Promoted evidence/authority.
    Promoted,
    /// Demoted evidence/authority.
    Demoted,
    /// No promotion movement.
    Unchanged,
    /// Promotion does not apply.
    NotApplicable,
}

impl PromotionEffect {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Promoted => "promoted",
            Self::Demoted => "demoted",
            Self::Unchanged => "unchanged",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// One scoped verification outcome (i94v.7.3.1): the shared core every
/// outcome-bearing event kind embeds. Not-applicable states are explicit
/// variants, never absent-by-convention fields.
#[derive(Debug, Clone, PartialEq)]
pub struct ScopedReceiptOutcome {
    /// Scope the receipt covers.
    pub scope: ReceiptScope,
    /// Receipt identity (hex) at that scope.
    pub receipt: String,
    /// How execution ended.
    pub disposition: ExecutionDisposition,
    /// What the requested predicate concluded.
    pub predicate: PredicateOutcome,
    /// Comma-joined registered evidence-method names ("" = none).
    pub evidence_methods: String,
    /// Epistemic grade of the produced evidence.
    pub grade: EpistemicGrade,
    /// Domain-applicability classification.
    pub applicability: DomainApplicability,
    /// Operational-support classification.
    pub support: OperationalSupport,
    /// Evidence completeness (separate from integrity).
    pub completeness: EvidenceCompleteness,
    /// Evidence integrity (separate from completeness).
    pub integrity: EvidenceIntegrity,
    /// Promotion effect.
    pub promotion: PromotionEffect,
    /// Diagnostic detail; REQUIRED non-empty for failed/refused/
    /// indeterminate dispositions (lint-enforced).
    pub detail: String,
}

/// Severity ladder. `Error` events MUST satisfy [`lint_failure_record`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// High-volume diagnostics (tile completions).
    Trace,
    /// Normal progress (solver iterations, case verdicts).
    Info,
    /// Something degraded but the run continues.
    Warn,
    /// A failure record — must be reproducible from its own payload.
    Error,
}

impl Severity {
    fn name(self) -> &'static str {
        match self {
            Severity::Trace => "trace",
            Severity::Info => "info",
            Severity::Warn => "warn",
            Severity::Error => "error",
        }
    }
}

/// Typed payload registry, v1. Additive evolution only: new kinds may be
/// added; existing fields may never change meaning (validator-enforced by
/// the golden lines in the conformance suite).
#[derive(Debug, Clone, PartialEq)]
pub enum EventKind {
    /// One iteration of an iterative solver.
    SolverResidual {
        /// Solver identity (e.g. "cg", "p-mg").
        solver: String,
        /// Iteration index.
        iter: u64,
        /// Residual norm.
        residual: f64,
    },
    /// A tile finished executing.
    TileComplete {
        /// Logical tile identity (the SAME identity that keys RNG streams
        /// and deterministic reductions — plan §5.2).
        tile: u64,
        /// Kernel name.
        kernel: String,
    },
    /// A scope was cancelled.
    Cancellation {
        /// Why (kill-handle, budget, panic containment, ...).
        reason: String,
    },
    /// Budget accounting delta (P4: budgets first).
    BudgetDelta {
        /// Resource name ("wall_s", "mem_bytes", "core_s", "energy_j").
        resource: String,
        /// Amount spent in this step.
        spent: f64,
        /// Remaining grant.
        remaining: f64,
    },
    /// A gradient verification outcome (the merge-gate evidence, §8.7).
    GradientCheck {
        /// Operator under test.
        op: String,
        /// Max relative error across probed directions.
        max_rel_err: f64,
        /// Verdict.
        pass: bool,
    },
    /// One conformance-suite case verdict (plan §13.3).
    ConformanceCase {
        /// Suite id (e.g. "fs-qty/conformance").
        suite: String,
        /// Case id (e.g. "qty-001/0.12Pa*s").
        case: String,
        /// Verdict.
        pass: bool,
        /// Human/agent-readable detail; REQUIRED non-empty when `pass=false`
        /// (lint-enforced) so failures reproduce from the log alone.
        detail: String,
        /// Replay seed when the case is randomized (lint: required on fail
        /// for randomized cases; 0 = not randomized).
        seed: u64,
    },
    /// A performance measurement (roofline harness rows).
    BenchmarkResult {
        /// Kernel name.
        kernel: String,
        /// Metric name ("glups", "gflops", "bandwidth_gbs", "mrays").
        metric: String,
        /// Measured value.
        value: f64,
        /// Machine fingerprint hash (fs-substrate probe).
        machine: u64,
    },
    /// A chaos/storm assertion outcome (G4).
    StormAssertion {
        /// Assertion name ("no-arena-leak", "cancel-latency-p99").
        name: String,
        /// Verdict.
        pass: bool,
        /// Storm seed for replay.
        seed: u64,
    },
    /// Manifest + source-snapshot selection for one verification run
    /// (i94v.7.3.1 F1).
    ManifestSelection {
        /// Selected manifest identity (hex).
        manifest: String,
        /// Source-snapshot identity (hex) the manifest was resolved against.
        source_snapshot: String,
    },
    /// Core/Max stratum and campaign-profile expansion (i94v.7.3.1 F1).
    StratumExpansion {
        /// Stratum name ("core", "max").
        stratum: String,
        /// Campaign profile the stratum expanded under.
        profile: String,
        /// Number of expanded cases.
        cases: u64,
    },
    /// One DSR run outcome carrying the full scoped receipt core
    /// (i94v.7.3.1 F1).
    DsrRun {
        /// DSR run identity.
        run: String,
        /// The scoped verification outcome.
        outcome: ScopedReceiptOutcome,
    },
    /// A concurrency-race observation from a G4 race/storm harness.
    RaceRecord {
        /// Contended resource or invariant ("arena-slot", "tune-row", ...).
        resource: String,
        /// Interleaving/schedule descriptor sufficient to replay the race.
        schedule: String,
        /// Verdict.
        pass: bool,
        /// Replay seed for the schedule (lint: required non-zero on fail).
        seed: u64,
    },
    /// A graceful-degradation decision at a bounded resource.
    DegradationEvent {
        /// Resource that hit its bound ("retained_bytes_per_scope", ...).
        resource: String,
        /// Bound in force.
        limit: u64,
        /// Observed demand when the bound fired.
        observed: u64,
        /// Action taken ("rotate-lane", "flush", "refuse", ...).
        action: String,
    },
    /// An external-artifact import decision (STEP/mesh/dataset quarantine).
    ImportReceipt {
        /// Import format ("step", "gltf", "csv-dataset", ...).
        format: String,
        /// Content hash of the imported artifact (hex).
        artifact: String,
        /// Whether the importer admitted the artifact.
        accepted: bool,
        /// Refusal/quarantine detail; REQUIRED non-empty when
        /// `accepted=false` (lint-enforced).
        detail: String,
    },
    /// A certificate check outcome (enclosure/trim/conversion verdicts).
    CertificateVerdict {
        /// Certificate family ("ivl-enclosure", "trim-cert", ...).
        certificate: String,
        /// Verdict.
        pass: bool,
        /// Certified numeric bound (0.0 when the certificate is non-numeric).
        bound: f64,
        /// Detail; REQUIRED non-empty when `pass=false` (lint-enforced).
        detail: String,
    },
    /// Escape hatch for kinds not yet in the registry. `json` is exact opaque
    /// UTF-8: whitespace and object-member order are semantic identity bytes;
    /// fs-obs never claims to canonicalize unchecked JSON. The caller must
    /// still supply one valid pre-serialized JSON object for `to_jsonl`.
    Custom {
        /// Kind name (kebab-case).
        name: String,
        /// Exact opaque UTF-8 bytes of the pre-serialized JSON object.
        json: String,
    },
}

impl EventKind {
    /// Stable wire/ledger kind name. This is the documented export surface
    /// for ledger `events`-table sinks (huq.16): sinks store
    /// (`Event::to_jsonl` bytes, `kind_name`) without re-encoding.
    #[must_use]
    pub fn kind_name(&self) -> &'static str {
        match self {
            EventKind::SolverResidual { .. } => "solver_residual",
            EventKind::TileComplete { .. } => "tile_complete",
            EventKind::Cancellation { .. } => "cancellation",
            EventKind::BudgetDelta { .. } => "budget_delta",
            EventKind::GradientCheck { .. } => "gradient_check",
            EventKind::ConformanceCase { .. } => "conformance_case",
            EventKind::BenchmarkResult { .. } => "benchmark_result",
            EventKind::StormAssertion { .. } => "storm_assertion",
            EventKind::ManifestSelection { .. } => "manifest_selection",
            EventKind::StratumExpansion { .. } => "stratum_expansion",
            EventKind::DsrRun { .. } => "dsr_run",
            EventKind::RaceRecord { .. } => "race_record",
            EventKind::DegradationEvent { .. } => "degradation_event",
            EventKind::ImportReceipt { .. } => "import_receipt",
            EventKind::CertificateVerdict { .. } => "certificate_verdict",
            EventKind::Custom { .. } => "custom",
        }
    }
}

/// One observability event: envelope + typed payload.
#[derive(Debug, Clone, PartialEq)]
pub struct Event {
    /// Session identity (study/run scope; content-hashed).
    pub session: String,
    /// Slash-separated scope path mirroring the asupersync scope tree
    /// (e.g. "study-x/op-3/kernel-lbm/tile-42"); content-hashed.
    pub scope: String,
    /// Per-emitter monotone sequence number; content-hashed (gives logs a
    /// deterministic total order per scope without wall-clock).
    pub seq: u64,
    /// Severity.
    pub severity: Severity,
    /// Typed payload.
    pub kind: EventKind,
    /// Wall-clock nanoseconds since the unix epoch. ENVELOPE ONLY: excluded
    /// from `content_hash`, always serialized LAST, `None` in deterministic
    /// replay comparisons.
    pub wall_ns: Option<u64>,
}

#[allow(dead_code)]
fn classify_event_identity_fields(event: &Event, receipt: &EventIdentityReceipt) {
    let Event {
        session: _,
        scope: _,
        seq: _,
        severity: _,
        kind: _,
        wall_ns: _,
    } = event;
    let EventIdentityReceipt {
        declared_identity_version: _,
        canonical_bytes: _,
        root: _,
    } = receipt;
}

fn esc(out: &mut String, s: &str) {
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
}

fn push_str_field(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, "\"{key}\":\"");
    esc(out, val);
    out.push('"');
}

/// Serialize a float for the wire: finite → shortest-round-trip; non-finite
/// → tagged string (JSON has no NaN/Inf; readers must handle both shapes).
fn push_f64(out: &mut String, key: &str, v: f64) {
    if v.is_finite() {
        let _ = write!(out, "\"{key}\":{v}");
    } else {
        let _ = write!(out, "\"{key}\":\"non-finite:{v}\"");
    }
}

impl Event {
    /// Serialize the CONTENT portion (everything except `wall_ns`) in
    /// canonical field order.
    fn content_json(&self) -> String {
        let mut s = String::with_capacity(160);
        let _ = write!(s, "{{\"v\":{SCHEMA_VERSION},");
        push_str_field(&mut s, "session", &self.session);
        s.push(',');
        push_str_field(&mut s, "scope", &self.scope);
        let _ = write!(s, ",\"seq\":{},", self.seq);
        push_str_field(&mut s, "severity", self.severity.name());
        s.push(',');
        push_str_field(&mut s, "kind", self.kind.kind_name());
        s.push_str(",\"payload\":{");
        match &self.kind {
            EventKind::SolverResidual {
                solver,
                iter,
                residual,
            } => {
                push_str_field(&mut s, "solver", solver);
                let _ = write!(s, ",\"iter\":{iter},");
                push_f64(&mut s, "residual", *residual);
            }
            EventKind::TileComplete { tile, kernel } => {
                let _ = write!(s, "\"tile\":{tile},");
                push_str_field(&mut s, "kernel", kernel);
            }
            EventKind::Cancellation { reason } => {
                push_str_field(&mut s, "reason", reason);
            }
            EventKind::BudgetDelta {
                resource,
                spent,
                remaining,
            } => {
                push_str_field(&mut s, "resource", resource);
                s.push(',');
                push_f64(&mut s, "spent", *spent);
                s.push(',');
                push_f64(&mut s, "remaining", *remaining);
            }
            EventKind::GradientCheck {
                op,
                max_rel_err,
                pass,
            } => {
                push_str_field(&mut s, "op", op);
                s.push(',');
                push_f64(&mut s, "max_rel_err", *max_rel_err);
                let _ = write!(s, ",\"pass\":{pass}");
            }
            EventKind::ConformanceCase {
                suite,
                case,
                pass,
                detail,
                seed,
            } => {
                push_str_field(&mut s, "suite", suite);
                s.push(',');
                push_str_field(&mut s, "case", case);
                let _ = write!(s, ",\"pass\":{pass},");
                push_str_field(&mut s, "detail", detail);
                let _ = write!(s, ",\"seed\":{seed}");
            }
            EventKind::BenchmarkResult {
                kernel,
                metric,
                value,
                machine,
            } => {
                push_str_field(&mut s, "kernel", kernel);
                s.push(',');
                push_str_field(&mut s, "metric", metric);
                s.push(',');
                push_f64(&mut s, "value", *value);
                let _ = write!(s, ",\"machine\":{machine}");
            }
            EventKind::StormAssertion { name, pass, seed } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"pass\":{pass},\"seed\":{seed}");
            }
            EventKind::ManifestSelection {
                manifest,
                source_snapshot,
            } => {
                push_str_field(&mut s, "manifest", manifest);
                s.push(',');
                push_str_field(&mut s, "source_snapshot", source_snapshot);
            }
            EventKind::StratumExpansion {
                stratum,
                profile,
                cases,
            } => {
                push_str_field(&mut s, "stratum", stratum);
                s.push(',');
                push_str_field(&mut s, "profile", profile);
                let _ = write!(s, ",\"cases\":{cases}");
            }
            EventKind::DsrRun { run, outcome } => {
                push_str_field(&mut s, "run", run);
                s.push_str(",\"outcome\":{");
                push_str_field(&mut s, "scope", outcome.scope.name());
                s.push(',');
                push_str_field(&mut s, "receipt", &outcome.receipt);
                s.push(',');
                push_str_field(&mut s, "disposition", outcome.disposition.name());
                s.push(',');
                push_str_field(&mut s, "predicate", outcome.predicate.name());
                s.push(',');
                push_str_field(&mut s, "evidence_methods", &outcome.evidence_methods);
                s.push(',');
                push_str_field(&mut s, "grade", outcome.grade.name());
                s.push(',');
                push_str_field(&mut s, "applicability", outcome.applicability.name());
                s.push(',');
                push_str_field(&mut s, "support", outcome.support.name());
                s.push(',');
                push_str_field(&mut s, "completeness", outcome.completeness.name());
                s.push(',');
                push_str_field(&mut s, "integrity", outcome.integrity.name());
                s.push(',');
                push_str_field(&mut s, "promotion", outcome.promotion.name());
                s.push(',');
                push_str_field(&mut s, "detail", &outcome.detail);
                s.push('}');
            }
            EventKind::RaceRecord {
                resource,
                schedule,
                pass,
                seed,
            } => {
                push_str_field(&mut s, "resource", resource);
                s.push(',');
                push_str_field(&mut s, "schedule", schedule);
                let _ = write!(s, ",\"pass\":{pass},\"seed\":{seed}");
            }
            EventKind::DegradationEvent {
                resource,
                limit,
                observed,
                action,
            } => {
                push_str_field(&mut s, "resource", resource);
                let _ = write!(s, ",\"limit\":{limit},\"observed\":{observed},");
                push_str_field(&mut s, "action", action);
            }
            EventKind::ImportReceipt {
                format,
                artifact,
                accepted,
                detail,
            } => {
                push_str_field(&mut s, "format", format);
                s.push(',');
                push_str_field(&mut s, "artifact", artifact);
                let _ = write!(s, ",\"accepted\":{accepted},");
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::CertificateVerdict {
                certificate,
                pass,
                bound,
                detail,
            } => {
                push_str_field(&mut s, "certificate", certificate);
                let _ = write!(s, ",\"pass\":{pass},");
                push_f64(&mut s, "bound", *bound);
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::Custom { name, json } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"data\":{json}");
            }
        }
        s.push('}');
        s.push('}');
        s
    }

    /// Serialize one JSON-line: content with `wall_ns` spliced in as the
    /// LAST field before the closing brace (envelope position).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let mut content = self.content_json();
        if let Some(ns) = self.wall_ns {
            content.pop(); // strip trailing '}'
            let _ = write!(content, ",\"wall_ns\":{ns}}}");
        }
        content
    }

    fn content_identity_with_versions(
        &self,
        event_identity_version: u32,
        event_wire_schema_version: u32,
    ) -> ident::ReplayIdentity {
        self.content_identity_with_schema(
            EVENT_CONTENT_IDENTITY_DOMAIN,
            event_identity_version,
            event_wire_schema_version,
        )
    }

    fn content_identity_with_schema(
        &self,
        artifact_domain: &str,
        event_identity_version: u32,
        event_wire_schema_version: u32,
    ) -> ident::ReplayIdentity {
        let Event {
            session,
            scope,
            seq,
            severity,
            kind,
            wall_ns: _,
        } = self;
        let builder = ident::IdentityBuilder::new(artifact_domain)
            .u64(
                "event_content_identity_version",
                u64::from(event_identity_version),
            )
            .u64(
                "event_wire_schema_version",
                u64::from(event_wire_schema_version),
            )
            .str("session", session)
            .str("scope", scope)
            .u64("seq", *seq)
            .str("severity", severity.name())
            .str("kind", kind.kind_name());
        let builder = match kind {
            EventKind::SolverResidual {
                solver,
                iter,
                residual,
            } => builder
                .str("solver", solver)
                .u64("iter", *iter)
                .f64_bits("residual", *residual),
            EventKind::TileComplete { tile, kernel } => {
                builder.u64("tile", *tile).str("kernel", kernel)
            }
            EventKind::Cancellation { reason } => builder.str("reason", reason),
            EventKind::BudgetDelta {
                resource,
                spent,
                remaining,
            } => builder
                .str("resource", resource)
                .f64_bits("spent", *spent)
                .f64_bits("remaining", *remaining),
            EventKind::GradientCheck {
                op,
                max_rel_err,
                pass,
            } => builder
                .str("op", op)
                .f64_bits("max_rel_err", *max_rel_err)
                .flag("pass", *pass),
            EventKind::ConformanceCase {
                suite,
                case,
                pass,
                detail,
                seed,
            } => builder
                .str("suite", suite)
                .str("case", case)
                .flag("pass", *pass)
                .str("detail", detail)
                .u64("seed", *seed),
            EventKind::BenchmarkResult {
                kernel,
                metric,
                value,
                machine,
            } => builder
                .str("kernel", kernel)
                .str("metric", metric)
                .f64_bits("value", *value)
                .u64("machine", *machine),
            EventKind::StormAssertion { name, pass, seed } => builder
                .str("name", name)
                .flag("pass", *pass)
                .u64("seed", *seed),
            EventKind::ManifestSelection {
                manifest,
                source_snapshot,
            } => builder
                .str("manifest", manifest)
                .str("source_snapshot", source_snapshot),
            EventKind::StratumExpansion {
                stratum,
                profile,
                cases,
            } => builder
                .str("stratum", stratum)
                .str("profile", profile)
                .u64("cases", *cases),
            EventKind::DsrRun { run, outcome } => builder
                .str("run", run)
                .str("outcome_scope", outcome.scope.name())
                .str("outcome_receipt", &outcome.receipt)
                .str("outcome_disposition", outcome.disposition.name())
                .str("outcome_predicate", outcome.predicate.name())
                .str("outcome_evidence_methods", &outcome.evidence_methods)
                .str("outcome_grade", outcome.grade.name())
                .str("outcome_applicability", outcome.applicability.name())
                .str("outcome_support", outcome.support.name())
                .str("outcome_completeness", outcome.completeness.name())
                .str("outcome_integrity", outcome.integrity.name())
                .str("outcome_promotion", outcome.promotion.name())
                .str("outcome_detail", &outcome.detail),
            EventKind::RaceRecord {
                resource,
                schedule,
                pass,
                seed,
            } => builder
                .str("resource", resource)
                .str("schedule", schedule)
                .flag("pass", *pass)
                .u64("seed", *seed),
            EventKind::DegradationEvent {
                resource,
                limit,
                observed,
                action,
            } => builder
                .str("resource", resource)
                .u64("limit", *limit)
                .u64("observed", *observed)
                .str("action", action),
            EventKind::ImportReceipt {
                format,
                artifact,
                accepted,
                detail,
            } => builder
                .str("format", format)
                .str("artifact", artifact)
                .flag("accepted", *accepted)
                .str("detail", detail),
            EventKind::CertificateVerdict {
                certificate,
                pass,
                bound,
                detail,
            } => builder
                .str("certificate", certificate)
                .flag("pass", *pass)
                .f64_bits("bound", *bound)
                .str("detail", detail),
            EventKind::Custom { name, json } => builder
                .str("name", name)
                .bytes("custom_json_opaque_utf8", json.as_bytes()),
        };
        builder
            .exclude(
                "wall_ns",
                "wall-clock is observability envelope, not replay identity",
            )
            .finish()
    }

    /// Canonical typed identity for the deterministic event content.
    ///
    /// This is deliberately independent of the display-oriented JSON line:
    /// floats bind by their exact bit patterns, every payload variant has a
    /// closed typed encoder, and wall-clock remains an explicit exclusion.
    #[must_use]
    pub fn content_identity(&self) -> ident::ReplayIdentity {
        self.content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION)
    }

    /// Capture a retained receipt containing the declared event-identity
    /// version, exact canonical bytes, and their root.
    #[must_use]
    pub fn content_identity_receipt(&self) -> EventIdentityReceipt {
        let identity = self.content_identity();
        EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            identity.canonical_bytes().to_vec(),
            identity.root(),
        )
    }

    /// Admit a retained identity receipt only when all three proof surfaces
    /// agree: declared semantics, root-over-retained-bytes, and the exact
    /// canonical identity of this event.
    ///
    /// # Errors
    /// [`EventIdentityAdmissionError`] if any proof surface differs.
    pub fn admit_content_identity(
        &self,
        receipt: &EventIdentityReceipt,
    ) -> Result<(), EventIdentityAdmissionError> {
        check_event_content_identity_version(receipt.declared_identity_version)
            .map_err(EventIdentityAdmissionError::UnsupportedVersion)?;

        let computed = fnv1a64(&receipt.canonical_bytes);
        if computed != receipt.root {
            return Err(EventIdentityAdmissionError::RootMismatch {
                declared: receipt.root,
                computed,
            });
        }

        let expected = self.content_identity();
        if receipt.canonical_bytes.as_slice() != expected.canonical_bytes()
            || receipt.root != expected.root()
        {
            return Err(EventIdentityAdmissionError::CanonicalBytesMismatch {
                declared_root: receipt.root,
                expected_root: expected.root(),
            });
        }
        Ok(())
    }

    /// Deterministic FNV-1a root over [`Event::content_identity`]'s exact
    /// typed bytes. Not cryptographic; ledger-grade content addressing uses
    /// the same canonical identity bytes under a stronger digest.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        self.content_identity().root()
    }
}

/// FNV-1a 64-bit (in-house, deterministic across platforms).
#[must_use]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// A per-scope emitter handing out monotone sequence numbers.
#[derive(Debug)]
pub struct Emitter {
    session: String,
    scope: String,
    /// Byte length of the immutable base scope; `exit_scope` never
    /// truncates past it.
    base_scope_len: usize,
    seq: u64,
}

impl Emitter {
    /// Create an emitter for one (session, scope) pair.
    #[must_use]
    pub fn new(session: impl Into<String>, scope: impl Into<String>) -> Self {
        let scope = scope.into();
        let base_scope_len = scope.len();
        Emitter {
            session: session.into(),
            scope,
            base_scope_len,
            seq: 0,
        }
    }

    /// Enter a child scope segment: the asupersync scope tree IS the trace
    /// tree (huq.16), so the emitter mirrors it with explicit enter/exit
    /// rather than thread-local magic. Entering emits nothing by itself;
    /// `seq` stays one monotone stream per emitter so interleaved child
    /// scopes replay in exact emission order.
    ///
    /// # Errors
    /// Refuses empty segments and segments containing `/` (the path
    /// separator) or control characters — a forged segment must not be able
    /// to impersonate a different tree position.
    pub fn enter_scope(&mut self, segment: &str) -> Result<(), SchemaError> {
        if segment.is_empty() {
            return Err(SchemaError {
                at: 0,
                message: "scope segment must be non-empty".to_string(),
            });
        }
        if segment.contains('/') || segment.chars().any(char::is_control) {
            return Err(SchemaError {
                at: 0,
                message: format!(
                    "scope segment {segment:?} must not contain '/' or control characters \
                     (path forgery guard)"
                ),
            });
        }
        self.scope.push('/');
        self.scope.push_str(segment);
        Ok(())
    }

    /// Exit the innermost child scope.
    ///
    /// # Errors
    /// Refuses at the base scope: exits must pair with entries, and a
    /// mismatched exit is a harness bug worth surfacing, not masking.
    pub fn exit_scope(&mut self) -> Result<(), SchemaError> {
        let Some(cut) = self.scope[self.base_scope_len..].rfind('/') else {
            return Err(SchemaError {
                at: 0,
                message: "exit_scope at the base scope: unbalanced enter/exit".to_string(),
            });
        };
        self.scope.truncate(self.base_scope_len + cut);
        Ok(())
    }

    /// The live slash-separated scope path events are currently stamped with.
    #[must_use]
    pub fn current_scope(&self) -> &str {
        &self.scope
    }

    /// Depth of entered child segments below the base scope.
    #[must_use]
    pub fn scope_depth(&self) -> usize {
        self.scope[self.base_scope_len..]
            .chars()
            .filter(|&c| c == '/')
            .count()
    }

    /// Build the next event (seq auto-increments). `wall_ns` is supplied by
    /// the caller because THIS crate must stay deterministic and I/O-free;
    /// runtime layers pass real clocks, tests pass `None`.
    pub fn emit(&mut self, severity: Severity, kind: EventKind, wall_ns: Option<u64>) -> Event {
        let e = Event {
            session: self.session.clone(),
            scope: self.scope.clone(),
            seq: self.seq,
            severity,
            kind,
            wall_ns,
        };
        self.seq += 1;
        e
    }
}

/// Validation failure for a JSON-line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaError {
    /// Byte position (approximate for structural errors).
    pub at: usize,
    /// What is wrong and how to fix it.
    pub message: String,
}

impl fmt::Display for SchemaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "event schema violation at byte {}: {}",
            self.at, self.message
        )
    }
}

impl core::error::Error for SchemaError {}

/// Registry of known kinds (validator uses it; keep in sync with EventKind).
pub const KNOWN_KINDS: &[&str] = &[
    "solver_residual",
    "tile_complete",
    "cancellation",
    "budget_delta",
    "gradient_check",
    "conformance_case",
    "benchmark_result",
    "storm_assertion",
    "manifest_selection",
    "stratum_expansion",
    "dsr_run",
    "race_record",
    "degradation_event",
    "import_receipt",
    "certificate_verdict",
    "custom",
];

/// Strict structural validation of one JSON-line: required envelope keys in
/// canonical order, known kind, balanced payload object. (The writer is
/// ours; deviation means corruption, so the checks are cheap and strict
/// rather than a full JSON parser.)
///
/// # Errors
/// Returns [`SchemaError`] naming the first violated requirement.
pub fn validate_line(line: &str) -> Result<(), SchemaError> {
    let need = |cond: bool, at: usize, msg: &str| -> Result<(), SchemaError> {
        if cond {
            Ok(())
        } else {
            Err(SchemaError {
                at,
                message: msg.to_string(),
            })
        }
    };
    need(
        line.starts_with('{') && line.ends_with('}'),
        0,
        "line must be one JSON object",
    )?;
    let ver = format!("{{\"v\":{SCHEMA_VERSION},");
    need(
        line.starts_with(&ver),
        0,
        "first field must be the schema version \"v\"",
    )?;
    for key in [
        "\"session\":",
        "\"scope\":",
        "\"seq\":",
        "\"severity\":",
        "\"kind\":",
        "\"payload\":{",
    ] {
        need(
            line.contains(key),
            0,
            &format!("missing required field {key}"),
        )?;
    }
    let kind_pos = line.find("\"kind\":\"").ok_or(SchemaError {
        at: 0,
        message: "kind must be a string".to_string(),
    })?;
    let kind_rest = &line[kind_pos + 8..];
    let kind_end = kind_rest.find('"').unwrap_or(0);
    let kind = &kind_rest[..kind_end];
    need(
        KNOWN_KINDS.contains(&kind),
        kind_pos,
        &format!("unknown kind {kind:?}; register new kinds in the payload registry"),
    )?;
    // Braces balance (escaped quotes handled by the writer's escaping).
    let mut depth = 0i64;
    let mut in_str = false;
    let mut prev_backslash = false;
    for (i, c) in line.char_indices() {
        match c {
            '"' if !prev_backslash => in_str = !in_str,
            '{' if !in_str => depth += 1,
            '}' if !in_str => depth -= 1,
            _ => {}
        }
        prev_backslash = c == '\\' && !prev_backslash;
        if depth < 0 {
            return Err(SchemaError {
                at: i,
                message: "unbalanced braces".to_string(),
            });
        }
    }
    need(
        depth == 0 && !in_str,
        line.len(),
        "unbalanced braces or unterminated string",
    )?;
    Ok(())
}

/// The failure-record completeness lint, v1 (the "log-replay lint" of the
/// observability bead): failure verdicts must carry enough to reproduce.
///
/// # Errors
/// Returns [`SchemaError`] describing the missing reproduction ingredient.
pub fn lint_failure_record(event: &Event) -> Result<(), SchemaError> {
    match &event.kind {
        EventKind::ConformanceCase {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing conformance case must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::StormAssertion {
            pass: false, seed, ..
        } if *seed == 0 => Err(SchemaError {
            at: 0,
            message: "failing storm assertion must carry its replay seed".to_string(),
        }),
        EventKind::GradientCheck {
            pass: false, op, ..
        } if op.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing gradient check must name its operator".to_string(),
        }),
        EventKind::DsrRun { outcome, .. }
            if matches!(
                outcome.disposition,
                ExecutionDisposition::Failed
                    | ExecutionDisposition::Refused
                    | ExecutionDisposition::Indeterminate
            ) && outcome.detail.is_empty() =>
        {
            Err(SchemaError {
                at: 0,
                message: "failed/refused/indeterminate scoped outcome must carry a \
                          non-empty diagnostic detail (reproduce-from-log-alone doctrine)"
                    .to_string(),
            })
        }
        EventKind::RaceRecord {
            pass: false, seed, ..
        } if *seed == 0 => Err(SchemaError {
            at: 0,
            message: "failing race record must carry its schedule replay seed".to_string(),
        }),
        EventKind::ImportReceipt {
            accepted: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "refused import must carry a non-empty quarantine detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::CertificateVerdict {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing certificate verdict must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_events() -> Vec<Event> {
        let mut em = Emitter::new("study-x", "op-1/kernel-cg");
        vec![
            em.emit(
                Severity::Info,
                EventKind::SolverResidual {
                    solver: "cg".into(),
                    iter: 3,
                    residual: 1.5e-7,
                },
                Some(1_000),
            ),
            em.emit(
                Severity::Trace,
                EventKind::TileComplete {
                    tile: 42,
                    kernel: "lbm_d3q19".into(),
                },
                None,
            ),
            em.emit(
                Severity::Warn,
                EventKind::Cancellation {
                    reason: "budget".into(),
                },
                Some(1_500),
            ),
            em.emit(
                Severity::Warn,
                EventKind::BudgetDelta {
                    resource: "wall_s".into(),
                    spent: 12.5,
                    remaining: 7187.5,
                },
                Some(2_000),
            ),
            em.emit(
                Severity::Info,
                EventKind::GradientCheck {
                    op: "poisson".into(),
                    max_rel_err: 2.5e-8,
                    pass: true,
                },
                None,
            ),
            em.emit(
                Severity::Error,
                EventKind::ConformanceCase {
                    suite: "fs-qty/conformance".into(),
                    case: "qty-001".into(),
                    pass: false,
                    detail: "value mismatch: got 0.13, want 0.12".into(),
                    seed: 7,
                },
                Some(3_000),
            ),
            em.emit(
                Severity::Info,
                EventKind::BenchmarkResult {
                    kernel: "gemm".into(),
                    metric: "gflops".into(),
                    value: 123.5,
                    machine: 0x1234,
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::StormAssertion {
                    name: "no-arena-leak".into(),
                    pass: true,
                    seed: 99,
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::Custom {
                    name: "regime-report".into(),
                    json: r#"{"re":100.5,"we":0.3}"#.into(),
                },
                None,
            ),
        ]
    }

    fn read_identity_len(bytes: &[u8], cursor: &mut usize) -> usize {
        let end = cursor
            .checked_add(core::mem::size_of::<u64>())
            .expect("test identity length offset fits usize");
        let encoded: [u8; 8] = bytes[*cursor..end]
            .try_into()
            .expect("identity length has eight bytes");
        *cursor = end;
        usize::try_from(u64::from_le_bytes(encoded)).expect("test identity length fits usize")
    }

    fn identity_field<'a>(canonical: &'a [u8], wanted: &str) -> (u8, &'a [u8]) {
        assert!(
            canonical.starts_with(ident::REPLAY_IDENTITY_DOMAIN.as_bytes()),
            "identity frame must start with its declared replay domain"
        );
        let mut cursor = ident::REPLAY_IDENTITY_DOMAIN.len() + core::mem::size_of::<u32>();
        let kind_len = read_identity_len(canonical, &mut cursor);
        cursor += kind_len;
        while cursor < canonical.len() {
            let tag = canonical[cursor];
            cursor += 1;
            let key_len = read_identity_len(canonical, &mut cursor);
            let key_end = cursor + key_len;
            let key = core::str::from_utf8(&canonical[cursor..key_end])
                .expect("identity field keys are UTF-8");
            cursor = key_end;
            let value_len = read_identity_len(canonical, &mut cursor);
            let value_end = cursor + value_len;
            if key == wanted {
                return (tag, &canonical[cursor..value_end]);
            }
            cursor = value_end;
        }
        panic!("identity field {wanted:?} was not encoded");
    }

    fn event_with_kind(kind: EventKind) -> Event {
        Event {
            session: "identity-session".into(),
            scope: "identity-scope".into(),
            seq: 17,
            severity: Severity::Info,
            kind,
            wall_ns: None,
        }
    }

    fn assert_payload_mutations(
        base: EventKind,
        mutations: Vec<(&'static str, EventKind)>,
        observed: &mut Vec<&'static str>,
    ) {
        let base_root = event_with_kind(base).content_hash();
        for (field, mutation) in mutations {
            assert_ne!(
                event_with_kind(mutation).content_hash(),
                base_root,
                "mutating {field} must move the exact event identity"
            );
            observed.push(field);
        }
    }

    #[test]
    fn every_kind_serializes_and_validates() {
        for e in sample_events() {
            let line = e.to_jsonl();
            validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        }
    }

    #[test]
    fn wall_clock_is_envelope_only() {
        let mut a = sample_events().remove(0);
        let h1 = a.content_hash();
        a.wall_ns = Some(999_999_999);
        let h2 = a.content_hash();
        assert_eq!(h1, h2, "content hash must exclude wall clock");
        // ...but the serialized line DOES carry it, as the last field.
        assert!(a.to_jsonl().ends_with(",\"wall_ns\":999999999}"));
    }

    #[test]
    fn content_hash_is_sensitive_to_content() {
        let events = sample_events();
        let mut hashes: Vec<u64> = events.iter().map(Event::content_hash).collect();
        hashes.sort_unstable();
        hashes.dedup();
        assert_eq!(
            hashes.len(),
            events.len(),
            "distinct events must hash distinctly"
        );
    }

    #[test]
    fn content_identity_preserves_bits_that_display_json_collapses() {
        let event = |residual| Event {
            session: "same-session".into(),
            scope: "same-scope".into(),
            seq: 7,
            severity: Severity::Info,
            kind: EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 3,
                residual,
            },
            wall_ns: None,
        };
        let first = event(f64::from_bits(0x7ff8_0000_0000_0001));
        let second = event(f64::from_bits(0x7ff8_0000_0000_0002));
        assert_eq!(
            first.to_jsonl(),
            second.to_jsonl(),
            "tagged display JSON intentionally does not expose NaN payload bits"
        );
        assert_ne!(first.content_hash(), second.content_hash());
        assert_ne!(
            first.content_identity().canonical_bytes(),
            second.content_identity().canonical_bytes()
        );
        assert_eq!(
            first.content_hash(),
            fnv1a64(first.content_identity().canonical_bytes()),
            "the stored root is derived from the exact canonical bytes"
        );
    }

    #[test]
    fn every_event_kind_payload_field_moves_identity() {
        let mut observed = Vec::new();

        assert_payload_mutations(
            EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 3,
                residual: 0.25,
            },
            vec![
                (
                    "solver_residual.solver",
                    EventKind::SolverResidual {
                        solver: "gmres".into(),
                        iter: 3,
                        residual: 0.25,
                    },
                ),
                (
                    "solver_residual.iter",
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 4,
                        residual: 0.25,
                    },
                ),
                (
                    "solver_residual.residual",
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 3,
                        residual: 0.5,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::TileComplete {
                tile: 7,
                kernel: "lbm".into(),
            },
            vec![
                (
                    "tile_complete.tile",
                    EventKind::TileComplete {
                        tile: 8,
                        kernel: "lbm".into(),
                    },
                ),
                (
                    "tile_complete.kernel",
                    EventKind::TileComplete {
                        tile: 7,
                        kernel: "gemm".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::Cancellation {
                reason: "budget".into(),
            },
            vec![(
                "cancellation.reason",
                EventKind::Cancellation {
                    reason: "panic".into(),
                },
            )],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::BudgetDelta {
                resource: "wall_s".into(),
                spent: 1.0,
                remaining: 9.0,
            },
            vec![
                (
                    "budget_delta.resource",
                    EventKind::BudgetDelta {
                        resource: "mem_bytes".into(),
                        spent: 1.0,
                        remaining: 9.0,
                    },
                ),
                (
                    "budget_delta.spent",
                    EventKind::BudgetDelta {
                        resource: "wall_s".into(),
                        spent: 2.0,
                        remaining: 9.0,
                    },
                ),
                (
                    "budget_delta.remaining",
                    EventKind::BudgetDelta {
                        resource: "wall_s".into(),
                        spent: 1.0,
                        remaining: 8.0,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 0.125,
                pass: true,
            },
            vec![
                (
                    "gradient_check.op",
                    EventKind::GradientCheck {
                        op: "elasticity".into(),
                        max_rel_err: 0.125,
                        pass: true,
                    },
                ),
                (
                    "gradient_check.max_rel_err",
                    EventKind::GradientCheck {
                        op: "poisson".into(),
                        max_rel_err: 0.25,
                        pass: true,
                    },
                ),
                (
                    "gradient_check.pass",
                    EventKind::GradientCheck {
                        op: "poisson".into(),
                        max_rel_err: 0.125,
                        pass: false,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ConformanceCase {
                suite: "suite-a".into(),
                case: "case-a".into(),
                pass: true,
                detail: "ok".into(),
                seed: 11,
            },
            vec![
                (
                    "conformance_case.suite",
                    EventKind::ConformanceCase {
                        suite: "suite-b".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.case",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-b".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.pass",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: false,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.detail",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "different".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.seed",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 12,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::BenchmarkResult {
                kernel: "gemm".into(),
                metric: "gflops".into(),
                value: 100.0,
                machine: 7,
            },
            vec![
                (
                    "benchmark_result.kernel",
                    EventKind::BenchmarkResult {
                        kernel: "spmv".into(),
                        metric: "gflops".into(),
                        value: 100.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.metric",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "bandwidth_gbs".into(),
                        value: 100.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.value",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "gflops".into(),
                        value: 101.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.machine",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "gflops".into(),
                        value: 100.0,
                        machine: 8,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: true,
                seed: 23,
            },
            vec![
                (
                    "storm_assertion.name",
                    EventKind::StormAssertion {
                        name: "cancel-latency".into(),
                        pass: true,
                        seed: 23,
                    },
                ),
                (
                    "storm_assertion.pass",
                    EventKind::StormAssertion {
                        name: "no-leak".into(),
                        pass: false,
                        seed: 23,
                    },
                ),
                (
                    "storm_assertion.seed",
                    EventKind::StormAssertion {
                        name: "no-leak".into(),
                        pass: true,
                        seed: 24,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ManifestSelection {
                manifest: "ab12".into(),
                source_snapshot: "cd34".into(),
            },
            vec![
                (
                    "manifest_selection.manifest",
                    EventKind::ManifestSelection {
                        manifest: "ef56".into(),
                        source_snapshot: "cd34".into(),
                    },
                ),
                (
                    "manifest_selection.source_snapshot",
                    EventKind::ManifestSelection {
                        manifest: "ab12".into(),
                        source_snapshot: "0789".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::StratumExpansion {
                stratum: "core".into(),
                profile: "nightly".into(),
                cases: 40,
            },
            vec![
                (
                    "stratum_expansion.stratum",
                    EventKind::StratumExpansion {
                        stratum: "max".into(),
                        profile: "nightly".into(),
                        cases: 40,
                    },
                ),
                (
                    "stratum_expansion.profile",
                    EventKind::StratumExpansion {
                        stratum: "core".into(),
                        profile: "release".into(),
                        cases: 40,
                    },
                ),
                (
                    "stratum_expansion.cases",
                    EventKind::StratumExpansion {
                        stratum: "core".into(),
                        profile: "nightly".into(),
                        cases: 41,
                    },
                ),
            ],
            &mut observed,
        );
        {
            let base_outcome = || ScopedReceiptOutcome {
                scope: ReceiptScope::Job,
                receipt: "ab12".into(),
                disposition: ExecutionDisposition::Completed,
                predicate: PredicateOutcome::Satisfied,
                evidence_methods: "oracle,enclosure".into(),
                grade: EpistemicGrade::Verified,
                applicability: DomainApplicability::InDomain,
                support: OperationalSupport::Supported,
                completeness: EvidenceCompleteness::Complete,
                integrity: EvidenceIntegrity::Intact,
                promotion: PromotionEffect::Unchanged,
                detail: "ok".into(),
            };
            let with = |mutate: &dyn Fn(&mut ScopedReceiptOutcome)| {
                let mut o = base_outcome();
                mutate(&mut o);
                EventKind::DsrRun {
                    run: "run-1".into(),
                    outcome: o,
                }
            };
            assert_payload_mutations(
                EventKind::DsrRun {
                    run: "run-1".into(),
                    outcome: base_outcome(),
                },
                vec![
                    (
                        "dsr_run.run",
                        EventKind::DsrRun {
                            run: "run-2".into(),
                            outcome: base_outcome(),
                        },
                    ),
                    (
                        "dsr_run.outcome.scope",
                        with(&|o| o.scope = ReceiptScope::Attempt),
                    ),
                    (
                        "dsr_run.outcome.receipt",
                        with(&|o| o.receipt = "cd34".into()),
                    ),
                    (
                        "dsr_run.outcome.disposition",
                        with(&|o| o.disposition = ExecutionDisposition::Failed),
                    ),
                    (
                        "dsr_run.outcome.predicate",
                        with(&|o| o.predicate = PredicateOutcome::Refuted),
                    ),
                    (
                        "dsr_run.outcome.evidence_methods",
                        with(&|o| o.evidence_methods = "oracle".into()),
                    ),
                    (
                        "dsr_run.outcome.grade",
                        with(&|o| o.grade = EpistemicGrade::Estimated),
                    ),
                    (
                        "dsr_run.outcome.applicability",
                        with(&|o| o.applicability = DomainApplicability::OutOfDomain),
                    ),
                    (
                        "dsr_run.outcome.support",
                        with(&|o| o.support = OperationalSupport::Degraded),
                    ),
                    (
                        "dsr_run.outcome.completeness",
                        with(&|o| o.completeness = EvidenceCompleteness::Incomplete),
                    ),
                    (
                        "dsr_run.outcome.integrity",
                        with(&|o| o.integrity = EvidenceIntegrity::Unchecked),
                    ),
                    (
                        "dsr_run.outcome.promotion",
                        with(&|o| o.promotion = PromotionEffect::Promoted),
                    ),
                    (
                        "dsr_run.outcome.detail",
                        with(&|o| o.detail = "different".into()),
                    ),
                ],
                &mut observed,
            );
        }
        assert_payload_mutations(
            EventKind::RaceRecord {
                resource: "arena-slot".into(),
                schedule: "a-b-a".into(),
                pass: true,
                seed: 31,
            },
            vec![
                (
                    "race_record.resource",
                    EventKind::RaceRecord {
                        resource: "tune-row".into(),
                        schedule: "a-b-a".into(),
                        pass: true,
                        seed: 31,
                    },
                ),
                (
                    "race_record.schedule",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "b-a-b".into(),
                        pass: true,
                        seed: 31,
                    },
                ),
                (
                    "race_record.pass",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "a-b-a".into(),
                        pass: false,
                        seed: 31,
                    },
                ),
                (
                    "race_record.seed",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "a-b-a".into(),
                        pass: true,
                        seed: 32,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::DegradationEvent {
                resource: "retained_bytes_per_scope".into(),
                limit: 4096,
                observed: 5000,
                action: "rotate-lane".into(),
            },
            vec![
                (
                    "degradation_event.resource",
                    EventKind::DegradationEvent {
                        resource: "events_per_scope".into(),
                        limit: 4096,
                        observed: 5000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.limit",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 8192,
                        observed: 5000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.observed",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 4096,
                        observed: 6000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.action",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 4096,
                        observed: 5000,
                        action: "refuse".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ImportReceipt {
                format: "step".into(),
                artifact: "ab12".into(),
                accepted: true,
                detail: "clean".into(),
            },
            vec![
                (
                    "import_receipt.format",
                    EventKind::ImportReceipt {
                        format: "gltf".into(),
                        artifact: "ab12".into(),
                        accepted: true,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.artifact",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "cd34".into(),
                        accepted: true,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.accepted",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "ab12".into(),
                        accepted: false,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.detail",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "ab12".into(),
                        accepted: true,
                        detail: "quarantined".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::CertificateVerdict {
                certificate: "ivl-enclosure".into(),
                pass: true,
                bound: 0.125,
                detail: "ok".into(),
            },
            vec![
                (
                    "certificate_verdict.certificate",
                    EventKind::CertificateVerdict {
                        certificate: "trim-cert".into(),
                        pass: true,
                        bound: 0.125,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.pass",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: false,
                        bound: 0.125,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.bound",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: true,
                        bound: 0.25,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.detail",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: true,
                        bound: 0.125,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::Custom {
                name: "opaque".into(),
                json: r#"{"a":1,"b":2}"#.into(),
            },
            vec![
                (
                    "custom.name",
                    EventKind::Custom {
                        name: "other".into(),
                        json: r#"{"a":1,"b":2}"#.into(),
                    },
                ),
                (
                    "custom.json",
                    EventKind::Custom {
                        name: "opaque".into(),
                        json: r#"{"b":2,"a":1}"#.into(),
                    },
                ),
            ],
            &mut observed,
        );

        let expected = [
            "solver_residual.solver",
            "solver_residual.iter",
            "solver_residual.residual",
            "tile_complete.tile",
            "tile_complete.kernel",
            "cancellation.reason",
            "budget_delta.resource",
            "budget_delta.spent",
            "budget_delta.remaining",
            "gradient_check.op",
            "gradient_check.max_rel_err",
            "gradient_check.pass",
            "conformance_case.suite",
            "conformance_case.case",
            "conformance_case.pass",
            "conformance_case.detail",
            "conformance_case.seed",
            "benchmark_result.kernel",
            "benchmark_result.metric",
            "benchmark_result.value",
            "benchmark_result.machine",
            "storm_assertion.name",
            "storm_assertion.pass",
            "storm_assertion.seed",
            "manifest_selection.manifest",
            "manifest_selection.source_snapshot",
            "stratum_expansion.stratum",
            "stratum_expansion.profile",
            "stratum_expansion.cases",
            "dsr_run.run",
            "dsr_run.outcome.scope",
            "dsr_run.outcome.receipt",
            "dsr_run.outcome.disposition",
            "dsr_run.outcome.predicate",
            "dsr_run.outcome.evidence_methods",
            "dsr_run.outcome.grade",
            "dsr_run.outcome.applicability",
            "dsr_run.outcome.support",
            "dsr_run.outcome.completeness",
            "dsr_run.outcome.integrity",
            "dsr_run.outcome.promotion",
            "dsr_run.outcome.detail",
            "race_record.resource",
            "race_record.schedule",
            "race_record.pass",
            "race_record.seed",
            "degradation_event.resource",
            "degradation_event.limit",
            "degradation_event.observed",
            "degradation_event.action",
            "import_receipt.format",
            "import_receipt.artifact",
            "import_receipt.accepted",
            "import_receipt.detail",
            "certificate_verdict.certificate",
            "certificate_verdict.pass",
            "certificate_verdict.bound",
            "certificate_verdict.detail",
            "custom.name",
            "custom.json",
        ];
        assert_eq!(
            observed.as_slice(),
            expected.as_slice(),
            "all 60 payload fields stay enumerated"
        );
    }

    #[test]
    fn scope_tree_mirrors_enter_exit_and_refuses_forgery_and_imbalance() {
        let mut em = Emitter::new("study-x", "op-1");
        assert_eq!(em.current_scope(), "op-1");
        assert_eq!(em.scope_depth(), 0);

        em.enter_scope("kernel-cg").expect("child scope");
        em.enter_scope("tile-7").expect("grandchild scope");
        assert_eq!(em.current_scope(), "op-1/kernel-cg/tile-7");
        assert_eq!(em.scope_depth(), 2);
        let deep = em.emit(
            Severity::Info,
            EventKind::TileComplete {
                tile: 7,
                kernel: "cg".into(),
            },
            None,
        );
        assert_eq!(deep.scope, "op-1/kernel-cg/tile-7");

        em.exit_scope().expect("exit grandchild");
        assert_eq!(em.current_scope(), "op-1/kernel-cg");
        let mid = em.emit(
            Severity::Info,
            EventKind::Cancellation {
                reason: "budget".into(),
            },
            None,
        );
        assert_eq!(mid.scope, "op-1/kernel-cg");
        assert_eq!(
            mid.seq,
            deep.seq + 1,
            "seq is one monotone stream across scope moves"
        );

        em.exit_scope().expect("exit child");
        assert_eq!(em.current_scope(), "op-1");
        assert!(
            em.exit_scope().is_err(),
            "exit at the base scope must refuse (unbalanced enter/exit)"
        );

        assert!(em.enter_scope("").is_err(), "empty segment must refuse");
        assert!(
            em.enter_scope("a/b").is_err(),
            "path separator inside a segment is scope forgery"
        );
        assert!(
            em.enter_scope("a\u{7}b").is_err(),
            "control characters inside a segment must refuse"
        );
        assert_eq!(
            em.current_scope(),
            "op-1",
            "refused entries must not move the scope"
        );
    }

    #[test]
    fn scope_tree_replay_is_deterministic_and_wire_valid() {
        let run = || {
            let mut em = Emitter::new("study-x", "op-1");
            let mut lines = Vec::new();
            em.enter_scope("solve").expect("scope");
            lines.push(
                em.emit(
                    Severity::Info,
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 0,
                        residual: 0.5,
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.enter_scope("tile-3").expect("scope");
            lines.push(
                em.emit(
                    Severity::Trace,
                    EventKind::TileComplete {
                        tile: 3,
                        kernel: "cg".into(),
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.exit_scope().expect("exit");
            em.exit_scope().expect("exit");
            lines
        };
        let first = run();
        let second = run();
        assert_eq!(
            first, second,
            "identical scope walks replay bit-identically"
        );
        for line in &first {
            validate_line(line).expect("scope-tree events stay wire-valid");
        }
    }

    #[test]
    fn verification_outcome_kinds_are_wire_valid_and_lint_enforced() {
        let mut em = Emitter::new("study-x", "dsr-1");
        let selection = em.emit(
            Severity::Info,
            EventKind::ManifestSelection {
                manifest: "ab12".into(),
                source_snapshot: "cd34".into(),
            },
            None,
        );
        let expansion = em.emit(
            Severity::Info,
            EventKind::StratumExpansion {
                stratum: "core".into(),
                profile: "nightly".into(),
                cases: 40,
            },
            None,
        );
        let mut outcome = ScopedReceiptOutcome {
            scope: ReceiptScope::Job,
            receipt: "ab12".into(),
            disposition: ExecutionDisposition::Completed,
            predicate: PredicateOutcome::Satisfied,
            evidence_methods: "oracle,enclosure".into(),
            grade: EpistemicGrade::Verified,
            applicability: DomainApplicability::InDomain,
            support: OperationalSupport::Supported,
            completeness: EvidenceCompleteness::Complete,
            integrity: EvidenceIntegrity::Intact,
            promotion: PromotionEffect::Unchanged,
            detail: "ok".into(),
        };
        let run = em.emit(
            Severity::Info,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome: outcome.clone(),
            },
            Some(9),
        );
        for event in [&selection, &expansion, &run] {
            validate_line(&event.to_jsonl()).expect("verification kinds stay wire-valid");
            lint_failure_record(event).expect("conclusive outcomes pass the lint");
        }

        // A failed outcome without diagnostic detail refuses.
        outcome.disposition = ExecutionDisposition::Failed;
        outcome.detail = String::new();
        let bad = em.emit(
            Severity::Error,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome: outcome.clone(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad).is_err(),
            "failed outcome without detail must refuse"
        );
        outcome.detail = "worker lost lease".into();
        let repaired = em.emit(
            Severity::Error,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome,
            },
            None,
        );
        lint_failure_record(&repaired).expect("detailed failure passes");
        validate_line(&repaired.to_jsonl()).expect("failure line stays wire-valid");
    }

    #[test]
    fn custom_payload_identity_is_exact_opaque_utf8() {
        let opaque = r#"{ "b":2, "a":1 }"#;
        let event = event_with_kind(EventKind::Custom {
            name: "opaque-json".into(),
            json: opaque.into(),
        });
        let identity = event.content_identity();
        let (tag, retained) = identity_field(identity.canonical_bytes(), "custom_json_opaque_utf8");
        assert_eq!(
            tag, 0x05,
            "custom JSON is bound as exact bytes, not text claiming canonical JSON"
        );
        assert_eq!(
            retained,
            opaque.as_bytes(),
            "opaque UTF-8 round-trips byte-for-byte through the identity frame"
        );
        assert!(event.to_jsonl().contains(&format!("\"data\":{opaque}")));
        validate_line(&event.to_jsonl())
            .expect("the valid opaque object remains a valid event line");

        let reordered = event_with_kind(EventKind::Custom {
            name: "opaque-json".into(),
            json: r#"{"a":1,"b":2}"#.into(),
        });
        assert_ne!(
            event.content_hash(),
            reordered.content_hash(),
            "whitespace and member order are honestly semantic under opaque-byte identity"
        );
    }

    #[test]
    fn event_content_identity_mutation_battery() {
        let base = Event {
            session: "session-a".into(),
            scope: "scope-a".into(),
            seq: 11,
            severity: Severity::Info,
            kind: EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 0.25,
                pass: true,
            },
            wall_ns: None,
        };
        let base_hash = base.content_hash();
        let mutations = [
            Event {
                session: "session-b".into(),
                ..base.clone()
            },
            Event {
                scope: "scope-b".into(),
                ..base.clone()
            },
            Event {
                seq: 12,
                ..base.clone()
            },
            Event {
                severity: Severity::Warn,
                ..base.clone()
            },
            Event {
                kind: EventKind::Cancellation {
                    reason: "budget".into(),
                },
                ..base.clone()
            },
            Event {
                kind: EventKind::GradientCheck {
                    op: "poisson".into(),
                    max_rel_err: 0.25_f64.next_up(),
                    pass: true,
                },
                ..base.clone()
            },
        ];
        assert!(
            mutations
                .iter()
                .all(|mutation| mutation.content_hash() != base_hash),
            "every mutable semantic event field must move the typed identity"
        );
        let mut envelope = base;
        envelope.wall_ns = Some(u64::MAX);
        assert_eq!(
            envelope.content_hash(),
            base_hash,
            "the declared wall-clock exclusion must not move identity"
        );
    }

    #[test]
    fn event_content_identity_domain_version_and_wire_schema_bytes_are_independent() {
        let event = event_with_kind(EventKind::Cancellation {
            reason: "budget".into(),
        });
        let current =
            event.content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION);
        let domain_mutation = event.content_identity_with_schema(
            "org.frankensim.fs-obs.event-content.v2.alternate",
            EVENT_CONTENT_IDENTITY_VERSION,
            SCHEMA_VERSION,
        );
        let identity_version_mutation = event
            .content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION + 1, SCHEMA_VERSION);
        let wire_schema_mutation = event
            .content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION + 1);

        let (identity_tag, current_identity_version) =
            identity_field(current.canonical_bytes(), "event_content_identity_version");
        let (schema_tag, current_wire_schema) =
            identity_field(current.canonical_bytes(), "event_wire_schema_version");
        assert_eq!(identity_tag, 0x02);
        assert_eq!(schema_tag, 0x02);
        // Canonical bytes carry the shared replay-identity framing —
        // "fsid" + u32 ident-schema version + u64 kind length — and bind the
        // event artifact domain as the length-framed identity kind, so the
        // domain begins after the fixed 16-byte frame header rather than at
        // byte zero (length framing is what makes prefix collisions
        // unrepresentable; see ident_003).
        const FRAME_HEADER_LEN: usize = ident::REPLAY_IDENTITY_DOMAIN.len()
            + core::mem::size_of::<u32>()
            + core::mem::size_of::<u64>();
        assert!(
            current
                .canonical_bytes()
                .starts_with(ident::REPLAY_IDENTITY_DOMAIN.as_bytes())
        );
        assert!(
            current.canonical_bytes()[FRAME_HEADER_LEN..]
                .starts_with(EVENT_CONTENT_IDENTITY_DOMAIN.as_bytes())
        );
        assert!(
            domain_mutation.canonical_bytes()[FRAME_HEADER_LEN..]
                .starts_with(b"org.frankensim.fs-obs.event-content.v2.alternate")
        );
        assert_eq!(
            current_identity_version,
            u64::from(EVENT_CONTENT_IDENTITY_VERSION)
                .to_le_bytes()
                .as_slice()
        );
        assert_eq!(
            current_wire_schema,
            u64::from(SCHEMA_VERSION).to_le_bytes().as_slice()
        );
        assert_eq!(
            identity_field(
                domain_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version,
            "artifact-domain mutation must leave identity-version bytes unchanged"
        );
        assert_eq!(
            identity_field(
                domain_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema,
            "artifact-domain mutation must leave wire-schema bytes unchanged"
        );

        assert_eq!(
            identity_field(
                identity_version_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema,
            "identity-version mutation must leave the wire-schema bytes unchanged"
        );
        assert_eq!(
            identity_field(
                wire_schema_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version,
            "wire-schema mutation must leave the identity-version bytes unchanged"
        );
        assert_ne!(
            identity_field(
                identity_version_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version
        );
        assert_ne!(
            identity_field(
                wire_schema_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema
        );
        assert_ne!(
            current.canonical_bytes(),
            identity_version_mutation.canonical_bytes()
        );
        assert_ne!(current.canonical_bytes(), domain_mutation.canonical_bytes());
        assert_ne!(current.root(), domain_mutation.root());
        assert_ne!(current.root(), identity_version_mutation.root());
        assert_ne!(
            current.canonical_bytes(),
            wire_schema_mutation.canonical_bytes()
        );
        assert_ne!(current.root(), wire_schema_mutation.root());
    }

    #[test]
    fn retained_event_identity_receipts_admit_exactly_or_fail_closed() {
        let event = event_with_kind(EventKind::StormAssertion {
            name: "no-leak".into(),
            pass: true,
            seed: 23,
        });
        let captured = event.content_identity_receipt();
        assert_eq!(
            captured.declared_identity_version(),
            EVENT_CONTENT_IDENTITY_VERSION
        );
        assert_eq!(captured.root(), fnv1a64(captured.canonical_bytes()));

        let retained = EventIdentityReceipt::from_retained_parts(
            captured.declared_identity_version(),
            captured.canonical_bytes().to_vec(),
            captured.root(),
        );
        assert_eq!(event.admit_content_identity(&retained), Ok(()));

        let stale = EventIdentityReceipt::from_retained_parts(
            0,
            captured.canonical_bytes().to_vec(),
            captured.root(),
        );
        assert_eq!(
            event.admit_content_identity(&stale),
            Err(EventIdentityAdmissionError::UnsupportedVersion(
                EventIdentityVersionError {
                    declared: 0,
                    supported: EVENT_CONTENT_IDENTITY_VERSION,
                }
            ))
        );

        let wrong_root = EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            captured.canonical_bytes().to_vec(),
            captured.root() ^ 1,
        );
        assert!(matches!(
            event.admit_content_identity(&wrong_root),
            Err(EventIdentityAdmissionError::RootMismatch { .. })
        ));

        let mut foreign_bytes = captured.canonical_bytes().to_vec();
        let last = foreign_bytes
            .last_mut()
            .expect("event identity canonical bytes are non-empty");
        *last ^= 1;
        let foreign_root = fnv1a64(&foreign_bytes);
        let self_consistent_but_foreign = EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            foreign_bytes,
            foreign_root,
        );
        assert_eq!(
            event.admit_content_identity(&self_consistent_but_foreign),
            Err(EventIdentityAdmissionError::CanonicalBytesMismatch {
                declared_root: foreign_root,
                expected_root: captured.root(),
            })
        );
    }

    #[test]
    fn event_content_identity_versions_fail_closed() {
        assert!(check_event_content_identity_version(EVENT_CONTENT_IDENTITY_VERSION).is_ok());
        for declared in [0, EVENT_CONTENT_IDENTITY_VERSION + 1] {
            assert_eq!(
                check_event_content_identity_version(declared),
                Err(EventIdentityVersionError {
                    declared,
                    supported: EVENT_CONTENT_IDENTITY_VERSION,
                })
            );
        }
    }

    #[test]
    fn sequence_numbers_are_monotone_per_emitter() {
        let events = sample_events();
        for (i, e) in events.iter().enumerate() {
            assert_eq!(e.seq, i as u64);
        }
    }

    #[test]
    fn golden_line_shape_is_stable() {
        // Schema evolution is additive-only; this golden line is the contract.
        // Changing it requires a SCHEMA_VERSION bump and a semantic justification
        // (golden-evidence policy, AGENTS.md).
        let e = Event {
            session: "s".into(),
            scope: "a/b".into(),
            seq: 5,
            severity: Severity::Info,
            kind: EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 1e-9,
                pass: true,
            },
            wall_ns: None,
        };
        // Note: Rust's shortest-round-trip float Display never uses scientific
        // notation, so 1e-9 serializes as 0.000000001 — that IS the contract.
        assert_eq!(
            e.to_jsonl(),
            "{\"v\":1,\"session\":\"s\",\"scope\":\"a/b\",\"seq\":5,\"severity\":\"info\",\
             \"kind\":\"gradient_check\",\"payload\":{\"op\":\"poisson\",\
             \"max_rel_err\":0.000000001,\"pass\":true}}"
        );
    }

    #[test]
    fn validator_rejects_corruption() {
        let good = sample_events()[0].to_jsonl();
        for bad in [
            String::new(),
            "not json".to_string(),
            good.replace("\"v\":1", "\"v\":99"),
            good.replace("solver_residual", "mystery_kind"),
            good.replace("\"session\"", "\"sesion\""),
            good[..good.len() - 1].to_string(),
        ] {
            assert!(validate_line(&bad).is_err(), "should reject: {bad}");
        }
    }

    #[test]
    fn failure_lint_demands_reproduction_ingredients() {
        let mut em = Emitter::new("s", "x");
        let bad = em.emit(
            Severity::Error,
            EventKind::ConformanceCase {
                suite: "s".into(),
                case: "c".into(),
                pass: false,
                detail: String::new(),
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad).is_err());
        let good = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 42,
            },
            None,
        );
        assert!(lint_failure_record(&good).is_ok());
        let bad_storm = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad_storm).is_err());
    }

    #[test]
    fn escaping_handles_hostile_strings() {
        let mut em = Emitter::new("s\"es\\sion\n", "sc\tope");
        let e = em.emit(
            Severity::Error,
            EventKind::Cancellation {
                reason: "quote\" backslash\\ newline\n tab\t".into(),
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        assert!(!line.contains('\n'), "JSONL lines must be single-line");
    }

    #[test]
    fn non_finite_floats_are_tagged_not_invalid() {
        let mut em = Emitter::new("s", "x");
        let e = em.emit(
            Severity::Info,
            EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 1,
                residual: f64::NAN,
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).expect("tagged non-finite must stay valid");
        assert!(line.contains("non-finite:NaN"));
    }

    #[test]
    fn fnv_matches_known_answers() {
        // Published FNV-1a 64 test vectors.
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x85944171f73967e8);
    }
}
