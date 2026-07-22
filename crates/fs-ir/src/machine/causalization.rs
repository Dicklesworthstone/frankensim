//! Versioned equation-variable hypergraph and causalization-receipt schemas.
//!
//! This module is the I02.1 structural boundary between admitted Machine IR
//! and later equation extraction, matching, index reduction, and execution.
//! It records what an analyzer is allowed to consume and what a future
//! analyzer is allowed to report; it deliberately implements none of those
//! algorithms. In particular, an admitted graph or receipt does not prove
//! numerical rank, numerical nonsingularity, physical cause and effect,
//! solvability, a DAE index, or correctness of an opaque operator reference.
//!
//! Every equation and variable has a nominal, domain-separated identity
//! derived from source lineage rather than a display name or vector position.
//! Graph identity excludes diagnostic labels, canonicalizes all collections,
//! and binds the exact admitted Machine graph. Opaque sources are admitted only
//! through an explicit source identity plus a separate audit-receipt identity.

use core::fmt;
use core::hash::{Hash, Hasher};
use core::num::NonZeroU64;

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use fs_blake3::identity::{
    CancellationProbe, CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema,
    ChildSpec, EntityId, EvidenceNodeId, Field, FieldSpec, IdentityReceipt, LimitKind, NeverCancel,
    OrderedBytesStreamDiagnostic, OrderedBytesStreamError, ProblemSemanticId, StrongIdentity,
    WireType,
};
use fs_exec::Cx;
use fs_qty::{DIMENSION_COUNT, Dims};

use super::semantics::{AdmittedMachineBehavior, MachineBehaviorIdV1, StateSlotContract};
use super::{
    AdmittedMachineGraph, ClockId, FrameBinding, InterfaceId, MachineElementId, MachineGraphIdV1,
    MachineReferenceError, PortId, RelationId, StateSlotId, SubsystemId, TerminalId,
    TerminalQuantitySpec, TerminalShape,
};

/// Shared candidate version for equation, variable, and incidence entity-ID
/// auxiliaries; this is not the top-level causal-structure identity version.
pub const CAUSAL_GRAPH_SCHEMA_VERSION_V1: u32 = 1;
/// Shared candidate version for matching-set and conditional-outcome-set
/// auxiliaries; this is not the complete causalization-receipt version.
pub const CAUSALIZATION_RECEIPT_SCHEMA_VERSION_V1: u32 = 1;
/// Schema version of normalized causal-structure identities.
pub const CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1: u32 = 1;
/// Schema version of provenance-bearing causal-graph artifacts.
pub const CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1: u32 = 1;
/// Schema version of producer-independent causal outcomes.
pub const CAUSAL_OUTCOME_IDENTITY_SCHEMA_VERSION_V1: u32 = 1;
/// Schema version of complete causalization evidence receipts.
pub const CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1: u32 = 1;
/// Domain of normalized equation-variable structure identities.
pub const CAUSAL_STRUCTURE_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-ir.machine.causal-structure.v1";
/// Domain of provenance-bearing causal graph artifacts.
pub const CAUSAL_GRAPH_ARTIFACT_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-ir.machine.causal-graph-artifact.v1";
/// Domain of producer-independent normalized causal outcomes.
pub const CAUSAL_OUTCOME_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-ir.machine.causal-outcome.v1";
/// Domain of complete provenance-bearing causalization receipts.
pub const CAUSALIZATION_RECEIPT_IDENTITY_DOMAIN_V1: &str =
    "org.frankensim.fs-ir.machine.causalization-receipt.v1";

/// Maximum equations in one graph draft.
pub const MAX_CAUSAL_EQUATIONS: usize = 65_536;
/// Maximum variables in one graph draft.
pub const MAX_CAUSAL_VARIABLES: usize = 65_536;
/// Maximum structural incidences in one graph draft.
pub const MAX_CAUSAL_INCIDENCES: usize = 1_048_576;
/// Maximum distinct activation-condition tables.
pub const MAX_CAUSAL_CONDITIONS: usize = 4_096;
/// Maximum aggregate activation dependencies.
pub const MAX_CAUSAL_CONDITION_DEPENDENCIES: usize = 65_536;
/// Maximum aggregate branches/modes across activation conditions.
pub const MAX_CAUSAL_CONDITION_BRANCHES: usize = 65_536;
/// Maximum aggregate activation-conjunction selections across graph rows.
pub const MAX_CAUSAL_ACTIVATION_SELECTIONS: usize = 262_144;
/// Maximum aggregate DNF cubes across graph rows.
pub const MAX_CAUSAL_ACTIVATION_CUBES: usize = 262_144;
/// Maximum cubes admitted in one DNF activation domain.
pub const MAX_CAUSAL_CUBES_PER_ACTIVATION: usize = 4_096;
/// Maximum selections admitted in one DNF cube.
pub const MAX_CAUSAL_SELECTIONS_PER_CUBE: usize = 4_096;
/// Maximum charged symbolic work units consumed by exact activation
/// implication, coverage, and overlap checks in one graph admission.
///
/// A work unit is charged for every explored partial assignment, cube
/// compatibility check, selection comparison, and finite-domain branch-table
/// lookup. This is intentionally an operation budget rather than merely a cap
/// on the number of materialized search states.
pub const MAX_CAUSAL_ACTIVATION_PROOF_STATES: usize = 1_048_576;
/// Maximum aggregate support references across all nodes.
pub const MAX_CAUSAL_SUPPORT_REFERENCES: usize = 262_144;
/// Maximum aggregate parent references across all derived node lineages.
pub const MAX_CAUSAL_DERIVATION_REFERENCES: usize = 262_144;
/// Maximum direct parents in one derived node lineage.
pub const MAX_CAUSAL_DERIVATION_PARENTS: usize = 4_096;
/// Maximum derivative order represented by the baseline structural schema.
pub const MAX_CAUSAL_DERIVATIVE_ORDER: u16 = 32;
/// Maximum UTF-8 bytes in a non-semantic diagnostic label.
pub const MAX_CAUSAL_DIAGNOSTIC_LABEL_BYTES: usize = 256;
/// Maximum matching pairs in one causalization receipt.
pub const MAX_CAUSAL_MATCHING_PAIRS: usize = 65_536;
/// Maximum distinct derivative-variable vertices representable by a receipt.
pub const MAX_CAUSAL_DERIVATIVE_VERTICES: usize = 131_072;
/// Maximum condition-specific child outcomes in one receipt.
pub const MAX_CAUSAL_CONDITIONAL_OUTCOMES: usize = 4_096;
/// Maximum receipt-wide condition-to-branch selections across a mode domain
/// or hybrid child outcomes (mutually exclusive in an admitted receipt).
pub const MAX_CAUSAL_CONDITIONAL_SELECTIONS: usize = 65_536;

const NODE_IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(16_384, 8_192, 1, 1, 256);
// One incidence may carry the complete aggregate activation-selection envelope:
// 262,144 pairs of independently namespaced condition/branch references are
// below 96 MiB at the Machine key bound. Keep the standalone receipt large
// enough to make the public semantic cap real rather than an 8 KiB accident.
const INCIDENCE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(128 * 1_024 * 1_024, 112 * 1_024 * 1_024, 1, 1, 256);
const IDENTITY_RECEIPT_ADJUDICATION_BYTES: usize = 3 * 32 + 8 + 4 + 8;
const CAUSAL_OWNER_CANONICAL_BYTES: usize = 1 + 32;
// A max-key audited bridge contributes a complete incidence adjudication tuple
// plus two canonical references to the provenance-only artifact row. Including
// row-length frames, all 1,048,576 publicly admitted incidences require just
// under 478 MiB in that ordered field. Keep the field and total limits above
// the declared graph cap instead of silently lowering it to an encoder accident.
const MAX_CAUSAL_REFERENCE_CANONICAL_BYTES: usize =
    8 + super::MAX_MACHINE_ENTITY_KEY_BYTES + 8 + 32;
const MAX_CAUSAL_INCIDENCE_ARTIFACT_ROW_BYTES: usize =
    IDENTITY_RECEIPT_ADJUDICATION_BYTES + 1 + 2 * MAX_CAUSAL_REFERENCE_CANONICAL_BYTES;
const MAX_CAUSAL_INCIDENCE_ARTIFACT_FIELD_BYTES: usize =
    8 + MAX_CAUSAL_INCIDENCES * (8 + MAX_CAUSAL_INCIDENCE_ARTIFACT_ROW_BYTES);
const CAUSAL_GRAPH_MAX_FIELD_BYTES: u64 = 512 * 1_024 * 1_024;
// Structure identity charges equation, variable, condition, typed-incidence,
// and incidence-adjudication collections independently at their public caps.
const CAUSAL_GRAPH_MAX_COLLECTION_ITEMS: u64 = 2_232_320;
#[allow(clippy::cast_possible_truncation)]
const _: () =
    assert!(MAX_CAUSAL_INCIDENCE_ARTIFACT_FIELD_BYTES <= CAUSAL_GRAPH_MAX_FIELD_BYTES as usize);
#[allow(clippy::cast_possible_truncation)]
const _: () = assert!(
    MAX_CAUSAL_EQUATIONS + MAX_CAUSAL_VARIABLES + MAX_CAUSAL_CONDITIONS + 2 * MAX_CAUSAL_INCIDENCES
        == CAUSAL_GRAPH_MAX_COLLECTION_ITEMS as usize
);
const CAUSAL_GRAPH_IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(
    768 * 1_024 * 1_024,
    CAUSAL_GRAPH_MAX_FIELD_BYTES,
    10,
    CAUSAL_GRAPH_MAX_COLLECTION_ITEMS,
    4_096,
);
const CAUSAL_RECEIPT_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(64 * 1_024 * 1_024, 32 * 1_024 * 1_024, 24, 300_000, 4_096);
const TIME_DERIVATIVE_AXIS: usize = 2;

/// Origin of one failure inside fallible ordered-row identity publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CausalIdentityStreamOrigin {
    /// The canonical encoder refused framing, limits, accounting, or
    /// cancellation.
    Canonical,
    /// The independent length planner or the admitted row producer refused.
    Producer,
    /// Exact row storage could not be reserved after row admission.
    Allocation,
}

/// Deterministic evidence retained when fallible causal-row storage refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CausalIdentityAllocationRefusal {
    requested_bytes: u64,
}

impl CausalIdentityAllocationRefusal {
    /// Exact admitted logical row bytes passed to `Vec::try_reserve_exact`.
    #[must_use]
    pub const fn requested_bytes(self) -> u64 {
        self.requested_bytes
    }
}

impl fmt::Display for CausalIdentityAllocationRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "causal identity row allocation refused {} admitted bytes",
            self.requested_bytes
        )
    }
}

impl std::error::Error for CausalIdentityAllocationRefusal {}

/// Stable source of a causal identity publication refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CausalIdentityFailure {
    /// Canonical framing, planning, production, or cancellation refused.
    Canonical(CanonicalError),
    /// Fallible storage reservation refused after exact row admission.
    Allocation(CausalIdentityAllocationRefusal),
}

impl fmt::Display for CausalIdentityFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Canonical(error) => error.fmt(f),
            Self::Allocation(refusal) => refusal.fmt(f),
        }
    }
}

/// Exact ordered-row progress retained by a causal identity refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CausalIdentityStreamRefusal {
    origin: CausalIdentityStreamOrigin,
    diagnostic: OrderedBytesStreamDiagnostic,
}

impl CausalIdentityStreamRefusal {
    /// Whether canonical work, caller production, or fallible row reservation
    /// originated the refusal.
    #[must_use]
    pub const fn origin(self) -> CausalIdentityStreamOrigin {
        self.origin
    }

    /// Schema, field, phase, row, completed-row, and byte progress at refusal.
    #[must_use]
    pub const fn diagnostic(self) -> OrderedBytesStreamDiagnostic {
        self.diagnostic
    }
}

/// Causal identity failure with optional ordered-row progress evidence.
///
/// `canonical_error` preserves canonical sources while `allocation_refusal`
/// preserves a distinct fallible-reservation source. When the failure occurred
/// in an ordered byte collection, `stream_refusal` additionally proves its
/// origin, exact phase, row, completed prefix, byte progress, and no-publication
/// disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CausalIdentityError {
    failure: CausalIdentityFailure,
    stream_refusal: Option<CausalIdentityStreamRefusal>,
}

impl CausalIdentityError {
    /// Stable canonical source, absent only for a distinct allocation refusal.
    #[must_use]
    pub const fn canonical_error(self) -> Option<CanonicalError> {
        match self.failure {
            CausalIdentityFailure::Canonical(error) => Some(error),
            CausalIdentityFailure::Allocation(_) => None,
        }
    }

    /// Exact fallible-reservation refusal, when row storage was unavailable.
    #[must_use]
    pub const fn allocation_refusal(self) -> Option<CausalIdentityAllocationRefusal> {
        match self.failure {
            CausalIdentityFailure::Canonical(_) => None,
            CausalIdentityFailure::Allocation(refusal) => Some(refusal),
        }
    }

    /// Ordered-row progress evidence, when the refusal came from that adapter.
    #[must_use]
    pub const fn stream_refusal(self) -> Option<CausalIdentityStreamRefusal> {
        self.stream_refusal
    }

    const fn is_cancelled(self) -> bool {
        matches!(
            self.failure,
            CausalIdentityFailure::Canonical(CanonicalError::Cancelled { .. })
        )
    }
}

impl From<CanonicalError> for CausalIdentityError {
    fn from(canonical_error: CanonicalError) -> Self {
        Self {
            failure: CausalIdentityFailure::Canonical(canonical_error),
            stream_refusal: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CausalIdentityRowError {
    Canonical(CanonicalError),
    Allocation(CausalIdentityAllocationRefusal),
}

impl From<CanonicalError> for CausalIdentityRowError {
    fn from(error: CanonicalError) -> Self {
        Self::Canonical(error)
    }
}

impl From<CausalIdentityAllocationRefusal> for CausalIdentityRowError {
    fn from(refusal: CausalIdentityAllocationRefusal) -> Self {
        Self::Allocation(refusal)
    }
}

impl From<OrderedBytesStreamError<CausalIdentityRowError>> for CausalIdentityError {
    fn from(error: OrderedBytesStreamError<CausalIdentityRowError>) -> Self {
        match error {
            OrderedBytesStreamError::Canonical { source, diagnostic } => Self {
                failure: CausalIdentityFailure::Canonical(source),
                stream_refusal: Some(CausalIdentityStreamRefusal {
                    origin: CausalIdentityStreamOrigin::Canonical,
                    diagnostic,
                }),
            },
            OrderedBytesStreamError::Producer { source, diagnostic } => match source {
                CausalIdentityRowError::Canonical(source) => Self {
                    failure: CausalIdentityFailure::Canonical(source),
                    stream_refusal: Some(CausalIdentityStreamRefusal {
                        origin: CausalIdentityStreamOrigin::Producer,
                        diagnostic,
                    }),
                },
                CausalIdentityRowError::Allocation(refusal) => Self {
                    failure: CausalIdentityFailure::Allocation(refusal),
                    stream_refusal: Some(CausalIdentityStreamRefusal {
                        origin: CausalIdentityStreamOrigin::Allocation,
                        diagnostic,
                    }),
                },
            },
        }
    }
}

impl fmt::Display for CausalIdentityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(stream) = self.stream_refusal {
            let diagnostic = stream.diagnostic();
            write!(
                f,
                "{} (ordered field `{}`, {:?} origin, phase {:?}, row {:?}, {}/{} rows complete)",
                self.failure,
                diagnostic.field_name(),
                stream.origin(),
                diagnostic.phase(),
                diagnostic.row_index(),
                diagnostic.completed_rows(),
                diagnostic.declared_rows(),
            )
        } else {
            self.failure.fmt(f)
        }
    }
}

impl std::error::Error for CausalIdentityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.failure {
            CausalIdentityFailure::Canonical(error) => Some(error),
            CausalIdentityFailure::Allocation(refusal) => Some(refusal),
        }
    }
}
const CAUSAL_CANCELLATION_POLL_STRIDE: usize = 256;
// Invalid drafts receive bounded, deterministic diagnostics rather than an
// attacker-controlled multi-million-row finding allocation. Crossing this
// budget yields the ordinary ResourceLimit sentinel and publishes no identity.
/// Maximum detailed findings retained before graph admission returns a single
/// deterministic [`CausalGraphRule::ResourceLimit`] sentinel.
pub const MAX_CAUSAL_GRAPH_FINDINGS: usize = 65_536;
/// Maximum detailed findings retained before receipt admission returns a
/// single deterministic [`CausalReceiptRule::ResourceLimit`] sentinel.
pub const MAX_CAUSAL_RECEIPT_FINDINGS: usize = 65_536;
const INCIDENCE_FIXED_CANONICAL_CAPACITY: usize = 1_024;
const NODE_FIXED_CANONICAL_CAPACITY: usize = 1_024;

// `IdentityReceipt::limits()` records the encoder policy that admitted a
// canonical preimage; it is evidence about construction, not part of that
// preimage's meaning. Collision adjudication therefore compares and composes
// every independently derived semantic/preimage/schema/metric coordinate while
// deliberately excluding only the limits record. Keep Eq, Ord, Hash, parent
// composition, and admission revalidation on this exact tuple.
fn identity_receipt_adjudication_cmp<I: StrongIdentity>(
    left: IdentityReceipt<I>,
    right: IdentityReceipt<I>,
) -> core::cmp::Ordering {
    left.id()
        .cmp(&right.id())
        .then_with(|| {
            left.canonical_preimage()
                .as_bytes()
                .cmp(right.canonical_preimage().as_bytes())
        })
        .then_with(|| {
            left.schema_id()
                .as_bytes()
                .cmp(right.schema_id().as_bytes())
        })
        .then_with(|| left.canonical_bytes().cmp(&right.canonical_bytes()))
        .then_with(|| left.field_count().cmp(&right.field_count()))
        .then_with(|| left.collection_items().cmp(&right.collection_items()))
}

fn identity_receipt_adjudication_eq<I: StrongIdentity>(
    left: IdentityReceipt<I>,
    right: IdentityReceipt<I>,
) -> bool {
    identity_receipt_adjudication_cmp(left, right) == core::cmp::Ordering::Equal
}

fn optional_identity_receipt_adjudication_eq<I: StrongIdentity>(
    left: Option<IdentityReceipt<I>>,
    right: Option<IdentityReceipt<I>>,
) -> bool {
    match (left, right) {
        (Some(left), Some(right)) => identity_receipt_adjudication_eq(left, right),
        (None, None) => true,
        (Some(_), None) | (None, Some(_)) => false,
    }
}

fn hash_identity_receipt_adjudication<I: StrongIdentity, H: Hasher>(
    receipt: IdentityReceipt<I>,
    state: &mut H,
) {
    receipt.id().hash(state);
    receipt.canonical_preimage().as_bytes().hash(state);
    receipt.schema_id().as_bytes().hash(state);
    receipt.canonical_bytes().hash(state);
    receipt.field_count().hash(state);
    receipt.collection_items().hash(state);
}

/// Minimal append-only byte surface shared by eager canonicalization and the
/// cap-enforcing streamed-row writer below.
///
/// Keeping this surface deliberately smaller than `Vec<u8>` prevents a
/// streamed producer from reserving, resizing, or otherwise bypassing its
/// admitted row plan. The `Vec<u8>` implementation preserves the existing
/// eager encoders; streamed producers receive only `CausalIdentityRowWriter`.
trait CausalCanonicalBytes {
    fn len(&self) -> usize;
    fn push(&mut self, byte: u8);
    fn extend_from_slice(&mut self, bytes: &[u8]);

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl CausalCanonicalBytes for Vec<u8> {
    fn len(&self) -> usize {
        Vec::len(self)
    }

    fn push(&mut self, byte: u8) {
        Vec::push(self, byte);
    }

    fn extend_from_slice(&mut self, bytes: &[u8]) {
        Vec::extend_from_slice(self, bytes);
    }
}

fn push_identity_receipt_adjudication<I: StrongIdentity, O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    receipt: IdentityReceipt<I>,
) {
    out.extend_from_slice(receipt.id().as_bytes());
    out.extend_from_slice(receipt.canonical_preimage().as_bytes());
    out.extend_from_slice(receipt.schema_id().as_bytes());
    out.extend_from_slice(&receipt.canonical_bytes().to_le_bytes());
    out.extend_from_slice(&receipt.field_count().to_le_bytes());
    out.extend_from_slice(&receipt.collection_items().to_le_bytes());
}

macro_rules! causal_ref {
    ($(#[$meta:meta])* $name:ident, $role:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            namespace: Box<str>,
            schema_version: NonZeroU64,
            semantic_digest: [u8; 32],
        }

        impl $name {
            /// Construct a versioned reference owned by another semantic domain.
            ///
            /// The digest is bound exactly; this constructor does not inspect,
            /// authenticate, execute, or strengthen the referenced artifact.
            ///
            /// # Errors
            /// Refuses a noncanonical namespace or an all-zero digest.
            pub fn new(
                namespace: impl Into<String>,
                schema_version: NonZeroU64,
                semantic_digest: [u8; 32],
            ) -> Result<Self, MachineReferenceError> {
                let namespace = namespace.into();
                super::validate_canonical_key($role, &namespace)
                    .map_err(MachineReferenceError::Namespace)?;
                if semantic_digest == [0; 32] {
                    return Err(MachineReferenceError::ZeroDigest { role: $role });
                }
                Ok(Self {
                    namespace: namespace.into_boxed_str(),
                    schema_version,
                    semantic_digest,
                })
            }

            /// External schema namespace.
            #[must_use]
            pub fn namespace(&self) -> &str {
                &self.namespace
            }

            /// Explicit external schema version.
            #[must_use]
            pub const fn schema_version(&self) -> NonZeroU64 {
                self.schema_version
            }

            /// Exact semantic digest supplied by the external owner.
            #[must_use]
            pub const fn semantic_digest(&self) -> [u8; 32] {
                self.semantic_digest
            }

            fn append_canonical<O: CausalCanonicalBytes + ?Sized>(&self, out: &mut O) {
                push_len_prefixed(out, self.namespace.as_bytes());
                out.extend_from_slice(&self.schema_version.get().to_le_bytes());
                out.extend_from_slice(&self.semantic_digest);
            }

        }
    };
}

causal_ref!(
    /// Exact source or hand-operator artifact named by an escape hatch.
    SourceArtifactRef,
    "causal-source-artifact-ref"
);
causal_ref!(
    /// Instance-local normalized equation/variable meaning shared by
    /// conformant generated and hand-authored producers.
    NormalizedNodeSemanticRef,
    "causal-normalized-node-semantic-ref"
);
causal_ref!(
    /// Separate audit receipt required before an opaque source may enter the graph.
    EscapeHatchAuditRef,
    "causal-escape-hatch-audit-ref"
);
causal_ref!(
    /// Semantic definition of a parameter, guard, or hybrid activation domain.
    ActivationConditionRef,
    "causal-activation-condition-ref"
);
causal_ref!(
    /// One named branch/mode inside an activation-condition domain.
    ActivationBranchRef,
    "causal-activation-branch-ref"
);
causal_ref!(
    /// Obligation governing simultaneous/root-solved guard evaluation.
    GuardSolveObligationRef,
    "causal-guard-solve-obligation-ref"
);
causal_ref!(
    /// Semantic definition of an incidence operator or transform.
    IncidenceOperatorRef,
    "causal-incidence-operator-ref"
);
causal_ref!(
    /// Semantic definition of a provenance-preserving derivation step.
    DerivationRuleRef,
    "causal-derivation-rule-ref"
);
causal_ref!(
    /// Checkable implication/equivalence obligation for a derived node.
    DerivationObligationRef,
    "causal-derivation-obligation-ref"
);
causal_ref!(
    /// Normalized transfer semantics between two logical clock domains.
    ClockBridgeRef,
    "causal-clock-bridge-ref"
);
causal_ref!(
    /// Audit receipt for one exact clock-bridge implementation.
    ClockBridgeAuditRef,
    "causal-clock-bridge-audit-ref"
);
causal_ref!(
    /// External geometric or topological support identity.
    ExternalSupportRef,
    "causal-external-support-ref"
);
causal_ref!(
    /// Normalized semantic crosswalk from one `fs-couple::PortSchema`
    /// projection into one dependency-neutral Machine-IR port.
    PortSchemaCrosswalkRef,
    "causal-port-schema-crosswalk-ref"
);
causal_ref!(
    /// Audit receipt for one exact PortSchema crosswalk implementation.
    PortSchemaCrosswalkAuditRef,
    "causal-port-schema-crosswalk-audit-ref"
);
causal_ref!(
    /// Analyzer implementation and configuration identity.
    CausalAnalyzerRef,
    "causal-analyzer-ref"
);
causal_ref!(
    /// Exact extractor implementation and configuration identity.
    CausalExtractorRef,
    "causal-extractor-ref"
);
causal_ref!(
    /// Exact assertion of source-domain coverage for one extraction.
    CausalExtractionCoverageRef,
    "causal-extraction-coverage-ref"
);
causal_ref!(
    /// Independent extraction-coverage checker artifact identity.
    CausalExtractionCheckerRef,
    "causal-extraction-checker-ref"
);
causal_ref!(
    /// Exact semantic boundary of an intentionally partial extraction.
    CausalPartialScopeRef,
    "causal-partial-scope-ref"
);
causal_ref!(
    /// Exact time, memory, and accuracy budget declaration.
    CausalBudgetRef,
    "causal-budget-ref"
);
causal_ref!(
    /// Exact admitted capability-set declaration.
    CausalCapabilityRef,
    "causal-capability-ref"
);
causal_ref!(
    /// Certificate that a retained matching is maximum for its declared graph.
    MaximumMatchingCertificateRef,
    "causal-maximum-matching-certificate-ref"
);
causal_ref!(
    /// Coverage/uniformity certificate for a declared hybrid-mode domain.
    ConditionalCoverageRef,
    "causal-conditional-coverage-ref"
);
causal_ref!(
    /// Deterministic analyzer checkpoint that may resume an incomplete result.
    CausalCheckpointRef,
    "causal-checkpoint-ref"
);
causal_ref!(
    /// Checker identity referenced by a receipt without granting authority.
    CausalCheckerRef,
    "causal-checker-ref"
);
causal_ref!(
    /// Explicit migration receipt for one predecessor artifact.
    CausalMigrationRef,
    "causal-migration-ref"
);
/// Effort/flow coordinate within one already admitted Machine-IR port.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PortCoordinate {
    /// Generalized effort coordinate.
    Effort,
    /// Generalized flow coordinate.
    Flow,
}

/// Explicit unit convention of every structural signal in a v1 graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalUnitConvention {
    /// Quantities are expressed as canonical SI-base dimensions; semantic
    /// quantity kinds remain distinct even when their dimensions coincide.
    SiBaseDimensions,
}

/// Determinism class declared by an extractor or analyzer invocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalDeterminism {
    /// Bit-stable output under collection permutations and worker-count changes.
    Deterministic,
    /// A deliberately relaxed mode, retained explicitly in artifact identity.
    Relaxed,
}

/// Explicit random-stream declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalSeedPolicy {
    /// The operation is structurally deterministic and consumes no RNG stream.
    NoRandomness,
    /// Counter-based stream key supplied by the caller. The key is semantic,
    /// not a worker or thread identifier.
    CounterBased { seed: u64, stream: u64 },
}

/// Semantic extent of one equation extraction.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalGraphScope {
    /// Extractor asserts coverage of its complete admitted Machine-model input.
    CompleteMachineModel,
    /// Extractor intentionally exposes only an exact named semantic boundary.
    Partial {
        /// Boundary definition whose contents are owned by the source domain.
        boundary: CausalPartialScopeRef,
    },
}

/// Honest evidence state for an extraction coverage assertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CausalExtractionEvidence {
    /// Coverage is a producer assertion only.
    Unverified,
    /// An exact independent checker artifact is retained without being
    /// authenticated or promoted by this schema boundary.
    CheckerReferenced(CausalExtractionCheckerRef),
}

/// Five-Explicits and provenance context of equation extraction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalExtractionContext {
    /// Exact extractor implementation/configuration and version.
    pub extractor: CausalExtractorRef,
    /// Exact source-domain coverage assertion; it is retained but not
    /// authenticated by this schema-only boundary.
    pub coverage: CausalExtractionCoverageRef,
    /// Explicit evidence maturity of the coverage assertion.
    pub evidence: CausalExtractionEvidence,
    /// Exact time/memory/accuracy budget.
    pub budget: CausalBudgetRef,
    /// Exact capability set under which extraction ran.
    pub capabilities: CausalCapabilityRef,
    /// Explicit seed or no-randomness declaration.
    pub seed_policy: CausalSeedPolicy,
    /// Declared determinism class.
    pub determinism: CausalDeterminism,
}

/// Five-Explicits and provenance context of causal analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalAnalysisContext {
    /// Exact analyzer implementation/configuration and version.
    pub analyzer: CausalAnalyzerRef,
    /// Exact time/memory/accuracy budget.
    pub budget: CausalBudgetRef,
    /// Exact capability set under which analysis ran.
    pub capabilities: CausalCapabilityRef,
    /// Explicit seed or no-randomness declaration.
    pub seed_policy: CausalSeedPolicy,
    /// Declared determinism class.
    pub determinism: CausalDeterminism,
}

/// Machine-owned source from which a node's semantic identity descends.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineNodeOrigin {
    /// Declarative subsystem/model boundary.
    Subsystem(SubsystemId),
    /// Directed Machine-IR relation.
    Relation(RelationId),
    /// Durable body, support, contact, terminal, port, or state element.
    Element(MachineElementId),
    /// Exact effort or flow coordinate of one port.
    PortCoordinate {
        /// Owning Machine-IR port.
        port: PortId,
        /// Coordinate selected within the port.
        coordinate: PortCoordinate,
    },
    /// Role-oriented Machine-IR interface.
    Interface(InterfaceId),
}

/// Audited escape hatch for a source not yet represented by generated Machine IR.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AuditedEscapeHatch {
    /// Exact opaque source identity.
    pub source: SourceArtifactRef,
    /// Separate audit/conformance receipt identity.
    pub audit: EscapeHatchAuditRef,
    /// Source identity the audit declares it checked. Admission requires exact
    /// equality with `source`, making the declared source coordinate explicit
    /// in the published artifact tuple. This structural equality check does
    /// not authenticate the opaque audit artifact or inspect its contents.
    pub audited_source: SourceArtifactRef,
}

/// Nominal parent identity in a provenance-preserving derivation.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ParentNodeId {
    /// Parent equation.
    Equation(EquationId),
    /// Parent base variable.
    Variable(VariableId),
}

/// Transitive lineage of a differentiated, eliminated, or transformed node.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedNodeLineage {
    /// Canonical nonempty parent set.
    pub parents: Vec<ParentNodeId>,
    /// Exact transformation semantics.
    pub transformation: DerivationRuleRef,
    /// Time-differentiation order applied by this derivation step.
    pub differentiation_order: u16,
    /// Separately checkable semantic-equivalence or implication obligation.
    pub obligation: DerivationObligationRef,
}

/// Closed node-origin vocabulary.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NodeOrigin {
    /// Traceable to an entity in the exact admitted Machine graph.
    Machine(MachineNodeOrigin),
    /// Opaque source with a separately bound audit receipt.
    AuditedEscapeHatch(AuditedEscapeHatch),
    /// Provenance-preserving child of existing equation/variable identities.
    Derived(DerivedNodeLineage),
}

/// Identity-bearing lineage of one equation or variable.
///
/// `semantic` is produced by the source-domain extractor and must be stable
/// under presentation-only renaming, collection reordering, and conformant
/// generated-versus-hand implementations. `instance` prevents two occurrences
/// of the same reusable component template from aliasing. A display label is
/// intentionally absent.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NodeLineage {
    /// Machine-owned or audited-escape source.
    pub origin: NodeOrigin,
    /// Accountable Machine instance/port/interface for this occurrence.
    pub instance: CausalOwner,
    /// Versioned normalized meaning, never an array index or display hash.
    pub semantic: NormalizedNodeSemanticRef,
}

impl NodeLineage {
    /// Construct node lineage.
    #[must_use]
    pub fn new(
        origin: NodeOrigin,
        instance: CausalOwner,
        semantic: NormalizedNodeSemanticRef,
    ) -> Self {
        Self {
            origin,
            instance,
            semantic,
        }
    }

    fn append_canonical_row_cancellable<O: CausalCanonicalBytes + ?Sized>(
        &self,
        out: &mut O,
        cx: &Cx<'_>,
    ) -> Result<(), CanonicalError> {
        push_node_origin_cancellable(out, &self.origin, cx)?;
        push_owner(out, &self.instance);
        self.semantic.append_canonical(out);
        identity_materialization_checkpoint(cx, out.len())
    }

    fn canonical_row_len_cancellable(&self, cx: &Cx<'_>) -> Result<usize, CanonicalError> {
        checked_row_len_sum([
            node_origin_canonical_len_cancellable(&self.origin, cx)?,
            CAUSAL_OWNER_CANONICAL_BYTES,
            causal_ref_canonical_len(self.semantic.namespace())?,
        ])
    }

    fn normalized_row(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(160);
        push_owner(&mut out, &self.instance);
        self.semantic.append_canonical(&mut out);
        out
    }
}

/// Canonical schema marker for equation identities.
pub enum EquationIdentitySchemaV1 {}

impl CanonicalSchema for EquationIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.causal-equation.v1";
    const NAME: &'static str = "causal-equation-id";
    const VERSION: u32 = CAUSAL_GRAPH_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one normalized equation meaning independent of producer provenance, presentation, and collection order";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "instance-qualified-normalized-meaning",
        WireType::Bytes,
    )];
}

/// Typed digest of one equation identity.
pub type EquationEntityIdV1 = EntityId<EquationIdentitySchemaV1>;

/// Nominal identity of one structural equation.
#[derive(Clone)]
pub struct EquationId {
    receipt: Arc<IdentityReceipt<EquationEntityIdV1>>,
}

impl EquationId {
    /// Derive a normalized equation identity from the lineage's source-local
    /// semantic digest. Producer provenance is retained separately by the graph
    /// artifact identity.
    ///
    /// # Errors
    /// Returns a bounded canonical-identity refusal.
    pub fn derive(lineage: &NodeLineage) -> Result<Self, CanonicalError> {
        let normalized = lineage.normalized_row();
        let receipt =
            CanonicalEncoder::<EquationEntityIdV1, _>::new(NODE_IDENTITY_LIMITS, NeverCancel)?
                .bytes(
                    Field::new(0, "instance-qualified-normalized-meaning"),
                    &normalized,
                )?
                .finish()?;
        Ok(Self {
            receipt: Arc::new(receipt),
        })
    }

    /// Typed identity digest.
    #[must_use]
    pub fn identity(&self) -> EquationEntityIdV1 {
        self.receipt.id()
    }

    /// Complete standalone identity receipt for collision adjudication.
    #[must_use]
    pub fn identity_receipt(&self) -> IdentityReceipt<EquationEntityIdV1> {
        *self.receipt
    }
}

/// Canonical schema marker for variable identities.
pub enum VariableIdentitySchemaV1 {}

impl CanonicalSchema for VariableIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.causal-variable.v1";
    const NAME: &'static str = "causal-variable-id";
    const VERSION: u32 = CAUSAL_GRAPH_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one normalized base-variable meaning independent of producer provenance, presentation, and collection order";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "instance-qualified-normalized-meaning",
        WireType::Bytes,
    )];
}

/// Typed digest of one variable identity.
pub type VariableEntityIdV1 = EntityId<VariableIdentitySchemaV1>;

/// Nominal identity of one structural variable.
#[derive(Clone)]
pub struct VariableId {
    receipt: Arc<IdentityReceipt<VariableEntityIdV1>>,
}

impl VariableId {
    /// Derive a normalized variable identity from the lineage's source-local
    /// semantic digest. Derivative order belongs to incidence keys.
    ///
    /// # Errors
    /// Returns a bounded canonical-identity refusal.
    pub fn derive(lineage: &NodeLineage) -> Result<Self, CanonicalError> {
        let normalized = lineage.normalized_row();
        let receipt =
            CanonicalEncoder::<VariableEntityIdV1, _>::new(NODE_IDENTITY_LIMITS, NeverCancel)?
                .bytes(
                    Field::new(0, "instance-qualified-normalized-meaning"),
                    &normalized,
                )?
                .finish()?;
        Ok(Self {
            receipt: Arc::new(receipt),
        })
    }

    /// Typed identity digest.
    #[must_use]
    pub fn identity(&self) -> VariableEntityIdV1 {
        self.receipt.id()
    }

    /// Complete standalone identity receipt for collision adjudication.
    #[must_use]
    pub fn identity_receipt(&self) -> IdentityReceipt<VariableEntityIdV1> {
        *self.receipt
    }
}

/// Canonical schema marker for normalized structural incidences.
pub enum IncidenceIdentitySchemaV1 {}

impl CanonicalSchema for IncidenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.causal-incidence.v1";
    const NAME: &'static str = "causal-incidence-id";
    const VERSION: u32 = CAUSAL_GRAPH_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact normalized equation-to-derivative-variable structural edge, including activation and transfer semantics";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "normalized-incidence-meaning",
        WireType::Bytes,
    )];
}

/// Typed digest of one normalized structural incidence.
pub type IncidenceEntityIdV1 = EntityId<IncidenceIdentitySchemaV1>;

/// Nominal identity of one exact structural incidence.
#[derive(Clone)]
pub struct IncidenceId {
    receipt: Arc<IdentityReceipt<IncidenceEntityIdV1>>,
}

impl IncidenceId {
    /// Typed incidence identity digest.
    #[must_use]
    pub fn identity(&self) -> IncidenceEntityIdV1 {
        self.receipt.id()
    }

    /// Complete standalone identity receipt for collision adjudication.
    #[must_use]
    pub fn identity_receipt(&self) -> IdentityReceipt<IncidenceEntityIdV1> {
        *self.receipt
    }
}

macro_rules! id_semantics {
    ($name:ident) => {
        impl PartialEq for $name {
            fn eq(&self, other: &Self) -> bool {
                identity_receipt_adjudication_eq(self.identity_receipt(), other.identity_receipt())
            }
        }

        impl Eq for $name {}

        impl PartialOrd for $name {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Ord for $name {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                identity_receipt_adjudication_cmp(self.identity_receipt(), other.identity_receipt())
            }
        }

        impl Hash for $name {
            fn hash<H: Hasher>(&self, state: &mut H) {
                hash_identity_receipt_adjudication(self.identity_receipt(), state);
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let receipt = self.identity_receipt();
                f.debug_struct(stringify!($name))
                    .field("identity", &receipt.id())
                    .field("canonical_preimage", &receipt.canonical_preimage())
                    .field("schema_id", &receipt.schema_id())
                    .field("canonical_bytes", &receipt.canonical_bytes())
                    .field("field_count", &receipt.field_count())
                    .field("collection_items", &receipt.collection_items())
                    .finish()
            }
        }
    };
}

id_semantics!(EquationId);
id_semantics!(VariableId);
id_semantics!(IncidenceId);

/// Explicit computational owner of an equation or variable.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalOwner {
    /// Machine subsystem/component.
    Subsystem(SubsystemId),
    /// Exact effort/flow port.
    Port(PortId),
    /// Cross-component interface.
    Interface(InterfaceId),
}

/// Spatial/topological support of an equation or variable.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalSupport {
    /// Lumped quantity with no stronger support claim.
    Lumped,
    /// Durable element in the exact admitted Machine graph.
    MachineElement(MachineElementId),
    /// Versioned support owned by a geometry/topology domain.
    External(ExternalSupportRef),
}

/// Exact quantity, shape, clock, frame, and orientation contract.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SignalContract {
    /// Six-base dimensions and optional stronger quantity semantics.
    pub quantity: TerminalQuantitySpec,
    /// Scalar, vector, tensor, or field-trace shape.
    pub shape: TerminalShape,
    /// Logical clock domain.
    pub clock: ClockId,
    /// Frame and orientation binding.
    pub frame: FrameBinding,
}

/// Exact normalized PortSchema projection plus a source-bound audit receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortSchemaCrosswalkBinding {
    /// Normalized projection semantics used by structural identity.
    pub projection: PortSchemaCrosswalkRef,
    /// Audit artifact identity retained only in provenance identity.
    pub audit: PortSchemaCrosswalkAuditRef,
    /// Projection identity the audit declares it checked. Admission requires
    /// exact equality with `projection`.
    pub audited_projection: PortSchemaCrosswalkRef,
}

/// Structural role of an equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EquationRole {
    /// Conservation or balance relation.
    Balance,
    /// Constitutive relation.
    Constitutive,
    /// Algebraic or geometric constraint.
    Constraint,
    /// Prescribed source relation.
    Source,
    /// State-update relation.
    StateUpdate,
    /// Effort/flow or interface closure relation.
    PortClosure,
    /// Guard or mode-selection relation.
    Guard,
}

/// Whether an equation belongs to the left-hand matching vertex set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EquationParticipation {
    /// Residual equation that must match one structural unknown occurrence.
    Matching,
    /// Prescribed source/closure row retained structurally but not matched.
    KnownClosure,
    /// Guard/mode relation used only to define activation.
    ConditionOnly,
}

/// Structural role of a variable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VariableRole {
    /// Algebraic unknown.
    Algebraic,
    /// Owned state value.
    State,
    /// Prescribed source value.
    Source,
    /// Parameter that is not solved by this graph.
    Parameter,
    /// Effort coordinate of a Machine port.
    PortEffort,
    /// Flow coordinate of a Machine port.
    PortFlow,
    /// Discrete mode or guard value.
    DiscreteMode,
}

/// Whether a variable belongs to the right-hand matching vertex set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SolveParticipation {
    /// Solved structural unknown. Algebraic/port base values contribute order
    /// zero explicitly; a retained [`VariableRole::State`] is represented by
    /// its active derivative occurrence when one exists.
    Unknown,
    /// Known source/boundary value that is incident but never matched.
    KnownInput,
    /// Condition dependency excluded from matching.
    ConditionOnly,
    /// Participation is resolved by mutually exclusive activation cells on
    /// incidences; the base declaration itself is never matched.
    ModeDependent,
}

/// Activation domain of a node or incidence.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ConditionBranchSelection {
    /// Condition definition.
    pub condition: ActivationConditionRef,
    /// Selected branch in that condition's finite domain.
    pub branch: ActivationBranchRef,
}

/// One conjunction inside a canonical disjunctive-normal-form activation.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ActivationCube {
    /// Nonempty, duplicate-condition-free conjunction.
    pub selections: Vec<ConditionBranchSelection>,
}

/// Activation domain of a node or incidence.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ActivationDomain {
    /// Always structurally present.
    Always,
    /// Present in a canonical disjunction of conjunctions.
    Conditional {
        /// Nonempty, bounded, duplicate-free DNF cubes. Admission sorts cubes
        /// and their selections; exact implication is decided separately by a
        /// budgeted finite-domain symbolic search.
        cubes: Vec<ActivationCube>,
    },
}

/// One equation declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EquationSpec {
    /// Identity derived from `lineage`.
    pub id: EquationId,
    /// Presentation-only diagnostic text, excluded from graph identity.
    pub diagnostic_label: Box<str>,
    /// Complete source lineage.
    pub lineage: NodeLineage,
    /// Computational owner.
    pub owner: CausalOwner,
    /// Spatial/topological support.
    pub supports: Vec<CausalSupport>,
    /// Residual quantity contract shared by every admitted term.
    pub residual: SignalContract,
    /// Source, constraint, state-update, or other structural role.
    pub role: EquationRole,
    /// Explicit participation in matching and DM decomposition.
    pub solve_participation: EquationParticipation,
    /// Always-active or condition-bound structure.
    pub activation: ActivationDomain,
}

/// One variable declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariableSpec {
    /// Identity derived from `lineage`.
    pub id: VariableId,
    /// Presentation-only diagnostic text, excluded from graph identity.
    pub diagnostic_label: Box<str>,
    /// Complete source lineage.
    pub lineage: NodeLineage,
    /// Computational owner.
    pub owner: CausalOwner,
    /// Spatial/topological support.
    pub supports: Vec<CausalSupport>,
    /// Value quantity contract before incidence-local differentiation.
    pub value: SignalContract,
    /// State, algebraic, source, parameter, port, or mode role.
    pub role: VariableRole,
    /// Explicit participation in matching and DM decomposition.
    pub solve_participation: SolveParticipation,
    /// Exact audited PortSchema-to-Machine projection for port coordinates.
    /// Required only for `PortEffort` and `PortFlow`; logical timestamp ticks
    /// are execution state and must not enter this structural identity.
    pub port_schema_crosswalk: Option<PortSchemaCrosswalkBinding>,
    /// Always-active or condition-bound structure.
    pub activation: ActivationDomain,
}

/// Checkable dependency table for one parameter/guard activation condition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ActivationConditionSource {
    /// One explicit Guard equation and its simultaneous/root-solve obligation.
    GuardEquation {
        /// Exact guard equation.
        equation: EquationId,
        /// Obligation preventing circular mode selection from being treated as
        /// an ordinary upstream Boolean evaluation.
        obligation: GuardSolveObligationRef,
    },
    /// Audited external predicate whose audit is bound to its exact source.
    AuditedPredicate(AuditedEscapeHatch),
}

/// Checkable dependency table for one parameter/guard activation condition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ActivationConditionSpec {
    /// Exact condition semantics used by activation rows.
    pub condition: ActivationConditionRef,
    /// Exact predicate/guard source.
    pub source: ActivationConditionSource,
    /// Canonical nonempty branch/mode set.
    pub branches: Vec<ActivationBranchRef>,
    /// Canonical nonempty set of variables read by the condition.
    pub dependencies: Vec<VariableId>,
}

/// Logical-clock relationship applied by one incidence.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IncidenceClockRelation {
    /// Variable and residual inhabit the exact same clock domain.
    SameClock,
    /// A separately audited bridge maps distinct clock domains.
    AuditedBridge {
        /// Exact source clock; must equal the variable contract clock.
        source: ClockId,
        /// Exact destination clock; must equal the equation term clock.
        target: ClockId,
        /// Transfer/resampling semantics and audit receipt.
        bridge: ClockBridgeRef,
        /// Separate audit receipt for the bridge implementation.
        audit: ClockBridgeAuditRef,
        /// Transfer identity the audit declares it checked. Admission requires
        /// exact equality with `bridge`.
        audited_bridge: ClockBridgeRef,
    },
}

/// One structural incidence between an equation and a variable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IncidenceSpec {
    /// Identity derived from the complete normalized incidence meaning.
    pub id: IncidenceId,
    /// Incident equation.
    pub equation: EquationId,
    /// Incident variable.
    pub variable: VariableId,
    /// Number of time derivatives applied to the variable at this occurrence.
    pub derivative_order: u16,
    /// Matching participation of this exact derivative-variable occurrence.
    /// This is independent of the retained base/state value at order zero.
    pub solve_participation: SolveParticipation,
    /// Dimension delta contributed by the coefficient/operator, before the
    /// time-derivative delta is applied.
    pub coefficient_dimensions: Dims,
    /// Exact term contract after derivative/operator application. Admission
    /// requires this to equal the equation residual contract.
    pub term: SignalContract,
    /// Opaque operator semantics when shape/frame/semantic-kind conversion is
    /// stronger than a direct structural occurrence.
    pub operator: Option<IncidenceOperatorRef>,
    /// Same-clock declaration or audited cross-clock bridge.
    pub clock_relation: IncidenceClockRelation,
    /// Always-active or condition-bound structural occurrence.
    pub activation: ActivationDomain,
}

impl IncidenceSpec {
    /// Construct an incidence and derive its exact normalized identity.
    ///
    /// # Errors
    /// Returns a bounded canonical-identity refusal, including cancellation
    /// before any incidence identity is published.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        equation: EquationId,
        variable: VariableId,
        derivative_order: u16,
        solve_participation: SolveParticipation,
        coefficient_dimensions: Dims,
        term: SignalContract,
        operator: Option<IncidenceOperatorRef>,
        clock_relation: IncidenceClockRelation,
        mut activation: ActivationDomain,
        cx: &Cx<'_>,
    ) -> Result<Self, CanonicalError> {
        preflight_constructor_activation(&activation, cx)?;
        canonicalize_activation_for_identity(&mut activation, cx)?;
        let meaning = incidence_meaning_row_parts_cancellable(
            &equation,
            &variable,
            derivative_order,
            solve_participation,
            coefficient_dimensions,
            &term,
            operator.as_ref(),
            &clock_relation,
            &activation,
            cx,
        )?;
        let receipt = incidence_receipt_from_meaning_cancellable(&meaning, cx)?;
        let id = IncidenceId {
            receipt: Arc::new(receipt),
        };
        Ok(Self {
            id,
            equation,
            variable,
            derivative_order,
            solve_participation,
            coefficient_dimensions,
            term,
            operator,
            clock_relation,
            activation,
        })
    }
}

fn preflight_constructor_activation(
    activation: &ActivationDomain,
    cx: &Cx<'_>,
) -> Result<(), CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let cubes = activation_cubes(activation);
    if cubes.len() > MAX_CAUSAL_CUBES_PER_ACTIVATION {
        return Err(collection_limit_error(
            cubes.len(),
            MAX_CAUSAL_CUBES_PER_ACTIVATION,
        ));
    }
    let mut selections = 0usize;
    for (cube_index, cube) in cubes.iter().enumerate() {
        identity_materialization_poll(cx, cube_index, 0)?;
        if cube.selections.len() > MAX_CAUSAL_SELECTIONS_PER_CUBE {
            return Err(collection_limit_error(
                cube.selections.len(),
                MAX_CAUSAL_SELECTIONS_PER_CUBE,
            ));
        }
        selections = selections.saturating_add(cube.selections.len());
        if selections > MAX_CAUSAL_ACTIVATION_SELECTIONS {
            return Err(collection_limit_error(
                selections,
                MAX_CAUSAL_ACTIVATION_SELECTIONS,
            ));
        }
    }
    identity_materialization_checkpoint(cx, 0)?;
    Ok(())
}

fn collection_limit_error(requested: usize, limit: usize) -> CanonicalError {
    CanonicalError::LimitExceeded {
        kind: LimitKind::CollectionItems,
        requested: u64::try_from(requested).unwrap_or(u64::MAX),
        limit: u64::try_from(limit).unwrap_or(u64::MAX),
    }
}

/// Mutable, authority-free equation-variable graph draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalGraphDraft {
    /// Explicit SI-base dimensional convention.
    pub units: CausalUnitConvention,
    /// Complete-model or explicitly bounded partial extraction.
    pub scope: CausalGraphScope,
    /// Producer, coverage, budget, seed, capability, and determinism context.
    pub extraction: CausalExtractionContext,
    /// Equations in arbitrary caller order.
    pub equations: Vec<EquationSpec>,
    /// Variables in arbitrary caller order.
    pub variables: Vec<VariableSpec>,
    /// Activation-condition dependency tables in arbitrary caller order.
    pub conditions: Vec<ActivationConditionSpec>,
    /// Structural incidences in arbitrary caller order.
    pub incidences: Vec<IncidenceSpec>,
}

const MACHINE_GRAPH_CHILD: ChildSpec = ChildSpec::for_identity::<MachineGraphIdV1>();
const INCIDENCE_CHILD: ChildSpec = ChildSpec::for_identity::<IncidenceEntityIdV1>();

/// Canonical schema marker for normalized equation-variable structure.
pub enum CausalStructureIdentitySchemaV1 {}

impl CanonicalSchema for CausalStructureIdentitySchemaV1 {
    const DOMAIN: &'static str = CAUSAL_STRUCTURE_IDENTITY_DOMAIN_V1;
    const NAME: &'static str = "normalized-causal-equation-variable-hypergraph";
    const VERSION: u32 = CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "normalized structural incidence bound to one admitted Machine graph; producer provenance, numerical rank, and physical causality excluded";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("causal-structure-schema-version", WireType::U64),
        FieldSpec::child_of("machine-graph-id", &MACHINE_GRAPH_CHILD),
        FieldSpec::required("machine-graph-receipt-adjudication", WireType::Bytes),
        FieldSpec::required("unit-convention", WireType::Variant),
        FieldSpec::required("extraction-scope", WireType::Bytes),
        FieldSpec::required("equations", WireType::OrderedBytes),
        FieldSpec::required("variables", WireType::OrderedBytes),
        FieldSpec::required("activation-conditions", WireType::OrderedBytes),
        FieldSpec::ordered_children_of("incidences", &INCIDENCE_CHILD),
        FieldSpec::required("incidence-receipt-adjudications", WireType::OrderedBytes),
    ];
}

/// Producer-independent structural identity used for conformance comparison.
pub type CausalStructureIdV1 = ProblemSemanticId<CausalStructureIdentitySchemaV1>;

const CAUSAL_STRUCTURE_CHILD: ChildSpec = ChildSpec::for_identity::<CausalStructureIdV1>();

/// Canonical schema marker for the full provenance-bearing graph artifact.
pub enum CausalGraphArtifactIdentitySchemaV1 {}

impl CanonicalSchema for CausalGraphArtifactIdentitySchemaV1 {
    const DOMAIN: &'static str = CAUSAL_GRAPH_ARTIFACT_IDENTITY_DOMAIN_V1;
    const NAME: &'static str = "provenance-bearing-causal-graph-artifact";
    const VERSION: u32 = CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "one normalized causal structure plus exact generated, derived, or audited-escape lineage";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("causal-graph-artifact-schema-version", WireType::U64),
        FieldSpec::child_of("causal-structure-id", &CAUSAL_STRUCTURE_CHILD),
        FieldSpec::required("causal-structure-receipt-adjudication", WireType::Bytes),
        FieldSpec::optional_bytes("machine-behavior-receipt-adjudication"),
        FieldSpec::required("extraction-context", WireType::Bytes),
        FieldSpec::required("equation-lineage", WireType::OrderedBytes),
        FieldSpec::required("variable-lineage", WireType::OrderedBytes),
        FieldSpec::required("activation-condition-provenance", WireType::OrderedBytes),
        FieldSpec::required("incidence-provenance", WireType::OrderedBytes),
    ];
}

/// Full provenance-bearing identity of one admitted graph artifact.
pub type CausalGraphArtifactIdV1 = EvidenceNodeId<CausalGraphArtifactIdentitySchemaV1>;

const CAUSAL_GRAPH_ARTIFACT_CHILD: ChildSpec = ChildSpec::for_identity::<CausalGraphArtifactIdV1>();

/// Stable graph-admission rule vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CausalGraphRule {
    /// A public collection, nested label, or aggregate count exceeded its cap.
    ResourceLimit = 1,
    /// Two equations used the same nominal identity.
    DuplicateEquation = 2,
    /// Two variables used the same nominal identity.
    DuplicateVariable = 3,
    /// Two incidences used the same nominal identity.
    DuplicateIncidence = 4,
    /// Two activation tables used the same condition identity.
    DuplicateCondition = 5,
    /// Equation identity did not derive from its attached lineage.
    EquationIdentityMismatch = 6,
    /// Variable identity did not derive from its attached lineage.
    VariableIdentityMismatch = 7,
    /// A presentation label was empty, too large, or contained control text.
    InvalidDiagnosticLabel = 8,
    /// Node owner did not exist in the admitted Machine graph.
    UnknownOwner = 9,
    /// Node origin did not exist in the admitted Machine graph.
    UnknownOrigin = 10,
    /// Declared owner was inconsistent with the source origin.
    OwnerOriginMismatch = 11,
    /// Machine-element support did not exist in the admitted Machine graph.
    UnknownSupport = 12,
    /// Support set was empty, duplicated, or mixed lumped/strong supports.
    InvalidSupportSet = 13,
    /// Signal quantity used an unsupported semantic value form.
    UnsupportedQuantityForm = 14,
    /// Signal contract named an unknown Machine clock.
    UnknownClock = 15,
    /// State or port role was not backed by the corresponding Machine origin.
    RoleOriginMismatch = 16,
    /// Port role omitted or spuriously carried its audited PortSchema crosswalk.
    PortCrosswalkMismatch = 17,
    /// Role contradicted explicit known/unknown/condition-only participation.
    SolveParticipationMismatch = 18,
    /// A state-bearing graph omitted or contradicted the Machine behavior overlay.
    StateBehaviorMismatch = 19,
    /// Derived lineage was dangling, cyclic, duplicated, or out of bounds.
    InvalidDerivedLineage = 20,
    /// Activation table was missing, empty, dangling, or malformed.
    InvalidActivationCondition = 21,
    /// Incidence named an unknown equation.
    UnknownIncidenceEquation = 22,
    /// Incidence named an unknown variable.
    UnknownIncidenceVariable = 23,
    /// Derivative order exceeded the current candidate-v1 public bound.
    DerivativeOrderLimit = 24,
    /// Incidence term did not equal the equation's residual contract.
    ResidualTermMismatch = 25,
    /// Variable, derivative, coefficient, and term dimensions did not close.
    IncidenceUnitMismatch = 26,
    /// Same-clock or bridged-clock declaration was inconsistent.
    IncidenceClockMismatch = 27,
    /// A semantic/shape/frame transform omitted its operator identity.
    MissingOperatorSemantics = 28,
    /// Incidence activation was not contained in both endpoint domains.
    ActivationMismatch = 29,
    /// Bounded canonical identity publication refused.
    Identity = 30,
    /// Admission observed cancellation and published no identity.
    Cancelled = 31,
    /// Incidence identity did not derive from its normalized meaning.
    IncidenceIdentityMismatch = 32,
    /// More than one variable claimed one exclusive state/port coordinate.
    DuplicateRoleBinding = 34,
    /// Escape-hatch audit declared a different source artifact.
    EscapeAuditMismatch = 35,
    /// Derivative-key solve participation was inconsistent across incidences
    /// or with the base order-zero variable.
    DerivativeParticipationMismatch = 36,
}

impl CausalGraphRule {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ResourceLimit => "CausalGraphResourceLimit",
            Self::DuplicateEquation => "CausalGraphDuplicateEquation",
            Self::DuplicateVariable => "CausalGraphDuplicateVariable",
            Self::DuplicateIncidence => "CausalGraphDuplicateIncidence",
            Self::DuplicateCondition => "CausalGraphDuplicateCondition",
            Self::EquationIdentityMismatch => "CausalGraphEquationIdentityMismatch",
            Self::VariableIdentityMismatch => "CausalGraphVariableIdentityMismatch",
            Self::InvalidDiagnosticLabel => "CausalGraphInvalidDiagnosticLabel",
            Self::UnknownOwner => "CausalGraphUnknownOwner",
            Self::UnknownOrigin => "CausalGraphUnknownOrigin",
            Self::OwnerOriginMismatch => "CausalGraphOwnerOriginMismatch",
            Self::UnknownSupport => "CausalGraphUnknownSupport",
            Self::InvalidSupportSet => "CausalGraphInvalidSupportSet",
            Self::UnsupportedQuantityForm => "CausalGraphUnsupportedQuantityForm",
            Self::UnknownClock => "CausalGraphUnknownClock",
            Self::RoleOriginMismatch => "CausalGraphRoleOriginMismatch",
            Self::PortCrosswalkMismatch => "CausalGraphPortCrosswalkMismatch",
            Self::SolveParticipationMismatch => "CausalGraphSolveParticipationMismatch",
            Self::StateBehaviorMismatch => "CausalGraphStateBehaviorMismatch",
            Self::InvalidDerivedLineage => "CausalGraphInvalidDerivedLineage",
            Self::InvalidActivationCondition => "CausalGraphInvalidActivationCondition",
            Self::UnknownIncidenceEquation => "CausalGraphUnknownIncidenceEquation",
            Self::UnknownIncidenceVariable => "CausalGraphUnknownIncidenceVariable",
            Self::DerivativeOrderLimit => "CausalGraphDerivativeOrderLimit",
            Self::ResidualTermMismatch => "CausalGraphResidualTermMismatch",
            Self::IncidenceUnitMismatch => "CausalGraphIncidenceUnitMismatch",
            Self::IncidenceClockMismatch => "CausalGraphIncidenceClockMismatch",
            Self::MissingOperatorSemantics => "CausalGraphMissingOperatorSemantics",
            Self::ActivationMismatch => "CausalGraphActivationMismatch",
            Self::Identity => "CausalGraphIdentity",
            Self::Cancelled => "CausalGraphCancelled",
            Self::IncidenceIdentityMismatch => "CausalGraphIncidenceIdentityMismatch",
            Self::DuplicateRoleBinding => "CausalGraphDuplicateRoleBinding",
            Self::EscapeAuditMismatch => "CausalGraphEscapeAuditMismatch",
            Self::DerivativeParticipationMismatch => "CausalGraphDerivativeParticipationMismatch",
        }
    }
}

/// Stable subject of one graph-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalGraphSubject {
    /// Complete graph.
    Graph,
    /// One equation identity.
    Equation(EquationId),
    /// One variable identity.
    Variable(VariableId),
    /// One structural incidence key.
    Incidence {
        /// Exact normalized incidence identity.
        incidence: IncidenceId,
        /// Incident equation.
        equation: EquationId,
        /// Incident variable.
        variable: VariableId,
        /// Differentiation order.
        derivative_order: u16,
    },
    /// One duplicated nominal incidence identity whose conflicting endpoints
    /// are deliberately not selected from caller order.
    IncidenceIdentity(IncidenceId),
    /// One declared owner.
    Owner(CausalOwner),
    /// One activation condition.
    Condition(ActivationConditionRef),
}

/// One deterministic, localized graph-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CausalGraphFinding {
    rule: CausalGraphRule,
    subject: CausalGraphSubject,
}

impl CausalGraphFinding {
    fn new(rule: CausalGraphRule, subject: CausalGraphSubject) -> Self {
        Self { rule, subject }
    }

    /// Stable rule.
    #[must_use]
    pub const fn rule(&self) -> CausalGraphRule {
        self.rule
    }

    /// Localized subject.
    #[must_use]
    pub const fn subject(&self) -> &CausalGraphSubject {
        &self.subject
    }

    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }
}

/// Complete deterministic refusal from graph admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalGraphRefusal {
    findings: Vec<CausalGraphFinding>,
    identity_error: Option<CanonicalError>,
    identity_allocation_refusal: Option<CausalIdentityAllocationRefusal>,
    identity_stream_refusal: Option<CausalIdentityStreamRefusal>,
}

impl CausalGraphRefusal {
    /// Sorted, duplicate-free findings.
    #[must_use]
    pub fn findings(&self) -> &[CausalGraphFinding] {
        &self.findings
    }

    /// Canonical identity error, absent for the distinct fallible-allocation
    /// refusal exposed by [`Self::identity_allocation_refusal`].
    #[must_use]
    pub const fn identity_error(&self) -> Option<&CanonicalError> {
        self.identity_error.as_ref()
    }

    /// Exact fallible row-storage refusal, when identity publication could not
    /// allocate the already-admitted logical row.
    #[must_use]
    pub const fn identity_allocation_refusal(&self) -> Option<CausalIdentityAllocationRefusal> {
        self.identity_allocation_refusal
    }

    /// Exact ordered-row progress when identity publication refused inside a
    /// streamed collection.
    #[must_use]
    pub const fn identity_stream_refusal(&self) -> Option<CausalIdentityStreamRefusal> {
        self.identity_stream_refusal
    }
}

impl fmt::Display for CausalGraphRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "causal graph refused with {} deterministic finding(s)",
            self.findings.len()
        )
    }
}

impl std::error::Error for CausalGraphRefusal {}

/// Submitted counts retained for resource diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CausalGraphSubmittedCounts {
    /// Whether the complete telemetry pass reached its end without
    /// cancellation or an early resource-cap refusal. Aggregate fields are
    /// prefix observations and must not be treated as exact when this is
    /// `false`; top-level collection lengths remain exact.
    pub complete: bool,
    /// Submitted equations.
    pub equations: usize,
    /// Submitted variables.
    pub variables: usize,
    /// Submitted activation-condition tables.
    pub conditions: usize,
    /// Submitted incidences.
    pub incidences: usize,
    /// Aggregate condition dependency references.
    pub condition_dependencies: usize,
    /// Aggregate declared branches/modes.
    pub condition_branches: usize,
    /// Aggregate node support references.
    pub supports: usize,
    /// Aggregate selections across node/incidence activation conjunctions.
    pub activation_selections: usize,
    /// Aggregate DNF cubes across node/incidence activations.
    pub activation_cubes: usize,
    /// Aggregate parent references across derived node lineages.
    pub derivation_references: usize,
    /// Aggregate diagnostic-label bytes.
    pub diagnostic_label_bytes: usize,
    /// Largest direct-parent list on one derived node.
    pub max_derivation_parents: usize,
    /// Largest DNF cube count on one node or incidence activation.
    pub max_activation_cubes_per_row: usize,
    /// Largest selection count on one activation cube.
    pub max_activation_selections_per_cube: usize,
}

/// Canonically ordered, admitted structural graph.
#[derive(Debug, Clone)]
pub struct AdmittedCausalGraph {
    machine_graph: IdentityReceipt<MachineGraphIdV1>,
    machine_behavior: Option<IdentityReceipt<MachineBehaviorIdV1>>,
    units: CausalUnitConvention,
    scope: CausalGraphScope,
    extraction: CausalExtractionContext,
    equations: Vec<EquationSpec>,
    variables: Vec<VariableSpec>,
    conditions: Vec<ActivationConditionSpec>,
    incidences: Vec<IncidenceSpec>,
    structure_receipt: IdentityReceipt<CausalStructureIdV1>,
    artifact_receipt: IdentityReceipt<CausalGraphArtifactIdV1>,
}

impl PartialEq for AdmittedCausalGraph {
    fn eq(&self, other: &Self) -> bool {
        identity_receipt_adjudication_eq(self.machine_graph, other.machine_graph)
            && optional_identity_receipt_adjudication_eq(
                self.machine_behavior,
                other.machine_behavior,
            )
            && self.units == other.units
            && self.scope == other.scope
            && self.extraction == other.extraction
            && self.equations == other.equations
            && self.variables == other.variables
            && self.conditions == other.conditions
            && self.incidences == other.incidences
            && identity_receipt_adjudication_eq(self.structure_receipt, other.structure_receipt)
            && identity_receipt_adjudication_eq(self.artifact_receipt, other.artifact_receipt)
    }
}

impl Eq for AdmittedCausalGraph {}

impl AdmittedCausalGraph {
    /// Producer-independent identity of normalized structure.
    #[must_use]
    pub const fn structure_identity(&self) -> CausalStructureIdV1 {
        self.structure_receipt.id()
    }

    /// Complete normalized-structure identity receipt.
    #[must_use]
    pub const fn structure_identity_receipt(&self) -> IdentityReceipt<CausalStructureIdV1> {
        self.structure_receipt
    }

    /// Full generated/derived/escape provenance-bearing artifact identity.
    #[must_use]
    pub const fn artifact_identity(&self) -> CausalGraphArtifactIdV1 {
        self.artifact_receipt.id()
    }

    /// Complete provenance-artifact identity receipt.
    #[must_use]
    pub const fn artifact_identity_receipt(&self) -> IdentityReceipt<CausalGraphArtifactIdV1> {
        self.artifact_receipt
    }

    /// Exact admitted Machine graph identity.
    #[must_use]
    pub const fn machine_graph(&self) -> MachineGraphIdV1 {
        self.machine_graph.id()
    }

    /// Complete Machine-graph identity receipt retained for collision
    /// adjudication.
    #[must_use]
    pub const fn machine_graph_identity_receipt(&self) -> IdentityReceipt<MachineGraphIdV1> {
        self.machine_graph
    }

    /// Exact behavior overlay identity, absent only for an explicitly
    /// state-free graph-only admission.
    #[must_use]
    pub const fn machine_behavior(&self) -> Option<MachineBehaviorIdV1> {
        match self.machine_behavior {
            Some(receipt) => Some(receipt.id()),
            None => None,
        }
    }

    /// Complete behavior-overlay identity receipt, when graph admission was
    /// behavior-aware.
    #[must_use]
    pub const fn machine_behavior_identity_receipt(
        &self,
    ) -> Option<IdentityReceipt<MachineBehaviorIdV1>> {
        self.machine_behavior
    }

    /// Explicit unit convention of every structural signal.
    #[must_use]
    pub const fn units(&self) -> CausalUnitConvention {
        self.units
    }

    /// Complete-model or explicitly bounded partial extraction scope.
    #[must_use]
    pub const fn scope(&self) -> &CausalGraphScope {
        &self.scope
    }

    /// Exact producer and Five-Explicits extraction context.
    #[must_use]
    pub const fn extraction(&self) -> &CausalExtractionContext {
        &self.extraction
    }

    /// Canonically ordered equations.
    #[must_use]
    pub fn equations(&self) -> &[EquationSpec] {
        &self.equations
    }

    /// Canonically ordered variables.
    #[must_use]
    pub fn variables(&self) -> &[VariableSpec] {
        &self.variables
    }

    /// Canonically ordered activation-condition tables.
    #[must_use]
    pub fn conditions(&self) -> &[ActivationConditionSpec] {
        &self.conditions
    }

    /// Canonically ordered incidences.
    #[must_use]
    pub fn incidences(&self) -> &[IncidenceSpec] {
        &self.incidences
    }
}

/// Bounded graph-admission outcome for structured tracing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalGraphAdmissionDecision {
    submitted: CausalGraphSubmittedCounts,
    result: Result<AdmittedCausalGraph, CausalGraphRefusal>,
}

impl CausalGraphAdmissionDecision {
    /// Counts observed before canonicalization.
    #[must_use]
    pub const fn submitted_counts(&self) -> CausalGraphSubmittedCounts {
        self.submitted
    }

    /// Stable top-level decision code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "CausalGraphAdmitted",
            Err(_) => "CausalGraphRefused",
        }
    }

    /// Borrow the conventional result.
    #[must_use]
    pub fn result(&self) -> Result<&AdmittedCausalGraph, &CausalGraphRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision.
    #[must_use]
    pub fn into_result(self) -> Result<AdmittedCausalGraph, CausalGraphRefusal> {
        self.result
    }
}

impl CausalGraphDraft {
    /// Admit this structural draft against one exact admitted Machine graph.
    ///
    /// # Errors
    /// Returns every bounded deterministic finding and publishes no identity
    /// when ownership, lineage, units, clocks, resources, or canonical identity
    /// fail closed.
    pub fn admit_against(
        self,
        machine: &AdmittedMachineGraph,
        cx: &Cx<'_>,
    ) -> Result<AdmittedCausalGraph, CausalGraphRefusal> {
        self.admit_with_decision(machine, cx).into_result()
    }

    /// Admit while retaining pre-canonicalization resource counts.
    #[must_use]
    pub fn admit_with_decision(
        self,
        machine: &AdmittedMachineGraph,
        cx: &Cx<'_>,
    ) -> CausalGraphAdmissionDecision {
        match submitted_graph_counts(&self, cx) {
            Ok(submitted) => CausalGraphAdmissionDecision {
                submitted,
                result: admit_causal_graph(self, machine, None, submitted, cx),
            },
            Err((submitted, refusal)) => CausalGraphAdmissionDecision {
                submitted,
                result: Err(refusal),
            },
        }
    }

    /// Admit against an exact graph-plus-behavior chain.
    ///
    /// This is mandatory when any variable has [`VariableRole::State`].
    ///
    /// # Errors
    /// Refuses a behavior bound to another graph or a state contract mismatch,
    /// in addition to every ordinary graph resource, semantic, cancellation,
    /// and identity-admission refusal.
    pub fn admit_against_behavior(
        self,
        machine: &AdmittedMachineGraph,
        behavior: &AdmittedMachineBehavior,
        cx: &Cx<'_>,
    ) -> Result<AdmittedCausalGraph, CausalGraphRefusal> {
        self.admit_with_behavior_decision(machine, behavior, cx)
            .into_result()
    }

    /// Admit against graph plus behavior while retaining submitted counts.
    #[must_use]
    pub fn admit_with_behavior_decision(
        self,
        machine: &AdmittedMachineGraph,
        behavior: &AdmittedMachineBehavior,
        cx: &Cx<'_>,
    ) -> CausalGraphAdmissionDecision {
        match submitted_graph_counts(&self, cx) {
            Ok(submitted) => CausalGraphAdmissionDecision {
                submitted,
                result: admit_causal_graph(self, machine, Some(behavior), submitted, cx),
            },
            Err((submitted, refusal)) => CausalGraphAdmissionDecision {
                submitted,
                result: Err(refusal),
            },
        }
    }
}

fn submitted_graph_counts(
    draft: &CausalGraphDraft,
    cx: &Cx<'_>,
) -> Result<CausalGraphSubmittedCounts, (CausalGraphSubmittedCounts, CausalGraphRefusal)> {
    let mut counts = CausalGraphSubmittedCounts {
        complete: false,
        equations: draft.equations.len(),
        variables: draft.variables.len(),
        conditions: draft.conditions.len(),
        incidences: draft.incidences.len(),
        condition_dependencies: 0,
        condition_branches: 0,
        supports: 0,
        activation_selections: 0,
        activation_cubes: 0,
        derivation_references: 0,
        diagnostic_label_bytes: 0,
        max_derivation_parents: 0,
        max_activation_cubes_per_row: 0,
        max_activation_selections_per_cube: 0,
    };
    if let Err(refusal) = graph_checkpoint(cx) {
        return Err((counts, refusal));
    }
    if graph_counts_exceed_limits(&counts) {
        return Err((counts, resource_graph_refusal()));
    }
    let mut work = 0usize;
    for equation in &draft.equations {
        if let Err(refusal) = graph_poll(cx, work) {
            return Err((counts, refusal));
        }
        work = work.saturating_add(1);
        counts.diagnostic_label_bytes = counts
            .diagnostic_label_bytes
            .saturating_add(equation.diagnostic_label.len());
        counts.supports = counts.supports.saturating_add(equation.supports.len());
        let parents = derived_parent_count(&equation.lineage);
        counts.derivation_references = counts.derivation_references.saturating_add(parents);
        counts.max_derivation_parents = counts.max_derivation_parents.max(parents);
        if graph_counts_exceed_limits(&counts) {
            return Err((counts, resource_graph_refusal()));
        }
        if let Err(refusal) = count_activation(&equation.activation, &mut counts, &mut work, cx) {
            return Err((counts, refusal));
        }
    }
    for variable in &draft.variables {
        if let Err(refusal) = graph_poll(cx, work) {
            return Err((counts, refusal));
        }
        work = work.saturating_add(1);
        counts.diagnostic_label_bytes = counts
            .diagnostic_label_bytes
            .saturating_add(variable.diagnostic_label.len());
        counts.supports = counts.supports.saturating_add(variable.supports.len());
        let parents = derived_parent_count(&variable.lineage);
        counts.derivation_references = counts.derivation_references.saturating_add(parents);
        counts.max_derivation_parents = counts.max_derivation_parents.max(parents);
        if graph_counts_exceed_limits(&counts) {
            return Err((counts, resource_graph_refusal()));
        }
        if let Err(refusal) = count_activation(&variable.activation, &mut counts, &mut work, cx) {
            return Err((counts, refusal));
        }
    }
    for condition in &draft.conditions {
        if let Err(refusal) = graph_poll(cx, work) {
            return Err((counts, refusal));
        }
        work = work.saturating_add(1);
        counts.condition_dependencies = counts
            .condition_dependencies
            .saturating_add(condition.dependencies.len());
        counts.condition_branches = counts
            .condition_branches
            .saturating_add(condition.branches.len());
        if graph_counts_exceed_limits(&counts) {
            return Err((counts, resource_graph_refusal()));
        }
    }
    for incidence in &draft.incidences {
        if let Err(refusal) = graph_poll(cx, work) {
            return Err((counts, refusal));
        }
        work = work.saturating_add(1);
        if let Err(refusal) = count_activation(&incidence.activation, &mut counts, &mut work, cx) {
            return Err((counts, refusal));
        }
    }
    counts.complete = true;
    Ok(counts)
}

fn count_activation(
    activation: &ActivationDomain,
    counts: &mut CausalGraphSubmittedCounts,
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let cubes = activation_cubes(activation);
    counts.activation_cubes = counts.activation_cubes.saturating_add(cubes.len());
    counts.max_activation_cubes_per_row = counts.max_activation_cubes_per_row.max(cubes.len());
    if graph_counts_exceed_limits(counts) {
        return Err(resource_graph_refusal());
    }
    for cube in cubes {
        graph_poll(cx, *work)?;
        *work = (*work).saturating_add(1);
        counts.activation_selections = counts
            .activation_selections
            .saturating_add(cube.selections.len());
        counts.max_activation_selections_per_cube = counts
            .max_activation_selections_per_cube
            .max(cube.selections.len());
        if graph_counts_exceed_limits(counts) {
            return Err(resource_graph_refusal());
        }
    }
    Ok(())
}

fn graph_counts_exceed_limits(counts: &CausalGraphSubmittedCounts) -> bool {
    let label_envelope = counts
        .equations
        .saturating_add(counts.variables)
        .saturating_mul(MAX_CAUSAL_DIAGNOSTIC_LABEL_BYTES);
    counts.equations > MAX_CAUSAL_EQUATIONS
        || counts.variables > MAX_CAUSAL_VARIABLES
        || counts.conditions > MAX_CAUSAL_CONDITIONS
        || counts.incidences > MAX_CAUSAL_INCIDENCES
        || counts.condition_dependencies > MAX_CAUSAL_CONDITION_DEPENDENCIES
        || counts.condition_branches > MAX_CAUSAL_CONDITION_BRANCHES
        || counts.supports > MAX_CAUSAL_SUPPORT_REFERENCES
        || counts.activation_selections > MAX_CAUSAL_ACTIVATION_SELECTIONS
        || counts.activation_cubes > MAX_CAUSAL_ACTIVATION_CUBES
        || counts.derivation_references > MAX_CAUSAL_DERIVATION_REFERENCES
        || counts.max_derivation_parents > MAX_CAUSAL_DERIVATION_PARENTS
        || counts.max_activation_cubes_per_row > MAX_CAUSAL_CUBES_PER_ACTIVATION
        || counts.max_activation_selections_per_cube > MAX_CAUSAL_SELECTIONS_PER_CUBE
        || counts.diagnostic_label_bytes > label_envelope
}

fn graph_refusal(
    mut findings: Vec<CausalGraphFinding>,
    identity_error: Option<CanonicalError>,
) -> CausalGraphRefusal {
    findings.sort();
    findings.dedup();
    debug_assert!(!findings.is_empty());
    CausalGraphRefusal {
        findings,
        identity_error,
        identity_allocation_refusal: None,
        identity_stream_refusal: None,
    }
}

fn graph_refusal_cancellable(
    mut findings: Vec<CausalGraphFinding>,
    cx: &Cx<'_>,
) -> Result<CausalGraphRefusal, CausalGraphRefusal> {
    cancellable_sort(&mut findings, || graph_checkpoint(cx))?;
    cancellable_dedup(&mut findings, || graph_checkpoint(cx))?;
    debug_assert!(!findings.is_empty());
    Ok(CausalGraphRefusal {
        findings,
        identity_error: None,
        identity_allocation_refusal: None,
        identity_stream_refusal: None,
    })
}

fn cancelled_graph_refusal() -> CausalGraphRefusal {
    graph_refusal(
        vec![CausalGraphFinding::new(
            CausalGraphRule::Cancelled,
            CausalGraphSubject::Graph,
        )],
        None,
    )
}

fn resource_graph_refusal() -> CausalGraphRefusal {
    graph_refusal(
        vec![CausalGraphFinding::new(
            CausalGraphRule::ResourceLimit,
            CausalGraphSubject::Graph,
        )],
        None,
    )
}

fn enforce_graph_finding_budget(findings: &[CausalGraphFinding]) -> Result<(), CausalGraphRefusal> {
    if findings.len() > MAX_CAUSAL_GRAPH_FINDINGS {
        Err(resource_graph_refusal())
    } else {
        Ok(())
    }
}

fn identity_graph_refusal(error: CausalIdentityError) -> CausalGraphRefusal {
    let rule = if error.is_cancelled() {
        CausalGraphRule::Cancelled
    } else {
        CausalGraphRule::Identity
    };
    let mut refusal = graph_refusal(
        vec![CausalGraphFinding::new(rule, CausalGraphSubject::Graph)],
        error.canonical_error(),
    );
    refusal.identity_allocation_refusal = error.allocation_refusal();
    refusal.identity_stream_refusal = error.stream_refusal();
    refusal
}

fn graph_checkpoint(cx: &Cx<'_>) -> Result<(), CausalGraphRefusal> {
    cx.checkpoint().map_err(|_| cancelled_graph_refusal())
}

fn graph_poll(cx: &Cx<'_>, index: usize) -> Result<(), CausalGraphRefusal> {
    if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
        graph_checkpoint(cx)?;
    }
    Ok(())
}

#[derive(Debug)]
struct MachineView {
    clocks: BTreeSet<ClockId>,
    subsystems: BTreeSet<SubsystemId>,
    elements: BTreeMap<MachineElementId, SubsystemId>,
    terminal_ports: BTreeMap<TerminalId, PortId>,
    ports: BTreeMap<PortId, SubsystemId>,
    port_contracts: BTreeMap<(PortId, PortCoordinate), SignalContract>,
    relations: BTreeMap<RelationId, (SubsystemId, SubsystemId)>,
    interfaces: BTreeSet<InterfaceId>,
}

impl MachineView {
    #[allow(clippy::too_many_lines)]
    fn new(machine: &AdmittedMachineGraph, cx: &Cx<'_>) -> Result<Self, CausalGraphRefusal> {
        let mut clocks = BTreeSet::new();
        for (index, clock) in machine.clocks().iter().enumerate() {
            graph_poll(cx, index)?;
            clocks.insert(clock.id.clone());
        }
        let mut subsystems = BTreeSet::new();
        let mut elements = BTreeMap::new();
        let mut work = 0usize;
        for subsystem in machine.subsystems() {
            graph_poll(cx, work)?;
            work = work.saturating_add(1);
            subsystems.insert(subsystem.id.clone());
            for body in &subsystem.bodies {
                graph_poll(cx, work)?;
                work = work.saturating_add(1);
                elements.insert(MachineElementId::Body(body.clone()), subsystem.id.clone());
            }
            for patch in &subsystem.surface_patches {
                graph_poll(cx, work)?;
                work = work.saturating_add(1);
                elements.insert(
                    MachineElementId::SurfacePatch(patch.clone()),
                    subsystem.id.clone(),
                );
            }
            for feature in &subsystem.contact_features {
                graph_poll(cx, work)?;
                work = work.saturating_add(1);
                elements.insert(
                    MachineElementId::ContactFeature(feature.clone()),
                    subsystem.id.clone(),
                );
            }
            for state in &subsystem.state_slots {
                graph_poll(cx, work)?;
                work = work.saturating_add(1);
                elements.insert(
                    MachineElementId::StateSlot(state.clone()),
                    subsystem.id.clone(),
                );
            }
        }
        let mut terminals = BTreeMap::new();
        let mut terminal_specs = BTreeMap::new();
        for (index, terminal) in machine.terminals().iter().enumerate() {
            graph_poll(cx, index)?;
            terminals.insert(terminal.id.clone(), terminal.owner.clone());
            terminal_specs.insert(terminal.id.clone(), terminal);
            elements.insert(
                MachineElementId::Terminal(terminal.id.clone()),
                terminal.owner.clone(),
            );
        }
        let mut ports = BTreeMap::new();
        let mut terminal_ports = BTreeMap::new();
        let mut port_contracts = BTreeMap::new();
        for (index, port) in machine.ports().iter().enumerate() {
            graph_poll(cx, index)?;
            ports.insert(port.id.clone(), port.owner.clone());
            terminal_ports.insert(port.effort.clone(), port.id.clone());
            terminal_ports.insert(port.flow.clone(), port.id.clone());
            elements.insert(MachineElementId::Port(port.id.clone()), port.owner.clone());
            for (coordinate, terminal_id) in [
                (PortCoordinate::Effort, &port.effort),
                (PortCoordinate::Flow, &port.flow),
            ] {
                if let Some(terminal) = terminal_specs.get(terminal_id) {
                    port_contracts.insert(
                        (port.id.clone(), coordinate),
                        SignalContract {
                            quantity: terminal.quantity,
                            shape: terminal.shape,
                            clock: terminal.clock.clone(),
                            frame: terminal.frame.clone(),
                        },
                    );
                }
            }
        }
        let mut relations = BTreeMap::new();
        for (index, relation) in machine.relations().iter().enumerate() {
            graph_poll(cx, index)?;
            if let (Some(source), Some(target)) = (
                terminals.get(&relation.source),
                terminals.get(&relation.target),
            ) {
                relations.insert(relation.id.clone(), (source.clone(), target.clone()));
            }
        }
        let mut interfaces = BTreeSet::new();
        for (index, interface) in machine.interfaces().iter().enumerate() {
            graph_poll(cx, index)?;
            interfaces.insert(interface.id.clone());
        }
        graph_checkpoint(cx)?;
        Ok(Self {
            clocks,
            subsystems,
            elements,
            terminal_ports,
            ports,
            port_contracts,
            relations,
            interfaces,
        })
    }

    fn owner_exists(&self, owner: &CausalOwner) -> bool {
        match owner {
            CausalOwner::Subsystem(id) => self.subsystems.contains(id),
            CausalOwner::Port(id) => self.ports.contains_key(id),
            CausalOwner::Interface(id) => self.interfaces.contains(id),
        }
    }

    fn origin_exists(&self, origin: &NodeOrigin) -> bool {
        match origin {
            NodeOrigin::AuditedEscapeHatch(_) | NodeOrigin::Derived(_) => true,
            NodeOrigin::Machine(machine) => match machine {
                MachineNodeOrigin::Subsystem(id) => self.subsystems.contains(id),
                MachineNodeOrigin::Relation(id) => self.relations.contains_key(id),
                MachineNodeOrigin::Element(id) => self.elements.contains_key(id),
                MachineNodeOrigin::PortCoordinate { port, .. } => self.ports.contains_key(port),
                MachineNodeOrigin::Interface(id) => self.interfaces.contains(id),
            },
        }
    }

    fn owner_matches_origin(&self, owner: &CausalOwner, origin: &NodeOrigin) -> bool {
        match origin {
            NodeOrigin::AuditedEscapeHatch(_) | NodeOrigin::Derived(_) => true,
            NodeOrigin::Machine(MachineNodeOrigin::Subsystem(id)) => {
                owner == &CausalOwner::Subsystem(id.clone())
            }
            NodeOrigin::Machine(MachineNodeOrigin::Relation(id)) => {
                self.relations.get(id).is_some_and(|(source, target)| {
                    owner == &CausalOwner::Subsystem(source.clone())
                        || owner == &CausalOwner::Subsystem(target.clone())
                })
            }
            NodeOrigin::Machine(MachineNodeOrigin::Element(element)) => {
                let subsystem_match = self
                    .elements
                    .get(element)
                    .is_some_and(|subsystem| owner == &CausalOwner::Subsystem(subsystem.clone()));
                let port_match = match element {
                    MachineElementId::Port(port) => owner == &CausalOwner::Port(port.clone()),
                    MachineElementId::Terminal(terminal) => self
                        .terminal_ports
                        .get(terminal)
                        .is_some_and(|port| owner == &CausalOwner::Port(port.clone())),
                    _ => false,
                };
                subsystem_match || port_match
            }
            NodeOrigin::Machine(MachineNodeOrigin::PortCoordinate { port, .. }) => {
                owner == &CausalOwner::Port(port.clone())
            }
            NodeOrigin::Machine(MachineNodeOrigin::Interface(id)) => {
                owner == &CausalOwner::Interface(id.clone())
            }
        }
    }

    fn support_exists(&self, support: &CausalSupport) -> bool {
        match support {
            CausalSupport::Lumped | CausalSupport::External(_) => true,
            CausalSupport::MachineElement(element) => self.elements.contains_key(element),
        }
    }

    fn port_contract(&self, port: &PortId, coordinate: PortCoordinate) -> Option<&SignalContract> {
        self.port_contracts.get(&(port.clone(), coordinate))
    }
}

fn equation_identity_matches_lineage(id: &EquationId, lineage: &NodeLineage) -> bool {
    EquationId::derive(lineage).is_ok_and(|derived| {
        identity_receipt_adjudication_eq(derived.identity_receipt(), id.identity_receipt())
    })
}

fn variable_identity_matches_lineage(id: &VariableId, lineage: &NodeLineage) -> bool {
    VariableId::derive(lineage).is_ok_and(|derived| {
        identity_receipt_adjudication_eq(derived.identity_receipt(), id.identity_receipt())
    })
}

#[allow(clippy::too_many_lines)]
fn admit_causal_graph(
    mut draft: CausalGraphDraft,
    machine: &AdmittedMachineGraph,
    behavior: Option<&AdmittedMachineBehavior>,
    counts: CausalGraphSubmittedCounts,
    cx: &Cx<'_>,
) -> Result<AdmittedCausalGraph, CausalGraphRefusal> {
    graph_checkpoint(cx)?;
    debug_assert!(counts.complete);
    if graph_counts_exceed_limits(&counts) {
        return Err(resource_graph_refusal());
    }

    graph_checkpoint(cx)?;
    canonicalize_graph_draft(&mut draft, cx)?;
    graph_checkpoint(cx)?;
    let view = MachineView::new(machine, cx)?;
    let mut findings = Vec::new();

    let behavior_contracts = if let Some(overlay) = behavior {
        let mut contracts = BTreeMap::new();
        for (index, contract) in overlay.state_contracts().iter().enumerate() {
            graph_poll(cx, index)?;
            contracts.insert(contract.id.clone(), contract);
        }
        graph_checkpoint(cx)?;
        Some(contracts)
    } else {
        None
    };
    if behavior.is_some_and(|overlay| overlay.base_graph() != machine.identity()) {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::StateBehaviorMismatch,
            CausalGraphSubject::Graph,
        ));
    }

    for (index, pair) in draft.equations.windows(2).enumerate() {
        graph_poll(cx, index)?;
        if pair[0].id.identity() == pair[1].id.identity() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateEquation,
                CausalGraphSubject::Equation(pair[1].id.clone()),
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }
    for (index, pair) in draft.variables.windows(2).enumerate() {
        graph_poll(cx, index)?;
        if pair[0].id.identity() == pair[1].id.identity() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateVariable,
                CausalGraphSubject::Variable(pair[1].id.clone()),
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }
    for (index, pair) in draft.conditions.windows(2).enumerate() {
        graph_poll(cx, index)?;
        if pair[0].condition == pair[1].condition {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateCondition,
                CausalGraphSubject::Condition(pair[1].condition.clone()),
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }
    for (index, pair) in draft.incidences.windows(2).enumerate() {
        graph_poll(cx, index)?;
        if pair[0].id.identity() == pair[1].id.identity() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateIncidence,
                CausalGraphSubject::IncidenceIdentity(pair[1].id.clone()),
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }
    // Duplicate nominal identities make every downstream lookup ambiguous.
    // Refuse that class before building ID-keyed maps so caller order cannot
    // choose which conflicting payload subsequent diagnostics observe.
    if !findings.is_empty() {
        return Err(graph_refusal_cancellable(findings, cx)?);
    }

    let mut equation_ids = BTreeSet::new();
    let mut equation_specs = BTreeMap::new();
    for (index, equation) in draft.equations.iter().enumerate() {
        graph_poll(cx, index)?;
        equation_ids.insert(equation.id.clone());
        equation_specs.insert(equation.id.clone(), equation);
    }
    let mut variable_ids = BTreeSet::new();
    let mut variable_specs = BTreeMap::new();
    for (index, variable) in draft.variables.iter().enumerate() {
        graph_poll(cx, index)?;
        variable_ids.insert(variable.id.clone());
        variable_specs.insert(variable.id.clone(), variable);
    }
    let mut condition_ids = BTreeSet::new();
    let mut condition_branches = BTreeSet::new();
    let mut condition_domains = BTreeMap::new();
    for (condition_index, condition) in draft.conditions.iter().enumerate() {
        graph_poll(cx, condition_index)?;
        condition_ids.insert(condition.condition.clone());
        let mut branches = BTreeSet::new();
        for (branch_index, branch) in condition.branches.iter().enumerate() {
            graph_poll(cx, branch_index)?;
            branches.insert(branch.clone());
            condition_branches.insert((condition.condition.clone(), branch.clone()));
        }
        condition_domains.insert(condition.condition.clone(), branches);
    }

    for (index, equation) in draft.equations.iter().enumerate() {
        graph_poll(cx, index)?;
        let subject = CausalGraphSubject::Equation(equation.id.clone());
        if !equation_identity_matches_lineage(&equation.id, &equation.lineage) {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::EquationIdentityMismatch,
                subject.clone(),
            ));
        }
        validate_label(&equation.diagnostic_label, &subject, &mut findings);
        validate_node_common(
            &equation.lineage,
            &equation.owner,
            &equation.supports,
            &equation.residual,
            &equation.activation,
            &condition_branches,
            &view,
            &subject,
            &mut findings,
            cx,
        )?;
        enforce_graph_finding_budget(&findings)?;
    }

    let mut exclusive_role_bindings = BTreeMap::<MachineNodeOrigin, VariableId>::new();
    for (index, variable) in draft.variables.iter().enumerate() {
        graph_poll(cx, index)?;
        let subject = CausalGraphSubject::Variable(variable.id.clone());
        if !variable_identity_matches_lineage(&variable.id, &variable.lineage) {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::VariableIdentityMismatch,
                subject.clone(),
            ));
        }
        validate_label(&variable.diagnostic_label, &subject, &mut findings);
        validate_node_common(
            &variable.lineage,
            &variable.owner,
            &variable.supports,
            &variable.value,
            &variable.activation,
            &condition_branches,
            &view,
            &subject,
            &mut findings,
            cx,
        )?;
        validate_variable_role(
            variable,
            behavior_contracts.as_ref(),
            &view,
            &subject,
            &mut findings,
        );
        let exclusive_origin = match (&variable.role, &variable.lineage.origin) {
            (
                VariableRole::State | VariableRole::PortEffort | VariableRole::PortFlow,
                NodeOrigin::Machine(origin),
            ) => Some(origin),
            _ => None,
        };
        if let Some(origin) = exclusive_origin
            && let Some(previous) =
                exclusive_role_bindings.insert(origin.clone(), variable.id.clone())
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateRoleBinding,
                CausalGraphSubject::Variable(previous),
            ));
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DuplicateRoleBinding,
                subject,
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }

    validate_derived_lineage(&draft, &equation_ids, &variable_ids, &mut findings, cx)?;
    enforce_graph_finding_budget(&findings)?;
    let mut equation_incidence_dependencies = BTreeMap::<EquationId, BTreeSet<VariableId>>::new();
    let mut equation_incidences_always_available = BTreeMap::<EquationId, bool>::new();
    for (index, incidence) in draft.incidences.iter().enumerate() {
        graph_poll(cx, index)?;
        equation_incidence_dependencies
            .entry(incidence.equation.clone())
            .or_default()
            .insert(incidence.variable.clone());
        equation_incidences_always_available
            .entry(incidence.equation.clone())
            .and_modify(|available| {
                *available &= incidence.activation == ActivationDomain::Always;
            })
            .or_insert(incidence.activation == ActivationDomain::Always);
    }
    validate_conditions(
        &draft,
        &equation_specs,
        &equation_incidence_dependencies,
        &equation_incidences_always_available,
        &variable_specs,
        &condition_ids,
        &condition_branches,
        &mut findings,
        cx,
    )?;
    enforce_graph_finding_budget(&findings)?;
    let mut derivative_participation = BTreeMap::<
        DerivativeVariableKey,
        BTreeMap<SolveParticipation, Vec<&ActivationDomain>>,
    >::new();
    let mut resolving_domain_indices = BTreeMap::<VariableId, Vec<usize>>::new();
    let mut activation_proof_states = 0usize;
    for (index, incidence) in draft.incidences.iter().enumerate() {
        graph_poll(cx, index)?;
        if !incidence_identity_matches_cancellable(incidence, cx)? {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::IncidenceIdentityMismatch,
                incidence_subject(incidence),
            ));
        }
        let matching_equation = equation_specs
            .get(&incidence.equation)
            .is_some_and(|equation| {
                equation.solve_participation == EquationParticipation::Matching
            });
        if matching_equation
            && matches!(
                incidence.solve_participation,
                SolveParticipation::Unknown | SolveParticipation::KnownInput
            )
        {
            resolving_domain_indices
                .entry(incidence.variable.clone())
                .or_default()
                .push(index);
        }
        let participation_conflict = if matching_equation
            && incidence.solve_participation != SolveParticipation::ConditionOnly
        {
            let derivative_key = DerivativeVariableKey {
                variable: incidence.variable.clone(),
                derivative_order: incidence.derivative_order,
            };
            let cases_by_participation =
                derivative_participation.entry(derivative_key).or_default();
            let mut conflict = false;
            for (participation, cases) in cases_by_participation.iter() {
                if *participation == incidence.solve_participation {
                    continue;
                }
                for activation in cases {
                    if activation_domains_overlap_bounded(
                        activation,
                        &incidence.activation,
                        &mut activation_proof_states,
                        cx,
                    )? {
                        conflict = true;
                        break;
                    }
                }
                if conflict {
                    break;
                }
            }
            cases_by_participation
                .entry(incidence.solve_participation)
                .or_default()
                .push(&incidence.activation);
            conflict
        } else {
            false
        };
        if incidence.solve_participation == SolveParticipation::ModeDependent
            || participation_conflict
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DerivativeParticipationMismatch,
                incidence_subject(incidence),
            ));
        }
        validate_incidence(
            incidence,
            &equation_specs,
            &variable_specs,
            &condition_branches,
            &condition_domains,
            &view,
            &mut findings,
            &mut activation_proof_states,
            cx,
        )?;
        enforce_graph_finding_budget(&findings)?;
    }

    for (index, variable) in draft.variables.iter().enumerate() {
        graph_poll(cx, index)?;
        if variable.solve_participation != SolveParticipation::ModeDependent {
            continue;
        }
        let resolving_domains = resolving_domain_indices
            .get(&variable.id)
            .into_iter()
            .flatten()
            .map(|&incidence_index| &draft.incidences[incidence_index].activation);
        if !activation_implies_union(
            &variable.activation,
            resolving_domains,
            &condition_domains,
            &mut activation_proof_states,
            cx,
        )? {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::DerivativeParticipationMismatch,
                CausalGraphSubject::Variable(variable.id.clone()),
            ));
        }
        enforce_graph_finding_budget(&findings)?;
    }

    let mut possible_unknown_vertices = BTreeSet::new();
    let mut always_unknown_state_derivatives = BTreeSet::new();
    for (index, incidence) in draft.incidences.iter().enumerate() {
        graph_poll(cx, index)?;
        if incidence.solve_participation != SolveParticipation::Unknown {
            continue;
        }
        if possible_unknown_vertices.insert(DerivativeVariableKey {
            variable: incidence.variable.clone(),
            derivative_order: incidence.derivative_order,
        }) && possible_unknown_vertices.len() > MAX_CAUSAL_DERIVATIVE_VERTICES
        {
            return Err(resource_graph_refusal());
        }
        if incidence.derivative_order > 0 && incidence.activation == ActivationDomain::Always {
            always_unknown_state_derivatives.insert(incidence.variable.clone());
        }
    }
    for (index, variable) in draft.variables.iter().enumerate() {
        graph_poll(cx, index)?;
        if variable.solve_participation == SolveParticipation::Unknown
            && (variable.role != VariableRole::State
                || !always_unknown_state_derivatives.contains(&variable.id))
        {
            if possible_unknown_vertices.insert(DerivativeVariableKey {
                variable: variable.id.clone(),
                derivative_order: 0,
            }) && possible_unknown_vertices.len() > MAX_CAUSAL_DERIVATIVE_VERTICES
            {
                return Err(resource_graph_refusal());
            }
        }
    }

    if !findings.is_empty() {
        return Err(graph_refusal_cancellable(findings, cx)?);
    }

    let behavior_identity = behavior.map(AdmittedMachineBehavior::identity_receipt);
    graph_checkpoint(cx)?;
    let structure_receipt = match causal_structure_identity(&draft, machine, cx) {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(identity_graph_refusal(error));
        }
    };
    let artifact_receipt =
        match causal_artifact_identity(&draft, structure_receipt, behavior_identity, cx) {
            Ok(receipt) => receipt,
            Err(error) => {
                return Err(identity_graph_refusal(error));
            }
        };
    graph_checkpoint(cx)?;

    Ok(AdmittedCausalGraph {
        machine_graph: machine.identity_receipt(),
        machine_behavior: behavior_identity,
        units: draft.units,
        scope: draft.scope,
        extraction: draft.extraction,
        equations: draft.equations,
        variables: draft.variables,
        conditions: draft.conditions,
        incidences: draft.incidences,
        structure_receipt,
        artifact_receipt,
    })
}

fn derived_parent_count(lineage: &NodeLineage) -> usize {
    match &lineage.origin {
        NodeOrigin::Derived(derived) => derived.parents.len(),
        NodeOrigin::Machine(_) | NodeOrigin::AuditedEscapeHatch(_) => 0,
    }
}

fn canonicalize_graph_draft(
    draft: &mut CausalGraphDraft,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let mut work = 0usize;
    for equation in &mut draft.equations {
        graph_poll(cx, work)?;
        work = work.saturating_add(1);
        cancellable_sort(&mut equation.supports, || graph_checkpoint(cx))?;
        canonicalize_lineage_cancellable(&mut equation.lineage, cx)?;
        canonicalize_activation_cancellable(&mut equation.activation, cx)?;
    }
    for variable in &mut draft.variables {
        graph_poll(cx, work)?;
        work = work.saturating_add(1);
        cancellable_sort(&mut variable.supports, || graph_checkpoint(cx))?;
        canonicalize_lineage_cancellable(&mut variable.lineage, cx)?;
        canonicalize_activation_cancellable(&mut variable.activation, cx)?;
    }
    for condition in &mut draft.conditions {
        graph_poll(cx, work)?;
        work = work.saturating_add(1);
        cancellable_sort(&mut condition.branches, || graph_checkpoint(cx))?;
        cancellable_sort(&mut condition.dependencies, || graph_checkpoint(cx))?;
    }
    for incidence in &mut draft.incidences {
        graph_poll(cx, work)?;
        work = work.saturating_add(1);
        canonicalize_activation_cancellable(&mut incidence.activation, cx)?;
    }
    cancellable_sort_by(
        &mut draft.equations,
        |left, right| left.id.cmp(&right.id),
        || graph_checkpoint(cx),
    )?;
    cancellable_sort_by(
        &mut draft.variables,
        |left, right| left.id.cmp(&right.id),
        || graph_checkpoint(cx),
    )?;
    cancellable_sort_by(
        &mut draft.conditions,
        |left, right| left.condition.cmp(&right.condition),
        || graph_checkpoint(cx),
    )?;
    cancellable_sort_by(
        &mut draft.incidences,
        |left, right| left.id.cmp(&right.id),
        || graph_checkpoint(cx),
    )?;
    graph_checkpoint(cx)
}

fn canonicalize_lineage_cancellable(
    lineage: &mut NodeLineage,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    if let NodeOrigin::Derived(derived) = &mut lineage.origin {
        cancellable_sort(&mut derived.parents, || graph_checkpoint(cx))?;
    }
    Ok(())
}

fn canonicalize_activation_for_identity(
    activation: &mut ActivationDomain,
    cx: &Cx<'_>,
) -> Result<(), CanonicalError> {
    if let ActivationDomain::Conditional { cubes } = activation {
        for (index, cube) in cubes.iter_mut().enumerate() {
            identity_materialization_poll(cx, index, 0)?;
            cancellable_sort(&mut cube.selections, || {
                identity_materialization_checkpoint(cx, 0)
            })?;
        }
        cancellable_sort_by_fallible(
            cubes,
            |left, right| {
                compare_activation_cubes_cancellable(left, right, || {
                    identity_materialization_checkpoint(cx, 0)
                })
            },
            || identity_materialization_checkpoint(cx, 0),
        )?;
    }
    identity_materialization_checkpoint(cx, 0)
}

fn canonicalize_activation_cancellable(
    activation: &mut ActivationDomain,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    if let ActivationDomain::Conditional { cubes } = activation {
        for (index, cube) in cubes.iter_mut().enumerate() {
            graph_poll(cx, index)?;
            cancellable_sort(&mut cube.selections, || graph_checkpoint(cx))?;
        }
        cancellable_sort_by_fallible(
            cubes,
            |left, right| {
                compare_activation_cubes_cancellable(left, right, || graph_checkpoint(cx))
            },
            || graph_checkpoint(cx),
        )?;
    }
    Ok(())
}

fn cancellable_sort<T: Ord, E>(
    values: &mut [T],
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    cancellable_sort_by(values, |left, right| left.cmp(right), checkpoint)
}

fn cancellable_sort_by<T, E>(
    values: &mut [T],
    mut compare: impl FnMut(&T, &T) -> core::cmp::Ordering,
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    cancellable_sort_by_fallible(values, |left, right| Ok(compare(left, right)), checkpoint)
}

fn cancellable_sort_by_fallible<T, E>(
    values: &mut [T],
    compare: impl FnMut(&T, &T) -> Result<core::cmp::Ordering, E>,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    cancellable_sort_by_fallible_observed(values, compare, |_| checkpoint())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CausalSortPhase {
    Entry,
    IndexInitialization,
    Merge,
    InverseInitialization,
    InverseMap,
    PayloadPosition,
    PayloadSwap,
    Complete,
}

fn cancellable_sort_by_fallible_observed<T, E>(
    values: &mut [T],
    mut compare: impl FnMut(&T, &T) -> Result<core::cmp::Ordering, E>,
    mut checkpoint: impl FnMut(CausalSortPhase) -> Result<(), E>,
) -> Result<(), E> {
    checkpoint(CausalSortPhase::Entry)?;
    if values.len() < 2 {
        return Ok(());
    }

    let mut order = Vec::with_capacity(values.len());
    let mut buffer = Vec::with_capacity(values.len());
    for index in 0..values.len() {
        if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint(CausalSortPhase::IndexInitialization)?;
        }
        order.push(index);
        buffer.push(0usize);
    }
    {
        let immutable_values: &[T] = values;
        let mut width = 1usize;
        let mut merge_work = 0usize;
        while width < order.len() {
            let run_width = width.saturating_mul(2);
            for start in (0..order.len()).step_by(run_width) {
                let middle = start.saturating_add(width).min(order.len());
                let end = start.saturating_add(run_width).min(order.len());
                let mut left = start;
                let mut right = middle;
                let mut output = start;
                while output < end {
                    if merge_work.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
                        checkpoint(CausalSortPhase::Merge)?;
                    }
                    let take_left = right == end
                        || (left < middle
                            && compare(
                                &immutable_values[order[left]],
                                &immutable_values[order[right]],
                            )?
                            .then_with(|| order[left].cmp(&order[right]))
                                != core::cmp::Ordering::Greater);
                    if take_left {
                        buffer[output] = order[left];
                        left += 1;
                    } else {
                        buffer[output] = order[right];
                        right += 1;
                    }
                    output += 1;
                    merge_work = merge_work.saturating_add(1);
                }
            }
            core::mem::swap(&mut order, &mut buffer);
            width = run_width;
        }
    }

    let mut target_position = Vec::with_capacity(order.len());
    for index in 0..order.len() {
        if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint(CausalSortPhase::InverseInitialization)?;
        }
        target_position.push(0usize);
    }
    for (new_position, old_position) in order.into_iter().enumerate() {
        if new_position.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint(CausalSortPhase::InverseMap)?;
        }
        target_position[old_position] = new_position;
    }
    let mut swaps = 0usize;
    for position in 0..values.len() {
        if position.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint(CausalSortPhase::PayloadPosition)?;
        }
        while target_position[position] != position {
            if swaps.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
                checkpoint(CausalSortPhase::PayloadSwap)?;
            }
            let target = target_position[position];
            values.swap(position, target);
            target_position.swap(position, target);
            swaps = swaps.saturating_add(1);
        }
    }
    checkpoint(CausalSortPhase::Complete)
}

fn cancellable_lexicographic_cmp<T: Ord, E>(
    left: &[T],
    right: &[T],
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<core::cmp::Ordering, E> {
    checkpoint()?;
    for (index, (left_item, right_item)) in left.iter().zip(right).enumerate() {
        if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint()?;
        }
        let ordering = left_item.cmp(right_item);
        if ordering != core::cmp::Ordering::Equal {
            return Ok(ordering);
        }
    }
    Ok(left.len().cmp(&right.len()))
}

fn cancellable_slice_eq<T: PartialEq, E>(
    left: &[T],
    right: &[T],
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<bool, E> {
    checkpoint()?;
    if left.len() != right.len() {
        return Ok(false);
    }
    for (index, (left_item, right_item)) in left.iter().zip(right).enumerate() {
        if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint()?;
        }
        if left_item != right_item {
            return Ok(false);
        }
    }
    Ok(true)
}

fn cancellable_set_eq<T: Ord, E>(
    left: &BTreeSet<T>,
    right: &BTreeSet<T>,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<bool, E> {
    checkpoint()?;
    if left.len() != right.len() {
        return Ok(false);
    }
    for (index, (left_item, right_item)) in left.iter().zip(right).enumerate() {
        if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint()?;
        }
        if left_item != right_item {
            return Ok(false);
        }
    }
    Ok(true)
}

fn compare_activation_cubes_cancellable<E>(
    left: &ActivationCube,
    right: &ActivationCube,
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<core::cmp::Ordering, E> {
    cancellable_lexicographic_cmp(&left.selections, &right.selections, checkpoint)
}

fn compare_conditional_outcomes_cancellable<E>(
    left: &ConditionalCausalOutcome,
    right: &ConditionalCausalOutcome,
    checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<core::cmp::Ordering, E> {
    let assignment =
        cancellable_lexicographic_cmp(&left.assignment, &right.assignment, checkpoint)?;
    Ok(assignment
        .then_with(|| identity_receipt_adjudication_cmp(left.structure, right.structure))
        .then_with(|| identity_receipt_adjudication_cmp(left.artifact, right.artifact))
        .then_with(|| left.determination.cmp(&right.determination))
        .then_with(|| left.structural_rank.cmp(&right.structural_rank))
        .then_with(|| left.unknown_axes.cmp(&right.unknown_axes))
        .then_with(|| identity_receipt_adjudication_cmp(left.outcome, right.outcome))
        .then_with(|| identity_receipt_adjudication_cmp(left.receipt, right.receipt)))
}

fn cancellable_dedup<T: PartialEq, E>(
    values: &mut Vec<T>,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<(), E> {
    checkpoint()?;
    if values.len() < 2 {
        return Ok(());
    }
    let mut write = 1usize;
    for read in 1..values.len() {
        if read.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
            checkpoint()?;
        }
        if values[read] != values[write - 1] {
            values.swap(write, read);
            write += 1;
        }
    }
    values.truncate(write);
    checkpoint()
}

fn validate_label(
    label: &str,
    subject: &CausalGraphSubject,
    findings: &mut Vec<CausalGraphFinding>,
) {
    if label.is_empty()
        || label.len() > MAX_CAUSAL_DIAGNOSTIC_LABEL_BYTES
        || label.chars().any(char::is_control)
    {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidDiagnosticLabel,
            subject.clone(),
        ));
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_node_common(
    lineage: &NodeLineage,
    owner: &CausalOwner,
    supports: &[CausalSupport],
    signal: &SignalContract,
    activation: &ActivationDomain,
    condition_branches: &BTreeSet<(ActivationConditionRef, ActivationBranchRef)>,
    view: &MachineView,
    subject: &CausalGraphSubject,
    findings: &mut Vec<CausalGraphFinding>,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    if !view.owner_exists(owner) {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownOwner,
            CausalGraphSubject::Owner(owner.clone()),
        ));
    }
    if lineage.instance != *owner
        || (view.origin_exists(&lineage.origin)
            && !view.owner_matches_origin(owner, &lineage.origin))
    {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::OwnerOriginMismatch,
            subject.clone(),
        ));
    }
    if !view.origin_exists(&lineage.origin) {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownOrigin,
            subject.clone(),
        ));
    }
    if let NodeOrigin::AuditedEscapeHatch(escape) = &lineage.origin
        && escape.source != escape.audited_source
    {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::EscapeAuditMismatch,
            subject.clone(),
        ));
    }
    validate_supports(supports, view, subject, findings, cx)?;
    if !signal.quantity.is_admitted() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnsupportedQuantityForm,
            subject.clone(),
        ));
    }
    if !view.clocks.contains(&signal.clock) {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownClock,
            subject.clone(),
        ));
    }
    if !activation_is_valid(activation, condition_branches, cx)? {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidActivationCondition,
            subject.clone(),
        ));
    }
    Ok(())
}

fn validate_supports(
    supports: &[CausalSupport],
    view: &MachineView,
    subject: &CausalGraphSubject,
    findings: &mut Vec<CausalGraphFinding>,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let mut duplicate = false;
    let mut invalid_lumped = false;
    let mut unknown_support = false;
    for (index, support) in supports.iter().enumerate() {
        graph_poll(cx, index)?;
        duplicate |= index > 0 && &supports[index - 1] == support;
        invalid_lumped |= supports.len() > 1 && *support == CausalSupport::Lumped;
        unknown_support |= !view.support_exists(support);
    }
    if supports.is_empty() || duplicate || invalid_lumped {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidSupportSet,
            subject.clone(),
        ));
    }
    if unknown_support {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownSupport,
            subject.clone(),
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn validate_variable_role(
    variable: &VariableSpec,
    behavior_contracts: Option<&BTreeMap<StateSlotId, &StateSlotContract>>,
    view: &MachineView,
    subject: &CausalGraphSubject,
    findings: &mut Vec<CausalGraphFinding>,
) {
    let expected_port_coordinate = match variable.role {
        VariableRole::PortEffort => Some(PortCoordinate::Effort),
        VariableRole::PortFlow => Some(PortCoordinate::Flow),
        _ => None,
    };
    if let Some(expected) = expected_port_coordinate {
        let port = match &variable.lineage.origin {
            NodeOrigin::Machine(MachineNodeOrigin::PortCoordinate { port, coordinate })
                if *coordinate == expected =>
            {
                Some(port)
            }
            _ => None,
        };
        if port.is_none() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::RoleOriginMismatch,
                subject.clone(),
            ));
        }
        if !port.is_some_and(|port| view.port_contract(port, expected) == Some(&variable.value)) {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::RoleOriginMismatch,
                subject.clone(),
            ));
        }
        if variable.port_schema_crosswalk.is_none() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::PortCrosswalkMismatch,
                subject.clone(),
            ));
        }
    } else if variable.port_schema_crosswalk.is_some() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::PortCrosswalkMismatch,
            subject.clone(),
        ));
    }
    if variable
        .port_schema_crosswalk
        .as_ref()
        .is_some_and(|binding| binding.projection != binding.audited_projection)
    {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::PortCrosswalkMismatch,
            subject.clone(),
        ));
    }

    match &variable.lineage.origin {
        NodeOrigin::Machine(MachineNodeOrigin::PortCoordinate { coordinate, .. }) => {
            let expected_role = match coordinate {
                PortCoordinate::Effort => VariableRole::PortEffort,
                PortCoordinate::Flow => VariableRole::PortFlow,
            };
            if variable.role != expected_role {
                findings.push(CausalGraphFinding::new(
                    CausalGraphRule::RoleOriginMismatch,
                    subject.clone(),
                ));
            }
        }
        NodeOrigin::Machine(MachineNodeOrigin::Element(MachineElementId::StateSlot(_)))
            if variable.role != VariableRole::State =>
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::RoleOriginMismatch,
                subject.clone(),
            ));
        }
        _ => {}
    }

    let participation_admitted = match variable.role {
        VariableRole::Algebraic => matches!(
            variable.solve_participation,
            SolveParticipation::Unknown | SolveParticipation::ModeDependent
        ),
        VariableRole::State => matches!(
            variable.solve_participation,
            SolveParticipation::Unknown
                | SolveParticipation::KnownInput
                | SolveParticipation::ModeDependent
        ),
        VariableRole::Source => variable.solve_participation == SolveParticipation::KnownInput,
        VariableRole::Parameter | VariableRole::DiscreteMode => {
            matches!(
                variable.solve_participation,
                SolveParticipation::KnownInput | SolveParticipation::ConditionOnly
            )
        }
        VariableRole::PortEffort | VariableRole::PortFlow => {
            matches!(
                variable.solve_participation,
                SolveParticipation::Unknown
                    | SolveParticipation::KnownInput
                    | SolveParticipation::ModeDependent
            )
        }
    };
    if !participation_admitted {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::SolveParticipationMismatch,
            subject.clone(),
        ));
    }

    if variable.role == VariableRole::State {
        let state_id = match &variable.lineage.origin {
            NodeOrigin::Machine(MachineNodeOrigin::Element(MachineElementId::StateSlot(id))) => {
                Some(id)
            }
            _ => None,
        };
        let contract = state_id.and_then(|id| behavior_contracts.and_then(|map| map.get(id)));
        let contract_matches = contract.is_some_and(|contract| {
            variable.owner == CausalOwner::Subsystem(contract.owner.clone())
                && variable.value
                    == (SignalContract {
                        quantity: contract.quantity,
                        shape: contract.shape,
                        clock: contract.clock.clone(),
                        frame: contract.frame.clone(),
                    })
        });
        if state_id.is_none() {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::RoleOriginMismatch,
                subject.clone(),
            ));
        }
        if !contract_matches {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::StateBehaviorMismatch,
                subject.clone(),
            ));
        }
    }
}

#[allow(clippy::too_many_lines)]
fn validate_derived_lineage(
    draft: &CausalGraphDraft,
    equation_ids: &BTreeSet<EquationId>,
    variable_ids: &BTreeSet<VariableId>,
    findings: &mut Vec<CausalGraphFinding>,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let mut derived_nodes = BTreeSet::new();
    let mut edges = Vec::<(ParentNodeId, ParentNodeId)>::new();
    let mut parents_by_child = BTreeMap::<ParentNodeId, Vec<ParentNodeId>>::new();
    let mut local_orders = BTreeMap::<ParentNodeId, u16>::new();
    let mut subjects = BTreeMap::<ParentNodeId, CausalGraphSubject>::new();
    for (index, (node, lineage, subject)) in draft
        .equations
        .iter()
        .map(|equation| {
            (
                ParentNodeId::Equation(equation.id.clone()),
                &equation.lineage,
                CausalGraphSubject::Equation(equation.id.clone()),
            )
        })
        .chain(draft.variables.iter().map(|variable| {
            (
                ParentNodeId::Variable(variable.id.clone()),
                &variable.lineage,
                CausalGraphSubject::Variable(variable.id.clone()),
            )
        }))
        .enumerate()
    {
        graph_poll(cx, index)?;
        let NodeOrigin::Derived(derived) = &lineage.origin else {
            continue;
        };
        derived_nodes.insert(node.clone());
        local_orders.insert(node.clone(), derived.differentiation_order);
        subjects.insert(node.clone(), subject.clone());
        let mut canonical_parents = Vec::with_capacity(derived.parents.len());
        let mut duplicate = false;
        let mut dangling = false;
        let mut self_parent = false;
        for (parent_index, parent) in derived.parents.iter().enumerate() {
            graph_poll(cx, parent_index)?;
            duplicate |= parent_index > 0 && &derived.parents[parent_index - 1] == parent;
            dangling |= match parent {
                ParentNodeId::Equation(id) => !equation_ids.contains(id),
                ParentNodeId::Variable(id) => !variable_ids.contains(id),
            };
            self_parent |= parent == &node;
            canonical_parents.push(parent.clone());
            edges.push((parent.clone(), node.clone()));
        }
        parents_by_child.insert(node.clone(), canonical_parents);
        if derived.parents.is_empty()
            || derived.parents.len() > MAX_CAUSAL_DERIVATION_PARENTS
            || duplicate
            || dangling
            || derived.differentiation_order > MAX_CAUSAL_DERIVATIVE_ORDER
            || self_parent
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::InvalidDerivedLineage,
                subject.clone(),
            ));
        }
        enforce_graph_finding_budget(findings)?;
    }

    let mut indegree = BTreeMap::<ParentNodeId, usize>::new();
    for (index, node) in derived_nodes.iter().enumerate() {
        graph_poll(cx, index)?;
        indegree.insert(node.clone(), 0);
    }
    let mut outgoing = BTreeMap::<ParentNodeId, Vec<ParentNodeId>>::new();
    for (index, (parent, child)) in edges.into_iter().enumerate() {
        graph_poll(cx, index)?;
        if derived_nodes.contains(&parent) {
            *indegree.entry(child.clone()).or_default() += 1;
            outgoing.entry(parent).or_default().push(child);
        }
    }
    let mut ready = BTreeSet::new();
    for (index, (node, degree)) in indegree.iter().enumerate() {
        graph_poll(cx, index)?;
        if *degree == 0 {
            ready.insert(node.clone());
        }
    }
    let mut visited = 0usize;
    let mut cumulative_orders = BTreeMap::<ParentNodeId, u32>::new();
    while let Some(node) = ready.pop_first() {
        graph_poll(cx, visited)?;
        visited += 1;
        let mut parent_order = 0u32;
        for (parent_index, parent) in parents_by_child
            .get(&node)
            .into_iter()
            .flatten()
            .enumerate()
        {
            graph_poll(cx, parent_index)?;
            if let Some(order) = cumulative_orders.get(parent) {
                parent_order = parent_order.max(*order);
            }
        }
        let cumulative_order = parent_order.saturating_add(u32::from(
            *local_orders
                .get(&node)
                .expect("derived node has a local differentiation order"),
        ));
        if cumulative_order > u32::from(MAX_CAUSAL_DERIVATIVE_ORDER) {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::InvalidDerivedLineage,
                subjects
                    .get(&node)
                    .expect("derived node has a diagnostic subject")
                    .clone(),
            ));
        }
        enforce_graph_finding_budget(findings)?;
        cumulative_orders.insert(node.clone(), cumulative_order);
        if let Some(children) = outgoing.get(&node) {
            for (child_index, child) in children.iter().enumerate() {
                graph_poll(cx, child_index)?;
                let degree = indegree
                    .get_mut(child)
                    .expect("derived child has an indegree row");
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(child.clone());
                }
            }
        }
    }
    if visited != derived_nodes.len() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidDerivedLineage,
            CausalGraphSubject::Graph,
        ));
    }
    Ok(())
}

fn graph_has_adjacent_duplicate<T: PartialEq>(
    values: &[T],
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    for (index, pair) in values.windows(2).enumerate() {
        graph_poll(cx, index)?;
        if pair[0] == pair[1] {
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn validate_conditions(
    draft: &CausalGraphDraft,
    equation_specs: &BTreeMap<EquationId, &EquationSpec>,
    equation_incidence_dependencies: &BTreeMap<EquationId, BTreeSet<VariableId>>,
    equation_incidences_always_available: &BTreeMap<EquationId, bool>,
    variable_specs: &BTreeMap<VariableId, &VariableSpec>,
    condition_ids: &BTreeSet<ActivationConditionRef>,
    condition_branches: &BTreeSet<(ActivationConditionRef, ActivationBranchRef)>,
    findings: &mut Vec<CausalGraphFinding>,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let mut used_dependencies = BTreeSet::new();
    let mut sourced_guard_equations = BTreeSet::new();
    for (condition_index, condition) in draft.conditions.iter().enumerate() {
        graph_poll(cx, condition_index)?;
        let subject = CausalGraphSubject::Condition(condition.condition.clone());
        let duplicate_branches = graph_has_adjacent_duplicate(&condition.branches, cx)?;
        let duplicate_dependencies = graph_has_adjacent_duplicate(&condition.dependencies, cx)?;
        let mut dependency_set = BTreeSet::new();
        let mut dangling = false;
        let mut audited_dependencies_valid = true;
        let mut dependencies_always_available = true;
        for (dependency_index, dependency) in condition.dependencies.iter().enumerate() {
            graph_poll(cx, dependency_index)?;
            let variable = variable_specs.get(dependency);
            dangling |= variable.is_none();
            audited_dependencies_valid &= variable.is_some_and(|variable| {
                matches!(
                    variable.solve_participation,
                    SolveParticipation::KnownInput | SolveParticipation::ConditionOnly
                )
            });
            dependencies_always_available &=
                variable.is_some_and(|variable| variable.activation == ActivationDomain::Always);
            dependency_set.insert(dependency.clone());
            used_dependencies.insert(dependency.clone());
        }
        let mut audit_mismatch = false;
        let source_valid = match &condition.source {
            ActivationConditionSource::GuardEquation {
                equation,
                obligation: _,
            } => {
                sourced_guard_equations.insert(equation.clone());
                if let Some(guard) = equation_specs.get(equation) {
                    let dependencies_match = match equation_incidence_dependencies.get(equation) {
                        Some(incidence_dependencies) => {
                            cancellable_set_eq(incidence_dependencies, &dependency_set, || {
                                graph_checkpoint(cx)
                            })?
                        }
                        None => false,
                    };
                    guard.role == EquationRole::Guard
                        && guard.solve_participation == EquationParticipation::ConditionOnly
                        && guard.activation == ActivationDomain::Always
                        && dependencies_always_available
                        && equation_incidences_always_available
                            .get(equation)
                            .copied()
                            .unwrap_or(false)
                        && dependencies_match
                } else {
                    false
                }
            }
            ActivationConditionSource::AuditedPredicate(escape) => {
                audit_mismatch = escape.source != escape.audited_source;
                audited_dependencies_valid && dependencies_always_available
            }
        };
        if audit_mismatch {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::EscapeAuditMismatch,
                subject.clone(),
            ));
        }
        if condition.branches.is_empty()
            || condition.dependencies.is_empty()
            || duplicate_branches
            || duplicate_dependencies
            || dangling
            || !source_valid
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::InvalidActivationCondition,
                subject,
            ));
        }
        enforce_graph_finding_budget(findings)?;
    }
    for (index, equation) in equation_specs.values().enumerate() {
        graph_poll(cx, index)?;
        if equation.role == EquationRole::Guard
            && (equation.solve_participation != EquationParticipation::ConditionOnly
                || !sourced_guard_equations.contains(&equation.id))
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::InvalidActivationCondition,
                CausalGraphSubject::Equation(equation.id.clone()),
            ));
        }
        enforce_graph_finding_budget(findings)?;
    }
    for (index, variable) in variable_specs.values().enumerate() {
        graph_poll(cx, index)?;
        if variable.solve_participation == SolveParticipation::ConditionOnly
            && !used_dependencies.contains(&variable.id)
        {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::InvalidActivationCondition,
                CausalGraphSubject::Variable(variable.id.clone()),
            ));
        }
        enforce_graph_finding_budget(findings)?;
    }
    let activations = draft
        .equations
        .iter()
        .map(|equation| &equation.activation)
        .chain(draft.variables.iter().map(|variable| &variable.activation))
        .chain(
            draft
                .incidences
                .iter()
                .map(|incidence| &incidence.activation),
        );
    let mut used_condition_ids = BTreeSet::new();
    let mut invalid_used_condition = false;
    let mut work = 0usize;
    for activation in activations {
        graph_poll(cx, work)?;
        work = work.saturating_add(1);
        for cube in activation_cubes(activation) {
            graph_poll(cx, work)?;
            work = work.saturating_add(1);
            for selection in &cube.selections {
                graph_poll(cx, work)?;
                work = work.saturating_add(1);
                invalid_used_condition |= !condition_ids.contains(&selection.condition)
                    || !condition_branches
                        .contains(&(selection.condition.clone(), selection.branch.clone()));
                used_condition_ids.insert(selection.condition.clone());
            }
        }
    }
    if invalid_used_condition {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidActivationCondition,
            CausalGraphSubject::Graph,
        ));
        enforce_graph_finding_budget(findings)?;
    }
    for (index, condition) in condition_ids.difference(&used_condition_ids).enumerate() {
        graph_poll(cx, index)?;
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::InvalidActivationCondition,
            CausalGraphSubject::Condition(condition.clone()),
        ));
        enforce_graph_finding_budget(findings)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn validate_incidence(
    incidence: &IncidenceSpec,
    equation_specs: &BTreeMap<EquationId, &EquationSpec>,
    variable_specs: &BTreeMap<VariableId, &VariableSpec>,
    condition_branches: &BTreeSet<(ActivationConditionRef, ActivationBranchRef)>,
    condition_domains: &BTreeMap<ActivationConditionRef, BTreeSet<ActivationBranchRef>>,
    view: &MachineView,
    findings: &mut Vec<CausalGraphFinding>,
    activation_proof_states: &mut usize,
    cx: &Cx<'_>,
) -> Result<(), CausalGraphRefusal> {
    let subject = incidence_subject(incidence);
    let equation = equation_specs.get(&incidence.equation).copied();
    let variable = variable_specs.get(&incidence.variable).copied();
    if equation.is_none() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownIncidenceEquation,
            subject.clone(),
        ));
    }
    if variable.is_none() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownIncidenceVariable,
            subject.clone(),
        ));
    }
    if incidence.derivative_order > MAX_CAUSAL_DERIVATIVE_ORDER {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::DerivativeOrderLimit,
            subject.clone(),
        ));
    }
    if !incidence.term.quantity.is_admitted() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnsupportedQuantityForm,
            subject.clone(),
        ));
    }
    if !view.clocks.contains(&incidence.term.clock) {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::UnknownClock,
            subject.clone(),
        ));
    }
    let (Some(equation), Some(variable)) = (equation, variable) else {
        return Ok(());
    };
    let equation_participation_valid = match equation.solve_participation {
        EquationParticipation::Matching => true,
        EquationParticipation::KnownClosure | EquationParticipation::ConditionOnly => !matches!(
            incidence.solve_participation,
            SolveParticipation::Unknown | SolveParticipation::ModeDependent
        ),
    };
    let condition_read = equation.solve_participation == EquationParticipation::ConditionOnly
        && incidence.solve_participation == SolveParticipation::ConditionOnly;
    let role_order_valid =
        variable.role != VariableRole::DiscreteMode || incidence.derivative_order == 0;
    let variable_participation_valid = role_order_valid
        && if condition_read {
            true
        } else if incidence.derivative_order == 0 {
            variable.solve_participation == SolveParticipation::ModeDependent
                || incidence.solve_participation == variable.solve_participation
        } else {
            match variable.role {
                VariableRole::State => matches!(
                    incidence.solve_participation,
                    SolveParticipation::Unknown | SolveParticipation::KnownInput
                ),
                VariableRole::Algebraic => {
                    incidence.solve_participation == SolveParticipation::Unknown
                }
                VariableRole::Source => {
                    incidence.solve_participation == SolveParticipation::KnownInput
                }
                VariableRole::Parameter => {
                    incidence.solve_participation == variable.solve_participation
                }
                VariableRole::PortEffort | VariableRole::PortFlow => {
                    variable.solve_participation == SolveParticipation::ModeDependent
                        || incidence.solve_participation == variable.solve_participation
                }
                VariableRole::DiscreteMode => false,
            }
        };
    let participation_valid = incidence.solve_participation != SolveParticipation::ModeDependent
        && equation_participation_valid
        && variable_participation_valid;
    if !participation_valid {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::DerivativeParticipationMismatch,
            subject.clone(),
        ));
    }
    if incidence.term != equation.residual {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::ResidualTermMismatch,
            subject.clone(),
        ));
    }
    if incidence.derivative_order <= MAX_CAUSAL_DERIVATIVE_ORDER {
        let mut derivative = [0i8; 6];
        derivative[TIME_DERIVATIVE_AXIS] = -(incidence.derivative_order as i8);
        let actual = variable
            .value
            .quantity
            .dims()
            .checked_plus(Dims(derivative))
            .and_then(|dims| dims.checked_plus(incidence.coefficient_dimensions));
        if actual != Some(incidence.term.quantity.dims()) {
            findings.push(CausalGraphFinding::new(
                CausalGraphRule::IncidenceUnitMismatch,
                subject.clone(),
            ));
        }
    }
    match &incidence.clock_relation {
        IncidenceClockRelation::SameClock => {
            if variable.value.clock != incidence.term.clock {
                findings.push(CausalGraphFinding::new(
                    CausalGraphRule::IncidenceClockMismatch,
                    subject.clone(),
                ));
            }
        }
        IncidenceClockRelation::AuditedBridge {
            source,
            target,
            bridge,
            audited_bridge,
            ..
        } => {
            if source != &variable.value.clock
                || target != &incidence.term.clock
                || source == target
                || !view.clocks.contains(source)
                || !view.clocks.contains(target)
                || bridge != audited_bridge
            {
                findings.push(CausalGraphFinding::new(
                    CausalGraphRule::IncidenceClockMismatch,
                    subject.clone(),
                ));
            }
        }
    }
    let semantic_transform =
        variable.value.quantity.semantic_type() != incidence.term.quantity.semantic_type();
    let structural_transform = variable.value.shape != incidence.term.shape
        || variable.value.frame != incidence.term.frame
        || incidence.coefficient_dimensions != Dims::NONE;
    if (semantic_transform || structural_transform) && incidence.operator.is_none() {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::MissingOperatorSemantics,
            subject.clone(),
        ));
    }
    let equation_contains = activation_implies_union(
        &incidence.activation,
        core::iter::once(&equation.activation),
        condition_domains,
        activation_proof_states,
        cx,
    )?;
    let variable_contains = activation_implies_union(
        &incidence.activation,
        core::iter::once(&variable.activation),
        condition_domains,
        activation_proof_states,
        cx,
    )?;
    if !equation_contains
        || !variable_contains
        || !activation_is_valid(&incidence.activation, condition_branches, cx)?
    {
        findings.push(CausalGraphFinding::new(
            CausalGraphRule::ActivationMismatch,
            subject,
        ));
    }
    Ok(())
}

fn activation_cubes(activation: &ActivationDomain) -> &[ActivationCube] {
    match activation {
        ActivationDomain::Always => &[],
        ActivationDomain::Conditional { cubes } => cubes,
    }
}

fn cubes_overlap_bounded(
    left: &ActivationCube,
    right: &ActivationCube,
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    let mut left_index = 0usize;
    let mut right_index = 0usize;
    while left_index < left.selections.len() && right_index < right.selections.len() {
        activation_proof_step(work, cx)?;
        let left_selection = &left.selections[left_index];
        let right_selection = &right.selections[right_index];
        match left_selection.condition.cmp(&right_selection.condition) {
            core::cmp::Ordering::Less => left_index += 1,
            core::cmp::Ordering::Greater => right_index += 1,
            core::cmp::Ordering::Equal => {
                if left_selection.branch != right_selection.branch {
                    return Ok(false);
                }
                left_index += 1;
                right_index += 1;
            }
        }
    }
    Ok(true)
}

fn activation_is_valid(
    activation: &ActivationDomain,
    condition_branches: &BTreeSet<(ActivationConditionRef, ActivationBranchRef)>,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    match activation {
        ActivationDomain::Always => Ok(true),
        ActivationDomain::Conditional { cubes } => {
            if cubes.is_empty() || cubes.len() > MAX_CAUSAL_CUBES_PER_ACTIVATION {
                return Ok(false);
            }
            for (cube_index, cube) in cubes.iter().enumerate() {
                graph_poll(cx, cube_index)?;
                let duplicate = cube_index > 0
                    && cancellable_slice_eq(
                        &cubes[cube_index - 1].selections,
                        &cube.selections,
                        || graph_checkpoint(cx),
                    )?;
                if duplicate
                    || cube.selections.is_empty()
                    || cube.selections.len() > MAX_CAUSAL_SELECTIONS_PER_CUBE
                {
                    return Ok(false);
                }
                for (selection_index, selection) in cube.selections.iter().enumerate() {
                    graph_poll(cx, selection_index)?;
                    if (selection_index > 0
                        && cube.selections[selection_index - 1].condition == selection.condition)
                        || !condition_branches
                            .contains(&(selection.condition.clone(), selection.branch.clone()))
                    {
                        return Ok(false);
                    }
                }
            }
            Ok(true)
        }
    }
}

fn activation_proof_step(states: &mut usize, cx: &Cx<'_>) -> Result<(), CausalGraphRefusal> {
    *states = states.saturating_add(1);
    if *states > MAX_CAUSAL_ACTIVATION_PROOF_STATES {
        return Err(graph_refusal(
            vec![CausalGraphFinding::new(
                CausalGraphRule::ResourceLimit,
                CausalGraphSubject::Graph,
            )],
            None,
        ));
    }
    graph_poll(cx, *states)
}

fn activation_domains_overlap_bounded(
    left: &ActivationDomain,
    right: &ActivationDomain,
    states: &mut usize,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    // Charge the comparison itself so Always and malformed empty-DNF cases
    // cannot bypass the global proof-work/cancellation envelope in a hostile
    // participation-conflict draft.
    activation_proof_step(states, cx)?;
    match (left, right) {
        (ActivationDomain::Always, _) | (_, ActivationDomain::Always) => Ok(true),
        (
            ActivationDomain::Conditional { cubes: left_cubes },
            ActivationDomain::Conditional { cubes: right_cubes },
        ) => {
            for left_cube in left_cubes {
                for right_cube in right_cubes {
                    activation_proof_step(states, cx)?;
                    if cubes_overlap_bounded(left_cube, right_cube, states, cx)? {
                        return Ok(true);
                    }
                }
            }
            Ok(false)
        }
    }
}

#[derive(Debug)]
struct ActivationProofDomain {
    condition_indices: BTreeMap<ActivationConditionRef, usize>,
    branches: Vec<Vec<ActivationBranchRef>>,
}

fn activation_proof_domain(
    domains: &BTreeMap<ActivationConditionRef, BTreeSet<ActivationBranchRef>>,
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<ActivationProofDomain, CausalGraphRefusal> {
    let mut condition_indices = BTreeMap::new();
    let mut branch_rows = Vec::with_capacity(domains.len());
    for (index, (condition, branches)) in domains.iter().enumerate() {
        activation_proof_step(work, cx)?;
        condition_indices.insert(condition.clone(), index);
        let mut branch_row = Vec::with_capacity(branches.len());
        for branch in branches {
            activation_proof_step(work, cx)?;
            branch_row.push(branch.clone());
        }
        branch_rows.push(branch_row);
    }
    Ok(ActivationProofDomain {
        condition_indices,
        branches: branch_rows,
    })
}

fn activation_branch_index(
    branches: &[ActivationBranchRef],
    target: &ActivationBranchRef,
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<Option<usize>, CausalGraphRefusal> {
    let mut lower = 0usize;
    let mut upper = branches.len();
    while lower < upper {
        activation_proof_step(work, cx)?;
        let middle = lower + (upper - lower) / 2;
        match branches[middle].cmp(target) {
            core::cmp::Ordering::Less => lower = middle + 1,
            core::cmp::Ordering::Greater => upper = middle,
            core::cmp::Ordering::Equal => return Ok(Some(middle)),
        }
    }
    Ok(None)
}

fn compile_activation_cube(
    cube: &ActivationCube,
    domain: &ActivationProofDomain,
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<Option<Vec<(usize, usize)>>, CausalGraphRefusal> {
    if cube.selections.is_empty() {
        return Ok(None);
    }
    let mut compiled = Vec::with_capacity(cube.selections.len());
    for selection in &cube.selections {
        activation_proof_step(work, cx)?;
        let Some(condition_index) = domain.condition_indices.get(&selection.condition).copied()
        else {
            return Ok(None);
        };
        let Some(branch_index) = activation_branch_index(
            &domain.branches[condition_index],
            &selection.branch,
            work,
            cx,
        )?
        else {
            return Ok(None);
        };
        compiled.push((condition_index, branch_index));
    }
    Ok(Some(compiled))
}

fn compiled_cube_has_uncovered_assignment(
    antecedent: &[(usize, usize)],
    consequents: &[Vec<(usize, usize)>],
    branch_counts: &[usize],
    work: &mut usize,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    activation_proof_step(work, cx)?;
    let mut assignment = vec![None; branch_counts.len()];
    for &(condition, branch) in antecedent {
        activation_proof_step(work, cx)?;
        let Some(slot) = assignment.get_mut(condition) else {
            return Ok(true);
        };
        if slot.is_some_and(|assigned| assigned != branch) {
            return Ok(true);
        }
        *slot = Some(branch);
    }

    // Each frame stores the condition introduced at one DFS depth and the
    // next branch to explore. The current assignment and stack are both
    // O(number of declared conditions); no search-state clone is retained.
    let mut frames = Vec::<(usize, usize)>::new();
    'search: loop {
        activation_proof_step(work, cx)?;
        let mut next_condition: Option<usize> = None;
        let mut any_compatible = false;
        let mut covered = false;
        for cube in consequents {
            activation_proof_step(work, cx)?;
            let mut compatible = true;
            let mut fully_satisfied = true;
            let mut first_missing = None;
            for &(condition, branch) in cube {
                activation_proof_step(work, cx)?;
                match assignment[condition] {
                    Some(assigned) if assigned != branch => {
                        compatible = false;
                        break;
                    }
                    Some(_) => {}
                    None => {
                        fully_satisfied = false;
                        first_missing.get_or_insert(condition);
                    }
                }
            }
            if !compatible {
                continue;
            }
            any_compatible = true;
            if fully_satisfied {
                covered = true;
                break;
            }
            if let Some(condition) = first_missing
                && next_condition.is_none_or(|current| condition < current)
            {
                next_condition = Some(condition);
            }
        }
        if !covered && !any_compatible {
            return Ok(true);
        }
        if !covered {
            let Some(condition) = next_condition else {
                return Ok(true);
            };
            let Some(branch_count) = branch_counts.get(condition).copied() else {
                return Ok(true);
            };
            if branch_count == 0 || assignment[condition].is_some() {
                return Ok(true);
            }
            activation_proof_step(work, cx)?;
            assignment[condition] = Some(0);
            frames.push((condition, 1));
            continue;
        }

        loop {
            let Some((condition, next_branch)) = frames.last_mut() else {
                return Ok(false);
            };
            activation_proof_step(work, cx)?;
            if *next_branch < branch_counts[*condition] {
                assignment[*condition] = Some(*next_branch);
                *next_branch += 1;
                continue 'search;
            }
            assignment[*condition] = None;
            frames.pop();
        }
    }
}

fn activation_implies_union<'a>(
    antecedent: &ActivationDomain,
    consequents: impl IntoIterator<Item = &'a ActivationDomain>,
    domains: &BTreeMap<ActivationConditionRef, BTreeSet<ActivationBranchRef>>,
    states: &mut usize,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    let mut consequent_activations = Vec::new();
    for activation in consequents {
        activation_proof_step(states, cx)?;
        if matches!(activation, ActivationDomain::Always) {
            return Ok(true);
        }
        consequent_activations.push(activation);
    }
    let domain = activation_proof_domain(domains, states, cx)?;
    let mut branch_counts = Vec::with_capacity(domain.branches.len());
    for branches in &domain.branches {
        activation_proof_step(states, cx)?;
        branch_counts.push(branches.len());
    }
    let mut consequent_cubes = Vec::new();
    for activation in consequent_activations {
        for cube in activation_cubes(activation) {
            activation_proof_step(states, cx)?;
            if let Some(compiled) = compile_activation_cube(cube, &domain, states, cx)? {
                consequent_cubes.push(compiled);
            }
        }
    }
    match antecedent {
        ActivationDomain::Always => compiled_cube_has_uncovered_assignment(
            &[],
            &consequent_cubes,
            &branch_counts,
            states,
            cx,
        )
        .map(|uncovered| !uncovered),
        ActivationDomain::Conditional { cubes } => {
            for cube in cubes {
                let Some(compiled) = compile_activation_cube(cube, &domain, states, cx)? else {
                    return Ok(false);
                };
                if compiled_cube_has_uncovered_assignment(
                    &compiled,
                    &consequent_cubes,
                    &branch_counts,
                    states,
                    cx,
                )? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
    }
}

fn assignment_is_valid_cancellable(
    assignment: &[ConditionBranchSelection],
    domain: &BTreeMap<ActivationConditionRef, BTreeSet<ActivationBranchRef>>,
    cx: &Cx<'_>,
) -> Result<bool, CausalReceiptRefusal> {
    if assignment.is_empty() || assignment.len() != domain.len() {
        return Ok(false);
    }
    for (index, selection) in assignment.iter().enumerate() {
        receipt_poll(cx, index)?;
        if (index > 0 && assignment[index - 1].condition == selection.condition)
            || !domain
                .get(&selection.condition)
                .is_some_and(|branches| branches.contains(&selection.branch))
        {
            return Ok(false);
        }
    }
    Ok(true)
}

fn activation_active_in_assignment_cancellable(
    activation: &ActivationDomain,
    assignment: &[ConditionBranchSelection],
    cx: &Cx<'_>,
) -> Result<bool, CausalReceiptRefusal> {
    match activation {
        ActivationDomain::Always => Ok(true),
        ActivationDomain::Conditional { cubes } => {
            for (cube_index, cube) in cubes.iter().enumerate() {
                receipt_poll(cx, cube_index)?;
                let mut active = true;
                for (selection_index, selection) in cube.selections.iter().enumerate() {
                    receipt_poll(cx, selection_index)?;
                    if assignment.binary_search(selection).is_err() {
                        active = false;
                        break;
                    }
                }
                if active {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}

fn receipt_row_is_active(
    activation: &ActivationDomain,
    assignment: Option<&[ConditionBranchSelection]>,
    cx: &Cx<'_>,
) -> Result<bool, CausalReceiptRefusal> {
    match assignment {
        Some(assignment) => activation_active_in_assignment_cancellable(activation, assignment, cx),
        None => Ok(matches!(activation, ActivationDomain::Always)),
    }
}

fn incidence_subject(incidence: &IncidenceSpec) -> CausalGraphSubject {
    CausalGraphSubject::Incidence {
        incidence: incidence.id.clone(),
        equation: incidence.equation.clone(),
        variable: incidence.variable.clone(),
        derivative_order: incidence.derivative_order,
    }
}

fn causal_structure_identity(
    draft: &CausalGraphDraft,
    machine: &AdmittedMachineGraph,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalStructureIdV1>, CausalIdentityError> {
    let mut machine_receipt_row = Vec::with_capacity(IDENTITY_RECEIPT_ADJUDICATION_BYTES);
    push_identity_receipt_adjudication(&mut machine_receipt_row, machine.identity_receipt());
    // Typed children retain recursive role/schema binding; their sibling
    // receipt rows keep digest collisions adjudicable across independently
    // admitted graphs without materializing the potentially huge child
    // canonical preimages themselves.
    let encoder =
        CanonicalEncoder::<CausalStructureIdV1, _>::new(CAUSAL_GRAPH_IDENTITY_LIMITS, || {
            cx.checkpoint().is_err()
        })?
        .u64(
            Field::new(0, "causal-structure-schema-version"),
            u64::from(CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1),
        )?
        .child(Field::new(1, "machine-graph-id"), machine.identity())?
        .bytes(
            Field::new(2, "machine-graph-receipt-adjudication"),
            &machine_receipt_row,
        )?
        .variant(
            Field::new(3, "unit-convention"),
            u32::from(unit_convention_tag(draft.units)),
            &[],
        )?
        .bytes(
            Field::new(4, "extraction-scope"),
            &extraction_scope_row(&draft.scope),
        )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(5, "equations"),
        &draft.equations,
        cx,
        equation_structure_row_len_cancellable,
        equation_structure_row,
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(6, "variables"),
        &draft.variables,
        cx,
        variable_structure_row_len_cancellable,
        variable_structure_row,
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(7, "activation-conditions"),
        &draft.conditions,
        cx,
        condition_canonical_len_cancellable,
        condition_row,
    )?;
    let encoder = encoder.ordered_children(
        Field::new(8, "incidences"),
        u64::try_from(draft.incidences.len()).map_err(|_| CanonicalError::LengthOverflow)?,
        draft
            .incidences
            .iter()
            .map(|incidence| incidence.id.identity()),
    )?;
    stream_identity_rows(
        encoder,
        Field::new(9, "incidence-receipt-adjudications"),
        &draft.incidences,
        cx,
        |_, _| Ok(IDENTITY_RECEIPT_ADJUDICATION_BYTES),
        |incidence, _, row| {
            push_identity_receipt_adjudication(row, incidence.id.identity_receipt());
            Ok(())
        },
    )?
    .finish()
    .map_err(Into::into)
}

fn causal_artifact_identity(
    draft: &CausalGraphDraft,
    structure: IdentityReceipt<CausalStructureIdV1>,
    behavior: Option<IdentityReceipt<MachineBehaviorIdV1>>,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalGraphArtifactIdV1>, CausalIdentityError> {
    let extraction = extraction_context_row(&draft.extraction);
    let behavior_row = behavior.map(machine_behavior_identity_row);
    let mut structure_row = Vec::with_capacity(IDENTITY_RECEIPT_ADJUDICATION_BYTES);
    push_identity_receipt_adjudication(&mut structure_row, structure);
    let encoder =
        CanonicalEncoder::<CausalGraphArtifactIdV1, _>::new(CAUSAL_GRAPH_IDENTITY_LIMITS, || {
            cx.checkpoint().is_err()
        })?
        .u64(
            Field::new(0, "causal-graph-artifact-schema-version"),
            u64::from(CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1),
        )?
        .child(Field::new(1, "causal-structure-id"), structure.id())?
        .bytes(
            Field::new(2, "causal-structure-receipt-adjudication"),
            &structure_row,
        )?
        .optional_bytes(
            Field::new(3, "machine-behavior-receipt-adjudication"),
            behavior_row.as_deref(),
        )?
        .bytes(Field::new(4, "extraction-context"), &extraction)?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(5, "equation-lineage"),
        &draft.equations,
        cx,
        |equation, cx| node_artifact_row_equation_len_cancellable(&equation.lineage, cx),
        |equation, cx, row| node_artifact_row_equation(&equation.id, &equation.lineage, cx, row),
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(6, "variable-lineage"),
        &draft.variables,
        cx,
        node_artifact_row_variable_len_cancellable,
        node_artifact_row_variable,
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(7, "activation-condition-provenance"),
        &draft.conditions,
        cx,
        |condition, _| condition_artifact_row_len(condition),
        |condition, _, row| {
            condition_artifact_row(condition, row);
            Ok(())
        },
    )?;
    stream_identity_rows(
        encoder,
        Field::new(8, "incidence-provenance"),
        &draft.incidences,
        cx,
        |incidence, _| incidence_artifact_row_len(incidence),
        |incidence, _, row| {
            incidence_artifact_row(incidence, row);
            Ok(())
        },
    )?
    .finish()
    .map_err(Into::into)
}

fn identity_materialization_checkpoint(
    cx: &Cx<'_>,
    absorbed_bytes: usize,
) -> Result<(), CanonicalError> {
    cx.checkpoint().map_err(|_| CanonicalError::Cancelled {
        absorbed_bytes: u64::try_from(absorbed_bytes).unwrap_or(u64::MAX),
    })
}

fn identity_materialization_poll(
    cx: &Cx<'_>,
    index: usize,
    absorbed_bytes: usize,
) -> Result<(), CanonicalError> {
    if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
        identity_materialization_checkpoint(cx, absorbed_bytes)?;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CausalIdentityRowPlan {
    row_index: u64,
    native_bytes: usize,
    canonical_bytes: u64,
}

impl CausalIdentityRowPlan {
    fn new(row_index: usize, native_bytes: usize) -> Result<Self, CanonicalError> {
        Ok(Self {
            row_index: u64::try_from(row_index).map_err(|_| CanonicalError::LengthOverflow)?,
            native_bytes,
            canonical_bytes: u64::try_from(native_bytes)
                .map_err(|_| CanonicalError::LengthOverflow)?,
        })
    }
}

fn checked_row_len_add(left: usize, right: usize) -> Result<usize, CanonicalError> {
    left.checked_add(right)
        .ok_or(CanonicalError::LengthOverflow)
}

fn checked_row_len_sum(parts: impl IntoIterator<Item = usize>) -> Result<usize, CanonicalError> {
    parts.into_iter().try_fold(0usize, checked_row_len_add)
}

fn push_row_count<O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    count: usize,
) -> Result<(), CanonicalError> {
    let count = u64::try_from(count).map_err(|_| CanonicalError::LengthOverflow)?;
    out.extend_from_slice(&count.to_le_bytes());
    Ok(())
}

fn reserve_identity_row(
    row: &mut Vec<u8>,
    plan: CausalIdentityRowPlan,
) -> Result<(), CausalIdentityAllocationRefusal> {
    row.try_reserve_exact(plan.native_bytes)
        .map_err(|_| CausalIdentityAllocationRefusal {
            requested_bytes: plan.canonical_bytes,
        })
}

/// Opaque producer surface for one already-admitted ordered-byte row.
///
/// The backing `Vec` is reserved fallibly before construction and never
/// exposed to the producer. Every attempted append is checked against the
/// logical row plan before mutation. An over-plan append records one stable
/// typed mismatch and becomes a no-op, so neither the row length nor its
/// allocation can grow because of the refused write.
struct CausalIdentityRowWriter {
    bytes: Vec<u8>,
    plan: CausalIdentityRowPlan,
    admitted_capacity: usize,
    refusal: Option<CanonicalError>,
}

impl CausalIdentityRowWriter {
    fn from_reserved(bytes: Vec<u8>, plan: CausalIdentityRowPlan) -> Result<Self, CanonicalError> {
        if !bytes.is_empty() {
            return Err(CanonicalError::InvalidSchemaDescriptor(
                "ordered row reservation returned nonempty storage",
            ));
        }
        if bytes.capacity() < plan.native_bytes {
            return Err(CanonicalError::InvalidSchemaDescriptor(
                "ordered row reservation returned less than the admitted plan",
            ));
        }
        let admitted_capacity = bytes.capacity();
        Ok(Self {
            bytes,
            plan,
            admitted_capacity,
            refusal: None,
        })
    }

    fn record_append(&mut self, additional_bytes: usize) -> bool {
        if self.refusal.is_some() {
            return false;
        }
        let Some(observed_native) = self.bytes.len().checked_add(additional_bytes) else {
            self.refusal = Some(CanonicalError::LengthOverflow);
            return false;
        };
        let Ok(observed) = u64::try_from(observed_native) else {
            self.refusal = Some(CanonicalError::LengthOverflow);
            return false;
        };
        if observed > self.plan.canonical_bytes {
            self.refusal = Some(CanonicalError::DeclaredLengthMismatch {
                declared: self.plan.canonical_bytes,
                observed,
            });
            return false;
        }
        true
    }

    fn refusal(&self) -> Option<CanonicalError> {
        self.refusal
    }

    fn into_bytes(self) -> Result<Vec<u8>, CanonicalError> {
        if let Some(refusal) = self.refusal {
            return Err(refusal);
        }
        debug_assert_eq!(
            self.bytes.capacity(),
            self.admitted_capacity,
            "cap-enforced row writes cannot reallocate admitted storage"
        );
        Ok(self.bytes)
    }

    #[cfg(test)]
    fn capacity(&self) -> usize {
        self.bytes.capacity()
    }
}

impl CausalCanonicalBytes for CausalIdentityRowWriter {
    fn len(&self) -> usize {
        self.bytes.len()
    }

    fn push(&mut self, byte: u8) {
        if self.record_append(1) {
            self.bytes.push(byte);
        }
    }

    fn extend_from_slice(&mut self, bytes: &[u8]) {
        if self.record_append(bytes.len()) {
            self.bytes.extend_from_slice(bytes);
        }
    }
}

fn stream_identity_rows<I, C, T, L, W>(
    encoder: CanonicalEncoder<I, C>,
    field: Field,
    values: &[T],
    cx: &Cx<'_>,
    row_len: L,
    write_row: W,
) -> Result<CanonicalEncoder<I, C>, CausalIdentityError>
where
    I: StrongIdentity,
    C: CancellationProbe,
    L: FnMut(&T, &Cx<'_>) -> Result<usize, CanonicalError>,
    W: FnMut(&T, &Cx<'_>, &mut CausalIdentityRowWriter) -> Result<(), CanonicalError>,
{
    stream_identity_rows_with_reserve(
        encoder,
        field,
        values,
        cx,
        row_len,
        write_row,
        reserve_identity_row,
    )
}

fn stream_identity_rows_with_reserve<I, C, T, L, W, R>(
    encoder: CanonicalEncoder<I, C>,
    field: Field,
    values: &[T],
    cx: &Cx<'_>,
    mut row_len: L,
    mut write_row: W,
    mut reserve_row: R,
) -> Result<CanonicalEncoder<I, C>, CausalIdentityError>
where
    I: StrongIdentity,
    C: CancellationProbe,
    L: FnMut(&T, &Cx<'_>) -> Result<usize, CanonicalError>,
    W: FnMut(&T, &Cx<'_>, &mut CausalIdentityRowWriter) -> Result<(), CanonicalError>,
    R: FnMut(&mut Vec<u8>, CausalIdentityRowPlan) -> Result<(), CausalIdentityAllocationRefusal>,
{
    let declared_count = u64::try_from(values.len()).map_err(|_| CanonicalError::LengthOverflow)?;
    // Length planning is deliberately independent from row production. The
    // encoder admits each exact declaration before the producer is allowed to
    // allocate or serialize that row.
    let pending_plan = core::cell::Cell::new(None);
    let row_lengths = values.iter().enumerate().map(|(row_index, value)| {
        row_len(value, cx)
            .and_then(|native_bytes| CausalIdentityRowPlan::new(row_index, native_bytes))
            .map(|plan| {
                pending_plan.set(Some(plan));
                plan.canonical_bytes
            })
            .map_err(CausalIdentityRowError::from)
    });
    encoder
        .ordered_bytes_stream(field, declared_count, row_lengths, |row_index, mut sink| {
            let plan = pending_plan
                .take()
                .ok_or(CausalIdentityRowError::Canonical(
                    CanonicalError::InvalidSchemaDescriptor(
                        "ordered row adapter lost its admitted length plan",
                    ),
                ))?;
            if plan.row_index != row_index {
                return Err(CausalIdentityRowError::Canonical(
                    CanonicalError::InvalidSchemaDescriptor(
                        "ordered row adapter plan index disagrees with producer index",
                    ),
                ));
            }
            let index = usize::try_from(row_index)
                .map_err(|_| CausalIdentityRowError::Canonical(CanonicalError::LengthOverflow))?;
            let value = values.get(index).ok_or(CausalIdentityRowError::Canonical(
                CanonicalError::InvalidSchemaDescriptor(
                    "ordered row adapter received an out-of-range row index",
                ),
            ))?;
            let mut bytes = Vec::new();
            reserve_row(&mut bytes, plan).map_err(CausalIdentityRowError::Allocation)?;
            let mut writer = CausalIdentityRowWriter::from_reserved(bytes, plan)
                .map_err(CausalIdentityRowError::Canonical)?;
            let producer_result = write_row(value, cx, &mut writer);
            if let Some(refusal) = writer.refusal() {
                return Err(CausalIdentityRowError::Canonical(refusal));
            }
            producer_result.map_err(CausalIdentityRowError::Canonical)?;
            let bytes = writer
                .into_bytes()
                .map_err(CausalIdentityRowError::Canonical)?;
            sink.write(&bytes)
                .map_err(CausalIdentityRowError::Canonical)
        })
        .map_err(Into::into)
}

fn terminal_quantity_canonical_len(quantity: TerminalQuantitySpec) -> usize {
    match quantity {
        TerminalQuantitySpec::Dimensional(_) => 1 + DIMENSION_COUNT,
        TerminalQuantitySpec::Semantic(semantic) => {
            let kind_bytes = match semantic.kind() {
                super::QuantityKind::Angle(_)
                | super::QuantityKind::AngularVelocity(_)
                | super::QuantityKind::Composition(_) => 2,
                super::QuantityKind::Strain { .. } => 3,
                super::QuantityKind::AbsoluteTemperature
                | super::QuantityKind::TemperatureDifference
                | super::QuantityKind::Torque
                | super::QuantityKind::Energy
                | super::QuantityKind::Pressure
                | super::QuantityKind::Stress
                | super::QuantityKind::Mass
                | super::QuantityKind::Amount
                | super::QuantityKind::MolarMass
                | super::QuantityKind::MassConcentration
                | super::QuantityKind::AmountConcentration
                | super::QuantityKind::Entropy
                | super::QuantityKind::HeatCapacity
                | super::QuantityKind::AcousticPressure
                | super::QuantityKind::AcousticPower => 1,
            };
            1 + kind_bytes + 1 + DIMENSION_COUNT
        }
    }
}

const fn terminal_shape_canonical_len(shape: TerminalShape) -> usize {
    match shape {
        TerminalShape::Scalar => 1,
        TerminalShape::Vector { .. } | TerminalShape::FieldTrace { .. } => 1 + 8,
        TerminalShape::Tensor { .. } => 1 + 16,
    }
}

fn signal_canonical_len(signal: &SignalContract) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        terminal_quantity_canonical_len(signal.quantity),
        terminal_shape_canonical_len(signal.shape),
        32,
        checked_row_len_add(8, signal.frame.canonical_key().len())?,
        1,
    ])
}

fn optional_causal_ref_canonical_len(namespace: Option<&str>) -> Result<usize, CanonicalError> {
    checked_row_len_add(
        1,
        match namespace {
            Some(namespace) => causal_ref_canonical_len(namespace)?,
            None => 0,
        },
    )
}

fn equation_structure_row_len_cancellable(
    equation: &EquationSpec,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        IDENTITY_RECEIPT_ADJUDICATION_BYTES,
        CAUSAL_OWNER_CANONICAL_BYTES,
        supports_canonical_len_cancellable(&equation.supports, cx)?,
        signal_canonical_len(&equation.residual)?,
        2,
        activation_canonical_len_cancellable(&equation.activation, cx)?,
    ])
}

fn variable_structure_row_len_cancellable(
    variable: &VariableSpec,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        IDENTITY_RECEIPT_ADJUDICATION_BYTES,
        CAUSAL_OWNER_CANONICAL_BYTES,
        supports_canonical_len_cancellable(&variable.supports, cx)?,
        signal_canonical_len(&variable.value)?,
        2,
        optional_causal_ref_canonical_len(
            variable
                .port_schema_crosswalk
                .as_ref()
                .map(|binding| binding.projection.namespace()),
        )?,
        activation_canonical_len_cancellable(&variable.activation, cx)?,
    ])
}

fn machine_behavior_identity_row(receipt: IdentityReceipt<MachineBehaviorIdV1>) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_len_prefixed(
        &mut out,
        b"org.frankensim.fs-ir.machine.behavior.v1/problem-semantic",
    );
    push_identity_receipt_adjudication(&mut out, receipt);
    out
}

fn extraction_scope_row(scope: &CausalGraphScope) -> Vec<u8> {
    let mut out = Vec::with_capacity(160);
    match scope {
        CausalGraphScope::CompleteMachineModel => out.push(1),
        CausalGraphScope::Partial { boundary } => {
            out.push(2);
            boundary.append_canonical(&mut out);
        }
    }
    out
}

fn extraction_context_row(context: &CausalExtractionContext) -> Vec<u8> {
    let mut out = Vec::with_capacity(640);
    context.extractor.append_canonical(&mut out);
    context.coverage.append_canonical(&mut out);
    match &context.evidence {
        CausalExtractionEvidence::Unverified => out.push(1),
        CausalExtractionEvidence::CheckerReferenced(checker) => {
            out.push(2);
            checker.append_canonical(&mut out);
        }
    }
    context.budget.append_canonical(&mut out);
    context.capabilities.append_canonical(&mut out);
    push_seed_policy(&mut out, context.seed_policy);
    out.push(determinism_tag(context.determinism));
    out
}

fn equation_structure_row(
    equation: &EquationSpec,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    push_identity_receipt_adjudication(out, equation.id.identity_receipt());
    push_owner(out, &equation.owner);
    push_supports_cancellable(out, &equation.supports, cx)?;
    push_signal(out, &equation.residual)?;
    out.push(equation_role_tag(equation.role));
    out.push(equation_participation_tag(equation.solve_participation));
    identity_materialization_checkpoint(cx, out.len())?;
    push_activation_cancellable(out, &equation.activation, cx)?;
    identity_materialization_checkpoint(cx, out.len())
}

fn variable_structure_row(
    variable: &VariableSpec,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    push_identity_receipt_adjudication(out, variable.id.identity_receipt());
    push_owner(out, &variable.owner);
    push_supports_cancellable(out, &variable.supports, cx)?;
    push_signal(out, &variable.value)?;
    out.push(variable_role_tag(variable.role));
    out.push(solve_participation_tag(variable.solve_participation));
    push_optional_ref(
        out,
        variable
            .port_schema_crosswalk
            .as_ref()
            .map(|binding| &binding.projection),
    );
    identity_materialization_checkpoint(cx, out.len())?;
    push_activation_cancellable(out, &variable.activation, cx)?;
    identity_materialization_checkpoint(cx, out.len())
}

fn condition_row(
    condition: &ActivationConditionSpec,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    condition.condition.append_canonical(out);
    match &condition.source {
        ActivationConditionSource::GuardEquation {
            equation,
            obligation,
        } => {
            out.push(1);
            push_identity_receipt_adjudication(out, equation.identity_receipt());
            obligation.append_canonical(out);
        }
        // `condition`, branches, and dependencies carry the normalized
        // predicate semantics. The exact opaque implementation and its audit
        // belong exclusively to `condition_artifact_row`, preserving the
        // structure/artifact identity split.
        ActivationConditionSource::AuditedPredicate(_) => out.push(2),
    }
    push_row_count(out, condition.branches.len())?;
    for (index, branch) in condition.branches.iter().enumerate() {
        identity_materialization_poll(cx, index, out.len())?;
        branch.append_canonical(out);
    }
    push_row_count(out, condition.dependencies.len())?;
    for (index, dependency) in condition.dependencies.iter().enumerate() {
        identity_materialization_poll(cx, index, out.len())?;
        push_identity_receipt_adjudication(out, dependency.identity_receipt());
    }
    identity_materialization_checkpoint(cx, out.len())
}

fn condition_artifact_row(condition: &ActivationConditionSpec, out: &mut CausalIdentityRowWriter) {
    debug_assert!(out.is_empty());
    condition.condition.append_canonical(out);
    match &condition.source {
        ActivationConditionSource::GuardEquation { .. } => out.push(1),
        ActivationConditionSource::AuditedPredicate(escape) => {
            out.push(2);
            escape.audit.append_canonical(out);
            escape.audited_source.append_canonical(out);
        }
    }
}

fn condition_artifact_row_len(
    condition: &ActivationConditionSpec,
) -> Result<usize, CanonicalError> {
    let source_bytes = match &condition.source {
        ActivationConditionSource::GuardEquation { .. } => 0,
        ActivationConditionSource::AuditedPredicate(escape) => checked_row_len_sum([
            causal_ref_canonical_len(escape.audit.namespace())?,
            causal_ref_canonical_len(escape.audited_source.namespace())?,
        ])?,
    };
    checked_row_len_sum([
        causal_ref_canonical_len(condition.condition.namespace())?,
        1,
        source_bytes,
    ])
}

fn incidence_meaning_row_cancellable(
    incidence: &IncidenceSpec,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CanonicalError> {
    incidence_meaning_row_parts_cancellable(
        &incidence.equation,
        &incidence.variable,
        incidence.derivative_order,
        incidence.solve_participation,
        incidence.coefficient_dimensions,
        &incidence.term,
        incidence.operator.as_ref(),
        &incidence.clock_relation,
        &incidence.activation,
        cx,
    )
}

#[allow(clippy::too_many_arguments)]
fn incidence_meaning_row_parts_cancellable(
    equation: &EquationId,
    variable: &VariableId,
    derivative_order: u16,
    solve_participation: SolveParticipation,
    coefficient_dimensions: Dims,
    term: &SignalContract,
    operator: Option<&IncidenceOperatorRef>,
    clock_relation: &IncidenceClockRelation,
    activation: &ActivationDomain,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CanonicalError> {
    let activation_bytes = activation_canonical_len_cancellable(activation, cx)?;
    let canonical_bytes =
        checked_row_len_add(INCIDENCE_FIXED_CANONICAL_CAPACITY, activation_bytes)?;
    let mut out = Vec::with_capacity(canonical_bytes);
    push_incidence_meaning_fixed(
        &mut out,
        equation,
        variable,
        derivative_order,
        solve_participation,
        coefficient_dimensions,
        term,
        operator,
        clock_relation,
    )?;
    debug_assert!(out.len() <= INCIDENCE_FIXED_CANONICAL_CAPACITY);
    identity_materialization_checkpoint(cx, out.len())?;
    push_activation_cancellable(&mut out, activation, cx)?;
    identity_materialization_checkpoint(cx, out.len())?;
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn push_incidence_meaning_fixed(
    out: &mut Vec<u8>,
    equation: &EquationId,
    variable: &VariableId,
    derivative_order: u16,
    solve_participation: SolveParticipation,
    coefficient_dimensions: Dims,
    term: &SignalContract,
    operator: Option<&IncidenceOperatorRef>,
    clock_relation: &IncidenceClockRelation,
) -> Result<(), CanonicalError> {
    push_identity_receipt_adjudication(out, equation.identity_receipt());
    push_identity_receipt_adjudication(out, variable.identity_receipt());
    out.extend_from_slice(&derivative_order.to_le_bytes());
    out.push(solve_participation_tag(solve_participation));
    push_dims(out, coefficient_dimensions);
    push_signal(out, term)?;
    push_optional_ref(out, operator);
    match clock_relation {
        IncidenceClockRelation::SameClock => out.push(1),
        IncidenceClockRelation::AuditedBridge {
            source,
            target,
            bridge,
            audit: _,
            audited_bridge: _,
        } => {
            out.push(2);
            out.extend_from_slice(source.identity().as_bytes());
            out.extend_from_slice(target.identity().as_bytes());
            bridge.append_canonical(out);
        }
    }
    Ok(())
}

fn incidence_receipt_from_meaning_cancellable(
    meaning: &[u8],
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<IncidenceEntityIdV1>, CanonicalError> {
    CanonicalEncoder::<IncidenceEntityIdV1, _>::new(INCIDENCE_IDENTITY_LIMITS, || {
        cx.checkpoint().is_err()
    })?
    .bytes(Field::new(0, "normalized-incidence-meaning"), meaning)?
    .finish()
}

fn incidence_identity_matches_cancellable(
    incidence: &IncidenceSpec,
    cx: &Cx<'_>,
) -> Result<bool, CausalGraphRefusal> {
    let recomputed = incidence_meaning_row_cancellable(incidence, cx)
        .and_then(|meaning| incidence_receipt_from_meaning_cancellable(&meaning, cx));
    match recomputed {
        Ok(recomputed) => Ok(identity_receipt_adjudication_eq(
            recomputed,
            incidence.id.identity_receipt(),
        )),
        Err(CanonicalError::Cancelled { .. }) => Err(cancelled_graph_refusal()),
        Err(_) => Ok(false),
    }
}

fn node_artifact_row_equation(
    id: &EquationId,
    lineage: &NodeLineage,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    let lineage_bytes = lineage.canonical_row_len_cancellable(cx)?;
    push_identity_receipt_adjudication(out, id.identity_receipt());
    push_row_count(out, lineage_bytes)?;
    lineage.append_canonical_row_cancellable(out, cx)?;
    identity_materialization_checkpoint(cx, out.len())
}

fn node_artifact_row_equation_len_cancellable(
    lineage: &NodeLineage,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        IDENTITY_RECEIPT_ADJUDICATION_BYTES,
        8,
        lineage.canonical_row_len_cancellable(cx)?,
    ])
}

fn node_artifact_row_variable(
    variable: &VariableSpec,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    let lineage_bytes = variable.lineage.canonical_row_len_cancellable(cx)?;
    push_identity_receipt_adjudication(out, variable.id.identity_receipt());
    push_row_count(out, lineage_bytes)?;
    variable.lineage.append_canonical_row_cancellable(out, cx)?;
    push_optional_ref(
        out,
        variable
            .port_schema_crosswalk
            .as_ref()
            .map(|binding| &binding.audit),
    );
    identity_materialization_checkpoint(cx, out.len())
}

fn node_artifact_row_variable_len_cancellable(
    variable: &VariableSpec,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        IDENTITY_RECEIPT_ADJUDICATION_BYTES,
        8,
        variable.lineage.canonical_row_len_cancellable(cx)?,
        optional_causal_ref_canonical_len(
            variable
                .port_schema_crosswalk
                .as_ref()
                .map(|binding| binding.audit.namespace()),
        )?,
    ])
}

fn incidence_artifact_row(incidence: &IncidenceSpec, out: &mut CausalIdentityRowWriter) {
    debug_assert!(out.is_empty());
    push_identity_receipt_adjudication(out, incidence.id.identity_receipt());
    match &incidence.clock_relation {
        IncidenceClockRelation::SameClock => out.push(1),
        IncidenceClockRelation::AuditedBridge {
            audit,
            audited_bridge,
            ..
        } => {
            out.push(2);
            audit.append_canonical(out);
            audited_bridge.append_canonical(out);
        }
    }
}

fn incidence_artifact_row_len(incidence: &IncidenceSpec) -> Result<usize, CanonicalError> {
    let relation_bytes = match &incidence.clock_relation {
        IncidenceClockRelation::SameClock => 1,
        IncidenceClockRelation::AuditedBridge {
            audit,
            audited_bridge,
            ..
        } => checked_row_len_sum([
            1,
            causal_ref_canonical_len(audit.namespace())?,
            causal_ref_canonical_len(audited_bridge.namespace())?,
        ])?,
    };
    checked_row_len_add(IDENTITY_RECEIPT_ADJUDICATION_BYTES, relation_bytes)
}

fn node_origin_canonical_len_cancellable(
    origin: &NodeOrigin,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let bytes = match origin {
        NodeOrigin::Machine(machine) => match machine {
            MachineNodeOrigin::Subsystem(_)
            | MachineNodeOrigin::Relation(_)
            | MachineNodeOrigin::Interface(_) => 2 + 32,
            MachineNodeOrigin::Element(_) | MachineNodeOrigin::PortCoordinate { .. } => 2 + 33,
        },
        NodeOrigin::AuditedEscapeHatch(escape) => checked_row_len_sum([
            1,
            causal_ref_canonical_len(escape.source.namespace())?,
            causal_ref_canonical_len(escape.audit.namespace())?,
            causal_ref_canonical_len(escape.audited_source.namespace())?,
        ])?,
        NodeOrigin::Derived(derived) => {
            let mut bytes = checked_row_len_add(1, 8)?;
            for (index, _) in derived.parents.iter().enumerate() {
                identity_materialization_poll(cx, index, bytes)?;
                bytes = checked_row_len_add(bytes, 1 + IDENTITY_RECEIPT_ADJUDICATION_BYTES)?;
            }
            checked_row_len_sum([
                bytes,
                causal_ref_canonical_len(derived.transformation.namespace())?,
                2,
                causal_ref_canonical_len(derived.obligation.namespace())?,
            ])?
        }
    };
    identity_materialization_checkpoint(cx, bytes)?;
    Ok(bytes)
}

fn push_node_origin_cancellable<O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    origin: &NodeOrigin,
    cx: &Cx<'_>,
) -> Result<(), CanonicalError> {
    match origin {
        NodeOrigin::Machine(machine) => {
            out.push(1);
            push_machine_origin(out, machine);
        }
        NodeOrigin::AuditedEscapeHatch(escape) => {
            out.push(2);
            escape.source.append_canonical(out);
            escape.audit.append_canonical(out);
            escape.audited_source.append_canonical(out);
        }
        NodeOrigin::Derived(derived) => {
            out.push(3);
            push_row_count(out, derived.parents.len())?;
            for (index, parent) in derived.parents.iter().enumerate() {
                identity_materialization_poll(cx, index, out.len())?;
                match parent {
                    ParentNodeId::Equation(id) => {
                        out.push(1);
                        push_identity_receipt_adjudication(out, id.identity_receipt());
                    }
                    ParentNodeId::Variable(id) => {
                        out.push(2);
                        push_identity_receipt_adjudication(out, id.identity_receipt());
                    }
                }
            }
            derived.transformation.append_canonical(out);
            out.extend_from_slice(&derived.differentiation_order.to_le_bytes());
            derived.obligation.append_canonical(out);
        }
    }
    identity_materialization_checkpoint(cx, out.len())?;
    Ok(())
}

fn push_machine_origin<O: CausalCanonicalBytes + ?Sized>(out: &mut O, origin: &MachineNodeOrigin) {
    match origin {
        MachineNodeOrigin::Subsystem(id) => {
            out.push(1);
            out.extend_from_slice(id.identity().as_bytes());
        }
        MachineNodeOrigin::Relation(id) => {
            out.push(2);
            out.extend_from_slice(id.identity().as_bytes());
        }
        MachineNodeOrigin::Element(element) => {
            out.push(3);
            out.extend_from_slice(&element.canonical_row());
        }
        MachineNodeOrigin::PortCoordinate { port, coordinate } => {
            out.push(4);
            out.extend_from_slice(port.identity().as_bytes());
            out.push(port_coordinate_tag(*coordinate));
        }
        MachineNodeOrigin::Interface(id) => {
            out.push(5);
            out.extend_from_slice(id.identity().as_bytes());
        }
    }
}

fn push_owner<O: CausalCanonicalBytes + ?Sized>(out: &mut O, owner: &CausalOwner) {
    match owner {
        CausalOwner::Subsystem(id) => {
            out.push(1);
            out.extend_from_slice(id.identity().as_bytes());
        }
        CausalOwner::Port(id) => {
            out.push(2);
            out.extend_from_slice(id.identity().as_bytes());
        }
        CausalOwner::Interface(id) => {
            out.push(3);
            out.extend_from_slice(id.identity().as_bytes());
        }
    }
}

fn push_supports_cancellable<O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    supports: &[CausalSupport],
    cx: &Cx<'_>,
) -> Result<(), CanonicalError> {
    push_row_count(out, supports.len())?;
    for (index, support) in supports.iter().enumerate() {
        identity_materialization_poll(cx, index, out.len())?;
        match support {
            CausalSupport::Lumped => out.push(1),
            CausalSupport::MachineElement(element) => {
                out.push(2);
                out.extend_from_slice(&element.canonical_row());
            }
            CausalSupport::External(reference) => {
                out.push(3);
                reference.append_canonical(out);
            }
        }
    }
    Ok(())
}

fn supports_canonical_len_cancellable(
    supports: &[CausalSupport],
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let mut bytes = 8usize;
    for (index, support) in supports.iter().enumerate() {
        identity_materialization_poll(cx, index, bytes)?;
        let support_bytes = match support {
            CausalSupport::Lumped => 1,
            CausalSupport::MachineElement(_) => 1 + 33,
            CausalSupport::External(reference) => {
                checked_row_len_add(1, causal_ref_canonical_len(reference.namespace())?)?
            }
        };
        bytes = checked_row_len_add(bytes, support_bytes)?;
    }
    identity_materialization_checkpoint(cx, bytes)?;
    Ok(bytes)
}

fn condition_canonical_len_cancellable(
    condition: &ActivationConditionSpec,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let mut bytes = causal_ref_canonical_len(condition.condition.namespace())?;
    let source_bytes = match &condition.source {
        ActivationConditionSource::GuardEquation { obligation, .. } => checked_row_len_sum([
            1,
            IDENTITY_RECEIPT_ADJUDICATION_BYTES,
            causal_ref_canonical_len(obligation.namespace())?,
        ])?,
        ActivationConditionSource::AuditedPredicate(_) => 1,
    };
    bytes = checked_row_len_add(bytes, source_bytes)?;
    bytes = checked_row_len_add(bytes, 8)?;
    for (index, branch) in condition.branches.iter().enumerate() {
        identity_materialization_poll(cx, index, bytes)?;
        bytes = checked_row_len_add(bytes, causal_ref_canonical_len(branch.namespace())?)?;
    }
    bytes = checked_row_len_add(bytes, 8)?;
    for (index, _) in condition.dependencies.iter().enumerate() {
        identity_materialization_poll(cx, index, bytes)?;
        bytes = checked_row_len_add(bytes, IDENTITY_RECEIPT_ADJUDICATION_BYTES)?;
    }
    identity_materialization_checkpoint(cx, bytes)?;
    Ok(bytes)
}

fn push_signal<O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    signal: &SignalContract,
) -> Result<(), CanonicalError> {
    match signal.quantity {
        TerminalQuantitySpec::Dimensional(dims) => {
            out.push(1);
            push_dims(out, dims);
        }
        TerminalQuantitySpec::Semantic(semantic_type) => {
            out.push(2);
            match semantic_type.kind() {
                super::QuantityKind::AbsoluteTemperature => out.push(1),
                super::QuantityKind::TemperatureDifference => out.push(2),
                super::QuantityKind::Angle(domain) => {
                    out.push(3);
                    out.push(super::angle_domain_tag(domain));
                }
                super::QuantityKind::AngularVelocity(domain) => {
                    out.push(4);
                    out.push(super::angle_domain_tag(domain));
                }
                super::QuantityKind::Torque => out.push(5),
                super::QuantityKind::Energy => out.push(6),
                super::QuantityKind::Pressure => out.push(7),
                super::QuantityKind::Stress => out.push(8),
                super::QuantityKind::Strain { basis, component } => {
                    out.push(9);
                    out.push(super::strain_basis_tag(basis));
                    out.push(super::strain_component_tag(component));
                }
                super::QuantityKind::Composition(basis) => {
                    out.push(10);
                    out.push(super::composition_basis_tag(basis));
                }
                super::QuantityKind::Mass => out.push(11),
                super::QuantityKind::Amount => out.push(12),
                super::QuantityKind::MolarMass => out.push(13),
                super::QuantityKind::MassConcentration => out.push(14),
                super::QuantityKind::AmountConcentration => out.push(15),
                super::QuantityKind::Entropy => out.push(16),
                super::QuantityKind::HeatCapacity => out.push(17),
                super::QuantityKind::AcousticPressure => out.push(18),
                super::QuantityKind::AcousticPower => out.push(19),
            }
            out.push(super::value_form_tag(semantic_type.form()));
            push_dims(out, semantic_type.expected_dims());
        }
    }
    match signal.shape {
        TerminalShape::Scalar => out.push(1),
        TerminalShape::Vector { components } => {
            out.push(2);
            out.extend_from_slice(&components.get().to_le_bytes());
        }
        TerminalShape::Tensor { rows, columns } => {
            out.push(3);
            out.extend_from_slice(&rows.get().to_le_bytes());
            out.extend_from_slice(&columns.get().to_le_bytes());
        }
        TerminalShape::FieldTrace { components } => {
            out.push(4);
            out.extend_from_slice(&components.get().to_le_bytes());
        }
    }
    out.extend_from_slice(signal.clock.identity().as_bytes());
    push_row_count(out, signal.frame.canonical_key().len())?;
    out.extend_from_slice(signal.frame.canonical_key().as_bytes());
    out.push(match signal.frame.orientation() {
        super::OrientationParity::Preserving => 1,
        super::OrientationParity::Reversing => 2,
    });
    Ok(())
}

fn push_dims<O: CausalCanonicalBytes + ?Sized>(out: &mut O, dims: Dims) {
    out.extend_from_slice(&dims.0.map(|exponent| exponent as u8));
}

fn push_activation_cancellable<O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    activation: &ActivationDomain,
    cx: &Cx<'_>,
) -> Result<(), CanonicalError> {
    match activation {
        ActivationDomain::Always => out.push(1),
        ActivationDomain::Conditional { cubes } => {
            out.push(2);
            push_row_count(out, cubes.len())?;
            for (cube_index, cube) in cubes.iter().enumerate() {
                identity_materialization_poll(cx, cube_index, out.len())?;
                push_row_count(out, cube.selections.len())?;
                for (selection_index, selection) in cube.selections.iter().enumerate() {
                    identity_materialization_poll(cx, selection_index, out.len())?;
                    selection.condition.append_canonical(out);
                    selection.branch.append_canonical(out);
                }
            }
        }
    }
    Ok(())
}

fn activation_canonical_len_cancellable(
    activation: &ActivationDomain,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let mut bytes = 1usize;
    if let ActivationDomain::Conditional { cubes } = activation {
        bytes = checked_row_len_add(bytes, 8)?;
        for (cube_index, cube) in cubes.iter().enumerate() {
            identity_materialization_poll(cx, cube_index, bytes)?;
            bytes = checked_row_len_add(bytes, 8)?;
            for (selection_index, selection) in cube.selections.iter().enumerate() {
                identity_materialization_poll(cx, selection_index, bytes)?;
                bytes = checked_row_len_sum([
                    bytes,
                    causal_ref_canonical_len(selection.condition.namespace())?,
                    causal_ref_canonical_len(selection.branch.namespace())?,
                ])?;
            }
        }
    }
    identity_materialization_checkpoint(cx, bytes)?;
    Ok(bytes)
}

fn assignment_canonical_len_cancellable(
    assignment: &[ConditionBranchSelection],
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let mut bytes = 8usize;
    for (index, selection) in assignment.iter().enumerate() {
        identity_materialization_poll(cx, index, bytes)?;
        bytes = checked_row_len_sum([
            bytes,
            causal_ref_canonical_len(selection.condition.namespace())?,
            causal_ref_canonical_len(selection.branch.namespace())?,
        ])?;
    }
    identity_materialization_checkpoint(cx, bytes)?;
    Ok(bytes)
}

fn causal_ref_canonical_len(namespace: &str) -> Result<usize, CanonicalError> {
    checked_row_len_sum([8, namespace.len(), 8, 32])
}

trait AppendCausalRef {
    fn append_to<O: CausalCanonicalBytes + ?Sized>(&self, out: &mut O);
}

macro_rules! append_ref_impl {
    ($($name:ident),+ $(,)?) => {
        $(
            impl AppendCausalRef for $name {
                fn append_to<O: CausalCanonicalBytes + ?Sized>(&self, out: &mut O) {
                    self.append_canonical(out);
                }
            }
        )+
    };
}

append_ref_impl!(
    PortSchemaCrosswalkRef,
    PortSchemaCrosswalkAuditRef,
    IncidenceOperatorRef,
    MaximumMatchingCertificateRef,
    ConditionalCoverageRef,
    CausalCheckpointRef,
);

fn push_optional_ref<T: AppendCausalRef, O: CausalCanonicalBytes + ?Sized>(
    out: &mut O,
    reference: Option<&T>,
) {
    match reference {
        None => out.push(0),
        Some(reference) => {
            out.push(1);
            reference.append_to(out);
        }
    }
}

fn push_len_prefixed<O: CausalCanonicalBytes + ?Sized>(out: &mut O, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

const fn port_coordinate_tag(coordinate: PortCoordinate) -> u8 {
    match coordinate {
        PortCoordinate::Effort => 1,
        PortCoordinate::Flow => 2,
    }
}

const fn equation_role_tag(role: EquationRole) -> u8 {
    match role {
        EquationRole::Balance => 1,
        EquationRole::Constitutive => 2,
        EquationRole::Constraint => 3,
        EquationRole::Source => 4,
        EquationRole::StateUpdate => 5,
        EquationRole::PortClosure => 6,
        EquationRole::Guard => 7,
    }
}

const fn equation_participation_tag(participation: EquationParticipation) -> u8 {
    match participation {
        EquationParticipation::Matching => 1,
        EquationParticipation::KnownClosure => 2,
        EquationParticipation::ConditionOnly => 3,
    }
}

const fn variable_role_tag(role: VariableRole) -> u8 {
    match role {
        VariableRole::Algebraic => 1,
        VariableRole::State => 2,
        VariableRole::Source => 3,
        VariableRole::Parameter => 4,
        VariableRole::PortEffort => 5,
        VariableRole::PortFlow => 6,
        VariableRole::DiscreteMode => 7,
    }
}

const fn solve_participation_tag(participation: SolveParticipation) -> u8 {
    match participation {
        SolveParticipation::Unknown => 1,
        SolveParticipation::KnownInput => 2,
        SolveParticipation::ConditionOnly => 3,
        SolveParticipation::ModeDependent => 4,
    }
}

const fn unit_convention_tag(units: CausalUnitConvention) -> u8 {
    match units {
        CausalUnitConvention::SiBaseDimensions => 1,
    }
}

const fn determinism_tag(determinism: CausalDeterminism) -> u8 {
    match determinism {
        CausalDeterminism::Deterministic => 1,
        CausalDeterminism::Relaxed => 2,
    }
}

fn push_seed_policy(out: &mut Vec<u8>, seed_policy: CausalSeedPolicy) {
    match seed_policy {
        CausalSeedPolicy::NoRandomness => out.push(1),
        CausalSeedPolicy::CounterBased { seed, stream } => {
            out.push(2);
            out.extend_from_slice(&seed.to_le_bytes());
            out.extend_from_slice(&stream.to_le_bytes());
        }
    }
}

/// Matching vertex for one base variable at one explicit derivative order.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivativeVariableKey {
    /// Base variable identity.
    pub variable: VariableId,
    /// Explicit time-derivative order.
    pub derivative_order: u16,
}

/// One selected structural equation-to-derivative-variable edge.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CausalMatchingPair {
    /// Exact incidence selected. This disambiguates branch/operator/clock
    /// variants with the same endpoint and derivative key.
    pub incidence: IncidenceId,
    /// Matched equation.
    pub equation: EquationId,
    /// Matched right-hand structural vertex.
    pub variable: DerivativeVariableKey,
}

/// Cardinality/deficiency classification. Under- and over-determined regions
/// may coexist and therefore have an explicit `Mixed` state.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DeterminationClass {
    /// No unmatched equations or unknown vertices in the declared witness.
    WellDetermined,
    /// Unmatched unknown vertices only.
    UnderDetermined,
    /// Unmatched equations only.
    OverDetermined,
    /// Both unmatched equations and unknown vertices.
    Mixed,
    /// Analyzer did not make a determination claim.
    Unknown,
    /// Both structural bipartition sides are empty in the exact analysis
    /// domain. This is a concrete domain result, not a vacuous full-rank or
    /// solvability claim; it is required for honest off/disengaged mode cells.
    /// Declared last to preserve the established tag/order of earlier states.
    EmptyProjection,
}

/// Generic structural-rank state, orthogonal to graph cardinality.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StructuralRankState {
    /// Matching reaches the smaller bipartition side.
    FullRelativeToMinSide,
    /// A checker-supported structural deficiency remains.
    Deficient,
    /// At least one bipartition side is empty, so min-side rank is vacuous
    /// rather than an informative full-rank claim.
    NotApplicable,
    /// No structural-rank claim.
    Unknown,
}

/// Whether one result is unconditional or condition/mode dependent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Conditionality {
    /// One result applies across the analyzer's declared domain.
    Unconditional,
    /// Separate branch/mode receipts are bound below.
    Conditional,
    /// Analyzer did not establish conditional coverage.
    Unknown,
}

/// Exact structural domain analyzed by one receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CausalReceiptDomain {
    /// Entire graph, admitted only when the graph has no activation conditions.
    UnconditionalGraph,
    /// One complete joint-mode assignment. Only rows active in this cell may
    /// participate in matching and complements.
    ModeCell {
        /// Exactly one selection for every declared graph condition.
        assignment: Vec<ConditionBranchSelection>,
    },
    /// Cross-mode summary with no union-graph matching witness. Exact
    /// mode-cell results live in `conditional_outcomes`.
    HybridSummary,
}

/// Typed reason why one or more receipt axes remain unknown.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalUnknownReason {
    /// Analysis was not attempted.
    NotAnalyzed,
    /// Request was cancelled and fully drained before publication.
    Cancelled,
    /// Structure class is outside the analyzer's admitted capability.
    UnsupportedStructure,
    /// Explicit time/memory/capability budget was exhausted.
    BudgetExhausted,
    /// Required source metadata was unavailable.
    IncompleteMetadata,
    /// At least two bound, concretely analyzed mode cells disagree. This does
    /// not by itself claim complete Cartesian mode coverage.
    NonUniformAcrossModes,
}

/// Orthogonal receipt axis whose claim remains unknown.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalOutcomeAxis {
    /// Empty-projection/under/over/well/mixed determination.
    Determination,
    /// Generic block-structural rank.
    StructuralRank,
    /// Hybrid-mode conditionality/coverage.
    Conditionality,
}

/// Axis-local unknown reason and optional deterministic resume point.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CausalUnknownAxisState {
    /// Exact axis left unknown.
    pub axis: CausalOutcomeAxis,
    /// Honest reason this axis is unknown.
    pub reason: CausalUnknownReason,
    /// Required for cancelled/budget-exhausted work and forbidden otherwise.
    pub resume_checkpoint: Option<CausalCheckpointRef>,
}

/// Joint-mode-specific child outcome constructed only from an already admitted
/// mode-cell receipt.
#[derive(Debug, Clone)]
pub struct ConditionalCausalOutcome {
    /// Canonical total assignment: exactly one branch for every declared
    /// condition in the graph.
    assignment: Vec<ConditionBranchSelection>,
    /// Exact normalized graph analyzed by the child.
    structure: IdentityReceipt<CausalStructureIdV1>,
    /// Exact provenance-bearing graph analyzed by the child.
    artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    /// Branch-local cardinality classification.
    determination: DeterminationClass,
    /// Branch-local structural-rank state.
    structural_rank: StructuralRankState,
    /// Branch-local progress reasons and deterministic resume points for any
    /// axes left unknown by an incomplete mode-cell analysis.
    unknown_axes: Vec<CausalUnknownAxisState>,
    /// Producer-independent child outcome identity.
    outcome: IdentityReceipt<CausalOutcomeIdV1>,
    /// Exact child receipt identity.
    receipt: IdentityReceipt<CausalizationReceiptIdV1>,
}

impl PartialEq for ConditionalCausalOutcome {
    fn eq(&self, other: &Self) -> bool {
        self.assignment == other.assignment
            && identity_receipt_adjudication_eq(self.structure, other.structure)
            && identity_receipt_adjudication_eq(self.artifact, other.artifact)
            && self.determination == other.determination
            && self.structural_rank == other.structural_rank
            && self.unknown_axes == other.unknown_axes
            && identity_receipt_adjudication_eq(self.outcome, other.outcome)
            && identity_receipt_adjudication_eq(self.receipt, other.receipt)
    }
}

impl Eq for ConditionalCausalOutcome {}

impl PartialOrd for ConditionalCausalOutcome {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ConditionalCausalOutcome {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.assignment
            .cmp(&other.assignment)
            .then_with(|| identity_receipt_adjudication_cmp(self.structure, other.structure))
            .then_with(|| identity_receipt_adjudication_cmp(self.artifact, other.artifact))
            .then_with(|| self.determination.cmp(&other.determination))
            .then_with(|| self.structural_rank.cmp(&other.structural_rank))
            .then_with(|| self.unknown_axes.cmp(&other.unknown_axes))
            .then_with(|| identity_receipt_adjudication_cmp(self.outcome, other.outcome))
            .then_with(|| identity_receipt_adjudication_cmp(self.receipt, other.receipt))
    }
}

impl Hash for ConditionalCausalOutcome {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.assignment.hash(state);
        hash_identity_receipt_adjudication(self.structure, state);
        hash_identity_receipt_adjudication(self.artifact, state);
        self.determination.hash(state);
        self.structural_rank.hash(state);
        self.unknown_axes.hash(state);
        hash_identity_receipt_adjudication(self.outcome, state);
        hash_identity_receipt_adjudication(self.receipt, state);
    }
}

/// Canonical schema marker for one exact set of mode-cell child receipts.
pub enum ConditionalOutcomeSetIdentitySchemaV1 {}

impl CanonicalSchema for ConditionalOutcomeSetIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.conditional-causal-outcome-set.v1";
    const NAME: &'static str = "conditional-causal-outcome-set";
    const VERSION: u32 = CAUSALIZATION_RECEIPT_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "canonical set of admitted mode-cell receipt identities and their exact assignments/axes";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "mode-cell-outcomes",
        WireType::OrderedBytes,
    )];
}

/// Strong identity of one exact conditional child-outcome set.
pub type ConditionalOutcomeSetIdV1 = EvidenceNodeId<ConditionalOutcomeSetIdentitySchemaV1>;

/// Canonical schema marker for an exact retained structural matching.
pub enum CausalMatchingSetIdentitySchemaV1 {}

impl CanonicalSchema for CausalMatchingSetIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.causal-matching-set.v1";
    const NAME: &'static str = "causal-matching-set";
    const VERSION: u32 = CAUSALIZATION_RECEIPT_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "canonical exact incidence/equation/derivative-variable matching witness";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "matching-pairs",
        WireType::OrderedBytes,
    )];
}

/// Strong identity of one exact structural matching witness.
pub type CausalMatchingSetIdV1 = EvidenceNodeId<CausalMatchingSetIdentitySchemaV1>;

/// Graph-, domain-, and witness-bound maximum-matching theorem commitment.
#[derive(Debug, Clone)]
pub struct MaximumMatchingBinding {
    structure: IdentityReceipt<CausalStructureIdV1>,
    artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    domain: CausalReceiptDomain,
    matching_set: IdentityReceipt<CausalMatchingSetIdV1>,
    certificate: MaximumMatchingCertificateRef,
    checker: CausalCheckerRef,
}

impl PartialEq for MaximumMatchingBinding {
    fn eq(&self, other: &Self) -> bool {
        identity_receipt_adjudication_eq(self.structure, other.structure)
            && identity_receipt_adjudication_eq(self.artifact, other.artifact)
            && self.domain == other.domain
            && identity_receipt_adjudication_eq(self.matching_set, other.matching_set)
            && self.certificate == other.certificate
            && self.checker == other.checker
    }
}

impl Eq for MaximumMatchingBinding {}

/// Exact reason a maximum-matching domain cannot describe the supplied graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaximumMatchingDomainError {
    /// An unconditional domain was supplied for a graph with conditions.
    ConditionalGraph,
    /// A mode cell was supplied for a graph without conditions.
    ConditionFreeGraph,
    /// The mode assignment was not exactly one selection per graph condition.
    AssignmentCardinality {
        /// Submitted selections.
        submitted: usize,
        /// Exact number of graph conditions.
        expected: usize,
        /// Public assignment cap.
        max: usize,
    },
    /// A canonical assignment entry named the wrong condition or an unknown branch.
    InvalidSelection {
        /// Canonical assignment position.
        index: usize,
    },
}

impl fmt::Display for MaximumMatchingDomainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ConditionalGraph => {
                f.write_str("an unconditional domain cannot describe a conditional graph")
            }
            Self::ConditionFreeGraph => {
                f.write_str("a mode cell cannot describe a condition-free graph")
            }
            Self::AssignmentCardinality {
                submitted,
                expected,
                max,
            } => write!(
                f,
                "mode assignment has {submitted} selections; expected {expected} within cap {max}"
            ),
            Self::InvalidSelection { index } => write!(
                f,
                "mode assignment entry {index} names the wrong condition or an unknown branch"
            ),
        }
    }
}

impl std::error::Error for MaximumMatchingDomainError {}

/// Exact reason a proposed matching witness is not an inhabited matching of
/// the supplied graph projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaximumMatchingWitnessError {
    /// The witness exceeded the public pair cap.
    PairLimit {
        /// Submitted pairs.
        submitted: usize,
        /// Public cap.
        max: usize,
    },
    /// Two canonical rows named the same incidence/endpoints/order tuple.
    DuplicatePair {
        /// Second row in canonical order.
        index: usize,
    },
    /// A row named an incidence outside the exact graph artifact.
    ForeignIncidence {
        /// Canonical witness row.
        index: usize,
    },
    /// The incidence is not an unknown-bearing matching edge.
    NonUnknownIncidence {
        /// Canonical witness row.
        index: usize,
    },
    /// The incidence's equation is not a matching equation.
    NonMatchingEquation {
        /// Canonical witness row.
        index: usize,
    },
    /// The row's equation, variable, or derivative order disagrees with its incidence.
    EndpointMismatch {
        /// Canonical witness row.
        index: usize,
    },
    /// The incidence is inactive in the bound mode cell.
    InactiveIncidence {
        /// Canonical witness row.
        index: usize,
    },
    /// Two rows match the same equation.
    DuplicateEquation {
        /// Second canonical witness row.
        index: usize,
    },
    /// Two rows match the same derivative-variable vertex.
    DuplicateVariable {
        /// Second canonical witness row.
        index: usize,
    },
}

impl fmt::Display for MaximumMatchingWitnessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PairLimit { submitted, max } => {
                write!(f, "matching has {submitted} pairs above cap {max}")
            }
            Self::DuplicatePair { index } => {
                write!(f, "canonical matching row {index} duplicates a prior pair")
            }
            Self::ForeignIncidence { index } => {
                write!(
                    f,
                    "matching row {index} names an incidence outside the graph"
                )
            }
            Self::NonUnknownIncidence { index } => {
                write!(f, "matching row {index} does not name an unknown incidence")
            }
            Self::NonMatchingEquation { index } => {
                write!(f, "matching row {index} does not name a matching equation")
            }
            Self::EndpointMismatch { index } => {
                write!(
                    f,
                    "matching row {index} disagrees with its incidence endpoints"
                )
            }
            Self::InactiveIncidence { index } => {
                write!(f, "matching row {index} is inactive in the bound mode cell")
            }
            Self::DuplicateEquation { index } => {
                write!(f, "matching row {index} reuses an equation")
            }
            Self::DuplicateVariable { index } => {
                write!(
                    f,
                    "matching row {index} reuses a derivative-variable vertex"
                )
            }
        }
    }
}

impl std::error::Error for MaximumMatchingWitnessError {}

/// Refusal from constructing a malformed maximum-matching binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaximumMatchingBindingError {
    /// Proposed witness is not a bounded, inhabited graph matching.
    InvalidMatchingSet(MaximumMatchingWitnessError),
    /// Proposed graph projection is not an exact inhabitable domain.
    InvalidDomain(MaximumMatchingDomainError),
    /// Hybrid summaries have no union-graph matching witness.
    HybridSummaryDomain,
    /// Canonical matching-set identity publication refused.
    Identity(CanonicalError),
    /// Ordered matching-row identity publication refused with exact progress.
    IdentityStream(CausalIdentityError),
}

impl fmt::Display for MaximumMatchingBindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMatchingSet(problem) => {
                write!(f, "maximum-matching witness is invalid: {problem}")
            }
            Self::InvalidDomain(problem) => {
                write!(f, "maximum-matching domain is invalid: {problem}")
            }
            Self::HybridSummaryDomain => {
                f.write_str("hybrid summaries cannot bind a union-graph maximum matching")
            }
            Self::Identity(error) => write!(f, "matching-set identity refused: {error}"),
            Self::IdentityStream(error) => {
                write!(f, "matching-set streamed identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for MaximumMatchingBindingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidMatchingSet(problem) => Some(problem),
            Self::InvalidDomain(problem) => Some(problem),
            Self::Identity(error) => Some(error),
            Self::IdentityStream(error) => Some(error),
            Self::HybridSummaryDomain => None,
        }
    }
}

fn maximum_matching_identity_error(error: CausalIdentityError) -> MaximumMatchingBindingError {
    if error.stream_refusal().is_some() {
        MaximumMatchingBindingError::IdentityStream(error)
    } else if let Some(canonical_error) = error.canonical_error() {
        MaximumMatchingBindingError::Identity(canonical_error)
    } else {
        MaximumMatchingBindingError::IdentityStream(error)
    }
}

impl MaximumMatchingBinding {
    /// Bind a maximum-matching theorem to one exact graph projection and one
    /// exact matching witness.
    ///
    /// # Errors
    /// Refuses a hybrid-summary domain, an oversized/invalid mode-cell domain,
    /// an oversized/invalid matching set, or a bounded canonical identity
    /// error, including cancellation before theorem-binding publication.
    pub fn new(
        graph: &AdmittedCausalGraph,
        mut domain: CausalReceiptDomain,
        matching: &[CausalMatchingPair],
        certificate: MaximumMatchingCertificateRef,
        checker: CausalCheckerRef,
        cx: &Cx<'_>,
    ) -> Result<Self, MaximumMatchingBindingError> {
        identity_materialization_checkpoint(cx, 0)
            .map_err(MaximumMatchingBindingError::Identity)?;
        canonicalize_maximum_matching_domain(graph, &mut domain, cx)?;
        let canonical_matching =
            canonicalize_maximum_matching_witness(graph, &domain, matching, cx)?;
        let matching_set = causal_matching_set_identity_cancellable(&canonical_matching, cx)
            .map_err(maximum_matching_identity_error)?;
        Ok(Self {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            domain,
            matching_set,
            certificate,
            checker,
        })
    }

    /// Exact maximum-matching certificate artifact reference.
    #[must_use]
    pub const fn certificate(&self) -> &MaximumMatchingCertificateRef {
        &self.certificate
    }

    /// Exact checker attestation bound to this theorem commitment.
    #[must_use]
    pub const fn checker(&self) -> &CausalCheckerRef {
        &self.checker
    }

    /// Exact graph structure bound by the theorem.
    #[must_use]
    pub const fn structure(&self) -> CausalStructureIdV1 {
        self.structure.id()
    }

    /// Exact provenance artifact bound by the theorem.
    #[must_use]
    pub const fn artifact(&self) -> CausalGraphArtifactIdV1 {
        self.artifact.id()
    }

    /// Complete normalized-structure identity receipt retained for collision
    /// adjudication as well as semantic-ID comparison.
    #[must_use]
    pub const fn structure_identity_receipt(&self) -> IdentityReceipt<CausalStructureIdV1> {
        self.structure
    }

    /// Complete provenance-artifact identity receipt retained for collision
    /// adjudication as well as semantic-ID comparison.
    #[must_use]
    pub const fn artifact_identity_receipt(&self) -> IdentityReceipt<CausalGraphArtifactIdV1> {
        self.artifact
    }

    /// Exact graph projection bound by the theorem.
    #[must_use]
    pub const fn domain(&self) -> &CausalReceiptDomain {
        &self.domain
    }

    /// Exact matching-witness identity bound by the theorem.
    #[must_use]
    pub const fn matching_set(&self) -> CausalMatchingSetIdV1 {
        self.matching_set.id()
    }

    /// Complete matching-set identity receipt retained for collision
    /// adjudication as well as semantic-ID comparison.
    #[must_use]
    pub const fn matching_set_identity_receipt(&self) -> IdentityReceipt<CausalMatchingSetIdV1> {
        self.matching_set
    }
}

fn canonicalize_maximum_matching_domain(
    graph: &AdmittedCausalGraph,
    domain: &mut CausalReceiptDomain,
    cx: &Cx<'_>,
) -> Result<(), MaximumMatchingBindingError> {
    match domain {
        CausalReceiptDomain::HybridSummary => {
            return Err(MaximumMatchingBindingError::HybridSummaryDomain);
        }
        CausalReceiptDomain::UnconditionalGraph if !graph.conditions().is_empty() => {
            return Err(MaximumMatchingBindingError::InvalidDomain(
                MaximumMatchingDomainError::ConditionalGraph,
            ));
        }
        CausalReceiptDomain::ModeCell { .. } if graph.conditions().is_empty() => {
            return Err(MaximumMatchingBindingError::InvalidDomain(
                MaximumMatchingDomainError::ConditionFreeGraph,
            ));
        }
        CausalReceiptDomain::UnconditionalGraph | CausalReceiptDomain::ModeCell { .. } => {}
    }
    let CausalReceiptDomain::ModeCell { assignment } = domain else {
        return Ok(());
    };
    if assignment.len() > MAX_CAUSAL_CONDITIONS || assignment.len() != graph.conditions().len() {
        return Err(MaximumMatchingBindingError::InvalidDomain(
            MaximumMatchingDomainError::AssignmentCardinality {
                submitted: assignment.len(),
                expected: graph.conditions().len(),
                max: MAX_CAUSAL_CONDITIONS,
            },
        ));
    }
    cancellable_sort(assignment, || identity_materialization_checkpoint(cx, 0))
        .map_err(MaximumMatchingBindingError::Identity)?;
    for (index, (selection, condition)) in assignment.iter().zip(graph.conditions()).enumerate() {
        identity_materialization_poll(cx, index, 0)
            .map_err(MaximumMatchingBindingError::Identity)?;
        if selection.condition != condition.condition
            || condition.branches.binary_search(&selection.branch).is_err()
        {
            return Err(MaximumMatchingBindingError::InvalidDomain(
                MaximumMatchingDomainError::InvalidSelection { index },
            ));
        }
    }
    Ok(())
}

fn canonicalize_maximum_matching_witness(
    graph: &AdmittedCausalGraph,
    domain: &CausalReceiptDomain,
    matching: &[CausalMatchingPair],
    cx: &Cx<'_>,
) -> Result<Vec<CausalMatchingPair>, MaximumMatchingBindingError> {
    if matching.len() > MAX_CAUSAL_MATCHING_PAIRS {
        return Err(MaximumMatchingBindingError::InvalidMatchingSet(
            MaximumMatchingWitnessError::PairLimit {
                submitted: matching.len(),
                max: MAX_CAUSAL_MATCHING_PAIRS,
            },
        ));
    }
    let mut canonical_matching = Vec::with_capacity(matching.len());
    for (index, pair) in matching.iter().enumerate() {
        identity_materialization_poll(cx, index, 0)
            .map_err(MaximumMatchingBindingError::Identity)?;
        canonical_matching.push(pair.clone());
    }
    cancellable_sort_by(
        &mut canonical_matching,
        compare_causal_matching_pairs_nominal,
        || identity_materialization_checkpoint(cx, 0),
    )
    .map_err(MaximumMatchingBindingError::Identity)?;
    for (index, pair) in canonical_matching.windows(2).enumerate() {
        identity_materialization_poll(cx, index, 0)
            .map_err(MaximumMatchingBindingError::Identity)?;
        if causal_matching_pair_nominal_eq(&pair[0], &pair[1]) {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::DuplicatePair { index: index + 1 },
            ));
        }
    }
    validate_maximum_matching_witness(graph, domain, &canonical_matching, cx)?;
    Ok(canonical_matching)
}

fn validate_maximum_matching_witness(
    graph: &AdmittedCausalGraph,
    domain: &CausalReceiptDomain,
    canonical_matching: &[CausalMatchingPair],
    cx: &Cx<'_>,
) -> Result<(), MaximumMatchingBindingError> {
    let mode_assignment = match domain {
        CausalReceiptDomain::ModeCell { assignment } => Some(assignment.as_slice()),
        CausalReceiptDomain::UnconditionalGraph => None,
        CausalReceiptDomain::HybridSummary => unreachable!("hybrid domain refused above"),
    };
    let mut matched_equations = BTreeSet::new();
    let mut matched_variables = BTreeSet::new();
    for (index, pair) in canonical_matching.iter().enumerate() {
        identity_materialization_poll(cx, index, 0)
            .map_err(MaximumMatchingBindingError::Identity)?;
        let Ok(incidence_index) = graph
            .incidences()
            .binary_search_by(|incidence| incidence.id.cmp(&pair.incidence))
        else {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::ForeignIncidence { index },
            ));
        };
        let incidence = &graph.incidences()[incidence_index];
        if incidence.solve_participation != SolveParticipation::Unknown {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::NonUnknownIncidence { index },
            ));
        }
        let equation_is_matching = graph
            .equations()
            .binary_search_by(|equation| equation.id.cmp(&incidence.equation))
            .ok()
            .is_some_and(|equation_index| {
                graph.equations()[equation_index].solve_participation
                    == EquationParticipation::Matching
            });
        if !equation_is_matching {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::NonMatchingEquation { index },
            ));
        }
        if pair.equation != incidence.equation
            || pair.variable.variable != incidence.variable
            || pair.variable.derivative_order != incidence.derivative_order
        {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::EndpointMismatch { index },
            ));
        }
        if !activation_active_for_binding(&incidence.activation, mode_assignment, cx)
            .map_err(MaximumMatchingBindingError::Identity)?
        {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::InactiveIncidence { index },
            ));
        }
        if !matched_equations.insert(pair.equation.clone()) {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::DuplicateEquation { index },
            ));
        }
        if !matched_variables.insert(pair.variable.clone()) {
            return Err(MaximumMatchingBindingError::InvalidMatchingSet(
                MaximumMatchingWitnessError::DuplicateVariable { index },
            ));
        }
    }
    Ok(())
}

fn activation_active_for_binding(
    activation: &ActivationDomain,
    assignment: Option<&[ConditionBranchSelection]>,
    cx: &Cx<'_>,
) -> Result<bool, CanonicalError> {
    match (activation, assignment) {
        (ActivationDomain::Always, _) => Ok(true),
        (ActivationDomain::Conditional { .. }, None) => Ok(false),
        (ActivationDomain::Conditional { cubes }, Some(assignment)) => {
            for (cube_index, cube) in cubes.iter().enumerate() {
                identity_materialization_poll(cx, cube_index, 0)?;
                let mut active = true;
                for (selection_index, selection) in cube.selections.iter().enumerate() {
                    identity_materialization_poll(cx, selection_index, 0)?;
                    if assignment.binary_search(selection).is_err() {
                        active = false;
                        break;
                    }
                }
                if active {
                    return Ok(true);
                }
            }
            Ok(false)
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum ConditionalCoverageClaim {
    ModeCells(IdentityReceipt<ConditionalOutcomeSetIdV1>),
    UniformTheorem {
        determination: DeterminationClass,
        structural_rank: StructuralRankState,
    },
}

impl PartialEq for ConditionalCoverageClaim {
    fn eq(&self, other: &Self) -> bool {
        conditional_coverage_claim_eq(*self, *other)
    }
}

impl Eq for ConditionalCoverageClaim {}

/// Source-bound conditional coverage or uniformity theorem commitment.
#[derive(Debug, Clone)]
pub struct ConditionalCoverageBinding {
    structure: IdentityReceipt<CausalStructureIdV1>,
    artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    claim: ConditionalCoverageClaim,
    certificate: ConditionalCoverageRef,
    checker: CausalCheckerRef,
}

impl PartialEq for ConditionalCoverageBinding {
    fn eq(&self, other: &Self) -> bool {
        identity_receipt_adjudication_eq(self.structure, other.structure)
            && identity_receipt_adjudication_eq(self.artifact, other.artifact)
            && conditional_coverage_claim_eq(self.claim, other.claim)
            && self.certificate == other.certificate
            && self.checker == other.checker
    }
}

impl Eq for ConditionalCoverageBinding {}

fn conditional_coverage_claim_eq(
    left: ConditionalCoverageClaim,
    right: ConditionalCoverageClaim,
) -> bool {
    match (left, right) {
        (ConditionalCoverageClaim::ModeCells(left), ConditionalCoverageClaim::ModeCells(right)) => {
            identity_receipt_adjudication_eq(left, right)
        }
        (
            ConditionalCoverageClaim::UniformTheorem {
                determination: left_determination,
                structural_rank: left_rank,
            },
            ConditionalCoverageClaim::UniformTheorem {
                determination: right_determination,
                structural_rank: right_rank,
            },
        ) => left_determination == right_determination && left_rank == right_rank,
        (
            ConditionalCoverageClaim::ModeCells(_),
            ConditionalCoverageClaim::UniformTheorem { .. },
        )
        | (
            ConditionalCoverageClaim::UniformTheorem { .. },
            ConditionalCoverageClaim::ModeCells(_),
        ) => false,
    }
}

/// Refusal from constructing a malformed coverage binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConditionalCoverageBindingError {
    /// Coverage requires at least one admitted mode-cell child.
    EmptyOutcomeSet,
    /// The submitted child count exceeded the public cap.
    OutcomeSetLimit {
        /// Submitted children.
        submitted: usize,
        /// Public cap.
        max: usize,
    },
    /// One child assignment or the aggregate assignments exceeded a public cap.
    SelectionLimit {
        /// Child in canonical order, or `None` for aggregate overflow.
        outcome_index: Option<usize>,
        /// Observed selection count.
        submitted: usize,
        /// Applicable public cap.
        max: usize,
    },
    /// A child receipt analyzed another graph.
    ForeignGraph {
        /// Child position after canonical ordering.
        outcome_index: usize,
    },
    /// A child left determination or rank unknown and therefore cannot witness a cell theorem.
    NonConcreteChild {
        /// Child position after canonical ordering.
        outcome_index: usize,
    },
    /// Two children describe the same exact mode assignment.
    DuplicateAssignment {
        /// Second child in canonical order.
        index: usize,
    },
    /// The graph's Cartesian mode domain overflowed addressable cardinality.
    CartesianDomainOverflow,
    /// The complete Cartesian mode domain exceeds the explicit-child envelope.
    ExplicitDomainTooLarge {
        /// Required mode cells.
        required_outcomes: usize,
        /// Public child cap.
        max_outcomes: usize,
        /// Required condition selections across all cells, saturated on overflow.
        required_selections: usize,
        /// Public aggregate selection cap.
        max_selections: usize,
    },
    /// Child count omitted or added cells relative to the exact Cartesian domain.
    IncompleteCartesianCover {
        /// Submitted unique mode cells.
        submitted: usize,
        /// Exact Cartesian cell count.
        expected: usize,
    },
    /// A canonical child did not equal the Cartesian cell at its ordinal.
    WrongCartesianCell {
        /// Canonical child position.
        index: usize,
    },
    /// The graph has no activation conditions and therefore no hybrid domain
    /// over which a coverage or uniformity theorem can range.
    InvalidGraphDomain,
    /// A uniform theorem must state concrete determination and rank axes.
    NonConcreteUniformClaim,
    /// Concrete determination and structural-rank axes contradicted each
    /// other's exact bipartition semantics.
    IncompatibleUniformClaim,
    /// Canonical child-set identity publication refused.
    Identity(CanonicalError),
    /// Ordered child-row identity publication refused with exact progress.
    IdentityStream(CausalIdentityError),
}

impl fmt::Display for ConditionalCoverageBindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyOutcomeSet => f.write_str("conditional coverage requires a child set"),
            Self::OutcomeSetLimit { submitted, max } => write!(
                f,
                "conditional coverage submitted {submitted} children above the {max}-child cap"
            ),
            Self::SelectionLimit {
                outcome_index,
                submitted,
                max,
            } => match outcome_index {
                Some(index) => write!(
                    f,
                    "conditional coverage canonical child {index} has {submitted} selections above the {max}-selection cap"
                ),
                None => write!(
                    f,
                    "conditional coverage aggregate has {submitted} selections above the {max}-selection cap"
                ),
            },
            Self::ForeignGraph { outcome_index } => write!(
                f,
                "conditional coverage canonical child {outcome_index} belongs to another causal graph"
            ),
            Self::NonConcreteChild { outcome_index } => write!(
                f,
                "conditional coverage canonical child {outcome_index} has a non-concrete outcome axis"
            ),
            Self::DuplicateAssignment { index } => write!(
                f,
                "conditional coverage canonical child {index} duplicates a mode assignment"
            ),
            Self::CartesianDomainOverflow => {
                f.write_str("conditional coverage Cartesian mode cardinality overflowed usize")
            }
            Self::ExplicitDomainTooLarge {
                required_outcomes,
                max_outcomes,
                required_selections,
                max_selections,
            } => write!(
                f,
                "conditional coverage requires {required_outcomes} cells/{required_selections} selections above caps {max_outcomes}/{max_selections}"
            ),
            Self::IncompleteCartesianCover {
                submitted,
                expected,
            } => write!(
                f,
                "conditional coverage submitted {submitted} cells but exact Cartesian coverage requires {expected}"
            ),
            Self::WrongCartesianCell { index } => write!(
                f,
                "conditional coverage canonical child {index} is not the expected Cartesian cell"
            ),
            Self::InvalidGraphDomain => {
                f.write_str("conditional coverage requires a graph with activation conditions")
            }
            Self::NonConcreteUniformClaim => {
                f.write_str("a uniform theorem requires concrete determination and rank axes")
            }
            Self::IncompatibleUniformClaim => {
                f.write_str("uniform determination and structural-rank axes are incompatible")
            }
            Self::Identity(error) => write!(f, "conditional outcome-set identity refused: {error}"),
            Self::IdentityStream(error) => {
                write!(
                    f,
                    "conditional outcome-set streamed identity refused: {error}"
                )
            }
        }
    }
}

impl std::error::Error for ConditionalCoverageBindingError {}

fn conditional_coverage_identity_error(
    error: CausalIdentityError,
) -> ConditionalCoverageBindingError {
    if error.stream_refusal().is_some() {
        ConditionalCoverageBindingError::IdentityStream(error)
    } else if let Some(canonical_error) = error.canonical_error() {
        ConditionalCoverageBindingError::Identity(canonical_error)
    } else {
        ConditionalCoverageBindingError::IdentityStream(error)
    }
}

impl ConditionalCoverageBinding {
    /// Bind a coverage certificate to an exact graph and exact canonical set
    /// of already-admitted mode-cell child receipts.
    ///
    /// # Errors
    /// Refuses a condition-free graph; an empty or oversized child set; child
    /// or aggregate selection overflow; foreign or non-concrete children;
    /// duplicate assignments; Cartesian cardinality overflow; a Cartesian
    /// domain larger than the explicit-child envelope; an omitted, extra, or
    /// non-Cartesian cell; or a bounded canonical-identity error, including
    /// cancellation before coverage-binding publication.
    pub fn for_mode_cells(
        graph: &AdmittedCausalGraph,
        outcomes: &[ConditionalCausalOutcome],
        certificate: ConditionalCoverageRef,
        checker: CausalCheckerRef,
        cx: &Cx<'_>,
    ) -> Result<Self, ConditionalCoverageBindingError> {
        identity_materialization_checkpoint(cx, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if graph.conditions().is_empty() {
            return Err(ConditionalCoverageBindingError::InvalidGraphDomain);
        }
        if outcomes.is_empty() {
            return Err(ConditionalCoverageBindingError::EmptyOutcomeSet);
        }
        if outcomes.len() > MAX_CAUSAL_CONDITIONAL_OUTCOMES {
            return Err(ConditionalCoverageBindingError::OutcomeSetLimit {
                submitted: outcomes.len(),
                max: MAX_CAUSAL_CONDITIONAL_OUTCOMES,
            });
        }
        let (strides, expected_outcomes) = explicit_cartesian_mode_domain(graph, cx)?;
        let canonical_outcomes = canonicalize_coverage_outcomes(graph, outcomes, cx)?;
        validate_cartesian_mode_cover(graph, &canonical_outcomes, &strides, expected_outcomes, cx)?;
        let outcome_set = conditional_outcome_set_identity_cancellable(&canonical_outcomes, cx)
            .map_err(conditional_coverage_identity_error)?;
        Ok(Self {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            claim: ConditionalCoverageClaim::ModeCells(outcome_set),
            certificate,
            checker,
        })
    }

    /// Bind a cross-mode uniformity theorem to one exact graph and concrete
    /// determination/rank claim.
    ///
    /// # Errors
    /// Returns [`ConditionalCoverageBindingError::InvalidGraphDomain`] when
    /// the graph has no activation condition,
    /// [`ConditionalCoverageBindingError::NonConcreteUniformClaim`] when
    /// either axis is unknown, or
    /// [`ConditionalCoverageBindingError::IncompatibleUniformClaim`] when the
    /// two concrete axes contradict bipartition semantics.
    pub fn for_uniform_theorem(
        graph: &AdmittedCausalGraph,
        determination: DeterminationClass,
        structural_rank: StructuralRankState,
        certificate: ConditionalCoverageRef,
        checker: CausalCheckerRef,
    ) -> Result<Self, ConditionalCoverageBindingError> {
        if graph.conditions().is_empty() {
            return Err(ConditionalCoverageBindingError::InvalidGraphDomain);
        }
        if determination == DeterminationClass::Unknown
            || structural_rank == StructuralRankState::Unknown
        {
            return Err(ConditionalCoverageBindingError::NonConcreteUniformClaim);
        }
        if !causal_axes_compatible(determination, structural_rank) {
            return Err(ConditionalCoverageBindingError::IncompatibleUniformClaim);
        }
        Ok(Self {
            structure: graph.structure_identity_receipt(),
            artifact: graph.artifact_identity_receipt(),
            claim: ConditionalCoverageClaim::UniformTheorem {
                determination,
                structural_rank,
            },
            certificate,
            checker,
        })
    }

    /// Exact certificate artifact reference.
    #[must_use]
    pub const fn certificate(&self) -> &ConditionalCoverageRef {
        &self.certificate
    }

    /// Exact checker attestation bound to this theorem commitment.
    #[must_use]
    pub const fn checker(&self) -> &CausalCheckerRef {
        &self.checker
    }

    /// Exact graph structure bound by the theorem.
    #[must_use]
    pub const fn structure(&self) -> CausalStructureIdV1 {
        self.structure.id()
    }

    /// Exact provenance artifact bound by the theorem.
    #[must_use]
    pub const fn artifact(&self) -> CausalGraphArtifactIdV1 {
        self.artifact.id()
    }

    /// Complete normalized-structure identity receipt retained for collision
    /// adjudication as well as semantic-ID comparison.
    #[must_use]
    pub const fn structure_identity_receipt(&self) -> IdentityReceipt<CausalStructureIdV1> {
        self.structure
    }

    /// Complete provenance-artifact identity receipt retained for collision
    /// adjudication as well as semantic-ID comparison.
    #[must_use]
    pub const fn artifact_identity_receipt(&self) -> IdentityReceipt<CausalGraphArtifactIdV1> {
        self.artifact
    }

    /// Canonical child-outcome set identity for a coverage claim.
    #[must_use]
    pub const fn outcome_set(&self) -> Option<ConditionalOutcomeSetIdV1> {
        match self.claim {
            ConditionalCoverageClaim::ModeCells(outcome_set) => Some(outcome_set.id()),
            ConditionalCoverageClaim::UniformTheorem { .. } => None,
        }
    }

    /// Complete child-outcome-set identity receipt retained for collision
    /// adjudication, or `None` for a uniform-theorem claim.
    #[must_use]
    pub const fn outcome_set_identity_receipt(
        &self,
    ) -> Option<IdentityReceipt<ConditionalOutcomeSetIdV1>> {
        match self.claim {
            ConditionalCoverageClaim::ModeCells(outcome_set) => Some(outcome_set),
            ConditionalCoverageClaim::UniformTheorem { .. } => None,
        }
    }

    /// Concrete axes for a cross-mode uniformity theorem.
    #[must_use]
    pub const fn uniform_axes(&self) -> Option<(DeterminationClass, StructuralRankState)> {
        match self.claim {
            ConditionalCoverageClaim::UniformTheorem {
                determination,
                structural_rank,
            } => Some((determination, structural_rank)),
            ConditionalCoverageClaim::ModeCells(_) => None,
        }
    }
}

fn explicit_cartesian_mode_domain(
    graph: &AdmittedCausalGraph,
    cx: &Cx<'_>,
) -> Result<(Vec<usize>, usize), ConditionalCoverageBindingError> {
    let mut strides = vec![0usize; graph.conditions().len()];
    let mut required_outcomes = 1usize;
    for reverse_offset in 0..graph.conditions().len() {
        identity_materialization_poll(cx, reverse_offset, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        let condition_index = graph.conditions().len() - 1 - reverse_offset;
        strides[condition_index] = required_outcomes;
        required_outcomes = required_outcomes
            .checked_mul(graph.conditions()[condition_index].branches.len())
            .ok_or(ConditionalCoverageBindingError::CartesianDomainOverflow)?;
    }
    let required_selections = required_outcomes.saturating_mul(graph.conditions().len());
    if required_outcomes > MAX_CAUSAL_CONDITIONAL_OUTCOMES
        || required_selections > MAX_CAUSAL_CONDITIONAL_SELECTIONS
    {
        return Err(ConditionalCoverageBindingError::ExplicitDomainTooLarge {
            required_outcomes,
            max_outcomes: MAX_CAUSAL_CONDITIONAL_OUTCOMES,
            required_selections,
            max_selections: MAX_CAUSAL_CONDITIONAL_SELECTIONS,
        });
    }
    Ok((strides, required_outcomes))
}

fn canonicalize_coverage_outcomes(
    graph: &AdmittedCausalGraph,
    outcomes: &[ConditionalCausalOutcome],
    cx: &Cx<'_>,
) -> Result<Vec<ConditionalCausalOutcome>, ConditionalCoverageBindingError> {
    let mut conditional_selections = 0usize;
    for (outcome_index, outcome) in outcomes.iter().enumerate() {
        identity_materialization_poll(cx, outcome_index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        conditional_selections = conditional_selections.saturating_add(outcome.assignment.len());
    }
    if conditional_selections > MAX_CAUSAL_CONDITIONAL_SELECTIONS {
        return Err(ConditionalCoverageBindingError::SelectionLimit {
            outcome_index: None,
            submitted: conditional_selections,
            max: MAX_CAUSAL_CONDITIONAL_SELECTIONS,
        });
    }

    // Canonicalize references before reporting any child-local defect. The
    // same invalid multiset must therefore choose the same rule and index
    // regardless of caller order, without cloning untrusted assignments first.
    let mut ordered_outcomes = outcomes.iter().collect::<Vec<_>>();
    cancellable_sort_by_fallible(
        &mut ordered_outcomes,
        |left, right| {
            compare_conditional_outcomes_cancellable(left, right, || {
                identity_materialization_checkpoint(cx, 0)
            })
        },
        || identity_materialization_checkpoint(cx, 0),
    )
    .map_err(ConditionalCoverageBindingError::Identity)?;

    for (outcome_index, outcome) in ordered_outcomes.iter().enumerate() {
        identity_materialization_poll(cx, outcome_index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if outcome.assignment.len() > MAX_CAUSAL_CONDITIONS {
            return Err(ConditionalCoverageBindingError::SelectionLimit {
                outcome_index: Some(outcome_index),
                submitted: outcome.assignment.len(),
                max: MAX_CAUSAL_CONDITIONS,
            });
        }
    }
    for (outcome_index, outcome) in ordered_outcomes.iter().enumerate() {
        identity_materialization_poll(cx, outcome_index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if !identity_receipt_adjudication_eq(outcome.structure, graph.structure_identity_receipt())
            || !identity_receipt_adjudication_eq(
                outcome.artifact,
                graph.artifact_identity_receipt(),
            )
        {
            return Err(ConditionalCoverageBindingError::ForeignGraph { outcome_index });
        }
    }
    for (outcome_index, outcome) in ordered_outcomes.iter().enumerate() {
        identity_materialization_poll(cx, outcome_index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if outcome.determination == DeterminationClass::Unknown
            || outcome.structural_rank == StructuralRankState::Unknown
        {
            return Err(ConditionalCoverageBindingError::NonConcreteChild { outcome_index });
        }
    }

    let mut canonical_outcomes = Vec::with_capacity(outcomes.len());
    for outcome in ordered_outcomes {
        let mut assignment = Vec::with_capacity(outcome.assignment.len());
        for (selection_index, selection) in outcome.assignment.iter().enumerate() {
            identity_materialization_poll(cx, selection_index, 0)
                .map_err(ConditionalCoverageBindingError::Identity)?;
            assignment.push(selection.clone());
        }
        canonical_outcomes.push(ConditionalCausalOutcome {
            assignment,
            structure: outcome.structure,
            artifact: outcome.artifact,
            determination: outcome.determination,
            structural_rank: outcome.structural_rank,
            unknown_axes: outcome.unknown_axes.clone(),
            outcome: outcome.outcome,
            receipt: outcome.receipt,
        });
    }
    for (index, pair) in canonical_outcomes.windows(2).enumerate() {
        identity_materialization_poll(cx, index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if cancellable_slice_eq(&pair[0].assignment, &pair[1].assignment, || {
            identity_materialization_checkpoint(cx, 0)
        })
        .map_err(ConditionalCoverageBindingError::Identity)?
        {
            return Err(ConditionalCoverageBindingError::DuplicateAssignment { index: index + 1 });
        }
    }
    Ok(canonical_outcomes)
}

fn validate_cartesian_mode_cover(
    graph: &AdmittedCausalGraph,
    canonical_outcomes: &[ConditionalCausalOutcome],
    strides: &[usize],
    expected_outcomes: usize,
    cx: &Cx<'_>,
) -> Result<(), ConditionalCoverageBindingError> {
    if canonical_outcomes.len() != expected_outcomes {
        return Err(ConditionalCoverageBindingError::IncompleteCartesianCover {
            submitted: canonical_outcomes.len(),
            expected: expected_outcomes,
        });
    }
    for (outcome_index, outcome) in canonical_outcomes.iter().enumerate() {
        identity_materialization_poll(cx, outcome_index, 0)
            .map_err(ConditionalCoverageBindingError::Identity)?;
        if outcome.assignment.len() != graph.conditions().len() {
            return Err(ConditionalCoverageBindingError::WrongCartesianCell {
                index: outcome_index,
            });
        }
        for (condition_index, (selection, condition)) in outcome
            .assignment
            .iter()
            .zip(graph.conditions())
            .enumerate()
        {
            identity_materialization_poll(cx, condition_index, 0)
                .map_err(ConditionalCoverageBindingError::Identity)?;
            let expected_branch = &condition.branches
                [(outcome_index / strides[condition_index]) % condition.branches.len()];
            if selection.condition != condition.condition || selection.branch != *expected_branch {
                return Err(ConditionalCoverageBindingError::WrongCartesianCell {
                    index: outcome_index,
                });
            }
        }
    }
    Ok(())
}

const fn causal_axes_compatible(
    determination: DeterminationClass,
    structural_rank: StructuralRankState,
) -> bool {
    match (determination, structural_rank) {
        (DeterminationClass::Unknown, _) | (_, StructuralRankState::Unknown) => true,
        (DeterminationClass::EmptyProjection, StructuralRankState::NotApplicable) => true,
        (
            DeterminationClass::WellDetermined
            | DeterminationClass::UnderDetermined
            | DeterminationClass::OverDetermined,
            StructuralRankState::FullRelativeToMinSide,
        )
        | (DeterminationClass::Mixed, StructuralRankState::Deficient)
        | (
            DeterminationClass::UnderDetermined | DeterminationClass::OverDetermined,
            StructuralRankState::NotApplicable,
        ) => true,
        _ => false,
    }
}

/// Explicit evidence state of a structurally admitted receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CausalReceiptEvidence {
    /// Local runtime/admission closure only; no independent checker or
    /// recursive schema-governance authority is claimed.
    Unverified,
    /// Exact checker artifact is referenced, but its presence alone does not
    /// authenticate or admit checker authority.
    CheckerReferenced(CausalCheckerRef),
}

/// Mutable, authority-free causalization receipt draft.
#[derive(Debug, Clone)]
pub struct CausalizationReceiptDraft {
    /// Complete normalized-structure receipt analyzed.
    pub structure: IdentityReceipt<CausalStructureIdV1>,
    /// Complete provenance-bearing graph-artifact receipt analyzed.
    pub artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    /// Analyzer implementation, budget, capability, seed, and determinism context.
    pub analysis: CausalAnalysisContext,
    /// Unconditional graph, exact mode cell, or cross-mode summary domain.
    pub domain: CausalReceiptDomain,
    /// Empty-projection/under/over/mixed/well/unknown axis.
    pub determination: DeterminationClass,
    /// Structural-rank axis.
    pub structural_rank: StructuralRankState,
    /// Conditionality axis.
    pub conditionality: Conditionality,
    /// Selected structural matching in arbitrary caller order.
    pub matching: Vec<CausalMatchingPair>,
    /// Exact unmatched equation complement.
    pub unmatched_equations: Vec<EquationId>,
    /// Exact unmatched unknown-vertex complement.
    pub unmatched_variables: Vec<DerivativeVariableKey>,
    /// Branch/mode child receipts.
    pub conditional_outcomes: Vec<ConditionalCausalOutcome>,
    /// Required before a non-min-side-saturating matching may assert maximum
    /// matching, deficiency, or a Mixed determination.
    pub maximum_matching_certificate: Option<MaximumMatchingBinding>,
    /// Required for complete conditional coverage or a mode-uniform theorem.
    pub conditional_coverage: Option<ConditionalCoverageBinding>,
    /// Canonical one-row-per-unknown-axis explanations and resume points.
    pub unknown_axes: Vec<CausalUnknownAxisState>,
    /// Honest checker/evidence state.
    pub evidence: CausalReceiptEvidence,
}

impl PartialEq for CausalizationReceiptDraft {
    fn eq(&self, other: &Self) -> bool {
        identity_receipt_adjudication_eq(self.structure, other.structure)
            && identity_receipt_adjudication_eq(self.artifact, other.artifact)
            && self.analysis == other.analysis
            && self.domain == other.domain
            && self.determination == other.determination
            && self.structural_rank == other.structural_rank
            && self.conditionality == other.conditionality
            && self.matching == other.matching
            && self.unmatched_equations == other.unmatched_equations
            && self.unmatched_variables == other.unmatched_variables
            && self.conditional_outcomes == other.conditional_outcomes
            && self.maximum_matching_certificate == other.maximum_matching_certificate
            && self.conditional_coverage == other.conditional_coverage
            && self.unknown_axes == other.unknown_axes
            && self.evidence == other.evidence
    }
}

impl Eq for CausalizationReceiptDraft {}

/// Canonical schema marker for producer-independent causal outcomes.
pub enum CausalOutcomeIdentitySchemaV1 {}

impl CanonicalSchema for CausalOutcomeIdentitySchemaV1 {
    const DOMAIN: &'static str = CAUSAL_OUTCOME_IDENTITY_DOMAIN_V1;
    const NAME: &'static str = "normalized-causal-outcome";
    const VERSION: u32 = CAUSAL_OUTCOME_IDENTITY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "producer-independent structural outcome over one normalized graph; analyzer, provenance artifact, certificates, and progress metadata excluded";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("causal-outcome-schema-version", WireType::U64),
        FieldSpec::child_of("causal-structure-id", &CAUSAL_STRUCTURE_CHILD),
        FieldSpec::required("causal-structure-receipt-adjudication", WireType::Bytes),
        FieldSpec::required("analysis-domain", WireType::Bytes),
        FieldSpec::required("outcome-axes", WireType::Bytes),
        FieldSpec::required("matching", WireType::OrderedBytes),
        FieldSpec::required("unmatched-equations", WireType::OrderedBytes),
        FieldSpec::required("unmatched-variables", WireType::OrderedBytes),
        FieldSpec::required("conditional-outcome-semantics", WireType::OrderedBytes),
    ];
}

/// Producer-independent identity of one normalized causal outcome.
pub type CausalOutcomeIdV1 = ProblemSemanticId<CausalOutcomeIdentitySchemaV1>;

const CAUSAL_OUTCOME_CHILD: ChildSpec = ChildSpec::for_identity::<CausalOutcomeIdV1>();

/// Canonical schema marker for causalization receipts.
///
/// This candidate grammar must be replaced by the recursive SCC-bundle
/// authority before that tracked blocker becomes terminal.
pub enum CausalizationReceiptIdentitySchemaV1 {}

impl CanonicalSchema for CausalizationReceiptIdentitySchemaV1 {
    const DOMAIN: &'static str = CAUSALIZATION_RECEIPT_IDENTITY_DOMAIN_V1;
    const NAME: &'static str = "causalization-structural-receipt";
    const VERSION: u32 = CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "locally admission-closed structural outcome; recursive schema governance and numerical-rank, solvability, DAE-index, or physical-causality authority excluded";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("causalization-receipt-schema-version", WireType::U64),
        FieldSpec::child_of("causal-structure-id", &CAUSAL_STRUCTURE_CHILD),
        FieldSpec::required("causal-structure-receipt-adjudication", WireType::Bytes),
        FieldSpec::child_of("causal-graph-artifact-id", &CAUSAL_GRAPH_ARTIFACT_CHILD),
        FieldSpec::required(
            "causal-graph-artifact-receipt-adjudication",
            WireType::Bytes,
        ),
        FieldSpec::required("analysis-context", WireType::Bytes),
        FieldSpec::required("analysis-domain", WireType::Bytes),
        FieldSpec::required("outcome-axes", WireType::Bytes),
        FieldSpec::required("matching", WireType::OrderedBytes),
        FieldSpec::required("unmatched-equations", WireType::OrderedBytes),
        FieldSpec::required("unmatched-variables", WireType::OrderedBytes),
        FieldSpec::required("conditional-outcomes", WireType::OrderedBytes),
        FieldSpec::optional_bytes("maximum-matching-certificate"),
        FieldSpec::optional_bytes("conditional-coverage"),
        FieldSpec::required("unknown-axes", WireType::OrderedBytes),
        FieldSpec::required("evidence-state", WireType::Variant),
        FieldSpec::child_of("normalized-causal-outcome-id", &CAUSAL_OUTCOME_CHILD),
        FieldSpec::required(
            "normalized-causal-outcome-receipt-adjudication",
            WireType::Bytes,
        ),
    ];
}

/// Strong identity of one causalization receipt artifact.
pub type CausalizationReceiptIdV1 = EvidenceNodeId<CausalizationReceiptIdentitySchemaV1>;

/// Stable receipt-admission rule vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CausalReceiptRule {
    /// A receipt collection exceeded its explicit cap.
    ResourceLimit = 1,
    /// Receipt targeted a different structure or graph artifact.
    GraphIdentityMismatch = 2,
    /// Matching pair named a missing or known/condition-only vertex.
    UnknownMatchingEndpoint = 3,
    /// Matching pair was not an admitted incidence.
    NonIncidenceMatch = 4,
    /// One equation or derivative-variable vertex was matched twice.
    DuplicateMatchingEndpoint = 5,
    /// Supplied unmatched sets were not exact matching complements.
    UnmatchedSetMismatch = 6,
    /// Determination/rank axes contradicted the retained witness.
    OutcomeAxisMismatch = 7,
    /// Conditional child outcomes were missing, duplicated, or foreign.
    ConditionalCoverageMismatch = 8,
    /// Unknown reason presence contradicted the three outcome axes.
    UnknownReasonMismatch = 9,
    /// Bounded canonical identity publication refused.
    Identity = 10,
    /// Admission observed cancellation and published no receipt identity.
    Cancelled = 11,
}

impl CausalReceiptRule {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ResourceLimit => "CausalReceiptResourceLimit",
            Self::GraphIdentityMismatch => "CausalReceiptGraphIdentityMismatch",
            Self::UnknownMatchingEndpoint => "CausalReceiptUnknownMatchingEndpoint",
            Self::NonIncidenceMatch => "CausalReceiptNonIncidenceMatch",
            Self::DuplicateMatchingEndpoint => "CausalReceiptDuplicateMatchingEndpoint",
            Self::UnmatchedSetMismatch => "CausalReceiptUnmatchedSetMismatch",
            Self::OutcomeAxisMismatch => "CausalReceiptOutcomeAxisMismatch",
            Self::ConditionalCoverageMismatch => "CausalReceiptConditionalCoverageMismatch",
            Self::UnknownReasonMismatch => "CausalReceiptUnknownReasonMismatch",
            Self::Identity => "CausalReceiptIdentity",
            Self::Cancelled => "CausalReceiptCancelled",
        }
    }
}

/// Local subject of one receipt finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalReceiptSubject {
    /// Complete receipt.
    Receipt,
    /// Normalized structure identity supplied by the receipt.
    StructureIdentity,
    /// Provenance-bearing graph artifact identity supplied by the receipt.
    GraphArtifactIdentity,
    /// Submitted unmatched-equation complement.
    UnmatchedEquations,
    /// Submitted unmatched-variable complement.
    UnmatchedVariables,
    /// One matching pair.
    Matching(CausalMatchingPair),
    /// One equation.
    Equation(EquationId),
    /// One derivative-variable vertex.
    Variable(DerivativeVariableKey),
    /// One condition branch.
    Condition {
        /// Condition definition.
        condition: ActivationConditionRef,
        /// Branch/mode identity.
        branch: ActivationBranchRef,
    },
}

/// One sorted receipt-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CausalReceiptFinding {
    rule: CausalReceiptRule,
    subject: CausalReceiptSubject,
}

impl CausalReceiptFinding {
    fn new(rule: CausalReceiptRule, subject: CausalReceiptSubject) -> Self {
        Self { rule, subject }
    }

    /// Stable rule.
    #[must_use]
    pub const fn rule(&self) -> CausalReceiptRule {
        self.rule
    }

    /// Local subject.
    #[must_use]
    pub const fn subject(&self) -> &CausalReceiptSubject {
        &self.subject
    }

    /// Stable code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }
}

/// Complete deterministic receipt refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalReceiptRefusal {
    findings: Vec<CausalReceiptFinding>,
    identity_error: Option<CanonicalError>,
    identity_allocation_refusal: Option<CausalIdentityAllocationRefusal>,
    identity_stream_refusal: Option<CausalIdentityStreamRefusal>,
}

impl CausalReceiptRefusal {
    /// Sorted, duplicate-free findings.
    #[must_use]
    pub fn findings(&self) -> &[CausalReceiptFinding] {
        &self.findings
    }

    /// Canonical publication error, absent for the distinct fallible-allocation
    /// refusal exposed by [`Self::identity_allocation_refusal`].
    #[must_use]
    pub const fn identity_error(&self) -> Option<&CanonicalError> {
        self.identity_error.as_ref()
    }

    /// Exact fallible row-storage refusal, when identity publication could not
    /// allocate the already-admitted logical row.
    #[must_use]
    pub const fn identity_allocation_refusal(&self) -> Option<CausalIdentityAllocationRefusal> {
        self.identity_allocation_refusal
    }

    /// Exact ordered-row progress when identity publication refused inside a
    /// streamed collection.
    #[must_use]
    pub const fn identity_stream_refusal(&self) -> Option<CausalIdentityStreamRefusal> {
        self.identity_stream_refusal
    }
}

impl fmt::Display for CausalReceiptRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "causalization receipt refused with {} deterministic finding(s)",
            self.findings.len()
        )
    }
}

impl std::error::Error for CausalReceiptRefusal {}

/// Structurally admitted causalization receipt.
#[derive(Debug, Clone)]
pub struct AdmittedCausalizationReceipt {
    outcome_receipt: IdentityReceipt<CausalOutcomeIdV1>,
    structure: IdentityReceipt<CausalStructureIdV1>,
    artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    analysis: CausalAnalysisContext,
    domain: CausalReceiptDomain,
    determination: DeterminationClass,
    structural_rank: StructuralRankState,
    conditionality: Conditionality,
    matching: Vec<CausalMatchingPair>,
    unmatched_equations: Vec<EquationId>,
    unmatched_variables: Vec<DerivativeVariableKey>,
    conditional_outcomes: Vec<ConditionalCausalOutcome>,
    maximum_matching_certificate: Option<MaximumMatchingBinding>,
    conditional_coverage: Option<ConditionalCoverageBinding>,
    unknown_axes: Vec<CausalUnknownAxisState>,
    evidence: CausalReceiptEvidence,
    receipt: IdentityReceipt<CausalizationReceiptIdV1>,
}

impl PartialEq for AdmittedCausalizationReceipt {
    fn eq(&self, other: &Self) -> bool {
        identity_receipt_adjudication_eq(self.outcome_receipt, other.outcome_receipt)
            && identity_receipt_adjudication_eq(self.structure, other.structure)
            && identity_receipt_adjudication_eq(self.artifact, other.artifact)
            && self.analysis == other.analysis
            && self.domain == other.domain
            && self.determination == other.determination
            && self.structural_rank == other.structural_rank
            && self.conditionality == other.conditionality
            && self.matching == other.matching
            && self.unmatched_equations == other.unmatched_equations
            && self.unmatched_variables == other.unmatched_variables
            && self.conditional_outcomes == other.conditional_outcomes
            && self.maximum_matching_certificate == other.maximum_matching_certificate
            && self.conditional_coverage == other.conditional_coverage
            && self.unknown_axes == other.unknown_axes
            && self.evidence == other.evidence
            && identity_receipt_adjudication_eq(self.receipt, other.receipt)
    }
}

impl Eq for AdmittedCausalizationReceipt {}

impl AdmittedCausalizationReceipt {
    /// Producer-independent normalized outcome identity.
    #[must_use]
    pub const fn outcome_identity(&self) -> CausalOutcomeIdV1 {
        self.outcome_receipt.id()
    }

    /// Complete normalized-outcome identity receipt.
    #[must_use]
    pub const fn outcome_identity_receipt(&self) -> IdentityReceipt<CausalOutcomeIdV1> {
        self.outcome_receipt
    }

    /// Strong receipt artifact identity.
    #[must_use]
    pub const fn identity(&self) -> CausalizationReceiptIdV1 {
        self.receipt.id()
    }

    /// Complete canonical identity receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<CausalizationReceiptIdV1> {
        self.receipt
    }

    /// Exact normalized structure.
    #[must_use]
    pub const fn structure(&self) -> CausalStructureIdV1 {
        self.structure.id()
    }

    /// Exact provenance-bearing graph artifact.
    #[must_use]
    pub const fn artifact(&self) -> CausalGraphArtifactIdV1 {
        self.artifact.id()
    }

    /// Complete normalized-structure receipt for collision adjudication.
    #[must_use]
    pub const fn structure_identity_receipt(&self) -> IdentityReceipt<CausalStructureIdV1> {
        self.structure
    }

    /// Complete provenance-artifact receipt for collision adjudication.
    #[must_use]
    pub const fn artifact_identity_receipt(&self) -> IdentityReceipt<CausalGraphArtifactIdV1> {
        self.artifact
    }

    /// Analyzer identity.
    #[must_use]
    pub const fn analysis(&self) -> &CausalAnalysisContext {
        &self.analysis
    }

    /// Exact graph projection analyzed.
    #[must_use]
    pub const fn domain(&self) -> &CausalReceiptDomain {
        &self.domain
    }

    /// Determination axis.
    #[must_use]
    pub const fn determination(&self) -> DeterminationClass {
        self.determination
    }

    /// Structural-rank axis.
    #[must_use]
    pub const fn structural_rank(&self) -> StructuralRankState {
        self.structural_rank
    }

    /// Conditionality axis.
    #[must_use]
    pub const fn conditionality(&self) -> Conditionality {
        self.conditionality
    }

    /// Canonical matching pairs.
    #[must_use]
    pub fn matching(&self) -> &[CausalMatchingPair] {
        &self.matching
    }

    /// Exact unmatched equation complement.
    #[must_use]
    pub fn unmatched_equations(&self) -> &[EquationId] {
        &self.unmatched_equations
    }

    /// Exact unmatched unknown-vertex complement.
    #[must_use]
    pub fn unmatched_variables(&self) -> &[DerivativeVariableKey] {
        &self.unmatched_variables
    }

    /// Canonical condition-specific child outcomes.
    #[must_use]
    pub fn conditional_outcomes(&self) -> &[ConditionalCausalOutcome] {
        &self.conditional_outcomes
    }

    /// Maximum-matching certificate reference, if asserted.
    #[must_use]
    pub const fn maximum_matching_certificate(&self) -> Option<&MaximumMatchingBinding> {
        self.maximum_matching_certificate.as_ref()
    }

    /// Conditional-domain coverage/uniformity certificate reference.
    #[must_use]
    pub const fn conditional_coverage(&self) -> Option<&ConditionalCoverageBinding> {
        self.conditional_coverage.as_ref()
    }

    /// Canonical explanations/checkpoints for every unknown axis.
    #[must_use]
    pub fn unknown_axes(&self) -> &[CausalUnknownAxisState] {
        &self.unknown_axes
    }

    /// Honest evidence/checker state.
    #[must_use]
    pub const fn evidence(&self) -> &CausalReceiptEvidence {
        &self.evidence
    }
}

/// Refusal from projecting an admitted receipt into a hybrid-summary child.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConditionalOutcomeError {
    /// Only an exact [`CausalReceiptDomain::ModeCell`] receipt may become a
    /// conditional child.
    NotModeCell,
    /// Construction was cancelled before all child assignment/progress state
    /// was copied and the completed object passed its final checkpoint.
    Cancelled,
}

impl fmt::Display for ConditionalOutcomeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotModeCell => {
                f.write_str("conditional outcomes require an admitted mode-cell receipt")
            }
            Self::Cancelled => f.write_str("conditional outcome construction was cancelled"),
        }
    }
}

impl std::error::Error for ConditionalOutcomeError {}

impl ConditionalCausalOutcome {
    /// Construct a non-forgeable child commitment from an admitted mode-cell
    /// receipt. Graph identity, assignment, outcome axes, unknown-axis resume
    /// state, and typed receipt identity are copied together and cannot be
    /// independently substituted later.
    ///
    /// # Errors
    /// Refuses an unconditional-graph or hybrid-summary receipt, or
    /// cancellation before the assignment/progress-state copy completes.
    pub fn from_mode_cell(
        child: &AdmittedCausalizationReceipt,
        cx: &Cx<'_>,
    ) -> Result<Self, ConditionalOutcomeError> {
        cx.checkpoint()
            .map_err(|_| ConditionalOutcomeError::Cancelled)?;
        let CausalReceiptDomain::ModeCell { assignment } = child.domain() else {
            return Err(ConditionalOutcomeError::NotModeCell);
        };
        let mut copied_assignment = Vec::with_capacity(assignment.len());
        for (index, selection) in assignment.iter().enumerate() {
            if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
                cx.checkpoint()
                    .map_err(|_| ConditionalOutcomeError::Cancelled)?;
            }
            copied_assignment.push(selection.clone());
        }
        let mut copied_unknown_axes = Vec::with_capacity(child.unknown_axes().len());
        for (index, state) in child.unknown_axes().iter().enumerate() {
            if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
                cx.checkpoint()
                    .map_err(|_| ConditionalOutcomeError::Cancelled)?;
            }
            copied_unknown_axes.push(state.clone());
        }
        cx.checkpoint()
            .map_err(|_| ConditionalOutcomeError::Cancelled)?;
        Ok(Self {
            assignment: copied_assignment,
            structure: child.structure_identity_receipt(),
            artifact: child.artifact_identity_receipt(),
            determination: child.determination(),
            structural_rank: child.structural_rank(),
            unknown_axes: copied_unknown_axes,
            outcome: child.outcome_identity_receipt(),
            receipt: child.identity_receipt(),
        })
    }

    /// Exact total mode assignment analyzed by the child.
    #[must_use]
    pub fn assignment(&self) -> &[ConditionBranchSelection] {
        &self.assignment
    }

    /// Exact normalized graph analyzed by the child.
    #[must_use]
    pub const fn structure(&self) -> CausalStructureIdV1 {
        self.structure.id()
    }

    /// Complete normalized-structure receipt retained for collision
    /// adjudication.
    #[must_use]
    pub const fn structure_identity_receipt(&self) -> IdentityReceipt<CausalStructureIdV1> {
        self.structure
    }

    /// Exact provenance-bearing graph analyzed by the child.
    #[must_use]
    pub const fn artifact(&self) -> CausalGraphArtifactIdV1 {
        self.artifact.id()
    }

    /// Complete provenance-artifact receipt retained for collision
    /// adjudication.
    #[must_use]
    pub const fn artifact_identity_receipt(&self) -> IdentityReceipt<CausalGraphArtifactIdV1> {
        self.artifact
    }

    /// Child determination axis.
    #[must_use]
    pub const fn determination(&self) -> DeterminationClass {
        self.determination
    }

    /// Child structural-rank axis.
    #[must_use]
    pub const fn structural_rank(&self) -> StructuralRankState {
        self.structural_rank
    }

    /// Branch-local explanations and resume points for unknown axes.
    #[must_use]
    pub fn unknown_axes(&self) -> &[CausalUnknownAxisState] {
        &self.unknown_axes
    }

    /// Exact typed child receipt identity.
    #[must_use]
    pub const fn receipt(&self) -> CausalizationReceiptIdV1 {
        self.receipt.id()
    }

    /// Producer-independent child outcome identity.
    #[must_use]
    pub const fn outcome(&self) -> CausalOutcomeIdV1 {
        self.outcome.id()
    }

    /// Complete child evidence-receipt identity metadata.
    #[must_use]
    pub const fn receipt_identity_receipt(&self) -> IdentityReceipt<CausalizationReceiptIdV1> {
        self.receipt
    }

    /// Complete producer-independent child outcome identity metadata.
    #[must_use]
    pub const fn outcome_identity_receipt(&self) -> IdentityReceipt<CausalOutcomeIdV1> {
        self.outcome
    }
}

impl CausalizationReceiptDraft {
    /// Admit this receipt against the exact graph it claims to analyze.
    ///
    /// Admission checks bounded internal closure only. In particular, it does
    /// not establish maximum matching, numerical rank, physical causality, or
    /// checker authenticity.
    ///
    /// # Errors
    /// Refuses foreign identities, non-incidence matches, duplicate endpoints,
    /// inexact unmatched complements, inconsistent axes, malformed conditional
    /// coverage, resource/diagnostic-budget overflow, cancellation, or bounded
    /// identity publication.
    pub fn admit_against(
        self,
        graph: &AdmittedCausalGraph,
        cx: &Cx<'_>,
    ) -> Result<AdmittedCausalizationReceipt, CausalReceiptRefusal> {
        self.admit_with_decision(graph, cx).into_result()
    }

    /// Admit while retaining pre-canonicalization collection counts and a
    /// stable top-level outcome code for tracing and ledger adapters.
    #[must_use]
    pub fn admit_with_decision(
        self,
        graph: &AdmittedCausalGraph,
        cx: &Cx<'_>,
    ) -> CausalReceiptAdmissionDecision {
        match submitted_receipt_counts(&self, cx) {
            Ok(submitted) => CausalReceiptAdmissionDecision {
                submitted,
                result: admit_causalization_receipt(self, graph, submitted, cx),
            },
            Err((submitted, refusal)) => CausalReceiptAdmissionDecision {
                submitted,
                result: Err(refusal),
            },
        }
    }
}

/// Submitted receipt collection counts retained on refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CausalReceiptSubmittedCounts {
    /// Whether the complete telemetry pass reached its end without
    /// cancellation or an early resource-cap refusal. Top-level collection
    /// lengths and `domain_selections` remain exact when this is `false`, but
    /// `conditional_selections` is then only a prefix observation.
    pub complete: bool,
    /// Submitted matching rows.
    pub matching: usize,
    /// Submitted unmatched-equation rows.
    pub unmatched_equations: usize,
    /// Submitted unmatched derivative-variable rows.
    pub unmatched_variables: usize,
    /// Submitted mode-specific child outcomes.
    pub conditional_outcomes: usize,
    /// Aggregate condition-to-branch selections.
    pub conditional_selections: usize,
    /// Condition selections in a mode-cell analysis domain.
    pub domain_selections: usize,
    /// Submitted axis-local unknown explanations.
    pub unknown_axes: usize,
}

fn submitted_receipt_counts(
    draft: &CausalizationReceiptDraft,
    cx: &Cx<'_>,
) -> Result<CausalReceiptSubmittedCounts, (CausalReceiptSubmittedCounts, CausalReceiptRefusal)> {
    let domain_selections = match &draft.domain {
        CausalReceiptDomain::ModeCell { assignment } => assignment.len(),
        CausalReceiptDomain::UnconditionalGraph | CausalReceiptDomain::HybridSummary => 0,
    };
    let mut counts = CausalReceiptSubmittedCounts {
        complete: false,
        matching: draft.matching.len(),
        unmatched_equations: draft.unmatched_equations.len(),
        unmatched_variables: draft.unmatched_variables.len(),
        conditional_outcomes: draft.conditional_outcomes.len(),
        conditional_selections: 0,
        domain_selections,
        unknown_axes: draft.unknown_axes.len(),
    };
    if let Err(refusal) = receipt_checkpoint(cx) {
        return Err((counts, refusal));
    }
    if counts.matching > MAX_CAUSAL_MATCHING_PAIRS
        || counts.unmatched_equations > MAX_CAUSAL_EQUATIONS
        || counts.unmatched_variables > MAX_CAUSAL_DERIVATIVE_VERTICES
        || counts.conditional_outcomes > MAX_CAUSAL_CONDITIONAL_OUTCOMES
        || counts.domain_selections > MAX_CAUSAL_CONDITIONS
        || counts.unknown_axes > 3
    {
        return Err((counts, resource_receipt_refusal()));
    }
    for (index, outcome) in draft.conditional_outcomes.iter().enumerate() {
        if let Err(refusal) = receipt_poll(cx, index) {
            return Err((counts, refusal));
        }
        if outcome.assignment.len() > MAX_CAUSAL_CONDITIONS {
            return Err((counts, resource_receipt_refusal()));
        }
        counts.conditional_selections = counts
            .conditional_selections
            .saturating_add(outcome.assignment.len());
        if counts
            .conditional_selections
            .saturating_add(counts.domain_selections)
            > MAX_CAUSAL_CONDITIONAL_SELECTIONS
        {
            return Err((counts, resource_receipt_refusal()));
        }
    }
    counts.complete = true;
    Ok(counts)
}

/// Cancellation-correct, structured receipt-admission outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalReceiptAdmissionDecision {
    submitted: CausalReceiptSubmittedCounts,
    result: Result<AdmittedCausalizationReceipt, CausalReceiptRefusal>,
}

impl CausalReceiptAdmissionDecision {
    /// Counts observed before canonicalization.
    #[must_use]
    pub const fn submitted_counts(&self) -> CausalReceiptSubmittedCounts {
        self.submitted
    }

    /// Stable top-level decision code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "CausalReceiptAdmitted",
            Err(_) => "CausalReceiptRefused",
        }
    }

    /// Borrow the conventional result.
    #[must_use]
    pub fn result(&self) -> Result<&AdmittedCausalizationReceipt, &CausalReceiptRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision.
    #[must_use]
    pub fn into_result(self) -> Result<AdmittedCausalizationReceipt, CausalReceiptRefusal> {
        self.result
    }
}

fn receipt_refusal(
    mut findings: Vec<CausalReceiptFinding>,
    identity_error: Option<CanonicalError>,
) -> CausalReceiptRefusal {
    findings.sort();
    findings.dedup();
    debug_assert!(!findings.is_empty());
    CausalReceiptRefusal {
        findings,
        identity_error,
        identity_allocation_refusal: None,
        identity_stream_refusal: None,
    }
}

fn identity_receipt_refusal(error: CausalIdentityError) -> CausalReceiptRefusal {
    let rule = if error.is_cancelled() {
        CausalReceiptRule::Cancelled
    } else {
        CausalReceiptRule::Identity
    };
    let mut refusal = receipt_refusal(
        vec![CausalReceiptFinding::new(
            rule,
            CausalReceiptSubject::Receipt,
        )],
        error.canonical_error(),
    );
    refusal.identity_allocation_refusal = error.allocation_refusal();
    refusal.identity_stream_refusal = error.stream_refusal();
    refusal
}

fn receipt_refusal_cancellable(
    mut findings: Vec<CausalReceiptFinding>,
    cx: &Cx<'_>,
) -> Result<CausalReceiptRefusal, CausalReceiptRefusal> {
    cancellable_sort(&mut findings, || receipt_checkpoint(cx))?;
    cancellable_dedup(&mut findings, || receipt_checkpoint(cx))?;
    debug_assert!(!findings.is_empty());
    Ok(CausalReceiptRefusal {
        findings,
        identity_error: None,
        identity_allocation_refusal: None,
        identity_stream_refusal: None,
    })
}

fn cancelled_receipt_refusal() -> CausalReceiptRefusal {
    receipt_refusal(
        vec![CausalReceiptFinding::new(
            CausalReceiptRule::Cancelled,
            CausalReceiptSubject::Receipt,
        )],
        None,
    )
}

fn resource_receipt_refusal() -> CausalReceiptRefusal {
    receipt_refusal(
        vec![CausalReceiptFinding::new(
            CausalReceiptRule::ResourceLimit,
            CausalReceiptSubject::Receipt,
        )],
        None,
    )
}

fn enforce_receipt_finding_budget(
    findings: &[CausalReceiptFinding],
) -> Result<(), CausalReceiptRefusal> {
    if findings.len() > MAX_CAUSAL_RECEIPT_FINDINGS {
        Err(resource_receipt_refusal())
    } else {
        Ok(())
    }
}

fn receipt_checkpoint(cx: &Cx<'_>) -> Result<(), CausalReceiptRefusal> {
    cx.checkpoint().map_err(|_| cancelled_receipt_refusal())
}

fn receipt_poll(cx: &Cx<'_>, index: usize) -> Result<(), CausalReceiptRefusal> {
    if index.is_multiple_of(CAUSAL_CANCELLATION_POLL_STRIDE) {
        receipt_checkpoint(cx)?;
    }
    Ok(())
}

fn receipt_condition_domain(
    graph: &AdmittedCausalGraph,
    cx: &Cx<'_>,
) -> Result<BTreeMap<ActivationConditionRef, BTreeSet<ActivationBranchRef>>, CausalReceiptRefusal> {
    let mut domain = BTreeMap::new();
    for (condition_index, condition) in graph.conditions.iter().enumerate() {
        receipt_poll(cx, condition_index)?;
        let mut branches = BTreeSet::new();
        for (branch_index, branch) in condition.branches.iter().enumerate() {
            receipt_poll(cx, branch_index)?;
            branches.insert(branch.clone());
        }
        domain.insert(condition.condition.clone(), branches);
    }
    receipt_checkpoint(cx)?;
    Ok(domain)
}

fn receipt_domains_equal_cancellable(
    left: &CausalReceiptDomain,
    right: &CausalReceiptDomain,
    cx: &Cx<'_>,
) -> Result<bool, CausalReceiptRefusal> {
    match (left, right) {
        (CausalReceiptDomain::UnconditionalGraph, CausalReceiptDomain::UnconditionalGraph)
        | (CausalReceiptDomain::HybridSummary, CausalReceiptDomain::HybridSummary) => Ok(true),
        (
            CausalReceiptDomain::ModeCell { assignment: left },
            CausalReceiptDomain::ModeCell { assignment: right },
        ) => cancellable_slice_eq(left, right, || receipt_checkpoint(cx)),
        _ => Ok(false),
    }
}

#[allow(clippy::too_many_lines)]
fn admit_causalization_receipt(
    mut draft: CausalizationReceiptDraft,
    graph: &AdmittedCausalGraph,
    counts: CausalReceiptSubmittedCounts,
    cx: &Cx<'_>,
) -> Result<AdmittedCausalizationReceipt, CausalReceiptRefusal> {
    // frankensim-unratified-candidate-identity:frankensim-leapfrog-2026-program-i94v.1.2.12:admit_causalization_receipt
    receipt_checkpoint(cx)?;
    debug_assert!(counts.complete);
    cancellable_sort_by(
        &mut draft.matching,
        compare_causal_matching_pairs_nominal,
        || receipt_checkpoint(cx),
    )?;
    cancellable_sort(&mut draft.unmatched_equations, || receipt_checkpoint(cx))?;
    cancellable_sort(&mut draft.unmatched_variables, || receipt_checkpoint(cx))?;
    for (index, outcome) in draft.conditional_outcomes.iter_mut().enumerate() {
        receipt_poll(cx, index)?;
        cancellable_sort(&mut outcome.assignment, || receipt_checkpoint(cx))?;
    }
    if let CausalReceiptDomain::ModeCell { assignment } = &mut draft.domain {
        cancellable_sort(assignment, || receipt_checkpoint(cx))?;
    }
    cancellable_sort(&mut draft.unknown_axes, || receipt_checkpoint(cx))?;
    cancellable_sort_by_fallible(
        &mut draft.conditional_outcomes,
        |left, right| {
            compare_conditional_outcomes_cancellable(left, right, || receipt_checkpoint(cx))
        },
        || receipt_checkpoint(cx),
    )?;

    let mut findings = Vec::new();
    if !identity_receipt_adjudication_eq(draft.structure, graph.structure_identity_receipt()) {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::GraphIdentityMismatch,
            CausalReceiptSubject::StructureIdentity,
        ));
    }
    if !identity_receipt_adjudication_eq(draft.artifact, graph.artifact_identity_receipt()) {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::GraphIdentityMismatch,
            CausalReceiptSubject::GraphArtifactIdentity,
        ));
    }

    let graph_conditions = receipt_condition_domain(graph, cx)?;
    let is_summary = matches!(draft.domain, CausalReceiptDomain::HybridSummary);
    let domain_valid = match &draft.domain {
        CausalReceiptDomain::UnconditionalGraph => graph_conditions.is_empty(),
        CausalReceiptDomain::ModeCell { assignment } => {
            !graph_conditions.is_empty()
                && assignment_is_valid_cancellable(assignment, &graph_conditions, cx)?
        }
        CausalReceiptDomain::HybridSummary => !graph_conditions.is_empty(),
    };
    if !domain_valid {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::ConditionalCoverageMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }
    let mode_assignment = match &draft.domain {
        CausalReceiptDomain::ModeCell { assignment } => Some(assignment.as_slice()),
        CausalReceiptDomain::UnconditionalGraph | CausalReceiptDomain::HybridSummary => None,
    };
    let mut equation_ids = BTreeSet::new();
    for (index, equation) in graph.equations().iter().enumerate() {
        receipt_poll(cx, index)?;
        if !is_summary
            && equation.solve_participation == EquationParticipation::Matching
            && receipt_row_is_active(&equation.activation, mode_assignment, cx)?
        {
            equation_ids.insert(equation.id.clone());
        }
    }
    let mut active_variables = BTreeSet::new();
    for (index, variable) in graph.variables().iter().enumerate() {
        receipt_poll(cx, index)?;
        if !is_summary && receipt_row_is_active(&variable.activation, mode_assignment, cx)? {
            active_variables.insert(variable.id.clone());
        }
    }
    let mut active_incidences = Vec::new();
    for (index, incidence) in graph.incidences().iter().enumerate() {
        receipt_poll(cx, index)?;
        if !is_summary
            && equation_ids.contains(&incidence.equation)
            && active_variables.contains(&incidence.variable)
            && receipt_row_is_active(&incidence.activation, mode_assignment, cx)?
        {
            active_incidences.push(incidence);
        }
    }
    let mut unknown_vertices = BTreeSet::new();
    let mut states_with_active_unknown_derivatives = BTreeSet::new();
    let mut incidence_keys = BTreeSet::new();
    for (index, incidence) in active_incidences.iter().enumerate() {
        receipt_poll(cx, index)?;
        if incidence.solve_participation == SolveParticipation::Unknown {
            if incidence.derivative_order > 0 {
                states_with_active_unknown_derivatives.insert(incidence.variable.clone());
            }
            unknown_vertices.insert(DerivativeVariableKey {
                variable: incidence.variable.clone(),
                derivative_order: incidence.derivative_order,
            });
            incidence_keys.insert(CausalMatchingPair {
                incidence: incidence.id.clone(),
                equation: incidence.equation.clone(),
                variable: DerivativeVariableKey {
                    variable: incidence.variable.clone(),
                    derivative_order: incidence.derivative_order,
                },
            });
        }
    }
    for (index, variable) in graph.variables().iter().enumerate() {
        receipt_poll(cx, index)?;
        if active_variables.contains(&variable.id)
            && variable.solve_participation == SolveParticipation::Unknown
            && (variable.role != VariableRole::State
                || !states_with_active_unknown_derivatives.contains(&variable.id))
        {
            unknown_vertices.insert(DerivativeVariableKey {
                variable: variable.id.clone(),
                derivative_order: 0,
            });
        }
    }
    let mut matched_equations = BTreeSet::new();
    let mut matched_variables = BTreeSet::new();
    let mut accepted_matching_rows = 0usize;
    for (index, pair) in draft.matching.iter().enumerate() {
        receipt_poll(cx, index)?;
        let endpoints_valid =
            equation_ids.contains(&pair.equation) && unknown_vertices.contains(&pair.variable);
        let incidence_valid = incidence_keys.contains(pair);
        if !endpoints_valid {
            findings.push(CausalReceiptFinding::new(
                CausalReceiptRule::UnknownMatchingEndpoint,
                CausalReceiptSubject::Matching(pair.clone()),
            ));
        }
        if !incidence_valid {
            findings.push(CausalReceiptFinding::new(
                CausalReceiptRule::NonIncidenceMatch,
                CausalReceiptSubject::Matching(pair.clone()),
            ));
        }
        if endpoints_valid && incidence_valid {
            if matched_equations.contains(&pair.equation)
                || matched_variables.contains(&pair.variable)
            {
                findings.push(CausalReceiptFinding::new(
                    CausalReceiptRule::DuplicateMatchingEndpoint,
                    CausalReceiptSubject::Matching(pair.clone()),
                ));
            } else {
                let equation_was_new = matched_equations.insert(pair.equation.clone());
                let variable_was_new = matched_variables.insert(pair.variable.clone());
                debug_assert!(equation_was_new && variable_was_new);
                accepted_matching_rows += 1;
            }
        }
        enforce_receipt_finding_budget(&findings)?;
    }

    let mut expected_unmatched_equations = Vec::new();
    for (index, equation) in equation_ids.difference(&matched_equations).enumerate() {
        receipt_poll(cx, index)?;
        expected_unmatched_equations.push(equation.clone());
    }
    let mut expected_unmatched_variables = Vec::new();
    for (index, variable) in unknown_vertices.difference(&matched_variables).enumerate() {
        receipt_poll(cx, index)?;
        expected_unmatched_variables.push(variable.clone());
    }
    if !cancellable_slice_eq(
        &draft.unmatched_equations,
        &expected_unmatched_equations,
        || receipt_checkpoint(cx),
    )? {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::UnmatchedSetMismatch,
            CausalReceiptSubject::UnmatchedEquations,
        ));
    }
    if !cancellable_slice_eq(
        &draft.unmatched_variables,
        &expected_unmatched_variables,
        || receipt_checkpoint(cx),
    )? {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::UnmatchedSetMismatch,
            CausalReceiptSubject::UnmatchedVariables,
        ));
    }

    let vacuous_empty_graph = equation_ids.is_empty() && unknown_vertices.is_empty();
    let witness_determination = match (
        vacuous_empty_graph,
        expected_unmatched_equations.is_empty(),
        expected_unmatched_variables.is_empty(),
    ) {
        (true, _, _) => DeterminationClass::EmptyProjection,
        (false, true, true) => DeterminationClass::WellDetermined,
        (false, true, false) => DeterminationClass::UnderDetermined,
        (false, false, true) => DeterminationClass::OverDetermined,
        (false, false, false) => DeterminationClass::Mixed,
    };
    let smaller_side = equation_ids.len().min(unknown_vertices.len());
    let saturates_min_side = !vacuous_empty_graph && accepted_matching_rows == smaller_side;
    let evidence_checker = match &draft.evidence {
        CausalReceiptEvidence::Unverified => None,
        CausalReceiptEvidence::CheckerReferenced(checker) => Some(checker),
    };
    let actual_matching_set = if draft.maximum_matching_certificate.is_some() {
        match causal_matching_set_identity_cancellable(&draft.matching, cx) {
            Ok(identity) => Some(identity),
            Err(error) => return Err(identity_receipt_refusal(error)),
        }
    } else {
        None
    };
    let maximum_binding_valid = if let Some(binding) = &draft.maximum_matching_certificate {
        identity_receipt_adjudication_eq(binding.structure, graph.structure_identity_receipt())
            && identity_receipt_adjudication_eq(binding.artifact, graph.artifact_identity_receipt())
            && receipt_domains_equal_cancellable(&binding.domain, &draft.domain, cx)?
            && actual_matching_set.is_some_and(|actual| {
                identity_receipt_adjudication_eq(actual, binding.matching_set)
            })
            && evidence_checker == Some(&binding.checker)
    } else {
        false
    };
    let maximum_supported = saturates_min_side || maximum_binding_valid;
    let expected_rank = if smaller_side == 0 {
        StructuralRankState::NotApplicable
    } else if saturates_min_side {
        StructuralRankState::FullRelativeToMinSide
    } else {
        StructuralRankState::Deficient
    };
    let mut uniform_child_determination = None;
    let mut uniform_child_rank = None;
    let mut child_determination_nonuniform = false;
    let mut child_rank_nonuniform = false;

    let mut graph_condition_ids = BTreeSet::new();
    for (index, condition) in graph_conditions.keys().enumerate() {
        receipt_poll(cx, index)?;
        graph_condition_ids.insert(condition.clone());
    }
    let mut duplicate_outcome = false;
    let mut children_bound_valid = true;
    let mut children_all_concrete = true;
    for (outcome_index, outcome) in draft.conditional_outcomes.iter().enumerate() {
        receipt_poll(cx, outcome_index)?;
        let duplicate_assignment = if outcome_index > 0 {
            cancellable_slice_eq(
                &draft.conditional_outcomes[outcome_index - 1].assignment,
                &outcome.assignment,
                || receipt_checkpoint(cx),
            )?
        } else {
            false
        };
        duplicate_outcome |= duplicate_assignment;
        let mut duplicate_selection = false;
        let mut branches_valid = true;
        let mut selected_conditions = BTreeSet::new();
        for (selection_index, selection) in outcome.assignment.iter().enumerate() {
            receipt_poll(cx, selection_index)?;
            duplicate_selection |= selection_index > 0
                && outcome.assignment[selection_index - 1].condition == selection.condition;
            selected_conditions.insert(selection.condition.clone());
            branches_valid &= graph_conditions
                .get(&selection.condition)
                .is_some_and(|branches| branches.contains(&selection.branch));
        }
        let condition_set_valid =
            cancellable_set_eq(&selected_conditions, &graph_condition_ids, || {
                receipt_checkpoint(cx)
            })?;
        let child_bound_valid = !outcome.assignment.is_empty()
            && identity_receipt_adjudication_eq(
                outcome.structure,
                graph.structure_identity_receipt(),
            )
            && identity_receipt_adjudication_eq(
                outcome.artifact,
                graph.artifact_identity_receipt(),
            )
            && !duplicate_selection
            && condition_set_valid
            && branches_valid;
        children_bound_valid &= child_bound_valid;
        children_all_concrete &= outcome.determination != DeterminationClass::Unknown
            && outcome.structural_rank != StructuralRankState::Unknown;
        if child_bound_valid && !duplicate_assignment {
            if outcome.determination != DeterminationClass::Unknown
                && !child_determination_nonuniform
            {
                match uniform_child_determination {
                    None => uniform_child_determination = Some(outcome.determination),
                    Some(existing) if existing != outcome.determination => {
                        uniform_child_determination = None;
                        child_determination_nonuniform = true;
                    }
                    Some(_) => {}
                }
            }
            if outcome.structural_rank != StructuralRankState::Unknown && !child_rank_nonuniform {
                match uniform_child_rank {
                    None => uniform_child_rank = Some(outcome.structural_rank),
                    Some(existing) if existing != outcome.structural_rank => {
                        uniform_child_rank = None;
                        child_rank_nonuniform = true;
                    }
                    Some(_) => {}
                }
            }
        }
    }
    let determination_valid = if is_summary {
        match draft.conditionality {
            Conditionality::Conditional => {
                draft.determination == DeterminationClass::Unknown
                    || uniform_child_determination == Some(draft.determination)
            }
            Conditionality::Unconditional => true,
            Conditionality::Unknown => draft.determination == DeterminationClass::Unknown,
        }
    } else if vacuous_empty_graph {
        matches!(
            draft.determination,
            DeterminationClass::EmptyProjection | DeterminationClass::Unknown
        )
    } else if maximum_supported {
        matches!(draft.determination, DeterminationClass::Unknown)
            || draft.determination == witness_determination
    } else {
        draft.determination == DeterminationClass::Unknown
    };
    let rank_valid = if is_summary {
        match draft.conditionality {
            Conditionality::Conditional => {
                draft.structural_rank == StructuralRankState::Unknown
                    || uniform_child_rank == Some(draft.structural_rank)
            }
            Conditionality::Unconditional => true,
            Conditionality::Unknown => draft.structural_rank == StructuralRankState::Unknown,
        }
    } else if maximum_supported || smaller_side == 0 {
        matches!(draft.structural_rank, StructuralRankState::Unknown)
            || draft.structural_rank == expected_rank
    } else {
        draft.structural_rank == StructuralRankState::Unknown
    };
    if !determination_valid
        || !rank_valid
        || !causal_axes_compatible(draft.determination, draft.structural_rank)
    {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::OutcomeAxisMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }
    if draft.maximum_matching_certificate.is_some() && !maximum_binding_valid {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::OutcomeAxisMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }
    if is_summary && draft.maximum_matching_certificate.is_some() {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::OutcomeAxisMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }
    let actual_conditional_outcome_set = match draft
        .conditional_coverage
        .as_ref()
        .map(|binding| binding.claim)
    {
        Some(ConditionalCoverageClaim::ModeCells(_)) => {
            match conditional_outcome_set_identity_cancellable(&draft.conditional_outcomes, cx) {
                Ok(identity) => Some(identity),
                Err(error) => return Err(identity_receipt_refusal(error)),
            }
        }
        Some(ConditionalCoverageClaim::UniformTheorem { .. }) | None => None,
    };
    let coverage_binding_valid = draft.conditional_coverage.as_ref().is_some_and(|binding| {
        if !identity_receipt_adjudication_eq(binding.structure, graph.structure_identity_receipt())
            || !identity_receipt_adjudication_eq(
                binding.artifact,
                graph.artifact_identity_receipt(),
            )
        {
            return false;
        }
        let claim_matches = match (&draft.domain, draft.conditionality, binding.claim) {
            (
                CausalReceiptDomain::HybridSummary,
                Conditionality::Conditional,
                ConditionalCoverageClaim::ModeCells(bound_outcomes),
            ) => actual_conditional_outcome_set
                .is_some_and(|actual| identity_receipt_adjudication_eq(actual, bound_outcomes)),
            (
                CausalReceiptDomain::HybridSummary,
                Conditionality::Unconditional,
                ConditionalCoverageClaim::UniformTheorem {
                    determination,
                    structural_rank,
                },
            ) => determination == draft.determination && structural_rank == draft.structural_rank,
            _ => false,
        };
        claim_matches && evidence_checker == Some(&binding.checker)
    });
    let conditional_valid = match (&draft.domain, draft.conditionality) {
        (
            CausalReceiptDomain::UnconditionalGraph | CausalReceiptDomain::ModeCell { .. },
            Conditionality::Unconditional,
        ) => draft.conditional_outcomes.is_empty() && draft.conditional_coverage.is_none(),
        (CausalReceiptDomain::HybridSummary, Conditionality::Conditional) => {
            !draft.conditional_outcomes.is_empty()
                && !duplicate_outcome
                && children_bound_valid
                && children_all_concrete
                && coverage_binding_valid
        }
        (CausalReceiptDomain::HybridSummary, Conditionality::Unconditional) => {
            draft.conditional_outcomes.is_empty() && coverage_binding_valid
        }
        (CausalReceiptDomain::HybridSummary, Conditionality::Unknown) => {
            !duplicate_outcome && children_bound_valid && draft.conditional_coverage.is_none()
        }
        _ => false,
    };
    if !conditional_valid {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::ConditionalCoverageMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }

    let expected_unknown_axes: BTreeSet<_> = [
        (draft.determination == DeterminationClass::Unknown)
            .then_some(CausalOutcomeAxis::Determination),
        (draft.structural_rank == StructuralRankState::Unknown)
            .then_some(CausalOutcomeAxis::StructuralRank),
        (draft.conditionality == Conditionality::Unknown)
            .then_some(CausalOutcomeAxis::Conditionality),
    ]
    .into_iter()
    .flatten()
    .collect();
    let supplied_unknown_axes: BTreeSet<_> =
        draft.unknown_axes.iter().map(|state| state.axis).collect();
    let duplicate_unknown_axis = draft
        .unknown_axes
        .windows(2)
        .any(|pair| pair[0].axis == pair[1].axis);
    let checkpoint_mismatch = draft.unknown_axes.iter().any(|state| {
        let checkpoint_required = matches!(
            state.reason,
            CausalUnknownReason::Cancelled | CausalUnknownReason::BudgetExhausted
        );
        checkpoint_required != state.resume_checkpoint.is_some()
    });
    let nonuniform_reason_mismatch = draft.unknown_axes.iter().any(|state| {
        let children_are_nonuniform = match state.axis {
            CausalOutcomeAxis::Determination => child_determination_nonuniform,
            CausalOutcomeAxis::StructuralRank => child_rank_nonuniform,
            CausalOutcomeAxis::Conditionality => false,
        };
        (children_are_nonuniform && state.reason != CausalUnknownReason::NonUniformAcrossModes)
            || (state.reason == CausalUnknownReason::NonUniformAcrossModes
                && (!is_summary || !children_are_nonuniform))
    });
    if duplicate_unknown_axis
        || checkpoint_mismatch
        || nonuniform_reason_mismatch
        || supplied_unknown_axes != expected_unknown_axes
    {
        findings.push(CausalReceiptFinding::new(
            CausalReceiptRule::UnknownReasonMismatch,
            CausalReceiptSubject::Receipt,
        ));
    }

    enforce_receipt_finding_budget(&findings)?;
    if !findings.is_empty() {
        return Err(receipt_refusal_cancellable(findings, cx)?);
    }
    receipt_checkpoint(cx)?;
    let structure_receipt = graph.structure_identity_receipt();
    let artifact_receipt = graph.artifact_identity_receipt();
    let outcome_receipt = match causal_outcome_identity(&draft, structure_receipt, cx) {
        Ok(receipt) => receipt,
        Err(error) => return Err(identity_receipt_refusal(error)),
    };
    let receipt = match causalization_receipt_identity(
        &draft,
        structure_receipt,
        artifact_receipt,
        outcome_receipt,
        cx,
    ) {
        Ok(receipt) => receipt,
        Err(error) => return Err(identity_receipt_refusal(error)),
    };
    receipt_checkpoint(cx)?;
    Ok(AdmittedCausalizationReceipt {
        outcome_receipt,
        structure: structure_receipt,
        artifact: artifact_receipt,
        analysis: draft.analysis,
        domain: draft.domain,
        determination: draft.determination,
        structural_rank: draft.structural_rank,
        conditionality: draft.conditionality,
        matching: draft.matching,
        unmatched_equations: draft.unmatched_equations,
        unmatched_variables: draft.unmatched_variables,
        conditional_outcomes: draft.conditional_outcomes,
        maximum_matching_certificate: draft.maximum_matching_certificate,
        conditional_coverage: draft.conditional_coverage,
        unknown_axes: draft.unknown_axes,
        evidence: draft.evidence,
        receipt,
    })
}

fn causal_outcome_identity(
    draft: &CausalizationReceiptDraft,
    structure: IdentityReceipt<CausalStructureIdV1>,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalOutcomeIdV1>, CausalIdentityError> {
    let domain = receipt_domain_row_cancellable(&draft.domain, cx)?;
    let axes = [
        determination_tag(draft.determination),
        structural_rank_tag(draft.structural_rank),
        conditionality_tag(draft.conditionality),
    ];
    let mut structure_row = Vec::with_capacity(116);
    push_identity_receipt_adjudication(&mut structure_row, structure);
    let encoder =
        CanonicalEncoder::<CausalOutcomeIdV1, _>::new(CAUSAL_RECEIPT_IDENTITY_LIMITS, || {
            cx.checkpoint().is_err()
        })?
        .u64(
            Field::new(0, "causal-outcome-schema-version"),
            u64::from(CAUSAL_OUTCOME_IDENTITY_SCHEMA_VERSION_V1),
        )?
        .child(Field::new(1, "causal-structure-id"), structure.id())?
        .bytes(
            Field::new(2, "causal-structure-receipt-adjudication"),
            &structure_row,
        )?
        .bytes(Field::new(3, "analysis-domain"), &domain)?
        .bytes(Field::new(4, "outcome-axes"), &axes)?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(5, "matching"),
        &draft.matching,
        cx,
        |_, _| Ok(3 * IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2),
        |pair, _, row| {
            matching_row(pair, row);
            Ok(())
        },
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(6, "unmatched-equations"),
        &draft.unmatched_equations,
        cx,
        |_, _| Ok(IDENTITY_RECEIPT_ADJUDICATION_BYTES),
        |equation, _, row| {
            push_identity_receipt_adjudication(row, equation.identity_receipt());
            Ok(())
        },
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(7, "unmatched-variables"),
        &draft.unmatched_variables,
        cx,
        |_, _| Ok(IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2),
        |variable, _, row| {
            derivative_variable_row(variable, row);
            Ok(())
        },
    )?;
    stream_identity_rows(
        encoder,
        Field::new(8, "conditional-outcome-semantics"),
        &draft.conditional_outcomes,
        cx,
        normalized_conditional_outcome_row_len_cancellable,
        normalized_conditional_outcome_row_cancellable,
    )?
    .finish()
    .map_err(Into::into)
}

fn causalization_receipt_identity(
    draft: &CausalizationReceiptDraft,
    structure: IdentityReceipt<CausalStructureIdV1>,
    artifact: IdentityReceipt<CausalGraphArtifactIdV1>,
    outcome: IdentityReceipt<CausalOutcomeIdV1>,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalizationReceiptIdV1>, CausalIdentityError> {
    let axes = [
        determination_tag(draft.determination),
        structural_rank_tag(draft.structural_rank),
        conditionality_tag(draft.conditionality),
    ];
    let analysis = analysis_context_row(&draft.analysis);
    let domain = receipt_domain_row_cancellable(&draft.domain, cx)?;
    let maximum_matching_certificate = match draft.maximum_matching_certificate.as_ref() {
        Some(binding) => Some(maximum_matching_binding_row_cancellable(binding, cx)?),
        None => None,
    };
    let conditional_coverage = draft
        .conditional_coverage
        .as_ref()
        .map(conditional_coverage_row);
    let (evidence_tag, evidence_payload) = match &draft.evidence {
        CausalReceiptEvidence::Unverified => (1, Vec::new()),
        CausalReceiptEvidence::CheckerReferenced(checker) => {
            let mut row = Vec::with_capacity(128);
            checker.append_canonical(&mut row);
            (2, row)
        }
    };
    let mut structure_row = Vec::with_capacity(116);
    push_identity_receipt_adjudication(&mut structure_row, structure);
    let mut artifact_row = Vec::with_capacity(116);
    push_identity_receipt_adjudication(&mut artifact_row, artifact);
    let mut outcome_row = Vec::with_capacity(116);
    push_identity_receipt_adjudication(&mut outcome_row, outcome);

    let encoder = CanonicalEncoder::<CausalizationReceiptIdV1, _>::new(
        CAUSAL_RECEIPT_IDENTITY_LIMITS,
        || cx.checkpoint().is_err(),
    )?
    .u64(
        Field::new(0, "causalization-receipt-schema-version"),
        u64::from(CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1),
    )?
    .child(Field::new(1, "causal-structure-id"), structure.id())?
    .bytes(
        Field::new(2, "causal-structure-receipt-adjudication"),
        &structure_row,
    )?
    .child(Field::new(3, "causal-graph-artifact-id"), artifact.id())?
    .bytes(
        Field::new(4, "causal-graph-artifact-receipt-adjudication"),
        &artifact_row,
    )?
    .bytes(Field::new(5, "analysis-context"), &analysis)?
    .bytes(Field::new(6, "analysis-domain"), &domain)?
    .bytes(Field::new(7, "outcome-axes"), &axes)?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(8, "matching"),
        &draft.matching,
        cx,
        |_, _| Ok(3 * IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2),
        |pair, _, row| {
            matching_row(pair, row);
            Ok(())
        },
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(9, "unmatched-equations"),
        &draft.unmatched_equations,
        cx,
        |_, _| Ok(IDENTITY_RECEIPT_ADJUDICATION_BYTES),
        |equation, _, row| {
            push_identity_receipt_adjudication(row, equation.identity_receipt());
            Ok(())
        },
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(10, "unmatched-variables"),
        &draft.unmatched_variables,
        cx,
        |_, _| Ok(IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2),
        |variable, _, row| {
            derivative_variable_row(variable, row);
            Ok(())
        },
    )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(11, "conditional-outcomes"),
        &draft.conditional_outcomes,
        cx,
        conditional_outcome_row_len_cancellable,
        conditional_outcome_row_cancellable,
    )?;
    let encoder = encoder
        .optional_bytes(
            Field::new(12, "maximum-matching-certificate"),
            maximum_matching_certificate.as_deref(),
        )?
        .optional_bytes(
            Field::new(13, "conditional-coverage"),
            conditional_coverage.as_deref(),
        )?;
    let encoder = stream_identity_rows(
        encoder,
        Field::new(14, "unknown-axes"),
        &draft.unknown_axes,
        cx,
        |state, _| unknown_axis_row_len(state),
        |state, _, row| {
            unknown_axis_row(state, row);
            Ok(())
        },
    )?;
    encoder
        .variant(
            Field::new(15, "evidence-state"),
            evidence_tag,
            &evidence_payload,
        )?
        .child(Field::new(16, "normalized-causal-outcome-id"), outcome.id())?
        .bytes(
            Field::new(17, "normalized-causal-outcome-receipt-adjudication"),
            &outcome_row,
        )?
        .finish()
        .map_err(Into::into)
}

fn receipt_domain_row_cancellable(
    domain: &CausalReceiptDomain,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CanonicalError> {
    let assignment_bytes = match domain {
        CausalReceiptDomain::ModeCell { assignment } => {
            assignment_canonical_len_cancellable(assignment, cx)?
        }
        CausalReceiptDomain::UnconditionalGraph | CausalReceiptDomain::HybridSummary => 0,
    };
    let mut out = Vec::with_capacity(checked_row_len_add(1, assignment_bytes)?);
    match domain {
        CausalReceiptDomain::UnconditionalGraph => out.push(1),
        CausalReceiptDomain::ModeCell { assignment } => {
            out.push(2);
            push_row_count(&mut out, assignment.len())?;
            for (index, selection) in assignment.iter().enumerate() {
                identity_materialization_poll(cx, index, out.len())?;
                selection.condition.append_canonical(&mut out);
                selection.branch.append_canonical(&mut out);
            }
        }
        CausalReceiptDomain::HybridSummary => out.push(3),
    }
    debug_assert_eq!(
        out.len(),
        1 + if assignment_bytes == 0 {
            0
        } else {
            assignment_bytes
        }
    );
    identity_materialization_checkpoint(cx, out.len())?;
    Ok(out)
}

fn unknown_axis_row(state: &CausalUnknownAxisState, out: &mut CausalIdentityRowWriter) {
    debug_assert!(out.is_empty());
    out.push(outcome_axis_tag(state.axis));
    out.push(unknown_reason_tag(state.reason));
    push_optional_ref(out, state.resume_checkpoint.as_ref());
}

fn unknown_axis_row_len(state: &CausalUnknownAxisState) -> Result<usize, CanonicalError> {
    checked_row_len_add(
        2,
        optional_causal_ref_canonical_len(
            state
                .resume_checkpoint
                .as_ref()
                .map(CausalCheckpointRef::namespace),
        )?,
    )
}

fn analysis_context_row(context: &CausalAnalysisContext) -> Vec<u8> {
    let mut out = Vec::with_capacity(512);
    context.analyzer.append_canonical(&mut out);
    context.budget.append_canonical(&mut out);
    context.capabilities.append_canonical(&mut out);
    push_seed_policy(&mut out, context.seed_policy);
    out.push(determinism_tag(context.determinism));
    out
}

fn matching_row(pair: &CausalMatchingPair, out: &mut CausalIdentityRowWriter) {
    debug_assert!(out.is_empty());
    push_identity_receipt_adjudication(out, pair.incidence.identity_receipt());
    push_identity_receipt_adjudication(out, pair.equation.identity_receipt());
    push_identity_receipt_adjudication(out, pair.variable.variable.identity_receipt());
    out.extend_from_slice(&pair.variable.derivative_order.to_le_bytes());
}

fn causal_matching_pair_nominal_eq(left: &CausalMatchingPair, right: &CausalMatchingPair) -> bool {
    left.incidence.identity() == right.incidence.identity()
        && left.equation.identity() == right.equation.identity()
        && left.variable.variable.identity() == right.variable.variable.identity()
        && left.variable.derivative_order == right.variable.derivative_order
}

fn compare_causal_matching_pairs_nominal(
    left: &CausalMatchingPair,
    right: &CausalMatchingPair,
) -> core::cmp::Ordering {
    left.incidence
        .identity()
        .cmp(&right.incidence.identity())
        .then_with(|| left.equation.identity().cmp(&right.equation.identity()))
        .then_with(|| {
            left.variable
                .variable
                .identity()
                .cmp(&right.variable.variable.identity())
        })
        .then_with(|| {
            left.variable
                .derivative_order
                .cmp(&right.variable.derivative_order)
        })
}

fn causal_matching_set_identity_cancellable(
    canonical_matching: &[CausalMatchingPair],
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalMatchingSetIdV1>, CausalIdentityError> {
    let encoder =
        CanonicalEncoder::<CausalMatchingSetIdV1, _>::new(CAUSAL_RECEIPT_IDENTITY_LIMITS, || {
            cx.checkpoint().is_err()
        })?;
    stream_identity_rows(
        encoder,
        Field::new(0, "matching-pairs"),
        canonical_matching,
        cx,
        |_, _| Ok(3 * IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2),
        |pair, _, row| {
            matching_row(pair, row);
            Ok(())
        },
    )?
    .finish()
    .map_err(Into::into)
}

fn maximum_matching_binding_row_cancellable(
    binding: &MaximumMatchingBinding,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CanonicalError> {
    let domain = receipt_domain_row_cancellable(&binding.domain, cx)?;
    let mut out = Vec::with_capacity(1_024usize.saturating_add(domain.len()));
    push_identity_receipt_adjudication(&mut out, binding.structure);
    push_identity_receipt_adjudication(&mut out, binding.artifact);
    push_len_prefixed(&mut out, &domain);
    out.extend_from_slice(binding.matching_set.id().as_bytes());
    out.extend_from_slice(binding.matching_set.canonical_preimage().as_bytes());
    out.extend_from_slice(binding.matching_set.schema_id().as_bytes());
    out.extend_from_slice(&binding.matching_set.canonical_bytes().to_le_bytes());
    out.extend_from_slice(&binding.matching_set.field_count().to_le_bytes());
    out.extend_from_slice(&binding.matching_set.collection_items().to_le_bytes());
    binding.certificate.append_canonical(&mut out);
    binding.checker.append_canonical(&mut out);
    identity_materialization_checkpoint(cx, out.len())?;
    Ok(out)
}

fn derivative_variable_row(variable: &DerivativeVariableKey, out: &mut CausalIdentityRowWriter) {
    debug_assert!(out.is_empty());
    push_identity_receipt_adjudication(out, variable.variable.identity_receipt());
    out.extend_from_slice(&variable.derivative_order.to_le_bytes());
}

fn conditional_outcome_row_cancellable(
    outcome: &ConditionalCausalOutcome,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    // The full child receipt below already commits its ordered unknown-axis
    // reasons and checkpoints. Keep those progress coordinates out of the
    // normalized child outcome while retaining them transitively in this
    // provenance-bearing parent row and directly on the typed child object.
    debug_assert!(out.is_empty());
    push_row_count(out, outcome.assignment.len())?;
    for (index, selection) in outcome.assignment.iter().enumerate() {
        identity_materialization_poll(cx, index, out.len())?;
        selection.condition.append_canonical(out);
        selection.branch.append_canonical(out);
    }
    out.push(determination_tag(outcome.determination));
    out.push(structural_rank_tag(outcome.structural_rank));
    push_identity_receipt_adjudication(out, outcome.outcome);
    push_identity_receipt_adjudication(out, outcome.receipt);
    identity_materialization_checkpoint(cx, out.len())
}

fn conditional_outcome_row_len_cancellable(
    outcome: &ConditionalCausalOutcome,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        assignment_canonical_len_cancellable(&outcome.assignment, cx)?,
        2,
        2 * IDENTITY_RECEIPT_ADJUDICATION_BYTES,
    ])
}

fn normalized_conditional_outcome_row_cancellable(
    outcome: &ConditionalCausalOutcome,
    cx: &Cx<'_>,
    out: &mut CausalIdentityRowWriter,
) -> Result<(), CanonicalError> {
    debug_assert!(out.is_empty());
    push_row_count(out, outcome.assignment.len())?;
    for (index, selection) in outcome.assignment.iter().enumerate() {
        identity_materialization_poll(cx, index, out.len())?;
        selection.condition.append_canonical(out);
        selection.branch.append_canonical(out);
    }
    out.push(determination_tag(outcome.determination));
    out.push(structural_rank_tag(outcome.structural_rank));
    push_identity_receipt_adjudication(out, outcome.outcome);
    identity_materialization_checkpoint(cx, out.len())
}

fn normalized_conditional_outcome_row_len_cancellable(
    outcome: &ConditionalCausalOutcome,
    cx: &Cx<'_>,
) -> Result<usize, CanonicalError> {
    checked_row_len_sum([
        assignment_canonical_len_cancellable(&outcome.assignment, cx)?,
        2,
        IDENTITY_RECEIPT_ADJUDICATION_BYTES,
    ])
}

fn conditional_outcome_set_identity_cancellable(
    canonical_outcomes: &[ConditionalCausalOutcome],
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<ConditionalOutcomeSetIdV1>, CausalIdentityError> {
    let encoder = CanonicalEncoder::<ConditionalOutcomeSetIdV1, _>::new(
        CAUSAL_RECEIPT_IDENTITY_LIMITS,
        || cx.checkpoint().is_err(),
    )?;
    stream_identity_rows(
        encoder,
        Field::new(0, "mode-cell-outcomes"),
        canonical_outcomes,
        cx,
        conditional_outcome_row_len_cancellable,
        conditional_outcome_row_cancellable,
    )?
    .finish()
    .map_err(Into::into)
}

fn conditional_coverage_row(binding: &ConditionalCoverageBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(560);
    push_identity_receipt_adjudication(&mut out, binding.structure);
    push_identity_receipt_adjudication(&mut out, binding.artifact);
    match binding.claim {
        ConditionalCoverageClaim::ModeCells(outcome_set) => {
            out.push(1);
            out.extend_from_slice(outcome_set.id().as_bytes());
            out.extend_from_slice(outcome_set.canonical_preimage().as_bytes());
            out.extend_from_slice(outcome_set.schema_id().as_bytes());
            out.extend_from_slice(&outcome_set.canonical_bytes().to_le_bytes());
            out.extend_from_slice(&outcome_set.field_count().to_le_bytes());
            out.extend_from_slice(&outcome_set.collection_items().to_le_bytes());
        }
        ConditionalCoverageClaim::UniformTheorem {
            determination,
            structural_rank,
        } => {
            out.push(2);
            out.push(determination_tag(determination));
            out.push(structural_rank_tag(structural_rank));
        }
    }
    binding.certificate.append_canonical(&mut out);
    binding.checker.append_canonical(&mut out);
    out
}

const fn determination_tag(state: DeterminationClass) -> u8 {
    match state {
        DeterminationClass::WellDetermined => 1,
        DeterminationClass::UnderDetermined => 2,
        DeterminationClass::OverDetermined => 3,
        DeterminationClass::Mixed => 4,
        DeterminationClass::Unknown => 5,
        DeterminationClass::EmptyProjection => 6,
    }
}

const fn structural_rank_tag(state: StructuralRankState) -> u8 {
    match state {
        StructuralRankState::FullRelativeToMinSide => 1,
        StructuralRankState::Deficient => 2,
        StructuralRankState::NotApplicable => 3,
        StructuralRankState::Unknown => 4,
    }
}

const fn conditionality_tag(state: Conditionality) -> u8 {
    match state {
        Conditionality::Unconditional => 1,
        Conditionality::Conditional => 2,
        Conditionality::Unknown => 3,
    }
}

const fn unknown_reason_tag(reason: CausalUnknownReason) -> u8 {
    match reason {
        CausalUnknownReason::NotAnalyzed => 1,
        CausalUnknownReason::Cancelled => 2,
        CausalUnknownReason::UnsupportedStructure => 3,
        CausalUnknownReason::BudgetExhausted => 4,
        CausalUnknownReason::IncompleteMetadata => 5,
        CausalUnknownReason::NonUniformAcrossModes => 6,
    }
}

const fn outcome_axis_tag(axis: CausalOutcomeAxis) -> u8 {
    match axis {
        CausalOutcomeAxis::Determination => 1,
        CausalOutcomeAxis::StructuralRank => 2,
        CausalOutcomeAxis::Conditionality => 3,
    }
}

/// Causal artifact family migrated by an explicit migration receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CausalMigrationArtifactKind {
    /// Producer-independent normalized structure.
    Structure,
    /// Provenance-bearing graph artifact.
    GraphArtifact,
    /// Causalization outcome receipt.
    CausalizationReceipt,
}

/// Exact bounded identity receipt retained from a predecessor schema.
///
/// Unlike a bare digest, this keeps the predecessor schema descriptor root and
/// complete canonical-preimage root so future readers can adjudicate whether a
/// migration actually consumed the recorded historical artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HistoricalCausalIdentityReceipt {
    /// Historical artifact family.
    artifact_kind: CausalMigrationArtifactKind,
    /// Historical schema version; zero is an explicit pre-canonical legacy lane.
    schema_version: u32,
    /// Historical typed semantic/evidence identity bytes.
    semantic_identity: [u8; 32],
    /// Root of the complete historical canonical preimage.
    canonical_preimage: [u8; 32],
    /// Historical schema descriptor root.
    schema_identity: [u8; 32],
    /// Historical canonical frame byte count.
    canonical_bytes: u64,
    /// Historical top-level field count.
    field_count: u32,
    /// Historical encoded collection-item count.
    collection_items: u64,
}

/// Refusal from constructing a migration draft around an incomplete historical
/// receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistoricalReceiptError {
    /// One mandatory digest was all zero.
    ZeroDigest,
    /// A canonical predecessor omitted mandatory frame/field metadata.
    IncompleteCanonicalMetadata,
}

impl fmt::Display for HistoricalReceiptError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDigest => f.write_str(
                "historical semantic, canonical-preimage, and schema digests must be nonzero",
            ),
            Self::IncompleteCanonicalMetadata => f.write_str(
                "canonical historical schemas require nonzero frame bytes and field count",
            ),
        }
    }
}

impl std::error::Error for HistoricalReceiptError {}

impl HistoricalCausalIdentityReceipt {
    /// Construct an exact historical identity receipt.
    ///
    /// # Errors
    /// Refuses any all-zero mandatory digest, or incomplete canonical frame
    /// metadata for a nonzero historical schema version.
    pub fn new(
        artifact_kind: CausalMigrationArtifactKind,
        schema_version: u32,
        semantic_identity: [u8; 32],
        canonical_preimage: [u8; 32],
        schema_identity: [u8; 32],
        canonical_bytes: u64,
        field_count: u32,
        collection_items: u64,
    ) -> Result<Self, HistoricalReceiptError> {
        if semantic_identity == [0; 32]
            || canonical_preimage == [0; 32]
            || schema_identity == [0; 32]
        {
            return Err(HistoricalReceiptError::ZeroDigest);
        }
        if schema_version != 0 && (canonical_bytes == 0 || field_count == 0) {
            return Err(HistoricalReceiptError::IncompleteCanonicalMetadata);
        }
        Ok(Self {
            artifact_kind,
            schema_version,
            semantic_identity,
            canonical_preimage,
            schema_identity,
            canonical_bytes,
            field_count,
            collection_items,
        })
    }

    /// Historical artifact family.
    #[must_use]
    pub const fn artifact_kind(&self) -> CausalMigrationArtifactKind {
        self.artifact_kind
    }

    /// Historical schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Historical typed identity bytes.
    #[must_use]
    pub const fn semantic_identity(&self) -> [u8; 32] {
        self.semantic_identity
    }

    /// Historical canonical-preimage root.
    #[must_use]
    pub const fn canonical_preimage(&self) -> [u8; 32] {
        self.canonical_preimage
    }

    /// Historical schema-descriptor root.
    #[must_use]
    pub const fn schema_identity(&self) -> [u8; 32] {
        self.schema_identity
    }

    /// Historical canonical-frame byte count.
    #[must_use]
    pub const fn canonical_bytes(&self) -> u64 {
        self.canonical_bytes
    }

    /// Historical top-level field count.
    #[must_use]
    pub const fn field_count(&self) -> u32 {
        self.field_count
    }

    /// Historical encoded collection-item count.
    #[must_use]
    pub const fn collection_items(&self) -> u64 {
        self.collection_items
    }
}

#[derive(Debug, Clone, Copy)]
enum NativeCausalTargetReceipt {
    Structure(IdentityReceipt<CausalStructureIdV1>),
    GraphArtifact(IdentityReceipt<CausalGraphArtifactIdV1>),
    Causalization(IdentityReceipt<CausalizationReceiptIdV1>),
}

impl PartialEq for NativeCausalTargetReceipt {
    fn eq(&self, other: &Self) -> bool {
        match (*self, *other) {
            (Self::Structure(left), Self::Structure(right)) => {
                identity_receipt_adjudication_eq(left, right)
            }
            (Self::GraphArtifact(left), Self::GraphArtifact(right)) => {
                identity_receipt_adjudication_eq(left, right)
            }
            (Self::Causalization(left), Self::Causalization(right)) => {
                identity_receipt_adjudication_eq(left, right)
            }
            (Self::Structure(_), Self::GraphArtifact(_) | Self::Causalization(_))
            | (Self::GraphArtifact(_), Self::Structure(_) | Self::Causalization(_))
            | (Self::Causalization(_), Self::Structure(_) | Self::GraphArtifact(_)) => false,
        }
    }
}

impl Eq for NativeCausalTargetReceipt {}

impl NativeCausalTargetReceipt {
    const fn kind(self) -> CausalMigrationArtifactKind {
        match self {
            Self::Structure(_) => CausalMigrationArtifactKind::Structure,
            Self::GraphArtifact(_) => CausalMigrationArtifactKind::GraphArtifact,
            Self::Causalization(_) => CausalMigrationArtifactKind::CausalizationReceipt,
        }
    }

    const fn schema_version(self) -> u32 {
        match self {
            Self::Structure(_) => CAUSAL_STRUCTURE_IDENTITY_SCHEMA_VERSION_V1,
            Self::GraphArtifact(_) => CAUSAL_GRAPH_ARTIFACT_IDENTITY_SCHEMA_VERSION_V1,
            Self::Causalization(_) => CAUSALIZATION_RECEIPT_IDENTITY_SCHEMA_VERSION_V1,
        }
    }

    fn identity(self) -> [u8; 32] {
        match self {
            Self::Structure(receipt) => *receipt.id().as_bytes(),
            Self::GraphArtifact(receipt) => *receipt.id().as_bytes(),
            Self::Causalization(receipt) => *receipt.id().as_bytes(),
        }
    }

    fn canonical_preimage(self) -> [u8; 32] {
        match self {
            Self::Structure(receipt) => *receipt.canonical_preimage().as_bytes(),
            Self::GraphArtifact(receipt) => *receipt.canonical_preimage().as_bytes(),
            Self::Causalization(receipt) => *receipt.canonical_preimage().as_bytes(),
        }
    }

    fn schema_identity(self) -> [u8; 32] {
        match self {
            Self::Structure(receipt) => *receipt.schema_id().as_bytes(),
            Self::GraphArtifact(receipt) => *receipt.schema_id().as_bytes(),
            Self::Causalization(receipt) => *receipt.schema_id().as_bytes(),
        }
    }

    const fn canonical_bytes(self) -> u64 {
        match self {
            Self::Structure(receipt) => receipt.canonical_bytes(),
            Self::GraphArtifact(receipt) => receipt.canonical_bytes(),
            Self::Causalization(receipt) => receipt.canonical_bytes(),
        }
    }

    const fn field_count(self) -> u32 {
        match self {
            Self::Structure(receipt) => receipt.field_count(),
            Self::GraphArtifact(receipt) => receipt.field_count(),
            Self::Causalization(receipt) => receipt.field_count(),
        }
    }

    const fn collection_items(self) -> u64 {
        match self {
            Self::Structure(receipt) => receipt.collection_items(),
            Self::GraphArtifact(receipt) => receipt.collection_items(),
            Self::Causalization(receipt) => receipt.collection_items(),
        }
    }
}

/// Explicit migration draft. Migration lineage is deliberately separate from
/// normalized causal structure identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CausalSchemaMigrationDraft {
    /// Complete predecessor identity receipt.
    predecessor: HistoricalCausalIdentityReceipt,
    /// Exact native typed target receipt. This private enum makes cross-family
    /// digest substitution unrepresentable through the public API.
    target: NativeCausalTargetReceipt,
    /// Audited migration implementation/receipt identity.
    migration: CausalMigrationRef,
}

/// Canonical schema marker for migration receipts.
pub enum CausalSchemaMigrationIdentitySchemaV1 {}

impl CanonicalSchema for CausalSchemaMigrationIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.causal-schema-migration.v1";
    const NAME: &'static str = "causal-schema-migration-receipt";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "exact predecessor identity/preimage/schema receipt bound to one native target without changing target semantics";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("migration-receipt-schema-version", WireType::U64),
        FieldSpec::required("artifact-kind", WireType::Variant),
        FieldSpec::required("predecessor-receipt", WireType::Bytes),
        FieldSpec::required("target-receipt", WireType::Bytes),
        FieldSpec::required("migration", WireType::Bytes),
    ];
}

/// Strong identity of one schema migration receipt.
pub type CausalSchemaMigrationIdV1 = EvidenceNodeId<CausalSchemaMigrationIdentitySchemaV1>;

/// Admitted migration receipt retaining both sides exactly.
#[derive(Debug, Clone)]
pub struct AdmittedCausalSchemaMigration {
    kind: CausalMigrationArtifactKind,
    predecessor: HistoricalCausalIdentityReceipt,
    target_identity: [u8; 32],
    target_canonical_preimage: [u8; 32],
    target_schema_identity: [u8; 32],
    target_canonical_bytes: u64,
    target_field_count: u32,
    target_collection_items: u64,
    migration: CausalMigrationRef,
    receipt: IdentityReceipt<CausalSchemaMigrationIdV1>,
}

impl PartialEq for AdmittedCausalSchemaMigration {
    fn eq(&self, other: &Self) -> bool {
        self.kind == other.kind
            && self.predecessor == other.predecessor
            && self.target_identity == other.target_identity
            && self.target_canonical_preimage == other.target_canonical_preimage
            && self.target_schema_identity == other.target_schema_identity
            && self.target_canonical_bytes == other.target_canonical_bytes
            && self.target_field_count == other.target_field_count
            && self.target_collection_items == other.target_collection_items
            && self.migration == other.migration
            && identity_receipt_adjudication_eq(self.receipt, other.receipt)
    }
}

impl Eq for AdmittedCausalSchemaMigration {}

impl AdmittedCausalSchemaMigration {
    /// Migration receipt identity.
    #[must_use]
    pub const fn identity(&self) -> CausalSchemaMigrationIdV1 {
        self.receipt.id()
    }

    /// Complete identity receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<CausalSchemaMigrationIdV1> {
        self.receipt
    }

    /// Migrated artifact family.
    #[must_use]
    pub const fn kind(&self) -> CausalMigrationArtifactKind {
        self.kind
    }

    /// Complete predecessor receipt.
    #[must_use]
    pub const fn predecessor(&self) -> &HistoricalCausalIdentityReceipt {
        &self.predecessor
    }

    /// Exact target identity bytes.
    #[must_use]
    pub const fn target_identity(&self) -> [u8; 32] {
        self.target_identity
    }

    /// Exact target canonical-preimage root.
    #[must_use]
    pub const fn target_canonical_preimage(&self) -> [u8; 32] {
        self.target_canonical_preimage
    }

    /// Exact target schema descriptor root.
    #[must_use]
    pub const fn target_schema_identity(&self) -> [u8; 32] {
        self.target_schema_identity
    }

    /// Exact target canonical-frame byte count.
    #[must_use]
    pub const fn target_canonical_bytes(&self) -> u64 {
        self.target_canonical_bytes
    }

    /// Exact target top-level field count.
    #[must_use]
    pub const fn target_field_count(&self) -> u32 {
        self.target_field_count
    }

    /// Exact target encoded collection-item count.
    #[must_use]
    pub const fn target_collection_items(&self) -> u64 {
        self.target_collection_items
    }

    /// Audited migration implementation/receipt reference.
    #[must_use]
    pub const fn migration(&self) -> &CausalMigrationRef {
        &self.migration
    }
}

/// Migration admission failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CausalMigrationError {
    /// Historical receipt family did not equal the native target family.
    ArtifactKindMismatch,
    /// Predecessor is not strictly older than the target schema.
    PredecessorNotOlder {
        /// Historical version.
        predecessor: u32,
        /// Target version.
        target: u32,
    },
    /// Target identity, preimage, or schema digest was all zero.
    InvalidTargetReceipt,
    /// Canonical migration identity publication refused.
    Identity(CanonicalError),
}

impl fmt::Display for CausalMigrationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ArtifactKindMismatch => {
                f.write_str("historical and native target artifact families differ")
            }
            Self::PredecessorNotOlder {
                predecessor,
                target,
            } => write!(
                f,
                "predecessor schema version {predecessor} must be older than target {target}"
            ),
            Self::InvalidTargetReceipt => f.write_str("target migration receipt is incomplete"),
            Self::Identity(error) => write!(f, "migration identity publication refused: {error}"),
        }
    }
}

impl std::error::Error for CausalMigrationError {}

impl CausalSchemaMigrationDraft {
    /// Bind a historical receipt to a native normalized-structure target.
    #[must_use]
    pub fn for_structure(
        predecessor: HistoricalCausalIdentityReceipt,
        target: &AdmittedCausalGraph,
        migration: CausalMigrationRef,
    ) -> Self {
        let receipt = target.structure_identity_receipt();
        Self {
            predecessor,
            target: NativeCausalTargetReceipt::Structure(receipt),
            migration,
        }
    }

    /// Bind a historical receipt to a native provenance-artifact target.
    #[must_use]
    pub fn for_graph_artifact(
        predecessor: HistoricalCausalIdentityReceipt,
        target: &AdmittedCausalGraph,
        migration: CausalMigrationRef,
    ) -> Self {
        let receipt = target.artifact_identity_receipt();
        Self {
            predecessor,
            target: NativeCausalTargetReceipt::GraphArtifact(receipt),
            migration,
        }
    }

    /// Bind a historical receipt to a native causalization-receipt target.
    #[must_use]
    pub fn for_causalization_receipt(
        predecessor: HistoricalCausalIdentityReceipt,
        target: &AdmittedCausalizationReceipt,
        migration: CausalMigrationRef,
    ) -> Self {
        let receipt = target.identity_receipt();
        Self {
            predecessor,
            target: NativeCausalTargetReceipt::Causalization(receipt),
            migration,
        }
    }

    /// Admit exact migration lineage without changing the target identity.
    ///
    /// # Errors
    /// Refuses a predecessor from another artifact family, a non-older
    /// predecessor, an incomplete target receipt, or bounded canonical
    /// publication error, including cancellation before publication.
    pub fn admit(self, cx: &Cx<'_>) -> Result<AdmittedCausalSchemaMigration, CausalMigrationError> {
        identity_materialization_checkpoint(cx, 0).map_err(CausalMigrationError::Identity)?;
        let kind = self.target.kind();
        let target_version = self.target.schema_version();
        if self.predecessor.artifact_kind != kind {
            return Err(CausalMigrationError::ArtifactKindMismatch);
        }
        if self.predecessor.schema_version >= target_version {
            return Err(CausalMigrationError::PredecessorNotOlder {
                predecessor: self.predecessor.schema_version,
                target: target_version,
            });
        }
        let target_identity = self.target.identity();
        let target_canonical_preimage = self.target.canonical_preimage();
        let target_schema_identity = self.target.schema_identity();
        if target_identity == [0; 32]
            || target_canonical_preimage == [0; 32]
            || target_schema_identity == [0; 32]
        {
            return Err(CausalMigrationError::InvalidTargetReceipt);
        }
        let receipt =
            causal_schema_migration_identity(&self, cx).map_err(CausalMigrationError::Identity)?;
        identity_materialization_checkpoint(cx, receipt.canonical_bytes() as usize)
            .map_err(CausalMigrationError::Identity)?;
        Ok(AdmittedCausalSchemaMigration {
            kind,
            predecessor: self.predecessor,
            target_identity,
            target_canonical_preimage,
            target_schema_identity,
            target_canonical_bytes: self.target.canonical_bytes(),
            target_field_count: self.target.field_count(),
            target_collection_items: self.target.collection_items(),
            migration: self.migration,
            receipt,
        })
    }
}

fn causal_schema_migration_identity(
    draft: &CausalSchemaMigrationDraft,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<CausalSchemaMigrationIdV1>, CanonicalError> {
    identity_materialization_checkpoint(cx, 0)?;
    let predecessor = historical_receipt_row(&draft.predecessor);
    let mut target = Vec::with_capacity(116);
    target.extend_from_slice(&draft.target.identity());
    target.extend_from_slice(&draft.target.canonical_preimage());
    target.extend_from_slice(&draft.target.schema_identity());
    target.extend_from_slice(&draft.target.canonical_bytes().to_le_bytes());
    target.extend_from_slice(&draft.target.field_count().to_le_bytes());
    target.extend_from_slice(&draft.target.collection_items().to_le_bytes());
    let mut migration = Vec::with_capacity(128);
    draft.migration.append_canonical(&mut migration);
    CanonicalEncoder::<CausalSchemaMigrationIdV1, _>::new(CAUSAL_RECEIPT_IDENTITY_LIMITS, || {
        cx.checkpoint().is_err()
    })?
    .u64(Field::new(0, "migration-receipt-schema-version"), 1)?
    .variant(
        Field::new(1, "artifact-kind"),
        u32::from(migration_artifact_kind_tag(draft.target.kind())),
        &[],
    )?
    .bytes(Field::new(2, "predecessor-receipt"), &predecessor)?
    .bytes(Field::new(3, "target-receipt"), &target)?
    .bytes(Field::new(4, "migration"), &migration)?
    .finish()
}

fn historical_receipt_row(receipt: &HistoricalCausalIdentityReceipt) -> Vec<u8> {
    let mut out = Vec::with_capacity(129);
    out.push(migration_artifact_kind_tag(receipt.artifact_kind));
    out.extend_from_slice(&receipt.schema_version.to_le_bytes());
    out.extend_from_slice(&receipt.semantic_identity);
    out.extend_from_slice(&receipt.canonical_preimage);
    out.extend_from_slice(&receipt.schema_identity);
    out.extend_from_slice(&receipt.canonical_bytes.to_le_bytes());
    out.extend_from_slice(&receipt.field_count.to_le_bytes());
    out.extend_from_slice(&receipt.collection_items.to_le_bytes());
    out
}

const fn migration_artifact_kind_tag(kind: CausalMigrationArtifactKind) -> u8 {
    match kind {
        CausalMigrationArtifactKind::Structure => 1,
        CausalMigrationArtifactKind::GraphArtifact => 2,
        CausalMigrationArtifactKind::CausalizationReceipt => 3,
    }
}

#[cfg(test)]
mod internal_tests {
    use super::*;
    use fs_blake3::identity::{OrderedBytesStreamDisposition, OrderedBytesStreamPhase};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    enum StreamParitySchemaV1 {}

    impl CanonicalSchema for StreamParitySchemaV1 {
        const DOMAIN: &'static str = "org.frankensim.fs-ir.test.causal-stream-parity.v1";
        const NAME: &'static str = "causal-stream-parity";
        const VERSION: u32 = 1;
        const CONTEXT: &'static str = "G3 eager-versus-streamed causal row identity parity";
        const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("rows", WireType::OrderedBytes)];
    }

    fn with_internal_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x57ea_0001,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn materialize_planned_test_row(
        planned_bytes: usize,
        write: impl FnOnce(&mut CausalIdentityRowWriter) -> Result<(), CanonicalError>,
    ) -> Vec<u8> {
        let plan = CausalIdentityRowPlan::new(0, planned_bytes).expect("test row plan");
        let mut row = Vec::new();
        reserve_identity_row(&mut row, plan).expect("test row reserve");
        let mut writer =
            CausalIdentityRowWriter::from_reserved(row, plan).expect("test checked row writer");
        let admitted_capacity = writer.capacity();
        write(&mut writer).expect("test row producer");
        assert_eq!(
            writer.capacity(),
            admitted_capacity,
            "test row producer must not grow past its admitted reservation"
        );
        let row = writer.into_bytes().expect("test row remains within plan");
        assert_eq!(row.len(), planned_bytes, "test row plan must be exact");
        row
    }

    fn injected_identity_reserve_failure(
        produced: &core::cell::Cell<usize>,
    ) -> CausalIdentityError {
        let rows = [b"planned".as_slice()];
        let limits = CanonicalLimits::new(4 * 1_024, 8, 1, 1, 16);
        with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows_with_reserve(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |row, _| Ok(row.len()),
                |row, _, out| {
                    produced.set(produced.get() + 1);
                    out.extend_from_slice(row);
                    Ok(())
                },
                |_, plan| {
                    Err(CausalIdentityAllocationRefusal {
                        requested_bytes: plan.canonical_bytes,
                    })
                },
            )
            .expect_err("injected reserve failure must consume without publication")
        })
    }

    #[test]
    fn g3_streamed_identity_rows_match_eager_canonical_bytes_and_roots() {
        let rows = vec![
            Vec::new(),
            b"one".to_vec(),
            vec![0, 1, 2, 0xff, 0x80],
            (0_u16..1_025)
                .flat_map(u16::to_le_bytes)
                .collect::<Vec<_>>(),
        ];
        let limits = CanonicalLimits::new(32 * 1_024, 8 * 1_024, 1, 4, 16);
        let eager =
            CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(limits, NeverCancel)
                .expect("valid eager encoder")
                .ordered_bytes(
                    Field::new(0, "rows"),
                    u64::try_from(rows.len()).expect("fixture cardinality fits u64"),
                    rows.iter().map(Vec::as_slice),
                )
                .expect("eager rows admit")
                .finish()
                .expect("eager identity publishes");

        let streamed = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid streamed encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |row, _| Ok(row.len()),
                |row, _, out| {
                    out.extend_from_slice(row);
                    Ok(())
                },
            )
            .expect("streamed rows admit")
            .finish()
            .expect("streamed identity publishes")
        });

        assert_eq!(streamed, eager);
        assert_eq!(streamed.id(), eager.id());
        assert_eq!(streamed.canonical_preimage(), eager.canonical_preimage());
        assert_eq!(streamed.canonical_bytes(), eager.canonical_bytes());
        assert_eq!(streamed.collection_items(), eager.collection_items());
    }

    #[test]
    fn g0_checked_row_terminal_encoding_matches_machine_wire_vocabulary() {
        use fs_qty::semantic::{
            AngleDomain, CompositionBasis, QuantityKind, SemanticType, StrainBasis,
            StrainComponent, ValueForm,
        };

        let quantities = [
            TerminalQuantitySpec::Dimensional(Dims([1, -2, 3, -4, 5, -6])),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Instantaneous,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Angle(AngleDomain::Mechanical),
                ValueForm::Peak,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Angle(AngleDomain::Electrical),
                ValueForm::Rms,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AngularVelocity(AngleDomain::Mechanical),
                ValueForm::Instantaneous,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AngularVelocity(AngleDomain::Electrical),
                ValueForm::Peak,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(QuantityKind::Torque, ValueForm::Rms)),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Energy,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Instantaneous,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Stress,
                ValueForm::Peak,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Tensor,
                    component: StrainComponent::Normal,
                },
                ValueForm::Rms,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Strain {
                    basis: StrainBasis::Engineering,
                    component: StrainComponent::Shear,
                },
                ValueForm::Instantaneous,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Composition(CompositionBasis::MassFraction),
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Composition(CompositionBasis::MoleFraction),
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Composition(CompositionBasis::VolumeFraction),
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Mass,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Amount,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::MolarMass,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::MassConcentration,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AmountConcentration,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::Entropy,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::HeatCapacity,
                ValueForm::Static,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AcousticPressure,
                ValueForm::Rms,
            )),
            TerminalQuantitySpec::Semantic(SemanticType::new(
                QuantityKind::AcousticPower,
                ValueForm::Static,
            )),
        ];
        let clock = ClockId::new("clock/i02-terminal-wire-parity").expect("valid parity clock");
        let frame = FrameBinding::new(
            "frame/i02-terminal-wire-parity",
            crate::machine::OrientationParity::Reversing,
        )
        .expect("valid parity frame");
        let append_signal_tail = |out: &mut Vec<u8>, signal: &SignalContract| {
            out.extend_from_slice(signal.clock.identity().as_bytes());
            let frame_bytes = u64::try_from(signal.frame.canonical_key().len())
                .expect("bounded frame key length fits u64");
            out.extend_from_slice(&frame_bytes.to_le_bytes());
            out.extend_from_slice(signal.frame.canonical_key().as_bytes());
            out.push(match signal.frame.orientation() {
                crate::machine::OrientationParity::Preserving => 1,
                crate::machine::OrientationParity::Reversing => 2,
            });
        };
        for quantity in quantities {
            let signal = SignalContract {
                quantity,
                shape: TerminalShape::Scalar,
                clock: clock.clone(),
                frame: frame.clone(),
            };
            let mut expected = Vec::new();
            super::super::push_terminal_quantity(&mut expected, quantity);
            assert_eq!(expected.len(), terminal_quantity_canonical_len(quantity));
            super::super::push_terminal_shape(&mut expected, signal.shape);
            append_signal_tail(&mut expected, &signal);
            let mut observed = Vec::new();
            push_signal(&mut observed, &signal).expect("checked signal row");
            assert_eq!(observed, expected, "quantity wire vocabulary drifted");
            assert_eq!(
                observed.len(),
                signal_canonical_len(&signal).expect("signal length")
            );
        }

        let shapes = [
            TerminalShape::Scalar,
            TerminalShape::Vector {
                components: NonZeroU64::new(3).expect("nonzero vector components"),
            },
            TerminalShape::Tensor {
                rows: NonZeroU64::new(2).expect("nonzero tensor rows"),
                columns: NonZeroU64::new(5).expect("nonzero tensor columns"),
            },
            TerminalShape::FieldTrace {
                components: NonZeroU64::new(7).expect("nonzero trace components"),
            },
        ];
        for shape in shapes {
            let signal = SignalContract {
                quantity: TerminalQuantitySpec::Dimensional(Dims([1, -2, 3, -4, 5, -6])),
                shape,
                clock: clock.clone(),
                frame: frame.clone(),
            };
            let mut expected = Vec::new();
            super::super::push_terminal_quantity(&mut expected, signal.quantity);
            let shape_start = expected.len();
            super::super::push_terminal_shape(&mut expected, shape);
            assert_eq!(
                expected.len() - shape_start,
                terminal_shape_canonical_len(shape)
            );
            append_signal_tail(&mut expected, &signal);
            let mut observed = Vec::new();
            push_signal(&mut observed, &signal).expect("checked signal row");
            assert_eq!(observed, expected, "terminal-shape wire vocabulary drifted");
            assert_eq!(
                observed.len(),
                signal_canonical_len(&signal).expect("signal length")
            );
        }
    }

    #[test]
    fn g3_actual_causal_row_length_plans_match_payload_families() {
        let version = NonZeroU64::new(1).expect("fixture version is nonzero");
        let owner = SubsystemId::new("subsystem/i02-stream-plan").expect("valid owner");
        let clock = ClockId::new("clock/i02-stream-plan").expect("valid clock");
        let condition_ref = ActivationConditionRef::new("i02/condition", version, [1; 32])
            .expect("valid condition");
        let branch =
            ActivationBranchRef::new("i02/branch/on", version, [2; 32]).expect("valid branch");
        let activation = ActivationDomain::Conditional {
            cubes: vec![ActivationCube {
                selections: vec![ConditionBranchSelection {
                    condition: condition_ref.clone(),
                    branch: branch.clone(),
                }],
            }],
        };
        let equation_lineage = NodeLineage::new(
            NodeOrigin::Machine(MachineNodeOrigin::Subsystem(owner.clone())),
            CausalOwner::Subsystem(owner.clone()),
            NormalizedNodeSemanticRef::new("i02/equation", version, [3; 32])
                .expect("valid equation meaning"),
        );
        let equation_id = EquationId::derive(&equation_lineage).expect("equation identity");
        let variable_lineage = NodeLineage::new(
            NodeOrigin::Machine(MachineNodeOrigin::Subsystem(owner.clone())),
            CausalOwner::Subsystem(owner.clone()),
            NormalizedNodeSemanticRef::new("i02/variable", version, [4; 32])
                .expect("valid variable meaning"),
        );
        let variable_id = VariableId::derive(&variable_lineage).expect("variable identity");
        let signal = SignalContract {
            quantity: TerminalQuantitySpec::Semantic(fs_qty::semantic::SemanticType::new(
                fs_qty::semantic::QuantityKind::Strain {
                    basis: fs_qty::semantic::StrainBasis::Tensor,
                    component: fs_qty::semantic::StrainComponent::Shear,
                },
                fs_qty::semantic::ValueForm::Static,
            )),
            shape: TerminalShape::Tensor {
                rows: NonZeroU64::new(2).expect("fixture row count is nonzero"),
                columns: NonZeroU64::new(3).expect("fixture column count is nonzero"),
            },
            clock,
            frame: FrameBinding::new(
                "frame/i02-stream-plan",
                crate::machine::OrientationParity::Preserving,
            )
            .expect("valid frame"),
        };
        let mut signal_row = Vec::new();
        push_signal(&mut signal_row, &signal).expect("signal row");
        assert_eq!(
            signal_canonical_len(&signal).expect("signal length"),
            signal_row.len()
        );
        let equation = EquationSpec {
            id: equation_id.clone(),
            diagnostic_label: "equation".into(),
            lineage: equation_lineage,
            owner: CausalOwner::Subsystem(owner.clone()),
            supports: vec![
                CausalSupport::Lumped,
                CausalSupport::External(
                    ExternalSupportRef::new("i02/support", version, [5; 32])
                        .expect("valid support"),
                ),
            ],
            residual: signal.clone(),
            role: EquationRole::Constraint,
            solve_participation: EquationParticipation::Matching,
            activation: activation.clone(),
        };
        let projection = PortSchemaCrosswalkRef::new("i02/projection", version, [6; 32])
            .expect("valid projection");
        let variable = VariableSpec {
            id: variable_id.clone(),
            diagnostic_label: "variable".into(),
            lineage: variable_lineage,
            owner: CausalOwner::Subsystem(owner),
            supports: vec![CausalSupport::Lumped],
            value: signal.clone(),
            role: VariableRole::Algebraic,
            solve_participation: SolveParticipation::Unknown,
            port_schema_crosswalk: Some(PortSchemaCrosswalkBinding {
                projection: projection.clone(),
                audit: PortSchemaCrosswalkAuditRef::new("i02/projection-audit", version, [7; 32])
                    .expect("valid projection audit"),
                audited_projection: projection,
            }),
            activation: activation.clone(),
        };
        let condition = ActivationConditionSpec {
            condition: condition_ref,
            source: ActivationConditionSource::GuardEquation {
                equation: equation_id.clone(),
                obligation: GuardSolveObligationRef::new("i02/guard-obligation", version, [8; 32])
                    .expect("valid guard obligation"),
            },
            branches: vec![branch],
            dependencies: vec![variable_id.clone()],
        };

        with_internal_cx(|cx| {
            let equation_len =
                equation_structure_row_len_cancellable(&equation, cx).expect("equation length");
            materialize_planned_test_row(equation_len, |row| {
                equation_structure_row(&equation, cx, row)
            });
            let variable_len =
                variable_structure_row_len_cancellable(&variable, cx).expect("variable length");
            materialize_planned_test_row(variable_len, |row| {
                variable_structure_row(&variable, cx, row)
            });
            let condition_len =
                condition_canonical_len_cancellable(&condition, cx).expect("condition length");
            materialize_planned_test_row(condition_len, |row| condition_row(&condition, cx, row));
            let condition_artifact_len =
                condition_artifact_row_len(&condition).expect("condition artifact length");
            materialize_planned_test_row(condition_artifact_len, |row| {
                condition_artifact_row(&condition, row);
                Ok(())
            });
            let equation_artifact_len =
                node_artifact_row_equation_len_cancellable(&equation.lineage, cx)
                    .expect("equation lineage length");
            materialize_planned_test_row(equation_artifact_len, |row| {
                node_artifact_row_equation(&equation.id, &equation.lineage, cx, row)
            });
            let variable_artifact_len = node_artifact_row_variable_len_cancellable(&variable, cx)
                .expect("variable lineage length");
            materialize_planned_test_row(variable_artifact_len, |row| {
                node_artifact_row_variable(&variable, cx, row)
            });

            let incidence = IncidenceSpec::new(
                equation_id.clone(),
                variable_id.clone(),
                0,
                SolveParticipation::Unknown,
                Dims::NONE,
                signal,
                None,
                IncidenceClockRelation::SameClock,
                activation,
                cx,
            )
            .expect("incidence identity");
            let incidence_artifact_len =
                incidence_artifact_row_len(&incidence).expect("incidence artifact length");
            materialize_planned_test_row(incidence_artifact_len, |row| {
                incidence_artifact_row(&incidence, row);
                Ok(())
            });
            let matching = CausalMatchingPair {
                incidence: incidence.id,
                equation: equation_id,
                variable: DerivativeVariableKey {
                    variable: variable_id,
                    derivative_order: 0,
                },
            };
            materialize_planned_test_row(3 * IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2, |row| {
                matching_row(&matching, row);
                Ok(())
            });
            materialize_planned_test_row(IDENTITY_RECEIPT_ADJUDICATION_BYTES + 2, |row| {
                derivative_variable_row(&matching.variable, row);
                Ok(())
            });
            let unknown = CausalUnknownAxisState {
                axis: CausalOutcomeAxis::Determination,
                reason: CausalUnknownReason::Cancelled,
                resume_checkpoint: Some(
                    CausalCheckpointRef::new("i02/checkpoint", version, [9; 32])
                        .expect("valid checkpoint"),
                ),
            };
            let unknown_len = unknown_axis_row_len(&unknown).expect("unknown-axis length");
            materialize_planned_test_row(unknown_len, |row| {
                unknown_axis_row(&unknown, row);
                Ok(())
            });
        });
    }

    #[test]
    fn g0_streamed_identity_rows_refuse_one_over_before_row_production() {
        let rows = [
            b"zero".to_vec(),
            b"one".to_vec(),
            b"two".to_vec(),
            b"three".to_vec(),
        ];
        let limits = CanonicalLimits::new(4 * 1_024, 1_024, 1, 3, 16);
        let planned = core::cell::Cell::new(0usize);
        let produced = core::cell::Cell::new(0usize);

        let refusal = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |row, _| {
                    planned.set(planned.get() + 1);
                    Ok(row.len())
                },
                |row, _, out| {
                    produced.set(produced.get() + 1);
                    out.extend_from_slice(row);
                    Ok(())
                },
            )
            .expect_err("four rows must refuse a three-row envelope")
        });

        assert_eq!(
            refusal.canonical_error(),
            Some(CanonicalError::LimitExceeded {
                kind: LimitKind::CollectionItems,
                requested: 4,
                limit: 3,
            })
        );
        let stream = refusal
            .stream_refusal()
            .expect("field refusal retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Canonical);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::FieldAdmission
        );
        assert_eq!(stream.diagnostic().row_index(), None);
        assert_eq!(stream.diagnostic().completed_rows(), 0);
        assert_eq!(planned.get(), 0, "field admission precedes length planning");
        assert_eq!(produced.get(), 0, "field admission precedes row allocation");
    }

    #[test]
    fn g0_oversized_stream_row_refuses_before_payload_allocation() {
        let rows = [()];
        let limits = CanonicalLimits::new(4 * 1_024, 8, 1, 1, 16);
        let produced = core::cell::Cell::new(0usize);

        let refusal = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |(), _| Ok(9),
                |(), _, out| {
                    produced.set(produced.get() + 1);
                    out.extend_from_slice(&[0; 9]);
                    Ok(())
                },
            )
            .expect_err("one-over row must refuse before its producer")
        });

        assert_eq!(
            refusal.canonical_error(),
            Some(CanonicalError::LimitExceeded {
                kind: LimitKind::FieldBytes,
                requested: 9,
                limit: 8,
            })
        );
        let stream = refusal
            .stream_refusal()
            .expect("row admission refusal retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Canonical);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::RowAdmission
        );
        assert_eq!(stream.diagnostic().row_index(), Some(0));
        assert_eq!(stream.diagnostic().declared_row_bytes(), Some(9));
        assert_eq!(stream.diagnostic().completed_rows(), 0);
        assert_eq!(
            produced.get(),
            0,
            "row admission precedes payload allocation"
        );
    }

    #[test]
    fn g0_streamed_identity_rows_admit_empty_and_exact_reduced_limit() {
        let rows = [Vec::new(), vec![0x5a; 8]];
        let limits = CanonicalLimits::new(4 * 1_024, 8, 1, 2, 16);

        let receipt = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid reduced-limit encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |row, _| Ok(row.len()),
                |row, _, out| {
                    out.extend_from_slice(row);
                    Ok(())
                },
            )
            .expect("empty and exact-cap rows admit")
            .finish()
            .expect("exact-cap identity publishes")
        });

        assert_eq!(receipt.collection_items(), 2);
    }

    #[test]
    fn g0_one_over_producer_refuses_before_mutation_growth_or_publication() {
        let rows = [()];
        let limits = CanonicalLimits::new(4 * 1_024, 8, 1, 1, 16);
        let capacity_before = core::cell::Cell::new(0usize);
        let capacity_after = core::cell::Cell::new(usize::MAX);
        let length_after = core::cell::Cell::new(usize::MAX);

        let refusal = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |(), _| Ok(4),
                |(), _, out| {
                    capacity_before.set(out.capacity());
                    out.extend_from_slice(&[1, 2, 3, 4, 5]);
                    capacity_after.set(out.capacity());
                    length_after.set(out.len());
                    Ok(())
                },
            )
            .expect_err("one-over producer must refuse without publishing")
        });

        assert_eq!(
            refusal.canonical_error(),
            Some(CanonicalError::DeclaredLengthMismatch {
                declared: 4,
                observed: 5,
            })
        );
        assert_eq!(refusal.allocation_refusal(), None);
        let stream = refusal
            .stream_refusal()
            .expect("one-over producer retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Producer);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::RowProducer
        );
        assert_eq!(stream.diagnostic().row_index(), Some(0));
        assert_eq!(stream.diagnostic().declared_row_bytes(), Some(4));
        assert_eq!(stream.diagnostic().written_row_bytes(), 0);
        assert_eq!(stream.diagnostic().completed_rows(), 0);
        assert_eq!(
            stream.diagnostic().disposition(),
            OrderedBytesStreamDisposition::EncoderConsumedNoPublication
        );
        assert_eq!(
            capacity_after.get(),
            capacity_before.get(),
            "refused one-over write cannot grow the admitted allocation"
        );
        assert_eq!(
            length_after.get(),
            0,
            "refused one-over write is checked before row mutation"
        );
    }

    #[test]
    fn g0_checked_row_single_byte_append_cannot_cross_admitted_cap() {
        let plan = CausalIdentityRowPlan::new(0, 4).expect("test row plan");
        let mut bytes = Vec::new();
        reserve_identity_row(&mut bytes, plan).expect("test row reserve");
        let mut writer =
            CausalIdentityRowWriter::from_reserved(bytes, plan).expect("checked row writer");
        for byte in [1, 2, 3, 4] {
            writer.push(byte);
        }
        let admitted_capacity = writer.capacity();
        writer.push(5);

        assert_eq!(writer.len(), 4, "refused byte is not appended");
        assert_eq!(
            writer.capacity(),
            admitted_capacity,
            "refused byte cannot grow storage"
        );
        assert_eq!(
            writer.refusal(),
            Some(CanonicalError::DeclaredLengthMismatch {
                declared: 4,
                observed: 5,
            })
        );
        assert_eq!(
            writer.into_bytes(),
            Err(CanonicalError::DeclaredLengthMismatch {
                declared: 4,
                observed: 5,
            }),
            "sticky over-plan refusal prevents the row from reaching a sink"
        );
    }

    #[test]
    fn g0_streamed_identity_rows_distinguish_injected_reserve_failure() {
        let produced = core::cell::Cell::new(0usize);
        let refusal = injected_identity_reserve_failure(&produced);

        assert_eq!(refusal.canonical_error(), None);
        assert_eq!(
            refusal.allocation_refusal(),
            Some(CausalIdentityAllocationRefusal { requested_bytes: 7 })
        );
        let stream = refusal
            .stream_refusal()
            .expect("reserve refusal retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Allocation);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::RowProducer
        );
        assert_eq!(stream.diagnostic().row_index(), Some(0));
        assert_eq!(stream.diagnostic().declared_row_bytes(), Some(7));
        assert_eq!(stream.diagnostic().written_row_bytes(), 0);
        assert_eq!(produced.get(), 0, "reserve refusal precedes serialization");
    }

    #[test]
    fn g0_graph_and_receipt_refusals_preserve_typed_identity_allocation_data() {
        let produced = core::cell::Cell::new(0usize);
        let identity_refusal = injected_identity_reserve_failure(&produced);
        let expected_allocation = CausalIdentityAllocationRefusal { requested_bytes: 7 };
        let expected_stream = identity_refusal
            .stream_refusal()
            .expect("injected allocation refusal carries stream context");

        let graph_refusal = identity_graph_refusal(identity_refusal);
        assert_eq!(graph_refusal.identity_error(), None);
        assert_eq!(
            graph_refusal.identity_allocation_refusal(),
            Some(expected_allocation)
        );
        assert_eq!(
            graph_refusal.identity_stream_refusal(),
            Some(expected_stream)
        );

        let receipt_refusal = identity_receipt_refusal(identity_refusal);
        assert_eq!(receipt_refusal.identity_error(), None);
        assert_eq!(
            receipt_refusal.identity_allocation_refusal(),
            Some(expected_allocation)
        );
        assert_eq!(
            receipt_refusal.identity_stream_refusal(),
            Some(expected_stream)
        );
        assert_eq!(
            produced.get(),
            0,
            "allocation refusal never invokes a producer"
        );
    }

    #[test]
    fn g0_streamed_identity_rows_refuse_length_payload_disagreement() {
        let rows = [()];
        let limits = CanonicalLimits::new(4 * 1_024, 8, 1, 1, 16);

        let refusal = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |(), _| Ok(4),
                |(), _, out| {
                    out.extend_from_slice(&[1, 2, 3]);
                    Ok(())
                },
            )
            .expect_err("short payload must refuse exact row completion")
        });

        assert_eq!(
            refusal.canonical_error(),
            Some(CanonicalError::DeclaredLengthMismatch {
                declared: 4,
                observed: 3,
            })
        );
        let stream = refusal
            .stream_refusal()
            .expect("length mismatch retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Canonical);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::RowCompletion
        );
        assert_eq!(stream.diagnostic().written_row_bytes(), 3);
    }

    #[test]
    fn g4_stream_row_producer_cancellation_retains_completed_prefix() {
        let rows = [b"first".as_slice(), b"second".as_slice()];
        let limits = CanonicalLimits::new(4 * 1_024, 1_024, 1, 2, 16);
        let produced = core::cell::Cell::new(0usize);

        let refusal = with_internal_cx(|cx| {
            let encoder = CanonicalEncoder::<EvidenceNodeId<StreamParitySchemaV1>, _>::new(
                limits,
                NeverCancel,
            )
            .expect("valid bounded encoder");
            stream_identity_rows(
                encoder,
                Field::new(0, "rows"),
                &rows,
                cx,
                |row, _| Ok(row.len()),
                |row, _, out| {
                    let index = produced.get();
                    produced.set(index + 1);
                    if index == 1 {
                        Err(CanonicalError::Cancelled { absorbed_bytes: 0 })
                    } else {
                        out.extend_from_slice(row);
                        Ok(())
                    }
                },
            )
            .expect_err("second producer cancellation must consume without publication")
        });

        assert_eq!(
            refusal.canonical_error(),
            Some(CanonicalError::Cancelled { absorbed_bytes: 0 })
        );
        let stream = refusal
            .stream_refusal()
            .expect("producer cancellation retains stream context");
        assert_eq!(stream.origin(), CausalIdentityStreamOrigin::Producer);
        assert_eq!(
            stream.diagnostic().phase(),
            OrderedBytesStreamPhase::RowProducer
        );
        assert_eq!(stream.diagnostic().row_index(), Some(1));
        assert_eq!(stream.diagnostic().completed_rows(), 1);
        assert!(stream.diagnostic().canonical_bytes() > 0);
        assert_eq!(produced.get(), 2);
    }

    #[test]
    fn g0_determination_order_tracks_preserved_wire_tags() {
        let states = [
            DeterminationClass::WellDetermined,
            DeterminationClass::UnderDetermined,
            DeterminationClass::OverDetermined,
            DeterminationClass::Mixed,
            DeterminationClass::Unknown,
            DeterminationClass::EmptyProjection,
        ];
        assert_eq!(states.map(determination_tag), [1, 2, 3, 4, 5, 6]);
        assert!(
            states.windows(2).all(|pair| pair[0] < pair[1]),
            "public Ord must remain aligned with canonical wire-tag order"
        );
    }

    /// Test-only constructor capability for the nominal identity roles used by
    /// these law fixtures. `StrongIdentity` intentionally has no generic
    /// constructor, so a `StrongIdentity` bound alone must not pretend that
    /// `CanonicalEncoder::<I, _>::new` exists for every sealed role.
    trait TestIdentityEncoder: StrongIdentity {
        fn test_encoder(
            limits: CanonicalLimits,
        ) -> Result<CanonicalEncoder<Self, NeverCancel>, CanonicalError>;
    }

    impl<D: CanonicalSchema> TestIdentityEncoder for EntityId<D> {
        fn test_encoder(
            limits: CanonicalLimits,
        ) -> Result<CanonicalEncoder<Self, NeverCancel>, CanonicalError> {
            CanonicalEncoder::<EntityId<D>, _>::new(limits, NeverCancel)
        }
    }

    impl<D: CanonicalSchema> TestIdentityEncoder for EvidenceNodeId<D> {
        fn test_encoder(
            limits: CanonicalLimits,
        ) -> Result<CanonicalEncoder<Self, NeverCancel>, CanonicalError> {
            CanonicalEncoder::<EvidenceNodeId<D>, _>::new(limits, NeverCancel)
        }
    }

    fn single_bytes_receipt<I: TestIdentityEncoder>(
        field_name: &'static str,
        value: &[u8],
        limits: CanonicalLimits,
    ) -> IdentityReceipt<I> {
        I::test_encoder(limits)
            .expect("test encoder")
            .bytes(Field::new(0, field_name), value)
            .expect("test identity field")
            .finish()
            .expect("test identity receipt")
    }

    fn single_ordered_bytes_receipt<I: TestIdentityEncoder>(
        field_name: &'static str,
        value: &[u8],
        limits: CanonicalLimits,
    ) -> IdentityReceipt<I> {
        I::test_encoder(limits)
            .expect("test encoder")
            .ordered_bytes(Field::new(0, field_name), 1, core::iter::once(value))
            .expect("test ordered identity field")
            .finish()
            .expect("test identity receipt")
    }

    fn hash_value<T: Hash>(value: &T) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        value.hash(&mut hasher);
        hasher.finish()
    }

    fn hash_adjudicated_receipt<I: StrongIdentity>(receipt: IdentityReceipt<I>) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        hash_identity_receipt_adjudication(receipt, &mut hasher);
        hasher.finish()
    }

    fn assert_adjudicated_id_laws<T>(left: &T, same: &T, different: &T)
    where
        T: Eq + Ord + Hash + fmt::Debug,
    {
        assert_eq!(left, same, "evidence-only encoder limits changed equality");
        assert_eq!(
            left.cmp(same),
            core::cmp::Ordering::Equal,
            "equality and ordering disagree"
        );
        assert_eq!(
            hash_value(left),
            hash_value(same),
            "equal adjudicated IDs must hash identically"
        );
        assert_ne!(left, different, "different canonical meaning collapsed");
        assert_ne!(
            left.cmp(different),
            core::cmp::Ordering::Equal,
            "different canonical meaning compared equal"
        );
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    struct StableItem {
        key: usize,
        insertion_order: usize,
    }

    #[test]
    fn g0_adjudicated_id_eq_ord_hash_ignore_only_evidence_limits() {
        let limits_a = CanonicalLimits::new(16_384, 8_192, 1, 1, 64);
        let limits_b = CanonicalLimits::new(32_768, 16_384, 1, 1, 512);

        let owner = SubsystemId::new("test/adjudication-owner").expect("test owner");
        let lineage = NodeLineage::new(
            NodeOrigin::Machine(MachineNodeOrigin::Subsystem(owner.clone())),
            CausalOwner::Subsystem(owner),
            NormalizedNodeSemanticRef::new(
                "test/adjudication-node",
                NonZeroU64::new(1).expect("nonzero schema version"),
                [41; 32],
            )
            .expect("test node semantics"),
        );
        let normalized = lineage.normalized_row();
        let equation_reencoded = EquationId {
            receipt: Arc::new(single_bytes_receipt::<EquationEntityIdV1>(
                "instance-qualified-normalized-meaning",
                &normalized,
                limits_b,
            )),
        };
        let variable_reencoded = VariableId {
            receipt: Arc::new(single_bytes_receipt::<VariableEntityIdV1>(
                "instance-qualified-normalized-meaning",
                &normalized,
                limits_b,
            )),
        };
        assert!(
            equation_identity_matches_lineage(&equation_reencoded, &lineage),
            "graph admission must accept an adjudication-equivalent equation receipt"
        );
        assert!(
            variable_identity_matches_lineage(&variable_reencoded, &lineage),
            "graph admission must accept an adjudication-equivalent variable receipt"
        );

        let equation_a = EquationId {
            receipt: Arc::new(single_bytes_receipt::<EquationEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"equation-a",
                limits_a,
            )),
        };
        let equation_same = EquationId {
            receipt: Arc::new(single_bytes_receipt::<EquationEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"equation-a",
                limits_b,
            )),
        };
        let equation_different = EquationId {
            receipt: Arc::new(single_bytes_receipt::<EquationEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"equation-b",
                limits_a,
            )),
        };
        assert_ne!(
            equation_a.identity_receipt(),
            equation_same.identity_receipt(),
            "fixture must differ in evidence-only limits"
        );
        assert!(identity_receipt_adjudication_eq(
            equation_a.identity_receipt(),
            equation_same.identity_receipt()
        ));
        assert_eq!(
            identity_receipt_adjudication_cmp(
                equation_a.identity_receipt(),
                equation_same.identity_receipt()
            ),
            core::cmp::Ordering::Equal
        );
        assert_eq!(
            hash_adjudicated_receipt(equation_a.identity_receipt()),
            hash_adjudicated_receipt(equation_same.identity_receipt())
        );
        assert_adjudicated_id_laws(&equation_a, &equation_same, &equation_different);

        let variable_a = VariableId {
            receipt: Arc::new(single_bytes_receipt::<VariableEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"variable-a",
                limits_a,
            )),
        };
        let variable_same = VariableId {
            receipt: Arc::new(single_bytes_receipt::<VariableEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"variable-a",
                limits_b,
            )),
        };
        let variable_different = VariableId {
            receipt: Arc::new(single_bytes_receipt::<VariableEntityIdV1>(
                "instance-qualified-normalized-meaning",
                b"variable-b",
                limits_a,
            )),
        };
        assert_adjudicated_id_laws(&variable_a, &variable_same, &variable_different);

        let incidence_a = IncidenceId {
            receipt: Arc::new(single_bytes_receipt::<IncidenceEntityIdV1>(
                "normalized-incidence-meaning",
                b"incidence-a",
                limits_a,
            )),
        };
        let incidence_same = IncidenceId {
            receipt: Arc::new(single_bytes_receipt::<IncidenceEntityIdV1>(
                "normalized-incidence-meaning",
                b"incidence-a",
                limits_b,
            )),
        };
        let incidence_different = IncidenceId {
            receipt: Arc::new(single_bytes_receipt::<IncidenceEntityIdV1>(
                "normalized-incidence-meaning",
                b"incidence-b",
                limits_a,
            )),
        };
        assert_adjudicated_id_laws(&incidence_a, &incidence_same, &incidence_different);

        let outcome_set_a = single_ordered_bytes_receipt::<ConditionalOutcomeSetIdV1>(
            "mode-cell-outcomes",
            b"outcome-a",
            limits_a,
        );
        let outcome_set_same = single_ordered_bytes_receipt::<ConditionalOutcomeSetIdV1>(
            "mode-cell-outcomes",
            b"outcome-a",
            limits_b,
        );
        assert_eq!(
            ConditionalCoverageClaim::ModeCells(outcome_set_a),
            ConditionalCoverageClaim::ModeCells(outcome_set_same),
            "private theorem claims must not reintroduce evidence-limit-sensitive equality"
        );
    }

    #[test]
    fn g0_receipt_finding_budget_is_inclusive_and_fails_closed_afterward() {
        let finding = CausalReceiptFinding::new(
            CausalReceiptRule::OutcomeAxisMismatch,
            CausalReceiptSubject::Receipt,
        );
        let mut findings = vec![finding.clone(); MAX_CAUSAL_RECEIPT_FINDINGS];
        enforce_receipt_finding_budget(&findings).expect("inclusive finding cap");
        findings.push(finding);
        let refusal = enforce_receipt_finding_budget(&findings)
            .expect_err("one finding above the cap must fail closed");
        assert_eq!(refusal.findings().len(), 1);
        assert_eq!(
            refusal.findings()[0].rule(),
            CausalReceiptRule::ResourceLimit
        );
    }

    #[test]
    fn g0_graph_finding_budget_is_inclusive_and_fails_closed_afterward() {
        let finding = CausalGraphFinding::new(
            CausalGraphRule::InvalidDiagnosticLabel,
            CausalGraphSubject::Graph,
        );
        let mut findings = vec![finding.clone(); MAX_CAUSAL_GRAPH_FINDINGS];
        enforce_graph_finding_budget(&findings).expect("inclusive finding cap");
        findings.push(finding);
        let refusal = enforce_graph_finding_budget(&findings)
            .expect_err("one finding above the cap must fail closed");
        assert_eq!(refusal.findings().len(), 1);
        assert_eq!(refusal.findings()[0].rule(), CausalGraphRule::ResourceLimit);
    }

    #[test]
    fn g0_cancellable_sort_matches_stable_standard_sort_across_boundaries() {
        let mut cases: Vec<Vec<StableItem>> = vec![
            Vec::new(),
            vec![StableItem {
                key: 7,
                insertion_order: 0,
            }],
            (0..257)
                .rev()
                .map(|value| StableItem {
                    key: value,
                    insertion_order: 256 - value,
                })
                .collect(),
            (0..1_025)
                .map(|insertion_order| StableItem {
                    key: (insertion_order * 37 + 11) % 23,
                    insertion_order,
                })
                .collect(),
        ];

        for (case_index, values) in cases.iter_mut().enumerate() {
            let mut expected = values.clone();
            expected.sort_by_key(|item| item.key);
            let mut checkpoints = 0usize;
            cancellable_sort_by(
                values,
                |left, right| left.key.cmp(&right.key),
                || {
                    checkpoints = checkpoints.saturating_add(1);
                    Ok::<(), &'static str>(())
                },
            )
            .unwrap_or_else(|error| {
                panic!(
                    "case={case_index} unexpectedly interrupted; checkpoints={checkpoints}; error={error}"
                )
            });
            let actual_order: Vec<_> = values
                .iter()
                .map(|item| (item.key, item.insertion_order))
                .collect();
            let expected_order: Vec<_> = expected
                .iter()
                .map(|item| (item.key, item.insertion_order))
                .collect();
            assert_eq!(
                actual_order, expected_order,
                "case={case_index}; checkpoints={checkpoints}; stable order diverged"
            );
            assert!(
                checkpoints > 0,
                "case={case_index} executed without a cancellation checkpoint"
            );
        }

        let mut natural = vec![5_i32, 3, 3, -1, 8, 0, 5];
        cancellable_sort(&mut natural, || Ok::<(), &'static str>(()))
            .expect("natural-order wrapper completes");
        assert_eq!(natural, vec![-1, 0, 3, 3, 5, 5, 8]);
    }

    #[test]
    fn g4_cancellable_sort_interrupts_each_observed_phase() {
        let phase_failures = [
            (CausalSortPhase::Entry, 1usize),
            (CausalSortPhase::IndexInitialization, 2usize),
            (CausalSortPhase::Merge, 2usize),
            (CausalSortPhase::InverseInitialization, 2usize),
            (CausalSortPhase::InverseMap, 2usize),
            (CausalSortPhase::PayloadPosition, 2usize),
            (CausalSortPhase::PayloadSwap, 2usize),
            (CausalSortPhase::Complete, 1usize),
        ];
        for (phase, fail_at) in phase_failures {
            let original: Vec<usize> = (0..4_096).rev().collect();
            let mut values = original.clone();
            let mut phase_checkpoints = 0usize;
            let result = cancellable_sort_by_fallible_observed(
                &mut values,
                |left, right| Ok::<_, (CausalSortPhase, usize)>(left.cmp(right)),
                |observed_phase| {
                    if observed_phase == phase {
                        phase_checkpoints = phase_checkpoints.saturating_add(1);
                    }
                    if observed_phase == phase && phase_checkpoints == fail_at {
                        Err((phase, phase_checkpoints))
                    } else {
                        Ok(())
                    }
                },
            );
            assert_eq!(
                result,
                Err((phase, fail_at)),
                "phase={phase:?}; expected injected phase checkpoint={fail_at}; observed={phase_checkpoints}"
            );

            let mut observed_multiset = values.clone();
            observed_multiset.sort_unstable();
            let mut expected_multiset = original.clone();
            expected_multiset.sort_unstable();
            assert_eq!(
                observed_multiset, expected_multiset,
                "phase={phase:?}; interruption changed the payload multiset"
            );
            if phase == CausalSortPhase::Complete {
                assert!(
                    values.windows(2).all(|pair| pair[0] <= pair[1]),
                    "complete-phase interruption must retain the fully sorted payload"
                );
            } else if matches!(
                phase,
                CausalSortPhase::PayloadPosition | CausalSortPhase::PayloadSwap
            ) {
                assert_ne!(
                    values, original,
                    "phase={phase:?}; interruption did not exercise partial payload movement"
                );
            } else {
                assert_eq!(
                    values, original,
                    "phase={phase:?}; index-only interruption moved the payload"
                );
            }
        }

        let mut identity_permutation: Vec<usize> = (0..4_096).collect();
        let identity_before = identity_permutation.clone();
        let mut position_checkpoints = 0usize;
        let result = cancellable_sort_by_fallible_observed(
            &mut identity_permutation,
            |left, right| Ok::<_, (CausalSortPhase, usize)>(left.cmp(right)),
            |phase| {
                if phase == CausalSortPhase::PayloadPosition {
                    position_checkpoints = position_checkpoints.saturating_add(1);
                }
                if phase == CausalSortPhase::PayloadPosition && position_checkpoints == 3 {
                    Err((phase, position_checkpoints))
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(
            result,
            Err((CausalSortPhase::PayloadPosition, 3)),
            "identity permutation must remain cancellable across long fixed-point runs"
        );
        assert_eq!(
            identity_permutation, identity_before,
            "fixed-point interruption unexpectedly moved an already canonical payload"
        );
    }

    #[test]
    fn g4_fallible_long_row_comparator_interrupts_before_payload_publication() {
        let condition = ActivationConditionRef::new(
            "test/sort-comparator-condition",
            NonZeroU64::new(1).expect("nonzero test version"),
            [1; 32],
        )
        .expect("condition reference");
        let common_branch = ActivationBranchRef::new(
            "test/sort-comparator-common-branch",
            NonZeroU64::new(1).expect("nonzero test version"),
            [2; 32],
        )
        .expect("branch reference");
        let differing_branch = ActivationBranchRef::new(
            "test/sort-comparator-differing-branch",
            NonZeroU64::new(1).expect("nonzero test version"),
            [3; 32],
        )
        .expect("branch reference");
        let common = ConditionBranchSelection {
            condition: condition.clone(),
            branch: common_branch,
        };
        let left = ActivationCube {
            selections: vec![common.clone(); 513],
        };
        let mut right = left.clone();
        right.selections[512] = ConditionBranchSelection {
            condition,
            branch: differing_branch,
        };
        let mut rows = vec![right, left];
        let original = rows.clone();
        let mut comparator_checkpoints = 0usize;
        let result = cancellable_sort_by_fallible(
            &mut rows,
            |left, right| {
                compare_activation_cubes_cancellable(left, right, || {
                    comparator_checkpoints = comparator_checkpoints.saturating_add(1);
                    if comparator_checkpoints == 3 {
                        Err("comparator-cancelled")
                    } else {
                        Ok(())
                    }
                })
            },
            || Ok::<(), &'static str>(()),
        );
        assert_eq!(result, Err("comparator-cancelled"));
        assert_eq!(
            rows, original,
            "fallible index comparison must not publish a partial payload order"
        );
    }
}
