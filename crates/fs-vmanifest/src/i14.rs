//! The I14 (multirung EMC and harness) VerificationManifest draft
//! (bead frankensim-leapfrog-2026-program-i94v.2.4.9.1).
//!
//! The independently promotable core covers native harness identity, typed
//! RLGC/MTL and PEEC models, grounding/bonding/shielding and bearing-current
//! paths, separately owned source/probe and victim semantics, fixed-rung
//! routing, adjoints, UQ inference mechanics, and independently promotable
//! full-wave schema admission. The maximal lattice separately owns full-wave
//! FEEC, BEM formulation correctness, FMM acceleration envelopes, certified
//! cross-rung error descent, blind robust mitigation, synthetic safety-case
//! traceability, governed EMC reliability, passive/causal sheaf-cosheaf
//! composition, hypercohomology-obstruction localization, cover-refinement
//! naturality, a KYP/sheaf passivity bridge, maximal falsification, governed
//! standards crosswalks, calibrated laboratory validation, and production
//! bearing-population reliability. A weaker receipt closes only its own lattice
//! element and is never relabeled as a stronger physical, theorem, standards,
//! population, safety, or regulatory claim.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};
use fs_blake3::{ContentHash, hash_domain};
use std::collections::BTreeSet;

const CAMPAIGN_POLICY_FIXTURE: &str = "i14-campaign-policy-v1";
const ACCEPTANCE_POLICY_FIXTURE: &str = "i14-acceptance-arithmetic-policy-v1";
const EM_CONVENTION_FIXTURE: &str = "i14-em-convention-card-v1";
const THEOREM_POLICY_FIXTURE: &str = "i14-theorem-formalization-policy-v1";

/// Domain separator for the exhaustive I14 terminal-status truth table.
const TERMINAL_STATUS_TABLE_DOMAIN_V1: &str = "org.frankensim.i14.terminal-status-truth-table.v1";

/// Domain separator for canonical terminal-result identity.
const TERMINAL_RESULT_DOMAIN_V1: &str = "org.frankensim.i14.terminal-result.v1";

/// Domain separator for receipt-bound, noncanonical timing telemetry.
const TELEMETRY_ENVELOPE_DOMAIN_V1: &str = "org.frankensim.i14.telemetry-envelope.v1";

/// Domain separator for validated I14 cancellation-card identity.
const CANCELLATION_CARD_DOMAIN_V2: &str = "org.frankensim.i14.cancellation-card.v2";

/// Domain separator for the authoritative first-terminal boundary prefix.
const TERMINAL_TRACE_DOMAIN_V2: &str = "org.frankensim.i14.terminal-trace.v2";

/// Domain separator for the recomputed clock-free logical execution stream.
const LOGICAL_EXECUTION_TRACE_DOMAIN_V2: &str = "org.frankensim.i14.logical-execution-trace.v2";

/// Domain separator for authoritative canonical terminal-result identity.
const TERMINAL_RESULT_DOMAIN_V2: &str = "org.frankensim.i14.terminal-result.v2";

/// Domain separator for receipt-bound, explicitly noncanonical timing telemetry.
const TELEMETRY_ENVELOPE_DOMAIN_V2: &str = "org.frankensim.i14.telemetry-envelope.v2";

/// Domain separator for a complete raw watchdog observation trace.
const WATCHDOG_RAW_TRACE_DOMAIN_V2: &str = "org.frankensim.i14.watchdog-raw-trace.v2";

/// Number of raw tuples in the exhaustive I14 terminal-status truth table.
pub const I14_TERMINAL_STATUS_TABLE_V1_TUPLES: usize = 3_600;

/// Maximum cancellation requests accepted by one cause-selector invocation.
pub const I14_MAX_CANCELLATION_REQUESTS_V1: usize = 16_384;

/// Maximum root-to-leaf scope ancestry accepted by the cause selector.
pub const I14_MAX_SCOPE_ANCESTRY_V1: usize = 256;

/// Maximum observer tiles admitted at one terminal boundary.
pub const I14_MAX_OBSERVER_TILES_V1: usize = 128;

/// Maximum watchdog observations admitted to one telemetry envelope.
pub const I14_MAX_WATCHDOG_OBSERVATIONS_V1: usize = 4_096;

/// Maximum genesis-to-terminal boundary records in one authoritative trace.
pub const I14_MAX_TERMINAL_BOUNDARIES_V2: usize = 4_096;

/// Maximum boundary/request arbitration pairs admitted to one V2 trace.
///
/// Both independent collection caps remain available, but their Cartesian
/// product is explicitly budgeted so the legacy per-boundary reference
/// selector cannot become an allocation-and-sort denial of service. A future
/// incremental selector may raise this versioned budget after equivalence and
/// performance evidence are frozen.
pub const I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2: usize = 1_048_576;

/// Genesis ordinal required by an authoritative I14 terminal trace.
pub const I14_TERMINAL_BOUNDARY_GENESIS_ORDINAL_V2: u64 = 0;

/// Known-answer digest of the canonical I14 terminal-status table.
///
/// The canonical payload starts with `I14_TERMINAL_STATUS_TRUTH_TABLE_V1\0`,
/// followed by the tuple count as little-endian `u32`. Every tuple follows in
/// lexicographic axis-tag order and contributes its eight raw tags, eight
/// normalized tags, two explicit normalization-action bits, and primary exit
/// code. [`i14_terminal_status_table_digest_v1`] hashes that payload under
/// `TERMINAL_STATUS_TABLE_DOMAIN_V1`.
pub const I14_TERMINAL_STATUS_TABLE_V1_DIGEST_HEX: &str =
    "ea696d41442ee10ebf44da163a90f3b6f55269721c3fde9ec106e4cf1e5fe6b0";

/// Known-answer digest for the canonical terminal-result fixture in the I14
/// conformance suite. This pins the version-1 canonical byte layout in
/// addition to its relational invariants.
pub const I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX: &str =
    "8bdedaea44b6d2b746d3ffe280ec590bf092fc32a7810a183a70b446cadc764f";

/// Known-answer digest for the receipt-bound telemetry fixture paired with
/// [`I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX`].
pub const I14_TELEMETRY_ENVELOPE_V1_KAT_HEX: &str =
    "da64b4c869861e616933c1f84b9595d9eda0e0e0c865037914e0965ac440f735";

/// Known-answer digest for the Core card used by the V2 conformance fixture.
/// Independent encoder reproduction remains a separate promotion gate.
pub const I14_CANCELLATION_CARD_V2_KAT_HEX: &str =
    "ff109fdef04188dcf47edf7569fe5c363b113d1a113181c0650a5803bc785f2d";

/// Known-answer digest for the selected request-inclusive terminal prefix
/// paired with the V2 conformance fixture.
pub const I14_TERMINAL_PREFIX_V2_KAT_HEX: &str =
    "797687334c795ad8b10497040a404ed5bff9fed260189c2cccd160427b5a0b48";

/// Exhaustive exact-byte KAT for all four I14DrainTriggerV2 wire variants
/// concatenated in normative tag order. The union uses one tag byte followed,
/// for tags 1 through 3, by one unframed little-endian u64.
pub const I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX: &str =
    "00010102030405060708021112131415161718032122232425262728";

/// Exact-byte KAT for a present canonical infrastructure-failure onset witness:
/// presence byte 1, sequence 0x0807060504030201, source tag Supervisor=5,
/// then 32 receipt bytes equal to 0xa7. Calibrated onset time is telemetry.
pub const I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX: &str =
    "01010203040506070805a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7a7";

/// Known-answer digest for the complete, multi-kind raw-watchdog trace paired
/// with the V2 conformance fixture. The fixture contains 33 poll observations
/// plus one external-heartbeat and one deadline observation; its exact encoded
/// payload is 629 bytes.
pub const I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX: &str =
    "a5035cc151cc1ad0952e7939a4b380b0514b1cff9def0507f9bbd0756de8b200";

/// Known-answer digest for the authoritative V2 canonical fixture.
pub const I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX: &str =
    "4108e92b5dbee524898a5c457dd08e90fd774bf75ea810710622a2c5917de4cb";

/// Known-answer digest for the paired noncanonical V2 telemetry.
pub const I14_TELEMETRY_ENVELOPE_V2_KAT_HEX: &str =
    "8b8d52877994335ae38314567343519243b228ff2ba88b65a66fa66ac14047af";

/// Whether the campaign compute reached a normal terminal boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14ExecutionDisposition {
    /// All scheduled work completed before any competing terminal cause.
    Completed = 0,
    /// A recorded cancellation request won terminal-cause adjudication.
    Cancelled = 1,
    /// A wall budget or cancellation/drain/finalization deadline expired.
    TimedOut = 2,
    /// A declared total resource ceiling was reached before completion.
    BudgetExhausted = 3,
    /// The supervisor, authentication, drain, or publication protocol failed.
    InfrastructureFailed = 4,
}

impl I14ExecutionDisposition {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 5] = [
        Self::Completed,
        Self::Cancelled,
        Self::TimedOut,
        Self::BudgetExhausted,
        Self::InfrastructureFailed,
    ];
}

/// Scientific adjudication of the requested claim, independent of execution.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14ClaimAdjudication {
    /// Evidence collection or adjudication has not reached a scientific verdict.
    Pending = 0,
    /// The admitted evidence supports the frozen claim.
    Supported = 1,
    /// The admitted evidence fails the frozen acceptance predicate.
    Failed = 2,
    /// A strength-matched admitted counterexample refutes the exact claim revision.
    Refuted = 3,
    /// The available admitted information does not determine a verdict.
    Unknown = 4,
}

impl I14ClaimAdjudication {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 5] = [
        Self::Pending,
        Self::Supported,
        Self::Failed,
        Self::Refuted,
        Self::Unknown,
    ];
}

/// Completeness of evidence required by the frozen acceptance card.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14EvidenceCompleteness {
    /// Every required component and supporting artifact is present.
    CompleteEvidence = 0,
    /// At least one required component is missing but some evidence is present.
    PartialEvidence = 1,
    /// No admissible evidence for the requested claim is present.
    NoEvidence = 2,
}

impl I14EvidenceCompleteness {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 3] = [
        Self::CompleteEvidence,
        Self::PartialEvidence,
        Self::NoEvidence,
    ];
}

/// Integrity of evidence bytes, custody, lineage, and checker execution.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14EvidenceIntegrity {
    /// All applicable integrity predicates were verified.
    IntegrityVerified = 0,
    /// At least one integrity predicate failed or could not be authenticated.
    IntegrityFailed = 1,
}

impl I14EvidenceIntegrity {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 2] = [Self::IntegrityVerified, Self::IntegrityFailed];
}

/// Structural and semantic validity of the submitted input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14InputValidity {
    /// The input satisfies its frozen structural and semantic schema.
    WellFormedInput = 0,
    /// The input is structurally or semantically malformed.
    MalformedInput = 1,
}

impl I14InputValidity {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 2] = [Self::WellFormedInput, Self::MalformedInput];
}

/// Applicability of the frozen claim domain to a well-formed input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14DomainApplicability {
    /// The input is inside the frozen claim domain.
    Admitted = 0,
    /// The input is well formed but outside the frozen claim domain.
    OutOfDomain = 1,
    /// Applicability cannot be determined, including for malformed input.
    Indeterminate = 2,
}

impl I14DomainApplicability {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 3] = [Self::Admitted, Self::OutOfDomain, Self::Indeterminate];
}

/// Whether the requested operation is implemented in the admitted capability set.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14OperationalSupport {
    /// The requested operation is supported by the admitted capability set.
    SupportedOperation = 0,
    /// The requested operation is outside the admitted capability set.
    UnsupportedOperation = 1,
}

impl I14OperationalSupport {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 2] = [Self::SupportedOperation, Self::UnsupportedOperation];
}

/// Schema, causal-event, and cross-axis validity of the terminal receipt.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14ReceiptValidity {
    /// The receipt is schema-valid and its axes are mutually consistent.
    WellFormedReceipt = 0,
    /// The receipt is malformed or contains a forbidden axis combination.
    MalformedReceipt = 1,
}

impl I14ReceiptValidity {
    /// Canonical enumeration order used by the version-1 truth table.
    pub const ALL: [Self; 2] = [Self::WellFormedReceipt, Self::MalformedReceipt];
}

/// The eight orthogonal I14 terminal receipt axes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct I14TerminalStatusV1 {
    /// Compute/supervisor disposition.
    pub execution: I14ExecutionDisposition,
    /// Scientific claim adjudication.
    pub claim: I14ClaimAdjudication,
    /// Required-evidence completeness.
    pub completeness: I14EvidenceCompleteness,
    /// Evidence/custody/checker integrity.
    pub integrity: I14EvidenceIntegrity,
    /// Submitted-input validity.
    pub input: I14InputValidity,
    /// Frozen-domain applicability.
    pub domain: I14DomainApplicability,
    /// Admitted-capability support.
    pub support: I14OperationalSupport,
    /// Receipt and cross-axis validity.
    pub receipt: I14ReceiptValidity,
}

impl I14TerminalStatusV1 {
    const fn tags(self) -> [u8; 8] {
        [
            self.execution as u8,
            self.claim as u8,
            self.completeness as u8,
            self.integrity as u8,
            self.input as u8,
            self.domain as u8,
            self.support as u8,
            self.receipt as u8,
        ]
    }
}

/// Fail-closed normalization and primary CLI projection of one terminal status.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14TerminalEvaluationV1 {
    /// Raw producer status retained for provenance and contradiction audit.
    pub raw: I14TerminalStatusV1,
    /// Authoritative normalized status; no input axis is discarded.
    pub normalized: I14TerminalStatusV1,
    /// Explicit normalization actions applied by the fail-closed evaluator.
    pub normalization: I14TerminalNormalizationV1,
    /// Deterministic primary CLI exit code selected by the frozen precedence.
    pub exit_code: u8,
}

/// Explicit, receipt-bound changes made while normalizing a raw terminal tuple.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14TerminalNormalizationV1 {
    /// Malformed input forced domain applicability to `Indeterminate`.
    pub domain_forced_indeterminate: bool,
    /// A raw cross-axis contradiction changed receipt validity to malformed.
    pub receipt_marked_malformed: bool,
}

/// Normalize and adjudicate one raw I14 terminal-status tuple.
///
/// Malformed input always normalizes domain applicability to
/// [`I14DomainApplicability::Indeterminate`]. Any forbidden combination marks
/// the receipt malformed; it never silently rewrites a scientific verdict.
#[must_use]
pub fn i14_evaluate_terminal_status_v1(mut status: I14TerminalStatusV1) -> I14TerminalEvaluationV1 {
    let raw = status;
    let conflicting_malformed_input_domain = status.input == I14InputValidity::MalformedInput
        && status.domain != I14DomainApplicability::Indeterminate;
    if status.input == I14InputValidity::MalformedInput {
        status.domain = I14DomainApplicability::Indeterminate;
    }

    let non_scientific_claim = matches!(
        status.claim,
        I14ClaimAdjudication::Unknown | I14ClaimAdjudication::Pending
    );
    let context_requires_non_scientific_claim = status.execution
        != I14ExecutionDisposition::Completed
        || status.completeness != I14EvidenceCompleteness::CompleteEvidence
        || status.integrity != I14EvidenceIntegrity::IntegrityVerified
        || status.input != I14InputValidity::WellFormedInput
        || status.domain != I14DomainApplicability::Admitted
        || status.support != I14OperationalSupport::SupportedOperation
        || status.receipt != I14ReceiptValidity::WellFormedReceipt;
    let pending_is_valid = status.claim != I14ClaimAdjudication::Pending
        || status.execution != I14ExecutionDisposition::Completed
        || status.completeness != I14EvidenceCompleteness::CompleteEvidence;

    if conflicting_malformed_input_domain
        || (context_requires_non_scientific_claim && !non_scientific_claim)
        || !pending_is_valid
    {
        status.receipt = I14ReceiptValidity::MalformedReceipt;
    }

    let exit_code = if status.receipt == I14ReceiptValidity::MalformedReceipt
        || status.input == I14InputValidity::MalformedInput
        || status.integrity == I14EvidenceIntegrity::IntegrityFailed
    {
        60
    } else if status.execution == I14ExecutionDisposition::InfrastructureFailed {
        70
    } else if status.execution == I14ExecutionDisposition::Cancelled {
        20
    } else if status.execution == I14ExecutionDisposition::TimedOut {
        21
    } else if status.execution == I14ExecutionDisposition::BudgetExhausted {
        22
    } else if status.claim == I14ClaimAdjudication::Failed {
        10
    } else if status.claim == I14ClaimAdjudication::Refuted {
        40
    } else if status.domain != I14DomainApplicability::Admitted
        || status.support == I14OperationalSupport::UnsupportedOperation
    {
        30
    } else if status.completeness != I14EvidenceCompleteness::CompleteEvidence {
        50
    } else if status.claim == I14ClaimAdjudication::Unknown {
        30
    } else {
        0
    };

    I14TerminalEvaluationV1 {
        raw,
        normalized: status,
        normalization: I14TerminalNormalizationV1 {
            domain_forced_indeterminate: raw.input == I14InputValidity::MalformedInput
                && raw.domain != I14DomainApplicability::Indeterminate,
            receipt_marked_malformed: raw.receipt == I14ReceiptValidity::WellFormedReceipt
                && status.receipt == I14ReceiptValidity::MalformedReceipt,
        },
        exit_code,
    }
}

/// One observed cancellation event paired with its request.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct I14CancellationObservationV1 {
    /// Coordinator-assigned logical event sequence.
    pub logical_sequence: u64,
    /// Receipt-bound monotonic telemetry timestamp in nanoseconds.
    pub monotonic_ns: u64,
    /// Logical tile that observed the request.
    pub observing_tile_id: u64,
    /// Latest completed logical boundary when observation occurred, or `None`
    /// when observation preceded the first completed boundary.
    pub latest_completed_boundary_ordinal: Option<u64>,
}

/// One cancellation request visible to terminal-boundary arbitration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14CancellationRequestV1 {
    /// Unique request identity within the execution scope.
    pub request_id: u64,
    /// Scope root cancelled by this request.
    pub scope_root: u64,
    /// Coordinator-assigned logical event sequence.
    pub logical_sequence: u64,
    /// Receipt-bound monotonic request timestamp in nanoseconds.
    pub requested_monotonic_ns: u64,
    /// Inclusive monotonic deadline for observing this request.
    pub observation_deadline_ns: u64,
    /// Observation, if it exists in the recorded event trace.
    pub observation: Option<I14CancellationObservationV1>,
}

/// One deterministic logical boundary considered for terminal selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct I14TerminalBoundaryV1<'a> {
    /// Logical boundary ordinal within the execution receipt.
    pub boundary_ordinal: u64,
    /// Coordinator sequence of the boundary-selection event.
    pub logical_sequence: u64,
    /// Recorded monotonic time at which this boundary was adjudicated.
    pub monotonic_ns: u64,
    /// Root-to-leaf scope ancestry, including the candidate scope itself.
    pub scope_ancestry: &'a [u64],
    /// Tile identities admitted to observe cancellation at this boundary.
    pub admitted_observer_tile_ids: &'a [u64],
    /// A same-boundary infrastructure failure is present.
    pub infrastructure_failed: bool,
    /// A same-boundary wall/deadline timeout independent of cancellation is present.
    pub timed_out: bool,
    /// The applicable total resource ceiling was reached.
    pub budget_exhausted: bool,
    /// All scheduled work completed normally.
    pub completed: bool,
}

/// Fail-closed refusal of malformed terminal-cause arbitration input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14TerminalCauseRefusalV1 {
    /// The request collection exceeds its admitted bound.
    TooManyRequests {
        /// Supplied request count.
        count: usize,
        /// Maximum admitted request count.
        cap: usize,
    },
    /// The candidate has no scope identity or ancestry.
    EmptyScopeAncestry,
    /// Scope ancestry exceeds its admitted bound.
    ScopeAncestryTooDeep {
        /// Supplied ancestry depth.
        depth: usize,
        /// Maximum admitted ancestry depth.
        cap: usize,
    },
    /// The root-to-leaf scope path repeats one identity and is cyclic.
    DuplicateScopeId {
        /// Repeated scope identity.
        scope_id: u64,
    },
    /// The admitted observer-tile catalog exceeds its bound.
    TooManyObserverTiles {
        /// Supplied observer-tile count.
        count: usize,
        /// Maximum admitted observer-tile count.
        cap: usize,
    },
    /// The observer-tile catalog repeats one tile identity.
    DuplicateObserverTileId {
        /// Duplicated observer-tile identity.
        observing_tile_id: u64,
    },
    /// Two requests reuse one identity.
    DuplicateRequestId {
        /// Duplicated request identity.
        request_id: u64,
    },
    /// Two causal events reuse one coordinator logical sequence.
    DuplicateLogicalSequence {
        /// Duplicated logical sequence.
        logical_sequence: u64,
    },
    /// A request deadline precedes its request timestamp.
    DeadlineBeforeRequest {
        /// Malformed request identity.
        request_id: u64,
    },
    /// An observation precedes its paired request in logical or monotonic order.
    ObservationBeforeRequest {
        /// Malformed request identity.
        request_id: u64,
    },
    /// Coordinator logical order and calibrated monotonic order disagree.
    NonMonotonicLogicalTimestamp {
        /// Earlier coordinator sequence.
        earlier_logical_sequence: u64,
        /// Timestamp carried by the earlier sequence.
        earlier_monotonic_ns: u64,
        /// Later coordinator sequence.
        later_logical_sequence: u64,
        /// Timestamp carried by the later sequence.
        later_monotonic_ns: u64,
    },
    /// An already-logical observation claims the candidate or a future boundary
    /// was already complete.
    ObservationBoundaryNotBeforeCandidate {
        /// Malformed request identity.
        request_id: u64,
    },
    /// A post-boundary observation claims that the candidate boundary had not
    /// yet completed.
    ObservationBoundaryBehindCandidate {
        /// Malformed request identity.
        request_id: u64,
        /// Latest completed boundary claimed by the observation.
        latest_completed_boundary_ordinal: Option<u64>,
        /// Candidate boundary that necessarily completed earlier.
        candidate_boundary_ordinal: u64,
    },
    /// An observation names a tile absent from the boundary's admitted catalog.
    UnknownObserverTile {
        /// Malformed request identity.
        request_id: u64,
        /// Unadmitted observing tile identity.
        observing_tile_id: u64,
    },
}

/// Deterministic result of scoped terminal-boundary cause arbitration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14TerminalBoundaryDecisionV1 {
    /// A terminal disposition was selected at this boundary.
    Selected {
        /// Selected disposition under the frozen cause priority.
        disposition: I14ExecutionDisposition,
        /// Winning cancellation request, when cancellation/deadline caused selection.
        request_id: Option<u64>,
    },
    /// A normal completion/budget candidate is held for an earlier scoped request.
    DeferredByCancellation {
        /// Deterministically earliest blocking request.
        request_id: u64,
    },
    /// No terminal cause or normal terminal candidate is present at this boundary.
    NotTerminal,
}

/// Semantic state of one cancellation request at the candidate boundary.
///
/// This state deliberately excludes raw clock values. A timing change that
/// crosses a request/deadline/boundary relation changes this state and thus the
/// canonical result; a timing change that preserves the relation affects only
/// the telemetry envelope.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14CancellationRequestStateV1 {
    /// The prior request targets no scope in the candidate ancestry.
    OutOfScope = 0,
    /// A prior in-scope request was observed no later than its deadline.
    ObservedWithinDeadline = 1,
    /// A prior in-scope request missed its observation deadline.
    MissedObservationDeadline = 2,
    /// A prior in-scope request is still awaiting an unexpired observation.
    PendingObservation = 3,
}

/// Clock-free canonical projection of one cancellation observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct I14CanonicalCancellationObservationV1 {
    /// Coordinator-assigned logical event sequence.
    pub logical_sequence: u64,
    /// Logical tile that observed the request.
    pub observing_tile_id: u64,
    /// Latest completed boundary, or `None` before the first completion.
    pub latest_completed_boundary_ordinal: Option<u64>,
}

/// Clock-free canonical projection of one cancellation request and its state.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct I14CanonicalCancellationRequestV1 {
    /// Unique request identity within the execution scope.
    pub request_id: u64,
    /// Scope root cancelled by the request.
    pub scope_root: u64,
    /// Coordinator-assigned request sequence.
    pub logical_sequence: u64,
    /// Semantic request state at the candidate boundary.
    pub state: I14CancellationRequestStateV1,
    /// Clock-free observation projection, when any observation was recorded.
    pub observation: Option<I14CanonicalCancellationObservationV1>,
}

/// Raw terminal-cause candidates bound into canonical result identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct I14TerminalCauseCandidatesV1 {
    /// Infrastructure failure candidate.
    pub infrastructure_failed: bool,
    /// Independent wall/deadline timeout candidate.
    pub timed_out: bool,
    /// Total resource ceiling candidate.
    pub budget_exhausted: bool,
    /// Normal completion candidate.
    pub completed: bool,
}

/// Inputs frozen into one canonical terminal-result identity.
#[derive(Clone, Copy, Debug)]
pub struct I14CanonicalTerminalResultInputV1<'a> {
    /// Candidate boundary and raw cause candidates.
    pub boundary: I14TerminalBoundaryV1<'a>,
    /// Recorded request trace; the canonical result takes its immutable
    /// logical cut strictly before `boundary.logical_sequence`.
    pub requests: &'a [I14CancellationRequestV1],
    /// Raw terminal-status tuple retained and fail-closed-normalized.
    pub terminal_status: I14TerminalStatusV1,
    /// Content identity of the scientific/operational result payload.
    pub semantic_payload_digest: ContentHash,
}

/// Fail-closed canonical terminal-result construction refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14CanonicalResultRefusalV1 {
    /// The immutable at-boundary trace is malformed.
    TerminalCause(I14TerminalCauseRefusalV1),
    /// This boundary has not selected a terminal disposition.
    BoundaryNotTerminal {
        /// Nonterminal decision returned by the cause selector.
        decision: I14TerminalBoundaryDecisionV1,
    },
    /// The terminal receipt disagrees with the selected execution cause.
    ExecutionDispositionMismatch {
        /// Disposition selected from the causal trace.
        selected: I14ExecutionDisposition,
        /// Disposition carried by the normalized terminal receipt.
        receipt: I14ExecutionDisposition,
    },
}

/// Validated, clock-free terminal result whose digest is G5-comparable.
///
/// Construction is available only through
/// [`i14_canonical_terminal_result_v1`], which invokes the same fail-closed
/// cause selector used for terminal adjudication.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct I14CanonicalTerminalResultV1 {
    boundary_ordinal: u64,
    boundary_logical_sequence: u64,
    scope_ancestry: Vec<u64>,
    admitted_observer_tile_ids: Vec<u64>,
    cause_candidates: I14TerminalCauseCandidatesV1,
    decision: I14TerminalBoundaryDecisionV1,
    requests: Vec<I14CanonicalCancellationRequestV1>,
    terminal_evaluation: I14TerminalEvaluationV1,
    semantic_payload_digest: ContentHash,
}

impl I14CanonicalTerminalResultV1 {
    /// Logical boundary ordinal selected for this result.
    #[must_use]
    pub const fn boundary_ordinal(&self) -> u64 {
        self.boundary_ordinal
    }

    /// Coordinator logical sequence of the boundary-selection event.
    #[must_use]
    pub const fn boundary_logical_sequence(&self) -> u64 {
        self.boundary_logical_sequence
    }

    /// Root-to-leaf acyclic scope path.
    #[must_use]
    pub fn scope_ancestry(&self) -> &[u64] {
        &self.scope_ancestry
    }

    /// Canonically sorted observer-tile catalog.
    #[must_use]
    pub fn admitted_observer_tile_ids(&self) -> &[u64] {
        &self.admitted_observer_tile_ids
    }

    /// Raw semantic cause candidates at the selected boundary.
    #[must_use]
    pub const fn cause_candidates(&self) -> I14TerminalCauseCandidatesV1 {
        self.cause_candidates
    }

    /// Deterministically selected terminal decision.
    #[must_use]
    pub const fn decision(&self) -> I14TerminalBoundaryDecisionV1 {
        self.decision
    }

    /// Canonically sequence-ordered request/observation trace.
    #[must_use]
    pub fn requests(&self) -> &[I14CanonicalCancellationRequestV1] {
        &self.requests
    }

    /// Raw/normalized terminal axes, normalization actions, and CLI projection.
    #[must_use]
    pub const fn terminal_evaluation(&self) -> I14TerminalEvaluationV1 {
        self.terminal_evaluation
    }

    /// Content identity of the scientific/operational result payload.
    #[must_use]
    pub const fn semantic_payload_digest(&self) -> ContentHash {
        self.semantic_payload_digest
    }

    /// Domain-separated clock-free canonical identity.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut payload = b"I14_CANONICAL_TERMINAL_RESULT_V1\0".to_vec();
        i14_push_u64(&mut payload, self.boundary_ordinal);
        i14_push_u64(&mut payload, self.boundary_logical_sequence);
        i14_push_u64_slice(&mut payload, &self.scope_ancestry);
        i14_push_u64_slice(&mut payload, &self.admitted_observer_tile_ids);
        payload.extend_from_slice(&[
            u8::from(self.cause_candidates.infrastructure_failed),
            u8::from(self.cause_candidates.timed_out),
            u8::from(self.cause_candidates.budget_exhausted),
            u8::from(self.cause_candidates.completed),
        ]);
        i14_push_terminal_decision(&mut payload, self.decision);
        i14_push_len(&mut payload, self.requests.len());
        for request in &self.requests {
            i14_push_u64(&mut payload, request.request_id);
            i14_push_u64(&mut payload, request.scope_root);
            i14_push_u64(&mut payload, request.logical_sequence);
            payload.push(request.state as u8);
            match request.observation {
                None => payload.push(0),
                Some(observation) => {
                    payload.push(1);
                    i14_push_u64(&mut payload, observation.logical_sequence);
                    i14_push_u64(&mut payload, observation.observing_tile_id);
                    i14_push_optional_u64(
                        &mut payload,
                        observation.latest_completed_boundary_ordinal,
                    );
                }
            }
        }
        payload.extend_from_slice(&self.terminal_evaluation.raw.tags());
        payload.extend_from_slice(&self.terminal_evaluation.normalized.tags());
        payload.extend_from_slice(&[
            u8::from(
                self.terminal_evaluation
                    .normalization
                    .domain_forced_indeterminate,
            ),
            u8::from(
                self.terminal_evaluation
                    .normalization
                    .receipt_marked_malformed,
            ),
            self.terminal_evaluation.exit_code,
        ]);
        payload.extend_from_slice(self.semantic_payload_digest.as_bytes());
        hash_domain(TERMINAL_RESULT_DOMAIN_V1, &payload)
    }
}

/// Kind of raw watchdog observation retained in the telemetry envelope.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14WatchdogObservationKindV1 {
    /// Local asynchronous watchdog poll.
    Poll = 0,
    /// Heartbeat from an admitted external child.
    ExternalHeartbeat = 1,
    /// Explicit deadline/expiry sample.
    DeadlineSample = 2,
}

/// One raw, noncanonical watchdog timing observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14WatchdogObservationV1 {
    /// Unique identity within the telemetry envelope.
    pub observation_id: u64,
    /// Observation kind.
    pub kind: I14WatchdogObservationKindV1,
    /// Calibrated monotonic timestamp in nanoseconds.
    pub monotonic_ns: u64,
}

/// Inputs whose raw timing/calibration fields form one telemetry envelope.
#[derive(Clone, Copy, Debug)]
pub struct I14TelemetryEnvelopeInputV1<'a> {
    /// Validated candidate boundary, including its raw monotonic time.
    pub boundary: I14TerminalBoundaryV1<'a>,
    /// Recorded requests/observations, including raw monotonic times/deadlines.
    pub requests: &'a [I14CancellationRequestV1],
    /// Raw terminal-status tuple bound into the canonical result.
    pub terminal_status: I14TerminalStatusV1,
    /// Content identity of the scientific/operational result payload.
    pub semantic_payload_digest: ContentHash,
    /// Raw watchdog observations; presentation order carries no identity.
    pub watchdog_observations: &'a [I14WatchdogObservationV1],
    /// Content identity of the clock-calibration artifact.
    pub clock_calibration_artifact: ContentHash,
}

/// Fail-closed telemetry-envelope construction refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14TelemetryEnvelopeRefusalV1 {
    /// The complete raw trace, including post-boundary events, is malformed.
    RecordedTrace(I14TerminalCauseRefusalV1),
    /// The immutable at-boundary terminal result cannot be constructed.
    CanonicalResult(I14CanonicalResultRefusalV1),
    /// The watchdog collection exceeds its admitted bound.
    TooManyWatchdogObservations {
        /// Supplied watchdog observation count.
        count: usize,
        /// Maximum admitted observation count.
        cap: usize,
    },
    /// Two watchdog observations reuse one identity.
    DuplicateWatchdogObservationId {
        /// Reused observation identity.
        observation_id: u64,
    },
}

fn i14_request_precedes(left: &I14CancellationRequestV1, right: &I14CancellationRequestV1) -> bool {
    left.logical_sequence < right.logical_sequence
}

fn i14_keep_earliest_request<'a>(
    slot: &mut Option<&'a I14CancellationRequestV1>,
    candidate: &'a I14CancellationRequestV1,
) {
    if slot.is_none_or(|current| i14_request_precedes(candidate, current)) {
        *slot = Some(candidate);
    }
}

/// Select a terminal cause at one caller-supplied scoped boundary.
///
/// This V1 operation is retained for local arbitration and legacy artifact
/// interpretation only. It proves neither that the supplied boundary is the
/// first terminal-eligible boundary nor that tier cancellation/drain SLOs were
/// met, and therefore has no promotion authority. New campaigns must use
/// [`i14_select_first_terminal_boundary_v2`] and the V2 canonical constructor.
///
/// Only requests with a logical sequence before the candidate boundary and a
/// scope root in `boundary.scope_ancestry` participate. Multiple requests are
/// canonicalized before validation and resolved by their globally unique
/// logical sequence. Observer-tile catalogs are likewise canonicalized before
/// duplicate diagnosis, so valid decisions and malformed refusals are input-
/// order invariant. A missed relevant observation deadline outranks an observed
/// cancellation; an unobserved, unexpired request defers only a normal
/// `Completed`/`BudgetExhausted` candidate. Every participating pre-boundary
/// observation must name an observer admitted at this boundary. A retained
/// later observation uses its later boundary's observer catalog, but must
/// acknowledge this candidate in its latest-completed boundary ordinal. The
/// full order is infrastructure, timeout, cancellation, budget exhaustion,
/// completion.
#[allow(clippy::too_many_lines)]
pub fn i14_select_terminal_boundary_v1(
    boundary: I14TerminalBoundaryV1<'_>,
    requests: &[I14CancellationRequestV1],
) -> Result<I14TerminalBoundaryDecisionV1, I14TerminalCauseRefusalV1> {
    if requests.len() > I14_MAX_CANCELLATION_REQUESTS_V1 {
        return Err(I14TerminalCauseRefusalV1::TooManyRequests {
            count: requests.len(),
            cap: I14_MAX_CANCELLATION_REQUESTS_V1,
        });
    }
    if boundary.scope_ancestry.is_empty() {
        return Err(I14TerminalCauseRefusalV1::EmptyScopeAncestry);
    }
    if boundary.scope_ancestry.len() > I14_MAX_SCOPE_ANCESTRY_V1 {
        return Err(I14TerminalCauseRefusalV1::ScopeAncestryTooDeep {
            depth: boundary.scope_ancestry.len(),
            cap: I14_MAX_SCOPE_ANCESTRY_V1,
        });
    }
    let mut scope_ids = BTreeSet::new();
    let mut duplicate_scope_ids = BTreeSet::new();
    for scope_id in boundary.scope_ancestry {
        if !scope_ids.insert(*scope_id) {
            duplicate_scope_ids.insert(*scope_id);
        }
    }
    if let Some(scope_id) = duplicate_scope_ids.first().copied() {
        return Err(I14TerminalCauseRefusalV1::DuplicateScopeId { scope_id });
    }
    if boundary.admitted_observer_tile_ids.len() > I14_MAX_OBSERVER_TILES_V1 {
        return Err(I14TerminalCauseRefusalV1::TooManyObserverTiles {
            count: boundary.admitted_observer_tile_ids.len(),
            cap: I14_MAX_OBSERVER_TILES_V1,
        });
    }
    let mut ordered_observer_tile_ids = boundary.admitted_observer_tile_ids.to_vec();
    ordered_observer_tile_ids.sort_unstable();
    let mut observer_tile_ids = BTreeSet::new();
    for observing_tile_id in ordered_observer_tile_ids {
        if !observer_tile_ids.insert(observing_tile_id) {
            return Err(I14TerminalCauseRefusalV1::DuplicateObserverTileId { observing_tile_id });
        }
    }
    let mut ordered_requests = requests.iter().collect::<Vec<_>>();
    ordered_requests.sort_unstable_by_key(|request| {
        (
            request.logical_sequence,
            request.request_id,
            request.scope_root,
            request.requested_monotonic_ns,
            request.observation_deadline_ns,
            request.observation,
        )
    });
    let mut request_ids = BTreeSet::new();
    let mut logical_sequences = BTreeSet::from([boundary.logical_sequence]);
    let mut logical_timestamps = vec![(boundary.logical_sequence, boundary.monotonic_ns)];
    for request in &ordered_requests {
        let request = *request;
        if !request_ids.insert(request.request_id) {
            return Err(I14TerminalCauseRefusalV1::DuplicateRequestId {
                request_id: request.request_id,
            });
        }
        if !logical_sequences.insert(request.logical_sequence) {
            return Err(I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
                logical_sequence: request.logical_sequence,
            });
        }
        logical_timestamps.push((request.logical_sequence, request.requested_monotonic_ns));
        if request.observation_deadline_ns < request.requested_monotonic_ns {
            return Err(I14TerminalCauseRefusalV1::DeadlineBeforeRequest {
                request_id: request.request_id,
            });
        }
        if let Some(observation) = request.observation {
            if observation.logical_sequence <= request.logical_sequence
                || observation.monotonic_ns < request.requested_monotonic_ns
            {
                return Err(I14TerminalCauseRefusalV1::ObservationBeforeRequest {
                    request_id: request.request_id,
                });
            }
            if !logical_sequences.insert(observation.logical_sequence) {
                return Err(I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
                    logical_sequence: observation.logical_sequence,
                });
            }
            logical_timestamps.push((observation.logical_sequence, observation.monotonic_ns));
            if observation.logical_sequence < boundary.logical_sequence {
                if scope_ids.contains(&request.scope_root)
                    && !observer_tile_ids.contains(&observation.observing_tile_id)
                {
                    return Err(I14TerminalCauseRefusalV1::UnknownObserverTile {
                        request_id: request.request_id,
                        observing_tile_id: observation.observing_tile_id,
                    });
                }
                if observation
                    .latest_completed_boundary_ordinal
                    .is_some_and(|ordinal| ordinal >= boundary.boundary_ordinal)
                {
                    return Err(
                        I14TerminalCauseRefusalV1::ObservationBoundaryNotBeforeCandidate {
                            request_id: request.request_id,
                        },
                    );
                }
            } else if observation.logical_sequence > boundary.logical_sequence
                && observation
                    .latest_completed_boundary_ordinal
                    .is_none_or(|ordinal| ordinal < boundary.boundary_ordinal)
            {
                return Err(
                    I14TerminalCauseRefusalV1::ObservationBoundaryBehindCandidate {
                        request_id: request.request_id,
                        latest_completed_boundary_ordinal: observation
                            .latest_completed_boundary_ordinal,
                        candidate_boundary_ordinal: boundary.boundary_ordinal,
                    },
                );
            }
        }
    }
    logical_timestamps.sort_unstable_by_key(|(logical_sequence, _)| *logical_sequence);
    for pair in logical_timestamps.windows(2) {
        let (earlier_logical_sequence, earlier_monotonic_ns) = pair[0];
        let (later_logical_sequence, later_monotonic_ns) = pair[1];
        if earlier_monotonic_ns > later_monotonic_ns {
            return Err(I14TerminalCauseRefusalV1::NonMonotonicLogicalTimestamp {
                earlier_logical_sequence,
                earlier_monotonic_ns,
                later_logical_sequence,
                later_monotonic_ns,
            });
        }
    }

    let mut missed = None;
    let mut observed = None;
    let mut pending = None;
    for request in ordered_requests {
        if request.logical_sequence >= boundary.logical_sequence
            || !scope_ids.contains(&request.scope_root)
        {
            continue;
        }
        match request.observation {
            Some(observation) if observation.logical_sequence < boundary.logical_sequence => {
                if observation.monotonic_ns <= request.observation_deadline_ns {
                    i14_keep_earliest_request(&mut observed, request);
                } else {
                    i14_keep_earliest_request(&mut missed, request);
                }
            }
            _ if boundary.monotonic_ns > request.observation_deadline_ns => {
                i14_keep_earliest_request(&mut missed, request);
            }
            _ => i14_keep_earliest_request(&mut pending, request),
        }
    }

    if boundary.infrastructure_failed {
        return Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::InfrastructureFailed,
            request_id: None,
        });
    }
    if boundary.timed_out || missed.is_some() {
        return Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::TimedOut,
            request_id: missed.map(|request| request.request_id),
        });
    }
    if let Some(request) = observed {
        return Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Cancelled,
            request_id: Some(request.request_id),
        });
    }
    if (boundary.budget_exhausted || boundary.completed)
        && let Some(request) = pending
    {
        return Ok(I14TerminalBoundaryDecisionV1::DeferredByCancellation {
            request_id: request.request_id,
        });
    }
    if boundary.budget_exhausted {
        return Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::BudgetExhausted,
            request_id: None,
        });
    }
    if boundary.completed {
        return Ok(I14TerminalBoundaryDecisionV1::Selected {
            disposition: I14ExecutionDisposition::Completed,
            request_id: None,
        });
    }
    Ok(I14TerminalBoundaryDecisionV1::NotTerminal)
}

fn i14_push_u64(payload: &mut Vec<u8>, value: u64) {
    payload.extend_from_slice(&value.to_le_bytes());
}

fn i14_push_len(payload: &mut Vec<u8>, length: usize) {
    let length = u64::try_from(length).expect("I14 admitted collection length fits u64");
    i14_push_u64(payload, length);
}

fn i14_push_u64_slice(payload: &mut Vec<u8>, values: &[u64]) {
    i14_push_len(payload, values.len());
    for value in values {
        i14_push_u64(payload, *value);
    }
}

fn i14_push_optional_u64(payload: &mut Vec<u8>, value: Option<u64>) {
    match value {
        None => payload.push(0),
        Some(value) => {
            payload.push(1);
            i14_push_u64(payload, value);
        }
    }
}

fn i14_push_drain_trigger(payload: &mut Vec<u8>, trigger: I14DrainTriggerV2) {
    match trigger {
        I14DrainTriggerV2::NonCancellationDrain => payload.push(0),
        I14DrainTriggerV2::CancellationObserved { request_id } => {
            payload.push(1);
            i14_push_u64(payload, request_id);
        }
        I14DrainTriggerV2::ObservationTimeoutDrain { request_id } => {
            payload.push(2);
            i14_push_u64(payload, request_id);
        }
        I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence,
        } => {
            payload.push(3);
            i14_push_u64(payload, onset_logical_sequence);
        }
    }
}

/// Return the normative exact wire bytes for one V2 drain-trigger union.
///
/// NonCancellationDrain is 00; the remaining variants are tag 01, 02,
/// or 03 followed by one unframed little-endian u64. The containing
/// canonical-result grammar supplies position and framing.
#[must_use]
pub fn i14_drain_trigger_encoding_v2(trigger: I14DrainTriggerV2) -> Vec<u8> {
    let mut payload = Vec::with_capacity(9);
    i14_push_drain_trigger(&mut payload, trigger);
    payload
}

fn i14_push_infrastructure_failure_onset_v2(
    payload: &mut Vec<u8>,
    onset: Option<I14InfrastructureFailureOnsetV2>,
) {
    i14_push_canonical_infrastructure_failure_onset_v2(
        payload,
        onset.map(|onset| onset.event.logical_sequence),
        onset.map(|onset| onset.source),
        onset.map(|onset| onset.verification_receipt_digest),
    );
}

fn i14_push_canonical_infrastructure_failure_onset_v2(
    payload: &mut Vec<u8>,
    logical_sequence: Option<u64>,
    source: Option<I14InfrastructureFailureSourceV2>,
    verification_receipt_digest: Option<ContentHash>,
) {
    match (logical_sequence, source, verification_receipt_digest) {
        (None, None, None) => payload.push(0),
        (Some(logical_sequence), Some(source), Some(verification_receipt_digest)) => {
            payload.push(1);
            i14_push_u64(payload, logical_sequence);
            payload.push(source.wire_tag());
            payload.extend_from_slice(verification_receipt_digest.as_bytes());
        }
        _ => unreachable!("validated canonical infrastructure onset fields agree on presence"),
    }
}

/// Return canonical onset-witness bytes, excluding calibrated telemetry time.
#[must_use]
pub fn i14_infrastructure_failure_onset_encoding_v2(
    onset: Option<I14InfrastructureFailureOnsetV2>,
) -> Vec<u8> {
    let mut payload = Vec::with_capacity(42);
    i14_push_infrastructure_failure_onset_v2(&mut payload, onset);
    payload
}

fn i14_push_terminal_decision(payload: &mut Vec<u8>, decision: I14TerminalBoundaryDecisionV1) {
    match decision {
        I14TerminalBoundaryDecisionV1::Selected {
            disposition,
            request_id,
        } => {
            payload.extend_from_slice(&[0, disposition as u8]);
            i14_push_optional_u64(payload, request_id);
        }
        I14TerminalBoundaryDecisionV1::DeferredByCancellation { request_id } => {
            payload.push(1);
            i14_push_u64(payload, request_id);
        }
        I14TerminalBoundaryDecisionV1::NotTerminal => payload.push(2),
    }
}

fn i14_request_state_v1(
    boundary: I14TerminalBoundaryV1<'_>,
    scope_ids: &BTreeSet<u64>,
    request: I14CancellationRequestV1,
) -> I14CancellationRequestStateV1 {
    debug_assert!(request.logical_sequence < boundary.logical_sequence);
    if !scope_ids.contains(&request.scope_root) {
        return I14CancellationRequestStateV1::OutOfScope;
    }
    match request.observation {
        Some(observation) if observation.logical_sequence < boundary.logical_sequence => {
            if observation.monotonic_ns <= request.observation_deadline_ns {
                I14CancellationRequestStateV1::ObservedWithinDeadline
            } else {
                I14CancellationRequestStateV1::MissedObservationDeadline
            }
        }
        _ if boundary.monotonic_ns > request.observation_deadline_ns => {
            I14CancellationRequestStateV1::MissedObservationDeadline
        }
        _ => I14CancellationRequestStateV1::PendingObservation,
    }
}

/// Validate and construct a legacy single-boundary canonical terminal result.
///
/// This V1 layout is frozen for old-ledger readability and has no
/// first-terminal or cancellation-lifecycle promotion authority. New campaign
/// receipts must use [`i14_canonical_terminal_result_v2`].
///
/// Request and observer-catalog presentation order are normalized. Events at
/// the candidate's logical sequence are first included in selector validation
/// so they fail as identity collisions; only the strict `<` cut enters a
/// successful canonical result. Raw monotonic timestamps, deadlines, watchdog
/// arrivals, and clock calibration are excluded, but every timing-derived
/// semantic request state and selected cause is retained.
pub fn i14_canonical_terminal_result_v1(
    input: I14CanonicalTerminalResultInputV1<'_>,
) -> Result<I14CanonicalTerminalResultV1, I14CanonicalResultRefusalV1> {
    let boundary = input.boundary;
    if input.requests.len() > I14_MAX_CANCELLATION_REQUESTS_V1 {
        return Err(I14CanonicalResultRefusalV1::TerminalCause(
            I14TerminalCauseRefusalV1::TooManyRequests {
                count: input.requests.len(),
                cap: I14_MAX_CANCELLATION_REQUESTS_V1,
            },
        ));
    }
    let validation_requests = input
        .requests
        .iter()
        .filter(|request| request.logical_sequence <= boundary.logical_sequence)
        .map(|request| I14CancellationRequestV1 {
            observation: request
                .observation
                .filter(|observation| observation.logical_sequence <= boundary.logical_sequence),
            ..*request
        })
        .collect::<Vec<_>>();
    let decision = i14_select_terminal_boundary_v1(boundary, &validation_requests)
        .map_err(I14CanonicalResultRefusalV1::TerminalCause)?;
    let I14TerminalBoundaryDecisionV1::Selected {
        disposition: selected,
        ..
    } = decision
    else {
        return Err(I14CanonicalResultRefusalV1::BoundaryNotTerminal { decision });
    };
    let terminal_evaluation = i14_evaluate_terminal_status_v1(input.terminal_status);
    if terminal_evaluation.normalized.execution != selected {
        return Err(I14CanonicalResultRefusalV1::ExecutionDispositionMismatch {
            selected,
            receipt: terminal_evaluation.normalized.execution,
        });
    }
    let scope_ids = boundary
        .scope_ancestry
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let mut admitted_observer_tile_ids = boundary.admitted_observer_tile_ids.to_vec();
    admitted_observer_tile_ids.sort_unstable();
    let mut ordered_requests = validation_requests
        .into_iter()
        .filter(|request| request.logical_sequence < boundary.logical_sequence)
        .map(|request| I14CancellationRequestV1 {
            observation: request
                .observation
                .filter(|observation| observation.logical_sequence < boundary.logical_sequence),
            ..request
        })
        .collect::<Vec<_>>();
    ordered_requests.sort_unstable_by_key(|request| (request.logical_sequence, request.request_id));
    let requests = ordered_requests
        .into_iter()
        .map(|request| I14CanonicalCancellationRequestV1 {
            request_id: request.request_id,
            scope_root: request.scope_root,
            logical_sequence: request.logical_sequence,
            state: i14_request_state_v1(boundary, &scope_ids, request),
            observation: request.observation.map(|observation| {
                I14CanonicalCancellationObservationV1 {
                    logical_sequence: observation.logical_sequence,
                    observing_tile_id: observation.observing_tile_id,
                    latest_completed_boundary_ordinal: observation
                        .latest_completed_boundary_ordinal,
                }
            }),
        })
        .collect();
    Ok(I14CanonicalTerminalResultV1 {
        boundary_ordinal: boundary.boundary_ordinal,
        boundary_logical_sequence: boundary.logical_sequence,
        scope_ancestry: boundary.scope_ancestry.to_vec(),
        admitted_observer_tile_ids,
        cause_candidates: I14TerminalCauseCandidatesV1 {
            infrastructure_failed: boundary.infrastructure_failed,
            timed_out: boundary.timed_out,
            budget_exhausted: boundary.budget_exhausted,
            completed: boundary.completed,
        },
        decision,
        requests,
        terminal_evaluation,
        semantic_payload_digest: input.semantic_payload_digest,
    })
}

/// Validate and hash the clock-free canonical terminal result.
pub fn i14_canonical_terminal_result_digest_v1(
    input: I14CanonicalTerminalResultInputV1<'_>,
) -> Result<ContentHash, I14CanonicalResultRefusalV1> {
    Ok(i14_canonical_terminal_result_v1(input)?.digest())
}

/// Validate and hash a legacy single-boundary noncanonical telemetry envelope.
///
/// This V1 layout does not prove first-terminal selection or drain/finalization
/// SLOs. New campaign telemetry must use
/// [`i14_telemetry_envelope_digest_v2`].
///
/// The envelope binds the canonical result digest plus every raw request,
/// observation, boundary, watchdog, and clock-calibration timing input. Request
/// and watchdog presentation order carry no identity.
pub fn i14_telemetry_envelope_digest_v1(
    input: I14TelemetryEnvelopeInputV1<'_>,
) -> Result<ContentHash, I14TelemetryEnvelopeRefusalV1> {
    if input.watchdog_observations.len() > I14_MAX_WATCHDOG_OBSERVATIONS_V1 {
        return Err(I14TelemetryEnvelopeRefusalV1::TooManyWatchdogObservations {
            count: input.watchdog_observations.len(),
            cap: I14_MAX_WATCHDOG_OBSERVATIONS_V1,
        });
    }
    i14_select_terminal_boundary_v1(input.boundary, input.requests)
        .map_err(I14TelemetryEnvelopeRefusalV1::RecordedTrace)?;
    let canonical_result = i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
        boundary: input.boundary,
        requests: input.requests,
        terminal_status: input.terminal_status,
        semantic_payload_digest: input.semantic_payload_digest,
    })
    .map_err(I14TelemetryEnvelopeRefusalV1::CanonicalResult)?;
    let mut watchdog_observations = input.watchdog_observations.to_vec();
    watchdog_observations
        .sort_unstable_by_key(|observation| (observation.observation_id, observation.kind));
    for pair in watchdog_observations.windows(2) {
        if pair[0].observation_id == pair[1].observation_id {
            return Err(
                I14TelemetryEnvelopeRefusalV1::DuplicateWatchdogObservationId {
                    observation_id: pair[0].observation_id,
                },
            );
        }
    }

    let mut ordered_requests = input.requests.to_vec();
    ordered_requests.sort_unstable_by_key(|request| (request.logical_sequence, request.request_id));
    let mut payload = b"I14_NONCANONICAL_TELEMETRY_ENVELOPE_V1\0".to_vec();
    payload.extend_from_slice(canonical_result.digest().as_bytes());
    i14_push_u64(&mut payload, input.boundary.monotonic_ns);
    i14_push_len(&mut payload, ordered_requests.len());
    for request in ordered_requests {
        i14_push_u64(&mut payload, request.request_id);
        i14_push_u64(&mut payload, request.scope_root);
        i14_push_u64(&mut payload, request.logical_sequence);
        i14_push_u64(&mut payload, request.requested_monotonic_ns);
        i14_push_u64(&mut payload, request.observation_deadline_ns);
        match request.observation {
            None => payload.push(0),
            Some(observation) => {
                payload.push(1);
                i14_push_u64(&mut payload, observation.logical_sequence);
                i14_push_u64(&mut payload, observation.monotonic_ns);
                i14_push_u64(&mut payload, observation.observing_tile_id);
                i14_push_optional_u64(&mut payload, observation.latest_completed_boundary_ordinal);
            }
        }
    }
    i14_push_len(&mut payload, watchdog_observations.len());
    for observation in watchdog_observations {
        i14_push_u64(&mut payload, observation.observation_id);
        payload.push(observation.kind as u8);
        i14_push_u64(&mut payload, observation.monotonic_ns);
    }
    payload.extend_from_slice(input.clock_calibration_artifact.as_bytes());
    Ok(hash_domain(TELEMETRY_ENVELOPE_DOMAIN_V1, &payload))
}

/// Unit carried by the single total-resource ledger governed by a V2 card.
///
/// The unit is part of canonical identity; consumers must not compare or add
/// ceilings expressed on different axes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14TotalResourceUnitV2 {
    /// Native harness-graph or event-trace records.
    GraphTraceRecords = 0,
    /// Field unknowns or quadrature points.
    FieldQuadratureUnknowns = 1,
    /// Search-tree nodes.
    SearchTreeNodes = 2,
    /// Formal declarations admitted to a theorem campaign.
    FormalDeclarations = 3,
}

/// Authoritative cancellation-contract tier for an I14 Core or Max campaign.
///
/// Smoke campaigns intentionally have no promotion-authoritative terminal
/// trace: they may exercise the V1 local arbiter, but cannot construct a V2
/// cancellation card or V2 terminal result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum I14CancellationTierV2 {
    /// Core campaign contract.
    Core = 0,
    /// Max campaign contract.
    Max = 1,
    /// Max theorem/falsifier campaign with the separately frozen 24-hour cap.
    MaxTheoremFalsifier = 2,
}

impl I14CancellationTierV2 {
    /// Maximum observer tiles and in-flight child scopes.
    #[must_use]
    pub const fn child_and_observer_cap(self) -> usize {
        match self {
            Self::Core => 32,
            Self::Max | Self::MaxTheoremFalsifier => 128,
        }
    }

    /// Required request-to-observation deadline in nanoseconds.
    #[must_use]
    pub const fn request_to_observation_ns(self) -> u64 {
        match self {
            Self::Core => 250_000_000,
            Self::Max | Self::MaxTheoremFalsifier => 1_000_000_000,
        }
    }

    /// Maximum independent watchdog interval in nanoseconds.
    #[must_use]
    pub const fn watchdog_quantum_ns(self) -> u64 {
        match self {
            Self::Core => 25_000_000,
            Self::Max | Self::MaxTheoremFalsifier => 100_000_000,
        }
    }

    /// Maximum drain-trigger-to-drained interval in nanoseconds.
    #[must_use]
    pub const fn trigger_to_drained_ns(self) -> u64 {
        match self {
            Self::Core => 2_000_000_000,
            Self::Max | Self::MaxTheoremFalsifier => 8_000_000_000,
        }
    }

    /// Maximum drained-to-finalized interval in nanoseconds.
    #[must_use]
    pub const fn drained_to_finalized_ns(self) -> u64 {
        match self {
            Self::Core => 2_000_000_000,
            Self::Max | Self::MaxTheoremFalsifier => 8_000_000_000,
        }
    }

    /// Maximum admitted campaign wall budget in nanoseconds.
    #[must_use]
    pub const fn campaign_wall_budget_cap_ns(self) -> u64 {
        match self {
            Self::Core => 5_400_000_000_000,
            Self::Max => 64_800_000_000_000,
            Self::MaxTheoremFalsifier => 86_400_000_000_000,
        }
    }

    /// Maximum admitted logical-memory ceiling in bytes.
    #[must_use]
    pub const fn logical_memory_cap_bytes(self) -> u64 {
        match self {
            Self::Core => 32 * 1_073_741_824,
            Self::Max | Self::MaxTheoremFalsifier => 96 * 1_073_741_824,
        }
    }
}

impl I14TotalResourceUnitV2 {
    /// Frozen total-count cap for this work kind and tier.
    #[must_use]
    pub const fn total_ceiling_cap(self, tier: I14CancellationTierV2) -> u64 {
        let core = match self {
            Self::GraphTraceRecords => 4_096,
            Self::FieldQuadratureUnknowns => 16_384,
            Self::SearchTreeNodes => 1_024,
            Self::FormalDeclarations => 256,
        };
        match tier {
            I14CancellationTierV2::Core => core,
            I14CancellationTierV2::Max | I14CancellationTierV2::MaxTheoremFalsifier => core * 4,
        }
    }

    /// Frozen maximum logical poll-tile item count for this work kind and tier.
    #[must_use]
    pub const fn tile_quantum_cap(self, tier: I14CancellationTierV2) -> u64 {
        let core = match self {
            Self::GraphTraceRecords => 64,
            Self::FieldQuadratureUnknowns => 256,
            Self::SearchTreeNodes => 16,
            Self::FormalDeclarations => 4,
        };
        match tier {
            I14CancellationTierV2::Core => core,
            I14CancellationTierV2::Max | I14CancellationTierV2::MaxTheoremFalsifier => core * 4,
        }
    }
}

/// Raw fields proposed for one schema-authoritative cancellation card.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14CancellationCardInputV2 {
    /// Core, ordinary Max, or theorem/falsifier Max authority tier.
    pub tier: I14CancellationTierV2,
    /// Content identity of the semantic work unit governed by this card.
    pub semantic_work_unit_digest: ContentHash,
    /// Admitted campaign wall budget in nanoseconds.
    pub campaign_wall_budget_ns: u64,
    /// Hard allocator-governed logical-memory ceiling in bytes.
    pub logical_memory_ceiling_bytes: u64,
    /// Unit of the independently ledgered total-resource ceiling.
    pub total_resource_unit: I14TotalResourceUnitV2,
    /// Nonzero total-resource ceiling in `total_resource_unit`.
    pub total_resource_ceiling: u64,
    /// Largest indivisible resource quantum in the same unit.
    pub maximum_resource_tile_quantum: u64,
    /// Identity of the frozen budget definition and metering authority.
    pub resource_budget_authority_digest: ContentHash,
    /// Identity of the deterministic logical work partition.
    pub logical_partition_spec_digest: ContentHash,
    /// Host/toolchain/scheduler fingerprint on which response bounds were shown.
    pub execution_environment_digest: ContentHash,
    /// Demonstrated worst-case request-to-logical-boundary response bound.
    pub logical_tile_response_bound_ns: u64,
    /// Demonstrated worst-case indivisible-item response bound.
    pub indivisible_item_response_bound_ns: u64,
    /// Demonstrated maximum heartbeat interval for external children.
    pub external_heartbeat_bound_ns: u64,
    /// Exact external-child catalog, or `None` when no external child is
    /// admitted and the heartbeat bound must be zero.
    pub external_child_catalog_digest: Option<ContentHash>,
}

/// Fail-closed cancellation-card admission refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14CancellationCardRefusalV2 {
    /// The admitted campaign wall budget is zero.
    ZeroCampaignWallBudget,
    /// The hard logical-memory ceiling is zero.
    ZeroLogicalMemoryCeiling,
    /// The hard logical-memory ceiling exceeds the tier envelope.
    LogicalMemoryCeilingExceedsTier {
        /// Proposed byte ceiling.
        declared_bytes: u64,
        /// Tier maximum in bytes.
        cap_bytes: u64,
    },
    /// The total-resource ceiling is zero and therefore cannot define a gate.
    ZeroTotalResourceCeiling,
    /// The count-defined indivisible resource quantum is zero.
    ZeroResourceTileQuantum,
    /// One indivisible resource quantum is wider than the entire ceiling.
    ResourceTileQuantumExceedsCeiling {
        /// Proposed largest indivisible quantum.
        declared: u64,
        /// Proposed total-resource ceiling.
        ceiling: u64,
    },
    /// The work-kind total ceiling exceeds its frozen tier cap.
    TotalResourceCeilingExceedsTier {
        /// Proposed ceiling.
        declared: u64,
        /// Work-kind/tier maximum.
        cap: u64,
    },
    /// The work-kind tile quantum exceeds its frozen tier cap.
    ResourceTileQuantumExceedsTier {
        /// Proposed quantum.
        declared: u64,
        /// Work-kind/tier maximum.
        cap: u64,
    },
    /// The demonstrated logical-tile response bound is zero.
    ZeroLogicalTileResponseBound,
    /// The demonstrated indivisible-item response bound is zero.
    ZeroIndivisibleItemResponseBound,
    /// External-child catalog and heartbeat-bound presence disagree.
    ExternalChildPolicyMismatch,
    /// The campaign wall budget exceeds the tier envelope.
    CampaignWallBudgetExceedsTier {
        /// Proposed budget.
        declared_ns: u64,
        /// Tier maximum.
        cap_ns: u64,
    },
    /// A logical tile cannot meet the tier request-observation contract.
    LogicalTileResponseExceedsTier {
        /// Demonstrated bound.
        declared_ns: u64,
        /// Tier maximum.
        cap_ns: u64,
    },
    /// One indivisible item cannot meet the independent watchdog contract.
    IndivisibleItemResponseExceedsTier {
        /// Demonstrated bound.
        declared_ns: u64,
        /// Tier maximum.
        cap_ns: u64,
    },
    /// An external child cannot meet the independent heartbeat contract.
    ExternalHeartbeatExceedsTier {
        /// Demonstrated bound.
        declared_ns: u64,
        /// Tier maximum.
        cap_ns: u64,
    },
}

/// Validated, immutable cancellation contract used by V2 schema paths.
///
/// Validation proves internal bounds and canonical identity, not that a caller
/// actually measured the claimed response envelope. Promotion additionally
/// requires HELM/ledger to authenticate the verification receipt whose identity
/// is bound by a terminal trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14CancellationCardV2 {
    tier: I14CancellationTierV2,
    semantic_work_unit_digest: ContentHash,
    campaign_wall_budget_ns: u64,
    logical_memory_ceiling_bytes: u64,
    total_resource_unit: I14TotalResourceUnitV2,
    total_resource_ceiling: u64,
    maximum_resource_tile_quantum: u64,
    resource_budget_authority_digest: ContentHash,
    logical_partition_spec_digest: ContentHash,
    execution_environment_digest: ContentHash,
    logical_tile_response_bound_ns: u64,
    indivisible_item_response_bound_ns: u64,
    external_heartbeat_bound_ns: u64,
    external_child_catalog_digest: Option<ContentHash>,
}

impl I14CancellationCardV2 {
    /// Core, ordinary Max, or theorem/falsifier Max authority tier.
    #[must_use]
    pub const fn tier(self) -> I14CancellationTierV2 {
        self.tier
    }

    /// Content identity of the governed semantic work unit.
    #[must_use]
    pub const fn semantic_work_unit_digest(self) -> ContentHash {
        self.semantic_work_unit_digest
    }

    /// Admitted campaign wall budget in nanoseconds.
    #[must_use]
    pub const fn campaign_wall_budget_ns(self) -> u64 {
        self.campaign_wall_budget_ns
    }

    /// Hard allocator-governed logical-memory ceiling in bytes.
    #[must_use]
    pub const fn logical_memory_ceiling_bytes(self) -> u64 {
        self.logical_memory_ceiling_bytes
    }

    /// Unit of the independently ledgered total-resource ceiling.
    #[must_use]
    pub const fn total_resource_unit(self) -> I14TotalResourceUnitV2 {
        self.total_resource_unit
    }

    /// Frozen total-resource ceiling.
    #[must_use]
    pub const fn total_resource_ceiling(self) -> u64 {
        self.total_resource_ceiling
    }

    /// Largest admitted indivisible resource quantum.
    #[must_use]
    pub const fn maximum_resource_tile_quantum(self) -> u64 {
        self.maximum_resource_tile_quantum
    }

    /// Identity of the budget definition and metering authority.
    #[must_use]
    pub const fn resource_budget_authority_digest(self) -> ContentHash {
        self.resource_budget_authority_digest
    }

    /// Deterministic logical-partition specification identity.
    #[must_use]
    pub const fn logical_partition_spec_digest(self) -> ContentHash {
        self.logical_partition_spec_digest
    }

    /// Host/toolchain/scheduler response-evidence identity.
    #[must_use]
    pub const fn execution_environment_digest(self) -> ContentHash {
        self.execution_environment_digest
    }

    /// Demonstrated request-to-boundary response bound in nanoseconds.
    #[must_use]
    pub const fn logical_tile_response_bound_ns(self) -> u64 {
        self.logical_tile_response_bound_ns
    }

    /// Demonstrated indivisible-item response bound in nanoseconds.
    #[must_use]
    pub const fn indivisible_item_response_bound_ns(self) -> u64 {
        self.indivisible_item_response_bound_ns
    }

    /// Demonstrated external heartbeat interval in nanoseconds.
    #[must_use]
    pub const fn external_heartbeat_bound_ns(self) -> u64 {
        self.external_heartbeat_bound_ns
    }

    /// Exact external-child catalog, if this card admits external execution.
    #[must_use]
    pub const fn external_child_catalog_digest(self) -> Option<ContentHash> {
        self.external_child_catalog_digest
    }

    /// Domain-separated identity of every admitted card field and tier limit.
    #[must_use]
    pub fn digest(self) -> ContentHash {
        let mut payload = b"I14_CANCELLATION_CARD_V2\0".to_vec();
        payload.push(self.tier as u8);
        payload.extend_from_slice(self.semantic_work_unit_digest.as_bytes());
        i14_push_u64(&mut payload, self.campaign_wall_budget_ns);
        i14_push_u64(&mut payload, self.logical_memory_ceiling_bytes);
        payload.push(self.total_resource_unit as u8);
        i14_push_u64(&mut payload, self.total_resource_ceiling);
        i14_push_u64(&mut payload, self.maximum_resource_tile_quantum);
        payload.extend_from_slice(self.resource_budget_authority_digest.as_bytes());
        payload.extend_from_slice(self.logical_partition_spec_digest.as_bytes());
        payload.extend_from_slice(self.execution_environment_digest.as_bytes());
        i14_push_u64(&mut payload, self.logical_tile_response_bound_ns);
        i14_push_u64(&mut payload, self.indivisible_item_response_bound_ns);
        i14_push_u64(&mut payload, self.external_heartbeat_bound_ns);
        match self.external_child_catalog_digest {
            None => payload.push(0),
            Some(digest) => {
                payload.push(1);
                payload.extend_from_slice(digest.as_bytes());
            }
        }
        i14_push_u64(
            &mut payload,
            u64::try_from(self.tier.child_and_observer_cap())
                .expect("I14 child/observer cap fits u64"),
        );
        i14_push_u64(&mut payload, self.tier.request_to_observation_ns());
        i14_push_u64(&mut payload, self.tier.watchdog_quantum_ns());
        i14_push_u64(&mut payload, self.tier.trigger_to_drained_ns());
        i14_push_u64(&mut payload, self.tier.drained_to_finalized_ns());
        i14_push_u64(&mut payload, self.tier.campaign_wall_budget_cap_ns());
        i14_push_u64(&mut payload, self.tier.logical_memory_cap_bytes());
        i14_push_u64(
            &mut payload,
            self.total_resource_unit.total_ceiling_cap(self.tier),
        );
        i14_push_u64(
            &mut payload,
            self.total_resource_unit.tile_quantum_cap(self.tier),
        );
        hash_domain(CANCELLATION_CARD_DOMAIN_V2, &payload)
    }
}

/// Admit a raw Core/ordinary-Max/theorem-Max card into the V2 schema path.
pub fn i14_admit_cancellation_card_v2(
    input: I14CancellationCardInputV2,
) -> Result<I14CancellationCardV2, I14CancellationCardRefusalV2> {
    if input.campaign_wall_budget_ns == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroCampaignWallBudget);
    }
    if input.logical_memory_ceiling_bytes == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroLogicalMemoryCeiling);
    }
    let memory_cap_bytes = input.tier.logical_memory_cap_bytes();
    if input.logical_memory_ceiling_bytes > memory_cap_bytes {
        return Err(
            I14CancellationCardRefusalV2::LogicalMemoryCeilingExceedsTier {
                declared_bytes: input.logical_memory_ceiling_bytes,
                cap_bytes: memory_cap_bytes,
            },
        );
    }
    if input.total_resource_ceiling == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroTotalResourceCeiling);
    }
    if input.maximum_resource_tile_quantum == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroResourceTileQuantum);
    }
    if input.maximum_resource_tile_quantum > input.total_resource_ceiling {
        return Err(
            I14CancellationCardRefusalV2::ResourceTileQuantumExceedsCeiling {
                declared: input.maximum_resource_tile_quantum,
                ceiling: input.total_resource_ceiling,
            },
        );
    }
    let total_resource_cap = input.total_resource_unit.total_ceiling_cap(input.tier);
    if input.total_resource_ceiling > total_resource_cap {
        return Err(
            I14CancellationCardRefusalV2::TotalResourceCeilingExceedsTier {
                declared: input.total_resource_ceiling,
                cap: total_resource_cap,
            },
        );
    }
    let tile_quantum_cap = input.total_resource_unit.tile_quantum_cap(input.tier);
    if input.maximum_resource_tile_quantum > tile_quantum_cap {
        return Err(
            I14CancellationCardRefusalV2::ResourceTileQuantumExceedsTier {
                declared: input.maximum_resource_tile_quantum,
                cap: tile_quantum_cap,
            },
        );
    }
    if input.logical_tile_response_bound_ns == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroLogicalTileResponseBound);
    }
    if input.indivisible_item_response_bound_ns == 0 {
        return Err(I14CancellationCardRefusalV2::ZeroIndivisibleItemResponseBound);
    }
    if (input.external_child_catalog_digest.is_none() && input.external_heartbeat_bound_ns != 0)
        || (input.external_child_catalog_digest.is_some() && input.external_heartbeat_bound_ns == 0)
    {
        return Err(I14CancellationCardRefusalV2::ExternalChildPolicyMismatch);
    }
    let campaign_cap = input.tier.campaign_wall_budget_cap_ns();
    if input.campaign_wall_budget_ns > campaign_cap {
        return Err(
            I14CancellationCardRefusalV2::CampaignWallBudgetExceedsTier {
                declared_ns: input.campaign_wall_budget_ns,
                cap_ns: campaign_cap,
            },
        );
    }
    let observation_cap = input.tier.request_to_observation_ns();
    if input.logical_tile_response_bound_ns > observation_cap {
        return Err(
            I14CancellationCardRefusalV2::LogicalTileResponseExceedsTier {
                declared_ns: input.logical_tile_response_bound_ns,
                cap_ns: observation_cap,
            },
        );
    }
    let watchdog_cap = input.tier.watchdog_quantum_ns();
    if input.indivisible_item_response_bound_ns > watchdog_cap {
        return Err(
            I14CancellationCardRefusalV2::IndivisibleItemResponseExceedsTier {
                declared_ns: input.indivisible_item_response_bound_ns,
                cap_ns: watchdog_cap,
            },
        );
    }
    if input.external_child_catalog_digest.is_some()
        && input.external_heartbeat_bound_ns > watchdog_cap
    {
        return Err(I14CancellationCardRefusalV2::ExternalHeartbeatExceedsTier {
            declared_ns: input.external_heartbeat_bound_ns,
            cap_ns: watchdog_cap,
        });
    }
    Ok(I14CancellationCardV2 {
        tier: input.tier,
        semantic_work_unit_digest: input.semantic_work_unit_digest,
        campaign_wall_budget_ns: input.campaign_wall_budget_ns,
        logical_memory_ceiling_bytes: input.logical_memory_ceiling_bytes,
        total_resource_unit: input.total_resource_unit,
        total_resource_ceiling: input.total_resource_ceiling,
        maximum_resource_tile_quantum: input.maximum_resource_tile_quantum,
        resource_budget_authority_digest: input.resource_budget_authority_digest,
        logical_partition_spec_digest: input.logical_partition_spec_digest,
        execution_environment_digest: input.execution_environment_digest,
        logical_tile_response_bound_ns: input.logical_tile_response_bound_ns,
        indivisible_item_response_bound_ns: input.indivisible_item_response_bound_ns,
        external_heartbeat_bound_ns: input.external_heartbeat_bound_ns,
        external_child_catalog_digest: input.external_child_catalog_digest,
    })
}

/// One boundary record in a schema-authoritative genesis-to-terminal trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14TerminalBoundaryRecordV2<'a> {
    /// Existing cause-arbitration boundary.
    pub boundary: I14TerminalBoundaryV1<'a>,
    /// Child scopes still in flight at this boundary.
    pub in_flight_children: u16,
    /// Most recent independent watchdog poll visible at this boundary.
    pub last_watchdog_poll_monotonic_ns: u64,
    /// Cumulative consumption on the card's total-resource axis.
    pub total_resource_consumed: u64,
    /// Frozen work items not yet completed at this boundary.
    pub work_items_remaining: u64,
    /// Cost of the receipt-bound next work item rejected by the governor, or
    /// `None` when no next item was rejected at this boundary.
    pub rejected_next_work_resource_cost: Option<u64>,
    /// Identity of the resource-ledger prefix through this boundary.
    pub resource_ledger_prefix_digest: ContentHash,
    /// Identity of the scheduled-work frontier through this boundary.
    pub work_frontier_prefix_digest: ContentHash,
}

/// Complete bounded trace supplied to authoritative first-terminal selection.
#[derive(Clone, Copy, Debug)]
pub struct I14TerminalBoundaryTraceV2<'a> {
    /// Independent signature/adjudication verification-receipt identity asserted
    /// for the logical stream. This crate binds but does not authenticate it.
    pub logical_execution_verification_receipt_digest: ContentHash,
    /// Validated cancellation contract for this exact semantic work unit.
    pub cancellation_card: I14CancellationCardV2,
    /// Calibrated campaign-start timestamp.
    pub campaign_started_monotonic_ns: u64,
    /// Genesis-to-frontier or genesis-to-first-terminal boundary records.
    pub boundaries: &'a [I14TerminalBoundaryRecordV2<'a>],
    /// Complete request/observation trace no later than the frontier boundary.
    pub requests: &'a [I14CancellationRequestV1],
}

/// Clock-free projection of one validated boundary prefix record.
#[derive(Clone, Debug, PartialEq, Eq)]
struct I14CanonicalBoundaryRecordV2 {
    boundary_ordinal: u64,
    logical_sequence: u64,
    scope_ancestry: Vec<u64>,
    admitted_observer_tile_ids: Vec<u64>,
    in_flight_children: u16,
    total_resource_consumed: u64,
    work_items_remaining: u64,
    rejected_next_work_resource_cost: Option<u64>,
    resource_ledger_prefix_digest: ContentHash,
    work_frontier_prefix_digest: ContentHash,
    cause_candidates: I14TerminalCauseCandidatesV1,
    decision: I14TerminalBoundaryDecisionV1,
}

/// Clock-free semantic projection of one cancellation request at the supplied
/// terminal/frontier cut. Raw calibrated times remain telemetry-only, while
/// every timing-derived state and logical observation is prefix-bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct I14CanonicalTraceRequestV2 {
    request_id: u64,
    scope_root: u64,
    logical_sequence: u64,
    state: I14CancellationRequestStateV1,
    observation: Option<I14CanonicalCancellationObservationV1>,
}

/// Opaque proof that the final supplied boundary is the first terminal one.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct I14FirstTerminalSelectionV2 {
    logical_execution_receipt_digest: ContentHash,
    logical_execution_verification_receipt_digest: ContentHash,
    cancellation_card: I14CancellationCardV2,
    selected_index: usize,
    boundary_prefix: Vec<I14CanonicalBoundaryRecordV2>,
    decision: I14TerminalBoundaryDecisionV1,
    prefix_digest: ContentHash,
}

impl I14FirstTerminalSelectionV2 {
    /// Recomputed clock-free content identity of the complete trace.
    #[must_use]
    pub const fn logical_execution_receipt_digest(&self) -> ContentHash {
        self.logical_execution_receipt_digest
    }

    /// Independent verification-receipt identity bound to the logical stream.
    ///
    /// Signature and authority validation occur in the HELM/ledger gate; this
    /// lower-layer value alone is deliberately not an authentication theorem.
    #[must_use]
    pub const fn logical_execution_verification_receipt_digest(&self) -> ContentHash {
        self.logical_execution_verification_receipt_digest
    }

    /// Cancellation contract validated for this trace.
    #[must_use]
    pub const fn cancellation_card(&self) -> I14CancellationCardV2 {
        self.cancellation_card
    }

    /// Zero-based position of the first selected terminal boundary.
    #[must_use]
    pub const fn selected_index(&self) -> usize {
        self.selected_index
    }

    /// Number of validated boundaries from genesis through selection.
    #[must_use]
    pub fn boundary_count(&self) -> usize {
        self.boundary_prefix.len()
    }

    /// Selected terminal decision.
    #[must_use]
    pub const fn decision(&self) -> I14TerminalBoundaryDecisionV1 {
        self.decision
    }

    /// Digest of the complete clock-free boundary prefix.
    #[must_use]
    pub const fn prefix_digest(&self) -> ContentHash {
        self.prefix_digest
    }
}

/// Opaque proof that a complete bounded prefix contains no terminal boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct I14TerminalFrontierCertificateV2 {
    logical_execution_receipt_digest: ContentHash,
    logical_execution_verification_receipt_digest: ContentHash,
    cancellation_card: I14CancellationCardV2,
    boundary_count: usize,
    last_boundary_ordinal: u64,
    last_boundary_logical_sequence: u64,
    prefix_digest: ContentHash,
}

impl I14TerminalFrontierCertificateV2 {
    /// Recomputed clock-free content identity of the complete trace.
    #[must_use]
    pub const fn logical_execution_receipt_digest(&self) -> ContentHash {
        self.logical_execution_receipt_digest
    }

    /// Independent verification-receipt identity bound to the logical stream.
    #[must_use]
    pub const fn logical_execution_verification_receipt_digest(&self) -> ContentHash {
        self.logical_execution_verification_receipt_digest
    }

    /// Cancellation contract validated for this frontier.
    #[must_use]
    pub const fn cancellation_card(&self) -> I14CancellationCardV2 {
        self.cancellation_card
    }

    /// Number of validated nonterminal boundaries.
    #[must_use]
    pub const fn boundary_count(&self) -> usize {
        self.boundary_count
    }

    /// Last validated frontier ordinal.
    #[must_use]
    pub const fn last_boundary_ordinal(&self) -> u64 {
        self.last_boundary_ordinal
    }

    /// Coordinator sequence of the last validated frontier boundary.
    #[must_use]
    pub const fn last_boundary_logical_sequence(&self) -> u64 {
        self.last_boundary_logical_sequence
    }

    /// Digest committing to every validated nonterminal boundary and decision.
    #[must_use]
    pub const fn prefix_digest(&self) -> ContentHash {
        self.prefix_digest
    }
}

/// Authoritative terminal-trace result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum I14TerminalTraceOutcomeV2 {
    /// The final supplied record is exactly the first terminal boundary.
    Selected(I14FirstTerminalSelectionV2),
    /// The complete supplied prefix is valid but remains nonterminal.
    Frontier(I14TerminalFrontierCertificateV2),
}

/// Fail-closed refusal of an authoritative first-terminal trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14TerminalTraceRefusalV2 {
    /// The trace contains no boundary, so neither genesis nor a frontier exists.
    EmptyBoundaryTrace,
    /// The trace exceeds its admitted bound.
    TooManyBoundaries {
        /// Supplied record count.
        count: usize,
        /// Maximum admitted record count.
        cap: usize,
    },
    /// A boundary's root-to-leaf scope ancestry exceeds the versioned schema
    /// cap. This is refused before cloning or comparing the ancestry.
    ScopeAncestryTooDeep {
        /// Boundary carrying the excess ancestry.
        boundary_ordinal: u64,
        /// Supplied ancestry depth.
        depth: usize,
        /// Versioned ancestry cap.
        cap: usize,
    },
    /// The boundary/request Cartesian work estimate exceeds the versioned
    /// arbitration budget of the reference selector.
    ArbitrationWorkBudgetExceeded {
        /// Supplied boundary count.
        boundary_count: usize,
        /// Supplied request count.
        request_count: usize,
        /// Exact multiplication of the two counts.
        pair_count: usize,
        /// Maximum admitted pair count.
        cap: usize,
    },
    /// The first boundary does not carry the required genesis ordinal.
    FirstBoundaryNotGenesis {
        /// Supplied first ordinal.
        found: u64,
    },
    /// A boundary ordinal is not the exact successor of its predecessor.
    NonContiguousBoundaryOrdinal {
        /// Position of the malformed record.
        index: usize,
        /// Required ordinal.
        expected: u64,
        /// Supplied ordinal.
        found: u64,
    },
    /// The predecessor ordinal cannot be incremented representably.
    BoundaryOrdinalOverflow {
        /// Predecessor ordinal.
        previous: u64,
    },
    /// Boundary coordinator sequences are not strictly increasing.
    NonIncreasingBoundaryLogicalSequence {
        /// Position of the malformed record.
        index: usize,
        /// Previous sequence.
        previous: u64,
        /// Supplied sequence.
        found: u64,
    },
    /// Calibrated boundary timestamps move backward.
    NonMonotonicBoundaryTimestamp {
        /// Position of the malformed record.
        index: usize,
        /// Previous or campaign-start timestamp.
        previous_ns: u64,
        /// Supplied timestamp.
        found_ns: u64,
    },
    /// Scope ancestry changed within one semantic-work-unit trace.
    ScopeAncestryChanged {
        /// Boundary where the frozen path changed.
        boundary_ordinal: u64,
    },
    /// Observer membership changed within one semantic-work-unit trace.
    ObserverCatalogChanged {
        /// Boundary where the frozen catalog changed.
        boundary_ordinal: u64,
    },
    /// A boundary exceeds its tier's in-flight-child cap.
    TooManyInFlightChildren {
        /// Boundary carrying the excess.
        boundary_ordinal: u64,
        /// Supplied child count.
        count: usize,
        /// Tier cap.
        cap: usize,
    },
    /// A boundary exceeds its tier's observer-catalog cap.
    TooManyObserverTilesForTier {
        /// Boundary carrying the excess.
        boundary_ordinal: u64,
        /// Supplied observer count.
        count: usize,
        /// Tier cap.
        cap: usize,
    },
    /// The claimed latest watchdog poll occurs after its boundary.
    WatchdogPollAfterBoundary {
        /// Malformed boundary.
        boundary_ordinal: u64,
    },
    /// A boundary claims a watchdog poll before campaign admission.
    WatchdogPollBeforeCampaign {
        /// Malformed boundary.
        boundary_ordinal: u64,
    },
    /// Latest-watchdog-poll time moves backward across boundary order.
    NonMonotonicWatchdogPoll {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Prior boundary's latest poll.
        previous_ns: u64,
        /// Supplied latest poll.
        found_ns: u64,
    },
    /// A watchdog lapse was not represented as an infrastructure failure.
    WatchdogFailureNotReflected {
        /// Boundary carrying the lapse.
        boundary_ordinal: u64,
        /// Measured interval since the latest poll.
        elapsed_ns: u64,
        /// Tier maximum.
        cap_ns: u64,
    },
    /// Campaign start plus its admitted wall budget overflowed `u64`.
    CampaignDeadlineOverflow,
    /// Wall-budget expiry was not represented as a timeout candidate.
    CampaignTimeoutNotReflected {
        /// Boundary at or beyond expiry.
        boundary_ordinal: u64,
        /// Derived immutable campaign deadline.
        deadline_ns: u64,
        /// Boundary timestamp.
        boundary_ns: u64,
    },
    /// A caller-provided request deadline differs from the card's exact SLO.
    ObservationDeadlineSloMismatch {
        /// Request carrying the malformed deadline.
        request_id: u64,
        /// Supplied deadline delta.
        declared_delta_ns: u64,
        /// Required deadline delta.
        required_delta_ns: u64,
    },
    /// A cancellation request predates campaign admission.
    RequestBeforeCampaign {
        /// Malformed request identity.
        request_id: u64,
        /// Supplied request timestamp.
        requested_ns: u64,
        /// Frozen campaign start.
        campaign_started_ns: u64,
    },
    /// A cancellation observation predates campaign admission.
    ObservationBeforeCampaign {
        /// Owning request identity.
        request_id: u64,
        /// Supplied observation timestamp.
        observed_ns: u64,
        /// Frozen campaign start.
        campaign_started_ns: u64,
    },
    /// A request or observation lies after the supplied complete frontier.
    EventAfterTraceFrontier {
        /// Owning request identity.
        request_id: u64,
        /// Offending logical sequence.
        logical_sequence: u64,
        /// Final supplied boundary sequence.
        frontier_logical_sequence: u64,
    },
    /// Cumulative resource consumption moved backward.
    NonMonotonicResourceConsumption {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Prior cumulative consumption.
        previous: u64,
        /// Supplied cumulative consumption.
        found: u64,
    },
    /// Remaining frozen work increased at a later boundary.
    NonMonotonicWorkRemaining {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Prior remaining count.
        previous: u64,
        /// Supplied remaining count.
        found: u64,
    },
    /// Cumulative resource consumption exceeds the hard ceiling.
    ResourceConsumptionExceedsCeiling {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Supplied cumulative consumption.
        consumed: u64,
        /// Frozen ceiling.
        ceiling: u64,
    },
    /// A rejected next-work witness is zero or exceeds the admitted quantum.
    RejectedWorkQuantumInvalid {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Rejected next-item cost.
        rejected_cost: u64,
        /// Largest admitted indivisible quantum.
        tile_quantum: u64,
    },
    /// The governor claims to reject work that still fits below the ceiling.
    RejectedWorkFitsWithinCeiling {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Cumulative committed consumption.
        consumed: u64,
        /// Rejected next-item cost.
        rejected_cost: u64,
        /// Frozen ceiling.
        ceiling: u64,
    },
    /// A completed work frontier also carries a rejected-next-work witness.
    RejectedWorkAfterCompletion {
        /// Malformed boundary.
        boundary_ordinal: u64,
    },
    /// The `budget_exhausted` candidate disagrees with the resource ledger.
    BudgetExhaustionMismatch {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Derivation from cumulative consumption and ceiling.
        expected: bool,
        /// Caller-supplied cause bit.
        found: bool,
    },
    /// The `completed` candidate disagrees with the frozen work frontier.
    CompletionMismatch {
        /// Malformed boundary.
        boundary_ordinal: u64,
        /// Derivation from remaining work items.
        expected: bool,
        /// Caller-supplied cause bit.
        found: bool,
    },
    /// A selected terminal boundary still reports live child scopes.
    TerminalBoundaryHasInFlightChildren {
        /// Malformed terminal boundary.
        boundary_ordinal: u64,
        /// Live child count.
        count: usize,
    },
    /// One local boundary/request cut is malformed.
    BoundaryCause {
        /// Boundary at which validation failed.
        boundary_ordinal: u64,
        /// Local fail-closed diagnosis.
        source: I14TerminalCauseRefusalV1,
    },
    /// The caller supplied records after the first selected terminal boundary.
    BoundaryAfterFirstTerminal {
        /// First selected terminal ordinal.
        first_terminal_ordinal: u64,
        /// First forbidden trailing ordinal.
        trailing_ordinal: u64,
    },
}

fn i14_push_canonical_boundary_v2(payload: &mut Vec<u8>, record: &I14CanonicalBoundaryRecordV2) {
    i14_push_u64(payload, record.boundary_ordinal);
    i14_push_u64(payload, record.logical_sequence);
    i14_push_u64_slice(payload, &record.scope_ancestry);
    i14_push_u64_slice(payload, &record.admitted_observer_tile_ids);
    i14_push_u64(payload, u64::from(record.in_flight_children));
    i14_push_u64(payload, record.total_resource_consumed);
    i14_push_u64(payload, record.work_items_remaining);
    i14_push_optional_u64(payload, record.rejected_next_work_resource_cost);
    payload.extend_from_slice(record.resource_ledger_prefix_digest.as_bytes());
    payload.extend_from_slice(record.work_frontier_prefix_digest.as_bytes());
    payload.extend_from_slice(&[
        u8::from(record.cause_candidates.infrastructure_failed),
        u8::from(record.cause_candidates.timed_out),
        u8::from(record.cause_candidates.budget_exhausted),
        u8::from(record.cause_candidates.completed),
    ]);
    i14_push_terminal_decision(payload, record.decision);
}

fn i14_push_canonical_trace_request_v2(payload: &mut Vec<u8>, request: I14CanonicalTraceRequestV2) {
    i14_push_u64(payload, request.request_id);
    i14_push_u64(payload, request.scope_root);
    i14_push_u64(payload, request.logical_sequence);
    payload.push(request.state as u8);
    match request.observation {
        None => payload.push(0),
        Some(observation) => {
            payload.push(1);
            i14_push_u64(payload, observation.logical_sequence);
            i14_push_u64(payload, observation.observing_tile_id);
            i14_push_optional_u64(payload, observation.latest_completed_boundary_ordinal);
        }
    }
}

fn i14_logical_execution_trace_digest_v2(
    cancellation_card: I14CancellationCardV2,
    records: &[I14CanonicalBoundaryRecordV2],
    requests: &[I14CanonicalTraceRequestV2],
) -> ContentHash {
    let mut payload = b"I14_LOGICAL_EXECUTION_TRACE_V2\0".to_vec();
    payload.extend_from_slice(cancellation_card.digest().as_bytes());
    i14_push_len(&mut payload, records.len());
    for record in records {
        i14_push_canonical_boundary_v2(&mut payload, record);
    }
    i14_push_len(&mut payload, requests.len());
    for request in requests {
        i14_push_canonical_trace_request_v2(&mut payload, *request);
    }
    hash_domain(LOGICAL_EXECUTION_TRACE_DOMAIN_V2, &payload)
}

fn i14_canonical_trace_requests_v2(
    frontier: I14TerminalBoundaryV1<'_>,
    requests: &[I14CancellationRequestV1],
) -> Vec<I14CanonicalTraceRequestV2> {
    let scope_ids = frontier
        .scope_ancestry
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let mut requests = requests
        .iter()
        .copied()
        .map(|request| I14CanonicalTraceRequestV2 {
            request_id: request.request_id,
            scope_root: request.scope_root,
            logical_sequence: request.logical_sequence,
            state: i14_request_state_v1(frontier, &scope_ids, request),
            observation: request.observation.map(|observation| {
                I14CanonicalCancellationObservationV1 {
                    logical_sequence: observation.logical_sequence,
                    observing_tile_id: observation.observing_tile_id,
                    latest_completed_boundary_ordinal: observation
                        .latest_completed_boundary_ordinal,
                }
            }),
        })
        .collect::<Vec<_>>();
    requests.sort_unstable_by_key(|request| (request.logical_sequence, request.request_id));
    requests
}

fn i14_terminal_prefix_digest_v2(
    logical_execution_receipt_digest: ContentHash,
    logical_execution_verification_receipt_digest: ContentHash,
    cancellation_card: I14CancellationCardV2,
    records: &[I14CanonicalBoundaryRecordV2],
    requests: &[I14CanonicalTraceRequestV2],
) -> ContentHash {
    let mut payload = b"I14_TERMINAL_BOUNDARY_PREFIX_V2\0".to_vec();
    payload.extend_from_slice(logical_execution_receipt_digest.as_bytes());
    payload.extend_from_slice(logical_execution_verification_receipt_digest.as_bytes());
    payload.extend_from_slice(cancellation_card.digest().as_bytes());
    i14_push_len(&mut payload, records.len());
    for record in records {
        i14_push_canonical_boundary_v2(&mut payload, record);
    }
    i14_push_len(&mut payload, requests.len());
    for request in requests {
        i14_push_canonical_trace_request_v2(&mut payload, *request);
    }
    hash_domain(TERMINAL_TRACE_DOMAIN_V2, &payload)
}

fn i14_sorted_observers(boundary: I14TerminalBoundaryV1<'_>) -> Vec<u64> {
    let mut observers = boundary.admitted_observer_tile_ids.to_vec();
    observers.sort_unstable();
    observers
}

/// Validate a complete boundary prefix and select only its first terminal.
///
/// Unlike [`i14_select_terminal_boundary_v1`], this is an authority-bearing
/// trace operation. It requires genesis, contiguous immutable boundary order,
/// an exact tier cancellation card, a fixed scope/observer catalog, bounded
/// watchdog/child state, and no record after the first selected terminal. The
/// returned proof is schema-authoritative only: a consuming HELM/ledger gate
/// must authenticate the verification receipt identified by the bound digest
/// before promotion.
#[allow(clippy::too_many_lines)]
pub fn i14_select_first_terminal_boundary_v2(
    trace: I14TerminalBoundaryTraceV2<'_>,
) -> Result<I14TerminalTraceOutcomeV2, I14TerminalTraceRefusalV2> {
    if trace.boundaries.is_empty() {
        return Err(I14TerminalTraceRefusalV2::EmptyBoundaryTrace);
    }
    if trace.boundaries.len() > I14_MAX_TERMINAL_BOUNDARIES_V2 {
        return Err(I14TerminalTraceRefusalV2::TooManyBoundaries {
            count: trace.boundaries.len(),
            cap: I14_MAX_TERMINAL_BOUNDARIES_V2,
        });
    }
    if trace.requests.len() > I14_MAX_CANCELLATION_REQUESTS_V1 {
        return Err(I14TerminalTraceRefusalV2::BoundaryCause {
            boundary_ordinal: trace.boundaries[0].boundary.boundary_ordinal,
            source: I14TerminalCauseRefusalV1::TooManyRequests {
                count: trace.requests.len(),
                cap: I14_MAX_CANCELLATION_REQUESTS_V1,
            },
        });
    }
    let first = trace.boundaries[0].boundary;
    if first.boundary_ordinal != I14_TERMINAL_BOUNDARY_GENESIS_ORDINAL_V2 {
        return Err(I14TerminalTraceRefusalV2::FirstBoundaryNotGenesis {
            found: first.boundary_ordinal,
        });
    }
    let campaign_deadline_ns = trace
        .campaign_started_monotonic_ns
        .checked_add(trace.cancellation_card.campaign_wall_budget_ns())
        .ok_or(I14TerminalTraceRefusalV2::CampaignDeadlineOverflow)?;
    let tier = trace.cancellation_card.tier();
    let child_cap = tier.child_and_observer_cap();
    let watchdog_cap_ns = tier.watchdog_quantum_ns();
    // Reject attacker-sized variable-width fields across the complete trace
    // before cloning/sorting the first catalog or comparing any full ancestry.
    // The later semantic pass may therefore assume these allocations and
    // comparisons are bounded by the frozen schema/tier contracts.
    for record in trace.boundaries {
        let boundary = record.boundary;
        if boundary.scope_ancestry.len() > I14_MAX_SCOPE_ANCESTRY_V1 {
            return Err(I14TerminalTraceRefusalV2::ScopeAncestryTooDeep {
                boundary_ordinal: boundary.boundary_ordinal,
                depth: boundary.scope_ancestry.len(),
                cap: I14_MAX_SCOPE_ANCESTRY_V1,
            });
        }
        if boundary.admitted_observer_tile_ids.len() > child_cap {
            return Err(I14TerminalTraceRefusalV2::TooManyObserverTilesForTier {
                boundary_ordinal: boundary.boundary_ordinal,
                count: boundary.admitted_observer_tile_ids.len(),
                cap: child_cap,
            });
        }
    }
    let expected_scope = first.scope_ancestry;
    let expected_observers = i14_sorted_observers(first);
    let mut previous_ordinal: Option<u64> = None;
    let mut previous_logical_sequence: Option<u64> = None;
    let mut previous_monotonic_ns = trace.campaign_started_monotonic_ns;
    let mut previous_watchdog_poll_ns = trace.campaign_started_monotonic_ns;
    let mut previous_resource_consumed: Option<u64> = None;
    let mut previous_work_items_remaining: Option<u64> = None;

    for (index, record) in trace.boundaries.iter().enumerate() {
        let boundary = record.boundary;
        if let Some(previous) = previous_ordinal {
            let expected = previous
                .checked_add(1)
                .ok_or(I14TerminalTraceRefusalV2::BoundaryOrdinalOverflow { previous })?;
            if boundary.boundary_ordinal != expected {
                return Err(I14TerminalTraceRefusalV2::NonContiguousBoundaryOrdinal {
                    index,
                    expected,
                    found: boundary.boundary_ordinal,
                });
            }
        }
        if let Some(previous) = previous_logical_sequence
            && boundary.logical_sequence <= previous
        {
            return Err(
                I14TerminalTraceRefusalV2::NonIncreasingBoundaryLogicalSequence {
                    index,
                    previous,
                    found: boundary.logical_sequence,
                },
            );
        }
        if boundary.monotonic_ns < previous_monotonic_ns {
            return Err(I14TerminalTraceRefusalV2::NonMonotonicBoundaryTimestamp {
                index,
                previous_ns: previous_monotonic_ns,
                found_ns: boundary.monotonic_ns,
            });
        }
        if boundary.scope_ancestry != expected_scope {
            return Err(I14TerminalTraceRefusalV2::ScopeAncestryChanged {
                boundary_ordinal: boundary.boundary_ordinal,
            });
        }
        if i14_sorted_observers(boundary) != expected_observers {
            return Err(I14TerminalTraceRefusalV2::ObserverCatalogChanged {
                boundary_ordinal: boundary.boundary_ordinal,
            });
        }
        if usize::from(record.in_flight_children) > child_cap {
            return Err(I14TerminalTraceRefusalV2::TooManyInFlightChildren {
                boundary_ordinal: boundary.boundary_ordinal,
                count: usize::from(record.in_flight_children),
                cap: child_cap,
            });
        }
        if let Some(previous) = previous_resource_consumed
            && record.total_resource_consumed < previous
        {
            return Err(I14TerminalTraceRefusalV2::NonMonotonicResourceConsumption {
                boundary_ordinal: boundary.boundary_ordinal,
                previous,
                found: record.total_resource_consumed,
            });
        }
        if let Some(previous) = previous_work_items_remaining
            && record.work_items_remaining > previous
        {
            return Err(I14TerminalTraceRefusalV2::NonMonotonicWorkRemaining {
                boundary_ordinal: boundary.boundary_ordinal,
                previous,
                found: record.work_items_remaining,
            });
        }
        let resource_ceiling = trace.cancellation_card.total_resource_ceiling();
        let tile_quantum = trace.cancellation_card.maximum_resource_tile_quantum();
        if record.total_resource_consumed > resource_ceiling {
            return Err(
                I14TerminalTraceRefusalV2::ResourceConsumptionExceedsCeiling {
                    boundary_ordinal: boundary.boundary_ordinal,
                    consumed: record.total_resource_consumed,
                    ceiling: resource_ceiling,
                },
            );
        }
        if record.work_items_remaining == 0 && record.rejected_next_work_resource_cost.is_some() {
            return Err(I14TerminalTraceRefusalV2::RejectedWorkAfterCompletion {
                boundary_ordinal: boundary.boundary_ordinal,
            });
        }
        if let Some(rejected_cost) = record.rejected_next_work_resource_cost {
            if rejected_cost == 0 || rejected_cost > tile_quantum {
                return Err(I14TerminalTraceRefusalV2::RejectedWorkQuantumInvalid {
                    boundary_ordinal: boundary.boundary_ordinal,
                    rejected_cost,
                    tile_quantum,
                });
            }
            if record
                .total_resource_consumed
                .checked_add(rejected_cost)
                .is_some_and(|next| next <= resource_ceiling)
            {
                return Err(I14TerminalTraceRefusalV2::RejectedWorkFitsWithinCeiling {
                    boundary_ordinal: boundary.boundary_ordinal,
                    consumed: record.total_resource_consumed,
                    rejected_cost,
                    ceiling: resource_ceiling,
                });
            }
        }
        let expected_budget_exhaustion =
            record.work_items_remaining > 0 && record.rejected_next_work_resource_cost.is_some();
        if boundary.budget_exhausted != expected_budget_exhaustion {
            return Err(I14TerminalTraceRefusalV2::BudgetExhaustionMismatch {
                boundary_ordinal: boundary.boundary_ordinal,
                expected: expected_budget_exhaustion,
                found: boundary.budget_exhausted,
            });
        }
        let expected_completion = record.work_items_remaining == 0;
        if boundary.completed != expected_completion {
            return Err(I14TerminalTraceRefusalV2::CompletionMismatch {
                boundary_ordinal: boundary.boundary_ordinal,
                expected: expected_completion,
                found: boundary.completed,
            });
        }
        if record.last_watchdog_poll_monotonic_ns > boundary.monotonic_ns {
            return Err(I14TerminalTraceRefusalV2::WatchdogPollAfterBoundary {
                boundary_ordinal: boundary.boundary_ordinal,
            });
        }
        if record.last_watchdog_poll_monotonic_ns < trace.campaign_started_monotonic_ns {
            return Err(I14TerminalTraceRefusalV2::WatchdogPollBeforeCampaign {
                boundary_ordinal: boundary.boundary_ordinal,
            });
        }
        if record.last_watchdog_poll_monotonic_ns < previous_watchdog_poll_ns {
            return Err(I14TerminalTraceRefusalV2::NonMonotonicWatchdogPoll {
                boundary_ordinal: boundary.boundary_ordinal,
                previous_ns: previous_watchdog_poll_ns,
                found_ns: record.last_watchdog_poll_monotonic_ns,
            });
        }
        let watchdog_elapsed_ns = boundary.monotonic_ns - record.last_watchdog_poll_monotonic_ns;
        if watchdog_elapsed_ns > watchdog_cap_ns && !boundary.infrastructure_failed {
            return Err(I14TerminalTraceRefusalV2::WatchdogFailureNotReflected {
                boundary_ordinal: boundary.boundary_ordinal,
                elapsed_ns: watchdog_elapsed_ns,
                cap_ns: watchdog_cap_ns,
            });
        }
        if boundary.monotonic_ns >= campaign_deadline_ns && !boundary.timed_out {
            return Err(I14TerminalTraceRefusalV2::CampaignTimeoutNotReflected {
                boundary_ordinal: boundary.boundary_ordinal,
                deadline_ns: campaign_deadline_ns,
                boundary_ns: boundary.monotonic_ns,
            });
        }
        previous_ordinal = Some(boundary.boundary_ordinal);
        previous_logical_sequence = Some(boundary.logical_sequence);
        previous_monotonic_ns = boundary.monotonic_ns;
        previous_watchdog_poll_ns = record.last_watchdog_poll_monotonic_ns;
        previous_resource_consumed = Some(record.total_resource_consumed);
        previous_work_items_remaining = Some(record.work_items_remaining);
    }

    let frontier_boundary = trace
        .boundaries
        .last()
        .expect("nonempty trace checked above")
        .boundary;
    let frontier_sequence = frontier_boundary.logical_sequence;
    let mut ordered_requests = trace.requests.iter().collect::<Vec<_>>();
    ordered_requests.sort_unstable_by_key(|request| {
        (
            request.logical_sequence,
            request.request_id,
            request.scope_root,
            request.requested_monotonic_ns,
            request.observation_deadline_ns,
            request.observation,
        )
    });
    for request in ordered_requests {
        if request.requested_monotonic_ns < trace.campaign_started_monotonic_ns {
            return Err(I14TerminalTraceRefusalV2::RequestBeforeCampaign {
                request_id: request.request_id,
                requested_ns: request.requested_monotonic_ns,
                campaign_started_ns: trace.campaign_started_monotonic_ns,
            });
        }
        if let Some(observation) = request.observation
            && observation.monotonic_ns < trace.campaign_started_monotonic_ns
        {
            return Err(I14TerminalTraceRefusalV2::ObservationBeforeCampaign {
                request_id: request.request_id,
                observed_ns: observation.monotonic_ns,
                campaign_started_ns: trace.campaign_started_monotonic_ns,
            });
        }
        if request.observation_deadline_ns < request.requested_monotonic_ns {
            return Err(I14TerminalTraceRefusalV2::BoundaryCause {
                boundary_ordinal: first.boundary_ordinal,
                source: I14TerminalCauseRefusalV1::DeadlineBeforeRequest {
                    request_id: request.request_id,
                },
            });
        }
        let declared_delta_ns = request.observation_deadline_ns - request.requested_monotonic_ns;
        let required_delta_ns = tier.request_to_observation_ns();
        if declared_delta_ns != required_delta_ns {
            return Err(I14TerminalTraceRefusalV2::ObservationDeadlineSloMismatch {
                request_id: request.request_id,
                declared_delta_ns,
                required_delta_ns,
            });
        }
        if request.logical_sequence > frontier_sequence {
            return Err(I14TerminalTraceRefusalV2::EventAfterTraceFrontier {
                request_id: request.request_id,
                logical_sequence: request.logical_sequence,
                frontier_logical_sequence: frontier_sequence,
            });
        }
        if request.logical_sequence == frontier_sequence {
            return Err(I14TerminalTraceRefusalV2::BoundaryCause {
                boundary_ordinal: frontier_boundary.boundary_ordinal,
                source: I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
                    logical_sequence: frontier_sequence,
                },
            });
        }
        if let Some(observation) = request.observation
            && observation.logical_sequence > frontier_sequence
        {
            return Err(I14TerminalTraceRefusalV2::EventAfterTraceFrontier {
                request_id: request.request_id,
                logical_sequence: observation.logical_sequence,
                frontier_logical_sequence: frontier_sequence,
            });
        }
        if let Some(observation) = request.observation
            && observation.logical_sequence == frontier_sequence
        {
            return Err(I14TerminalTraceRefusalV2::BoundaryCause {
                boundary_ordinal: frontier_boundary.boundary_ordinal,
                source: I14TerminalCauseRefusalV1::DuplicateLogicalSequence {
                    logical_sequence: frontier_sequence,
                },
            });
        }
    }
    let canonical_requests = i14_canonical_trace_requests_v2(frontier_boundary, trace.requests);
    let arbitration_pair_count = trace.boundaries.len().saturating_mul(trace.requests.len());
    if arbitration_pair_count > I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2 {
        return Err(I14TerminalTraceRefusalV2::ArbitrationWorkBudgetExceeded {
            boundary_count: trace.boundaries.len(),
            request_count: trace.requests.len(),
            pair_count: arbitration_pair_count,
            cap: I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2,
        });
    }

    let mut canonical_prefix = Vec::with_capacity(trace.boundaries.len());
    for (index, record) in trace.boundaries.iter().enumerate() {
        let boundary = record.boundary;
        let decision =
            i14_select_terminal_boundary_v1(boundary, trace.requests).map_err(|source| {
                I14TerminalTraceRefusalV2::BoundaryCause {
                    boundary_ordinal: boundary.boundary_ordinal,
                    source,
                }
            })?;
        canonical_prefix.push(I14CanonicalBoundaryRecordV2 {
            boundary_ordinal: boundary.boundary_ordinal,
            logical_sequence: boundary.logical_sequence,
            scope_ancestry: boundary.scope_ancestry.to_vec(),
            admitted_observer_tile_ids: i14_sorted_observers(boundary),
            in_flight_children: record.in_flight_children,
            total_resource_consumed: record.total_resource_consumed,
            work_items_remaining: record.work_items_remaining,
            rejected_next_work_resource_cost: record.rejected_next_work_resource_cost,
            resource_ledger_prefix_digest: record.resource_ledger_prefix_digest,
            work_frontier_prefix_digest: record.work_frontier_prefix_digest,
            cause_candidates: I14TerminalCauseCandidatesV1 {
                infrastructure_failed: boundary.infrastructure_failed,
                timed_out: boundary.timed_out,
                budget_exhausted: boundary.budget_exhausted,
                completed: boundary.completed,
            },
            decision,
        });
        if matches!(decision, I14TerminalBoundaryDecisionV1::Selected { .. }) {
            if record.in_flight_children != 0 {
                return Err(
                    I14TerminalTraceRefusalV2::TerminalBoundaryHasInFlightChildren {
                        boundary_ordinal: boundary.boundary_ordinal,
                        count: usize::from(record.in_flight_children),
                    },
                );
            }
            if index + 1 != trace.boundaries.len() {
                return Err(I14TerminalTraceRefusalV2::BoundaryAfterFirstTerminal {
                    first_terminal_ordinal: boundary.boundary_ordinal,
                    trailing_ordinal: trace.boundaries[index + 1].boundary.boundary_ordinal,
                });
            }
            let logical_execution_receipt_digest = i14_logical_execution_trace_digest_v2(
                trace.cancellation_card,
                &canonical_prefix,
                &canonical_requests,
            );
            let prefix_digest = i14_terminal_prefix_digest_v2(
                logical_execution_receipt_digest,
                trace.logical_execution_verification_receipt_digest,
                trace.cancellation_card,
                &canonical_prefix,
                &canonical_requests,
            );
            return Ok(I14TerminalTraceOutcomeV2::Selected(
                I14FirstTerminalSelectionV2 {
                    logical_execution_receipt_digest,
                    logical_execution_verification_receipt_digest: trace
                        .logical_execution_verification_receipt_digest,
                    cancellation_card: trace.cancellation_card,
                    selected_index: index,
                    boundary_prefix: canonical_prefix,
                    decision,
                    prefix_digest,
                },
            ));
        }
    }

    let last = trace
        .boundaries
        .last()
        .expect("nonempty trace checked above");
    let logical_execution_receipt_digest = i14_logical_execution_trace_digest_v2(
        trace.cancellation_card,
        &canonical_prefix,
        &canonical_requests,
    );
    let prefix_digest = i14_terminal_prefix_digest_v2(
        logical_execution_receipt_digest,
        trace.logical_execution_verification_receipt_digest,
        trace.cancellation_card,
        &canonical_prefix,
        &canonical_requests,
    );
    Ok(I14TerminalTraceOutcomeV2::Frontier(
        I14TerminalFrontierCertificateV2 {
            logical_execution_receipt_digest,
            logical_execution_verification_receipt_digest: trace
                .logical_execution_verification_receipt_digest,
            cancellation_card: trace.cancellation_card,
            boundary_count: canonical_prefix.len(),
            last_boundary_ordinal: last.boundary.boundary_ordinal,
            last_boundary_logical_sequence: last.boundary.logical_sequence,
            prefix_digest,
        },
    ))
}

/// One lifecycle event with deterministic logical and calibrated time axes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct I14TimedLogicalEventV2 {
    /// Coordinator-assigned logical event sequence.
    pub logical_sequence: u64,
    /// Receipt-bound monotonic timestamp in nanoseconds.
    pub monotonic_ns: u64,
}

/// Closed cause taxonomy for a receipt-bound infrastructure-failure onset.
///
/// The first five sources are locally cross-checked against lifecycle evidence.
/// The final four bind protocol witnesses whose issuer, signature, capability,
/// trust policy, and revocation status are authenticated only by HELM/ledger.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14InfrastructureFailureSourceV2 {
    /// Independent watchdog coverage violated its admitted contract.
    WatchdogCoverage = 0,
    /// Scheduler before/after-tile polling evidence is incomplete.
    TilePollCoverage = 1,
    /// An admitted external child lacks heartbeat/termination coverage.
    ExternalHeartbeatCoverage = 2,
    /// One or more active descendants did not reach the drained cut.
    DescendantDrain = 3,
    /// A child was spawned at or after the frozen spawn frontier.
    SpawnAfterFrontierClosure = 4,
    /// Supervisor infrastructure failed outside the locally derived classes.
    Supervisor = 5,
    /// Authentication, trust-policy, or revocation verification failed.
    Authentication = 6,
    /// Drain protocol failed independently of descendant-count evidence.
    DrainProtocol = 7,
    /// Atomic publication protocol failed.
    PublicationProtocol = 8,
}

impl I14InfrastructureFailureSourceV2 {
    const fn wire_tag(self) -> u8 {
        match self {
            Self::WatchdogCoverage => 0,
            Self::TilePollCoverage => 1,
            Self::ExternalHeartbeatCoverage => 2,
            Self::DescendantDrain => 3,
            Self::SpawnAfterFrontierClosure => 4,
            Self::Supervisor => 5,
            Self::Authentication => 6,
            Self::DrainProtocol => 7,
            Self::PublicationProtocol => 8,
        }
    }
}

/// Typed, receipt-bound first infrastructure-failure onset.
///
/// The event's logical sequence, source tag, and receipt identity are canonical.
/// Its calibrated monotonic time is receipt-bound telemetry and participates in
/// causal validation and SLO arithmetic without entering canonical identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14InfrastructureFailureOnsetV2 {
    /// Closed infrastructure-failure source.
    pub source: I14InfrastructureFailureSourceV2,
    /// First structurally admitted, receipt-bound onset event on both axes.
    pub event: I14TimedLogicalEventV2,
    /// Independent verification receipt identity; local code binds but does not
    /// authenticate this digest.
    pub verification_receipt_digest: ContentHash,
}

/// Deterministic event class that caused terminal drain to begin.
///
/// The variant is canonical semantics; calibrated observation, deadline, and
/// failure-onset times remain in telemetry. An observation-timeout trigger
/// references the request whose calibrated deadline expired rather than
/// inventing a logical event from wall time.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum I14DrainTriggerV2 {
    /// An in-scope request was observed no later than its deadline.
    CancellationObserved {
        /// Request whose observation was the earliest effective drain cut.
        request_id: u64,
    },
    /// An in-scope request passed its inclusive observation deadline without
    /// an on-time observation.
    ObservationTimeoutDrain {
        /// Request whose missed deadline caused timeout drain to begin.
        request_id: u64,
    },
    /// A structurally admitted, receipt-bound lifecycle infrastructure failure
    /// caused or coalesced with drain before any earlier effective cancellation
    /// observation or timeout deadline. HELM/ledger authenticates its receipt.
    InfrastructureFailure {
        /// Logical sequence of the first receipt-bound infrastructure onset.
        onset_logical_sequence: u64,
    },
    /// Drain began for a terminal cause unrelated to cancellation.
    NonCancellationDrain,
}

/// Scheduler evidence that cancellation observation or terminal drain closed
/// the child-spawn frontier.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14SpawnFrontierEvidenceV2 {
    /// Cancellation request governing an earlier observation cut, or `None`
    /// when the unconditional drain-start cut governs the audit.
    pub request_id: Option<u64>,
    /// Clock-free identity of scheduler child membership and lifecycle
    /// semantics. Raw child timing has a separate telemetry root.
    pub scheduler_semantic_trace_digest: ContentHash,
    /// Last child-spawn event, if any child was ever spawned.
    pub last_child_spawn: Option<I14TimedLogicalEventV2>,
    /// Count of child spawns at or after the governed observation/drain cut.
    pub post_frontier_spawn_count: u64,
}

/// Bounded summary of the complete independent watchdog trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14WatchdogCoverageV2 {
    /// Number of polls represented by the complete watchdog trace.
    pub poll_count: u64,
    /// First calibrated watchdog poll.
    pub first_poll_monotonic_ns: u64,
    /// Last calibrated watchdog poll.
    pub last_poll_monotonic_ns: u64,
    /// Maximum consecutive poll gap measured over the complete trace.
    pub maximum_poll_gap_ns: u64,
    /// Clock-free watchdog coverage semantics used by canonical identity.
    pub watchdog_semantic_trace_digest: ContentHash,
    /// Identity of the complete multi-kind raw watchdog observation stream
    /// retained in telemetry. Poll records within it derive this summary.
    pub watchdog_raw_trace_digest: ContentHash,
    /// Identity of the independent verification/adjudication receipt for the
    /// clock-free watchdog semantics. Signature verification is a HELM/ledger
    /// gate.
    pub watchdog_verification_receipt_digest: ContentHash,
}

/// Receipt-bound coverage of scheduler cancellation polls around logical tiles.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14TilePollCoverageV2 {
    /// Logical tiles admitted by the frozen partition.
    pub admitted_tile_count: u64,
    /// Tiles carrying both required before/after boundary polls.
    pub fully_bracketed_tile_count: u64,
    /// Whether the required poll before item zero was recorded.
    pub before_item_zero_poll_observed: bool,
    /// Clock-free tile-membership and poll-placement trace identity.
    pub tile_poll_semantic_trace_digest: ContentHash,
    /// Identity of the complete raw tile-poll trace retained in telemetry.
    pub tile_poll_raw_trace_digest: ContentHash,
    /// Independent verification/adjudication receipt identity for the
    /// clock-free tile-poll semantic projection.
    pub tile_poll_verification_receipt_digest: ContentHash,
}

/// Summary and independent receipt identity for admitted external children.
///
/// Caller-supplied counts and the maximum-gap summary are checked here only for
/// internal consistency. The receipt digest is bound into canonical identity,
/// but the consuming HELM/ledger promotion gate must authenticate its signature
/// and issuer authority and independently verify completeness, membership,
/// acknowledgements, publication atomicity, and the derived gap verdict.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14ExternalHeartbeatCoverageV2 {
    /// Number of external children admitted during the campaign.
    pub admitted_external_children: u32,
    /// Number whose entire admitted lifetime is present in the heartbeat trace.
    pub fully_covered_external_children: u32,
    /// Largest simultaneously active external-child population.
    pub maximum_concurrent_external_children: u16,
    /// Number of raw heartbeat observations represented by the evidence.
    pub heartbeat_count: u64,
    /// Largest consecutive heartbeat gap for any admitted external child.
    pub maximum_heartbeat_gap_ns: u64,
    /// Clock-free identity of child membership and derived coverage semantics.
    pub heartbeat_semantic_trace_digest: ContentHash,
    /// Identity of the complete raw heartbeat trace retained in telemetry.
    pub heartbeat_raw_trace_digest: ContentHash,
    /// Independent verification/adjudication receipt identity for the
    /// clock-free heartbeat semantic projection.
    pub heartbeat_verification_receipt_digest: ContentHash,
    /// Whether every external child produced its interrupt/terminate-drain ack.
    pub all_termination_acks_observed: bool,
    /// Caller-supplied structural assertion of no partial publication outside
    /// the atomic finalization transaction; HELM/ledger verifies its provenance
    /// and truth before the projection gains authority.
    pub atomic_publication_verified: bool,
}

/// Complete drain/finalize evidence paired with a selected V2 terminal trace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14TerminalLifecycleTraceV2 {
    /// Execution start event.
    pub execution_started: I14TimedLogicalEventV2,
    /// Deterministic event class that caused drain to begin.
    pub drain_trigger: I14DrainTriggerV2,
    /// Drain-start event.
    pub drain_started: I14TimedLogicalEventV2,
    /// Descendant-drained event.
    pub drained: I14TimedLogicalEventV2,
    /// Finalization event, strictly before the terminal receipt boundary.
    pub finalized: I14TimedLogicalEventV2,
    /// Child scopes active when drain began.
    pub active_children_at_drain_start: u16,
    /// Child scopes accounted for by the drained event.
    pub drained_children: u16,
    /// Unconditional clock-free identity of child membership, causal order,
    /// outcomes, and losing-race semantics, including non-cancellation paths.
    pub child_lifecycle_semantic_trace_digest: ContentHash,
    /// Identity of the complete raw child lifecycle and timing trace retained
    /// in telemetry, including spawn, interrupt, acknowledgement, and drain
    /// times. It is deliberately excluded from canonical identity.
    pub child_lifecycle_raw_trace_digest: ContentHash,
    /// Independent verification/adjudication receipt identity for the clock-free
    /// child semantic projection. HELM/ledger authenticates the corresponding
    /// receipt and separately verifies raw-trace derivation; this layer only
    /// binds the semantic, raw-trace, and receipt identities.
    pub child_lifecycle_verification_receipt_digest: ContentHash,
    /// Mandatory spawn-frontier audit for every terminal path.
    pub spawn_frontier_audit: Option<I14SpawnFrontierEvidenceV2>,
    /// Coverage summary for the independent watchdog.
    pub watchdog_coverage: I14WatchdogCoverageV2,
    /// Coverage of before-item-zero and before/after-tile scheduler polls.
    pub tile_poll_coverage: I14TilePollCoverageV2,
    /// Coverage summary for admitted external-child heartbeats.
    pub external_heartbeat_coverage: I14ExternalHeartbeatCoverageV2,
    /// Earliest typed, receipt-bound infrastructure-failure onset. This is
    /// mandatory for an InfrastructureFailed terminal boundary, including
    /// generic supervisor/authentication/drain/publication protocol failures.
    pub first_infrastructure_failure_onset: Option<I14InfrastructureFailureOnsetV2>,
    /// Logical sequence of the first terminal-boundary record at or after the
    /// locally derived timeout onset. This is a latch-boundary identity, not a
    /// caller-invented logical event at the calibrated onset instant.
    pub first_timeout_failure_onset_logical_sequence: Option<u64>,
    /// Locally derived first nanosecond outside the earliest violated inclusive
    /// timeout SLO. Presence must match the latch-boundary identity exactly;
    /// the calibrated value is telemetry.
    pub first_timeout_failure_onset_monotonic_ns: Option<u64>,
}

/// Terminal cause class whose asynchronous onset must latch at the first
/// subsequent logical boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14LifecycleCauseClassV2 {
    /// Watchdog, polling, external-child, descendant, or spawn-frontier failure.
    InfrastructureFailed,
    /// Drain-trigger-to-drained or drained-to-finalized deadline failure.
    TimedOut,
}

/// Failure classes that must be reflected in the selected terminal cause.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14LifecycleFailureV2 {
    /// The complete watchdog trace violated its tier quantum.
    WatchdogCoverage,
    /// Required before-item-zero or before/after-tile polls are incomplete.
    TilePollCoverage,
    /// One or more external children lack complete bounded-heartbeat coverage.
    ExternalHeartbeatCoverage,
    /// Not every active descendant reached the drained cut.
    DescendantDrain,
    /// A child spawn occurred at or after the observation/drain frontier cut.
    SpawnAfterFrontierClosure,
    /// The selected drain trigger did not reach drained within its tier SLO.
    TriggerToDrainedDeadline,
    /// Finalization did not complete within its tier SLO.
    DrainedToFinalizedDeadline,
}

/// Fail-closed refusal of drain/finalize/watchdog evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14LifecycleRefusalV2 {
    /// One trace event reuses a coordinator logical sequence.
    DuplicateTraceLogicalSequence {
        /// Reused sequence.
        logical_sequence: u64,
    },
    /// Execution/drain/finalization logical order is invalid.
    LifecycleLogicalOrder,
    /// Calibrated lifecycle time order is inconsistent with logical order.
    LifecycleTimestampOrder,
    /// Finalization is not strictly before the selected terminal boundary.
    FinalizedNotBeforeBoundary,
    /// The drain trigger differs from the deterministic participating request.
    DrainTriggerMismatch {
        /// Deterministically required trigger.
        expected: I14DrainTriggerV2,
        /// Supplied trigger.
        found: I14DrainTriggerV2,
    },
    /// No mandatory spawn-frontier audit was supplied.
    MissingSpawnAudit {
        /// Request governing an earlier observation cut, if any.
        request_id: Option<u64>,
    },
    /// The spawn audit names a different cancellation request.
    SpawnAuditRequestMismatch {
        /// Required request.
        expected: Option<u64>,
        /// Supplied request.
        found: Option<u64>,
    },
    /// The last-spawn sequence and post-cut count contradict one another.
    SpawnAuditInconsistent,
    /// The mandatory spawn-frontier audit does not bind the unconditional
    /// child-lifecycle semantic trace used by the terminal result.
    SpawnAuditTraceDigestMismatch,
    /// A retained last-child-spawn event falls outside execution/final-drain
    /// order or collides with another logical event.
    SpawnEventOrder,
    /// Lifecycle child count exceeds the card's admitted cap.
    TooManyActiveChildren {
        /// Supplied active-child count.
        count: usize,
        /// Tier cap.
        cap: usize,
    },
    /// More descendants are reported drained than were active at drain start.
    DrainedChildrenExceedActive {
        /// Child scopes active when drain began.
        active: u16,
        /// Child scopes claimed by the drained cut.
        drained: u16,
    },
    /// The final terminal boundary still reports live child scopes.
    TerminalBoundaryHasInFlightChildren {
        /// Live child count at the terminal boundary.
        count: usize,
    },
    /// A boundary after the descendant-drained cut still reports live child
    /// scopes, contradicting the cut's finality.
    BoundaryAfterDrainedHasInFlightChildren {
        /// Boundary whose child count contradicts the drained cut.
        boundary_ordinal: u64,
        /// Live child count reported at that boundary.
        count: usize,
    },
    /// No nonempty watchdog coverage summary was supplied.
    MissingWatchdogCoverage,
    /// Watchdog coverage endpoints are outside campaign/terminal order.
    WatchdogCoverageOrder,
    /// Poll count, span, and claimed maximum gap cannot coexist.
    WatchdogCoverageInconsistent {
        /// Number of polls in the claimed complete trace.
        poll_count: u64,
        /// Time from first to last poll.
        span_ns: u64,
        /// Claimed maximum consecutive gap.
        maximum_poll_gap_ns: u64,
    },
    /// Tile count and fully bracketed count contradict one another.
    TilePollCoverageInconsistent,
    /// External-child coverage counts or empty-trace semantics contradict.
    ExternalHeartbeatCoverageInconsistent,
    /// Logical and calibrated axes for one failure onset disagree on presence.
    FailureOnsetAxisPresenceMismatch {
        /// Cause class governed by the onset.
        cause: I14LifecycleCauseClassV2,
        /// Supplied earliest logical onset.
        onset_logical_sequence: Option<u64>,
        /// Supplied earliest calibrated onset time.
        onset_monotonic_ns: Option<u64>,
    },
    /// Presence of a lifecycle-derived failure and its required onset/latch
    /// witness disagree.
    FailureOnsetPresenceMismatch {
        /// Cause class governed by the onset.
        cause: I14LifecycleCauseClassV2,
        /// Whether summary validation derived at least one such failure.
        failure_present: bool,
        /// Supplied earliest onset.
        onset_logical_sequence: Option<u64>,
        /// Supplied calibrated onset time.
        onset_monotonic_ns: Option<u64>,
    },
    /// A caller-supplied timeout onset does not equal the first nanosecond
    /// outside the inclusive admitted SLO.
    FailureOnsetTimestampMismatch {
        /// Cause class governed by the onset.
        cause: I14LifecycleCauseClassV2,
        /// Exactly derived first failing nanosecond, if a breach exists.
        expected_monotonic_ns: Option<u64>,
        /// Supplied calibrated onset time.
        found_monotonic_ns: Option<u64>,
    },
    /// A timeout witness names something other than the first boundary at or
    /// after the locally derived calibrated onset.
    FailureOnsetLogicalSequenceMismatch {
        /// Cause class governed by the latch.
        cause: I14LifecycleCauseClassV2,
        /// Exactly derived first latch-boundary logical sequence, if a breach
        /// exists.
        expected_logical_sequence: Option<u64>,
        /// Supplied latch-boundary logical sequence.
        found_logical_sequence: Option<u64>,
    },
    /// Infrastructure onset names a locally derived source whose evidence is
    /// not actually failed.
    InfrastructureFailureSourceMismatch {
        /// Unsupported source/evidence pairing.
        source: I14InfrastructureFailureSourceV2,
    },
    /// Infrastructure witness uses an all-zero verification-receipt identity.
    InfrastructureFailureVerificationReceiptDigestZero {
        /// Source whose receipt identity is absent.
        source: I14InfrastructureFailureSourceV2,
    },
    /// A supplied failure onset lies outside execution-to-terminal order.
    FailureOnsetOrder {
        /// Cause class governed by the onset.
        cause: I14LifecycleCauseClassV2,
        /// Malformed onset sequence.
        onset_logical_sequence: u64,
        /// Malformed calibrated onset time.
        onset_monotonic_ns: u64,
    },
    /// A lifecycle failure was not latched at the first boundary on or after
    /// its onset.
    FailureNotLatchedAtFirstBoundary {
        /// Cause class governed by the onset.
        cause: I14LifecycleCauseClassV2,
        /// Bound infrastructure-onset or derived timeout-latch sequence.
        onset_logical_sequence: u64,
        /// Boundary that first had to expose the cause.
        required_boundary_ordinal: u64,
        /// Boundary actually selected by the supplied trace.
        selected_boundary_ordinal: u64,
    },
    /// Maximum concurrently active external children exceeds the tier cap.
    TooManyConcurrentExternalChildren {
        /// Supplied maximum concurrency.
        count: usize,
        /// Tier cap.
        cap: usize,
    },
    /// A real lifecycle failure was hidden by a favorable terminal cause.
    LifecycleFailureNotReflected {
        /// Failure proven by the lifecycle trace.
        failure: I14LifecycleFailureV2,
        /// Cause selected by the terminal trace.
        selected: I14ExecutionDisposition,
    },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct I14CanonicalLifecycleProjectionV2 {
    execution_started_logical_sequence: u64,
    drain_trigger: I14DrainTriggerV2,
    drain_started_logical_sequence: u64,
    drained_logical_sequence: u64,
    finalized_logical_sequence: u64,
    active_children_at_drain_start: u16,
    drained_children: u16,
    child_lifecycle_semantic_trace_digest: ContentHash,
    child_lifecycle_verification_receipt_digest: ContentHash,
    spawn_frontier_request_id: Option<u64>,
    last_child_spawn_logical_sequence: Option<u64>,
    post_frontier_spawn_count: u64,
    watchdog_semantic_trace_digest: ContentHash,
    watchdog_verification_receipt_digest: ContentHash,
    admitted_tile_count: u64,
    fully_bracketed_tile_count: u64,
    before_item_zero_poll_observed: bool,
    tile_poll_semantic_trace_digest: ContentHash,
    tile_poll_verification_receipt_digest: ContentHash,
    admitted_external_children: u32,
    fully_covered_external_children: u32,
    maximum_concurrent_external_children: u16,
    heartbeat_semantic_trace_digest: ContentHash,
    heartbeat_verification_receipt_digest: ContentHash,
    all_external_termination_acks_observed: bool,
    external_atomic_publication_verified: bool,
    first_infrastructure_failure_onset_logical_sequence: Option<u64>,
    first_infrastructure_failure_source: Option<I14InfrastructureFailureSourceV2>,
    infrastructure_failure_verification_receipt_digest: Option<ContentHash>,
    first_timeout_failure_onset_logical_sequence: Option<u64>,
    watchdog_slo_breached: bool,
    tile_poll_coverage_failed: bool,
    external_heartbeat_slo_breached: bool,
    descendant_drain_failed: bool,
    post_frontier_spawned: bool,
    trigger_to_drained_slo_breached: bool,
    drained_to_finalized_slo_breached: bool,
}

impl I14CanonicalLifecycleProjectionV2 {
    /// Logical execution-start sequence.
    #[must_use]
    pub const fn execution_started_logical_sequence(self) -> u64 {
        self.execution_started_logical_sequence
    }

    /// Mandatory typed event class that deterministically triggered drain.
    #[must_use]
    pub const fn drain_trigger(self) -> I14DrainTriggerV2 {
        self.drain_trigger
    }

    /// Logical drain-start sequence.
    #[must_use]
    pub const fn drain_started_logical_sequence(self) -> u64 {
        self.drain_started_logical_sequence
    }

    /// Logical descendant-drained sequence.
    #[must_use]
    pub const fn drained_logical_sequence(self) -> u64 {
        self.drained_logical_sequence
    }

    /// Logical finalization sequence.
    #[must_use]
    pub const fn finalized_logical_sequence(self) -> u64 {
        self.finalized_logical_sequence
    }

    /// Child scopes active when drain began.
    #[must_use]
    pub const fn active_children_at_drain_start(self) -> u16 {
        self.active_children_at_drain_start
    }

    /// Child scopes accounted for by the drained cut.
    #[must_use]
    pub const fn drained_children(self) -> u16 {
        self.drained_children
    }

    /// Clock-free child/losing-race lifecycle semantic identity.
    #[must_use]
    pub const fn child_lifecycle_semantic_trace_digest(self) -> ContentHash {
        self.child_lifecycle_semantic_trace_digest
    }

    /// Independent verification-receipt identity for the child trace.
    #[must_use]
    pub const fn child_lifecycle_verification_receipt_digest(self) -> ContentHash {
        self.child_lifecycle_verification_receipt_digest
    }

    /// Cancellation request governing the spawn cut, if observation precedes
    /// unconditional drain start.
    #[must_use]
    pub const fn spawn_frontier_request_id(self) -> Option<u64> {
        self.spawn_frontier_request_id
    }

    /// Last child-spawn logical sequence, if the trace contains one.
    #[must_use]
    pub const fn last_child_spawn_logical_sequence(self) -> Option<u64> {
        self.last_child_spawn_logical_sequence
    }

    /// Number of spawns at or after the observation/drain frontier cut.
    #[must_use]
    pub const fn post_frontier_spawn_count(self) -> u64 {
        self.post_frontier_spawn_count
    }

    /// Clock-free watchdog semantic trace identity.
    #[must_use]
    pub const fn watchdog_semantic_trace_digest(self) -> ContentHash {
        self.watchdog_semantic_trace_digest
    }

    /// Independent watchdog verification-receipt identity.
    #[must_use]
    pub const fn watchdog_verification_receipt_digest(self) -> ContentHash {
        self.watchdog_verification_receipt_digest
    }

    /// Logical tiles represented by scheduler-poll evidence.
    #[must_use]
    pub const fn admitted_tile_count(self) -> u64 {
        self.admitted_tile_count
    }

    /// Logical tiles carrying both required boundary polls.
    #[must_use]
    pub const fn fully_bracketed_tile_count(self) -> u64 {
        self.fully_bracketed_tile_count
    }

    /// Whether the scheduler recorded its required poll before item zero.
    #[must_use]
    pub const fn before_item_zero_poll_observed(self) -> bool {
        self.before_item_zero_poll_observed
    }

    /// Clock-free tile-poll semantic trace identity.
    #[must_use]
    pub const fn tile_poll_semantic_trace_digest(self) -> ContentHash {
        self.tile_poll_semantic_trace_digest
    }

    /// Independent tile-poll verification-receipt identity.
    #[must_use]
    pub const fn tile_poll_verification_receipt_digest(self) -> ContentHash {
        self.tile_poll_verification_receipt_digest
    }

    /// Admitted external-child population represented by heartbeat evidence.
    #[must_use]
    pub const fn admitted_external_children(self) -> u32 {
        self.admitted_external_children
    }

    /// External children with complete admitted-lifetime heartbeat coverage.
    #[must_use]
    pub const fn fully_covered_external_children(self) -> u32 {
        self.fully_covered_external_children
    }

    /// Largest simultaneous external-child population.
    #[must_use]
    pub const fn maximum_concurrent_external_children(self) -> u16 {
        self.maximum_concurrent_external_children
    }

    /// Clock-free external-heartbeat semantic trace identity.
    #[must_use]
    pub const fn heartbeat_semantic_trace_digest(self) -> ContentHash {
        self.heartbeat_semantic_trace_digest
    }

    /// Independent external-heartbeat verification-receipt identity.
    #[must_use]
    pub const fn heartbeat_verification_receipt_digest(self) -> ContentHash {
        self.heartbeat_verification_receipt_digest
    }

    /// Whether every external child produced its termination/drain ack.
    #[must_use]
    pub const fn all_external_termination_acks_observed(self) -> bool {
        self.all_external_termination_acks_observed
    }

    /// Caller-supplied structural assertion that external execution preserved
    /// atomic no-partial publication; HELM/ledger verifies it for authority.
    #[must_use]
    pub const fn external_atomic_publication_verified(self) -> bool {
        self.external_atomic_publication_verified
    }

    /// Earliest structurally admitted receipt-bound infrastructure onset.
    #[must_use]
    pub const fn first_infrastructure_failure_onset_logical_sequence(self) -> Option<u64> {
        self.first_infrastructure_failure_onset_logical_sequence
    }

    /// Closed source of the first infrastructure-failure witness.
    #[must_use]
    pub const fn first_infrastructure_failure_source(
        self,
    ) -> Option<I14InfrastructureFailureSourceV2> {
        self.first_infrastructure_failure_source
    }

    /// Bound independent verification-receipt identity for that witness.
    #[must_use]
    pub const fn infrastructure_failure_verification_receipt_digest(self) -> Option<ContentHash> {
        self.infrastructure_failure_verification_receipt_digest
    }

    /// First boundary logical sequence that latched the derived timeout onset.
    #[must_use]
    pub const fn first_timeout_failure_onset_logical_sequence(self) -> Option<u64> {
        self.first_timeout_failure_onset_logical_sequence
    }

    /// Whether watchdog coverage violated the admitted quantum.
    #[must_use]
    pub const fn watchdog_slo_breached(self) -> bool {
        self.watchdog_slo_breached
    }

    /// Whether required scheduler tile-poll coverage was incomplete.
    #[must_use]
    pub const fn tile_poll_coverage_failed(self) -> bool {
        self.tile_poll_coverage_failed
    }

    /// Whether external-heartbeat coverage was incomplete or too sparse.
    #[must_use]
    pub const fn external_heartbeat_slo_breached(self) -> bool {
        self.external_heartbeat_slo_breached
    }

    /// Whether descendant drainage failed.
    #[must_use]
    pub const fn descendant_drain_failed(self) -> bool {
        self.descendant_drain_failed
    }

    /// Whether a child spawned at or after the spawn frontier closed.
    #[must_use]
    pub const fn post_frontier_spawned(self) -> bool {
        self.post_frontier_spawned
    }

    /// Whether the selected trigger missed its drain deadline.
    #[must_use]
    pub const fn trigger_to_drained_slo_breached(self) -> bool {
        self.trigger_to_drained_slo_breached
    }

    /// Whether finalization missed its post-drain deadline.
    #[must_use]
    pub const fn drained_to_finalized_slo_breached(self) -> bool {
        self.drained_to_finalized_slo_breached
    }
}

fn i14_selected_disposition(decision: I14TerminalBoundaryDecisionV1) -> I14ExecutionDisposition {
    match decision {
        I14TerminalBoundaryDecisionV1::Selected { disposition, .. } => disposition,
        I14TerminalBoundaryDecisionV1::DeferredByCancellation { .. }
        | I14TerminalBoundaryDecisionV1::NotTerminal => {
            unreachable!("selected V2 proof always carries a terminal decision")
        }
    }
}

fn i14_expected_drain_trigger_v2(
    trace: I14TerminalBoundaryTraceV2<'_>,
    selected_boundary: I14TerminalBoundaryV1<'_>,
    drain_started: I14TimedLogicalEventV2,
    infrastructure_failure_onset: Option<I14TimedLogicalEventV2>,
) -> I14DrainTriggerV2 {
    let scope_ids = selected_boundary
        .scope_ancestry
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    // Tuple field 1 is a causal tie rank, deliberately distinct from the
    // canonical drain-trigger wire tag: infrastructure=0, observation=1,
    // timeout=2. Effective time precedes this rank; logical sequence and stable
    // identity order candidates only within one simultaneous cause class.
    let cancellation_candidates = trace.requests.iter().filter_map(|request| {
        if request.logical_sequence >= drain_started.logical_sequence
            || !scope_ids.contains(&request.scope_root)
        {
            return None;
        }
        match request.observation {
            Some(observation)
                if observation.logical_sequence < drain_started.logical_sequence
                    && observation.monotonic_ns <= request.observation_deadline_ns =>
            {
                Some((
                    observation.monotonic_ns,
                    1_u8,
                    observation.logical_sequence,
                    request.request_id,
                    I14DrainTriggerV2::CancellationObserved {
                        request_id: request.request_id,
                    },
                ))
            }
            _ if drain_started.monotonic_ns > request.observation_deadline_ns => {
                let timeout_onset_ns = request
                    .observation_deadline_ns
                    .checked_add(1)
                    .expect("a later drain-start timestamp proves deadline+1 is representable");
                Some((
                    timeout_onset_ns,
                    2_u8,
                    request.logical_sequence,
                    request.request_id,
                    I14DrainTriggerV2::ObservationTimeoutDrain {
                        request_id: request.request_id,
                    },
                ))
            }
            _ => None,
        }
    });
    let infrastructure_candidate = infrastructure_failure_onset
        .filter(|onset| {
            onset.logical_sequence <= drain_started.logical_sequence
                && onset.monotonic_ns <= drain_started.monotonic_ns
        })
        .map(|onset| {
            (
                onset.monotonic_ns,
                0_u8,
                onset.logical_sequence,
                0_u64,
                I14DrainTriggerV2::InfrastructureFailure {
                    onset_logical_sequence: onset.logical_sequence,
                },
            )
        });
    cancellation_candidates
        .chain(infrastructure_candidate)
        .min_by_key(|candidate| (candidate.0, candidate.1, candidate.2, candidate.3))
        .map_or(I14DrainTriggerV2::NonCancellationDrain, |candidate| {
            candidate.4
        })
}

fn i14_expected_spawn_audit_request_v2(
    trace: I14TerminalBoundaryTraceV2<'_>,
    selected_boundary: I14TerminalBoundaryV1<'_>,
    drain_started_logical_sequence: u64,
) -> Option<u64> {
    let scope_ids = selected_boundary
        .scope_ancestry
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    trace
        .requests
        .iter()
        .filter_map(|request| {
            let observation = request.observation?;
            (request.logical_sequence < drain_started_logical_sequence
                && observation.logical_sequence < drain_started_logical_sequence
                && scope_ids.contains(&request.scope_root))
            .then_some((observation.logical_sequence, request.request_id))
        })
        .min()
        .map(|(_, request_id)| request_id)
}

fn i14_failure_reflected_v2(
    failure: I14LifecycleFailureV2,
    selected_boundary: I14TerminalBoundaryV1<'_>,
    selected: I14ExecutionDisposition,
) -> bool {
    match failure {
        I14LifecycleFailureV2::WatchdogCoverage
        | I14LifecycleFailureV2::TilePollCoverage
        | I14LifecycleFailureV2::ExternalHeartbeatCoverage
        | I14LifecycleFailureV2::DescendantDrain
        | I14LifecycleFailureV2::SpawnAfterFrontierClosure => {
            selected_boundary.infrastructure_failed
                && selected == I14ExecutionDisposition::InfrastructureFailed
        }
        I14LifecycleFailureV2::TriggerToDrainedDeadline
        | I14LifecycleFailureV2::DrainedToFinalizedDeadline => {
            selected_boundary.timed_out
                && matches!(
                    selected,
                    I14ExecutionDisposition::TimedOut
                        | I14ExecutionDisposition::InfrastructureFailed
                )
        }
    }
}

fn i14_first_ns_outside_inclusive_slo(
    started_ns: u64,
    inclusive_cap_ns: u64,
    completed_ns: u64,
) -> Option<u64> {
    let cutoff_ns = started_ns.checked_add(inclusive_cap_ns)?;
    (completed_ns > cutoff_ns).then(|| {
        cutoff_ns
            .checked_add(1)
            .expect("a later observed completion proves cutoff+1 is representable")
    })
}

fn i14_validate_failure_onset_v2(
    trace: I14TerminalBoundaryTraceV2<'_>,
    selected_index: usize,
    execution_started: I14TimedLogicalEventV2,
    cause: I14LifecycleCauseClassV2,
    failure_present: bool,
    onset_logical_sequence: Option<u64>,
    onset_monotonic_ns: Option<u64>,
) -> Result<(), I14LifecycleRefusalV2> {
    if onset_logical_sequence.is_some() != onset_monotonic_ns.is_some() {
        return Err(I14LifecycleRefusalV2::FailureOnsetAxisPresenceMismatch {
            cause,
            onset_logical_sequence,
            onset_monotonic_ns,
        });
    }
    if failure_present != onset_logical_sequence.is_some() {
        return Err(I14LifecycleRefusalV2::FailureOnsetPresenceMismatch {
            cause,
            failure_present,
            onset_logical_sequence,
            onset_monotonic_ns,
        });
    }
    let (Some(onset_logical_sequence), Some(onset_monotonic_ns)) =
        (onset_logical_sequence, onset_monotonic_ns)
    else {
        return Ok(());
    };
    let selected_boundary = trace.boundaries[selected_index].boundary;
    if onset_logical_sequence <= execution_started.logical_sequence
        || onset_logical_sequence > selected_boundary.logical_sequence
        || onset_monotonic_ns < execution_started.monotonic_ns
        || onset_monotonic_ns > selected_boundary.monotonic_ns
    {
        return Err(I14LifecycleRefusalV2::FailureOnsetOrder {
            cause,
            onset_logical_sequence,
            onset_monotonic_ns,
        });
    }
    let first_boundary_index = match cause {
        I14LifecycleCauseClassV2::InfrastructureFailed => trace
            .boundaries
            .iter()
            .position(|record| record.boundary.logical_sequence >= onset_logical_sequence),
        I14LifecycleCauseClassV2::TimedOut => trace
            .boundaries
            .iter()
            .position(|record| record.boundary.monotonic_ns >= onset_monotonic_ns),
    }
    .expect("selected boundary is at or after validated onset");
    if cause == I14LifecycleCauseClassV2::TimedOut
        && trace.boundaries[first_boundary_index]
            .boundary
            .logical_sequence
            < onset_logical_sequence
    {
        // Deadline expiry is ordered before a boundary stamped at the same
        // nanosecond. A later logical onset cannot retroactively label it.
        return Err(I14LifecycleRefusalV2::FailureOnsetOrder {
            cause,
            onset_logical_sequence,
            onset_monotonic_ns,
        });
    }
    if first_boundary_index != selected_index {
        return Err(I14LifecycleRefusalV2::FailureNotLatchedAtFirstBoundary {
            cause,
            onset_logical_sequence,
            required_boundary_ordinal: trace.boundaries[first_boundary_index]
                .boundary
                .boundary_ordinal,
            selected_boundary_ordinal: selected_boundary.boundary_ordinal,
        });
    }
    Ok(())
}

#[allow(clippy::too_many_lines)]
fn i14_validate_lifecycle_v2(
    trace: I14TerminalBoundaryTraceV2<'_>,
    selection: &I14FirstTerminalSelectionV2,
    lifecycle: I14TerminalLifecycleTraceV2,
) -> Result<I14CanonicalLifecycleProjectionV2, I14LifecycleRefusalV2> {
    let selected_record = trace.boundaries[selection.selected_index];
    let selected_boundary = selected_record.boundary;
    let genesis_boundary = trace.boundaries[0].boundary;
    if lifecycle.execution_started.logical_sequence >= genesis_boundary.logical_sequence {
        return Err(I14LifecycleRefusalV2::LifecycleLogicalOrder);
    }
    if lifecycle.execution_started.monotonic_ns > genesis_boundary.monotonic_ns {
        return Err(I14LifecycleRefusalV2::LifecycleTimestampOrder);
    }
    if !(lifecycle.execution_started.logical_sequence < lifecycle.drain_started.logical_sequence
        && lifecycle.drain_started.logical_sequence < lifecycle.drained.logical_sequence
        && lifecycle.drained.logical_sequence < lifecycle.finalized.logical_sequence)
    {
        return Err(I14LifecycleRefusalV2::LifecycleLogicalOrder);
    }
    if !(lifecycle.execution_started.monotonic_ns <= lifecycle.drain_started.monotonic_ns
        && lifecycle.drain_started.monotonic_ns <= lifecycle.drained.monotonic_ns
        && lifecycle.drained.monotonic_ns <= lifecycle.finalized.monotonic_ns)
    {
        return Err(I14LifecycleRefusalV2::LifecycleTimestampOrder);
    }
    if lifecycle.finalized.logical_sequence >= selected_boundary.logical_sequence
        || lifecycle.finalized.monotonic_ns > selected_boundary.monotonic_ns
    {
        return Err(I14LifecycleRefusalV2::FinalizedNotBeforeBoundary);
    }
    if selected_record.in_flight_children != 0 {
        return Err(I14LifecycleRefusalV2::TerminalBoundaryHasInFlightChildren {
            count: usize::from(selected_record.in_flight_children),
        });
    }
    if let Some(record) = trace.boundaries.iter().find(|record| {
        record.boundary.logical_sequence > lifecycle.drained.logical_sequence
            && record.in_flight_children != 0
    }) {
        return Err(
            I14LifecycleRefusalV2::BoundaryAfterDrainedHasInFlightChildren {
                boundary_ordinal: record.boundary.boundary_ordinal,
                count: usize::from(record.in_flight_children),
            },
        );
    }

    let mut trace_events = Vec::new();
    for record in trace.boundaries {
        trace_events.push((
            record.boundary.logical_sequence,
            record.boundary.monotonic_ns,
        ));
    }
    for request in trace.requests {
        trace_events.push((request.logical_sequence, request.requested_monotonic_ns));
        if let Some(observation) = request.observation {
            trace_events.push((observation.logical_sequence, observation.monotonic_ns));
        }
    }
    for event in [
        lifecycle.execution_started,
        lifecycle.drain_started,
        lifecycle.drained,
        lifecycle.finalized,
    ] {
        trace_events.push((event.logical_sequence, event.monotonic_ns));
    }
    if let Some(last_spawn) = lifecycle
        .spawn_frontier_audit
        .and_then(|audit| audit.last_child_spawn)
    {
        if !(lifecycle.execution_started.logical_sequence < last_spawn.logical_sequence
            && last_spawn.logical_sequence < lifecycle.drained.logical_sequence)
            || !(lifecycle.execution_started.monotonic_ns <= last_spawn.monotonic_ns
                && last_spawn.monotonic_ns <= lifecycle.drained.monotonic_ns)
        {
            return Err(I14LifecycleRefusalV2::SpawnEventOrder);
        }
        trace_events.push((last_spawn.logical_sequence, last_spawn.monotonic_ns));
    }
    let infrastructure_failure_onset = lifecycle
        .first_infrastructure_failure_onset
        .map(|witness| witness.event);
    if let Some(witness) = lifecycle.first_infrastructure_failure_onset {
        if witness
            .verification_receipt_digest
            .as_bytes()
            .iter()
            .all(|byte| *byte == 0)
        {
            return Err(
                I14LifecycleRefusalV2::InfrastructureFailureVerificationReceiptDigestZero {
                    source: witness.source,
                },
            );
        }
    }
    let timeout_failure_onset = match (
        lifecycle.first_timeout_failure_onset_logical_sequence,
        lifecycle.first_timeout_failure_onset_monotonic_ns,
    ) {
        (None, None) => None,
        (Some(logical_sequence), Some(monotonic_ns)) => Some(I14TimedLogicalEventV2 {
            logical_sequence,
            monotonic_ns,
        }),
        (onset_logical_sequence, onset_monotonic_ns) => {
            return Err(I14LifecycleRefusalV2::FailureOnsetAxisPresenceMismatch {
                cause: I14LifecycleCauseClassV2::TimedOut,
                onset_logical_sequence,
                onset_monotonic_ns,
            });
        }
    };
    // Infrastructure onset is a real receipt-bound logical event. It may share
    // identity only with the exact drain-start event when failure onset and
    // drain start are one coalesced transaction. Equal timestamps never make
    // an alias with a request, observation, boundary, spawn, drained, finalized
    // or execution-started event valid. Timeout onset is calibrated arithmetic
    // between events; its logical field names the first boundary that must
    // latch it and therefore must not be inserted as a second event colliding
    // with that boundary.
    for onset in infrastructure_failure_onset {
        if onset.logical_sequence == lifecycle.drain_started.logical_sequence {
            if onset.monotonic_ns != lifecycle.drain_started.monotonic_ns {
                return Err(I14LifecycleRefusalV2::LifecycleTimestampOrder);
            }
        } else if trace_events
            .iter()
            .any(|(logical_sequence, _)| *logical_sequence == onset.logical_sequence)
        {
            return Err(I14LifecycleRefusalV2::DuplicateTraceLogicalSequence {
                logical_sequence: onset.logical_sequence,
            });
        } else {
            trace_events.push((onset.logical_sequence, onset.monotonic_ns));
        }
    }
    trace_events.sort_unstable_by_key(|(logical_sequence, _)| *logical_sequence);
    for pair in trace_events.windows(2) {
        if pair[0].0 == pair[1].0 {
            return Err(I14LifecycleRefusalV2::DuplicateTraceLogicalSequence {
                logical_sequence: pair[0].0,
            });
        }
        if pair[0].1 > pair[1].1 {
            return Err(I14LifecycleRefusalV2::LifecycleTimestampOrder);
        }
    }

    let expected_trigger = i14_expected_drain_trigger_v2(
        trace,
        selected_boundary,
        lifecycle.drain_started,
        infrastructure_failure_onset,
    );
    if lifecycle.drain_trigger != expected_trigger {
        return Err(I14LifecycleRefusalV2::DrainTriggerMismatch {
            expected: expected_trigger,
            found: lifecycle.drain_trigger,
        });
    }
    let expected_spawn_audit_request = i14_expected_spawn_audit_request_v2(
        trace,
        selected_boundary,
        lifecycle.drain_started.logical_sequence,
    );
    let audit = lifecycle
        .spawn_frontier_audit
        .ok_or(I14LifecycleRefusalV2::MissingSpawnAudit {
            request_id: expected_spawn_audit_request,
        })?;
    if audit.request_id != expected_spawn_audit_request {
        return Err(I14LifecycleRefusalV2::SpawnAuditRequestMismatch {
            expected: expected_spawn_audit_request,
            found: audit.request_id,
        });
    }
    if audit.scheduler_semantic_trace_digest != lifecycle.child_lifecycle_semantic_trace_digest {
        return Err(I14LifecycleRefusalV2::SpawnAuditTraceDigestMismatch);
    }
    let spawn_cut_observation = expected_spawn_audit_request
        .and_then(|request_id| {
            trace
                .requests
                .iter()
                .find(|request| request.request_id == request_id)
        })
        .and_then(|request| request.observation)
        .filter(|observation| {
            observation.logical_sequence < lifecycle.drain_started.logical_sequence
        });
    // The spawn frontier closes at the first actual observation (on time or
    // late) before drain, otherwise at drain start. Request submission and a
    // calibrated deadline alone never synthesize a logical spawn cut.
    let spawn_frontier_closed_at = spawn_cut_observation
        .map_or(lifecycle.drain_started.logical_sequence, |observation| {
            observation.logical_sequence
        });
    let post_frontier_spawned = audit
        .last_child_spawn
        .is_some_and(|event| event.logical_sequence >= spawn_frontier_closed_at);
    if post_frontier_spawned != (audit.post_frontier_spawn_count > 0) {
        return Err(I14LifecycleRefusalV2::SpawnAuditInconsistent);
    }
    let last_child_spawn_logical_sequence =
        audit.last_child_spawn.map(|event| event.logical_sequence);
    let post_frontier_spawn_count = audit.post_frontier_spawn_count;

    let child_cap = trace.cancellation_card.tier().child_and_observer_cap();
    let largest_lifecycle_child_count = usize::from(
        lifecycle
            .active_children_at_drain_start
            .max(lifecycle.drained_children),
    );
    if largest_lifecycle_child_count > child_cap {
        return Err(I14LifecycleRefusalV2::TooManyActiveChildren {
            count: largest_lifecycle_child_count,
            cap: child_cap,
        });
    }
    if lifecycle.drained_children > lifecycle.active_children_at_drain_start {
        return Err(I14LifecycleRefusalV2::DrainedChildrenExceedActive {
            active: lifecycle.active_children_at_drain_start,
            drained: lifecycle.drained_children,
        });
    }
    let descendant_drain_failed =
        lifecycle.drained_children < lifecycle.active_children_at_drain_start;

    let tile_polls = lifecycle.tile_poll_coverage;
    if tile_polls.fully_bracketed_tile_count > tile_polls.admitted_tile_count
        || (tile_polls.admitted_tile_count == 0 && tile_polls.before_item_zero_poll_observed)
    {
        return Err(I14LifecycleRefusalV2::TilePollCoverageInconsistent);
    }
    let tile_poll_coverage_failed = tile_polls.admitted_tile_count > 0
        && (tile_polls.fully_bracketed_tile_count != tile_polls.admitted_tile_count
            || !tile_polls.before_item_zero_poll_observed);

    let external = lifecycle.external_heartbeat_coverage;
    if usize::from(external.maximum_concurrent_external_children) > child_cap {
        return Err(I14LifecycleRefusalV2::TooManyConcurrentExternalChildren {
            count: usize::from(external.maximum_concurrent_external_children),
            cap: child_cap,
        });
    }
    let empty_external_trace = external.admitted_external_children == 0;
    let external_counts_inconsistent = external.fully_covered_external_children
        > external.admitted_external_children
        || u32::from(external.maximum_concurrent_external_children)
            > external.admitted_external_children
        || (trace
            .cancellation_card
            .external_child_catalog_digest()
            .is_none()
            != empty_external_trace)
        || (empty_external_trace
            && (external.fully_covered_external_children != 0
                || external.maximum_concurrent_external_children != 0
                || external.heartbeat_count != 0
                || external.maximum_heartbeat_gap_ns != 0
                || !external.all_termination_acks_observed
                || !external.atomic_publication_verified))
        || (!empty_external_trace && external.maximum_concurrent_external_children == 0)
        || u64::from(external.fully_covered_external_children) > external.heartbeat_count;
    if external_counts_inconsistent {
        return Err(I14LifecycleRefusalV2::ExternalHeartbeatCoverageInconsistent);
    }
    let external_heartbeat_slo_breached = external.fully_covered_external_children
        != external.admitted_external_children
        || external.heartbeat_count < u64::from(external.admitted_external_children)
        || external.maximum_heartbeat_gap_ns
            > trace.cancellation_card.external_heartbeat_bound_ns()
        || (!empty_external_trace
            && (!external.all_termination_acks_observed || !external.atomic_publication_verified));

    let coverage = lifecycle.watchdog_coverage;
    if coverage.poll_count == 0 {
        return Err(I14LifecycleRefusalV2::MissingWatchdogCoverage);
    }
    if lifecycle.execution_started.monotonic_ns < trace.campaign_started_monotonic_ns
        || coverage.first_poll_monotonic_ns < trace.campaign_started_monotonic_ns
        || coverage.last_poll_monotonic_ns < coverage.first_poll_monotonic_ns
        || coverage.last_poll_monotonic_ns > selected_boundary.monotonic_ns
        || coverage.last_poll_monotonic_ns != selected_record.last_watchdog_poll_monotonic_ns
    {
        return Err(I14LifecycleRefusalV2::WatchdogCoverageOrder);
    }
    let watchdog_span_ns = coverage.last_poll_monotonic_ns - coverage.first_poll_monotonic_ns;
    let represented_span_ns = u128::from(coverage.poll_count.saturating_sub(1))
        * u128::from(coverage.maximum_poll_gap_ns);
    if (coverage.poll_count == 1 && (watchdog_span_ns != 0 || coverage.maximum_poll_gap_ns != 0))
        || (coverage.poll_count > 1 && coverage.maximum_poll_gap_ns > watchdog_span_ns)
        || u128::from(watchdog_span_ns) > represented_span_ns
    {
        return Err(I14LifecycleRefusalV2::WatchdogCoverageInconsistent {
            poll_count: coverage.poll_count,
            span_ns: watchdog_span_ns,
            maximum_poll_gap_ns: coverage.maximum_poll_gap_ns,
        });
    }
    let watchdog_cap_ns = trace.cancellation_card.tier().watchdog_quantum_ns();
    let watchdog_slo_breached =
        coverage.first_poll_monotonic_ns - trace.campaign_started_monotonic_ns > watchdog_cap_ns
            || coverage.maximum_poll_gap_ns > watchdog_cap_ns
            || selected_boundary.monotonic_ns - coverage.last_poll_monotonic_ns > watchdog_cap_ns;

    let trigger_monotonic_ns = match lifecycle.drain_trigger {
        I14DrainTriggerV2::CancellationObserved { request_id } => {
            trace
                .requests
                .iter()
                .find(|request| request.request_id == request_id)
                .and_then(|request| request.observation)
                .expect("validated observed trigger references its trace observation")
                .monotonic_ns
        }
        I14DrainTriggerV2::ObservationTimeoutDrain { request_id } => trace
            .requests
            .iter()
            .find(|request| request.request_id == request_id)
            .expect("validated timeout trigger references its trace request")
            .observation_deadline_ns
            .checked_add(1)
            .expect("validated timeout drain proves deadline+1 is representable"),
        I14DrainTriggerV2::InfrastructureFailure {
            onset_logical_sequence,
        } => {
            let onset = infrastructure_failure_onset
                .expect("validated infrastructure trigger references its onset");
            debug_assert_eq!(onset.logical_sequence, onset_logical_sequence);
            onset.monotonic_ns
        }
        I14DrainTriggerV2::NonCancellationDrain => lifecycle.drain_started.monotonic_ns,
    };
    let trigger_to_drained_failure_onset_ns = i14_first_ns_outside_inclusive_slo(
        trigger_monotonic_ns,
        trace.cancellation_card.tier().trigger_to_drained_ns(),
        lifecycle.drained.monotonic_ns,
    );
    let drained_to_finalized_failure_onset_ns = i14_first_ns_outside_inclusive_slo(
        lifecycle.drained.monotonic_ns,
        trace.cancellation_card.tier().drained_to_finalized_ns(),
        lifecycle.finalized.monotonic_ns,
    );
    let trigger_to_drained_slo_breached = trigger_to_drained_failure_onset_ns.is_some();
    let drained_to_finalized_slo_breached = drained_to_finalized_failure_onset_ns.is_some();
    let expected_timeout_failure_onset_ns = trigger_to_drained_failure_onset_ns
        .into_iter()
        .chain(drained_to_finalized_failure_onset_ns)
        .min();
    let found_timeout_failure_onset_ns = timeout_failure_onset.map(|onset| onset.monotonic_ns);
    if found_timeout_failure_onset_ns != expected_timeout_failure_onset_ns {
        return Err(I14LifecycleRefusalV2::FailureOnsetTimestampMismatch {
            cause: I14LifecycleCauseClassV2::TimedOut,
            expected_monotonic_ns: expected_timeout_failure_onset_ns,
            found_monotonic_ns: found_timeout_failure_onset_ns,
        });
    }
    let expected_timeout_latch_logical_sequence =
        expected_timeout_failure_onset_ns.and_then(|onset_ns| {
            trace
                .boundaries
                .iter()
                .find(|record| record.boundary.monotonic_ns >= onset_ns)
                .map(|record| record.boundary.logical_sequence)
        });
    let found_timeout_latch_logical_sequence =
        timeout_failure_onset.map(|onset| onset.logical_sequence);
    if found_timeout_latch_logical_sequence != expected_timeout_latch_logical_sequence {
        return Err(I14LifecycleRefusalV2::FailureOnsetLogicalSequenceMismatch {
            cause: I14LifecycleCauseClassV2::TimedOut,
            expected_logical_sequence: expected_timeout_latch_logical_sequence,
            found_logical_sequence: found_timeout_latch_logical_sequence,
        });
    }

    let derived_infrastructure_failure_present = watchdog_slo_breached
        || tile_poll_coverage_failed
        || external_heartbeat_slo_breached
        || descendant_drain_failed
        || post_frontier_spawned;
    let infrastructure_failure_present =
        derived_infrastructure_failure_present || selected_boundary.infrastructure_failed;
    if let Some(witness) = lifecycle.first_infrastructure_failure_onset {
        let source_matches = match witness.source {
            I14InfrastructureFailureSourceV2::WatchdogCoverage => watchdog_slo_breached,
            I14InfrastructureFailureSourceV2::TilePollCoverage => tile_poll_coverage_failed,
            I14InfrastructureFailureSourceV2::ExternalHeartbeatCoverage => {
                external_heartbeat_slo_breached
            }
            I14InfrastructureFailureSourceV2::DescendantDrain => descendant_drain_failed,
            I14InfrastructureFailureSourceV2::SpawnAfterFrontierClosure => post_frontier_spawned,
            I14InfrastructureFailureSourceV2::Supervisor
            | I14InfrastructureFailureSourceV2::Authentication
            | I14InfrastructureFailureSourceV2::DrainProtocol
            | I14InfrastructureFailureSourceV2::PublicationProtocol => true,
        };
        if !source_matches {
            return Err(I14LifecycleRefusalV2::InfrastructureFailureSourceMismatch {
                source: witness.source,
            });
        }
    }
    let timeout_failure_present = expected_timeout_failure_onset_ns.is_some();
    i14_validate_failure_onset_v2(
        trace,
        selection.selected_index,
        lifecycle.execution_started,
        I14LifecycleCauseClassV2::InfrastructureFailed,
        infrastructure_failure_present,
        lifecycle
            .first_infrastructure_failure_onset
            .map(|witness| witness.event.logical_sequence),
        lifecycle
            .first_infrastructure_failure_onset
            .map(|witness| witness.event.monotonic_ns),
    )?;
    i14_validate_failure_onset_v2(
        trace,
        selection.selected_index,
        lifecycle.execution_started,
        I14LifecycleCauseClassV2::TimedOut,
        timeout_failure_present,
        lifecycle.first_timeout_failure_onset_logical_sequence,
        lifecycle.first_timeout_failure_onset_monotonic_ns,
    )?;
    let selected = i14_selected_disposition(selection.decision);
    for (failure, present) in [
        (
            I14LifecycleFailureV2::WatchdogCoverage,
            watchdog_slo_breached,
        ),
        (
            I14LifecycleFailureV2::TilePollCoverage,
            tile_poll_coverage_failed,
        ),
        (
            I14LifecycleFailureV2::ExternalHeartbeatCoverage,
            external_heartbeat_slo_breached,
        ),
        (
            I14LifecycleFailureV2::DescendantDrain,
            descendant_drain_failed,
        ),
        (
            I14LifecycleFailureV2::SpawnAfterFrontierClosure,
            post_frontier_spawned,
        ),
        (
            I14LifecycleFailureV2::TriggerToDrainedDeadline,
            trigger_to_drained_slo_breached,
        ),
        (
            I14LifecycleFailureV2::DrainedToFinalizedDeadline,
            drained_to_finalized_slo_breached,
        ),
    ] {
        if present && !i14_failure_reflected_v2(failure, selected_boundary, selected) {
            return Err(I14LifecycleRefusalV2::LifecycleFailureNotReflected { failure, selected });
        }
    }

    Ok(I14CanonicalLifecycleProjectionV2 {
        execution_started_logical_sequence: lifecycle.execution_started.logical_sequence,
        drain_trigger: lifecycle.drain_trigger,
        drain_started_logical_sequence: lifecycle.drain_started.logical_sequence,
        drained_logical_sequence: lifecycle.drained.logical_sequence,
        finalized_logical_sequence: lifecycle.finalized.logical_sequence,
        active_children_at_drain_start: lifecycle.active_children_at_drain_start,
        drained_children: lifecycle.drained_children,
        child_lifecycle_semantic_trace_digest: lifecycle.child_lifecycle_semantic_trace_digest,
        child_lifecycle_verification_receipt_digest: lifecycle
            .child_lifecycle_verification_receipt_digest,
        spawn_frontier_request_id: audit.request_id,
        last_child_spawn_logical_sequence,
        post_frontier_spawn_count,
        watchdog_semantic_trace_digest: coverage.watchdog_semantic_trace_digest,
        watchdog_verification_receipt_digest: coverage.watchdog_verification_receipt_digest,
        admitted_tile_count: tile_polls.admitted_tile_count,
        fully_bracketed_tile_count: tile_polls.fully_bracketed_tile_count,
        before_item_zero_poll_observed: tile_polls.before_item_zero_poll_observed,
        tile_poll_semantic_trace_digest: tile_polls.tile_poll_semantic_trace_digest,
        tile_poll_verification_receipt_digest: tile_polls.tile_poll_verification_receipt_digest,
        admitted_external_children: external.admitted_external_children,
        fully_covered_external_children: external.fully_covered_external_children,
        maximum_concurrent_external_children: external.maximum_concurrent_external_children,
        heartbeat_semantic_trace_digest: external.heartbeat_semantic_trace_digest,
        heartbeat_verification_receipt_digest: external.heartbeat_verification_receipt_digest,
        all_external_termination_acks_observed: external.all_termination_acks_observed,
        external_atomic_publication_verified: external.atomic_publication_verified,
        first_infrastructure_failure_onset_logical_sequence: lifecycle
            .first_infrastructure_failure_onset
            .map(|witness| witness.event.logical_sequence),
        first_infrastructure_failure_source: lifecycle
            .first_infrastructure_failure_onset
            .map(|witness| witness.source),
        infrastructure_failure_verification_receipt_digest: lifecycle
            .first_infrastructure_failure_onset
            .map(|witness| witness.verification_receipt_digest),
        first_timeout_failure_onset_logical_sequence: lifecycle
            .first_timeout_failure_onset_logical_sequence,
        watchdog_slo_breached,
        tile_poll_coverage_failed,
        external_heartbeat_slo_breached,
        descendant_drain_failed,
        post_frontier_spawned,
        trigger_to_drained_slo_breached,
        drained_to_finalized_slo_breached,
    })
}

/// Inputs to the canonical terminal-result constructor used by promotion gates.
#[derive(Clone, Copy, Debug)]
pub struct I14CanonicalTerminalResultInputV2<'a> {
    /// Complete genesis-to-first-terminal trace.
    pub trace: I14TerminalBoundaryTraceV2<'a>,
    /// Drain, finalization, spawn-frontier, and watchdog evidence.
    pub lifecycle: I14TerminalLifecycleTraceV2,
    /// Raw terminal-status tuple retained and fail-closed-normalized.
    pub terminal_status: I14TerminalStatusV1,
    /// Content identity of the scientific/operational result payload.
    pub semantic_payload_digest: ContentHash,
}

/// Fail-closed authoritative canonical-result refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14CanonicalResultRefusalV2 {
    /// The complete first-terminal trace is malformed.
    TerminalTrace(I14TerminalTraceRefusalV2),
    /// The supplied trace is a valid nonterminal frontier, not a result.
    TraceAtNonterminalFrontier {
        /// Number of validated frontier boundaries.
        boundary_count: usize,
        /// Digest of the complete nonterminal prefix.
        prefix_digest: ContentHash,
    },
    /// Drain/finalize/watchdog evidence is malformed or hides a failure.
    Lifecycle(I14LifecycleRefusalV2),
    /// The selected boundary cannot form its V1 local-result projection.
    LocalResult(I14CanonicalResultRefusalV1),
}

/// Validated clock-free I14 terminal projection.
///
/// This type is not by itself a signed promotion receipt. It becomes an
/// authority-bearing gate input only after the consuming HELM/ledger verifier
/// authenticates the receipt identified by
/// `logical_execution_verification_receipt_digest` against the frozen issuer,
/// capability, revocation, card, and semantic-work policy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct I14CanonicalTerminalResultV2 {
    logical_execution_receipt_digest: ContentHash,
    logical_execution_verification_receipt_digest: ContentHash,
    terminal_prefix_digest: ContentHash,
    cancellation_card: I14CancellationCardV2,
    selected_index: usize,
    boundary_prefix: Vec<I14CanonicalBoundaryRecordV2>,
    lifecycle: I14CanonicalLifecycleProjectionV2,
    local_result: I14CanonicalTerminalResultV1,
}

impl I14CanonicalTerminalResultV2 {
    /// Recomputed clock-free content identity of the complete trace.
    #[must_use]
    pub const fn logical_execution_receipt_digest(&self) -> ContentHash {
        self.logical_execution_receipt_digest
    }

    /// Independent verification-receipt identity for the logical stream.
    #[must_use]
    pub const fn logical_execution_verification_receipt_digest(&self) -> ContentHash {
        self.logical_execution_verification_receipt_digest
    }

    /// Digest of the complete request-inclusive first-terminal prefix.
    #[must_use]
    pub const fn terminal_prefix_digest(&self) -> ContentHash {
        self.terminal_prefix_digest
    }

    /// Validated cancellation contract.
    #[must_use]
    pub const fn cancellation_card(&self) -> I14CancellationCardV2 {
        self.cancellation_card
    }

    /// Zero-based first-terminal boundary position.
    #[must_use]
    pub const fn selected_index(&self) -> usize {
        self.selected_index
    }

    /// Count of canonical boundaries from genesis through first terminal.
    #[must_use]
    pub fn boundary_count(&self) -> usize {
        self.boundary_prefix.len()
    }

    /// Validated single-boundary semantic result at the selected frontier.
    #[must_use]
    pub const fn local_result(&self) -> &I14CanonicalTerminalResultV1 {
        &self.local_result
    }

    /// Validated clock-free lifecycle projection used by this result.
    #[must_use]
    pub const fn lifecycle(&self) -> &I14CanonicalLifecycleProjectionV2 {
        &self.lifecycle
    }

    /// Domain-separated V2 identity of card, prefix, lifecycle, and result.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut payload = b"I14_CANONICAL_TERMINAL_RESULT_V2\0".to_vec();
        payload.extend_from_slice(self.logical_execution_receipt_digest.as_bytes());
        payload.extend_from_slice(
            self.logical_execution_verification_receipt_digest
                .as_bytes(),
        );
        payload.extend_from_slice(self.terminal_prefix_digest.as_bytes());
        payload.extend_from_slice(self.cancellation_card.digest().as_bytes());
        i14_push_len(&mut payload, self.boundary_prefix.len());
        i14_push_u64(
            &mut payload,
            u64::try_from(self.selected_index).expect("I14 selected index fits u64"),
        );
        for record in &self.boundary_prefix {
            i14_push_canonical_boundary_v2(&mut payload, record);
        }
        i14_push_u64(
            &mut payload,
            self.lifecycle.execution_started_logical_sequence,
        );
        i14_push_drain_trigger(&mut payload, self.lifecycle.drain_trigger);
        i14_push_u64(&mut payload, self.lifecycle.drain_started_logical_sequence);
        i14_push_u64(&mut payload, self.lifecycle.drained_logical_sequence);
        i14_push_u64(&mut payload, self.lifecycle.finalized_logical_sequence);
        i14_push_u64(
            &mut payload,
            u64::from(self.lifecycle.active_children_at_drain_start),
        );
        i14_push_u64(&mut payload, u64::from(self.lifecycle.drained_children));
        payload.extend_from_slice(
            self.lifecycle
                .child_lifecycle_semantic_trace_digest
                .as_bytes(),
        );
        payload.extend_from_slice(
            self.lifecycle
                .child_lifecycle_verification_receipt_digest
                .as_bytes(),
        );
        i14_push_optional_u64(&mut payload, self.lifecycle.spawn_frontier_request_id);
        i14_push_optional_u64(
            &mut payload,
            self.lifecycle.last_child_spawn_logical_sequence,
        );
        i14_push_u64(&mut payload, self.lifecycle.post_frontier_spawn_count);
        payload.extend_from_slice(self.lifecycle.watchdog_semantic_trace_digest.as_bytes());
        payload.extend_from_slice(
            self.lifecycle
                .watchdog_verification_receipt_digest
                .as_bytes(),
        );
        i14_push_u64(&mut payload, self.lifecycle.admitted_tile_count);
        i14_push_u64(&mut payload, self.lifecycle.fully_bracketed_tile_count);
        payload.push(u8::from(self.lifecycle.before_item_zero_poll_observed));
        payload.extend_from_slice(self.lifecycle.tile_poll_semantic_trace_digest.as_bytes());
        payload.extend_from_slice(
            self.lifecycle
                .tile_poll_verification_receipt_digest
                .as_bytes(),
        );
        i14_push_u64(
            &mut payload,
            u64::from(self.lifecycle.admitted_external_children),
        );
        i14_push_u64(
            &mut payload,
            u64::from(self.lifecycle.fully_covered_external_children),
        );
        i14_push_u64(
            &mut payload,
            u64::from(self.lifecycle.maximum_concurrent_external_children),
        );
        payload.extend_from_slice(self.lifecycle.heartbeat_semantic_trace_digest.as_bytes());
        payload.extend_from_slice(
            self.lifecycle
                .heartbeat_verification_receipt_digest
                .as_bytes(),
        );
        i14_push_canonical_infrastructure_failure_onset_v2(
            &mut payload,
            self.lifecycle
                .first_infrastructure_failure_onset_logical_sequence,
            self.lifecycle.first_infrastructure_failure_source,
            self.lifecycle
                .infrastructure_failure_verification_receipt_digest,
        );
        i14_push_optional_u64(
            &mut payload,
            self.lifecycle.first_timeout_failure_onset_logical_sequence,
        );
        payload.extend_from_slice(&[
            u8::from(self.lifecycle.all_external_termination_acks_observed),
            u8::from(self.lifecycle.external_atomic_publication_verified),
            u8::from(self.lifecycle.watchdog_slo_breached),
            u8::from(self.lifecycle.tile_poll_coverage_failed),
            u8::from(self.lifecycle.external_heartbeat_slo_breached),
            u8::from(self.lifecycle.descendant_drain_failed),
            u8::from(self.lifecycle.post_frontier_spawned),
            u8::from(self.lifecycle.trigger_to_drained_slo_breached),
            u8::from(self.lifecycle.drained_to_finalized_slo_breached),
        ]);
        payload.extend_from_slice(self.local_result.digest().as_bytes());
        hash_domain(TERMINAL_RESULT_DOMAIN_V2, &payload)
    }
}

/// Validate and construct an authoritative clock-free terminal result.
pub fn i14_canonical_terminal_result_v2(
    input: I14CanonicalTerminalResultInputV2<'_>,
) -> Result<I14CanonicalTerminalResultV2, I14CanonicalResultRefusalV2> {
    let selection = match i14_select_first_terminal_boundary_v2(input.trace)
        .map_err(I14CanonicalResultRefusalV2::TerminalTrace)?
    {
        I14TerminalTraceOutcomeV2::Selected(selection) => selection,
        I14TerminalTraceOutcomeV2::Frontier(frontier) => {
            return Err(I14CanonicalResultRefusalV2::TraceAtNonterminalFrontier {
                boundary_count: frontier.boundary_count,
                prefix_digest: frontier.prefix_digest,
            });
        }
    };
    let lifecycle = i14_validate_lifecycle_v2(input.trace, &selection, input.lifecycle)
        .map_err(I14CanonicalResultRefusalV2::Lifecycle)?;
    let selected_boundary = input.trace.boundaries[selection.selected_index].boundary;
    let local_result = i14_canonical_terminal_result_v1(I14CanonicalTerminalResultInputV1 {
        boundary: selected_boundary,
        requests: input.trace.requests,
        terminal_status: input.terminal_status,
        semantic_payload_digest: input.semantic_payload_digest,
    })
    .map_err(I14CanonicalResultRefusalV2::LocalResult)?;
    Ok(I14CanonicalTerminalResultV2 {
        logical_execution_receipt_digest: selection.logical_execution_receipt_digest,
        logical_execution_verification_receipt_digest: selection
            .logical_execution_verification_receipt_digest,
        terminal_prefix_digest: selection.prefix_digest,
        cancellation_card: selection.cancellation_card,
        selected_index: selection.selected_index,
        boundary_prefix: selection.boundary_prefix,
        lifecycle,
        local_result,
    })
}

/// Validate and hash an authoritative clock-free terminal result.
pub fn i14_canonical_terminal_result_digest_v2(
    input: I14CanonicalTerminalResultInputV2<'_>,
) -> Result<ContentHash, I14CanonicalResultRefusalV2> {
    Ok(i14_canonical_terminal_result_v2(input)?.digest())
}

/// Receipt-bound summary of logical events arriving after terminal selection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct I14LateEventTailV2 {
    /// Number of verification-receipt-bound logical events retained after
    /// finalization. Receipt authentication is deferred to HELM/ledger.
    pub event_count: u64,
    /// First post-terminal logical sequence, or `None` for an empty tail.
    pub first_logical_sequence: Option<u64>,
    /// Last post-terminal logical sequence, or `None` for an empty tail.
    pub last_logical_sequence: Option<u64>,
    /// Clock-free identity of the retained post-terminal event tail.
    pub semantic_trace_digest: ContentHash,
    /// Independent verification/adjudication receipt identity for the tail.
    pub verification_receipt_digest: ContentHash,
}

/// Raw timing/calibration fields paired with an authoritative V2 result.
#[derive(Clone, Copy, Debug)]
pub struct I14TelemetryEnvelopeInputV2<'a> {
    /// Authoritative terminal input; it is revalidated, never trusted by digest.
    pub terminal: I14CanonicalTerminalResultInputV2<'a>,
    /// Raw watchdog samples retained for diagnosis; order carries no identity.
    /// They may be a bounded subset because the canonical lifecycle separately
    /// binds a semantic-trace digest and independent verification-receipt
    /// identity whose receipt HELM/ledger authenticates.
    pub watchdog_observations: &'a [I14WatchdogObservationV1],
    /// Whether `watchdog_observations` is the complete multi-kind watchdog
    /// observation stream rather than a diagnostic subset. A complete stream
    /// is reconciled with the poll summary and must derive the bound raw root.
    pub watchdog_samples_complete: bool,
    /// Verification-receipt-bound telemetry-only events arriving after final
    /// disposition. They are retained but can never rewrite the canonical
    /// result; receipt authentication is deferred to HELM/ledger.
    pub late_event_tail: I14LateEventTailV2,
    /// Content identity of the clock-calibration artifact.
    pub clock_calibration_artifact: ContentHash,
}

/// Fail-closed noncanonical telemetry-envelope refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum I14TelemetryEnvelopeRefusalV2 {
    /// The authoritative canonical result cannot be constructed.
    CanonicalResult(I14CanonicalResultRefusalV2),
    /// The watchdog sample collection exceeds its admitted bound.
    TooManyWatchdogObservations {
        /// Supplied observation count.
        count: usize,
        /// Maximum admitted observation count.
        cap: usize,
    },
    /// Two watchdog samples reuse one identity.
    DuplicateWatchdogObservationId {
        /// Reused observation identity.
        observation_id: u64,
    },
    /// A sample collection marked complete contradicts the lifecycle summary.
    CompleteWatchdogSampleMismatch,
    /// A boundary's cached latest-poll snapshot is not derivable from the
    /// complete supplied poll stream.
    CompleteWatchdogBoundarySnapshotMismatch {
        /// Boundary whose cached snapshot is inconsistent.
        boundary_ordinal: u64,
        /// Actual latest poll at or before the boundary, if one exists.
        expected_last_poll_monotonic_ns: Option<u64>,
        /// Latest-poll timestamp claimed by the boundary record.
        found_last_poll_monotonic_ns: u64,
    },
    /// A complete sample collection does not derive the bound raw-trace root.
    CompleteWatchdogRawTraceDigestMismatch {
        /// Digest derived from the complete supplied sample trace.
        expected: ContentHash,
        /// Raw trace digest declared by the lifecycle evidence.
        found: ContentHash,
    },
    /// Late-event count and sequence endpoints contradict one another.
    LateEventTailInconsistent,
    /// A purported late-event tail starts no later than terminal selection.
    LateEventTailNotAfterTerminal {
        /// First purported late logical sequence.
        first_logical_sequence: u64,
        /// Selected terminal logical sequence.
        terminal_logical_sequence: u64,
    },
}

fn i14_sorted_watchdog_observations_v2(
    observations: &[I14WatchdogObservationV1],
) -> Result<Vec<I14WatchdogObservationV1>, I14TelemetryEnvelopeRefusalV2> {
    if observations.len() > I14_MAX_WATCHDOG_OBSERVATIONS_V1 {
        return Err(I14TelemetryEnvelopeRefusalV2::TooManyWatchdogObservations {
            count: observations.len(),
            cap: I14_MAX_WATCHDOG_OBSERVATIONS_V1,
        });
    }
    let mut observations = observations.to_vec();
    observations.sort_unstable_by_key(|observation| (observation.observation_id, observation.kind));
    for pair in observations.windows(2) {
        if pair[0].observation_id == pair[1].observation_id {
            return Err(
                I14TelemetryEnvelopeRefusalV2::DuplicateWatchdogObservationId {
                    observation_id: pair[0].observation_id,
                },
            );
        }
    }
    Ok(observations)
}

fn i14_watchdog_raw_trace_digest_from_sorted_v2(
    observations: &[I14WatchdogObservationV1],
) -> ContentHash {
    let mut payload = b"I14_WATCHDOG_RAW_TRACE_V2\0".to_vec();
    i14_push_len(&mut payload, observations.len());
    for observation in observations {
        i14_push_u64(&mut payload, observation.observation_id);
        payload.push(observation.kind as u8);
        i14_push_u64(&mut payload, observation.monotonic_ns);
    }
    hash_domain(WATCHDOG_RAW_TRACE_DOMAIN_V2, &payload)
}

/// Validate and derive the versioned identity of a complete raw watchdog
/// observation trace. Presentation order carries no identity.
pub fn i14_watchdog_raw_trace_digest_v2(
    observations: &[I14WatchdogObservationV1],
) -> Result<ContentHash, I14TelemetryEnvelopeRefusalV2> {
    let observations = i14_sorted_watchdog_observations_v2(observations)?;
    Ok(i14_watchdog_raw_trace_digest_from_sorted_v2(&observations))
}

fn i14_push_raw_request_v2(payload: &mut Vec<u8>, request: I14CancellationRequestV1) {
    i14_push_u64(payload, request.request_id);
    i14_push_u64(payload, request.scope_root);
    i14_push_u64(payload, request.logical_sequence);
    i14_push_u64(payload, request.requested_monotonic_ns);
    i14_push_u64(payload, request.observation_deadline_ns);
    match request.observation {
        None => payload.push(0),
        Some(observation) => {
            payload.push(1);
            i14_push_u64(payload, observation.logical_sequence);
            i14_push_u64(payload, observation.monotonic_ns);
            i14_push_u64(payload, observation.observing_tile_id);
            i14_push_optional_u64(payload, observation.latest_completed_boundary_ordinal);
        }
    }
}

/// Validate and hash all raw timing fields paired with a V2 canonical result.
pub fn i14_telemetry_envelope_digest_v2(
    input: I14TelemetryEnvelopeInputV2<'_>,
) -> Result<ContentHash, I14TelemetryEnvelopeRefusalV2> {
    let watchdogs = i14_sorted_watchdog_observations_v2(input.watchdog_observations)?;
    let canonical = i14_canonical_terminal_result_v2(input.terminal)
        .map_err(I14TelemetryEnvelopeRefusalV2::CanonicalResult)?;
    let terminal_logical_sequence = canonical.local_result().boundary_logical_sequence();
    let tail = input.late_event_tail;
    match (
        tail.event_count,
        tail.first_logical_sequence,
        tail.last_logical_sequence,
    ) {
        (0, None, None) => {}
        (0, _, _) | (_, None, _) | (_, _, None) => {
            return Err(I14TelemetryEnvelopeRefusalV2::LateEventTailInconsistent);
        }
        (count, Some(first), Some(last)) => {
            let represented_slots = last.checked_sub(first).and_then(|span| span.checked_add(1));
            if (count == 1 && first != last) || represented_slots.is_none_or(|slots| slots < count)
            {
                return Err(I14TelemetryEnvelopeRefusalV2::LateEventTailInconsistent);
            }
            if first <= terminal_logical_sequence {
                return Err(
                    I14TelemetryEnvelopeRefusalV2::LateEventTailNotAfterTerminal {
                        first_logical_sequence: first,
                        terminal_logical_sequence,
                    },
                );
            }
        }
    }
    if input.watchdog_samples_complete {
        let mut polls = watchdogs
            .iter()
            .filter(|observation| observation.kind == I14WatchdogObservationKindV1::Poll)
            .copied()
            .collect::<Vec<_>>();
        polls.sort_unstable_by_key(|observation| {
            (observation.monotonic_ns, observation.observation_id)
        });
        let coverage = input.terminal.lifecycle.watchdog_coverage;
        let poll_count_matches = u64::try_from(polls.len()).ok() == Some(coverage.poll_count);
        let endpoints_and_gap_match = match (polls.first(), polls.last()) {
            (Some(first), Some(last)) => {
                let maximum_gap_ns = polls
                    .windows(2)
                    .map(|pair| pair[1].monotonic_ns - pair[0].monotonic_ns)
                    .max()
                    .unwrap_or(0);
                first.monotonic_ns == coverage.first_poll_monotonic_ns
                    && last.monotonic_ns == coverage.last_poll_monotonic_ns
                    && maximum_gap_ns == coverage.maximum_poll_gap_ns
            }
            (None, None) => false,
            _ => unreachable!("first/last presence agrees"),
        };
        if !poll_count_matches || !endpoints_and_gap_match {
            return Err(I14TelemetryEnvelopeRefusalV2::CompleteWatchdogSampleMismatch);
        }
        let mut next_poll_index = 0;
        let mut latest_poll_monotonic_ns = None;
        for record in input.terminal.trace.boundaries {
            while next_poll_index < polls.len()
                && polls[next_poll_index].monotonic_ns <= record.boundary.monotonic_ns
            {
                latest_poll_monotonic_ns = Some(polls[next_poll_index].monotonic_ns);
                next_poll_index += 1;
            }
            if latest_poll_monotonic_ns != Some(record.last_watchdog_poll_monotonic_ns) {
                return Err(
                    I14TelemetryEnvelopeRefusalV2::CompleteWatchdogBoundarySnapshotMismatch {
                        boundary_ordinal: record.boundary.boundary_ordinal,
                        expected_last_poll_monotonic_ns: latest_poll_monotonic_ns,
                        found_last_poll_monotonic_ns: record.last_watchdog_poll_monotonic_ns,
                    },
                );
            }
        }
        let expected = i14_watchdog_raw_trace_digest_from_sorted_v2(&watchdogs);
        let found = coverage.watchdog_raw_trace_digest;
        if found != expected {
            return Err(
                I14TelemetryEnvelopeRefusalV2::CompleteWatchdogRawTraceDigestMismatch {
                    expected,
                    found,
                },
            );
        }
    }
    let trace = input.terminal.trace;
    let lifecycle = input.terminal.lifecycle;
    let mut payload = b"I14_NONCANONICAL_TELEMETRY_ENVELOPE_V2\0".to_vec();
    payload.extend_from_slice(canonical.digest().as_bytes());
    payload.push(u8::from(input.watchdog_samples_complete));
    i14_push_u64(&mut payload, tail.event_count);
    i14_push_optional_u64(&mut payload, tail.first_logical_sequence);
    i14_push_optional_u64(&mut payload, tail.last_logical_sequence);
    payload.extend_from_slice(tail.semantic_trace_digest.as_bytes());
    payload.extend_from_slice(tail.verification_receipt_digest.as_bytes());
    i14_push_u64(&mut payload, trace.campaign_started_monotonic_ns);
    i14_push_len(&mut payload, trace.boundaries.len());
    for record in trace.boundaries {
        i14_push_u64(&mut payload, record.boundary.monotonic_ns);
        i14_push_u64(&mut payload, record.last_watchdog_poll_monotonic_ns);
    }
    let mut ordered_requests = trace.requests.to_vec();
    ordered_requests.sort_unstable_by_key(|request| (request.logical_sequence, request.request_id));
    i14_push_len(&mut payload, ordered_requests.len());
    for request in ordered_requests {
        i14_push_raw_request_v2(&mut payload, request);
    }
    for event in [
        lifecycle.execution_started,
        lifecycle.drain_started,
        lifecycle.drained,
        lifecycle.finalized,
    ] {
        i14_push_u64(&mut payload, event.monotonic_ns);
    }
    // Presence is already bound by the canonical onset-sequence fields, so
    // only present calibrated values are appended. This preserves the frozen
    // encoding of lifecycle records with no failure onset while keeping each
    // supplied raw onset time in telemetry identity.
    if let Some(onset) = lifecycle.first_infrastructure_failure_onset {
        i14_push_u64(&mut payload, onset.event.monotonic_ns);
    }
    if let Some(monotonic_ns) = lifecycle.first_timeout_failure_onset_monotonic_ns {
        i14_push_u64(&mut payload, monotonic_ns);
    }
    match lifecycle
        .spawn_frontier_audit
        .and_then(|audit| audit.last_child_spawn)
    {
        None => payload.push(0),
        Some(event) => {
            payload.push(1);
            i14_push_u64(&mut payload, event.monotonic_ns);
        }
    }
    payload.extend_from_slice(lifecycle.child_lifecycle_raw_trace_digest.as_bytes());
    let coverage = lifecycle.watchdog_coverage;
    i14_push_u64(&mut payload, coverage.poll_count);
    i14_push_u64(&mut payload, coverage.first_poll_monotonic_ns);
    i14_push_u64(&mut payload, coverage.last_poll_monotonic_ns);
    i14_push_u64(&mut payload, coverage.maximum_poll_gap_ns);
    payload.extend_from_slice(coverage.watchdog_raw_trace_digest.as_bytes());
    payload.extend_from_slice(
        lifecycle
            .tile_poll_coverage
            .tile_poll_raw_trace_digest
            .as_bytes(),
    );
    let external = lifecycle.external_heartbeat_coverage;
    payload.extend_from_slice(external.heartbeat_raw_trace_digest.as_bytes());
    i14_push_u64(&mut payload, external.heartbeat_count);
    i14_push_u64(&mut payload, external.maximum_heartbeat_gap_ns);
    i14_push_len(&mut payload, watchdogs.len());
    for observation in watchdogs {
        i14_push_u64(&mut payload, observation.observation_id);
        payload.push(observation.kind as u8);
        i14_push_u64(&mut payload, observation.monotonic_ns);
    }
    payload.extend_from_slice(input.clock_calibration_artifact.as_bytes());
    Ok(hash_domain(TELEMETRY_ENVELOPE_DOMAIN_V2, &payload))
}

/// Artifact categories governed by the I14 retention policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum I14ArtifactCategoryV1 {
    /// Structured event after schema-aware redaction.
    Event = 0,
    /// Frozen/successor manifest after schema-aware redaction.
    Manifest = 1,
    /// Adjudication receipt after schema-aware redaction.
    AdjudicationReceipt = 2,
    /// Structured log after schema-aware redaction.
    Log = 3,
    /// Oracle output after schema-aware redaction.
    OracleOutput = 4,
    /// Replay capsule after schema-aware redaction.
    ReplayCapsule = 5,
    /// Sanitized minimized counterexample.
    MinimizedCounterexample = 6,
    /// Sanitized scientific refutation artifact.
    Refutation = 7,
    /// Sanitized immutable failure bundle.
    FailureBundle = 8,
    /// Raw licensed bytes.
    LicensedBytes = 9,
    /// Raw secret or credential bytes.
    SecretBytes = 10,
    /// Raw specimen identity or re-identification material.
    SpecimenIdentity = 11,
    /// Raw governed holdout/validation/population bytes before or after controlled reveal.
    GovernedHoldoutBytes = 12,
    /// Raw derived slice whose content remains sensitive.
    DerivedSensitiveSlice = 13,
    /// Raw unredacted diagnostic slice.
    UnredactedDiagnosticSlice = 14,
}

impl I14ArtifactCategoryV1 {
    /// Canonical order used by retention conformance tests.
    pub const ALL: [Self; 15] = [
        Self::Event,
        Self::Manifest,
        Self::AdjudicationReceipt,
        Self::Log,
        Self::OracleOutput,
        Self::ReplayCapsule,
        Self::MinimizedCounterexample,
        Self::Refutation,
        Self::FailureBundle,
        Self::LicensedBytes,
        Self::SecretBytes,
        Self::SpecimenIdentity,
        Self::GovernedHoldoutBytes,
        Self::DerivedSensitiveSlice,
        Self::UnredactedDiagnosticSlice,
    ];
}

/// Durable retention class assigned to an I14 artifact category.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum I14RetentionClassV1 {
    /// Sanitized manifest/evidence/replay material retained durably.
    EvidenceDurable,
    /// Sanitized counterexamples/refutations/failures retained permanently.
    FailurePermanent,
    /// Raw sensitive material retained only under encrypted governed capability.
    GovernedRestricted,
}

/// Complete typed retention rule for one artifact category.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(clippy::struct_excessive_bools)]
pub struct I14RetentionRuleV1 {
    /// Assigned durable retention class.
    pub class: I14RetentionClassV1,
    /// Schema-aware sanitization is mandatory before retention.
    pub sanitize_before_retention: bool,
    /// Storage must be encrypted and capability controlled.
    pub encrypted_capability_controlled: bool,
    /// Every access must enter the complete access ledger.
    pub complete_access_ledger: bool,
    /// The class-specific retention/erasure decision must be retained.
    pub retain_retention_or_erasure_decision: bool,
}

/// Return the fail-closed retention rule for an I14 artifact category.
#[must_use]
pub const fn i14_retention_rule_v1(category: I14ArtifactCategoryV1) -> I14RetentionRuleV1 {
    let (class, sanitize_before_retention, encrypted_capability_controlled) = match category {
        I14ArtifactCategoryV1::Event
        | I14ArtifactCategoryV1::Manifest
        | I14ArtifactCategoryV1::AdjudicationReceipt
        | I14ArtifactCategoryV1::Log
        | I14ArtifactCategoryV1::OracleOutput
        | I14ArtifactCategoryV1::ReplayCapsule => {
            (I14RetentionClassV1::EvidenceDurable, true, false)
        }
        I14ArtifactCategoryV1::MinimizedCounterexample
        | I14ArtifactCategoryV1::Refutation
        | I14ArtifactCategoryV1::FailureBundle => {
            (I14RetentionClassV1::FailurePermanent, true, false)
        }
        I14ArtifactCategoryV1::LicensedBytes
        | I14ArtifactCategoryV1::SecretBytes
        | I14ArtifactCategoryV1::SpecimenIdentity
        | I14ArtifactCategoryV1::GovernedHoldoutBytes
        | I14ArtifactCategoryV1::DerivedSensitiveSlice
        | I14ArtifactCategoryV1::UnredactedDiagnosticSlice => {
            (I14RetentionClassV1::GovernedRestricted, false, true)
        }
    };
    I14RetentionRuleV1 {
        class,
        sanitize_before_retention,
        encrypted_capability_controlled,
        complete_access_ledger: true,
        retain_retention_or_erasure_decision: true,
    }
}

/// Compute the canonical exhaustive version-1 I14 terminal-status table digest.
#[must_use]
pub fn i14_terminal_status_table_digest_v1() -> ContentHash {
    let mut payload = b"I14_TERMINAL_STATUS_TRUTH_TABLE_V1\0".to_vec();
    payload.extend_from_slice(&(I14_TERMINAL_STATUS_TABLE_V1_TUPLES as u32).to_le_bytes());
    let mut tuples = 0usize;
    for execution in I14ExecutionDisposition::ALL {
        for claim in I14ClaimAdjudication::ALL {
            for completeness in I14EvidenceCompleteness::ALL {
                for integrity in I14EvidenceIntegrity::ALL {
                    for input in I14InputValidity::ALL {
                        for domain in I14DomainApplicability::ALL {
                            for support in I14OperationalSupport::ALL {
                                for receipt in I14ReceiptValidity::ALL {
                                    let raw = I14TerminalStatusV1 {
                                        execution,
                                        claim,
                                        completeness,
                                        integrity,
                                        input,
                                        domain,
                                        support,
                                        receipt,
                                    };
                                    let evaluated = i14_evaluate_terminal_status_v1(raw);
                                    payload.extend_from_slice(&raw.tags());
                                    payload.extend_from_slice(&evaluated.normalized.tags());
                                    payload.push(u8::from(
                                        evaluated.normalization.domain_forced_indeterminate,
                                    ));
                                    payload.push(u8::from(
                                        evaluated.normalization.receipt_marked_malformed,
                                    ));
                                    payload.push(evaluated.exit_code);
                                    tuples += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    debug_assert_eq!(tuples, I14_TERMINAL_STATUS_TABLE_V1_TUPLES);
    hash_domain(TERMINAL_STATUS_TABLE_DOMAIN_V1, &payload)
}

/// Build the I14 multirung EMC/harness draft. Consumers freeze it themselves;
/// freezing preregisters authority but does not prove that any named solver,
/// checker, benchmark, script, theorem, or laboratory artifact exists.
#[must_use]
pub fn i14_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I14",
        title: "Multirung EMC/harness gate: native topology, MTL/RLGC, PEEC, grounding, shielding, bearing currents, full-wave escalation, source-to-victim evidence, robust design, and passive-causal composition theorems",
        version: 1,
        explicits: FiveExplicits {
            units: "six-base SI quantity system L M T I Theta N; conductor and shield geometry in metres; per-unit-length R,L,G,C in ohm-per-metre, henry-per-metre, siemens-per-metre, and farad-per-metre. RMS phasors use Re{alpha X exp(+i omega t)}, where alpha is the canonical positive algebraic root alpha^2=2 encoded by minimal polynomial x^2-2 plus a rational isolating interval, and complex power is V I* or integral(E cross H*) dot n; radiation power is an outward-oriented boundary flux rather than an interior loss, while peak-phasor adapters apply the exact symbolic one-half factor. Numeric adapters evaluate algebraic constants only through outward-rounded enclosures and never promise bit-exact multiply-then-divide recovery. Laplace transforms use X(s)=integral_0^infinity x(t) exp(-s t) dt. Strict/asymptotically stable rational realizations have poles in Re(s)<0; generalized positive-real passive immittances are analytic on Re(s)>0 and may have only admitted simple boundary-axis poles (including zero) with positive-semidefinite Hermitian residue matrices plus separately checked behavior at infinity. Frequency response is the non-tangential Re(s) downarrow 0 boundary value where it exists, not an unconditional substitution at a pole. Every terminal, field, spectrum, uncertainty, normalized error, cost, mass, temperature, and safety QoI declares its unit and a pre-candidate AcceptanceCard with exact comparator arithmetic",
            seeds: "Philox 4x32-10 deterministic fixture streams use aliases i14/<fixture-id>/<purpose>; d=BLAKE3::derive_key('org.frankensim.i14.fixture-stream.v1', alias_utf8); k0=LE32(d[0..4]); k1=LE32(d[4..8]); c0=low32(semantic_case_index); c1=high32(semantic_case_index); c2=low32(output_block_ordinal); c3=high32(output_block_ordinal). Output lanes serialize exactly as LE32(r0)||LE32(r1)||LE32(r2)||LE32(r3); lane selection is the Philox output-word index and is never folded into the counter, and native-endian casts are forbidden. Development indices 0..=4095, PublicReplayCore 65536..=69631, PublicReplayMax 131072..=135167, and falsifier indices 196608..=212991 are disjoint inclusive ranges. PublicReplay deterministic holdouts provide replay/conformance mechanics only and no IID or untouched-data authority. GovernedBlindSlot and GovernedPopulationSlot authority require their separately scoped candidate-before-protected-access slot/envelope/discharge transaction in i14-campaign-policy-v1. GovernedPhysicalSlot permits logged least-privilege AuthorizedCalibration and AsBuiltModelInstantiation access only to frozen calibration/model-input strata, then requires CandidateFrozen before candidate-side validation access or commitment opening and the joined-root/envelope/discharge transaction before adjudication. GovernedStandards alone permits least-privilege licensed-input AuthorizedConstruction after GovernanceCommitted and before CandidateFrozen; it remains NoPromotionAuthority until its post-construction, pre-adjudication same-ID envelope/discharge transaction commits",
            budgets: "smoke <= 120 s and <= 8 GiB on one host; core <= 90 min and <= 32 GiB; max <= 18 h and <= 96 GiB on a quiet campaign host; theorem/falsifier max <= 24 h and <= 96 GiB. Core total resource ceilings are 4096 graph/trace records, 16384 field/quadrature unknowns, 1024 search/tree nodes, or 256 formal declarations; Core logical poll tiles contain at most 64 graph/trace records, 256 field/quadrature unknowns, 16 search/tree nodes, or 4 formal declarations, while an asynchronous cancellation watchdog polls at intervals <=25 ms without ending or repartitioning a logical tile. Core measured request-to-observation <=250 ms, drain-trigger-to-drained <=2 s, and drained-to-finalized <=2 s on the admitted host. Max total resource ceilings are four times the Core resource counts; Max logical poll tiles contain at most four times the Core item counts, while the asynchronous watchdog polls at intervals <=100 ms without changing logical tile membership or order. Max request-to-observation <=1 s, drain-trigger-to-drained <=8 s, and drained-to-finalized <=8 s. The drain-trigger clock starts at the referenced on-time observation, first nanosecond after the missed inclusive request deadline, structurally admitted receipt-bound infrastructure-failure onset whose receipt HELM/ledger authenticates, or drain-start timestamp for CancellationObserved, ObservationTimeoutDrain, InfrastructureFailure, or NonCancellationDrain respectively. Overall campaign wall-budget expiry is TimedOut. External solvers/checkers must expose interrupt, bounded terminate-drain, authenticated no-partial-publication and a heartbeat within the tier watchdog quantum; admission refuses an indivisible operation or logical tile lacking a demonstrated tier response bound. p50/p95/p99 time and peak memory are unclaimed until measured against a pinned machine fingerprint; numerical and decision accuracy are fixed only by per-claim directed tolerances and fixture-local budgets",
            versions: "fs-vmanifest schema v2; HarnessGraph schema v1; PortSchema electrical RMS-phasor convention v1; MtlOperator/RLGC schema v1; PEEC quotient-and-gauge schema v1; FidelityRegion/Crosswalk schema v1; FullWaveProblem/Receipt schema v1; source-probe-victim ModeLedger schema v1; I14 TerminalStatusTruthTable schema v1 with code-pinned domain-separated digest; I14 legacy local scoped terminal-cause selector schema v1; I14 legacy single-boundary canonical terminal-result schema v1 under org.frankensim.i14.terminal-result.v1; I14 legacy noncanonical telemetry-envelope schema v1 under org.frankensim.i14.telemetry-envelope.v1; I14 validated cancellation-card schema v2 under org.frankensim.i14.cancellation-card.v2; I14 recomputed clock-free logical execution trace schema v2 under org.frankensim.i14.logical-execution-trace.v2; I14 authoritative genesis-to-first-terminal trace schema v2 under org.frankensim.i14.terminal-trace.v2; I14 authoritative canonical terminal-result schema v2 under org.frankensim.i14.terminal-result.v2; I14 receipt-bound noncanonical telemetry-envelope schema v2 under org.frankensim.i14.telemetry-envelope.v2; I14 complete raw-watchdog trace schema v2 under org.frankensim.i14.watchdog-raw-trace.v2; I14 typed four-variant drain-trigger/lifecycle schema v2; I14 typed artifact-retention policy schema v1; uncertainty and LossOwnershipReceipt schema v1; theorem-card target schema v1; toolchain pinned by rust-toolchain.toml and sibling revisions by constellation.lock",
            capabilities: "core: harness-graph, synthetic-ap242-adapter-mechanics, rlgc-operator, mtl-propagation, peec, ground-bond-shield, bearing-current-path, source-probe, immunity-victim, fixed-rung-router, fixed-regime-adjoint, emc-uq-inference, fullwave-schema; maximal feature gates: fullwave-feec, maxwell-bem, maxwell-fmm, adaptive-fidelity, robust-emc-design, safety-case-link, passive-causal-composition-theorem, hypercohomology-obstruction-theorem, cover-refinement-naturality-theorem, kyp-sheaf-bridge-theorem, governed-standards-crosswalk, governed-laboratory-validation, governed-emc-reliability, governed-bearing-population; no network or production FFI; deterministic mode mandatory for G5; external licensed standards, laboratory data, production populations, proprietary as-built geometry, and blind custody enter only through separately scoped fs-vvreg-governed artifacts",
        },
        claims: i14_claims(),
        fixtures: i14_fixtures(),
        obligations: i14_obligations(),
        waivers: i14_waivers(),
        amendment_rules: "After campaign start every semantic change forks a successor through FrozenManifest::amend. Changes to topology or occurrence identity, EmConventionCard, port/reference/gauge/phasor/Laplace convention, material state, loss ownership, fidelity region, crosswalk, AcceptanceCard arithmetic, fixture partition or realized holdout root, tolerance, oracle/TCB, theorem proposition, standard edition, budget, capability, custody stage, or discharge receipt invalidate exactly the affected predecessor claim and producer-leaf authorities. No result may edit this version in place, select a favorable component set, widen a band after evidence, skip a custody stage, or resolve an AuthorityContradiction without an independently audited successor",
    }
}

#[allow(clippy::too_many_lines)]
fn i14_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i14-harnessgraph-identity-connectivity",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "A native HarnessGraph gives conductors, shield layers, connector occurrences and pins, splices, bundles, routes, frames, terminals, nets, chassis/reference conductors, material-card links, and as-built lineage stable typed identities. Electrical connectivity follows only explicit incidence and pin/splice semantics; geometric proximity never creates a net. Ambiguous, missing, duplicated, or orientation-inconsistent native connectivity refuses",
            hypotheses: &[
                "every entity has distinct lineage EntityId, exact-source identity where imported, representation identity, and semantic occurrence identity; repeated connector definitions do not alias occurrences",
                "route curves carry declared frames, orientation, length/error bounds, bundle membership intervals, conductor order, shield/drain termination, and topology independent of electrical property values",
                "connector pin, splice, terminal, net, ground/chassis, open-circuit, and shield-continuity relations are explicit graph edges; neither geometric contact nor name similarity is an electrical bond",
            ],
            qoi: "native_harness_identity_connectivity_and_orientation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i14.oracle.harnessgraph.v1 at fs-vmanifest-oracles/i14/harnessgraph.rs::reconstruct_and_compare_incidence",
                independent: true,
                tcb_overlap: "shares canonical native source bytes and identity hash primitive only; independently reconstructs incidence, occurrences and route transforms",
            },
            activation: "the native HarnessGraph schema is frozen before production implementation",
            kill: "one proximity-created net, swapped or aliased occurrence pin, lost open/shield-drain edge, nondeterministic graph identity, or silently ambiguous native edge kills this claim",
            fallback: "retain the native partial graph with structured Unknown edges and ranked repair",
            no_claim: "native topology identity does not determine adapter correctness, RLGC, current, shielding, EMC, manufacturing correctness, or licensed AP242 conformance",
        },
        ClaimSpec {
            id: "i14-synthetic-ap242-adapter-mechanics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For one legally authored AP242-shaped supported subset, a one-way synthetic adapter binds exact source bytes, units, frames, assembly and occurrence reuse, connector/pin identity, open/short/splice/shield/drain/bond semantics, route and material lineage, and every information loss to a native HarnessGraph result. Unsupported, ambiguous, duplicated, unresolved, or orientation-inconsistent records quarantine or refuse rather than disappear",
            hypotheses: &[
                "the authored subset grammar, source-byte root, edition-neutral record semantics, unit/frame transforms, occurrence map and loss taxonomy freeze before adapter execution",
                "every source construct maps to exactly one native identity, explicit loss record, quarantined artifact, or structured refusal, with repeated definitions kept distinct from semantic occurrences",
                "the independent oracle reconstructs source occurrence and transform graphs without production importer code; neither deck contains licensed standards text",
            ],
            qoi: "synthetic_ap242_subset_adapter_identity_transform_loss_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i14.oracle.ap242_adapter.v1 at fs-vmanifest-oracles/i14/ap242_adapter.rs::reconstruct_occurrences_and_losses",
                independent: true,
                tcb_overlap: "shares authored source bytes and quantity dimensions only; independently reconstructs occurrence, transform, connectivity and loss edges",
            },
            activation: "the synthetic adapter subset and native HarnessGraph target schema are frozen",
            kill: "one aliased occurrence, wrong unit/frame transform, invented connectivity, omitted supported record, unaccounted loss, or favorable treatment of ambiguity kills the adapter claim",
            fallback: "retain the authenticated source and any independently valid native partial graph while quarantining unsupported or ambiguous records",
            no_claim: "synthetic subset mechanics are not licensed-edition AP242 conformance, a universal STEP importer, physical harness validation, or EMC evidence",
        },
        ClaimSpec {
            id: "i14-fullwave-problem-convention-admission",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "A FullWaveProblem/Receipt contract admits no implicit electromagnetic convention. A structured EmConventionCard binds time versus frequency domain, RMS exp(+i omega t) phasors, symbolic algebraic RMS/peak scaling, real-block complex encoding, primal/twisted-dual FEEC field degrees, orientation, source amplitudes and relative phases, port normal, terminal reference, material frequency/temperature state, constitutive/passivity signs, outgoing-wave branch, boundary/exterior category, topology, gauge, source/storage/dissipation/radiation power ownership, QoI, adjoint convention, and budgets. Only a global reference-phase shift is presentation-equivalent; physical relative phases and power-owner categories are never normalized away. A versioned crosswalk transforms the entire coupled convention card or refuses, with no field-by-field partial conjugation",
            hypotheses: &[
                "frequency-domain Maxwell uses curl E=-i omega B and curl H=J_impressed+(sigma+i omega epsilon)E under the frozen exp(+i omega t) convention; other conventions enter only through a declared conjugation/sign transform",
                "RMS phasors define average complex terminal power V I* and outward complex Poynting flux integral(E cross H*) dot n; peak adapters use the symbolic algebraic alpha^2=2 and exact one-half power factor, with outward-rounded numeric evaluation rather than false bit-exact amplitude round trips; port current is positive into the modeled subsystem",
                "for exp(+i omega t), the passive outgoing scalar Green branch is exp(-i k r)/(4 pi r) with a frozen passive wavenumber branch; PML, BEM/FMM kernels, impedance conditions, radiation traces and adjoints use the same branch and orientation",
                "converting to the conjugate phasor convention simultaneously conjugates sources, relative phases, constitutive functions, field/port traces, impedances/admittances, outgoing kernels, PML stretches, power forms and adjoint operators; a global phase quotient never identifies different relative source phases",
                "the complete conductivity/dielectric/magnetic/bianisotropic constitutive or admittance block freezes the Hermitian dissipative inequality appropriate to exp(+i omega t); conductivity, dielectric loss, radiation flux, PML absorption, source work and storage are distinct owners and cannot be double counted or generically relabeled loss",
                "field degrees, primal/dual orientation, relative boundary complex, harmonic sector, gauge/divergence treatment, material passivity class, source compatibility, and PEC/PMC/impedance/port/PML/exterior-BEM boundary ownership are serialized",
                "admission states are orthogonal to numerical, validation, theorem, standards, and safety evidence and cannot mint those authorities",
            ],
            qoi: "fullwave_semantic_admission_and_explicit_crosswalk_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i14.oracle.fullwave_schema.v1 at fs-vmanifest-oracles/i14/fullwave_schema.rs::check_conventions",
                independent: true,
                tcb_overlap: "shares quantity dimensions and canonical bytes; independently applies sign, conjugation, orientation, and RMS/peak transformations",
            },
            activation: "before any full-wave solver, boundary, BEM, or FMM result is inspected",
            kill: "one invented default, wrong factor-of-two, conjugation/sign ambiguity, incompatible field degree, missing boundary owner, or passive-path admission of active/noncausal media kills the contract",
            fallback: "return a structured semantic diff or Refused problem; a lower admitted rung remains available only under its own contract",
            no_claim: "well-typed conventions prove neither well-posedness, absence of spurious modes, solver accuracy, open-boundary quality, nor physical validation",
        },
        ClaimSpec {
            id: "i14-mtl-rlgc-operator-admission",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A versioned MtlOperator binds ordered conductors and explicit reference/return topology to dimensioned frequency- and temperature-dependent per-unit-length R,L,G,C data and to the frozen series/shunt operators Z(s) and Y(s), whose frequency-boundary values are Z(i omega)=R+i omega L and Y(i omega)=G+i omega C only for the admitted decomposition. Source evidence, covariance or certified bounds, reciprocity class, and a QuasiTemValidityReceipt are inseparable. Potential coefficients and Maxwell capacitance are related only on their declared gauge quotient. Every accepted passive class has compatible realness/conjugate symmetry, terminal continuity, explicit conductor/reference transforms, and either full analytic positive-real Z(s)/Y(s) authority or an explicitly weaker band dissipativity receipt",
            hypotheses: &[
                "each matrix sample binds conductor order, reference conductor or quotient, terminal orientation, geometry/material lineage, frequency, temperature, phasor convention, and measured/computed/fitted evidence class",
                "for a quasi-static or nondispersive decomposition, R and G dissipative Hermitian parts are PSD and L and C storage blocks are positive on the admitted physical quotient under their declared reciprocal/passive hypotheses; this statement is not transferred blindly to dispersive fitted matrices",
                "dispersive authority belongs to the full matrix-valued analytic positive-real Z(s)/Y(s) on Re(s)>0 or to a passive internal-state/convolution realization: F(conj(s))=conj(F(s)), Hermitian(F(s)) is PSD in the open right half-plane, quotient regularity is explicit, and any admitted boundary-axis pole is simple with PSD Hermitian residue plus checked zero/infinity behavior. Strict positive-real/asymptotically stable subclasses keep all poles in Re(s)<0. Storage uses realization states; frequency-derivative energy formulas require their separately proved lossless or weak-loss regime. Sampled nonnegative Hermitian parts on a bounded omega set provide only band-limited dissipativity evidence",
                "geometry-derived entries consume their own I03/PEEC/MQS extraction receipts; measured data retain calibration/covariance, and interpolation cannot upgrade evidence or cross a material domain",
                "quasi-TEM validity names transverse/longitudinal scale, higher-mode cutoff margin, discontinuity severity, skin/proximity and material-dispersion indicators, source spectrum, and requested QoI",
            ],
            qoi: "maximum_preregistered_normalized_rlgc_structure_source_and_validity_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.rlgc.v1 at fs-vmanifest-oracles/i14/rlgc.rs::check_operator_and_quotient",
                independent: true,
                tcb_overlap: "shares canonical geometry/material/source records; analytic coax/two-wire and independent matrix/quotient checks do not call production extraction or interpolation",
            },
            activation: "HarnessGraph identity and at least one admitted extraction or measurement source are green",
            kill: "missing reference/gauge, incompatible terminal order, unit or quotient mismatch, negative passive loss, unsupported interpolation, evidence laundering, or invalid quasi-TEM admission kills the operator claim",
            fallback: "use admitted tabulated data with weaker evidence, route to PEEC/full wave, or return Unknown with the violated metric",
            no_claim: "an RLGC table or sampled j-omega dissipativity check is not a causal positive-real transient theorem, a universal common/differential basis, connector field solution, higher-mode solution, or full-wave validation",
        },
        ClaimSpec {
            id: "i14-mtl-passive-causal-propagation",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Within a pinned quasi-TEM region and the RMS exp(+i omega t) convention, frequency-domain propagation satisfies dV/dz=-(R+i omega L)I and dI/dz=-(G+i omega C)V for the admitted decomposition. Frequency- and time-domain multiconductor telegrapher solutions over admitted causal operators reproduce propagation, attenuation, reflection, crosstalk, and mode conversion while closing terminal power, stored energy, transmitted flux, and conductor/dielectric loss exactly once. Stable positive-real rational models are used only where their analytic domain and asymptotics are valid; diffusive, delayed, convolutional, fractional, or tabulated causal classes remain distinct. Repeated modes carry deterministic invariant-subspace identities rather than invented individual labels",
            hypotheses: &[
                "with +z conductor current, dV/dz=-(R+i omega L)I and dI/dz=-(G+i omega C)V under the RMS exp(+i omega t) convention; endpoint port currents are separately oriented positive into each connected subsystem",
                "the fit proves or encloses its exact strict or generalized positive-real class, causality/analytic-domain or class-specific supply inequality, real time response, open-half-plane poles or admitted simple boundary poles with PSD residues, zero/high-frequency asymptotics, conditioning, truncation, and any passivity-repair perturbation",
                "terminations, connectors, descriptor circuits, source bandwidth, convolution/state initialization, timestep, multirate boundary, and reference conductor are admitted and have unique loss/storage owners",
                "clustered/crossing modes use power-normalized invariant subspaces with deterministic transport; common/differential coordinates are application-declared and invertible, never universally canonical",
            ],
            qoi: "maximum_preregistered_normalized_mtl_wave_port_power_and_time_frequency_discrepancy",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.mtl_wave.v1 at fs-vmanifest-oracles/i14/mtl_wave.rs::solve_analytic_and_balance",
                independent: true,
                tcb_overlap: "shares frozen line/termination/source bytes; analytic matrix-exponential and independent convolution/dense descriptor routes do not call production fitter or propagator",
            },
            activation: "the MtlOperator and passive-causal model class are admitted",
            kill: "unstable/noncausal pole or delay, active pocket hidden between samples, mode-label discontinuity, wrong reflection sign/delay, double-counted loss, nonclosing power, or silent higher-mode/discontinuity extrapolation kills the lane",
            fallback: "retain the raw frequency operator, choose a stronger causal operator class, escalate to PEEC/full wave, or return Unknown",
            no_claim: "small fit residual or sampled j-omega positivity is not a global passivity/causality theorem; compact connector and source models are not 3-D field/device solutions",
        },
        ClaimSpec {
            id: "i14-peec-extraction-power-mor",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Quasistatic and explicitly admitted retarded PEEC extraction keeps partial inductance, potential coefficient P, quotient capacitance C, resistance, skin/proximity operators, charge continuity, Lorenz-gauge/reference data, and retardation errors distinct. Circuit stamping closes discrete charge and port-power identities, separates storage, conductor/dielectric dissipation, source work, and exterior radiation, and assigns positive-real/causal authority only to an analytic model class. Structure-preserving reduction retains declared nullspaces, ports, reciprocity, passivity/causality status, and held-out transfer/time error or refuses",
            hypotheses: &[
                "canonical filament/surface/volume unknowns, conductor components, incidence/continuity maps, terminal quotient, gauge, orientation, singular self-term treatment, geometry/discretization and material class are frozen",
                "P-to-C inversion or pseudoinversion occurs only on the declared reference quotient with explicit nullspace; reciprocal quasistatic symmetry/PSD statements do not cross into nonreciprocal, active, or unsupported dispersive media",
                "retarded scalar and vector potential/delay terms are co-discretized with discrete continuity and Lorenz-gauge residuals, causal kernels, initial history, delay/quadrature bounds, and wavelength admission",
                "MNA stamping uses RMS effort/flow orientation, exact algebraic power where claimed, named numerical/model defects, analytic-domain passivity or supply witnesses, and radiation only as boundary Poynting flux",
                "MOR binds full/reduced hashes, projection maps, descriptor constraints, port/tangent subspaces, validity band, error indicator, and an independent full-order/dense held-out adjudication",
            ],
            qoi: "maximum_preregistered_normalized_peec_charge_gauge_power_passivity_and_reduction_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.peec.v1 at fs-vmanifest-oracles/i14/peec.rs::extract_stamp_and_adjudicate",
                independent: true,
                tcb_overlap: "shares geometry/material and terminal bytes; analytic loop/plate/coax plus dense independently assembled potential/inductance and descriptor checks do not share production quadrature, compression, or MOR",
            },
            activation: "HarnessGraph identity and PEEC quotient/gauge schema are green",
            kill: "P/C conflation, missing nullspace/reference, scalar/vector retardation inconsistency, sampled-only passivity promotion, radiation-as-local-loss, active or noncausal reduced pocket, missing descriptor mode, or held-out error above one kills the affected PEEC authority",
            fallback: "use the admitted full-order dense model, narrower quasistatic path, MTL/full-wave rung, or Unknown with the failed identity named",
            no_claim: "a passive stamp does not validate extraction geometry; good port fit does not prove internal fields; dense agreement does not establish physical or regulatory validity",
        },
        ClaimSpec {
            id: "i14-ground-bond-shield-current-closure",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Grounding, bonding, reference-conductor, shield, aperture/seam, pigtail and 360-degree termination, transfer-impedance/admittance, and common-impedance models form one oriented current/power graph with explicit frequency, construction, contact, temperature, environment, and rung validity. Every source and return current closes, every conductor/dielectric/shield/contact loss and every radiation boundary flux has exactly one owner, and geometric touch never substitutes for a declared bond",
            hypotheses: &[
                "chassis/reference semantics, bond/joint topology, contact material/state, conductor/shield orientation, drain/termination geometry, aperture/seam model, current direction, gauge, and validity domain are explicit",
                "transfer impedance/admittance and distributed shield loss bind construction, frequency, temperature, mechanical/environment state, calibration/covariance or certified bounds, passivity/causality class, and termination geometry",
                "MTL, PEEC, circuit, and full-wave adapters preserve terminal/reference identities and RMS port power; compatibility and conservation are independently checked and model agreement does not excuse a shared wrong net graph",
            ],
            qoi: "maximum_preregistered_normalized_return_current_shield_transfer_and_loss_ownership_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.ground_shield.v1 at fs-vmanifest-oracles/i14/ground_shield.rs::reconcile_paths",
                independent: true,
                tcb_overlap: "shares frozen HarnessGraph and property records; independently assembled net/field balance, analytic coax/ground-loop, and retained measurement routes do not share production path reduction",
            },
            activation: "HarnessGraph plus at least one admitted MTL/PEEC/circuit field path is green",
            kill: "hidden/missing bond, floating or swapped reference, sign-wrong return, unowned/double-counted loss, noncausal fit, silent pigtail/aperture extrapolation, or normalized closure above one kills the path claim",
            fallback: "return Unknown for affected path and expose the open return/loss owner; escalate shield/aperture region to full wave",
            no_claim: "no universal shielding coefficient, construction-independent transfer impedance, legal EMC compliance, or physical validation outside the named specimen and state",
        },
        ClaimSpec {
            id: "i14-bearing-current-hybrid-path",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A motor-drive common-mode path couples converter source, winding/stator/rotor/shaft capacitances, frame/brush grounds, bearings, lubricant-film state, and return network with explicit voltage/current orientation and hybrid breakdown/recovery modes. Current, stored energy, dissipated event work, and thermal credits close under the declared bearing-film card; event ordering, hysteresis, speed/load/temperature dependence, uncertainty, and held-out measurements remain visible",
            hypotheses: &[
                "machine geometry and capacitance/impedance source, shaft/frame/bearing orientation, lubricant/material/history state, speed, load, temperature, brush/ground topology, controller clock, and common-mode source spectrum are frozen",
                "film conduction, displacement, breakdown, arc/discharge, recovery, and any erosion proxy use versioned hybrid constitutive/event cards with explicit guards, reset ordering, storage/dissipation and calibration domains",
                "current-path evidence is distinct from tribological damage or service-life evidence; simultaneous contacts, chatter, unresolved event roots, or missing telemetry return set-valued/Unknown rather than an arbitrary path",
            ],
            qoi: "maximum_preregistered_normalized_shaft_voltage_bearing_current_event_and_energy_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.bearing_current.v1 at fs-vmanifest-oracles/i14/bearing_current.rs::check_hybrid_path",
                independent: true,
                tcb_overlap: "shares machine/path and material-card bytes; analytic capacitive dividers, independent hybrid circuit integration, and governed measurements do not share production event solver",
            },
            activation: "ground/shield path plus source and bearing-film cards are admitted",
            kill: "missing path owner, wrong bearing-voltage ratio/current sign, unclosed event energy, hidden retry/order dependence, Supported result after unresolved event, or held-out score above one kills the path claim",
            fallback: "linear bounded capacitive path with breakdown/damage Unknown, or full refusal when the return graph is incomplete",
            no_claim: "bearing-current simulation alone is not a fluting, erosion, lubricant-life, bearing-life, product-safety, or physical-validation theorem",
        },
        ClaimSpec {
            id: "i14-core-fidelity-crosswalk-routing",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A deterministic core router compares and composes EQS, MTL, and PEEC only through content-addressed FidelityRegion and Crosswalk artifacts carrying topology, reference conductors, ports, RMS work variables, state projections, modeled/omitted loss, validity metrics, discrepancy/error bounds, proof state, cost, hysteresis, and fallback. Adjacent-rung power and charge close, evidence color never strengthens by routing, unknown eligibility returns Unknown or an admitted conservative rung, and a fixed-rung baseline remains replayable",
            hypotheses: &[
                "electrical size, retardation, transverse/higher-mode scale, skin/proximity, discontinuity, material dispersion/thermal state, source spectrum, topology, QoI, reference conductor, and error/cost budgets are explicit",
                "one-way and bidirectional crosswalks name exact state/port maps, units, energy/coenergy sign, loss ownership, approximation theorem or empirical discrepancy, naturality/power defect, and TCB overlap",
                "routing uses deterministic total order and route hysteresis under identical evidence/budget; unresolved overlap, nonpassive coupling, stale calibration, or missing loss owner cannot choose the cheapest rung",
                "the core authority ends at fixed/admitted EQS-MTL-PEEC routing; adaptive full-wave composition and global cross-rung error theorems remain separate maximal claims",
            ],
            qoi: "core_route_admit_escalate_unknown_and_crosswalk_power_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.core_router.v1 at fs-vmanifest-oracles/i14/core_router.rs::adjudicate_crosswalks",
                independent: true,
                tcb_overlap: "shares exact rung outputs and crosswalk declarations; independently maps ports, recomputes charge/power/loss and compares held-out adjacent-rung QoIs without production routing logic",
            },
            activation: "EQS, MTL, and PEEC terminal receipts needed by the requested QoI are available",
            kill: "reference/port mismatch, hidden loss, nonpassive or noncausal composition, holdout-bound failure, evidence strengthening, route chatter, cost-model drift outside its band, or silent Unknown eligibility kills the route",
            fallback: "fixed-rung execution at the highest admitted conservative core rung or structured Unknown/refusal",
            no_claim: "cross-rung agreement is numerical discrepancy evidence, not physical validation, theorem proof, full-wave adequacy, or standards compliance",
        },
        ClaimSpec {
            id: "i14-switching-source-probe-semantics",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Switching sources and conducted/radiated probes are versioned measurement operators, not generic spectra: device/common/differential source waveforms and impedances, clocks/phases/events, reference conductors, LISN/current/voltage/field probes, loading, orientation, bandwidth, transfer function, window, detector and near/far-field semantics map into admitted rungs with exact lineage. Time traces reconstruct their declared spectra and injected/source/field/circuit power closes without hidden sources or probe energy",
            hypotheses: &[
                "source waveform, edge model, impedance, clock/phase, common/differential decomposition, event identity, converter/motor state, and aliasing/bandwidth limits are explicit",
                "probe calibration, loading, orientation/frame, transfer function, bandwidth/window/detector, RMS/peak and near/far-field convention, numerical uncertainty, and valid spatial/frequency region are frozen",
                "source, coupling, probe, material, measurement, and model-rung uncertainties are separate owners; a source/probe adapter cannot absorb victim or standards authority",
            ],
            qoi: "maximum_preregistered_normalized_source_spectrum_probe_and_power_chain_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.source_probe.v1 at fs-vmanifest-oracles/i14/source_probe.rs::reconstruct_and_balance",
                independent: true,
                tcb_overlap: "shares raw waveform and calibration bytes; independent exact DFT/window/detector, circuit energy, and analytic probe routes do not share production source/probe compiler",
            },
            activation: "a source card, probe card, and admitted coupling rung exist",
            kill: "wrong phase/reference/orientation, aliasing or detector mismatch, unmodeled probe loading, spectrum reconstruction or power score above one, hidden source, or silent out-of-band extrapolation kills the scoped chain",
            fallback: "retain raw waveform/field with probe result Unknown and ranked calibration/bandwidth/rung remediation",
            no_claim: "source/probe agreement is scoped measurement semantics, not a universal hardware emissions model, product compliance, laboratory validation, or causal victim attribution",
        },
        ClaimSpec {
            id: "i14-immunity-victim-mode-ledger",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Conducted, radiated, BCI, EFT/surge, and explicitly scoped ESD-like injection cards couple through oriented ports to a versioned victim state machine whose functional-performance class, upset/latch/reset/recovery guards, hysteresis, superdense simultaneous-event order, telemetry, and safety-aware abort/finalize policy are frozen. An independent source-to-probe-to-victim replay reconciles injected energy/current, field/circuit power, interventions, state transitions and recovery; missed, ambiguous, silent, or noncausal events become PossibleEvent/Unknown/Refuted rather than a favorable trace",
            hypotheses: &[
                "injection waveform/source impedance/coupling clamp or field port, calibration, location/orientation, clock and event order, uncertainty, and admitted rung are explicit",
                "victim inputs, internal state, functional class, thresholds, hysteresis, sampling/quantization/delay, upset/latch/reset/recovery maps, telemetry and fault-containment boundary are content-addressed",
                "event-time and reset claims use the admitted hybrid/event class; simultaneous, grazing, chatter, Zeno, missing telemetry, or nonlinear model gaps retain set-valued/Unknown semantics",
                "causal attribution requires the typed source-coupling-probe-victim chain and intervention/counterfactual assumptions; observational correlation alone is never relabeled causal",
            ],
            qoi: "source_to_victim_power_event_upset_and_recovery_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.victim_chain.v1 at fs-vmanifest-oracles/i14/victim_chain.rs::replay_interventions",
                independent: true,
                tcb_overlap: "shares frozen waveform, port, and victim automaton bytes; independent event ledger, power integration, and intervention replay do not share production compiler/runtime",
            },
            activation: "source/probe and victim/injection contracts plus one admitted coupling rung are green",
            kill: "hidden source, double-counted power, missed or order-dependent upset, silent victim, noncausal recovery, ambiguous event promoted to Pass, or stale victim version kills the chain",
            fallback: "report bounded probe exposure and ModeLedger Unknown/PossibleEvent; safety owner receives the unresolved transition",
            no_claim: "field amplitude or simulated recovery alone is not product immunity, causal hardware validation, a safety conclusion, or regulatory certification",
        },
        ClaimSpec {
            id: "i14-fixed-regime-adjoint-closure",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For one admitted fixed EMC rung, topology, mesh, active set, event mode, material branch, reference/gauge, and router decision, TraceAdjoint and ImplicitResidualAdjoint remain distinct typed derivatives. Complex states use a frozen real-variable or Wirtinger convention, conjugate-transpose metric, and real scalar QoI. Tangent-adjoint duality, multiscale Taylor remainder, independent differences/duals, primal/adjoint residual and inexact-solve terms close for routes, spacing, shield/bond/filter/ground, aperture, material and smooth-geometry variables; a switch returns NoGradient, SetValued, or Unknown",
            hypotheses: &[
                "the exact recorded algorithm or converged residual, state/operator hashes, real objective, complex derivative convention, inner products, units/scales, solver tolerances, preconditioner/stopping state, held variables and objective count are frozen",
                "TraceAdjoint differentiates the exact recorded branch/iteration trace; ImplicitResidualAdjoint separately proves residual differentiability, state-Jacobian invertibility/stability, and inexact primal/adjoint correction",
                "finite differences or forward duals use an independently assembled path and conditioning/truncation enclosure; complex-step is forbidden unless the entire geometry/material/solve/QoI path is holomorphic and preserves the perturbation",
                "topology, remesh, active-set, event, router/rung, eigen-subspace, absolute-value/norm kink, or constitutive branch changes have no classical gradient without a separately proved generalized derivative",
            ],
            qoi: "maximum_preregistered_normalized_tangent_adjoint_taylor_and_independent_derivative_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.adjoint.v1 at fs-vmanifest-oracles/i14/adjoint.rs::check_fixed_regime_derivative",
                independent: true,
                tcb_overlap: "shares frozen primal values and objectives; separately assembled tangent/dual and interval-controlled differences do not call production reverse path",
            },
            activation: "one core or maximal rung has an admitted differentiable fixed-regime residual or exact trace",
            kill: "wrong conjugation/metric, trace/implicit authority conflation, missing residual/inexact term, nonholomorphic complex-step, branch leakage, or normalized derivative score above one kills gradient promotion",
            fallback: "use derivative-free or mixed-discrete robust search with the discontinuity/rung boundary explicit",
            no_claim: "exact discrete derivative means exact only for the frozen trace or residual theorem; it is not model truth, cross-rung differentiability, topology derivative, or parameter-independent end-to-end cost",
        },
        ClaimSpec {
            id: "i14-emc-uq-inference-mechanics",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "EMC inference mechanics keep source variability, harness geometry/tolerance, contact/bond state, material dispersion, calibration/measurement uncertainty, spatial/process dependence, numerical error, model-form discrepancy, and fidelity-router uncertainty separate. Deterministic ensembles, QMC, interval/affine enclosures, multifidelity control variates, rare-event estimators, robust/DRO/chance/worst-case objectives, and anytime-valid e-process/confidence procedures activate only under their exact synthetic or declared assumptions; tail or discrepancy insufficiency escalates or returns Unknown/conservative bounds",
            hypotheses: &[
                "every random/set-valued input declares population or ambiguity set, support, units, dependence/correlation, calibration lineage, sampling measure, importance weights, exchangeability/martingale/mixing assumptions, and drift/revalidation trigger",
                "candidate/model/toolchain, objectives, guards, stopping/multiplicity rule, sample identities, fidelity policy, holdout custody and artifact roots freeze before adaptive evaluation; public deterministic fixtures provide replay but no IID authority",
                "optional stopping and adaptive allocation use the declared anytime-valid procedure; tail effective sample size, importance-weight stability, discrepancy, extrapolation and numerical bounds are monitored and cannot be replaced by nominal percentiles",
                "epistemic model-rung disagreement and measurement uncertainty remain separate from aleatory variability and never average into a categorical evidence color",
            ],
            qoi: "uq_coverage_optional_stopping_tail_and_escalation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.uq.v1 at fs-vmanifest-oracles/i14/uq.rs::audit_sampling_and_coverage",
                independent: true,
                tcb_overlap: "shares sealed observations and distribution declarations; independent analytic cases, null simulations, interval containment, stopping and weight audits do not share production sampler/planner",
            },
            activation: "core rung receipts and the public replay inference-mechanics decks are frozen; no governed population is required or implied",
            kill: "holdout leakage, undeclared dependence, invalid optional stopping/multiplicity, unstable weights, tail or discrepancy budget failure, evidence laundering, or a nominal result after required escalation kills decision-grade UQ",
            fallback: "return conservative interval/worst-case bound or Unknown and rank the missing data/model/rung evidence",
            no_claim: "this Core authority proves inference, stopping, containment, and escalation mechanics only: no distribution-free guarantee for arbitrary dependent/nonstationary processes, physical validation, compliance probability, physical-population reliability, or truth of a population model",
        },
        ClaimSpec {
            id: "i14-fullwave-feec-stability-energy",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Inside an admitted full-wave frequency- or time-domain class, mixed FEEC Maxwell fields on a discrete de Rham complex satisfy exact incidence separately from bounded commuting-projection/discrete-compactness or formulation-specific stability, gauge/harmonic and boundary-complex treatment, charge/divergence/source compatibility, and h/p spectral-pollution obligations. Frequency-domain driven/eigen problems close complex power/Poynting ledgers. Time-domain lossless/passive-dispersive problems close fully discrete Gauss, charge continuity, source/boundary/PML-aware Poynting, auxiliary-state storage/dissipation, nonlinear/iterative defect, CFL/implicit stability, and numerical-dispersion receipts",
            hypotheses: &[
                "FullWaveProblem convention admission is green and binds mesh regularity, element family/order, material regularity/passivity/dispersion, topology, relative boundary complex, harmonic sector, gauge/divergence constraint, source class and QoI",
                "exact incidence is not used as a no-spurious-mode theorem; the admitted family proves a bounded commuting projection and discrete compactness or its exact coercivity/inf-sup substitute plus source/boundary compatibility",
                "frequency-domain real-block signs and RMS power are independently checked; time-domain ADE/state initialization, passive pole/residue convention, timestep/nonlinear solve, boundary/PML energy, and fully discrete rather than merely semidiscrete inequality are explicit",
                "PML qualification covers thickness/grading, mesh/p order, corners/edges, grazing, evanescent, anisotropic/dispersive and discontinuous media, late-time growth and reflection; failure demotes without hiding the boundary energy",
            ],
            qoi: "maximum_preregistered_normalized_fullwave_stability_gauss_spectral_power_and_dispersion_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.fullwave_feec.v1 at fs-vmanifest-oracles/i14/fullwave_feec.rs::adjudicate_stability_and_balance",
                independent: true,
                tcb_overlap: "shares mesh/material/source bytes and exact incidence only; analytic cavity/waveguide/interface fields, independent mixed operators, spectra and energy integration do not share production assembly/solver/PML",
            },
            activation: "the full-wave convention contract, FEEC family, material card, boundary and solver capabilities are enabled behind maximal gates",
            kill: "spurious retained mode, failed commuting/stability/gauge/source premise, divergence/Gauss/charge error, unclosed power/energy, nonpassive ADE, unresolved late PML growth, or normalized score above one kills the affected full-wave class",
            fallback: "lower admitted EQS/MTL/PEEC rung, explicit bounded-domain solve, refined boundary/PML, or Unknown/refusal",
            no_claim: "no arbitrary broadband, high-frequency, full-vehicle, active/hysteretic/nonlinear-media, exact-PML, or physical-validation authority outside the pinned class and scale",
        },
        ClaimSpec {
            id: "i14-exterior-bem-formulation-correctness",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Exterior Maxwell BEM admits three formulation classes without conflation: closed PEC conductors use explicitly selected EFIE/MFIE/CFIE with typed tangential traces; closed penetrable dielectric interfaces use an explicitly selected PMCHWT, Mueller, or JMCFIE transmission formulation with both material-side traces; separately admitted open PEC screens use a screen EFIE with screen trace spaces and edge-singularity treatment. Each class has orientation-sensitive singular/near-singular quadrature, resonance/low-frequency routing, and explicit FEEC-BEM continuity/power coupling. Independently assembled dense BEM on pinned small problems is the numerical correctness oracle; acceleration authority is separate",
            hypotheses: &[
                "each surface is orientable, manifold and inside the regularity/corner class required by its formulation; closed PEC, closed dielectric-transmission, and admitted open-screen domains are distinct, and normal, sided trace, surface-current basis, edge space, kernel/sign, wavenumber/branch, frequency/material and source/incident field are frozen",
                "singular and near-singular quadrature, selected PEC or dielectric or screen formulation, interior-resonance treatment, low-frequency conditioning, FEEC trace coupling and outward Poynting convention have independent receipts",
                "dense assembly, solve, trace, near/far field and power receipts separate formulation, geometry, discretization, singular quadrature and algebraic error",
            ],
            qoi: "maximum_preregistered_normalized_bem_formulation_trace_power_and_dense_qoi_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.maxwell_bem.v1 at fs-vmanifest-oracles/i14/maxwell_bem.rs::dense_trace_and_farfield_check",
                independent: true,
                tcb_overlap: "shares surface/source bytes and scalar low-level arithmetic; independent dense kernel/quadrature, analytic PEC sphere/dipole and power checker do not share production BEM assembly",
            },
            activation: "full-wave convention admission is green and the selected BEM formulation, trace spaces and independent dense reference are available",
            kill: "trace/orientation/power mismatch, hidden resonance/low-frequency failure, wrong-side dielectric trace, open screen admitted without screen spaces/edge treatment, nonmanifold surface, or dense correctness miss kills the formulation authority",
            fallback: "a different admitted dense formulation, bounded-domain FEEC/PML, lower rung, or refusal",
            no_claim: "no arbitrary open screens beyond the admitted screen-EFIE class, unadmitted corners/material interfaces, universal low-frequency stability, FMM acceleration, exact exterior truncation, or physical validation",
        },
        ClaimSpec {
            id: "i14-maxwell-fmm-acceleration-envelope",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A deterministic Maxwell FMM accelerator receives authority only when its formulation-bound vector/tangential Helmholtz kernel, near/far split, complex representation, tree, expansion/interpolation order, frequency/electrical-size regime, adjoint, cancellation tiling, and error allocation independently enclose dense BEM matvecs and solved near/far-field QoIs. Acceleration, crossover, and performance evidence remain orthogonal to BEM formulation correctness and physical validation",
            hypotheses: &[
                "the exact admitted dense BEM operator, EmConventionCard, surface/basis order, kernel branch, trace spaces and requested field/port QoIs are frozen",
                "tree construction, expansion or interpolation family/order, near/far criterion, vector/tangential coupling, tolerance allocation, dense crossover and deterministic traversal are versioned",
                "every accelerated matvec and solved QoI carries an independent dense baseline and certified outward error envelope separating FMM, iteration, quadrature, geometry and discretization error",
                "performance promotion binds p50/p95/p99 time and peak memory, machine fingerprint, scale, achieved QoI error and measured crossover; a tighter request need not yield monotonically smaller observed error after an algorithmic route change",
            ],
            qoi: "maximum_preregistered_normalized_fmm_dense_matvec_solved_qoi_and_performance_envelope_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i14.oracle.maxwell_fmm.v1 at fs-vmanifest-oracles/i14/maxwell_fmm.rs::dense_envelope_and_crossover_check",
                independent: true,
                tcb_overlap: "shares frozen surface, basis and kernel card only; independent dense matvec/solve and measurement harness do not call the production hierarchy, traversal or expansion code",
            },
            activation: "the exact BEM formulation authority is green and the FMM gate plus independent dense envelope are available",
            kill: "one dense-envelope escape, wrong vector/tangential coupling, adjoint mismatch, nondeterministic tree, cancellation-bound violation, unmeasured crossover, or performance claim outside its fingerprint kills FMM authority",
            fallback: "use the admitted dense BEM operator, a narrower frequency/scale/order envelope, or refusal",
            no_claim: "FMM agreement does not repair a wrong BEM formulation or discretization and proves no universal high-frequency superiority, exterior exactness, or physical validation",
        },
        ClaimSpec {
            id: "i14-certified-fidelity-descent-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked EMC fidelity-descent theorem composes local one-way or bidirectional EQS/MTL/PEEC/full-wave Crosswalk morphisms into a global requested-QoI error majorant when each state/port map, approximation remainder, stability constant, reference/gauge transform, loss owner, and model/domain overlap satisfies the frozen hypotheses. The theorem preserves an evidence preorder: routing may change cost and representation but cannot strengthen epistemic authority. Shared TCB and correlated defects remain joint terms rather than root-sum-square fiction; refinement/escalation is monotone only in the declared bound preorder",
            hypotheses: &[
                "each rung is an object with typed state/port spaces, reference conductor, topology, material/source/QoI domain, storage/dissipation/radiation ownership, validity predicate, cost and evidence state",
                "each crosswalk supplies commuting charge/power diagrams, exact orientation/unit maps, projection/lift, local directed naturality defect, stability/reliability constant, discrepancy dependence graph and no-claim boundary; no missing higher-mode/radiation/loss term is treated as zero",
                "the composed graph is finite and well posed; local majorants and correlation/joint-enclosure semantics satisfy the theorem's exact hypotheses, with positive denominators and no signed net work as a scale",
                "a canonical binding receipt proves byte/semantic equivalence among manifest claim, machine theorem-card AST, generated formal proposition, elaborated declaration and runtime-premise schema; the pre-proof successor freezes exact translation and axiom closure before theorem authority exists",
            ],
            qoi: "independent_fidelity_descent_theorem_and_runtime_majorant_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.fidelity_descent.v1 composite proofs/i14/FidelityDescent.lean::emcFidelityDescent plus fs-vmanifest-oracles/i14/fidelity_descent.rs::check_runtime_majorant",
                independent: true,
                tcb_overlap: "shares bound theorem AST, crosswalk bytes and rational constants; pinned proof kernel and separately implemented interval/dependence majorant checker do not share production router",
            },
            activation: "all used rung/crosswalk receipts exist and a pre-candidate successor has frozen the complete proposition/definition/runtime-premise AST, deterministic translation, declaration identity, axiom policy and nonvacuity family required by i14-theorem-formalization-policy",
            kill: "proof/declaration binding rejection, inadmissible runtime premise, observed QoI error outside the joint bound, hidden/correlated omitted term, evidence strengthening, or one independently admitted countermodel refutes exactly the bound theorem revision",
            fallback: "retain pairwise discrepancies and fixed-rung execution with no global descent authority",
            no_claim: "adjacent-rung agreement or finite benchmark survival is not this theorem, not physical validation, and not a proof that the runtime instance satisfies its premises",
        },
        ClaimSpec {
            id: "i14-robust-mitigation-heldout",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "Against frozen catalog, ablation, and fixed-budget baselines, at least one preregistered harness route, shield/termination, common-mode filter, bond/ground, aperture/enclosure, or mixed mitigation improves every directed governed-blind high-fidelity EMC decision QoI required by its objective under the declared uncertainty evidence while satisfying mass, cost, thermal, mechanical, manufacturing, serviceability, signal-integrity, safety and system-function guards. All design/rung/topology switches, total compute, feasible incumbents and Pareto lineage are retained, and an independent adjudicator reconstructs and reruns the sealed simulation cases. Governed physical/laboratory mitigation validation is a separate authority that additionally requires laboratory calibration and as-built specimen packs",
            hypotheses: &[
                "GovernanceCommitted freezes design variables, topology choices, fixed-regime gradient domains, objectives, directed comparison arithmetic, guards, uncertainty/dependence sets, baselines, evaluation budget, optimizer family, fidelity escalation, activation/kill rules, blind generator and custody protocol before candidate execution",
                "after CandidateFrozen, one atomic RealizationCommitted authority transaction replaces i14-mitigation-max-holdout with the typed inaccessible realized root, installs i14-external-blind-mitigation-custody-pack as its same-ID External discharge envelope, removes its Waiver row through a verified fs-vvreg DischargeReceipt, verifies the AmendmentRecord and advances the authority head; waiver-only, slot-only, envelope-only, split, public-deterministic or raw-digest changes cannot discharge this statistical/untouched obligation",
                "classical gradients are consumed only inside their fixed-rung regular domains; mixed discrete/topology/router/event changes use derivative-free, one-sided or set-valued methods with no false gradient crossing",
                "winner's curse, adaptive selection and multiplicity use a declared nested-holdout/anytime-valid adjudication; every guard uses canonical cross-domain artifacts and remains distinct from objective improvement",
            ],
            qoi: "blind_robust_mitigation_improvement_all_guards_and_independent_reproduction_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.robust_mitigation.v1 at fs-vmanifest-oracles/i14/robust_mitigation.rs::blind_reconstruct_and_compare",
                independent: true,
                tcb_overlap: "shares frozen objective/guard definitions and sealed final designs; independent geometry reconstruction, baseline runner, high-fidelity simulation selection and statistical adjudication do not share production optimizer/surrogate",
            },
            activation: "core receipt, UQ inference mechanics, fixed-regime adjoint or declared derivative-free path, guards, and the atomic slot plus same-ID envelope plus waiver-discharge RealizationCommitted transaction are green",
            kill: "holdout leak, post-hoc objective/guard/band, invalid gradient/rung crossing, baseline unfairness, one guard failure, missing compute/Pareto lineage, or failed independent improvement/refusal kills promotion of the mitigation",
            fallback: "retain best feasible baseline and every negative result/counterexample; no-improvement is an honest terminal state",
            no_claim: "optimization and simulation do not establish product compliance, laboratory immunity/emissions, manufacturing feasibility beyond the checked guards, or regulatory approval",
        },
        ClaimSpec {
            id: "i14-emc-safety-case-integration",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "One flagship source-to-coupling-to-probe-to-victim-to-upset/recovery-to-mitigation chain binds scoped EMC evidence into an existing assurance case without laundering it into a safety or regulatory conclusion. Every HazardId assumption, functional requirement, operating envelope, fault-containment region, monitor, violation effect, owner, expiry and revalidation trigger is explicit; bearing-current/insulation/electric-shock evidence, source/coupling uncertainty, experimental validation, standards conformance and residual common-mode risk remain separate typed edges",
            hypotheses: &[
                "the exact source, harness/as-built, rung/crosswalk, probe, victim/controller, event, mitigation, uncertainty, operating envelope and hazard/requirement revisions are content-addressed",
                "scientific Verified, empirical Validated, Estimated, StandardConformance, safety-process, and regulatory-approval states are orthogonal; the assurance join is monotone and invalid/refuted evidence dominates favorable observations",
                "the source-to-victim causal chain names interventions and model assumptions; hazard conclusion is owned by the safety program and simulation evidence can support, weaken, or leave it Unknown but cannot self-certify it",
                "all external standard editions, lab setups, instrumentation uncertainty, physical specimens, calibration/validation split, signatures/custody and license constraints are pinned and their corresponding governed producer receipts independently adjudicated, or the edge stays waived/Unknown/NoPromotionAuthority",
            ],
            qoi: "scoped_emc_evidence_to_hazard_traceability_and_no_laundering_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.assurance.v1 at fs-vmanifest-oracles/i14/assurance.rs::audit_hazard_edges",
                independent: true,
                tcb_overlap: "shares claim/receipt and HazardId schemas; independent assurance graph reconstruction and authority-color checker do not share production simulator/router or safety adjudicator",
            },
            activation: "a scoped core or maximal EMC chain and the target assurance-case interface exist; each external physical or standards edge activates only after its complete atomic discharge transaction and the corresponding independently adjudicated governed producer receipt are both green, otherwise that edge remains Unknown/NoPromotionAuthority",
            kill: "stale/missing hazard edge, assumption without owner/monitor/expiry, scope or evidence-color laundering, untraceable source-to-victim link, hidden common-mode dependence, or simulation-only compliance/safety verdict kills the integration claim",
            fallback: "publish the scoped EMC evidence package with the safety/regulatory conclusion Unknown and ranked missing edges",
            no_claim: "this is traceable assurance evidence, never legal compliance, product certification, regulatory approval, or proof that the physical hazard model/population is true",
        },
        ClaimSpec {
            id: "i14-governed-standards-crosswalk",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For exact governed editions and licensed clause subsets, a standards crosswalk binds native HarnessGraph/AP242 adapter semantics and scoped EMC or rotating-machine test requirements to content-addressed requirement, evidence, exception, and loss-accounting edges. Every edition, corrigendum, option, unit/frame transform, occurrence map, applicability predicate, test configuration, clause owner, and unimplemented construct is explicit; a machine-verifiable clause graph returns ConformantToScopedSubset, Nonconformant, Unknown, or Refused without embedding or paraphrasing restricted text into public artifacts",
            hypotheses: &[
                "fs-vvreg supplies authenticated licensed bytes or a legally distributable clause digest graph for each exact edition/corrigendum and records license, custodian, access, signature and clause scope",
                "the AP242 subset maps occurrence reuse, assembly transforms, units, frames, connectivity, open/short/shield semantics, material/property lineage and every information loss; unsupported or ambiguous constructs refuse rather than disappear",
                "EMC and machine-test crosswalks bind exact configuration, detector/bandwidth, calibration, operating mode, limit arithmetic, uncertainty treatment, evidence color and exclusions; a schema-shaped public deck cannot discharge licensed-edition authority",
            ],
            qoi: "scoped_exact_edition_clause_crosswalk_and_loss_accounting_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.standards_crosswalk.v1 at fs-vmanifest-oracles/i14/standards_crosswalk.rs::reconstruct_clause_and_adapter_edges",
                independent: true,
                tcb_overlap: "shares authenticated edition/clause and native artifact bytes; independently reconstructs clause applicability, AP242 occurrence/unit/frame maps, exclusions and loss receipts without production importer or assurance traversal",
            },
            activation: "the GovernedStandards GovernanceCommitted, AuthorizedConstruction, CandidateFrozen and StandardsAuthorityCommitted transaction path closes before independent adjudication, and every referenced native synthetic authority is green",
            kill: "edition or corrigendum ambiguity, missing licensed digest, unmapped supported clause, silent AP242 loss, configuration mismatch, favorable treatment of Unknown, or cross-edition result reuse kills the scoped crosswalk",
            fallback: "retain synthetic adapter mechanics and an explicit Unknown external edge with the missing edition, clause, mapping or evidence requirement",
            no_claim: "a scoped machine crosswalk is not legal advice, blanket AP242 implementation, laboratory validation, product certification, regulatory approval, or conformance outside the authenticated clauses and configurations",
        },
        ClaimSpec {
            id: "i14-governed-laboratory-emc-validation",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a frozen as-built specimen and operating envelope, an independently governed laboratory campaign compares preregistered source, terminal, field, current, spectrum, victim-mode, energy/loss, and uncertainty QoIs against the exact corresponding simulation artifacts. Calibration chains, fixture/cable/chamber state, environmental history, instrumentation transfer functions, censoring, missingness, repeats, exclusions, model updates, and validation-versus-calibration separation are retained; the result is a scoped physical-validation vector and never a generic green badge",
            hypotheses: &[
                "GovernanceCommitted freezes calibration/model-input/validation stratum rules, laboratory procedure, instrument/calibration policy, uncertainty budget, environment, custody, AcceptanceCards and adjudicator; an independent custodian/fs-vvreg capability first creates salted or equivalently hiding calibration/model-input content commitments, a hiding validation source-universe/frame commitment, disjoint-membership commitment, exact validation-selection algorithm and pre-candidate secret-seed/VRF or equivalent non-adaptive selection commitment. Candidate builders, fitters, checker/threshold owners and their transitive capabilities receive only opaque commitment identities, never validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material. AuthorizedCalibration exposes only committed calibration and as-built model-input strata under a complete read/output ledger; AsBuiltModelInstantiation performs only preregistered fitting/model selection; CandidateFrozen then binds the resulting immutable model/toolchain and every validation rule before any candidate-side protected validation access. The one RealizationCommitted transaction binds separately addressable calibration, model-input and validation roots, membership proofs to pre-access commitments, mutual disjoint-membership proof, non-adaptive selection proof, a contamination receipt naming every audited principal and transitive capability, same-ID discharge envelopes and distinct receipts",
                "calibration data, parameter fitting, model selection and threshold choice are disjoint from validation cases or use a preregistered nested design with multiplicity/optional-stopping control",
                "every comparison consumes an AcceptanceCard and reports ClaimAdjudication, DomainApplicability, EvidenceCompleteness and EvidenceIntegrity independently; missing or censored traces cannot be silently dropped and IntegrityFailed is never a scientific claim verdict",
            ],
            qoi: "scoped_asbuilt_laboratory_emc_validation_vector_and_integrity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.lab_validation.v1 at fs-vmanifest-oracles/i14/lab_validation.rs::reconstruct_calibration_and_qoi_vector",
                independent: true,
                tcb_overlap: "shares governed raw/specimen/configuration bytes and frozen AcceptanceCards; independent calibration propagation, trace reconstruction and comparison do not call production simulator, fitter or threshold owner",
            },
            activation: "one exact RealizationCommitted authority transaction replaces the laboratory commitment slot, installs both same-ID discharge envelopes, retires both calibrated-laboratory and proprietary-asbuilt-geometry waivers, and verifies their distinct receipts before joint reveal",
            kill: "specimen/configuration drift, calibration or custody gap, calibration-validation leakage, undeclared exclusion, unpropagated uncertainty, post-hoc QoI/scale, or favorable collapse of a failed/unknown component kills the scoped validation vector",
            fallback: "publish calibrated traces and a componentwise Unknown/Failed validation vector with ranked missing authority; synthetic solver authority remains orthogonal",
            no_claim: "one laboratory campaign is not population reliability, standards conformance unless separately crosswalked, universal physical truth, product certification, safety proof, or regulatory approval",
        },
        ClaimSpec {
            id: "i14-production-bearing-population-reliability",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a preregistered production bearing, lubricant, material, machine, drive, installation, duty-cycle, environment, maintenance, and censoring population, a governed campaign estimates the declared bearing-current event, shaft-voltage, discharge-energy, damage-proxy, and time-to-event reliability function while preserving aleatory, epistemic, model-discrepancy and selection uncertainty. Simulation is used only as a frozen covariate or prior component; population conclusions require independently governed observations, coverage, missingness and anytime-valid stopping evidence",
            hypotheses: &[
                "GovernanceCommitted freezes the exact inclusion/exclusion frame, sampling design, lot and configuration lineage, duty/environment exposure, telemetry, inspection/damage definitions, censoring/competing-risk rules, missingness model, multiplicity, stopping rule and target estimand, and CandidateFrozen closes before protected observation or adjudication access",
                "a separately governed bearing-metrology pack establishes bearing-specific shaft-voltage/current sensor, discharge-event classifier, inspection/damage proxy, timing and uncertainty behavior for the claimed range; generic EMC chamber/cable calibration, laboratory traces and PublicReplay fixtures do not substitute for this metrology or for production-population sampling",
                "coverage and tail claims use an independently checked finite-sample, e-process, conformal, survival, partial-identification or other explicitly admitted method whose assumptions and sensitivity bounds remain visible",
            ],
            qoi: "production_bearing_population_reliability_coverage_and_integrity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.bearing_population.v1 at fs-vmanifest-oracles/i14/bearing_population.rs::audit_frame_events_and_coverage",
                independent: true,
                tcb_overlap: "shares governed population records and frozen estimand only; independent frame/censoring/event/coverage audit does not call production hybrid solver, optimizer or reliability estimator",
            },
            activation: "one exact RealizationCommitted authority transaction replaces the bearing-population commitment slot, installs both same-ID discharge envelopes, retires the production-population and bearing-metrology waivers, and verifies their distinct receipts; the synthetic bearing-path claim is green for the declared covariate range",
            kill: "convenience-sample laundering, lot/configuration ambiguity, unmodeled censoring or competing risk, event-definition drift, optional-stopping/multiplicity defect, invalid coverage, or extrapolation outside the population frame kills the reliability claim",
            fallback: "retain specimen- or campaign-scoped observations and synthetic bearing-current paths with population reliability Unknown",
            no_claim: "population reliability does not itself prove causal damage, universal bearing life, product safety, standards conformance, warranty fitness, certification, or regulatory approval",
        },
        ClaimSpec {
            id: "i14-governed-emc-reliability-validation",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a preregistered source, harness, bond/shield, victim, operating-envelope and deployment population, a governed campaign estimates declared emissions, susceptibility, upset/recovery and protection reliability functionals while keeping aleatory variation, epistemic uncertainty, measurement error, model discrepancy, numerical error, missingness and selection effects separate. Simulation is a frozen covariate, control variate or prior component only; population authority comes from independently governed observations and exact finite-sample, survival, conformal, e-process or partial-identification assumptions",
            hypotheses: &[
                "GovernanceCommitted freezes the sampling frame, inclusion/exclusion, configuration and lot lineage, operating and environmental exposure, event/outcome definitions, measurement chain, censoring/missingness/dependence, target estimand, multiplicity and stopping rule, and CandidateFrozen closes before protected observation or adjudication access",
                "one RealizationCommitted authority transaction replaces the GovernedPopulationSlot, installs the population waiver subject's same-ID discharge envelope, retires its Waiver row and verifies the exact realized root, typed DischargeReceipt and AmendmentRecord under i14-campaign-policy-v1 before joint reveal",
                "calibration, model selection, discrepancy fitting and threshold choice are disjoint from validation observations or use a preregistered nested design whose dependence and optional stopping remain valid",
                "every AcceptanceCard component reports ClaimAdjudication, DomainApplicability, EvidenceCompleteness and EvidenceIntegrity independently; public replay, a synthetic UQ deck, or one laboratory specimen cannot substitute for the target population, and IntegrityFailed is never a scientific claim verdict",
            ],
            qoi: "governed_emc_population_reliability_coverage_and_integrity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.emc_reliability.v1 at fs-vmanifest-oracles/i14/emc_reliability.rs::audit_frame_events_and_anytime_coverage",
                independent: true,
                tcb_overlap: "shares governed observations, frozen estimand and AcceptanceCards only; independent frame, missingness, stopping, event and coverage audit does not call production simulator or reliability estimator",
            },
            activation: "the governed EMC population slot, same-ID discharge envelope, Waiver row, typed receipt and AmendmentRecord close in one RealizationCommitted pre-reveal authority transaction and the relevant synthetic source/path/victim authorities are green",
            kill: "population-frame drift, access before candidate freeze, calibration-validation leakage, undeclared dependence or exclusion, missingness/censoring defect, invalid optional stopping/multiplicity, miscoverage, or out-of-frame extrapolation kills reliability authority",
            fallback: "retain synthetic inference-mechanics and specimen-scoped observations while reporting population reliability Unknown",
            no_claim: "this scoped population result is not standards conformance, universal immunity, product safety, certification, warranty fitness, causal attribution, or authority outside the frozen frame",
        },
        ClaimSpec {
            id: "i14-passive-causal-sheaf-composition-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked passive-causal composition theorem treats local harness, field, MTL, PEEC, circuit, shield and victim models as a paired compatibility sheaf and balance cosheaf over a finite cellular assembly/regime site. Under exact typed hypotheses, restriction/crosswalk maps preserve oriented port variables and relative-hypercohomology winding/cut representatives, aggregation maps conserve charge/flow, higher-overlap cocycles cohere, and a chain-level natural power pairing into the orientation/dualizing cosheaf satisfies restriction/corestriction adjointness as a discrete Green-Stokes identity. Causal dissipative component relations interconnected by a maximally isotropic lossless Dirac relation with respect to a frozen nondegenerate ambient split-signature power pairing and a well-posed feedback/descriptor closure yield a quantified finite- and infinite-horizon global storage/supply inequality. Zero unaccounted defect is required for an exact loss identity; with bounded accounted gluing/time-discretization defects, non-strict dominance of certified integrated dissipation over the dependency-aware worst-case violation bound yields passivity, while strict dominance by a frozen positive signal margin yields quantified robust strict passivity",
            hypotheses: &[
                "each local model has versioned causal operator/relation semantics on a declared signal space, state/storage functional, supply rate, dissipation/loss owners, initialization/history, reference/gauge and validity domain; active sources/reservoirs are separate boundary ports",
                "the cellular site/cover, coefficient rings and torsion policy, stalk/state/effort/flow objects, restriction and corestriction maps, pairwise and higher-overlap cocycle/coherence laws, orientations/dualizing cosheaf, reference conductor, relative total chain/cochain complex, induced maps, kernels/cokernels, winding/cut classes, trace maps and power metrics are explicit; gauge/tree-cotree/cut representative changes are related by checked chain homotopies that preserve complete physical cochains, sources, supply and observables",
                "relative-boundary compatibility is checked on terminals, exterior radiation boundaries, interfaces and excised sources; the interconnection is a maximally isotropic power-conserving Dirac relation with respect to a nondegenerate ambient split-signature pairing, with clean composition, no unpaired latent port and a separately proved regular/well-posed descriptor closure free of unowned impulsive modes; delays, sampling, switching, multirate windows, numerical iteration and time discretization carry their own stability/passivity or defect premises",
                "for every admitted trajectory in the frozen signal/history space and every finite T>=0, storage S is lower bounded and the theorem proves S(x(T))-S(x(0)) <= integral_[0,T] <e_ext,f_ext> dt - integral_[0,T] d(t) dt + Delta_accounted(T), with d>=0 and a dependency-preserving outward enclosure on Delta_accounted. Exact global loss identity requires Delta_unaccounted=0 and equality under its stated owners. Defect-tolerant passivity requires m(T)=integral d-Delta_accounted >=0 over all admitted trajectories/horizons; robust strict passivity additionally requires m(T)>=mu*norm_signal_[0,T]^2 for one frozen mu>0. The infinite-horizon conclusion separately freezes coercivity/detectability and limit/integrability premises",
                "a matrix/frequency-domain positive-real component enters this time-domain storage theorem only through a checked KYP/passive-realization bridge that covers generalized boundary poles, quotient modes and initialization; sampled-band dissipativity alone is inadmissible",
                "cover refinement and representative invariance are proved by a checked nerve/cofinality or equivalent comparison theorem on the relative total complex; nonzero hypercohomology or crosswalk-holonomy obstruction classes force escalation rather than being averaged into a defect budget",
                "a canonical binding receipt proves byte/semantic equivalence among the manifest claim, machine theorem-card AST, definitions, generated formal proposition, elaborated declaration, runtime-premise schema, proof term and complete transitive axiom closure under the exact policy in i14-theorem-formalization-policy",
            ],
            qoi: "independent_passive_causal_sheaf_cosheaf_composition_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.passive_composition.v1 composite proofs/i14/PassiveCausalComposition.lean::assemblyPassivity plus fs-vmanifest-oracles/i14/passive_composition.rs::check_runtime_diagrams",
                independent: true,
                tcb_overlap: "shares bound theorem AST, rational finite baselines, incidence and declared port metrics; pinned proof kernel and independent diagram/power/runtime-premise checker do not share production composition/router",
            },
            activation: "the theorem target card is frozen and a pre-candidate successor freezes complete proposition/definition/runtime-premise ASTs, deterministic translation, exact declaration and axiom closure, and strength-matched nonvacuity families before proof/search/runtime promotion",
            kill: "proof/declaration binding rejection, failed higher-overlap/refinement or relative-boundary coherence, noncommuting compatibility/balance/power diagram, non-maximal Dirac relation, degenerate ambient power pairing, unclean composition, ill-posed or noncausal closure, missing KYP/passive-realization bridge, unowned source/loss/impulse, negative trajectory/horizon margin after dependency-aware debit, unresolved hypercohomology/holonomy obstruction, chain-homotopy-dependent physical observable, or an AuthorityContradiction blocks the exact bound revision",
            fallback: "retain component/rung receipts and direct global charge/power audit with no compositional theorem authority",
            no_claim: "component passivity, a Dirac interconnection, exact incidence, bare sheaf/cohomology equivalence, or sampled cross-rung agreement alone proves neither global passivity, causality, representative invariance, time-discrete stability, nor physical validity",
        },
        ClaimSpec {
            id: "i14-hypercohomology-obstruction-localization-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked relative-hypercohomology obstruction theorem constructs the total compatibility/balance complex for a finite typed harness assembly and proves that a nonzero obstruction class certifies the impossibility of globally gluing the requested local port/field/cochain data under the frozen crosswalks. When the exact acyclicity, descent and boundary hypotheses hold, vanishing plus a constructive witness yields a global section; persistent representatives localize a minimal certified family of seams, ports, winding/cut classes or model transitions whose repair can remove the obstruction",
            hypotheses: &[
                "the finite cellular site or cover, coefficient ring and torsion policy, relative boundary, local complexes, restriction/corestriction maps, crosswalks, total differential, sign convention and requested observable are canonical and executable",
                "the proposition distinguishes necessary obstruction vanishing from sufficient gluing hypotheses; it never promotes bare zero cohomology without exactness, acyclicity or constructive witness premises",
                "representative localization is invariant under admitted chain homotopy and reports all tied minimal supports under deterministic order rather than inventing one unique physical cause",
                "a canonical theorem binding, axiom audit, runtime-premise checker, nonvacuity family and adversarial nonzero classes satisfy i14-theorem-formalization-policy-v1",
            ],
            qoi: "hypercohomology_obstruction_gluing_witness_and_minimal_support_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.hypercohomology_obstruction.v1 composite proofs/i14/HypercohomologyObstruction.lean::relativeGluingObstruction plus fs-vmanifest-oracles/i14/hypercohomology.rs::check_total_complex_and_support",
                independent: true,
                tcb_overlap: "shares canonical finite complexes and theorem AST; pinned proof kernel and independent Smith/field-linear total-complex checker do not call production router or repair search",
            },
            activation: "the obstruction target card and theorem policy are frozen and a successor binds the complete total-complex proposition, runtime schema and nonvacuity families",
            kill: "d_total squared nonzero, wrong relative boundary or torsion treatment, false vanishing/sufficiency direction, noninvariant support, missed tied repair, theorem binding failure, or AuthorityContradiction blocks the theorem revision",
            fallback: "retain explicit seam/crosswalk residuals and direct global solve with obstruction authority Unknown",
            no_claim: "vanishing obstruction alone is not physical validity, global passivity, unique repair, standards conformance, or proof outside the admitted coefficient and cover class",
        },
        ClaimSpec {
            id: "i14-cover-refinement-naturality-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked cover-refinement theorem proves that admitted common refinements, nerve/cofinal comparison maps and canonical reaggregation preserve physical global sections, relative obstruction classes, port power and the passive-composition verdict while transporting every defect and ownership term exactly once. It also characterizes the hypotheses under which coarsening loses authority and returns an explicit comparison obstruction instead of silently changing the theorem",
            hypotheses: &[
                "both covers/sites, refinement functor, relative boundaries, orientation data, local complexes, comparison natural transformations, chain homotopies, coefficient/torsion policy and defect-owner maps are frozen",
                "a Leray, acyclicity, cofinality, quasi-isomorphism or exact alternative premise is named for each claimed invariant; arbitrary cover change receives no equivalence authority",
                "the comparison commutes with restriction/corestriction, total differential, trace, power pairing and external boundary supply, and maps shared defects without duplication or cancellation laundering",
                "the formal proposition, runtime comparison checker, axiom closure and refinement/coarsening negative twins satisfy the theorem policy",
            ],
            qoi: "cover_refinement_section_obstruction_power_and_passivity_naturality_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.cover_refinement.v1 composite proofs/i14/CoverRefinement.lean::passivityNaturality plus fs-vmanifest-oracles/i14/cover_refinement.rs::check_comparison_diagrams",
                independent: true,
                tcb_overlap: "shares canonical cover and comparison ASTs; pinned proof kernel and independent diagram/total-complex checker do not call production refinement or assembly code",
            },
            activation: "the refinement target card, both cover identities and exact comparison hypotheses are bound in a pre-candidate theorem successor",
            kill: "noncommuting comparison, changed physical section or obstruction class, duplicated/lost supply or defect, unjustified coarsening equivalence, theorem binding failure, or AuthorityContradiction blocks the revision",
            fallback: "retain cover-specific receipts with no cross-cover theorem authority and escalate to a common admitted refinement",
            no_claim: "cover refinement is not universally invariant; no equivalence is claimed without the exact comparison, topology, coefficient, boundary and regularity premises",
        },
        ClaimSpec {
            id: "i14-kyp-sheaf-passivity-bridge-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked KYP-to-sheaf bridge maps each admitted finite-dimensional rational or descriptor generalized-positive-real MTL, PEEC, circuit or reduced component to the exact storage, supply, dissipation, quotient, initialization and boundary-port relation consumed by passive sheaf-cosheaf composition. The bridge preserves simple imaginary-axis/zero/infinity storage poles, semidefinite storage, algebraic constraints and lossless modes without falsely upgrading generalized PR to strict PR or sampled-band dissipativity to a time-domain theorem",
            hypotheses: &[
                "the real or realified complex descriptor realization, regular pencil, impulse-free/consistent initialization class, controllability/observability or exact nonminimal substitute, quotient modes, port orientation and supply matrix are frozen",
                "generalized positive-realness means open-right-half-plane analyticity and PSD Hermitian part with separately admitted simple boundary-axis poles and PSD residues plus checked infinity term; strict PR uses its own shifted or coercive definition",
                "the KYP/LMI or storage witness is independently checked, transported through the exact port/reference crosswalk and bound to the local relation used by the composition theorem",
                "delay, fractional, convolutional, switching or infinite-dimensional components remain outside this bridge unless a successor proves a strength-matched extension",
            ],
            qoi: "generalized_pr_descriptor_kyp_storage_supply_and_local_relation_bridge_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.kyp_sheaf_bridge.v1 composite proofs/i14/KypSheafBridge.lean::generalizedPrToLocalDissipativity plus fs-vmanifest-oracles/i14/kyp_bridge.rs::check_descriptor_witness",
                independent: true,
                tcb_overlap: "shares canonical realization and rational witness bytes; pinned proof kernel and independent exact-rational descriptor/KYP checker do not call production fitter, MOR or composition code",
            },
            activation: "the KYP bridge card, exact realization class and theorem policy bind in a pre-candidate successor",
            kill: "unadmitted boundary pole, indefinite residue/storage, irregular or impulsive descriptor mode, inconsistent initialization, supply/port sign mismatch, sampled-only promotion, theorem binding failure, or AuthorityContradiction blocks the bridge",
            fallback: "retain sampled frequency evidence or a directly proved time-domain local supply inequality without KYP bridge authority",
            no_claim: "generalized PR or a feasible numerical LMI alone is not strict stability, a universal causal material theorem, global passivity, or physical validation",
        },
        ClaimSpec {
            id: "i14-maximal-counterexample-search",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: exhaust a future cardinality-proved finite rational microgrammar and search a separately seeded adversarial supergrammar for full-wave stability failures, BEM/FMM dense-envelope escapes, fidelity-bound violations, robust mitigations that violate the frozen safety envelope, assurance laundering, and countermodels to the fidelity-descent, passive-causal composition, hypercohomology-obstruction, cover-refinement-naturality and KYP-sheaf-bridge theorem revisions. Every candidate is independently classified, every nonvacuity/coverage floor is checked, and every genuine in-domain countermodel is minimized and retained",
            hypotheses: &[
                "the pre-search successor replaces target prose with complete canonical grammar ASTs, coefficient/degree/unit/reference/topology/signal-space semantics, validity and theorem-premise predicates, exact enumeration/exclusion order, rank/unrank/sharding algorithms, cardinality/bijection proof, cost preflight and Merkle completeness root",
                "the exhaustive microgrammar includes nontrivial passive rational MTL/PEEC/circuit components, at least one multiply connected relative-complex/winding case, one crosswalk/gluing defect, one well-posed feedback interconnection and declared negative twins; the larger grammar covers dispersive, delayed, nonnormal, topology/rank/event, PML, BEM/FMM and hostile uncertainty cases",
                "candidate states remain ClassificationPending until independent adjudication returns GenuineCountermodel, OutOfDomain, SpecificationDefect, AdmissionCheckerDefect, ImplementationCheckerDefect, ProofKernelOrTcbDefect, AuthorityContradiction, or reasoned RejectedOrIndeterminate. A GenuineCountermodel admitted against every immutable premise refutes an unproved universal theorem revision; if the exact proposition/axiom digest already has an admitted kernel proof, the impossible pair becomes AuthorityContradiction and quarantines both authorities pending deterministic binding/premise/arithmetic/oracle/axiom/kernel audit",
                "finite falsifier survival never promotes mathematics to Proved, physical validation, standards conformance or safety; formal kernel acceptance and runtime-premise admission remain separately mandatory",
            ],
            qoi: "exact_nonvacuity_coverage_zero_genuine_countermodel_and_integrity_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i14.oracle.maximal_falsifier.v1 at fs-vmanifest-oracles/i14/maximal_falsifier.rs::verify_enumeration_membership_and_minimize",
                independent: true,
                tcb_overlap: "shares canonical candidate and declaration bytes only; independent decoder/enumerator, premise/admission checker, formal kernels, high-precision evaluator, dense Maxwell/BEM paths and minimizer are separately version-pinned",
            },
            activation: "maximal declaration identities are frozen and a pre-candidate successor discharges the complete executable grammar and formal-projection gates in i14-theorem-falsifier-grammar and i14-theorem-formalization-policy",
            kill: "the first independently admitted in-domain countermodel refutes its exact unproved claim revision; an exact-digest proof/countermodel pair is AuthorityContradiction, not simultaneous proof and refutation, and blocks all affected authority while retaining both artifacts and auditing AST binding, assumptions, admission, arithmetic, oracle independence, axiom closure and proof kernel. Grammar escape, rank/unrank/canonicalization/cardinality defect, missed nonvacuity floor, incomplete budget, correlated checker, or artifact-integrity failure fails the campaign rather than passing it",
            fallback: "narrow or replace the exact claim through an authenticated amendment while retaining every candidate, failure and countermodel",
            no_claim: "version-1 prose has no exhaustive-search authority; a bounded empty search is not theorem proof, physical validation, safety evidence, or regulatory compliance",
        },
    ]
}

const fn authored_fixture(
    id: &'static str,
    spec: &'static str,
    partition: Partition,
) -> FixturePin {
    FixturePin {
        id,
        source: FixtureSource::AuthoredSpec { spec },
        partition,
    }
}

#[allow(clippy::too_many_lines)]
fn i14_fixtures() -> Vec<FixturePin> {
    vec![
        authored_fixture(
            CAMPAIGN_POLICY_FIXTURE,
            r#"POLICY: EVIDENCE_AUTHORITY_CLASS={PublicReplayCore,PublicReplayMax,GovernedBlindSlot,GovernedPhysicalSlot,GovernedPopulationSlot,GovernedStandards}. PublicReplay development and held-out AuthoredSpec generators may promote only deterministic schema, algebra, numerical-conformance and replay mechanics under their exact public grammar; they never imply IID sampling, untouched-data status, physical validation, standards conformance or population validity. GovernedStandards is an access-controlled authorized-input class, not a statistical holdout class. Partition and FixtureSource syntax alone carry no custody, discharge or promotion authority.
GOVERNED_EXTERNAL_BINDING: every governed external waiver remains NoPromotionAuthority until fs-vvreg verifies a typed DischargeReceipt binding the exact waiver subject and predicate, artifact kind, acquisition/generator/procedure/edition/version and validity rules, protected artifact roots, custodian signatures, access-ledger prefix, candidate/model/toolchain/checker/AcceptanceCard roots as applicable, claim/leaf scope, predecessor manifest digest, successor version and domain-separated transaction-intent digest. Local manifest structure binds receipt identities but does not authenticate issuer, capability, signature, trust policy or revocation state; HELM/ledger performs that authority check.
TRANSACTION_INTENT_V1: transaction_intent=BLAKE3::derive_key('org.frankensim.i14.governed-transaction-intent.v1',P), where P is the canonical length-framed successor-intent projection. P begins with the 25 exact ASCII bytes I14_TRANSACTION_INTENT_V1 followed by one byte 0x00, then U16LE(23), then exactly twenty-three U16LE(tag)||U64LE(payload_byte_len)||payload fields in strictly increasing tag order: 0x0001 initiative_id Utf8; 0x0002 schema_identity Utf8; 0x0003 operation_kind U16LE where BlindPopulationRealization=1,PhysicalValidationRealization=2,StandardsAuthorityCommit=3,CoupledExternalRealization=4; 0x0004 successor_version U64LE; 0x0005 predecessor_manifest_digest Digest; 0x0006 expected_authority_head AuthorityHead; 0x0007 predecessor_stage_receipt_digest Digest; 0x0008 candidate_freeze_commitment_digest Digest; 0x0009 mutation_fence MutationFence; 0x000a governance_stage U16LE where RealizationCommitted=1,StandardsAuthorityCommitted=2; 0x000b authority_scope Utf8; 0x000c coupled_transaction_group Utf8; 0x000d retired_waiver_subjects Utf8List; 0x000e protected_bindings ProtectedBindingList; 0x000f governed_slot_replacements FutureArtifactList; 0x0010 discharge_envelope_slots FutureArtifactList; 0x0011 final_successor_slot FutureDigest; 0x0012 amendment_record_slot FutureDigest; 0x0013 discharge_receipt_schema_digest Digest; 0x0014 authority_commit_receipt_schema_digest Digest; 0x0015 governance_protocol_schema_digest Digest; 0x0016 access_ledger_prefix_digest Digest; 0x0017 redaction_policy_digest Digest. Cross-field validity requires initiative_id='I14', schema_identity='i14-governed-transaction-intent-v1', checked successor_version=predecessor_version+1, one nonempty coupled_transaction_group, the exact operation-specific authority scope committed by governance_protocol_schema_digest, BlindPopulationRealization/PhysicalValidationRealization/CoupledExternalRealization with RealizationCommitted, StandardsAuthorityCommit with StandardsAuthorityCommitted, a nonempty exact predecessor-minus-successor retired-waiver set, and one MutationFence binding nonzero Digest(idempotency_key)||U64LE(attempt_id)||U64LE(capability_epoch). The retired-waiver subjects are exactly the waiver subjects represented by protected_bindings and biject the singleton related-waiver sets of same-ID discharge-envelope slots. Each governed-slot replacement's related-waiver set is exactly the sorted set of protected bindings naming its target slot/role/schema, which permits an indivisible joined output to discharge multiple coupled waivers without duplicating or role-swapping the output. StandardsAuthorityCommit alone requires an empty governed-slot-replacement list; every other operation requires the exact nonempty role-addressed image committed by the governance protocol. Every output, final-successor and AmendmentRecord slot id is pairwise distinct, and every future digest union is Pending. Coupled operations bind every component in the one group and forbid partial membership.
TRANSACTION_INTENT_ENCODING_V1: Digest is exactly 32 nonzero raw bytes. Utf8 is U64LE(byte_len)||exact nonempty case-sensitive bytes without normalization and has at most 65536 payload bytes. AuthorityHead is Digest||U64LE(generation), exactly 40 bytes. MutationFence is Digest(idempotency_key)||U64LE(attempt_id)||U64LE(capability_epoch), exactly 48 bytes. FRAME_BYTES(x)=U64LE(byte_len)||x. Utf8List=U64LE(count)||concat(FRAME_BYTES(Utf8)); it is lexicographically ordered by each complete framed element and duplicate-free. ProtectedBinding=FRAME_BYTES(Utf8(waiver_subject))||FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||Digest(protected_root), and ProtectedBindingList=U64LE(count)||concat(FRAME_BYTES(ProtectedBinding)); the list is ordered by complete encoded record, duplicate-free, and unique by waiver-subject/slot/role relationship. FutureDigest=FRAME_BYTES(Utf8(slot_id))||digest_union. FutureArtifact=FRAME_BYTES(Utf8(target_slot_id))||FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||Utf8List(related_waiver_subjects)||digest_union, and FutureArtifactList=U64LE(count)||concat(FRAME_BYTES(FutureArtifact)); the list is ordered by complete encoded record, duplicate-free, and unique by target-slot/role pair. digest_union is exactly byte 0=Pending with no following bytes or byte 1=Digest followed by 32 raw bytes; tags 0x000f..0x0012 require Pending. Pending envelope slots carry the I14_ENVELOPE_DIGEST_PENDING_V1 semantics and the pending AmendmentRecord/final-successor slots carry I14_AMENDMENT_RECORD_PENDING_V1 semantics without serializing a fabricated digest. Every scalar tag occurs exactly once and every fixed-width payload has exactly the stated length. Every compound list contains at most 4096 records, P is at most 16777216 bytes, every count/length/allocation uses checked arithmetic against remaining input before allocation, and missing, extra, duplicate, out-of-order, unknown, nonminimal, over-cap or trailing bytes are refused. Tag 0x0007 binds the already-existing predecessor-stage receipt; P never serializes a newly created discharge-receipt digest, authority-commit-receipt digest, realized envelope digest, final successor digest, realized AmendmentRecord digest or new authority head.
TRANSACTION_OPERATION_MATRIX_V1: BlindPopulationRealization requires governance_stage=RealizationCommitted, authority_scope='i14.blind-population-realization', at least one governed-slot replacement and at least one retired waiver. PhysicalValidationRealization requires governance_stage=RealizationCommitted, authority_scope='i14.physical-validation-realization', at least one governed-slot replacement and at least one retired waiver. StandardsAuthorityCommit requires governance_stage=StandardsAuthorityCommitted, authority_scope='i14.standards-authority-commit', exactly zero governed-slot replacements and at least one retired waiver. CoupledExternalRealization requires governance_stage=RealizationCommitted, authority_scope='i14.coupled-external-realization', at least one governed-slot replacement and at least two retired waivers. In every row discharge-envelope count equals retired-waiver count, protected bindings cover every and only retired waiver under the role-addressed relationships above, and coupled_transaction_group byte-equals the nonempty group identity authenticated by predecessor_stage_receipt_digest. No other operation/stage/scope/cardinality tuple is valid.
TRANSACTION_SCHEMA_REGISTRY_V1: SchemaDigestV1(role,schema_bytes)=BLAKE3::derive_key('org.frankensim.i14.governed-schema-artifact.v1',FRAME_BYTES(Utf8(role))||FRAME_BYTES(schema_bytes)). The three exact case-sensitive roles are DischargeReceipt, AuthorityCommitReceipt and GovernanceProtocol, and tags 0x0013,0x0014,0x0015 respectively must equal their SchemaDigestV1 values. GovernanceSchemaSetRootV1=BLAKE3::derive_key('org.frankensim.i14.governance-schema-set.v1',U64LE(3)||FRAME_BYTES(Utf8('DischargeReceipt')||raw32(discharge_receipt_schema_digest))||FRAME_BYTES(Utf8('AuthorityCommitReceipt')||raw32(authority_commit_receipt_schema_digest))||FRAME_BYTES(Utf8('GovernanceProtocol')||raw32(governance_protocol_schema_digest))) in exactly that role order. The authenticated predecessor-stage receipt binds this exact root, each role/digest, the operation-matrix row, group identity and a closed role-addressed artifact-schema table; transaction-time callers cannot substitute another schema or scope. ArtifactSchemaDigestV1(artifact_role,schema_bytes)=BLAKE3::derive_key('org.frankensim.i14.governed-artifact-schema.v1',FRAME_BYTES(Utf8(artifact_role))||FRAME_BYTES(schema_bytes)); every ProtectedBinding/FutureArtifact artifact_schema must equal the table's exact ArtifactSchemaDigestV1 for that role. Every fetched schema is at most 1048576 bytes, has nesting depth at most 64 and at most 4096 field/constraint records; every DischargeReceipt is at most 1048576 bytes with nesting depth at most 64. Counts, lengths and remaining input are checked before allocation, and cap+1, depth+1, truncation, trailing bytes, role alias and schema swap are IntegrityFailed before mutation.
TRANSACTION_DECODER_DEPTH_V1: The top-level schema or receipt value begins at nesting depth 0. Entering the payload of each framed nested record, list, union or constraint-AST node increments depth exactly once; fixed-width scalars and the outer tag/length header do not increment it. Depth 64 is admitted, any attempted entry to depth 65 is refused before allocating or decoding that child, and sibling traversal never accumulates depth.
TRANSACTION_SCHEMA_ROLE_CONFORMANCE_V1: The decoded DischargeReceipt schema must require every receipt binding named by TRANSACTION_INTENT_AUTHORITY_V1 and its exact role/digest/GovernanceSchemaSetRootV1 membership fields; the decoded AuthorityCommitReceipt schema must describe exactly the closed 361-byte AuthorityCommitReceiptBytesV1 layout with no optional or extension field; and the decoded GovernanceProtocol schema must describe exactly TRANSACTION_OPERATION_MATRIX_V1, the relationship/bijection rules and the closed artifact-role/schema table. The AuthorityCommitReceipt proves schema membership by binding transaction_intent; verification resolves the exact P, recomputes all three schema digests and GovernanceSchemaSetRootV1 and decodes the receipt before accepting its digest. A schema that omits, relaxes, renames, retypes or extends any required field/constraint is a schema mismatch even when its bytes hash correctly.
TRANSACTION_INTENT_AUTHORITY_V1: The exact bytes of every content-addressed transaction, receipt and artifact schema must be available, rehashed under its exact role-bound domain and used to decode the corresponding receipt, governed-slot wrapper, discharge-envelope wrapper and realized output before mutation. A verified DischargeReceipt binds its exact schema role/digest and GovernanceSchemaSetRootV1 membership, transaction_intent, predecessor manifest/head/stage, its exact protected binding and schema commitments, capability epoch and idempotency fence, but never the final successor digest. discharge_receipt_set_root=BLAKE3::derive_key('org.frankensim.i14.discharge-receipt-set.v1',U64LE(receipt_count)||concat(FRAME_BYTES(raw32(receipt_digest)) in lexicographic raw-digest order)); receipt_count equals the retired-waiver count, is at most 4096, and duplicate, missing, extra, wrong-schema or relationship-mismatched receipts are IntegrityFailed. RealizedArtifactBytesV1=FRAME_BYTES(Utf8(target_slot_id))||FRAME_BYTES(Utf8(artifact_role))||raw32(artifact_schema_digest)||Utf8List(related_waiver_subjects)||raw32(realized_artifact_digest), and realized_output_set_root=BLAKE3::derive_key('org.frankensim.i14.realized-output-set.v1',U64LE(output_count)||concat(FRAME_BYTES(RealizedArtifactBytesV1) in lexicographic full-record order)); it is the exact realized image of every Pending governed-slot and discharge-envelope record and no other output. The complete successor embeds each realized output at its bound target slot; it embeds neither the AmendmentRecord nor the new authority head. The AmendmentRecord binds predecessor and final successor digests. successor_authority_generation=checked_add(expected_authority_generation,1) and successor_authority_head_digest=BLAKE3::derive_key('org.frankensim.i14.authority-head.v1',raw32(expected_authority_head_digest)||U64LE(expected_authority_generation)||raw32(transaction_intent)||raw32(discharge_receipt_set_root)||raw32(realized_output_set_root)||raw32(final_successor_digest)||raw32(amendment_record_digest)||U64LE(successor_version)). The new head and commit receipt are external authority-ledger outputs and are not embedded in the successor. AuthorityCommitReceiptBytesV1 is the 31 exact ASCII bytes I14_AUTHORITY_COMMIT_RECEIPT_V1 plus NUL, U64LE(attempt_id), U64LE(capability_epoch), U64LE(successor_version), raw32(idempotency_key_digest), raw32(transaction_intent), raw32(discharge_receipt_set_root), raw32(realized_output_set_root), raw32(final_successor_digest), raw32(amendment_record_digest), raw32(expected_authority_head_digest), U64LE(expected_authority_generation), raw32(successor_authority_head_digest), U64LE(successor_authority_generation), one byte cas_result where Committed=1, and raw32(authority_capability_digest), exactly 361 bytes with every digest nonzero; authority_commit_receipt_digest=BLAKE3::derive_key('org.frankensim.i14.authority-commit-receipt.v1',AuthorityCommitReceiptBytesV1). Every repeated receipt value must byte-equal its P/MutationFence value or the independently recomputed root, digest or head, and authority_capability_digest must byte-equal the authenticated capability that performs the commit; any per-field mismatch is IntegrityFailed. One indivisible authority transaction requires the current scope head to byte-equal expected_authority_head and the current capability epoch to equal the bound epoch; validates every schema, receipt, list relationship, root, capability and access-ledger/redaction commitment; installs every realized slot/envelope; removes exactly the retired waivers; freezes the complete successor; verifies the AmendmentRecord; compare-and-swaps the exact expected head/generation to the exact successor head/generation uniquely derived above; durably records the idempotency-key-to-intent-and-receipt mapping; and commits the receipt. A stale head, caller-proposed alternative head, generation overflow, ABA capability epoch, count/root/bijection/schema mismatch, split successor or partial durable state is IntegrityFailed with no authority advance or other effect. A crash before the single durable commit point leaves no effect; a crash after it preserves the exact committed head, successor, waiver/output changes, idempotency mapping and receipt, and recovery or exact replay returns that byte-identical receipt without a second effect. Exact-key byte-identical replay performs no second mutation or capability consumption; a conflicting payload under one key is IntegrityFailed. To avoid a content-hash cycle, the authority DAG is P -> DischargeReceipts -> realized outputs -> final successor -> AmendmentRecord -> successor head/AuthorityCommitReceipt. Before any discharge or promotion authority, two independently implemented encoders must reproduce published exact P, schema-set, receipt-set, output-set, successor-head and commit-receipt KATs plus header/NUL/tag/order/type/length/cardinality/cap/depth, role-swap, schema-swap, per-field receipt mismatch, output-add/drop/swap, stale-head, alternative-head, generation-overflow, ABA, replay, crash-before-commit, crash-after-commit and partial-state mutation twins. Version 1 does not claim those independent encoders already exist: until that proof gate and HELM/ledger authentication pass, this transaction path remains NoPromotionAuthority. A raw digest, receipt alone, same-ID fixture alone, slot replacement alone, waiver removal alone, split successors, or a structurally frozen successor lacking the verified fs-vvreg transaction never proves discharge.
TRANSACTION_PARTIAL_STATE_RECOVERY_V1: In TRANSACTION_INTENT_AUTHORITY_V1, rejecting a proposed partial mutation means validation aborts before the durable commit point with no effect. If recovery ever observes only a strict subset of the promised head/successor/waiver/output/idempotency/receipt state, that observation is itself durable corruption: authority and promotion fail closed, the state is quarantined for authenticated reconciliation, and it must never be reported as a successful zero-effect refusal or silently completed from unauthenticated evidence.
GOVERNED_BLIND_POPULATION_PATH: this path applies only to GovernedBlindSlot and GovernedPopulationSlot. GovernanceCommitted freezes generator/acquisition procedure, validity, applicability and exclusion rules, custodian, membership/case-order policy, thresholds, AcceptanceCards, checker, candidate-input permissions and access protocol before protected outcome, label, aggregate or adjudication access. CandidateFrozen binds the immutable candidate/model/toolchain and every decision rule before that protected access. The independent custodian then realizes or selects the inaccessible cases. One RealizationCommitted authority transaction must both replace every governed AuthoredSpec commitment slot consumed by the leaf with its typed content-addressed External realized root and install every waiver subject's same-ID typed External discharge-envelope root while retiring every corresponding waiver through its distinct verified DischargeReceipt. RevealedForAdjudication begins only after that transaction commits; Closed records every access, exclusion, retry, optional-stopping and multiplicity decision, missingness event, attempted candidate and immutable FailureBundle. Skipping, reordering or splitting a stage, partial multi-pack discharge, arbitrary-digest substitution, slot-only replacement, same-ID-envelope-only replacement, post-reveal amendment, or waiver-only retirement is IntegrityFailed.
GOVERNED_PHYSICAL_VALIDATION_PATH: GovernedPhysicalSlot separates protected calibration/model-input strata from an untouched validation stratum. Before either candidate-accessible stratum is read, an independent custodian/fs-vvreg capability creates salted or equivalently hiding content/Merkle commitments for calibration and as-built model-input, a hiding validation source-universe/frame commitment, a disjoint-membership commitment, the exact validation-selection algorithm, and a pre-candidate secret-seed/VRF or equivalently non-adaptive selection commitment. GovernanceCommitted freezes those opaque identities, stratum membership/selection rules, calibration and as-built model-instantiation procedure, permitted parameter/model updates, contamination predicates, custodian, AcceptanceCards, checker and every validation decision. Candidate builders, fitters, checker/threshold owners and their transitive capabilities receive only opaque commitment identities; before CandidateFrozen they receive no validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material. AuthorizedCalibration grants named builders least-privilege access only to the committed calibration and as-built model-input strata under a complete read/output/derivation ledger. AsBuiltModelInstantiation may fit parameters, select among preregistered models and construct the as-built simulation only under that frozen procedure. CandidateFrozen then binds the resulting model, parameters, toolchain, QoIs, thresholds and decision rules; an independent contamination receipt names and audits every candidate-side principal and transitive capability and proves the forbidden validation information was not accessed. The custodian independently realizes or selects the still-inaccessible validation stratum through the committed mechanism. One RealizationCommitted authority transaction replaces the physical AuthoredSpec commitment slot with a joined typed External root whose separately addressable calibration, model-input and validation roots carry membership proofs to the pre-access commitments, a mutual disjoint-membership proof, and a non-adaptive selection proof; installs every coupled waiver subject's same-ID discharge envelope; removes all coupled Waiver rows through distinct receipts; verifies the AmendmentRecord; and advances the authority head. RevealedForAdjudication exposes only the frozen validation stratum after that transaction commits. Closed retains calibration/model lineage, contamination audit, every access/exclusion/retry and immutable FailureBundle. Candidate-side validation access or opening before freeze, calibration-to-validation leakage, hidden model update, adaptive validation selection, partial multi-pack discharge, stratum aliasing, split retirement, post-reveal amendment or waiver-only retirement is IntegrityFailed.
GOVERNED_STANDARDS_PATH: GovernedStandards has no HeldOut AuthoredSpec commitment slot in manifest version 1. GovernanceCommitted freezes exact publisher, edition, corrigenda, license and source-root identities; required clause-inventory and applicability/exclusion grammar; authorized principals and sandbox; derivation/crosswalk procedure; checker, AcceptanceCards, disclosure filter and access protocol before licensed-input access. If an authenticated legally distributable clause graph is insufficient, AuthorizedConstruction may grant least-privilege licensed-byte access only to named builders under a complete read/output ledger and dual control. CandidateFrozen occurs after that authorized construction and binds the derived content-addressed crosswalk, toolchain, every mapping/exclusion/loss, restricted-text-free public projection and decision rule before independent adjudication against the same authenticated authority. One StandardsAuthorityCommitted FrozenManifest::amend transaction installs the waiver subject's same-ID typed External discharge-envelope root, removes the Waiver row through its verified DischargeReceipt, verifies the AmendmentRecord and atomically advances the authority head before adjudication. StandardsAdjudicated and Closed retain independent reconstruction, disagreement, access, disclosure and FailureBundle receipts. This path grants exact-edition scoped crosswalk authority only, never blind, untouched-data, IID, legal, certification or regulatory authority. Unauthorized access, edition shopping, omitted in-scope clauses, unledgered derived output, restricted-text leakage, split retirement, cross-root result reuse or post-adjudication change is IntegrityFailed.
HOLDOUT_REALIZATION: version-1 PublicReplayCore/PublicReplayMax HeldOut AuthoredSpecs are fully identified deterministic generator contracts and need no secret replacement for their narrowly scoped mechanics authority. Version-1 GovernedBlindSlot/GovernedPopulationSlot HeldOut AuthoredSpecs are commitment schemas only and remain NoPromotionAuthority until their blind/population slot replacement, every same-ID discharge-envelope binding and every waiver retirement close atomically. GovernedPhysicalSlot follows the calibration/model-instantiation/untouched-validation path and remains NoPromotionAuthority until its joined stratum root, contamination receipt, coupled envelopes and waivers close atomically. GovernedStandards has no version-1 AuthoredSpec slot and follows the authorized-input construction path instead.
JOINT_EXTERNAL_TRANSACTION: if one leaf consumes multiple protected strata or external packs, one atomic signed transaction binds the complete candidate, models, checkers, thresholds, AcceptanceCards, every slot root, every same-ID discharge envelope and every DischargeReceipt before protected adjudication access. Sequential adaptation, selective opening, partial retirement, retries across strata, threshold edits after partial reveal or revealing one component while another remains waived are IntegrityFailed.
SEED DOMAINS: development 0..=4095; Core PublicReplay holdout 65536..=69631; Max PublicReplay holdout 131072..=135167; falsifier 196608..=212991. These ranges are disjoint inclusive intervals under the Five Explicits derivation and confer replay only."#,
            Partition::Development,
        ),
        authored_fixture(
            ACCEPTANCE_POLICY_FIXTURE,
            r#"ACCEPTANCE_ARITHMETIC_V2: manifest version 1 fixes the schema and refuses numeric promotion until a pre-candidate FrozenManifest::amend successor supplies one content-addressed AcceptanceCard per numeric claim and per independently scored component. Each card freezes physical unit, comparator class, exact reference artifact, certified point or enclosure representation, field/port/support, norm, quadrature or sampling measure, complex-value convention, aggregation order, applicability predicate, missing/nonfinite behavior, and every uncertainty/coverage rule. HARD COMPARATORS have no tolerance laundering: for a certified result enclosure X=[x_lo,x_hi], HardUpper(u) is Supported iff x_hi<=u; HardLower(l) iff x_lo>=l; HardInterval(l,u) iff finite dimension-compatible l<=u and l<=x_lo and x_hi<=u. Their diagnostic violation magnitudes are h_upper=max(0,x_hi-u), h_lower=max(0,l-x_lo), and h_interval=max(0,l-x_lo,x_hi-u); Supported iff the applicable h=0 exactly under outward rounding, and h divided by any scale is never adjudicated by <=1. Any overlap that fails to establish containment is Unknown or Failed under the frozen one-sided policy, never a normalized pass. SOFT COMPARATORS explicitly freeze absolute tolerance a>0, relative coefficient rho>=0, candidate-independent reference scale r>=0, and tau=max(a,rho*r), with a and r in the physical unit and rho dimensionless. SoftEquality uses d_eq=sup_{x in X,y in Y} norm(x-y); SoftUpper uses d_upper=max(0,x_hi-u); SoftLower uses d_lower=max(0,l-x_lo); SoftInterval uses d_interval=max(0,l-x_lo,x_hi-u). For the explicitly authorized soft comparator e=d/tau and Supported iff its outward upper bound e_hi<=1. Every hard predicate and every preregistered applicable soft component are conjoined; boundary equality is inclusive. A claim score is max_i e_i only for the soft-score vector after every component is admitted and never masks a hard failure. The card freezes IEEE-754 format, round-to-nearest-ties-to-even point arithmetic where used, operation/reduction order, FMA contraction policy, signed-zero/canonical-NaN policy, subnormal/overflow behavior, and directed outward rounding for interval and algebraic-constant evaluation. Exact algebraic constants use minimal polynomial plus rational isolating interval; rational or IEEE bit patterns remain allowed. Measurement uncertainty, covariance, equivalence tests, confidence/e-process coverage and censoring are separately typed comparators whose confidence/coverage level and uncertainty budget freeze before evidence and cannot be inflated to create acceptance. Missing, undefined, nonfinite, invalid enclosure, zero/negative tau input, unit mismatch, unbound denominator, post-hoc component selection, changed aggregation, or applicability ambiguity is IntegrityFailed or Unknown and cannot promote. Exact-bit claims remain exact and cannot be converted into a forgiving numeric score. Every row consumes this policy, and the adjudication receipt binds every AcceptanceCard digest, hard-comparator result and full soft component/enclosure vector rather than only a maximum."#,
            Partition::Development,
        ),
        authored_fixture(
            EM_CONVENTION_FIXTURE,
            r#"EM_CONVENTION_CARD_V1: alpha is the unique positive algebraic root alpha^2=2, encoded by minimal polynomial x^2-2 plus a rational isolating interval; x(t)=Re{alpha X_+ exp(+i omega t)}. A global time-origin/reference-phase shift transforms every coupled source, field, port, probe and material response together and preserves relative phase; relative source/probe phase is physical frozen data. The opposite exp(-i omega t) presentation uses X_-=conj(X_+) for every coupled quantity at once, never partial conjugation. Under +i omega t, curl E=-i omega B, curl H=J_impressed+(sigma+i omega epsilon)E, RMS terminal complex power is V I*, port current is positive into the modeled subsystem, and outward complex Poynting flux is integral(E cross H*) dot n. The passive outgoing scalar branch is G_+(r)=exp(-i k r)/(4 pi r) with Re(k)>=0 and Im(k)<=0 and the corresponding Sommerfeld sign. Passive local conductivity has PSD Hermitian part; for an admitted local frequency-domain dielectric/magnetic representation at omega>0 the Hermitian imaginary parts of epsilon and mu are negative semidefinite under this convention, or the material consumes a stronger causal passive realization. PML stretch signs, impedance conditions, BEM kernels and sided traces, radiation ownership, source/storage/dissipation balance, and complex-adjoint conjugate transpose all derive from this card. RMS/peak conversion transports alpha symbolically; every floating adapter uses a certified outward enclosure and never a decimal approximation of alpha or a false bit-exact amplitude round trip."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-harnessgraph-canonical",
            r#"GENERATOR: typed HarnessGraph assemblies with occurrence-stable connectors, contacts, conductors, shields, drains, splices, backshells, chassis/bond nodes, reference conductors, sensor/actuator/victim pins, geometric route segments, material/temperature states, port orientations, and explicit open/unmated/short/fault states. Include repeated part numbers, permutation-isomorphic assemblies, dangling shields, multi-ground loops, occurrence-preserving edits, and invalid alias/connectivity twins. ORACLE: independently canonicalize the typed incidence/occurrence relation and compare connectivity, orientation, units, and graph digest. SHRINK: remove leaves/segments while preserving the first identity or conservation fault. SEEDS: i14/i14-harnessgraph-canonical/topology, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-harnessgraph-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: the HarnessGraph grammar of i14-harnessgraph-canonical with unseen connector multiplicities, braided shield terminations, repeated subassemblies, floating references, and topology/fault twins. Sole consumer i14-harnessgraph-core. Public deterministic replay only, withheld from ordinary development execution until adjudication; it has no IID or untouched-data authority. SEEDS: i14/i14-harnessgraph-core-holdout/topology, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-ap242-adapter-mechanics",
            r#"DECK: legally authored AP242-shaped harness exchange records for the declared synthetic adapter subset. Cover occurrence reuse, nested assemblies, product-definition relationships, units, local/global frames, transforms, connector/pin occurrence identity, open/short/splice/shield/drain/bond semantics, route and material-property lineage, exact source-byte roots, unsupported constructs, ambiguity, duplicate identifiers, missing references, and one-way loss receipts. This deck tests adapter mechanics and embeds no licensed standard text. ORACLE: independently reconstruct occurrence and transform graphs and reconcile every source record to one target identity, explicit loss, quarantine, or refusal. SEEDS: i14/i14-ap242-adapter-mechanics/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-ap242-adapter-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen legally authored AP242-shaped adapter cases with repeated occurrences, mixed unit/frame transforms, nested reuse, loss-accounting traps, ambiguous semantics, unsupported records and refusal twins. Sole consumer i14-ap242-adapter-core. Public deterministic Core replay holdout only; it grants no licensed-edition conformance. SEEDS: i14/i14-ap242-adapter-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-mtl-coax-two-wire",
            r#"DECK: analytic and high-precision coaxial, twin-lead, two-wire-above-ground, and short multiconductor transmission-line cases. Bind reference conductor, conductor ordering, cross-section/material state, RMS exp(+i omega t) phasor convention, terminal orientation, and per-unit-length R,L,G,C units. Include lossless, skin/proximity-loss, dielectric-loss, dispersive rational, nearly degenerate modal, open/short/matched/mismatched, and invalid nonreciprocal/passivity twins. ORACLES: closed-form TEM subsets; independent high-precision matrix exponential/Riccati propagation; positive-real and Kramers-Kronig residual checks. SEEDS: i14/i14-mtl-coax-two-wire/sweep, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-rlgc-operator-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen ordered multiconductor RLGC/operator records with reference-quotient changes, quasi-static and dispersive classes, generalized-PR boundary poles, bad residues, active pockets between samples, nonminimal passive realizations, zero/infinity asymptotics, covariance/source traps and invalid interpolation twins. Sole consumer i14-rlgc-operator-core. Public deterministic Core replay only; no IID, population or universal causal authority. SEEDS: i14/i14-rlgc-operator-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-mtl-bundle-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: 3..16-conductor routed segments with admitted RLGC inputs, unseen conductor permutations, reference changes, nearly repeated modes, dispersive convolution/state initialization, connector discontinuities, open/short/matched loads, delayed reflections and power-balance negative twins. Sole consumer i14-mtl-propagation-core. Public deterministic Core replay; no IID or statistical-population authority. SEEDS: i14/i14-mtl-bundle-core-holdout/sweep, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-rlgc-fit-adversaries",
            r#"GENERATOR: frequency samples and rational/state-space RLGC candidates spanning strict positive-real/asymptotically stable, generalized positive-real lossless-boundary-pole, passive convolutional, weakly nonpassive, unstable-pole, non-Hermitian, bad-residue, zero/infinity-asymptotic, reference-singular, aliased, underresolved, and extrapolation-hostile cases. Preserve raw sample covariance and provenance. CONVENTION: X(s)=integral_0^infinity x(t) exp(-s t) dt; analytic positive-real domain Re(s)>0; frequency response is the non-tangential boundary value where it exists. Strict PR rational realizations have poles Re(s)<0. Generalized PR permits only admitted simple imaginary-axis poles, including zero, with PSD Hermitian residue matrices and separately checked infinity behavior. LAWS: conjugate realness, open-half-plane matrix Hermitian inequality, quotient regularity, dimensional/Hermitian structure; quasi-static L/C positivity only under its premises; sampled-band Hermitian dissipativity kept distinct from analytic positive-real Z(s)/Y(s); passive internal-state/convolution realization and its storage; stable causal realization; Kramers-Kronig/causality residual; passivity-preserving interpolation; refusal outside support or at an unhandled pole. SHRINK: delete frequencies/poles/states while preserving the first authority defect. SEEDS: i14/i14-rlgc-fit-adversaries/fit, development and falsifier ranges as declared by the consuming lane."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-peec-analytic",
            r#"DECK: partial-element equivalent-circuit cells with analytic self/mutual partial inductance and potential-coefficient P references, separately derived capacitance C after gauge/reference quotient, retardation/delay cases, conductor-loss owners, radiation ports, charge/current continuity, and Lorenz-gauge residuals. Include near-singular separation, disconnected/reference-rank defects, P-versus-C substitution traps, double-counted loss, nonpassive reduction, and invalid gauge twins. ORACLES: independent high-precision quadrature/dense solve plus exact small-network power audit. SEEDS: i14/i14-peec-analytic/cell, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-peec-mtl-crossover-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: routed geometries intentionally spanning electrically short, transition, and retarded regimes where MTL and PEEC validity overlap or separate. Bind the common physical ports, reference/gauge transforms, source normalization, loss owners, and independent fine PEEC/MTL or analytic baselines. Sole consumer i14-fidelity-routing-core. Public deterministic holdout, not an IID sample and not proof of universal cross-rung validity. SEEDS: i14/i14-peec-mtl-crossover-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-peec-network-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen quasistatic and retarded PEEC conductors, disconnected/reference quotients, near-singular cells, charge-continuity and Lorenz-gauge traps, passive/nonpassive reductions, delay/history cases and P-versus-C substitution twins. Sole consumer i14-peec-network. Public deterministic Core replay holdout only; no universal geometry, physical-validation or population authority. SEEDS: i14/i14-peec-network-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-ground-bond-shield",
            r#"DECK: DC-through-RF bond/ground/shield networks with finite contact impedance, surface-transfer impedance, braid coverage, pigtail/360-degree termination, aperture/backshell coupling, common-impedance coupling, chassis return, ground loops, and explicit shield internal/external surface currents. Include open, floating, corroded, nonlinear contact, duplicated-loss, omitted-return, and orientation-negative twins. ORACLES: exact small networks, analytic coax shield-transfer subsets, and independent charge/current/power closure. SEEDS: i14/i14-ground-bond-shield/network, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-shield-ground-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen shield/bond/backshell/chassis assemblies with multiple apertures, nested shields, asymmetric grounds, uncertain contact states, and coupled common/differential paths. Sole consumer i14-ground-shield-core. Public deterministic holdout with replay authority only. SEEDS: i14/i14-shield-ground-core-holdout/network, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-motor-bearing-current",
            r#"DECK: inverter-fed motor common-mode source, winding-frame capacitance, rotor/stator/frame/bearing paths, shaft voltage, bearing EDM threshold/event model, circulating current, grounding brush, cable shield, and frequency-dependent bearing impedance. Bind temperature/lubricant/material state and uniquely own dielectric/contact/arc losses. Include open-ground, insulated-bearing, raceway-discharge, threshold-grazing, event-order, and invalid energy-creation twins. ORACLES: exact reduced circuits, independent event integration, and charge/power budgets. SEEDS: i14/i14-motor-bearing-current/transient, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-bearing-current-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen motor-drive common-mode networks, winding/frame/shaft/bearing capacitance ratios, ground-brush and cable-shield variants, lubricant-film histories, threshold-grazing breakdown/recovery, simultaneous contacts, chatter, censoring and energy-creation twins. Sole consumer i14-bearing-current-core. Public deterministic Core replay holdout for path/event semantics only; no bearing-life or population authority. SEEDS: i14/i14-bearing-current-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-switching-source-probe",
            r#"DECK: PWM/edge source models with bus impedance, dead time, common/differential-mode decomposition, jitter and spectrum windows; typed voltage/current/field probes with bandwidth, transfer function, loading, placement/orientation, calibration uncertainty, anti-alias filtering, window/normalization convention, and RMS/peak adapters. Include ideal-step bandwidth escape, probe saturation/loading, missing calibration, alias/leakage, sign/conjugation, and double-normalization twins. ORACLES: analytic pulses/spectra and independent sampled-instrument pipeline. SEEDS: i14/i14-switching-source-probe/trace, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-source-probe-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen switching patterns, cable terminations, instrument transfer functions, placements, jitter/window phases, bandwidth edges, and saturation cases. Sole consumer i14-source-probe-core. Public deterministic replay holdout; it is not blind laboratory evidence. SEEDS: i14/i14-source-probe-core-holdout/trace, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-immunity-victim",
            r#"DECK: victim transfer and upset-mode ledgers connecting typed terminal/field exposure to susceptibility mode, state, threshold distribution, dwell/debounce, recovery, hysteresis, uncertainty, and outcome. Include linear susceptibility, nonlinear latch, metastable/event-order, multiple-cause, censored/no-upset, ambiguous-mode, and threshold-unit negative twins. Distinguish simulated exposure, component characterization, laboratory observation, and standards evidence. ORACLES: analytic reduced victims and independent event/mode classifier. SEEDS: i14/i14-immunity-victim/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-victim-upset-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen victim modes, correlated thresholds, temporal windows, recovery paths, multi-source attribution, and adversarial near-boundary outcomes. Sole consumer i14-victim-mode-core. Public deterministic holdout only; statistical reliability requires a separate GovernedPopulation campaign. SEEDS: i14/i14-victim-upset-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-core-crosswalk-ladder",
            r#"DECK: HarnessGraph-to-MTL, MTL-to-PEEC, PEEC-to-reduced-circuit, source/probe, and victim crosswalks over declared fidelity regions. Each case binds common physical ports, orientations, gauge/reference transforms, phasor/time convention, state/material snapshot, loss ownership, projection/lift maps, conditioning, residual/error budget, and valid QoIs. Include overlapping and disjoint regions, hysteresis/chatter boundaries, topology/rank changes, lossy noncommuting paths, and misleading agreement twins. ORACLES: independent fine-rung or analytic baselines plus triangle/diagram residual audit. SEEDS: i14/i14-core-crosswalk-ladder/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-adjoint-uq-mitigation",
            r#"DECK: fixed-regime complex-valued adjoint cases and uncertainty/mitigation candidates across MTL, PEEC, grounding, shielding, probe, and victim models. Bind Wirtinger/real objective convention, conjugate transpose, parameter units/transforms, frozen topology/event/routing regime, covariance/dependence model, aleatory/epistemic labels, censoring/missingness, multiplicity/stopping policy, fabrication/cost/thermal constraints, and candidate history. ORACLES: complex-step or finite-difference directional derivatives away from events, independent sampling/bounds, and exhaustive small design grids. SEEDS: i14/i14-adjoint-uq-mitigation/study, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-adjoint-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen fixed-rung complex adjoint cases across MTL, PEEC, grounding, shielding, source/probe and smooth victim regimes, with parameter-unit transforms, nonnormal systems, near conditioning limits, inexact solves, event/topology boundary refusals and intentionally wrong conjugate/trace derivatives. Sole consumer i14-fixed-regime-adjoint-core. Public deterministic Core replay holdout. SEEDS: i14/i14-adjoint-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-uq-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen aleatory/epistemic/dependence/censoring/missingness and optional-stopping cases, including correlated tails, rare-event stress, public-deterministic no-population traps, miscoverage and certify-or-escalate boundaries. Sole consumer i14-uq-inference-core. Public deterministic Core replay holdout for inference mechanics only; no physical-population authority. SEEDS: i14/i14-uq-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-mitigation-max-holdout",
            r#"HOLDOUT_KIND: GovernedBlindSlot. COMMITMENT SLOT: robust-mitigation cases span routed geometry, termination, filtering, grounding, shielding, source spectrum, victim state, material/contact uncertainty and correlated corner populations. No public seed or authored bytes confer blind authority. GovernanceCommitted freezes the generator/procedure, validity/exclusion rules, custodian, case-count/order policy, checker, thresholds, AcceptanceCards and access protocol. CandidateFrozen binds the immutable candidate/model/toolchain without protected holdout root, byte, label, aggregate or derived-statistic access; public development fixtures and other explicitly authorized candidate inputs remain available under the frozen candidate-input permissions. The independent custodian then realizes the inaccessible population and one atomic RealizationCommitted FrozenManifest::amend authority transaction replaces this slot with its typed External content/Merkle root, installs i14-external-blind-mitigation-custody-pack as its same-ID typed External discharge-envelope root, removes that Waiver row through a verified fs-vvreg DischargeReceipt, verifies the AmendmentRecord and advances the authority head. Only then may one-shot RevealedForAdjudication begin. Sole consumer i14-robust-mitigation-max. Slot-only replacement, envelope-only replacement, waiver-only retirement, arbitrary digest substitution, split successors, sequential reveal or any post-reveal successor is IntegrityFailed. This authored slot provides schema validation only until the complete authority transaction closes."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-fullwave-cavity-waveguide",
            r#"DECK: manufactured and canonical Maxwell cavity, waveguide, coax discontinuity, resonator, multiply connected domain, dispersive medium, conductive wall, port, and radiation cases. Bind relative/absolute boundary conditions, de Rham complex, orientations, gauge/nullspaces, source compatibility, RMS exp(+i omega t) convention, complex power, material passivity/causality, mesh/geometry error, and exact eigen/manufactured fields. Include spurious-mode, topology-rank, sign/conjugation, anisotropy, high-contrast, and underresolved twins. SEEDS: i14/i14-fullwave-cavity-waveguide/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fullwave-schema-crosswalk",
            r#"DECK: typed FullWaveProblem/Receipt records spanning time/frequency domain, RMS/peak, exp(+i omega t)/conjugate convention, source phase, port-current and normal orientation, primal/dual FEEC degrees, topology/harmonic/gauge data, passive/active/dispersive material state, PEC/PMC/impedance/port/PML/BEM boundary ownership, and source/storage/dissipation/radiation loss owners. Exact positive twins are paired with missing, conflicting, wrong-half-factor, conjugation, boundary-owner and loss-double-count twins. ORACLE: independent canonical crosswalk and refusal checker; no PDE solve. SEEDS: i14/i14-fullwave-schema-crosswalk/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fullwave-schema-core-holdout",
            r#"HOLDOUT_KIND: PublicReplayCore. GENERATOR: unseen FullWaveProblem convention/crosswalk records with sign, source, field-degree, relative-boundary, gauge, material-class, port, exterior, RMS/peak and loss-owner traps. Sole consumer i14-fullwave-schema-core. Public deterministic Core replay holdout; it proves schema admission mechanics only, never PDE stability or accuracy. SEEDS: i14/i14-fullwave-schema-core-holdout/case, core holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-fullwave-pml-dispersion",
            r#"DECK: open-domain Maxwell cases with analytic radiation or oversized-domain references, PML thickness/stretch sweeps, grazing/evanescent incidence, corners, anisotropic/dispersive/negative-index-adjacent passive media, late-time transients, pole proximity, and mesh/time refinement. Explicitly separate physical material dissipation from PML absorption and radiation boundary flux. Include unstable, reflecting, energy-creating, and deceptively short-window twins. SEEDS: i14/i14-fullwave-pml-dispersion/case, development and falsifier ranges."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fullwave-bem-scattering",
            r#"DECK: separately typed (1) closed PEC EFIE/MFIE/CFIE spheres and canonical conductors, (2) closed penetrable dielectric PMCHWT/Mueller/JMCFIE transmission scatterers with sided traces, and (3) admitted open PEC screen-EFIE cases with screen trace and edge-singularity spaces. Include interior-resonance traps, near-touching surfaces, multiply connected conductors, low-frequency and dense-mesh regimes, wrong-formulation/open-surface refusal twins, and cross-formulation cases only where mathematics permits. Small problems carry independent dense high-precision boundary-integral matrices/solves; larger cases remain formulation/consistency evidence only until independently certified. SEEDS: i14/i14-fullwave-bem-scattering/case, development and falsifier ranges."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fmm-dense-envelope-adversaries",
            r#"GENERATOR: formulation-bound Maxwell dense/FMM operator pairs spanning tree shapes, near/far boundaries, expansion/interpolation orders, vector/tangential coupling, complex outgoing kernels, low/high electrical size, nonuniform panels, near touching, adjoints, cancellation frontiers and measured dense crossover. Include false monotonic-tolerance assumptions, nondeterministic traversal, underestimated translation error, wrong kernel branch and formulation-error laundering twins. ORACLE: independently assemble dense matvecs/solves and separate FMM, iteration, quadrature, geometry and discretization envelopes. SEEDS: i14/i14-fmm-dense-envelope-adversaries/case, development and falsifier ranges."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fullwave-crosswalk-max",
            r#"DECK: maximal MTL/PEEC/full-wave overlap cases with common port and field-trace observables, topology-aware winding/cut representatives, source normalization, unique loss ownership, independent refinement ladders, and signed error decompositions. Include resonance, electrically long, radiation-dominant, retardation, rank/topology change, near-cancellation, and validity-gap cases. Agreement alone is never an error bound. SEEDS: i14/i14-fullwave-crosswalk-max/case, development and falsifier ranges."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fullwave-max-holdout",
            r#"HOLDOUT_KIND: PublicReplayMax. GENERATOR: unseen full-wave FEEC/PML cases emphasizing multiply connected topology, resonance, high contrast, dispersive loss, grazing radiation, relative boundaries, late-time growth and validity-region escape. Sole consumer i14-fullwave-feec-max. Public deterministic maximal replay holdout only; no population, laboratory, or universal theorem authority. SEEDS: i14/i14-fullwave-max-holdout/case, maximal holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-bem-formulation-max-holdout",
            r#"HOLDOUT_KIND: PublicReplayMax. GENERATOR: unseen closed-PEC, penetrable-dielectric and admitted open-screen scattering cases with resonance, low-frequency, near-touching, topology, sided-trace, edge-singularity and formulation/trace traps. Sole consumer i14-bem-formulation-max. Public deterministic Max replay only; no physical-validation or universal formulation authority. SEEDS: i14/i14-bem-formulation-max-holdout/case, maximal holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-fmm-acceleration-max-holdout",
            r#"HOLDOUT_KIND: PublicReplayMax. GENERATOR: unseen formulation-bound dense/FMM comparisons with hostile tree geometry, vector/tangential kernels, near/far boundary cases, adjoints, algorithm-route changes, nonmonotone observed error, deterministic traversal, cancellation and crossover traps. Sole consumer i14-fmm-acceleration-max. Public deterministic Max replay only; no physical-validation or universal performance authority. SEEDS: i14/i14-fmm-acceleration-max-holdout/case, maximal holdout range."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-safety-case-source-victim",
            r#"DECK: synthetic assurance graphs connecting source, coupling path, probe/model evidence, victim mode, hazard, control, assumption, monitor, owner, expiry, evidence color, and residual uncertainty. Include stale assumption, missing owner, hidden common cause, scope mismatch, superseded standard, simulation-only compliance claim, and untraceable edge twins. No licensed standard text and no legal/regulatory verdict is embedded. SEEDS: i14/i14-safety-case-source-victim/graph, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-standards-crosswalk-adversaries",
            r#"DECK: synthetic clause/digest and AP242/EMC/machine crosswalk graphs with stale editions, corrigendum mismatch, missing mandatory clauses, wrong applicability, occurrence/unit/frame loss, configuration mismatch, favorable Unknown omission, borrowed evidence and license/custody defects. No restricted text or actual conformance authority. ORACLE: independently reconstruct clause coverage, adapter loss and evidence-color graph. SEEDS: i14/i14-standards-crosswalk-adversaries/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-laboratory-validation-adversaries",
            r#"DECK: synthetic specimen/laboratory packages with expired or wrong-range calibration, cable/fixture/chamber drift, calibration-validation leakage, trace truncation, retry cherry-picking, covariance omission, censored/nonfinite samples, coordinate mismatch, post-hoc QoI/scale selection and favorable aggregation. No real laboratory or proprietary bytes. ORACLE: independent custody, calibration propagation and AcceptanceCard audit. SEEDS: i14/i14-laboratory-validation-adversaries/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-bearing-population-adversaries",
            r#"DECK: synthetic bearing-population records with lot/supplier/material/lubricant confounding, duty/environment shift, convenience sampling, left/right/interval censoring, competing risks, missingness, inspection misclassification, repeated-unit dependence, optional stopping, sparse-tail failure and estimand drift. No actual production-population authority. ORACLE: independent frame, censoring, event and coverage audit. SEEDS: i14/i14-bearing-population-adversaries/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-emc-reliability-adversaries",
            r#"DECK: synthetic source/harness/bond/shield/victim population records with configuration and lot shift, environmental/duty confounding, repeated-unit dependence, calibration drift, left/right/interval censoring, competing upset/recovery events, MNAR missingness, convenience sampling, optional stopping, sparse tails, invalid importance weights, post-hoc subgrouping and estimand drift. No physical-population authority. ORACLE: independent frame, event, missingness, stopping and coverage audit. SEEDS: i14/i14-emc-reliability-adversaries/case, development range."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-emc-reliability-max-holdout",
            r#"HOLDOUT_KIND: GovernedPopulationSlot. COMMITMENT SLOT: a governed source/harness/bond/shield/victim deployment population with frozen frame, outcome, censoring/missingness, dependence and estimand. No public seed or authored bytes confer population authority. One atomic RealizationCommitted FrozenManifest::amend authority transaction must replace this slot with its typed inaccessible External population root, install i14-external-emc-reliability-population-pack as its same-ID typed External discharge-envelope root, remove that Waiver row through its verified DischargeReceipt, verify the AmendmentRecord and advance the authority head before joint reveal. Sole consumer i14-emc-reliability-validation-max. Slot-only replacement, envelope-only replacement, waiver-only retirement, split or sequential access, or post-reveal amendment is IntegrityFailed."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-laboratory-validation-max-holdout",
            r#"HOLDOUT_KIND: GovernedPhysicalSlot. COMMITMENT SLOT: preregistered disjoint calibration, as-built model-input and untouched validation strata for EMC specimen configurations, calibrated raw traces and componentwise validation cases. No public seed or authored bytes confer physical-validation authority. Before candidate-side access, an independent custodian/fs-vvreg capability creates salted or equivalently hiding calibration/model-input content commitments, a hiding validation source-universe/frame commitment, disjoint-membership commitment, exact validation-selection algorithm and pre-candidate secret-seed/VRF or equivalent non-adaptive selection commitment. GovernanceCommitted freezes those opaque identities, stratum rules, calibration/model-instantiation procedure, permitted updates, contamination predicates, AcceptanceCards and every validation rule. Candidate builders, fitters, checker/threshold owners and their transitive capabilities receive only opaque commitment identities and no validation bytes, membership witnesses, labels, aggregates, derived statistics, selection outputs or commitment-opening material. AuthorizedCalibration exposes only committed calibration and as-built model-input strata under logged least privilege; AsBuiltModelInstantiation may fit/select only under the frozen procedure. CandidateFrozen binds the resulting model, parameters and rules, and an independent contamination receipt names and audits every candidate-side principal and transitive capability. One atomic RealizationCommitted FrozenManifest::amend authority transaction must replace this slot with its joined typed External root carrying separately addressable calibration, model-input and validation roots, membership proofs to the pre-access commitments, mutual disjoint-membership proof and non-adaptive selection proof; install i14-external-emc-laboratory-calibration-pack and i14-external-asbuilt-specimen-geometry-pack as distinct same-ID typed External discharge-envelope roots; remove both Waiver rows through distinct verified DischargeReceipts; verify the AmendmentRecord; and advance the authority head before untouched validation reveal. Sole consumer i14-laboratory-validation-max. Candidate-side validation access or opening before freeze, calibration-to-validation leakage, hidden model update, adaptive validation selection, partial discharge, slot-only or envelope-only replacement, stratum aliasing, split retirement, sample reordering without timestamp identity, sequential reveal or retry is IntegrityFailed."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-bearing-population-max-holdout",
            r#"HOLDOUT_KIND: GovernedPopulationSlot. COMMITMENT SLOT: a governed production bearing population with frozen frame, per-unit lineage/exposure, event and inspection definitions, censoring/missingness/competing-risk rules and estimand. No public seed or authored bytes confer population authority. One atomic RealizationCommitted FrozenManifest::amend authority transaction must replace this slot with its joined typed inaccessible External population root, install i14-external-bearing-population-reliability-pack and i14-external-bearing-population-metrology-pack as distinct same-ID typed External discharge-envelope roots, remove both Waiver rows through distinct verified DischargeReceipts, verify the AmendmentRecord and advance the authority head before joint reveal. Sole consumer i14-bearing-population-validation-max. Partial discharge, slot-only or envelope-only replacement, split retirement, sequential access or post-reveal amendment is IntegrityFailed."#,
            Partition::HeldOut,
        ),
        authored_fixture(
            "i14-passive-composition-theorem-card",
            r#"TARGET CARD ONLY: formalize a finite typed paired compatibility-sheaf/balance-cosheaf assembly over a cellular site with local causal dissipative relations, lower-bounded storage, supply/dissipation rates, frozen signal/history spaces, oriented relative boundary ports, a chain-level natural power pairing into the orientation/dualizing cosheaf, relative total complexes/hypercohomology, restriction/corestriction/crosswalk maps, pairwise and higher-overlap cocycle/coherence laws, relative-boundary compatibility, checked chain-homotopy representative invariance, and a maximally isotropic lossless Dirac relation with respect to a nondegenerate ambient split-signature power pairing. Require clean composition, no latent unpaired port, well-posed descriptor/feedback closure, initialization/history, a KYP/passive-realization bridge for frequency-domain positive-real components, checked cover-refinement comparison, and explicit gluing/solver/quadrature/time defects. For every admitted trajectory and finite horizon T, target S(x(T))-S(x(0)) <= integral_[0,T] <e_ext,f_ext> dt - integral_[0,T] d(t) dt + Delta_accounted(T), d>=0, with exact units/measure and dependency-aware outward bounds. Target conclusions distinguish (1) an exact global loss identity requiring zero unaccounted defect and exact owner equality, (2) defect-tolerant global passivity when integral d>=Delta_accounted for every admitted trajectory/horizon, and (3) quantified robust strict passivity when integral d-Delta_accounted>=mu*norm_signal^2 for one frozen mu>0; infinite-horizon authority adds coercivity/detectability and limit/integrability premises. Nonzero relative-hypercohomology or fidelity-holonomy obstruction forces escalation. This version freezes intent, not a theorem AST, formal proposition, proof, translation, axiom closure, or nonvacuity result."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-hypercohomology-obstruction-card",
            r#"TARGET CARD ONLY: formalize the relative total compatibility/balance complex of a finite typed cellular harness site, including coefficient/torsion policy, boundary, restriction/corestriction and crosswalk maps, total-differential signs, observable sheaf and independent executable canonicalization. Prove the exact necessary obstruction, the additional hypotheses under which vanishing plus a constructive witness is sufficient for global gluing, and chain-homotopy-invariant localization of every tied minimal obstruction support. Include nonzero, torsion, relative-boundary, false-vanishing and nonunique-repair twins. This card is not a proposition AST or proof."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-cover-refinement-naturality-card",
            r#"TARGET CARD ONLY: formalize two finite covers/sites, their common refinement/cofinal comparison, relative boundaries, coefficient/torsion policy, local complexes, comparison transformations, chain homotopies, trace/power pairing and defect-owner transport. State exact Leray/acyclicity/cofinality/quasi-isomorphism alternatives sufficient to preserve sections, obstruction classes, boundary supply and passive verdict, plus explicit comparison obstructions for invalid coarsening. Include duplicated-loss and false-universal-invariance twins. This card is not a proposition AST or proof."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-kyp-sheaf-bridge-card",
            r#"TARGET CARD ONLY: formalize a finite-dimensional real or realified-complex rational/descriptor generalized-positive-real realization with regular pencil, impulse-free consistent initialization, quotient/lossless modes, simple imaginary-axis/zero/infinity storage poles and PSD residues, exact port orientation and supply matrix. Prove the strength-matched KYP/storage witness consumed by the local sheaf relation without upgrading generalized PR to strict PR or sampled-band dissipativity to a time-domain theorem. Include irregular, indefinite-residue, wrong-supply-sign and nonminimal twins. This card is not a proposition AST or proof."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-fidelity-descent-theorem-card",
            r#"TARGET CARD ONLY: formalize conditions under which a certified fine-rung QoI enclosure, stable projection/lift and trace maps, topology/reference/phasor-compatible crosswalks, unique loss ownership, and signed local residual/geometry/material/solver/time/model-form budgets descend to a coarser rung while preserving a declared observable and total error bound. Include countertargets for resonance, topology/rank change, invalid asymptotics, cancellation, discontinuous events, and dependence-laundered budgets. This card is not a theorem AST or proof."#,
            Partition::Development,
        ),
        authored_fixture(
            THEOREM_POLICY_FIXTURE,
            r#"POLICY: theorem authority begins only in a pre-candidate successor that freezes (1) a canonical proposition-and-definition AST with units, quantifiers, domains, exclusions, validity predicates, signal/history spaces, boundary/reference/gauge/topology semantics, defect arithmetic, and no-claim boundary; (2) deterministic total AST-to-formal translation with binding receipts to the exact generated declaration; (3) runtime-premise schemas and sound admission checkers bound to the same AST; (4) complete transitive axiom closure; (5) strength-matched nonvacuity witnesses, negative twins, and independent countermodel/admission checkers. The only permitted ambient Lean axioms are exactly {propext, Quot.sound, Classical.choice}; `sorryAx`, unbound declarations, hidden native_decide-style executable trust, proof by production boolean, or a proposition weaker than the manifest claim refuse activation. A kernel proof is necessary but not sufficient for physical validity. AUTHORITY_CONTRADICTION: an independently admitted GenuineCountermodel against the exact same proposition/definition/axiom digest as an admitted kernel proof cannot be resolved by voting or immediate refutation. It quarantines proof and runtime authority, preserves both artifacts, and executes a deterministic fault tree over AST/proposition binding, assumptions, runtime-premise admission, arithmetic/enclosures, oracle independence, transitive axiom report and proof-kernel/TCB integrity; only a versioned defect resolution may classify the revision Proved or Refuted. Bare cohomology equivalence, exact d*d=0 incidence, component passivity, or naming a Dirac structure has no authority for global causal passivity."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-theorem-falsifier-grammar",
            r#"VERSION-1 TARGET GRAMMAR: candidate families for rational finite MTL/PEEC/circuit components, paired compatibility/balance diagrams, relative-complex winding/cut cases, crosswalk defects, well-posed feedback, dispersive/delayed/nonnormal/time-discrete/PML/BEM/FMM cases, uncertainty dependence, and negative twins. A successor must replace this prose with canonical executable grammar ASTs, exact validity and theorem-premise predicates, coefficient/degree/cardinality bounds, rank/unrank/sharding, canonicalization, exclusion order, cardinality/bijection proof, resource preflight, Merkle completeness root, independent classification, and minimization. This version has no exhaustive-search or theorem-survival authority."#,
            Partition::Development,
        ),
        authored_fixture(
            "i14-maximal-adversaries",
            r#"GENERATOR: cross-domain adversaries combining topology/reference changes, MTL/PEEC/full-wave validity escape, resonance, PML late-time instability, BEM/FMM envelope error, source/probe aliasing, victim event discontinuity, correlated uncertainty, mitigation violating the frozen safety envelope, assurance laundering, theorem-premise boundary cases, and reason-coded negative controls. Every candidate binds its grammar revision, rank, validity decision, model fingerprints, and independent adjudication state. SEEDS: i14/i14-maximal-adversaries/falsifier, falsifier range."#,
            Partition::Development,
        ),
    ]
}

#[allow(clippy::too_many_lines)]
fn i14_obligations() -> Vec<ObligationRow> {
    const UNIT_CASES: &[&str] = &[
        "happy",
        "empty",
        "boundary",
        "max",
        "error",
        "unit-dimension",
        "tie-break",
        "cancellation",
        "migration",
    ];
    const G5_MATRIX: &str = "threads {1,2,7} x shards {1,4} x mode {deterministic} x ISA {Apple-aarch64,x86_64}; bitwise comparison is authorized only for the authoritative V2 canonical terminal-result digest under an identical recorded logical event/cause trace, identical bound verification-receipt identities and identical topology/ISA/toolchain fingerprint. That canonical terminal-result digest binds the canonicalized logical event/cause trace from genesis through the first terminal boundary, including the cancellation-card identity, every prefix decision, selected disposition/cause, request ids, logical sequences, boundary ordinals, lifecycle projection and semantic payload, while excluding calibrated monotonic timestamps, deadlines, live watchdog arrival and clock calibration; the separate receipt-bound telemetry-envelope digest adds those noncanonical timing/calibration fields and is intentionally not bitwise compared. Cross-fingerprint comparison uses canonical discrete verdicts and preregistered numeric bands; accessibility/agent parity requires every campaign and replay through the documented noninteractive CLI with no privileged UI or hidden operator step; performance evidence binds p50/p95/p99 time, peak memory, quiet-host state, machine fingerprint and scale to the manifest budget, and smoke evidence never promotes a core/max performance claim";
    macro_rules! g0_contract {
        ($specific:literal) => {
            concat!(
                $specific,
                "; manifest binding: this obligation row is the exact per-cluster manifest slice; its row decks are the G1/G2 deck identities, each covered claim supplies the only authorized QoI/unit/tolerance/oracle arithmetic, and cross-crate/IR/API round trips must preserve the frozen schema, units, references, evidence color and digest; drift, an unbound metric, or an unnamed skip refuses"
            )
        };
    }
    macro_rules! g4_contract {
        ($specific:literal) => {
            concat!(
                $specific,
                "; numeric cancellation contract: every CancellationCard names its semantic work unit, total resource ceiling, count-defined logical poll-tile quantum, admitted worst-case indivisible-item time, asynchronous watchdog quantum and external heartbeat bound. The logical tile partition and boundary ordinals are immutable functions of admitted input and logical work identity. The scheduler polls before item 0 and immediately before and after every logical tile; the watchdog polls independently without ending or repartitioning a tile and supervises external heartbeats.",
                " Core total ceilings are 4096 graph/trace records, 16384 field/quadrature unknowns, 1024 search/tree nodes, or 256 formal declarations; a Core logical poll tile contains at most 64 graph/trace records, 256 field/quadrature unknowns, 16 search/tree nodes, or 4 formal declarations. The Core asynchronous watchdog polls at intervals <=25 ms without changing logical tile membership or order. Core request-to-observation <=250 ms, at most 32 in-flight children, drain-trigger-to-drained <=2000 ms and drained-to-finalized <=2000 ms.",
                " Max total ceilings are 16384 graph/trace records, 65536 field/quadrature unknowns, 4096 search/tree nodes, or 1024 formal declarations; a Max logical poll tile contains at most 256 graph/trace records, 1024 field/quadrature unknowns, 64 search/tree nodes, or 16 formal declarations. The Max asynchronous watchdog polls at intervals <=100 ms without changing logical tile membership or order. Max request-to-observation <=1000 ms, at most 128 in-flight children, drain-trigger-to-drained <=8000 ms and drained-to-finalized <=8000 ms. Admission refuses any indivisible item or external heartbeat whose demonstrated response bound exceeds its tier watchdog quantum, and any logical tile whose admitted worst-case request-to-boundary bound exceeds its tier request-to-observation SLO.",
                " execution.requested means job admission, not cancellation. cancellation.requested binds request id, scope root, coordinator-assigned globally unique deterministic logical sequence and calibrated monotonic timestamp; cancellation.observed binds the same request id, its own globally unique sequence, observing tile id, optional latest-completed-boundary ordinal and calibrated monotonic timestamp; execution.cancelled is the exactly-once terminal-selection event binding that request, terminal boundary and selected Cancelled cause. request-to-observation is their timestamp delta; drain-trigger-to-drained starts at the referenced on-time observation, first nanosecond after the missed inclusive request deadline, structurally admitted receipt-bound infrastructure-failure onset whose receipt HELM/ledger authenticates, or drain-start timestamp for CancellationObserved, ObservationTimeoutDrain, InfrastructureFailure, or NonCancellationDrain respectively and ends at the last execution.drained descendant; drained-to-finalized begins there and ends at execution.finalized. Only a request whose logical sequence precedes the candidate boundary and whose scope root occurs in that candidate's acyclic root-to-leaf scope ancestry participates; request slices and observer-tile catalogs are canonicalized before validation, so valid decisions and malformed refusals are input-order invariant, and multiple relevant requests are totally ordered by their unique logical sequence. I14_MAX_CANCELLATION_REQUESTS_V1=16384, I14_MAX_SCOPE_ANCESTRY_V1=256 and I14_MAX_OBSERVER_TILES_V1=128 bound the legacy local schema. i14_select_terminal_boundary_v1 validates those predicates, unique scope identities, strict request-before-observation causality with strict-pre-boundary participation in the frozen cut, admitted observing-tile identity, and that any participating observation's latest-completed boundary strictly precedes the candidate; a retained observation later than the candidate must name a latest-completed boundary ordinal at least the candidate ordinal. Calibrated monotonic timestamps must be nondecreasing in globally unique logical-sequence order across requests, observations and the candidate boundary. It also validates deadlines, caps and duplicate ids/sequences, then implements the same-boundary InfrastructureFailed > TimedOut > Cancelled > BudgetExhausted > Completed selector, but V1 proves only local arbitration at a caller-supplied boundary and has zero promotion authority. I14CancellationCardV2 admits Core, ordinary Max, or the explicit MaxTheoremFalsifier subtype, binds the semantic work-unit digest, campaign wall budget, hard logical-memory ceiling, exact count resource kind/ceiling/tile quantum, resource authority, deterministic partition, execution-environment fingerprint and external-child catalog, preserves the frozen 90-minute/18-hour/24-hour envelopes, and refuses logical-tile, indivisible-item or external-heartbeat response bounds wider than the frozen tier SLO. It refuses zero bounds, over-tier memory/count/tile contracts and inconsistent external-child policy. The V2 schema path enforces the exact 250 ms/1000 ms request-deadline delta, the tier-specific 32/128 observer and in-flight-child caps, campaign expiry, resource/work-frontier monotonicity and watchdog coverage. An actual observation closes the spawn frontier at its logical event; without one, request/deadline alone does not synthesize a logical cut and timeout/failure drain start closes it. All descendants and losing races drain. A multi-kind campaign freezes one card per independently governed work leaf; no card silently changes units.",
                " The coordinator selects the first terminal-eligible logical tile boundary in deterministic logical tile/event order. i14_select_first_terminal_boundary_v2 requires a nonempty trace beginning at genesis ordinal 0, at most I14_MAX_TERMINAL_BOUNDARIES_V2=4096 contiguous ordinals, strictly increasing logical sequences, nondecreasing calibrated boundary times, an immutable scope path and observer catalog, and no supplied record after the first Selected decision; its opaque nonterminal frontier certificate cannot be converted into a terminal result. The reference selector additionally refuses more than I14_MAX_TERMINAL_ARBITRATION_PAIRS_V2=1048576 boundary/request pairs before repeated legacy arbitration. Temporal precedence across boundaries is distinct from same-boundary cause precedence: an earlier eligible boundary always wins over a later nominally higher-priority cause. A cancellation request whose coordinator-assigned logical sequence precedes a normal Completed or BudgetExhausted candidate makes that candidate ineligible until the same request is observed; observation within its inclusive SLO selects Cancelled when no higher-priority cause occurs, while the first nanosecond after a missed observation deadline selects TimedOut. Within one boundary the cause priority is InfrastructureFailed > TimedOut > Cancelled > BudgetExhausted > Completed. A supervisor/authentication/drain protocol failure yields InfrastructureFailed; any overall campaign wall-clock expiry or request-observation/drain/finalization deadline miss yields TimedOut. BudgetExhausted requires nonempty frozen work plus a receipt-bound next-work quantum that the governor rejected because it would cross the exact ceiling; completion exactly at the ceiling remains Completed, and only normal completion of all scheduled work before every deadline/cancellation with an empty work frontier yields Completed. I14TerminalLifecycleTraceV2 makes execution-started, drain-started, drained and finalized events mandatory and globally sequence-unique; binds active/drained descendant counts; and carries I14DrainTriggerV2={CancellationObserved,ObservationTimeoutDrain,InfrastructureFailure,NonCancellationDrain}. At the drain-start cut the trigger is the earliest effective in-scope on-time observation, first nanosecond after a missed inclusive observation deadline or structurally admitted receipt-bound infrastructure-failure onset under the frozen calibrated order, with effective time first, causal tie rank second (Infrastructure=0, Observation=1, Timeout=2; distinct from wire tags), causal logical sequence third and stable identity fourth; absent an earlier candidate, drain start is the non-cancellation trigger. HELM/ledger authenticates the onset receipt. Canonical bytes bind the derived variant plus request id or infrastructure-onset logical sequence and bind the closed infrastructure source plus independent verification-receipt identity, while raw observation/deadline/onset times remain telemetry. Every terminal path requires I14SpawnFrontierEvidenceV2: request_id names the earliest actual pre-drain observation cut, may differ from the drain trigger, and is None when drain start itself closes the frontier. The audit binds the unconditional child semantic root; a separate child raw root is telemetry-only; and an independent child semantic-verification receipt identity is canonical. Watchdog, tile-poll and heartbeat evidence likewise split canonical semantic roots from telemetry-only raw roots and independent semantic-verification receipt identities. Local code validates only structural/internal consistency, caps, trigger arithmetic and failure reflection; the HELM/ledger promotion gate must authenticate issuers and independently verify membership, completeness, child/external acknowledgements, raw-root derivations and publication atomicity. When a complete multi-kind watchdog stream is supplied, this layer additionally derives and checks its versioned raw root. A typed infrastructure-onset witness binds source={WatchdogCoverage,TilePollCoverage,ExternalHeartbeatCoverage,DescendantDrain,SpawnAfterFrontierClosure,Supervisor,Authentication,DrainProtocol,PublicationProtocol}, event and independent verification-receipt identity; local derived source tags must match the corresponding failed evidence, while HELM/ledger authenticates every receipt and all generic protocol sources. Timeout onset is derived, never caller-selected, as the earliest first nanosecond outside the inclusive trigger-to-drained or drained-to-finalized cap; deadline expiry precedes a boundary at the same nanosecond, and the first boundary at or after that derived time must latch it. Real watchdog, descendant-drain or spawn-after-frontier failure remains hashable only when selected InfrastructureFailed reports it; tile-poll or external-heartbeat/termination/publication failure follows the same rule, and real trigger-to-drained or drained-to-finalized deadline failure remains hashable only when selected TimedOut or higher-priority InfrastructureFailed reports it. Watchdog and calibrated wall-time observations are receipt-bound telemetry: their raw timing fields never change logical tile membership or boundary ordinal and are excluded from canonical-result bytes, but when their recorded semantics select a different trigger or terminal cause, the canonical disposition and digest honestly change. Late triggers are retained but never rewrite an already finalized boundary through a telemetry-only verification-receipt-bound late-event tail whose authentication is deferred to HELM/ledger, while a lifecycle failure before final receipt upgrades the cause to TimedOut or InfrastructureFailed. Each non-Completed disposition has zero partial promotion authority. The schema-valid projection itself has no promotion authority until the consuming HELM/ledger verifier authenticates the bound issuer, capability, trust-policy, revocation, card and semantic-work receipt.",
                " The timeout logical field is not a synthetic onset event: it must equal the exact first latch-boundary logical sequence at or after the locally derived calibrated onset, eliminating caller-selected canonical identity.",
                " The canonical lifecycle projection also binds the child semantic trace and independent child semantic-verification receipt identity, mandatory spawn-frontier audit, watchdog/tile-poll/heartbeat semantic traces and independent semantic-verification receipt identities, the receipt-bound infrastructure-onset latch, the locally derived timeout latch and every derived lifecycle-failure bit. HELM/ledger authenticates the corresponding receipts; local code only validates structural binding and consistency. Its paired telemetry binds the distinct child/watchdog/tile-poll/heartbeat raw roots and calibrated failure-onset times. I14_DRAIN_TRIGGER_ENCODING_V2_KAT_HEX pins all four exact trigger tags and U64LE payloads; I14_INFRASTRUCTURE_FAILURE_ONSET_ENCODING_V2_KAT_HEX pins the present canonical presence/sequence/source/receipt form; and I14_WATCHDOG_RAW_TRACE_V2_KAT_HEX pins the complete multi-kind raw-watchdog byte layout. These evidence-layout KATs are distinct from the four pre-existing V2 schema-layout KATs named above; a checked-in KAT alone is not evidence that an independent encoder reproduced it.",
                " Terminal receipts preserve eight orthogonal axes exactly: ExecutionDisposition={Completed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed}; ClaimAdjudication={Pending,Supported,Failed,Refuted,Unknown}; EvidenceCompleteness={CompleteEvidence,PartialEvidence,NoEvidence}; EvidenceIntegrity={IntegrityVerified,IntegrityFailed}; InputValidity={WellFormedInput,MalformedInput}; DomainApplicability={Admitted,OutOfDomain,Indeterminate}; OperationalSupport={SupportedOperation,UnsupportedOperation}; ReceiptValidity={WellFormedReceipt,MalformedReceipt}; one axis never substitutes for another. EvidenceIntegrity covers evidence bytes, custody and checker integrity; InputValidity covers structural input validity; DomainApplicability is evaluated only for well-formed input; ReceiptValidity covers schema, causal-event and cross-axis consistency.",
                " Combination validity is fail-closed: MalformedInput forces DomainApplicability=Indeterminate before any other combination check. ClaimAdjudication in {Supported,Failed,Refuted} requires ExecutionDisposition=Completed, EvidenceCompleteness=CompleteEvidence, EvidenceIntegrity=IntegrityVerified, InputValidity=WellFormedInput, DomainApplicability=Admitted, OperationalSupport=SupportedOperation and ReceiptValidity=WellFormedReceipt; ClaimAdjudication=Pending requires a non-Completed disposition or non-complete evidence; MalformedInput, DomainApplicability in {OutOfDomain,Indeterminate}, OperationalSupport=UnsupportedOperation, non-Completed execution, incomplete evidence, failed integrity or malformed receipt permits only ClaimAdjudication in {Unknown,Pending}; every forbidden combination is ReceiptValidity=MalformedReceipt, never evidence corruption by definition.",
                " I14TerminalStatusV1 and i14_evaluate_terminal_status_v1 implement TerminalStatusTruthTableV1, which exhaustively enumerates all 3600 Cartesian tuples, retains the raw producer tuple, records every normalization action, normalizes each combination and maps it to one primary exit; a malformed input paired with a raw non-Indeterminate domain is a cross-axis contradiction and marks the receipt malformed. G0 enumerates the typed evaluator and compares i14_terminal_status_table_digest_v1 against I14_TERMINAL_STATUS_TABLE_V1_DIGEST_HEX. The deterministic primary CLI exit code is selected without discarding any receipt axis, in this exact fail-closed precedence: ReceiptValidity=MalformedReceipt or InputValidity=MalformedInput or EvidenceIntegrity=IntegrityFailed -> 60; ExecutionDisposition=InfrastructureFailed -> 70; ExecutionDisposition=Cancelled -> 20; ExecutionDisposition=TimedOut -> 21; ExecutionDisposition=BudgetExhausted -> 22; ClaimAdjudication=Failed -> 10; ClaimAdjudication=Refuted -> 40; DomainApplicability in {OutOfDomain,Indeterminate} or OperationalSupport=UnsupportedOperation -> 30; EvidenceCompleteness in {PartialEvidence,NoEvidence} -> 50; ClaimAdjudication=Unknown -> 30; 0 only for ExecutionDisposition=Completed with every requested ClaimAdjudication=Supported, EvidenceCompleteness=CompleteEvidence, EvidenceIntegrity=IntegrityVerified, InputValidity=WellFormedInput, DomainApplicability=Admitted, OperationalSupport=SupportedOperation and ReceiptValidity=WellFormedReceipt. Pending has no standalone primary exit: a well-formed Pending receipt is selected earlier by its execution, domain/support or incomplete-evidence cause. Every declared exit code has a well-formed-receipt witness, and every precedence disjunct has an explicit witness including the necessarily malformed-receipt branch. IntegrityFailed, MalformedInput, MalformedReceipt and InfrastructureFailed are evidentiary/admission/operational states, never scientific refutation.",
                " Logical event identity, causal order and semantic payload are deterministic. I14CanonicalTerminalResultInputV2 plus i14_canonical_terminal_result_v2 and i14_canonical_terminal_result_digest_v2 revalidate the complete genesis-to-first-terminal prefix, validated cancellation card and terminal lifecycle before accepting a normalized receipt whose execution matches the causal disposition; they recompute the clock-free logical execution root from canonical boundary/resource/work/request semantics and bind its independent verification-receipt identity, request-inclusive terminal-prefix digest, cancellation-card identity, every boundary prefix decision, raw and normalized terminal axes, normalization actions, exit projection, semantic payload digest, cause candidates, immutable strict-pre-boundary request cut, lifecycle logical sequences, unconditional child trace, tile/watchdog/external semantic roots and verification-receipt identities, and every derived lifecycle-failure bit. A nonterminal frontier has no result constructor, and post-terminal boundaries are refused. The authoritative V2 canonical terminal-result digest binds the canonicalized logical event/cause trace, including selected disposition/cause, request ids, logical sequences, boundary ordinals and semantic payload, while excluding calibrated monotonic timestamps, deadlines, live watchdog arrival and clock calibration. Here authoritative means the frozen V2 byte/schema contract; only an independently authenticated HELM/ledger adjudication receipt can authorize promotion. i14_telemetry_envelope_digest_v2 revalidates the authoritative result through the same local implementation and binds its digest plus campaign start, every raw boundary/watchdog timestamp, every raw request/observation timestamp and deadline, lifecycle times, complete-versus-subset watchdog-sample status, raw watchdog/external counts, at most I14_MAX_WATCHDOG_OBSERVATIONS_V1=4096 canonically ordered watchdog samples, a telemetry-only verification-receipt-bound late-event tail whose authentication is deferred to HELM/ledger, and the clock-calibration artifact. I14_CANCELLATION_CARD_V2_KAT_HEX, I14_TERMINAL_PREFIX_V2_KAT_HEX, I14_CANONICAL_TERMINAL_RESULT_V2_KAT_HEX and I14_TELEMETRY_ENVELOPE_V2_KAT_HEX pin the four V2 byte layouts as checked-in known answers; they do not by themselves prove an independent encoder exists or agrees. The V1 canonical and telemetry layouts and I14_CANONICAL_TERMINAL_RESULT_V1_KAT_HEX/I14_TELEMETRY_ENVELOPE_V1_KAT_HEX remain only for legacy ledger readability and have no first-terminal or lifecycle promotion authority. G5 bitwise comparison applies only to the V2 canonical result digest under an identical recorded logical event/cause trace, identical bound verification-receipt identities and identical topology/ISA/toolchain fingerprint, never to a live telemetry envelope. A different live watchdog/expiry trace may honestly select a different disposition and has no cross-run bit-stability claim. Every encoded event is bounded to 64 KiB after schema-aware redaction; an oversized field becomes a governed digest/slot plus typed explicit-truncation metadata. Licensed text, secrets, specimen identities, governed holdout/validation/population bytes before or after controlled reveal, and derived sensitive slices never enter a public artifact. Sanitized manifest, adjudication receipt, logs, oracle output and replay capsule map exactly to I14_EVIDENCE_DURABLE; sanitized minimized counterexamples, refutations and FailureBundles map exactly to I14_FAILURE_PERMANENT; raw licensed/secret/specimen/governed-holdout bytes, derived sensitive slices and unredacted diagnostic slices map exactly to encrypted capability-controlled I14_GOVERNED_RESTRICTED with complete access ledger and class-specific retention/erasure policy. I14ArtifactCategoryV1 and i14_retention_rule_v1 exhaustively map every named category to its retention class, sanitization requirement, encryption/capability rule, complete-access-ledger requirement and retained class-specific retention/erasure decision. Schema-aware sanitization applies before retention to events, manifests, adjudication receipts, logs, oracle output, replay capsules, minimized counterexamples, refutations and FailureBundles, and the retention tail must preserve the complete access ledger and class-specific retention/erasure decision; independent no-orphan/waiver/drift lint and manifest-adjudication receipt are mandatory. A Core row is consumed only by I14.G4 after I14.G2/G3 and a Max row only by I14.G7 after I14.G5/G6; any stale/missing artifact or lifecycle state fails closed"
            )
        };
    }

    vec![
        ObligationRow {
            leaf: "i14-harnessgraph-core",
            claims_covered: &["i14-harnessgraph-identity-connectivity"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: typed native HarnessGraph assemblies, invalid occurrence aliases, open/unmated/fault states, repeated part-number and permutation-isomorphic twins; validity predicates: typed incidence, occurrence identity, connector/contact capacity, orientation/reference and unit closure; laws: canonical identity is permutation-stable, an occurrence edit changes only its lineage, connectivity and port signs agree with an independent incidence traversal, charge/current has no undeclared path; shrinkers: delete graph leaves, route segments and attributes while preserving the first identity/connectivity fault; replay seeds and ranges follow the Five Explicits and fixture specs"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                "i14-harnessgraph-canonical",
                "i14-harnessgraph-core-holdout",
            ],
            g3_relations: &[
                "renaming and insertion-order permutation preserve canonical occurrence/connectivity semantics and digest",
                "rigid route transforms preserve graph identity while changing only geometry-dependent descendants",
                "splitting then losslessly rejoining a route segment preserves terminal connectivity and conservation",
                "an explicit open/unmated transform removes exactly the declared path and never aliases a sibling occurrence",
            ],
            g4_schedule: g4_contract!(
                "request cancellation at graph ingestion, canonicalization, connectivity traversal, evidence serialization and heldout-boundary tiles; drain every child and oracle, finalize exactly one terminal receipt, retain the smallest failure bundle, then checkpoint, migrate, resume and fork at each stable frontier; resumed and uninterrupted canonical graph/evidence digests must agree and no partial authority may publish"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_harnessgraph_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-harnessgraph-core dsr quality --tool frankensim",
            obs_events: &[
                "harness.ingested",
                "harness.identity_checked",
                "harness.connectivity_checked",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_harnessgraph_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-ap242-adapter-core",
            claims_covered: &["i14-synthetic-ap242-adapter-mechanics"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: legally authored AP242-shaped subset records, repeated/nested occurrences, mixed unit/frame transforms, supported and unsupported connectivity semantics, loss-accounting and ambiguity twins; validity predicates: exact source root, authored subset membership, one-way occurrence/transform/connectivity map and exhaustive source disposition; laws: every record maps to one native identity, explicit loss, quarantine or refusal, occurrence reuse never aliases instances, and no unsupported construct disappears; shrinkers: remove records and assembly levels while preserving the first adapter/loss fault; replay seeds follow the adapter fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                "i14-ap242-adapter-mechanics",
                "i14-ap242-adapter-core-holdout",
            ],
            g3_relations: &[
                "source record and insertion-order permutation preserves the canonical native occurrence/loss result",
                "equivalent unit/frame presentation preserves physical coordinates only through the complete frozen transform",
                "reusing one product definition for two occurrences preserves two semantic occurrence identities",
                "adding an unsupported source record adds exactly one quarantine/loss/refusal disposition and cannot change a supported sibling",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during source decode, occurrence/transform reconstruction, native mapping, loss audit, oracle and heldout tiles; drain every parser/checker, finalize one adapter/loss receipt with no partial import authority, retain quarantined source records, then checkpoint, migrate, resume and fork at stable record/assembly frontiers; resumed and uninterrupted mapping and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_ap242_adapter_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-ap242-adapter-core dsr quality --tool frankensim",
            obs_events: &[
                "adapter.source_bound",
                "adapter.occurrence_mapped",
                "adapter.loss_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_ap242_adapter_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-rlgc-operator-core",
            claims_covered: &["i14-mtl-rlgc-operator-admission"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: analytic and fitted ordered multiconductor RLGC operators, quotient/reference transforms, quasi-static, generalized-PR, strict-PR, passive-realization, sampled-only and hostile pole/residue/interpolation twins; validity predicates: EmConventionCard, dimensions, conductor/reference order, quotient rank, source support, quasi-static storage premises and exact generalized-PR or weaker band authority class; laws: analytic cases and independent matrix/quotient checks agree, complete basis/reference transforms preserve physical ports and supply, boundary poles/residues obey their typed class, and sampled data never become analytic authority; shrinkers: reduce conductors, samples, poles and states while preserving the first operator/passivity defect; replay seeds follow the operator fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-mtl-coax-two-wire",
                "i14-rlgc-fit-adversaries",
                "i14-rlgc-operator-core-holdout",
            ],
            g3_relations: &[
                "consistent conductor permutation conjugates RLGC and ports while preserving physical terminal supply",
                "reference change with the complete quotient/source map preserves gauge-invariant port observables",
                "metre-to-millimetre rescaling with exact per-length unit transport preserves the physical operator",
                "adding a sampled frequency cannot upgrade band dissipativity to analytic generalized-PR authority",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during extraction ingestion, rational fitting, generalized-PR/KYP checking, interpolation, oracle and operator-heldout tiles; drain every speculative fit/checker, finalize one operator/source/passivity receipt, persist pole/residue/quotient state, then checkpoint, migrate, resume and fork at stable frequency/state frontiers; resumed and uninterrupted operator and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_rlgc_operator_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-rlgc-operator-core dsr quality --tool frankensim",
            obs_events: &[
                "rlgc.operator_admitted",
                "rlgc.passivity_checked",
                "rlgc.fit_refused",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_rlgc_operator_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-mtl-propagation-core",
            claims_covered: &["i14-mtl-passive-causal-propagation"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: admitted line operators, multiconductor routes, terminations, nearly repeated modes, convolution/state initializations, lossless/lossy/dispersive and reflection/delay/power twins; validity predicates: green RLGC operator, EmConventionCard, declared telegrapher signs, stable causal realization or typed causal operator, terminal orientation, support and unique storage/loss owners; laws: dV/dz=-(R+i omega L)I and dI/dz=-(G+i omega C)V, matched/open/short analytic limits and terminal power close, repeated-mode subspaces are deterministic, and invalid extrapolation refuses; shrinkers: reduce conductors, segments, modes, states and frequencies while retaining the first propagation fault; replay seeds follow the propagation fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-mtl-coax-two-wire",
                "i14-mtl-bundle-core-holdout",
            ],
            g3_relations: &[
                "segment subdivision and reassociation preserve propagation inside the composed error budget",
                "complete conductor/reference basis transport preserves gauge-invariant terminal voltage, current and power",
                "reciprocal lossless reversal with exchanged oriented ports satisfies the declared reciprocity twin",
                "mode-basis rotation inside a repeated invariant subspace preserves reconstructed physical fields and ports",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during modal/Schur propagation, convolution/state stepping, frequency sweep, oracle and propagation-heldout tiles; drain every solver/oracle, finalize one wave/terminal-power receipt with no partial propagation authority, persist segment/mode/history state, then checkpoint, migrate, resume and fork at stable frequency/segment frontiers; resumed and uninterrupted terminal and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_mtl_propagation_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-mtl-propagation-core dsr quality --tool frankensim",
            obs_events: &[
                "mtl.propagation_checked",
                "mtl.mode_subspace_bound",
                "mtl.power_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_mtl_propagation_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-peec-network",
            claims_covered: &["i14-peec-extraction-power-mor"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: analytic PEEC cells, connected/disconnected conductors, retardation states, gauge/reference transforms, dense and reduced networks, near-singular and invalid P-versus-C/loss/passivity twins; validity predicates: geometry/material/state support, P potential-coefficient semantics distinct from derived C, Lorenz-gauge and reference quotient rank, continuity, delay causality and unique loss ownership; laws: partial elements agree with independent dense high-precision baselines, charge/current continuity and complex-power closure hold, complete gauge transforms preserve observables, and reduction retains passivity/ports/error band or refuses; shrinkers: remove cells, basis functions, ports and reduction states while preserving the first quotient/gauge/power fault; replay seeds per fixture"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-peec-analytic",
                "i14-peec-network-core-holdout",
            ],
            g3_relations: &[
                "rigid transform preserves partial elements and terminal observables within geometry/quadrature budgets",
                "complete tree-cotree/reference change preserves physical cochains, terminal QoIs and power",
                "a pinned nested stable refinement family tightens its certified discretization envelope or records the exact nonnested/conditioning refusal rather than assuming stepwise monotonicity",
                "passivity-preserving congruence and port permutation preserve the reduced network supply inequality",
                "short-line PEEC-to-MTL crosswalk agrees only inside the independently admitted overlap region",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during geometry integration, dense assembly, delay evaluation, quotient solve, model reduction, oracle and holdout tiles; drain quadrature/solver/reduction branches, finalize one loss/power receipt and immutable failure bundle, save resumable Krylov/reduction state, then checkpoint, migrate, resume and fork; resumed and uninterrupted network, projection and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_peec_network.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-peec-network dsr quality --tool frankensim",
            obs_events: &[
                "peec.extracted",
                "peec.gauge_checked",
                "peec.power_audited",
                "peec.reduction_checked",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_peec_network.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-ground-shield-core",
            claims_covered: &["i14-ground-bond-shield-current-closure"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: bond/ground/shield/backshell/chassis networks, contact and shield states, DC-to-RF frequency sweeps, and missing-return/double-loss/orientation twins; validity predicates: typed orientation/reference, frequency/material/temperature support, shield surface and aperture semantics, unique conductor/dielectric/contact loss owners, and a distinct radiation-boundary-flux owner; laws: KCL/charge and complex/transient power close, analytic shield-transfer cases agree, isolated components do not perturb connected QoIs, and no current, loss, or radiation flux appears without an owner; shrinkers: prune bonds/shields/apertures while retaining the first closure fault; replay seeds per fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-ground-bond-shield",
                "i14-shield-ground-core-holdout",
            ],
            g3_relations: &[
                "consistent port/orientation reversal changes signed currents but preserves scalar loss and closure",
                "adding an electrically isolated ground component cannot alter admitted connected-component QoIs",
                "360-degree termination to pigtail degradation cannot be presented as monotone outside the deck's validity band",
                "moving a uniquely owned loss between reporting groups preserves total power but invalidates the ownership receipt",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during frequency sweep, nonlinear contact solve, shield/path assembly, oracle and holdout tiles; drain all transient/speculative paths, finalize one terminal closure/loss receipt, retain current and loss ledgers, then checkpoint, migrate, resume and fork at stable time/frequency frontiers; resumed and uninterrupted verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_ground_shield_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-ground-shield-core dsr quality --tool frankensim",
            obs_events: &[
                "ground.current_closed",
                "shield.transfer_checked",
                "loss.ownership_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_ground_shield_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-bearing-current-core",
            claims_covered: &["i14-bearing-current-hybrid-path"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: motor-drive common-mode/bearing paths, capacitance/impedance cards, lubricant-film histories, nonlinear discharge/recovery events and missing-path/energy-creation twins; validity predicates: typed orientation/reference, source/material/temperature/speed/load support, event ordering, unique dielectric/contact/arc loss owners and explicit Unknown for unresolved roots; laws: charge and transient power close, analytic capacitive-divider and reduced hybrid cases agree, open/insulated transforms remove only declared paths, and event refinement preserves ordering outside simultaneous ambiguity; shrinkers: prune capacitances, contacts and event modes while retaining the first path or energy fault; replay seeds per fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-motor-bearing-current",
                "i14-bearing-current-core-holdout",
            ],
            g3_relations: &[
                "consistent shaft/bearing orientation reversal changes signed voltage/current but preserves dissipated event work",
                "removing an admitted parallel return path cannot increase that path's current contribution without a classified mode change",
                "time-step refinement preserves event ordering away from declared simultaneous-event ambiguity and tightens charge/power residual",
                "relabeling bearing occurrences preserves aggregate results only when every state/history and orientation follows the occurrence map",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during common-mode solve, switching/event localization, bearing-discharge classification, oracle and holdout tiles; drain all event/speculative paths, finalize one path/event/energy receipt without partial authority, retain hybrid states and event ledger, then checkpoint, migrate, resume and fork at stable event frontiers; resumed and uninterrupted event/verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_bearing_current_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-bearing-current-core dsr quality --tool frankensim",
            obs_events: &[
                "bearing.path_closed",
                "bearing.event_classified",
                "bearing.energy_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_bearing_current_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-source-probe-core",
            claims_covered: &["i14-switching-source-probe-semantics"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: PWM/common/differential sources and probe transfer/loading/calibration/window chains with alias, saturation, phase, reference and normalization twins; validity predicates: EmConventionCard, RMS/peak and time/frequency normalization, source impedance/support, probe placement/orientation/loading, bandwidth and calibration uncertainty; laws: analytic pulses/spectra and an independent instrument pipeline agree, calibrated equivalent probes agree after transfer correction, probe energy is not invented, and out-of-band exposure refuses; shrinkers: remove samples, source terms and probe stages while preserving the first semantic/power fault; replay seeds follow the source/probe fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-switching-source-probe",
                "i14-source-probe-core-holdout",
            ],
            g3_relations: &[
                "time shift with the consistently shifted window preserves shift-invariant spectral QoIs",
                "sampling-rate increase with compatible anti-alias filtering preserves in-band calibrated QoIs",
                "symbolic RMS-to-peak transport uses the exact algebraic alpha and leaves physical power invariant",
                "inserting an ideal zero-loading identity probe preserves exposure and lineage while changing artifact identity",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during source generation, convolution/sampling, spectrum estimation, transfer correction, oracle and source-heldout tiles; drain instrument tasks, finalize exactly one calibrated exposure/power receipt, retain saturation/alias failures, then checkpoint, migrate, resume and fork at stable sample/window frontiers; resumed and uninterrupted spectra and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_source_probe_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-source-probe-core dsr quality --tool frankensim",
            obs_events: &[
                "source.waveform_bound",
                "probe.calibration_applied",
                "probe.power_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_source_probe_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-victim-mode-core",
            claims_covered: &["i14-immunity-victim-mode-ledger"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: typed injection/exposure records, victim susceptibility states, uncertain thresholds, dwell/debounce, hysteresis, simultaneous-event and recovery/attribution twins; validity predicates: green typed exposure, EmConventionCard, victim schema/state/version, event order, telemetry, evidence color and intervention assumptions; laws: source-to-exposure-to-mode lineage is complete, event classification is stable away from declared boundaries, relative phase and event order remain physical, and ambiguous/censored cases stay PossibleEvent or Unknown; shrinkers: remove injection edges and victim modes while preserving the first event/attribution fault; replay seeds follow the victim fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-immunity-victim",
                "i14-victim-upset-core-holdout",
            ],
            g3_relations: &[
                "victim mode relabeling preserves outcomes only when every ledger edge and classifier label transforms consistently",
                "increasing temporal resolution preserves event order away from a declared simultaneous/grazing boundary",
                "removing required telemetry weakens the result to Unknown and never manufactures NoUpset",
                "a counterfactual intervention may change causal attribution while preserving the observed exposure trace",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during injection binding, victim event localization, mode attribution, recovery, oracle and victim-heldout tiles; drain every classifier/event task, finalize one ModeLedger receipt, retain censored/ambiguous failures, then checkpoint, migrate, resume and fork at stable event frontiers; resumed and uninterrupted mode-ledger and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_victim_mode_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-victim-mode-core dsr quality --tool frankensim",
            obs_events: &[
                "victim.mode_classified",
                "victim.outcome_unknown",
                "victim.recovery_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_victim_mode_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fidelity-routing-core",
            claims_covered: &["i14-core-fidelity-crosswalk-routing"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: fixed-rung and MTL/PEEC overlap-region cases, competing routes, validity/topology/reference/conditioning negative twins and explicit crosswalk diagrams; validity predicates: common physical port/QoI/reference/phasor/loss semantics, admitted fidelity region, typed error budget, route hysteresis and evidence preorder; laws: deterministic tie-break prevents chatter, independent fine-rung checks bound declared routed QoIs, missing ownership or Unknown eligibility never selects the cheapest route, and MTL/PEEC crossover failures block only routing authority; shrinkers: simplify route graph, ports and overlap domains while preserving the first admission fault; replay seeds per decks"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-core-crosswalk-ladder",
                "i14-peec-mtl-crossover-core-holdout",
            ],
            g3_relations: &[
                "equivalent common-port basis changes preserve routed physical QoIs",
                "tightening a valid error budget cannot make an inadmissible route admissible",
                "reassociating an admitted route preserves the final receipt only when dependency and signed loss/error ownership are retained",
                "permuting candidate enumeration preserves the selected route under deterministic tie-break",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during route admission/race, crosswalk, escalation, oracle and holdout tiles; cancel and fully drain losing routes, finalize one authority-colored routing receipt, save route hysteresis and error-owner state, then checkpoint, migrate, resume and fork; resumed and uninterrupted selections and digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fidelity_routing_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fidelity-routing-core dsr quality --tool frankensim",
            obs_events: &[
                "fidelity.route_admitted",
                "fidelity.escalated",
                "fidelity.crosswalk_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fidelity_routing_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fixed-regime-adjoint-core",
            claims_covered: &["i14-fixed-regime-adjoint-closure"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: fixed-rung complex adjoint studies, parameter transforms, inexact solves, nonnormal operators, event/topology boundary refusals and wrong-conjugate/trace negative twins; validity predicates: fixed topology/mesh/active set/event/routing regime, declared real or Wirtinger convention, parameter units, primal/adjoint residual budgets and independent derivative route; laws: tangent-adjoint duality, Taylor remainder, unit covariance and independent directional differences close only inside the frozen regime, while switches return NoGradient/SetValued/Unknown; shrinkers: reduce states, parameters and residual terms while retaining the first derivative fault; replay seeds per decks"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-adjoint-uq-mitigation",
                "i14-adjoint-core-holdout",
            ],
            g3_relations: &[
                "parameter unit rescaling transforms gradients inversely and leaves predicted physical delta invariant",
                "equivalent common-port basis changes transform adjoints covariantly",
                "tightening primal or adjoint solve tolerance cannot excuse a larger independently recomputed derivative defect",
                "crossing a topology, event, active-set or route boundary converts classical authority to NoGradient or SetValued",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during primal/tangent/adjoint solve, independent difference, Taylor sweep, oracle and holdout tiles; drain every solver child, finalize one derivative receipt with no partial gradient authority, save Krylov/tape/parameter/regime state, then checkpoint, migrate, resume and fork; resumed and uninterrupted derivative vectors and digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fixed_regime_adjoint_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fixed-regime-adjoint-core dsr quality --tool frankensim",
            obs_events: &[
                "adjoint.regime_admitted",
                "adjoint.checked",
                "adjoint.no_gradient",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fixed_regime_adjoint_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-uq-inference-core",
            claims_covered: &["i14-emc-uq-inference-mechanics"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: aleatory/epistemic/dependence models, bounded and sampled reliability studies, correlated tails, censoring/missingness/optional-stopping cases and miscoverage/public-holdout authority twins; validity predicates: population or bounded-set scope, dependence, weights, censoring, multiplicity/stopping, QoI direction, escalation and evidence color; laws: sample/worker permutation preserves deterministic estimates, dependence cannot be broken to manufacture evidence, public replay never gains population authority, and certify-or-escalate/refuse is monotone under lost evidence; shrinkers: reduce variables, samples and tail events while preserving the first inference fault; replay seeds per decks"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-adjoint-uq-mitigation",
                "i14-uq-core-holdout",
            ],
            g3_relations: &[
                "permuting uncertainty samples and deterministic worker assignment preserves estimates and decisions",
                "splitting correlated uncertainty into duplicate variables without preserving dependence is detected rather than credited as more evidence",
                "widening an epistemic set cannot improve a worst-case guard without a classified applicability change",
                "removing population or custody authority turns the corresponding conclusion Unknown without erasing deterministic replay evidence",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during sampling/bounding, tail estimation, escalation, oracle and holdout tiles; drain all estimator children, finalize one evidence-colored receipt without optional-stopping reuse, save counters/weights/bounds/stopping state, then checkpoint, migrate, resume and fork; resumed and uninterrupted estimates, decisions and digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_uq_inference_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-uq-inference-core dsr quality --tool frankensim",
            obs_events: &[
                "uq.authority_routed",
                "uq.coverage_checked",
                "uq.escalated",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_uq_inference_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fullwave-schema-core",
            claims_covered: &["i14-fullwave-problem-convention-admission"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: typed FullWaveProblem/Receipt convention and crosswalk records plus missing/conflicting sign, source, RMS/peak, field-degree, boundary, material, gauge, topology and loss-owner twins; validity predicates: exact units/phasor/Laplace/source/port/boundary/material/loss semantics and versioned crosswalk applicability; laws: schema refuses every ambiguity, exact crosswalk roundtrips preserve real fields and physical power, and admission never mints solver or physical authority; shrinkers: delete fields and boundaries while retaining the first semantic conflict; replay seeds per fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-fullwave-schema-crosswalk",
                "i14-fullwave-schema-core-holdout",
            ],
            g3_relations: &[
                "consistent orientation and conjugate phasor transform preserves real physical fields and reverses only declared signed complex quantities",
                "RMS-to-peak conversion followed by its inverse preserves terminal and Poynting power with the exact one-half factor",
                "renaming topology/gauge representatives preserves admission only when every bound source, boundary, field degree and loss owner is transported",
                "removing a required semantic field can only preserve Refused/Unknown or weaken admission, never manufacture success",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during schema decoding, crosswalk, independent reconstruction, oracle and holdout tiles; drain all checker children, finalize one admission receipt with no PDE authority, persist canonical problem/crosswalk state, then checkpoint, migrate, resume and fork; resumed and uninterrupted schema/verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fullwave_schema_core.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fullwave-schema-core dsr quality --tool frankensim",
            obs_events: &[
                "fullwave.problem_admitted",
                "fullwave.crosswalk_checked",
                "fullwave.problem_refused",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fullwave_schema_core.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fullwave-feec-max",
            claims_covered: &["i14-fullwave-feec-stability-energy"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: cavity/waveguide/manufactured fields, multiply connected complexes, dispersive/high-contrast media, PML/radiation cases and hostile spurious-mode/gauge/resonance/late-time twins; validity predicates: green fullwave schema, de Rham spaces and topology/nullspaces, inf-sup/coercivity or declared stabilization premise, passive causal dispersion, PML/radiation and independent error decomposition; laws: exact discrete sequence is necessary but never treated as sufficient stability, manufactured/convergence order and Gauss/energy/complex-power balances close, and PML remains stable over the declared late-time window; shrinkers: coarsen complex/material/port/boundary and basis while retaining the first stability fault; replay seeds per fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-fullwave-cavity-waveguide",
                "i14-fullwave-pml-dispersion",
                "i14-fullwave-max-holdout",
            ],
            g3_relations: &[
                "rigid transform preserves admitted eigen and field QoIs with transformed ports and fields",
                "h/p/time refinement follows the manufactured convergence band and separates geometry, discretization, algebraic and time errors",
                "PML thickness/stretch sweep exposes nonmonotone or unstable regimes rather than asserting universal monotonic improvement",
                "equivalent harmonic-basis and gauge representatives preserve complete physical fields, sources and power",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during complex/mesh construction, assembly/apply, eigensolve/time integration, PML monitoring, oracle and maximal-heldout tiles; drain all kernels and speculative formulations, finalize one stability/energy/loss receipt with no partial authority, persist mesh/Krylov/time/material state, then checkpoint, migrate, resume and fork at stable tile/frontier boundaries; resumed and uninterrupted verdict/evidence digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fullwave_feec_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fullwave-feec-max dsr quality --tool frankensim",
            obs_events: &[
                "fullwave.stability_checked",
                "fullwave.energy_audited",
                "pml.late_time_checked",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fullwave_feec_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-bem-formulation-max",
            claims_covered: &["i14-exterior-bem-formulation-correctness"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: separately typed closed-PEC, closed-dielectric-transmission and admitted open-screen scatterers, resonance/low-frequency/near-touching/edge and wrong-formulation twins; validity predicates: green fullwave schema and EmConventionCard, formulation-specific trace and regularity spaces, orientation, singular quadrature, resonance/conditioning route, FEEC coupling and independent dense error allocation; laws: each class refuses incompatible formulations and dense small-problem matrices, solved QoIs, trace and power are the correctness oracle; shrinkers: simplify surface, material sides, edge and basis while retaining the first trace/formulation fault; replay seeds per fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-fullwave-bem-scattering",
                "i14-bem-formulation-max-holdout",
            ],
            g3_relations: &[
                "rigid transform preserves admitted scattering QoIs with transformed currents, normals and fields",
                "swapping dielectric sides with the complete material/normal/trace transformation preserves the physical solution",
                "an open-screen case cannot be reclassified as a closed PEC or transmission problem by adding a zero-area closure",
                "a change of equivalent supported trace basis preserves physical currents and fields after the complete dual transport",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during surface/basis construction, singular quadrature, dense assembly/solve, oracle and formulation-heldout tiles; drain all kernels and speculative formulations, finalize one trace/power/dense receipt with no partial authority, persist surface/Krylov state, then checkpoint, migrate, resume and fork at stable tile/frontier boundaries; resumed and uninterrupted verdict/evidence digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_bem_formulation_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-bem-formulation-max dsr quality --tool frankensim",
            obs_events: &[
                "bem.formulation_admitted",
                "bem.dense_checked",
                "bem.power_audited",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_bem_formulation_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fmm-acceleration-max",
            claims_covered: &["i14-maxwell-fmm-acceleration-envelope"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: formulation-bound dense/FMM operator pairs, hostile trees, near/far boundaries, expansion orders, vector/tangential kernels, adjoints, route changes, cancellation and crossover twins; validity predicates: green dense BEM formulation, EmConventionCard, exact kernel/tree/order/tolerance and independent dense references with separated error owners; laws: every accelerated matvec and solved QoI lies in its certified dense envelope, deterministic traversal is stable, and a tighter request never grants pointwise monotonicity across algorithm changes; shrinkers: reduce panels, tree levels, orders and QoIs while retaining the first envelope or traversal fault; replay seeds per FMM fixtures"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-fullwave-bem-scattering",
                "i14-fmm-dense-envelope-adversaries",
                "i14-fmm-acceleration-max-holdout",
            ],
            g3_relations: &[
                "rigid transform preserves accelerated scattering QoIs with transformed currents, normals and fields",
                "permuting input panel order preserves the deterministic canonical tree and physical result",
                "tightening the requested tolerance cannot widen the declared certified envelope, but observed point error may be nonmonotone and every route is checked independently",
                "switching across the measured dense/FMM crossover changes the route receipt and cannot inherit performance authority from the other side",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during tree construction, near/far partition, expansion/apply, adjoint, dense oracle, performance measurement and FMM-heldout tiles; drain all kernels and losing route tasks, finalize one dense-envelope/crossover receipt, persist tree/order/frontier state, then checkpoint, migrate, resume and fork at stable tree/frontier boundaries; resumed and uninterrupted tree, matvec, QoI and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fmm_acceleration_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fmm-acceleration-max dsr quality --tool frankensim",
            obs_events: &[
                "fmm.tree_bound",
                "fmm.dense_envelope_checked",
                "fmm.crossover_measured",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fmm_acceleration_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-fidelity-descent-theorem-max",
            claims_covered: &["i14-certified-fidelity-descent-theorem"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: maximal cross-rung ladders, theorem target instances and resonance/topology/asymptotic/dependence negative twins; validity predicates: typed common QoI/port/topology/reference/loss maps, signed dependency-aware error budgets, theorem binding, runtime premises, axiom closure and nonvacuity; laws: the formal/runtime checker either encloses the fine QoI or refuses, invalid premise cases cannot inherit bounds, evidence never strengthens, and a genuine countermodel refutes only its exact theorem revision; shrinkers: minimize rung graph, error owner and theorem premise while retaining the first counterexample; replay per theorem policy"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-core-crosswalk-ladder",
                "i14-fullwave-crosswalk-max",
                "i14-fidelity-descent-theorem-card",
                THEOREM_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "refining the fine rung while preserving all descent hypotheses cannot widen the independently recomputed fine discretization component without a reasoned failure",
                "rung-path reassociation preserves the final enclosure only when dependency and signed-error ownership are retained",
                "equivalent reference/gauge/port representation changes preserve the physical enclosure after complete map transport",
                "removing a theorem premise can only refuse the runtime instance, never preserve favorable authority by default",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during AST translation, proof/kernel checking, runtime-premise admission, crosswalk/majorant assembly and oracle tiles; drain every proof/checker child, finalize one theorem/runtime receipt, retain declarations and premise failures, then checkpoint, migrate, resume and fork at canonical proof/graph frontiers; resumed and uninterrupted declarations, bounds and digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_fidelity_descent_theorem_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-fidelity-descent-theorem-max dsr quality --tool frankensim",
            obs_events: &[
                "fidelity.theorem_bound",
                "fidelity.descent_checked",
                "theorem.binding_checked",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_fidelity_descent_theorem_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-robust-mitigation-max",
            claims_covered: &["i14-robust-mitigation-heldout"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: constrained mitigation candidates, fixed catalog/ablation baselines, inaccessible heldout slots and synthetic pre-governance access, stage reordering, adaptive realization, partial discharge, prereveal, retry, objective/guard/leakage/multiplicity twins; validity predicates: GovernanceCommitted generator/acquisition/access protocol and AcceptanceCards, immutable CandidateFrozen candidate/model/toolchain/checker roots, independent custodian realization, typed slot replacement plus same-ID discharge envelope and verified receipt in one atomic RealizationCommitted authority-head transaction, RevealedForAdjudication one-shot access, complete adjudication/access/attempt ledger and Closed receipt; laws: optimization and every candidate-side transitive capability never read a protected root, byte, label, aggregate or derived statistic before CandidateFrozen and committed realization, the custodian follows the frozen non-adaptive generator/selection mechanism, one immutable candidate must improve every directed objective and satisfy every guard or remain Unknown/Failed, partial/split retirement grants no authority, and public replay never acquires blind authority; shrinkers: minimize lifecycle stages, principals/capabilities, design variables, uncertainty and guard sets while retaining the first authority leak, adaptation or regression; replay/custody uses every governed stage receipt"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-adjoint-uq-mitigation",
                "i14-mitigation-max-holdout",
                "i14-external-blind-mitigation-custody-pack",
            ],
            g3_relations: &[
                "renaming/reordering feasible design variables preserves the selected physical candidate under deterministic tie-break",
                "tightening a safety or fabrication constraint cannot admit a previously infeasible candidate",
                "adding a required directed guard can only preserve or weaken promotion, never improve it",
                "revealing any heldout root, label or aggregate before the joint commitment makes the campaign IntegrityFailed",
                "reordering GovernanceCommitted, CandidateFrozen, custodian realization, RealizationCommitted, RevealedForAdjudication and Closed changes the campaign to IntegrityFailed",
                "changing the frozen generator, selection mechanism, candidate, checker, AcceptanceCard, retry, exclusion, multiplicity or stopping rule requires a successor and cannot preserve blind authority",
                "slot-only replacement, same-ID-envelope-only installation, waiver-only retirement or a split authority-head transition cannot preserve the atomic realization receipt",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during GovernanceCommitted freeze, design search, adjoint/UQ evaluation, fidelity escalation, CandidateFrozen sealing, independent custodian realization, slot/envelope/receipt verification, atomic RealizationCommitted authority-head commit, RevealedForAdjudication access, heldout adjudication, guard checking and Closed finalization; drain every losing design/rung/oracle/custodian child, forbid post-cancel protected access, reveal and authority-head advancement unless their exact stage transaction had already committed, finalize a single candidate and authority-colored receipt without threshold/holdout reuse, retain all attempted candidates, access/exclusion/retry/multiplicity/stopping decisions, protected-root commitments and failures, then checkpoint, migrate, resume and fork at every stable lifecycle/search/adjudication frontier; resumed and uninterrupted stage, authority-head, selected-candidate and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_robust_mitigation_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-robust-mitigation-max dsr quality --tool frankensim",
            obs_events: &[
                "mitigation.governance_committed",
                "mitigation.candidate_frozen",
                "mitigation.custodian_realized",
                "mitigation.authority_committed",
                "mitigation.revealed_for_adjudication",
                "mitigation.holdout_adjudicated",
                "mitigation.guard_checked",
                "mitigation.closed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_robust_mitigation_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-safety-case-integration-max",
            claims_covered: &["i14-emc-safety-case-integration"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: synthetic source-to-victim-to-hazard assurance graphs with stale premise, missing owner, common-cause, evidence-color, scope, expiry and simulation-as-compliance twins; validity predicates: exact artifact revisions, monotone evidence colors, source-to-victim causal assumptions, HazardId owners/monitors/expiry and explicit Unknown external edges; laws: graph traversal preserves owners/scope/expiry, invalid/refuted evidence dominates favorable observations, removing a premise makes every dependent conclusion Unknown, and simulation never becomes compliance or safety approval; shrinkers: minimize assurance nodes/edges while retaining the first laundering defect; replay seeds per deck"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-safety-case-source-victim",
            ],
            g3_relations: &[
                "removing or expiring an assurance premise turns dependent conclusions Unknown and never leaves a green compliance edge",
                "renaming graph nodes preserves the verdict only with every content-addressed owner and edge transported",
                "adding a weaker evidence edge cannot upgrade a stronger unresolved or refuted dependency",
                "external standards/laboratory/population edges may remain Unknown while synthetic traceability independently promotes",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during assurance reconstruction, causal-edge audit, dependency traversal and adjudication tiles; drain every traversal/checker child, finalize one scoped traceability receipt with external edges explicitly Unknown, retain all missing/stale edges, then checkpoint, migrate, resume and fork at stable graph frontiers; resumed and uninterrupted assurance graph and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_safety_case_integration_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-safety-case-integration-max dsr quality --tool frankensim",
            obs_events: &[
                "safety.edge_audited",
                "safety.laundering_refused",
                "authority.unknown_routed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_safety_case_integration_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-passive-composition-theorem-max",
            claims_covered: &["i14-passive-causal-sheaf-composition-theorem"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: strength-matched finite rational compatibility-sheaf/balance-cosheaf assemblies, relative-complex boundaries and winding/cut representatives, local causal dissipative relations, maximal-isotropic Dirac interconnections in a nondegenerate ambient split pairing, well/ill-posed feedback and gluing/time/solver defects; validity predicates: exact policy-bound proposition/definition/runtime-premise identity and axiom closure, typed signal/history/storage/supply spaces, pairwise and higher-overlap coherence, relative-boundary compatibility, chain-homotopy representative maps, owned losses/sources and dependency-aware margin debit; laws: kernel/rebinding/nonvacuity/admission all agree before theorem authority, zero unaccounted defect is required for exact loss identity, integral dissipation greater than or equal to the accounted violation yields passivity, a frozen positive signal margin yields robust strict passivity, and representative transforms preserve physical observables; shrinkers: premise-preserving diagram/component/topology minimization; replay per theorem policy"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-passive-composition-theorem-card",
                THEOREM_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "cover refinement followed by canonical reaggregation preserves the global inequality when every pairwise/higher-overlap map and defect budget is transported",
                "complete chain-homotopic gauge/tree-cotree/cut representative change preserves physical cochains, supply and theorem verdict",
                "adding a disconnected zero-power component preserves the original verdict but changes the bound assembly identity",
                "for an identical admitted trajectory, increasing an owned nonnegative dissipation term cannot reduce the certified margin, while increasing the violation bound cannot improve it",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during AST translation/elaboration, axiom audit, runtime-premise admission, proof/nonvacuity and diagram-checker tiles; drain proof/checker children, finalize one theorem receipt with every premise state, persist declaration/proof identities, then checkpoint, migrate, resume and fork at every canonical frontier; resumed and uninterrupted theorem and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_passive_composition_theorem_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-passive-composition-theorem-max dsr quality --tool frankensim",
            obs_events: &[
                "theorem.binding_checked",
                "theorem.axioms_audited",
                "theorem.nonvacuity_checked",
                "theorem.runtime_premises_checked",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_passive_composition_theorem_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-hypercohomology-obstruction-max",
            claims_covered: &["i14-hypercohomology-obstruction-localization-theorem"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: canonical finite relative total complexes, torsion/field coefficient cases, gluable/non-gluable local data, nonzero obstruction classes, false-vanishing and tied-minimal-support twins; validity predicates: exact site/cover, boundary, coefficient policy, maps, total signs, theorem binding, runtime admission and nonvacuity; laws: d_total squared is zero, nonzero admitted class blocks gluing, vanishing yields a witness only under the exact sufficiency premises, and every tied minimal support is retained; shrinkers: premise-preserving complex, overlap and support minimization; replay per theorem policy"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-passive-composition-theorem-card",
                "i14-hypercohomology-obstruction-card",
                THEOREM_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "chain-homotopy-equivalent representatives preserve the obstruction class and physical observable",
                "adding a contractible zero-data patch preserves the original class after canonical comparison",
                "coefficient change never transports a vanishing or torsion verdict without the exact base-change theorem",
                "permuting tied minimal repair supports changes only deterministic presentation order",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during total-complex construction, exact reduction, proof/kernel checking, runtime admission, support minimization and oracle tiles; drain every algebra/proof/checker child, finalize one obstruction/witness/support receipt, persist reduction and proof frontiers, then checkpoint, migrate, resume and fork at canonical frontiers; resumed and uninterrupted classes, supports and digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_hypercohomology_obstruction_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-hypercohomology-obstruction-max dsr quality --tool frankensim",
            obs_events: &[
                "hypercohomology.complex_bound",
                "hypercohomology.obstruction_checked",
                "hypercohomology.support_minimized",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_hypercohomology_obstruction_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-cover-refinement-naturality-max",
            claims_covered: &["i14-cover-refinement-naturality-theorem"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: pairs of covers/sites, common refinements, valid and invalid cofinal/Leray/quasi-isomorphic comparisons, relative boundaries, port/power and duplicated-defect twins; validity predicates: exact refinement functor, comparison maps/homotopies, coefficient policy, theorem binding, runtime admission and no duplicated owner; laws: admitted comparison preserves sections, obstruction classes, boundary supply and theorem verdict, while invalid coarsening returns a comparison obstruction; shrinkers: minimize covers, overlaps and owner maps while retaining the first naturality fault; replay per theorem policy"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-passive-composition-theorem-card",
                "i14-cover-refinement-naturality-card",
                THEOREM_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "refinement followed by canonical reaggregation preserves admitted global sections and supply",
                "two refinement paths through a common refinement agree only when the frozen comparison coherence holds",
                "adding a redundant acyclic patch preserves authority through its explicit contraction",
                "invalid coarsening weakens authority to a comparison obstruction and never silently inherits the fine verdict",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during cover/nerve construction, comparison-map assembly, proof/kernel checking, runtime admission and diagram oracle tiles; drain every topology/proof/checker child, finalize one naturality/comparison receipt, persist cover/map/proof frontiers, then checkpoint, migrate, resume and fork at canonical frontiers; resumed and uninterrupted comparison and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_cover_refinement_naturality_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-cover-refinement-naturality-max dsr quality --tool frankensim",
            obs_events: &[
                "refinement.comparison_bound",
                "refinement.naturality_checked",
                "refinement.obstruction_reported",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_cover_refinement_naturality_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-kyp-sheaf-bridge-max",
            claims_covered: &["i14-kyp-sheaf-passivity-bridge-theorem"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: exact rational/state-space and descriptor generalized-PR components, simple boundary poles, quotient/lossless modes, KYP witnesses and irregular/impulsive/indefinite/wrong-supply twins; validity predicates: EmConventionCard, regular impulse-free descriptor class, consistent initialization, exact generalized-PR/strict-PR distinction, theorem binding and runtime relation identity; laws: the checked witness yields exactly the local storage/supply relation, boundary storage poles remain typed, and sampled-band evidence refuses; shrinkers: reduce states, constraints, poles and ports while retaining the first bridge fault; replay per theorem policy"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-rlgc-fit-adversaries",
                "i14-kyp-sheaf-bridge-card",
                THEOREM_POLICY_FIXTURE,
            ],
            g3_relations: &[
                "an exact similarity transform transports storage and supply witnesses without changing the physical relation",
                "complete port congruence preserves the KYP inequality and oriented supply",
                "adding an unreachable lossless mode preserves external behavior but changes the nonminimal witness identity",
                "deleting the descriptor consistency premise forces refusal rather than strict-stability promotion",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during realization canonicalization, exact KYP checking, descriptor admission, proof/kernel checking and local-relation oracle tiles; drain every algebra/proof/checker child, finalize one witness/bridge receipt, persist exact-factorization and proof frontiers, then checkpoint, migrate, resume and fork at canonical frontiers; resumed and uninterrupted witness, relation and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_kyp_sheaf_bridge_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-kyp-sheaf-bridge-max dsr quality --tool frankensim",
            obs_events: &[
                "kyp.realization_admitted",
                "kyp.witness_checked",
                "kyp.local_relation_bound",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_kyp_sheaf_bridge_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-maximal-falsifier-max",
            claims_covered: &["i14-maximal-counterexample-search"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: formal theorem ASTs, executable finite rational microgrammar candidates and separately seeded cross-domain adversaries over RLGC, PML, BEM/FMM, routing, mitigation and assurance; validity predicates: exact target-card identity, theorem-premise admission, executable grammar membership/rank/cardinality, independent classification and nonvacuity floors; laws: rank/unrank and independent cardinality/Merkle roots close, every candidate receives one reason-coded classification, and only GenuineCountermodel refutes its exact immutable revision; shrinkers: premise-preserving grammar/candidate/countermodel minimization; replay seeds per falsifier range"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-passive-composition-theorem-card",
                "i14-hypercohomology-obstruction-card",
                "i14-cover-refinement-naturality-card",
                "i14-kyp-sheaf-bridge-card",
                "i14-fidelity-descent-theorem-card",
                THEOREM_POLICY_FIXTURE,
                "i14-theorem-falsifier-grammar",
                "i14-maximal-adversaries",
                "i14-rlgc-fit-adversaries",
                "i14-fullwave-pml-dispersion",
                "i14-fullwave-bem-scattering",
            ],
            g3_relations: &[
                "rank then unrank is identity for every exhaustive microgrammar candidate",
                "shard permutation preserves the independently recomputed Merkle completeness root",
                "adding out-of-domain or checker-defect candidates cannot be counted as theorem countermodels",
                "minimization preserves the exact target revision, admission premises and GenuineCountermodel classification",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during grammar enumeration/sharding, theorem-premise admission, independent classification, countermodel minimization and artifact assembly; drain search/checker children, finalize each candidate as a reason-coded state and one terminal campaign receipt, persist enumeration rank/Merkle frontier and target identities, then checkpoint, migrate, resume and fork at every canonical frontier; resumed and uninterrupted roots, classifications, minimized countermodels and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_maximal_falsifier_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-maximal-falsifier-max dsr quality --tool frankensim",
            obs_events: &[
                "falsifier.grammar_bound",
                "falsifier.candidate_classified",
                "falsifier.countermodel_minimized",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_maximal_falsifier_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-standards-crosswalk-max",
            claims_covered: &["i14-governed-standards-crosswalk"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: governed exact-edition clause graphs plus legally authored AP242/EMC/machine synthetic pre-governance access, unauthorized principal, stage reordering, edition shopping, omitted scope, stale-edition, occurrence/unit/frame, applicability, configuration, restricted-text leakage, partial discharge and favorable-Unknown twins; validity predicates: GovernanceCommitted publisher/edition/corrigendum/license/source-root/scope identities and authorized principals/sandbox/procedure/disclosure filter, least-privilege AuthorizedConstruction read/output ledger, complete derived crosswalk/loss mapping and restricted-text-free projection, CandidateFrozen crosswalk/toolchain/checker/AcceptanceCard roots, same-ID discharge envelope plus verified receipt and atomic StandardsAuthorityCommitted authority-head transaction, independent StandardsAdjudicated reconstruction and Closed receipt; laws: no licensed byte or derived output is accessed before GovernanceCommitted or outside AuthorizedConstruction, every in-scope clause and source construct maps to one supported edge, explicit loss, Unknown, Nonconformant or refusal, CandidateFrozen precedes independent adjudication, partial/split retirement grants no authority, restricted text never enters public artifacts, and cross-edition or schema-shaped public evidence never promotes exact-edition authority; shrinkers: minimize lifecycle stages, principals/access edges, clauses, occurrences and mapping edges while retaining the first authority, omission or disclosure fault; replay uses governed digests and stage receipts and emits no restricted text"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-ap242-adapter-mechanics",
                "i14-standards-crosswalk-adversaries",
                "i14-safety-case-source-victim",
                "i14-external-standards-edition-clause-pack",
            ],
            g3_relations: &[
                "adding a mandatory scoped clause cannot preserve completeness if it lacks an independently adjudicated edge",
                "renumbering clauses or occurrences preserves the verdict only through authenticated identity maps",
                "edition or corrigendum substitution turns dependent edges Unknown until a successor rebinds them",
                "removing licensed authority leaves synthetic adapter mechanics intact but removes exact-edition promotion",
                "reordering GovernanceCommitted, AuthorizedConstruction, CandidateFrozen, StandardsAuthorityCommitted, StandardsAdjudicated and Closed changes the campaign to IntegrityFailed",
                "changing publisher, edition, corrigendum, scope, principal, sandbox, procedure or disclosure filter requires a successor and cannot preserve exact-edition authority",
                "same-ID-envelope-only installation, waiver-only retirement or a split authority-head transition cannot preserve the atomic standards receipt",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during GovernanceCommitted freeze, AuthorizedConstruction access, exact-edition admission, clause/occurrence reconstruction, loss/disclosure audit, CandidateFrozen sealing, same-ID envelope/receipt verification, atomic StandardsAuthorityCommitted authority-head commit, independent StandardsAdjudicated reconstruction and Closed finalization; drain every builder/parser/checker/adjudicator child, forbid post-cancel licensed access, publication, adjudication and authority-head advancement unless their exact stage transaction had already committed, finalize one scoped crosswalk receipt with restricted content represented only by governed digests, retain complete read/output/disclosure ledgers, omissions and refusals, then checkpoint, migrate, resume and fork at every stable lifecycle/graph frontier; resumed and uninterrupted stage, authority-head, mapping and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_standards_crosswalk_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-standards-crosswalk-max dsr quality --tool frankensim",
            obs_events: &[
                "standards.governance_committed",
                "standards.construction_authorized",
                "standards.edition_admitted",
                "standards.clause_mapped",
                "adapter.loss_audited",
                "standards.candidate_frozen",
                "standards.authority_committed",
                "standards.adjudicated",
                "standards.closed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_standards_crosswalk_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-laboratory-validation-max",
            claims_covered: &["i14-governed-laboratory-emc-validation"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: governed as-built specimen/laboratory packages plus synthetic pre-access commitment substitution, unauthorized read, commitment-opening leak, stage reordering, adaptive validation selection, calibration expiry, fixture/chamber drift, contamination, hidden model update, partial multi-pack discharge, validation prereveal, retry, exclusion, covariance, censoring, nonfinite and post-hoc-normalization twins; validity predicates: root-only/no-capability pack-commitment registration, GovernanceCommitted hiding calibration/model-input commitments and validation source-universe/selection commitment, AuthorizedCalibration least-privilege ledger, AsBuiltModelInstantiation derivations, CandidateFrozen candidate/checker/AcceptanceCard roots, independent contamination receipt naming audited principals and transitive capabilities, separately addressable realized-root membership/disjointness/non-adaptive-selection proofs, distinct same-ID discharge envelopes and receipts, atomic RealizationCommitted authority-head advance, validation-only reveal and Closed receipt; laws: authorized lifecycle order is strict, pre-governance registration exposes only opaque hiding commitment identities, no candidate builder/fitter/checker/threshold principal or transitive capability accesses validation bytes, membership witnesses, labels, aggregates, derived statistics, selection output or commitment-opening material before CandidateFrozen and contamination check, every realized stratum proves membership in its pre-access commitment and mutual disjointness, validation follows the committed non-adaptive mechanism, the coupled packs commit atomically before validation reveal, independent trace reconstruction and calibration propagation close, every preregistered component is reported, no missing/nonfinite/failed QoI is dropped, and synthetic solver authority remains orthogonal; shrinkers: minimize lifecycle stages, principals/capabilities, access/derivation edges, commitment proofs, traces, instruments and components while retaining the first authority, integrity or validation fault; replay uses governed raw roots and stage receipts"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-switching-source-probe",
                "i14-immunity-victim",
                "i14-laboratory-validation-adversaries",
                "i14-laboratory-validation-max-holdout",
                "i14-external-emc-laboratory-calibration-pack",
                "i14-external-asbuilt-specimen-geometry-pack",
            ],
            g3_relations: &[
                "permuting storage or worker order preserves componentwise results only when immutable timestamps, channel ids, sequence numbers and acquisition-clock metadata travel with every sample and canonical reconstruction restores the identical physical trace order",
                "removing or expiring a calibration can only turn dependent components Unknown or IntegrityFailed",
                "adding a preregistered QoI cannot improve the all-components validation vector",
                "changing specimen/configuration identity invalidates the scoped result even if numeric traces happen to match",
                "equivalent chunking of a frozen stratum preserves its realized root only through verified membership, mutual-disjointness and canonical-reconstruction proofs",
                "giving any candidate builder/fitter/checker/threshold principal or transitive capability validation bytes, membership witnesses, labels, aggregates, derived statistics, selection output or commitment-opening material before CandidateFrozen and the contamination check changes the result to IntegrityFailed regardless of numeric agreement",
                "changing the validation source universe, committed selection algorithm or pre-candidate secret-seed/VRF commitment requires a successor and cannot preserve physical-validation authority",
                "reordering, splitting or partially committing the two discharge envelopes and receipts cannot preserve the atomic authority-head transition",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during root-only/no-capability opaque pack-commitment registration, GovernanceCommitted commitment freeze, AuthorizedCalibration access, calibration checking, AsBuiltModelInstantiation, CandidateFrozen sealing, contamination checking, validation realization, separately addressable root/proof verification, atomic RealizationCommitted authority-head commit, validation-only reveal, trace reconstruction, component comparison, adjudication and Closed finalization; drain every builder/parser/oracle/custodian child, forbid post-cancel access, reveal and authority-head advancement unless their exact stage transaction had already committed, finalize one stage-complete componentwise validation vector without favorable omission, retain pre-access commitments, audited-principal/capability and access/derivation ledgers, raw roots, membership/disjointness/non-adaptive-selection proofs, both discharge receipts and all failures, then checkpoint, migrate, resume and fork at every stable lifecycle/trace/component frontier; resumed and uninterrupted stage, authority-head, component and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_laboratory_validation_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-laboratory-validation-max dsr quality --tool frankensim",
            obs_events: &[
                "laboratory.pack_commitment_registered",
                "laboratory.governance_committed",
                "laboratory.calibration_authorized",
                "laboratory.calibration_checked",
                "laboratory.model_instantiated",
                "laboratory.candidate_frozen",
                "laboratory.contamination_checked",
                "laboratory.validation_realized",
                "laboratory.authority_committed",
                "laboratory.validation_revealed",
                "laboratory.qoi_adjudicated",
                "laboratory.closed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_laboratory_validation_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-emc-reliability-validation-max",
            claims_covered: &["i14-governed-emc-reliability-validation"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: governed EMC deployment populations plus synthetic pre-freeze access, lifecycle reordering, partial-discharge, post-reveal amendment, frame/configuration/lot shift, duty/environment confounding, dependence, censoring, missingness, calibration, optional-stopping, sparse-tail and estimand twins; validity predicates: GovernanceCommitted population frame, outcome/event, EmConventionCard, metrology, inclusion/exclusion, dependence, missingness/censoring, multiplicity/stopping, AcceptanceCards, checker, candidate-input permissions and access protocol; immutable CandidateFrozen model/toolchain/checker/AcceptanceCard roots; independent custodian realization; typed slot replacement plus same-ID discharge envelope and verified receipt in one atomic RealizationCommitted authority-head transaction; one-shot RevealedForAdjudication and complete access/attempt/adjudication ledger in Closed; laws: no candidate-side principal or transitive capability reads a protected population root, byte, label, aggregate or derived statistic before CandidateFrozen and committed realization, the custodian follows the frozen non-adaptive acquisition/selection mechanism, independent frame/event reconstruction and anytime-valid coverage close, PublicReplay or one specimen never substitutes for population evidence, partial/split retirement grants no authority, and out-of-frame extrapolation refuses; shrinkers: reduce lifecycle stages, principals/capabilities, units, events, strata and covariates while retaining the first authority, integrity or coverage fault; replay uses every governed stage receipt and population root"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-adjoint-uq-mitigation",
                "i14-emc-reliability-adversaries",
                "i14-emc-reliability-max-holdout",
                "i14-external-emc-reliability-population-pack",
            ],
            g3_relations: &[
                "permuting unit records and deterministic workers preserves keyed estimates and coverage decisions",
                "adding censored, missing or competing outcomes cannot be treated as event-free exposure",
                "narrowing or shifting the target population requires a new estimand/frame binding",
                "reordering GovernanceCommitted, CandidateFrozen, custodian realization, RealizationCommitted, RevealedForAdjudication and Closed changes the campaign to IntegrityFailed",
                "changing the frozen population frame, acquisition/selection mechanism, candidate, checker, AcceptanceCard, stopping rule or estimand requires a successor before protected access",
                "slot-only replacement, same-ID-envelope-only installation, waiver-only retirement, split authority-head transition, sequential reveal or post-reveal amendment grants no population authority",
                "removing governed population authority leaves Core inference mechanics intact but reliability Unknown",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during GovernanceCommitted freeze, CandidateFrozen sealing, independent custodian realization, same-ID envelope/receipt verification, atomic RealizationCommitted authority-head commit, RevealedForAdjudication access, frame/event reconstruction, survival/e-process/coverage computation, adjudication and Closed finalization; drain every candidate/custodian/estimator/oracle child, forbid post-cancel protected access, reveal, adjudication and authority-head advancement unless the exact atomic stage had already committed, finalize one population-scoped component vector with no favorable missingness, retain complete frame/event/censoring/access/attempt/stopping ledgers, protected-root commitment, receipt and every failure, then checkpoint, migrate, resume and fork at every stable lifecycle/population frontier; resumed and uninterrupted stage, authority-head, estimates, coverage and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_emc_reliability_validation_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-emc-reliability-validation-max dsr quality --tool frankensim",
            obs_events: &[
                "emc.governance_committed",
                "emc.candidate_frozen",
                "emc.custodian_realized",
                "emc.authority_committed",
                "emc.revealed_for_adjudication",
                "emc.reliability_event_audited",
                "emc.reliability_coverage_checked",
                "emc.reliability_adjudicated",
                "emc.closed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_emc_reliability_validation_max.sh --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i14-bearing-population-validation-max",
            claims_covered: &["i14-production-bearing-population-reliability"],
            unit_cases: UNIT_CASES,
            g0: g0_contract!(
                "generators: governed production bearing populations plus synthetic pre-freeze access, lifecycle reordering, partial two-pack discharge, sequential reveal, post-reveal amendment, lot/supplier/material/lubricant confounding, duty/environment shift, censoring, competing-risk, missingness, dependence, optional-stopping, sparse-tail and estimand-drift twins; validity predicates: GovernanceCommitted sampling frame/estimand/configuration/event/inspection definitions, calibrated metrology, inclusion/exclusion, censoring/missingness, dependence, multiplicity/stopping, AcceptanceCards, checker, candidate-input permissions and access protocol; immutable CandidateFrozen model/toolchain/checker/AcceptanceCard roots; independent custodian joint realization; joined typed population/metrology slot root, two distinct same-ID discharge envelopes and verified receipts, both waiver retirements and authority-head advance in one atomic RealizationCommitted coupled transaction; one-shot joint RevealedForAdjudication and complete access/attempt/adjudication ledger in Closed; laws: no candidate-side principal or transitive capability reads a protected population/metrology root, byte, outcome, aggregate or derived statistic before CandidateFrozen and committed realization, the custodian follows the frozen non-adaptive acquisition/selection mechanism, both packs commit or neither grants authority, frame/event/metrology reconstruction and anytime-valid coverage close, public or laboratory samples never substitute for production population, and out-of-frame extrapolation refuses; shrinkers: reduce lifecycle stages, principals/capabilities, lots, units, events and covariates while retaining the first authority, frame, metrology or coverage fault; replay uses every governed stage receipt and separately addressable joined-root component"
            ),
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                ACCEPTANCE_POLICY_FIXTURE,
                EM_CONVENTION_FIXTURE,
                "i14-motor-bearing-current",
                "i14-bearing-population-adversaries",
                "i14-bearing-population-max-holdout",
                "i14-external-bearing-population-reliability-pack",
                "i14-external-bearing-population-metrology-pack",
            ],
            g3_relations: &[
                "permuting units and deterministic analysis workers preserves estimates and coverage decisions",
                "adding censored or competing-risk outcomes cannot be treated as event-free exposure",
                "narrowing the target population cannot silently reuse a broader estimand without an authenticated amendment",
                "reordering GovernanceCommitted, CandidateFrozen, custodian joint realization, coupled RealizationCommitted, joint RevealedForAdjudication and Closed changes the campaign to IntegrityFailed",
                "changing the frozen frame, metrology, acquisition/selection mechanism, candidate, checker, AcceptanceCard, stopping rule or estimand requires a successor before protected access",
                "slot-only replacement, either same-ID envelope alone, either waiver retirement alone, split authority-head transitions, sequential pack reveal or post-reveal amendment grants no bearing-population authority",
                "removing production-population authority leaves synthetic bearing-path semantics intact but population reliability Unknown",
            ],
            g4_schedule: g4_contract!(
                "request cancellation during GovernanceCommitted freeze, CandidateFrozen sealing, independent custodian joint realization, both same-ID envelope/receipt verifications, atomic coupled RealizationCommitted authority-head commit, joint RevealedForAdjudication access, metrology audit, frame reconstruction, event classification, survival/e-process/coverage computation, adjudication and Closed finalization; drain every candidate/custodian/metrology/estimator/oracle child, forbid post-cancel protected access, either-pack reveal, adjudication and authority-head advancement unless the exact coupled atomic stage had already committed, finalize one population-scoped receipt with no favorable missingness, retain complete frame/event/metrology/censoring/access/attempt/stopping ledgers, joined-root components, both receipts and every failure, then checkpoint, migrate, resume and fork at every stable lifecycle/population frontier; resumed and uninterrupted stage, authority-head, metrology, estimates, coverage and verdict digests must agree"
            ),
            g5_matrix: G5_MATRIX,
            entry_point: "scripts/e2e/leapfrog/i14_bearing_population_validation_max.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i14-bearing-population-validation-max dsr quality --tool frankensim",
            obs_events: &[
                "bearing.governance_committed",
                "bearing.candidate_frozen",
                "bearing.custodian_realized",
                "bearing.coupled_authority_committed",
                "bearing.revealed_for_adjudication",
                "bearing.metrology_audited",
                "bearing.event_audited",
                "bearing.coverage_checked",
                "bearing.reliability_adjudicated",
                "bearing.closed",
                "execution.requested",
                "cancellation.requested",
                "cancellation.observed",
                "execution.cancelled",
                "execution.drained",
                "execution.finalized",
                "checkpoint.saved",
                "checkpoint.resumed",
                "checkpoint.forked",
                "evidence.adjudication_receipt",
                "evidence.failure_bundle",
            ],
            replay_command: "scripts/e2e/leapfrog/i14_bearing_population_validation_max.sh --replay <artifact-id>",
        },
    ]
}

fn i14_waivers() -> Vec<Waiver> {
    vec![
        Waiver {
            subject: "i14-external-standards-edition-clause-pack",
            reason: "licensed AP242, EMC, and rotating-machine standard editions, corrigenda, clause inventories, and restricted crosswalk evidence cannot be embedded in the repository; public schema-shaped fixtures are deliberately not substitutes",
            owner: "I14 standards-crosswalk campaign owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-standards-edition-clause-pack after GovernanceCommitted, authorized least-privilege AuthorizedConstruction and CandidateFrozen, binding exact publisher/edition/corrigenda/license/source, authenticated licensed and complete realized inventory/applicability roots, every derived mapping/exclusion/loss, native dependency, builder/toolchain/checker/AcceptanceCard root, disclosure audit, dual-control signatures and access/output ledger. The same atomic StandardsAuthorityCommitted FrozenManifest::amend transaction installs i14-external-standards-edition-clause-pack as its same-ID typed External discharge-envelope root, removes this Waiver row, verifies the AmendmentRecord and advances the authority head before independent adjudication; no AuthoredSpec slot is replaced. Discharge alone grants no blind, untouched, IID, crosswalk, conformance, legal, certification or regulatory authority; omitted scope, restricted-text leakage, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before any ungoverned licensed-input access and before i14-standards-crosswalk-max independent adjudication or promotion; AuthorizedConstruction access is permitted only through the frozen governed capability, and review is mandatory on every edition, corrigendum, clause-scope, license, crosswalk, candidate, builder/toolchain/checker/AcceptanceCard, principal, sandbox, disclosure-filter, envelope/receipt, authority-head, or signature change",
            promotion_effect: "only i14-governed-standards-crosswalk and its producer remain NoPromotionAuthority while live; native synthetic HarnessGraph/AP242 adapter mechanics and synthetic safety traceability may promote independently, but no public result may be relabeled exact-edition conformance, legal compliance, certification, or regulatory approval",
        },
        Waiver {
            subject: "i14-external-emc-laboratory-calibration-pack",
            reason: "calibrated laboratory procedures, instrument/calibration records, environments, raw traces, uncertainty/covariance, specimen identities, exclusions, retries, and custody cannot be embedded publicly",
            owner: "I14 laboratory-validation campaign owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-emc-laboratory-calibration-pack after GovernanceCommitted, AuthorizedCalibration, AsBuiltModelInstantiation and CandidateFrozen, binding procedure/setup, instruments and calibration certificates/windows, fixture/cable/chamber/environment state, immutable pre-access calibration/model-input content commitments, validation source-universe root, disjoint-membership commitment, exact validation-selection algorithm, pre-candidate secret-seed/VRF or equivalent non-adaptive selection commitment, separately addressable calibration/model-input/untouched-validation roots and membership proofs to those commitments, mutual disjoint-membership proof, non-adaptive selection proof, uncertainty/covariance, retry/exclusion/missingness ledger, complete calibration/model derivation and access ledger, contamination receipt naming every audited candidate-side principal and transitive capability and proving no forbidden validation payload/opening access, candidate/checker/AcceptanceCard roots, custodian and independent adjudicator. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-laboratory-validation-max-holdout with its joined typed External realized root, installs i14-external-emc-laboratory-calibration-pack and i14-external-asbuilt-specimen-geometry-pack as distinct same-ID typed External discharge-envelope roots, removes this Waiver row and the coupled Waiver row through distinct receipts, verifies the AmendmentRecord and advances the authority head before untouched validation reveal; candidate-side validation access or opening before freeze, calibration-to-validation leakage, hidden update, adaptive selection, partial discharge, stratum aliasing, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before any ungoverned calibration/as-built access and before i14-laboratory-validation-max untouched-validation reveal, independent adjudication or promotion; AuthorizedCalibration and AsBuiltModelInstantiation access is permitted only through the frozen governed capability, and review is mandatory on every stratum, pre-access content commitment, validation source-universe/frame commitment, disjoint-membership commitment, validation-selection algorithm, seed/VRF or equivalent non-adaptive commitment, membership/disjointness/non-adaptive proof scheme, audited principal/transitive capability, procedure, setup, instrument, calibration, environment, trace, exclusion, model update, contamination predicate, candidate/model/toolchain/checker/AcceptanceCard identity, discharge envelope or receipt, authority head, custodian signature, independent-adjudicator signature or custody change",
            promotion_effect: "only i14-governed-laboratory-emc-validation and i14-laboratory-validation-max remain NoPromotionAuthority while live; synthetic solver, source/probe/victim, bearing-path, safety-traceability, standards-crosswalk and bearing-population mechanics remain orthogonal",
        },
        Waiver {
            subject: "i14-external-bearing-population-metrology-pack",
            reason: "bearing-specific shaft-voltage/current sensors, discharge-event classifiers, synchronized acquisition, inspection/damage metrology, calibration records, uncertainty/covariance, classifier validation and custody cannot be embedded publicly",
            owner: "I14 bearing-population metrology owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-bearing-population-metrology-pack after GovernanceCommitted and CandidateFrozen, binding bearing-specific sensor and acquisition identities, calibration certificates/windows/ranges, event-classifier version/confusion bounds, synchronized timing, inspection/damage-proxy definitions, raw validation roots, uncertainty/covariance, calibration-validation split, candidate/checker/AcceptanceCard roots, custodian, independent adjudicator and access ledger. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-bearing-population-max-holdout with its joined typed External realized root, installs i14-external-bearing-population-metrology-pack and i14-external-bearing-population-reliability-pack as distinct same-ID typed External discharge-envelope roots, removes this Waiver row and the coupled Waiver row through distinct receipts, verifies the AmendmentRecord and advances the authority head before joint reveal; partial discharge, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before any protected bearing metrology/population access and before i14-bearing-population-validation-max joint RevealedForAdjudication, independent adjudication or promotion; review on every population/frame commitment, acquisition/selection mechanism, sensor, acquisition clock, calibration, event classifier, inspection method, damage proxy, uncertainty, joined root/component, candidate/model/toolchain/checker/AcceptanceCard, envelope/receipt, audited principal/access capability, authority-head or custody change",
            promotion_effect: "only i14-production-bearing-population-reliability and i14-bearing-population-validation-max remain NoPromotionAuthority while live; generic EMC laboratory validation and synthetic bearing-current/event mechanics remain independently scoped",
        },
        Waiver {
            subject: "i14-external-asbuilt-specimen-geometry-pack",
            reason: "proprietary as-built occurrence geometry, tolerances, materials, assembly/aging state, placement, transforms, uncertainty, and specimen provenance cannot be embedded publicly",
            owner: "I14 laboratory-validation campaign owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-asbuilt-specimen-geometry-pack after GovernanceCommitted, AuthorizedCalibration, AsBuiltModelInstantiation and CandidateFrozen, binding exact specimen/configuration and occurrence identities, topology, dimensions/tolerances, materials, assembly/aging/temperature state, source/probe/victim placement, coordinate transforms, uncertainty, provenance, immutable pre-access calibration/model-input content commitments, validation source-universe root, disjoint-membership commitment, exact validation-selection algorithm, pre-candidate secret-seed/VRF or equivalent non-adaptive selection commitment, separately addressable calibration/model-input/untouched-validation roots and membership proofs to those commitments, mutual disjoint-membership proof, non-adaptive selection proof, every authorized model-input derivation/access, contamination receipt naming every audited candidate-side principal and transitive capability and proving no forbidden validation payload/opening access, candidate/checker/AcceptanceCard roots and custodian ledger. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-laboratory-validation-max-holdout with its joined typed External realized root, installs i14-external-asbuilt-specimen-geometry-pack and i14-external-emc-laboratory-calibration-pack as distinct same-ID typed External discharge-envelope roots, removes this Waiver row and the coupled Waiver row through distinct receipts, verifies the AmendmentRecord and advances the authority head before untouched validation reveal; candidate-side validation access or opening before freeze, calibration-to-validation leakage, hidden update, adaptive selection, partial discharge, stratum aliasing, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before any ungoverned as-built/calibration access and before i14-laboratory-validation-max untouched-validation reveal, independent adjudication or promotion; AuthorizedCalibration and AsBuiltModelInstantiation access is permitted only through the frozen governed capability, and review is mandatory on every stratum, pre-access content commitment, validation source-universe/frame commitment, disjoint-membership commitment, validation-selection algorithm, seed/VRF or equivalent non-adaptive commitment, membership/disjointness/non-adaptive proof scheme, audited principal/transitive capability, specimen, configuration, occurrence, geometry, material, placement, state, uncertainty, model derivation, contamination predicate, candidate/model/toolchain/checker/AcceptanceCard identity, discharge envelope or receipt, authority head, custodian signature, independent-adjudicator signature or provenance change",
            promotion_effect: "only i14-governed-laboratory-emc-validation and its producer remain NoPromotionAuthority while live; synthetic geometry, harness, solver and assurance-graph claims may promote without being relabeled as as-built physical validation",
        },
        Waiver {
            subject: "i14-external-bearing-population-reliability-pack",
            reason: "production bearing/material/lubricant populations, supplier/lot histories, duty/environment exposure, censored outcomes, raw traces, inspections, and independent population adjudication cannot be embedded publicly",
            owner: "I14 bearing-population campaign owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-bearing-population-reliability-pack after GovernanceCommitted and CandidateFrozen, binding the sampling frame, lots/suppliers/materials/lubricants, per-unit configuration and exposure lineage, inclusion/exclusion, censoring/missingness/competing-risk policy, outcome definitions, raw observation roots, estimand/stopping/multiplicity rules, candidate/checker/AcceptanceCard roots, custodian, independent adjudicator and access ledger. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-bearing-population-max-holdout with its joined typed External realized root, installs i14-external-bearing-population-reliability-pack and i14-external-bearing-population-metrology-pack as distinct same-ID typed External discharge-envelope roots, removes this Waiver row and the coupled Waiver row through distinct receipts, verifies the AmendmentRecord and advances the authority head before joint reveal; partial discharge, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before any protected bearing population/metrology access and before i14-bearing-population-validation-max joint RevealedForAdjudication, independent adjudication or promotion; review on every population/frame commitment, acquisition/selection mechanism, supplier, lot, material, lubricant, duty, environment, event/outcome definition, observation, censoring/missingness/competing-risk rule, estimand, multiplicity/stopping rule, joined root/component, candidate/model/toolchain/checker/AcceptanceCard, envelope/receipt, audited principal/access capability, authority-head or custody change",
            promotion_effect: "only i14-production-bearing-population-reliability and its producer remain NoPromotionAuthority while live; synthetic bearing-current path/event semantics may promote but cannot be relabeled fluting, erosion, life, warranty, safety, or population reliability",
        },
        Waiver {
            subject: "i14-external-blind-mitigation-custody-pack",
            reason: "the blind robust-mitigation generator commitment, inaccessible realized population, candidate-isolation/access ledger, one-shot reveal, retries/exclusions, multiplicity/stopping evidence, and independent adjudication cannot be embedded before the campaign",
            owner: "I14 robust-mitigation custody owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-blind-mitigation-custody-pack after GovernanceCommitted and CandidateFrozen, binding generator/population commitment, inaccessible realized External/Merkle root, candidate/model/checker/AcceptanceCard digests, candidate isolation/access log, immutable final candidate, one-shot joint reveal, every attempt/exclusion/retry, multiplicity, stopping, missingness, custodian and independent adjudicator. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-mitigation-max-holdout with its typed External realized root, installs i14-external-blind-mitigation-custody-pack as its same-ID typed External discharge-envelope root, removes this Waiver row, verifies the AmendmentRecord and advances the authority head before joint reveal; slot-only replacement, envelope-only replacement, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before i14-robust-mitigation-max accesses any protected root, byte, label, aggregate or derived statistic and before one-shot RevealedForAdjudication, independent adjudication or promotion; review on every generator/population commitment, acquisition/selection mechanism, candidate/model/toolchain/checker/AcceptanceCard, threshold/guard, realized root, envelope/receipt, access capability, exclusion/retry/multiplicity/stopping rule, authority-head or custody change",
            promotion_effect: "only i14-robust-mitigation-heldout and its producer remain NoPromotionAuthority while live; deterministic design mechanics, fidelity theorem and synthetic safety traceability remain independently promotable, but no public holdout may be relabeled blind or untouched",
        },
        Waiver {
            subject: "i14-external-emc-reliability-population-pack",
            reason: "deployment-population source, harness, bond/shield, victim, configuration, duty/environment, outcome, censoring/missingness, calibration and raw observation records cannot be embedded publicly",
            owner: "I14 governed EMC reliability campaign owner (frankensim-leapfrog-2026-program-i94v.2.4.9.1)",
            predicate: "fs-vvreg issues a typed DischargeReceipt for i14-external-emc-reliability-population-pack after GovernanceCommitted and CandidateFrozen, binding the preregistered population frame and estimand, per-unit configuration/exposure lineage, event/outcome definitions, metrology/calibration, inclusion/exclusion, dependence, censoring/missingness, raw realized root, multiplicity/stopping rules, candidate/checker/AcceptanceCard roots, custodian, access ledger and independent adjudicator. The same atomic RealizationCommitted FrozenManifest::amend transaction replaces i14-emc-reliability-max-holdout with its typed External realized root, installs i14-external-emc-reliability-population-pack as its same-ID typed External discharge-envelope root, removes this Waiver row, verifies the AmendmentRecord and advances the authority head before joint reveal; slot-only replacement, envelope-only replacement, split retirement or an unverified transaction is IntegrityFailed",
            expiry: "before i14-emc-reliability-validation-max accesses any protected root, byte, label, aggregate or derived statistic and before one-shot RevealedForAdjudication, independent adjudication or promotion; review on every population/frame commitment, acquisition/selection mechanism, configuration, event/outcome definition, measurement/metrology, censoring/missingness/dependence rule, estimand, multiplicity/stopping rule, candidate/model/toolchain/checker/AcceptanceCard, realized root, envelope/receipt, audited principal/access capability, authority-head or custody change",
            promotion_effect: "only i14-governed-emc-reliability-validation and i14-emc-reliability-validation-max remain NoPromotionAuthority while live; Core UQ inference mechanics and synthetic source/path/victim claims remain independently scoped",
        },
    ]
}
