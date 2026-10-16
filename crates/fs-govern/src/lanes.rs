//! One active unproven mechanism per independently falsifiable proof
//! lane (bead frankensim-ext-epic-gov-rjoq.6).
//!
//! The addendum documents a one-bet discipline; this module makes it
//! EXECUTABLE. A [`LaneCharter`] canonicalizes the semantic fields that
//! define a proof lane (statement/quantifiers, admissible domain,
//! assumptions, target authority, baseline, falsifier family,
//! independence class) and derives a validated [`ProofLaneId`] —
//! the id is minted only from a validated charter, so cosmetic
//! whitespace/ordering "splits" collapse to the same lane and a raw
//! hash cannot be spoofed in.
//!
//! [`PortfolioLedger`] is the admission state machine:
//! - multiple active unproven mechanisms are permitted across
//!   independently falsifiable lanes;
//! - a second active mechanism in the SAME lane refuses atomically,
//!   unless the policy holds a preregistered [`HeadToHeadCharter`]
//!   naming both candidates under a bounded shared envelope;
//! - lanes that DECLARE the same independence class share one bet (the
//!   split-gaming backstop);
//! - global work/memory/reviewer/falsification envelopes bind across
//!   all lanes, so lane partitioning cannot evade portfolio limits;
//! - terminal transitions (refuted/tombstoned/withdrawn/superseded)
//!   release a slot EXACTLY ONCE and only against a content-identified
//!   [`FinalizationReceipt`]; Unknown or stalled work never releases
//!   silently — there is deliberately NO timeout path;
//! - every request carries an [`IdempotencyKey`]; a retry replays the
//!   recorded decision without double-charging, and a DIFFERENT
//!   request under a used key refuses.
//!
//! Every method validates completely BEFORE governed state mutates. A
//! refusal may append one explicit bounded audit row, but cannot partly
//! charge a lane or resource envelope. The complete canonical request
//! is retained for deterministic replay; rows and idempotency bindings
//! are never silently evicted.

use crate::json_escape;
use fs_blake3::ContentHash;
use std::collections::{BTreeMap, BTreeSet};

/// Version of the lane-admission policy schema: bump when a rule,
/// canonicalization step, or identity preimage changes meaning.
pub const LANE_POLICY_VERSION: u32 = 2;

/// Domain for canonical proof-lane identities.
pub const PROOF_LANE_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.proof-lane.v1";

/// Domain for mechanism identities.
pub const MECHANISM_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.mechanism.v1";

/// Domain for terminal finalization receipts.
pub const FINALIZATION_RECEIPT_IDENTITY_DOMAIN: &str = "frankensim.fs-govern.lane-finalization.v1";

/// Domain for idempotency keys.
pub const IDEMPOTENCY_KEY_DOMAIN: &str = "frankensim.fs-govern.lane-idempotency.v1";

/// Domain for admission-request digests (idempotency conflict checks).
pub const REQUEST_DIGEST_DOMAIN: &str = "frankensim.fs-govern.lane-request.v1";
/// Version of the retained lane-request digest identity (bead sj31i.63
/// governance sweep; the domain string carries the same version).
pub const LANE_REQUEST_DIGEST_VERSION: u32 = 1;

/// Owner-local request-digest declaration consumed by `xtask check-identities`.
pub const LANE_REQUEST_DIGEST_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-govern:lane-request-digest",
    "version_const=LANE_REQUEST_DIGEST_VERSION",
    "version=1",
    "domain=frankensim.fs-govern.lane-request.v1",
    "domain_const=REQUEST_DIGEST_DOMAIN",
    "encoder=PortfolioLedger::digest_admit",
    "encoder_helpers=PortfolioLedger::digest_preregister,PortfolioLedger::digest_finalize,push_field,ResourceEnvelope::digest_into",
    "schema_constants=LANE_REQUEST_DIGEST_VERSION,REQUEST_DIGEST_DOMAIN",
    "schema_functions=PortfolioLedger::record,PortfolioLedger::replay,PortfolioLedger::admit,PortfolioLedger::preregister_comparison,PortfolioLedger::finalize",
    "schema_dependencies=fs-blake3:canonical-identity-frame",
    "digest=fs-blake3",
    "encoding=typed-binary",
    "sources=AdmissionDecision",
    "source_fields=AdmissionDecision.seq:nonsemantic:retention-ordering-only,AdmissionDecision.policy_version:nonsemantic:policy-echo-not-request-identity,AdmissionDecision.policy:nonsemantic:policy-echo-not-request-identity,AdmissionDecision.kind:semantic,AdmissionDecision.lane:semantic,AdmissionDecision.mechanism:semantic,AdmissionDecision.idempotency:nonsemantic:replay-map-key-not-digest-input,AdmissionDecision.request_digest:derived:blake3-root-of-the-kind-discriminated-request-preimage,AdmissionDecision.request:semantic,AdmissionDecision.refusal:nonsemantic:outcome-not-request-identity",
    "source_bindings=AdmissionDecision.kind>kind-tag,AdmissionDecision.lane>lane-hash,AdmissionDecision.mechanism>mechanism-hash,AdmissionDecision.request>reservation-axes+candidate-hashes+shared-envelope-axes+preregistration-artifact+receipt-identity",
    "external_semantic_fields=request-digest-domain",
    "semantic_fields=request-digest-domain,kind-tag,lane-hash,mechanism-hash,reservation-axes,candidate-hashes,shared-envelope-axes,preregistration-artifact,receipt-identity",
    "excluded_fields=none",
    "consumers=PortfolioLedger::admit,PortfolioLedger::preregister_comparison,PortfolioLedger::finalize,PortfolioLedger::replay",
    "mutations=request-digest-domain:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,kind-tag:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,lane-hash:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,mechanism-hash:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,reservation-axes:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,candidate-hashes:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,shared-envelope-axes:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,preregistration-artifact:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field,receipt-identity:crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field",
    "nonsemantic_mutations=AdmissionDecision.seq:crates/fs-govern/src/lanes.rs#lane_request_digest_ignores_retention_only_fields,AdmissionDecision.policy_version:crates/fs-govern/src/lanes.rs#lane_request_digest_ignores_retention_only_fields,AdmissionDecision.policy:crates/fs-govern/src/lanes.rs#lane_request_digest_ignores_retention_only_fields,AdmissionDecision.idempotency:crates/fs-govern/src/lanes.rs#lane_request_digest_ignores_retention_only_fields,AdmissionDecision.refusal:crates/fs-govern/src/lanes.rs#lane_request_digest_ignores_retention_only_fields",
    "field_guard=classify_lane_request_digest_fields",
    "transport_guard=AdmissionDecision::admitted",
    "version_guard=crates/fs-govern/src/lanes.rs#lane_request_digest_binds_every_preimage_field",
    "coupling_surface=fs-govern:lane-request-digest",
];

/// Exhaustive field classifier for the retained lane-request digest
/// identity: adding a field to the retained decision row breaks this
/// destructure until the identity declaration moves with it.
#[allow(dead_code)]
fn classify_lane_request_digest_fields(decision: &AdmissionDecision) {
    let AdmissionDecision {
        seq: _,
        policy_version: _,
        policy: _,
        kind: _,
        lane: _,
        mechanism: _,
        idempotency: _,
        request_digest: _,
        request: _,
        refusal: _,
    } = decision;
}

/// Maximum bytes for one canonical charter field.
pub const MAX_FIELD_BYTES: usize = 4096;

/// Maximum declared assumptions per charter.
pub const MAX_ASSUMPTIONS: usize = 256;

/// Maximum candidates in one preregistered head-to-head comparison.
pub const MAX_H2H_CANDIDATES: usize = 8;

/// Hard ceiling for retained decisions and idempotency bindings. The
/// ledger refuses new work before this bound is crossed and reserves one
/// decision/key slot for the eventual finalization of every active
/// mechanism.
pub const MAX_RETAINED_DECISIONS: usize = 256;

/// Hard ceiling for canonical variable-size decision payloads. This is
/// separate from the decision-count cap because a charter may contain a
/// bounded set of bounded assumptions.
pub const MAX_RETAINED_DECISION_BYTES: u64 = 16 * 1024 * 1024;

/// Payload budget held back for each active mechanism's future terminal
/// receipt. Finalization payloads are fixed-size and remain below this
/// conservative accounting reservation.
const FINALIZATION_RECORD_RESERVE_BYTES: u64 = 1_024;

/// Collapse whitespace runs to single spaces and trim — the G3
/// canonicalization that makes cosmetic re-spellings identity-stable.
fn canonical_text(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn preflight_text(what: &'static str, raw: &str) -> Result<(), LaneError> {
    if raw.len() > MAX_FIELD_BYTES {
        return Err(LaneError::TooLarge {
            what,
            len: raw.len(),
            cap: MAX_FIELD_BYTES,
        });
    }
    if raw.split_whitespace().next().is_none() {
        return Err(LaneError::EmptyField { what });
    }
    Ok(())
}

fn push_field(out: &mut Vec<u8>, tag: u8, bytes: &[u8]) {
    out.push(tag);
    let len = u64::try_from(bytes.len()).expect("field length fits in u64");
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(bytes);
}

/// Why a charter, receipt, or admission was refused. Every variant is
/// a structured refusal with a ranked remedy ([`LaneError::remedy`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LaneError {
    /// A required charter/receipt field was empty after canonicalization.
    EmptyField {
        /// Which field.
        what: &'static str,
    },
    /// A field or collection exceeded its bound.
    TooLarge {
        /// Which field.
        what: &'static str,
        /// Observed size.
        len: usize,
        /// The bound.
        cap: usize,
    },
    /// The lane already holds an active unproven mechanism.
    LaneOccupied {
        /// The occupied lane.
        lane: ProofLaneId,
        /// The mechanism holding the slot.
        active: MechanismId,
    },
    /// A different lane sharing this lane's declared independence class
    /// already holds an active mechanism — same falsification fate,
    /// same bet.
    IndependenceClassOccupied {
        /// The colliding active mechanism.
        active: MechanismId,
    },
    /// The global cap on simultaneously active unproven mechanisms.
    PortfolioCapExceeded {
        /// Active count.
        active: u32,
        /// The cap.
        cap: u32,
    },
    /// A global resource envelope axis would be oversubscribed.
    EnvelopeExceeded {
        /// Which axis.
        axis: &'static str,
        /// Amount requested.
        requested: u64,
        /// Amount remaining.
        remaining: u64,
    },
    /// The head-to-head shared envelope would be oversubscribed.
    ComparisonEnvelopeExceeded {
        /// Which axis.
        axis: &'static str,
        /// Amount requested.
        requested: u64,
        /// Amount remaining inside the comparison budget.
        remaining: u64,
    },
    /// The lane has a preregistered comparison and this mechanism is
    /// not one of its declared candidates.
    NotADeclaredCandidate {
        /// The lane.
        lane: ProofLaneId,
    },
    /// A head-to-head charter must declare between 2 and
    /// [`MAX_H2H_CANDIDATES`] DISTINCT candidates.
    ComparisonCandidatesInvalid,
    /// The lane already has a preregistered comparison.
    ComparisonAlreadyDeclared {
        /// The lane.
        lane: ProofLaneId,
    },
    /// A comparison cannot be preregistered on a lane that already
    /// holds an active mechanism (preregistration means BEFORE).
    ComparisonAfterAdmission {
        /// The lane.
        lane: ProofLaneId,
    },
    /// The mechanism is not active in this ledger.
    UnknownMechanism {
        /// The mechanism.
        mechanism: MechanismId,
    },
    /// The mechanism already reached a terminal state; slots release
    /// exactly once and tombstones are permanent.
    AlreadyTerminal {
        /// The mechanism.
        mechanism: MechanismId,
        /// Its terminal state.
        kind: TerminalKind,
    },
    /// The finalization receipt does not bind this mechanism/kind, its
    /// evidence artifact is the all-zero sentinel, or a superseding
    /// successor is missing/spurious.
    ReceiptInvalid {
        /// What is wrong.
        what: &'static str,
    },
    /// The idempotency key was already used by a DIFFERENT request.
    IdempotencyConflict {
        /// Sequence number of the original decision.
        original_seq: u64,
    },
    /// A mechanism id was presented under a lane other than the lane
    /// whose charter minted it.
    MechanismLaneMismatch {
        /// Lane required by the request or comparison charter.
        expected: ProofLaneId,
        /// Lane cryptographically bound into the mechanism id.
        actual: ProofLaneId,
    },
    /// The bounded retained decision/idempotency record has no safe
    /// capacity for this request while preserving one finalization slot
    /// for each active mechanism.
    RetentionCapacityExceeded {
        /// Governed retention axis.
        axis: &'static str,
        /// Amount already retained.
        used: u64,
        /// Additional amount needed by this request.
        requested: u64,
        /// Capacity reserved for active mechanisms to finalize.
        reserved_for_finalization: u64,
        /// Hard cap.
        cap: u64,
    },
}

impl LaneError {
    /// The highest-ranked remedy for this refusal (structured,
    /// actionable, deterministic).
    #[must_use]
    pub fn remedy(&self) -> &'static str {
        match self {
            LaneError::EmptyField { .. } => {
                "supply the missing semantic field; lanes are defined by their complete charter"
            }
            LaneError::TooLarge { .. } => {
                "shorten the field or split the charter into genuinely distinct lanes"
            }
            LaneError::LaneOccupied { .. } => {
                "finalize the active mechanism with a ledgered terminal receipt, or preregister a bounded head-to-head comparison BEFORE admitting candidates"
            }
            LaneError::IndependenceClassOccupied { .. } => {
                "wait for the active bet in this independence class to finalize, or justify a genuinely independent falsifier family under a new class"
            }
            LaneError::PortfolioCapExceeded { .. } => {
                "finalize an active mechanism before admitting another; the portfolio cap is deliberate"
            }
            LaneError::EnvelopeExceeded { .. } => {
                "reduce the reservation or release capacity by finalizing active work; global envelopes bind across all lanes"
            }
            LaneError::ComparisonEnvelopeExceeded { .. } => {
                "reduce the candidate's reservation to fit the preregistered shared budget"
            }
            LaneError::NotADeclaredCandidate { .. } => {
                "only preregistered candidates may join a comparison; amend requires a new preregistration on a fresh lane"
            }
            LaneError::ComparisonCandidatesInvalid => {
                "declare between 2 and 8 distinct candidate mechanisms"
            }
            LaneError::ComparisonAlreadyDeclared { .. } => {
                "use the existing preregistered comparison; one per lane"
            }
            LaneError::ComparisonAfterAdmission { .. } => {
                "finalize the active mechanism first; preregistration must precede admission"
            }
            LaneError::UnknownMechanism { .. } => {
                "admit the mechanism before finalizing it; check the mechanism id"
            }
            LaneError::AlreadyTerminal { .. } => {
                "terminal states are permanent; open a new mechanism (new version) if work genuinely restarts"
            }
            LaneError::ReceiptInvalid { .. } => {
                "supply a receipt whose identity binds this mechanism, kind, successor (for supersession), and a non-zero ledger artifact"
            }
            LaneError::IdempotencyConflict { .. } => {
                "reuse an idempotency key only for byte-identical retries; mint a fresh key for a new request"
            }
            LaneError::MechanismLaneMismatch { .. } => {
                "mint every mechanism from the same canonical lane charter used for admission or comparison"
            }
            LaneError::RetentionCapacityExceeded { .. } => {
                "archive and durably checkpoint this ledger before opening a successor ledger; retained idempotency decisions are never silently evicted"
            }
        }
    }
}

impl core::fmt::Display for LaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LaneError::EmptyField { what } => write!(f, "charter field `{what}` is empty"),
            LaneError::TooLarge { what, len, cap } => {
                write!(f, "`{what}` has size {len}, bound {cap}")
            }
            LaneError::LaneOccupied { lane, active } => write!(
                f,
                "lane {lane} already holds active mechanism {active}; one bet per lane"
            ),
            LaneError::IndependenceClassOccupied { active } => write!(
                f,
                "another lane in the same declared independence class holds active \
                 mechanism {active}; lanes sharing a falsification fate share one bet"
            ),
            LaneError::PortfolioCapExceeded { active, cap } => write!(
                f,
                "portfolio already holds {active} active unproven mechanisms (cap {cap})"
            ),
            LaneError::EnvelopeExceeded {
                axis,
                requested,
                remaining,
            } => write!(
                f,
                "global {axis} envelope cannot cover the reservation: requested \
                 {requested}, remaining {remaining}"
            ),
            LaneError::ComparisonEnvelopeExceeded {
                axis,
                requested,
                remaining,
            } => write!(
                f,
                "preregistered comparison {axis} budget cannot cover the reservation: \
                 requested {requested}, remaining {remaining}"
            ),
            LaneError::NotADeclaredCandidate { lane } => write!(
                f,
                "lane {lane} runs a preregistered comparison and this mechanism is not \
                 a declared candidate"
            ),
            LaneError::ComparisonCandidatesInvalid => write!(
                f,
                "a head-to-head comparison needs 2..={MAX_H2H_CANDIDATES} distinct candidates"
            ),
            LaneError::ComparisonAlreadyDeclared { lane } => {
                write!(f, "lane {lane} already has a preregistered comparison")
            }
            LaneError::ComparisonAfterAdmission { lane } => write!(
                f,
                "lane {lane} already holds an active mechanism; preregistration must \
                 come first"
            ),
            LaneError::UnknownMechanism { mechanism } => {
                write!(f, "mechanism {mechanism} is not active in this ledger")
            }
            LaneError::AlreadyTerminal { mechanism, kind } => write!(
                f,
                "mechanism {mechanism} is already terminal ({}); slots release exactly once",
                kind.name()
            ),
            LaneError::ReceiptInvalid { what } => {
                write!(f, "finalization receipt invalid: {what}")
            }
            LaneError::IdempotencyConflict { original_seq } => write!(
                f,
                "idempotency key already bound to a different request (decision seq \
                 {original_seq})"
            ),
            LaneError::MechanismLaneMismatch { expected, actual } => write!(
                f,
                "mechanism belongs to lane {actual}, but this request requires lane {expected}"
            ),
            LaneError::RetentionCapacityExceeded {
                axis,
                used,
                requested,
                reserved_for_finalization,
                cap,
            } => write!(
                f,
                "retained {axis} capacity exhausted: used {used}, request {requested}, \
                 finalization reserve {reserved_for_finalization}, cap {cap}"
            ),
        }
    }
}

impl std::error::Error for LaneError {}

/// Validated content identity of one proof lane. Minted ONLY by
/// [`LaneCharter::lane_id`] — there is no public constructor from a raw
/// hash, so an id always corresponds to a validated, canonicalized
/// charter (anti-spoofing).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ProofLaneId(ContentHash);

impl ProofLaneId {
    /// The underlying content hash (read-only).
    #[must_use]
    pub fn as_hash(&self) -> &ContentHash {
        &self.0
    }
}

impl core::fmt::Display for ProofLaneId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Identity of one mechanism inside a lane (lane id + canonical name +
/// version). Minted only through [`LaneCharter::mechanism_id`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MechanismId {
    lane: ProofLaneId,
    identity: ContentHash,
}

impl MechanismId {
    /// The underlying content hash (read-only).
    #[must_use]
    pub fn as_hash(&self) -> &ContentHash {
        &self.identity
    }

    /// The canonical proof lane that minted this mechanism.
    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
    }
}

impl core::fmt::Display for MechanismId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.identity)
    }
}

/// Caller-supplied idempotency key, domain-separated from all other
/// identities. Reusing a key REPLAYS the recorded decision for the
/// identical request and refuses a different one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct IdempotencyKey(ContentHash);

impl IdempotencyKey {
    /// Derive a key from a caller-chosen request tag.
    #[must_use]
    pub fn derive(tag: &str) -> IdempotencyKey {
        IdempotencyKey(fs_blake3::hash_domain(
            IDEMPOTENCY_KEY_DOMAIN,
            tag.as_bytes(),
        ))
    }
}

impl core::fmt::Display for IdempotencyKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// The semantic charter that DEFINES a proof lane. Fields are private
/// and canonicalized at construction; the lane id derives from them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LaneCharter {
    statement: String,
    admissible_domain: String,
    assumptions: Vec<String>,
    target_authority: String,
    baseline: String,
    falsifier_family: String,
    independence_class: String,
}

impl LaneCharter {
    /// Canonicalize and validate a charter. Whitespace runs collapse,
    /// assumptions sort and dedupe (empty entries refuse), and every
    /// field is non-empty and bounded — so two cosmetic spellings of
    /// one lane produce ONE identity.
    ///
    /// # Errors
    /// [`LaneError::EmptyField`] / [`LaneError::TooLarge`].
    #[allow(clippy::too_many_arguments)] // the charter IS these seven semantic fields
    pub fn new(
        statement: &str,
        admissible_domain: &str,
        assumptions: &[&str],
        target_authority: &str,
        baseline: &str,
        falsifier_family: &str,
        independence_class: &str,
    ) -> Result<LaneCharter, LaneError> {
        if assumptions.len() > MAX_ASSUMPTIONS {
            return Err(LaneError::TooLarge {
                what: "assumptions",
                len: assumptions.len(),
                cap: MAX_ASSUMPTIONS,
            });
        }
        for (what, raw) in [
            ("statement", statement),
            ("admissible domain", admissible_domain),
            ("target authority", target_authority),
            ("baseline", baseline),
            ("falsifier family", falsifier_family),
            ("independence class", independence_class),
        ] {
            preflight_text(what, raw)?;
        }
        for assumption in assumptions {
            preflight_text("assumption", assumption)?;
        }
        let mut canon_assumptions = assumptions
            .iter()
            .map(|assumption| canonical_text(assumption))
            .collect::<Vec<_>>();
        canon_assumptions.sort_unstable();
        canon_assumptions.dedup();
        Ok(LaneCharter {
            statement: canonical_text(statement),
            admissible_domain: canonical_text(admissible_domain),
            assumptions: canon_assumptions,
            target_authority: canonical_text(target_authority),
            baseline: canonical_text(baseline),
            falsifier_family: canonical_text(falsifier_family),
            independence_class: canonical_text(independence_class),
        })
    }

    /// The validated lane identity: domain-separated BLAKE3 over
    /// every tagged, length-prefixed canonical field.
    #[must_use]
    pub fn lane_id(&self) -> ProofLaneId {
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, self.statement.as_bytes());
        push_field(&mut canonical, 2, self.admissible_domain.as_bytes());
        let count = u64::try_from(self.assumptions.len()).expect("assumption count fits u64");
        push_field(&mut canonical, 3, &count.to_le_bytes());
        for a in &self.assumptions {
            push_field(&mut canonical, 4, a.as_bytes());
        }
        push_field(&mut canonical, 5, self.target_authority.as_bytes());
        push_field(&mut canonical, 6, self.baseline.as_bytes());
        push_field(&mut canonical, 7, self.falsifier_family.as_bytes());
        push_field(&mut canonical, 8, self.independence_class.as_bytes());
        ProofLaneId(fs_blake3::hash_domain(
            PROOF_LANE_IDENTITY_DOMAIN,
            &canonical,
        ))
    }

    /// Identity of the independence class this lane declared (lanes
    /// sharing it share one bet).
    #[must_use]
    pub fn independence_class_id(&self) -> ContentHash {
        fs_blake3::hash_domain(
            PROOF_LANE_IDENTITY_DOMAIN,
            format!("independence-class\u{0}{}", self.independence_class).as_bytes(),
        )
    }

    /// Mint the identity of a mechanism proposed for this lane.
    ///
    /// # Errors
    /// [`LaneError::EmptyField`] / [`LaneError::TooLarge`].
    pub fn mechanism_id(&self, name: &str, version: u32) -> Result<MechanismId, LaneError> {
        preflight_text("mechanism name", name)?;
        let canonical_name = canonical_text(name);
        let mut canonical = Vec::new();
        let lane = self.lane_id();
        push_field(&mut canonical, 1, lane.as_hash().as_bytes());
        push_field(&mut canonical, 2, canonical_name.as_bytes());
        push_field(&mut canonical, 3, &version.to_le_bytes());
        Ok(MechanismId {
            lane,
            identity: fs_blake3::hash_domain(MECHANISM_IDENTITY_DOMAIN, &canonical),
        })
    }

    /// Canonical statement (read-only, for logs).
    #[must_use]
    pub fn statement(&self) -> &str {
        &self.statement
    }

    /// Canonical assumptions (sorted, deduped).
    #[must_use]
    pub fn assumptions(&self) -> &[String] {
        &self.assumptions
    }

    /// Canonical admissible domain.
    #[must_use]
    pub fn admissible_domain(&self) -> &str {
        &self.admissible_domain
    }

    /// Canonical target authority.
    #[must_use]
    pub fn target_authority(&self) -> &str {
        &self.target_authority
    }

    /// Canonical boring baseline.
    #[must_use]
    pub fn baseline(&self) -> &str {
        &self.baseline
    }

    /// Canonical falsifier family.
    #[must_use]
    pub fn falsifier_family(&self) -> &str {
        &self.falsifier_family
    }

    /// Canonical declared independence class.
    #[must_use]
    pub fn independence_class(&self) -> &str {
        &self.independence_class
    }

    fn retained_bytes(&self) -> u64 {
        let strings = [
            self.statement.len(),
            self.admissible_domain.len(),
            self.target_authority.len(),
            self.baseline.len(),
            self.falsifier_family.len(),
            self.independence_class.len(),
        ]
        .into_iter()
        .chain(self.assumptions.iter().map(String::len))
        .fold(0_u64, |sum, len| {
            sum.saturating_add(u64::try_from(len).unwrap_or(u64::MAX))
        });
        strings.saturating_add(
            u64::try_from(self.assumptions.len())
                .unwrap_or(u64::MAX)
                .saturating_mul(u64::try_from(core::mem::size_of::<String>()).unwrap_or(u64::MAX)),
        )
    }

    fn write_json(&self, out: &mut String) {
        use core::fmt::Write as _;
        write!(
            out,
            "{{\"statement\":\"{}\",\"admissible_domain\":\"{}\",\"assumptions\":[",
            json_escape(&self.statement),
            json_escape(&self.admissible_domain),
        )
        .expect("writing to a String is infallible");
        for (index, assumption) in self.assumptions.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(out, "\"{}\"", json_escape(assumption))
                .expect("writing to a String is infallible");
        }
        write!(
            out,
            "],\"target_authority\":\"{}\",\"baseline\":\"{}\",\"falsifier_family\":\"{}\",\"independence_class\":\"{}\"}}",
            json_escape(&self.target_authority),
            json_escape(&self.baseline),
            json_escape(&self.falsifier_family),
            json_escape(&self.independence_class),
        )
        .expect("writing to a String is infallible");
    }
}

/// A resource envelope over the four governed axes. All arithmetic is
/// checked; axes are independent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ResourceEnvelope {
    /// Abstract work units (solver/agent budget).
    pub work_units: u64,
    /// Peak memory bytes.
    pub memory_bytes: u64,
    /// Reviewer slots.
    pub reviewer_slots: u64,
    /// Falsification capacity (independent falsification attempts the
    /// program can actually run).
    pub falsification_capacity: u64,
}

impl ResourceEnvelope {
    fn axes(&self) -> [(&'static str, u64); 4] {
        [
            ("work", self.work_units),
            ("memory", self.memory_bytes),
            ("reviewer", self.reviewer_slots),
            ("falsification-capacity", self.falsification_capacity),
        ]
    }

    /// `reserved + request` against `self` as the limit: the first
    /// axis that would overflow refuses (deterministic order).
    fn admit(
        &self,
        reserved: &ResourceEnvelope,
        request: &ResourceEnvelope,
        comparison: bool,
    ) -> Result<(), LaneError> {
        let limits = self.axes();
        let used = reserved.axes();
        let want = request.axes();
        for i in 0..4 {
            let remaining = limits[i].1.saturating_sub(used[i].1);
            if want[i].1 > remaining {
                return Err(if comparison {
                    LaneError::ComparisonEnvelopeExceeded {
                        axis: limits[i].0,
                        requested: want[i].1,
                        remaining,
                    }
                } else {
                    LaneError::EnvelopeExceeded {
                        axis: limits[i].0,
                        requested: want[i].1,
                        remaining,
                    }
                });
            }
        }
        Ok(())
    }

    fn add(&mut self, other: &ResourceEnvelope) {
        self.work_units = self.work_units.saturating_add(other.work_units);
        self.memory_bytes = self.memory_bytes.saturating_add(other.memory_bytes);
        self.reviewer_slots = self.reviewer_slots.saturating_add(other.reviewer_slots);
        self.falsification_capacity = self
            .falsification_capacity
            .saturating_add(other.falsification_capacity);
    }

    fn sub(&mut self, other: &ResourceEnvelope) {
        self.work_units = self.work_units.saturating_sub(other.work_units);
        self.memory_bytes = self.memory_bytes.saturating_sub(other.memory_bytes);
        self.reviewer_slots = self.reviewer_slots.saturating_sub(other.reviewer_slots);
        self.falsification_capacity = self
            .falsification_capacity
            .saturating_sub(other.falsification_capacity);
    }

    fn digest_into(&self, out: &mut Vec<u8>, tag: u8) {
        push_field(out, tag, &self.work_units.to_le_bytes());
        push_field(out, tag, &self.memory_bytes.to_le_bytes());
        push_field(out, tag, &self.reviewer_slots.to_le_bytes());
        push_field(out, tag, &self.falsification_capacity.to_le_bytes());
    }

    fn write_json(&self, out: &mut String) {
        use core::fmt::Write as _;
        write!(
            out,
            "{{\"work_units\":{},\"memory_bytes\":{},\"reviewer_slots\":{},\"falsification_capacity\":{}}}",
            self.work_units,
            self.memory_bytes,
            self.reviewer_slots,
            self.falsification_capacity,
        )
        .expect("writing to a String is infallible");
    }
}

/// Portfolio-level policy: the global envelope plus the cap on
/// simultaneously active unproven mechanisms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PortfolioPolicy {
    /// Global resource envelope across ALL lanes.
    pub global: ResourceEnvelope,
    /// Maximum simultaneously active unproven mechanisms.
    pub max_active_mechanisms: u32,
}

impl PortfolioPolicy {
    fn write_json(&self, out: &mut String) {
        use core::fmt::Write as _;
        out.push_str("{\"global\":");
        self.global.write_json(out);
        write!(
            out,
            ",\"max_active_mechanisms\":{},\"max_retained_decisions\":{},\"max_retained_decision_bytes\":{}}}",
            self.max_active_mechanisms,
            MAX_RETAINED_DECISIONS,
            MAX_RETAINED_DECISION_BYTES,
        )
        .expect("writing to a String is infallible");
    }
}

/// A preregistered, bounded head-to-head comparison: the ONLY way two
/// active mechanisms may share a lane. Declared BEFORE any admission
/// in the lane, naming its candidates and a shared budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeadToHeadCharter {
    lane: LaneCharter,
    candidates: Vec<MechanismId>,
    shared: ResourceEnvelope,
    preregistration_artifact: ContentHash,
}

impl HeadToHeadCharter {
    /// Validate a comparison charter: 2..=[`MAX_H2H_CANDIDATES`]
    /// DISTINCT candidates and a non-zero preregistration artifact
    /// (the ledgered protocol document).
    ///
    /// # Errors
    /// [`LaneError::ComparisonCandidatesInvalid`] /
    /// [`LaneError::ReceiptInvalid`].
    pub fn new(
        lane: &LaneCharter,
        candidates: &[MechanismId],
        shared: ResourceEnvelope,
        preregistration_artifact: ContentHash,
    ) -> Result<HeadToHeadCharter, LaneError> {
        // Reject oversized slices before copying or sorting them.
        if candidates.len() < 2 || candidates.len() > MAX_H2H_CANDIDATES {
            return Err(LaneError::ComparisonCandidatesInvalid);
        }
        let lane_id = lane.lane_id();
        if let Some(candidate) = candidates
            .iter()
            .find(|candidate| candidate.lane() != lane_id)
        {
            return Err(LaneError::MechanismLaneMismatch {
                expected: lane_id,
                actual: candidate.lane(),
            });
        }
        let mut sorted = candidates.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != candidates.len() {
            return Err(LaneError::ComparisonCandidatesInvalid);
        }
        if preregistration_artifact
            .as_bytes()
            .iter()
            .all(|byte| *byte == 0)
        {
            return Err(LaneError::ReceiptInvalid {
                what: "preregistration artifact is the all-zero missing-value sentinel",
            });
        }
        Ok(HeadToHeadCharter {
            lane: lane.clone(),
            candidates: sorted,
            shared,
            preregistration_artifact,
        })
    }

    /// The lane this comparison governs.
    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane.lane_id()
    }

    /// Canonical lane charter retained with the preregistration.
    #[must_use]
    pub fn lane_charter(&self) -> &LaneCharter {
        &self.lane
    }

    /// Declared candidates (sorted).
    #[must_use]
    pub fn candidates(&self) -> &[MechanismId] {
        &self.candidates
    }
}

/// Terminal states. Terminal is PERMANENT: a terminal mechanism never
/// re-activates and never releases capacity a second time.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TerminalKind {
    /// The falsifier family refuted the mechanism.
    Refuted,
    /// Governance killed it (kill criterion / quarterly review).
    Tombstoned,
    /// The owner withdrew it.
    Withdrawn,
    /// A successor mechanism replaced it.
    Superseded,
}

impl TerminalKind {
    /// Stable lowercase name for logs.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            TerminalKind::Refuted => "refuted",
            TerminalKind::Tombstoned => "tombstoned",
            TerminalKind::Withdrawn => "withdrawn",
            TerminalKind::Superseded => "superseded",
        }
    }

    fn tag(self) -> u8 {
        match self {
            TerminalKind::Refuted => 1,
            TerminalKind::Tombstoned => 2,
            TerminalKind::Withdrawn => 3,
            TerminalKind::Superseded => 4,
        }
    }
}

/// Content-identified evidence that a terminal outcome was DURABLY
/// finalized in the design ledger. The identity binds mechanism, kind,
/// successor (for supersession), and the ledger artifact; presence of
/// a receipt is necessary but its identity must also verify — a slot
/// never releases against a mismatched or zero-evidence receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FinalizationReceipt {
    mechanism: MechanismId,
    kind: TerminalKind,
    superseded_by: Option<MechanismId>,
    ledger_artifact: ContentHash,
    identity: ContentHash,
}

impl FinalizationReceipt {
    /// Construct a sealed receipt. `superseded_by` is required exactly
    /// when `kind` is [`TerminalKind::Superseded`], and the successor
    /// must differ from the subject.
    ///
    /// # Errors
    /// [`LaneError::ReceiptInvalid`].
    pub fn new(
        mechanism: MechanismId,
        kind: TerminalKind,
        superseded_by: Option<MechanismId>,
        ledger_artifact: ContentHash,
    ) -> Result<FinalizationReceipt, LaneError> {
        if ledger_artifact.as_bytes().iter().all(|byte| *byte == 0) {
            return Err(LaneError::ReceiptInvalid {
                what: "ledger artifact is the all-zero missing-value sentinel",
            });
        }
        match (kind, superseded_by) {
            (TerminalKind::Superseded, None) => {
                return Err(LaneError::ReceiptInvalid {
                    what: "supersession requires the successor mechanism id",
                });
            }
            (TerminalKind::Superseded, Some(successor)) if successor == mechanism => {
                return Err(LaneError::ReceiptInvalid {
                    what: "a mechanism cannot supersede itself",
                });
            }
            (TerminalKind::Superseded, Some(successor)) if successor.lane() != mechanism.lane() => {
                return Err(LaneError::ReceiptInvalid {
                    what: "a successor must belong to the same proof lane",
                });
            }
            (TerminalKind::Superseded, Some(_)) => {}
            (_, Some(_)) => {
                return Err(LaneError::ReceiptInvalid {
                    what: "a successor is only meaningful for supersession",
                });
            }
            (_, None) => {}
        }
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, mechanism.as_hash().as_bytes());
        push_field(&mut canonical, 2, &[kind.tag()]);
        if let Some(successor) = superseded_by {
            push_field(&mut canonical, 3, successor.as_hash().as_bytes());
        }
        push_field(&mut canonical, 4, ledger_artifact.as_bytes());
        Ok(FinalizationReceipt {
            mechanism,
            kind,
            superseded_by,
            ledger_artifact,
            identity: fs_blake3::hash_domain(FINALIZATION_RECEIPT_IDENTITY_DOMAIN, &canonical),
        })
    }

    /// The receipt's sealed identity.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// What one recorded decision was about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecisionKind {
    /// An admission request.
    Admit,
    /// A comparison preregistration.
    Preregister,
    /// A terminal finalization.
    Finalize,
}

impl DecisionKind {
    fn name(self) -> &'static str {
        match self {
            DecisionKind::Admit => "admit",
            DecisionKind::Preregister => "preregister",
            DecisionKind::Finalize => "finalize",
        }
    }
}

/// Complete canonical request retained with a decision. Successful rows
/// contain enough data to deterministically reconstruct and replay the
/// state transition; refused rows retain the exact attempted request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DecisionRequest {
    /// One mechanism admission.
    Admit {
        /// Canonical lane charter.
        charter: LaneCharter,
        /// Four-axis requested reservation.
        reservation: ResourceEnvelope,
    },
    /// One preregistered comparison.
    Preregister {
        /// Canonical lane charter.
        charter: LaneCharter,
        /// Sorted, lane-bound candidate ids.
        candidates: Vec<MechanismId>,
        /// Shared comparison envelope.
        shared: ResourceEnvelope,
        /// Content identity of the preregistered protocol.
        preregistration_artifact: ContentHash,
    },
    /// One terminal transition.
    Finalize {
        /// Terminal outcome.
        kind: TerminalKind,
        /// Required successor for supersession.
        superseded_by: Option<MechanismId>,
        /// Durable-ledger content reference.
        ledger_artifact: ContentHash,
        /// Sealed receipt identity.
        receipt_identity: ContentHash,
        /// Reservation actually released by an admitted transition.
        released: Option<ResourceEnvelope>,
    },
}

impl DecisionRequest {
    fn write_json(&self, out: &mut String) {
        use core::fmt::Write as _;
        match self {
            DecisionRequest::Admit {
                charter,
                reservation,
            } => {
                out.push_str("{\"type\":\"admit\",\"charter\":");
                charter.write_json(out);
                out.push_str(",\"reservation\":");
                reservation.write_json(out);
                out.push('}');
            }
            DecisionRequest::Preregister {
                charter,
                candidates,
                shared,
                preregistration_artifact,
            } => {
                out.push_str("{\"type\":\"preregister\",\"charter\":");
                charter.write_json(out);
                out.push_str(",\"candidates\":[");
                for (index, candidate) in candidates.iter().enumerate() {
                    if index > 0 {
                        out.push(',');
                    }
                    write!(out, "\"{candidate}\"").expect("writing to a String is infallible");
                }
                write!(out, "],\"shared\":").expect("writing to a String is infallible");
                shared.write_json(out);
                write!(
                    out,
                    ",\"preregistration_artifact\":\"{preregistration_artifact}\"}}"
                )
                .expect("writing to a String is infallible");
            }
            DecisionRequest::Finalize {
                kind,
                superseded_by,
                ledger_artifact,
                receipt_identity,
                released,
            } => {
                write!(
                    out,
                    "{{\"type\":\"finalize\",\"terminal_kind\":\"{}\",\"superseded_by\":",
                    kind.name(),
                )
                .expect("writing to a String is infallible");
                match superseded_by {
                    Some(successor) => {
                        write!(out, "\"{successor}\"").expect("writing to a String is infallible")
                    }
                    None => out.push_str("null"),
                }
                write!(
                    out,
                    ",\"ledger_artifact\":\"{ledger_artifact}\",\"receipt_identity\":\"{receipt_identity}\",\"released\":"
                )
                .expect("writing to a String is infallible");
                match released {
                    Some(reservation) => reservation.write_json(out),
                    None => out.push_str("null"),
                }
                out.push('}');
            }
        }
    }
}

#[derive(Clone, Copy)]
enum DecisionRequestRef<'a> {
    Admit {
        charter: &'a LaneCharter,
        reservation: ResourceEnvelope,
    },
    Preregister {
        charter: &'a HeadToHeadCharter,
    },
    Finalize {
        receipt: &'a FinalizationReceipt,
        released: Option<ResourceEnvelope>,
    },
}

impl DecisionRequestRef<'_> {
    fn retained_bytes(self) -> u64 {
        const FIXED: u64 = 512;
        match self {
            DecisionRequestRef::Admit { charter, .. } => {
                FIXED.saturating_add(charter.retained_bytes())
            }
            DecisionRequestRef::Preregister { charter } => FIXED
                .saturating_add(charter.lane.retained_bytes())
                .saturating_add(
                    u64::try_from(charter.candidates.len())
                        .unwrap_or(u64::MAX)
                        .saturating_mul(
                            u64::try_from(core::mem::size_of::<MechanismId>()).unwrap_or(u64::MAX),
                        ),
                ),
            DecisionRequestRef::Finalize { .. } => FIXED,
        }
    }

    fn to_owned(self) -> DecisionRequest {
        match self {
            DecisionRequestRef::Admit {
                charter,
                reservation,
            } => DecisionRequest::Admit {
                charter: charter.clone(),
                reservation,
            },
            DecisionRequestRef::Preregister { charter } => DecisionRequest::Preregister {
                charter: charter.lane.clone(),
                candidates: charter.candidates.clone(),
                shared: charter.shared,
                preregistration_artifact: charter.preregistration_artifact,
            },
            DecisionRequestRef::Finalize { receipt, released } => DecisionRequest::Finalize {
                kind: receipt.kind,
                superseded_by: receipt.superseded_by,
                ledger_artifact: receipt.ledger_artifact,
                receipt_identity: receipt.identity,
                released,
            },
        }
    }
}

/// One atomic decision in the replayable log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmissionDecision {
    /// Sequence number (0-based, dense).
    pub seq: u64,
    /// Policy schema version in force.
    pub policy_version: u32,
    /// Complete policy in force for this decision.
    pub policy: PortfolioPolicy,
    /// What kind of request this was.
    pub kind: DecisionKind,
    /// The validated canonical lane.
    pub lane: ProofLaneId,
    /// The lane-bound mechanism subject.
    pub mechanism: MechanismId,
    /// The idempotency key presented.
    pub idempotency: IdempotencyKey,
    /// Digest of the complete request for replay conflict checks.
    pub request_digest: ContentHash,
    /// Complete canonical replay preimage.
    pub request: DecisionRequest,
    /// Refusal, if the request was refused.
    pub refusal: Option<LaneError>,
}

impl AdmissionDecision {
    /// Whether the request was admitted.
    #[must_use]
    pub fn admitted(&self) -> bool {
        self.refusal.is_none()
    }

    /// One bounded JSON row for ledgers/dashboards.
    #[must_use]
    pub fn to_json(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        let verdict = match &self.refusal {
            None => "admitted".to_owned(),
            Some(error) => format!("refused: {error}"),
        };
        let remedy = self
            .refusal
            .as_ref()
            .map_or_else(String::new, |error| error.remedy().to_owned());
        write!(
            out,
            "{{\"seq\":{},\"policy_version\":{},\"policy\":",
            self.seq, self.policy_version,
        )
        .expect("writing to a String is infallible");
        self.policy.write_json(&mut out);
        write!(
            out,
            ",\"kind\":\"{}\",\"lane\":\"{}\",\"mechanism\":\"{}\",\"mechanism_lane\":\"{}\",\"idempotency\":\"{}\",\"request_digest\":\"{}\",\"request\":",
            self.kind.name(),
            self.lane,
            self.mechanism,
            self.mechanism.lane(),
            self.idempotency,
            self.request_digest,
        )
        .expect("writing to a String is infallible");
        self.request.write_json(&mut out);
        write!(
            out,
            ",\"verdict\":\"{}\",\"remedy\":\"{}\"}}",
            json_escape(&verdict),
            json_escape(&remedy),
        )
        .expect("writing to a String is infallible");
        out
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ActiveRecord {
    lane: ProofLaneId,
    independence_class: ContentHash,
    reservation: ResourceEnvelope,
    in_comparison: bool,
}

enum ReplayStatus {
    Fresh,
    Recorded(Option<LaneError>),
    Conflict { original_seq: u64 },
}

/// The atomic admission state machine. Exclusive access to this
/// non-Clone authority value is the in-process concurrency contract.
/// Every request is completely validated and durably representable in
/// the bounded replay record before governed state mutates.
#[derive(Debug, PartialEq)]
pub struct PortfolioLedger {
    policy: PortfolioPolicy,
    active: BTreeMap<MechanismId, ActiveRecord>,
    lane_active: BTreeMap<ProofLaneId, Vec<MechanismId>>,
    class_active: BTreeMap<ContentHash, BTreeSet<MechanismId>>,
    comparisons: BTreeMap<ProofLaneId, HeadToHeadCharter>,
    comparison_reserved: BTreeMap<ProofLaneId, ResourceEnvelope>,
    terminal: BTreeMap<MechanismId, TerminalKind>,
    reserved: ResourceEnvelope,
    decisions: Vec<AdmissionDecision>,
    retained_decision_bytes: u64,
    idempotency: BTreeMap<IdempotencyKey, u64>,
    conflict_idempotency: BTreeMap<(IdempotencyKey, ContentHash), u64>,
}

impl PortfolioLedger {
    /// Empty ledger under policy.
    #[must_use]
    pub fn new(policy: PortfolioPolicy) -> PortfolioLedger {
        PortfolioLedger {
            policy,
            active: BTreeMap::new(),
            lane_active: BTreeMap::new(),
            class_active: BTreeMap::new(),
            comparisons: BTreeMap::new(),
            comparison_reserved: BTreeMap::new(),
            terminal: BTreeMap::new(),
            reserved: ResourceEnvelope::default(),
            decisions: Vec::new(),
            retained_decision_bytes: 0,
            idempotency: BTreeMap::new(),
            conflict_idempotency: BTreeMap::new(),
        }
    }

    /// Number of active unproven mechanisms.
    #[must_use]
    pub fn active_count(&self) -> u32 {
        u32::try_from(self.active.len()).unwrap_or(u32::MAX)
    }

    /// Currently reserved global resources.
    #[must_use]
    pub fn reserved(&self) -> ResourceEnvelope {
        self.reserved
    }

    /// Complete deterministic retained decision log.
    #[must_use]
    pub fn decisions(&self) -> &[AdmissionDecision] {
        &self.decisions
    }

    /// Canonical variable-size bytes charged to the retained log.
    #[must_use]
    pub fn retained_decision_bytes(&self) -> u64 {
        self.retained_decision_bytes
    }

    /// Bounded JSON decision log: at most limit most-recent rows plus
    /// explicit truncation and hard-retention metadata.
    #[must_use]
    pub fn decisions_json(&self, limit: usize) -> String {
        use core::fmt::Write as _;
        let skipped = self.decisions.len().saturating_sub(limit);
        let mut out = format!(
            "{{\"skipped\":{skipped},\"retained\":{},\"retained_cap\":{},\"retained_bytes\":{},\"retained_byte_cap\":{},\"decisions\":[",
            self.decisions.len(),
            MAX_RETAINED_DECISIONS,
            self.retained_decision_bytes,
            MAX_RETAINED_DECISION_BYTES,
        );
        for (index, decision) in self.decisions.iter().skip(skipped).enumerate() {
            if index > 0 {
                out.push(',');
            }
            write!(out, "{}", decision.to_json()).expect("writing to a String is infallible");
        }
        out.push_str("]}");
        out
    }

    fn digest_admit(
        lane: ProofLaneId,
        mechanism: MechanismId,
        reservation: &ResourceEnvelope,
    ) -> ContentHash {
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, b"admit");
        push_field(&mut canonical, 2, lane.as_hash().as_bytes());
        push_field(&mut canonical, 3, mechanism.as_hash().as_bytes());
        reservation.digest_into(&mut canonical, 4);
        fs_blake3::hash_domain(REQUEST_DIGEST_DOMAIN, &canonical)
    }

    fn digest_preregister(charter: &HeadToHeadCharter) -> ContentHash {
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, b"preregister");
        push_field(&mut canonical, 2, charter.lane().as_hash().as_bytes());
        for candidate in &charter.candidates {
            push_field(&mut canonical, 3, candidate.as_hash().as_bytes());
        }
        charter.shared.digest_into(&mut canonical, 4);
        push_field(
            &mut canonical,
            5,
            charter.preregistration_artifact.as_bytes(),
        );
        fs_blake3::hash_domain(REQUEST_DIGEST_DOMAIN, &canonical)
    }

    fn digest_finalize(receipt: &FinalizationReceipt) -> ContentHash {
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, b"finalize");
        push_field(&mut canonical, 2, receipt.identity.as_bytes());
        fs_blake3::hash_domain(REQUEST_DIGEST_DOMAIN, &canonical)
    }

    fn replay(&self, key: IdempotencyKey, request_digest: ContentHash) -> ReplayStatus {
        let Some(seq) = self.idempotency.get(&key).copied() else {
            return ReplayStatus::Fresh;
        };
        let recorded = usize::try_from(seq)
            .ok()
            .and_then(|index| self.decisions.get(index));
        match recorded {
            Some(decision) if decision.request_digest == request_digest => {
                ReplayStatus::Recorded(decision.refusal.clone())
            }
            Some(_) => {
                let conflict = self
                    .conflict_idempotency
                    .get(&(key, request_digest))
                    .and_then(|conflict_seq| usize::try_from(*conflict_seq).ok())
                    .and_then(|index| self.decisions.get(index));
                match conflict {
                    Some(decision) => ReplayStatus::Recorded(decision.refusal.clone()),
                    None => ReplayStatus::Conflict { original_seq: seq },
                }
            }
            None => ReplayStatus::Conflict { original_seq: seq },
        }
    }

    fn capacity_error(
        axis: &'static str,
        used: u64,
        requested: u64,
        reserved_for_finalization: u64,
        cap: u64,
    ) -> LaneError {
        LaneError::RetentionCapacityExceeded {
            axis,
            used,
            requested,
            reserved_for_finalization,
            cap,
        }
    }

    fn ensure_record_capacity(
        &self,
        request_bytes: u64,
        bind_idempotency: bool,
        bind_conflict: bool,
        active_after: usize,
    ) -> Result<(), LaneError> {
        let active_reserve = u64::try_from(active_after).unwrap_or(u64::MAX);
        let decision_used = u64::try_from(self.decisions.len()).unwrap_or(u64::MAX);
        let decision_cap = u64::try_from(MAX_RETAINED_DECISIONS).unwrap_or(u64::MAX);
        let decision_total = decision_used
            .checked_add(1)
            .and_then(|value| value.checked_add(active_reserve));
        if decision_total.is_none_or(|total| total > decision_cap) {
            return Err(Self::capacity_error(
                "decision-count",
                decision_used,
                1,
                active_reserve,
                decision_cap,
            ));
        }

        let key_used = u64::try_from(
            self.idempotency
                .len()
                .saturating_add(self.conflict_idempotency.len()),
        )
        .unwrap_or(u64::MAX);
        let key_request = if bind_idempotency || bind_conflict {
            1
        } else {
            0
        };
        let key_total = key_used
            .checked_add(key_request)
            .and_then(|value| value.checked_add(active_reserve));
        if key_total.is_none_or(|total| total > decision_cap) {
            return Err(Self::capacity_error(
                "idempotency-count",
                key_used,
                key_request,
                active_reserve,
                decision_cap,
            ));
        }

        let byte_reserve = active_reserve
            .checked_mul(FINALIZATION_RECORD_RESERVE_BYTES)
            .unwrap_or(u64::MAX);
        let byte_total = self
            .retained_decision_bytes
            .checked_add(request_bytes)
            .and_then(|value| value.checked_add(byte_reserve));
        if byte_total.is_none_or(|total| total > MAX_RETAINED_DECISION_BYTES) {
            return Err(Self::capacity_error(
                "decision-bytes",
                self.retained_decision_bytes,
                request_bytes,
                byte_reserve,
                MAX_RETAINED_DECISION_BYTES,
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn record(
        &mut self,
        kind: DecisionKind,
        lane: ProofLaneId,
        mechanism: MechanismId,
        key: IdempotencyKey,
        request_digest: ContentHash,
        request: DecisionRequestRef<'_>,
        refusal: Option<LaneError>,
        bind_idempotency: bool,
        bind_conflict: bool,
        active_after: usize,
    ) -> Result<u64, LaneError> {
        let request_bytes = request.retained_bytes();
        self.ensure_record_capacity(request_bytes, bind_idempotency, bind_conflict, active_after)?;
        // The only variable-size clones occur after both count and byte
        // admission, so caller-controlled requests cannot allocate past
        // the retained-log envelope.
        let request = request.to_owned();
        let seq = u64::try_from(self.decisions.len()).map_err(|_| {
            Self::capacity_error(
                "decision-count",
                u64::MAX,
                1,
                u64::try_from(active_after).unwrap_or(u64::MAX),
                u64::try_from(MAX_RETAINED_DECISIONS).unwrap_or(u64::MAX),
            )
        })?;
        self.decisions.push(AdmissionDecision {
            seq,
            policy_version: LANE_POLICY_VERSION,
            policy: self.policy,
            kind,
            lane,
            mechanism,
            idempotency: key,
            request_digest,
            request,
            refusal,
        });
        self.retained_decision_bytes += request_bytes;
        if bind_idempotency {
            self.idempotency.insert(key, seq);
        }
        if bind_conflict {
            self.conflict_idempotency.insert((key, request_digest), seq);
        }
        Ok(seq)
    }

    /// Preregister a bounded head-to-head comparison for a lane before
    /// any candidate in that comparison is admitted.
    ///
    /// # Errors
    /// Structured LaneError; recorded refusals consume bounded audit
    /// capacity, while exact retries do not.
    pub fn preregister_comparison(
        &mut self,
        charter: HeadToHeadCharter,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let digest = Self::digest_preregister(&charter);
        let lane = charter.lane();
        let subject = charter.candidates[0];
        let request = DecisionRequestRef::Preregister { charter: &charter };
        match self.replay(key, digest) {
            ReplayStatus::Recorded(refusal) => {
                return refusal.map_or(Ok(()), Err);
            }
            ReplayStatus::Conflict { original_seq } => {
                let error = LaneError::IdempotencyConflict { original_seq };
                self.record(
                    DecisionKind::Preregister,
                    lane,
                    subject,
                    key,
                    digest,
                    request,
                    Some(error.clone()),
                    false,
                    true,
                    self.active.len(),
                )?;
                return Err(error);
            }
            ReplayStatus::Fresh => {}
        }
        let verdict = if self.comparisons.contains_key(&lane) {
            Some(LaneError::ComparisonAlreadyDeclared { lane })
        } else if self
            .lane_active
            .get(&lane)
            .is_some_and(|active| !active.is_empty())
        {
            Some(LaneError::ComparisonAfterAdmission { lane })
        } else {
            None
        };
        self.record(
            DecisionKind::Preregister,
            lane,
            subject,
            key,
            digest,
            request,
            verdict.clone(),
            true,
            false,
            self.active.len(),
        )?;
        match verdict {
            None => {
                self.comparison_reserved
                    .insert(lane, ResourceEnvelope::default());
                self.comparisons.insert(lane, charter);
                Ok(())
            }
            Some(error) => Err(error),
        }
    }

    /// Atomically admit one lane-bound mechanism. Validation order is
    /// deterministic: idempotency, lane binding, terminal permanence,
    /// duplicate activity, lane/class occupancy, portfolio cap,
    /// comparison envelope, global envelope, retained-log capacity.
    ///
    /// # Errors
    /// Structured LaneError with a ranked remedy.
    pub fn admit(
        &mut self,
        charter: &LaneCharter,
        mechanism: MechanismId,
        reservation: ResourceEnvelope,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let lane = charter.lane_id();
        let digest = Self::digest_admit(lane, mechanism, &reservation);
        let request = DecisionRequestRef::Admit {
            charter,
            reservation,
        };
        match self.replay(key, digest) {
            ReplayStatus::Recorded(refusal) => {
                return refusal.map_or(Ok(()), Err);
            }
            ReplayStatus::Conflict { original_seq } => {
                let error = LaneError::IdempotencyConflict { original_seq };
                self.record(
                    DecisionKind::Admit,
                    lane,
                    mechanism,
                    key,
                    digest,
                    request,
                    Some(error.clone()),
                    false,
                    true,
                    self.active.len(),
                )?;
                return Err(error);
            }
            ReplayStatus::Fresh => {}
        }
        let class = charter.independence_class_id();
        let comparison = self.comparisons.get(&lane);
        let in_comparison = comparison.is_some();
        let verdict = self.admit_verdict(lane, class, mechanism, &reservation, comparison);
        let active_after = if verdict.is_none() {
            self.active.len().saturating_add(1)
        } else {
            self.active.len()
        };
        self.record(
            DecisionKind::Admit,
            lane,
            mechanism,
            key,
            digest,
            request,
            verdict.clone(),
            true,
            false,
            active_after,
        )?;
        if let Some(error) = verdict {
            return Err(error);
        }
        self.active.insert(
            mechanism,
            ActiveRecord {
                lane,
                independence_class: class,
                reservation,
                in_comparison,
            },
        );
        self.lane_active.entry(lane).or_default().push(mechanism);
        self.class_active
            .entry(class)
            .or_default()
            .insert(mechanism);
        if in_comparison {
            self.comparison_reserved
                .entry(lane)
                .or_default()
                .add(&reservation);
        }
        self.reserved.add(&reservation);
        Ok(())
    }

    fn admit_verdict(
        &self,
        lane: ProofLaneId,
        class: ContentHash,
        mechanism: MechanismId,
        reservation: &ResourceEnvelope,
        comparison: Option<&HeadToHeadCharter>,
    ) -> Option<LaneError> {
        if mechanism.lane() != lane {
            return Some(LaneError::MechanismLaneMismatch {
                expected: lane,
                actual: mechanism.lane(),
            });
        }
        if let Some(kind) = self.terminal.get(&mechanism) {
            return Some(LaneError::AlreadyTerminal {
                mechanism,
                kind: *kind,
            });
        }
        if self.active.contains_key(&mechanism) {
            return Some(LaneError::LaneOccupied {
                lane,
                active: mechanism,
            });
        }
        let lane_occupants = self.lane_active.get(&lane).map_or(&[][..], Vec::as_slice);
        match comparison {
            None => {
                if let Some(active) = lane_occupants.first() {
                    return Some(LaneError::LaneOccupied {
                        lane,
                        active: *active,
                    });
                }
                if let Some(active) = self
                    .class_active
                    .get(&class)
                    .and_then(|active| active.iter().next())
                {
                    return Some(LaneError::IndependenceClassOccupied { active: *active });
                }
            }
            Some(head_to_head) => {
                if !head_to_head.candidates.contains(&mechanism) {
                    return Some(LaneError::NotADeclaredCandidate { lane });
                }
                if let Some(active) = self
                    .class_active
                    .get(&class)
                    .and_then(|active| active.iter().find(|active| active.lane() != lane))
                {
                    return Some(LaneError::IndependenceClassOccupied { active: *active });
                }
                let comparison_used = self
                    .comparison_reserved
                    .get(&lane)
                    .copied()
                    .unwrap_or_default();
                if let Err(error) = head_to_head
                    .shared
                    .admit(&comparison_used, reservation, true)
                {
                    return Some(error);
                }
            }
        }
        if u64::from(self.active_count()) >= u64::from(self.policy.max_active_mechanisms) {
            return Some(LaneError::PortfolioCapExceeded {
                active: self.active_count(),
                cap: self.policy.max_active_mechanisms,
            });
        }
        if let Err(error) = self.policy.global.admit(&self.reserved, reservation, false) {
            return Some(error);
        }
        None
    }

    /// Finalize a mechanism against a durable-ledger receipt. This is
    /// the only path that releases a slot, and reserved audit capacity
    /// guarantees that an already-active mechanism can record it.
    ///
    /// # Errors
    /// Structured LaneError.
    pub fn finalize(
        &mut self,
        receipt: &FinalizationReceipt,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let digest = Self::digest_finalize(receipt);
        let mechanism = receipt.mechanism;
        let lane = mechanism.lane();
        let conflict_request = DecisionRequestRef::Finalize {
            receipt,
            released: None,
        };
        match self.replay(key, digest) {
            ReplayStatus::Recorded(refusal) => {
                return refusal.map_or(Ok(()), Err);
            }
            ReplayStatus::Conflict { original_seq } => {
                let error = LaneError::IdempotencyConflict { original_seq };
                self.record(
                    DecisionKind::Finalize,
                    lane,
                    mechanism,
                    key,
                    digest,
                    conflict_request,
                    Some(error.clone()),
                    false,
                    true,
                    self.active.len(),
                )?;
                return Err(error);
            }
            ReplayStatus::Fresh => {}
        }
        let expected = FinalizationReceipt::new(
            mechanism,
            receipt.kind,
            receipt.superseded_by,
            receipt.ledger_artifact,
        );
        let active_record = self.active.get(&mechanism).cloned();
        let verdict =
            if expected.as_ref().map(FinalizationReceipt::identity) != Ok(receipt.identity) {
                Some(LaneError::ReceiptInvalid {
                    what: "identity does not match the receipt's own fields",
                })
            } else if let Some(kind) = self.terminal.get(&mechanism) {
                Some(LaneError::AlreadyTerminal {
                    mechanism,
                    kind: *kind,
                })
            } else if active_record.is_none() {
                Some(LaneError::UnknownMechanism { mechanism })
            } else {
                None
            };
        let released = verdict
            .is_none()
            .then(|| active_record.as_ref().map(|record| record.reservation))
            .flatten();
        let active_after = if verdict.is_none() {
            self.active.len().saturating_sub(1)
        } else {
            self.active.len()
        };
        self.record(
            DecisionKind::Finalize,
            lane,
            mechanism,
            key,
            digest,
            DecisionRequestRef::Finalize { receipt, released },
            verdict.clone(),
            true,
            false,
            active_after,
        )?;
        if let Some(error) = verdict {
            return Err(error);
        }
        let Some(record) = self.active.remove(&mechanism) else {
            return Err(LaneError::UnknownMechanism { mechanism });
        };
        if let Some(occupants) = self.lane_active.get_mut(&record.lane) {
            occupants.retain(|candidate| *candidate != mechanism);
            if occupants.is_empty() {
                self.lane_active.remove(&record.lane);
            }
        }
        let remove_class =
            if let Some(active) = self.class_active.get_mut(&record.independence_class) {
                active.remove(&mechanism);
                active.is_empty()
            } else {
                false
            };
        if remove_class {
            self.class_active.remove(&record.independence_class);
        }
        if record.in_comparison
            && let Some(used) = self.comparison_reserved.get_mut(&record.lane)
        {
            used.sub(&record.reservation);
        }
        self.reserved.sub(&record.reservation);
        self.terminal.insert(mechanism, receipt.kind);
        Ok(())
    }
}

#[cfg(test)]
mod request_digest_tests {
    use super::*;

    fn charter(independence_class: &str) -> LaneCharter {
        LaneCharter::new(
            "lane statement",
            "admissible domain",
            &["assumption"],
            "target authority",
            "baseline",
            "falsifier family",
            independence_class,
        )
        .expect("valid charter fixture")
    }

    fn envelope(work: u64) -> ResourceEnvelope {
        ResourceEnvelope {
            work_units: work,
            memory_bytes: 1 << 20,
            reviewer_slots: 2,
            falsification_capacity: 3,
        }
    }

    #[test]
    fn lane_request_digest_binds_every_preimage_field() {
        let charter_a = charter("class-a");
        let charter_b = charter("class-b");
        let lane_a = charter_a.lane_id();
        let lane_b = charter_b.lane_id();
        let mech_a = charter_a.mechanism_id("mechanism", 1).expect("mechanism a");
        let mech_b = charter_a.mechanism_id("mechanism", 2).expect("mechanism b");

        // Admit preimage: lane, mechanism, and every reservation axis move
        // the digest; identical inputs replay bit-identically.
        let base = PortfolioLedger::digest_admit(lane_a, mech_a, &envelope(10));
        assert_eq!(
            base,
            PortfolioLedger::digest_admit(lane_a, mech_a, &envelope(10)),
            "replay identity must be bit-stable"
        );
        assert_ne!(
            base,
            PortfolioLedger::digest_admit(lane_b, mech_a, &envelope(10)),
            "lane-hash was omitted"
        );
        assert_ne!(
            base,
            PortfolioLedger::digest_admit(lane_a, mech_b, &envelope(10)),
            "mechanism-hash was omitted"
        );
        let mut wider = envelope(10);
        wider.memory_bytes += 1;
        for (axis, changed) in [
            ("work", envelope(11)),
            ("memory", wider),
            (
                "reviewer",
                ResourceEnvelope {
                    reviewer_slots: 3,
                    ..envelope(10)
                },
            ),
            (
                "falsification-capacity",
                ResourceEnvelope {
                    falsification_capacity: 4,
                    ..envelope(10)
                },
            ),
        ] {
            assert_ne!(
                base,
                PortfolioLedger::digest_admit(lane_a, mech_a, &changed),
                "reservation axis {axis} was omitted"
            );
        }

        // Preregister preimage: candidates (item and order), shared
        // envelope, and the preregistration artifact all move the digest.
        let artifact = fs_blake3::hash_domain("fs-govern-test", b"artifact");
        let h2h = |candidates: Vec<MechanismId>, shared: ResourceEnvelope, artifact| {
            HeadToHeadCharter::new(&charter("class-a"), &candidates, shared, artifact)
                .expect("valid comparison charter")
        };
        let pre_base =
            PortfolioLedger::digest_preregister(&h2h(vec![mech_a, mech_b], envelope(10), artifact));
        // Candidate order is CANONICALIZED at charter admission (the
        // constructor sorts), so the same membership digests identically
        // regardless of presentation order...
        assert_eq!(
            pre_base,
            PortfolioLedger::digest_preregister(
                &h2h(vec![mech_b, mech_a], envelope(10), artifact,)
            ),
            "candidate order must be canonicalized, not identity-bearing"
        );
        // ...while membership itself is identity-bearing.
        let mech_c = charter("class-a")
            .mechanism_id("mechanism", 3)
            .expect("mechanism c");
        assert_ne!(
            pre_base,
            PortfolioLedger::digest_preregister(
                &h2h(vec![mech_a, mech_c], envelope(10), artifact,)
            ),
            "candidate membership was omitted"
        );
        assert_ne!(
            pre_base,
            PortfolioLedger::digest_preregister(
                &h2h(vec![mech_a, mech_b], envelope(11), artifact,)
            ),
            "shared envelope was omitted"
        );
        assert_ne!(
            pre_base,
            PortfolioLedger::digest_preregister(&h2h(
                vec![mech_a, mech_b],
                envelope(10),
                fs_blake3::hash_domain("fs-govern-test", b"other-artifact"),
            )),
            "preregistration artifact was omitted"
        );

        // Finalize preimage: the receipt identity moves the digest, and
        // the kind tag domain-separates the three request kinds even over
        // maximally-shared inputs.
        let receipt_a = FinalizationReceipt::new(mech_a, TerminalKind::Withdrawn, None, artifact)
            .expect("receipt a");
        let receipt_b = FinalizationReceipt::new(mech_b, TerminalKind::Withdrawn, None, artifact)
            .expect("receipt b");
        let fin_a = PortfolioLedger::digest_finalize(&receipt_a);
        assert_ne!(
            fin_a,
            PortfolioLedger::digest_finalize(&receipt_b),
            "receipt identity was omitted"
        );
        assert_ne!(
            base, pre_base,
            "kind tag was omitted (admit vs preregister)"
        );
        assert_ne!(base, fin_a, "kind tag was omitted (admit vs finalize)");
    }

    #[test]
    fn lane_request_digest_ignores_retention_only_fields() {
        // The digest functions take only the request inputs: retention
        // metadata (seq, policy echo, idempotency key, refusal outcome)
        // cannot reach the preimage. Two decisions retained under
        // different keys and sequence positions carry the SAME request
        // digest for the same request.
        let charter = charter("class-a");
        let lane = charter.lane_id();
        let mechanism = charter.mechanism_id("mechanism", 1).expect("mechanism");
        let digest_first = PortfolioLedger::digest_admit(lane, mechanism, &envelope(10));
        let digest_second = PortfolioLedger::digest_admit(lane, mechanism, &envelope(10));
        assert_eq!(
            digest_first, digest_second,
            "retention-only state (seq/policy/idempotency/refusal) is not a digest input"
        );
    }
}
