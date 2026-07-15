//! One active unproven mechanism per independently falsifiable proof
//! lane (bead frankensim-ext-epic-gov-rjoq.6).
//!
//! The addendum documents a one-bet discipline; this module makes it
//! EXECUTABLE. A [`LaneCharter`] canonicalizes the semantic fields that
//! define a proof lane (statement/quantifiers, admissible domain,
//! assumptions, target authority, baseline, falsifier family,
//! independence class) and derives an authenticated [`ProofLaneId`] —
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
//! Every method validates completely BEFORE mutating, so a refusal
//! leaves the ledger observably unchanged (no partial admission), and
//! the decision log is a deterministic, bounded, replayable record.

use crate::json_escape;
use fs_blake3::ContentHash;
use std::collections::BTreeMap;

/// Version of the lane-admission policy schema: bump when a rule,
/// canonicalization step, or identity preimage changes meaning.
pub const LANE_POLICY_VERSION: u32 = 1;

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

/// Maximum bytes for one canonical charter field.
pub const MAX_FIELD_BYTES: usize = 4096;

/// Maximum declared assumptions per charter.
pub const MAX_ASSUMPTIONS: usize = 256;

/// Maximum candidates in one preregistered head-to-head comparison.
pub const MAX_H2H_CANDIDATES: usize = 8;

/// Collapse whitespace runs to single spaces and trim — the G3
/// canonicalization that makes cosmetic re-spellings identity-stable.
fn canonical_text(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
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
        }
    }
}

impl std::error::Error for LaneError {}

/// Authenticated identity of one proof lane. Minted ONLY by
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
pub struct MechanismId(ContentHash);

impl MechanismId {
    /// The underlying content hash (read-only).
    #[must_use]
    pub fn as_hash(&self) -> &ContentHash {
        &self.0
    }
}

impl core::fmt::Display for MechanismId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}", self.0)
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
        let field = |what: &'static str, raw: &str| -> Result<String, LaneError> {
            let canonical = canonical_text(raw);
            if canonical.is_empty() {
                return Err(LaneError::EmptyField { what });
            }
            if canonical.len() > MAX_FIELD_BYTES {
                return Err(LaneError::TooLarge {
                    what,
                    len: canonical.len(),
                    cap: MAX_FIELD_BYTES,
                });
            }
            Ok(canonical)
        };
        if assumptions.len() > MAX_ASSUMPTIONS {
            return Err(LaneError::TooLarge {
                what: "assumptions",
                len: assumptions.len(),
                cap: MAX_ASSUMPTIONS,
            });
        }
        let mut canon_assumptions = assumptions
            .iter()
            .map(|a| field("assumption", a))
            .collect::<Result<Vec<_>, _>>()?;
        canon_assumptions.sort_unstable();
        canon_assumptions.dedup();
        Ok(LaneCharter {
            statement: field("statement", statement)?,
            admissible_domain: field("admissible domain", admissible_domain)?,
            assumptions: canon_assumptions,
            target_authority: field("target authority", target_authority)?,
            baseline: field("baseline", baseline)?,
            falsifier_family: field("falsifier family", falsifier_family)?,
            independence_class: field("independence class", independence_class)?,
        })
    }

    /// The authenticated lane identity: domain-separated BLAKE3 over
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
        let canonical_name = canonical_text(name);
        if canonical_name.is_empty() {
            return Err(LaneError::EmptyField {
                what: "mechanism name",
            });
        }
        if canonical_name.len() > MAX_FIELD_BYTES {
            return Err(LaneError::TooLarge {
                what: "mechanism name",
                len: canonical_name.len(),
                cap: MAX_FIELD_BYTES,
            });
        }
        let mut canonical = Vec::new();
        push_field(&mut canonical, 1, self.lane_id().as_hash().as_bytes());
        push_field(&mut canonical, 2, canonical_name.as_bytes());
        push_field(&mut canonical, 3, &version.to_le_bytes());
        Ok(MechanismId(fs_blake3::hash_domain(
            MECHANISM_IDENTITY_DOMAIN,
            &canonical,
        )))
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

/// A preregistered, bounded head-to-head comparison: the ONLY way two
/// active mechanisms may share a lane. Declared BEFORE any admission
/// in the lane, naming its candidates and a shared budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeadToHeadCharter {
    lane: ProofLaneId,
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
        lane: ProofLaneId,
        candidates: &[MechanismId],
        shared: ResourceEnvelope,
        preregistration_artifact: ContentHash,
    ) -> Result<HeadToHeadCharter, LaneError> {
        let mut sorted = candidates.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != candidates.len() || sorted.len() < 2 || sorted.len() > MAX_H2H_CANDIDATES
        {
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
            lane,
            candidates: sorted,
            shared,
            preregistration_artifact,
        })
    }

    /// The lane this comparison governs.
    #[must_use]
    pub fn lane(&self) -> ProofLaneId {
        self.lane
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

/// One atomic decision in the replayable log.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmissionDecision {
    /// Sequence number (0-based, dense).
    pub seq: u64,
    /// Policy schema version in force.
    pub policy_version: u32,
    /// What kind of request this was.
    pub kind: DecisionKind,
    /// The lane.
    pub lane: ProofLaneId,
    /// The mechanism (subject).
    pub mechanism: MechanismId,
    /// The idempotency key presented.
    pub idempotency: IdempotencyKey,
    /// Digest of the complete request (for replay conflict checks).
    pub request_digest: ContentHash,
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
            Some(e) => format!("refused: {e}"),
        };
        let remedy = self
            .refusal
            .as_ref()
            .map_or_else(String::new, |e| e.remedy().to_owned());
        write!(
            out,
            "{{\"seq\":{},\"policy_version\":{},\"kind\":\"{}\",\"lane\":\"{}\",\"mechanism\":\"{}\",\"idempotency\":\"{}\",\"request_digest\":\"{}\",\"verdict\":\"{}\",\"remedy\":\"{}\"}}",
            self.seq,
            self.policy_version,
            self.kind.name(),
            self.lane,
            self.mechanism,
            self.idempotency,
            self.request_digest,
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

/// The atomic admission state machine. Exclusive access (`&mut self`)
/// is the concurrency contract: every method validates COMPLETELY
/// before mutating, so refused or replayed requests leave state
/// unchanged and interleaved retries cannot oversubscribe.
#[derive(Debug, Clone, PartialEq)]
pub struct PortfolioLedger {
    policy: PortfolioPolicy,
    active: BTreeMap<MechanismId, ActiveRecord>,
    lane_active: BTreeMap<ProofLaneId, Vec<MechanismId>>,
    class_active: BTreeMap<ContentHash, MechanismId>,
    comparisons: BTreeMap<ProofLaneId, HeadToHeadCharter>,
    comparison_reserved: BTreeMap<ProofLaneId, ResourceEnvelope>,
    terminal: BTreeMap<MechanismId, TerminalKind>,
    reserved: ResourceEnvelope,
    decisions: Vec<AdmissionDecision>,
    idempotency: BTreeMap<IdempotencyKey, u64>,
}

impl PortfolioLedger {
    /// Empty ledger under `policy`.
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
            idempotency: BTreeMap::new(),
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

    /// The full deterministic decision log.
    #[must_use]
    pub fn decisions(&self) -> &[AdmissionDecision] {
        &self.decisions
    }

    /// Bounded JSON decision log: at most `limit` most-recent rows plus
    /// an explicit truncation count (never a silent cap).
    #[must_use]
    pub fn decisions_json(&self, limit: usize) -> String {
        use core::fmt::Write as _;
        let skipped = self.decisions.len().saturating_sub(limit);
        let mut out = format!("{{\"skipped\":{skipped},\"decisions\":[");
        for (i, d) in self.decisions.iter().skip(skipped).enumerate() {
            if i > 0 {
                out.push(',');
            }
            write!(out, "{}", d.to_json()).expect("writing to a String is infallible");
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
        push_field(&mut canonical, 2, charter.lane.as_hash().as_bytes());
        for c in &charter.candidates {
            push_field(&mut canonical, 3, c.as_hash().as_bytes());
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

    /// Idempotency gate: `Ok(Some(..))` replays the recorded verdict
    /// for a byte-identical request; `Err` refuses a different request
    /// under a used key; `Ok(None)` means the key is fresh.
    fn replay(
        &self,
        key: IdempotencyKey,
        request_digest: ContentHash,
    ) -> Result<Option<&AdmissionDecision>, LaneError> {
        match self.idempotency.get(&key) {
            None => Ok(None),
            Some(seq) => {
                let recorded = &self.decisions[usize::try_from(*seq).expect("seq fits usize")];
                if recorded.request_digest == request_digest {
                    Ok(Some(recorded))
                } else {
                    Err(LaneError::IdempotencyConflict { original_seq: *seq })
                }
            }
        }
    }

    fn record(
        &mut self,
        kind: DecisionKind,
        lane: ProofLaneId,
        mechanism: MechanismId,
        key: IdempotencyKey,
        request_digest: ContentHash,
        refusal: Option<LaneError>,
    ) -> u64 {
        let seq = u64::try_from(self.decisions.len()).expect("decision count fits u64");
        self.decisions.push(AdmissionDecision {
            seq,
            policy_version: LANE_POLICY_VERSION,
            kind,
            lane,
            mechanism,
            idempotency: key,
            request_digest,
            refusal,
        });
        self.idempotency.insert(key, seq);
        seq
    }

    /// Preregister a bounded head-to-head comparison for a lane —
    /// BEFORE any admission in that lane, at most one per lane.
    ///
    /// # Errors
    /// Structured [`LaneError`]; the refusal is also recorded in the
    /// decision log under the presented idempotency key.
    pub fn preregister_comparison(
        &mut self,
        charter: HeadToHeadCharter,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let digest = Self::digest_preregister(&charter);
        if let Some(recorded) = self.replay(key, digest)? {
            return match &recorded.refusal {
                None => Ok(()),
                Some(e) => Err(e.clone()),
            };
        }
        let lane = charter.lane;
        let subject = charter.candidates[0];
        let verdict = if self.comparisons.contains_key(&lane) {
            Some(LaneError::ComparisonAlreadyDeclared { lane })
        } else if self.lane_active.get(&lane).is_some_and(|v| !v.is_empty()) {
            Some(LaneError::ComparisonAfterAdmission { lane })
        } else {
            None
        };
        let refused = verdict.clone();
        self.record(
            DecisionKind::Preregister,
            lane,
            subject,
            key,
            digest,
            verdict,
        );
        match refused {
            None => {
                self.comparison_reserved
                    .insert(lane, ResourceEnvelope::default());
                self.comparisons.insert(lane, charter);
                Ok(())
            }
            Some(e) => Err(e),
        }
    }

    /// Atomically admit one mechanism into its lane. Validation order
    /// is deterministic: idempotency, terminal permanence, duplicate
    /// activity, lane occupancy (with the preregistered-comparison
    /// carve-out), independence-class collision, portfolio cap,
    /// comparison envelope, global envelope. Any refusal leaves the
    /// ledger unchanged except for the recorded decision.
    ///
    /// # Errors
    /// Structured [`LaneError`] with a ranked remedy.
    pub fn admit(
        &mut self,
        charter: &LaneCharter,
        mechanism: MechanismId,
        reservation: ResourceEnvelope,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let lane = charter.lane_id();
        let digest = Self::digest_admit(lane, mechanism, &reservation);
        if let Some(recorded) = self.replay(key, digest)? {
            return match &recorded.refusal {
                None => Ok(()),
                Some(e) => Err(e.clone()),
            };
        }
        let class = charter.independence_class_id();
        let comparison = self.comparisons.get(&lane);
        let in_comparison = comparison.is_some();
        let verdict = self.admit_verdict(lane, class, mechanism, &reservation, comparison);
        let refused = verdict.clone();
        self.record(DecisionKind::Admit, lane, mechanism, key, digest, verdict);
        if let Some(e) = refused {
            return Err(e);
        }
        // Commit — every check passed; mutations are now unconditional.
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
        self.class_active.entry(class).or_insert(mechanism);
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
                // The split-gaming backstop: an ACTIVE bet in the same
                // declared independence class (necessarily another
                // lane here) blocks this one.
                if let Some(active) = self.class_active.get(&class) {
                    return Some(LaneError::IndependenceClassOccupied { active: *active });
                }
            }
            Some(h2h) => {
                if !h2h.candidates.contains(&mechanism) {
                    return Some(LaneError::NotADeclaredCandidate { lane });
                }
                // A comparison licenses multiple bets INSIDE its lane
                // only — an active same-class bet in a DIFFERENT lane
                // still blocks (the backstop cannot be evaded by
                // preregistering a comparison elsewhere).
                if let Some(active) = self.class_active.get(&class)
                    && self.active.get(active).is_some_and(|r| r.lane != lane)
                {
                    return Some(LaneError::IndependenceClassOccupied { active: *active });
                }
                let comparison_used = self
                    .comparison_reserved
                    .get(&lane)
                    .copied()
                    .unwrap_or_default();
                if let Err(e) = h2h.shared.admit(&comparison_used, reservation, true) {
                    return Some(e);
                }
            }
        }
        if u64::from(self.active_count()) >= u64::from(self.policy.max_active_mechanisms) {
            return Some(LaneError::PortfolioCapExceeded {
                active: self.active_count(),
                cap: self.policy.max_active_mechanisms,
            });
        }
        if let Err(e) = self.policy.global.admit(&self.reserved, reservation, false) {
            return Some(e);
        }
        None
    }

    /// Finalize a mechanism against a durable ledger receipt: the ONLY
    /// path that releases a slot, and it releases exactly once. There
    /// is deliberately no timeout/stall path — Unknown never releases.
    ///
    /// # Errors
    /// Structured [`LaneError`].
    pub fn finalize(
        &mut self,
        receipt: &FinalizationReceipt,
        key: IdempotencyKey,
    ) -> Result<(), LaneError> {
        let digest = Self::digest_finalize(receipt);
        if let Some(recorded) = self.replay(key, digest)? {
            return match &recorded.refusal {
                None => Ok(()),
                Some(e) => Err(e.clone()),
            };
        }
        let mechanism = receipt.mechanism;
        let expected = FinalizationReceipt::new(
            mechanism,
            receipt.kind,
            receipt.superseded_by,
            receipt.ledger_artifact,
        );
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
            } else if !self.active.contains_key(&mechanism) {
                Some(LaneError::UnknownMechanism { mechanism })
            } else {
                None
            };
        // For a finalize refused before any lane is known (unknown
        // mechanism), the log row's lane column carries the receipt
        // identity — a deterministic placeholder that cannot collide
        // with a real lane id (different hash domain).
        let lane = self
            .active
            .get(&mechanism)
            .map_or_else(|| ProofLaneId(receipt.identity), |r| r.lane);
        let refused = verdict.clone();
        self.record(
            DecisionKind::Finalize,
            lane,
            mechanism,
            key,
            digest,
            verdict,
        );
        if let Some(e) = refused {
            return Err(e);
        }
        let record = self
            .active
            .remove(&mechanism)
            .expect("presence checked before commit");
        if let Some(occupants) = self.lane_active.get_mut(&record.lane) {
            occupants.retain(|m| *m != mechanism);
            if occupants.is_empty() {
                self.lane_active.remove(&record.lane);
            }
        }
        if self.class_active.get(&record.independence_class) == Some(&mechanism) {
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
