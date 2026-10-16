//! VerificationManifest v1 identity core (bead
//! frankensim-leapfrog-2026-program-i94v.7.1.1): stable semantic identity
//! and source authority, so renamed tests, revised Beads, reorganized
//! scripts, or changed contracts can never sever evidence from claims or
//! retroactively reinterpret an old receipt.
//!
//! [`ClaimId`] names a continuing conceptual lineage; [`ClaimRevisionId`]
//! is CONTENT-ADDRESSED over one exact statement (kind, quantifiers,
//! units/conventions, hypotheses, domain, code/contract surface, no-claim
//! boundary), so distinct revisions cannot collide and identical content
//! is idempotent. Supersession appends a new revision pointing at the
//! old; nothing ever mutates an old revision.
//!
//! [`ClaimRelationReceipt`] types the admitted claim-graph edges
//! (implication, refinement, restriction, counterexample, certified
//! equivalence) with checker/TCB, quantifier variance, and policy version.
//! Directed kinds preserve their endpoint orientation; certified equivalence
//! normalizes its bidirectional endpoints into revision-id order. Promotion
//! transfers ONLY across admitted relations;
//! directed cycles refuse unless every member is joined by certified
//! equivalence, in which case the strongly connected component
//! canonicalizes to one representative WITHOUT erasing members.
//!
//! [`NormalizedGraph`] is the one canonical normalized form: equivalent
//! semantic manifests normalize to the same domain-separated digest, and
//! the human/JSON/ledger renderings are tested SEMANTIC projections of
//! it, never independent truths.
//!
//! No-claims: the manifest is metadata and obligation authority, not
//! proof — compiling one never invents a scientific adjudication; the
//! frozen inventory compiler is V.1.2 scope, the lint battery V.1.3, and
//! ledger persistence is fs-obs/fs-ledger scope.

use core::fmt::{self, Write as _};
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, DomainHasher};

/// Version stamped into every v1 record and digest domain.
pub const MANIFEST_V1_SCHEMA_VERSION: u32 = 1;

const REVISION_DOMAIN: &str = "org.frankensim.fs-vmanifest.claim-revision.v1";
const GRAPH_DOMAIN: &str = "org.frankensim.fs-vmanifest.manifest-graph.v1";
/// Maximum UTF-8 bytes in one v1 semantic text field.
pub const MAX_V1_TEXT_BYTES: usize = 4096;
/// Maximum ASCII bytes in one admitted v1 lineage/case/journey id.
pub const MAX_V1_ID_BYTES: usize = 128;
const V1_LOGICAL_BYTES_PER_REVISION: u64 = 256;
const V1_LOGICAL_BYTES_PER_RECEIPT: u64 = 512;
const V1_PROJECTION_FIXED_BYTES_PER_ROW: u64 = 1024;
const V1_PROJECTION_ESCAPE_EXPANSION: u64 = 6;

/// Version-1 resource ceilings for claim-graph admission.
///
/// Callers may tighten any field through [`admit_graph_with_limits`], but may
/// not loosen these protocol ceilings. Byte fields are deterministic logical
/// charges, not allocator metadata or resident-set measurements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct V1AdmissionLimits {
    /// Maximum exact claim revisions.
    pub max_revisions: u32,
    /// Maximum relation receipts.
    pub max_receipts: u32,
    /// Maximum UTF-8 bytes across claim ids, revision semantics, and receipt
    /// checker/TCB/domain text.
    pub max_semantic_bytes: u64,
    /// Maximum conservative logical graph-construction bytes.
    pub max_working_bytes: u64,
    /// Maximum rows in any complete graph projection, including its header.
    pub max_projection_rows: u32,
    /// Maximum conservative bytes in any one human, JSON-lines, or ledger
    /// projection.
    pub max_projection_bytes: u64,
}

impl V1AdmissionLimits {
    /// The fixed VerificationManifest-v1 protocol envelope.
    pub const DEFAULT: Self = Self {
        max_revisions: 4096,
        max_receipts: 8192,
        max_semantic_bytes: 16 * 1024 * 1024,
        max_working_bytes: 32 * 1024 * 1024,
        max_projection_rows: 1 + 4096 + 8192,
        max_projection_bytes: 128 * 1024 * 1024,
    };

    fn validate_protocol_ceiling(self) -> Result<(), V1Error> {
        for (quantity, required, admitted, unit) in [
            (
                "configured revision ceiling",
                u64::from(self.max_revisions),
                u64::from(Self::DEFAULT.max_revisions),
                "revisions",
            ),
            (
                "configured receipt ceiling",
                u64::from(self.max_receipts),
                u64::from(Self::DEFAULT.max_receipts),
                "receipts",
            ),
            (
                "configured semantic-byte ceiling",
                self.max_semantic_bytes,
                Self::DEFAULT.max_semantic_bytes,
                "bytes",
            ),
            (
                "configured working-byte ceiling",
                self.max_working_bytes,
                Self::DEFAULT.max_working_bytes,
                "logical bytes",
            ),
            (
                "configured projection-row ceiling",
                u64::from(self.max_projection_rows),
                u64::from(Self::DEFAULT.max_projection_rows),
                "rows",
            ),
            (
                "configured projection-byte ceiling",
                self.max_projection_bytes,
                Self::DEFAULT.max_projection_bytes,
                "logical bytes",
            ),
        ] {
            if required > admitted {
                return Err(resource_refusal(quantity, required, admitted, unit));
            }
        }
        Ok(())
    }
}

impl Default for V1AdmissionLimits {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Deterministic resource charges retained by an admitted v1 graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct V1ResourceEnvelope {
    revision_count: u32,
    receipt_count: u32,
    semantic_bytes: u64,
    logical_working_bytes: u64,
    projection_rows: u32,
    projection_bytes_upper_bound: u64,
}

impl V1ResourceEnvelope {
    /// Exact admitted revision count.
    #[must_use]
    pub const fn revision_count(self) -> u32 {
        self.revision_count
    }

    /// Exact admitted receipt count.
    #[must_use]
    pub const fn receipt_count(self) -> u32 {
        self.receipt_count
    }

    /// Exact aggregate semantic UTF-8 bytes charged at admission.
    #[must_use]
    pub const fn semantic_bytes(self) -> u64 {
        self.semantic_bytes
    }

    /// Conservative logical graph-construction bytes.
    #[must_use]
    pub const fn logical_working_bytes(self) -> u64 {
        self.logical_working_bytes
    }

    /// Exact rows in each complete projection, including its header.
    #[must_use]
    pub const fn projection_rows(self) -> u32 {
        self.projection_rows
    }

    /// Conservative bytes for any one complete projection.
    #[must_use]
    pub const fn projection_bytes_upper_bound(self) -> u64 {
        self.projection_bytes_upper_bound
    }
}

/// Machine-readable resource-limit evidence carried by a v1 refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct V1ResourceRefusal {
    /// Stable quantity name.
    pub quantity: &'static str,
    /// Exact amount required by the rejected input or configured ceiling.
    pub required: u64,
    /// Amount admitted by the active v1 envelope.
    pub admitted: u64,
    /// Stable unit label for both numeric values.
    pub unit: &'static str,
}

/// A typed refusal with a stable rule slug and RANKED candidate fixes —
/// the "stable diagnostics and ranked fixes" the success criteria demand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct V1Error {
    rule: &'static str,
    detail: String,
    fixes: Vec<String>,
    resource: Option<V1ResourceRefusal>,
}

impl V1Error {
    fn new(rule: &'static str, detail: impl Into<String>) -> V1Error {
        V1Error {
            rule,
            detail: detail.into(),
            fixes: Vec::new(),
            resource: None,
        }
    }

    fn with_fix(mut self, fix: impl Into<String>) -> V1Error {
        self.fixes.push(fix.into());
        self
    }

    fn with_resource(mut self, resource: V1ResourceRefusal) -> V1Error {
        self.resource = Some(resource);
        self
    }

    /// Crate-internal constructor with a ranked fix list.
    pub(crate) fn with_fixes(
        rule: &'static str,
        detail: impl Into<String>,
        fixes: &[&str],
    ) -> V1Error {
        V1Error {
            rule,
            detail: detail.into(),
            fixes: fixes.iter().map(|f| (*f).to_owned()).collect(),
            resource: None,
        }
    }

    /// The stable rule slug.
    #[must_use]
    pub const fn rule(&self) -> &'static str {
        self.rule
    }

    /// Human-readable detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }

    /// Ranked candidate fixes, most likely first.
    #[must_use]
    pub fn ranked_fixes(&self) -> &[String] {
        &self.fixes
    }

    /// Machine-readable required/admitted/unit evidence for resource refusals.
    #[must_use]
    pub const fn resource_refusal(&self) -> Option<V1ResourceRefusal> {
        self.resource
    }
}

impl fmt::Display for V1Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule, self.detail)?;
        for (i, fix) in self.fixes.iter().enumerate() {
            write!(f, "\n  fix[{i}]: {fix}")?;
        }
        Ok(())
    }
}

impl std::error::Error for V1Error {}

fn resource_refusal(
    quantity: &'static str,
    required: u64,
    admitted: u64,
    unit: &'static str,
) -> V1Error {
    V1Error::new(
        "v1-resource-limit",
        format!("{quantity} requires {required} {unit}, admitted {admitted} {unit}"),
    )
    .with_resource(V1ResourceRefusal {
        quantity,
        required,
        admitted,
        unit,
    })
    .with_fix(format!(
        "reduce {quantity} to at most {admitted} {unit} or split the graph"
    ))
}

fn resource_overflow(quantity: &'static str, unit: &'static str) -> V1Error {
    V1Error::new(
        "v1-resource-overflow",
        format!("{quantity} cannot be represented in v1 {unit}"),
    )
    .with_fix("reduce the graph before admission; wrapped resource accounting is forbidden")
}

fn allocation_refusal(payload: &'static str, required_bytes: u64) -> V1Error {
    V1Error::new(
        "v1-allocation-refused",
        format!(
            "{payload} requires {required_bytes} admitted logical bytes, but the allocator refused the reservation"
        ),
    )
    .with_fix("release memory pressure and retry the same bounded input")
}

fn checked_add_resource(left: u64, right: u64, quantity: &'static str) -> Result<u64, V1Error> {
    left.checked_add(right)
        .ok_or_else(|| resource_overflow(quantity, "u64 bytes"))
}

fn checked_mul_resource(left: u64, right: u64, quantity: &'static str) -> Result<u64, V1Error> {
    left.checked_mul(right)
        .ok_or_else(|| resource_overflow(quantity, "u64 bytes"))
}

fn checked_u32_length(field: &'static str, length: usize) -> Result<u32, V1Error> {
    u32::try_from(length).map_err(|_| {
        V1Error::new(
            "v1-identity-length-overflow",
            format!("{field} length {length} cannot be encoded as a v1 u32 length"),
        )
        .with_fix("admit a bounded value before requesting authoritative identity")
    })
}

fn hash_framed_text(
    hasher: &mut DomainHasher,
    field: &'static str,
    value: &str,
) -> Result<(), V1Error> {
    hasher.update(&checked_u32_length(field, value.len())?.to_be_bytes());
    hasher.update(value.as_bytes());
    Ok(())
}

fn bounded_id(kind: &'static str, value: &str) -> Result<String, V1Error> {
    if value.is_empty() || value.len() > MAX_V1_ID_BYTES {
        return Err(V1Error::new(
            "v1-id-bounds",
            format!(
                "{kind} id length {} outside 1..={MAX_V1_ID_BYTES}",
                value.len()
            ),
        )
        .with_fix(format!("shorten or supply a non-empty {kind} id")));
    }
    if !value
        .bytes()
        .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b"-_.:/".contains(&b))
    {
        return Err(V1Error::new(
            "v1-id-bounds",
            format!("{kind} id {value:?} outside [a-z0-9-_.:/]"),
        )
        .with_fix("kebab-case the id; identity must not depend on display casing"));
    }
    Ok(value.to_owned())
}

fn bounded_text(field: &'static str, value: &str) -> Result<(), V1Error> {
    if value.is_empty() || value.len() > MAX_V1_TEXT_BYTES {
        return Err(V1Error::new(
            "v1-field-bounds",
            format!(
                "{field} length {} outside 1..={MAX_V1_TEXT_BYTES}",
                value.len()
            ),
        )
        .with_fix(format!(
            "populate {field}; empty semantic fields sever evidence from claims"
        )));
    }
    Ok(())
}

fn escape_json_string(value: &str) -> String {
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
                write!(&mut escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to a String cannot fail");
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

/// The continuing conceptual lineage of one claim.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ClaimId(String);

impl ClaimId {
    /// Admit a claim lineage id.
    pub fn new(value: &str) -> Result<ClaimId, V1Error> {
        Ok(ClaimId(bounded_id("claim", value)?))
    }

    /// The id text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Stable identity of one conformance case.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CaseId(String);

impl CaseId {
    /// Admit a case id.
    pub fn new(value: &str) -> Result<CaseId, V1Error> {
        Ok(CaseId(bounded_id("case", value)?))
    }

    /// The id text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Stable identity of one verification journey.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct JourneyId(String);

impl JourneyId {
    /// Admit a journey id.
    pub fn new(value: &str) -> Result<JourneyId, V1Error> {
        Ok(JourneyId(bounded_id("journey", value)?))
    }

    /// The id text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// The kind of statement one revision binds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum ClaimKind {
    /// A behavioral/functional property.
    Behavioral,
    /// A quantitative bound with units.
    QuantitativeBound,
    /// A determinism/replay property.
    Determinism,
    /// A refusal/fail-closed property.
    Refusal,
    /// A theorem-backed property.
    Theorem,
}

const fn claim_kind_code(kind: ClaimKind) -> u8 {
    match kind {
        ClaimKind::Behavioral => 1,
        ClaimKind::QuantitativeBound => 2,
        ClaimKind::Determinism => 3,
        ClaimKind::Refusal => 4,
        ClaimKind::Theorem => 5,
    }
}

/// One exact, immutable claim statement. The identity is content-
/// addressed over every semantic field, so revising ANY of them is a new
/// revision, and an old receipt can never be retroactively reinterpreted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimRevision {
    /// The lineage this revision belongs to.
    pub claim: ClaimId,
    /// The statement kind.
    pub kind: ClaimKind,
    /// The exact statement text.
    pub statement: String,
    /// Quantifier structure ("for all admitted meshes", "exists a
    /// certified witness", ...).
    pub quantifiers: String,
    /// Units and conventions the statement is read under.
    pub units_conventions: String,
    /// Hypotheses/preconditions.
    pub hypotheses: String,
    /// The validity domain.
    pub domain: String,
    /// The exact code/contract surface the claim binds
    /// (crate::path/CONTRACT section).
    pub surface: String,
    /// The no-claim boundary: what this revision deliberately does NOT
    /// claim.
    pub no_claim: String,
    /// The revision this one supersedes, if any. Supersession appends;
    /// the old revision is never mutated or erased.
    pub supersedes: Option<ClaimRevisionId>,
}

/// Content-addressed identity of one exact claim revision.
pub type ClaimRevisionId = ContentHash;

impl ClaimRevision {
    /// Validate every semantic field (empty fields sever evidence).
    pub fn validate(&self) -> Result<(), V1Error> {
        bounded_text("statement", &self.statement)?;
        bounded_text("quantifiers", &self.quantifiers)?;
        bounded_text("units_conventions", &self.units_conventions)?;
        bounded_text("hypotheses", &self.hypotheses)?;
        bounded_text("domain", &self.domain)?;
        bounded_text("surface", &self.surface)?;
        bounded_text("no_claim", &self.no_claim)?;
        Ok(())
    }

    fn revision_id_after_validation(&self) -> Result<ClaimRevisionId, V1Error> {
        let mut hasher = DomainHasher::new(REVISION_DOMAIN);
        hasher.update(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
        hash_framed_text(&mut hasher, "claim", self.claim.as_str())?;
        hasher.update(&[claim_kind_code(self.kind)]);
        hash_framed_text(&mut hasher, "statement", &self.statement)?;
        hash_framed_text(&mut hasher, "quantifiers", &self.quantifiers)?;
        hash_framed_text(&mut hasher, "units_conventions", &self.units_conventions)?;
        hash_framed_text(&mut hasher, "hypotheses", &self.hypotheses)?;
        hash_framed_text(&mut hasher, "domain", &self.domain)?;
        hash_framed_text(&mut hasher, "surface", &self.surface)?;
        hash_framed_text(&mut hasher, "no_claim", &self.no_claim)?;
        match &self.supersedes {
            None => hasher.update(&[0]),
            Some(prev) => {
                hasher.update(&[1]);
                hasher.update(prev.as_bytes());
            }
        }
        Ok(hasher.finalize())
    }

    /// The authoritative content-addressed revision identity: a streamed,
    /// domain-separated hash over every validated semantic field in fixed
    /// order. Invalid raw drafts refuse and cannot mint a revision id.
    pub fn revision_id(&self) -> Result<ClaimRevisionId, V1Error> {
        self.validate()?;
        self.revision_id_after_validation()
    }
}

/// The source-authority lattice: when sources conflict, the HIGHER
/// authority wins only by producing a NEW snapshot/revision — never by
/// silently reinterpreting the lower one. Order is total and explicit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum SourceAuthority {
    /// A script or generated artifact (lowest: fully derived).
    GeneratedArtifact,
    /// A test file or fixture in the tree.
    TestSource,
    /// The crate CONTRACT.md.
    Contract,
    /// The Bead obligation record.
    BeadObligation,
    /// A frozen manifest snapshot (highest: immutable, hash-pinned).
    FrozenSnapshot,
}

/// One pinned source: what it is, its authority, and the exact snapshot
/// hash it was read at. A mutated source is a NEW snapshot, never an
/// in-place change.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourcePin {
    /// Repo-relative path or logical name.
    pub source: String,
    /// The authority level.
    pub authority: SourceAuthority,
    /// The exact content hash the record was compiled against.
    pub snapshot: ContentHash,
}

/// The typed relation kinds the claim graph admits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum RelationKind {
    /// `from` implies `to` (promotion transfers forward).
    Implication,
    /// `from` refines `to` (tighter statement, same lineage direction).
    Refinement,
    /// `from` restricted to a subdomain yields `to`.
    Restriction,
    /// `from` is a counterexample against `to` (transfers REFUTATION,
    /// never promotion).
    Counterexample,
    /// Certified equivalence (bidirectional; the only admitted cycle
    /// former).
    CertifiedEquivalence,
}

const fn relation_code(kind: RelationKind) -> u8 {
    match kind {
        RelationKind::Implication => 1,
        RelationKind::Refinement => 2,
        RelationKind::Restriction => 3,
        RelationKind::Counterexample => 4,
        RelationKind::CertifiedEquivalence => 5,
    }
}

const fn relation_name(kind: RelationKind) -> &'static str {
    match kind {
        RelationKind::Implication => "Implication",
        RelationKind::Refinement => "Refinement",
        RelationKind::Restriction => "Restriction",
        RelationKind::Counterexample => "Counterexample",
        RelationKind::CertifiedEquivalence => "CertifiedEquivalence",
    }
}

const fn relation_arrow(kind: RelationKind) -> &'static str {
    match kind {
        RelationKind::CertifiedEquivalence => "<->",
        RelationKind::Implication
        | RelationKind::Refinement
        | RelationKind::Restriction
        | RelationKind::Counterexample => "->",
    }
}

/// How the quantifier structure varies across a relation edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum QuantifierVariance {
    /// Identical quantifier structure.
    Preserved,
    /// The target claim is quantifier-weaker (admissible for transfer).
    Weakened,
    /// The target claim is quantifier-STRONGER: promotion transfer along
    /// this edge is unsound and refuses.
    Strengthened,
}

const fn variance_code(variance: QuantifierVariance) -> u8 {
    match variance {
        QuantifierVariance::Preserved => 1,
        QuantifierVariance::Weakened => 2,
        QuantifierVariance::Strengthened => 3,
    }
}

const fn variance_name(variance: QuantifierVariance) -> &'static str {
    match variance {
        QuantifierVariance::Preserved => "Preserved",
        QuantifierVariance::Weakened => "Weakened",
        QuantifierVariance::Strengthened => "Strengthened",
    }
}

/// One typed relation draft between exact claim revisions. Structural equality
/// compares the caller-authored draft; it is NOT canonical graph identity.
/// Admission clones the value and normalizes certified-equivalence endpoints.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimRelationReceipt {
    /// The edge kind.
    pub kind: RelationKind,
    /// Source revision. Normalized certified-equivalence receipts expose the
    /// smaller revision id here; directed kinds preserve the supplied source.
    pub from: ClaimRevisionId,
    /// Target revision. Normalized certified-equivalence receipts expose the
    /// larger revision id here; directed kinds preserve the supplied target.
    pub to: ClaimRevisionId,
    /// The proof/checker identity that certified this edge (1..=4096
    /// UTF-8 bytes at admission).
    pub checker: String,
    /// The trusted computing base the checker ran under (1..=4096 UTF-8
    /// bytes at admission).
    pub tcb: String,
    /// Quantifier variance across the edge.
    pub variance: QuantifierVariance,
    /// Domain relationship note (subdomain, identical, disjoint —
    /// nonempty free text bound into the receipt identity, capped at 4096
    /// UTF-8 bytes at admission).
    pub domain_note: String,
    /// The relation-policy version this edge was admitted under.
    pub policy_version: u32,
}

impl ClaimRelationReceipt {
    /// Whether PROMOTION may transfer along this edge (fallback/refutation
    /// rules are the mirror image and live with the consuming gate).
    /// Counterexample edges never transfer promotion; strengthened
    /// quantifiers never transfer promotion. This raw-draft convenience does
    /// not imply that the receipt will admit; certified equivalence refuses
    /// unless its variance is [`QuantifierVariance::Preserved`].
    #[must_use]
    pub fn promotion_transfers(&self) -> bool {
        !matches!(self.kind, RelationKind::Counterexample)
            && !matches!(self.variance, QuantifierVariance::Strengthened)
    }
}

/// The canonical normalized claim graph: sorted revisions, sorted typed
/// edges, certified-equivalence components canonicalized to
/// representatives WITHOUT erasing members. Equivalent semantic
/// manifests normalize to the same digest.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedGraph {
    /// Every admitted revision id, ascending.
    revisions: Vec<ClaimRevisionId>,
    /// Every admitted edge, canonically sorted.
    edges: Vec<ClaimRelationReceipt>,
    /// Certified-equivalence component representative per member
    /// (identity for revisions outside any component).
    representatives: BTreeMap<ClaimRevisionId, ClaimRevisionId>,
    /// Cached streamed identity over the canonical graph.
    digest: ContentHash,
    /// Deterministic admission charges (not identity-forming).
    resources: V1ResourceEnvelope,
}

impl NormalizedGraph {
    /// Every admitted revision id in canonical ascending order.
    #[must_use]
    pub fn revisions(&self) -> &[ClaimRevisionId] {
        &self.revisions
    }

    /// Every admitted relation receipt in its complete canonical order.
    #[must_use]
    pub fn edges(&self) -> &[ClaimRelationReceipt] {
        &self.edges
    }

    /// The immutable representative map for certified-equivalence
    /// components. Every admitted revision has exactly one entry.
    #[must_use]
    pub fn representatives(&self) -> &BTreeMap<ClaimRevisionId, ClaimRevisionId> {
        &self.representatives
    }

    /// Return the canonical representative for an admitted revision.
    #[must_use]
    pub fn representative_of(&self, revision: &ClaimRevisionId) -> Option<ClaimRevisionId> {
        self.representatives.get(revision).copied()
    }

    /// The resource envelope retained from successful admission.
    #[must_use]
    pub const fn resource_envelope(&self) -> V1ResourceEnvelope {
        self.resources
    }

    /// The cached domain-separated digest of the normalized graph.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        self.digest
    }

    /// Human projection: one physical line per revision and edge. Every
    /// digest-forming field is present; untrusted text is escaped so it
    /// cannot forge a line or hide a field.
    #[must_use]
    pub fn render_human(&self) -> String {
        let digest = self.digest().to_hex();
        let mut out = format!(
            "VerificationManifest schema_version={} graph revisions={} edges={} digest={} ordinal=0\n",
            MANIFEST_V1_SCHEMA_VERSION,
            self.revisions.len(),
            self.edges.len(),
            digest
        );
        for (index, revision) in self.revisions.iter().enumerate() {
            let representative = self.representatives[revision];
            out.push_str(&format!(
                "  revision schema_version={} graph_digest={} ordinal={} revision={} representative={}\n",
                MANIFEST_V1_SCHEMA_VERSION,
                digest,
                index + 1,
                revision.to_hex(),
                representative.to_hex(),
            ));
        }
        for (index, edge) in self.edges.iter().enumerate() {
            out.push_str(&format!(
                "  edge schema_version={} graph_digest={} ordinal={} kind={} from={} {} to={} checker=\"{}\" tcb=\"{}\" variance={} domain_note=\"{}\" policy_version={}\n",
                MANIFEST_V1_SCHEMA_VERSION,
                digest,
                self.revisions.len() + index + 1,
                relation_name(edge.kind),
                edge.from.to_hex(),
                relation_arrow(edge.kind),
                edge.to.to_hex(),
                escape_json_string(&edge.checker),
                escape_json_string(&edge.tcb),
                variance_name(edge.variance),
                escape_json_string(&edge.domain_note),
                edge.policy_version,
            ));
        }
        out
    }

    /// Strict JSON-lines projection. The header, every revision, and every
    /// edge carry the schema version, graph digest, and a global canonical
    /// ordinal. Revision records expose the representative map; edge records
    /// expose every digest-forming receipt field.
    #[must_use]
    pub fn render_json_lines(&self) -> Vec<String> {
        let digest = self.digest().to_hex();
        let mut out = vec![format!(
            "{{\"vmanifest\":\"graph\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":0,\"revisions\":{},\"edges\":{}}}",
            MANIFEST_V1_SCHEMA_VERSION,
            digest,
            self.revisions.len(),
            self.edges.len(),
        )];
        for (index, revision) in self.revisions.iter().enumerate() {
            let representative = self.representatives[revision];
            out.push(format!(
                "{{\"vmanifest\":\"revision\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":{},\"revision\":\"{}\",\"representative\":\"{}\"}}",
                MANIFEST_V1_SCHEMA_VERSION,
                digest,
                index + 1,
                revision.to_hex(),
                representative.to_hex(),
            ));
        }
        for (index, edge) in self.edges.iter().enumerate() {
            out.push(format!(
                "{{\"vmanifest\":\"edge\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":{},\"kind\":\"{}\",\"from\":\"{}\",\"to\":\"{}\",\"checker\":\"{}\",\"tcb\":\"{}\",\"variance\":\"{}\",\"domain_note\":\"{}\",\"policy_version\":{}}}",
                MANIFEST_V1_SCHEMA_VERSION,
                digest,
                self.revisions.len() + index + 1,
                relation_name(edge.kind),
                edge.from.to_hex(),
                edge.to.to_hex(),
                escape_json_string(&edge.checker),
                escape_json_string(&edge.tcb),
                variance_name(edge.variance),
                escape_json_string(&edge.domain_note),
                edge.policy_version,
            ));
        }
        out
    }

    /// Ledger-row projection: (row kind, digest-bound payload) tuples.
    #[must_use]
    pub fn render_ledger_rows(&self) -> Vec<(String, String)> {
        let digest = self.digest().to_hex();
        let mut rows = vec![(
            "graph".to_owned(),
            format!(
                "schema_version={} ordinal=0 revisions={} edges={} graph_digest={}",
                MANIFEST_V1_SCHEMA_VERSION,
                self.revisions.len(),
                self.edges.len(),
                digest,
            ),
        )];
        for (ordinal, revision) in self.revisions.iter().enumerate() {
            rows.push((
                "revision".to_owned(),
                format!(
                    "schema_version={} ordinal={} revision={} representative={} graph_digest={}",
                    MANIFEST_V1_SCHEMA_VERSION,
                    ordinal + 1,
                    revision.to_hex(),
                    self.representatives[revision].to_hex(),
                    digest,
                ),
            ));
        }
        for (index, e) in self.edges.iter().enumerate() {
            rows.push((
                format!("edge/{}", relation_name(e.kind)),
                format!(
                    "schema_version={} ordinal={} kind={} from={} to={} checker=\"{}\" tcb=\"{}\" variance={} domain_note=\"{}\" policy_version={} graph_digest={}",
                    MANIFEST_V1_SCHEMA_VERSION,
                    self.revisions.len() + index + 1,
                    relation_name(e.kind),
                    e.from.to_hex(),
                    e.to.to_hex(),
                    escape_json_string(&e.checker),
                    escape_json_string(&e.tcb),
                    variance_name(e.variance),
                    escape_json_string(&e.domain_note),
                    e.policy_version,
                    digest,
                ),
            ));
        }
        rows
    }
}

fn directed(kind: RelationKind) -> bool {
    !matches!(kind, RelationKind::CertifiedEquivalence)
}

fn canonical_endpoints(
    kind: RelationKind,
    mut from: ClaimRevisionId,
    mut to: ClaimRevisionId,
) -> (ClaimRevisionId, ClaimRevisionId) {
    if kind == RelationKind::CertifiedEquivalence && to < from {
        core::mem::swap(&mut from, &mut to);
    }
    (from, to)
}

fn canonicalize_equivalence_endpoints(receipt: &mut ClaimRelationReceipt) {
    (receipt.from, receipt.to) = canonical_endpoints(receipt.kind, receipt.from, receipt.to);
}

fn charge_semantic_bytes(total: &mut u64, value: &str) -> Result<(), V1Error> {
    let bytes = u64::try_from(value.len())
        .map_err(|_| resource_overflow("aggregate semantic bytes", "u64 bytes"))?;
    *total = checked_add_resource(*total, bytes, "aggregate semantic bytes")?;
    Ok(())
}

fn preflight_graph(
    revisions: &[ClaimRevision],
    receipts: &[ClaimRelationReceipt],
    limits: V1AdmissionLimits,
) -> Result<V1ResourceEnvelope, V1Error> {
    limits.validate_protocol_ceiling()?;

    let revision_count_required = u64::try_from(revisions.len())
        .map_err(|_| resource_overflow("revision count", "u64 revisions"))?;
    if revision_count_required > u64::from(limits.max_revisions) {
        return Err(resource_refusal(
            "revision count",
            revision_count_required,
            u64::from(limits.max_revisions),
            "revisions",
        ));
    }
    let receipt_count_required = u64::try_from(receipts.len())
        .map_err(|_| resource_overflow("receipt count", "u64 receipts"))?;
    if receipt_count_required > u64::from(limits.max_receipts) {
        return Err(resource_refusal(
            "receipt count",
            receipt_count_required,
            u64::from(limits.max_receipts),
            "receipts",
        ));
    }
    let revision_count = u32::try_from(revision_count_required)
        .map_err(|_| resource_overflow("revision count", "u32 revisions"))?;
    let receipt_count = u32::try_from(receipt_count_required)
        .map_err(|_| resource_overflow("receipt count", "u32 receipts"))?;

    let projection_rows_required = checked_add_resource(
        1,
        checked_add_resource(
            revision_count_required,
            receipt_count_required,
            "projection row count",
        )?,
        "projection row count",
    )?;
    if projection_rows_required > u64::from(limits.max_projection_rows) {
        return Err(resource_refusal(
            "projection row count",
            projection_rows_required,
            u64::from(limits.max_projection_rows),
            "rows",
        ));
    }
    let projection_rows = u32::try_from(projection_rows_required)
        .map_err(|_| resource_overflow("projection row count", "u32 rows"))?;

    let mut semantic_bytes = 0u64;
    for revision in revisions {
        for value in [
            revision.claim.as_str(),
            &revision.statement,
            &revision.quantifiers,
            &revision.units_conventions,
            &revision.hypotheses,
            &revision.domain,
            &revision.surface,
            &revision.no_claim,
        ] {
            charge_semantic_bytes(&mut semantic_bytes, value)?;
        }
    }

    let mut receipt_semantic_bytes = 0u64;
    for receipt in receipts {
        for value in [
            receipt.checker.as_str(),
            receipt.tcb.as_str(),
            receipt.domain_note.as_str(),
        ] {
            charge_semantic_bytes(&mut semantic_bytes, value)?;
            charge_semantic_bytes(&mut receipt_semantic_bytes, value)?;
        }
    }
    if semantic_bytes > limits.max_semantic_bytes {
        return Err(resource_refusal(
            "aggregate semantic bytes",
            semantic_bytes,
            limits.max_semantic_bytes,
            "UTF-8 bytes",
        ));
    }

    let revision_slots = checked_mul_resource(
        revision_count_required,
        V1_LOGICAL_BYTES_PER_REVISION,
        "logical graph working bytes",
    )?;
    let receipt_slots = checked_mul_resource(
        receipt_count_required,
        V1_LOGICAL_BYTES_PER_RECEIPT,
        "logical graph working bytes",
    )?;
    let logical_working_bytes = checked_add_resource(
        checked_add_resource(revision_slots, receipt_slots, "logical graph working bytes")?,
        receipt_semantic_bytes,
        "logical graph working bytes",
    )?;
    if logical_working_bytes > limits.max_working_bytes {
        return Err(resource_refusal(
            "logical graph working bytes",
            logical_working_bytes,
            limits.max_working_bytes,
            "logical bytes",
        ));
    }

    let projection_bytes_upper_bound = checked_add_resource(
        checked_mul_resource(
            projection_rows_required,
            V1_PROJECTION_FIXED_BYTES_PER_ROW,
            "projection byte upper bound",
        )?,
        checked_mul_resource(
            receipt_semantic_bytes,
            V1_PROJECTION_ESCAPE_EXPANSION,
            "projection byte upper bound",
        )?,
        "projection byte upper bound",
    )?;
    if projection_bytes_upper_bound > limits.max_projection_bytes {
        return Err(resource_refusal(
            "projection byte upper bound",
            projection_bytes_upper_bound,
            limits.max_projection_bytes,
            "logical bytes",
        ));
    }

    // Only after every aggregate/count envelope is proven do we scan the
    // individual semantic contracts. This keeps cap failure precedence stable
    // for hostile batches while still refusing empty or oversized fields.
    for revision in revisions {
        revision.validate()?;
    }
    for receipt in receipts {
        bounded_text("checker", &receipt.checker)?;
        bounded_text("tcb", &receipt.tcb)?;
        bounded_text("domain_note", &receipt.domain_note)?;
    }

    Ok(V1ResourceEnvelope {
        revision_count,
        receipt_count,
        semantic_bytes,
        logical_working_bytes,
        projection_rows,
        projection_bytes_upper_bound,
    })
}

fn try_clone_text(value: &str, payload: &'static str) -> Result<String, V1Error> {
    let required_bytes = u64::try_from(value.len())
        .map_err(|_| resource_overflow("text clone length", "u64 bytes"))?;
    let mut cloned = String::new();
    cloned
        .try_reserve_exact(value.len())
        .map_err(|_| allocation_refusal(payload, required_bytes))?;
    cloned.push_str(value);
    Ok(cloned)
}

fn try_clone_receipts(
    receipts: &[ClaimRelationReceipt],
    resources: V1ResourceEnvelope,
) -> Result<Vec<ClaimRelationReceipt>, V1Error> {
    let mut cloned = Vec::new();
    cloned.try_reserve_exact(receipts.len()).map_err(|_| {
        allocation_refusal("relation receipt index", resources.logical_working_bytes)
    })?;
    for receipt in receipts {
        cloned.push(ClaimRelationReceipt {
            kind: receipt.kind,
            from: receipt.from,
            to: receipt.to,
            checker: try_clone_text(&receipt.checker, "relation checker")?,
            tcb: try_clone_text(&receipt.tcb, "relation TCB")?,
            variance: receipt.variance,
            domain_note: try_clone_text(&receipt.domain_note, "relation domain note")?,
            policy_version: receipt.policy_version,
        });
    }
    Ok(cloned)
}

fn normalized_graph_digest(
    revisions: &[ClaimRevisionId],
    edges: &[ClaimRelationReceipt],
    representatives: &BTreeMap<ClaimRevisionId, ClaimRevisionId>,
    resources: V1ResourceEnvelope,
) -> Result<ContentHash, V1Error> {
    let mut hasher = DomainHasher::new(GRAPH_DOMAIN);
    hasher.update(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
    hasher.update(&resources.revision_count.to_be_bytes());
    for revision in revisions {
        hasher.update(revision.as_bytes());
    }
    hasher.update(&resources.receipt_count.to_be_bytes());
    for edge in edges {
        hasher.update(&[relation_code(edge.kind)]);
        hasher.update(edge.from.as_bytes());
        hasher.update(edge.to.as_bytes());
        hash_framed_text(&mut hasher, "checker", &edge.checker)?;
        hash_framed_text(&mut hasher, "tcb", &edge.tcb)?;
        hasher.update(&[variance_code(edge.variance)]);
        hash_framed_text(&mut hasher, "domain_note", &edge.domain_note)?;
        hasher.update(&edge.policy_version.to_be_bytes());
    }
    for (member, representative) in representatives {
        hasher.update(member.as_bytes());
        hasher.update(representative.as_bytes());
    }
    Ok(hasher.finalize())
}

/// Admit a claim graph: every endpoint must be a known revision and every
/// semantic text field must satisfy its bound. Certified-equivalence endpoints
/// normalize to ascending revision-id order; canonical duplicate receipts
/// refuse. DIRECTED cycles refuse unless every cycle member is joined into one
/// certified-equivalence component, which canonicalizes to its smallest member
/// as representative without erasing anyone.
pub fn admit_graph(
    revisions: &[ClaimRevision],
    receipts: &[ClaimRelationReceipt],
) -> Result<NormalizedGraph, V1Error> {
    admit_graph_with_limits(revisions, receipts, V1AdmissionLimits::DEFAULT)
}

/// Admit a graph under a caller-tightened v1 resource envelope. Limits above
/// [`V1AdmissionLimits::DEFAULT`] refuse; limits are policy metadata and never
/// change the identity of an otherwise identical admitted graph.
pub fn admit_graph_with_limits(
    revisions: &[ClaimRevision],
    receipts: &[ClaimRelationReceipt],
    limits: V1AdmissionLimits,
) -> Result<NormalizedGraph, V1Error> {
    // Cardinality, semantic, logical-working, and projection envelopes are all
    // proven before graph-owned trees, indexes, or receipt text are allocated.
    let resources = preflight_graph(revisions, receipts, limits)?;

    let mut id_list = Vec::new();
    id_list.try_reserve_exact(revisions.len()).map_err(|_| {
        allocation_refusal("revision identity index", resources.logical_working_bytes)
    })?;
    for revision in revisions {
        id_list.push(revision.revision_id_after_validation()?);
    }
    id_list.sort_unstable();
    if id_list.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(V1Error::new(
            "v1-duplicate-revision",
            "identical revision content supplied twice",
        )
        .with_fix("content-addressed identity is idempotent: deduplicate the input"));
    }

    // Validate endpoint semantics against canonical endpoint copies before
    // cloning any receipt text. This preserves bidirectional diagnostics while
    // refusing invalid graphs before graph-owned payload allocation.
    for receipt in receipts {
        let (from, to) = canonical_endpoints(receipt.kind, receipt.from, receipt.to);
        for endpoint in [&from, &to] {
            if id_list.binary_search(endpoint).is_err() {
                return Err(V1Error::new(
                    "v1-dangling-relation",
                    format!(
                        "relation endpoint {} is not an admitted revision",
                        endpoint.to_hex()
                    ),
                )
                .with_fix("admit the endpoint revision first")
                .with_fix("or drop the stale receipt (renames never rebind identities)"));
            }
        }
        if from == to {
            return Err(V1Error::new(
                "v1-self-relation",
                "a revision cannot relate to itself",
            ));
        }
        if receipt.kind == RelationKind::CertifiedEquivalence
            && receipt.variance != QuantifierVariance::Preserved
        {
            return Err(V1Error::new(
                "v1-equivalence-variance",
                "certified equivalence is bidirectional and requires preserved quantifiers",
            )
            .with_fix(
                "use Preserved only with orientation-neutral bidirectional evidence, or choose the directed relation kind matching the variance",
            ));
        }
    }

    // Equivalence is semantically bidirectional: normalize the fallibly cloned
    // private receipts before sorting and all subsequent graph operations.
    // Caller-owned drafts stay untouched.
    let mut edges = try_clone_receipts(receipts, resources)?;
    for receipt in &mut edges {
        canonicalize_equivalence_endpoints(receipt);
    }

    // The canonical receipt order includes EVERY field later fed to the graph
    // digest. Canonical duplicate receipts refuse rather than silently changing
    // cardinality or being presentation-order dependent.
    edges.sort_by(|a, b| {
        relation_code(a.kind)
            .cmp(&relation_code(b.kind))
            .then_with(|| a.from.cmp(&b.from))
            .then_with(|| a.to.cmp(&b.to))
            .then_with(|| a.checker.cmp(&b.checker))
            .then_with(|| a.tcb.cmp(&b.tcb))
            .then_with(|| variance_code(a.variance).cmp(&variance_code(b.variance)))
            .then_with(|| a.domain_note.cmp(&b.domain_note))
            .then_with(|| a.policy_version.cmp(&b.policy_version))
    });
    if edges.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(V1Error::new(
            "v1-duplicate-relation",
            "the same canonical relation receipt was supplied more than once",
        )
        .with_fix(
            "deduplicate the receipt input; reversing certified-equivalence endpoints grants no additional authority",
        ));
    }

    // Union-find over certified-equivalence edges.
    let index: BTreeMap<ClaimRevisionId, usize> =
        id_list.iter().enumerate().map(|(i, id)| (*id, i)).collect();
    let mut parent = Vec::new();
    parent.try_reserve_exact(id_list.len()).map_err(|_| {
        allocation_refusal("equivalence union-find", resources.logical_working_bytes)
    })?;
    parent.extend(0..id_list.len());
    // Iterative path-halving find (explicit-stack doctrine: no recursion).
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for receipt in &edges {
        if receipt.kind == RelationKind::CertifiedEquivalence {
            let a = find(&mut parent, index[&receipt.from]);
            let b = find(&mut parent, index[&receipt.to]);
            parent[a.max(b)] = a.min(b);
        }
    }
    let mut representative = BTreeMap::new();
    for (i, id) in id_list.iter().enumerate() {
        let root = find(&mut parent, i);
        representative.insert(*id, id_list[root]);
    }

    // Directed-cycle check on the QUOTIENT graph (equivalence components
    // collapsed to representatives): any remaining directed cycle refuses.
    let mut adjacency: BTreeMap<ClaimRevisionId, BTreeSet<ClaimRevisionId>> = BTreeMap::new();
    for receipt in &edges {
        if directed(receipt.kind) {
            let from = representative[&receipt.from];
            let to = representative[&receipt.to];
            if from != to {
                adjacency.entry(from).or_default().insert(to);
            } else if receipt.kind != RelationKind::Counterexample {
                // A directed promotion-bearing edge INSIDE an equivalence
                // component is redundant but sound; a counterexample
                // inside one is a contradiction.
            } else {
                return Err(V1Error::new(
                    "v1-domain-contradiction",
                    "a counterexample edge joins certified-equivalent revisions",
                )
                .with_fix("revoke the equivalence certification or the counterexample"));
            }
        }
    }
    // Iterative DFS three-color cycle detection (deterministic order).
    let mut color: BTreeMap<ClaimRevisionId, u8> = BTreeMap::new();
    let stack_slots = id_list
        .len()
        .checked_add(edges.len())
        .and_then(|slots| slots.checked_add(1))
        .ok_or_else(|| resource_overflow("cycle stack slots", "usize slots"))?;
    let mut stack = Vec::new();
    stack
        .try_reserve_exact(stack_slots)
        .map_err(|_| allocation_refusal("cycle detector stack", resources.logical_working_bytes))?;
    for &start in adjacency.keys() {
        if color.get(&start).copied().unwrap_or(0) != 0 {
            continue;
        }
        stack.clear();
        stack.push((start, false));
        while let Some((node, children_done)) = stack.pop() {
            if children_done {
                color.insert(node, 2);
                continue;
            }
            match color.get(&node).copied().unwrap_or(0) {
                1 | 2 => continue,
                _ => {}
            }
            color.insert(node, 1);
            stack.push((node, true));
            if let Some(next) = adjacency.get(&node) {
                for &n in next {
                    match color.get(&n).copied().unwrap_or(0) {
                        1 => {
                            return Err(V1Error::new(
                                "v1-relation-cycle",
                                format!(
                                    "directed relation cycle through {} — cycles are admitted \
                                     only as certified equivalence",
                                    n.to_hex()
                                ),
                            )
                            .with_fix("certify the equivalence (both directions, checker + TCB)")
                            .with_fix("or break the cycle by revoking one implication"));
                        }
                        0 => stack.push((n, false)),
                        _ => {}
                    }
                }
            }
        }
    }

    let digest = normalized_graph_digest(&id_list, &edges, &representative, resources)?;
    Ok(NormalizedGraph {
        revisions: id_list,
        edges,
        representatives: representative,
        digest,
        resources,
    })
}

/// Migration classification between two record field sets: additive
/// (old readers still sound) versus breaking (a lossy report is
/// REQUIRED, naming every dropped field).
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum Migration {
    /// New fields only; old readers remain sound.
    Additive {
        /// The added field names.
        added: Vec<String>,
    },
    /// Fields were removed or re-typed: breaking, with the explicit
    /// lossy report.
    Breaking {
        /// The lossy-migration report: every dropped field and why.
        report: LossyMigrationReport,
    },
}

/// The explicit record of what a breaking migration loses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LossyMigrationReport {
    /// Fields present before and absent after.
    pub dropped: Vec<String>,
    /// The stated reason.
    pub reason: String,
}

/// Classify a migration from `old` to `new` field sets. Removals demand
/// a reason; without one the migration refuses rather than silently
/// dropping semantics.
pub fn classify_migration(
    old: &BTreeSet<String>,
    new: &BTreeSet<String>,
    lossy_reason: Option<&str>,
) -> Result<Migration, V1Error> {
    let dropped: Vec<String> = old.difference(new).cloned().collect();
    let added: Vec<String> = new.difference(old).cloned().collect();
    if dropped.is_empty() {
        return Ok(Migration::Additive { added });
    }
    match lossy_reason {
        Some(reason) => Ok(Migration::Breaking {
            report: LossyMigrationReport {
                dropped,
                reason: reason.to_owned(),
            },
        }),
        None => Err(V1Error::new(
            "v1-lossy-migration-undeclared",
            format!("migration drops {dropped:?} without a lossy report"),
        )
        .with_fix("declare the lossy migration report with a reason")
        .with_fix("or keep the fields and deprecate additively")),
    }
}

/// Per-field metadata for the manifest record: units, cardinality,
/// authority, default visibility, and migration semantics — declared in
/// data so the success criterion is checkable, not aspirational.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldSpec {
    /// The field name.
    pub name: &'static str,
    /// Units or "identity"/"text"/"hash" for unitless fields.
    pub units: &'static str,
    /// Cardinality ("1", "0..1", "0..n", "1..n").
    pub cardinality: &'static str,
    /// The authority that owns the field's truth.
    pub authority: SourceAuthority,
    /// Whether the field is visible by default in human projections.
    pub default_visible: bool,
    /// Migration semantics ("additive", "breaking-if-removed",
    /// "identity-forming").
    pub migration: &'static str,
}

/// The declared v1 manifest-record field registry. Identity-forming
/// fields participate in digests; removing any field here is a breaking
/// migration by definition.
pub const MANIFEST_RECORD_FIELDS: &[FieldSpec] = &[
    FieldSpec {
        name: "source-snapshots",
        units: "hash",
        cardinality: "1..n",
        authority: SourceAuthority::FrozenSnapshot,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "bead-obligation",
        units: "identity",
        cardinality: "1",
        authority: SourceAuthority::BeadObligation,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "claim-revision",
        units: "hash",
        cardinality: "1",
        authority: SourceAuthority::FrozenSnapshot,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "stratum",
        units: "enum(core|max)",
        cardinality: "1",
        authority: SourceAuthority::BeadObligation,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "campaign-profiles",
        units: "identity",
        cardinality: "1..n",
        authority: SourceAuthority::BeadObligation,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "ambition",
        units: "enum(S|F|M)",
        cardinality: "1",
        authority: SourceAuthority::BeadObligation,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "public-surface",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::Contract,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "case-ids",
        units: "identity",
        cardinality: "1..n",
        authority: SourceAuthority::TestSource,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "journey-ids",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::TestSource,
        default_visible: true,
        migration: "additive",
    },
    FieldSpec {
        name: "ownership",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::BeadObligation,
        default_visible: true,
        migration: "additive",
    },
    FieldSpec {
        name: "fixture-ids",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::TestSource,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "oracle-ids",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::TestSource,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "checker-ids",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::TestSource,
        default_visible: true,
        migration: "identity-forming",
    },
    FieldSpec {
        name: "tcb-overlap",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::Contract,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "tolerance-derivation",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::Contract,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "budgets",
        units: "mixed(bytes|ns|units)",
        cardinality: "0..n",
        authority: SourceAuthority::BeadObligation,
        default_visible: false,
        migration: "additive",
    },
    FieldSpec {
        name: "capabilities",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::BeadObligation,
        default_visible: false,
        migration: "additive",
    },
    FieldSpec {
        name: "event-kinds",
        units: "identity",
        cardinality: "0..n",
        authority: SourceAuthority::Contract,
        default_visible: false,
        migration: "additive",
    },
    FieldSpec {
        name: "retention",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::BeadObligation,
        default_visible: false,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "replay-command",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::GeneratedArtifact,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "dsr-lane",
        units: "text",
        cardinality: "1",
        authority: SourceAuthority::GeneratedArtifact,
        default_visible: true,
        migration: "breaking-if-removed",
    },
    FieldSpec {
        name: "receipt-expectations",
        units: "text",
        cardinality: "1..n",
        authority: SourceAuthority::Contract,
        default_visible: true,
        migration: "breaking-if-removed",
    },
];

/// Resolve a source-authority conflict: two pins for the same source at
/// different snapshots. The higher authority wins ONLY by minting the
/// record against its snapshot explicitly; equal-authority conflicts
/// refuse with ranked fixes (never a silent pick).
pub fn resolve_authority<'a>(a: &'a SourcePin, b: &'a SourcePin) -> Result<&'a SourcePin, V1Error> {
    if a.source != b.source {
        return Err(V1Error::new(
            "v1-authority-mismatch",
            "conflict resolution requires pins of the same source",
        ));
    }
    if a.snapshot == b.snapshot {
        return Ok(a);
    }
    match a.authority.cmp(&b.authority) {
        core::cmp::Ordering::Greater => Ok(a),
        core::cmp::Ordering::Less => Ok(b),
        core::cmp::Ordering::Equal => Err(V1Error::new(
            "v1-authority-conflict",
            format!(
                "source {:?} pinned at two snapshots with equal authority {:?}",
                a.source, a.authority
            ),
        )
        .with_fix("re-freeze the source and re-pin both records against the new snapshot")
        .with_fix("or split the records so each names its own source")),
    }
}
