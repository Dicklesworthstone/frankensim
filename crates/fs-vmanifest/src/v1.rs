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
//! equivalence) with direction, checker/TCB, quantifier variance, and
//! policy version. Promotion transfers ONLY across admitted relations;
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

use fs_blake3::{ContentHash, hash_domain};

/// Version stamped into every v1 record and digest domain.
pub const MANIFEST_V1_SCHEMA_VERSION: u32 = 1;

const REVISION_DOMAIN: &str = "org.frankensim.fs-vmanifest.claim-revision.v1";
const GRAPH_DOMAIN: &str = "org.frankensim.fs-vmanifest.manifest-graph.v1";
const MAX_V1_TEXT_BYTES: usize = 4096;
const MAX_V1_ID_BYTES: usize = 128;

/// A typed refusal with a stable rule slug and RANKED candidate fixes —
/// the "stable diagnostics and ranked fixes" the success criteria demand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct V1Error {
    rule: &'static str,
    detail: String,
    fixes: Vec<String>,
}

impl V1Error {
    fn new(rule: &'static str, detail: impl Into<String>) -> V1Error {
        V1Error {
            rule,
            detail: detail.into(),
            fixes: Vec::new(),
        }
    }

    fn with_fix(mut self, fix: impl Into<String>) -> V1Error {
        self.fixes.push(fix.into());
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
            '\u{0000}'..='\u{001f}' => {
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

fn push_field(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u32).to_be_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

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

    /// The content-addressed revision identity: a domain-separated hash
    /// over every semantic field in fixed order. Distinct content cannot
    /// collide; identical content is idempotent.
    #[must_use]
    pub fn revision_id(&self) -> ClaimRevisionId {
        let mut bytes = Vec::with_capacity(256);
        bytes.extend_from_slice(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
        push_field(&mut bytes, self.claim.as_str());
        bytes.push(claim_kind_code(self.kind));
        push_field(&mut bytes, &self.statement);
        push_field(&mut bytes, &self.quantifiers);
        push_field(&mut bytes, &self.units_conventions);
        push_field(&mut bytes, &self.hypotheses);
        push_field(&mut bytes, &self.domain);
        push_field(&mut bytes, &self.surface);
        push_field(&mut bytes, &self.no_claim);
        match &self.supersedes {
            None => bytes.push(0),
            Some(prev) => {
                bytes.push(1);
                bytes.extend_from_slice(prev.as_bytes());
            }
        }
        hash_domain(REVISION_DOMAIN, &bytes)
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

/// One typed, checkable relation edge between exact claim revisions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClaimRelationReceipt {
    /// The edge kind.
    pub kind: RelationKind,
    /// Source revision.
    pub from: ClaimRevisionId,
    /// Target revision.
    pub to: ClaimRevisionId,
    /// The proof/checker identity that certified this edge.
    pub checker: String,
    /// The trusted computing base the checker ran under.
    pub tcb: String,
    /// Quantifier variance across the edge.
    pub variance: QuantifierVariance,
    /// Domain relationship note (subdomain, identical, disjoint —
    /// free text bound into the receipt identity).
    pub domain_note: String,
    /// The relation-policy version this edge was admitted under.
    pub policy_version: u32,
}

impl ClaimRelationReceipt {
    /// Whether PROMOTION may transfer along this edge (fallback/refutation
    /// rules are the mirror image and live with the consuming gate).
    /// Counterexample edges never transfer promotion; strengthened
    /// quantifiers never transfer promotion.
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
    pub revisions: Vec<ClaimRevisionId>,
    /// Every admitted edge, canonically sorted.
    pub edges: Vec<ClaimRelationReceipt>,
    /// Certified-equivalence component representative per member
    /// (identity for revisions outside any component).
    pub representative: BTreeMap<ClaimRevisionId, ClaimRevisionId>,
}

impl NormalizedGraph {
    /// The domain-separated digest of the normalized graph.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
        bytes.extend_from_slice(&(self.revisions.len() as u32).to_be_bytes());
        for r in &self.revisions {
            bytes.extend_from_slice(r.as_bytes());
        }
        bytes.extend_from_slice(&(self.edges.len() as u32).to_be_bytes());
        for e in &self.edges {
            bytes.push(relation_code(e.kind));
            bytes.extend_from_slice(e.from.as_bytes());
            bytes.extend_from_slice(e.to.as_bytes());
            push_field(&mut bytes, &e.checker);
            push_field(&mut bytes, &e.tcb);
            bytes.push(match e.variance {
                QuantifierVariance::Preserved => 1,
                QuantifierVariance::Weakened => 2,
                QuantifierVariance::Strengthened => 3,
            });
            push_field(&mut bytes, &e.domain_note);
            bytes.extend_from_slice(&e.policy_version.to_be_bytes());
        }
        for (member, repr) in &self.representative {
            bytes.extend_from_slice(member.as_bytes());
            bytes.extend_from_slice(repr.as_bytes());
        }
        hash_domain(GRAPH_DOMAIN, &bytes)
    }

    /// Human projection: one line per revision and edge, digest-stamped.
    #[must_use]
    pub fn render_human(&self) -> String {
        let mut out = format!(
            "VerificationManifest v1 graph — {} revisions, {} edges, digest {}\n",
            self.revisions.len(),
            self.edges.len(),
            self.digest().to_hex()
        );
        for e in &self.edges {
            out.push_str(&format!(
                "  {:?}: {} -> {} [{}]\n",
                e.kind,
                &e.from.to_hex()[..12],
                &e.to.to_hex()[..12],
                e.checker
            ));
        }
        out
    }

    /// JSON-lines projection, digest-stamped.
    #[must_use]
    pub fn render_json_lines(&self) -> Vec<String> {
        let mut out = vec![format!(
            "{{\"vmanifest\":\"graph\",\"revisions\":{},\"edges\":{},\"digest\":\"{}\"}}",
            self.revisions.len(),
            self.edges.len(),
            self.digest().to_hex()
        )];
        for e in &self.edges {
            out.push(format!(
                "{{\"vmanifest\":\"edge\",\"kind\":\"{:?}\",\"from\":\"{}\",\"to\":\"{}\",\"checker\":\"{}\",\"policy\":{}}}",
                e.kind,
                e.from.to_hex(),
                e.to.to_hex(),
                escape_json_string(&e.checker),
                e.policy_version
            ));
        }
        out
    }

    /// Ledger-row projection: (row kind, digest-bound payload) tuples.
    #[must_use]
    pub fn render_ledger_rows(&self) -> Vec<(String, String)> {
        let digest = self.digest().to_hex();
        let mut rows = vec![("graph-digest".to_owned(), digest.clone())];
        for e in &self.edges {
            rows.push((
                format!("edge/{:?}", e.kind),
                format!("{}->{}@{}", e.from.to_hex(), e.to.to_hex(), digest),
            ));
        }
        rows
    }
}

fn directed(kind: RelationKind) -> bool {
    !matches!(kind, RelationKind::CertifiedEquivalence)
}

/// Admit a claim graph: every endpoint must be a known revision; DIRECTED
/// cycles refuse unless every cycle member is joined into one certified-
/// equivalence component, which canonicalizes to its smallest member as
/// representative without erasing anyone.
pub fn admit_graph(
    revisions: &[ClaimRevision],
    receipts: &[ClaimRelationReceipt],
) -> Result<NormalizedGraph, V1Error> {
    let mut ids = BTreeSet::new();
    for revision in revisions {
        revision.validate()?;
        if !ids.insert(revision.revision_id()) {
            return Err(V1Error::new(
                "v1-duplicate-revision",
                "identical revision content supplied twice",
            )
            .with_fix("content-addressed identity is idempotent: deduplicate the input"));
        }
    }
    for receipt in receipts {
        bounded_text("checker", &receipt.checker)?;
        bounded_text("tcb", &receipt.tcb)?;
        for endpoint in [&receipt.from, &receipt.to] {
            if !ids.contains(endpoint) {
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
        if receipt.from == receipt.to {
            return Err(V1Error::new(
                "v1-self-relation",
                "a revision cannot relate to itself",
            ));
        }
    }

    // Union-find over certified-equivalence edges.
    let id_list: Vec<ClaimRevisionId> = ids.iter().copied().collect();
    let index: BTreeMap<ClaimRevisionId, usize> =
        id_list.iter().enumerate().map(|(i, id)| (*id, i)).collect();
    let mut parent: Vec<usize> = (0..id_list.len()).collect();
    // Iterative path-halving find (explicit-stack doctrine: no recursion).
    fn find(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for receipt in receipts {
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
    for receipt in receipts {
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
    for &start in adjacency.keys() {
        if color.get(&start).copied().unwrap_or(0) != 0 {
            continue;
        }
        let mut stack = vec![(start, false)];
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

    let mut edges = receipts.to_vec();
    edges.sort_by(|a, b| {
        (relation_code(a.kind), a.from, a.to, &a.checker).cmp(&(
            relation_code(b.kind),
            b.from,
            b.to,
            &b.checker,
        ))
    });
    Ok(NormalizedGraph {
        revisions: id_list,
        edges,
        representative,
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
