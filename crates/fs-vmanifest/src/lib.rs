//! fs-vmanifest — the typed VerificationManifest schema for leapfrog G1
//! claim/evidence freezes (first instance: bead
//! frankensim-leapfrog-2026-program-i94v.1.1.8.1, initiative I01). Layer:
//! UTIL.
//!
//! A G1 freeze preregisters a campaign's claims, fixtures, obligations,
//! and waivers BEFORE any implementation result is inspected, so
//! tolerances, claim wording, capability scope, and failure policy cannot
//! drift toward favorable outcomes. The schema enforces that discipline
//! structurally:
//!
//! - Freezing is fail-closed: [`ManifestDraft::freeze`] refuses blank
//!   fields, over-cap collections (checked BEFORE any deep scan),
//!   duplicate ids, non-independent oracles, orphan deck references,
//!   uncovered claims, invalid tolerances, and malformed digests — each
//!   with a typed [`FreezeRefusal`] naming the gate.
//! - Authority is sealed: [`FrozenManifest`] has no public constructor
//!   and no mutating API; holding one proves the gates ran. Post-freeze
//!   "alteration" is impossible by construction — change happens only
//!   through [`FrozenManifest::amend`], which requires the successor
//!   version and names exactly the invalidated descendants.
//! - Identity is canonical: components sort into one total order with
//!   content tie-breaks, and the manifest digest is a domain-separated,
//!   length-framed BLAKE3 hash, byte-stable across runs on the same ISA.
//!
//! Preregistration is not proof: a frozen manifest asserts nothing about
//! implementation correctness, and no color or promotion authority is
//! minted here.

pub use fs_blake3::ContentHash;

use fs_blake3::hash_domain;
use std::fmt;

mod i01;
mod i02;
mod i04;
mod i08;
mod i12;
mod i15;

pub use i01::i01_draft;
pub use i02::i02_draft;
pub use i04::i04_draft;
pub use i08::i08_draft;
pub use i12::i12_draft;
pub use i15::i15_draft;

/// Manifest schema version (canonical bytes are comparable only within it).
pub const VMANIFEST_SCHEMA_VERSION: u32 = 1;
/// Maximum claims per manifest (checked before any per-claim scan).
pub const MAX_CLAIMS: usize = 256;
/// Maximum fixture pins per manifest.
pub const MAX_FIXTURES: usize = 512;
/// Maximum obligation rows per manifest.
pub const MAX_OBLIGATIONS: usize = 256;
/// Maximum waivers per manifest.
pub const MAX_WAIVERS: usize = 128;
/// Maximum items in any per-row string list (hypotheses, decks, events...).
pub const MAX_ROW_ITEMS: usize = 64;

const MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.manifest.v1";
const CLAIM_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.claim.v1";
const FIXTURE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.fixture.v1";
const OBLIGATION_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.obligation.v1";
const WAIVER_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.waiver.v1";
const SPEC_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.fixture-spec.v1";

/// Ambition lattice element ([S]/[F]/[M]); a weaker receipt closes its own
/// element and is never relabeled as the stronger theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ambition {
    /// Established mathematics/engineering.
    Solid,
    /// Research-backed, high-upside engineering risk.
    Frontier,
    /// Novel synthesis; stays behind flags until proven.
    Moonshot,
}

impl Ambition {
    /// Stable single-letter tag.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            Self::Solid => "S",
            Self::Frontier => "F",
            Self::Moonshot => "M",
        }
    }

    const fn byte(self) -> u8 {
        match self {
            Self::Solid => 1,
            Self::Frontier => 2,
            Self::Moonshot => 3,
        }
    }
}

/// Gauntlet evidence tier a claim targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GauntletTier {
    /// Property tests and algebraic laws.
    G0,
    /// Manufactured solutions and convergence-order verification.
    G1,
    /// Canonical benchmarks.
    G2,
    /// Metamorphic tests.
    G3,
    /// Chaos, cancellation storms, leak/deadlock checks.
    G4,
    /// Determinism audits.
    G5,
}

impl GauntletTier {
    /// Stable tag.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            Self::G0 => "G0",
            Self::G1 => "G1",
            Self::G2 => "G2",
            Self::G3 => "G3",
            Self::G4 => "G4",
            Self::G5 => "G5",
        }
    }

    const fn byte(self) -> u8 {
        match self {
            Self::G0 => 0,
            Self::G1 => 1,
            Self::G2 => 2,
            Self::G3 => 3,
            Self::G4 => 4,
            Self::G5 => 5,
        }
    }
}

/// Whether the claim asserts a property or preregisters a refutation route.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClaimPolarity {
    /// The campaign tries to establish the stated property.
    Affirmative,
    /// The campaign tries to refute the stated property (falsifier lane).
    Refutation,
}

impl ClaimPolarity {
    const fn byte(self) -> u8 {
        match self {
            Self::Affirmative => 1,
            Self::Refutation => 2,
        }
    }
}

/// Acceptance arithmetic for a claim's QoI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ToleranceSemantics {
    /// `|candidate - reference| <= atol` (finite, `> 0`).
    Absolute {
        /// Absolute tolerance.
        atol: f64,
    },
    /// `|candidate - reference| <= rtol * |reference|` (finite, `> 0`).
    Relative {
        /// Relative tolerance.
        rtol: f64,
    },
    /// Combined absolute + relative acceptance (each finite, `>= 0`, not
    /// both zero).
    AbsRel {
        /// Absolute part.
        atol: f64,
        /// Relative part.
        rtol: f64,
    },
    /// Candidate must lie in `[lo, hi]` (finite, `lo <= hi`).
    Interval {
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
    /// Exact boolean/bitwise verdict (replay equality, admission verdicts).
    Exact,
}

/// The independent checker route for a claim. Reusing the production code
/// path as its own oracle is a freeze refusal, not a waivable style issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OracleRoute {
    /// Stable oracle identity (crate/path/method).
    pub identity: &'static str,
    /// Declared independence from the production path under test.
    pub independent: bool,
    /// Declared trusted-computing-base overlap with the production path
    /// (honesty about shared kernels; `"none"` if disjoint).
    pub tcb_overlap: &'static str,
}

/// One preregistered claim: hypotheses, QoI, acceptance, oracle,
/// activation/kill criteria, fallback, and the Unknown/no-claim boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClaimSpec {
    /// Stable claim id (kebab-case, initiative-prefixed).
    pub id: &'static str,
    /// Lattice element this claim closes.
    pub ambition: Ambition,
    /// Affirmative or refutation lane.
    pub polarity: ClaimPolarity,
    /// The claim statement.
    pub statement: &'static str,
    /// Explicit hypotheses (must be non-empty).
    pub hypotheses: &'static [&'static str],
    /// Quantity of interest.
    pub qoi: &'static str,
    /// QoI unit (`"1"` for dimensionless, `"bit"` for exact verdicts).
    pub unit: &'static str,
    /// Acceptance arithmetic.
    pub tolerance: ToleranceSemantics,
    /// Targeted Gauntlet evidence tier.
    pub evidence_tier: GauntletTier,
    /// Independent checker route.
    pub oracle: OracleRoute,
    /// Activation criterion (when the claim's campaign lane starts).
    pub activation: &'static str,
    /// Kill criterion (when the lane is abandoned, with its receipt).
    pub kill: &'static str,
    /// Fallback when the claim fails or stays Unknown.
    pub fallback: &'static str,
    /// The explicit Unknown/no-claim boundary.
    pub no_claim: &'static str,
}

/// Which partition a fixture belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Partition {
    /// Visible during development.
    Development,
    /// Withheld until adjudication (held-out).
    HeldOut,
}

impl Partition {
    const fn byte(self) -> u8 {
        match self {
            Self::Development => 1,
            Self::HeldOut => 2,
        }
    }
}

/// How a fixture's bytes are identified.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixtureSource {
    /// A canonical generator/spec text authored in this crate; the digest
    /// is computed from these exact bytes.
    AuthoredSpec {
        /// The complete canonical spec text.
        spec: &'static str,
    },
    /// External artifact pinned by its BLAKE3 digest (64 hex chars).
    External {
        /// Lowercase 64-char hex of the artifact digest.
        digest_hex: &'static str,
    },
}

/// One pinned fixture corpus element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FixturePin {
    /// Stable fixture id (referenced by obligation deck lists).
    pub id: &'static str,
    /// Byte identity.
    pub source: FixtureSource,
    /// Development or held-out partition.
    pub partition: Partition,
}

impl FixturePin {
    /// The fixture's content identity, if well-formed.
    #[must_use]
    pub fn digest(&self) -> Option<ContentHash> {
        match self.source {
            FixtureSource::AuthoredSpec { spec } => {
                if spec.trim().is_empty() {
                    None
                } else {
                    Some(hash_domain(SPEC_IDENTITY_DOMAIN, spec.as_bytes()))
                }
            }
            FixtureSource::External { digest_hex } => ContentHash::from_hex(digest_hex),
        }
    }
}

/// Campaign tier an obligation row runs at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CampaignTier {
    /// Cheap always-on smoke slice.
    Smoke,
    /// The core evidence run.
    Core,
    /// The maximal fidelity run.
    Max,
}

impl CampaignTier {
    const fn byte(self) -> u8 {
        match self {
            Self::Smoke => 1,
            Self::Core => 2,
            Self::Max => 3,
        }
    }
}

/// One execution leaf's complete verification obligation: the mapping the
/// bead calls "no unnamed skips".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObligationRow {
    /// Execution-leaf/cluster id.
    pub leaf: &'static str,
    /// Claim ids this row's evidence feeds (must exist in the manifest).
    pub claims_covered: &'static [&'static str],
    /// Required unit-case classes (happy/empty/boundary/max/error/unit/
    /// tie/cancellation/migration).
    pub unit_cases: &'static [&'static str],
    /// G0 generators, validity predicates, laws, shrinkers, replay seeds.
    pub g0: &'static str,
    /// G1/G2 deck ids (each must resolve to a fixture pin or a waiver).
    pub decks: &'static [&'static str],
    /// G3 metamorphic relations.
    pub g3_relations: &'static [&'static str],
    /// G4 fault/cancellation/request-drain-finalize/checkpoint schedule.
    pub g4_schedule: &'static str,
    /// G5 thread/shard/mode/ISA determinism matrix.
    pub g5_matrix: &'static str,
    /// The named campaign entry point (`scripts/e2e/leapfrog/*.sh`).
    pub entry_point: &'static str,
    /// Smoke/core/max tier.
    pub tier: CampaignTier,
    /// DSR lane that owns the run.
    pub dsr_lane: &'static str,
    /// Required fs-obs event kinds.
    pub obs_events: &'static [&'static str],
    /// Exact replay command.
    pub replay_command: &'static str,
}

/// A named skip: narrow reason, owner, predicate, expiry, and its explicit
/// effect on promotion. There are no unnamed skips.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Waiver {
    /// The claim/leaf/deck id this waiver covers.
    pub subject: &'static str,
    /// Narrow reason for the skip.
    pub reason: &'static str,
    /// Owner responsible for discharging it.
    pub owner: &'static str,
    /// Predicate whose truth retires the waiver.
    pub predicate: &'static str,
    /// Expiry/review point.
    pub expiry: &'static str,
    /// Explicit effect on promotion while the waiver is live.
    pub promotion_effect: &'static str,
}

/// The Five Explicits carried by every manifest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FiveExplicits {
    /// Unit system declaration.
    pub units: &'static str,
    /// Seed policy (streams, ranges, partition split).
    pub seeds: &'static str,
    /// Budget declaration (time/memory/accuracy).
    pub budgets: &'static str,
    /// Version pins (schema, toolchain, constellation).
    pub versions: &'static str,
    /// Capability flags in force.
    pub capabilities: &'static str,
}

/// A mutable manifest under assembly. Freezing consumes it.
#[derive(Debug, Clone, PartialEq)]
pub struct ManifestDraft {
    /// Initiative id (e.g. `"I01"`).
    pub initiative: &'static str,
    /// Human-readable campaign title.
    pub title: &'static str,
    /// Manifest version (`>= 1`; amendments increment it).
    pub version: u32,
    /// The Five Explicits.
    pub explicits: FiveExplicits,
    /// Preregistered claims.
    pub claims: Vec<ClaimSpec>,
    /// Pinned fixture corpus.
    pub fixtures: Vec<FixturePin>,
    /// Execution-leaf obligations.
    pub obligations: Vec<ObligationRow>,
    /// Named skips.
    pub waivers: Vec<Waiver>,
    /// The rules under which this manifest may be amended.
    pub amendment_rules: &'static str,
}

/// Why a draft cannot freeze. Variants name the first failing gate in the
/// documented order: caps, blank fields, versions, duplicates, oracle
/// independence, tolerance validity, fixture well-formedness, coverage,
/// orphan references.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FreezeRefusal {
    /// A collection exceeds its cap (checked before any deep scan).
    OverCap {
        /// Which collection.
        what: &'static str,
        /// Offered length.
        len: usize,
        /// The cap.
        cap: usize,
    },
    /// Version must be `>= 1`.
    ZeroVersion,
    /// A required text field is blank.
    BlankField {
        /// Owning component id (manifest initiative when top-level).
        id: String,
        /// Which field.
        field: &'static str,
    },
    /// Two components share an id.
    DuplicateId {
        /// Component kind.
        kind: &'static str,
        /// The colliding id.
        id: String,
    },
    /// A claim's oracle is not independent of the production path.
    ProductionOracleReuse {
        /// The claim id.
        claim: String,
    },
    /// A claim's tolerance arithmetic is invalid.
    InvalidTolerance {
        /// The claim id.
        claim: String,
        /// What is wrong.
        reason: &'static str,
    },
    /// A fixture pin is malformed (blank spec or bad digest hex).
    MalformedFixture {
        /// The fixture id.
        fixture: String,
    },
    /// An obligation references a deck id with no fixture pin and no
    /// waiver (the no-orphan lint).
    OrphanDeck {
        /// The obligation leaf.
        leaf: String,
        /// The unresolved deck id.
        deck: String,
    },
    /// An obligation covers a claim id that does not exist.
    OrphanClaimRef {
        /// The obligation leaf.
        leaf: String,
        /// The unknown claim id.
        claim: String,
    },
    /// A claim is covered by no obligation row and no waiver.
    UncoveredClaim {
        /// The uncovered claim id.
        claim: String,
    },
}

impl fmt::Display for FreezeRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OverCap { what, len, cap } => {
                write!(f, "{what}: {len} exceeds the cap of {cap}")
            }
            Self::ZeroVersion => f.write_str("manifest version must be >= 1"),
            Self::BlankField { id, field } => {
                write!(f, "'{id}': required field '{field}' is blank")
            }
            Self::DuplicateId { kind, id } => {
                write!(f, "duplicate {kind} id '{id}'; freeze fails closed")
            }
            Self::ProductionOracleReuse { claim } => write!(
                f,
                "claim '{claim}': oracle is not independent of the production path"
            ),
            Self::InvalidTolerance { claim, reason } => {
                write!(f, "claim '{claim}': tolerance invalid: {reason}")
            }
            Self::MalformedFixture { fixture } => {
                write!(f, "fixture '{fixture}': blank spec or malformed digest")
            }
            Self::OrphanDeck { leaf, deck } => write!(
                f,
                "obligation '{leaf}': deck '{deck}' has no fixture pin and no waiver"
            ),
            Self::OrphanClaimRef { leaf, claim } => write!(
                f,
                "obligation '{leaf}': covered claim '{claim}' does not exist"
            ),
            Self::UncoveredClaim { claim } => write!(
                f,
                "claim '{claim}' is covered by no obligation row and no waiver"
            ),
        }
    }
}

impl std::error::Error for FreezeRefusal {}

/// Why an amendment was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AmendmentRefusal {
    /// The successor draft must carry `version == predecessor + 1`.
    WrongVersion {
        /// Expected successor version.
        expected: u32,
        /// Offered version.
        offered: u32,
    },
    /// The successor draft itself failed its freeze gates.
    SuccessorRefused(FreezeRefusal),
}

impl fmt::Display for AmendmentRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongVersion { expected, offered } => {
                write!(f, "amendment must carry version {expected}, got {offered}")
            }
            Self::SuccessorRefused(refusal) => {
                write!(f, "successor draft refused: {refusal}")
            }
        }
    }
}

impl std::error::Error for AmendmentRefusal {}

/// The record an amendment produces: which version replaced which, and
/// exactly which evidence descendants are invalidated (claims whose
/// content identity changed or vanished, plus obligations likewise).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AmendmentRecord {
    /// Predecessor version.
    pub from_version: u32,
    /// Successor version.
    pub to_version: u32,
    /// Predecessor manifest digest.
    pub from_digest: ContentHash,
    /// Successor manifest digest.
    pub to_digest: ContentHash,
    /// Ids (claims and obligation leaves) whose evidence is invalidated.
    pub invalidated: Vec<String>,
}

/// A sealed, immutable, frozen manifest.
///
/// SEALED: no public constructor and no mutating API. The only producers
/// are [`ManifestDraft::freeze`] and [`FrozenManifest::amend`], so holding
/// one proves the fail-closed gates ran on exactly this content.
#[derive(Debug, Clone, PartialEq)]
pub struct FrozenManifest {
    draft: ManifestDraft,
    digest: ContentHash,
}

impl ManifestDraft {
    /// Freeze the draft, consuming it.
    ///
    /// Gate order (documented so refusal tests are stable): collection
    /// caps (before any deep scan), version, top-level blank fields,
    /// per-component blank fields and list caps, duplicate ids, oracle
    /// independence, tolerance validity, fixture well-formedness, orphan
    /// claim references, orphan decks, uncovered claims.
    ///
    /// # Errors
    ///
    /// The typed [`FreezeRefusal`] for the first failing gate.
    pub fn freeze(self) -> Result<FrozenManifest, FreezeRefusal> {
        check_caps(&self)?;
        if self.version == 0 {
            return Err(FreezeRefusal::ZeroVersion);
        }
        check_top_level_text(&self)?;
        check_components(&self)?;
        check_duplicates(&self)?;
        check_claims(&self)?;
        check_fixtures(&self)?;
        check_references(&self)?;
        let digest = manifest_digest(&self);
        Ok(FrozenManifest {
            draft: canonicalize(self),
            digest,
        })
    }
}

impl FrozenManifest {
    /// Initiative id.
    #[must_use]
    pub fn initiative(&self) -> &'static str {
        self.draft.initiative
    }

    /// Manifest version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.draft.version
    }

    /// The canonical content identity this freeze is bound to.
    #[must_use]
    pub const fn digest(&self) -> ContentHash {
        self.digest
    }

    /// The frozen claims, in canonical order.
    #[must_use]
    pub fn claims(&self) -> &[ClaimSpec] {
        &self.draft.claims
    }

    /// The frozen fixture pins, in canonical order.
    #[must_use]
    pub fn fixtures(&self) -> &[FixturePin] {
        &self.draft.fixtures
    }

    /// The frozen obligations, in canonical order.
    #[must_use]
    pub fn obligations(&self) -> &[ObligationRow] {
        &self.draft.obligations
    }

    /// The frozen waivers, in canonical order.
    #[must_use]
    pub fn waivers(&self) -> &[Waiver] {
        &self.draft.waivers
    }

    /// The Five Explicits.
    #[must_use]
    pub const fn explicits(&self) -> FiveExplicits {
        self.draft.explicits
    }

    /// Look up one claim.
    #[must_use]
    pub fn claim(&self, id: &str) -> Option<&ClaimSpec> {
        self.draft.claims.iter().find(|c| c.id == id)
    }

    /// Amend into a successor version. The successor must pass every
    /// freeze gate and carry `version == self.version() + 1`; the record
    /// names exactly the invalidated descendants (claims and obligation
    /// leaves whose content identity changed or vanished).
    ///
    /// # Errors
    ///
    /// [`AmendmentRefusal::WrongVersion`] on a version skip/reuse, or the
    /// successor's own [`FreezeRefusal`].
    pub fn amend(
        &self,
        successor: ManifestDraft,
    ) -> Result<(FrozenManifest, AmendmentRecord), AmendmentRefusal> {
        let expected = self.draft.version + 1;
        if successor.version != expected {
            return Err(AmendmentRefusal::WrongVersion {
                expected,
                offered: successor.version,
            });
        }
        let frozen = successor
            .freeze()
            .map_err(AmendmentRefusal::SuccessorRefused)?;
        let mut invalidated = Vec::new();
        for old in &self.draft.claims {
            let survives = frozen
                .draft
                .claims
                .iter()
                .any(|new| new.id == old.id && claim_digest(new) == claim_digest(old));
            if !survives {
                invalidated.push(old.id.to_string());
            }
        }
        for old in &self.draft.obligations {
            let survives = frozen.draft.obligations.iter().any(|new| {
                new.leaf == old.leaf && obligation_digest(new) == obligation_digest(old)
            });
            if !survives {
                invalidated.push(old.leaf.to_string());
            }
        }
        invalidated.sort();
        invalidated.dedup();
        let record = AmendmentRecord {
            from_version: self.draft.version,
            to_version: frozen.draft.version,
            from_digest: self.digest,
            to_digest: frozen.digest,
            invalidated,
        };
        Ok((frozen, record))
    }
}

fn check_caps(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    let caps: [(&'static str, usize, usize); 4] = [
        ("claims", draft.claims.len(), MAX_CLAIMS),
        ("fixtures", draft.fixtures.len(), MAX_FIXTURES),
        ("obligations", draft.obligations.len(), MAX_OBLIGATIONS),
        ("waivers", draft.waivers.len(), MAX_WAIVERS),
    ];
    for (what, len, cap) in caps {
        if len > cap {
            return Err(FreezeRefusal::OverCap { what, len, cap });
        }
    }
    // Per-row list caps, still before any content scan.
    for claim in &draft.claims {
        if claim.hypotheses.len() > MAX_ROW_ITEMS {
            return Err(FreezeRefusal::OverCap {
                what: "claim.hypotheses",
                len: claim.hypotheses.len(),
                cap: MAX_ROW_ITEMS,
            });
        }
    }
    for row in &draft.obligations {
        let lists: [(&'static str, usize); 4] = [
            ("obligation.claims_covered", row.claims_covered.len()),
            ("obligation.unit_cases", row.unit_cases.len()),
            ("obligation.decks", row.decks.len()),
            ("obligation.g3_relations", row.g3_relations.len()),
        ];
        for (what, len) in lists {
            if len > MAX_ROW_ITEMS {
                return Err(FreezeRefusal::OverCap {
                    what,
                    len,
                    cap: MAX_ROW_ITEMS,
                });
            }
        }
        if row.obs_events.len() > MAX_ROW_ITEMS {
            return Err(FreezeRefusal::OverCap {
                what: "obligation.obs_events",
                len: row.obs_events.len(),
                cap: MAX_ROW_ITEMS,
            });
        }
    }
    Ok(())
}

fn blank(value: &str) -> bool {
    value.trim().is_empty()
}

fn check_top_level_text(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    let id = draft.initiative.to_string();
    let checks: [(&str, &'static str); 8] = [
        (draft.initiative, "initiative"),
        (draft.title, "title"),
        (draft.explicits.units, "explicits.units"),
        (draft.explicits.seeds, "explicits.seeds"),
        (draft.explicits.budgets, "explicits.budgets"),
        (draft.explicits.versions, "explicits.versions"),
        (draft.explicits.capabilities, "explicits.capabilities"),
        (draft.amendment_rules, "amendment_rules"),
    ];
    for (value, field) in checks {
        if blank(value) {
            return Err(FreezeRefusal::BlankField {
                id: id.clone(),
                field,
            });
        }
    }
    Ok(())
}

fn check_components(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    for claim in &draft.claims {
        let checks: [(&str, &'static str); 9] = [
            (claim.id, "claim.id"),
            (claim.statement, "claim.statement"),
            (claim.qoi, "claim.qoi"),
            (claim.unit, "claim.unit"),
            (claim.oracle.identity, "claim.oracle.identity"),
            (claim.activation, "claim.activation"),
            (claim.kill, "claim.kill"),
            (claim.fallback, "claim.fallback"),
            (claim.no_claim, "claim.no_claim"),
        ];
        for (value, field) in checks {
            if blank(value) {
                return Err(FreezeRefusal::BlankField {
                    id: claim.id.to_string(),
                    field,
                });
            }
        }
        if blank(claim.oracle.tcb_overlap) {
            return Err(FreezeRefusal::BlankField {
                id: claim.id.to_string(),
                field: "claim.oracle.tcb_overlap",
            });
        }
        if claim.hypotheses.is_empty() {
            return Err(FreezeRefusal::BlankField {
                id: claim.id.to_string(),
                field: "claim.hypotheses",
            });
        }
        for hypothesis in claim.hypotheses {
            if blank(hypothesis) {
                return Err(FreezeRefusal::BlankField {
                    id: claim.id.to_string(),
                    field: "claim.hypotheses[]",
                });
            }
        }
    }
    for row in &draft.obligations {
        let checks: [(&str, &'static str); 7] = [
            (row.leaf, "obligation.leaf"),
            (row.g0, "obligation.g0"),
            (row.g4_schedule, "obligation.g4_schedule"),
            (row.g5_matrix, "obligation.g5_matrix"),
            (row.entry_point, "obligation.entry_point"),
            (row.dsr_lane, "obligation.dsr_lane"),
            (row.replay_command, "obligation.replay_command"),
        ];
        for (value, field) in checks {
            if blank(value) {
                return Err(FreezeRefusal::BlankField {
                    id: row.leaf.to_string(),
                    field,
                });
            }
        }
        if row.claims_covered.is_empty() {
            return Err(FreezeRefusal::BlankField {
                id: row.leaf.to_string(),
                field: "obligation.claims_covered",
            });
        }
        if row.unit_cases.is_empty() {
            return Err(FreezeRefusal::BlankField {
                id: row.leaf.to_string(),
                field: "obligation.unit_cases",
            });
        }
    }
    for waiver in &draft.waivers {
        let checks: [(&str, &'static str); 6] = [
            (waiver.subject, "waiver.subject"),
            (waiver.reason, "waiver.reason"),
            (waiver.owner, "waiver.owner"),
            (waiver.predicate, "waiver.predicate"),
            (waiver.expiry, "waiver.expiry"),
            (waiver.promotion_effect, "waiver.promotion_effect"),
        ];
        for (value, field) in checks {
            if blank(value) {
                return Err(FreezeRefusal::BlankField {
                    id: waiver.subject.to_string(),
                    field,
                });
            }
        }
    }
    for fixture in &draft.fixtures {
        if blank(fixture.id) {
            return Err(FreezeRefusal::BlankField {
                id: fixture.id.to_string(),
                field: "fixture.id",
            });
        }
    }
    Ok(())
}

fn find_duplicate(mut ids: Vec<&str>) -> Option<String> {
    ids.sort_unstable();
    for pair in ids.windows(2) {
        if pair[0] == pair[1] {
            return Some(pair[0].to_string());
        }
    }
    None
}

fn check_duplicates(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    if let Some(id) = find_duplicate(draft.claims.iter().map(|c| c.id).collect()) {
        return Err(FreezeRefusal::DuplicateId { kind: "claim", id });
    }
    if let Some(id) = find_duplicate(draft.fixtures.iter().map(|x| x.id).collect()) {
        return Err(FreezeRefusal::DuplicateId {
            kind: "fixture",
            id,
        });
    }
    if let Some(id) = find_duplicate(draft.obligations.iter().map(|o| o.leaf).collect()) {
        return Err(FreezeRefusal::DuplicateId {
            kind: "obligation",
            id,
        });
    }
    if let Some(id) = find_duplicate(draft.waivers.iter().map(|w| w.subject).collect()) {
        return Err(FreezeRefusal::DuplicateId { kind: "waiver", id });
    }
    Ok(())
}

fn check_claims(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    for claim in &draft.claims {
        if !claim.oracle.independent {
            return Err(FreezeRefusal::ProductionOracleReuse {
                claim: claim.id.to_string(),
            });
        }
        check_tolerance(claim)?;
    }
    Ok(())
}

fn check_tolerance(claim: &ClaimSpec) -> Result<(), FreezeRefusal> {
    let refuse = |reason: &'static str| FreezeRefusal::InvalidTolerance {
        claim: claim.id.to_string(),
        reason,
    };
    match claim.tolerance {
        ToleranceSemantics::Absolute { atol } => {
            if !atol.is_finite() || atol <= 0.0 {
                return Err(refuse("absolute tolerance must be finite and > 0"));
            }
        }
        ToleranceSemantics::Relative { rtol } => {
            if !rtol.is_finite() || rtol <= 0.0 {
                return Err(refuse("relative tolerance must be finite and > 0"));
            }
        }
        ToleranceSemantics::AbsRel { atol, rtol } => {
            if !atol.is_finite() || !rtol.is_finite() {
                return Err(refuse("non-finite tolerance"));
            }
            if atol < 0.0 || rtol < 0.0 {
                return Err(refuse("negative tolerance"));
            }
            if atol == 0.0 && rtol == 0.0 {
                return Err(refuse("zero-width tolerance (use Exact for bit verdicts)"));
            }
        }
        ToleranceSemantics::Interval { lo, hi } => {
            if !lo.is_finite() || !hi.is_finite() {
                return Err(refuse("non-finite interval bound"));
            }
            if lo > hi {
                return Err(refuse("inverted interval"));
            }
        }
        ToleranceSemantics::Exact => {}
    }
    Ok(())
}

fn check_fixtures(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    for fixture in &draft.fixtures {
        if fixture.digest().is_none() {
            return Err(FreezeRefusal::MalformedFixture {
                fixture: fixture.id.to_string(),
            });
        }
    }
    Ok(())
}

fn check_references(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    let claim_ids: Vec<&str> = draft.claims.iter().map(|c| c.id).collect();
    let fixture_ids: Vec<&str> = draft.fixtures.iter().map(|x| x.id).collect();
    let waiver_subjects: Vec<&str> = draft.waivers.iter().map(|w| w.subject).collect();
    for row in &draft.obligations {
        for claim in row.claims_covered {
            if !claim_ids.contains(claim) {
                return Err(FreezeRefusal::OrphanClaimRef {
                    leaf: row.leaf.to_string(),
                    claim: (*claim).to_string(),
                });
            }
        }
        for deck in row.decks {
            if !fixture_ids.contains(deck) && !waiver_subjects.contains(deck) {
                return Err(FreezeRefusal::OrphanDeck {
                    leaf: row.leaf.to_string(),
                    deck: (*deck).to_string(),
                });
            }
        }
    }
    for claim in &draft.claims {
        let covered = draft
            .obligations
            .iter()
            .any(|row| row.claims_covered.contains(&claim.id));
        let waived = waiver_subjects.contains(&claim.id);
        if !covered && !waived {
            return Err(FreezeRefusal::UncoveredClaim {
                claim: claim.id.to_string(),
            });
        }
    }
    Ok(())
}

/// Sort every collection into one total, content-tie-broken order so
/// assembly order can never move the digest.
fn canonicalize(mut draft: ManifestDraft) -> ManifestDraft {
    draft
        .claims
        .sort_by_cached_key(|c| (c.id, claim_digest(c).0));
    draft
        .fixtures
        .sort_by_cached_key(|x| (x.id, fixture_digest(x).0));
    draft
        .obligations
        .sort_by_cached_key(|o| (o.leaf, obligation_digest(o).0));
    draft
        .waivers
        .sort_by_cached_key(|w| (w.subject, waiver_digest(w).0));
    draft
}

/// Length-frame `bytes` into `out` (u64 LE length + bytes).
fn frame(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

fn frame_list(out: &mut Vec<u8>, items: &[&str]) {
    out.extend_from_slice(&(items.len() as u64).to_le_bytes());
    for item in items {
        frame(out, item.as_bytes());
    }
}

fn tolerance_bytes(out: &mut Vec<u8>, tolerance: ToleranceSemantics) {
    match tolerance {
        ToleranceSemantics::Absolute { atol } => {
            out.push(1);
            out.extend_from_slice(&atol.to_bits().to_le_bytes());
        }
        ToleranceSemantics::Relative { rtol } => {
            out.push(2);
            out.extend_from_slice(&rtol.to_bits().to_le_bytes());
        }
        ToleranceSemantics::AbsRel { atol, rtol } => {
            out.push(3);
            out.extend_from_slice(&atol.to_bits().to_le_bytes());
            out.extend_from_slice(&rtol.to_bits().to_le_bytes());
        }
        ToleranceSemantics::Interval { lo, hi } => {
            out.push(4);
            out.extend_from_slice(&lo.to_bits().to_le_bytes());
            out.extend_from_slice(&hi.to_bits().to_le_bytes());
        }
        ToleranceSemantics::Exact => out.push(5),
    }
}

/// The canonical identity of one claim (every semantic field is
/// mutation-sensitive; floats contribute exact IEEE-754 bits).
#[must_use]
pub fn claim_digest(claim: &ClaimSpec) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, claim.id.as_bytes());
    payload.push(claim.ambition.byte());
    payload.push(claim.polarity.byte());
    frame(&mut payload, claim.statement.as_bytes());
    frame_list(&mut payload, claim.hypotheses);
    frame(&mut payload, claim.qoi.as_bytes());
    frame(&mut payload, claim.unit.as_bytes());
    tolerance_bytes(&mut payload, claim.tolerance);
    payload.push(claim.evidence_tier.byte());
    frame(&mut payload, claim.oracle.identity.as_bytes());
    payload.push(u8::from(claim.oracle.independent));
    frame(&mut payload, claim.oracle.tcb_overlap.as_bytes());
    frame(&mut payload, claim.activation.as_bytes());
    frame(&mut payload, claim.kill.as_bytes());
    frame(&mut payload, claim.fallback.as_bytes());
    frame(&mut payload, claim.no_claim.as_bytes());
    hash_domain(CLAIM_IDENTITY_DOMAIN, &payload)
}

/// The canonical identity of one fixture pin.
#[must_use]
pub fn fixture_digest(fixture: &FixturePin) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, fixture.id.as_bytes());
    payload.push(fixture.partition.byte());
    match fixture.source {
        FixtureSource::AuthoredSpec { spec } => {
            payload.push(1);
            let digest = hash_domain(SPEC_IDENTITY_DOMAIN, spec.as_bytes());
            payload.extend_from_slice(digest.as_bytes());
        }
        // Well-formed external hex normalizes to raw bytes; malformed hex
        // keeps its raw text under a DISTINCT tag (the malformed state is
        // itself mutation-sensitive; freeze refuses it anyway).
        FixtureSource::External { digest_hex } => match ContentHash::from_hex(digest_hex) {
            Some(digest) => {
                payload.push(2);
                payload.extend_from_slice(digest.as_bytes());
            }
            None => {
                payload.push(3);
                frame(&mut payload, digest_hex.as_bytes());
            }
        },
    }
    hash_domain(FIXTURE_IDENTITY_DOMAIN, &payload)
}

/// The canonical identity of one obligation row.
#[must_use]
pub fn obligation_digest(row: &ObligationRow) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, row.leaf.as_bytes());
    frame_list(&mut payload, row.claims_covered);
    frame_list(&mut payload, row.unit_cases);
    frame(&mut payload, row.g0.as_bytes());
    frame_list(&mut payload, row.decks);
    frame_list(&mut payload, row.g3_relations);
    frame(&mut payload, row.g4_schedule.as_bytes());
    frame(&mut payload, row.g5_matrix.as_bytes());
    frame(&mut payload, row.entry_point.as_bytes());
    payload.push(row.tier.byte());
    frame(&mut payload, row.dsr_lane.as_bytes());
    frame_list(&mut payload, row.obs_events);
    frame(&mut payload, row.replay_command.as_bytes());
    hash_domain(OBLIGATION_IDENTITY_DOMAIN, &payload)
}

/// The canonical identity of one waiver.
#[must_use]
pub fn waiver_digest(waiver: &Waiver) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, waiver.subject.as_bytes());
    frame(&mut payload, waiver.reason.as_bytes());
    frame(&mut payload, waiver.owner.as_bytes());
    frame(&mut payload, waiver.predicate.as_bytes());
    frame(&mut payload, waiver.expiry.as_bytes());
    frame(&mut payload, waiver.promotion_effect.as_bytes());
    hash_domain(WAIVER_IDENTITY_DOMAIN, &payload)
}

/// The whole-manifest canonical identity: schema version, top-level
/// fields, then the canonically sorted component digests.
fn manifest_digest(draft: &ManifestDraft) -> ContentHash {
    let mut claim_hashes: Vec<[u8; 32]> = draft.claims.iter().map(|c| claim_digest(c).0).collect();
    let mut fixture_hashes: Vec<[u8; 32]> =
        draft.fixtures.iter().map(|x| fixture_digest(x).0).collect();
    let mut obligation_hashes: Vec<[u8; 32]> = draft
        .obligations
        .iter()
        .map(|o| obligation_digest(o).0)
        .collect();
    let mut waiver_hashes: Vec<[u8; 32]> =
        draft.waivers.iter().map(|w| waiver_digest(w).0).collect();
    claim_hashes.sort_unstable();
    fixture_hashes.sort_unstable();
    obligation_hashes.sort_unstable();
    waiver_hashes.sort_unstable();

    let mut payload = Vec::new();
    payload.extend_from_slice(&VMANIFEST_SCHEMA_VERSION.to_le_bytes());
    frame(&mut payload, draft.initiative.as_bytes());
    frame(&mut payload, draft.title.as_bytes());
    payload.extend_from_slice(&draft.version.to_le_bytes());
    frame(&mut payload, draft.explicits.units.as_bytes());
    frame(&mut payload, draft.explicits.seeds.as_bytes());
    frame(&mut payload, draft.explicits.budgets.as_bytes());
    frame(&mut payload, draft.explicits.versions.as_bytes());
    frame(&mut payload, draft.explicits.capabilities.as_bytes());
    frame(&mut payload, draft.amendment_rules.as_bytes());
    for group in [
        &claim_hashes,
        &fixture_hashes,
        &obligation_hashes,
        &waiver_hashes,
    ] {
        payload.extend_from_slice(&(group.len() as u64).to_le_bytes());
        for hash in group {
            payload.extend_from_slice(hash);
        }
    }
    hash_domain(MANIFEST_IDENTITY_DOMAIN, &payload)
}
