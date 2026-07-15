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
//! - Freezing is fail-closed: [`ManifestDraft::freeze`] refuses empty claim
//!   authority, blank fields, over-cap collections/lists/cumulative text
//!   (checked BEFORE any semantic scan), duplicate or cross-kind ambiguous
//!   evidence ids, non-independent oracles, orphan references, uncovered
//!   claims, orphan waivers, invalid tolerances, and malformed digests — each
//!   with a typed [`FreezeRefusal`] naming the gate.
//! - Authority is sealed: [`FrozenManifest`] has no public constructor
//!   and no mutating API; holding one proves the gates ran. Post-freeze
//!   "alteration" is impossible by construction — change happens only
//!   through [`FrozenManifest::amend`], which preserves initiative
//!   identity, requires the representable successor version, and follows
//!   reverse dependencies to name exactly the invalidated predecessor
//!   claim and obligation evidence.
//! - Identity is canonical: components sort into one total order with
//!   content tie-breaks, and the manifest digest is a domain-separated,
//!   length-framed BLAKE3 hash, byte-stable across runs on the same ISA.
//!
//! Preregistration is not proof: a frozen manifest asserts nothing about
//! implementation correctness, and no color or promotion authority is
//! minted here.

pub use fs_blake3::ContentHash;

use fs_blake3::hash_domain;
use std::{collections::BTreeSet, fmt};

mod cp;
mod i01;
mod i02;
mod i03;
mod i04;
mod i05;
mod i06;
mod i08;
mod i09;
mod i10;
mod i11;
mod i12;
mod i15;
mod pd;
mod rl;

pub use cp::cp_draft;
pub use i01::i01_draft;
pub use i02::i02_draft;
pub use i03::i03_draft;
pub use i04::i04_draft;
pub use i05::i05_draft;
pub use i06::i06_draft;
pub use i08::i08_draft;
pub use i09::i09_draft;
pub use i10::i10_draft;
pub use i11::i11_draft;
pub use i12::i12_draft;
pub use i15::i15_draft;
pub use pd::pd_draft;
pub use rl::rl_draft;

/// Manifest schema version (canonical bytes are comparable only within it).
pub const VMANIFEST_SCHEMA_VERSION: u32 = 2;
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
/// Maximum cumulative UTF-8 bytes carried by one manifest. Checked in the
/// cap phase before blank scans, sorting, or hashing.
pub const MAX_MANIFEST_TEXT_BYTES: usize = 16 * 1024 * 1024;

const MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.manifest.v2";
const CLAIM_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.claim.v1";
const FIXTURE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.fixture.v1";
const OBLIGATION_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vmanifest.obligation.v2";
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
#[derive(Debug, Clone, Copy)]
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

/// Tolerance identity uses exact IEEE-754 encodings, matching
/// [`claim_digest`]. In particular, `-0.0` and `+0.0` are distinct authored
/// bounds even though ordinary floating-point comparison treats them as equal.
impl PartialEq for ToleranceSemantics {
    fn eq(&self, other: &Self) -> bool {
        match (*self, *other) {
            (Self::Absolute { atol: left }, Self::Absolute { atol: right })
            | (Self::Relative { rtol: left }, Self::Relative { rtol: right }) => {
                left.to_bits() == right.to_bits()
            }
            (
                Self::AbsRel {
                    atol: left_atol,
                    rtol: left_rtol,
                },
                Self::AbsRel {
                    atol: right_atol,
                    rtol: right_rtol,
                },
            ) => {
                left_atol.to_bits() == right_atol.to_bits()
                    && left_rtol.to_bits() == right_rtol.to_bits()
            }
            (
                Self::Interval {
                    lo: left_lo,
                    hi: left_hi,
                },
                Self::Interval {
                    lo: right_lo,
                    hi: right_hi,
                },
            ) => left_lo.to_bits() == right_lo.to_bits() && left_hi.to_bits() == right_hi.to_bits(),
            (Self::Exact, Self::Exact) => true,
            _ => false,
        }
    }
}

impl Eq for ToleranceSemantics {}

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
#[derive(Debug, Clone, Copy)]
pub enum FixtureSource {
    /// A canonical generator/spec text authored in this crate; the digest
    /// is computed from these exact bytes.
    AuthoredSpec {
        /// The complete canonical spec text.
        spec: &'static str,
    },
    /// External artifact pinned by its BLAKE3 digest (64 hex chars).
    External {
        /// Case-insensitive 64-char hex of the artifact digest. Canonical
        /// identity is computed from the decoded bytes, so hex case is only
        /// presentation.
        digest_hex: &'static str,
    },
}

impl PartialEq for FixtureSource {
    fn eq(&self, other: &Self) -> bool {
        match (*self, *other) {
            (Self::AuthoredSpec { spec: left }, Self::AuthoredSpec { spec: right }) => {
                left == right
            }
            (Self::External { digest_hex: left }, Self::External { digest_hex: right }) => {
                match (ContentHash::from_hex(left), ContentHash::from_hex(right)) {
                    (Some(left), Some(right)) => left == right,
                    (None, None) => left == right,
                    _ => false,
                }
            }
            _ => false,
        }
    }
}

impl Eq for FixtureSource {}

/// One pinned fixture corpus element.
#[derive(Debug, Clone, Copy)]
pub struct FixturePin {
    /// Stable fixture id (referenced by obligation deck lists).
    pub id: &'static str,
    /// Byte identity.
    pub source: FixtureSource,
    /// Development or held-out partition.
    pub partition: Partition,
}

impl PartialEq for FixturePin {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id && self.source == other.source && self.partition == other.partition
    }
}

impl Eq for FixturePin {}

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
#[derive(Debug, Clone, Copy)]
pub struct ObligationRow {
    /// Execution-leaf/cluster id.
    pub leaf: &'static str,
    /// Claim ids this row's evidence feeds (must exist and be unique in the
    /// manifest). This field is a set for identity purposes: presentation
    /// order does not change the obligation digest.
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

impl ObligationRow {
    /// Covered claims in canonical lexical order. Raw field order is authored
    /// presentation only and must not drive execution or serialization.
    #[must_use]
    pub fn canonical_claims_covered(&self) -> Vec<&'static str> {
        canonical_string_set(self.claims_covered)
    }

    /// Unit-case classes in canonical lexical order.
    #[must_use]
    pub fn canonical_unit_cases(&self) -> Vec<&'static str> {
        canonical_string_set(self.unit_cases)
    }

    /// Fixture/deck ids in canonical lexical order.
    #[must_use]
    pub fn canonical_decks(&self) -> Vec<&'static str> {
        canonical_string_set(self.decks)
    }

    /// Metamorphic relations in canonical lexical order.
    #[must_use]
    pub fn canonical_g3_relations(&self) -> Vec<&'static str> {
        canonical_string_set(self.g3_relations)
    }

    /// Observation event kinds in canonical lexical order.
    #[must_use]
    pub fn canonical_obs_events(&self) -> Vec<&'static str> {
        canonical_string_set(self.obs_events)
    }
}

impl PartialEq for ObligationRow {
    fn eq(&self, other: &Self) -> bool {
        same_string_set(self.claims_covered, other.claims_covered)
            && same_obligation_execution_semantics(self, other)
    }
}

impl Eq for ObligationRow {}

/// Canonical, immutable projection of one accepted obligation row.
///
/// Draft rows retain authored presentation order so freeze can diagnose the
/// exact submitted bytes. Frozen rows own lexically sorted set fields, making
/// the only public frozen execution/serialization view canonical by
/// construction rather than by caller convention.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrozenObligationRow {
    /// Execution-leaf/cluster id.
    leaf: &'static str,
    /// Covered claim ids in canonical lexical order.
    claims_covered: Vec<&'static str>,
    /// Required unit-case classes in canonical lexical order.
    unit_cases: Vec<&'static str>,
    /// G0 generators, predicates, laws, shrinkers, and replay seeds.
    g0: &'static str,
    /// Fixture/deck ids in canonical lexical order.
    decks: Vec<&'static str>,
    /// G3 relations in canonical lexical order.
    g3_relations: Vec<&'static str>,
    /// G4 fault/cancellation/checkpoint schedule.
    g4_schedule: &'static str,
    /// G5 determinism matrix.
    g5_matrix: &'static str,
    /// Named campaign entry point.
    entry_point: &'static str,
    /// Smoke/core/max tier.
    tier: CampaignTier,
    /// DSR lane that owns the run.
    dsr_lane: &'static str,
    /// Required observation events in canonical lexical order.
    obs_events: Vec<&'static str>,
    /// Exact replay command.
    replay_command: &'static str,
    /// Canonical authored-row component identity.
    digest: ContentHash,
}

impl FrozenObligationRow {
    fn from_accepted(row: &ObligationRow) -> Self {
        Self {
            leaf: row.leaf,
            claims_covered: row.canonical_claims_covered(),
            unit_cases: row.canonical_unit_cases(),
            g0: row.g0,
            decks: row.canonical_decks(),
            g3_relations: row.canonical_g3_relations(),
            g4_schedule: row.g4_schedule,
            g5_matrix: row.g5_matrix,
            entry_point: row.entry_point,
            tier: row.tier,
            dsr_lane: row.dsr_lane,
            obs_events: row.canonical_obs_events(),
            replay_command: row.replay_command,
            digest: obligation_digest(row),
        }
    }

    /// Execution-leaf/cluster id.
    #[must_use]
    pub const fn leaf(&self) -> &'static str {
        self.leaf
    }

    /// Covered claim ids in canonical lexical order.
    #[must_use]
    pub fn claims_covered(&self) -> &[&'static str] {
        &self.claims_covered
    }

    /// Required unit-case classes in canonical lexical order.
    #[must_use]
    pub fn unit_cases(&self) -> &[&'static str] {
        &self.unit_cases
    }

    /// G0 generator/predicate/law/shrinker/seed contract.
    #[must_use]
    pub const fn g0(&self) -> &'static str {
        self.g0
    }

    /// Fixture/deck ids in canonical lexical order.
    #[must_use]
    pub fn decks(&self) -> &[&'static str] {
        &self.decks
    }

    /// G3 relations in canonical lexical order.
    #[must_use]
    pub fn g3_relations(&self) -> &[&'static str] {
        &self.g3_relations
    }

    /// G4 fault/cancellation/checkpoint schedule.
    #[must_use]
    pub const fn g4_schedule(&self) -> &'static str {
        self.g4_schedule
    }

    /// G5 determinism matrix.
    #[must_use]
    pub const fn g5_matrix(&self) -> &'static str {
        self.g5_matrix
    }

    /// Named campaign entry point.
    #[must_use]
    pub const fn entry_point(&self) -> &'static str {
        self.entry_point
    }

    /// Smoke/core/max tier.
    #[must_use]
    pub const fn tier(&self) -> CampaignTier {
        self.tier
    }

    /// DSR lane that owns the run.
    #[must_use]
    pub const fn dsr_lane(&self) -> &'static str {
        self.dsr_lane
    }

    /// Required observation events in canonical lexical order.
    #[must_use]
    pub fn obs_events(&self) -> &[&'static str] {
        &self.obs_events
    }

    /// Exact replay command.
    #[must_use]
    pub const fn replay_command(&self) -> &'static str {
        self.replay_command
    }

    /// Canonical component identity, equal to the accepted draft row's
    /// [`obligation_digest`].
    #[must_use]
    pub const fn digest(&self) -> ContentHash {
        self.digest
    }
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
    /// Opaque schema, toolchain, dependency, and data-contract version pins.
    ///
    /// [`ManifestDraft::version`] is the only machine-interpreted manifest
    /// instance revision. This string is hashed as provenance but deliberately
    /// not parsed as a second revision authority. Authored instances have a
    /// deliberately narrow, non-exhaustive conformance lint for known legacy
    /// semicolon-field spellings of the numeric revision; arbitrary prose is
    /// not assigned revision semantics. This separation lets a targeted
    /// amendment advance the instance revision without pretending that the
    /// toolchain or every evidence dependency changed.
    pub versions: &'static str,
    /// Capability flags in force.
    pub capabilities: &'static str,
}

/// A mutable manifest under assembly. Freezing consumes it.
#[derive(Debug, Clone, PartialEq)]
pub struct ManifestDraft {
    /// Initiative id (e.g. `"I01"`).
    pub initiative: &'static str,
    /// Identity-bearing campaign-authority title. Display-only labels belong
    /// outside the manifest because changing this field invalidates all
    /// predecessor evidence.
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
/// documented order: all collection/list/cumulative-text caps, version,
/// top-level blanks, required nonempty claim authority and component blanks,
/// duplicates, oracle independence, tolerance validity, fixture
/// well-formedness, orphan references, coverage, and waiver subjects.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FreezeRefusal {
    /// A collection, per-row list, or cumulative text budget exceeds its cap
    /// (checked before any semantic scan).
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
    /// A semantically required collection is empty.
    EmptyCollection {
        /// Required collection name.
        what: &'static str,
    },
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
    /// A waiver subject names no claim, referenced deck slot, or obligation
    /// leaf in the same manifest.
    OrphanWaiver {
        /// The unused or misspelled waiver subject.
        subject: String,
    },
    /// A waiver's untyped subject collides across claim, referenced-deck,
    /// and obligation-leaf namespaces.
    AmbiguousWaiverSubject {
        /// The multiply-resolved waiver subject.
        subject: String,
    },
}

impl fmt::Display for FreezeRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OverCap { what, len, cap } => {
                write!(f, "{what}: {len} exceeds the cap of {cap}")
            }
            Self::ZeroVersion => f.write_str("manifest version must be >= 1"),
            Self::EmptyCollection { what } => {
                write!(f, "required manifest collection '{what}' is empty")
            }
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
            Self::OrphanWaiver { subject } => write!(
                f,
                "waiver subject '{subject}' names no claim, referenced deck slot, or obligation leaf"
            ),
            Self::AmbiguousWaiverSubject { subject } => write!(
                f,
                "waiver subject '{subject}' resolves in more than one of the claim, referenced-deck, and obligation-leaf namespaces"
            ),
        }
    }
}

impl std::error::Error for FreezeRefusal {}

/// Why an amendment was refused.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AmendmentRefusal {
    /// An amendment cannot cross an initiative identity boundary.
    InitiativeChanged {
        /// Initiative carried by the frozen predecessor.
        expected: String,
        /// Initiative offered by the successor draft.
        offered: String,
    },
    /// The predecessor is already at `u32::MAX`, so no successor version
    /// can be represented.
    VersionExhausted {
        /// Unincrementable predecessor version.
        version: u32,
    },
    /// The successor draft must carry `version == predecessor + 1`.
    WrongVersion {
        /// Expected successor version.
        expected: u32,
        /// Offered version.
        offered: u32,
    },
    /// A predecessor claim id cannot be reused as a successor obligation leaf,
    /// or vice versa: invalidation ids would otherwise alias authority kinds.
    EvidenceKindChanged {
        /// Reused evidence authority id.
        id: String,
        /// Predecessor kind.
        from_kind: &'static str,
        /// Successor kind.
        to_kind: &'static str,
    },
    /// The successor draft itself failed its freeze gates.
    SuccessorRefused(FreezeRefusal),
}

impl fmt::Display for AmendmentRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InitiativeChanged { expected, offered } => write!(
                f,
                "amendment cannot change initiative from '{expected}' to '{offered}'"
            ),
            Self::VersionExhausted { version } => write!(
                f,
                "manifest version {version} cannot be incremented for an amendment"
            ),
            Self::WrongVersion { expected, offered } => {
                write!(f, "amendment must carry version {expected}, got {offered}")
            }
            Self::EvidenceKindChanged {
                id,
                from_kind,
                to_kind,
            } => write!(
                f,
                "amendment cannot reuse evidence id '{id}' as {to_kind}; it was a predecessor {from_kind}"
            ),
            Self::SuccessorRefused(refusal) => {
                write!(f, "successor draft refused: {refusal}")
            }
        }
    }
}

impl std::error::Error for AmendmentRefusal {}

/// The record an amendment produces: which version replaced which, and
/// exactly which predecessor claim and obligation-leaf evidence is
/// invalidated after reverse-dependency propagation.
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
#[derive(Debug, Clone)]
pub struct FrozenManifest {
    draft: ManifestDraft,
    obligations: Vec<FrozenObligationRow>,
    digest: ContentHash,
}

/// Frozen-manifest equality is canonical content identity, not the incidental
/// presentation order retained by borrowed string slices inside the draft.
impl PartialEq for FrozenManifest {
    fn eq(&self, other: &Self) -> bool {
        self.digest == other.digest
    }
}

impl Eq for FrozenManifest {}

impl ManifestDraft {
    /// Freeze the draft, consuming it.
    ///
    /// Gate order (documented so refusal tests are stable): all collection,
    /// per-row list, and cumulative-text caps (before any semantic scan),
    /// version, top-level blank fields, required nonempty claim authority,
    /// per-component blank fields, duplicate ids, oracle
    /// independence, tolerance validity, fixture well-formedness, orphan
    /// claim references, orphan decks, uncovered claims, orphan waivers.
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
        let draft = canonicalize(self);
        let obligations = draft
            .obligations
            .iter()
            .map(FrozenObligationRow::from_accepted)
            .collect();
        Ok(FrozenManifest {
            draft,
            obligations,
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
    pub fn obligations(&self) -> &[FrozenObligationRow] {
        &self.obligations
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
    /// freeze gate, preserve the initiative, and carry
    /// `version == self.version() + 1`; the record names exactly the
    /// invalidated predecessor claims and obligation leaves, including
    /// reverse-dependency propagation through claims, fixture decks,
    /// waivers, obligations, and global campaign policy. A title change is
    /// global authority; a numeric version-only successor leaves the set
    /// empty because the revision itself does not falsify component evidence.
    ///
    /// # Errors
    ///
    /// [`AmendmentRefusal::InitiativeChanged`] on a cross-initiative
    /// successor, [`AmendmentRefusal::VersionExhausted`] when the
    /// predecessor version cannot be incremented,
    /// [`AmendmentRefusal::WrongVersion`] on a version skip/reuse,
    /// [`AmendmentRefusal::EvidenceKindChanged`] if a predecessor claim id is
    /// reused as a successor leaf (or vice versa), or the successor's own
    /// [`FreezeRefusal`].
    #[allow(clippy::too_many_lines)]
    pub fn amend(
        &self,
        successor: ManifestDraft,
    ) -> Result<(FrozenManifest, AmendmentRecord), AmendmentRefusal> {
        if successor.initiative != self.draft.initiative {
            return Err(AmendmentRefusal::InitiativeChanged {
                expected: self.draft.initiative.to_string(),
                offered: successor.initiative.to_string(),
            });
        }
        let expected =
            self.draft
                .version
                .checked_add(1)
                .ok_or(AmendmentRefusal::VersionExhausted {
                    version: self.draft.version,
                })?;
        if successor.version != expected {
            return Err(AmendmentRefusal::WrongVersion {
                expected,
                offered: successor.version,
            });
        }
        let frozen = successor
            .freeze()
            .map_err(AmendmentRefusal::SuccessorRefused)?;
        for claim in &self.draft.claims {
            if frozen
                .draft
                .obligations
                .iter()
                .any(|row| row.leaf == claim.id)
            {
                return Err(AmendmentRefusal::EvidenceKindChanged {
                    id: claim.id.to_string(),
                    from_kind: "claim",
                    to_kind: "obligation leaf",
                });
            }
        }
        for row in &self.draft.obligations {
            if frozen.draft.claims.iter().any(|claim| claim.id == row.leaf) {
                return Err(AmendmentRefusal::EvidenceKindChanged {
                    id: row.leaf.to_string(),
                    from_kind: "obligation leaf",
                    to_kind: "claim",
                });
            }
        }
        let mut invalidated = BTreeSet::new();

        if self.draft.title != frozen.draft.title
            || self.draft.explicits != frozen.draft.explicits
            || self.draft.amendment_rules != frozen.draft.amendment_rules
        {
            invalidate_all_predecessor_evidence(&self.draft, &mut invalidated);
        }

        for old in &self.draft.claims {
            let survives = frozen
                .draft
                .claims
                .iter()
                .any(|new| new.id == old.id && claim_digest(new) == claim_digest(old));
            if !survives {
                invalidate_predecessor_claim_change(&self.draft, old.id, &mut invalidated);
            }
        }
        for old in &self.draft.fixtures {
            let survives = frozen
                .draft
                .fixtures
                .iter()
                .any(|new| new.id == old.id && fixture_digest(new) == fixture_digest(old));
            if !survives {
                invalidate_predecessor_deck(&self.draft, old.id, &mut invalidated);
            }
        }
        // An added fixture can replace a predecessor's waived deck slot
        // without changing the obligation row or the waiver. That newly
        // available evidence source still changes every predecessor claim
        // and leaf that named the slot, so additions need the same reverse
        // dependency walk as edits/removals.
        for new in &frozen.draft.fixtures {
            let existed_unchanged = self
                .draft
                .fixtures
                .iter()
                .any(|old| old.id == new.id && fixture_digest(old) == fixture_digest(new));
            if !existed_unchanged {
                invalidate_predecessor_deck(&self.draft, new.id, &mut invalidated);
            }
        }
        for old in &self.draft.obligations {
            let successor_row = frozen
                .draft
                .obligations
                .iter()
                .find(|new| new.leaf == old.leaf);
            let survives =
                successor_row.is_some_and(|new| obligation_digest(new) == obligation_digest(old));
            if !survives {
                invalidated.insert(old.leaf.to_string());
                match successor_row {
                    Some(new) if same_obligation_execution_semantics(old, new) => {
                        // A mapping-only edit invalidates the predecessor
                        // leaf's mapping-bound authority and claims removed
                        // from that producer. It does not revoke unchanged
                        // sibling adjudications: their execution payload can
                        // be rebound only through the amendment lineage.
                        for claim in old
                            .claims_covered
                            .iter()
                            .filter(|claim| !new.claims_covered.contains(claim))
                        {
                            invalidated.insert((*claim).to_string());
                        }
                    }
                    Some(_) | None => {
                        invalidate_predecessor_leaf(&self.draft, old.leaf, &mut invalidated);
                    }
                }
            }
        }
        for new in &frozen.draft.obligations {
            let predecessor_row = self
                .draft
                .obligations
                .iter()
                .find(|old| old.leaf == new.leaf);
            let existed_unchanged =
                predecessor_row.is_some_and(|old| obligation_digest(old) == obligation_digest(new));
            if !existed_unchanged {
                let mapping_only = predecessor_row
                    .is_some_and(|old| same_obligation_execution_semantics(old, new));
                for claim in new.claims_covered {
                    if mapping_only
                        && predecessor_row.is_some_and(|old| old.claims_covered.contains(claim))
                    {
                        continue;
                    }
                    if self.draft.claims.iter().any(|old| old.id == *claim) {
                        invalidated.insert((*claim).to_string());
                    }
                }
            }
        }
        for old in &self.draft.waivers {
            let survives =
                frozen.draft.waivers.iter().any(|new| {
                    new.subject == old.subject && waiver_digest(new) == waiver_digest(old)
                });
            if !survives {
                invalidate_predecessor_subject(&self.draft, old.subject, &mut invalidated);
            }
        }
        for new in &frozen.draft.waivers {
            let existed_unchanged =
                self.draft.waivers.iter().any(|old| {
                    old.subject == new.subject && waiver_digest(old) == waiver_digest(new)
                });
            if !existed_unchanged {
                invalidate_predecessor_subject(&self.draft, new.subject, &mut invalidated);
            }
        }
        let record = AmendmentRecord {
            from_version: self.draft.version,
            to_version: frozen.draft.version,
            from_digest: self.digest,
            to_digest: frozen.digest,
            invalidated: invalidated.into_iter().collect(),
        };
        Ok((frozen, record))
    }
}

fn canonical_string_set<'a>(items: &[&'a str]) -> Vec<&'a str> {
    let mut sorted = items.to_vec();
    sorted.sort_unstable();
    sorted
}

fn same_string_set(left: &[&str], right: &[&str]) -> bool {
    if left.len() != right.len() {
        return false;
    }
    canonical_string_set(left) == canonical_string_set(right)
}

/// Whether two rows differ only in which claim ids consume the same
/// execution payload. The list-valued obligation fields are mathematical or
/// policy sets: freeze rejects duplicates, and presentation order carries no
/// identity. Keeping this distinction prevents a claim removal or rename from
/// needlessly revoking an unchanged sibling adjudication, while any
/// executable, deck, oracle-adjacent, policy, or observability change still
/// invalidates every claim fed by the predecessor leaf.
fn same_obligation_execution_semantics(left: &ObligationRow, right: &ObligationRow) -> bool {
    left.leaf == right.leaf
        && same_string_set(left.unit_cases, right.unit_cases)
        && left.g0 == right.g0
        && same_string_set(left.decks, right.decks)
        && same_string_set(left.g3_relations, right.g3_relations)
        && left.g4_schedule == right.g4_schedule
        && left.g5_matrix == right.g5_matrix
        && left.entry_point == right.entry_point
        && left.tier == right.tier
        && left.dsr_lane == right.dsr_lane
        && same_string_set(left.obs_events, right.obs_events)
        && left.replay_command == right.replay_command
}

fn invalidate_all_predecessor_evidence(draft: &ManifestDraft, invalidated: &mut BTreeSet<String>) {
    invalidated.extend(draft.claims.iter().map(|claim| claim.id.to_string()));
    invalidated.extend(
        draft
            .obligations
            .iter()
            .map(|obligation| obligation.leaf.to_string()),
    );
}

/// A claim-content change invalidates that claim and each predecessor leaf
/// producing its evidence, but deliberately not unrelated sibling claims
/// covered by the same leaf.
fn invalidate_predecessor_claim_change(
    draft: &ManifestDraft,
    claim: &str,
    invalidated: &mut BTreeSet<String>,
) {
    if !draft.claims.iter().any(|candidate| candidate.id == claim) {
        return;
    }
    invalidated.insert(claim.to_string());
    invalidated.extend(
        draft
            .obligations
            .iter()
            .filter(|row| row.claims_covered.contains(&claim))
            .map(|row| row.leaf.to_string()),
    );
}

/// An obligation-content change invalidates the predecessor leaf and every
/// predecessor claim whose evidence that leaf produced.
fn invalidate_predecessor_leaf(
    draft: &ManifestDraft,
    leaf: &str,
    invalidated: &mut BTreeSet<String>,
) {
    let Some(row) = draft.obligations.iter().find(|row| row.leaf == leaf) else {
        return;
    };
    invalidated.insert(row.leaf.to_string());
    for claim in row.claims_covered {
        if draft.claims.iter().any(|candidate| candidate.id == *claim) {
            invalidated.insert((*claim).to_string());
        }
    }
}

/// A changed fixture/deck invalidates every predecessor consumer leaf and
/// every predecessor claim fed by those leaves.
fn invalidate_predecessor_deck(
    draft: &ManifestDraft,
    deck: &str,
    invalidated: &mut BTreeSet<String>,
) {
    for leaf in draft
        .obligations
        .iter()
        .filter(|row| row.decks.contains(&deck))
        .map(|row| row.leaf)
    {
        invalidate_predecessor_leaf(draft, leaf, invalidated);
    }
}

/// Waiver subjects may name claims, deck slots, or obligation leaves. General
/// ids can overlap across those namespaces, but freeze requires each untyped
/// waiver subject to resolve in exactly one of them.
fn invalidate_predecessor_subject(
    draft: &ManifestDraft,
    subject: &str,
    invalidated: &mut BTreeSet<String>,
) {
    invalidate_predecessor_claim_change(draft, subject, invalidated);
    invalidate_predecessor_deck(draft, subject, invalidated);
    invalidate_predecessor_leaf(draft, subject, invalidated);
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
    let text_bytes = manifest_text_bytes(draft).unwrap_or(usize::MAX);
    if text_bytes > MAX_MANIFEST_TEXT_BYTES {
        return Err(FreezeRefusal::OverCap {
            what: "manifest text bytes",
            len: text_bytes,
            cap: MAX_MANIFEST_TEXT_BYTES,
        });
    }
    Ok(())
}

fn manifest_text_bytes(draft: &ManifestDraft) -> Option<usize> {
    fn add(total: &mut usize, value: &str) -> Option<()> {
        *total = total.checked_add(value.len())?;
        Some(())
    }

    let mut total = 0usize;
    for value in [
        draft.initiative,
        draft.title,
        draft.explicits.units,
        draft.explicits.seeds,
        draft.explicits.budgets,
        draft.explicits.versions,
        draft.explicits.capabilities,
        draft.amendment_rules,
    ] {
        add(&mut total, value)?;
    }
    for claim in &draft.claims {
        for value in [
            claim.id,
            claim.statement,
            claim.qoi,
            claim.unit,
            claim.oracle.identity,
            claim.oracle.tcb_overlap,
            claim.activation,
            claim.kill,
            claim.fallback,
            claim.no_claim,
        ] {
            add(&mut total, value)?;
        }
        for hypothesis in claim.hypotheses {
            add(&mut total, hypothesis)?;
        }
    }
    for fixture in &draft.fixtures {
        add(&mut total, fixture.id)?;
        match fixture.source {
            FixtureSource::AuthoredSpec { spec } => add(&mut total, spec)?,
            FixtureSource::External { digest_hex } => add(&mut total, digest_hex)?,
        }
    }
    for row in &draft.obligations {
        for value in [
            row.leaf,
            row.g0,
            row.g4_schedule,
            row.g5_matrix,
            row.entry_point,
            row.dsr_lane,
            row.replay_command,
        ] {
            add(&mut total, value)?;
        }
        for list in [
            row.claims_covered,
            row.unit_cases,
            row.decks,
            row.g3_relations,
            row.obs_events,
        ] {
            for value in list {
                add(&mut total, value)?;
            }
        }
    }
    for waiver in &draft.waivers {
        for value in [
            waiver.subject,
            waiver.reason,
            waiver.owner,
            waiver.predicate,
            waiver.expiry,
            waiver.promotion_effect,
        ] {
            add(&mut total, value)?;
        }
    }
    Some(total)
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

#[allow(clippy::too_many_lines)]
fn check_components(draft: &ManifestDraft) -> Result<(), FreezeRefusal> {
    if draft.claims.is_empty() {
        return Err(FreezeRefusal::EmptyCollection { what: "claims" });
    }
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
        let required_lists: [(&[&str], &'static str); 5] = [
            (row.claims_covered, "obligation.claims_covered"),
            (row.unit_cases, "obligation.unit_cases"),
            (row.decks, "obligation.decks"),
            (row.g3_relations, "obligation.g3_relations"),
            (row.obs_events, "obligation.obs_events"),
        ];
        for (items, field) in required_lists {
            if items.is_empty() {
                return Err(FreezeRefusal::BlankField {
                    id: row.leaf.to_string(),
                    field,
                });
            }
            if items.iter().any(|item| blank(item)) {
                return Err(FreezeRefusal::BlankField {
                    id: row.leaf.to_string(),
                    field,
                });
            }
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
    if let Some(id) = find_duplicate(
        draft
            .claims
            .iter()
            .map(|claim| claim.id)
            .chain(draft.obligations.iter().map(|row| row.leaf))
            .collect(),
    ) {
        return Err(FreezeRefusal::DuplicateId {
            kind: "claim/obligation evidence",
            id,
        });
    }
    for row in &draft.obligations {
        let set_fields: [(&[&str], &'static str); 5] = [
            (row.claims_covered, "obligation claim mapping"),
            (row.unit_cases, "obligation unit case"),
            (row.decks, "obligation deck"),
            (row.g3_relations, "obligation G3 relation"),
            (row.obs_events, "obligation observation event"),
        ];
        for (items, kind) in set_fields {
            if let Some(id) = find_duplicate(items.to_vec()) {
                return Err(FreezeRefusal::DuplicateId { kind, id });
            }
        }
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
    }
    for claim in &draft.claims {
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
    let leaf_ids: Vec<&str> = draft.obligations.iter().map(|row| row.leaf).collect();
    let deck_ids: Vec<&str> = draft
        .obligations
        .iter()
        .flat_map(|row| row.decks.iter().copied())
        .collect();
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
    }
    for row in &draft.obligations {
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
    for subject in waiver_subjects {
        let namespace_count = claim_ids.contains(&subject) as usize
            + leaf_ids.contains(&subject) as usize
            + deck_ids.contains(&subject) as usize;
        if namespace_count == 0 {
            return Err(FreezeRefusal::OrphanWaiver {
                subject: subject.to_string(),
            });
        }
        if namespace_count > 1 {
            return Err(FreezeRefusal::AmbiguousWaiverSubject {
                subject: subject.to_string(),
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

fn frame_string_set(out: &mut Vec<u8>, items: &[&str]) {
    let sorted = canonical_string_set(items);
    frame_list(out, &sorted);
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

/// The canonical identity of one obligation row. Required-case, deck,
/// metamorphic-relation, observation-event, and claim-mapping lists are sets;
/// their presentation order does not change identity.
#[must_use]
pub fn obligation_digest(row: &ObligationRow) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, row.leaf.as_bytes());
    frame_string_set(&mut payload, row.claims_covered);
    frame_string_set(&mut payload, row.unit_cases);
    frame(&mut payload, row.g0.as_bytes());
    frame_string_set(&mut payload, row.decks);
    frame_string_set(&mut payload, row.g3_relations);
    frame(&mut payload, row.g4_schedule.as_bytes());
    frame(&mut payload, row.g5_matrix.as_bytes());
    frame(&mut payload, row.entry_point.as_bytes());
    payload.push(row.tier.byte());
    frame(&mut payload, row.dsr_lane.as_bytes());
    frame_string_set(&mut payload, row.obs_events);
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
