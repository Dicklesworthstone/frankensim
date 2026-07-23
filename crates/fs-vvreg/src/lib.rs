//! fs-vvreg — the Gauntlet G1/G2 benchmark & V&V registry (bead
//! frankensim-ext-benchmark-vv-registry-f1gq; epic E0c). Layer: UTIL.
//!
//! A family name (TEAM, NAFEMS, CFR, IFToMM, ECN) is NOT an executable
//! benchmark. Before any solver claims a G1/G2 result against a registry
//! entry, that entry must pin: exact version/edition, source, license,
//! input-deck identity, oracle binding, QoIs, and acceptance envelopes. Citation is
//! fail-closed: [`Registry::cite`] returns a typed [`CitationRefusal`]
//! naming the first missing field instead of letting an unpinned family
//! name act as an oracle, and only the seeded registry behind
//! [`registry()`] can mint a [`CitationReceipt`] — caller-built rows and
//! registries are validation/lint-only ([`validate_entry`]).
//!
//! The registry also seeds the primary-reference table: literature and
//! standards anchors for definitions and benchmark provenance. References
//! are never authority-by-citation — no API on [`PrimaryReference`] mints
//! an evidence color; each consuming bead records its own
//! [`ConsumptionStatus`] (Appendix-D discipline) and pins the exact
//! artifact version through the entry digest.
//!
//! Color rule (encoded on [`CitationReceipt`]): a standards or benchmark
//! calculation earns at most numerical `Verified` for the exact edition and
//! scope; its physical prediction is `Estimated` unless independent
//! held-out evidence for the named QoI and population earns `Validated`.
//! No color is inherited from a publisher's name.

pub use fs_blake3::ContentHash;
pub use fs_evidence::ColorRank;

use fs_blake3::hash_domain;
use std::fmt;
use std::sync::LazyLock;

mod seed;

/// Versioned, evidence-bearing validation datasets with fail-closed partition
/// and context-of-use queries.
pub mod corpus;

/// Purpose-typed corpus access, transitive calibration taint, and versioned
/// repartition receipts.
pub mod partition;

/// Claim-scoped external-evidence axes and fail-closed portfolio admission.
pub mod portfolio;

/// Reference-only Level-A thermal analytic values and G1 order targets.
pub mod thermal_level_a;

/// Versioned standards-edition, source-lineage, and derived-rule manifest.
pub mod standards;

/// The registry payload version (rows are only comparable within a version).
pub const VVREG_VERSION: u32 = 1;
/// Maximum QoIs accepted on one entry (checked before the duplicate-name
/// scan, so gate cost is bounded even for hostile rows).
pub const MAX_QOIS_PER_ENTRY: usize = 64;
/// Canonical length-framed registry-identity schema version.
pub const REGISTRY_IDENTITY_SCHEMA_VERSION: u32 = 1;
/// Canonical executable-envelope diagnostic schema version.
pub const ENVELOPE_VERDICT_SCHEMA_VERSION: u32 = 1;

const REGISTRY_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.registry.v1";
const ENTRY_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.entry.v1";
const DECK_IDENTITY_DOMAIN: &str = "org.frankensim.fs-vvreg.deck.v1";

/// Which Gauntlet tier an entry serves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegistryTier {
    /// G1: analytic/manufactured oracle with a closed-form or exactly
    /// parameterized reference.
    G1Analytic,
    /// G2: canonical external benchmark with a pinned deck artifact.
    G2Benchmark,
}

impl RegistryTier {
    /// Stable tag for canonical rows.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            Self::G1Analytic => "G1",
            Self::G2Benchmark => "G2",
        }
    }

    const fn byte(self) -> u8 {
        match self {
            Self::G1Analytic => 1,
            Self::G2Benchmark => 2,
        }
    }
}

/// Exact version/edition pin. A family or problem number alone is not an
/// edition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Edition {
    /// The exact edition/revision/version string the deck was taken from.
    Exact {
        /// Exact edition identifier (e.g. `"ISO 6336-1:2019"`).
        version: &'static str,
    },
    /// Not yet pinned — blocks citation.
    Unpinned,
}

/// License/terms state for the deck artifact. Determines whether the deck
/// bytes may live in-repo or must live in quarantined external storage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LicenseState {
    /// SPDX-identified license; artifact bytes may be stored in-repo.
    Spdx {
        /// SPDX identifier (e.g. `"MIT"`, `"CC-BY-4.0"`).
        id: &'static str,
    },
    /// Named restricted terms; bytes stored out-of-repo at the named
    /// location, only the digest is registered.
    Restricted {
        /// Human-readable terms name (e.g. `"NAFEMS member license"`).
        terms: &'static str,
        /// Storage convention for the licensed bytes.
        storage: &'static str,
    },
    /// Not yet established — blocks citation.
    Unpinned,
}

/// The pinned input-deck identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeckPin {
    /// A canonical spec authored in this crate: the exact parameterization,
    /// closed form, and assumption set, embedded as bytes and
    /// domain-hashed. Used by G1 analytic entries.
    AuthoredSpec {
        /// The complete canonical spec text (the hashed deck bytes).
        spec: &'static str,
    },
    /// An external artifact pinned by its BLAKE3 content digest (64 hex
    /// chars); the bytes live wherever the license permits.
    External {
        /// Lowercase 64-char hex of the artifact's `ContentHash`.
        digest_hex: &'static str,
    },
    /// Not yet pinned — blocks citation.
    Unpinned,
}

impl DeckPin {
    /// The deck's content identity, if pinned and well-formed. A blank
    /// authored spec and a malformed external hex are NOT well-formed and
    /// yield `None`, matching the admission gates.
    #[must_use]
    pub fn digest(&self) -> Option<ContentHash> {
        match self {
            Self::AuthoredSpec { spec } => {
                if spec.trim().is_empty() {
                    None
                } else {
                    Some(hash_domain(DECK_IDENTITY_DOMAIN, spec.as_bytes()))
                }
            }
            Self::External { digest_hex } => ContentHash::from_hex(digest_hex),
            Self::Unpinned => None,
        }
    }
}

/// Acceptance envelope for one quantity of interest.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AcceptanceEnvelope {
    /// Accept when `|candidate - reference| <= atol + rtol * |reference|`.
    Tolerance {
        /// Absolute tolerance (finite, `>= 0`).
        atol: f64,
        /// Relative tolerance (finite, `>= 0`).
        rtol: f64,
    },
    /// Accept when the candidate lies in `[lo, hi]` (finite, `lo <= hi`).
    Interval {
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
    /// Not yet pinned — blocks citation.
    Unpinned,
}

/// Whether the deck's oracle is executable from the deck text alone.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OracleBinding {
    /// The oracle identity or executable comparison procedure has not yet
    /// been pinned. This is a target declaration, never a citable oracle.
    Unpinned,
    /// The deck contains the complete closed form or uniquely determined
    /// procedure; a consumer can implement it from the deck alone.
    SelfContained,
    /// The deck deliberately delegates a load-bearing derivation to the
    /// consumer (the anti-mnemonic discipline). Such an entry is a
    /// registered target but stays NON-CITABLE until a future derivation
    /// receipt mechanism binds the obligation.
    DerivationRequired {
        /// The delegated obligation, verbatim.
        obligation: &'static str,
    },
}

/// A named quantity of interest with its unit and acceptance envelope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Qoi {
    /// Stable QoI name (snake_case).
    pub name: &'static str,
    /// Unit expression (SI symbol string; `"1"` for dimensionless).
    pub unit: &'static str,
    /// Acceptance envelope for this QoI.
    pub envelope: AcceptanceEnvelope,
}

/// One computed quantity offered to an [`AcceptanceEnvelope`]. The variant is
/// explicit so a tolerance gate cannot silently treat an interval bound as an
/// oracle reference (or vice versa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnvelopeObservation {
    /// Compare a computed value with an independently supplied reference.
    AgainstReference {
        /// Reference/oracle value for this exact QoI.
        reference: f64,
        /// Computed value being gated.
        computed: f64,
    },
    /// Compare a computed value directly with the registered interval.
    AgainstInterval {
        /// Computed value being gated.
        computed: f64,
    },
}

impl EnvelopeObservation {
    const fn tag(self) -> &'static str {
        match self {
            Self::AgainstReference { .. } => "tolerance",
            Self::AgainstInterval { .. } => "interval",
        }
    }
}

/// Fully evaluated envelope rule retained in an [`EnvelopeVerdict`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnvelopeEvaluation {
    /// Absolute-relative comparison and its derived allowed deviation.
    Tolerance {
        /// Independent reference value.
        reference: f64,
        /// Registered absolute tolerance.
        atol: f64,
        /// Registered relative tolerance.
        rtol: f64,
        /// `atol + rtol * abs(reference)`.
        allowed: f64,
        /// `abs(computed - reference)`.
        deviation: f64,
    },
    /// Inclusive registered interval.
    Interval {
        /// Inclusive lower bound.
        lo: f64,
        /// Inclusive upper bound.
        hi: f64,
    },
}

/// Registry-bound inputs retained before arithmetic evaluation.
///
/// The fields are sealed: callers can inspect an attempt returned in a typed
/// refusal, but cannot manufacture a registry identity through this type.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeAttempt {
    entry_id: &'static str,
    qoi: &'static str,
    unit: &'static str,
    envelope: AcceptanceEnvelope,
    observation: EnvelopeObservation,
    entry_digest: ContentHash,
    registry_digest: ContentHash,
    registry_version: u32,
}

impl EnvelopeAttempt {
    /// Bound registry entry id.
    #[must_use]
    pub const fn entry_id(&self) -> &'static str {
        self.entry_id
    }

    /// Bound QoI name.
    #[must_use]
    pub const fn qoi(&self) -> &'static str {
        self.qoi
    }

    /// Bound QoI unit.
    #[must_use]
    pub const fn unit(&self) -> &'static str {
        self.unit
    }

    /// Exact envelope stored on the registry QoI.
    #[must_use]
    pub const fn envelope(&self) -> AcceptanceEnvelope {
        self.envelope
    }

    /// Exact caller observation, including non-finite bit patterns.
    #[must_use]
    pub const fn observation(&self) -> EnvelopeObservation {
        self.observation
    }

    /// Computed scalar offered by the observation.
    #[must_use]
    pub const fn computed(&self) -> f64 {
        match self.observation {
            EnvelopeObservation::AgainstReference { computed, .. }
            | EnvelopeObservation::AgainstInterval { computed } => computed,
        }
    }

    /// Content identity of the bound registry entry.
    #[must_use]
    pub const fn entry_digest(&self) -> ContentHash {
        self.entry_digest
    }

    /// Content identity of the seeded registry used for lookup.
    #[must_use]
    pub const fn registry_digest(&self) -> ContentHash {
        self.registry_digest
    }

    /// Registry payload version used for lookup.
    #[must_use]
    pub const fn registry_version(&self) -> u32 {
        self.registry_version
    }

    fn push_identity_json(&self, line: &mut String, kind: &str) {
        use fmt::Write as _;

        line.push_str("{\"vvreg\":");
        push_json_str(line, kind);
        let _ = write!(
            line,
            ",\"schema\":{ENVELOPE_VERDICT_SCHEMA_VERSION},\"registry_version\":{}",
            self.registry_version
        );
        line.push_str(",\"registry_digest\":");
        push_json_str(line, &self.registry_digest.to_hex());
        line.push_str(",\"entry\":");
        push_json_str(line, self.entry_id);
        line.push_str(",\"entry_digest\":");
        push_json_str(line, &self.entry_digest.to_hex());
        line.push_str(",\"qoi\":");
        push_json_str(line, self.qoi);
        line.push_str(",\"unit\":");
        push_json_str(line, self.unit);
    }

    fn push_definition_json(&self, line: &mut String) {
        match self.envelope {
            AcceptanceEnvelope::Tolerance { atol, rtol } => {
                line.push_str(",\"envelope_mode\":\"tolerance\",\"atol\":");
                push_f64_bits(line, atol);
                line.push_str(",\"rtol\":");
                push_f64_bits(line, rtol);
            }
            AcceptanceEnvelope::Interval { lo, hi } => {
                line.push_str(",\"envelope_mode\":\"interval\",\"lo\":");
                push_f64_bits(line, lo);
                line.push_str(",\"hi\":");
                push_f64_bits(line, hi);
            }
            AcceptanceEnvelope::Unpinned => {
                line.push_str(",\"envelope_mode\":\"unpinned\"");
            }
        }
        match self.observation {
            EnvelopeObservation::AgainstReference {
                reference,
                computed,
            } => {
                line.push_str(",\"observation_mode\":\"reference\",\"reference\":");
                push_f64_bits(line, reference);
                line.push_str(",\"computed\":");
                push_f64_bits(line, computed);
            }
            EnvelopeObservation::AgainstInterval { computed } => {
                line.push_str(",\"observation_mode\":\"interval\",\"computed\":");
                push_f64_bits(line, computed);
            }
        }
    }

    /// Canonical replay record with exact IEEE-754 bit tokens.
    #[must_use]
    pub fn json_line(&self) -> String {
        let mut line = String::new();
        self.push_identity_json(&mut line, "acceptance-envelope-attempt");
        self.push_definition_json(&mut line);
        line.push('}');
        line
    }
}

/// Deterministic arithmetic verdict for one seeded-registry QoI.
///
/// This record is diagnostic data, not an evidence color or citation receipt.
/// It is sealed and can only be returned by
/// [`Registry::check_acceptance_envelope`]. A failing verdict is retained
/// inside [`EnvelopeGateError::Violation`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeVerdict {
    attempt: EnvelopeAttempt,
    evaluation: EnvelopeEvaluation,
    margin: f64,
    pass: bool,
}

impl EnvelopeVerdict {
    /// Registry-bound inputs and identities used by the evaluation.
    #[must_use]
    pub const fn attempt(&self) -> &EnvelopeAttempt {
        &self.attempt
    }

    /// Fully evaluated tolerance or interval.
    #[must_use]
    pub const fn evaluation(&self) -> EnvelopeEvaluation {
        self.evaluation
    }

    /// Signed admission margin: zero or positive passes, negative refuses.
    #[must_use]
    pub const fn margin(&self) -> f64 {
        self.margin
    }

    /// Whether the candidate is inside the inclusive envelope.
    #[must_use]
    pub const fn passed(&self) -> bool {
        self.pass
    }

    /// Canonical diagnostic JSON line with exact IEEE-754 bit tokens.
    #[must_use]
    pub fn json_line(&self) -> String {
        let mut line = String::new();
        self.attempt
            .push_identity_json(&mut line, "acceptance-envelope-verdict");
        match self.evaluation {
            EnvelopeEvaluation::Tolerance {
                reference,
                atol,
                rtol,
                allowed,
                deviation,
            } => {
                line.push_str(",\"mode\":\"tolerance\",\"reference\":");
                push_f64_bits(&mut line, reference);
                line.push_str(",\"computed\":");
                push_f64_bits(&mut line, self.attempt.computed());
                line.push_str(",\"atol\":");
                push_f64_bits(&mut line, atol);
                line.push_str(",\"rtol\":");
                push_f64_bits(&mut line, rtol);
                line.push_str(",\"allowed\":");
                push_f64_bits(&mut line, allowed);
                line.push_str(",\"deviation\":");
                push_f64_bits(&mut line, deviation);
            }
            EnvelopeEvaluation::Interval { lo, hi } => {
                line.push_str(",\"mode\":\"interval\",\"reference\":null,\"computed\":");
                push_f64_bits(&mut line, self.attempt.computed());
                line.push_str(",\"lo\":");
                push_f64_bits(&mut line, lo);
                line.push_str(",\"hi\":");
                push_f64_bits(&mut line, hi);
            }
        }
        line.push_str(",\"margin\":");
        push_f64_bits(&mut line, self.margin);
        line.push_str(if self.pass {
            ",\"pass\":true}"
        } else {
            ",\"pass\":false}"
        });
        line
    }
}

/// Why an executable acceptance-envelope comparison refused.
#[derive(Debug, Clone, PartialEq)]
pub enum EnvelopeGateError {
    /// The seeded-registry lookup or stored envelope definition refused.
    Registry(CitationRefusal),
    /// The named QoI is absent from the uniquely selected registry entry.
    UnknownQoi {
        /// Registry entry id.
        id: &'static str,
        /// Requested QoI name.
        qoi: String,
    },
    /// The observation mode does not match the stored envelope variant.
    ModeMismatch {
        /// Complete registry-bound attempted evaluation.
        attempt: Box<EnvelopeAttempt>,
        /// Envelope mode required by the registry row.
        expected: &'static str,
        /// Mode the caller supplied.
        got: &'static str,
    },
    /// A reference or computed value was non-finite.
    NonFiniteInput {
        /// Complete registry-bound attempted evaluation.
        attempt: Box<EnvelopeAttempt>,
        /// Non-finite input field.
        field: &'static str,
    },
    /// Finite inputs overflowed a derived comparison quantity.
    ArithmeticOverflow {
        /// Complete registry-bound attempted evaluation.
        attempt: Box<EnvelopeAttempt>,
        /// Derived operation that overflowed.
        operation: &'static str,
    },
    /// The finite comparison completed outside the inclusive envelope.
    Violation {
        /// Complete failing diagnostic record.
        verdict: Box<EnvelopeVerdict>,
    },
}

impl EnvelopeGateError {
    /// Canonical refusal/verdict JSON. Arithmetic refusals retain the exact
    /// registry definition and caller observation needed for replay.
    #[must_use]
    pub fn json_line(&self) -> String {
        let (attempt, outcome, detail_key, detail) = match self {
            Self::ModeMismatch {
                attempt,
                expected,
                got,
            } => {
                let detail = format!("expected {expected}; got {got}");
                (attempt, "mode-mismatch", "detail", detail)
            }
            Self::NonFiniteInput { attempt, field } => {
                (attempt, "non-finite-input", "field", (*field).to_string())
            }
            Self::ArithmeticOverflow { attempt, operation } => (
                attempt,
                "arithmetic-overflow",
                "operation",
                (*operation).to_string(),
            ),
            Self::Violation { verdict } => return verdict.json_line(),
            Self::Registry(refusal) => {
                let mut line =
                    String::from("{\"vvreg\":\"acceptance-envelope-refusal\",\"schema\":");
                use fmt::Write as _;
                let _ = write!(line, "{ENVELOPE_VERDICT_SCHEMA_VERSION}");
                line.push_str(",\"outcome\":\"registry\",\"reason\":");
                push_json_str(&mut line, &refusal.to_string());
                line.push('}');
                return line;
            }
            Self::UnknownQoi { id, qoi } => {
                let mut line =
                    String::from("{\"vvreg\":\"acceptance-envelope-refusal\",\"schema\":");
                use fmt::Write as _;
                let _ = write!(line, "{ENVELOPE_VERDICT_SCHEMA_VERSION}");
                line.push_str(",\"outcome\":\"unknown-qoi\",\"entry\":");
                push_json_str(&mut line, id);
                line.push_str(",\"qoi\":");
                push_json_str(&mut line, qoi);
                line.push('}');
                return line;
            }
        };
        let mut line = String::new();
        attempt.push_identity_json(&mut line, "acceptance-envelope-refusal");
        attempt.push_definition_json(&mut line);
        line.push_str(",\"outcome\":");
        push_json_str(&mut line, outcome);
        line.push(',');
        push_json_str(&mut line, detail_key);
        line.push(':');
        push_json_str(&mut line, &detail);
        line.push('}');
        line
    }
}

impl fmt::Display for EnvelopeGateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.json_line())
    }
}

impl std::error::Error for EnvelopeGateError {}

/// One registry row: an executable (or not-yet-executable) benchmark
/// definition. Incomplete rows may exist in the registry — recording a
/// known target is honest — but they refuse citation until every
/// load-bearing field is pinned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegistryEntry {
    /// Stable registry slug (e.g. `"g1-hertz-sphere-plane"`).
    pub id: &'static str,
    /// Gauntlet tier this entry serves.
    pub tier: RegistryTier,
    /// Family name (e.g. `"TEAM"`, `"Hertz contact"`). Never citable alone.
    pub family: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Exact edition pin.
    pub edition: Edition,
    /// Citation/locator for the defining source.
    pub source: &'static str,
    /// License state governing deck storage.
    pub license: LicenseState,
    /// Pinned input-deck identity.
    pub deck: DeckPin,
    /// Whether the deck's oracle is executable from the deck alone.
    pub oracle: OracleBinding,
    /// Quantities of interest with acceptance envelopes.
    pub qois: &'static [Qoi],
    /// Assumption and no-claim boundary notes.
    pub notes: &'static str,
}

/// Why an entry cannot be cited by a Gauntlet test. Each variant names the
/// first gate that failed in the documented entry-validation or lookup
/// order: id shape/size, QoI-count cap, required fields, edition, license,
/// deck, oracle, QoI presence/uniqueness, and per-QoI envelope.
#[derive(Debug, Clone, PartialEq)]
pub enum CitationRefusal {
    /// Citation was attempted against a caller-built registry. Only the
    /// seeded workspace registry behind [`registry()`] carries authority;
    /// caller-built registries are lint-only.
    UnauthoritativeRegistry,
    /// The requested or declared id exceeds [`MAX_LOOKUP_ID_LEN`] bytes. The
    /// refusal deliberately carries only the length so hostile input is
    /// rejected before it is copied into an error value.
    OversizedLookupId {
        /// The offered length in bytes.
        len: usize,
    },
    /// The requested or declared id is not a lowercase ASCII slug. This refusal
    /// deliberately carries no copy of the hostile input.
    InvalidLookupId,
    /// The requested id is not in the registry.
    UnknownEntry {
        /// The id that was requested.
        id: String,
    },
    /// More than one registry row carries this id; admission fails closed
    /// rather than picking one of the conflicting definitions.
    DuplicateEntry {
        /// The ambiguous id.
        id: String,
    },
    /// A required text field is blank.
    EmptyField {
        /// Entry id (empty string when the id itself is blank).
        id: &'static str,
        /// Which field is blank.
        field: &'static str,
    },
    /// The exact version/edition is not pinned.
    UnpinnedEdition {
        /// Entry id.
        id: &'static str,
    },
    /// The license/terms state is not pinned.
    UnpinnedLicense {
        /// Entry id.
        id: &'static str,
    },
    /// The input deck is not pinned.
    UnpinnedDeck {
        /// Entry id.
        id: &'static str,
    },
    /// The external deck digest is not 64 valid hex chars.
    MalformedDeckDigest {
        /// Entry id.
        id: &'static str,
    },
    /// The oracle identity or executable comparison procedure is not pinned.
    UnpinnedOracle {
        /// Entry id.
        id: &'static str,
    },
    /// The deck delegates a load-bearing derivation to the consumer and no
    /// derivation receipt mechanism binds it yet.
    UnboundOracle {
        /// Entry id.
        id: &'static str,
        /// The delegated obligation.
        obligation: &'static str,
    },
    /// The entry declares no quantities of interest.
    MissingQois {
        /// Entry id.
        id: &'static str,
    },
    /// The entry declares more QoIs than [`MAX_QOIS_PER_ENTRY`].
    TooManyQois {
        /// Entry id.
        id: &'static str,
        /// The declared QoI count.
        count: usize,
    },
    /// Two QoIs on the entry share a name.
    DuplicateQoi {
        /// Entry id.
        id: &'static str,
        /// The duplicated QoI name.
        qoi: &'static str,
    },
    /// A QoI has no acceptance envelope.
    UnpinnedEnvelope {
        /// Entry id.
        id: &'static str,
        /// QoI name.
        qoi: &'static str,
    },
    /// An envelope bound is non-finite, negative-tolerance, or inverted.
    InvalidEnvelope {
        /// Entry id.
        id: &'static str,
        /// QoI name.
        qoi: &'static str,
        /// What is wrong with the bound.
        reason: &'static str,
    },
}

impl fmt::Display for CitationRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnauthoritativeRegistry => f.write_str(
                "caller-built registries are lint-only; cite through the seeded registry()",
            ),
            Self::OversizedLookupId { len } => write!(
                f,
                "registry citation id is {len} bytes; the cap is {MAX_LOOKUP_ID_LEN}"
            ),
            Self::InvalidLookupId => {
                f.write_str("registry citation id must be a non-empty lowercase ASCII slug")
            }
            Self::UnknownEntry { id } => {
                write!(f, "registry has no entry '{id}'")
            }
            Self::DuplicateEntry { id } => write!(
                f,
                "registry has conflicting rows for '{id}'; admission fails closed"
            ),
            Self::EmptyField { id, field } => {
                write!(f, "entry '{id}': required field '{field}' is blank")
            }
            Self::UnpinnedEdition { id } => write!(
                f,
                "entry '{id}': exact version/edition is not pinned; a family or problem number alone is not a deck"
            ),
            Self::UnpinnedLicense { id } => {
                write!(f, "entry '{id}': license/terms are not pinned")
            }
            Self::UnpinnedDeck { id } => {
                write!(f, "entry '{id}': input deck is not pinned")
            }
            Self::MalformedDeckDigest { id } => {
                write!(f, "entry '{id}': external deck digest is not 64 hex chars")
            }
            Self::UnpinnedOracle { id } => {
                write!(
                    f,
                    "entry '{id}': oracle identity or procedure is not pinned"
                )
            }
            Self::UnboundOracle { id, obligation } => write!(
                f,
                "entry '{id}': deck delegates a load-bearing derivation ({obligation}); non-citable until a derivation receipt binds it"
            ),
            Self::MissingQois { id } => {
                write!(f, "entry '{id}': no quantities of interest are declared")
            }
            Self::TooManyQois { id, count } => write!(
                f,
                "entry '{id}': {count} QoIs exceeds the cap of {MAX_QOIS_PER_ENTRY}"
            ),
            Self::DuplicateQoi { id, qoi } => {
                write!(f, "entry '{id}': QoI '{qoi}' is declared more than once")
            }
            Self::UnpinnedEnvelope { id, qoi } => {
                write!(f, "entry '{id}': QoI '{qoi}' has no acceptance envelope")
            }
            Self::InvalidEnvelope { id, qoi, reason } => {
                write!(f, "entry '{id}': QoI '{qoi}' envelope invalid: {reason}")
            }
        }
    }
}

impl std::error::Error for CitationRefusal {}

/// Proof that a specific, fully pinned entry was admitted for citation at a
/// specific registry version. The receipt carries the color caps of the
/// registry rule; it never mints a color itself.
///
/// SEALED: all fields are private and there is no public constructor —
/// the only way to obtain a receipt is [`Registry::cite`] on the seeded
/// registry behind [`registry()`], so holding one proves admission ran
/// against that exact registry state (see `registry_digest`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CitationReceipt {
    entry_id: &'static str,
    tier: RegistryTier,
    edition: &'static str,
    deck_digest: ContentHash,
    entry_digest: ContentHash,
    registry_digest: ContentHash,
    registry_version: u32,
}

impl CitationReceipt {
    /// The admitted entry id.
    #[must_use]
    pub const fn entry_id(&self) -> &'static str {
        self.entry_id
    }

    /// Gauntlet tier of the admitted entry.
    #[must_use]
    pub const fn tier(&self) -> RegistryTier {
        self.tier
    }

    /// The exact pinned edition string.
    #[must_use]
    pub const fn edition(&self) -> &'static str {
        self.edition
    }

    /// The pinned deck identity.
    #[must_use]
    pub const fn deck_digest(&self) -> ContentHash {
        self.deck_digest
    }

    /// The full entry identity (pins the artifact version consumers record).
    #[must_use]
    pub const fn entry_digest(&self) -> ContentHash {
        self.entry_digest
    }

    /// The content identity of the seeded registry this receipt was minted
    /// against (binds the receipt to the exact registry state).
    #[must_use]
    pub const fn registry_digest(&self) -> ContentHash {
        self.registry_digest
    }

    /// Registry payload version the admission ran under.
    #[must_use]
    pub const fn registry_version(&self) -> u32 {
        self.registry_version
    }

    /// RULE: a standards/benchmark calculation earns at most numerical
    /// `Verified` (or standard-conformance) for the exact edition and scope.
    // The cap is deliberately reachable only THROUGH an admitted receipt.
    #[allow(clippy::unused_self)]
    #[must_use]
    pub const fn numerical_claim_cap(&self) -> ColorRank {
        ColorRank::Verified
    }

    /// RULE: the physical prediction cap is `Estimated`, unconditionally,
    /// in this slice. Upgrading to `Validated` requires a TYPED binding of
    /// independent held-out evidence for the named QoI and population — a
    /// bare caller-asserted flag is forgeable and is deliberately not
    /// offered. No color is inherited from a publisher's name.
    // The cap is deliberately reachable only THROUGH an admitted receipt.
    #[allow(clippy::unused_self)]
    #[must_use]
    pub const fn physical_claim_cap(&self) -> ColorRank {
        ColorRank::Estimated
    }
}

/// Validate one row against the citation gates WITHOUT minting authority:
/// no receipt and no color API can be reached through this function. It
/// exists so seed authors and downstream pin work can probe gates.
///
/// Check order (documented so refusal tests are stable): id shape/size,
/// QoI-count cap, blank required fields (family, title, source, notes, QoI
/// names/units), edition, license, deck, oracle binding, QoI presence,
/// duplicate QoI names, per-QoI envelope pin, per-QoI envelope validity.
///
/// # Errors
///
/// Returns the typed [`CitationRefusal`] for the first unpinned, blank,
/// ambiguous, or invalid load-bearing field.
pub fn validate_entry(entry: &RegistryEntry) -> Result<(), CitationRefusal> {
    gate_entry(entry).map(|_pinned| ())
}

/// The full gate chain; on success returns the pinned edition and deck
/// digest. Private: receipts are minted only through [`Registry::cite`]
/// on the authoritative seeded registry.
fn gate_entry(entry: &RegistryEntry) -> Result<(&'static str, ContentHash), CitationRefusal> {
    validate_registry_id(entry.id)?;
    // Bound every later QoI traversal before touching caller-supplied rows.
    if entry.qois.len() > MAX_QOIS_PER_ENTRY {
        return Err(CitationRefusal::TooManyQois {
            id: entry.id,
            count: entry.qois.len(),
        });
    }
    check_text_fields(entry)?;
    let Edition::Exact { version } = entry.edition else {
        return Err(CitationRefusal::UnpinnedEdition { id: entry.id });
    };
    if version.trim().is_empty() {
        return Err(CitationRefusal::EmptyField {
            id: entry.id,
            field: "edition.version",
        });
    }
    check_license(entry)?;
    let deck_digest = check_deck(entry)?;
    match entry.oracle {
        OracleBinding::Unpinned => {
            return Err(CitationRefusal::UnpinnedOracle { id: entry.id });
        }
        OracleBinding::DerivationRequired { obligation } => {
            return Err(CitationRefusal::UnboundOracle {
                id: entry.id,
                obligation,
            });
        }
        OracleBinding::SelfContained => {}
    }
    if entry.qois.is_empty() {
        return Err(CitationRefusal::MissingQois { id: entry.id });
    }
    for (i, qoi) in entry.qois.iter().enumerate() {
        if entry.qois[..i].iter().any(|prior| prior.name == qoi.name) {
            return Err(CitationRefusal::DuplicateQoi {
                id: entry.id,
                qoi: qoi.name,
            });
        }
    }
    for qoi in entry.qois {
        check_envelope(entry.id, qoi)?;
    }
    Ok((version, deck_digest))
}

fn check_text_fields(entry: &RegistryEntry) -> Result<(), CitationRefusal> {
    let checks: [(&str, &'static str); 5] = [
        (entry.id, "id"),
        (entry.family, "family"),
        (entry.title, "title"),
        (entry.source, "source"),
        (entry.notes, "notes"),
    ];
    for (value, field) in checks {
        if value.trim().is_empty() {
            return Err(CitationRefusal::EmptyField {
                id: entry.id,
                field,
            });
        }
    }
    for qoi in entry.qois {
        if qoi.name.trim().is_empty() {
            return Err(CitationRefusal::EmptyField {
                id: entry.id,
                field: "qoi.name",
            });
        }
        if qoi.unit.trim().is_empty() {
            return Err(CitationRefusal::EmptyField {
                id: entry.id,
                field: "qoi.unit",
            });
        }
    }
    Ok(())
}

fn check_license(entry: &RegistryEntry) -> Result<(), CitationRefusal> {
    match entry.license {
        LicenseState::Spdx { id } => {
            if id.trim().is_empty() {
                return Err(CitationRefusal::EmptyField {
                    id: entry.id,
                    field: "license.spdx",
                });
            }
        }
        LicenseState::Restricted { terms, storage } => {
            if terms.trim().is_empty() || storage.trim().is_empty() {
                return Err(CitationRefusal::EmptyField {
                    id: entry.id,
                    field: "license.restricted",
                });
            }
        }
        LicenseState::Unpinned => {
            return Err(CitationRefusal::UnpinnedLicense { id: entry.id });
        }
    }
    Ok(())
}

fn check_deck(entry: &RegistryEntry) -> Result<ContentHash, CitationRefusal> {
    match entry.deck {
        DeckPin::AuthoredSpec { spec } => {
            if spec.trim().is_empty() {
                return Err(CitationRefusal::EmptyField {
                    id: entry.id,
                    field: "deck.spec",
                });
            }
            Ok(hash_domain(DECK_IDENTITY_DOMAIN, spec.as_bytes()))
        }
        DeckPin::External { digest_hex } => ContentHash::from_hex(digest_hex)
            .ok_or(CitationRefusal::MalformedDeckDigest { id: entry.id }),
        DeckPin::Unpinned => Err(CitationRefusal::UnpinnedDeck { id: entry.id }),
    }
}

fn check_envelope(id: &'static str, qoi: &Qoi) -> Result<(), CitationRefusal> {
    match qoi.envelope {
        AcceptanceEnvelope::Tolerance { atol, rtol } => {
            if !atol.is_finite() || !rtol.is_finite() {
                return Err(CitationRefusal::InvalidEnvelope {
                    id,
                    qoi: qoi.name,
                    reason: "non-finite tolerance",
                });
            }
            if atol < 0.0 || rtol < 0.0 {
                return Err(CitationRefusal::InvalidEnvelope {
                    id,
                    qoi: qoi.name,
                    reason: "negative tolerance",
                });
            }
            if atol == 0.0 && rtol == 0.0 {
                return Err(CitationRefusal::InvalidEnvelope {
                    id,
                    qoi: qoi.name,
                    reason: "zero-width tolerance (declare an Interval for exact claims)",
                });
            }
        }
        AcceptanceEnvelope::Interval { lo, hi } => {
            if !lo.is_finite() || !hi.is_finite() {
                return Err(CitationRefusal::InvalidEnvelope {
                    id,
                    qoi: qoi.name,
                    reason: "non-finite bound",
                });
            }
            if lo > hi {
                return Err(CitationRefusal::InvalidEnvelope {
                    id,
                    qoi: qoi.name,
                    reason: "inverted interval",
                });
            }
        }
        AcceptanceEnvelope::Unpinned => {
            return Err(CitationRefusal::UnpinnedEnvelope { id, qoi: qoi.name });
        }
    }
    Ok(())
}

fn evaluate_envelope_attempt(
    attempt: EnvelopeAttempt,
) -> Result<EnvelopeVerdict, EnvelopeGateError> {
    match (attempt.envelope, attempt.observation) {
        (
            AcceptanceEnvelope::Tolerance { atol, rtol },
            EnvelopeObservation::AgainstReference {
                reference,
                computed,
            },
        ) => evaluate_tolerance_attempt(attempt, reference, computed, atol, rtol),
        (
            AcceptanceEnvelope::Interval { lo, hi },
            EnvelopeObservation::AgainstInterval { computed },
        ) => evaluate_interval_attempt(attempt, computed, lo, hi),
        (AcceptanceEnvelope::Tolerance { .. }, observation) => {
            Err(EnvelopeGateError::ModeMismatch {
                attempt: Box::new(attempt),
                expected: "tolerance",
                got: observation.tag(),
            })
        }
        (AcceptanceEnvelope::Interval { .. }, observation) => {
            Err(EnvelopeGateError::ModeMismatch {
                attempt: Box::new(attempt),
                expected: "interval",
                got: observation.tag(),
            })
        }
        (AcceptanceEnvelope::Unpinned, _) => Err(EnvelopeGateError::Registry(
            CitationRefusal::UnpinnedEnvelope {
                id: attempt.entry_id,
                qoi: attempt.qoi,
            },
        )),
    }
}

fn evaluate_tolerance_attempt(
    attempt: EnvelopeAttempt,
    reference: f64,
    computed: f64,
    atol: f64,
    rtol: f64,
) -> Result<EnvelopeVerdict, EnvelopeGateError> {
    if !reference.is_finite() {
        return Err(EnvelopeGateError::NonFiniteInput {
            attempt: Box::new(attempt),
            field: "reference",
        });
    }
    if !computed.is_finite() {
        return Err(EnvelopeGateError::NonFiniteInput {
            attempt: Box::new(attempt),
            field: "computed",
        });
    }
    let relative = rtol * reference.abs();
    if !relative.is_finite() {
        return Err(EnvelopeGateError::ArithmeticOverflow {
            attempt: Box::new(attempt),
            operation: "rtol * abs(reference)",
        });
    }
    let allowed = atol + relative;
    if !allowed.is_finite() {
        return Err(EnvelopeGateError::ArithmeticOverflow {
            attempt: Box::new(attempt),
            operation: "atol + relative tolerance",
        });
    }
    let delta = computed - reference;
    if !delta.is_finite() {
        return Err(EnvelopeGateError::ArithmeticOverflow {
            attempt: Box::new(attempt),
            operation: "computed - reference",
        });
    }
    let deviation = delta.abs();
    let margin = allowed - deviation;
    if !margin.is_finite() {
        return Err(EnvelopeGateError::ArithmeticOverflow {
            attempt: Box::new(attempt),
            operation: "allowed - deviation",
        });
    }
    finish_envelope_verdict(EnvelopeVerdict {
        attempt,
        evaluation: EnvelopeEvaluation::Tolerance {
            reference,
            atol,
            rtol,
            allowed,
            deviation,
        },
        margin,
        pass: deviation <= allowed,
    })
}

fn evaluate_interval_attempt(
    attempt: EnvelopeAttempt,
    computed: f64,
    lo: f64,
    hi: f64,
) -> Result<EnvelopeVerdict, EnvelopeGateError> {
    if !computed.is_finite() {
        return Err(EnvelopeGateError::NonFiniteInput {
            attempt: Box::new(attempt),
            field: "computed",
        });
    }
    let lower_margin = computed - lo;
    let upper_margin = hi - computed;
    if !lower_margin.is_finite() || !upper_margin.is_finite() {
        return Err(EnvelopeGateError::ArithmeticOverflow {
            attempt: Box::new(attempt),
            operation: "interval signed margin",
        });
    }
    finish_envelope_verdict(EnvelopeVerdict {
        attempt,
        evaluation: EnvelopeEvaluation::Interval { lo, hi },
        margin: lower_margin.min(upper_margin),
        pass: computed >= lo && computed <= hi,
    })
}

fn finish_envelope_verdict(verdict: EnvelopeVerdict) -> Result<EnvelopeVerdict, EnvelopeGateError> {
    if verdict.pass {
        Ok(verdict)
    } else {
        Err(EnvelopeGateError::Violation {
            verdict: Box::new(verdict),
        })
    }
}

/// Appendix-D consumption discipline: how a consuming bead has engaged with
/// a registry artifact. All five states are recordable — honesty about
/// `Unread` is the point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsumptionStatus {
    /// Cited but not yet read.
    Unread,
    /// Read; definitions understood.
    Read,
    /// Result re-derived from the source.
    Derived,
    /// Result reproduced executably.
    Reproduced,
    /// Independently falsified/checked through a second route.
    IndependentlyFalsified,
}

impl ConsumptionStatus {
    /// Stable tag for rows and logs.
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self {
            Self::Unread => "unread",
            Self::Read => "read",
            Self::Derived => "derived",
            Self::Reproduced => "reproduced",
            Self::IndependentlyFalsified => "independently_falsified",
        }
    }

    /// All states, in fixed order (for exhaustive row/round-trip tests).
    #[must_use]
    pub const fn all() -> [Self; 5] {
        [
            Self::Unread,
            Self::Read,
            Self::Derived,
            Self::Reproduced,
            Self::IndependentlyFalsified,
        ]
    }
}

/// Maximum accepted consuming-bead id length in bytes (validated BEFORE
/// the string is copied into a record).
pub const MAX_BEAD_ID_LEN: usize = 256;

/// Maximum accepted registry entry or lookup id length in bytes. IDs are
/// also restricted to lowercase ASCII slugs; lookup checks run before a
/// missing id can be copied into [`CitationRefusal::UnknownEntry`].
pub const MAX_LOOKUP_ID_LEN: usize = 256;

fn is_registry_slug(id: &str) -> bool {
    let bytes = id.as_bytes();
    let Some((&first, rest)) = bytes.split_first() else {
        return false;
    };
    first.is_ascii_lowercase()
        && bytes
            .last()
            .is_some_and(|last| last.is_ascii_alphanumeric())
        && rest
            .iter()
            .all(|byte| byte.is_ascii_lowercase() || byte.is_ascii_digit() || *byte == b'-')
        && !bytes.windows(2).any(|pair| pair == b"--")
}

fn validate_registry_id(id: &str) -> Result<(), CitationRefusal> {
    if id.len() > MAX_LOOKUP_ID_LEN {
        return Err(CitationRefusal::OversizedLookupId { len: id.len() });
    }
    if !is_registry_slug(id) {
        return Err(CitationRefusal::InvalidLookupId);
    }
    Ok(())
}

/// Why a consumption record could not be bound.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConsumptionRefusal {
    /// The consuming bead id is blank.
    EmptyBead,
    /// The consuming bead id exceeds [`MAX_BEAD_ID_LEN`] bytes.
    OversizedBead {
        /// The offered length in bytes.
        len: usize,
    },
}

impl fmt::Display for ConsumptionRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyBead => f.write_str("consumption record requires a non-blank bead id"),
            Self::OversizedBead { len } => write!(
                f,
                "consumption record bead id is {len} bytes; the cap is {MAX_BEAD_ID_LEN}"
            ),
        }
    }
}

impl std::error::Error for ConsumptionRefusal {}

/// A consuming bead's recorded engagement with one admitted registry entry,
/// pinned to the exact artifact version (the entry digest).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsumptionRecord {
    /// Consuming bead id.
    pub bead: String,
    /// The cited entry id.
    pub entry_id: &'static str,
    /// The exact entry identity at citation time.
    pub entry_digest: ContentHash,
    /// Registry payload version at citation time.
    pub registry_version: u32,
    /// Recorded engagement state.
    pub status: ConsumptionStatus,
}

impl ConsumptionRecord {
    /// Bind a record to an admitted citation.
    ///
    /// # Errors
    ///
    /// Refuses a blank bead id, and refuses an oversized bead id BEFORE
    /// copying it (validate-before-allocate).
    pub fn bind(
        receipt: &CitationReceipt,
        bead: &str,
        status: ConsumptionStatus,
    ) -> Result<Self, ConsumptionRefusal> {
        if bead.len() > MAX_BEAD_ID_LEN {
            return Err(ConsumptionRefusal::OversizedBead { len: bead.len() });
        }
        if bead.trim().is_empty() {
            return Err(ConsumptionRefusal::EmptyBead);
        }
        Ok(Self {
            bead: bead.to_string(),
            entry_id: receipt.entry_id,
            entry_digest: receipt.entry_digest,
            registry_version: receipt.registry_version,
            status,
        })
    }
}

/// A primary literature/standards reference: anchors definitions and
/// benchmark provenance. Deliberately mints no color and grants no
/// authority — the type has no path to `Color`/`ColorRank`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimaryReference {
    /// Stable 1-based index (matches the bead's seed list).
    pub index: u32,
    /// Stable slug (e.g. `"feec-stability-afw"`).
    pub key: &'static str,
    /// Human citation.
    pub citation: &'static str,
    /// Exact locator (arXiv id, DOI, report number).
    pub locator: &'static str,
    /// What the reference anchors (definitions, not authority).
    pub anchors: &'static str,
    /// What the reference does NOT establish for this workspace.
    pub boundary: &'static str,
}

/// Registry-wide integrity finding (seed-authoring defects, not citation
/// gates).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IntegrityFinding {
    /// Two rows share an id.
    DuplicateEntryId {
        /// The colliding id.
        id: &'static str,
    },
    /// Two references share a key or index.
    DuplicateReference {
        /// The colliding key.
        key: &'static str,
    },
    /// A reference field is blank.
    BlankReferenceField {
        /// Reference key.
        key: &'static str,
        /// Which field is blank.
        field: &'static str,
    },
}

/// The registry lint verdict: which entries are citable, which refuse and
/// why, and any seed-integrity defects.
#[derive(Debug, Clone, PartialEq)]
pub struct RegistryLint {
    /// Ids that admit citation today.
    pub citable: Vec<&'static str>,
    /// Refusals for every non-citable entry.
    pub refused: Vec<CitationRefusal>,
    /// Seed-authoring defects.
    pub integrity: Vec<IntegrityFinding>,
}

/// Whether a registry instance carries citation authority. Private on
/// purpose: only the seeded constructor behind [`registry()`] can set
/// `Seeded`, so a caller-built registry can never mint receipts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RegistryAuthority {
    /// The workspace-seeded registry: `cite` may mint receipts.
    Seeded,
    /// A caller-built registry: lint/serialization only.
    Unauthoritative,
}

/// The benchmark & V&V registry: entries sorted by id plus the
/// primary-reference table sorted by index.
#[derive(Debug, Clone, PartialEq)]
pub struct Registry {
    entries: Vec<RegistryEntry>,
    references: Vec<PrimaryReference>,
    authority: RegistryAuthority,
}

static REGISTRY: LazyLock<Registry> = LazyLock::new(|| {
    let mut seeded = Registry::build(seed::entries(), seed::references());
    seeded.authority = RegistryAuthority::Seeded;
    seeded
});

/// The seeded workspace registry.
#[must_use]
pub fn registry() -> &'static Registry {
    &REGISTRY
}

impl Registry {
    /// Build a caller-owned registry from rows (sorted deterministically
    /// by id/index). Building never panics; integrity defects surface
    /// through [`Registry::lint`]. The result is UNTRUSTED: it supports
    /// lint, lookup, serialization, and digests, but [`Registry::cite`]
    /// refuses to mint receipts from it — only the seeded registry behind
    /// [`registry()`] carries citation authority.
    #[must_use]
    pub fn build(mut entries: Vec<RegistryEntry>, mut references: Vec<PrimaryReference>) -> Self {
        // Total canonical order: id first, then the full content identity
        // as the tie-break, so even conflicting duplicate-id rows land in
        // one input-order-independent arrangement (the digest contract).
        entries.sort_by_cached_key(|e| (e.id, entry_digest(e).0));
        references
            .sort_by_cached_key(|r| (r.index, r.key, r.citation, r.locator, r.anchors, r.boundary));
        Self {
            entries,
            references,
            authority: RegistryAuthority::Unauthoritative,
        }
    }

    /// All entries, sorted by id.
    #[must_use]
    pub fn entries(&self) -> &[RegistryEntry] {
        &self.entries
    }

    /// All primary references, sorted by index.
    #[must_use]
    pub fn references(&self) -> &[PrimaryReference] {
        &self.references
    }

    /// Look up one entry.
    #[must_use]
    pub fn entry(&self, id: &str) -> Option<&RegistryEntry> {
        self.entries.iter().find(|e| e.id == id)
    }

    /// Execute the exact acceptance envelope stored on one uniquely selected
    /// QoI in the seeded workspace registry.
    ///
    /// The returned diagnostic binds both the entry and registry content
    /// digests. Caller-built registries refuse before lookup, and callers name
    /// only the QoI: they cannot substitute a looser envelope or different
    /// unit. This is still an arithmetic gate, not a citation receipt or
    /// evidence color; reference/run provenance remains the consuming lane's
    /// responsibility.
    ///
    /// # Errors
    ///
    /// [`EnvelopeGateError::Registry`] for an unauthoritative registry,
    /// malformed/unknown/duplicate entry, duplicate QoI, or invalid/unpinned
    /// stored envelope; [`EnvelopeGateError::UnknownQoi`] for an absent QoI;
    /// the remaining variants retain the complete bound attempt or verdict.
    pub fn check_acceptance_envelope(
        &self,
        id: &str,
        qoi_name: &str,
        observation: EnvelopeObservation,
    ) -> Result<EnvelopeVerdict, EnvelopeGateError> {
        if self.authority != RegistryAuthority::Seeded {
            return Err(EnvelopeGateError::Registry(
                CitationRefusal::UnauthoritativeRegistry,
            ));
        }
        validate_registry_id(id).map_err(EnvelopeGateError::Registry)?;
        let mut entries = self.entries.iter().filter(|entry| entry.id == id);
        let entry = entries.next().ok_or_else(|| {
            EnvelopeGateError::Registry(CitationRefusal::UnknownEntry { id: id.to_string() })
        })?;
        if entries.next().is_some() {
            return Err(EnvelopeGateError::Registry(
                CitationRefusal::DuplicateEntry { id: id.to_string() },
            ));
        }

        let mut qois = entry.qois.iter().filter(|qoi| qoi.name == qoi_name);
        let qoi = qois.next().ok_or_else(|| EnvelopeGateError::UnknownQoi {
            id: entry.id,
            qoi: qoi_name.to_string(),
        })?;
        if qois.next().is_some() {
            return Err(EnvelopeGateError::Registry(CitationRefusal::DuplicateQoi {
                id: entry.id,
                qoi: qoi.name,
            }));
        }
        if qoi.name.trim().is_empty() {
            return Err(EnvelopeGateError::Registry(CitationRefusal::EmptyField {
                id: entry.id,
                field: "qoi.name",
            }));
        }
        if qoi.unit.trim().is_empty() {
            return Err(EnvelopeGateError::Registry(CitationRefusal::EmptyField {
                id: entry.id,
                field: "qoi.unit",
            }));
        }
        check_envelope(entry.id, qoi).map_err(EnvelopeGateError::Registry)?;

        evaluate_envelope_attempt(EnvelopeAttempt {
            entry_id: entry.id,
            qoi: qoi.name,
            unit: qoi.unit,
            envelope: qoi.envelope,
            observation,
            entry_digest: entry_digest(entry),
            registry_digest: self.digest(),
            registry_version: VVREG_VERSION,
        })
    }

    /// Look up one primary reference by key.
    #[must_use]
    pub fn reference(&self, key: &str) -> Option<&PrimaryReference> {
        self.references.iter().find(|r| r.key == key)
    }

    /// Cite an entry for a Gauntlet test: fail-closed admission. The
    /// receipt binds the seeded registry's content digest, so it names
    /// exactly which registry state admitted the entry.
    ///
    /// # Errors
    ///
    /// [`CitationRefusal::UnauthoritativeRegistry`] on any caller-built
    /// registry (synthetic rows can never mint receipts),
    /// [`CitationRefusal::OversizedLookupId`] or
    /// [`CitationRefusal::InvalidLookupId`] before copying hostile input,
    /// [`CitationRefusal::UnknownEntry`] for a missing id,
    /// [`CitationRefusal::DuplicateEntry`] when the registry holds
    /// conflicting rows for the id (admission never picks one), otherwise
    /// the first failing gate from [`validate_entry`].
    pub fn cite(&self, id: &str) -> Result<CitationReceipt, CitationRefusal> {
        if self.authority != RegistryAuthority::Seeded {
            return Err(CitationRefusal::UnauthoritativeRegistry);
        }
        validate_registry_id(id)?;
        let mut matches = self.entries.iter().filter(|e| e.id == id);
        let entry = matches
            .next()
            .ok_or_else(|| CitationRefusal::UnknownEntry { id: id.to_string() })?;
        if matches.next().is_some() {
            return Err(CitationRefusal::DuplicateEntry { id: id.to_string() });
        }
        let (edition, deck_digest) = gate_entry(entry)?;
        Ok(CitationReceipt {
            entry_id: entry.id,
            tier: entry.tier,
            edition,
            deck_digest,
            entry_digest: entry_digest(entry),
            registry_digest: self.digest(),
            registry_version: VVREG_VERSION,
        })
    }

    /// Lint the whole registry: citability of every row plus seed
    /// integrity (duplicate ids/keys/indices, blank reference fields).
    /// A duplicated id is NEVER citable: every row carrying it is refused
    /// with [`CitationRefusal::DuplicateEntry`] regardless of pin state.
    #[must_use]
    pub fn lint(&self) -> RegistryLint {
        let mut duplicate_ids = Vec::new();
        let mut last_duplicate = None;
        for pair in self.entries.windows(2) {
            if pair[0].id == pair[1].id && last_duplicate != Some(pair[0].id) {
                duplicate_ids.push(pair[0].id);
                last_duplicate = Some(pair[0].id);
            }
        }
        let mut citable = Vec::new();
        let mut refused = Vec::new();
        for entry in &self.entries {
            // duplicate_ids is sorted (built from id-sorted entries).
            if duplicate_ids.binary_search(&entry.id).is_ok() {
                refused.push(CitationRefusal::DuplicateEntry {
                    id: entry.id.to_string(),
                });
                continue;
            }
            match validate_entry(entry) {
                Ok(()) => citable.push(entry.id),
                Err(refusal) => refused.push(refusal),
            }
        }
        let mut integrity = Vec::new();
        for &id in &duplicate_ids {
            integrity.push(IntegrityFinding::DuplicateEntryId { id });
        }
        // Same index (rows are sorted by (index, key), so index duplicates
        // are adjacent) ...
        for pair in self.references.windows(2) {
            if pair[0].index == pair[1].index {
                integrity.push(IntegrityFinding::DuplicateReference { key: pair[1].key });
            }
        }
        // ... and same key ANYWHERE, including at different indices, where
        // sorting by (index, key) does not make the collision adjacent.
        let mut keys: Vec<&'static str> = self.references.iter().map(|r| r.key).collect();
        keys.sort_unstable();
        for pair in keys.windows(2) {
            if pair[0] == pair[1] {
                integrity.push(IntegrityFinding::DuplicateReference { key: pair[1] });
            }
        }
        for reference in &self.references {
            let checks: [(&str, &'static str); 5] = [
                (reference.key, "key"),
                (reference.citation, "citation"),
                (reference.locator, "locator"),
                (reference.anchors, "anchors"),
                (reference.boundary, "boundary"),
            ];
            for (value, field) in checks {
                if value.trim().is_empty() {
                    integrity.push(IntegrityFinding::BlankReferenceField {
                        key: reference.key,
                        field,
                    });
                }
            }
        }
        RegistryLint {
            citable,
            refused,
            integrity,
        }
    }

    /// Canonical serialization: one deterministic row string per entry,
    /// sorted by id. Floats are rendered as IEEE-754 bit patterns
    /// (`0x`-prefixed 16-hex-digit tokens) so the row text is
    /// formatting-independent.
    #[must_use]
    pub fn canonical_rows(&self) -> Vec<String> {
        self.entries.iter().map(canonical_row).collect()
    }

    /// The registry's canonical content identity: a domain-separated
    /// BLAKE3 hash over the length-framed schema version, payload version,
    /// per-entry digests, and the primary-reference table.
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(&REGISTRY_IDENTITY_SCHEMA_VERSION.to_le_bytes());
        payload.extend_from_slice(&VVREG_VERSION.to_le_bytes());
        payload.extend_from_slice(&(self.entries.len() as u64).to_le_bytes());
        for entry in &self.entries {
            payload.extend_from_slice(entry_digest(entry).as_bytes());
        }
        payload.extend_from_slice(&(self.references.len() as u64).to_le_bytes());
        for reference in &self.references {
            payload.extend_from_slice(&reference.index.to_le_bytes());
            frame(&mut payload, reference.key.as_bytes());
            frame(&mut payload, reference.citation.as_bytes());
            frame(&mut payload, reference.locator.as_bytes());
            frame(&mut payload, reference.anchors.as_bytes());
            frame(&mut payload, reference.boundary.as_bytes());
        }
        hash_domain(REGISTRY_IDENTITY_DOMAIN, &payload)
    }
}

/// Length-frame `bytes` into `out` (u64 LE length + bytes).
fn frame(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

/// The canonical identity of one entry: every semantic field is
/// mutation-sensitive. Strings are length-framed; enums carry variant
/// tags; floats contribute their exact IEEE-754 bits (LE).
#[must_use]
pub fn entry_digest(entry: &RegistryEntry) -> ContentHash {
    let mut payload = Vec::new();
    frame(&mut payload, entry.id.as_bytes());
    payload.push(entry.tier.byte());
    match entry.edition {
        Edition::Exact { version } => {
            payload.push(1);
            frame(&mut payload, version.as_bytes());
        }
        Edition::Unpinned => payload.push(0),
    }
    frame(&mut payload, entry.family.as_bytes());
    frame(&mut payload, entry.title.as_bytes());
    frame(&mut payload, entry.source.as_bytes());
    match entry.license {
        LicenseState::Spdx { id } => {
            payload.push(1);
            frame(&mut payload, id.as_bytes());
        }
        LicenseState::Restricted { terms, storage } => {
            payload.push(2);
            frame(&mut payload, terms.as_bytes());
            frame(&mut payload, storage.as_bytes());
        }
        LicenseState::Unpinned => payload.push(0),
    }
    match entry.deck {
        DeckPin::AuthoredSpec { spec } => {
            payload.push(1);
            let digest = hash_domain(DECK_IDENTITY_DOMAIN, spec.as_bytes());
            payload.extend_from_slice(digest.as_bytes());
        }
        // A well-formed external digest is normalized to its raw 32 bytes
        // (one canonical spelling: hex case cannot fork the identity); a
        // malformed one keeps its raw text under a DISTINCT tag so the
        // malformed state itself is mutation-sensitive.
        DeckPin::External { digest_hex } => match ContentHash::from_hex(digest_hex) {
            Some(digest) => {
                payload.push(2);
                payload.extend_from_slice(digest.as_bytes());
            }
            None => {
                payload.push(3);
                frame(&mut payload, digest_hex.as_bytes());
            }
        },
        DeckPin::Unpinned => payload.push(0),
    }
    match entry.oracle {
        OracleBinding::Unpinned => payload.push(0),
        OracleBinding::SelfContained => payload.push(1),
        OracleBinding::DerivationRequired { obligation } => {
            payload.push(2);
            frame(&mut payload, obligation.as_bytes());
        }
    }
    payload.extend_from_slice(&(entry.qois.len() as u64).to_le_bytes());
    for qoi in entry.qois {
        frame(&mut payload, qoi.name.as_bytes());
        frame(&mut payload, qoi.unit.as_bytes());
        match qoi.envelope {
            AcceptanceEnvelope::Tolerance { atol, rtol } => {
                payload.push(1);
                payload.extend_from_slice(&atol.to_bits().to_le_bytes());
                payload.extend_from_slice(&rtol.to_bits().to_le_bytes());
            }
            AcceptanceEnvelope::Interval { lo, hi } => {
                payload.push(2);
                payload.extend_from_slice(&lo.to_bits().to_le_bytes());
                payload.extend_from_slice(&hi.to_bits().to_le_bytes());
            }
            AcceptanceEnvelope::Unpinned => payload.push(0),
        }
    }
    frame(&mut payload, entry.notes.as_bytes());
    hash_domain(ENTRY_IDENTITY_DOMAIN, &payload)
}

fn escape_json_into(out: &mut String, s: &str) {
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                use fmt::Write as _;
                let code = c as u32;
                let _ = write!(out, "\\u{code:04x}");
            }
            c => out.push(c),
        }
    }
}

fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    escape_json_into(out, s);
    out.push('"');
}

fn push_f64_bits(out: &mut String, value: f64) {
    use fmt::Write as _;
    let bits = value.to_bits();
    let _ = write!(out, "\"0x{bits:016x}\"");
}

/// The canonical row string for one entry (deterministic, sorted-field,
/// bit-exact float tokens). This is the serialization the goldens pin.
#[must_use]
pub fn canonical_row(entry: &RegistryEntry) -> String {
    let mut row = String::new();
    row.push_str("{\"id\":");
    push_json_str(&mut row, entry.id);
    row.push_str(",\"tier\":");
    push_json_str(&mut row, entry.tier.tag());
    row.push_str(",\"family\":");
    push_json_str(&mut row, entry.family);
    row.push_str(",\"title\":");
    push_json_str(&mut row, entry.title);
    row.push_str(",\"edition\":");
    match entry.edition {
        Edition::Exact { version } => push_json_str(&mut row, version),
        Edition::Unpinned => row.push_str("null"),
    }
    row.push_str(",\"source\":");
    push_json_str(&mut row, entry.source);
    row.push_str(",\"license\":");
    match entry.license {
        LicenseState::Spdx { id } => {
            row.push_str("{\"spdx\":");
            push_json_str(&mut row, id);
            row.push('}');
        }
        LicenseState::Restricted { terms, storage } => {
            row.push_str("{\"terms\":");
            push_json_str(&mut row, terms);
            row.push_str(",\"storage\":");
            push_json_str(&mut row, storage);
            row.push('}');
        }
        LicenseState::Unpinned => row.push_str("null"),
    }
    row.push_str(",\"deck\":");
    // The row preserves the deck VARIANT and state, mirroring the identity
    // encoding: authored vs external vs malformed-external vs unpinned are
    // four distinct spellings, and valid external hex is normalized
    // through `ContentHash` to one canonical lowercase form.
    match entry.deck {
        // A blank authored spec has no well-formed digest; the row keeps
        // the authored VARIANT visible with a null digest (a state
        // distinct from unpinned and from malformed-external).
        DeckPin::AuthoredSpec { .. } => match entry.deck.digest() {
            Some(digest) => {
                row.push_str("{\"authored\":");
                push_json_str(&mut row, &digest.to_hex());
                row.push('}');
            }
            None => row.push_str("{\"authored\":null}"),
        },
        DeckPin::External { digest_hex } => match ContentHash::from_hex(digest_hex) {
            Some(digest) => {
                row.push_str("{\"external\":");
                push_json_str(&mut row, &digest.to_hex());
                row.push('}');
            }
            None => {
                row.push_str("{\"malformed\":");
                push_json_str(&mut row, digest_hex);
                row.push('}');
            }
        },
        DeckPin::Unpinned => row.push_str("null"),
    }
    row.push_str(",\"oracle\":");
    match entry.oracle {
        OracleBinding::Unpinned => row.push_str("null"),
        OracleBinding::SelfContained => row.push_str("\"self-contained\""),
        OracleBinding::DerivationRequired { obligation } => {
            row.push_str("{\"derivation_required\":");
            push_json_str(&mut row, obligation);
            row.push('}');
        }
    }
    row.push_str(",\"qois\":[");
    for (i, qoi) in entry.qois.iter().enumerate() {
        if i > 0 {
            row.push(',');
        }
        row.push_str("{\"name\":");
        push_json_str(&mut row, qoi.name);
        row.push_str(",\"unit\":");
        push_json_str(&mut row, qoi.unit);
        row.push_str(",\"envelope\":");
        match qoi.envelope {
            AcceptanceEnvelope::Tolerance { atol, rtol } => {
                row.push_str("{\"kind\":\"tolerance\",\"atol\":");
                push_f64_bits(&mut row, atol);
                row.push_str(",\"rtol\":");
                push_f64_bits(&mut row, rtol);
                row.push('}');
            }
            AcceptanceEnvelope::Interval { lo, hi } => {
                row.push_str("{\"kind\":\"interval\",\"lo\":");
                push_f64_bits(&mut row, lo);
                row.push_str(",\"hi\":");
                push_f64_bits(&mut row, hi);
                row.push('}');
            }
            AcceptanceEnvelope::Unpinned => row.push_str("null"),
        }
        row.push('}');
    }
    row.push_str("],\"notes\":");
    push_json_str(&mut row, entry.notes);
    row.push('}');
    row
}
