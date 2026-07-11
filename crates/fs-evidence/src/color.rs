//! The THREE-COLOR epistemic schema (Proposal 3): every ledger quantity
//! carries a COLOR — **verified** (interval-certified numerics),
//! **validated** (anchored to experimental data in a stated REGIME), or
//! **estimated** (cross-model probes, surrogates) — and the composition
//! algebra is CONSERVATIVE BY TYPE: an estimate can never be laundered
//! into a certificate, and validation is a REGIONAL property that
//! auto-demotes the moment execution leaves its regime.
//!
//! Why: a pipeline can be interval-certified end-to-end and still be
//! precisely WRONG, because model-form error does not compose
//! algebraically. The colors keep the boundary between numerical
//! certainty and modeling uncertainty visible in the type system.
//!
//! Layer discipline (bead qmao.1 polish notes): the enum, payloads, and
//! pairwise algebra live HERE in fs-evidence so any layer can color a
//! `Certified<T>` without touching HELM; write-time enforcement lives
//! in fs-ledger over already-colored values.

use crate::{ModelEvidence, NumericalCertificate, NumericalKind, ValidityDomain};
use std::collections::BTreeMap;

/// Maximum byte length of a machine-readable color-provenance identity.
///
/// Composition replaces longer human-readable chains with a domain-separated
/// digest so conservative evidence propagation remains total and bounded.
pub const MAX_COLOR_IDENTITY_BYTES: usize = 256;
/// Semantic version of color composition, identity grammar, and canonical
/// color bytes. Durable derived receipts must bind this value.
pub const COLOR_ALGEBRA_VERSION: u32 = 2;
const _: () = assert!(COLOR_ALGEBRA_VERSION <= u8::MAX as u32);

const COMPOSED_IDENTITY_DOMAIN: &str = "org.frankensim.fs-evidence.composed-identity.v2";
const INVALID_CARD_SET_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-evidence.invalid-card-set-identity.v2";
const DERIVED_IDENTITY_PREFIX: &str = "derived:v2:";

fn is_placeholder_token(value: &str) -> bool {
    [
        "-",
        "?",
        "n/a",
        "na",
        "none",
        "not run",
        "pending",
        "placeholder",
        "tbd",
        "todo",
        "unknown",
    ]
    .iter()
    .any(|placeholder| value.eq_ignore_ascii_case(placeholder))
}

/// Structural reason a machine-readable color identity is unusable.
///
/// Derived identities are allowed here because derived color nodes must pass
/// structural replay. Leaf admission additionally calls
/// [`color_leaf_identity_reason`] to require lineage for the reserved derived
/// namespace.
#[must_use]
pub fn color_identity_reason(value: &str) -> Option<&'static str> {
    if value.len() > MAX_COLOR_IDENTITY_BYTES {
        return Some("too-long");
    }
    let trimmed = value.trim();
    if trimmed.is_empty() {
        Some("blank")
    } else if trimmed != value {
        Some("surrounding-whitespace")
    } else if value.chars().any(char::is_control) {
        Some("control-character")
    } else if !value.bytes().all(|byte| {
        byte.is_ascii_alphanumeric()
            || matches!(byte, b'-' | b'_' | b'.' | b'/' | b':' | b'@' | b'+' | b'=')
    }) {
        Some("invalid-character")
    } else if is_placeholder_token(value) {
        Some("placeholder")
    } else {
        None
    }
}

/// Structural reason an identity cannot enter as a provenance leaf.
#[must_use]
pub fn color_leaf_identity_reason(value: &str) -> Option<&'static str> {
    color_identity_reason(value).or_else(|| {
        value
            .starts_with("derived:")
            .then_some("derived-identity-requires-lineage")
    })
}

fn compact_identity(hash_label: &str, value: &str) -> String {
    let label_hash = fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, hash_label.as_bytes());
    let value_hash = fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, value.as_bytes());
    let mut identity = [0_u8; 64];
    identity[..32].copy_from_slice(label_hash.as_bytes());
    identity[32..].copy_from_slice(value_hash.as_bytes());
    format!(
        "{DERIVED_IDENTITY_PREFIX}{hash_label}:{}",
        fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, &identity)
    )
}

fn bounded_pair_identity(hash_label: &str, left: &str, separator: &str, right: &str) -> String {
    if color_identity_reason(left).is_none() && color_identity_reason(right).is_none() {
        let readable = format!(
            "{DERIVED_IDENTITY_PREFIX}{hash_label}:{}:{left}{separator}{}:{right}",
            left.len(),
            right.len(),
        );
        if readable.len() <= MAX_COLOR_IDENTITY_BYTES {
            return readable;
        }
    }

    let label_hash = fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, hash_label.as_bytes());
    let left_hash = fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, left.as_bytes());
    let right_hash = fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, right.as_bytes());
    let mut identity = [0_u8; 96];
    identity[..32].copy_from_slice(label_hash.as_bytes());
    identity[32..64].copy_from_slice(left_hash.as_bytes());
    identity[64..].copy_from_slice(right_hash.as_bytes());
    format!(
        "{DERIVED_IDENTITY_PREFIX}{hash_label}:{}",
        fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, &identity)
    )
}

fn composed_estimator_identity(left: &str, right: &str) -> String {
    bounded_pair_identity("composed", left, "+", right)
}

fn verified_estimator_identity(estimator: &str) -> String {
    bounded_pair_identity("composed-verified", estimator, "+", "verified")
}

fn validated_estimator_identity(estimator: &str, dataset: &str) -> String {
    bounded_pair_identity("composed-validated", estimator, "+", dataset)
}

fn composed_dataset_identity(left: &str, right: &str) -> String {
    bounded_pair_identity("datasets", left, "+", right)
}

fn verified_dataset_identity(dataset: &str) -> String {
    bounded_pair_identity("datasets-verified", dataset, "+", "verified")
}

fn disjoint_regime_identity(left: &str, right: &str) -> String {
    bounded_pair_identity("disjoint-regimes", left, "+", right)
}

/// Bounded canonical estimator identity for a validated-color regime exit.
///
/// Exposed so replay engines reconstruct the exact same derived color rather
/// than duplicating the identity grammar.
#[must_use]
pub fn demotion_estimator_identity(dataset: &str, axis: &str) -> String {
    bounded_pair_identity("regime-exit", dataset, "@", axis)
}

fn compact_invalid_card_set_identity(cards: &[&String]) -> String {
    fn update_field(hasher: &mut fs_blake3::Blake3, bytes: &[u8]) {
        let len = u64::try_from(bytes.len()).expect("a Rust allocation length fits u64");
        hasher.update(&len.to_le_bytes());
        hasher.update(bytes);
    }

    let mut hasher = fs_blake3::Blake3::new();
    update_field(&mut hasher, INVALID_CARD_SET_IDENTITY_DOMAIN.as_bytes());
    let count = u64::try_from(cards.len()).expect("a Rust allocation length fits u64");
    hasher.update(&count.to_le_bytes());
    for card in cards {
        let reason = color_leaf_identity_reason(card).unwrap_or("valid");
        update_field(&mut hasher, reason.as_bytes());
        update_field(&mut hasher, card.as_bytes());
    }
    let set_hash = fs_blake3::hash_domain(
        INVALID_CARD_SET_IDENTITY_DOMAIN,
        hasher.finalize().as_bytes(),
    );
    compact_identity("invalid-card-set", &set_hash.to_string())
}

fn joined_card_identity(cards: &[String]) -> Result<String, String> {
    let mut canonical = cards.iter().collect::<Vec<_>>();
    canonical.sort_unstable();
    canonical.dedup();
    if canonical
        .iter()
        .any(|card| color_leaf_identity_reason(card).is_some())
    {
        return Err(compact_invalid_card_set_identity(&canonical));
    }
    let mut cards = canonical.into_iter();
    let Some(first) = cards.next() else {
        return Ok("uncarded-numerics".to_string());
    };
    cards.try_fold(first.clone(), |identity, card| {
        Ok(bounded_pair_identity("model-cards", &identity, "+", card))
    })
}

/// The epistemic color, with its color-specific payload.
#[derive(Debug, Clone, PartialEq)]
pub enum Color {
    /// Interval-certified numerics: the payload IS the bound.
    Verified {
        /// Certified lower bound.
        lo: f64,
        /// Certified upper bound.
        hi: f64,
    },
    /// Anchored to experimental data INSIDE a stated regime.
    Validated {
        /// The region of feature space where the anchoring holds
        /// (Reynolds range, strain range, …).
        regime: ValidityDomain,
        /// Identity of the anchoring dataset.
        dataset: String,
    },
    /// Cross-model discrepancy probes, surrogates, heuristics.
    Estimated {
        /// The estimator's identity.
        estimator: String,
        /// The estimator's own dispersion (∞ = no spread claim).
        dispersion: f64,
    },
}

/// A malformed color payload. Claim-strength exceptions cannot authorize any
/// of these structural defects.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColorPayloadError {
    /// A machine-readable identity is malformed.
    InvalidIdentity {
        /// Payload field (`dataset`, `axis`, or `estimator`).
        field: &'static str,
        /// Offending value.
        value: String,
        /// Stable reason.
        reason: &'static str,
    },
    /// A Verified interval contains NaN or is inverted. Ordered infinities are
    /// valid, possibly vacuous enclosures.
    InvalidVerifiedInterval {
        /// Stable reason.
        reason: &'static str,
    },
    /// A Validated regime is missing or unusable.
    InvalidValidatedRegime {
        /// Offending axis, or empty when no axis was declared.
        axis: String,
        /// Stable reason.
        reason: &'static str,
    },
    /// Estimated dispersion is NaN or negative.
    InvalidEstimatedDispersion {
        /// Stable reason.
        reason: &'static str,
    },
}

impl core::fmt::Display for ColorPayloadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidIdentity {
                field,
                value,
                reason,
            } => write!(f, "invalid {field} identity {value:?}: {reason}"),
            Self::InvalidVerifiedInterval { reason } => {
                write!(f, "invalid Verified interval: {reason}")
            }
            Self::InvalidValidatedRegime { axis, reason } if axis.is_empty() => {
                write!(f, "invalid Validated regime: {reason}")
            }
            Self::InvalidValidatedRegime { axis, reason } => {
                write!(f, "invalid Validated regime axis {axis:?}: {reason}")
            }
            Self::InvalidEstimatedDispersion { reason } => {
                write!(f, "invalid Estimated dispersion: {reason}")
            }
        }
    }
}

impl std::error::Error for ColorPayloadError {}

/// Validate the shared structural invariant for one color payload.
///
/// Ordered infinite Verified endpoints are allowed: `[-inf,+inf]` is a sound
/// but vacuous enclosure. NaN and inverted intervals are never valid.
///
/// # Errors
/// Returns the exact malformed field and stable reason.
pub fn validate_color_payload(color: &Color) -> Result<(), ColorPayloadError> {
    match color {
        Color::Verified { lo, hi } => {
            if lo.is_nan() || hi.is_nan() {
                Err(ColorPayloadError::InvalidVerifiedInterval {
                    reason: "bounds contain NaN",
                })
            } else if lo > hi {
                Err(ColorPayloadError::InvalidVerifiedInterval {
                    reason: "lower bound exceeds upper bound",
                })
            } else {
                Ok(())
            }
        }
        Color::Validated { regime, dataset } => {
            if let Some(reason) = color_identity_reason(dataset) {
                return Err(ColorPayloadError::InvalidIdentity {
                    field: "dataset",
                    value: dataset.clone(),
                    reason,
                });
            }
            if regime.bounds().is_empty() {
                return Err(ColorPayloadError::InvalidValidatedRegime {
                    axis: String::new(),
                    reason: "at least one bounded axis is required",
                });
            }
            for (axis, (lo, hi)) in regime.bounds() {
                if let Some(reason) = color_identity_reason(axis) {
                    return Err(ColorPayloadError::InvalidIdentity {
                        field: "axis",
                        value: axis.clone(),
                        reason,
                    });
                }
                if !lo.is_finite() || !hi.is_finite() {
                    return Err(ColorPayloadError::InvalidValidatedRegime {
                        axis: axis.clone(),
                        reason: "bounds must be finite",
                    });
                }
                if lo > hi {
                    return Err(ColorPayloadError::InvalidValidatedRegime {
                        axis: axis.clone(),
                        reason: "lower bound exceeds upper bound",
                    });
                }
            }
            Ok(())
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            if let Some(reason) = color_identity_reason(estimator) {
                return Err(ColorPayloadError::InvalidIdentity {
                    field: "estimator",
                    value: estimator.clone(),
                    reason,
                });
            }
            if dispersion.is_nan() || *dispersion < 0.0 {
                Err(ColorPayloadError::InvalidEstimatedDispersion {
                    reason: "value is NaN or negative",
                })
            } else {
                Ok(())
            }
        }
    }
}

fn compact_malformed_color_identity(label: &str, color: &Color) -> String {
    let mut hasher = fs_blake3::Blake3::new();
    hasher.update(fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, label.as_bytes()).as_bytes());
    match color {
        Color::Verified { lo, hi } => {
            hasher.update(&[0]);
            hasher.update(&lo.to_bits().to_le_bytes());
            hasher.update(&hi.to_bits().to_le_bytes());
        }
        Color::Validated { regime, dataset } => {
            hasher.update(&[1]);
            hasher.update(
                fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, dataset.as_bytes()).as_bytes(),
            );
            hasher.update(
                &u64::try_from(regime.bounds().len())
                    .expect("a Rust allocation length fits u64")
                    .to_le_bytes(),
            );
            for (axis, (lo, hi)) in regime.bounds() {
                hasher.update(
                    fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, axis.as_bytes()).as_bytes(),
                );
                hasher.update(&lo.to_bits().to_le_bytes());
                hasher.update(&hi.to_bits().to_le_bytes());
            }
        }
        Color::Estimated {
            estimator,
            dispersion,
        } => {
            hasher.update(&[2]);
            hasher.update(
                fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, estimator.as_bytes()).as_bytes(),
            );
            hasher.update(&dispersion.to_bits().to_le_bytes());
        }
    }
    let payload_hash = hasher.finalize();
    format!(
        "{DERIVED_IDENTITY_PREFIX}{label}:{}",
        fs_blake3::hash_domain(COMPOSED_IDENTITY_DOMAIN, payload_hash.as_bytes())
    )
}

/// Lattice rank (higher = stronger claim). The composition result can
/// never OUTRANK the weakest operand — that is the no-laundering law.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ColorRank {
    /// Estimated.
    Estimated,
    /// Validated.
    Validated,
    /// Verified.
    Verified,
}

impl Color {
    /// The lattice rank.
    #[must_use]
    pub fn rank(&self) -> ColorRank {
        match self {
            Color::Verified { .. } => ColorRank::Verified,
            Color::Validated { .. } => ColorRank::Validated,
            Color::Estimated { .. } => ColorRank::Estimated,
        }
    }

    /// Stable name for ledger rows.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Color::Verified { .. } => "verified",
            Color::Validated { .. } => "validated",
            Color::Estimated { .. } => "estimated",
        }
    }

    /// Canonical JSON payload for ledger rows.
    #[must_use]
    pub fn payload_json(&self) -> String {
        match self {
            Color::Verified { lo, hi } => {
                format!("{{\"interval\":[{},{}]}}", json_f64(*lo), json_f64(*hi))
            }
            Color::Validated { regime, dataset } => {
                use core::fmt::Write as _;
                let mut axes = String::new();
                for (k, (lo, hi)) in regime.bounds() {
                    if !axes.is_empty() {
                        axes.push(',');
                    }
                    let _ = write!(
                        axes,
                        "{}:[{},{}]",
                        json_string(k),
                        json_f64(*lo),
                        json_f64(*hi)
                    );
                }
                format!(
                    "{{\"dataset\":{},\"regime\":{{{axes}}}}}",
                    json_string(dataset)
                )
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => format!(
                "{{\"estimator\":{},\"dispersion\":{}}}",
                json_string(estimator),
                json_f64(*dispersion)
            ),
        }
    }

    /// Versioned canonical identity bytes for hashing and authorization.
    ///
    /// Encoding v2 starts with a version byte and a variant byte, then uses
    /// u64-LE length prefixes for every variable-width field. Floating-point
    /// values are encoded as their exact IEEE-754 bit patterns in little-endian
    /// order, so values that render identically in [`Self::payload_json`] (and
    /// `0.0` versus `-0.0`) remain distinct. Validated-regime axes follow the
    /// deterministic [`BTreeMap`] order exposed by [`ValidityDomain::bounds`].
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        const VERSION: u8 = COLOR_ALGEBRA_VERSION as u8;
        const VERIFIED: u8 = 0;
        const VALIDATED: u8 = 1;
        const ESTIMATED: u8 = 2;

        let mut out = Vec::new();
        out.push(VERSION);
        match self {
            Color::Verified { lo, hi } => {
                out.push(VERIFIED);
                push_canonical_field(&mut out, &lo.to_bits().to_le_bytes());
                push_canonical_field(&mut out, &hi.to_bits().to_le_bytes());
            }
            Color::Validated { regime, dataset } => {
                out.push(VALIDATED);
                push_canonical_field(&mut out, dataset.as_bytes());
                push_canonical_len(&mut out, regime.bounds().len());
                for (axis, (lo, hi)) in regime.bounds() {
                    push_canonical_field(&mut out, axis.as_bytes());
                    push_canonical_field(&mut out, &lo.to_bits().to_le_bytes());
                    push_canonical_field(&mut out, &hi.to_bits().to_le_bytes());
                }
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => {
                out.push(ESTIMATED);
                push_canonical_field(&mut out, estimator.as_bytes());
                push_canonical_field(&mut out, &dispersion.to_bits().to_le_bytes());
            }
        }
        out
    }
}

fn push_canonical_len(out: &mut Vec<u8>, len: usize) {
    let len = u64::try_from(len).expect("a Rust allocation length always fits in u64");
    out.extend_from_slice(&len.to_le_bytes());
}

fn push_canonical_field(out: &mut Vec<u8>, bytes: &[u8]) {
    push_canonical_len(out, bytes.len());
    out.extend_from_slice(bytes);
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        format!("\"non-finite:{value}\"")
    }
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// How verified intervals combine under the ledger operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalOp {
    /// Sum of quantities: bounds add.
    Add,
    /// Product: the four-corner hull.
    Mul,
    /// Anything else: the conservative hull of both bounds.
    Hull,
}

/// A demotion record (validated → estimated on regime exit).
#[derive(Debug, Clone, PartialEq)]
pub struct Demotion {
    /// The dataset whose regime was exited.
    pub dataset: String,
    /// The axis that left its range.
    pub axis: String,
    /// The offending state value.
    pub value: f64,
}

/// Teaching errors for the color algebra.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorError {
    /// A `Verified` claim was attempted from a non-enclosure
    /// certificate — the laundering refusal.
    LaunderingRefused {
        /// What the certificate actually was.
        actual: &'static str,
    },
    /// A regime axis referenced by the state is missing from the tag.
    IncompleteState {
        /// The axis the regime declares but the state omits.
        axis: String,
    },
}

impl core::fmt::Display for ColorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ColorError::LaunderingRefused { actual } => write!(
                f,
                "refusing to mark this quantity `verified`: its certificate is \
                 `{actual}`, not an enclosure — estimates cannot be laundered into \
                 certificates (attach a signed waiver to the LEDGER write if a \
                 human accepts responsibility, and the waiver will travel in \
                 provenance)"
            ),
            ColorError::IncompleteState { axis } => write!(
                f,
                "the execution state does not report regime axis `{axis}`; a \
                 validated claim cannot be checked against an unreported axis \
                 (report it, or the value demotes)"
            ),
        }
    }
}

impl std::error::Error for ColorError {}

/// Color a quantity `verified` FROM an interval-certified numerical
/// certificate — the only door in. Anything weaker refuses.
///
/// # Errors
/// [`ColorError::LaunderingRefused`] for estimate/no-claim kinds.
pub fn verified_from(cert: &NumericalCertificate) -> Result<Color, ColorError> {
    match cert.kind {
        NumericalKind::Exact | NumericalKind::Enclosure => {
            // A `Verified` color asserts `[lo, hi]` is a valid enclosure; a NaN
            // or inverted (`lo > hi`) interval encloses NOTHING, so fail closed
            // rather than mint a false certificate from garbage bounds — e.g.
            // `exact(NaN)` or a hand-built inverted cert (bead wa8i E4).
            // Infinite bounds are a valid (if loose) enclosure and pass.
            if cert.lo.is_nan() || cert.hi.is_nan() || cert.lo > cert.hi {
                return Err(ColorError::LaunderingRefused {
                    actual: "non-enclosure (NaN or inverted bounds)",
                });
            }
            Ok(Color::Verified {
                lo: cert.lo,
                hi: cert.hi,
            })
        }
        NumericalKind::Estimate => Err(ColorError::LaunderingRefused { actual: "estimate" }),
        NumericalKind::NoClaim => Err(ColorError::LaunderingRefused { actual: "no-claim" }),
    }
}

/// Derive the honest color of an existing [`crate::Evidence`] receipt.
/// Model-free enclosure-grade numerics can become `Verified`; any plain
/// [`ModelEvidence`] remains `Estimated` because model cards, simulation
/// discrepancy pairs, and validity boxes do not authenticate an experimental
/// anchor. A future typed anchored-source receipt is the admission path for
/// `Validated`.
#[must_use]
pub fn color_of(numerical: &NumericalCertificate, model: &ModelEvidence) -> Color {
    let model_is_absent = model.cards.is_empty()
        && model.assumptions.is_empty()
        && model.validity.bounds().is_empty()
        && model.discrepancy_rel.to_bits() == 0.0_f64.to_bits()
        && model.in_domain;
    let cards = joined_card_identity(&model.cards);
    let numerical_is_usable = match numerical.kind {
        NumericalKind::Exact | NumericalKind::Enclosure => verified_from(numerical).is_ok(),
        NumericalKind::Estimate => {
            numerical.lo.is_finite() && numerical.hi.is_finite() && numerical.lo <= numerical.hi
        }
        NumericalKind::NoClaim => false,
    };
    let model_is_usable = (model_is_absent || !model.cards.is_empty())
        && model.in_domain
        && !model.discrepancy_rel.is_nan()
        && model.discrepancy_rel >= 0.0
        && model.validity.bounds().iter().all(|(axis, (lo, hi))| {
            color_identity_reason(axis).is_none() && lo.is_finite() && hi.is_finite() && lo <= hi
        });

    if model_is_absent {
        // Route through the guarded door so uncarded numerics with NaN/inverted
        // bounds fall through to Estimated instead of minting a false Verified
        // (bead wa8i E4). Valid Exact/Enclosure bounds are unchanged.
        if let Ok(verified) = verified_from(numerical) {
            return verified;
        }
    }
    let estimator = match cards {
        Ok(identity) => identity,
        Err(invalid_identity) => {
            return Color::Estimated {
                estimator: invalid_identity,
                dispersion: f64::INFINITY,
            };
        }
    };
    if !numerical_is_usable || !model_is_usable {
        return Color::Estimated {
            estimator,
            dispersion: f64::INFINITY,
        };
    }
    let numerical_dispersion = if numerical.kind == NumericalKind::Estimate {
        let reference = numerical.lo / 2.0 + numerical.hi / 2.0;
        numerical.rel_half_width(reference)
    } else {
        0.0
    };
    Color::Estimated {
        estimator,
        dispersion: model.discrepancy_rel + numerical_dispersion,
    }
}

fn normalize_composition_input(color: &Color) -> Color {
    let Err(error) = validate_color_payload(color) else {
        return color.clone();
    };
    let label = match error {
        ColorPayloadError::InvalidIdentity { .. } => "invalid-color-identity",
        ColorPayloadError::InvalidVerifiedInterval { .. } => "invalid-verified-interval",
        ColorPayloadError::InvalidValidatedRegime { .. } => "invalid-validity-regime",
        ColorPayloadError::InvalidEstimatedDispersion { .. } => "invalid-estimated-dispersion",
    };
    Color::Estimated {
        estimator: compact_malformed_color_identity(label, color),
        dispersion: f64::INFINITY,
    }
}

/// Outward-round an arithmetic interval endpoint pair (one ulp each way) so a
/// composed enclosure stays a TRUE enclosure under round-to-nearest — otherwise
/// the lower bound could round UP (or the upper bound DOWN) and EXCLUDE the
/// exact result, laundering the very precision a `Verified` color certifies. A
/// NaN endpoint (∞ + −∞) degrades to the whole real line, the only sound
/// enclosure. Mirrors `NumericalCertificate::compose` in `lib.rs`.
fn outward_round(lo: f64, hi: f64) -> (f64, f64) {
    if lo.is_nan() || hi.is_nan() {
        (f64::NEG_INFINITY, f64::INFINITY)
    } else {
        (lo.next_down(), hi.next_up())
    }
}

/// The TOTAL, conservative pairwise composition: the result's rank is
/// the MINIMUM of the operands' ranks (no laundering), verified bounds
/// combine per `op` (outward-rounded to preserve enclosure), validated
/// regimes INTERSECT (both anchors must hold), and estimated dispersions
/// add (conservative).
#[must_use]
pub fn compose(a: &Color, b: &Color, op: IntervalOp) -> Color {
    let a = normalize_composition_input(a);
    let b = normalize_composition_input(b);
    match (&a, &b) {
        (Color::Verified { lo: a0, hi: a1 }, Color::Verified { lo: b0, hi: b1 }) => {
            let (lo, hi) = match op {
                IntervalOp::Add => outward_round(a0 + b0, a1 + b1),
                IntervalOp::Mul => {
                    let c = [a0 * b0, a0 * b1, a1 * b0, a1 * b1];
                    if c.iter().any(|x| x.is_nan()) {
                        // A NaN corner (0×∞) is indeterminate AND the min/max
                        // fold below would SILENTLY DROP it (f64::min/max ignore
                        // NaN), reporting a bogus tight interval. The whole real
                        // line is the only sound enclosure.
                        (f64::NEG_INFINITY, f64::INFINITY)
                    } else {
                        outward_round(
                            c.iter().copied().fold(f64::INFINITY, f64::min),
                            c.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                        )
                    }
                }
                // Hull is exact endpoint selection (no arithmetic → no rounding).
                IntervalOp::Hull => (a0.min(*b0), a1.max(*b1)),
            };
            Color::Verified { lo, hi }
        }
        // Validated ⊕ verified stays validated (the weaker anchor);
        // validated ⊕ validated intersects regimes.
        (Color::Validated { regime, dataset }, Color::Verified { .. })
        | (Color::Verified { .. }, Color::Validated { regime, dataset }) => Color::Validated {
            regime: regime.clone(),
            dataset: verified_dataset_identity(dataset),
        },
        (
            Color::Validated {
                regime: r1,
                dataset: d1,
            },
            Color::Validated {
                regime: r2,
                dataset: d2,
            },
        ) => {
            let regime = intersect_domains(r1, r2);
            if regime.is_empty() {
                // Mutually unsatisfiable regimes: there is NO state where both
                // anchors hold, so the composition cannot honestly stay
                // Validated. There is also no defensible spread claim, so the
                // demotion carries infinite dispersion just like a regime exit.
                Color::Estimated {
                    estimator: disjoint_regime_identity(d1, d2),
                    dispersion: f64::INFINITY,
                }
            } else {
                Color::Validated {
                    regime,
                    dataset: composed_dataset_identity(d1, d2),
                }
            }
        }
        // Anything ⊕ estimated → estimated. No exceptions here; the
        // waiver door lives at the LEDGER, in provenance.
        (
            Color::Estimated {
                estimator,
                dispersion,
            },
            other,
        )
        | (
            other,
            Color::Estimated {
                estimator,
                dispersion,
            },
        ) => {
            let (other_disp, composed_identity) = match other {
                Color::Estimated {
                    estimator: e2,
                    dispersion: v,
                } => (*v, composed_estimator_identity(estimator, e2)),
                Color::Verified { .. } => (0.0, verified_estimator_identity(estimator)),
                Color::Validated { dataset, .. } => {
                    (0.0, validated_estimator_identity(estimator, dataset))
                }
            };
            Color::Estimated {
                estimator: composed_identity,
                dispersion: dispersion + other_disp,
            }
        }
    }
}

/// Intersection of two validity domains (axis-wise; an axis present in
/// either constrains the result — both anchors must hold). Delegates to
/// [`ValidityDomain::intersect`], which PRESERVES emptiness (`lo > hi`) on a
/// disjoint axis so callers can detect an unsatisfiable regime via
/// [`ValidityDomain::is_empty`]. (A previous `hi.max(lo)` clamp here collapsed a
/// disjoint intersection into a phantom single point, silently claiming
/// validity at a state where neither anchor holds.)
#[must_use]
pub fn intersect_domains(a: &ValidityDomain, b: &ValidityDomain) -> ValidityDomain {
    a.intersect(b)
}

/// Check a validated color against the CURRENT execution state:
/// inside the regime → unchanged; outside (or unreported axis) →
/// AUTOMATIC DEMOTION to estimated, with the flag returned. Verified
/// and estimated colors pass through untouched.
#[must_use]
pub fn check_regime(color: &Color, state: &BTreeMap<String, f64>) -> (Color, Option<Demotion>) {
    let Some(reason) = regime_demotion(color, state) else {
        return (color.clone(), None);
    };
    (
        Color::Estimated {
            estimator: demotion_estimator_identity(&reason.dataset, &reason.axis),
            dispersion: f64::INFINITY,
        },
        Some(reason),
    )
}

/// Determine whether a validated color exits its regime without cloning the
/// color payload. This is the borrowed preflight used by admission layers that
/// must bound a multi-parent fold before constructing the derived regime.
#[must_use]
pub fn regime_demotion(color: &Color, state: &BTreeMap<String, f64>) -> Option<Demotion> {
    let Color::Validated { regime, dataset } = color else {
        return None;
    };
    if regime.bounds().is_empty() {
        return Some(Demotion {
            dataset: dataset.clone(),
            axis: "<undeclared-regime>".to_string(),
            value: f64::NAN,
        });
    }
    for (axis, &(lo, hi)) in regime.bounds() {
        let Some(&v) = state.get(axis) else {
            return Some(Demotion {
                dataset: dataset.clone(),
                axis: axis.clone(),
                value: f64::NAN,
            });
        };
        if !lo.is_finite() || !hi.is_finite() || lo > hi || !v.is_finite() || v < lo || v > hi {
            return Some(Demotion {
                dataset: dataset.clone(),
                axis: axis.clone(),
                value: v,
            });
        }
    }
    None
}
