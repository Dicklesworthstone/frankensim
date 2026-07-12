//! EXPLANATION OBJECTS (addendum Proposal B, bead knh1.5; [F] — behind
//! the `explanation-objects` feature): when a QoI moves, return not a
//! NUMBER but a CAUSAL DECOMPOSITION along named physical channels,
//! each term with a bound, as a first-class integrity-checked artifact whose
//! aggregate evidence color states its claim strength. Built-in engines bind
//! exact inputs into derivation digests; without a separately retained root
//! and evidence resolver those digests do not claim external authority. The
//! difference between confabulation and understanding is
//! whether the system can CHECK the story — so the explanation is an
//! OBJECT the system checks, and the natural-language rendering on top
//! is explicitly NON-AUTHORITATIVE.
//!
//! Three attribution engines feed one tree:
//! 1. ADJOINT attribution — for elliptic compliance the bilinear trick
//!    gives the EXACT identity `J₁ − J₀ = −∫ Δa·u₀′·u₁′`, so channel
//!    masks decompose ΔJ exactly (quadrature-level bounds).
//! 2. PROVENANCE attribution — which EDIT moved the number, by exact
//!    telescoping over replayed ledger states.
//! 3. PHYSICAL decomposition — the far-field drag FLAGSHIP: induced
//!    drag via the Trefftz-plane wake integral on a lifting-line
//!    fixture (reconciling with the analytic `C_L²/(π·AR)`), a
//!    viscous strip channel, and the wave channel DECLARED zero in the
//!    subsonic regime rather than silently omitted.
//!
//! THE HONESTY GATE (a permanent runtime invariant): if the
//! unattributed residual exceeds its threshold the system REFUSES to
//! explain rather than smearing the residual across plausible
//! channels. A partial explanation with a declared gap beats a
//! complete story with a hidden one.

use fs_evidence::{COLOR_ALGEBRA_VERSION, Color, IntervalOp, compose, validate_color_payload};
use std::collections::BTreeSet;

/// Canonical explanation-node fingerprint semantics.
pub const EXPLANATION_FINGERPRINT_VERSION: u32 = 2;
const EXPLANATION_FINGERPRINT_DOMAIN: &str = "org.frankensim.fs-adjoint.explanation-node.v2";
/// Canonical finalized-explanation receipt semantics.
pub const EXPLANATION_RECEIPT_VERSION: u32 = 1;
const EXPLANATION_RECEIPT_DOMAIN: &str = "org.frankensim.fs-adjoint.explanation-receipt.v1";
const UNRETAINED_DERIVATION_DOMAIN: &str =
    "org.frankensim.fs-adjoint.explanation-unretained-derivation.v1";
const ADJOINT_DERIVATION_DOMAIN: &str =
    "org.frankensim.fs-adjoint.explanation-adjoint-derivation.v1";
const PROVENANCE_DERIVATION_DOMAIN: &str =
    "org.frankensim.fs-adjoint.explanation-provenance-derivation.v1";
const LIFTING_LINE_DERIVATION_DOMAIN: &str =
    "org.frankensim.fs-adjoint.explanation-lifting-line-derivation.v1";
const DRAG_DERIVATION_DOMAIN: &str = "org.frankensim.fs-adjoint.explanation-drag-derivation.v1";
const MAX_CHANNEL_BYTES: usize = 256;
const MAX_EVIDENCE_BYTES: usize = 512;
const MAX_EVIDENCE_LINKS: usize = 64;
const MAX_EXPLANATION_NODES: usize = 1_024;
const MAX_LIFTING_LINE_STATIONS: usize = 4_096;
/// Maximum interior nodes admitted by the fixture-scale elliptic explanation
/// engine before it allocates solver or channel-mask storage.
pub const MAX_ELLIPTIC_INTERIOR_NODES: usize = 65_536;
const MAX_DECLARED_SUBSONIC_MACH: f64 = 0.8;

/// Deterministic refusal from an explanation constructor or engine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExplanationError {
    /// A bounded collection or fixture dimension is outside its contract.
    InvalidCount {
        /// Count being checked.
        field: &'static str,
        /// Caller-provided count.
        value: usize,
        /// Smallest accepted count.
        min: usize,
        /// Largest accepted count.
        max: usize,
    },
    /// A text identity is blank, untrimmed, oversized, or control-bearing.
    InvalidText {
        /// Identity field.
        field: &'static str,
        /// Element index for a collection field.
        index: Option<usize>,
        /// Stable reason.
        reason: &'static str,
    },
    /// A numerical input or derived value is unusable.
    InvalidNumber {
        /// Numerical field or operation.
        field: &'static str,
        /// Element/index location when applicable.
        index: Option<usize>,
        /// Stable reason.
        reason: &'static str,
    },
    /// Related vectors have incompatible dimensions.
    LengthMismatch {
        /// Vector or relation being checked.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Actual length.
        actual: usize,
    },
    /// An index would be outside its admitted collection.
    InvalidIndex {
        /// Indexed collection.
        field: &'static str,
        /// Rejected index.
        index: usize,
        /// Exclusive upper bound.
        upper_bound: usize,
    },
    /// A uniqueness invariant was violated.
    DuplicateIdentity {
        /// Duplicated identity class.
        field: &'static str,
    },
    /// Two channel masks claim the same finite-element index.
    OverlappingChannelElement {
        /// Repeated element index.
        element: usize,
    },
    /// Adjacent provenance edits do not telescope bit-exactly.
    DisconnectedHistory {
        /// Index of the later edit whose input does not match.
        edit_index: usize,
    },
    /// A color payload cannot support an explanation node.
    InvalidColor {
        /// Stable reason without echoing unbounded caller payloads.
        reason: &'static str,
    },
    /// Private authority metadata or a retained fingerprint is inconsistent.
    IntegrityMismatch {
        /// Integrity field.
        field: &'static str,
        /// Node index when applicable.
        index: Option<usize>,
    },
    /// The tridiagonal fixture solve encountered an unusable pivot.
    SingularPivot {
        /// Pivot row.
        index: usize,
    },
}

impl core::fmt::Display for ExplanationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidCount {
                field,
                value,
                min,
                max,
            } => write!(f, "{field} must be in {min}..={max}, got {value}"),
            Self::InvalidText {
                field,
                index,
                reason,
            }
            | Self::InvalidNumber {
                field,
                index,
                reason,
            } => match index {
                Some(index) => write!(f, "invalid {field}[{index}]: {reason}"),
                None => write!(f, "invalid {field}: {reason}"),
            },
            Self::LengthMismatch {
                field,
                expected,
                actual,
            } => write!(f, "{field} length must be {expected}, got {actual}"),
            Self::InvalidIndex {
                field,
                index,
                upper_bound,
            } => write!(f, "{field} index {index} is outside 0..{upper_bound}"),
            Self::DuplicateIdentity { field } => write!(f, "duplicate {field}"),
            Self::OverlappingChannelElement { element } => {
                write!(f, "channel element {element} is claimed more than once")
            }
            Self::DisconnectedHistory { edit_index } => write!(
                f,
                "provenance edit {edit_index} does not start from the preceding edit's output"
            ),
            Self::InvalidColor { reason } => write!(f, "invalid explanation color: {reason}"),
            Self::IntegrityMismatch { field, index } => match index {
                Some(index) => write!(f, "explanation {field} mismatch at node {index}"),
                None => write!(f, "explanation {field} mismatch"),
            },
            Self::SingularPivot { index } => {
                write!(
                    f,
                    "elliptic tridiagonal pivot {index} is non-finite or non-positive"
                )
            }
        }
    }
}

impl std::error::Error for ExplanationError {}

fn push_len(out: &mut Vec<u8>, len: usize) {
    let len = u64::try_from(len).expect("a Rust allocation length fits u64");
    out.extend_from_slice(&len.to_le_bytes());
}

fn push_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    push_len(out, bytes.len());
    out.extend_from_slice(bytes);
}

fn push_str(out: &mut Vec<u8>, value: &str) {
    push_bytes(out, value.as_bytes());
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    out.extend_from_slice(&value.to_bits().to_le_bytes());
}

fn push_usize(out: &mut Vec<u8>, value: usize) {
    let value = u64::try_from(value).expect("a Rust allocation index fits u64");
    out.extend_from_slice(&value.to_le_bytes());
}

fn bounded_text_is_valid(value: &str, max_bytes: usize) -> bool {
    !value.is_empty()
        && value.trim() == value
        && value.len() <= max_bytes
        && !value.chars().any(char::is_control)
}

fn evidence_is_valid(evidence: &[String]) -> bool {
    !evidence.is_empty()
        && evidence.len() <= MAX_EVIDENCE_LINKS
        && evidence
            .iter()
            .all(|item| bounded_text_is_valid(item, MAX_EVIDENCE_BYTES))
}

fn digest_is_valid(digest: &str) -> bool {
    digest.len() == 64
        && digest
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn derivation_digest(domain: &str, payload: &[u8]) -> String {
    fs_blake3::hash_domain(domain, payload).to_string()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExplanationNodeOrigin {
    BuiltIn,
    Unretained,
}

impl ExplanationNodeOrigin {
    fn tag(self) -> u8 {
        match self {
            Self::BuiltIn => 1,
            Self::Unretained => 2,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ExplanationNodeAuthority {
    origin: ExplanationNodeOrigin,
    derivation_digest: String,
    batch_digest: String,
    batch_index: usize,
    batch_size: usize,
}

impl ExplanationNodeAuthority {
    fn unretained(derivation_digest: String) -> Self {
        Self {
            origin: ExplanationNodeOrigin::Unretained,
            batch_digest: derivation_digest.clone(),
            derivation_digest,
            batch_index: 0,
            batch_size: 1,
        }
    }

    fn built_in(
        derivation_digest: String,
        batch_digest: String,
        batch_index: usize,
        batch_size: usize,
    ) -> Self {
        Self {
            origin: ExplanationNodeOrigin::BuiltIn,
            derivation_digest,
            batch_digest,
            batch_index,
            batch_size,
        }
    }

    fn is_valid(&self) -> bool {
        digest_is_valid(&self.derivation_digest)
            && digest_is_valid(&self.batch_digest)
            && (1..=MAX_EXPLANATION_NODES).contains(&self.batch_size)
            && self.batch_index < self.batch_size
    }
}

fn validate_node_payload(
    channel: &str,
    contribution: f64,
    bound: f64,
    color: &Color,
    evidence: &[String],
) -> Result<(), ExplanationError> {
    if !bounded_text_is_valid(channel, MAX_CHANNEL_BYTES) {
        return Err(ExplanationError::InvalidText {
            field: "channel",
            index: None,
            reason: "must be bounded, trimmed, non-empty, and control-free",
        });
    }
    if !contribution.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "contribution",
            index: None,
            reason: "must be finite",
        });
    }
    if !bound.is_finite() || bound < 0.0 {
        return Err(ExplanationError::InvalidNumber {
            field: "bound",
            index: None,
            reason: "must be finite and non-negative",
        });
    }
    if !(1..=MAX_EVIDENCE_LINKS).contains(&evidence.len()) {
        return Err(ExplanationError::InvalidCount {
            field: "evidence links",
            value: evidence.len(),
            min: 1,
            max: MAX_EVIDENCE_LINKS,
        });
    }
    for (index, item) in evidence.iter().enumerate() {
        if !bounded_text_is_valid(item, MAX_EVIDENCE_BYTES) {
            return Err(ExplanationError::InvalidText {
                field: "evidence link",
                index: Some(index),
                reason: "must be bounded, trimmed, non-empty, and control-free",
            });
        }
    }
    if validate_color_payload(color).is_err() {
        return Err(ExplanationError::InvalidColor {
            reason: "shared color payload is structurally malformed",
        });
    }
    if matches!(color, Color::Validated { .. }) {
        return Err(ExplanationError::InvalidColor {
            reason: "Validated requires a retained regime-membership witness",
        });
    }
    if let Color::Verified { lo, hi } = color {
        let covered_lo = contribution - bound;
        let covered_hi = contribution + bound;
        if !covered_lo.is_finite() || !covered_hi.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "verified contribution envelope",
                index: None,
                reason: "contribution plus/minus bound must remain finite",
            });
        }
        if !lo.is_finite() || !hi.is_finite() || lo > hi || *lo > covered_lo || *hi < covered_hi {
            return Err(ExplanationError::InvalidColor {
                reason: "Verified interval must cover contribution plus/minus bound",
            });
        }
    }
    Ok(())
}

fn node_payload_is_valid(
    channel: &str,
    contribution: f64,
    bound: f64,
    color: &Color,
    evidence: &[String],
) -> bool {
    bounded_text_is_valid(channel, MAX_CHANNEL_BYTES)
        && contribution.is_finite()
        && bound.is_finite()
        && bound >= 0.0
        && evidence_is_valid(evidence)
        && validate_color_payload(color).is_ok()
        && match color {
            Color::Verified { lo, hi } => {
                lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && *lo <= contribution - bound
                    && *hi >= contribution + bound
            }
            Color::Estimated { .. } => true,
            Color::Validated { .. } => false,
        }
}

fn node_fingerprint(
    channel: &str,
    contribution: f64,
    bound: f64,
    color: &Color,
    evidence: &[String],
    authority: &ExplanationNodeAuthority,
) -> String {
    let mut canon = Vec::new();
    canon.extend_from_slice(&EXPLANATION_FINGERPRINT_VERSION.to_le_bytes());
    canon.extend_from_slice(&COLOR_ALGEBRA_VERSION.to_le_bytes());
    canon.push(authority.origin.tag());
    push_str(&mut canon, channel);
    push_f64(&mut canon, contribution);
    push_f64(&mut canon, bound);
    push_bytes(&mut canon, &color.canonical_bytes());
    push_len(&mut canon, evidence.len());
    for item in evidence {
        push_str(&mut canon, item);
    }
    push_str(&mut canon, &authority.derivation_digest);
    push_str(&mut canon, &authority.batch_digest);
    push_usize(&mut canon, authority.batch_index);
    push_usize(&mut canon, authority.batch_size);
    fs_blake3::hash_domain(EXPLANATION_FINGERPRINT_DOMAIN, &canon).to_string()
}

/// One attribution node: a named channel's contribution with its
/// bound, evidence links, a derivation digest, and a payload fingerprint.
#[derive(Debug, Clone, PartialEq)]
pub struct ExplanationNode {
    /// The named physical/provenance channel.
    channel: String,
    /// Signed contribution to the observed ΔQoI.
    contribution: f64,
    /// Declared half-width on the contribution. Only a built-in-origin
    /// `Verified` node's bound contributes to certified residual coverage.
    bound: f64,
    /// The evidence color of this term.
    color: Color,
    /// Evidence-link identities claimed to back the term.
    evidence: Vec<String>,
    authority: ExplanationNodeAuthority,
    /// Deterministic payload-integrity fingerprint. This is not an authority
    /// signature and becomes externally meaningful only when a trusted owner
    /// retains the expected digest and derivation inputs.
    fingerprint: String,
    /// Version of the canonical fingerprint semantics.
    fingerprint_version: u32,
}

impl ExplanationNode {
    /// Build a caller-owned node with an unretained derivation digest derived
    /// from the payload. This detects later mutation; it does not authenticate
    /// the caller's evidence claims.
    ///
    /// # Errors
    /// Returns [`ExplanationError`] for malformed text, numbers, colors, or
    /// evidence links.
    pub fn new(
        channel: &str,
        contribution: f64,
        bound: f64,
        color: Color,
        evidence: Vec<String>,
    ) -> Result<ExplanationNode, ExplanationError> {
        validate_node_payload(channel, contribution, bound, &color, &evidence)?;
        let mut derivation = Vec::new();
        push_str(&mut derivation, channel);
        push_f64(&mut derivation, contribution);
        push_f64(&mut derivation, bound);
        push_bytes(&mut derivation, &color.canonical_bytes());
        push_len(&mut derivation, evidence.len());
        for item in &evidence {
            push_str(&mut derivation, item);
        }
        let derivation_digest = derivation_digest(UNRETAINED_DERIVATION_DOMAIN, &derivation);
        Self::new_with_authority(
            channel,
            contribution,
            bound,
            color,
            evidence,
            ExplanationNodeAuthority::unretained(derivation_digest),
        )
    }

    fn new_with_authority(
        channel: &str,
        contribution: f64,
        bound: f64,
        color: Color,
        evidence: Vec<String>,
        authority: ExplanationNodeAuthority,
    ) -> Result<ExplanationNode, ExplanationError> {
        validate_node_payload(channel, contribution, bound, &color, &evidence)?;
        if !authority.is_valid() {
            return Err(ExplanationError::IntegrityMismatch {
                field: "authority digest/batch coordinates",
                index: None,
            });
        }
        let fingerprint =
            node_fingerprint(channel, contribution, bound, &color, &evidence, &authority);
        Ok(ExplanationNode {
            channel: channel.to_string(),
            contribution,
            bound,
            color,
            evidence,
            authority,
            fingerprint,
            fingerprint_version: EXPLANATION_FINGERPRINT_VERSION,
        })
    }

    /// Named attribution channel.
    #[must_use]
    pub fn channel(&self) -> &str {
        &self.channel
    }

    /// Signed nominal contribution.
    #[must_use]
    pub fn contribution(&self) -> f64 {
        self.contribution
    }

    /// Declared contribution half-width.
    #[must_use]
    pub fn bound(&self) -> f64 {
        self.bound
    }

    /// Node evidence color.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// Claimed evidence-link identities.
    #[must_use]
    pub fn evidence(&self) -> &[String] {
        &self.evidence
    }

    /// Domain-separated exact-input or unretained-payload digest.
    #[must_use]
    pub fn derivation_digest(&self) -> &str {
        &self.authority.derivation_digest
    }

    /// Payload-integrity fingerprint.
    #[must_use]
    pub fn fingerprint(&self) -> &str {
        &self.fingerprint
    }

    /// Fingerprint schema version.
    #[must_use]
    pub fn fingerprint_version(&self) -> u32 {
        self.fingerprint_version
    }

    fn aggregate_color(&self) -> Color {
        match (self.authority.origin, &self.color) {
            (ExplanationNodeOrigin::BuiltIn, color)
            | (ExplanationNodeOrigin::Unretained, color @ Color::Estimated { .. }) => color.clone(),
            (ExplanationNodeOrigin::Unretained, Color::Verified { .. }) => Color::Estimated {
                estimator: "unretained-verified-explanation-node".to_string(),
                dispersion: f64::INFINITY,
            },
            (ExplanationNodeOrigin::Unretained, Color::Validated { .. }) => Color::Estimated {
                estimator: "unretained-validated-explanation-node".to_string(),
                dispersion: f64::INFINITY,
            },
        }
    }

    fn verifies(&self) -> bool {
        self.fingerprint_version == EXPLANATION_FINGERPRINT_VERSION
            && node_payload_is_valid(
                &self.channel,
                self.contribution,
                self.bound,
                &self.color,
                &self.evidence,
            )
            && self.authority.is_valid()
            && self.fingerprint
                == node_fingerprint(
                    &self.channel,
                    self.contribution,
                    self.bound,
                    &self.color,
                    &self.evidence,
                    &self.authority,
                )
    }
}

#[derive(Debug, Clone, Copy)]
struct SumEnvelope {
    rounded: f64,
    lo: f64,
    hi: f64,
}

impl SumEnvelope {
    fn roundoff(self) -> f64 {
        (self.rounded - self.lo)
            .abs()
            .max((self.hi - self.rounded).abs())
    }
}

fn enclosed_sum(values: &[f64]) -> SumEnvelope {
    let Some((&first, rest)) = values.split_first() else {
        return SumEnvelope {
            rounded: 0.0,
            lo: 0.0,
            hi: 0.0,
        };
    };
    let mut envelope = SumEnvelope {
        rounded: first,
        lo: first,
        hi: first,
    };
    for &value in rest {
        envelope.rounded += value;
        envelope.lo = (envelope.lo + value).next_down();
        envelope.hi = (envelope.hi + value).next_up();
    }
    envelope
}

fn aggregate_color(nodes: &[ExplanationNode]) -> Color {
    let mut colors = nodes.iter().map(ExplanationNode::aggregate_color);
    let Some(first) = colors.next() else {
        return Color::Estimated {
            estimator: "empty-explanation".to_string(),
            dispersion: f64::INFINITY,
        };
    };
    colors.fold(first, |aggregate, color| {
        compose(&aggregate, &color, IntervalOp::Add)
    })
}

fn certified_coverage(nodes: &[ExplanationNode], aggregation_roundoff: f64) -> f64 {
    let verified_bounds = nodes
        .iter()
        .filter_map(|node| match (node.authority.origin, &node.color) {
            (ExplanationNodeOrigin::BuiltIn, Color::Verified { .. }) => Some(node.bound),
            (ExplanationNodeOrigin::BuiltIn, Color::Validated { .. } | Color::Estimated { .. })
            | (ExplanationNodeOrigin::Unretained, _) => None,
        })
        .collect::<Vec<_>>();
    let verified_upper = enclosed_sum(&verified_bounds).hi;
    if verified_upper == 0.0 {
        aggregation_roundoff
    } else if aggregation_roundoff == 0.0 {
        verified_upper
    } else {
        (verified_upper + aggregation_roundoff).next_up()
    }
}

fn nodes_are_unique(nodes: &[ExplanationNode]) -> bool {
    let mut fingerprints = BTreeSet::new();
    let mut derivations = BTreeSet::new();
    nodes.iter().all(|node| {
        fingerprints.insert(node.fingerprint.as_str())
            && derivations.insert(node.authority.derivation_digest.as_str())
    })
}

fn built_in_batch_is_coherent(nodes: &[ExplanationNode]) -> bool {
    let built_in = nodes
        .iter()
        .filter(|node| node.authority.origin == ExplanationNodeOrigin::BuiltIn)
        .collect::<Vec<_>>();
    let Some(first) = built_in.first() else {
        return true;
    };
    if built_in.len() != first.authority.batch_size {
        return false;
    }
    let mut indices = BTreeSet::new();
    built_in.iter().all(|node| {
        node.authority.batch_digest == first.authority.batch_digest
            && node.authority.batch_size == first.authority.batch_size
            && indices.insert(node.authority.batch_index)
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExplanationVariant {
    Explained,
    Refused,
}

impl ExplanationVariant {
    fn tag(self) -> u8 {
        match self {
            Self::Explained => 1,
            Self::Refused => 2,
        }
    }
}

/// Replay and integrity metadata for a finalized explanation.
#[derive(Debug, Clone, PartialEq)]
pub struct ExplanationReceipt {
    /// Receipt schema version.
    pub version: u32,
    /// Color-algebra version used for aggregate composition.
    pub color_algebra_version: u32,
    /// Caller-requested residual threshold.
    pub requested_threshold: f64,
    /// Sum of built-in `Verified` bounds plus explicit aggregation roundoff.
    pub certified_coverage: f64,
    /// Stricter of requested threshold and certified coverage.
    pub effective_limit: f64,
    /// Outward-rounded error envelope for summing node contributions.
    pub aggregation_roundoff: f64,
    /// Conservative Add-composition of node colors.
    pub aggregate_color: Color,
    /// Domain-separated payload-integrity root, not an authority signature.
    pub root: String,
}

fn explanation_root(
    variant: ExplanationVariant,
    nodes: &[ExplanationNode],
    observed: f64,
    residual: f64,
    receipt: &ExplanationReceipt,
) -> String {
    let mut canon = Vec::new();
    canon.extend_from_slice(&receipt.version.to_le_bytes());
    canon.push(variant.tag());
    push_len(&mut canon, nodes.len());
    for node in nodes {
        canon.extend_from_slice(&node.fingerprint_version.to_le_bytes());
        push_str(&mut canon, &node.authority.derivation_digest);
        push_str(&mut canon, &node.authority.batch_digest);
        push_usize(&mut canon, node.authority.batch_index);
        push_usize(&mut canon, node.authority.batch_size);
        push_str(&mut canon, &node.fingerprint);
    }
    push_f64(&mut canon, observed);
    push_f64(&mut canon, residual);
    push_f64(&mut canon, receipt.requested_threshold);
    push_f64(&mut canon, receipt.certified_coverage);
    push_f64(&mut canon, receipt.effective_limit);
    push_f64(&mut canon, receipt.aggregation_roundoff);
    canon.extend_from_slice(&receipt.color_algebra_version.to_le_bytes());
    push_bytes(&mut canon, &receipt.aggregate_color.canonical_bytes());
    fs_blake3::hash_domain(EXPLANATION_RECEIPT_DOMAIN, &canon).to_string()
}

fn build_receipt(
    variant: ExplanationVariant,
    nodes: &[ExplanationNode],
    observed: f64,
    residual: f64,
    mut receipt: ExplanationReceipt,
) -> ExplanationReceipt {
    receipt.root = explanation_root(variant, nodes, observed, residual, &receipt);
    receipt
}

fn receipt_payload(
    requested_threshold: f64,
    certified_coverage: f64,
    aggregation_roundoff: f64,
    aggregate_color: Color,
) -> ExplanationReceipt {
    ExplanationReceipt {
        version: EXPLANATION_RECEIPT_VERSION,
        color_algebra_version: COLOR_ALGEBRA_VERSION,
        requested_threshold,
        certified_coverage,
        effective_limit: requested_threshold.min(certified_coverage),
        aggregation_roundoff,
        aggregate_color,
        root: String::new(),
    }
}

fn same_f64(left: f64, right: f64) -> bool {
    left.to_bits() == right.to_bits()
}

fn replay_values_match(left: f64, right: f64) -> bool {
    const SIGNLESS_BITS: u64 = 0x7fff_ffff_ffff_ffff;
    same_f64(left, right)
        || (left.to_bits() & SIGNLESS_BITS == 0 && right.to_bits() & SIGNLESS_BITS == 0)
}

fn common_structure_is_valid(
    variant: ExplanationVariant,
    nodes: &[ExplanationNode],
    observed: f64,
    residual: f64,
    receipt: &ExplanationReceipt,
) -> bool {
    if nodes.is_empty()
        || nodes.len() > MAX_EXPLANATION_NODES
        || !nodes_are_unique(nodes)
        || !built_in_batch_is_coherent(nodes)
        || !observed.is_finite()
        || !residual.is_finite()
        || receipt.version != EXPLANATION_RECEIPT_VERSION
        || receipt.color_algebra_version != COLOR_ALGEBRA_VERSION
        || !receipt.requested_threshold.is_finite()
        || receipt.requested_threshold < 0.0
        || !receipt.certified_coverage.is_finite()
        || receipt.certified_coverage < 0.0
        || !receipt.effective_limit.is_finite()
        || receipt.effective_limit < 0.0
        || !receipt.aggregation_roundoff.is_finite()
        || receipt.aggregation_roundoff < 0.0
        || validate_color_payload(&receipt.aggregate_color).is_err()
        || !digest_is_valid(&receipt.root)
        || !nodes.iter().all(ExplanationNode::verifies)
    {
        return false;
    }

    let contributions = nodes
        .iter()
        .map(|node| node.contribution)
        .collect::<Vec<_>>();
    let sum = enclosed_sum(&contributions);
    if !sum.rounded.is_finite() || !sum.lo.is_finite() || !sum.hi.is_finite() {
        return false;
    }
    let expected_roundoff = sum.roundoff();
    let expected_residual = observed - sum.rounded;
    let expected_coverage = certified_coverage(nodes, expected_roundoff);
    let expected_color = aggregate_color(nodes);
    let expected_limit = receipt.requested_threshold.min(expected_coverage);

    same_f64(receipt.aggregation_roundoff, expected_roundoff)
        && same_f64(receipt.certified_coverage, expected_coverage)
        && same_f64(receipt.effective_limit, expected_limit)
        && same_f64(residual, expected_residual)
        && receipt.aggregate_color.canonical_bytes() == expected_color.canonical_bytes()
        && receipt.root == explanation_root(variant, nodes, observed, residual, receipt)
}

/// A finalized explanation, or the refusal that keeps it honest.
#[derive(Debug, Clone, PartialEq)]
pub enum Explanation {
    /// The tree reconciles within certified coverage.
    Explained {
        /// Channel nodes.
        nodes: Vec<ExplanationNode>,
        /// The observed ΔQoI being explained.
        observed: f64,
        /// The declared unattributed residual.
        residual: f64,
        /// Versioned replay and integrity receipt.
        receipt: ExplanationReceipt,
    },
    /// The residual exceeded the effective limit: no explanation is issued.
    Refused {
        /// The partial non-authoritative attribution.
        partial: Vec<ExplanationNode>,
        /// The observed ΔQoI that could not be explained.
        observed: f64,
        /// The unattributed residual that triggered refusal.
        residual: f64,
        /// Versioned replay and integrity receipt.
        receipt: ExplanationReceipt,
    },
}

impl Explanation {
    /// THE PERMANENT INVARIANT (the Proposal-B kill criterion): the
    /// certified channel bounds must COVER the unattributed residual.
    /// An engine failing this on any case is lying and ships nothing.
    ///
    /// The earlier form checked `attributed + residual == observed`,
    /// which `finalize` makes true BY CONSTRUCTION (residual is
    /// defined as observed − attributed) — a kill criterion that can
    /// never fire is not a kill criterion (bead 9sf6 F3). Coverage is
    /// the non-vacuous claim: bad channel math leaves a residual the
    /// bounds cannot absorb, and THIS fires.
    #[must_use]
    pub fn reconciles(&self) -> bool {
        match self {
            Explanation::Explained {
                nodes,
                observed,
                residual,
                receipt,
            } => {
                common_structure_is_valid(
                    ExplanationVariant::Explained,
                    nodes,
                    *observed,
                    *residual,
                    receipt,
                ) && residual.abs() <= receipt.effective_limit
            }
            Explanation::Refused { .. } => false,
        }
    }

    /// Validate payloads, node and receipt digests, admission arithmetic, and
    /// success/refusal semantics without treating a refusal as reconciled.
    #[must_use]
    pub fn is_structurally_valid(&self) -> bool {
        match self {
            Explanation::Explained {
                nodes,
                observed,
                residual,
                receipt,
            } => {
                common_structure_is_valid(
                    ExplanationVariant::Explained,
                    nodes,
                    *observed,
                    *residual,
                    receipt,
                ) && residual.abs() <= receipt.effective_limit
            }
            Explanation::Refused {
                partial,
                observed,
                residual,
                receipt,
            } => {
                common_structure_is_valid(
                    ExplanationVariant::Refused,
                    partial,
                    *observed,
                    *residual,
                    receipt,
                ) && residual.abs() > receipt.effective_limit
            }
        }
    }

    /// The versioned receipt retained by either outcome.
    #[must_use]
    pub fn receipt(&self) -> &ExplanationReceipt {
        match self {
            Explanation::Explained { receipt, .. } | Explanation::Refused { receipt, .. } => {
                receipt
            }
        }
    }

    /// NON-AUTHORITATIVE natural-language rendering. The TREE is the
    /// artifact; this string is for humans skimming and says so.
    #[must_use]
    pub fn render_narrative(&self) -> String {
        use std::fmt::Write as _;
        let mut out =
            String::from("NON-AUTHORITATIVE RENDERING (the explanation tree is the artifact):\n");
        match self {
            Explanation::Explained {
                nodes,
                observed,
                residual,
                receipt,
            } => {
                let _ = writeln!(out, "observed change {observed:+.6e}");
                for n in nodes {
                    let _ = writeln!(
                        out,
                        "  {} contributed {:+.6e} (± {:.1e})",
                        n.channel, n.contribution, n.bound
                    );
                }
                let _ = writeln!(out, "  unattributed residual {residual:+.6e}");
                let _ = writeln!(out, "  aggregate color {}", receipt.aggregate_color.name());
            }
            Explanation::Refused {
                residual, receipt, ..
            } => {
                let _ = writeln!(
                    out,
                    "REFUSED: unattributed residual {residual:.3e} exceeds the honesty \
                     limit {:.3e}; no causal story is issued.",
                    receipt.effective_limit
                );
            }
        }
        out
    }
}

fn validate_explanation_inputs(
    nodes: &[ExplanationNode],
    observed: f64,
    threshold: f64,
) -> Result<(), ExplanationError> {
    if !observed.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "observed change",
            index: None,
            reason: "must be finite",
        });
    }
    if !threshold.is_finite() || threshold < 0.0 {
        return Err(ExplanationError::InvalidNumber {
            field: "explanation threshold",
            index: None,
            reason: "must be finite and non-negative",
        });
    }
    if !(1..=MAX_EXPLANATION_NODES).contains(&nodes.len()) {
        return Err(ExplanationError::InvalidCount {
            field: "explanation nodes",
            value: nodes.len(),
            min: 1,
            max: MAX_EXPLANATION_NODES,
        });
    }
    if !nodes_are_unique(nodes) {
        return Err(ExplanationError::DuplicateIdentity {
            field: "node fingerprint or derivation digest",
        });
    }
    if !built_in_batch_is_coherent(nodes) {
        return Err(ExplanationError::IntegrityMismatch {
            field: "built-in attribution batch",
            index: None,
        });
    }
    for (index, node) in nodes.iter().enumerate() {
        validate_node_payload(
            &node.channel,
            node.contribution,
            node.bound,
            &node.color,
            &node.evidence,
        )?;
        if !node.verifies() {
            return Err(ExplanationError::IntegrityMismatch {
                field: "node fingerprint/version",
                index: Some(index),
            });
        }
    }
    Ok(())
}

/// Assemble + gate: compute the residual against the observed change
/// and REFUSE when it exceeds `threshold`.
///
/// # Errors
/// Returns [`ExplanationError`] when inputs, retained node integrity, or
/// derived arithmetic are unusable.
pub fn finalize(
    nodes: Vec<ExplanationNode>,
    observed: f64,
    threshold: f64,
) -> Result<Explanation, ExplanationError> {
    validate_explanation_inputs(&nodes, observed, threshold)?;
    let contributions = nodes
        .iter()
        .map(|node| node.contribution)
        .collect::<Vec<_>>();
    let sum = enclosed_sum(&contributions);
    if !sum.rounded.is_finite() || !sum.lo.is_finite() || !sum.hi.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "aggregate contribution enclosure",
            index: None,
            reason: "summation must remain finite",
        });
    }
    let aggregation_roundoff = sum.roundoff();
    let certified_coverage = certified_coverage(&nodes, aggregation_roundoff);
    if !certified_coverage.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "certified coverage",
            index: None,
            reason: "must remain finite",
        });
    }
    let residual = observed - sum.rounded;
    if !residual.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "explanation residual",
            index: None,
            reason: "subtraction must remain finite",
        });
    }
    let aggregate_color = aggregate_color(&nodes);
    let effective_limit = threshold.min(certified_coverage);
    if residual.abs() > effective_limit {
        let receipt = build_receipt(
            ExplanationVariant::Refused,
            &nodes,
            observed,
            residual,
            receipt_payload(
                threshold,
                certified_coverage,
                aggregation_roundoff,
                aggregate_color,
            ),
        );
        Ok(Explanation::Refused {
            partial: nodes,
            observed,
            residual,
            receipt,
        })
    } else {
        let receipt = build_receipt(
            ExplanationVariant::Explained,
            &nodes,
            observed,
            residual,
            receipt_payload(
                threshold,
                certified_coverage,
                aggregation_roundoff,
                aggregate_color,
            ),
        );
        Ok(Explanation::Explained {
            nodes,
            observed,
            residual,
            receipt,
        })
    }
}

// ---------------------------------------------------------------------------
// Engine 1: ADJOINT attribution on the elliptic compliance fixture.
// ---------------------------------------------------------------------------

/// The 1-D elliptic fixture: `−(a u′)′ = 1`, u(0)=u(1)=0, P1 elements;
/// compliance `J = ∫ u`. Channel masks are disjoint; they may intentionally
/// omit elements so the honesty gate can expose an unattributed residual.
#[derive(Debug, Clone)]
pub struct Elliptic1d {
    n: usize,
}

impl Elliptic1d {
    /// Construct a bounded 1-D elliptic fixture.
    ///
    /// # Errors
    /// Refuses zero or oversized interior-node counts before allocation.
    pub fn new(n: usize) -> Result<Self, ExplanationError> {
        if !(1..=MAX_ELLIPTIC_INTERIOR_NODES).contains(&n) {
            return Err(ExplanationError::InvalidCount {
                field: "Elliptic1d interior nodes",
                value: n,
                min: 1,
                max: MAX_ELLIPTIC_INTERIOR_NODES,
            });
        }
        Ok(Self { n })
    }

    /// Number of interior nodes.
    #[must_use]
    pub fn interior_nodes(&self) -> usize {
        self.n
    }

    fn validate(&self) -> Result<(), ExplanationError> {
        if !(1..=MAX_ELLIPTIC_INTERIOR_NODES).contains(&self.n) {
            return Err(ExplanationError::InvalidCount {
                field: "Elliptic1d interior nodes",
                value: self.n,
                min: 1,
                max: MAX_ELLIPTIC_INTERIOR_NODES,
            });
        }
        Ok(())
    }

    fn assemble_stiffness(&self, a: &[f64]) -> Result<(f64, Vec<f64>, Vec<f64>), ExplanationError> {
        self.validate()?;
        let n = self.n;
        let element_count = n + 1;
        if a.len() != element_count {
            return Err(ExplanationError::LengthMismatch {
                field: "conductivity",
                expected: element_count,
                actual: a.len(),
            });
        }
        for (index, value) in a.iter().enumerate() {
            if !value.is_finite() || *value <= 0.0 {
                return Err(ExplanationError::InvalidNumber {
                    field: "conductivity",
                    index: Some(index),
                    reason: "must be finite and positive",
                });
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut diag = vec![0.0f64; n];
        let mut off = vec![0.0f64; n.saturating_sub(1)];
        for (e, &ae) in a.iter().enumerate() {
            let w = ae / h;
            if !w.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "assembled stiffness weight",
                    index: Some(e),
                    reason: "division overflowed",
                });
            }
            if e < n {
                let value = diag[e] + w;
                if !value.is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "stiffness diagonal",
                        index: Some(e),
                        reason: "assembly overflowed",
                    });
                }
                diag[e] = value;
            }
            if e > 0 {
                let index = e - 1;
                let value = diag[index] + w;
                if !value.is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "stiffness diagonal",
                        index: Some(index),
                        reason: "assembly overflowed",
                    });
                }
                diag[index] = value;
            }
            if e > 0 && e < n {
                let index = e - 1;
                let value = off[index] - w;
                if !value.is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "stiffness off-diagonal",
                        index: Some(index),
                        reason: "assembly overflowed",
                    });
                }
                off[index] = value;
            }
        }
        Ok((h, diag, off))
    }

    /// Solve with per-element conductivity `a` (length n+1).
    ///
    /// # Errors
    /// Refuses malformed dimensions/conductivities and non-finite assembly or
    /// pivot arithmetic.
    pub fn solve(&self, a: &[f64]) -> Result<Vec<f64>, ExplanationError> {
        let n = self.n;
        let (h, diag, off) = self.assemble_stiffness(a)?;
        let mut c = off.clone();
        let mut d = vec![h; n];
        if !diag[0].is_finite() || diag[0] <= 0.0 {
            return Err(ExplanationError::SingularPivot { index: 0 });
        }
        if n > 1 {
            c[0] /= diag[0];
            if !c[0].is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "tridiagonal forward coefficient",
                    index: Some(0),
                    reason: "division overflowed",
                });
            }
        }
        d[0] /= diag[0];
        if !d[0].is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "tridiagonal forward state",
                index: Some(0),
                reason: "division overflowed",
            });
        }
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if !m.is_finite() || m <= 0.0 {
                return Err(ExplanationError::SingularPivot { index: i });
            }
            if i < n - 1 {
                c[i] = off[i] / m;
                if !c[i].is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "tridiagonal forward coefficient",
                        index: Some(i),
                        reason: "division overflowed",
                    });
                }
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
            if !d[i].is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "tridiagonal forward state",
                    index: Some(i),
                    reason: "elimination overflowed",
                });
            }
        }
        for i in (0..n - 1).rev() {
            let t = c[i] * d[i + 1];
            d[i] -= t;
            if !t.is_finite() || !d[i].is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "tridiagonal back substitution",
                    index: Some(i),
                    reason: "arithmetic overflowed",
                });
            }
        }
        Ok(d)
    }

    /// Compliance `J = h Σ u`.
    ///
    /// # Errors
    /// Refuses a malformed/non-finite state or non-finite accumulation.
    pub fn compliance(&self, u: &[f64]) -> Result<f64, ExplanationError> {
        self.validate()?;
        if u.len() != self.n {
            return Err(ExplanationError::LengthMismatch {
                field: "elliptic state",
                expected: self.n,
                actual: u.len(),
            });
        }
        let mut sum = 0.0f64;
        for (index, value) in u.iter().enumerate() {
            if !value.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "elliptic state",
                    index: Some(index),
                    reason: "must be finite",
                });
            }
            sum += value;
            if !sum.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "compliance accumulation",
                    index: Some(index),
                    reason: "summation overflowed",
                });
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (self.n as f64 + 1.0);
        let compliance = h * sum;
        if !compliance.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "compliance",
                index: None,
                reason: "derived value is non-finite",
            });
        }
        Ok(compliance)
    }

    /// Element slope of the P1 solution.
    fn slope(&self, u: &[f64], e: usize) -> Result<f64, ExplanationError> {
        if e > self.n {
            return Err(ExplanationError::InvalidIndex {
                field: "elliptic element",
                index: e,
                upper_bound: self.n + 1,
            });
        }
        let n = self.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let lo = if e == 0 { 0.0 } else { u[e - 1] };
        let hi = if e == n { 0.0 } else { u[e] };
        let slope = (hi - lo) / h;
        if !slope.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "elliptic element slope",
                index: Some(e),
                reason: "derived value is non-finite",
            });
        }
        Ok(slope)
    }
}

fn validate_adjoint_channels(
    fixture: &Elliptic1d,
    channels: &[(&str, Vec<usize>)],
) -> Result<(), ExplanationError> {
    fixture.validate()?;
    if !(1..=MAX_EXPLANATION_NODES).contains(&channels.len()) {
        return Err(ExplanationError::InvalidCount {
            field: "adjoint channels",
            value: channels.len(),
            min: 1,
            max: MAX_EXPLANATION_NODES,
        });
    }
    let element_count = fixture.n + 1;
    let mut claimed_elements = vec![false; element_count];
    let mut channel_names = BTreeSet::new();
    for (name, elements) in channels {
        if !bounded_text_is_valid(name, MAX_CHANNEL_BYTES) {
            return Err(ExplanationError::InvalidText {
                field: "adjoint channel",
                index: None,
                reason: "must be bounded, trimmed, non-empty, and control-free",
            });
        }
        if !channel_names.insert(*name) {
            return Err(ExplanationError::DuplicateIdentity {
                field: "adjoint channel name",
            });
        }
        if !(1..=element_count).contains(&elements.len()) {
            return Err(ExplanationError::InvalidCount {
                field: "adjoint channel mask elements",
                value: elements.len(),
                min: 1,
                max: element_count,
            });
        }
        for &element in elements {
            if element >= element_count {
                return Err(ExplanationError::InvalidIndex {
                    field: "adjoint channel element",
                    index: element,
                    upper_bound: element_count,
                });
            }
            if claimed_elements[element] {
                return Err(ExplanationError::OverlappingChannelElement { element });
            }
            claimed_elements[element] = true;
        }
    }
    Ok(())
}

fn adjoint_derivation_digests(
    fixture: &Elliptic1d,
    a0: &[f64],
    a1: &[f64],
    channels: &[(&str, Vec<usize>)],
) -> (String, String) {
    let mut problem_payload = Vec::new();
    push_usize(&mut problem_payload, fixture.n);
    push_len(&mut problem_payload, a0.len());
    for &value in a0 {
        push_f64(&mut problem_payload, value);
    }
    push_len(&mut problem_payload, a1.len());
    for &value in a1 {
        push_f64(&mut problem_payload, value);
    }
    let problem_digest = derivation_digest(ADJOINT_DERIVATION_DOMAIN, &problem_payload);
    let mut batch_payload = Vec::new();
    push_str(&mut batch_payload, &problem_digest);
    push_len(&mut batch_payload, channels.len());
    for (name, elements) in channels {
        push_str(&mut batch_payload, name);
        push_len(&mut batch_payload, elements.len());
        for &element in elements {
            push_usize(&mut batch_payload, element);
        }
    }
    let batch_digest = derivation_digest(ADJOINT_DERIVATION_DOMAIN, &batch_payload);
    (problem_digest, batch_digest)
}

/// ADJOINT attribution of a conductivity edit `a0 → a1` over named
/// disjoint channel masks (element index sets). Uses the EXACT bilinear
/// identity `J(a1) − J(a0) = −∫ Δa · u0′ · u1′` (compliance is
/// self-adjoint; both states enter, no linearization error), so each supplied
/// channel follows the exact discrete identity before floating-point solve and
/// accumulation error. A complete mask set sums to the observed change; an
/// intentionally partial set leaves a residual for finalize to refuse. The v0
/// rounding allowance is heuristic, so every node remains `Estimated` until an
/// outward-rounded solver and accumulation proof is implemented.
///
/// # Errors
/// Refuses malformed fixtures, conductivities, channel masks, solver failures,
/// or non-finite derived attribution arithmetic.
pub fn adjoint_attribution(
    fixture: &Elliptic1d,
    a0: &[f64],
    a1: &[f64],
    channels: &[(&str, Vec<usize>)],
) -> Result<Vec<ExplanationNode>, ExplanationError> {
    fixture.validate()?;
    let expected = fixture.n + 1;
    if a0.len() != expected {
        return Err(ExplanationError::LengthMismatch {
            field: "initial conductivity",
            expected,
            actual: a0.len(),
        });
    }
    if a1.len() != expected {
        return Err(ExplanationError::LengthMismatch {
            field: "final conductivity",
            expected,
            actual: a1.len(),
        });
    }
    validate_adjoint_channels(fixture, channels)?;
    let u0 = fixture.solve(a0)?;
    let u1 = fixture.solve(a1)?;
    let (problem_digest, batch_digest) = adjoint_derivation_digests(fixture, a0, a1, channels);
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / (fixture.n as f64 + 1.0);
    let mut nodes = Vec::with_capacity(channels.len());
    for (index, (name, elems)) in channels.iter().enumerate() {
        let mut acc = 0.0f64;
        for &element in elems {
            let da = a1[element] - a0[element];
            let slope0 = fixture.slope(&u0, element)?;
            let slope1 = fixture.slope(&u1, element)?;
            let term = da * slope0 * slope1 * h;
            if !da.is_finite() || !term.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "adjoint channel term",
                    index: Some(element),
                    reason: "derived product is non-finite",
                });
            }
            acc -= term;
            if !acc.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "adjoint channel accumulation",
                    index: Some(index),
                    reason: "summation overflowed",
                });
            }
        }
        // Diagnostic only: this allowance has not been proved against the
        // full solve conditioning and accumulation path. It is retained for
        // prioritization, but its Estimated color cannot certify coverage.
        let bound = 1e-12 * acc.abs().max(1.0);
        if !bound.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "adjoint rounding diagnostic",
                index: Some(index),
                reason: "derived value is non-finite",
            });
        }
        let mut channel_payload = Vec::new();
        push_str(&mut channel_payload, &problem_digest);
        push_str(&mut channel_payload, name);
        push_len(&mut channel_payload, elems.len());
        for &element in elems {
            push_usize(&mut channel_payload, element);
        }
        let channel_digest = derivation_digest(ADJOINT_DERIVATION_DOMAIN, &channel_payload);
        nodes.push(ExplanationNode::new_with_authority(
            name,
            acc,
            bound,
            Color::Estimated {
                estimator: format!("derived:v2:adjoint-rounding:{channel_digest}"),
                dispersion: bound,
            },
            vec![
                format!("problem:{problem_digest}"),
                format!("derivation:{channel_digest}"),
                format!("mask:{name}"),
            ],
            ExplanationNodeAuthority::built_in(
                channel_digest,
                batch_digest.clone(),
                index,
                channels.len(),
            ),
        )?);
    }
    Ok(nodes)
}

// ---------------------------------------------------------------------------
// Engine 2: PROVENANCE attribution (which edit moved the number).
// ---------------------------------------------------------------------------

/// Telescoping edit attribution from caller-provided state values. The history
/// is input-bound and deterministic, but this API does not authenticate a
/// ledger replay. Each rounded subtraction carries a one-ulp diagnostic
/// envelope and remains `Estimated`; an authenticated replay-receipt overload
/// is required before provenance nodes may become built-in `Verified` evidence.
///
/// # Errors
/// Refuses empty/oversized, malformed, duplicated, disconnected, non-finite,
/// or arithmetically overflowing edit histories.
pub fn provenance_attribution(
    edits: &[(String, f64, f64)],
) -> Result<Vec<ExplanationNode>, ExplanationError> {
    if !(1..=MAX_EXPLANATION_NODES).contains(&edits.len()) {
        return Err(ExplanationError::InvalidCount {
            field: "provenance edits",
            value: edits.len(),
            min: 1,
            max: MAX_EXPLANATION_NODES,
        });
    }
    let mut history_payload = Vec::new();
    push_len(&mut history_payload, edits.len());
    let mut names = BTreeSet::new();
    for (index, (name, before, after)) in edits.iter().enumerate() {
        if !bounded_text_is_valid(name, MAX_CHANNEL_BYTES - "edit:".len()) {
            return Err(ExplanationError::InvalidText {
                field: "provenance edit",
                index: Some(index),
                reason: "must be bounded, trimmed, non-empty, and control-free",
            });
        }
        if !names.insert(name.as_str()) {
            return Err(ExplanationError::DuplicateIdentity {
                field: "provenance edit name",
            });
        }
        if !before.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "provenance before state",
                index: Some(index),
                reason: "must be finite",
            });
        }
        if !after.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "provenance after state",
                index: Some(index),
                reason: "must be finite",
            });
        }
        push_str(&mut history_payload, name);
        push_f64(&mut history_payload, *before);
        push_f64(&mut history_payload, *after);
    }
    for (index, (previous, next)) in edits.iter().zip(edits.iter().skip(1)).enumerate() {
        if !replay_values_match(previous.2, next.1) {
            return Err(ExplanationError::DisconnectedHistory {
                edit_index: index + 1,
            });
        }
    }
    let history_digest = derivation_digest(PROVENANCE_DERIVATION_DOMAIN, &history_payload);
    let mut nodes = Vec::with_capacity(edits.len());
    for (index, (name, before, after)) in edits.iter().enumerate() {
        let mut edit_payload = Vec::new();
        push_str(&mut edit_payload, &history_digest);
        push_usize(&mut edit_payload, index);
        let edit_digest = derivation_digest(PROVENANCE_DERIVATION_DOMAIN, &edit_payload);
        let contribution = after - before;
        if !contribution.is_finite() {
            return Err(ExplanationError::InvalidNumber {
                field: "provenance contribution",
                index: Some(index),
                reason: "subtraction overflowed",
            });
        }
        let lo = contribution.next_down();
        let hi = contribution.next_up();
        let lower_gap = (contribution - lo).abs();
        let upper_gap = (hi - contribution).abs();
        if !lo.is_finite()
            || !hi.is_finite()
            || !lower_gap.is_finite()
            || !upper_gap.is_finite()
            || lower_gap <= 0.0
            || upper_gap <= 0.0
        {
            return Err(ExplanationError::InvalidNumber {
                field: "provenance rounding envelope",
                index: Some(index),
                reason: "both adjacent values and both one-ulp gaps must be finite and positive",
            });
        }
        let bound = lower_gap.max(upper_gap);
        nodes.push(ExplanationNode::new(
            &format!("edit:{name}"),
            contribution,
            bound,
            Color::Estimated {
                estimator: format!("derived:v2:provenance-input:{history_digest}"),
                dispersion: bound,
            },
            vec![
                format!("history:{history_digest}"),
                format!("derivation:{edit_digest}"),
                format!("replay:{name}"),
            ],
        )?);
    }
    Ok(nodes)
}

// ---------------------------------------------------------------------------
// Engine 3: PHYSICAL decomposition — the far-field drag flagship.
// ---------------------------------------------------------------------------

/// Lifting-line wing fixture: span-stations with circulation Γ(y) on a
/// span `b`, freestream `v_inf`, reference area `s_ref`.
#[derive(Debug, Clone)]
pub struct LiftingLine {
    gamma: Vec<f64>,
    b: f64,
    v_inf: f64,
    s_ref: f64,
    derivation_digest: String,
}

impl LiftingLine {
    /// Elliptic distribution `Γ = Γ0 √(1 − (2y/b)²)` at `n` stations.
    ///
    /// # Errors
    /// Refuses invalid counts/parameters or non-finite derived circulation
    /// before returning a fixture.
    pub fn elliptic(
        gamma0: f64,
        b: f64,
        v_inf: f64,
        s_ref: f64,
        n: usize,
    ) -> Result<LiftingLine, ExplanationError> {
        if !(1..=MAX_LIFTING_LINE_STATIONS).contains(&n) {
            return Err(ExplanationError::InvalidCount {
                field: "lifting-line stations",
                value: n,
                min: 1,
                max: MAX_LIFTING_LINE_STATIONS,
            });
        }
        for (field, value, positive) in [
            ("elliptic circulation amplitude", gamma0, false),
            ("lifting-line span", b, true),
            ("freestream speed", v_inf, true),
            ("reference area", s_ref, true),
        ] {
            if !value.is_finite() || (positive && value <= 0.0) {
                return Err(ExplanationError::InvalidNumber {
                    field,
                    index: None,
                    reason: if positive {
                        "must be finite and positive"
                    } else {
                        "must be finite"
                    },
                });
            }
        }
        let mut gamma = Vec::with_capacity(n);
        for i in 0..n {
            #[allow(clippy::cast_precision_loss)]
            let y = -0.5 + (i as f64 + 0.5) / n as f64; // 2y/b in (−1,1)
            let value = gamma0 * (1.0 - (2.0 * y) * (2.0 * y)).max(0.0).sqrt();
            if !value.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "lifting-line circulation",
                    index: Some(i),
                    reason: "derived value is non-finite",
                });
            }
            gamma.push(value);
        }
        let mut derivation_payload = Vec::new();
        push_f64(&mut derivation_payload, gamma0);
        push_f64(&mut derivation_payload, b);
        push_f64(&mut derivation_payload, v_inf);
        push_f64(&mut derivation_payload, s_ref);
        push_usize(&mut derivation_payload, n);
        let derivation_digest =
            derivation_digest(LIFTING_LINE_DERIVATION_DOMAIN, &derivation_payload);
        Ok(LiftingLine {
            gamma,
            b,
            v_inf,
            s_ref,
            derivation_digest,
        })
    }

    fn validate(&self) -> Result<(), ExplanationError> {
        if !(1..=MAX_LIFTING_LINE_STATIONS).contains(&self.gamma.len()) {
            return Err(ExplanationError::InvalidCount {
                field: "lifting-line stations",
                value: self.gamma.len(),
                min: 1,
                max: MAX_LIFTING_LINE_STATIONS,
            });
        }
        if self.gamma.iter().any(|value| !value.is_finite())
            || !self.b.is_finite()
            || self.b <= 0.0
            || !self.v_inf.is_finite()
            || self.v_inf <= 0.0
            || !self.s_ref.is_finite()
            || self.s_ref <= 0.0
        {
            return Err(ExplanationError::InvalidNumber {
                field: "lifting-line state",
                index: None,
                reason: "circulation must be finite and dimensions must be finite and positive",
            });
        }
        if !digest_is_valid(&self.derivation_digest) {
            return Err(ExplanationError::IntegrityMismatch {
                field: "lifting-line derivation digest",
                index: None,
            });
        }
        Ok(())
    }

    /// Lift coefficient from the circulation integral (KJ theorem).
    ///
    /// # Errors
    /// Refuses non-finite derived arithmetic.
    pub fn cl(&self) -> Result<f64, ExplanationError> {
        self.validate()?;
        let mut gamma_sum = 0.0f64;
        for (index, gamma) in self.gamma.iter().enumerate() {
            gamma_sum += gamma;
            if !gamma_sum.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "circulation accumulation",
                    index: Some(index),
                    reason: "summation overflowed",
                });
            }
        }
        #[allow(clippy::cast_precision_loss)]
        let dy = self.b / self.gamma.len() as f64;
        let lift_per_rho = self.v_inf * gamma_sum * dy;
        let dynamic_reference = 0.5 * self.v_inf * self.v_inf * self.s_ref;
        let cl = lift_per_rho / dynamic_reference;
        if !dy.is_finite()
            || !lift_per_rho.is_finite()
            || !dynamic_reference.is_finite()
            || dynamic_reference <= 0.0
            || !cl.is_finite()
        {
            return Err(ExplanationError::InvalidNumber {
                field: "lift coefficient",
                index: None,
                reason: "derived arithmetic is non-finite or degenerate",
            });
        }
        Ok(cl)
    }

    /// INDUCED drag by the TREFFTZ-PLANE wake integral: the shed
    /// vorticity sheet's kinetic energy,
    /// `D_i/ρ = (1/4π) ΣΣ γ_i γ_j ln|y_i − y_j|`-free discrete form via
    /// downwash: `w(y_i) = Σ_j γ'_j / (4π (y_i − y_j))`,
    /// `D_i/ρ = Σ_i Γ_i w_i dy`. Deterministic midpoint discretization.
    ///
    /// # Errors
    /// Refuses non-finite wake assembly, accumulation, or normalization.
    pub fn induced_drag_coefficient(&self) -> Result<f64, ExplanationError> {
        self.validate()?;
        let n = self.gamma.len();
        #[allow(clippy::cast_precision_loss)]
        let dy = self.b / n as f64;
        // Shed vorticity between stations: γ_shed = −dΓ/dy at panel
        // edges (n+1 trailing vortices including tips).
        // Downwash convention (Katz & Plotkin): w(y) =
        // −(1/4π)∫(dΓ/dy′)/(y−y′) dy′; the discrete jump ΔΓ_j enters
        // NEGATED, i.e. as (Γ_j − Γ_{j−1}) with the sign folded here —
        // the original left−right form double-negated and produced a
        // wake integral of the right magnitude and wrong sign (caught
        // by the analytic envelope in conformance).
        let mut shed = Vec::with_capacity(n + 1);
        for j in 0..=n {
            let left = if j == 0 { 0.0 } else { self.gamma[j - 1] };
            let right = if j == n { 0.0 } else { self.gamma[j] };
            let value = right - left;
            if !value.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "shed circulation",
                    index: Some(j),
                    reason: "subtraction overflowed",
                });
            }
            shed.push(value);
        }
        let mut drag_per_rho = 0.0f64;
        for (i, &g) in self.gamma.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let yi = (i as f64 + 0.5) * dy;
            // Downwash at station i from every trailing vortex edge.
            let mut w = 0.0f64;
            for (j, &sv) in shed.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let yj = j as f64 * dy;
                let r = yi - yj;
                let term = sv / (4.0 * std::f64::consts::PI * r);
                if !term.is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "Trefftz downwash term",
                        index: Some(j),
                        reason: "derived arithmetic is non-finite",
                    });
                }
                w += term;
                if !w.is_finite() {
                    return Err(ExplanationError::InvalidNumber {
                        field: "Trefftz downwash accumulation",
                        index: Some(i),
                        reason: "summation overflowed",
                    });
                }
            }
            let term = g * w * dy;
            drag_per_rho += term;
            if !term.is_finite() || !drag_per_rho.is_finite() {
                return Err(ExplanationError::InvalidNumber {
                    field: "Trefftz drag accumulation",
                    index: Some(i),
                    reason: "derived arithmetic is non-finite",
                });
            }
        }
        let dynamic_reference = 0.5 * self.v_inf * self.v_inf * self.s_ref;
        let coefficient = drag_per_rho / dynamic_reference;
        if !dy.is_finite()
            || !dynamic_reference.is_finite()
            || dynamic_reference <= 0.0
            || !coefficient.is_finite()
        {
            return Err(ExplanationError::InvalidNumber {
                field: "induced drag coefficient",
                index: None,
                reason: "derived arithmetic is non-finite or degenerate",
            });
        }
        Ok(coefficient)
    }

    /// Aspect ratio.
    ///
    /// # Errors
    /// Refuses non-finite derived geometry.
    pub fn aspect_ratio(&self) -> Result<f64, ExplanationError> {
        self.validate()?;
        let aspect_ratio = self.b * self.b / self.s_ref;
        if !aspect_ratio.is_finite() || aspect_ratio <= 0.0 {
            return Err(ExplanationError::InvalidNumber {
                field: "aspect ratio",
                index: None,
                reason: "derived value must be finite and positive",
            });
        }
        Ok(aspect_ratio)
    }
}

fn drag_channel_digests(
    wing: &LiftingLine,
    cf_strip: f64,
    wetted_over_sref: f64,
    mach: f64,
    subsonic_evidence: &str,
) -> (String, String, String, String) {
    let mut drag_payload = Vec::new();
    push_str(&mut drag_payload, &wing.derivation_digest);
    push_f64(&mut drag_payload, cf_strip);
    push_f64(&mut drag_payload, wetted_over_sref);
    push_f64(&mut drag_payload, mach);
    push_str(&mut drag_payload, subsonic_evidence);
    let drag_digest = derivation_digest(DRAG_DERIVATION_DOMAIN, &drag_payload);
    let channel_digest = |channel: &str| {
        let mut payload = Vec::new();
        push_str(&mut payload, &drag_digest);
        push_str(&mut payload, channel);
        derivation_digest(DRAG_DERIVATION_DOMAIN, &payload)
    };
    let induced = channel_digest("induced");
    let viscous = channel_digest("viscous");
    let wave = channel_digest("wave");
    (drag_digest, induced, viscous, wave)
}

fn wave_explanation_node(
    mach: f64,
    subsonic_evidence: &str,
    wave_digest: String,
    batch_digest: String,
) -> Result<ExplanationNode, ExplanationError> {
    let (channel, estimator, regime) = if mach < MAX_DECLARED_SUBSONIC_MACH {
        (
            "wave (declared zero: subsonic regime)",
            format!("derived:v2:wave-subsonic-declaration:{wave_digest}"),
            format!("regime:mach<{MAX_DECLARED_SUBSONIC_MACH}"),
        )
    } else {
        (
            "wave (unresolved outside declared subsonic regime)",
            format!("derived:v2:wave-outside-regime:{wave_digest}"),
            format!("regime:mach>={MAX_DECLARED_SUBSONIC_MACH}"),
        )
    };
    ExplanationNode::new_with_authority(
        channel,
        0.0,
        0.0,
        Color::Estimated {
            estimator,
            dispersion: f64::INFINITY,
        },
        vec![
            regime,
            subsonic_evidence.to_string(),
            format!("mach-bits:{:016x}", mach.to_bits()),
            format!("derivation:{wave_digest}"),
        ],
        ExplanationNodeAuthority::built_in(wave_digest, batch_digest, 2, 3),
    )
}

/// The FLAGSHIP: decompose total drag into (induced, viscous, wave)
/// with measured/heuristic diagnostics, via the wake integral + a strip-friction model + the
/// declared-zero subsonic wave channel. `mach` and `subsonic_evidence`
/// explicitly state the regime supporting that zero; at Mach 0.8 or above the
/// wave channel demotes to unbounded `Estimated`. `cd_total_observed` is the
/// near-field measurement being explained.
///
/// # Errors
/// Refuses malformed public inputs, invalid retained wing state, non-finite
/// derived drag diagnostics, or explanation finalization failures.
pub fn drag_decomposition(
    wing: &LiftingLine,
    cf_strip: f64,
    wetted_over_sref: f64,
    mach: f64,
    subsonic_evidence: &str,
    cd_total_observed: f64,
    threshold: f64,
) -> Result<Explanation, ExplanationError> {
    for (field, value) in [
        ("strip friction coefficient", cf_strip),
        ("wetted/reference area ratio", wetted_over_sref),
        ("Mach number", mach),
        ("observed total drag", cd_total_observed),
        ("drag explanation threshold", threshold),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(ExplanationError::InvalidNumber {
                field,
                index: None,
                reason: "must be finite and non-negative",
            });
        }
    }
    if !bounded_text_is_valid(subsonic_evidence, MAX_EVIDENCE_BYTES) {
        return Err(ExplanationError::InvalidText {
            field: "subsonic evidence identity",
            index: None,
            reason: "must be bounded, trimmed, non-empty, and control-free",
        });
    }
    wing.validate()?;
    let n_stations = wing.gamma.len();
    let cdi = wing.induced_drag_coefficient()?;
    if cdi < 0.0 {
        return Err(ExplanationError::InvalidNumber {
            field: "induced drag diagnostic",
            index: None,
            reason: "must be non-negative",
        });
    }
    // Empirical midpoint-discretization diagnostic. The observed O(1/n)
    // convergence on the elliptic fixture is not an outward-rounded theorem
    // over the full public parameter domain, so this value cannot discharge
    // certified coverage or mint `Verified`.
    #[allow(clippy::cast_precision_loss)]
    let cdi_bound = cdi.abs() / n_stations as f64 * 4.0;
    let cdv = cf_strip * wetted_over_sref;
    if !cdi_bound.is_finite() || !cdv.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "drag diagnostic contribution/dispersion",
            index: None,
            reason: "derived arithmetic is non-finite",
        });
    }
    let viscous_dispersion = 0.15 * cdv;
    if !viscous_dispersion.is_finite() {
        return Err(ExplanationError::InvalidNumber {
            field: "viscous drag diagnostic dispersion",
            index: None,
            reason: "derived arithmetic is non-finite",
        });
    }
    let (batch_digest, induced_digest, viscous_digest, wave_digest) =
        drag_channel_digests(wing, cf_strip, wetted_over_sref, mach, subsonic_evidence);
    let nodes = vec![
        ExplanationNode::new_with_authority(
            "induced (Trefftz wake integral)",
            cdi,
            cdi_bound,
            Color::Estimated {
                estimator: format!("derived:v2:trefftz-midpoint-o1-diagnostic:{induced_digest}"),
                dispersion: cdi_bound,
            },
            vec![
                "wake-integral".to_string(),
                format!("stations:{n_stations}"),
                format!("derivation:{induced_digest}"),
            ],
            ExplanationNodeAuthority::built_in(induced_digest, batch_digest.clone(), 0, 3),
        )?,
        ExplanationNode::new_with_authority(
            "viscous (strip friction)",
            cdv,
            viscous_dispersion,
            Color::Estimated {
                estimator: "strip-friction".to_string(),
                dispersion: viscous_dispersion,
            },
            vec![
                format!("cf:{cf_strip}"),
                format!("swet/sref:{wetted_over_sref}"),
                format!("derivation:{viscous_digest}"),
            ],
            ExplanationNodeAuthority::built_in(viscous_digest, batch_digest.clone(), 1, 3),
        )?,
        wave_explanation_node(mach, subsonic_evidence, wave_digest, batch_digest)?,
    ];
    finalize(nodes, cd_total_observed, threshold)
}
