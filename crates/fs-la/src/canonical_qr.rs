//! Canonical rank-deficient QR: typed policy and result surface.
//!
//! Bead frankensim-epic-bedrock-6ys.5.1.2. This module freezes the CONTRACT
//! VOCABULARY that the canonical-gauge implementation ([`6ys.5.1.3`]) and the
//! independent certificate checker ([`6ys.5.1.4`]) must speak. It defines NO
//! algorithm and NO certificate checker.
//!
//! Claim vocabulary is bound to the theorem tiers frozen by
//! `tests/tsqr_rank_deficient.rs` (bead 6ys.5.1.1, commit a5d99bd6):
//!
//! * **T0** exact-arithmetic reconstruction authority,
//! * **T1** same-ISA bitwise rerun determinism,
//! * **T2** full-column-rank cross-schedule agreement + positive-diagonal
//!   uniqueness,
//! * **T3** is deliberately NOT a [`ClaimTier`] variant: it is the *absence*
//!   of a claim, represented by [`Authority::NoClaim`] with a typed reason.
//!   Encoding T3 as a claimable tier would let callers forge authority by
//!   naming it.
//!
//! # Non-forgeability law
//!
//! An [`OutcomeAuthority::Certified`] value cannot be constructed attached
//! to an outcome unless a checker receipt reference is present
//! ([`ReplayIdentity::certificate_ref`]); [`CanonicalQrOutcome::checked`]
//! refuses the pairing otherwise. Producer self-assessment never mints a
//! certified tier.
//!
//! # Dimensionful-tolerance ban
//!
//! [`RankTolerance`] has a private `f64` reachable only through
//! [`RankTolerance::relative`], which accepts scale-*relative* factors only.
//! There is no constructor from an absolute pivot cutoff, so an absolute
//! threshold literal cannot be smuggled into a policy — the failure the
//! scale-sweep fixture (`threshold_boundary_is_scale_dependent_no_claim`)
//! demonstrated.

use std::fmt;

use fs_blake3::{ContentHash, DomainHasher};

/// Semantic version of the frozen theorem tiers this surface speaks.
pub const CANONICAL_QR_THEOREM_VERSION: u32 = 1;

/// Semantic version of the wire/codec encoding below.
pub const CANONICAL_QR_SCHEMA_VERSION: u32 = 1;

/// Implementation revision of the producing code path (bumped by 6ys.5.1.3
/// when the algorithmic surface changes).
pub const CANONICAL_QR_IMPLEMENTATION_VERSION: u32 = 0;

/// Domain separating canonical-QR replay identities from every other hash
/// space in the repository.
pub const CANONICAL_QR_IDENTITY_DOMAIN: &str = "frankensim.fs-la.canonical-qr-replay.v1";

/// f64 unit roundoff, 2⁻⁵³ (same constant as `mixed.rs`).
const EPS64: f64 = 1.110_223_024_625_156_5e-16;

/// Typed refusals for every checked constructor and codec in this module.
/// Failure is DATA: no data condition panics (fs-la contract law).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolicyError {
    /// A budget/tolerance field was NaN, infinite, zero, or negative where
    /// a positive finite scale-relative factor was required.
    InvalidScaleRelativeFactor,
    /// A budget field was outside its stated admissible window.
    BudgetOutOfRange,
    /// Wire bytes name a schema version this decoder does not understand.
    UnknownSchemaVersion(u32),
    /// Wire payload is truncated, has trailing bytes, or misframes a field.
    MalformedEncoding,
    /// The R storage length does not equal the exactly-checked n*n product.
    ShapeMismatch { expected: usize, got: usize },
    /// R carries a nonzero entry below the diagonal.
    NotUpperTriangular { row: usize, col: usize },
    /// R has a strictly-negative computed diagonal (violates the flip law:
    /// strictly-negative diagonals are flipped before admission).
    StrictlyNegativeDiagonal { index: usize },
    /// A declared-certified tier was paired with no checker receipt.
    UncertifiedClaim,
    /// Identity versions are stale or cross-domain relative to this build.
    StaleIdentity { field: &'static str },
}

impl fmt::Display for PolicyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidScaleRelativeFactor => {
                write!(f, "tolerance/budget must be a positive finite scale-relative factor")
            }
            Self::BudgetOutOfRange => write!(f, "budget outside its admissible window"),
            Self::UnknownSchemaVersion(v) => write!(f, "unknown canonical-qr schema version {v}"),
            Self::MalformedEncoding => write!(f, "malformed canonical-qr encoding"),
            Self::ShapeMismatch { expected, got } => {
                write!(f, "R shape mismatch: expected {expected} entries, got {got}")
            }
            Self::NotUpperTriangular { row, col } => {
                write!(f, "R[{row}][{col}] below diagonal must be exactly zero")
            }
            Self::StrictlyNegativeDiagonal { index } => {
                write!(f, "R diagonal [{index}] strictly negative; flip law violated")
            }
            Self::UncertifiedClaim => write!(
                f,
                "certified tier requires an independent checker receipt reference"
            ),
            Self::StaleIdentity { field } => write!(f, "stale/cross-domain identity: {field}"),
        }
    }
}

impl std::error::Error for PolicyError {}

/// Scale-aware numerical rank threshold. Classification divides a candidate
/// pivot by its column-scale reference and compares against `relative`
/// times machine epsilon; ABSOLUTE cutoffs are unrepresentable (private
/// field, no absolute constructor).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RankTolerance {
    relative: f64,
}

impl RankTolerance {
    /// Admissible window for scale-relative factors: (ε₆₄, 1e6]. Below the
    /// floor the threshold is finer than rounding can resolve; above the
    /// ceiling it would classify obviously-nonzero pivots as zero.
    pub const WINDOW_UPPER: f64 = 1.0e6;

    /// Checked construction from a scale-relative factor (multiples of the
    /// column-norm reference, not of unity).
    pub fn relative(scale_relative: f64) -> Result<Self, PolicyError> {
        if !(scale_relative.is_finite() && scale_relative > EPS64 && scale_relative <= Self::WINDOW_UPPER)
        {
            return Err(PolicyError::InvalidScaleRelativeFactor);
        }
        Ok(Self { relative: scale_relative })
    }

    /// Documented default: √ε₆₄, the classical square-root-of-roundoff rank
    /// heuristic expressed RELATIVELY.
    #[must_use]
    pub fn default_f64() -> Self {
        Self { relative: EPS64.sqrt() }
    }

    /// The stored scale-relative factor (inspection only; there is no way to
    /// read an absolute cutoff back because none exists).
    #[must_use]
    pub fn factor(&self) -> f64 {
        self.relative
    }
}

/// Relative residual/error budget carried with the policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrorBudget {
    residual_relative: f64,
}

impl ErrorBudget {
    /// Admissible window: (0, 1].
    pub fn relative(residual_relative: f64) -> Result<Self, PolicyError> {
        if !(residual_relative.is_finite() && residual_relative > 0.0 && residual_relative <= 1.0)
        {
            return Err(PolicyError::BudgetOutOfRange);
        }
        Ok(Self { residual_relative })
    }

    #[must_use]
    pub fn factor(&self) -> f64 {
        self.residual_relative
    }
}

/// Claim tiers bound to the frozen theorem statement. T3 is intentionally
/// absent (see module docs): absence of claim lives in [`Authority::NoClaim`]
/// so it can never be mistaken for an affirmative result class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ClaimTier {
    /// T0: reconstruction authority (`RᵀR = AᵀA`, rank preservation).
    ExactReconstruction,
    /// T1: same-ISA bitwise rerun stability for a fixed schedule.
    SameIsaDeterministic,
    /// T2: full-column-rank cross-schedule agreement + uniqueness.
    FullRankTreeAgreement,
}

impl ClaimTier {
    /// Stable wire tag (u8) used by the codec and diagnostics.
    #[must_use]
    pub fn tag(self) -> u8 {
        match self {
            Self::ExactReconstruction => 0,
            Self::SameIsaDeterministic => 1,
            Self::FullRankTreeAgreement => 2,
        }
    }

    /// Inverse of [`ClaimTier::tag`]; `None` on unknown tags (fail-closed).
    #[must_use]
    pub fn from_tag(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::ExactReconstruction),
            1 => Some(Self::SameIsaDeterministic),
            2 => Some(Self::FullRankTreeAgreement),
            _ => None,
        }
    }
}

/// Typed reasons a result carries NO claim. Exhaustive on purpose: every
/// honest refusal in the pipeline must name one of these.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NoClaimReason {
    /// Rank-deficient input: cross-schedule factor equality is not claimed.
    RankDeficientCrossScheduleEquality,
    /// Near-boundary pivot: no rank verdict at any absolute threshold.
    AmbiguousRankBoundary,
    /// Input contained non-finite entries (outside every tier; typed refusal
    /// per the 6ys.5.1.1 status-quo pin).
    NonFiniteInput,
    /// Requested arithmetic mode is not admitted by this implementation.
    UnsupportedArithmeticMode,
}

/// The complete honest authority statement about one outcome.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OutcomeAuthority {
    /// An affirmative claim at the named tier (requires a checker receipt;
    /// enforced by [`CanonicalQrOutcome::checked`]).
    Certified(ClaimTier),
    /// Explicit, typed absence of claim.
    NoClaim(NoClaimReason),
}

impl OutcomeAuthority {
    /// Stable wire tag pair for the codec.
    #[must_use]
    pub fn tag(self) -> [u8; 2] {
        match self {
            Self::Certified(t) => [0, t.tag()],
            // Reason tags continue after the tier tag space.
            Self::NoClaim(r) => [
                1,
                match r {
                    NoClaimReason::RankDeficientCrossScheduleEquality => 0,
                    NoClaimReason::AmbiguousRankBoundary => 1,
                    NoClaimReason::NonFiniteInput => 2,
                    NoClaimReason::UnsupportedArithmeticMode => 3,
                },
            ],
        }
    }

    /// Fail-closed inverse of [`OutcomeAuthority::tag`].
    #[must_use]
    pub fn from_tag(bytes: [u8; 2]) -> Option<Self> {
        match bytes {
            [0, t] => ClaimTier::from_tag(t).map(Self::Certified),
            [1, 0] => Some(Self::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality)),
            [1, 1] => Some(Self::NoClaim(NoClaimReason::AmbiguousRankBoundary)),
            [1, 2] => Some(Self::NoClaim(NoClaimReason::NonFiniteInput)),
            [1, 3] => Some(Self::NoClaim(NoClaimReason::UnsupportedArithmeticMode)),
            _ => None,
        }
    }
}

/// Determinism class admitted by the frozen tiers. One variant today: the
/// same-ISA bitwise-stable class. Cross-ISA claims stay out until a G5
/// audit earns them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DeterminismClass {
    SameIsaBitStable,
}

/// Arithmetic mode requested of the producer. Only plain binary64 is
/// admitted at this implementation revision; directed-rounded and
/// extended-precision modes are named so refusals are typed rather than
/// silent ([`NoClaimReason::UnsupportedArithmeticMode`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ArithmeticMode {
    /// Round-to-nearest binary64 (the current `tsqr_r` reality).
    Binary64RoundToNearest,
    /// Reserved for 6ys.5.1.3's directed/outward kernel work.
    Binary64DirectedOutward,
}

impl ArithmeticMode {
    /// Whether this build admits the mode for production results.
    #[must_use]
    pub fn admitted(self) -> bool {
        matches!(self, Self::Binary64RoundToNearest)
    }

    /// Stable wire tag.
    #[must_use]
    pub fn tag(self) -> u8 {
        match self {
            Self::Binary64RoundToNearest => 0,
            Self::Binary64DirectedOutward => 1,
        }
    }

    /// Fail-closed inverse of [`ArithmeticMode::tag`].
    #[must_use]
    pub fn from_tag(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::Binary64RoundToNearest),
            1 => Some(Self::Binary64DirectedOutward),
            _ => None,
        }
    }
}

/// Deterministic tie policy (mirrors the LU lowest-index convention).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TiePolicy {
    LowestIndexFirst,
}

/// The full checked policy object. Construction is total-failure-safe: every
/// invalid combination returns a typed [`PolicyError`].
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalQrPolicy {
    rank_tolerance: RankTolerance,
    error_budget: ErrorBudget,
    determinism: DeterminismClass,
    mode: ArithmeticMode,
    ties: TiePolicy,
    theorem_version: u32,
    schema_version: u32,
}

impl CanonicalQrPolicy {
    /// Checked construction. Rejects malformed budgets and pins the current
    /// theorem/schema versions; older identities decode but do not validate
    /// against this build (see [`ReplayIdentity::validate`]).
    pub fn new(
        rank_tolerance: RankTolerance,
        error_budget: ErrorBudget,
        determinism: DeterminismClass,
        mode: ArithmeticMode,
        ties: TiePolicy,
    ) -> Result<Self, PolicyError> {
        // Field values were validated by their own constructors; assemble.
        Ok(Self {
            rank_tolerance,
            error_budget,
            determinism,
            mode,
            ties,
            theorem_version: CANONICAL_QR_THEOREM_VERSION,
            schema_version: CANONICAL_QR_SCHEMA_VERSION,
        })
    }

    #[must_use]
    pub fn rank_tolerance(&self) -> RankTolerance {
        self.rank_tolerance
    }

    #[must_use]
    pub fn error_budget(&self) -> ErrorBudget {
        self.error_budget
    }

    #[must_use]
    pub fn determinism(&self) -> DeterminismClass {
        self.determinism
    }

    #[must_use]
    pub fn arithmetic_mode(&self) -> ArithmeticMode {
        self.mode
    }

    #[must_use]
    pub fn tie_policy(&self) -> TiePolicy {
        self.ties
    }

    #[must_use]
    pub fn theorem_version(&self) -> u32 {
        self.theorem_version
    }

    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Canonical little-endian encoding (fixed field order; the codec's
    /// canonical form is also the digest preimage).
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(64);
        out.extend_from_slice(&self.schema_version.to_le_bytes());
        out.extend_from_slice(&self.theorem_version.to_le_bytes());
        out.extend_from_slice(&self.rank_tolerance.relative.to_le_bytes());
        out.extend_from_slice(&self.error_budget.residual_relative.to_le_bytes());
        out.push(self.determinism.tag_u8());
        out.push(self.mode.tag());
        out.push(self.ties.tag_u8());
        out
    }

    /// Decode with fail-closed framing: exact length, known versions,
    /// revalidated field windows (a bit-flipped payload must refuse, not
    /// silently construct nonsense).
    pub fn decode(bytes: &[u8]) -> Result<Self, PolicyError> {
        let expected = 4 + 4 + 8 + 8 + 3;
        if bytes.len() != expected {
            return Err(PolicyError::MalformedEncoding);
        }
        let schema = u32::from_le_bytes(bytes[0..4].try_into().expect("framed"));
        if schema != CANONICAL_QR_SCHEMA_VERSION {
            return Err(PolicyError::UnknownSchemaVersion(schema));
        }
        let theorem = u32::from_le_bytes(bytes[4..8].try_into().expect("framed"));
        if theorem != CANONICAL_QR_THEOREM_VERSION {
            return Err(PolicyError::StaleIdentity { field: "theorem_version" });
        }
        let tol = f64::from_le_bytes(bytes[8..16].try_into().expect("framed"));
        let budget = f64::from_le_bytes(bytes[16..24].try_into().expect("framed"));
        let determinism = DeterminismClass::from_tag_u8(bytes[24]).ok_or(PolicyError::MalformedEncoding)?;
        let mode = ArithmeticMode::from_tag(bytes[25]).ok_or(PolicyError::MalformedEncoding)?;
        let ties = TiePolicy::from_tag_u8(bytes[26]).ok_or(PolicyError::MalformedEncoding)?;
        Ok(Self {
            rank_tolerance: RankTolerance::relative(tol)?,
            error_budget: ErrorBudget::relative(budget)?,
            determinism,
            mode,
            ties,
            theorem_version: theorem,
            schema_version: schema,
        })
    }

    /// Content identity of the policy alone (domain-separated).
    #[must_use]
    pub fn digest(&self) -> ContentHash {
        let mut h = DomainHasher::new(CANONICAL_QR_IDENTITY_DOMAIN);
        h.update(b"policy:");
        h.update(&self.encode());
        h.finalize()
    }
}

impl DeterminismClass {
    fn tag_u8(self) -> u8 {
        match self {
            Self::SameIsaBitStable => 0,
        }
    }

    fn from_tag_u8(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::SameIsaBitStable),
            _ => None,
        }
    }
}

impl TiePolicy {
    fn tag_u8(self) -> u8 {
        match self {
            Self::LowestIndexFirst => 0,
        }
    }

    fn from_tag_u8(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::LowestIndexFirst),
            _ => None,
        }
    }
}

/// Per-pivot classification produced by the (future) rank profiler. Every
/// class names what may be claimed about that diagonal position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PivotClass {
    /// Clearly separated from zero under the scale-aware threshold.
    Nonzero,
    /// Structurally/exactly zero (e.g. dependent-column elimination).
    Zero,
    /// Inside the ambiguity band: NO verdict (feeds
    /// [`NoClaimReason::AmbiguousRankBoundary`]).
    Ambiguous,
}

impl PivotClass {
    fn tag(self) -> u8 {
        match self {
            Self::Nonzero => 0,
            Self::Zero => 1,
            Self::Ambiguous => 2,
        }
    }

    fn from_tag(tag: u8) -> Option<Self> {
        match tag {
            0 => Some(Self::Nonzero),
            1 => Some(Self::Zero),
            2 => Some(Self::Ambiguous),
            _ => None,
        }
    }
}

/// Typed rank profile: the count and per-position classes. `rank` MUST equal
/// the number of [`PivotClass::Nonzero`] entries — enforced, because a
/// profile whose headline number disagrees with its own detail rows is a
/// forged certificate waiting to happen.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CertifiedRankProfile {
    rank: usize,
    pivots: Vec<PivotClass>,
}

impl CertifiedRankProfile {
    /// Checked construction enforcing rank == count(Nonzero).
    pub fn checked(pivots: Vec<PivotClass>) -> Result<Self, PolicyError> {
        let rank = pivots.iter().filter(|p| **p == PivotClass::Nonzero).count();
        Ok(Self { rank, pivots })
    }

    #[must_use]
    pub fn rank(&self) -> usize {
        self.rank
    }

    #[must_use]
    pub fn pivots(&self) -> &[PivotClass] {
        &self.pivots
    }

    /// True iff any position refused a verdict — callers must then treat the
    /// whole profile as rank-unresolved for certification purposes.
    #[must_use]
    pub fn has_ambiguity(&self) -> bool {
        self.pivots.contains(&PivotClass::Ambiguous)
    }

    fn encode_into(&self, out: &mut Vec<u8>) {
        out.extend_from_slice(&(self.pivots.len() as u64).to_le_bytes());
        out.extend(self.pivots.iter().map(|p| p.tag()));
    }

    fn decode_from(bytes: &[u8], cursor: &mut usize) -> Result<Self, PolicyError> {
        let end = *cursor + 8;
        let len_bytes = bytes.get(*cursor..end).ok_or(PolicyError::MalformedEncoding)?;
        let len = u64::from_le_bytes(len_bytes.try_into().expect("framed")) as usize;
        // Refuse giant infallible allocations: length must fit the payload.
        if bytes.len() < end + len {
            return Err(PolicyError::MalformedEncoding);
        }
        let mut pivots = Vec::with_capacity(len);
        for b in &bytes[end..end + len] {
            pivots.push(PivotClass::from_tag(*b).ok_or(PolicyError::MalformedEncoding)?);
        }
        *cursor = end + len;
        Self::checked(pivots)
    }
}

/// Replay identity binding input, policy, logical tree, arithmetic mode,
/// versions, result, and (when present) the independent checker receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplayIdentity {
    /// Caller-supplied digest over raw input bits (row-major f64 le bytes).
    pub input_digest: ContentHash,
    /// Digest over the leaf-partition/tree schedule actually scheduled.
    pub tree_digest: ContentHash,
    /// Digest over the canonical outcome encoding.
    pub result_digest: ContentHash,
    /// Independent checker receipt reference; `None` caps authority at
    /// no-claim (see the non-forgeability law).
    pub certificate_ref: Option<ContentHash>,
    /// Mode actually executed (may be a refused-mode echo in a NoClaim
    /// outcome).
    pub arithmetic_mode: ArithmeticMode,
}

impl ReplayIdentity {
    /// Fail-closed coherence check against this build's versions.
    pub fn validate(&self, policy: &CanonicalQrPolicy) -> Result<(), PolicyError> {
        if policy.theorem_version() != CANONICAL_QR_THEOREM_VERSION {
            return Err(PolicyError::StaleIdentity { field: "theorem_version" });
        }
        if policy.schema_version() != CANONICAL_QR_SCHEMA_VERSION {
            return Err(PolicyError::StaleIdentity { field: "schema_version" });
        }
        if policy.arithmetic_mode() != self.arithmetic_mode {
            return Err(PolicyError::StaleIdentity { field: "arithmetic_mode" });
        }
        Ok(())
    }

    /// Composite domain-separated identity (what a ledger row keys on).
    #[must_use]
    pub fn composite_digest(&self) -> ContentHash {
        let mut h = DomainHasher::new(CANONICAL_QR_IDENTITY_DOMAIN);
        h.update(b"replay:");
        h.update(self.input_digest.as_bytes());
        h.update(self.tree_digest.as_bytes());
        h.update(self.result_digest.as_bytes());
        const NO_CERT: [u8; 32] = [0u8; 32];
        h.update(match self.certificate_ref.as_ref() {
            Some(c) => c.as_bytes(),
            None => &NO_CERT,
        });
        h.update(&[self.arithmetic_mode.tag()]);
        h.update(&CANONICAL_QR_IMPLEMENTATION_VERSION.to_le_bytes());
        h.finalize()
    }
}

/// The typed outcome surface. Construct ONLY through
/// [`CanonicalQrOutcome::checked`]; every field invariant is enforced there.
#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalQrOutcome {
    n: usize,
    r_factor: Vec<f64>,
    rank_profile: CertifiedRankProfile,
    authority: OutcomeAuthority,
    replay: ReplayIdentity,
}

impl CanonicalQrOutcome {
    /// Checked construction enforcing, in order: exactly-checked shape
    /// (`n*n`, refusing overflow and wrong lengths), upper-triangularity
    /// (exact zeros below the diagonal), the strictly-negative-diagonal flip
    /// law, finiteness (non-finite entries are a typed refusal, never a
    /// silent carrier), rank-profile consistency, and the non-forgeability
    /// law binding certified tiers to checker receipts.
    pub fn checked(
        r_factor: Vec<f64>,
        n: usize,
        rank_profile: CertifiedRankProfile,
        authority: OutcomeAuthority,
        replay: ReplayIdentity,
    ) -> Result<Self, PolicyError> {
        let Some(expected) = n.checked_mul(n) else {
            return Err(PolicyError::ShapeMismatch { expected: usize::MAX, got: r_factor.len() });
        };
        if r_factor.len() != expected {
            return Err(PolicyError::ShapeMismatch { expected, got: r_factor.len() });
        }
        if rank_profile.pivots().len() != n {
            return Err(PolicyError::ShapeMismatch { expected: n, got: rank_profile.pivots().len() });
        }
        for i in 0..n {
            if r_factor[i * n + i] < 0.0 {
                return Err(PolicyError::StrictlyNegativeDiagonal { index: i });
            }
        }
        for i in 1..n {
            for j in 0..i {
                if r_factor[i * n + j] != 0.0 {
                    return Err(PolicyError::NotUpperTriangular { row: i, col: j });
                }
            }
        }
        if r_factor.iter().any(|v| !v.is_finite()) {
            // A non-finite factor cannot carry ANY tier; the caller asked to
            // admit it as an outcome, which the contract refuses outright.
            return Err(PolicyError::InvalidScaleRelativeFactor);
        }
        // Non-forgeability: certified tiers demand an independent receipt.
        if matches!(authority, OutcomeAuthority::Certified(_)) && replay.certificate_ref.is_none()
        {
            return Err(PolicyError::UncertifiedClaim);
        }
        Ok(Self { n, r_factor, rank_profile, authority, replay })
    }

    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Row-major n×n upper-triangular R.
    #[must_use]
    pub fn r_factor(&self) -> &[f64] {
        &self.r_factor
    }

    #[must_use]
    pub fn rank_profile(&self) -> &CertifiedRankProfile {
        &self.rank_profile
    }

    #[must_use]
    pub fn authority(&self) -> OutcomeAuthority {
        self.authority
    }

    #[must_use]
    pub fn replay(&self) -> &ReplayIdentity {
        &self.replay
    }

    /// Canonical outcome encoding (digest preimage): dims, factor bits in
    /// row-major order, profile, authority tags, identity fields.
    #[must_use]
    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(16 + self.r_factor.len() * 8 + self.rank_profile.pivots().len() + 96);
        out.extend_from_slice(&(self.n as u64).to_le_bytes());
        for v in &self.r_factor {
            out.extend_from_slice(&v.to_le_bytes());
        }
        self.rank_profile.encode_into(&mut out);
        out.extend_from_slice(&self.authority.tag());
        out.extend_from_slice(self.replay.input_digest.as_bytes());
        out.extend_from_slice(self.replay.tree_digest.as_bytes());
        out.extend_from_slice(&[self.replay.arithmetic_mode.tag()]);
        out
    }

    /// Result digest over the canonical encoding (domain-separated).
    #[must_use]
    pub fn result_digest(&self) -> ContentHash {
        let mut h = DomainHasher::new(CANONICAL_QR_IDENTITY_DOMAIN);
        h.update(b"result:");
        h.update(&self.encode());
        h.finalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_tolerance_rejects_absolute_and_degenerate_values() {
        assert_eq!(RankTolerance::relative(f64::NAN), Err(PolicyError::InvalidScaleRelativeFactor));
        assert_eq!(RankTolerance::relative(f64::INFINITY), Err(PolicyError::InvalidScaleRelativeFactor));
        assert_eq!(RankTolerance::relative(0.0), Err(PolicyError::InvalidScaleRelativeFactor));
        assert_eq!(RankTolerance::relative(-1e-9), Err(PolicyError::InvalidScaleRelativeFactor));
        assert_eq!(RankTolerance::relative(EPS64), Err(PolicyError::InvalidScaleRelativeFactor));
        assert!(RankTolerance::relative(RankTolerance::WINDOW_UPPER).is_ok());
        assert!(RankTolerance::relative(10.0).is_ok());
    }

    #[test]
    fn policy_codec_roundtrip_and_failures() {
        let policy = CanonicalQrPolicy::new(
            RankTolerance::default_f64(),
            ErrorBudget::relative(1e-12).expect("in window"),
            DeterminismClass::SameIsaBitStable,
            ArithmeticMode::Binary64RoundToNearest,
            TiePolicy::LowestIndexFirst,
        )
        .expect("valid");
        let bytes = policy.encode();
        assert_eq!(CanonicalQrPolicy::decode(&bytes), Ok(policy.clone()));
        // Truncation and trailing bytes refuse.
        assert_eq!(CanonicalQrPolicy::decode(&bytes[..bytes.len() - 1]), Err(PolicyError::MalformedEncoding));
        let mut trailing = bytes.clone();
        trailing.push(0);
        assert_eq!(CanonicalQrPolicy::decode(&trailing), Err(PolicyError::MalformedEncoding));
        // Bit-flipped tolerance revalidates through the checked constructor.
        let mut flipped = bytes.clone();
        flipped[8] ^= 0x80; // sign bit of the tolerance field
        assert_eq!(
            CanonicalQrPolicy::decode(&flipped),
            Err(PolicyError::InvalidScaleRelativeFactor)
        );
        // Unknown schema version fails closed.
        let mut foreign = bytes;
        foreign[0] = 0x7F;
        assert!(matches!(
            CanonicalQrPolicy::decode(&foreign),
            Err(PolicyError::UnknownSchemaVersion(_))
        ));
    }

    #[test]
    fn policy_digest_moves_with_every_field() {
        let base = CanonicalQrPolicy::new(
            RankTolerance::default_f64(),
            ErrorBudget::relative(1e-12).expect("in window"),
            DeterminismClass::SameIsaBitStable,
            ArithmeticMode::Binary64RoundToNearest,
            TiePolicy::LowestIndexFirst,
        )
        .expect("valid");
        let d0 = base.digest();
        let variants = [
            CanonicalQrPolicy::new(
                RankTolerance::relative(1e-6).expect("in window"),
                ErrorBudget::relative(1e-12).expect("in window"),
                DeterminismClass::SameIsaBitStable,
                ArithmeticMode::Binary64RoundToNearest,
                TiePolicy::LowestIndexFirst,
            )
            .expect("valid"),
            CanonicalQrPolicy::new(
                RankTolerance::default_f64(),
                ErrorBudget::relative(1e-10).expect("in window"),
                DeterminismClass::SameIsaBitStable,
                ArithmeticMode::Binary64RoundToNearest,
                TiePolicy::LowestIndexFirst,
            )
            .expect("valid"),
            CanonicalQrPolicy::new(
                RankTolerance::default_f64(),
                ErrorBudget::relative(1e-12).expect("in window"),
                DeterminismClass::SameIsaBitStable,
                ArithmeticMode::Binary64DirectedOutward,
                TiePolicy::LowestIndexFirst,
            )
            .expect("valid"),
        ];
        for v in &variants {
            assert_ne!(v.digest(), d0, "digest collision across policy mutation");
        }
    }

    #[test]
    fn outcome_checked_enforces_shape_triangularity_and_law() {
        let n = 2usize;
        let ok = vec![1.0, 0.5, 0.0, 0.0];
        let profile = CertifiedRankProfile::checked(vec![PivotClass::Nonzero, PivotClass::Zero])
            .expect("consistent");
        let identity = ReplayIdentity {
            input_digest: fs_blake3::hash_bytes(b"a"),
            tree_digest: fs_blake3::hash_bytes(b"t"),
            result_digest: fs_blake3::hash_bytes(b"r"),
            certificate_ref: None,
            arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
        };
        let policy = CanonicalQrPolicy::new(
            RankTolerance::default_f64(),
            ErrorBudget::relative(1e-12).expect("in window"),
            DeterminismClass::SameIsaBitStable,
            ArithmeticMode::Binary64RoundToNearest,
            TiePolicy::LowestIndexFirst,
        )
        .expect("valid");
        identity.validate(&policy).expect("coherent");

        // Wrong length refuses.
        assert!(matches!(
            CanonicalQrOutcome::checked(vec![1.0], n, profile.clone(), OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality), identity.clone()),
            Err(PolicyError::ShapeMismatch { .. })
        ));
        // Lower-triangular dirt refuses.
        let dirty = vec![1.0, 0.5, 1e-30, 0.0];
        assert_eq!(
            CanonicalQrOutcome::checked(dirty, n, profile.clone(), OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality), identity.clone()),
            Err(PolicyError::NotUpperTriangular { row: 1, col: 0 })
        );
        // Strictly-negative diagonal refuses.
        let neg = vec![-1.0, 0.5, 0.0, -2.0];
        assert_eq!(
            CanonicalQrOutcome::checked(neg, n, profile.clone(), OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality), identity.clone()),
            Err(PolicyError::StrictlyNegativeDiagonal { index: 0 })
        );
        // Non-finite entries refuse.
        let nan = vec![1.0, f64::NAN, 0.0, 0.0];
        assert!(matches!(
            CanonicalQrOutcome::checked(nan, n, profile.clone(), OutcomeAuthority::NoClaim(NoClaimReason::NonFiniteInput), identity.clone()),
            Err(PolicyError::InvalidScaleRelativeFactor)
        ));
        // Valid no-claim outcome admits.
        let outcome = CanonicalQrOutcome::checked(
            ok,
            n,
            profile,
            OutcomeAuthority::NoClaim(NoClaimReason::RankDeficientCrossScheduleEquality),
            identity,
        )
        .expect("valid");
        assert_eq!(outcome.result_digest(), outcome.result_digest());
    }

    #[test]
    fn authority_non_forgeability_requires_certificate() {
        let profile = CertifiedRankProfile::checked(vec![PivotClass::Nonzero]).expect("consistent");
        let identity = ReplayIdentity {
            input_digest: fs_blake3::hash_bytes(b"a"),
            tree_digest: fs_blake3::hash_bytes(b"t"),
            result_digest: fs_blake3::hash_bytes(b"r"),
            certificate_ref: None,
            arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
        };
        // Certified WITHOUT a receipt refuses: producers cannot mint T2.
        assert_eq!(
            CanonicalQrOutcome::checked(
                vec![1.0],
                1,
                profile.clone(),
                OutcomeAuthority::Certified(ClaimTier::FullRankTreeAgreement),
                identity.clone()
            ),
            Err(PolicyError::UncertifiedClaim)
        );
        // With a receipt reference the pairing admits (the CHECKER decides
        // whether the receipt itself is valid — that is 6ys.5.1.4's job).
        let certified_identity = ReplayIdentity {
            certificate_ref: Some(fs_blake3::hash_bytes(b"checker-receipt")),
            ..identity
        };
        assert!(
            CanonicalQrOutcome::checked(
                vec![1.0],
                1,
                profile,
                OutcomeAuthority::Certified(ClaimTier::FullRankTreeAgreement),
                certified_identity,
            )
            .is_ok()
        );
    }

    #[test]
    fn rank_profile_consistency_and_ambiguity() {
        let bad = vec![PivotClass::Nonzero, PivotClass::Nonzero];
        let p = CertifiedRankProfile::checked(bad).expect("consistent");
        assert_eq!(p.rank(), 2);
        let mixed = vec![PivotClass::Nonzero, PivotClass::Ambiguous, PivotClass::Zero];
        let p2 = CertifiedRankProfile::checked(mixed).expect("consistent");
        assert_eq!(p2.rank(), 1);
        assert!(p2.has_ambiguity());
    }

    #[test]
    fn composite_digest_binds_every_field() {
        let mk = |cert: Option<ContentHash>, tree: &str| ReplayIdentity {
            input_digest: fs_blake3::hash_bytes(b"input"),
            tree_digest: fs_blake3::hash_bytes(tree.as_bytes()),
            result_digest: fs_blake3::hash_bytes(b"result"),
            certificate_ref: cert,
            arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
        };
        let base = mk(Some(fs_blake3::hash_bytes(b"c")), "tree");
        let d0 = base.composite_digest();
        // Certificate removal moves the identity (absent ≠ present).
        assert_ne!(mk(None, "tree").composite_digest(), d0);
        // Tree change moves the identity.
        assert_ne!(mk(Some(fs_blake3::hash_bytes(b"c")), "other").composite_digest(), d0);
    }
}
