//! Orthogonal truth semantics for spectral results.
//!
//! Coverage, epistemic authority, cluster resolution, and termination are not
//! one total order. A certified partial enclosure and an estimated complete
//! candidate list are incomparable. This module keeps those axes separate and
//! models repeated eigenvalues as set-valued clusters rather than unstable
//! individual-vector identities.

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, LimitKind, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::admission::{
    AdmittedSpectralWitnessV1, CompletenessScopeV1, DescriptorRoleV1, InfiniteEigenvaluePolicyV1,
    RegionBoundaryPolicyV1, SpectralNormId, SpectralOrderingV1, SpectralProblemId,
    SpectralPropositionId, SpectralPropositionKindV1, ValidatedSpectralProblemV1,
    spectral_proposition_receipt,
};

const TRUTH_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 20, 1 << 20, 16, 8192, 8192);
/// Maximum cluster records accepted before sorting, hashing, or pair checks.
pub const MAX_SPECTRAL_CLUSTERS_V1: usize = 4096;
/// Maximum cluster references accepted in a region-boundary classification.
pub const MAX_REGION_BOUNDARY_REFERENCES_V1: usize = 4096;

/// Domain-separated identity schema for a canonical set-valued spectral
/// result before whole-result evidence is attached.
pub enum SpectralResultSetIdentitySchemaV1 {}

impl CanonicalSchema for SpectralResultSetIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.result-set.v1";
    const NAME: &'static str = "spectral-result-set";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str =
        "canonical spectral clusters, enclosures, multiplicities, and internal states";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("clusters", WireType::CanonicalSet)];
}

/// Typed semantic identity of a canonical result set. Identity is not result
/// authority.
pub type SpectralResultSetIdV1 = ProblemSemanticId<SpectralResultSetIdentitySchemaV1>;

fn push_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn canonical_f64_bits(value: f64) -> u64 {
    if value.to_bits() == (-0.0f64).to_bits() {
        0.0f64.to_bits()
    } else {
        value.to_bits()
    }
}

/// Canonical membership/lineage identity of one set-valued spectral cluster.
/// Evidence producers must mint a new identity whenever the represented member
/// set changes; enclosure refinement may retain it only while membership is
/// unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpectralClusterIdV1([u8; 32]);

impl SpectralClusterIdV1 {
    /// Construct from exact lineage digest bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Exact lineage digest bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Finite closed interval with canonical signed-zero handling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteIntervalV1 {
    lower: f64,
    upper: f64,
}

impl FiniteIntervalV1 {
    /// Construct `[lower, upper]`; both endpoints must be finite and ordered.
    ///
    /// # Errors
    ///
    /// Returns [`SpectralTruthErrorV1::NonFiniteInterval`] for a NaN or
    /// infinite endpoint and [`SpectralTruthErrorV1::ReversedInterval`] when
    /// `lower > upper`.
    pub fn new(lower: f64, upper: f64) -> Result<Self, SpectralTruthErrorV1> {
        if !(lower.is_finite() && upper.is_finite()) {
            return Err(SpectralTruthErrorV1::NonFiniteInterval);
        }
        if lower > upper {
            return Err(SpectralTruthErrorV1::ReversedInterval);
        }
        Ok(Self {
            lower: f64::from_bits(canonical_f64_bits(lower)),
            upper: f64::from_bits(canonical_f64_bits(upper)),
        })
    }

    /// Lower endpoint.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Upper endpoint.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    /// Whether this interval intersects another closed interval.
    #[must_use]
    pub fn intersects(self, other: Self) -> bool {
        self.lower <= other.upper && other.lower <= self.upper
    }
}

/// Set-valued enclosure for a real, complex, or infinite eigenvalue cluster.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectralEnclosureV1 {
    /// Real interval.
    Real(FiniteIntervalV1),
    /// Axis-aligned complex rectangle.
    ComplexBox {
        /// Real-part interval.
        real: FiniteIntervalV1,
        /// Imaginary-part interval.
        imag: FiniteIntervalV1,
    },
    /// Projective point at infinity for a pencil/polynomial problem.
    ProjectiveInfinity,
}

fn push_interval(out: &mut Vec<u8>, interval: FiniteIntervalV1) {
    push_u64(out, canonical_f64_bits(interval.lower));
    push_u64(out, canonical_f64_bits(interval.upper));
}

fn push_enclosure(out: &mut Vec<u8>, enclosure: SpectralEnclosureV1) {
    match enclosure {
        SpectralEnclosureV1::Real(interval) => {
            out.push(0);
            push_interval(out, interval);
        }
        SpectralEnclosureV1::ComplexBox { real, imag } => {
            out.push(1);
            push_interval(out, real);
            push_interval(out, imag);
        }
        SpectralEnclosureV1::ProjectiveInfinity => out.push(2),
    }
}

/// Unvalidated algebraic or geometric multiplicity claim carried by a draft.
/// Favorable semantics become inspectable as result truth only through
/// [`ValidatedSpectralClusterV1`] after every witness binding has passed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable multiplicity variants intentionally retain their admitted evidence receipt inline"
)]
pub enum MultiplicityClaimV1 {
    /// Multiplicity unresolved.
    Unknown,
    /// Admitted positive lower bound only.
    LowerBound {
        /// Proven lower bound.
        value: u32,
        /// Evidence bound to this exact cluster and statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Admitted inclusive positive lower/upper bounds.
    Bounds {
        /// Lower bound.
        lower: u32,
        /// Upper bound.
        upper: u32,
        /// Evidence bound to this exact cluster and statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Exact positive multiplicity with supporting evidence.
    Exact {
        /// Exact value.
        value: u32,
        /// Evidence bound to this exact cluster and statement.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// Witness-free multiplicity semantics embedded into another proposition.
/// Internal separation and degeneracy receipts bind both multiplicity axes
/// without recursively embedding their authority records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiplicityStatementV1 {
    /// Multiplicity unresolved.
    Unknown,
    /// Positive lower bound with no finite upper bound.
    LowerBound {
        /// Proven lower bound.
        value: u32,
    },
    /// Inclusive positive lower and upper bounds.
    Bounds {
        /// Proven lower bound.
        lower: u32,
        /// Proven upper bound.
        upper: u32,
    },
    /// One exact positive value.
    Exact {
        /// Exact multiplicity.
        value: u32,
    },
}

impl MultiplicityClaimV1 {
    fn validate(self) -> Result<(), SpectralTruthErrorV1> {
        match self {
            Self::LowerBound { value: 0, .. }
            | Self::Bounds { lower: 0, .. }
            | Self::Bounds { upper: 0, .. }
            | Self::Exact { value: 0, .. } => Err(SpectralTruthErrorV1::ZeroMultiplicity),
            Self::Bounds { lower, upper, .. } if lower > upper => {
                Err(SpectralTruthErrorV1::ReversedMultiplicityBounds)
            }
            Self::Unknown | Self::LowerBound { .. } | Self::Bounds { .. } | Self::Exact { .. } => {
                Ok(())
            }
        }
    }

    const fn minimum(self) -> Option<u32> {
        match self {
            Self::Unknown => None,
            Self::LowerBound { value, .. } | Self::Exact { value, .. } => Some(value),
            Self::Bounds { lower, .. } => Some(lower),
        }
    }

    const fn maximum(self) -> Option<u32> {
        match self {
            Self::Unknown | Self::LowerBound { .. } => None,
            Self::Bounds { upper: value, .. } | Self::Exact { value, .. } => Some(value),
        }
    }

    const fn exact(self) -> Option<u32> {
        match self {
            Self::Exact { value, .. } => Some(value),
            Self::Unknown | Self::LowerBound { .. } | Self::Bounds { .. } => None,
        }
    }

    const fn statement(self) -> MultiplicityStatementV1 {
        match self {
            Self::Unknown => MultiplicityStatementV1::Unknown,
            Self::LowerBound { value, .. } => MultiplicityStatementV1::LowerBound { value },
            Self::Bounds { lower, upper, .. } => MultiplicityStatementV1::Bounds { lower, upper },
            Self::Exact { value, .. } => MultiplicityStatementV1::Exact { value },
        }
    }

    fn push_semantics(self, out: &mut Vec<u8>) {
        match self {
            Self::Unknown => out.push(0),
            Self::LowerBound { value, .. } => {
                out.push(1);
                push_u32(out, value);
            }
            Self::Bounds { lower, upper, .. } => {
                out.push(2);
                push_u32(out, lower);
                push_u32(out, upper);
            }
            Self::Exact { value, .. } => {
                out.push(3);
                push_u32(out, value);
            }
        }
    }

    const fn witness(self) -> Option<AdmittedSpectralWitnessV1> {
        match self {
            Self::Unknown => None,
            Self::LowerBound { witness, .. }
            | Self::Bounds { witness, .. }
            | Self::Exact { witness, .. } => Some(witness),
        }
    }
}

fn push_multiplicity_statement(out: &mut Vec<u8>, statement: MultiplicityStatementV1) {
    match statement {
        MultiplicityStatementV1::Unknown => out.push(0),
        MultiplicityStatementV1::LowerBound { value } => {
            out.push(1);
            push_u32(out, value);
        }
        MultiplicityStatementV1::Bounds { lower, upper } => {
            out.push(2);
            push_u32(out, lower);
            push_u32(out, upper);
        }
        MultiplicityStatementV1::Exact { value } => {
            out.push(3);
            push_u32(out, value);
        }
    }
}

/// Algebraic versus geometric multiplicity proposition target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiplicityKindV1 {
    /// Algebraic multiplicity.
    Algebraic,
    /// Geometric multiplicity.
    Geometric,
}

/// Shape of a multiplicity proposition; exact `m` cannot be replayed as a
/// generic interval `[m,m]` or vice versa.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiplicityAssertionV1 {
    /// A positive lower bound with no finite upper bound.
    LowerBound,
    /// Inclusive positive lower and upper bounds.
    Bounds,
    /// One exact positive multiplicity.
    Exact,
}

impl MultiplicityAssertionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::LowerBound => 0,
            Self::Bounds => 1,
            Self::Exact => 2,
        }
    }
}

impl MultiplicityKindV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Algebraic => 0,
            Self::Geometric => 1,
        }
    }
}

/// What can be concluded about defectivity from validated exact
/// multiplicities. Favorable variants are emitted only by
/// [`ValidatedSpectralClusterV1::defectivity`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DefectivityStateV1 {
    /// Exact geometric multiplicity is smaller than algebraic multiplicity.
    ProvenDefective,
    /// Exact algebraic and geometric multiplicities are equal.
    ProvenSemisimple,
    /// At least one multiplicity remains non-exact.
    Unknown,
}

/// Epistemic status of one cluster localization. Favorable variants retain an
/// admitted proposition-bound witness; candidate localization is explicitly
/// non-authoritative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalizationAuthorityV1 {
    /// Raw candidate region.
    Candidate,
    /// Model-based estimate admitted under a named verifier/policy.
    Estimated,
    /// Rigorous enclosure admitted under a named verifier/policy.
    Enclosed,
}

/// Unvalidated set-valued cluster-localization draft.
///
/// Construction records a claimed authority and witness, but does not prove
/// that the witness names this problem, cluster, enclosure, or authority
/// family. Favorable fields are intentionally not inspectable until validation
/// produces a [`ValidatedSpectralClusterV1`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralLocalizationV1 {
    enclosure: SpectralEnclosureV1,
    authority: LocalizationAuthorityV1,
    witness: Option<AdmittedSpectralWitnessV1>,
}

impl SpectralLocalizationV1 {
    /// Construct a non-authoritative candidate localization.
    #[must_use]
    pub const fn candidate(enclosure: SpectralEnclosureV1) -> Self {
        Self {
            enclosure,
            authority: LocalizationAuthorityV1::Candidate,
            witness: None,
        }
    }

    /// Construct a draft model-based estimate claim. Exact witness binding is
    /// checked only when the containing truth draft is validated.
    #[must_use]
    pub const fn estimated(
        enclosure: SpectralEnclosureV1,
        witness: AdmittedSpectralWitnessV1,
    ) -> Self {
        Self {
            enclosure,
            authority: LocalizationAuthorityV1::Estimated,
            witness: Some(witness),
        }
    }

    /// Construct a draft rigorous-enclosure claim. Exact witness binding is
    /// checked only when the containing truth draft is validated.
    #[must_use]
    pub const fn enclosed(
        enclosure: SpectralEnclosureV1,
        witness: AdmittedSpectralWitnessV1,
    ) -> Self {
        Self {
            enclosure,
            authority: LocalizationAuthorityV1::Enclosed,
            witness: Some(witness),
        }
    }

    /// Enclosed/candidate set.
    #[must_use]
    pub const fn enclosure(self) -> SpectralEnclosureV1 {
        self.enclosure
    }
}

/// Why a well-formed separation proposition remains unresolved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnknownSeparationReasonV1 {
    /// Enclosures remain too wide to resolve overlap.
    EnclosureTooWide,
    /// Iteration/budget stopped before the requested bound was established.
    BudgetExhausted,
    /// Required residual or metric evidence is unavailable.
    MissingEvidence,
}

impl UnknownSeparationReasonV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::EnclosureTooWide => 0,
            Self::BudgetExhausted => 1,
            Self::MissingEvidence => 2,
        }
    }
}

/// Unvalidated internal-resolution claim carried by a cluster draft.
/// Favorable variants become inspectable as truth only through
/// [`ValidatedSpectralClusterV1::internal`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable internal-resolution variants intentionally retain their admitted evidence receipt inline"
)]
pub enum InternalClusterStateV1 {
    /// No internal-resolution proposition.
    NoClaim,
    /// Proposition is meaningful but unresolved.
    Unknown {
        /// Why the internal-resolution proposition could not be decided.
        reason: UnknownSeparationReasonV1,
    },
    /// Explicit no-claim state used when the producer declines to assign
    /// internal-separation semantics. This makes no proposition that
    /// separation is mathematically inapplicable; [`Self::Unknown`] instead
    /// says a chosen separation proposition is meaningful but unresolved.
    NoClaimUndefined,
    /// Exact algebraic multiplicity one makes internal separation vacuous.
    Simple,
    /// Repetition/zero internal separation is positively established.
    ProvenDegenerate {
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Members represented inside this cluster have a positive lower split.
    Resolved {
        /// Positive dimensionless lower bound.
        lower: PositiveFiniteV1,
        /// Norm/chart model identity.
        norm: SpectralNormId,
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// One unvalidated set-valued spectral cluster draft.
///
/// The draft exposes only neutral lineage and enclosure data. Claimed
/// localization authority, multiplicity, internal resolution, and defectivity
/// remain sealed until [`SpectralTruthV1::new`] returns a
/// [`ValidatedSpectralClusterV1`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralClusterV1 {
    id: SpectralClusterIdV1,
    localization: SpectralLocalizationV1,
    algebraic_multiplicity: MultiplicityClaimV1,
    geometric_multiplicity: MultiplicityClaimV1,
    internal: InternalClusterStateV1,
}

impl SpectralClusterV1 {
    /// Construct a cluster draft and reject locally impossible multiplicity
    /// statements. This does not validate any evidence binding.
    ///
    /// # Errors
    ///
    /// Returns a structured truth error when either multiplicity claim is
    /// malformed, geometric multiplicity can exceed algebraic multiplicity,
    /// or the declared internal state is impossible for the multiplicities.
    pub fn new(
        id: SpectralClusterIdV1,
        localization: SpectralLocalizationV1,
        algebraic_multiplicity: MultiplicityClaimV1,
        geometric_multiplicity: MultiplicityClaimV1,
        internal: InternalClusterStateV1,
    ) -> Result<Self, SpectralTruthErrorV1> {
        algebraic_multiplicity.validate()?;
        geometric_multiplicity.validate()?;
        if let (Some(geometric_min), Some(algebraic_max)) = (
            geometric_multiplicity.minimum(),
            algebraic_multiplicity.maximum(),
        ) && geometric_min > algebraic_max
        {
            return Err(SpectralTruthErrorV1::GeometricExceedsAlgebraic);
        }
        if let (Some(algebraic), Some(geometric)) = (
            algebraic_multiplicity.exact(),
            geometric_multiplicity.exact(),
        ) && geometric > algebraic
        {
            return Err(SpectralTruthErrorV1::GeometricExceedsAlgebraic);
        }
        match internal {
            InternalClusterStateV1::Simple if algebraic_multiplicity.exact() != Some(1) => {
                return Err(SpectralTruthErrorV1::InvalidInternalClusterState);
            }
            InternalClusterStateV1::ProvenDegenerate { .. }
            | InternalClusterStateV1::Resolved { .. }
                if algebraic_multiplicity
                    .minimum()
                    .is_none_or(|value| value < 2) =>
            {
                return Err(SpectralTruthErrorV1::InvalidInternalClusterState);
            }
            InternalClusterStateV1::NoClaim
            | InternalClusterStateV1::Unknown { .. }
            | InternalClusterStateV1::NoClaimUndefined
            | InternalClusterStateV1::Simple
            | InternalClusterStateV1::ProvenDegenerate { .. }
            | InternalClusterStateV1::Resolved { .. } => {}
        }
        Ok(Self {
            id,
            localization,
            algebraic_multiplicity,
            geometric_multiplicity,
            internal,
        })
    }

    /// Stable cluster lineage identity.
    #[must_use]
    pub const fn id(self) -> SpectralClusterIdV1 {
        self.id
    }

    /// Set-valued enclosure.
    #[must_use]
    pub const fn enclosure(self) -> SpectralEnclosureV1 {
        self.localization.enclosure
    }
}

/// Non-forgeable view of one cluster whose complete favorable evidence has
/// been checked against the enclosing admitted problem and result semantics.
/// Instances are minted only as part of [`SpectralTruthV1`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValidatedSpectralClusterV1 {
    cluster: SpectralClusterV1,
}

impl ValidatedSpectralClusterV1 {
    const fn from_validated(cluster: SpectralClusterV1) -> Self {
        Self { cluster }
    }

    /// Stable cluster lineage/membership identity.
    #[must_use]
    pub const fn id(&self) -> SpectralClusterIdV1 {
        self.cluster.id
    }

    /// Validated set-valued enclosure.
    #[must_use]
    pub const fn enclosure(&self) -> SpectralEnclosureV1 {
        self.cluster.localization.enclosure
    }

    /// Validated localization authority. This enum is not a total scientific
    /// evidence order.
    #[must_use]
    pub const fn localization_authority(&self) -> LocalizationAuthorityV1 {
        self.cluster.localization.authority
    }

    /// Proposition-bound localization witness, if the validated localization
    /// is favorable. Candidate localization returns `None`.
    #[must_use]
    pub const fn localization_witness(&self) -> Option<AdmittedSpectralWitnessV1> {
        self.cluster.localization.witness
    }

    /// Validated algebraic-multiplicity claim.
    #[must_use]
    pub const fn algebraic_multiplicity(&self) -> MultiplicityClaimV1 {
        self.cluster.algebraic_multiplicity
    }

    /// Validated geometric-multiplicity claim.
    #[must_use]
    pub const fn geometric_multiplicity(&self) -> MultiplicityClaimV1 {
        self.cluster.geometric_multiplicity
    }

    /// Validated per-cluster internal-resolution state.
    #[must_use]
    pub const fn internal(&self) -> InternalClusterStateV1 {
        self.cluster.internal
    }

    /// Infer defectivity only from two validated exact multiplicities, never
    /// from a repeated numerical value alone.
    #[must_use]
    pub const fn defectivity(&self) -> DefectivityStateV1 {
        match (
            self.cluster.algebraic_multiplicity.exact(),
            self.cluster.geometric_multiplicity.exact(),
        ) {
            (Some(algebraic), Some(geometric)) if geometric < algebraic => {
                DefectivityStateV1::ProvenDefective
            }
            (Some(algebraic), Some(geometric)) if geometric == algebraic => {
                DefectivityStateV1::ProvenSemisimple
            }
            _ => DefectivityStateV1::Unknown,
        }
    }
}

/// Finite strictly positive dimensionless value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PositiveFiniteV1(f64);

impl PositiveFiniteV1 {
    /// Validate and construct a positive finite value.
    ///
    /// # Errors
    ///
    /// Returns [`SpectralTruthErrorV1::InvalidSeparationBound`] when `value`
    /// is non-finite or not strictly positive.
    pub fn new(value: f64) -> Result<Self, SpectralTruthErrorV1> {
        if !(value.is_finite() && value > 0.0) {
            return Err(SpectralTruthErrorV1::InvalidSeparationBound);
        }
        Ok(Self(value))
    }

    /// Raw dimensionless value.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

/// Separation at the boundary of a partial ordered request.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable partial-boundary variants intentionally retain their admitted evidence receipt inline"
)]
pub enum PartialBoundaryStateV1 {
    /// Boundary remains unresolved.
    Unknown {
        /// Why separation at the requested prefix boundary is unresolved.
        reason: UnknownSeparationReasonV1,
    },
    /// Returned prefix is positively separated from the unreturned suffix.
    Separated {
        /// Positive lower bound on separation from the unreturned suffix.
        lower: PositiveFiniteV1,
        /// Norm/chart model in which `lower` is stated.
        norm: SpectralNormId,
        /// Evidence bound to this exact result-set boundary statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Request boundary cut a repeated cluster, which was returned whole.
    ClusterClosed {
        /// Repeated boundary cluster returned in its entirety.
        cluster: SpectralClusterIdV1,
        /// Evidence bound to this exact cluster-closure statement.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// Membership resolution at a named region boundary.
#[derive(Debug, Clone, PartialEq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable region-boundary variants intentionally retain their admitted evidence receipt inline"
)]
pub enum RegionBoundaryStateV1 {
    /// Boundary remains unresolved.
    Unknown {
        /// Why membership at the named region boundary is unresolved.
        reason: UnknownSeparationReasonV1,
    },
    /// Returned set is positively separated from the region boundary.
    Separated {
        /// Positive lower bound on separation from the region boundary.
        lower: PositiveFiniteV1,
        /// Norm/chart model in which `lower` is stated.
        norm: SpectralNormId,
        /// Evidence bound to this exact region-separation statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Boundary intersections were explicitly classified under the admitted
    /// open/closed policy.
    IntersectionsResolved {
        /// Boundary-intersecting clusters classified as included.
        included: Vec<SpectralClusterIdV1>,
        /// Exact algebraic multiplicity classified as excluded.
        excluded_algebraic: u32,
        /// Evidence bound to this exact intersection classification.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// Scope-boundary truth. Internal cluster resolution lives on each cluster;
/// this axis describes only the returned set versus its requested complement.
#[derive(Debug, Clone, PartialEq)]
#[allow(
    clippy::large_enum_variant,
    reason = "scope boundary variants preserve their evidence-bearing subordinate state inline"
)]
pub enum ScopeBoundaryStateV1 {
    /// No boundary proposition.
    NoClaim,
    /// Boundary of a partial ordered request.
    Partial(PartialBoundaryStateV1),
    /// Boundary of a named region request.
    Region(RegionBoundaryStateV1),
    /// Full-spectrum accounting proves there is no external complement.
    FullSpectrum,
}

/// Evidence-bearing authority of the returned result set.
///
/// This is deliberately not `Ord` and exposes no `min`/`max` lattice:
/// estimates, residual bounds, and rigorous enclosures are different
/// propositions and are incomparable unless their exact models are related by
/// an additional admitted theorem.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable result-authority variants intentionally retain their admitted evidence receipt inline"
)]
pub enum SpectralResultAuthorityV1 {
    /// No proposition is made.
    NoClaim,
    /// A proposition is well-defined but unresolved.
    Unknown,
    /// Raw candidate output without an error authority.
    Candidate,
    /// Estimated result with non-rigorous error information.
    Estimated {
        /// Evidence for the exact result-set estimate proposition.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Residual bounded in the named norm/model.
    ResidualBounded {
        /// Finite nonnegative residual upper bound.
        upper: f64,
        /// Norm/model identity.
        norm: SpectralNormId,
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Independently admitted enclosure of the exact result set.
    CertifiedEnclosure {
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// Status of an algebraic-cardinality prefix request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialCoverageStatusV1 {
    /// Fewer algebraic values were returned than requested.
    Incomplete,
    /// Exactly the requested algebraic cardinality was returned.
    Satisfied,
    /// The request cut through a repeated cluster, so that cluster was
    /// returned whole and the algebraic count legitimately exceeded `k`.
    ClusterClosureOverrun {
        /// Repeated cluster that straddled the requested cardinality boundary.
        boundary_cluster: SpectralClusterIdV1,
        /// Exact algebraic cardinality preceding `boundary_cluster`.
        preceding_algebraic: u32,
    },
}

/// Explicit projective/infinite accounting for a full finite-dimensional
/// spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::large_enum_variant,
    reason = "infinity-accounting variants intentionally retain their admitted evidence receipt inline"
)]
pub enum InfinityAccountingV1 {
    /// Ordinary problem with no projective infinity semantics.
    NotApplicable,
    /// Infinite multiplicity is returned as an optional projective cluster.
    Included {
        /// Exact algebraic multiplicity at projective infinity.
        algebraic: u32,
        /// Returned projective cluster, absent only when `algebraic` is zero.
        cluster: Option<SpectralClusterIdV1>,
        /// Evidence bound to this exact included-infinity statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Infinite multiplicity is omitted but exactly counted.
    ExcludedWithCount {
        /// Exact omitted algebraic multiplicity at projective infinity.
        algebraic: u32,
        /// Evidence bound to this exact excluded-infinity statement.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// Achieved coverage. Requested scope is taken only from the bound validated
/// problem; callers cannot repeat or silently replace it here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(
    clippy::large_enum_variant,
    reason = "favorable coverage variants intentionally retain their admitted evidence receipts inline"
)]
pub enum SpectralCoverageV1 {
    /// No result exists.
    NoResult,
    /// Candidate clusters with no completeness proposition.
    Candidates,
    /// Progress toward `CompletenessScopeV1::Partial` measured in algebraic
    /// cardinality, not cluster-record count.
    Partial {
        /// Exact algebraic cardinality represented by the returned clusters.
        returned_algebraic: u32,
        /// Relationship between that cardinality and the requested prefix.
        status: PartialCoverageStatusV1,
        /// Evidence bound to this exact partial-completeness statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Complete membership accounting for the exact admitted named region.
    RegionComplete {
        /// Exact algebraic cardinality inside the admitted named region.
        algebraic_cardinality: u32,
        /// Evidence bound to this exact region-completeness statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Complete finite/projective accounting for the admitted full scope.
    FullFinite {
        /// Exact algebraic cardinality represented by finite clusters.
        finite_algebraic: u32,
        /// Explicit accounting for the projective point at infinity.
        infinity: InfinityAccountingV1,
        /// Evidence bound to this exact full-completeness statement.
        witness: AdmittedSpectralWitnessV1,
    },
}

/// How the computation terminated; this is not an epistemic verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralTerminationV1 {
    /// No computation started.
    NotStarted,
    /// Declared stopping criterion completed.
    Completed,
    /// Budget exhausted with retained partial evidence.
    BudgetExhausted,
    /// Cancellation drained and finalized with retained partial evidence.
    Cancelled,
    /// Numerical failure. Independently validated candidates or partial
    /// evidence may still be retained; termination does not erase truth axes.
    NumericalFailure,
    /// Admission refused before execution.
    Refused,
}

/// Witness-free exact statement used to mint truth authority receipts. Each
/// variant owns one canonical payload dialect; callers cannot relabel a
/// localization witness as multiplicity, separation, or completeness.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectralTruthPropositionV1 {
    /// One cluster is localized with the stated authority and enclosure.
    ClusterLocalization {
        /// Stable lineage identity of the localized cluster.
        cluster: SpectralClusterIdV1,
        /// Epistemic class of the localization claim.
        authority: LocalizationAuthorityV1,
        /// Exact set-valued enclosure named by the claim.
        enclosure: SpectralEnclosureV1,
    },
    /// One algebraic or geometric multiplicity statement for a cluster.
    Multiplicity {
        /// Stable lineage identity of the cluster being counted.
        cluster: SpectralClusterIdV1,
        /// Exact cluster enclosure to which the count applies.
        enclosure: SpectralEnclosureV1,
        /// Whether algebraic or geometric multiplicity is asserted.
        kind: MultiplicityKindV1,
        /// Shape of the lower/bounded/exact assertion.
        assertion: MultiplicityAssertionV1,
        /// Positive lower bound, or exact value for `Exact`.
        lower: u32,
        /// Inclusive upper bound; absent only for a lower-bound-only claim.
        upper: Option<u32>,
    },
    /// Repetition or zero internal separation is established for one cluster.
    InternalDegeneracy {
        /// Stable lineage identity of the degenerate cluster.
        cluster: SpectralClusterIdV1,
        /// Exact cluster enclosure to which the assertion applies.
        enclosure: SpectralEnclosureV1,
        /// Algebraic multiplicity semantics at verification time.
        algebraic: MultiplicityStatementV1,
        /// Geometric multiplicity semantics at verification time.
        geometric: MultiplicityStatementV1,
    },
    /// Members represented within one repeated cluster have a positive split.
    InternalResolution {
        /// Stable lineage identity of the resolved cluster.
        cluster: SpectralClusterIdV1,
        /// Exact cluster enclosure to which the assertion applies.
        enclosure: SpectralEnclosureV1,
        /// Algebraic multiplicity semantics at verification time.
        algebraic: MultiplicityStatementV1,
        /// Geometric multiplicity semantics at verification time.
        geometric: MultiplicityStatementV1,
        /// Positive lower bound on internal separation.
        lower: PositiveFiniteV1,
        /// Norm/chart model in which `lower` is stated.
        norm: SpectralNormId,
    },
    /// The complete canonical result set has model-based estimate authority.
    ResultEstimate {
        /// Canonical result-set identity receiving the authority claim.
        result_set: SpectralResultSetIdV1,
    },
    /// The complete canonical result set has the stated residual upper bound.
    ResultResidualBound {
        /// Canonical result-set identity receiving the authority claim.
        result_set: SpectralResultSetIdV1,
        /// Finite nonnegative residual upper bound.
        upper: f64,
        /// Norm/model identity in which `upper` is stated.
        norm: SpectralNormId,
    },
    /// The complete canonical result set has rigorous enclosure authority.
    ResultCertifiedEnclosure {
        /// Canonical result-set identity receiving the authority claim.
        result_set: SpectralResultSetIdV1,
    },
    /// The result set accounts for the stated partial algebraic prefix.
    PartialCoverage {
        /// Canonical result-set identity receiving the completeness claim.
        result_set: SpectralResultSetIdV1,
        /// Exact algebraic cardinality represented by the result set.
        returned_algebraic: u32,
        /// Relationship between that cardinality and the requested prefix.
        status: PartialCoverageStatusV1,
    },
    /// The result set completely accounts for an admitted named region.
    RegionCompleteness {
        /// Canonical result-set identity receiving the completeness claim.
        result_set: SpectralResultSetIdV1,
        /// Exact algebraic cardinality in the named region.
        algebraic_cardinality: u32,
    },
    /// The result set completely accounts for the finite/projective spectrum.
    FullCompleteness {
        /// Canonical result-set identity receiving the completeness claim.
        result_set: SpectralResultSetIdV1,
        /// Exact algebraic cardinality represented by finite clusters.
        finite_algebraic: u32,
        /// Witness-free projective-infinity accounting embedded in the claim.
        infinity: InfinityAccountingStatementV1,
    },
    /// Projective-infinity multiplicity is included in the returned result.
    IncludedInfinity {
        /// Canonical result-set identity receiving the multiplicity claim.
        result_set: SpectralResultSetIdV1,
        /// Exact algebraic multiplicity at projective infinity.
        algebraic: u32,
        /// Returned projective cluster, absent only for zero multiplicity.
        cluster: Option<SpectralClusterIdV1>,
    },
    /// Projective-infinity multiplicity is omitted but exactly counted.
    ExcludedInfinity {
        /// Canonical result-set identity receiving the multiplicity claim.
        result_set: SpectralResultSetIdV1,
        /// Exact omitted algebraic multiplicity at projective infinity.
        algebraic: u32,
    },
    /// A partial result is positively separated from its unreturned suffix.
    PartialBoundarySeparated {
        /// Canonical result-set identity receiving the separation claim.
        result_set: SpectralResultSetIdV1,
        /// Positive lower bound on separation from the suffix.
        lower: PositiveFiniteV1,
        /// Norm/chart model in which `lower` is stated.
        norm: SpectralNormId,
    },
    /// A repeated cluster crossing the partial boundary was returned whole.
    PartialBoundaryClusterClosed {
        /// Canonical result-set identity receiving the closure claim.
        result_set: SpectralResultSetIdV1,
        /// Repeated boundary cluster returned in its entirety.
        cluster: SpectralClusterIdV1,
    },
    /// A region result is positively separated from the region boundary.
    RegionBoundarySeparated {
        /// Canonical result-set identity receiving the separation claim.
        result_set: SpectralResultSetIdV1,
        /// Positive lower bound on separation from the region boundary.
        lower: PositiveFiniteV1,
        /// Norm/chart model in which `lower` is stated.
        norm: SpectralNormId,
    },
    /// Intersections with a named region boundary are explicitly classified.
    RegionBoundaryIntersections {
        /// Canonical result-set identity receiving the classification claim.
        result_set: SpectralResultSetIdV1,
        /// Boundary-intersecting clusters classified as included.
        included: Vec<SpectralClusterIdV1>,
        /// Exact algebraic multiplicity classified as excluded.
        excluded_algebraic: u32,
    },
}

/// Witness-free infinity summary embedded into the whole-completeness
/// proposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InfinityAccountingStatementV1 {
    /// Ordinary problem with no projective-infinity semantics.
    NotApplicable,
    /// Infinite multiplicity is included in an optional projective cluster.
    Included {
        /// Exact algebraic multiplicity at projective infinity.
        algebraic: u32,
        /// Returned projective cluster, absent only for zero multiplicity.
        cluster: Option<SpectralClusterIdV1>,
    },
    /// Infinite multiplicity is omitted but exactly counted.
    ExcludedWithCount {
        /// Exact omitted algebraic multiplicity at projective infinity.
        algebraic: u32,
    },
}

impl InfinityAccountingV1 {
    const fn statement(self) -> InfinityAccountingStatementV1 {
        match self {
            Self::NotApplicable => InfinityAccountingStatementV1::NotApplicable,
            Self::Included {
                algebraic, cluster, ..
            } => InfinityAccountingStatementV1::Included { algebraic, cluster },
            Self::ExcludedWithCount { algebraic, .. } => {
                InfinityAccountingStatementV1::ExcludedWithCount { algebraic }
            }
        }
    }
}

fn push_result_set_id(out: &mut Vec<u8>, id: SpectralResultSetIdV1) {
    out.extend_from_slice(id.as_bytes());
}

fn push_optional_cluster(out: &mut Vec<u8>, cluster: Option<SpectralClusterIdV1>) {
    match cluster {
        Some(cluster) => {
            out.push(1);
            out.extend_from_slice(cluster.as_bytes());
        }
        None => out.push(0),
    }
}

fn push_partial_status(out: &mut Vec<u8>, status: PartialCoverageStatusV1) {
    match status {
        PartialCoverageStatusV1::Incomplete => out.push(0),
        PartialCoverageStatusV1::Satisfied => out.push(1),
        PartialCoverageStatusV1::ClusterClosureOverrun {
            boundary_cluster,
            preceding_algebraic,
        } => {
            out.push(2);
            out.extend_from_slice(boundary_cluster.as_bytes());
            push_u32(out, preceding_algebraic);
        }
    }
}

fn push_infinity_statement(out: &mut Vec<u8>, statement: InfinityAccountingStatementV1) {
    match statement {
        InfinityAccountingStatementV1::NotApplicable => out.push(0),
        InfinityAccountingStatementV1::Included { algebraic, cluster } => {
            out.push(1);
            push_u32(out, algebraic);
            push_optional_cluster(out, cluster);
        }
        InfinityAccountingStatementV1::ExcludedWithCount { algebraic } => {
            out.push(2);
            push_u32(out, algebraic);
        }
    }
}

/// Canonical proposition receipt for a truth claim under one exact validated
/// problem identity.
///
/// # Errors
///
/// Returns [`CanonicalError`] when canonical field, collection, or total-byte
/// limits are exceeded or a payload length cannot be represented.
#[allow(
    clippy::too_many_lines,
    reason = "the exhaustive match is the single canonical payload dialect table for all v1 truth propositions"
)]
pub fn truth_proposition_receipt(
    problem: SpectralProblemId,
    proposition: &SpectralTruthPropositionV1,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    let mut payload = Vec::with_capacity(256);
    payload.extend_from_slice(problem.as_bytes());
    let kind = match proposition {
        SpectralTruthPropositionV1::ClusterLocalization {
            cluster,
            authority,
            enclosure,
        } => {
            payload.push(0);
            payload.extend_from_slice(cluster.as_bytes());
            payload.push(match authority {
                LocalizationAuthorityV1::Candidate => 0,
                LocalizationAuthorityV1::Estimated => 1,
                LocalizationAuthorityV1::Enclosed => 2,
            });
            push_enclosure(&mut payload, *enclosure);
            SpectralPropositionKindV1::ResultAuthority
        }
        SpectralTruthPropositionV1::Multiplicity {
            cluster,
            enclosure,
            kind,
            assertion,
            lower,
            upper,
        } => {
            payload.push(0);
            payload.extend_from_slice(cluster.as_bytes());
            push_enclosure(&mut payload, *enclosure);
            payload.push(kind.tag());
            payload.push(assertion.tag());
            push_u32(&mut payload, *lower);
            match upper {
                Some(upper) => {
                    payload.push(1);
                    push_u32(&mut payload, *upper);
                }
                None => payload.push(0),
            }
            SpectralPropositionKindV1::Multiplicity
        }
        SpectralTruthPropositionV1::InternalDegeneracy {
            cluster,
            enclosure,
            algebraic,
            geometric,
        } => {
            payload.push(0);
            payload.extend_from_slice(cluster.as_bytes());
            push_enclosure(&mut payload, *enclosure);
            push_multiplicity_statement(&mut payload, *algebraic);
            push_multiplicity_statement(&mut payload, *geometric);
            SpectralPropositionKindV1::Separation
        }
        SpectralTruthPropositionV1::InternalResolution {
            cluster,
            enclosure,
            algebraic,
            geometric,
            lower,
            norm,
        } => {
            payload.push(1);
            payload.extend_from_slice(cluster.as_bytes());
            push_enclosure(&mut payload, *enclosure);
            push_multiplicity_statement(&mut payload, *algebraic);
            push_multiplicity_statement(&mut payload, *geometric);
            push_u64(&mut payload, canonical_f64_bits(lower.get()));
            payload.extend_from_slice(norm.as_bytes());
            SpectralPropositionKindV1::Separation
        }
        SpectralTruthPropositionV1::ResultEstimate { result_set } => {
            payload.push(1);
            push_result_set_id(&mut payload, *result_set);
            SpectralPropositionKindV1::ResultAuthority
        }
        SpectralTruthPropositionV1::ResultResidualBound {
            result_set,
            upper,
            norm,
        } => {
            if !upper.is_finite() {
                return Err(CanonicalError::NonFiniteFloat {
                    bits: upper.to_bits(),
                });
            }
            payload.push(2);
            push_result_set_id(&mut payload, *result_set);
            push_u64(&mut payload, canonical_f64_bits(*upper));
            payload.extend_from_slice(norm.as_bytes());
            SpectralPropositionKindV1::ResultAuthority
        }
        SpectralTruthPropositionV1::ResultCertifiedEnclosure { result_set } => {
            payload.push(3);
            push_result_set_id(&mut payload, *result_set);
            SpectralPropositionKindV1::ResultAuthority
        }
        SpectralTruthPropositionV1::PartialCoverage {
            result_set,
            returned_algebraic,
            status,
        } => {
            payload.push(0);
            push_result_set_id(&mut payload, *result_set);
            push_u32(&mut payload, *returned_algebraic);
            push_partial_status(&mut payload, *status);
            SpectralPropositionKindV1::Completeness
        }
        SpectralTruthPropositionV1::RegionCompleteness {
            result_set,
            algebraic_cardinality,
        } => {
            payload.push(1);
            push_result_set_id(&mut payload, *result_set);
            push_u32(&mut payload, *algebraic_cardinality);
            SpectralPropositionKindV1::Completeness
        }
        SpectralTruthPropositionV1::FullCompleteness {
            result_set,
            finite_algebraic,
            infinity,
        } => {
            payload.push(2);
            push_result_set_id(&mut payload, *result_set);
            push_u32(&mut payload, *finite_algebraic);
            push_infinity_statement(&mut payload, *infinity);
            SpectralPropositionKindV1::Completeness
        }
        SpectralTruthPropositionV1::IncludedInfinity {
            result_set,
            algebraic,
            cluster,
        } => {
            payload.push(1);
            push_result_set_id(&mut payload, *result_set);
            push_u32(&mut payload, *algebraic);
            push_optional_cluster(&mut payload, *cluster);
            SpectralPropositionKindV1::Multiplicity
        }
        SpectralTruthPropositionV1::ExcludedInfinity {
            result_set,
            algebraic,
        } => {
            payload.push(2);
            push_result_set_id(&mut payload, *result_set);
            push_u32(&mut payload, *algebraic);
            SpectralPropositionKindV1::Multiplicity
        }
        SpectralTruthPropositionV1::PartialBoundarySeparated {
            result_set,
            lower,
            norm,
        } => {
            payload.push(2);
            push_result_set_id(&mut payload, *result_set);
            push_u64(&mut payload, canonical_f64_bits(lower.get()));
            payload.extend_from_slice(norm.as_bytes());
            SpectralPropositionKindV1::Separation
        }
        SpectralTruthPropositionV1::RegionBoundarySeparated {
            result_set,
            lower,
            norm,
        } => {
            payload.push(4);
            push_result_set_id(&mut payload, *result_set);
            push_u64(&mut payload, canonical_f64_bits(lower.get()));
            payload.extend_from_slice(norm.as_bytes());
            SpectralPropositionKindV1::Separation
        }
        SpectralTruthPropositionV1::PartialBoundaryClusterClosed {
            result_set,
            cluster,
        } => {
            payload.push(3);
            push_result_set_id(&mut payload, *result_set);
            payload.extend_from_slice(cluster.as_bytes());
            SpectralPropositionKindV1::Separation
        }
        SpectralTruthPropositionV1::RegionBoundaryIntersections {
            result_set,
            included,
            excluded_algebraic,
        } => {
            if included.len() > MAX_REGION_BOUNDARY_REFERENCES_V1 {
                return Err(CanonicalError::LimitExceeded {
                    kind: LimitKind::CollectionItems,
                    requested: u64::try_from(included.len())
                        .map_err(|_| CanonicalError::LengthOverflow)?,
                    limit: u64::try_from(MAX_REGION_BOUNDARY_REFERENCES_V1)
                        .map_err(|_| CanonicalError::LengthOverflow)?,
                });
            }
            payload.push(5);
            push_result_set_id(&mut payload, *result_set);
            let mut included = included.clone();
            included.sort_unstable();
            push_u32(
                &mut payload,
                u32::try_from(included.len()).map_err(|_| CanonicalError::LengthOverflow)?,
            );
            for cluster in included {
                payload.extend_from_slice(cluster.as_bytes());
            }
            push_u32(&mut payload, *excluded_algebraic);
            SpectralPropositionKindV1::Separation
        }
    };
    spectral_proposition_receipt(kind, &payload)
}

fn push_cluster_semantics(out: &mut Vec<u8>, cluster: &SpectralClusterV1) {
    out.extend_from_slice(cluster.id.as_bytes());
    out.push(match cluster.localization.authority {
        LocalizationAuthorityV1::Candidate => 0,
        LocalizationAuthorityV1::Estimated => 1,
        LocalizationAuthorityV1::Enclosed => 2,
    });
    push_enclosure(out, cluster.localization.enclosure);
    cluster.algebraic_multiplicity.push_semantics(out);
    cluster.geometric_multiplicity.push_semantics(out);
    match cluster.internal {
        InternalClusterStateV1::NoClaim => out.push(0),
        InternalClusterStateV1::Unknown { reason } => {
            out.push(1);
            out.push(reason.tag());
        }
        InternalClusterStateV1::NoClaimUndefined => out.push(2),
        InternalClusterStateV1::Simple => out.push(3),
        InternalClusterStateV1::ProvenDegenerate { .. } => out.push(4),
        InternalClusterStateV1::Resolved { lower, norm, .. } => {
            out.push(5);
            push_u64(out, canonical_f64_bits(lower.get()));
            out.extend_from_slice(norm.as_bytes());
        }
    }
}

fn cluster_semantics(cluster: &SpectralClusterV1) -> Vec<u8> {
    let mut payload = Vec::with_capacity(192);
    push_cluster_semantics(&mut payload, cluster);
    payload
}

/// Canonical result-set receipt, independent of evidence order and authority
/// anchors. Favorable evidence later binds this identity.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the cluster count or canonical payload
/// exceeds the v1 identity envelope, or when a length cannot be represented.
pub fn spectral_result_set_receipt(
    clusters: &[SpectralClusterV1],
) -> Result<IdentityReceipt<SpectralResultSetIdV1>, CanonicalError> {
    if clusters.len() > MAX_SPECTRAL_CLUSTERS_V1 {
        return Err(CanonicalError::LimitExceeded {
            kind: LimitKind::CollectionItems,
            requested: u64::try_from(clusters.len()).map_err(|_| CanonicalError::LengthOverflow)?,
            limit: u64::try_from(MAX_SPECTRAL_CLUSTERS_V1)
                .map_err(|_| CanonicalError::LengthOverflow)?,
        });
    }
    let mut payloads: Vec<Vec<u8>> = clusters.iter().map(cluster_semantics).collect();
    payloads.sort();
    CanonicalEncoder::<SpectralResultSetIdV1, _>::new(TRUTH_IDENTITY_LIMITS, NeverCancel)?
        .canonical_set(
            Field::new(0, "clusters"),
            u64::try_from(payloads.len()).map_err(|_| CanonicalError::LengthOverflow)?,
            payloads.iter().map(Vec::as_slice),
        )?
        .finish()
}

/// Raw result draft. It cannot mint truth until validated against the complete
/// admitted problem descriptor and every favorable witness receipt.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralTruthDraftV1 {
    authority: SpectralResultAuthorityV1,
    coverage: SpectralCoverageV1,
    clusters: Vec<SpectralClusterV1>,
    boundary: ScopeBoundaryStateV1,
    termination: SpectralTerminationV1,
}

impl SpectralTruthDraftV1 {
    /// Assemble an untrusted truth draft for validation against an admitted
    /// problem; construction alone grants no authority.
    #[must_use]
    pub fn new(
        authority: SpectralResultAuthorityV1,
        coverage: SpectralCoverageV1,
        clusters: Vec<SpectralClusterV1>,
        boundary: ScopeBoundaryStateV1,
        termination: SpectralTerminationV1,
    ) -> Self {
        Self {
            authority,
            coverage,
            clusters,
            boundary,
            termination,
        }
    }
}

/// Validated product truth bound to one full problem descriptor and canonical
/// result-set identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralTruthV1 {
    problem_id: SpectralProblemId,
    result_set_id: SpectralResultSetIdV1,
    authority: SpectralResultAuthorityV1,
    coverage: SpectralCoverageV1,
    clusters: Vec<ValidatedSpectralClusterV1>,
    boundary: ScopeBoundaryStateV1,
    termination: SpectralTerminationV1,
}

impl SpectralTruthV1 {
    /// Validate a draft against the full admitted problem.
    ///
    /// # Errors
    ///
    /// Returns [`SpectralTruthReportV1`] containing every deterministically
    /// ranked cross-axis, evidence, cardinality, boundary, and identity
    /// violation found in the bounded draft.
    pub fn new(
        problem: &ValidatedSpectralProblemV1,
        draft: SpectralTruthDraftV1,
    ) -> Result<Self, SpectralTruthReportV1> {
        validate_truth_v1(problem, draft)
    }

    /// Full admitted problem identity to which this truth is bound.
    #[must_use]
    pub const fn problem_id(&self) -> SpectralProblemId {
        self.problem_id
    }

    /// Canonical identity of the validated set-valued cluster collection.
    #[must_use]
    pub const fn result_set_id(&self) -> SpectralResultSetIdV1 {
        self.result_set_id
    }

    /// Evidence-bearing authority of the complete returned result set.
    #[must_use]
    pub const fn authority(&self) -> SpectralResultAuthorityV1 {
        self.authority
    }

    /// Achieved algebraic coverage relative to the admitted request.
    #[must_use]
    pub const fn coverage(&self) -> SpectralCoverageV1 {
        self.coverage
    }

    /// Canonically ordered set-valued spectral clusters.
    #[must_use]
    pub fn clusters(&self) -> &[ValidatedSpectralClusterV1] {
        &self.clusters
    }

    /// Validated state of the returned set versus its requested complement.
    #[must_use]
    pub const fn boundary(&self) -> &ScopeBoundaryStateV1 {
        &self.boundary
    }

    /// Computation termination state, orthogonal to epistemic authority.
    #[must_use]
    pub const fn termination(&self) -> SpectralTerminationV1 {
        self.termination
    }
}

fn validate_truth_witness(
    witness: AdmittedSpectralWitnessV1,
    problem: SpectralProblemId,
    proposition: SpectralTruthPropositionV1,
    issues: &mut Vec<SpectralTruthErrorV1>,
) {
    match truth_proposition_receipt(problem, &proposition) {
        Ok(expected) if witness.proposition() != expected.id() => {
            issues.push(SpectralTruthErrorV1::WitnessPropositionMismatch {
                expected: expected.id(),
                found: witness.proposition(),
            });
        }
        Ok(expected) if !witness.matches_receipt(expected) => {
            issues.push(SpectralTruthErrorV1::WitnessObservationMismatch {
                proposition: expected.id(),
            });
        }
        Ok(_) => {}
        Err(error) => issues.push(SpectralTruthErrorV1::Identity(error)),
    }
}

fn validate_multiplicity_witness(
    problem: SpectralProblemId,
    cluster: SpectralClusterIdV1,
    enclosure: SpectralEnclosureV1,
    kind: MultiplicityKindV1,
    claim: MultiplicityClaimV1,
    issues: &mut Vec<SpectralTruthErrorV1>,
) {
    let Some(witness) = claim.witness() else {
        return;
    };
    let (assertion, lower, upper) = match claim {
        MultiplicityClaimV1::Unknown => return,
        MultiplicityClaimV1::LowerBound { value, .. } => {
            (MultiplicityAssertionV1::LowerBound, value, None)
        }
        MultiplicityClaimV1::Bounds { lower, upper, .. } => {
            (MultiplicityAssertionV1::Bounds, lower, Some(upper))
        }
        MultiplicityClaimV1::Exact { value, .. } => {
            (MultiplicityAssertionV1::Exact, value, Some(value))
        }
    };
    validate_truth_witness(
        witness,
        problem,
        SpectralTruthPropositionV1::Multiplicity {
            cluster,
            enclosure,
            kind,
            assertion,
            lower,
            upper,
        },
        issues,
    );
}

fn exact_algebraic_sum(clusters: &[SpectralClusterV1]) -> Option<u32> {
    clusters.iter().try_fold(0_u32, |sum, cluster| {
        sum.checked_add(cluster.algebraic_multiplicity.exact()?)
    })
}

fn inferred_algebraic_minimum_sum(clusters: &[SpectralClusterV1]) -> Option<u64> {
    clusters.iter().try_fold(0_u64, |sum, cluster| {
        let algebraic = cluster.algebraic_multiplicity.minimum().unwrap_or(0);
        let geometric = cluster.geometric_multiplicity.minimum().unwrap_or(0);
        // Geometric multiplicity cannot exceed algebraic multiplicity, so an
        // admitted geometric lower bound is also an algebraic lower bound even
        // when the algebraic axis itself is reported as unknown.
        sum.checked_add(u64::from(algebraic.max(geometric)))
    })
}

fn cluster_exists(clusters: &[SpectralClusterV1], id: SpectralClusterIdV1) -> bool {
    clusters
        .binary_search_by_key(&id, |cluster| cluster.id)
        .is_ok()
}

fn projective_clusters(clusters: &[SpectralClusterV1]) -> Vec<SpectralClusterV1> {
    clusters
        .iter()
        .copied()
        .filter(|cluster| {
            matches!(
                cluster.localization.enclosure,
                SpectralEnclosureV1::ProjectiveInfinity
            )
        })
        .collect()
}

fn projective_cluster_has_favorable_claim(cluster: &SpectralClusterV1) -> bool {
    cluster.localization.authority != LocalizationAuthorityV1::Candidate
        || !matches!(cluster.algebraic_multiplicity, MultiplicityClaimV1::Unknown)
        || !matches!(cluster.geometric_multiplicity, MultiplicityClaimV1::Unknown)
        || matches!(
            cluster.internal,
            InternalClusterStateV1::Simple
                | InternalClusterStateV1::ProvenDegenerate { .. }
                | InternalClusterStateV1::Resolved { .. }
        )
}

fn enclosure_can_contain_real_spectrum(enclosure: SpectralEnclosureV1) -> bool {
    match enclosure {
        SpectralEnclosureV1::Real(_) | SpectralEnclosureV1::ProjectiveInfinity => true,
        SpectralEnclosureV1::ComplexBox { imag, .. } => imag.lower() <= 0.0 && imag.upper() >= 0.0,
    }
}

fn validate_cluster_evidence(
    problem: SpectralProblemId,
    clusters: &[SpectralClusterV1],
    issues: &mut Vec<SpectralTruthErrorV1>,
) {
    for cluster in clusters {
        match (cluster.localization.authority, cluster.localization.witness) {
            (LocalizationAuthorityV1::Candidate, None) => {}
            (LocalizationAuthorityV1::Candidate, Some(_))
            | (LocalizationAuthorityV1::Estimated | LocalizationAuthorityV1::Enclosed, None) => {
                issues.push(SpectralTruthErrorV1::AuthorityEvidenceMismatch);
            }
            (authority, Some(witness)) => validate_truth_witness(
                witness,
                problem,
                SpectralTruthPropositionV1::ClusterLocalization {
                    cluster: cluster.id,
                    authority,
                    enclosure: cluster.localization.enclosure,
                },
                issues,
            ),
        }
        validate_multiplicity_witness(
            problem,
            cluster.id,
            cluster.localization.enclosure,
            MultiplicityKindV1::Algebraic,
            cluster.algebraic_multiplicity,
            issues,
        );
        validate_multiplicity_witness(
            problem,
            cluster.id,
            cluster.localization.enclosure,
            MultiplicityKindV1::Geometric,
            cluster.geometric_multiplicity,
            issues,
        );
        match cluster.internal {
            InternalClusterStateV1::ProvenDegenerate { witness } => validate_truth_witness(
                witness,
                problem,
                SpectralTruthPropositionV1::InternalDegeneracy {
                    cluster: cluster.id,
                    enclosure: cluster.localization.enclosure,
                    algebraic: cluster.algebraic_multiplicity.statement(),
                    geometric: cluster.geometric_multiplicity.statement(),
                },
                issues,
            ),
            InternalClusterStateV1::Resolved {
                lower,
                norm,
                witness,
            } => validate_truth_witness(
                witness,
                problem,
                SpectralTruthPropositionV1::InternalResolution {
                    cluster: cluster.id,
                    enclosure: cluster.localization.enclosure,
                    algebraic: cluster.algebraic_multiplicity.statement(),
                    geometric: cluster.geometric_multiplicity.statement(),
                    lower,
                    norm,
                },
                issues,
            ),
            InternalClusterStateV1::NoClaim
            | InternalClusterStateV1::Unknown { .. }
            | InternalClusterStateV1::NoClaimUndefined
            | InternalClusterStateV1::Simple => {}
        }
    }
}

/// Validate every truth axis, canonicalize ordering, and mint no partial token
/// on any failure.
///
/// # Errors
///
/// Returns [`SpectralTruthReportV1`] containing every deterministically ranked
/// violation found after resource-cap admission. No validated cluster or truth
/// token is returned on failure.
#[allow(
    clippy::too_many_lines,
    reason = "the validator deliberately accumulates all bounded cross-axis v1 diagnostics before minting truth"
)]
pub fn validate_truth_v1(
    problem: &ValidatedSpectralProblemV1,
    mut draft: SpectralTruthDraftV1,
) -> Result<SpectralTruthV1, SpectralTruthReportV1> {
    let mut issues = Vec::new();
    if draft.clusters.len() > MAX_SPECTRAL_CLUSTERS_V1 {
        issues.push(SpectralTruthErrorV1::TooManyClusters {
            found: draft.clusters.len(),
            limit: MAX_SPECTRAL_CLUSTERS_V1,
        });
    }
    if let ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
        included,
        ..
    }) = &draft.boundary
        && included.len() > MAX_REGION_BOUNDARY_REFERENCES_V1
    {
        issues.push(SpectralTruthErrorV1::TooManyBoundaryReferences {
            found: included.len(),
            limit: MAX_REGION_BOUNDARY_REFERENCES_V1,
        });
    }
    if !issues.is_empty() {
        return Err(SpectralTruthReportV1::new(issues));
    }
    if let SpectralResultAuthorityV1::ResidualBounded { upper, .. } = &mut draft.authority
        && *upper == 0.0
    {
        *upper = 0.0;
    }
    // Cluster IDs are membership identities, but malformed drafts may repeat
    // one ID with conflicting semantics. The full evidence-free semantic key
    // makes every later `.find` and diagnostic independent of caller order.
    draft.clusters.sort_by_cached_key(cluster_semantics);
    if draft
        .clusters
        .windows(2)
        .any(|pair| pair[0].id == pair[1].id)
    {
        issues.push(SpectralTruthErrorV1::DuplicateCluster);
    }
    if let ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
        included,
        ..
    }) = &mut draft.boundary
    {
        included.sort_unstable();
        if included.windows(2).any(|pair| pair[0] == pair[1]) {
            issues.push(SpectralTruthErrorV1::DuplicateSeparationReference);
        }
        if included
            .iter()
            .any(|id| !cluster_exists(&draft.clusters, *id))
        {
            issues.push(SpectralTruthErrorV1::DanglingClusterReference);
        }
    }
    let result_receipt = spectral_result_set_receipt(&draft.clusters);
    let result_set_id = match result_receipt {
        Ok(receipt) => receipt.id(),
        Err(error) => {
            issues.push(SpectralTruthErrorV1::Identity(error));
            return Err(SpectralTruthReportV1::new(issues));
        }
    };
    let problem_id = problem.problem_id();
    validate_cluster_evidence(problem_id, &draft.clusters, &mut issues);
    if problem.requires_real_spectrum_truth()
        || matches!(
            problem.spec().ordering(),
            SpectralOrderingV1::RealAscending | SpectralOrderingV1::RealDescending
        )
    {
        let whole_set_enclosed = matches!(
            draft.authority,
            SpectralResultAuthorityV1::CertifiedEnclosure { .. }
        );
        if draft.clusters.iter().any(|cluster| {
            (whole_set_enclosed
                || cluster.localization.authority == LocalizationAuthorityV1::Enclosed)
                && !enclosure_can_contain_real_spectrum(cluster.localization.enclosure)
        }) {
            issues.push(SpectralTruthErrorV1::RealSpectrumEnclosureConflict);
        }
    }
    let projective = projective_clusters(&draft.clusters);
    if !projective.is_empty()
        && !matches!(
            problem.spec().class().descriptor(),
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective
            }
        )
    {
        issues.push(SpectralTruthErrorV1::ProjectiveClusterNotAdmitted);
    }
    if projective.len() > 1 {
        issues.push(SpectralTruthErrorV1::MultipleProjectiveClusters);
    }
    let whole_result_favorable = matches!(
        draft.authority,
        SpectralResultAuthorityV1::Estimated { .. }
            | SpectralResultAuthorityV1::ResidualBounded { .. }
            | SpectralResultAuthorityV1::CertifiedEnclosure { .. }
    );
    let coverage_asserts_membership = matches!(
        draft.coverage,
        SpectralCoverageV1::Partial { .. }
            | SpectralCoverageV1::RegionComplete { .. }
            | SpectralCoverageV1::FullFinite { .. }
    );
    if problem.projective_infinity_is_excluded()
        && !projective.is_empty()
        && (whole_result_favorable
            || coverage_asserts_membership
            || projective
                .iter()
                .any(projective_cluster_has_favorable_claim))
    {
        issues.push(SpectralTruthErrorV1::ProjectiveInfinityExcludedByRegularity);
    }

    match draft.authority {
        SpectralResultAuthorityV1::NoClaim
        | SpectralResultAuthorityV1::Unknown
        | SpectralResultAuthorityV1::Candidate => {}
        SpectralResultAuthorityV1::Estimated { witness } => validate_truth_witness(
            witness,
            problem_id,
            SpectralTruthPropositionV1::ResultEstimate {
                result_set: result_set_id,
            },
            &mut issues,
        ),
        SpectralResultAuthorityV1::ResidualBounded {
            upper,
            norm,
            witness,
        } => {
            if !(upper.is_finite() && upper >= 0.0) {
                issues.push(SpectralTruthErrorV1::InvalidResidualBound);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::ResultResidualBound {
                    result_set: result_set_id,
                    upper,
                    norm,
                },
                &mut issues,
            );
        }
        SpectralResultAuthorityV1::CertifiedEnclosure { witness } => validate_truth_witness(
            witness,
            problem_id,
            SpectralTruthPropositionV1::ResultCertifiedEnclosure {
                result_set: result_set_id,
            },
            &mut issues,
        ),
    }

    let exact_sum = exact_algebraic_sum(&draft.clusters);
    if let (Some(returned), Some(total)) = (exact_sum, problem.known_algebraic_cardinality())
        && returned > total
    {
        issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
    }
    if let (Some(returned_minimum), Some(total)) = (
        inferred_algebraic_minimum_sum(&draft.clusters),
        problem.known_algebraic_cardinality(),
    ) && returned_minimum > u64::from(total)
    {
        issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
    }
    if let (
        ScopeBoundaryStateV1::Region(RegionBoundaryStateV1::IntersectionsResolved {
            excluded_algebraic,
            ..
        }),
        Some(total),
    ) = (&draft.boundary, problem.known_algebraic_cardinality())
    {
        let total = u64::from(total);
        let excluded = u64::from(*excluded_algebraic);
        let returned_minimum = inferred_algebraic_minimum_sum(&draft.clusters);
        if excluded > total
            || returned_minimum
                .and_then(|returned| returned.checked_add(excluded))
                .is_none_or(|accounted| accounted > total)
        {
            issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
        }
        if let SpectralCoverageV1::RegionComplete {
            algebraic_cardinality,
            ..
        } = draft.coverage
            && u64::from(algebraic_cardinality)
                .checked_add(excluded)
                .is_none_or(|accounted| accounted > total)
        {
            issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
        }
    }
    // Metrics name the operator spaces after any admitted reduction. Gauge
    // nullity is pre-reduction provenance and must not be subtracted again.
    let operator_dimension = problem.spec().spaces().domain().dimension();
    for cluster in &draft.clusters {
        if let Some(minimum) = cluster.geometric_multiplicity.minimum()
            && minimum > operator_dimension
        {
            issues.push(SpectralTruthErrorV1::GeometricCapacityExceeded {
                cluster: cluster.id,
                minimum,
                dimension: operator_dimension,
            });
        }
    }
    let requested_scope = problem.spec().requested_scope();
    let complete = matches!(
        draft.coverage,
        SpectralCoverageV1::RegionComplete { .. } | SpectralCoverageV1::FullFinite { .. }
    );
    match draft.coverage {
        SpectralCoverageV1::NoResult => {
            if !draft.clusters.is_empty() {
                issues.push(SpectralTruthErrorV1::NoResultHasClusters);
            }
            if !matches!(draft.authority, SpectralResultAuthorityV1::NoClaim)
                || !matches!(&draft.boundary, ScopeBoundaryStateV1::NoClaim)
            {
                issues.push(SpectralTruthErrorV1::NoResultHasClaims);
            }
        }
        SpectralCoverageV1::Candidates => {}
        SpectralCoverageV1::Partial {
            returned_algebraic,
            status,
            witness,
        } => {
            if problem.known_algebraic_cardinality().is_none() {
                issues.push(SpectralTruthErrorV1::DiscreteSpectrumRegularityNotEstablished);
            }
            if let CompletenessScopeV1::Partial { requested } = requested_scope {
                if exact_sum != Some(returned_algebraic) {
                    issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
                }
                match status {
                    PartialCoverageStatusV1::Incomplete
                        if returned_algebraic > 0
                            && returned_algebraic < requested
                            && !matches!(
                                &draft.boundary,
                                ScopeBoundaryStateV1::Partial(
                                    PartialBoundaryStateV1::ClusterClosed { .. }
                                )
                            ) => {}
                    PartialCoverageStatusV1::Satisfied
                        if returned_algebraic == requested
                            && matches!(
                                &draft.boundary,
                                ScopeBoundaryStateV1::Partial(
                                    PartialBoundaryStateV1::Separated { .. }
                                )
                            ) => {}
                    PartialCoverageStatusV1::ClusterClosureOverrun {
                        boundary_cluster,
                        preceding_algebraic,
                    } if preceding_algebraic < requested
                        && returned_algebraic > requested
                        && cluster_exists(&draft.clusters, boundary_cluster)
                        && draft
                            .clusters
                            .iter()
                            .find(|cluster| cluster.id == boundary_cluster)
                            .and_then(|cluster| cluster.algebraic_multiplicity.exact())
                            .and_then(|multiplicity| {
                                preceding_algebraic.checked_add(multiplicity)
                            })
                            == Some(returned_algebraic)
                        && matches!(
                            &draft.boundary,
                            ScopeBoundaryStateV1::Partial(
                                PartialBoundaryStateV1::ClusterClosed { cluster, .. }
                            ) if *cluster == boundary_cluster
                        ) => {}
                    PartialCoverageStatusV1::Incomplete
                    | PartialCoverageStatusV1::Satisfied
                    | PartialCoverageStatusV1::ClusterClosureOverrun { .. } => {
                        issues.push(SpectralTruthErrorV1::InvalidPartialCoverage);
                    }
                }
            } else {
                issues.push(SpectralTruthErrorV1::CoverageScopeMismatch);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::PartialCoverage {
                    result_set: result_set_id,
                    returned_algebraic,
                    status,
                },
                &mut issues,
            );
        }
        SpectralCoverageV1::RegionComplete {
            algebraic_cardinality,
            witness,
        } => {
            if problem.known_algebraic_cardinality().is_none() {
                issues.push(SpectralTruthErrorV1::DiscreteSpectrumRegularityNotEstablished);
            }
            if !matches!(requested_scope, CompletenessScopeV1::Region { .. }) {
                issues.push(SpectralTruthErrorV1::CoverageScopeMismatch);
            }
            if exact_sum != Some(algebraic_cardinality)
                || (algebraic_cardinality == 0) != draft.clusters.is_empty()
            {
                issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
            }
            let boundary_compatible = if let ScopeBoundaryStateV1::Region(state) = &draft.boundary {
                region_boundary_state_is_compatible(requested_scope, state, false)
            } else {
                false
            };
            if !boundary_compatible {
                issues.push(SpectralTruthErrorV1::BoundaryScopeMismatch);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::RegionCompleteness {
                    result_set: result_set_id,
                    algebraic_cardinality,
                },
                &mut issues,
            );
        }
        SpectralCoverageV1::FullFinite {
            finite_algebraic,
            infinity,
            witness,
        } => {
            if let CompletenessScopeV1::FullFinite {
                algebraic_cardinality,
                infinity_policy,
            } = requested_scope
            {
                validate_full_accounting(
                    problem,
                    problem_id,
                    result_set_id,
                    &draft.clusters,
                    finite_algebraic,
                    infinity,
                    algebraic_cardinality,
                    infinity_policy,
                    &mut issues,
                );
                if problem.known_algebraic_cardinality() != Some(algebraic_cardinality) {
                    issues.push(
                        SpectralTruthErrorV1::FullCompletenessRegularityNotEstablished {
                            requested: algebraic_cardinality,
                            established: problem.known_algebraic_cardinality(),
                        },
                    );
                }
                if !matches!(&draft.boundary, ScopeBoundaryStateV1::FullSpectrum) {
                    issues.push(SpectralTruthErrorV1::BoundaryCoverageMismatch);
                }
            } else {
                issues.push(SpectralTruthErrorV1::CoverageScopeMismatch);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::FullCompleteness {
                    result_set: result_set_id,
                    finite_algebraic,
                    infinity: infinity.statement(),
                },
                &mut issues,
            );
        }
    }

    if let ScopeBoundaryStateV1::Partial(PartialBoundaryStateV1::ClusterClosed {
        cluster, ..
    }) = &draft.boundary
    {
        let exact_repeated = draft
            .clusters
            .iter()
            .find(|candidate| candidate.id == *cluster)
            .and_then(|candidate| candidate.algebraic_multiplicity.exact())
            .is_some_and(|multiplicity| multiplicity >= 2);
        let matching_overrun = matches!(
            draft.coverage,
            SpectralCoverageV1::Partial {
                status: PartialCoverageStatusV1::ClusterClosureOverrun {
                    boundary_cluster,
                    ..
                },
                ..
            } if boundary_cluster == *cluster
        );
        if !(exact_repeated && matching_overrun) {
            issues.push(SpectralTruthErrorV1::BoundaryCoverageMismatch);
        }
    }

    validate_boundary_evidence(
        problem_id,
        result_set_id,
        requested_scope,
        &draft.clusters,
        &draft.boundary,
        &mut issues,
    );
    if matches!(&draft.boundary, ScopeBoundaryStateV1::FullSpectrum)
        && !matches!(draft.coverage, SpectralCoverageV1::FullFinite { .. })
    {
        issues.push(SpectralTruthErrorV1::BoundaryCoverageMismatch);
    }

    if matches!(
        draft.termination,
        SpectralTerminationV1::NotStarted | SpectralTerminationV1::Refused
    ) && (!matches!(draft.authority, SpectralResultAuthorityV1::NoClaim)
        || !matches!(draft.coverage, SpectralCoverageV1::NoResult)
        || !draft.clusters.is_empty()
        || !matches!(&draft.boundary, ScopeBoundaryStateV1::NoClaim))
    {
        issues.push(SpectralTruthErrorV1::InactiveComputationHasClaims);
    }
    if complete && draft.termination != SpectralTerminationV1::Completed {
        issues.push(SpectralTruthErrorV1::CompleteCoverageRequiresCompletion);
    }
    finish_truth(problem_id, result_set_id, draft, issues)
}

fn finish_truth(
    problem_id: SpectralProblemId,
    result_set_id: SpectralResultSetIdV1,
    draft: SpectralTruthDraftV1,
    issues: Vec<SpectralTruthErrorV1>,
) -> Result<SpectralTruthV1, SpectralTruthReportV1> {
    if issues.is_empty() {
        let SpectralTruthDraftV1 {
            authority,
            coverage,
            clusters,
            boundary,
            termination,
        } = draft;
        Ok(SpectralTruthV1 {
            problem_id,
            result_set_id,
            authority,
            coverage,
            clusters: clusters
                .into_iter()
                .map(ValidatedSpectralClusterV1::from_validated)
                .collect(),
            boundary,
            termination,
        })
    } else {
        Err(SpectralTruthReportV1::new(issues))
    }
}

fn region_boundary_state_is_compatible(
    requested_scope: CompletenessScopeV1,
    state: &RegionBoundaryStateV1,
    allow_unknown: bool,
) -> bool {
    let CompletenessScopeV1::Region { boundary, .. } = requested_scope else {
        return false;
    };
    match state {
        RegionBoundaryStateV1::Unknown { .. } => allow_unknown,
        RegionBoundaryStateV1::Separated { .. } => true,
        RegionBoundaryStateV1::IntersectionsResolved {
            included,
            excluded_algebraic,
            ..
        } => match boundary {
            RegionBoundaryPolicyV1::Closed => *excluded_algebraic == 0,
            RegionBoundaryPolicyV1::Open => included.is_empty(),
            RegionBoundaryPolicyV1::RefuseIntersection => false,
        },
    }
}

#[allow(
    clippy::too_many_arguments,
    clippy::too_many_lines,
    reason = "full-spectrum accounting keeps all bound identities, cardinalities, policy, and diagnostic sink explicit"
)]
fn validate_full_accounting(
    problem: &ValidatedSpectralProblemV1,
    problem_id: SpectralProblemId,
    result_set_id: SpectralResultSetIdV1,
    clusters: &[SpectralClusterV1],
    finite_algebraic: u32,
    infinity: InfinityAccountingV1,
    total_algebraic: u32,
    policy: InfiniteEigenvaluePolicyV1,
    issues: &mut Vec<SpectralTruthErrorV1>,
) {
    let projective = projective_clusters(clusters);
    let finite_sum = clusters
        .iter()
        .filter(|cluster| {
            !matches!(
                cluster.localization.enclosure,
                SpectralEnclosureV1::ProjectiveInfinity
            )
        })
        .try_fold(0_u32, |sum, cluster| {
            sum.checked_add(cluster.algebraic_multiplicity.exact()?)
        });
    if finite_sum != Some(finite_algebraic) {
        issues.push(SpectralTruthErrorV1::CoverageCardinalityMismatch);
    }
    match infinity {
        InfinityAccountingV1::NotApplicable => {
            if policy != InfiniteEigenvaluePolicyV1::NoClaim
                || !matches!(
                    problem.spec().class().descriptor(),
                    DescriptorRoleV1::Ordinary
                )
                || !projective.is_empty()
                || finite_algebraic != total_algebraic
            {
                issues.push(SpectralTruthErrorV1::InfinityPolicyMismatch);
            }
        }
        InfinityAccountingV1::Included {
            algebraic,
            cluster,
            witness,
        } => {
            if policy != InfiniteEigenvaluePolicyV1::IncludeProjective
                || !matches!(
                    problem.spec().class().descriptor(),
                    DescriptorRoleV1::Descriptor {
                        infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective
                    }
                )
                || finite_algebraic.checked_add(algebraic) != Some(total_algebraic)
            {
                issues.push(SpectralTruthErrorV1::InfinityPolicyMismatch);
            }
            let projective_ok = match (algebraic, cluster, projective.as_slice()) {
                (0, None, []) => true,
                (value, Some(expected), [actual]) if value > 0 && actual.id == expected => {
                    actual.algebraic_multiplicity.exact() == Some(value)
                }
                _ => false,
            };
            if !projective_ok {
                issues.push(SpectralTruthErrorV1::InfinityAccountingMismatch);
            }
            if problem.projective_infinity_is_excluded() && algebraic != 0 {
                issues.push(SpectralTruthErrorV1::InfinityAccountingMismatch);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::IncludedInfinity {
                    result_set: result_set_id,
                    algebraic,
                    cluster,
                },
                issues,
            );
        }
        InfinityAccountingV1::ExcludedWithCount { algebraic, witness } => {
            if policy != InfiniteEigenvaluePolicyV1::ExcludeWithCount
                || !matches!(
                    problem.spec().class().descriptor(),
                    DescriptorRoleV1::Descriptor {
                        infinity_policy: InfiniteEigenvaluePolicyV1::ExcludeWithCount
                    }
                )
                || !projective.is_empty()
                || finite_algebraic.checked_add(algebraic) != Some(total_algebraic)
            {
                issues.push(SpectralTruthErrorV1::InfinityPolicyMismatch);
            }
            if problem.projective_infinity_is_excluded() && algebraic != 0 {
                issues.push(SpectralTruthErrorV1::InfinityAccountingMismatch);
            }
            validate_truth_witness(
                witness,
                problem_id,
                SpectralTruthPropositionV1::ExcludedInfinity {
                    result_set: result_set_id,
                    algebraic,
                },
                issues,
            );
        }
    }
}

#[allow(
    clippy::too_many_lines,
    reason = "one exhaustive match keeps every scope-boundary witness dialect and compatibility rule auditable together"
)]
fn validate_boundary_evidence(
    problem_id: SpectralProblemId,
    result_set_id: SpectralResultSetIdV1,
    requested_scope: CompletenessScopeV1,
    clusters: &[SpectralClusterV1],
    boundary: &ScopeBoundaryStateV1,
    issues: &mut Vec<SpectralTruthErrorV1>,
) {
    match boundary {
        ScopeBoundaryStateV1::NoClaim => {}
        ScopeBoundaryStateV1::Partial(state) => {
            if !matches!(requested_scope, CompletenessScopeV1::Partial { .. }) {
                issues.push(SpectralTruthErrorV1::BoundaryScopeMismatch);
            }
            match state {
                PartialBoundaryStateV1::Unknown { .. } => {}
                PartialBoundaryStateV1::Separated {
                    lower,
                    norm,
                    witness,
                } => validate_truth_witness(
                    *witness,
                    problem_id,
                    SpectralTruthPropositionV1::PartialBoundarySeparated {
                        result_set: result_set_id,
                        lower: *lower,
                        norm: *norm,
                    },
                    issues,
                ),
                PartialBoundaryStateV1::ClusterClosed { cluster, witness } => {
                    if !cluster_exists(clusters, *cluster) {
                        issues.push(SpectralTruthErrorV1::DanglingClusterReference);
                    }
                    validate_truth_witness(
                        *witness,
                        problem_id,
                        SpectralTruthPropositionV1::PartialBoundaryClusterClosed {
                            result_set: result_set_id,
                            cluster: *cluster,
                        },
                        issues,
                    );
                }
            }
        }
        ScopeBoundaryStateV1::Region(state) => {
            if !region_boundary_state_is_compatible(requested_scope, state, true) {
                issues.push(SpectralTruthErrorV1::BoundaryScopeMismatch);
            }
            match state {
                RegionBoundaryStateV1::Unknown { .. } => {}
                RegionBoundaryStateV1::Separated {
                    lower,
                    norm,
                    witness,
                } => validate_truth_witness(
                    *witness,
                    problem_id,
                    SpectralTruthPropositionV1::RegionBoundarySeparated {
                        result_set: result_set_id,
                        lower: *lower,
                        norm: *norm,
                    },
                    issues,
                ),
                RegionBoundaryStateV1::IntersectionsResolved {
                    included,
                    excluded_algebraic,
                    witness,
                } => validate_truth_witness(
                    *witness,
                    problem_id,
                    SpectralTruthPropositionV1::RegionBoundaryIntersections {
                        result_set: result_set_id,
                        included: included.clone(),
                        excluded_algebraic: *excluded_algebraic,
                    },
                    issues,
                ),
            }
        }
        ScopeBoundaryStateV1::FullSpectrum => {
            if !matches!(requested_scope, CompletenessScopeV1::FullFinite { .. }) {
                issues.push(SpectralTruthErrorV1::BoundaryScopeMismatch);
            }
        }
    }
}

/// Structured malformed-truth refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpectralTruthErrorV1 {
    /// Untrusted input exceeded the cluster-record preprocessing budget.
    TooManyClusters {
        /// Number of cluster records supplied.
        found: usize,
        /// Maximum cluster records admitted by v1.
        limit: usize,
    },
    /// Untrusted input exceeded the region-boundary reference budget.
    TooManyBoundaryReferences {
        /// Number of boundary cluster references supplied.
        found: usize,
        /// Maximum boundary references admitted by v1.
        limit: usize,
    },
    /// Interval endpoint was NaN or infinite.
    NonFiniteInterval,
    /// Interval lower endpoint exceeded upper endpoint.
    ReversedInterval,
    /// Multiplicity zero is invalid.
    ZeroMultiplicity,
    /// Multiplicity lower bound exceeded upper bound.
    ReversedMultiplicityBounds,
    /// Geometric multiplicity cannot exceed algebraic multiplicity.
    GeometricExceedsAlgebraic,
    /// A cluster's admitted geometric-multiplicity minimum exceeds the
    /// operator-space dimension.
    GeometricCapacityExceeded {
        /// Cluster whose eigenspace lower bound is impossible.
        cluster: SpectralClusterIdV1,
        /// Admitted geometric-multiplicity minimum.
        minimum: u32,
        /// Maximum eigenspace dimension of the operator space.
        dimension: u32,
    },
    /// Separation lower bound was nonpositive or nonfinite.
    InvalidSeparationBound,
    /// Residual upper bound was negative or nonfinite.
    InvalidResidualBound,
    /// Per-cluster internal state conflicts with multiplicity semantics.
    InvalidInternalClusterState,
    /// Cluster lineage appears twice.
    DuplicateCluster,
    /// Overlap set contains the same cluster twice.
    DuplicateSeparationReference,
    /// Separation state references an absent cluster.
    DanglingClusterReference,
    /// Partial coverage count is zero/inverted.
    InvalidPartialCoverage,
    /// Achieved coverage conflicts with the bound admitted request.
    CoverageScopeMismatch,
    /// Favorable prefix or region completeness lacks admitted equation and
    /// descriptor regularity establishing a discrete finite-dimensional
    /// spectrum.
    DiscreteSpectrumRegularityNotEstablished,
    /// Full finite/projective completeness lacks the admitted equation and
    /// descriptor regularity needed to make its total cardinality meaningful.
    FullCompletenessRegularityNotEstablished {
        /// Algebraic cardinality requested by the admitted scope.
        requested: u32,
        /// Cardinality established by admitted regularity theorem closure, if
        /// any.
        established: Option<u32>,
    },
    /// Cluster multiplicity accounting conflicts with achieved coverage.
    CoverageCardinalityMismatch,
    /// Scope-boundary state conflicts with the admitted request.
    BoundaryScopeMismatch,
    /// A favorable boundary state is unsupported by achieved coverage.
    BoundaryCoverageMismatch,
    /// Descriptor/projective policy conflicts with result accounting.
    InfinityPolicyMismatch,
    /// Returned projective cluster/count conflicts with infinity accounting.
    InfinityAccountingMismatch,
    /// `NoResult` retained one or more clusters.
    NoResultHasClusters,
    /// `NoResult` retained authority or scope-boundary claims.
    NoResultHasClaims,
    /// A projective-infinity cluster was returned without an admitted include
    /// policy for the descriptor problem.
    ProjectiveClusterNotAdmitted,
    /// Favorable projective-infinity membership conflicts with admitted
    /// invertibility of the pencil weight or polynomial leading coefficient.
    ProjectiveInfinityExcludedByRegularity,
    /// More than one projective-infinity cluster was returned for one result.
    MultipleProjectiveClusters,
    /// A certified enclosure under an admitted real-spectrum requirement
    /// cannot contain any real spectral value.
    RealSpectrumEnclosureConflict,
    /// Authority class and presence/absence of evidence disagree.
    AuthorityEvidenceMismatch,
    /// Admitted evidence names another exact proposition.
    WitnessPropositionMismatch {
        /// Proposition identity required by the truth statement.
        expected: SpectralPropositionId,
        /// Proposition identity carried by the supplied witness.
        found: SpectralPropositionId,
    },
    /// Typed digest matched but canonical-preimage observation differed.
    WitnessObservationMismatch {
        /// Proposition identity whose canonical observations disagreed.
        proposition: SpectralPropositionId,
    },
    /// Refused/not-started computation carried result claims.
    InactiveComputationHasClaims,
    /// Complete region/full coverage requires successful completion in v1.
    CompleteCoverageRequiresCompletion,
    /// Canonical typed identity construction failed closed.
    Identity(CanonicalError),
}

impl SpectralTruthErrorV1 {
    fn sort_key(&self) -> u16 {
        match self {
            Self::TooManyClusters { .. } => 0,
            Self::TooManyBoundaryReferences { .. } => 1,
            Self::NonFiniteInterval => 2,
            Self::ReversedInterval => 3,
            Self::ZeroMultiplicity => 4,
            Self::ReversedMultiplicityBounds => 5,
            Self::GeometricExceedsAlgebraic | Self::GeometricCapacityExceeded { .. } => 6,
            Self::InvalidSeparationBound => 7,
            Self::InvalidResidualBound => 8,
            Self::InvalidInternalClusterState => 9,
            Self::DuplicateCluster => 10,
            Self::DuplicateSeparationReference => 11,
            Self::DanglingClusterReference => 12,
            Self::InvalidPartialCoverage => 13,
            Self::CoverageScopeMismatch => 14,
            Self::DiscreteSpectrumRegularityNotEstablished => 15,
            Self::FullCompletenessRegularityNotEstablished { .. } => 16,
            Self::CoverageCardinalityMismatch => 17,
            Self::BoundaryScopeMismatch => 18,
            Self::BoundaryCoverageMismatch => 19,
            Self::InfinityPolicyMismatch => 20,
            Self::InfinityAccountingMismatch => 21,
            Self::NoResultHasClusters => 22,
            Self::NoResultHasClaims => 23,
            Self::ProjectiveClusterNotAdmitted => 24,
            Self::ProjectiveInfinityExcludedByRegularity => 25,
            Self::MultipleProjectiveClusters => 26,
            Self::RealSpectrumEnclosureConflict => 27,
            Self::AuthorityEvidenceMismatch => 28,
            Self::WitnessPropositionMismatch { .. } => 29,
            Self::WitnessObservationMismatch { .. } => 30,
            Self::InactiveComputationHasClaims => 31,
            Self::CompleteCoverageRequiresCompletion => 32,
            Self::Identity(_) => 33,
        }
    }
}

impl fmt::Display for SpectralTruthErrorV1 {
    #[allow(clippy::too_many_lines)] // One exhaustive, auditable diagnostic mapping for the public truth-error enum.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooManyClusters { found, limit } => write!(
                f,
                "spectral truth has {found} clusters; v1 admits at most {limit}"
            ),
            Self::TooManyBoundaryReferences { found, limit } => write!(
                f,
                "region boundary has {found} cluster references; v1 admits at most {limit}"
            ),
            Self::NonFiniteInterval => f.write_str("spectral interval endpoints must be finite"),
            Self::ReversedInterval => {
                f.write_str("spectral interval lower endpoint exceeds upper endpoint")
            }
            Self::ZeroMultiplicity => f.write_str("spectral multiplicity must be positive"),
            Self::ReversedMultiplicityBounds => {
                f.write_str("multiplicity lower bound exceeds upper bound")
            }
            Self::GeometricExceedsAlgebraic => {
                f.write_str("geometric multiplicity cannot exceed algebraic multiplicity")
            }
            Self::GeometricCapacityExceeded {
                cluster,
                minimum,
                dimension,
            } => write!(
                f,
                "cluster {cluster:?} has geometric multiplicity at least {minimum}, exceeding operator dimension {dimension}"
            ),
            Self::InvalidSeparationBound => {
                f.write_str("separation lower bound must be finite and strictly positive")
            }
            Self::InvalidResidualBound => {
                f.write_str("residual upper bound must be finite and nonnegative")
            }
            Self::InvalidInternalClusterState => {
                f.write_str("internal cluster state conflicts with multiplicity evidence")
            }
            Self::DuplicateCluster => f.write_str("duplicate spectral cluster lineage"),
            Self::DuplicateSeparationReference => {
                f.write_str("duplicate cluster in separation overlap set")
            }
            Self::DanglingClusterReference => {
                f.write_str("separation state references a missing cluster")
            }
            Self::InvalidPartialCoverage => f.write_str(
                "partial algebraic coverage/status conflicts with its requested boundary",
            ),
            Self::CoverageScopeMismatch => {
                f.write_str("coverage state conflicts with the admitted requested scope")
            }
            Self::DiscreteSpectrumRegularityNotEstablished => f.write_str(
                "prefix or region completeness lacks admitted discrete-spectrum regularity",
            ),
            Self::FullCompletenessRegularityNotEstablished {
                requested,
                established,
            } => write!(
                f,
                "full spectral coverage requested algebraic cardinality {requested}, but admitted regularity establishes {established:?}"
            ),
            Self::CoverageCardinalityMismatch => {
                f.write_str("coverage cardinality conflicts with exact cluster multiplicities")
            }
            Self::BoundaryScopeMismatch => {
                f.write_str("scope-boundary state conflicts with the admitted request")
            }
            Self::BoundaryCoverageMismatch => {
                f.write_str("scope-boundary claim conflicts with achieved spectral coverage")
            }
            Self::InfinityPolicyMismatch => {
                f.write_str("projective/infinity accounting conflicts with admitted policy")
            }
            Self::InfinityAccountingMismatch => {
                f.write_str("projective cluster multiplicity conflicts with infinity accounting")
            }
            Self::NoResultHasClusters => f.write_str("NoResult cannot retain clusters"),
            Self::NoResultHasClaims => {
                f.write_str("NoResult cannot retain authority or scope-boundary claims")
            }
            Self::ProjectiveClusterNotAdmitted => f.write_str(
                "projective-infinity clusters require an admitted descriptor include policy",
            ),
            Self::ProjectiveInfinityExcludedByRegularity => f.write_str(
                "projective-infinity membership conflicts with admitted weight/leading-coefficient invertibility",
            ),
            Self::MultipleProjectiveClusters => {
                f.write_str("a result may contain at most one projective-infinity cluster")
            }
            Self::RealSpectrumEnclosureConflict => f.write_str(
                "certified enclosure under admitted real-spectrum semantics excludes the real axis",
            ),
            Self::AuthorityEvidenceMismatch => {
                f.write_str("result authority and evidence presence disagree")
            }
            Self::WitnessPropositionMismatch { expected, found } => write!(
                f,
                "truth witness proposition {} does not match required {}",
                found.to_hex(),
                expected.to_hex()
            ),
            Self::WitnessObservationMismatch { proposition } => write!(
                f,
                "truth proposition {} has conflicting canonical observations",
                proposition.to_hex()
            ),
            Self::InactiveComputationHasClaims => {
                f.write_str("refused/not-started computation cannot carry result claims")
            }
            Self::CompleteCoverageRequiresCompletion => {
                f.write_str("complete spectral coverage requires Completed termination")
            }
            Self::Identity(error) => write!(f, "spectral truth identity failed: {error}"),
        }
    }
}

impl core::error::Error for SpectralTruthErrorV1 {}

/// Deterministically ranked complete truth-refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpectralTruthReportV1 {
    issues: Vec<SpectralTruthErrorV1>,
}

impl SpectralTruthReportV1 {
    fn new(mut issues: Vec<SpectralTruthErrorV1>) -> Self {
        issues.sort_by_cached_key(|issue| (issue.sort_key(), format!("{issue:?}")));
        issues.dedup();
        Self { issues }
    }

    /// Deterministically ranked, deduplicated refusal diagnostics.
    #[must_use]
    pub fn issues(&self) -> &[SpectralTruthErrorV1] {
        &self.issues
    }
}

impl fmt::Display for SpectralTruthReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "spectral truth refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for SpectralTruthReportV1 {}
