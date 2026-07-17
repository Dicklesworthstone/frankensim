//! Versioned, fail-closed spectral problem admission.
//!
//! Spectral representation, semantic role, and mathematical structure are
//! deliberately orthogonal. A descriptor polynomial may be Hamiltonian; a
//! monodromy operator may be symplectic; and normality is meaningful only in
//! a named metric. This module therefore rejects the tempting but incorrect
//! design of one mutually-exclusive "problem class" enum.
//!
//! The types here classify and bind evidence. They do not prove the supplied
//! witnesses, select an implementation, or make a spectral-completeness
//! claim. Iterative algorithms consume only [`ValidatedSpectralProblemV1`],
//! while [`assess_method_class`] is a pure obligation check for the later
//! routing layer.

use core::fmt;

use fs_blake3::hash_domain;
use fs_blake3::identity::{
    AuthorityRef, CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field,
    FieldSpec, IdentityAuditRecord, IdentityReceipt, KeyPolicyId, NeverCancel, ObservedIdentity,
    PolicyRelativeAdmitted, Presented, ProblemSemanticId, PromotionAuditRecord, PromotionRefusal,
    PromotionTrustRoot, PromotionWitness, StrongIdentity, Verified, VerifierId, WireType,
};
use fs_qty::{Angle, Dims, QtyAny, Time};

/// Current version of the spectral problem admission and identity schema.
///
/// Version 2 binds the root-relative promotion audit carried by every
/// favorable witness. Version-1 descriptors are refused rather than silently
/// reinterpreted under the stronger identity encoding.
pub const SPECTRAL_PROBLEM_SCHEMA_VERSION: u32 = 2;
/// Maximum structure claims admitted by the v1 Rust descriptor layout before
/// any sorting or pair checks. Identity schema v2 leaves this envelope intact.
pub const MAX_STRUCTURE_CLAIMS_V1: usize = 256;
/// Maximum mutually compatible regularity families in the v1 Rust descriptor
/// layout. This is also the pre-sort input cap; adding a new family requires
/// revisiting the canonical envelope and this bound together.
pub const MAX_REGULARITY_CLAIMS_V1: usize = 5;

const IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(1 << 18, 1 << 16, 16, 4096, 4096);
// A schema-v2 problem can carry 256 promotion-bearing structure claims. Their
// one canonical-set field is larger than the tighter verifier/policy field
// envelope. Keep its field bound at the complete 256 KiB frame cap and isolate
// that broader allowance to problem identities.
const PROBLEM_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 18, 1 << 18, 16, 4096, 4096);
// Truth propositions may carry the full bounded region-boundary reference set
// (4096 typed 32-byte IDs plus framing). Keep this distinct from the tighter
// verifier/policy descriptor limit so the public truth cap is actually
// encodable without broadening unrelated identity inputs.
const PROPOSITION_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 18, 1 << 18, 16, 4096, 4096);
const SPECTRAL_PROMOTION_WITNESS_ENCODING_DOMAIN_V1: &[u8] =
    b"org.frankensim.fs-spectral.promotion-witness-encoding.v1";

macro_rules! opaque_digest_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct from exact typed digest bytes.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact digest bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }
    };
}

opaque_digest_id!(
    /// Identity of the operator, pencil, polynomial, map, or operator function
    /// whose spectrum is requested. Identity is not an authority claim.
    SpectralSubjectId
);
opaque_digest_id!(
    /// Identity of a metric artifact. Domain and codomain metrics remain
    /// distinct by their role in [`SpectralSpaceContextV1`]; equal bytes mean
    /// the same metric artifact and conflicting descriptors are refused.
    SpectralMetricId
);
opaque_digest_id!(
    /// Identity of a non-metric form such as a symplectic two-form,
    /// Krein/J form, or conjugation/real-structure artifact.
    SpectralFormId
);
opaque_digest_id!(
    /// Identity of the norm/model used to interpret a structure-defect
    /// tolerance. Equal numeric tolerances in different norms are not
    /// interchangeable.
    SpectralNormId
);
opaque_digest_id!(
    /// Identity of the complete left/right/operator/inverse scaling bundle.
    SpectralScalingId
);
opaque_digest_id!(
    /// Identity of one scaling map within a scaling bundle.
    SpectralScalingMapId
);
opaque_digest_id!(
    /// Identity of an analytic operator function or branch definition.
    SpectralFunctionId
);
opaque_digest_id!(
    /// Identity of a continuation path/mode-lineage artifact used to track a
    /// continuous Floquet logarithm branch through crossings.
    SpectralContinuationId
);
opaque_digest_id!(
    /// Identity of a requested spectral region.
    SpectralRegionId
);
opaque_digest_id!(
    /// Identity of the projective chart used to order finite coordinates and
    /// distinguish the projective point at infinity.
    SpectralProjectiveChartId
);
opaque_digest_id!(
    /// Identity of exact gauge-fixing conditions and their lineage artifact.
    SpectralGaugeArtifactId
);
opaque_digest_id!(
    /// Identity of the quotient/reduction map whose target is the declared
    /// post-reduction operator space.
    SpectralQuotientMapId
);

/// Static identity schema for a fully validated spectral problem descriptor.
pub enum SpectralProblemIdentitySchemaV2 {}

impl CanonicalSchema for SpectralProblemIdentitySchemaV2 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.problem-semantic.v2";
    const NAME: &'static str = "spectral-problem-semantic";
    const VERSION: u32 = SPECTRAL_PROBLEM_SCHEMA_VERSION;
    const CONTEXT: &'static str = "operator, form, role, structure witnesses, units, scaling, metrics, regularity, ordering, and requested scope";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("subject", WireType::Bytes),
        FieldSpec::required("scalar-field", WireType::Variant),
        FieldSpec::required("class", WireType::Bytes),
        FieldSpec::required("structure-claims", WireType::CanonicalSet),
        FieldSpec::required("scaling", WireType::Bytes),
        FieldSpec::required("spaces", WireType::Bytes),
        FieldSpec::required("regularity", WireType::CanonicalSet),
        FieldSpec::required("ordering", WireType::Variant),
        FieldSpec::required("requested-scope", WireType::Variant),
    ];
}

/// Domain-separated identity schema for one exact spectral proposition.
///
/// A proposition identity is narrower than a problem identity: it binds the
/// exact statement an evidence artifact was asked to support, including its
/// form/metric, polarity, bounds, and semantic context. It is not itself
/// evidence or authority.
pub enum SpectralPropositionIdentitySchemaV1 {}

impl CanonicalSchema for SpectralPropositionIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.proposition-semantic.v1";
    const NAME: &'static str = "spectral-proposition-semantic";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str =
        "exact spectral proposition, bound problem/form context, polarity, and parameters";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("kind", WireType::Variant),
        FieldSpec::required("payload", WireType::Bytes),
    ];
}

/// Typed identity of an exact proposition. Equality is not proof that the
/// proposition is true.
pub type SpectralPropositionId = ProblemSemanticId<SpectralPropositionIdentitySchemaV1>;

/// Identity namespace for the verifier implementation/configuration that
/// checked a spectral proposition.
pub enum SpectralVerifierIdentitySchemaV1 {}

impl CanonicalSchema for SpectralVerifierIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.verifier.v1";
    const NAME: &'static str = "spectral-verifier";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "spectral proposition verifier implementation and configuration";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("descriptor", WireType::Bytes)];
}

/// Identity namespace for the policy that admits a verifier decision for use
/// by a method-family or truth constructor.
pub enum SpectralAuthorityPolicySchemaV1 {}

impl CanonicalSchema for SpectralAuthorityPolicySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.authority-policy.v1";
    const NAME: &'static str = "spectral-authority-policy";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "policy admitting verified spectral propositions";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("descriptor", WireType::Bytes)];
}

/// The only authority typestate from which a favorable spectral witness token
/// can be minted. Presented or merely verified material cannot cross this
/// boundary.
pub type SpectralAuthorityVerifierIdV1 = VerifierId<SpectralVerifierIdentitySchemaV1>;
/// Typed identity of the policy admitting a spectral verifier decision.
pub type SpectralAuthorityPolicyIdV1 = KeyPolicyId<SpectralAuthorityPolicySchemaV1>;
/// Configured trust root required to promote policy-relative spectral
/// admission into a favorable witness.
pub type SpectralPromotionTrustRootV1 =
    PromotionTrustRoot<SpectralVerifierIdentitySchemaV1, SpectralAuthorityPolicySchemaV1>;
/// Opaque promotion decision minted by a [`SpectralPromotionTrustRootV1`].
pub type SpectralPromotionWitnessV1 = PromotionWitness<
    SpectralPropositionId,
    SpectralVerifierIdentitySchemaV1,
    SpectralAuthorityPolicySchemaV1,
>;
/// Stable context bound into every spectral promotion decision created via
/// [`spectral_promotion_trust_root`].
pub const SPECTRAL_PROMOTION_CONTEXT_V1: &str = "org.frankensim.fs-spectral.promotion.v1";
/// Untrusted presented spectral authority.
pub type PresentedSpectralAuthorityV1 = AuthorityRef<
    SpectralPropositionId,
    SpectralVerifierIdentitySchemaV1,
    SpectralAuthorityPolicySchemaV1,
    Presented,
>;
/// Verifier-accepted but not yet policy-admitted spectral authority.
pub type VerifiedSpectralAuthorityV1 = AuthorityRef<
    SpectralPropositionId,
    SpectralVerifierIdentitySchemaV1,
    SpectralAuthorityPolicySchemaV1,
    Verified,
>;
/// Verifier-accepted and separately policy-admitted spectral authority.
///
/// This state is explicitly policy-relative. It is not sufficient to mint a
/// favorable spectral witness and must first pass a configured
/// [`SpectralPromotionTrustRootV1`].
pub type PolicyRelativeSpectralAuthorityV1 = AuthorityRef<
    SpectralPropositionId,
    SpectralVerifierIdentitySchemaV1,
    SpectralAuthorityPolicySchemaV1,
    PolicyRelativeAdmitted,
>;
/// Compatibility name for policy-relative spectral admission. Despite the
/// historical name, this type is not promotion authority.
pub type AdmittedSpectralAuthorityV1 = PolicyRelativeSpectralAuthorityV1;

/// Configure a spectral promotion root from independently retained canonical
/// verifier and policy receipts.
///
/// The caller remains responsible for choosing meaningful verifier and policy
/// configurations. The root binds both typed identities and their exact
/// canonical-byte observations; digest-only configuration is not accepted.
/// This v1 root is configuration-relative and does not authenticate an owner
/// or root instance.
///
/// # Errors
///
/// Returns [`PromotionRefusal`] if the fixed spectral promotion context is
/// invalid.
pub const fn spectral_promotion_trust_root(
    verifier: IdentityReceipt<SpectralAuthorityVerifierIdV1>,
    policy: IdentityReceipt<SpectralAuthorityPolicyIdV1>,
) -> Result<SpectralPromotionTrustRootV1, PromotionRefusal> {
    SpectralPromotionTrustRootV1::configure(
        ObservedIdentity::from_receipt(verifier),
        ObservedIdentity::from_receipt(policy),
        SPECTRAL_PROMOTION_CONTEXT_V1,
    )
}

/// Exact binding mismatch between a policy-relative authority and an opaque
/// spectral promotion decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralPromotionBindingErrorV1 {
    /// The subject receipt, preimage, or bounded receipt metadata differs.
    Subject,
    /// The external anchor differs.
    Anchor,
    /// The typed verifier identity differs.
    Verifier,
    /// The typed key-policy identity differs.
    KeyPolicy,
    /// The promotion witness was minted for another configured context.
    Context,
}

impl fmt::Display for SpectralPromotionBindingErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Subject => f.write_str("spectral promotion subject receipt differs"),
            Self::Anchor => f.write_str("spectral promotion anchor differs"),
            Self::Verifier => f.write_str("spectral promotion verifier differs"),
            Self::KeyPolicy => f.write_str("spectral promotion key policy differs"),
            Self::Context => f.write_str("spectral promotion context differs"),
        }
    }
}

impl core::error::Error for SpectralPromotionBindingErrorV1 {}

/// Copyable favorable token produced only from a configured spectral
/// promotion witness.
///
/// Generic verifier/admitter capabilities yield only
/// [`PolicyRelativeSpectralAuthorityV1`]. They cannot be converted directly
/// into this type; a separately configured [`SpectralPromotionTrustRootV1`]
/// must first mint an opaque [`SpectralPromotionWitnessV1`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmittedSpectralWitnessV1 {
    proposition: SpectralPropositionId,
    audit: IdentityAuditRecord,
    promotion: PromotionAuditRecord,
}

impl AdmittedSpectralWitnessV1 {
    /// Inspect a policy-relative authority together with the opaque,
    /// root-relative promotion decision for that exact binding.
    ///
    /// Policy-relative admission alone is deliberately insufficient:
    ///
    /// ```compile_fail,E0061
    /// use fs_spectral::admission::{
    ///     AdmittedSpectralAuthorityV1, AdmittedSpectralWitnessV1,
    /// };
    /// fn foreign_policy_cannot_promote(admitted: AdmittedSpectralAuthorityV1) {
    ///     let _ = AdmittedSpectralWitnessV1::from_authority(&admitted);
    /// }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`SpectralPromotionBindingErrorV1`] unless the authority and
    /// promotion decision bind the exact same subject receipt, anchor,
    /// verifier, key policy, and fixed spectral promotion context.
    pub fn from_authority(
        authority: &AdmittedSpectralAuthorityV1,
        promotion: SpectralPromotionWitnessV1,
    ) -> Result<Self, SpectralPromotionBindingErrorV1> {
        if authority.receipt() != promotion.subject() {
            return Err(SpectralPromotionBindingErrorV1::Subject);
        }
        if authority.anchor() != promotion.anchor() {
            return Err(SpectralPromotionBindingErrorV1::Anchor);
        }
        if authority.verifier() != promotion.verifier().id() {
            return Err(SpectralPromotionBindingErrorV1::Verifier);
        }
        if authority.key_policy() != promotion.key_policy().id() {
            return Err(SpectralPromotionBindingErrorV1::KeyPolicy);
        }
        if promotion.context() != SPECTRAL_PROMOTION_CONTEXT_V1 {
            return Err(SpectralPromotionBindingErrorV1::Context);
        }
        Ok(Self {
            proposition: authority.receipt().id(),
            audit: authority.audit_record(),
            promotion: promotion.audit(),
        })
    }

    /// Exact proposition identity checked by the admitted verifier/policy.
    #[must_use]
    pub const fn proposition(&self) -> SpectralPropositionId {
        self.proposition
    }

    /// Bounded subject, anchor, verifier/policy identity, trust-state, and
    /// no-claim audit record retained from policy-relative admission.
    #[must_use]
    pub const fn audit(&self) -> IdentityAuditRecord {
        self.audit
    }

    /// Bounded verifier/policy namespace, canonical-byte observations, and
    /// configured promotion context.
    #[must_use]
    pub const fn promotion_audit(&self) -> PromotionAuditRecord {
        self.promotion
    }

    /// Whether this token is bound to the typed digest, independently retained
    /// canonical-preimage root, and exact canonical byte length of the expected
    /// proposition receipt.
    #[must_use]
    pub fn matches_receipt(&self, expected: IdentityReceipt<SpectralPropositionId>) -> bool {
        self.proposition == expected.id()
            && self.audit.canonical_preimage() == expected.canonical_preimage()
            && self.audit.canonical_bytes() == expected.canonical_bytes()
    }
}

/// Typed semantic identity of one validated problem descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpectralProblemId(ProblemSemanticId<SpectralProblemIdentitySchemaV2>);

impl SpectralProblemId {
    /// Exact typed digest bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        self.0.as_bytes()
    }

    /// Lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_hex(self) -> String {
        self.0.to_hex()
    }
}

/// Complete producer receipt for one canonical spectral-problem identity.
///
/// The public [`SpectralProblemId`] remains a narrow semantic identifier. This
/// receipt additionally retains the independently adjudicable canonical-frame
/// root, exact byte length, field count, collection cardinality, and limits
/// that produced those digest bytes.
pub type SpectralProblemIdentityReceiptV2 =
    IdentityReceipt<ProblemSemanticId<SpectralProblemIdentitySchemaV2>>;

/// Scalar field of the admitted operator spaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SpectralScalarFieldV1 {
    /// Real scalar field.
    Real,
    /// Complex scalar field.
    Complex,
}

impl SpectralScalarFieldV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Real => 0,
            Self::Complex => 1,
        }
    }
}

/// Algebraic representation of the spectral equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SpectralRepresentationV1 {
    /// Standard linear problem `A x = lambda x`.
    StandardLinear,
    /// Generalized pencil `A x = lambda B x`.
    GeneralizedPencil,
    /// Matrix polynomial `sum_k P_k lambda^k`, with declared grade.
    MatrixPolynomial {
        /// Declared polynomial grade, including a possibly singular leading
        /// coefficient whose treatment is governed by regularity evidence.
        grade: u32,
    },
}

impl SpectralRepresentationV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::StandardLinear => 0,
            Self::GeneralizedPencil => 1,
            Self::MatrixPolynomial { .. } => 2,
        }
    }
}

/// Whether generalized/polynomial semantics include a descriptor system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DescriptorRoleV1 {
    /// Ordinary, non-descriptor problem.
    Ordinary,
    /// Descriptor semantics, including explicit infinite-eigenvalue policy.
    Descriptor {
        /// How infinite eigenvalues are represented in the requested result.
        infinity_policy: InfiniteEigenvaluePolicyV1,
    },
}

impl DescriptorRoleV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Ordinary => 0,
            Self::Descriptor { .. } => 1,
        }
    }
}

/// Explicit handling of projective/infinite eigenvalues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InfiniteEigenvaluePolicyV1 {
    /// Infinite eigenvalues are included as projective result clusters.
    IncludeProjective,
    /// They are excluded, with count and exclusion recorded by the solver.
    ExcludeWithCount,
    /// The problem declares no claim about infinite eigenvalues.
    NoClaim,
}

impl InfiniteEigenvaluePolicyV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::IncludeProjective => 0,
            Self::ExcludeWithCount => 1,
            Self::NoClaim => 2,
        }
    }
}

/// Whether a Floquet request returns multipliers or logarithmic exponents.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FloquetParameterV1 {
    /// Dimensionless monodromy multipliers.
    Multiplier,
    /// Floquet exponents with inverse-time dimensions.
    Exponent,
}

impl FloquetParameterV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Multiplier => 0,
            Self::Exponent => 1,
        }
    }
}

/// Explicit logarithm/phase convention for Floquet results.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FloquetBranchConventionV1 {
    /// No logarithm is taken; only multipliers are requested.
    MultipliersOnly,
    /// Principal complex logarithm.
    PrincipalLog,
    /// A continuous branch tracked by a retained continuation/mode-lineage
    /// artifact and anchored at a finite phase in radians.
    ContinuousFrom {
        /// Exact continuation/lineage artifact.
        continuation: SpectralContinuationId,
        /// Dimensionless phase anchor in radians.
        anchor_phase: Angle,
    },
}

impl FloquetBranchConventionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::MultipliersOnly => 0,
            Self::PrincipalLog => 1,
            Self::ContinuousFrom { .. } => 2,
        }
    }
}

/// Branch policy for an analytic operator function.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperatorFunctionBranchPolicyV1 {
    /// The admitted domain is single-valued and has no branch cut.
    SingleValued,
    /// A branch is bound by the supplied function identity.
    ExplicitBranch,
    /// No branch-compatibility claim is made.
    NoClaim,
}

impl OperatorFunctionBranchPolicyV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::SingleValued => 0,
            Self::ExplicitBranch => 1,
            Self::NoClaim => 2,
        }
    }
}

/// Provenance/semantic origin of the operator whose spectrum is requested.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectralOperatorOriginV1 {
    /// Direct algebraic operator/pencil/polynomial.
    Direct,
    /// Monodromy map of a periodic evolution.
    MonodromyFloquet {
        /// Physical period in coherent SI seconds, enforced by the type.
        period: Time,
        /// Multiplier versus exponent semantics.
        parameter: FloquetParameterV1,
        /// Explicit logarithm/phase convention.
        branch: FloquetBranchConventionV1,
    },
    /// Analytic operator-function problem.
    AnalyticOperatorFunction {
        /// Identity of the exact function/domain artifact.
        function: SpectralFunctionId,
        /// Branch policy admitted by that artifact.
        branch_policy: OperatorFunctionBranchPolicyV1,
    },
}

impl SpectralOperatorOriginV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Direct => 0,
            Self::MonodromyFloquet { .. } => 1,
            Self::AnalyticOperatorFunction { .. } => 2,
        }
    }
}

/// Product classification of equation form, descriptor semantics, and origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralProblemClassV1 {
    representation: SpectralRepresentationV1,
    descriptor: DescriptorRoleV1,
    origin: SpectralOperatorOriginV1,
}

impl SpectralProblemClassV1 {
    /// Construct an orthogonal class product. Validation of cross-field
    /// obligations happens in [`validate_problem`].
    #[must_use]
    pub const fn new(
        representation: SpectralRepresentationV1,
        descriptor: DescriptorRoleV1,
        origin: SpectralOperatorOriginV1,
    ) -> Self {
        Self {
            representation,
            descriptor,
            origin,
        }
    }

    /// Algebraic representation.
    #[must_use]
    pub const fn representation(&self) -> SpectralRepresentationV1 {
        self.representation
    }

    /// Descriptor role.
    #[must_use]
    pub const fn descriptor(&self) -> DescriptorRoleV1 {
        self.descriptor
    }

    /// Operator origin.
    #[must_use]
    pub const fn origin(&self) -> SpectralOperatorOriginV1 {
        self.origin
    }
}

/// Mathematical property attached to a named metric and an evidence witness.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StructurePropertyV1 {
    /// Adjoint equality in the witness-bound metric.
    SelfAdjoint,
    /// Normality in the witness-bound metric.
    Normal,
    /// Nonnormality in the witness-bound metric.
    Nonnormal,
    /// Hamiltonian generator identity.
    Hamiltonian,
    /// Symplectic map identity.
    Symplectic,
    /// J-self-adjoint identity.
    JSelfAdjoint,
    /// Composite Hermitian-definite pencil proposition: the pencil weight is
    /// Hermitian positive definite and the induced operator is self-adjoint in
    /// that weight. This is deliberately stronger than relabeling an
    /// undifferentiated generalized pencil `SelfAdjoint`.
    HermitianDefinitePencil,
    /// Gyroscopic polynomial structure.
    Gyroscopic,
    /// Palindromic polynomial structure under the named involution/parity.
    Palindromic {
        /// Sign under coefficient reversal.
        parity: PalindromicParityV1,
        /// Transpose versus conjugate-transpose coefficient involution.
        involution: PolynomialInvolutionV1,
    },
    /// Entire requested spectral set is real. This is not implied by a real
    /// coefficient field or by conjugate-pair symmetry.
    RealSpectrum,
    /// Real-coefficient conjugate-pair symmetry.
    RealConjugatePairs,
}

impl StructurePropertyV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::SelfAdjoint => 0,
            Self::Normal => 1,
            Self::Nonnormal => 2,
            Self::Hamiltonian => 3,
            Self::Symplectic => 4,
            Self::JSelfAdjoint => 5,
            Self::HermitianDefinitePencil => 6,
            Self::Gyroscopic => 7,
            Self::Palindromic { .. } => 8,
            Self::RealSpectrum => 9,
            Self::RealConjugatePairs => 10,
        }
    }
}

/// Reversal sign for a structured matrix polynomial.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PalindromicParityV1 {
    /// `P(lambda)` equals its reversed involution.
    Palindromic,
    /// `P(lambda)` is the negative of its reversed involution.
    AntiPalindromic,
}

impl PalindromicParityV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Palindromic => 0,
            Self::AntiPalindromic => 1,
        }
    }
}

/// Coefficient involution used by a palindromic polynomial proposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PolynomialInvolutionV1 {
    /// Plain transpose.
    Transpose,
    /// Complex conjugate transpose.
    ConjugateTranspose,
}

impl PolynomialInvolutionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Transpose => 0,
            Self::ConjugateTranspose => 1,
        }
    }
}

/// Exact geometric object with respect to which a structure proposition is
/// stated. Equal digest bytes in different variants never cross-satisfy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StructureSupportV1 {
    /// Named nonsingular metric/inner-product model; this may be positive,
    /// indefinite, or unresolved, but never an admitted singular metric.
    InnerProduct(SpectralMetricId),
    /// Nondegenerate skew symplectic form.
    SymplecticForm(SpectralFormId),
    /// Nondegenerate indefinite Hermitian/Krein form.
    KreinForm(SpectralFormId),
    /// Conjugation/real-structure artifact.
    Conjugation(SpectralFormId),
    /// Proposition is independent of an auxiliary form.
    FormFree,
}

impl StructureSupportV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::InnerProduct(_) => 0,
            Self::SymplecticForm(_) => 1,
            Self::KreinForm(_) => 2,
            Self::Conjugation(_) => 3,
            Self::FormFree => 4,
        }
    }
}

/// Support requirement reported by a method-family obligation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructureSupportRequirementV1 {
    /// One exact metric/form artifact.
    Exact(StructureSupportV1),
    /// Any proposition-bound nondegenerate symplectic form.
    SymplecticForm,
    /// Any proposition-bound nondegenerate Krein/J form.
    KreinForm,
}

impl StructureSupportRequirementV1 {
    fn accepts(self, support: StructureSupportV1) -> bool {
        match self {
            Self::Exact(expected) => support == expected,
            Self::SymplecticForm => matches!(support, StructureSupportV1::SymplecticForm(_)),
            Self::KreinForm => matches!(support, StructureSupportV1::KreinForm(_)),
        }
    }
}

/// Whether a witness supports or contradicts its named proposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum WitnessDispositionV1 {
    /// Evidence supports the proposition within its declared tolerance/scope.
    Witnessed,
    /// Evidence contradicts the proposition within its declared scope.
    Contradicted,
}

impl WitnessDispositionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Witnessed => 0,
            Self::Contradicted => 1,
        }
    }
}

/// One witnessed or contradicted structure proposition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StructureClaimV1 {
    property: StructurePropertyV1,
    support: StructureSupportV1,
    disposition: WitnessDispositionV1,
    tolerance: f64,
    norm: SpectralNormId,
    witness: AdmittedSpectralWitnessV1,
}

impl StructureClaimV1 {
    /// Construct one explicit property claim.
    #[must_use]
    pub const fn new(
        property: StructurePropertyV1,
        support: StructureSupportV1,
        disposition: WitnessDispositionV1,
        tolerance: f64,
        norm: SpectralNormId,
        witness: AdmittedSpectralWitnessV1,
    ) -> Self {
        Self {
            property,
            support,
            disposition,
            tolerance,
            norm,
            witness,
        }
    }

    /// Property being assessed.
    #[must_use]
    pub const fn property(&self) -> StructurePropertyV1 {
        self.property
    }

    /// Metric/form/non-form support of this exact proposition.
    #[must_use]
    pub const fn support(&self) -> StructureSupportV1 {
        self.support
    }

    /// Supporting versus contradicting disposition.
    #[must_use]
    pub const fn disposition(&self) -> WitnessDispositionV1 {
        self.disposition
    }

    /// Finite nonnegative defect tolerance in the named norm.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Norm/model identity in which the tolerance is interpreted.
    #[must_use]
    pub const fn norm(&self) -> SpectralNormId {
        self.norm
    }

    /// Admitted evidence for the exact proposition.
    #[must_use]
    pub const fn witness(&self) -> &AdmittedSpectralWitnessV1 {
        &self.witness
    }
}

/// Deterministically ordered structure claims. Properties may repeat across
/// distinct supports, norms, or tolerances, but one normalized proposition
/// cannot be duplicated or both witnessed and contradicted.
#[derive(Debug, Clone, PartialEq)]
pub struct StructureProfileV1 {
    claims: Vec<StructureClaimV1>,
}

impl StructureProfileV1 {
    /// Construct a raw profile. Canonical ordering and contradiction checks
    /// occur in [`validate_problem`].
    #[must_use]
    pub fn new(claims: Vec<StructureClaimV1>) -> Self {
        Self { claims }
    }

    /// Claims in caller order. A validated problem exposes canonical order.
    #[must_use]
    pub fn claims(&self) -> &[StructureClaimV1] {
        &self.claims
    }
}

/// Definiteness/nondegeneracy state of a named metric.
#[allow(clippy::large_enum_variant)] // Admitted audit evidence stays inline and replay-complete.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MetricDefinitenessV1 {
    /// Canonical Euclidean metric.
    Euclidean,
    /// Positive-definite metric with finite spectral bounds and evidence.
    PositiveDefinite {
        /// Strictly positive lower eigenvalue bound.
        lower: f64,
        /// Finite upper eigenvalue bound, `>= lower`.
        upper: f64,
        /// Evidence binding the bounds to the metric artifact.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Nondegenerate indefinite/Krein metric with declared signature.
    Indefinite {
        /// Positive signature count.
        positive: u32,
        /// Negative signature count.
        negative: u32,
        /// Evidence binding signature and nondegeneracy.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Singular/semidefinite metric. It is representable but cannot authorize
    /// ordinary self-adjoint Lanczos.
    Singular {
        /// Proven rank.
        rank: u32,
        /// Evidence binding the rank claim.
        witness: AdmittedSpectralWitnessV1,
    },
    /// No definiteness claim.
    Unknown,
}

/// Witness-free metric proposition used to mint a proposition receipt before
/// the admitted witness is embedded into [`MetricDefinitenessV1`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MetricDefinitenessPropositionV1 {
    /// Strictly positive finite spectral bounds.
    PositiveDefinite {
        /// Strictly positive lower eigenvalue bound.
        lower: f64,
        /// Finite upper eigenvalue bound, no smaller than `lower`.
        upper: f64,
    },
    /// Exact nondegenerate indefinite signature.
    Indefinite {
        /// Positive signature count.
        positive: u32,
        /// Negative signature count.
        negative: u32,
    },
    /// Exact singular rank.
    Singular {
        /// Proven rank of the metric artifact.
        rank: u32,
    },
}

impl MetricDefinitenessPropositionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::PositiveDefinite { .. } => 0,
            Self::Indefinite { .. } => 1,
            Self::Singular { .. } => 2,
        }
    }
}

impl MetricDefinitenessV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Euclidean => 0,
            Self::PositiveDefinite { .. } => 1,
            Self::Indefinite { .. } => 2,
            Self::Singular { .. } => 3,
            Self::Unknown => 4,
        }
    }

    const fn is_positive_definite(self) -> bool {
        matches!(self, Self::Euclidean | Self::PositiveDefinite { .. })
    }

    const fn is_adjoint_compatible(self) -> bool {
        matches!(
            self,
            Self::Euclidean | Self::PositiveDefinite { .. } | Self::Indefinite { .. }
        )
    }
}

/// One domain or codomain metric.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralMetricV1 {
    id: SpectralMetricId,
    dimension: u32,
    definiteness: MetricDefinitenessV1,
}

impl SpectralMetricV1 {
    /// Construct the unique built-in Euclidean metric for a dimension. The
    /// identity is domain-separated from all caller-supplied metric artifacts,
    /// so a fabricated digest cannot acquire witness-free definiteness.
    #[must_use]
    pub fn euclidean(dimension: u32) -> Self {
        let id = SpectralMetricId::from_bytes(
            hash_domain(
                "org.frankensim.fs-spectral.euclidean-metric.v1",
                &dimension.to_le_bytes(),
            )
            .0,
        );
        Self {
            id,
            dimension,
            definiteness: MetricDefinitenessV1::Euclidean,
        }
    }

    /// Construct a metric descriptor.
    #[must_use]
    pub const fn new(
        id: SpectralMetricId,
        dimension: u32,
        definiteness: MetricDefinitenessV1,
    ) -> Self {
        Self {
            id,
            dimension,
            definiteness,
        }
    }

    /// Metric identity.
    #[must_use]
    pub const fn id(&self) -> SpectralMetricId {
        self.id
    }

    /// Space dimension.
    #[must_use]
    pub const fn dimension(&self) -> u32 {
        self.dimension
    }

    /// Definiteness state.
    #[must_use]
    pub const fn definiteness(&self) -> MetricDefinitenessV1 {
        self.definiteness
    }
}

/// Gauge/nullspace convention attached to the operator spaces.
#[allow(clippy::large_enum_variant)] // Admitted audit evidence stays inline and replay-complete.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaugeConventionV1 {
    /// A checker admitted that no structural gauge/nullspace is present.
    CertifiedNone {
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Gauge fixed by explicit conditions.
    Fixed {
        /// Certified nullity before fixing.
        nullity: u32,
        /// Exact gauge-condition/constraint artifact.
        gauge: SpectralGaugeArtifactId,
        /// Evidence for the gauge/nullity statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Declared domain/codomain metrics describe an already-induced quotient
    /// operator space.
    Quotiented {
        /// Certified nullity in the pre-quotient space.
        nullity: u32,
        /// Exact reduction map targeting the declared quotient space.
        quotient: SpectralQuotientMapId,
        /// Evidence for the nullspace/quotient statement.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Gauge/nullity convention is unresolved.
    Unknown,
}

/// Witness-free gauge proposition used when producing admitted evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaugePropositionV1 {
    /// No structural gauge/nullspace is present.
    None,
    /// Gauge is fixed and had the stated nullity before fixing.
    Fixed {
        /// Certified nullity before applying the gauge-fixing conditions.
        nullity: u32,
        /// Exact gauge-condition/constraint artifact.
        gauge: SpectralGaugeArtifactId,
    },
    /// Declared operator spaces are induced by an exact quotient map.
    Quotiented {
        /// Certified dimension of the pre-quotient nullspace.
        nullity: u32,
        /// Exact quotient/reduction map targeting the declared spaces.
        quotient: SpectralQuotientMapId,
    },
}

/// Witness-free identity context binding serialization/nullspace evidence to
/// the exact gauge or reduction lineage active for the declared spaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GaugeContextV1 {
    /// No structural gauge/nullspace is present.
    CertifiedNone,
    /// Exact gauge-fixing lineage.
    Fixed {
        /// Certified nullity before gauge fixing.
        nullity: u32,
        /// Exact gauge-condition/constraint artifact.
        gauge: SpectralGaugeArtifactId,
    },
    /// Exact quotient/reduction lineage.
    Quotiented {
        /// Certified pre-reduction nullity.
        nullity: u32,
        /// Exact quotient/reduction map.
        quotient: SpectralQuotientMapId,
    },
    /// Gauge/nullspace lineage is unresolved. Evidence bound here cannot be
    /// replayed after a later gauge or quotient choice.
    Unknown,
}

impl GaugeContextV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::CertifiedNone => 0,
            Self::Fixed { .. } => 1,
            Self::Quotiented { .. } => 2,
            Self::Unknown => 3,
        }
    }
}

impl GaugePropositionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::None => 0,
            Self::Fixed { .. } => 1,
            Self::Quotiented { .. } => 2,
        }
    }
}

impl GaugeConventionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::CertifiedNone { .. } => 0,
            Self::Fixed { .. } => 1,
            Self::Quotiented { .. } => 2,
            Self::Unknown => 3,
        }
    }

    /// Witness-free semantic context for binding related serialization
    /// evidence to this exact gauge/reduction lineage.
    #[must_use]
    pub const fn context(self) -> GaugeContextV1 {
        match self {
            Self::CertifiedNone { .. } => GaugeContextV1::CertifiedNone,
            Self::Fixed { nullity, gauge, .. } => GaugeContextV1::Fixed { nullity, gauge },
            Self::Quotiented {
                nullity, quotient, ..
            } => GaugeContextV1::Quotiented { nullity, quotient },
            Self::Unknown => GaugeContextV1::Unknown,
        }
    }
}

/// Whether nullspace zeros are present in serialized spectral samples.
#[allow(clippy::large_enum_variant)] // Admitted audit evidence stays inline and replay-complete.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZeroPaddingConventionV1 {
    /// A checker admitted that no structural null zeros are present.
    CertifiedNonePresent {
        /// Exact serialization/nullity proposition evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Exact number of zeros explicitly included.
    ExplicitlyPadded {
        /// Included structural-zero count. For an already-quotiented problem
        /// this describes the pre-reduction serialization lineage and may
        /// exceed the declared target-space dimension.
        count: u32,
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Exact number of structural zeros intentionally omitted.
    Omitted {
        /// Omitted structural-zero count. For an already-quotiented problem
        /// this describes the pre-reduction serialization lineage and may
        /// exceed the declared target-space dimension.
        count: u32,
        /// Exact proposition-bound evidence.
        witness: AdmittedSpectralWitnessV1,
    },
    /// Convention is unresolved; gap methods must fail closed.
    Unknown,
}

/// Witness-free serialization/nullspace proposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZeroPaddingPropositionV1 {
    /// No structural null zeros are present.
    NonePresent,
    /// Exact number of zeros serialized into the samples.
    ExplicitlyPadded {
        /// Exact number of structural zeros included in serialized samples.
        count: u32,
    },
    /// Exact number of structural zeros intentionally omitted.
    Omitted {
        /// Exact number of structural zeros omitted from serialized samples.
        count: u32,
    },
}

impl ZeroPaddingPropositionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::NonePresent => 0,
            Self::ExplicitlyPadded { .. } => 1,
            Self::Omitted { .. } => 2,
        }
    }
}

impl ZeroPaddingConventionV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::CertifiedNonePresent { .. } => 0,
            Self::ExplicitlyPadded { .. } => 1,
            Self::Omitted { .. } => 2,
            Self::Unknown => 3,
        }
    }
}

/// Domain/codomain spaces, metrics, and nullspace conventions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralSpaceContextV1 {
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
    gauge: GaugeConventionV1,
    zero_padding: ZeroPaddingConventionV1,
}

impl SpectralSpaceContextV1 {
    /// Construct an explicit space context.
    #[must_use]
    pub const fn new(
        domain: SpectralMetricV1,
        codomain: SpectralMetricV1,
        gauge: GaugeConventionV1,
        zero_padding: ZeroPaddingConventionV1,
    ) -> Self {
        Self {
            domain,
            codomain,
            gauge,
            zero_padding,
        }
    }

    /// Domain metric.
    #[must_use]
    pub const fn domain(&self) -> SpectralMetricV1 {
        self.domain
    }

    /// Codomain metric.
    #[must_use]
    pub const fn codomain(&self) -> SpectralMetricV1 {
        self.codomain
    }

    /// Gauge convention.
    #[must_use]
    pub const fn gauge(&self) -> GaugeConventionV1 {
        self.gauge
    }

    /// Zero-padding convention.
    #[must_use]
    pub const fn zero_padding(&self) -> ZeroPaddingConventionV1 {
        self.zero_padding
    }
}

/// Unit and normalization boundary used before numerical kernels consume
/// dimensionless samples.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralScalingContextV1 {
    id: SpectralScalingId,
    spectral_dims: Dims,
    spectral_scale_si: f64,
    left_map: SpectralScalingMapId,
    right_map: SpectralScalingMapId,
    operator_map: SpectralScalingMapId,
    inverse_map: SpectralScalingMapId,
}

impl SpectralScalingContextV1 {
    /// Construct a complete scaling descriptor. The scale is the positive SI
    /// value represented by one dimensionless spectral unit.
    #[allow(clippy::too_many_arguments)] // One complete normalization bundle; no ambient map defaults.
    #[must_use]
    pub const fn new(
        id: SpectralScalingId,
        spectral_dims: Dims,
        spectral_scale_si: f64,
        left_map: SpectralScalingMapId,
        right_map: SpectralScalingMapId,
        operator_map: SpectralScalingMapId,
        inverse_map: SpectralScalingMapId,
    ) -> Self {
        Self {
            id,
            spectral_dims,
            spectral_scale_si,
            left_map,
            right_map,
            operator_map,
            inverse_map,
        }
    }

    /// Scaling-bundle identity.
    #[must_use]
    pub const fn id(&self) -> SpectralScalingId {
        self.id
    }

    /// Left-space normalization map identity.
    #[must_use]
    pub const fn left_map(&self) -> SpectralScalingMapId {
        self.left_map
    }

    /// Right-space normalization map identity.
    #[must_use]
    pub const fn right_map(&self) -> SpectralScalingMapId {
        self.right_map
    }

    /// Operator normalization map identity.
    #[must_use]
    pub const fn operator_map(&self) -> SpectralScalingMapId {
        self.operator_map
    }

    /// Inverse/result normalization map identity.
    #[must_use]
    pub const fn inverse_map(&self) -> SpectralScalingMapId {
        self.inverse_map
    }

    /// Six-base SI dimensions `[m, kg, s, K, A, mol]` of spectral values.
    #[must_use]
    pub const fn spectral_dims(&self) -> Dims {
        self.spectral_dims
    }

    /// Positive canonical SI normalization scale.
    #[must_use]
    pub const fn spectral_scale_si(&self) -> f64 {
        self.spectral_scale_si
    }

    /// Normalize a finite runtime-dimensioned SI value for a dimensionless
    /// numerical kernel.
    ///
    /// # Errors
    ///
    /// Returns a structured issue when the value has the wrong dimensions,
    /// either operand is non-finite or non-positive where required, or the
    /// normalization loses a nonzero value to underflow.
    #[must_use = "normalization failures must be handled"]
    pub fn normalize(&self, value_si: QtyAny) -> Result<f64, SpectralAdmissionIssueV1> {
        if value_si.dims != self.spectral_dims {
            return Err(SpectralAdmissionIssueV1::UnitMismatch {
                expected: self.spectral_dims,
                found: value_si.dims,
            });
        }
        if !value_si.value.is_finite() {
            return Err(SpectralAdmissionIssueV1::NonFinite {
                field: AdmissionFieldV1::SpectralValue,
            });
        }
        if !(self.spectral_scale_si.is_finite() && self.spectral_scale_si > 0.0) {
            return Err(SpectralAdmissionIssueV1::NonPositive {
                field: AdmissionFieldV1::SpectralScale,
            });
        }
        let normalized = value_si.value / self.spectral_scale_si;
        if !normalized.is_finite() {
            return Err(SpectralAdmissionIssueV1::NonFinite {
                field: AdmissionFieldV1::NormalizedSpectralValue,
            });
        }
        if value_si.value != 0.0 && normalized == 0.0 {
            return Err(SpectralAdmissionIssueV1::Underflow {
                field: AdmissionFieldV1::NormalizedSpectralValue,
            });
        }
        Ok(if normalized == 0.0 { 0.0 } else { normalized })
    }

    /// Map a dimensionless numerical result back to a runtime-dimensioned
    /// coherent-SI value.
    ///
    /// # Errors
    ///
    /// Returns a structured issue when the normalized value or configured
    /// scale is invalid, the SI result is non-finite, or a nonzero value
    /// underflows to zero.
    #[must_use = "denormalization failures must be handled"]
    pub fn denormalize(&self, value: f64) -> Result<QtyAny, SpectralAdmissionIssueV1> {
        if !value.is_finite() {
            return Err(SpectralAdmissionIssueV1::NonFinite {
                field: AdmissionFieldV1::NormalizedSpectralValue,
            });
        }
        if !(self.spectral_scale_si.is_finite() && self.spectral_scale_si > 0.0) {
            return Err(SpectralAdmissionIssueV1::NonPositive {
                field: AdmissionFieldV1::SpectralScale,
            });
        }
        let value_si = value * self.spectral_scale_si;
        if !value_si.is_finite() {
            return Err(SpectralAdmissionIssueV1::NonFinite {
                field: AdmissionFieldV1::SpectralValue,
            });
        }
        if value != 0.0 && value_si == 0.0 {
            return Err(SpectralAdmissionIssueV1::Underflow {
                field: AdmissionFieldV1::SpectralValue,
            });
        }
        Ok(QtyAny::new(
            if value_si == 0.0 { 0.0 } else { value_si },
            self.spectral_dims,
        ))
    }
}

/// Regularity proposition required by an algorithm family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RegularityClassV1 {
    /// Finite-dimensional bounded operator.
    FiniteDimensional,
    /// Regular generalized pencil.
    RegularPencil,
    /// The generalized-pencil weight is invertible, excluding projective
    /// eigenvalues from an `Ordinary` pencil classification.
    InvertiblePencilWeight,
    /// Regular matrix polynomial of the declared grade.
    RegularPolynomial {
        /// Exact polynomial grade to which the regularity evidence applies.
        grade: u32,
    },
    /// The exact-grade polynomial leading coefficient is invertible, excluding
    /// projective roots from an `Ordinary` polynomial classification.
    InvertiblePolynomialLeadingCoefficient {
        /// Polynomial grade whose leading coefficient was checked.
        grade: u32,
    },
    /// Descriptor problem with an admitted regular/singular decomposition.
    RegularDescriptor,
    /// Well-posed periodic evolution and monodromy construction.
    WellPosedMonodromy,
    /// Analytic operator function on the admitted parameter domain/branch.
    AnalyticOperatorFunction,
}

impl RegularityClassV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::FiniteDimensional => 0,
            Self::RegularPencil => 1,
            Self::InvertiblePencilWeight => 2,
            Self::RegularPolynomial { .. } => 3,
            Self::InvertiblePolynomialLeadingCoefficient { .. } => 4,
            Self::RegularDescriptor => 5,
            Self::WellPosedMonodromy => 6,
            Self::AnalyticOperatorFunction => 7,
        }
    }
}

/// One supported or contradicted regularity proposition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegularityClaimV1 {
    class: RegularityClassV1,
    disposition: WitnessDispositionV1,
    witness: AdmittedSpectralWitnessV1,
}

impl RegularityClaimV1 {
    /// Construct one explicit regularity claim.
    #[must_use]
    pub const fn new(
        class: RegularityClassV1,
        disposition: WitnessDispositionV1,
        witness: AdmittedSpectralWitnessV1,
    ) -> Self {
        Self {
            class,
            disposition,
            witness,
        }
    }

    /// Named proposition.
    #[must_use]
    pub const fn class(&self) -> RegularityClassV1 {
        self.class
    }

    /// Supporting versus contradicting disposition.
    #[must_use]
    pub const fn disposition(&self) -> WitnessDispositionV1 {
        self.disposition
    }

    /// Bound evidence.
    #[must_use]
    pub const fn witness(&self) -> &AdmittedSpectralWitnessV1 {
        &self.witness
    }
}

/// Product regularity profile. An empty profile is explicit `Unknown` at the
/// method-obligation boundary. Multiple independent claims are required for
/// combinations such as a regular descriptor matrix polynomial.
#[derive(Debug, Clone, PartialEq)]
pub struct RegularityProfileV1 {
    claims: Vec<RegularityClaimV1>,
}

impl RegularityProfileV1 {
    /// Construct a raw regularity profile.
    #[must_use]
    pub fn new(claims: Vec<RegularityClaimV1>) -> Self {
        Self { claims }
    }

    /// Raw claims. Validated descriptors expose canonical order separately.
    #[must_use]
    pub fn claims(&self) -> &[RegularityClaimV1] {
        &self.claims
    }
}

/// Deterministic secondary order for equal primary magnitude/shift distance.
/// The final cluster-lineage comparison makes the order total even when both
/// coordinates coincide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComplexTieBreakV1 {
    /// Compare normalized real coordinate, then imaginary coordinate, then
    /// stable cluster lineage.
    RealThenImagThenLineage,
    /// Compare normalized imaginary coordinate, then real coordinate, then
    /// stable cluster lineage.
    ImagThenRealThenLineage,
}

impl ComplexTieBreakV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::RealThenImagThenLineage => 0,
            Self::ImagThenRealThenLineage => 1,
        }
    }
}

/// Placement of the projective point at infinity in an explicitly charted
/// total order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProjectiveInfinityPlacementV1 {
    /// Infinity precedes all finite chart coordinates.
    First,
    /// Infinity follows all finite chart coordinates.
    Last,
}

impl ProjectiveInfinityPlacementV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::First => 0,
            Self::Last => 1,
        }
    }
}

/// Deterministic spectral ordering/target convention. Every nontrivial complex
/// prefix order includes its tie semantics in the problem identity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpectralOrderingV1 {
    /// No total order is claimed; results are set-valued clusters.
    SetValued,
    /// Ascending real value; valid only under a real-spectrum obligation.
    RealAscending,
    /// Descending real value; valid only under a real-spectrum obligation.
    RealDescending,
    /// Ascending complex magnitude with an explicit total secondary order.
    MagnitudeAscending {
        /// Tie policy for equal magnitudes.
        tie_break: ComplexTieBreakV1,
    },
    /// Nearest to a finite dimensionless complex shift.
    NearestShift {
        /// Real part of the normalized shift.
        real: f64,
        /// Imaginary part of the normalized shift.
        imag: f64,
        /// Tie policy for equal distances.
        tie_break: ComplexTieBreakV1,
    },
    /// Membership in a named region; boundary policy lives in the scope.
    NamedRegion {
        /// Identity of the target region artifact.
        region: SpectralRegionId,
    },
    /// Explicitly charted projective ordering including infinity.
    Projective {
        /// Content identity of the admitted projective chart/order artifact.
        chart: SpectralProjectiveChartId,
        /// Where the projective point at infinity appears.
        infinity: ProjectiveInfinityPlacementV1,
        /// Secondary order for equal finite chart radii/coordinates.
        tie_break: ComplexTieBreakV1,
    },
}

impl SpectralOrderingV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::SetValued => 0,
            Self::RealAscending => 1,
            Self::RealDescending => 2,
            Self::MagnitudeAscending { .. } => 3,
            Self::NearestShift { .. } => 4,
            Self::NamedRegion { .. } => 5,
            Self::Projective { .. } => 6,
        }
    }
}

/// Boundary membership convention for region-completeness claims.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RegionBoundaryPolicyV1 {
    /// Include the boundary.
    Closed,
    /// Exclude the boundary.
    Open,
    /// Boundary intersections make membership unresolved.
    RefuseIntersection,
}

impl RegionBoundaryPolicyV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Closed => 0,
            Self::Open => 1,
            Self::RefuseIntersection => 2,
        }
    }
}

/// Requested coverage/completeness scope. This is a request, never evidence
/// that the returned spectrum is complete.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompletenessScopeV1 {
    /// Candidate generation only; no completeness request.
    CandidateOnly,
    /// A bounded algebraic prefix under the declared ordering/target. Repeated
    /// boundary clusters are returned whole, so the achieved result may
    /// contain more than this requested algebraic cardinality.
    Partial {
        /// Ordered algebraic cardinality requested from the solver.
        requested: u32,
    },
    /// All eigenvalues in a named region.
    Region {
        /// Region identity.
        region: SpectralRegionId,
        /// Deterministic boundary membership policy.
        boundary: RegionBoundaryPolicyV1,
    },
    /// Complete finite/projective spectrum of the declared dimension.
    FullFinite {
        /// Expected total algebraic cardinality, including projective roots.
        /// This is `n` for a linear `n x n` problem and normally `grade * n`
        /// for a regular grade-`grade` matrix polynomial. When a quotient has
        /// been applied, `n` is the declared dimension of the already-induced
        /// post-reduction metric space; the pre-reduction nullity is provenance
        /// and is not subtracted again.
        algebraic_cardinality: u32,
        /// Infinite-eigenvalue handling.
        infinity_policy: InfiniteEigenvaluePolicyV1,
    },
}

impl CompletenessScopeV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::CandidateOnly => 0,
            Self::Partial { .. } => 1,
            Self::Region { .. } => 2,
            Self::FullFinite { .. } => 3,
        }
    }
}

/// Raw versioned problem descriptor. It is intentionally not an admitted
/// token: retained or agent-produced input must pass [`validate_problem`].
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralProblemSpecV1 {
    schema_version: u32,
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    structures: StructureProfileV1,
    scaling: SpectralScalingContextV1,
    spaces: SpectralSpaceContextV1,
    regularity: RegularityProfileV1,
    ordering: SpectralOrderingV1,
    requested_scope: CompletenessScopeV1,
}

impl SpectralProblemSpecV1 {
    /// Construct a current-version raw descriptor.
    #[allow(clippy::too_many_arguments)] // The Five Explicits and all orthogonal problem axes stay visible.
    #[must_use]
    pub fn new(
        subject: SpectralSubjectId,
        scalar_field: SpectralScalarFieldV1,
        class: SpectralProblemClassV1,
        structures: StructureProfileV1,
        scaling: SpectralScalingContextV1,
        spaces: SpectralSpaceContextV1,
        regularity: RegularityProfileV1,
        ordering: SpectralOrderingV1,
        requested_scope: CompletenessScopeV1,
    ) -> Self {
        Self::with_schema_version(
            SPECTRAL_PROBLEM_SCHEMA_VERSION,
            subject,
            scalar_field,
            class,
            structures,
            scaling,
            spaces,
            regularity,
            ordering,
            requested_scope,
        )
    }

    /// Construct decoded versioned input. Unknown versions remain raw and
    /// fail closed in [`validate_problem`].
    #[allow(clippy::too_many_arguments)] // Decoded schema version plus every identity-bearing problem axis.
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        subject: SpectralSubjectId,
        scalar_field: SpectralScalarFieldV1,
        class: SpectralProblemClassV1,
        structures: StructureProfileV1,
        scaling: SpectralScalingContextV1,
        spaces: SpectralSpaceContextV1,
        regularity: RegularityProfileV1,
        ordering: SpectralOrderingV1,
        requested_scope: CompletenessScopeV1,
    ) -> Self {
        Self {
            schema_version,
            subject,
            scalar_field,
            class,
            structures,
            scaling,
            spaces,
            regularity,
            ordering,
            requested_scope,
        }
    }

    /// Declared schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact subject identity.
    #[must_use]
    pub const fn subject(&self) -> SpectralSubjectId {
        self.subject
    }

    /// Scalar field of the operator spaces.
    #[must_use]
    pub const fn scalar_field(&self) -> SpectralScalarFieldV1 {
        self.scalar_field
    }

    /// Class product.
    #[must_use]
    pub const fn class(&self) -> SpectralProblemClassV1 {
        self.class
    }

    /// Structure profile.
    #[must_use]
    pub const fn structures(&self) -> &StructureProfileV1 {
        &self.structures
    }

    /// Unit/scaling context.
    #[must_use]
    pub const fn scaling(&self) -> SpectralScalingContextV1 {
        self.scaling
    }

    /// Space context.
    #[must_use]
    pub const fn spaces(&self) -> SpectralSpaceContextV1 {
        self.spaces
    }

    /// Regularity state.
    #[must_use]
    pub const fn regularity(&self) -> &RegularityProfileV1 {
        &self.regularity
    }

    /// Requested deterministic ordering/target convention.
    #[must_use]
    pub const fn ordering(&self) -> SpectralOrderingV1 {
        self.ordering
    }

    /// Requested completeness scope.
    #[must_use]
    pub const fn requested_scope(&self) -> CompletenessScopeV1 {
        self.requested_scope
    }
}

/// Field associated with a structured admission issue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AdmissionFieldV1 {
    /// Schema version.
    SchemaVersion,
    /// Polynomial grade.
    PolynomialGrade,
    /// Domain dimension.
    DomainDimension,
    /// Codomain dimension.
    CodomainDimension,
    /// Spectral normalization scale.
    SpectralScale,
    /// Physical spectral value.
    SpectralValue,
    /// Dimensionless normalized spectral value.
    NormalizedSpectralValue,
    /// Floquet period.
    FloquetPeriod,
    /// Floquet branch anchor.
    FloquetBranch,
    /// Nearest-shift target.
    OrderingShift,
    /// Metric bounds/signature.
    Metric,
    /// Witness verifier version.
    WitnessVersion,
    /// Witness tolerance.
    WitnessTolerance,
    /// Requested completeness count/dimension.
    CompletenessScope,
    /// Certified gauge/nullity count.
    GaugeNullity,
    /// Structural zero-padding count.
    ZeroPaddingCount,
}

impl AdmissionFieldV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::SchemaVersion => 0,
            Self::PolynomialGrade => 1,
            Self::DomainDimension => 2,
            Self::CodomainDimension => 3,
            Self::SpectralScale => 4,
            Self::SpectralValue => 5,
            Self::NormalizedSpectralValue => 6,
            Self::FloquetPeriod => 7,
            Self::FloquetBranch => 8,
            Self::OrderingShift => 9,
            Self::Metric => 10,
            Self::WitnessVersion => 11,
            Self::WitnessTolerance => 12,
            Self::CompletenessScope => 13,
            Self::GaugeNullity => 14,
            Self::ZeroPaddingCount => 15,
        }
    }
}

/// Bounded raw profile whose item limit was exceeded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ClaimProfileV1 {
    /// Mathematical structure claims.
    Structure,
    /// Regularity claims.
    Regularity,
}

impl ClaimProfileV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Structure => 0,
            Self::Regularity => 1,
        }
    }
}

/// Structured fail-closed admission or method-obligation issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpectralAdmissionIssueV1 {
    /// Unsupported problem schema.
    UnsupportedSchemaVersion {
        /// Schema version supplied by the raw descriptor.
        found: u32,
        /// Sole schema version understood by this admission implementation.
        supported: u32,
    },
    /// Untrusted input exceeded its pre-processing item budget.
    TooManyClaims {
        /// Bounded profile whose limit was exceeded.
        profile: ClaimProfileV1,
        /// Number of claims supplied by the caller.
        found: usize,
        /// Maximum number admitted before sorting or pairwise checks.
        limit: usize,
    },
    /// Zero is invalid for the named dimension/count field.
    Zero {
        /// Field whose value was zero.
        field: AdmissionFieldV1,
    },
    /// Non-finite value.
    NonFinite {
        /// Field containing a NaN or infinity.
        field: AdmissionFieldV1,
    },
    /// A nonzero finite value rounded to zero at a normalization boundary.
    Underflow {
        /// Field whose nonzero value underflowed.
        field: AdmissionFieldV1,
    },
    /// Value must be finite and strictly positive.
    NonPositive {
        /// Field that failed the finite, strictly-positive requirement.
        field: AdmissionFieldV1,
    },
    /// Runtime physical dimensions do not match the declared spectral units.
    UnitMismatch {
        /// Dimensions declared by the spectral scaling context.
        expected: Dims,
        /// Dimensions carried by the supplied runtime quantity.
        found: Dims,
    },
    /// Domain/codomain or scope dimension mismatch.
    DimensionMismatch {
        /// First dimension or count participating in the failed equality or
        /// upper-bound check.
        left: u32,
        /// Required or comparison dimension/count.
        right: u32,
    },
    /// Polynomial grade times operator dimension exceeded the schema range.
    AlgebraicCardinalityOverflow {
        /// Operator-space dimension.
        dimension: u32,
        /// Polynomial grade used in the cardinality product.
        grade: u32,
    },
    /// Metric bounds/signature are internally inconsistent.
    InvalidMetric {
        /// Identity of the invalid metric descriptor.
        metric: SpectralMetricId,
    },
    /// One content identity was attached to two distinct metric descriptors.
    MetricIdentityConflict {
        /// Rebound metric identity.
        metric: SpectralMetricId,
    },
    /// Representation is incompatible with another class axis.
    RepresentationConflict,
    /// Floquet parameter dimensions or branch semantics are inconsistent.
    FloquetSemanticMismatch,
    /// Ordering is unavailable for the declared scalar/class context.
    OrderingUnavailable,
    /// Ordering target and requested completeness scope disagree.
    ScopeOrderingMismatch,
    /// Descriptor and requested-result infinity policies disagree.
    InfinityPolicyMismatch,
    /// A partial prefix that can contain projective infinity requires an
    /// explicitly charted projective ordering and infinity placement.
    ProjectivePrefixOrderingRequired,
    /// Structure claim appears twice with the same normalized property,
    /// support, norm, and tolerance key and the same disposition.
    DuplicateStructure {
        /// Duplicated mathematical property.
        property: StructurePropertyV1,
        /// Duplicated metric/form support.
        support: StructureSupportV1,
    },
    /// The same normalized proposition is both supported and contradicted.
    ContradictoryStructure {
        /// Mathematical property with opposed dispositions.
        property: StructurePropertyV1,
        /// Metric/form support on which the contradiction occurs.
        support: StructureSupportV1,
    },
    /// Exact normality and exact nonnormality are logical complements, yet
    /// both were explicitly contradicted on the same admitted support.
    ComplementaryStructureConflict {
        /// Shared inner product on which both complement propositions were
        /// contradicted.
        support: StructureSupportV1,
    },
    /// Exact theorem closure from one admitted structure proposition conflicts
    /// with another admitted proposition.
    StructureTheoremConflict {
        /// Exact witnessed premise.
        premise: StructurePropertyV1,
        /// Consequence contradicted by the other claim.
        consequence: StructurePropertyV1,
        /// Geometric support on which the implication is valid.
        support: StructureSupportV1,
    },
    /// An exact structure theorem implies a regularity proposition that the
    /// same admitted profile explicitly contradicts.
    StructureRegularityTheoremConflict {
        /// Exact witnessed structure premise.
        premise: StructurePropertyV1,
        /// Contradicted regularity consequence.
        consequence: RegularityClassV1,
        /// Geometric support on which the implication is valid.
        support: StructureSupportV1,
    },
    /// One admitted regularity proposition implies another that the same
    /// profile explicitly contradicts.
    RegularityTheoremConflict {
        /// Witnessed regularity premise.
        premise: RegularityClassV1,
        /// Contradicted regularity consequence.
        consequence: RegularityClassV1,
    },
    /// Structure proposition was attached to the wrong kind of metric/form.
    InvalidStructureSupport {
        /// Property whose support contract was violated.
        property: StructurePropertyV1,
        /// Supplied metric/form support.
        support: StructureSupportV1,
    },
    /// Structure proposition is undefined for the declared equation
    /// representation (for example, a palindromic polynomial claim on a
    /// standard-linear problem).
    StructureRepresentationMismatch {
        /// Inapplicable mathematical property.
        property: StructurePropertyV1,
        /// Declared equation representation.
        representation: SpectralRepresentationV1,
    },
    /// Structure tolerance was negative or non-finite.
    InvalidStructureTolerance {
        /// Property carrying the invalid tolerance.
        property: StructurePropertyV1,
    },
    /// Admitted evidence names a different exact proposition.
    WitnessPropositionMismatch {
        /// Proposition required by the enclosing claim.
        expected: SpectralPropositionId,
        /// Proposition actually admitted by the verifier/policy.
        found: SpectralPropositionId,
    },
    /// Typed proposition digest matched but the independently retained
    /// canonical-preimage observation differed.
    WitnessObservationMismatch {
        /// Digest whose retained canonical-preimage root did not match.
        proposition: SpectralPropositionId,
    },
    /// Regularity class conflicts with the explicit finite-space schema,
    /// representation, role, or origin.
    RegularityMismatch,
    /// An ordinary pencil/polynomial classification lacks proposition-backed
    /// exclusion of projective/infinite eigenvalues.
    OrdinaryFiniteSpectrumWitnessRequired {
        /// Exact invertibility proposition required by the representation.
        required: RegularityClassV1,
    },
    /// Method requires an equation representation not supplied by the problem.
    MethodRepresentationMismatch {
        /// Method family whose representation obligation failed.
        method: SpectralMethodClassV1,
    },
    /// Method requires a semantic origin not supplied by the problem.
    MethodOriginMismatch {
        /// Method family whose origin obligation failed.
        method: SpectralMethodClassV1,
    },
    /// Method cannot accept the descriptor role.
    MethodDescriptorMismatch {
        /// Method family whose descriptor-role obligation failed.
        method: SpectralMethodClassV1,
    },
    /// Method requires identical domain/codomain spaces and metrics.
    MethodSpaceMismatch {
        /// Method family whose space obligation failed.
        method: SpectralMethodClassV1,
    },
    /// A symplectic/Hamiltonian method requires even dimension.
    EvenDimensionRequired {
        /// Structure-preserving method family requiring even dimension.
        method: SpectralMethodClassV1,
    },
    /// Descriptor execution requires an explicit projective/infinity policy.
    DescriptorInfinityPolicyRequired {
        /// Descriptor-capable method family requiring the policy.
        method: SpectralMethodClassV1,
    },
    /// Required structure witness is absent in the method's metric.
    MissingStructureWitness {
        /// Method family whose structure obligation is unsatisfied.
        method: SpectralMethodClassV1,
        /// Required mathematical property.
        property: StructurePropertyV1,
        /// Required exact or form-family support.
        support: StructureSupportRequirementV1,
    },
    /// Available witness explicitly contradicts the method obligation.
    ContradictedMethodObligation {
        /// Refused method family.
        method: SpectralMethodClassV1,
        /// Explicitly contradicted property.
        property: StructurePropertyV1,
    },
    /// V1 structure-preserving methods require a zero-tolerance proposition;
    /// approximate claims remain available to routers but cannot authorize the
    /// specialized kernel family without an explicit error budget API.
    ExactStructureWitnessRequired {
        /// Method family requiring an exact proposition.
        method: SpectralMethodClassV1,
        /// Property available only with a nonzero tolerance or otherwise
        /// lacking an exact uncontradicted witness.
        property: StructurePropertyV1,
    },
    /// Method requires a positive-definite metric.
    PositiveDefiniteMetricRequired {
        /// Method family requiring the positive-definite metric.
        method: SpectralMethodClassV1,
    },
    /// Method requires a nondegenerate indefinite/Krein metric.
    IndefiniteMetricRequired {
        /// Method family requiring the indefinite metric.
        method: SpectralMethodClassV1,
    },
    /// Required regularity proposition is absent, contradicted, or mismatched.
    RegularityRequired {
        /// Method family whose regularity obligation failed.
        method: SpectralMethodClassV1,
        /// Exact regularity proposition required by that method.
        required: RegularityClassV1,
    },
    /// Structural-zero-sensitive gap interpretation requires an explicit gauge
    /// convention, independently of eigensolver method admission.
    GapGaugeConventionRequired,
    /// Structural-zero-sensitive gap interpretation requires an explicit
    /// zero-padding/omission convention.
    GapZeroPaddingConventionRequired,
    /// The serialized structural-zero count contradicts the certified gauge
    /// nullity. A fixed/quotiented gauge may serialize no padding at all, but
    /// any explicit padded/omitted count must equal the certified nullity.
    GapStructuralZeroCountMismatch {
        /// Nullity certified by the gauge proposition.
        gauge_nullity: u32,
        /// Count declared by the padding/omission proposition.
        declared_zero_count: u32,
    },
    /// Canonical typed identity construction failed closed.
    Identity(CanonicalError),
}

impl SpectralAdmissionIssueV1 {
    fn sort_key(&self) -> (u16, u16, u16) {
        match self {
            Self::UnsupportedSchemaVersion { .. } => (0, 0, 0),
            Self::TooManyClaims { profile, .. } => (0, 1, u16::from(profile.tag())),
            Self::Zero { field } => (1, u16::from(field.tag()), 0),
            Self::NonFinite { field } => (2, u16::from(field.tag()), 0),
            Self::Underflow { field } => (2, u16::from(field.tag()), 1),
            Self::NonPositive { field } => (3, u16::from(field.tag()), 0),
            Self::UnitMismatch { .. } => (4, 0, 0),
            Self::DimensionMismatch { .. } => (5, 0, 0),
            Self::AlgebraicCardinalityOverflow { .. } => (6, 0, 0),
            Self::InvalidMetric { .. } => (7, 0, 0),
            Self::MetricIdentityConflict { .. } => (8, 0, 0),
            Self::RepresentationConflict => (9, 0, 0),
            Self::FloquetSemanticMismatch => (10, 0, 0),
            Self::OrderingUnavailable => (11, 0, 0),
            Self::ScopeOrderingMismatch => (12, 0, 0),
            Self::InfinityPolicyMismatch => (13, 0, 0),
            Self::ProjectivePrefixOrderingRequired => (13, 0, 1),
            Self::DuplicateStructure { property, .. } => (14, u16::from(property.tag()), 0),
            Self::ContradictoryStructure { property, .. } => (15, u16::from(property.tag()), 0),
            Self::ComplementaryStructureConflict { .. } => {
                (15, u16::from(StructurePropertyV1::Normal.tag()), 1)
            }
            Self::StructureTheoremConflict {
                premise,
                consequence,
                ..
            } => (16, u16::from(premise.tag()), u16::from(consequence.tag())),
            Self::StructureRegularityTheoremConflict {
                premise,
                consequence,
                ..
            } => (
                16,
                64 + u16::from(premise.tag()),
                u16::from(consequence.tag()),
            ),
            Self::RegularityTheoremConflict {
                premise,
                consequence,
            } => (
                16,
                128 + u16::from(premise.tag()),
                u16::from(consequence.tag()),
            ),
            Self::InvalidStructureSupport { property, .. } => (17, u16::from(property.tag()), 0),
            Self::StructureRepresentationMismatch { property, .. } => {
                (17, u16::from(property.tag()), 1)
            }
            Self::InvalidStructureTolerance { property } => (18, u16::from(property.tag()), 0),
            Self::WitnessPropositionMismatch { .. } => (19, 0, 0),
            Self::WitnessObservationMismatch { .. } => (20, 0, 0),
            Self::RegularityMismatch => (21, 0, 0),
            Self::OrdinaryFiniteSpectrumWitnessRequired { required } => {
                (22, u16::from(required.tag()), 0)
            }
            Self::MethodRepresentationMismatch { method } => (23, u16::from(method.tag()), 0),
            Self::MethodOriginMismatch { method } => (24, u16::from(method.tag()), 0),
            Self::MethodDescriptorMismatch { method } => (25, u16::from(method.tag()), 0),
            Self::MethodSpaceMismatch { method } => (26, u16::from(method.tag()), 0),
            Self::EvenDimensionRequired { method } => (27, u16::from(method.tag()), 0),
            Self::DescriptorInfinityPolicyRequired { method } => (28, u16::from(method.tag()), 0),
            Self::MissingStructureWitness {
                method, property, ..
            } => (29, u16::from(method.tag()), u16::from(property.tag())),
            Self::ContradictedMethodObligation { method, property } => {
                (30, u16::from(method.tag()), u16::from(property.tag()))
            }
            Self::ExactStructureWitnessRequired { method, property } => {
                (31, u16::from(method.tag()), u16::from(property.tag()))
            }
            Self::PositiveDefiniteMetricRequired { method } => (32, u16::from(method.tag()), 0),
            Self::IndefiniteMetricRequired { method } => (33, u16::from(method.tag()), 0),
            Self::RegularityRequired { method, required } => {
                (34, u16::from(method.tag()), u16::from(required.tag()))
            }
            Self::GapGaugeConventionRequired => (35, 0, 0),
            Self::GapZeroPaddingConventionRequired => (36, 0, 0),
            Self::GapStructuralZeroCountMismatch { .. } => (37, 0, 0),
            Self::Identity(_) => (38, 0, 0),
        }
    }
}

impl fmt::Display for SpectralAdmissionIssueV1 {
    #[allow(clippy::too_many_lines)] // One exhaustive, auditable diagnostic mapping for the public issue enum.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSchemaVersion { found, supported } => write!(
                f,
                "unsupported spectral schema version {found}; supported version is {supported}"
            ),
            Self::TooManyClaims {
                profile,
                found,
                limit,
            } => write!(
                f,
                "{profile:?} profile has {found} claims; v1 admits at most {limit}"
            ),
            Self::Zero { field } => write!(f, "{field:?} must be nonzero"),
            Self::NonFinite { field } => write!(f, "{field:?} must be finite"),
            Self::Underflow { field } => write!(
                f,
                "nonzero {field:?} underflowed to zero at the spectral normalization boundary"
            ),
            Self::NonPositive { field } => {
                write!(f, "{field:?} must be finite and strictly positive")
            }
            Self::UnitMismatch { expected, found } => write!(
                f,
                "spectral value units {} do not match expected {}",
                found.unit_string(),
                expected.unit_string()
            ),
            Self::DimensionMismatch { left, right } => {
                write!(f, "spectral dimensions differ: {left} versus {right}")
            }
            Self::AlgebraicCardinalityOverflow { dimension, grade } => write!(
                f,
                "algebraic cardinality {dimension} * {grade} exceeds the v1 schema range"
            ),
            Self::InvalidMetric { metric } => write!(f, "metric {metric:?} is inconsistent"),
            Self::MetricIdentityConflict { metric } => write!(
                f,
                "metric identity {metric:?} was rebound to a different descriptor"
            ),
            Self::RepresentationConflict => {
                f.write_str("spectral representation conflicts with descriptor/origin semantics")
            }
            Self::FloquetSemanticMismatch => {
                f.write_str("Floquet period, output units, or logarithm branch is inconsistent")
            }
            Self::OrderingUnavailable => {
                f.write_str("requested ordering is undefined for this admitted problem")
            }
            Self::ScopeOrderingMismatch => {
                f.write_str("spectral ordering target conflicts with requested coverage scope")
            }
            Self::InfinityPolicyMismatch => f.write_str(
                "descriptor, projective ordering, and requested infinity policies disagree",
            ),
            Self::ProjectivePrefixOrderingRequired => f.write_str(
                "a partial spectrum including projective infinity requires an explicit projective chart and infinity placement",
            ),
            Self::DuplicateStructure {
                property, support, ..
            } => write!(
                f,
                "duplicate {property:?} claim on {support:?} with the same normalized norm/tolerance key; combine evidence first"
            ),
            Self::ContradictoryStructure { property, support } => write!(
                f,
                "{property:?} is both supported and contradicted on {support:?}"
            ),
            Self::ComplementaryStructureConflict { support } => write!(
                f,
                "exact Normal and exact Nonnormal were both contradicted on {support:?}; logical complements cannot both be false"
            ),
            Self::StructureTheoremConflict {
                premise,
                consequence,
                support,
            } => write!(
                f,
                "exact {premise:?} on {support:?} implies {consequence:?}, but the admitted profile contradicts that consequence"
            ),
            Self::StructureRegularityTheoremConflict {
                premise,
                consequence,
                support,
            } => write!(
                f,
                "exact {premise:?} on {support:?} implies {consequence:?}, but the admitted regularity profile contradicts that consequence"
            ),
            Self::RegularityTheoremConflict {
                premise,
                consequence,
            } => write!(
                f,
                "witnessed {premise:?} implies {consequence:?}, but the admitted regularity profile contradicts that consequence"
            ),
            Self::InvalidStructureSupport { property, support } => write!(
                f,
                "{property:?} proposition cannot use {support:?} as its geometric support"
            ),
            Self::StructureRepresentationMismatch {
                property,
                representation,
            } => write!(
                f,
                "{property:?} proposition is not defined for {representation:?} equations"
            ),
            Self::InvalidStructureTolerance { property } => write!(
                f,
                "{property:?} structure tolerance must be finite and nonnegative"
            ),
            Self::WitnessPropositionMismatch { expected, found } => write!(
                f,
                "admitted witness proposition {} does not match required proposition {}",
                found.to_hex(),
                expected.to_hex()
            ),
            Self::WitnessObservationMismatch { proposition } => write!(
                f,
                "proposition {} was presented for different canonical observations",
                proposition.to_hex()
            ),
            Self::RegularityMismatch => f.write_str(
                "regularity state conflicts with the finite-space schema, representation, role, or origin",
            ),
            Self::OrdinaryFiniteSpectrumWitnessRequired { required } => write!(
                f,
                "ordinary pencil/polynomial semantics require witnessed {required:?} to exclude projective roots"
            ),
            Self::MethodRepresentationMismatch { method } => {
                write!(f, "{method:?} does not admit this equation representation")
            }
            Self::MethodOriginMismatch { method } => {
                write!(f, "{method:?} does not admit this operator origin")
            }
            Self::MethodDescriptorMismatch { method } => {
                write!(f, "{method:?} does not admit this descriptor role")
            }
            Self::MethodSpaceMismatch { method } => {
                write!(f, "{method:?} requires identical admitted operator spaces")
            }
            Self::EvenDimensionRequired { method } => {
                write!(f, "{method:?} requires an even-dimensional operator space")
            }
            Self::DescriptorInfinityPolicyRequired { method } => write!(
                f,
                "{method:?} requires an explicit descriptor infinity policy"
            ),
            Self::MissingStructureWitness {
                method,
                property,
                support,
            } => write!(
                f,
                "{method:?} requires a {property:?} witness on {support:?}"
            ),
            Self::ContradictedMethodObligation { method, property } => write!(
                f,
                "{method:?} is refused because {property:?} is contradicted"
            ),
            Self::ExactStructureWitnessRequired { method, property } => write!(
                f,
                "{method:?} requires a zero-tolerance {property:?} proposition in v1"
            ),
            Self::PositiveDefiniteMetricRequired { method } => {
                write!(f, "{method:?} requires a positive-definite metric")
            }
            Self::IndefiniteMetricRequired { method } => {
                write!(f, "{method:?} requires a nondegenerate indefinite metric")
            }
            Self::RegularityRequired { method, required } => {
                write!(f, "{method:?} requires {required:?} regularity evidence")
            }
            Self::GapGaugeConventionRequired => {
                f.write_str("spectral-gap interpretation requires an explicit gauge convention")
            }
            Self::GapZeroPaddingConventionRequired => f.write_str(
                "spectral-gap interpretation requires an explicit zero-padding convention",
            ),
            Self::GapStructuralZeroCountMismatch {
                gauge_nullity,
                declared_zero_count,
            } => write!(
                f,
                "spectral-gap zero convention declares {declared_zero_count} structural zeros but the gauge certifies nullity {gauge_nullity}"
            ),
            Self::Identity(error) => write!(f, "spectral identity construction failed: {error}"),
        }
    }
}

impl core::error::Error for SpectralAdmissionIssueV1 {}

/// Deterministically ranked collection of admission issues.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpectralAdmissionReportV1 {
    issues: Vec<SpectralAdmissionIssueV1>,
}

impl SpectralAdmissionReportV1 {
    fn new(mut issues: Vec<SpectralAdmissionIssueV1>) -> Self {
        issues.sort_by_cached_key(|issue| (issue.sort_key(), format!("{issue:?}")));
        issues.dedup();
        Self { issues }
    }

    /// Ranked actionable issues.
    #[must_use]
    pub fn issues(&self) -> &[SpectralAdmissionIssueV1] {
        &self.issues
    }
}

impl fmt::Display for SpectralAdmissionReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "spectral admission refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for SpectralAdmissionReportV1 {}

/// Algorithm family whose prerequisites can be assessed without selecting a
/// concrete implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SpectralMethodClassV1 {
    /// Standard positive-metric self-adjoint Lanczos.
    SelfAdjointLanczos,
    /// Generalized self-adjoint-definite Lanczos.
    GeneralizedSelfAdjointLanczos,
    /// General nonnormal Arnoldi.
    GeneralArnoldi,
    /// Polynomial/QEP Krylov method.
    PolynomialKrylov,
    /// Hamiltonian structure-preserving method.
    HamiltonianStructurePreserving,
    /// Symplectic map method.
    SymplecticStructurePreserving,
    /// J/Krein self-adjoint method.
    KreinJOrthogonal,
    /// Monodromy/Floquet Arnoldi.
    MonodromyArnoldi,
    /// Descriptor pencil method with projective semantics.
    DescriptorPencil,
    /// Analytic operator-function Krylov method.
    OperatorFunctionKrylov,
}

impl SpectralMethodClassV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::SelfAdjointLanczos => 0,
            Self::GeneralizedSelfAdjointLanczos => 1,
            Self::GeneralArnoldi => 2,
            Self::PolynomialKrylov => 3,
            Self::HamiltonianStructurePreserving => 4,
            Self::SymplecticStructurePreserving => 5,
            Self::KreinJOrthogonal => 6,
            Self::MonodromyArnoldi => 7,
            Self::DescriptorPencil => 8,
            Self::OperatorFunctionKrylov => 9,
        }
    }
}

/// Sealed, validated descriptor. Fields remain private so raw input cannot be
/// mistaken for admitted input.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedSpectralProblemV1 {
    spec: SpectralProblemSpecV1,
    canonical_claims: Vec<StructureClaimV1>,
    canonical_regularity: Vec<RegularityClaimV1>,
    problem_receipt: SpectralProblemIdentityReceiptV2,
}

impl ValidatedSpectralProblemV1 {
    /// Observational view of the validated raw descriptor. Detached fields are
    /// not authorization capabilities; downstream authority-bearing APIs must
    /// consume this complete validated problem token.
    #[must_use]
    pub const fn spec(&self) -> &SpectralProblemSpecV1 {
        &self.spec
    }

    /// Canonically ordered observational structure claims. A detached claim is
    /// not an authorization capability.
    #[must_use]
    pub fn structure_claims(&self) -> &[StructureClaimV1] {
        &self.canonical_claims
    }

    /// Canonically ordered observational regularity claims. A detached claim
    /// is not an authorization capability.
    #[must_use]
    pub fn regularity_claims(&self) -> &[RegularityClaimV1] {
        &self.canonical_regularity
    }

    /// Typed semantic identity binding every authority-bearing field.
    #[must_use]
    pub const fn problem_id(&self) -> SpectralProblemId {
        SpectralProblemId(self.problem_receipt.id())
    }

    /// Complete canonical producer receipt retained for ledger adjudication.
    ///
    /// Digest equality alone is not collision adjudication. Durable consumers
    /// should retain this observation (or an equivalent ledger record) when
    /// the admitted descriptor crosses a trust boundary.
    #[must_use]
    pub const fn identity_receipt(&self) -> SpectralProblemIdentityReceiptV2 {
        self.problem_receipt
    }

    /// Exact algebraic cardinality implied by the admitted finite-dimensional
    /// equation class, when regularity makes that number meaningful. Ordinary
    /// pencils/polynomials are finite by admission; descriptor problems expose
    /// a count only after both descriptor and equation regularity are admitted.
    #[must_use]
    pub fn known_algebraic_cardinality(&self) -> Option<u32> {
        admitted_algebraic_cardinality(
            &self.spec,
            &self.canonical_claims,
            &self.canonical_regularity,
        )
    }

    /// Whether admitted exact theorem closure excludes nonzero projective
    /// infinity multiplicity for this equation representation. For a
    /// generalized pencil this consumes an invertible weight directly or via
    /// an exact Hermitian-definite-pencil proposition; for a matrix polynomial
    /// it consumes an invertible exact-grade leading coefficient.
    #[must_use]
    pub fn projective_infinity_is_excluded(&self) -> bool {
        projective_infinity_is_excluded_by_theorem(
            &self.spec,
            &self.canonical_claims,
            &self.canonical_regularity,
        )
    }

    /// Whether admitted exact structure semantics require every certified
    /// finite enclosure to intersect the real axis. This includes an explicit
    /// real-spectrum proposition and the currently registered solid theorem
    /// consequences; it is independent of the requested display ordering.
    #[must_use]
    pub fn requires_real_spectrum_truth(&self) -> bool {
        real_spectrum_is_admitted(&self.spec, &self.canonical_claims)
    }
}

/// Proof that one method class's schema obligations were satisfied. This is
/// not a routing decision or a claim that an implementation converged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmittedSpectralMethodClassV1 {
    problem_id: SpectralProblemId,
    method: SpectralMethodClassV1,
    selected_support: Option<StructureSupportV1>,
}

impl AdmittedSpectralMethodClassV1 {
    /// Bound problem identity.
    #[must_use]
    pub const fn problem_id(&self) -> SpectralProblemId {
        self.problem_id
    }

    /// Admitted method family.
    #[must_use]
    pub const fn method(&self) -> SpectralMethodClassV1 {
        self.method
    }

    /// Exact metric/form selected for a structure-preserving method. Generic
    /// methods return `None`; a token never leaves an any-form obligation
    /// ambiguous for the numerical kernel.
    #[must_use]
    pub const fn selected_support(&self) -> Option<StructureSupportV1> {
        self.selected_support
    }
}

/// Proof that structural-zero-sensitive gap interpretation has explicit,
/// proposition-validated gauge and serialized-zero conventions. This token is
/// deliberately separate from eigensolver method admission: an algorithm can
/// compute candidates without being authorized to interpret the first gap as
/// a physical/mechanism separation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmittedSpectralGapSemanticsV1 {
    problem_id: SpectralProblemId,
    gauge: GaugeConventionV1,
    zero_padding: ZeroPaddingConventionV1,
}

impl AdmittedSpectralGapSemanticsV1 {
    /// Bound problem identity.
    #[must_use]
    pub const fn problem_id(&self) -> SpectralProblemId {
        self.problem_id
    }

    /// Exact admitted gauge/nullspace convention.
    #[must_use]
    pub const fn gauge(&self) -> GaugeConventionV1 {
        self.gauge
    }

    /// Exact admitted structural-zero serialization convention.
    #[must_use]
    pub const fn zero_padding(&self) -> ZeroPaddingConventionV1 {
        self.zero_padding
    }
}

/// Stable proposition-family tags shared with the result-truth module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SpectralPropositionKindV1 {
    MetricDefiniteness,
    Structure,
    Regularity,
    Gauge,
    ZeroPadding,
    ResultAuthority,
    Multiplicity,
    Separation,
    Completeness,
}

impl SpectralPropositionKindV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::MetricDefiniteness => 0,
            Self::Structure => 1,
            Self::Regularity => 2,
            Self::Gauge => 3,
            Self::ZeroPadding => 4,
            Self::ResultAuthority => 5,
            Self::Multiplicity => 6,
            Self::Separation => 7,
            Self::Completeness => 8,
        }
    }
}

pub(crate) fn spectral_proposition_receipt(
    kind: SpectralPropositionKindV1,
    payload: &[u8],
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    CanonicalEncoder::<SpectralPropositionId, _>::new(PROPOSITION_IDENTITY_LIMITS, NeverCancel)?
        .variant(Field::new(0, "kind"), kind.tag(), &[])?
        .bytes(Field::new(1, "payload"), payload)?
        .finish()
}

/// Canonical verifier identity receipt. The descriptor should identify an
/// immutable implementation, configuration, and checker contract.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the descriptor exceeds the bounded identity
/// schema or cannot be canonically encoded.
#[must_use = "the verifier receipt or canonicalization failure must be handled"]
pub fn spectral_verifier_receipt(
    descriptor: &[u8],
) -> Result<IdentityReceipt<SpectralAuthorityVerifierIdV1>, CanonicalError> {
    CanonicalEncoder::<SpectralAuthorityVerifierIdV1, _>::new(IDENTITY_LIMITS, NeverCancel)?
        .bytes(Field::new(0, "descriptor"), descriptor)?
        .finish()
}

/// Canonical authority-policy receipt. The descriptor should identify the
/// exact policy manifest, not a display label.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the descriptor exceeds the bounded identity
/// schema or cannot be canonically encoded.
#[must_use = "the authority-policy receipt or canonicalization failure must be handled"]
pub fn spectral_authority_policy_receipt(
    descriptor: &[u8],
) -> Result<IdentityReceipt<SpectralAuthorityPolicyIdV1>, CanonicalError> {
    CanonicalEncoder::<SpectralAuthorityPolicyIdV1, _>::new(IDENTITY_LIMITS, NeverCancel)?
        .bytes(Field::new(0, "descriptor"), descriptor)?
        .finish()
}

fn push_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn canonical_f64_bits(value: f64) -> u64 {
    if value == 0.0 {
        0.0f64.to_bits()
    } else {
        value.to_bits()
    }
}

fn push_gauge_context(out: &mut Vec<u8>, context: GaugeContextV1) {
    out.push(context.tag());
    match context {
        GaugeContextV1::Fixed { nullity, gauge } => {
            push_u32(out, nullity);
            out.extend_from_slice(gauge.as_bytes());
        }
        GaugeContextV1::Quotiented { nullity, quotient } => {
            push_u32(out, nullity);
            out.extend_from_slice(quotient.as_bytes());
        }
        GaugeContextV1::CertifiedNone | GaugeContextV1::Unknown => {}
    }
}

fn witness_bytes(witness: &AdmittedSpectralWitnessV1) -> Vec<u8> {
    let audit = witness.audit;
    let promotion = witness.promotion;
    let mut out = Vec::with_capacity(32 * 8 + 64);
    push_u64(
        &mut out,
        u64::try_from(SPECTRAL_PROMOTION_WITNESS_ENCODING_DOMAIN_V1.len())
            .expect("static witness encoding domain length fits u64"),
    );
    out.extend_from_slice(SPECTRAL_PROMOTION_WITNESS_ENCODING_DOMAIN_V1);
    out.extend_from_slice(witness.proposition.as_bytes());
    out.extend_from_slice(audit.canonical_preimage().as_bytes());
    push_u64(&mut out, audit.canonical_bytes());
    match audit.anchor() {
        Some(anchor) => {
            out.push(1);
            out.extend_from_slice(anchor.as_bytes());
        }
        None => out.push(0),
    }
    match audit.verifier() {
        Some(verifier) => {
            out.push(1);
            out.extend_from_slice(&verifier);
        }
        None => out.push(0),
    }
    match audit.key_policy() {
        Some(policy) => {
            out.push(1);
            out.extend_from_slice(&policy);
        }
        None => out.push(0),
    }
    push_u64(
        &mut out,
        u64::try_from(promotion.verifier_domain.len())
            .expect("static verifier domain length fits u64"),
    );
    out.extend_from_slice(promotion.verifier_domain.as_bytes());
    out.extend_from_slice(promotion.verifier_observation.content_id().as_bytes());
    push_u64(&mut out, promotion.verifier_observation.length());
    push_u64(
        &mut out,
        u64::try_from(promotion.key_policy_domain.len())
            .expect("static key-policy domain length fits u64"),
    );
    out.extend_from_slice(promotion.key_policy_domain.as_bytes());
    out.extend_from_slice(promotion.key_policy_observation.content_id().as_bytes());
    push_u64(&mut out, promotion.key_policy_observation.length());
    push_u64(
        &mut out,
        u64::try_from(promotion.context.len()).expect("static promotion context length fits u64"),
    );
    out.extend_from_slice(promotion.context.as_bytes());
    out
}

fn class_bytes(class: SpectralProblemClassV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(96);
    out.push(class.representation.tag());
    if let SpectralRepresentationV1::MatrixPolynomial { grade } = class.representation {
        push_u32(&mut out, grade);
    }
    out.push(class.descriptor.tag());
    if let DescriptorRoleV1::Descriptor { infinity_policy } = class.descriptor {
        out.push(infinity_policy.tag());
    }
    out.push(class.origin.tag());
    match class.origin {
        SpectralOperatorOriginV1::Direct => {}
        SpectralOperatorOriginV1::MonodromyFloquet {
            period,
            parameter,
            branch,
        } => {
            push_u64(&mut out, canonical_f64_bits(period.value()));
            out.push(parameter.tag());
            out.push(branch.tag());
            if let FloquetBranchConventionV1::ContinuousFrom {
                continuation,
                anchor_phase,
            } = branch
            {
                out.extend_from_slice(continuation.as_bytes());
                push_u64(&mut out, canonical_f64_bits(anchor_phase.value()));
            }
        }
        SpectralOperatorOriginV1::AnalyticOperatorFunction {
            function,
            branch_policy,
        } => {
            out.extend_from_slice(function.as_bytes());
            out.push(branch_policy.tag());
        }
    }
    out
}

fn push_structure_property(out: &mut Vec<u8>, property: StructurePropertyV1) {
    out.push(property.tag());
    if let StructurePropertyV1::Palindromic { parity, involution } = property {
        out.push(parity.tag());
        out.push(involution.tag());
    }
}

fn push_structure_support(out: &mut Vec<u8>, support: StructureSupportV1) {
    out.push(support.tag());
    match support {
        StructureSupportV1::InnerProduct(metric) => out.extend_from_slice(metric.as_bytes()),
        StructureSupportV1::SymplecticForm(form)
        | StructureSupportV1::KreinForm(form)
        | StructureSupportV1::Conjugation(form) => out.extend_from_slice(form.as_bytes()),
        StructureSupportV1::FormFree => {}
    }
}

fn push_space_signature(out: &mut Vec<u8>, domain: SpectralMetricV1, codomain: SpectralMetricV1) {
    out.extend_from_slice(domain.id.as_bytes());
    push_u32(out, domain.dimension);
    out.extend_from_slice(codomain.id.as_bytes());
    push_u32(out, codomain.dimension);
}

/// Canonical receipt for one exact structure proposition. Evidence producers
/// must verify this receipt's preimage and external anchor before admitting an
/// authority reference.
///
/// # Errors
///
/// Returns [`CanonicalError`] when `tolerance` is non-finite or the fully bound
/// proposition exceeds the canonical identity limits.
#[allow(clippy::too_many_arguments)] // Every argument is an independent part of the exact proposition identity.
#[must_use = "the structure-proposition receipt or canonicalization failure must be handled"]
pub fn structure_proposition_receipt(
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
    property: StructurePropertyV1,
    support: StructureSupportV1,
    disposition: WitnessDispositionV1,
    tolerance: f64,
    norm: SpectralNormId,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    if !tolerance.is_finite() {
        return Err(CanonicalError::NonFiniteFloat {
            bits: tolerance.to_bits(),
        });
    }
    let mut payload = Vec::with_capacity(256);
    payload.extend_from_slice(subject.as_bytes());
    payload.push(scalar_field.tag() as u8);
    let class = class_bytes(class);
    push_u32(&mut payload, class.len() as u32);
    payload.extend_from_slice(&class);
    let scaling = scaling_bytes(scaling);
    push_u32(&mut payload, scaling.len() as u32);
    payload.extend_from_slice(&scaling);
    push_space_signature(&mut payload, domain, codomain);
    push_structure_property(&mut payload, property);
    push_structure_support(&mut payload, support);
    payload.push(disposition.tag());
    push_u64(&mut payload, canonical_f64_bits(tolerance));
    payload.extend_from_slice(norm.as_bytes());
    spectral_proposition_receipt(SpectralPropositionKindV1::Structure, &payload)
}

/// Canonical receipt for an exact regularity proposition.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the fully bound proposition exceeds the
/// canonical identity limits or cannot be encoded.
#[allow(clippy::too_many_arguments)] // Every argument is an independent part of the exact proposition identity.
#[must_use = "the regularity-proposition receipt or canonicalization failure must be handled"]
pub fn regularity_proposition_receipt(
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
    regularity: RegularityClassV1,
    disposition: WitnessDispositionV1,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    let mut payload = Vec::with_capacity(192);
    payload.extend_from_slice(subject.as_bytes());
    payload.push(scalar_field.tag() as u8);
    let class = class_bytes(class);
    push_u32(&mut payload, class.len() as u32);
    payload.extend_from_slice(&class);
    let scaling = scaling_bytes(scaling);
    push_u32(&mut payload, scaling.len() as u32);
    payload.extend_from_slice(&scaling);
    push_space_signature(&mut payload, domain, codomain);
    payload.push(regularity.tag());
    if let RegularityClassV1::RegularPolynomial { grade }
    | RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade } = regularity
    {
        push_u32(&mut payload, grade);
    }
    payload.push(disposition.tag());
    spectral_proposition_receipt(SpectralPropositionKindV1::Regularity, &payload)
}

/// Canonical receipt for exact metric bounds, signature, or rank.
///
/// # Errors
///
/// Returns [`CanonicalError`] when positive-definite bounds are non-finite or
/// the proposition exceeds the canonical identity limits.
#[must_use = "the metric-proposition receipt or canonicalization failure must be handled"]
pub fn metric_proposition_receipt(
    metric: SpectralMetricId,
    dimension: u32,
    proposition: MetricDefinitenessPropositionV1,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    if let MetricDefinitenessPropositionV1::PositiveDefinite { lower, upper } = proposition {
        if !lower.is_finite() {
            return Err(CanonicalError::NonFiniteFloat {
                bits: lower.to_bits(),
            });
        }
        if !upper.is_finite() {
            return Err(CanonicalError::NonFiniteFloat {
                bits: upper.to_bits(),
            });
        }
    }
    let mut payload = Vec::with_capacity(96);
    payload.extend_from_slice(metric.as_bytes());
    push_u32(&mut payload, dimension);
    payload.push(proposition.tag());
    match proposition {
        MetricDefinitenessPropositionV1::PositiveDefinite { lower, upper } => {
            push_u64(&mut payload, canonical_f64_bits(lower));
            push_u64(&mut payload, canonical_f64_bits(upper));
        }
        MetricDefinitenessPropositionV1::Indefinite { positive, negative } => {
            push_u32(&mut payload, positive);
            push_u32(&mut payload, negative);
        }
        MetricDefinitenessPropositionV1::Singular { rank } => push_u32(&mut payload, rank),
    }
    spectral_proposition_receipt(SpectralPropositionKindV1::MetricDefiniteness, &payload)
}

fn nullspace_proposition_payload(
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
) -> Vec<u8> {
    let mut payload = Vec::with_capacity(144);
    payload.extend_from_slice(subject.as_bytes());
    payload.push(scalar_field.tag() as u8);
    let class = class_bytes(class);
    push_u32(&mut payload, class.len() as u32);
    payload.extend_from_slice(&class);
    let scaling = scaling_bytes(scaling);
    push_u32(&mut payload, scaling.len() as u32);
    payload.extend_from_slice(&scaling);
    push_space_signature(&mut payload, domain, codomain);
    payload
}

/// Canonical receipt for the exact gauge/nullspace convention.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the fully bound proposition exceeds the
/// canonical identity limits or cannot be encoded.
#[must_use = "the gauge-proposition receipt or canonicalization failure must be handled"]
pub fn gauge_proposition_receipt(
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
    proposition: GaugePropositionV1,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    let mut payload =
        nullspace_proposition_payload(subject, scalar_field, class, scaling, domain, codomain);
    payload.push(proposition.tag());
    match proposition {
        GaugePropositionV1::Fixed { nullity, gauge } => {
            push_u32(&mut payload, nullity);
            payload.extend_from_slice(gauge.as_bytes());
        }
        GaugePropositionV1::Quotiented { nullity, quotient } => {
            push_u32(&mut payload, nullity);
            payload.extend_from_slice(quotient.as_bytes());
        }
        GaugePropositionV1::None => {}
    }
    spectral_proposition_receipt(SpectralPropositionKindV1::Gauge, &payload)
}

/// Canonical receipt for the exact zero-padding serialization convention.
///
/// # Errors
///
/// Returns [`CanonicalError`] when the fully bound proposition exceeds the
/// canonical identity limits or cannot be encoded.
#[must_use = "the zero-padding receipt or canonicalization failure must be handled"]
pub fn zero_padding_proposition_receipt(
    subject: SpectralSubjectId,
    scalar_field: SpectralScalarFieldV1,
    class: SpectralProblemClassV1,
    scaling: SpectralScalingContextV1,
    domain: SpectralMetricV1,
    codomain: SpectralMetricV1,
    gauge: GaugeContextV1,
    proposition: ZeroPaddingPropositionV1,
) -> Result<IdentityReceipt<SpectralPropositionId>, CanonicalError> {
    let mut payload =
        nullspace_proposition_payload(subject, scalar_field, class, scaling, domain, codomain);
    push_gauge_context(&mut payload, gauge);
    payload.push(proposition.tag());
    match proposition {
        ZeroPaddingPropositionV1::ExplicitlyPadded { count }
        | ZeroPaddingPropositionV1::Omitted { count } => push_u32(&mut payload, count),
        ZeroPaddingPropositionV1::NonePresent => {}
    }
    spectral_proposition_receipt(SpectralPropositionKindV1::ZeroPadding, &payload)
}

fn claim_bytes(claim: &StructureClaimV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    push_structure_property(&mut out, claim.property);
    push_structure_support(&mut out, claim.support);
    out.push(claim.disposition.tag());
    push_u64(&mut out, canonical_f64_bits(claim.tolerance));
    out.extend_from_slice(claim.norm.as_bytes());
    out.extend_from_slice(&witness_bytes(&claim.witness));
    out
}

fn scaling_bytes(scaling: SpectralScalingContextV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(32 * 5 + 32);
    out.extend_from_slice(scaling.id.as_bytes());
    for exponent in scaling.spectral_dims.0 {
        out.push(exponent as u8);
    }
    push_u64(&mut out, canonical_f64_bits(scaling.spectral_scale_si));
    out.extend_from_slice(scaling.left_map.as_bytes());
    out.extend_from_slice(scaling.right_map.as_bytes());
    out.extend_from_slice(scaling.operator_map.as_bytes());
    out.extend_from_slice(scaling.inverse_map.as_bytes());
    out
}

fn metric_bytes(metric: SpectralMetricV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    out.extend_from_slice(metric.id.as_bytes());
    push_u32(&mut out, metric.dimension);
    out.push(metric.definiteness.tag());
    match metric.definiteness {
        MetricDefinitenessV1::Euclidean | MetricDefinitenessV1::Unknown => {}
        MetricDefinitenessV1::PositiveDefinite {
            lower,
            upper,
            witness,
        } => {
            push_u64(&mut out, canonical_f64_bits(lower));
            push_u64(&mut out, canonical_f64_bits(upper));
            out.extend_from_slice(&witness_bytes(&witness));
        }
        MetricDefinitenessV1::Indefinite {
            positive,
            negative,
            witness,
        } => {
            push_u32(&mut out, positive);
            push_u32(&mut out, negative);
            out.extend_from_slice(&witness_bytes(&witness));
        }
        MetricDefinitenessV1::Singular { rank, witness } => {
            push_u32(&mut out, rank);
            out.extend_from_slice(&witness_bytes(&witness));
        }
    }
    out
}

fn spaces_bytes(spaces: SpectralSpaceContextV1) -> Vec<u8> {
    let domain = metric_bytes(spaces.domain);
    let codomain = metric_bytes(spaces.codomain);
    let mut out = Vec::with_capacity(domain.len() + codomain.len() + 80);
    push_u32(&mut out, domain.len() as u32);
    out.extend_from_slice(&domain);
    push_u32(&mut out, codomain.len() as u32);
    out.extend_from_slice(&codomain);
    out.push(spaces.gauge.tag());
    match spaces.gauge {
        GaugeConventionV1::CertifiedNone { witness } => {
            out.extend_from_slice(&witness_bytes(&witness));
        }
        GaugeConventionV1::Unknown => {}
        GaugeConventionV1::Fixed {
            nullity,
            gauge,
            witness,
        } => {
            push_u32(&mut out, nullity);
            out.extend_from_slice(gauge.as_bytes());
            out.extend_from_slice(&witness_bytes(&witness));
        }
        GaugeConventionV1::Quotiented {
            nullity,
            quotient,
            witness,
        } => {
            push_u32(&mut out, nullity);
            out.extend_from_slice(quotient.as_bytes());
            out.extend_from_slice(&witness_bytes(&witness));
        }
    }
    out.push(spaces.zero_padding.tag());
    match spaces.zero_padding {
        ZeroPaddingConventionV1::CertifiedNonePresent { witness } => {
            out.extend_from_slice(&witness_bytes(&witness));
        }
        ZeroPaddingConventionV1::ExplicitlyPadded { count, witness }
        | ZeroPaddingConventionV1::Omitted { count, witness } => {
            push_u32(&mut out, count);
            out.extend_from_slice(&witness_bytes(&witness));
        }
        ZeroPaddingConventionV1::Unknown => {}
    }
    out
}

fn regularity_claim_bytes(claim: &RegularityClaimV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    out.push(claim.class.tag());
    if let RegularityClassV1::RegularPolynomial { grade }
    | RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade } = claim.class
    {
        push_u32(&mut out, grade);
    }
    out.push(claim.disposition.tag());
    out.extend_from_slice(&witness_bytes(&claim.witness));
    out
}

fn ordering_payload(ordering: SpectralOrderingV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(72);
    match ordering {
        SpectralOrderingV1::MagnitudeAscending { tie_break } => {
            out.push(tie_break.tag());
        }
        SpectralOrderingV1::NearestShift {
            real,
            imag,
            tie_break,
        } => {
            push_u64(&mut out, canonical_f64_bits(real));
            push_u64(&mut out, canonical_f64_bits(imag));
            out.push(tie_break.tag());
        }
        SpectralOrderingV1::NamedRegion { region } => {
            out.extend_from_slice(region.as_bytes());
        }
        SpectralOrderingV1::Projective {
            chart,
            infinity,
            tie_break,
        } => {
            out.extend_from_slice(chart.as_bytes());
            out.push(infinity.tag());
            out.push(tie_break.tag());
        }
        SpectralOrderingV1::SetValued
        | SpectralOrderingV1::RealAscending
        | SpectralOrderingV1::RealDescending => {}
    }
    out
}

fn scope_payload(scope: CompletenessScopeV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(72);
    match scope {
        CompletenessScopeV1::CandidateOnly => {}
        CompletenessScopeV1::Partial { requested } => push_u32(&mut out, requested),
        CompletenessScopeV1::Region { region, boundary } => {
            out.extend_from_slice(region.as_bytes());
            out.push(boundary.tag());
        }
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality,
            infinity_policy,
        } => {
            push_u32(&mut out, algebraic_cardinality);
            out.push(infinity_policy.tag());
        }
    }
    out
}

fn canonical_problem_receipt(
    spec: &SpectralProblemSpecV1,
    claims: &[StructureClaimV1],
    regularity_claims: &[RegularityClaimV1],
) -> Result<SpectralProblemIdentityReceiptV2, CanonicalError> {
    let class = class_bytes(spec.class);
    let scaling = scaling_bytes(spec.scaling);
    let spaces = spaces_bytes(spec.spaces);
    let ordering = ordering_payload(spec.ordering);
    let scope = scope_payload(spec.requested_scope);
    let mut claim_payloads: Vec<Vec<u8>> = claims.iter().map(claim_bytes).collect();
    claim_payloads.sort();
    let mut regularity_payloads: Vec<Vec<u8>> = regularity_claims
        .iter()
        .map(regularity_claim_bytes)
        .collect();
    regularity_payloads.sort();

    CanonicalEncoder::<ProblemSemanticId<SpectralProblemIdentitySchemaV2>, _>::new(
        PROBLEM_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .bytes(Field::new(0, "subject"), spec.subject.as_bytes())?
    .variant(Field::new(1, "scalar-field"), spec.scalar_field.tag(), &[])?
    .bytes(Field::new(2, "class"), &class)?
    .canonical_set(
        Field::new(3, "structure-claims"),
        claim_payloads.len() as u64,
        claim_payloads.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(4, "scaling"), &scaling)?
    .bytes(Field::new(5, "spaces"), &spaces)?
    .canonical_set(
        Field::new(6, "regularity"),
        regularity_payloads.len() as u64,
        regularity_payloads.iter().map(Vec::as_slice),
    )?
    .variant(Field::new(7, "ordering"), spec.ordering.tag(), &ordering)?
    .variant(
        Field::new(8, "requested-scope"),
        spec.requested_scope.tag(),
        &scope,
    )?
    .finish()
}

fn validate_witness_receipt(
    witness: &AdmittedSpectralWitnessV1,
    expected: Result<IdentityReceipt<SpectralPropositionId>, CanonicalError>,
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) {
    let expected = match expected {
        Ok(expected) => expected,
        Err(error) => {
            issues.push(SpectralAdmissionIssueV1::Identity(error));
            return;
        }
    };
    if witness.proposition != expected.id() {
        issues.push(SpectralAdmissionIssueV1::WitnessPropositionMismatch {
            expected: expected.id(),
            found: witness.proposition,
        });
    } else if !witness.matches_receipt(expected) {
        issues.push(SpectralAdmissionIssueV1::WitnessObservationMismatch {
            proposition: expected.id(),
        });
    }
}

fn validate_metric(metric: SpectralMetricV1, issues: &mut Vec<SpectralAdmissionIssueV1>) {
    match metric.definiteness {
        MetricDefinitenessV1::Euclidean => {
            if metric != SpectralMetricV1::euclidean(metric.dimension) {
                issues.push(SpectralAdmissionIssueV1::InvalidMetric { metric: metric.id });
            }
        }
        MetricDefinitenessV1::Unknown => {}
        MetricDefinitenessV1::PositiveDefinite {
            lower,
            upper,
            witness,
        } => {
            if !(lower.is_finite() && upper.is_finite() && lower > 0.0 && upper >= lower) {
                issues.push(SpectralAdmissionIssueV1::InvalidMetric { metric: metric.id });
            }
            validate_witness_receipt(
                &witness,
                metric_proposition_receipt(
                    metric.id,
                    metric.dimension,
                    MetricDefinitenessPropositionV1::PositiveDefinite { lower, upper },
                ),
                issues,
            );
        }
        MetricDefinitenessV1::Indefinite {
            positive,
            negative,
            witness,
        } => {
            let total = positive.checked_add(negative);
            if positive == 0 || negative == 0 || total != Some(metric.dimension) {
                issues.push(SpectralAdmissionIssueV1::InvalidMetric { metric: metric.id });
            }
            validate_witness_receipt(
                &witness,
                metric_proposition_receipt(
                    metric.id,
                    metric.dimension,
                    MetricDefinitenessPropositionV1::Indefinite { positive, negative },
                ),
                issues,
            );
        }
        MetricDefinitenessV1::Singular { rank, witness } => {
            if rank >= metric.dimension {
                issues.push(SpectralAdmissionIssueV1::InvalidMetric { metric: metric.id });
            }
            validate_witness_receipt(
                &witness,
                metric_proposition_receipt(
                    metric.id,
                    metric.dimension,
                    MetricDefinitenessPropositionV1::Singular { rank },
                ),
                issues,
            );
        }
    }
}

fn regularity_matches(spec: &SpectralProblemSpecV1, class: RegularityClassV1) -> bool {
    match class {
        RegularityClassV1::FiniteDimensional => true,
        RegularityClassV1::RegularPencil | RegularityClassV1::InvertiblePencilWeight => matches!(
            spec.class.representation,
            SpectralRepresentationV1::GeneralizedPencil
        ),
        RegularityClassV1::RegularPolynomial { grade }
        | RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade } => matches!(
            spec.class.representation,
            SpectralRepresentationV1::MatrixPolynomial { grade: declared } if grade == declared
        ),
        RegularityClassV1::RegularDescriptor => {
            matches!(spec.class.descriptor, DescriptorRoleV1::Descriptor { .. })
        }
        RegularityClassV1::WellPosedMonodromy => matches!(
            spec.class.origin,
            SpectralOperatorOriginV1::MonodromyFloquet { .. }
        ),
        RegularityClassV1::AnalyticOperatorFunction => matches!(
            spec.class.origin,
            SpectralOperatorOriginV1::AnalyticOperatorFunction { .. }
        ),
    }
}

fn support_matches_property(property: StructurePropertyV1, support: StructureSupportV1) -> bool {
    match property {
        StructurePropertyV1::SelfAdjoint
        | StructurePropertyV1::Normal
        | StructurePropertyV1::Nonnormal
        | StructurePropertyV1::HermitianDefinitePencil => {
            matches!(support, StructureSupportV1::InnerProduct(_))
        }
        StructurePropertyV1::Hamiltonian | StructurePropertyV1::Symplectic => {
            matches!(support, StructureSupportV1::SymplecticForm(_))
        }
        StructurePropertyV1::JSelfAdjoint => {
            matches!(support, StructureSupportV1::KreinForm(_))
        }
        StructurePropertyV1::RealConjugatePairs => {
            matches!(support, StructureSupportV1::Conjugation(_))
        }
        StructurePropertyV1::Gyroscopic
        | StructurePropertyV1::Palindromic { .. }
        | StructurePropertyV1::RealSpectrum => {
            matches!(support, StructureSupportV1::FormFree)
        }
    }
}

fn support_is_positive_definite(spec: &SpectralProblemSpecV1, support: StructureSupportV1) -> bool {
    shared_inner_product(spec, support)
        .is_some_and(|metric| metric.definiteness.is_positive_definite())
}

fn support_is_adjoint_compatible(
    spec: &SpectralProblemSpecV1,
    support: StructureSupportV1,
) -> bool {
    shared_inner_product(spec, support)
        .is_some_and(|metric| metric.definiteness.is_adjoint_compatible())
}

fn shared_inner_product(
    spec: &SpectralProblemSpecV1,
    support: StructureSupportV1,
) -> Option<SpectralMetricV1> {
    let StructureSupportV1::InnerProduct(metric_id) = support else {
        return None;
    };
    let domain = spec.spaces.domain;
    let codomain = spec.spaces.codomain;
    (domain == codomain && domain.id == metric_id).then_some(domain)
}

fn property_matches_representation(
    property: StructurePropertyV1,
    representation: SpectralRepresentationV1,
) -> bool {
    match property {
        StructurePropertyV1::HermitianDefinitePencil => {
            matches!(representation, SpectralRepresentationV1::GeneralizedPencil)
        }
        StructurePropertyV1::Gyroscopic => matches!(
            representation,
            SpectralRepresentationV1::MatrixPolynomial { grade: 2 }
        ),
        StructurePropertyV1::Palindromic { .. } => {
            matches!(
                representation,
                SpectralRepresentationV1::MatrixPolynomial { .. }
            )
        }
        StructurePropertyV1::SelfAdjoint
        | StructurePropertyV1::Normal
        | StructurePropertyV1::Nonnormal
        | StructurePropertyV1::Hamiltonian
        | StructurePropertyV1::Symplectic
        | StructurePropertyV1::JSelfAdjoint
        | StructurePropertyV1::RealSpectrum
        | StructurePropertyV1::RealConjugatePairs => true,
    }
}

#[allow(clippy::too_many_lines)] // One local theorem registry keeps each premise and contradiction visibly adjacent.
fn validate_exact_structure_theorems(
    spec: &SpectralProblemSpecV1,
    claims: &[StructureClaimV1],
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) {
    let exact_witness = |claim: &StructureClaimV1| {
        claim.disposition == WitnessDispositionV1::Witnessed
            && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
    };
    let contradicts = |property: StructurePropertyV1, support: Option<StructureSupportV1>| {
        claims.iter().any(|claim| {
            claim.property == property
                && support.is_none_or(|expected| claim.support == expected)
                && match property {
                    StructurePropertyV1::Nonnormal => {
                        claim.disposition == WitnessDispositionV1::Witnessed
                    }
                    StructurePropertyV1::SelfAdjoint
                    | StructurePropertyV1::Normal
                    | StructurePropertyV1::RealSpectrum => {
                        claim.disposition == WitnessDispositionV1::Contradicted
                    }
                    StructurePropertyV1::Hamiltonian
                    | StructurePropertyV1::Symplectic
                    | StructurePropertyV1::JSelfAdjoint
                    | StructurePropertyV1::HermitianDefinitePencil
                    | StructurePropertyV1::Gyroscopic
                    | StructurePropertyV1::Palindromic { .. }
                    | StructurePropertyV1::RealConjugatePairs => false,
                }
        })
    };
    for normal_refutation in claims.iter().filter(|claim| {
        claim.property == StructurePropertyV1::Normal
            && claim.disposition == WitnessDispositionV1::Contradicted
            && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
    }) {
        let support = normal_refutation.support;
        let nonnormal_also_refuted = claims.iter().any(|claim| {
            claim.property == StructurePropertyV1::Nonnormal
                && claim.support == support
                && claim.disposition == WitnessDispositionV1::Contradicted
                && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
        });
        if nonnormal_also_refuted && support_is_adjoint_compatible(spec, support) {
            issues.push(SpectralAdmissionIssueV1::ComplementaryStructureConflict { support });
        }
    }
    for premise in claims.iter().filter(|claim| exact_witness(claim)) {
        let standard_linear = matches!(
            spec.class.representation,
            SpectralRepresentationV1::StandardLinear
        );
        let valid_adjoint_support = support_is_adjoint_compatible(spec, premise.support);
        if valid_adjoint_support
            && premise.property == StructurePropertyV1::Normal
            && contradicts(StructurePropertyV1::Nonnormal, Some(premise.support))
        {
            issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                premise: StructurePropertyV1::Normal,
                consequence: StructurePropertyV1::Normal,
                support: premise.support,
            });
        }
        if standard_linear
            && valid_adjoint_support
            && premise.property == StructurePropertyV1::SelfAdjoint
        {
            if contradicts(StructurePropertyV1::Normal, Some(premise.support))
                || contradicts(StructurePropertyV1::Nonnormal, Some(premise.support))
            {
                issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                    premise: StructurePropertyV1::SelfAdjoint,
                    consequence: StructurePropertyV1::Normal,
                    support: premise.support,
                });
            }
            if support_is_positive_definite(spec, premise.support)
                && contradicts(StructurePropertyV1::RealSpectrum, None)
            {
                issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                    premise: StructurePropertyV1::SelfAdjoint,
                    consequence: StructurePropertyV1::RealSpectrum,
                    support: premise.support,
                });
            }
        }
        if premise.property == StructurePropertyV1::HermitianDefinitePencil
            && matches!(
                spec.class.representation,
                SpectralRepresentationV1::GeneralizedPencil
            )
            && support_is_positive_definite(spec, premise.support)
        {
            if contradicts(StructurePropertyV1::SelfAdjoint, Some(premise.support)) {
                issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                    premise: StructurePropertyV1::HermitianDefinitePencil,
                    consequence: StructurePropertyV1::SelfAdjoint,
                    support: premise.support,
                });
            }
            if contradicts(StructurePropertyV1::Normal, Some(premise.support))
                || contradicts(StructurePropertyV1::Nonnormal, Some(premise.support))
            {
                issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                    premise: StructurePropertyV1::HermitianDefinitePencil,
                    consequence: StructurePropertyV1::Normal,
                    support: premise.support,
                });
            }
            if contradicts(StructurePropertyV1::RealSpectrum, None) {
                issues.push(SpectralAdmissionIssueV1::StructureTheoremConflict {
                    premise: StructurePropertyV1::HermitianDefinitePencil,
                    consequence: StructurePropertyV1::RealSpectrum,
                    support: premise.support,
                });
            }
        }
    }
}

fn validate_regularity_theorems(
    spec: &SpectralProblemSpecV1,
    structures: &[StructureClaimV1],
    regularity: &[RegularityClaimV1],
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) {
    let contradicted = |consequence: RegularityClassV1| {
        regularity.iter().any(|claim| {
            claim.class == consequence && claim.disposition == WitnessDispositionV1::Contradicted
        })
    };

    if matches!(
        spec.class.representation,
        SpectralRepresentationV1::GeneralizedPencil
    ) {
        for premise in structures.iter().filter(|claim| {
            claim.property == StructurePropertyV1::HermitianDefinitePencil
                && claim.disposition == WitnessDispositionV1::Witnessed
                && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
                && support_is_positive_definite(spec, claim.support)
        }) {
            for consequence in [
                RegularityClassV1::InvertiblePencilWeight,
                RegularityClassV1::RegularPencil,
            ] {
                if contradicted(consequence) {
                    issues.push(
                        SpectralAdmissionIssueV1::StructureRegularityTheoremConflict {
                            premise: StructurePropertyV1::HermitianDefinitePencil,
                            consequence,
                            support: premise.support,
                        },
                    );
                }
            }
        }
    }

    for premise in regularity
        .iter()
        .filter(|claim| claim.disposition == WitnessDispositionV1::Witnessed)
    {
        let consequence = match premise.class {
            RegularityClassV1::InvertiblePencilWeight => Some(RegularityClassV1::RegularPencil),
            RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade } => {
                Some(RegularityClassV1::RegularPolynomial { grade })
            }
            RegularityClassV1::FiniteDimensional
            | RegularityClassV1::RegularPencil
            | RegularityClassV1::RegularPolynomial { .. }
            | RegularityClassV1::RegularDescriptor
            | RegularityClassV1::WellPosedMonodromy
            | RegularityClassV1::AnalyticOperatorFunction => None,
        };
        if let Some(consequence) = consequence
            && contradicted(consequence)
        {
            issues.push(SpectralAdmissionIssueV1::RegularityTheoremConflict {
                premise: premise.class,
                consequence,
            });
        }
    }
}

fn validate_ordinary_finite_spectrum(
    spec: &SpectralProblemSpecV1,
    claims: &[StructureClaimV1],
    regularity: &[RegularityClaimV1],
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) {
    if !matches!(spec.class.descriptor, DescriptorRoleV1::Ordinary) {
        return;
    }
    let required = match spec.class.representation {
        SpectralRepresentationV1::StandardLinear => return,
        SpectralRepresentationV1::GeneralizedPencil => RegularityClassV1::InvertiblePencilWeight,
        SpectralRepresentationV1::MatrixPolynomial { grade } => {
            RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade }
        }
    };
    let contradicted = regularity.iter().any(|claim| {
        claim.class == required && claim.disposition == WitnessDispositionV1::Contradicted
    });
    let witnessed = regularity.iter().any(|claim| {
        claim.class == required && claim.disposition == WitnessDispositionV1::Witnessed
    });
    let theorem_backed_pencil = matches!(required, RegularityClassV1::InvertiblePencilWeight)
        && claims.iter().any(|claim| {
            claim.property == StructurePropertyV1::HermitianDefinitePencil
                && claim.disposition == WitnessDispositionV1::Witnessed
                && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
                && support_is_positive_definite(spec, claim.support)
        });
    if contradicted || !(witnessed || theorem_backed_pencil) {
        issues.push(SpectralAdmissionIssueV1::OrdinaryFiniteSpectrumWitnessRequired { required });
    }
}

fn equation_algebraic_cardinality(spec: &SpectralProblemSpecV1) -> Option<u32> {
    let grade = match spec.class.representation {
        SpectralRepresentationV1::StandardLinear | SpectralRepresentationV1::GeneralizedPencil => 1,
        SpectralRepresentationV1::MatrixPolynomial { grade } => grade,
    };
    spec.spaces.domain.dimension.checked_mul(grade)
}

fn regularity_is_admitted(
    spec: &SpectralProblemSpecV1,
    structures: &[StructureClaimV1],
    regularity: &[RegularityClaimV1],
    required: RegularityClassV1,
) -> bool {
    if regularity.iter().any(|claim| {
        claim.class == required && claim.disposition == WitnessDispositionV1::Witnessed
    }) {
        return true;
    }

    match required {
        RegularityClassV1::InvertiblePencilWeight
            if matches!(
                spec.class.representation,
                SpectralRepresentationV1::GeneralizedPencil
            ) =>
        {
            structures.iter().any(|claim| {
                claim.property == StructurePropertyV1::HermitianDefinitePencil
                    && claim.disposition == WitnessDispositionV1::Witnessed
                    && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
                    && support_is_positive_definite(spec, claim.support)
            })
        }
        RegularityClassV1::RegularPencil
            if matches!(
                spec.class.representation,
                SpectralRepresentationV1::GeneralizedPencil
            ) =>
        {
            regularity_is_admitted(
                spec,
                structures,
                regularity,
                RegularityClassV1::InvertiblePencilWeight,
            )
        }
        RegularityClassV1::RegularPolynomial { grade }
            if matches!(
                spec.class.representation,
                SpectralRepresentationV1::MatrixPolynomial {
                    grade: representation_grade
                } if representation_grade == grade
            ) =>
        {
            regularity.iter().any(|claim| {
                claim.class == RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade }
                    && claim.disposition == WitnessDispositionV1::Witnessed
            })
        }
        RegularityClassV1::FiniteDimensional
        | RegularityClassV1::RegularPencil
        | RegularityClassV1::InvertiblePencilWeight
        | RegularityClassV1::RegularPolynomial { .. }
        | RegularityClassV1::InvertiblePolynomialLeadingCoefficient { .. }
        | RegularityClassV1::RegularDescriptor
        | RegularityClassV1::WellPosedMonodromy
        | RegularityClassV1::AnalyticOperatorFunction => false,
    }
}

fn real_spectrum_is_admitted(
    spec: &SpectralProblemSpecV1,
    structures: &[StructureClaimV1],
) -> bool {
    structures.iter().any(|claim| {
        claim.disposition == WitnessDispositionV1::Witnessed
            && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
            && match claim.property {
                StructurePropertyV1::RealSpectrum => claim.support == StructureSupportV1::FormFree,
                StructurePropertyV1::SelfAdjoint => {
                    matches!(
                        spec.class.representation,
                        SpectralRepresentationV1::StandardLinear
                    ) && support_is_positive_definite(spec, claim.support)
                }
                StructurePropertyV1::HermitianDefinitePencil => {
                    matches!(
                        spec.class.representation,
                        SpectralRepresentationV1::GeneralizedPencil
                    ) && support_is_positive_definite(spec, claim.support)
                }
                StructurePropertyV1::Normal
                | StructurePropertyV1::Nonnormal
                | StructurePropertyV1::Hamiltonian
                | StructurePropertyV1::Symplectic
                | StructurePropertyV1::JSelfAdjoint
                | StructurePropertyV1::Gyroscopic
                | StructurePropertyV1::Palindromic { .. }
                | StructurePropertyV1::RealConjugatePairs => false,
            }
    })
}

fn projective_infinity_is_excluded_by_theorem(
    spec: &SpectralProblemSpecV1,
    structures: &[StructureClaimV1],
    regularity: &[RegularityClaimV1],
) -> bool {
    match spec.class.representation {
        SpectralRepresentationV1::StandardLinear => true,
        SpectralRepresentationV1::GeneralizedPencil => regularity_is_admitted(
            spec,
            structures,
            regularity,
            RegularityClassV1::InvertiblePencilWeight,
        ),
        SpectralRepresentationV1::MatrixPolynomial { grade } => regularity_is_admitted(
            spec,
            structures,
            regularity,
            RegularityClassV1::InvertiblePolynomialLeadingCoefficient { grade },
        ),
    }
}

fn admitted_algebraic_cardinality(
    spec: &SpectralProblemSpecV1,
    structures: &[StructureClaimV1],
    regularity: &[RegularityClaimV1],
) -> Option<u32> {
    let equation_regular = match spec.class.representation {
        SpectralRepresentationV1::StandardLinear => true,
        SpectralRepresentationV1::GeneralizedPencil => regularity_is_admitted(
            spec,
            structures,
            regularity,
            RegularityClassV1::RegularPencil,
        ),
        SpectralRepresentationV1::MatrixPolynomial { grade } => regularity_is_admitted(
            spec,
            structures,
            regularity,
            RegularityClassV1::RegularPolynomial { grade },
        ),
    };
    let class_regular = match spec.class.descriptor {
        DescriptorRoleV1::Ordinary => true,
        DescriptorRoleV1::Descriptor { .. } => {
            equation_regular
                && regularity.iter().any(|claim| {
                    claim.class == RegularityClassV1::RegularDescriptor
                        && claim.disposition == WitnessDispositionV1::Witnessed
                })
        }
    };
    if class_regular {
        equation_algebraic_cardinality(spec)
    } else {
        None
    }
}

/// Validate raw input and mint its typed semantic identity. All detected
/// defects are returned in deterministic rank order; no partial admitted token
/// escapes on failure.
///
/// # Errors
///
/// Returns [`SpectralAdmissionReportV1`] when any schema, unit, dimension,
/// metric, witness, structure, regularity, ordering, scope, or canonical
/// identity obligation fails.
#[allow(clippy::too_many_lines)] // One ordered fail-closed admission matrix keeps cross-field obligations auditable.
#[must_use = "the admission result must be handled before the problem can be used"]
pub fn validate_problem(
    mut spec: SpectralProblemSpecV1,
) -> Result<ValidatedSpectralProblemV1, SpectralAdmissionReportV1> {
    let mut issues = Vec::new();
    if spec.structures.claims.len() > MAX_STRUCTURE_CLAIMS_V1 {
        issues.push(SpectralAdmissionIssueV1::TooManyClaims {
            profile: ClaimProfileV1::Structure,
            found: spec.structures.claims.len(),
            limit: MAX_STRUCTURE_CLAIMS_V1,
        });
    }
    if spec.regularity.claims.len() > MAX_REGULARITY_CLAIMS_V1 {
        issues.push(SpectralAdmissionIssueV1::TooManyClaims {
            profile: ClaimProfileV1::Regularity,
            found: spec.regularity.claims.len(),
            limit: MAX_REGULARITY_CLAIMS_V1,
        });
    }
    if !issues.is_empty() {
        return Err(SpectralAdmissionReportV1::new(issues));
    }
    if spec.schema_version != SPECTRAL_PROBLEM_SCHEMA_VERSION {
        issues.push(SpectralAdmissionIssueV1::UnsupportedSchemaVersion {
            found: spec.schema_version,
            supported: SPECTRAL_PROBLEM_SCHEMA_VERSION,
        });
    }
    for claim in &mut spec.structures.claims {
        if claim.tolerance == 0.0 {
            claim.tolerance = 0.0;
        }
    }
    if let SpectralOrderingV1::NearestShift { real, imag, .. } = &mut spec.ordering {
        if *real == 0.0 {
            *real = 0.0;
        }
        if *imag == 0.0 {
            *imag = 0.0;
        }
    }
    if let SpectralOperatorOriginV1::MonodromyFloquet {
        branch: FloquetBranchConventionV1::ContinuousFrom { anchor_phase, .. },
        ..
    } = &mut spec.class.origin
        && anchor_phase.value() == 0.0
    {
        *anchor_phase = Angle::new(0.0);
    }
    if let SpectralRepresentationV1::MatrixPolynomial { grade } = spec.class.representation
        && grade == 0
    {
        issues.push(SpectralAdmissionIssueV1::Zero {
            field: AdmissionFieldV1::PolynomialGrade,
        });
    }
    if matches!(spec.class.descriptor, DescriptorRoleV1::Descriptor { .. })
        && matches!(
            spec.class.representation,
            SpectralRepresentationV1::StandardLinear
        )
    {
        issues.push(SpectralAdmissionIssueV1::RepresentationConflict);
    }
    if !(spec.scaling.spectral_scale_si.is_finite() && spec.scaling.spectral_scale_si > 0.0) {
        issues.push(SpectralAdmissionIssueV1::NonPositive {
            field: AdmissionFieldV1::SpectralScale,
        });
    }
    if spec.spaces.domain.dimension == 0 {
        issues.push(SpectralAdmissionIssueV1::Zero {
            field: AdmissionFieldV1::DomainDimension,
        });
    }
    if spec.spaces.codomain.dimension == 0 {
        issues.push(SpectralAdmissionIssueV1::Zero {
            field: AdmissionFieldV1::CodomainDimension,
        });
    }
    if spec.spaces.domain.dimension != spec.spaces.codomain.dimension {
        issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
            left: spec.spaces.domain.dimension,
            right: spec.spaces.codomain.dimension,
        });
    }
    if equation_algebraic_cardinality(&spec).is_none() {
        let grade = match spec.class.representation {
            SpectralRepresentationV1::MatrixPolynomial { grade } => grade,
            SpectralRepresentationV1::StandardLinear
            | SpectralRepresentationV1::GeneralizedPencil => 1,
        };
        issues.push(SpectralAdmissionIssueV1::AlgebraicCardinalityOverflow {
            dimension: spec.spaces.domain.dimension,
            grade,
        });
    }
    if spec.spaces.domain.id == spec.spaces.codomain.id
        && spec.spaces.domain != spec.spaces.codomain
    {
        issues.push(SpectralAdmissionIssueV1::MetricIdentityConflict {
            metric: spec.spaces.domain.id,
        });
    }
    validate_metric(spec.spaces.domain, &mut issues);
    validate_metric(spec.spaces.codomain, &mut issues);
    match spec.spaces.gauge {
        GaugeConventionV1::CertifiedNone { witness } => validate_witness_receipt(
            &witness,
            gauge_proposition_receipt(
                spec.subject,
                spec.scalar_field,
                spec.class,
                spec.scaling,
                spec.spaces.domain,
                spec.spaces.codomain,
                GaugePropositionV1::None,
            ),
            &mut issues,
        ),
        GaugeConventionV1::Fixed {
            nullity,
            gauge,
            witness,
        } => {
            if nullity == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::GaugeNullity,
                });
            } else if nullity > spec.spaces.domain.dimension {
                issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
                    left: nullity,
                    right: spec.spaces.domain.dimension,
                });
            }
            validate_witness_receipt(
                &witness,
                gauge_proposition_receipt(
                    spec.subject,
                    spec.scalar_field,
                    spec.class,
                    spec.scaling,
                    spec.spaces.domain,
                    spec.spaces.codomain,
                    GaugePropositionV1::Fixed { nullity, gauge },
                ),
                &mut issues,
            );
        }
        GaugeConventionV1::Quotiented {
            nullity,
            quotient,
            witness,
        } => {
            if nullity == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::GaugeNullity,
                });
            }
            validate_witness_receipt(
                &witness,
                gauge_proposition_receipt(
                    spec.subject,
                    spec.scalar_field,
                    spec.class,
                    spec.scaling,
                    spec.spaces.domain,
                    spec.spaces.codomain,
                    GaugePropositionV1::Quotiented { nullity, quotient },
                ),
                &mut issues,
            );
        }
        GaugeConventionV1::Unknown => {}
    }
    match spec.spaces.zero_padding {
        ZeroPaddingConventionV1::CertifiedNonePresent { witness } => validate_witness_receipt(
            &witness,
            zero_padding_proposition_receipt(
                spec.subject,
                spec.scalar_field,
                spec.class,
                spec.scaling,
                spec.spaces.domain,
                spec.spaces.codomain,
                spec.spaces.gauge.context(),
                ZeroPaddingPropositionV1::NonePresent,
            ),
            &mut issues,
        ),
        ZeroPaddingConventionV1::ExplicitlyPadded { count, witness } => {
            if count == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::ZeroPaddingCount,
                });
            } else if !matches!(spec.spaces.gauge, GaugeConventionV1::Quotiented { .. })
                && count > spec.spaces.domain.dimension
            {
                issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
                    left: count,
                    right: spec.spaces.domain.dimension,
                });
            }
            validate_witness_receipt(
                &witness,
                zero_padding_proposition_receipt(
                    spec.subject,
                    spec.scalar_field,
                    spec.class,
                    spec.scaling,
                    spec.spaces.domain,
                    spec.spaces.codomain,
                    spec.spaces.gauge.context(),
                    ZeroPaddingPropositionV1::ExplicitlyPadded { count },
                ),
                &mut issues,
            );
        }
        ZeroPaddingConventionV1::Omitted { count, witness } => {
            if count == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::ZeroPaddingCount,
                });
            } else if !matches!(spec.spaces.gauge, GaugeConventionV1::Quotiented { .. })
                && count > spec.spaces.domain.dimension
            {
                issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
                    left: count,
                    right: spec.spaces.domain.dimension,
                });
            }
            validate_witness_receipt(
                &witness,
                zero_padding_proposition_receipt(
                    spec.subject,
                    spec.scalar_field,
                    spec.class,
                    spec.scaling,
                    spec.spaces.domain,
                    spec.spaces.codomain,
                    spec.spaces.gauge.context(),
                    ZeroPaddingPropositionV1::Omitted { count },
                ),
                &mut issues,
            );
        }
        ZeroPaddingConventionV1::Unknown => {}
    }

    match spec.class.origin {
        SpectralOperatorOriginV1::Direct
        // `NoClaim` is representable for analytic-function classification,
        // but the specialized method obligation below refuses it.
        | SpectralOperatorOriginV1::AnalyticOperatorFunction { .. } => {}
        SpectralOperatorOriginV1::MonodromyFloquet {
            period,
            parameter,
            branch,
        } => {
            if !(period.value().is_finite() && period.value() > 0.0) {
                issues.push(SpectralAdmissionIssueV1::NonPositive {
                    field: AdmissionFieldV1::FloquetPeriod,
                });
            }
            let expected_dims = match parameter {
                FloquetParameterV1::Multiplier => Dims::NONE,
                FloquetParameterV1::Exponent => Dims([0, 0, -1, 0, 0, 0]),
            };
            let branch_ok = match (parameter, branch) {
                (FloquetParameterV1::Multiplier, FloquetBranchConventionV1::MultipliersOnly)
                | (FloquetParameterV1::Exponent, FloquetBranchConventionV1::PrincipalLog) => true,
                (
                    FloquetParameterV1::Exponent,
                    FloquetBranchConventionV1::ContinuousFrom { anchor_phase, .. },
                ) => anchor_phase.value().is_finite(),
                _ => false,
            };
            if spec.scaling.spectral_dims != expected_dims || !branch_ok {
                issues.push(SpectralAdmissionIssueV1::FloquetSemanticMismatch);
            }
        }
    }

    if let SpectralOrderingV1::NearestShift { real, imag, .. } = spec.ordering
        && !(real.is_finite() && imag.is_finite())
    {
        issues.push(SpectralAdmissionIssueV1::NonFinite {
            field: AdmissionFieldV1::OrderingShift,
        });
    }
    if matches!(
        spec.ordering,
        SpectralOrderingV1::RealAscending | SpectralOrderingV1::RealDescending
    ) {
        let real_spectrum_admitted = real_spectrum_is_admitted(&spec, &spec.structures.claims);
        let real_spectrum_contradicted = spec.structures.claims.iter().any(|claim| {
            claim.disposition == WitnessDispositionV1::Contradicted
                && claim.property == StructurePropertyV1::RealSpectrum
                && claim.support == StructureSupportV1::FormFree
        });
        let includes_projective_infinity = matches!(
            spec.class.descriptor,
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective
            }
        );
        let projective_infinity_excluded = projective_infinity_is_excluded_by_theorem(
            &spec,
            &spec.structures.claims,
            &spec.regularity.claims,
        );
        if !real_spectrum_admitted
            || real_spectrum_contradicted
            || (includes_projective_infinity && !projective_infinity_excluded)
        {
            issues.push(SpectralAdmissionIssueV1::OrderingUnavailable);
        }
    }
    match (spec.ordering, spec.requested_scope) {
        (
            SpectralOrderingV1::NamedRegion { region: ordered },
            CompletenessScopeV1::Region {
                region: requested, ..
            },
        ) if ordered == requested => {}
        (SpectralOrderingV1::NamedRegion { .. }, _) | (_, CompletenessScopeV1::Region { .. })
            if !matches!(spec.ordering, SpectralOrderingV1::SetValued) =>
        {
            issues.push(SpectralAdmissionIssueV1::ScopeOrderingMismatch);
        }
        _ => {}
    }
    if matches!(
        (spec.ordering, spec.requested_scope),
        (
            SpectralOrderingV1::SetValued,
            CompletenessScopeV1::Partial { .. }
        )
    ) {
        issues.push(SpectralAdmissionIssueV1::ScopeOrderingMismatch);
    }
    if matches!(spec.ordering, SpectralOrderingV1::Projective { .. })
        && !matches!(
            spec.class.descriptor,
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective
            }
        )
    {
        issues.push(SpectralAdmissionIssueV1::InfinityPolicyMismatch);
    }
    if matches!(spec.requested_scope, CompletenessScopeV1::Partial { .. })
        && matches!(
            spec.class.descriptor,
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::IncludeProjective
            }
        )
        && !matches!(spec.ordering, SpectralOrderingV1::Projective { .. })
        && !projective_infinity_is_excluded_by_theorem(
            &spec,
            &spec.structures.claims,
            &spec.regularity.claims,
        )
    {
        issues.push(SpectralAdmissionIssueV1::ProjectivePrefixOrderingRequired);
    }
    if matches!(spec.requested_scope, CompletenessScopeV1::Partial { .. })
        && matches!(
            spec.class.descriptor,
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim
                    | InfiniteEigenvaluePolicyV1::ExcludeWithCount
            }
        )
    {
        // V1 `Partial` is a prefix of the mathematical spectrum, not a
        // finite-only serialized view. Excluding or declining to classify
        // projective infinity would make that prefix unauditable; a future
        // finite-only scope must carry its own infinity accounting.
        issues.push(SpectralAdmissionIssueV1::InfinityPolicyMismatch);
    }
    if let (
        DescriptorRoleV1::Descriptor {
            infinity_policy: declared,
        },
        CompletenessScopeV1::FullFinite {
            infinity_policy: requested,
            ..
        },
    ) = (spec.class.descriptor, spec.requested_scope)
        && declared != requested
    {
        issues.push(SpectralAdmissionIssueV1::InfinityPolicyMismatch);
    }
    if let CompletenessScopeV1::FullFinite {
        infinity_policy: requested,
        ..
    } = spec.requested_scope
    {
        match spec.class.descriptor {
            DescriptorRoleV1::Ordinary if requested != InfiniteEigenvaluePolicyV1::NoClaim => {
                issues.push(SpectralAdmissionIssueV1::InfinityPolicyMismatch);
            }
            DescriptorRoleV1::Descriptor {
                infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim,
            } => issues.push(SpectralAdmissionIssueV1::InfinityPolicyMismatch),
            DescriptorRoleV1::Ordinary | DescriptorRoleV1::Descriptor { .. } => {}
        }
    }
    match spec.requested_scope {
        CompletenessScopeV1::Partial { requested } => {
            if requested == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::CompletenessScope,
                });
            }
        }
        CompletenessScopeV1::FullFinite {
            algebraic_cardinality,
            ..
        } => {
            if algebraic_cardinality == 0 {
                issues.push(SpectralAdmissionIssueV1::Zero {
                    field: AdmissionFieldV1::CompletenessScope,
                });
            } else {
                let expected = equation_algebraic_cardinality(&spec);
                match expected {
                    Some(expected) if expected != algebraic_cardinality => {
                        issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
                            left: algebraic_cardinality,
                            right: expected,
                        });
                    }
                    None => {
                        let grade = match spec.class.representation {
                            SpectralRepresentationV1::MatrixPolynomial { grade } => grade,
                            SpectralRepresentationV1::StandardLinear
                            | SpectralRepresentationV1::GeneralizedPencil => 1,
                        };
                        issues.push(SpectralAdmissionIssueV1::AlgebraicCardinalityOverflow {
                            dimension: spec.spaces.domain.dimension,
                            grade,
                        });
                    }
                    Some(_) => {}
                }
            }
        }
        CompletenessScopeV1::CandidateOnly | CompletenessScopeV1::Region { .. } => {}
    }

    let mut regularity_claims = spec.regularity.claims.clone();
    regularity_claims.sort_by_cached_key(regularity_claim_bytes);
    for claim in &regularity_claims {
        validate_witness_receipt(
            &claim.witness,
            regularity_proposition_receipt(
                spec.subject,
                spec.scalar_field,
                spec.class,
                spec.scaling,
                spec.spaces.domain,
                spec.spaces.codomain,
                claim.class,
                claim.disposition,
            ),
            &mut issues,
        );
        if !regularity_matches(&spec, claim.class)
            || (claim.class == RegularityClassV1::FiniteDimensional
                && claim.disposition == WitnessDispositionV1::Contradicted)
        {
            issues.push(SpectralAdmissionIssueV1::RegularityMismatch);
        }
    }
    for pair in regularity_claims.windows(2) {
        if pair[0].class == pair[1].class {
            issues.push(SpectralAdmissionIssueV1::RegularityMismatch);
        }
    }

    let mut claims = spec.structures.claims.clone();
    claims.sort_by_key(|claim| {
        (
            claim.property,
            claim.support,
            claim.norm,
            canonical_f64_bits(claim.tolerance),
            claim.disposition,
            claim.witness.proposition,
        )
    });
    for claim in &claims {
        validate_witness_receipt(
            &claim.witness,
            structure_proposition_receipt(
                spec.subject,
                spec.scalar_field,
                spec.class,
                spec.scaling,
                spec.spaces.domain,
                spec.spaces.codomain,
                claim.property,
                claim.support,
                claim.disposition,
                claim.tolerance,
                claim.norm,
            ),
            &mut issues,
        );
        let support_known = match claim.support {
            StructureSupportV1::InnerProduct(metric) => {
                metric == spec.spaces.domain.id || metric == spec.spaces.codomain.id
            }
            StructureSupportV1::SymplecticForm(_)
            | StructureSupportV1::KreinForm(_)
            | StructureSupportV1::Conjugation(_)
            | StructureSupportV1::FormFree => true,
        };
        if !(claim.tolerance.is_finite() && claim.tolerance >= 0.0) {
            issues.push(SpectralAdmissionIssueV1::InvalidStructureTolerance {
                property: claim.property,
            });
        }
        let definiteness_matches = match claim.property {
            StructurePropertyV1::HermitianDefinitePencil => {
                support_is_positive_definite(&spec, claim.support)
            }
            StructurePropertyV1::SelfAdjoint
            | StructurePropertyV1::Normal
            | StructurePropertyV1::Nonnormal => support_is_adjoint_compatible(&spec, claim.support),
            StructurePropertyV1::Hamiltonian
            | StructurePropertyV1::Symplectic
            | StructurePropertyV1::JSelfAdjoint
            | StructurePropertyV1::Gyroscopic
            | StructurePropertyV1::Palindromic { .. }
            | StructurePropertyV1::RealSpectrum
            | StructurePropertyV1::RealConjugatePairs => true,
        };
        if !support_known
            || !support_matches_property(claim.property, claim.support)
            || !definiteness_matches
        {
            issues.push(SpectralAdmissionIssueV1::InvalidStructureSupport {
                property: claim.property,
                support: claim.support,
            });
        }
        if !property_matches_representation(claim.property, spec.class.representation) {
            issues.push(SpectralAdmissionIssueV1::StructureRepresentationMismatch {
                property: claim.property,
                representation: spec.class.representation,
            });
        }
    }
    for pair in claims.windows(2) {
        let left_key = (
            pair[0].property,
            pair[0].support,
            pair[0].norm,
            canonical_f64_bits(pair[0].tolerance),
        );
        let right_key = (
            pair[1].property,
            pair[1].support,
            pair[1].norm,
            canonical_f64_bits(pair[1].tolerance),
        );
        if left_key == right_key {
            if pair[0].disposition == pair[1].disposition {
                issues.push(SpectralAdmissionIssueV1::DuplicateStructure {
                    property: pair[0].property,
                    support: pair[0].support,
                });
            } else {
                issues.push(SpectralAdmissionIssueV1::ContradictoryStructure {
                    property: pair[0].property,
                    support: pair[0].support,
                });
            }
        }
    }
    for (index, left) in claims.iter().enumerate() {
        for right in &claims[index + 1..] {
            if left.property != right.property
                || left.support != right.support
                || left.disposition == right.disposition
                || !(left.tolerance.is_finite() && right.tolerance.is_finite())
            {
                continue;
            }
            let (witnessed_tolerance, contradicted_tolerance) =
                if left.disposition == WitnessDispositionV1::Witnessed {
                    (left.tolerance, right.tolerance)
                } else {
                    (right.tolerance, left.tolerance)
                };
            let witnessed_is_exact = canonical_f64_bits(witnessed_tolerance) == 0.0f64.to_bits();
            if left.norm != right.norm && !witnessed_is_exact {
                continue;
            }
            if witnessed_tolerance <= contradicted_tolerance {
                issues.push(SpectralAdmissionIssueV1::ContradictoryStructure {
                    property: left.property,
                    support: left.support,
                });
            }
        }
    }
    validate_exact_structure_theorems(&spec, &claims, &mut issues);
    validate_regularity_theorems(&spec, &claims, &regularity_claims, &mut issues);
    validate_ordinary_finite_spectrum(&spec, &claims, &regularity_claims, &mut issues);

    if issues.is_empty()
        && let CompletenessScopeV1::Partial { requested } = spec.requested_scope
        && let Some(total) = admitted_algebraic_cardinality(&spec, &claims, &regularity_claims)
        && requested > total
    {
        issues.push(SpectralAdmissionIssueV1::DimensionMismatch {
            left: requested,
            right: total,
        });
    }

    if !issues.is_empty() {
        return Err(SpectralAdmissionReportV1::new(issues));
    }
    spec.structures = StructureProfileV1::new(claims.clone());
    spec.regularity = RegularityProfileV1::new(regularity_claims.clone());
    let problem_receipt =
        canonical_problem_receipt(&spec, &claims, &regularity_claims).map_err(|error| {
            SpectralAdmissionReportV1::new(vec![SpectralAdmissionIssueV1::Identity(error)])
        })?;
    Ok(ValidatedSpectralProblemV1 {
        spec,
        canonical_claims: claims,
        canonical_regularity: regularity_claims,
        problem_receipt,
    })
}

fn require_structure(
    problem: &ValidatedSpectralProblemV1,
    method: SpectralMethodClassV1,
    property: StructurePropertyV1,
    support: StructureSupportRequirementV1,
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) -> Option<StructureSupportV1> {
    let matching = || {
        problem
            .canonical_claims
            .iter()
            .filter(|claim| claim.property == property && support.accepts(claim.support))
    };
    if matching().next().is_none() {
        issues.push(SpectralAdmissionIssueV1::MissingStructureWitness {
            method,
            property,
            support,
        });
        return None;
    }

    for candidate in matching() {
        let candidate_support = candidate.support;
        let on_support = || matching().filter(move |claim| claim.support == candidate_support);
        let contradicted =
            on_support().any(|claim| claim.disposition == WitnessDispositionV1::Contradicted);
        let exactly_witnessed = on_support().any(|claim| {
            claim.disposition == WitnessDispositionV1::Witnessed
                && canonical_f64_bits(claim.tolerance) == 0.0f64.to_bits()
        });
        if exactly_witnessed && !contradicted {
            return Some(candidate_support);
        }
    }

    let has_uncontradicted_support = matching().any(|candidate| {
        let candidate_support = candidate.support;
        !matching().any(|claim| {
            claim.support == candidate_support
                && claim.disposition == WitnessDispositionV1::Contradicted
        })
    });
    if has_uncontradicted_support {
        issues.push(SpectralAdmissionIssueV1::ExactStructureWitnessRequired { method, property });
    } else {
        issues.push(SpectralAdmissionIssueV1::ContradictedMethodObligation { method, property });
    }
    None
}

fn require_regularity(
    problem: &ValidatedSpectralProblemV1,
    method: SpectralMethodClassV1,
    required: RegularityClassV1,
    issues: &mut Vec<SpectralAdmissionIssueV1>,
) {
    let satisfied = regularity_is_admitted(
        &problem.spec,
        &problem.canonical_claims,
        &problem.canonical_regularity,
        required,
    );
    if !satisfied {
        issues.push(SpectralAdmissionIssueV1::RegularityRequired { method, required });
    }
}

/// Assess one algorithm family's structural obligations. This function never
/// chooses between implementations and never asserts convergence; it returns
/// a sealed method-class token only when every minimum obligation is present.
///
/// # Errors
///
/// Returns [`SpectralAdmissionReportV1`] when the validated problem lacks or
/// contradicts any representation, origin, descriptor, space, structure, or
/// regularity obligation of `method`.
#[allow(clippy::too_many_lines)] // The exhaustive method-obligation matrix is clearest and safest in one match.
#[must_use = "the method-obligation decision must be handled"]
pub fn assess_method_class(
    problem: &ValidatedSpectralProblemV1,
    method: SpectralMethodClassV1,
) -> Result<AdmittedSpectralMethodClassV1, SpectralAdmissionReportV1> {
    let mut issues = Vec::new();
    let representation = problem.spec.class.representation;
    let descriptor = problem.spec.class.descriptor;
    let origin = problem.spec.class.origin;
    let metric = problem.spec.spaces.domain;
    let ordinary_direct = matches!(descriptor, DescriptorRoleV1::Ordinary)
        && matches!(origin, SpectralOperatorOriginV1::Direct);
    let mut selected_support = None;

    match method {
        SpectralMethodClassV1::SelfAdjointLanczos => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if !metric.definiteness.is_positive_definite() {
                issues.push(SpectralAdmissionIssueV1::PositiveDefiniteMetricRequired { method });
            }
            if problem.spec.spaces.domain.id != problem.spec.spaces.codomain.id {
                issues.push(SpectralAdmissionIssueV1::MethodSpaceMismatch { method });
            }
            selected_support = require_structure(
                problem,
                method,
                StructurePropertyV1::SelfAdjoint,
                StructureSupportRequirementV1::Exact(StructureSupportV1::InnerProduct(metric.id)),
                &mut issues,
            );
            require_regularity(
                problem,
                method,
                RegularityClassV1::FiniteDimensional,
                &mut issues,
            );
        }
        SpectralMethodClassV1::GeneralizedSelfAdjointLanczos => {
            if representation != SpectralRepresentationV1::GeneralizedPencil {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !ordinary_direct {
                if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                    issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
                }
                if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                    issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
                }
            }
            if !metric.definiteness.is_positive_definite() {
                issues.push(SpectralAdmissionIssueV1::PositiveDefiniteMetricRequired { method });
            }
            if problem.spec.spaces.domain.id != problem.spec.spaces.codomain.id {
                issues.push(SpectralAdmissionIssueV1::MethodSpaceMismatch { method });
            }
            selected_support = require_structure(
                problem,
                method,
                StructurePropertyV1::HermitianDefinitePencil,
                StructureSupportRequirementV1::Exact(StructureSupportV1::InnerProduct(metric.id)),
                &mut issues,
            );
            require_regularity(
                problem,
                method,
                RegularityClassV1::RegularPencil,
                &mut issues,
            );
        }
        SpectralMethodClassV1::GeneralArnoldi => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !ordinary_direct {
                if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                    issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
                }
                if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                    issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
                }
            }
            require_regularity(
                problem,
                method,
                RegularityClassV1::FiniteDimensional,
                &mut issues,
            );
        }
        SpectralMethodClassV1::PolynomialKrylov => {
            let grade = match representation {
                SpectralRepresentationV1::MatrixPolynomial { grade } => Some(grade),
                _ => None,
            };
            if grade.is_none() {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if let DescriptorRoleV1::Descriptor { infinity_policy } = descriptor {
                if infinity_policy == InfiniteEigenvaluePolicyV1::NoClaim {
                    issues.push(SpectralAdmissionIssueV1::DescriptorInfinityPolicyRequired {
                        method,
                    });
                }
                require_regularity(
                    problem,
                    method,
                    RegularityClassV1::RegularDescriptor,
                    &mut issues,
                );
            }
            if let Some(grade) = grade {
                require_regularity(
                    problem,
                    method,
                    RegularityClassV1::RegularPolynomial { grade },
                    &mut issues,
                );
            }
        }
        SpectralMethodClassV1::HamiltonianStructurePreserving => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if !metric.dimension.is_multiple_of(2) {
                issues.push(SpectralAdmissionIssueV1::EvenDimensionRequired { method });
            }
            if problem.spec.spaces.domain.id != problem.spec.spaces.codomain.id {
                issues.push(SpectralAdmissionIssueV1::MethodSpaceMismatch { method });
            }
            selected_support = require_structure(
                problem,
                method,
                StructurePropertyV1::Hamiltonian,
                StructureSupportRequirementV1::SymplecticForm,
                &mut issues,
            );
            require_regularity(
                problem,
                method,
                RegularityClassV1::FiniteDimensional,
                &mut issues,
            );
        }
        SpectralMethodClassV1::SymplecticStructurePreserving => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if !matches!(
                origin,
                SpectralOperatorOriginV1::Direct
                    | SpectralOperatorOriginV1::MonodromyFloquet { .. }
            ) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if !metric.dimension.is_multiple_of(2) {
                issues.push(SpectralAdmissionIssueV1::EvenDimensionRequired { method });
            }
            if problem.spec.spaces.domain.id != problem.spec.spaces.codomain.id {
                issues.push(SpectralAdmissionIssueV1::MethodSpaceMismatch { method });
            }
            selected_support = require_structure(
                problem,
                method,
                StructurePropertyV1::Symplectic,
                StructureSupportRequirementV1::SymplecticForm,
                &mut issues,
            );
            if matches!(origin, SpectralOperatorOriginV1::MonodromyFloquet { .. }) {
                require_regularity(
                    problem,
                    method,
                    RegularityClassV1::WellPosedMonodromy,
                    &mut issues,
                );
            } else {
                require_regularity(
                    problem,
                    method,
                    RegularityClassV1::FiniteDimensional,
                    &mut issues,
                );
            }
        }
        SpectralMethodClassV1::KreinJOrthogonal => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if !matches!(metric.definiteness, MetricDefinitenessV1::Indefinite { .. }) {
                issues.push(SpectralAdmissionIssueV1::IndefiniteMetricRequired { method });
            }
            if problem.spec.spaces.domain.id != problem.spec.spaces.codomain.id {
                issues.push(SpectralAdmissionIssueV1::MethodSpaceMismatch { method });
            }
            selected_support = require_structure(
                problem,
                method,
                StructurePropertyV1::JSelfAdjoint,
                StructureSupportRequirementV1::KreinForm,
                &mut issues,
            );
            require_regularity(
                problem,
                method,
                RegularityClassV1::FiniteDimensional,
                &mut issues,
            );
        }
        SpectralMethodClassV1::MonodromyArnoldi => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::MonodromyFloquet { .. }) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            require_regularity(
                problem,
                method,
                RegularityClassV1::WellPosedMonodromy,
                &mut issues,
            );
        }
        SpectralMethodClassV1::DescriptorPencil => {
            if !matches!(descriptor, DescriptorRoleV1::Descriptor { .. }) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            if matches!(representation, SpectralRepresentationV1::StandardLinear) {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(origin, SpectralOperatorOriginV1::Direct) {
                issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method });
            }
            if matches!(
                descriptor,
                DescriptorRoleV1::Descriptor {
                    infinity_policy: InfiniteEigenvaluePolicyV1::NoClaim
                }
            ) {
                issues.push(SpectralAdmissionIssueV1::DescriptorInfinityPolicyRequired { method });
            }
            require_regularity(
                problem,
                method,
                RegularityClassV1::RegularDescriptor,
                &mut issues,
            );
            match representation {
                SpectralRepresentationV1::GeneralizedPencil => require_regularity(
                    problem,
                    method,
                    RegularityClassV1::RegularPencil,
                    &mut issues,
                ),
                SpectralRepresentationV1::MatrixPolynomial { grade } => require_regularity(
                    problem,
                    method,
                    RegularityClassV1::RegularPolynomial { grade },
                    &mut issues,
                ),
                SpectralRepresentationV1::StandardLinear => {}
            }
        }
        SpectralMethodClassV1::OperatorFunctionKrylov => {
            if representation != SpectralRepresentationV1::StandardLinear {
                issues.push(SpectralAdmissionIssueV1::MethodRepresentationMismatch { method });
            }
            if !matches!(descriptor, DescriptorRoleV1::Ordinary) {
                issues.push(SpectralAdmissionIssueV1::MethodDescriptorMismatch { method });
            }
            match origin {
                SpectralOperatorOriginV1::AnalyticOperatorFunction {
                    branch_policy:
                        OperatorFunctionBranchPolicyV1::SingleValued
                        | OperatorFunctionBranchPolicyV1::ExplicitBranch,
                    ..
                } => {}
                _ => issues.push(SpectralAdmissionIssueV1::MethodOriginMismatch { method }),
            }
            require_regularity(
                problem,
                method,
                RegularityClassV1::AnalyticOperatorFunction,
                &mut issues,
            );
        }
    }

    if issues.is_empty() {
        Ok(AdmittedSpectralMethodClassV1 {
            problem_id: problem.problem_id(),
            method,
            selected_support,
        })
    } else {
        Err(SpectralAdmissionReportV1::new(issues))
    }
}

/// Admit interpretation of structural-zero-sensitive spectral gaps only when
/// the validated problem carries explicit proposition-bound gauge and
/// zero-padding semantics. This does not claim that a numerical gap is healthy
/// or that a cluster is separated.
///
/// # Errors
///
/// Returns [`SpectralAdmissionReportV1`] when either convention is unresolved
/// or its declared structural-zero count contradicts the certified nullity.
#[must_use = "the gap-semantics decision must be handled"]
pub fn assess_gap_semantics(
    problem: &ValidatedSpectralProblemV1,
) -> Result<AdmittedSpectralGapSemanticsV1, SpectralAdmissionReportV1> {
    let gauge = problem.spec.spaces.gauge;
    let zero_padding = problem.spec.spaces.zero_padding;
    let mut issues = Vec::new();
    if matches!(gauge, GaugeConventionV1::Unknown) {
        issues.push(SpectralAdmissionIssueV1::GapGaugeConventionRequired);
    }
    if matches!(zero_padding, ZeroPaddingConventionV1::Unknown) {
        issues.push(SpectralAdmissionIssueV1::GapZeroPaddingConventionRequired);
    }
    let gauge_nullity = match gauge {
        GaugeConventionV1::CertifiedNone { .. } => Some(0),
        GaugeConventionV1::Fixed { nullity, .. }
        | GaugeConventionV1::Quotiented { nullity, .. } => Some(nullity),
        GaugeConventionV1::Unknown => None,
    };
    let declared_zero_count = match zero_padding {
        ZeroPaddingConventionV1::CertifiedNonePresent { .. } | ZeroPaddingConventionV1::Unknown => {
            None
        }
        ZeroPaddingConventionV1::ExplicitlyPadded { count, .. }
        | ZeroPaddingConventionV1::Omitted { count, .. } => Some(count),
    };
    if let (Some(gauge_nullity), Some(declared_zero_count)) = (gauge_nullity, declared_zero_count)
        && gauge_nullity != declared_zero_count
    {
        issues.push(SpectralAdmissionIssueV1::GapStructuralZeroCountMismatch {
            gauge_nullity,
            declared_zero_count,
        });
    }
    if issues.is_empty() {
        Ok(AdmittedSpectralGapSemanticsV1 {
            problem_id: problem.problem_id(),
            gauge,
            zero_padding,
        })
    } else {
        Err(SpectralAdmissionReportV1::new(issues))
    }
}
