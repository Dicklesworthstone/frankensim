//! Versioned one-way physical-domain to spectral-problem crosswalks (RB.8a).
//!
//! Domain crates may adapt physical models DOWNWARD into [`crate::admission`]
//! problem semantics, but `fs-spectral` never imports those higher layers or
//! mutates their objects. This module records the exact source identity and
//! every convention needed to interpret the induced problem and its results.
//! Unknown, ambiguous, or lossy mappings refuse before an adapter identity is
//! minted.

use core::fmt;
use std::collections::BTreeSet;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, WireType,
};

use crate::admission::{
    DescriptorRoleV1, SpectralMetricId, SpectralNormId, SpectralOperatorOriginV1,
    SpectralProblemId, SpectralRepresentationV1, SpectralScalingId, ValidatedSpectralProblemV1,
};

/// Current wire/semantic version of the physical-domain adapter schema.
pub const SPECTRAL_ADAPTER_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum QoI crosswalks accepted before sorting or identity work.
pub const MAX_SPECTRAL_ADAPTER_QOIS_V1: usize = 256;

const ADAPTER_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 18, 1 << 18, 16, 4096, 4096);

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
}

macro_rules! opaque_adapter_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct from exact typed digest bytes. Identity is not
            /// scientific authority.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact typed digest bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }

        impl DigestBytes for $name {
            fn digest_bytes(&self) -> &[u8; 32] {
                self.as_bytes()
            }
        }
    };
}

opaque_adapter_id!(
    /// Content identity of the physical source artifact.
    PhysicalSourceArtifactIdV1
);
opaque_adapter_id!(
    /// Immutable version identity of the physical source model.
    PhysicalModelVersionIdV1
);
opaque_adapter_id!(
    /// Identity of the physical primal/state space.
    PhysicalStateSpaceIdV1
);
opaque_adapter_id!(
    /// Identity of the physical dual/test space.
    PhysicalDualSpaceIdV1
);
opaque_adapter_id!(
    /// Identity of the source unit and semantic-quantity convention.
    PhysicalUnitSystemIdV1
);
opaque_adapter_id!(
    /// Identity of a source or target frame convention.
    PhysicalFrameIdV1
);
opaque_adapter_id!(
    /// Identity of the physical metric artifact.
    PhysicalMetricIdV1
);
opaque_adapter_id!(
    /// Identity of the physical norm artifact.
    PhysicalNormIdV1
);
opaque_adapter_id!(
    /// Identity of equality/inequality constraint data.
    PhysicalConstraintSetIdV1
);
opaque_adapter_id!(
    /// Identity of a physical nullspace or gauge artifact.
    PhysicalNullspaceIdV1
);
opaque_adapter_id!(
    /// Identity of a source parameter schema and held-variable convention.
    PhysicalParameterSchemaIdV1
);
opaque_adapter_id!(
    /// Identity of the physical boundary-condition schema.
    PhysicalBoundarySchemaIdV1
);
opaque_adapter_id!(
    /// Identity of a frozen physical linearization point/state.
    PhysicalLinearizationPointIdV1
);
opaque_adapter_id!(
    /// Identity of a periodic phase, Poincare section, or event-word artifact.
    PhysicalPhaseSectionIdV1
);
opaque_adapter_id!(
    /// Identity of source-domain structure evidence retained by the adapter.
    PhysicalStructureWitnessIdV1
);
opaque_adapter_id!(
    /// Identity of a source-domain QoI.
    PhysicalQoiIdV1
);
opaque_adapter_id!(
    /// Identity of the corresponding spectral QoI/interpretation target.
    SpectralQoiIdV1
);
opaque_adapter_id!(
    /// Identity of an exact forward, inverse, quotient, or interpretation map.
    SpectralAdapterMapIdV1
);
opaque_adapter_id!(
    /// Identity of an explicit no-claim statement or non-applicability proof.
    SpectralAdapterNoClaimIdV1
);

/// Domain-separated identity schema for one admitted adapter crosswalk.
pub enum SpectralAdapterIdentitySchemaV1 {}

impl CanonicalSchema for SpectralAdapterIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-spectral.physical-adapter.v1";
    const NAME: &'static str = "physical-spectral-adapter";
    const VERSION: u32 = SPECTRAL_ADAPTER_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one-way physical source identity, conventions, exact maps, target problem, and result interpretation";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source-domain", WireType::Variant),
        FieldSpec::required("source-artifact", WireType::Bytes),
        FieldSpec::required("model-version", WireType::Bytes),
        FieldSpec::required("source-context", WireType::Bytes),
        FieldSpec::required("target-problem", WireType::Bytes),
        FieldSpec::required("crosswalk", WireType::Bytes),
        FieldSpec::required("qoi-crosswalks", WireType::CanonicalSet),
    ];
}

/// Typed identity of one canonical crosswalk. Equality is not proof that its
/// source physics is scientifically correct.
pub type SpectralAdapterIdV1 = ProblemSemanticId<SpectralAdapterIdentitySchemaV1>;

/// Admitted physical source families. Concrete domain types remain upstream.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PhysicalSourceDomainV1 {
    /// Rotordynamic mass/damping/gyroscopic/stiffness models.
    Rotordynamics,
    /// State-space, descriptor, transfer, or closed-loop control models.
    Control,
    /// Acoustic pressure/velocity or boundary-operator models.
    Acoustics,
    /// Electromagnetic field/circuit/harness models.
    Electromagnetics,
    /// Periodic smooth or hybrid orbit/monodromy models.
    PeriodicDynamics,
}

impl PhysicalSourceDomainV1 {
    const fn tag(self) -> u32 {
        match self {
            Self::Rotordynamics => 0,
            Self::Control => 1,
            Self::Acoustics => 2,
            Self::Electromagnetics => 3,
            Self::PeriodicDynamics => 4,
        }
    }
}

/// Source equation class, independently cross-checked against the validated
/// spectral problem representation and descriptor role.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicalOperatorClassV1 {
    /// Standard linear operator.
    StandardLinear,
    /// Generalized pencil, with ordinary versus descriptor semantics explicit.
    GeneralizedPencil {
        /// Whether the source equation has descriptor/infinite semantics.
        descriptor: bool,
    },
    /// Exact-grade matrix polynomial.
    MatrixPolynomial {
        /// Declared grade.
        grade: u32,
        /// Whether descriptor/infinite semantics are retained.
        descriptor: bool,
    },
}

impl PhysicalOperatorClassV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::StandardLinear => 0,
            Self::GeneralizedPencil { .. } => 1,
            Self::MatrixPolynomial { .. } => 2,
        }
    }
}

/// Explicit state of a physical semantic component.
///
/// `NotApplicable` requires a content-addressed justification. `Unknown`
/// remains representable for decoded/draft input but is never admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterBindingV1<I> {
    /// Artifact is retained through an exact named map.
    Retained {
        /// Exact source artifact.
        artifact: I,
        /// Exact source-to-spectral map.
        map: SpectralAdapterMapIdV1,
    },
    /// Component provably does not apply to this source problem.
    NotApplicable {
        /// Explicit non-applicability/no-claim artifact.
        justification: SpectralAdapterNoClaimIdV1,
    },
    /// Ambiguous or omitted draft state; validation refuses.
    Unknown,
}

/// Frozen operating-point semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearizationContextV1 {
    /// Exact physical state/parameter point about which the operator was
    /// linearized.
    Frozen {
        /// Linearization-point artifact.
        point: PhysicalLinearizationPointIdV1,
        /// Exact physical model version at which the point was defined.
        /// Validation rejects replay against any other source version.
        model_version: PhysicalModelVersionIdV1,
    },
    /// No linearization is mathematically involved, with explicit evidence.
    NotApplicable {
        /// Non-applicability artifact.
        justification: SpectralAdapterNoClaimIdV1,
    },
    /// Missing or ambiguous operating point; validation refuses.
    Unknown,
}

/// Semantic class of the physical-frame map.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameMapKindV1 {
    /// No coordinate change. Source and target frame identities must match.
    Identity,
    /// A named exact coordinate transform connects distinct frames.
    ExactTransform,
    /// The frame mapping drops information and is never admitted.
    Lossy,
    /// The relationship is omitted or ambiguous and is never admitted.
    Unknown,
}

/// Periodic phase/section semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseSectionContextV1 {
    /// Periodic source with an exact phase/section and induced map.
    Periodic {
        /// Physical phase/section/event-word artifact.
        artifact: PhysicalPhaseSectionIdV1,
        /// Exact map into the spectral origin convention.
        map: SpectralAdapterMapIdV1,
    },
    /// Source is not periodic, with explicit evidence.
    NotPeriodic {
        /// Non-periodicity/non-applicability artifact.
        justification: SpectralAdapterNoClaimIdV1,
    },
    /// Missing or ambiguous phase convention; validation refuses.
    Unknown,
}

/// Whether the forward physical-to-spectral map loses information.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdapterFidelityV1 {
    /// Exact bijective change of representation.
    ExactIsomorphism {
        /// Forward map.
        forward: SpectralAdapterMapIdV1,
        /// Exact inverse map.
        inverse: SpectralAdapterMapIdV1,
    },
    /// Exact one-way map; no inverse is promised.
    ExactOneWay {
        /// Forward map.
        forward: SpectralAdapterMapIdV1,
        /// Explicit reason an inverse is absent or intentionally unclaimed.
        no_inverse: SpectralAdapterNoClaimIdV1,
    },
    /// Exact quotient/reduction with its kernel and absent inverse explicit.
    ExactQuotient {
        /// Forward quotient/reduction map.
        forward: SpectralAdapterMapIdV1,
        /// Exact eliminated/null kernel.
        kernel: PhysicalNullspaceIdV1,
        /// Explicit no-inverse statement.
        no_inverse: SpectralAdapterNoClaimIdV1,
    },
    /// Information is discarded without an exact quotient theorem; refused.
    Lossy {
        /// Exact loss/no-claim artifact.
        reason: SpectralAdapterNoClaimIdV1,
    },
    /// Multiple inequivalent forward interpretations remain; refused.
    Ambiguous,
}

impl AdapterFidelityV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::ExactIsomorphism { .. } => 0,
            Self::ExactOneWay { .. } => 1,
            Self::ExactQuotient { .. } => 2,
            Self::Lossy { .. } => 3,
            Self::Ambiguous => 4,
        }
    }
}

/// Reverse interpretation of spectral modes or QoIs in the source domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReverseInterpretationV1 {
    /// Exact reconstruction/interpretation map.
    Exact {
        /// Exact reverse map.
        map: SpectralAdapterMapIdV1,
    },
    /// Partial interpretation with an exact no-claim boundary.
    Partial {
        /// Partial reverse map.
        map: SpectralAdapterMapIdV1,
        /// What the map does not reconstruct or claim.
        no_claim: SpectralAdapterNoClaimIdV1,
    },
    /// No reverse map exists; this is explicit rather than a fake inverse.
    Unavailable {
        /// Exact no-reverse/no-claim artifact.
        no_claim: SpectralAdapterNoClaimIdV1,
    },
}

impl ReverseInterpretationV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Exact { .. } => 0,
            Self::Partial { .. } => 1,
            Self::Unavailable { .. } => 2,
        }
    }
}

/// Exact source/target metric and norm crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MetricNormCrosswalkV1 {
    /// Source metric identity.
    pub source_metric: PhysicalMetricIdV1,
    /// Source norm identity.
    pub source_norm: PhysicalNormIdV1,
    /// Exact target domain metric.
    pub target_domain_metric: SpectralMetricId,
    /// Exact target codomain metric.
    pub target_codomain_metric: SpectralMetricId,
    /// Metric transport map.
    pub metric_map: SpectralAdapterMapIdV1,
    /// Norm transport map.
    pub norm_map: SpectralAdapterMapIdV1,
    /// Target norm/model used by downstream evidence.
    pub target_norm: SpectralNormId,
}

/// Exact unit/normalization crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UnitCrosswalkV1 {
    /// Source unit/quantity-kind schema.
    pub source_units: PhysicalUnitSystemIdV1,
    /// Target spectral scaling bundle.
    pub target_scaling: SpectralScalingId,
    /// Exact source-to-target unit map.
    pub map: SpectralAdapterMapIdV1,
}

/// Exact frame crosswalk. The target frame remains explicit because v1
/// spectral problem admission does not yet carry physical frames itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrameCrosswalkV1 {
    /// Source frame.
    pub source_frame: PhysicalFrameIdV1,
    /// Target spectral-coordinate frame.
    pub target_frame: PhysicalFrameIdV1,
    /// Exact frame map.
    pub map: SpectralAdapterMapIdV1,
    /// Whether the named map is identity, exact, lossy, or unresolved.
    pub kind: FrameMapKindV1,
}

/// Complete source-side physical context.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PhysicalSourceContextV1 {
    domain: PhysicalSourceDomainV1,
    source_artifact: PhysicalSourceArtifactIdV1,
    model_version: PhysicalModelVersionIdV1,
    operator_class: PhysicalOperatorClassV1,
    state_space: PhysicalStateSpaceIdV1,
    dual_space: PhysicalDualSpaceIdV1,
    state_dimension: u32,
    dual_dimension: u32,
    units: UnitCrosswalkV1,
    frame: FrameCrosswalkV1,
    metric_norm: MetricNormCrosswalkV1,
    constraints: AdapterBindingV1<PhysicalConstraintSetIdV1>,
    nullspace: AdapterBindingV1<PhysicalNullspaceIdV1>,
    parameters: AdapterBindingV1<PhysicalParameterSchemaIdV1>,
    boundaries: AdapterBindingV1<PhysicalBoundarySchemaIdV1>,
    linearization: LinearizationContextV1,
    phase_section: PhaseSectionContextV1,
    structure: AdapterBindingV1<PhysicalStructureWitnessIdV1>,
}

impl PhysicalSourceContextV1 {
    /// Construct one explicit source context. Validation remains separate.
    #[allow(clippy::too_many_arguments)] // Every physical convention stays explicit at the boundary.
    #[must_use]
    pub fn new(
        domain: PhysicalSourceDomainV1,
        source_artifact: PhysicalSourceArtifactIdV1,
        model_version: PhysicalModelVersionIdV1,
        operator_class: PhysicalOperatorClassV1,
        state_space: PhysicalStateSpaceIdV1,
        dual_space: PhysicalDualSpaceIdV1,
        state_dimension: u32,
        dual_dimension: u32,
        units: UnitCrosswalkV1,
        frame: FrameCrosswalkV1,
        metric_norm: MetricNormCrosswalkV1,
        constraints: AdapterBindingV1<PhysicalConstraintSetIdV1>,
        nullspace: AdapterBindingV1<PhysicalNullspaceIdV1>,
        parameters: AdapterBindingV1<PhysicalParameterSchemaIdV1>,
        boundaries: AdapterBindingV1<PhysicalBoundarySchemaIdV1>,
        linearization: LinearizationContextV1,
        phase_section: PhaseSectionContextV1,
        structure: AdapterBindingV1<PhysicalStructureWitnessIdV1>,
    ) -> Self {
        Self {
            domain,
            source_artifact,
            model_version,
            operator_class,
            state_space,
            dual_space,
            state_dimension,
            dual_dimension,
            units,
            frame,
            metric_norm,
            constraints,
            nullspace,
            parameters,
            boundaries,
            linearization,
            phase_section,
            structure,
        }
    }

    /// Source domain family.
    #[must_use]
    pub const fn domain(&self) -> PhysicalSourceDomainV1 {
        self.domain
    }

    /// Exact source artifact identity.
    #[must_use]
    pub const fn source_artifact(&self) -> PhysicalSourceArtifactIdV1 {
        self.source_artifact
    }

    /// Exact source model version.
    #[must_use]
    pub const fn model_version(&self) -> PhysicalModelVersionIdV1 {
        self.model_version
    }

    /// Source operator class.
    #[must_use]
    pub const fn operator_class(&self) -> PhysicalOperatorClassV1 {
        self.operator_class
    }

    /// Exact physical state-space identity.
    #[must_use]
    pub const fn state_space(&self) -> PhysicalStateSpaceIdV1 {
        self.state_space
    }

    /// Exact physical dual/test-space identity.
    #[must_use]
    pub const fn dual_space(&self) -> PhysicalDualSpaceIdV1 {
        self.dual_space
    }

    /// Source state dimension.
    #[must_use]
    pub const fn state_dimension(&self) -> u32 {
        self.state_dimension
    }

    /// Source dual dimension.
    #[must_use]
    pub const fn dual_dimension(&self) -> u32 {
        self.dual_dimension
    }

    /// Unit and scaling crosswalk.
    #[must_use]
    pub const fn units(&self) -> UnitCrosswalkV1 {
        self.units
    }

    /// Physical-to-spectral frame crosswalk.
    #[must_use]
    pub const fn frame(&self) -> FrameCrosswalkV1 {
        self.frame
    }

    /// Metric and norm crosswalk.
    #[must_use]
    pub const fn metric_norm(&self) -> MetricNormCrosswalkV1 {
        self.metric_norm
    }

    /// Physical constraint binding.
    #[must_use]
    pub const fn constraints(&self) -> AdapterBindingV1<PhysicalConstraintSetIdV1> {
        self.constraints
    }

    /// Physical nullspace/gauge binding.
    #[must_use]
    pub const fn nullspace(&self) -> AdapterBindingV1<PhysicalNullspaceIdV1> {
        self.nullspace
    }

    /// Parameter and held-variable schema binding.
    #[must_use]
    pub const fn parameters(&self) -> AdapterBindingV1<PhysicalParameterSchemaIdV1> {
        self.parameters
    }

    /// Boundary-condition schema binding.
    #[must_use]
    pub const fn boundaries(&self) -> AdapterBindingV1<PhysicalBoundarySchemaIdV1> {
        self.boundaries
    }

    /// Frozen linearization or explicit non-applicability state.
    #[must_use]
    pub const fn linearization(&self) -> LinearizationContextV1 {
        self.linearization
    }

    /// Periodic phase/section convention.
    #[must_use]
    pub const fn phase_section(&self) -> PhaseSectionContextV1 {
        self.phase_section
    }

    /// Source structure-witness binding.
    #[must_use]
    pub const fn structure(&self) -> AdapterBindingV1<PhysicalStructureWitnessIdV1> {
        self.structure
    }
}

/// One source QoI to spectral interpretation crosswalk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpectralQoiCrosswalkV1 {
    source: PhysicalQoiIdV1,
    target: SpectralQoiIdV1,
    forward: SpectralAdapterMapIdV1,
    reverse: ReverseInterpretationV1,
}

impl SpectralQoiCrosswalkV1 {
    /// Construct an explicit QoI crosswalk.
    #[must_use]
    pub const fn new(
        source: PhysicalQoiIdV1,
        target: SpectralQoiIdV1,
        forward: SpectralAdapterMapIdV1,
        reverse: ReverseInterpretationV1,
    ) -> Self {
        Self {
            source,
            target,
            forward,
            reverse,
        }
    }

    /// Source QoI identity.
    #[must_use]
    pub const fn source(&self) -> PhysicalQoiIdV1 {
        self.source
    }

    /// Target spectral QoI identity.
    #[must_use]
    pub const fn target(&self) -> SpectralQoiIdV1 {
        self.target
    }

    /// Exact source-to-spectral QoI map.
    #[must_use]
    pub const fn forward(&self) -> SpectralAdapterMapIdV1 {
        self.forward
    }

    /// Reverse interpretation and its no-claim boundary.
    #[must_use]
    pub const fn reverse(&self) -> ReverseInterpretationV1 {
        self.reverse
    }
}

/// Raw versioned adapter descriptor. It has no authority until validated
/// against the complete target problem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpectralAdapterSpecV1 {
    schema_version: u32,
    source: PhysicalSourceContextV1,
    target_problem: SpectralProblemId,
    fidelity: AdapterFidelityV1,
    mode_interpretation: ReverseInterpretationV1,
    qois: Vec<SpectralQoiCrosswalkV1>,
}

impl SpectralAdapterSpecV1 {
    /// Construct a current-version raw adapter descriptor.
    #[must_use]
    pub fn new(
        source: PhysicalSourceContextV1,
        target_problem: SpectralProblemId,
        fidelity: AdapterFidelityV1,
        mode_interpretation: ReverseInterpretationV1,
        qois: Vec<SpectralQoiCrosswalkV1>,
    ) -> Self {
        Self::with_schema_version(
            SPECTRAL_ADAPTER_SCHEMA_VERSION_V1,
            source,
            target_problem,
            fidelity,
            mode_interpretation,
            qois,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        source: PhysicalSourceContextV1,
        target_problem: SpectralProblemId,
        fidelity: AdapterFidelityV1,
        mode_interpretation: ReverseInterpretationV1,
        qois: Vec<SpectralQoiCrosswalkV1>,
    ) -> Self {
        Self {
            schema_version,
            source,
            target_problem,
            fidelity,
            mode_interpretation,
            qois,
        }
    }

    /// Declared adapter schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Physical source context.
    #[must_use]
    pub const fn source(&self) -> &PhysicalSourceContextV1 {
        &self.source
    }

    /// Exact validated target-problem identity expected by this draft.
    #[must_use]
    pub const fn target_problem(&self) -> SpectralProblemId {
        self.target_problem
    }

    /// Forward fidelity declaration.
    #[must_use]
    pub const fn fidelity(&self) -> AdapterFidelityV1 {
        self.fidelity
    }

    /// Mode reverse-interpretation convention.
    #[must_use]
    pub const fn mode_interpretation(&self) -> ReverseInterpretationV1 {
        self.mode_interpretation
    }

    /// Raw QoI crosswalks.
    #[must_use]
    pub fn qois(&self) -> &[SpectralQoiCrosswalkV1] {
        &self.qois
    }
}

/// Adapter field associated with an unresolved semantic binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AdapterFieldV1 {
    /// Equality/inequality constraints.
    Constraints,
    /// Nullspace/gauge data.
    Nullspace,
    /// Parameter and held-variable schema.
    Parameters,
    /// Boundary-condition schema.
    Boundaries,
    /// Linearization point.
    Linearization,
    /// Periodic phase/section.
    PhaseSection,
    /// Physical structure evidence.
    Structure,
}

/// Structured fail-closed adapter issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpectralAdapterIssueV1 {
    /// Unsupported adapter schema.
    UnsupportedSchemaVersion {
        /// Version supplied.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// Untrusted QoI collection exceeded its pre-sort cap.
    TooManyQois {
        /// Items supplied.
        found: usize,
        /// Maximum admitted.
        limit: usize,
    },
    /// Draft names a different target problem than the validated token.
    TargetProblemMismatch,
    /// Source state/dual dimension is zero or differs from the target space.
    DimensionMismatch {
        /// State or dual dimension.
        field: &'static str,
        /// Source dimension.
        source: u32,
        /// Target dimension.
        target: u32,
    },
    /// Source equation class is incompatible with the target problem class.
    OperatorClassMismatch,
    /// Target metric IDs do not match the validated problem spaces.
    MetricMismatch,
    /// Target scaling bundle does not match the validated problem.
    ScalingMismatch,
    /// An alleged identity frame map names different frames.
    FrameMismatch,
    /// Frame transport is lossy or unresolved.
    InadmissibleFrameMap,
    /// Frozen linearization belongs to a different source-model version.
    StaleLinearization,
    /// A descriptor problem omitted its retained physical constraints.
    DescriptorConstraintsMissing,
    /// A mandatory physical semantic component is unknown.
    UnknownBinding {
        /// Missing/ambiguous field.
        field: AdapterFieldV1,
    },
    /// Periodic-domain and monodromy/Floquet semantics disagree.
    PhaseOriginMismatch,
    /// Target structure claims exist but source structure was not retained.
    StructureNotRetained,
    /// Mapping is lossy or ambiguous rather than an exact admitted relation.
    InadmissibleFidelity,
    /// A one-way or quotient map claimed a whole-mode exact reverse map.
    ReverseInterpretationOverclaim,
    /// Exact quotient fidelity does not bind the same retained nullspace.
    QuotientKernelMismatch,
    /// Duplicate source or target QoI identity would make interpretation
    /// ambiguous.
    DuplicateQoi,
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Complete deterministic refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpectralAdapterReportV1 {
    issues: Vec<SpectralAdapterIssueV1>,
}

impl SpectralAdapterReportV1 {
    fn new(issues: Vec<SpectralAdapterIssueV1>) -> Self {
        Self { issues }
    }

    /// Deterministically ordered issues.
    #[must_use]
    pub fn issues(&self) -> &[SpectralAdapterIssueV1] {
        &self.issues
    }
}

impl fmt::Display for SpectralAdapterReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "physical spectral adapter refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for SpectralAdapterReportV1 {}

/// Sealed adapter validated against one complete spectral problem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedSpectralAdapterV1 {
    spec: SpectralAdapterSpecV1,
    canonical_qois: Vec<SpectralQoiCrosswalkV1>,
    receipt: IdentityReceipt<SpectralAdapterIdV1>,
}

impl ValidatedSpectralAdapterV1 {
    /// Observational raw descriptor view; it is not detachable authority.
    #[must_use]
    pub const fn spec(&self) -> &SpectralAdapterSpecV1 {
        &self.spec
    }

    /// Canonically ordered QoI crosswalks.
    #[must_use]
    pub fn qois(&self) -> &[SpectralQoiCrosswalkV1] {
        &self.canonical_qois
    }

    /// Typed semantic adapter identity.
    #[must_use]
    pub const fn adapter_id(&self) -> SpectralAdapterIdV1 {
        self.receipt.id()
    }

    /// Exact typed identity and canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<SpectralAdapterIdV1> {
        self.receipt
    }
}

/// Validate a raw physical-domain crosswalk against the exact target problem.
///
/// Resource caps are checked before sorting or identity work. No partial token
/// escapes on failure.
///
/// # Errors
/// Returns [`SpectralAdapterReportV1`] for any missing convention, target
/// mismatch, lossy/ambiguous mapping, duplicate QoI, or identity failure.
#[must_use = "the adapter admission result must be handled before use"]
pub fn validate_adapter_v1(
    mut spec: SpectralAdapterSpecV1,
    target: &ValidatedSpectralProblemV1,
) -> Result<ValidatedSpectralAdapterV1, SpectralAdapterReportV1> {
    if spec.qois.len() > MAX_SPECTRAL_ADAPTER_QOIS_V1 {
        return Err(SpectralAdapterReportV1::new(vec![
            SpectralAdapterIssueV1::TooManyQois {
                found: spec.qois.len(),
                limit: MAX_SPECTRAL_ADAPTER_QOIS_V1,
            },
        ]));
    }

    let mut issues = Vec::new();
    if spec.schema_version != SPECTRAL_ADAPTER_SCHEMA_VERSION_V1 {
        issues.push(SpectralAdapterIssueV1::UnsupportedSchemaVersion {
            found: spec.schema_version,
            supported: SPECTRAL_ADAPTER_SCHEMA_VERSION_V1,
        });
    }
    if spec.target_problem != target.problem_id() {
        issues.push(SpectralAdapterIssueV1::TargetProblemMismatch);
    }
    let target_domain = target.spec().spaces().domain();
    let target_codomain = target.spec().spaces().codomain();
    for (field, source, target_dimension) in [
        (
            "state",
            spec.source.state_dimension,
            target_domain.dimension(),
        ),
        (
            "dual",
            spec.source.dual_dimension,
            target_codomain.dimension(),
        ),
    ] {
        if source == 0 || source != target_dimension {
            issues.push(SpectralAdapterIssueV1::DimensionMismatch {
                field,
                source,
                target: target_dimension,
            });
        }
    }
    if !operator_class_matches(spec.source.operator_class, target) {
        issues.push(SpectralAdapterIssueV1::OperatorClassMismatch);
    }
    if spec.source.metric_norm.target_domain_metric != target_domain.id()
        || spec.source.metric_norm.target_codomain_metric != target_codomain.id()
    {
        issues.push(SpectralAdapterIssueV1::MetricMismatch);
    }
    if spec.source.units.target_scaling != target.spec().scaling().id() {
        issues.push(SpectralAdapterIssueV1::ScalingMismatch);
    }
    match spec.source.frame.kind {
        FrameMapKindV1::Identity
            if spec.source.frame.source_frame != spec.source.frame.target_frame =>
        {
            issues.push(SpectralAdapterIssueV1::FrameMismatch);
        }
        FrameMapKindV1::Lossy | FrameMapKindV1::Unknown => {
            issues.push(SpectralAdapterIssueV1::InadmissibleFrameMap);
        }
        FrameMapKindV1::Identity | FrameMapKindV1::ExactTransform => {}
    }

    check_binding(
        spec.source.constraints,
        AdapterFieldV1::Constraints,
        &mut issues,
    );
    check_binding(
        spec.source.nullspace,
        AdapterFieldV1::Nullspace,
        &mut issues,
    );
    check_binding(
        spec.source.parameters,
        AdapterFieldV1::Parameters,
        &mut issues,
    );
    check_binding(
        spec.source.boundaries,
        AdapterFieldV1::Boundaries,
        &mut issues,
    );
    match spec.source.linearization {
        LinearizationContextV1::Frozen { model_version, .. }
            if model_version != spec.source.model_version =>
        {
            issues.push(SpectralAdapterIssueV1::StaleLinearization);
        }
        LinearizationContextV1::Unknown => {
            issues.push(SpectralAdapterIssueV1::UnknownBinding {
                field: AdapterFieldV1::Linearization,
            });
        }
        LinearizationContextV1::Frozen { .. } | LinearizationContextV1::NotApplicable { .. } => {}
    }
    if matches!(spec.source.phase_section, PhaseSectionContextV1::Unknown) {
        issues.push(SpectralAdapterIssueV1::UnknownBinding {
            field: AdapterFieldV1::PhaseSection,
        });
    }
    check_binding(
        spec.source.structure,
        AdapterFieldV1::Structure,
        &mut issues,
    );

    if matches!(
        target.spec().class().descriptor(),
        DescriptorRoleV1::Descriptor { .. }
    ) && !matches!(spec.source.constraints, AdapterBindingV1::Retained { .. })
    {
        issues.push(SpectralAdapterIssueV1::DescriptorConstraintsMissing);
    }

    let target_is_periodic = matches!(
        target.spec().class().origin(),
        SpectralOperatorOriginV1::MonodromyFloquet { .. }
    );
    let source_is_periodic = spec.source.domain == PhysicalSourceDomainV1::PeriodicDynamics;
    let phase_is_periodic = matches!(
        spec.source.phase_section,
        PhaseSectionContextV1::Periodic { .. }
    );
    if target_is_periodic != phase_is_periodic || source_is_periodic != target_is_periodic {
        issues.push(SpectralAdapterIssueV1::PhaseOriginMismatch);
    }
    if !target.structure_claims().is_empty()
        && matches!(
            spec.source.structure,
            AdapterBindingV1::NotApplicable { .. }
        )
    {
        issues.push(SpectralAdapterIssueV1::StructureNotRetained);
    }
    if matches!(
        spec.fidelity,
        AdapterFidelityV1::Lossy { .. } | AdapterFidelityV1::Ambiguous
    ) {
        issues.push(SpectralAdapterIssueV1::InadmissibleFidelity);
    }
    if matches!(
        spec.fidelity,
        AdapterFidelityV1::ExactOneWay { .. } | AdapterFidelityV1::ExactQuotient { .. }
    ) && matches!(
        spec.mode_interpretation,
        ReverseInterpretationV1::Exact { .. }
    ) {
        issues.push(SpectralAdapterIssueV1::ReverseInterpretationOverclaim);
    }
    if let AdapterFidelityV1::ExactQuotient { kernel, .. } = spec.fidelity
        && !matches!(
            spec.source.nullspace,
            AdapterBindingV1::Retained { artifact, .. } if artifact == kernel
        )
    {
        issues.push(SpectralAdapterIssueV1::QuotientKernelMismatch);
    }

    let mut source_qois = BTreeSet::new();
    let mut target_qois = BTreeSet::new();
    for qoi in &spec.qois {
        if !source_qois.insert(*qoi.source.as_bytes())
            || !target_qois.insert(*qoi.target.as_bytes())
        {
            issues.push(SpectralAdapterIssueV1::DuplicateQoi);
            break;
        }
    }
    if !issues.is_empty() {
        return Err(SpectralAdapterReportV1::new(issues));
    }

    spec.qois.sort_by_key(qoi_bytes);
    let receipt = match adapter_receipt(&spec) {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(SpectralAdapterReportV1::new(vec![
                SpectralAdapterIssueV1::Identity(error),
            ]));
        }
    };
    Ok(ValidatedSpectralAdapterV1 {
        canonical_qois: spec.qois.clone(),
        spec,
        receipt,
    })
}

fn operator_class_matches(
    source: PhysicalOperatorClassV1,
    target: &ValidatedSpectralProblemV1,
) -> bool {
    let class = target.spec().class();
    let descriptor = matches!(class.descriptor(), DescriptorRoleV1::Descriptor { .. });
    match (source, class.representation()) {
        (PhysicalOperatorClassV1::StandardLinear, SpectralRepresentationV1::StandardLinear) => {
            !descriptor
        }
        (
            PhysicalOperatorClassV1::GeneralizedPencil {
                descriptor: source_descriptor,
            },
            SpectralRepresentationV1::GeneralizedPencil,
        ) => source_descriptor == descriptor,
        (
            PhysicalOperatorClassV1::MatrixPolynomial {
                grade: source_grade,
                descriptor: source_descriptor,
            },
            SpectralRepresentationV1::MatrixPolynomial { grade },
        ) => source_grade == grade && source_descriptor == descriptor,
        _ => false,
    }
}

fn check_binding<I: Copy>(
    binding: AdapterBindingV1<I>,
    field: AdapterFieldV1,
    issues: &mut Vec<SpectralAdapterIssueV1>,
) {
    if matches!(binding, AdapterBindingV1::Unknown) {
        issues.push(SpectralAdapterIssueV1::UnknownBinding { field });
    }
}

fn push_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_binding<I: DigestBytes>(out: &mut Vec<u8>, binding: AdapterBindingV1<I>) {
    match binding {
        AdapterBindingV1::Retained { artifact, map } => {
            out.push(0);
            out.extend_from_slice(artifact.digest_bytes());
            out.extend_from_slice(map.as_bytes());
        }
        AdapterBindingV1::NotApplicable { justification } => {
            out.push(1);
            out.extend_from_slice(justification.as_bytes());
        }
        AdapterBindingV1::Unknown => out.push(2),
    }
}

fn push_reverse(out: &mut Vec<u8>, reverse: ReverseInterpretationV1) {
    out.push(reverse.tag());
    match reverse {
        ReverseInterpretationV1::Exact { map } => out.extend_from_slice(map.as_bytes()),
        ReverseInterpretationV1::Partial { map, no_claim } => {
            out.extend_from_slice(map.as_bytes());
            out.extend_from_slice(no_claim.as_bytes());
        }
        ReverseInterpretationV1::Unavailable { no_claim } => {
            out.extend_from_slice(no_claim.as_bytes());
        }
    }
}

fn source_context_bytes(source: &PhysicalSourceContextV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(1024);
    out.extend_from_slice(source.state_space.as_bytes());
    out.extend_from_slice(source.dual_space.as_bytes());
    push_u32(&mut out, source.state_dimension);
    push_u32(&mut out, source.dual_dimension);
    out.push(source.operator_class.tag());
    match source.operator_class {
        PhysicalOperatorClassV1::StandardLinear => {}
        PhysicalOperatorClassV1::GeneralizedPencil { descriptor } => {
            out.push(u8::from(descriptor));
        }
        PhysicalOperatorClassV1::MatrixPolynomial { grade, descriptor } => {
            push_u32(&mut out, grade);
            out.push(u8::from(descriptor));
        }
    }
    out.extend_from_slice(source.units.source_units.as_bytes());
    out.extend_from_slice(source.units.target_scaling.as_bytes());
    out.extend_from_slice(source.units.map.as_bytes());
    out.extend_from_slice(source.frame.source_frame.as_bytes());
    out.extend_from_slice(source.frame.target_frame.as_bytes());
    out.extend_from_slice(source.frame.map.as_bytes());
    out.push(match source.frame.kind {
        FrameMapKindV1::Identity => 0,
        FrameMapKindV1::ExactTransform => 1,
        FrameMapKindV1::Lossy => 2,
        FrameMapKindV1::Unknown => 3,
    });
    out.extend_from_slice(source.metric_norm.source_metric.as_bytes());
    out.extend_from_slice(source.metric_norm.source_norm.as_bytes());
    out.extend_from_slice(source.metric_norm.target_domain_metric.as_bytes());
    out.extend_from_slice(source.metric_norm.target_codomain_metric.as_bytes());
    out.extend_from_slice(source.metric_norm.metric_map.as_bytes());
    out.extend_from_slice(source.metric_norm.norm_map.as_bytes());
    out.extend_from_slice(source.metric_norm.target_norm.as_bytes());
    push_binding(&mut out, source.constraints);
    push_binding(&mut out, source.nullspace);
    push_binding(&mut out, source.parameters);
    push_binding(&mut out, source.boundaries);
    match source.linearization {
        LinearizationContextV1::Frozen {
            point,
            model_version,
        } => {
            out.push(0);
            out.extend_from_slice(point.as_bytes());
            out.extend_from_slice(model_version.as_bytes());
        }
        LinearizationContextV1::NotApplicable { justification } => {
            out.push(1);
            out.extend_from_slice(justification.as_bytes());
        }
        LinearizationContextV1::Unknown => out.push(2),
    }
    match source.phase_section {
        PhaseSectionContextV1::Periodic { artifact, map } => {
            out.push(0);
            out.extend_from_slice(artifact.as_bytes());
            out.extend_from_slice(map.as_bytes());
        }
        PhaseSectionContextV1::NotPeriodic { justification } => {
            out.push(1);
            out.extend_from_slice(justification.as_bytes());
        }
        PhaseSectionContextV1::Unknown => out.push(2),
    }
    push_binding(&mut out, source.structure);
    out
}

fn crosswalk_bytes(spec: &SpectralAdapterSpecV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(256);
    out.push(spec.fidelity.tag());
    match spec.fidelity {
        AdapterFidelityV1::ExactIsomorphism { forward, inverse } => {
            out.extend_from_slice(forward.as_bytes());
            out.extend_from_slice(inverse.as_bytes());
        }
        AdapterFidelityV1::ExactOneWay {
            forward,
            no_inverse,
        } => {
            out.extend_from_slice(forward.as_bytes());
            out.extend_from_slice(no_inverse.as_bytes());
        }
        AdapterFidelityV1::ExactQuotient {
            forward,
            kernel,
            no_inverse,
        } => {
            out.extend_from_slice(forward.as_bytes());
            out.extend_from_slice(kernel.as_bytes());
            out.extend_from_slice(no_inverse.as_bytes());
        }
        AdapterFidelityV1::Lossy { reason } => out.extend_from_slice(reason.as_bytes()),
        AdapterFidelityV1::Ambiguous => {}
    }
    push_reverse(&mut out, spec.mode_interpretation);
    out
}

fn qoi_bytes(qoi: &SpectralQoiCrosswalkV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(160);
    out.extend_from_slice(qoi.source.as_bytes());
    out.extend_from_slice(qoi.target.as_bytes());
    out.extend_from_slice(qoi.forward.as_bytes());
    push_reverse(&mut out, qoi.reverse);
    out
}

fn adapter_receipt(
    spec: &SpectralAdapterSpecV1,
) -> Result<IdentityReceipt<SpectralAdapterIdV1>, CanonicalError> {
    let source_context = source_context_bytes(&spec.source);
    let crosswalk = crosswalk_bytes(spec);
    let qois: Vec<Vec<u8>> = spec.qois.iter().map(qoi_bytes).collect();
    CanonicalEncoder::<SpectralAdapterIdV1, _>::new(ADAPTER_IDENTITY_LIMITS, NeverCancel)?
        .variant(
            Field::new(0, "source-domain"),
            spec.source.domain.tag(),
            &[],
        )?
        .bytes(
            Field::new(1, "source-artifact"),
            spec.source.source_artifact.as_bytes(),
        )?
        .bytes(
            Field::new(2, "model-version"),
            spec.source.model_version.as_bytes(),
        )?
        .bytes(Field::new(3, "source-context"), &source_context)?
        .bytes(
            Field::new(4, "target-problem"),
            spec.target_problem.as_bytes(),
        )?
        .bytes(Field::new(5, "crosswalk"), &crosswalk)?
        .canonical_set(
            Field::new(6, "qoi-crosswalks"),
            qois.len() as u64,
            qois.iter().map(Vec::as_slice),
        )?
        .finish()
}
