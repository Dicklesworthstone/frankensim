//! Graph-bound mating-feature fit and clearance requirements.
//!
//! Version one binds caller-declared internal/external contact-feature pairs to
//! explicit basic size and signed diametral-gap envelopes. Machine-IR proves
//! graph existence and subsystem co-ownership for each declared body/feature
//! endpoint. It does not prove cylindrical geometry, coaxiality, physical
//! containment, assembly feasibility, or standards-table interpretation.

use core::fmt;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphIdV1, MachineIdError, SubsystemId,
};
use super::ManufacturingArtifactRefV1;

/// Identity/admission schema version for fit requirements.
pub const MACHINE_FIT_CLEARANCE_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum fit requirements retained by version one.
pub const MAX_MACHINE_FIT_REQUIREMENTS_V1: usize = 4_096;

const FIT_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(12 * 1_024 * 1_024, 8 * 1_024 * 1_024, 4, 4_096, 4_096);

/// Stable identity of one mating-feature fit requirement.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FitRequirementIdV1(Box<str>);

impl FitRequirementIdV1 {
    /// Admit one bounded canonical requirement key.
    ///
    /// # Errors
    /// Refuses text outside the Machine-IR canonical key grammar.
    pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
        let key = key.into();
        super::super::validate_canonical_key("fit-requirement-id", &key)?;
        Ok(Self(key.into_boxed_str()))
    }

    /// Exact canonical key retained in aggregate identity rows.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for FitRequirementIdV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.canonical_key())
    }
}

macro_rules! artifact_role {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ManufacturingArtifactRefV1);

        impl $name {
            /// Assign an admitted manufacturing artifact coordinate to this role.
            #[must_use]
            pub const fn new(artifact: ManufacturingArtifactRefV1) -> Self {
                Self(artifact)
            }

            /// Exact nominal coordinate retained in aggregate identity.
            #[must_use]
            pub const fn artifact(&self) -> &ManufacturingArtifactRefV1 {
                &self.0
            }
        }
    };
}

artifact_role!(
    /// Exact standard/model coordinate used to interpret a fit requirement.
    FitSpecificationRefV1
);
artifact_role!(
    /// Machine-readable semantic requirement source.
    FitSemanticSourceRefV1
);
artifact_role!(
    /// Presentation-only graphical coordinate carrying no semantic authority.
    FitPresentationRefV1
);

/// Explicit submitted unit for fit dimensions and diametral gaps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum FitLengthUnitV1 {
    /// Metres.
    Metre = 1,
    /// Millimetres.
    Millimetre = 2,
    /// Micrometres.
    Micrometre = 3,
    /// Nanometres.
    Nanometre = 4,
    /// International inches, exactly 0.0254 metres at the unit-definition level.
    Inch = 5,
}

impl FitLengthUnitV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Binary64 multiplier used to normalize to coherent-SI metres.
    #[must_use]
    pub const fn metres_per_unit(self) -> f64 {
        match self {
            Self::Metre => 1.0,
            Self::Millimetre => 1.0e-3,
            Self::Micrometre => 1.0e-6,
            Self::Nanometre => 1.0e-9,
            Self::Inch => 0.0254,
        }
    }

    /// Stable unit spelling.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::Metre => "m",
            Self::Millimetre => "mm",
            Self::Micrometre => "um",
            Self::Nanometre => "nm",
            Self::Inch => "in",
        }
    }
}

/// Refusal from constructing a signed unit-bearing fit length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SignedFitLengthErrorV1 {
    /// NaN or infinity was supplied.
    NonFinite,
    /// Unit normalization produced a non-finite value.
    SiNonFinite,
    /// A nonzero submitted value vanished during SI normalization.
    SiUnderflow,
}

impl SignedFitLengthErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonFinite => "FitLengthNonFinite",
            Self::SiNonFinite => "FitLengthSiNonFinite",
            Self::SiUnderflow => "FitLengthSiUnderflow",
        }
    }
}

impl fmt::Display for SignedFitLengthErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "fit length must be finite",
            Self::SiNonFinite => "fit-length SI normalization must remain finite",
            Self::SiUnderflow => "nonzero fit length vanished in SI normalization",
        })
    }
}

impl std::error::Error for SignedFitLengthErrorV1 {}

/// Canonical signed binary64 length retaining source unit and coherent-SI bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SignedFitLengthV1 {
    submitted_bits: u64,
    unit: FitLengthUnitV1,
    metres_bits: u64,
}

impl SignedFitLengthV1 {
    /// Validate, canonicalize signed zero, and normalize one signed length.
    ///
    /// # Errors
    /// Refuses a non-finite source or a nonrepresentable nonzero SI result.
    pub fn try_new(value: f64, unit: FitLengthUnitV1) -> Result<Self, SignedFitLengthErrorV1> {
        if !value.is_finite() {
            return Err(SignedFitLengthErrorV1::NonFinite);
        }
        let submitted = if value == 0.0 { 0.0 } else { value };
        let metres = submitted * unit.metres_per_unit();
        if !metres.is_finite() {
            return Err(SignedFitLengthErrorV1::SiNonFinite);
        }
        if submitted != 0.0 && metres == 0.0 {
            return Err(SignedFitLengthErrorV1::SiUnderflow);
        }
        Ok(Self {
            submitted_bits: submitted.to_bits(),
            unit,
            metres_bits: metres.to_bits(),
        })
    }

    /// Canonical value in the submitted unit.
    #[must_use]
    pub fn submitted_value(self) -> f64 {
        f64::from_bits(self.submitted_bits)
    }

    /// Exact submitted unit.
    #[must_use]
    pub const fn unit(self) -> FitLengthUnitV1 {
        self.unit
    }

    /// Coherent-SI binary64 value in metres.
    #[must_use]
    pub fn metres(self) -> f64 {
        f64::from_bits(self.metres_bits)
    }

    /// Canonical coherent-SI bits.
    #[must_use]
    pub const fn metres_bits(self) -> u64 {
        self.metres_bits
    }

    /// Whether this value is canonical positive zero.
    #[must_use]
    pub const fn is_zero(self) -> bool {
        self.submitted_bits == 0
    }

    fn append_canonical(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.submitted_bits.to_le_bytes());
        out.push(self.unit.tag());
        out.extend_from_slice(&self.metres_bits.to_le_bytes());
    }
}

/// Refusal from constructing a non-positive basic fit size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PositiveFitLengthErrorV1 {
    /// Signed-length construction failed.
    Length(SignedFitLengthErrorV1),
    /// Basic size was zero or negative.
    NonPositive,
}

impl fmt::Display for PositiveFitLengthErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Length(error) => write!(formatter, "invalid positive fit length: {error}"),
            Self::NonPositive => formatter.write_str("basic fit size must be strictly positive"),
        }
    }
}

impl std::error::Error for PositiveFitLengthErrorV1 {}

impl From<SignedFitLengthErrorV1> for PositiveFitLengthErrorV1 {
    fn from(error: SignedFitLengthErrorV1) -> Self {
        Self::Length(error)
    }
}

/// Strictly positive basic fit size.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PositiveFitLengthV1(SignedFitLengthV1);

impl PositiveFitLengthV1 {
    /// Construct one strictly positive unit-bearing basic size.
    ///
    /// # Errors
    /// Refuses non-finite, non-positive, or nonrepresentable input.
    pub fn try_new(value: f64, unit: FitLengthUnitV1) -> Result<Self, PositiveFitLengthErrorV1> {
        let length = SignedFitLengthV1::try_new(value, unit)?;
        if length.metres() <= 0.0 {
            return Err(PositiveFitLengthErrorV1::NonPositive);
        }
        Ok(Self(length))
    }

    /// Canonical value in the submitted unit.
    #[must_use]
    pub fn submitted_value(self) -> f64 {
        self.0.submitted_value()
    }

    /// Exact submitted unit.
    #[must_use]
    pub const fn unit(self) -> FitLengthUnitV1 {
        self.0.unit()
    }

    /// Coherent-SI binary64 value in metres.
    #[must_use]
    pub fn metres(self) -> f64 {
        self.0.metres()
    }

    fn append_canonical(self, out: &mut Vec<u8>) {
        self.0.append_canonical(out);
    }
}

/// Fit regime derived from the signed diametral-gap envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum FitRegimeV1 {
    /// Every admitted gap is nonnegative and at least one bound is positive.
    Clearance = 1,
    /// The admitted envelope spans negative and positive gap.
    Transition = 2,
    /// Every admitted gap is nonpositive and at least one bound is negative.
    Interference = 3,
}

impl FitRegimeV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// Refusal from constructing an invalid signed gap envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitAllowanceErrorV1 {
    /// Minimum gap exceeded maximum gap in coherent-SI metres.
    Inverted,
    /// Both bounds were exactly zero, conveying no fit regime.
    DegenerateZero,
}

impl fmt::Display for FitAllowanceErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Inverted => "minimum fit gap must not exceed maximum fit gap",
            Self::DegenerateZero => "fit gap envelope cannot be exactly [0, 0]",
        })
    }
}

impl std::error::Error for FitAllowanceErrorV1 {}

/// Ordered signed diametral-gap envelope with a derived, non-contradictory regime.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitAllowanceV1 {
    minimum_diametral_gap: SignedFitLengthV1,
    maximum_diametral_gap: SignedFitLengthV1,
    regime: FitRegimeV1,
}

impl FitAllowanceV1 {
    /// Validate an ordered gap envelope and derive its only compatible regime.
    ///
    /// # Errors
    /// Refuses inverted bounds or a degenerate all-zero envelope.
    pub fn try_new(
        minimum_diametral_gap: SignedFitLengthV1,
        maximum_diametral_gap: SignedFitLengthV1,
    ) -> Result<Self, FitAllowanceErrorV1> {
        let minimum = minimum_diametral_gap.metres();
        let maximum = maximum_diametral_gap.metres();
        if minimum > maximum {
            return Err(FitAllowanceErrorV1::Inverted);
        }
        if minimum == 0.0 && maximum == 0.0 {
            return Err(FitAllowanceErrorV1::DegenerateZero);
        }
        let regime = if minimum >= 0.0 {
            FitRegimeV1::Clearance
        } else if maximum <= 0.0 {
            FitRegimeV1::Interference
        } else {
            FitRegimeV1::Transition
        };
        Ok(Self {
            minimum_diametral_gap,
            maximum_diametral_gap,
            regime,
        })
    }

    /// Minimum signed diametral gap.
    #[must_use]
    pub const fn minimum_diametral_gap(&self) -> SignedFitLengthV1 {
        self.minimum_diametral_gap
    }

    /// Maximum signed diametral gap.
    #[must_use]
    pub const fn maximum_diametral_gap(&self) -> SignedFitLengthV1 {
        self.maximum_diametral_gap
    }

    /// Regime derived from the two gap bounds.
    #[must_use]
    pub const fn regime(&self) -> FitRegimeV1 {
        self.regime
    }

    fn append_canonical(&self, out: &mut Vec<u8>) {
        self.minimum_diametral_gap.append_canonical(out);
        self.maximum_diametral_gap.append_canonical(out);
        out.push(self.regime.tag());
    }
}

/// Caller-declared body/contact-feature endpoint of a mating pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitFeatureSelectorV1 {
    declared_body: BodyId,
    feature: ContactFeatureId,
}

impl FitFeatureSelectorV1 {
    /// Construct one authority-free feature selector.
    #[must_use]
    pub const fn new(declared_body: BodyId, feature: ContactFeatureId) -> Self {
        Self {
            declared_body,
            feature,
        }
    }

    /// Caller-declared body; containment is not proved by version one.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable contact feature selected as one fit endpoint.
    #[must_use]
    pub const fn feature(&self) -> &ContactFeatureId {
        &self.feature
    }

    fn append_canonical(&self, out: &mut Vec<u8>, role: u8) {
        out.push(role);
        append_bytes(out, self.declared_body.identity().as_bytes());
        append_bytes(out, self.declared_body.canonical_key().as_bytes());
        append_bytes(out, self.feature.identity().as_bytes());
        append_bytes(out, self.feature.canonical_key().as_bytes());
    }
}

/// Ordered internal/external endpoints of one fit requirement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitPairTargetV1 {
    internal: FitFeatureSelectorV1,
    external: FitFeatureSelectorV1,
}

impl FitPairTargetV1 {
    /// Construct one role-ordered pair. Admission rejects aliased endpoints.
    #[must_use]
    pub const fn new(internal: FitFeatureSelectorV1, external: FitFeatureSelectorV1) -> Self {
        Self { internal, external }
    }

    /// Caller-declared internal-feature role.
    #[must_use]
    pub const fn internal(&self) -> &FitFeatureSelectorV1 {
        &self.internal
    }

    /// Caller-declared external-feature role.
    #[must_use]
    pub const fn external(&self) -> &FitFeatureSelectorV1 {
        &self.external
    }
}

/// One semantic fit requirement over a role-ordered mating-feature pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FitRequirementV1 {
    id: FitRequirementIdV1,
    target: FitPairTargetV1,
    basic_size: PositiveFitLengthV1,
    allowance: FitAllowanceV1,
    specification: FitSpecificationRefV1,
    semantic_source: FitSemanticSourceRefV1,
    presentation: Option<FitPresentationRefV1>,
}

impl FitRequirementV1 {
    /// Construct one authority-free semantic fit requirement.
    #[must_use]
    pub const fn new(
        id: FitRequirementIdV1,
        target: FitPairTargetV1,
        basic_size: PositiveFitLengthV1,
        allowance: FitAllowanceV1,
        specification: FitSpecificationRefV1,
        semantic_source: FitSemanticSourceRefV1,
        presentation: Option<FitPresentationRefV1>,
    ) -> Self {
        Self {
            id,
            target,
            basic_size,
            allowance,
            specification,
            semantic_source,
            presentation,
        }
    }

    /// Stable requirement identity.
    #[must_use]
    pub const fn id(&self) -> &FitRequirementIdV1 {
        &self.id
    }

    /// Ordered internal/external target roles.
    #[must_use]
    pub const fn target(&self) -> &FitPairTargetV1 {
        &self.target
    }

    /// Strictly positive basic size.
    #[must_use]
    pub const fn basic_size(&self) -> PositiveFitLengthV1 {
        self.basic_size
    }

    /// Signed gap envelope and derived regime.
    #[must_use]
    pub const fn allowance(&self) -> &FitAllowanceV1 {
        &self.allowance
    }

    /// Exact standard/model interpretation coordinate.
    #[must_use]
    pub const fn specification(&self) -> &FitSpecificationRefV1 {
        &self.specification
    }

    /// Exact machine-readable semantic source coordinate.
    #[must_use]
    pub const fn semantic_source(&self) -> &FitSemanticSourceRefV1 {
        &self.semantic_source
    }

    /// Optional presentation-only coordinate.
    #[must_use]
    pub const fn presentation(&self) -> Option<&FitPresentationRefV1> {
        self.presentation.as_ref()
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(512);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        self.target.internal.append_canonical(&mut row, 1);
        self.target.external.append_canonical(&mut row, 2);
        self.basic_size.append_canonical(&mut row);
        self.allowance.append_canonical(&mut row);
        append_artifact(&mut row, self.specification.artifact());
        append_artifact(&mut row, self.semantic_source.artifact());
        match &self.presentation {
            Some(presentation) => {
                row.push(1);
                append_artifact(&mut row, presentation.artifact());
            }
            None => row.push(0),
        }
        row
    }
}

/// Mutable-by-construction fit-requirement draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineFitClearanceDraftV1 {
    /// Semantic requirements in non-semantic caller order.
    pub requirements: Vec<FitRequirementV1>,
}

impl MachineFitClearanceDraftV1 {
    /// Admit, canonicalize, and bind this fit catalog to one exact graph.
    ///
    /// # Errors
    /// Refuses empty/resource-overflow input, aliases, invalid graph selectors,
    /// or bounded identity publication failure.
    #[allow(clippy::result_large_err)] // Preserve exact owned IDs in refusals.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineFitClearanceV1, MachineFitClearanceAdmissionErrorV1> {
        admit_fit_clearance(self, graph)
    }
}

/// Canonical identity schema for one graph-bound fit catalog.
pub enum MachineFitClearanceIdentitySchemaV1 {}

impl CanonicalSchema for MachineFitClearanceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.fit-clearance.v1";
    const NAME: &'static str = "admitted-machine-fit-clearance";
    const VERSION: u32 = MACHINE_FIT_CLEARANCE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "one exact Machine graph plus canonical role-ordered mating-feature fit requirements";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("fit-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("requirements", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted fit catalog.
pub type MachineFitClearanceIdV1 = ProblemSemanticId<MachineFitClearanceIdentitySchemaV1>;

/// Canonically ordered graph-bound fit requirements.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineFitClearanceV1 {
    graph: MachineGraphIdV1,
    requirements: Vec<FitRequirementV1>,
    receipt: IdentityReceipt<MachineFitClearanceIdV1>,
}

impl AdmittedMachineFitClearanceV1 {
    /// Exact Machine graph extended by this catalog.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Requirements in canonical ID order.
    #[must_use]
    pub fn requirements(&self) -> &[FitRequirementV1] {
        &self.requirements
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineFitClearanceIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineFitClearanceIdV1> {
        self.receipt
    }
}

/// Which role contained an invalid graph selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitEndpointRoleV1 {
    /// Internal-feature endpoint.
    Internal,
    /// External-feature endpoint.
    External,
}

/// Structured refusal from fit-catalog admission.
#[allow(clippy::large_enum_variant)] // Preserve exact owned IDs in rich refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineFitClearanceAdmissionErrorV1 {
    /// At least one semantic fit requirement is required.
    NoRequirements,
    /// Raw submissions exceeded the fixed cap.
    RequirementLimit {
        /// Submitted requirement count.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One requirement identity appeared more than once.
    DuplicateRequirement {
        /// Repeated requirement identity.
        requirement: FitRequirementIdV1,
    },
    /// Two IDs selected the same unordered contact-feature pair.
    DuplicateFeaturePair {
        /// First requirement in canonical ID order.
        first: FitRequirementIdV1,
        /// Later requirement selecting the same pair.
        duplicate: FitRequirementIdV1,
    },
    /// Both endpoint roles named the same body.
    SameBody {
        /// Invalid requirement identity.
        requirement: FitRequirementIdV1,
        /// Reused body identity.
        body: BodyId,
    },
    /// Both endpoint roles named the same contact feature.
    SameFeature {
        /// Invalid requirement identity.
        requirement: FitRequirementIdV1,
        /// Reused feature identity.
        feature: ContactFeatureId,
    },
    /// One endpoint named a body absent from the graph.
    UnknownBody {
        /// Requirement containing the missing body.
        requirement: FitRequirementIdV1,
        /// Invalid endpoint role.
        role: FitEndpointRoleV1,
        /// Missing body identity.
        body: BodyId,
    },
    /// One endpoint named a contact feature absent from the graph.
    UnknownFeature {
        /// Requirement containing the missing feature.
        requirement: FitRequirementIdV1,
        /// Invalid endpoint role.
        role: FitEndpointRoleV1,
        /// Missing feature identity.
        feature: ContactFeatureId,
    },
    /// One declared body and feature have different subsystem owners.
    FeatureOwnerMismatch {
        /// Requirement containing the cross-owner selector.
        requirement: FitRequirementIdV1,
        /// Invalid endpoint role.
        role: FitEndpointRoleV1,
        /// Declared body identity.
        body: BodyId,
        /// Selected contact-feature identity.
        feature: ContactFeatureId,
        /// Graph owner of the body.
        body_owner: SubsystemId,
        /// Graph owner of the feature.
        feature_owner: SubsystemId,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineFitClearanceAdmissionErrorV1 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoRequirements => "MachineFitNoRequirements",
            Self::RequirementLimit { .. } => "MachineFitRequirementLimit",
            Self::DuplicateRequirement { .. } => "MachineFitDuplicateRequirement",
            Self::DuplicateFeaturePair { .. } => "MachineFitDuplicateFeaturePair",
            Self::SameBody { .. } => "MachineFitSameBody",
            Self::SameFeature { .. } => "MachineFitSameFeature",
            Self::UnknownBody { .. } => "MachineFitUnknownBody",
            Self::UnknownFeature { .. } => "MachineFitUnknownFeature",
            Self::FeatureOwnerMismatch { .. } => "MachineFitFeatureOwnerMismatch",
            Self::Identity(_) => "MachineFitIdentity",
        }
    }
}

impl fmt::Display for MachineFitClearanceAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoRequirements => formatter.write_str("fit catalog requires a requirement"),
            Self::RequirementLimit { actual, max } => {
                write!(
                    formatter,
                    "fit catalog has {actual} requirements; maximum is {max}"
                )
            }
            Self::DuplicateRequirement { requirement } => {
                write!(
                    formatter,
                    "fit requirement {requirement} appears more than once"
                )
            }
            Self::DuplicateFeaturePair { first, duplicate } => write!(
                formatter,
                "fit requirement {duplicate} aliases the unordered feature pair of {first}"
            ),
            Self::SameBody { requirement, body } => write!(
                formatter,
                "fit requirement {requirement} uses body {body} in both endpoint roles"
            ),
            Self::SameFeature {
                requirement,
                feature,
            } => write!(
                formatter,
                "fit requirement {requirement} uses feature {feature} in both endpoint roles"
            ),
            Self::UnknownBody {
                requirement,
                role,
                body,
            } => write!(
                formatter,
                "fit requirement {requirement} has unknown {role:?} body {body}"
            ),
            Self::UnknownFeature {
                requirement,
                role,
                feature,
            } => write!(
                formatter,
                "fit requirement {requirement} has unknown {role:?} feature {feature}"
            ),
            Self::FeatureOwnerMismatch {
                requirement,
                role,
                body,
                feature,
                body_owner,
                feature_owner,
            } => write!(
                formatter,
                "fit requirement {requirement} {role:?} body {body} is owned by {body_owner}, but feature {feature} is owned by {feature_owner}"
            ),
            Self::Identity(error) => write!(formatter, "fit identity refused: {error}"),
        }
    }
}

impl std::error::Error for MachineFitClearanceAdmissionErrorV1 {}

impl From<CanonicalError> for MachineFitClearanceAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::result_large_err)] // Preserve exact owned IDs in structured refusals.
fn admit_fit_clearance(
    draft: MachineFitClearanceDraftV1,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineFitClearanceV1, MachineFitClearanceAdmissionErrorV1> {
    if draft.requirements.is_empty() {
        return Err(MachineFitClearanceAdmissionErrorV1::NoRequirements);
    }
    if draft.requirements.len() > MAX_MACHINE_FIT_REQUIREMENTS_V1 {
        return Err(MachineFitClearanceAdmissionErrorV1::RequirementLimit {
            actual: draft.requirements.len(),
            max: MAX_MACHINE_FIT_REQUIREMENTS_V1,
        });
    }

    let mut requirements = draft.requirements;
    requirements.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = requirements
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(MachineFitClearanceAdmissionErrorV1::DuplicateRequirement {
            requirement: pair[0].id.clone(),
        });
    }

    let body_owners = graph
        .subsystems()
        .iter()
        .flat_map(|subsystem| {
            subsystem
                .bodies
                .iter()
                .cloned()
                .map(move |body| (body, subsystem.id.clone()))
        })
        .collect::<BTreeMap<_, _>>();
    let feature_owners = graph
        .subsystems()
        .iter()
        .flat_map(|subsystem| {
            subsystem
                .contact_features
                .iter()
                .cloned()
                .map(move |feature| (feature, subsystem.id.clone()))
        })
        .collect::<BTreeMap<_, _>>();

    let mut unordered_pairs =
        BTreeMap::<(ContactFeatureId, ContactFeatureId), FitRequirementIdV1>::new();
    for requirement in &requirements {
        let internal = requirement.target.internal();
        let external = requirement.target.external();
        if internal.declared_body == external.declared_body {
            return Err(MachineFitClearanceAdmissionErrorV1::SameBody {
                requirement: requirement.id.clone(),
                body: internal.declared_body.clone(),
            });
        }
        if internal.feature == external.feature {
            return Err(MachineFitClearanceAdmissionErrorV1::SameFeature {
                requirement: requirement.id.clone(),
                feature: internal.feature.clone(),
            });
        }
        let pair = if internal.feature < external.feature {
            (internal.feature.clone(), external.feature.clone())
        } else {
            (external.feature.clone(), internal.feature.clone())
        };
        if let Some(first) = unordered_pairs.insert(pair, requirement.id.clone()) {
            return Err(MachineFitClearanceAdmissionErrorV1::DuplicateFeaturePair {
                first,
                duplicate: requirement.id.clone(),
            });
        }
        validate_selector(
            requirement.id(),
            FitEndpointRoleV1::Internal,
            internal,
            &body_owners,
            &feature_owners,
        )?;
        validate_selector(
            requirement.id(),
            FitEndpointRoleV1::External,
            external,
            &body_owners,
            &feature_owners,
        )?;
    }

    let rows = requirements
        .iter()
        .map(FitRequirementV1::canonical_row)
        .collect::<Vec<_>>();
    let graph_id = graph.identity();
    let receipt =
        CanonicalEncoder::<MachineFitClearanceIdV1, _>::new(FIT_IDENTITY_LIMITS, NeverCancel)?
            .u64(
                Field::new(0, "fit-schema-version"),
                u64::from(MACHINE_FIT_CLEARANCE_SCHEMA_VERSION_V1),
            )?
            .u64(
                Field::new(1, "frankenscript-ir-version"),
                u64::from(IR_VERSION),
            )?
            .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
            .ordered_bytes(
                Field::new(3, "requirements"),
                rows.len() as u64,
                rows.iter().map(Vec::as_slice),
            )?
            .finish()?;

    Ok(AdmittedMachineFitClearanceV1 {
        graph: graph_id,
        requirements,
        receipt,
    })
}

#[allow(clippy::result_large_err)] // Preserve exact owned IDs in structured refusals.
fn validate_selector(
    requirement: &FitRequirementIdV1,
    role: FitEndpointRoleV1,
    selector: &FitFeatureSelectorV1,
    body_owners: &BTreeMap<BodyId, SubsystemId>,
    feature_owners: &BTreeMap<ContactFeatureId, SubsystemId>,
) -> Result<(), MachineFitClearanceAdmissionErrorV1> {
    let Some(body_owner) = body_owners.get(selector.declared_body()) else {
        return Err(MachineFitClearanceAdmissionErrorV1::UnknownBody {
            requirement: requirement.clone(),
            role,
            body: selector.declared_body().clone(),
        });
    };
    let Some(feature_owner) = feature_owners.get(selector.feature()) else {
        return Err(MachineFitClearanceAdmissionErrorV1::UnknownFeature {
            requirement: requirement.clone(),
            role,
            feature: selector.feature().clone(),
        });
    };
    if body_owner != feature_owner {
        return Err(MachineFitClearanceAdmissionErrorV1::FeatureOwnerMismatch {
            requirement: requirement.clone(),
            role,
            body: selector.declared_body().clone(),
            feature: selector.feature().clone(),
            body_owner: body_owner.clone(),
            feature_owner: feature_owner.clone(),
        });
    }
    Ok(())
}

fn append_artifact(out: &mut Vec<u8>, artifact: &ManufacturingArtifactRefV1) {
    let row = artifact.canonical_row();
    append_bytes(out, &row);
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
