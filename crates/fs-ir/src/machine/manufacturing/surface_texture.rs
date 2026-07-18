//! Graph-bound surface-texture requirements and as-built observations.
//!
//! Design requirements and measured observations are different public types.
//! Both retain explicit metric, unit, filter/evaluation, artifact, and graph
//! context; neither can silently promote a graphical annotation or measurement
//! into a verified manufacturing-conformance claim.

use core::fmt;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, FrameBinding, MachineGraphIdV1, MachineIdError,
    OrientationParity, SubsystemId, SurfacePatchId,
};
use super::ManufacturingArtifactRefV1;

/// Identity/admission schema version for surface-texture state.
pub const MACHINE_SURFACE_TEXTURE_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum design requirements retained by version one.
pub const MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1: usize = 4_096;
/// Maximum measured observations retained by version one.
pub const MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1: usize = 4_096;

const SURFACE_TEXTURE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(16 * 1_024 * 1_024, 8 * 1_024 * 1_024, 5, 12_288, 8_192);

macro_rules! texture_key {
    ($(#[$meta:meta])* $name:ident, $role:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(Box<str>);

        impl $name {
            /// Admit one bounded canonical key.
            ///
            /// # Errors
            /// Refuses text outside the Machine-IR canonical key grammar.
            pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
                let key = key.into();
                super::super::validate_canonical_key($role, &key)?;
                Ok(Self(key.into_boxed_str()))
            }

            /// Exact canonical key retained in aggregate identity rows.
            #[must_use]
            pub fn canonical_key(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str(self.canonical_key())
            }
        }
    };
}

texture_key!(
    /// Stable identity of one semantic surface-texture design requirement.
    SurfaceTextureRequirementIdV1,
    "surface-texture-requirement-id"
);
texture_key!(
    /// Stable identity of one measured as-built surface-texture observation.
    SurfaceTextureObservationIdV1,
    "surface-texture-observation-id"
);

macro_rules! artifact_role {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ManufacturingArtifactRefV1);

        impl $name {
            /// Assign an already-admitted manufacturing artifact coordinate to
            /// this exact semantic role.
            #[must_use]
            pub const fn new(artifact: ManufacturingArtifactRefV1) -> Self {
                Self(artifact)
            }

            /// Exact nominal artifact coordinate retained in aggregate identity.
            #[must_use]
            pub const fn artifact(&self) -> &ManufacturingArtifactRefV1 {
                &self.0
            }
        }
    };
}

artifact_role!(
    /// Filter/evaluation interpretation used by both requirement and observation.
    SurfaceTextureFilterRefV1
);
artifact_role!(
    /// Standard/model coordinate used to interpret a design requirement.
    SurfaceTextureStandardRefV1
);
artifact_role!(
    /// Machine-readable semantic source; presentation links cannot inhabit this role.
    SurfaceTextureSemanticSourceRefV1
);
artifact_role!(
    /// Presentation-only graphical coordinate carrying no semantic authority.
    SurfaceTexturePresentationRefV1
);
artifact_role!(
    /// Measurement-result coordinate carrying no automatic metrology authority.
    SurfaceTextureMeasurementRefV1
);
artifact_role!(
    /// Calibration/context coordinate carrying no automatic traceability authority.
    SurfaceTextureCalibrationRefV1
);

/// Explicit submitted unit for one surface-texture length scalar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SurfaceTextureLengthUnitV1 {
    /// Metres.
    Metre = 1,
    /// Millimetres.
    Millimetre = 2,
    /// Micrometres.
    Micrometre = 3,
    /// Nanometres.
    Nanometre = 4,
}

impl SurfaceTextureLengthUnitV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Coherent-SI metre multiplier.
    #[must_use]
    pub const fn metres_per_unit(self) -> f64 {
        match self {
            Self::Metre => 1.0,
            Self::Millimetre => 1.0e-3,
            Self::Micrometre => 1.0e-6,
            Self::Nanometre => 1.0e-9,
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
        }
    }
}

/// Refusal from constructing a nonnegative unit-bearing length.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceTextureLengthErrorV1 {
    /// NaN or infinity was supplied.
    NonFinite,
    /// Negative lengths are outside the admitted domain.
    Negative,
    /// Unit normalization overflowed to infinity.
    SiNonFinite,
    /// A positive submitted length vanished during SI normalization.
    SiUnderflow,
}

impl fmt::Display for SurfaceTextureLengthErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "surface-texture length must be finite",
            Self::Negative => "surface-texture length must be nonnegative",
            Self::SiNonFinite => "surface-texture SI normalization must remain finite",
            Self::SiUnderflow => "positive surface-texture length vanished in SI normalization",
        })
    }
}

impl SurfaceTextureLengthErrorV1 {
    /// Stable machine-actionable refusal code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonFinite => "SurfaceTextureLengthNonFinite",
            Self::Negative => "SurfaceTextureLengthNegative",
            Self::SiNonFinite => "SurfaceTextureLengthSiNonFinite",
            Self::SiUnderflow => "SurfaceTextureLengthSiUnderflow",
        }
    }
}

impl std::error::Error for SurfaceTextureLengthErrorV1 {}

/// Canonical nonnegative binary64 length with submitted and coherent-SI value.
///
/// The submitted unit/value are retained for round trip. `metres` is derived
/// once at construction and retained for dimensionally closed comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SurfaceTextureLengthV1 {
    submitted_bits: u64,
    unit: SurfaceTextureLengthUnitV1,
    metres_bits: u64,
}

impl SurfaceTextureLengthV1 {
    /// Validate, canonicalize signed zero, and normalize one length.
    ///
    /// # Errors
    /// Refuses non-finite/negative input or a nonrepresentable positive SI
    /// value.
    pub fn try_new(
        value: f64,
        unit: SurfaceTextureLengthUnitV1,
    ) -> Result<Self, SurfaceTextureLengthErrorV1> {
        if !value.is_finite() {
            return Err(SurfaceTextureLengthErrorV1::NonFinite);
        }
        if value < 0.0 {
            return Err(SurfaceTextureLengthErrorV1::Negative);
        }
        let submitted = if value == 0.0 { 0.0 } else { value };
        let metres = submitted * unit.metres_per_unit();
        if !metres.is_finite() {
            return Err(SurfaceTextureLengthErrorV1::SiNonFinite);
        }
        if submitted != 0.0 && metres == 0.0 {
            return Err(SurfaceTextureLengthErrorV1::SiUnderflow);
        }
        Ok(Self {
            submitted_bits: submitted.to_bits(),
            unit,
            metres_bits: metres.to_bits(),
        })
    }

    /// Canonical value expressed in the submitted unit.
    #[must_use]
    pub fn submitted_value(self) -> f64 {
        f64::from_bits(self.submitted_bits)
    }

    /// Exact submitted unit.
    #[must_use]
    pub const fn unit(self) -> SurfaceTextureLengthUnitV1 {
        self.unit
    }

    /// Coherent-SI value in metres.
    #[must_use]
    pub fn metres(self) -> f64 {
        f64::from_bits(self.metres_bits)
    }

    /// Whether this length is exact positive zero after canonicalization.
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

/// Supported profile-height statistic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SurfaceTextureMetricV1 {
    /// Arithmetical mean profile height Ra.
    Ra = 1,
    /// Root-mean-square profile height Rq.
    Rq = 2,
    /// Maximum profile height Rz under the referenced evaluation procedure.
    Rz = 3,
    /// Total profile height Rt under the referenced evaluation procedure.
    Rt = 4,
}

impl SurfaceTextureMetricV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Stable diagnostic symbol.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::Ra => "Ra",
            Self::Rq => "Rq",
            Self::Rz => "Rz",
            Self::Rt => "Rt",
        }
    }
}

/// Semantic limit applied to one profile-height metric.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SurfaceTextureLimitV1 {
    /// Metric must not exceed the supplied positive limit.
    Maximum(SurfaceTextureLengthV1),
    /// Metric must not fall below the supplied positive limit.
    Minimum(SurfaceTextureLengthV1),
    /// Metric must lie between the supplied inclusive bounds.
    Range {
        /// Nonnegative lower bound.
        minimum: SurfaceTextureLengthV1,
        /// Strictly greater upper bound.
        maximum: SurfaceTextureLengthV1,
    },
}

impl SurfaceTextureLimitV1 {
    fn append_canonical(&self, out: &mut Vec<u8>) {
        match self {
            Self::Maximum(value) => {
                out.push(1);
                value.append_canonical(out);
            }
            Self::Minimum(value) => {
                out.push(2);
                value.append_canonical(out);
            }
            Self::Range { minimum, maximum } => {
                out.push(3);
                minimum.append_canonical(out);
                maximum.append_canonical(out);
            }
        }
    }
}

/// Declared dominant lay direction in the named lay frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SurfaceLayV1 {
    /// No directional lay claim.
    Unspecified = 1,
    /// Nominally parallel lay.
    Parallel = 2,
    /// Nominally perpendicular lay.
    Perpendicular = 3,
    /// Crossed lay.
    Crossed = 4,
    /// Multi-directional lay.
    Multidirectional = 5,
    /// Approximately circular lay.
    Circular = 6,
    /// Approximately radial lay.
    Radial = 7,
    /// Non-directional or particulate lay.
    Particulate = 8,
}

impl SurfaceLayV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Whether version one requires a nominal frame for this lay tag.
    #[must_use]
    pub const fn requires_frame(self) -> bool {
        !matches!(self, Self::Unspecified | Self::Particulate)
    }
}

/// Declared material-removal policy attached to a requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum SurfaceProductionRuleV1 {
    /// No material-removal restriction is asserted.
    Unspecified = 1,
    /// Material removal is required.
    MaterialRemovalRequired = 2,
    /// Material removal is prohibited.
    MaterialRemovalProhibited = 3,
}

impl SurfaceProductionRuleV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// One semantic design requirement for a graph-owned surface patch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurfaceTextureRequirementV1 {
    id: SurfaceTextureRequirementIdV1,
    declared_body: BodyId,
    patch: SurfacePatchId,
    metric: SurfaceTextureMetricV1,
    limit: SurfaceTextureLimitV1,
    filter_cutoff: SurfaceTextureLengthV1,
    evaluation_length: SurfaceTextureLengthV1,
    filter_specification: SurfaceTextureFilterRefV1,
    lay: SurfaceLayV1,
    lay_frame: Option<FrameBinding>,
    production_rule: SurfaceProductionRuleV1,
    standard_specification: SurfaceTextureStandardRefV1,
    semantic_source: SurfaceTextureSemanticSourceRefV1,
    presentation: Option<SurfaceTexturePresentationRefV1>,
}

impl SurfaceTextureRequirementV1 {
    /// Construct one authority-free semantic requirement declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        id: SurfaceTextureRequirementIdV1,
        declared_body: BodyId,
        patch: SurfacePatchId,
        metric: SurfaceTextureMetricV1,
        limit: SurfaceTextureLimitV1,
        filter_cutoff: SurfaceTextureLengthV1,
        evaluation_length: SurfaceTextureLengthV1,
        filter_specification: SurfaceTextureFilterRefV1,
        lay: SurfaceLayV1,
        lay_frame: Option<FrameBinding>,
        production_rule: SurfaceProductionRuleV1,
        standard_specification: SurfaceTextureStandardRefV1,
        semantic_source: SurfaceTextureSemanticSourceRefV1,
        presentation: Option<SurfaceTexturePresentationRefV1>,
    ) -> Self {
        Self {
            id,
            declared_body,
            patch,
            metric,
            limit,
            filter_cutoff,
            evaluation_length,
            filter_specification,
            lay,
            lay_frame,
            production_rule,
            standard_specification,
            semantic_source,
            presentation,
        }
    }

    /// Stable requirement identity.
    #[must_use]
    pub const fn id(&self) -> &SurfaceTextureRequirementIdV1 {
        &self.id
    }

    /// Caller-declared body; co-ownership is checked, containment is not proved.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable surface patch carrying the requirement.
    #[must_use]
    pub const fn patch(&self) -> &SurfacePatchId {
        &self.patch
    }

    /// Explicit profile-height metric.
    #[must_use]
    pub const fn metric(&self) -> SurfaceTextureMetricV1 {
        self.metric
    }

    /// Typed limit retaining its submitted units.
    #[must_use]
    pub const fn limit(&self) -> &SurfaceTextureLimitV1 {
        &self.limit
    }

    /// Filter cutoff length.
    #[must_use]
    pub const fn filter_cutoff(&self) -> SurfaceTextureLengthV1 {
        self.filter_cutoff
    }

    /// Evaluation length.
    #[must_use]
    pub const fn evaluation_length(&self) -> SurfaceTextureLengthV1 {
        self.evaluation_length
    }

    /// Exact external filter interpretation.
    #[must_use]
    pub const fn filter_specification(&self) -> &SurfaceTextureFilterRefV1 {
        &self.filter_specification
    }

    /// Declared lay category.
    #[must_use]
    pub const fn lay(&self) -> SurfaceLayV1 {
        self.lay
    }

    /// Named frame and orientation for the lay declaration.
    #[must_use]
    pub const fn lay_frame(&self) -> Option<&FrameBinding> {
        self.lay_frame.as_ref()
    }

    /// Material-removal policy.
    #[must_use]
    pub const fn production_rule(&self) -> SurfaceProductionRuleV1 {
        self.production_rule
    }

    /// Exact standard/model coordinate used to interpret the requirement.
    #[must_use]
    pub const fn standard_specification(&self) -> &SurfaceTextureStandardRefV1 {
        &self.standard_specification
    }

    /// Exact semantic source record, distinct from optional presentation.
    #[must_use]
    pub const fn semantic_source(&self) -> &SurfaceTextureSemanticSourceRefV1 {
        &self.semantic_source
    }

    /// Optional graphical presentation link carrying no semantic authority.
    #[must_use]
    pub const fn presentation(&self) -> Option<&SurfaceTexturePresentationRefV1> {
        self.presentation.as_ref()
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(512);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        append_bytes(&mut row, self.declared_body.identity().as_bytes());
        append_bytes(&mut row, self.declared_body.canonical_key().as_bytes());
        append_bytes(&mut row, self.patch.identity().as_bytes());
        append_bytes(&mut row, self.patch.canonical_key().as_bytes());
        row.push(self.metric.tag());
        self.limit.append_canonical(&mut row);
        self.filter_cutoff.append_canonical(&mut row);
        self.evaluation_length.append_canonical(&mut row);
        append_artifact(&mut row, self.filter_specification.artifact());
        row.push(self.lay.tag());
        match &self.lay_frame {
            Some(frame) => {
                row.push(1);
                append_bytes(&mut row, frame.canonical_key().as_bytes());
                row.push(orientation_tag(frame.orientation()));
            }
            None => row.push(0),
        }
        row.push(self.production_rule.tag());
        append_artifact(&mut row, self.standard_specification.artifact());
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

/// One measured as-built observation tied to an exact design requirement.
///
/// This record does not compute pass/fail or promote metrology authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurfaceTextureObservationV1 {
    id: SurfaceTextureObservationIdV1,
    requirement: SurfaceTextureRequirementIdV1,
    metric: SurfaceTextureMetricV1,
    filter_cutoff: SurfaceTextureLengthV1,
    evaluation_length: SurfaceTextureLengthV1,
    filter_specification: SurfaceTextureFilterRefV1,
    measured_value: SurfaceTextureLengthV1,
    standard_uncertainty: SurfaceTextureLengthV1,
    measurement_artifact: SurfaceTextureMeasurementRefV1,
    calibration_artifact: SurfaceTextureCalibrationRefV1,
}

impl SurfaceTextureObservationV1 {
    /// Construct one authority-free measured-observation declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        id: SurfaceTextureObservationIdV1,
        requirement: SurfaceTextureRequirementIdV1,
        metric: SurfaceTextureMetricV1,
        filter_cutoff: SurfaceTextureLengthV1,
        evaluation_length: SurfaceTextureLengthV1,
        filter_specification: SurfaceTextureFilterRefV1,
        measured_value: SurfaceTextureLengthV1,
        standard_uncertainty: SurfaceTextureLengthV1,
        measurement_artifact: SurfaceTextureMeasurementRefV1,
        calibration_artifact: SurfaceTextureCalibrationRefV1,
    ) -> Self {
        Self {
            id,
            requirement,
            metric,
            filter_cutoff,
            evaluation_length,
            filter_specification,
            measured_value,
            standard_uncertainty,
            measurement_artifact,
            calibration_artifact,
        }
    }

    /// Stable observation identity.
    #[must_use]
    pub const fn id(&self) -> &SurfaceTextureObservationIdV1 {
        &self.id
    }

    /// Exact design requirement this observation addresses.
    #[must_use]
    pub const fn requirement(&self) -> &SurfaceTextureRequirementIdV1 {
        &self.requirement
    }

    /// Explicit observed metric, checked against the referenced requirement.
    #[must_use]
    pub const fn metric(&self) -> SurfaceTextureMetricV1 {
        self.metric
    }

    /// Actual filter cutoff used by the observation.
    #[must_use]
    pub const fn filter_cutoff(&self) -> SurfaceTextureLengthV1 {
        self.filter_cutoff
    }

    /// Actual evaluation length used by the observation.
    #[must_use]
    pub const fn evaluation_length(&self) -> SurfaceTextureLengthV1 {
        self.evaluation_length
    }

    /// Exact observed filter/evaluation interpretation.
    #[must_use]
    pub const fn filter_specification(&self) -> &SurfaceTextureFilterRefV1 {
        &self.filter_specification
    }

    /// Observed metric value in its submitted unit.
    #[must_use]
    pub const fn measured_value(&self) -> SurfaceTextureLengthV1 {
        self.measured_value
    }

    /// Caller-supplied standard uncertainty in its submitted unit.
    #[must_use]
    pub const fn standard_uncertainty(&self) -> SurfaceTextureLengthV1 {
        self.standard_uncertainty
    }

    /// Exact measurement-result artifact coordinate.
    #[must_use]
    pub const fn measurement_artifact(&self) -> &SurfaceTextureMeasurementRefV1 {
        &self.measurement_artifact
    }

    /// Exact calibration/context artifact coordinate.
    #[must_use]
    pub const fn calibration_artifact(&self) -> &SurfaceTextureCalibrationRefV1 {
        &self.calibration_artifact
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(256);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        append_bytes(&mut row, self.requirement.canonical_key().as_bytes());
        row.push(self.metric.tag());
        self.filter_cutoff.append_canonical(&mut row);
        self.evaluation_length.append_canonical(&mut row);
        append_artifact(&mut row, self.filter_specification.artifact());
        self.measured_value.append_canonical(&mut row);
        self.standard_uncertainty.append_canonical(&mut row);
        append_artifact(&mut row, self.measurement_artifact.artifact());
        append_artifact(&mut row, self.calibration_artifact.artifact());
        row
    }
}

/// Relationship error found within one requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceTextureRequirementIssueV1 {
    /// A maximum/minimum limit was exactly zero.
    NonPositiveLimit,
    /// A range maximum was not strictly greater than its minimum in SI units.
    RangeNotIncreasing,
    /// The filter cutoff was exactly zero.
    NonPositiveFilterCutoff,
    /// The evaluation length was exactly zero.
    NonPositiveEvaluationLength,
    /// Evaluation length was shorter than the filter cutoff in SI units.
    EvaluationShorterThanCutoff,
    /// A directional lay omitted its nominal frame.
    MissingLayFrame,
    /// A lay tag that requires no frame supplied one.
    UnexpectedLayFrame,
}

impl SurfaceTextureRequirementIssueV1 {
    /// Stable subdiagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonPositiveLimit => "SurfaceTextureNonPositiveLimit",
            Self::RangeNotIncreasing => "SurfaceTextureRangeNotIncreasing",
            Self::NonPositiveFilterCutoff => "SurfaceTextureNonPositiveFilterCutoff",
            Self::NonPositiveEvaluationLength => "SurfaceTextureNonPositiveEvaluationLength",
            Self::EvaluationShorterThanCutoff => "SurfaceTextureEvaluationShorterThanCutoff",
            Self::MissingLayFrame => "SurfaceTextureMissingLayFrame",
            Self::UnexpectedLayFrame => "SurfaceTextureUnexpectedLayFrame",
        }
    }
}

/// Exact observation context that disagreed with its design requirement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceTextureObservationContextIssueV1 {
    /// The observed height statistic differs.
    Metric,
    /// The coherent-SI filter cutoff differs.
    FilterCutoff,
    /// The coherent-SI evaluation length differs.
    EvaluationLength,
    /// The exact filter/evaluation interpretation artifact differs.
    FilterSpecification,
}

impl SurfaceTextureObservationContextIssueV1 {
    /// Stable subdiagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Metric => "SurfaceTextureObservationMetricMismatch",
            Self::FilterCutoff => "SurfaceTextureObservationFilterCutoffMismatch",
            Self::EvaluationLength => "SurfaceTextureObservationEvaluationLengthMismatch",
            Self::FilterSpecification => "SurfaceTextureObservationFilterSpecificationMismatch",
        }
    }
}

/// Structured refusal from surface-texture aggregate admission.
#[allow(clippy::large_enum_variant)] // Preserve exact owned IDs in rich refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineSurfaceTextureAdmissionErrorV1 {
    /// At least one semantic design requirement is required.
    NoRequirements,
    /// Raw requirement submissions exceeded the fixed cap.
    RequirementLimit {
        /// Submitted requirement count.
        actual: usize,
        /// Maximum admitted requirement count.
        max: usize,
    },
    /// Raw observation submissions exceeded the fixed cap.
    ObservationLimit {
        /// Submitted observation count.
        actual: usize,
        /// Maximum admitted observation count.
        max: usize,
    },
    /// A requirement identity appeared more than once.
    DuplicateRequirement {
        /// Repeated requirement identity.
        requirement: SurfaceTextureRequirementIdV1,
    },
    /// Two identities selected the same surface and metric.
    DuplicateRequirementSelector {
        /// First requirement in canonical ID order.
        first: SurfaceTextureRequirementIdV1,
        /// Later requirement selecting the same surface and metric.
        duplicate: SurfaceTextureRequirementIdV1,
    },
    /// An observation identity appeared more than once.
    DuplicateObservation {
        /// Repeated observation identity.
        observation: SurfaceTextureObservationIdV1,
    },
    /// A requirement named a body absent from the graph.
    UnknownBody {
        /// Requirement containing the missing body.
        requirement: SurfaceTextureRequirementIdV1,
        /// Missing body identity.
        body: BodyId,
    },
    /// A requirement named a surface patch absent from the graph.
    UnknownSurfacePatch {
        /// Requirement containing the missing patch.
        requirement: SurfaceTextureRequirementIdV1,
        /// Missing surface-patch identity.
        patch: SurfacePatchId,
    },
    /// Body and surface exist but have different subsystem owners.
    SurfaceOwnerMismatch {
        /// Requirement containing the cross-owner selector.
        requirement: SurfaceTextureRequirementIdV1,
        /// Declared body identity.
        body: BodyId,
        /// Selected surface-patch identity.
        patch: SurfacePatchId,
        /// Graph owner of the declared body.
        body_owner: SubsystemId,
        /// Graph owner of the selected patch.
        patch_owner: SubsystemId,
    },
    /// Requirement-local length relationships were invalid.
    InvalidRequirement {
        /// Invalid requirement identity.
        requirement: SurfaceTextureRequirementIdV1,
        /// Exact invalid relationship.
        issue: SurfaceTextureRequirementIssueV1,
    },
    /// An observation referenced an undeclared design requirement.
    MissingObservationRequirement {
        /// Observation containing the missing reference.
        observation: SurfaceTextureObservationIdV1,
        /// Missing requirement identity.
        requirement: SurfaceTextureRequirementIdV1,
    },
    /// An observation's actual evaluation context disagreed with its requirement.
    ObservationContextMismatch {
        /// Observation containing the mismatched context.
        observation: SurfaceTextureObservationIdV1,
        /// Referenced requirement identity.
        requirement: SurfaceTextureRequirementIdV1,
        /// First deterministic context mismatch.
        issue: SurfaceTextureObservationContextIssueV1,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineSurfaceTextureAdmissionErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoRequirements => "MachineSurfaceTextureNoRequirements",
            Self::RequirementLimit { .. } => "MachineSurfaceTextureRequirementLimit",
            Self::ObservationLimit { .. } => "MachineSurfaceTextureObservationLimit",
            Self::DuplicateRequirement { .. } => "MachineSurfaceTextureDuplicateRequirement",
            Self::DuplicateRequirementSelector { .. } => {
                "MachineSurfaceTextureDuplicateRequirementSelector"
            }
            Self::DuplicateObservation { .. } => "MachineSurfaceTextureDuplicateObservation",
            Self::UnknownBody { .. } => "MachineSurfaceTextureUnknownBody",
            Self::UnknownSurfacePatch { .. } => "MachineSurfaceTextureUnknownPatch",
            Self::SurfaceOwnerMismatch { .. } => "MachineSurfaceTextureOwnerMismatch",
            Self::InvalidRequirement { .. } => "MachineSurfaceTextureInvalidRequirement",
            Self::MissingObservationRequirement { .. } => {
                "MachineSurfaceTextureMissingObservationRequirement"
            }
            Self::ObservationContextMismatch { .. } => {
                "MachineSurfaceTextureObservationContextMismatch"
            }
            Self::Identity(_) => "MachineSurfaceTextureIdentity",
        }
    }
}

impl fmt::Display for MachineSurfaceTextureAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoRequirements => {
                formatter.write_str("surface-texture state requires a design requirement")
            }
            Self::RequirementLimit { actual, max } => write!(
                formatter,
                "surface-texture state has {actual} requirements; maximum is {max}"
            ),
            Self::ObservationLimit { actual, max } => write!(
                formatter,
                "surface-texture state has {actual} observations; maximum is {max}"
            ),
            Self::DuplicateRequirement { requirement } => write!(
                formatter,
                "surface-texture requirement {requirement} appears more than once"
            ),
            Self::DuplicateRequirementSelector { first, duplicate } => write!(
                formatter,
                "surface-texture requirement {duplicate} aliases the surface/metric selector of {first}"
            ),
            Self::DuplicateObservation { observation } => write!(
                formatter,
                "surface-texture observation {observation} appears more than once"
            ),
            Self::UnknownBody { requirement, body } => write!(
                formatter,
                "surface-texture requirement {requirement} names unknown body {body}"
            ),
            Self::UnknownSurfacePatch { requirement, patch } => write!(
                formatter,
                "surface-texture requirement {requirement} names unknown surface {patch}"
            ),
            Self::SurfaceOwnerMismatch {
                requirement,
                body,
                patch,
                body_owner,
                patch_owner,
            } => write!(
                formatter,
                "surface-texture requirement {requirement} declares body {body} owned by \
                 {body_owner}, but patch {patch} is owned by {patch_owner}"
            ),
            Self::InvalidRequirement { requirement, issue } => write!(
                formatter,
                "surface-texture requirement {requirement} has invalid relationship {}",
                issue.code()
            ),
            Self::MissingObservationRequirement {
                observation,
                requirement,
            } => write!(
                formatter,
                "surface-texture observation {observation} names missing requirement {requirement}"
            ),
            Self::ObservationContextMismatch {
                observation,
                requirement,
                issue,
            } => write!(
                formatter,
                "surface-texture observation {observation} disagrees with requirement {requirement}: {}",
                issue.code()
            ),
            Self::Identity(error) => {
                write!(formatter, "surface-texture identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for MachineSurfaceTextureAdmissionErrorV1 {}

impl From<CanonicalError> for MachineSurfaceTextureAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

/// Mutable-by-construction surface-texture state draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineSurfaceTextureDraftV1 {
    /// Semantic design requirements in non-semantic caller order.
    pub requirements: Vec<SurfaceTextureRequirementV1>,
    /// Measured as-built observations in non-semantic caller order.
    pub observations: Vec<SurfaceTextureObservationV1>,
}

impl MachineSurfaceTextureDraftV1 {
    /// Admit, canonicalize, and bind this state to one exact graph.
    ///
    /// # Errors
    /// Refuses raw resource overflow, duplicates, unknown/cross-owner targets,
    /// invalid length relationships, missing requirement references, or bounded
    /// identity publication failure.
    #[allow(clippy::result_large_err)] // Preserve exact owned IDs in structured refusals.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineSurfaceTextureV1, MachineSurfaceTextureAdmissionErrorV1> {
        admit_surface_texture(self, graph)
    }
}

/// Canonical identity schema for graph-bound surface-texture state.
pub enum MachineSurfaceTextureIdentitySchemaV1 {}

impl CanonicalSchema for MachineSurfaceTextureIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.surface-texture.v1";
    const NAME: &'static str = "admitted-machine-surface-texture";
    const VERSION: u32 = MACHINE_SURFACE_TEXTURE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact Machine graph plus canonical semantic surface-texture requirements and separately typed measured observations";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("surface-texture-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("requirements", WireType::OrderedBytes),
        FieldSpec::required("observations", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of admitted surface-texture state.
pub type MachineSurfaceTextureIdV1 = ProblemSemanticId<MachineSurfaceTextureIdentitySchemaV1>;

/// Canonically ordered graph-bound requirements and observations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineSurfaceTextureV1 {
    graph: MachineGraphIdV1,
    requirements: Vec<SurfaceTextureRequirementV1>,
    observations: Vec<SurfaceTextureObservationV1>,
    receipt: IdentityReceipt<MachineSurfaceTextureIdV1>,
}

impl AdmittedMachineSurfaceTextureV1 {
    /// Exact Machine graph extended by this state.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Semantic requirements in canonical ID order.
    #[must_use]
    pub fn requirements(&self) -> &[SurfaceTextureRequirementV1] {
        &self.requirements
    }

    /// Measured observations in canonical ID order.
    #[must_use]
    pub fn observations(&self) -> &[SurfaceTextureObservationV1] {
        &self.observations
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineSurfaceTextureIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineSurfaceTextureIdV1> {
        self.receipt
    }
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::result_large_err)] // Preserve exact owned IDs in structured refusals.
fn admit_surface_texture(
    draft: MachineSurfaceTextureDraftV1,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineSurfaceTextureV1, MachineSurfaceTextureAdmissionErrorV1> {
    if draft.requirements.is_empty() {
        return Err(MachineSurfaceTextureAdmissionErrorV1::NoRequirements);
    }
    if draft.requirements.len() > MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1 {
        return Err(MachineSurfaceTextureAdmissionErrorV1::RequirementLimit {
            actual: draft.requirements.len(),
            max: MAX_MACHINE_SURFACE_TEXTURE_REQUIREMENTS_V1,
        });
    }
    if draft.observations.len() > MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1 {
        return Err(MachineSurfaceTextureAdmissionErrorV1::ObservationLimit {
            actual: draft.observations.len(),
            max: MAX_MACHINE_SURFACE_TEXTURE_OBSERVATIONS_V1,
        });
    }

    let mut requirements = draft.requirements;
    requirements.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = requirements
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(
            MachineSurfaceTextureAdmissionErrorV1::DuplicateRequirement {
                requirement: pair[0].id.clone(),
            },
        );
    }
    let mut selectors =
        BTreeMap::<(SurfacePatchId, SurfaceTextureMetricV1), SurfaceTextureRequirementIdV1>::new();
    for requirement in &requirements {
        let selector = (requirement.patch.clone(), requirement.metric);
        if let Some(first) = selectors.insert(selector, requirement.id.clone()) {
            return Err(
                MachineSurfaceTextureAdmissionErrorV1::DuplicateRequirementSelector {
                    first,
                    duplicate: requirement.id.clone(),
                },
            );
        }
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
    let patch_owners = graph
        .subsystems()
        .iter()
        .flat_map(|subsystem| {
            subsystem
                .surface_patches
                .iter()
                .cloned()
                .map(move |patch| (patch, subsystem.id.clone()))
        })
        .collect::<BTreeMap<_, _>>();

    for requirement in &requirements {
        let Some(body_owner) = body_owners.get(&requirement.declared_body) else {
            return Err(MachineSurfaceTextureAdmissionErrorV1::UnknownBody {
                requirement: requirement.id.clone(),
                body: requirement.declared_body.clone(),
            });
        };
        let Some(patch_owner) = patch_owners.get(&requirement.patch) else {
            return Err(MachineSurfaceTextureAdmissionErrorV1::UnknownSurfacePatch {
                requirement: requirement.id.clone(),
                patch: requirement.patch.clone(),
            });
        };
        if body_owner != patch_owner {
            return Err(
                MachineSurfaceTextureAdmissionErrorV1::SurfaceOwnerMismatch {
                    requirement: requirement.id.clone(),
                    body: requirement.declared_body.clone(),
                    patch: requirement.patch.clone(),
                    body_owner: body_owner.clone(),
                    patch_owner: patch_owner.clone(),
                },
            );
        }
        validate_requirement(requirement)?;
    }

    let requirements_by_id = requirements
        .iter()
        .map(|requirement| (requirement.id.clone(), requirement))
        .collect::<BTreeMap<_, _>>();

    let mut observations = draft.observations;
    observations.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = observations
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(
            MachineSurfaceTextureAdmissionErrorV1::DuplicateObservation {
                observation: pair[0].id.clone(),
            },
        );
    }
    for observation in &observations {
        let Some(requirement) = requirements_by_id.get(&observation.requirement) else {
            return Err(
                MachineSurfaceTextureAdmissionErrorV1::MissingObservationRequirement {
                    observation: observation.id.clone(),
                    requirement: observation.requirement.clone(),
                },
            );
        };
        validate_observation_context(observation, requirement)?;
    }

    let requirement_rows = requirements
        .iter()
        .map(SurfaceTextureRequirementV1::canonical_row)
        .collect::<Vec<_>>();
    let observation_rows = observations
        .iter()
        .map(SurfaceTextureObservationV1::canonical_row)
        .collect::<Vec<_>>();
    let graph_id = graph.identity();
    let receipt = CanonicalEncoder::<MachineSurfaceTextureIdV1, _>::new(
        SURFACE_TEXTURE_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "surface-texture-schema-version"),
        u64::from(MACHINE_SURFACE_TEXTURE_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
    .ordered_bytes(
        Field::new(3, "requirements"),
        requirement_rows.len() as u64,
        requirement_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(4, "observations"),
        observation_rows.len() as u64,
        observation_rows.iter().map(Vec::as_slice),
    )?
    .finish()?;

    Ok(AdmittedMachineSurfaceTextureV1 {
        graph: graph_id,
        requirements,
        observations,
        receipt,
    })
}

#[allow(clippy::result_large_err)] // Preserve the aggregate refusal type end to end.
fn validate_requirement(
    requirement: &SurfaceTextureRequirementV1,
) -> Result<(), MachineSurfaceTextureAdmissionErrorV1> {
    let issue = match &requirement.limit {
        SurfaceTextureLimitV1::Maximum(value) | SurfaceTextureLimitV1::Minimum(value)
            if value.is_zero() =>
        {
            Some(SurfaceTextureRequirementIssueV1::NonPositiveLimit)
        }
        SurfaceTextureLimitV1::Range { minimum, maximum }
            if maximum.metres() <= minimum.metres() =>
        {
            Some(SurfaceTextureRequirementIssueV1::RangeNotIncreasing)
        }
        _ if requirement.filter_cutoff.is_zero() => {
            Some(SurfaceTextureRequirementIssueV1::NonPositiveFilterCutoff)
        }
        _ if requirement.evaluation_length.is_zero() => {
            Some(SurfaceTextureRequirementIssueV1::NonPositiveEvaluationLength)
        }
        _ if requirement.evaluation_length.metres() < requirement.filter_cutoff.metres() => {
            Some(SurfaceTextureRequirementIssueV1::EvaluationShorterThanCutoff)
        }
        _ if !requirement.lay.requires_frame() && requirement.lay_frame.is_some() => {
            Some(SurfaceTextureRequirementIssueV1::UnexpectedLayFrame)
        }
        _ if requirement.lay.requires_frame() && requirement.lay_frame.is_none() => {
            Some(SurfaceTextureRequirementIssueV1::MissingLayFrame)
        }
        _ => None,
    };
    if let Some(issue) = issue {
        return Err(MachineSurfaceTextureAdmissionErrorV1::InvalidRequirement {
            requirement: requirement.id.clone(),
            issue,
        });
    }
    Ok(())
}

#[allow(clippy::result_large_err)] // Preserve the aggregate refusal type end to end.
fn validate_observation_context(
    observation: &SurfaceTextureObservationV1,
    requirement: &SurfaceTextureRequirementV1,
) -> Result<(), MachineSurfaceTextureAdmissionErrorV1> {
    let issue = if observation.metric != requirement.metric {
        Some(SurfaceTextureObservationContextIssueV1::Metric)
    } else if observation.filter_cutoff.metres_bits != requirement.filter_cutoff.metres_bits {
        Some(SurfaceTextureObservationContextIssueV1::FilterCutoff)
    } else if observation.evaluation_length.metres_bits != requirement.evaluation_length.metres_bits
    {
        Some(SurfaceTextureObservationContextIssueV1::EvaluationLength)
    } else if observation.filter_specification != requirement.filter_specification {
        Some(SurfaceTextureObservationContextIssueV1::FilterSpecification)
    } else {
        None
    };
    if let Some(issue) = issue {
        return Err(
            MachineSurfaceTextureAdmissionErrorV1::ObservationContextMismatch {
                observation: observation.id.clone(),
                requirement: requirement.id.clone(),
                issue,
            },
        );
    }
    Ok(())
}

fn orientation_tag(orientation: OrientationParity) -> u8 {
    match orientation {
        OrientationParity::Preserving => 1,
        OrientationParity::Reversing => 2,
    }
}

fn append_artifact(out: &mut Vec<u8>, artifact: &ManufacturingArtifactRefV1) {
    let row = artifact.canonical_row();
    append_bytes(out, &row);
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
