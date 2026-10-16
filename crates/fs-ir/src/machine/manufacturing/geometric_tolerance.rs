//! Datum-backed geometric-tolerance control admission.
//!
//! Version one intentionally supports only flatness, parallelism, and
//! perpendicularity over caller-declared body/surface-patch selectors. It
//! binds exact graph and datum-catalog identities and validates structural
//! reference rules. It does not construct planes/axes, measure geometry, or
//! establish ASME/ISO conformance.

use core::fmt;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, MachineGraphIdV1, MachineIdError, SubsystemId, SurfacePatchId,
};
use super::ManufacturingArtifactRefV1;
use super::datum_system::{
    AdmittedMachineDatumSystemV1, DatumFeatureIdV1, DatumReferenceFrameIdV1,
};

/// Identity/admission schema version for geometric-tolerance controls.
pub const MACHINE_GEOMETRIC_TOLERANCE_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum controls retained by version one.
pub const MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1: usize = 4_096;

const GEOMETRIC_TOLERANCE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(16 * 1_024 * 1_024, 8 * 1_024 * 1_024, 5, 4_096, 4_096);

/// Stable identity of one geometric-tolerance control.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeometricToleranceControlIdV1(Box<str>);

impl GeometricToleranceControlIdV1 {
    /// Admit one bounded canonical control key.
    ///
    /// # Errors
    /// Refuses text outside the Machine-IR canonical key grammar.
    pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
        let key = key.into();
        super::super::validate_canonical_key("geometric-tolerance-control-id", &key)?;
        Ok(Self(key.into_boxed_str()))
    }

    /// Exact canonical key retained in identity rows.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for GeometricToleranceControlIdV1 {
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
            /// Assign an admitted artifact coordinate to this exact role.
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
    /// Exact standard/model coordinate used to interpret one control.
    GeometricToleranceSpecificationRefV1
);
artifact_role!(
    /// Machine-readable semantic control source.
    GeometricToleranceSemanticSourceRefV1
);
artifact_role!(
    /// Presentation-only graphical coordinate carrying no semantic authority.
    GeometricTolerancePresentationRefV1
);

/// Explicit submitted unit for one tolerance-zone width.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum GeometricToleranceLengthUnitV1 {
    /// Metres.
    Metre = 1,
    /// Millimetres.
    Millimetre = 2,
    /// Micrometres.
    Micrometre = 3,
    /// Nanometres.
    Nanometre = 4,
}

impl GeometricToleranceLengthUnitV1 {
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

/// Refusal from constructing a positive tolerance-zone width.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometricToleranceLengthErrorV1 {
    /// NaN or infinity was supplied.
    NonFinite,
    /// Zone width was zero or negative.
    NonPositive,
    /// Unit normalization produced a non-finite value.
    SiNonFinite,
    /// A positive submitted value vanished during SI normalization.
    SiUnderflow,
}

impl GeometricToleranceLengthErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonFinite => "GeometricToleranceLengthNonFinite",
            Self::NonPositive => "GeometricToleranceLengthNonPositive",
            Self::SiNonFinite => "GeometricToleranceLengthSiNonFinite",
            Self::SiUnderflow => "GeometricToleranceLengthSiUnderflow",
        }
    }
}

impl fmt::Display for GeometricToleranceLengthErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "geometric-tolerance zone width must be finite",
            Self::NonPositive => "geometric-tolerance zone width must be positive",
            Self::SiNonFinite => "geometric-tolerance SI normalization must remain finite",
            Self::SiUnderflow => "positive geometric-tolerance width vanished in SI normalization",
        })
    }
}

impl std::error::Error for GeometricToleranceLengthErrorV1 {}

/// Strictly positive binary64 length retaining source unit and coherent-SI bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GeometricToleranceLengthV1 {
    submitted_bits: u64,
    unit: GeometricToleranceLengthUnitV1,
    metres_bits: u64,
}

impl GeometricToleranceLengthV1 {
    /// Validate and normalize one positive tolerance-zone width.
    ///
    /// # Errors
    /// Refuses non-finite, non-positive, or nonrepresentable input.
    pub fn try_new(
        value: f64,
        unit: GeometricToleranceLengthUnitV1,
    ) -> Result<Self, GeometricToleranceLengthErrorV1> {
        if !value.is_finite() {
            return Err(GeometricToleranceLengthErrorV1::NonFinite);
        }
        if value <= 0.0 {
            return Err(GeometricToleranceLengthErrorV1::NonPositive);
        }
        let metres = value * unit.metres_per_unit();
        if !metres.is_finite() {
            return Err(GeometricToleranceLengthErrorV1::SiNonFinite);
        }
        if metres == 0.0 {
            return Err(GeometricToleranceLengthErrorV1::SiUnderflow);
        }
        Ok(Self {
            submitted_bits: value.to_bits(),
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
    pub const fn unit(self) -> GeometricToleranceLengthUnitV1 {
        self.unit
    }

    /// Coherent-SI binary64 value in metres.
    #[must_use]
    pub fn metres(self) -> f64 {
        f64::from_bits(self.metres_bits)
    }

    fn append_canonical(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.submitted_bits.to_le_bytes());
        out.push(self.unit.tag());
        out.extend_from_slice(&self.metres_bits.to_le_bytes());
    }
}

/// Closed geometric-characteristic subset admitted by version one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum GeometricCharacteristicV1 {
    /// Flatness of one controlled surface; no datum frame is permitted.
    Flatness = 1,
    /// Parallelism to one admitted datum frame.
    Parallelism = 2,
    /// Perpendicularity to one admitted datum frame.
    Perpendicularity = 3,
}

impl GeometricCharacteristicV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Whether the characteristic structurally requires a datum frame.
    #[must_use]
    pub const fn requires_datum_frame(self) -> bool {
        !matches!(self, Self::Flatness)
    }
}

/// One semantic geometric-tolerance control declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeometricToleranceControlV1 {
    id: GeometricToleranceControlIdV1,
    declared_body: BodyId,
    controlled_patch: SurfacePatchId,
    characteristic: GeometricCharacteristicV1,
    zone_width: GeometricToleranceLengthV1,
    datum_frame: Option<DatumReferenceFrameIdV1>,
    specification: GeometricToleranceSpecificationRefV1,
    semantic_source: GeometricToleranceSemanticSourceRefV1,
    presentation: Option<GeometricTolerancePresentationRefV1>,
}

impl GeometricToleranceControlV1 {
    /// Construct one authority-free semantic control declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub const fn new(
        id: GeometricToleranceControlIdV1,
        declared_body: BodyId,
        controlled_patch: SurfacePatchId,
        characteristic: GeometricCharacteristicV1,
        zone_width: GeometricToleranceLengthV1,
        datum_frame: Option<DatumReferenceFrameIdV1>,
        specification: GeometricToleranceSpecificationRefV1,
        semantic_source: GeometricToleranceSemanticSourceRefV1,
        presentation: Option<GeometricTolerancePresentationRefV1>,
    ) -> Self {
        Self {
            id,
            declared_body,
            controlled_patch,
            characteristic,
            zone_width,
            datum_frame,
            specification,
            semantic_source,
            presentation,
        }
    }

    /// Stable control identity.
    #[must_use]
    pub const fn id(&self) -> &GeometricToleranceControlIdV1 {
        &self.id
    }

    /// Caller-declared body; containment is not proved by version one.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable controlled surface patch.
    #[must_use]
    pub const fn controlled_patch(&self) -> &SurfacePatchId {
        &self.controlled_patch
    }

    /// Closed characteristic tag.
    #[must_use]
    pub const fn characteristic(&self) -> GeometricCharacteristicV1 {
        self.characteristic
    }

    /// Strictly positive tolerance-zone width.
    #[must_use]
    pub const fn zone_width(&self) -> GeometricToleranceLengthV1 {
        self.zone_width
    }

    /// Datum frame required by orientation controls and forbidden for flatness.
    #[must_use]
    pub const fn datum_frame(&self) -> Option<&DatumReferenceFrameIdV1> {
        self.datum_frame.as_ref()
    }

    /// Exact standard/model interpretation coordinate.
    #[must_use]
    pub const fn specification(&self) -> &GeometricToleranceSpecificationRefV1 {
        &self.specification
    }

    /// Exact machine-readable semantic source coordinate.
    #[must_use]
    pub const fn semantic_source(&self) -> &GeometricToleranceSemanticSourceRefV1 {
        &self.semantic_source
    }

    /// Optional presentation-only coordinate.
    #[must_use]
    pub const fn presentation(&self) -> Option<&GeometricTolerancePresentationRefV1> {
        self.presentation.as_ref()
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(512);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        append_bytes(&mut row, self.declared_body.identity().as_bytes());
        append_bytes(&mut row, self.declared_body.canonical_key().as_bytes());
        append_bytes(&mut row, self.controlled_patch.identity().as_bytes());
        append_bytes(&mut row, self.controlled_patch.canonical_key().as_bytes());
        row.push(self.characteristic.tag());
        self.zone_width.append_canonical(&mut row);
        match &self.datum_frame {
            Some(frame) => {
                row.push(1);
                append_bytes(&mut row, frame.canonical_key().as_bytes());
            }
            None => row.push(0),
        }
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

/// Mutable-by-construction geometric-tolerance control draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineGeometricToleranceDraftV1 {
    /// Semantic controls in non-semantic caller order.
    pub controls: Vec<GeometricToleranceControlV1>,
}

impl MachineGeometricToleranceDraftV1 {
    /// Admit and bind controls to one exact graph and datum catalog.
    ///
    /// # Errors
    /// Refuses graph/datum mismatch, empty/resource-overflow input, aliases,
    /// invalid graph selectors, invalid datum use, or identity failure.
    #[allow(clippy::result_large_err)] // Preserve exact owned IDs in refusals.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
        datum: &AdmittedMachineDatumSystemV1,
    ) -> Result<AdmittedMachineGeometricToleranceV1, MachineGeometricToleranceAdmissionErrorV1>
    {
        admit_geometric_tolerance(self, graph, datum)
    }
}

/// Canonical identity schema for one graph/datum-bound control catalog.
pub enum MachineGeometricToleranceIdentitySchemaV1 {}

impl CanonicalSchema for MachineGeometricToleranceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.geometric-tolerance.v1";
    const NAME: &'static str = "admitted-machine-geometric-tolerance";
    const VERSION: u32 = MACHINE_GEOMETRIC_TOLERANCE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "one exact Machine graph and datum catalog plus canonical geometric-tolerance controls";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("geometric-tolerance-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("datum-system", WireType::Bytes),
        FieldSpec::required("controls", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted control catalog.
pub type MachineGeometricToleranceIdV1 =
    ProblemSemanticId<MachineGeometricToleranceIdentitySchemaV1>;

/// Canonically ordered graph/datum-bound controls and complete receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineGeometricToleranceV1 {
    graph: MachineGraphIdV1,
    datum: super::datum_system::MachineDatumSystemIdV1,
    controls: Vec<GeometricToleranceControlV1>,
    receipt: IdentityReceipt<MachineGeometricToleranceIdV1>,
}

impl AdmittedMachineGeometricToleranceV1 {
    /// Exact Machine graph extended by this catalog.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact admitted datum catalog used by orientation controls.
    #[must_use]
    pub const fn datum(&self) -> super::datum_system::MachineDatumSystemIdV1 {
        self.datum
    }

    /// Controls in canonical ID order.
    #[must_use]
    pub fn controls(&self) -> &[GeometricToleranceControlV1] {
        &self.controls
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineGeometricToleranceIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineGeometricToleranceIdV1> {
        self.receipt
    }
}

/// Structural mismatch between a characteristic and its datum-frame field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometricToleranceDatumUseIssueV1 {
    /// Flatness supplied a datum frame.
    FlatnessHasDatum,
    /// Parallelism or perpendicularity omitted a datum frame.
    OrientationMissingDatum,
}

/// Structured refusal from geometric-tolerance admission.
#[allow(clippy::large_enum_variant)] // Preserve exact owned IDs in rich refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineGeometricToleranceAdmissionErrorV1 {
    /// Datum catalog and requested graph identities differ.
    GraphDatumMismatch {
        /// Graph supplied directly to admission.
        graph: MachineGraphIdV1,
        /// Graph bound by the admitted datum catalog.
        datum_graph: MachineGraphIdV1,
    },
    /// At least one semantic control is required.
    NoControls,
    /// Raw controls exceeded the fixed cap.
    ControlLimit {
        /// Submitted control count.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One control identity appeared more than once.
    DuplicateControl {
        /// Repeated control identity.
        control: GeometricToleranceControlIdV1,
    },
    /// Two identities selected the same effective structural control.
    DuplicateControlSelector {
        /// First control in canonical ID order.
        first: GeometricToleranceControlIdV1,
        /// Later control selecting the same patch/characteristic/frame.
        duplicate: GeometricToleranceControlIdV1,
    },
    /// A control named a body absent from the graph.
    UnknownBody {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Missing body identity.
        body: BodyId,
    },
    /// A control named a surface patch absent from the graph.
    UnknownSurfacePatch {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Missing surface-patch identity.
        patch: SurfacePatchId,
    },
    /// Body and surface exist but have different subsystem owners.
    SurfaceOwnerMismatch {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Declared body identity.
        body: BodyId,
        /// Controlled surface-patch identity.
        patch: SurfacePatchId,
        /// Graph owner of the body.
        body_owner: SubsystemId,
        /// Graph owner of the patch.
        patch_owner: SubsystemId,
    },
    /// Characteristic and datum-frame presence disagree.
    InvalidDatumUse {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Exact structural mismatch.
        issue: GeometricToleranceDatumUseIssueV1,
    },
    /// An orientation control referenced an absent datum frame.
    UnknownDatumFrame {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Missing frame identity.
        frame: DatumReferenceFrameIdV1,
    },
    /// An admitted frame's primary datum could not be resolved.
    DatumInvariantGap {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Referenced frame identity.
        frame: DatumReferenceFrameIdV1,
        /// Missing primary datum feature.
        primary: DatumFeatureIdV1,
    },
    /// Controlled body and datum-frame body differ.
    DatumBodyMismatch {
        /// Invalid control identity.
        control: GeometricToleranceControlIdV1,
        /// Referenced frame identity.
        frame: DatumReferenceFrameIdV1,
        /// Body declared by the control.
        controlled_body: BodyId,
        /// Body declared by the frame's primary datum feature.
        datum_body: BodyId,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineGeometricToleranceAdmissionErrorV1 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::GraphDatumMismatch { .. } => "MachineGeometricToleranceGraphDatumMismatch",
            Self::NoControls => "MachineGeometricToleranceNoControls",
            Self::ControlLimit { .. } => "MachineGeometricToleranceControlLimit",
            Self::DuplicateControl { .. } => "MachineGeometricToleranceDuplicateControl",
            Self::DuplicateControlSelector { .. } => "MachineGeometricToleranceDuplicateSelector",
            Self::UnknownBody { .. } => "MachineGeometricToleranceUnknownBody",
            Self::UnknownSurfacePatch { .. } => "MachineGeometricToleranceUnknownPatch",
            Self::SurfaceOwnerMismatch { .. } => "MachineGeometricToleranceOwnerMismatch",
            Self::InvalidDatumUse { .. } => "MachineGeometricToleranceInvalidDatumUse",
            Self::UnknownDatumFrame { .. } => "MachineGeometricToleranceUnknownDatumFrame",
            Self::DatumInvariantGap { .. } => "MachineGeometricToleranceDatumInvariantGap",
            Self::DatumBodyMismatch { .. } => "MachineGeometricToleranceDatumBodyMismatch",
            Self::Identity(_) => "MachineGeometricToleranceIdentity",
        }
    }
}

impl fmt::Display for MachineGeometricToleranceAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GraphDatumMismatch { graph, datum_graph } => write!(
                formatter,
                "geometric-tolerance graph {graph} differs from datum graph {datum_graph}"
            ),
            Self::NoControls => {
                formatter.write_str("geometric-tolerance catalog requires a control")
            }
            Self::ControlLimit { actual, max } => write!(
                formatter,
                "geometric-tolerance catalog has {actual} controls; maximum is {max}"
            ),
            Self::DuplicateControl { control } => {
                write!(
                    formatter,
                    "geometric-tolerance control {control} is repeated"
                )
            }
            Self::DuplicateControlSelector { first, duplicate } => write!(
                formatter,
                "geometric-tolerance control {duplicate} aliases the selector of {first}"
            ),
            Self::UnknownBody { control, body } => {
                write!(formatter, "control {control} names unknown body {body}")
            }
            Self::UnknownSurfacePatch { control, patch } => {
                write!(formatter, "control {control} names unknown patch {patch}")
            }
            Self::SurfaceOwnerMismatch {
                control,
                body,
                patch,
                body_owner,
                patch_owner,
            } => write!(
                formatter,
                "control {control} body {body} is owned by {body_owner}, but patch {patch} is owned by {patch_owner}"
            ),
            Self::InvalidDatumUse { control, issue } => {
                write!(
                    formatter,
                    "control {control} has invalid datum use {issue:?}"
                )
            }
            Self::UnknownDatumFrame { control, frame } => {
                write!(
                    formatter,
                    "control {control} names unknown datum frame {frame}"
                )
            }
            Self::DatumInvariantGap {
                control,
                frame,
                primary,
            } => write!(
                formatter,
                "control {control} frame {frame} has unresolved primary datum {primary}"
            ),
            Self::DatumBodyMismatch {
                control,
                frame,
                controlled_body,
                datum_body,
            } => write!(
                formatter,
                "control {control} body {controlled_body} differs from frame {frame} body {datum_body}"
            ),
            Self::Identity(error) => {
                write!(formatter, "geometric-tolerance identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for MachineGeometricToleranceAdmissionErrorV1 {}

impl From<CanonicalError> for MachineGeometricToleranceAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::result_large_err)] // Preserve exact owned IDs in rich refusals.
fn admit_geometric_tolerance(
    draft: MachineGeometricToleranceDraftV1,
    graph: &AdmittedMachineGraph,
    datum: &AdmittedMachineDatumSystemV1,
) -> Result<AdmittedMachineGeometricToleranceV1, MachineGeometricToleranceAdmissionErrorV1> {
    let graph_id = graph.identity();
    if datum.graph() != graph_id {
        return Err(
            MachineGeometricToleranceAdmissionErrorV1::GraphDatumMismatch {
                graph: graph_id,
                datum_graph: datum.graph(),
            },
        );
    }
    if draft.controls.is_empty() {
        return Err(MachineGeometricToleranceAdmissionErrorV1::NoControls);
    }
    if draft.controls.len() > MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1 {
        return Err(MachineGeometricToleranceAdmissionErrorV1::ControlLimit {
            actual: draft.controls.len(),
            max: MAX_MACHINE_GEOMETRIC_TOLERANCE_CONTROLS_V1,
        });
    }

    let mut controls = draft.controls;
    controls.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = controls.windows(2).find(|pair| pair[0].id == pair[1].id) {
        return Err(
            MachineGeometricToleranceAdmissionErrorV1::DuplicateControl {
                control: pair[0].id.clone(),
            },
        );
    }
    let mut selectors = BTreeMap::<
        (
            SurfacePatchId,
            GeometricCharacteristicV1,
            Option<DatumReferenceFrameIdV1>,
        ),
        GeometricToleranceControlIdV1,
    >::new();
    for control in &controls {
        let selector = (
            control.controlled_patch.clone(),
            control.characteristic,
            control.datum_frame.clone(),
        );
        if let Some(first) = selectors.insert(selector, control.id.clone()) {
            return Err(
                MachineGeometricToleranceAdmissionErrorV1::DuplicateControlSelector {
                    first,
                    duplicate: control.id.clone(),
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
    let datum_features = datum
        .datum_features()
        .iter()
        .map(|feature| (feature.id().clone(), feature))
        .collect::<BTreeMap<_, _>>();
    let datum_frames = datum
        .reference_frames()
        .iter()
        .map(|frame| (frame.id().clone(), frame))
        .collect::<BTreeMap<_, _>>();

    for control in &controls {
        let Some(body_owner) = body_owners.get(&control.declared_body) else {
            return Err(MachineGeometricToleranceAdmissionErrorV1::UnknownBody {
                control: control.id.clone(),
                body: control.declared_body.clone(),
            });
        };
        let Some(patch_owner) = patch_owners.get(&control.controlled_patch) else {
            return Err(
                MachineGeometricToleranceAdmissionErrorV1::UnknownSurfacePatch {
                    control: control.id.clone(),
                    patch: control.controlled_patch.clone(),
                },
            );
        };
        if body_owner != patch_owner {
            return Err(
                MachineGeometricToleranceAdmissionErrorV1::SurfaceOwnerMismatch {
                    control: control.id.clone(),
                    body: control.declared_body.clone(),
                    patch: control.controlled_patch.clone(),
                    body_owner: body_owner.clone(),
                    patch_owner: patch_owner.clone(),
                },
            );
        }

        match (
            control.characteristic.requires_datum_frame(),
            &control.datum_frame,
        ) {
            (false, Some(_)) => {
                return Err(MachineGeometricToleranceAdmissionErrorV1::InvalidDatumUse {
                    control: control.id.clone(),
                    issue: GeometricToleranceDatumUseIssueV1::FlatnessHasDatum,
                });
            }
            (true, None) => {
                return Err(MachineGeometricToleranceAdmissionErrorV1::InvalidDatumUse {
                    control: control.id.clone(),
                    issue: GeometricToleranceDatumUseIssueV1::OrientationMissingDatum,
                });
            }
            (false, None) => {}
            (true, Some(frame_id)) => {
                let Some(frame) = datum_frames.get(frame_id) else {
                    return Err(
                        MachineGeometricToleranceAdmissionErrorV1::UnknownDatumFrame {
                            control: control.id.clone(),
                            frame: frame_id.clone(),
                        },
                    );
                };
                let Some(primary) = datum_features.get(frame.primary()) else {
                    return Err(
                        MachineGeometricToleranceAdmissionErrorV1::DatumInvariantGap {
                            control: control.id.clone(),
                            frame: frame_id.clone(),
                            primary: frame.primary().clone(),
                        },
                    );
                };
                if primary.declared_body() != &control.declared_body {
                    return Err(
                        MachineGeometricToleranceAdmissionErrorV1::DatumBodyMismatch {
                            control: control.id.clone(),
                            frame: frame_id.clone(),
                            controlled_body: control.declared_body.clone(),
                            datum_body: primary.declared_body().clone(),
                        },
                    );
                }
            }
        }
    }

    let rows = controls
        .iter()
        .map(GeometricToleranceControlV1::canonical_row)
        .collect::<Vec<_>>();
    let datum_id = datum.identity();
    let receipt = CanonicalEncoder::<MachineGeometricToleranceIdV1, _>::new(
        GEOMETRIC_TOLERANCE_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "geometric-tolerance-schema-version"),
        u64::from(MACHINE_GEOMETRIC_TOLERANCE_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
    .bytes(Field::new(3, "datum-system"), datum_id.as_bytes())?
    .ordered_bytes(
        Field::new(4, "controls"),
        rows.len() as u64,
        rows.iter().map(Vec::as_slice),
    )?
    .finish()?;

    Ok(AdmittedMachineGeometricToleranceV1 {
        graph: graph_id,
        datum: datum_id,
        controls,
        receipt,
    })
}

fn append_artifact(out: &mut Vec<u8>, artifact: &ManufacturingArtifactRefV1) {
    let row = artifact.canonical_row();
    append_bytes(out, &row);
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
