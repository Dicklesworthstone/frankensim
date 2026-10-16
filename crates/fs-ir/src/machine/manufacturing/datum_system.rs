//! Bounded graph-bound datum-feature and reference-frame admission.
//!
//! Version one binds caller-declared bodies to durable surface or contact
//! features and records primary/secondary/tertiary precedence. Machine-IR can
//! prove that both IDs exist and share one subsystem owner; it cannot prove
//! that the selected feature is geometrically contained on the declared body
//! or construct a physical datum reference frame. Those assertions remain
//! explicit caller responsibility.

use core::fmt;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphIdV1, MachineIdError, SubsystemId,
    SurfacePatchId,
};

/// Identity/admission schema version for graph-bound datum systems.
pub const MACHINE_DATUM_SYSTEM_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum datum-feature declarations retained by version one.
pub const MAX_MACHINE_DATUM_FEATURES_V1: usize = 4_096;
/// Maximum datum reference frames retained by version one.
pub const MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1: usize = 2_048;

const DATUM_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(8 * 1_024 * 1_024, 4 * 1_024 * 1_024, 5, 8_192, 4_096);

macro_rules! datum_key {
    ($(#[$meta:meta])* $name:ident, $role:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(Box<str>);

        impl $name {
            /// Admit one bounded canonical datum key.
            ///
            /// # Errors
            /// Refuses a key outside the Machine-IR canonical key grammar.
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

datum_key!(
    /// Stable caller-facing identity of one declared datum feature.
    DatumFeatureIdV1,
    "datum-feature-id"
);
datum_key!(
    /// Stable caller-facing identity of one ordered datum reference frame.
    DatumReferenceFrameIdV1,
    "datum-reference-frame-id"
);

/// Closed durable target class admitted as a datum feature in version one.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DatumFeatureTargetV1 {
    /// Named surface support independent of discretization indices.
    SurfacePatch(SurfacePatchId),
    /// Named contact attachment feature.
    ContactFeature(ContactFeatureId),
}

impl DatumFeatureTargetV1 {
    /// Stable selector role tag used by the aggregate identity.
    #[must_use]
    pub const fn tag(&self) -> u8 {
        match self {
            Self::SurfacePatch(_) => 1,
            Self::ContactFeature(_) => 2,
        }
    }

    /// Exact canonical target key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        match self {
            Self::SurfacePatch(patch) => patch.canonical_key(),
            Self::ContactFeature(feature) => feature.canonical_key(),
        }
    }
}

/// One datum feature bound to a caller-declared body and durable target.
///
/// The declared body is identity-semantic. Admission proves only graph
/// existence and subsystem co-ownership, not geometric containment.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DatumFeatureBindingV1 {
    id: DatumFeatureIdV1,
    declared_body: BodyId,
    target: DatumFeatureTargetV1,
}

impl DatumFeatureBindingV1 {
    /// Construct one authority-free datum-feature declaration.
    #[must_use]
    pub fn new(id: DatumFeatureIdV1, declared_body: BodyId, target: DatumFeatureTargetV1) -> Self {
        Self {
            id,
            declared_body,
            target,
        }
    }

    /// Stable datum-feature identity.
    #[must_use]
    pub const fn id(&self) -> &DatumFeatureIdV1 {
        &self.id
    }

    /// Caller-declared body associated with this feature.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable selected surface or contact feature.
    #[must_use]
    pub const fn target(&self) -> &DatumFeatureTargetV1 {
        &self.target
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(
            self.id.canonical_key().len()
                + self.declared_body.canonical_key().len()
                + self.target.canonical_key().len()
                + 96,
        );
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        append_bytes(&mut row, self.declared_body.identity().as_bytes());
        append_bytes(&mut row, self.declared_body.canonical_key().as_bytes());
        row.push(self.target.tag());
        match &self.target {
            DatumFeatureTargetV1::SurfacePatch(patch) => {
                append_bytes(&mut row, patch.identity().as_bytes());
            }
            DatumFeatureTargetV1::ContactFeature(feature) => {
                append_bytes(&mut row, feature.identity().as_bytes());
            }
        }
        append_bytes(&mut row, self.target.canonical_key().as_bytes());
        row
    }
}

/// Ordered constraint role of one datum reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum DatumPrecedenceV1 {
    /// First reference in the declared frame.
    Primary = 1,
    /// Second reference in the declared frame.
    Secondary = 2,
    /// Third reference in the declared frame.
    Tertiary = 3,
}

impl DatumPrecedenceV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Primary => "primary",
            Self::Secondary => "secondary",
            Self::Tertiary => "tertiary",
        }
    }
}

/// Caller-declared primary/secondary/tertiary datum precedence.
///
/// This is a canonical semantic record, not a geometrically constructed
/// coordinate frame or a 3-2-1 constraint-sufficiency certificate.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DatumReferenceFrameV1 {
    id: DatumReferenceFrameIdV1,
    primary: DatumFeatureIdV1,
    secondary: Option<DatumFeatureIdV1>,
    tertiary: Option<DatumFeatureIdV1>,
}

impl DatumReferenceFrameV1 {
    /// Construct an authority-free frame declaration.
    ///
    /// Invalid gaps, duplicates, and references remain representable so the
    /// aggregate admission boundary can return one structured refusal.
    #[must_use]
    pub fn new(
        id: DatumReferenceFrameIdV1,
        primary: DatumFeatureIdV1,
        secondary: Option<DatumFeatureIdV1>,
        tertiary: Option<DatumFeatureIdV1>,
    ) -> Self {
        Self {
            id,
            primary,
            secondary,
            tertiary,
        }
    }

    /// Stable reference-frame identity.
    #[must_use]
    pub const fn id(&self) -> &DatumReferenceFrameIdV1 {
        &self.id
    }

    /// Required primary datum feature.
    #[must_use]
    pub const fn primary(&self) -> &DatumFeatureIdV1 {
        &self.primary
    }

    /// Optional secondary datum feature.
    #[must_use]
    pub const fn secondary(&self) -> Option<&DatumFeatureIdV1> {
        self.secondary.as_ref()
    }

    /// Optional tertiary datum feature.
    #[must_use]
    pub const fn tertiary(&self) -> Option<&DatumFeatureIdV1> {
        self.tertiary.as_ref()
    }

    /// Number of explicitly declared precedence tiers.
    #[must_use]
    pub const fn reference_count(&self) -> usize {
        1 + self.secondary.is_some() as usize + self.tertiary.is_some() as usize
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(64);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        row.push(DatumPrecedenceV1::Primary.tag());
        append_bytes(&mut row, self.primary.canonical_key().as_bytes());
        match &self.secondary {
            Some(secondary) => {
                row.push(1);
                row.push(DatumPrecedenceV1::Secondary.tag());
                append_bytes(&mut row, secondary.canonical_key().as_bytes());
            }
            None => row.push(0),
        }
        match &self.tertiary {
            Some(tertiary) => {
                row.push(1);
                row.push(DatumPrecedenceV1::Tertiary.tag());
                append_bytes(&mut row, tertiary.canonical_key().as_bytes());
            }
            None => row.push(0),
        }
        row
    }
}

/// Mutable-by-construction graph-bound datum-system draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineDatumSystemDraftV1 {
    /// Datum-feature bindings in non-semantic caller order.
    pub datum_features: Vec<DatumFeatureBindingV1>,
    /// Reference-frame declarations in non-semantic caller order.
    pub reference_frames: Vec<DatumReferenceFrameV1>,
}

impl MachineDatumSystemDraftV1 {
    /// Admit, canonicalize, and bind the declarations to one exact graph.
    ///
    /// # Errors
    /// Refuses resource overflow, duplicate or unknown declarations, ownership
    /// mismatch, malformed precedence, mixed-body frames, unused features, or
    /// bounded identity publication failure.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineDatumSystemV1, MachineDatumAdmissionErrorV1> {
        admit_datum_system(self, graph)
    }
}

/// Structured refusal from graph-bound datum-system admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineDatumAdmissionErrorV1 {
    /// At least one datum feature is required.
    NoDatumFeatures,
    /// Raw datum-feature submissions exceeded the fixed cap.
    DatumFeatureLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Fixed schema cap.
        max: usize,
    },
    /// At least one reference frame is required.
    NoReferenceFrames,
    /// Raw reference-frame submissions exceeded the fixed cap.
    ReferenceFrameLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Fixed schema cap.
        max: usize,
    },
    /// One datum-feature identity appeared more than once.
    DuplicateDatumFeature {
        /// Repeated datum-feature identifier.
        feature: DatumFeatureIdV1,
    },
    /// Two datum IDs tried to alias one typed body/feature selector.
    DuplicateDatumSelector {
        /// Lexically first datum ID selecting the target.
        first: DatumFeatureIdV1,
        /// Later datum ID selecting the same target.
        duplicate: DatumFeatureIdV1,
    },
    /// A selector named a body absent from the admitted graph.
    UnknownBody {
        /// Datum declaration being admitted.
        feature: DatumFeatureIdV1,
        /// Missing caller-declared body.
        body: BodyId,
    },
    /// A selector named a surface/contact feature absent from the graph.
    UnknownFeatureTarget {
        /// Datum declaration being admitted.
        feature: DatumFeatureIdV1,
        /// Missing typed target.
        target: DatumFeatureTargetV1,
    },
    /// Body and feature exist but have different subsystem owners.
    FeatureOwnerMismatch {
        /// Datum declaration being admitted.
        feature: DatumFeatureIdV1,
        /// Caller-declared body.
        body: BodyId,
        /// Selected durable target.
        target: DatumFeatureTargetV1,
        /// Subsystem owning the body.
        body_owner: SubsystemId,
        /// Subsystem owning the selected feature.
        target_owner: SubsystemId,
    },
    /// One reference-frame identity appeared more than once.
    DuplicateReferenceFrame {
        /// Repeated reference-frame identifier.
        frame: DatumReferenceFrameIdV1,
    },
    /// A tertiary datum was supplied without a secondary datum.
    TertiaryWithoutSecondary {
        /// Reference frame supplying a tertiary datum without a secondary datum.
        frame: DatumReferenceFrameIdV1,
    },
    /// A frame referenced an undeclared datum feature.
    MissingDatumReference {
        /// Frame containing the missing reference.
        frame: DatumReferenceFrameIdV1,
        /// Exact precedence tier.
        precedence: DatumPrecedenceV1,
        /// Missing datum identity.
        feature: DatumFeatureIdV1,
    },
    /// The same datum feature occupied more than one tier in a frame.
    RepeatedDatumReference {
        /// Frame containing the duplicate.
        frame: DatumReferenceFrameIdV1,
        /// Repeated feature.
        feature: DatumFeatureIdV1,
        /// First tier using the feature.
        first: DatumPrecedenceV1,
        /// Later tier repeating it.
        repeated: DatumPrecedenceV1,
    },
    /// Version one admits only single-declared-body reference frames.
    MixedBodyReferenceFrame {
        /// Frame containing mixed body declarations.
        frame: DatumReferenceFrameIdV1,
        /// Body declared by the first resolved tier.
        first_body: BodyId,
        /// Conflicting body declared by a later tier.
        conflicting_body: BodyId,
    },
    /// A declared datum feature was not used by any frame.
    UnusedDatumFeature {
        /// Declared datum feature not referenced by any frame.
        feature: DatumFeatureIdV1,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineDatumAdmissionErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoDatumFeatures => "MachineDatumNoFeatures",
            Self::DatumFeatureLimit { .. } => "MachineDatumFeatureLimit",
            Self::NoReferenceFrames => "MachineDatumNoReferenceFrames",
            Self::ReferenceFrameLimit { .. } => "MachineDatumReferenceFrameLimit",
            Self::DuplicateDatumFeature { .. } => "MachineDatumDuplicateFeature",
            Self::DuplicateDatumSelector { .. } => "MachineDatumDuplicateSelector",
            Self::UnknownBody { .. } => "MachineDatumUnknownBody",
            Self::UnknownFeatureTarget { .. } => "MachineDatumUnknownFeatureTarget",
            Self::FeatureOwnerMismatch { .. } => "MachineDatumFeatureOwnerMismatch",
            Self::DuplicateReferenceFrame { .. } => "MachineDatumDuplicateReferenceFrame",
            Self::TertiaryWithoutSecondary { .. } => "MachineDatumTertiaryWithoutSecondary",
            Self::MissingDatumReference { .. } => "MachineDatumMissingReference",
            Self::RepeatedDatumReference { .. } => "MachineDatumRepeatedReference",
            Self::MixedBodyReferenceFrame { .. } => "MachineDatumMixedBodyReferenceFrame",
            Self::UnusedDatumFeature { .. } => "MachineDatumUnusedFeature",
            Self::Identity(_) => "MachineDatumIdentity",
        }
    }
}

impl fmt::Display for MachineDatumAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoDatumFeatures => formatter.write_str("datum system must declare a feature"),
            Self::DatumFeatureLimit { actual, max } => {
                write!(
                    formatter,
                    "datum system has {actual} features; maximum is {max}"
                )
            }
            Self::NoReferenceFrames => {
                formatter.write_str("datum system must declare a reference frame")
            }
            Self::ReferenceFrameLimit { actual, max } => write!(
                formatter,
                "datum system has {actual} reference frames; maximum is {max}"
            ),
            Self::DuplicateDatumFeature { feature } => {
                write!(formatter, "datum feature {feature} appears more than once")
            }
            Self::DuplicateDatumSelector { first, duplicate } => write!(
                formatter,
                "datum features {first} and {duplicate} select the same body feature"
            ),
            Self::UnknownBody { feature, body } => {
                write!(
                    formatter,
                    "datum feature {feature} names unknown body {body}"
                )
            }
            Self::UnknownFeatureTarget { feature, target } => write!(
                formatter,
                "datum feature {feature} names unknown target {}:{}",
                target.tag(),
                target.canonical_key()
            ),
            Self::FeatureOwnerMismatch {
                feature,
                body,
                target,
                body_owner,
                target_owner,
            } => write!(
                formatter,
                "datum feature {feature} declares body {body} owned by {body_owner}, but target \
                 {}:{} is owned by {target_owner}",
                target.tag(),
                target.canonical_key()
            ),
            Self::DuplicateReferenceFrame { frame } => {
                write!(
                    formatter,
                    "datum reference frame {frame} appears more than once"
                )
            }
            Self::TertiaryWithoutSecondary { frame } => write!(
                formatter,
                "datum reference frame {frame} supplies tertiary without secondary"
            ),
            Self::MissingDatumReference {
                frame,
                precedence,
                feature,
            } => write!(
                formatter,
                "datum reference frame {frame} has missing {} feature {feature}",
                precedence.name()
            ),
            Self::RepeatedDatumReference {
                frame,
                feature,
                first,
                repeated,
            } => write!(
                formatter,
                "datum reference frame {frame} repeats feature {feature} at {} after {}",
                repeated.name(),
                first.name()
            ),
            Self::MixedBodyReferenceFrame {
                frame,
                first_body,
                conflicting_body,
            } => write!(
                formatter,
                "datum reference frame {frame} mixes declared bodies {first_body} and \
                 {conflicting_body}"
            ),
            Self::UnusedDatumFeature { feature } => {
                write!(
                    formatter,
                    "datum feature {feature} is not used by any frame"
                )
            }
            Self::Identity(error) => write!(formatter, "datum identity refused: {error}"),
        }
    }
}

impl std::error::Error for MachineDatumAdmissionErrorV1 {}

impl From<CanonicalError> for MachineDatumAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

/// Canonical identity schema for one admitted graph-bound datum catalog.
pub enum MachineDatumSystemIdentitySchemaV1 {}

impl CanonicalSchema for MachineDatumSystemIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.manufacturing-datum-system.v1";
    const NAME: &'static str = "admitted-machine-datum-system";
    const VERSION: u32 = MACHINE_DATUM_SYSTEM_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact Machine graph plus canonical datum-feature selectors and single-body primary-secondary-tertiary reference frames";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("datum-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("datum-features", WireType::OrderedBytes),
        FieldSpec::required("reference-frames", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of an admitted graph-bound datum catalog.
pub type MachineDatumSystemIdV1 = ProblemSemanticId<MachineDatumSystemIdentitySchemaV1>;

/// Canonically ordered graph-bound datum declarations plus identity receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineDatumSystemV1 {
    graph: MachineGraphIdV1,
    datum_features: Vec<DatumFeatureBindingV1>,
    reference_frames: Vec<DatumReferenceFrameV1>,
    receipt: IdentityReceipt<MachineDatumSystemIdV1>,
}

impl AdmittedMachineDatumSystemV1 {
    /// Exact Machine graph extended by these declarations.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Datum features in canonical feature-ID order.
    #[must_use]
    pub fn datum_features(&self) -> &[DatumFeatureBindingV1] {
        &self.datum_features
    }

    /// Reference frames in canonical frame-ID order.
    #[must_use]
    pub fn reference_frames(&self) -> &[DatumReferenceFrameV1] {
        &self.reference_frames
    }

    /// Domain-separated aggregate semantic identity.
    #[must_use]
    pub const fn identity(&self) -> MachineDatumSystemIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineDatumSystemIdV1> {
        self.receipt
    }
}

fn admit_datum_system(
    draft: MachineDatumSystemDraftV1,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineDatumSystemV1, MachineDatumAdmissionErrorV1> {
    if draft.datum_features.is_empty() {
        return Err(MachineDatumAdmissionErrorV1::NoDatumFeatures);
    }
    if draft.datum_features.len() > MAX_MACHINE_DATUM_FEATURES_V1 {
        return Err(MachineDatumAdmissionErrorV1::DatumFeatureLimit {
            actual: draft.datum_features.len(),
            max: MAX_MACHINE_DATUM_FEATURES_V1,
        });
    }
    if draft.reference_frames.is_empty() {
        return Err(MachineDatumAdmissionErrorV1::NoReferenceFrames);
    }
    if draft.reference_frames.len() > MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1 {
        return Err(MachineDatumAdmissionErrorV1::ReferenceFrameLimit {
            actual: draft.reference_frames.len(),
            max: MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1,
        });
    }

    let mut datum_features = draft.datum_features;
    datum_features.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = datum_features
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(MachineDatumAdmissionErrorV1::DuplicateDatumFeature {
            feature: pair[0].id.clone(),
        });
    }

    let mut selected_targets = BTreeMap::<DatumFeatureTargetV1, DatumFeatureIdV1>::new();
    for binding in &datum_features {
        if let Some(first) = selected_targets.insert(binding.target.clone(), binding.id.clone()) {
            return Err(MachineDatumAdmissionErrorV1::DuplicateDatumSelector {
                first,
                duplicate: binding.id.clone(),
            });
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
    let surface_owners = graph
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
    let contact_owners = graph
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

    for binding in &datum_features {
        let Some(body_owner) = body_owners.get(&binding.declared_body) else {
            return Err(MachineDatumAdmissionErrorV1::UnknownBody {
                feature: binding.id.clone(),
                body: binding.declared_body.clone(),
            });
        };
        let target_owner = match &binding.target {
            DatumFeatureTargetV1::SurfacePatch(patch) => surface_owners.get(patch),
            DatumFeatureTargetV1::ContactFeature(feature) => contact_owners.get(feature),
        };
        let Some(target_owner) = target_owner else {
            return Err(MachineDatumAdmissionErrorV1::UnknownFeatureTarget {
                feature: binding.id.clone(),
                target: binding.target.clone(),
            });
        };
        if body_owner != target_owner {
            return Err(MachineDatumAdmissionErrorV1::FeatureOwnerMismatch {
                feature: binding.id.clone(),
                body: binding.declared_body.clone(),
                target: binding.target.clone(),
                body_owner: body_owner.clone(),
                target_owner: target_owner.clone(),
            });
        }
    }

    let mut reference_frames = draft.reference_frames;
    reference_frames.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = reference_frames
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(MachineDatumAdmissionErrorV1::DuplicateReferenceFrame {
            frame: pair[0].id.clone(),
        });
    }

    let feature_by_id = datum_features
        .iter()
        .map(|binding| (binding.id.clone(), binding))
        .collect::<BTreeMap<_, _>>();
    let mut referenced = BTreeSet::<DatumFeatureIdV1>::new();
    for frame in &reference_frames {
        if frame.tertiary.is_some() && frame.secondary.is_none() {
            return Err(MachineDatumAdmissionErrorV1::TertiaryWithoutSecondary {
                frame: frame.id.clone(),
            });
        }

        let references = [
            (DatumPrecedenceV1::Primary, Some(&frame.primary)),
            (DatumPrecedenceV1::Secondary, frame.secondary.as_ref()),
            (DatumPrecedenceV1::Tertiary, frame.tertiary.as_ref()),
        ];
        let mut seen = BTreeMap::<DatumFeatureIdV1, DatumPrecedenceV1>::new();
        let mut frame_body: Option<&BodyId> = None;
        for (precedence, feature) in references {
            let Some(feature) = feature else {
                continue;
            };
            let Some(binding) = feature_by_id.get(feature) else {
                return Err(MachineDatumAdmissionErrorV1::MissingDatumReference {
                    frame: frame.id.clone(),
                    precedence,
                    feature: feature.clone(),
                });
            };
            if let Some(first) = seen.insert(feature.clone(), precedence) {
                return Err(MachineDatumAdmissionErrorV1::RepeatedDatumReference {
                    frame: frame.id.clone(),
                    feature: feature.clone(),
                    first,
                    repeated: precedence,
                });
            }
            if let Some(first_body) = frame_body {
                if first_body != &binding.declared_body {
                    return Err(MachineDatumAdmissionErrorV1::MixedBodyReferenceFrame {
                        frame: frame.id.clone(),
                        first_body: first_body.clone(),
                        conflicting_body: binding.declared_body.clone(),
                    });
                }
            } else {
                frame_body = Some(&binding.declared_body);
            }
            referenced.insert(feature.clone());
        }
    }

    if let Some(feature) = datum_features
        .iter()
        .map(|binding| &binding.id)
        .find(|feature| !referenced.contains(*feature))
    {
        return Err(MachineDatumAdmissionErrorV1::UnusedDatumFeature {
            feature: feature.clone(),
        });
    }

    let feature_rows = datum_features
        .iter()
        .map(DatumFeatureBindingV1::canonical_row)
        .collect::<Vec<_>>();
    let frame_rows = reference_frames
        .iter()
        .map(DatumReferenceFrameV1::canonical_row)
        .collect::<Vec<_>>();
    let graph_id = graph.identity();
    let receipt =
        CanonicalEncoder::<MachineDatumSystemIdV1, _>::new(DATUM_IDENTITY_LIMITS, NeverCancel)?
            .u64(
                Field::new(0, "datum-schema-version"),
                u64::from(MACHINE_DATUM_SYSTEM_SCHEMA_VERSION_V1),
            )?
            .u64(
                Field::new(1, "frankenscript-ir-version"),
                u64::from(IR_VERSION),
            )?
            .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
            .ordered_bytes(
                Field::new(3, "datum-features"),
                feature_rows.len() as u64,
                feature_rows.iter().map(Vec::as_slice),
            )?
            .ordered_bytes(
                Field::new(4, "reference-frames"),
                frame_rows.len() as u64,
                frame_rows.iter().map(Vec::as_slice),
            )?
            .finish()?;

    Ok(AdmittedMachineDatumSystemV1 {
        graph: graph_id,
        datum_features,
        reference_frames,
        receipt,
    })
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
