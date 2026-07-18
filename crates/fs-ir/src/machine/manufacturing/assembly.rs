//! Ordered graph-bound assembly and joining admission.
//!
//! Version one records a total assembly order over one growing assembly rooted
//! at an initial body. `AttachIncoming` introduces a body; `ContinueExisting`
//! adds another joint between present bodies. Operations retain durable
//! body/contact-feature selectors, closed joining-family tags, an explicit
//! preloaded-bolt force target, and exact procedure/path/evidence coordinates.
//! Admission proves only structural graph ownership and order availability. It
//! does not execute a process, validate a path, or establish joint physics.

use core::fmt;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphIdV1, MachineIdError, SubsystemId,
};
use super::ManufacturingArtifactRefV1;

/// Identity/admission schema version for ordered assembly declarations.
pub const MACHINE_ASSEMBLY_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum operations retained by one version-one assembly receipt.
pub const MAX_MACHINE_ASSEMBLY_OPERATIONS_V1: usize = 4_096;

const ASSEMBLY_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(12 * 1_024 * 1_024, 8 * 1_024 * 1_024, 5, 4_096, 4_096);

/// Stable identity of one assembly operation.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AssemblyOperationIdV1(Box<str>);

impl AssemblyOperationIdV1 {
    /// Admit one bounded canonical operation key.
    ///
    /// # Errors
    /// Refuses text outside the Machine-IR canonical key grammar.
    pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
        let key = key.into();
        super::super::validate_canonical_key("assembly-operation-id", &key)?;
        Ok(Self(key.into_boxed_str()))
    }

    /// Exact canonical key retained in aggregate identity rows.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for AssemblyOperationIdV1 {
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
    /// Exact joining procedure/specification coordinate for one operation.
    AssemblyProcedureRefV1
);
artifact_role!(
    /// Exact nominal insertion/approach path coordinate for one operation.
    AssemblyPathRefV1
);
artifact_role!(
    /// Exact nominal as-built execution/evidence coordinate for one operation.
    AssemblyExecutionEvidenceRefV1
);

/// Explicit source unit for one declared preloaded-bolt force target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AssemblyPreloadUnitV1 {
    /// Newtons.
    Newton = 1,
    /// Kilonewtons.
    Kilonewton = 2,
}

impl AssemblyPreloadUnitV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Binary64 multiplier used to normalize to coherent-SI newtons.
    #[must_use]
    pub const fn newtons_per_unit(self) -> f64 {
        match self {
            Self::Newton => 1.0,
            Self::Kilonewton => 1_000.0,
        }
    }

    /// Stable unit spelling.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::Newton => "N",
            Self::Kilonewton => "kN",
        }
    }
}

/// Refusal from constructing a positive preloaded-bolt force target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyPreloadErrorV1 {
    /// NaN or infinity was supplied.
    NonFinite,
    /// Force target was zero or negative.
    NonPositive,
    /// Unit normalization overflowed binary64.
    SiNonFinite,
}

impl AssemblyPreloadErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::NonFinite => "AssemblyPreloadNonFinite",
            Self::NonPositive => "AssemblyPreloadNonPositive",
            Self::SiNonFinite => "AssemblyPreloadSiNonFinite",
        }
    }
}

impl fmt::Display for AssemblyPreloadErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "assembly preload target must be finite",
            Self::NonPositive => "assembly preload target must be positive",
            Self::SiNonFinite => "assembly preload SI normalization must remain finite",
        })
    }
}

impl std::error::Error for AssemblyPreloadErrorV1 {}

/// Strictly positive preload target retaining source and coherent-SI bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AssemblyPreloadV1 {
    submitted_bits: u64,
    unit: AssemblyPreloadUnitV1,
    newtons_bits: u64,
}

impl AssemblyPreloadV1 {
    /// Validate and normalize one positive preload target.
    ///
    /// # Errors
    /// Refuses non-finite, non-positive, or SI-overflowing input.
    pub fn try_new(
        value: f64,
        unit: AssemblyPreloadUnitV1,
    ) -> Result<Self, AssemblyPreloadErrorV1> {
        if !value.is_finite() {
            return Err(AssemblyPreloadErrorV1::NonFinite);
        }
        if value <= 0.0 {
            return Err(AssemblyPreloadErrorV1::NonPositive);
        }
        let newtons = value * unit.newtons_per_unit();
        if !newtons.is_finite() {
            return Err(AssemblyPreloadErrorV1::SiNonFinite);
        }
        Ok(Self {
            submitted_bits: value.to_bits(),
            unit,
            newtons_bits: newtons.to_bits(),
        })
    }

    /// Canonical submitted value.
    #[must_use]
    pub fn submitted_value(self) -> f64 {
        f64::from_bits(self.submitted_bits)
    }

    /// Exact submitted unit.
    #[must_use]
    pub const fn unit(self) -> AssemblyPreloadUnitV1 {
        self.unit
    }

    /// Coherent-SI binary64 force in newtons.
    #[must_use]
    pub fn newtons(self) -> f64 {
        f64::from_bits(self.newtons_bits)
    }

    fn append_canonical(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.submitted_bits.to_le_bytes());
        out.push(self.unit.tag());
        out.extend_from_slice(&self.newtons_bits.to_le_bytes());
    }
}

/// Closed joining-family vocabulary admitted by version one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AssemblyJoinKindV1 {
    /// Bolted joint with one explicit positive preload target.
    PreloadedBolt = 1,
    /// Welded joint declaration.
    Weld = 2,
    /// Adhesive-bonded joint declaration.
    AdhesiveBond = 3,
    /// Keyed joint declaration.
    Key = 4,
    /// Splined joint declaration.
    Spline = 5,
    /// Interference-fit joining declaration.
    InterferenceFit = 6,
}

impl AssemblyJoinKindV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Whether this family structurally requires a preload target.
    #[must_use]
    pub const fn requires_preload(self) -> bool {
        matches!(self, Self::PreloadedBolt)
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::PreloadedBolt => "preloaded-bolt",
            Self::Weld => "weld",
            Self::AdhesiveBond => "adhesive-bond",
            Self::Key => "key",
            Self::Spline => "spline",
            Self::InterferenceFit => "interference-fit",
        }
    }
}

/// Availability transition performed by one ordered operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AssemblyOperationModeV1 {
    /// Introduce an incoming body not yet in the growing assembly.
    AttachIncoming = 1,
    /// Add another joint between two bodies already in the growing assembly.
    ContinueExisting = 2,
}

impl AssemblyOperationModeV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// Ordered endpoint role within one assembly operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AssemblyEndpointRoleV1 {
    /// Feature already present in the partially assembled body set.
    Base,
    /// Second endpoint: introduced by attach mode, already present by continue mode.
    Incoming,
}

/// One caller-declared body/contact-feature selector.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AssemblyFeatureSelectorV1 {
    declared_body: BodyId,
    contact_feature: ContactFeatureId,
}

impl AssemblyFeatureSelectorV1 {
    /// Construct one authority-free assembly endpoint selector.
    #[must_use]
    pub fn new(declared_body: BodyId, contact_feature: ContactFeatureId) -> Self {
        Self {
            declared_body,
            contact_feature,
        }
    }

    /// Caller-declared body; feature containment is not proved by version one.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable contact attachment feature.
    #[must_use]
    pub const fn contact_feature(&self) -> &ContactFeatureId {
        &self.contact_feature
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(336);
        append_bytes(&mut row, self.declared_body.identity().as_bytes());
        append_bytes(&mut row, self.declared_body.canonical_key().as_bytes());
        append_bytes(&mut row, self.contact_feature.identity().as_bytes());
        append_bytes(&mut row, self.contact_feature.canonical_key().as_bytes());
        row
    }
}

/// One explicitly ordered structural assembly operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyOperationV1 {
    id: AssemblyOperationIdV1,
    ordinal: u32,
    mode: AssemblyOperationModeV1,
    base: AssemblyFeatureSelectorV1,
    incoming: AssemblyFeatureSelectorV1,
    join_kind: AssemblyJoinKindV1,
    preload: Option<AssemblyPreloadV1>,
    procedure: AssemblyProcedureRefV1,
    path: AssemblyPathRefV1,
    execution_evidence: AssemblyExecutionEvidenceRefV1,
}

impl AssemblyOperationV1 {
    /// Construct one authority-free assembly operation declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        id: AssemblyOperationIdV1,
        ordinal: u32,
        mode: AssemblyOperationModeV1,
        base: AssemblyFeatureSelectorV1,
        incoming: AssemblyFeatureSelectorV1,
        join_kind: AssemblyJoinKindV1,
        preload: Option<AssemblyPreloadV1>,
        procedure: AssemblyProcedureRefV1,
        path: AssemblyPathRefV1,
        execution_evidence: AssemblyExecutionEvidenceRefV1,
    ) -> Self {
        Self {
            id,
            ordinal,
            mode,
            base,
            incoming,
            join_kind,
            preload,
            procedure,
            path,
            execution_evidence,
        }
    }

    /// Stable operation identity.
    #[must_use]
    pub const fn id(&self) -> &AssemblyOperationIdV1 {
        &self.id
    }

    /// Zero-based total-order position.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Availability transition performed at this ordinal.
    #[must_use]
    pub const fn mode(&self) -> AssemblyOperationModeV1 {
        self.mode
    }

    /// Already-available endpoint.
    #[must_use]
    pub const fn base(&self) -> &AssemblyFeatureSelectorV1 {
        &self.base
    }

    /// Incoming-role endpoint, whose required availability depends on mode.
    #[must_use]
    pub const fn incoming(&self) -> &AssemblyFeatureSelectorV1 {
        &self.incoming
    }

    /// Closed joining-family declaration.
    #[must_use]
    pub const fn join_kind(&self) -> AssemblyJoinKindV1 {
        self.join_kind
    }

    /// Preload required only for the preloaded-bolt family.
    #[must_use]
    pub const fn preload(&self) -> Option<AssemblyPreloadV1> {
        self.preload
    }

    /// Exact procedure/specification coordinate.
    #[must_use]
    pub const fn procedure(&self) -> &AssemblyProcedureRefV1 {
        &self.procedure
    }

    /// Exact nominal approach-path coordinate.
    #[must_use]
    pub const fn path(&self) -> &AssemblyPathRefV1 {
        &self.path
    }

    /// Exact nominal execution/evidence coordinate.
    #[must_use]
    pub const fn execution_evidence(&self) -> &AssemblyExecutionEvidenceRefV1 {
        &self.execution_evidence
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(1_536);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        row.extend_from_slice(&self.ordinal.to_le_bytes());
        row.push(self.mode.tag());
        append_bytes(&mut row, &self.base.canonical_row());
        append_bytes(&mut row, &self.incoming.canonical_row());
        row.push(self.join_kind.tag());
        match self.preload {
            Some(preload) => {
                row.push(1);
                preload.append_canonical(&mut row);
            }
            None => row.push(0),
        }
        append_artifact(&mut row, self.procedure.artifact());
        append_artifact(&mut row, self.path.artifact());
        append_artifact(&mut row, self.execution_evidence.artifact());
        row
    }
}

/// Mutable-by-construction ordered assembly draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssemblyDraftV1 {
    /// Body available before ordinal zero.
    pub initial_body: BodyId,
    /// Operations in non-semantic caller collection order.
    pub operations: Vec<AssemblyOperationV1>,
}

impl MachineAssemblyDraftV1 {
    /// Admit and bind one total structural assembly order to an exact graph.
    ///
    /// # Errors
    /// Refuses resource overflow, invalid order, unavailable bodies, invalid
    /// endpoint ownership or reuse, preload mismatch, and identity failure.
    #[allow(clippy::result_large_err)] // Preserve exact owned IDs in refusals.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineAssemblyV1, MachineAssemblyAdmissionErrorV1> {
        admit_assembly(self, graph)
    }
}

/// Canonical identity schema for one graph-bound assembly order.
pub enum MachineAssemblyIdentitySchemaV1 {}

impl CanonicalSchema for MachineAssemblyIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.manufacturing-assembly.v1";
    const NAME: &'static str = "admitted-machine-assembly";
    const VERSION: u32 = MACHINE_ASSEMBLY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one exact Machine graph, initial body, and canonical ordered structural assembly operations";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("assembly-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("initial-body", WireType::Bytes),
        FieldSpec::required("operations", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted assembly order.
pub type MachineAssemblyIdV1 = ProblemSemanticId<MachineAssemblyIdentitySchemaV1>;

/// Canonically ordered graph-bound assembly operations plus complete receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineAssemblyV1 {
    graph: MachineGraphIdV1,
    initial_body: BodyId,
    operations: Vec<AssemblyOperationV1>,
    receipt: IdentityReceipt<MachineAssemblyIdV1>,
}

impl AdmittedMachineAssemblyV1 {
    /// Exact Machine graph extended by this assembly record.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Body available before ordinal zero.
    #[must_use]
    pub const fn initial_body(&self) -> &BodyId {
        &self.initial_body
    }

    /// Operations in checked ordinal order.
    #[must_use]
    pub fn operations(&self) -> &[AssemblyOperationV1] {
        &self.operations
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineAssemblyIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineAssemblyIdV1> {
        self.receipt
    }
}

/// Structural mismatch between a join family and its preload field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyPreloadUseIssueV1 {
    /// A preloaded-bolt operation omitted its target.
    PreloadedBoltMissing,
    /// A non-bolted operation supplied a preload target.
    NonBoltHasPreload,
}

/// Structured refusal from assembly admission.
#[allow(clippy::large_enum_variant)] // Preserve exact owned IDs in rich refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineAssemblyAdmissionErrorV1 {
    /// At least one assembly operation is required.
    NoOperations,
    /// Raw operations exceeded the fixed cap.
    OperationLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// Initial body is absent from the admitted graph.
    UnknownInitialBody {
        /// Missing initial body.
        body: BodyId,
    },
    /// One operation identity appeared more than once.
    DuplicateOperation {
        /// Repeated operation identity.
        operation: AssemblyOperationIdV1,
    },
    /// Two operations declared one ordinal.
    DuplicateOrdinal {
        /// Repeated ordinal.
        ordinal: u32,
        /// Lexically first operation at that ordinal.
        first: AssemblyOperationIdV1,
        /// Later operation at that ordinal.
        duplicate: AssemblyOperationIdV1,
    },
    /// Sorted ordinals were not exactly zero through N minus one.
    OrdinalGap {
        /// Operation exposing the gap.
        operation: AssemblyOperationIdV1,
        /// Required zero-based ordinal.
        expected: u32,
        /// Submitted ordinal.
        actual: u32,
    },
    /// An endpoint named a body absent from the graph.
    UnknownBody {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Endpoint role.
        role: AssemblyEndpointRoleV1,
        /// Missing body.
        body: BodyId,
    },
    /// An endpoint named a contact feature absent from the graph.
    UnknownFeature {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Endpoint role.
        role: AssemblyEndpointRoleV1,
        /// Missing contact feature.
        feature: ContactFeatureId,
    },
    /// Endpoint body and feature exist under different subsystem owners.
    FeatureOwnerMismatch {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Endpoint role.
        role: AssemblyEndpointRoleV1,
        /// Caller-declared body.
        body: BodyId,
        /// Selected contact feature.
        feature: ContactFeatureId,
        /// Graph owner of the body.
        body_owner: SubsystemId,
        /// Graph owner of the feature.
        feature_owner: SubsystemId,
    },
    /// Both endpoint roles selected one durable contact feature.
    SameFeature {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Reused feature.
        feature: ContactFeatureId,
    },
    /// Both endpoint roles declared one body.
    SameBody {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Self-joined body.
        body: BodyId,
    },
    /// Join family and preload presence disagree.
    InvalidPreloadUse {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Exact mismatch.
        issue: AssemblyPreloadUseIssueV1,
    },
    /// Base body was not yet present at this ordinal.
    BaseUnavailable {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Unavailable base body.
        body: BodyId,
    },
    /// Incoming body had already entered the assembly.
    IncomingAlreadyAttached {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Reintroduced body.
        body: BodyId,
    },
    /// A continuation named an incoming body not yet in the assembly.
    ContinuationIncomingUnavailable {
        /// Invalid operation.
        operation: AssemblyOperationIdV1,
        /// Unavailable incoming body.
        body: BodyId,
    },
    /// One durable contact feature was consumed by multiple operations.
    FeatureReuse {
        /// Reused contact feature.
        feature: ContactFeatureId,
        /// First operation consuming it.
        first: AssemblyOperationIdV1,
        /// Later operation consuming it.
        duplicate: AssemblyOperationIdV1,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineAssemblyAdmissionErrorV1 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoOperations => "MachineAssemblyNoOperations",
            Self::OperationLimit { .. } => "MachineAssemblyOperationLimit",
            Self::UnknownInitialBody { .. } => "MachineAssemblyUnknownInitialBody",
            Self::DuplicateOperation { .. } => "MachineAssemblyDuplicateOperation",
            Self::DuplicateOrdinal { .. } => "MachineAssemblyDuplicateOrdinal",
            Self::OrdinalGap { .. } => "MachineAssemblyOrdinalGap",
            Self::UnknownBody { .. } => "MachineAssemblyUnknownBody",
            Self::UnknownFeature { .. } => "MachineAssemblyUnknownFeature",
            Self::FeatureOwnerMismatch { .. } => "MachineAssemblyFeatureOwnerMismatch",
            Self::SameFeature { .. } => "MachineAssemblySameFeature",
            Self::SameBody { .. } => "MachineAssemblySameBody",
            Self::InvalidPreloadUse { .. } => "MachineAssemblyInvalidPreloadUse",
            Self::BaseUnavailable { .. } => "MachineAssemblyBaseUnavailable",
            Self::IncomingAlreadyAttached { .. } => "MachineAssemblyIncomingAlreadyAttached",
            Self::ContinuationIncomingUnavailable { .. } => {
                "MachineAssemblyContinuationIncomingUnavailable"
            }
            Self::FeatureReuse { .. } => "MachineAssemblyFeatureReuse",
            Self::Identity(_) => "MachineAssemblyIdentity",
        }
    }
}

impl fmt::Display for MachineAssemblyAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoOperations => formatter.write_str("assembly requires an operation"),
            Self::OperationLimit { actual, max } => {
                write!(
                    formatter,
                    "assembly has {actual} operations; maximum is {max}"
                )
            }
            Self::UnknownInitialBody { body } => {
                write!(
                    formatter,
                    "assembly initial body {body} is absent from the graph"
                )
            }
            Self::DuplicateOperation { operation } => {
                write!(formatter, "assembly operation {operation} is repeated")
            }
            Self::DuplicateOrdinal {
                ordinal,
                first,
                duplicate,
            } => write!(
                formatter,
                "assembly ordinal {ordinal} is shared by {first} and {duplicate}"
            ),
            Self::OrdinalGap {
                operation,
                expected,
                actual,
            } => write!(
                formatter,
                "assembly operation {operation} has ordinal {actual}; expected {expected}"
            ),
            Self::UnknownBody {
                operation,
                role,
                body,
            } => write!(
                formatter,
                "assembly operation {operation} {role:?} endpoint names unknown body {body}"
            ),
            Self::UnknownFeature {
                operation,
                role,
                feature,
            } => write!(
                formatter,
                "assembly operation {operation} {role:?} endpoint names unknown feature {feature}"
            ),
            Self::FeatureOwnerMismatch {
                operation,
                role,
                body,
                feature,
                body_owner,
                feature_owner,
            } => write!(
                formatter,
                "assembly operation {operation} {role:?} body {body} is owned by {body_owner}, but feature {feature} is owned by {feature_owner}"
            ),
            Self::SameFeature { operation, feature } => write!(
                formatter,
                "assembly operation {operation} uses feature {feature} in both endpoint roles"
            ),
            Self::SameBody { operation, body } => write!(
                formatter,
                "assembly operation {operation} uses body {body} in both endpoint roles"
            ),
            Self::InvalidPreloadUse { operation, issue } => {
                write!(
                    formatter,
                    "assembly operation {operation} has invalid preload use {issue:?}"
                )
            }
            Self::BaseUnavailable { operation, body } => write!(
                formatter,
                "assembly operation {operation} base body {body} is not yet available"
            ),
            Self::IncomingAlreadyAttached { operation, body } => write!(
                formatter,
                "assembly operation {operation} incoming body {body} is already attached"
            ),
            Self::ContinuationIncomingUnavailable { operation, body } => write!(
                formatter,
                "assembly continuation {operation} incoming body {body} is not yet available"
            ),
            Self::FeatureReuse {
                feature,
                first,
                duplicate,
            } => write!(
                formatter,
                "assembly feature {feature} is reused by {duplicate} after {first}"
            ),
            Self::Identity(error) => write!(formatter, "assembly identity refused: {error}"),
        }
    }
}

impl std::error::Error for MachineAssemblyAdmissionErrorV1 {}

impl From<CanonicalError> for MachineAssemblyAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::result_large_err)] // Preserve exact owned IDs in rich refusals.
fn admit_assembly(
    draft: MachineAssemblyDraftV1,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineAssemblyV1, MachineAssemblyAdmissionErrorV1> {
    if draft.operations.is_empty() {
        return Err(MachineAssemblyAdmissionErrorV1::NoOperations);
    }
    if draft.operations.len() > MAX_MACHINE_ASSEMBLY_OPERATIONS_V1 {
        return Err(MachineAssemblyAdmissionErrorV1::OperationLimit {
            actual: draft.operations.len(),
            max: MAX_MACHINE_ASSEMBLY_OPERATIONS_V1,
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

    if !body_owners.contains_key(&draft.initial_body) {
        return Err(MachineAssemblyAdmissionErrorV1::UnknownInitialBody {
            body: draft.initial_body,
        });
    }

    let mut operations = draft.operations;
    operations.sort_by(|left, right| {
        left.ordinal
            .cmp(&right.ordinal)
            .then_with(|| left.id.cmp(&right.id))
    });

    let mut operation_ids = BTreeSet::<AssemblyOperationIdV1>::new();
    for operation in &operations {
        if !operation_ids.insert(operation.id.clone()) {
            return Err(MachineAssemblyAdmissionErrorV1::DuplicateOperation {
                operation: operation.id.clone(),
            });
        }
    }
    if let Some(pair) = operations
        .windows(2)
        .find(|pair| pair[0].ordinal == pair[1].ordinal)
    {
        return Err(MachineAssemblyAdmissionErrorV1::DuplicateOrdinal {
            ordinal: pair[0].ordinal,
            first: pair[0].id.clone(),
            duplicate: pair[1].id.clone(),
        });
    }
    for (index, operation) in operations.iter().enumerate() {
        let expected =
            u32::try_from(index).map_err(|_| MachineAssemblyAdmissionErrorV1::OperationLimit {
                actual: operations.len(),
                max: MAX_MACHINE_ASSEMBLY_OPERATIONS_V1,
            })?;
        if operation.ordinal != expected {
            return Err(MachineAssemblyAdmissionErrorV1::OrdinalGap {
                operation: operation.id.clone(),
                expected,
                actual: operation.ordinal,
            });
        }
    }

    let mut available_bodies = BTreeSet::<BodyId>::from([draft.initial_body.clone()]);
    let mut used_features = BTreeMap::<ContactFeatureId, AssemblyOperationIdV1>::new();
    for operation in &operations {
        for (role, endpoint) in [
            (AssemblyEndpointRoleV1::Base, &operation.base),
            (AssemblyEndpointRoleV1::Incoming, &operation.incoming),
        ] {
            let Some(body_owner) = body_owners.get(&endpoint.declared_body) else {
                return Err(MachineAssemblyAdmissionErrorV1::UnknownBody {
                    operation: operation.id.clone(),
                    role,
                    body: endpoint.declared_body.clone(),
                });
            };
            let Some(feature_owner) = feature_owners.get(&endpoint.contact_feature) else {
                return Err(MachineAssemblyAdmissionErrorV1::UnknownFeature {
                    operation: operation.id.clone(),
                    role,
                    feature: endpoint.contact_feature.clone(),
                });
            };
            if body_owner != feature_owner {
                return Err(MachineAssemblyAdmissionErrorV1::FeatureOwnerMismatch {
                    operation: operation.id.clone(),
                    role,
                    body: endpoint.declared_body.clone(),
                    feature: endpoint.contact_feature.clone(),
                    body_owner: body_owner.clone(),
                    feature_owner: feature_owner.clone(),
                });
            }
        }

        if operation.base.contact_feature == operation.incoming.contact_feature {
            return Err(MachineAssemblyAdmissionErrorV1::SameFeature {
                operation: operation.id.clone(),
                feature: operation.base.contact_feature.clone(),
            });
        }
        if operation.base.declared_body == operation.incoming.declared_body {
            return Err(MachineAssemblyAdmissionErrorV1::SameBody {
                operation: operation.id.clone(),
                body: operation.base.declared_body.clone(),
            });
        }

        match (operation.join_kind.requires_preload(), operation.preload) {
            (true, None) => {
                return Err(MachineAssemblyAdmissionErrorV1::InvalidPreloadUse {
                    operation: operation.id.clone(),
                    issue: AssemblyPreloadUseIssueV1::PreloadedBoltMissing,
                });
            }
            (false, Some(_)) => {
                return Err(MachineAssemblyAdmissionErrorV1::InvalidPreloadUse {
                    operation: operation.id.clone(),
                    issue: AssemblyPreloadUseIssueV1::NonBoltHasPreload,
                });
            }
            (true, Some(_)) | (false, None) => {}
        }

        if !available_bodies.contains(&operation.base.declared_body) {
            return Err(MachineAssemblyAdmissionErrorV1::BaseUnavailable {
                operation: operation.id.clone(),
                body: operation.base.declared_body.clone(),
            });
        }
        match operation.mode {
            AssemblyOperationModeV1::AttachIncoming => {
                if available_bodies.contains(&operation.incoming.declared_body) {
                    return Err(MachineAssemblyAdmissionErrorV1::IncomingAlreadyAttached {
                        operation: operation.id.clone(),
                        body: operation.incoming.declared_body.clone(),
                    });
                }
            }
            AssemblyOperationModeV1::ContinueExisting => {
                if !available_bodies.contains(&operation.incoming.declared_body) {
                    return Err(
                        MachineAssemblyAdmissionErrorV1::ContinuationIncomingUnavailable {
                            operation: operation.id.clone(),
                            body: operation.incoming.declared_body.clone(),
                        },
                    );
                }
            }
        }

        for feature in [
            &operation.base.contact_feature,
            &operation.incoming.contact_feature,
        ] {
            if let Some(first) = used_features.insert(feature.clone(), operation.id.clone()) {
                return Err(MachineAssemblyAdmissionErrorV1::FeatureReuse {
                    feature: feature.clone(),
                    first,
                    duplicate: operation.id.clone(),
                });
            }
        }
        if operation.mode == AssemblyOperationModeV1::AttachIncoming {
            available_bodies.insert(operation.incoming.declared_body.clone());
        }
    }

    let initial_body = draft.initial_body;
    let mut initial_body_row = Vec::with_capacity(176);
    append_bytes(&mut initial_body_row, initial_body.identity().as_bytes());
    append_bytes(
        &mut initial_body_row,
        initial_body.canonical_key().as_bytes(),
    );
    let rows = operations
        .iter()
        .map(AssemblyOperationV1::canonical_row)
        .collect::<Vec<_>>();
    let graph_id = graph.identity();
    let receipt =
        CanonicalEncoder::<MachineAssemblyIdV1, _>::new(ASSEMBLY_IDENTITY_LIMITS, NeverCancel)?
            .u64(
                Field::new(0, "assembly-schema-version"),
                u64::from(MACHINE_ASSEMBLY_SCHEMA_VERSION_V1),
            )?
            .u64(
                Field::new(1, "frankenscript-ir-version"),
                u64::from(IR_VERSION),
            )?
            .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
            .bytes(Field::new(3, "initial-body"), &initial_body_row)?
            .ordered_bytes(
                Field::new(4, "operations"),
                rows.len() as u64,
                rows.iter().map(Vec::as_slice),
            )?
            .finish()?;

    Ok(AdmittedMachineAssemblyV1 {
        graph: graph_id,
        initial_body,
        operations,
        receipt,
    })
}

fn append_artifact(out: &mut Vec<u8>, artifact: &ManufacturingArtifactRefV1) {
    append_bytes(out, &artifact.canonical_row());
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
