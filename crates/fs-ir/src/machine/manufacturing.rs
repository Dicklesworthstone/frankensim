//! Graph-bound as-built manufacturing process lineage.
//!
//! This additive E7 seed records a bounded, linear process history for each
//! durable [`BodyId`] in one admitted Machine graph. Every step binds the
//! process specification and the resulting microstructure, residual-stress,
//! and property-state artifact coordinates, together with the exact external
//! correlation model used for process tolerances. Caller order is not
//! semantic: admission reconstructs each body's predecessor chain before
//! publishing a receipt.
//!
//! The receipt is structural. It does not prove that a process ran, that an
//! artifact hash names authentic measurements, that a correlation model is
//! valid, or that the derived material state is physically correct. Process
//! steps are exposed as manufacturing lineage dependents so the existing
//! Machine-IR split/remesh/wear law rebinds them only across an unambiguous
//! durable-body morphism and invalidates them on an ambiguous split.

use core::fmt;
use core::num::NonZeroU64;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;
use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::{
    AdmittedMachineGraph, BodyId, DependentBinding, DependentKind, MachineGraphIdV1, MachineIdError,
};

/// Ordered, graph-bound assembly/joining declarations.
pub mod assembly;
/// Bounded graph-bound datum features and ordered reference frames.
pub mod datum_system;
/// Typed mating-feature size envelopes and derived clearance state.
pub mod fit_clearance;
/// Explicit applicability links between fit endpoints and GD&T controls.
pub mod fit_gdt_crosswalk;
/// Datum-backed flatness, parallelism, and perpendicularity controls.
pub mod geometric_tolerance;
/// Typed surface-texture requirements and as-built observations.
pub mod surface_texture;
/// Exact positional binding from Machine-IR tolerance axes to an admitted
/// `fs-toleralloc` correlated-stack receipt.
pub mod tolerance_axis;

/// Identity-schema version for graph-bound manufacturing state.
pub const MACHINE_MANUFACTURING_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum process steps retained by one manufacturing-state receipt.
pub const MAX_MANUFACTURING_PROCESS_STEPS_V1: usize = 4_096;

const MANUFACTURING_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(16 * 1_024 * 1_024, 8 * 1_024 * 1_024, 5, 8_192, 4_096);

/// Refusal from constructing an exact manufacturing artifact coordinate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ManufacturingReferenceErrorV1 {
    /// The namespace violates the bounded Machine-IR key grammar.
    Namespace(MachineIdError),
    /// An all-zero digest cannot identify external artifact content.
    ZeroDigest,
}

impl ManufacturingReferenceErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Namespace(_) => "ManufacturingReferenceNamespace",
            Self::ZeroDigest => "ManufacturingReferenceZeroDigest",
        }
    }
}

impl fmt::Display for ManufacturingReferenceErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Namespace(error) => write!(formatter, "invalid manufacturing reference: {error}"),
            Self::ZeroDigest => {
                formatter.write_str("manufacturing artifact digest must not be all zero")
            }
        }
    }
}

impl std::error::Error for ManufacturingReferenceErrorV1 {}

/// Bounded, versioned, content-addressed manufacturing artifact coordinate.
///
/// The coordinate is nominal. This module binds it exactly but does not
/// authenticate its owner, bytes, schema implementation, or scientific claim.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ManufacturingArtifactRefV1 {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    content_hash: ContentHash,
}

impl ManufacturingArtifactRefV1 {
    /// Construct one exact external artifact coordinate.
    ///
    /// # Errors
    /// Refuses a noncanonical namespace or all-zero digest.
    pub fn new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        content_hash: ContentHash,
    ) -> Result<Self, ManufacturingReferenceErrorV1> {
        let namespace = namespace.into();
        super::validate_canonical_key("manufacturing-artifact-ref", &namespace)
            .map_err(ManufacturingReferenceErrorV1::Namespace)?;
        if content_hash.as_bytes() == &[0; 32] {
            return Err(ManufacturingReferenceErrorV1::ZeroDigest);
        }
        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            content_hash,
        })
    }

    /// External schema namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Explicit nonzero external schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact caller-supplied content hash.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(self.namespace.len() + 48);
        append_bytes(&mut row, self.namespace.as_bytes());
        row.extend_from_slice(&self.schema_version.get().to_le_bytes());
        row.extend_from_slice(self.content_hash.as_bytes());
        row
    }
}

/// Stable, human-auditable identity of one manufacturing process step.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ManufacturingStepIdV1(Box<str>);

impl ManufacturingStepIdV1 {
    /// Admit a canonical hierarchical step key.
    ///
    /// # Errors
    /// Refuses a key outside the Machine-IR bounded key grammar.
    pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
        let key = key.into();
        super::validate_canonical_key("manufacturing-step-id", &key)?;
        Ok(Self(key.into_boxed_str()))
    }

    /// Exact canonical key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for ManufacturingStepIdV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.canonical_key())
    }
}

/// Declared process family for one as-built lineage transition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ManufacturingProcessKindV1 {
    /// Casting from feedstock into a body occurrence.
    Casting = 1,
    /// Forging or forming operation.
    Forging = 2,
    /// Material-removal machining operation.
    Machining = 3,
    /// Additive manufacturing operation.
    AdditiveManufacturing = 4,
    /// Heat-treatment operation.
    HeatTreatment = 5,
    /// Coating or surface-treatment operation.
    Coating = 6,
    /// Assembly/joining operation represented at body-state granularity.
    Assembly = 7,
}

impl ManufacturingProcessKindV1 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u64 {
        self as u64
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Casting => "casting",
            Self::Forging => "forging",
            Self::Machining => "machining",
            Self::AdditiveManufacturing => "additive-manufacturing",
            Self::HeatTreatment => "heat-treatment",
            Self::Coating => "coating",
            Self::Assembly => "assembly",
        }
    }
}

/// One explicit process-to-material-state transition for a durable body.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ManufacturingProcessStepV1 {
    id: ManufacturingStepIdV1,
    body: BodyId,
    predecessor: Option<ManufacturingStepIdV1>,
    process: ManufacturingProcessKindV1,
    process_specification: ManufacturingArtifactRefV1,
    input_material_state: ManufacturingArtifactRefV1,
    microstructure_state: ManufacturingArtifactRefV1,
    residual_stress_state: ManufacturingArtifactRefV1,
    property_state: ManufacturingArtifactRefV1,
}

impl ManufacturingProcessStepV1 {
    /// Construct one authority-free process-lineage declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        id: ManufacturingStepIdV1,
        body: BodyId,
        predecessor: Option<ManufacturingStepIdV1>,
        process: ManufacturingProcessKindV1,
        process_specification: ManufacturingArtifactRefV1,
        input_material_state: ManufacturingArtifactRefV1,
        microstructure_state: ManufacturingArtifactRefV1,
        residual_stress_state: ManufacturingArtifactRefV1,
        property_state: ManufacturingArtifactRefV1,
    ) -> Self {
        Self {
            id,
            body,
            predecessor,
            process,
            process_specification,
            input_material_state,
            microstructure_state,
            residual_stress_state,
            property_state,
        }
    }

    /// Stable process-step key.
    #[must_use]
    pub const fn id(&self) -> &ManufacturingStepIdV1 {
        &self.id
    }

    /// Durable body whose as-built state this step advances.
    #[must_use]
    pub const fn body(&self) -> &BodyId {
        &self.body
    }

    /// Previous step in the same body's linear history, or `None` for its root.
    #[must_use]
    pub const fn predecessor(&self) -> Option<&ManufacturingStepIdV1> {
        self.predecessor.as_ref()
    }

    /// Declared process family.
    #[must_use]
    pub const fn process(&self) -> ManufacturingProcessKindV1 {
        self.process
    }

    /// Exact process-specification artifact coordinate.
    #[must_use]
    pub const fn process_specification(&self) -> &ManufacturingArtifactRefV1 {
        &self.process_specification
    }

    /// Exact material-state artifact consumed by the process.
    #[must_use]
    pub const fn input_material_state(&self) -> &ManufacturingArtifactRefV1 {
        &self.input_material_state
    }

    /// Exact resulting microstructure artifact coordinate.
    #[must_use]
    pub const fn microstructure_state(&self) -> &ManufacturingArtifactRefV1 {
        &self.microstructure_state
    }

    /// Exact resulting residual-stress artifact coordinate.
    #[must_use]
    pub const fn residual_stress_state(&self) -> &ManufacturingArtifactRefV1 {
        &self.residual_stress_state
    }

    /// Exact resulting property-state artifact coordinate.
    #[must_use]
    pub const fn property_state(&self) -> &ManufacturingArtifactRefV1 {
        &self.property_state
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(512);
        append_bytes(&mut row, self.id.canonical_key().as_bytes());
        row.extend_from_slice(self.body.identity().as_bytes());
        match &self.predecessor {
            Some(predecessor) => {
                row.push(1);
                append_bytes(&mut row, predecessor.canonical_key().as_bytes());
            }
            None => row.push(0),
        }
        row.extend_from_slice(&self.process.tag().to_le_bytes());
        for artifact in [
            &self.process_specification,
            &self.input_material_state,
            &self.microstructure_state,
            &self.residual_stress_state,
            &self.property_state,
        ] {
            append_bytes(&mut row, &artifact.canonical_row());
        }
        row
    }
}

/// Mutable-by-construction, authority-free manufacturing-state draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineManufacturingDraftV1 {
    /// Exact external correlation model used by declared process tolerances.
    pub correlation_model: ManufacturingArtifactRefV1,
    /// Process transitions in arbitrary caller order.
    pub process_steps: Vec<ManufacturingProcessStepV1>,
}

impl MachineManufacturingDraftV1 {
    /// Admit body-bound linear process histories against one exact graph.
    ///
    /// # Errors
    /// Refuses empty/oversized histories, unknown bodies, duplicate or dangling
    /// step IDs, cross-body predecessors, forks, cycles, disconnected chains,
    /// or canonical identity publication failure.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineManufacturingStateV1, ManufacturingAdmissionErrorV1> {
        admit_manufacturing_state(self, graph)
    }
}

/// Structured refusal from manufacturing-state admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ManufacturingAdmissionErrorV1 {
    /// A manufacturing-state receipt must contain at least one process step.
    NoProcessSteps,
    /// The public process-step cap was exceeded.
    ProcessStepLimit {
        /// Supplied step count.
        actual: usize,
        /// Versioned maximum.
        max: usize,
    },
    /// One process step names a body absent from the admitted graph.
    UnknownBody {
        /// Offending step.
        step: ManufacturingStepIdV1,
        /// Missing body.
        body: BodyId,
    },
    /// A stable step ID appeared more than once.
    DuplicateStep {
        /// Duplicate step ID.
        step: ManufacturingStepIdV1,
    },
    /// One predecessor key is absent from the submitted state.
    MissingPredecessor {
        /// Dependent step.
        step: ManufacturingStepIdV1,
        /// Missing predecessor.
        predecessor: ManufacturingStepIdV1,
    },
    /// A process chain attempted to cross durable body identities.
    CrossBodyPredecessor {
        /// Dependent step.
        step: ManufacturingStepIdV1,
        /// Predecessor on a different body.
        predecessor: ManufacturingStepIdV1,
    },
    /// Two steps claim the same immediate predecessor.
    ForkedProcessChain {
        /// Fork point.
        predecessor: ManufacturingStepIdV1,
        /// First successor in canonical step-ID order.
        first: ManufacturingStepIdV1,
        /// Second successor in canonical step-ID order.
        second: ManufacturingStepIdV1,
    },
    /// A body history has zero or multiple roots.
    RootCardinality {
        /// Body whose history is malformed.
        body: BodyId,
        /// Observed root count.
        roots: usize,
    },
    /// A body history contains an unreachable cycle or disconnected segment.
    DisconnectedProcessChain {
        /// Body whose history is malformed.
        body: BodyId,
    },
    /// Canonical identity publication failed.
    Identity(CanonicalError),
}

impl ManufacturingAdmissionErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoProcessSteps => "ManufacturingNoProcessSteps",
            Self::ProcessStepLimit { .. } => "ManufacturingProcessStepLimit",
            Self::UnknownBody { .. } => "ManufacturingUnknownBody",
            Self::DuplicateStep { .. } => "ManufacturingDuplicateStep",
            Self::MissingPredecessor { .. } => "ManufacturingMissingPredecessor",
            Self::CrossBodyPredecessor { .. } => "ManufacturingCrossBodyPredecessor",
            Self::ForkedProcessChain { .. } => "ManufacturingForkedProcessChain",
            Self::RootCardinality { .. } => "ManufacturingRootCardinality",
            Self::DisconnectedProcessChain { .. } => "ManufacturingDisconnectedProcessChain",
            Self::Identity(_) => "ManufacturingIdentity",
        }
    }
}

impl fmt::Display for ManufacturingAdmissionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoProcessSteps => {
                formatter.write_str("manufacturing state must contain at least one process step")
            }
            Self::ProcessStepLimit { actual, max } => write!(
                formatter,
                "manufacturing state has {actual} process steps; maximum is {max}"
            ),
            Self::UnknownBody { step, body } => {
                write!(
                    formatter,
                    "manufacturing step {step} names unknown body {body}"
                )
            }
            Self::DuplicateStep { step } => {
                write!(
                    formatter,
                    "manufacturing step {step} appears more than once"
                )
            }
            Self::MissingPredecessor { step, predecessor } => write!(
                formatter,
                "manufacturing step {step} names missing predecessor {predecessor}"
            ),
            Self::CrossBodyPredecessor { step, predecessor } => write!(
                formatter,
                "manufacturing step {step} and predecessor {predecessor} belong to different bodies"
            ),
            Self::ForkedProcessChain {
                predecessor,
                first,
                second,
            } => write!(
                formatter,
                "manufacturing predecessor {predecessor} forks to {first} and {second}"
            ),
            Self::RootCardinality { body, roots } => write!(
                formatter,
                "manufacturing history for {body} has {roots} roots; exactly one is required"
            ),
            Self::DisconnectedProcessChain { body } => write!(
                formatter,
                "manufacturing history for {body} is cyclic or disconnected"
            ),
            Self::Identity(error) => write!(formatter, "manufacturing identity refused: {error}"),
        }
    }
}

impl std::error::Error for ManufacturingAdmissionErrorV1 {}

impl From<CanonicalError> for ManufacturingAdmissionErrorV1 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

/// Canonical identity schema for one graph-bound manufacturing state.
pub enum MachineManufacturingIdentitySchemaV1 {}

impl CanonicalSchema for MachineManufacturingIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.manufacturing-state.v1";
    const NAME: &'static str = "admitted-machine-manufacturing-state";
    const VERSION: u32 = MACHINE_MANUFACTURING_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one admitted Machine graph, one exact external process-tolerance correlation model, and canonical per-body process-to-material-state lineages";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("manufacturing-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("correlation-model", WireType::Bytes),
        FieldSpec::required("process-steps", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of admitted as-built manufacturing state.
pub type MachineManufacturingStateIdV1 = ProblemSemanticId<MachineManufacturingIdentitySchemaV1>;

/// Canonically ordered, graph-bound manufacturing state plus semantic receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineManufacturingStateV1 {
    graph: MachineGraphIdV1,
    correlation_model: ManufacturingArtifactRefV1,
    process_steps: Vec<ManufacturingProcessStepV1>,
    receipt: IdentityReceipt<MachineManufacturingStateIdV1>,
}

impl AdmittedMachineManufacturingStateV1 {
    /// Exact admitted Machine graph extended by this state.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact external process-tolerance correlation-model coordinate.
    #[must_use]
    pub const fn correlation_model(&self) -> &ManufacturingArtifactRefV1 {
        &self.correlation_model
    }

    /// Process steps in canonical body-chain order.
    #[must_use]
    pub fn process_steps(&self) -> &[ManufacturingProcessStepV1] {
        &self.process_steps
    }

    /// Domain-separated semantic identity.
    #[must_use]
    pub const fn identity(&self) -> MachineManufacturingStateIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineManufacturingStateIdV1> {
        self.receipt
    }

    /// Existing Machine-lineage dependents for every retained process step.
    ///
    /// Supplying these to [`super::LineageRecord::admit`] moves the complete
    /// history across a one-to-one body morphism and invalidates it across an
    /// ambiguous split/remesh rather than guessing a descendant.
    #[must_use]
    pub fn lineage_dependents(&self) -> Vec<DependentBinding> {
        self.process_steps
            .iter()
            .map(|step| {
                DependentBinding::new(
                    DependentKind::ManufacturingState,
                    step.id.canonical_key(),
                    step.body.clone().into(),
                )
                .expect("admitted manufacturing step IDs remain canonical")
            })
            .collect()
    }
}

fn admit_manufacturing_state(
    draft: MachineManufacturingDraftV1,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineManufacturingStateV1, ManufacturingAdmissionErrorV1> {
    if draft.process_steps.is_empty() {
        return Err(ManufacturingAdmissionErrorV1::NoProcessSteps);
    }
    if draft.process_steps.len() > MAX_MANUFACTURING_PROCESS_STEPS_V1 {
        return Err(ManufacturingAdmissionErrorV1::ProcessStepLimit {
            actual: draft.process_steps.len(),
            max: MAX_MANUFACTURING_PROCESS_STEPS_V1,
        });
    }

    let known_bodies: BTreeSet<BodyId> = graph
        .subsystems()
        .iter()
        .flat_map(|subsystem| subsystem.bodies.iter().cloned())
        .collect();

    let mut seen_ids = BTreeSet::<ManufacturingStepIdV1>::new();
    let mut duplicate_ids = BTreeSet::<ManufacturingStepIdV1>::new();
    let mut unknown_bodies = BTreeSet::<(ManufacturingStepIdV1, BodyId)>::new();
    for step in &draft.process_steps {
        if !seen_ids.insert(step.id.clone()) {
            duplicate_ids.insert(step.id.clone());
        }
        if !known_bodies.contains(&step.body) {
            unknown_bodies.insert((step.id.clone(), step.body.clone()));
        }
    }
    if let Some(step) = duplicate_ids.into_iter().next() {
        return Err(ManufacturingAdmissionErrorV1::DuplicateStep { step });
    }
    if let Some((step, body)) = unknown_bodies.into_iter().next() {
        return Err(ManufacturingAdmissionErrorV1::UnknownBody { step, body });
    }

    let mut by_id = BTreeMap::<ManufacturingStepIdV1, ManufacturingProcessStepV1>::new();
    for step in draft.process_steps {
        let id = step.id.clone();
        let previous = by_id.insert(id, step);
        debug_assert!(previous.is_none(), "duplicate IDs were preflighted");
    }

    let mut successors = BTreeMap::<ManufacturingStepIdV1, ManufacturingStepIdV1>::new();
    let mut roots = BTreeMap::<BodyId, Vec<ManufacturingStepIdV1>>::new();
    let mut steps_by_body = BTreeMap::<BodyId, usize>::new();
    for step in by_id.values() {
        *steps_by_body.entry(step.body.clone()).or_default() += 1;
        let Some(predecessor_id) = &step.predecessor else {
            roots
                .entry(step.body.clone())
                .or_default()
                .push(step.id.clone());
            continue;
        };
        let Some(predecessor) = by_id.get(predecessor_id) else {
            return Err(ManufacturingAdmissionErrorV1::MissingPredecessor {
                step: step.id.clone(),
                predecessor: predecessor_id.clone(),
            });
        };
        if predecessor.body != step.body {
            return Err(ManufacturingAdmissionErrorV1::CrossBodyPredecessor {
                step: step.id.clone(),
                predecessor: predecessor_id.clone(),
            });
        }
        if let Some(first) = successors.insert(predecessor_id.clone(), step.id.clone()) {
            let second = step.id.clone();
            let (first, second) = if first <= second {
                (first, second)
            } else {
                (second, first)
            };
            return Err(ManufacturingAdmissionErrorV1::ForkedProcessChain {
                predecessor: predecessor_id.clone(),
                first,
                second,
            });
        }
    }

    let mut ordered = Vec::with_capacity(by_id.len());
    for (body, count) in steps_by_body {
        let body_roots = roots.get(&body).map_or(&[][..], Vec::as_slice);
        if body_roots.len() != 1 {
            return Err(ManufacturingAdmissionErrorV1::RootCardinality {
                body,
                roots: body_roots.len(),
            });
        }
        let mut current = body_roots[0].clone();
        let mut visited = 0usize;
        loop {
            let step = by_id
                .get(&current)
                .expect("root and successors are drawn from admitted step IDs");
            ordered.push(step.clone());
            visited += 1;
            let Some(next) = successors.get(&current) else {
                break;
            };
            current = next.clone();
            if visited > count {
                return Err(ManufacturingAdmissionErrorV1::DisconnectedProcessChain { body });
            }
        }
        if visited != count {
            return Err(ManufacturingAdmissionErrorV1::DisconnectedProcessChain { body });
        }
    }

    let rows: Vec<Vec<u8>> = ordered
        .iter()
        .map(ManufacturingProcessStepV1::canonical_row)
        .collect();
    let correlation_row = draft.correlation_model.canonical_row();
    let graph_id = graph.identity();
    let receipt = CanonicalEncoder::<MachineManufacturingStateIdV1, _>::new(
        MANUFACTURING_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "manufacturing-schema-version"),
        u64::from(MACHINE_MANUFACTURING_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
    .bytes(Field::new(3, "correlation-model"), &correlation_row)?
    .ordered_bytes(
        Field::new(4, "process-steps"),
        rows.len() as u64,
        rows.iter().map(Vec::as_slice),
    )?
    .finish()?;

    Ok(AdmittedMachineManufacturingStateV1 {
        graph: graph_id,
        correlation_model: draft.correlation_model,
        process_steps: ordered,
        receipt,
    })
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
