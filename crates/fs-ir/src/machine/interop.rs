//! Checked engineering-workflow and quarantined interoperability receipts.
//!
//! This module binds the ordered E7 engineering workflow to one exact admitted
//! Machine graph. It also admits receipts for FMI 3.0.2 or SSP 2.0 artifact
//! coordinates supplied with a nominal adapter and isolation-receipt
//! coordinate. Foreign execution is never performed here, and its outputs can
//! enter only at [`ColorRank::Estimated`] with unbounded dispersion. Attempts
//! to claim `Validated` or `Verified` evidence refuse instead of being silently
//! demoted or laundered.
//!
//! Admission proves bounded structure, stable-ID closure, exact external
//! artifact coordinates, and deterministic receipt identity. It does not parse
//! FMI/SSP XML, authenticate an adapter or isolation receipt, execute an FMU,
//! establish unit compatibility, or validate a foreign model.

use core::fmt;
use core::num::NonZeroU64;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::{COLOR_ALGEBRA_VERSION, ColorRank};

use crate::IR_VERSION;

use super::{AdmittedMachineGraph, MachineGraphIdV1, MachineIdError, StateSlotId, TerminalId};

/// Identity-schema version for the admitted engineering workflow.
pub const MACHINE_WORKFLOW_SCHEMA_VERSION_V1: u32 = 1;
/// Identity-schema version for one quarantined foreign execution receipt.
pub const FOREIGN_EXECUTION_SCHEMA_VERSION_V1: u32 = 1;
/// Exact number of stages in the version-one engineering workflow.
pub const MACHINE_WORKFLOW_STAGE_COUNT_V1: usize = 10;
/// Maximum foreign output bindings retained by one execution receipt.
pub const MAX_FOREIGN_OUTPUT_BINDINGS_V1: usize = 4_096;
/// Mandatory marker bound into every foreign execution identity.
pub const FOREIGN_EXECUTION_NO_AUTHORITY_POLICY_V1: &str =
    "foreign-execution-estimated-only-no-native-authority";

// This is an interop-schema tag, not ColorRank's Rust discriminant.
const FOREIGN_ESTIMATED_RANK_CODE_V1: u64 = 1;

const WORKFLOW_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(256 * 1_024, 128 * 1_024, 4, 64, 64);
const FOREIGN_EXECUTION_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(8 * 1_024 * 1_024, 4 * 1_024 * 1_024, 4, 64, 8_192);

/// Refusal from constructing an exact external interoperability artifact ref.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InteropReferenceErrorV1 {
    /// The namespace violates the bounded Machine-IR key grammar.
    Namespace(MachineIdError),
    /// An all-zero digest cannot identify external artifact content.
    ZeroDigest,
}

impl InteropReferenceErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Namespace(_) => "InteropReferenceNamespace",
            Self::ZeroDigest => "InteropReferenceZeroDigest",
        }
    }
}

impl fmt::Display for InteropReferenceErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Namespace(error) => write!(formatter, "invalid interop namespace: {error}"),
            Self::ZeroDigest => formatter.write_str("interop artifact digest must not be all zero"),
        }
    }
}

impl std::error::Error for InteropReferenceErrorV1 {}

/// Bounded, versioned, content-addressed external artifact coordinate.
///
/// The coordinate is nominal. This module binds it exactly but does not
/// authenticate the owner, bytes, schema implementation, or process that
/// supplied it.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InteropArtifactRefV1 {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    content_hash: ContentHash,
}

impl InteropArtifactRefV1 {
    /// Construct one exact external artifact coordinate.
    ///
    /// # Errors
    /// Refuses a noncanonical namespace or all-zero digest.
    pub fn new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        content_hash: ContentHash,
    ) -> Result<Self, InteropReferenceErrorV1> {
        let namespace = namespace.into();
        super::validate_canonical_key("interop-artifact-ref", &namespace)
            .map_err(InteropReferenceErrorV1::Namespace)?;
        if content_hash.as_bytes() == &[0; 32] {
            return Err(InteropReferenceErrorV1::ZeroDigest);
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

    /// Exact caller-supplied artifact content hash.
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

/// Ordered stage in the version-one end-to-end engineering workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineWorkflowStageV1 {
    /// Import or create parameterized geometry.
    ImportOrCreateGeometry,
    /// Assemble bodies and declare joints/interfaces.
    AssembleAndJoint,
    /// Bind materials, interfaces, and as-built state.
    BindMachineState,
    /// Declare IC/BC, drives, loads, tolerances, faults, and hazards.
    DeclareScenario,
    /// Choose ContextOfUse, QoIs, acceptance criteria, and budgets.
    ChooseDecisionContext,
    /// Admit fidelity, solver, and escalation policy.
    AdmitExecutionPolicy,
    /// Calibrate only on declared training data.
    CalibrateTrainingData,
    /// Verify and validate against declared holdouts.
    ValidateHoldouts,
    /// Simulate, quantify uncertainty, or optimize under the admitted policy.
    SimulateUqOptimize,
    /// Export replay, evidence, and the decision package.
    ExportDecisionPackage,
}

impl MachineWorkflowStageV1 {
    /// The one canonical version-one stage order.
    pub const ORDERED: [Self; MACHINE_WORKFLOW_STAGE_COUNT_V1] = [
        Self::ImportOrCreateGeometry,
        Self::AssembleAndJoint,
        Self::BindMachineState,
        Self::DeclareScenario,
        Self::ChooseDecisionContext,
        Self::AdmitExecutionPolicy,
        Self::CalibrateTrainingData,
        Self::ValidateHoldouts,
        Self::SimulateUqOptimize,
        Self::ExportDecisionPackage,
    ];

    /// Stable machine-readable stage name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::ImportOrCreateGeometry => "import-or-create-geometry",
            Self::AssembleAndJoint => "assemble-and-joint",
            Self::BindMachineState => "bind-machine-state",
            Self::DeclareScenario => "declare-scenario",
            Self::ChooseDecisionContext => "choose-decision-context",
            Self::AdmitExecutionPolicy => "admit-execution-policy",
            Self::CalibrateTrainingData => "calibrate-training-data",
            Self::ValidateHoldouts => "validate-holdouts",
            Self::SimulateUqOptimize => "simulate-uq-optimize",
            Self::ExportDecisionPackage => "export-decision-package",
        }
    }

    const fn code(self) -> u64 {
        match self {
            Self::ImportOrCreateGeometry => 1,
            Self::AssembleAndJoint => 2,
            Self::BindMachineState => 3,
            Self::DeclareScenario => 4,
            Self::ChooseDecisionContext => 5,
            Self::AdmitExecutionPolicy => 6,
            Self::CalibrateTrainingData => 7,
            Self::ValidateHoldouts => 8,
            Self::SimulateUqOptimize => 9,
            Self::ExportDecisionPackage => 10,
        }
    }
}

/// One asserted stage artifact in the checked workflow spine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineWorkflowStepV1 {
    /// Exact stage represented by this artifact.
    pub stage: MachineWorkflowStageV1,
    /// Caller-asserted content-addressed artifact coordinate for this stage.
    pub artifact: InteropArtifactRefV1,
}

impl MachineWorkflowStepV1 {
    fn canonical_row(&self) -> Vec<u8> {
        let artifact = self.artifact.canonical_row();
        let mut row = Vec::with_capacity(16 + artifact.len());
        row.extend_from_slice(&self.stage.code().to_le_bytes());
        append_bytes(&mut row, &artifact);
        row
    }
}

/// Draft of the complete checked engineering workflow spine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineWorkflowDraftV1 {
    /// Exact Machine graph the workflow claims to extend.
    pub machine_graph_id: MachineGraphIdV1,
    /// Caller order must equal [`MachineWorkflowStageV1::ORDERED`].
    pub steps: Vec<MachineWorkflowStepV1>,
}

/// Canonical identity schema for one admitted engineering workflow.
pub enum MachineWorkflowIdentitySchemaV1 {}

impl CanonicalSchema for MachineWorkflowIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.workflow.v1";
    const NAME: &'static str = "admitted-machine-workflow";
    const VERSION: u32 = MACHINE_WORKFLOW_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "exact Machine graph and complete ordered content-addressed engineering workflow";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("workflow-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("ordered-stage-artifacts", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted engineering workflow.
pub type MachineWorkflowIdV1 = ProblemSemanticId<MachineWorkflowIdentitySchemaV1>;

/// Refusal from admitting a complete engineering workflow spine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineWorkflowRefusalV1 {
    /// The draft names a different Machine graph than the supplied authority.
    MachineGraphMismatch {
        /// Graph named by the draft.
        declared: MachineGraphIdV1,
        /// Exact supplied admitted graph.
        supplied: MachineGraphIdV1,
    },
    /// The workflow omitted or added a stage.
    StageCount {
        /// Submitted stage count.
        actual: usize,
        /// Exact version-one stage count.
        expected: usize,
    },
    /// One position does not contain the required stage.
    StageOrder {
        /// Zero-based position.
        index: usize,
        /// Required stage at this position.
        expected: MachineWorkflowStageV1,
        /// Submitted stage at this position.
        actual: MachineWorkflowStageV1,
    },
    /// Canonical identity publication was refused.
    Identity(CanonicalError),
}

impl MachineWorkflowRefusalV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::MachineGraphMismatch { .. } => "MachineWorkflowGraphMismatch",
            Self::StageCount { .. } => "MachineWorkflowStageCount",
            Self::StageOrder { .. } => "MachineWorkflowStageOrder",
            Self::Identity(_) => "MachineWorkflowIdentity",
        }
    }

    /// Actionable deterministic repair hint.
    #[must_use]
    pub const fn fix(&self) -> &'static str {
        match self {
            Self::MachineGraphMismatch { .. } => {
                "bind the workflow draft to the exact supplied admitted Machine graph"
            }
            Self::StageCount { .. } | Self::StageOrder { .. } => {
                "supply every version-one workflow stage exactly once in canonical order"
            }
            Self::Identity(_) => "reduce the bounded workflow payload and retry admission",
        }
    }
}

impl fmt::Display for MachineWorkflowRefusalV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MachineGraphMismatch { .. } => {
                formatter.write_str("workflow draft is bound to another Machine graph")
            }
            Self::StageCount { actual, expected } => write!(
                formatter,
                "workflow has {actual} stages; version one requires {expected}"
            ),
            Self::StageOrder {
                index,
                expected,
                actual,
            } => write!(
                formatter,
                "workflow stage {index} is {}; expected {}",
                actual.name(),
                expected.name()
            ),
            Self::Identity(error) => write!(formatter, "workflow identity refused: {error}"),
        }
    }
}

impl std::error::Error for MachineWorkflowRefusalV1 {}

/// Admitted exact-graph, exact-order engineering workflow spine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineWorkflowV1 {
    machine_graph_id: MachineGraphIdV1,
    steps: Box<[MachineWorkflowStepV1]>,
    receipt: IdentityReceipt<MachineWorkflowIdV1>,
}

impl MachineWorkflowDraftV1 {
    /// Admit the complete workflow against one exact Machine graph.
    ///
    /// # Errors
    /// Refuses a graph rebind, missing/extra/reordered stage, or bounded
    /// identity-publication failure.
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineWorkflowV1, MachineWorkflowRefusalV1> {
        let supplied = graph.identity();
        if self.machine_graph_id != supplied {
            return Err(MachineWorkflowRefusalV1::MachineGraphMismatch {
                declared: self.machine_graph_id,
                supplied,
            });
        }
        if self.steps.len() != MACHINE_WORKFLOW_STAGE_COUNT_V1 {
            return Err(MachineWorkflowRefusalV1::StageCount {
                actual: self.steps.len(),
                expected: MACHINE_WORKFLOW_STAGE_COUNT_V1,
            });
        }
        for (index, (step, expected)) in self
            .steps
            .iter()
            .zip(MachineWorkflowStageV1::ORDERED)
            .enumerate()
        {
            if step.stage != expected {
                return Err(MachineWorkflowRefusalV1::StageOrder {
                    index,
                    expected,
                    actual: step.stage,
                });
            }
        }

        let rows = self
            .steps
            .iter()
            .map(MachineWorkflowStepV1::canonical_row)
            .collect::<Vec<_>>();
        let receipt = encode_workflow_identity(supplied, &rows)
            .map_err(MachineWorkflowRefusalV1::Identity)?;
        Ok(AdmittedMachineWorkflowV1 {
            machine_graph_id: supplied,
            steps: self.steps.into_boxed_slice(),
            receipt,
        })
    }
}

impl AdmittedMachineWorkflowV1 {
    /// Exact admitted Machine graph identity.
    #[must_use]
    pub const fn machine_graph_id(&self) -> MachineGraphIdV1 {
        self.machine_graph_id
    }

    /// Complete canonical workflow stages and their artifact coordinates.
    #[must_use]
    pub fn steps(&self) -> &[MachineWorkflowStepV1] {
        &self.steps
    }

    /// Strong semantic identity of this workflow spine.
    #[must_use]
    pub const fn identity(&self) -> MachineWorkflowIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineWorkflowIdV1> {
        self.receipt
    }
}

/// Closed caller-declared interchange standard/version coordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InterchangeStandardV1 {
    /// Functional Mock-up Interface 3.0.2.
    Fmi302,
    /// System Structure and Parameterization 2.0.
    Ssp20,
}

impl InterchangeStandardV1 {
    /// Stable standard family name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Fmi302 => "fmi",
            Self::Ssp20 => "ssp",
        }
    }

    /// Exact admitted external standard version.
    #[must_use]
    pub const fn version(self) -> &'static str {
        match self {
            Self::Fmi302 => "3.0.2",
            Self::Ssp20 => "2.0",
        }
    }

    const fn code(self) -> u64 {
        match self {
            Self::Fmi302 => 1,
            Self::Ssp20 => 2,
        }
    }
}

/// Stable Machine-IR target associated with one quarantined foreign output.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ForeignOutputTargetV1 {
    /// A declared Machine graph terminal.
    Terminal(TerminalId),
    /// A state slot directly owned by one declared subsystem.
    StateSlot(StateSlotId),
}

impl ForeignOutputTargetV1 {
    /// Stable target role.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        match self {
            Self::Terminal(_) => "terminal",
            Self::StateSlot(_) => "state-slot",
        }
    }

    /// Human-auditable canonical Machine key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        match self {
            Self::Terminal(id) => id.canonical_key(),
            Self::StateSlot(id) => id.canonical_key(),
        }
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(self.canonical_key().len() + 48);
        match self {
            Self::Terminal(id) => {
                row.extend_from_slice(&1_u64.to_le_bytes());
                append_bytes(&mut row, id.canonical_key().as_bytes());
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::StateSlot(id) => {
                row.extend_from_slice(&2_u64.to_le_bytes());
                append_bytes(&mut row, id.canonical_key().as_bytes());
                row.extend_from_slice(id.identity().as_bytes());
            }
        }
        row
    }

    fn exists_in(&self, graph: &AdmittedMachineGraph) -> bool {
        match self {
            Self::Terminal(id) => graph.terminals().iter().any(|terminal| &terminal.id == id),
            Self::StateSlot(id) => graph.subsystems().iter().any(|subsystem| {
                subsystem
                    .state_slots
                    .iter()
                    .any(|candidate| candidate == id)
            }),
        }
    }
}

/// One opaque foreign output artifact bound to a stable Machine target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ForeignOutputBindingV1 {
    /// Bounded canonical output name, unique within one execution.
    pub name: String,
    /// Exact Machine terminal or state slot associated with the output.
    pub target: ForeignOutputTargetV1,
    /// Opaque content-addressed output artifact.
    pub artifact: InteropArtifactRefV1,
}

impl ForeignOutputBindingV1 {
    fn canonical_row(&self) -> Vec<u8> {
        let target = self.target.canonical_row();
        let artifact = self.artifact.canonical_row();
        let mut row = Vec::with_capacity(self.name.len() + target.len() + artifact.len() + 24);
        append_bytes(&mut row, self.name.as_bytes());
        append_bytes(&mut row, &target);
        append_bytes(&mut row, &artifact);
        row
    }
}

/// Draft receipt for output coordinates claimed from outside the trusted graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ForeignExecutionDraftV1 {
    /// Exact Machine graph the outputs claim to target.
    pub machine_graph_id: MachineGraphIdV1,
    /// Exact admitted workflow associated with the boundary crossing.
    pub workflow_id: MachineWorkflowIdV1,
    /// Caller-declared FMI 3.0.2 or SSP 2.0 coordinate.
    pub standard: InterchangeStandardV1,
    /// Caller-supplied coordinate for imported model-description/system-structure bytes.
    pub model_description: InteropArtifactRefV1,
    /// Caller-supplied adapter implementation coordinate.
    pub adapter: InteropArtifactRefV1,
    /// Caller-supplied nominal isolation/execution receipt coordinate.
    pub isolation_receipt: InteropArtifactRefV1,
    /// Evidence strength requested by the caller. Only `Estimated` is legal.
    pub claimed_rank: ColorRank,
    /// Opaque output artifacts and their stable Machine targets.
    pub outputs: Vec<ForeignOutputBindingV1>,
}

/// Canonical identity schema for one quarantined foreign execution receipt.
pub enum ForeignExecutionIdentitySchemaV1 {}

impl CanonicalSchema for ForeignExecutionIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.foreign-execution.v1";
    const NAME: &'static str = "admitted-foreign-execution";
    const VERSION: u32 = FOREIGN_EXECUTION_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "exact Machine workflow and bounded FMI or SSP output quarantine with Estimated-only evidence";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("foreign-execution-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("machine-workflow", WireType::Bytes),
        FieldSpec::required("interchange-standard", WireType::U64),
        FieldSpec::required("model-description", WireType::Bytes),
        FieldSpec::required("adapter", WireType::Bytes),
        FieldSpec::required("isolation-receipt", WireType::Bytes),
        FieldSpec::required("color-algebra-version", WireType::U64),
        FieldSpec::required("evidence-rank", WireType::U64),
        FieldSpec::required("estimator-identity", WireType::Utf8),
        FieldSpec::required("estimated-dispersion-bits", WireType::U64),
        FieldSpec::required("no-authority-policy", WireType::Utf8),
        FieldSpec::required("canonical-output-bindings", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one quarantined foreign execution.
pub type ForeignExecutionIdV1 = ProblemSemanticId<ForeignExecutionIdentitySchemaV1>;

/// Fail-closed refusal from admitting foreign outputs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ForeignExecutionRefusalV1 {
    /// The supplied workflow was admitted against another graph.
    WorkflowGraphMismatch,
    /// The draft names another graph.
    MachineGraphMismatch,
    /// The draft names another workflow.
    WorkflowIdentityMismatch,
    /// An execution with no outputs has no consumable boundary artifact.
    NoOutputs,
    /// The public output-binding cap was exceeded.
    TooManyOutputs {
        /// Submitted bindings.
        actual: usize,
        /// Versioned maximum.
        max: usize,
    },
    /// Foreign output tried to claim native certificate authority.
    EvidenceLaundering {
        /// Caller-requested rank.
        claimed: ColorRank,
        /// Strongest rank this boundary can admit.
        maximum: ColorRank,
    },
    /// One output name violates the bounded canonical key grammar.
    InvalidOutputName {
        /// Submitted position before canonicalization.
        index: usize,
        /// Exact naming refusal.
        error: MachineIdError,
    },
    /// Two outputs share a name.
    DuplicateOutputName {
        /// First submitted position.
        first_index: usize,
        /// Duplicate submitted position.
        duplicate_index: usize,
        /// Bounded canonical name.
        name: String,
    },
    /// Two outputs compete for the same stable target.
    DuplicateOutputTarget {
        /// First submitted position.
        first_index: usize,
        /// Duplicate submitted position.
        duplicate_index: usize,
        /// Repeated stable target.
        target: ForeignOutputTargetV1,
    },
    /// An output target is absent from the supplied admitted graph.
    UnknownOutputTarget {
        /// Submitted position.
        index: usize,
        /// Foreign target.
        target: ForeignOutputTargetV1,
    },
    /// Canonical identity publication was refused.
    Identity(CanonicalError),
}

impl ForeignExecutionRefusalV1 {
    /// Stable diagnostic rule code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::WorkflowGraphMismatch => "ForeignExecutionWorkflowGraphMismatch",
            Self::MachineGraphMismatch => "ForeignExecutionGraphMismatch",
            Self::WorkflowIdentityMismatch => "ForeignExecutionWorkflowMismatch",
            Self::NoOutputs => "ForeignExecutionNoOutputs",
            Self::TooManyOutputs { .. } => "ForeignExecutionOutputLimit",
            Self::EvidenceLaundering { .. } => "ForeignExecutionEvidenceLaundering",
            Self::InvalidOutputName { .. } => "ForeignExecutionOutputName",
            Self::DuplicateOutputName { .. } => "ForeignExecutionDuplicateOutputName",
            Self::DuplicateOutputTarget { .. } => "ForeignExecutionDuplicateOutputTarget",
            Self::UnknownOutputTarget { .. } => "ForeignExecutionUnknownOutputTarget",
            Self::Identity(_) => "ForeignExecutionIdentity",
        }
    }

    /// Actionable deterministic repair hint.
    #[must_use]
    pub const fn fix(&self) -> &'static str {
        match self {
            Self::WorkflowGraphMismatch
            | Self::MachineGraphMismatch
            | Self::WorkflowIdentityMismatch => {
                "bind the draft, admitted workflow, and supplied Machine graph to the same exact identities"
            }
            Self::NoOutputs => "supply at least one content-addressed foreign output binding",
            Self::TooManyOutputs { .. } => {
                "split the external execution into bounded independently receipted batches"
            }
            Self::EvidenceLaundering { .. } => {
                "declare foreign outputs Estimated and retain native validation or verification as a separate downstream operation"
            }
            Self::InvalidOutputName { .. } => {
                "use a bounded lowercase slash-separated canonical output name"
            }
            Self::DuplicateOutputName { .. } => "give every foreign output a unique stable name",
            Self::DuplicateOutputTarget { .. } => {
                "merge competing foreign values before binding one output to the Machine target"
            }
            Self::UnknownOutputTarget { .. } => {
                "bind the output to a terminal or state slot in the supplied admitted Machine graph"
            }
            Self::Identity(_) => "reduce the bounded execution payload and retry admission",
        }
    }
}

impl fmt::Display for ForeignExecutionRefusalV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WorkflowGraphMismatch => {
                formatter.write_str("workflow belongs to another Machine graph")
            }
            Self::MachineGraphMismatch => {
                formatter.write_str("foreign execution draft names another Machine graph")
            }
            Self::WorkflowIdentityMismatch => {
                formatter.write_str("foreign execution draft names another workflow")
            }
            Self::NoOutputs => formatter.write_str("foreign execution contains no outputs"),
            Self::TooManyOutputs { actual, max } => {
                write!(
                    formatter,
                    "foreign execution has {actual} outputs; maximum is {max}"
                )
            }
            Self::EvidenceLaundering { claimed, maximum } => write!(
                formatter,
                "foreign execution claimed {claimed:?}; maximum admissible rank is {maximum:?}"
            ),
            Self::InvalidOutputName { index, error } => {
                write!(
                    formatter,
                    "foreign output {index} has an invalid name: {error}"
                )
            }
            Self::DuplicateOutputName {
                first_index,
                duplicate_index,
                name,
            } => write!(
                formatter,
                "foreign outputs {first_index} and {duplicate_index} share name {name:?}"
            ),
            Self::DuplicateOutputTarget {
                first_index,
                duplicate_index,
                target,
            } => write!(
                formatter,
                "foreign outputs {first_index} and {duplicate_index} share target {}:{}",
                target.kind(),
                target.canonical_key()
            ),
            Self::UnknownOutputTarget { index, target } => write!(
                formatter,
                "foreign output {index} targets unknown {}:{}",
                target.kind(),
                target.canonical_key()
            ),
            Self::Identity(error) => {
                write!(formatter, "foreign execution identity refused: {error}")
            }
        }
    }
}

impl std::error::Error for ForeignExecutionRefusalV1 {}

/// Admitted foreign execution receipt whose outputs remain Estimated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedForeignExecutionV1 {
    machine_graph_id: MachineGraphIdV1,
    workflow_id: MachineWorkflowIdV1,
    standard: InterchangeStandardV1,
    model_description: InteropArtifactRefV1,
    adapter: InteropArtifactRefV1,
    isolation_receipt: InteropArtifactRefV1,
    outputs: Box<[ForeignOutputBindingV1]>,
    estimator_identity: Box<str>,
    receipt: IdentityReceipt<ForeignExecutionIdV1>,
}

impl ForeignExecutionDraftV1 {
    /// Admit foreign outputs against one exact workflow and Machine graph.
    ///
    /// # Errors
    /// Refuses identity rebinding, empty/oversized/ambiguous output sets,
    /// targets outside the graph, evidence stronger than `Estimated`, or a
    /// bounded canonical-identity failure.
    pub fn admit_against(
        self,
        workflow: &AdmittedMachineWorkflowV1,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedForeignExecutionV1, ForeignExecutionRefusalV1> {
        if workflow.machine_graph_id() != graph.identity() {
            return Err(ForeignExecutionRefusalV1::WorkflowGraphMismatch);
        }
        if self.machine_graph_id != graph.identity() {
            return Err(ForeignExecutionRefusalV1::MachineGraphMismatch);
        }
        if self.workflow_id != workflow.identity() {
            return Err(ForeignExecutionRefusalV1::WorkflowIdentityMismatch);
        }
        if self.outputs.is_empty() {
            return Err(ForeignExecutionRefusalV1::NoOutputs);
        }
        if self.outputs.len() > MAX_FOREIGN_OUTPUT_BINDINGS_V1 {
            return Err(ForeignExecutionRefusalV1::TooManyOutputs {
                actual: self.outputs.len(),
                max: MAX_FOREIGN_OUTPUT_BINDINGS_V1,
            });
        }
        if self.claimed_rank != ColorRank::Estimated {
            return Err(ForeignExecutionRefusalV1::EvidenceLaundering {
                claimed: self.claimed_rank,
                maximum: ColorRank::Estimated,
            });
        }

        let mut names = BTreeMap::<String, usize>::new();
        let mut targets = BTreeMap::<ForeignOutputTargetV1, usize>::new();
        for (index, output) in self.outputs.iter().enumerate() {
            super::validate_canonical_key("foreign-output", &output.name)
                .map_err(|error| ForeignExecutionRefusalV1::InvalidOutputName { index, error })?;
            if let Some(&first_index) = names.get(&output.name) {
                return Err(ForeignExecutionRefusalV1::DuplicateOutputName {
                    first_index,
                    duplicate_index: index,
                    name: output.name.clone(),
                });
            }
            if !output.target.exists_in(graph) {
                return Err(ForeignExecutionRefusalV1::UnknownOutputTarget {
                    index,
                    target: output.target.clone(),
                });
            }
            if let Some(&first_index) = targets.get(&output.target) {
                return Err(ForeignExecutionRefusalV1::DuplicateOutputTarget {
                    first_index,
                    duplicate_index: index,
                    target: output.target.clone(),
                });
            }
            names.insert(output.name.clone(), index);
            targets.insert(output.target.clone(), index);
        }

        let mut outputs = self.outputs;
        outputs.sort_by(|left, right| {
            left.name
                .cmp(&right.name)
                .then_with(|| left.target.cmp(&right.target))
        });
        let rows = outputs
            .iter()
            .map(ForeignOutputBindingV1::canonical_row)
            .collect::<Vec<_>>();
        let model_row = self.model_description.canonical_row();
        let adapter_row = self.adapter.canonical_row();
        let isolation_row = self.isolation_receipt.canonical_row();
        let estimator_hash = hash_domain(
            "org.frankensim.fs-ir.machine.external-adapter-estimator.v1",
            &adapter_row,
        );
        let estimator_identity = format!("external-adapter:{}", estimator_hash.to_hex());
        let receipt = encode_foreign_execution_identity(
            graph.identity(),
            workflow.identity(),
            self.standard,
            &model_row,
            &adapter_row,
            &isolation_row,
            &estimator_identity,
            &rows,
        )
        .map_err(ForeignExecutionRefusalV1::Identity)?;
        Ok(AdmittedForeignExecutionV1 {
            machine_graph_id: graph.identity(),
            workflow_id: workflow.identity(),
            standard: self.standard,
            model_description: self.model_description,
            adapter: self.adapter,
            isolation_receipt: self.isolation_receipt,
            outputs: outputs.into_boxed_slice(),
            estimator_identity: estimator_identity.into_boxed_str(),
            receipt,
        })
    }
}

impl AdmittedForeignExecutionV1 {
    /// Exact Machine graph whose stable targets the outputs reference.
    #[must_use]
    pub const fn machine_graph_id(&self) -> MachineGraphIdV1 {
        self.machine_graph_id
    }

    /// Exact checked workflow authorizing the boundary crossing.
    #[must_use]
    pub const fn workflow_id(&self) -> MachineWorkflowIdV1 {
        self.workflow_id
    }

    /// Exact caller-declared external standard/version coordinate.
    #[must_use]
    pub const fn standard(&self) -> InterchangeStandardV1 {
        self.standard
    }

    /// Caller-supplied imported model-description coordinate.
    #[must_use]
    pub const fn model_description(&self) -> &InteropArtifactRefV1 {
        &self.model_description
    }

    /// Caller-supplied adapter coordinate.
    #[must_use]
    pub const fn adapter(&self) -> &InteropArtifactRefV1 {
        &self.adapter
    }

    /// Caller-supplied nominal isolation receipt coordinate.
    #[must_use]
    pub const fn isolation_receipt(&self) -> &InteropArtifactRefV1 {
        &self.isolation_receipt
    }

    /// Canonically ordered opaque foreign outputs and their stable targets.
    #[must_use]
    pub fn outputs(&self) -> &[ForeignOutputBindingV1] {
        &self.outputs
    }

    /// Forced evidence rank. No stronger variant can be stored in this type.
    #[must_use]
    pub const fn evidence_rank(&self) -> ColorRank {
        ColorRank::Estimated
    }

    /// Deterministic external-adapter estimator identity.
    #[must_use]
    pub fn estimator_identity(&self) -> &str {
        &self.estimator_identity
    }

    /// Foreign execution carries no finite dispersion claim.
    #[must_use]
    pub const fn estimated_dispersion(&self) -> f64 {
        f64::INFINITY
    }

    /// Mandatory policy marker denying native authority inheritance.
    #[must_use]
    pub const fn no_authority_policy(&self) -> &'static str {
        FOREIGN_EXECUTION_NO_AUTHORITY_POLICY_V1
    }

    /// Strong semantic identity of this complete quarantine receipt.
    #[must_use]
    pub const fn identity(&self) -> ForeignExecutionIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<ForeignExecutionIdV1> {
        self.receipt
    }
}

fn encode_workflow_identity(
    graph: MachineGraphIdV1,
    rows: &[Vec<u8>],
) -> Result<IdentityReceipt<MachineWorkflowIdV1>, CanonicalError> {
    CanonicalEncoder::<MachineWorkflowIdV1, _>::new(WORKFLOW_IDENTITY_LIMITS, NeverCancel)?
        .u64(
            Field::new(0, "workflow-schema-version"),
            u64::from(MACHINE_WORKFLOW_SCHEMA_VERSION_V1),
        )?
        .u64(
            Field::new(1, "frankenscript-ir-version"),
            u64::from(IR_VERSION),
        )?
        .bytes(Field::new(2, "machine-graph"), graph.as_bytes())?
        .ordered_bytes(
            Field::new(3, "ordered-stage-artifacts"),
            rows.len() as u64,
            rows.iter().map(Vec::as_slice),
        )?
        .finish()
}

#[allow(clippy::too_many_arguments)]
fn encode_foreign_execution_identity(
    graph: MachineGraphIdV1,
    workflow: MachineWorkflowIdV1,
    standard: InterchangeStandardV1,
    model_description: &[u8],
    adapter: &[u8],
    isolation_receipt: &[u8],
    estimator_identity: &str,
    outputs: &[Vec<u8>],
) -> Result<IdentityReceipt<ForeignExecutionIdV1>, CanonicalError> {
    CanonicalEncoder::<ForeignExecutionIdV1, _>::new(
        FOREIGN_EXECUTION_IDENTITY_LIMITS,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "foreign-execution-schema-version"),
        u64::from(FOREIGN_EXECUTION_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph.as_bytes())?
    .bytes(Field::new(3, "machine-workflow"), workflow.as_bytes())?
    .u64(Field::new(4, "interchange-standard"), standard.code())?
    .bytes(Field::new(5, "model-description"), model_description)?
    .bytes(Field::new(6, "adapter"), adapter)?
    .bytes(Field::new(7, "isolation-receipt"), isolation_receipt)?
    .u64(
        Field::new(8, "color-algebra-version"),
        u64::from(COLOR_ALGEBRA_VERSION),
    )?
    .u64(
        Field::new(9, "evidence-rank"),
        FOREIGN_ESTIMATED_RANK_CODE_V1,
    )?
    .utf8(Field::new(10, "estimator-identity"), estimator_identity)?
    .u64(
        Field::new(11, "estimated-dispersion-bits"),
        f64::INFINITY.to_bits(),
    )?
    .utf8(
        Field::new(12, "no-authority-policy"),
        FOREIGN_EXECUTION_NO_AUTHORITY_POLICY_V1,
    )?
    .ordered_bytes(
        Field::new(13, "canonical-output-bindings"),
        outputs.len() as u64,
        outputs.iter().map(Vec::as_slice),
    )?
    .finish()
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
