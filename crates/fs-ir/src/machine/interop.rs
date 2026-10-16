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
//! artifact coordinates, and deterministic receipt identity. The additive
//! model-description preflight checks bounded UTF-8 bytes, the required root
//! coordinate, and exact raw-byte identity before a typed checked draft can be
//! admitted. It is not full XML/XSD or archive validation and does not
//! authenticate an adapter or isolation receipt, execute an FMU, establish unit
//! compatibility, or validate a foreign model.

use core::fmt;
use core::num::NonZeroU64;

use std::collections::BTreeMap;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_blake3::{ContentHash, hash_bytes, hash_domain};
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
/// Maximum raw FMI model-description or SSP system-structure bytes preflighted.
pub const MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1: usize = 4 * 1_024 * 1_024;
/// Maximum bytes inspected before and through the root start tag.
pub const MAX_INTEROP_ROOT_PREFIX_BYTES_V1: usize = 64 * 1_024;
/// Maximum attributes accepted on an interop document root.
pub const MAX_INTEROP_ROOT_ATTRIBUTES_V1: usize = 64;
/// Mandatory marker bound into every foreign execution identity.
pub const FOREIGN_EXECUTION_NO_AUTHORITY_POLICY_V1: &str =
    "foreign-execution-estimated-only-no-native-authority";

const SSP20_SYSTEM_STRUCTURE_NAMESPACE_V1: &str =
    "http://ssp-standard.org/SSP1/SystemStructureDescription";

// This is an interop-schema tag, not ColorRank's Rust discriminant.
const FOREIGN_ESTIMATED_RANK_CODE_V1: u64 = 1;

const WORKFLOW_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(256 * 1_024, 128 * 1_024, 4, 64, 64);
// max_fields must admit every ForeignExecutionIdentitySchemaV1 field (14).
const FOREIGN_EXECUTION_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(8 * 1_024 * 1_024, 4 * 1_024 * 1_024, 14, 64, 8_192);

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

    /// Root element required by this coordinate's model-description document.
    #[must_use]
    pub const fn model_description_root(self) -> &'static str {
        match self {
            Self::Fmi302 => "fmiModelDescription",
            Self::Ssp20 => "SystemStructureDescription",
        }
    }

    /// Version written in this coordinate's model-description root attribute.
    ///
    /// FMI patch releases share the `fmiVersion="3.0"` XML schema coordinate,
    /// so this deliberately differs from [`Self::version`] for FMI 3.0.2.
    #[must_use]
    pub const fn model_description_version(self) -> &'static str {
        match self {
            Self::Fmi302 => "3.0",
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

/// Fail-closed refusal from bounded FMI/SSP root-coordinate preflight.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelDescriptionPreflightRefusalV1 {
    /// Raw model-description bytes exceeded the public input cap.
    TooLarge {
        /// Submitted byte length.
        actual: usize,
        /// Versioned maximum byte length.
        max: usize,
    },
    /// Input was not valid UTF-8.
    InvalidUtf8 {
        /// First byte not known to be valid UTF-8.
        valid_up_to: usize,
    },
    /// XML-forbidden NUL appeared in the input.
    NulByte {
        /// Exact byte offset.
        offset: usize,
    },
    /// A DTD or entity declaration was refused by the structural preflight.
    ForbiddenDeclaration {
        /// Stable declaration kind.
        kind: &'static str,
        /// Exact byte offset.
        offset: usize,
    },
    /// FMI modelDescription.xml omitted its required leading XML declaration.
    MissingXmlDeclaration,
    /// The root start tag was not reached within the public prefix bound.
    RootPrefixTooLarge {
        /// Versioned maximum prefix length.
        max: usize,
    },
    /// The XML declaration or root start tag violated the bounded grammar.
    Malformed {
        /// Exact byte offset where the refusal was detected.
        offset: usize,
        /// Stable bounded grammar rule.
        rule: &'static str,
    },
    /// The root start tag exceeded the public attribute-count cap.
    TooManyAttributes {
        /// Submitted attribute count at refusal.
        actual: usize,
        /// Versioned maximum attribute count.
        max: usize,
    },
    /// Two root attributes used the same exact qualified name.
    DuplicateAttribute {
        /// Repeated bounded attribute name.
        name: String,
    },
    /// The first element was not the standard's required root.
    WrongRoot {
        /// Required unqualified or local root name.
        expected: &'static str,
        /// Submitted bounded qualified root name.
        actual: String,
    },
    /// A required root attribute was absent.
    MissingAttribute {
        /// Required attribute name.
        name: &'static str,
    },
    /// A required identifying attribute was empty.
    EmptyAttribute {
        /// Required attribute name.
        name: &'static str,
    },
    /// A required attribute did not have the exact standard value.
    WrongAttributeValue {
        /// Checked attribute name.
        name: String,
        /// Required exact value.
        expected: &'static str,
        /// Submitted bounded value.
        actual: String,
    },
    /// The caller-declared artifact digest did not identify the inspected bytes.
    ContentHashMismatch {
        /// Digest bound by the draft's artifact coordinate.
        declared: ContentHash,
        /// Digest of the exact inspected bytes.
        actual: ContentHash,
    },
}

impl ModelDescriptionPreflightRefusalV1 {
    /// Stable diagnostic rule code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::TooLarge { .. } => "InteropModelDescriptionSize",
            Self::InvalidUtf8 { .. } => "InteropModelDescriptionUtf8",
            Self::NulByte { .. } => "InteropModelDescriptionNul",
            Self::ForbiddenDeclaration { .. } => "InteropModelDescriptionDeclaration",
            Self::MissingXmlDeclaration => "InteropModelDescriptionXmlDeclaration",
            Self::RootPrefixTooLarge { .. } => "InteropModelDescriptionRootPrefix",
            Self::Malformed { .. } => "InteropModelDescriptionMalformed",
            Self::TooManyAttributes { .. } => "InteropModelDescriptionAttributeLimit",
            Self::DuplicateAttribute { .. } => "InteropModelDescriptionDuplicateAttribute",
            Self::WrongRoot { .. } => "InteropModelDescriptionRoot",
            Self::MissingAttribute { .. } => "InteropModelDescriptionMissingAttribute",
            Self::EmptyAttribute { .. } => "InteropModelDescriptionEmptyAttribute",
            Self::WrongAttributeValue { .. } => "InteropModelDescriptionAttributeValue",
            Self::ContentHashMismatch { .. } => "InteropModelDescriptionContentHash",
        }
    }

    /// Actionable deterministic repair hint.
    #[must_use]
    pub const fn fix(&self) -> &'static str {
        match self {
            Self::TooLarge { .. } => {
                "supply the bounded model-description XML separately from archive payloads"
            }
            Self::InvalidUtf8 { .. } | Self::NulByte { .. } => {
                "encode the model-description document as valid NUL-free UTF-8"
            }
            Self::ForbiddenDeclaration { .. } => {
                "remove DTD and entity declarations before structural preflight"
            }
            Self::MissingXmlDeclaration => {
                "start FMI modelDescription.xml with an XML 1.0 UTF-8 declaration"
            }
            Self::RootPrefixTooLarge { .. } | Self::TooManyAttributes { .. } => {
                "reduce prolog or root-tag complexity and retry"
            }
            Self::Malformed { .. } | Self::DuplicateAttribute { .. } => {
                "supply a single bounded root start tag with unique quoted attributes"
            }
            Self::WrongRoot { .. }
            | Self::MissingAttribute { .. }
            | Self::EmptyAttribute { .. }
            | Self::WrongAttributeValue { .. } => {
                "supply the exact root, namespace, version, and identifying attributes required by the declared standard"
            }
            Self::ContentHashMismatch { .. } => {
                "bind the model-description artifact coordinate to the exact inspected raw bytes"
            }
        }
    }
}

impl fmt::Display for ModelDescriptionPreflightRefusalV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TooLarge { actual, max } => {
                write!(
                    formatter,
                    "model-description has {actual} bytes; maximum is {max}"
                )
            }
            Self::InvalidUtf8 { valid_up_to } => {
                write!(
                    formatter,
                    "model-description is not UTF-8 at byte {valid_up_to}"
                )
            }
            Self::NulByte { offset } => {
                write!(formatter, "model-description contains NUL at byte {offset}")
            }
            Self::ForbiddenDeclaration { kind, offset } => write!(
                formatter,
                "model-description contains forbidden {kind} declaration at byte {offset}"
            ),
            Self::MissingXmlDeclaration => formatter.write_str(
                "FMI modelDescription.xml is missing its leading XML 1.0 UTF-8 declaration",
            ),
            Self::RootPrefixTooLarge { max } => write!(
                formatter,
                "model-description root was not reached within {max} bytes"
            ),
            Self::Malformed { offset, rule } => {
                write!(
                    formatter,
                    "malformed model-description at byte {offset}: {rule}"
                )
            }
            Self::TooManyAttributes { actual, max } => write!(
                formatter,
                "model-description root has at least {actual} attributes; maximum is {max}"
            ),
            Self::DuplicateAttribute { name } => {
                write!(
                    formatter,
                    "model-description root repeats attribute {name:?}"
                )
            }
            Self::WrongRoot { expected, actual } => write!(
                formatter,
                "model-description root is {actual:?}; expected {expected:?}"
            ),
            Self::MissingAttribute { name } => {
                write!(
                    formatter,
                    "model-description root is missing attribute {name:?}"
                )
            }
            Self::EmptyAttribute { name } => {
                write!(
                    formatter,
                    "model-description root attribute {name:?} is empty"
                )
            }
            Self::WrongAttributeValue {
                name,
                expected,
                actual,
            } => write!(
                formatter,
                "model-description root attribute {name:?} is {actual:?}; expected {expected:?}"
            ),
            Self::ContentHashMismatch { declared, actual } => write!(
                formatter,
                "model-description digest {} does not match inspected bytes {}",
                declared.to_hex(),
                actual.to_hex()
            ),
        }
    }
}

impl std::error::Error for ModelDescriptionPreflightRefusalV1 {}

/// Opaque evidence that exact bounded bytes passed structural root preflight.
///
/// This receipt proves only the checks named by this module. It is not an XML
/// well-formedness, XSD-conformance, archive-safety, or executable-model claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InteropModelDescriptionPreflightV1 {
    standard: InterchangeStandardV1,
    content_hash: ContentHash,
    byte_len: usize,
    root_attribute_count: usize,
}

impl InteropModelDescriptionPreflightV1 {
    /// Exact declared standard whose root coordinate was checked.
    #[must_use]
    pub const fn standard(&self) -> InterchangeStandardV1 {
        self.standard
    }

    /// Hash of the exact raw bytes inspected.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    /// Exact bounded raw byte length inspected.
    #[must_use]
    pub const fn byte_len(&self) -> usize {
        self.byte_len
    }

    /// Number of unique attributes parsed from the root start tag.
    #[must_use]
    pub const fn root_attribute_count(&self) -> usize {
        self.root_attribute_count
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

/// Foreign-execution draft whose exact model-description bytes were preflighted.
///
/// The inner draft is deliberately private so its standard and content hash
/// cannot be rebound after the preflight receipt is minted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreflightedForeignExecutionDraftV1 {
    draft: ForeignExecutionDraftV1,
    preflight: InteropModelDescriptionPreflightV1,
}

impl PreflightedForeignExecutionDraftV1 {
    /// Retained bounded structural-preflight evidence.
    #[must_use]
    pub const fn preflight(&self) -> InteropModelDescriptionPreflightV1 {
        self.preflight
    }

    /// Admit the checked draft against one exact workflow and Machine graph.
    ///
    /// # Errors
    /// Returns the same graph, workflow, output, evidence, or identity refusal
    /// as [`ForeignExecutionDraftV1::admit_against`].
    pub fn admit_against(
        self,
        workflow: &AdmittedMachineWorkflowV1,
        graph: &AdmittedMachineGraph,
    ) -> Result<PreflightedAdmittedForeignExecutionV1, ForeignExecutionRefusalV1> {
        let execution = self.draft.admit_against(workflow, graph)?;
        Ok(PreflightedAdmittedForeignExecutionV1 {
            execution,
            preflight: self.preflight,
        })
    }
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

/// Admitted Estimated-only foreign execution retaining byte-preflight evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreflightedAdmittedForeignExecutionV1 {
    execution: AdmittedForeignExecutionV1,
    preflight: InteropModelDescriptionPreflightV1,
}

impl PreflightedAdmittedForeignExecutionV1 {
    /// The ordinary quarantined execution receipt.
    #[must_use]
    pub const fn execution(&self) -> &AdmittedForeignExecutionV1 {
        &self.execution
    }

    /// Retained bounded structural-preflight evidence for the exact model bytes.
    #[must_use]
    pub const fn preflight(&self) -> InteropModelDescriptionPreflightV1 {
        self.preflight
    }

    /// Explicitly discard the extra preflight evidence and retain the execution.
    #[must_use]
    pub fn into_execution(self) -> AdmittedForeignExecutionV1 {
        self.execution
    }
}

impl ForeignExecutionDraftV1 {
    /// Check exact bounded FMI/SSP model-description bytes and consume the draft.
    ///
    /// The preflight validates UTF-8 and a bounded root start tag, refuses DTD
    /// and entity declarations, enforces the standard-specific root/version/
    /// namespace/identity attributes, and checks the raw-byte content hash. It
    /// does not validate the complete XML document, XSD, FMU/SSP archive, or
    /// executable model.
    ///
    /// # Errors
    /// Refuses oversized or malformed input, a mismatched standard coordinate,
    /// unsafe declarations, or bytes not identified by `model_description`.
    pub fn preflight_model_description(
        self,
        bytes: &[u8],
    ) -> Result<PreflightedForeignExecutionDraftV1, ModelDescriptionPreflightRefusalV1> {
        let preflight = preflight_model_description_bytes(self.standard, bytes)?;
        if preflight.content_hash != self.model_description.content_hash() {
            return Err(ModelDescriptionPreflightRefusalV1::ContentHashMismatch {
                declared: self.model_description.content_hash(),
                actual: preflight.content_hash,
            });
        }
        Ok(PreflightedForeignExecutionDraftV1 {
            draft: self,
            preflight,
        })
    }

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

#[derive(Debug)]
struct XmlRootAttribute<'a> {
    name: &'a str,
    value: &'a str,
}

#[derive(Debug)]
struct XmlRootStartTag<'a> {
    name: &'a str,
    attributes: Vec<XmlRootAttribute<'a>>,
}

fn preflight_model_description_bytes(
    standard: InterchangeStandardV1,
    bytes: &[u8],
) -> Result<InteropModelDescriptionPreflightV1, ModelDescriptionPreflightRefusalV1> {
    if bytes.len() > MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1 {
        return Err(ModelDescriptionPreflightRefusalV1::TooLarge {
            actual: bytes.len(),
            max: MAX_INTEROP_MODEL_DESCRIPTION_BYTES_V1,
        });
    }
    let source = core::str::from_utf8(bytes).map_err(|error| {
        ModelDescriptionPreflightRefusalV1::InvalidUtf8 {
            valid_up_to: error.valid_up_to(),
        }
    })?;
    if let Some(offset) = bytes.iter().position(|byte| *byte == 0) {
        return Err(ModelDescriptionPreflightRefusalV1::NulByte { offset });
    }
    for (needle, kind) in [
        (b"<!DOCTYPE".as_slice(), "DOCTYPE"),
        (b"<!ENTITY".as_slice(), "ENTITY"),
    ] {
        if let Some(offset) = find_bytes(bytes, needle) {
            return Err(ModelDescriptionPreflightRefusalV1::ForbiddenDeclaration { kind, offset });
        }
    }

    let root = scan_xml_root_start_tag(source, standard)?;
    validate_model_description_root(standard, &root)?;
    Ok(InteropModelDescriptionPreflightV1 {
        standard,
        content_hash: hash_bytes(bytes),
        byte_len: bytes.len(),
        root_attribute_count: root.attributes.len(),
    })
}

fn scan_xml_root_start_tag(
    source: &str,
    standard: InterchangeStandardV1,
) -> Result<XmlRootStartTag<'_>, ModelDescriptionPreflightRefusalV1> {
    let bytes = source.as_bytes();
    let mut cursor = usize::from(bytes.starts_with(&[0xef, 0xbb, 0xbf])) * 3;
    let has_declaration = bytes.get(cursor..).is_some_and(|tail| {
        tail.starts_with(b"<?xml") && tail.get(5).is_some_and(|byte| is_xml_space(*byte))
    });

    if standard == InterchangeStandardV1::Fmi302 && !has_declaration {
        return Err(ModelDescriptionPreflightRefusalV1::MissingXmlDeclaration);
    }
    if has_declaration {
        let close = find_before_root_limit(bytes, cursor + 5, b"?>")?;
        if standard == InterchangeStandardV1::Fmi302
            && bytes[cursor..close]
                .iter()
                .any(|byte| matches!(byte, b'\r' | b'\n'))
        {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor,
                rule: "FMI XML declaration must remain on the first line",
            });
        }
        let declaration = parse_xml_attributes(source, cursor + 5, close)?;
        require_exact_attribute(&declaration, "version", "1.0")?;
        let encoding = declaration
            .iter()
            .find(|attribute| attribute.name == "encoding")
            .map(|attribute| attribute.value);
        if standard == InterchangeStandardV1::Fmi302 && encoding.is_none() {
            return Err(ModelDescriptionPreflightRefusalV1::MissingAttribute { name: "encoding" });
        }
        if let Some(encoding) = encoding.filter(|encoding| !encoding.eq_ignore_ascii_case("UTF-8"))
        {
            return Err(ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
                name: "encoding".to_owned(),
                expected: "UTF-8",
                actual: encoding.to_owned(),
            });
        }
        cursor = close + 2;
    }

    loop {
        skip_xml_space(bytes, &mut cursor);
        if cursor >= MAX_INTEROP_ROOT_PREFIX_BYTES_V1 {
            return Err(ModelDescriptionPreflightRefusalV1::RootPrefixTooLarge {
                max: MAX_INTEROP_ROOT_PREFIX_BYTES_V1,
            });
        }
        if bytes
            .get(cursor..)
            .is_some_and(|tail| tail.starts_with(b"<!--"))
        {
            let close = find_before_root_limit(bytes, cursor + 4, b"-->")?;
            cursor = close + 3;
            continue;
        }
        if bytes
            .get(cursor..)
            .is_some_and(|tail| tail.starts_with(b"<?"))
        {
            if bytes
                .get(cursor..)
                .is_some_and(|tail| tail.starts_with(b"<?xml"))
            {
                return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                    offset: cursor,
                    rule: "XML declaration must precede whitespace and processing instructions",
                });
            }
            let close = find_before_root_limit(bytes, cursor + 2, b"?>")?;
            cursor = close + 2;
            continue;
        }
        break;
    }

    if bytes.get(cursor) != Some(&b'<') {
        return Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: cursor.min(bytes.len()),
            rule: "expected root element start tag",
        });
    }
    let close = find_root_tag_close(bytes, cursor)?;
    let name_start = cursor + 1;
    let name_end = parse_xml_name_end(bytes, name_start, close)?;
    if name_end < close && !is_xml_space(bytes[name_end]) {
        return Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: name_end,
            rule: "expected whitespace after root element name",
        });
    }
    let attributes = parse_xml_attributes(source, name_end, close)?;
    Ok(XmlRootStartTag {
        name: &source[name_start..name_end],
        attributes,
    })
}

fn validate_model_description_root(
    standard: InterchangeStandardV1,
    root: &XmlRootStartTag<'_>,
) -> Result<(), ModelDescriptionPreflightRefusalV1> {
    match standard {
        InterchangeStandardV1::Fmi302 => {
            if root.name != standard.model_description_root() {
                return Err(ModelDescriptionPreflightRefusalV1::WrongRoot {
                    expected: standard.model_description_root(),
                    actual: root.name.to_owned(),
                });
            }
            if let Some(namespace) = root
                .attributes
                .iter()
                .find(|attribute| attribute.name == "xmlns")
                .map(|attribute| attribute.value)
                .filter(|namespace| !namespace.is_empty())
            {
                return Err(ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
                    name: "xmlns".to_owned(),
                    expected: "",
                    actual: namespace.to_owned(),
                });
            }
            require_exact_attribute(
                &root.attributes,
                "fmiVersion",
                standard.model_description_version(),
            )?;
            require_nonempty_attribute(&root.attributes, "modelName")?;
            require_nonempty_attribute(&root.attributes, "instantiationToken")?;
        }
        InterchangeStandardV1::Ssp20 => {
            let namespace_attribute = match root.name.split_once(':') {
                None if root.name == standard.model_description_root() => "xmlns".to_owned(),
                Some((prefix, local))
                    if !prefix.is_empty()
                        && !matches!(prefix, "xml" | "xmlns")
                        && !local.contains(':')
                        && local == standard.model_description_root() =>
                {
                    format!("xmlns:{prefix}")
                }
                _ => {
                    return Err(ModelDescriptionPreflightRefusalV1::WrongRoot {
                        expected: standard.model_description_root(),
                        actual: root.name.to_owned(),
                    });
                }
            };
            let namespace = root
                .attributes
                .iter()
                .find(|attribute| attribute.name == namespace_attribute)
                .map(|attribute| attribute.value)
                .ok_or(ModelDescriptionPreflightRefusalV1::MissingAttribute {
                    name: "root namespace declaration",
                })?;
            if namespace != SSP20_SYSTEM_STRUCTURE_NAMESPACE_V1 {
                return Err(ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
                    name: namespace_attribute,
                    expected: SSP20_SYSTEM_STRUCTURE_NAMESPACE_V1,
                    actual: namespace.to_owned(),
                });
            }
            require_exact_attribute(
                &root.attributes,
                "version",
                standard.model_description_version(),
            )?;
            require_nonempty_attribute(&root.attributes, "name")?;
        }
    }
    Ok(())
}

fn require_attribute<'a>(
    attributes: &[XmlRootAttribute<'a>],
    name: &'static str,
) -> Result<&'a str, ModelDescriptionPreflightRefusalV1> {
    attributes
        .iter()
        .find(|attribute| attribute.name == name)
        .map(|attribute| attribute.value)
        .ok_or(ModelDescriptionPreflightRefusalV1::MissingAttribute { name })
}

fn require_nonempty_attribute(
    attributes: &[XmlRootAttribute<'_>],
    name: &'static str,
) -> Result<(), ModelDescriptionPreflightRefusalV1> {
    if require_attribute(attributes, name)?.is_empty() {
        return Err(ModelDescriptionPreflightRefusalV1::EmptyAttribute { name });
    }
    Ok(())
}

fn require_exact_attribute(
    attributes: &[XmlRootAttribute<'_>],
    name: &'static str,
    expected: &'static str,
) -> Result<(), ModelDescriptionPreflightRefusalV1> {
    let actual = require_attribute(attributes, name)?;
    if actual != expected {
        return Err(ModelDescriptionPreflightRefusalV1::WrongAttributeValue {
            name: name.to_owned(),
            expected,
            actual: actual.to_owned(),
        });
    }
    Ok(())
}

fn parse_xml_attributes<'a>(
    source: &'a str,
    mut cursor: usize,
    end: usize,
) -> Result<Vec<XmlRootAttribute<'a>>, ModelDescriptionPreflightRefusalV1> {
    let bytes = source.as_bytes();
    let mut attributes = Vec::new();
    while cursor < end {
        let before_space = cursor;
        skip_xml_space_until(bytes, &mut cursor, end);
        if cursor == end {
            break;
        }
        if cursor == before_space {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor,
                rule: "expected whitespace before attribute",
            });
        }
        if bytes[cursor] == b'/' {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor,
                rule: "model-description root cannot be self-closing",
            });
        }
        let name_start = cursor;
        let name_end = parse_xml_name_end(bytes, name_start, end)?;
        cursor = name_end;
        skip_xml_space_until(bytes, &mut cursor, end);
        if bytes.get(cursor) != Some(&b'=') {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor.min(end),
                rule: "expected equals sign after attribute name",
            });
        }
        cursor += 1;
        skip_xml_space_until(bytes, &mut cursor, end);
        let quote = *bytes
            .get(cursor)
            .ok_or(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor.min(end),
                rule: "expected quoted attribute value",
            })?;
        if !matches!(quote, b'\'' | b'"') {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor,
                rule: "expected single- or double-quoted attribute value",
            });
        }
        cursor += 1;
        let value_start = cursor;
        while cursor < end && bytes[cursor] != quote {
            if bytes[cursor] == b'<' {
                return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                    offset: cursor,
                    rule: "raw less-than sign is forbidden in an attribute value",
                });
            }
            cursor += 1;
        }
        if cursor == end {
            return Err(ModelDescriptionPreflightRefusalV1::Malformed {
                offset: cursor,
                rule: "unterminated attribute value",
            });
        }
        let name = &source[name_start..name_end];
        if attributes
            .iter()
            .any(|attribute: &XmlRootAttribute<'_>| attribute.name == name)
        {
            return Err(ModelDescriptionPreflightRefusalV1::DuplicateAttribute {
                name: name.to_owned(),
            });
        }
        if attributes.len() == MAX_INTEROP_ROOT_ATTRIBUTES_V1 {
            return Err(ModelDescriptionPreflightRefusalV1::TooManyAttributes {
                actual: attributes.len() + 1,
                max: MAX_INTEROP_ROOT_ATTRIBUTES_V1,
            });
        }
        attributes.push(XmlRootAttribute {
            name,
            value: &source[value_start..cursor],
        });
        cursor += 1;
    }
    Ok(attributes)
}

fn parse_xml_name_end(
    bytes: &[u8],
    start: usize,
    end: usize,
) -> Result<usize, ModelDescriptionPreflightRefusalV1> {
    let Some(&first) = bytes.get(start).filter(|_| start < end) else {
        return Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: start.min(end),
            rule: "missing XML name",
        });
    };
    if !is_xml_name_start(first) {
        return Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: start,
            rule: "XML name must start with an ASCII letter, underscore, or colon",
        });
    }
    let mut cursor = start + 1;
    while cursor < end && is_xml_name_continue(bytes[cursor]) {
        cursor += 1;
    }
    Ok(cursor)
}

const fn is_xml_name_start(byte: u8) -> bool {
    byte.is_ascii_alphabetic() || matches!(byte, b'_' | b':')
}

const fn is_xml_name_continue(byte: u8) -> bool {
    is_xml_name_start(byte) || byte.is_ascii_digit() || matches!(byte, b'-' | b'.')
}

fn skip_xml_space(bytes: &[u8], cursor: &mut usize) {
    let limit = bytes.len().min(MAX_INTEROP_ROOT_PREFIX_BYTES_V1);
    while *cursor < limit && is_xml_space(bytes[*cursor]) {
        *cursor += 1;
    }
}

fn skip_xml_space_until(bytes: &[u8], cursor: &mut usize, end: usize) {
    while *cursor < end && is_xml_space(bytes[*cursor]) {
        *cursor += 1;
    }
}

const fn is_xml_space(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\r' | b'\n')
}

fn find_root_tag_close(
    bytes: &[u8],
    start: usize,
) -> Result<usize, ModelDescriptionPreflightRefusalV1> {
    let limit = bytes.len().min(MAX_INTEROP_ROOT_PREFIX_BYTES_V1);
    let mut cursor = start + 1;
    let mut quote = None;
    while cursor < limit {
        match (quote, bytes[cursor]) {
            (Some(open), byte) if byte == open => quote = None,
            (None, byte @ (b'\'' | b'"')) => quote = Some(byte),
            (None, b'>') => return Ok(cursor),
            _ => {}
        }
        cursor += 1;
    }
    if bytes.len() >= MAX_INTEROP_ROOT_PREFIX_BYTES_V1 {
        Err(ModelDescriptionPreflightRefusalV1::RootPrefixTooLarge {
            max: MAX_INTEROP_ROOT_PREFIX_BYTES_V1,
        })
    } else {
        Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: bytes.len(),
            rule: "unterminated root start tag",
        })
    }
}

fn find_before_root_limit(
    bytes: &[u8],
    start: usize,
    needle: &[u8],
) -> Result<usize, ModelDescriptionPreflightRefusalV1> {
    let limit = bytes.len().min(MAX_INTEROP_ROOT_PREFIX_BYTES_V1);
    if start <= limit
        && needle.len() <= limit.saturating_sub(start)
        && let Some(relative) = find_bytes(&bytes[start..limit], needle)
    {
        return Ok(start + relative);
    }
    if bytes.len() >= MAX_INTEROP_ROOT_PREFIX_BYTES_V1 {
        Err(ModelDescriptionPreflightRefusalV1::RootPrefixTooLarge {
            max: MAX_INTEROP_ROOT_PREFIX_BYTES_V1,
        })
    } else {
        Err(ModelDescriptionPreflightRefusalV1::Malformed {
            offset: bytes.len(),
            rule: "unterminated XML declaration, comment, or processing instruction",
        })
    }
}

fn find_bytes(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
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
