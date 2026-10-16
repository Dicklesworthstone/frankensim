//! One-way Machine-IR projection into scenario and external-domain artifacts.
//!
//! The admitted Machine graph, behavior, and assurance overlays intentionally
//! contain opaque model, value, motion, reset, distribution, evidence, and
//! policy references. This module therefore does not guess a concrete
//! scenario from those declarations. A caller must supply the already
//! materialized [`Scenario`], every externally owned domain artifact, and an
//! explicit crosswalk for every durable Machine identity. Admission validates
//! that package, canonicalizes it, and publishes a projection receipt.
//!
//! The direction is deliberately one way: Machine IR -> domain artifacts.
//! There is no reverse constructor, equivalence claim, or evidence-strength
//! transfer API here.

use core::fmt;
use core::hash::Hash;
use core::num::NonZeroU64;

use std::collections::BTreeSet;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EntityId, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};
use fs_blake3::{ContentHash, hash_bytes};
use fs_exec::Cx;
use fs_scenario::ir::{SCENARIO_IR_VERSION, parse_ir, write_ir};
use fs_scenario::{
    Scenario, ScenarioError, ValidationBudget, ValidationError, ValidationPlan, Violation,
};

use crate::{IR_VERSION, VERSION as FS_IR_VERSION};

use super::assurance::{
    AccountingWindowId, AdmittedMachineAssurance, ExperimentId, FaultId, FidelityRungId, HazardId,
    MachineAssuranceIdV1, SensorId,
};
use super::semantics::{AdmittedMachineBehavior, EventId, MachineBehaviorIdV1, ToleranceId};
use super::{
    AdmittedMachineGraph, ClockId, InterfaceId, MachineElementId, MachineGraphIdV1, MachineIdError,
    RelationId, SubsystemId,
};

/// Canonical schema version of the PR-5 Machine-domain projection.
pub const MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum external artifacts admitted by one projection.
pub const MAX_MACHINE_DOMAIN_ARTIFACTS: usize = 4_096;
/// Maximum complete Machine-source crosswalk rows admitted by one projection.
pub const MAX_MACHINE_DOMAIN_CROSSWALKS: usize = 65_536;
/// Maximum aggregate targets across all crosswalk rows.
pub const MAX_MACHINE_DOMAIN_TARGETS: usize = 131_072;
/// Maximum bytes in one external selector or scenario locator component.
pub const MAX_MACHINE_DOMAIN_SELECTOR_BYTES: usize = 4_096;
/// Maximum canonical manifest bytes retained by one admitted projection.
pub const MAX_MACHINE_DOMAIN_MANIFEST_BYTES: usize = 32 * 1_024 * 1_024;
/// Maximum framed manifest-plus-scenario bytes retained for replay transport.
pub const MAX_MACHINE_DOMAIN_PORTABLE_PAYLOAD_BYTES: usize = 128 * 1_024 * 1_024;
/// Hard canonical scenario-IR ceiling inherited from the current parser.
pub const MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES: usize = 16 * 1_024 * 1_024;
/// Maximum aggregate fixed-shape records admitted to the bounded legacy plan scan.
pub const MAX_MACHINE_DOMAIN_SCENARIO_RECORDS: usize =
    (MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES - 4_096) / 2_048;

const MACHINE_DOMAIN_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(64 * 1_024 * 1_024, 16 * 1_024 * 1_024, 16, 262_144, 131_072);
const PORTABLE_PAYLOAD_MAGIC: &[u8] = b"FSMACHINE-DOMAIN-LOWERING-V1\0";

/// Canonical schema marker for [`MachineDomainArtifactId`].
pub enum MachineDomainArtifactIdSchemaV1 {}

impl CanonicalSchema for MachineDomainArtifactIdSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.domain-artifact-id.v1";
    const NAME: &'static str = "machine-domain-artifact-id";
    const VERSION: u32 = MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "one locally durable handle for an exact externally owned domain artifact";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required("canonical-key", WireType::Utf8)];
}

/// Typed digest of one locally named external domain artifact.
pub type MachineDomainArtifactEntityIdV1 = EntityId<MachineDomainArtifactIdSchemaV1>;

/// Locally durable identity of one externally owned domain artifact.
#[derive(Clone)]
pub struct MachineDomainArtifactId {
    canonical_key: Box<str>,
    receipt: IdentityReceipt<MachineDomainArtifactEntityIdV1>,
}

impl MachineDomainArtifactId {
    /// Admit a canonical human-auditable artifact key.
    ///
    /// # Errors
    /// Refuses a noncanonical key or bounded identity-construction failure.
    pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
        let key = key.into();
        super::validate_canonical_key("machine-domain-artifact-id", &key)?;
        let receipt = CanonicalEncoder::<MachineDomainArtifactEntityIdV1, _>::new(
            super::MACHINE_IDENTITY_LIMITS,
            NeverCancel,
        )?
        .utf8(Field::new(0, "canonical-key"), &key)?
        .finish()?;
        Ok(Self {
            canonical_key: key.into_boxed_str(),
            receipt,
        })
    }

    /// Exact canonical key retained for diagnostics and manifests.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.canonical_key
    }

    /// Domain-separated durable identity.
    #[must_use]
    pub const fn identity(&self) -> MachineDomainArtifactEntityIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineDomainArtifactEntityIdV1> {
        self.receipt
    }
}

impl PartialEq for MachineDomainArtifactId {
    fn eq(&self, other: &Self) -> bool {
        self.identity() == other.identity()
    }
}

impl Eq for MachineDomainArtifactId {}

impl PartialOrd for MachineDomainArtifactId {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MachineDomainArtifactId {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.identity().cmp(&other.identity())
    }
}

impl core::hash::Hash for MachineDomainArtifactId {
    fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
        self.identity().hash(state);
    }
}

impl fmt::Display for MachineDomainArtifactId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.canonical_key)
    }
}

impl fmt::Debug for MachineDomainArtifactId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MachineDomainArtifactId")
            .field("canonical_key", &self.canonical_key)
            .field("identity", &self.identity())
            .finish()
    }
}

/// Closed coarse domain class for one external artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MachineDomainArtifactKindV1 {
    /// Motion/path semantics.
    Motion = 1,
    /// Logical or physical time semantics.
    Time = 2,
    /// Uncertainty, distribution, or correlation semantics.
    Uncertainty = 3,
    /// Controller, reset, or hybrid execution semantics.
    Control = 4,
    /// Geometry or representation semantics.
    Geometry = 5,
    /// Physics, material, or constitutive semantics.
    Physics = 6,
    /// Evidence, policy, or assurance semantics.
    Assurance = 7,
    /// Another explicitly namespaced externally owned domain.
    External = 255,
}

impl MachineDomainArtifactKindV1 {
    const fn tag(self) -> u8 {
        self as u8
    }
}

/// Structured refusal while constructing an exact versioned reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineDomainReferenceError {
    /// A namespace or selector was empty.
    Empty {
        /// Reference component role.
        role: &'static str,
    },
    /// A namespace or selector exceeded its public byte bound.
    TooLong {
        /// Reference component role.
        role: &'static str,
        /// Submitted bytes.
        bytes: usize,
        /// Maximum admitted bytes.
        max: usize,
    },
    /// A retained content digest was all zeroes.
    ZeroDigest {
        /// Reference role.
        role: &'static str,
    },
    /// The locally durable ID was malformed.
    Id(MachineIdError),
}

impl fmt::Display for MachineDomainReferenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty { role } => write!(formatter, "{role} must not be empty"),
            Self::TooLong { role, bytes, max } => {
                write!(formatter, "{role} has {bytes} bytes; maximum is {max}")
            }
            Self::ZeroDigest { role } => write!(formatter, "{role} digest must not be all zero"),
            Self::Id(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl core::error::Error for MachineDomainReferenceError {}

impl From<MachineIdError> for MachineDomainReferenceError {
    fn from(error: MachineIdError) -> Self {
        Self::Id(error)
    }
}

fn validate_reference_text(
    role: &'static str,
    value: &str,
) -> Result<(), MachineDomainReferenceError> {
    if value.is_empty() {
        return Err(MachineDomainReferenceError::Empty { role });
    }
    if value.len() > MAX_MACHINE_DOMAIN_SELECTOR_BYTES {
        return Err(MachineDomainReferenceError::TooLong {
            role,
            bytes: value.len(),
            max: MAX_MACHINE_DOMAIN_SELECTOR_BYTES,
        });
    }
    Ok(())
}

fn validate_reference_namespace(
    role: &'static str,
    value: &str,
) -> Result<(), MachineDomainReferenceError> {
    validate_reference_text(role, value)?;
    super::validate_canonical_key(role, value).map_err(MachineDomainReferenceError::Id)
}

fn validate_reference_digest(
    role: &'static str,
    digest: ContentHash,
) -> Result<(), MachineDomainReferenceError> {
    if digest.as_bytes().iter().all(|byte| *byte == 0) {
        Err(MachineDomainReferenceError::ZeroDigest { role })
    } else {
        Ok(())
    }
}

/// Exact versioned externally owned domain artifact.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MachineDomainArtifactRefV1 {
    id: MachineDomainArtifactId,
    kind: MachineDomainArtifactKindV1,
    namespace: Box<str>,
    schema_version: NonZeroU64,
    content_hash: ContentHash,
}

impl MachineDomainArtifactRefV1 {
    /// Construct an exact external artifact reference.
    ///
    /// # Errors
    /// Refuses an empty/oversized namespace or an all-zero content hash.
    pub fn new(
        id: MachineDomainArtifactId,
        kind: MachineDomainArtifactKindV1,
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        content_hash: ContentHash,
    ) -> Result<Self, MachineDomainReferenceError> {
        let namespace = namespace.into();
        validate_reference_namespace("machine-domain-artifact-namespace", &namespace)?;
        validate_reference_digest("machine-domain-artifact", content_hash)?;
        Ok(Self {
            id,
            kind,
            namespace: namespace.into_boxed_str(),
            schema_version,
            content_hash,
        })
    }

    /// Locally durable artifact handle.
    #[must_use]
    pub const fn id(&self) -> &MachineDomainArtifactId {
        &self.id
    }

    /// Coarse external domain class.
    #[must_use]
    pub const fn kind(&self) -> MachineDomainArtifactKindV1 {
        self.kind
    }

    /// Exact externally owned schema namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Exact nonzero external schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact content hash supplied by the owning domain.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    fn canonical_row_len(&self) -> Option<usize> {
        32usize
            .checked_add(8)?
            .checked_add(self.id.canonical_key().len())?
            .checked_add(1)?
            .checked_add(8)?
            .checked_add(self.namespace.len())?
            .checked_add(8)?
            .checked_add(32)
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        row.extend_from_slice(self.id.identity().as_bytes());
        append_bytes(row, self.id.canonical_key().as_bytes());
        row.push(self.kind.tag());
        append_bytes(row, self.namespace.as_bytes());
        row.extend_from_slice(&self.schema_version.get().to_le_bytes());
        row.extend_from_slice(self.content_hash.as_bytes());
    }
}

/// Exact versioned law explaining one Machine-to-domain mapping row.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MachineCrosswalkLawRefV1 {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    content_hash: ContentHash,
}

impl MachineCrosswalkLawRefV1 {
    /// Construct an exact mapping-law reference.
    ///
    /// # Errors
    /// Refuses an empty/oversized namespace or all-zero digest.
    pub fn new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        content_hash: ContentHash,
    ) -> Result<Self, MachineDomainReferenceError> {
        let namespace = namespace.into();
        validate_reference_namespace("machine-crosswalk-law-namespace", &namespace)?;
        validate_reference_digest("machine-crosswalk-law", content_hash)?;
        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            content_hash,
        })
    }

    /// Exact mapping-law namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Exact mapping-law schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact mapping-law artifact hash.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    fn canonical_row_len(&self) -> Option<usize> {
        8usize
            .checked_add(self.namespace.len())?
            .checked_add(8)?
            .checked_add(32)
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        append_bytes(row, self.namespace.as_bytes());
        row.extend_from_slice(&self.schema_version.get().to_le_bytes());
        row.extend_from_slice(self.content_hash.as_bytes());
    }
}

/// Exact versioned policy governing the complete projection.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MachineLoweringPolicyRefV1 {
    namespace: Box<str>,
    schema_version: NonZeroU64,
    content_hash: ContentHash,
}

impl MachineLoweringPolicyRefV1 {
    /// Construct an exact projection-policy reference.
    ///
    /// # Errors
    /// Refuses an empty/oversized namespace or all-zero digest.
    pub fn new(
        namespace: impl Into<String>,
        schema_version: NonZeroU64,
        content_hash: ContentHash,
    ) -> Result<Self, MachineDomainReferenceError> {
        let namespace = namespace.into();
        validate_reference_namespace("machine-lowering-policy-namespace", &namespace)?;
        validate_reference_digest("machine-lowering-policy", content_hash)?;
        Ok(Self {
            namespace: namespace.into_boxed_str(),
            schema_version,
            content_hash,
        })
    }

    /// Exact policy namespace.
    #[must_use]
    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    /// Exact policy schema version.
    #[must_use]
    pub const fn schema_version(&self) -> NonZeroU64 {
        self.schema_version
    }

    /// Exact policy artifact hash.
    #[must_use]
    pub const fn content_hash(&self) -> ContentHash {
        self.content_hash
    }

    fn canonical_row_len(&self) -> Option<usize> {
        8usize
            .checked_add(self.namespace.len())?
            .checked_add(8)?
            .checked_add(32)
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        append_bytes(row, self.namespace.as_bytes());
        row.extend_from_slice(&self.schema_version.get().to_le_bytes());
        row.extend_from_slice(self.content_hash.as_bytes());
    }
}

/// Typed stable Machine source that must appear exactly once in a crosswalk.
///
/// The three aggregate variants retain the complete graph/behavior/assurance
/// identity. They transitively bind declarations which intentionally have no
/// separate local ID, such as condition rows, material bindings, Context QoI
/// rows, and V&V receipts. Every independently durable local ID is also
/// enumerated explicitly; vector positions are never used.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineLoweringSourceV1 {
    /// Complete admitted Machine graph.
    Graph(MachineGraphIdV1),
    /// Complete admitted Machine behavior overlay.
    Behavior(MachineBehaviorIdV1),
    /// Complete admitted Machine assurance overlay.
    Assurance(MachineAssuranceIdV1),
    /// Named logical clock.
    Clock(ClockId),
    /// Named subsystem.
    Subsystem(SubsystemId),
    /// One of the six durable Machine element roles.
    Element(MachineElementId),
    /// Directed graph relation.
    Relation(RelationId),
    /// Role-oriented interface.
    Interface(InterfaceId),
    /// Guarded event/reset declaration.
    Event(EventId),
    /// Parameter-tolerance declaration.
    Tolerance(ToleranceId),
    /// Sensor declaration.
    Sensor(SensorId),
    /// Local experiment binding.
    Experiment(ExperimentId),
    /// Hazard declaration.
    Hazard(HazardId),
    /// Fault declaration.
    Fault(FaultId),
    /// Accounting-window declaration.
    AccountingWindow(AccountingWindowId),
    /// Fidelity-rung declaration.
    FidelityRung(FidelityRungId),
}

impl MachineLoweringSourceV1 {
    const fn canonical_row_len(&self) -> usize {
        match self {
            Self::Element(_) => 34,
            _ => 33,
        }
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        match self {
            Self::Graph(id) => {
                row.push(1);
                row.extend_from_slice(id.as_bytes());
            }
            Self::Behavior(id) => {
                row.push(2);
                row.extend_from_slice(id.as_bytes());
            }
            Self::Assurance(id) => {
                row.push(3);
                row.extend_from_slice(id.as_bytes());
            }
            Self::Clock(id) => {
                row.push(4);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Subsystem(id) => {
                row.push(5);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Element(id) => {
                row.push(6);
                row.push(id.kind().tag());
                row.extend_from_slice(&id.digest_bytes());
            }
            Self::Relation(id) => {
                row.push(7);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Interface(id) => {
                row.push(8);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Event(id) => {
                row.push(9);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Tolerance(id) => {
                row.push(10);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Sensor(id) => {
                row.push(11);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Experiment(id) => {
                row.push(12);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Hazard(id) => {
                row.push(13);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::Fault(id) => {
                row.push(14);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::AccountingWindow(id) => {
                row.push(15);
                row.extend_from_slice(id.identity().as_bytes());
            }
            Self::FidelityRung(id) => {
                row.push(16);
                row.extend_from_slice(id.identity().as_bytes());
            }
        }
    }
}

/// Stable semantic locator inside the concrete scenario artifact.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineScenarioLocatorV1 {
    /// Complete scenario artifact.
    Root,
    /// Implicit world frame.
    WorldFrame,
    /// Exact non-world frame ID and name.
    Frame {
        /// Numeric scenario frame identity.
        id: u32,
        /// Exact retained frame name.
        name: Box<str>,
    },
    /// Exact region name appearing in a BC or contact declaration.
    Region(Box<str>),
    /// Exact unique load-case name.
    LoadCase(Box<str>),
    /// Exact unique combination name.
    Combination(Box<str>),
    /// Exact unique stochastic-ensemble name.
    Ensemble(Box<str>),
    /// Exact ordered contact-region pair.
    ContactPair {
        /// First declared endpoint.
        region_a: Box<str>,
        /// Second declared endpoint.
        region_b: Box<str>,
    },
}

impl MachineScenarioLocatorV1 {
    fn canonical_row_len(&self) -> Option<usize> {
        match self {
            Self::Root | Self::WorldFrame => Some(1),
            Self::Frame { name, .. } => 1usize
                .checked_add(4)?
                .checked_add(8)?
                .checked_add(name.len()),
            Self::Region(name)
            | Self::LoadCase(name)
            | Self::Combination(name)
            | Self::Ensemble(name) => 1usize.checked_add(8)?.checked_add(name.len()),
            Self::ContactPair { region_a, region_b } => 1usize
                .checked_add(8)?
                .checked_add(region_a.len())?
                .checked_add(8)?
                .checked_add(region_b.len()),
        }
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        match self {
            Self::Root => row.push(1),
            Self::WorldFrame => row.push(2),
            Self::Frame { id, name } => {
                row.push(3);
                row.extend_from_slice(&id.to_le_bytes());
                append_bytes(row, name.as_bytes());
            }
            Self::Region(name) => {
                row.push(4);
                append_bytes(row, name.as_bytes());
            }
            Self::LoadCase(name) => {
                row.push(5);
                append_bytes(row, name.as_bytes());
            }
            Self::Combination(name) => {
                row.push(6);
                append_bytes(row, name.as_bytes());
            }
            Self::Ensemble(name) => {
                row.push(7);
                append_bytes(row, name.as_bytes());
            }
            Self::ContactPair { region_a, region_b } => {
                row.push(8);
                append_bytes(row, region_a.as_bytes());
                append_bytes(row, region_b.as_bytes());
            }
        }
    }
}

/// Exact target of one Machine crosswalk row.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineDomainTargetV1 {
    /// Stable semantic locator inside the supplied scenario.
    Scenario(MachineScenarioLocatorV1),
    /// Exact selector interpreted by an externally owned domain artifact.
    ExternalArtifact {
        /// Locally durable external artifact handle.
        artifact: MachineDomainArtifactId,
        /// Exact selector owned by that external artifact's schema.
        selector: Box<str>,
    },
}

impl MachineDomainTargetV1 {
    /// Construct a whole-artifact or sub-artifact external target.
    ///
    /// # Errors
    /// Refuses an empty or oversized selector.
    pub fn external(
        artifact: MachineDomainArtifactId,
        selector: impl Into<String>,
    ) -> Result<Self, MachineDomainReferenceError> {
        let selector = selector.into();
        validate_reference_text("machine-domain target selector", &selector)?;
        Ok(Self::ExternalArtifact {
            artifact,
            selector: selector.into_boxed_str(),
        })
    }

    fn canonical_row_len(&self) -> Option<usize> {
        match self {
            Self::Scenario(locator) => 1usize
                .checked_add(8)?
                .checked_add(locator.canonical_row_len()?),
            Self::ExternalArtifact { selector, .. } => 1usize
                .checked_add(32)?
                .checked_add(8)?
                .checked_add(selector.len()),
        }
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        match self {
            Self::Scenario(locator) => {
                row.push(1);
                append_length(
                    row,
                    locator
                        .canonical_row_len()
                        .expect("bounded scenario locator row length"),
                );
                locator.append_canonical_row(row);
            }
            Self::ExternalArtifact { artifact, selector } => {
                row.push(2);
                row.extend_from_slice(artifact.identity().as_bytes());
                append_bytes(row, selector.as_bytes());
            }
        }
    }
}

/// Construction refusal for one individual crosswalk row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineCrosswalkEntryError {
    /// The source had no declared target.
    EmptyTargets,
    /// The target set repeated an exact target.
    DuplicateTarget {
        /// Repeated exact target.
        target: MachineDomainTargetV1,
    },
    /// Aggregate target count exceeded its public bound.
    TargetLimit {
        /// Submitted target count.
        count: usize,
        /// Maximum admitted target count.
        max: usize,
    },
}

impl fmt::Display for MachineCrosswalkEntryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyTargets => formatter.write_str("Machine crosswalk source has no target"),
            Self::DuplicateTarget { target } => {
                write!(formatter, "Machine crosswalk repeats target {target:?}")
            }
            Self::TargetLimit { count, max } => write!(
                formatter,
                "Machine crosswalk has {count} targets; maximum is {max}"
            ),
        }
    }
}

impl core::error::Error for MachineCrosswalkEntryError {}

/// One complete stable-ID to domain-target mapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineDomainCrosswalkEntryV1 {
    source: MachineLoweringSourceV1,
    targets: Vec<MachineDomainTargetV1>,
    law: MachineCrosswalkLawRefV1,
}

impl MachineDomainCrosswalkEntryV1 {
    /// Construct and canonicalize one source mapping.
    ///
    /// # Errors
    /// Refuses an empty, duplicate, or oversized target set.
    pub fn new(
        source: MachineLoweringSourceV1,
        mut targets: Vec<MachineDomainTargetV1>,
        law: MachineCrosswalkLawRefV1,
    ) -> Result<Self, MachineCrosswalkEntryError> {
        if targets.is_empty() {
            return Err(MachineCrosswalkEntryError::EmptyTargets);
        }
        if targets.len() > MAX_MACHINE_DOMAIN_TARGETS {
            return Err(MachineCrosswalkEntryError::TargetLimit {
                count: targets.len(),
                max: MAX_MACHINE_DOMAIN_TARGETS,
            });
        }
        targets.sort();
        if let Some(pair) = targets.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(MachineCrosswalkEntryError::DuplicateTarget {
                target: pair[0].clone(),
            });
        }
        Ok(Self {
            source,
            targets,
            law,
        })
    }

    /// Exact typed Machine source.
    #[must_use]
    pub const fn source(&self) -> &MachineLoweringSourceV1 {
        &self.source
    }

    /// Canonically ordered nonempty target set.
    #[must_use]
    pub fn targets(&self) -> &[MachineDomainTargetV1] {
        &self.targets
    }

    /// Exact mapping law.
    #[must_use]
    pub const fn law(&self) -> &MachineCrosswalkLawRefV1 {
        &self.law
    }

    fn canonical_row_len(&self) -> Option<usize> {
        let mut len = 8usize
            .checked_add(self.source.canonical_row_len())?
            .checked_add(8)?;
        for target in &self.targets {
            len = len
                .checked_add(8)?
                .checked_add(target.canonical_row_len()?)?;
        }
        len.checked_add(8)?
            .checked_add(self.law.canonical_row_len()?)
    }

    fn append_canonical_row(&self, row: &mut Vec<u8>) {
        append_length(row, self.source.canonical_row_len());
        self.source.append_canonical_row(row);
        append_length(row, self.targets.len());
        for target in &self.targets {
            append_length(
                row,
                target
                    .canonical_row_len()
                    .expect("bounded crosswalk target row length"),
            );
            target.append_canonical_row(row);
        }
        append_length(
            row,
            self.law
                .canonical_row_len()
                .expect("bounded crosswalk law row length"),
        );
        self.law.append_canonical_row(row);
    }
}

/// Authority-free caller-supplied projection draft.
#[derive(Debug, Clone, PartialEq)]
pub struct MachineDomainLoweringDraftV1 {
    /// Concrete scenario materialized by an explicit caller-owned resolver.
    pub scenario: Scenario,
    /// Exact external motion/time/UQ/control/domain artifacts.
    pub external_artifacts: Vec<MachineDomainArtifactRefV1>,
    /// Complete mapping for every required durable Machine source.
    pub crosswalks: Vec<MachineDomainCrosswalkEntryV1>,
    /// Exact policy governing the complete projection.
    pub policy: MachineLoweringPolicyRefV1,
}

/// Submitted projection shape retained even when admission refuses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineDomainLoweringSubmittedCountsV1 {
    external_artifacts: usize,
    crosswalks: usize,
    crosswalk_targets: Option<usize>,
}

impl MachineDomainLoweringSubmittedCountsV1 {
    /// Number of exact external artifacts submitted.
    #[must_use]
    pub const fn external_artifacts(self) -> usize {
        self.external_artifacts
    }

    /// Number of stable-ID crosswalk rows submitted.
    #[must_use]
    pub const fn crosswalks(self) -> usize {
        self.crosswalks
    }

    /// Aggregate target count when its cancellation-polled scan completed.
    #[must_use]
    pub const fn crosswalk_targets(self) -> Option<usize> {
        self.crosswalk_targets
    }
}

fn submitted_counts(
    draft: &MachineDomainLoweringDraftV1,
) -> MachineDomainLoweringSubmittedCountsV1 {
    MachineDomainLoweringSubmittedCountsV1 {
        external_artifacts: draft.external_artifacts.len(),
        crosswalks: draft.crosswalks.len(),
        crosswalk_targets: None,
    }
}

fn submitted_target_count(
    draft: &MachineDomainLoweringDraftV1,
    cx: &Cx<'_>,
) -> Result<usize, MachineDomainLoweringRefusal> {
    checkpoint(cx, "initial")?;
    enforce_limit(
        "external artifact count",
        draft.external_artifacts.len(),
        MAX_MACHINE_DOMAIN_ARTIFACTS,
    )?;
    enforce_limit(
        "crosswalk count",
        draft.crosswalks.len(),
        MAX_MACHINE_DOMAIN_CROSSWALKS,
    )?;
    let mut targets = 0usize;
    for (index, entry) in draft.crosswalks.iter().enumerate() {
        poll_stride(cx, "submitted crosswalk-target counting", index)?;
        targets = targets.checked_add(entry.targets.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "crosswalk target count",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_TARGETS,
            },
        )?;
        enforce_limit(
            "crosswalk target count",
            targets,
            MAX_MACHINE_DOMAIN_TARGETS,
        )?;
    }
    Ok(targets)
}

/// Canonical identity schema for one admitted Machine-domain projection.
pub enum MachineDomainLoweringIdentitySchemaV1 {}

impl CanonicalSchema for MachineDomainLoweringIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.domain-lowering.v1";
    const NAME: &'static str = "admitted-machine-domain-lowering";
    const VERSION: u32 = MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "exact Machine graph, behavior, assurance, validated scenario, external artifacts, and complete one-way stable-ID crosswalk";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("lowering-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("machine-behavior", WireType::Bytes),
        FieldSpec::required("machine-assurance", WireType::Bytes),
        FieldSpec::required("scenario-ir-version", WireType::U64),
        FieldSpec::required("scenario-ir-hash", WireType::Bytes),
        FieldSpec::required("scenario-validation-budget", WireType::Bytes),
        FieldSpec::required("scenario-validation-plan", WireType::Bytes),
        FieldSpec::required("lowering-policy", WireType::Bytes),
        FieldSpec::required("external-artifacts", WireType::OrderedBytes),
        FieldSpec::required("crosswalks", WireType::OrderedBytes),
        FieldSpec::required("canonical-manifest-hash", WireType::Bytes),
    ];
}

/// Strong semantic identity of one admitted one-way projection.
pub type MachineDomainLoweringIdV1 = ProblemSemanticId<MachineDomainLoweringIdentitySchemaV1>;

/// Exact reason a scenario locator failed to resolve uniquely.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MachineScenarioLocatorProblemV1 {
    /// A text component was empty or exceeded the public bound.
    InvalidText,
    /// No declared scenario object matched.
    Missing,
    /// More than one scenario object matched a locator requiring uniqueness.
    Ambiguous,
}

/// Typed fail-closed refusal from projection admission.
#[derive(Debug, Clone, PartialEq)]
pub enum MachineDomainLoweringRefusal {
    /// Behavior was admitted against another graph.
    BehaviorGraphMismatch,
    /// Assurance was admitted against another graph.
    AssuranceGraphMismatch,
    /// Assurance was admitted against another behavior.
    AssuranceBehaviorMismatch,
    /// Cancellation was observed before publication.
    Cancelled {
        /// Deterministic admission phase.
        phase: &'static str,
    },
    /// One public or aggregate resource exceeded its bound.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Submitted amount.
        requested: usize,
        /// Maximum admitted amount.
        limit: usize,
    },
    /// A bounded scratch allocation was refused.
    AllocationRefused {
        /// Stable scratch-resource name.
        resource: &'static str,
        /// Requested elements.
        requested: usize,
    },
    /// An external artifact ID was declared more than once.
    DuplicateArtifact {
        /// Repeated local artifact ID.
        artifact: MachineDomainArtifactId,
    },
    /// A required Machine source had more than one crosswalk row.
    DuplicateCrosswalk {
        /// Repeated source.
        source: MachineLoweringSourceV1,
    },
    /// One exact admitted source had no crosswalk row.
    MissingCrosswalk {
        /// Unrepresented required source.
        source: MachineLoweringSourceV1,
    },
    /// A crosswalk row named a source outside this admitted stack.
    UnexpectedCrosswalk {
        /// Foreign or otherwise unexpected source.
        source: MachineLoweringSourceV1,
    },
    /// A scenario locator failed exact resolution.
    ScenarioTarget {
        /// Exact failed locator.
        locator: MachineScenarioLocatorV1,
        /// Stable failure class.
        problem: MachineScenarioLocatorProblemV1,
    },
    /// A crosswalk named an undeclared external artifact.
    UnknownExternalTarget {
        /// Missing artifact handle.
        artifact: MachineDomainArtifactId,
    },
    /// A directly constructed external target used an empty selector.
    EmptyExternalSelector {
        /// Exact artifact whose selector was empty.
        artifact: MachineDomainArtifactId,
    },
    /// An external artifact was not reachable from any crosswalk row.
    OrphanExternalArtifact {
        /// Unused external artifact handle.
        artifact: MachineDomainArtifactId,
    },
    /// Scenario validation could not be admitted under the explicit budget.
    ScenarioValidation(ValidationError),
    /// Scenario semantic validation produced one or more findings.
    ScenarioFindings(Vec<Violation>),
    /// The canonical scenario writer emitted bytes its parser refused.
    ScenarioParse(ScenarioError),
    /// Canonical scenario bytes reported a non-current version.
    ScenarioVersion {
        /// Version found by the parser.
        found: u32,
        /// Required current version.
        expected: u32,
    },
    /// Canonical current-version bytes unexpectedly required migration.
    UnexpectedScenarioMigration,
    /// Parse/reprint changed the scenario value or canonical bytes.
    ScenarioRoundTripMismatch,
    /// Canonical identity publication failed.
    Identity(CanonicalError),
}

impl MachineDomainLoweringRefusal {
    /// Stable machine-readable top-level rule code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::BehaviorGraphMismatch => "MachineLoweringBehaviorGraphMismatch",
            Self::AssuranceGraphMismatch => "MachineLoweringAssuranceGraphMismatch",
            Self::AssuranceBehaviorMismatch => "MachineLoweringAssuranceBehaviorMismatch",
            Self::Cancelled { .. } => "MachineLoweringCancelled",
            Self::ResourceLimit { .. } => "MachineLoweringResourceLimit",
            Self::AllocationRefused { .. } => "MachineLoweringAllocationRefused",
            Self::DuplicateArtifact { .. } => "MachineLoweringDuplicateArtifact",
            Self::DuplicateCrosswalk { .. } => "MachineLoweringDuplicateCrosswalk",
            Self::MissingCrosswalk { .. } => "MachineLoweringMissingCrosswalk",
            Self::UnexpectedCrosswalk { .. } => "MachineLoweringUnexpectedCrosswalk",
            Self::ScenarioTarget { .. } => "MachineLoweringScenarioTarget",
            Self::UnknownExternalTarget { .. } => "MachineLoweringUnknownExternalTarget",
            Self::EmptyExternalSelector { .. } => "MachineLoweringEmptyExternalSelector",
            Self::OrphanExternalArtifact { .. } => "MachineLoweringOrphanExternalArtifact",
            Self::ScenarioValidation(_) => "MachineLoweringScenarioValidation",
            Self::ScenarioFindings(_) => "MachineLoweringScenarioFindings",
            Self::ScenarioParse(_) => "MachineLoweringScenarioParse",
            Self::ScenarioVersion { .. } => "MachineLoweringScenarioVersion",
            Self::UnexpectedScenarioMigration => "MachineLoweringUnexpectedScenarioMigration",
            Self::ScenarioRoundTripMismatch => "MachineLoweringScenarioRoundTripMismatch",
            Self::Identity(_) => "MachineLoweringIdentity",
        }
    }
}

impl fmt::Display for MachineDomainLoweringRefusal {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::BehaviorGraphMismatch => {
                formatter.write_str("Machine behavior belongs to another graph")
            }
            Self::AssuranceGraphMismatch => {
                formatter.write_str("Machine assurance belongs to another graph")
            }
            Self::AssuranceBehaviorMismatch => {
                formatter.write_str("Machine assurance belongs to another behavior")
            }
            Self::Cancelled { phase } => {
                write!(formatter, "Machine lowering cancelled during {phase}")
            }
            Self::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "Machine lowering {resource} request {requested} exceeds limit {limit}"
            ),
            Self::AllocationRefused {
                resource,
                requested,
            } => write!(
                formatter,
                "Machine lowering could not allocate {requested} {resource} elements"
            ),
            Self::DuplicateArtifact { artifact } => {
                write!(formatter, "duplicate external artifact {artifact}")
            }
            Self::DuplicateCrosswalk { source } => {
                write!(formatter, "duplicate Machine crosswalk source {source:?}")
            }
            Self::MissingCrosswalk { source } => {
                write!(formatter, "missing Machine crosswalk source {source:?}")
            }
            Self::UnexpectedCrosswalk { source } => {
                write!(formatter, "unexpected Machine crosswalk source {source:?}")
            }
            Self::ScenarioTarget { locator, problem } => {
                write!(
                    formatter,
                    "scenario locator {locator:?} failed: {problem:?}"
                )
            }
            Self::UnknownExternalTarget { artifact } => {
                write!(
                    formatter,
                    "crosswalk names undeclared external artifact {artifact}"
                )
            }
            Self::EmptyExternalSelector { artifact } => {
                write!(
                    formatter,
                    "external target {artifact} has an empty selector"
                )
            }
            Self::OrphanExternalArtifact { artifact } => {
                write!(
                    formatter,
                    "external artifact {artifact} has no crosswalk target"
                )
            }
            Self::ScenarioValidation(error) => fmt::Display::fmt(error, formatter),
            Self::ScenarioFindings(findings) => write!(
                formatter,
                "scenario lowering refused with {} semantic finding(s)",
                findings.len()
            ),
            Self::ScenarioParse(error) => fmt::Display::fmt(error, formatter),
            Self::ScenarioVersion { found, expected } => write!(
                formatter,
                "scenario IR version {found} does not match required {expected}"
            ),
            Self::UnexpectedScenarioMigration => {
                formatter.write_str("current scenario IR unexpectedly required migration")
            }
            Self::ScenarioRoundTripMismatch => {
                formatter.write_str("scenario canonical parse/reprint was not byte-stable")
            }
            Self::Identity(error) => fmt::Display::fmt(error, formatter),
        }
    }
}

impl core::error::Error for MachineDomainLoweringRefusal {}

impl From<CanonicalError> for MachineDomainLoweringRefusal {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

/// Sealed, canonically ordered one-way Machine-domain projection.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedMachineDomainLoweringV1 {
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    assurance: MachineAssuranceIdV1,
    scenario: Scenario,
    canonical_scenario_ir: Box<str>,
    scenario_hash: ContentHash,
    validation_budget: ValidationBudget,
    validation_plan: ValidationPlan,
    external_artifacts: Vec<MachineDomainArtifactRefV1>,
    crosswalks: Vec<MachineDomainCrosswalkEntryV1>,
    policy: MachineLoweringPolicyRefV1,
    canonical_manifest: Box<[u8]>,
    portable_payload: Box<[u8]>,
    receipt: IdentityReceipt<MachineDomainLoweringIdV1>,
}

/// Bounded structured outcome for one Machine-domain projection attempt.
///
/// This record is suitable for structured tracing. It is not a canonical
/// digest or replay record for an early-refused caller draft.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub struct MachineDomainLoweringAdmissionDecisionV1 {
    submitted: MachineDomainLoweringSubmittedCountsV1,
    result: Result<AdmittedMachineDomainLoweringV1, MachineDomainLoweringRefusal>,
}

impl MachineDomainLoweringAdmissionDecisionV1 {
    /// Exact bounded collection counts observed before canonicalization.
    #[must_use]
    pub const fn submitted_counts(&self) -> MachineDomainLoweringSubmittedCountsV1 {
        self.submitted
    }

    /// Stable admitted/refused top-level outcome code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "MachineDomainLoweringAdmitted",
            Err(_) => "MachineDomainLoweringRefused",
        }
    }

    /// Stable detailed refusal code, if admission refused.
    #[must_use]
    pub fn refusal_code(&self) -> Option<&'static str> {
        self.result.as_ref().err().map(|refusal| refusal.code())
    }

    /// Borrow the admitted projection or complete typed refusal.
    #[must_use]
    pub fn result(
        &self,
    ) -> Result<&AdmittedMachineDomainLoweringV1, &MachineDomainLoweringRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision and recover the conventional result.
    #[must_use]
    pub fn into_result(
        self,
    ) -> Result<AdmittedMachineDomainLoweringV1, MachineDomainLoweringRefusal> {
        self.result
    }
}

impl AdmittedMachineDomainLoweringV1 {
    /// Exact source graph identity.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact source behavior identity.
    #[must_use]
    pub const fn behavior(&self) -> MachineBehaviorIdV1 {
        self.behavior
    }

    /// Exact source assurance identity.
    #[must_use]
    pub const fn assurance(&self) -> MachineAssuranceIdV1 {
        self.assurance
    }

    /// Strong identity of the complete projection.
    #[must_use]
    pub const fn identity(&self) -> MachineDomainLoweringIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineDomainLoweringIdV1> {
        self.receipt
    }

    /// Validated concrete scenario retained by the projection.
    #[must_use]
    pub const fn scenario(&self) -> &Scenario {
        &self.scenario
    }

    /// Byte-stable current-version scenario IR.
    #[must_use]
    pub fn canonical_scenario_ir(&self) -> &str {
        &self.canonical_scenario_ir
    }

    /// Plain content hash of the exact canonical scenario bytes.
    #[must_use]
    pub const fn scenario_hash(&self) -> ContentHash {
        self.scenario_hash
    }

    /// Explicit semantic-validation budget bound by projection identity.
    #[must_use]
    pub const fn validation_budget(&self) -> ValidationBudget {
        self.validation_budget
    }

    /// Preflighted scenario validation plan bound by projection identity.
    #[must_use]
    pub const fn validation_plan(&self) -> ValidationPlan {
        self.validation_plan
    }

    /// Canonically ordered exact external domain artifacts.
    #[must_use]
    pub fn external_artifacts(&self) -> &[MachineDomainArtifactRefV1] {
        &self.external_artifacts
    }

    /// Canonically ordered complete stable-ID crosswalk.
    #[must_use]
    pub fn crosswalks(&self) -> &[MachineDomainCrosswalkEntryV1] {
        &self.crosswalks
    }

    /// Exact projection policy.
    #[must_use]
    pub const fn policy(&self) -> &MachineLoweringPolicyRefV1 {
        &self.policy
    }

    /// Canonical manifest containing source IDs, artifact hashes, target rows,
    /// validation policy, and the canonical scenario hash.
    #[must_use]
    pub fn canonical_manifest(&self) -> &[u8] {
        &self.canonical_manifest
    }

    /// Portable exact replay payload: framed manifest plus canonical scenario
    /// IR. This is suitable as a package semantic witness when it fits the
    /// package's independent payload budget.
    #[must_use]
    pub fn portable_payload(&self) -> &[u8] {
        &self.portable_payload
    }

    /// Plain content hash of the exact portable payload.
    #[must_use]
    pub fn portable_payload_hash(&self) -> ContentHash {
        hash_bytes(&self.portable_payload)
    }

    /// Re-parse and re-derive every retained deterministic projection surface.
    ///
    /// This proves exact replay consistency only. It does not authenticate an
    /// external mapping law or strengthen scientific evidence.
    ///
    /// # Errors
    /// Returns a typed refusal if scenario, manifest, payload, or projection
    /// identity drifted from the admitted value.
    pub fn verify_replay(&self) -> Result<MachineDomainReplayReceiptV1, MachineDomainReplayError> {
        let decoded = parse_ir(&self.canonical_scenario_ir)
            .map_err(MachineDomainReplayError::ScenarioParse)?;
        if decoded.source_version() != SCENARIO_IR_VERSION || decoded.migration().is_some() {
            return Err(MachineDomainReplayError::ScenarioVersion);
        }
        if decoded.scenario() != &self.scenario
            || write_ir(decoded.scenario()) != self.canonical_scenario_ir.as_ref()
            || hash_bytes(self.canonical_scenario_ir.as_bytes()) != self.scenario_hash
        {
            return Err(MachineDomainReplayError::ScenarioDrift);
        }

        let artifact_rows = canonical_artifact_rows(&self.external_artifacts, None)
            .map_err(|_| MachineDomainReplayError::ManifestDrift)?;
        let crosswalk_rows = canonical_crosswalk_rows(&self.crosswalks, None)
            .map_err(|_| MachineDomainReplayError::ManifestDrift)?;
        let budget_row = validation_budget_row(self.validation_budget);
        let plan_row = validation_plan_row(self.validation_plan);
        let policy_row = canonical_policy_row(&self.policy)
            .map_err(|_| MachineDomainReplayError::ManifestDrift)?;
        let manifest_parts = ManifestParts {
            graph: self.graph,
            behavior: self.behavior,
            assurance: self.assurance,
            scenario_hash: self.scenario_hash,
            validation_budget: &budget_row,
            validation_plan: &plan_row,
            policy: &policy_row,
            artifacts: &artifact_rows,
            crosswalks: &crosswalk_rows,
        };
        let manifest = build_manifest(
            &manifest_parts,
            manifest_encoded_len(&manifest_parts).ok_or(MachineDomainReplayError::ManifestDrift)?,
        )
        .map_err(|_| MachineDomainReplayError::ManifestDrift)?;
        if manifest.as_slice() != self.canonical_manifest.as_ref() {
            return Err(MachineDomainReplayError::ManifestDrift);
        }
        let payload_len =
            portable_payload_encoded_len(manifest.len(), self.canonical_scenario_ir.len())
                .ok_or(MachineDomainReplayError::PayloadDrift)?;
        let payload = build_portable_payload(
            &manifest,
            self.canonical_scenario_ir.as_bytes(),
            payload_len,
        )
        .map_err(|_| MachineDomainReplayError::PayloadDrift)?;
        if payload.as_slice() != self.portable_payload.as_ref() {
            return Err(MachineDomainReplayError::PayloadDrift);
        }
        let receipt = encode_lowering_identity(
            self.graph,
            self.behavior,
            self.assurance,
            self.scenario_hash,
            &budget_row,
            &plan_row,
            &policy_row,
            &artifact_rows,
            &crosswalk_rows,
            hash_bytes(&manifest),
            NeverCancel,
        )
        .map_err(MachineDomainReplayError::Identity)?;
        if receipt.id() != self.identity()
            || receipt.canonical_preimage() != self.receipt.canonical_preimage()
        {
            return Err(MachineDomainReplayError::IdentityDrift);
        }

        Ok(MachineDomainReplayReceiptV1 {
            lowering: self.identity(),
            scenario_hash: self.scenario_hash,
            manifest_hash: hash_bytes(&manifest),
            portable_payload_hash: hash_bytes(&payload),
        })
    }
}

/// Successful exact replay-consistency receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineDomainReplayReceiptV1 {
    lowering: MachineDomainLoweringIdV1,
    scenario_hash: ContentHash,
    manifest_hash: ContentHash,
    portable_payload_hash: ContentHash,
}

impl MachineDomainReplayReceiptV1 {
    /// Replayed projection identity.
    #[must_use]
    pub const fn lowering(&self) -> MachineDomainLoweringIdV1 {
        self.lowering
    }

    /// Replayed canonical scenario hash.
    #[must_use]
    pub const fn scenario_hash(&self) -> ContentHash {
        self.scenario_hash
    }

    /// Replayed canonical manifest hash.
    #[must_use]
    pub const fn manifest_hash(&self) -> ContentHash {
        self.manifest_hash
    }

    /// Replayed portable payload hash.
    #[must_use]
    pub const fn portable_payload_hash(&self) -> ContentHash {
        self.portable_payload_hash
    }
}

/// Typed exact-replay failure.
#[derive(Debug, Clone, PartialEq)]
pub enum MachineDomainReplayError {
    /// Retained scenario bytes no longer parse.
    ScenarioParse(ScenarioError),
    /// Retained scenario bytes are not current-version unmigrated IR.
    ScenarioVersion,
    /// Retained scenario value, bytes, or hash changed.
    ScenarioDrift,
    /// Canonical manifest no longer re-derives exactly.
    ManifestDrift,
    /// Framed portable payload no longer re-derives exactly.
    PayloadDrift,
    /// Canonical identity could not be re-derived.
    Identity(CanonicalError),
    /// Re-derived typed identity or canonical preimage changed.
    IdentityDrift,
}

impl fmt::Display for MachineDomainReplayError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ScenarioParse(error) => fmt::Display::fmt(error, formatter),
            Self::ScenarioVersion => formatter.write_str("replay scenario version drifted"),
            Self::ScenarioDrift => formatter.write_str("replay scenario content drifted"),
            Self::ManifestDrift => formatter.write_str("replay manifest drifted"),
            Self::PayloadDrift => formatter.write_str("replay portable payload drifted"),
            Self::Identity(error) => fmt::Display::fmt(error, formatter),
            Self::IdentityDrift => formatter.write_str("replay projection identity drifted"),
        }
    }
}

impl core::error::Error for MachineDomainReplayError {}

/// Borrowed, structurally verified view of a portable PR-5 replay payload.
///
/// This decoder exposes the exact typed source identities and scenario/framed-
/// manifest hashes needed by a package or ledger consumer without
/// reconstructing authority-bearing external artifact objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineDomainPortablePayloadViewV1<'a> {
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    assurance: MachineAssuranceIdV1,
    manifest_bytes: &'a [u8],
    canonical_scenario_ir: &'a str,
    scenario_hash: ContentHash,
    manifest_hash: ContentHash,
    payload_hash: ContentHash,
    external_artifact_rows: usize,
    crosswalk_rows: usize,
}

impl<'a> MachineDomainPortablePayloadViewV1<'a> {
    /// Exact source graph identity decoded under its nominal schema.
    #[must_use]
    pub const fn graph(&self) -> MachineGraphIdV1 {
        self.graph
    }

    /// Exact source behavior identity decoded under its nominal schema.
    #[must_use]
    pub const fn behavior(&self) -> MachineBehaviorIdV1 {
        self.behavior
    }

    /// Exact source assurance identity decoded under its nominal schema.
    #[must_use]
    pub const fn assurance(&self) -> MachineAssuranceIdV1 {
        self.assurance
    }

    /// Exact opaque, structurally framed manifest bytes carried by the payload.
    ///
    /// The producer's canonical encoding is preserved byte-for-byte, but this
    /// structural decoder does not interpret or authenticate opaque row bodies.
    #[must_use]
    pub const fn manifest_bytes(&self) -> &'a [u8] {
        self.manifest_bytes
    }

    /// Exact current-version canonical scenario IR carried by the payload.
    #[must_use]
    pub const fn canonical_scenario_ir(&self) -> &'a str {
        self.canonical_scenario_ir
    }

    /// Plain content hash of the canonical scenario bytes.
    #[must_use]
    pub const fn scenario_hash(&self) -> ContentHash {
        self.scenario_hash
    }

    /// Plain content hash of the complete framed manifest bytes.
    #[must_use]
    pub const fn manifest_hash(&self) -> ContentHash {
        self.manifest_hash
    }

    /// Plain content hash of the complete framed payload.
    #[must_use]
    pub const fn payload_hash(&self) -> ContentHash {
        self.payload_hash
    }

    /// Number of exact external-artifact rows retained in the manifest.
    #[must_use]
    pub const fn external_artifact_rows(&self) -> usize {
        self.external_artifact_rows
    }

    /// Number of exact stable-ID crosswalk rows retained in the manifest.
    #[must_use]
    pub const fn crosswalk_rows(&self) -> usize {
        self.crosswalk_rows
    }
}

/// Typed structural refusal from portable replay-payload decoding.
#[derive(Debug, Clone, PartialEq)]
pub enum MachineDomainPortablePayloadError {
    /// Payload or nested field exceeded its hard byte/count bound.
    ResourceLimit {
        /// Stable field/resource name.
        resource: &'static str,
        /// Submitted amount.
        requested: usize,
        /// Maximum admitted amount.
        limit: usize,
    },
    /// A framed field ended before its declared size.
    Truncated {
        /// Stable field name.
        field: &'static str,
    },
    /// A frame had the wrong fixed magic bytes.
    Magic {
        /// Outer payload or nested manifest.
        field: &'static str,
    },
    /// A retained schema/language/scenario version is unsupported.
    Version {
        /// Versioned field.
        field: &'static str,
        /// Version found on the wire.
        found: u32,
        /// Version accepted by this decoder.
        expected: u32,
    },
    /// A retained crate-version stamp does not match this decoder build.
    BuildVersion {
        /// Crate-version field.
        field: &'static str,
    },
    /// A typed 32-byte source identity was malformed.
    Identity {
        /// Identity field.
        field: &'static str,
    },
    /// A text field was not valid UTF-8.
    Utf8 {
        /// Text field.
        field: &'static str,
    },
    /// A fixed-shape manifest field had the wrong byte count.
    FieldShape {
        /// Malformed field.
        field: &'static str,
    },
    /// A nested scenario failed canonical parsing.
    ScenarioParse(ScenarioError),
    /// Scenario bytes were not current-version, migration-free, and byte-stable.
    ScenarioRoundTrip,
    /// Retained scenario bytes did not match the manifest hash.
    ScenarioHash,
    /// One frame retained bytes after its complete declared shape.
    TrailingBytes {
        /// Frame with trailing bytes.
        field: &'static str,
    },
}

impl fmt::Display for MachineDomainPortablePayloadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                formatter,
                "portable Machine-domain {resource} request {requested} exceeds {limit}"
            ),
            Self::Truncated { field } => write!(formatter, "portable payload truncated at {field}"),
            Self::Magic { field } => {
                write!(formatter, "portable payload has invalid {field} magic")
            }
            Self::Version {
                field,
                found,
                expected,
            } => write!(
                formatter,
                "portable payload {field} version {found} does not match {expected}"
            ),
            Self::BuildVersion { field } => {
                write!(
                    formatter,
                    "portable payload {field} build version does not match"
                )
            }
            Self::Identity { field } => write!(formatter, "portable payload has invalid {field}"),
            Self::Utf8 { field } => write!(formatter, "portable payload {field} is not UTF-8"),
            Self::FieldShape { field } => {
                write!(
                    formatter,
                    "portable payload {field} has the wrong byte shape"
                )
            }
            Self::ScenarioParse(error) => fmt::Display::fmt(error, formatter),
            Self::ScenarioRoundTrip => {
                formatter.write_str("portable scenario is not canonical current-version IR")
            }
            Self::ScenarioHash => formatter.write_str("portable scenario hash does not match"),
            Self::TrailingBytes { field } => {
                write!(formatter, "portable payload {field} has trailing bytes")
            }
        }
    }
}

impl core::error::Error for MachineDomainPortablePayloadError {}

/// Decode and structurally verify an exact portable PR-5 replay payload.
///
/// This boundary independently checks framing, hard resource caps, build and
/// schema versions, nominal Machine identity widths, complete manifest row
/// framing, canonical scenario parse/reprint, and scenario hash binding. It
/// proves transport integrity only; external law semantics remain owned by
/// their named domains.
///
/// # Errors
/// Returns a typed refusal for malformed framing, resource excess, version or
/// identity mismatch, scenario drift, or trailing bytes.
#[allow(clippy::too_many_lines)]
pub fn parse_machine_domain_portable_payload_v1(
    bytes: &[u8],
) -> Result<MachineDomainPortablePayloadViewV1<'_>, MachineDomainPortablePayloadError> {
    enforce_portable_limit(
        "payload bytes",
        bytes.len(),
        MAX_MACHINE_DOMAIN_PORTABLE_PAYLOAD_BYTES,
    )?;
    let mut outer = PayloadCursor::new(bytes);
    outer.expect_magic("payload")?;
    let manifest = outer.length_prefixed("canonical manifest")?;
    enforce_portable_limit(
        "canonical manifest bytes",
        manifest.len(),
        MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
    )?;
    let scenario_bytes = outer.length_prefixed("canonical scenario IR")?;
    enforce_portable_limit(
        "canonical scenario IR bytes",
        scenario_bytes.len(),
        MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES,
    )?;
    outer.finish("payload")?;

    let mut cursor = PayloadCursor::new(manifest);
    cursor.expect_magic("manifest")?;
    cursor.expect_version("lowering schema", MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1)?;
    cursor.expect_version("FrankenScript IR", IR_VERSION)?;
    if cursor.length_prefixed("fs-ir version")? != FS_IR_VERSION.as_bytes() {
        return Err(MachineDomainPortablePayloadError::BuildVersion { field: "fs-ir" });
    }
    if cursor.length_prefixed("fs-scenario version")? != fs_scenario::VERSION.as_bytes() {
        return Err(MachineDomainPortablePayloadError::BuildVersion {
            field: "fs-scenario",
        });
    }
    let graph = MachineGraphIdV1::parse_slice(cursor.take("machine graph", 32)?).ok_or(
        MachineDomainPortablePayloadError::Identity {
            field: "machine graph",
        },
    )?;
    let behavior = MachineBehaviorIdV1::parse_slice(cursor.take("machine behavior", 32)?).ok_or(
        MachineDomainPortablePayloadError::Identity {
            field: "machine behavior",
        },
    )?;
    let assurance = MachineAssuranceIdV1::parse_slice(cursor.take("machine assurance", 32)?)
        .ok_or(MachineDomainPortablePayloadError::Identity {
            field: "machine assurance",
        })?;
    cursor.expect_version("scenario IR", SCENARIO_IR_VERSION)?;
    let scenario_hash = ContentHash::from_slice(cursor.take("scenario hash", 32)?).ok_or(
        MachineDomainPortablePayloadError::FieldShape {
            field: "scenario hash",
        },
    )?;
    if cursor.length_prefixed("validation budget")?.len() != 14 * 16 {
        return Err(MachineDomainPortablePayloadError::FieldShape {
            field: "validation budget",
        });
    }
    if cursor.length_prefixed("validation plan")?.len() != 14 * 16 {
        return Err(MachineDomainPortablePayloadError::FieldShape {
            field: "validation plan",
        });
    }
    if cursor.length_prefixed("lowering policy")?.is_empty() {
        return Err(MachineDomainPortablePayloadError::FieldShape {
            field: "lowering policy",
        });
    }
    let external_artifact_rows =
        cursor.skip_rows("external artifact rows", MAX_MACHINE_DOMAIN_ARTIFACTS)?;
    let crosswalk_rows = cursor.skip_rows("crosswalk rows", MAX_MACHINE_DOMAIN_CROSSWALKS)?;
    cursor.finish("manifest")?;

    let scenario_ir = core::str::from_utf8(scenario_bytes).map_err(|_| {
        MachineDomainPortablePayloadError::Utf8 {
            field: "canonical scenario IR",
        }
    })?;
    if hash_bytes(scenario_bytes) != scenario_hash {
        return Err(MachineDomainPortablePayloadError::ScenarioHash);
    }
    let decoded =
        parse_ir(scenario_ir).map_err(MachineDomainPortablePayloadError::ScenarioParse)?;
    if decoded.source_version() != SCENARIO_IR_VERSION
        || decoded.migration().is_some()
        || write_ir(decoded.scenario()) != scenario_ir
    {
        return Err(MachineDomainPortablePayloadError::ScenarioRoundTrip);
    }

    Ok(MachineDomainPortablePayloadViewV1 {
        graph,
        behavior,
        assurance,
        manifest_bytes: manifest,
        canonical_scenario_ir: scenario_ir,
        scenario_hash,
        manifest_hash: hash_bytes(manifest),
        payload_hash: hash_bytes(bytes),
        external_artifact_rows,
        crosswalk_rows,
    })
}

fn enforce_portable_limit(
    resource: &'static str,
    requested: usize,
    limit: usize,
) -> Result<(), MachineDomainPortablePayloadError> {
    if requested > limit {
        Err(MachineDomainPortablePayloadError::ResourceLimit {
            resource,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

struct PayloadCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> PayloadCursor<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(
        &mut self,
        field: &'static str,
        amount: usize,
    ) -> Result<&'a [u8], MachineDomainPortablePayloadError> {
        let end = self
            .offset
            .checked_add(amount)
            .ok_or(MachineDomainPortablePayloadError::Truncated { field })?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(MachineDomainPortablePayloadError::Truncated { field })?;
        self.offset = end;
        Ok(value)
    }

    fn expect_magic(
        &mut self,
        field: &'static str,
    ) -> Result<(), MachineDomainPortablePayloadError> {
        if self.take(field, PORTABLE_PAYLOAD_MAGIC.len())? == PORTABLE_PAYLOAD_MAGIC {
            Ok(())
        } else {
            Err(MachineDomainPortablePayloadError::Magic { field })
        }
    }

    fn read_u32(&mut self, field: &'static str) -> Result<u32, MachineDomainPortablePayloadError> {
        let raw: [u8; 4] = self
            .take(field, 4)?
            .try_into()
            .map_err(|_| MachineDomainPortablePayloadError::Truncated { field })?;
        Ok(u32::from_le_bytes(raw))
    }

    fn read_u64(&mut self, field: &'static str) -> Result<u64, MachineDomainPortablePayloadError> {
        let raw: [u8; 8] = self
            .take(field, 8)?
            .try_into()
            .map_err(|_| MachineDomainPortablePayloadError::Truncated { field })?;
        Ok(u64::from_le_bytes(raw))
    }

    fn expect_version(
        &mut self,
        field: &'static str,
        expected: u32,
    ) -> Result<(), MachineDomainPortablePayloadError> {
        let found = self.read_u32(field)?;
        if found == expected {
            Ok(())
        } else {
            Err(MachineDomainPortablePayloadError::Version {
                field,
                found,
                expected,
            })
        }
    }

    fn length_prefixed(
        &mut self,
        field: &'static str,
    ) -> Result<&'a [u8], MachineDomainPortablePayloadError> {
        let available = self.bytes.len();
        let encoded_length = self.read_u64(field)?;
        let length = usize::try_from(encoded_length).map_err(|_| {
            MachineDomainPortablePayloadError::ResourceLimit {
                resource: field,
                requested: usize::MAX,
                limit: available,
            }
        })?;
        self.take(field, length)
    }

    fn skip_rows(
        &mut self,
        field: &'static str,
        limit: usize,
    ) -> Result<usize, MachineDomainPortablePayloadError> {
        let encoded_count = self.read_u64(field)?;
        let count = usize::try_from(encoded_count).map_err(|_| {
            MachineDomainPortablePayloadError::ResourceLimit {
                resource: field,
                requested: usize::MAX,
                limit,
            }
        })?;
        enforce_portable_limit(field, count, limit)?;
        for _ in 0..count {
            if self.length_prefixed(field)?.is_empty() {
                return Err(MachineDomainPortablePayloadError::FieldShape { field });
            }
        }
        Ok(count)
    }

    fn finish(&self, field: &'static str) -> Result<(), MachineDomainPortablePayloadError> {
        if self.offset == self.bytes.len() {
            Ok(())
        } else {
            Err(MachineDomainPortablePayloadError::TrailingBytes { field })
        }
    }
}

/// Enumerate the complete canonical source vocabulary required for one exact
/// admitted Machine stack.
///
/// # Errors
/// Refuses a graph/behavior/assurance identity mismatch.
#[allow(clippy::too_many_lines)]
pub fn required_machine_lowering_sources_v1(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
) -> Result<Vec<MachineLoweringSourceV1>, MachineDomainLoweringRefusal> {
    required_machine_lowering_sources_impl(graph, behavior, assurance, None)
}

#[allow(clippy::too_many_lines)]
fn required_machine_lowering_sources_impl(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
    cx: Option<&Cx<'_>>,
) -> Result<Vec<MachineLoweringSourceV1>, MachineDomainLoweringRefusal> {
    validate_stack(graph, behavior, assurance)?;
    let mut sources = RequiredSourceCollector::new(cx);
    sources.insert(MachineLoweringSourceV1::Graph(graph.identity()))?;
    sources.insert(MachineLoweringSourceV1::Behavior(behavior.identity()))?;
    sources.insert(MachineLoweringSourceV1::Assurance(assurance.identity()))?;

    for clock in graph.clocks() {
        sources.insert(MachineLoweringSourceV1::Clock(clock.id.clone()))?;
    }
    for subsystem in graph.subsystems() {
        sources.insert(MachineLoweringSourceV1::Subsystem(subsystem.id.clone()))?;
        sources.extend(
            subsystem
                .bodies
                .iter()
                .cloned()
                .map(|id| MachineLoweringSourceV1::Element(id.into())),
        )?;
        sources.extend(
            subsystem
                .surface_patches
                .iter()
                .cloned()
                .map(|id| MachineLoweringSourceV1::Element(id.into())),
        )?;
        sources.extend(
            subsystem
                .contact_features
                .iter()
                .cloned()
                .map(|id| MachineLoweringSourceV1::Element(id.into())),
        )?;
        sources.extend(
            subsystem
                .state_slots
                .iter()
                .cloned()
                .map(|id| MachineLoweringSourceV1::Element(id.into())),
        )?;
    }
    sources.extend(
        graph
            .terminals()
            .iter()
            .map(|terminal| MachineLoweringSourceV1::Element(terminal.id.clone().into())),
    )?;
    sources.extend(
        graph
            .ports()
            .iter()
            .map(|port| MachineLoweringSourceV1::Element(port.id.clone().into())),
    )?;
    sources.extend(
        graph
            .relations()
            .iter()
            .map(|relation| MachineLoweringSourceV1::Relation(relation.id.clone())),
    )?;
    sources.extend(
        graph
            .interfaces()
            .iter()
            .map(|interface| MachineLoweringSourceV1::Interface(interface.id.clone())),
    )?;
    sources.extend(
        behavior
            .events()
            .iter()
            .map(|event| MachineLoweringSourceV1::Event(event.id.clone())),
    )?;
    sources.extend(
        behavior
            .tolerances()
            .iter()
            .map(|tolerance| MachineLoweringSourceV1::Tolerance(tolerance.id.clone())),
    )?;
    sources.extend(
        assurance
            .sensors()
            .iter()
            .map(|sensor| MachineLoweringSourceV1::Sensor(sensor.id.clone())),
    )?;
    sources.extend(
        assurance
            .experiments()
            .iter()
            .map(|experiment| MachineLoweringSourceV1::Experiment(experiment.id.clone())),
    )?;
    sources.extend(
        assurance
            .hazards()
            .iter()
            .map(|hazard| MachineLoweringSourceV1::Hazard(hazard.id.clone())),
    )?;
    sources.extend(
        assurance
            .faults()
            .iter()
            .map(|fault| MachineLoweringSourceV1::Fault(fault.id.clone())),
    )?;
    sources.extend(
        assurance
            .accounting_windows()
            .iter()
            .map(|window| MachineLoweringSourceV1::AccountingWindow(window.id.clone())),
    )?;
    sources.extend(
        assurance
            .fidelity()
            .rungs
            .iter()
            .map(|rung| MachineLoweringSourceV1::FidelityRung(rung.id.clone())),
    )?;

    sources.finish()
}

struct RequiredSourceCollector<'cx, 'arena> {
    sources: BTreeSet<MachineLoweringSourceV1>,
    cx: Option<&'cx Cx<'arena>>,
    visited: usize,
}

impl<'cx, 'arena> RequiredSourceCollector<'cx, 'arena> {
    fn new(cx: Option<&'cx Cx<'arena>>) -> Self {
        Self {
            sources: BTreeSet::new(),
            cx,
            visited: 0,
        }
    }

    fn insert(
        &mut self,
        source: MachineLoweringSourceV1,
    ) -> Result<(), MachineDomainLoweringRefusal> {
        if let Some(cx) = self.cx {
            poll_stride(cx, "required Machine-source collection", self.visited)?;
        }
        self.visited =
            self.visited
                .checked_add(1)
                .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
                    resource: "required source count",
                    requested: usize::MAX,
                    limit: MAX_MACHINE_DOMAIN_CROSSWALKS,
                })?;
        enforce_limit(
            "required source count",
            self.visited,
            MAX_MACHINE_DOMAIN_CROSSWALKS,
        )?;
        self.sources.insert(source);
        Ok(())
    }

    fn extend(
        &mut self,
        sources: impl IntoIterator<Item = MachineLoweringSourceV1>,
    ) -> Result<(), MachineDomainLoweringRefusal> {
        for source in sources {
            self.insert(source)?;
        }
        Ok(())
    }

    fn finish(self) -> Result<Vec<MachineLoweringSourceV1>, MachineDomainLoweringRefusal> {
        let mut sources = reserved_vec(self.sources.len(), "required Machine sources")?;
        sources.extend(self.sources);
        Ok(sources)
    }
}

/// Admit a caller-materialized one-way Machine-domain projection.
///
/// The supplied scenario is semantically validated and round-tripped through
/// current canonical scenario IR. Crosswalk rows are required for every exact
/// aggregate and durable local Machine identity. External artifacts and
/// mapping laws remain opaque, exact, versioned references; admission does not
/// authenticate or execute them.
///
/// # Errors
/// Returns a typed refusal for source-chain mismatch, resource excess,
/// incomplete/foreign/ambiguous crosswalks, invalid scenario semantics,
/// cancellation, canonical drift, or identity-publication failure.
pub fn admit_machine_domain_lowering_v1(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
    draft: MachineDomainLoweringDraftV1,
    validation_budget: ValidationBudget,
    cx: &Cx<'_>,
) -> Result<AdmittedMachineDomainLoweringV1, MachineDomainLoweringRefusal> {
    admit_machine_domain_lowering_inner_v1(graph, behavior, assurance, draft, validation_budget, cx)
}

/// Attempt one-way projection admission while retaining submitted counts and
/// a stable structured outcome code.
#[must_use]
pub fn admit_machine_domain_lowering_with_decision_v1(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
    draft: MachineDomainLoweringDraftV1,
    validation_budget: ValidationBudget,
    cx: &Cx<'_>,
) -> MachineDomainLoweringAdmissionDecisionV1 {
    let mut submitted = submitted_counts(&draft);
    let result = match submitted_target_count(&draft, cx) {
        Ok(targets) => {
            submitted.crosswalk_targets = Some(targets);
            admit_machine_domain_lowering_inner_v1(
                graph,
                behavior,
                assurance,
                draft,
                validation_budget,
                cx,
            )
        }
        Err(refusal) => Err(refusal),
    };
    MachineDomainLoweringAdmissionDecisionV1 { submitted, result }
}

#[allow(clippy::too_many_lines)]
fn admit_machine_domain_lowering_inner_v1(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
    mut draft: MachineDomainLoweringDraftV1,
    validation_budget: ValidationBudget,
    cx: &Cx<'_>,
) -> Result<AdmittedMachineDomainLoweringV1, MachineDomainLoweringRefusal> {
    checkpoint(cx, "initial")?;
    validate_stack(graph, behavior, assurance)?;
    preflight_draft(&draft, cx)?;
    preflight_scenario_records(&draft.scenario, cx)?;
    checkpoint(cx, "post-preflight")?;

    let validation_plan = draft
        .scenario
        .validation_plan(validation_budget)
        .map_err(MachineDomainLoweringRefusal::ScenarioValidation)?;
    enforce_scenario_serialization_bound(validation_plan)?;
    let findings = draft
        .scenario
        .validate_with_budget(validation_budget, cx)
        .map_err(|error| match error {
            ValidationError::Cancelled { phase, .. } => {
                MachineDomainLoweringRefusal::Cancelled { phase }
            }
            other => MachineDomainLoweringRefusal::ScenarioValidation(other),
        })?;
    if !findings.is_empty() {
        return Err(MachineDomainLoweringRefusal::ScenarioFindings(findings));
    }
    checkpoint(cx, "post-scenario-validation")?;

    checkpoint(cx, "pre-projection-sort")?;
    draft.external_artifacts.sort();
    if let Some(pair) = draft
        .external_artifacts
        .windows(2)
        .find(|pair| pair[0].id == pair[1].id)
    {
        return Err(MachineDomainLoweringRefusal::DuplicateArtifact {
            artifact: pair[0].id.clone(),
        });
    }
    draft
        .crosswalks
        .sort_by(|left, right| left.source.cmp(&right.source));
    checkpoint(cx, "post-projection-sort")?;
    if let Some(pair) = draft
        .crosswalks
        .windows(2)
        .find(|pair| pair[0].source == pair[1].source)
    {
        return Err(MachineDomainLoweringRefusal::DuplicateCrosswalk {
            source: pair[0].source.clone(),
        });
    }

    checkpoint(cx, "pre-required-source-collection")?;
    let required = required_machine_lowering_sources_impl(graph, behavior, assurance, Some(cx))?;
    checkpoint(cx, "post-required-source-collection")?;
    let mut actual_index = 0usize;
    let mut required_index = 0usize;
    while actual_index < draft.crosswalks.len() && required_index < required.len() {
        poll_stride(
            cx,
            "crosswalk source-set comparison",
            actual_index.saturating_add(required_index),
        )?;
        match draft.crosswalks[actual_index]
            .source
            .cmp(&required[required_index])
        {
            core::cmp::Ordering::Less => {
                return Err(MachineDomainLoweringRefusal::UnexpectedCrosswalk {
                    source: draft.crosswalks[actual_index].source.clone(),
                });
            }
            core::cmp::Ordering::Greater => {
                return Err(MachineDomainLoweringRefusal::MissingCrosswalk {
                    source: required[required_index].clone(),
                });
            }
            core::cmp::Ordering::Equal => {
                actual_index += 1;
                required_index += 1;
            }
        }
    }
    if let Some(entry) = draft.crosswalks.get(actual_index) {
        return Err(MachineDomainLoweringRefusal::UnexpectedCrosswalk {
            source: entry.source.clone(),
        });
    }
    if let Some(source) = required.get(required_index) {
        return Err(MachineDomainLoweringRefusal::MissingCrosswalk {
            source: source.clone(),
        });
    }

    let scenario_catalog = ScenarioLocatorCatalog::build(&draft.scenario, validation_plan, cx)?;
    validate_crosswalk_targets(
        &scenario_catalog,
        &draft.external_artifacts,
        &draft.crosswalks,
        cx,
    )?;
    checkpoint(cx, "post-crosswalk")?;

    let canonical_scenario_ir = write_ir(&draft.scenario);
    enforce_limit(
        "canonical scenario IR bytes",
        canonical_scenario_ir.len(),
        MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES,
    )?;
    let decoded =
        parse_ir(&canonical_scenario_ir).map_err(MachineDomainLoweringRefusal::ScenarioParse)?;
    if decoded.source_version() != SCENARIO_IR_VERSION {
        return Err(MachineDomainLoweringRefusal::ScenarioVersion {
            found: decoded.source_version(),
            expected: SCENARIO_IR_VERSION,
        });
    }
    if decoded.migration().is_some() {
        return Err(MachineDomainLoweringRefusal::UnexpectedScenarioMigration);
    }
    if decoded.scenario() != &draft.scenario
        || write_ir(decoded.scenario()) != canonical_scenario_ir
    {
        return Err(MachineDomainLoweringRefusal::ScenarioRoundTripMismatch);
    }

    let scenario_hash = hash_bytes(canonical_scenario_ir.as_bytes());
    let budget_row = validation_budget_row(validation_budget);
    let plan_row = validation_plan_row(validation_plan);
    let manifest_len = manifest_encoded_len_from_values(
        &budget_row,
        &plan_row,
        &draft.policy,
        &draft.external_artifacts,
        &draft.crosswalks,
    )
    .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
        resource: "canonical manifest bytes",
        requested: usize::MAX,
        limit: MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
    })?;
    enforce_limit(
        "canonical manifest bytes",
        manifest_len,
        MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
    )?;
    let portable_payload_len =
        portable_payload_encoded_len(manifest_len, canonical_scenario_ir.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "portable replay payload bytes",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_PORTABLE_PAYLOAD_BYTES,
            },
        )?;
    enforce_limit(
        "portable replay payload bytes",
        portable_payload_len,
        MAX_MACHINE_DOMAIN_PORTABLE_PAYLOAD_BYTES,
    )?;
    let policy_row = canonical_policy_row(&draft.policy)?;
    let artifact_rows = canonical_artifact_rows(&draft.external_artifacts, Some(cx))?;
    let crosswalk_rows = canonical_crosswalk_rows(&draft.crosswalks, Some(cx))?;
    let manifest_parts = ManifestParts {
        graph: graph.identity(),
        behavior: behavior.identity(),
        assurance: assurance.identity(),
        scenario_hash,
        validation_budget: &budget_row,
        validation_plan: &plan_row,
        policy: &policy_row,
        artifacts: &artifact_rows,
        crosswalks: &crosswalk_rows,
    };
    debug_assert_eq!(manifest_encoded_len(&manifest_parts), Some(manifest_len));
    let manifest = build_manifest(&manifest_parts, manifest_len)?;
    let portable_payload = build_portable_payload(
        &manifest,
        canonical_scenario_ir.as_bytes(),
        portable_payload_len,
    )?;
    checkpoint(cx, "pre-identity")?;

    let receipt = encode_lowering_identity(
        graph.identity(),
        behavior.identity(),
        assurance.identity(),
        scenario_hash,
        &budget_row,
        &plan_row,
        &policy_row,
        &artifact_rows,
        &crosswalk_rows,
        hash_bytes(&manifest),
        || cx.is_cancel_requested(),
    )
    .map_err(|error| match error {
        CanonicalError::Cancelled { .. } => {
            MachineDomainLoweringRefusal::Cancelled { phase: "identity" }
        }
        other => MachineDomainLoweringRefusal::Identity(other),
    })?;
    checkpoint(cx, "pre-publication")?;

    Ok(AdmittedMachineDomainLoweringV1 {
        graph: graph.identity(),
        behavior: behavior.identity(),
        assurance: assurance.identity(),
        scenario: draft.scenario,
        canonical_scenario_ir: canonical_scenario_ir.into_boxed_str(),
        scenario_hash,
        validation_budget,
        validation_plan,
        external_artifacts: draft.external_artifacts,
        crosswalks: draft.crosswalks,
        policy: draft.policy,
        canonical_manifest: manifest.into_boxed_slice(),
        portable_payload: portable_payload.into_boxed_slice(),
        receipt,
    })
}

fn validate_stack(
    graph: &AdmittedMachineGraph,
    behavior: &AdmittedMachineBehavior,
    assurance: &AdmittedMachineAssurance,
) -> Result<(), MachineDomainLoweringRefusal> {
    if behavior.base_graph() != graph.identity() {
        return Err(MachineDomainLoweringRefusal::BehaviorGraphMismatch);
    }
    if assurance.base_graph() != graph.identity() {
        return Err(MachineDomainLoweringRefusal::AssuranceGraphMismatch);
    }
    if assurance.base_behavior() != behavior.identity() {
        return Err(MachineDomainLoweringRefusal::AssuranceBehaviorMismatch);
    }
    Ok(())
}

fn preflight_draft(
    draft: &MachineDomainLoweringDraftV1,
    cx: &Cx<'_>,
) -> Result<(), MachineDomainLoweringRefusal> {
    enforce_limit(
        "external artifact count",
        draft.external_artifacts.len(),
        MAX_MACHINE_DOMAIN_ARTIFACTS,
    )?;
    enforce_limit(
        "crosswalk count",
        draft.crosswalks.len(),
        MAX_MACHINE_DOMAIN_CROSSWALKS,
    )?;
    let mut targets = 0usize;
    let mut bytes = draft.policy.namespace.len();
    for (artifact_index, artifact) in draft.external_artifacts.iter().enumerate() {
        poll_stride(cx, "external artifact preflight", artifact_index)?;
        bytes = bytes
            .checked_add(artifact.id.canonical_key.len())
            .and_then(|value| value.checked_add(artifact.namespace.len()))
            .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
                resource: "reference bytes",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
            })?;
        enforce_limit("reference bytes", bytes, MAX_MACHINE_DOMAIN_MANIFEST_BYTES)?;
    }
    let mut target_index = 0usize;
    for (entry_index, entry) in draft.crosswalks.iter().enumerate() {
        poll_stride(cx, "crosswalk row preflight", entry_index)?;
        targets = targets.checked_add(entry.targets.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "crosswalk target count",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_TARGETS,
            },
        )?;
        enforce_limit(
            "crosswalk target count",
            targets,
            MAX_MACHINE_DOMAIN_TARGETS,
        )?;
        bytes = bytes.checked_add(entry.law.namespace.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "reference bytes",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
            },
        )?;
        enforce_limit("reference bytes", bytes, MAX_MACHINE_DOMAIN_MANIFEST_BYTES)?;
        for target in &entry.targets {
            poll_stride(cx, "crosswalk target preflight", target_index)?;
            target_index += 1;
            bytes = bytes.checked_add(target_text_bytes(target)).ok_or(
                MachineDomainLoweringRefusal::ResourceLimit {
                    resource: "reference bytes",
                    requested: usize::MAX,
                    limit: MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
                },
            )?;
            enforce_limit("reference bytes", bytes, MAX_MACHINE_DOMAIN_MANIFEST_BYTES)?;
        }
    }
    Ok(())
}

fn target_text_bytes(target: &MachineDomainTargetV1) -> usize {
    match target {
        MachineDomainTargetV1::Scenario(locator) => match locator {
            MachineScenarioLocatorV1::Root | MachineScenarioLocatorV1::WorldFrame => 0,
            MachineScenarioLocatorV1::Frame { name, .. }
            | MachineScenarioLocatorV1::Region(name)
            | MachineScenarioLocatorV1::LoadCase(name)
            | MachineScenarioLocatorV1::Combination(name)
            | MachineScenarioLocatorV1::Ensemble(name) => name.len(),
            MachineScenarioLocatorV1::ContactPair { region_a, region_b } => {
                region_a.len().saturating_add(region_b.len())
            }
        },
        MachineDomainTargetV1::ExternalArtifact { selector, .. } => selector.len(),
    }
}

fn enforce_limit(
    resource: &'static str,
    requested: usize,
    limit: usize,
) -> Result<(), MachineDomainLoweringRefusal> {
    if requested > limit {
        Err(MachineDomainLoweringRefusal::ResourceLimit {
            resource,
            requested,
            limit,
        })
    } else {
        Ok(())
    }
}

fn checkpoint(cx: &Cx<'_>, phase: &'static str) -> Result<(), MachineDomainLoweringRefusal> {
    cx.checkpoint()
        .map_err(|_| MachineDomainLoweringRefusal::Cancelled { phase })
}

fn preflight_scenario_records(
    scenario: &Scenario,
    cx: &Cx<'_>,
) -> Result<(), MachineDomainLoweringRefusal> {
    let mut records = [
        scenario.frames.frames.len(),
        scenario.base_bcs.len(),
        scenario.cases.len(),
        scenario.combinations.len(),
        scenario.ensembles.len(),
        scenario.contacts.len(),
    ]
    .into_iter()
    .try_fold(0usize, usize::checked_add)
    .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
        resource: "scenario structural records",
        requested: usize::MAX,
        limit: MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
    })?;
    enforce_limit(
        "scenario structural records",
        records,
        MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
    )?;

    for (index, case) in scenario.cases.iter().enumerate() {
        poll_stride(cx, "scenario case-record preflight", index)?;
        records = records.checked_add(case.bcs.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "scenario structural records",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
            },
        )?;
        enforce_limit(
            "scenario structural records",
            records,
            MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
        )?;
    }
    for (index, combination) in scenario.combinations.iter().enumerate() {
        poll_stride(cx, "scenario combination-record preflight", index)?;
        records = records.checked_add(combination.terms.len()).ok_or(
            MachineDomainLoweringRefusal::ResourceLimit {
                resource: "scenario structural records",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
            },
        )?;
        enforce_limit(
            "scenario structural records",
            records,
            MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
        )?;
    }
    checkpoint(cx, "post-scenario-record-preflight")
}

fn enforce_scenario_serialization_bound(
    plan: ValidationPlan,
) -> Result<(), MachineDomainLoweringRefusal> {
    let records = [
        plan.frames,
        plan.base_bcs,
        plan.cases,
        plan.case_bcs,
        plan.combinations,
        plan.combination_terms,
        plan.ensembles,
        plan.contacts,
    ]
    .into_iter()
    .try_fold(0u128, |total, count| total.checked_add(count as u128))
    .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
        resource: "scenario serialization estimate bytes",
        requested: usize::MAX,
        limit: MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES,
    })?;
    // Scenario strings can expand by at most six bytes per input byte under
    // escaped IR text; an f64 needs fewer than 32 bytes including separator;
    // 2 KiB per fixed-shape record is deliberately conservative. This hard
    // pre-write bound keeps the non-cancellable legacy writer inside the
    // current parser's 16 MiB safety envelope.
    let estimate = 4_096u128
        .checked_add((plan.identity_bytes as u128).saturating_mul(6))
        .and_then(|value| value.checked_add((plan.signal_scalars as u128).saturating_mul(32)))
        .and_then(|value| value.checked_add(records.saturating_mul(2_048)))
        .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
            resource: "scenario serialization estimate bytes",
            requested: usize::MAX,
            limit: MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES,
        })?;
    if estimate > MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES as u128 {
        return Err(MachineDomainLoweringRefusal::ResourceLimit {
            resource: "scenario serialization estimate bytes",
            requested: usize::try_from(estimate).unwrap_or(usize::MAX),
            limit: MAX_MACHINE_DOMAIN_SCENARIO_IR_BYTES,
        });
    }
    Ok(())
}

struct ScenarioLocatorCatalog<'a> {
    frames: Vec<(u32, &'a str)>,
    regions: Vec<&'a str>,
    cases: Vec<&'a str>,
    combinations: Vec<&'a str>,
    ensembles: Vec<&'a str>,
    contacts: Vec<(&'a str, &'a str)>,
}

impl<'a> ScenarioLocatorCatalog<'a> {
    fn build(
        scenario: &'a Scenario,
        plan: ValidationPlan,
        cx: &Cx<'_>,
    ) -> Result<Self, MachineDomainLoweringRefusal> {
        let region_capacity = plan
            .base_bcs
            .checked_add(plan.case_bcs)
            .and_then(|value| value.checked_add(plan.contacts.saturating_mul(2)))
            .ok_or(MachineDomainLoweringRefusal::ResourceLimit {
                resource: "scenario locator region count",
                requested: usize::MAX,
                limit: MAX_MACHINE_DOMAIN_TARGETS,
            })?;
        let mut catalog = Self {
            frames: reserved_vec(plan.frames, "scenario frame locators")?,
            regions: reserved_vec(region_capacity, "scenario region locators")?,
            cases: reserved_vec(plan.cases, "scenario case locators")?,
            combinations: reserved_vec(plan.combinations, "scenario combination locators")?,
            ensembles: reserved_vec(plan.ensembles, "scenario ensemble locators")?,
            contacts: reserved_vec(plan.contacts, "scenario contact locators")?,
        };

        for (index, frame) in scenario.frames.frames.iter().enumerate() {
            poll_stride(cx, "scenario frame locator indexing", index)?;
            catalog.frames.push((frame.id.0, &frame.name));
        }
        for (index, boundary) in scenario.base_bcs.iter().enumerate() {
            poll_stride(cx, "scenario base-region locator indexing", index)?;
            catalog.regions.push(&boundary.region);
        }
        let mut boundary_index = scenario.base_bcs.len();
        for (case_index, case) in scenario.cases.iter().enumerate() {
            poll_stride(cx, "scenario case locator indexing", case_index)?;
            catalog.cases.push(&case.name);
            for boundary in &case.bcs {
                poll_stride(cx, "scenario case-region locator indexing", boundary_index)?;
                catalog.regions.push(&boundary.region);
                boundary_index += 1;
            }
        }
        for (index, combination) in scenario.combinations.iter().enumerate() {
            poll_stride(cx, "scenario combination locator indexing", index)?;
            catalog.combinations.push(&combination.name);
        }
        for (index, ensemble) in scenario.ensembles.iter().enumerate() {
            poll_stride(cx, "scenario ensemble locator indexing", index)?;
            catalog.ensembles.push(&ensemble.name);
        }
        for (index, contact) in scenario.contacts.iter().enumerate() {
            poll_stride(cx, "scenario contact locator indexing", index)?;
            catalog.regions.push(&contact.region_a);
            catalog.regions.push(&contact.region_b);
            catalog
                .contacts
                .push((&contact.region_a, &contact.region_b));
        }
        checkpoint(cx, "pre-scenario-locator-sort")?;
        catalog.frames.sort_unstable();
        catalog.regions.sort_unstable();
        catalog.cases.sort_unstable();
        catalog.combinations.sort_unstable();
        catalog.ensembles.sort_unstable();
        catalog.contacts.sort_unstable();
        checkpoint(cx, "post-scenario-locator-sort")?;
        Ok(catalog)
    }
}

fn reserved_vec<T>(
    requested: usize,
    resource: &'static str,
) -> Result<Vec<T>, MachineDomainLoweringRefusal> {
    let mut values = Vec::new();
    values.try_reserve_exact(requested).map_err(|_| {
        MachineDomainLoweringRefusal::AllocationRefused {
            resource,
            requested,
        }
    })?;
    Ok(values)
}

fn poll_stride(
    cx: &Cx<'_>,
    phase: &'static str,
    index: usize,
) -> Result<(), MachineDomainLoweringRefusal> {
    if index.is_multiple_of(256) {
        checkpoint(cx, phase)?;
    }
    Ok(())
}

fn validate_crosswalk_targets(
    scenario: &ScenarioLocatorCatalog<'_>,
    artifacts: &[MachineDomainArtifactRefV1],
    crosswalks: &[MachineDomainCrosswalkEntryV1],
    cx: &Cx<'_>,
) -> Result<(), MachineDomainLoweringRefusal> {
    let mut used_artifacts = reserved_vec(artifacts.len(), "used external-artifact flags")?;
    used_artifacts.resize(artifacts.len(), false);
    let mut target_index = 0usize;
    for entry in crosswalks {
        for target in &entry.targets {
            poll_stride(cx, "Machine crosswalk target validation", target_index)?;
            target_index += 1;
            match target {
                MachineDomainTargetV1::Scenario(locator) => {
                    validate_scenario_locator(scenario, locator)?;
                }
                MachineDomainTargetV1::ExternalArtifact { artifact, selector } => {
                    if selector.is_empty() {
                        return Err(MachineDomainLoweringRefusal::EmptyExternalSelector {
                            artifact: artifact.clone(),
                        });
                    }
                    if selector.len() > MAX_MACHINE_DOMAIN_SELECTOR_BYTES {
                        return Err(MachineDomainLoweringRefusal::ResourceLimit {
                            resource: "external target selector bytes",
                            requested: selector.len(),
                            limit: MAX_MACHINE_DOMAIN_SELECTOR_BYTES,
                        });
                    }
                    let artifact_index = artifacts
                        .binary_search_by(|candidate| candidate.id.cmp(artifact))
                        .map_err(|_| MachineDomainLoweringRefusal::UnknownExternalTarget {
                            artifact: artifact.clone(),
                        })?;
                    used_artifacts[artifact_index] = true;
                }
            }
        }
    }
    if let Some((index, _)) = used_artifacts.iter().enumerate().find(|(_, used)| !**used) {
        return Err(MachineDomainLoweringRefusal::OrphanExternalArtifact {
            artifact: artifacts[index].id.clone(),
        });
    }
    Ok(())
}

fn validate_scenario_locator(
    scenario: &ScenarioLocatorCatalog<'_>,
    locator: &MachineScenarioLocatorV1,
) -> Result<(), MachineDomainLoweringRefusal> {
    let text_valid =
        |value: &str| !value.is_empty() && value.len() <= MAX_MACHINE_DOMAIN_SELECTOR_BYTES;
    let count = match locator {
        MachineScenarioLocatorV1::Root | MachineScenarioLocatorV1::WorldFrame => return Ok(()),
        MachineScenarioLocatorV1::Frame { id, name } => {
            if !text_valid(name) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            sorted_count(&scenario.frames, &(*id, name.as_ref()))
        }
        MachineScenarioLocatorV1::Region(region) => {
            if !text_valid(region) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            usize::from(scenario.regions.binary_search(&region.as_ref()).is_ok())
        }
        MachineScenarioLocatorV1::LoadCase(name) => {
            if !text_valid(name) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            sorted_count(&scenario.cases, &name.as_ref())
        }
        MachineScenarioLocatorV1::Combination(name) => {
            if !text_valid(name) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            sorted_count(&scenario.combinations, &name.as_ref())
        }
        MachineScenarioLocatorV1::Ensemble(name) => {
            if !text_valid(name) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            sorted_count(&scenario.ensembles, &name.as_ref())
        }
        MachineScenarioLocatorV1::ContactPair { region_a, region_b } => {
            if !text_valid(region_a) || !text_valid(region_b) {
                return Err(MachineDomainLoweringRefusal::ScenarioTarget {
                    locator: locator.clone(),
                    problem: MachineScenarioLocatorProblemV1::InvalidText,
                });
            }
            sorted_count(&scenario.contacts, &(region_a.as_ref(), region_b.as_ref()))
        }
    };
    match count {
        1 => Ok(()),
        0 => Err(MachineDomainLoweringRefusal::ScenarioTarget {
            locator: locator.clone(),
            problem: MachineScenarioLocatorProblemV1::Missing,
        }),
        _ => Err(MachineDomainLoweringRefusal::ScenarioTarget {
            locator: locator.clone(),
            problem: MachineScenarioLocatorProblemV1::Ambiguous,
        }),
    }
}

fn sorted_count<T: Ord>(values: &[T], needle: &T) -> usize {
    let first = values.partition_point(|value| value < needle);
    let after = values.partition_point(|value| value <= needle);
    after - first
}

fn canonical_policy_row(
    policy: &MachineLoweringPolicyRefV1,
) -> Result<Vec<u8>, MachineDomainLoweringRefusal> {
    let length = encoded_len_or_limit(policy.canonical_row_len(), "lowering policy row bytes")?;
    let mut row = reserved_bytes(length, "lowering policy row bytes")?;
    policy.append_canonical_row(&mut row);
    debug_assert_eq!(row.len(), length);
    Ok(row)
}

fn canonical_artifact_rows(
    artifacts: &[MachineDomainArtifactRefV1],
    cx: Option<&Cx<'_>>,
) -> Result<Vec<Vec<u8>>, MachineDomainLoweringRefusal> {
    let mut rows = reserved_vec(artifacts.len(), "canonical external-artifact rows")?;
    for (index, artifact) in artifacts.iter().enumerate() {
        if let Some(cx) = cx {
            poll_stride(cx, "external-artifact row encoding", index)?;
        }
        let length =
            encoded_len_or_limit(artifact.canonical_row_len(), "external-artifact row bytes")?;
        let mut row = reserved_bytes(length, "external-artifact row bytes")?;
        artifact.append_canonical_row(&mut row);
        debug_assert_eq!(row.len(), length);
        rows.push(row);
    }
    Ok(rows)
}

fn canonical_crosswalk_rows(
    crosswalks: &[MachineDomainCrosswalkEntryV1],
    cx: Option<&Cx<'_>>,
) -> Result<Vec<Vec<u8>>, MachineDomainLoweringRefusal> {
    let mut rows = reserved_vec(crosswalks.len(), "canonical crosswalk rows")?;
    for (index, crosswalk) in crosswalks.iter().enumerate() {
        if let Some(cx) = cx {
            poll_stride(cx, "crosswalk row encoding", index)?;
        }
        let length = encoded_len_or_limit(crosswalk.canonical_row_len(), "crosswalk row bytes")?;
        let mut row = reserved_bytes(length, "crosswalk row bytes")?;
        crosswalk.append_canonical_row(&mut row);
        debug_assert_eq!(row.len(), length);
        rows.push(row);
    }
    Ok(rows)
}

fn encoded_len_or_limit(
    length: Option<usize>,
    resource: &'static str,
) -> Result<usize, MachineDomainLoweringRefusal> {
    length.ok_or(MachineDomainLoweringRefusal::ResourceLimit {
        resource,
        requested: usize::MAX,
        limit: MAX_MACHINE_DOMAIN_MANIFEST_BYTES,
    })
}

fn reserved_bytes(
    requested: usize,
    resource: &'static str,
) -> Result<Vec<u8>, MachineDomainLoweringRefusal> {
    let mut bytes = Vec::new();
    bytes.try_reserve_exact(requested).map_err(|_| {
        MachineDomainLoweringRefusal::AllocationRefused {
            resource,
            requested,
        }
    })?;
    Ok(bytes)
}

fn validation_budget_row(budget: ValidationBudget) -> Vec<u8> {
    let mut row = Vec::with_capacity(14 * 16);
    for value in [
        budget.max_frames,
        budget.max_base_bcs,
        budget.max_cases,
        budget.max_case_bcs,
        budget.max_combinations,
        budget.max_combination_terms,
        budget.max_ensembles,
        budget.max_contacts,
        budget.max_signal_scalars,
        budget.max_flux_checkpoints,
        budget.max_identity_bytes,
        budget.max_identity_component_bytes,
        budget.max_findings,
    ] {
        row.extend_from_slice(&(value as u128).to_le_bytes());
    }
    row.extend_from_slice(&budget.max_work.to_le_bytes());
    row
}

fn validation_plan_row(plan: ValidationPlan) -> Vec<u8> {
    let mut row = Vec::with_capacity(14 * 16);
    for value in [
        plan.frames,
        plan.base_bcs,
        plan.cases,
        plan.case_bcs,
        plan.combinations,
        plan.combination_terms,
        plan.ensembles,
        plan.contacts,
        plan.signal_scalars,
        plan.flux_checkpoints,
        plan.identity_bytes,
        plan.identity_component_bytes,
        plan.finding_capacity,
    ] {
        row.extend_from_slice(&(value as u128).to_le_bytes());
    }
    row.extend_from_slice(&plan.planned_work.to_le_bytes());
    row
}

struct ManifestParts<'a> {
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    assurance: MachineAssuranceIdV1,
    scenario_hash: ContentHash,
    validation_budget: &'a [u8],
    validation_plan: &'a [u8],
    policy: &'a [u8],
    artifacts: &'a [Vec<u8>],
    crosswalks: &'a [Vec<u8>],
}

fn manifest_encoded_len(parts: &ManifestParts<'_>) -> Option<usize> {
    let mut len = PORTABLE_PAYLOAD_MAGIC
        .len()
        .checked_add(core::mem::size_of::<u32>() * 3)?
        .checked_add(32 * 4)?;
    for bytes in [
        FS_IR_VERSION.as_bytes(),
        fs_scenario::VERSION.as_bytes(),
        parts.validation_budget,
        parts.validation_plan,
        parts.policy,
    ] {
        len = len.checked_add(8)?.checked_add(bytes.len())?;
    }
    for rows in [parts.artifacts, parts.crosswalks] {
        len = len.checked_add(8)?;
        for row in rows {
            len = len.checked_add(8)?.checked_add(row.len())?;
        }
    }
    Some(len)
}

fn manifest_encoded_len_from_values(
    validation_budget: &[u8],
    validation_plan: &[u8],
    policy: &MachineLoweringPolicyRefV1,
    artifacts: &[MachineDomainArtifactRefV1],
    crosswalks: &[MachineDomainCrosswalkEntryV1],
) -> Option<usize> {
    let mut len = PORTABLE_PAYLOAD_MAGIC
        .len()
        .checked_add(core::mem::size_of::<u32>() * 3)?
        .checked_add(32 * 4)?;
    for length in [
        FS_IR_VERSION.len(),
        fs_scenario::VERSION.len(),
        validation_budget.len(),
        validation_plan.len(),
        policy.canonical_row_len()?,
    ] {
        len = len.checked_add(8)?.checked_add(length)?;
    }
    len = len.checked_add(8)?;
    for artifact in artifacts {
        len = len
            .checked_add(8)?
            .checked_add(artifact.canonical_row_len()?)?;
    }
    len = len.checked_add(8)?;
    for crosswalk in crosswalks {
        len = len
            .checked_add(8)?
            .checked_add(crosswalk.canonical_row_len()?)?;
    }
    Some(len)
}

fn build_manifest(
    parts: &ManifestParts<'_>,
    capacity: usize,
) -> Result<Vec<u8>, MachineDomainLoweringRefusal> {
    let mut out = reserved_bytes(capacity, "canonical manifest bytes")?;
    out.extend_from_slice(PORTABLE_PAYLOAD_MAGIC);
    out.extend_from_slice(&MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1.to_le_bytes());
    out.extend_from_slice(&IR_VERSION.to_le_bytes());
    append_bytes(&mut out, FS_IR_VERSION.as_bytes());
    append_bytes(&mut out, fs_scenario::VERSION.as_bytes());
    out.extend_from_slice(parts.graph.as_bytes());
    out.extend_from_slice(parts.behavior.as_bytes());
    out.extend_from_slice(parts.assurance.as_bytes());
    out.extend_from_slice(&SCENARIO_IR_VERSION.to_le_bytes());
    out.extend_from_slice(parts.scenario_hash.as_bytes());
    append_bytes(&mut out, parts.validation_budget);
    append_bytes(&mut out, parts.validation_plan);
    append_bytes(&mut out, parts.policy);
    append_rows(&mut out, parts.artifacts);
    append_rows(&mut out, parts.crosswalks);
    debug_assert_eq!(out.len(), capacity);
    Ok(out)
}

fn portable_payload_encoded_len(manifest: usize, scenario_ir: usize) -> Option<usize> {
    PORTABLE_PAYLOAD_MAGIC
        .len()
        .checked_add(16)?
        .checked_add(manifest)?
        .checked_add(scenario_ir)
}

fn build_portable_payload(
    manifest: &[u8],
    scenario_ir: &[u8],
    capacity: usize,
) -> Result<Vec<u8>, MachineDomainLoweringRefusal> {
    let mut out = reserved_bytes(capacity, "portable replay payload bytes")?;
    out.extend_from_slice(PORTABLE_PAYLOAD_MAGIC);
    append_bytes(&mut out, manifest);
    append_bytes(&mut out, scenario_ir);
    debug_assert_eq!(out.len(), capacity);
    Ok(out)
}

fn append_length(out: &mut Vec<u8>, length: usize) {
    out.extend_from_slice(&(length as u64).to_le_bytes());
}

fn append_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    append_length(out, bytes.len());
    out.extend_from_slice(bytes);
}

fn append_rows(out: &mut Vec<u8>, rows: &[Vec<u8>]) {
    append_length(out, rows.len());
    for row in rows {
        append_bytes(out, row);
    }
}

#[allow(clippy::too_many_arguments)]
fn encode_lowering_identity<C>(
    graph: MachineGraphIdV1,
    behavior: MachineBehaviorIdV1,
    assurance: MachineAssuranceIdV1,
    scenario_hash: ContentHash,
    validation_budget: &[u8],
    validation_plan: &[u8],
    policy: &[u8],
    artifacts: &[Vec<u8>],
    crosswalks: &[Vec<u8>],
    manifest_hash: ContentHash,
    cancellation: C,
) -> Result<IdentityReceipt<MachineDomainLoweringIdV1>, CanonicalError>
where
    C: fs_blake3::identity::CancellationProbe,
{
    CanonicalEncoder::<MachineDomainLoweringIdV1, _>::new(
        MACHINE_DOMAIN_IDENTITY_LIMITS,
        cancellation,
    )?
    .u64(
        Field::new(0, "lowering-schema-version"),
        u64::from(MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph.as_bytes())?
    .bytes(Field::new(3, "machine-behavior"), behavior.as_bytes())?
    .bytes(Field::new(4, "machine-assurance"), assurance.as_bytes())?
    .u64(
        Field::new(5, "scenario-ir-version"),
        u64::from(SCENARIO_IR_VERSION),
    )?
    .bytes(Field::new(6, "scenario-ir-hash"), scenario_hash.as_bytes())?
    .bytes(
        Field::new(7, "scenario-validation-budget"),
        validation_budget,
    )?
    .bytes(Field::new(8, "scenario-validation-plan"), validation_plan)?
    .bytes(Field::new(9, "lowering-policy"), policy)?
    .ordered_bytes(
        Field::new(10, "external-artifacts"),
        artifacts.len() as u64,
        artifacts.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(11, "crosswalks"),
        crosswalks.len() as u64,
        crosswalks.iter().map(Vec::as_slice),
    )?
    .bytes(
        Field::new(12, "canonical-manifest-hash"),
        manifest_hash.as_bytes(),
    )?
    .finish()
}
