//! Graph-bound assembly chronology and family-specific joint topology.
//!
//! Version two deliberately separates chronological body availability from
//! physical joint occurrences. Durable contact features may be reused through
//! distinct feature-use identities, unordered hyperedge members canonicalize
//! independently of chronology, and directional semantics exist only as
//! family-specific physical roles. Admission proves structural graph ownership
//! and a deterministic availability history; it does not prove containment,
//! path feasibility, execution, or joint physics.

use core::fmt;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, ContentId, Field,
    FieldSpec, IdentityReceipt, NeverCancel, ProblemSemanticId, StrongIdentity, WireType,
};

use crate::IR_VERSION;

use super::super::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphIdV1, MachineIdError, SubsystemId,
};
use super::ManufacturingArtifactRefV1;

/// Identity and admission schema version for the corrected assembly model.
pub const MACHINE_ASSEMBLY_SCHEMA_VERSION_V2: u32 = 2;
/// Version of the authenticated availability-transition commitment carried by V2.
pub const MACHINE_ASSEMBLY_AVAILABILITY_COMMITMENT_VERSION_V2: u32 = 1;
/// Maximum chronological steps retained by one assembly receipt.
pub const MAX_MACHINE_ASSEMBLY_STEPS_V2: usize = 4_096;
/// Maximum physical joint occurrences retained by one assembly receipt.
pub const MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2: usize = 4_096;
/// Maximum bodies available before chronological step zero.
pub const MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2: usize = 4_096;
/// Maximum bodies atomically introduced by one chronological step.
pub const MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2: usize = 64;
/// Maximum physical occurrences atomically admitted by one chronological step.
pub const MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2: usize = 64;
/// Maximum typed physical participants in one joint occurrence.
pub const MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2: usize = 64;

/// Explicit canonical resource envelope for V2 assembly identity publication.
///
/// The 160 MiB field ceiling covers the independently derived 138,235,912-byte
/// field containing 4,096 true maximum-width `ExecutionClaimed` preloaded-bolt
/// rows. The 256 MiB aggregate ceiling covers the independently derived
/// 226,037,784-byte three-collection payload plus fixed frame overhead. The
/// availability evidence hashes the initial set once and then only each
/// bounded transition, so no complete growing set is serialized per step.
pub const MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2: CanonicalLimits =
    CanonicalLimits::new(256 * 1_024 * 1_024, 160 * 1_024 * 1_024, 6, 4_096, 4_096);

const AVAILABILITY_INITIAL_ROOT_DOMAIN_V2: &str = "org.frankensim.fs-ir.machine.manufacturing-assembly.v2/availability-transition-chain.v1/initial";
const AVAILABILITY_STEP_ROOT_DOMAIN_V2: &str =
    "org.frankensim.fs-ir.machine.manufacturing-assembly.v2/availability-transition-chain.v1/step";

macro_rules! assembly_key_id_v2 {
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

            /// Exact canonical key retained in aggregate identity.
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

assembly_key_id_v2!(
    /// Stable identity of one chronological availability step.
    AssemblyStepIdV2,
    "assembly-step-id"
);
assembly_key_id_v2!(
    /// Stable identity of one physical joint occurrence.
    JointOccurrenceIdV2,
    "assembly-joint-occurrence-id"
);
assembly_key_id_v2!(
    /// Assembly-declaration-local identity of one occurrence-local durable-feature use.
    JointFeatureUseIdV2,
    "assembly-joint-feature-use-id"
);

macro_rules! assembly_artifact_role_v2 {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(ManufacturingArtifactRefV1);

        impl $name {
            /// Assign an admitted artifact coordinate to this nominal role.
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

assembly_artifact_role_v2!(
    /// Exact joining procedure/specification coordinate.
    AssemblyProcedureRefV2
);
assembly_artifact_role_v2!(
    /// Exact nominal insertion/approach-path coordinate.
    AssemblyPathRefV2
);
assembly_artifact_role_v2!(
    /// Exact nominal evidence coordinate accompanying an execution claim.
    AssemblyExecutionEvidenceRefV2
);

/// Closed truthful lifecycle for one physical joint occurrence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssemblyLifecycleV2 {
    /// A declared procedure and path that have not been claimed as executed.
    Planned {
        /// Exact nominal procedure coordinate.
        procedure: AssemblyProcedureRefV2,
        /// Exact nominal path coordinate.
        path: AssemblyPathRefV2,
    },
    /// A caller claim of execution with nominal evidence, not verification.
    ExecutionClaimed {
        /// Exact nominal procedure coordinate.
        procedure: AssemblyProcedureRefV2,
        /// Exact nominal path coordinate.
        path: AssemblyPathRefV2,
        /// Exact nominal evidence coordinate supporting the caller claim.
        evidence: AssemblyExecutionEvidenceRefV2,
    },
}

impl AssemblyLifecycleV2 {
    /// Stable lifecycle discriminant.
    #[must_use]
    pub const fn tag(&self) -> u8 {
        match self {
            Self::Planned { .. } => 1,
            Self::ExecutionClaimed { .. } => 2,
        }
    }

    /// Exact nominal procedure coordinate in either lifecycle state.
    #[must_use]
    pub const fn procedure(&self) -> &AssemblyProcedureRefV2 {
        match self {
            Self::Planned { procedure, .. } | Self::ExecutionClaimed { procedure, .. } => procedure,
        }
    }

    /// Exact nominal path coordinate in either lifecycle state.
    #[must_use]
    pub const fn path(&self) -> &AssemblyPathRefV2 {
        match self {
            Self::Planned { path, .. } | Self::ExecutionClaimed { path, .. } => path,
        }
    }

    /// Nominal execution evidence only when execution is explicitly claimed.
    #[must_use]
    pub const fn execution_evidence(&self) -> Option<&AssemblyExecutionEvidenceRefV2> {
        match self {
            Self::Planned { .. } => None,
            Self::ExecutionClaimed { evidence, .. } => Some(evidence),
        }
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(384);
        row.push(self.tag());
        append_artifact_v2(&mut row, self.procedure().artifact());
        append_artifact_v2(&mut row, self.path().artifact());
        if let Some(evidence) = self.execution_evidence() {
            append_artifact_v2(&mut row, evidence.artifact());
        }
        row
    }
}

/// Explicit source unit for a preloaded-bolt force target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AssemblyPreloadUnitV2 {
    /// Newtons.
    Newton = 1,
    /// Kilonewtons.
    Kilonewton = 2,
}

impl AssemblyPreloadUnitV2 {
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

/// Refusal from constructing a strictly positive finite preload target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyPreloadErrorV2 {
    /// NaN or infinity was supplied.
    NonFinite,
    /// Force target was zero or negative.
    NonPositive,
    /// Unit normalization overflowed binary64.
    SiNonFinite,
}

impl AssemblyPreloadErrorV2 {
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

impl fmt::Display for AssemblyPreloadErrorV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::NonFinite => "assembly preload target must be finite",
            Self::NonPositive => "assembly preload target must be positive",
            Self::SiNonFinite => "assembly preload SI normalization must remain finite",
        })
    }
}

impl std::error::Error for AssemblyPreloadErrorV2 {}

/// Strictly positive preload target retaining source and coherent-SI bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AssemblyPreloadV2 {
    submitted_bits: u64,
    unit: AssemblyPreloadUnitV2,
    newtons_bits: u64,
}

impl AssemblyPreloadV2 {
    /// Validate and normalize one positive preload target.
    ///
    /// # Errors
    /// Refuses non-finite, non-positive, or SI-overflowing input.
    pub fn try_new(
        value: f64,
        unit: AssemblyPreloadUnitV2,
    ) -> Result<Self, AssemblyPreloadErrorV2> {
        if !value.is_finite() {
            return Err(AssemblyPreloadErrorV2::NonFinite);
        }
        if value <= 0.0 {
            return Err(AssemblyPreloadErrorV2::NonPositive);
        }
        let newtons = value * unit.newtons_per_unit();
        if !newtons.is_finite() {
            return Err(AssemblyPreloadErrorV2::SiNonFinite);
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
    pub const fn unit(self) -> AssemblyPreloadUnitV2 {
        self.unit
    }

    /// Coherent-SI binary64 force in newtons.
    #[must_use]
    pub fn newtons(self) -> f64 {
        f64::from_bits(self.newtons_bits)
    }

    /// Exact submitted binary64 bits.
    #[must_use]
    pub const fn submitted_bits(self) -> u64 {
        self.submitted_bits
    }

    /// Exact coherent-SI binary64 bits.
    #[must_use]
    pub const fn newtons_bits(self) -> u64 {
        self.newtons_bits
    }

    fn append_canonical(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.submitted_bits.to_le_bytes());
        out.push(self.unit.tag());
        out.extend_from_slice(&self.newtons_bits.to_le_bytes());
    }
}

/// Explicit physical-feature reuse policy carried by each occurrence-local use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum PhysicalFeatureUsePolicyV2 {
    /// The durable physical feature may support other separately identified uses.
    Reusable = 1,
    /// No other use in this assembly declaration may select the feature.
    ExclusiveWithinAssembly = 2,
}

impl PhysicalFeatureUsePolicyV2 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// Caller-declared body/contact-feature selector.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AssemblyFeatureSelectorV2 {
    declared_body: BodyId,
    contact_feature: ContactFeatureId,
}

impl AssemblyFeatureSelectorV2 {
    /// Construct one authority-free selector.
    #[must_use]
    pub fn new(declared_body: BodyId, contact_feature: ContactFeatureId) -> Self {
        Self {
            declared_body,
            contact_feature,
        }
    }

    /// Caller-declared body; physical containment is not proved by V2.
    #[must_use]
    pub const fn declared_body(&self) -> &BodyId {
        &self.declared_body
    }

    /// Durable physical contact feature.
    #[must_use]
    pub const fn contact_feature(&self) -> &ContactFeatureId {
        &self.contact_feature
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(352);
        append_body_v2(&mut row, &self.declared_body);
        append_feature_v2(&mut row, &self.contact_feature);
        row
    }
}

/// One occurrence-local use of a reusable or explicitly exclusive feature.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct JointFeatureUseV2 {
    id: JointFeatureUseIdV2,
    selector: AssemblyFeatureSelectorV2,
    policy: PhysicalFeatureUsePolicyV2,
}

impl JointFeatureUseV2 {
    /// Construct one authority-free physical-feature use.
    #[must_use]
    pub fn new(
        id: JointFeatureUseIdV2,
        selector: AssemblyFeatureSelectorV2,
        policy: PhysicalFeatureUsePolicyV2,
    ) -> Self {
        Self {
            id,
            selector,
            policy,
        }
    }

    /// Stable identity, unique within one assembly declaration, for this use.
    #[must_use]
    pub const fn id(&self) -> &JointFeatureUseIdV2 {
        &self.id
    }

    /// Physical body/feature selector.
    #[must_use]
    pub const fn selector(&self) -> &AssemblyFeatureSelectorV2 {
        &self.selector
    }

    /// Explicit physical-feature reuse policy.
    #[must_use]
    pub const fn policy(&self) -> PhysicalFeatureUsePolicyV2 {
        self.policy
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(512);
        append_bytes_v2(&mut row, self.id.canonical_key().as_bytes());
        row.push(self.policy.tag());
        append_bytes_v2(&mut row, &self.selector.canonical_row());
        row
    }
}

/// Typed participant in the physically ordered fastener stack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum BoltStackRoleV2 {
    /// Bolt or screw body, exactly one per occurrence.
    Bolt = 1,
    /// Nut body, at most one per occurrence.
    Nut = 2,
    /// Washer body.
    Washer = 3,
    /// Spacer or sleeve body.
    Spacer = 4,
    /// Locking element body.
    LockingElement = 5,
}

impl BoltStackRoleV2 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }
}

/// One position in a bolt-head-to-thread-end fastener stack.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoltStackParticipantV2 {
    position: u16,
    role: BoltStackRoleV2,
    feature_use: JointFeatureUseV2,
}

impl BoltStackParticipantV2 {
    /// Construct one typed fastener-stack participant.
    #[must_use]
    pub fn new(position: u16, role: BoltStackRoleV2, feature_use: JointFeatureUseV2) -> Self {
        Self {
            position,
            role,
            feature_use,
        }
    }

    /// Zero-based physical position from bolt head toward threaded end.
    #[must_use]
    pub const fn position(&self) -> u16 {
        self.position
    }

    /// Closed fastener-stack role.
    #[must_use]
    pub const fn role(&self) -> BoltStackRoleV2 {
        self.role
    }

    /// Occurrence-local physical-feature use.
    #[must_use]
    pub const fn feature_use(&self) -> &JointFeatureUseV2 {
        &self.feature_use
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(520);
        row.extend_from_slice(&self.position.to_le_bytes());
        row.push(self.role.tag());
        append_bytes_v2(&mut row, &self.feature_use.canonical_row());
        row
    }
}

/// Closed physical joint-family vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum AssemblyJointFamilyV2 {
    /// Preloaded bolt hyperedge.
    PreloadedBolt = 1,
    /// Weld member hyperedge.
    Weld = 2,
    /// Adhesive adherend hyperedge.
    AdhesiveBond = 3,
    /// Directed shaft/hub/key topology.
    Key = 4,
    /// Directed external/internal spline topology.
    Spline = 5,
    /// Directed external/internal interference-fit topology.
    InterferenceFit = 6,
}

impl AssemblyJointFamilyV2 {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
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

/// Family-specific physical topology for one durable joint occurrence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JointTopologyV2 {
    /// A bolt hyperedge with unordered clamped members and an ordered stack.
    PreloadedBolt {
        /// Unordered clamped-member participants; at least two are required.
        clamped_members: Vec<JointFeatureUseV2>,
        /// Bolt-head-to-thread-end typed participants; exactly one bolt and at
        /// most one nut are required.
        fastener_stack: Vec<BoltStackParticipantV2>,
        /// Positive finite force declaration retained only for this family.
        preload: AssemblyPreloadV2,
    },
    /// An unordered welded-member hyperedge with at least two members.
    Weld {
        /// Unordered welded physical members.
        members: Vec<JointFeatureUseV2>,
    },
    /// An unordered adhesive-adherend hyperedge with at least two adherends.
    AdhesiveBond {
        /// Unordered adherend physical members.
        adherends: Vec<JointFeatureUseV2>,
    },
    /// A physically directed three-body keyed topology.
    Key {
        /// Shaft-side participant.
        shaft: JointFeatureUseV2,
        /// Hub-side participant.
        hub: JointFeatureUseV2,
        /// Distinct key-body participant.
        key: JointFeatureUseV2,
    },
    /// A physically directed spline topology.
    Spline {
        /// External spline participant.
        external: JointFeatureUseV2,
        /// Internal spline participant.
        internal: JointFeatureUseV2,
    },
    /// A physically directed interference-fit topology.
    InterferenceFit {
        /// External fit participant.
        external: JointFeatureUseV2,
        /// Internal fit participant.
        internal: JointFeatureUseV2,
    },
}

impl JointTopologyV2 {
    /// Closed family represented by this exact payload.
    #[must_use]
    pub const fn family(&self) -> AssemblyJointFamilyV2 {
        match self {
            Self::PreloadedBolt { .. } => AssemblyJointFamilyV2::PreloadedBolt,
            Self::Weld { .. } => AssemblyJointFamilyV2::Weld,
            Self::AdhesiveBond { .. } => AssemblyJointFamilyV2::AdhesiveBond,
            Self::Key { .. } => AssemblyJointFamilyV2::Key,
            Self::Spline { .. } => AssemblyJointFamilyV2::Spline,
            Self::InterferenceFit { .. } => AssemblyJointFamilyV2::InterferenceFit,
        }
    }

    fn canonicalize_unordered_members(&mut self) {
        match self {
            Self::PreloadedBolt {
                clamped_members,
                fastener_stack,
                ..
            } => {
                clamped_members.sort_by_cached_key(JointFeatureUseV2::canonical_row);
                fastener_stack.sort_by(|left, right| {
                    left.position
                        .cmp(&right.position)
                        .then_with(|| left.canonical_row().cmp(&right.canonical_row()))
                });
            }
            Self::Weld { members } => {
                members.sort_by_cached_key(JointFeatureUseV2::canonical_row);
            }
            Self::AdhesiveBond { adherends } => {
                adherends.sort_by_cached_key(JointFeatureUseV2::canonical_row);
            }
            Self::Key { .. } | Self::Spline { .. } | Self::InterferenceFit { .. } => {}
        }
    }

    /// Typed participants in canonical family order.
    #[must_use]
    pub fn participants(&self) -> Vec<(AssemblyParticipantRoleV2, &JointFeatureUseV2)> {
        match self {
            Self::PreloadedBolt {
                clamped_members,
                fastener_stack,
                ..
            } => clamped_members
                .iter()
                .enumerate()
                .map(|(index, feature_use)| {
                    (
                        AssemblyParticipantRoleV2::BoltClampedMember {
                            canonical_index: index,
                        },
                        feature_use,
                    )
                })
                .chain(fastener_stack.iter().map(|participant| {
                    (
                        AssemblyParticipantRoleV2::BoltStack {
                            position: participant.position,
                            role: participant.role,
                        },
                        &participant.feature_use,
                    )
                }))
                .collect(),
            Self::Weld { members } => members
                .iter()
                .enumerate()
                .map(|(index, feature_use)| {
                    (
                        AssemblyParticipantRoleV2::WeldMember {
                            canonical_index: index,
                        },
                        feature_use,
                    )
                })
                .collect(),
            Self::AdhesiveBond { adherends } => adherends
                .iter()
                .enumerate()
                .map(|(index, feature_use)| {
                    (
                        AssemblyParticipantRoleV2::AdhesiveAdherend {
                            canonical_index: index,
                        },
                        feature_use,
                    )
                })
                .collect(),
            Self::Key { shaft, hub, key } => vec![
                (AssemblyParticipantRoleV2::KeyShaft, shaft),
                (AssemblyParticipantRoleV2::KeyHub, hub),
                (AssemblyParticipantRoleV2::KeyBody, key),
            ],
            Self::Spline { external, internal } => vec![
                (AssemblyParticipantRoleV2::SplineExternal, external),
                (AssemblyParticipantRoleV2::SplineInternal, internal),
            ],
            Self::InterferenceFit { external, internal } => vec![
                (AssemblyParticipantRoleV2::InterferenceExternal, external),
                (AssemblyParticipantRoleV2::InterferenceInternal, internal),
            ],
        }
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(1_024);
        row.push(self.family().tag());
        match self {
            Self::PreloadedBolt {
                clamped_members,
                fastener_stack,
                preload,
            } => {
                append_rows_v2(
                    &mut row,
                    clamped_members.iter().map(JointFeatureUseV2::canonical_row),
                );
                append_rows_v2(
                    &mut row,
                    fastener_stack
                        .iter()
                        .map(BoltStackParticipantV2::canonical_row),
                );
                preload.append_canonical(&mut row);
            }
            Self::Weld { members } => {
                append_rows_v2(
                    &mut row,
                    members.iter().map(JointFeatureUseV2::canonical_row),
                );
            }
            Self::AdhesiveBond { adherends } => {
                append_rows_v2(
                    &mut row,
                    adherends.iter().map(JointFeatureUseV2::canonical_row),
                );
            }
            Self::Key { shaft, hub, key } => {
                append_bytes_v2(&mut row, &shaft.canonical_row());
                append_bytes_v2(&mut row, &hub.canonical_row());
                append_bytes_v2(&mut row, &key.canonical_row());
            }
            Self::Spline { external, internal } | Self::InterferenceFit { external, internal } => {
                append_bytes_v2(&mut row, &external.canonical_row());
                append_bytes_v2(&mut row, &internal.canonical_row());
            }
        }
        row
    }
}

/// Exact typed role of a participant within one family payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AssemblyParticipantRoleV2 {
    /// Canonically indexed unordered clamped member.
    BoltClampedMember {
        /// Index after canonical unordered-member sorting.
        canonical_index: usize,
    },
    /// Physically ordered fastener-stack element.
    BoltStack {
        /// Zero-based bolt-head-to-thread-end position.
        position: u16,
        /// Closed stack role.
        role: BoltStackRoleV2,
    },
    /// Canonically indexed unordered weld member.
    WeldMember {
        /// Index after canonical unordered-member sorting.
        canonical_index: usize,
    },
    /// Canonically indexed unordered adhesive adherend.
    AdhesiveAdherend {
        /// Index after canonical unordered-member sorting.
        canonical_index: usize,
    },
    /// Shaft side of a keyed joint.
    KeyShaft,
    /// Hub side of a keyed joint.
    KeyHub,
    /// Distinct key body.
    KeyBody,
    /// External spline side.
    SplineExternal,
    /// Internal spline side.
    SplineInternal,
    /// External interference-fit side.
    InterferenceExternal,
    /// Internal interference-fit side.
    InterferenceInternal,
}

/// One family-specific physical joint occurrence, independent of chronology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JointOccurrenceV2 {
    id: JointOccurrenceIdV2,
    topology: JointTopologyV2,
    lifecycle: AssemblyLifecycleV2,
}

impl JointOccurrenceV2 {
    /// Construct one authority-free occurrence declaration.
    #[must_use]
    pub fn new(
        id: JointOccurrenceIdV2,
        topology: JointTopologyV2,
        lifecycle: AssemblyLifecycleV2,
    ) -> Self {
        Self {
            id,
            topology,
            lifecycle,
        }
    }

    /// Stable physical occurrence identity.
    #[must_use]
    pub const fn id(&self) -> &JointOccurrenceIdV2 {
        &self.id
    }

    /// Closed family-specific physical topology.
    #[must_use]
    pub const fn topology(&self) -> &JointTopologyV2 {
        &self.topology
    }

    /// Truthful planned or execution-claimed lifecycle.
    #[must_use]
    pub const fn lifecycle(&self) -> &AssemblyLifecycleV2 {
        &self.lifecycle
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(2_048);
        append_bytes_v2(&mut row, self.id.canonical_key().as_bytes());
        append_bytes_v2(&mut row, &self.topology.canonical_row());
        append_bytes_v2(&mut row, &self.lifecycle.canonical_row());
        row
    }
}

/// One chronological, atomic availability transition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyStepV2 {
    id: AssemblyStepIdV2,
    ordinal: u32,
    introduced_bodies: Vec<BodyId>,
    occurrence_ids: Vec<JointOccurrenceIdV2>,
}

impl AssemblyStepV2 {
    /// Construct one authority-free chronological step.
    #[must_use]
    pub fn new(
        id: AssemblyStepIdV2,
        ordinal: u32,
        introduced_bodies: Vec<BodyId>,
        occurrence_ids: Vec<JointOccurrenceIdV2>,
    ) -> Self {
        Self {
            id,
            ordinal,
            introduced_bodies,
            occurrence_ids,
        }
    }

    /// Stable chronological step identity.
    #[must_use]
    pub const fn id(&self) -> &AssemblyStepIdV2 {
        &self.id
    }

    /// Zero-based total-order position.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Bodies introduced atomically after this whole step validates.
    #[must_use]
    pub fn introduced_bodies(&self) -> &[BodyId] {
        &self.introduced_bodies
    }

    /// Physical occurrences scheduled atomically at this step.
    #[must_use]
    pub fn occurrence_ids(&self) -> &[JointOccurrenceIdV2] {
        &self.occurrence_ids
    }
}

/// Mutable-by-construction V2 assembly draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineAssemblyDraftV2 {
    /// Bodies available before chronological step zero; caller order is not semantic.
    pub initial_available_bodies: Vec<BodyId>,
    /// Chronological steps; caller collection order is not semantic.
    pub steps: Vec<AssemblyStepV2>,
    /// Physical occurrences; caller collection order is not semantic.
    pub occurrences: Vec<JointOccurrenceV2>,
}

impl MachineAssemblyDraftV2 {
    /// Admit topology and chronology against one exact Machine graph.
    ///
    /// # Errors
    /// Refuses resource overflow, graph/ownership gaps, family invariant
    /// violations, duplicate occurrence/use identities, ambiguous exclusivity,
    /// and non-atomic availability histories.
    #[allow(clippy::result_large_err)]
    pub fn admit_against(
        self,
        graph: &AdmittedMachineGraph,
    ) -> Result<AdmittedMachineAssemblyV2, MachineAssemblyAdmissionErrorV2> {
        admit_assembly_v2(self, graph)
    }
}

/// Canonical identity schema for a graph-bound V2 assembly declaration.
pub enum MachineAssemblyIdentitySchemaV2 {}

impl CanonicalSchema for MachineAssemblyIdentitySchemaV2 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.manufacturing-assembly.v2";
    const NAME: &'static str = "admitted-machine-assembly";
    const VERSION: u32 = MACHINE_ASSEMBLY_SCHEMA_VERSION_V2;
    const CONTEXT: &'static str = "one exact Machine graph, initial availability set, family-specific physical occurrences, and versioned authenticated chronological availability-transition chain";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("assembly-schema-version", WireType::U64),
        FieldSpec::required("frankenscript-ir-version", WireType::U64),
        FieldSpec::required("machine-graph", WireType::Bytes),
        FieldSpec::required("initial-available-bodies", WireType::OrderedBytes),
        FieldSpec::required("joint-occurrences", WireType::OrderedBytes),
        FieldSpec::required("assembly-steps", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of one admitted V2 assembly declaration.
pub type MachineAssemblyIdV2 = ProblemSemanticId<MachineAssemblyIdentitySchemaV2>;

/// One validated transition in the replayable authenticated availability chain.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedAssemblyStepV2 {
    step: AssemblyStepV2,
    available_before_count: usize,
    available_after_count: usize,
    availability_before_root: ContentId,
    availability_after_root: ContentId,
}

impl AdmittedAssemblyStepV2 {
    /// Original step after deterministic canonical sorting.
    #[must_use]
    pub const fn step(&self) -> &AssemblyStepV2 {
        &self.step
    }

    /// Number of available bodies before the step.
    #[must_use]
    pub const fn available_before_count(&self) -> usize {
        self.available_before_count
    }

    /// Number of available bodies after full step validation.
    #[must_use]
    pub const fn available_after_count(&self) -> usize {
        self.available_after_count
    }

    /// History-bound root immediately before this step.
    ///
    /// For step zero this is the canonical initial-set root. For every later
    /// step it is exactly the preceding step's `availability_after_root`.
    #[must_use]
    pub const fn availability_before_root(&self) -> ContentId {
        self.availability_before_root
    }

    /// History-bound root after authenticating this exact transition.
    ///
    /// This is not a history-independent set digest. A verifier recomputes it
    /// from the prior root/count, step ID/ordinal, sorted introductions, and
    /// after-count under the versioned V2 transition domain.
    #[must_use]
    pub const fn availability_after_root(&self) -> ContentId {
        self.availability_after_root
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(2_048);
        append_bytes_v2(&mut row, self.step.id.canonical_key().as_bytes());
        row.extend_from_slice(&self.step.ordinal.to_le_bytes());
        append_rows_v2(
            &mut row,
            self.step.introduced_bodies.iter().map(body_row_v2),
        );
        append_rows_v2(
            &mut row,
            self.step
                .occurrence_ids
                .iter()
                .map(|id| id_key_row_v2(id.canonical_key())),
        );
        row.extend_from_slice(&(self.available_before_count as u64).to_le_bytes());
        append_bytes_v2(&mut row, self.availability_before_root.as_bytes());
        row.extend_from_slice(&(self.available_after_count as u64).to_le_bytes());
        append_bytes_v2(&mut row, self.availability_after_root.as_bytes());
        row
    }
}

/// Canonically ordered graph-bound topology and chronology plus full receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineAssemblyV2 {
    admitted_against_graph: MachineGraphIdV1,
    initial_available_bodies: Vec<BodyId>,
    initial_availability_root: ContentId,
    occurrences: Vec<JointOccurrenceV2>,
    steps: Vec<AdmittedAssemblyStepV2>,
    receipt: IdentityReceipt<MachineAssemblyIdV2>,
}

impl AdmittedMachineAssemblyV2 {
    /// Exact Machine graph this declaration was admitted against.
    #[must_use]
    pub const fn admitted_against_graph(&self) -> MachineGraphIdV1 {
        self.admitted_against_graph
    }

    /// Canonically sorted bodies available before step zero.
    #[must_use]
    pub fn initial_available_bodies(&self) -> &[BodyId] {
        &self.initial_available_bodies
    }

    /// Canonical-set root that seeds the history-bound transition chain.
    #[must_use]
    pub const fn initial_availability_root(&self) -> ContentId {
        self.initial_availability_root
    }

    /// Physical occurrences in stable occurrence-ID order.
    #[must_use]
    pub fn occurrences(&self) -> &[JointOccurrenceV2] {
        &self.occurrences
    }

    /// Chronological steps in checked ordinal order.
    #[must_use]
    pub fn steps(&self) -> &[AdmittedAssemblyStepV2] {
        &self.steps
    }

    /// Domain-separated aggregate identity.
    #[must_use]
    pub const fn identity(&self) -> MachineAssemblyIdV2 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineAssemblyIdV2> {
        self.receipt
    }
}

/// Exact family-specific topology defect.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssemblyTopologyIssueV2 {
    /// A family-specific unordered member set was too small.
    TooFewMembers {
        /// Closed joint family.
        family: AssemblyJointFamilyV2,
        /// Submitted member count.
        actual: usize,
        /// Minimum admitted member count.
        min: usize,
    },
    /// Total typed participants exceeded the per-occurrence cap.
    ParticipantLimit {
        /// Submitted count before nested validation.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// A preloaded-bolt topology did not contain exactly one bolt.
    BoltCount {
        /// Submitted number of bolt-role stack participants.
        actual: usize,
    },
    /// A preloaded-bolt topology contained more than one nut.
    NutCount {
        /// Submitted number of nut-role stack participants.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// Two fastener-stack participants declared one physical position.
    DuplicateStackPosition {
        /// Repeated zero-based position.
        position: u16,
    },
    /// Fastener-stack positions were not exactly zero through N minus one.
    StackPositionGap {
        /// Required zero-based position.
        expected: u16,
        /// Submitted position at the sorted index.
        actual: u16,
    },
    /// Two typed roles in one occurrence declared one body.
    DuplicateParticipantBody {
        /// Repeated body.
        body: BodyId,
        /// First canonical role.
        first_role: AssemblyParticipantRoleV2,
        /// Later canonical role.
        duplicate_role: AssemblyParticipantRoleV2,
    },
    /// Two typed roles in one occurrence selected one physical feature.
    DuplicatePhysicalFeature {
        /// Repeated physical feature.
        feature: ContactFeatureId,
        /// First canonical role.
        first_role: AssemblyParticipantRoleV2,
        /// Later canonical role.
        duplicate_role: AssemblyParticipantRoleV2,
    },
}

/// Structured refusal from V2 assembly admission.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineAssemblyAdmissionErrorV2 {
    /// No body was declared available before step zero.
    EmptyInitialAvailability,
    /// Raw initial availability exceeded the fixed cap.
    InitialBodyLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One initial body appeared more than once.
    DuplicateInitialBody {
        /// Repeated body.
        body: BodyId,
    },
    /// An initially available body is absent from the admitted graph.
    UnknownInitialBody {
        /// Missing body.
        body: BodyId,
    },
    /// At least one chronological step is required.
    NoSteps,
    /// Raw chronological steps exceeded the fixed cap.
    StepLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One step identity appeared more than once.
    DuplicateStep {
        /// Repeated step identity.
        step: AssemblyStepIdV2,
    },
    /// Two steps declared one ordinal.
    DuplicateOrdinal {
        /// Repeated ordinal.
        ordinal: u32,
        /// Lexically first step at that ordinal.
        first: AssemblyStepIdV2,
        /// Later step at that ordinal.
        duplicate: AssemblyStepIdV2,
    },
    /// Sorted ordinals were not exactly zero through N minus one.
    OrdinalGap {
        /// Step exposing the gap.
        step: AssemblyStepIdV2,
        /// Required zero-based ordinal.
        expected: u32,
        /// Submitted ordinal.
        actual: u32,
    },
    /// One step introduced too many bodies.
    StepIntroductionLimit {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Submitted count before nested validation.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One step scheduled no physical occurrence.
    StepWithoutOccurrence {
        /// Invalid step.
        step: AssemblyStepIdV2,
    },
    /// One step scheduled too many physical occurrences.
    StepOccurrenceLimit {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Submitted count before nested validation.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One body appeared twice in a step's introduction set.
    DuplicateIntroducedBody {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Repeated body.
        body: BodyId,
    },
    /// A step introduced a body absent from the admitted graph.
    UnknownIntroducedBody {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Missing body.
        body: BodyId,
    },
    /// A step attempted to reintroduce an already available body.
    BodyAlreadyAvailable {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Reintroduced body.
        body: BodyId,
    },
    /// At least one physical occurrence is required.
    NoOccurrences,
    /// Raw physical occurrences exceeded the fixed cap.
    OccurrenceLimit {
        /// Submitted count before sorting or deduplication.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One physical occurrence identity appeared more than once.
    DuplicateOccurrence {
        /// Repeated occurrence identity.
        occurrence: JointOccurrenceIdV2,
    },
    /// A family payload violated a closed topology invariant.
    InvalidTopology {
        /// Invalid physical occurrence.
        occurrence: JointOccurrenceIdV2,
        /// Exact family-specific defect.
        issue: AssemblyTopologyIssueV2,
    },
    /// One occurrence-local feature-use identity appeared more than once.
    DuplicateFeatureUse {
        /// Repeated use identity.
        feature_use: JointFeatureUseIdV2,
        /// First physical occurrence containing it.
        first: JointOccurrenceIdV2,
        /// Later physical occurrence containing it.
        duplicate: JointOccurrenceIdV2,
    },
    /// A participant named a body absent from the admitted graph.
    UnknownParticipantBody {
        /// Invalid occurrence.
        occurrence: JointOccurrenceIdV2,
        /// Exact family role.
        role: AssemblyParticipantRoleV2,
        /// Missing body.
        body: BodyId,
    },
    /// A participant named a feature absent from the admitted graph.
    UnknownParticipantFeature {
        /// Invalid occurrence.
        occurrence: JointOccurrenceIdV2,
        /// Exact family role.
        role: AssemblyParticipantRoleV2,
        /// Missing feature.
        feature: ContactFeatureId,
    },
    /// Participant body and feature exist under different subsystem owners.
    ParticipantOwnerMismatch {
        /// Invalid occurrence.
        occurrence: JointOccurrenceIdV2,
        /// Exact family role.
        role: AssemblyParticipantRoleV2,
        /// Caller-declared body.
        body: BodyId,
        /// Selected physical feature.
        feature: ContactFeatureId,
        /// Graph owner of the body.
        body_owner: SubsystemId,
        /// Graph owner of the feature.
        feature_owner: SubsystemId,
    },
    /// One durable feature was associated with conflicting declared bodies.
    ConflictingFeatureBodyDeclaration {
        /// Physical feature with inconsistent declarations.
        feature: ContactFeatureId,
        /// First canonical declared body.
        first_body: BodyId,
        /// Later conflicting declared body.
        conflicting_body: BodyId,
    },
    /// A use declared exclusivity while another use selected the same feature.
    ExclusiveFeatureReuse {
        /// Exclusively claimed physical feature.
        feature: ContactFeatureId,
        /// First use selecting it.
        first: JointFeatureUseIdV2,
        /// Conflicting later use.
        duplicate: JointFeatureUseIdV2,
    },
    /// One occurrence reference appeared twice within a step.
    DuplicateOccurrenceReference {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Repeated occurrence reference.
        occurrence: JointOccurrenceIdV2,
    },
    /// A step referenced no declared physical occurrence.
    UnknownOccurrenceReference {
        /// Invalid step.
        step: AssemblyStepIdV2,
        /// Missing occurrence.
        occurrence: JointOccurrenceIdV2,
    },
    /// One physical occurrence was scheduled by multiple steps.
    OccurrenceScheduledTwice {
        /// Repeated physical occurrence.
        occurrence: JointOccurrenceIdV2,
        /// First chronological step.
        first: AssemblyStepIdV2,
        /// Later chronological step.
        duplicate: AssemblyStepIdV2,
    },
    /// A scheduled occurrence selected a body unavailable at that atomic step.
    ParticipantBodyUnavailable {
        /// Invalid chronological step.
        step: AssemblyStepIdV2,
        /// Scheduled physical occurrence.
        occurrence: JointOccurrenceIdV2,
        /// Exact family role.
        role: AssemblyParticipantRoleV2,
        /// Unavailable body.
        body: BodyId,
    },
    /// A newly introduced body did not structurally participate in the step.
    IntroducedBodyDoesNotParticipate {
        /// Invalid chronological step.
        step: AssemblyStepIdV2,
        /// Unused introduced body.
        body: BodyId,
    },
    /// A declared physical occurrence was never scheduled.
    UnscheduledOccurrence {
        /// Unscheduled occurrence.
        occurrence: JointOccurrenceIdV2,
    },
    /// Canonical aggregate identity publication failed.
    Identity(CanonicalError),
}

impl MachineAssemblyAdmissionErrorV2 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::EmptyInitialAvailability => "MachineAssemblyEmptyInitialAvailability",
            Self::InitialBodyLimit { .. } => "MachineAssemblyInitialBodyLimit",
            Self::DuplicateInitialBody { .. } => "MachineAssemblyDuplicateInitialBody",
            Self::UnknownInitialBody { .. } => "MachineAssemblyUnknownInitialBody",
            Self::NoSteps => "MachineAssemblyNoSteps",
            Self::StepLimit { .. } => "MachineAssemblyStepLimit",
            Self::DuplicateStep { .. } => "MachineAssemblyDuplicateStep",
            Self::DuplicateOrdinal { .. } => "MachineAssemblyDuplicateOrdinal",
            Self::OrdinalGap { .. } => "MachineAssemblyOrdinalGap",
            Self::StepIntroductionLimit { .. } => "MachineAssemblyStepIntroductionLimit",
            Self::StepWithoutOccurrence { .. } => "MachineAssemblyStepWithoutOccurrence",
            Self::StepOccurrenceLimit { .. } => "MachineAssemblyStepOccurrenceLimit",
            Self::DuplicateIntroducedBody { .. } => "MachineAssemblyDuplicateIntroducedBody",
            Self::UnknownIntroducedBody { .. } => "MachineAssemblyUnknownIntroducedBody",
            Self::BodyAlreadyAvailable { .. } => "MachineAssemblyBodyAlreadyAvailable",
            Self::NoOccurrences => "MachineAssemblyNoOccurrences",
            Self::OccurrenceLimit { .. } => "MachineAssemblyOccurrenceLimit",
            Self::DuplicateOccurrence { .. } => "MachineAssemblyDuplicateOccurrence",
            Self::InvalidTopology { .. } => "MachineAssemblyInvalidTopology",
            Self::DuplicateFeatureUse { .. } => "MachineAssemblyDuplicateFeatureUse",
            Self::UnknownParticipantBody { .. } => "MachineAssemblyUnknownParticipantBody",
            Self::UnknownParticipantFeature { .. } => "MachineAssemblyUnknownParticipantFeature",
            Self::ParticipantOwnerMismatch { .. } => "MachineAssemblyParticipantOwnerMismatch",
            Self::ConflictingFeatureBodyDeclaration { .. } => {
                "MachineAssemblyConflictingFeatureBodyDeclaration"
            }
            Self::ExclusiveFeatureReuse { .. } => "MachineAssemblyExclusiveFeatureReuse",
            Self::DuplicateOccurrenceReference { .. } => {
                "MachineAssemblyDuplicateOccurrenceReference"
            }
            Self::UnknownOccurrenceReference { .. } => "MachineAssemblyUnknownOccurrenceReference",
            Self::OccurrenceScheduledTwice { .. } => "MachineAssemblyOccurrenceScheduledTwice",
            Self::ParticipantBodyUnavailable { .. } => "MachineAssemblyParticipantBodyUnavailable",
            Self::IntroducedBodyDoesNotParticipate { .. } => {
                "MachineAssemblyIntroducedBodyDoesNotParticipate"
            }
            Self::UnscheduledOccurrence { .. } => "MachineAssemblyUnscheduledOccurrence",
            Self::Identity(_) => "MachineAssemblyIdentity",
        }
    }
}

impl fmt::Display for MachineAssemblyAdmissionErrorV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for MachineAssemblyAdmissionErrorV2 {}

impl From<CanonicalError> for MachineAssemblyAdmissionErrorV2 {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

#[derive(Debug)]
struct GraphOwnersV2 {
    body: BTreeMap<BodyId, SubsystemId>,
    feature: BTreeMap<ContactFeatureId, SubsystemId>,
}

impl GraphOwnersV2 {
    fn from_graph(graph: &AdmittedMachineGraph) -> Self {
        let body = graph
            .subsystems()
            .iter()
            .flat_map(|subsystem| {
                subsystem
                    .bodies
                    .iter()
                    .cloned()
                    .map(move |body| (body, subsystem.id.clone()))
            })
            .collect();
        let feature = graph
            .subsystems()
            .iter()
            .flat_map(|subsystem| {
                subsystem
                    .contact_features
                    .iter()
                    .cloned()
                    .map(move |feature| (feature, subsystem.id.clone()))
            })
            .collect();
        Self { body, feature }
    }
}

fn validate_topology_shape_v2(
    occurrence: &JointOccurrenceV2,
) -> Result<(), MachineAssemblyAdmissionErrorV2> {
    let participants = occurrence.topology.participants();
    if participants.len() > MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2 {
        return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
            occurrence: occurrence.id.clone(),
            issue: AssemblyTopologyIssueV2::ParticipantLimit {
                actual: participants.len(),
                max: MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2,
            },
        });
    }

    match &occurrence.topology {
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            ..
        } => {
            if clamped_members.len() < 2 {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::TooFewMembers {
                        family: AssemblyJointFamilyV2::PreloadedBolt,
                        actual: clamped_members.len(),
                        min: 2,
                    },
                });
            }
            let bolt_count = fastener_stack
                .iter()
                .filter(|participant| participant.role == BoltStackRoleV2::Bolt)
                .count();
            if bolt_count != 1 {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::BoltCount { actual: bolt_count },
                });
            }
            let nut_count = fastener_stack
                .iter()
                .filter(|participant| participant.role == BoltStackRoleV2::Nut)
                .count();
            if nut_count > 1 {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::NutCount {
                        actual: nut_count,
                        max: 1,
                    },
                });
            }
            if let Some(pair) = fastener_stack
                .windows(2)
                .find(|pair| pair[0].position == pair[1].position)
            {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::DuplicateStackPosition {
                        position: pair[0].position,
                    },
                });
            }
            for (index, participant) in fastener_stack.iter().enumerate() {
                let expected = u16::try_from(index).map_err(|_| {
                    MachineAssemblyAdmissionErrorV2::InvalidTopology {
                        occurrence: occurrence.id.clone(),
                        issue: AssemblyTopologyIssueV2::ParticipantLimit {
                            actual: participants.len(),
                            max: MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2,
                        },
                    }
                })?;
                if participant.position != expected {
                    return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                        occurrence: occurrence.id.clone(),
                        issue: AssemblyTopologyIssueV2::StackPositionGap {
                            expected,
                            actual: participant.position,
                        },
                    });
                }
            }
        }
        JointTopologyV2::Weld { members } => {
            if members.len() < 2 {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::TooFewMembers {
                        family: AssemblyJointFamilyV2::Weld,
                        actual: members.len(),
                        min: 2,
                    },
                });
            }
        }
        JointTopologyV2::AdhesiveBond { adherends } => {
            if adherends.len() < 2 {
                return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                    occurrence: occurrence.id.clone(),
                    issue: AssemblyTopologyIssueV2::TooFewMembers {
                        family: AssemblyJointFamilyV2::AdhesiveBond,
                        actual: adherends.len(),
                        min: 2,
                    },
                });
            }
        }
        JointTopologyV2::Key { .. }
        | JointTopologyV2::Spline { .. }
        | JointTopologyV2::InterferenceFit { .. } => {}
    }

    let mut bodies = BTreeMap::<BodyId, AssemblyParticipantRoleV2>::new();
    let mut features = BTreeMap::<ContactFeatureId, AssemblyParticipantRoleV2>::new();
    for (role, feature_use) in participants {
        let selector = feature_use.selector();
        if let Some(first_role) = bodies.insert(selector.declared_body().clone(), role) {
            return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                occurrence: occurrence.id.clone(),
                issue: AssemblyTopologyIssueV2::DuplicateParticipantBody {
                    body: selector.declared_body().clone(),
                    first_role,
                    duplicate_role: role,
                },
            });
        }
        if let Some(first_role) = features.insert(selector.contact_feature().clone(), role) {
            return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                occurrence: occurrence.id.clone(),
                issue: AssemblyTopologyIssueV2::DuplicatePhysicalFeature {
                    feature: selector.contact_feature().clone(),
                    first_role,
                    duplicate_role: role,
                },
            });
        }
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct FirstPhysicalUseV2 {
    declared_body: BodyId,
    use_id: JointFeatureUseIdV2,
    policy: PhysicalFeatureUsePolicyV2,
}

fn validate_occurrences_v2(
    occurrences: &[JointOccurrenceV2],
    owners: &GraphOwnersV2,
) -> Result<(), MachineAssemblyAdmissionErrorV2> {
    let mut use_ids = BTreeMap::<JointFeatureUseIdV2, JointOccurrenceIdV2>::new();
    let mut physical_uses = BTreeMap::<ContactFeatureId, FirstPhysicalUseV2>::new();

    for occurrence in occurrences {
        validate_topology_shape_v2(occurrence)?;
        for (role, feature_use) in occurrence.topology.participants() {
            if let Some(first) = use_ids.insert(feature_use.id.clone(), occurrence.id.clone()) {
                return Err(MachineAssemblyAdmissionErrorV2::DuplicateFeatureUse {
                    feature_use: feature_use.id.clone(),
                    first,
                    duplicate: occurrence.id.clone(),
                });
            }

            let selector = feature_use.selector();
            let Some(body_owner) = owners.body.get(selector.declared_body()) else {
                return Err(MachineAssemblyAdmissionErrorV2::UnknownParticipantBody {
                    occurrence: occurrence.id.clone(),
                    role,
                    body: selector.declared_body().clone(),
                });
            };
            let Some(feature_owner) = owners.feature.get(selector.contact_feature()) else {
                return Err(MachineAssemblyAdmissionErrorV2::UnknownParticipantFeature {
                    occurrence: occurrence.id.clone(),
                    role,
                    feature: selector.contact_feature().clone(),
                });
            };
            if body_owner != feature_owner {
                return Err(MachineAssemblyAdmissionErrorV2::ParticipantOwnerMismatch {
                    occurrence: occurrence.id.clone(),
                    role,
                    body: selector.declared_body().clone(),
                    feature: selector.contact_feature().clone(),
                    body_owner: body_owner.clone(),
                    feature_owner: feature_owner.clone(),
                });
            }

            match physical_uses.get(selector.contact_feature()) {
                Some(first) if first.declared_body != *selector.declared_body() => {
                    return Err(
                        MachineAssemblyAdmissionErrorV2::ConflictingFeatureBodyDeclaration {
                            feature: selector.contact_feature().clone(),
                            first_body: first.declared_body.clone(),
                            conflicting_body: selector.declared_body().clone(),
                        },
                    );
                }
                Some(first)
                    if first.policy == PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly
                        || feature_use.policy
                            == PhysicalFeatureUsePolicyV2::ExclusiveWithinAssembly =>
                {
                    return Err(MachineAssemblyAdmissionErrorV2::ExclusiveFeatureReuse {
                        feature: selector.contact_feature().clone(),
                        first: first.use_id.clone(),
                        duplicate: feature_use.id.clone(),
                    });
                }
                Some(_) => {}
                None => {
                    physical_uses.insert(
                        selector.contact_feature().clone(),
                        FirstPhysicalUseV2 {
                            declared_body: selector.declared_body().clone(),
                            use_id: feature_use.id.clone(),
                            policy: feature_use.policy,
                        },
                    );
                }
            }
        }
    }
    Ok(())
}

fn raw_participant_count_v2(topology: &JointTopologyV2) -> usize {
    match topology {
        JointTopologyV2::PreloadedBolt {
            clamped_members,
            fastener_stack,
            ..
        } => clamped_members.len().saturating_add(fastener_stack.len()),
        JointTopologyV2::Weld { members } => members.len(),
        JointTopologyV2::AdhesiveBond { adherends } => adherends.len(),
        JointTopologyV2::Key { .. } => 3,
        JointTopologyV2::Spline { .. } | JointTopologyV2::InterferenceFit { .. } => 2,
    }
}

#[allow(clippy::too_many_lines)]
#[allow(clippy::result_large_err)]
fn admit_assembly_v2(
    draft: MachineAssemblyDraftV2,
    graph: &AdmittedMachineGraph,
) -> Result<AdmittedMachineAssemblyV2, MachineAssemblyAdmissionErrorV2> {
    if draft.initial_available_bodies.is_empty() {
        return Err(MachineAssemblyAdmissionErrorV2::EmptyInitialAvailability);
    }
    if draft.initial_available_bodies.len() > MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2 {
        return Err(MachineAssemblyAdmissionErrorV2::InitialBodyLimit {
            actual: draft.initial_available_bodies.len(),
            max: MAX_MACHINE_ASSEMBLY_INITIAL_BODIES_V2,
        });
    }
    if draft.steps.is_empty() {
        return Err(MachineAssemblyAdmissionErrorV2::NoSteps);
    }
    if draft.steps.len() > MAX_MACHINE_ASSEMBLY_STEPS_V2 {
        return Err(MachineAssemblyAdmissionErrorV2::StepLimit {
            actual: draft.steps.len(),
            max: MAX_MACHINE_ASSEMBLY_STEPS_V2,
        });
    }
    if draft.occurrences.is_empty() {
        return Err(MachineAssemblyAdmissionErrorV2::NoOccurrences);
    }
    if draft.occurrences.len() > MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2 {
        return Err(MachineAssemblyAdmissionErrorV2::OccurrenceLimit {
            actual: draft.occurrences.len(),
            max: MAX_MACHINE_ASSEMBLY_OCCURRENCES_V2,
        });
    }
    let mut initial_available_bodies = draft.initial_available_bodies;
    initial_available_bodies.sort();
    if let Some(pair) = initial_available_bodies
        .windows(2)
        .find(|pair| pair[0] == pair[1])
    {
        return Err(MachineAssemblyAdmissionErrorV2::DuplicateInitialBody {
            body: pair[0].clone(),
        });
    }
    let mut occurrences = draft.occurrences;
    occurrences.sort_by(|left, right| left.id.cmp(&right.id));
    if let Some(pair) = occurrences.windows(2).find(|pair| pair[0].id == pair[1].id) {
        return Err(MachineAssemblyAdmissionErrorV2::DuplicateOccurrence {
            occurrence: pair[0].id.clone(),
        });
    }
    for occurrence in &occurrences {
        let raw_count = raw_participant_count_v2(&occurrence.topology);
        if raw_count > MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2 {
            return Err(MachineAssemblyAdmissionErrorV2::InvalidTopology {
                occurrence: occurrence.id.clone(),
                issue: AssemblyTopologyIssueV2::ParticipantLimit {
                    actual: raw_count,
                    max: MAX_MACHINE_ASSEMBLY_PARTICIPANTS_PER_OCCURRENCE_V2,
                },
            });
        }
    }
    let mut steps = draft.steps;
    steps.sort_by(|left, right| {
        left.ordinal
            .cmp(&right.ordinal)
            .then_with(|| left.id.cmp(&right.id))
    });
    let mut step_ids = BTreeSet::<AssemblyStepIdV2>::new();
    for step in &steps {
        if !step_ids.insert(step.id.clone()) {
            return Err(MachineAssemblyAdmissionErrorV2::DuplicateStep {
                step: step.id.clone(),
            });
        }
    }
    if let Some(pair) = steps
        .windows(2)
        .find(|pair| pair[0].ordinal == pair[1].ordinal)
    {
        return Err(MachineAssemblyAdmissionErrorV2::DuplicateOrdinal {
            ordinal: pair[0].ordinal,
            first: pair[0].id.clone(),
            duplicate: pair[1].id.clone(),
        });
    }
    for (index, step) in steps.iter().enumerate() {
        let expected =
            u32::try_from(index).map_err(|_| MachineAssemblyAdmissionErrorV2::StepLimit {
                actual: steps.len(),
                max: MAX_MACHINE_ASSEMBLY_STEPS_V2,
            })?;
        if step.ordinal != expected {
            return Err(MachineAssemblyAdmissionErrorV2::OrdinalGap {
                step: step.id.clone(),
                expected,
                actual: step.ordinal,
            });
        }
    }
    for step in &steps {
        if step.introduced_bodies.len() > MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2 {
            return Err(MachineAssemblyAdmissionErrorV2::StepIntroductionLimit {
                step: step.id.clone(),
                actual: step.introduced_bodies.len(),
                max: MAX_MACHINE_ASSEMBLY_INTRODUCTIONS_PER_STEP_V2,
            });
        }
        if step.occurrence_ids.is_empty() {
            return Err(MachineAssemblyAdmissionErrorV2::StepWithoutOccurrence {
                step: step.id.clone(),
            });
        }
        if step.occurrence_ids.len() > MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2 {
            return Err(MachineAssemblyAdmissionErrorV2::StepOccurrenceLimit {
                step: step.id.clone(),
                actual: step.occurrence_ids.len(),
                max: MAX_MACHINE_ASSEMBLY_OCCURRENCES_PER_STEP_V2,
            });
        }
    }

    // Complete every raw nested-cardinality preflight before traversing the
    // admitted graph or allocating canonical participant rows. Sorting only
    // durable IDs/ordinals first makes the selected refusal caller-order
    // invariant when several submitted collections are independently invalid.
    let owners = GraphOwnersV2::from_graph(graph);
    for body in &initial_available_bodies {
        if !owners.body.contains_key(body) {
            return Err(MachineAssemblyAdmissionErrorV2::UnknownInitialBody { body: body.clone() });
        }
    }
    for occurrence in &mut occurrences {
        occurrence.topology.canonicalize_unordered_members();
    }
    validate_occurrences_v2(&occurrences, &owners)?;

    let occurrence_index = occurrences
        .iter()
        .enumerate()
        .map(|(index, occurrence)| (occurrence.id.clone(), index))
        .collect::<BTreeMap<_, _>>();

    for step in &mut steps {
        step.introduced_bodies.sort();
        step.occurrence_ids.sort();
    }

    let mut available = initial_available_bodies
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    let initial_availability_root = initial_availability_root_v2(&initial_available_bodies);
    let mut availability_root = initial_availability_root;
    let mut scheduled = BTreeMap::<JointOccurrenceIdV2, AssemblyStepIdV2>::new();
    let mut admitted_steps = Vec::with_capacity(steps.len());

    for step in steps {
        if let Some(pair) = step
            .introduced_bodies
            .windows(2)
            .find(|pair| pair[0] == pair[1])
        {
            return Err(MachineAssemblyAdmissionErrorV2::DuplicateIntroducedBody {
                step: step.id.clone(),
                body: pair[0].clone(),
            });
        }
        if let Some(pair) = step
            .occurrence_ids
            .windows(2)
            .find(|pair| pair[0] == pair[1])
        {
            return Err(
                MachineAssemblyAdmissionErrorV2::DuplicateOccurrenceReference {
                    step: step.id.clone(),
                    occurrence: pair[0].clone(),
                },
            );
        }

        let available_before_count = available.len();
        let introduced = step
            .introduced_bodies
            .iter()
            .cloned()
            .collect::<BTreeSet<_>>();
        for body in &step.introduced_bodies {
            if !owners.body.contains_key(body) {
                return Err(MachineAssemblyAdmissionErrorV2::UnknownIntroducedBody {
                    step: step.id.clone(),
                    body: body.clone(),
                });
            }
            if available.contains(body) {
                return Err(MachineAssemblyAdmissionErrorV2::BodyAlreadyAvailable {
                    step: step.id.clone(),
                    body: body.clone(),
                });
            }
        }

        let mut participating_bodies = BTreeSet::<BodyId>::new();
        let mut pending_occurrences = Vec::<JointOccurrenceIdV2>::new();
        for occurrence_id in &step.occurrence_ids {
            let Some(index) = occurrence_index.get(occurrence_id) else {
                return Err(
                    MachineAssemblyAdmissionErrorV2::UnknownOccurrenceReference {
                        step: step.id.clone(),
                        occurrence: occurrence_id.clone(),
                    },
                );
            };
            if let Some(first) = scheduled.get(occurrence_id) {
                return Err(MachineAssemblyAdmissionErrorV2::OccurrenceScheduledTwice {
                    occurrence: occurrence_id.clone(),
                    first: first.clone(),
                    duplicate: step.id.clone(),
                });
            }
            let occurrence = &occurrences[*index];
            for (role, feature_use) in occurrence.topology.participants() {
                let body = feature_use.selector().declared_body();
                if !available.contains(body) && !introduced.contains(body) {
                    return Err(
                        MachineAssemblyAdmissionErrorV2::ParticipantBodyUnavailable {
                            step: step.id.clone(),
                            occurrence: occurrence_id.clone(),
                            role,
                            body: body.clone(),
                        },
                    );
                }
                participating_bodies.insert(body.clone());
            }
            pending_occurrences.push(occurrence_id.clone());
        }
        for body in &step.introduced_bodies {
            if !participating_bodies.contains(body) {
                return Err(
                    MachineAssemblyAdmissionErrorV2::IntroducedBodyDoesNotParticipate {
                        step: step.id.clone(),
                        body: body.clone(),
                    },
                );
            }
        }

        for occurrence_id in pending_occurrences {
            scheduled.insert(occurrence_id, step.id.clone());
        }
        let available_after_count = available_before_count + step.introduced_bodies.len();
        let availability_after_root = availability_transition_root_v2(
            availability_root,
            available_before_count,
            &step,
            available_after_count,
        );
        for body in &step.introduced_bodies {
            let inserted = available.insert(body.clone());
            debug_assert!(inserted, "validated introductions are new");
        }
        debug_assert_eq!(available.len(), available_after_count);
        admitted_steps.push(AdmittedAssemblyStepV2 {
            step,
            available_before_count,
            available_after_count,
            availability_before_root: availability_root,
            availability_after_root,
        });
        availability_root = availability_after_root;
    }

    if let Some(occurrence) = occurrences
        .iter()
        .find(|occurrence| !scheduled.contains_key(&occurrence.id))
    {
        return Err(MachineAssemblyAdmissionErrorV2::UnscheduledOccurrence {
            occurrence: occurrence.id.clone(),
        });
    }

    let initial_rows = initial_available_bodies
        .iter()
        .map(body_row_v2)
        .collect::<Vec<_>>();
    let occurrence_rows = occurrences
        .iter()
        .map(JointOccurrenceV2::canonical_row)
        .collect::<Vec<_>>();
    let step_rows = admitted_steps
        .iter()
        .map(AdmittedAssemblyStepV2::canonical_row)
        .collect::<Vec<_>>();
    let graph_id = graph.identity();
    let receipt = CanonicalEncoder::<MachineAssemblyIdV2, _>::new(
        MACHINE_ASSEMBLY_IDENTITY_LIMITS_V2,
        NeverCancel,
    )?
    .u64(
        Field::new(0, "assembly-schema-version"),
        u64::from(MACHINE_ASSEMBLY_SCHEMA_VERSION_V2),
    )?
    .u64(
        Field::new(1, "frankenscript-ir-version"),
        u64::from(IR_VERSION),
    )?
    .bytes(Field::new(2, "machine-graph"), graph_id.as_bytes())?
    .ordered_bytes(
        Field::new(3, "initial-available-bodies"),
        initial_rows.len() as u64,
        initial_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(4, "joint-occurrences"),
        occurrence_rows.len() as u64,
        occurrence_rows.iter().map(Vec::as_slice),
    )?
    .ordered_bytes(
        Field::new(5, "assembly-steps"),
        step_rows.len() as u64,
        step_rows.iter().map(Vec::as_slice),
    )?
    .finish()?;

    Ok(AdmittedMachineAssemblyV2 {
        admitted_against_graph: graph_id,
        initial_available_bodies,
        initial_availability_root,
        occurrences,
        steps: admitted_steps,
        receipt,
    })
}

fn body_row_v2(body: &BodyId) -> Vec<u8> {
    let mut row = Vec::with_capacity(176);
    append_body_v2(&mut row, body);
    row
}

fn initial_availability_root_v2(bodies: &[BodyId]) -> ContentId {
    let mut preimage = Vec::with_capacity(128 + bodies.len().saturating_mul(184));
    append_bytes_v2(
        &mut preimage,
        AVAILABILITY_INITIAL_ROOT_DOMAIN_V2.as_bytes(),
    );
    preimage.extend_from_slice(&MACHINE_ASSEMBLY_AVAILABILITY_COMMITMENT_VERSION_V2.to_le_bytes());
    append_rows_v2(&mut preimage, bodies.iter().map(body_row_v2));
    ContentId::of_bytes(&preimage)
}

fn availability_transition_root_v2(
    prior_root: ContentId,
    available_before_count: usize,
    step: &AssemblyStepV2,
    available_after_count: usize,
) -> ContentId {
    let mut preimage = Vec::with_capacity(256 + step.introduced_bodies.len().saturating_mul(184));
    append_bytes_v2(&mut preimage, AVAILABILITY_STEP_ROOT_DOMAIN_V2.as_bytes());
    preimage.extend_from_slice(&MACHINE_ASSEMBLY_AVAILABILITY_COMMITMENT_VERSION_V2.to_le_bytes());
    append_bytes_v2(&mut preimage, prior_root.as_bytes());
    preimage.extend_from_slice(&(available_before_count as u64).to_le_bytes());
    append_bytes_v2(&mut preimage, step.id.canonical_key().as_bytes());
    preimage.extend_from_slice(&step.ordinal.to_le_bytes());
    append_rows_v2(
        &mut preimage,
        step.introduced_bodies.iter().map(body_row_v2),
    );
    preimage.extend_from_slice(&(available_after_count as u64).to_le_bytes());
    ContentId::of_bytes(&preimage)
}

fn id_key_row_v2(key: &str) -> Vec<u8> {
    let mut row = Vec::with_capacity(136);
    append_bytes_v2(&mut row, key.as_bytes());
    row
}

fn append_body_v2(out: &mut Vec<u8>, body: &BodyId) {
    append_bytes_v2(out, body.identity().as_bytes());
    append_bytes_v2(out, body.canonical_key().as_bytes());
}

fn append_feature_v2(out: &mut Vec<u8>, feature: &ContactFeatureId) {
    append_bytes_v2(out, feature.identity().as_bytes());
    append_bytes_v2(out, feature.canonical_key().as_bytes());
}

fn append_artifact_v2(out: &mut Vec<u8>, artifact: &ManufacturingArtifactRefV1) {
    let mut row = Vec::with_capacity(176);
    append_bytes_v2(&mut row, artifact.namespace().as_bytes());
    row.extend_from_slice(&artifact.schema_version().get().to_le_bytes());
    row.extend_from_slice(artifact.content_hash().as_bytes());
    append_bytes_v2(out, &row);
}

fn append_rows_v2(out: &mut Vec<u8>, rows: impl std::iter::ExactSizeIterator<Item = Vec<u8>>) {
    out.extend_from_slice(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        append_bytes_v2(out, &row);
    }
}

fn append_bytes_v2(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}
