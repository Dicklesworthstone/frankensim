//! Durable Machine-IR entity identity and topology-lineage kernel.
//!
//! PR-1 establishes the durable entity and topology-lineage law: array
//! positions are never identity, entity roles are nominally distinct, and a
//! topology change may rebind an attachment only when its source has one
//! unambiguous target. PR-2 adds a dependency-neutral, versioned machine graph
//! whose subsystem, clock, terminal, port, relation, material, and interface
//! declarations are admitted before a semantic identity is published. PR-3's
//! [`semantics`] module admits a separate behavior overlay for state contracts,
//! conditions, motion, events/resets, tolerances, and dependence declarations
//! without changing the already-published graph-v1 identity grammar. PR-4's
//! [`assurance`] module then binds sensors and experiment lineage, operational
//! context and QoIs, hazards/faults, explicit accounting boundaries, and a
//! terminating fidelity-escalation policy to one exact admitted behavior.
//! PR-5's [`lowering`] module admits a one-way projection into one concrete,
//! validated `fs-scenario` artifact plus exact externally owned domain
//! artifacts and a complete stable-ID crosswalk. It deliberately does not
//! infer executable meaning from the opaque references admitted by PR-2--4.
//!
//! Runtime coupling, executable material/interface cards, controllers,
//! safety adjudication, physical validation, balance proofs, adaptive routing,
//! and inverse/equivalence claims remain outside this module's authority.

pub mod assurance;
pub mod lowering;
pub mod semantics;

use core::fmt;
use core::hash::{Hash, Hasher};
use core::num::NonZeroU64;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EntityId, Field, FieldSpec,
    IdentityReceipt, NeverCancel, ProblemSemanticId, SemanticId, StrongIdentity, WireType,
};
use fs_qty::Dims;
use fs_qty::semantic::{
    AngleDomain, CompositionBasis, QuantityKind, SemanticType, StrainBasis, StrainComponent,
    ValueForm,
};

/// Version of every durable Machine-IR entity-key schema in this module.
pub const MACHINE_ENTITY_ID_SCHEMA_VERSION_V1: u32 = 1;
/// Version of the canonical lineage-record and invalidation schemas.
pub const MACHINE_LINEAGE_SCHEMA_VERSION_V1: u32 = 1;
/// Version of the admitted Machine-IR graph schema.
pub const MACHINE_GRAPH_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum canonical bytes in one human-auditable entity key.
pub const MAX_MACHINE_ENTITY_KEY_BYTES: usize = 128;
/// Maximum source relations in one synchronous lineage record.
pub const MAX_LINEAGE_RELATIONS: usize = 4_096;
/// Maximum targets declared for one source entity.
pub const MAX_LINEAGE_TARGETS_PER_SOURCE: usize = 4_096;
/// Maximum total source-to-target endpoints in one lineage record.
pub const MAX_LINEAGE_ENDPOINTS: usize = 8_192;
/// Maximum dependent attachments considered by one lineage admission.
pub const MAX_LINEAGE_DEPENDENTS: usize = 4_096;
/// Maximum clock declarations in one graph draft.
pub const MAX_MACHINE_GRAPH_CLOCKS: usize = 1_024;
/// Maximum subsystem declarations in one graph draft.
pub const MAX_MACHINE_GRAPH_SUBSYSTEMS: usize = 1_024;
/// Maximum terminal declarations in one graph draft.
pub const MAX_MACHINE_GRAPH_TERMINALS: usize = 4_096;
/// Maximum port declarations in one graph draft.
pub const MAX_MACHINE_GRAPH_PORTS: usize = 2_048;
/// Maximum relation declarations in one graph draft.
pub const MAX_MACHINE_GRAPH_RELATIONS: usize = 8_192;
/// Maximum material bindings in one graph draft.
pub const MAX_MACHINE_GRAPH_MATERIALS: usize = 4_096;
/// Maximum interface bindings in one graph draft.
pub const MAX_MACHINE_GRAPH_INTERFACES: usize = 2_048;
/// Maximum durable topology/state elements owned across all subsystems.
pub const MAX_MACHINE_GRAPH_OWNED_ELEMENTS: usize = 16_384;

const MACHINE_IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(4_096, 128, 1, 1, 256);
const LINEAGE_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(4 * 1_024 * 1_024, 1_024 * 1_024, 5, 16_384, 4_096);
const MACHINE_GRAPH_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(16 * 1_024 * 1_024, 4 * 1_024 * 1_024, 8, 64_000, 16_384);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);

/// Nominal role of a durable Machine-IR entity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum MachineElementKind {
    /// Physical or logical body occurrence.
    Body = 1,
    /// Named surface support that may survive remeshing.
    SurfacePatch = 2,
    /// Contact attachment feature.
    ContactFeature = 3,
    /// Typed subsystem terminal.
    Terminal = 4,
    /// Coupling port.
    Port = 5,
    /// Owned state slot.
    StateSlot = 6,
}

impl MachineElementKind {
    /// Stable binary tag used by lineage schemas.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Body => "body",
            Self::SurfacePatch => "surface-patch",
            Self::ContactFeature => "contact-feature",
            Self::Terminal => "terminal",
            Self::Port => "port",
            Self::StateSlot => "state-slot",
        }
    }
}

/// Structured refusal from canonical Machine-IR key construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MachineIdError {
    /// The identifier was empty.
    Empty {
        /// Nominal identifier role.
        role: &'static str,
    },
    /// The identifier exceeded its public byte envelope.
    TooLong {
        /// Nominal identifier role.
        role: &'static str,
        /// Supplied UTF-8 byte count.
        bytes: usize,
        /// Maximum admitted byte count.
        max: usize,
    },
    /// A slash-delimited key segment was empty.
    EmptySegment {
        /// Nominal identifier role.
        role: &'static str,
        /// Zero-based segment index.
        segment: usize,
    },
    /// A segment did not begin with an ASCII lowercase letter.
    InvalidSegmentStart {
        /// Nominal identifier role.
        role: &'static str,
        /// Zero-based segment index.
        segment: usize,
        /// Byte offset in the complete key.
        at: usize,
        /// Offending byte.
        byte: u8,
    },
    /// A byte was outside the canonical `[a-z][a-z0-9-]*` segment grammar.
    InvalidByte {
        /// Nominal identifier role.
        role: &'static str,
        /// Byte offset in the complete key.
        at: usize,
        /// Offending byte.
        byte: u8,
    },
    /// A segment ended in `-`.
    TrailingSeparator {
        /// Nominal identifier role.
        role: &'static str,
        /// Zero-based segment index.
        segment: usize,
    },
    /// Adjacent `-` separators would create a noncanonical spelling.
    RepeatedSeparator {
        /// Nominal identifier role.
        role: &'static str,
        /// Byte offset of the second separator.
        at: usize,
    },
    /// The bounded canonical identity encoder refused publication.
    Identity(CanonicalError),
}

impl MachineIdError {
    /// Stable rule code for structured admission diagnostics.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Empty { .. } => "MachineIdEmpty",
            Self::TooLong { .. } => "MachineIdTooLong",
            Self::EmptySegment { .. } => "MachineIdEmptySegment",
            Self::InvalidSegmentStart { .. } => "MachineIdInvalidSegmentStart",
            Self::InvalidByte { .. } => "MachineIdInvalidByte",
            Self::TrailingSeparator { .. } => "MachineIdTrailingSeparator",
            Self::RepeatedSeparator { .. } => "MachineIdRepeatedSeparator",
            Self::Identity(_) => "MachineIdIdentity",
        }
    }
}

impl fmt::Display for MachineIdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty { role } => write!(f, "{role} id must not be empty"),
            Self::TooLong { role, bytes, max } => {
                write!(f, "{role} id has {bytes} bytes; maximum is {max}")
            }
            Self::EmptySegment { role, segment } => {
                write!(f, "{role} id segment {segment} must not be empty")
            }
            Self::InvalidSegmentStart {
                role,
                segment,
                at,
                byte,
            } => write!(
                f,
                "{role} id segment {segment} must start with an ASCII lowercase letter; byte \
                 0x{byte:02x} at offset {at} is invalid"
            ),
            Self::InvalidByte { role, at, byte } => write!(
                f,
                "{role} id byte 0x{byte:02x} at offset {at} is outside the canonical \
                 [a-z][a-z0-9-]* segment grammar"
            ),
            Self::TrailingSeparator { role, segment } => {
                write!(f, "{role} id segment {segment} must not end in '-'")
            }
            Self::RepeatedSeparator { role, at } => write!(
                f,
                "{role} id contains repeated '-' separators ending at offset {at}"
            ),
            Self::Identity(error) => write!(f, "Machine-IR identity refused: {error}"),
        }
    }
}

impl core::error::Error for MachineIdError {}

impl From<CanonicalError> for MachineIdError {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

fn validate_canonical_key(role: &'static str, value: &str) -> Result<(), MachineIdError> {
    if value.is_empty() {
        return Err(MachineIdError::Empty { role });
    }
    if value.len() > MAX_MACHINE_ENTITY_KEY_BYTES {
        return Err(MachineIdError::TooLong {
            role,
            bytes: value.len(),
            max: MAX_MACHINE_ENTITY_KEY_BYTES,
        });
    }

    let mut base = 0usize;
    for (segment_index, segment) in value.as_bytes().split(|byte| *byte == b'/').enumerate() {
        if segment.is_empty() {
            return Err(MachineIdError::EmptySegment {
                role,
                segment: segment_index,
            });
        }
        if !segment[0].is_ascii_lowercase() {
            return Err(MachineIdError::InvalidSegmentStart {
                role,
                segment: segment_index,
                at: base,
                byte: segment[0],
            });
        }
        let mut previous_separator = false;
        for (offset, byte) in segment.iter().copied().enumerate() {
            let at = base + offset;
            match byte {
                b'a'..=b'z' | b'0'..=b'9' => previous_separator = false,
                b'-' if previous_separator => {
                    return Err(MachineIdError::RepeatedSeparator { role, at });
                }
                b'-' => previous_separator = true,
                _ => return Err(MachineIdError::InvalidByte { role, at, byte }),
            }
        }
        if previous_separator {
            return Err(MachineIdError::TrailingSeparator {
                role,
                segment: segment_index,
            });
        }
        base += segment.len() + 1;
    }
    Ok(())
}

macro_rules! durable_machine_id {
    (
        $(#[$meta:meta])*
        $name:ident,
        $schema:ident,
        $identity:ident,
        $kind:expr,
        $role:literal,
        $domain:literal,
        $context:literal
    ) => {
        #[doc = concat!("Canonical schema marker for `", stringify!($name), "`.")]
        pub enum $schema {}

        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str = $domain;
            const NAME: &'static str = $role;
            const VERSION: u32 = MACHINE_ENTITY_ID_SCHEMA_VERSION_V1;
            const CONTEXT: &'static str = $context;
            const FIELDS: &'static [FieldSpec] =
                &[FieldSpec::required("canonical-key", WireType::Utf8)];
        }

        #[doc = concat!("Typed durable entity digest for `", stringify!($name), "`.")]
        pub type $identity = EntityId<$schema>;

        $(#[$meta])*
        #[derive(Clone)]
        pub struct $name {
            canonical_key: Box<str>,
            receipt: IdentityReceipt<$identity>,
        }

        impl $name {
            /// Admit a canonical, human-auditable durable key.
            ///
            /// # Errors
            /// Refuses noncanonical text or a bounded identity-encoding error.
            pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
                let key = key.into();
                validate_canonical_key($role, &key)?;
                let receipt = CanonicalEncoder::<$identity, _>::new(
                    MACHINE_IDENTITY_LIMITS,
                    NeverCancel,
                )?
                .utf8(Field::new(0, "canonical-key"), &key)?
                .finish()?;
                Ok(Self {
                    canonical_key: key.into_boxed_str(),
                    receipt,
                })
            }

            /// Nominal entity kind.
            #[must_use]
            pub const fn kind(&self) -> MachineElementKind {
                $kind
            }

            /// Exact canonical key retained for diagnostics and wire lowering.
            #[must_use]
            pub fn canonical_key(&self) -> &str {
                &self.canonical_key
            }

            /// Domain- and role-separated durable entity identity.
            #[must_use]
            pub const fn identity(&self) -> $identity {
                self.receipt.id()
            }

            /// Complete canonical-preimage receipt for collision adjudication.
            #[must_use]
            pub const fn identity_receipt(&self) -> IdentityReceipt<$identity> {
                self.receipt
            }

            fn digest_bytes(&self) -> [u8; 32] {
                *self.identity().as_bytes()
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &Self) -> bool {
                self.identity() == other.identity()
            }
        }

        impl Eq for $name {}

        impl PartialOrd for $name {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Ord for $name {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                self.identity().cmp(&other.identity())
            }
        }

        impl Hash for $name {
            fn hash<H: Hasher>(&self, state: &mut H) {
                self.identity().hash(state);
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.canonical_key)
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("canonical_key", &self.canonical_key)
                    .field("identity", &self.identity())
                    .finish()
            }
        }
    };
}

durable_machine_id!(
    /// Durable identity of a Machine-IR body occurrence.
    BodyId,
    BodyIdSchemaV1,
    BodyEntityIdV1,
    MachineElementKind::Body,
    "body-id",
    "org.frankensim.fs-ir.machine.body-id.v1",
    "one durable Machine-IR body occurrence named by a canonical hierarchical key"
);
durable_machine_id!(
    /// Durable identity of a surface support independent of mesh indices.
    SurfacePatchId,
    SurfacePatchIdSchemaV1,
    SurfacePatchEntityIdV1,
    MachineElementKind::SurfacePatch,
    "surface-patch-id",
    "org.frankensim.fs-ir.machine.surface-patch-id.v1",
    "one durable Machine-IR surface patch independent of discretization indices"
);
durable_machine_id!(
    /// Durable identity of a contact attachment feature.
    ContactFeatureId,
    ContactFeatureIdSchemaV1,
    ContactFeatureEntityIdV1,
    MachineElementKind::ContactFeature,
    "contact-feature-id",
    "org.frankensim.fs-ir.machine.contact-feature-id.v1",
    "one durable Machine-IR contact feature independent of solver ordering"
);
durable_machine_id!(
    /// Durable identity of a typed subsystem terminal.
    TerminalId,
    TerminalIdSchemaV1,
    TerminalEntityIdV1,
    MachineElementKind::Terminal,
    "terminal-id",
    "org.frankensim.fs-ir.machine.terminal-id.v1",
    "one durable Machine-IR terminal independent of graph serialization order"
);
durable_machine_id!(
    /// Durable identity of a coupling port.
    PortId,
    PortIdSchemaV1,
    PortEntityIdV1,
    MachineElementKind::Port,
    "port-id",
    "org.frankensim.fs-ir.machine.port-id.v1",
    "one durable Machine-IR coupling port independent of graph serialization order"
);
durable_machine_id!(
    /// Durable identity of an owned state slot.
    StateSlotId,
    StateSlotIdSchemaV1,
    StateSlotEntityIdV1,
    MachineElementKind::StateSlot,
    "state-slot-id",
    "org.frankensim.fs-ir.machine.state-slot-id.v1",
    "one durable Machine-IR state slot independent of vector position"
);

/// One of the six non-confusable durable Machine-IR entity roles.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineElementId {
    /// Body occurrence.
    Body(BodyId),
    /// Surface support.
    SurfacePatch(SurfacePatchId),
    /// Contact feature.
    ContactFeature(ContactFeatureId),
    /// Terminal.
    Terminal(TerminalId),
    /// Port.
    Port(PortId),
    /// State slot.
    StateSlot(StateSlotId),
}

impl MachineElementId {
    /// Nominal role.
    #[must_use]
    pub const fn kind(&self) -> MachineElementKind {
        match self {
            Self::Body(_) => MachineElementKind::Body,
            Self::SurfacePatch(_) => MachineElementKind::SurfacePatch,
            Self::ContactFeature(_) => MachineElementKind::ContactFeature,
            Self::Terminal(_) => MachineElementKind::Terminal,
            Self::Port(_) => MachineElementKind::Port,
            Self::StateSlot(_) => MachineElementKind::StateSlot,
        }
    }

    /// Human-auditable canonical key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        match self {
            Self::Body(id) => id.canonical_key(),
            Self::SurfacePatch(id) => id.canonical_key(),
            Self::ContactFeature(id) => id.canonical_key(),
            Self::Terminal(id) => id.canonical_key(),
            Self::Port(id) => id.canonical_key(),
            Self::StateSlot(id) => id.canonical_key(),
        }
    }

    fn digest_bytes(&self) -> [u8; 32] {
        match self {
            Self::Body(id) => id.digest_bytes(),
            Self::SurfacePatch(id) => id.digest_bytes(),
            Self::ContactFeature(id) => id.digest_bytes(),
            Self::Terminal(id) => id.digest_bytes(),
            Self::Port(id) => id.digest_bytes(),
            Self::StateSlot(id) => id.digest_bytes(),
        }
    }

    fn canonical_row(&self) -> [u8; 33] {
        let mut row = [0u8; 33];
        row[0] = self.kind().tag();
        row[1..].copy_from_slice(&self.digest_bytes());
        row
    }
}

impl fmt::Display for MachineElementId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.kind().name(), self.canonical_key())
    }
}

macro_rules! element_conversion {
    ($id:ident, $variant:ident) => {
        impl From<$id> for MachineElementId {
            fn from(value: $id) -> Self {
                Self::$variant(value)
            }
        }
    };
}

element_conversion!(BodyId, Body);
element_conversion!(SurfacePatchId, SurfacePatch);
element_conversion!(ContactFeatureId, ContactFeature);
element_conversion!(TerminalId, Terminal);
element_conversion!(PortId, Port);
element_conversion!(StateSlotId, StateSlot);

/// Topology or discretization event represented by a lineage morphism.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum LineageEvent {
    /// One entity becomes multiple descendants.
    Split = 1,
    /// Multiple entities become one descendant.
    Merge = 2,
    /// Discretization changes while physical support may be preserved or split.
    Remesh = 3,
    /// Wear advances an entity to a successor support.
    Wear = 4,
    /// Fracture creates multiple descendants.
    Fracture = 5,
}

impl LineageEvent {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u64 {
        self as u64
    }
}

/// One source entity and its explicitly declared target set.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LineageRelation {
    source: MachineElementId,
    targets: Vec<MachineElementId>,
}

impl LineageRelation {
    /// Construct and canonicalize one source relation.
    ///
    /// # Errors
    /// Refuses an empty/oversized target set, cross-role targets, or duplicates.
    pub fn new(
        source: MachineElementId,
        mut targets: Vec<MachineElementId>,
    ) -> Result<Self, LineageRefusal> {
        if targets.is_empty() {
            return Err(LineageRefusal::NoTargets { source });
        }
        if targets.len() > MAX_LINEAGE_TARGETS_PER_SOURCE {
            return Err(LineageRefusal::TargetLimit {
                source,
                count: targets.len(),
                max: MAX_LINEAGE_TARGETS_PER_SOURCE,
            });
        }
        targets.sort();
        if let Some(target) = targets.iter().find(|target| target.kind() != source.kind()) {
            return Err(LineageRefusal::TargetKindMismatch {
                source,
                target: target.clone(),
            });
        }
        if let Some(pair) = targets.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(LineageRefusal::DuplicateTarget {
                source,
                target: pair[0].clone(),
            });
        }
        Ok(Self { source, targets })
    }

    /// Source entity.
    #[must_use]
    pub const fn source(&self) -> &MachineElementId {
        &self.source
    }

    /// Canonically ordered target set.
    #[must_use]
    pub fn targets(&self) -> &[MachineElementId] {
        &self.targets
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(33 + 8 + self.targets.len() * 33);
        row.extend_from_slice(&self.source.canonical_row());
        row.extend_from_slice(&(self.targets.len() as u64).to_le_bytes());
        for target in &self.targets {
            row.extend_from_slice(&target.canonical_row());
        }
        row
    }
}

/// Downstream artifact class whose attachment must follow or invalidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum DependentKind {
    /// Derived cache keyed by the source entity.
    Cache = 1,
    /// Contact relation attached to the source entity.
    Contact = 2,
    /// Winding or coil attachment.
    Winding = 3,
    /// Adjoint/tape dependency.
    Adjoint = 4,
}

impl DependentKind {
    /// Stable identity tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Cache => "cache",
            Self::Contact => "contact",
            Self::Winding => "winding",
            Self::Adjoint => "adjoint",
        }
    }
}

/// One named downstream attachment to a durable source entity.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DependentBinding {
    kind: DependentKind,
    canonical_key: Box<str>,
    source: MachineElementId,
}

impl DependentBinding {
    /// Construct a checked dependent attachment.
    ///
    /// # Errors
    /// Refuses a noncanonical dependent key.
    pub fn new(
        kind: DependentKind,
        key: impl Into<String>,
        source: MachineElementId,
    ) -> Result<Self, MachineIdError> {
        let key = key.into();
        validate_canonical_key(kind.name(), &key)?;
        Ok(Self {
            kind,
            canonical_key: key.into_boxed_str(),
            source,
        })
    }

    /// Downstream artifact class.
    #[must_use]
    pub const fn kind(&self) -> DependentKind {
        self.kind
    }

    /// Canonical dependent identity key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.canonical_key
    }

    /// Entity to which this dependent is attached.
    #[must_use]
    pub const fn source(&self) -> &MachineElementId {
        &self.source
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = Vec::with_capacity(1 + 8 + self.canonical_key.len() + 33);
        row.push(self.kind.tag());
        row.extend_from_slice(&(self.canonical_key.len() as u64).to_le_bytes());
        row.extend_from_slice(self.canonical_key.as_bytes());
        row.extend_from_slice(&self.source.canonical_row());
        row
    }
}

/// One admitted deterministic attachment move.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LineageRebinding {
    dependent: DependentBinding,
    target: MachineElementId,
}

impl LineageRebinding {
    /// Dependent whose attachment moves.
    #[must_use]
    pub const fn dependent(&self) -> &DependentBinding {
        &self.dependent
    }

    /// Unique target selected by the admitted relation.
    #[must_use]
    pub const fn target(&self) -> &MachineElementId {
        &self.target
    }

    fn canonical_row(&self) -> Vec<u8> {
        let mut row = self.dependent.canonical_row();
        row.extend_from_slice(&self.target.canonical_row());
        row
    }
}

/// Canonical schema marker for one admitted lineage record.
pub enum LineageRecordSchemaV1 {}

impl CanonicalSchema for LineageRecordSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.lineage-record.v1";
    const NAME: &'static str = "machine-lineage-record";
    const VERSION: u32 = MACHINE_LINEAGE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one canonical topology/discretization event, its explicit source-target relations, and only deterministic dependent rebindings";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("event", WireType::U64),
        FieldSpec::required("relations", WireType::OrderedBytes),
        FieldSpec::required("rebindings", WireType::OrderedBytes),
    ];
}

/// Semantic identity of one admitted lineage record.
pub type LineageRecordIdV1 = SemanticId<LineageRecordSchemaV1>;

/// Canonical schema marker for a fail-closed lineage invalidation.
pub enum LineageInvalidationSchemaV1 {}

impl CanonicalSchema for LineageInvalidationSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.lineage-invalidation.v1";
    const NAME: &'static str = "machine-lineage-invalidation";
    const VERSION: u32 = MACHINE_LINEAGE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "one refused ambiguous lineage event, its complete source-target relation set, every considered dependent binding, the one-to-many subset that caused ambiguity, and every dependent invalidated instead of silently rebound";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("event", WireType::U64),
        FieldSpec::required("relations", WireType::OrderedBytes),
        FieldSpec::required("considered-dependents", WireType::OrderedBytes),
        FieldSpec::required("ambiguous-relations", WireType::OrderedBytes),
        FieldSpec::required("invalidated-dependents", WireType::OrderedBytes),
    ];
}

/// Semantic identity of one fail-closed invalidation set.
pub type LineageInvalidationIdV1 = SemanticId<LineageInvalidationSchemaV1>;

/// Published lineage event containing only unambiguous rebindings.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LineageRecord {
    event: LineageEvent,
    relations: Vec<LineageRelation>,
    rebindings: Vec<LineageRebinding>,
    receipt: IdentityReceipt<LineageRecordIdV1>,
}

impl LineageRecord {
    /// Admit a bounded lineage event and its dependent attachments.
    ///
    /// Inputs are canonicalized by durable identity, so caller order is not
    /// semantic. A dependent on a one-to-many source never reaches the success
    /// type: admission returns [`LineageRefusal::Ambiguous`] with a complete
    /// [`LineageInvalidation`] receipt.
    ///
    /// # Errors
    /// Refuses malformed event shapes, duplicates, missing sources, resource
    /// limits, identity publication failure, or ambiguous live attachments.
    pub fn admit(
        event: LineageEvent,
        relations: Vec<LineageRelation>,
        dependents: Vec<DependentBinding>,
    ) -> Result<Self, LineageRefusal> {
        Self::admit_with_decision(event, relations, dependents).into_result()
    }

    /// Execute admission and retain its deterministic outcome summary.
    ///
    /// The summary is suitable for structured tracing: it carries the attempted
    /// event and input counts together with either the admitted lineage receipt
    /// or the exact typed refusal. A replay ledger must separately retain the
    /// canonical attempted inputs; early refusals intentionally discard them.
    #[must_use]
    pub fn admit_with_decision(
        event: LineageEvent,
        relations: Vec<LineageRelation>,
        dependents: Vec<DependentBinding>,
    ) -> LineageAdmissionDecision {
        let submitted_relations = relations.len();
        let submitted_dependents = dependents.len();
        let result = Self::admit_inner(event, relations, dependents);
        LineageAdmissionDecision {
            event,
            submitted_relations,
            submitted_dependents,
            result,
        }
    }

    fn admit_inner(
        event: LineageEvent,
        mut relations: Vec<LineageRelation>,
        mut dependents: Vec<DependentBinding>,
    ) -> Result<Self, LineageRefusal> {
        if relations.is_empty() {
            return Err(LineageRefusal::NoRelations);
        }
        if relations.len() > MAX_LINEAGE_RELATIONS {
            return Err(LineageRefusal::RelationLimit {
                count: relations.len(),
                max: MAX_LINEAGE_RELATIONS,
            });
        }
        if dependents.len() > MAX_LINEAGE_DEPENDENTS {
            return Err(LineageRefusal::DependentLimit {
                count: dependents.len(),
                max: MAX_LINEAGE_DEPENDENTS,
            });
        }
        let endpoints = relations.iter().try_fold(0usize, |total, relation| {
            total.checked_add(relation.targets.len())
        });
        let Some(endpoints) = endpoints else {
            return Err(LineageRefusal::EndpointLimit {
                count: usize::MAX,
                max: MAX_LINEAGE_ENDPOINTS,
            });
        };
        if endpoints > MAX_LINEAGE_ENDPOINTS {
            return Err(LineageRefusal::EndpointLimit {
                count: endpoints,
                max: MAX_LINEAGE_ENDPOINTS,
            });
        }

        relations.sort_by(|left, right| left.source.cmp(&right.source));
        if let Some(pair) = relations
            .windows(2)
            .find(|pair| pair[0].source == pair[1].source)
        {
            return Err(LineageRefusal::DuplicateSource {
                source: pair[0].source.clone(),
            });
        }
        validate_event_shape(event, &relations)?;

        dependents.sort();
        if let Some(pair) = dependents.windows(2).find(|pair| {
            pair[0].kind == pair[1].kind && pair[0].canonical_key == pair[1].canonical_key
        }) {
            return Err(LineageRefusal::DuplicateDependent {
                kind: pair[0].kind,
                key: pair[0].canonical_key.clone(),
            });
        }
        for dependent in &dependents {
            if relations
                .binary_search_by(|relation| relation.source.cmp(&dependent.source))
                .is_err()
            {
                return Err(LineageRefusal::UnknownDependentSource {
                    kind: dependent.kind,
                    key: dependent.canonical_key.clone(),
                    source: dependent.source.clone(),
                });
            }
        }

        let ambiguous_relations: Vec<_> = relations
            .iter()
            .filter(|relation| relation.targets.len() != 1)
            .cloned()
            .collect();
        let invalidated_dependents: Vec<_> = dependents
            .iter()
            .filter(|dependent| {
                ambiguous_relations
                    .binary_search_by(|relation| relation.source.cmp(&dependent.source))
                    .is_ok()
            })
            .cloned()
            .collect();
        if !invalidated_dependents.is_empty() {
            return Err(LineageRefusal::Ambiguous(LineageInvalidation::new(
                event,
                relations,
                dependents,
                ambiguous_relations,
                invalidated_dependents,
            )?));
        }

        let mut rebindings = Vec::with_capacity(dependents.len());
        for dependent in dependents {
            let relation_index = relations
                .binary_search_by(|relation| relation.source.cmp(&dependent.source))
                .expect("dependent source membership checked above");
            let relation = &relations[relation_index];
            if let [target] = relation.targets.as_slice() {
                rebindings.push(LineageRebinding {
                    dependent,
                    target: target.clone(),
                });
            }
        }
        rebindings.sort();
        let receipt = lineage_record_identity(event, &relations, &rebindings)?;
        Ok(Self {
            event,
            relations,
            rebindings,
            receipt,
        })
    }

    /// Event kind.
    #[must_use]
    pub const fn event(&self) -> LineageEvent {
        self.event
    }

    /// Canonically ordered complete source-target relations.
    #[must_use]
    pub fn relations(&self) -> &[LineageRelation] {
        &self.relations
    }

    /// Canonically ordered deterministic attachment moves.
    #[must_use]
    pub fn rebindings(&self) -> &[LineageRebinding] {
        &self.rebindings
    }

    /// Domain-separated semantic identity.
    #[must_use]
    pub const fn identity(&self) -> LineageRecordIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<LineageRecordIdV1> {
        self.receipt
    }
}

/// Exact fail-closed payload for ambiguous live attachments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LineageInvalidation {
    event: LineageEvent,
    relations: Vec<LineageRelation>,
    considered_dependents: Vec<DependentBinding>,
    ambiguous_relations: Vec<LineageRelation>,
    invalidated_dependents: Vec<DependentBinding>,
    receipt: IdentityReceipt<LineageInvalidationIdV1>,
}

impl LineageInvalidation {
    fn new(
        event: LineageEvent,
        relations: Vec<LineageRelation>,
        considered_dependents: Vec<DependentBinding>,
        ambiguous_relations: Vec<LineageRelation>,
        invalidated_dependents: Vec<DependentBinding>,
    ) -> Result<Self, CanonicalError> {
        let receipt = lineage_invalidation_identity(
            event,
            &relations,
            &considered_dependents,
            &ambiguous_relations,
            &invalidated_dependents,
        )?;
        Ok(Self {
            event,
            relations,
            considered_dependents,
            ambiguous_relations,
            invalidated_dependents,
            receipt,
        })
    }

    /// Refused event kind.
    #[must_use]
    pub const fn event(&self) -> LineageEvent {
        self.event
    }

    /// Complete canonically ordered relation set of the refused event.
    #[must_use]
    pub fn relations(&self) -> &[LineageRelation] {
        &self.relations
    }

    /// Every canonically ordered dependent considered by the refused event.
    #[must_use]
    pub fn considered_dependents(&self) -> &[DependentBinding] {
        &self.considered_dependents
    }

    /// Every one-to-many relation in the refused event.
    #[must_use]
    pub fn ambiguous_relations(&self) -> &[LineageRelation] {
        &self.ambiguous_relations
    }

    /// Complete deterministic invalidation set; nothing was rebound.
    #[must_use]
    pub fn invalidated_dependents(&self) -> &[DependentBinding] {
        &self.invalidated_dependents
    }

    /// Domain-separated invalidation identity.
    #[must_use]
    pub const fn identity(&self) -> LineageInvalidationIdV1 {
        self.receipt.id()
    }

    /// Complete canonical-preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<LineageInvalidationIdV1> {
        self.receipt
    }
}

/// Structured fail-closed lineage refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LineageRefusal {
    /// No source relation was supplied.
    NoRelations,
    /// Too many relations were supplied.
    RelationLimit {
        /// Supplied relation count.
        count: usize,
        /// Public maximum.
        max: usize,
    },
    /// One relation had no target.
    NoTargets {
        /// Offending source.
        source: MachineElementId,
    },
    /// One source had too many targets.
    TargetLimit {
        /// Offending source.
        source: MachineElementId,
        /// Supplied target count.
        count: usize,
        /// Public maximum.
        max: usize,
    },
    /// The whole event had too many target endpoints.
    EndpointLimit {
        /// Supplied or saturated count.
        count: usize,
        /// Public maximum.
        max: usize,
    },
    /// Source and target roles differed.
    TargetKindMismatch {
        /// Source entity.
        source: MachineElementId,
        /// Incompatible target entity.
        target: MachineElementId,
    },
    /// One source listed the same target twice.
    DuplicateTarget {
        /// Source entity.
        source: MachineElementId,
        /// Duplicate target.
        target: MachineElementId,
    },
    /// The event declared the same source more than once.
    DuplicateSource {
        /// Duplicate source.
        source: MachineElementId,
    },
    /// The relation cardinality did not match the named event.
    EventShape {
        /// Event kind.
        event: LineageEvent,
        /// Stable violated rule.
        rule: &'static str,
    },
    /// Too many dependent attachments were supplied.
    DependentLimit {
        /// Supplied dependent count.
        count: usize,
        /// Public maximum.
        max: usize,
    },
    /// A dependent named no source relation in this event.
    UnknownDependentSource {
        /// Dependent class.
        kind: DependentKind,
        /// Dependent key.
        key: Box<str>,
        /// Missing source.
        source: MachineElementId,
    },
    /// A dependent identity appeared more than once.
    DuplicateDependent {
        /// Dependent class.
        kind: DependentKind,
        /// Duplicate key.
        key: Box<str>,
    },
    /// Canonical identity publication failed.
    Identity(CanonicalError),
    /// One or more live dependents had multiple possible targets. The payload
    /// is the exact invalidation record; no success value was published.
    Ambiguous(LineageInvalidation),
}

impl LineageRefusal {
    /// Stable rule code for structured admission diagnostics.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoRelations => "LineageNoRelations",
            Self::RelationLimit { .. } => "LineageRelationLimit",
            Self::NoTargets { .. } => "LineageNoTargets",
            Self::TargetLimit { .. } => "LineageTargetLimit",
            Self::EndpointLimit { .. } => "LineageEndpointLimit",
            Self::TargetKindMismatch { .. } => "LineageTargetKindMismatch",
            Self::DuplicateTarget { .. } => "LineageDuplicateTarget",
            Self::DuplicateSource { .. } => "LineageDuplicateSource",
            Self::EventShape { .. } => "LineageEventShape",
            Self::DependentLimit { .. } => "LineageDependentLimit",
            Self::UnknownDependentSource { .. } => "LineageUnknownDependentSource",
            Self::DuplicateDependent { .. } => "LineageDuplicateDependent",
            Self::Identity(_) => "LineageIdentity",
            Self::Ambiguous(_) => "LineageAmbiguous",
        }
    }

    /// Exact invalidation payload when this is an ambiguity refusal.
    #[must_use]
    pub const fn invalidation(&self) -> Option<&LineageInvalidation> {
        match self {
            Self::Ambiguous(invalidation) => Some(invalidation),
            _ => None,
        }
    }
}

/// Bounded deterministic outcome summary for one lineage-admission attempt.
///
/// The core library does not print or select a tracing backend. Callers can
/// inspect this summary without reconstructing a refusal from display text.
/// It is not a digest or replay record of the complete submitted graph.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LineageAdmissionDecision {
    event: LineageEvent,
    submitted_relations: usize,
    submitted_dependents: usize,
    result: Result<LineageRecord, LineageRefusal>,
}

impl LineageAdmissionDecision {
    /// Attempted lineage-event kind, including for early refusals.
    #[must_use]
    pub const fn event(&self) -> LineageEvent {
        self.event
    }

    /// Number of relations submitted before canonicalization or refusal.
    #[must_use]
    pub const fn submitted_relation_count(&self) -> usize {
        self.submitted_relations
    }

    /// Number of dependent bindings submitted before canonicalization.
    #[must_use]
    pub const fn submitted_dependent_count(&self) -> usize {
        self.submitted_dependents
    }

    /// Stable decision/rule code for structured event fields.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "LineageAdmitted",
            Err(refusal) => refusal.code(),
        }
    }

    /// Borrow the admitted record, or the complete refusal.
    #[must_use]
    pub fn result(&self) -> Result<&LineageRecord, &LineageRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision record and recover the conventional result.
    #[must_use]
    pub fn into_result(self) -> Result<LineageRecord, LineageRefusal> {
        self.result
    }
}

impl From<CanonicalError> for LineageRefusal {
    fn from(error: CanonicalError) -> Self {
        Self::Identity(error)
    }
}

impl fmt::Display for LineageRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoRelations => f.write_str("lineage event must contain at least one relation"),
            Self::RelationLimit { count, max } => {
                write!(f, "lineage event has {count} relations; maximum is {max}")
            }
            Self::NoTargets { source } => write!(f, "lineage source {source} has no target"),
            Self::TargetLimit { source, count, max } => write!(
                f,
                "lineage source {source} has {count} targets; maximum is {max}"
            ),
            Self::EndpointLimit { count, max } => {
                write!(f, "lineage event has {count} endpoints; maximum is {max}")
            }
            Self::TargetKindMismatch { source, target } => write!(
                f,
                "lineage target {target} is not the same nominal entity kind as source {source}"
            ),
            Self::DuplicateTarget { source, target } => {
                write!(f, "lineage source {source} repeats target {target}")
            }
            Self::DuplicateSource { source } => {
                write!(f, "lineage event repeats source {source}")
            }
            Self::EventShape { event, rule } => {
                write!(f, "{event:?} lineage shape refused: {rule}")
            }
            Self::DependentLimit { count, max } => write!(
                f,
                "lineage event has {count} dependent attachments; maximum is {max}"
            ),
            Self::UnknownDependentSource { kind, key, source } => write!(
                f,
                "{} dependent `{key}` names source {source}, which is absent from the event",
                kind.name()
            ),
            Self::DuplicateDependent { kind, key } => {
                write!(
                    f,
                    "{} dependent `{key}` appears more than once",
                    kind.name()
                )
            }
            Self::Identity(error) => write!(f, "lineage identity refused: {error}"),
            Self::Ambiguous(invalidation) => write!(
                f,
                "ambiguous {:?} lineage invalidated {} dependent attachment(s); no automatic \
                 rebinding was published",
                invalidation.event,
                invalidation.invalidated_dependents.len()
            ),
        }
    }
}

impl core::error::Error for LineageRefusal {}

fn validate_event_shape(
    event: LineageEvent,
    relations: &[LineageRelation],
) -> Result<(), LineageRefusal> {
    match event {
        LineageEvent::Split | LineageEvent::Fracture => {
            if relations.len() != 1 || relations[0].targets.len() < 2 {
                return Err(LineageRefusal::EventShape {
                    event,
                    rule: "split/fracture requires exactly one source and at least two targets",
                });
            }
        }
        LineageEvent::Merge => {
            if relations.len() < 2 || relations.iter().any(|relation| relation.targets.len() != 1) {
                return Err(LineageRefusal::EventShape {
                    event,
                    rule: "merge requires at least two sources, each with exactly one target",
                });
            }
            let target = &relations[0].targets[0];
            if relations
                .iter()
                .skip(1)
                .any(|relation| relation.targets[0] != *target)
            {
                return Err(LineageRefusal::EventShape {
                    event,
                    rule: "every merge source must name the same target",
                });
            }
        }
        LineageEvent::Wear => {
            if relations.iter().any(|relation| relation.targets.len() != 1) {
                return Err(LineageRefusal::EventShape {
                    event,
                    rule: "wear requires exactly one successor for every source",
                });
            }
        }
        LineageEvent::Remesh => {}
    }
    Ok(())
}

fn lineage_record_identity(
    event: LineageEvent,
    relations: &[LineageRelation],
    rebindings: &[LineageRebinding],
) -> Result<IdentityReceipt<LineageRecordIdV1>, CanonicalError> {
    let relation_rows: Vec<_> = relations
        .iter()
        .map(LineageRelation::canonical_row)
        .collect();
    let rebinding_rows: Vec<_> = rebindings
        .iter()
        .map(LineageRebinding::canonical_row)
        .collect();
    CanonicalEncoder::<LineageRecordIdV1, _>::new(LINEAGE_IDENTITY_LIMITS, NeverCancel)?
        .u64(Field::new(0, "event"), event.tag())?
        .ordered_bytes(
            Field::new(1, "relations"),
            relation_rows.len() as u64,
            relation_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(2, "rebindings"),
            rebinding_rows.len() as u64,
            rebinding_rows.iter().map(Vec::as_slice),
        )?
        .finish()
}

fn lineage_invalidation_identity(
    event: LineageEvent,
    relations: &[LineageRelation],
    considered_dependents: &[DependentBinding],
    ambiguous_relations: &[LineageRelation],
    dependents: &[DependentBinding],
) -> Result<IdentityReceipt<LineageInvalidationIdV1>, CanonicalError> {
    let relation_rows: Vec<_> = relations
        .iter()
        .map(LineageRelation::canonical_row)
        .collect();
    let considered_rows: Vec<_> = considered_dependents
        .iter()
        .map(DependentBinding::canonical_row)
        .collect();
    let ambiguous_rows: Vec<_> = ambiguous_relations
        .iter()
        .map(LineageRelation::canonical_row)
        .collect();
    let dependent_rows: Vec<_> = dependents
        .iter()
        .map(DependentBinding::canonical_row)
        .collect();
    CanonicalEncoder::<LineageInvalidationIdV1, _>::new(LINEAGE_IDENTITY_LIMITS, NeverCancel)?
        .u64(Field::new(0, "event"), event.tag())?
        .ordered_bytes(
            Field::new(1, "relations"),
            relation_rows.len() as u64,
            relation_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(2, "considered-dependents"),
            considered_rows.len() as u64,
            considered_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(3, "ambiguous-relations"),
            ambiguous_rows.len() as u64,
            ambiguous_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(4, "invalidated-dependents"),
            dependent_rows.len() as u64,
            dependent_rows.iter().map(Vec::as_slice),
        )?
        .finish()
}

// ---------------------------------------------------------------------------
// Machine graph schema and admission (E0 PR-2).
// ---------------------------------------------------------------------------

macro_rules! durable_graph_id {
    (
        $(#[$meta:meta])*
        $name:ident,
        $schema:ident,
        $identity:ident,
        $role:literal,
        $domain:literal,
        $context:literal
    ) => {
        #[doc = concat!("Canonical schema marker for `", stringify!($name), "`.")]
        pub enum $schema {}

        impl CanonicalSchema for $schema {
            const DOMAIN: &'static str = $domain;
            const NAME: &'static str = $role;
            const VERSION: u32 = MACHINE_GRAPH_SCHEMA_VERSION_V1;
            const CONTEXT: &'static str = $context;
            const FIELDS: &'static [FieldSpec] =
                &[FieldSpec::required("canonical-key", WireType::Utf8)];
        }

        #[doc = concat!("Typed durable graph identity for `", stringify!($name), "`.")]
        pub type $identity = EntityId<$schema>;

        $(#[$meta])*
        #[derive(Clone)]
        pub struct $name {
            canonical_key: Box<str>,
            receipt: IdentityReceipt<$identity>,
        }

        impl $name {
            /// Admit a canonical, human-auditable graph key.
            ///
            /// # Errors
            /// Refuses noncanonical text or a bounded identity-encoding error.
            pub fn new(key: impl Into<String>) -> Result<Self, MachineIdError> {
                let key = key.into();
                validate_canonical_key($role, &key)?;
                let receipt = CanonicalEncoder::<$identity, _>::new(
                    MACHINE_IDENTITY_LIMITS,
                    NeverCancel,
                )?
                .utf8(Field::new(0, "canonical-key"), &key)?
                .finish()?;
                Ok(Self {
                    canonical_key: key.into_boxed_str(),
                    receipt,
                })
            }

            /// Exact canonical key retained for diagnostics and lowering.
            #[must_use]
            pub fn canonical_key(&self) -> &str {
                &self.canonical_key
            }

            /// Domain- and role-separated durable identity.
            #[must_use]
            pub const fn identity(&self) -> $identity {
                self.receipt.id()
            }

            /// Complete canonical-preimage receipt for collision adjudication.
            #[must_use]
            pub const fn identity_receipt(&self) -> IdentityReceipt<$identity> {
                self.receipt
            }

            fn digest_bytes(&self) -> [u8; 32] {
                *self.identity().as_bytes()
            }
        }

        impl PartialEq for $name {
            fn eq(&self, other: &Self) -> bool {
                self.identity() == other.identity()
            }
        }

        impl Eq for $name {}

        impl PartialOrd for $name {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl Ord for $name {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                self.identity().cmp(&other.identity())
            }
        }

        impl Hash for $name {
            fn hash<H: Hasher>(&self, state: &mut H) {
                self.identity().hash(state);
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.canonical_key)
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.debug_struct(stringify!($name))
                    .field("canonical_key", &self.canonical_key)
                    .field("identity", &self.identity())
                    .finish()
            }
        }
    };
}

durable_graph_id!(
    /// Durable identity of a declarative machine subsystem.
    SubsystemId,
    SubsystemIdSchemaV1,
    SubsystemEntityIdV1,
    "subsystem-id",
    "org.frankensim.fs-ir.machine.subsystem-id.v1",
    "one declarative Machine-IR subsystem independent of graph serialization order"
);
durable_graph_id!(
    /// Durable identity of a machine relation.
    RelationId,
    RelationIdSchemaV1,
    RelationEntityIdV1,
    "relation-id",
    "org.frankensim.fs-ir.machine.relation-id.v1",
    "one typed Machine-IR relation independent of graph serialization order"
);
durable_graph_id!(
    /// Durable identity of a logical machine clock domain.
    ClockId,
    ClockIdSchemaV1,
    ClockEntityIdV1,
    "clock-id",
    "org.frankensim.fs-ir.machine.clock-id.v1",
    "one declared Machine-IR clock domain rather than a runtime scheduler clock"
);
durable_graph_id!(
    /// Durable identity of a machine interface binding.
    InterfaceId,
    InterfaceIdSchemaV1,
    InterfaceEntityIdV1,
    "interface-id",
    "org.frankensim.fs-ir.machine.interface-id.v1",
    "one role-oriented Machine-IR interface binding between two declared ports"
);

/// Refusal from constructing a dependency-neutral external schema reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MachineReferenceError {
    /// The reference namespace was not a canonical Machine-IR key.
    Namespace(MachineIdError),
    /// An all-zero digest cannot name an external semantic artifact.
    ZeroDigest {
        /// Nominal reference role.
        role: &'static str,
    },
}

impl MachineReferenceError {
    /// Stable diagnostic rule code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Namespace(_) => "MachineReferenceNamespace",
            Self::ZeroDigest { .. } => "MachineReferenceZeroDigest",
        }
    }
}

impl fmt::Display for MachineReferenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Namespace(error) => write!(f, "invalid external reference namespace: {error}"),
            Self::ZeroDigest { role } => write!(f, "{role} semantic digest must not be all zero"),
        }
    }
}

impl std::error::Error for MachineReferenceError {}

macro_rules! versioned_machine_ref {
    ($(#[$meta:meta])* $name:ident, $role:literal) => {
        $(#[$meta])*
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            namespace: Box<str>,
            schema_version: NonZeroU64,
            semantic_digest: [u8; 32],
        }

        impl $name {
            /// Construct an opaque, versioned semantic reference.
            ///
            /// `namespace` identifies the external schema family; the digest
            /// is retained exactly and is not rehashed or authenticated here.
            ///
            /// # Errors
            /// Refuses a noncanonical namespace or an all-zero digest.
            pub fn new(
                namespace: impl Into<String>,
                schema_version: NonZeroU64,
                semantic_digest: [u8; 32],
            ) -> Result<Self, MachineReferenceError> {
                let namespace = namespace.into();
                validate_canonical_key($role, &namespace)
                    .map_err(MachineReferenceError::Namespace)?;
                if semantic_digest == [0; 32] {
                    return Err(MachineReferenceError::ZeroDigest { role: $role });
                }
                Ok(Self {
                    namespace: namespace.into_boxed_str(),
                    schema_version,
                    semantic_digest,
                })
            }

            /// External schema namespace.
            #[must_use]
            pub fn namespace(&self) -> &str {
                &self.namespace
            }

            /// Explicit external schema version.
            #[must_use]
            pub const fn schema_version(&self) -> NonZeroU64 {
                self.schema_version
            }

            /// Exact semantic digest supplied by the external owner.
            #[must_use]
            pub const fn semantic_digest(&self) -> [u8; 32] {
                self.semantic_digest
            }

            fn append_canonical(&self, out: &mut Vec<u8>) {
                push_len_prefixed(out, self.namespace.as_bytes());
                out.extend_from_slice(&self.schema_version.get().to_le_bytes());
                out.extend_from_slice(&self.semantic_digest);
            }
        }
    };
}

versioned_machine_ref!(
    /// Opaque versioned reference to an executable subsystem model schema.
    ModelRef,
    "model-ref"
);
versioned_machine_ref!(
    /// Opaque versioned reference to an immutable material-card schema.
    MaterialCardRef,
    "material-card-ref"
);
versioned_machine_ref!(
    /// Opaque versioned reference to an interface-system card schema.
    InterfaceCardRef,
    "interface-card-ref"
);
versioned_machine_ref!(
    /// Opaque versioned reference to an algebraic-loop solve-policy schema.
    SolvePolicyRef,
    "solve-policy-ref"
);

/// Logical clock behavior declared by a machine graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MachineClock {
    /// Continuous-time semantics; this is not a wall-clock guarantee.
    Continuous,
    /// Periodic logical sampling in integer nanoseconds.
    Periodic {
        /// Strictly positive logical period.
        period_ns: NonZeroU64,
        /// Phase in `[0, period_ns)`.
        phase_ns: u64,
    },
    /// Event-driven logical time; event/reset semantics land in PR-3.
    EventDriven,
}

/// One named logical clock declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClockSpec {
    /// Stable clock-domain identity.
    pub id: ClockId,
    /// Declared logical behavior.
    pub clock: MachineClock,
}

/// Quantity contract carried by a terminal.
///
/// `Dimensional` is an explicit no-semantic-kind claim, not a wildcard. It is
/// compatible only with the exact same dimensional declaration. `Semantic`
/// additionally distinguishes equal-dimensional meanings such as pressure
/// versus stress and absolute versus difference temperature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TerminalQuantitySpec {
    /// Six-base dimension vector with no stronger semantic-kind claim.
    Dimensional(Dims),
    /// Sealed semantic quantity kind and scalar value form.
    Semantic(SemanticType),
}

impl TerminalQuantitySpec {
    /// Exact six-base dimension vector required by the terminal.
    #[must_use]
    pub const fn dims(self) -> Dims {
        match self {
            Self::Dimensional(dims) => dims,
            Self::Semantic(semantic_type) => semantic_type.expected_dims(),
        }
    }

    /// Strong semantic type, when one was declared.
    #[must_use]
    pub const fn semantic_type(self) -> Option<SemanticType> {
        match self {
            Self::Dimensional(_) => None,
            Self::Semantic(semantic_type) => Some(semantic_type),
        }
    }

    fn is_admitted(self) -> bool {
        match self {
            Self::Dimensional(_) => true,
            Self::Semantic(semantic_type) => semantic_type
                .kind()
                .admits_scalar_form(semantic_type.form()),
        }
    }
}

/// Declared terminal value shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TerminalShape {
    /// One real scalar.
    Scalar,
    /// Fixed-width real vector.
    Vector {
        /// Number of components.
        components: NonZeroU64,
    },
    /// Fixed rectangular real tensor.
    Tensor {
        /// Row count.
        rows: NonZeroU64,
        /// Column count.
        columns: NonZeroU64,
    },
    /// Fixed-component trace field on an external geometric support.
    FieldTrace {
        /// Number of value components at each trace point.
        components: NonZeroU64,
    },
}

/// Orientation of terminal coordinates relative to their named frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OrientationParity {
    /// Coordinate orientation preserves the named frame orientation.
    Preserving,
    /// Coordinate orientation reverses the named frame orientation.
    Reversing,
}

/// Canonical frame and orientation declaration for a terminal.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FrameBinding {
    canonical_key: Box<str>,
    orientation: OrientationParity,
}

impl FrameBinding {
    /// Construct a frame binding from a canonical frame key.
    ///
    /// # Errors
    /// Refuses a noncanonical key before it can enter graph identity.
    pub fn new(
        canonical_key: impl Into<String>,
        orientation: OrientationParity,
    ) -> Result<Self, MachineIdError> {
        let canonical_key = canonical_key.into();
        validate_canonical_key("frame-binding", &canonical_key)?;
        Ok(Self {
            canonical_key: canonical_key.into_boxed_str(),
            orientation,
        })
    }

    /// Canonical frame key.
    #[must_use]
    pub fn canonical_key(&self) -> &str {
        &self.canonical_key
    }

    /// Explicit orientation parity within the named frame.
    #[must_use]
    pub const fn orientation(&self) -> OrientationParity {
        self.orientation
    }
}

/// Causality role of a subsystem terminal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TerminalCausality {
    /// Value must be closed by exactly one graph relation.
    Input,
    /// Value is produced by the owning subsystem model.
    Output,
    /// Value is supplied across the admitted graph boundary.
    ExternalInput,
}

/// One typed subsystem terminal declaration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TerminalSpec {
    /// Stable terminal identity.
    pub id: TerminalId,
    /// Owning subsystem.
    pub owner: SubsystemId,
    /// Quantity dimensions and optional stronger semantic kind.
    pub quantity: TerminalQuantitySpec,
    /// Scalar/vector/tensor/trace shape.
    pub shape: TerminalShape,
    /// Input/output boundary role.
    pub causality: TerminalCausality,
    /// Logical clock domain.
    pub clock: ClockId,
    /// Coordinate frame and orientation.
    pub frame: FrameBinding,
}

/// Declarative subsystem and directly owned durable elements.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubsystemSpec {
    /// Stable subsystem identity.
    pub id: SubsystemId,
    /// Opaque versioned external model reference.
    pub model: ModelRef,
    /// Directly owned body occurrences.
    pub bodies: Vec<BodyId>,
    /// Directly owned named surface supports.
    pub surface_patches: Vec<SurfacePatchId>,
    /// Directly owned contact features.
    pub contact_features: Vec<ContactFeatureId>,
    /// Directly owned state slots.
    pub state_slots: Vec<StateSlotId>,
}

/// Port-side energy orientation relative to the owning subsystem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PortEnergyRole {
    /// Positive effort-flow product is directed into the subsystem.
    IntoSubsystem,
    /// Positive effort-flow product is directed out of the subsystem.
    OutOfSubsystem,
}

/// Neutral effort/flow port declaration.
///
/// This is graph data only. Executable transport, interpolation, buffering,
/// and synchronization remain owned by `fs-couple`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PortSpec {
    /// Stable coupling-port identity.
    pub id: PortId,
    /// Owning subsystem.
    pub owner: SubsystemId,
    /// Terminal carrying the effort variable.
    pub effort: TerminalId,
    /// Terminal carrying the flow variable.
    pub flow: TerminalId,
    /// Explicit sign convention for power accounting.
    pub energy_role: PortEnergyRole,
}

/// Timing/solve boundary of one directed relation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RelationMode {
    /// Instantaneous dependency. A named policy cuts the dependency edge for
    /// structural algebraic-loop admission; it does not prove convergence.
    Algebraic {
        /// Explicit external solve-policy reference, when present.
        solve_policy: Option<SolvePolicyRef>,
    },
    /// Stateful dependency whose declared state slot breaks feedthrough.
    Stateful {
        /// State written by this relation and owned by the target subsystem.
        state_slot: StateSlotId,
    },
}

/// One directed typed relation from an output to an input terminal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RelationSpec {
    /// Stable relation identity.
    pub id: RelationId,
    /// Source output terminal.
    pub source: TerminalId,
    /// Target input terminal.
    pub target: TerminalId,
    /// Algebraic or state-breaking dependency semantics.
    pub mode: RelationMode,
}

/// Durable target of a material-card binding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MaterialTarget {
    /// Complete body occurrence.
    Body(BodyId),
    /// Named surface support with a distinct surface material declaration.
    SurfacePatch(SurfacePatchId),
}

/// Opaque material-card binding to one owned durable target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaterialBinding {
    /// Owned body or surface target.
    pub target: MaterialTarget,
    /// Versioned immutable-card semantic reference.
    pub material: MaterialCardRef,
}

/// Physical orientation of an interface's negative and positive port roles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum InterfaceOrientation {
    /// Endpoint coordinate orientations must match.
    Aligned,
    /// Endpoint coordinate orientations must be opposite.
    Opposed,
}

/// Role-oriented interface binding between two ports.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InterfaceBinding {
    /// Stable interface identity.
    pub id: InterfaceId,
    /// Negative-role port; endpoint order is identity-semantic.
    pub negative: PortId,
    /// Positive-role port; endpoint order is identity-semantic.
    pub positive: PortId,
    /// Versioned external interface-system reference.
    pub interface: InterfaceCardRef,
    /// Explicit endpoint orientation law.
    pub orientation: InterfaceOrientation,
}

/// Canonical identity schema for one admitted machine graph.
pub enum MachineGraphIdentitySchemaV1 {}

impl CanonicalSchema for MachineGraphIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-ir.machine.graph.v1";
    const NAME: &'static str = "admitted-machine-graph";
    const VERSION: u32 = MACHINE_GRAPH_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "complete canonical subsystem, clock, terminal, port, relation, material, and interface declarations";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("graph-schema-version", WireType::U64),
        FieldSpec::required("clocks", WireType::OrderedBytes),
        FieldSpec::required("subsystems", WireType::OrderedBytes),
        FieldSpec::required("terminals", WireType::OrderedBytes),
        FieldSpec::required("ports", WireType::OrderedBytes),
        FieldSpec::required("relations", WireType::OrderedBytes),
        FieldSpec::required("materials", WireType::OrderedBytes),
        FieldSpec::required("interfaces", WireType::OrderedBytes),
    ];
}

/// Strong semantic identity of an admitted Machine-IR graph.
pub type MachineGraphIdV1 = ProblemSemanticId<MachineGraphIdentitySchemaV1>;

/// Closed admission rule vocabulary for deterministic graph diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum MachineGraphRule {
    /// A public collection or aggregate ownership limit was exceeded.
    ResourceLimit = 1,
    /// A clock identity was declared more than once.
    DuplicateClock = 2,
    /// A periodic clock phase was not strictly below its period.
    InvalidClockPhase = 3,
    /// A subsystem identity was declared more than once.
    DuplicateSubsystem = 4,
    /// A durable topology/state element had more than one owner.
    DuplicateOwnership = 5,
    /// A terminal identity was declared more than once.
    DuplicateTerminal = 6,
    /// A terminal named an undeclared subsystem owner.
    UnknownTerminalOwner = 7,
    /// A terminal named an undeclared clock.
    UnknownTerminalClock = 8,
    /// A semantic terminal kind used a forbidden scalar form.
    UnsupportedTerminalForm = 9,
    /// A port identity was declared more than once.
    DuplicatePort = 10,
    /// A port named an undeclared subsystem owner.
    UnknownPortOwner = 11,
    /// A port named an undeclared terminal.
    UnknownPortTerminal = 12,
    /// A port terminal belonged to a different subsystem.
    PortOwnerMismatch = 13,
    /// A port reused one terminal for both energy roles or across ports.
    PortTerminalConflict = 14,
    /// Effort and flow terminals used different clock domains.
    PortClockMismatch = 15,
    /// Effort and flow terminals used different frames.
    PortFrameMismatch = 16,
    /// Effort and flow terminal orientations disagreed.
    PortOrientationMismatch = 17,
    /// Effort and flow dimensions did not multiply to physical power.
    PortPowerDimensionMismatch = 18,
    /// A relation identity was declared more than once.
    DuplicateRelation = 19,
    /// A relation named an undeclared terminal.
    UnknownRelationTerminal = 20,
    /// Relation endpoints did not form output-to-input causality.
    RelationCausalityGap = 21,
    /// Relation endpoints had different quantity contracts.
    RelationQuantityGap = 22,
    /// Relation endpoints had different value shapes.
    RelationShapeGap = 23,
    /// Relation endpoints had different clocks.
    RelationClockGap = 24,
    /// Relation endpoints had different frames.
    RelationFrameGap = 25,
    /// Relation endpoints had different orientation declarations.
    RelationOrientationGap = 26,
    /// One input terminal had more than one producer relation.
    MultipleInputSources = 27,
    /// One internal input terminal had no producer relation.
    MissingSourceClosure = 28,
    /// A stateful relation named an undeclared state slot.
    UnknownStateSlot = 29,
    /// The named state slot was owned outside the target subsystem.
    StateOwnerMismatch = 30,
    /// More than one relation wrote the same state slot.
    MultipleStateWriters = 31,
    /// A declared state slot had no explicit stateful writer.
    UnaccountedState = 32,
    /// Instantaneous ungoverned dependencies contained a directed cycle.
    AlgebraicLoopWithoutSolvePolicy = 33,
    /// More than one material binding named the same target.
    DuplicateMaterialBinding = 34,
    /// A material binding named an unowned target.
    UnknownMaterialTarget = 35,
    /// An owned body had no material-card binding.
    MissingBodyMaterial = 36,
    /// An interface identity was declared more than once.
    DuplicateInterface = 37,
    /// An interface named an undeclared port.
    UnknownInterfacePort = 38,
    /// An interface used the same port for both endpoint roles.
    SameInterfaceEndpoint = 39,
    /// One port was bound into more than one interface.
    PortInterfaceConflict = 40,
    /// Interface endpoint effort/flow quantity contracts disagreed.
    InterfaceQuantityGap = 41,
    /// Interface endpoint effort/flow shapes disagreed.
    InterfaceShapeGap = 42,
    /// Interface endpoint clocks disagreed.
    InterfaceClockGap = 43,
    /// Interface endpoint frames disagreed.
    InterfaceFrameGap = 44,
    /// Interface endpoint orientation did not satisfy the declared law.
    InterfaceOrientationGap = 45,
    /// Interface endpoint power-direction roles were not complementary.
    InterfaceEnergyRoleGap = 46,
    /// Canonical semantic identity construction refused publication.
    Identity = 47,
    /// Effort and flow terminals used incompatible value shapes.
    PortShapeMismatch = 48,
    /// An interface terminal pair did not contain exactly one output and one input.
    InterfaceCausalityGap = 49,
    /// Interface endpoints lacked their exact directed terminal relations.
    InterfaceRelationGap = 50,
}

impl MachineGraphRule {
    /// Stable rule code for structured admission logging.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ResourceLimit => "MachineGraphResourceLimit",
            Self::DuplicateClock => "MachineGraphDuplicateClock",
            Self::InvalidClockPhase => "MachineGraphInvalidClockPhase",
            Self::DuplicateSubsystem => "MachineGraphDuplicateSubsystem",
            Self::DuplicateOwnership => "MachineGraphDuplicateOwnership",
            Self::DuplicateTerminal => "MachineGraphDuplicateTerminal",
            Self::UnknownTerminalOwner => "MachineGraphUnknownTerminalOwner",
            Self::UnknownTerminalClock => "MachineGraphUnknownTerminalClock",
            Self::UnsupportedTerminalForm => "MachineGraphUnsupportedTerminalForm",
            Self::DuplicatePort => "MachineGraphDuplicatePort",
            Self::UnknownPortOwner => "MachineGraphUnknownPortOwner",
            Self::UnknownPortTerminal => "MachineGraphUnknownPortTerminal",
            Self::PortOwnerMismatch => "MachineGraphPortOwnerMismatch",
            Self::PortTerminalConflict => "MachineGraphPortTerminalConflict",
            Self::PortClockMismatch => "MachineGraphPortClockMismatch",
            Self::PortFrameMismatch => "MachineGraphPortFrameMismatch",
            Self::PortOrientationMismatch => "MachineGraphPortOrientationMismatch",
            Self::PortPowerDimensionMismatch => "MachineGraphPortPowerDimensionMismatch",
            Self::DuplicateRelation => "MachineGraphDuplicateRelation",
            Self::UnknownRelationTerminal => "MachineGraphUnknownRelationTerminal",
            Self::RelationCausalityGap => "MachineGraphRelationCausalityGap",
            Self::RelationQuantityGap => "MachineGraphRelationQuantityGap",
            Self::RelationShapeGap => "MachineGraphRelationShapeGap",
            Self::RelationClockGap => "MachineGraphRelationClockGap",
            Self::RelationFrameGap => "MachineGraphRelationFrameGap",
            Self::RelationOrientationGap => "MachineGraphRelationOrientationGap",
            Self::MultipleInputSources => "MachineGraphMultipleInputSources",
            Self::MissingSourceClosure => "MachineGraphMissingSourceClosure",
            Self::UnknownStateSlot => "MachineGraphUnknownStateSlot",
            Self::StateOwnerMismatch => "MachineGraphStateOwnerMismatch",
            Self::MultipleStateWriters => "MachineGraphMultipleStateWriters",
            Self::UnaccountedState => "MachineGraphUnaccountedState",
            Self::AlgebraicLoopWithoutSolvePolicy => "MachineGraphAlgebraicLoopWithoutSolvePolicy",
            Self::DuplicateMaterialBinding => "MachineGraphDuplicateMaterialBinding",
            Self::UnknownMaterialTarget => "MachineGraphUnknownMaterialTarget",
            Self::MissingBodyMaterial => "MachineGraphMissingBodyMaterial",
            Self::DuplicateInterface => "MachineGraphDuplicateInterface",
            Self::UnknownInterfacePort => "MachineGraphUnknownInterfacePort",
            Self::SameInterfaceEndpoint => "MachineGraphSameInterfaceEndpoint",
            Self::PortInterfaceConflict => "MachineGraphPortInterfaceConflict",
            Self::InterfaceQuantityGap => "MachineGraphInterfaceQuantityGap",
            Self::InterfaceShapeGap => "MachineGraphInterfaceShapeGap",
            Self::InterfaceClockGap => "MachineGraphInterfaceClockGap",
            Self::InterfaceFrameGap => "MachineGraphInterfaceFrameGap",
            Self::InterfaceOrientationGap => "MachineGraphInterfaceOrientationGap",
            Self::InterfaceEnergyRoleGap => "MachineGraphInterfaceEnergyRoleGap",
            Self::Identity => "MachineGraphIdentity",
            Self::PortShapeMismatch => "MachineGraphPortShapeMismatch",
            Self::InterfaceCausalityGap => "MachineGraphInterfaceCausalityGap",
            Self::InterfaceRelationGap => "MachineGraphInterfaceRelationGap",
        }
    }
}

/// Typed subject named by a graph-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineGraphSubject {
    /// Complete submitted graph.
    Graph,
    /// Clock declaration.
    Clock(ClockId),
    /// Subsystem declaration.
    Subsystem(SubsystemId),
    /// Durable topology or state element.
    MachineElement(MachineElementId),
    /// Terminal declaration.
    Terminal(TerminalId),
    /// Port declaration.
    Port(PortId),
    /// Relation declaration.
    Relation(RelationId),
    /// Material binding target.
    MaterialTarget(MaterialTarget),
    /// Interface binding.
    Interface(InterfaceId),
}

/// One deterministic, structured graph-admission finding.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MachineGraphFinding {
    rule: MachineGraphRule,
    subject: MachineGraphSubject,
    related: Option<MachineGraphSubject>,
}

impl MachineGraphFinding {
    fn new(
        rule: MachineGraphRule,
        subject: MachineGraphSubject,
        related: Option<MachineGraphSubject>,
    ) -> Self {
        Self {
            rule,
            subject,
            related,
        }
    }

    /// Closed rule that produced this finding.
    #[must_use]
    pub const fn rule(&self) -> MachineGraphRule {
        self.rule
    }

    /// Stable rule code for structured traces.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        self.rule.code()
    }

    /// Exact offending graph subject.
    #[must_use]
    pub const fn subject(&self) -> &MachineGraphSubject {
        &self.subject
    }

    /// Optional second subject needed to explain the conflict.
    #[must_use]
    pub const fn related(&self) -> Option<&MachineGraphSubject> {
        self.related.as_ref()
    }
}

/// Complete, deterministically sorted refusal from graph admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineGraphRefusal {
    findings: Vec<MachineGraphFinding>,
    identity_error: Option<CanonicalError>,
}

impl MachineGraphRefusal {
    /// Sorted, duplicate-free findings. This slice is always nonempty.
    #[must_use]
    pub fn findings(&self) -> &[MachineGraphFinding] {
        &self.findings
    }

    /// Canonical identity error, only for `MachineGraphIdentity`.
    #[must_use]
    pub const fn identity_error(&self) -> Option<&CanonicalError> {
        self.identity_error.as_ref()
    }

    /// Stable outcome code; inspect individual findings for specific rules.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        "MachineGraphRefused"
    }
}

impl fmt::Display for MachineGraphRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "machine graph refused with {} finding(s); first rule is {}",
            self.findings.len(),
            self.findings[0].code()
        )
    }
}

impl std::error::Error for MachineGraphRefusal {}

/// Collection sizes retained for every graph-admission attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MachineGraphSubmittedCounts {
    /// Submitted clock count.
    pub clocks: usize,
    /// Submitted subsystem count.
    pub subsystems: usize,
    /// Submitted terminal count.
    pub terminals: usize,
    /// Submitted port count.
    pub ports: usize,
    /// Submitted relation count.
    pub relations: usize,
    /// Submitted material-binding count.
    pub materials: usize,
    /// Submitted interface-binding count.
    pub interfaces: usize,
    /// Submitted owned-element count across subsystems.
    pub owned_elements: usize,
}

/// Mutable-by-construction, authority-free machine graph draft.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineGraphDraft {
    /// Logical clock declarations.
    pub clocks: Vec<ClockSpec>,
    /// Subsystem/model declarations and durable ownership.
    pub subsystems: Vec<SubsystemSpec>,
    /// Typed terminal declarations.
    pub terminals: Vec<TerminalSpec>,
    /// Neutral effort/flow port declarations.
    pub ports: Vec<PortSpec>,
    /// Directed typed relations.
    pub relations: Vec<RelationSpec>,
    /// Material-card bindings.
    pub materials: Vec<MaterialBinding>,
    /// Role-oriented interface bindings.
    pub interfaces: Vec<InterfaceBinding>,
}

impl MachineGraphDraft {
    /// Attempt graph admission and publish no semantic identity on refusal.
    ///
    /// # Errors
    /// Returns all deterministic findings discovered within the bounded draft.
    pub fn admit(self) -> Result<AdmittedMachineGraph, MachineGraphRefusal> {
        self.admit_with_decision().into_result()
    }

    /// Attempt admission while retaining structured submitted counts.
    #[must_use]
    pub fn admit_with_decision(self) -> MachineGraphAdmissionDecision {
        let submitted = submitted_counts(&self);
        MachineGraphAdmissionDecision {
            submitted,
            result: admit_machine_graph(self),
        }
    }
}

/// Canonically ordered admitted machine graph plus its semantic receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedMachineGraph {
    clocks: Vec<ClockSpec>,
    subsystems: Vec<SubsystemSpec>,
    terminals: Vec<TerminalSpec>,
    ports: Vec<PortSpec>,
    relations: Vec<RelationSpec>,
    materials: Vec<MaterialBinding>,
    interfaces: Vec<InterfaceBinding>,
    receipt: IdentityReceipt<MachineGraphIdV1>,
}

impl AdmittedMachineGraph {
    /// Strong semantic identity of the complete admitted graph.
    #[must_use]
    pub const fn identity(&self) -> MachineGraphIdV1 {
        self.receipt.id()
    }

    /// Complete canonical identity receipt for collision adjudication.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<MachineGraphIdV1> {
        self.receipt
    }

    /// Canonically ordered clocks.
    #[must_use]
    pub fn clocks(&self) -> &[ClockSpec] {
        &self.clocks
    }

    /// Canonically ordered subsystems.
    #[must_use]
    pub fn subsystems(&self) -> &[SubsystemSpec] {
        &self.subsystems
    }

    /// Canonically ordered terminals.
    #[must_use]
    pub fn terminals(&self) -> &[TerminalSpec] {
        &self.terminals
    }

    /// Canonically ordered ports.
    #[must_use]
    pub fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    /// Canonically ordered relations.
    #[must_use]
    pub fn relations(&self) -> &[RelationSpec] {
        &self.relations
    }

    /// Canonically ordered material bindings.
    #[must_use]
    pub fn materials(&self) -> &[MaterialBinding] {
        &self.materials
    }

    /// Canonically ordered interface bindings.
    #[must_use]
    pub fn interfaces(&self) -> &[InterfaceBinding] {
        &self.interfaces
    }
}

/// Bounded deterministic outcome summary for one graph-admission attempt.
///
/// This summary is suitable for structured tracing but is not itself a digest
/// or replay record of an early-refused draft.
#[must_use]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MachineGraphAdmissionDecision {
    submitted: MachineGraphSubmittedCounts,
    result: Result<AdmittedMachineGraph, MachineGraphRefusal>,
}

impl MachineGraphAdmissionDecision {
    /// Exact collection counts observed before canonicalization.
    #[must_use]
    pub const fn submitted_counts(&self) -> MachineGraphSubmittedCounts {
        self.submitted
    }

    /// Stable top-level decision code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match &self.result {
            Ok(_) => "MachineGraphAdmitted",
            Err(_) => "MachineGraphRefused",
        }
    }

    /// Borrow the admitted graph or complete refusal.
    #[must_use]
    pub fn result(&self) -> Result<&AdmittedMachineGraph, &MachineGraphRefusal> {
        self.result.as_ref()
    }

    /// Consume the decision and recover the conventional result.
    #[must_use]
    pub fn into_result(self) -> Result<AdmittedMachineGraph, MachineGraphRefusal> {
        self.result
    }
}

fn submitted_counts(draft: &MachineGraphDraft) -> MachineGraphSubmittedCounts {
    let owned_elements = draft.subsystems.iter().fold(0usize, |count, subsystem| {
        count
            .saturating_add(subsystem.bodies.len())
            .saturating_add(subsystem.surface_patches.len())
            .saturating_add(subsystem.contact_features.len())
            .saturating_add(subsystem.state_slots.len())
    });
    MachineGraphSubmittedCounts {
        clocks: draft.clocks.len(),
        subsystems: draft.subsystems.len(),
        terminals: draft.terminals.len(),
        ports: draft.ports.len(),
        relations: draft.relations.len(),
        materials: draft.materials.len(),
        interfaces: draft.interfaces.len(),
        owned_elements,
    }
}

fn resource_limit_findings(counts: MachineGraphSubmittedCounts) -> Vec<MachineGraphFinding> {
    let over_limit = counts.clocks > MAX_MACHINE_GRAPH_CLOCKS
        || counts.subsystems > MAX_MACHINE_GRAPH_SUBSYSTEMS
        || counts.terminals > MAX_MACHINE_GRAPH_TERMINALS
        || counts.ports > MAX_MACHINE_GRAPH_PORTS
        || counts.relations > MAX_MACHINE_GRAPH_RELATIONS
        || counts.materials > MAX_MACHINE_GRAPH_MATERIALS
        || counts.interfaces > MAX_MACHINE_GRAPH_INTERFACES
        || counts.owned_elements > MAX_MACHINE_GRAPH_OWNED_ELEMENTS;
    if over_limit {
        vec![MachineGraphFinding::new(
            MachineGraphRule::ResourceLimit,
            MachineGraphSubject::Graph,
            None,
        )]
    } else {
        Vec::new()
    }
}

fn graph_refusal(
    mut findings: Vec<MachineGraphFinding>,
    identity_error: Option<CanonicalError>,
) -> MachineGraphRefusal {
    findings.sort();
    findings.dedup();
    debug_assert!(!findings.is_empty());
    MachineGraphRefusal {
        findings,
        identity_error,
    }
}

// The admission sequence intentionally mirrors the closed rule vocabulary so
// every finding remains local, ordered, and auditable in one pass.
#[allow(clippy::too_many_lines)]
fn admit_machine_graph(
    mut draft: MachineGraphDraft,
) -> Result<AdmittedMachineGraph, MachineGraphRefusal> {
    let counts = submitted_counts(&draft);
    let resource_findings = resource_limit_findings(counts);
    if !resource_findings.is_empty() {
        return Err(graph_refusal(resource_findings, None));
    }

    canonicalize_graph_draft(&mut draft);
    let mut findings = Vec::new();

    for pair in draft.clocks.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateClock,
                MachineGraphSubject::Clock(pair[1].id.clone()),
                None,
            ));
        }
    }
    for clock in &draft.clocks {
        if let MachineClock::Periodic {
            period_ns,
            phase_ns,
        } = clock.clock
            && phase_ns >= period_ns.get()
        {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::InvalidClockPhase,
                MachineGraphSubject::Clock(clock.id.clone()),
                None,
            ));
        }
    }
    let clock_ids: BTreeSet<ClockId> = draft.clocks.iter().map(|clock| clock.id.clone()).collect();

    for pair in draft.subsystems.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateSubsystem,
                MachineGraphSubject::Subsystem(pair[1].id.clone()),
                None,
            ));
        }
    }
    let subsystem_ids: BTreeSet<SubsystemId> = draft
        .subsystems
        .iter()
        .map(|subsystem| subsystem.id.clone())
        .collect();
    let subsystem_indices: BTreeMap<SubsystemId, usize> = draft
        .subsystems
        .iter()
        .enumerate()
        .map(|(index, subsystem)| (subsystem.id.clone(), index))
        .collect();

    let mut element_owners = BTreeMap::<MachineElementId, SubsystemId>::new();
    let mut state_owners = BTreeMap::<StateSlotId, SubsystemId>::new();
    let mut body_ids = BTreeSet::<BodyId>::new();
    for subsystem in &draft.subsystems {
        for body in &subsystem.bodies {
            body_ids.insert(body.clone());
            record_element_owner(
                &mut element_owners,
                MachineElementId::Body(body.clone()),
                &subsystem.id,
                &mut findings,
            );
        }
        for patch in &subsystem.surface_patches {
            record_element_owner(
                &mut element_owners,
                MachineElementId::SurfacePatch(patch.clone()),
                &subsystem.id,
                &mut findings,
            );
        }
        for feature in &subsystem.contact_features {
            record_element_owner(
                &mut element_owners,
                MachineElementId::ContactFeature(feature.clone()),
                &subsystem.id,
                &mut findings,
            );
        }
        for state_slot in &subsystem.state_slots {
            record_element_owner(
                &mut element_owners,
                MachineElementId::StateSlot(state_slot.clone()),
                &subsystem.id,
                &mut findings,
            );
            if let Some(first_owner) = state_owners.insert(state_slot.clone(), subsystem.id.clone())
            {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::DuplicateOwnership,
                    MachineGraphSubject::MachineElement(MachineElementId::StateSlot(
                        state_slot.clone(),
                    )),
                    Some(MachineGraphSubject::Subsystem(first_owner)),
                ));
            }
        }
    }

    for pair in draft.terminals.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateTerminal,
                MachineGraphSubject::Terminal(pair[1].id.clone()),
                None,
            ));
        }
    }
    let terminal_indices: BTreeMap<TerminalId, usize> = draft
        .terminals
        .iter()
        .enumerate()
        .map(|(index, terminal)| (terminal.id.clone(), index))
        .collect();
    for terminal in &draft.terminals {
        record_element_owner(
            &mut element_owners,
            MachineElementId::Terminal(terminal.id.clone()),
            &terminal.owner,
            &mut findings,
        );
        if !subsystem_ids.contains(&terminal.owner) {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownTerminalOwner,
                MachineGraphSubject::Terminal(terminal.id.clone()),
                Some(MachineGraphSubject::Subsystem(terminal.owner.clone())),
            ));
        }
        if !clock_ids.contains(&terminal.clock) {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownTerminalClock,
                MachineGraphSubject::Terminal(terminal.id.clone()),
                Some(MachineGraphSubject::Clock(terminal.clock.clone())),
            ));
        }
        if !terminal.quantity.is_admitted() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnsupportedTerminalForm,
                MachineGraphSubject::Terminal(terminal.id.clone()),
                None,
            ));
        }
    }

    for pair in draft.ports.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicatePort,
                MachineGraphSubject::Port(pair[1].id.clone()),
                None,
            ));
        }
    }
    let port_indices: BTreeMap<PortId, usize> = draft
        .ports
        .iter()
        .enumerate()
        .map(|(index, port)| (port.id.clone(), index))
        .collect();
    let mut terminal_ports = BTreeMap::<TerminalId, PortId>::new();
    for port in &draft.ports {
        record_element_owner(
            &mut element_owners,
            MachineElementId::Port(port.id.clone()),
            &port.owner,
            &mut findings,
        );
        if !subsystem_ids.contains(&port.owner) {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownPortOwner,
                MachineGraphSubject::Port(port.id.clone()),
                Some(MachineGraphSubject::Subsystem(port.owner.clone())),
            ));
        }
        if port.effort == port.flow {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::PortTerminalConflict,
                MachineGraphSubject::Port(port.id.clone()),
                Some(MachineGraphSubject::Terminal(port.effort.clone())),
            ));
        }
        for terminal_id in [&port.effort, &port.flow] {
            if let Some(first_port) = terminal_ports.insert(terminal_id.clone(), port.id.clone()) {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortTerminalConflict,
                    MachineGraphSubject::Port(port.id.clone()),
                    Some(MachineGraphSubject::Port(first_port)),
                ));
            }
        }
        let effort = terminal_indices
            .get(&port.effort)
            .map(|index| &draft.terminals[*index]);
        let flow = terminal_indices
            .get(&port.flow)
            .map(|index| &draft.terminals[*index]);
        if effort.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownPortTerminal,
                MachineGraphSubject::Port(port.id.clone()),
                Some(MachineGraphSubject::Terminal(port.effort.clone())),
            ));
        }
        if flow.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownPortTerminal,
                MachineGraphSubject::Port(port.id.clone()),
                Some(MachineGraphSubject::Terminal(port.flow.clone())),
            ));
        }
        if let (Some(effort), Some(flow)) = (effort, flow) {
            for terminal in [effort, flow] {
                if terminal.owner != port.owner {
                    findings.push(MachineGraphFinding::new(
                        MachineGraphRule::PortOwnerMismatch,
                        MachineGraphSubject::Port(port.id.clone()),
                        Some(MachineGraphSubject::Terminal(terminal.id.clone())),
                    ));
                }
            }
            if effort.clock != flow.clock {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortClockMismatch,
                    MachineGraphSubject::Port(port.id.clone()),
                    None,
                ));
            }
            if effort.frame.canonical_key() != flow.frame.canonical_key() {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortFrameMismatch,
                    MachineGraphSubject::Port(port.id.clone()),
                    None,
                ));
            }
            if effort.frame.orientation() != flow.frame.orientation() {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortOrientationMismatch,
                    MachineGraphSubject::Port(port.id.clone()),
                    None,
                ));
            }
            if effort.shape != flow.shape {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortShapeMismatch,
                    MachineGraphSubject::Port(port.id.clone()),
                    None,
                ));
            }
            if effort.quantity.dims().checked_plus(flow.quantity.dims()) != Some(POWER_DIMS) {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortPowerDimensionMismatch,
                    MachineGraphSubject::Port(port.id.clone()),
                    None,
                ));
            }
        }
    }

    let oriented_interface_relations = declared_interface_relation_orientations(
        &draft.interfaces,
        &draft.ports,
        &port_indices,
        &draft.terminals,
        &terminal_indices,
    );

    for pair in draft.relations.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateRelation,
                MachineGraphSubject::Relation(pair[1].id.clone()),
                None,
            ));
        }
    }
    let mut target_source_counts = BTreeMap::<TerminalId, usize>::new();
    let mut state_writer_counts = BTreeMap::<StateSlotId, usize>::new();
    let mut algebraic_adjacency = vec![Vec::<(usize, RelationId)>::new(); draft.subsystems.len()];
    for relation in &draft.relations {
        let source = terminal_indices
            .get(&relation.source)
            .map(|index| &draft.terminals[*index]);
        let target = terminal_indices
            .get(&relation.target)
            .map(|index| &draft.terminals[*index]);
        if source.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownRelationTerminal,
                MachineGraphSubject::Relation(relation.id.clone()),
                Some(MachineGraphSubject::Terminal(relation.source.clone())),
            ));
        }
        if target.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownRelationTerminal,
                MachineGraphSubject::Relation(relation.id.clone()),
                Some(MachineGraphSubject::Terminal(relation.target.clone())),
            ));
        }
        if let Some(target) = target {
            *target_source_counts.entry(target.id.clone()).or_default() += 1;
        }
        if let (Some(source), Some(target)) = (source, target) {
            if source.causality != TerminalCausality::Output
                || target.causality != TerminalCausality::Input
            {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationCausalityGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            if source.quantity != target.quantity {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationQuantityGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            if source.shape != target.shape {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationShapeGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            if source.clock != target.clock {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationClockGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            if source.frame.canonical_key() != target.frame.canonical_key() {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationFrameGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            let orientations_match = source.frame.orientation() == target.frame.orientation();
            let orientation_admitted =
                match oriented_interface_relations.get(&(source.id.clone(), target.id.clone())) {
                    Some(InterfaceOrientation::Aligned) => orientations_match,
                    Some(InterfaceOrientation::Opposed) => !orientations_match,
                    None => orientations_match,
                };
            if !orientation_admitted {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::RelationOrientationGap,
                    MachineGraphSubject::Relation(relation.id.clone()),
                    Some(MachineGraphSubject::Terminal(target.id.clone())),
                ));
            }
            match &relation.mode {
                RelationMode::Stateful { state_slot } => {
                    *state_writer_counts.entry(state_slot.clone()).or_default() += 1;
                    match state_owners.get(state_slot) {
                        None => findings.push(MachineGraphFinding::new(
                            MachineGraphRule::UnknownStateSlot,
                            MachineGraphSubject::Relation(relation.id.clone()),
                            Some(MachineGraphSubject::MachineElement(
                                MachineElementId::StateSlot(state_slot.clone()),
                            )),
                        )),
                        Some(owner) if *owner != target.owner => {
                            findings.push(MachineGraphFinding::new(
                                MachineGraphRule::StateOwnerMismatch,
                                MachineGraphSubject::Relation(relation.id.clone()),
                                Some(MachineGraphSubject::Subsystem(owner.clone())),
                            ));
                        }
                        Some(_) => {}
                    }
                }
                RelationMode::Algebraic { solve_policy: None } => {
                    if let (Some(source_owner), Some(target_owner)) = (
                        subsystem_indices.get(&source.owner),
                        subsystem_indices.get(&target.owner),
                    ) {
                        algebraic_adjacency[*source_owner]
                            .push((*target_owner, relation.id.clone()));
                    }
                }
                RelationMode::Algebraic {
                    solve_policy: Some(_),
                } => {}
            }
        }
    }
    for terminal in &draft.terminals {
        if terminal.causality == TerminalCausality::Input {
            match target_source_counts.get(&terminal.id).copied().unwrap_or(0) {
                0 => findings.push(MachineGraphFinding::new(
                    MachineGraphRule::MissingSourceClosure,
                    MachineGraphSubject::Terminal(terminal.id.clone()),
                    None,
                )),
                1 => {}
                _ => findings.push(MachineGraphFinding::new(
                    MachineGraphRule::MultipleInputSources,
                    MachineGraphSubject::Terminal(terminal.id.clone()),
                    None,
                )),
            }
        }
    }
    for state_slot in state_owners.keys() {
        match state_writer_counts.get(state_slot).copied().unwrap_or(0) {
            0 => findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnaccountedState,
                MachineGraphSubject::MachineElement(MachineElementId::StateSlot(
                    state_slot.clone(),
                )),
                None,
            )),
            1 => {}
            _ => findings.push(MachineGraphFinding::new(
                MachineGraphRule::MultipleStateWriters,
                MachineGraphSubject::MachineElement(MachineElementId::StateSlot(
                    state_slot.clone(),
                )),
                None,
            )),
        }
    }
    if let Some(relation) = first_algebraic_cycle_relation(&algebraic_adjacency) {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::AlgebraicLoopWithoutSolvePolicy,
            MachineGraphSubject::Relation(relation),
            None,
        ));
    }

    let mut material_targets = BTreeSet::<MaterialTarget>::new();
    for material in &draft.materials {
        if !material_targets.insert(material.target.clone()) {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateMaterialBinding,
                MachineGraphSubject::MaterialTarget(material.target.clone()),
                None,
            ));
        }
        let owned = match &material.target {
            MaterialTarget::Body(body) => {
                element_owners.contains_key(&MachineElementId::Body(body.clone()))
            }
            MaterialTarget::SurfacePatch(patch) => {
                element_owners.contains_key(&MachineElementId::SurfacePatch(patch.clone()))
            }
        };
        if !owned {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownMaterialTarget,
                MachineGraphSubject::MaterialTarget(material.target.clone()),
                None,
            ));
        }
    }
    for body in body_ids {
        let target = MaterialTarget::Body(body);
        if !material_targets.contains(&target) {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::MissingBodyMaterial,
                MachineGraphSubject::MaterialTarget(target),
                None,
            ));
        }
    }

    for pair in draft.interfaces.windows(2) {
        if pair[0].id == pair[1].id {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::DuplicateInterface,
                MachineGraphSubject::Interface(pair[1].id.clone()),
                None,
            ));
        }
    }
    let relation_pairs: BTreeSet<(TerminalId, TerminalId)> = draft
        .relations
        .iter()
        .map(|relation| (relation.source.clone(), relation.target.clone()))
        .collect();
    let mut port_interfaces = BTreeMap::<PortId, InterfaceId>::new();
    for interface in &draft.interfaces {
        if interface.negative == interface.positive {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::SameInterfaceEndpoint,
                MachineGraphSubject::Interface(interface.id.clone()),
                Some(MachineGraphSubject::Port(interface.negative.clone())),
            ));
        }
        for port_id in [&interface.negative, &interface.positive] {
            if let Some(first_interface) =
                port_interfaces.insert(port_id.clone(), interface.id.clone())
            {
                findings.push(MachineGraphFinding::new(
                    MachineGraphRule::PortInterfaceConflict,
                    MachineGraphSubject::Interface(interface.id.clone()),
                    Some(MachineGraphSubject::Interface(first_interface)),
                ));
            }
        }
        let negative = port_indices
            .get(&interface.negative)
            .map(|index| &draft.ports[*index]);
        let positive = port_indices
            .get(&interface.positive)
            .map(|index| &draft.ports[*index]);
        if negative.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownInterfacePort,
                MachineGraphSubject::Interface(interface.id.clone()),
                Some(MachineGraphSubject::Port(interface.negative.clone())),
            ));
        }
        if positive.is_none() {
            findings.push(MachineGraphFinding::new(
                MachineGraphRule::UnknownInterfacePort,
                MachineGraphSubject::Interface(interface.id.clone()),
                Some(MachineGraphSubject::Port(interface.positive.clone())),
            ));
        }
        if let (Some(negative), Some(positive)) = (negative, positive) {
            check_interface_compatibility(
                interface,
                negative,
                positive,
                &draft.terminals,
                &terminal_indices,
                &relation_pairs,
                &mut findings,
            );
        }
    }

    if !findings.is_empty() {
        return Err(graph_refusal(findings, None));
    }

    let receipt = match machine_graph_identity(&draft) {
        Ok(receipt) => receipt,
        Err(error) => {
            return Err(graph_refusal(
                vec![MachineGraphFinding::new(
                    MachineGraphRule::Identity,
                    MachineGraphSubject::Graph,
                    None,
                )],
                Some(error),
            ));
        }
    };
    Ok(AdmittedMachineGraph {
        clocks: draft.clocks,
        subsystems: draft.subsystems,
        terminals: draft.terminals,
        ports: draft.ports,
        relations: draft.relations,
        materials: draft.materials,
        interfaces: draft.interfaces,
        receipt,
    })
}

fn canonicalize_graph_draft(draft: &mut MachineGraphDraft) {
    draft.clocks.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| clock_row(left).cmp(&clock_row(right)))
    });
    for subsystem in &mut draft.subsystems {
        subsystem.bodies.sort();
        subsystem.surface_patches.sort();
        subsystem.contact_features.sort();
        subsystem.state_slots.sort();
    }
    draft.subsystems.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| subsystem_row(left).cmp(&subsystem_row(right)))
    });
    draft.terminals.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| terminal_row(left).cmp(&terminal_row(right)))
    });
    draft.ports.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| port_row(left).cmp(&port_row(right)))
    });
    draft.relations.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| relation_row(left).cmp(&relation_row(right)))
    });
    draft.materials.sort_by(|left, right| {
        left.target
            .cmp(&right.target)
            .then_with(|| left.material.cmp(&right.material))
    });
    draft.interfaces.sort_by(|left, right| {
        left.id
            .cmp(&right.id)
            .then_with(|| interface_row(left).cmp(&interface_row(right)))
    });
}

fn record_element_owner(
    owners: &mut BTreeMap<MachineElementId, SubsystemId>,
    element: MachineElementId,
    owner: &SubsystemId,
    findings: &mut Vec<MachineGraphFinding>,
) {
    if let Some(first_owner) = owners.insert(element.clone(), owner.clone()) {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::DuplicateOwnership,
            MachineGraphSubject::MachineElement(element),
            Some(MachineGraphSubject::Subsystem(first_owner)),
        ));
    }
}

fn first_algebraic_cycle_relation(adjacency: &[Vec<(usize, RelationId)>]) -> Option<RelationId> {
    fn visit(
        node: usize,
        adjacency: &[Vec<(usize, RelationId)>],
        colors: &mut [u8],
    ) -> Option<RelationId> {
        colors[node] = 1;
        for (next, relation) in &adjacency[node] {
            match colors[*next] {
                1 => return Some(relation.clone()),
                0 => {
                    if let Some(relation) = visit(*next, adjacency, colors) {
                        return Some(relation);
                    }
                }
                _ => {}
            }
        }
        colors[node] = 2;
        None
    }

    let mut colors = vec![0u8; adjacency.len()];
    for node in 0..adjacency.len() {
        if colors[node] == 0
            && let Some(relation) = visit(node, adjacency, &mut colors)
        {
            return Some(relation);
        }
    }
    None
}

fn declared_interface_relation_orientations(
    interfaces: &[InterfaceBinding],
    ports: &[PortSpec],
    port_indices: &BTreeMap<PortId, usize>,
    terminals: &[TerminalSpec],
    terminal_indices: &BTreeMap<TerminalId, usize>,
) -> BTreeMap<(TerminalId, TerminalId), InterfaceOrientation> {
    let mut relations = BTreeMap::new();
    for interface in interfaces {
        let Some(negative) = port_indices
            .get(&interface.negative)
            .map(|index| &ports[*index])
        else {
            continue;
        };
        let Some(positive) = port_indices
            .get(&interface.positive)
            .map(|index| &ports[*index])
        else {
            continue;
        };
        for (negative_terminal, positive_terminal) in [
            (&negative.effort, &positive.effort),
            (&negative.flow, &positive.flow),
        ] {
            let Some(negative_terminal) = terminal_indices
                .get(negative_terminal)
                .map(|index| &terminals[*index])
            else {
                continue;
            };
            let Some(positive_terminal) = terminal_indices
                .get(positive_terminal)
                .map(|index| &terminals[*index])
            else {
                continue;
            };
            if let Some(pair) = directed_terminal_pair(negative_terminal, positive_terminal) {
                relations.insert(pair, interface.orientation);
            }
        }
    }
    relations
}

fn directed_terminal_pair(
    first: &TerminalSpec,
    second: &TerminalSpec,
) -> Option<(TerminalId, TerminalId)> {
    match (first.causality, second.causality) {
        (TerminalCausality::Output, TerminalCausality::Input) => {
            Some((first.id.clone(), second.id.clone()))
        }
        (TerminalCausality::Input, TerminalCausality::Output) => {
            Some((second.id.clone(), first.id.clone()))
        }
        _ => None,
    }
}

// Interface checks deliberately keep all paired declarations in one view so
// a finding can name the interface instead of leaking an incidental loop order.
#[allow(clippy::too_many_lines)]
fn check_interface_compatibility(
    interface: &InterfaceBinding,
    negative: &PortSpec,
    positive: &PortSpec,
    terminals: &[TerminalSpec],
    terminal_indices: &BTreeMap<TerminalId, usize>,
    relation_pairs: &BTreeSet<(TerminalId, TerminalId)>,
    findings: &mut Vec<MachineGraphFinding>,
) {
    let Some(negative_effort) = terminal_indices
        .get(&negative.effort)
        .map(|index| &terminals[*index])
    else {
        return;
    };
    let Some(negative_flow) = terminal_indices
        .get(&negative.flow)
        .map(|index| &terminals[*index])
    else {
        return;
    };
    let Some(positive_effort) = terminal_indices
        .get(&positive.effort)
        .map(|index| &terminals[*index])
    else {
        return;
    };
    let Some(positive_flow) = terminal_indices
        .get(&positive.flow)
        .map(|index| &terminals[*index])
    else {
        return;
    };

    let subject = || MachineGraphSubject::Interface(interface.id.clone());
    if negative_effort.quantity != positive_effort.quantity
        || negative_flow.quantity != positive_flow.quantity
    {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceQuantityGap,
            subject(),
            None,
        ));
    }
    if negative_effort.shape != positive_effort.shape || negative_flow.shape != positive_flow.shape
    {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceShapeGap,
            subject(),
            None,
        ));
    }
    if negative_effort.clock != positive_effort.clock || negative_flow.clock != positive_flow.clock
    {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceClockGap,
            subject(),
            None,
        ));
    }
    if negative_effort.frame.canonical_key() != positive_effort.frame.canonical_key()
        || negative_flow.frame.canonical_key() != positive_flow.frame.canonical_key()
    {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceFrameGap,
            subject(),
            None,
        ));
    }
    let effort_orientations_match =
        negative_effort.frame.orientation() == positive_effort.frame.orientation();
    let flow_orientations_match =
        negative_flow.frame.orientation() == positive_flow.frame.orientation();
    let orientation_admitted = match interface.orientation {
        InterfaceOrientation::Aligned => effort_orientations_match && flow_orientations_match,
        InterfaceOrientation::Opposed => !effort_orientations_match && !flow_orientations_match,
    };
    if !orientation_admitted {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceOrientationGap,
            subject(),
            None,
        ));
    }
    if negative.energy_role == positive.energy_role {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceEnergyRoleGap,
            subject(),
            None,
        ));
    }
    let effort_relation = directed_terminal_pair(negative_effort, positive_effort);
    let flow_relation = directed_terminal_pair(negative_flow, positive_flow);
    if effort_relation.is_none() || flow_relation.is_none() {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceCausalityGap,
            subject(),
            None,
        ));
    }
    if effort_relation
        .as_ref()
        .is_some_and(|pair| !relation_pairs.contains(pair))
        || flow_relation
            .as_ref()
            .is_some_and(|pair| !relation_pairs.contains(pair))
    {
        findings.push(MachineGraphFinding::new(
            MachineGraphRule::InterfaceRelationGap,
            subject(),
            None,
        ));
    }
}

fn machine_graph_identity(
    draft: &MachineGraphDraft,
) -> Result<IdentityReceipt<MachineGraphIdV1>, CanonicalError> {
    let clock_rows: Vec<Vec<u8>> = draft.clocks.iter().map(clock_row).collect();
    let subsystem_rows: Vec<Vec<u8>> = draft.subsystems.iter().map(subsystem_row).collect();
    let terminal_rows: Vec<Vec<u8>> = draft.terminals.iter().map(terminal_row).collect();
    let port_rows: Vec<Vec<u8>> = draft.ports.iter().map(port_row).collect();
    let relation_rows: Vec<Vec<u8>> = draft.relations.iter().map(relation_row).collect();
    let material_rows: Vec<Vec<u8>> = draft.materials.iter().map(material_row).collect();
    let interface_rows: Vec<Vec<u8>> = draft.interfaces.iter().map(interface_row).collect();

    CanonicalEncoder::<MachineGraphIdV1, _>::new(MACHINE_GRAPH_IDENTITY_LIMITS, NeverCancel)?
        .u64(
            Field::new(0, "graph-schema-version"),
            u64::from(MACHINE_GRAPH_SCHEMA_VERSION_V1),
        )?
        .ordered_bytes(
            Field::new(1, "clocks"),
            clock_rows.len() as u64,
            clock_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(2, "subsystems"),
            subsystem_rows.len() as u64,
            subsystem_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(3, "terminals"),
            terminal_rows.len() as u64,
            terminal_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(4, "ports"),
            port_rows.len() as u64,
            port_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(5, "relations"),
            relation_rows.len() as u64,
            relation_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(6, "materials"),
            material_rows.len() as u64,
            material_rows.iter().map(Vec::as_slice),
        )?
        .ordered_bytes(
            Field::new(7, "interfaces"),
            interface_rows.len() as u64,
            interface_rows.iter().map(Vec::as_slice),
        )?
        .finish()
}

fn push_len_prefixed(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

fn push_identity(out: &mut Vec<u8>, identity: &[u8; 32]) {
    out.extend_from_slice(identity);
}

fn push_identity_collection<T>(
    out: &mut Vec<u8>,
    values: &[T],
    mut identity: impl FnMut(&T) -> [u8; 32],
) {
    out.extend_from_slice(&(values.len() as u64).to_le_bytes());
    for value in values {
        push_identity(out, &identity(value));
    }
}

fn clock_row(spec: &ClockSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(56);
    push_identity(&mut out, &spec.id.digest_bytes());
    match spec.clock {
        MachineClock::Continuous => out.push(1),
        MachineClock::Periodic {
            period_ns,
            phase_ns,
        } => {
            out.push(2);
            out.extend_from_slice(&period_ns.get().to_le_bytes());
            out.extend_from_slice(&phase_ns.to_le_bytes());
        }
        MachineClock::EventDriven => out.push(3),
    }
    out
}

fn subsystem_row(spec: &SubsystemSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(
        96 + 32
            * (spec.bodies.len()
                + spec.surface_patches.len()
                + spec.contact_features.len()
                + spec.state_slots.len()),
    );
    push_identity(&mut out, &spec.id.digest_bytes());
    spec.model.append_canonical(&mut out);
    push_identity_collection(&mut out, &spec.bodies, BodyId::digest_bytes);
    push_identity_collection(
        &mut out,
        &spec.surface_patches,
        SurfacePatchId::digest_bytes,
    );
    push_identity_collection(
        &mut out,
        &spec.contact_features,
        ContactFeatureId::digest_bytes,
    );
    push_identity_collection(&mut out, &spec.state_slots, StateSlotId::digest_bytes);
    out
}

fn terminal_row(spec: &TerminalSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_identity(&mut out, &spec.id.digest_bytes());
    push_identity(&mut out, &spec.owner.digest_bytes());
    push_terminal_quantity(&mut out, spec.quantity);
    push_terminal_shape(&mut out, spec.shape);
    out.push(terminal_causality_tag(spec.causality));
    push_identity(&mut out, &spec.clock.digest_bytes());
    push_len_prefixed(&mut out, spec.frame.canonical_key().as_bytes());
    out.push(orientation_parity_tag(spec.frame.orientation()));
    out
}

fn port_row(spec: &PortSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(129);
    push_identity(&mut out, &spec.id.digest_bytes());
    push_identity(&mut out, &spec.owner.digest_bytes());
    push_identity(&mut out, &spec.effort.digest_bytes());
    push_identity(&mut out, &spec.flow.digest_bytes());
    out.push(port_energy_role_tag(spec.energy_role));
    out
}

fn relation_row(spec: &RelationSpec) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_identity(&mut out, &spec.id.digest_bytes());
    push_identity(&mut out, &spec.source.digest_bytes());
    push_identity(&mut out, &spec.target.digest_bytes());
    match &spec.mode {
        RelationMode::Algebraic { solve_policy: None } => out.push(1),
        RelationMode::Algebraic {
            solve_policy: Some(policy),
        } => {
            out.push(2);
            policy.append_canonical(&mut out);
        }
        RelationMode::Stateful { state_slot } => {
            out.push(3);
            push_identity(&mut out, &state_slot.digest_bytes());
        }
    }
    out
}

fn material_row(binding: &MaterialBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(128);
    match &binding.target {
        MaterialTarget::Body(body) => {
            out.push(1);
            push_identity(&mut out, &body.digest_bytes());
        }
        MaterialTarget::SurfacePatch(patch) => {
            out.push(2);
            push_identity(&mut out, &patch.digest_bytes());
        }
    }
    binding.material.append_canonical(&mut out);
    out
}

fn interface_row(binding: &InterfaceBinding) -> Vec<u8> {
    let mut out = Vec::with_capacity(192);
    push_identity(&mut out, &binding.id.digest_bytes());
    push_identity(&mut out, &binding.negative.digest_bytes());
    push_identity(&mut out, &binding.positive.digest_bytes());
    binding.interface.append_canonical(&mut out);
    out.push(interface_orientation_tag(binding.orientation));
    out
}

fn push_terminal_quantity(out: &mut Vec<u8>, quantity: TerminalQuantitySpec) {
    match quantity {
        TerminalQuantitySpec::Dimensional(dims) => {
            out.push(1);
            push_dims(out, dims);
        }
        TerminalQuantitySpec::Semantic(semantic_type) => {
            out.push(2);
            push_quantity_kind(out, semantic_type.kind());
            out.push(value_form_tag(semantic_type.form()));
            push_dims(out, semantic_type.expected_dims());
        }
    }
}

fn push_dims(out: &mut Vec<u8>, dims: Dims) {
    out.extend(dims.0.map(|exponent| exponent as u8));
}

fn push_quantity_kind(out: &mut Vec<u8>, kind: QuantityKind) {
    match kind {
        QuantityKind::AbsoluteTemperature => out.push(1),
        QuantityKind::TemperatureDifference => out.push(2),
        QuantityKind::Angle(domain) => {
            out.push(3);
            out.push(angle_domain_tag(domain));
        }
        QuantityKind::AngularVelocity(domain) => {
            out.push(4);
            out.push(angle_domain_tag(domain));
        }
        QuantityKind::Torque => out.push(5),
        QuantityKind::Energy => out.push(6),
        QuantityKind::Pressure => out.push(7),
        QuantityKind::Stress => out.push(8),
        QuantityKind::Strain { basis, component } => {
            out.push(9);
            out.push(strain_basis_tag(basis));
            out.push(strain_component_tag(component));
        }
        QuantityKind::Composition(basis) => {
            out.push(10);
            out.push(composition_basis_tag(basis));
        }
        QuantityKind::Mass => out.push(11),
        QuantityKind::Amount => out.push(12),
        QuantityKind::MolarMass => out.push(13),
        QuantityKind::MassConcentration => out.push(14),
        QuantityKind::AmountConcentration => out.push(15),
        QuantityKind::Entropy => out.push(16),
        QuantityKind::HeatCapacity => out.push(17),
        QuantityKind::AcousticPressure => out.push(18),
        QuantityKind::AcousticPower => out.push(19),
    }
}

fn push_terminal_shape(out: &mut Vec<u8>, shape: TerminalShape) {
    match shape {
        TerminalShape::Scalar => out.push(1),
        TerminalShape::Vector { components } => {
            out.push(2);
            out.extend_from_slice(&components.get().to_le_bytes());
        }
        TerminalShape::Tensor { rows, columns } => {
            out.push(3);
            out.extend_from_slice(&rows.get().to_le_bytes());
            out.extend_from_slice(&columns.get().to_le_bytes());
        }
        TerminalShape::FieldTrace { components } => {
            out.push(4);
            out.extend_from_slice(&components.get().to_le_bytes());
        }
    }
}

const fn angle_domain_tag(domain: AngleDomain) -> u8 {
    match domain {
        AngleDomain::Mechanical => 1,
        AngleDomain::Electrical => 2,
    }
}

const fn strain_basis_tag(basis: StrainBasis) -> u8 {
    match basis {
        StrainBasis::Tensor => 1,
        StrainBasis::Engineering => 2,
    }
}

const fn strain_component_tag(component: StrainComponent) -> u8 {
    match component {
        StrainComponent::Normal => 1,
        StrainComponent::Shear => 2,
    }
}

const fn composition_basis_tag(basis: CompositionBasis) -> u8 {
    match basis {
        CompositionBasis::MassFraction => 1,
        CompositionBasis::MoleFraction => 2,
        CompositionBasis::VolumeFraction => 3,
    }
}

const fn value_form_tag(form: ValueForm) -> u8 {
    match form {
        ValueForm::Static => 1,
        ValueForm::Instantaneous => 2,
        ValueForm::Peak => 3,
        ValueForm::Rms => 4,
    }
}

const fn terminal_causality_tag(causality: TerminalCausality) -> u8 {
    match causality {
        TerminalCausality::Input => 1,
        TerminalCausality::Output => 2,
        TerminalCausality::ExternalInput => 3,
    }
}

const fn orientation_parity_tag(orientation: OrientationParity) -> u8 {
    match orientation {
        OrientationParity::Preserving => 1,
        OrientationParity::Reversing => 2,
    }
}

const fn port_energy_role_tag(role: PortEnergyRole) -> u8 {
    match role {
        PortEnergyRole::IntoSubsystem => 1,
        PortEnergyRole::OutOfSubsystem => 2,
    }
}

const fn interface_orientation_tag(orientation: InterfaceOrientation) -> u8 {
    match orientation {
        InterfaceOrientation::Aligned => 1,
        InterfaceOrientation::Opposed => 2,
    }
}
