//! Run-scoped operation memory lease (bead wf9.16).
//!
//! One operation — a TilePool run and everything it allocates — obtains a
//! single [`OperationMemoryLease`] and every constituent mechanism charges
//! it: executor root metadata before worker launch, and every tile arena's
//! chunks while the operation holds them. The lease and [`ArenaPool`]'s
//! process-wide `limit_bytes` are DIFFERENT ledgers with different
//! lifetimes: the pool counts OS-reserved bytes (in-use + free-listed,
//! across operations), the lease counts one operation's live set. A chunk
//! recycled from the pool free list charges the acquiring operation's lease
//! exactly while held and never twice; free-list inventory belongs to no
//! operation. Both gates must admit; a refusal names whichever refused.
//!
//! Receipts have a canonical structure and exact values for the observed
//! admission trace. Identical plans with identical pool-cache state have
//! deterministic cumulative demand; cache history and near-limit concurrent
//! refusals/peaks are intentionally visible and can change the receipt.
//! Thread stacks and allocator overhead are explicitly NOT claimed (CONTRACT
//! no-claims).
//!
//! [`ArenaPool`]: crate::ArenaPool

use core::fmt;
use core::marker::PhantomData;
use std::sync::{Arc, Mutex};

/// One refused lease reservation (the FIRST refusal is retained verbatim in
/// the receipt; later refusals only count).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseRefusal {
    /// Component that requested the bytes (e.g. `"tilepool-root-metadata"`,
    /// `"arena-chunk"`).
    pub what: &'static str,
    /// Bytes the component asked for.
    pub requested_bytes: u64,
    /// Lease bytes in use at refusal time.
    pub used_bytes: u64,
    /// The lease limit in force.
    pub limit_bytes: u64,
    reason: LeaseRefusalReason,
    sequence: u128,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LeaseRefusalReason {
    Capacity,
    Sealed,
    CounterOverflow,
}

impl LeaseRefusal {
    /// Stable reason captured in the same serialized transition as the
    /// refusal.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            LeaseRefusalReason::Capacity => "capacity",
            LeaseRefusalReason::Sealed => "sealed",
            LeaseRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Root-ledger sequence at which this refusal linearized.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for LeaseRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "operation memory lease refused {} B for `{}` with {} B of the {} B lease in use",
            self.requested_bytes, self.what, self.used_bytes, self.limit_bytes
        )
    }
}

/// Canonically serialized lease accounting snapshot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseReceipt {
    /// The limit in force (`None` = unbounded legacy wrapper).
    pub limit_bytes: Option<u64>,
    /// Exact cumulative bytes of granted reservations.
    pub requested_bytes: u128,
    /// Conservative logical high-water of concurrently held bytes.
    pub peak_bytes: u64,
    /// Bytes still held when the snapshot was taken.
    pub used_bytes: u64,
    /// Exact number of refused reservations.
    pub refusals: u128,
    /// Internal release attempts that did not match a live reservation.
    /// The counter remains fail-closed (used bytes are not changed) and makes
    /// an invariant violation visible without panicking from `Drop`.
    pub release_invariant_violations: u128,
    /// The first refusal, verbatim.
    pub first_refusal: Option<LeaseRefusal>,
}

impl LeaseReceipt {
    /// Canonical JSON object (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-operation-lease-v2\"");
        match self.limit_bytes {
            Some(limit) => {
                let _ = write!(out, ",\"limit_bytes\":{limit}");
            }
            None => out.push_str(",\"limit_bytes\":null"),
        }
        let _ = write!(
            out,
            ",\"requested_bytes\":{},\"peak_bytes\":{},\"used_bytes\":{},\"refusals\":{},\"release_invariant_violations\":{}",
            self.requested_bytes,
            self.peak_bytes,
            self.used_bytes,
            self.refusals,
            self.release_invariant_violations
        );
        match &self.first_refusal {
            Some(refusal) => {
                let what = json_escape(refusal.what);
                let _ = write!(
                    out,
                    ",\"first_refusal\":{{\"what\":\"{}\",\"requested_bytes\":{},\"used_bytes\":{},\"limit_bytes\":{},\"reason\":\"{}\",\"sequence\":{}}}",
                    what,
                    refusal.requested_bytes,
                    refusal.used_bytes,
                    refusal.limit_bytes,
                    refusal.reason(),
                    refusal.sequence()
                );
            }
            None => out.push_str(",\"first_refusal\":null"),
        }
        out.push('}');
        out
    }
}

const VERIFIED_RECEIPT_SCHEMA_VERSION: u16 = 2;
const MAX_DELEGATION_RECORDS: usize = 4096;
const MAX_LOGICAL_ID_BYTES: usize = 256;
const LEASE_IDENTITY_SCHEMA_VERSION: u16 = 1;
const LEASE_IDENTITY_DOMAIN_BYTES: usize = 8;
const LEASE_IDENTITY_SUBJECT_BYTES: usize = 32;
const LEASE_IDENTITY_MAX_PATH_COMPONENTS: usize = 16;
const LEASE_IDENTITY_ENCODED_BYTES: usize = 2
    + LEASE_IDENTITY_DOMAIN_BYTES
    + LEASE_IDENTITY_SUBJECT_BYTES
    + LEASE_IDENTITY_SUBJECT_BYTES
    + LEASE_IDENTITY_SUBJECT_BYTES
    + 1
    + LEASE_IDENTITY_MAX_PATH_COMPONENTS * size_of::<u64>();
const PUBLISHED_TRANSFER_BINDING_SCHEMA_VERSION: u16 = 1;
const PUBLISHED_TRANSFER_BINDING_FIELD_BYTES: usize = 32;
const PUBLISHED_TRANSFER_BINDING_ENCODED_BYTES: usize =
    2 + PUBLISHED_TRANSFER_BINDING_FIELD_BYTES * 4;
const PUBLISHED_TRANSFER_ENVELOPE_SCHEMA_VERSION: u16 = 1;
const PUBLISHED_TRANSFER_ENVELOPE_ENCODED_BYTES: usize = 2 + size_of::<u64>() * 3;
const PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION: u16 = 1;

/// Allocation-free, versioned authority identity for one root or delegated
/// memory envelope.
///
/// `domain` separates protocols, `root_subject` binds every descendant to
/// the same invocation (or other root authority), `parent_subject` proves the
/// immediate transferor, `owner_subject` binds the current envelope to its
/// exact child authority, and `path` records deterministic subdivision
/// ordinals. The fixed-size canonical encoding can be hashed or logged without
/// formatting a dynamic identifier or allocating.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LeaseIdentity {
    domain: [u8; LEASE_IDENTITY_DOMAIN_BYTES],
    root_subject: [u8; LEASE_IDENTITY_SUBJECT_BYTES],
    parent_subject: [u8; LEASE_IDENTITY_SUBJECT_BYTES],
    owner_subject: [u8; LEASE_IDENTITY_SUBJECT_BYTES],
    path_len: u8,
    path: [u64; LEASE_IDENTITY_MAX_PATH_COMPONENTS],
}

impl LeaseIdentity {
    /// Canonical encoding schema version.
    pub const SCHEMA_VERSION: u16 = LEASE_IDENTITY_SCHEMA_VERSION;
    /// Required fixed domain-tag width.
    pub const DOMAIN_BYTES: usize = LEASE_IDENTITY_DOMAIN_BYTES;
    /// Required root/owner subject width.
    pub const SUBJECT_BYTES: usize = LEASE_IDENTITY_SUBJECT_BYTES;
    /// Maximum deterministic delegation depth.
    pub const MAX_PATH_COMPONENTS: usize = LEASE_IDENTITY_MAX_PATH_COMPONENTS;
    /// Fixed canonical encoding width.
    pub const ENCODED_BYTES: usize = LEASE_IDENTITY_ENCODED_BYTES;

    /// Construct a root identity. The root owns itself and has an empty path.
    #[must_use]
    pub const fn root(
        domain: [u8; LEASE_IDENTITY_DOMAIN_BYTES],
        root_subject: [u8; LEASE_IDENTITY_SUBJECT_BYTES],
    ) -> Self {
        Self {
            domain,
            root_subject,
            parent_subject: root_subject,
            owner_subject: root_subject,
            path_len: 0,
            path: [0; LEASE_IDENTITY_MAX_PATH_COMPONENTS],
        }
    }

    /// Construct one direct child by appending a deterministic path component
    /// and binding the child authority's exact subject.
    ///
    /// # Errors
    ///
    /// Returns [`LeaseIdentityPathError`] when the fixed path capacity is
    /// exhausted. No mutation or allocation occurs on refusal.
    pub fn child(
        self,
        owner_subject: [u8; LEASE_IDENTITY_SUBJECT_BYTES],
        path_component: u64,
    ) -> Result<Self, LeaseIdentityPathError> {
        let index = usize::from(self.path_len);
        if index >= LEASE_IDENTITY_MAX_PATH_COMPONENTS {
            return Err(LeaseIdentityPathError {
                maximum_components: LEASE_IDENTITY_MAX_PATH_COMPONENTS,
            });
        }
        let mut child = self;
        child.parent_subject = self.owner_subject;
        child.owner_subject = owner_subject;
        child.path[index] = path_component;
        child.path_len += 1;
        Ok(child)
    }

    /// Fixed protocol-domain tag.
    #[must_use]
    pub const fn domain(self) -> [u8; LEASE_IDENTITY_DOMAIN_BYTES] {
        self.domain
    }

    /// Root authority subject shared by the complete delegation tree.
    #[must_use]
    pub const fn root_subject(self) -> [u8; LEASE_IDENTITY_SUBJECT_BYTES] {
        self.root_subject
    }

    /// Exact subject that transferred this envelope. A root names itself.
    #[must_use]
    pub const fn parent_subject(self) -> [u8; LEASE_IDENTITY_SUBJECT_BYTES] {
        self.parent_subject
    }

    /// Exact subject that owns this envelope.
    #[must_use]
    pub const fn owner_subject(self) -> [u8; LEASE_IDENTITY_SUBJECT_BYTES] {
        self.owner_subject
    }

    /// Deterministic subdivision path.
    #[must_use]
    pub fn path(&self) -> &[u64] {
        &self.path[..usize::from(self.path_len)]
    }

    /// Versioned fixed-size canonical binary encoding.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; LEASE_IDENTITY_ENCODED_BYTES] {
        let mut encoded = [0_u8; LEASE_IDENTITY_ENCODED_BYTES];
        let mut cursor = 0;
        encoded[cursor..cursor + 2].copy_from_slice(&LEASE_IDENTITY_SCHEMA_VERSION.to_le_bytes());
        cursor += 2;
        encoded[cursor..cursor + LEASE_IDENTITY_DOMAIN_BYTES].copy_from_slice(&self.domain);
        cursor += LEASE_IDENTITY_DOMAIN_BYTES;
        encoded[cursor..cursor + LEASE_IDENTITY_SUBJECT_BYTES].copy_from_slice(&self.root_subject);
        cursor += LEASE_IDENTITY_SUBJECT_BYTES;
        encoded[cursor..cursor + LEASE_IDENTITY_SUBJECT_BYTES]
            .copy_from_slice(&self.parent_subject);
        cursor += LEASE_IDENTITY_SUBJECT_BYTES;
        encoded[cursor..cursor + LEASE_IDENTITY_SUBJECT_BYTES].copy_from_slice(&self.owner_subject);
        cursor += LEASE_IDENTITY_SUBJECT_BYTES;
        encoded[cursor] = self.path_len;
        cursor += 1;
        for component in self.path {
            encoded[cursor..cursor + size_of::<u64>()].copy_from_slice(&component.to_le_bytes());
            cursor += size_of::<u64>();
        }
        encoded
    }

    /// Canonical JSON object with fixed field order and no display-derived
    /// identity.
    #[must_use]
    pub fn to_json(self) -> String {
        use fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-lease-identity-v1\"");
        let _ = write!(
            out,
            ",\"schema_version\":{},\"domain\":\"{}\",\"root_subject\":\"{}\",\"parent_subject\":\"{}\",\"owner_subject\":\"{}\",\"path\":[",
            LEASE_IDENTITY_SCHEMA_VERSION,
            byte_slice_hex(&self.domain),
            byte_slice_hex(&self.root_subject),
            byte_slice_hex(&self.parent_subject),
            byte_slice_hex(&self.owner_subject)
        );
        for (index, component) in self.path().iter().enumerate() {
            if index != 0 {
                out.push(',');
            }
            let _ = write!(out, "{component}");
        }
        out.push_str("]}");
        out
    }

    fn is_root(self) -> bool {
        self.path_len == 0
            && self.parent_subject == self.root_subject
            && self.owner_subject == self.root_subject
    }

    fn is_direct_child_of(self, parent: Self) -> bool {
        if self.domain != parent.domain
            || self.root_subject != parent.root_subject
            || self.parent_subject != parent.owner_subject
            || usize::from(self.path_len) != usize::from(parent.path_len) + 1
        {
            return false;
        }
        let parent_len = usize::from(parent.path_len);
        self.path[..parent_len] == parent.path[..parent_len]
    }
}

impl fmt::Debug for LeaseIdentity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LeaseIdentity")
            .field("domain", &byte_slice_hex(&self.domain))
            .field("root_subject", &byte_slice_hex(&self.root_subject))
            .field("parent_subject", &byte_slice_hex(&self.parent_subject))
            .field("owner_subject", &byte_slice_hex(&self.owner_subject))
            .field("path_len", &self.path_len)
            .field("path", &self.path())
            .finish()
    }
}

/// Fixed-path exhaustion while deriving a child identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LeaseIdentityPathError {
    maximum_components: usize,
}

impl LeaseIdentityPathError {
    /// Maximum supported deterministic path depth.
    #[must_use]
    pub const fn maximum_components(self) -> usize {
        self.maximum_components
    }
}

impl fmt::Display for LeaseIdentityPathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "lease identity path exceeds {} components",
            self.maximum_components
        )
    }
}

impl std::error::Error for LeaseIdentityPathError {}

/// Fixed-width identity of one planned publication into one destination.
///
/// The four fields deliberately remain distinct. A plan may execute the same
/// occurrence more than once, an occurrence may produce multiple outputs, and
/// one output identity may be routed to different destinations. Treating the
/// tuple as a single caller-formatted label would lose those authority
/// boundaries and make duplicate detection ambiguous.
#[allow(clippy::struct_field_names)] // the four *_identity fields ARE the point: distinct authority boundaries, per the doc comment above
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PublishedTransferBinding {
    plan_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
    occurrence_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
    output_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
    destination_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
}

impl PublishedTransferBinding {
    /// Canonical encoding schema version.
    pub const SCHEMA_VERSION: u16 = PUBLISHED_TRANSFER_BINDING_SCHEMA_VERSION;
    /// Width of each identity field.
    pub const FIELD_BYTES: usize = PUBLISHED_TRANSFER_BINDING_FIELD_BYTES;
    /// Width of the fixed canonical encoding.
    pub const ENCODED_BYTES: usize = PUBLISHED_TRANSFER_BINDING_ENCODED_BYTES;

    /// Construct one exact publication identity tuple.
    #[must_use]
    pub const fn new(
        plan_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
        occurrence_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
        output_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
        destination_identity: [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES],
    ) -> Self {
        Self {
            plan_identity,
            occurrence_identity,
            output_identity,
            destination_identity,
        }
    }

    /// Identity of the plan that authorized publication.
    #[must_use]
    pub const fn plan_identity(self) -> [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES] {
        self.plan_identity
    }

    /// Identity of the exact plan occurrence.
    #[must_use]
    pub const fn occurrence_identity(self) -> [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES] {
        self.occurrence_identity
    }

    /// Identity of the exact produced output.
    #[must_use]
    pub const fn output_identity(self) -> [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES] {
        self.output_identity
    }

    /// Identity of the destination that accepts ownership.
    #[must_use]
    pub const fn destination_identity(self) -> [u8; PUBLISHED_TRANSFER_BINDING_FIELD_BYTES] {
        self.destination_identity
    }

    /// Fixed, versioned canonical binary encoding.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; PUBLISHED_TRANSFER_BINDING_ENCODED_BYTES] {
        let mut encoded = [0_u8; PUBLISHED_TRANSFER_BINDING_ENCODED_BYTES];
        let mut cursor = 0;
        encoded[cursor..cursor + 2]
            .copy_from_slice(&PUBLISHED_TRANSFER_BINDING_SCHEMA_VERSION.to_le_bytes());
        cursor += 2;
        for field in [
            self.plan_identity,
            self.occurrence_identity,
            self.output_identity,
            self.destination_identity,
        ] {
            encoded[cursor..cursor + PUBLISHED_TRANSFER_BINDING_FIELD_BYTES]
                .copy_from_slice(&field);
            cursor += PUBLISHED_TRANSFER_BINDING_FIELD_BYTES;
        }
        encoded
    }

    /// Canonical JSON object with fixed field order.
    #[must_use]
    pub fn to_json(self) -> String {
        format!(
            "{{\"schema\":\"fs-alloc-published-transfer-binding-v1\",\"schema_version\":{},\"plan_identity\":\"{}\",\"occurrence_identity\":\"{}\",\"output_identity\":\"{}\",\"destination_identity\":\"{}\"}}",
            PUBLISHED_TRANSFER_BINDING_SCHEMA_VERSION,
            byte_slice_hex(&self.plan_identity),
            byte_slice_hex(&self.occurrence_identity),
            byte_slice_hex(&self.output_identity),
            byte_slice_hex(&self.destination_identity)
        )
    }
}

/// Exact byte composition transferred with a published output.
///
/// Payload, layout/indexing storage, and ownership overhead stay separate in
/// evidence even though admission enforces their checked total against one
/// already-live staging charge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[allow(clippy::struct_field_names)] // payload/layout/overhead stay separate byte accounts in evidence, per the doc comment above
pub struct PublishedTransferEnvelope {
    payload_bytes: u64,
    layout_bytes: u64,
    overhead_bytes: u64,
}

impl PublishedTransferEnvelope {
    /// Canonical encoding schema version.
    pub const SCHEMA_VERSION: u16 = PUBLISHED_TRANSFER_ENVELOPE_SCHEMA_VERSION;
    /// Width of the fixed canonical encoding.
    pub const ENCODED_BYTES: usize = PUBLISHED_TRANSFER_ENVELOPE_ENCODED_BYTES;

    /// Construct an exact three-part output envelope.
    #[must_use]
    pub const fn new(payload_bytes: u64, layout_bytes: u64, overhead_bytes: u64) -> Self {
        Self {
            payload_bytes,
            layout_bytes,
            overhead_bytes,
        }
    }

    /// Construct an envelope whose complete charge is payload.
    #[must_use]
    pub const fn payload_only(payload_bytes: u64) -> Self {
        Self::new(payload_bytes, 0, 0)
    }

    /// Payload bytes.
    #[must_use]
    pub const fn payload_bytes(self) -> u64 {
        self.payload_bytes
    }

    /// Layout/indexing bytes.
    #[must_use]
    pub const fn layout_bytes(self) -> u64 {
        self.layout_bytes
    }

    /// Ownership overhead bytes.
    #[must_use]
    pub const fn overhead_bytes(self) -> u64 {
        self.overhead_bytes
    }

    /// Checked total bytes represented by the envelope.
    #[must_use]
    pub const fn total_bytes(self) -> Option<u64> {
        let Some(payload_and_layout) = self.payload_bytes.checked_add(self.layout_bytes) else {
            return None;
        };
        payload_and_layout.checked_add(self.overhead_bytes)
    }

    /// Fixed, versioned canonical binary encoding.
    #[must_use]
    pub fn canonical_bytes(self) -> [u8; PUBLISHED_TRANSFER_ENVELOPE_ENCODED_BYTES] {
        let mut encoded = [0_u8; PUBLISHED_TRANSFER_ENVELOPE_ENCODED_BYTES];
        encoded[..2].copy_from_slice(&PUBLISHED_TRANSFER_ENVELOPE_SCHEMA_VERSION.to_le_bytes());
        let mut cursor = 2;
        for value in [self.payload_bytes, self.layout_bytes, self.overhead_bytes] {
            encoded[cursor..cursor + size_of::<u64>()].copy_from_slice(&value.to_le_bytes());
            cursor += size_of::<u64>();
        }
        encoded
    }

    /// Canonical JSON object with fixed field order.
    #[must_use]
    pub fn to_json(self) -> String {
        format!(
            "{{\"schema\":\"fs-alloc-published-transfer-envelope-v1\",\"schema_version\":{},\"payload_bytes\":{},\"layout_bytes\":{},\"overhead_bytes\":{},\"total_bytes\":{}}}",
            PUBLISHED_TRANSFER_ENVELOPE_SCHEMA_VERSION,
            self.payload_bytes,
            self.layout_bytes,
            self.overhead_bytes,
            self.total_bytes()
                .map_or_else(|| String::from("null"), |total| total.to_string())
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfigurationRefusalReason {
    UnboundedRoot,
    NotPristine,
    AlreadyConfigured,
    InvalidRootIdentity,
    MetadataLimit,
    MetadataAllocation,
    Sealed,
}

/// Refusal to enable verified delegation on an existing root lease.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseConfigurationRefusal {
    root_identity: LeaseIdentity,
    root_id: &'static str,
    reason: ConfigurationRefusalReason,
    sequence: u128,
}

impl LeaseConfigurationRefusal {
    /// Typed root identity requested by the caller.
    #[must_use]
    pub fn root_identity(&self) -> LeaseIdentity {
        self.root_identity
    }

    /// Static diagnostic label requested by the caller.
    #[must_use]
    pub fn root_id(&self) -> &'static str {
        self.root_id
    }

    /// Stable refusal code captured while the root-state mutex was held.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            ConfigurationRefusalReason::UnboundedRoot => "unbounded_root",
            ConfigurationRefusalReason::NotPristine => "root_not_pristine",
            ConfigurationRefusalReason::AlreadyConfigured => "already_configured",
            ConfigurationRefusalReason::InvalidRootIdentity => "invalid_root_identity",
            ConfigurationRefusalReason::MetadataLimit => "metadata_limit",
            ConfigurationRefusalReason::MetadataAllocation => "metadata_allocation",
            ConfigurationRefusalReason::Sealed => "sealed",
        }
    }

    /// Root-ledger sequence at which the refusal linearized.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for LeaseConfigurationRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "operation memory lease delegation configuration for `{}` refused: {} at sequence {}",
            self.root_id,
            self.reason(),
            self.sequence
        )
    }
}

impl std::error::Error for LeaseConfigurationRefusal {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DelegationRefusalReason {
    UnboundedParent,
    UnconfiguredRoot,
    RootSealed,
    InvalidIdentityRelationship,
    InvalidLogicalPath,
    DuplicateIdentity,
    DuplicatePath,
    MetadataExhausted,
    Capacity,
    ParentReturned,
    CounterOverflow,
}

/// Refusal to transfer capacity to a logical child.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseDelegationRefusal {
    root_identity: Option<LeaseIdentity>,
    parent_identity: Option<LeaseIdentity>,
    identity: LeaseIdentity,
    root_id: Option<&'static str>,
    parent_path: Option<&'static str>,
    logical_path: &'static str,
    requested_bytes: u64,
    used_bytes: u64,
    limit_bytes: Option<u64>,
    reason: DelegationRefusalReason,
    sequence: u128,
}

impl LeaseDelegationRefusal {
    /// Configured typed root identity, when one exists.
    #[must_use]
    pub fn root_identity(&self) -> Option<LeaseIdentity> {
        self.root_identity
    }

    /// Typed parent identity; `None` denotes the root.
    #[must_use]
    pub fn parent_identity(&self) -> Option<LeaseIdentity> {
        self.parent_identity
    }

    /// Exact typed identity requested for the child.
    #[must_use]
    pub fn identity(&self) -> LeaseIdentity {
        self.identity
    }

    /// Configured root identity, when one exists.
    #[must_use]
    pub fn root_id(&self) -> Option<&'static str> {
        self.root_id
    }

    /// Logical parent path; `None` denotes the root.
    #[must_use]
    pub fn parent_path(&self) -> Option<&'static str> {
        self.parent_path
    }

    /// Full caller-supplied child path.
    #[must_use]
    pub fn logical_path(&self) -> &'static str {
        self.logical_path
    }

    /// Child capacity requested in bytes.
    #[must_use]
    pub fn requested_bytes(&self) -> u64 {
        self.requested_bytes
    }

    /// Parent bytes occupied at the serialized refusal point.
    #[must_use]
    pub fn used_bytes(&self) -> u64 {
        self.used_bytes
    }

    /// Parent capacity, or `None` for an unbounded compatibility root.
    #[must_use]
    pub fn limit_bytes(&self) -> Option<u64> {
        self.limit_bytes
    }

    /// Capacity still available at the serialized refusal point.
    #[must_use]
    pub fn available_bytes(&self) -> Option<u64> {
        self.limit_bytes
            .and_then(|limit| limit.checked_sub(self.used_bytes))
    }

    /// Stable refusal code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            DelegationRefusalReason::UnboundedParent => "unbounded_parent",
            DelegationRefusalReason::UnconfiguredRoot => "unconfigured_root",
            DelegationRefusalReason::RootSealed => "root_sealed",
            DelegationRefusalReason::InvalidIdentityRelationship => "invalid_identity_relationship",
            DelegationRefusalReason::InvalidLogicalPath => "invalid_logical_path",
            DelegationRefusalReason::DuplicateIdentity => "duplicate_identity",
            DelegationRefusalReason::DuplicatePath => "duplicate_path",
            DelegationRefusalReason::MetadataExhausted => "metadata_exhausted",
            DelegationRefusalReason::Capacity => "capacity",
            DelegationRefusalReason::ParentReturned => "parent_returned",
            DelegationRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Root-ledger sequence at which the refusal linearized.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for LeaseDelegationRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "operation memory lease refused {} B transfer to `{}`: {} at sequence {}",
            self.requested_bytes,
            self.logical_path,
            self.reason(),
            self.sequence
        )
    }
}

impl std::error::Error for LeaseDelegationRefusal {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DelegatedReservationRefusalReason {
    RootSealed,
    ChildReturned,
    Capacity,
    CounterOverflow,
}

/// Refusal of an allocation reservation inside a delegated envelope.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DelegatedLeaseRefusal {
    root_identity: LeaseIdentity,
    identity: LeaseIdentity,
    root_id: &'static str,
    logical_path: &'static str,
    site: &'static str,
    requested_bytes: u64,
    used_bytes: u64,
    limit_bytes: u64,
    reason: DelegatedReservationRefusalReason,
    sequence: u128,
}

impl DelegatedLeaseRefusal {
    /// Typed root authority identity.
    #[must_use]
    pub fn root_identity(&self) -> LeaseIdentity {
        self.root_identity
    }

    /// Exact typed delegated identity.
    #[must_use]
    pub fn identity(&self) -> LeaseIdentity {
        self.identity
    }

    /// Stable refusal code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            DelegatedReservationRefusalReason::RootSealed => "root_sealed",
            DelegatedReservationRefusalReason::ChildReturned => "child_returned",
            DelegatedReservationRefusalReason::Capacity => "capacity",
            DelegatedReservationRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Root identity.
    #[must_use]
    pub fn root_id(&self) -> &'static str {
        self.root_id
    }

    /// Exact delegated logical path.
    #[must_use]
    pub fn logical_path(&self) -> &'static str {
        self.logical_path
    }

    /// Allocation site.
    #[must_use]
    pub fn site(&self) -> &'static str {
        self.site
    }

    /// Requested bytes.
    #[must_use]
    pub fn requested_bytes(&self) -> u64 {
        self.requested_bytes
    }

    /// Child bytes occupied at refusal.
    #[must_use]
    pub fn used_bytes(&self) -> u64 {
        self.used_bytes
    }

    /// Child capacity.
    #[must_use]
    pub fn limit_bytes(&self) -> u64 {
        self.limit_bytes
    }

    /// Capacity still available at the serialized refusal point.
    #[must_use]
    pub fn available_bytes(&self) -> Option<u64> {
        self.limit_bytes.checked_sub(self.used_bytes)
    }

    /// Serialized root-ledger sequence.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for DelegatedLeaseRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "delegated memory lease `{}` refused {} B for `{}`: {} at sequence {}",
            self.logical_path,
            self.requested_bytes,
            self.site,
            self.reason(),
            self.sequence
        )
    }
}

impl std::error::Error for DelegatedLeaseRefusal {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PublishedTransferOperation {
    Prepare,
    Publish,
    Rollback,
    CloseDestination,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PublishedTransferRefusalReason {
    RootSealed,
    ChildReturned,
    ZeroBytes,
    EnvelopeMismatch,
    DuplicateBinding,
    MetadataExhausted,
    TransferUnavailable,
    ConservationMismatch,
    CounterOverflow,
}

/// Structured refusal for one publication ownership transition.
///
/// The diagnostic captures the complete fixed-width publication binding and
/// typed lease authority at the same mutex-serialized sequence as the failed
/// transition. No dynamic diagnostic allocation occurs on the refusal path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublishedTransferRefusal {
    root_identity: LeaseIdentity,
    parent_identity: Option<LeaseIdentity>,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    operation: PublishedTransferOperation,
    reason: PublishedTransferRefusalReason,
    sequence: u128,
}

impl PublishedTransferRefusal {
    /// Typed root authority identity.
    #[must_use]
    pub fn root_identity(&self) -> LeaseIdentity {
        self.root_identity
    }

    /// Typed delegated parent, or `None` for a direct root child.
    #[must_use]
    pub fn parent_identity(&self) -> Option<LeaseIdentity> {
        self.parent_identity
    }

    /// Exact child that owned the staging allocation.
    #[must_use]
    pub fn child_identity(&self) -> LeaseIdentity {
        self.child_identity
    }

    /// Exact plan/occurrence/output/destination tuple.
    #[must_use]
    pub fn binding(&self) -> PublishedTransferBinding {
        self.binding
    }

    /// Exact payload/layout/overhead byte composition.
    #[must_use]
    pub fn envelope(&self) -> PublishedTransferEnvelope {
        self.envelope
    }

    /// Bytes affected by the refused transition.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Stable transition name.
    #[must_use]
    pub fn operation(&self) -> &'static str {
        match self.operation {
            PublishedTransferOperation::Prepare => "prepare",
            PublishedTransferOperation::Publish => "publish",
            PublishedTransferOperation::Rollback => "rollback",
            PublishedTransferOperation::CloseDestination => "close_destination",
        }
    }

    /// Stable refusal code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            PublishedTransferRefusalReason::RootSealed => "root_sealed",
            PublishedTransferRefusalReason::ChildReturned => "child_returned",
            PublishedTransferRefusalReason::ZeroBytes => "zero_bytes",
            PublishedTransferRefusalReason::EnvelopeMismatch => "envelope_mismatch",
            PublishedTransferRefusalReason::DuplicateBinding => "duplicate_binding",
            PublishedTransferRefusalReason::MetadataExhausted => "metadata_exhausted",
            PublishedTransferRefusalReason::TransferUnavailable => "transfer_unavailable",
            PublishedTransferRefusalReason::ConservationMismatch => "conservation_mismatch",
            PublishedTransferRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Serialized root-ledger sequence.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for PublishedTransferRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "published transfer {} of {} B refused: {} at sequence {}",
            self.operation(),
            self.bytes,
            self.reason(),
            self.sequence
        )
    }
}

impl std::error::Error for PublishedTransferRefusal {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CloseRefusalReason {
    LiveAllocation,
    LiveChild,
    AlreadyReturned,
    ConservationMismatch,
    CounterOverflow,
}

/// Fail-closed child-return diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DelegatedLeaseCloseRefusal {
    identity: LeaseIdentity,
    logical_path: &'static str,
    reason: CloseRefusalReason,
    used_bytes: u64,
    live_children: u64,
    sequence: u128,
}

impl DelegatedLeaseCloseRefusal {
    /// Exact typed delegated identity.
    #[must_use]
    pub fn identity(&self) -> LeaseIdentity {
        self.identity
    }

    /// Stable refusal code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            CloseRefusalReason::LiveAllocation => "live_allocation",
            CloseRefusalReason::LiveChild => "live_child",
            CloseRefusalReason::AlreadyReturned => "already_returned",
            CloseRefusalReason::ConservationMismatch => "conservation_mismatch",
            CloseRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Child path.
    #[must_use]
    pub fn logical_path(&self) -> &'static str {
        self.logical_path
    }

    /// Bytes still occupied.
    #[must_use]
    pub fn used_bytes(&self) -> u64 {
        self.used_bytes
    }

    /// Direct children not yet returned.
    #[must_use]
    pub fn live_children(&self) -> u64 {
        self.live_children
    }

    /// Serialized root-ledger sequence.
    #[must_use]
    pub fn sequence(&self) -> u128 {
        self.sequence
    }
}

impl fmt::Display for DelegatedLeaseCloseRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "delegated memory lease `{}` could not return: {} at sequence {}",
            self.logical_path,
            self.reason(),
            self.sequence
        )
    }
}

impl std::error::Error for DelegatedLeaseCloseRefusal {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SealRefusalReason {
    UnverifiedRoot,
    LiveCapacity,
    ConservationMismatch,
    ReleaseInvariant,
    CounterOverflow,
}

/// Fail-closed reason a root could not yet produce a verified close receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseSealRefusal {
    root_identity: Option<LeaseIdentity>,
    root_id: Option<&'static str>,
    reason: SealRefusalReason,
    used_bytes: u64,
    active_delegations: u64,
    release_invariant_violations: u128,
    seal_sequence: u128,
    observation_sequence: u128,
}

impl LeaseSealRefusal {
    /// Typed configured root identity, when configuration completed.
    #[must_use]
    pub fn root_identity(&self) -> Option<LeaseIdentity> {
        self.root_identity
    }

    /// Static root diagnostic label, when configuration completed.
    #[must_use]
    pub fn root_id(&self) -> Option<&'static str> {
        self.root_id
    }

    /// Stable refusal code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        match self.reason {
            SealRefusalReason::UnverifiedRoot => "unverified_root",
            SealRefusalReason::LiveCapacity => "live_capacity",
            SealRefusalReason::ConservationMismatch => "conservation_mismatch",
            SealRefusalReason::ReleaseInvariant => "release_invariant",
            SealRefusalReason::CounterOverflow => "counter_overflow",
        }
    }

    /// Root bytes still live.
    #[must_use]
    pub fn used_bytes(&self) -> u64 {
        self.used_bytes
    }

    /// Child identities not yet returned.
    #[must_use]
    pub fn active_delegations(&self) -> u64 {
        self.active_delegations
    }

    /// Internal exact-return invariant failures.
    #[must_use]
    pub fn release_invariant_violations(&self) -> u128 {
        self.release_invariant_violations
    }

    /// Sequence that permanently froze admission.
    #[must_use]
    pub fn seal_sequence(&self) -> u128 {
        self.seal_sequence
    }

    /// Sequence of this close observation.
    #[must_use]
    pub fn observation_sequence(&self) -> u128 {
        self.observation_sequence
    }
}

impl fmt::Display for LeaseSealRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "operation memory lease could not close: {} ({} B live, {} child transfer(s), sequence {})",
            self.reason(),
            self.used_bytes,
            self.active_delegations,
            self.observation_sequence
        )
    }
}

impl std::error::Error for LeaseSealRefusal {}

/// Diagnostic returned when a close receipt does not verify.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeaseReceiptVerificationError {
    reason: &'static str,
}

impl LeaseReceiptVerificationError {
    /// Stable verification failure code.
    #[must_use]
    pub fn reason(&self) -> &'static str {
        self.reason
    }
}

impl fmt::Display for LeaseReceiptVerificationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "lease close receipt verification failed: {}",
            self.reason
        )
    }
}

impl std::error::Error for LeaseReceiptVerificationError {}

/// Verified evidence that one staging allocation transferred exactly once to
/// its bound destination.
///
/// This receipt is copyable evidence; the associated [`PublishedTransfer`]
/// owner is deliberately affine. Public fields permit serialization and
/// adversarial verification of decoded values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublishedTransferReceipt {
    /// Schema version.
    pub schema_version: u16,
    /// Typed root authority.
    pub root_identity: LeaseIdentity,
    /// Typed parent of the staging child.
    pub parent_identity: Option<LeaseIdentity>,
    /// Exact child that owned the staging allocation.
    pub child_identity: LeaseIdentity,
    /// Exact plan/occurrence/output/destination identity tuple.
    pub binding: PublishedTransferBinding,
    /// Exact payload/layout/overhead byte composition.
    pub envelope: PublishedTransferEnvelope,
    /// Transferred bytes.
    pub bytes: u64,
    /// Sequence that admitted the prepared transition.
    pub prepared_sequence: u128,
    /// Sequence that transferred ownership to the destination.
    pub published_sequence: u128,
    /// Deterministic root over every field above.
    pub receipt_root: [u8; 32],
}

impl PublishedTransferReceipt {
    /// Recompute the deterministic receipt root.
    #[must_use]
    pub fn recompute_root(&self) -> [u8; 32] {
        let mut hash = ReceiptHasher::new(b"fs-alloc-published-transfer-v1");
        hash.u16(self.schema_version);
        hash.identity(self.root_identity);
        hash.optional_identity(self.parent_identity);
        hash.identity(self.child_identity);
        hash.published_binding(self.binding);
        hash.published_envelope(self.envelope);
        hash.u64(self.bytes);
        hash.u128(self.prepared_sequence);
        hash.u128(self.published_sequence);
        hash.finish()
    }

    /// Verify the external authority context, sequencing, and receipt root.
    ///
    /// # Errors
    ///
    /// Returns a stable code for the first failed invariant.
    pub fn verify_for(
        &self,
        root_identity: LeaseIdentity,
        parent_identity: Option<LeaseIdentity>,
        child_identity: LeaseIdentity,
        binding: PublishedTransferBinding,
        envelope: PublishedTransferEnvelope,
    ) -> Result<(), LeaseReceiptVerificationError> {
        if self.schema_version != PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION {
            return verification_error("schema_version");
        }
        if self.root_identity != root_identity
            || self.parent_identity != parent_identity
            || self.child_identity != child_identity
            || self.binding != binding
            || self.envelope != envelope
        {
            return verification_error("identity");
        }
        if !self.root_identity.is_root()
            || self.parent_identity.is_none()
                && !self.child_identity.is_direct_child_of(self.root_identity)
            || self
                .parent_identity
                .is_some_and(|parent| !self.child_identity.is_direct_child_of(parent))
        {
            return verification_error("identity_relationship");
        }
        if self.bytes == 0 {
            return verification_error("zero_bytes");
        }
        if self.envelope.total_bytes() != Some(self.bytes) {
            return verification_error("envelope");
        }
        if self.published_sequence <= self.prepared_sequence {
            return verification_error("sequence");
        }
        if self.recompute_root() != self.receipt_root {
            return verification_error("receipt_root");
        }
        Ok(())
    }

    /// Canonical one-line JSON with deterministic field order.
    #[must_use]
    pub fn to_json(&self) -> String {
        use fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-published-transfer-v1\"");
        let _ = write!(
            out,
            ",\"schema_version\":{},\"root_identity\":{},\"parent_identity\":",
            self.schema_version,
            self.root_identity.to_json()
        );
        match self.parent_identity {
            Some(identity) => out.push_str(&identity.to_json()),
            None => out.push_str("null"),
        }
        let _ = write!(
            out,
            ",\"child_identity\":{},\"binding\":{},\"envelope\":{},\"bytes\":{},\"prepared_sequence\":{},\"published_sequence\":{},\"receipt_root\":\"{}\"}}",
            self.child_identity.to_json(),
            self.binding.to_json(),
            self.envelope.to_json(),
            self.bytes,
            self.prepared_sequence,
            self.published_sequence,
            hash_hex(self.receipt_root)
        );
        out
    }
}

/// Verified evidence that prepared staging ownership rolled back to its child.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublishedTransferRollbackReceipt {
    /// Schema version.
    pub schema_version: u16,
    /// Typed root authority.
    pub root_identity: LeaseIdentity,
    /// Typed parent of the staging child.
    pub parent_identity: Option<LeaseIdentity>,
    /// Exact child that regained the staging bytes.
    pub child_identity: LeaseIdentity,
    /// Exact attempted publication identity tuple.
    pub binding: PublishedTransferBinding,
    /// Exact payload/layout/overhead byte composition.
    pub envelope: PublishedTransferEnvelope,
    /// Returned staging bytes.
    pub bytes: u64,
    /// Sequence that admitted the prepared transition.
    pub prepared_sequence: u128,
    /// Sequence that returned staging ownership.
    pub rolled_back_sequence: u128,
    /// Whether `Drop` performed rollback.
    pub implicit_rollback: bool,
    /// Deterministic root over every field above.
    pub receipt_root: [u8; 32],
}

impl PublishedTransferRollbackReceipt {
    /// Recompute the deterministic receipt root.
    #[must_use]
    pub fn recompute_root(&self) -> [u8; 32] {
        let mut hash = ReceiptHasher::new(b"fs-alloc-published-rollback-v1");
        hash.u16(self.schema_version);
        hash.identity(self.root_identity);
        hash.optional_identity(self.parent_identity);
        hash.identity(self.child_identity);
        hash.published_binding(self.binding);
        hash.published_envelope(self.envelope);
        hash.u64(self.bytes);
        hash.u128(self.prepared_sequence);
        hash.u128(self.rolled_back_sequence);
        hash.boolean(self.implicit_rollback);
        hash.finish()
    }

    /// Verify exact authority context, sequencing, and root.
    ///
    /// # Errors
    ///
    /// Returns a stable code for the first failed invariant.
    pub fn verify_for(
        &self,
        root_identity: LeaseIdentity,
        parent_identity: Option<LeaseIdentity>,
        child_identity: LeaseIdentity,
        binding: PublishedTransferBinding,
        envelope: PublishedTransferEnvelope,
    ) -> Result<(), LeaseReceiptVerificationError> {
        if self.schema_version != PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION {
            return verification_error("schema_version");
        }
        if self.root_identity != root_identity
            || self.parent_identity != parent_identity
            || self.child_identity != child_identity
            || self.binding != binding
            || self.envelope != envelope
        {
            return verification_error("identity");
        }
        if !self.root_identity.is_root()
            || self.parent_identity.is_none()
                && !self.child_identity.is_direct_child_of(self.root_identity)
            || self
                .parent_identity
                .is_some_and(|parent| !self.child_identity.is_direct_child_of(parent))
        {
            return verification_error("identity_relationship");
        }
        if self.bytes == 0 {
            return verification_error("zero_bytes");
        }
        if self.envelope.total_bytes() != Some(self.bytes) {
            return verification_error("envelope");
        }
        if self.rolled_back_sequence <= self.prepared_sequence {
            return verification_error("sequence");
        }
        if self.recompute_root() != self.receipt_root {
            return verification_error("receipt_root");
        }
        Ok(())
    }

    /// Canonical one-line JSON with deterministic field order.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"fs-alloc-published-rollback-v1\",\"schema_version\":{},\"root_identity\":{},\"parent_identity\":{},\"child_identity\":{},\"binding\":{},\"envelope\":{},\"bytes\":{},\"prepared_sequence\":{},\"rolled_back_sequence\":{},\"implicit_rollback\":{},\"receipt_root\":\"{}\"}}",
            self.schema_version,
            self.root_identity.to_json(),
            self.parent_identity
                .map_or_else(|| String::from("null"), LeaseIdentity::to_json),
            self.child_identity.to_json(),
            self.binding.to_json(),
            self.envelope.to_json(),
            self.bytes,
            self.prepared_sequence,
            self.rolled_back_sequence,
            self.implicit_rollback,
            hash_hex(self.receipt_root)
        )
    }
}

/// Verified evidence that destination ownership ended after publication.
///
/// A root may seal before this receipt exists. Destination closure is a
/// separate lifecycle boundary and never rewrites the already-issued root
/// terminal receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PublishedTransferCloseReceipt {
    /// Schema version.
    pub schema_version: u16,
    /// Typed root authority.
    pub root_identity: LeaseIdentity,
    /// Typed parent of the staging child.
    pub parent_identity: Option<LeaseIdentity>,
    /// Exact child that staged the output.
    pub child_identity: LeaseIdentity,
    /// Exact publication identity tuple.
    pub binding: PublishedTransferBinding,
    /// Exact payload/layout/overhead byte composition.
    pub envelope: PublishedTransferEnvelope,
    /// Destination-owned bytes closed.
    pub bytes: u64,
    /// Root of the successful publication receipt.
    pub published_receipt_root: [u8; 32],
    /// Successful publication sequence.
    pub published_sequence: u128,
    /// Destination close sequence.
    pub closed_sequence: u128,
    /// Whether `Drop` performed destination closure.
    pub implicit_close: bool,
    /// Deterministic root over every field above.
    pub receipt_root: [u8; 32],
}

impl PublishedTransferCloseReceipt {
    /// Recompute the deterministic receipt root.
    #[must_use]
    pub fn recompute_root(&self) -> [u8; 32] {
        let mut hash = ReceiptHasher::new(b"fs-alloc-published-close-v1");
        hash.u16(self.schema_version);
        hash.identity(self.root_identity);
        hash.optional_identity(self.parent_identity);
        hash.identity(self.child_identity);
        hash.published_binding(self.binding);
        hash.published_envelope(self.envelope);
        hash.u64(self.bytes);
        hash.bytes(&self.published_receipt_root);
        hash.u128(self.published_sequence);
        hash.u128(self.closed_sequence);
        hash.boolean(self.implicit_close);
        hash.finish()
    }

    /// Verify this close against the successful publication receipt.
    ///
    /// # Errors
    ///
    /// Returns a stable code for the first failed invariant.
    pub fn verify_for(
        &self,
        published: &PublishedTransferReceipt,
    ) -> Result<(), LeaseReceiptVerificationError> {
        published.verify_for(
            published.root_identity,
            published.parent_identity,
            published.child_identity,
            published.binding,
            published.envelope,
        )?;
        if self.schema_version != PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION {
            return verification_error("schema_version");
        }
        if self.root_identity != published.root_identity
            || self.parent_identity != published.parent_identity
            || self.child_identity != published.child_identity
            || self.binding != published.binding
            || self.envelope != published.envelope
            || self.bytes != published.bytes
            || self.published_receipt_root != published.receipt_root
        {
            return verification_error("identity");
        }
        if self.published_sequence != published.published_sequence
            || self.closed_sequence <= self.published_sequence
        {
            return verification_error("sequence");
        }
        if self.recompute_root() != self.receipt_root {
            return verification_error("receipt_root");
        }
        Ok(())
    }

    /// Canonical one-line JSON with deterministic field order.
    #[must_use]
    pub fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"fs-alloc-published-close-v1\",\"schema_version\":{},\"root_identity\":{},\"parent_identity\":{},\"child_identity\":{},\"binding\":{},\"envelope\":{},\"bytes\":{},\"published_receipt_root\":\"{}\",\"published_sequence\":{},\"closed_sequence\":{},\"implicit_close\":{},\"receipt_root\":\"{}\"}}",
            self.schema_version,
            self.root_identity.to_json(),
            self.parent_identity
                .map_or_else(|| String::from("null"), LeaseIdentity::to_json),
            self.child_identity.to_json(),
            self.binding.to_json(),
            self.envelope.to_json(),
            self.bytes,
            hash_hex(self.published_receipt_root),
            self.published_sequence,
            self.closed_sequence,
            self.implicit_close,
            hash_hex(self.receipt_root)
        )
    }
}

/// Versioned terminal receipt for one returned delegated envelope.
///
/// Fields are public because serialized receipts are untrusted data. Callers
/// must use [`Self::verify_for`] before treating a decoded or modified value as
/// evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DelegatedLeaseCloseReceipt {
    /// Schema version.
    pub schema_version: u16,
    /// Typed root authority identity.
    pub root_identity: LeaseIdentity,
    /// Typed parent identity, or `None` for a direct root child.
    pub parent_identity: Option<LeaseIdentity>,
    /// Exact typed identity of this returned child.
    pub identity: LeaseIdentity,
    /// Exact root identity.
    pub root_id: &'static str,
    /// Logical parent path, or `None` for a direct root child.
    pub parent_path: Option<&'static str>,
    /// Exact full logical child path.
    pub logical_path: &'static str,
    /// Delegated envelope.
    pub capacity_bytes: u64,
    /// Exact cumulative child allocation grants.
    pub allocation_granted_bytes: u128,
    /// Exact cumulative child allocation returns.
    pub allocation_returned_bytes: u128,
    /// Exact cumulative child allocation bytes transferred to destinations.
    pub allocation_published_bytes: u128,
    /// Exact cumulative nested transfers.
    pub delegated_bytes: u128,
    /// Exact cumulative nested transfer returns.
    pub returned_delegated_bytes: u128,
    /// Child occupied high-water.
    pub peak_used_bytes: u64,
    /// Final occupied bytes; verified receipts require zero.
    pub final_used_bytes: u64,
    /// Exact refused request count.
    pub refused_requests: u128,
    /// Exact refused requested bytes.
    pub refused_bytes: u128,
    /// Number of prepared publication records retained for this child.
    pub publication_record_count: usize,
    /// Number of successful destination transfers for this child.
    pub published_transfer_count: usize,
    /// Number of prepared transfers rolled back to staging.
    pub rolled_back_transfer_count: usize,
    /// Deterministic root over publication records owned by this child.
    pub publication_root: [u8; 32],
    /// Sequence of the transfer.
    pub created_sequence: u128,
    /// Sequence of exact return.
    pub returned_sequence: u128,
    /// Whether `Drop` performed the return.
    pub implicit_return: bool,
    /// Deterministic receipt root over every field above.
    pub receipt_root: [u8; 32],
}

impl DelegatedLeaseCloseReceipt {
    /// Recompute the deterministic root from the receipt fields.
    #[must_use]
    pub fn recompute_root(&self) -> [u8; 32] {
        let mut hash = ReceiptHasher::new(b"fs-alloc-child-close-v2");
        hash.u16(self.schema_version);
        hash.identity(self.root_identity);
        hash.optional_identity(self.parent_identity);
        hash.identity(self.identity);
        hash.text(self.root_id);
        hash.optional_text(self.parent_path);
        hash.text(self.logical_path);
        hash.u64(self.capacity_bytes);
        hash.u128(self.allocation_granted_bytes);
        hash.u128(self.allocation_returned_bytes);
        hash.u128(self.allocation_published_bytes);
        hash.u128(self.delegated_bytes);
        hash.u128(self.returned_delegated_bytes);
        hash.u64(self.peak_used_bytes);
        hash.u64(self.final_used_bytes);
        hash.u128(self.refused_requests);
        hash.u128(self.refused_bytes);
        hash.usize(self.publication_record_count);
        hash.usize(self.published_transfer_count);
        hash.usize(self.rolled_back_transfer_count);
        hash.bytes(&self.publication_root);
        hash.u128(self.created_sequence);
        hash.u128(self.returned_sequence);
        hash.boolean(self.implicit_return);
        hash.finish()
    }

    /// Verify schema, identities, exact conservation, and the receipt root.
    ///
    /// # Errors
    ///
    /// Returns a stable mismatch code for the first failed invariant.
    pub fn verify_for(
        &self,
        root_identity: LeaseIdentity,
        parent_identity: Option<LeaseIdentity>,
        identity: LeaseIdentity,
    ) -> Result<(), LeaseReceiptVerificationError> {
        if self.schema_version != VERIFIED_RECEIPT_SCHEMA_VERSION {
            return verification_error("schema_version");
        }
        if self.root_identity != root_identity
            || self.parent_identity != parent_identity
            || self.identity != identity
        {
            return verification_error("identity");
        }
        if !self.root_identity.is_root()
            || self.identity.root_subject() != self.root_identity.root_subject()
            || self.identity.domain() != self.root_identity.domain()
            || self.parent_identity.is_none()
                && !self.identity.is_direct_child_of(self.root_identity)
            || self
                .parent_identity
                .is_some_and(|parent| !self.identity.is_direct_child_of(parent))
        {
            return verification_error("identity_relationship");
        }
        if self.final_used_bytes != 0 {
            return verification_error("live_capacity");
        }
        if self
            .allocation_returned_bytes
            .checked_add(self.allocation_published_bytes)
            != Some(self.allocation_granted_bytes)
            || self.delegated_bytes != self.returned_delegated_bytes
        {
            return verification_error("conservation");
        }
        if self
            .published_transfer_count
            .checked_add(self.rolled_back_transfer_count)
            != Some(self.publication_record_count)
        {
            return verification_error("publication_count");
        }
        if self.allocation_published_bytes == 0 && self.published_transfer_count != 0
            || self.allocation_published_bytes != 0 && self.published_transfer_count == 0
        {
            return verification_error("publication_count");
        }
        if self.peak_used_bytes > self.capacity_bytes {
            return verification_error("peak");
        }
        if self.returned_sequence <= self.created_sequence {
            return verification_error("sequence");
        }
        if self.recompute_root() != self.receipt_root {
            return verification_error("receipt_root");
        }
        Ok(())
    }

    /// Canonical one-line JSON with deterministic field order.
    #[must_use]
    pub fn to_json(&self) -> String {
        use fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-delegated-close-v2\"");
        let _ = write!(
            out,
            ",\"schema_version\":{},\"root_identity\":{},\"parent_identity\":",
            self.schema_version,
            self.root_identity.to_json()
        );
        match self.parent_identity {
            Some(identity) => out.push_str(&identity.to_json()),
            None => out.push_str("null"),
        }
        let _ = write!(
            out,
            ",\"identity\":{},\"root_id\":\"{}\",\"parent_path\":",
            self.identity.to_json(),
            json_escape(self.root_id)
        );
        match self.parent_path {
            Some(path) => {
                let _ = write!(out, "\"{}\"", json_escape(path));
            }
            None => out.push_str("null"),
        }
        let _ = write!(
            out,
            ",\"logical_path\":\"{}\",\"capacity_bytes\":{},\"allocation_granted_bytes\":{},\"allocation_returned_bytes\":{},\"allocation_published_bytes\":{},\"delegated_bytes\":{},\"returned_delegated_bytes\":{},\"peak_used_bytes\":{},\"final_used_bytes\":{},\"refused_requests\":{},\"refused_bytes\":{},\"publication_record_count\":{},\"published_transfer_count\":{},\"rolled_back_transfer_count\":{},\"publication_root\":\"{}\",\"created_sequence\":{},\"returned_sequence\":{},\"implicit_return\":{},\"receipt_root\":\"{}\"}}",
            json_escape(self.logical_path),
            self.capacity_bytes,
            self.allocation_granted_bytes,
            self.allocation_returned_bytes,
            self.allocation_published_bytes,
            self.delegated_bytes,
            self.returned_delegated_bytes,
            self.peak_used_bytes,
            self.final_used_bytes,
            self.refused_requests,
            self.refused_bytes,
            self.publication_record_count,
            self.published_transfer_count,
            self.rolled_back_transfer_count,
            hash_hex(self.publication_root),
            self.created_sequence,
            self.returned_sequence,
            self.implicit_return,
            hash_hex(self.receipt_root)
        );
        out
    }
}

/// Immutable successful terminal snapshot of a permanently sealed root.
///
/// The receipt covers the exact root identity, the frozen admission cut, all
/// direct and delegated conservation totals, the canonical child ledger root,
/// refusal totals, and the final deterministic receipt root.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SealedLeaseReceipt {
    /// Schema version.
    pub schema_version: u16,
    /// Typed configured root authority identity.
    pub root_identity: LeaseIdentity,
    /// Exact configured root identity.
    pub root_id: &'static str,
    /// Root byte limit.
    pub limit_bytes: u64,
    /// Pre-admitted maximum for each child and publication record ledger.
    pub metadata_limit: usize,
    /// Number of unique child identities ever transferred.
    pub delegation_count: usize,
    /// Exact direct allocation grants.
    pub direct_granted_bytes: u128,
    /// Exact direct allocation returns.
    pub direct_returned_bytes: u128,
    /// Exact top-level transferred capacity.
    pub delegated_bytes: u128,
    /// Exact top-level returned capacity.
    pub returned_delegated_bytes: u128,
    /// Exact allocation grants inside all child envelopes.
    pub child_granted_bytes: u128,
    /// Exact allocation returns inside all child envelopes.
    pub child_returned_bytes: u128,
    /// Exact allocation bytes transferred from child staging to destinations.
    pub child_published_bytes: u128,
    /// Number of prepared publication records retained by the root.
    pub publication_record_count: usize,
    /// Number of successful destination transfers.
    pub published_transfer_count: usize,
    /// Number of prepared transfers rolled back to staging.
    pub rolled_back_transfer_count: usize,
    /// Exact semantic refusal count before the terminal receipt.
    pub refused_requests: u128,
    /// Exact semantic refused bytes before the terminal receipt.
    pub refused_bytes: u128,
    /// Root occupied high-water.
    pub peak_used_bytes: u64,
    /// Final root live bytes; verified receipts require zero.
    pub final_used_bytes: u64,
    /// Final live child identities; verified receipts require zero.
    pub active_delegations: u64,
    /// Exact internal return failures; verified receipts require zero.
    pub release_invariant_violations: u128,
    /// Whether any exact counter overflowed.
    pub counter_overflowed: bool,
    /// Sequence that froze all new root and child admission.
    pub seal_sequence: u128,
    /// Sequence that produced the terminal receipt.
    pub close_sequence: u128,
    /// Deterministic accumulator covering all pre-close semantic refusals.
    pub refusal_root: [u8; 32],
    /// Deterministic root over child records sorted by logical path.
    pub delegation_root: [u8; 32],
    /// Deterministic root over publication records sorted by binding.
    pub publication_root: [u8; 32],
    /// Deterministic root over this entire receipt.
    pub receipt_root: [u8; 32],
    receipt: LeaseReceipt,
}

impl SealedLeaseReceipt {
    /// Compatibility accounting view captured at close.
    #[must_use]
    pub fn receipt(&self) -> &LeaseReceipt {
        &self.receipt
    }

    /// Recompute the deterministic receipt root.
    #[must_use]
    pub fn recompute_root(&self) -> [u8; 32] {
        let mut hash = ReceiptHasher::new(b"fs-alloc-root-close-v2");
        hash.u16(self.schema_version);
        hash.identity(self.root_identity);
        hash.text(self.root_id);
        hash.u64(self.limit_bytes);
        hash.usize(self.metadata_limit);
        hash.usize(self.delegation_count);
        hash.u128(self.direct_granted_bytes);
        hash.u128(self.direct_returned_bytes);
        hash.u128(self.delegated_bytes);
        hash.u128(self.returned_delegated_bytes);
        hash.u128(self.child_granted_bytes);
        hash.u128(self.child_returned_bytes);
        hash.u128(self.child_published_bytes);
        hash.usize(self.publication_record_count);
        hash.usize(self.published_transfer_count);
        hash.usize(self.rolled_back_transfer_count);
        hash.u128(self.refused_requests);
        hash.u128(self.refused_bytes);
        hash.u64(self.peak_used_bytes);
        hash.u64(self.final_used_bytes);
        hash.u64(self.active_delegations);
        hash.u128(self.release_invariant_violations);
        hash.boolean(self.counter_overflowed);
        hash.u128(self.seal_sequence);
        hash.u128(self.close_sequence);
        hash.bytes(&self.refusal_root);
        hash.bytes(&self.delegation_root);
        hash.bytes(&self.publication_root);
        hash.finish()
    }

    /// Verify the expected root identity, terminal state, conservation, and
    /// deterministic root.
    ///
    /// # Errors
    ///
    /// Returns a stable mismatch code for the first failed invariant.
    pub fn verify_for(
        &self,
        root_identity: LeaseIdentity,
    ) -> Result<(), LeaseReceiptVerificationError> {
        if self.schema_version != VERIFIED_RECEIPT_SCHEMA_VERSION {
            return verification_error("schema_version");
        }
        if self.root_identity != root_identity || !self.root_identity.is_root() {
            return verification_error("identity");
        }
        if self.final_used_bytes != 0 || self.active_delegations != 0 {
            return verification_error("live_capacity");
        }
        if self.direct_granted_bytes != self.direct_returned_bytes
            || self.delegated_bytes != self.returned_delegated_bytes
            || self
                .child_returned_bytes
                .checked_add(self.child_published_bytes)
                != Some(self.child_granted_bytes)
        {
            return verification_error("conservation");
        }
        if self
            .published_transfer_count
            .checked_add(self.rolled_back_transfer_count)
            != Some(self.publication_record_count)
        {
            return verification_error("publication_count");
        }
        if self.child_published_bytes == 0 && self.published_transfer_count != 0
            || self.child_published_bytes != 0 && self.published_transfer_count == 0
        {
            return verification_error("publication_count");
        }
        if self.release_invariant_violations != 0 {
            return verification_error("release_invariant");
        }
        if self.counter_overflowed {
            return verification_error("counter_overflow");
        }
        if self.peak_used_bytes > self.limit_bytes {
            return verification_error("peak");
        }
        if self.delegation_count > self.metadata_limit
            || self.publication_record_count > self.metadata_limit
        {
            return verification_error("metadata_limit");
        }
        if self.close_sequence < self.seal_sequence {
            return verification_error("sequence");
        }
        if self.recompute_root() != self.receipt_root {
            return verification_error("receipt_root");
        }
        Ok(())
    }

    /// Canonical one-line JSON with deterministic field order.
    #[must_use]
    pub fn to_json(&self) -> String {
        use fmt::Write as _;
        let mut out = String::from("{\"schema\":\"fs-alloc-root-close-v2\"");
        let _ = write!(
            out,
            ",\"schema_version\":{},\"root_identity\":{},\"root_id\":\"{}\",\"limit_bytes\":{},\"metadata_limit\":{},\"delegation_count\":{},\"direct_granted_bytes\":{},\"direct_returned_bytes\":{},\"delegated_bytes\":{},\"returned_delegated_bytes\":{},\"child_granted_bytes\":{},\"child_returned_bytes\":{},\"child_published_bytes\":{},\"publication_record_count\":{},\"published_transfer_count\":{},\"rolled_back_transfer_count\":{},\"refused_requests\":{},\"refused_bytes\":{},\"peak_used_bytes\":{},\"final_used_bytes\":{},\"active_delegations\":{},\"release_invariant_violations\":{},\"counter_overflowed\":{},\"seal_sequence\":{},\"close_sequence\":{},\"refusal_root\":\"{}\",\"delegation_root\":\"{}\",\"publication_root\":\"{}\",\"receipt_root\":\"{}\"}}",
            self.schema_version,
            self.root_identity.to_json(),
            json_escape(self.root_id),
            self.limit_bytes,
            self.metadata_limit,
            self.delegation_count,
            self.direct_granted_bytes,
            self.direct_returned_bytes,
            self.delegated_bytes,
            self.returned_delegated_bytes,
            self.child_granted_bytes,
            self.child_returned_bytes,
            self.child_published_bytes,
            self.publication_record_count,
            self.published_transfer_count,
            self.rolled_back_transfer_count,
            self.refused_requests,
            self.refused_bytes,
            self.peak_used_bytes,
            self.final_used_bytes,
            self.active_delegations,
            self.release_invariant_violations,
            self.counter_overflowed,
            self.seal_sequence,
            self.close_sequence,
            hash_hex(self.refusal_root),
            hash_hex(self.delegation_root),
            hash_hex(self.publication_root),
            hash_hex(self.receipt_root)
        );
        out
    }
}

fn verification_error<T>(reason: &'static str) -> Result<T, LeaseReceiptVerificationError> {
    Err(LeaseReceiptVerificationError { reason })
}

pub(crate) fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control <= '\u{1f}' => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}

#[derive(Clone, Copy)]
struct ReceiptHasher {
    lanes: [u64; 4],
}

impl ReceiptHasher {
    fn new(domain: &[u8]) -> Self {
        let mut hash = Self {
            lanes: [
                0xcbf2_9ce4_8422_2325,
                0x9e37_79b9_7f4a_7c15,
                0x243f_6a88_85a3_08d3,
                0x1319_8a2e_0370_7344,
            ],
        };
        hash.bytes(domain);
        hash
    }

    fn bytes(&mut self, bytes: &[u8]) {
        self.raw_bytes(&(bytes.len() as u128).to_le_bytes());
        self.raw_bytes(bytes);
    }

    fn raw_bytes(&mut self, bytes: &[u8]) {
        for (index, byte) in bytes.iter().copied().enumerate() {
            let lane = index & 3;
            self.lanes[lane] ^= u64::from(byte);
            self.lanes[lane] = self.lanes[lane]
                .wrapping_mul(0x0000_0100_0000_01b3)
                .rotate_left(13 + u32::try_from(lane).expect("lane fits") * 7);
            let next = (lane + 1) & 3;
            self.lanes[next] ^= self.lanes[lane].rotate_right(17);
        }
    }

    fn text(&mut self, value: &str) {
        self.bytes(value.as_bytes());
    }

    fn optional_text(&mut self, value: Option<&str>) {
        match value {
            Some(value) => {
                self.boolean(true);
                self.text(value);
            }
            None => self.boolean(false),
        }
    }

    fn identity(&mut self, value: LeaseIdentity) {
        self.bytes(&value.canonical_bytes());
    }

    fn optional_identity(&mut self, value: Option<LeaseIdentity>) {
        match value {
            Some(value) => {
                self.boolean(true);
                self.identity(value);
            }
            None => self.boolean(false),
        }
    }

    fn published_binding(&mut self, value: PublishedTransferBinding) {
        self.bytes(&value.canonical_bytes());
    }

    fn published_envelope(&mut self, value: PublishedTransferEnvelope) {
        self.bytes(&value.canonical_bytes());
    }

    fn boolean(&mut self, value: bool) {
        self.bytes(&[u8::from(value)]);
    }

    fn u16(&mut self, value: u16) {
        self.bytes(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes(&value.to_le_bytes());
    }

    fn u128(&mut self, value: u128) {
        self.bytes(&value.to_le_bytes());
    }

    fn usize(&mut self, value: usize) {
        self.u128(value as u128);
    }

    fn finish(mut self) -> [u8; 32] {
        for round in 0..12 {
            let lane = round & 3;
            let next = (lane + 1) & 3;
            self.lanes[lane] ^= self.lanes[next].rotate_left(19);
            self.lanes[lane] = self.lanes[lane]
                .wrapping_mul(0x9e37_79b1_85eb_ca87)
                .rotate_left(23);
        }
        let mut out = [0_u8; 32];
        for (index, lane) in self.lanes.into_iter().enumerate() {
            out[index * 8..(index + 1) * 8].copy_from_slice(&lane.to_le_bytes());
        }
        out
    }
}

fn hash_hex(hash: [u8; 32]) -> String {
    byte_slice_hex(&hash)
}

fn byte_slice_hex(bytes: &[u8]) -> String {
    use fmt::Write as _;
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        let _ = write!(out, "{byte:02x}");
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DelegationDisposition {
    Live,
    Returned,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PublicationDisposition {
    Prepared,
    Published,
    RolledBack,
}

#[derive(Debug, Clone)]
struct PublicationRecord {
    parent_identity: Option<LeaseIdentity>,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    prepared_sequence: u128,
    resolved_sequence: Option<u128>,
    disposition: PublicationDisposition,
    implicit_rollback: bool,
    published_receipt_root: Option<[u8; 32]>,
    destination_closed_sequence: Option<u128>,
    implicit_destination_close: bool,
    destination_close_root: Option<[u8; 32]>,
}

#[derive(Debug, Clone)]
struct DelegationRecord {
    parent_identity: Option<LeaseIdentity>,
    identity: LeaseIdentity,
    parent_path: Option<&'static str>,
    logical_path: &'static str,
    capacity_bytes: u64,
    used_bytes: u64,
    peak_used_bytes: u64,
    allocation_granted_bytes: u128,
    allocation_returned_bytes: u128,
    allocation_published_bytes: u128,
    delegated_bytes: u128,
    returned_delegated_bytes: u128,
    refused_requests: u128,
    refused_bytes: u128,
    live_children: u64,
    created_sequence: u128,
    returned_sequence: Option<u128>,
    implicit_return: bool,
    close_root: Option<[u8; 32]>,
    disposition: DelegationDisposition,
}

struct LeaseState {
    used_bytes: u64,
    peak_bytes: u64,
    requested_bytes: u128,
    refusals: u128,
    release_invariant_violations: u128,
    first_refusal: Option<LeaseRefusal>,
    sequence: u128,
    counter_overflowed: bool,
    sealed: bool,
    seal_sequence: Option<u128>,
    terminal_receipt: Option<SealedLeaseReceipt>,
    root_identity: Option<LeaseIdentity>,
    root_id: Option<&'static str>,
    metadata_limit: usize,
    delegations: Vec<DelegationRecord>,
    publications: Vec<PublicationRecord>,
    direct_granted_bytes: u128,
    direct_returned_bytes: u128,
    delegated_bytes: u128,
    returned_delegated_bytes: u128,
    child_granted_bytes: u128,
    child_returned_bytes: u128,
    child_published_bytes: u128,
    semantic_refusals: u128,
    semantic_refused_bytes: u128,
    refusal_root: [u8; 32],
    post_terminal_refusals: u128,
}

impl LeaseState {
    fn new() -> Self {
        Self {
            used_bytes: 0,
            peak_bytes: 0,
            requested_bytes: 0,
            refusals: 0,
            release_invariant_violations: 0,
            first_refusal: None,
            sequence: 0,
            counter_overflowed: false,
            sealed: false,
            seal_sequence: None,
            terminal_receipt: None,
            root_identity: None,
            root_id: None,
            metadata_limit: 0,
            delegations: Vec::new(),
            publications: Vec::new(),
            direct_granted_bytes: 0,
            direct_returned_bytes: 0,
            delegated_bytes: 0,
            returned_delegated_bytes: 0,
            child_granted_bytes: 0,
            child_returned_bytes: 0,
            child_published_bytes: 0,
            semantic_refusals: 0,
            semantic_refused_bytes: 0,
            refusal_root: ReceiptHasher::new(b"fs-alloc-refusals-v1").finish(),
            post_terminal_refusals: 0,
        }
    }
}

struct LeaseShared {
    limit_bytes: Option<u64>,
    state: Mutex<LeaseState>,
}

/// Cloneable run-scoped root memory lease. Every reserve, release, transfer,
/// and seal transition is serialized by one mutex-protected state machine.
#[derive(Clone)]
pub struct OperationMemoryLease {
    shared: Arc<LeaseShared>,
}

impl fmt::Debug for OperationMemoryLease {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OperationMemoryLease")
            .field("receipt", &self.receipt().to_json())
            .finish_non_exhaustive()
    }
}

impl OperationMemoryLease {
    /// A lease enforcing a hard byte limit.
    #[must_use]
    pub fn bounded(limit_bytes: u64) -> Self {
        Self::with_limit(Some(limit_bytes))
    }

    /// The compatibility lease: accounting only, with no bounded-memory close
    /// authority.
    #[must_use]
    pub fn unbounded() -> Self {
        Self::with_limit(None)
    }

    fn with_limit(limit_bytes: Option<u64>) -> Self {
        Self {
            shared: Arc::new(LeaseShared {
                limit_bytes,
                state: Mutex::new(LeaseState::new()),
            }),
        }
    }

    /// Enable verified affine delegation on a pristine bounded root.
    ///
    /// This configures an already allocated root, so the fallible metadata
    /// pre-admission cannot be followed by a hidden `Arc` allocation. All
    /// child and publication records are reserved up front. Zero-capacity
    /// children remain real affine control authorities: they consume metadata,
    /// block seal while live, and must return exactly once.
    ///
    /// # Errors
    ///
    /// Returns a typed refusal for an unbounded, used, sealed, invalid,
    /// already-configured, or unreservable root.
    #[allow(clippy::result_large_err)]
    pub fn enable_delegation(
        &self,
        root_identity: LeaseIdentity,
        root_id: &'static str,
        metadata_limit: usize,
    ) -> Result<(), LeaseConfigurationRefusal> {
        let mut state = self.lock_state();
        let sequence = next_sequence(&mut state);
        let reason = if self.shared.limit_bytes.is_none() {
            Some(ConfigurationRefusalReason::UnboundedRoot)
        } else if state.sealed {
            Some(ConfigurationRefusalReason::Sealed)
        } else if state.root_identity.is_some() {
            Some(ConfigurationRefusalReason::AlreadyConfigured)
        } else if state.used_bytes != 0
            || state.requested_bytes != 0
            || state.refusals != 0
            || state.release_invariant_violations != 0
        {
            Some(ConfigurationRefusalReason::NotPristine)
        } else if !root_identity.is_root() || !valid_root_id(root_id) {
            Some(ConfigurationRefusalReason::InvalidRootIdentity)
        } else if metadata_limit == 0 || metadata_limit > MAX_DELEGATION_RECORDS {
            Some(ConfigurationRefusalReason::MetadataLimit)
        } else {
            None
        };
        if let Some(reason) = reason {
            record_semantic_refusal(
                &mut state,
                reason_code_configuration(reason),
                Some(root_identity),
                root_id,
                None,
                0,
                0,
                sequence,
            );
            return Err(LeaseConfigurationRefusal {
                root_identity,
                root_id,
                reason,
                sequence,
            });
        }
        if state.delegations.try_reserve_exact(metadata_limit).is_err()
            || state
                .publications
                .try_reserve_exact(metadata_limit)
                .is_err()
        {
            let reason = ConfigurationRefusalReason::MetadataAllocation;
            record_semantic_refusal(
                &mut state,
                reason_code_configuration(reason),
                Some(root_identity),
                root_id,
                None,
                0,
                0,
                sequence,
            );
            return Err(LeaseConfigurationRefusal {
                root_identity,
                root_id,
                reason,
                sequence,
            });
        }
        state.root_identity = Some(root_identity);
        state.root_id = Some(root_id);
        state.metadata_limit = metadata_limit;
        Ok(())
    }

    /// The limit in force.
    #[must_use]
    pub fn limit_bytes(&self) -> Option<u64> {
        self.shared.limit_bytes
    }

    /// Number of child identities not yet returned.
    #[must_use]
    pub fn active_delegations(&self) -> u64 {
        let state = self.lock_state();
        u64::try_from(
            state
                .delegations
                .iter()
                .filter(|record| record.disposition == DelegationDisposition::Live)
                .count(),
        )
        .unwrap_or(u64::MAX)
    }

    /// Whether the root admission cut has been frozen.
    #[must_use]
    pub fn is_sealed(&self) -> bool {
        self.lock_state().sealed
    }

    /// Transfer capacity to one exact, caller-supplied logical path.
    ///
    /// A zero-capacity child still owns a unique affine control token, consumes
    /// one retained metadata record, and blocks root seal until returned.
    ///
    /// The returned owner is affine: it is neither `Clone` nor `Copy`, cannot
    /// be substituted for an ordinary root lease, and is lifetime-bound to
    /// this root handle.
    ///
    /// ```compile_fail
    /// use fs_alloc::{LeaseIdentity, OperationMemoryLease};
    ///
    /// let root = OperationMemoryLease::bounded(8);
    /// let root_id = LeaseIdentity::root(*b"example1", [1; 32]);
    /// let child_id = root_id.child([2; 32], 0).unwrap();
    /// root.enable_delegation(root_id, "run", 1).unwrap();
    /// let child = root.delegate_capacity(child_id, "run/child", 8).unwrap();
    /// let _duplicate_owner = child.clone();
    /// ```
    ///
    /// ```compile_fail
    /// use fs_alloc::{LeaseIdentity, LeasedVec, OperationMemoryLease};
    ///
    /// let root = OperationMemoryLease::bounded(8);
    /// let root_id = LeaseIdentity::root(*b"example2", [1; 32]);
    /// let child_id = root_id.child([2; 32], 0).unwrap();
    /// root.enable_delegation(root_id, "run", 1).unwrap();
    /// let child = root.delegate_capacity(child_id, "run/child", 8).unwrap();
    /// let _ = LeasedVec::<u8>::with_capacity(&child, "payload", 8);
    /// ```
    ///
    /// ```compile_fail
    /// use fs_alloc::{LeaseIdentity, OperationMemoryLease};
    ///
    /// let child = {
    ///     let root = OperationMemoryLease::bounded(8);
    ///     let root_id = LeaseIdentity::root(*b"example3", [1; 32]);
    ///     let child_id = root_id.child([2; 32], 0).unwrap();
    ///     root.enable_delegation(root_id, "run", 1).unwrap();
    ///     root.delegate_capacity(child_id, "run/child", 8).unwrap()
    /// };
    /// let _ = child.capacity_bytes();
    /// ```
    #[allow(clippy::result_large_err)]
    pub fn delegate_capacity<'root>(
        &'root self,
        identity: LeaseIdentity,
        logical_path: &'static str,
        capacity_bytes: u64,
    ) -> Result<DelegatedMemoryLease<'root>, LeaseDelegationRefusal> {
        self.delegate_from(None, None, identity, logical_path, capacity_bytes)?;
        Ok(DelegatedMemoryLease {
            root: self.clone(),
            parent_identity: None,
            identity,
            parent_path: None,
            logical_path,
            capacity_bytes,
            closed: false,
            _parent: PhantomData,
        })
    }

    // Keeping this transition together makes the single-lock admission order
    // and its refusal snapshot directly auditable. Refusals deliberately carry
    // fixed-width identities and complete allocation-free diagnostic facts.
    #[allow(clippy::result_large_err, clippy::too_many_lines)]
    fn delegate_from(
        &self,
        parent_identity: Option<LeaseIdentity>,
        parent_path: Option<&'static str>,
        identity: LeaseIdentity,
        logical_path: &'static str,
        capacity_bytes: u64,
    ) -> Result<(), LeaseDelegationRefusal> {
        let mut state = self.lock_state();
        let sequence = next_sequence(&mut state);
        let root_identity = state.root_identity;
        let root_id = state.root_id;
        let parent_snapshot = parent_identity.and_then(|identity| {
            state
                .delegations
                .iter()
                .find(|record| record.identity == identity)
                .map(|record| (record.used_bytes, record.capacity_bytes, record.disposition))
        });
        let (used_bytes, limit_bytes) = match parent_snapshot {
            Some((used, limit, _)) => (used, Some(limit)),
            None if parent_identity.is_some() => (0, None),
            None => (state.used_bytes, self.shared.limit_bytes),
        };

        let reason = if self.shared.limit_bytes.is_none() {
            Some(DelegationRefusalReason::UnboundedParent)
        } else if root_identity.is_none() || root_id.is_none() {
            Some(DelegationRefusalReason::UnconfiguredRoot)
        } else if state.sealed {
            Some(DelegationRefusalReason::RootSealed)
        } else if !identity.is_direct_child_of(
            parent_identity
                .unwrap_or_else(|| root_identity.expect("configured root identity checked")),
        ) {
            Some(DelegationRefusalReason::InvalidIdentityRelationship)
        } else if !valid_logical_path(root_id.expect("checked"), parent_path, logical_path) {
            Some(DelegationRefusalReason::InvalidLogicalPath)
        } else if state
            .delegations
            .iter()
            .any(|record| record.identity == identity)
        {
            Some(DelegationRefusalReason::DuplicateIdentity)
        } else if state.delegations.iter().any(|record| {
            record.parent_identity == parent_identity && record.identity.path() == identity.path()
        }) {
            Some(DelegationRefusalReason::DuplicatePath)
        } else if state.delegations.len() >= state.metadata_limit {
            Some(DelegationRefusalReason::MetadataExhausted)
        } else if parent_identity.is_some()
            && parent_snapshot
                .is_none_or(|(_, _, disposition)| disposition == DelegationDisposition::Returned)
        {
            Some(DelegationRefusalReason::ParentReturned)
        } else if used_bytes
            .checked_add(capacity_bytes)
            .is_none_or(|next| limit_bytes.is_none_or(|limit| next > limit))
        {
            Some(DelegationRefusalReason::Capacity)
        } else if state.counter_overflowed
            || state
                .requested_bytes
                .checked_add(u128::from(capacity_bytes))
                .is_none()
            || (parent_identity.is_none()
                && state
                    .delegated_bytes
                    .checked_add(u128::from(capacity_bytes))
                    .is_none())
            || parent_identity.is_some_and(|identity| {
                let record = state
                    .delegations
                    .iter()
                    .find(|record| record.identity == identity)
                    .expect("live parent checked");
                record
                    .delegated_bytes
                    .checked_add(u128::from(capacity_bytes))
                    .is_none()
                    || record.live_children.checked_add(1).is_none()
            })
        {
            Some(DelegationRefusalReason::CounterOverflow)
        } else {
            None
        };

        if let Some(reason) = reason {
            if reason == DelegationRefusalReason::CounterOverflow {
                state.counter_overflowed = true;
            }
            if let Some(parent_identity) = parent_identity
                && let Some(index) = state
                    .delegations
                    .iter()
                    .position(|record| record.identity == parent_identity)
            {
                let next_refused_requests =
                    state.delegations[index].refused_requests.checked_add(1);
                let next_refused_bytes = state.delegations[index]
                    .refused_bytes
                    .checked_add(u128::from(capacity_bytes));
                if let (Some(refused_requests), Some(refused_bytes)) =
                    (next_refused_requests, next_refused_bytes)
                {
                    let record = &mut state.delegations[index];
                    record.refused_requests = refused_requests;
                    record.refused_bytes = refused_bytes;
                } else {
                    state.counter_overflowed = true;
                }
            }
            record_semantic_refusal(
                &mut state,
                reason_code_delegation(reason),
                Some(identity),
                logical_path,
                Some(capacity_bytes),
                used_bytes,
                limit_bytes.unwrap_or(u64::MAX),
                sequence,
            );
            return Err(LeaseDelegationRefusal {
                root_identity,
                parent_identity,
                identity,
                root_id,
                parent_path,
                logical_path,
                requested_bytes: capacity_bytes,
                used_bytes,
                limit_bytes,
                reason,
                sequence,
            });
        }

        if let Some(parent_identity) = parent_identity {
            let parent = state
                .delegations
                .iter_mut()
                .find(|record| record.identity == parent_identity)
                .expect("validated parent");
            parent.used_bytes += capacity_bytes;
            parent.peak_used_bytes = parent.peak_used_bytes.max(parent.used_bytes);
            parent.delegated_bytes += u128::from(capacity_bytes);
            parent.live_children += 1;
        } else {
            state.used_bytes += capacity_bytes;
            state.peak_bytes = state.peak_bytes.max(state.used_bytes);
            state.requested_bytes += u128::from(capacity_bytes);
            state.delegated_bytes += u128::from(capacity_bytes);
        }
        state.delegations.push(DelegationRecord {
            parent_identity,
            identity,
            parent_path,
            logical_path,
            capacity_bytes,
            used_bytes: 0,
            peak_used_bytes: 0,
            allocation_granted_bytes: 0,
            allocation_returned_bytes: 0,
            allocation_published_bytes: 0,
            delegated_bytes: 0,
            returned_delegated_bytes: 0,
            refused_requests: 0,
            refused_bytes: 0,
            live_children: 0,
            created_sequence: sequence,
            returned_sequence: None,
            implicit_return: false,
            close_root: None,
            disposition: DelegationDisposition::Live,
        });
        Ok(())
    }

    /// Permanently freeze all root and child admission, then produce a
    /// verified terminal receipt once existing charges and transfers drain.
    ///
    /// Seal, reserve, and delegate are one mutex-serialized transition family:
    /// a race has exactly one winner, and the first seal sequence is immutable.
    #[allow(clippy::too_many_lines)] // one mutex-serialized transition family with an immutable first-seal sequence; splitting would scatter the race invariants xdu4's proof lane reasons about
    #[allow(clippy::result_large_err)]
    pub fn seal(&self) -> Result<SealedLeaseReceipt, LeaseSealRefusal> {
        let mut state = self.lock_state();
        if let Some(receipt) = state.terminal_receipt.as_ref() {
            return Ok(receipt.clone());
        }
        let observation_sequence = next_sequence(&mut state);
        let seal_sequence = *state.seal_sequence.get_or_insert(observation_sequence);
        state.sealed = true;
        let active_delegations = count_live_delegations(&state);
        let reason = if state.root_identity.is_none()
            || state.root_id.is_none()
            || self.shared.limit_bytes.is_none()
        {
            Some(SealRefusalReason::UnverifiedRoot)
        } else if state.counter_overflowed {
            Some(SealRefusalReason::CounterOverflow)
        } else if state.release_invariant_violations != 0 {
            Some(SealRefusalReason::ReleaseInvariant)
        } else if state.used_bytes != 0 || active_delegations != 0 {
            Some(SealRefusalReason::LiveCapacity)
        } else if state.direct_granted_bytes != state.direct_returned_bytes
            || state.delegated_bytes != state.returned_delegated_bytes
            || state
                .child_returned_bytes
                .checked_add(state.child_published_bytes)
                != Some(state.child_granted_bytes)
            || state.delegations.iter().any(|record| {
                record.disposition != DelegationDisposition::Returned
                    || record.used_bytes != 0
                    || record.live_children != 0
                    || record
                        .allocation_returned_bytes
                        .checked_add(record.allocation_published_bytes)
                        != Some(record.allocation_granted_bytes)
                    || record.delegated_bytes != record.returned_delegated_bytes
            })
            || state
                .publications
                .iter()
                .any(|record| record.disposition == PublicationDisposition::Prepared)
        {
            Some(SealRefusalReason::ConservationMismatch)
        } else {
            None
        };
        if let Some(reason) = reason {
            let refusal_site = state.root_id.unwrap_or("unverified-root");
            let refusal_identity = state.root_identity;
            let refusal_used = state.used_bytes;
            let refusal_active = u64::try_from(active_delegations).unwrap_or(u64::MAX);
            record_semantic_refusal(
                &mut state,
                reason_code_seal(reason),
                refusal_identity,
                refusal_site,
                None,
                refusal_used,
                refusal_active,
                observation_sequence,
            );
            return Err(LeaseSealRefusal {
                root_identity: refusal_identity,
                root_id: state.root_id,
                reason,
                used_bytes: refusal_used,
                active_delegations: refusal_active,
                release_invariant_violations: state.release_invariant_violations,
                seal_sequence,
                observation_sequence,
            });
        }

        state
            .delegations
            .sort_unstable_by_key(|record| record.identity);
        state.publications.sort_unstable_by_key(|record| {
            (
                record.binding,
                record.child_identity,
                record.prepared_sequence,
            )
        });
        let delegation_root = delegation_root(&state.delegations);
        let publication_root = publication_root(&state.publications);
        let published_transfer_count = state
            .publications
            .iter()
            .filter(|record| record.disposition == PublicationDisposition::Published)
            .count();
        let rolled_back_transfer_count = state
            .publications
            .iter()
            .filter(|record| record.disposition == PublicationDisposition::RolledBack)
            .count();
        let receipt = lease_receipt_locked(self.shared.limit_bytes, &state);
        let mut terminal = SealedLeaseReceipt {
            schema_version: VERIFIED_RECEIPT_SCHEMA_VERSION,
            root_identity: state.root_identity.expect("verified root identity"),
            root_id: state.root_id.expect("verified root"),
            limit_bytes: self.shared.limit_bytes.expect("verified bounded root"),
            metadata_limit: state.metadata_limit,
            delegation_count: state.delegations.len(),
            direct_granted_bytes: state.direct_granted_bytes,
            direct_returned_bytes: state.direct_returned_bytes,
            delegated_bytes: state.delegated_bytes,
            returned_delegated_bytes: state.returned_delegated_bytes,
            child_granted_bytes: state.child_granted_bytes,
            child_returned_bytes: state.child_returned_bytes,
            child_published_bytes: state.child_published_bytes,
            publication_record_count: state.publications.len(),
            published_transfer_count,
            rolled_back_transfer_count,
            refused_requests: state.semantic_refusals,
            refused_bytes: state.semantic_refused_bytes,
            peak_used_bytes: state.peak_bytes,
            final_used_bytes: state.used_bytes,
            active_delegations: 0,
            release_invariant_violations: state.release_invariant_violations,
            counter_overflowed: state.counter_overflowed,
            seal_sequence,
            close_sequence: observation_sequence,
            refusal_root: state.refusal_root,
            delegation_root,
            publication_root,
            receipt_root: [0; 32],
            receipt,
        };
        terminal.receipt_root = terminal.recompute_root();
        terminal
            .verify_for(terminal.root_identity)
            .expect("constructed terminal receipt verifies");
        state.terminal_receipt = Some(terminal.clone());
        Ok(terminal)
    }

    /// Reserve root bytes, returning an exact RAII charge.
    pub fn reserve(&self, what: &'static str, bytes: u64) -> Result<LeaseCharge, LeaseRefusal> {
        let mut state = self.lock_state();
        let sequence = next_sequence(&mut state);
        let used = state.used_bytes;
        let reason = if state.sealed {
            Some(LeaseRefusalReason::Sealed)
        } else if state.counter_overflowed
            || state
                .requested_bytes
                .checked_add(u128::from(bytes))
                .is_none()
            || (state.root_id.is_some()
                && state
                    .direct_granted_bytes
                    .checked_add(u128::from(bytes))
                    .is_none())
        {
            Some(LeaseRefusalReason::CounterOverflow)
        } else if used
            .checked_add(bytes)
            .is_none_or(|next| self.shared.limit_bytes.is_some_and(|limit| next > limit))
        {
            Some(LeaseRefusalReason::Capacity)
        } else {
            None
        };
        if let Some(reason) = reason {
            if reason == LeaseRefusalReason::CounterOverflow {
                state.counter_overflowed = true;
            }
            let refusal = record_root_refusal(
                &mut state,
                self.shared.limit_bytes,
                what,
                bytes,
                used,
                reason,
                sequence,
            );
            return Err(refusal);
        }
        state.used_bytes += bytes;
        state.peak_bytes = state.peak_bytes.max(state.used_bytes);
        state.requested_bytes += u128::from(bytes);
        if state.root_id.is_some() {
            state.direct_granted_bytes += u128::from(bytes);
        }
        drop(state);
        Ok(LeaseCharge {
            lease: self.clone(),
            bytes,
        })
    }

    /// Side-effect-free admission hint. The subsequent reserve remains the
    /// serialized authority.
    pub(crate) fn can_reserve_now(&self, bytes: u64) -> bool {
        let state = self.lock_state();
        !state.sealed
            && !state.counter_overflowed
            && state
                .used_bytes
                .checked_add(bytes)
                .is_some_and(|next| self.shared.limit_bytes.is_none_or(|limit| next <= limit))
    }

    /// Release bytes whose RAII guard was transferred to an aggregate owner.
    pub(crate) fn release_raw(&self, bytes: u64) -> bool {
        let mut state = self.lock_state();
        let _sequence = next_sequence(&mut state);
        let Some(next_used) = state.used_bytes.checked_sub(bytes) else {
            increment_invariant_violation(&mut state);
            return false;
        };
        let next_returned = if state.root_id.is_some() {
            let Some(next) = state.direct_returned_bytes.checked_add(u128::from(bytes)) else {
                state.counter_overflowed = true;
                increment_invariant_violation(&mut state);
                return false;
            };
            Some(next)
        } else {
            None
        };
        state.used_bytes = next_used;
        if let Some(next_returned) = next_returned {
            state.direct_returned_bytes = next_returned;
        }
        true
    }

    fn lock_state(&self) -> std::sync::MutexGuard<'_, LeaseState> {
        self.shared
            .state
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Accounting snapshot with canonical serialization.
    #[must_use]
    pub fn receipt(&self) -> LeaseReceipt {
        lease_receipt_locked(self.shared.limit_bytes, &self.lock_state())
    }
}

/// Non-cloneable owner of one exact delegated capacity envelope.
///
/// Allocation charges borrow this owner, and nested owners borrow their
/// parent. Safe Rust therefore prevents owner return before live charges or
/// children. Runtime checks remain fail-closed for deliberately leaked values.
pub struct DelegatedMemoryLease<'parent> {
    root: OperationMemoryLease,
    parent_identity: Option<LeaseIdentity>,
    identity: LeaseIdentity,
    parent_path: Option<&'static str>,
    logical_path: &'static str,
    capacity_bytes: u64,
    closed: bool,
    _parent: PhantomData<&'parent ()>,
}

impl fmt::Debug for DelegatedMemoryLease<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("DelegatedMemoryLease")
            .field("logical_path", &self.logical_path)
            .field("capacity_bytes", &self.capacity_bytes)
            .field("used_bytes", &self.used_bytes())
            .finish_non_exhaustive()
    }
}

impl DelegatedMemoryLease<'_> {
    /// Exact typed authority identity.
    #[must_use]
    pub fn identity(&self) -> LeaseIdentity {
        self.identity
    }

    /// Typed parent identity, or `None` for a direct root child.
    #[must_use]
    pub fn parent_identity(&self) -> Option<LeaseIdentity> {
        self.parent_identity
    }

    /// Exact full logical path.
    #[must_use]
    pub fn logical_path(&self) -> &'static str {
        self.logical_path
    }

    /// Parent path, or `None` for a direct root child.
    #[must_use]
    pub fn parent_path(&self) -> Option<&'static str> {
        self.parent_path
    }

    /// Capacity owned by this affine envelope.
    #[must_use]
    pub fn capacity_bytes(&self) -> u64 {
        self.capacity_bytes
    }

    /// Current child occupancy, including live nested envelopes.
    #[must_use]
    pub fn used_bytes(&self) -> u64 {
        let state = self.root.lock_state();
        state
            .delegations
            .iter()
            .find(|record| record.identity == self.identity)
            .map_or(0, |record| record.used_bytes)
    }

    /// Child occupied high-water.
    #[must_use]
    pub fn peak_used_bytes(&self) -> u64 {
        let state = self.root.lock_state();
        state
            .delegations
            .iter()
            .find(|record| record.identity == self.identity)
            .map_or(0, |record| record.peak_used_bytes)
    }

    /// Number of direct child envelopes still live.
    #[must_use]
    pub fn active_delegations(&self) -> u64 {
        let state = self.root.lock_state();
        state
            .delegations
            .iter()
            .find(|record| record.identity == self.identity)
            .map_or(0, |record| record.live_children)
    }

    /// Reserve allocation bytes inside this envelope.
    #[allow(clippy::result_large_err)]
    pub fn reserve<'lease>(
        &'lease self,
        site: &'static str,
        bytes: u64,
    ) -> Result<DelegatedLeaseCharge<'lease>, DelegatedLeaseRefusal> {
        reserve_delegated_charge(
            &self.root,
            self.identity,
            self.capacity_bytes,
            self.logical_path,
            site,
            bytes,
        )
    }

    /// Transfer capacity to a deterministic full descendant path.
    #[allow(clippy::result_large_err)]
    pub fn delegate_capacity<'child>(
        &'child self,
        identity: LeaseIdentity,
        logical_path: &'static str,
        capacity_bytes: u64,
    ) -> Result<DelegatedMemoryLease<'child>, LeaseDelegationRefusal> {
        self.root.delegate_from(
            Some(self.identity),
            Some(self.logical_path),
            identity,
            logical_path,
            capacity_bytes,
        )?;
        Ok(DelegatedMemoryLease {
            root: self.root.clone(),
            parent_identity: Some(self.identity),
            identity,
            parent_path: Some(self.logical_path),
            logical_path,
            capacity_bytes,
            closed: false,
            _parent: PhantomData,
        })
    }

    /// Explicitly return the envelope and obtain its verified close receipt.
    ///
    /// ```compile_fail
    /// use fs_alloc::{LeaseIdentity, OperationMemoryLease};
    ///
    /// let root = OperationMemoryLease::bounded(8);
    /// let root_id = LeaseIdentity::root(*b"example4", [1; 32]);
    /// let child_id = root_id.child([2; 32], 0).unwrap();
    /// root.enable_delegation(root_id, "run", 1).unwrap();
    /// let child = root.delegate_capacity(child_id, "run/child", 8).unwrap();
    /// let _receipt = child.close().unwrap();
    /// let _use_after_return = child.capacity_bytes();
    /// ```
    #[allow(clippy::result_large_err)]
    pub fn close(mut self) -> Result<DelegatedLeaseCloseReceipt, DelegatedLeaseCloseRefusal> {
        let result = return_delegation(&self.root, self.identity, false);
        self.closed = true;
        result
    }
}

impl Drop for DelegatedMemoryLease<'_> {
    fn drop(&mut self) {
        if !self.closed {
            let _ = return_delegation(&self.root, self.identity, true);
            self.closed = true;
        }
    }
}

/// Borrowed RAII charge inside one exact delegated owner.
#[must_use = "dropping the charge returns its reserved bytes"]
pub struct DelegatedLeaseCharge<'owner> {
    root: OperationMemoryLease,
    identity: LeaseIdentity,
    bytes: u64,
    active: bool,
    _owner: PhantomData<&'owner ()>,
}

impl fmt::Debug for DelegatedLeaseCharge<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("DelegatedLeaseCharge")
            .field("identity", &self.identity)
            .field("bytes", &self.bytes)
            .finish_non_exhaustive()
    }
}

impl<'owner> DelegatedLeaseCharge<'owner> {
    /// Bytes held by this charge.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Consume staging ownership and admit one exact prepared publication.
    ///
    /// Preparation retains the bytes inside the child envelope. Dropping the
    /// returned affine handle rolls those bytes back exactly once. If
    /// preparation is refused, consuming this charge returns its staging bytes
    /// through the ordinary charge-drop path.
    ///
    /// # Errors
    ///
    /// Refuses zero-byte, duplicate, post-seal, exhausted-metadata, unavailable
    /// child, and counter-overflow transitions without allocating.
    #[allow(clippy::result_large_err)]
    pub fn prepare_published_transfer(
        self,
        binding: PublishedTransferBinding,
    ) -> Result<PreparedPublishedTransfer<'owner>, PublishedTransferRefusal> {
        let envelope = PublishedTransferEnvelope::payload_only(self.bytes);
        self.prepare_published_transfer_with_envelope(binding, envelope)
    }

    /// Consume staging ownership and prepare an explicit byte composition.
    ///
    /// # Errors
    ///
    /// In addition to ordinary preparation refusals, the checked sum of
    /// payload, layout, and overhead must equal this charge exactly.
    #[allow(clippy::result_large_err)]
    pub fn prepare_published_transfer_with_envelope(
        mut self,
        binding: PublishedTransferBinding,
        envelope: PublishedTransferEnvelope,
    ) -> Result<PreparedPublishedTransfer<'owner>, PublishedTransferRefusal> {
        let prepared_sequence =
            prepare_published_transfer(&self.root, self.identity, self.bytes, binding, envelope)?;
        self.active = false;
        Ok(PreparedPublishedTransfer {
            root: self.root.clone(),
            child_identity: self.identity,
            binding,
            envelope,
            bytes: self.bytes,
            prepared_sequence,
            resolved: false,
            _owner: PhantomData,
        })
    }
}

impl Drop for DelegatedLeaseCharge<'_> {
    fn drop(&mut self) {
        if self.active {
            release_delegated(&self.root, self.identity, self.bytes);
            self.active = false;
        }
    }
}

/// Affine prepared publication that still owns staging bytes.
///
/// The handle borrows the delegated owner. It must resolve as either exactly
/// one successful publication or one rollback before that child can return.
///
/// ```compile_fail
/// use fs_alloc::{LeaseIdentity, OperationMemoryLease, PublishedTransferBinding};
///
/// let root = OperationMemoryLease::bounded(8);
/// let root_id = LeaseIdentity::root(*b"pubprep1", [1; 32]);
/// let child_id = root_id.child([2; 32], 0).unwrap();
/// root.enable_delegation(root_id, "run", 1).unwrap();
/// let child = root.delegate_capacity(child_id, "run/child", 8).unwrap();
/// let charge = child.reserve("output", 8).unwrap();
/// let prepared = charge
///     .prepare_published_transfer(PublishedTransferBinding::new(
///         [3; 32], [4; 32], [5; 32], [6; 32],
///     ))
///     .unwrap();
/// let _duplicate = prepared.clone();
/// ```
#[must_use = "prepared output must publish or roll back"]
pub struct PreparedPublishedTransfer<'owner> {
    root: OperationMemoryLease,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    prepared_sequence: u128,
    resolved: bool,
    _owner: PhantomData<&'owner ()>,
}

impl fmt::Debug for PreparedPublishedTransfer<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PreparedPublishedTransfer")
            .field("child_identity", &self.child_identity)
            .field("binding", &self.binding)
            .field("envelope", &self.envelope)
            .field("bytes", &self.bytes)
            .field("prepared_sequence", &self.prepared_sequence)
            .finish_non_exhaustive()
    }
}

impl PreparedPublishedTransfer<'_> {
    /// Exact staging child identity.
    #[must_use]
    pub fn child_identity(&self) -> LeaseIdentity {
        self.child_identity
    }

    /// Exact publication identity tuple.
    #[must_use]
    pub fn binding(&self) -> PublishedTransferBinding {
        self.binding
    }

    /// Exact payload/layout/overhead byte composition.
    #[must_use]
    pub fn envelope(&self) -> PublishedTransferEnvelope {
        self.envelope
    }

    /// Bytes retained in prepared staging.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Serialized preparation sequence.
    #[must_use]
    pub fn prepared_sequence(&self) -> u128 {
        self.prepared_sequence
    }

    /// Transfer staging ownership exactly once to the bound destination.
    ///
    /// A publish transition already prepared before the root admission cut may
    /// complete after sealing begins; it drains staging rather than admitting
    /// new work. On refusal, `Drop` attempts a fail-closed rollback.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal if the prepared record cannot be resolved
    /// or exact counters cannot advance.
    #[allow(clippy::result_large_err)]
    pub fn publish(mut self) -> Result<PublishedTransfer, PublishedTransferRefusal> {
        let receipt = publish_prepared_transfer(
            &self.root,
            self.child_identity,
            self.binding,
            self.envelope,
            self.bytes,
            self.prepared_sequence,
        )?;
        self.resolved = true;
        Ok(PublishedTransfer {
            root: self.root.clone(),
            receipt,
            closed: false,
        })
    }

    /// Explicitly return prepared staging bytes to the child.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal if the exact prepared record cannot be
    /// rolled back or exact counters cannot advance.
    #[allow(clippy::result_large_err)]
    pub fn rollback(
        mut self,
    ) -> Result<PublishedTransferRollbackReceipt, PublishedTransferRefusal> {
        let result = rollback_prepared_transfer(
            &self.root,
            self.child_identity,
            self.binding,
            self.envelope,
            self.bytes,
            self.prepared_sequence,
            false,
        );
        self.resolved = true;
        result
    }
}

impl Drop for PreparedPublishedTransfer<'_> {
    fn drop(&mut self) {
        if !self.resolved {
            let _ = rollback_prepared_transfer(
                &self.root,
                self.child_identity,
                self.binding,
                self.envelope,
                self.bytes,
                self.prepared_sequence,
                true,
            );
            self.resolved = true;
        }
    }
}

/// Affine ownership of bytes successfully transferred to one destination.
///
/// This handle is intentionally independent of the staging child lifetime.
/// The child may return and the root may seal while destination ownership
/// remains live. Destination `Drop` is recorded separately and never mutates a
/// previously issued root terminal receipt.
///
/// ```compile_fail
/// use fs_alloc::{LeaseIdentity, OperationMemoryLease, PublishedTransferBinding};
///
/// let root = OperationMemoryLease::bounded(8);
/// let root_id = LeaseIdentity::root(*b"pubdest1", [1; 32]);
/// let child_id = root_id.child([2; 32], 0).unwrap();
/// root.enable_delegation(root_id, "run", 1).unwrap();
/// let child = root.delegate_capacity(child_id, "run/child", 8).unwrap();
/// let published = child
///     .reserve("output", 8)
///     .unwrap()
///     .prepare_published_transfer(PublishedTransferBinding::new(
///         [3; 32], [4; 32], [5; 32], [6; 32],
///     ))
///     .unwrap()
///     .publish()
///     .unwrap();
/// let _duplicate = published.clone();
/// ```
#[must_use = "published destination ownership must be explicitly closed or dropped"]
pub struct PublishedTransfer {
    root: OperationMemoryLease,
    receipt: PublishedTransferReceipt,
    closed: bool,
}

impl fmt::Debug for PublishedTransfer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PublishedTransfer")
            .field("binding", &self.receipt.binding)
            .field("bytes", &self.receipt.bytes)
            .field("published_sequence", &self.receipt.published_sequence)
            .finish_non_exhaustive()
    }
}

impl PublishedTransfer {
    /// Verified successful publication receipt.
    #[must_use]
    pub fn receipt(&self) -> &PublishedTransferReceipt {
        &self.receipt
    }

    /// Destination-owned bytes.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.receipt.bytes
    }

    /// Exact destination identity tuple.
    #[must_use]
    pub fn binding(&self) -> PublishedTransferBinding {
        self.receipt.binding
    }

    /// Exact payload/layout/overhead byte composition.
    #[must_use]
    pub fn envelope(&self) -> PublishedTransferEnvelope {
        self.receipt.envelope
    }

    /// Explicitly close destination ownership exactly once.
    ///
    /// # Errors
    ///
    /// Returns a structured refusal if the exact published record is
    /// unavailable or its close counters cannot advance.
    #[allow(clippy::result_large_err)]
    pub fn close(mut self) -> Result<PublishedTransferCloseReceipt, PublishedTransferRefusal> {
        let result = close_published_transfer(&self.root, &self.receipt, false);
        self.closed = true;
        result
    }
}

impl Drop for PublishedTransfer {
    fn drop(&mut self) {
        if !self.closed {
            let _ = close_published_transfer(&self.root, &self.receipt, true);
            self.closed = true;
        }
    }
}

/// RAII root-lease charge: releases its bytes on drop (including unwinds).
#[derive(Debug)]
pub struct LeaseCharge {
    lease: OperationMemoryLease,
    bytes: u64,
}

impl LeaseCharge {
    /// Bytes held by this charge.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.bytes
    }

    /// Transfer this charge to a crate-internal owner that will release the
    /// same bytes manually.
    pub(crate) fn commit_to_manual_release(mut self) {
        self.bytes = 0;
    }
}

impl Drop for LeaseCharge {
    fn drop(&mut self) {
        if self.bytes > 0 {
            let _ = self.lease.release_raw(self.bytes);
        }
    }
}

#[allow(clippy::result_large_err)]
fn prepare_published_transfer(
    root: &OperationMemoryLease,
    child_identity: LeaseIdentity,
    bytes: u64,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
) -> Result<u128, PublishedTransferRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let root_identity = state
        .root_identity
        .expect("delegated charges have a configured root");
    let child_index = state
        .delegations
        .iter()
        .position(|record| record.identity == child_identity)
        .expect("delegated charge retains its child record");
    let child = &state.delegations[child_index];
    let parent_identity = child.parent_identity;
    let child_used = child.used_bytes;
    let reason = if state.sealed {
        Some(PublishedTransferRefusalReason::RootSealed)
    } else if child.disposition != DelegationDisposition::Live {
        Some(PublishedTransferRefusalReason::ChildReturned)
    } else if bytes == 0 {
        Some(PublishedTransferRefusalReason::ZeroBytes)
    } else if envelope.total_bytes() != Some(bytes) {
        Some(PublishedTransferRefusalReason::EnvelopeMismatch)
    } else if state
        .publications
        .iter()
        .any(|record| record.binding == binding)
    {
        Some(PublishedTransferRefusalReason::DuplicateBinding)
    } else if state.publications.len() >= state.metadata_limit {
        Some(PublishedTransferRefusalReason::MetadataExhausted)
    } else if child_used < bytes {
        Some(PublishedTransferRefusalReason::ConservationMismatch)
    } else if state.counter_overflowed {
        Some(PublishedTransferRefusalReason::CounterOverflow)
    } else {
        None
    };
    if let Some(reason) = reason {
        let refusal = record_published_transfer_refusal(
            &mut state,
            root_identity,
            parent_identity,
            child_identity,
            binding,
            envelope,
            bytes,
            PublishedTransferOperation::Prepare,
            reason,
            sequence,
        );
        return Err(refusal);
    }
    state.publications.push(PublicationRecord {
        parent_identity,
        child_identity,
        binding,
        envelope,
        bytes,
        prepared_sequence: sequence,
        resolved_sequence: None,
        disposition: PublicationDisposition::Prepared,
        implicit_rollback: false,
        published_receipt_root: None,
        destination_closed_sequence: None,
        implicit_destination_close: false,
        destination_close_root: None,
    });
    Ok(sequence)
}

// Publication is an allocation disposition, so its staging release and both
// sides of the granted = returned + published equation must linearize together.
#[allow(clippy::result_large_err, clippy::too_many_lines)]
fn publish_prepared_transfer(
    root: &OperationMemoryLease,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    prepared_sequence: u128,
) -> Result<PublishedTransferReceipt, PublishedTransferRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let root_identity = state
        .root_identity
        .expect("prepared transfers have a configured root");
    let child_index = state
        .delegations
        .iter()
        .position(|record| record.identity == child_identity)
        .expect("prepared transfer retains its child record");
    let parent_identity = state.delegations[child_index].parent_identity;
    let child_used = state.delegations[child_index].used_bytes;
    let publication_index = state.publications.iter().position(|record| {
        record.child_identity == child_identity
            && record.binding == binding
            && record.envelope == envelope
            && record.bytes == bytes
            && record.prepared_sequence == prepared_sequence
    });
    let reason = if state.delegations[child_index].disposition != DelegationDisposition::Live {
        Some(PublishedTransferRefusalReason::ChildReturned)
    } else if publication_index.is_none_or(|index| {
        state.publications[index].disposition != PublicationDisposition::Prepared
    }) {
        Some(PublishedTransferRefusalReason::TransferUnavailable)
    } else if child_used < bytes {
        Some(PublishedTransferRefusalReason::ConservationMismatch)
    } else if state.counter_overflowed
        || state.delegations[child_index]
            .allocation_published_bytes
            .checked_add(u128::from(bytes))
            .is_none()
        || state
            .child_published_bytes
            .checked_add(u128::from(bytes))
            .is_none()
    {
        Some(PublishedTransferRefusalReason::CounterOverflow)
    } else {
        None
    };
    if let Some(reason) = reason {
        let refusal = record_published_transfer_refusal(
            &mut state,
            root_identity,
            parent_identity,
            child_identity,
            binding,
            envelope,
            bytes,
            PublishedTransferOperation::Publish,
            reason,
            sequence,
        );
        return Err(refusal);
    }
    let publication_index = publication_index.expect("prepared record checked");
    let mut receipt = PublishedTransferReceipt {
        schema_version: PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION,
        root_identity,
        parent_identity,
        child_identity,
        binding,
        envelope,
        bytes,
        prepared_sequence,
        published_sequence: sequence,
        receipt_root: [0; 32],
    };
    receipt.receipt_root = receipt.recompute_root();
    receipt
        .verify_for(
            root_identity,
            parent_identity,
            child_identity,
            binding,
            envelope,
        )
        .expect("constructed publication receipt verifies");

    state.child_published_bytes += u128::from(bytes);
    let child = &mut state.delegations[child_index];
    child.used_bytes -= bytes;
    child.allocation_published_bytes += u128::from(bytes);
    let publication = &mut state.publications[publication_index];
    publication.disposition = PublicationDisposition::Published;
    publication.resolved_sequence = Some(sequence);
    publication.published_receipt_root = Some(receipt.receipt_root);
    Ok(receipt)
}

// Rollback is the alternate terminal outcome for one prepared record. The
// exact returned counters and prepared disposition change under one lock.
#[allow(clippy::result_large_err, clippy::too_many_lines)]
fn rollback_prepared_transfer(
    root: &OperationMemoryLease,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    prepared_sequence: u128,
    implicit_rollback: bool,
) -> Result<PublishedTransferRollbackReceipt, PublishedTransferRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let root_identity = state
        .root_identity
        .expect("prepared transfers have a configured root");
    let child_index = state
        .delegations
        .iter()
        .position(|record| record.identity == child_identity)
        .expect("prepared transfer retains its child record");
    let parent_identity = state.delegations[child_index].parent_identity;
    let child_used = state.delegations[child_index].used_bytes;
    let publication_index = state.publications.iter().position(|record| {
        record.child_identity == child_identity
            && record.binding == binding
            && record.envelope == envelope
            && record.bytes == bytes
            && record.prepared_sequence == prepared_sequence
    });
    let reason = if state.delegations[child_index].disposition != DelegationDisposition::Live {
        Some(PublishedTransferRefusalReason::ChildReturned)
    } else if publication_index.is_none_or(|index| {
        state.publications[index].disposition != PublicationDisposition::Prepared
    }) {
        Some(PublishedTransferRefusalReason::TransferUnavailable)
    } else if child_used < bytes {
        Some(PublishedTransferRefusalReason::ConservationMismatch)
    } else if state.counter_overflowed
        || state.delegations[child_index]
            .allocation_returned_bytes
            .checked_add(u128::from(bytes))
            .is_none()
        || state
            .child_returned_bytes
            .checked_add(u128::from(bytes))
            .is_none()
    {
        Some(PublishedTransferRefusalReason::CounterOverflow)
    } else {
        None
    };
    if let Some(reason) = reason {
        if implicit_rollback {
            increment_invariant_violation(&mut state);
        }
        let refusal = record_published_transfer_refusal(
            &mut state,
            root_identity,
            parent_identity,
            child_identity,
            binding,
            envelope,
            bytes,
            PublishedTransferOperation::Rollback,
            reason,
            sequence,
        );
        return Err(refusal);
    }
    let publication_index = publication_index.expect("prepared record checked");
    let mut receipt = PublishedTransferRollbackReceipt {
        schema_version: PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION,
        root_identity,
        parent_identity,
        child_identity,
        binding,
        envelope,
        bytes,
        prepared_sequence,
        rolled_back_sequence: sequence,
        implicit_rollback,
        receipt_root: [0; 32],
    };
    receipt.receipt_root = receipt.recompute_root();
    receipt
        .verify_for(
            root_identity,
            parent_identity,
            child_identity,
            binding,
            envelope,
        )
        .expect("constructed rollback receipt verifies");

    state.child_returned_bytes += u128::from(bytes);
    let child = &mut state.delegations[child_index];
    child.used_bytes -= bytes;
    child.allocation_returned_bytes += u128::from(bytes);
    let publication = &mut state.publications[publication_index];
    publication.disposition = PublicationDisposition::RolledBack;
    publication.resolved_sequence = Some(sequence);
    publication.implicit_rollback = implicit_rollback;
    Ok(receipt)
}

#[allow(clippy::result_large_err)]
fn close_published_transfer(
    root: &OperationMemoryLease,
    published: &PublishedTransferReceipt,
    implicit_close: bool,
) -> Result<PublishedTransferCloseReceipt, PublishedTransferRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let publication_index = state.publications.iter().position(|record| {
        record.child_identity == published.child_identity
            && record.binding == published.binding
            && record.envelope == published.envelope
            && record.bytes == published.bytes
            && record.prepared_sequence == published.prepared_sequence
    });
    let reason = if publication_index.is_none_or(|index| {
        let record = &state.publications[index];
        record.disposition != PublicationDisposition::Published
            || record.resolved_sequence != Some(published.published_sequence)
            || record.published_receipt_root != Some(published.receipt_root)
            || record.destination_closed_sequence.is_some()
    }) {
        Some(PublishedTransferRefusalReason::TransferUnavailable)
    } else if state.counter_overflowed {
        Some(PublishedTransferRefusalReason::CounterOverflow)
    } else {
        None
    };
    if let Some(reason) = reason {
        if implicit_close {
            increment_invariant_violation(&mut state);
        }
        let refusal = record_published_transfer_refusal(
            &mut state,
            published.root_identity,
            published.parent_identity,
            published.child_identity,
            published.binding,
            published.envelope,
            published.bytes,
            PublishedTransferOperation::CloseDestination,
            reason,
            sequence,
        );
        return Err(refusal);
    }
    let publication_index = publication_index.expect("published record checked");
    let mut receipt = PublishedTransferCloseReceipt {
        schema_version: PUBLISHED_TRANSFER_RECEIPT_SCHEMA_VERSION,
        root_identity: published.root_identity,
        parent_identity: published.parent_identity,
        child_identity: published.child_identity,
        binding: published.binding,
        envelope: published.envelope,
        bytes: published.bytes,
        published_receipt_root: published.receipt_root,
        published_sequence: published.published_sequence,
        closed_sequence: sequence,
        implicit_close,
        receipt_root: [0; 32],
    };
    receipt.receipt_root = receipt.recompute_root();
    receipt
        .verify_for(published)
        .expect("constructed destination close receipt verifies");
    let publication = &mut state.publications[publication_index];
    publication.destination_closed_sequence = Some(sequence);
    publication.implicit_destination_close = implicit_close;
    publication.destination_close_root = Some(receipt.receipt_root);
    Ok(receipt)
}

/// One serialized delegated-reservation admission, shared by the owner
/// method and typed restaging after a rollback. Keeping the refusal table,
/// refusal bookkeeping, and grant counters under one lock here means the
/// owner method and the typed layer cannot drift.
#[allow(clippy::result_large_err)]
fn reserve_delegated_charge<'owner>(
    root: &OperationMemoryLease,
    identity: LeaseIdentity,
    capacity_bytes: u64,
    logical_path: &'static str,
    site: &'static str,
    bytes: u64,
) -> Result<DelegatedLeaseCharge<'owner>, DelegatedLeaseRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let root_identity = state
        .root_identity
        .expect("delegated roots have typed identities");
    let root_id = state.root_id.expect("delegated roots are configured");
    let index = state
        .delegations
        .iter()
        .position(|record| record.identity == identity)
        .expect("delegated record retained for root lifetime");
    let record = &state.delegations[index];
    let used = record.used_bytes;
    let reason = if state.sealed {
        Some(DelegatedReservationRefusalReason::RootSealed)
    } else if record.disposition == DelegationDisposition::Returned {
        Some(DelegatedReservationRefusalReason::ChildReturned)
    } else if state.counter_overflowed
        || record
            .allocation_granted_bytes
            .checked_add(u128::from(bytes))
            .is_none()
        || state
            .child_granted_bytes
            .checked_add(u128::from(bytes))
            .is_none()
    {
        Some(DelegatedReservationRefusalReason::CounterOverflow)
    } else if used
        .checked_add(bytes)
        .is_none_or(|next| next > capacity_bytes)
    {
        Some(DelegatedReservationRefusalReason::Capacity)
    } else {
        None
    };
    if let Some(reason) = reason {
        if reason == DelegatedReservationRefusalReason::CounterOverflow {
            state.counter_overflowed = true;
        }
        let next_refused_requests = state.delegations[index].refused_requests.checked_add(1);
        let next_refused_bytes = state.delegations[index]
            .refused_bytes
            .checked_add(u128::from(bytes));
        if let (Some(requests), Some(refused_bytes)) = (next_refused_requests, next_refused_bytes) {
            let record = &mut state.delegations[index];
            record.refused_requests = requests;
            record.refused_bytes = refused_bytes;
        } else {
            state.counter_overflowed = true;
        }
        record_semantic_refusal(
            &mut state,
            reason_code_delegated_reservation(reason),
            Some(identity),
            site,
            Some(bytes),
            used,
            capacity_bytes,
            sequence,
        );
        return Err(DelegatedLeaseRefusal {
            root_identity,
            identity,
            root_id,
            logical_path,
            site,
            requested_bytes: bytes,
            used_bytes: used,
            limit_bytes: capacity_bytes,
            reason,
            sequence,
        });
    }
    state.child_granted_bytes += u128::from(bytes);
    let record = &mut state.delegations[index];
    record.used_bytes += bytes;
    record.peak_used_bytes = record.peak_used_bytes.max(record.used_bytes);
    record.allocation_granted_bytes += u128::from(bytes);
    drop(state);
    Ok(DelegatedLeaseCharge {
        root: root.clone(),
        identity,
        bytes,
        active: true,
        _owner: PhantomData,
    })
}

#[allow(clippy::too_many_arguments)]
fn record_published_transfer_refusal(
    state: &mut LeaseState,
    root_identity: LeaseIdentity,
    parent_identity: Option<LeaseIdentity>,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    operation: PublishedTransferOperation,
    reason: PublishedTransferRefusalReason,
    sequence: u128,
) -> PublishedTransferRefusal {
    let observed_used = state
        .delegations
        .iter()
        .find(|record| record.identity == child_identity)
        .map_or(0, |record| record.used_bytes);
    let publication_count = u64::try_from(state.publications.len()).unwrap_or(u64::MAX);
    record_semantic_refusal(
        state,
        reason_code_published_transfer(reason),
        Some(child_identity),
        operation_code_published_transfer(operation),
        Some(bytes),
        observed_used,
        publication_count,
        sequence,
    );
    PublishedTransferRefusal {
        root_identity,
        parent_identity,
        child_identity,
        binding,
        envelope,
        bytes,
        operation,
        reason,
        sequence,
    }
}

// This is the exact inverse transition for delegation admission. Keeping its
// conservation checks and receipt snapshot under one lock avoids stale facts.
#[allow(clippy::result_large_err, clippy::too_many_lines)]
fn return_delegation(
    root: &OperationMemoryLease,
    identity: LeaseIdentity,
    implicit_return: bool,
) -> Result<DelegatedLeaseCloseReceipt, DelegatedLeaseCloseRefusal> {
    let mut state = root.lock_state();
    let sequence = next_sequence(&mut state);
    let root_identity = state
        .root_identity
        .expect("delegated roots have typed identities");
    let root_id = state.root_id.expect("delegated roots are configured");
    let index = state
        .delegations
        .iter()
        .position(|record| record.identity == identity)
        .expect("delegated record retained for root lifetime");
    let record = &state.delegations[index];
    let reason = if record.disposition == DelegationDisposition::Returned {
        Some(CloseRefusalReason::AlreadyReturned)
    } else if record.live_children != 0 {
        Some(CloseRefusalReason::LiveChild)
    } else if record.used_bytes != 0 {
        Some(CloseRefusalReason::LiveAllocation)
    } else if record
        .allocation_returned_bytes
        .checked_add(record.allocation_published_bytes)
        != Some(record.allocation_granted_bytes)
        || record.delegated_bytes != record.returned_delegated_bytes
    {
        Some(CloseRefusalReason::ConservationMismatch)
    } else if state.counter_overflowed {
        Some(CloseRefusalReason::CounterOverflow)
    } else {
        None
    };
    if let Some(reason) = reason {
        let used_bytes = state.delegations[index].used_bytes;
        let live_children = state.delegations[index].live_children;
        let logical_path = state.delegations[index].logical_path;
        if implicit_return {
            increment_invariant_violation(&mut state);
        }
        record_semantic_refusal(
            &mut state,
            reason_code_close(reason),
            Some(identity),
            logical_path,
            None,
            used_bytes,
            live_children,
            sequence,
        );
        return Err(DelegatedLeaseCloseRefusal {
            identity,
            logical_path,
            reason,
            used_bytes,
            live_children,
            sequence,
        });
    }

    state.publications.sort_unstable_by_key(|record| {
        (
            record.binding,
            record.child_identity,
            record.prepared_sequence,
        )
    });
    let snapshot = state.delegations[index].clone();
    let publication_record_count = state
        .publications
        .iter()
        .filter(|record| record.child_identity == identity)
        .count();
    let published_transfer_count = state
        .publications
        .iter()
        .filter(|record| {
            record.child_identity == identity
                && record.disposition == PublicationDisposition::Published
        })
        .count();
    let rolled_back_transfer_count = state
        .publications
        .iter()
        .filter(|record| {
            record.child_identity == identity
                && record.disposition == PublicationDisposition::RolledBack
        })
        .count();
    if published_transfer_count.checked_add(rolled_back_transfer_count)
        != Some(publication_record_count)
    {
        increment_invariant_violation(&mut state);
        return Err(DelegatedLeaseCloseRefusal {
            identity,
            logical_path: snapshot.logical_path,
            reason: CloseRefusalReason::ConservationMismatch,
            used_bytes: snapshot.used_bytes,
            live_children: snapshot.live_children,
            sequence,
        });
    }
    let publication_root = publication_root_for_child(&state.publications, identity);
    let capacity = snapshot.capacity_bytes;
    let parent_identity = snapshot.parent_identity;
    if let Some(parent_identity) = parent_identity {
        let parent_index = state
            .delegations
            .iter()
            .position(|record| record.identity == parent_identity)
            .expect("parent record retained");
        let parent = &state.delegations[parent_index];
        if parent.used_bytes < capacity
            || parent.live_children == 0
            || parent
                .returned_delegated_bytes
                .checked_add(u128::from(capacity))
                .is_none()
        {
            increment_invariant_violation(&mut state);
            return Err(DelegatedLeaseCloseRefusal {
                identity,
                logical_path: snapshot.logical_path,
                reason: CloseRefusalReason::ConservationMismatch,
                used_bytes: snapshot.used_bytes,
                live_children: snapshot.live_children,
                sequence,
            });
        }
    } else if state.used_bytes < capacity
        || state
            .returned_delegated_bytes
            .checked_add(u128::from(capacity))
            .is_none()
    {
        increment_invariant_violation(&mut state);
        return Err(DelegatedLeaseCloseRefusal {
            identity,
            logical_path: snapshot.logical_path,
            reason: CloseRefusalReason::ConservationMismatch,
            used_bytes: snapshot.used_bytes,
            live_children: snapshot.live_children,
            sequence,
        });
    }

    if let Some(parent_identity) = parent_identity {
        let parent = state
            .delegations
            .iter_mut()
            .find(|record| record.identity == parent_identity)
            .expect("validated parent");
        parent.used_bytes -= capacity;
        parent.live_children -= 1;
        parent.returned_delegated_bytes += u128::from(capacity);
    } else {
        state.used_bytes -= capacity;
        state.returned_delegated_bytes += u128::from(capacity);
    }

    let mut receipt = DelegatedLeaseCloseReceipt {
        schema_version: VERIFIED_RECEIPT_SCHEMA_VERSION,
        root_identity,
        parent_identity: snapshot.parent_identity,
        identity: snapshot.identity,
        root_id,
        parent_path: snapshot.parent_path,
        logical_path: snapshot.logical_path,
        capacity_bytes: snapshot.capacity_bytes,
        allocation_granted_bytes: snapshot.allocation_granted_bytes,
        allocation_returned_bytes: snapshot.allocation_returned_bytes,
        allocation_published_bytes: snapshot.allocation_published_bytes,
        delegated_bytes: snapshot.delegated_bytes,
        returned_delegated_bytes: snapshot.returned_delegated_bytes,
        peak_used_bytes: snapshot.peak_used_bytes,
        final_used_bytes: snapshot.used_bytes,
        refused_requests: snapshot.refused_requests,
        refused_bytes: snapshot.refused_bytes,
        publication_record_count,
        published_transfer_count,
        rolled_back_transfer_count,
        publication_root,
        created_sequence: snapshot.created_sequence,
        returned_sequence: sequence,
        implicit_return,
        receipt_root: [0; 32],
    };
    receipt.receipt_root = receipt.recompute_root();
    receipt
        .verify_for(root_identity, snapshot.parent_identity, identity)
        .expect("constructed child receipt verifies");
    let record = &mut state.delegations[index];
    record.disposition = DelegationDisposition::Returned;
    record.returned_sequence = Some(sequence);
    record.implicit_return = implicit_return;
    record.close_root = Some(receipt.receipt_root);
    Ok(receipt)
}

fn release_delegated(root: &OperationMemoryLease, identity: LeaseIdentity, bytes: u64) {
    let mut state = root.lock_state();
    let _sequence = next_sequence(&mut state);
    let Some(index) = state
        .delegations
        .iter()
        .position(|record| record.identity == identity)
    else {
        increment_invariant_violation(&mut state);
        return;
    };
    let record = &state.delegations[index];
    if record.disposition != DelegationDisposition::Live || record.used_bytes < bytes {
        increment_invariant_violation(&mut state);
        return;
    }
    let Some(record_returned) = record
        .allocation_returned_bytes
        .checked_add(u128::from(bytes))
    else {
        state.counter_overflowed = true;
        increment_invariant_violation(&mut state);
        return;
    };
    let Some(global_returned) = state.child_returned_bytes.checked_add(u128::from(bytes)) else {
        state.counter_overflowed = true;
        increment_invariant_violation(&mut state);
        return;
    };
    state.child_returned_bytes = global_returned;
    let record = &mut state.delegations[index];
    record.used_bytes -= bytes;
    record.allocation_returned_bytes = record_returned;
}

fn next_sequence(state: &mut LeaseState) -> u128 {
    if let Some(next) = state.sequence.checked_add(1) {
        state.sequence = next;
        next
    } else {
        state.counter_overflowed = true;
        u128::MAX
    }
}

fn increment_invariant_violation(state: &mut LeaseState) {
    if let Some(next) = state.release_invariant_violations.checked_add(1) {
        state.release_invariant_violations = next;
    } else {
        state.counter_overflowed = true;
    }
}

fn record_root_refusal(
    state: &mut LeaseState,
    limit_bytes: Option<u64>,
    what: &'static str,
    bytes: u64,
    used: u64,
    reason: LeaseRefusalReason,
    sequence: u128,
) -> LeaseRefusal {
    let root_identity = state.root_identity;
    let refusal = LeaseRefusal {
        what,
        requested_bytes: bytes,
        used_bytes: used,
        limit_bytes: limit_bytes.unwrap_or(u64::MAX),
        reason,
        sequence,
    };
    if let Some(next) = state.refusals.checked_add(1) {
        state.refusals = next;
    } else {
        state.counter_overflowed = true;
    }
    if state.first_refusal.is_none() {
        state.first_refusal = Some(refusal.clone());
    }
    record_semantic_refusal(
        state,
        reason_code_root_reservation(reason),
        root_identity,
        what,
        Some(bytes),
        used,
        limit_bytes.unwrap_or(u64::MAX),
        sequence,
    );
    refusal
}

// The arguments are the fixed canonical refusal-event fields; grouping them in
// an allocated container would weaken the fail-closed refusal path.
#[allow(clippy::too_many_arguments)]
fn record_semantic_refusal(
    state: &mut LeaseState,
    reason: &'static str,
    identity: Option<LeaseIdentity>,
    site: &'static str,
    requested_bytes: Option<u64>,
    observed_used_bytes: u64,
    observed_auxiliary: u64,
    sequence: u128,
) {
    if state.terminal_receipt.is_some() {
        if let Some(next) = state.post_terminal_refusals.checked_add(1) {
            state.post_terminal_refusals = next;
        } else {
            state.counter_overflowed = true;
        }
        return;
    }
    if let Some(next) = state.semantic_refusals.checked_add(1) {
        state.semantic_refusals = next;
    } else {
        state.counter_overflowed = true;
    }
    if let Some(requested_bytes) = requested_bytes {
        if let Some(next) = state
            .semantic_refused_bytes
            .checked_add(u128::from(requested_bytes))
        {
            state.semantic_refused_bytes = next;
        } else {
            state.counter_overflowed = true;
        }
    }
    let mut hash = ReceiptHasher::new(b"fs-alloc-refusal-event-v1");
    hash.bytes(&state.refusal_root);
    hash.text(reason);
    hash.optional_identity(identity);
    hash.text(site);
    match requested_bytes {
        Some(requested_bytes) => {
            hash.boolean(true);
            hash.u64(requested_bytes);
        }
        None => hash.boolean(false),
    }
    hash.u64(observed_used_bytes);
    hash.u64(observed_auxiliary);
    hash.u128(sequence);
    state.refusal_root = hash.finish();
}

fn count_live_delegations(state: &LeaseState) -> usize {
    state
        .delegations
        .iter()
        .filter(|record| record.disposition == DelegationDisposition::Live)
        .count()
}

fn lease_receipt_locked(limit_bytes: Option<u64>, state: &LeaseState) -> LeaseReceipt {
    LeaseReceipt {
        limit_bytes,
        requested_bytes: state.requested_bytes,
        peak_bytes: state.peak_bytes,
        used_bytes: state.used_bytes,
        refusals: state.refusals,
        release_invariant_violations: state.release_invariant_violations,
        first_refusal: state.first_refusal.clone(),
    }
}

fn valid_root_id(root_id: &str) -> bool {
    !root_id.is_empty()
        && root_id.len() <= MAX_LOGICAL_ID_BYTES
        && !root_id.contains('/')
        && root_id.bytes().all(valid_identity_byte)
}

fn valid_logical_path(root_id: &str, parent_path: Option<&str>, path: &str) -> bool {
    if path.is_empty()
        || path.len() > MAX_LOGICAL_ID_BYTES
        || path.starts_with('/')
        || path.ends_with('/')
        || path.split('/').any(str::is_empty)
        || !path
            .bytes()
            .all(|byte| byte == b'/' || valid_identity_byte(byte))
    {
        return false;
    }
    let prefix = parent_path.unwrap_or(root_id);
    path.strip_prefix(prefix)
        .is_some_and(|suffix| suffix.starts_with('/') && suffix.len() > 1)
}

fn valid_identity_byte(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-')
}

fn delegation_root(records: &[DelegationRecord]) -> [u8; 32] {
    let mut hash = ReceiptHasher::new(b"fs-alloc-delegation-ledger-v2");
    hash.usize(records.len());
    for record in records {
        hash.optional_identity(record.parent_identity);
        hash.identity(record.identity);
        hash.optional_text(record.parent_path);
        hash.text(record.logical_path);
        hash.u64(record.capacity_bytes);
        hash.u64(record.peak_used_bytes);
        hash.u128(record.allocation_granted_bytes);
        hash.u128(record.allocation_returned_bytes);
        hash.u128(record.allocation_published_bytes);
        hash.u128(record.delegated_bytes);
        hash.u128(record.returned_delegated_bytes);
        hash.u128(record.refused_requests);
        hash.u128(record.refused_bytes);
        hash.u128(record.created_sequence);
        hash.u128(record.returned_sequence.unwrap_or(0));
        hash.boolean(record.implicit_return);
        hash.bytes(&record.close_root.unwrap_or([0; 32]));
    }
    hash.finish()
}

fn publication_root(records: &[PublicationRecord]) -> [u8; 32] {
    publication_root_filtered(records, None)
}

fn publication_root_for_child(
    records: &[PublicationRecord],
    child_identity: LeaseIdentity,
) -> [u8; 32] {
    publication_root_filtered(records, Some(child_identity))
}

fn publication_root_filtered(
    records: &[PublicationRecord],
    child_identity: Option<LeaseIdentity>,
) -> [u8; 32] {
    let mut hash = ReceiptHasher::new(b"fs-alloc-publication-ledger-v1");
    let count = records
        .iter()
        .filter(|record| child_identity.is_none_or(|identity| record.child_identity == identity))
        .count();
    hash.usize(count);
    for record in records
        .iter()
        .filter(|record| child_identity.is_none_or(|identity| record.child_identity == identity))
    {
        hash.optional_identity(record.parent_identity);
        hash.identity(record.child_identity);
        hash.published_binding(record.binding);
        hash.published_envelope(record.envelope);
        hash.u64(record.bytes);
        hash.u128(record.prepared_sequence);
        hash.u128(record.resolved_sequence.unwrap_or(0));
        hash.u16(match record.disposition {
            PublicationDisposition::Prepared => 0,
            PublicationDisposition::Published => 1,
            PublicationDisposition::RolledBack => 2,
        });
        hash.boolean(record.implicit_rollback);
        hash.bytes(&record.published_receipt_root.unwrap_or([0; 32]));
    }
    hash.finish()
}

fn reason_code_configuration(reason: ConfigurationRefusalReason) -> &'static str {
    match reason {
        ConfigurationRefusalReason::UnboundedRoot => "unbounded_root",
        ConfigurationRefusalReason::NotPristine => "root_not_pristine",
        ConfigurationRefusalReason::AlreadyConfigured => "already_configured",
        ConfigurationRefusalReason::InvalidRootIdentity => "invalid_root_identity",
        ConfigurationRefusalReason::MetadataLimit => "metadata_limit",
        ConfigurationRefusalReason::MetadataAllocation => "metadata_allocation",
        ConfigurationRefusalReason::Sealed => "sealed",
    }
}

fn reason_code_delegation(reason: DelegationRefusalReason) -> &'static str {
    match reason {
        DelegationRefusalReason::UnboundedParent => "unbounded_parent",
        DelegationRefusalReason::UnconfiguredRoot => "unconfigured_root",
        DelegationRefusalReason::RootSealed => "root_sealed",
        DelegationRefusalReason::InvalidIdentityRelationship => "invalid_identity_relationship",
        DelegationRefusalReason::InvalidLogicalPath => "invalid_logical_path",
        DelegationRefusalReason::DuplicateIdentity => "duplicate_identity",
        DelegationRefusalReason::DuplicatePath => "duplicate_path",
        DelegationRefusalReason::MetadataExhausted => "metadata_exhausted",
        DelegationRefusalReason::Capacity => "capacity",
        DelegationRefusalReason::ParentReturned => "parent_returned",
        DelegationRefusalReason::CounterOverflow => "counter_overflow",
    }
}

fn reason_code_delegated_reservation(reason: DelegatedReservationRefusalReason) -> &'static str {
    match reason {
        DelegatedReservationRefusalReason::RootSealed => "root_sealed",
        DelegatedReservationRefusalReason::ChildReturned => "child_returned",
        DelegatedReservationRefusalReason::Capacity => "capacity",
        DelegatedReservationRefusalReason::CounterOverflow => "counter_overflow",
    }
}

fn reason_code_published_transfer(reason: PublishedTransferRefusalReason) -> &'static str {
    match reason {
        PublishedTransferRefusalReason::RootSealed => "root_sealed",
        PublishedTransferRefusalReason::ChildReturned => "child_returned",
        PublishedTransferRefusalReason::ZeroBytes => "zero_bytes",
        PublishedTransferRefusalReason::EnvelopeMismatch => "envelope_mismatch",
        PublishedTransferRefusalReason::DuplicateBinding => "duplicate_binding",
        PublishedTransferRefusalReason::MetadataExhausted => "metadata_exhausted",
        PublishedTransferRefusalReason::TransferUnavailable => "transfer_unavailable",
        PublishedTransferRefusalReason::ConservationMismatch => "conservation_mismatch",
        PublishedTransferRefusalReason::CounterOverflow => "counter_overflow",
    }
}

fn operation_code_published_transfer(operation: PublishedTransferOperation) -> &'static str {
    match operation {
        PublishedTransferOperation::Prepare => "prepare_published_transfer",
        PublishedTransferOperation::Publish => "publish_prepared_transfer",
        PublishedTransferOperation::Rollback => "rollback_prepared_transfer",
        PublishedTransferOperation::CloseDestination => "close_published_destination",
    }
}

fn reason_code_close(reason: CloseRefusalReason) -> &'static str {
    match reason {
        CloseRefusalReason::LiveAllocation => "live_allocation",
        CloseRefusalReason::LiveChild => "live_child",
        CloseRefusalReason::AlreadyReturned => "already_returned",
        CloseRefusalReason::ConservationMismatch => "conservation_mismatch",
        CloseRefusalReason::CounterOverflow => "counter_overflow",
    }
}

fn reason_code_seal(reason: SealRefusalReason) -> &'static str {
    match reason {
        SealRefusalReason::UnverifiedRoot => "unverified_root",
        SealRefusalReason::LiveCapacity => "live_capacity",
        SealRefusalReason::ConservationMismatch => "conservation_mismatch",
        SealRefusalReason::ReleaseInvariant => "release_invariant",
        SealRefusalReason::CounterOverflow => "counter_overflow",
    }
}

fn reason_code_root_reservation(reason: LeaseRefusalReason) -> &'static str {
    match reason {
        LeaseRefusalReason::Capacity => "capacity",
        LeaseRefusalReason::Sealed => "sealed",
        LeaseRefusalReason::CounterOverflow => "counter_overflow",
    }
}

// ---------------------------------------------------------------------------
// Typed two-phase publication guard (bead 6ys.21.1.3.2)
//
// The byte-level transfer surface above proves the ledger transitions; this
// layer binds exactly one live `T` to the same authority so a committed value
// can never escape beside raw bytes. Refusals hand staging and value back
// UNCHANGED by calling the module-free transition functions directly instead
// of the consuming byte-level wrappers, whose error paths release staging.

/// Charged authority bytes for one `T`.
///
/// A zero-sized value still occupies one counted authority byte: every
/// publication remains a counted authority even when the payload occupies
/// none, and the ledger's zero-byte preparation refusal stays intact.
#[must_use]
pub fn authority_bytes_for<T>() -> u64 {
    size_of::<T>().max(1) as u64
}

/// Allocator-originated staging for exactly one live `T`.
///
/// Constructed only through [`DelegatedMemoryLease::allocate`]; the charged
/// staging envelope is exactly [`authority_bytes_for::<T>`] declared as pure
/// payload. Dropping the allocation releases the charge and destroys the
/// value normally; preparing consumes it into an affine two-phase handle.
pub struct LeasedAllocation<'owner, T> {
    inner: Option<StagedAllocation<'owner, T>>,
}

struct StagedAllocation<'owner, T> {
    charge: DelegatedLeaseCharge<'owner>,
    capacity_bytes: u64,
    logical_path: &'static str,
    value: T,
}

impl<T: fmt::Debug> fmt::Debug for LeasedAllocation<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.inner.as_ref() {
            Some(staged) => f
                .debug_struct("LeasedAllocation")
                .field("bytes", &staged.charge.bytes)
                .field("value", &staged.value)
                .finish_non_exhaustive(),
            None => f.debug_struct("LeasedAllocation").finish_non_exhaustive(),
        }
    }
}

impl<'owner> DelegatedMemoryLease<'owner> {
    /// Reserve exact staging authority for and bind one live `T`.
    ///
    /// The charged envelope is [`authority_bytes_for::<T>`] as pure payload;
    /// a zero-sized value charges one counted authority byte. On reservation
    /// refusal the candidate value is dropped — reserve capacity first when
    /// losing the value would matter.
    ///
    /// # Errors
    ///
    /// Returns the ordinary delegated-reservation refusal unchanged.
    #[allow(clippy::result_large_err)]
    pub fn allocate<T>(
        &'owner self,
        site: &'static str,
        value: T,
    ) -> Result<LeasedAllocation<'owner, T>, DelegatedLeaseRefusal> {
        let charge = self.reserve(site, authority_bytes_for::<T>())?;
        Ok(LeasedAllocation {
            inner: Some(StagedAllocation {
                charge,
                capacity_bytes: self.capacity_bytes,
                logical_path: self.logical_path,
                value,
            }),
        })
    }
}

impl<'owner, T> LeasedAllocation<'owner, T> {
    /// Exact charged staging bytes.
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.inner.as_ref().map_or(0, |staged| staged.charge.bytes)
    }

    /// Consume staging ownership and admit one prepared publication of the
    /// bound value under one exact binding tuple.
    ///
    /// Attempt generations are caller discipline: after a non-mutating
    /// rollback the same value and staging may be re-prepared under a fresh
    /// binding (typically a regenerated occurrence identity); the ledger's
    /// duplicate-binding refusal rejects reused tuples deterministically.
    ///
    /// # Errors
    ///
    /// Every preparation refusal returns the allocation UNCHANGED — staging
    /// charge still active, value untouched — inside the rejection, so a
    /// corrected retry loses nothing.
    #[allow(clippy::result_large_err)]
    pub fn prepare(
        mut self,
        binding: PublishedTransferBinding,
    ) -> Result<TypedPreparedPublication<'owner, T>, TypedPrepareRejection<'owner, T>> {
        let mut staged = self
            .inner
            .take()
            .expect("affine staging consumed at most once");
        let bytes = staged.charge.bytes;
        let envelope = PublishedTransferEnvelope::payload_only(bytes);
        let prepared_sequence = match prepare_published_transfer(
            &staged.charge.root,
            staged.charge.identity,
            bytes,
            binding,
            envelope,
        ) {
            Ok(prepared_sequence) => prepared_sequence,
            Err(refusal) => {
                return Err(TypedPrepareRejection {
                    allocation: Self {
                        inner: Some(staged),
                    },
                    refusal,
                });
            }
        };
        // The prepared record owns the staging disposition now; the charge
        // must not release its bytes on drop.
        staged.charge.active = false;
        Ok(TypedPreparedPublication {
            root: staged.charge.root.clone(),
            child_identity: staged.charge.identity,
            binding,
            envelope,
            bytes,
            child_capacity: staged.capacity_bytes,
            child_path: staged.logical_path,
            prepared_sequence,
            resolved: false,
            value: Some(staged.value),
            _owner: PhantomData,
        })
    }
}
/// Affine prepared publication that carries the staged value itself.
///
/// Must resolve as exactly one successful commit or one rollback; dropping
/// without resolution performs the fail-closed implicit rollback (with the
/// same invariant-violation accounting as the byte-level surface) and
/// destroys the value.
#[must_use = "prepared output must commit or roll back"]
pub struct TypedPreparedPublication<'owner, T> {
    root: OperationMemoryLease,
    child_identity: LeaseIdentity,
    binding: PublishedTransferBinding,
    envelope: PublishedTransferEnvelope,
    bytes: u64,
    child_capacity: u64,
    child_path: &'static str,
    prepared_sequence: u128,
    resolved: bool,
    value: Option<T>,
    _owner: PhantomData<&'owner ()>,
}

impl<T> fmt::Debug for TypedPreparedPublication<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TypedPreparedPublication")
            .field("child_identity", &self.child_identity)
            .field("binding", &self.binding)
            .field("envelope", &self.envelope)
            .field("bytes", &self.bytes)
            .field("prepared_sequence", &self.prepared_sequence)
            .finish_non_exhaustive()
    }
}

impl<'owner, T> TypedPreparedPublication<'owner, T> {
    /// Exact publication identity tuple.
    #[must_use]
    pub fn binding(&self) -> PublishedTransferBinding {
        self.binding
    }

    /// Serialized preparation sequence.
    #[must_use]
    pub fn prepared_sequence(&self) -> u128 {
        self.prepared_sequence
    }

    /// Commit exactly once: transfer staging authority to the destination
    /// and bind the live value into the resulting published allocation.
    ///
    /// A publish already prepared before the seal cut may complete during
    /// drain rather than admitting new work.
    ///
    /// # Errors
    ///
    /// Every refusal returns the prepared publication UNCHANGED — record
    /// still prepared, value intact — inside the rejection, so commit may be
    /// retried after diagnosing or replaced by rollback.
    #[allow(clippy::result_large_err)]
    pub fn commit(mut self) -> Result<PublishedAllocation<T>, TypedCommitRejection<'owner, T>> {
        let value = self
            .value
            .take()
            .expect("carried value lives until resolution");
        match publish_prepared_transfer(
            &self.root,
            self.child_identity,
            self.binding,
            self.envelope,
            self.bytes,
            self.prepared_sequence,
        ) {
            Ok(receipt) => {
                self.resolved = true;
                Ok(PublishedAllocation {
                    root: self.root.clone(),
                    receipt,
                    value: Some(value),
                    destroyed: false,
                    closed: false,
                })
            }
            Err(refusal) => {
                self.value = Some(value);
                Err(TypedCommitRejection {
                    prepared: self,
                    refusal,
                })
            }
        }
    }

    /// Roll back exactly once, then genuinely re-reserve the freed staging
    /// so the returned allocation carries live authority again — never a
    /// fabricated charge over bytes the ledger no longer counts.
    ///
    /// # Errors
    ///
    /// If the ledger refuses the rollback itself, the prepared publication
    /// is returned UNCHANGED inside the rejection. If the rollback succeeded
    /// but re-reservation was refused (for example a concurrent seal cut),
    /// the untouched value and the verified rollback receipt are returned
    /// instead; the staging authority is already back with the child pool.
    #[allow(clippy::result_large_err)]
    pub fn rollback(
        mut self,
    ) -> Result<TypedRollback<'owner, T>, TypedRollbackRejection<'owner, T>> {
        let value = self
            .value
            .take()
            .expect("carried value lives until resolution");
        let receipt = match rollback_prepared_transfer(
            &self.root,
            self.child_identity,
            self.binding,
            self.envelope,
            self.bytes,
            self.prepared_sequence,
            false,
        ) {
            Ok(receipt) => receipt,
            Err(refusal) => {
                self.value = Some(value);
                return Err(TypedRollbackRejection::Prepared(self, refusal));
            }
        };
        self.resolved = true;
        let charge = match reserve_delegated_charge(
            &self.root,
            self.child_identity,
            self.child_capacity,
            self.child_path,
            "restaged-output",
            self.bytes,
        ) {
            Ok(charge) => charge,
            Err(refusal) => {
                return Err(TypedRollbackRejection::Released(value, receipt, refusal));
            }
        };
        Ok(TypedRollback {
            allocation: LeasedAllocation {
                inner: Some(StagedAllocation {
                    charge,
                    capacity_bytes: self.child_capacity,
                    logical_path: self.child_path,
                    value,
                }),
            },
            receipt,
        })
    }
}

impl<T> Drop for TypedPreparedPublication<'_, T> {
    fn drop(&mut self) {
        drop(self.value.take());
        if !self.resolved {
            let _ = rollback_prepared_transfer(
                &self.root,
                self.child_identity,
                self.binding,
                self.envelope,
                self.bytes,
                self.prepared_sequence,
                true,
            );
            self.resolved = true;
        }
    }
}

/// Affine ownership of one successfully published value.
///
/// Supports shared observation and consuming authority-preserving handoff
/// only: no raw parts, no mutable access, no growth, no detached value
/// serialization, no clone. Close destroys the value first and only then
/// records the destination close, so a destructor panic can never be
/// promoted to a successful close.
///
/// ```compile_fail
/// use fs_alloc::{LeaseIdentity, OperationMemoryLease, PublishedAllocation};
///
/// let forged = PublishedAllocation {
///     root_identity_placeholder: (),
/// };
/// ```
#[must_use = "published ownership must be explicitly closed or dropped"]
pub struct PublishedAllocation<T> {
    root: OperationMemoryLease,
    receipt: PublishedTransferReceipt,
    value: Option<T>,
    /// Set only after value destruction completed WITHOUT panicking; a
    /// destructor panic leaves it false so no close is ever recorded.
    destroyed: bool,
    closed: bool,
}

impl<T> fmt::Debug for PublishedAllocation<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PublishedAllocation")
            .field("binding", &self.receipt.binding)
            .field("bytes", &self.receipt.bytes)
            .field("published_sequence", &self.receipt.published_sequence)
            .finish_non_exhaustive()
    }
}

impl<T> PublishedAllocation<T> {
    /// Shared observation of the published value.
    ///
    /// # Panics
    ///
    /// Never on a live guard: the value lives until close consumes it.
    #[must_use]
    pub fn observe(&self) -> &T {
        self.value
            .as_ref()
            .expect("published value lives until close")
    }

    /// Verified successful publication receipt.
    #[must_use]
    pub fn receipt(&self) -> &PublishedTransferReceipt {
        &self.receipt
    }

    /// Exact publication identity tuple.
    #[must_use]
    pub fn binding(&self) -> PublishedTransferBinding {
        self.receipt.binding
    }

    /// Charged authority bytes (one counted byte for zero-sized values).
    #[must_use]
    pub fn bytes(&self) -> u64 {
        self.receipt.bytes
    }

    /// Explicitly close exactly once: destroy the value first, then record
    /// the destination close against the exact published record.
    ///
    /// If the value destructor panics, unwinding destroys this guard with
    /// `destroyed` still false, so no close is recorded: the published record
    /// stays open and root seal refuses — fail-closed no-liveness, never a
    /// fabricated success.
    ///
    /// # Errors
    ///
    /// A refused close returns the guard with the value already destroyed
    /// (close semantics are destroy-then-record); its eventual drop retries
    /// the close implicitly.
    #[allow(clippy::result_large_err)]
    pub fn close(mut self) -> Result<PublishedTransferCloseReceipt, TypedCloseRejection<T>> {
        let value = self.value.take();
        if let Some(value) = value {
            drop(value);
        }
        self.destroyed = true;
        match close_published_transfer(&self.root, &self.receipt, false) {
            Ok(receipt) => {
                self.closed = true;
                Ok(receipt)
            }
            Err(refusal) => Err(TypedCloseRejection {
                published: self,
                refusal,
            }),
        }
    }
}

impl<T> Drop for PublishedAllocation<T> {
    fn drop(&mut self) {
        if let Some(value) = self.value.take() {
            drop(value);
            self.destroyed = true;
        }
        if !self.closed && self.destroyed {
            let _ = close_published_transfer(&self.root, &self.receipt, true);
            self.closed = true;
        }
        // A destructor panic mid-close leaves `destroyed == false`: recording
        // any close here would fabricate success after incomplete value
        // destruction. The published record stays open — fail-closed
        // no-liveness, exactly like deliberate forget.
    }
}

/// Preparation refusal carrying the allocation back unchanged.
impl<T> fmt::Debug for TypedPrepareRejection<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TypedPrepareRejection")
            .field("refusal", &self.refusal)
            .finish_non_exhaustive()
    }
}

/// Prepared-publication refusal: the structured ledger refusal plus the
/// prepared publication returned unchanged (record still prepared, value
/// intact) for a corrected retry.
pub struct TypedPrepareRejection<'owner, T> {
    allocation: LeasedAllocation<'owner, T>,
    refusal: PublishedTransferRefusal,
}

impl<'owner, T> TypedPrepareRejection<'owner, T> {
    /// Structured ledger refusal.
    #[must_use]
    pub fn refusal(&self) -> &PublishedTransferRefusal {
        &self.refusal
    }

    /// Recover the unchanged allocation for a corrected retry.
    #[must_use]
    pub fn into_allocation(self) -> LeasedAllocation<'owner, T> {
        self.allocation
    }
}

/// Commit refusal carrying the prepared publication back unchanged.
impl<T> fmt::Debug for TypedCommitRejection<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TypedCommitRejection")
            .field("refusal", &self.refusal)
            .finish_non_exhaustive()
    }
}

/// Commit refusal: the structured ledger refusal plus the prepared
/// publication returned unchanged for a corrected commit attempt.
pub struct TypedCommitRejection<'owner, T> {
    prepared: TypedPreparedPublication<'owner, T>,
    refusal: PublishedTransferRefusal,
}

impl<'owner, T> TypedCommitRejection<'owner, T> {
    /// Structured ledger refusal.
    #[must_use]
    pub fn refusal(&self) -> &PublishedTransferRefusal {
        &self.refusal
    }

    /// Recover the unchanged prepared publication.
    pub fn into_prepared(self) -> TypedPreparedPublication<'owner, T> {
        self.prepared
    }
}

/// Rollback refusal, split by the phase that refused.
pub enum TypedRollbackRejection<'owner, T> {
    /// The ledger refused the rollback itself; the prepared publication is
    /// returned unchanged — record still prepared, value intact.
    Prepared(
        TypedPreparedPublication<'owner, T>,
        PublishedTransferRefusal,
    ),
    /// The rollback committed but re-reservation was refused (for example a
    /// concurrent seal cut); the untouched value and the verified receipt
    /// are returned, with the staging authority back in the child pool.
    Released(T, PublishedTransferRollbackReceipt, DelegatedLeaseRefusal),
}

impl<T> fmt::Debug for TypedRollbackRejection<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Prepared(_, refusal) => f
                .debug_struct("TypedRollbackRejection::Prepared")
                .field("refusal", refusal)
                .finish_non_exhaustive(),
            Self::Released(_, _, refusal) => f
                .debug_struct("TypedRollbackRejection::Released")
                .field("refusal", refusal)
                .finish_non_exhaustive(),
        }
    }
}

/// Successful non-mutating rollback: staging plus untouched value.
impl<T> fmt::Debug for TypedRollback<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TypedRollback")
            .field("receipt", &self.receipt)
            .finish_non_exhaustive()
    }
}

/// Successful non-mutating rollback: the restaged allocation plus the
/// verified rollback receipt.
pub struct TypedRollback<'owner, T> {
    allocation: LeasedAllocation<'owner, T>,
    receipt: PublishedTransferRollbackReceipt,
}

impl<'owner, T> TypedRollback<'owner, T> {
    /// Verified rollback receipt.
    #[must_use]
    pub fn receipt(&self) -> &PublishedTransferRollbackReceipt {
        &self.receipt
    }

    /// Restaged allocation with the untouched value, ready to re-prepare
    /// under a fresh attempt generation or to drop normally.
    #[must_use]
    pub fn into_allocation(self) -> LeasedAllocation<'owner, T> {
        self.allocation
    }
}

/// Close refusal carrying the guard whose value was already destroyed.
impl<T> fmt::Debug for TypedCloseRejection<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TypedCloseRejection")
            .field("refusal", &self.refusal)
            .finish_non_exhaustive()
    }
}

/// Close refusal carrying the published guard whose value was already
/// destroyed; the eventual drop retries the close implicitly.
pub struct TypedCloseRejection<T> {
    published: PublishedAllocation<T>,
    refusal: PublishedTransferRefusal,
}

impl<T> TypedCloseRejection<T> {
    /// Structured ledger refusal.
    #[must_use]
    pub fn refusal(&self) -> &PublishedTransferRefusal {
        &self.refusal
    }

    /// Recover the guard; its value is destroyed and its eventual drop
    /// retries the close implicitly.
    pub fn into_published(self) -> PublishedAllocation<T> {
        self.published
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reserve_release_and_receipt_account_exactly() {
        let lease = OperationMemoryLease::bounded(1000);
        let a = lease.reserve("root", 600).expect("fits");
        assert_eq!(lease.receipt().used_bytes, 600);
        let refusal = lease.reserve("chunk", 500).expect_err("over limit");
        assert_eq!(refusal.what, "chunk");
        assert_eq!(refusal.used_bytes, 600);
        assert_eq!(refusal.limit_bytes, 1000);
        let b = lease.reserve("chunk", 400).expect("exactly fits");
        drop(a);
        drop(b);
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0);
        assert_eq!(receipt.requested_bytes, 1000);
        assert_eq!(receipt.peak_bytes, 1000);
        assert_eq!(receipt.refusals, 1);
        assert_eq!(receipt.release_invariant_violations, 0);
        let first = receipt.first_refusal.as_ref().expect("recorded");
        assert_eq!(first.what, "chunk");
        assert_eq!(first.requested_bytes, 500);
        assert!(receipt.to_json().contains("\"refusals\":1"));
    }

    #[test]
    fn unbounded_lease_accounts_within_the_representable_live_set() {
        let lease = OperationMemoryLease::unbounded();
        let charge = lease
            .reserve("huge", u64::MAX / 2)
            .expect("unbounded admits");
        assert_eq!(lease.receipt().peak_bytes, u64::MAX / 2);
        drop(charge);
        assert_eq!(lease.receipt().used_bytes, 0);
        assert_eq!(lease.receipt().refusals, 0);
    }

    #[test]
    fn charges_release_on_unwind() {
        let lease = OperationMemoryLease::bounded(100);
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _held = lease.reserve("tile", 80).expect("fits");
            panic!("tile body panicked");
        }));
        assert!(result.is_err());
        assert_eq!(
            lease.receipt().used_bytes,
            0,
            "an unwinding holder must release its charge"
        );
    }

    #[test]
    fn concurrent_reservations_never_exceed_the_limit() {
        let lease = OperationMemoryLease::bounded(64);
        std::thread::scope(|s| {
            for _ in 0..8 {
                let lease = lease.clone();
                s.spawn(move || {
                    for _ in 0..200 {
                        if let Ok(charge) = lease.reserve("hammer", 16) {
                            assert!(lease.receipt().used_bytes <= 64);
                            drop(charge);
                        }
                    }
                });
            }
        });
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0);
        assert!(receipt.peak_bytes <= 64);
    }

    #[test]
    fn cumulative_counters_are_exact_beyond_u64_and_refusal_json_is_escaped() {
        let lease = OperationMemoryLease::unbounded();
        drop(
            lease
                .reserve("first", u64::MAX - 1)
                .expect("representable live set"),
        );
        drop(
            lease
                .reserve("second", 2)
                .expect("sequential reservation fits"),
        );
        assert_eq!(
            lease.receipt().requested_bytes,
            u128::from(u64::MAX) + 1,
            "sequential reuse must not silently saturate cumulative demand"
        );

        let refusing = OperationMemoryLease::bounded(0);
        let hostile = "chunk\"}\n{\"forged\":true";
        refusing
            .reserve(hostile, 1)
            .expect_err("zero-byte limit refuses");
        let json = refusing.receipt().to_json();
        assert!(!json.contains('\n'));
        assert!(json.contains("chunk\\\"}\\n{\\\"forged\\\":true"));
    }

    #[test]
    fn unmatched_release_is_visible_and_fail_closed() {
        let lease = OperationMemoryLease::bounded(16);
        assert!(!lease.release_raw(1));
        let receipt = lease.receipt();
        assert_eq!(receipt.used_bytes, 0, "underflow must not wrap");
        assert_eq!(receipt.release_invariant_violations, 1);
        assert!(
            receipt
                .to_json()
                .contains("\"release_invariant_violations\":1")
        );
    }

    #[test]
    fn publication_ledger_commits_binding_and_outcome_but_not_later_destination_close() {
        let root_identity = LeaseIdentity::root(*b"unitpub1", [1; 32]);
        let child_identity = root_identity.child([2; 32], 0).expect("depth fits");
        let binding = PublishedTransferBinding::new([3; 32], [4; 32], [5; 32], [6; 32]);
        let record = PublicationRecord {
            parent_identity: None,
            child_identity,
            binding,
            envelope: PublishedTransferEnvelope::payload_only(8),
            bytes: 8,
            prepared_sequence: 4,
            resolved_sequence: Some(5),
            disposition: PublicationDisposition::Published,
            implicit_rollback: false,
            published_receipt_root: Some([7; 32]),
            destination_closed_sequence: None,
            implicit_destination_close: false,
            destination_close_root: None,
        };
        let original = publication_root(core::slice::from_ref(&record));

        let mut later_closed = record.clone();
        later_closed.destination_closed_sequence = Some(9);
        later_closed.implicit_destination_close = true;
        later_closed.destination_close_root = Some([8; 32]);
        assert_eq!(
            original,
            publication_root(core::slice::from_ref(&later_closed)),
            "destination lifetime is separate from the frozen publication outcome"
        );

        let mut destination_substitution = record.clone();
        destination_substitution.binding =
            PublishedTransferBinding::new([3; 32], [4; 32], [5; 32], [9; 32]);
        assert_ne!(
            original,
            publication_root(core::slice::from_ref(&destination_substitution)),
            "the destination identity is committed"
        );

        let mut rollback_substitution = record;
        rollback_substitution.disposition = PublicationDisposition::RolledBack;
        rollback_substitution.published_receipt_root = None;
        assert_ne!(
            original,
            publication_root(core::slice::from_ref(&rollback_substitution)),
            "publish and rollback are distinct terminal outcomes"
        );
    }
}
