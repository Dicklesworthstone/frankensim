//! Persistent entity identity: `Assembly -> Part -> (Region | Surface |
//! Interface)` as typed identities instead of strings.
//!
//! A scenario that binds a boundary condition to the STRING `"inlet"` is one
//! CAD rename away from silently meaning something else — or nothing at all.
//! This module makes the binding target a value: an [`EntityId`] derived from
//! the entity's kind, its parent path, its declared name, and (when the
//! importer supplies one) a geometry fingerprint. Display names live beside
//! identity and never enter it, so renaming is a receipt, not an orphaning.
//!
//! Every event that can move a binding — declaration, rename, import, re-mesh,
//! revision migration, retirement, legacy migration — appends an
//! [`IdentityReceipt`] to a hash-chained log. Resolution follows supersession
//! links and reports the WEAKEST evidence tier along the chain, so a strong
//! link never launders a weak one.
//!
//! What a content-derived id proves is narrow and stated in
//! [`GeometryFingerprint`]: equal ids mean the derivation inputs were
//! byte-equal. They are not a claim that two revisions describe the same
//! physical part.

use crate::frame::{FrameId, WORLD};
use crate::scenario::{Scenario, Violation};
use fs_blake3::{ContentHash, DomainHasher};
use fs_qty::{Dims, QtyAny};
use std::fmt;

/// Domain separation for entity identity derivation (wire-stable).
const ENTITY_ID_DOMAIN: &str = "fs-scenario/entity-id/v1";
/// Domain separation for geometry fingerprints (wire-stable).
const GEOMETRY_FINGERPRINT_DOMAIN: &str = "fs-scenario/geometry-fingerprint/v1";
/// Domain separation for identity receipts (wire-stable).
const IDENTITY_RECEIPT_DOMAIN: &str = "fs-scenario/identity-receipt/v1";
/// Domain separation for datum identity derivation (wire-stable).
const DATUM_ID_DOMAIN: &str = "fs-scenario/datum-id/v1";
/// Genesis payload of the receipt hash chain.
const RECEIPT_CHAIN_GENESIS: &[u8] = b"fs-scenario/identity-receipt-chain/genesis";

/// Maximum source bytes copied from one entity name into a diagnostic.
const DIAGNOSTIC_NAME_PREVIEW_BYTES: usize = 128;

/// Length dimensions (m): the only tolerance magnitude admitted by V1.
const LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Maximum datums in one tolerance's ordered datum reference frame
/// (primary, secondary, tertiary — ASME Y14.5 practice).
pub const MAX_DATUM_FRAME_LEN: usize = 3;

/// Bounded name rendering for diagnostics: long names are truncated at a UTF-8
/// boundary and the exact original byte length is appended.
#[derive(Clone, Copy)]
struct NamePreview<'a>(&'a str);

impl fmt::Debug for NamePreview<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.len() <= DIAGNOSTIC_NAME_PREVIEW_BYTES {
            return fmt::Debug::fmt(self.0, formatter);
        }
        let mut end = DIAGNOSTIC_NAME_PREVIEW_BYTES;
        while !self.0.is_char_boundary(end) {
            end -= 1;
        }
        fmt::Debug::fmt(&self.0[..end], formatter)?;
        write!(formatter, "…<{} bytes total>", self.0.len())
    }
}

fn preview(name: &str) -> String {
    format!("{:?}", NamePreview(name))
}

fn absorb_bytes(hasher: &mut DomainHasher, bytes: &[u8]) {
    hasher.update(&(bytes.len() as u64).to_le_bytes());
    hasher.update(bytes);
}

fn absorb_entity_id(hasher: &mut DomainHasher, id: Option<EntityId>) {
    match id {
        None => hasher.update(&[0u8]),
        Some(id) => {
            hasher.update(&[1u8, id.kind.tag()]);
            hasher.update(id.digest.as_bytes());
        }
    }
}

fn absorb_fingerprint(hasher: &mut DomainHasher, fingerprint: Option<GeometryFingerprint>) {
    match fingerprint {
        None => hasher.update(&[0u8]),
        Some(fingerprint) => {
            hasher.update(&[1u8]);
            hasher.update(fingerprint.0.as_bytes());
        }
    }
}

fn push_checked<T>(
    values: &mut Vec<T>,
    value: T,
    resource: &'static str,
) -> Result<(), EntityError> {
    values
        .try_reserve(1)
        .map_err(|_| EntityError::AllocationRefused { resource })?;
    values.push(value);
    Ok(())
}

// ---------------------------------------------------------------------------
// Kinds and identity
// ---------------------------------------------------------------------------

/// The closed set of entity kinds.
///
/// The containment rule is `Assembly -> (Assembly | Part)`,
/// `Part -> (Region | Surface)`, and `Interface` under either the assembly or
/// the part that owns both of its sides.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntityKind {
    /// A container of parts and sub-assemblies.
    Assembly,
    /// One manufactured or modelled body.
    Part,
    /// A volumetric subdomain of a part.
    Region,
    /// A boundary patch of a part.
    Surface,
    /// A first-class pairing of two sides.
    Interface,
}

impl EntityKind {
    /// Wire-stable tag absorbed into identity derivation.
    const fn tag(self) -> u8 {
        match self {
            EntityKind::Assembly => 1,
            EntityKind::Part => 2,
            EntityKind::Region => 3,
            EntityKind::Surface => 4,
            EntityKind::Interface => 5,
        }
    }

    /// Stable lowercase label used in diagnostics and identity tokens.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            EntityKind::Assembly => "assembly",
            EntityKind::Part => "part",
            EntityKind::Region => "region",
            EntityKind::Surface => "surface",
            EntityKind::Interface => "interface",
        }
    }

    /// Whether this kind admits `parent` as its containing kind.
    #[must_use]
    pub const fn admits_parent(self, parent: Option<EntityKind>) -> bool {
        matches!(
            (self, parent),
            (EntityKind::Assembly, None | Some(EntityKind::Assembly))
                | (EntityKind::Part, Some(EntityKind::Assembly))
                | (
                    EntityKind::Region | EntityKind::Surface,
                    Some(EntityKind::Part)
                )
                | (
                    EntityKind::Interface,
                    Some(EntityKind::Assembly | EntityKind::Part)
                )
        )
    }

    /// Whether this kind may stand on one side of an interface.
    #[must_use]
    pub const fn admits_interface_side(self) -> bool {
        matches!(self, EntityKind::Surface | EntityKind::Region)
    }
}

impl fmt::Display for EntityKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// A 32-byte content root over geometry bytes the IMPORTER chose to hash.
///
/// # What this proves
/// Two equal fingerprints prove the supplied byte strings were equal under
/// BLAKE3. That is all.
///
/// # What this does not prove
/// It is not a claim that the two revisions describe the same physical part,
/// that either byte string is a valid geometry, that a re-export of an
/// unchanged part reproduces the same bytes (it usually does not), or that
/// different fingerprints mean different geometry. This crate never parses,
/// canonicalizes, or meshes geometry; it hashes exactly what it is given.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeometryFingerprint(ContentHash);

impl GeometryFingerprint {
    /// Fingerprint the exact supplied bytes under the fingerprint domain.
    #[must_use]
    pub fn of_bytes(bytes: &[u8]) -> Self {
        let mut hasher = DomainHasher::new(GEOMETRY_FINGERPRINT_DOMAIN);
        absorb_bytes(&mut hasher, bytes);
        GeometryFingerprint(hasher.finalize())
    }

    /// Adopt a root an upstream layer already computed.
    #[must_use]
    pub const fn from_hash(hash: ContentHash) -> Self {
        GeometryFingerprint(hash)
    }

    /// The underlying 32-byte root.
    #[must_use]
    pub const fn hash(self) -> ContentHash {
        self.0
    }
}

impl fmt::Display for GeometryFingerprint {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, formatter)
    }
}

/// A stable, content-derived entity identity.
///
/// The digest covers, in a length-prefixed and therefore unambiguous preimage:
/// the kind tag, the parent identity (which transitively encodes the whole
/// parent path), the declared name, the geometry fingerprint when present, and
/// the interface pairing when present. Display names are NOT part of the
/// preimage.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityId {
    kind: EntityKind,
    digest: ContentHash,
}

impl EntityId {
    /// Reconstruct an already-derived identity carried by a checked crate wire
    /// envelope. This does not derive or validate a declaration preimage.
    pub(crate) const fn from_wire(kind: EntityKind, digest: ContentHash) -> Self {
        Self { kind, digest }
    }

    /// The kind this identity was derived for.
    #[must_use]
    pub const fn kind(self) -> EntityKind {
        self.kind
    }

    /// The 32-byte identity root.
    #[must_use]
    pub const fn digest(self) -> ContentHash {
        self.digest
    }

    /// Stable `kind:hex` token used in diagnostics and reports.
    #[must_use]
    pub fn token(self) -> String {
        format!("{}:{}", self.kind.label(), self.digest.to_hex())
    }

    /// Short `kind:hex16` token for human-facing tables.
    #[must_use]
    pub fn short_token(self) -> String {
        let hex = self.digest.to_hex();
        format!("{}:{}", self.kind.label(), &hex[..16])
    }
}

impl fmt::Debug for EntityId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.short_token())
    }
}

impl fmt::Display for EntityId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.token())
    }
}

/// How the two sides of an interface relate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InterfacePairing {
    /// Physics distinguishes the sides: `from` is the side a one-sided
    /// treatment (a TIM layer, a coating, a film coefficient) is applied to.
    Ordered,
    /// Physics does not distinguish the sides; the pair is canonicalized so
    /// `(a, b)` and `(b, a)` are ONE identity.
    Unordered,
}

impl InterfacePairing {
    const fn tag(self) -> u8 {
        match self {
            InterfacePairing::Ordered => 1,
            InterfacePairing::Unordered => 2,
        }
    }
}

/// A first-class interface pair with its own ordering semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InterfacePair {
    from: EntityId,
    to: EntityId,
    pairing: InterfacePairing,
}

impl InterfacePair {
    /// An ordered pair: `from` is the side the one-sided treatment applies to.
    /// `(a, b)` and `(b, a)` are DIFFERENT identities.
    #[must_use]
    pub const fn ordered(from: EntityId, to: EntityId) -> Self {
        InterfacePair {
            from,
            to,
            pairing: InterfacePairing::Ordered,
        }
    }

    /// An unordered pair, canonicalized so declaration order cannot create two
    /// identities for one physical interface.
    #[must_use]
    pub fn unordered(a: EntityId, b: EntityId) -> Self {
        let (from, to) = if a <= b { (a, b) } else { (b, a) };
        InterfacePair {
            from,
            to,
            pairing: InterfacePairing::Unordered,
        }
    }

    /// First side (the applied-to side when ordered).
    #[must_use]
    pub const fn from(self) -> EntityId {
        self.from
    }

    /// Second side.
    #[must_use]
    pub const fn to(self) -> EntityId {
        self.to
    }

    /// Declared pairing semantics.
    #[must_use]
    pub const fn pairing(self) -> InterfacePairing {
        self.pairing
    }

    /// The side a one-sided treatment is applied to, or `None` for an
    /// unordered pair — a refusal to answer, never a silent `from`.
    #[must_use]
    pub const fn applied_side(self) -> Option<EntityId> {
        match self.pairing {
            InterfacePairing::Ordered => Some(self.from),
            InterfacePairing::Unordered => None,
        }
    }
}

/// The declaration an identity is derived from.
///
/// Construct with [`EntityDeclaration::assembly`] and friends, then attach the
/// optional geometry fingerprint and display name. [`EntityDeclaration::identity`]
/// is pure: it needs no catalog.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntityDeclaration {
    kind: EntityKind,
    parent: Option<EntityId>,
    declared_name: String,
    display_name: String,
    fingerprint: Option<GeometryFingerprint>,
    pair: Option<InterfacePair>,
    legacy: bool,
}

impl EntityDeclaration {
    fn new(kind: EntityKind, parent: Option<EntityId>, declared_name: &str) -> Self {
        EntityDeclaration {
            kind,
            parent,
            declared_name: declared_name.to_string(),
            display_name: declared_name.to_string(),
            fingerprint: None,
            pair: None,
            legacy: false,
        }
    }

    /// A root assembly.
    #[must_use]
    pub fn assembly(declared_name: &str) -> Self {
        Self::new(EntityKind::Assembly, None, declared_name)
    }

    /// A sub-assembly of `parent`.
    #[must_use]
    pub fn sub_assembly(parent: EntityId, declared_name: &str) -> Self {
        Self::new(EntityKind::Assembly, Some(parent), declared_name)
    }

    /// A part inside an assembly.
    #[must_use]
    pub fn part(parent: EntityId, declared_name: &str) -> Self {
        Self::new(EntityKind::Part, Some(parent), declared_name)
    }

    /// A volumetric region of a part.
    #[must_use]
    pub fn region(parent: EntityId, declared_name: &str) -> Self {
        Self::new(EntityKind::Region, Some(parent), declared_name)
    }

    /// A boundary surface of a part.
    #[must_use]
    pub fn surface(parent: EntityId, declared_name: &str) -> Self {
        Self::new(EntityKind::Surface, Some(parent), declared_name)
    }

    /// An interface owned by `parent` between the two sides of `pair`.
    #[must_use]
    pub fn interface(parent: EntityId, declared_name: &str, pair: InterfacePair) -> Self {
        let mut declaration = Self::new(EntityKind::Interface, Some(parent), declared_name);
        declaration.pair = Some(pair);
        declaration
    }

    /// Attach a geometry fingerprint (identity-bearing).
    #[must_use]
    pub fn with_fingerprint(mut self, fingerprint: GeometryFingerprint) -> Self {
        self.fingerprint = Some(fingerprint);
        self
    }

    /// Set the display name (NOT identity-bearing).
    #[must_use]
    pub fn with_display_name(mut self, display_name: &str) -> Self {
        self.display_name = display_name.to_string();
        self
    }

    /// Mark this declaration as mechanically migrated from a bare string.
    /// The marker is metadata, not identity, so migration is idempotent.
    #[must_use]
    pub fn with_legacy_marker(mut self) -> Self {
        self.legacy = true;
        self
    }

    /// Declared kind.
    #[must_use]
    pub const fn kind(&self) -> EntityKind {
        self.kind
    }

    /// Parent identity, or `None` for a root assembly.
    #[must_use]
    pub const fn parent(&self) -> Option<EntityId> {
        self.parent
    }

    /// The identity-bearing declared name.
    #[must_use]
    pub fn declared_name(&self) -> &str {
        &self.declared_name
    }

    /// The display name at declaration time.
    #[must_use]
    pub fn display_name(&self) -> &str {
        &self.display_name
    }

    /// Geometry fingerprint, when the importer supplied one.
    #[must_use]
    pub const fn fingerprint(&self) -> Option<GeometryFingerprint> {
        self.fingerprint
    }

    /// Interface pairing, for `Interface` declarations.
    #[must_use]
    pub const fn pair(&self) -> Option<InterfacePair> {
        self.pair
    }

    /// Whether this declaration carries the legacy-string marker.
    #[must_use]
    pub const fn is_legacy(&self) -> bool {
        self.legacy
    }

    /// Derive the identity of this declaration. Pure: no catalog is consulted,
    /// so callers can compute the identity of an entity before declaring it.
    #[must_use]
    pub fn identity(&self) -> EntityId {
        let mut hasher = DomainHasher::new(ENTITY_ID_DOMAIN);
        hasher.update(&[self.kind.tag()]);
        absorb_entity_id(&mut hasher, self.parent);
        absorb_bytes(&mut hasher, self.declared_name.as_bytes());
        absorb_fingerprint(&mut hasher, self.fingerprint);
        match self.pair {
            None => hasher.update(&[0u8]),
            Some(pair) => {
                hasher.update(&[1u8, pair.pairing.tag()]);
                absorb_entity_id(&mut hasher, Some(pair.from));
                absorb_entity_id(&mut hasher, Some(pair.to));
            }
        }
        EntityId {
            kind: self.kind,
            digest: hasher.finalize(),
        }
    }

    /// Whether two declarations have the same identity PREIMAGE — the test
    /// that separates an honest duplicate declaration from a digest collision.
    #[must_use]
    pub fn same_identity_preimage(&self, other: &EntityDeclaration) -> bool {
        self.kind == other.kind
            && self.parent == other.parent
            && self.declared_name == other.declared_name
            && self.fingerprint == other.fingerprint
            && self.pair == other.pair
    }
}

/// Lifecycle status of one catalog entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityStatus {
    /// Current.
    Active,
    /// Replaced by a later identity; references resolve forward.
    Superseded {
        /// The identity that replaced this one.
        successor: EntityId,
        /// Sequence of the receipt that recorded the supersession.
        receipt: u64,
    },
    /// Removed by a complete revision; references do NOT resolve forward.
    Retired {
        /// Sequence of the receipt that recorded the retirement.
        receipt: u64,
    },
}

/// One catalog entry: an immutable identity plus mutable display state.
#[derive(Debug, Clone, PartialEq)]
pub struct Entity {
    id: EntityId,
    declaration: EntityDeclaration,
    display_name: String,
    status: EntityStatus,
    row: usize,
}

impl Entity {
    /// Immutable identity.
    #[must_use]
    pub const fn id(&self) -> EntityId {
        self.id
    }

    /// The declaration this identity was derived from.
    #[must_use]
    pub const fn declaration(&self) -> &EntityDeclaration {
        &self.declaration
    }

    /// Entity kind.
    #[must_use]
    pub const fn kind(&self) -> EntityKind {
        self.id.kind
    }

    /// Parent identity, or `None` for a root assembly.
    #[must_use]
    pub const fn parent(&self) -> Option<EntityId> {
        self.declaration.parent
    }

    /// The identity-bearing declared name.
    #[must_use]
    pub fn declared_name(&self) -> &str {
        &self.declaration.declared_name
    }

    /// The CURRENT display name (renames move this, never the identity).
    #[must_use]
    pub fn display_name(&self) -> &str {
        &self.display_name
    }

    /// Geometry fingerprint, when one was supplied.
    #[must_use]
    pub const fn fingerprint(&self) -> Option<GeometryFingerprint> {
        self.declaration.fingerprint
    }

    /// Interface pairing, for interfaces.
    #[must_use]
    pub const fn pair(&self) -> Option<InterfacePair> {
        self.declaration.pair
    }

    /// Lifecycle status.
    #[must_use]
    pub const fn status(&self) -> EntityStatus {
        self.status
    }

    /// Whether this entity was mechanically migrated from a bare string.
    #[must_use]
    pub const fn is_legacy(&self) -> bool {
        self.declaration.legacy
    }

    /// Declaration row (deterministic declaration order).
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

// ---------------------------------------------------------------------------
// Receipts
// ---------------------------------------------------------------------------

/// The event that produced an identity receipt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RebindEvent {
    /// An entity was authored directly.
    Declaration,
    /// A display name changed; identity did not.
    Rename,
    /// A geometry import revision was applied.
    Import,
    /// A re-mesh revision was applied.
    Remesh,
    /// A revision (design-history) migration was applied.
    RevisionMigration,
    /// A bare string was mechanically adopted as a declared name.
    LegacyMigration,
    /// A complete revision did not contain a previously active entity.
    Retirement,
}

impl RebindEvent {
    const fn tag(self) -> u8 {
        match self {
            RebindEvent::Declaration => 1,
            RebindEvent::Rename => 2,
            RebindEvent::Import => 3,
            RebindEvent::Remesh => 4,
            RebindEvent::RevisionMigration => 5,
            RebindEvent::LegacyMigration => 6,
            RebindEvent::Retirement => 7,
        }
    }

    /// Stable label used in receipts and reports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            RebindEvent::Declaration => "declaration",
            RebindEvent::Rename => "rename",
            RebindEvent::Import => "import",
            RebindEvent::Remesh => "remesh",
            RebindEvent::RevisionMigration => "revision-migration",
            RebindEvent::LegacyMigration => "legacy-migration",
            RebindEvent::Retirement => "retirement",
        }
    }

    /// Whether this event may be carried by an [`ImportRevision`].
    #[must_use]
    pub const fn is_revision_event(self) -> bool {
        matches!(
            self,
            RebindEvent::Import | RebindEvent::Remesh | RebindEvent::RevisionMigration
        )
    }
}

impl fmt::Display for RebindEvent {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// How strong the evidence behind a correspondence is. Ordered weakest first.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceTier {
    /// A human or upstream tool asserted it; this crate verified nothing.
    Asserted,
    /// A declared name was adopted; no structural or content evidence.
    Declared,
    /// The declared name and parent path corresponded.
    PathMatched,
    /// The geometry fingerprints were byte-equal.
    ContentMatched,
    /// The identity itself did not change.
    Identical,
}

impl EvidenceTier {
    /// Stable label used in receipts and reports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            EvidenceTier::Asserted => "asserted",
            EvidenceTier::Declared => "declared",
            EvidenceTier::PathMatched => "path-matched",
            EvidenceTier::ContentMatched => "content-matched",
            EvidenceTier::Identical => "identical",
        }
    }
}

impl fmt::Display for EvidenceTier {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// Why a receipt claims one identity corresponds to another.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MatchBasis {
    /// A new identity with no claimed predecessor.
    NewIdentity,
    /// The identity did not change across the event.
    Unchanged,
    /// Geometry fingerprints were byte-equal; the declared path changed.
    GeometryFingerprint,
    /// Declared name and parent path corresponded; geometry bytes differ or
    /// are absent on one side.
    DeclaredPath,
    /// A bare legacy string was adopted as the declared name.
    LegacyName,
    /// The caller asserted the correspondence.
    Asserted,
    /// A complete revision did not contain the entity.
    Absent,
}

impl MatchBasis {
    const fn tag(self) -> u8 {
        match self {
            MatchBasis::NewIdentity => 1,
            MatchBasis::Unchanged => 2,
            MatchBasis::GeometryFingerprint => 3,
            MatchBasis::DeclaredPath => 4,
            MatchBasis::LegacyName => 5,
            MatchBasis::Asserted => 6,
            MatchBasis::Absent => 7,
        }
    }

    /// Stable label used in receipts and reports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            MatchBasis::NewIdentity => "new-identity",
            MatchBasis::Unchanged => "unchanged",
            MatchBasis::GeometryFingerprint => "geometry-fingerprint",
            MatchBasis::DeclaredPath => "declared-path",
            MatchBasis::LegacyName => "legacy-name",
            MatchBasis::Asserted => "asserted",
            MatchBasis::Absent => "absent",
        }
    }

    /// The evidence tier of a supersession carrying this basis, or `None` when
    /// the basis does not link two identities.
    #[must_use]
    pub const fn tier(self) -> Option<EvidenceTier> {
        match self {
            MatchBasis::NewIdentity | MatchBasis::Absent => None,
            MatchBasis::Unchanged => Some(EvidenceTier::Identical),
            MatchBasis::GeometryFingerprint => Some(EvidenceTier::ContentMatched),
            MatchBasis::DeclaredPath => Some(EvidenceTier::PathMatched),
            MatchBasis::LegacyName => Some(EvidenceTier::Declared),
            MatchBasis::Asserted => Some(EvidenceTier::Asserted),
        }
    }

    /// Whether this basis establishes that geometry bytes matched.
    ///
    /// True ONLY for [`MatchBasis::GeometryFingerprint`], and even then the
    /// claim is byte equality of what the importer hashed — see
    /// [`GeometryFingerprint`].
    #[must_use]
    pub const fn proves_geometry_bytes_matched(self) -> bool {
        matches!(self, MatchBasis::GeometryFingerprint)
    }
}

impl fmt::Display for MatchBasis {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// Immutable evidence for one identity event, chained to its predecessor.
///
/// The digest binds the previous chain root, the sequence, the event, the
/// subject and predecessor identities, the basis, and both display names.
/// [`IdentityReceipt::verifies`] recomputes it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityReceipt {
    sequence: u64,
    event: RebindEvent,
    subject: EntityId,
    predecessor: Option<EntityId>,
    basis: MatchBasis,
    display_before: Option<String>,
    display_after: Option<String>,
    previous: ContentHash,
    digest: ContentHash,
}

impl IdentityReceipt {
    // The digest binds every retained field; bundling them into a struct here
    // would only move the same arity behind a constructor.
    #[allow(clippy::too_many_arguments)]
    fn compute_digest(
        previous: ContentHash,
        sequence: u64,
        event: RebindEvent,
        subject: EntityId,
        predecessor: Option<EntityId>,
        basis: MatchBasis,
        display_before: Option<&str>,
        display_after: Option<&str>,
    ) -> ContentHash {
        let mut hasher = DomainHasher::new(IDENTITY_RECEIPT_DOMAIN);
        hasher.update(previous.as_bytes());
        hasher.update(&sequence.to_le_bytes());
        hasher.update(&[event.tag(), basis.tag()]);
        absorb_entity_id(&mut hasher, Some(subject));
        absorb_entity_id(&mut hasher, predecessor);
        match display_before {
            None => hasher.update(&[0u8]),
            Some(name) => {
                hasher.update(&[1u8]);
                absorb_bytes(&mut hasher, name.as_bytes());
            }
        }
        match display_after {
            None => hasher.update(&[0u8]),
            Some(name) => {
                hasher.update(&[1u8]);
                absorb_bytes(&mut hasher, name.as_bytes());
            }
        }
        hasher.finalize()
    }

    /// Position in the receipt log (also the chain index).
    #[must_use]
    pub const fn sequence(&self) -> u64 {
        self.sequence
    }

    /// The event that produced this receipt.
    #[must_use]
    pub const fn event(&self) -> RebindEvent {
        self.event
    }

    /// The identity the receipt is about (the successor, for supersessions).
    #[must_use]
    pub const fn subject(&self) -> EntityId {
        self.subject
    }

    /// The superseded identity, when this receipt links two identities.
    #[must_use]
    pub const fn predecessor(&self) -> Option<EntityId> {
        self.predecessor
    }

    /// Why the correspondence is claimed.
    #[must_use]
    pub const fn basis(&self) -> MatchBasis {
        self.basis
    }

    /// Display name before the event, when the event changed it.
    #[must_use]
    pub fn display_before(&self) -> Option<&str> {
        self.display_before.as_deref()
    }

    /// Display name after the event.
    #[must_use]
    pub fn display_after(&self) -> Option<&str> {
        self.display_after.as_deref()
    }

    /// Chain root this receipt extends.
    #[must_use]
    pub const fn previous(&self) -> ContentHash {
        self.previous
    }

    /// This receipt's chain root.
    #[must_use]
    pub const fn digest(&self) -> ContentHash {
        self.digest
    }

    /// Recompute the digest from the retained fields.
    ///
    /// This detects mutation of any recorded field. It does not prove the
    /// receipt describes a real geometry operation, and — on its own — it does
    /// not detect truncation of the log's tail; that needs the chain root
    /// pinned outside this catalog.
    #[must_use]
    pub fn verifies(&self) -> bool {
        Self::compute_digest(
            self.previous,
            self.sequence,
            self.event,
            self.subject,
            self.predecessor,
            self.basis,
            self.display_before.as_deref(),
            self.display_after.as_deref(),
        ) == self.digest
    }
}

// ---------------------------------------------------------------------------
// Datums, tolerances, placements
// ---------------------------------------------------------------------------

/// The reference-feature kind of a datum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DatumFeature {
    /// A datum plane.
    Plane,
    /// A datum axis.
    Axis,
    /// A datum point.
    Point,
}

impl DatumFeature {
    const fn tag(self) -> u8 {
        match self {
            DatumFeature::Plane => 1,
            DatumFeature::Axis => 2,
            DatumFeature::Point => 3,
        }
    }

    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            DatumFeature::Plane => "plane",
            DatumFeature::Axis => "axis",
            DatumFeature::Point => "point",
        }
    }
}

impl fmt::Display for DatumFeature {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// A content-derived datum identity.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DatumId(ContentHash);

impl DatumId {
    /// The 32-byte identity root.
    #[must_use]
    pub const fn digest(self) -> ContentHash {
        self.0
    }

    /// Stable `datum:hex16` token for reports.
    #[must_use]
    pub fn short_token(self) -> String {
        let hex = self.0.to_hex();
        format!("datum:{}", &hex[..16])
    }
}

impl fmt::Debug for DatumId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.short_token())
    }
}

impl fmt::Display for DatumId {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "datum:{}", self.0.to_hex())
    }
}

/// A named reference feature attached to an entity.
///
/// The identity derives from the owner, the declared name, the feature kind,
/// and the referenced datums — so the datum hierarchy is a DAG by
/// construction: a datum cannot name a datum that does not exist yet, and its
/// identity depends on the identities it names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Datum {
    id: DatumId,
    owner: EntityId,
    declared_name: String,
    feature: DatumFeature,
    references: Vec<DatumId>,
    row: usize,
}

impl Datum {
    /// Content-derived identity.
    #[must_use]
    pub const fn id(&self) -> DatumId {
        self.id
    }

    /// The entity this datum is attached to.
    #[must_use]
    pub const fn owner(&self) -> EntityId {
        self.owner
    }

    /// Identity-bearing declared name (drawing letter, typically).
    #[must_use]
    pub fn declared_name(&self) -> &str {
        &self.declared_name
    }

    /// Reference-feature kind.
    #[must_use]
    pub const fn feature(&self) -> DatumFeature {
        self.feature
    }

    /// Datums that establish this one, in declared order.
    #[must_use]
    pub fn references(&self) -> &[DatumId] {
        &self.references
    }

    /// Declaration row.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// The declared geometric characteristic a tolerance controls.
///
/// V1 admits only length-dimensioned characteristics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ToleranceKind {
    /// Position (requires a datum reference frame).
    Position,
    /// Flatness (a form control: no datums).
    Flatness,
    /// Parallelism (requires a datum reference frame).
    Parallelism,
    /// Perpendicularity (requires a datum reference frame).
    Perpendicularity,
    /// Surface profile (requires a datum reference frame).
    Profile,
    /// Declared assembly gap at an interface (requires a datum frame).
    InterfaceGap,
}

impl ToleranceKind {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            ToleranceKind::Position => "position",
            ToleranceKind::Flatness => "flatness",
            ToleranceKind::Parallelism => "parallelism",
            ToleranceKind::Perpendicularity => "perpendicularity",
            ToleranceKind::Profile => "profile",
            ToleranceKind::InterfaceGap => "interface-gap",
        }
    }

    /// Whether this characteristic requires a datum reference frame.
    ///
    /// Form controls (flatness) take no datums; orientation, location, and
    /// profile controls require at least one.
    #[must_use]
    pub const fn requires_datum_frame(self) -> bool {
        !matches!(self, ToleranceKind::Flatness)
    }

    /// The SI dimensions the magnitude must carry.
    #[must_use]
    pub const fn expected_dims(self) -> Dims {
        LENGTH_DIMS
    }
}

impl fmt::Display for ToleranceKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// Where a declared tolerance came from. There is no unsourced tolerance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ToleranceSource {
    /// A drawing callout.
    Drawing {
        /// Sheet identifier.
        sheet: String,
        /// Note or feature-control-frame identifier.
        note: String,
    },
    /// A standard clause.
    Standard {
        /// Clause identifier, e.g. `ASME Y14.5-2018 §7.3`.
        clause: String,
    },
    /// An engineering assumption, with the rationale that justifies it.
    Assumed {
        /// Why this value was assumed.
        rationale: String,
    },
}

impl ToleranceSource {
    /// The free-text field that must be nonempty for this source.
    fn nonempty_fields(&self) -> [(&'static str, &str); 2] {
        match self {
            ToleranceSource::Drawing { sheet, note } => [("sheet", sheet.as_str()), ("note", note)],
            ToleranceSource::Standard { clause } => {
                [("clause", clause.as_str()), ("clause", clause)]
            }
            ToleranceSource::Assumed { rationale } => {
                [("rationale", rationale.as_str()), ("rationale", rationale)]
            }
        }
    }

    /// Stable label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        match self {
            ToleranceSource::Drawing { .. } => "drawing",
            ToleranceSource::Standard { .. } => "standard",
            ToleranceSource::Assumed { .. } => "assumed",
        }
    }
}

/// A declared geometric tolerance on an entity.
#[derive(Debug, Clone, PartialEq)]
pub struct Tolerance {
    subject: EntityId,
    kind: ToleranceKind,
    magnitude: QtyAny,
    datum_frame: Vec<DatumId>,
    source: ToleranceSource,
    row: usize,
}

impl Tolerance {
    /// Declare a tolerance. Structural admission happens in
    /// [`EntityCatalog::declare_tolerance`]; dimensional and datum-frame
    /// checks are reported by [`EntityCatalog::validate`].
    #[must_use]
    pub fn new(
        subject: EntityId,
        kind: ToleranceKind,
        magnitude: QtyAny,
        datum_frame: Vec<DatumId>,
        source: ToleranceSource,
    ) -> Self {
        Tolerance {
            subject,
            kind,
            magnitude,
            datum_frame,
            source,
            row: 0,
        }
    }

    /// The entity or interface this tolerance controls.
    #[must_use]
    pub const fn subject(&self) -> EntityId {
        self.subject
    }

    /// The controlled characteristic.
    #[must_use]
    pub const fn kind(&self) -> ToleranceKind {
        self.kind
    }

    /// The dimensioned magnitude.
    #[must_use]
    pub const fn magnitude(&self) -> QtyAny {
        self.magnitude
    }

    /// The ordered datum reference frame (primary, secondary, tertiary).
    #[must_use]
    pub fn datum_frame(&self) -> &[DatumId] {
        &self.datum_frame
    }

    /// Declared source.
    #[must_use]
    pub const fn source(&self) -> &ToleranceSource {
        &self.source
    }

    /// Declaration row.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

/// What a placement record claims about where an occurrence sits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PlacementBasis {
    /// The nominal (as-designed) placement.
    Nominal,
    /// The measured (as-built) placement, cited by the content identity of a
    /// calibrated registration record.
    ///
    /// The citation is a bare `ContentHash` on purpose. This crate is L3 and
    /// the registration record lives in L2 `fs-asbuilt`; carrying the record
    /// itself would either invert the layer direction or duplicate a 6-dof
    /// covariance that already has exactly one owner. What this variant
    /// asserts is therefore narrow and checkable: *this occurrence is placed
    /// on the authority of the artifact with this identity*.
    ///
    /// It asserts nothing about that artifact's existence, authenticity, or
    /// fitness. Resolving the identity back to a
    /// `fs_asbuilt::rigid3::CalibratedRigid3Registration` — and refusing when
    /// it does not resolve — is the product layer's obligation. See the
    /// no-claim boundary in `CONTRACT.md`.
    AsBuilt {
        /// Content identity of the calibrated registration record, as
        /// published by `CalibratedRigid3Registration::model_identity()`.
        registration_ref: ContentHash,
    },
}

impl PlacementBasis {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            PlacementBasis::Nominal => "nominal",
            PlacementBasis::AsBuilt { .. } => "as-built",
        }
    }

    /// The cited registration identity, when this basis is measured.
    ///
    /// `None` for [`PlacementBasis::Nominal`]: a nominal placement cites no
    /// measurement, which is a different statement from citing one that
    /// happens to be unavailable.
    #[must_use]
    pub const fn registration_ref(self) -> Option<ContentHash> {
        match self {
            PlacementBasis::Nominal => None,
            PlacementBasis::AsBuilt { registration_ref } => Some(registration_ref),
        }
    }
}

/// A typed placement record binding one part occurrence to a scenario frame.
///
/// The transform itself is the scenario's existing [`crate::frame::FrameTree`]
/// entry: this record does not introduce a parallel frame system, it names
/// which frame carries the occurrence's placement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    occurrence: EntityId,
    frame: FrameId,
    basis: PlacementBasis,
    row: usize,
}

impl Placement {
    /// The placed part occurrence.
    #[must_use]
    pub const fn occurrence(&self) -> EntityId {
        self.occurrence
    }

    /// The frame carrying the placement transform.
    #[must_use]
    pub const fn frame(&self) -> FrameId {
        self.frame
    }

    /// What the placement claims.
    #[must_use]
    pub const fn basis(&self) -> PlacementBasis {
        self.basis
    }

    /// Declaration row.
    #[must_use]
    pub const fn row(&self) -> usize {
        self.row
    }
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Explicit admission limits for one catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EntityBudget {
    /// Maximum entities in the catalog.
    pub max_entities: usize,
    /// Maximum receipts in the log.
    pub max_receipts: usize,
    /// Maximum bytes in any declared or display name.
    pub max_name_bytes: usize,
    /// Maximum supersession hops followed during resolution.
    pub max_supersession_hops: usize,
    /// Maximum parent hops followed during containment checks.
    pub max_hierarchy_depth: usize,
    /// Maximum datums in the catalog.
    pub max_datums: usize,
    /// Maximum tolerances in the catalog.
    pub max_tolerances: usize,
    /// Maximum placements in the catalog.
    pub max_placements: usize,
    /// Maximum entries in one binding table.
    pub max_bindings: usize,
}

/// The default entity budget: explicit, and named at the call site.
pub const DEFAULT_ENTITY_BUDGET: EntityBudget = EntityBudget {
    max_entities: 65_536,
    max_receipts: 262_144,
    max_name_bytes: 4_096,
    max_supersession_hops: 64,
    max_hierarchy_depth: 64,
    max_datums: 16_384,
    max_tolerances: 16_384,
    max_placements: 65_536,
    max_bindings: 262_144,
};

impl Default for EntityBudget {
    fn default() -> Self {
        DEFAULT_ENTITY_BUDGET
    }
}

/// Structured refusals from the entity catalog.
///
/// Every variant carries a fix hint through [`EntityError::fix`] and converts
/// into the crate's `Violation` shape through [`EntityError::into_violation`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EntityError {
    /// A declared or display name was empty.
    EmptyName {
        /// The kind being declared.
        kind: EntityKind,
        /// Which name field was empty.
        field: &'static str,
    },
    /// A name exceeded the budget's byte cap.
    NameTooLong {
        /// The kind being declared.
        kind: EntityKind,
        /// Which name field was too long.
        field: &'static str,
        /// Supplied byte length.
        bytes: usize,
        /// Admitted byte length.
        limit: usize,
    },
    /// The declared parent is not in the catalog.
    UnknownParent {
        /// The kind being declared.
        kind: EntityKind,
        /// The unresolved parent identity.
        parent: EntityId,
    },
    /// The declared parent is superseded or retired.
    InactiveParent {
        /// The parent identity.
        parent: EntityId,
        /// Its status.
        status: EntityStatus,
    },
    /// The containment rule forbids this child/parent kind pair.
    ParentKind {
        /// The child kind.
        child: EntityKind,
        /// The parent kind, or `None` for a root declaration.
        parent: Option<EntityKind>,
    },
    /// An `Interface` declaration carried no pair, or a non-interface did.
    InterfacePairShape {
        /// The declared kind.
        kind: EntityKind,
        /// Whether a pair was supplied.
        supplied: bool,
    },
    /// An interface side is not in the catalog.
    UnknownInterfaceSide {
        /// The unresolved side.
        side: EntityId,
    },
    /// An interface side is not a region or surface, or is inactive.
    InterfaceSideKind {
        /// The offending side.
        side: EntityId,
        /// Its kind.
        kind: EntityKind,
    },
    /// Both interface sides are the same entity.
    SelfInterface {
        /// The repeated side.
        side: EntityId,
    },
    /// An interface side is not contained in the interface's parent.
    InterfaceSideScope {
        /// The offending side.
        side: EntityId,
        /// The interface's declared parent.
        parent: EntityId,
    },
    /// The identical declaration was already made.
    DuplicateDeclaration {
        /// The already-present identity.
        id: EntityId,
        /// Its declaration row.
        row: usize,
    },
    /// Two DIFFERENT declarations derived the same identity.
    IdentityCollision {
        /// The colliding identity.
        id: EntityId,
        /// Bounded preview of the stored declared name.
        existing: String,
        /// Bounded preview of the incoming declared name.
        incoming: String,
    },
    /// A revision re-declared an identity that is already superseded/retired.
    InactiveIdentityRedeclared {
        /// The identity.
        id: EntityId,
        /// Its status.
        status: EntityStatus,
    },
    /// The referenced entity is not in the catalog.
    UnknownEntity {
        /// The unresolved identity.
        id: EntityId,
    },
    /// Automatic import matching found more than one candidate at the
    /// strongest available tier. Never resolved by picking the first.
    AmbiguousImportMatch {
        /// The incoming identity.
        incoming: EntityId,
        /// First candidate, in declaration order.
        first: EntityId,
        /// Second candidate, in declaration order.
        second: EntityId,
        /// Total candidates at this tier.
        total: usize,
        /// The tier at which the ambiguity occurred.
        tier: EvidenceTier,
    },
    /// A predecessor already has a successor.
    AmbiguousSupersession {
        /// The predecessor.
        predecessor: EntityId,
        /// The already-recorded successor.
        existing: EntityId,
        /// The proposed successor.
        proposed: EntityId,
    },
    /// A revision asserted a predecessor that is superseded or retired.
    InactivePredecessor {
        /// The predecessor.
        predecessor: EntityId,
        /// Its status.
        status: EntityStatus,
    },
    /// An import revision carried a non-revision event.
    NonRevisionEvent {
        /// The supplied event.
        event: RebindEvent,
    },
    /// A declared-name lookup matched more than one entity.
    AmbiguousDeclaredName {
        /// Bounded preview of the name.
        name: String,
        /// First match, in declaration order.
        first: EntityId,
        /// Second match, in declaration order.
        second: EntityId,
        /// Total matches.
        total: usize,
    },
    /// A declared-name lookup matched nothing.
    UnknownDeclaredName {
        /// Bounded preview of the name.
        name: String,
    },
    /// A datum reference is not in the catalog.
    UnknownDatum {
        /// The unresolved datum identity.
        id: DatumId,
    },
    /// The identical datum was already declared.
    DuplicateDatum {
        /// The already-present datum identity.
        id: DatumId,
        /// Its declaration row.
        row: usize,
    },
    /// An interface cannot carry a datum feature.
    DatumOwnerKind {
        /// The proposed owner.
        owner: EntityId,
        /// Its kind.
        kind: EntityKind,
    },
    /// Only a part occurrence can be placed.
    PlacementOccurrenceKind {
        /// The proposed occurrence.
        occurrence: EntityId,
        /// Its kind.
        kind: EntityKind,
    },
    /// An as-built placement cited the all-zero registration identity.
    PlacementRegistrationUnbound {
        /// The occurrence whose citation is unbound.
        occurrence: EntityId,
    },
    /// A named admission limit was exceeded.
    CapacityExceeded {
        /// Which resource.
        resource: &'static str,
        /// Requested count.
        requested: usize,
        /// Admitted count.
        limit: usize,
    },
    /// A fallible reservation was refused before mutation.
    AllocationRefused {
        /// Which collection.
        resource: &'static str,
    },
    /// A containment or supersession walk exceeded its hop budget.
    HopBudgetExceeded {
        /// Which walk.
        walk: &'static str,
        /// Where the walk started.
        start: EntityId,
        /// The admitted hop budget.
        limit: usize,
    },
}

impl EntityError {
    /// Stable machine-readable code, shared with the `Violation` shape.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            EntityError::EmptyName { .. } => "entity-empty-name",
            EntityError::NameTooLong { .. } => "entity-name-too-long",
            EntityError::UnknownParent { .. } => "entity-unknown-parent",
            EntityError::InactiveParent { .. } => "entity-inactive-parent",
            EntityError::ParentKind { .. } => "entity-parent-kind",
            EntityError::InterfacePairShape { .. } => "entity-interface-pair-shape",
            EntityError::UnknownInterfaceSide { .. } => "entity-unknown-interface-side",
            EntityError::InterfaceSideKind { .. } => "entity-interface-side-kind",
            EntityError::SelfInterface { .. } => "entity-self-interface",
            EntityError::InterfaceSideScope { .. } => "entity-interface-side-scope",
            EntityError::DuplicateDeclaration { .. } => "entity-duplicate-declaration",
            EntityError::IdentityCollision { .. } => "entity-identity-collision",
            EntityError::InactiveIdentityRedeclared { .. } => "entity-inactive-redeclared",
            EntityError::UnknownEntity { .. } => "entity-unknown",
            EntityError::AmbiguousImportMatch { .. } => "entity-ambiguous-import-match",
            EntityError::AmbiguousSupersession { .. } => "entity-ambiguous-supersession",
            EntityError::InactivePredecessor { .. } => "entity-inactive-predecessor",
            EntityError::NonRevisionEvent { .. } => "entity-non-revision-event",
            EntityError::AmbiguousDeclaredName { .. } => "entity-ambiguous-declared-name",
            EntityError::UnknownDeclaredName { .. } => "entity-unknown-declared-name",
            EntityError::UnknownDatum { .. } => "entity-unknown-datum",
            EntityError::DuplicateDatum { .. } => "entity-duplicate-datum",
            EntityError::DatumOwnerKind { .. } => "entity-datum-owner-kind",
            EntityError::PlacementOccurrenceKind { .. } => "entity-placement-occurrence-kind",
            EntityError::PlacementRegistrationUnbound { .. } => {
                "entity-placement-registration-unbound"
            }
            EntityError::CapacityExceeded { .. } => "entity-capacity-exceeded",
            EntityError::AllocationRefused { .. } => "entity-allocation-refused",
            EntityError::HopBudgetExceeded { .. } => "entity-hop-budget-exceeded",
        }
    }

    /// The ranked repair for this refusal.
    #[must_use]
    pub fn fix(&self) -> String {
        match self {
            EntityError::EmptyName { field, .. } => {
                format!("supply a nonempty {field}; identity is derived from it")
            }
            EntityError::NameTooLong { limit, .. } => {
                format!("shorten the name to at most {limit} bytes")
            }
            EntityError::UnknownParent { .. } => {
                "declare the parent before its children, or bind to the parent's current identity"
                    .to_string()
            }
            EntityError::InactiveParent { .. } => {
                "declare the child under the parent's CURRENT identity (resolve the parent first)"
                    .to_string()
            }
            EntityError::ParentKind { .. } => {
                "follow the containment rule: assembly -> (assembly | part), part -> (region | surface), interface under the assembly or part owning both sides"
                    .to_string()
            }
            EntityError::InterfacePairShape { .. } => {
                "declare interfaces with InterfacePair::ordered/unordered and no other kind with a pair"
                    .to_string()
            }
            EntityError::UnknownInterfaceSide { .. } => {
                "declare both sides before the interface".to_string()
            }
            EntityError::InterfaceSideKind { .. } => {
                "use active region or surface entities as interface sides".to_string()
            }
            EntityError::SelfInterface { .. } => {
                "an interface needs two distinct sides; declare the second side".to_string()
            }
            EntityError::InterfaceSideScope { .. } => {
                "declare the interface under an ancestor that contains both sides".to_string()
            }
            EntityError::DuplicateDeclaration { .. } => {
                "this exact declaration is already in the catalog; reuse its identity".to_string()
            }
            EntityError::IdentityCollision { .. } => {
                "two different declarations derived one digest: report this as a collision and change one declared name to proceed"
                    .to_string()
            }
            EntityError::InactiveIdentityRedeclared { .. } => {
                "do not revive a superseded or retired identity; declare a new one and assert the correspondence"
                    .to_string()
            }
            EntityError::UnknownEntity { .. } => {
                "declare the entity, or resolve the reference to its current identity".to_string()
            }
            EntityError::AmbiguousImportMatch { .. } => {
                "the automatic match is ambiguous; supply Correspondence::Asserted(predecessor) or Correspondence::New for this entity"
                    .to_string()
            }
            EntityError::AmbiguousSupersession { .. } => {
                "one identity may be superseded once; retire the extra revision or assert a different predecessor"
                    .to_string()
            }
            EntityError::InactivePredecessor { .. } => {
                "assert the CURRENT identity as the predecessor".to_string()
            }
            EntityError::NonRevisionEvent { .. } => {
                "use RebindEvent::Import, Remesh, or RevisionMigration for an import revision"
                    .to_string()
            }
            EntityError::AmbiguousDeclaredName { total, .. } => {
                format!(
                    "the name matches {total} entities: reference the intended entity by identity"
                )
            }
            EntityError::UnknownDeclaredName { .. } => {
                "declare the entity before referencing it by name".to_string()
            }
            EntityError::UnknownDatum { .. } => {
                "declare the datum before referencing it".to_string()
            }
            EntityError::DuplicateDatum { .. } => {
                "this exact datum is already in the catalog; reuse its identity".to_string()
            }
            EntityError::DatumOwnerKind { .. } => {
                "attach the datum to the assembly, part, region, or surface that carries the reference feature"
                    .to_string()
            }
            EntityError::PlacementOccurrenceKind { .. } => {
                "place a part occurrence; regions, surfaces, and interfaces move with their part"
                    .to_string()
            }
            EntityError::PlacementRegistrationUnbound { .. } => {
                "cite the identity of a real calibrated registration record, or declare the placement nominal"
                    .to_string()
            }
            EntityError::CapacityExceeded { resource, .. } => {
                format!("raise the explicit {resource} limit in the EntityBudget, or split the model")
            }
            EntityError::AllocationRefused { .. } => {
                "reduce the model size; the allocation was refused before any mutation".to_string()
            }
            EntityError::HopBudgetExceeded { walk, .. } => {
                format!("flatten the {walk} chain or raise its explicit hop budget")
            }
        }
    }

    /// Convert into the crate's structured `Violation` shape.
    #[must_use]
    pub fn into_violation(self) -> Violation {
        Violation {
            code: self.code(),
            what: self.to_string(),
            fix: self.fix(),
        }
    }
}

impl fmt::Display for EntityError {
    #[allow(clippy::too_many_lines)]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EntityError::EmptyName { kind, field } => {
                write!(f, "{kind} declaration has an empty {field}")
            }
            EntityError::NameTooLong {
                kind,
                field,
                bytes,
                limit,
            } => write!(
                f,
                "{kind} {field} is {bytes} bytes, above the admitted {limit}"
            ),
            EntityError::UnknownParent { kind, parent } => {
                write!(f, "{kind} declares unknown parent {parent}")
            }
            EntityError::InactiveParent { parent, status } => {
                write!(f, "parent {parent} is not active ({status:?})")
            }
            EntityError::ParentKind { child, parent } => match parent {
                Some(parent) => write!(f, "{child} cannot be contained by {parent}"),
                None => write!(f, "{child} cannot be declared without a parent"),
            },
            EntityError::InterfacePairShape { kind, supplied } => write!(
                f,
                "{kind} declaration {} an interface pair",
                if *supplied {
                    "must not carry"
                } else {
                    "must carry"
                }
            ),
            EntityError::UnknownInterfaceSide { side } => {
                write!(f, "interface side {side} is not in the catalog")
            }
            EntityError::InterfaceSideKind { side, kind } => {
                write!(
                    f,
                    "interface side {side} is a {kind}, not an active region or surface"
                )
            }
            EntityError::SelfInterface { side } => {
                write!(f, "interface declares {side} on both sides")
            }
            EntityError::InterfaceSideScope { side, parent } => {
                write!(f, "interface side {side} is not contained in {parent}")
            }
            EntityError::DuplicateDeclaration { id, row } => {
                write!(f, "{id} was already declared at row {row}")
            }
            EntityError::IdentityCollision {
                id,
                existing,
                incoming,
            } => write!(
                f,
                "identity collision on {id}: stored declaration {existing} and incoming {incoming} differ"
            ),
            EntityError::InactiveIdentityRedeclared { id, status } => {
                write!(f, "{id} is {status:?} and cannot be redeclared")
            }
            EntityError::UnknownEntity { id } => write!(f, "{id} is not in the catalog"),
            EntityError::AmbiguousImportMatch {
                incoming,
                first,
                second,
                total,
                tier,
            } => write!(
                f,
                "import match for {incoming} is ambiguous at tier {tier}: {total} candidates including {first} and {second}"
            ),
            EntityError::AmbiguousSupersession {
                predecessor,
                existing,
                proposed,
            } => write!(
                f,
                "{predecessor} is already superseded by {existing}; {proposed} would be a second successor"
            ),
            EntityError::InactivePredecessor {
                predecessor,
                status,
            } => write!(f, "asserted predecessor {predecessor} is {status:?}"),
            EntityError::NonRevisionEvent { event } => {
                write!(f, "{event} is not an import-revision event")
            }
            EntityError::AmbiguousDeclaredName {
                name,
                first,
                second,
                total,
            } => write!(
                f,
                "declared name {name} matches {total} entities including {first} and {second}"
            ),
            EntityError::UnknownDeclaredName { name } => {
                write!(f, "declared name {name} matches no entity")
            }
            EntityError::UnknownDatum { id } => write!(f, "{id} is not in the catalog"),
            EntityError::DuplicateDatum { id, row } => {
                write!(f, "{id} was already declared at row {row}")
            }
            EntityError::DatumOwnerKind { owner, kind } => {
                write!(f, "{owner} is a {kind} and cannot carry a datum feature")
            }
            EntityError::PlacementOccurrenceKind { occurrence, kind } => {
                write!(
                    f,
                    "{occurrence} is a {kind}, not a placeable part occurrence"
                )
            }
            EntityError::PlacementRegistrationUnbound { occurrence } => {
                write!(
                    f,
                    "the as-built placement of {occurrence} cites the all-zero registration identity, which names no artifact"
                )
            }
            EntityError::CapacityExceeded {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "{resource} request {requested} exceeds the admitted limit {limit}"
            ),
            EntityError::AllocationRefused { resource } => {
                write!(f, "{resource} reservation was refused before mutation")
            }
            EntityError::HopBudgetExceeded { walk, start, limit } => {
                write!(f, "{walk} walk from {start} exceeded {limit} hops")
            }
        }
    }
}

impl core::error::Error for EntityError {}

// ---------------------------------------------------------------------------
// Resolution
// ---------------------------------------------------------------------------

/// What kind an entity reference is required to designate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum KindExpectation {
    /// Exactly one kind.
    Exact(EntityKind),
    /// A volumetric domain: a part or a region.
    Domain,
    /// A boundary: a surface or an interface.
    Boundary,
    /// Any entity kind.
    Any,
}

impl KindExpectation {
    /// Whether `kind` satisfies this expectation.
    #[must_use]
    pub const fn admits(self, kind: EntityKind) -> bool {
        match self {
            KindExpectation::Exact(expected) => expected.tag() == kind.tag(),
            KindExpectation::Domain => matches!(kind, EntityKind::Part | EntityKind::Region),
            KindExpectation::Boundary => {
                matches!(kind, EntityKind::Surface | EntityKind::Interface)
            }
            KindExpectation::Any => true,
        }
    }

    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            KindExpectation::Exact(kind) => kind.label(),
            KindExpectation::Domain => "part|region",
            KindExpectation::Boundary => "surface|interface",
            KindExpectation::Any => "any",
        }
    }
}

impl fmt::Display for KindExpectation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.label())
    }
}

/// A typed reference to an entity: an identity plus the kind the referring
/// site requires.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EntityRef {
    target: EntityId,
    expect: KindExpectation,
}

impl EntityRef {
    /// A reference to `target` that must designate a kind satisfying `expect`.
    #[must_use]
    pub const fn new(target: EntityId, expect: KindExpectation) -> Self {
        EntityRef { target, expect }
    }

    /// The referenced identity as authored.
    #[must_use]
    pub const fn target(self) -> EntityId {
        self.target
    }

    /// The required kind.
    #[must_use]
    pub const fn expect(self) -> KindExpectation {
        self.expect
    }
}

/// A successful reference resolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Resolution {
    requested: EntityId,
    current: EntityId,
    hops: usize,
    tier: EvidenceTier,
}

impl Resolution {
    /// The identity the caller asked for.
    #[must_use]
    pub const fn requested(self) -> EntityId {
        self.requested
    }

    /// The current identity the reference resolves to.
    #[must_use]
    pub const fn current(self) -> EntityId {
        self.current
    }

    /// Supersession hops followed (0 when the identity is current).
    #[must_use]
    pub const fn hops(self) -> usize {
        self.hops
    }

    /// The WEAKEST evidence tier along the followed chain.
    ///
    /// A chain of one asserted hop and one content-matched hop resolves at
    /// `Asserted`: composition never launders a weak link into a strong one.
    #[must_use]
    pub const fn tier(self) -> EvidenceTier {
        self.tier
    }
}

/// Why a reference did not resolve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolutionFault {
    /// No entity with that identity.
    Dangling {
        /// The unresolved identity.
        id: EntityId,
    },
    /// The chain ends at a retired entity.
    Retired {
        /// The retired identity.
        id: EntityId,
        /// The retirement receipt sequence.
        receipt: u64,
    },
    /// The supersession chain exceeded the hop budget.
    DepthExceeded {
        /// Where the walk started.
        id: EntityId,
        /// The admitted hop budget.
        limit: usize,
    },
    /// The resolved entity is the wrong kind for the referring site.
    KindMismatch {
        /// The resolved identity.
        id: EntityId,
        /// Its actual kind.
        actual: EntityKind,
        /// The required kind.
        expect: KindExpectation,
    },
}

impl ResolutionFault {
    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            ResolutionFault::Dangling { .. } => "entity-dangling-reference",
            ResolutionFault::Retired { .. } => "entity-retired-reference",
            ResolutionFault::DepthExceeded { .. } => "entity-supersession-depth",
            ResolutionFault::KindMismatch { .. } => "entity-kind-mismatch",
        }
    }

    /// The ranked repair.
    #[must_use]
    pub fn fix(self) -> String {
        match self {
            ResolutionFault::Dangling { .. } => {
                "declare the entity, or re-import the revision that introduced it so a supersession receipt exists"
                    .to_string()
            }
            ResolutionFault::Retired { .. } => {
                "the entity was retired by a complete revision: rebind this site to a current entity or restore the entity in a new revision"
                    .to_string()
            }
            ResolutionFault::DepthExceeded { limit, .. } => format!(
                "collapse the supersession chain or raise the {limit}-hop resolution budget"
            ),
            ResolutionFault::KindMismatch { expect, .. } => {
                format!("bind this site to a {expect} entity")
            }
        }
    }
}

impl fmt::Display for ResolutionFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ResolutionFault::Dangling { id } => write!(f, "{id} does not resolve to any entity"),
            ResolutionFault::Retired { id, receipt } => {
                write!(f, "{id} was retired by receipt {receipt}")
            }
            ResolutionFault::DepthExceeded { id, limit } => {
                write!(f, "supersession chain from {id} exceeded {limit} hops")
            }
            ResolutionFault::KindMismatch { id, actual, expect } => {
                write!(f, "{id} is a {actual} where {expect} is required")
            }
        }
    }
}

/// The outcome of a name lookup: strings are not identities, and this type is
/// how that shows up in the API.
///
/// Every variant quantifies over ACTIVE entities only. A name carried by one
/// active entity and one superseded entity is `Unique`; the superseded entity
/// is still reachable by identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NameLookup {
    /// No active entity carries the name.
    Missing,
    /// Exactly one active entity carries the name.
    Unique(EntityId),
    /// More than one active entity carries the name.
    Ambiguous {
        /// First match in declaration order.
        first: EntityId,
        /// Second match in declaration order.
        second: EntityId,
        /// Total matches.
        total: usize,
    },
}

// ---------------------------------------------------------------------------
// Import revisions
// ---------------------------------------------------------------------------

/// What an import revision claims about entities it does NOT contain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportScope {
    /// The revision is a partial update: omitted entities are untouched.
    Partial,
    /// The revision is the complete content of `root`: active strict
    /// descendants of `root` that it does not contain are RETIRED, each with a
    /// receipt. The root itself is never auto-retired.
    Complete {
        /// The subtree the revision claims to cover completely.
        root: EntityId,
    },
}

/// How an incoming declaration corresponds to an existing identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Correspondence {
    /// Let the catalog match on fingerprint, then on declared path. An
    /// ambiguous match is a typed refusal, never a first-match guess.
    Auto,
    /// The caller asserts the predecessor. This crate verifies nothing beyond
    /// the predecessor being a known, active identity.
    Asserted(EntityId),
    /// The caller asserts there is no predecessor.
    New,
}

/// One entity in an import revision.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportedEntity {
    /// The incoming declaration.
    pub declaration: EntityDeclaration,
    /// How it corresponds to the existing catalog.
    pub correspondence: Correspondence,
}

/// A geometry import, re-mesh, or revision migration applied atomically.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportRevision {
    /// Revision label, retained for reports.
    pub label: String,
    /// Which rebinding event this revision records.
    pub event: RebindEvent,
    /// What the revision claims about entities it omits.
    pub scope: ImportScope,
    /// The declarations, parents before children.
    pub entities: Vec<ImportedEntity>,
}

/// What one entity in a revision did to the catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ImportStep {
    /// The identity already existed and did not change.
    Unchanged {
        /// The identity.
        id: EntityId,
    },
    /// A new identity was declared with no predecessor.
    Declared {
        /// The new identity.
        id: EntityId,
    },
    /// A new identity superseded an existing one.
    Superseded {
        /// The new identity.
        id: EntityId,
        /// The superseded identity.
        predecessor: EntityId,
        /// Why the correspondence is claimed.
        basis: MatchBasis,
    },
    /// An active entity was retired by a complete revision.
    Retired {
        /// The retired identity.
        id: EntityId,
    },
}

/// The deterministic result of applying one revision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImportOutcome {
    steps: Vec<ImportStep>,
    first_receipt: u64,
    receipts: u64,
}

impl ImportOutcome {
    /// Per-entity steps, in revision order followed by retirements in
    /// declaration order.
    #[must_use]
    pub fn steps(&self) -> &[ImportStep] {
        &self.steps
    }

    /// Sequence of the first receipt this revision appended.
    #[must_use]
    pub const fn first_receipt(&self) -> u64 {
        self.first_receipt
    }

    /// How many receipts this revision appended.
    #[must_use]
    pub const fn receipts(&self) -> u64 {
        self.receipts
    }
}

// ---------------------------------------------------------------------------
// Catalog
// ---------------------------------------------------------------------------

/// The entity catalog: identities, their lifecycle, and the receipt log.
#[derive(Debug, Clone, PartialEq)]
pub struct EntityCatalog {
    budget: EntityBudget,
    entities: Vec<Entity>,
    index: Vec<(EntityId, usize)>,
    receipts: Vec<IdentityReceipt>,
    chain: ContentHash,
    datums: Vec<Datum>,
    datum_index: Vec<(DatumId, usize)>,
    tolerances: Vec<Tolerance>,
    placements: Vec<Placement>,
}

impl Default for EntityCatalog {
    fn default() -> Self {
        Self::new()
    }
}

impl EntityCatalog {
    /// An empty catalog under [`DEFAULT_ENTITY_BUDGET`].
    #[must_use]
    pub fn new() -> Self {
        Self::with_budget(DEFAULT_ENTITY_BUDGET)
    }

    /// An empty catalog under an explicit budget.
    #[must_use]
    pub fn with_budget(budget: EntityBudget) -> Self {
        EntityCatalog {
            budget,
            entities: Vec::new(),
            index: Vec::new(),
            receipts: Vec::new(),
            chain: fs_blake3::hash_domain(IDENTITY_RECEIPT_DOMAIN, RECEIPT_CHAIN_GENESIS),
            datums: Vec::new(),
            datum_index: Vec::new(),
            tolerances: Vec::new(),
            placements: Vec::new(),
        }
    }

    /// The admission limits in force.
    #[must_use]
    pub const fn budget(&self) -> EntityBudget {
        self.budget
    }

    /// Every entity in declaration order.
    #[must_use]
    pub fn entities(&self) -> &[Entity] {
        &self.entities
    }

    /// Every receipt in append order.
    #[must_use]
    pub fn receipts(&self) -> &[IdentityReceipt] {
        &self.receipts
    }

    /// The current receipt-chain root.
    ///
    /// Pin this outside the catalog to detect tail truncation; on its own the
    /// chain proves only that the retained log is internally consistent.
    #[must_use]
    pub const fn receipt_root(&self) -> ContentHash {
        self.chain
    }

    /// Every declared datum in declaration order.
    #[must_use]
    pub fn datums(&self) -> &[Datum] {
        &self.datums
    }

    /// Every declared tolerance in declaration order.
    #[must_use]
    pub fn tolerances(&self) -> &[Tolerance] {
        &self.tolerances
    }

    /// Every declared placement in declaration order.
    #[must_use]
    pub fn placements(&self) -> &[Placement] {
        &self.placements
    }

    /// Recompute the receipt chain and verify every link.
    #[must_use]
    pub fn verify_receipts(&self) -> bool {
        let mut previous = fs_blake3::hash_domain(IDENTITY_RECEIPT_DOMAIN, RECEIPT_CHAIN_GENESIS);
        for (position, receipt) in self.receipts.iter().enumerate() {
            if receipt.sequence != position as u64
                || receipt.previous != previous
                || !receipt.verifies()
            {
                return false;
            }
            previous = receipt.digest;
        }
        previous == self.chain
    }

    fn row_of(&self, id: EntityId) -> Option<usize> {
        let position = self.index.partition_point(|(candidate, _)| *candidate < id);
        self.index
            .get(position)
            .and_then(|(candidate, row)| (*candidate == id).then_some(*row))
    }

    /// The entity with this identity, whatever its status.
    #[must_use]
    pub fn get(&self, id: EntityId) -> Option<&Entity> {
        self.row_of(id).map(|row| &self.entities[row])
    }

    /// Whether the identity is present and active.
    #[must_use]
    pub fn is_active(&self, id: EntityId) -> bool {
        self.get(id)
            .is_some_and(|entity| entity.status == EntityStatus::Active)
    }

    /// Look up entities by their identity-bearing declared name.
    ///
    /// A name is not an identity: this returns [`NameLookup::Ambiguous`]
    /// rather than a first match whenever more than one entity carries it.
    #[must_use]
    pub fn lookup_declared_name(&self, name: &str) -> NameLookup {
        self.lookup_by(|entity| entity.declared_name() == name)
    }

    /// Look up entities by their CURRENT display name.
    #[must_use]
    pub fn lookup_display_name(&self, name: &str) -> NameLookup {
        self.lookup_by(|entity| entity.display_name() == name)
    }

    fn lookup_by(&self, mut predicate: impl FnMut(&Entity) -> bool) -> NameLookup {
        let mut first = None;
        let mut second = None;
        let mut total = 0usize;
        for entity in &self.entities {
            if entity.status != EntityStatus::Active || !predicate(entity) {
                continue;
            }
            total += 1;
            if first.is_none() {
                first = Some(entity.id);
            } else if second.is_none() {
                second = Some(entity.id);
            }
        }
        match (first, second) {
            (None, _) => NameLookup::Missing,
            (Some(first), None) => NameLookup::Unique(first),
            (Some(first), Some(second)) => NameLookup::Ambiguous {
                first,
                second,
                total,
            },
        }
    }

    /// Resolve a declared name to exactly one identity, or refuse.
    ///
    /// # Errors
    /// [`EntityError::UnknownDeclaredName`] when nothing matches and
    /// [`EntityError::AmbiguousDeclaredName`] when more than one entity does.
    pub fn resolve_declared_name(&self, name: &str) -> Result<EntityId, EntityError> {
        match self.lookup_declared_name(name) {
            NameLookup::Missing => Err(EntityError::UnknownDeclaredName {
                name: preview(name),
            }),
            NameLookup::Unique(id) => Ok(id),
            NameLookup::Ambiguous {
                first,
                second,
                total,
            } => Err(EntityError::AmbiguousDeclaredName {
                name: preview(name),
                first,
                second,
                total,
            }),
        }
    }

    /// Resolve a typed reference through supersession.
    ///
    /// # Errors
    /// Returns a [`ResolutionFault`] for dangling, retired, over-deep, and
    /// wrong-kind references.
    pub fn resolve(&self, reference: EntityRef) -> Result<Resolution, ResolutionFault> {
        let requested = reference.target();
        let mut current = requested;
        let mut hops = 0usize;
        let mut tier = EvidenceTier::Identical;
        loop {
            let Some(entity) = self.get(current) else {
                return Err(ResolutionFault::Dangling { id: current });
            };
            match entity.status {
                EntityStatus::Active => break,
                EntityStatus::Retired { receipt } => {
                    return Err(ResolutionFault::Retired {
                        id: current,
                        receipt,
                    });
                }
                EntityStatus::Superseded { successor, receipt } => {
                    if hops >= self.budget.max_supersession_hops {
                        return Err(ResolutionFault::DepthExceeded {
                            id: requested,
                            limit: self.budget.max_supersession_hops,
                        });
                    }
                    let hop_tier = self
                        .receipts
                        .get(receipt as usize)
                        .and_then(|entry| entry.basis.tier())
                        .unwrap_or(EvidenceTier::Asserted);
                    tier = tier.min(hop_tier);
                    current = successor;
                    hops += 1;
                }
            }
        }
        let actual = current.kind();
        if !reference.expect().admits(actual) {
            return Err(ResolutionFault::KindMismatch {
                id: current,
                actual,
                expect: reference.expect(),
            });
        }
        Ok(Resolution {
            requested,
            current,
            hops,
            tier,
        })
    }

    fn contains_in_subtree(&self, node: EntityId, ancestor: EntityId) -> Result<bool, EntityError> {
        let mut current = node;
        let mut hops = 0usize;
        loop {
            if current == ancestor {
                return Ok(true);
            }
            let Some(entity) = self.get(current) else {
                return Ok(false);
            };
            let Some(parent) = entity.parent() else {
                return Ok(false);
            };
            if hops >= self.budget.max_hierarchy_depth {
                return Err(EntityError::HopBudgetExceeded {
                    walk: "containment",
                    start: node,
                    limit: self.budget.max_hierarchy_depth,
                });
            }
            current = parent;
            hops += 1;
        }
    }

    fn admit(&self, declaration: &EntityDeclaration) -> Result<EntityId, EntityError> {
        let kind = declaration.kind;
        for (field, name) in [
            ("declared name", declaration.declared_name.as_str()),
            ("display name", declaration.display_name.as_str()),
        ] {
            if name.is_empty() {
                return Err(EntityError::EmptyName { kind, field });
            }
            if name.len() > self.budget.max_name_bytes {
                return Err(EntityError::NameTooLong {
                    kind,
                    field,
                    bytes: name.len(),
                    limit: self.budget.max_name_bytes,
                });
            }
        }
        let parent_kind = match declaration.parent {
            None => None,
            Some(parent) => {
                let entity = self
                    .get(parent)
                    .ok_or(EntityError::UnknownParent { kind, parent })?;
                if entity.status != EntityStatus::Active {
                    return Err(EntityError::InactiveParent {
                        parent,
                        status: entity.status,
                    });
                }
                Some(entity.kind())
            }
        };
        if !kind.admits_parent(parent_kind) {
            return Err(EntityError::ParentKind {
                child: kind,
                parent: parent_kind,
            });
        }
        match (kind, declaration.pair) {
            (EntityKind::Interface, None) => {
                return Err(EntityError::InterfacePairShape {
                    kind,
                    supplied: false,
                });
            }
            (EntityKind::Interface, Some(pair)) => {
                self.admit_pair(pair, declaration.parent)?;
            }
            (_, Some(_)) => {
                return Err(EntityError::InterfacePairShape {
                    kind,
                    supplied: true,
                });
            }
            (_, None) => {}
        }
        if self.entities.len() >= self.budget.max_entities {
            return Err(EntityError::CapacityExceeded {
                resource: "entities",
                requested: self.entities.len() + 1,
                limit: self.budget.max_entities,
            });
        }
        if self.receipts.len() >= self.budget.max_receipts {
            return Err(EntityError::CapacityExceeded {
                resource: "receipts",
                requested: self.receipts.len() + 1,
                limit: self.budget.max_receipts,
            });
        }
        let id = declaration.identity();
        if let Some(row) = self.row_of(id) {
            let existing = &self.entities[row];
            if existing.declaration.same_identity_preimage(declaration) {
                return Err(EntityError::DuplicateDeclaration { id, row });
            }
            return Err(EntityError::IdentityCollision {
                id,
                existing: preview(existing.declared_name()),
                incoming: preview(&declaration.declared_name),
            });
        }
        Ok(id)
    }

    fn admit_pair(&self, pair: InterfacePair, parent: Option<EntityId>) -> Result<(), EntityError> {
        if pair.from == pair.to {
            return Err(EntityError::SelfInterface { side: pair.from });
        }
        for side in [pair.from, pair.to] {
            let entity = self
                .get(side)
                .ok_or(EntityError::UnknownInterfaceSide { side })?;
            if !entity.kind().admits_interface_side() || entity.status != EntityStatus::Active {
                return Err(EntityError::InterfaceSideKind {
                    side,
                    kind: entity.kind(),
                });
            }
            if let Some(parent) = parent
                && !self.contains_in_subtree(side, parent)?
            {
                return Err(EntityError::InterfaceSideScope { side, parent });
            }
        }
        Ok(())
    }

    fn reserve_declaration_slots(&mut self) -> Result<(), EntityError> {
        self.entities
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "entities",
            })?;
        self.index
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused { resource: "index" })?;
        self.receipts
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "receipts",
            })
    }

    fn insert_entity(&mut self, id: EntityId, declaration: EntityDeclaration) -> usize {
        let row = self.entities.len();
        let display_name = declaration.display_name.clone();
        self.entities.push(Entity {
            id,
            declaration,
            display_name,
            status: EntityStatus::Active,
            row,
        });
        let position = self.index.partition_point(|(candidate, _)| *candidate < id);
        self.index.insert(position, (id, row));
        row
    }

    fn append_receipt(
        &mut self,
        event: RebindEvent,
        subject: EntityId,
        predecessor: Option<EntityId>,
        basis: MatchBasis,
        display_before: Option<String>,
        display_after: Option<String>,
    ) -> u64 {
        let sequence = self.receipts.len() as u64;
        let previous = self.chain;
        let digest = IdentityReceipt::compute_digest(
            previous,
            sequence,
            event,
            subject,
            predecessor,
            basis,
            display_before.as_deref(),
            display_after.as_deref(),
        );
        self.chain = digest;
        self.receipts.push(IdentityReceipt {
            sequence,
            event,
            subject,
            predecessor,
            basis,
            display_before,
            display_after,
            previous,
            digest,
        });
        sequence
    }

    /// Declare a new entity and append its declaration receipt.
    ///
    /// # Errors
    /// Returns a typed [`EntityError`] for empty/oversized names, unknown or
    /// inactive parents, containment violations, malformed interface pairs,
    /// duplicate declarations, digest collisions, and capacity refusals.
    /// Nothing is mutated when the declaration is refused.
    pub fn declare(&mut self, declaration: EntityDeclaration) -> Result<EntityId, EntityError> {
        let id = self.admit(&declaration)?;
        self.reserve_declaration_slots()?;
        let (event, basis) = if declaration.legacy {
            (RebindEvent::LegacyMigration, MatchBasis::LegacyName)
        } else {
            (RebindEvent::Declaration, MatchBasis::NewIdentity)
        };
        let display = declaration.display_name.clone();
        self.insert_entity(id, declaration);
        self.append_receipt(event, id, None, basis, None, Some(display));
        Ok(id)
    }

    /// Declare an entity, or return the existing identity when the SAME
    /// declaration preimage is already present. Used by mechanical migration,
    /// which must be idempotent; authoring uses [`EntityCatalog::declare`],
    /// which refuses duplicates.
    ///
    /// # Errors
    /// Same as [`EntityCatalog::declare`], minus duplicate declarations.
    pub fn declare_or_existing(
        &mut self,
        declaration: EntityDeclaration,
    ) -> Result<EntityId, EntityError> {
        match self.declare(declaration) {
            Err(EntityError::DuplicateDeclaration { id, .. }) => Ok(id),
            other => other,
        }
    }

    /// Change an entity's DISPLAY name. The identity does not move, so every
    /// reference stays bound; the rename is recorded as a receipt.
    ///
    /// # Errors
    /// [`EntityError::UnknownEntity`] for an unknown identity,
    /// [`EntityError::InactiveIdentityRedeclared`] for a superseded or retired
    /// one, and the name-shape refusals.
    pub fn rename(&mut self, id: EntityId, display_name: &str) -> Result<u64, EntityError> {
        let Some(row) = self.row_of(id) else {
            return Err(EntityError::UnknownEntity { id });
        };
        let status = self.entities[row].status;
        if status != EntityStatus::Active {
            return Err(EntityError::InactiveIdentityRedeclared { id, status });
        }
        if display_name.is_empty() {
            return Err(EntityError::EmptyName {
                kind: id.kind(),
                field: "display name",
            });
        }
        if display_name.len() > self.budget.max_name_bytes {
            return Err(EntityError::NameTooLong {
                kind: id.kind(),
                field: "display name",
                bytes: display_name.len(),
                limit: self.budget.max_name_bytes,
            });
        }
        if self.receipts.len() >= self.budget.max_receipts {
            return Err(EntityError::CapacityExceeded {
                resource: "receipts",
                requested: self.receipts.len() + 1,
                limit: self.budget.max_receipts,
            });
        }
        self.receipts
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "receipts",
            })?;
        let before = self.entities[row].display_name.clone();
        self.entities[row].display_name = display_name.to_string();
        Ok(self.append_receipt(
            RebindEvent::Rename,
            id,
            None,
            MatchBasis::Unchanged,
            Some(before),
            Some(display_name.to_string()),
        ))
    }

    fn supersede(
        &mut self,
        predecessor: EntityId,
        successor: EntityId,
        event: RebindEvent,
        basis: MatchBasis,
    ) -> Result<u64, EntityError> {
        let Some(row) = self.row_of(predecessor) else {
            return Err(EntityError::UnknownEntity { id: predecessor });
        };
        match self.entities[row].status {
            EntityStatus::Active => {}
            EntityStatus::Superseded {
                successor: existing,
                ..
            } => {
                return Err(EntityError::AmbiguousSupersession {
                    predecessor,
                    existing,
                    proposed: successor,
                });
            }
            status @ EntityStatus::Retired { .. } => {
                return Err(EntityError::InactivePredecessor {
                    predecessor,
                    status,
                });
            }
        }
        if self.receipts.len() >= self.budget.max_receipts {
            return Err(EntityError::CapacityExceeded {
                resource: "receipts",
                requested: self.receipts.len() + 1,
                limit: self.budget.max_receipts,
            });
        }
        self.receipts
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "receipts",
            })?;
        let receipt = self.append_receipt(event, successor, Some(predecessor), basis, None, None);
        self.entities[row].status = EntityStatus::Superseded { successor, receipt };
        Ok(receipt)
    }

    fn retire(&mut self, id: EntityId, event: RebindEvent) -> Result<u64, EntityError> {
        let Some(row) = self.row_of(id) else {
            return Err(EntityError::UnknownEntity { id });
        };
        let status = self.entities[row].status;
        if status != EntityStatus::Active {
            return Err(EntityError::InactiveIdentityRedeclared { id, status });
        }
        if self.receipts.len() >= self.budget.max_receipts {
            return Err(EntityError::CapacityExceeded {
                resource: "receipts",
                requested: self.receipts.len() + 1,
                limit: self.budget.max_receipts,
            });
        }
        self.receipts
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "receipts",
            })?;
        let receipt = self.append_receipt(event, id, None, MatchBasis::Absent, None, None);
        self.entities[row].status = EntityStatus::Retired { receipt };
        Ok(receipt)
    }

    /// Retire an entity explicitly (a deleted part in a revision).
    ///
    /// # Errors
    /// [`EntityError::UnknownEntity`] or
    /// [`EntityError::InactiveIdentityRedeclared`].
    pub fn retire_entity(&mut self, id: EntityId) -> Result<u64, EntityError> {
        self.retire(id, RebindEvent::Retirement)
    }

    fn auto_match(
        &self,
        declaration: &EntityDeclaration,
    ) -> Result<Option<(EntityId, MatchBasis)>, EntityError> {
        let mut content: Vec<EntityId> = Vec::new();
        let mut path: Vec<EntityId> = Vec::new();
        for entity in &self.entities {
            if entity.status != EntityStatus::Active || entity.kind() != declaration.kind {
                continue;
            }
            let fingerprints_match = match (entity.fingerprint(), declaration.fingerprint) {
                (Some(existing), Some(incoming)) => existing == incoming,
                _ => false,
            };
            if fingerprints_match {
                push_checked(&mut content, entity.id, "import-match-candidates")?;
                continue;
            }
            if entity.declared_name() == declaration.declared_name
                && self.parent_corresponds(entity.parent(), declaration.parent)
            {
                push_checked(&mut path, entity.id, "import-match-candidates")?;
            }
        }
        for (candidates, basis, tier) in [
            (
                &content,
                MatchBasis::GeometryFingerprint,
                EvidenceTier::ContentMatched,
            ),
            (&path, MatchBasis::DeclaredPath, EvidenceTier::PathMatched),
        ] {
            match candidates.len() {
                0 => {}
                1 => return Ok(Some((candidates[0], basis))),
                total => {
                    return Err(EntityError::AmbiguousImportMatch {
                        incoming: declaration.identity(),
                        first: candidates[0],
                        second: candidates[1],
                        total,
                        tier,
                    });
                }
            }
        }
        Ok(None)
    }

    /// Whether a candidate's parent corresponds to an incoming parent, either
    /// exactly or through supersessions already applied in this revision. This
    /// is what lets a renamed sub-assembly keep its children's path match.
    fn parent_corresponds(
        &self,
        candidate_parent: Option<EntityId>,
        incoming_parent: Option<EntityId>,
    ) -> bool {
        match (candidate_parent, incoming_parent) {
            (None, None) => true,
            (Some(candidate), Some(incoming)) => {
                candidate == incoming
                    || self
                        .resolve(EntityRef::new(candidate, KindExpectation::Any))
                        .is_ok_and(|resolution| resolution.current() == incoming)
            }
            _ => false,
        }
    }

    /// Apply one import/re-mesh/revision-migration revision ATOMICALLY.
    ///
    /// Either every step is applied and receipted, or the catalog is left
    /// exactly as it was and a typed [`EntityError`] is returned.
    ///
    /// # Errors
    /// Any declaration refusal, an ambiguous automatic match, an asserted
    /// predecessor that is unknown or inactive, a revision that re-declares a
    /// superseded identity, or a capacity refusal.
    pub fn apply_import(
        &mut self,
        revision: &ImportRevision,
    ) -> Result<ImportOutcome, EntityError> {
        if !revision.event.is_revision_event() {
            return Err(EntityError::NonRevisionEvent {
                event: revision.event,
            });
        }
        let mut working = self.clone();
        let outcome = working.apply_import_in_place(revision)?;
        *self = working;
        Ok(outcome)
    }

    fn apply_import_in_place(
        &mut self,
        revision: &ImportRevision,
    ) -> Result<ImportOutcome, EntityError> {
        let first_receipt = self.receipts.len() as u64;
        let mut steps: Vec<ImportStep> = Vec::new();
        let mut present: Vec<EntityId> = Vec::new();
        for item in &revision.entities {
            let incoming = item.declaration.identity();
            if let Some(row) = self.row_of(incoming) {
                let existing = &self.entities[row];
                // The digest hit must be an honest re-declaration of the same
                // entity, not a different declaration aliasing onto it.
                if !existing
                    .declaration
                    .same_identity_preimage(&item.declaration)
                {
                    return Err(EntityError::IdentityCollision {
                        id: incoming,
                        existing: preview(existing.declared_name()),
                        incoming: preview(item.declaration.declared_name()),
                    });
                }
                let status = existing.status;
                if status != EntityStatus::Active {
                    return Err(EntityError::InactiveIdentityRedeclared {
                        id: incoming,
                        status,
                    });
                }
                if self.receipts.len() >= self.budget.max_receipts {
                    return Err(EntityError::CapacityExceeded {
                        resource: "receipts",
                        requested: self.receipts.len() + 1,
                        limit: self.budget.max_receipts,
                    });
                }
                self.receipts
                    .try_reserve(1)
                    .map_err(|_| EntityError::AllocationRefused {
                        resource: "receipts",
                    })?;
                self.append_receipt(
                    revision.event,
                    incoming,
                    None,
                    MatchBasis::Unchanged,
                    None,
                    None,
                );
                push_checked(&mut steps, ImportStep::Unchanged { id: incoming }, "steps")?;
                push_checked(&mut present, incoming, "import-present")?;
                continue;
            }
            let predecessor = match item.correspondence {
                Correspondence::New => None,
                Correspondence::Asserted(predecessor) => {
                    let entity = self
                        .get(predecessor)
                        .ok_or(EntityError::UnknownEntity { id: predecessor })?;
                    if entity.status != EntityStatus::Active {
                        return Err(EntityError::InactivePredecessor {
                            predecessor,
                            status: entity.status,
                        });
                    }
                    Some((predecessor, MatchBasis::Asserted))
                }
                Correspondence::Auto => self.auto_match(&item.declaration)?,
            };
            let id = self.declare_imported(item.declaration.clone(), revision.event)?;
            match predecessor {
                None => push_checked(&mut steps, ImportStep::Declared { id }, "steps")?,
                Some((predecessor, basis)) => {
                    self.supersede(predecessor, id, revision.event, basis)?;
                    push_checked(
                        &mut steps,
                        ImportStep::Superseded {
                            id,
                            predecessor,
                            basis,
                        },
                        "steps",
                    )?;
                }
            }
            push_checked(&mut present, id, "import-present")?;
        }
        if let ImportScope::Complete { root } = revision.scope {
            if self.get(root).is_none() {
                return Err(EntityError::UnknownEntity { id: root });
            }
            let mut retire_rows: Vec<EntityId> = Vec::new();
            for entity in &self.entities {
                if entity.status != EntityStatus::Active
                    || entity.id == root
                    || present.contains(&entity.id)
                {
                    continue;
                }
                if self.contains_in_subtree(entity.id, root)? {
                    push_checked(&mut retire_rows, entity.id, "retirements")?;
                }
            }
            for id in retire_rows {
                self.retire(id, revision.event)?;
                push_checked(&mut steps, ImportStep::Retired { id }, "steps")?;
            }
        }
        let receipts = self.receipts.len() as u64 - first_receipt;
        Ok(ImportOutcome {
            steps,
            first_receipt,
            receipts,
        })
    }

    fn declare_imported(
        &mut self,
        declaration: EntityDeclaration,
        event: RebindEvent,
    ) -> Result<EntityId, EntityError> {
        let id = self.admit(&declaration)?;
        self.reserve_declaration_slots()?;
        let display = declaration.display_name.clone();
        self.insert_entity(id, declaration);
        self.append_receipt(
            event,
            id,
            None,
            MatchBasis::NewIdentity,
            None,
            Some(display),
        );
        Ok(id)
    }

    // -- datums, tolerances, placements -------------------------------------

    /// Declare a datum feature on an entity.
    ///
    /// # Errors
    /// Unknown/inactive owner, an owner kind that cannot carry a datum, an
    /// unknown or repeated datum reference, a name refusal, or a capacity
    /// refusal.
    pub fn declare_datum(
        &mut self,
        owner: EntityId,
        declared_name: &str,
        feature: DatumFeature,
        references: &[DatumId],
    ) -> Result<DatumId, EntityError> {
        if declared_name.is_empty() {
            return Err(EntityError::EmptyName {
                kind: owner.kind(),
                field: "datum name",
            });
        }
        if declared_name.len() > self.budget.max_name_bytes {
            return Err(EntityError::NameTooLong {
                kind: owner.kind(),
                field: "datum name",
                bytes: declared_name.len(),
                limit: self.budget.max_name_bytes,
            });
        }
        let entity = self
            .get(owner)
            .ok_or(EntityError::UnknownEntity { id: owner })?;
        if entity.status != EntityStatus::Active {
            return Err(EntityError::InactiveParent {
                parent: owner,
                status: entity.status,
            });
        }
        if entity.kind() == EntityKind::Interface {
            return Err(EntityError::DatumOwnerKind {
                owner,
                kind: entity.kind(),
            });
        }
        if references.len() > MAX_DATUM_FRAME_LEN {
            return Err(EntityError::CapacityExceeded {
                resource: "datum-references",
                requested: references.len(),
                limit: MAX_DATUM_FRAME_LEN,
            });
        }
        for reference in references {
            if self.datum_row_of(*reference).is_none() {
                return Err(EntityError::UnknownDatum { id: *reference });
            }
        }
        if self.datums.len() >= self.budget.max_datums {
            return Err(EntityError::CapacityExceeded {
                resource: "datums",
                requested: self.datums.len() + 1,
                limit: self.budget.max_datums,
            });
        }
        let mut hasher = DomainHasher::new(DATUM_ID_DOMAIN);
        absorb_entity_id(&mut hasher, Some(owner));
        absorb_bytes(&mut hasher, declared_name.as_bytes());
        hasher.update(&[feature.tag()]);
        hasher.update(&(references.len() as u64).to_le_bytes());
        for reference in references {
            hasher.update(reference.0.as_bytes());
        }
        let id = DatumId(hasher.finalize());
        if let Some(row) = self.datum_row_of(id) {
            return Err(EntityError::DuplicateDatum { id, row });
        }
        let mut retained = Vec::new();
        retained
            .try_reserve(references.len())
            .map_err(|_| EntityError::AllocationRefused {
                resource: "datum-references",
            })?;
        retained.extend_from_slice(references);
        let row = self.datums.len();
        push_checked(
            &mut self.datums,
            Datum {
                id,
                owner,
                declared_name: declared_name.to_string(),
                feature,
                references: retained,
                row,
            },
            "datums",
        )?;
        let position = self
            .datum_index
            .partition_point(|(candidate, _)| *candidate < id);
        self.datum_index
            .try_reserve(1)
            .map_err(|_| EntityError::AllocationRefused {
                resource: "datum-index",
            })?;
        self.datum_index.insert(position, (id, row));
        Ok(id)
    }

    fn datum_row_of(&self, id: DatumId) -> Option<usize> {
        let position = self
            .datum_index
            .partition_point(|(candidate, _)| *candidate < id);
        self.datum_index
            .get(position)
            .and_then(|(candidate, row)| (*candidate == id).then_some(*row))
    }

    /// The datum with this identity.
    #[must_use]
    pub fn datum(&self, id: DatumId) -> Option<&Datum> {
        self.datum_row_of(id).map(|row| &self.datums[row])
    }

    /// Declare a geometric tolerance.
    ///
    /// Structural admission (subject present, capacity) happens here;
    /// dimensional and datum-frame findings are reported by
    /// [`EntityCatalog::validate`] as `Violation`s so an agent gets the whole
    /// repair list at once.
    ///
    /// # Errors
    /// Unknown subject or a capacity refusal.
    pub fn declare_tolerance(&mut self, tolerance: Tolerance) -> Result<usize, EntityError> {
        if self.get(tolerance.subject).is_none() {
            return Err(EntityError::UnknownEntity {
                id: tolerance.subject,
            });
        }
        if self.tolerances.len() >= self.budget.max_tolerances {
            return Err(EntityError::CapacityExceeded {
                resource: "tolerances",
                requested: self.tolerances.len() + 1,
                limit: self.budget.max_tolerances,
            });
        }
        let row = self.tolerances.len();
        let mut tolerance = tolerance;
        tolerance.row = row;
        push_checked(&mut self.tolerances, tolerance, "tolerances")?;
        Ok(row)
    }

    /// Declare the placement of one part occurrence in a scenario frame. The
    /// transform itself lives in the scenario's `FrameTree`; the basis says
    /// whether that transform is the as-designed one or is carried on the
    /// authority of a measured registration.
    ///
    /// # Errors
    /// Unknown or inactive occurrence, a non-part occurrence, an as-built
    /// basis citing the all-zero registration identity, or a capacity refusal.
    pub fn declare_placement(
        &mut self,
        occurrence: EntityId,
        frame: FrameId,
        basis: PlacementBasis,
    ) -> Result<usize, EntityError> {
        let entity = self
            .get(occurrence)
            .ok_or(EntityError::UnknownEntity { id: occurrence })?;
        if entity.status != EntityStatus::Active {
            return Err(EntityError::InactiveIdentityRedeclared {
                id: occurrence,
                status: entity.status,
            });
        }
        if entity.kind() != EntityKind::Part {
            return Err(EntityError::PlacementOccurrenceKind {
                occurrence,
                kind: entity.kind(),
            });
        }
        // The all-zero hash is what an uninitialized or defaulted citation
        // looks like. Admitting it would let an as-built basis — the stronger
        // claim — be declared while naming no artifact at all, which is the
        // one failure this variant exists to prevent. Refuse at declaration
        // rather than leaving it for a downstream resolver that may not run.
        if let PlacementBasis::AsBuilt { registration_ref } = basis
            && registration_ref.as_bytes() == &[0u8; 32]
        {
            return Err(EntityError::PlacementRegistrationUnbound { occurrence });
        }
        if self.placements.len() >= self.budget.max_placements {
            return Err(EntityError::CapacityExceeded {
                resource: "placements",
                requested: self.placements.len() + 1,
                limit: self.budget.max_placements,
            });
        }
        let row = self.placements.len();
        push_checked(
            &mut self.placements,
            Placement {
                occurrence,
                frame,
                basis,
                row,
            },
            "placements",
        )?;
        Ok(row)
    }

    /// Every as-built placement, as `(occurrence, registration_ref)` pairs in
    /// declaration order.
    ///
    /// This is the resolution seam for the product layer: it enumerates the
    /// registration identities a solve would have to resolve and authenticate
    /// before it may claim an as-built geometry term. Nominal placements are
    /// absent because they cite nothing.
    ///
    /// The order is declaration order, so the result is deterministic for a
    /// given catalog and is safe to fold into a downstream identity.
    #[must_use]
    pub fn as_built_registrations(&self) -> Vec<(EntityId, ContentHash)> {
        self.placements
            .iter()
            .filter_map(|placement| {
                placement
                    .basis
                    .registration_ref()
                    .map(|citation| (placement.occurrence, citation))
            })
            .collect()
    }

    /// The placement of one occurrence, when exactly one row names it.
    ///
    /// Returns `None` both when nothing places the occurrence and when more
    /// than one row does — an ambiguous placement is never resolved to a first
    /// match here. [`EntityCatalog::validate`] reports the duplicate as
    /// `placement-duplicate`.
    #[must_use]
    pub fn placement_of(&self, occurrence: EntityId) -> Option<&Placement> {
        let mut found = None;
        for placement in &self.placements {
            if placement.occurrence == occurrence {
                if found.is_some() {
                    return None;
                }
                found = Some(placement);
            }
        }
        found
    }

    /// Catalog-internal findings: inactive parents under active children,
    /// datum hierarchy, and tolerance typing.
    ///
    /// Findings are the crate's structured `Violation { code, what, fix }`.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        let mut out = Vec::new();
        for entity in &self.entities {
            if entity.status != EntityStatus::Active {
                continue;
            }
            let Some(parent) = entity.parent() else {
                continue;
            };
            match self.get(parent).map(Entity::status) {
                None => out.push(Violation {
                    code: "entity-dangling-parent",
                    what: format!(
                        "entity row {} {} declares parent {parent}, which is not in the catalog",
                        entity.row,
                        entity.id
                    ),
                    fix: "declare the parent, or re-import the subtree so a supersession receipt exists"
                        .to_string(),
                }),
                Some(EntityStatus::Active) => {}
                Some(status) => out.push(Violation {
                    code: "entity-inactive-parent",
                    what: format!(
                        "active entity row {} {} hangs under {parent}, which is {status:?}",
                        entity.row, entity.id
                    ),
                    fix: "re-import the child under the parent's current identity, or retire the child"
                        .to_string(),
                }),
            }
        }
        for datum in &self.datums {
            if let Err(fault) = self.resolve(EntityRef::new(datum.owner, KindExpectation::Any)) {
                out.push(Violation {
                    code: "datum-unknown-owner",
                    what: format!(
                        "datum row {} {} is attached to {}: {fault}",
                        datum.row,
                        preview(&datum.declared_name),
                        datum.owner
                    ),
                    fix: fault.fix(),
                });
            }
            for (position, reference) in datum.references.iter().enumerate() {
                if self.datum_row_of(*reference).is_none() {
                    out.push(Violation {
                        code: "datum-dangling-reference",
                        what: format!(
                            "datum row {} reference {position} {reference} is not in the catalog",
                            datum.row
                        ),
                        fix: "declare the referenced datum first".to_string(),
                    });
                }
                if datum.references[..position].contains(reference) {
                    out.push(Violation {
                        code: "datum-repeated-reference",
                        what: format!(
                            "datum row {} repeats reference {reference} at position {position}",
                            datum.row
                        ),
                        fix: "list each establishing datum once".to_string(),
                    });
                }
            }
        }
        for tolerance in &self.tolerances {
            self.check_tolerance(tolerance, &mut out);
        }
        let mut seen: Vec<EntityId> = Vec::new();
        for placement in &self.placements {
            // Duplicates are judged on the RESOLVED occurrence, so placing a
            // predecessor and its successor is one occurrence, not two.
            let resolved = self
                .resolve(EntityRef::new(placement.occurrence, KindExpectation::Any))
                .map(Resolution::current);
            let key = resolved.unwrap_or(placement.occurrence);
            if seen.contains(&key) {
                out.push(Violation {
                    code: "placement-duplicate",
                    what: format!(
                        "placement row {} re-places occurrence {}",
                        placement.row, placement.occurrence
                    ),
                    fix: "declare exactly one placement per occurrence".to_string(),
                });
            } else {
                seen.push(key);
            }
            if let Err(fault) = resolved {
                out.push(Violation {
                    code: "placement-unknown-occurrence",
                    what: format!(
                        "placement row {} places {}: {fault}",
                        placement.row, placement.occurrence
                    ),
                    fix: fault.fix(),
                });
            }
        }
        out
    }

    fn check_tolerance(&self, tolerance: &Tolerance, out: &mut Vec<Violation>) {
        let row = tolerance.row;
        if let Err(fault) = self.resolve(EntityRef::new(tolerance.subject, KindExpectation::Any)) {
            out.push(Violation {
                code: "tolerance-unknown-subject",
                what: format!(
                    "tolerance row {row} controls {}: {fault}",
                    tolerance.subject
                ),
                fix: fault.fix(),
            });
        }
        let expected = tolerance.kind.expected_dims();
        if tolerance.magnitude.dims != expected {
            out.push(Violation {
                code: "tolerance-dims",
                what: format!(
                    "tolerance row {row} ({}) has dimensions {:?}, not {:?}",
                    tolerance.kind, tolerance.magnitude.dims.0, expected.0
                ),
                fix: "express the tolerance magnitude as a length (SI exponents [1,0,0,0,0,0])"
                    .to_string(),
            });
        }
        if !(tolerance.magnitude.value.is_finite() && tolerance.magnitude.value > 0.0) {
            out.push(Violation {
                code: "tolerance-magnitude",
                what: format!(
                    "tolerance row {row} magnitude {} is not a finite positive value",
                    tolerance.magnitude.value
                ),
                fix: "declare a finite positive tolerance zone".to_string(),
            });
        }
        if tolerance.datum_frame.len() > MAX_DATUM_FRAME_LEN {
            out.push(Violation {
                code: "tolerance-datum-frame-arity",
                what: format!(
                    "tolerance row {row} names {} datums; at most {MAX_DATUM_FRAME_LEN} (primary, secondary, tertiary) are admitted",
                    tolerance.datum_frame.len()
                ),
                fix: "reduce the datum reference frame to primary, secondary, and tertiary"
                    .to_string(),
            });
        }
        if tolerance.kind.requires_datum_frame() && tolerance.datum_frame.is_empty() {
            out.push(Violation {
                code: "tolerance-datum-required",
                what: format!(
                    "tolerance row {row} ({}) is an orientation/location control with no datum reference frame",
                    tolerance.kind
                ),
                fix: "name at least a primary datum".to_string(),
            });
        }
        if !tolerance.kind.requires_datum_frame() && !tolerance.datum_frame.is_empty() {
            out.push(Violation {
                code: "tolerance-datum-forbidden",
                what: format!(
                    "tolerance row {row} ({}) is a form control and cannot reference datums",
                    tolerance.kind
                ),
                fix: "remove the datum reference frame from the form control".to_string(),
            });
        }
        for (position, reference) in tolerance.datum_frame.iter().enumerate() {
            if self.datum_row_of(*reference).is_none() {
                out.push(Violation {
                    code: "tolerance-datum-dangling",
                    what: format!(
                        "tolerance row {row} datum {position} {reference} is not in the catalog"
                    ),
                    fix: "declare the datum before referencing it".to_string(),
                });
            }
            if tolerance.datum_frame[..position].contains(reference) {
                out.push(Violation {
                    code: "tolerance-datum-repeated",
                    what: format!("tolerance row {row} repeats datum {reference}"),
                    fix: "name each datum once in the reference frame".to_string(),
                });
            }
        }
        for (field, value) in tolerance.source.nonempty_fields() {
            if value.is_empty() {
                out.push(Violation {
                    code: "tolerance-source-empty",
                    what: format!(
                        "tolerance row {row} declares a {} source with an empty {field}",
                        tolerance.source.label()
                    ),
                    fix: "name the drawing sheet/note, the standard clause, or the assumption's rationale"
                        .to_string(),
                });
            }
        }
    }

    #[cfg(test)]
    fn force_insert_for_collision_test(
        &mut self,
        id: EntityId,
        declaration: EntityDeclaration,
    ) -> usize {
        self.insert_entity(id, declaration)
    }
}

// ---------------------------------------------------------------------------
// Reference sites and binding table
// ---------------------------------------------------------------------------

/// Which side of a contact pair a reference addresses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ContactSide {
    /// `ContactLaw::region_a`.
    A,
    /// `ContactLaw::region_b`.
    B,
}

impl ContactSide {
    /// Stable label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            ContactSide::A => "a",
            ContactSide::B => "b",
        }
    }
}

/// A typed coordinate for the scenario object that holds a reference.
///
/// The four reserved families (`Load`, `MaterialBinding`, `Sensor`,
/// `Requirement`) exist so later scenario objects reuse ONE diagnostic shape.
/// This crate cannot yet enumerate those objects, so their rows are not
/// existence-checked — see the no-claim boundaries in `CONTRACT.md`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReferenceSite {
    /// `Scenario::base_bcs[row]`.
    BaseBoundaryCondition {
        /// Row in `base_bcs`.
        row: usize,
    },
    /// `Scenario::cases[case_row].bcs[bc_row]`.
    CaseBoundaryCondition {
        /// Row in `cases`.
        case_row: usize,
        /// Row in that case's `bcs`.
        bc_row: usize,
    },
    /// One side of `Scenario::contacts[row]`.
    Contact {
        /// Row in `contacts`.
        row: usize,
        /// Which side.
        side: ContactSide,
    },
    /// A load object (reserved).
    Load {
        /// Row in the load collection.
        row: usize,
    },
    /// A material binding (reserved).
    MaterialBinding {
        /// Row in the material-binding collection.
        row: usize,
    },
    /// A sensor (reserved).
    Sensor {
        /// Row in the sensor collection.
        row: usize,
    },
    /// A requirement (reserved).
    Requirement {
        /// Row in the requirement collection.
        row: usize,
    },
}

impl ReferenceSite {
    /// Whether this crate can enumerate the site's owning collection today.
    #[must_use]
    pub const fn is_enumerable(self) -> bool {
        matches!(
            self,
            ReferenceSite::BaseBoundaryCondition { .. }
                | ReferenceSite::CaseBoundaryCondition { .. }
                | ReferenceSite::Contact { .. }
        )
    }

    /// The kind this site requires of whatever it references.
    #[must_use]
    pub const fn required_kind(self) -> KindExpectation {
        match self {
            ReferenceSite::BaseBoundaryCondition { .. }
            | ReferenceSite::CaseBoundaryCondition { .. }
            | ReferenceSite::Contact { .. } => KindExpectation::Boundary,
            ReferenceSite::MaterialBinding { .. } => KindExpectation::Domain,
            ReferenceSite::Load { .. }
            | ReferenceSite::Sensor { .. }
            | ReferenceSite::Requirement { .. } => KindExpectation::Any,
        }
    }
}

impl fmt::Display for ReferenceSite {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReferenceSite::BaseBoundaryCondition { row } => write!(f, "base BC row {row}"),
            ReferenceSite::CaseBoundaryCondition { case_row, bc_row } => {
                write!(f, "case row {case_row} BC row {bc_row}")
            }
            ReferenceSite::Contact { row, side } => {
                write!(f, "contact row {row} side {}", side.label())
            }
            ReferenceSite::Load { row } => write!(f, "load row {row}"),
            ReferenceSite::MaterialBinding { row } => write!(f, "material binding row {row}"),
            ReferenceSite::Sensor { row } => write!(f, "sensor row {row}"),
            ReferenceSite::Requirement { row } => write!(f, "requirement row {row}"),
        }
    }
}

/// One site-to-identity binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Binding {
    site: ReferenceSite,
    reference: EntityRef,
}

impl Binding {
    /// The referring site.
    #[must_use]
    pub const fn site(self) -> ReferenceSite {
        self.site
    }

    /// The typed reference.
    #[must_use]
    pub const fn reference(self) -> EntityRef {
        self.reference
    }
}

/// One row of a rendered binding table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BindingRow {
    /// The referring site.
    pub site: ReferenceSite,
    /// The identity as authored.
    pub requested: EntityId,
    /// The current identity, when the reference resolves.
    pub current: Option<EntityId>,
    /// The current display name, when the reference resolves.
    pub display_name: Option<String>,
    /// Supersession hops followed.
    pub hops: usize,
    /// Weakest evidence tier along the chain, when the reference resolves.
    pub tier: Option<EvidenceTier>,
    /// Why the reference did not resolve, when it did not.
    pub fault: Option<ResolutionFault>,
}

/// The scenario's site-to-identity bindings.
///
/// This is the authoritative reference layer: scenario objects name entities
/// by identity here, while the legacy string fields on `BoundaryCondition` and
/// `ContactLaw` remain the v2 IR wire form.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BindingTable {
    bindings: Vec<Binding>,
}

impl BindingTable {
    /// An empty table.
    #[must_use]
    pub fn new() -> Self {
        BindingTable {
            bindings: Vec::new(),
        }
    }

    /// Every binding in insertion order.
    #[must_use]
    pub fn bindings(&self) -> &[Binding] {
        &self.bindings
    }

    /// Bind a site to an identity.
    ///
    /// Duplicate sites are admitted here and reported as
    /// `entity-duplicate-binding` by [`validate_bindings`], so an agent gets
    /// the whole repair list at once instead of one refusal at a time.
    ///
    /// # Errors
    /// [`EntityError::CapacityExceeded`] or [`EntityError::AllocationRefused`].
    pub fn bind(
        &mut self,
        site: ReferenceSite,
        reference: EntityRef,
        budget: EntityBudget,
    ) -> Result<(), EntityError> {
        if self.bindings.len() >= budget.max_bindings {
            return Err(EntityError::CapacityExceeded {
                resource: "bindings",
                requested: self.bindings.len() + 1,
                limit: budget.max_bindings,
            });
        }
        push_checked(&mut self.bindings, Binding { site, reference }, "bindings")
    }

    /// Bind a site by DECLARED NAME, refusing ambiguity.
    ///
    /// This is the bridge for callers that still hold strings. A name that
    /// matches two entities is a typed refusal naming both candidates — never
    /// a first match.
    ///
    /// # Errors
    /// [`EntityError::UnknownDeclaredName`],
    /// [`EntityError::AmbiguousDeclaredName`], or a capacity refusal.
    pub fn bind_by_declared_name(
        &mut self,
        site: ReferenceSite,
        name: &str,
        expect: KindExpectation,
        catalog: &EntityCatalog,
    ) -> Result<EntityId, EntityError> {
        let id = catalog.resolve_declared_name(name)?;
        self.bind(site, EntityRef::new(id, expect), catalog.budget())?;
        Ok(id)
    }

    /// The single binding for one site, or `None` when unbound or duplicated.
    #[must_use]
    pub fn binding_for(&self, site: ReferenceSite) -> Option<Binding> {
        let mut found = None;
        for binding in &self.bindings {
            if binding.site == site {
                if found.is_some() {
                    return None;
                }
                found = Some(*binding);
            }
        }
        found
    }

    /// Render the resolved binding table in insertion order.
    #[must_use]
    pub fn report(&self, catalog: &EntityCatalog) -> Vec<BindingRow> {
        let mut rows = Vec::new();
        for binding in &self.bindings {
            let row = match catalog.resolve(binding.reference) {
                Ok(resolution) => BindingRow {
                    site: binding.site,
                    requested: resolution.requested(),
                    current: Some(resolution.current()),
                    display_name: catalog
                        .get(resolution.current())
                        .map(|entity| entity.display_name().to_string()),
                    hops: resolution.hops(),
                    tier: Some(resolution.tier()),
                    fault: None,
                },
                Err(fault) => BindingRow {
                    site: binding.site,
                    requested: binding.reference.target(),
                    current: None,
                    display_name: None,
                    hops: 0,
                    tier: None,
                    fault: Some(fault),
                },
            };
            rows.push(row);
        }
        rows
    }
}

/// Every enumerable string-bearing reference site in a scenario, in a
/// deterministic order: base BCs, then case BCs, then contact sides.
#[must_use]
pub fn scenario_reference_sites(scenario: &Scenario) -> Vec<(ReferenceSite, &str)> {
    let mut sites = Vec::new();
    for (row, bc) in scenario.base_bcs.iter().enumerate() {
        sites.push((
            ReferenceSite::BaseBoundaryCondition { row },
            bc.region.as_str(),
        ));
    }
    for (case_row, case) in scenario.cases.iter().enumerate() {
        for (bc_row, bc) in case.bcs.iter().enumerate() {
            sites.push((
                ReferenceSite::CaseBoundaryCondition { case_row, bc_row },
                bc.region.as_str(),
            ));
        }
    }
    for (row, contact) in scenario.contacts.iter().enumerate() {
        sites.push((
            ReferenceSite::Contact {
                row,
                side: ContactSide::A,
            },
            contact.region_a.as_str(),
        ));
        sites.push((
            ReferenceSite::Contact {
                row,
                side: ContactSide::B,
            },
            contact.region_b.as_str(),
        ));
    }
    sites
}

/// Resolve EVERY reference in the binding table and report the failures.
///
/// Findings use the crate's `Violation { code, what, fix }` shape. Codes:
/// `entity-dangling-reference`, `entity-retired-reference`,
/// `entity-supersession-depth`, `entity-kind-mismatch`,
/// `entity-site-kind-mismatch`, `entity-duplicate-binding`,
/// `entity-orphan-site`, `entity-unbound-site`, `entity-ambiguous-name`,
/// `entity-legacy-name-drift`, plus the catalog and placement codes from
/// [`EntityCatalog::validate`].
#[must_use]
pub fn validate_bindings(
    scenario: &Scenario,
    catalog: &EntityCatalog,
    bindings: &BindingTable,
) -> Vec<Violation> {
    let mut out = catalog.validate();
    let sites = scenario_reference_sites(scenario);
    check_placement_frames(scenario, catalog, &mut out);
    check_bound_sites(catalog, bindings, &sites, &mut out);
    check_unbound_sites(catalog, bindings, &sites, &mut out);
    out
}

fn check_placement_frames(scenario: &Scenario, catalog: &EntityCatalog, out: &mut Vec<Violation>) {
    for placement in catalog.placements() {
        let frame = placement.frame();
        if frame == WORLD {
            continue;
        }
        let matches = scenario
            .frames
            .frames
            .iter()
            .filter(|candidate| candidate.id == frame)
            .count();
        if matches == 0 {
            out.push(Violation {
                code: "placement-unknown-frame",
                what: format!(
                    "placement row {} names frame {}, which the scenario does not declare",
                    placement.row(),
                    frame.0
                ),
                fix: "declare the placement frame in the scenario FrameTree, or place the occurrence in the world frame"
                    .to_string(),
            });
        } else if matches > 1 {
            out.push(Violation {
                code: "placement-ambiguous-frame",
                what: format!(
                    "placement row {} names frame {}, which the scenario declares {matches} times",
                    placement.row(),
                    frame.0
                ),
                fix: "give each scenario frame a unique id before placing occurrences in it"
                    .to_string(),
            });
        }
    }
}

fn check_bound_sites(
    catalog: &EntityCatalog,
    bindings: &BindingTable,
    sites: &[(ReferenceSite, &str)],
    out: &mut Vec<Violation>,
) {
    for (position, binding) in bindings.bindings().iter().enumerate() {
        let site = binding.site();
        if bindings.bindings()[..position]
            .iter()
            .any(|earlier| earlier.site() == site)
        {
            out.push(Violation {
                code: "entity-duplicate-binding",
                what: format!("{site} is bound more than once (binding row {position})"),
                fix: "bind each reference site exactly once".to_string(),
            });
        }
        if site.is_enumerable() && !sites.iter().any(|(candidate, _)| *candidate == site) {
            out.push(Violation {
                code: "entity-orphan-site",
                what: format!(
                    "binding row {position} names {site}, which the scenario does not have"
                ),
                fix: "remove the binding, or restore the scenario object it names".to_string(),
            });
            continue;
        }
        match catalog.resolve(binding.reference()) {
            Err(fault) => out.push(Violation {
                code: fault.code(),
                what: format!("{site}: {fault}"),
                fix: fault.fix(),
            }),
            Ok(resolution) => {
                let actual = resolution.current().kind();
                if !site.required_kind().admits(actual) {
                    out.push(Violation {
                        code: "entity-site-kind-mismatch",
                        what: format!(
                            "{site} resolves to {} (a {actual}) where {} is required",
                            resolution.current(),
                            site.required_kind()
                        ),
                        fix: format!("bind {site} to a {} entity", site.required_kind()),
                    });
                }
                if let Some((_, declared)) = sites.iter().find(|(candidate, _)| *candidate == site)
                    && let Some(entity) = catalog.get(resolution.current())
                    && entity.is_legacy()
                    && entity.declared_name() != *declared
                {
                    out.push(Violation {
                        code: "entity-legacy-name-drift",
                        what: format!(
                            "{site} carries the legacy string {} but is bound to legacy entity {} declared as {}",
                            preview(declared),
                            resolution.current(),
                            preview(entity.declared_name())
                        ),
                        fix: "re-run the legacy migration for this scenario, or replace the legacy entity with an imported identity"
                            .to_string(),
                    });
                }
            }
        }
    }
}

fn check_unbound_sites(
    catalog: &EntityCatalog,
    bindings: &BindingTable,
    sites: &[(ReferenceSite, &str)],
    out: &mut Vec<Violation>,
) {
    for (site, declared) in sites {
        if bindings
            .bindings()
            .iter()
            .any(|binding| binding.site() == *site)
        {
            continue;
        }
        match catalog.lookup_declared_name(declared) {
            NameLookup::Ambiguous {
                first,
                second,
                total,
            } => out.push(Violation {
                code: "entity-ambiguous-name",
                what: format!(
                    "{site} is unbound and its string {} names {total} entities including {first} and {second}",
                    preview(declared)
                ),
                fix: "bind this site to the intended identity; the string cannot choose between them"
                    .to_string(),
            }),
            NameLookup::Unique(id) => out.push(Violation {
                code: "entity-unbound-site",
                what: format!(
                    "{site} references geometry by the string {} only",
                    preview(declared)
                ),
                fix: format!("bind {site} to {}", id.token()),
            }),
            NameLookup::Missing => out.push(Violation {
                code: "entity-unbound-site",
                what: format!(
                    "{site} references geometry by the string {} only, and no entity carries that declared name",
                    preview(declared)
                ),
                fix: "declare the entity (or migrate the scenario with migrate_legacy_scenario), then bind this site"
                    .to_string(),
            }),
        }
    }
}

/// What a mechanical legacy migration produced.
#[derive(Debug, Clone, PartialEq)]
pub struct LegacyMigration {
    root: EntityId,
    part: EntityId,
    surfaces: Vec<(String, EntityId)>,
    bindings: BindingTable,
    first_receipt: u64,
    receipts: u64,
}

impl LegacyMigration {
    /// The assembly the migration declared (or reused).
    #[must_use]
    pub const fn root(&self) -> EntityId {
        self.root
    }

    /// The synthetic part every migrated surface hangs under.
    ///
    /// A bare string carries no part structure; inventing one would be a
    /// stronger claim than the input supports.
    #[must_use]
    pub const fn part(&self) -> EntityId {
        self.part
    }

    /// Distinct migrated strings and their identities, in first-use order.
    #[must_use]
    pub fn surfaces(&self) -> &[(String, EntityId)] {
        &self.surfaces
    }

    /// The complete binding table for the migrated scenario.
    #[must_use]
    pub const fn bindings(&self) -> &BindingTable {
        &self.bindings
    }

    /// Sequence of the first receipt this migration appended.
    #[must_use]
    pub const fn first_receipt(&self) -> u64 {
        self.first_receipt
    }

    /// How many receipts this migration appended (zero on a repeat run).
    #[must_use]
    pub const fn receipts(&self) -> u64 {
        self.receipts
    }
}

/// Mechanically migrate a string-referencing scenario into the entity model.
///
/// Every distinct region string becomes a `Surface` whose DECLARED NAME is
/// that exact string, carrying the legacy marker, under one synthetic part
/// inside `root_name`'s assembly. No geometry fingerprint is invented, and no
/// part structure is inferred: a string does not carry either.
///
/// The migration is idempotent — the legacy marker is metadata, not identity,
/// so re-running it re-derives the same identities and appends no receipts.
///
/// # Errors
/// Any declaration refusal, including an empty region string.
pub fn migrate_legacy_scenario(
    scenario: &Scenario,
    catalog: &mut EntityCatalog,
    root_name: &str,
) -> Result<LegacyMigration, EntityError> {
    let first_receipt = catalog.receipts().len() as u64;
    let root =
        catalog.declare_or_existing(EntityDeclaration::assembly(root_name).with_legacy_marker())?;
    let part = catalog
        .declare_or_existing(EntityDeclaration::part(root, "legacy").with_legacy_marker())?;
    let mut surfaces: Vec<(String, EntityId)> = Vec::new();
    let mut bindings = BindingTable::new();
    for (site, declared) in scenario_reference_sites(scenario) {
        let id = if let Some((_, id)) = surfaces.iter().find(|(name, _)| name == declared) {
            *id
        } else {
            let id = catalog.declare_or_existing(
                EntityDeclaration::surface(part, declared).with_legacy_marker(),
            )?;
            push_checked(&mut surfaces, (declared.to_string(), id), "legacy-surfaces")?;
            id
        };
        bindings.bind(
            site,
            EntityRef::new(id, KindExpectation::Boundary),
            catalog.budget(),
        )?;
    }
    let receipts = catalog.receipts().len() as u64 - first_receipt;
    Ok(LegacyMigration {
        root,
        part,
        surfaces,
        bindings,
        first_receipt,
        receipts,
    })
}

#[cfg(test)]
mod identity_tests {
    use super::*;

    fn stack() -> (EntityCatalog, EntityId, EntityId) {
        let mut catalog = EntityCatalog::new();
        let assembly = catalog
            .declare(EntityDeclaration::assembly("stack"))
            .expect("assembly");
        let part = catalog
            .declare(EntityDeclaration::part(assembly, "plate"))
            .expect("part");
        (catalog, assembly, part)
    }

    #[test]
    fn identity_ignores_display_name_but_tracks_declared_name() {
        let (_, assembly, _) = stack();
        let base = EntityDeclaration::part(assembly, "cold-plate");
        let renamed = EntityDeclaration::part(assembly, "cold-plate")
            .with_display_name("Cold Plate (rev B), anodized");
        assert_eq!(base.identity(), renamed.identity());
        let other = EntityDeclaration::part(assembly, "cold-plate-2");
        assert_ne!(base.identity(), other.identity());
    }

    #[test]
    fn identity_preimage_is_length_prefixed_and_injective() {
        let (mut catalog, assembly, _) = stack();
        // A naive concatenation of parent name and child name would map
        // ("ab", "c") and ("a", "bc") onto one preimage.
        let left_parent = catalog
            .declare(EntityDeclaration::sub_assembly(assembly, "ab"))
            .expect("left");
        let right_parent = catalog
            .declare(EntityDeclaration::sub_assembly(assembly, "a"))
            .expect("right");
        let left = EntityDeclaration::part(left_parent, "c").identity();
        let right = EntityDeclaration::part(right_parent, "bc").identity();
        assert_ne!(left, right);

        // Same declared name, different parent path.
        let same_name_left = EntityDeclaration::part(left_parent, "x").identity();
        let same_name_right = EntityDeclaration::part(right_parent, "x").identity();
        assert_ne!(same_name_left, same_name_right);

        // Same name and parent, different kind.
        let part = EntityDeclaration::part(left_parent, "y").identity();
        let sub = EntityDeclaration::sub_assembly(left_parent, "y").identity();
        assert_ne!(part.digest(), sub.digest());
        assert_eq!(part.kind(), EntityKind::Part);
        assert_eq!(sub.kind(), EntityKind::Assembly);
    }

    #[test]
    fn fingerprints_participate_in_identity() {
        let (_, assembly, _) = stack();
        let plain = EntityDeclaration::part(assembly, "die");
        let hashed = EntityDeclaration::part(assembly, "die")
            .with_fingerprint(GeometryFingerprint::of_bytes(b"solid-die-v1"));
        let rehashed = EntityDeclaration::part(assembly, "die")
            .with_fingerprint(GeometryFingerprint::of_bytes(b"solid-die-v2"));
        assert_ne!(plain.identity(), hashed.identity());
        assert_ne!(hashed.identity(), rehashed.identity());
        assert_eq!(
            hashed.identity(),
            EntityDeclaration::part(assembly, "die")
                .with_fingerprint(GeometryFingerprint::of_bytes(b"solid-die-v1"))
                .identity()
        );
        // An empty byte string is a fingerprint, not the absence of one.
        assert_ne!(
            plain.identity(),
            EntityDeclaration::part(assembly, "die")
                .with_fingerprint(GeometryFingerprint::of_bytes(b""))
                .identity()
        );
    }

    #[test]
    fn interface_ordering_is_identity_bearing() {
        let (mut catalog, assembly, part) = stack();
        let top = catalog
            .declare(EntityDeclaration::surface(part, "top"))
            .expect("top");
        let bottom = catalog
            .declare(EntityDeclaration::surface(part, "bottom"))
            .expect("bottom");
        let forward =
            EntityDeclaration::interface(assembly, "tim", InterfacePair::ordered(top, bottom));
        let backward =
            EntityDeclaration::interface(assembly, "tim", InterfacePair::ordered(bottom, top));
        assert_ne!(forward.identity(), backward.identity());
        assert_eq!(
            forward.pair().and_then(InterfacePair::applied_side),
            Some(top)
        );

        let one = EntityDeclaration::interface(
            assembly,
            "contact",
            InterfacePair::unordered(top, bottom),
        );
        let other = EntityDeclaration::interface(
            assembly,
            "contact",
            InterfacePair::unordered(bottom, top),
        );
        assert_eq!(one.identity(), other.identity());
        assert_eq!(
            one.pair().and_then(InterfacePair::applied_side),
            None,
            "an unordered pair refuses to name an applied side"
        );
    }

    #[test]
    fn interface_sides_are_checked() {
        let (mut catalog, assembly, part) = stack();
        let top = catalog
            .declare(EntityDeclaration::surface(part, "top"))
            .expect("top");
        let orphan_part = EntityDeclaration::part(assembly, "absent").identity();
        let orphan_surface = EntityDeclaration::surface(orphan_part, "ghost").identity();
        assert_eq!(
            catalog.declare(EntityDeclaration::interface(
                assembly,
                "self",
                InterfacePair::ordered(top, top)
            )),
            Err(EntityError::SelfInterface { side: top })
        );
        assert_eq!(
            catalog.declare(EntityDeclaration::interface(
                assembly,
                "ghost",
                InterfacePair::ordered(top, orphan_surface)
            )),
            Err(EntityError::UnknownInterfaceSide {
                side: orphan_surface
            })
        );
        assert_eq!(
            catalog.declare(EntityDeclaration::interface(
                assembly,
                "wrong-kind",
                InterfacePair::ordered(top, part)
            )),
            Err(EntityError::InterfaceSideKind {
                side: part,
                kind: EntityKind::Part
            })
        );
    }

    #[test]
    fn containment_rules_are_enforced() {
        let (mut catalog, assembly, part) = stack();
        assert_eq!(
            catalog.declare(EntityDeclaration::part(part, "nested")),
            Err(EntityError::ParentKind {
                child: EntityKind::Part,
                parent: Some(EntityKind::Part)
            })
        );
        assert_eq!(
            catalog.declare(EntityDeclaration::surface(assembly, "loose")),
            Err(EntityError::ParentKind {
                child: EntityKind::Surface,
                parent: Some(EntityKind::Assembly)
            })
        );
        let unknown = EntityDeclaration::assembly("elsewhere").identity();
        assert_eq!(
            catalog.declare(EntityDeclaration::part(unknown, "floating")),
            Err(EntityError::UnknownParent {
                kind: EntityKind::Part,
                parent: unknown
            })
        );
    }

    #[test]
    fn duplicate_declaration_and_digest_collision_are_distinct_refusals() {
        let (mut catalog, _, part) = stack();
        let surface = catalog
            .declare(EntityDeclaration::surface(part, "inlet"))
            .expect("surface");
        let row = catalog.get(surface).expect("row").row();
        assert_eq!(
            catalog.declare(EntityDeclaration::surface(part, "inlet")),
            Err(EntityError::DuplicateDeclaration { id: surface, row })
        );
        // The duplicate path must not be reached by a DIFFERENT preimage that
        // lands on the same digest. Force that state and prove it fails closed
        // rather than aliasing two declarations onto one identity.
        let victim = EntityDeclaration::surface(part, "outlet");
        let forged = EntityDeclaration::surface(part, "not-outlet");
        catalog.force_insert_for_collision_test(victim.identity(), forged);
        let error = catalog
            .declare(victim)
            .expect_err("a colliding declaration must be refused");
        assert_eq!(error.code(), "entity-identity-collision");
        assert!(matches!(error, EntityError::IdentityCollision { .. }));
    }

    #[test]
    fn rename_moves_the_display_name_only() {
        let (mut catalog, _, part) = stack();
        let before_root = catalog.receipt_root();
        let sequence = catalog.rename(part, "Cold Plate (rev B)").expect("rename");
        let entity = catalog.get(part).expect("entity");
        assert_eq!(entity.id(), part);
        assert_eq!(entity.declared_name(), "plate");
        assert_eq!(entity.display_name(), "Cold Plate (rev B)");
        let receipt = &catalog.receipts()[sequence as usize];
        assert_eq!(receipt.event(), RebindEvent::Rename);
        assert_eq!(receipt.basis(), MatchBasis::Unchanged);
        assert_eq!(receipt.display_before(), Some("plate"));
        assert_eq!(receipt.display_after(), Some("Cold Plate (rev B)"));
        assert_ne!(before_root, catalog.receipt_root());
        assert!(catalog.verify_receipts());
    }

    #[test]
    fn receipt_chain_verifies_and_detects_tampering() {
        let (mut catalog, _, part) = stack();
        catalog.rename(part, "Plate B").expect("rename");
        assert!(catalog.verify_receipts());
        assert!(catalog.receipts().iter().all(IdentityReceipt::verifies));

        let mut field_tampered = catalog.clone();
        field_tampered.receipts[1].display_after = Some("Plate C".to_string());
        assert!(!field_tampered.receipts[1].verifies());
        assert!(!field_tampered.verify_receipts());

        let mut reordered = catalog.clone();
        reordered.receipts.swap(0, 1);
        assert!(!reordered.verify_receipts());

        let mut truncated = catalog.clone();
        truncated.receipts.pop();
        assert!(
            !truncated.verify_receipts(),
            "the retained root still pins the dropped tail"
        );
    }

    #[test]
    fn resolution_reports_the_weakest_link_in_the_chain() {
        let mut catalog = EntityCatalog::new();
        let assembly = catalog
            .declare(EntityDeclaration::assembly("stack"))
            .expect("assembly");
        let part = catalog
            .declare(EntityDeclaration::part(assembly, "plate"))
            .expect("part");
        let first = catalog
            .declare(
                EntityDeclaration::surface(part, "inlet")
                    .with_fingerprint(GeometryFingerprint::of_bytes(b"patch-1")),
            )
            .expect("first");

        let second_declaration = EntityDeclaration::surface(part, "inlet-renamed")
            .with_fingerprint(GeometryFingerprint::of_bytes(b"patch-2"));
        let second = second_declaration.identity();
        catalog
            .apply_import(&ImportRevision {
                label: "rev-b".to_string(),
                event: RebindEvent::Import,
                scope: ImportScope::Partial,
                entities: vec![ImportedEntity {
                    declaration: second_declaration,
                    correspondence: Correspondence::Asserted(first),
                }],
            })
            .expect("asserted import");

        let third_declaration = EntityDeclaration::surface(part, "inlet-final")
            .with_fingerprint(GeometryFingerprint::of_bytes(b"patch-2"));
        let third = third_declaration.identity();
        catalog
            .apply_import(&ImportRevision {
                label: "rev-c".to_string(),
                event: RebindEvent::Import,
                scope: ImportScope::Partial,
                entities: vec![ImportedEntity {
                    declaration: third_declaration,
                    correspondence: Correspondence::Auto,
                }],
            })
            .expect("auto import");

        let resolution = catalog
            .resolve(EntityRef::new(first, KindExpectation::Boundary))
            .expect("resolves");
        assert_eq!(resolution.current(), third);
        assert_eq!(resolution.hops(), 2);
        assert_eq!(
            resolution.tier(),
            EvidenceTier::Asserted,
            "a content-matched hop must not launder an asserted hop"
        );

        let direct = catalog
            .resolve(EntityRef::new(second, KindExpectation::Boundary))
            .expect("resolves");
        assert_eq!(direct.tier(), EvidenceTier::ContentMatched);
        assert_eq!(direct.hops(), 1);
        assert!(catalog.verify_receipts());
    }

    #[test]
    fn supersession_depth_is_budgeted() {
        let budget = EntityBudget {
            max_supersession_hops: 1,
            ..DEFAULT_ENTITY_BUDGET
        };
        let mut catalog = EntityCatalog::with_budget(budget);
        let assembly = catalog
            .declare(EntityDeclaration::assembly("stack"))
            .expect("assembly");
        let part = catalog
            .declare(EntityDeclaration::part(assembly, "plate"))
            .expect("part");
        let first = catalog
            .declare(EntityDeclaration::surface(part, "s0"))
            .expect("s0");
        let mut previous = first;
        for step in 1..=2u32 {
            let declaration = EntityDeclaration::surface(part, &format!("s{step}"));
            let id = declaration.identity();
            catalog
                .apply_import(&ImportRevision {
                    label: format!("rev-{step}"),
                    event: RebindEvent::Import,
                    scope: ImportScope::Partial,
                    entities: vec![ImportedEntity {
                        declaration,
                        correspondence: Correspondence::Asserted(previous),
                    }],
                })
                .expect("import");
            previous = id;
        }
        assert_eq!(
            catalog.resolve(EntityRef::new(first, KindExpectation::Boundary)),
            Err(ResolutionFault::DepthExceeded {
                id: first,
                limit: 1
            })
        );
    }

    #[test]
    fn ambiguous_automatic_matching_is_refused_not_guessed() {
        let mut catalog = EntityCatalog::new();
        let assembly = catalog
            .declare(EntityDeclaration::assembly("stack"))
            .expect("assembly");
        let fastener = GeometryFingerprint::of_bytes(b"m3-screw");
        let left = catalog
            .declare(EntityDeclaration::part(assembly, "screw-left").with_fingerprint(fastener))
            .expect("left");
        let right = catalog
            .declare(EntityDeclaration::part(assembly, "screw-right").with_fingerprint(fastener))
            .expect("right");
        let incoming = EntityDeclaration::part(assembly, "screw-a").with_fingerprint(fastener);
        let error = catalog
            .apply_import(&ImportRevision {
                label: "rev-b".to_string(),
                event: RebindEvent::Import,
                scope: ImportScope::Partial,
                entities: vec![ImportedEntity {
                    declaration: incoming.clone(),
                    correspondence: Correspondence::Auto,
                }],
            })
            .expect_err("two identical fingerprints cannot be matched automatically");
        assert_eq!(
            error,
            EntityError::AmbiguousImportMatch {
                incoming: incoming.identity(),
                first: left,
                second: right,
                total: 2,
                tier: EvidenceTier::ContentMatched,
            }
        );
        assert!(
            catalog.get(incoming.identity()).is_none(),
            "a refused revision must not mutate the catalog"
        );
        assert!(catalog.verify_receipts());
    }

    #[test]
    fn a_refused_revision_leaves_the_catalog_untouched() {
        let (mut catalog, assembly, part) = stack();
        catalog
            .declare(EntityDeclaration::surface(part, "top"))
            .expect("top");
        let before = catalog.clone();
        let ghost_parent = EntityDeclaration::assembly("ghost").identity();
        let error = catalog
            .apply_import(&ImportRevision {
                label: "rev-b".to_string(),
                event: RebindEvent::Import,
                scope: ImportScope::Partial,
                entities: vec![
                    ImportedEntity {
                        declaration: EntityDeclaration::part(assembly, "new-part"),
                        correspondence: Correspondence::New,
                    },
                    ImportedEntity {
                        declaration: EntityDeclaration::part(ghost_parent, "orphan"),
                        correspondence: Correspondence::New,
                    },
                ],
            })
            .expect_err("the second entity must refuse");
        assert_eq!(error.code(), "entity-unknown-parent");
        assert_eq!(catalog, before);
    }

    #[test]
    fn import_events_must_be_revision_events() {
        let (mut catalog, assembly, _) = stack();
        assert_eq!(
            catalog.apply_import(&ImportRevision {
                label: "rev".to_string(),
                event: RebindEvent::Rename,
                scope: ImportScope::Partial,
                entities: vec![ImportedEntity {
                    declaration: EntityDeclaration::part(assembly, "x"),
                    correspondence: Correspondence::New,
                }],
            }),
            Err(EntityError::NonRevisionEvent {
                event: RebindEvent::Rename
            })
        );
    }

    #[test]
    fn datum_and_placement_owners_are_kind_checked() {
        let (mut catalog, assembly, part) = stack();
        let top = catalog
            .declare(EntityDeclaration::surface(part, "top"))
            .expect("top");
        let bottom = catalog
            .declare(EntityDeclaration::surface(part, "bottom"))
            .expect("bottom");
        let interface = catalog
            .declare(EntityDeclaration::interface(
                assembly,
                "seam",
                InterfacePair::unordered(top, bottom),
            ))
            .expect("interface");
        assert_eq!(
            catalog.declare_datum(interface, "A", DatumFeature::Plane, &[]),
            Err(EntityError::DatumOwnerKind {
                owner: interface,
                kind: EntityKind::Interface
            })
        );
        let datum = catalog
            .declare_datum(part, "A", DatumFeature::Plane, &[])
            .expect("datum");
        assert_eq!(
            catalog.declare_datum(part, "A", DatumFeature::Plane, &[]),
            Err(EntityError::DuplicateDatum { id: datum, row: 0 })
        );
        assert_eq!(
            catalog.declare_placement(top, WORLD, PlacementBasis::Nominal),
            Err(EntityError::PlacementOccurrenceKind {
                occurrence: top,
                kind: EntityKind::Surface
            })
        );
        catalog
            .declare_placement(part, WORLD, PlacementBasis::Nominal)
            .expect("placement");
        assert!(catalog.validate().is_empty());
    }

    #[test]
    fn an_import_cannot_alias_onto_a_colliding_identity() {
        let (mut catalog, _, part) = stack();
        let victim = EntityDeclaration::surface(part, "outlet");
        let forged = EntityDeclaration::surface(part, "not-outlet");
        catalog.force_insert_for_collision_test(victim.identity(), forged);
        let before = catalog.clone();
        let error = catalog
            .apply_import(&ImportRevision {
                label: "rev-b".to_string(),
                event: RebindEvent::Import,
                scope: ImportScope::Partial,
                entities: vec![ImportedEntity {
                    declaration: victim,
                    correspondence: Correspondence::New,
                }],
            })
            .expect_err("a colliding revision entity must be refused");
        assert_eq!(error.code(), "entity-identity-collision");
        assert_eq!(catalog, before);
    }

    #[test]
    fn declared_names_are_not_identities() {
        let mut catalog = EntityCatalog::new();
        let assembly = catalog
            .declare(EntityDeclaration::assembly("stack"))
            .expect("assembly");
        let left = catalog
            .declare(EntityDeclaration::part(assembly, "left"))
            .expect("left");
        let right = catalog
            .declare(EntityDeclaration::part(assembly, "right"))
            .expect("right");
        let left_inlet = catalog
            .declare(EntityDeclaration::surface(left, "inlet"))
            .expect("left inlet");
        let right_inlet = catalog
            .declare(EntityDeclaration::surface(right, "inlet"))
            .expect("right inlet");
        assert_ne!(left_inlet, right_inlet);
        assert_eq!(
            catalog.lookup_declared_name("inlet"),
            NameLookup::Ambiguous {
                first: left_inlet,
                second: right_inlet,
                total: 2
            }
        );
        let error = catalog
            .resolve_declared_name("inlet")
            .expect_err("a string cannot choose");
        assert_eq!(error.code(), "entity-ambiguous-declared-name");
        assert_eq!(catalog.lookup_declared_name("outlet"), NameLookup::Missing);
    }
}
