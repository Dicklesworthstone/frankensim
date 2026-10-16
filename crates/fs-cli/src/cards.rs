//! Admitted material and interface card-pack inputs for the `solve` verb
//! (bead frankensim-hp7tb, half (b)).
//!
//! The solve pipeline's `material-resolve` stage needs a
//! [`fs_project::CardLibrary`], which holds `fs_matdb` cards keyed by content
//! hash. Cards travel on the wire inside the normalized pack envelopes
//! `FSMCDPK\0` ([`NormalizedMaterialCardPack`]) and `FSINTPK\0`
//! ([`NormalizedInterfacePack`]). This module owns the boundary between
//! caller-supplied pack bytes and that library.
//!
//! # What admission establishes
//!
//! * Every pack decodes through its own canonical envelope. `from_bytes`
//!   re-encodes and compares, so semantically valid but non-canonically
//!   spelled bytes refuse rather than silently normalizing.
//! * The set is canonical by content: packs are ordered by pack root, so
//!   flag order cannot change [`CardPackSet::root`] or anything derived from
//!   it. Byte-identical duplicates are idempotent.
//! * Two distinct packs that reconstruct the *same* card refuse as a
//!   conflict. [`fs_project::CardLibrary`] is keyed by card hash, so
//!   admitting both would silently drop one pack's provenance; there is no
//!   last-one-wins path.
//!
//! # What admission does not establish
//!
//! Admission is a decoding and set-canonicality gate. It does not
//! authenticate a pack's producer, does not check that the packs cover the
//! project's declared bindings (that is [`fs_project::resolve_bindings`]),
//! and confers no scientific authority on the claims a pack carries. A pack
//! that decodes is admissible input, not validated data.
//!
//! # Source labels are diagnostics, never identity
//!
//! Caller paths are retained only for refusal text. They never enter
//! [`CardPackSet::root`] or any receipt, because the same content supplied
//! from a different path must produce the same run and the same stage
//! receipt — resume re-attests packs from the ledger and has no paths at all.

use fs_blake3::{ContentHash, hash_bytes, hash_domain};
use fs_matdb::{NormalizedInterfacePack, NormalizedMaterialCardPack, PackError};
use fs_project::CardLibrary;

/// Domain separating the canonical card-pack-set root from every other hash.
pub const CARD_PACK_SET_DOMAIN: &str = "org.frankensim.fs-cli.solve-card-pack-set.v1";

/// Per-file read cap for one card pack.
///
/// Decoding one pack is a single indivisible step inside `fs-matdb`, so this
/// cap is also the coarsest cancellation latency the admission loop can
/// offer: the driver observes cancellation at pack boundaries, not inside a
/// decode. Four mebibytes is far above any realistic claim table and matches
/// the retained-receipt read cap the resume path already uses.
pub const MAX_CARD_PACK_BYTES: u64 = 4 * 1024 * 1024;

/// Maximum packs one invocation may supply, across both kinds.
///
/// The stage operation links one ledger edge per pack plus its predecessor,
/// receipt, and checkpoint; the run-input operation additionally carries the
/// project source and every retained import edge. This ceiling keeps both
/// operations well inside the ledger's bounded edge scan.
pub const MAX_CARD_PACKS: usize = 128;

/// Longest caller-supplied source label retained for diagnostics.
pub const MAX_CARD_PACK_SOURCE_BYTES: usize = 4096;

/// Which normalized pack envelope a caller supplied.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CardPackKind {
    /// `FSMCDPK\0` — one material card and its complete claim pack.
    Material,
    /// `FSINTPK\0` — one interface-system card and its complete claim pack.
    Interface,
}

impl CardPackKind {
    /// The grammar flag that supplies this kind.
    #[must_use]
    pub const fn flag(self) -> &'static str {
        match self {
            CardPackKind::Material => "--materials",
            CardPackKind::Interface => "--interfaces",
        }
    }

    /// Stable lowercase label used in diagnostics and receipts.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            CardPackKind::Material => "material",
            CardPackKind::Interface => "interface",
        }
    }

    /// Ledger artifact kind retaining this pack's exact bytes.
    #[must_use]
    pub const fn artifact_kind(self) -> &'static str {
        match self {
            CardPackKind::Material => "solve-material-card-pack",
            CardPackKind::Interface => "solve-interface-card-pack",
        }
    }

    /// Recover the kind from a retained artifact kind.
    #[must_use]
    pub fn from_artifact_kind(kind: &str) -> Option<CardPackKind> {
        match kind {
            "solve-material-card-pack" => Some(CardPackKind::Material),
            "solve-interface-card-pack" => Some(CardPackKind::Interface),
            _ => None,
        }
    }
}

/// One caller-supplied pack file before admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawCardPack {
    /// Which envelope the invocation declared for these bytes.
    pub kind: CardPackKind,
    /// Caller-facing label (normally the path) retained only for diagnostics.
    pub source: String,
    /// Exact file bytes.
    pub bytes: Vec<u8>,
    /// Optional caller-pinned pack identity. When present the pack decodes
    /// through the verified path and a moved identity refuses.
    pub expect: Option<ContentHash>,
}

/// Structured card-pack admission refusal.
///
/// The shape mirrors the solve driver's refusal so the caller can lift it
/// without inventing new vocabulary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CardPackRefusal {
    /// Stable machine code.
    pub code: &'static str,
    /// What was refused, including the offending source label.
    pub what: String,
    /// Actionable fix.
    pub fix: String,
}

impl CardPackRefusal {
    fn new(code: &'static str, what: impl Into<String>, fix: impl Into<String>) -> CardPackRefusal {
        CardPackRefusal {
            code,
            what: what.into(),
            fix: fix.into(),
        }
    }
}

impl std::fmt::Display for CardPackRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.what)
    }
}

impl std::error::Error for CardPackRefusal {}

/// A decoded pack, retained with the exact bytes that produced it.
#[derive(Debug, Clone, PartialEq)]
enum DecodedPack {
    Material(Box<NormalizedMaterialCardPack>),
    Interface(Box<NormalizedInterfacePack>),
}

/// One admitted pack: canonical bytes, pack root, and the reconstructed card
/// identity that will key the [`CardLibrary`].
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedCardPack {
    pack: DecodedPack,
    root: ContentHash,
    artifact: ContentHash,
    card: ContentHash,
    identity: String,
    bytes: Vec<u8>,
}

impl AdmittedCardPack {
    /// Which envelope this pack used.
    #[must_use]
    pub fn kind(&self) -> CardPackKind {
        match self.pack {
            DecodedPack::Material(_) => CardPackKind::Material,
            DecodedPack::Interface(_) => CardPackKind::Interface,
        }
    }

    /// Semantic content root of the whole pack envelope, as `fs-matdb`
    /// defines it: a domain-separated hash, not a plain hash of the bytes.
    /// This is what the set root — and therefore the run identity — binds.
    #[must_use]
    pub fn root(&self) -> ContentHash {
        self.root
    }

    /// Ledger content address of the exact retained bytes.
    ///
    /// Deliberately distinct from [`AdmittedCardPack::root`]: the ledger
    /// addresses artifacts by plain content hash, so the semantic pack
    /// identity is not an artifact address and the two must never be
    /// substituted for one another.
    #[must_use]
    pub fn artifact(&self) -> ContentHash {
        self.artifact
    }

    /// Content hash of the reconstructed card. This is the [`CardLibrary`]
    /// key and the value a project binding's `:card` field must name.
    #[must_use]
    pub fn card(&self) -> ContentHash {
        self.card
    }

    /// Human-readable card state identity, retained in the stage receipt.
    ///
    /// For a material pack this is exactly the `MaterialStateId` rendering a
    /// project binding's `:state` field must match.
    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }

    /// Exact canonical bytes, retained as a ledger artifact.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

/// The admitted, canonically ordered set of card packs for one invocation.
#[derive(Debug, Clone, PartialEq)]
pub struct CardPackSet {
    materials: Vec<AdmittedCardPack>,
    interfaces: Vec<AdmittedCardPack>,
    root: ContentHash,
}

impl CardPackSet {
    /// The empty set. Its root is well defined, so a project that declares no
    /// bindings still derives a stable run identity.
    #[must_use]
    pub fn empty() -> CardPackSet {
        CardPackSetBuilder::new()
            .finish()
            .unwrap_or_else(|_| unreachable!("the empty pack set cannot conflict"))
    }

    /// Admit a complete set in one call.
    ///
    /// # Errors
    /// [`CardPackRefusal`] when any pack fails to decode, the count ceiling is
    /// exceeded, or two distinct packs reconstruct the same card.
    pub fn admit(raw: Vec<RawCardPack>) -> Result<CardPackSet, CardPackRefusal> {
        let mut builder = CardPackSetBuilder::new();
        for pack in raw {
            builder.push(pack)?;
        }
        builder.finish()
    }

    /// Canonical set root: the ordered pack roots under a separating domain.
    #[must_use]
    pub fn root(&self) -> ContentHash {
        self.root
    }

    /// Admitted material packs in canonical (ascending pack root) order.
    #[must_use]
    pub fn materials(&self) -> &[AdmittedCardPack] {
        &self.materials
    }

    /// Admitted interface packs in canonical (ascending pack root) order.
    #[must_use]
    pub fn interfaces(&self) -> &[AdmittedCardPack] {
        &self.interfaces
    }

    /// Every admitted pack, materials before interfaces, each in canonical
    /// order. This is the retention and edge-linking order.
    pub fn iter(&self) -> impl Iterator<Item = &AdmittedCardPack> {
        self.materials.iter().chain(self.interfaces.iter())
    }

    /// Total admitted packs.
    #[must_use]
    pub fn len(&self) -> usize {
        self.materials.len() + self.interfaces.len()
    }

    /// Whether the invocation supplied no packs at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Build the complete card library.
    ///
    /// Every admitted pack contributes its card with its complete claim set;
    /// no claim is selected or discarded here. Claim selection happens inside
    /// [`fs_project::resolve_bindings`], which leaves replayable usage
    /// receipts for the claims it actually used.
    #[must_use]
    pub fn library(&self) -> CardLibrary {
        let mut library = CardLibrary::new();
        for pack in self.iter() {
            match &pack.pack {
                DecodedPack::Material(material) => {
                    library.insert_material(material.card().clone());
                }
                DecodedPack::Interface(interface) => {
                    library.insert_interface(interface.card().clone());
                }
            }
        }
        library
    }
}

/// Incremental admission so the caller can observe cancellation and charge
/// work between packs.
#[derive(Debug, Default)]
pub struct CardPackSetBuilder {
    materials: Vec<AdmittedCardPack>,
    interfaces: Vec<AdmittedCardPack>,
}

impl CardPackSetBuilder {
    /// A builder holding no packs.
    #[must_use]
    pub fn new() -> CardPackSetBuilder {
        CardPackSetBuilder::default()
    }

    /// Packs admitted so far.
    #[must_use]
    pub fn len(&self) -> usize {
        self.materials.len() + self.interfaces.len()
    }

    /// Whether no pack has been admitted yet.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Decode and admit one pack.
    ///
    /// A byte-identical repeat of an already-admitted pack is idempotent and
    /// returns `Ok` without growing the set.
    ///
    /// # Errors
    /// [`CardPackRefusal`] for an oversized source label, a decode or pinned
    /// identity failure, the count ceiling, or a card conflict with a pack
    /// already admitted.
    pub fn push(&mut self, raw: RawCardPack) -> Result<(), CardPackRefusal> {
        if raw.source.len() > MAX_CARD_PACK_SOURCE_BYTES {
            return Err(CardPackRefusal::new(
                "cli-solve-card-pack-source",
                format!(
                    "a `{}` source label is {} bytes, above the {MAX_CARD_PACK_SOURCE_BYTES}-byte \
                     diagnostic ceiling",
                    raw.kind.flag(),
                    raw.source.len()
                ),
                "supply the pack through a shorter path",
            ));
        }
        if raw.bytes.len() as u64 > MAX_CARD_PACK_BYTES {
            return Err(CardPackRefusal::new(
                "cli-solve-card-pack-size",
                format!(
                    "{} pack `{}` is {} bytes, above the {MAX_CARD_PACK_BYTES}-byte admission \
                     ceiling",
                    raw.kind.label(),
                    raw.source,
                    raw.bytes.len()
                ),
                "split the claim table or raise the documented pack ceiling deliberately",
            ));
        }
        let admitted = decode(&raw)?;
        let existing = match raw.kind {
            CardPackKind::Material => &self.materials,
            CardPackKind::Interface => &self.interfaces,
        };
        // A repeat of exactly the same pack is idempotent. Content addressing
        // makes equal roots equal bytes, so no byte comparison is needed.
        if existing.iter().any(|pack| pack.root == admitted.root) {
            return Ok(());
        }
        // The library is keyed by card hash. Two distinct packs reconstructing
        // one card disagree about provenance under a single key, so neither
        // may win silently.
        if let Some(other) = existing.iter().find(|pack| pack.card == admitted.card) {
            let (first, second) = ordered_pair(other.root, admitted.root);
            return Err(CardPackRefusal::new(
                "cli-solve-card-pack-conflict",
                format!(
                    "{} packs {} and {} reconstruct the same card {}; the card library is keyed \
                     by card identity, so admitting both would silently drop one pack's \
                     provenance",
                    raw.kind.label(),
                    first.to_hex(),
                    second.to_hex(),
                    admitted.card.to_hex()
                ),
                "supply exactly one pack per card identity",
            ));
        }
        if self.len() >= MAX_CARD_PACKS {
            return Err(CardPackRefusal::new(
                "cli-solve-card-pack-count",
                format!("the invocation supplies more than {MAX_CARD_PACKS} distinct card packs"),
                "reduce the pack set to the cards the project actually binds",
            ));
        }
        match raw.kind {
            CardPackKind::Material => self.materials.push(admitted),
            CardPackKind::Interface => self.interfaces.push(admitted),
        }
        Ok(())
    }

    /// Canonicalize and close the set.
    ///
    /// # Errors
    /// [`CardPackRefusal`] is reserved for future whole-set rules; the
    /// current implementation refuses every conflict eagerly in
    /// [`CardPackSetBuilder::push`].
    pub fn finish(mut self) -> Result<CardPackSet, CardPackRefusal> {
        self.materials
            .sort_by(|a, b| a.root.as_bytes().cmp(b.root.as_bytes()));
        self.interfaces
            .sort_by(|a, b| a.root.as_bytes().cmp(b.root.as_bytes()));
        let mut preimage =
            Vec::with_capacity(16 + 32 * (self.materials.len() + self.interfaces.len()));
        push_roots(&mut preimage, &self.materials);
        push_roots(&mut preimage, &self.interfaces);
        let root = hash_domain(CARD_PACK_SET_DOMAIN, &preimage);
        Ok(CardPackSet {
            materials: self.materials,
            interfaces: self.interfaces,
            root,
        })
    }
}

fn push_roots(preimage: &mut Vec<u8>, packs: &[AdmittedCardPack]) {
    preimage.extend_from_slice(&(packs.len() as u64).to_le_bytes());
    for pack in packs {
        preimage.extend_from_slice(pack.root.as_bytes());
    }
}

fn ordered_pair(a: ContentHash, b: ContentHash) -> (ContentHash, ContentHash) {
    if a.as_bytes() <= b.as_bytes() {
        (a, b)
    } else {
        (b, a)
    }
}

fn decode(raw: &RawCardPack) -> Result<AdmittedCardPack, CardPackRefusal> {
    match raw.kind {
        CardPackKind::Material => {
            let pack = match raw.expect {
                Some(expected) => {
                    NormalizedMaterialCardPack::from_bytes_verified(expected, &raw.bytes)
                }
                None => NormalizedMaterialCardPack::from_bytes(&raw.bytes),
            }
            .map_err(|error| decode_refusal(raw, &error))?;
            let card = pack.card().content_hash();
            let identity = pack.card().id().to_string();
            Ok(AdmittedCardPack {
                root: pack.content_hash(),
                artifact: hash_bytes(&raw.bytes),
                pack: DecodedPack::Material(Box::new(pack)),
                card,
                identity,
                bytes: raw.bytes.clone(),
            })
        }
        CardPackKind::Interface => {
            let pack = match raw.expect {
                Some(expected) => {
                    NormalizedInterfacePack::from_bytes_verified(expected, &raw.bytes)
                }
                None => NormalizedInterfacePack::from_bytes(&raw.bytes),
            }
            .map_err(|error| decode_refusal(raw, &error))?;
            let card = pack.card().content_hash();
            let identity = format!(
                "{} | {} | {}",
                pack.card().surface_a().material,
                pack.card().surface_b().material,
                pack.card().medium()
            );
            Ok(AdmittedCardPack {
                root: pack.content_hash(),
                artifact: hash_bytes(&raw.bytes),
                pack: DecodedPack::Interface(Box::new(pack)),
                card,
                identity,
                bytes: raw.bytes.clone(),
            })
        }
    }
}

fn decode_refusal(raw: &RawCardPack, error: &PackError) -> CardPackRefusal {
    let code = match error {
        PackError::IdentityMismatch { .. } => "cli-solve-card-pack-identity",
        PackError::ResourceLimit { .. } => "cli-solve-card-pack-size",
        _ => "cli-solve-card-pack-decode",
    };
    CardPackRefusal::new(
        code,
        format!(
            "{} pack `{}` did not decode through the {} envelope: {error}",
            raw.kind.label(),
            raw.source,
            raw.kind.flag()
        ),
        "supply a pack produced by the fs-matdb normalized pack compiler at the schema version \
         this driver admits",
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedPack, ObservationDataset,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_qty::Dims;

    const FIXTURE_DOMAIN: &str = "org.frankensim.fs-cli.tests.card-pack-fixture.v1";
    const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

    fn provenance() -> Provenance {
        Provenance {
            source: "fixture guarded-hot-plate table".to_string(),
            license: "CC-BY-4.0; redistribution permitted with attribution".to_string(),
            artifact: Some(hash_domain(FIXTURE_DOMAIN, b"fixture-table")),
        }
    }

    /// Build one canonical `FSMCDPK\0` payload.
    ///
    /// `chemistry` selects the material state, so it changes the
    /// reconstructed card. `pack_id` lives only in the nested claim pack
    /// envelope, so changing it alone produces different pack bytes carrying
    /// the *same* card — exactly the conflict the library key cannot absorb.
    fn material_bytes(chemistry: &str, pack_id: &str) -> Vec<u8> {
        let state = MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "wrought".to_string(),
            process: "T6".to_string(),
            revision: 0,
        };
        let mut claims = ClaimSet::new();
        let observation = claims
            .register_observation(ObservationDataset {
                specimen: "fixture coupon".to_string(),
                method: "fixture campaign".to_string(),
                artifact: hash_domain(FIXTURE_DOMAIN, b"raw-observation"),
                caveats: "fixture value; not a seed-dataset authority".to_string(),
                provenance: provenance(),
            })
            .expect("licensed observation inserts");
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new("thermal-conductivity", CONDUCTIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value: 167.0,
                    dims: CONDUCTIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 450.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation],
                provenance: provenance(),
            })
            .expect("conductivity claim inserts");
        let claims_pack = NormalizedPack::new(
            pack_id,
            "frankensim-material-card-pack-compiler-v1",
            hash_domain(FIXTURE_DOMAIN, b"source-envelope"),
            "CC-BY-4.0: redistribution permitted with attribution",
            claims,
            Vec::new(),
            Vec::new(),
        )
        .expect("claim pack admits");
        NormalizedMaterialCardPack::new(state, claims_pack)
            .expect("material-card pack admits")
            .to_bytes()
    }

    fn raw(kind: CardPackKind, source: &str, bytes: Vec<u8>) -> RawCardPack {
        RawCardPack {
            kind,
            source: source.to_string(),
            bytes,
            expect: None,
        }
    }

    #[test]
    fn g0_flag_order_cannot_change_the_set_root() {
        let a = material_bytes("AA6061", "pack-a");
        let b = material_bytes("Cu-OFE", "pack-b");
        let forward = CardPackSet::admit(vec![
            raw(CardPackKind::Material, "a.fsmcdpk", a.clone()),
            raw(CardPackKind::Material, "b.fsmcdpk", b.clone()),
        ])
        .expect("both packs admit");
        let reverse = CardPackSet::admit(vec![
            raw(CardPackKind::Material, "b.fsmcdpk", b),
            raw(CardPackKind::Material, "a.fsmcdpk", a),
        ])
        .expect("both packs admit");
        assert_eq!(forward.root().to_hex(), reverse.root().to_hex());
        assert_eq!(forward.len(), 2);
        let forward_roots: Vec<String> = forward.iter().map(|p| p.root().to_hex()).collect();
        let reverse_roots: Vec<String> = reverse.iter().map(|p| p.root().to_hex()).collect();
        assert_eq!(forward_roots, reverse_roots, "retention order is canonical");
    }

    #[test]
    fn g0_byte_identical_duplicates_are_idempotent() {
        let a = material_bytes("AA6061", "pack-a");
        let once = CardPackSet::admit(vec![raw(CardPackKind::Material, "a.fsmcdpk", a.clone())])
            .expect("admits");
        let twice = CardPackSet::admit(vec![
            raw(CardPackKind::Material, "a.fsmcdpk", a.clone()),
            raw(CardPackKind::Material, "copy/a.fsmcdpk", a),
        ])
        .expect("the repeat is idempotent");
        assert_eq!(twice.len(), 1);
        assert_eq!(once.root().to_hex(), twice.root().to_hex());
    }

    #[test]
    fn g0_distinct_packs_for_one_card_refuse_without_last_one_wins() {
        // Same claims and same material state, different pack envelopes: the
        // reconstructed card is identical, so the library key collides.
        let first = material_bytes("AA6061", "pack-a");
        let second = material_bytes("AA6061", "pack-b");
        assert_ne!(first, second, "the fixture must differ in envelope bytes");
        let refusal = CardPackSet::admit(vec![
            raw(CardPackKind::Material, "a.fsmcdpk", first),
            raw(CardPackKind::Material, "b.fsmcdpk", second),
        ])
        .expect_err("one card cannot come from two packs");
        assert_eq!(refusal.code, "cli-solve-card-pack-conflict");
    }

    #[test]
    fn g0_a_moved_pinned_identity_refuses() {
        let bytes = material_bytes("AA6061", "pack-a");
        let wrong = hash_domain("org.frankensim.tests.wrong-pin.v1", b"not-this-pack");
        let refusal = CardPackSet::admit(vec![RawCardPack {
            kind: CardPackKind::Material,
            source: "a.fsmcdpk".to_string(),
            bytes: bytes.clone(),
            expect: Some(wrong),
        }])
        .expect_err("a moved identity refuses");
        assert_eq!(refusal.code, "cli-solve-card-pack-identity");

        let set = CardPackSet::admit(vec![RawCardPack {
            kind: CardPackKind::Material,
            source: "a.fsmcdpk".to_string(),
            bytes,
            expect: Some(
                CardPackSet::admit(vec![raw(
                    CardPackKind::Material,
                    "a.fsmcdpk",
                    material_bytes("AA6061", "pack-a"),
                )])
                .expect("admits")
                .materials()[0]
                    .root(),
            ),
        }])
        .expect("the exact pinned identity admits");
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn g0_wrong_magic_and_truncation_refuse_by_decode() {
        let bytes = material_bytes("AA6061", "pack-a");
        let mut wrong_magic = bytes.clone();
        wrong_magic[0] = b'X';
        let refusal =
            CardPackSet::admit(vec![raw(CardPackKind::Material, "a", wrong_magic)]).expect_err("");
        assert_eq!(refusal.code, "cli-solve-card-pack-decode");

        let truncated = bytes[..bytes.len() / 2].to_vec();
        let refusal =
            CardPackSet::admit(vec![raw(CardPackKind::Material, "a", truncated)]).expect_err("");
        assert_eq!(refusal.code, "cli-solve-card-pack-decode");

        // An interface flag over material bytes is a decode refusal, not a
        // silent kind coercion.
        let refusal = CardPackSet::admit(vec![raw(CardPackKind::Interface, "a", bytes)])
            .expect_err("kind is not inferred from content");
        assert_eq!(refusal.code, "cli-solve-card-pack-decode");
    }

    #[test]
    fn g0_the_source_label_ceiling_admits_the_cap_and_refuses_one_byte_past() {
        // The label is retained only for diagnostics, so the ceiling exists to
        // bound what a refusal message can echo back. Exactly at the cap is
        // still a legal label; the refusal starts one byte later.
        let bytes = material_bytes("AA6061", "pack-a");
        let at_cap = "s".repeat(MAX_CARD_PACK_SOURCE_BYTES);
        let set = CardPackSet::admit(vec![raw(CardPackKind::Material, &at_cap, bytes.clone())])
            .expect("a label exactly at the ceiling is still admitted");
        assert_eq!(set.len(), 1);

        let past_cap = "s".repeat(MAX_CARD_PACK_SOURCE_BYTES + 1);
        let refusal = CardPackSet::admit(vec![raw(CardPackKind::Material, &past_cap, bytes)])
            .expect_err("one byte past the ceiling refuses");
        assert_eq!(refusal.code, "cli-solve-card-pack-source");
        assert!(
            !refusal.what.contains(&past_cap),
            "a refusal about an oversized label must not echo the whole label back"
        );
    }

    #[test]
    fn g0_the_pack_byte_ceiling_refuses_before_decode_is_attempted() {
        // Both lengths are undecodable garbage, so the *code* is what
        // discriminates: at the cap the bytes reach the decoder and fail on
        // magic, one byte past the ceiling refuses first. That is the exact
        // boundary and the guard ordering in a single comparison.
        let at_cap = vec![0u8; MAX_CARD_PACK_BYTES as usize];
        let refusal = CardPackSet::admit(vec![raw(CardPackKind::Material, "at-cap", at_cap)])
            .expect_err("undecodable bytes never admit");
        assert_eq!(
            refusal.code, "cli-solve-card-pack-decode",
            "bytes exactly at the ceiling are inside admission and reach the decoder"
        );

        let past_cap = vec![0u8; MAX_CARD_PACK_BYTES as usize + 1];
        let refusal = CardPackSet::admit(vec![raw(CardPackKind::Material, "past-cap", past_cap)])
            .expect_err("one byte past the ceiling refuses");
        assert_eq!(
            refusal.code, "cli-solve-card-pack-size",
            "the size ceiling must fire before the decoder is handed the bytes"
        );
    }

    #[test]
    fn g0_the_pack_count_ceiling_admits_exactly_max_card_packs_and_refuses_the_next() {
        // Distinct chemistries give distinct material states, so every fixture
        // reconstructs a different card and none collapses as a duplicate or
        // refuses as a conflict.
        let mut packs = Vec::with_capacity(MAX_CARD_PACKS + 1);
        for index in 0..=MAX_CARD_PACKS {
            packs.push(raw(
                CardPackKind::Material,
                &format!("pack-{index}.fsmcdpk"),
                material_bytes(&format!("Alloy{index:04}"), "pack-a"),
            ));
        }

        let at_cap = CardPackSet::admit(packs[..MAX_CARD_PACKS].to_vec())
            .expect("exactly the ceiling admits");
        assert_eq!(at_cap.len(), MAX_CARD_PACKS);

        let refusal = CardPackSet::admit(packs).expect_err("one pack past the ceiling refuses");
        assert_eq!(refusal.code, "cli-solve-card-pack-count");
    }

    #[test]
    fn g0_idempotent_repeats_do_not_consume_the_pack_count_budget() {
        // The ceiling bounds the *admitted set*, not the number of times a
        // caller names a pack: the idempotent-repeat path returns before the
        // count check. This is the semantic difference between this ceiling
        // and the invocation-level one the CLI grammar applies, and it is why
        // the two carry different messages under the same code.
        let mut packs = Vec::new();
        for index in 0..MAX_CARD_PACKS {
            packs.push(raw(
                CardPackKind::Material,
                &format!("pack-{index}.fsmcdpk"),
                material_bytes(&format!("Alloy{index:04}"), "pack-a"),
            ));
        }
        let distinct = CardPackSet::admit(packs.clone()).expect("the full set admits");

        // Every pack named a second time: twice the declarations, the same
        // admitted set and the same canonical root.
        let mut repeated = packs.clone();
        repeated.extend(packs);
        let collapsed =
            CardPackSet::admit(repeated).expect("repeats collapse instead of exhausting the count");
        assert_eq!(collapsed.len(), MAX_CARD_PACKS);
        assert_eq!(distinct.root().to_hex(), collapsed.root().to_hex());
    }

    #[test]
    fn g0_the_empty_set_has_a_stable_root() {
        assert_eq!(
            CardPackSet::empty().root().to_hex(),
            CardPackSet::admit(Vec::new())
                .expect("admits")
                .root()
                .to_hex()
        );
        assert!(CardPackSet::empty().is_empty());
    }

    #[test]
    fn g0_the_library_carries_every_admitted_card() {
        let set = CardPackSet::admit(vec![
            raw(
                CardPackKind::Material,
                "a",
                material_bytes("AA6061", "pack-a"),
            ),
            raw(
                CardPackKind::Material,
                "b",
                material_bytes("Cu-OFE", "pack-b"),
            ),
        ])
        .expect("admits");
        let library = set.library();
        for pack in set.iter() {
            assert!(
                library.material(&pack.card().to_hex()).is_some(),
                "every admitted card is keyed by its own content hash"
            );
        }
    }
}
