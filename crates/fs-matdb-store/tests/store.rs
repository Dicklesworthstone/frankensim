//! Conformance battery for the FrankenSQLite material store
//! (bead frankensim-oecdy): discovery correctness, evaluation parity
//! by construction, staleness/tamper fail-closed, and bitwise rebuild.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MODEL_PACK_TARGET_BASIS, MatDbError, MaterialStateId, ModelNormalizationReceipt,
    ModelNormalizationTarget, NormalizedInterfacePack, NormalizedMaterialCardPack,
    NormalizedModelPack, NormalizedPack, NormalizedSpeciesPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SPECIES_MOLAR_MASS_DIMS,
    SPECIES_PACK_TARGET_BASIS, SPECIES_REFERENCE_PRESSURE_DIMS, SelectionPolicy,
    SpeciesAssociation, SpeciesNormalizationReceipt, SpeciesNormalizationTarget, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_matdb_store::{CatalogPack, MaterialStore, PackKind, STORE_SCHEMA_VERSION, StoreError};
use fs_qty::Dims;
use fsqlite::{AsyncConnection, SqliteValue};

fn dims_none() -> Dims {
    Dims::NONE
}

fn provenance() -> Provenance {
    Provenance {
        source: "Store battery synthetic source".to_string(),
        license: "CC-BY-4.0".to_string(),
        artifact: None,
    }
}

/// A synthetic pack with the given scalar properties, each valid on
/// `temperature` in [200, 400].
fn test_pack(pack_id: &str, properties: &[(&str, f64)]) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: format!("{pack_id} synthetic specimen"),
            method: "authored".to_string(),
            artifact: fs_blake3_hash(pack_id.as_bytes()),
            caveats: "synthetic battery data".to_string(),
            provenance: provenance(),
        })
        .expect("observation");
    for &(name, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims_none()),
                value: PropertyValue::Scalar {
                    value,
                    dims: dims_none(),
                },
                validity: ValidityDomain::unconstrained().with("temperature", 200.0, 400.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation],
                provenance: provenance(),
            })
            .expect("claim");
    }
    NormalizedPack::new(
        pack_id,
        "store-battery-v1",
        fs_blake3_hash(b"synthetic source artifact"),
        "synthetic redistribution permitted for tests",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .expect("pack")
}

fn fs_blake3_hash(bytes: &[u8]) -> fs_blake3::ContentHash {
    fs_blake3::hash_bytes(bytes)
}

fn scratch_path(tag: &str) -> String {
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let sequence = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "fs-matdb-store-{tag}-{}-{nonce}-{sequence}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).expect("fresh scratch dir without deleting previous evidence");
    dir.join("store.db").to_string_lossy().into_owned()
}

fn seeded_store(path: &str) -> MaterialStore {
    let store = MaterialStore::open(path).expect("open");
    store
        .ingest_pack(&test_pack(
            "steel-304-synth",
            &[("young_modulus", 193.0e9), ("melting_point", 1673.0)],
        ))
        .expect("ingest steel");
    store
        .ingest_pack(&test_pack(
            "balsa-synth",
            &[("young_modulus", 3.4e9), ("density", 160.0)],
        ))
        .expect("ingest balsa");
    store.seal_corpus().expect("seal");
    store
}

#[test]
fn discovery_surfaces_answer_the_owner_questions() {
    let path = scratch_path("discovery");
    let store = seeded_store(&path);
    // Everything about one material.
    let steel = store.properties_of("steel-304-synth").expect("props");
    assert_eq!(steel.len(), 2);
    assert!(steel.iter().any(|r| r.property == "melting_point"));
    // Every material with a property, range-filtered.
    let stiff = store
        .materials_with("young_modulus", Some((100.0e9, 300.0e9)))
        .expect("range");
    assert_eq!(stiff.len(), 1);
    assert_eq!(stiff[0].pack_id, "steel-304-synth");
    let all = store.materials_with("young_modulus", None).expect("all");
    assert_eq!(all.len(), 2);
    // Validity-window discovery: both packs valid at 300 K, none at
    // 500 K.
    assert_eq!(
        store
            .valid_at("young_modulus", "temperature", 300.0)
            .expect("valid")
            .len(),
        2
    );
    assert_eq!(
        store
            .valid_at("young_modulus", "temperature", 500.0)
            .expect("valid")
            .len(),
        0
    );
    // Unknown pack refuses by name.
    assert!(matches!(
        store.properties_of("unobtainium"),
        Err(StoreError::UnknownPack { .. })
    ));
    println!("{{\"suite\":\"fs-matdb-store\",\"case\":\"discovery\",\"verdict\":\"pass\"}}");
}

#[test]
fn evaluation_is_the_same_answer_as_direct_pack_use() {
    // Parity by construction, ASSERTED: the store's answer must be
    // bitwise the in-memory answer, receipt hash included.
    let path = scratch_path("parity");
    let store = seeded_store(&path);
    let pack = test_pack(
        "steel-304-synth",
        &[("young_modulus", 193.0e9), ("melting_point", 1673.0)],
    );
    let point = QueryPoint::new().with("temperature", 300.0).expect("point");
    let direct = pack
        .claims()
        .query("young_modulus", &point, SelectionPolicy::SingleClaimOnly)
        .expect("direct");
    let stored = store
        .evaluate(
            "steel-304-synth",
            "young_modulus",
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("stored");
    assert_eq!(
        stored.evidence.value.value.to_bits(),
        direct.evidence.value.value.to_bits()
    );
    assert_eq!(
        stored.receipt.content_hash().to_hex(),
        direct.receipt.content_hash().to_hex()
    );
    // Extrapolation refusal passes through UNCHANGED.
    let outside = QueryPoint::new().with("temperature", 500.0).expect("point");
    let refusal = store.evaluate(
        "steel-304-synth",
        "young_modulus",
        &outside,
        SelectionPolicy::SingleClaimOnly,
    );
    assert!(matches!(
        refusal,
        Err(StoreError::MatDb(MatDbError::NoClaimInDomain { .. }))
    ));
    // Absent data is a named refusal, never a fabricated row.
    let missing = store.evaluate(
        "balsa-synth",
        "melting_point",
        &point,
        SelectionPolicy::SingleClaimOnly,
    );
    assert!(matches!(
        missing,
        Err(StoreError::MatDb(MatDbError::UnknownProperty { .. }))
    ));
    println!(
        "{{\"suite\":\"fs-matdb-store\",\"case\":\"evaluation-parity\",\"verdict\":\"pass\"}}"
    );
}

#[test]
fn staleness_fails_closed_and_resealing_is_deliberate() {
    let path = scratch_path("staleness");
    let store = seeded_store(&path);
    assert!(store.require_sealed().is_ok());
    // A pack ingested after sealing makes EVERY surface refuse.
    store
        .ingest_pack(&test_pack("late-arrival", &[("density", 1000.0)]))
        .expect("late ingest");
    assert!(matches!(
        store.require_sealed(),
        Err(StoreError::CorpusChanged { .. })
    ));
    assert!(matches!(
        store.properties_of("steel-304-synth"),
        Err(StoreError::CorpusChanged { .. })
    ));
    assert!(matches!(
        store.evaluate(
            "steel-304-synth",
            "young_modulus",
            &QueryPoint::new().with("temperature", 300.0).expect("point"),
            SelectionPolicy::SingleClaimOnly,
        ),
        Err(StoreError::CorpusChanged { .. })
    ));
    // Re-sealing is an explicit act that restores service.
    store.seal_corpus().expect("re-seal");
    assert!(store.require_sealed().is_ok());
    // A fresh unsealed store refuses by name.
    let empty_path = scratch_path("unsealed");
    let empty = MaterialStore::open(&empty_path).expect("open");
    assert!(matches!(
        empty.materials_with("density", None),
        Err(StoreError::NotSealed)
    ));
    println!("{{\"suite\":\"fs-matdb-store\",\"case\":\"staleness\",\"verdict\":\"pass\"}}");
}

#[test]
fn index_tampering_is_detected_and_cannot_poison_answers() {
    let path = scratch_path("tamper");
    let store = seeded_store(&path);
    assert!(store.verify_index("steel-304-synth").is_ok());
    // Tamper with the DERIVED index through a raw connection.
    let mut raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_with_params_sync(
        "UPDATE claims SET scalar_value = ?1 WHERE property = 'young_modulus' \
         AND pack_id = 'steel-304-synth'",
        &[SqliteValue::Float(999.0)],
    )
    .expect("tamper");
    raw.close_sync().expect("close raw");
    // (a) The cross-check catches it by name.
    assert!(matches!(
        store.verify_index("steel-304-synth"),
        Err(StoreError::IndexMismatch {
            what: "scalar value",
            ..
        })
    ));
    // (b) THE LOAD-BEARING CLAIM: evaluation is untouched, because it
    // never reads the index — the canonical hash-verified bytes are
    // the only evaluation path.
    let answer = store
        .evaluate(
            "steel-304-synth",
            "young_modulus",
            &QueryPoint::new().with("temperature", 300.0).expect("point"),
            SelectionPolicy::SingleClaimOnly,
        )
        .expect("evaluate");
    assert_eq!(answer.evidence.value.value.to_bits(), 193.0e9f64.to_bits());
    // (b2) Review finding: valid_at is DRIVEN by the validity table,
    // so verify_index must cover it too.
    let mut raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_sync("UPDATE validity SET lo = 0.0, hi = 1.0e9 WHERE pack_id = 'steel-304-synth'")
        .expect("validity tamper");
    raw.close_sync().expect("close raw");
    assert!(matches!(
        store.verify_index("steel-304-synth"),
        Err(StoreError::IndexMismatch {
            what: "validity bounds" | "scalar value",
            ..
        })
    ));
    // (c) Tampering with the CANONICAL BYTES is caught by the content
    // hash at decode time.
    let mut raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_sync("UPDATE packs SET bytes = x'deadbeef' WHERE pack_id = 'steel-304-synth'")
        .expect("corrupt");
    raw.close_sync().expect("close raw");
    assert!(matches!(
        store.evaluate(
            "steel-304-synth",
            "young_modulus",
            &QueryPoint::new().with("temperature", 300.0).expect("point"),
            SelectionPolicy::SingleClaimOnly,
        ),
        Err(StoreError::PackCorrupt { .. })
    ));
    println!("{{\"suite\":\"fs-matdb-store\",\"case\":\"tamper\",\"verdict\":\"pass\"}}");
}

#[test]
fn rebuild_is_bitwise_identical() {
    let a = seeded_store(&scratch_path("rebuild-a"));
    let b = seeded_store(&scratch_path("rebuild-b"));
    assert_eq!(
        a.seal_corpus().expect("digest a").to_hex(),
        b.seal_corpus().expect("digest b").to_hex()
    );
    assert_eq!(
        a.canonical_dump().expect("dump a"),
        b.canonical_dump().expect("dump b")
    );
    println!("{{\"suite\":\"fs-matdb-store\",\"case\":\"bitwise-rebuild\",\"verdict\":\"pass\"}}");
}

#[test]
fn partial_ingest_rolls_back_atomically() {
    // Review finding follow-through, with an executed twist: an
    // unlicensed claim turned out to be UNREPRESENTABLE — fs-matdb's
    // own ClaimSet::insert_claim refuses it (MissingLicense), so the
    // store's admission pre-pass is defense-in-depth for future pack
    // paths, and this test pins the upstream gate...
    let mut claims = ClaimSet::new();
    let refused = claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("young_modulus", dims_none()),
        value: PropertyValue::Scalar {
            value: 1.0e9,
            dims: dims_none(),
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "unlicensed".to_string(),
            license: "  ".to_string(),
            artifact: None,
        },
    });
    assert!(
        refused.is_err(),
        "fs-matdb itself must refuse an unlicensed claim"
    );
    // ...and the TRANSACTIONAL machinery is exercised by an induced
    // mid-ingest SQL failure: with the validity table dropped, the
    // third insert of a fresh pack fails and the whole pack — packs
    // row and claims rows — must roll back, leaving a retry able to
    // succeed once the schema is restored.
    let path = scratch_path("rollback");
    let store = MaterialStore::open(&path).expect("open");
    let mut raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_sync("DROP TABLE validity")
        .expect("drop validity");
    raw.close_sync().expect("close raw");
    assert!(matches!(
        store.ingest_pack(&test_pack("rollback-pack", &[("density", 1000.0)])),
        Err(StoreError::Sql { .. })
    ));
    // Restore the schema; the retry must NOT hit DuplicatePack —
    // nothing of the failed ingest survived.
    let mut raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_sync(
        "CREATE TABLE validity(pack_id TEXT NOT NULL, claim_hash BLOB NOT NULL, \
         axis TEXT NOT NULL, lo REAL NOT NULL, hi REAL NOT NULL)",
    )
    .expect("recreate validity");
    raw.close_sync().expect("close raw");
    store
        .ingest_pack(&test_pack("rollback-pack", &[("density", 1000.0)]))
        .expect("retry succeeds — the failed ingest left no residue");
    store.seal_corpus().expect("seal");
    assert_eq!(
        store.properties_of("rollback-pack").expect("props").len(),
        1
    );
    store
        .verify_index("rollback-pack")
        .expect("index consistent");
    println!("{{\"suite\":\"fs-matdb-store\",\"case\":\"atomic-rollback\",\"verdict\":\"pass\"}}");
}

#[test]
fn unknown_property_refuses_by_name() {
    // Review finding: a typo'd property must be a NAMED refusal, not
    // an empty result set.
    let path = scratch_path("unknown-prop");
    let store = seeded_store(&path);
    assert!(matches!(
        store.materials_with("yuong_modulus", None),
        Err(StoreError::UnknownProperty { .. })
    ));
    assert!(matches!(
        store.valid_at("yuong_modulus", "temperature", 300.0),
        Err(StoreError::UnknownProperty { .. })
    ));
    // An empty RANGE result on a KNOWN property is a legitimate empty
    // set, not a refusal.
    assert!(
        store
            .materials_with("young_modulus", Some((1.0, 2.0)))
            .expect("known property, empty range")
            .is_empty()
    );
}

#[test]
fn duplicates_refuse() {
    let path = scratch_path("dup");
    let store = seeded_store(&path);
    assert!(matches!(
        store.ingest_pack(&test_pack("balsa-synth", &[("density", 160.0)])),
        Err(StoreError::DuplicatePack { .. })
    ));
}

fn named_state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "synthetic-test-specimen".into(),
        revision: 0,
    }
}

/// Real canonical codecs and real SQL storage, with explicitly synthetic data.
/// These fixtures verify transport/evaluation, not measured material accuracy.
fn family_bundle() -> Vec<CatalogPack> {
    let property = test_pack("unbound-claims", &[("density", 1000.0)]);
    let material = NormalizedMaterialCardPack::new(
        named_state("body"),
        test_pack("body-card", &[("young_modulus", 2.0e9)]),
    )
    .unwrap();
    let interface = NormalizedInterfacePack::new(
        SurfaceSpec {
            material: named_state("body"),
            texture_frame: "body/grain".into(),
        },
        SurfaceSpec {
            material: named_state("base"),
            texture_frame: "base/polish".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "air".into(),
            history: "virgin".into(),
        },
        test_pack("ordered-contact", &[("friction", 0.25)]),
    )
    .unwrap();
    let (models, species) = model_and_species();
    vec![
        CatalogPack::Properties(property),
        CatalogPack::MaterialCard(material),
        CatalogPack::Interface(interface),
        CatalogPack::Model(models),
        CatalogPack::Species(species),
    ]
}

fn model_and_species() -> (NormalizedModelPack, NormalizedSpeciesPack) {
    let source = fs_blake3_hash(b"synthetic-model-and-species-source");
    let provenance = Provenance {
        artifact: Some(source),
        ..provenance()
    };
    let model = ConstitutiveModelCard {
        law: LawId("synthetic-law".into()),
        law_version: 7,
        parameters: std::collections::BTreeMap::from([(
            "stiffness".into(),
            LawParameter {
                value: 3.0,
                dims: Dims([0, 1, -2, 0, 0, 0]),
            },
        )]),
        state_schema_version: 2,
        initial_state: InitialStatePolicy::RequiresDeclaredState,
        validity: ValidityDomain::unconstrained(),
        sources: vec![source],
        provenance: provenance.clone(),
    };
    let receipt = ModelNormalizationReceipt::new(
        ModelNormalizationTarget::Parameter {
            model: model.content_hash(),
            parameter: "stiffness".into(),
        },
        source,
        Dims([0, 1, -2, 0, 0, 0]),
        1.0,
        0.0,
        "N/m",
        MODEL_PACK_TARGET_BASIS,
        None,
        None,
    );
    let models = NormalizedModelPack::new(
        "constitutive-models",
        "test-compiler",
        source,
        "test redistribution",
        vec![model],
        vec![receipt],
    )
    .unwrap();
    let association = SpeciesAssociation::new(
        fs_qty::chemistry::SpeciesId::new("N2").unwrap(),
        0.028_013_4,
        "gas",
        "ideal-gas",
        100_000.0,
        "NASA-TP-2002-211556",
        vec![source],
        provenance,
    )
    .unwrap();
    let normalizations = vec![
        SpeciesNormalizationReceipt::new(
            SpeciesNormalizationTarget::MolarMass,
            source,
            SPECIES_MOLAR_MASS_DIMS,
            1.0,
            0.0,
            "kg/mol",
            SPECIES_PACK_TARGET_BASIS,
        ),
        SpeciesNormalizationReceipt::new(
            SpeciesNormalizationTarget::ReferencePressure,
            source,
            SPECIES_REFERENCE_PRESSURE_DIMS,
            1.0,
            0.0,
            "Pa",
            SPECIES_PACK_TARGET_BASIS,
        ),
    ];
    let species = NormalizedSpeciesPack::new(
        "N2",
        "test-compiler",
        source,
        "test redistribution",
        association,
        normalizations,
    )
    .unwrap();
    (models, species)
}

fn assert_family_round_trip(store: &MaterialStore, pack: &CatalogPack, point: &QueryPoint) {
    let discovered = store.packs(Some(pack.kind())).unwrap();
    assert_eq!(discovered.len(), 1);
    assert_eq!(discovered[0].content_hash, pack.content_hash());
    let loaded = store
        .load_by_hash(pack.kind(), pack.content_hash())
        .unwrap();
    assert_eq!(&loaded, pack);
    assert_eq!(loaded.to_bytes(), pack.to_bytes());
    store.verify_index(pack.pack_id()).unwrap();
    if let Some(claims) = pack.claims_pack() {
        let name = claims
            .claims()
            .claims_ordered()
            .next()
            .unwrap()
            .1
            .key
            .name();
        let direct = claims
            .claims()
            .query(name, point, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        let stored = store
            .evaluate(
                pack.pack_id(),
                name,
                point,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(
            stored.evidence.value.value.to_bits(),
            direct.evidence.value.value.to_bits()
        );
        assert_eq!(stored.receipt.content_hash(), direct.receipt.content_hash());
        assert_eq!(
            store.properties_of(pack.pack_id()).unwrap()[0].pack_kind,
            pack.kind()
        );
    } else {
        assert!(store.properties_of(pack.pack_id()).unwrap().is_empty());
        assert!(matches!(
            store.evaluate(
                pack.pack_id(),
                "stiffness",
                point,
                SelectionPolicy::SingleClaimOnly
            ),
            Err(StoreError::NoPropertyClaims { .. })
        ));
    }
}

#[test]
fn g0_all_families_round_trip_and_resolve_body_and_ordered_interface() {
    let path = scratch_path("families");
    let bundle = family_bundle();
    {
        let store = MaterialStore::open(&path).unwrap();
        store.ingest_bundle(&bundle).unwrap();
        store.seal_corpus().unwrap();
    }
    let store = MaterialStore::open(&path).unwrap();
    assert_eq!(store.packs(None).unwrap().len(), 5);
    let point = QueryPoint::new().with("temperature", 300.0).unwrap();
    for pack in &bundle {
        assert_family_round_trip(&store, pack, &point);
    }
    let CatalogPack::MaterialCard(body) = store
        .load_by_hash(PackKind::MaterialCard, bundle[1].content_hash())
        .unwrap()
    else {
        panic!("material family");
    };
    let CatalogPack::Interface(contact) = store
        .load_by_hash(PackKind::Interface, bundle[2].content_hash())
        .unwrap()
    else {
        panic!("interface family");
    };
    assert_eq!(contact.card().surface_a().material, *body.card().id());
    assert_eq!(contact.card().surface_b().material, named_state("base"));
    assert_eq!(contact.card().surface_a().texture_frame, "body/grain");
    assert_eq!(contact.card().history(), "virgin");
    let reversed = NormalizedInterfacePack::new(
        contact.card().surface_b().clone(),
        contact.card().surface_a().clone(),
        contact.card().context().clone(),
        contact.claims_pack().clone(),
    )
    .unwrap();
    assert_ne!(
        reversed.content_hash(),
        contact.content_hash(),
        "surface order remains identity-bearing"
    );
    assert!(matches!(
        store.load_by_hash(PackKind::Interface, reversed.content_hash()),
        Err(StoreError::UnknownContent { .. })
    ));
    assert!(matches!(
        store.load_by_hash(PackKind::Model, bundle[1].content_hash()),
        Err(StoreError::WrongPackKind { .. })
    ));
    assert!(matches!(
        store.load_pack("body-card"),
        Err(StoreError::WrongPackKind { .. })
    ));
}

#[test]
fn g4_bundle_failure_rolls_back_earlier_families_and_preserves_the_seal() {
    let store = MaterialStore::open(":memory:").unwrap();
    let baseline = test_pack("already-present", &[("density", 42.0)]);
    store.ingest_pack(&baseline).unwrap();
    let seal = store.seal_corpus().unwrap();
    let dump = store.canonical_dump().unwrap();
    let mut bundle = family_bundle();
    // The final existing-id conflict occurs after all five distinct families
    // have been written inside the transaction. It must roll all of them back.
    bundle.push(CatalogPack::Properties(baseline));
    assert!(matches!(
        store.ingest_bundle(&bundle),
        Err(StoreError::DuplicatePack { .. })
    ));
    assert_eq!(store.canonical_dump().unwrap(), dump);
    store.require_sealed().unwrap();
    assert_eq!(store.seal_corpus().unwrap(), seal);
    assert_eq!(store.packs(None).unwrap().len(), 1);
    assert!(matches!(
        store.load_catalog_pack("body-card"),
        Err(StoreError::UnknownPack { .. })
    ));
    store.ingest_bundle(&bundle[..5]).unwrap();
    assert!(matches!(
        store.packs(None),
        Err(StoreError::CorpusChanged { .. })
    ));
    store.seal_corpus().unwrap();
    assert_eq!(store.packs(None).unwrap().len(), 6);
}

#[test]
fn g5_family_bundle_rebuild_is_independent_of_ingest_order() {
    let a = MaterialStore::open(":memory:").unwrap();
    let b = MaterialStore::open(":memory:").unwrap();
    let mut bundle = family_bundle();
    a.ingest_bundle(&bundle).unwrap();
    bundle.reverse();
    b.ingest_bundle(&bundle).unwrap();
    assert_eq!(a.seal_corpus().unwrap(), b.seal_corpus().unwrap());
    assert_eq!(a.canonical_dump().unwrap(), b.canonical_dump().unwrap());
    a.ingest_bundle(&[]).unwrap();
    a.require_sealed().unwrap();
}

#[test]
fn g0_family_hash_and_wire_version_refuse_before_ingestion() {
    for pack in family_bundle() {
        assert!(
            CatalogPack::from_bytes_verified(
                pack.kind(),
                fs_blake3_hash(b"wrong hash"),
                &pack.to_bytes()
            )
            .is_err()
        );
        let wrong_kind = if pack.kind() == PackKind::Properties {
            PackKind::Model
        } else {
            PackKind::Properties
        };
        assert!(
            CatalogPack::from_bytes_verified(wrong_kind, pack.content_hash(), &pack.to_bytes())
                .is_err()
        );
        let mut bytes = pack.to_bytes();
        bytes[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
        // Recompute the exact family hash so rejection must come from the
        // version-aware decoder, rather than merely a mismatched old digest.
        let domain = match pack.kind() {
            PackKind::Properties => "org.frankensim.fs-matdb.normalized-pack.v1",
            PackKind::MaterialCard => "org.frankensim.fs-matdb.normalized-material-card-pack.v1",
            PackKind::Interface => "org.frankensim.fs-matdb.normalized-interface-pack.v1",
            PackKind::Model => "org.frankensim.fs-matdb.normalized-model-pack.v1",
            PackKind::Species => "org.frankensim.fs-matdb.normalized-species-pack.v1",
        };
        let hash = fs_blake3::hash_domain(domain, &bytes);
        assert!(
            matches!(
                CatalogPack::from_bytes_verified(pack.kind(), hash, &bytes),
                Err(fs_matdb::PackError::Malformed { .. })
            ),
            "{:?}",
            pack.kind()
        );
    }
}

#[test]
fn g0_v1_migration_preserves_canonical_bytes_and_requires_explicit_reseal() {
    let path = scratch_path("migration");
    let pack = test_pack("legacy", &[("density", 1000.0)]);
    let mut raw = AsyncConnection::open_sync(&path).unwrap();
    raw.execute_sync("CREATE TABLE packs(pack_id TEXT PRIMARY KEY, content_hash BLOB NOT NULL, bytes BLOB NOT NULL, compiler TEXT NOT NULL, redistribution TEXT NOT NULL)").unwrap();
    raw.execute_sync("CREATE TABLE claims(pack_id TEXT NOT NULL, claim_hash BLOB NOT NULL, property TEXT NOT NULL, kind TEXT NOT NULL, scalar_value REAL, observation_backed INTEGER NOT NULL, license TEXT NOT NULL, source TEXT NOT NULL)").unwrap();
    raw.execute_sync("CREATE TABLE validity(pack_id TEXT NOT NULL, claim_hash BLOB NOT NULL, axis TEXT NOT NULL, lo REAL NOT NULL, hi REAL NOT NULL)").unwrap();
    raw.execute_sync("CREATE TABLE corpus_seal(singleton INTEGER PRIMARY KEY, digest BLOB NOT NULL, pack_count INTEGER NOT NULL)").unwrap();
    raw.execute_with_params_sync(
        "INSERT INTO packs VALUES (?1, ?2, ?3, ?4, ?5)",
        &[
            SqliteValue::Text(pack.pack_id().into()),
            SqliteValue::Blob(pack.content_hash().as_bytes().to_vec().into()),
            SqliteValue::Blob(pack.to_bytes().into()),
            SqliteValue::Text(pack.compiler().into()),
            SqliteValue::Text(pack.redistribution_terms().into()),
        ],
    )
    .unwrap();
    let (claim_id, claim) = pack.claims().claims_ordered().next().unwrap();
    raw.execute_with_params_sync(
        "INSERT INTO claims VALUES ('legacy', ?1, 'density', 'scalar', 1000.0, 1, ?2, ?3)",
        &[
            SqliteValue::Blob(claim_id.0.as_bytes().to_vec().into()),
            SqliteValue::Text(claim.provenance.license.clone().into()),
            SqliteValue::Text(claim.provenance.source.clone().into()),
        ],
    )
    .unwrap();
    raw.execute_with_params_sync(
        "INSERT INTO validity VALUES ('legacy', ?1, 'temperature', 200.0, 400.0)",
        &[SqliteValue::Blob(claim_id.0.as_bytes().to_vec().into())],
    )
    .unwrap();
    let mut old = fs_blake3::DomainHasher::new("org.frankensim.fs-matdb-store.corpus.v1");
    old.update(&(pack.pack_id().len() as u64).to_le_bytes());
    old.update(pack.pack_id().as_bytes());
    old.update(&32_u64.to_le_bytes());
    old.update(pack.content_hash().as_bytes());
    raw.execute_with_params_sync(
        "INSERT INTO corpus_seal VALUES (1, ?1, 1)",
        &[SqliteValue::Blob(old.finalize().as_bytes().to_vec().into())],
    )
    .unwrap();
    raw.execute_sync("PRAGMA user_version = 1").unwrap();
    raw.close_sync().unwrap();
    let store = MaterialStore::open(&path).unwrap();
    assert!(matches!(
        store.load_pack("legacy"),
        Err(StoreError::CorpusChanged { .. })
    ));
    store.seal_corpus().unwrap();
    assert_eq!(store.load_pack("legacy").unwrap(), pack);
    store.verify_index("legacy").unwrap();
    assert_eq!(store.properties_of("legacy").unwrap().len(), 1);
    assert_eq!(
        store
            .valid_at("density", "temperature", 300.0)
            .unwrap()
            .len(),
        1
    );
    assert!(
        store
            .valid_at("density", "temperature", 500.0)
            .unwrap()
            .is_empty()
    );
    assert_eq!(store.packs(None).unwrap()[0].kind, PackKind::Properties);
    drop(store);
    let mut raw = AsyncConnection::open_sync(&path).unwrap();
    assert!(
        matches!(raw.query_sync("PRAGMA user_version").unwrap()[0].get(0), Some(SqliteValue::Integer(v)) if *v == STORE_SCHEMA_VERSION)
    );
    raw.execute_sync("PRAGMA user_version = 99").unwrap();
    raw.close_sync().unwrap();
    assert!(matches!(
        MaterialStore::open(&path),
        Err(StoreError::UnsupportedSchema { version: 99 })
    ));
}
