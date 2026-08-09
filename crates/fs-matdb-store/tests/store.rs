//! Conformance battery for the FrankenSQLite material store
//! (bead frankensim-oecdy): discovery correctness, evaluation parity
//! by construction, staleness/tamper fail-closed, and bitwise rebuild.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MatDbError, NormalizedPack, ObservationDataset, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_matdb_store::{MaterialStore, StoreError};
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
    let dir = std::env::temp_dir().join("fs-matdb-store-battery");
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir.join(format!("{tag}-{}.db", std::process::id()))
        .to_string_lossy()
        .into_owned()
}

fn seeded_store(path: &str) -> MaterialStore {
    let _ = std::fs::remove_file(path);
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
    let _ = std::fs::remove_file(&empty_path);
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
    let raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_with_params_sync(
        "UPDATE claims SET scalar_value = ?1 WHERE property = 'young_modulus' \
         AND pack_id = 'steel-304-synth'",
        &[SqliteValue::Float(999.0)],
    )
    .expect("tamper");
    drop(raw);
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
    // (c) Tampering with the CANONICAL BYTES is caught by the content
    // hash at decode time.
    let raw = AsyncConnection::open_sync(&path).expect("raw");
    raw.execute_sync("UPDATE packs SET bytes = x'deadbeef' WHERE pack_id = 'steel-304-synth'")
        .expect("corrupt");
    drop(raw);
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
fn duplicates_refuse() {
    let path = scratch_path("dup");
    let store = seeded_store(&path);
    assert!(matches!(
        store.ingest_pack(&test_pack("balsa-synth", &[("density", 160.0)])),
        Err(StoreError::DuplicatePack { .. })
    ));
}
