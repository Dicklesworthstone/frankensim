//! fs-ledger conformance suite (plan §13.3): any reimplementation must pass.
//!
//! Covers the bead's acceptance criteria: versioned schema migrations,
//! oracle-verified BLAKE3 content addressing, hash dedupe (same bytes → one
//! row), chunked multi-part storage round-trips, concurrent snapshot readers
//! during an appending write sweep, a kill -9 crash-recovery battery (G4
//! class), corruption failing loudly, and an events/sec throughput smoke
//! ledgered as its own metric.
//!
//! Every case emits a JSON-lines verdict (seeds and fixture data inline) so
//! failures are reproducible from the log alone (docs/CONVENTIONS.md).

mod common;

use std::sync::{
    Arc, Barrier,
    atomic::{AtomicU32, Ordering},
};

use common::SyncConnection;
use fs_ledger::{
    Blake3, CONTROLLED_ARTIFACT_TILE_LEN, CONTROLLED_OP_FIELD_TILE_LEN, ControlledOpRead,
    ControlledVisibleOpPage, EdgeRole, EventRow, ExecMode, FiveExplicits, Ledger, LedgerError,
    MAIN_BRANCH, OpOutcome, OpVariableField, SCHEMA_VERSION, STORAGE_CHUNK_LEN, hash_bytes,
};

static NEXT_DB: AtomicU32 = AtomicU32::new(0);

fn temp_db(tag: &str) -> String {
    let n = NEXT_DB.fetch_add(1, Ordering::Relaxed);
    std::env::temp_dir()
        .join(format!(
            "fs-ledger-conf-{tag}-{}-{n}.db",
            std::process::id()
        ))
        .display()
        .to_string()
}

/// Best-effort cleanup of a database and its WAL/shm sidecars.
fn cleanup_db(path: &str) {
    for suffix in [
        "",
        "-wal",
        "-shm",
        "-journal",
        ".fsqlite-wal",
        ".fsqlite-shm",
        ".fsqlite-journal",
    ] {
        let _ = std::fs::remove_file(format!("{path}{suffix}"));
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ledger/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

const FX: FiveExplicits<'static> = FiveExplicits {
    seed: &[0x5E, 0xED, 0x00, 0x01],
    versions: r#"{"constellation":"f92683cc4572a198"}"#,
    budget: r#"{"wall_s":30}"#,
    capability: r#"{"ops":["ledger.*"]}"#,
};

/// The official BLAKE3 test input: bytes cycling 0..=250.
fn official_pattern(len: usize) -> Vec<u8> {
    (0..len).map(|i| (i % 251) as u8).collect()
}

/// Official BLAKE3 test vectors (subset spanning sub-block, block-edge,
/// chunk-edge, multi-chunk, and multi-level-tree inputs), cross-generated
/// from an independent oracle implementation.
const BLAKE3_VECTORS: &[(usize, &str)] = &[
    (
        0,
        "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262",
    ),
    (
        1,
        "2d3adedff11b61f14c886e35afa036736dcd87a74d27b5c1510225d0f592e213",
    ),
    (
        2,
        "7b7015bb92cf0b318037702a6cdd81dee41224f734684c2c122cd6359cb1ee63",
    ),
    (
        3,
        "e1be4d7a8ab5560aa4199eea339849ba8e293d55ca0a81006726d184519e647f",
    ),
    (
        63,
        "e9bc37a594daad83be9470df7f7b3798297c3d834ce80ba85d6e207627b7db7b",
    ),
    (
        64,
        "4eed7141ea4a5cd4b788606bd23f46e212af9cacebacdc7d1f4c6dc7f2511b98",
    ),
    (
        65,
        "de1e5fa0be70df6d2be8fffd0e99ceaa8eb6e8c93a63f2d8d1c30ecb6b263dee",
    ),
    (
        127,
        "d81293fda863f008c09e92fc382a81f5a0b4a1251cba1634016a0f86a6bd640d",
    ),
    (
        128,
        "f17e570564b26578c33bb7f44643f539624b05df1a76c81f30acd548c44b45ef",
    ),
    (
        129,
        "683aaae9f3c5ba37eaaf072aed0f9e30bac0865137bae68b1fde4ca2aebdcb12",
    ),
    (
        1023,
        "10108970eeda3eb932baac1428c7a2163b0e924c9a9e25b35bba72b28f70bd11",
    ),
    (
        1024,
        "42214739f095a406f3fc83deb889744ac00df831c10daa55189b5d121c855af7",
    ),
    (
        1025,
        "d00278ae47eb27b34faecf67b4fe263f82d5412916c1ffd97c8cb7fb814b8444",
    ),
    (
        2048,
        "e776b6028c7cd22a4d0ba182a8bf62205d2ef576467e838ed6f2529b85fba24a",
    ),
    (
        2049,
        "5f4d72f40d7a5f82b15ca2b2e44b1de3c2ef86c426c95c1af0b6879522563030",
    ),
    (
        3072,
        "b98cb0ff3623be03326b373de6b9095218513e64f1ee2edd2525c7ad1e5cffd2",
    ),
    (
        3073,
        "7124b49501012f81cc7f11ca069ec9226cecb8a2c850cfe644e327d22d3e1cd3",
    ),
    (
        4096,
        "015094013f57a5277b59d8475c0501042c0b642e531b0a1c8f58d2163229e969",
    ),
    (
        4097,
        "9b4052b38f1c5fc8b1f9ff7ac7b27cd242487b3d890d15c96a1c25b8aa0fb995",
    ),
    (
        5120,
        "9cadc15fed8b5d854562b26a9536d9707cadeda9b143978f319ab34230535833",
    ),
    (
        5121,
        "628bd2cb2004694adaab7bbd778a25df25c47b9d4155a55f8fbd79f2fe154cff",
    ),
    (
        6144,
        "3e2e5b74e048f3add6d21faab3f83aa44d3b2278afb83b80b3c35164ebeca205",
    ),
    (
        6145,
        "f1323a8631446cc50536a9f705ee5cb619424d46887f3c376c695b70e0f0507f",
    ),
    (
        7168,
        "61da957ec2499a95d6b8023e2b0e604ec7f6b50e80a9678b89d2628e99ada77a",
    ),
    (
        7169,
        "a003fc7a51754a9b3c7fae0367ab3d782dccf28855a03d435f8cfe74605e7817",
    ),
    (
        8192,
        "aae792484c8efe4f19e2ca7d371d8c467ffb10748d8a5a1ae579948f718a2a63",
    ),
    (
        8193,
        "bab6c09cb8ce8cf459261398d2e7aef35700bf488116ceb94a36d0f5f1b7bc3b",
    ),
    (
        16384,
        "f875d6646de28985646f34ee13be9a576fd515f76b5b0a26bb324735041ddde4",
    ),
    (
        31744,
        "62b6960e1a44bcc1eb1a611a8d6235b6b4b78f32e7abc4fb4c6cdcce94895c47",
    ),
    (
        102400,
        "bc3e3d41a1146b069abffad3c0d44860cf664390afce4d9661f7902e7943e085",
    ),
    (
        1048576,
        "74cb441fd087764ca9c3694da742ebe30cbeb3060a17009ca81825c7a8d10343",
    ),
    (
        2097153,
        "52dc212cb4cc61cb94d25bd7b1d47b256e4c3a6d68956df50c235c37a2aeacd7",
    ),
];

#[test]
fn ledger_001_blake3_matches_official_vectors() {
    for &(len, expected) in BLAKE3_VECTORS {
        let got = hash_bytes(&official_pattern(len)).to_hex();
        assert_eq!(got, expected, "BLAKE3 mismatch at input length {len}");
    }
    verdict(
        "ledger-001",
        "32 official-pattern BLAKE3 vectors (0 B .. 2 MiB+1) matched",
    );
}

#[test]
fn ledger_002_blake3_streaming_split_property() {
    // Deterministic LCG; the seed is in the log line for replay.
    let seed: u64 = 0x5EED_1ED6_E200_0002;
    let mut x = seed;
    let mut lcg = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x
    };
    for round in 0..64 {
        let len = (lcg() % 20_000) as usize;
        let data: Vec<u8> = (0..len).map(|_| (lcg() >> 33) as u8).collect();
        let whole = hash_bytes(&data);
        let mut hasher = Blake3::new();
        let mut offset = 0;
        while offset < len {
            let piece = 1 + (lcg() % 4096) as usize;
            let end = (offset + piece).min(len);
            hasher.update(&data[offset..end]);
            offset = end;
        }
        assert_eq!(
            hasher.finalize(),
            whole,
            "streaming/one-shot divergence: seed={seed:#x} round={round} len={len}"
        );
    }
    verdict(
        "ledger-002",
        "64 random split patterns agree with one-shot hashing (seeded LCG)",
    );
}

#[test]
fn ledger_003_schema_migration_versioned() {
    let db = temp_db("migrate");
    {
        let l = Ledger::open(&db).expect("fresh open");
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
    }
    {
        // Reopen: idempotent, same version, tables intact.
        let l = Ledger::open(&db).expect("reopen");
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
        assert!(l.lint().unwrap().is_clean());
    }
    {
        // A file from the future is refused, not clobbered.
        let raw = SyncConnection::open(&db).expect("raw open");
        raw.execute(&format!("PRAGMA user_version = {}", SCHEMA_VERSION + 40))
            .unwrap();
        drop(raw);
        match Ledger::open(&db) {
            Err(err) => assert_eq!(err.code(), "LedgerFutureSchema"),
            Ok(_) => panic!("future schema must be refused"),
        }
    }
    cleanup_db(&db);
    verdict(
        "ledger-003",
        "fresh migrate + idempotent reopen + future-version refusal",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // One move/alias/replacement/migration identity scenario.
fn ledger_003b_instance_identity_survives_moves_aliases_and_migration() {
    let first_memory = Ledger::open(":memory:").expect("first memory ledger");
    let first_memory_id = first_memory.instance_id();
    assert_eq!(
        first_memory
            .checked_instance_id()
            .expect("checked identity"),
        first_memory_id
    );
    let uuid = first_memory_id.as_bytes();
    assert_eq!(uuid[6] & 0xf0, 0x40, "identity has UUID v4 shape");
    assert_eq!(uuid[8] & 0xc0, 0x80, "identity has RFC 4122 variant");
    let rendered = first_memory_id.to_uuid_string();
    assert_eq!(rendered.len(), 36);
    assert_eq!(
        rendered
            .bytes()
            .enumerate()
            .filter_map(|(index, byte)| (byte == b'-').then_some(index))
            .collect::<Vec<_>>(),
        vec![8, 13, 18, 23]
    );
    let moved_memory = Box::new(first_memory);
    assert_eq!(
        moved_memory.instance_id(),
        first_memory_id,
        "moving one handle cannot change its sink authority"
    );
    assert_ne!(
        moved_memory.instance_id(),
        Ledger::open(":memory:")
            .expect("independent memory ledger")
            .instance_id(),
        "independent memory databases must never alias by address reuse"
    );

    let db = temp_db("instance-id");
    let original_id = {
        let ledger = Ledger::open(&db).expect("create persistent ledger");
        ledger.instance_id()
    };
    assert_eq!(
        Ledger::open(&db).expect("reopen").instance_id(),
        original_id,
        "reopening the same file preserves identity"
    );
    let path = std::path::Path::new(&db);
    let alias = format!(
        "{}/./{}",
        path.parent().expect("temporary parent").display(),
        path.file_name()
            .expect("temporary filename")
            .to_string_lossy()
    );
    assert_eq!(
        Ledger::open(&alias).expect("path alias").instance_id(),
        original_id,
        "path spelling is not physical sink identity"
    );

    let archived = format!("{db}.archived");
    std::fs::rename(&db, &archived).expect("archive the original database without deleting it");
    // A SQLite WAL is part of the physical database, not disposable path
    // decoration. Move any surviving sidecars with the archived main file so
    // the original path truly denotes a fresh physical database.
    for suffix in ["-wal", "-shm"] {
        let sidecar = format!("{db}{suffix}");
        if std::path::Path::new(&sidecar).exists() {
            std::fs::rename(&sidecar, format!("{archived}{suffix}"))
                .expect("archive the original database sidecar");
        }
    }
    let replacement_id = Ledger::open(&db)
        .expect("replacement database at the same path")
        .instance_id();
    assert_ne!(
        replacement_id, original_id,
        "a replacement database cannot inherit path-based authority"
    );
    cleanup_db(&db);
    cleanup_db(&archived);

    let old = temp_db("instance-id-v3");
    {
        let raw = SyncConnection::open(&old).expect("raw v3 ledger");
        for batch in [
            fs_ledger::schema::V1,
            fs_ledger::schema::V2,
            fs_ledger::schema::V3,
        ] {
            for ddl in batch {
                raw.execute(ddl).expect("construct shipped v3 schema");
            }
        }
        raw.execute("PRAGMA user_version = 3")
            .expect("mark v3 schema");
    }
    let migrated = Ledger::open(&old).expect("v3 identity migration");
    assert_eq!(migrated.schema_version().unwrap(), SCHEMA_VERSION);
    assert_eq!(migrated.table_count("ledger_identity").unwrap(), 1);
    let migrated_id = migrated.instance_id();
    drop(migrated);
    assert_eq!(
        Ledger::open(&old)
            .expect("reopen migrated v3")
            .instance_id(),
        migrated_id,
        "migration seeds identity atomically and only once"
    );
    cleanup_db(&old);

    let old_v4 = temp_db("instance-id-v4");
    let expected_v4_id = [
        0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x46, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee,
        0xff,
    ];
    {
        let raw = SyncConnection::open(&old_v4).expect("raw v4 ledger");
        for batch in [
            fs_ledger::schema::V1,
            fs_ledger::schema::V2,
            fs_ledger::schema::V3,
            fs_ledger::schema::V4,
        ] {
            for ddl in batch {
                raw.execute(ddl).expect("construct shipped v4 schema");
            }
        }
        raw.execute(
            "INSERT INTO ledger_identity(singleton, instance_id) VALUES \
             (1, X'00112233445546778899AABBCCDDEEFF')",
        )
        .expect("seed v4 identity");
        raw.execute("PRAGMA user_version = 4")
            .expect("mark v4 schema");
    }
    let migrated_v4 = Ledger::open(&old_v4).expect("v4 guard migration");
    assert_eq!(migrated_v4.schema_version().unwrap(), SCHEMA_VERSION);
    assert_eq!(migrated_v4.instance_id().as_bytes(), expected_v4_id);
    let second_v5_handle = Ledger::open(&old_v4).expect("second v5 handle");
    assert_eq!(
        second_v5_handle
            .checked_instance_id()
            .expect("second checked identity"),
        migrated_v4.instance_id()
    );
    {
        let raw = SyncConnection::open(&old_v4).expect("raw immutable identity check");
        assert!(
            raw.execute(
                "UPDATE ledger_identity SET instance_id = \
                 X'102132435465467788A9BACBDCEDFE0F' WHERE singleton = 1",
            )
            .is_err(),
            "v5 must refuse even a valid UUID-shaped identity replacement"
        );
        assert!(
            raw.execute("DELETE FROM ledger_identity WHERE singleton = 1")
                .is_err(),
            "v5 must refuse identity deletion"
        );
        assert!(
            raw.execute(
                "INSERT OR REPLACE INTO ledger_identity(singleton, instance_id) VALUES \
                 (1, X'102132435465467788A9BACBDCEDFE0F')",
            )
            .is_err(),
            "v5 must refuse replacement through the insert conflict path"
        );
    }
    assert_eq!(
        migrated_v4
            .checked_instance_id()
            .expect("identity after refused mutations")
            .as_bytes(),
        expected_v4_id
    );
    drop(second_v5_handle);
    drop(migrated_v4);
    cleanup_db(&old_v4);
    verdict(
        "ledger-003b",
        "opaque identity survives moves, aliases, v3/v4 migration, and v5 mutation guards; replacement rotates",
    );
}

#[test]
fn ledger_003c_missing_persisted_identity_fails_closed() {
    let db = temp_db("missing-instance-id");
    {
        let raw = SyncConnection::open(&db).expect("raw missing-identity v4 ledger");
        for batch in [
            fs_ledger::schema::V1,
            fs_ledger::schema::V2,
            fs_ledger::schema::V3,
            fs_ledger::schema::V4,
        ] {
            for ddl in batch {
                raw.execute(ddl).expect("construct shipped v4 schema");
            }
        }
        raw.execute("PRAGMA user_version = 4")
            .expect("mark missing-identity v4 schema");
    }
    assert!(matches!(
        Ledger::open(&db),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));
    {
        let raw = SyncConnection::open(&db).expect("inspect refused v4 migration");
        assert_eq!(
            raw.query_row("PRAGMA user_version")
                .expect("version")
                .get(0),
            Some(&fsqlite::SqliteValue::Integer(4))
        );
    }
    cleanup_db(&db);

    let malformed = temp_db("malformed-instance-id");
    {
        let raw = SyncConnection::open(&malformed).expect("raw malformed v4 ledger");
        for batch in [
            fs_ledger::schema::V1,
            fs_ledger::schema::V2,
            fs_ledger::schema::V3,
            fs_ledger::schema::V4,
        ] {
            for ddl in batch {
                raw.execute(ddl).expect("construct shipped v4 schema");
            }
        }
        raw.execute(
            "INSERT INTO ledger_identity(singleton, instance_id) VALUES \
             (1, X'00000000000000000000000000000000')",
        )
        .expect("install invalid UUID bits");
        raw.execute("PRAGMA user_version = 4")
            .expect("mark malformed v4 schema");
    }
    assert!(matches!(
        Ledger::open(&malformed),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));
    cleanup_db(&malformed);

    let premature = temp_db("premature-malformed-instance-id");
    {
        let raw = SyncConnection::open(&premature).expect("raw premature v4 ledger");
        for batch in [
            fs_ledger::schema::V1,
            fs_ledger::schema::V2,
            fs_ledger::schema::V3,
            fs_ledger::schema::V4,
        ] {
            for ddl in batch {
                raw.execute(ddl).expect("construct premature v4 schema");
            }
        }
        raw.execute(
            "INSERT INTO ledger_identity(singleton, instance_id) VALUES \
             (1, X'00000000000000000000000000000000')",
        )
        .expect("install premature invalid UUID bits");
        raw.execute("PRAGMA user_version = 3")
            .expect("retain the old marker");
    }
    assert!(matches!(
        Ledger::open(&premature),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));
    {
        let raw = SyncConnection::open(&premature).expect("reopen refused migration");
        assert_eq!(
            raw.query_row("PRAGMA user_version")
                .expect("version")
                .get(0),
            Some(&fsqlite::SqliteValue::Integer(3)),
            "malformed premature identity must not receive the v4 marker"
        );
    }
    cleanup_db(&premature);

    verdict(
        "ledger-003c",
        "v4 and premature-v4 schemas missing or carrying malformed identity refuse before version advancement",
    );
}

#[test]
fn ledger_003d_checked_identity_detects_stale_open_handles() {
    let db = temp_db("stale-instance-id-handle");
    let old = Ledger::open(&db).expect("old handle");
    let peer = Ledger::open(&db).expect("peer old handle");
    let original = old.instance_id();
    assert_eq!(peer.instance_id(), original);

    {
        let raw = SyncConnection::open(&db).expect("raw DDL bypass fixture");
        raw.execute("DROP TRIGGER trg_ledger_identity_immutable_update")
            .expect("remove update guard for out-of-band fixture");
        assert!(matches!(
            Ledger::open(&db),
            Err(LedgerError::SchemaMismatch { .. })
        ));
        raw.execute(
            "UPDATE ledger_identity SET instance_id = \
             X'102132435465467788A9BACBDCEDFE0F' WHERE singleton = 1",
        )
        .expect("install valid replacement identity out of band");
        raw.execute(fs_ledger::schema::V5[0])
            .expect("restore exact shipped update guard");
    }

    assert!(matches!(
        old.checked_instance_id(),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));
    assert!(matches!(
        peer.checked_instance_id(),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));
    assert!(matches!(
        old.lint(),
        Err(LedgerError::InstanceIdentityCorrupt { .. })
    ));

    let current = Ledger::open(&db).expect("new handle after out-of-band replacement");
    assert_ne!(current.instance_id(), original);
    assert_eq!(
        current
            .checked_instance_id()
            .expect("new handle agrees with current row"),
        current.instance_id()
    );
    drop(current);
    drop(peer);
    drop(old);
    cleanup_db(&db);
    verdict(
        "ledger-003d",
        "checked identity and lint reject old handles after an out-of-band valid UUID replacement",
    );
}

#[test]
fn ledger_003e_checked_identity_re_attests_exact_mutation_guards() {
    let db = temp_db("identity-guard-attestation");
    let ledger = Ledger::open(&db).expect("guard-attested ledger");
    let original = ledger.instance_id();
    let guards = [
        "trg_ledger_identity_immutable_update",
        "trg_ledger_identity_immutable_delete",
        "trg_ledger_identity_immutable_reinsert",
    ];
    let raw = SyncConnection::open(&db).expect("raw identity-guard fixture");

    for (index, name) in guards.into_iter().enumerate() {
        raw.execute(&format!("DROP TRIGGER {name}"))
            .expect("drop one identity guard without changing the row");
        assert!(matches!(
            ledger.checked_instance_id(),
            Err(LedgerError::SchemaMismatch { .. })
        ));
        if index == 0 {
            assert!(matches!(
                ledger.lint(),
                Err(LedgerError::SchemaMismatch { .. })
            ));
        }
        raw.execute(fs_ledger::schema::V5[index])
            .expect("restore exact shipped identity guard");
        assert_eq!(
            ledger
                .checked_instance_id()
                .expect("restored guard and unchanged row attest"),
            original
        );
    }

    raw.execute("DROP TRIGGER trg_ledger_identity_immutable_update")
        .expect("drop update guard for changed-definition fixture");
    raw.execute(
        "CREATE TRIGGER trg_ledger_identity_immutable_update \
         BEFORE UPDATE ON ledger_identity WHEN 0 \
         BEGIN \
           SELECT RAISE(ABORT, 'ledger_identity is immutable'); \
         END",
    )
    .expect("install same-name weakened identity guard");
    assert!(matches!(
        ledger.checked_instance_id(),
        Err(LedgerError::SchemaMismatch { .. })
    ));
    raw.execute("DROP TRIGGER trg_ledger_identity_immutable_update")
        .expect("drop weakened identity guard");
    raw.execute(fs_ledger::schema::V5[0])
        .expect("restore exact update guard");
    assert_eq!(
        ledger
            .checked_instance_id()
            .expect("exact guards restore checked authority"),
        original
    );

    drop(raw);
    drop(ledger);
    cleanup_db(&db);
    verdict(
        "ledger-003e",
        "checked identity and lint re-attest all three exact mutation guards before trusting an unchanged row",
    );
}

#[test]
fn ledger_004_dedupe_and_chunked_round_trip() {
    let db = temp_db("chunk");
    let mut l = Ledger::open(&db).expect("open");
    // Strictly larger than one storage chunk → chunk rows.
    let big = official_pattern(STORAGE_CHUNK_LEN + STORAGE_CHUNK_LEN / 2 + 17);
    let r1 = l
        .put_artifact("field", &big, Some(r#"{"units":"Pa"}"#))
        .unwrap();
    assert!(r1.chunked && !r1.deduped);
    assert_eq!(l.table_count("artifact_chunks").unwrap(), 2);
    // Identical bytes through BOTH write paths dedupe to the same single row.
    let r2 = l.put_artifact("field", &big, None).unwrap();
    assert!(r2.deduped);
    let mut w = l.artifact_writer("field").unwrap();
    for piece in big.chunks(65_537) {
        w.write(piece).unwrap();
    }
    let r3 = w.finish(None).unwrap();
    assert!(r3.deduped);
    assert_eq!(r1.hash, r3.hash);
    assert_eq!(l.table_count("artifacts").unwrap(), 1);
    // Round trips: whole-value and streamed reads reproduce the bytes.
    assert_eq!(l.get_artifact(&r1.hash).unwrap().unwrap(), big);
    let mut streamed = Vec::new();
    let n = l
        .read_artifact_chunks(&r1.hash, &mut |c| streamed.extend_from_slice(c))
        .unwrap();
    assert_eq!(n, Some(big.len() as u64));
    assert_eq!(streamed, big);
    assert!(l.lint().unwrap().is_clean());
    assert!(l.verify_artifact_integrity().unwrap().is_clean());
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-004",
        "6 MiB artifact: chunked storage, dual-path dedupe, byte round trip",
    );
}

#[test]
fn ledger_004b_controlled_bounded_stream_stops_at_callback_boundary() {
    let db = temp_db("controlled-stream");
    let l = Ledger::open(&db).expect("open");

    let inline = l
        .put_artifact("bounded-inline", b"inline evidence", None)
        .expect("inline artifact");
    let mut inline_callbacks = 0usize;
    let inline_result = l
        .read_artifact_chunks_bounded_controlled(&inline.hash, inline.len, &mut |chunk| {
            inline_callbacks += 1;
            assert_eq!(chunk, b"inline evidence");
            std::ops::ControlFlow::Break("planned-inline-stop")
        })
        .expect("controlled inline read")
        .expect("inline artifact exists");
    assert_eq!(
        inline_result,
        std::ops::ControlFlow::Break("planned-inline-stop")
    );
    assert_eq!(inline_callbacks, 1);

    let exact_cap_bytes = official_pattern(STORAGE_CHUNK_LEN);
    let exact_cap = l
        .put_artifact("bounded-inline-exact-cap", &exact_cap_bytes, None)
        .expect("exact-cap inline artifact");
    assert!(!exact_cap.chunked, "the exact storage bound remains inline");
    let mut exact_cap_callbacks = 0usize;
    let exact_cap_stop = l
        .read_artifact_chunks_bounded_controlled(&exact_cap.hash, exact_cap.len, &mut |tile| {
            exact_cap_callbacks += 1;
            assert_eq!(tile.len(), CONTROLLED_ARTIFACT_TILE_LEN);
            assert_eq!(tile, &exact_cap_bytes[..CONTROLLED_ARTIFACT_TILE_LEN]);
            std::ops::ControlFlow::Break("planned-first-tile-stop")
        })
        .expect("controlled exact-cap read")
        .expect("exact-cap artifact exists");
    assert_eq!(
        exact_cap_stop,
        std::ops::ControlFlow::Break("planned-first-tile-stop")
    );
    assert_eq!(
        exact_cap_callbacks, 1,
        "an exact-cap inline row stops after its first controlled tile"
    );

    let chunked_bytes = official_pattern(STORAGE_CHUNK_LEN * 2 + 17);
    let chunked = l
        .put_artifact("bounded-chunked", &chunked_bytes, None)
        .expect("chunked artifact");
    assert!(chunked.chunked);
    let mut stopped_callbacks = 0usize;
    let stopped = l
        .read_artifact_chunks_bounded_controlled(&chunked.hash, chunked.len, &mut |chunk| {
            stopped_callbacks += 1;
            assert_eq!(chunk.len(), CONTROLLED_ARTIFACT_TILE_LEN);
            std::ops::ControlFlow::Break("planned-first-tile-stop")
        })
        .expect("controlled chunked read")
        .expect("chunked artifact exists");
    assert_eq!(
        stopped,
        std::ops::ControlFlow::Break("planned-first-tile-stop")
    );
    assert_eq!(
        stopped_callbacks, 1,
        "a callback break must not deliver a later controlled tile"
    );

    let mut over_cap_callbacks = 0usize;
    let over_cap =
        l.read_artifact_chunks_bounded_controlled(&chunked.hash, chunked.len - 1, &mut |_| {
            over_cap_callbacks += 1;
            std::ops::ControlFlow::<()>::Continue(())
        });
    assert!(matches!(
        over_cap,
        Err(fs_ledger::LedgerError::ArtifactReadLimit {
            limit,
            observed,
            ..
        }) if limit == chunked.len - 1 && observed == chunked.len
    ));
    assert_eq!(
        over_cap_callbacks, 0,
        "metadata preflight must reject an over-cap artifact before delivering bytes"
    );

    let mut completed_bytes = Vec::new();
    let mut completed_callbacks = 0usize;
    let completed = l
        .read_artifact_chunks_bounded_controlled(&chunked.hash, chunked.len, &mut |chunk| {
            completed_callbacks += 1;
            assert!(chunk.len() <= CONTROLLED_ARTIFACT_TILE_LEN);
            completed_bytes.extend_from_slice(chunk);
            std::ops::ControlFlow::<()>::Continue(())
        })
        .expect("complete controlled read")
        .expect("chunked artifact exists");
    assert_eq!(
        completed,
        std::ops::ControlFlow::Continue(chunked.len),
        "only a complete integrity-checked read returns the total length"
    );
    assert_eq!(completed_bytes, chunked_bytes);
    assert_eq!(
        completed_callbacks,
        chunked_bytes.len().div_ceil(CONTROLLED_ARTIFACT_TILE_LEN)
    );

    let mut infallible_row_lengths = Vec::new();
    let infallible_len = l
        .read_artifact_chunks_bounded(&chunked.hash, chunked.len, &mut |row| {
            infallible_row_lengths.push(row.len());
        })
        .expect("infallible bounded read")
        .expect("chunked artifact exists");
    assert_eq!(infallible_len, chunked.len);
    assert_eq!(
        infallible_row_lengths,
        vec![STORAGE_CHUNK_LEN, STORAGE_CHUNK_LEN, 17],
        "the established infallible callback retains storage-row granularity"
    );

    let mismatched = l
        .put_artifact("bounded-complete-mismatch", b"complete then compare", None)
        .expect("mismatch artifact");
    l.corrupt_artifact_for_test(&mismatched.hash)
        .expect("corrupt mismatch artifact");
    let mut mismatch_callbacks = 0usize;
    let mismatch =
        l.read_artifact_chunks_bounded_controlled(&mismatched.hash, mismatched.len, &mut |_| {
            mismatch_callbacks += 1;
            std::ops::ControlFlow::<()>::Continue(())
        });
    assert!(matches!(
        mismatch,
        Err(fs_ledger::LedgerError::Corrupt { .. })
    ));
    assert_eq!(
        mismatch_callbacks, 1,
        "a completing callback cannot suppress the final content-hash comparison"
    );

    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-004b",
        "controlled bounded reads stop after one hash-coupled tile and complete reads retain full integrity authority",
    );
}

#[test]
fn ledger_005_corruption_fails_loudly() {
    let db = temp_db("corrupt");
    let l = Ledger::open(&db).expect("open");
    let small = l
        .put_artifact("blob", b"inline precious bytes", None)
        .unwrap();
    let big = l
        .put_artifact("field", &official_pattern(STORAGE_CHUNK_LEN + 1), None)
        .unwrap();
    assert!(l.verify_artifact_integrity().unwrap().is_clean());
    l.corrupt_artifact_for_test(&small.hash).unwrap();
    l.corrupt_artifact_for_test(&big.hash).unwrap();
    assert!(matches!(
        l.get_artifact(&small.hash),
        Err(fs_ledger::LedgerError::Corrupt { .. })
    ));
    assert!(matches!(
        l.get_artifact(&big.hash),
        Err(fs_ledger::LedgerError::Corrupt { .. })
    ));
    let report = l.verify_artifact_integrity().unwrap();
    assert_eq!(report.checked, 2);
    assert!(report.corrupted.contains(&small.hash.to_hex()));
    assert!(report.corrupted.contains(&big.hash.to_hex()));
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-005",
        "ordinary reads and full scan reject inline and chunked corruption",
    );
}

#[test]
fn ledger_006_concurrent_snapshot_readers_during_sweep() {
    const UNITS: usize = 120;
    const READERS: usize = 3;
    let db = temp_db("stress");
    // Create the schema before anyone races.
    drop(Ledger::open(&db).expect("schema open"));

    let db_ref = &db;
    std::thread::scope(|scope| {
        let writer = scope.spawn(move || {
            let l = Ledger::open(db_ref).expect("writer open");
            for i in 0..UNITS {
                // Retry the whole atomic unit on contention.
                let mut attempts = 0;
                loop {
                    match write_unit(&l, i as u64) {
                        Ok(()) => break,
                        Err(e) if e.code() == "LedgerBusy" && attempts < 50 => {
                            attempts += 1;
                            let _ = l.rollback();
                            std::thread::sleep(std::time::Duration::from_millis(2));
                        }
                        Err(e) => panic!("writer unit {i}: {e}"),
                    }
                }
            }
        });
        let mut readers = Vec::new();
        for r in 0..READERS {
            readers.push(scope.spawn(move || {
                let l = Ledger::open(db_ref).expect("reader open");
                let mut last_ops = 0u64;
                let deadline = std::time::Instant::now() + std::time::Duration::from_secs(60);
                loop {
                    assert!(std::time::Instant::now() < deadline, "reader {r} timed out");
                    // One snapshot: counts + hygiene must be internally
                    // consistent no matter what the writer is doing.
                    match snapshot_check(&l) {
                        Ok((ops, edges)) => {
                            assert!(ops >= last_ops, "reader {r}: ops count went backwards");
                            assert_eq!(ops, edges, "reader {r}: edges out of sync inside snapshot");
                            last_ops = ops;
                            if ops >= UNITS as u64 {
                                break;
                            }
                        }
                        Err(e) => {
                            assert_eq!(e.code(), "LedgerBusy", "reader {r}: {e}");
                            let _ = l.rollback();
                        }
                    }
                    std::thread::sleep(std::time::Duration::from_millis(1));
                }
                last_ops
            }));
        }
        writer.join().expect("writer thread");
        for reader in readers {
            assert_eq!(reader.join().expect("reader thread"), UNITS as u64);
        }
    });

    let l = Ledger::open(&db).expect("final open");
    assert_eq!(l.table_count("ops").unwrap(), UNITS as u64);
    assert_eq!(l.table_count("edges").unwrap(), UNITS as u64);
    assert!(l.lint().unwrap().is_clean());
    assert!(l.verify_artifact_integrity().unwrap().is_clean());
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-006",
        "3 snapshot readers observed monotone, internally consistent state during a 120-op sweep",
    );
}

/// One atomic writer unit: op + artifact + edge + metric + event.
fn write_unit(l: &Ledger, i: u64) -> Result<(), fs_ledger::LedgerError> {
    let ti = i64::try_from(i).expect("test counter fits i64");
    l.begin()?;
    let result = (|| {
        let op = l.begin_op(
            Some(b"stress".as_slice()),
            &format!("{{\"unit\":{i}}}"),
            &FX,
            ti,
        )?;
        let bytes: Vec<u8> = (0..2048u64)
            .map(|j| ((i * 31 + j * 7) % 251) as u8)
            .collect();
        let receipt = l.put_artifact("stress-blob", &bytes, None)?;
        l.link(op, &receipt.hash, EdgeRole::Out)?;
        l.record_metric(op, 0, "unit", i as f64)?;
        l.append_event(&EventRow {
            session: Some(b"stress".as_slice()),
            t: ti,
            kind: "tile_complete",
            payload: None,
        })?;
        l.finish_op(op, OpOutcome::Ok, None, ti + 1)?;
        Ok(())
    })();
    match result {
        Ok(()) => l.commit(),
        Err(e) => {
            let _ = l.rollback();
            Err(e)
        }
    }
}

/// Read (ops, edges) counts plus hygiene inside one snapshot transaction.
fn snapshot_check(l: &Ledger) -> Result<(u64, u64), fs_ledger::LedgerError> {
    l.begin()?;
    let result = (|| {
        let ops = l.table_count("ops")?;
        let edges = l.table_count("edges")?;
        let lint = l.lint()?;
        assert!(lint.is_clean(), "snapshot saw dirty state: {lint:?}");
        Ok((ops, edges))
    })();
    match result {
        Ok(v) => {
            l.commit()?;
            Ok(v)
        }
        Err(e) => {
            let _ = l.rollback();
            Err(e)
        }
    }
}

/// Child-process entry for the kill -9 battery: loops atomic writer units
/// until killed. A no-op unless the parent set `FS_LEDGER_CRASH_DB`.
#[test]
fn ledger_007_crash_child_writer() {
    let Ok(db) = std::env::var("FS_LEDGER_CRASH_DB") else {
        return; // not a child invocation
    };
    let l = Ledger::open(&db).expect("child open");
    let mut i = 0u64;
    loop {
        write_unit(&l, i).expect("child write unit");
        i += 1;
    }
}

#[test]
fn ledger_007_crash_kill9_battery() {
    let exe = std::env::current_exe().expect("current exe");
    let seed: u64 = 0x5EED_C4A5_4000_0007;
    let mut x = seed;
    let mut lcg = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x
    };
    for round in 0..6 {
        let db = temp_db(&format!("crash{round}"));
        // Schema exists before the child starts, so the kill always lands in
        // write traffic rather than migration.
        drop(Ledger::open(&db).expect("pre-create"));
        let kill_ms = 40 + (lcg() % 260);
        let mut child = std::process::Command::new(&exe)
            .args(["--exact", "ledger_007_crash_child_writer", "--nocapture"])
            .env("FS_LEDGER_CRASH_DB", &db)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()
            .expect("spawn crash child");
        std::thread::sleep(std::time::Duration::from_millis(kill_ms));
        child.kill().expect("SIGKILL child");
        let _ = child.wait();

        // Recovery: reopen replays the WAL; state must be exactly the
        // committed prefix — complete ops, zero orphans, intact bytes.
        let l = Ledger::open(&db).expect("recovery open");
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
        let lint = l.lint().unwrap();
        assert!(
            lint.is_clean(),
            "round {round}: post-recovery lint dirty: {lint:?}"
        );
        let integrity = l.verify_artifact_integrity().unwrap();
        assert!(
            integrity.is_clean(),
            "round {round}: corrupted artifacts after recovery: {:?}",
            integrity.corrupted
        );
        let ops = l.table_count("ops").unwrap();
        let edges = l.table_count("edges").unwrap();
        let events = l.table_count("events").unwrap();
        assert_eq!(
            ops, edges,
            "round {round}: partial op groups survived the crash"
        );
        assert_eq!(
            ops, events,
            "round {round}: event stream out of sync with ops"
        );
        println!(
            "{{\"suite\":\"fs-ledger/conformance\",\"case\":\"ledger-007\",\"round\":{round},\
             \"seed\":\"{seed:#x}\",\"kill_ms\":{kill_ms},\"recovered_ops\":{ops},\
             \"lint_clean\":true,\"integrity_clean\":true}}"
        );
        drop(l);
        cleanup_db(&db);
    }
    verdict(
        "ledger-007",
        "6 kill -9 rounds recovered to lint-clean, integrity-clean state",
    );
}

#[test]
fn ledger_008_event_throughput_ledgered() {
    let db = temp_db("bench");
    let l = Ledger::open(&db).expect("open");
    let payload = r#"{"tile":42,"kernel":"bench"}"#;
    let batch: Vec<EventRow<'_>> = (0..1000)
        .map(|i| EventRow {
            session: Some(b"bench".as_slice()),
            t: i,
            kind: "tile_complete",
            payload: Some(payload),
        })
        .collect();
    let start = std::time::Instant::now();
    let mut appended = 0u64;
    while start.elapsed() < std::time::Duration::from_millis(500) {
        l.append_events(&batch).expect("batch append");
        appended += batch.len() as u64;
    }
    let events_per_sec = appended as f64 / start.elapsed().as_secs_f64();
    // Ledger the measurement in the ledger itself (plan §11.2: throughput
    // benchmarks are metrics rows). The op records the benchmark context.
    let op = l
        .begin_op(
            Some(b"bench".as_slice()),
            r#"{"op":"bench.append_events"}"#,
            &FX,
            0,
        )
        .expect("bench op");
    l.record_metric(op, 0, "events_per_sec", events_per_sec)
        .expect("metric");
    l.finish_op(op, OpOutcome::Ok, None, 1).expect("finish");
    // De-minimis floor: proves the batched path works at all under debug
    // builds + parallel test I/O. NOT a roofline claim (those need machine
    // fingerprints and acceptance bands; plan §14) — the real number lives
    // in the metrics row this test just wrote.
    assert!(
        events_per_sec > 100.0,
        "throughput smoke floor: {events_per_sec:.0}/s"
    );
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-008",
        &format!("sustained {events_per_sec:.0} events/s (batched, fsync-per-commit), ledgered"),
    );
}

#[test]
fn ledger_011_schema_attestation_gauntlet() {
    // gp3.18: `CREATE TABLE IF NOT EXISTS` must never launder an alien,
    // partial, or mangled schema into a version-stamped one. Six probes.

    // (1) valid-empty-database: a fresh file initializes ATOMICALLY to the
    // current version and re-attests clean on reopen.
    let db = temp_db("attest-empty");
    {
        let l = Ledger::open(&db).expect("fresh empty file initializes");
        assert_eq!(l.schema_version().unwrap(), SCHEMA_VERSION);
    }
    Ledger::open(&db).expect("reopen attests the schema it just wrote");
    cleanup_db(&db);

    // (2) conflicting-object at v0: ANY pre-existing user object in an
    // unversioned file refuses initialization (fail closed, file intact).
    let db = temp_db("attest-conflict");
    {
        let raw = SyncConnection::open(&db).expect("raw");
        raw.execute("CREATE TABLE ops(x TEXT)").expect("alien ops");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("alien object at v0 must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    assert!(
        err.to_string().contains("pre-existing table `ops`"),
        "{err}"
    );
    {
        // Fail closed means UNTOUCHED: still v0, alien table intact.
        let raw = SyncConnection::open(&db).expect("raw reopen");
        assert_eq!(
            raw.query_row("PRAGMA user_version").expect("v").get(0),
            Some(&fsqlite::SqliteValue::Integer(0))
        );
    }
    cleanup_db(&db);

    // (2b) an internal-looking but legal user name is still inventoried.
    // SQL LIKE treats `_` as a wildcard, so `NOT LIKE 'sqlite_%'` would hide
    // this table even though only the literal `sqlite_` prefix is reserved.
    let db = temp_db("attest-internal-lookalike");
    {
        let raw = SyncConnection::open(&db).expect("raw");
        raw.execute("CREATE TABLE sqlitex_hidden(x TEXT)")
            .expect("legal user table");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("an internal-looking user object at v0 must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    assert!(
        err.to_string()
            .contains("pre-existing table `sqlitex_hidden`"),
        "{err}"
    );
    cleanup_db(&db);

    // (3) partial-schema: half of v1 committed with the marker forced to 1
    // refuses with the missing objects named.
    let db = temp_db("attest-partial");
    {
        let raw = SyncConnection::open(&db).expect("raw");
        for ddl in fs_ledger::schema::V1.iter().take(4) {
            raw.execute(ddl).expect("partial v1");
        }
        raw.execute("PRAGMA user_version = 1").expect("marker");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("partial schema must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    assert!(err.to_string().contains("missing table"), "{err}");
    cleanup_db(&db);

    verdict(
        "ledger-011",
        "schema attestation: empty-init atomic + conflict/partial refused fail-closed",
    );
}

#[test]
fn ledger_011b_schema_attestation_mangled_definitions() {
    // gp3.18 continued: same-name objects with divergent definitions.

    // (4) wrong-column: a v1 table rebuilt with a different column set
    // (same name) is not the shipped definition.
    let db = temp_db("attest-wrongcol");
    {
        let raw = SyncConnection::open(&db).expect("raw");
        for ddl in fs_ledger::schema::V1 {
            if !ddl.contains("EXISTS metrics") {
                raw.execute(ddl).expect("v1 ddl");
            }
        }
        raw.execute("CREATE TABLE metrics(id INTEGER PRIMARY KEY) STRICT")
            .expect("wrong metrics");
        raw.execute("PRAGMA user_version = 1").expect("marker");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("wrong column set must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    assert!(
        err.to_string().contains("table `metrics` differs") || err.to_string().contains("column"),
        "{err}"
    );
    cleanup_db(&db);

    // (5) wrong-affinity: same column names, one declared type changed.
    let db = temp_db("attest-affinity");
    {
        let raw = SyncConnection::open(&db).expect("raw");
        for ddl in fs_ledger::schema::V1 {
            if ddl.contains("EXISTS edges") {
                raw.execute(&ddl.replace("op INTEGER NOT NULL", "op TEXT NOT NULL"))
                    .expect("mangled edges");
            } else {
                raw.execute(ddl).expect("v1 ddl");
            }
        }
        raw.execute("PRAGMA user_version = 1").expect("marker");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("wrong affinity must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    cleanup_db(&db);

    // (6) missing-index: a fully migrated file whose index was dropped is
    // not the schema its version claims.
    let db = temp_db("attest-noindex");
    {
        let l = Ledger::open(&db).expect("initialize");
        drop(l);
        let raw = SyncConnection::open(&db).expect("raw");
        raw.execute("DROP INDEX idx_ops_session")
            .expect("drop index");
    }
    let Err(err) = Ledger::open(&db) else {
        panic!("missing index must refuse")
    };
    assert_eq!(err.code(), "LedgerSchemaMismatch");
    assert!(
        err.to_string().contains("missing index `idx_ops_session`"),
        "{err}"
    );
    cleanup_db(&db);

    verdict(
        "ledger-011b",
        "schema attestation: wrong-column/wrong-affinity/missing-index refused",
    );
}

#[test]
fn ledger_012_artifact_envelope_conflicts() {
    // gp3.19: identical bytes dedupe ONLY under an agreeing envelope;
    // conflicting kind/meta refuses structurally instead of silently
    // retaining whichever envelope arrived first.
    let bytes = b"envelope gauntlet payload".as_slice();

    // Recommit: the exact same envelope dedupes idempotently.
    let db = temp_db("env-recommit");
    let mut l = Ledger::open(&db).expect("open");
    let first = l
        .put_artifact("field", bytes, Some(r#"{"units":"Pa"}"#))
        .expect("first put");
    assert!(!first.deduped);
    let again = l
        .put_artifact("field", bytes, Some(r#"{"units":"Pa"}"#))
        .expect("recommit with the identical envelope");
    assert!(again.deduped);
    // Whitespace-variant metadata is canonically equal through the engine.
    assert!(
        l.put_artifact("field", bytes, Some(r#"{ "units" : "Pa" }"#))
            .expect("canonically equal meta dedupes")
            .deduped
    );
    // No-claim metadata accepts the stored envelope (streaming contract).
    assert!(
        l.put_artifact("field", bytes, None)
            .expect("no-claim meta dedupes")
            .deduped
    );
    // Conflicting kind refuses.
    let Err(err) = l.put_artifact("mesh", bytes, Some(r#"{"units":"Pa"}"#)) else {
        panic!("conflicting kind must refuse")
    };
    assert_eq!(err.code(), "LedgerArtifactEnvelopeConflict");
    assert!(err.to_string().contains("kind=field"), "{err}");
    // Conflicting metadata refuses and names both sides.
    let Err(err) = l.put_artifact("field", bytes, Some(r#"{"units":"kPa"}"#)) else {
        panic!("conflicting meta must refuse")
    };
    assert_eq!(err.code(), "LedgerArtifactEnvelopeConflict");
    assert!(err.to_string().contains("kPa"), "{err}");
    // A claim against a row stored WITHOUT metadata also refuses.
    let bare = l
        .put_artifact("bare", b"no-meta bytes", None)
        .expect("bare");
    assert!(!bare.deduped);
    assert!(matches!(
        l.put_artifact("bare", b"no-meta bytes", Some("{}")),
        Err(LedgerError::ArtifactEnvelopeConflict { field: "meta", .. })
    ));
    // The streaming writer's dedupe path enforces the same gate.
    let mut w = l.artifact_writer("mesh").expect("writer");
    w.write(bytes).expect("stage");
    assert!(matches!(
        w.finish(None),
        Err(LedgerError::ArtifactEnvelopeConflict { field: "kind", .. })
    ));
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-012",
        "artifact envelopes: agree-or-refuse at every dedupe site",
    );
}

#[test]
fn ledger_012b_envelope_order_concurrency_replay() {
    let bytes = b"envelope gauntlet payload".as_slice();
    // Reversed-order: either arrival order yields ONE stored envelope and
    // one structured refusal — never a silent swallow.
    for (first_kind, second_kind) in [("field", "mesh"), ("mesh", "field")] {
        let db = temp_db("env-order");
        let l = Ledger::open(&db).expect("open");
        l.put_artifact(first_kind, bytes, None).expect("winner");
        let refused = l.put_artifact(second_kind, bytes, None);
        assert!(
            matches!(
                &refused,
                Err(LedgerError::ArtifactEnvelopeConflict { field: "kind", stored, .. })
                    if stored == first_kind
            ),
            "order {first_kind}->{second_kind}: {refused:?}"
        );
        drop(l);
        cleanup_db(&db);
    }

    // Concurrent handles: a second connection deduping through the
    // duplicate-key race path still attests the envelope.
    let db = temp_db("env-concurrent");
    let l1 = Ledger::open(&db).expect("open 1");
    let l2 = Ledger::open(&db).expect("open 2");
    l1.put_artifact("field", bytes, None).expect("writer one");
    assert!(
        l2.put_artifact("field", bytes, None)
            .expect("agreeing concurrent put dedupes")
            .deduped
    );
    assert!(matches!(
        l2.put_artifact("mesh", bytes, None),
        Err(LedgerError::ArtifactEnvelopeConflict { field: "kind", .. })
    ));
    drop((l1, l2));
    cleanup_db(&db);

    // Replay + migration: a file written before this gate (envelope rows
    // unchanged by construction — no schema change) replays the accepted
    // sequence identically after reopen, and the gate binds new puts.
    let db = temp_db("env-replay");
    {
        let l = Ledger::open(&db).expect("open");
        l.put_artifact("field", bytes, Some(r#"{"units":"Pa"}"#))
            .expect("seed");
    }
    {
        let l = Ledger::open(&db).expect("reopen (schema attests)");
        assert!(
            l.put_artifact("field", bytes, Some(r#"{"units":"Pa"}"#))
                .expect("replayed put dedupes identically")
                .deduped
        );
        assert!(matches!(
            l.put_artifact("field", bytes, Some(r#"{"units":"bar"}"#)),
            Err(LedgerError::ArtifactEnvelopeConflict { .. })
        ));
    }
    cleanup_db(&db);

    verdict(
        "ledger-012b",
        "artifact envelopes: order-independent verdicts, concurrent + replay",
    );
}

#[test]
fn ledger_013_file_backed_concurrent_tune_upsert_is_atomic_and_replayable() {
    let db = temp_db("tune-concurrent");
    {
        let seed = Ledger::open(&db).expect("open seed ledger");
        let payload = format!(r#"{{"seed":true,"padding":"{}"}}"#, "x".repeat(1024));
        for index in 0..160 {
            seed.tune_put(
                "concurrent-upsert",
                &format!("seed-{index:04}"),
                b"machine",
                &payload,
                &payload,
            )
            .expect("seed enough rows to split tune leaves");
        }
    }

    let barrier = Arc::new(Barrier::new(2));
    let mut workers = Vec::new();
    for writer in 0..2 {
        let path = db.clone();
        let barrier = Arc::clone(&barrier);
        workers.push(std::thread::spawn(move || {
            let ledger = Ledger::open(&path).expect("open worker ledger");
            barrier.wait();
            for sequence in 0..64 {
                let payload = format!(
                    r#"{{"writer":{writer},"sequence":{sequence},"padding":"{}"}}"#,
                    "y".repeat(1024)
                );
                let mut stored = false;
                for _attempt in 0..200 {
                    match ledger.tune_put(
                        "concurrent-upsert",
                        "shared-shape",
                        b"machine",
                        &payload,
                        &payload,
                    ) {
                        Ok(()) => {
                            stored = true;
                            break;
                        }
                        Err(LedgerError::Busy { .. }) => {
                            std::thread::sleep(std::time::Duration::from_millis(2));
                        }
                        Err(error) => panic!("writer {writer} sequence {sequence}: {error}"),
                    }
                }
                assert!(
                    stored,
                    "writer {writer} sequence {sequence} exhausted retries"
                );
            }
        }));
    }
    for worker in workers {
        worker.join().expect("tune writer thread");
    }

    let final_payload = {
        let ledger = Ledger::open(&db).expect("open verifier");
        assert_eq!(ledger.table_count("tune").expect("count tune rows"), 161);
        let row = ledger
            .tune_get("concurrent-upsert", "shared-shape", b"machine")
            .expect("read final shared row")
            .expect("shared row exists");
        assert_eq!(row.params, row.measured, "upsert fields must never tear");
        assert!(
            row.params.contains(r#""sequence":63"#),
            "each writer's final operation has sequence 63: {}",
            row.params
        );
        assert_eq!(
            ledger
                .tune_rows("concurrent-upsert")
                .expect("bounded scan")
                .len(),
            161
        );
        assert!(ledger.lint().expect("lint after racing upserts").is_clean());
        row.params
    };
    {
        let replay = Ledger::open(&db).expect("reopen tune ledger");
        let row = replay
            .tune_get("concurrent-upsert", "shared-shape", b"machine")
            .expect("replay final row")
            .expect("replayed row exists");
        assert_eq!(row.params, final_payload);
        assert_eq!(row.measured, final_payload);
        assert!(replay.lint().expect("lint replayed ledger").is_clean());
    }
    cleanup_db(&db);
    verdict(
        "ledger-013",
        "file-backed concurrent tune upserts remain single-row, untorn, bounded, and replayable",
    );
}

#[test]
fn ledger_009_version_is_stamped() {
    assert!(!fs_ledger::VERSION.is_empty());
    verdict("ledger-009", "crate version stamped");
}

#[test]
fn ledger_010_nightly_writer_records_run() {
    fn one_json_line<'a>(bytes: &'a [u8], stream: &str) -> &'a str {
        let text = std::str::from_utf8(bytes).expect("nightly_ledger emits UTF-8");
        let line = text
            .strip_suffix('\n')
            .unwrap_or_else(|| panic!("{stream} must end with exactly one newline: {text:?}"));
        assert!(
            !line.chars().any(|ch| ch < ' '),
            "{stream} contains a raw JSON control: {text:?}"
        );
        line
    }

    fn assert_no_db_artifacts(path: &str) {
        for suffix in [
            "",
            "-wal",
            "-shm",
            "-journal",
            ".fsqlite-wal",
            ".fsqlite-shm",
            ".fsqlite-journal",
        ] {
            let candidate = format!("{path}{suffix}");
            assert!(
                !std::path::Path::new(&candidate).exists(),
                "pure nightly admission unexpectedly created {candidate:?}"
            );
        }
    }

    const SUITE: &str = "nightly\",\"forged\":true,\"tail\":\"\\line\nnext\u{0001}";
    const SUITE_JSON: &str =
        "nightly\\\",\\\"forged\\\":true,\\\"tail\\\":\\\"\\\\line\\nnext\\u0001";
    const SHA: &str = "sha\"\\\n\u{0003}";
    const SHA_JSON: &str = "sha\\\"\\\\\\n\\u0003";
    const RUNNER: &str = "runner\"\\\r\t\u{0002}";
    const RUNNER_JSON: &str = "runner\\\"\\\\\\r\\t\\u0002";

    let db = temp_db("nightly");
    let exe = env!("CARGO_BIN_EXE_nightly_ledger");
    let out = std::process::Command::new(exe)
        .args([db.as_str(), "ok", SUITE, "1"])
        .env("GITHUB_SHA", SHA)
        .env("RUNNER_OS", RUNNER)
        .output()
        .expect("run nightly_ledger");
    assert!(
        out.status.success(),
        "nightly_ledger failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out.stderr.is_empty(), "success must not emit stderr");
    let stdout = one_json_line(&out.stdout, "stdout");
    assert!(stdout.contains(&format!("\"metric\":\"{SUITE_JSON}\"")));
    assert!(!stdout.contains("\"forged\":true"));

    let l = Ledger::open(&db).expect("open nightly db");
    assert_eq!(l.table_count("ops").unwrap(), 1);
    assert_eq!(l.table_count("metrics").unwrap(), 1);
    assert_eq!(l.table_count("events").unwrap(), 1);
    let op = l.op(1).unwrap().expect("op row");
    assert_eq!(op.outcome.as_deref(), Some("ok"));
    assert_eq!(
        op.ir,
        format!("{{\"op\":\"ci.nightly\",\"suite\":\"{SUITE_JSON}\"}}")
    );
    assert_eq!(
        op.versions,
        format!(
            "{{\"frankensim\":{{\"availability\":\"reported\",\"source\":\
             \"process_environment\",\"variable\":\"GITHUB_SHA\",\"value\":\"{SHA_JSON}\"}}}}"
        )
    );
    assert_eq!(
        op.capability,
        format!(
            "{{\"ops\":[\"ci.nightly\"],\"runner\":{{\"availability\":\"reported\",\
             \"source\":\"process_environment\",\"variable\":\"RUNNER_OS\",\
             \"value\":\"{RUNNER_JSON}\"}}}}"
        )
    );
    l.append_event(&EventRow {
        session: None,
        t: 1,
        kind: "nightly_writer_stdout_validation",
        payload: Some(stdout),
    })
    .expect("stdout is one valid JSON value");

    assert!(l.lint().unwrap().is_clean());
    // Bad arguments are structured refusals on stderr, never panics.
    const BAD_OUTCOME: &str = "bogus\",\"forged\":true,\"tail\":\"\\\n\u{0002}";
    let bad = std::process::Command::new(exe)
        .args([db.as_str(), BAD_OUTCOME, "x", "y"])
        .output()
        .expect("run nightly_ledger with bad args");
    assert!(!bad.status.success());
    assert!(bad.stdout.is_empty(), "refusal must not emit stdout");
    let bad_stderr = one_json_line(&bad.stderr, "bad-argument stderr");
    assert_eq!(
        bad_stderr,
        format!(
            "{{\"error\":\"NightlyLedger\",\"detail\":\"NightlyInput: field 'outcome' \
             rejected: UTF-8 encoding must contain 1..=5 bytes; observed {}\"}}",
            BAD_OUTCOME.len()
        )
    );
    assert!(!bad_stderr.contains("\"forged\":true"));
    l.append_event(&EventRow {
        session: None,
        t: 2,
        kind: "nightly_writer_refusal_validation",
        payload: Some(bad_stderr),
    })
    .expect("refusal is one valid JSON value");

    // Non-finite values must refuse before opening or partially mutating the
    // ledger; they are not JSON numbers and metrics require finite REALs.
    let nonfinite = std::process::Command::new(exe)
        .args([db.as_str(), "ok", "nonfinite", "NaN"])
        .output()
        .expect("run nightly_ledger with non-finite value");
    assert!(!nonfinite.status.success());
    assert!(nonfinite.stdout.is_empty(), "refusal must not emit stdout");
    let nonfinite_stderr = one_json_line(&nonfinite.stderr, "non-finite stderr");
    assert!(nonfinite_stderr.contains("field 'value' rejected: must be finite"));
    l.append_event(&EventRow {
        session: None,
        t: 3,
        kind: "nightly_writer_nonfinite_validation",
        payload: Some(nonfinite_stderr),
    })
    .expect("non-finite refusal is one valid JSON value");

    // Empty metric identity is a pure admission refusal, before database open.
    let empty_suite = std::process::Command::new(exe)
        .args([db.as_str(), "ok", "", "1"])
        .output()
        .expect("run nightly_ledger with an empty suite");
    assert!(!empty_suite.status.success());
    assert!(
        empty_suite.stdout.is_empty(),
        "refusal must not emit stdout"
    );
    let empty_suite_stderr = one_json_line(&empty_suite.stderr, "empty-suite stderr");
    assert!(empty_suite_stderr.contains(
        "field 'suite' rejected: UTF-8 encoding must contain 1..=65536 bytes; observed 0"
    ));
    l.append_event(&EventRow {
        session: None,
        t: 4,
        kind: "nightly_writer_admission_validation",
        payload: Some(empty_suite_stderr),
    })
    .expect("empty-suite refusal is one valid JSON value");

    // Pure admission must not open a fresh database or create any known
    // journal sidecar. Simultaneously-invalid later fields prove the fixed
    // diagnostic precedence: suite precedes value and provenance.
    let admission_db = temp_db("nightly-admission-refusal");
    cleanup_db(&admission_db);
    assert_no_db_artifacts(&admission_db);
    let admission = std::process::Command::new(exe)
        .args([admission_db.as_str(), "ok", "", "NaN"])
        .env("GITHUB_SHA", "")
        .env("RUNNER_OS", "")
        .output()
        .expect("run pure-admission refusal");
    assert!(!admission.status.success());
    assert!(admission.stdout.is_empty(), "refusal must not emit stdout");
    assert_eq!(
        one_json_line(&admission.stderr, "admission stderr"),
        "{\"error\":\"NightlyLedger\",\"detail\":\"NightlyInput: field 'suite' rejected: UTF-8 encoding must contain 1..=65536 bytes; observed 0\"}"
    );
    assert_no_db_artifacts(&admission_db);

    let cap_db = temp_db("nightly-suite-cap-refusal");
    cleanup_db(&cap_db);
    let oversized_suite = "s".repeat(64 * 1024 + 1);
    let cap = std::process::Command::new(exe)
        .args([cap_db.as_str(), "ok", &oversized_suite, "1"])
        .env_remove("GITHUB_SHA")
        .env_remove("RUNNER_OS")
        .output()
        .expect("run suite cap+1 refusal");
    assert!(!cap.status.success());
    assert!(cap.stdout.is_empty(), "refusal must not emit stdout");
    assert!(one_json_line(&cap.stderr, "cap+1 stderr").contains("observed 65537"));
    assert_no_db_artifacts(&cap_db);

    for (tag, outcome, suite, value, sha, runner, expected_field) in [
        (
            "outcome",
            "bogus",
            "suite",
            "1",
            "sha",
            "runner",
            "field 'outcome'",
        ),
        (
            "value",
            "ok",
            "suite",
            "NaN",
            "sha",
            "runner",
            "field 'value'",
        ),
        (
            "sha",
            "ok",
            "suite",
            "1",
            "",
            "runner",
            "field 'GITHUB_SHA'",
        ),
        ("runner", "ok", "suite", "1", "sha", "", "field 'RUNNER_OS'"),
    ] {
        let pure_db = temp_db(&format!("nightly-pure-{tag}"));
        cleanup_db(&pure_db);
        let refusal = std::process::Command::new(exe)
            .args([pure_db.as_str(), outcome, suite, value])
            .env("GITHUB_SHA", sha)
            .env("RUNNER_OS", runner)
            .output()
            .expect("run pure nightly refusal");
        assert!(!refusal.status.success());
        assert!(refusal.stdout.is_empty(), "refusal must not emit stdout");
        assert!(
            one_json_line(&refusal.stderr, "pure-refusal stderr").contains(expected_field),
            "wrong diagnostic precedence for {tag}"
        );
        assert_no_db_artifacts(&pure_db);
    }

    // Missing provenance is explicit typed data, never a fabricated commit,
    // machine, or ambiguous `local` sentinel.
    let unavailable_db = temp_db("nightly-unavailable-provenance");
    cleanup_db(&unavailable_db);
    let unavailable = std::process::Command::new(exe)
        .args([unavailable_db.as_str(), "error", "unavailable", "0"])
        .env_remove("GITHUB_SHA")
        .env_remove("RUNNER_OS")
        .output()
        .expect("run nightly_ledger without optional provenance variables");
    assert!(
        unavailable.status.success(),
        "typed-unavailable run failed: {}",
        String::from_utf8_lossy(&unavailable.stderr)
    );
    assert!(
        unavailable.stderr.is_empty(),
        "success must not emit stderr"
    );
    one_json_line(&unavailable.stdout, "typed-unavailable stdout");
    let unavailable_ledger = Ledger::open(&unavailable_db).expect("open unavailable ledger");
    let unavailable_op = unavailable_ledger.op(1).unwrap().expect("unavailable op");
    assert_eq!(
        unavailable_op.versions,
        "{\"frankensim\":{\"availability\":\"unavailable\",\"source\":\"process_environment\",\"variable\":\"GITHUB_SHA\"}}"
    );
    assert_eq!(
        unavailable_op.capability,
        "{\"ops\":[\"ci.nightly\"],\"runner\":{\"availability\":\"unavailable\",\"source\":\"process_environment\",\"variable\":\"RUNNER_OS\"}}"
    );
    assert!(!unavailable_op.versions.contains("local"));
    assert!(!unavailable_op.capability.contains("local"));
    assert!(unavailable_ledger.lint().unwrap().is_clean());
    drop(unavailable_ledger);
    cleanup_db(&unavailable_db);

    #[cfg(unix)]
    {
        // These bytes are valid POSIX filename data but reserved path syntax
        // on Windows. The resulting stdout must still be exactly one JSON line.
        let hostile_db = std::env::temp_dir().join(format!(
            "fs-ledger-nightly-path-{}-quote\"-slash\\-line\n-control\u{0001}.db",
            std::process::id()
        ));
        let hostile_db = hostile_db.to_str().expect("hostile path remains UTF-8");
        cleanup_db(hostile_db);
        let hostile_path = std::process::Command::new(exe)
            .args([hostile_db, "ok", "hostile-path", "1"])
            .env_remove("GITHUB_SHA")
            .env_remove("RUNNER_OS")
            .output()
            .expect("run nightly_ledger with hostile valid path");
        assert!(
            hostile_path.status.success(),
            "hostile path run failed: {}",
            String::from_utf8_lossy(&hostile_path.stderr)
        );
        assert!(
            hostile_path.stderr.is_empty(),
            "success must not emit stderr"
        );
        let hostile_stdout = one_json_line(&hostile_path.stdout, "hostile-path stdout");
        assert!(hostile_stdout.contains("quote\\\"-slash\\\\-line\\n-control\\u0001.db"));
        let hostile_ledger = Ledger::open(hostile_db).expect("reopen hostile path ledger");
        assert_eq!(hostile_ledger.table_count("ops").unwrap(), 1);
        assert!(hostile_ledger.lint().unwrap().is_clean());
        drop(hostile_ledger);
        cleanup_db(hostile_db);
    }

    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStringExt as _;

        let non_utf8_db = temp_db("nightly-non-utf8-suite");
        cleanup_db(&non_utf8_db);
        let non_utf8_suite = std::ffi::OsString::from_vec(vec![0xff]);
        let non_utf8 = std::process::Command::new(exe)
            .arg(non_utf8_db.as_str())
            .arg("ok")
            .arg(non_utf8_suite)
            .arg("1")
            .env_remove("GITHUB_SHA")
            .env_remove("RUNNER_OS")
            .output()
            .expect("run nightly_ledger with non-UTF-8 suite");
        assert!(!non_utf8.status.success());
        assert!(non_utf8.stdout.is_empty(), "refusal must not emit stdout");
        let non_utf8_stderr = one_json_line(&non_utf8.stderr, "non-UTF-8 stderr");
        assert!(non_utf8_stderr.contains("field 'suite' rejected: must be valid UTF-8"));
        assert_no_db_artifacts(&non_utf8_db);
        l.append_event(&EventRow {
            session: None,
            t: 5,
            kind: "nightly_writer_non_utf8_validation",
            payload: Some(non_utf8_stderr),
        })
        .expect("non-UTF-8 refusal is one valid JSON value");

        for (tag, variable) in [
            ("nightly-non-utf8-sha", "GITHUB_SHA"),
            ("nightly-non-utf8-runner", "RUNNER_OS"),
        ] {
            let env_db = temp_db(tag);
            cleanup_db(&env_db);
            let mut command = std::process::Command::new(exe);
            command.args([env_db.as_str(), "ok", "non-utf8-env", "1"]);
            if variable == "GITHUB_SHA" {
                command
                    .env("GITHUB_SHA", std::ffi::OsString::from_vec(vec![0xff]))
                    .env_remove("RUNNER_OS");
            } else {
                command
                    .env("GITHUB_SHA", "valid-sha")
                    .env("RUNNER_OS", std::ffi::OsString::from_vec(vec![0xfe]));
            }
            let invalid_env = command.output().expect("run with non-UTF-8 provenance");
            assert!(!invalid_env.status.success());
            assert!(
                invalid_env.stdout.is_empty(),
                "refusal must not emit stdout"
            );
            let stderr = one_json_line(&invalid_env.stderr, "non-UTF-8 provenance stderr");
            assert!(stderr.contains(&format!("field '{variable}' rejected: must be valid UTF-8")));
            assert_no_db_artifacts(&env_db);
        }
    }

    assert_eq!(l.table_count("ops").unwrap(), 1);
    assert_eq!(l.table_count("metrics").unwrap(), 1);
    #[cfg(unix)]
    assert_eq!(l.table_count("events").unwrap(), 6);
    #[cfg(not(unix))]
    assert_eq!(l.table_count("events").unwrap(), 5);
    assert!(l.lint().unwrap().is_clean());

    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-010",
        "nightly_ledger preserves hostile JSON, records typed provenance, and refuses every pure invalid input before database creation",
    );
}

#[test]
fn ledger_039_controlled_operation_discovery_is_bounded_and_replayable() {
    let db = temp_db("controlled-operation-discovery");
    let ledger = Ledger::open(&db).expect("open controlled discovery ledger");
    let first = ledger
        .begin_op_on(
            MAIN_BRANCH,
            ExecMode::Deterministic,
            None,
            r#"{"op":"prefix"}"#,
            &FX,
            1,
        )
        .expect("record shared prefix");
    let child = ledger
        .fork("controlled-conformance-child", MAIN_BRANCH)
        .expect("fork controlled child");
    let ir = format!(
        "{{\"payload\":\"{}\"}}",
        "x".repeat(CONTROLLED_OP_FIELD_TILE_LEN)
    );
    let selected = ledger
        .begin_op_on(
            child,
            ExecMode::Fast,
            Some(b"controlled-conformance-session"),
            &ir,
            &FX,
            2,
        )
        .expect("record selected child operation");
    ledger
        .finish_op(selected, OpOutcome::Ok, Some(r#"{"status":"complete"}"#), 3)
        .expect("finish selected operation");

    let first_page = ledger
        .visible_op_ids_page_controlled(child, None, 1, |_| {
            std::ops::ControlFlow::<()>::Continue(())
        })
        .expect("first bounded visible page");
    let ControlledVisibleOpPage::Complete(first_page) = first_page else {
        panic!("continue-only page controller cannot stop")
    };
    assert_eq!(first_page.op_ids, vec![selected]);
    let cursor = first_page.next_cursor.expect("shared prefix remains");

    let late = ledger
        .begin_op_on(
            child,
            ExecMode::Deterministic,
            None,
            r#"{"op":"late"}"#,
            &FX,
            4,
        )
        .expect("append after captured high-water");
    assert!(late > selected);
    let continuation = ledger
        .visible_op_ids_page_controlled(child, Some(cursor), 1, |_| {
            std::ops::ControlFlow::<()>::Continue(())
        })
        .expect("resume captured visible cohort");
    let ControlledVisibleOpPage::Complete(continuation) = continuation else {
        panic!("continue-only continuation cannot stop")
    };
    assert_eq!(continuation.op_ids, vec![first]);
    assert!(!continuation.truncated);

    let mut last_field = None;
    let mut next_offset = 0_u64;
    let mut maximum_tile = 0_usize;
    let read = ledger
        .read_op_fields_controlled(selected, |field, offset, tile| {
            if last_field == Some(field) {
                assert_eq!(offset, next_offset);
            } else {
                last_field = Some(field);
                next_offset = 0;
                assert_eq!(offset, 0);
            }
            maximum_tile = maximum_tile.max(tile.len());
            next_offset += u64::try_from(tile.len()).unwrap();
            std::ops::ControlFlow::<()>::Continue(())
        })
        .expect("controlled field read");
    let ControlledOpRead::Complete(complete) = read else {
        panic!("continue-only field controller cannot stop")
    };
    assert_eq!(complete.metadata.id, selected);
    assert_eq!(complete.metadata.branch, child);
    assert_eq!(complete.metadata.exec_mode, ExecMode::Fast);
    assert_eq!(complete.metadata.outcome, Some(OpOutcome::Ok));
    assert_eq!(last_field, Some(OpVariableField::Diagnostic));
    assert!(maximum_tile <= CONTROLLED_OP_FIELD_TILE_LEN);
    ledger
        .verify_op_content_identity_prehashed(&complete.prehashed_content)
        .expect("selected candidate sidecar agrees with delivered bytes");

    let generation_probe = ledger
        .begin_op_on(
            child,
            ExecMode::Deterministic,
            None,
            r#"{"op":"generation-probe"}"#,
            &FX,
            5,
        )
        .expect("record same-length source-change probe");
    let raw = SyncConnection::open(&db).expect("open independent fixture writer");
    let mut changed_source = false;
    let generation_error = ledger
        .read_op_fields_controlled::<()>(generation_probe, |field, _, _| {
            if !changed_source && field == OpVariableField::Ir {
                raw.execute_with_params(
                    "UPDATE ops SET seed = ?1 WHERE id = ?2",
                    &[
                        fsqlite::SqliteValue::Blob(vec![0xA5; FX.seed.len()].into()),
                        fsqlite::SqliteValue::Integer(generation_probe),
                    ],
                )
                .expect("replace seed through independent fixture writer");
                changed_source = true;
            }
            std::ops::ControlFlow::Continue(())
        })
        .expect_err("source generation change must refuse a mixed completion");
    assert!(changed_source);
    assert!(matches!(
        generation_error,
        LedgerError::OpCorrupt {
            op: rejected,
            detail,
        } if rejected == generation_probe
            && (detail.contains("source generation")
                || detail.contains("changed after bounded metadata preflight")
                || detail.contains("disappeared or changed"))
    ));

    let optional_probe = ledger
        .begin_op_on(
            child,
            ExecMode::Deterministic,
            None,
            r#"{"op":"optional-probe"}"#,
            &FX,
            6,
        )
        .expect("record optional-envelope probe");
    let mut finished_elsewhere = false;
    let optional_error = ledger
        .read_op_fields_controlled::<()>(optional_probe, |field, _, _| {
            if !finished_elsewhere && field == OpVariableField::Capability {
                raw.execute_with_params(
                    "UPDATE ops
                     SET t_end = 7, outcome = 'ok', diag = '{\"done\":true}'
                     WHERE id = ?1",
                    &[fsqlite::SqliteValue::Integer(optional_probe)],
                )
                .expect("finish operation through independent fixture writer");
                finished_elsewhere = true;
            }
            std::ops::ControlFlow::Continue(())
        })
        .expect_err("optional envelope change must refuse a stale completion");
    assert!(finished_elsewhere);
    assert!(matches!(
        optional_error,
        LedgerError::OpCorrupt {
            op: rejected,
            detail,
        } if rejected == optional_probe
            && detail.contains("envelope or source generation changed")
    ));

    drop(raw);
    drop(ledger);
    cleanup_db(&db);
    verdict(
        "ledger-039",
        "controlled visible-op pagination freezes a high-water cohort and controlled field reads verify only after complete bounded delivery",
    );
}
