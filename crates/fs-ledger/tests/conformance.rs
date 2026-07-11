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

use std::sync::atomic::{AtomicU32, Ordering};

use fs_ledger::{
    Blake3, EdgeRole, EventRow, FiveExplicits, Ledger, LedgerError, OpOutcome, SCHEMA_VERSION,
    STORAGE_CHUNK_LEN, hash_bytes,
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
    for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
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
        let raw = fsqlite::Connection::open(&db).expect("raw open");
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
fn ledger_004_dedupe_and_chunked_round_trip() {
    let db = temp_db("chunk");
    let l = Ledger::open(&db).expect("open");
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
    let report = l.verify_artifact_integrity().unwrap();
    assert_eq!(report.checked, 2);
    assert!(report.corrupted.contains(&small.hash.to_hex()));
    assert!(report.corrupted.contains(&big.hash.to_hex()));
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-005",
        "inline and chunked corruption both detected by re-hash",
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
        let raw = fsqlite::Connection::open(&db).expect("raw reopen");
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
        let raw = fsqlite::Connection::open(&db).expect("raw");
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
    let l = Ledger::open(&db).expect("open");
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
fn ledger_009_version_is_stamped() {
    assert!(!fs_ledger::VERSION.is_empty());
    verdict("ledger-009", "crate version stamped");
}

#[test]
fn ledger_010_nightly_writer_records_run() {
    let db = temp_db("nightly");
    let exe = env!("CARGO_BIN_EXE_nightly_ledger");
    let out = std::process::Command::new(exe)
        .args([db.as_str(), "ok", "nightly_gauntlet_pass", "1"])
        .output()
        .expect("run nightly_ledger");
    assert!(
        out.status.success(),
        "nightly_ledger failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let l = Ledger::open(&db).expect("open nightly db");
    assert_eq!(l.table_count("ops").unwrap(), 1);
    assert_eq!(l.table_count("metrics").unwrap(), 1);
    assert_eq!(l.table_count("events").unwrap(), 1);
    let op = l.op(1).unwrap().expect("op row");
    assert_eq!(op.outcome.as_deref(), Some("ok"));
    assert!(l.lint().unwrap().is_clean());
    // Bad arguments are structured refusals on stderr, never panics.
    let bad = std::process::Command::new(exe)
        .args([db.as_str(), "bogus", "x", "y"])
        .output()
        .expect("run nightly_ledger with bad args");
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("NightlyLedger"));
    drop(l);
    cleanup_db(&db);
    verdict(
        "ledger-010",
        "nightly_ledger records op+metric+event into a lint-clean ledger; bad args refuse structurally",
    );
}
