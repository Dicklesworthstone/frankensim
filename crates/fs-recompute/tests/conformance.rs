//! fs-recompute conformance suite (CONTRACT.md: any reimplementation
//! must pass). Node hashing stability and field sensitivity, the
//! determinism trip-wire, skip soundness with slack certificates, the
//! ACROSS-WORKER-COUNTS determinism certification (real threads,
//! adversarial completion orders — the G5-at-scale primitive), pinning
//! against eviction, and snapshot/fork stability. Completed cases emit
//! canonical aggregate fs-obs verdicts; assertions and expectations that
//! abort before that point remain ordinary Rust test diagnostics. The one
//! randomized case carries its actual input-generation seed, while
//! deterministic cases use zero. `NodeRecord::rng_seed` values describe
//! manifest records and are not treated as hidden test randomness.

use fs_exec::Reduce;
use fs_exec::reduce::{det_sum, pairwise_fold};
use fs_ledger::hash_bytes;
use fs_recompute::{
    ARTIFACT_CONTENT_IDENTITY_DOMAIN, ARTIFACT_CONTENT_IDENTITY_VERSION,
    NODE_RECORD_IDENTITY_DOMAIN, NODE_RECORD_IDENTITY_VERSION, NodeRecord, ParamValue, PinReason,
    PutOutcome, SNAPSHOT_FORMAT_DOMAIN, SNAPSHOT_FORMAT_VERSION, SkipDecision,
    SnapshotAdmissionError, SnapshotIdentity, Store, StoreError,
};
use std::sync::{Arc, Mutex};

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-recompute/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-recompute/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("recompute conformance verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line)
        .expect("recompute conformance verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }

    fn below(&mut self, n: u64) -> u64 {
        (self.next() >> 32) % n
    }
}

fn record(op: &str, seed: u64, achieved: f64, required: f64) -> NodeRecord {
    NodeRecord {
        op_id: op.to_string(),
        input_hashes: vec![hash_bytes(b"input-a"), hash_bytes(b"input-b")],
        params: vec![
            ("h".to_string(), ParamValue::f(0.05)),
            ("order".to_string(), ParamValue::Int(2)),
            ("scheme".to_string(), ParamValue::Str("weno".to_string())),
        ],
        code_version_hash: hash_bytes(b"code-v1"),
        rng_seed: seed,
        achieved_error: achieved,
        required_tolerance: required,
    }
}

const SNAPSHOT_V3_HEADER: &str = "fsrecompute v3\n\
node_identity_version=2\n\
node_identity_domain=fs-recompute-node-v2\n\
artifact_identity_version=1\n\
artifact_identity_domain=org.frankensim.fs-recompute.artifact-content.v1\n\
--\n";

#[test]
fn node_record_op_id_moves_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.op_id = "assemble-mass".to_string();
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn node_record_input_hashes_move_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.input_hashes[0] = hash_bytes(b"input-c");
    assert_ne!(base.content_hash(), changed.content_hash());
    let mut reordered = base.clone();
    reordered.input_hashes.reverse();
    assert_ne!(
        base.content_hash(),
        reordered.content_hash(),
        "input hashes are an ordered semantic sequence"
    );
}

#[test]
fn node_record_parameters_move_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.params[0].1 = ParamValue::f(0.051);
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn node_record_parameter_input_order_does_not_move_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut reordered = base.clone();
    reordered.params.reverse();
    assert_eq!(base.content_hash(), reordered.content_hash());
}

#[test]
fn node_record_code_version_moves_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.code_version_hash = hash_bytes(b"code-v2");
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn node_record_rng_seed_moves_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.rng_seed = 43;
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn node_record_achieved_error_moves_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.achieved_error = 2e-6;
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn node_record_required_tolerance_moves_identity() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let mut changed = base.clone();
    changed.required_tolerance = 2e-4;
    assert_ne!(base.content_hash(), changed.content_hash());
}

#[test]
fn artifact_content_bytes_move_identity() {
    let record = record("artifact-content", 42, 1e-6, 1e-4);
    let mut first = Store::new();
    let mut second = Store::new();
    first
        .put(record.clone(), b"artifact-a")
        .expect("first artifact must be admitted");
    second
        .put(record.clone(), b"artifact-b")
        .expect("second artifact must be admitted in an independent store");

    let first_hash = first
        .lookup(&record)
        .expect("first artifact stored")
        .artifact_hash;
    let second_hash = second
        .lookup(&record)
        .expect("second artifact stored")
        .artifact_hash;
    assert_ne!(first_hash, second_hash);
    assert_ne!(first_hash, hash_bytes(b"artifact-a"));
    assert_ne!(second_hash, hash_bytes(b"artifact-b"));
}

#[test]
fn snapshot_v3_admission_validates_identity_metadata_before_exposing_rows() {
    assert_eq!(SNAPSHOT_FORMAT_DOMAIN, "fsrecompute");
    assert_eq!(SNAPSHOT_FORMAT_VERSION, 3);
    assert_eq!(NODE_RECORD_IDENTITY_VERSION, 2);
    assert_eq!(NODE_RECORD_IDENTITY_DOMAIN, "fs-recompute-node-v2");
    assert_eq!(ARTIFACT_CONTENT_IDENTITY_VERSION, 1);
    assert_eq!(
        ARTIFACT_CONTENT_IDENTITY_DOMAIN,
        "org.frankensim.fs-recompute.artifact-content.v1"
    );

    let mut store = Store::new();
    store
        .put(record("snapshot-admission", 73, 1e-8, 1e-6), b"artifact")
        .expect("fixture row");
    let snapshot = store.snapshot();
    assert_eq!(
        snapshot,
        format!("{SNAPSHOT_V3_HEADER}{}\n", store.rows()[0])
    );
    let admitted = Store::admit_snapshot(&snapshot).expect("canonical v3 snapshot");
    assert_eq!(admitted.rows(), &[store.rows()[0].as_str()]);

    let opaque = format!("{SNAPSHOT_V3_HEADER}not-json\nstill-opaque");
    let admitted_opaque = Store::admit_snapshot(&opaque).expect("rows stay opaque");
    assert_eq!(admitted_opaque.rows(), &["not-json", "still-opaque"]);

    assert!(matches!(
        Store::admit_snapshot("fsrecompute v2\n"),
        Err(SnapshotAdmissionError::LegacyV2 { declared: 2 })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v1\n"),
        Err(SnapshotAdmissionError::StaleVersion {
            identity: SnapshotIdentity::SnapshotFormat,
            declared: 1,
            supported: 3,
        })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v4\n"),
        Err(SnapshotAdmissionError::FutureVersion {
            identity: SnapshotIdentity::SnapshotFormat,
            declared: 4,
            supported: 3,
        })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v3\nnode_identity_version=1\n"),
        Err(SnapshotAdmissionError::StaleVersion {
            identity: SnapshotIdentity::NodeRecord,
            declared: 1,
            supported: 2,
        })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v3\nnode_identity_version=3\n"),
        Err(SnapshotAdmissionError::FutureVersion {
            identity: SnapshotIdentity::NodeRecord,
            declared: 3,
            supported: 2,
        })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=fs-recompute-node-v2\n\
             artifact_identity_version=0\n"
        ),
        Err(SnapshotAdmissionError::StaleVersion {
            identity: SnapshotIdentity::ArtifactContent,
            declared: 0,
            supported: 1,
        })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=fs-recompute-node-v2\n\
             artifact_identity_version=2\n"
        ),
        Err(SnapshotAdmissionError::FutureVersion {
            identity: SnapshotIdentity::ArtifactContent,
            declared: 2,
            supported: 1,
        })
    ));

    assert!(matches!(
        Store::admit_snapshot("foreign-recompute v3\n"),
        Err(SnapshotAdmissionError::DomainMismatch {
            identity: SnapshotIdentity::SnapshotFormat,
            ..
        })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=shadow-node-domain\n"
        ),
        Err(SnapshotAdmissionError::DomainMismatch {
            identity: SnapshotIdentity::NodeRecord,
            ..
        })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=fs-recompute-node-v2\n\
             artifact_identity_version=1\n\
             artifact_identity_domain=shadow-artifact-domain\n"
        ),
        Err(SnapshotAdmissionError::DomainMismatch {
            identity: SnapshotIdentity::ArtifactContent,
            ..
        })
    ));

    assert!(matches!(
        Store::admit_snapshot("fsrecompute vthree\n"),
        Err(SnapshotAdmissionError::MalformedHeader { line: 1, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v3"),
        Err(SnapshotAdmissionError::MalformedHeader { line: 2, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_domain=fs-recompute-node-v2\n\
             node_identity_version=2\n"
        ),
        Err(SnapshotAdmissionError::NonCanonicalHeader { line: 2, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v03\n"),
        Err(SnapshotAdmissionError::NonCanonicalHeader { line: 1, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot("fsrecompute v3\nnode_identity_version=02\n"),
        Err(SnapshotAdmissionError::NonCanonicalHeader { line: 2, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=fs-recompute-node-v2\n\
             artifact_identity_version=01\n"
        ),
        Err(SnapshotAdmissionError::NonCanonicalHeader { line: 4, .. })
    ));
    assert!(matches!(
        Store::admit_snapshot(
            "fsrecompute v3\n\
             node_identity_version=2\n\
             node_identity_domain=fs-recompute-node-v2\n\
             artifact_identity_version=1\n\
             artifact_identity_domain=org.frankensim.fs-recompute.artifact-content.v1\n\
             extra-header-line\n\
             --\nrow\n"
        ),
        Err(SnapshotAdmissionError::NonCanonicalHeader { line: 6, .. })
    ));
}

/// rcs-001 — hashing: stable across repeats, sensitive to EVERY one of
/// the seven fields, param-order canonical, floats by bits; boundary
/// cases: negative slack representable, empty DAG, single node, and a
/// 1000-deep chain hash-stable.
#[test]
fn rcs_001_hashing_stability() {
    let base = record("assemble-stiffness", 42, 1e-6, 1e-4);
    let stable = base.content_hash() == base.content_hash();
    // Param order does not matter (canonicalized).
    let mut reordered = base.clone();
    reordered.params.reverse();
    let canonical = reordered.content_hash() == base.content_hash();
    // Every field matters.
    let mut probes = Vec::new();
    let mut m = base.clone();
    m.op_id = "assemble-mass".to_string();
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.input_hashes[0] = hash_bytes(b"input-c");
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.params[0].1 = ParamValue::f(0.1 / 2.0); // division by 2 is exact: same bits
    let same_bits = m.content_hash() == base.content_hash();
    let mut m = base.clone();
    m.params[0].1 = ParamValue::f(0.051);
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.code_version_hash = hash_bytes(b"code-v2");
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.rng_seed = 43;
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.achieved_error = 2e-6;
    probes.push(m.content_hash());
    let mut m = base.clone();
    m.required_tolerance = 2e-4;
    probes.push(m.content_hash());
    let all_differ = probes.iter().all(|p| *p != base.content_hash());
    // Delimiter injection must not alias a real input edge. The v1
    // line-oriented encoding serialized these two distinct records to
    // identical bytes.
    let injected_input = hash_bytes(b"injected-input");
    let mut injected = base.clone();
    injected.op_id = format!("assemble\nin:{}", injected_input.to_hex());
    injected.input_hashes.clear();
    let mut structured = base.clone();
    structured.op_id = "assemble".to_string();
    structured.input_hashes = vec![injected_input];
    let injection_safe = injected.content_hash() != structured.content_hash();
    // Negative slack is representable and first-class.
    let over_budget = record("expensive-op", 7, 1e-3, 1e-4);
    let over_budget_row = over_budget.to_row(b"x");
    let negative = over_budget.slack() < 0.0
        && over_budget_row.contains("\"slack\":-")
        && over_budget_row.contains("\"slack_bits\":");
    // Deep DAG: a 1000-node chain where each node's input is the
    // previous node's hash — stable across two builds.
    let deep = |n: u32| -> fs_ledger::ContentHash {
        let mut prev = hash_bytes(b"root");
        for k in 0..n {
            let r = NodeRecord {
                op_id: format!("step-{k}"),
                input_hashes: vec![prev],
                params: vec![],
                code_version_hash: hash_bytes(b"code-v1"),
                rng_seed: 0,
                achieved_error: 1e-9,
                required_tolerance: 1e-6,
            };
            prev = r.content_hash();
        }
        prev
    };
    let deep_stable = deep(1000) == deep(1000);
    // Empty and single-node stores behave.
    let empty = Store::new();
    let empty_ok = empty.is_empty() && empty.rows().is_empty();
    let mut single = Store::new();
    let _ = single.put(base.clone(), b"artifact").expect("put");
    let single_ok = single.len() == 1;
    verdict(
        "rcs-001",
        stable
            && canonical
            && same_bits
            && all_differ
            && injection_safe
            && negative
            && deep_stable
            && empty_ok
            && single_ok,
        &format!(
            "node hashes are repeat-stable, param-order canonical, and injective \
             across delimiter-looking caller strings; every field differs and \
             boundary stores behave [stable={stable} canonical={canonical} \
             same_bits={same_bits} all_differ={all_differ} injection_safe={injection_safe} \
             negative={negative} deep_stable={deep_stable} empty={empty_ok} \
             single={single_ok}]"
        ),
        0,
    );
}

/// rcs-002 — the determinism trip-wire's exact boundary: identical
/// SEVEN-FIELD record + identical artifact dedupes; identical record +
/// DIFFERENT artifact is a stop-the-line error naming both hashes. A
/// remeasured output changes the record identity and is stored separately,
/// so the diagnostic must not overstate this as a five-input-field check.
#[test]
fn rcs_002_determinism_tripwire() {
    let mut store = Store::new();
    let r = record("tile-reduce", 7, 1e-8, 1e-6);
    let first = store.put(r.clone(), b"bits-v1").expect("insert");
    let inserted = matches!(first, PutOutcome::Inserted(_));
    let again = store.put(r.clone(), b"bits-v1").expect("dedup");
    let deduped = matches!(again, PutOutcome::Deduped(_));
    let trip = store.put(r.clone(), b"bits-v2");
    let trip_message = trip
        .expect_err("same record with different bytes must trip")
        .to_string();
    let tripped = trip_message.contains("stop-the-line")
        && trip_message.contains("same seven-field node record")
        && trip_message.contains("achieved error, required tolerance")
        && !trip_message.contains("same (op, inputs, params, code, seed)");

    let mut remeasured = r;
    remeasured.achieved_error = 2e-8;
    let output_changed_identity = matches!(
        store
            .put(remeasured, b"bits-v2")
            .expect("a different measured-output field is a different seven-field node"),
        PutOutcome::Inserted(_)
    );
    verdict(
        "rcs-002",
        inserted && deduped && tripped && output_changed_identity,
        "identical (seven-field record, artifact) dedupes as a write-time memo hit; \
         the same full record with different bytes trips the DETERMINISM CONTRACT, \
         while a changed measured output is honestly a distinct node identity",
        0,
    );
}

/// rcs-003 — skip soundness: hits carry their slack certificate,
/// tightened tolerances name their deficit, unknown identities miss,
/// and the boundary (tolerance == achieved) is a zero-slack hit.
#[test]
fn rcs_003_skip_soundness() {
    let mut store = Store::new();
    let r = record("adapt-solve", 11, 1e-6, 1e-4);
    store.put(r.clone(), b"solution").expect("put");
    let hit = store.can_skip(&r, 1e-4);
    let hit_ok = matches!(hit, SkipDecision::Hit { slack, .. } if (slack - 9.9e-5).abs() < 1e-9);
    let boundary = store.can_skip(&r, 1e-6);
    let boundary_ok = matches!(boundary, SkipDecision::Hit { slack, .. } if slack == 0.0);
    let tightened = store.can_skip(&r, 1e-7);
    let tightened_ok = matches!(tightened, SkipDecision::ToleranceTightened { deficit }
        if (deficit - 9e-7).abs() < 1e-12);
    let invalid_tolerances = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0]
        .into_iter()
        .all(|tolerance| {
            matches!(
                store.can_skip(&r, tolerance),
                SkipDecision::InvalidTolerance { tolerance_bits }
                    if tolerance_bits == tolerance.to_bits()
            )
        });
    let mut other = r.clone();
    other.rng_seed = 999;
    let miss = store.can_skip(&other, 1e-4);
    let miss_ok = miss == SkipDecision::Miss;
    // The skip identity ignores the RECORDED tolerances (a node cached
    // under a looser requirement still hits when it achieved enough).
    let mut retolerated = r.clone();
    retolerated.required_tolerance = 5e-3;
    let still_hits = matches!(store.can_skip(&retolerated, 1e-4), SkipDecision::Hit { .. });
    verdict(
        "rcs-003",
        hit_ok && boundary_ok && tightened_ok && invalid_tolerances && miss_ok && still_hits,
        "skips carry slack certificates (9.9e-5 on the fixture), the exact boundary \
         is a zero-slack hit, tightened tolerances name their deficit (9e-7), \
         malformed tolerances fail closed, unknown identities miss, and the skip \
         identity correctly ignores recorded tolerances",
        0,
    );
}

/// rcs-004 — THE DETERMINISM CERTIFICATION (G5-at-scale primitive):
/// a fixture study (10k-element deterministic reduction over 64 fixed
/// tiles) produces BIT-IDENTICAL artifacts across {1,2,4,8} real
/// worker threads AND across adversarial permuted completion orders —
/// certified by the store accepting every re-put as a dedup.
#[test]
fn rcs_004_worker_count_certification() {
    const INPUT_SEED: u64 = 0x1001_2026_0707_0054;
    const COMPLETION_ORDER_SEEDS: [u64; 2] = [0xA1, 0xB2];

    #[derive(Clone, Copy)]
    struct Sum(f64);
    impl Reduce for Sum {
        fn identity() -> Self {
            Sum(0.0)
        }

        fn merge(self, other: Self) -> Self {
            Sum(self.0 + other.0)
        }
    }
    let data: Vec<f64> = {
        let mut rng = Lcg(INPUT_SEED);
        (0..10_240).map(|_| rng.unit() * 2.0 - 1.0).collect()
    };
    let tiles = 64usize;
    let tile_len = data.len() / tiles;
    // The study: tile partials via det_sum, global pairwise_fold in
    // FIXED tile order — the reduction tree never depends on which
    // worker computed what, or when it finished.
    let run = |workers: usize, permute_seed: Option<u64>| -> Vec<u8> {
        let partials: Arc<Mutex<Vec<Option<f64>>>> = Arc::new(Mutex::new(vec![None; tiles]));
        // Adversarial completion order: each worker processes its
        // tiles in a permuted order.
        let mut assignment: Vec<Vec<usize>> = vec![Vec::new(); workers];
        for t in 0..tiles {
            assignment[t % workers].push(t);
        }
        if let Some(seed) = permute_seed {
            let mut rng = Lcg(seed);
            for lane in &mut assignment {
                for i in (1..lane.len()).rev() {
                    let j = (rng.below(1 + i as u64)) as usize;
                    lane.swap(i, j);
                }
            }
        }
        std::thread::scope(|scope| {
            for lane in assignment {
                let partials = Arc::clone(&partials);
                let data = &data;
                scope.spawn(move || {
                    for t in lane {
                        let chunk = &data[t * tile_len..(t + 1) * tile_len];
                        let v = det_sum(chunk);
                        partials.lock().expect("lock")[t] = Some(v);
                    }
                });
            }
        });
        let finals: Vec<f64> = partials
            .lock()
            .expect("lock")
            .iter()
            .map(|v| v.expect("tile computed"))
            .collect();
        let total = pairwise_fold(finals.iter().map(|&v| Sum(v)).collect::<Vec<_>>()).0;
        total.to_le_bytes().to_vec()
    };
    let reference = run(1, None);
    let mut store = Store::new();
    let r = record("fixture-study-reduction", 0x54, 0.0, 1e-12);
    store.put(r.clone(), &reference).expect("reference put");
    let mut all_dedup = true;
    let mut runs = 0;
    for workers in [1usize, 2, 4, 8] {
        for permute in [
            None,
            Some(COMPLETION_ORDER_SEEDS[0]),
            Some(COMPLETION_ORDER_SEEDS[1]),
        ] {
            let artifact = run(workers, permute);
            runs += 1;
            match store.put(r.clone(), &artifact) {
                Ok(PutOutcome::Deduped(_)) => {}
                Ok(PutOutcome::Inserted(_)) => unreachable!("same record"),
                Err(_) => all_dedup = false,
            }
        }
    }
    verdict(
        "rcs-004",
        all_dedup && runs == 12,
        &format!(
            "the fixture study produced BIT-IDENTICAL artifacts across {{1,2,4,8}} \
             real worker threads x {{sequential, 2 adversarial permuted completion \
             orders}} = {runs} runs, every re-put accepted as a dedup by the \
             determinism contract (fixed tile partition + order-fixed pairwise \
             fold); fixture-input seed 0x1001_2026_0707_0054, completion-order \
             seeds [0xA1, 0xB2]; NodeRecord rng_seed 0x54 is manifest data"
        ),
        INPUT_SEED,
    );
}

/// rcs-005 — pinning: evidence-package/contract pins survive eviction;
/// eviction is deterministic (oldest unpinned first) and never touches
/// pinned nodes.
#[test]
fn rcs_005_pinning() {
    let mut store = Store::new();
    let mut hashes = Vec::new();
    for k in 0..10u64 {
        let r = record(&format!("op-{k}"), k, 1e-6, 1e-4);
        let PutOutcome::Inserted(h) = store.put(r, format!("art-{k}").as_bytes()).expect("put")
        else {
            unreachable!("fresh records");
        };
        hashes.push(h);
    }
    store
        .pin(&hashes[2], PinReason::EvidencePackage("EVP-7".to_string()))
        .expect("pin 2");
    store
        .pin(&hashes[5], PinReason::Contract("CTR-3".to_string()))
        .expect("pin 5");
    let evicted = store.evict_unpinned(3);
    let pinned_survive = store.get(&hashes[2]).is_some() && store.get(&hashes[5]).is_some();
    // Oldest unpinned first: 0, 1, 3, 4, 6 evicted (keep 3 unpinned).
    let expected_gone = [0usize, 1, 3, 4, 6]
        .iter()
        .all(|&i| store.get(&hashes[i]).is_none());
    let expected_kept = [7usize, 8, 9]
        .iter()
        .all(|&i| store.get(&hashes[i]).is_some());
    let unknown = store.pin(&hash_bytes(b"nope"), PinReason::Contract("x".to_string()));
    verdict(
        "rcs-005",
        evicted == 5 && pinned_survive && expected_gone && expected_kept && unknown.is_err(),
        "eviction removes exactly the 5 oldest UNPINNED nodes; evidence-package and \
         contract pins are untouchable; pinning unknown nodes teaches",
        0,
    );
}

/// rcs-006 — rows + fork stability: ledger rows carry all seven fields,
/// identity metadata, and slack; snapshots are bitwise-deterministic,
/// and the obs event ships the slack table.
#[test]
fn rcs_006_rows_and_fork() {
    let build = || -> (Store, Vec<String>) {
        let mut store = Store::new();
        for k in 0..5u64 {
            let r = record(&format!("op-{k}"), k, 1e-6 * (k + 1) as f64, 1e-4);
            store.put(r, format!("art-{k}").as_bytes()).expect("put");
        }
        let rows = store.rows().to_vec();
        (store, rows)
    };
    let (s1, r1) = build();
    let (s2, r2) = build();
    let rows_deterministic = r1 == r2;
    let snapshot1 = s1.snapshot();
    let snapshot2 = s2.snapshot();
    let fork_stable = snapshot1 == snapshot2
        && snapshot1.starts_with(SNAPSHOT_V3_HEADER)
        && Store::admit_snapshot(&snapshot1).is_ok();
    let has_fields = r1[0].contains("\"op\":")
        && r1[0].contains("\"node\":")
        && r1[0].contains("\"artifact\":")
        && r1[0].contains(
            "\"node_identity\":{\"version\":2,\"domain\":\"fs-recompute-node-v2\"}",
        )
        && r1[0].contains(
            "\"artifact_identity\":{\"version\":1,\"domain\":\"org.frankensim.fs-recompute.artifact-content.v1\"}",
        )
        && r1[0].contains("\"seed\":")
        && r1[0].contains("\"achieved\":")
        && r1[0].contains("\"required\":")
        && r1[0].contains("\"slack\":");
    let mut escaped_record = record("quoted \"op\"\nnext", 9, 1e-6, 1e-4);
    escaped_record.params.push((
        "quoted\"key".to_string(),
        ParamValue::Str("line1\nline2".to_string()),
    ));
    let escaped = escaped_record.to_row(b"escaped");
    let row_is_structured = escaped.contains("quoted \\\"op\\\"\\nnext")
        && escaped.contains("line1\\nline2")
        && !escaped.contains("line1\nline2")
        && escaped.contains("\"input_hashes\":[")
        && escaped.contains("\"params\":[")
        && escaped.contains("\"code_version\":");
    let mut em = fs_obs::Emitter::new("fs-recompute/conformance", "rcs-006/slack");
    let event = em.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: "recompute-slack-table".to_string(),
            json: format!("{{\"nodes\":{},\"rows\":{}}}", s1.len(), r1.len()),
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("slack table must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("slack table must use the fs-obs wire schema");
    println!("{line}");
    verdict(
        "rcs-006",
        rows_deterministic && fork_stable && has_fields && row_is_structured,
        "ledger rows carry all seven fields, exact node/artifact identity metadata, \
         and slack; identical builds give bitwise-identical rows, caller strings are \
         JSON-escaped, and canonical v3 snapshots are fork-stable and admissible",
        0,
    );
}

/// rcs-007 — malformed error magnitudes and slack burns fail closed.
#[test]
fn rcs_007_invalid_magnitudes_refuse() {
    let mut store = Store::new();
    let invalid_nan = store.put(record("nan", 1, f64::NAN, 1.0), b"nan");
    let invalid_negative = store.put(record("negative", 1, 0.0, -1.0), b"negative");
    let valid = record("valid", 1, 0.0, 1.0);
    let PutOutcome::Inserted(node) = store.put(valid, b"valid").expect("valid insert") else {
        unreachable!("fresh record");
    };
    let before = store.get(&node).expect("node").effective_slack();
    let burn = store.burn_slack(&node, -0.5);
    let after = store.get(&node).expect("node").effective_slack();
    let pass = matches!(invalid_nan, Err(StoreError::InvalidErrorBudget { .. }))
        && matches!(invalid_negative, Err(StoreError::InvalidErrorBudget { .. }))
        && matches!(burn, Err(StoreError::InvalidSlackBurn { .. }))
        && before.to_bits() == after.to_bits();
    verdict(
        "rcs-007",
        pass,
        "NaN/negative error budgets and negative slack burns refuse without mutating \
         the store; error magnitudes cannot manufacture reusable slack",
        0,
    );
}
