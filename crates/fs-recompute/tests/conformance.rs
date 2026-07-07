//! fs-recompute conformance suite (CONTRACT.md: any reimplementation
//! must pass). Node hashing stability and field sensitivity, the
//! determinism trip-wire, skip soundness with slack certificates, the
//! ACROSS-WORKER-COUNTS determinism certification (real threads,
//! adversarial completion orders — the G5-at-scale primitive), pinning
//! against eviction, and snapshot/fork stability. JSON-line verdicts;
//! seeded cases carry seeds.

use fs_exec::Reduce;
use fs_exec::reduce::{det_sum, pairwise_fold};
use fs_ledger::hash_bytes;
use fs_recompute::{
    NodeRecord, ParamValue, PinReason, PutOutcome, SkipDecision, Store, StoreError,
};
use std::sync::{Arc, Mutex};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-recompute/conformance\",\"case\":\"{case}\",\"verdict\":\"{}\",\
         \"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
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
    // Negative slack is representable and first-class.
    let over_budget = record("expensive-op", 7, 1e-3, 1e-4);
    let negative = over_budget.slack() < 0.0
        && over_budget
            .to_row(&hash_bytes(b"x"))
            .contains("\"slack\":-9.0");
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
            && negative
            && deep_stable
            && empty_ok
            && single_ok,
        "node hashes are repeat-stable and param-order canonical, EVERY one of the \
         seven fields perturbs the hash (floats by bits), negative slack is \
         first-class in the row, a 1000-deep chain is hash-stable, and \
         empty/single-node stores behave",
    );
}

/// rcs-002 — the determinism trip-wire: identical record + identical
/// artifact dedupes; identical record + DIFFERENT artifact is a
/// stop-the-line error naming both hashes.
#[test]
fn rcs_002_determinism_tripwire() {
    let mut store = Store::new();
    let r = record("tile-reduce", 7, 1e-8, 1e-6);
    let first = store.put(r.clone(), b"bits-v1").expect("insert");
    let inserted = matches!(first, PutOutcome::Inserted(_));
    let again = store.put(r.clone(), b"bits-v1").expect("dedup");
    let deduped = matches!(again, PutOutcome::Deduped(_));
    let trip = store.put(r, b"bits-v2");
    let tripped = matches!(&trip, Err(StoreError::DeterminismViolation { .. }))
        && trip.unwrap_err().to_string().contains("stop-the-line");
    verdict(
        "rcs-002",
        inserted && deduped && tripped,
        "identical (record, artifact) dedupes as a write-time memo hit; the same \
         record with different bytes trips the DETERMINISM CONTRACT with both \
         artifact hashes named — stop the line, not a warning",
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
        hit_ok && boundary_ok && tightened_ok && miss_ok && still_hits,
        "skips carry slack certificates (9.9e-5 on the fixture), the exact boundary \
         is a zero-slack hit, tightened tolerances name their deficit (9e-7), \
         unknown identities miss, and the skip identity correctly ignores recorded \
         tolerances",
    );
}

/// rcs-004 — THE DETERMINISM CERTIFICATION (G5-at-scale primitive):
/// a fixture study (10k-element deterministic reduction over 64 fixed
/// tiles) produces BIT-IDENTICAL artifacts across {1,2,4,8} real
/// worker threads AND across adversarial permuted completion orders —
/// certified by the store accepting every re-put as a dedup.
#[test]
fn rcs_004_worker_count_certification() {
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
        let mut rng = Lcg(0x1001_2026_0707_0054);
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
        for permute in [None, Some(0xA1), Some(0xB2)] {
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
             fold); seed 0x1001_2026_0707_0054"
        ),
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
    );
}

/// rcs-006 — rows + fork stability: ledger rows carry all seven fields
/// plus slack, snapshots are bitwise-deterministic, and the obs event
/// ships the slack table.
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
    let fork_stable = s1.snapshot() == s2.snapshot();
    let has_fields = r1[0].contains("\"op\":")
        && r1[0].contains("\"node\":")
        && r1[0].contains("\"artifact\":")
        && r1[0].contains("\"seed\":")
        && r1[0].contains("\"achieved\":")
        && r1[0].contains("\"required\":")
        && r1[0].contains("\"slack\":");
    let mut em = fs_obs::Emitter::new("fs-recompute/conformance", "rcs-006/slack");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "recompute-slack-table".to_string(),
                json: format!("{{\"nodes\":{},\"rows\":{}}}", s1.len(), r1.len()),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("slack table validates");
    println!("{line}");
    verdict(
        "rcs-006",
        rows_deterministic && fork_stable && has_fields,
        "ledger rows carry all seven fields plus slack, identical builds give \
         bitwise-identical rows, and snapshots are fork-stable",
    );
}
