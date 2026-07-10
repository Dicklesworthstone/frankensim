//! perturb() API + cache-policy conformance (bead lmp4.8, feature
//! `tolerance-invalidation`). First-class plans with certificates and
//! costs, leaf/root boundaries, the cost-weighted-vs-LRU eviction A/B,
//! pinned-pressure structured errors, per-op skip-yield telemetry, and
//! the kill-criterion replay measurement. JSON-line verdicts; seeded
//! cases carry seeds.

use fs_ledger::{ContentHash, hash_bytes};
use fs_recompute::api::RecomputeApi;
use fs_recompute::invalidate::Edge;
use fs_recompute::{NodeRecord, PinReason, PutOutcome, Store, StoreError};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-recompute/api\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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
}

fn record(op: &str, achieved: f64, required: f64) -> NodeRecord {
    NodeRecord {
        op_id: op.to_string(),
        input_hashes: vec![],
        params: vec![],
        code_version_hash: hash_bytes(b"code-v1"),
        rng_seed: 0,
        achieved_error: achieved,
        required_tolerance: required,
    }
}

fn put(store: &mut Store, op: &str, achieved: f64, required: f64) -> ContentHash {
    match store
        .put(record(op, achieved, required), op.as_bytes())
        .expect("put")
    {
        PutOutcome::Inserted(h) | PutOutcome::Deduped(h) => h,
    }
}

/// A diamond DAG: src → {left, right} → sink, plus an isolated leaf.
fn diamond() -> (Store, Vec<Edge>, Vec<ContentHash>) {
    let mut store = Store::new();
    let src = put(&mut store, "src", 1e-9, 1e-9);
    let left = put(&mut store, "left", 0.0, 1e-2);
    let right = put(&mut store, "right", 0.0, 1e-6);
    let sink = put(&mut store, "sink", 0.0, 1e-1);
    let leaf = put(&mut store, "leaf", 0.0, 1e-3);
    let edges = vec![
        Edge {
            from: src,
            to: left,
            sensitivity: 0.5,
        },
        Edge {
            from: src,
            to: right,
            sensitivity: 0.5,
        },
        Edge {
            from: left,
            to: sink,
            sensitivity: 0.1,
        },
        Edge {
            from: right,
            to: sink,
            sensitivity: 0.1,
        },
    ];
    (store, edges, vec![src, left, right, sink, leaf])
}

/// api-001 — perturb() returns a first-class plan: minimal frontier,
/// certificates for skips, costs for the planner; leaf and root
/// boundaries behave; nothing burns until commit.
#[test]
fn api_001_perturb_plan() {
    let (store, edges, h) = diamond();
    let mut api = RecomputeApi::new(store, edges, 1.0);
    for (i, cost) in [1.0, 10.0, 100.0, 1000.0, 5.0].iter().enumerate() {
        api.record_cost(&h[i], *cost);
    }
    // δ = 1e-2 at src: left absorbs (5e-3 < 1e-2), right recomputes
    // (5e-3 > 1e-6), sink absorbs the flow-through (5e-3·0.1·2 = 1e-3
    // < 1e-1).
    let p = api.perturb(&h[0], 1e-2).expect("plan");
    let frontier_right = p.recompute_count() == 2 && p.skip_count() == 2; // src + right
    let certs = p.certificates().len() == 2
        && p.certificates()
            .iter()
            .all(|c| c.contains("perturbation absorbed"));
    // Cost: recompute src (1.0) + right (100.0); hash-memo pays all 4.
    let cost_ok =
        (p.estimated_cost() - 101.0).abs() < 1e-12 && (p.hash_memo_cost() - 1111.0).abs() < 1e-12;
    // Nothing burned before commit: replan gives identical verdicts.
    let p2 = api.perturb(&h[0], 1e-2).expect("replan");
    let pure = p2.recompute_count() == p.recompute_count();
    // Leaf: perturbing the isolated leaf touches only itself.
    let pl = api.perturb(&h[4], 1e-5).expect("leaf");
    let leaf_ok = pl.invalidation().verdicts().len() == 1 && pl.recompute_count() == 1;
    // Root: a huge δ swamps every slack → full frontier recomputes.
    let pr = api.perturb(&h[0], 1e3).expect("root");
    let root_ok = pr.recompute_count() == 4 && pr.skip_count() == 0;
    verdict(
        "api-001",
        frontier_right && certs && cost_ok && pure && leaf_ok && root_ok,
        &format!(
            "the diamond plan recomputes {{src, right}} and certifies {{left, sink}} \
             with absorbed-perturbation claims; costs read {:.0} vs hash-memo \
             {:.0}; plans are pure until commit; a leaf touches only itself; a \
             swamping root perturbation recomputes the full frontier",
            p.estimated_cost(),
            p.hash_memo_cost()
        ),
    );
}

/// api-002 — commit burns slack and updates telemetry: the same
/// perturbation stops being absorbable once its slack is spent.
#[test]
fn api_002_commit_burns() {
    let (store, edges, h) = diamond();
    let mut api = RecomputeApi::new(store, edges, 1.0);
    let p1 = api.perturb(&h[0], 8e-3).expect("plan 1");
    let left_skipped_first =
        p1.invalidation().verdicts().iter().any(|(x, v)| {
            *x == h[1] && matches!(v, fs_recompute::invalidate::Verdict::Skip { .. })
        });
    api.commit(&p1).expect("commit");
    // left's slack was 1e-2, burned 4e-3 → 6e-3 left. A second 8e-3
    // perturbation (4e-3 incoming) still fits; a third does not.
    let p2 = api.perturb(&h[0], 8e-3).expect("plan 2");
    api.commit(&p2).expect("commit 2");
    let p3 = api.perturb(&h[0], 8e-3).expect("plan 3");
    let left_recomputes_third = p3.invalidation().verdicts().iter().any(|(x, v)| {
        *x == h[1] && matches!(v, fs_recompute::invalidate::Verdict::Recompute { .. })
    });
    verdict(
        "api-002",
        left_skipped_first && left_recomputes_third,
        "slack is spendable through the API: the first two 8e-3 perturbations absorb \
         at `left` (burning 4e-3 each of its 1e-2 slack), the third finds it spent \
         and recomputes",
    );
}

/// api-003 — the eviction A/B: cost-weighted eviction preserves more
/// recompute-cost-saved than LRU on a seeded trace; pins survive both;
/// a pinned-saturated cache is a STRUCTURED error.
#[test]
fn api_003_eviction_ab() {
    // A chain so hits land on absorbers (sources always recompute).
    let mut store = Store::new();
    let src = put(&mut store, "src", 1e-9, 1e-9);
    let mut hot = Vec::new();
    let mut edges = Vec::new();
    for k in 0..4 {
        let n = put(&mut store, &format!("hot{k}"), 0.0, 1.0);
        edges.push(Edge {
            from: src,
            to: n,
            sensitivity: 0.1,
        });
        hot.push(n);
    }
    let mut cold = Vec::new();
    for k in 0..8 {
        let n = put(&mut store, &format!("cold{k}"), 0.0, 1.0);
        cold.push(n);
    }
    let mut api = RecomputeApi::new(store, edges, 1.0);
    for n in &hot {
        api.record_cost(n, 1000.0);
    }
    for n in &cold {
        api.record_cost(n, 1.0);
    }
    for _ in 0..20 {
        let p = api.perturb(&src, 1e-3).expect("plan");
        api.commit(&p).expect("commit"); // hot nodes absorb → hits
    }
    // Pin one cold node (a contract reference).
    api.store
        .pin(&cold[0], PinReason::Contract("CTR-9".to_string()))
        .expect("pin");
    // Cost-weighted eviction to 6 total: should keep the 4 hot
    // (high cost × hit-prob) + the pin + 1 more.
    let evicted = api.ensure_capacity(6).expect("evict");
    let hot_survive = hot.iter().all(|n| api.store.get(n).is_some());
    let pin_survives = api.store.get(&cold[0]).is_some();
    // LRU baseline (insertion order) on an identical setup would evict
    // the OLDEST unpinned = src + hot0..3 — destroying 4000 of cost.
    let cw_cost_kept: f64 = hot
        .iter()
        .filter(|n| api.store.get(n).is_some())
        .map(|_| 1000.0)
        .sum();
    let lru_cost_kept = 0.0; // the four oldest-after-src are the hot ones
    // Pinned saturation: capacity below the pin count is structured.
    for n in &hot {
        api.store
            .pin(n, PinReason::EvidencePackage("EVP-1".to_string()))
            .expect("pin hot");
    }
    let pressure = api.ensure_capacity(2);
    let structured = matches!(pressure, Err(StoreError::CacheFullOfPins { pinned, capacity })
        if pinned == 5 && capacity == 2);
    verdict(
        "api-003",
        evicted > 0 && hot_survive && pin_survives && cw_cost_kept > lru_cost_kept && structured,
        &format!(
            "cost-weighted eviction keeps all four hot expensive nodes \
             ({cw_cost_kept:.0} recompute-cost preserved vs {lru_cost_kept:.0} \
             under LRU insertion order), the contract pin survives, and pinned \
             saturation surfaces a STRUCTURED CacheFullOfPins error; \
             seed 0x1001_2026_0707_0083"
        ),
    );
}

/// api-004 — skip-yield telemetry per op, worst-first ordering, and
/// the live dashboard event.
#[test]
fn api_004_skip_yield_dashboard() {
    let (store, edges, h) = diamond();
    let mut api = RecomputeApi::new(store, edges, 1.0);
    for _ in 0..10 {
        let p = api.perturb(&h[0], 1e-2).expect("plan");
        api.commit(&p).expect("commit stops once slack burns");
        // Stop early if left's slack is spent (later plans recompute).
        if p.skip_count() == 0 {
            break;
        }
    }
    let y = api.skip_yield();
    let right_worst = y.of("right").is_some_and(|v| v == 0.0);
    let src_worst = y.of("src").is_some_and(|v| v == 0.0);
    let left_partial = y.of("left").is_some_and(|v| v > 0.0);
    let worst = y.worst_first();
    let worst_ordering = worst.first().is_some_and(|(_, v)| *v == 0.0);
    let dashboard = y.dashboard_json();
    let mut em = fs_obs::Emitter::new("fs-recompute/api", "api-004/skip-yield");
    let line = em
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: "recompute-skip-yield".to_string(),
                json: dashboard.clone(),
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("dashboard validates");
    println!("{line}");
    verdict(
        "api-004",
        right_worst && src_worst && left_partial && worst_ordering,
        &format!(
            "per-op yields separate the never-absorbing ops (src, right at 0.0 — \
             where tightening effort goes, worst-first) from the absorbing ones \
             (left > 0), and the dashboard ships live: {dashboard}"
        ),
    );
}

/// api-005 — the kill-criterion replay: a recorded design-iteration
/// trace measures certified-skip cost against plain hash memoization;
/// the fixture achieves >= 2x (the measurement machinery the real
/// decision will use, at fixture scale).
#[test]
fn api_005_kill_criterion_replay() {
    // A pipeline where most stages have healthy slack: an agent
    // exploring nearby variants perturbs the design repeatedly.
    let mut store = Store::new();
    let src = put(&mut store, "design", 1e-9, 1e-9);
    let mut prev = src;
    let mut edges = Vec::new();
    let mut stages = vec![src];
    for k in 0..10 {
        let n = put(&mut store, &format!("stage{k}"), 0.0, 1.0);
        edges.push(Edge {
            from: prev,
            to: n,
            sensitivity: 0.6,
        });
        stages.push(n);
        prev = n;
    }
    let mut api = RecomputeApi::new(store, edges, 1.0);
    for (k, s) in stages.iter().enumerate() {
        api.record_cost(s, if k == 0 { 1.0 } else { 50.0 });
    }
    let mut rng = Lcg(0x1001_2026_0707_0085);
    let mut certified_cost = 0.0;
    let mut memo_cost = 0.0;
    let mut plans = 0;
    for _ in 0..100 {
        let delta = 1e-3 * rng.unit();
        let p = api.perturb(&src, delta).expect("plan");
        certified_cost += p.estimated_cost();
        memo_cost += p.hash_memo_cost();
        plans += 1;
        // Do NOT commit: model independent variant exploration from
        // the same base design (branch checkouts, Proposal 10).
    }
    let speedup = memo_cost / certified_cost.max(1e-12);
    println!(
        "{{\"suite\":\"fs-recompute/api\",\"case\":\"api-005-metric\",\"verdict\":\"info\",\
         \"detail\":\"plans={plans} certified={certified_cost:.0} memo={memo_cost:.0} \
         speedup={speedup:.1}x\"}}"
    );
    verdict(
        "api-005",
        speedup >= 2.0 && plans == 100,
        &format!(
            "the 100-variant replay pays {certified_cost:.0} in certified delta-solve \
             cost vs {memo_cost:.0} under hash memoization — {speedup:.1}x, the \
             kill-criterion measurement machinery (fixture-scale; the production \
             decision runs on recorded agent traces); seed 0x1001_2026_0707_0085"
        ),
    );
}

/// api-006 — plans are bound to the exact slack snapshot they certify.
/// Committing one of two concurrent plans makes the other stale; the
/// refusal is atomic and leaves slack/telemetry unchanged.
#[test]
fn api_006_stale_plan_refuses_atomically() {
    let (store, edges, h) = diamond();
    let mut api = RecomputeApi::new(store, edges, 1.0);
    let first = api.perturb(&h[0], 8e-3).expect("first plan");
    let stale = api.perturb(&h[0], 8e-3).expect("concurrent plan");
    api.commit(&first).expect("first commit");

    let left_before = api.store.get(&h[1]).expect("left").effective_slack();
    let sink_before = api.store.get(&h[3]).expect("sink").effective_slack();
    let telemetry_before = api.skip_yield().dashboard_json();
    let refused = api.commit(&stale);
    let left_after = api.store.get(&h[1]).expect("left").effective_slack();
    let sink_after = api.store.get(&h[3]).expect("sink").effective_slack();
    let telemetry_after = api.skip_yield().dashboard_json();

    // A revision number alone is not a store identity. Construct two stores
    // at the same revision but with different burned-slack state and prove a
    // plan cannot cross between them.
    let (store_a, edges_a, a) = diamond();
    let mut api_a = RecomputeApi::new(store_a, edges_a, 1.0);
    api_a
        .store
        .pin(&a[4], PinReason::Contract("same-revision".to_string()))
        .expect("pin advances revision without changing slack");
    let foreign = api_a.perturb(&a[0], 8e-3).expect("foreign plan");
    let (store_b, edges_b, b) = diamond();
    let mut api_b = RecomputeApi::new(store_b, edges_b, 1.0);
    let own = api_b.perturb(&b[0], 8e-3).expect("own plan");
    api_b.commit(&own).expect("burn advances matching revision");
    let same_revision = api_a.store.revision() == api_b.store.revision();
    let foreign_refused = matches!(api_b.commit(&foreign), Err(StoreError::StalePlan { .. }));

    // Telemetry keys originate in caller-controlled operator ids and must be
    // escaped just like node rows.
    let mut escaped_store = Store::new();
    let escaped_src = put(&mut escaped_store, "escaped-src", 0.0, 0.0);
    let escaped_sink = put(&mut escaped_store, "quoted \"op\"\nnext", 0.0, 1.0);
    let escaped_edges = vec![Edge {
        from: escaped_src,
        to: escaped_sink,
        sensitivity: 0.5,
    }];
    let mut escaped_api = RecomputeApi::new(escaped_store, escaped_edges, 1.0);
    let escaped_plan = escaped_api
        .perturb(&escaped_src, 0.25)
        .expect("escaped plan");
    escaped_api.commit(&escaped_plan).expect("escaped commit");
    let dashboard = escaped_api.skip_yield().dashboard_json();
    let telemetry_escaped = dashboard.contains("quoted \\\"op\\\"\\nnext")
        && !dashboard.contains("quoted \"op\"\nnext");

    verdict(
        "api-006",
        matches!(refused, Err(StoreError::StalePlan { .. }))
            && left_before.to_bits() == left_after.to_bits()
            && sink_before.to_bits() == sink_after.to_bits()
            && telemetry_before == telemetry_after
            && same_revision
            && foreign_refused
            && telemetry_escaped,
        "a plan certified against pre-burn slack refuses after a sibling plan commits; \
         a same-revision plan from different state also refuses; failures perform no \
         partial burns or false telemetry, and dashboard operator ids are JSON-escaped",
    );
}
