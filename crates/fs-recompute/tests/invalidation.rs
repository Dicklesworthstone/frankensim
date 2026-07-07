//! Tolerance-aware invalidation conformance (bead lmp4.7, feature
//! `tolerance-invalidation`). Propagation and absorption on chains,
//! the fail-closed hardening zoo, the G3 SOUNDNESS battery against a
//! REAL executable DAG (any unsound skip is Sev-0), the falsifier,
//! graceful degradation under loose bounds with skip-yield measured,
//! and slack burning. JSON-line verdicts; seeded cases carry seeds.

use fs_ledger::{ContentHash, hash_bytes};
use fs_recompute::invalidate::{Edge, RecomputeReason, Verdict, apply_plan, plan};
use fs_recompute::{NodeRecord, PutOutcome, Store};

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-recompute/invalidation\",\"case\":\"{case}\",\"verdict\":\"{}\",\
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

fn simple_record(op: &str, achieved: f64, required: f64) -> NodeRecord {
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

fn put(store: &mut Store, op: &str, achieved: f64, required: f64, artifact: &[u8]) -> ContentHash {
    match store.put(simple_record(op, achieved, required), artifact).expect("put") {
        PutOutcome::Inserted(h) | PutOutcome::Deduped(h) => h,
    }
}

/// inv-001 — propagation and absorption: on a decaying chain the
/// perturbation is absorbed exactly where the accumulated bound first
/// fits the slack, and propagation STOPS there; δ = 0 is an empty
/// frontier.
#[test]
fn inv_001_absorption() {
    let mut store = Store::new();
    // Chain a → b → c → d with L = 0.5 on each edge.
    // Slacks: a: 0 (source), b: 1e-3, c: 1e-1, d: 1e-6.
    let a = put(&mut store, "a", 1e-9, 1e-9, b"a");
    let b = put(&mut store, "b", 1e-6, 1.001e-3, b"b");
    let c = put(&mut store, "c", 1e-6, 1.001e-1, b"c");
    let d = put(&mut store, "d", 1e-9, 1.001e-6, b"d");
    let edges = vec![
        Edge { from: a, to: b, sensitivity: 0.5 },
        Edge { from: b, to: c, sensitivity: 0.5 },
        Edge { from: c, to: d, sensitivity: 0.5 },
    ];
    // δ = 0.02 at a: a recomputes (source); b sees 0.01 > slack 1e-3
    // → recompute, passes 0.01·0.5... b's bound = 0.01, b recomputes,
    // c sees 0.5·0.01 = 5e-3 < 1e-1 → SKIP, absorbed; d untouched.
    let p = plan(&store, &edges, &[(a, 0.02)]).expect("plan");
    let touched: Vec<_> = p.verdicts.iter().map(|(h, v)| (*h, v.clone())).collect();
    let a_recomputes = matches!(
        touched.iter().find(|(h, _)| *h == a),
        Some((_, Verdict::Recompute { reason: RecomputeReason::SourcePerturbed, .. }))
    );
    let b_recomputes = matches!(
        touched.iter().find(|(h, _)| *h == b),
        Some((_, Verdict::Recompute { reason: RecomputeReason::BoundExceedsSlack, .. }))
    );
    let c_skips = matches!(
        touched.iter().find(|(h, _)| *h == c),
        Some((_, Verdict::Skip { bound, .. })) if (bound - 5e-3).abs() < 1e-15
    );
    let d_untouched = !touched.iter().any(|(h, _)| *h == d);
    // δ = 0: empty frontier.
    let empty = plan(&store, &edges, &[(a, 0.0)]).expect("plan");
    let empty_ok = empty.verdicts.is_empty() && (empty.skip_yield - 1.0).abs() < 1e-15;
    verdict(
        "inv-001",
        a_recomputes && b_recomputes && c_skips && d_untouched && empty_ok,
        "the edit recomputes, the tight node recomputes, the slack-rich node ABSORBS \
         (bound 5e-3, propagation stops), the downstream node is untouched, and a \
         zero perturbation is an empty frontier",
    );
}

/// inv-002 — fail-closed hardening: exact ties recompute; negative
/// slack never skips; non-finite sensitivities force recompute;
/// non-topological edges refuse.
#[test]
fn inv_002_fail_closed() {
    let mut store = Store::new();
    let a = put(&mut store, "a2", 1e-9, 1e-9, b"a2");
    // Tie: slack exactly equals incoming bound 0.005.
    let b = put(&mut store, "b2", 0.0, 0.005, b"b2");
    // Negative slack.
    let c = put(&mut store, "c2", 0.01, 0.001, b"c2");
    // Downstream of a non-finite edge.
    let d = put(&mut store, "d2", 0.0, 1e9, b"d2");
    let edges = vec![
        Edge { from: a, to: b, sensitivity: 0.5 },
        Edge { from: a, to: c, sensitivity: 1e-9 },
        Edge { from: a, to: d, sensitivity: f64::INFINITY },
    ];
    let p = plan(&store, &edges, &[(a, 0.01)]).expect("plan");
    let find = |h: ContentHash| p.verdicts.iter().find(|(x, _)| *x == h).map(|(_, v)| v.clone());
    let tie = matches!(find(b),
        Some(Verdict::Recompute { reason: RecomputeReason::TieOnBoundary, .. }));
    let negative = matches!(find(c),
        Some(Verdict::Recompute { reason: RecomputeReason::NegativeSlack, .. }));
    let nonfinite = matches!(find(d),
        Some(Verdict::Recompute { reason: RecomputeReason::NonFiniteSensitivity, .. }));
    // Non-topological edge refuses.
    let bad = plan(
        &store,
        &[Edge { from: b, to: a, sensitivity: 1.0 }],
        &[(a, 0.01)],
    );
    let refuses = bad.is_err();
    verdict(
        "inv-002",
        tie && negative && nonfinite && refuses,
        "the exact tie recomputes (never skip on the boundary), negative slack never \
         skips, non-finite sensitivities FAIL CLOSED, and non-topological edges \
         refuse with teaching text",
    );
}

/// An executable fixture DAG: values + ops with TRUE Lipschitz bounds.
/// Layout: x0, x1 inputs; m = 0.5·x0 + 0.5·x1; t = tanh(m); s = 0.3·t;
/// q = t + s; r = 0.25·q.
struct Fixture {
    store: Store,
    hashes: Vec<ContentHash>,
    edges: Vec<Edge>,
}

fn eval_fixture(x0: f64, x1: f64) -> [f64; 7] {
    let m = 0.5 * x0 + 0.5 * x1;
    let t = m.tanh();
    let s = 0.3 * t;
    let q = t + s;
    let r = 0.25 * q;
    [x0, x1, m, t, s, q, r]
}

fn build_fixture(x0: f64, x1: f64, tolerances: &[f64; 7]) -> Fixture {
    let vals = eval_fixture(x0, x1);
    let mut store = Store::new();
    let names = ["x0", "x1", "m", "t", "s", "q", "r"];
    let mut hashes = Vec::new();
    for (i, name) in names.iter().enumerate() {
        let h = put(
            &mut store,
            name,
            0.0, // baseline achieves exactly
            tolerances[i],
            &vals[i].to_le_bytes(),
        );
        hashes.push(h);
    }
    // True Lipschitz bounds: linear coefficients exact; tanh' ≤ 1.
    let edges = vec![
        Edge { from: hashes[0], to: hashes[2], sensitivity: 0.5 },
        Edge { from: hashes[1], to: hashes[2], sensitivity: 0.5 },
        Edge { from: hashes[2], to: hashes[3], sensitivity: 1.0 },
        Edge { from: hashes[3], to: hashes[4], sensitivity: 0.3 },
        Edge { from: hashes[3], to: hashes[5], sensitivity: 1.0 },
        Edge { from: hashes[4], to: hashes[5], sensitivity: 1.0 },
        Edge { from: hashes[5], to: hashes[6], sensitivity: 0.25 },
    ];
    Fixture {
        store,
        hashes,
        edges,
    }
}

/// inv-003 — THE SOUNDNESS BATTERY (G3; any violation is Sev-0): over
/// seeded random perturbation traces on the executable DAG, every
/// SKIPPED node's cached value lies within its required tolerance of
/// the full-recompute ground truth; the falsifier force-recomputes
/// skipped nodes and asserts agreement within the certified bound.
#[test]
fn inv_003_soundness_battery() {
    let mut rng = Lcg(0x1001_2026_0707_0063);
    let mut skipped_total = 0u32;
    let mut violations = 0u32;
    let mut falsifier_checks = 0u32;
    let mut falsifier_failures = 0u32;
    for _ in 0..40 {
        let x0 = rng.unit() * 2.0 - 1.0;
        let x1 = rng.unit() * 2.0 - 1.0;
        // Tolerances: mixed tight/loose so skips and recomputes co-occur.
        let tolerances = [
            1e-12,
            1e-12,
            10f64.powf(-(2.0 + 4.0 * rng.unit())),
            10f64.powf(-(2.0 + 4.0 * rng.unit())),
            10f64.powf(-(2.0 + 4.0 * rng.unit())),
            10f64.powf(-(2.0 + 4.0 * rng.unit())),
            10f64.powf(-(2.0 + 4.0 * rng.unit())),
        ];
        let fx = build_fixture(x0, x1, &tolerances);
        let baseline = eval_fixture(x0, x1);
        // Perturb one input by a small δ.
        let which = usize::from(rng.unit() < 0.5);
        let delta = (rng.unit() * 2.0 - 1.0) * 1e-3;
        let (nx0, nx1) = if which == 0 {
            (x0 + delta, x1)
        } else {
            (x0, x1 + delta)
        };
        let truth = eval_fixture(nx0, nx1);
        let p = plan(&fx.store, &fx.edges, &[(fx.hashes[which], delta.abs())])
            .expect("plan");
        for (h, v) in &p.verdicts {
            let idx = fx.hashes.iter().position(|x| x == h).expect("fixture node");
            if let Verdict::Skip { bound, .. } = v {
                skipped_total += 1;
                // SOUNDNESS: cached (baseline) value within required
                // tolerance of the true new value.
                let err = (baseline[idx] - truth[idx]).abs();
                if err > fx.store.get(h).expect("node").record.required_tolerance {
                    violations += 1;
                }
                // Falsifier: force-recompute ~30% of skips; agreement
                // must be within the certified bound.
                if rng.unit() < 0.3 {
                    falsifier_checks += 1;
                    if err > *bound {
                        falsifier_failures += 1;
                    }
                }
            }
        }
    }
    verdict(
        "inv-003",
        violations == 0 && falsifier_failures == 0 && skipped_total > 30 && falsifier_checks > 5,
        &format!(
            "0 unsound skips over 40 seeded traces ({skipped_total} skips, every \
             cached value within tolerance of full-recompute truth — the Sev-0 \
             gate) and the falsifier's {falsifier_checks} forced recomputes all \
             agreed within their certified bounds; seed 0x1001_2026_0707_0063"
        ),
    );
}

/// inv-004 — graceful degradation (risk R4): a loose-bound op mid-DAG
/// balloons the frontier to hash-memoization behavior — still CORRECT
/// — and the skip-yield metric records the collapse.
#[test]
fn inv_004_graceful_degradation() {
    let tolerances = [1e-12, 1e-12, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2];
    let fx = build_fixture(0.3, -0.2, &tolerances);
    let healthy = plan(&fx.store, &fx.edges, &[(fx.hashes[0], 1e-4)]).expect("plan");
    // Loosen: replace m→t sensitivity with a pessimistic 1e6 (the
    // interval-derivative worst case for a nasty op).
    let mut loose_edges = fx.edges.clone();
    loose_edges[2].sensitivity = 1e6;
    let degraded = plan(&fx.store, &loose_edges, &[(fx.hashes[0], 1e-4)]).expect("plan");
    let healthy_skips = healthy.skip_yield;
    let degraded_skips = degraded.skip_yield;
    // Degraded: everything downstream of t recomputes (bounds blown).
    let all_downstream_recompute = degraded
        .verdicts
        .iter()
        .filter(|(h, _)| {
            let idx = fx.hashes.iter().position(|x| x == h).expect("node");
            idx >= 3
        })
        .all(|(_, v)| matches!(v, Verdict::Recompute { .. }));
    println!(
        "{{\"suite\":\"fs-recompute/invalidation\",\"case\":\"inv-004-yield\",\
         \"verdict\":\"info\",\"detail\":\"healthy={healthy_skips:.2} \
         degraded={degraded_skips:.2}\"}}"
    );
    verdict(
        "inv-004",
        healthy_skips > degraded_skips && all_downstream_recompute,
        &format!(
            "loose bounds collapse the frontier to hash-memoization behavior \
             (skip yield {healthy_skips:.2} -> {degraded_skips:.2}, everything \
             downstream of the pessimistic op recomputes) while remaining correct — \
             the R4 health metric in action"
        ),
    );
}

/// inv-005 — verified-color skip claims and SLACK BURNING: skip rows
/// carry the verified certificate; applying the plan burns absorbed
/// bounds, so a repeat perturbation finds the slack spent and
/// recomputes.
#[test]
fn inv_005_claims_and_burning() {
    let mut store = Store::new();
    let a = put(&mut store, "a5", 1e-9, 1e-9, b"a5");
    let b = put(&mut store, "b5", 0.0, 0.012, b"b5");
    let edges = vec![Edge { from: a, to: b, sensitivity: 1.0 }];
    let p1 = plan(&store, &edges, &[(a, 0.01)]).expect("plan 1");
    let skip_row = p1
        .rows
        .iter()
        .find(|r| r.contains("\"verdict\":\"skip\""))
        .expect("skip row");
    let claim_ok = skip_row.contains("perturbation absorbed")
        && skip_row.contains("\"color\":\"verified\"")
        && skip_row.contains("interval");
    apply_plan(&mut store, &p1).expect("apply");
    // Slack was 0.012, burned 0.01 → 0.002 left. The same perturbation
    // again cannot be absorbed.
    let p2 = plan(&store, &edges, &[(a, 0.01)]).expect("plan 2");
    let now_recomputes = p2.verdicts.iter().any(|(h, v)| {
        *h == b && matches!(v, Verdict::Recompute { reason: RecomputeReason::BoundExceedsSlack, .. })
    });
    // A smaller perturbation still fits the remainder.
    let p3 = plan(&store, &edges, &[(a, 0.001)]).expect("plan 3");
    let small_still_skips = p3
        .verdicts
        .iter()
        .any(|(h, v)| *h == b && matches!(v, Verdict::Skip { .. }));
    verdict(
        "inv-005",
        claim_ok && now_recomputes && small_still_skips,
        "skip rows carry the verified-color interval claim; applying the plan BURNS \
         slack (0.012 - 0.01), so the repeat 0.01 perturbation recomputes while a \
         0.001 one still fits the remainder — slack is a real, spendable resource",
    );
}
