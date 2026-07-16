//! FLYWHEEL CLOSES conformance (the lmp4.18 flagship; runs under the
//! `flywheel-e2e` feature). Acceptance: composed speedup EXCEEDS the
//! best isolated speedup by a stated margin across seeded replays with
//! variance reported; colors propagate through one retained graph
//! (estimated speculation stays estimated through cache → merge → query;
//! re-rooting and laundering are blocked); a mid-loop cancellation leaves
//! consistent event and evidence prefixes and a re-run completes clean (G4);
//! whole-loop replay reproduces the input-bound report and graph commitment
//! bit-for-bit (G5).
#![cfg(feature = "flywheel-e2e")]

use fs_evidence::{Color, IntervalOp, ValidityDomain};
use fs_flywheel_e2e::{LoopConfig, run_loop, speedups};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"flywheel-e2e\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The super-additivity margin: composed must beat the best isolated
/// speedup by at least this factor (review round 3: a stated margin,
/// not a single-run "exceeds").
const MARGIN: f64 = 1.15;

/// Test-only exact waiver capability. The signature is a public content hash,
/// so this exercises waiver taint and request binding without claiming
/// cryptographic authentication.
struct ExactFixtureWaiverVerifier;

impl fs_ledger::WaiverVerifier for ExactFixtureWaiverVerifier {
    fn verify(&self, key_id: &str, payload: &[u8], signature: &[u8]) -> fs_ledger::PolicyDecision {
        let policy = fs_ledger::hash_bytes(b"fs-flywheel-e2e/headline-waiver-policy/v1");
        let expected_signature = fs_ledger::hash_bytes(payload);
        let accepted =
            key_id == "flywheel-headline-fixture-key" && signature == expected_signature.as_bytes();
        if accepted {
            fs_ledger::PolicyDecision::accept(policy)
        } else {
            fs_ledger::PolicyDecision::reject(policy)
        }
    }
}

fn signed_fixture_source_grant(name: &str, color: &Color) -> fs_ledger::WaiverGrant {
    let mut grant = fs_ledger::WaiverGrant {
        annotation: fs_ledger::Waiver {
            id: "WVR-FLYWHEEL-HEADLINE-001".to_string(),
            signer: "flywheel-e2e-fixture".to_string(),
            reason: "exercise fail-closed scientific headline access".to_string(),
        },
        key_id: "flywheel-headline-fixture-key".to_string(),
        scope: fs_ledger::WAIVER_SCOPE_SOURCE_COLOR.to_string(),
        node_name: name.to_string(),
        claimed_color: color.canonical_bytes(),
        parent_hashes: Vec::new(),
        expires_day: 400,
        signature: Vec::new(),
    };
    grant.signature = fs_ledger::hash_bytes(&grant.signing_payload_source())
        .as_bytes()
        .to_vec();
    grant
}

#[test]
fn fw_001_compounding_measured_with_margin_and_variance() {
    let seeds = [11u64, 23, 47, 89, 173];
    let mut ratios = Vec::new();
    for &seed in &seeds {
        let (isolated, composed) = speedups(12, seed);
        let best = isolated.values().fold(0.0f64, |a, &b| a.max(b));
        let ratio = composed / best;
        assert!(
            ratio > MARGIN,
            "seed {seed}: composed {composed:.2}x vs best isolated {best:.2}x — the \
             loop must COMPOUND (ratio {ratio:.2} <= margin {MARGIN})"
        );
        ratios.push(ratio);
        println!(
            "{{\"metric\":\"compounding\",\"seed\":{seed},\"isolated\":{isolated:?},\
             \"composed\":{composed:.3},\"ratio\":{ratio:.3}}}"
        );
    }
    // Across-replay variance, reported (and sane).
    let mean = ratios.iter().sum::<f64>() / ratios.len() as f64;
    let var = ratios.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / ratios.len() as f64;
    let cv = var.sqrt() / mean;
    println!(
        "{{\"metric\":\"compounding-variance\",\"mean\":{mean:.3},\"variance\":{var:.5},\
         \"cv\":{cv:.3},\"replays\":5}}"
    );
    assert!(
        cv < 0.25,
        "the across-replay coefficient of variation is bounded: {cv}"
    );
    verdict(
        "fw-001",
        "5 seeded replays: composed speedup beats the best isolated by >1.15x every \
         time; variance reported",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // one auditable retained-lineage refusal fixture
fn fw_002_colors_survive_the_loop_no_laundering() {
    let mut report = run_loop(&LoopConfig::composed(), 12, 7);
    // Speculation fired, so the headline MUST be estimated (weakest
    // input) — a verified headline over accepted speculation would be
    // laundering.
    assert!(report.accept_rate > 0.0, "speculation was live");
    let headline = report
        .headline()
        .expect("a run-produced headline node resolves")
        .clone();
    assert!(
        matches!(headline, Color::Estimated { .. }),
        "estimated speculation results keep the headline estimated: {headline:?}"
    );
    let Color::Estimated { estimator, .. } = &headline else {
        unreachable!("the assertion above established the headline rank");
    };
    assert!(
        estimator.starts_with("derived:v2:"),
        "the first accepted heterogeneous step already has a lineage-only identity: {estimator}"
    );
    report
        .color_graph
        .verify_replay()
        .expect("the retained graph replays from its authenticated baseline");
    let baseline = report
        .color_graph
        .node(0)
        .expect("the baseline source is retained");
    assert!(
        matches!(
            baseline.origin(),
            Some(fs_ledger::SourceOrigin::Certificate { .. })
        ) && baseline.origin_policy_fingerprint().is_some(),
        "the verified baseline retains its typed certificate and admitting policy"
    );

    let proposer_identities: std::collections::BTreeSet<&str> = report
        .color_graph
        .nodes()
        .iter()
        .filter_map(|node| match node.declared_color_unverified() {
            Color::Estimated { estimator, .. } if node.parents().is_empty() => {
                Some(estimator.as_str())
            }
            _ => None,
        })
        .collect();
    assert!(
        proposer_identities.contains("wedge-proposer-v1/agent-0/dataset-agent-0")
            && proposer_identities.contains("wedge-proposer-v1/agent-1/dataset-agent-1")
            && proposer_identities.contains("wedge-validation-probe-v1/agent-0/dataset-holdout-1")
            && proposer_identities.contains("wedge-validation-probe-v1/agent-1/dataset-holdout-0"),
        "heterogeneous proposer and holdout identities survive as source leaves: \
         {proposer_identities:?}"
    );

    // A derived estimator identity is reserved for nodes with lineage. It
    // cannot be re-rooted as a fresh source, whether copied from this report
    // or forged directly.
    let reroot = report
        .color_graph
        .source("illicit-reroot", headline.clone());
    assert!(
        matches!(
            reroot,
            Err(fs_ledger::ColorWriteError::InvalidEstimatedSource {
                field: "estimator",
                why: "derived-identity-requires-lineage",
            })
        ),
        "a derived headline is refused specifically because re-rooting loses lineage"
    );
    let forged_reserved = report.color_graph.source(
        "forged-derived-source",
        Color::Estimated {
            estimator: "derived:v2:estimators:forged".to_string(),
            dispersion: 0.1,
        },
    );
    assert!(
        matches!(
            forged_reserved,
            Err(fs_ledger::ColorWriteError::InvalidEstimatedSource {
                field: "estimator",
                why: "derived-identity-requires-lineage",
            })
        ),
        "the reserved derived namespace is refused specifically without retained lineage"
    );

    // The downstream query acts on the retained headline node, rather than
    // manufacturing a new leaf, and the write gate refuses its upgrade.
    let upgrade = report.color_graph.derive(
        "query-answer",
        &[report.headline_node],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained().with("reynolds", 1.0e4, 1.0e6),
            dataset: "wishful".to_string(),
        }),
        &std::collections::BTreeMap::new(),
        None,
    );
    assert!(
        matches!(
            &upgrade,
            Err(fs_ledger::ColorWriteError::LaunderingRefused {
                claimed: fs_evidence::ColorRank::Validated,
                derived: fs_evidence::ColorRank::Estimated,
                offending_parents,
            }) if offending_parents == &[report.headline_node]
        ),
        "a structurally valid Validated claim is refused specifically as laundering: {upgrade:?}"
    );

    // Even an authenticated waiver is an explicit human-responsibility door,
    // not ordinary scientific evidence. Pointing the report at a waived node
    // must therefore fail closed through headline(), as must an unknown id.
    let waived_name = "waived-headline-fixture";
    let waived_color = Color::Validated {
        regime: ValidityDomain::unconstrained().with("reynolds", 1.0e4, 1.0e6),
        dataset: "fixture-waiver-only".to_string(),
    };
    let waived_node = report
        .color_graph
        .source_waived(
            waived_name,
            waived_color.clone(),
            signed_fixture_source_grant(waived_name, &waived_color),
            &ExactFixtureWaiverVerifier,
            200,
        )
        .expect("the exact fixture waiver is admitted and remains visibly tainted");
    assert!(
        report
            .color_graph
            .node(waived_node)
            .is_some_and(fs_ledger::ColorNode::depends_on_waiver),
        "the regression fixture itself must remain waiver-tainted"
    );
    report.headline_node = waived_node;
    assert!(
        report.headline().is_none(),
        "an authenticated waiver cannot silently become a scientific headline"
    );
    report.headline_node = u64::MAX;
    assert!(
        report.headline().is_none(),
        "an unresolved node id cannot become a scientific headline"
    );
    verdict(
        "fw-002",
        "authenticated baseline + heterogeneous proposal leaves retain one lineage; re-rooting, \
         reserved-source forgery, tail upgrades, and non-scientific headlines are refused",
    );
}

#[test]
fn fw_003_g5_whole_loop_determinism() {
    let a = run_loop(&LoopConfig::composed(), 12, 99);
    let b = run_loop(&LoopConfig::composed(), 12, 99);
    assert_eq!(a, b, "every semantic report and evidence field replays");
    assert_eq!(
        a.total_cost.to_bits(),
        b.total_cost.to_bits(),
        "bit-equal costs"
    );
    assert_eq!(a.trace_hash(), b.trace_hash(), "identical trace hashes");
    assert_eq!(a.events.len(), b.events.len());
    // A different seed genuinely changes the trace (the hash is not
    // vacuous).
    let c = run_loop(&LoopConfig::composed(), 12, 100);
    assert_ne!(a.trace_hash(), c.trace_hash(), "the hash sees the workload");

    // Field-sensitivity does not need the large workload above. A zero-step
    // report still carries the authenticated baseline graph and every report
    // field, keeping this exhaustive mutation battery cheap.
    let hash_fixture = run_loop(&LoopConfig::baseline(), 0, 991);
    let base_hash = hash_fixture.trace_hash();
    let mut negative_zero = run_loop(&LoopConfig::baseline(), 0, 991);
    negative_zero.total_cost = -0.0;
    assert_ne!(
        hash_fixture, negative_zero,
        "report equality preserves signed-zero bit identity"
    );
    assert_ne!(
        base_hash,
        negative_zero.trace_hash(),
        "the G5 hash preserves signed-zero bit identity"
    );
    let assert_sensitive = |label: &str, mutate: &dyn Fn(&mut fs_flywheel_e2e::LoopReport)| {
        let mut changed = run_loop(&LoopConfig::baseline(), 0, 991);
        mutate(&mut changed);
        assert_ne!(base_hash, changed.trace_hash(), "the G5 hash binds {label}");
    };
    assert_sensitive("speculation config", &|r| r.config.speculation = true);
    assert_sensitive("recompute config", &|r| r.config.recompute = true);
    assert_sensitive("merge config", &|r| r.config.merge = true);
    assert_sensitive("tombstone config", &|r| r.config.tombstones = true);
    assert_sensitive("cancellation config", &|r| {
        r.config.cancel_after_stages = Some(1);
    });
    assert_sensitive("requested iterations", &|r| r.requested_iterations += 1);
    assert_sensitive("seed", &|r| r.seed += 1);
    assert_sensitive("total_cost", &|r| r.total_cost += 1.0);
    assert_sensitive("iterations", &|r| r.iterations += 1);
    assert_sensitive("cancelled", &|r| r.cancelled = !r.cancelled);
    assert_sensitive("events", &|r| r.events.push("tamper".to_string()));
    assert_sensitive("accept_rate", &|r| r.accept_rate += 0.125);
    assert_sensitive("skips", &|r| r.skips += 1);
    assert_sensitive("resolved merges", &|r| r.merges.0 += 1);
    assert_sensitive("conflicted merges", &|r| r.merges.1 += 1);
    assert_sensitive("tombstone blocks", &|r| r.tombstone_blocks += 1);
    assert_sensitive("headline node", &|r| r.headline_node = u64::MAX);
    assert_sensitive("color graph", &|r| {
        r.color_graph
            .source(
                "post-report-tamper",
                Color::Estimated {
                    estimator: "post-report-tamper-v1".to_string(),
                    dispersion: f64::INFINITY,
                },
            )
            .expect("the graph mutation is structurally valid");
    });
    verdict(
        "fw-003",
        "whole-loop replay reproduces the complete report and graph; the G5 hash binds every \
         semantic report field and all retained evidence",
    );
}

#[test]
fn fw_004_g4_cancellation_storm_mid_loop() {
    // Cancel at escalating points through the loop: every partial trace
    // must be a clean PREFIX state (consistent events, no partial-stage
    // residue), and the full re-run completes.
    let full = run_loop(&LoopConfig::composed(), 12, 55);
    assert_eq!(
        full.iterations, full.requested_iterations,
        "every non-cancelled terminal outcome completes a design iteration"
    );
    for cancel_at in [1usize, 3, 7, 15, 31] {
        let partial = run_loop(
            &LoopConfig {
                cancel_after_stages: Some(cancel_at),
                ..LoopConfig::composed()
            },
            12,
            55,
        );
        assert!(partial.cancelled, "the storm fired at {cancel_at}");
        assert!(
            partial.iterations < full.iterations,
            "cancellation genuinely interrupted the loop"
        );
        // The partial trace is a strict prefix of the full trace: no
        // torn stages, no divergence before the cut.
        assert!(
            partial.events.len() <= full.events.len(),
            "no extra residue events"
        );
        for (p, f) in partial.events.iter().zip(&full.events) {
            assert_eq!(p, f, "the partial trace is a clean prefix");
        }
        partial
            .color_graph
            .verify_replay()
            .expect("a cancelled run retains a replayable evidence prefix");
        assert!(
            full.color_graph
                .rows()
                .starts_with(partial.color_graph.rows()),
            "cancellation leaves only a canonical graph-row prefix"
        );
        assert!(
            partial.color_graph.nodes().len() <= full.color_graph.nodes().len(),
            "a partial graph cannot contain nodes beyond the full replay"
        );
        assert!(
            full.color_graph
                .nodes()
                .iter()
                .map(fs_ledger::ColorNode::hash)
                .zip(
                    partial
                        .color_graph
                        .nodes()
                        .iter()
                        .map(fs_ledger::ColorNode::hash)
                )
                .all(|(full, partial)| full == partial),
            "cancellation leaves only a provenance-node prefix"
        );
        assert!(
            partial.headline().is_some(),
            "the partial headline always resolves inside its retained graph"
        );
    }
    // And a fresh full run after the storms is untouched (no leaked
    // state between runs — every run owns its stores).
    let again = run_loop(&LoopConfig::composed(), 12, 55);
    assert_eq!(full, again, "no semantic or evidence residue crosses runs");
    assert_eq!(
        full.trace_hash(),
        again.trace_hash(),
        "no cross-run residue"
    );
    verdict(
        "fw-004",
        "5 cancel points: every partial trace is a clean prefix; re-runs reproduce the \
         full trace exactly",
    );
}

#[test]
fn fw_005_telemetry_completeness() {
    let report = run_loop(&LoopConfig::composed(), 12, 7);
    // The whole flywheel's telemetry in one trace: every dial moved.
    assert!(
        report.accept_rate > 0.3,
        "accept rate: {}",
        report.accept_rate
    );
    assert!(report.skips > 0, "recompute skips: {}", report.skips);
    assert!(report.merges.0 > 0, "resolved merges: {:?}", report.merges);
    assert!(
        report.tombstone_blocks > 0,
        "tombstone blocks: {}",
        report.tombstone_blocks
    );
    assert!(!report.events.is_empty());
    assert!(
        report.color_graph.nodes().len() > 1,
        "accepted evidence lineage is retained in telemetry"
    );
    println!(
        "{{\"metric\":\"flywheel-telemetry\",\"accept_rate\":{:.2},\"skips\":{},\
         \"merges_ok\":{},\"merges_conflict\":{},\"tombstone_blocks\":{},\"events\":{},\
         \"evidence_nodes\":{}}}",
        report.accept_rate,
        report.skips,
        report.merges.0,
        report.merges.1,
        report.tombstone_blocks,
        report.events.len(),
        report.color_graph.nodes().len()
    );
    verdict(
        "fw-005",
        "accept rate, skip count, merge verdicts, tombstone blocks, event trace, and retained \
         evidence graph all live in one report",
    );
}

fn events_with_stage<'a>(report: &'a fs_flywheel_e2e::LoopReport, stage: &str) -> Vec<&'a str> {
    report
        .events
        .iter()
        .filter(|event| event.contains(stage))
        .map(String::as_str)
        .collect()
}

fn solve_event_costs(report: &fs_flywheel_e2e::LoopReport) -> Vec<f64> {
    events_with_stage(report, "stage=solve")
        .into_iter()
        .map(|event| {
            event
                .rsplit_once("cost=")
                .expect("solve event carries a cost")
                .1
                .parse::<f64>()
                .expect("solve event cost is finite decimal data")
        })
        .collect()
}

#[test]
fn fw_006_feature_toggles_do_not_perturb_common_workload_draws() {
    let anchor = LoopConfig {
        tombstones: true,
        ..LoopConfig::baseline()
    };
    let baseline = run_loop(&anchor, 18, 0x5eed);
    let baseline_candidates = events_with_stage(&baseline, "stage=tombstone");

    for config in [
        LoopConfig {
            speculation: true,
            ..anchor
        },
        LoopConfig {
            recompute: true,
            ..anchor
        },
        LoopConfig {
            merge: true,
            ..anchor
        },
        LoopConfig::composed(),
    ] {
        let toggled = run_loop(&config, 18, 0x5eed);
        assert_eq!(
            events_with_stage(&toggled, "stage=tombstone"),
            baseline_candidates,
            "proposal toggles cannot shift candidate-velocity draws"
        );
    }

    let merge_anchor = run_loop(
        &LoopConfig {
            merge: true,
            tombstones: true,
            ..LoopConfig::baseline()
        },
        18,
        0x5eed,
    );
    let merge_verdicts = events_with_stage(&merge_anchor, "stage=merge");
    for config in [
        LoopConfig {
            speculation: true,
            merge: true,
            tombstones: true,
            ..LoopConfig::baseline()
        },
        LoopConfig {
            recompute: true,
            merge: true,
            tombstones: true,
            ..LoopConfig::baseline()
        },
        LoopConfig::composed(),
    ] {
        assert_eq!(
            events_with_stage(&run_loop(&config, 18, 0x5eed), "stage=merge"),
            merge_verdicts,
            "unrelated toggles cannot shift merge-taint or gauge draws"
        );
    }
}

#[test]
fn fw_007_cancellation_charges_every_completed_branch_prefix() {
    let parallel_one = run_loop(
        &LoopConfig {
            cancel_after_stages: Some(2),
            ..LoopConfig::composed()
        },
        1,
        55,
    );
    let one_cost = solve_event_costs(&parallel_one);
    assert!(parallel_one.cancelled);
    assert_eq!(one_cost.len(), 1);
    assert!((parallel_one.total_cost - one_cost[0]).abs() < 1e-12);

    let parallel_two = run_loop(
        &LoopConfig {
            cancel_after_stages: Some(3),
            ..LoopConfig::composed()
        },
        1,
        55,
    );
    let two_costs = solve_event_costs(&parallel_two);
    assert!(parallel_two.cancelled);
    assert_eq!(two_costs.len(), 2);
    assert!(
        (parallel_two.total_cost - two_costs[0].max(two_costs[1])).abs() < 1e-12,
        "parallel work settles to the completed critical path before merge"
    );

    let serialized_two = run_loop(
        &LoopConfig {
            tombstones: true,
            cancel_after_stages: Some(3),
            ..LoopConfig::baseline()
        },
        1,
        55,
    );
    let serialized_costs = solve_event_costs(&serialized_two);
    assert!(serialized_two.cancelled);
    assert_eq!(serialized_costs.len(), 2);
    assert!(
        (serialized_two.total_cost - (serialized_costs[0] + serialized_costs[1])).abs() < 1e-12,
        "serialized work charges both completed branches before merge"
    );
}

#[test]
fn fw_008_zero_work_speedups_are_neutral_and_finite() {
    let (isolated, composed) = fs_flywheel_e2e::speedups(0, 0x5eed);
    assert_eq!(isolated.len(), 4);
    assert!(
        isolated
            .values()
            .all(|ratio| ratio.to_bits() == 1.0_f64.to_bits())
    );
    assert_eq!(composed.to_bits(), 1.0_f64.to_bits());
}

#[test]
fn fw_010_admitted_headline_carries_lineage_and_refuses_estimated() {
    use fs_flywheel_e2e::HeadlineAdmissionRefusal;
    use fs_ledger::colors::ColorAdmissionRefusal;

    // No speculation: the headline stays positive and must convert to
    // admitted evidence whose lineage names the exact headline node.
    let baseline = run_loop(&LoopConfig::baseline(), 12, 7);
    let headline = baseline
        .headline()
        .expect("baseline run resolves a scientific headline")
        .clone();
    assert!(
        !matches!(headline, Color::Estimated { .. }),
        "baseline headline is positive: {headline:?}"
    );
    let admitted = baseline
        .admitted_headline()
        .expect("a positive replayable headline admits");
    assert_eq!(admitted.admitted_color(), &headline);
    assert_eq!(
        admitted.receipt().node_hash(),
        baseline
            .color_graph
            .node(baseline.headline_node)
            .expect("headline node exists")
            .hash(),
        "the admission lineage names the exact headline node"
    );

    // Speculation fired: the headline is estimated (fw_002) and admission
    // must refuse NotPositive rather than laundering the rank upward.
    let composed = run_loop(&LoopConfig::composed(), 12, 7);
    let refusal = composed
        .admitted_headline()
        .expect_err("an estimated headline never acquires scientific admission");
    assert!(
        matches!(
            refusal,
            HeadlineAdmissionRefusal::Mint(ColorAdmissionRefusal::NotPositive { .. })
        ),
        "estimated headline refuses at the minting gate: {refusal}"
    );
}
