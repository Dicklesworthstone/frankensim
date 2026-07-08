//! FLYWHEEL CLOSES conformance (the lmp4.18 flagship; runs under the
//! `flywheel-e2e` feature). Acceptance: composed speedup EXCEEDS the
//! best isolated speedup by a stated margin across seeded replays with
//! variance reported; colors propagate (estimated speculation stays
//! estimated through cache → merge → query — laundering blocked); a
//! mid-loop cancellation leaves a consistent trace and a re-run
//! completes clean (G4); whole-loop replay reproduces the trace hash
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
fn fw_002_colors_survive_the_loop_no_laundering() {
    let report = run_loop(&LoopConfig::composed(), 12, 7);
    // Speculation fired, so the headline MUST be estimated (weakest
    // input) — a verified headline over accepted speculation would be
    // laundering.
    assert!(report.accept_rate > 0.0, "speculation was live");
    assert!(
        matches!(report.headline, Color::Estimated { .. }),
        "estimated speculation results keep the headline estimated: {:?}",
        report.headline
    );
    // And the write gate itself refuses the upgrade attempt downstream:
    // a query stage claiming Verified over this headline fails.
    let mut graph = fs_ledger::ColorGraph::new();
    let head = graph.source("loop-headline", report.headline.clone());
    let upgrade = graph.derive(
        "query-answer",
        &[head],
        IntervalOp::Hull,
        Some(Color::Validated {
            regime: ValidityDomain::unconstrained(),
            dataset: "wishful".to_string(),
        }),
        &std::collections::BTreeMap::new(),
        None,
    );
    assert!(
        upgrade.is_err(),
        "the loop's tail cannot upgrade the speculation color without a waiver"
    );
    verdict(
        "fw-002",
        "estimated stays estimated through cache -> merge -> query; the tail upgrade is \
         refused at write time",
    );
}

#[test]
fn fw_003_g5_whole_loop_determinism() {
    let a = run_loop(&LoopConfig::composed(), 12, 99);
    let b = run_loop(&LoopConfig::composed(), 12, 99);
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
    verdict(
        "fw-003",
        "whole-loop replay reproduces bit-equal cost and trace hash; different seeds \
         differ",
    );
}

#[test]
fn fw_004_g4_cancellation_storm_mid_loop() {
    // Cancel at escalating points through the loop: every partial trace
    // must be a clean PREFIX state (consistent events, no partial-stage
    // residue), and the full re-run completes.
    let full = run_loop(&LoopConfig::composed(), 12, 55);
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
            partial.iterations < 12,
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
    }
    // And a fresh full run after the storms is untouched (no leaked
    // state between runs — every run owns its stores).
    let again = run_loop(&LoopConfig::composed(), 12, 55);
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
    println!(
        "{{\"metric\":\"flywheel-telemetry\",\"accept_rate\":{:.2},\"skips\":{},\
         \"merges_ok\":{},\"merges_conflict\":{},\"tombstone_blocks\":{},\"events\":{}}}",
        report.accept_rate,
        report.skips,
        report.merges.0,
        report.merges.1,
        report.tombstone_blocks,
        report.events.len()
    );
    verdict(
        "fw-005",
        "accept rate, skip count, merge verdicts, tombstone blocks, and the event trace \
         all live in one report",
    );
}
