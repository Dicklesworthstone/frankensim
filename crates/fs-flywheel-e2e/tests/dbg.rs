//! Config-sweep smoke: every single-proposal config beats the baseline
//! on cost, and the composed loop beats them all (the quick sanity
//! check behind the measured fw-001 compounding assertion).
#![cfg(feature = "flywheel-e2e")]

use fs_flywheel_e2e::{LoopConfig, run_loop};

#[test]
fn config_sweep_smoke() {
    let base = run_loop(&LoopConfig::baseline(), 10, 11).total_cost;
    let composed = run_loop(&LoopConfig::composed(), 10, 11).total_cost;
    for (name, cfg) in [
        (
            "speculation",
            LoopConfig {
                speculation: true,
                ..LoopConfig::baseline()
            },
        ),
        (
            "recompute",
            LoopConfig {
                recompute: true,
                ..LoopConfig::baseline()
            },
        ),
        (
            "merge",
            LoopConfig {
                merge: true,
                ..LoopConfig::baseline()
            },
        ),
        (
            "tombstones",
            LoopConfig {
                tombstones: true,
                ..LoopConfig::baseline()
            },
        ),
    ] {
        let cost = run_loop(&cfg, 10, 11).total_cost;
        assert!(cost < base, "{name} alone helps: {cost} < {base}");
        assert!(
            composed < cost,
            "composed beats {name} alone: {composed} < {cost}"
        );
    }
}
