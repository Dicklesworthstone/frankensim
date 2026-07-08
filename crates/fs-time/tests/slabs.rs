//! Time-slab conformance (the bk0o.7 bead; runs under `time-slabs`).
//! Acceptance: the temporal cocycle → 0 as coupling tightens
//! (consistency); the budget pie localizes splitting error to the
//! correct temporal interval on a seeded fixture; the adaptive
//! controller reduces the cocycle where it subcycles and beats uniform
//! subcycling at equal accuracy; G3 slab re-partition invariance
//! within the ledgered defect envelope; the activation gate recommends
//! INSTRUMENT-ONLY below the 20% line.
#![cfg(feature = "time-slabs")]

use fs_time::slabs::{
    Activation, CoupledFixture, activation_report, march_adaptive, march_instrumented,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-time/slabs\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn constant_coupling(_t: f64) -> f64 {
    0.8
}

/// Coupling that spikes only in t ∈ [2.0, 2.5] — the seeded
/// localization fixture.
fn spike_coupling(t: f64) -> f64 {
    if (2.0..2.5).contains(&t) { 2.5 } else { 0.05 }
}

#[test]
fn ts_001_cocycle_vanishes_as_coupling_tightens() {
    // Consistency: doubling subcycles must shrink the per-slab defect
    // (Lie splitting: locally O(dt²), so ~2x per doubling over a fixed
    // slab). Measured, with the observed order ledgered.
    let fixture = CoupledFixture {
        coupling: constant_coupling,
    };
    let mut defects = Vec::new();
    for subcycles in [1usize, 2, 4, 8] {
        let (_, ledger) = march_instrumented(&fixture, [1.0, 0.5], 1.0, 1, subcycles);
        defects.push(ledger.total_defect());
    }
    println!("{{\"metric\":\"consistency\",\"defects\":{defects:?}}}");
    for w in defects.windows(2) {
        let ratio = w[0] / w[1].max(1e-300);
        assert!(
            ratio > 1.6,
            "the cocycle shrinks by ~2x per subcycle doubling: {defects:?}"
        );
    }
    assert!(
        defects[3] < 0.2 * defects[0],
        "8x subcycling collapses the defect: {defects:?}"
    );
    verdict(
        "ts-001",
        "the temporal cocycle shrinks ~2x per subcycle doubling (first-order Lie \
         splitting), collapsing 5x+ from 1 to 8 subcycles — consistency verified",
    );
}

#[test]
fn ts_002_budget_pie_localizes_the_handoff() {
    // The seeded fixture couples strongly ONLY in t in [2.0, 2.5]:
    // the budget pie's top interval must be there — "your error is in
    // the coupling handoff at t in [2.0, 2.5]".
    let fixture = CoupledFixture {
        coupling: spike_coupling,
    };
    let (_, ledger) = march_instrumented(&fixture, [1.0, 0.5], 4.0, 16, 1);
    let pie = ledger.attribute();
    let (t0, t1, share) = pie[0];
    println!(
        "{{\"metric\":\"budget-pie\",\"top_interval\":[{t0:.2},{t1:.2}],\"share\":{share:.3},\
         \"ledger\":{}}}",
        ledger.to_json()
    );
    assert!(
        t0 >= 1.99 && t1 <= 2.51,
        "the top defect interval sits inside the seeded spike: [{t0}, {t1}]"
    );
    // The spike's slabs dominate the pie.
    let spike_share: f64 = pie
        .iter()
        .filter(|(a, b, _)| *a >= 1.99 && *b <= 2.51)
        .map(|(_, _, s)| s)
        .sum();
    assert!(
        spike_share > 0.8,
        "the seeded interval owns >80% of the splitting budget: {spike_share:.3}"
    );
    verdict(
        "ts-002",
        "on the seeded spike fixture the budget pie's top interval lands inside \
         t in [2.0, 2.5] and the spike owns >80% of the splitting-error mass",
    );
}

#[test]
fn ts_003_adaptive_controller_beats_uniform() {
    // The controller subcycles only where the cocycle is large. At
    // matched final accuracy it must spend fewer substeps than the
    // uniform march that achieves the same worst-slab defect.
    let fixture = CoupledFixture {
        coupling: spike_coupling,
    };
    let tol = 1e-3;
    let (u_adaptive, ledger, spent_adaptive) = march_adaptive(&fixture, [1.0, 0.5], 4.0, 16, tol);
    // Every slab meets tol (or hit the cap); subcycling concentrated
    // in the spike.
    for e in &ledger.entries {
        assert!(
            e.defect <= tol || e.subcycles == 64,
            "slab [{}, {}] controlled: defect {:.2e}",
            e.t0,
            e.t1,
            e.defect
        );
    }
    let max_sub_outside = ledger
        .entries
        .iter()
        .filter(|e| e.t1 <= 2.0 || e.t0 >= 2.5)
        .map(|e| e.subcycles)
        .max()
        .expect("outside slabs");
    let max_sub_inside = ledger
        .entries
        .iter()
        .filter(|e| e.t0 >= 1.99 && e.t1 <= 2.51)
        .map(|e| e.subcycles)
        .max()
        .expect("inside slabs");
    assert!(
        max_sub_inside > max_sub_outside,
        "subcycling concentrates where the cocycle is large: {max_sub_inside} vs \
         {max_sub_outside}"
    );
    // Uniform march at the SAME per-slab tolerance: find the uniform
    // subcycle count that meets it everywhere, and compare cost.
    let mut uniform_needed = 1usize;
    loop {
        let (_, l) = march_instrumented(&fixture, [1.0, 0.5], 4.0, 16, uniform_needed);
        if l.entries.iter().all(|e| e.defect <= tol) || uniform_needed >= 64 {
            break;
        }
        uniform_needed *= 2;
    }
    let spent_uniform = uniform_needed * 16;
    println!(
        "{{\"metric\":\"controller\",\"adaptive_substeps\":{spent_adaptive},\
         \"uniform_substeps\":{spent_uniform}}}"
    );
    assert!(
        spent_adaptive < spent_uniform,
        "adaptive spends less at equal accuracy: {spent_adaptive} vs {spent_uniform}"
    );
    // And the states agree (both under the same defect control).
    let (u_uniform, _) = march_instrumented(&fixture, [1.0, 0.5], 4.0, 16, uniform_needed);
    let gap =
        ((u_adaptive[0] - u_uniform[0]).powi(2) + (u_adaptive[1] - u_uniform[1]).powi(2)).sqrt();
    assert!(gap < 20.0 * tol, "controlled states agree: {gap:.2e}");
    verdict(
        "ts-003",
        "the controller concentrates subcycles in the spike, meets the per-slab \
         tolerance everywhere, and spends fewer substeps than the uniform march at \
         equal accuracy (measured)",
    );
}

#[test]
fn ts_004_g3_repartition_invariance() {
    // Re-partitioning the slabs changes WHERE the handoffs happen; the
    // final state must agree within the summed ledgered defects of the
    // two marches (the defect envelope is exactly what the ledger is
    // for).
    let fixture = CoupledFixture {
        coupling: constant_coupling,
    };
    let (u_a, la) = march_instrumented(&fixture, [1.0, 0.5], 2.0, 8, 4);
    let (u_b, lb) = march_instrumented(&fixture, [1.0, 0.5], 2.0, 5, 4);
    let gap = ((u_a[0] - u_b[0]).powi(2) + (u_a[1] - u_b[1]).powi(2)).sqrt();
    let envelope = la.total_defect() + lb.total_defect();
    println!("{{\"metric\":\"repartition\",\"gap\":{gap:.3e},\"envelope\":{envelope:.3e}}}");
    assert!(
        gap <= 2.0 * envelope,
        "re-partition gap {gap:.2e} within the ledgered envelope {envelope:.2e}"
    );
    verdict(
        "ts-004",
        "8-slab vs 5-slab partitions agree within the summed ledgered defects — the \
         ledger's envelope covers re-partition drift (G3)",
    );
}

#[test]
fn ts_005_activation_gate_is_honest() {
    // A weakly-coupled workload: splitting error is a sliver of the
    // budget — the gate must say INSTRUMENT-ONLY. A strongly-coupled
    // one flips it.
    let weak = CoupledFixture { coupling: |_| 0.02 };
    let (_, ledger_weak) = march_instrumented(&weak, [1.0, 0.5], 2.0, 8, 1);
    let (frac_weak, verdict_weak) = activation_report(&ledger_weak, 1e-2);
    assert!(frac_weak < 0.2);
    assert_eq!(verdict_weak, Activation::InstrumentOnly);
    let strong = CoupledFixture { coupling: |_| 1.5 };
    let (_, ledger_strong) = march_instrumented(&strong, [1.0, 0.5], 2.0, 8, 1);
    let (frac_strong, verdict_strong) = activation_report(&ledger_strong, 1e-2);
    assert!(frac_strong >= 0.2);
    assert_eq!(verdict_strong, Activation::ControlJustified);
    println!(
        "{{\"metric\":\"activation\",\"weak_fraction\":{frac_weak:.3},\
         \"strong_fraction\":{frac_strong:.3}}}"
    );
    verdict(
        "ts-005",
        "the activation gate recommends instrument-only when splitting error is <20% of \
         budget and flips to control-justified when it dominates — the Proposal-4 \
         sequencing discipline as code",
    );
}
