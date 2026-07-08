//! Persistence-penalty conformance (the 7tv.15 bead; runs under
//! `moonshot-topo-persistence`). Acceptance: penalty zero iff the
//! diagram matches the target up to τ (G0); seeded violations
//! (enclosed void in a bracket, spurious island) are driven to target
//! topology by attribution-guided descent; attribution directions are
//! verified by perturbation; sub-persistence noise is ignored while
//! real features are enforced; the [M] promotion gate — the graded,
//! localized penalty beats the connected-component-labeling heuristic
//! (which has no gradient) on the fixture suite, measured and
//! ledgered; PH stability under field perturbation.
#![cfg(feature = "moonshot-topo-persistence")]

use fs_topo::cubical::VoxelField;
use fs_topo::penalty::{
    TopoSpec, apply_attribution_step, enclosed_voids, evaluate, heuristic_cc_penalty,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-topo/penalty\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A solid n³ bracket: value −1 inside a centered box, +1 outside
/// (solid = value < 0).
fn bracket(n: u32) -> VoxelField {
    let mut values = vec![1.0f64; (n * n * n) as usize];
    for z in 2..n - 2 {
        for y in 2..n - 2 {
            for x in 2..n - 2 {
                values[((z * n + y) * n + x) as usize] = -1.0;
            }
        }
    }
    VoxelField {
        dims: [n, n, n],
        values,
        h: 1.0,
    }
}

fn carve(field: &mut VoxelField, lo: [u32; 3], hi: [u32; 3], value: f64) {
    let [nx, ny, _] = field.dims;
    for z in lo[2]..hi[2] {
        for y in lo[1]..hi[1] {
            for x in lo[0]..hi[0] {
                field.values[((z * ny + y) * nx + x) as usize] = value;
            }
        }
    }
}

fn spec() -> TopoSpec {
    TopoSpec {
        components: 1,
        tunnels: 0,
        enclosed_voids: 0,
        tau: 0.3,
        level: 0.0,
    }
}

#[test]
fn tp_001_g0_zero_iff_target_matched() {
    // A clean single-component bracket meets (1, 0, 0): penalty 0.
    let clean = bracket(16);
    let report = evaluate(&clean, &spec());
    assert_eq!(report.betti, (1, 0, 0));
    assert!(report.total == 0.0, "clean bracket: {}", report.total);
    assert!(report.attributions.is_empty());
    // Seed an ENCLOSED VOID (castability violation): penalty > 0 with
    // a fill-direction attribution exactly on the void.
    let mut voided = bracket(16);
    carve(&mut voided, [6, 6, 6], [10, 10, 10], 1.0);
    let report = evaluate(&voided, &spec());
    assert!(report.total > 0.0, "void detected");
    let att = report
        .attributions
        .iter()
        .find(|a| a.channel == "enclosed-void")
        .expect("void attribution");
    assert_eq!(att.voxels.len(), 64, "the 4x4x4 void is fully attributed");
    assert!(att.direction > 0.0, "fill direction");
    // Seed a SPURIOUS ISLAND (excess component): carve-direction
    // attribution on the island.
    let mut island = bracket(16);
    carve(&mut island, [0, 0, 0], [2, 2, 2], -0.9);
    let report = evaluate(&island, &spec());
    assert!(report.total > 0.0);
    assert!(
        report
            .attributions
            .iter()
            .any(|a| a.channel == "excess-component" && a.direction < 0.0),
        "island attributed for carving: {:?}",
        report
            .attributions
            .iter()
            .map(|a| a.channel)
            .collect::<Vec<_>>()
    );
    verdict(
        "tp-001",
        "penalty is exactly zero on the compliant bracket; a seeded 4x4x4 void and a \
         spurious island each produce positive penalty with correctly-signed, exactly \
         localized attributions",
    );
}

#[test]
fn tp_002_attribution_perturbation_directions() {
    // Moving density WHERE the attribution says reduces the penalty;
    // moving it elsewhere does not.
    let mut voided = bracket(16);
    carve(&mut voided, [6, 6, 6], [10, 10, 10], 1.0);
    let report = evaluate(&voided, &spec());
    let before = report.total;
    // Fill AT the attribution.
    let mut guided = voided.clone();
    let touched = apply_attribution_step(&mut guided, &report, 0.6);
    assert!(touched >= 64);
    let after = evaluate(&guided, &spec()).total;
    assert!(
        after < before,
        "attributed step reduces the penalty: {before} -> {after}"
    );
    // Perturb AWAY from the attribution (a corner of the solid): the
    // penalty must not decrease.
    let mut misguided = voided.clone();
    carve(&mut misguided, [3, 3, 3], [5, 5, 5], -1.6);
    let off = evaluate(&misguided, &spec()).total;
    assert!(
        off >= before - 1e-12,
        "off-target perturbation does not reduce the penalty: {before} -> {off}"
    );
    verdict(
        "tp-002",
        "perturbation test: the attributed fill strictly reduces the void penalty; an \
         off-target perturbation of equal magnitude does not",
    );
}

#[test]
fn tp_003_threshold_is_the_feature_size_floor() {
    // A SHALLOW dimple (depth < tau) is noise: ignored. A DEEP void
    // (depth > tau) is real: enforced.
    let mut shallow = bracket(16);
    carve(&mut shallow, [6, 6, 6], [10, 10, 10], 0.2); // depth 0.2 < 0.3
    let report = evaluate(&shallow, &spec());
    assert!(
        report.void_depths.is_empty() && report.total == 0.0,
        "sub-tau dimple ignored: {:?}",
        report.void_depths
    );
    let mut deep = bracket(16);
    carve(&mut deep, [6, 6, 6], [10, 10, 10], 0.5); // depth 0.5 > 0.3
    let report = evaluate(&deep, &spec());
    assert!(
        report.total > 0.0 && report.void_depths == vec![0.5],
        "supra-tau void enforced: {:?}",
        report.void_depths
    );
    // PH STABILITY (inherited): a small uniform perturbation moves the
    // penalty by at most a comparable amount.
    let mut jittered = deep.clone();
    let mut lcg = 0x1234u64;
    for v in &mut jittered.values {
        lcg = lcg
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *v += (((lcg >> 11) as f64) / (1u64 << 53) as f64 - 0.5) * 0.02;
    }
    let (a, b) = (
        evaluate(&deep, &spec()).total,
        evaluate(&jittered, &spec()).total,
    );
    assert!(
        (a - b).abs() < 0.15 * a.max(1.0),
        "stability: {a} vs {b} under 0.01-scale jitter"
    );
    verdict(
        "tp-003",
        "0.2-deep dimple ignored, 0.5-deep void enforced at tau = 0.3 (the feature-size \
         floor); the penalty is stable under small field jitter",
    );
}

#[test]
fn tp_004_m_gate_beats_the_cc_heuristic() {
    // THE PROMOTION GATE: attribution-guided descent drives the seeded
    // bracket to castable topology; the CC-labeling heuristic has ZERO
    // gradient until the void vanishes, so the same step budget leaves
    // it stuck (its penalty is flat and it offers no direction).
    let mut voided = bracket(16);
    carve(&mut voided, [6, 6, 6], [10, 10, 10], 1.0);
    let s = spec();
    // Persistence route: iterate attribution steps.
    let mut field = voided.clone();
    let mut steps = 0usize;
    let mut trace = Vec::new();
    for _ in 0..12 {
        let report = evaluate(&field, &s);
        trace.push((report.total * 1e3).round() / 1e3);
        if report.total == 0.0 {
            break;
        }
        let _ = apply_attribution_step(&mut field, &report, 0.6);
        steps += 1;
    }
    let final_pen = evaluate(&field, &s);
    println!(
        "{{\"metric\":\"m-gate\",\"persistence_steps\":{steps},\
         \"persistence_final\":{:.4},\"trace\":{trace:?}}}",
        final_pen.total
    );
    assert!(
        final_pen.total == 0.0 && final_pen.betti == (1, 0, 0),
        "attribution descent reaches castable topology in {steps} steps"
    );
    // Heuristic route: the CC penalty gives count-only feedback — flat
    // over any local density change that does not remove the void
    // outright. Verify the flatness that makes it un-followable.
    let h0 = heuristic_cc_penalty(&voided, &s);
    let mut nudged = voided.clone();
    carve(&mut nudged, [7, 7, 7], [9, 9, 9], 0.4); // partial fill
    let h1 = heuristic_cc_penalty(&nudged, &s);
    // Counts are integers cast to f64: exact comparison is meaningful,
    // spelled with an epsilon to satisfy the float-eq lint.
    assert!(
        (h0 - h1).abs() < 1e-12 && h0 > 0.0,
        "the heuristic is FLAT under partial progress ({h0} == {h1}) — no gradient to \
         follow; the persistence penalty is graded (trace above)"
    );
    verdict(
        "tp-004",
        "the [M] gate, measured: attribution-guided persistence descent reaches castable \
         topology in a handful of steps while the CC-labeling heuristic is provably flat \
         under partial progress",
    );
}

#[test]
fn tp_005_duality_route_and_tunnel_counting() {
    // The enclosed-void finder agrees with betti's H2 on the same
    // field (the duality route's cross-check).
    let mut voided = bracket(12);
    carve(&mut voided, [5, 5, 5], [8, 8, 8], 1.0);
    let voids = enclosed_voids(&voided, 0.0);
    let b = fs_topo::cubical::betti(&voided, 0.0);
    assert_eq!(
        voids.len() as u32,
        b.2,
        "duality cross-check: |voids| == b2"
    );
    // Tunnel channel: a through-hole (donut) counts b1 = 1; the spec
    // demanding 1 tunnel is satisfied, demanding 0 penalizes.
    let mut donut = bracket(16);
    carve(&mut donut, [7, 7, 0], [9, 9, 16], 1.0); // full-z channel
    let b = fs_topo::cubical::betti(&donut, 0.0);
    assert_eq!(b.1, 1, "the through-channel is one tunnel");
    let mut want_tunnel = spec();
    want_tunnel.tunnels = 1;
    let ok = evaluate(&donut, &want_tunnel);
    assert!(ok.total == 0.0, "routing spec satisfied: {}", ok.total);
    let forbid = evaluate(&donut, &spec());
    assert!(
        forbid.total > 0.0
            && forbid
                .attributions
                .iter()
                .any(|a| a.channel == "tunnel-mismatch"),
        "excess tunnel penalized (counted, not localized — the documented no-claim)"
    );
    verdict(
        "tp-005",
        "enclosed-void finder cross-checks against betti H2 (duality); a through-channel \
         satisfies a 1-tunnel routing spec and is penalized under a 0-tunnel spec",
    );
}
