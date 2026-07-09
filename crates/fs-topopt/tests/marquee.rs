//! CutFEM-octree marquee conformance (the b7d0 bead; runs under
//! `cutfem-marquee`). Acceptance: the objective improves at fixed
//! volume with ZERO mesh rebuilds (the marquee property — asserted
//! from the run log); topology EVOLVES (void components change or
//! boundary-side nodes flip) without any remesh; DWR-driven splits
//! concentrate near the design boundary; the benchmark envelope is
//! recorded as evidence with per-iteration timing (the
//! interactive-cadence measurement, debug-build labeled); the
//! medial-axis thickness oracle audits the optimized geometry's
//! minimum length scale.
#![cfg(feature = "cutfem-marquee")]

use fs_topopt::marquee::{DensityDesign, run_marquee};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-topopt/marquee\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// The volume-to-point start: a solid plate with two seeded cooling
/// voids (the classic heat-sink layout germ).
fn seeded_design() -> DensityDesign {
    let n = 17;
    let mut d = DensityDesign::uniform(n, 0.75);
    #[allow(clippy::cast_precision_loss)]
    let scale = (n - 1) as f64;
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (i as f64 / scale, j as f64 / scale);
            for &(cx, cy) in &[(0.3, 0.5), (0.7, 0.5)] {
                if (x - cx).powi(2) + (y - cy).powi(2) < 0.012 {
                    d.rho[j * n + i] = 0.1;
                }
            }
        }
    }
    d
}

#[test]
fn tm_001_improves_with_zero_rebuilds() {
    // The marquee property is zero mesh rebuilds: the background
    // quadtree is built once, then the configured adaptation path may
    // refine it by octree splits. Split concentration is tm-003's
    // executable evidence.
    let report = run_marquee(seeded_design(), 4, 6, 6, 4).expect("marquee runs");
    let first = report.iterations.first().expect("first").compliance;
    let last = report.iterations.last().expect("last").compliance;
    println!(
        "{{\"metric\":\"marquee-objective\",\"first\":{first:.5},\"last\":{last:.5},\
         \"total_splits\":{},\"total_rebuilds\":{}}}",
        report.total_splits, report.total_rebuilds
    );
    assert!(last < first, "compliance improves: {first:.5} -> {last:.5}");
    // THE MARQUEE PROPERTY: zero mesh rebuilds across the whole run;
    // every refinement event is an octree split.
    assert_eq!(report.total_rebuilds, 0, "no mesh rebuild events, ever");
    assert!(
        report.iterations.iter().all(|it| it.rebuilds == 0),
        "no per-iteration rebuilds either"
    );
    assert!(report.total_splits > 0, "refinement happened as splits");
    // Volume held.
    for it in &report.iterations {
        assert!(
            (it.volume - 0.75).abs() < 0.06,
            "solid fraction held at iter {}: {}",
            it.iter,
            it.volume
        );
    }
    verdict(
        "tm-001",
        "compliance improves at fixed volume with ZERO mesh rebuilds; adaptation \
         enters only as octree splits, straight from the run log",
    );
}

#[test]
fn tm_002_topology_evolves_without_remeshing() {
    // Start from a SOLID plate (one void region: the exterior ring of
    // clamped-low boundary values doesn't exist here — fully solid,
    // zero interior voids) and let the optimizer act: redistribution
    // plus the volume constraint must change the void-component count
    // at SOME point in the run — new topology, same background grid.
    let start = seeded_design();
    let report = run_marquee(start.clone(), 4, 6, 8, 0).expect("runs");
    let void_counts: Vec<usize> = report.iterations.iter().map(|it| it.voids).collect();
    // The evolution witness: lattice nodes whose SIDE of the design
    // boundary flipped between start and finish (the void-count is a
    // coarser signal that only moves on merge/split events).
    let flips = start
        .rho
        .iter()
        .zip(&report.design.rho)
        .filter(|(a, b)| (**a > 0.5) != (**b > 0.5))
        .count();
    #[allow(clippy::cast_precision_loss)]
    let flip_frac = flips as f64 / start.rho.len() as f64;
    println!(
        "{{\"metric\":\"topology\",\"void_counts\":{void_counts:?},\
         \"boundary_flips\":{flips},\"flip_frac\":{flip_frac:.3}}}"
    );
    assert!(
        flips >= 3 || void_counts.windows(2).any(|w| w[0] != w[1]),
        "the design genuinely evolves: {flips} flips, voids {void_counts:?}"
    );
    assert_eq!(report.total_rebuilds, 0, "and still zero rebuilds");
    verdict(
        "tm-002",
        "the boundary-side witness or void-component count changes mid-run — topology \
         genuinely evolves on the same never-rebuilt background grid",
    );
}

#[test]
fn tm_003_dwr_splits_concentrate_at_the_boundary() {
    // The DWR x |grad rho| criterion must put its splits where the
    // design boundary lives, not uniformly.
    let design = seeded_design();
    let report = run_marquee(design, 4, 6, 3, 6).expect("runs");
    // Re-derive the split footprint from the final grid: leaves finer
    // than the base level must mostly sit in the design-boundary halo.
    assert!(report.total_splits >= 6, "enough splits to measure");
    let refined = report.refined_boundary_leaves + report.refined_off_boundary_leaves;
    assert!(refined > 0, "refined leaf footprint is recorded");
    #[allow(clippy::cast_precision_loss)]
    let boundary_frac = report.refined_boundary_leaves as f64 / refined as f64;
    println!(
        "{{\"metric\":\"dwr-boundary-concentration\",\"near\":{},\"far\":{},\
         \"near_frac\":{boundary_frac:.3},\"total_splits\":{}}}",
        report.refined_boundary_leaves, report.refined_off_boundary_leaves, report.total_splits
    );
    // The ghost-penalty contract (CutBandNotUniform) forces the
    // refined band to include the one-cell HALO around every cut cell
    // — the equal-level face-neighbor requirement — so the honest
    // concentration ceiling is band/(band+halo) ~ 2/3, not the 0.80 a
    // halo-free marker could reach.
    assert!(
        boundary_frac >= 0.60,
        "refined leaves concentrate near the design boundary: near {} far {}",
        report.refined_boundary_leaves,
        report.refined_off_boundary_leaves
    );
    verdict(
        "tm-003",
        "DWR-weighted refinement executes its split budget and the final refined \
         leaf footprint is concentrated in the design-boundary halo",
    );
}

#[test]
fn tm_004_benchmark_envelope_and_cadence() {
    // THE GOLDEN ENVELOPE, recorded as evidence (heat-conduction
    // topopt class — the volume-to-point fixture; elasticity envelopes
    // ride CutFEM elasticity when it lands, per the CONTRACT):
    // at 75% solid, 6 iterations, the final compliance must sit inside
    // the recorded band. The band was measured from this fixture and
    // is intentionally loose enough to survive floating-point drift
    // but tight enough to catch regressions.
    let report = run_marquee(seeded_design(), 4, 5, 4, 2).expect("runs");
    let last = report.iterations.last().expect("last").compliance;
    let envelope = (0.005, 0.035);
    println!(
        "{{\"metric\":\"envelope\",\"final\":{last:.5},\"band\":[{},{}],\
         \"timings_ms\":{:?}}}",
        envelope.0,
        envelope.1,
        report
            .iterations
            .iter()
            .map(|it| (it.wall_ms * 10.0).round() / 10.0)
            .collect::<Vec<_>>()
    );
    assert!(
        last > envelope.0 && last < envelope.1,
        "the converged compliance sits in the recorded envelope: {last:.5}"
    );
    // Interactive cadence: iterations don't blow up over the run
    // (splits add cost, but no remeshing spikes — the max/min ratio
    // stays bounded; debug-build measurement, labeled).
    let times: Vec<f64> = report.iterations.iter().map(|it| it.wall_ms).collect();
    let (lo, hi) = times
        .iter()
        .fold((f64::INFINITY, 0.0f64), |(l, h), &t| (l.min(t), h.max(t)));
    assert!(
        hi / lo.max(0.1) < 25.0,
        "per-iteration cost stays bounded (no remeshing cliffs): {times:?}"
    );
    verdict(
        "tm-004",
        "the converged compliance sits inside the ledgered golden envelope; \
         per-iteration wall times (debug label) show no remeshing cliffs",
    );
}

#[test]
fn tm_005_thickness_oracle_audits_the_result() {
    let report = run_marquee(seeded_design(), 4, 6, 6, 0).expect("runs");
    let min_feature = report.design.min_feature_cells();
    println!(
        "{{\"metric\":\"thickness-oracle\",\"min_feature_cells\":{min_feature},\
         \"floor\":2}}"
    );
    // The optimized geometry keeps a resolvable minimum length scale:
    // the medial-axis-class audit (geometry layer auditing optimizer
    // output) must report at least the 2-cell floor.
    assert!(
        min_feature >= 2,
        "the thickness oracle certifies the length-scale floor: {min_feature}"
    );
    verdict(
        "tm-005",
        "the medial-axis thickness oracle reports the optimized design's minimum \
         feature at or above the 2-cell floor — geometry auditing optimization output",
    );
}
