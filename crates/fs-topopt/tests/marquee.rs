//! CutFEM-quadtree marquee conformance (the 2D analogue of the planned
//! octree lane; the b7d0 bead; runs under
//! `cutfem-marquee`). Acceptance: the objective improves at fixed
//! volume with ZERO mesh rebuilds (the marquee property — asserted
//! from the run log); topology EVOLVES (void components change or
//! boundary-side nodes flip) without any remesh; the DWR cut-band gate
//! concentrates refinement near the design boundary; the benchmark envelope is
//! recorded as evidence with per-iteration timing (the
//! interactive-cadence measurement, debug-build labeled); the
//! medial-axis thickness oracle audits the optimized geometry's
//! minimum length scale.
#![cfg(feature = "cutfem-marquee")]

use fs_cutfem::{CellKey, Circle, CutElasticity, CutFemError, CutSdf, Quadtree};
use fs_dwr::estimate_elasticity_compliance;
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::marquee::{
    DWR_CUT_BAND_MASS_GATE, DWR_CUT_BAND_POLICY_VERSION, DensityDesign, DwrBandDecision,
    refine_dwr_cut_band, run_marquee,
};
use std::collections::{BTreeMap, BTreeSet};

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

struct InfiniteHalo<'a>(&'a DensityDesign);

impl CutSdf for InfiniteHalo<'_> {
    fn value(&self, point: [f64; 2]) -> f64 {
        self.0.value(point)
    }

    fn gradient(&self, point: [f64; 2]) -> [f64; 2] {
        self.0.gradient(point)
    }

    fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> Interval {
        if hi[0] - lo[0] > 0.25 || hi[1] - lo[1] > 0.25 {
            if lo[0] > 0.0 {
                Interval::WHOLE
            } else {
                Interval::new(-1.0, 1.0)
            }
        } else {
            self.0.enclose(lo, hi)
        }
    }
}

fn assert_grid_and_level_unchanged(
    grid: &Quadtree,
    leaves_before: &[CellKey],
    level: u32,
    expected_level: u32,
) {
    assert_eq!(grid.leaves().collect::<Vec<_>>(), leaves_before);
    assert_eq!(level, expected_level);
}

fn vector_dwr_problem<'a>(
    grid: &'a Quadtree,
    disk: &'a Circle,
    material: &'a IsotropicElastic,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf: disk,
        material,
        nitsche_beta: 100.0,
        ghost_gamma: 0.5,
        stabilization_scaling: fs_cutfem::CutStabilizationScaling::MuScaled,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-12,
        solver_max_iters: 60_000,
    }
}

#[test]
#[should_panic(expected = "density lattice needs at least 2 nodes per side")]
fn tm_000_rejects_degenerate_lattice() {
    let _ = DensityDesign::uniform(1, 0.5);
}

#[test]
#[should_panic(expected = "density lattice length must equal n*n")]
fn tm_000_rejects_manual_shape_mismatch() {
    let design = DensityDesign {
        n: 3,
        rho: vec![0.5; 8],
    };
    let _ = design.volume();
}

#[test]
fn tm_001_improves_with_zero_rebuilds() {
    // The marquee property is zero mesh rebuilds: the background
    // quadtree is built once, then the configured adaptation path may
    // refine it by quadtree splits. Split concentration is tm-003's
    // executable evidence.
    let report = run_marquee(seeded_design(), 4, 6, 6, true).expect("marquee runs");
    let first = report
        .iterations
        .first()
        .expect("first")
        .analysis_compliance;
    let final_compliance = report.final_compliance;
    println!(
        "{{\"metric\":\"marquee-objective\",\"first_analysis\":{first:.5},\
         \"final_returned_design\":{final_compliance:.5},\
         \"total_splits\":{},\"total_rebuilds\":{}}}",
        report.total_splits, report.total_rebuilds
    );
    assert!(
        final_compliance < first,
        "returned-design compliance improves: {first:.5} -> {final_compliance:.5}"
    );
    // THE MARQUEE PROPERTY: zero mesh rebuilds across the whole run;
    // every refinement event is a quadtree split.
    assert_eq!(report.total_rebuilds, 0, "no mesh rebuild events, ever");
    assert!(
        report.iterations.iter().all(|it| it.rebuilds == 0),
        "no per-iteration rebuilds either"
    );
    assert!(report.total_splits > 0, "refinement happened as splits");
    // Volume held.
    for it in &report.iterations {
        assert!(
            (it.analysis_volume - 0.75).abs() < 0.06 && (it.target_volume - 0.75).abs() < 0.06,
            "analysis/target solid fractions held at iter {}: {}/{}",
            it.iter,
            it.analysis_volume,
            it.target_volume
        );
    }
    verdict(
        "tm-001",
        "compliance improves at fixed volume with ZERO mesh rebuilds; adaptation \
         enters only as quadtree splits, straight from the run log",
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
    let report = run_marquee(start.clone(), 4, 6, 8, false).expect("runs");
    let void_counts: Vec<usize> = report.iterations.iter().map(|it| it.target_voids).collect();
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
#[allow(
    clippy::too_many_lines,
    reason = "the acceptance keeps its complete versioned gate, split-accounting, and final-footprint evidence in one trace"
)]
fn tm_003_dwr_splits_concentrate_at_the_boundary() {
    // The DWR cut-band mass gate must put its refinement where the
    // design boundary lives, not uniformly.
    let base_level = 4u32;
    let max_level = 6u32;
    let iteration_count = 3usize;
    let design = seeded_design();
    let report = run_marquee(design, base_level, max_level, iteration_count, true).expect("runs");
    assert_eq!(report.iterations.len(), iteration_count);
    let initial_grid_leaf_count = Quadtree::with_room(base_level, max_level).leaf_count();
    let mut expected_level = base_level;
    let mut policy_advances = 0usize;
    let mut policy_splits = 0usize;
    let mut summed_splits = 0usize;
    for iteration in &report.iterations {
        let dwr = &iteration.dwr_refinement;
        assert_eq!(dwr.policy_version, DWR_CUT_BAND_POLICY_VERSION);
        assert_eq!(
            dwr.previous_level, expected_level,
            "iteration {} consumes the prior iteration's band level",
            iteration.iter
        );
        assert!(
            dwr.cut_mass >= 0.0 && dwr.cut_mass <= dwr.total_mass,
            "cut mass is a nonnegative subset of total mass"
        );
        assert_eq!(
            dwr.leaves_after - dwr.leaves_before,
            3 * dwr.splits,
            "every reported quadtree split replaces one leaf by four"
        );
        assert_eq!(
            iteration.splits,
            dwr.splits + iteration.post_update_conformance_splits,
            "the forensic split categories exactly reconstruct the iteration total"
        );
        assert_eq!(dwr.analysis_max_level, max_level);
        assert_eq!(
            iteration.target_grid_leaf_count,
            dwr.leaves_after + 3 * iteration.post_update_conformance_splits,
            "motion-conformance accounting reaches the row's target grid"
        );
        let expected_decision = if dwr.total_mass == 0.0 {
            DwrBandDecision::ZeroMass
        } else if dwr.cut_mass <= DWR_CUT_BAND_MASS_GATE * dwr.total_mass {
            DwrBandDecision::GateNotMet
        } else if dwr.previous_level >= max_level {
            DwrBandDecision::LevelHeadroomExhausted
        } else {
            DwrBandDecision::Advanced
        };
        assert_eq!(
            dwr.decision, expected_decision,
            "the exhaustive decision explains the versioned policy outcome"
        );
        if dwr.decision.is_advanced() {
            policy_advances += 1;
            assert_eq!(dwr.band_level, dwr.previous_level + 1);
        } else {
            assert_eq!(dwr.band_level, dwr.previous_level);
            assert_eq!(dwr.splits, 0);
        }
        policy_splits += dwr.splits;
        expected_level = dwr.band_level;
        summed_splits += iteration.splits;
        let cut_fraction = if dwr.total_mass > 0.0 {
            dwr.cut_mass / dwr.total_mass
        } else {
            0.0
        };
        println!(
            "{{\"metric\":\"dwr-band-time-level\",\"iteration\":{},\
             \"analysis_design_witness\":\"{:016x}\",\
             \"target_design_witness\":\"{:016x}\",\
             \"analysis_compliance\":{:.10e},\"analysis_volume\":{:.8},\
             \"target_volume\":{:.8},\"analysis_voids\":{},\"target_voids\":{},\
             \"policy_version\":{},\
             \"cut_mass\":{:.10e},\"total_mass\":{:.10e},\
             \"cut_fraction\":{cut_fraction:.8},\"decision\":\"{:?}\",\
             \"advanced\":{},\
             \"level_before\":{},\"level_after\":{},\"leaves_before\":{},\
             \"policy_leaves_after\":{},\"target_grid_leaves\":{},\
             \"policy_splits\":{},\
             \"post_update_conformance_splits\":{},\"iteration_splits\":{}}}",
            iteration.iter,
            iteration.analysis_design_witness,
            iteration.target_design_witness,
            iteration.analysis_compliance,
            iteration.analysis_volume,
            iteration.target_volume,
            iteration.analysis_voids,
            iteration.target_voids,
            dwr.policy_version,
            dwr.cut_mass,
            dwr.total_mass,
            dwr.decision,
            dwr.decision.is_advanced(),
            dwr.previous_level,
            dwr.band_level,
            dwr.leaves_before,
            dwr.leaves_after,
            iteration.target_grid_leaf_count,
            dwr.splits,
            iteration.post_update_conformance_splits,
            iteration.splits,
        );
    }
    for pair in report.iterations.windows(2) {
        assert_eq!(
            pair[0].target_design_witness, pair[1].analysis_design_witness,
            "the forensic witnesses correlate a row's target with the next analysis design"
        );
        assert_eq!(
            pair[0].target_grid_leaf_count, pair[1].dwr_refinement.leaves_before,
            "split-only accounting permits no hidden leaf-count change before the next analysis"
        );
    }
    let last_iteration = report.iterations.last().expect("nonempty fixture run");
    assert_eq!(
        last_iteration.target_design_witness, report.final_design_witness,
        "the forensic witnesses correlate the last target with the returned design"
    );
    assert_eq!(
        last_iteration.target_grid_leaf_count, report.final_grid_leaf_count,
        "split-only accounting permits no hidden leaf-count change before the final re-solve"
    );
    assert_eq!(
        report.final_design_witness,
        report.design.state_witness(),
        "the report witness recomputes from the returned design"
    );
    assert!(report.final_compliance.is_finite());
    assert!(report.final_grid_leaf_count > 0);
    assert_eq!(
        report.total_splits, summed_splits,
        "per-iteration split evidence reconstructs the run total"
    );
    assert_eq!(
        report.final_grid_leaf_count,
        initial_grid_leaf_count + 3 * report.total_splits,
        "the final tree independently reconstructs every 1-to-4 split wave"
    );
    let level_headroom = usize::try_from(max_level - base_level).expect("level delta fits usize");
    let maximum_advances = iteration_count.min(level_headroom);
    assert!(policy_advances > 0, "the DWR gate authorizes refinement");
    assert!(
        policy_splits > 0,
        "the fixture's DWR authorizations perform structural splits"
    );
    assert!(
        policy_advances <= maximum_advances,
        "one-level policy advances stay within iteration and level headroom"
    );
    // Re-derive the split footprint from the final grid: leaves finer
    // than the base level must mostly sit in the design-boundary halo.
    assert!(
        report.total_splits > 0,
        "an authorized band advance creates a measurable footprint"
    );
    let refined = report.refined_boundary_leaves + report.refined_off_boundary_leaves;
    assert!(refined > 0, "refined leaf footprint is recorded");
    #[allow(clippy::cast_precision_loss)]
    let boundary_frac = report.refined_boundary_leaves as f64 / refined as f64;
    println!(
        "{{\"metric\":\"dwr-boundary-concentration\",\"near\":{},\"far\":{},\
         \"near_frac\":{boundary_frac:.3},\"total_splits\":{}}}",
        report.refined_boundary_leaves, report.refined_off_boundary_leaves, report.total_splits
    );
    // Fixture-golden final-footprint diagnostic: the ghost-penalty contract
    // (CutBandNotUniform) forces the
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
        "the DWR cut-band mass gate advances refinement and the final refined leaf \
         footprint is concentrated in the design-boundary halo",
    );
}

#[test]
fn tm_004_benchmark_envelope_and_cadence() {
    // THE GOLDEN ENVELOPE, recorded as evidence (heat-conduction
    // topopt class — the volume-to-point fixture; elasticity envelopes
    // ride CutFEM elasticity when it lands, per the CONTRACT):
    // at 75% solid, 4 updates, the returned design's re-solved compliance
    // must sit inside
    // the recorded band. The band was measured from this fixture and
    // is intentionally loose enough to survive floating-point drift
    // but tight enough to catch regressions.
    let report = run_marquee(seeded_design(), 4, 5, 4, true).expect("runs");
    let final_compliance = report.final_compliance;
    let envelope = (0.005, 0.035);
    println!(
        "{{\"metric\":\"envelope\",\"final\":{final_compliance:.5},\"band\":[{},{}],\
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
        final_compliance > envelope.0 && final_compliance < envelope.1,
        "the returned-design compliance sits in the recorded envelope: {final_compliance:.5}"
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
        "the returned design's re-solved compliance sits inside the ledgered golden envelope; \
         per-iteration wall times (debug label) show no remeshing cliffs",
    );
}

#[test]
fn tm_005_thickness_oracle_audits_the_result() {
    let report = run_marquee(seeded_design(), 4, 6, 6, false).expect("runs");
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

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "the shared-policy acceptance keeps deterministic replay, metamorphic variants, and disabled behavior together"
)]
fn tm_006_synthetic_indicators_drive_shared_band_policy_once() {
    // This is deliberately estimator-agnostic helper coverage. Vector
    // compliance DWR can supply the same CellKey map, but this test does not
    // claim a graded vector re-solve: it exercises one planning decision only.
    let design = seeded_design();
    let make_grid = || Quadtree::with_room(2, 4);
    let seed_grid = make_grid();
    let indicators: BTreeMap<_, _> = seed_grid
        .leaves()
        .map(|cell| {
            let (lo, hi) = seed_grid.rect(cell);
            let eta = if design.enclose(lo, hi).contains_zero() {
                1.0
            } else {
                1.0 / 64.0
            };
            (cell, eta)
        })
        .collect();

    let mut grid_a = make_grid();
    let mut grid_b = make_grid();
    let mut level_a = 2;
    let mut level_b = 2;
    let decision_a = refine_dwr_cut_band(&mut grid_a, &design, &indicators, &mut level_a, true)
        .expect("valid shared policy input");
    let decision_b = refine_dwr_cut_band(&mut grid_b, &design, &indicators, &mut level_b, true)
        .expect("valid shared policy input");
    let leaves_a: Vec<_> = grid_a.leaves().collect();
    let leaves_b: Vec<_> = grid_b.leaves().collect();

    assert_eq!(decision_a, decision_b, "policy evidence is deterministic");
    assert_eq!(
        leaves_a, leaves_b,
        "the structural split set is deterministic"
    );
    assert_eq!(
        decision_a.decision,
        DwrBandDecision::Advanced,
        "cut-band indicator mass passes the gate"
    );
    assert!(
        decision_a.splits > 0,
        "one planning step performs actual splits"
    );
    assert_eq!(decision_a.previous_level, 2);
    assert_eq!(decision_a.band_level, 3);
    assert_eq!(decision_a.policy_version, DWR_CUT_BAND_POLICY_VERSION);
    assert_eq!(
        decision_a.leaves_after - decision_a.leaves_before,
        3 * decision_a.splits
    );

    // G3: this policy uses absolute aggregate partition mass, not signs,
    // insertion history, or within-partition magnitude ranking. Negation,
    // exact power-of-two scaling, reverse insertion, and mass-preserving
    // redistribution preserve the authorization and target topology.
    let negated: BTreeMap<_, _> = indicators
        .iter()
        .map(|(&cell, &eta)| (cell, -eta))
        .collect();
    let scaled: BTreeMap<_, _> = indicators
        .iter()
        .map(|(&cell, &eta)| (cell, 2.0 * eta))
        .collect();
    let reverse_inserted: BTreeMap<_, _> = indicators
        .iter()
        .rev()
        .map(|(&cell, &eta)| (cell, eta))
        .collect();
    let mut cut_cells = Vec::new();
    let mut off_cut_cells = Vec::new();
    for &cell in indicators.keys() {
        let (lo, hi) = seed_grid.rect(cell);
        if design.enclose(lo, hi).contains_zero() {
            cut_cells.push(cell);
        } else {
            off_cut_cells.push(cell);
        }
    }
    assert!(cut_cells.len() >= 2 && off_cut_cells.len() >= 2);
    let mut redistributed = indicators.clone();
    assert_eq!(redistributed.insert(cut_cells[0], 1.5), Some(1.0));
    assert_eq!(redistributed.insert(cut_cells[1], 0.5), Some(1.0));
    assert_eq!(
        redistributed.insert(off_cut_cells[0], 3.0 / 128.0),
        Some(1.0 / 64.0)
    );
    assert_eq!(
        redistributed.insert(off_cut_cells[1], 1.0 / 128.0),
        Some(1.0 / 64.0)
    );
    for (name, variant, mass_scale) in [
        ("negated", negated, 1.0),
        ("power-of-two-scaled", scaled, 2.0),
        ("reverse-inserted", reverse_inserted, 1.0),
        ("within-partition-redistributed", redistributed, 1.0),
    ] {
        let mut variant_grid = make_grid();
        let mut variant_level = 2;
        let variant_decision = refine_dwr_cut_band(
            &mut variant_grid,
            &design,
            &variant,
            &mut variant_level,
            true,
        )
        .expect("metamorphic policy input remains valid");
        assert_eq!(variant_decision.decision, decision_a.decision, "{name}");
        assert_eq!(
            variant_decision.previous_level, decision_a.previous_level,
            "{name}"
        );
        assert_eq!(variant_decision.band_level, decision_a.band_level, "{name}");
        assert_eq!(variant_decision.splits, decision_a.splits, "{name}");
        assert_eq!(
            variant_decision.cut_mass.to_bits(),
            (mass_scale * decision_a.cut_mass).to_bits(),
            "{name} cut mass"
        );
        assert_eq!(
            variant_decision.total_mass.to_bits(),
            (mass_scale * decision_a.total_mass).to_bits(),
            "{name} total mass"
        );
        assert_eq!(
            variant_grid.leaves().collect::<Vec<_>>(),
            leaves_a,
            "{name}"
        );
    }

    let mut disabled_grid = make_grid();
    let disabled_leaf_count = disabled_grid.leaf_count();
    let mut disabled_level = 2;
    let disabled = refine_dwr_cut_band(
        &mut disabled_grid,
        &design,
        &indicators,
        &mut disabled_level,
        false,
    )
    .expect("valid disabled policy input");
    assert_eq!(
        disabled.decision,
        DwrBandDecision::Disabled,
        "disabled policy is a no-op"
    );
    assert_eq!(disabled.splits, 0);
    assert_eq!(disabled_grid.leaf_count(), disabled_leaf_count);
    assert_eq!(disabled_level, 2);

    let mut zero_grid = make_grid();
    let zero_before: Vec<_> = zero_grid.leaves().collect();
    let zero_indicators: BTreeMap<_, _> = zero_before
        .iter()
        .copied()
        .map(|cell| (cell, 0.0))
        .collect();
    let mut zero_level = 2;
    let zero = refine_dwr_cut_band(
        &mut zero_grid,
        &design,
        &zero_indicators,
        &mut zero_level,
        true,
    )
    .expect("zero mass is a valid deterministic no-op");
    assert_eq!(zero.decision, DwrBandDecision::ZeroMass);
    assert_eq!(zero.splits, 0);
    assert_grid_and_level_unchanged(&zero_grid, &zero_before, zero_level, 2);

    let mut threshold_grid = make_grid();
    let threshold_before: Vec<_> = threshold_grid.leaves().collect();
    let threshold_indicators = BTreeMap::from([
        (cut_cells[0], DWR_CUT_BAND_MASS_GATE),
        (off_cut_cells[0], 1.0 - DWR_CUT_BAND_MASS_GATE),
    ]);
    let mut threshold_level = 2;
    let threshold = refine_dwr_cut_band(
        &mut threshold_grid,
        &design,
        &threshold_indicators,
        &mut threshold_level,
        true,
    )
    .expect("exact threshold is a valid deterministic no-op");
    assert_eq!(
        threshold.cut_mass.to_bits(),
        (DWR_CUT_BAND_MASS_GATE * threshold.total_mass).to_bits(),
        "the fixture lands bit-exactly on the policy threshold"
    );
    assert_eq!(threshold.decision, DwrBandDecision::GateNotMet);
    assert_eq!(threshold.splits, 0);
    assert_grid_and_level_unchanged(&threshold_grid, &threshold_before, threshold_level, 2);

    let mut exhausted_grid = Quadtree::with_room(2, 2);
    let exhausted_before: Vec<_> = exhausted_grid.leaves().collect();
    let mut exhausted_level = 2;
    let exhausted = refine_dwr_cut_band(
        &mut exhausted_grid,
        &design,
        &indicators,
        &mut exhausted_level,
        true,
    )
    .expect("exhausted headroom is a valid deterministic no-op");
    assert_eq!(exhausted.decision, DwrBandDecision::LevelHeadroomExhausted);
    assert_eq!(exhausted.analysis_max_level, 2);
    assert_eq!(exhausted.splits, 0);
    assert_grid_and_level_unchanged(&exhausted_grid, &exhausted_before, exhausted_level, 2);

    println!(
        "{{\"metric\":\"dwr-shared-band-policy\",\"cut_mass\":{:.3},\
         \"total_mass\":{:.3},\"decision\":\"{:?}\",\"splits\":{},\"band_level\":{}}}",
        decision_a.cut_mass,
        decision_a.total_mass,
        decision_a.decision,
        decision_a.splits,
        decision_a.band_level
    );
    verdict(
        "tm-006",
        "synthetic CellKey indicators deterministically drive one shared cut-band \
         planning step; no graded vector re-solve is claimed",
    );
}

#[test]
fn tm_007_invalid_band_policy_inputs_refuse_without_mutation() {
    let design = seeded_design();
    let make_grid = || Quadtree::with_room(2, 4);

    let mut nonfinite_grid = make_grid();
    let nonfinite_before: Vec<_> = nonfinite_grid.leaves().collect();
    let nonfinite_cell = nonfinite_before[0];
    let mut nonfinite_level = 2;
    let nonfinite = BTreeMap::from([(nonfinite_cell, f64::NAN)]);
    let err = refine_dwr_cut_band(
        &mut nonfinite_grid,
        &design,
        &nonfinite,
        &mut nonfinite_level,
        true,
    )
    .expect_err("non-finite indicator must refuse");
    assert!(matches!(err, CutFemError::InvalidFemInput { .. }));
    assert_grid_and_level_unchanged(&nonfinite_grid, &nonfinite_before, nonfinite_level, 2);

    let mut overflow_grid = make_grid();
    let overflow_before: Vec<_> = overflow_grid.leaves().collect();
    let mut overflow_level = 2;
    let overflow = BTreeMap::from([
        (overflow_before[0], f64::MAX),
        (overflow_before[1], f64::MAX),
    ]);
    let err = refine_dwr_cut_band(
        &mut overflow_grid,
        &design,
        &overflow,
        &mut overflow_level,
        true,
    )
    .expect_err("non-finite accumulated mass must refuse");
    assert!(matches!(err, CutFemError::InvalidFemInput { .. }));
    assert_grid_and_level_unchanged(&overflow_grid, &overflow_before, overflow_level, 2);

    let mut nonleaf_grid = make_grid();
    let nonleaf_before: Vec<_> = nonleaf_grid.leaves().collect();
    let mut nonleaf_level = 2;
    let nonleaf = BTreeMap::from([((1, 0, 0), 1.0)]);
    let err = refine_dwr_cut_band(
        &mut nonleaf_grid,
        &design,
        &nonleaf,
        &mut nonleaf_level,
        true,
    )
    .expect_err("non-leaf indicator key must refuse");
    assert!(matches!(err, CutFemError::InvalidFemInput { .. }));
    assert_grid_and_level_unchanged(&nonleaf_grid, &nonleaf_before, nonleaf_level, 2);

    let mut enclosure_grid = make_grid();
    let enclosure_before: Vec<_> = enclosure_grid.leaves().collect();
    let mut enclosure_level = 2;
    let indicators: BTreeMap<_, _> = enclosure_before
        .iter()
        .copied()
        .map(|cell| {
            let (lo, hi) = enclosure_grid.rect(cell);
            let eta = if design.enclose(lo, hi).contains_zero() {
                1.0
            } else {
                0.01
            };
            (cell, eta)
        })
        .collect();
    let err = refine_dwr_cut_band(
        &mut enclosure_grid,
        &InfiniteHalo(&design),
        &indicators,
        &mut enclosure_level,
        true,
    )
    .expect_err("non-finite recursive halo enclosure must refuse");
    assert!(matches!(err, CutFemError::InvalidFemInput { .. }));
    assert_grid_and_level_unchanged(&enclosure_grid, &enclosure_before, enclosure_level, 2);

    let mut level_grid = make_grid();
    let level_before: Vec<_> = level_grid.leaves().collect();
    let level_cell = level_before[0];
    let mut invalid_level = 5;
    let indicators = BTreeMap::from([(level_cell, 1.0)]);
    let err = refine_dwr_cut_band(
        &mut level_grid,
        &design,
        &indicators,
        &mut invalid_level,
        true,
    )
    .expect_err("out-of-range band level must refuse");
    assert!(matches!(err, CutFemError::InvalidFemInput { .. }));
    assert_grid_and_level_unchanged(&level_grid, &level_before, invalid_level, 5);

    verdict(
        "tm-007",
        "non-finite indicators, accumulated masses, and recursive enclosures, \
         non-leaf keys, and invalid levels return structured refusals without \
         mutating grid or band level",
    );
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "the acceptance retains two complete estimate/refine/re-solve cycles and their replay evidence"
)]
fn tm_008_real_vector_dwr_advances_and_resolves_twice_on_graded_trees() {
    // This is the fs-dwr effectivity battery's manufactured all-embedded disk
    // family. The radius is large enough to retain coarse interior active
    // leaves while the DWR policy advances the cut band and halo, so every
    // post-refinement solve genuinely exercises a mixed-level active space.
    // The real estimator map is passed through unchanged on both cycles.
    let mut grid = Quadtree::with_room(3, 5);
    let disk = Circle {
        center: [0.47, 0.53],
        radius: 0.33,
    };
    let material = IsotropicElastic::new(1.0, 0.3, 10.0).expect("compressible material");
    let (lambda, mu) = material.lame();
    let displacement_scale = 0.01;
    let body = |_: f64, _: f64| [-2.0 * displacement_scale * (lambda + 3.0 * mu), 0.0];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let mut band_level = 3;
    let mut total_splits = 0usize;

    for cycle in 0..2u32 {
        let estimate = {
            let problem = vector_dwr_problem(&grid, &disk, &material);
            estimate_elasticity_compliance(&problem, &body, &zero)
                .expect("real vector DWR estimate on current grid")
        };
        assert!(
            estimate.eta_abs.is_finite() && estimate.eta_abs > 0.0,
            "cycle {cycle} must retain nondegenerate authentic marking mass"
        );
        let before_leaf_count = grid.leaf_count();
        let mut replay_grid = grid.clone();
        let mut replay_level = band_level;
        let decision = refine_dwr_cut_band(
            &mut grid,
            &disk,
            &estimate.indicators,
            &mut band_level,
            true,
        )
        .expect("real vector indicators satisfy shared policy input contract");
        let replay_decision = refine_dwr_cut_band(
            &mut replay_grid,
            &disk,
            &estimate.indicators,
            &mut replay_level,
            true,
        )
        .expect("real vector refinement replay");
        assert_eq!(decision, replay_decision, "cycle {cycle} policy replay");
        assert_eq!(band_level, replay_level, "cycle {cycle} band replay");
        assert_eq!(
            grid.leaves().collect::<Vec<_>>(),
            replay_grid.leaves().collect::<Vec<_>>(),
            "cycle {cycle} structural replay"
        );
        assert_eq!(
            decision.total_mass.to_bits(),
            estimate.eta_abs.to_bits(),
            "the helper consumes the estimator's exact marking mass"
        );
        let cut_mass_fraction = decision.cut_mass / decision.total_mass;
        assert!(
            cut_mass_fraction.is_finite() && cut_mass_fraction > DWR_CUT_BAND_MASS_GATE,
            "cycle {cycle} authentic boundary-dominated indicators must pass the shared gate: {cut_mass_fraction:.6}"
        );
        assert_eq!(decision.policy_version, DWR_CUT_BAND_POLICY_VERSION);
        assert_eq!(
            decision.decision,
            DwrBandDecision::Advanced,
            "cycle {cycle} must advance the band"
        );
        assert!(
            decision.splits > 0,
            "cycle {cycle} must perform real splits"
        );
        assert_eq!(decision.previous_level, 3 + cycle);
        assert_eq!(decision.band_level, 4 + cycle);
        assert_eq!(band_level, 4 + cycle);
        assert_eq!(
            grid.leaf_count() - before_leaf_count,
            3 * decision.splits,
            "each quadtree split replaces one leaf by four"
        );
        assert_eq!(decision.leaves_before, before_leaf_count);
        assert_eq!(decision.leaves_after, grid.leaf_count());
        total_splits += decision.splits;

        let adapted = vector_dwr_problem(&grid, &disk, &material)
            .solve(&body, &zero)
            .expect("graded vector re-solve after authentic refinement");
        let replay = vector_dwr_problem(&replay_grid, &disk, &material)
            .solve(&body, &zero)
            .expect("graded vector re-solve replay");
        let active_levels: BTreeSet<_> = adapted.active_cells().iter().map(|cell| cell.0).collect();
        assert!(
            active_levels.len() >= 2,
            "cycle {cycle} must re-solve on genuinely mixed active leaves: {active_levels:?}"
        );
        assert!(
            adapted.compliance().is_finite()
                && adapted.rel_residual.is_finite()
                && adapted.dof_count() > 0,
            "cycle {cycle} graded solve evidence must be finite and nonempty"
        );
        assert_eq!(adapted.active_cells(), replay.active_cells());
        assert_eq!(adapted.coefficients().len(), replay.coefficients().len());
        for (left, right) in adapted.coefficients().iter().zip(replay.coefficients()) {
            assert_eq!(
                left.to_bits(),
                right.to_bits(),
                "cycle {cycle} coefficient replay"
            );
        }
        assert_eq!(adapted.nodal().len(), replay.nodal().len());
        for ((left_node, left), (right_node, right)) in adapted.nodal().iter().zip(replay.nodal()) {
            assert_eq!(left_node, right_node);
            assert_eq!(left[0].to_bits(), right[0].to_bits());
            assert_eq!(left[1].to_bits(), right[1].to_bits());
        }
        assert_eq!(
            adapted.compliance().to_bits(),
            replay.compliance().to_bits(),
            "cycle {cycle} compliance replay"
        );

        let enriched_delta = estimate.j_enriched - estimate.j_primal;
        let eta_over_enriched_delta = estimate.eta_signed / enriched_delta;
        assert!(
            enriched_delta.is_finite()
                && enriched_delta.abs() > 0.0
                && eta_over_enriched_delta.is_finite(),
            "cycle {cycle} enriched-compliance diagnostic must be finite and non-degenerate"
        );
        let active_min_level = *active_levels.first().expect("nonempty active levels");
        let active_max_level = *active_levels.last().expect("nonempty active levels");
        println!(
            "{{\"metric\":\"real-vector-dwr-graded-resolve-cycle\",\"cycle\":{cycle},\"j_h\":{:.10e},\"j_h2\":{:.10e},\"eta_signed\":{:.10e},\"eta_abs\":{:.10e},\"eta_over_enriched_delta\":{eta_over_enriched_delta:.8},\"cut_mass_fraction\":{cut_mass_fraction:.8},\"splits\":{},\"band_level\":{},\"adapted_dofs\":{},\"adapted_compliance\":{:.10e},\"active_min_level\":{active_min_level},\"active_max_level\":{active_max_level}}}",
            estimate.j_primal,
            estimate.j_enriched,
            estimate.eta_signed,
            estimate.eta_abs,
            decision.splits,
            decision.band_level,
            adapted.dof_count(),
            adapted.compliance(),
        );
    }

    assert!(total_splits > 0);
    assert_eq!(band_level, 5);
    println!(
        "{{\"metric\":\"real-vector-dwr-graded-resolve\",\"cycles\":2,\"total_splits\":{total_splits},\"final_band_level\":{band_level}}}"
    );
    verdict(
        "tm-008",
        "two authentic vector-compliance DWR estimates advance the shared cut-band policy; \
         each positive split wave is followed by a deterministic mixed-level vector re-solve",
    );
}
