//! fs-grammar-e2e — GrammarForge: certified-fabricable geometric program
//! discovery. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! A CAD model is one hand-built artifact with no guarantees. This instead
//! ILLUMINATES the diverse family of CSG PROGRAMS that approximate a target
//! shape and are fabricable, composing crates never designed to meet:
//!
//! - **Programs as data** ([`fs_shapeprog`]): a candidate is a CSG program —
//!   `sphere(r₁)@−d ∪ (sphere(r₂)⊕o)@+d`. Its fidelity is the worst-case SDF
//!   discrepancy from the target over a sample grid.
//! - **Certificate-preserving simplification** ([`fs_shapeprog::simplify`]): the
//!   rewrite engine drops redundant tiny offsets and applies geometric
//!   identities, each with a fidelity certificate (`Exact` or `Approximate{bound}`),
//!   so the simplified program is provably within `max_error` of the original —
//!   and the campaign INDEPENDENTLY re-measures the discrepancy to confirm the
//!   certificate holds (certifying the certifier).
//! - **Manufacturability** ([`fs_fab`]): a minimum-feature-size constraint scores
//!   each program's smallest feature — the fabrication margin.
//! - **Illumination** ([`fs_archive`]): MAP-Elites over (program size × fab
//!   margin) keeps the best-matching program in every complexity/fabricability
//!   niche — the diverse atlas, not one model.
//! - **Honest colors** ([`fs_evidence`]): a program that matches within tolerance,
//!   is fab-satisfied, and simplifies soundly is `Verified`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_archive::MapElites;
use fs_evidence::Color;
use fs_fab::min_feature_size;
use fs_shapeprog::{Geom, max_sdf_discrepancy, simplify};

/// The target shape: a "peanut" — two unit spheres at `x = ±0.8`.
#[must_use]
pub fn target() -> Geom {
    Geom::sphere(1.0)
        .translate([-0.8, 0.0, 0.0])
        .union(Geom::sphere(1.0).translate([0.8, 0.0, 0.0]))
}

/// Build a candidate program from parameters `[r1, r2, d, o]`.
#[must_use]
pub fn build_program(r1: f64, r2: f64, d: f64, o: f64) -> Geom {
    let left = Geom::sphere(r1).translate([-d, 0.0, 0.0]);
    let right = Geom::sphere(r2).offset(o).translate([d, 0.0, 0.0]);
    left.union(right)
}

/// A deterministic 3-D sample grid over `[-2, 2]³` for SDF discrepancy.
#[must_use]
fn sample_points() -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    let n = 7;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let c = |t: usize| -2.0 + 4.0 * t as f64 / (n - 1) as f64;
                pts.push([c(i), c(j), c(k)]);
            }
        }
    }
    pts
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct GrammarReport {
    /// Fraction of (size × fab-margin) niches filled.
    pub coverage: f64,
    /// Quality-diversity score (Σ elite fitness = −Σ discrepancy).
    pub qd_score: f64,
    /// Number of filled niches.
    pub num_elites: usize,
    /// The best (lowest) SDF discrepancy from the target.
    pub best_discrepancy: f64,
    /// The best program's parameters `[r1, r2, d, o]`.
    pub best_params: [f64; 4],
    /// Elites whose program simplified to something smaller.
    pub simplified_count: usize,
    /// Total elite program size before simplification.
    pub size_before: usize,
    /// Total elite program size after simplification.
    pub size_after: usize,
    /// The largest certified simplification error over all elites (`≤ tol`).
    pub max_certified_error: f64,
    /// Did EVERY simplification's re-measured discrepancy stay within its
    /// certified bound? (certifying the rewrite engine).
    pub simplification_sound: bool,
    /// Elites that satisfy the minimum-feature-size fabrication constraint.
    pub fab_satisfied: usize,
    /// The headline color: `Verified` iff the best design matches within
    /// tolerance, is fab-satisfied, and simplifies soundly.
    pub headline_color: Color,
}

/// Run the GrammarForge campaign; `match_tol` is the SDF discrepancy under which
/// a program is deemed to match the target; `simplify_tol` bounds dropped offsets.
#[must_use]
pub fn run_campaign(match_tol: f64, simplify_tol: f64) -> GrammarReport {
    let target = target();
    let samples = sample_points();
    // A minimum-feature-size rule that actually discriminates: the thinnest
    // (r=0.7) spheres fail it, so not every program is fabricable.
    let fab = min_feature_size(0.8);

    let r_vals = [0.7, 0.9, 1.0, 1.1];
    let d_vals = [0.6, 0.8, 1.0];
    let o_vals = [0.0, 0.02, 0.05];

    // MAP-Elites over (total material `r1+r2`, dipole separation `d`) — the
    // behavioral axes of the peanut family.
    let mut archive = MapElites::new(vec![1.3, 0.5], vec![2.3, 1.1], vec![6, 4]);
    for &r1 in &r_vals {
        for &r2 in &r_vals {
            for &d in &d_vals {
                for &o in &o_vals {
                    let prog = build_program(r1, r2, d, o);
                    let disc = max_sdf_discrepancy(&prog, &target, &samples);
                    // Closeness score in (0, 1] — higher is a better match; the
                    // archive requires a non-negative fitness.
                    let fitness = 1.0 / (1.0 + disc);
                    let descriptor = vec![r1 + r2, d];
                    archive.add(vec![r1, r2, d, o], descriptor, fitness);
                }
            }
        }
    }

    // Post-process the elites: simplification soundness + fabrication tally.
    let (mut simplified_count, mut size_before, mut size_after, mut fab_satisfied) = (0, 0, 0, 0);
    let mut max_certified_error = 0.0_f64;
    let mut simplification_sound = true;
    for e in archive.elites() {
        let prog = build_program(e.solution[0], e.solution[1], e.solution[2], e.solution[3]);
        let before = prog.size();
        let simp = simplify(&prog, simplify_tol);
        let after = simp.program.size();
        size_before += before;
        size_after += after;
        if after < before {
            simplified_count += 1;
        }
        max_certified_error = max_certified_error.max(simp.max_error);
        // Independently re-measure: the simplified program must stay within its
        // OWN certified error bound of the original (certifying the certifier).
        let actual = max_sdf_discrepancy(&prog, &simp.program, &samples);
        if actual > simp.max_error + 1e-9 {
            simplification_sound = false;
        }
        if fab.satisfied(e.solution[0].min(e.solution[1])) {
            fab_satisfied += 1;
        }
    }

    let best = archive.best().expect("archive has at least one elite");
    let best_discrepancy = 1.0 / best.fitness - 1.0;
    let best_params = [
        best.solution[0],
        best.solution[1],
        best.solution[2],
        best.solution[3],
    ];
    let best_fab_ok = fab.satisfied(best.solution[0].min(best.solution[1]));
    let headline_color = if best_discrepancy <= match_tol && best_fab_ok && simplification_sound {
        // declared-color-ok: demo headline candidate from local discrepancy/fabricability checks; admitted only at a consumer's authority boundary (6pf9)
        Color::Verified {
            lo: 0.0,
            hi: best_discrepancy,
        }
    } else {
        Color::Estimated {
            estimator: "grammar-open".to_string(),
            dispersion: best_discrepancy,
        }
    };

    GrammarReport {
        coverage: archive.coverage(),
        qd_score: archive.qd_score(),
        num_elites: archive.num_elites(),
        best_discrepancy,
        best_params,
        simplified_count,
        size_before,
        size_after,
        max_certified_error,
        simplification_sound,
        fab_satisfied,
        headline_color,
    }
}
