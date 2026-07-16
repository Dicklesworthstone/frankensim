//! fs-dwr conformance battery (bead tfz.23).
//!
//! - dwr-001 G0+G3: adjoint-weighting consistency — effectivity
//!   indices near 1 against known-truth goals on two frontends (disk
//!   Nitsche, strip strong+Nitsche), and estimator monotonicity under
//!   refinement at the theoretical rate.
//! - dwr-002 P2: marking determinism (bitwise-identical marks across
//!   runs) and the Dörfler prefix property.
//! - dwr-003: goal-oriented beats uniform on the localized-QoI fixture
//!   (small window far from a sharp source) — accuracy-per-DOF,
//!   ledgered.
//! - dwr-004: anisotropic metric synthesis — alignment with the layer,
//!   complexity normalization, and a metric-instantiated graded mesh
//!   beating isotropic at equal DOF.
//! - dwr-005: DWR-weighted Haar thresholding — goal impact under
//!   budget at ≥5× compression, and measurably better than unweighted
//!   thresholding at matched compression.
//! - dwr-006: the h-vs-p decision signal — kinks route to h, smooth
//!   regions to p.

use fs_cutfem::quad::{cut_cell_rules, tensor_gauss};
use fs_cutfem::{Circle, CutSdf, FemParams, HalfPlane, Quadtree};
use fs_dwr::{
    Decision, GoalContext, adapt_loop, dorfler, estimate, h_vs_p, haar_threshold, synthesize_metric,
};
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

// ---------------------------------------------------------------- fixtures

fn mms_u(x: f64, y: f64) -> f64 {
    (PI * x).sin() * (PI * y).sin()
}
fn mms_f(x: f64, y: f64) -> f64 {
    2.0 * PI * PI * mms_u(x, y)
}

fn goal_weight(x: f64, y: f64) -> f64 {
    let r2 = (x - 0.55) * (x - 0.55) + (y - 0.5) * (y - 0.5);
    (-100.0 * r2).exp()
}

/// ∫_Ω jw·u_exact by high-depth certified quadrature (no solve).
fn exact_goal(sdf: &dyn CutSdf, u: &dyn Fn(f64, f64) -> f64, jw: &dyn Fn(f64, f64) -> f64) -> f64 {
    let grid = Quadtree::uniform(6);
    let mut j = 0.0;
    for c in grid.leaves() {
        let (lo, hi) = grid.rect(c);
        let iv = sdf.enclose(lo, hi);
        let rule = if iv.hi() < 0.0 {
            let mut v = Vec::with_capacity(9);
            tensor_gauss(lo, hi, &mut v);
            v
        } else if iv.lo() > 0.0 {
            Vec::new()
        } else {
            cut_cell_rules(sdf, lo, hi, 4).bulk
        };
        for (p, w) in rule {
            j += w * jw(p[0], p[1]) * u(p[0], p[1]);
        }
    }
    j
}

// ------------------------------------------------------------------ dwr-001

#[test]
fn dwr_001_effectivity_and_monotonicity() {
    let mut rows = String::new();
    let mut effectivities = Vec::new();
    let mut rates = Vec::new();
    // Frontend 1: all-embedded disk (Nitsche data from the MMS).
    let disk = Circle {
        center: [0.5, 0.5],
        radius: 0.35,
    };
    // Frontend 2: strip (strong outer + Nitsche embedded side).
    let strip = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.55,
    };
    let fixtures: [(&str, &dyn CutSdf, FemParams); 2] = [
        ("disk", &disk, FemParams::default()),
        (
            "strip",
            &strip,
            FemParams {
                strong_outer: true,
                ..FemParams::default()
            },
        ),
    ];
    let strip_goal = |x: f64, y: f64| {
        let r2 = (x - 0.3) * (x - 0.3) + (y - 0.5) * (y - 0.5);
        (-100.0 * r2).exp()
    };
    for (name, sdf, params) in fixtures {
        // The strip's goal sits well inside Ω (a bump ON Γ would put
        // half its mass outside the domain).
        let jw: &dyn Fn(f64, f64) -> f64 = if name == "strip" {
            &strip_goal
        } else {
            &goal_weight
        };
        let goal = GoalContext { weight: jw };
        let j_true = exact_goal(sdf, &mms_u, jw);
        let mut etas = Vec::new();
        for level in [4u32, 5] {
            let grid = Quadtree::uniform(level);
            let est = estimate(&grid, sdf, params, &mms_f, &mms_u, &goal).expect("estimates");
            let true_err = j_true - est.j_primal;
            let eff = est.eta_signed / true_err;
            let _ = write!(
                rows,
                "{{\"fixture\":\"{name}\",\"level\":{level},\"j_h\":{:.6e},\
                 \"true_err\":{true_err:.3e},\"eta_signed\":{:.3e},\
                 \"effectivity\":{eff:.3}}},",
                est.j_primal, est.eta_signed
            );
            effectivities.push(eff);
            etas.push(est.eta_abs);
        }
        rates.push((etas[0] / etas[1]).log2());
    }
    let eff_ok = effectivities.iter().all(|e| (0.5..=1.6).contains(e));
    let mono_ok = rates.iter().all(|r| *r > 1.2);
    verdict(
        "dwr-001",
        eff_ok && mono_ok,
        &format!(
            "\"detail\":\"enriched-adjoint effectivity + G3 monotonicity, two frontends\",\
             \"rows\":[{}],\"eta_rates\":{rates:?}",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ dwr-002

#[test]
fn dwr_002_marking_determinism_and_prefix() {
    let disk = Circle {
        center: [0.5, 0.5],
        radius: 0.35,
    };
    let goal = GoalContext {
        weight: &goal_weight,
    };
    let grid = Quadtree::uniform(4);
    let a = estimate(&grid, &disk, FemParams::default(), &mms_f, &mms_u, &goal).expect("est a");
    let b = estimate(&grid, &disk, FemParams::default(), &mms_f, &mms_u, &goal).expect("est b");
    let bitwise = a
        .indicators
        .iter()
        .zip(&b.indicators)
        .all(|((ka, va), (kb, vb))| ka == kb && va.to_bits() == vb.to_bits());
    let marked_a = dorfler(&a.indicators, 0.5);
    let marked_b = dorfler(&b.indicators, 0.5);
    let marks_equal = marked_a == marked_b;
    // Prefix property: marked mass ≥ θ·total, and dropping the last
    // marked cell falls below θ.
    let total: f64 = a.indicators.values().map(|v| v.abs()).sum();
    let mass: f64 = marked_a.iter().map(|c| a.indicators[c].abs()).sum();
    let last = marked_a.last().expect("nonempty marking");
    let minimal = mass - a.indicators[last].abs() < 0.5 * total;
    let pass = bitwise && marks_equal && mass >= 0.5 * total && minimal;
    verdict(
        "dwr-002",
        pass,
        &format!(
            "\"detail\":\"bitwise-deterministic Doerfler marking, minimal prefix\",\
             \"marked\":{},\"mass_fraction\":{:.3},\"bitwise\":{bitwise}",
            marked_a.len(),
            mass / total
        ),
    );
}

// ------------------------------------------------------------------ dwr-003

#[test]
fn dwr_003_goal_oriented_beats_uniform() {
    // Localized QoI far from a sharp source, inside a disk.
    let disk = Circle {
        center: [0.5, 0.5],
        radius: 0.42,
    };
    let src = |x: f64, y: f64| {
        let r2 = (x - 0.3) * (x - 0.3) + (y - 0.5) * (y - 0.5);
        300.0 * (-1000.0 * r2).exp()
    };
    let qoi = |x: f64, y: f64| {
        let r2 = (x - 0.72) * (x - 0.72) + (y - 0.5) * (y - 0.5);
        (-800.0 * r2).exp()
    };
    let zero = |_: f64, _: f64| 0.0;
    let goal = GoalContext { weight: &qoi };
    // Reference truth: uniform level-7 solve.
    let j_ref = {
        let grid = Quadtree::uniform(7);
        let space = fs_cutfem::Space::build(&grid, &disk, FemParams::default()).expect("ref");
        let sol = space.solve(&src, &zero).expect("ref solves");
        fs_dwr::goal_value(&space, &sol.nodal, &goal).expect("reference goal")
    };
    // Uniform ladder.
    let mut uni_rows = String::new();
    let mut uniform: Vec<(usize, f64)> = Vec::new();
    for level in [3u32, 4, 5] {
        let grid = Quadtree::uniform(level);
        let est = estimate(&grid, &disk, FemParams::default(), &src, &zero, &goal)
            .expect("uniform estimates");
        let err = (j_ref - est.j_primal).abs();
        uniform.push((est.dofs, err));
        let _ = write!(uni_rows, "{{\"dofs\":{},\"err\":{err:.3e}}},", est.dofs);
    }
    // Adaptive loop from base 3. `adapt_loop` reserves its final row for an
    // estimate-only receipt, so six rows attempt five refinement waves. The
    // unchanged uniform-5 DOF cap below prevents that extra wave from buying
    // accuracy through a larger space.
    let mut grid = Quadtree::with_room(3, 8);
    grid.refine_toward_interface(&disk, 3);
    let (steps, _) = adapt_loop(
        &mut grid,
        &disk,
        FemParams::default(),
        &src,
        &zero,
        &goal,
        0.5,
        6,
    )
    .expect("adaptive loop");
    let mut ad_rows = String::new();
    for s in &steps {
        let _ = write!(ad_rows, "{},", s.to_json());
    }
    let final_step = steps.last().expect("steps");
    let err_adapt = (j_ref - final_step.j).abs();
    let (dofs_u5, err_u5) = uniform[2];
    // Accuracy-per-DOF: strictly better accuracy at no more DOFs —
    // the measured margin (error ratio at the final step) is the
    // ledgered figure of merit.
    #[allow(clippy::cast_precision_loss)]
    let dof_fraction = final_step.dofs as f64 / dofs_u5 as f64;
    let pass = err_adapt <= 0.5 * err_u5 && dof_fraction <= 1.0;
    verdict(
        "dwr-003",
        pass,
        &format!(
            "\"detail\":\"localized QoI: adaptive vs uniform accuracy-per-DOF\",\
             \"j_ref\":{j_ref:.8e},\"uniform\":[{}],\"adaptive\":[{}],\
             \"err_adapt\":{err_adapt:.3e},\"err_uniform5\":{err_u5:.3e},\
             \"dof_fraction\":{dof_fraction:.3}",
            uni_rows.trim_end_matches(','),
            ad_rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ dwr-004

#[test]
#[allow(clippy::too_many_lines)] // synthesis + three property checks are one narrative
fn dwr_004_anisotropic_metric_synthesis() {
    // Layer field on a uniform grid.
    let grid = Quadtree::uniform(5);
    let layer = |x: f64, y: f64| {
        let _ = x;
        (20.0 * (y - 0.5)).tanh()
    };
    let mut nodal: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    let ext = grid.node_extent();
    for gi in 0..=ext {
        for gj in 0..=ext {
            let p = grid.node_pos((gi, gj));
            nodal.insert((gi, gj), layer(p[0], p[1]));
        }
    }
    let weight: BTreeMap<(u32, u32, u32), f64> = grid.leaves().map(|c| (c, 1.0)).collect();
    let target = 800.0;
    let metric = synthesize_metric(&grid, &nodal, &weight, target);
    // (a) Complexity normalization.
    let h = 1.0 / 32.0;
    let implied: f64 = metric
        .values()
        .map(|m| (m[0][0] * m[1][1] - m[0][1] * m[1][0]).max(0.0).sqrt() * h * h)
        .sum();
    // (b) Alignment + anisotropy in the layer band.
    let mut band_align = Vec::new();
    let mut band_ratio = Vec::new();
    for (c, m) in &metric {
        let (lo, hi) = grid.rect(*c);
        let yc = f64::midpoint(lo[1], hi[1]);
        if (yc - 0.5).abs() < 0.05 {
            // Dominant eigenvector of the 2×2 metric.
            let tr = m[0][0] + m[1][1];
            let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
            let disc = (0.25 * tr * tr - det).max(0.0).sqrt();
            let l1 = 0.5 * tr + disc;
            let l2 = 0.5 * tr - disc;
            let ey = if m[0][1].abs() > 1e-20 {
                let n = (l1 - m[1][1]).hypot(m[0][1]);
                m[0][1] / n
            } else if m[0][0] >= m[1][1] {
                0.0
            } else {
                1.0
            };
            band_align.push(ey.abs());
            band_ratio.push((l1 / l2.max(1e-30)).sqrt());
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let mean_align = band_align.iter().sum::<f64>() / band_align.len() as f64;
    #[allow(clippy::cast_precision_loss)]
    let mean_ratio = band_ratio.iter().sum::<f64>() / band_ratio.len() as f64;
    // (c) Metric-instantiated graded tensor mesh beats isotropic at
    // equal DOF: y-spacing ∝ M_yy^{-1/2} column-averaged.
    let n1d = 24usize;
    let mut density: Vec<f64> = Vec::new(); // per y-row of cells
    for row in 0..32u32 {
        let mut avg = 0.0;
        for col in 0..32u32 {
            avg += metric[&(5u32, col, row)][1][1].max(0.0).sqrt();
        }
        density.push(avg / 32.0);
    }
    let total_density: f64 = density.iter().sum();
    // Graded knots: CONTINUOUS inversion of the piecewise-constant
    // density (snapping to row boundaries would cap the dense-region
    // spacing at the recovery lattice and forfeit the win).
    let mut knots = vec![0.0f64];
    let mut acc = 0.0;
    #[allow(clippy::cast_precision_loss)]
    let quantum = total_density / n1d as f64;
    let mut next_target = quantum;
    for (row, d) in density.iter().enumerate() {
        #[allow(clippy::cast_precision_loss)]
        while acc + d >= next_target - 1e-12 && knots.len() < n1d {
            let frac = ((next_target - acc) / d).clamp(0.0, 1.0);
            knots.push((row as f64 + frac) / 32.0);
            next_target += quantum;
        }
        acc += d;
    }
    knots.push(1.0);
    let interp_err = |ys: &[f64]| -> f64 {
        // Max midpoint interpolation error of the layer in y.
        let mut worst = 0.0f64;
        for w in ys.windows(2) {
            let (a, b) = (w[0], w[1]);
            let mid = f64::midpoint(a, b);
            let lin = f64::midpoint(layer(0.5, a), layer(0.5, b));
            worst = worst.max((layer(0.5, mid) - lin).abs());
        }
        worst
    };
    #[allow(clippy::cast_precision_loss)]
    let iso: Vec<f64> = (0..=n1d).map(|k| k as f64 / n1d as f64).collect();
    let err_graded = interp_err(&knots);
    let err_iso = interp_err(&iso);
    let pass = (implied - target).abs() / target < 0.05
        && mean_align > 0.9
        && mean_ratio > 5.0
        && err_graded < 0.5 * err_iso;
    verdict(
        "dwr-004",
        pass,
        &format!(
            "\"detail\":\"metric synthesis: complexity, alignment, graded-vs-iso interpolation\",\
             \"implied_cells\":{implied:.1},\"target\":{target},\
             \"mean_layer_alignment\":{mean_align:.3},\"mean_anisotropy\":{mean_ratio:.1},\
             \"graded_interp_err\":{err_graded:.3e},\"iso_interp_err\":{err_iso:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ dwr-005

#[test]
fn dwr_005_weighted_tile_thresholding() {
    let n = 64usize;
    #[allow(clippy::cast_precision_loss)]
    let cell_center = |k: usize| ((k as f64) + 0.5) / n as f64;
    let field: Vec<f64> = (0..n * n)
        .map(|k| {
            let (x, y) = (cell_center(k % n), cell_center(k / n));
            (2.0 * PI * x).sin() * (2.0 * PI * y).sin()
                + 0.2 * (31.0 * PI * (x + 0.31)).sin() * (33.0 * PI * (y + 0.17)).sin()
        })
        .collect();
    // Goal lives on the left strip; the adjoint importance decays away
    // from it.
    // Non-separable positive window so dropped details cannot cancel.
    let jw = |x: f64| if x < 0.25 { 1.0 + 4.0 * x } else { 0.0 };
    let importance = |x: f64| (-8.0 * (x - 0.125).max(0.0)).exp();
    let goal_of = |f: &[f64]| -> f64 {
        let mut j = 0.0;
        for (k, v) in f.iter().enumerate() {
            let x = cell_center(k % n);
            #[allow(clippy::cast_precision_loss)]
            let da = 1.0 / (n * n) as f64;
            j += jw(x) * v * da;
        }
        j
    };
    let j0 = goal_of(&field);
    let eps = 0.02;
    let weighted = haar_threshold(&field, n, &|ci, _| eps / importance(cell_center(ci)));
    let dj_weighted = (goal_of(&weighted.field) - j0).abs();
    // Unweighted at MATCHED compression: bisect a flat budget until the
    // kept-count matches.
    let mut lo = 1e-6;
    let mut hi = 10.0;
    let mut flat = haar_threshold(&field, n, &|_, _| eps);
    for _ in 0..40 {
        let mid = f64::midpoint(lo, hi);
        flat = haar_threshold(&field, n, &|_, _| mid);
        if flat.kept > weighted.kept {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let dj_flat = (goal_of(&flat.field) - j0).abs();
    let pass =
        weighted.ratio() >= 5.0 && dj_weighted <= eps && dj_flat > 2.0 * dj_weighted.max(1e-12);
    verdict(
        "dwr-005",
        pass,
        &format!(
            "\"detail\":\"DWR-weighted Haar budgets: goal-safe compression\",\
             \"ratio\":{:.1},\"kept\":{}/{},\"dj_weighted\":{dj_weighted:.3e},\
             \"dj_flat_matched\":{dj_flat:.3e},\"budget\":{eps}",
            weighted.ratio(),
            weighted.kept,
            weighted.total
        ),
    );
}

// ------------------------------------------------------------------ dwr-006

#[test]
fn dwr_006_h_vs_p_decision_signal() {
    let grid = Quadtree::uniform(5);
    let f = |x: f64, y: f64| 0.2 * (PI * x).sin() * (PI * y).sin() + 2.0 * (x - 0.6).abs();
    let mut nodal: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    let ext = grid.node_extent();
    for gi in 0..=ext {
        for gj in 0..=ext {
            let p = grid.node_pos((gi, gj));
            nodal.insert((gi, gj), f(p[0], p[1]));
        }
    }
    let decisions = h_vs_p(&grid, &nodal, 0.3);
    let h = 1.0 / 32.0;
    let mut kink_h = 0usize;
    let mut kink_total = 0usize;
    let mut smooth_p = 0usize;
    let mut smooth_total = 0usize;
    for (c, d) in &decisions {
        let (lo, hi) = grid.rect(*c);
        let xc = f64::midpoint(lo[0], hi[0]);
        if (xc - 0.6).abs() < h {
            kink_total += 1;
            if *d == Decision::HRefine {
                kink_h += 1;
            }
        } else if (xc - 0.6).abs() > 0.15 {
            smooth_total += 1;
            if *d == Decision::PEnrich {
                smooth_p += 1;
            }
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let kink_frac = kink_h as f64 / kink_total.max(1) as f64;
    #[allow(clippy::cast_precision_loss)]
    let smooth_frac = smooth_p as f64 / smooth_total.max(1) as f64;
    let pass = kink_frac > 0.9 && smooth_frac > 0.9;
    verdict(
        "dwr-006",
        pass,
        &format!(
            "\"detail\":\"kinks route to h, smooth to p (execution awaits local-p spaces)\",\
             \"kink_h_fraction\":{kink_frac:.3},\"smooth_p_fraction\":{smooth_frac:.3},\
             \"cells\":{}",
            decisions.len()
        ),
    );
}
