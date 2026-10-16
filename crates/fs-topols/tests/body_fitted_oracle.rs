//! Independent body-fitted elasticity oracle for the level-set compliance study
//! (bead q61wp.75 item 4).
//!
//! The optimizer and its final evaluation solve on CutFEM: cut quadrature,
//! Nitsche clamping and ghost-penalty stabilization over a fixed lattice. This
//! oracle shares none of that. It clips every sub-lattice triangle by the
//! linearly interpolated level set into a conforming body-fitted P1 mesh and
//! solves plane-strain elasticity with fs-solid's standard displacement
//! elements, strong clamping and edge-integrated traction.
//!
//! The design's bilinear field reproduces exactly under dyadic refinement, so
//! both methods can be refined on the identical geometry. The acceptance band
//! is derived from each method's own refinement change, not chosen. A mutated
//! (eroded) design must fall outside it.

use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_topols::{Cantilever, GridSdf, OptimizeSettings, evaluate_compliance_design, optimize_compliance};
use std::collections::BTreeMap;

const FIXTURE: Cantilever = Cantilever { load: 1.0, band: 0.125 };

fn settings(level: u32) -> OptimizeSettings {
    OptimizeSettings { level, iterations: 12, ..OptimizeSettings::default() }
}

/// The same bilinear geometry on a lattice `factor` times finer.
fn refined(phi: &GridSdf, factor: usize) -> GridSdf {
    GridSdf::from_fn(phi.n() * factor, &|x, y| phi.value_at([x, y]))
}

/// Conforming body-fitted P1 mesh of `{phi < 0}` on an `n x n` sub-lattice.
fn body_fitted(phi: &GridSdf, n: usize) -> Mesh2 {
    let h = 1.0 / n as f64;
    let scale = phi.nodes().iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let value = |i: usize, j: usize| {
        let v = phi.value_at([i as f64 * h, j as f64 * h]);
        // An exact zero would make a zero-area element; treat it as void.
        if v == 0.0 { 1e-14 * scale } else { v }
    };
    let lattice: Vec<f64> = (0..=n).flat_map(|j| (0..=n).map(move |i| (i, j))).map(|(i, j)| value(i, j)).collect();
    let id = |i: usize, j: usize| j * (n + 1) + i;
    let point = |k: usize| [(k % (n + 1)) as f64 * h, (k / (n + 1)) as f64 * h];
    let mut nodes: Vec<[f64; 2]> = Vec::new();
    let mut lattice_node: BTreeMap<usize, usize> = BTreeMap::new();
    let mut cut_node: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut elems: Vec<Vec<usize>> = Vec::new();
    let mut vertex = |k: usize, nodes: &mut Vec<[f64; 2]>| {
        *lattice_node.entry(k).or_insert_with(|| {
            nodes.push(point(k));
            nodes.len() - 1
        })
    };
    for j in 0..n {
        for i in 0..n {
            let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
            for tri in [[a, b, c], [a, c, d]] {
                // Sutherland-Hodgman clip of a counterclockwise triangle.
                let mut polygon = Vec::with_capacity(4);
                for e in 0..3 {
                    let (p, q) = (tri[e], tri[(e + 1) % 3]);
                    let (fp, fq) = (lattice[p], lattice[q]);
                    if fp < 0.0 {
                        polygon.push(vertex(p, &mut nodes));
                    }
                    if (fp < 0.0) != (fq < 0.0) {
                        let key = (p.min(q), p.max(q));
                        let node = *cut_node.entry(key).or_insert_with(|| {
                            let t = fp / (fp - fq);
                            let (xp, xq) = (point(p), point(q));
                            nodes.push([xp[0] + t * (xq[0] - xp[0]), xp[1] + t * (xq[1] - xp[1])]);
                            nodes.len() - 1
                        });
                        polygon.push(node);
                    }
                }
                for k in 1..polygon.len().saturating_sub(1) {
                    let (p, q, r) = (nodes[polygon[0]], nodes[polygon[k]], nodes[polygon[k + 1]]);
                    let area = 0.5 * ((q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1]));
                    if area > 1e-12 * h * h {
                        elems.push(vec![polygon[0], polygon[k], polygon[k + 1]]);
                    }
                }
            }
        }
    }
    let mut patches = Vec::new();
    for (patch, x) in [(Patch::Left, 0.0), (Patch::Right, 1.0)] {
        let mut edges = Vec::new();
        for elem in &elems {
            for e in 0..3 {
                let (p, q) = (elem[e], elem[(e + 1) % 3]);
                if nodes[p][0] == x && nodes[q][0] == x {
                    edges.push((p, q));
                }
            }
        }
        patches.push((patch, edges));
    }
    Mesh2 { nodes, elems, patches }
}

/// Compliance of the body-fitted P1 solve: the load's work on its edges.
fn oracle(phi: &GridSdf, n: usize, material: (f64, f64)) -> f64 {
    let mesh = body_fitted(phi, n);
    let band = |y: f64| (y - 0.5).abs() <= FIXTURE.band + 1e-12;
    let traction = |_: f64, y: f64| if band(y) { [0.0, -FIXTURE.load] } else { [0.0, 0.0] };
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: material.0,
        poisson: material.1,
        plane: PlaneKind::Strain,
        formulation: Formulation::Standard,
        body_force: None,
        dirichlet: vec![(Patch::Left, &zero)],
        traction: vec![(Patch::Right, &traction)],
        symmetry: Vec::new(),
    };
    let u = problem.solve().expect("body-fitted solve");
    // P1 displacement is linear on each edge; the band ends fall on lattice
    // nodes, so the traction is constant per loaded edge and this is exact.
    mesh.patch_edges(Patch::Right)
        .unwrap()
        .iter()
        .filter(|&&(p, q)| band(0.5 * (mesh.nodes[p][1] + mesh.nodes[q][1])))
        .map(|&(p, q)| {
            let length = (mesh.nodes[q][1] - mesh.nodes[p][1]).abs();
            -FIXTURE.load * 0.5 * (u[p][1] + u[q][1]) * length
        })
        .sum()
}

/// Observed order (clamped to [0.5, 3]) and Richardson limit of three
/// successive dyadic refinements.
fn richardson(v: [f64; 3]) -> (f64, f64) {
    let order = ((v[1] - v[0]).abs() / (v[2] - v[1]).abs()).log2().clamp(0.5, 3.0);
    (order, v[2] + (v[2] - v[1]) / (2f64.powf(order) - 1.0))
}

fn cutfem(phi: &GridSdf, factor: usize, base: &OptimizeSettings) -> f64 {
    let fine = refined(phi, factor);
    let level = base.level + factor.trailing_zeros();
    evaluate_compliance_design(&fine, FIXTURE, OptimizeSettings { level, ..*base })
        .expect("CutFEM evaluation")
        .compliance
}

#[test]
fn body_fitted_p1_and_cutfem_agree_on_the_optimized_design_within_their_refinement_band() {
    let base = settings(4);
    let n = 1usize << base.level;
    let mut phi = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.35);
    optimize_compliance(&mut phi, FIXTURE, base).expect("optimizer runs");
    let material = (base.youngs, base.poisson);

    let cut = [cutfem(&phi, 1, &base), cutfem(&phi, 2, &base), cutfem(&phi, 4, &base)];
    let fit = [oracle(&phi, 2 * n, material), oracle(&phi, 4 * n, material), oracle(&phi, 8 * n, material)];
    assert!(fit.iter().chain(&cut).all(|c| c.is_finite() && *c > 0.0));
    let (cut_order, cut_limit) = richardson(cut);
    let (fit_order, fit_limit) = richardson(fit);
    // Each limit is uncertain by about its own extrapolation correction; the
    // band is their sum. Measured 2026-09-25: gap 0.0124 against band 0.064
    // (about 5x headroom). CutFEM converges at observed order ~1.0 here and
    // the body-fitted P1 solve at ~1.4.
    let band = (cut_limit - cut[2]).abs() + (fit_limit - fit[2]).abs();
    let gap = (cut_limit - fit_limit).abs();
    println!(
        "{{\"cutfem\":{cut:?},\"body_fitted\":{fit:?},\"orders\":[{cut_order:.3},{fit_order:.3}],\"limits\":[{cut_limit},{fit_limit}],\"gap\":{gap:e},\"band\":{band:e}}}"
    );
    assert!(gap <= band, "independent methods disagree: gap {gap:e} > band {band:e}");

    // Mutation: erode the design by a quarter lattice cell. Its extrapolated
    // body-fitted compliance must leave the band around the CutFEM limit.
    let eroded = GridSdf::from_fn(n, &|x, y| phi.value_at([x, y]) + 0.25 / n as f64);
    let (_, mutated) = richardson([
        oracle(&eroded, 2 * n, material),
        oracle(&eroded, 4 * n, material),
        oracle(&eroded, 8 * n, material),
    ]);
    println!("{{\"eroded_limit\":{mutated},\"separation\":{:e}}}", (cut_limit - mutated).abs());
    assert!(
        (cut_limit - mutated).abs() > band,
        "the band cannot tell a quarter-cell erosion apart: {:e} <= {band:e}",
        (cut_limit - mutated).abs()
    );
}
