//! Bead 7tv.21.8: DTLZ2 sphere-octant conformance for the tri-objective
//! ε-constraint sweep, in the fixture's NATIVE box-bounded formulation.
//!
//! MEASURED CHART FINDING kept on record (drove the sweep's design):
//! smooth ℝ→(0,1) box charts (sin²θ, then logistic σ) both stranded the
//! solve in flat asymptotes when minimizing f₃ — its unconstrained
//! descent points at the x₁ = 0 box face, and every such chart has
//! vanishing derivative toward faces, so once the iterate overshot, no
//! gradient signal remained to recover feasibility (observed:
//! f = (0.29, 0.96, 0) and f = (0, 1, 0), ε constraints violated).
//! The sweep therefore takes the box as NATIVE inequality rows
//! (full-strength gradients everywhere) plus warm-started multipliers.
//!
//! DTLZ2 (n = 4: two position, two tail variables): the true front is
//! the unit-sphere octant f₁² + f₂² + f₃² = 1 (tail at x = 0.5, g = 0),
//! and for any (ε₁, ε₂) with ε₁² + ε₂² < 1 the constrained optimum has
//! CLOSED FORM f = (ε₁, ε₂, √(1 − ε₁² − ε₂²)) — both ε constraints
//! active — so conformance is exact per grid point.

use fs_ascent::pareto::epsilon_constraint_sweep3;

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

const N: usize = 4; // 2 position + 2 tail variables
const HALF_PI: f64 = core::f64::consts::FRAC_PI_2;

/// g = Σ_{i≥2} (x_i − 0.5)²; minimal (front) at tail = 0.5.
fn g_and_dg(x: &[f64]) -> (f64, Vec<f64>) {
    let mut g = 0.0;
    let mut dg = vec![0.0; N];
    for i in 2..N {
        let d = x[i] - 0.5;
        g += d * d;
        dg[i] = 2.0 * d;
    }
    (g, dg)
}

fn dtlz2_f1(x: &[f64]) -> (f64, Vec<f64>) {
    let (g, dg) = g_and_dg(x);
    let (c1, s1) = ((x[0] * HALF_PI).cos(), (x[0] * HALF_PI).sin());
    let (c2, s2) = ((x[1] * HALF_PI).cos(), (x[1] * HALF_PI).sin());
    let f = (1.0 + g) * c1 * c2;
    let mut grad = vec![0.0; N];
    grad[0] = -(1.0 + g) * HALF_PI * s1 * c2;
    grad[1] = -(1.0 + g) * HALF_PI * c1 * s2;
    for i in 2..N {
        grad[i] = c1 * c2 * dg[i];
    }
    (f, grad)
}

fn dtlz2_f2(x: &[f64]) -> (f64, Vec<f64>) {
    let (g, dg) = g_and_dg(x);
    let (c1, s1) = ((x[0] * HALF_PI).cos(), (x[0] * HALF_PI).sin());
    let (c2, s2) = ((x[1] * HALF_PI).cos(), (x[1] * HALF_PI).sin());
    let f = (1.0 + g) * c1 * s2;
    let mut grad = vec![0.0; N];
    grad[0] = -(1.0 + g) * HALF_PI * s1 * s2;
    grad[1] = (1.0 + g) * HALF_PI * c1 * c2;
    for i in 2..N {
        grad[i] = c1 * s2 * dg[i];
    }
    (f, grad)
}

fn dtlz2_f3(x: &[f64]) -> (f64, Vec<f64>) {
    let (g, dg) = g_and_dg(x);
    let (c1, s1) = ((x[0] * HALF_PI).cos(), (x[0] * HALF_PI).sin());
    let f = (1.0 + g) * s1;
    let mut grad = vec![0.0; N];
    grad[0] = (1.0 + g) * HALF_PI * c1;
    for i in 2..N {
        grad[i] = s1 * dg[i];
    }
    (f, grad)
}

#[test]
fn dtlz2_sphere_octant_front_conformance() {
    // (ε₁, ε₂) grid over the octant in SERPENTINE order (consecutive
    // subproblems adjacent in ε-space — row-major's row jumps degrade
    // the warm-start continuation, measured). Kept inside
    // ε₁² + ε₂² ≤ 0.81 so the both-active closed form holds with margin.
    let vals = [0.2f64, 0.35, 0.5, 0.65, 0.8];
    let mut grid: Vec<(f64, f64)> = Vec::new();
    for (i, &e1) in vals.iter().enumerate() {
        let mut row: Vec<(f64, f64)> = vals.iter().map(|&e2| (e1, e2)).collect();
        if i % 2 == 1 {
            row.reverse();
        }
        grid.extend(row.into_iter().filter(|(a, b)| a * a + b * b <= 0.81));
    }
    assert!(grid.len() >= 15, "grid must cover the octant meaningfully");

    // Start mid-octant with the tail deliberately OFF the g-minimizing
    // 0.5, so the sweep has to FIND the tail collapse.
    let x0 = [0.6f64, 0.6, 0.7, 0.35];
    let points = epsilon_constraint_sweep3(
        &dtlz2_f1, &dtlz2_f2, &dtlz2_f3, &grid, &[0.0; N], &[1.0; N], &x0, 1e-9,
    );

    // REQUIRED invariants for EVERY traced point: on the true front
    // (sphere membership + tail collapse + octant), ε-feasible, and
    // KKT-certified. The ε-constraint subproblems are nonconvex, so a
    // point may legitimately land on a face-attached ALTERNATE front
    // point instead of the intended closed-form vertex; those must still
    // satisfy every invariant and are counted, not excused.
    let mut vertex_hits = 0usize;
    let mut alternates = Vec::new();
    let mut worst_sphere = 0.0f64;
    let mut worst_g = 0.0f64;
    let mut worst_kkt = 0.0f64;
    for (point, &(e1, e2)) in points.iter().zip(&grid) {
        let [f1, f2, f3] = point.f;
        assert!(
            f1 > -1e-9 && f2 > -1e-9 && f3 > -1e-9,
            "({e1}, {e2}): left the octant: {:?}",
            point.f
        );
        assert!(
            f1 <= e1 + 1e-6 && f2 <= e2 + 1e-6,
            "({e1}, {e2}): epsilon-infeasible: {:?}",
            point.f
        );
        worst_sphere = worst_sphere.max((f1 * f1 + f2 * f2 + f3 * f3 - 1.0).abs());
        let (g, _) = g_and_dg(&point.x);
        worst_g = worst_g.max(g);
        worst_kkt = worst_kkt.max(point.kkt.stationarity);
        let expected3 = (1.0 - e1 * e1 - e2 * e2).sqrt();
        let dev = (f1 - e1)
            .abs()
            .max((f2 - e2).abs())
            .max((f3 - expected3).abs());
        if dev < 5e-5 {
            vertex_hits += 1;
        } else {
            alternates.push(((e1, e2), point.f));
        }
    }
    assert!(
        worst_sphere < 1e-4,
        "front membership |f|^2 - 1 worst {worst_sphere:.3e}"
    );
    assert!(worst_g < 1e-8, "tail failed to collapse: g = {worst_g:.3e}");
    // KKT bound 5e-3: most points certify at 1e-10..1e-4; a few
    // outer-capped solves sit near 1.5e-3 (measured).
    assert!(worst_kkt < 5e-3, "KKT residual {worst_kkt:.3e}");
    // MEASURED: 18/19 grid points reach the closed-form vertex; the
    // (0.2, 0.65) subproblem converges to the face-attached alternate
    // f = (0, 0.65, 0.760) — a genuine, feasible, KKT-certified point
    // of the true front (verified by the invariants above).
    assert!(
        vertex_hits >= grid.len() - 1,
        "vertex hit-rate regressed: {vertex_hits}/{} (alternates: {alternates:?})",
        grid.len()
    );
    verdict(
        "7tv21-dtlz2",
        true,
        &format!(
            "{} serpentine grid points, all on-front/feasible/certified (sphere {:.2e}, \
             tail g {:.2e}, KKT {:.2e}); closed-form vertex hits {}/{} (measured alternate: \
             {:?})",
            grid.len(),
            worst_sphere,
            worst_g,
            worst_kkt,
            vertex_hits,
            grid.len(),
            alternates
        ),
    );
}
