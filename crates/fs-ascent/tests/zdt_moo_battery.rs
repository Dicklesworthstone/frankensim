//! Bead 7tv.21.4: MOO standards — ZDT1/ZDT2 front conformance for the
//! Pareto tracing sweeps.
//!
//! The sweeps are 2-objective gradient tracers with no box support, so
//! ZDT's [0,1]^n box is entered through the smooth chart
//! x_i = sin²(θ_i): surjective onto the box and front-geometry
//! preserving, with every composite differentiable along the sweep (the
//! |sin| kink of ZDT1's √x₁ sits at the measure-zero chart seams
//! θ₁ = kπ, which the ε grid keeps away from). n = 3 chart variables —
//! front conformance is n-independent; scalability is separate parent
//! scope.

use fs_ascent::pareto::{ParetoPoint, epsilon_constraint_sweep, weighted_sum_sweep};

const SUITE: &str = "fs-ascent/zdt-moo";

fn verdict(name: &str, pass: bool, details: &str) {
    // Canonical fs-obs verdict (51bbx suite adoption): linted, wire-validated,
    // emitted BEFORE the assertion so failures reproduce from the log alone.
    let mut emitter = fs_obs::Emitter::new(SUITE, name);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: name.to_string(),
            pass,
            detail: details.to_string(),
            seed: 0, // deterministic fixed-input batteries; embedded seeds live in detail
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("verdict must satisfy the failure-record lint");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "{name}: {details}");
}

const N: usize = 3;

/// Chart: x_i = sin²(θ_i) ∈ [0, 1]; dx/dθ = sin(2θ).
fn chart(theta: &[f64]) -> Vec<f64> {
    theta.iter().map(|t| t.sin() * t.sin()).collect()
}

fn chart_jac(theta: &[f64]) -> Vec<f64> {
    theta.iter().map(|t| (2.0 * t).sin()).collect()
}

/// g(x) = 1 + 9·Σ_{i≥2} x_i / (n−1); minimal (true front) at g = 1.
fn g_of(x: &[f64]) -> f64 {
    1.0 + 9.0 * x[1..].iter().sum::<f64>() / (N as f64 - 1.0)
}

/// ZDT1 in chart coordinates: f1 = x1, f2 = g·(1 − √(x1/g)).
fn zdt1_f1(theta: &[f64]) -> (f64, Vec<f64>) {
    let x = chart(theta);
    let j = chart_jac(theta);
    let mut grad = vec![0.0; N];
    grad[0] = j[0];
    (x[0], grad)
}

fn zdt1_f2(theta: &[f64]) -> (f64, Vec<f64>) {
    let x = chart(theta);
    let j = chart_jac(theta);
    let g = g_of(&x);
    let ratio = (x[0] / g).max(1e-300);
    let root = ratio.sqrt();
    let f2 = g * (1.0 - root);
    // From f2 = g − √x1·√g with ∂g/∂x_i = c (i ≥ 2):
    //   ∂f2/∂x1 = −√g/(2√x1);  ∂f2/∂x_i = c − √x1·c/(2√g).
    let c = 9.0 / (N as f64 - 1.0);
    let sx1 = x[0].max(1e-300).sqrt();
    let mut grad = vec![0.0; N];
    grad[0] = -(g.sqrt()) / (2.0 * sx1) * j[0];
    for i in 1..N {
        grad[i] = (c - sx1 * c / (2.0 * g.sqrt())) * j[i];
    }
    (f2, grad)
}

/// ZDT2 in chart coordinates: f1 = x1, f2 = g·(1 − (x1/g)²).
fn zdt2_f1(theta: &[f64]) -> (f64, Vec<f64>) {
    zdt1_f1(theta)
}

fn zdt2_f2(theta: &[f64]) -> (f64, Vec<f64>) {
    let x = chart(theta);
    let j = chart_jac(theta);
    let g = g_of(&x);
    let f2 = g - x[0] * x[0] / g;
    let c = 9.0 / (N as f64 - 1.0);
    let mut grad = vec![0.0; N];
    grad[0] = -2.0 * x[0] / g * j[0];
    for i in 1..N {
        grad[i] = (c + x[0] * x[0] / (g * g) * c) * j[i];
    }
    (f2, grad)
}

/// Discrete 2-D hypervolume against reference (1.1, 1.1) of a mutually
/// non-dominated set sorted by f1 ascending (f2 descending).
fn hypervolume(points: &[[f64; 2]]) -> f64 {
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let (ref1, ref2) = (1.1f64, 1.1f64);
    let mut hv = 0.0;
    for (k, p) in pts.iter().enumerate() {
        let next_x = pts.get(k + 1).map_or(ref1, |q| q[0]);
        hv += (next_x - p[0]) * (ref2 - p[1]);
    }
    hv
}

/// Shared conformance walk over one traced front.
struct FrontReport {
    worst_front_err: f64,
    worst_g_excess: f64,
    worst_kkt: f64,
    f1_min: f64,
    f1_max: f64,
    hv_traced: f64,
    hv_analytic: f64,
}

fn conform(points: &[ParetoPoint], analytic_f2: impl Fn(f64) -> f64) -> FrontReport {
    let mut worst_front_err = 0.0f64;
    let mut worst_g_excess = 0.0f64;
    let mut worst_kkt = 0.0f64;
    let mut f1_min = f64::INFINITY;
    let mut f1_max = f64::NEG_INFINITY;
    let mut traced = Vec::new();
    let mut analytic = Vec::new();
    for p in points {
        let [f1, f2] = p.f;
        worst_front_err = worst_front_err.max((f2 - analytic_f2(f1)).abs());
        let x = chart(&p.x);
        worst_g_excess = worst_g_excess.max(g_of(&x) - 1.0);
        if let Some(kkt) = &p.kkt {
            worst_kkt = worst_kkt.max(kkt.stationarity);
        }
        f1_min = f1_min.min(f1);
        f1_max = f1_max.max(f1);
        traced.push([f1, f2]);
        analytic.push([f1, analytic_f2(f1)]);
    }
    FrontReport {
        worst_front_err,
        worst_g_excess,
        worst_kkt,
        f1_min,
        f1_max,
        hv_traced: hypervolume(&traced),
        hv_analytic: hypervolume(&analytic),
    }
}

fn epsilons() -> Vec<f64> {
    (1..=19).map(|k| f64::from(k) * 0.05).collect()
}

#[test]
fn zdt1_convex_front_conformance() {
    // Start on the box interior, away from chart seams.
    let theta0 = [0.8f64, 0.3, 0.3];
    let points = epsilon_constraint_sweep(&zdt1_f1, &zdt1_f2, &epsilons(), &theta0, 1e-8);
    let r = conform(&points, |f1| 1.0 - f1.sqrt());
    // Tolerances: front error 5e-4 (auglag tol 1e-8 in θ maps to ~1e-4
    // in f near the √ singular end), g-collapse 1e-4, KKT 1e-4,
    // coverage must span [0.05, 0.95], HV parity 1e-3.
    assert!(
        r.worst_front_err < 5e-4,
        "zdt1: front error {:.3e}",
        r.worst_front_err
    );
    assert!(
        r.worst_g_excess < 1e-4,
        "zdt1: traced points did not collapse g to 1 (excess {:.3e})",
        r.worst_g_excess
    );
    assert!(r.worst_kkt < 1e-4, "zdt1: KKT residual {:.3e}", r.worst_kkt);
    assert!(
        r.f1_min < 0.06 && r.f1_max > 0.94,
        "zdt1: coverage [{:.3}, {:.3}] misses the sweep span",
        r.f1_min,
        r.f1_max
    );
    let hv_dev = (r.hv_traced - r.hv_analytic).abs();
    assert!(hv_dev < 1e-3, "zdt1: hypervolume deviation {hv_dev:.3e}");
    verdict(
        "7tv21-zdt1",
        true,
        &format!(
            "convex front 19 eps-points: front-err {:.2e} (tol 5e-4), g-excess {:.2e}, \
             KKT {:.2e}, coverage [{:.2},{:.2}], HV dev {:.2e} (tol 1e-3)",
            r.worst_front_err, r.worst_g_excess, r.worst_kkt, r.f1_min, r.f1_max, hv_dev
        ),
    );
}

#[test]
fn zdt2_concave_front_epsilon_covers_where_weighted_sum_collapses() {
    let theta0 = [0.8f64, 0.3, 0.3];
    let points = epsilon_constraint_sweep(&zdt2_f1, &zdt2_f2, &epsilons(), &theta0, 1e-9);
    let r = conform(&points, |f1| 1.0 - f1 * f1);
    assert!(
        r.worst_front_err < 5e-4,
        "zdt2: front error {:.3e}",
        r.worst_front_err
    );
    assert!(
        r.worst_g_excess < 1e-4,
        "zdt2: traced points did not collapse g to 1 (excess {:.3e})",
        r.worst_g_excess
    );
    assert!(r.worst_kkt < 1e-4, "zdt2: KKT residual {:.3e}", r.worst_kkt);
    assert!(
        r.f1_min < 0.06 && r.f1_max > 0.94,
        "zdt2: coverage [{:.3}, {:.3}] misses the sweep span",
        r.f1_min,
        r.f1_max
    );
    let hv_dev = (r.hv_traced - r.hv_analytic).abs();
    assert!(hv_dev < 1e-3, "zdt2: hypervolume deviation {hv_dev:.3e}");

    // The module's documented weighted-sum-collapse claim, exhibited on
    // a STANDARD concave fixture: interior weights land only at the
    // front's extremes (f1 ≈ 0 or f1 ≈ 1), never on the interior.
    let weights: Vec<f64> = (1..=9).map(|k| f64::from(k) * 0.1).collect();
    let ws = weighted_sum_sweep(&zdt2_f1, &zdt2_f2, &weights, &theta0);
    let interior = ws.iter().filter(|p| p.f[0] > 0.05 && p.f[0] < 0.95).count();
    assert_eq!(
        interior, 0,
        "zdt2: weighted sums must collapse to the extremes on a concave front"
    );
    verdict(
        "7tv21-zdt2",
        true,
        &format!(
            "concave front 19 eps-points: front-err {:.2e}, g-excess {:.2e}, KKT {:.2e}, \
             coverage [{:.2},{:.2}], HV dev {:.2e}; weighted-sum collapse exhibited \
             (0/9 interior weights reach the interior front)",
            r.worst_front_err, r.worst_g_excess, r.worst_kkt, r.f1_min, r.f1_max, hv_dev
        ),
    );
}
