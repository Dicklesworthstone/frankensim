//! Pareto-tracing battery (vcia lane c): the convex quadratic pair
//! with its CLOSED-FORM front (weighted-sum points match
//! analytically); the Fonseca–Fleming CONCAVE front where weighted
//! sums provably collapse to the extremes (measured cluster count)
//! while ε-constraint covers the middle with KKT certificates on
//! every point; bitwise replay; golden.

use fs_ascent::{epsilon_constraint_sweep, weighted_sum_sweep};
use std::panic::{AssertUnwindSafe, catch_unwind};

type ObjPair = (
    Box<dyn Fn(&[f64]) -> (f64, Vec<f64>)>,
    Box<dyn Fn(&[f64]) -> (f64, Vec<f64>)>,
);

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ascent-pareto\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

/// f₁ = ‖x − a‖², f₂ = ‖x − b‖² with a = 0, b = 1 (3-d):
/// x*(w) = (1−w)·b + w·a ⇒ f₁ = 3(1−w)², f₂ = 3w².
fn quad_pair() -> ObjPair {
    let f1 = |x: &[f64]| -> (f64, Vec<f64>) {
        let v: f64 = x.iter().map(|t| t * t).sum();
        (v, x.iter().map(|t| 2.0 * t).collect())
    };
    let f2 = |x: &[f64]| -> (f64, Vec<f64>) {
        let v: f64 = x.iter().map(|t| (t - 1.0) * (t - 1.0)).sum();
        (v, x.iter().map(|t| 2.0 * (t - 1.0)).collect())
    };
    (Box::new(f1), Box::new(f2))
}

/// Fonseca–Fleming (2-d): smooth, CONCAVE front.
fn ff_pair() -> ObjPair {
    let c = 1.0 / fs_math::det::sqrt(2.0);
    let f1 = move |x: &[f64]| -> (f64, Vec<f64>) {
        let s: f64 = x.iter().map(|t| (t - c) * (t - c)).sum();
        let e = fs_math::det::exp(-s);
        (1.0 - e, x.iter().map(|t| 2.0 * (t - c) * e).collect())
    };
    let f2 = move |x: &[f64]| -> (f64, Vec<f64>) {
        let s: f64 = x.iter().map(|t| (t + c) * (t + c)).sum();
        let e = fs_math::det::exp(-s);
        (1.0 - e, x.iter().map(|t| 2.0 * (t + c) * e).collect())
    };
    (Box::new(f1), Box::new(f2))
}

#[test]
fn weighted_sum_matches_closed_form_on_convex() {
    let (f1, f2) = quad_pair();
    let weights: Vec<f64> = (1..10).map(|k| f64::from(k) / 10.0).collect();
    let front = weighted_sum_sweep(&f1, &f2, &weights, &[0.5, 0.5, 0.5]);
    let mut worst = 0.0f64;
    for (p, &w) in front.iter().zip(&weights) {
        let f1_true = 3.0 * (1.0 - w) * (1.0 - w);
        let f2_true = 3.0 * w * w;
        worst = worst
            .max((p.f[0] - f1_true).abs())
            .max((p.f[1] - f2_true).abs());
        assert!(p.grad_norm < 1e-9, "scalarized solve not certified");
    }
    assert!(
        worst < 1e-7,
        "closed-form front missed: worst dev {worst:.2e}"
    );
    log(
        "convex-closed-form",
        "pass",
        &format!("worst dev {worst:.1e}"),
    );
}

#[test]
fn concave_front_epsilon_covers_weighted_collapses() {
    let (f1, f2) = ff_pair();
    // Weighted sums on the CONCAVE front: sweep 9 weights, count
    // distinct clusters in f-space (radius 0.05). The classic result:
    // everything lands at the two extremes.
    let weights: Vec<f64> = (1..10).map(|k| f64::from(k) / 10.0).collect();
    let ws = weighted_sum_sweep(&f1, &f2, &weights, &[0.1, -0.1]);
    let mut clusters: Vec<[f64; 2]> = Vec::new();
    for p in &ws {
        if !clusters
            .iter()
            .any(|c| (c[0] - p.f[0]).abs() + (c[1] - p.f[1]).abs() < 0.05)
        {
            clusters.push(p.f);
        }
    }
    assert!(
        clusters.len() <= 3,
        "weighted sums should collapse on a concave front: {} clusters",
        clusters.len()
    );
    // ε-constraint: covers the middle. Trace f₁ ∈ [0.15, 0.85].
    let epsilons: Vec<f64> = (0..8).map(|k| 0.1f64.mul_add(f64::from(k), 0.15)).collect();
    let ec = epsilon_constraint_sweep(&f1, &f2, &epsilons, &[0.0, 0.0], 1e-7);
    // On-front check: the FF Pareto set is x₁ = x₂ = t, t ∈ [−c, c].
    let mut worst_off = 0.0f64;
    let mut worst_kkt = 0.0f64;
    for p in &ec {
        worst_off = worst_off.max((p.x[0] - p.x[1]).abs());
        let k = p.kkt.as_ref().expect("certificate present");
        worst_kkt = worst_kkt
            .max(k.stationarity)
            .max(k.feasibility)
            .max(k.dual_feasibility)
            .max(k.complementarity);
    }
    assert!(
        worst_off < 1e-4,
        "traced points off the Pareto set (x1 != x2): {worst_off:.2e}"
    );
    assert!(
        worst_kkt < 1e-5,
        "KKT certificates too loose: {worst_kkt:.2e}"
    );
    // Coverage: f₁ spread across the concave middle.
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in &ec {
        lo = lo.min(p.f[0]);
        hi = hi.max(p.f[0]);
    }
    assert!(
        hi - lo > 0.6,
        "epsilon-constraint should cover the front: [{lo:.3}, {hi:.3}]"
    );
    log(
        "concave-coverage",
        "pass",
        &format!(
            "weighted clusters {}, eps spread [{lo:.2},{hi:.2}], worst KKT {worst_kkt:.1e}",
            clusters.len()
        ),
    );
}

const GOLDEN_HASH: u64 = 0x301b_04df_db91_3965; // recorded at vcia lane c, frozen

#[test]
fn pareto_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (f1, f2) = quad_pair();
    let front = weighted_sum_sweep(&f1, &f2, &[0.25, 0.5, 0.75], &[0.2, 0.8, 0.5]);
    for p in &front {
        feed(p.f[0]);
        feed(p.f[1]);
    }
    let (g1, g2) = ff_pair();
    let ec = epsilon_constraint_sweep(&g1, &g2, &[0.3, 0.5, 0.7], &[0.0, 0.0], 1e-7);
    for p in &ec {
        feed(p.f[0]);
        feed(p.f[1]);
        feed(p.grad_norm);
    }
    // Determinism: repeat must be bitwise.
    let ec2 = epsilon_constraint_sweep(&g1, &g2, &[0.3, 0.5, 0.7], &[0.0, 0.0], 1e-7);
    for (a, b) in ec.iter().zip(&ec2) {
        assert!(a.f[1].to_bits() == b.f[1].to_bits(), "sweep not replayable");
    }
    log("pareto-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "pareto bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn pareto_contract_guards_fail_fast() {
    let good = |x: &[f64]| -> (f64, Vec<f64>) {
        let f = x.iter().map(|v| v * v).sum();
        (f, x.iter().map(|v| 2.0 * v).collect())
    };
    let bad_grad = |_: &[f64]| -> (f64, Vec<f64>) { (0.0, vec![0.0]) };
    let bad_value = |x: &[f64]| -> (f64, Vec<f64>) { (f64::NAN, vec![0.0; x.len()]) };
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = weighted_sum_sweep(&good, &good, &[1.25], &[0.0, 1.0]);
        }))
        .is_err(),
        "weights outside [0, 1] must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = weighted_sum_sweep(&good, &good, &[0.5], &[0.0, f64::NAN]);
        }))
        .is_err(),
        "non-finite initial decision vectors must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = weighted_sum_sweep(&bad_grad, &good, &[0.5], &[0.0, 1.0]);
        }))
        .is_err(),
        "mismatched weighted-sum gradients must fail fast instead of truncating"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = epsilon_constraint_sweep(&good, &good, &[f64::NAN], &[0.0, 1.0], 1e-6);
        }))
        .is_err(),
        "non-finite epsilon constraints must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = epsilon_constraint_sweep(&good, &good, &[0.1], &[0.0, 1.0], 0.0);
        }))
        .is_err(),
        "non-positive epsilon-constraint tolerances must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = epsilon_constraint_sweep(&good, &bad_value, &[0.1], &[0.0, 1.0], 1e-6);
        }))
        .is_err(),
        "non-finite objective values must fail fast"
    );
    log(
        "contract-guards",
        "pass",
        "invalid weights, epsilons, tolerances, decisions, gradients, and objective values fail fast",
    );
}
