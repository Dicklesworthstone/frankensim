//! fs-bo battery (7tv.5 slices 1–3): kernel PSD-ness across the
//! Matérn family, POSTERIOR CONSISTENCY (conditioning on D₁∪D₂ =
//! conditioning on the union directly — the G0 law from the
//! acceptance), known-answer posteriors (noiseless interpolation +
//! hand-computed one-point closed form), Φ/Φ⁻¹ round-trip accuracy at
//! the documented tolerances, EI sanity laws, q-EI dominance over
//! single-point EI, BO-beats-QMC-random on Branin (matched budget,
//! fixed seed set — the ledgered evidence), bitwise replay, and the
//! golden hash.

use fs_bo::{
    BoConfig, Gp, Kernel, Matern, expected_improvement, minimize, normal_bank, phi_cdf, phi_inv,
    q_expected_improvement,
};
use fs_rand::StreamKey;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-bo\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_points(n: usize, d: usize, tile: u32) -> Vec<Vec<f64>> {
    let mut s = StreamKey {
        seed: 61,
        kernel: 0xB0B0,
        tile,
    }
    .stream();
    (0..n)
        .map(|_| (0..d).map(|_| s.next_f64()).collect())
        .collect()
}

fn branin(x: &[f64]) -> f64 {
    // Standard Branin on [-5,10]×[0,15], rescaled here from [0,1]².
    let x1 = 15.0f64.mul_add(x[0], -5.0);
    let x2 = 15.0 * x[1];
    let a = 1.0;
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let r = 6.0;
    let s = 10.0;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - r));
    let cosx1 = fs_math::det::cos(x1);
    a * inner * inner + s * (1.0 - t) * cosx1 + s
}

#[test]
fn kernels_are_psd() {
    // Cholesky must succeed (with tiny jitter) on random point sets
    // for every family — the G0 PSD gate.
    let x = rand_points(40, 3, 1);
    for family in [Matern::Half, Matern::ThreeHalves, Matern::FiveHalves] {
        let kernel = Kernel {
            family,
            signal: 1.3,
            lengthscales: vec![0.4, 0.7, 1.1],
        };
        let n = x.len();
        let mut k = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                k[i * n + j] = kernel.eval(&x[i], &x[j]);
            }
            k[i * n + i] += 1e-10;
        }
        assert!(
            fs_la::factor::cholesky(&k, n).is_ok(),
            "{family:?} kernel matrix not PSD"
        );
    }
    log("kernel-psd", "pass", "3 families x 40 random points");
}

#[test]
fn posterior_consistency_union_law() {
    // Fit on D1 ∪ D2 directly vs the same union assembled in a
    // different order: posteriors must agree at probe points (exact
    // inference has no order dependence beyond roundoff).
    let d1 = rand_points(12, 2, 2);
    let d2 = rand_points(9, 2, 3);
    let f = |p: &[f64]| fs_math::det::sin(3.0 * p[0]) + 0.5 * p[1] * p[1];
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.5, 0.5],
    };
    let mut xu: Vec<Vec<f64>> = d1.clone();
    xu.extend(d2.iter().cloned());
    let yu: Vec<f64> = xu.iter().map(|p| f(p)).collect();
    let mut xr: Vec<Vec<f64>> = d2.clone();
    xr.extend(d1.iter().cloned());
    let yr: Vec<f64> = xr.iter().map(|p| f(p)).collect();
    let g1 = Gp::fit(&xu, &yu, kernel.clone(), 1e-8);
    let g2 = Gp::fit(&xr, &yr, kernel, 1e-8);
    let probes = rand_points(15, 2, 4);
    let mut worst = 0.0f64;
    for p in &probes {
        let (m1, v1) = g1.predict(p);
        let (m2, v2) = g2.predict(p);
        worst = worst.max((m1 - m2).abs()).max((v1 - v2).abs());
    }
    assert!(
        worst < 1e-8,
        "union-order consistency violated: {worst:.3e}"
    );
    log(
        "posterior-consistency",
        "pass",
        &format!("worst dev {worst:.1e}"),
    );
}

#[test]
fn known_answer_posteriors() {
    // Noiseless interpolation: mean == y at data, var -> 0.
    let x = rand_points(10, 2, 5);
    let y: Vec<f64> = x.iter().map(|p| p[0] - 2.0 * p[1]).collect();
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.6, 0.6],
    };
    let gp = Gp::fit(&x, &y, kernel, 1e-10);
    for (xi, yi) in x.iter().zip(&y) {
        let (m, v) = gp.predict(xi);
        assert!((m - yi).abs() < 1e-5, "interpolation broke: {m} vs {yi}");
        assert!(v < 1e-5, "variance at datum not ~0: {v:.3e}");
    }
    // One-point closed form: μ(x*) = k(x*,x₁)/(k(x₁,x₁)+σ²)·y₁.
    let x1 = vec![vec![0.3f64, 0.4]];
    let y1 = vec![2.0f64];
    let kernel1 = Kernel {
        family: Matern::ThreeHalves,
        signal: 1.5,
        lengthscales: vec![0.5, 0.8],
    };
    let noise = 0.1;
    let gp1 = Gp::fit(&x1, &y1, kernel1.clone(), noise);
    let xs = vec![0.5f64, 0.1];
    let (m, v) = gp1.predict(&xs);
    let k11 = kernel1.eval(&x1[0], &x1[0]);
    let ks = kernel1.eval(&x1[0], &xs);
    let m_hand = ks / (k11 + noise) * y1[0];
    let v_hand = kernel1.eval(&xs, &xs) - ks * ks / (k11 + noise);
    assert!((m - m_hand).abs() < 1e-12, "{m} vs hand {m_hand}");
    assert!((v - v_hand).abs() < 1e-12, "{v} vs hand {v_hand}");
    log(
        "known-answer",
        "pass",
        "interpolation + 1-point closed form",
    );
}

#[test]
fn kernel_dimension_mismatches_fail_fast() {
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.5, 0.5],
    };
    let eval_mismatch = catch_unwind(AssertUnwindSafe(|| kernel.eval(&[0.1, 0.2], &[0.1])));
    assert!(
        eval_mismatch.is_err(),
        "kernel eval must reject mismatched point dimensions"
    );
    let scale_mismatch = Kernel {
        lengthscales: vec![0.5],
        ..kernel.clone()
    };
    let ard_mismatch = catch_unwind(AssertUnwindSafe(|| {
        scale_mismatch.eval(&[0.1, 0.2], &[0.3, 0.4])
    }));
    assert!(
        ard_mismatch.is_err(),
        "kernel eval must reject ARD lengthscale dimension mismatches"
    );
    let gp = Gp::fit(&[vec![0.2, 0.3]], &[1.0], kernel, 1e-6);
    let predict_mismatch = catch_unwind(AssertUnwindSafe(|| gp.predict(&[0.4])));
    assert!(
        predict_mismatch.is_err(),
        "GP prediction must reject probe dimension mismatches"
    );
    log("kernel-dimensions", "pass", "mismatches fail fast");
}

#[test]
fn normal_functions_accuracy() {
    // Φ at table values; Φ⁻¹ round trip at the documented accuracy.
    assert!((phi_cdf(0.0) - 0.5).abs() < 1e-9);
    assert!((phi_cdf(1.96) - 0.975_002_1).abs() < 1e-6);
    assert!((phi_cdf(-1.0) - 0.158_655_25).abs() < 1e-6);
    let mut worst = 0.0f64;
    for k in 1..100 {
        let p = f64::from(k) / 100.0;
        let rt = phi_cdf(phi_inv(p));
        worst = worst.max((rt - p).abs());
    }
    assert!(worst < 1e-6, "round-trip error {worst:.3e}");
    log("normal-accuracy", "pass", &format!("roundtrip {worst:.1e}"));
}

#[test]
fn ei_laws_and_qei_dominance() {
    let x = rand_points(15, 2, 6);
    let y: Vec<f64> = x.iter().map(|p| branin(p)).collect();
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 20.0,
        lengthscales: vec![0.3, 0.3],
    };
    let gp = Gp::fit(&x, &y, kernel, 1e-6);
    let f_best = y.iter().copied().fold(f64::INFINITY, f64::min);
    // EI ≥ 0 everywhere; EI ~ 0 at a well-observed datum.
    let probes = rand_points(20, 2, 7);
    for p in &probes {
        assert!(expected_improvement(&gp, p, f_best, 0.0) >= 0.0);
    }
    let at_datum = expected_improvement(&gp, &x[3], f_best, 0.0);
    assert!(
        at_datum < 1e-3,
        "EI at noiseless datum should be ~0: {at_datum:.3e}"
    );
    // q-EI(X ∪ {x}) ≥ q-EI(X) (monotone in the batch, fixed bank) and
    // q-EI({x}) ≈ EI(x) at matching sample banks (MC tolerance).
    let bank2 = normal_bank(4096, 2, 99);
    let bank1: Vec<f64> = (0..4096).map(|s| bank2[s * 2]).collect();
    let xa = vec![0.25f64, 0.6];
    let xb = vec![0.7f64, 0.2];
    let q1 = q_expected_improvement(&gp, core::slice::from_ref(&xa), f_best, &bank1);
    let q2 = q_expected_improvement(&gp, &[xa.clone(), xb], f_best, &bank2);
    assert!(
        q2 >= q1 - 1e-9,
        "batch q-EI must dominate its sub-batch: {q2} vs {q1}"
    );
    let ei_closed = expected_improvement(&gp, &xa, f_best, 0.0);
    let rel = (q1 - ei_closed).abs() / ei_closed.max(1e-12);
    assert!(
        rel < 0.05,
        "q-EI(1) vs closed-form EI: {q1:.6e} vs {ei_closed:.6e} (rel {rel:.3})"
    );
    log(
        "ei-qei",
        "pass",
        &format!("qEI(1) vs EI rel {rel:.3}, dominance ok"),
    );
}

#[test]
fn bo_beats_random_on_branin_and_replays() {
    // The acceptance's baseline comparison, ledgered: median best-found
    // after a matched budget over a fixed seed set, EI-BO vs
    // scrambled-Sobol random search. Branin's global minimum is
    // 0.397887.
    let budget_init = 8usize;
    let iters = 10usize;
    let seeds = [11u64, 23, 47];
    let config_for = |seed: u64| BoConfig {
        bounds: (0.0, 1.0),
        family: Matern::FiveHalves,
        log_box: (-2.5, 1.0),
        hyper_starts: 4,
        acq_starts: 3,
        acq_evals: 300,
        q: 1,
        mc_samples: 256,
        seed,
    };
    let mut bo_bests = Vec::new();
    let mut rand_bests = Vec::new();
    for &seed in &seeds {
        let mut f = |x: &[f64]| branin(x);
        let rep = minimize(&mut f, 2, budget_init, iters, &config_for(seed));
        bo_bests.push(*rep.best_trace.last().expect("trace"));
        // Random baseline: same TOTAL budget from the same generator
        // family (scrambled Sobol at a shifted seed).
        let total = budget_init + iters;
        let sobol = fs_rand::qmc::Sobol::scrambled(2, seed ^ 0xDEAD);
        let mut best = f64::INFINITY;
        let mut pt = [0.0f64; 2];
        for s in 0..total {
            sobol.point(u32::try_from(s + 1).expect("small"), &mut pt);
            best = best.min(branin(&pt));
        }
        rand_bests.push(best);
    }
    let med = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let bo_med = med(&mut bo_bests);
    let rand_med = med(&mut rand_bests);
    assert!(
        bo_med < rand_med,
        "BO must beat random at matched budget: {bo_med:.4} vs {rand_med:.4}"
    );
    assert!(
        bo_med < 1.0,
        "BO should approach Branin's optimum (0.3979): {bo_med:.4}"
    );
    // Bitwise replay of a whole BO run.
    let mut f1 = |x: &[f64]| branin(x);
    let r1 = minimize(&mut f1, 2, budget_init, 3, &config_for(7));
    let mut f2 = |x: &[f64]| branin(x);
    let r2 = minimize(&mut f2, 2, budget_init, 3, &config_for(7));
    assert!(
        r1.y.iter()
            .zip(&r2.y)
            .all(|(a, b)| a.to_bits() == b.to_bits()),
        "BO run not bitwise replayable"
    );
    log(
        "bo-vs-random",
        "pass",
        &format!("BO median {bo_med:.4} vs random {rand_med:.4} (opt 0.3979), replay bitwise"),
    );
}

// BUMPED once at tzeh lane-a landing (was 0x4f5a_0601_3cd1_6f46):
// predict_joint's flat 1e-10 jitter became adaptive escalation
// (rel*max_diag from 1e-10 up, diagonal fallback for numerically void
// joints) — REQUIRED for TuRBO's collapsed trust regions, where the
// flat jitter measurably failed on near-duplicate candidate sets.
// Semantic justification per the golden-evidence policy.
const GOLDEN_HASH: u64 = 0x5db8_f433_fd71_f738;

#[test]
fn bo_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // GP posterior fingerprint.
    let x = rand_points(12, 2, 8);
    let y: Vec<f64> = x.iter().map(|p| branin(p)).collect();
    let gp = Gp::fit(
        &x,
        &y,
        Kernel {
            family: Matern::FiveHalves,
            signal: 15.0,
            lengthscales: vec![0.4, 0.5],
        },
        1e-6,
    );
    feed(gp.lml);
    for p in rand_points(6, 2, 9) {
        let (m, v) = gp.predict(&p);
        feed(m);
        feed(v);
    }
    // Acquisition fingerprints.
    let f_best = y.iter().copied().fold(f64::INFINITY, f64::min);
    feed(expected_improvement(&gp, &[0.4, 0.4], f_best, 0.0));
    let bank = normal_bank(512, 3, 42);
    feed(q_expected_improvement(
        &gp,
        &[vec![0.2, 0.3], vec![0.6, 0.7], vec![0.9, 0.1]],
        f_best,
        &bank,
    ));
    // Short BO run fingerprint.
    let mut f = |p: &[f64]| branin(p);
    let rep = minimize(
        &mut f,
        2,
        6,
        2,
        &BoConfig {
            bounds: (0.0, 1.0),
            family: Matern::FiveHalves,
            log_box: (-2.0, 0.5),
            hyper_starts: 2,
            acq_starts: 2,
            acq_evals: 150,
            q: 1,
            mc_samples: 128,
            seed: 5,
        },
    );
    for v in &rep.y {
        feed(*v);
    }
    log("bo-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "bo bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
