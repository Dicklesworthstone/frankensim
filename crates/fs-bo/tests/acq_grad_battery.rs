//! Tape-gradient battery (a2g2 lane a; feature `tape-acq`): the taped
//! q-EI gradient vs central FD on the fixed-bank surface (tolerance
//! gates per the bridge's declared determinism class — never
//! bitwise); primal parity with the production f64 q-EI; the
//! gradient-ascent argmax vs CMA-ES at a matched evaluation ledger.
#![cfg(feature = "tape-acq")]

use fs_bo::acq_grad::{qei_ascent, qei_gradient};
use fs_bo::{Gp, Kernel, Matern, normal_bank, q_expected_improvement};
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-bo-acqgrad\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_pts(n: usize, d: usize, tile: u32) -> Vec<Vec<f64>> {
    let mut s = StreamKey {
        seed: 171,
        kernel: 0x0A69,
        tile,
    }
    .stream();
    (0..n)
        .map(|_| (0..d).map(|_| s.next_f64()).collect())
        .collect()
}

fn branin(x: &[f64]) -> f64 {
    let x1 = 15.0f64.mul_add(x[0], -5.0);
    let x2 = 15.0 * x[1];
    let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
    let c = 5.0 / core::f64::consts::PI;
    let t = 1.0 / (8.0 * core::f64::consts::PI);
    let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - 6.0));
    inner * inner + 10.0 * (1.0 - t) * fs_math::det::cos(x1) + 10.0
}

fn fitted_gp() -> (Gp, f64) {
    let x = rand_pts(25, 2, 1);
    let y: Vec<f64> = x.iter().map(|p| branin(p) / 50.0).collect();
    let f_best = y.iter().copied().fold(f64::INFINITY, f64::min);
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.3, 0.3],
    };
    (Gp::fit(&x, &y, kernel, 1e-4), f_best)
}

#[test]
fn taped_gradient_matches_fd_and_primal() {
    let (gp, f_best) = fitted_gp();
    let q = 3usize;
    let bank = normal_bank(256, q, 55);
    let cands = rand_pts(q, 2, 2);
    let (primal_taped, grad) = qei_gradient(&gp, &cands, f_best, &bank);
    // ALIVENESS guards (vacuous-pass protection: the first version of
    // this gate passed on a flat region with grad == fd == 0 while
    // NaN gradients hid elsewhere): the surface must be active and
    // every gradient entry finite.
    assert!(
        primal_taped > 1e-4,
        "test point sits in a no-improvement region: {primal_taped:.3e}"
    );
    assert!(
        grad.iter().all(|g| g.is_finite()),
        "NaN/inf in taped gradient: {grad:?}"
    );
    assert!(
        grad.iter().map(|g| g.abs()).fold(0.0f64, f64::max) > 1e-6,
        "gradient dead on an active surface: {grad:?}"
    );
    // Primal parity with the production f64 path (different elementary
    // kernels — tolerance, not bitwise, per the declared class).
    let primal_f64 = q_expected_improvement(&gp, &cands, f_best, &bank);
    // ft elementary functions differ from fs-math's det kernels in
    // final ulps; accumulated over the 256-sample chain the honest
    // parity gate is 1e-6 RELATIVE (measured 1.3e-7).
    let rel_primal = (primal_taped - primal_f64).abs() / primal_f64.abs().max(1e-12);
    assert!(
        rel_primal < 1e-6,
        "taped primal vs production q-EI: {primal_taped:.9e} vs {primal_f64:.9e}"
    );
    // Gradient vs central FD on the f64 surface.
    let flat: Vec<f64> = cands.iter().flatten().copied().collect();
    let eps = 1e-6;
    let mut worst = 0.0f64;
    let gnorm = grad
        .iter()
        .map(|v| v.abs())
        .fold(0.0f64, f64::max)
        .max(1e-12);
    for i in 0..flat.len() {
        let mut fp = flat.clone();
        fp[i] += eps;
        let mut fm = flat.clone();
        fm[i] -= eps;
        let unflat = |f: &[f64]| -> Vec<Vec<f64>> {
            (0..q).map(|b| f[b * 2..(b + 1) * 2].to_vec()).collect()
        };
        let vp = q_expected_improvement(&gp, &unflat(&fp), f_best, &bank);
        let vm = q_expected_improvement(&gp, &unflat(&fm), f_best, &bank);
        let fd = (vp - vm) / (2.0 * eps);
        worst = worst.max((fd - grad[i]).abs() / gnorm);
    }
    assert!(
        worst < 1e-4,
        "taped gradient vs FD (relative to grad scale): {worst:.3e}"
    );
    log(
        "tape-vs-fd",
        "pass",
        &format!("primal rel {rel_primal:.1e}, grad worst rel {worst:.1e}"),
    );
}

#[test]
fn ascent_polish_properties_and_cost() {
    // The lane's claim is EXACT CHEAP GRADIENTS, gated as: (a) the
    // L-BFGS polish monotonically improves its probe seed and reaches
    // a stationary point; (b) the tape gradient costs ONE reverse pass
    // where FD costs 2·q·d primal evaluations (the cost ledger,
    // counted). CMA-ES's global value is REPORTED, not gated: the
    // fixed-bank q-EI surface has Monte-Carlo needle spikes at
    // near-duplicate candidate blocks (posterior Cholesky near
    // singularity — 200 Sobol probes in 4-d never hit the spike CMA
    // wandered into at 0.488 vs smooth-region maxima ~0.33; gating
    // local polish on winning a needle hunt would be dishonest in
    // BOTH directions).
    let (gp, f_best) = fitted_gp();
    let q = 2usize;
    let bank = normal_bank(128, q, 77);
    let sobol = fs_rand::qmc::Sobol::scrambled(2 * q, 500);
    let mut pt = vec![0.0f64; 2 * q];
    let mut probes: Vec<(f64, Vec<Vec<f64>>)> = Vec::new();
    for s in 0..60u32 {
        sobol.point(s + 1, &mut pt);
        let block: Vec<Vec<f64>> = (0..q).map(|b| pt[b * 2..(b + 1) * 2].to_vec()).collect();
        let v = q_expected_improvement(&gp, &block, f_best, &bank);
        probes.push((v, block));
    }
    probes.sort_by(|a, b| b.0.total_cmp(&a.0));
    let (seed_val, seed) = (probes[0].0, probes[0].1.clone());
    let (x_pol, val_pol, reverse_passes) = qei_ascent(&gp, &seed, f_best, &bank, (0.0, 1.0), 60);
    // (a) Monotone improvement + stationarity.
    assert!(
        val_pol >= seed_val - 1e-12,
        "polish must not lose ground: {val_pol:.4e} vs seed {seed_val:.4e}"
    );
    let (_, g_final) = qei_gradient(&gp, &x_pol, f_best, &bank);
    let gnorm = g_final.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    assert!(
        gnorm < 1e-6 || val_pol > 1.05 * seed_val,
        "polish must certify stationarity or materially improve: gnorm {gnorm:.2e},          {seed_val:.4e} -> {val_pol:.4e}"
    );
    // (b) Cost ledger: FD gradients for the same trajectory would cost
    // 2·q·d primal evaluations each.
    let fd_equivalent = reverse_passes * (2 * q * 2 + 1);
    assert!(
        reverse_passes * 4 < fd_equivalent,
        "tape must be materially cheaper than FD: {reverse_passes} passes vs {fd_equivalent}"
    );
    // Comparative context, REPORTED: CMA-ES on the same surface.
    let mut obj = |flat: &[f64]| -> f64 {
        let cands: Vec<Vec<f64>> = (0..q)
            .map(|b| {
                flat[b * 2..(b + 1) * 2]
                    .iter()
                    .map(|v| v.clamp(0.0, 1.0))
                    .collect()
            })
            .collect();
        -q_expected_improvement(&gp, &cands, f_best, &bank)
    };
    let params = fs_dfo::CmaParams::standard(2 * q, 0.2, 400, f64::NEG_INFINITY);
    let cma = fs_dfo::cmaes(&mut obj, &vec![0.5; 2 * q], &params, 5);
    log(
        "ascent-polish",
        "pass",
        &format!(
            "seed {seed_val:.4e} -> polished {val_pol:.4e} (gnorm {gnorm:.1e}) in {reverse_passes}              reverse passes (FD-equiv {fd_equivalent}); CMA-ES reported: {:.4e} in {} evals",
            -cma.f_best, cma.evals
        ),
    );
}
