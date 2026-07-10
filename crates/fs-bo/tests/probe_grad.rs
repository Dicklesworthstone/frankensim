//! Regression probe (a2g2 lane a): the Σ-diagonal Matérn eval hits
//! √0, whose BACKWARD derivative is infinite — the first draft
//! produced NaN gradients that two FD gates missed by passing
//! VACUOUSLY on flat regions. This probe pins an ACTIVE surface and
//! asserts finite, FD-matching gradients forever.
#![cfg(feature = "tape-acq")]
use fs_bo::acq_grad::qei_gradient;
use fs_bo::{Gp, Kernel, Matern, normal_bank, q_expected_improvement};
use fs_rand::StreamKey;

#[test]
fn probe_gradient_values() {
    fn branin(x: &[f64]) -> f64 {
        let x1 = 15.0f64.mul_add(x[0], -5.0);
        let x2 = 15.0 * x[1];
        let b = 5.1 / (4.0 * core::f64::consts::PI * core::f64::consts::PI);
        let c = 5.0 / core::f64::consts::PI;
        let t = 1.0 / (8.0 * core::f64::consts::PI);
        let inner = b.mul_add(-(x1 * x1), c.mul_add(x1, x2 - 6.0));
        inner * inner + 10.0 * (1.0 - t) * fs_math::det::cos(x1) + 10.0
    }
    let mut s = StreamKey {
        seed: 171,
        kernel: 0x0A69,
        tile: 1,
    }
    .stream();
    let x: Vec<Vec<f64>> = (0..25)
        .map(|_| (0..2).map(|_| s.next_f64()).collect())
        .collect();
    let y: Vec<f64> = x.iter().map(|p| branin(p) / 50.0).collect();
    let f_best = y.iter().copied().fold(f64::INFINITY, f64::min);
    let kernel = Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.3, 0.3],
    };
    let gp = Gp::fit(&x, &y, kernel, 1e-4);
    let bank = normal_bank(256, 3, 55);
    let cands: Vec<Vec<f64>> = (0..3)
        .map(|_| (0..2).map(|_| s.next_f64()).collect())
        .collect();
    println!("cands = {cands:?}");
    let (val, grad) = qei_gradient(&gp, &cands, f_best, &bank);
    assert!(
        val > 1e-3,
        "probe surface must be ACTIVE (vacuous-pass guard): {val:.3e}"
    );
    assert!(
        grad.iter().all(|g| g.is_finite()),
        "NaN regression: {grad:?}"
    );
    let gnorm = grad.iter().map(|g| g.abs()).fold(0.0f64, f64::max);
    assert!(gnorm > 0.1, "gradient dead on the active probe: {grad:?}");
    let flat: Vec<f64> = cands.iter().flatten().copied().collect();
    for (i, &gi) in grad.iter().enumerate() {
        let mut fp = flat.clone();
        fp[i] += 1e-6;
        let mut fm = flat.clone();
        fm[i] -= 1e-6;
        let un = |f: &[f64]| vec![f[0..2].to_vec(), f[2..4].to_vec(), f[4..6].to_vec()];
        let vp = q_expected_improvement(&gp, &un(&fp), f_best, &bank);
        let vm = q_expected_improvement(&gp, &un(&fm), f_best, &bank);
        let fd = (vp - vm) / 2e-6;
        assert!(
            (fd - gi).abs() / gnorm < 1e-5,
            "tape/FD mismatch at {i}: {gi:.8e} vs {fd:.8e}"
        );
    }
}
