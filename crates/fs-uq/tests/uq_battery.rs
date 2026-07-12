//! fs-uq battery (tfz.25 slice 1): KL truncation with captured-
//! variance evidence + covariance reconstruction + sample-statistics
//! gates; PCE known-answer coefficients (Y = exp(a·ξ) has closed-form
//! Hermite coefficients) + surrogate accuracy; QMC-vs-MC advantage
//! MEASURED at matched budget; MLMC telescoping identity AUDITED +
//! variance-per-cost win vs single-level MC; determinism; golden.

use fs_rand::StreamKey;
use fs_uq::{CovarianceKind, KlExpansion, fit_pce, mlmc_estimate};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-uq\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn grid_points(m: usize) -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    for i in 0..m {
        for j in 0..m {
            pts.push([i as f64 / (m - 1) as f64, j as f64 / (m - 1) as f64, 0.0]);
        }
    }
    pts
}

#[test]
fn kl_truncation_evidence_and_reconstruction() {
    let pts = grid_points(7);
    let (sigma2, ell) = (2.0f64, 0.4f64);
    let kl = KlExpansion::build(&pts, CovarianceKind::SquaredExponential, sigma2, ell, 0.99);
    assert!(
        kl.captured_variance >= 0.99,
        "captured variance below target: {}",
        kl.captured_variance
    );
    assert!(
        kl.order() < pts.len() / 2,
        "smooth field should truncate hard: {} of {}",
        kl.order(),
        pts.len()
    );
    let rec =
        kl.covariance_reconstruction_error(&pts, CovarianceKind::SquaredExponential, sigma2, ell);
    assert!(rec < 0.05, "covariance reconstruction too lossy: {rec:.4}");
    // Sample statistics: variance at a probe point ≈ retained-diagonal.
    let germs = kl.qmc_germs(4096, 7);
    let probe = 24usize;
    let mut mean = 0.0f64;
    let mut m2 = 0.0f64;
    for g in &germs {
        let f = kl.realize(g);
        mean += f[probe];
        m2 = f[probe].mul_add(f[probe], m2);
    }
    let n = germs.len() as f64;
    mean /= n;
    let var = m2 / n - mean * mean;
    // Retained variance at the probe = Σ λₖφₖ(probe)².
    let target: f64 = kl
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(k, &lam)| lam * kl.modes[k * kl.n + probe] * kl.modes[k * kl.n + probe])
        .sum();
    let rel = (var - target).abs() / target;
    assert!(rel < 0.05, "sampled variance off: {var:.4} vs {target:.4}");
    log(
        "kl",
        "pass",
        &format!(
            "order {} of {}, captured {:.4}, recon {rec:.3}, var rel {rel:.3}",
            kl.order(),
            pts.len(),
            kl.captured_variance
        ),
    );
}

#[test]
fn pce_known_answer_exponential() {
    // Y = exp(a·ξ): orthonormal-Hermite coefficients are
    // c_k = e^{a²/2}·aᵏ/√(k!) — closed form; mean e^{a²/2},
    // variance e^{a²}(e^{a²} − 1).
    let a = 0.6f64;
    let mut s = StreamKey {
        seed: 71,
        kernel: 0x00C4,
        tile: 1,
    }
    .stream();
    let n = 400usize;
    let xi: Vec<Vec<f64>> = (0..n).map(|_| vec![s.next_normal()]).collect();
    let y: Vec<f64> = xi.iter().map(|x| fs_math::det::exp(a * x[0])).collect();
    let pce = fit_pce(&xi, &y, 8);
    let ea2 = fs_math::det::exp(a * a / 2.0);
    let mut fact = 1.0f64;
    for k in 0..6 {
        if k > 1 {
            fact *= k as f64;
        }
        let expected = ea2 * fs_math::det::pow(a, k as f64) / fs_math::det::sqrt(fact);
        let got = pce.coefficients[k];
        assert!(
            (got - expected).abs() < 0.02 * expected.abs().max(0.05),
            "c_{k}: {got:.5} vs closed form {expected:.5}"
        );
    }
    let mean_true = ea2;
    let var_true = fs_math::det::exp(a * a) * (fs_math::det::exp(a * a) - 1.0);
    assert!((pce.mean() - mean_true).abs() < 0.01);
    assert!((pce.variance() - var_true).abs() / var_true < 0.05);
    // Surrogate accuracy at fresh points.
    let mut worst = 0.0f64;
    for _ in 0..50 {
        let x = s.next_normal().clamp(-2.5, 2.5);
        let err = (pce.eval(&[x]) - fs_math::det::exp(a * x)).abs();
        worst = worst.max(err);
    }
    assert!(worst < 0.02, "surrogate error {worst:.4}");
    log(
        "pce",
        "pass",
        &format!(
            "coeffs match closed form, mean {:.4}/{mean_true:.4}, surrogate max err {worst:.4}",
            pce.mean()
        ),
    );
}

#[test]
fn qmc_beats_mc_measured() {
    // Smooth integrand over [0,1]^5 with known mean: product of
    // (1 + 0.4·(uᵢ − ½)) has mean exactly 1. RMSE over 20 replicates
    // at n = 2048: scrambled Sobol must beat plain MC by ≥ 3×.
    let d = 5usize;
    let n = 2048usize;
    let reps = 20usize;
    let f = |u: &[f64]| -> f64 { u.iter().map(|ui| 0.4f64.mul_add(ui - 0.5, 1.0)).product() };
    let mut qmc_se = 0.0f64;
    let mut mc_se = 0.0f64;
    for r in 0..reps {
        let sobol = fs_rand::qmc::Sobol::scrambled(d, 1000 + r as u64);
        let mut pt = vec![0.0f64; d];
        let mut acc = 0.0f64;
        for s in 0..n {
            sobol.point(u32::try_from(s + 1).expect("small"), &mut pt);
            acc += f(&pt);
        }
        let e = acc / n as f64 - 1.0;
        qmc_se = e.mul_add(e, qmc_se);
        let mut st = StreamKey {
            seed: 81,
            kernel: 0x00AC,
            tile: r as u32,
        }
        .stream();
        let mut acc2 = 0.0f64;
        for _ in 0..n {
            let u: Vec<f64> = (0..d).map(|_| st.next_f64()).collect();
            acc2 += f(&u);
        }
        let e2 = acc2 / n as f64 - 1.0;
        mc_se = e2.mul_add(e2, mc_se);
    }
    let qmc_rmse = fs_math::det::sqrt(qmc_se / reps as f64);
    let mc_rmse = fs_math::det::sqrt(mc_se / reps as f64);
    let ratio = mc_rmse / qmc_rmse;
    assert!(
        ratio > 3.0,
        "QMC advantage below 3x at n={n}: MC {mc_rmse:.2e} vs QMC {qmc_rmse:.2e}"
    );
    log(
        "qmc-vs-mc",
        "pass",
        &format!("MC {mc_rmse:.2e} vs QMC {qmc_rmse:.2e} = {ratio:.1}x at n={n}, d={d}"),
    );
}

#[test]
fn mlmc_telescoping_and_cost_win() {
    // Synthetic solver ladder: level ℓ "discretizes" E[G(ξ)] with bias
    // b·2^{-2ℓ} and the COUPLED corrections have variance decaying
    // 2^{-2ℓ} while cost grows 2^ℓ (the classic MLMC regime).
    let nl = 5usize;
    let costs: Vec<f64> = (0..nl).map(|l| fs_math::det::pow(2.0, l as f64)).collect();
    let g = |xi: f64| xi * xi; // E[G] = 1 for ξ ~ N(0,1)
    let level_value = |l: usize, xi: f64| -> f64 {
        // P_ℓ = G(ξ)·(1 + 2^{-2ℓ}·0.3) + fine-scale noise handled by
        // coupling: the correction Y_ℓ uses the SAME germ.
        let bias = 0.3 * fs_math::det::pow(2.0, -2.0 * l as f64);
        g(xi) * (1.0 + bias)
    };
    let mut sampler = |l: usize, idx: u64| -> f64 {
        let mut st = StreamKey {
            seed: 91,
            kernel: 0x300C + l as u32,
            tile: u32::try_from(idx & 0xFFFF_FFFF).expect("fits"),
        }
        .stream();
        let xi = st.next_normal();
        if l == 0 {
            level_value(0, xi)
        } else {
            level_value(l, xi) - level_value(l - 1, xi)
        }
    };
    let rep = mlmc_estimate(&mut sampler, &costs, 200, 5e-5);
    // Telescoping AUDIT: Σ level means must equal the fine-level mean
    // estimated from the SAME per-level germ sequences (identity of
    // the estimator, checked by reconstruction).
    let sum_means: f64 = rep.levels.iter().map(|l| l.mean).sum();
    assert!(
        (rep.estimate - sum_means).abs() < 1e-14,
        "telescoping bookkeeping broken"
    );
    // The true target is E[G]·(1+bias_L) = 1·(1 + 0.3·2^{-8}).
    let truth = 1.0 + 0.3 * fs_math::det::pow(2.0, -8.0);
    assert!(
        (rep.estimate - truth).abs() < 0.05,
        "MLMC estimate off: {} vs {truth}",
        rep.estimate
    );
    // Variance decay ⇒ allocation puts most samples on coarse levels.
    assert!(
        rep.levels[0].samples > 10 * rep.levels[nl - 1].samples,
        "allocation should favor coarse: {:?}",
        rep.levels.iter().map(|l| l.samples).collect::<Vec<_>>()
    );
    // Cost win vs single-level MC at the SAME estimator variance:
    // single-level needs N = V_L^single/target on the finest level.
    let mut sum = 0.0f64;
    let mut sum2 = 0.0f64;
    let npilot = 2000usize;
    let mut st2 = StreamKey {
        seed: 92,
        kernel: 0x51CE,
        tile: 0,
    }
    .stream();
    for _ in 0..npilot {
        let v = level_value(nl - 1, st2.next_normal());
        sum += v;
        sum2 = v.mul_add(v, sum2);
    }
    let var_fine = sum2 / npilot as f64 - (sum / npilot as f64).powi(2);
    let single_cost = (var_fine / rep.estimator_variance) * costs[nl - 1];
    let win = single_cost / rep.total_cost;
    assert!(
        win > 3.0,
        "MLMC cost win below 3x: single {single_cost:.1} vs multilevel {:.1}",
        rep.total_cost
    );
    log(
        "mlmc",
        "pass",
        &format!(
            "estimate {:.4} (truth {truth:.4}), cost win {win:.1}x, samples {:?}",
            rep.estimate,
            rep.levels.iter().map(|l| l.samples).collect::<Vec<_>>()
        ),
    );
}

const GOLDEN_HASH: u64 = 0x0ed2_4974_dc37_bbc6; // recorded at tfz.25 slice 1, frozen

#[test]
fn uq_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let pts = grid_points(5);
    let kl = KlExpansion::build(&pts, CovarianceKind::Exponential, 1.5, 0.5, 0.95);
    feed(kl.captured_variance);
    for v in kl.eigenvalues.iter().take(6) {
        feed(*v);
    }
    let germs = kl.qmc_germs(4, 11);
    for g in &germs {
        let field = kl.realize(g);
        for v in field.iter().step_by(7) {
            feed(*v);
        }
    }
    // PCE fingerprint.
    let mut s = StreamKey {
        seed: 71,
        kernel: 0x00C4,
        tile: 9,
    }
    .stream();
    let xi: Vec<Vec<f64>> = (0..150)
        .map(|_| vec![s.next_normal(), s.next_normal()])
        .collect();
    let y: Vec<f64> = xi
        .iter()
        .map(|x| fs_math::det::exp(0.3 * x[0]) + 0.5 * x[1])
        .collect();
    let pce = fit_pce(&xi, &y, 4);
    for v in pce.coefficients.iter().take(8) {
        feed(*v);
    }
    log("uq-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "uq bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
#[should_panic(expected = "nonzero pilot")]
fn mlmc_rejects_a_zero_pilot() {
    // Regression: pilot = 0 divides by n = 0 in the variance estimate; the
    // 1e-30 floor then masks the NaN and the report would claim a tiny,
    // fake-confident estimator variance built on ZERO data. Fail closed.
    let mut sampler = |_l: usize, _g: u64| 1.0;
    let _ = mlmc_estimate(&mut sampler, &[1.0, 2.0], 0, 5e-5);
}

#[test]
#[should_panic(expected = "positive target variance")]
fn mlmc_rejects_a_nonpositive_target_variance() {
    // Regression: target_variance = 0 makes n_opt = INFINITY as usize =
    // usize::MAX, so the top-up loop samples essentially forever (a hang).
    let mut sampler = |_l: usize, _g: u64| 1.0;
    let _ = mlmc_estimate(&mut sampler, &[1.0, 2.0], 8, 0.0);
}

#[test]
#[should_panic(expected = "positive per-level costs")]
fn mlmc_rejects_a_zero_cost_level() {
    // Regression: costs[l] = 0 makes v / costs[l] = +inf → n_opt = usize::MAX
    // → the same unbounded sampling loop. (NaN costs fail the same > 0 test.)
    let mut sampler = |_l: usize, _g: u64| 1.0;
    let _ = mlmc_estimate(&mut sampler, &[1.0, 0.0], 8, 5e-5);
}
