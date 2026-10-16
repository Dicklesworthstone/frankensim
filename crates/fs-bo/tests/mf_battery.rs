//! Multi-fidelity battery (tzeh lane b): correlation recovery +
//! variance reduction from cheap data; cost-aware ALLOCATION (low
//! fidelity dominates at 10:1 cost ratio); COST-TO-TARGET vs
//! single-fidelity EI-BO at matched total cost (median over fixed
//! seeds, ledgered); bitwise replay; golden.

use fs_bo::{BoConfig, Matern, MfConfig, MfGp, MfKernel, fit_mf, mf_minimize, minimize};
use fs_rand::StreamKey;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-bo-mf\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
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

/// Low fidelity: Branin plus a smooth known bias (the classic
/// two-fidelity test family) — HIGHLY correlated with the target.
fn branin_low(x: &[f64]) -> f64 {
    branin(x) + 3.0f64.mul_add(x[0], -2.0 * x[1]) + 5.0
}

fn rand_pts(n: usize, d: usize, tile: u32) -> Vec<Vec<f64>> {
    let mut s = StreamKey {
        seed: 121,
        kernel: 0x00F1,
        tile,
    }
    .stream();
    (0..n)
        .map(|_| (0..d).map(|_| s.next_f64()).collect())
        .collect()
}

#[test]
fn correlation_recovery_and_variance_reduction() {
    // Fit the joint GP on 12 high + 30 low observations; the learned
    // between-fidelity correlation must be high (the fixture's
    // empirical correlation is ~0.99), and adding the low data must
    // REDUCE high-fidelity posterior variance at held-out points.
    let xh = rand_pts(12, 2, 1);
    let xl = rand_pts(30, 2, 2);
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let mut fid = Vec::new();
    let mut ys = Vec::new();
    for x in &xh {
        ys.push(branin(x));
        xs.push(x.clone());
        fid.push(1);
    }
    for x in &xl {
        ys.push(branin_low(x));
        xs.push(x.clone());
        fid.push(0);
    }
    let n = ys.len() as f64;
    let mean = ys.iter().sum::<f64>() / n;
    let scale = fs_math::det::sqrt(ys.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n);
    let ys_std: Vec<f64> = ys.iter().map(|v| (v - mean) / scale).collect();
    let gp = fit_mf(&xs, &fid, &ys_std, (-2.0, 1.0), 3, 7);
    let rho = gp.kernel.correlation();
    assert!(
        rho > 0.8,
        "learned inter-fidelity correlation too low: {rho:.3}"
    );
    // High-only GP for the variance comparison.
    let ysh: Vec<f64> = ys_std[..12].to_vec();
    let fidh = vec![1usize; 12];
    let gp_high = fit_mf(&xh, &fidh, &ysh, (-2.0, 1.0), 3, 7);
    let probes = rand_pts(10, 2, 3);
    let mut reduced = 0usize;
    for p in &probes {
        let (_, v_joint) = gp.predict(p, 1);
        let (_, v_high) = gp_high.predict(p, 1);
        if v_joint < v_high {
            reduced += 1;
        }
    }
    assert!(
        reduced >= 7,
        "low-fidelity data must reduce high-fidelity variance broadly: {reduced}/10"
    );
    log(
        "correlation",
        "pass",
        &format!("rho {rho:.3}, variance reduced at {reduced}/10 probes"),
    );
}

#[test]
fn mf_kernel_dimension_mismatches_fail_fast() {
    let kernel = MfKernel {
        lengthscales: vec![0.4, 0.6],
        l_fid: [1.2, 0.9, 0.3],
        noise: 1e-6,
    };
    let point_mismatch = catch_unwind(AssertUnwindSafe(|| kernel.eval(&[0.1, 0.2], 0, &[0.3], 1)));
    assert!(
        point_mismatch.is_err(),
        "MF kernel eval must reject mismatched point dimensions"
    );
    let scale_mismatch = MfKernel {
        lengthscales: vec![0.4],
        ..kernel.clone()
    };
    let ard_mismatch = catch_unwind(AssertUnwindSafe(|| {
        scale_mismatch.eval(&[0.1, 0.2], 0, &[0.3, 0.4], 1)
    }));
    assert!(
        ard_mismatch.is_err(),
        "MF kernel eval must reject ARD lengthscale dimension mismatches"
    );
    let fidelity_mismatch = catch_unwind(AssertUnwindSafe(|| {
        kernel.eval(&[0.1, 0.2], 2, &[0.3, 0.4], 1)
    }));
    assert!(
        fidelity_mismatch.is_err(),
        "MF kernel eval must reject invalid fidelity indices"
    );
    log("mf-kernel-dimensions", "pass", "mismatches fail fast");
}

fn mf_config(seed: u64, budget: f64) -> MfConfig {
    MfConfig {
        bounds: (0.0, 1.0),
        log_box: (-2.0, 1.0),
        hyper_starts: 2,
        n_init_low: 10,
        n_init_high: 4,
        cost_low: 1.0,
        cost_high: 10.0,
        budget,
        acq_starts: 2,
        acq_evals: 200,
        seed,
    }
}

#[test]
fn allocation_cost_to_target_and_replay() {
    let seeds = [13u64, 31];
    let budget = 120.0f64;
    let mut mf_bests = Vec::new();
    let mut sf_bests = Vec::new();
    let mut low_share = Vec::new();
    for &seed in &seeds {
        let mut f = |x: &[f64], m: usize| -> f64 { if m == 0 { branin_low(x) } else { branin(x) } };
        let rep = mf_minimize(&mut f, 2, &mf_config(seed, budget));
        // ALLOCATION: at 10:1 costs with a highly-correlated cheap
        // fidelity, low evaluations must dominate the COUNT.
        assert!(
            rep.evals_low > rep.evals_high,
            "cheap fidelity should dominate the count: {} low vs {} high",
            rep.evals_low,
            rep.evals_high
        );
        low_share.push(rep.evals_low as f64 / (rep.evals_low + rep.evals_high) as f64);
        mf_bests.push(rep.f_best_high);
        // Single-fidelity EI-BO at the SAME total cost (budget/10
        // high evaluations).
        let sf_evals = (budget / 10.0) as usize;
        let mut fh = |x: &[f64]| branin(x);
        let sf = minimize(
            &mut fh,
            2,
            4,
            sf_evals - 4,
            &BoConfig {
                bounds: (0.0, 1.0),
                family: Matern::FiveHalves,
                log_box: (-2.5, 1.0),
                hyper_starts: 2,
                acq_starts: 2,
                acq_evals: 200,
                q: 1,
                mc_samples: 128,
                seed,
            },
        );
        sf_bests.push(*sf.best_trace.last().expect("trace"));
    }
    let med = |v: &mut Vec<f64>| -> f64 {
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let mf_med = med(&mut mf_bests);
    let sf_med = med(&mut sf_bests);
    // COST-TO-TARGET: at matched cost the MF run must be at least
    // competitive (<= 1.5x the single-fidelity best-found gap to the
    // 0.397887 optimum) — and the ledger records both numbers.
    let opt = 0.397_887f64;
    assert!(
        mf_med - opt <= 1.5 * (sf_med - opt) + 0.05,
        "MF-BO fell behind single-fidelity at matched cost: {mf_med:.4} vs {sf_med:.4}"
    );
    // Bitwise replay.
    let mut f1 = |x: &[f64], m: usize| if m == 0 { branin_low(x) } else { branin(x) };
    let r1 = mf_minimize(&mut f1, 2, &mf_config(5, 60.0));
    let mut f2 = |x: &[f64], m: usize| if m == 0 { branin_low(x) } else { branin(x) };
    let r2 = mf_minimize(&mut f2, 2, &mf_config(5, 60.0));
    assert_eq!(r1.trace.len(), r2.trace.len());
    for (a, b) in r1.trace.iter().zip(&r2.trace) {
        assert!(a.1.to_bits() == b.1.to_bits(), "MF run not replayable");
    }
    log(
        "allocation-cost",
        "pass",
        &format!(
            "MF {mf_med:.4} vs SF {sf_med:.4} (opt 0.3979) at cost {budget}, low share {:.0}%",
            100.0 * low_share.iter().sum::<f64>() / low_share.len() as f64
        ),
    );
}

const GOLDEN_HASH: u64 = 0x6411_f077_1d5e_9f88; // recorded at tzeh lane b, frozen

#[test]
fn mf_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // Kernel + fit fingerprint.
    let kernel = MfKernel {
        lengthscales: vec![0.4, 0.6],
        l_fid: [1.2, 0.9, 0.3],
        noise: 1e-6,
    };
    feed(kernel.correlation());
    feed(kernel.eval(&[0.2, 0.3], 0, &[0.5, 0.1], 1));
    let xs = rand_pts(10, 2, 8);
    let fid: Vec<usize> = (0..10).map(|i| i % 2).collect();
    let ys: Vec<f64> = xs
        .iter()
        .zip(&fid)
        .map(|(x, &m)| if m == 0 { branin_low(x) } else { branin(x) } / 50.0)
        .collect();
    let gp = MfGp::try_fit(&xs, &fid, &ys, kernel).expect("SPD");
    feed(gp.lml);
    let (m0, v0) = gp.predict(&[0.4, 0.4], 1);
    feed(m0);
    feed(v0);
    feed(gp.posterior_fid_correlation(&[0.4, 0.4]));
    // Short MF-BO fingerprint.
    let mut f = |x: &[f64], m: usize| if m == 0 { branin_low(x) } else { branin(x) };
    let rep = mf_minimize(&mut f, 2, &mf_config(9, 40.0));
    feed(rep.f_best_high);
    feed(rep.cost);
    feed(rep.learned_correlation);
    log("mf-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "mf bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
