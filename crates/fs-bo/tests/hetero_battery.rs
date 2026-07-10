//! Heteroscedastic-GP + anytime-stopped noisy BO battery (a2g2 lane
//! b): the COINCIDENT-POINT closed form (2×2 posterior algebra by
//! hand — precision-weighted pull toward the low-noise observation);
//! declared-noisy clusters do not drag predictions (measured against
//! the same data fit homoscedastically); an anytime-stopped noisy BO
//! loop composing fs-uq's e-process primitive — adaptive replication
//! measured against a fixed-replication baseline; golden.

use fs_bo::{Gp, Kernel, Matern};
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!("{{\"suite\":\"fs-bo-hetero\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}");
}

fn kernel() -> Kernel {
    Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.4],
    }
}

#[test]
fn coincident_point_closed_form() {
    // Two observations at the SAME x with noises σa² ≪ σb²: the
    // posterior mean at x is k₀·1ᵀ(K+S)⁻¹y with K = k₀·11ᵀ — a 2×2
    // system solvable by hand. The precision-weighted pull must favor
    // the low-noise observation.
    let x = vec![vec![0.5f64], vec![0.5f64]];
    let (ya, yb) = (1.0f64, -1.0f64);
    let (sa, sb) = (1e-4f64, 1e-1f64);
    let gp = Gp::try_fit_diag(&x, &[ya, yb], kernel(), &[sa, sb]).expect("SPD");
    let (mean, _) = gp.predict(&[0.5]);
    // Hand algebra: K+S = [[k0+sa, k0], [k0, k0+sb]], k0 = signal = 1.
    let k0 = 1.0f64;
    let det = (k0 + sa).mul_add(k0 + sb, -(k0 * k0));
    let inv = [
        (k0 + sb) / det,
        -k0 / det,
        -k0 / det,
        (k0 + sa) / det,
    ];
    let alpha = [
        inv[0].mul_add(ya, inv[1] * yb),
        inv[2].mul_add(ya, inv[3] * yb),
    ];
    let mean_hand = k0 * (alpha[0] + alpha[1]);
    assert!(
        (mean - mean_hand).abs() < 1e-10,
        "posterior mean {mean:.8} vs hand {mean_hand:.8}"
    );
    // Precision weighting: the mean must sit near ya (low noise).
    assert!(
        (mean - ya).abs() < 0.01,
        "low-noise observation must dominate: mean {mean:.4} vs ya {ya}"
    );
    log(
        "coincident",
        "pass",
        &format!("mean {mean:.6} == hand {mean_hand:.6}, pulled to low-noise"),
    );
}

#[test]
fn declared_noise_does_not_drag() {
    // Clean sine data plus a corrupted cluster DECLARED noisy: the
    // heteroscedastic fit must track the truth where the homoscedastic
    // fit (same data, uniform noise) gets dragged.
    let mut s = StreamKey { seed: 151, kernel: 0x11E7, tile: 0 }.stream();
    let truth = |x: f64| fs_math::det::sin(3.0 * x);
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let mut ys = Vec::new();
    let mut noises = Vec::new();
    for k in 0..20 {
        let x = f64::from(k) / 19.0;
        xs.push(vec![x]);
        ys.push(truth(x) + 0.01 * s.next_normal());
        noises.push(1e-4);
    }
    // Corrupted cluster near x = 0.5, declared with honest big noise.
    for _ in 0..6 {
        let x = 0.02f64.mul_add(s.next_normal(), 0.5);
        xs.push(vec![x]);
        ys.push(truth(x) + 2.0 + 0.5 * s.next_normal());
        noises.push(4.0);
    }
    let hetero = Gp::try_fit_diag(&xs, &ys, kernel(), &noises).expect("SPD");
    let homo = Gp::fit(&xs, &ys, kernel(), 1e-2);
    let probe = [0.5f64];
    let (mh, _) = hetero.predict(&probe);
    let (mo, _) = homo.predict(&probe);
    let err_hetero = (mh - truth(0.5)).abs();
    let err_homo = (mo - truth(0.5)).abs();
    assert!(
        err_hetero < 0.1,
        "heteroscedastic fit dragged by declared noise: err {err_hetero:.3}"
    );
    assert!(
        err_homo > 3.0 * err_hetero,
        "the comparison should show the drag: hetero {err_hetero:.3} vs homo {err_homo:.3}"
    );
    log(
        "no-drag",
        "pass",
        &format!("hetero err {err_hetero:.4} vs homo err {err_homo:.4} at the corrupted cluster"),
    );
}

#[test]
fn anytime_stopped_noisy_bo() {
    // E-RACING candidate elimination: sample each challenger until its
    // confidence sequence SEPARATES from the incumbent's (eliminate or
    // dethrone) or a half-width floor is reached. Clearly-bad
    // candidates separate FAST — adaptive allocation, measured against
    // uniform max replication. (First draft compared half-width-target
    // stopping across candidates — a fixed-σ CS half-width depends
    // only on n, so every candidate stopped at the SAME count and the
    // comparison was vacuous; racing is the real e-process pattern.)
    // Plus the Bet-5 validity claim: every stopped interval covers the
    // true mean AT the stopping time.
    let truth = |x: f64| 0.5f64.mul_add(-fs_math::det::exp(-30.0 * (x - 0.7) * (x - 0.7)), 0.6);
    let candidates: Vec<f64> = (0..12).map(|k| f64::from(k) / 11.0).collect();
    let noise = 0.15f64;
    let cap = 4000u64;
    let mut incumbent: Option<(usize, f64, f64)> = None; // (idx, lo, hi)
    let mut total_adaptive = 0u64;
    let mut per_candidate = Vec::new();
    let mut worst_miss = false;
    for (ci, &c) in candidates.iter().enumerate() {
        let mut s = StreamKey {
            seed: 161,
            kernel: 0x0A27,
            tile: u32::try_from(ci).expect("few"),
        }
        .stream();
        let mut cs = fs_eproc::GaussianMixtureCs::new(0.5, 1.0, 0.05);
        let mut n = 0u64;
        let (mut lo, mut hi) = (0.0f64, 1.0f64);
        while n < cap {
            let x = (truth(c) + noise * s.next_normal()).clamp(0.0, 1.0);
            cs.observe(x);
            n += 1;
            let (center, radius) = cs.interval().expect("data seen");
            lo = center - radius;
            hi = center + radius;
            if let Some((_, _, inc_hi)) = incumbent {
                // Challenger provably worse than the incumbent: eliminate.
                if lo > inc_hi {
                    break;
                }
                // Challenger provably better: dethrone.
                if hi < lo.max(0.0).mul_add(0.0, incumbent.expect("set").1) {
                    break;
                }
            }
            if 2.0 * radius < 0.04 {
                break;
            }
        }
        // Validity at the stopping time (the Bet-5 claim).
        if truth(c) < lo || truth(c) > hi {
            worst_miss = true;
        }
        total_adaptive += n;
        per_candidate.push(n);
        let dethrone = match incumbent {
            None => true,
            Some((_, _, inc_hi)) => hi < inc_hi,
        };
        if dethrone {
            incumbent = Some((ci, lo, hi));
        }
    }
    assert!(!worst_miss, "a stopped CS missed its true mean (anytime validity violated)");
    let best_idx = incumbent.expect("incumbent exists").0;
    let x_best = candidates[best_idx];
    assert!(
        (x_best - 0.7).abs() < 0.15,
        "racing missed the optimum region: x* {x_best:.2}"
    );
    // Adaptive allocation: total must be well under uniform max
    // replication, and bad candidates must stop early.
    let max_n = *per_candidate.iter().max().expect("nonempty");
    let total_fixed = max_n * candidates.len() as u64;
    assert!(
        total_adaptive * 2 < total_fixed,
        "racing should eliminate bad candidates early: {total_adaptive} vs uniform {total_fixed}"
    );
    log(
        "anytime-bo",
        "pass",
        &format!(
            "x* {x_best:.2}, adaptive {total_adaptive} vs uniform {total_fixed} samples, all CSs valid at stop"
        ),
    );
}

const GOLDEN_HASH: u64 = 0; // recorded on first run, then frozen

#[test]
fn hetero_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let xs: Vec<Vec<f64>> = (0..10).map(|k| vec![f64::from(k) / 9.0]).collect();
    let ys: Vec<f64> = xs.iter().map(|x| fs_math::det::sin(3.0 * x[0])).collect();
    let noises: Vec<f64> = (0..10).map(|k| if k % 3 == 0 { 0.1 } else { 1e-4 }).collect();
    let gp = Gp::try_fit_diag(&xs, &ys, kernel(), &noises).expect("SPD");
    feed(gp.lml);
    for k in 0..5 {
        let (m, v) = gp.predict(&[0.1f64.mul_add(f64::from(k), 0.05)]);
        feed(m);
        feed(v);
    }
    log("hetero-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "hetero bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
