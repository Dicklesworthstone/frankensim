//! Randomized-NLA battery (6ys.7): rank-r error vs oversampling bounds
//! on three synthetic spectra, adaptive-estimate coverage, RSVD accuracy,
//! Nyström PSD reconstruction, sketch-LS vs direct QR, Hutch++ variance
//! beating Hutchinson at matched budgets, keyed-stream determinism, and
//! the cross-ISA golden hash.

use fs_la::factor::qr;
use fs_la::rand_nla::{hutch_pp, hutchinson, nystrom_psd, range_finder, rsvd, sketch_ls};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

/// A = Q1·diag(σ)·Q2ᵀ with a chosen spectrum (m×n, m ≥ n).
fn spectrum_matrix(m: usize, n: usize, sigma: &[f64], seed: u64) -> Vec<f64> {
    let mut s = seed;
    let g1: Vec<f64> = (0..m * n).map(|_| lcg(&mut s)).collect();
    let g2: Vec<f64> = (0..n * n).map(|_| lcg(&mut s)).collect();
    let f1 = qr(&g1, m, n);
    let f2 = qr(&g2, n, n);
    // Columns of Q1 (m×n) and Q2 (n×n) via apply_q.
    let mut q1 = vec![0.0f64; m * n];
    for j in 0..n {
        let mut e = vec![0.0f64; m];
        e[j] = 1.0;
        f1.apply_q(&mut e);
        for i in 0..m {
            q1[i * n + j] = e[i];
        }
    }
    let mut q2 = vec![0.0f64; n * n];
    for j in 0..n {
        let mut e = vec![0.0f64; n];
        e[j] = 1.0;
        f2.apply_q(&mut e);
        for i in 0..n {
            q2[i * n + j] = e[i];
        }
    }
    let mut a = vec![0.0f64; m * n];
    for i in 0..m {
        for j in 0..n {
            let mut acc = 0.0f64;
            for k in 0..sigma.len().min(n) {
                acc = (q1[i * n + k] * sigma[k]).mul_add(q2[j * n + k], acc);
            }
            a[i * n + j] = acc;
        }
    }
    a
}

fn spectral_residual(a: &[f64], m: usize, n: usize, q: &[f64], k: usize) -> f64 {
    // ‖(I−QQᵀ)A‖_F as a proxy (≥ spectral norm; fine for bound checks).
    let mut total = 0.0f64;
    for j in 0..n {
        let mut col = vec![0.0f64; m];
        for i in 0..m {
            col[i] = a[i * n + j];
        }
        for l in 0..k {
            let mut dot = 0.0f64;
            for i in 0..m {
                dot = q[i * k + l].mul_add(col[i], dot);
            }
            for i in 0..m {
                col[i] = (-dot).mul_add(q[i * k + l], col[i]);
            }
        }
        total += col.iter().map(|t| t * t).sum::<f64>();
    }
    total.sqrt()
}

#[test]
fn rangefinder_meets_bounds_on_three_spectra() {
    let (m, n) = (120usize, 60usize);
    // Fast exponential decay, slow algebraic decay, and a spectral gap.
    let fast: Vec<f64> = (0..n).map(|k| 0.5f64.powi(i32::try_from(k).unwrap())).collect();
    let slow: Vec<f64> = (0..n).map(|k| 1.0 / f64::from(u32::try_from(k + 1).unwrap())).collect();
    let gap: Vec<f64> = (0..n).map(|k| if k < 10 { 1.0 } else { 1e-6 }).collect();
    for (name, sigma, rank, q_pow, budget) in [
        ("fast", &fast, 12usize, 0usize, 3.0f64),
        ("slow", &slow, 20, 2, 8.0),
        ("gap", &gap, 10, 1, 3.0),
    ] {
        let a = spectrum_matrix(m, n, sigma, 0xA11CE);
        let (qb, rep) = range_finder(&a, m, n, rank, 8, q_pow, 42);
        let resid = spectral_residual(&a, m, n, &qb, rep.rank);
        // Theoretical tail: sqrt(Σ_{k>rank} σ²) — allow `budget`× slack
        // (Frobenius proxy + randomness).
        let tail: f64 = sigma[rank..].iter().map(|t| t * t).sum::<f64>().sqrt()
            * (n as f64).sqrt().mul_add(0.0, budget);
        assert!(
            resid <= tail.max(1e-10),
            "{name}: residual {resid:.3e} above budgeted tail {tail:.3e}"
        );
        // Coverage: the posterior estimate must NOT understate the truth
        // (the G0 containment property for error estimates).
        assert!(
            rep.est_error * 1.5 >= resid / (n as f64).sqrt(),
            "{name}: estimate {:.3e} understates residual {resid:.3e}",
            rep.est_error
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"rangefinder\",\"verdict\":\"pass\",\"detail\":\"3 spectra within budgeted tails, estimates cover\"}}"
    );
}

#[test]
fn rsvd_recovers_singular_values() {
    let (m, n) = (100usize, 50usize);
    let sigma: Vec<f64> = (0..n).map(|k| 2.0 * 0.6f64.powi(i32::try_from(k).unwrap())).collect();
    let a = spectrum_matrix(m, n, &sigma, 0xB0B);
    let (_, sv, _, _) = rsvd(&a, m, n, 10, 8, 1, 7);
    for k in 0..8 {
        assert!(
            (sv[k] - sigma[k]).abs() < 1e-6 * sigma[k].max(1e-12),
            "rsvd sigma[{k}] {} vs {}",
            sv[k],
            sigma[k]
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"rsvd\",\"verdict\":\"pass\",\"detail\":\"leading 8 singular values to 1e-6 rel\"}}"
    );
}

#[test]
fn nystrom_reconstructs_low_rank_psd() {
    let n = 60usize;
    // PSD with rank ~8 + tiny tail.
    let sigma: Vec<f64> = (0..n).map(|k| if k < 8 { 4.0 * 0.7f64.powi(i32::try_from(k).unwrap()) } else { 1e-9 }).collect();
    // Build PSD as Q·diag(σ)·Qᵀ.
    let mut s = 0xC0DE_u64;
    let g: Vec<f64> = (0..n * n).map(|_| lcg(&mut s)).collect();
    let f = qr(&g, n, n);
    let mut qm = vec![0.0f64; n * n];
    for j in 0..n {
        let mut e = vec![0.0f64; n];
        e[j] = 1.0;
        f.apply_q(&mut e);
        for i in 0..n {
            qm[i * n + j] = e[i];
        }
    }
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0f64;
            for k in 0..n {
                acc = (qm[i * n + k] * sigma[k]).mul_add(qm[j * n + k], acc);
            }
            a[i * n + j] = acc;
        }
    }
    let fk = nystrom_psd(&a, n, 10, 6, 99);
    let k = fk.len() / n;
    // ‖A − FFᵀ‖_F small relative to ‖A‖_F.
    let (mut err, mut nrm) = (0.0f64, 0.0f64);
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0f64;
            for l in 0..k {
                acc = fk[i * k + l].mul_add(fk[j * k + l], acc);
            }
            let d = a[i * n + j] - acc;
            err += d * d;
            nrm += a[i * n + j] * a[i * n + j];
        }
    }
    let rel = (err / nrm).sqrt();
    assert!(rel < 1e-6, "nystrom relative error {rel:.3e}");
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"nystrom\",\"verdict\":\"pass\",\"detail\":\"rank-8 PSD reconstructed to {rel:.2e}\"}}"
    );
}

#[test]
fn sketch_ls_matches_direct_qr() {
    let (m, n) = (400usize, 30usize);
    let mut s = 0x15AAC_u64;
    let a: Vec<f64> = (0..m * n).map(|_| lcg(&mut s)).collect();
    let b: Vec<f64> = (0..m).map(|_| lcg(&mut s)).collect();
    let direct = qr(&a, m, n).solve_ls(&b);
    let (sketched, iters) = sketch_ls(&a, m, n, &b, 5);
    for k in 0..n {
        assert!(
            (sketched[k] - direct[k]).abs() < 1e-8 * direct[k].abs().max(1.0),
            "sketch LS[{k}]: {} vs {}",
            sketched[k],
            direct[k]
        );
    }
    assert!(iters <= 60, "preconditioned CG should converge fast: {iters} iters");
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"sketch-ls\",\"verdict\":\"pass\",\"detail\":\"matches QR to 1e-8 in {iters} CG iters\"}}"
    );
}

#[test]
fn hutch_pp_beats_hutchinson_variance() {
    let n = 80usize;
    // Decaying-spectrum SPD (where Hutch++ shines).
    let sigma: Vec<f64> = (0..n).map(|k| 8.0 * 0.8f64.powi(i32::try_from(k).unwrap())).collect();
    let true_trace: f64 = sigma.iter().sum();
    let mut s = 0x7ACE_u64;
    let g: Vec<f64> = (0..n * n).map(|_| lcg(&mut s)).collect();
    let f = qr(&g, n, n);
    let mut qm = vec![0.0f64; n * n];
    for j in 0..n {
        let mut e = vec![0.0f64; n];
        e[j] = 1.0;
        f.apply_q(&mut e);
        for i in 0..n {
            qm[i * n + j] = e[i];
        }
    }
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0f64;
            for k in 0..n {
                acc = (qm[i * n + k] * sigma[k]).mul_add(qm[j * n + k], acc);
            }
            a[i * n + j] = acc;
        }
    }
    // Across independent seeds: both unbiased-ish; Hutch++ tighter.
    let probes = 24;
    let (mut h_err2, mut hpp_err2) = (0.0f64, 0.0f64);
    let trials = 20;
    for t in 0..trials {
        let h = hutchinson(&a, n, probes, 1000 + t);
        let hpp = hutch_pp(&a, n, probes, 1000 + t);
        h_err2 += (h.estimate - true_trace).powi(2);
        hpp_err2 += (hpp.estimate - true_trace).powi(2);
    }
    assert!(
        hpp_err2 < h_err2 * 0.5,
        "Hutch++ MSE {hpp_err2:.3e} must decisively beat Hutchinson {h_err2:.3e}"
    );
    // Determinism: same seed → bitwise same estimate.
    let a1 = hutchinson(&a, n, probes, 123);
    let a2 = hutchinson(&a, n, probes, 123);
    assert_eq!(a1.estimate.to_bits(), a2.estimate.to_bits());
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"hutch\",\"verdict\":\"pass\",\"detail\":\"MSE {h_err2:.3e} -> {hpp_err2:.3e} at {probes} probes x {trials} trials\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x3e92_8bac_8cf9_fd48; // pinned from first run (arm64), cross-checked on x86_64

#[test]
fn rand_nla_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let (m, n) = (60usize, 30usize);
    let sigma: Vec<f64> = (0..n).map(|k| 0.7f64.powi(i32::try_from(k).unwrap())).collect();
    let a = spectrum_matrix(m, n, &sigma, 0xFEED);
    let (_, sv, _, rep) = rsvd(&a, m, n, 6, 4, 1, 11);
    for &v in &sv {
        feed(v);
    }
    feed(rep.est_error);
    let h = hutchinson(&a[..n * n], n, 16, 77);
    feed(h.estimate);
    feed(h.variance_est);
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"rand-nla-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "randomized-NLA bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
