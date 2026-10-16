//! Distribution battery (6ys.19): moment gates within CLT bands,
//! deterministic-consumption replay contracts, bitwise alias
//! construction, vMF geometry invariants, truncation guarantees, and the
//! cross-ISA golden hash.

use fs_rand::dist::AliasTable;
use fs_rand::{Stream, StreamCheckpoint, StreamKey};

const KEY: StreamKey = StreamKey {
    seed: 0xD157_0001,
    kernel: 11,
    tile: 3,
};

fn resume(index: u64) -> Stream {
    let retained = StreamCheckpoint::current(KEY, index).to_canonical_le_bytes();
    Stream::resume_retained(&retained)
        .expect("fixture checkpoints use the current canonical replay transport")
}

#[test]
fn gamma_moments_and_replay() {
    const N: usize = 200_000;
    for &alpha in &[0.5f64, 1.0, 2.5, 7.0] {
        let mut s = KEY.stream();
        let (mut m1, mut m2) = (0.0f64, 0.0f64);
        for _ in 0..N {
            let g = s.next_gamma(alpha);
            assert!(g > 0.0, "gamma must be positive");
            m1 += g;
            m2 += g * g;
        }
        let n = N as f64;
        let mean = m1 / n;
        let var = m2 / n - mean * mean;
        // CLT bands: sd(mean) = sqrt(var/n); allow 5σ.
        let mean_tol = 5.0 * (alpha / n).sqrt();
        assert!(
            (mean - alpha).abs() < mean_tol,
            "gamma({alpha}) mean {mean} vs {alpha} (tol {mean_tol})"
        );
        assert!(
            (var - alpha).abs() < 0.05 * alpha + mean_tol * 4.0,
            "gamma({alpha}) var {var} vs {alpha}"
        );
        // Deterministic-consumption replay: same key + index → same value
        // AND same post-index, even mid-stream.
        let mut a = resume(777);
        let _ = a.next_f64(); // interleave
        let idx_before = a.index();
        let va = a.next_gamma(alpha);
        let consumed = a.index() - idx_before;
        let mut b = resume(idx_before);
        let vb = b.next_gamma(alpha);
        assert_eq!(va.to_bits(), vb.to_bits(), "gamma replay value");
        assert_eq!(b.index() - idx_before, consumed, "gamma replay consumption");
    }
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"gamma\",\"verdict\":\"pass\",\"detail\":\"4 shapes: moments in CLT bands + replay contract\"}}"
    );
}

#[test]
fn beta_and_dirichlet_moments() {
    const N: usize = 100_000;
    let mut s = KEY.stream();
    let (a, b) = (2.0f64, 5.0f64);
    let mut m1 = 0.0f64;
    for _ in 0..N {
        let x = s.next_beta(a, b);
        assert!((0.0..=1.0).contains(&x));
        m1 += x;
    }
    let want = a / (a + b);
    assert!(
        (m1 / N as f64 - want).abs() < 0.003,
        "beta mean {} vs {want}",
        m1 / N as f64
    );
    // Dirichlet: components sum to 1; means proportional to alphas.
    let alphas = [1.0f64, 2.0, 4.0];
    let total: f64 = alphas.iter().sum();
    let mut sums = [0.0f64; 3];
    let mut out = [0.0f64; 3];
    for _ in 0..N {
        s.next_dirichlet(&alphas, &mut out);
        let sum: f64 = out.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "dirichlet must sum to 1: {sum}");
        for (acc, &v) in sums.iter_mut().zip(&out) {
            *acc += v;
        }
    }
    for (k, (&acc, &alpha)) in sums.iter().zip(&alphas).enumerate() {
        let want = alpha / total;
        assert!(
            (acc / N as f64 - want).abs() < 0.004,
            "dirichlet mean[{k}] {} vs {want}",
            acc / N as f64
        );
    }
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"beta-dirichlet\",\"verdict\":\"pass\",\"detail\":\"beta(2,5) + dirichlet(1,2,4) moments in band\"}}"
    );
}

#[test]
fn extreme_beta_and_dirichlet_shapes_are_total_and_replay_safe() {
    const TINY: f64 = f64::MIN_POSITIVE;
    const SMALLEST: f64 = f64::from_bits(1);
    const HUGE: f64 = f64::MAX;

    for (case, &(alpha, beta)) in [
        (SMALLEST, SMALLEST),
        (TINY, TINY),
        (TINY, 2.0 * TINY),
        (1.0e-300, TINY),
        (HUGE, HUGE),
    ]
    .iter()
    .enumerate()
    {
        let start = 120_000 + 100 * case as u64;
        let mut a = resume(start);
        let value = a.next_beta(alpha, beta);
        let end = a.index();

        assert!(value.is_finite(), "extreme-shape beta returned {value}");
        assert!((0.0..=1.0).contains(&value), "beta out of range: {value}");

        let mut replay = resume(start);
        let replayed = replay.next_beta(alpha, beta);
        assert_eq!(value.to_bits(), replayed.to_bits(), "beta replay bits");
        assert_eq!(end, replay.index(), "beta fallback consumed hidden draws");
    }

    for (case, alphas) in [
        [SMALLEST, SMALLEST, SMALLEST],
        [TINY, TINY, TINY],
        [TINY, 2.0 * TINY, 1.0e-300],
        [HUGE, HUGE, HUGE],
    ]
    .into_iter()
    .enumerate()
    {
        let start = 140_000 + 100 * case as u64;
        let mut a = resume(start);
        let mut values = [0.0; 3];
        a.next_dirichlet(&alphas, &mut values);
        let end = a.index();

        assert!(
            values.iter().all(|v| v.is_finite() && *v >= 0.0),
            "extreme-shape Dirichlet returned {values:?}"
        );
        let sum: f64 = values.iter().sum();
        assert!(sum > 0.0, "Dirichlet normalization was all zero");
        assert!((sum - 1.0).abs() <= f64::EPSILON, "simplex sum {sum}");

        let mut replay = resume(start);
        let mut replayed = [0.0; 3];
        replay.next_dirichlet(&alphas, &mut replayed);
        assert_eq!(
            values.map(f64::to_bits),
            replayed.map(f64::to_bits),
            "Dirichlet replay bits"
        );
        assert_eq!(
            end,
            replay.index(),
            "Dirichlet fallback consumed hidden draws"
        );
    }
}

#[test]
fn alias_table_bitwise_construction_and_chi_square() {
    const N: usize = 200_000;
    let weights = [1.0f64, 2.0, 3.0, 10.0, 0.5];
    let t1 = AliasTable::new(&weights);
    let t2 = AliasTable::new(&weights);
    // Bitwise-identical construction (P2 on setup).
    for i in 0..weights.len() {
        let mut s1 = resume(40_000 + i as u64);
        let mut s2 = resume(40_000 + i as u64);
        assert_eq!(
            t1.sample(&mut s1),
            t2.sample(&mut s2),
            "tables must behave identically"
        );
    }
    // Single-draw consumption contract.
    let mut s = KEY.stream();
    let before = s.index();
    let _ = t1.sample(&mut s);
    assert_eq!(
        s.index() - before,
        1,
        "alias sampling must consume exactly 1 draw"
    );
    // Chi-square against the pmf.
    let total: f64 = weights.iter().sum();
    let mut counts = [0u32; 5];
    let mut st = resume(90_000);
    for _ in 0..N {
        counts[t1.sample(&mut st)] += 1;
    }
    let mut chi2 = 0.0f64;
    for (c, &w) in counts.iter().zip(&weights) {
        let expect = N as f64 * w / total;
        chi2 += (f64::from(*c) - expect).powi(2) / expect;
    }
    // 4 dof: mean 4, sd ~2.8; accept generously (deterministic seed).
    assert!(
        chi2 < 25.0,
        "alias chi-square {chi2} out of band; counts {counts:?}"
    );
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"alias\",\"verdict\":\"pass\",\"detail\":\"bitwise construction, 1-draw contract, chi2 {chi2:.2}\"}}"
    );
}

#[test]
fn vmf_geometry_and_fixed_consumption() {
    const N: usize = 50_000;
    let mu = {
        // A deliberately non-axis mean direction, normalized.
        let raw = [0.3f64, -0.5, 0.81];
        let n = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
        [raw[0] / n, raw[1] / n, raw[2] / n]
    };
    for &kappa in &[1.0f64, 10.0, 100.0] {
        let mut s = KEY.stream();
        let mut resultant = [0.0f64; 3];
        for _ in 0..N {
            let before = s.index();
            let v = s.next_vmf3(mu, kappa);
            assert_eq!(s.index() - before, 2, "vMF must consume exactly 2 draws");
            let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            assert!(
                (norm - 1.0).abs() < 1e-12,
                "vMF output must be unit: {norm}"
            );
            for (r, &c) in resultant.iter_mut().zip(&v) {
                *r += c;
            }
        }
        let rlen = (resultant[0].powi(2) + resultant[1].powi(2) + resultant[2].powi(2)).sqrt();
        let dot = (resultant[0] * mu[0] + resultant[1] * mu[1] + resultant[2] * mu[2]) / rlen;
        assert!(
            dot > 0.999,
            "mean resultant must align with mu at kappa={kappa}: dot {dot}"
        );
        // Analytic mean resultant length for vMF on S²: coth(κ) − 1/κ.
        let coth = (fs_math::det::exp(kappa) + fs_math::det::exp(-kappa))
            / (fs_math::det::exp(kappa) - fs_math::det::exp(-kappa));
        let want = coth - 1.0 / kappa;
        let got = rlen / N as f64;
        assert!(
            (got - want).abs() < 0.02,
            "resultant length at kappa={kappa}: {got} vs analytic {want}"
        );
    }
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"vmf\",\"verdict\":\"pass\",\"detail\":\"unit norm, 2-draw contract, resultant matches coth(k)-1/k at 3 kappas\"}}"
    );
}

#[test]
fn truncated_variants_respect_bounds() {
    const N: usize = 100_000;
    let mut s = KEY.stream();
    // Truncated exponential: in [0, cap], exactly 1 draw each.
    let cap = 1.5f64;
    for _ in 0..N / 10 {
        let before = s.index();
        let x = s.next_truncated_exponential(cap);
        assert_eq!(s.index() - before, 1);
        assert!((0.0..=cap).contains(&x), "trunc-exp out of range: {x}");
    }
    // Truncated normal at lo = 2.0: all samples ≥ lo; mean matches the
    // analytic φ(lo)/(1−Φ(lo)) (computed via the landed erf).
    let lo = 2.0f64;
    let mut m1 = 0.0f64;
    for _ in 0..N {
        let z = s.next_truncated_normal(lo);
        assert!(z >= lo, "truncated normal below bound: {z}");
        m1 += z;
    }
    let phi = fs_math::det::exp(-0.5 * lo * lo) / (2.0 * std::f64::consts::PI).sqrt();
    let tail = 0.5 * fs_math::det::erfc(lo / std::f64::consts::SQRT_2);
    let want = phi / tail;
    let got = m1 / N as f64;
    assert!(
        (got - want).abs() < 0.01,
        "truncated-normal mean {got} vs analytic {want}"
    );
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"truncated\",\"verdict\":\"pass\",\"detail\":\"bounds hold; mean {got:.4} vs analytic {want:.4}\"}}"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const DIST_STREAM_V1_AGGREGATE_FNV1A64_GOLDEN: u64 = 0x4224_6e28_56de_673c;

#[test]
fn dist_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let mut s = KEY.stream();
    for _ in 0..2000 {
        feed(s.next_gamma(2.5));
        feed(s.next_beta(2.0, 5.0));
        feed(s.next_truncated_normal(1.0));
        let v = s.next_vmf3([0.0, 0.0, 1.0], 5.0);
        feed(v[0]);
        feed(v[1]);
        feed(v[2]);
    }
    let t = AliasTable::new(&[1.0, 2.0, 3.0, 10.0, 0.5]);
    for _ in 0..2000 {
        feed(t.sample(&mut s) as f64);
    }
    println!(
        "{{\"suite\":\"fs-rand\",\"case\":\"dist-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, DIST_STREAM_V1_AGGREGATE_FNV1A64_GOLDEN,
        "distribution bits changed: {acc:#018x} vs {DIST_STREAM_V1_AGGREGATE_FNV1A64_GOLDEN:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
