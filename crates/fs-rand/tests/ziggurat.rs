//! Ziggurat accuracy + determinism gates (bead frankensim-1za9). Box–Muller is
//! the strict oracle; the ziggurat perf path must match it in DISTRIBUTION and
//! be bit-deterministic (replay-safe) — the bead's own acceptance for the
//! fast-mode path.

use fs_rand::{Stream, StreamKey};

fn stream(seed: u64, tile: u32) -> Stream {
    StreamKey {
        seed,
        kernel: 7,
        tile,
    }
    .stream()
}

#[test]
fn ziggurat_matches_standard_normal_moments() {
    let mut s = stream(0xABCD, 3);
    let m = 2_000_000usize;
    let (mut s1, mut s2, mut s3, mut s4) = (0.0f64, 0.0, 0.0, 0.0);
    for _ in 0..m {
        let x = s.next_normal_ziggurat();
        s1 += x;
        s2 += x * x;
        s3 += x * x * x;
        s4 += x * x * x * x;
    }
    let n = m as f64;
    let mean = s1 / n;
    let var = s2 / n - mean * mean;
    let skew = (s3 / n - 3.0 * mean * var - mean.powi(3)) / var.powf(1.5);
    let exkurt = (s4 / n) / (var * var) - 3.0;
    // Tolerances are ≥ 7σ of the sampling error at m = 2e6 (won't flake).
    assert!(mean.abs() < 0.005, "mean {mean}");
    assert!((var - 1.0).abs() < 0.01, "var {var}");
    assert!(skew.abs() < 0.02, "skew {skew}");
    assert!(exkurt.abs() < 0.05, "exkurt {exkurt}");
    println!(
        "{{\"test\":\"zig-moments\",\"mean\":{mean:.5},\"var\":{var:.5},\
         \"skew\":{skew:.5},\"exkurt\":{exkurt:.5}}}"
    );
}

#[test]
fn ziggurat_is_bit_deterministic() {
    // Same key → bit-identical sequence, replayed (the deterministic-consumption
    // contract holds through the rejection loop).
    let mut a = stream(0x5151, 3);
    let mut b = stream(0x5151, 3);
    for _ in 0..100_000 {
        assert_eq!(
            a.next_normal_ziggurat().to_bits(),
            b.next_normal_ziggurat().to_bits()
        );
    }
    // Divergence is by IDENTITY (a different tile), not by shared state.
    let mut c = stream(0x5151, 4);
    let mut d = stream(0x5151, 3);
    let diff = (0..2000)
        .filter(|_| c.next_normal_ziggurat().to_bits() != d.next_normal_ziggurat().to_bits())
        .count();
    assert!(
        diff > 1900,
        "streams should diverge by identity: {diff}/2000"
    );
}

#[test]
fn ziggurat_agrees_with_box_muller_ks() {
    // Two-sample Kolmogorov–Smirnov: the ziggurat and the strict Box–Muller path
    // are both N(0,1), so their empirical CDFs must nearly coincide.
    let m = 200_000usize;
    let mut sz = stream(0x2222, 3);
    let mut sb = stream(0x3333, 3);
    let mut z: Vec<f64> = (0..m).map(|_| sz.next_normal_ziggurat()).collect();
    let mut bm: Vec<f64> = (0..m).map(|_| sb.next_normal()).collect();
    z.sort_by(|a, b| a.partial_cmp(b).unwrap());
    bm.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (mut i, mut j, mut ks) = (0usize, 0usize, 0.0f64);
    while i < m && j < m {
        if z[i] <= bm[j] {
            i += 1;
        } else {
            j += 1;
        }
        ks = ks.max((i as f64 - j as f64).abs() / m as f64);
    }
    // KS α≈0.001 critical value ≈ 1.95·√(2/m) ≈ 0.0062; gate generously at 0.02.
    assert!(
        ks < 0.02,
        "two-sample KS {ks} — ziggurat disagrees with Box–Muller"
    );
    println!("{{\"test\":\"zig-vs-boxmuller-ks\",\"ks\":{ks:.5},\"m\":{m}}}");
}
