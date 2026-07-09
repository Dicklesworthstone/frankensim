//! DEV-ONLY statistical battery for the Philox stream (bead frankensim-1za9,
//! item 3). A PractRand-class subset in pure Rust — a test target is isolated
//! from the production dependency graph by construction. Gates are set well
//! beyond the sampling noise (≥ ~4σ) so they flag real defects, not luck.

use fs_rand::{Stream, StreamKey};

fn stream(seed: u64, tile: u32) -> Stream {
    StreamKey {
        seed,
        kernel: 11,
        tile,
    }
    .stream()
}

/// χ² of `m` uniform draws over `k` equal bins — tests marginal uniformity.
#[test]
fn uniform_chi_square() {
    let (k, m) = (100usize, 2_000_000usize);
    let mut s = stream(0xC01D, 1);
    let mut bins = vec![0u64; k];
    for _ in 0..m {
        let b = (s.next_f64() * k as f64) as usize;
        bins[b.min(k - 1)] += 1;
    }
    let exp = m as f64 / k as f64;
    let chi2: f64 = bins.iter().map(|&o| (o as f64 - exp).powi(2) / exp).sum();
    // df = 99: mean 99, sd ≈ 14; the 99.9th percentile is ≈ 148. Gate at 180.
    assert!(chi2 < 180.0, "uniform χ² {chi2:.1} (df=99)");
    println!(
        "{{\"test\":\"uniform-chi2\",\"chi2\":{chi2:.2},\"df\":{}}}",
        k - 1
    );
}

/// Lag-1 serial (Pearson) correlation of the uniform stream — tests independence
/// of consecutive draws.
#[test]
fn serial_correlation_lag1() {
    let m = 1_000_000usize;
    let mut s = stream(0x5E71, 1);
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut prev = s.next_f64();
    for _ in 0..m {
        let cur = s.next_f64();
        sx += prev;
        sy += cur;
        sxx += prev * prev;
        syy += cur * cur;
        sxy += prev * cur;
        prev = cur;
    }
    let n = m as f64;
    let cov = sxy / n - (sx / n) * (sy / n);
    let vx = sxx / n - (sx / n).powi(2);
    let vy = syy / n - (sy / n).powi(2);
    let corr = cov / (vx * vy).sqrt();
    // sd ≈ 1/√m ≈ 0.001; gate at 0.01 (~10σ).
    assert!(corr.abs() < 0.01, "lag-1 correlation {corr:.5}");
    println!("{{\"test\":\"serial-corr\",\"corr\":{corr:.5}}}");
}

/// Monobit balance: over `m` 64-bit words, the fraction of set bits ≈ ½.
#[test]
fn monobit_balance() {
    let m = 200_000usize;
    let mut s = stream(0xB17B, 1);
    let ones: u64 = (0..m).map(|_| u64::from(s.next_u64().count_ones())).sum();
    let total = (m * 64) as f64;
    let frac = ones as f64 / total;
    // sd of the fraction ≈ 0.5/√(64m) ≈ 1.4e-4; gate at 1e-3 (~7σ).
    assert!((frac - 0.5).abs() < 1e-3, "set-bit fraction {frac:.6}");
    println!("{{\"test\":\"monobit\",\"set_fraction\":{frac:.6}}}");
}

/// Inter-stream INDEPENDENCE: streams that differ only by logical identity (tile)
/// must be uncorrelated — the whole point of counter-based, identity-keyed RNG.
#[test]
fn inter_stream_correlation_matrix() {
    let (streams, m) = (8usize, 100_000usize);
    // Draw each stream's sequence.
    let seqs: Vec<Vec<f64>> = (0..streams)
        .map(|tile| {
            let mut s = stream(0x1A9E, tile as u32);
            (0..m).map(|_| s.next_f64()).collect()
        })
        .collect();
    let mut worst = 0.0f64;
    for a in 0..streams {
        for b in (a + 1)..streams {
            let (xa, xb) = (&seqs[a], &seqs[b]);
            let (mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for i in 0..m {
                sx += xa[i];
                sy += xb[i];
                sxx += xa[i] * xa[i];
                syy += xb[i] * xb[i];
                sxy += xa[i] * xb[i];
            }
            let n = m as f64;
            let cov = sxy / n - (sx / n) * (sy / n);
            let vx = sxx / n - (sx / n).powi(2);
            let vy = syy / n - (sy / n).powi(2);
            worst = worst.max((cov / (vx * vy).sqrt()).abs());
        }
    }
    // Per-pair sd ≈ 1/√m ≈ 0.0032; 28 pairs; gate the max at 0.02 (~6σ).
    assert!(worst < 0.02, "worst inter-stream correlation {worst:.5}");
    println!("{{\"test\":\"inter-stream-corr\",\"pairs\":28,\"worst\":{worst:.5}}}");
}
