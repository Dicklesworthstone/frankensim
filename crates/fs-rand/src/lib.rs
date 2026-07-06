//! fs-rand — counter-based Philox streams keyed by LOGICAL work identity
//! (plan §6.7; Decalogue P2's seed pillar).
//!
//! The design that makes e-raced tournaments bit-reproducible and MC results
//! scheduling-independent: a draw is a pure function of
//! `(seed, kernel, tile, index)` — never of which thread ran when. Streams
//! support RANDOM ACCESS by index ([`Stream::at`]), so replay, forking, and
//! out-of-order tile execution cannot perturb randomness.
//!
//! Distributions are built on fs-math's STRICT functions (Box–Muller normals
//! via `det::{ln, sin, cos}`, exponentials via `det::ln`), so even sampled
//! values are cross-ISA bit-deterministic — a deliberate deviation from the
//! plan's ziggurat suggestion (recorded on the bead): ziggurat's table +
//! rejection layout is a PERF follow-up; determinism-first wins v1.
//!
//! Field widths (documented contract): seed 64 bits (Philox key), tile id
//! 32 bits, kernel id 32 bits (together counter words 2–3), draw index 64
//! bits (counter words 0–1). 2⁶⁴ draws per (seed, kernel, tile) stream.

pub mod philox;
pub mod qmc;

use fs_math::det;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The logical identity of a stream: the Cx-carried key (plan §5.2 —
/// "keyed by (seed, kernel_id, tile_id, iteration)", with the iteration as
/// the draw index).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamKey {
    /// Study seed (one of the Five Explicits).
    pub seed: u64,
    /// Kernel identity (registry-assigned; stable across runs).
    pub kernel: u32,
    /// Logical tile identity (NOT the worker/thread id — the whole point).
    pub tile: u32,
}

impl StreamKey {
    /// Open the stream at index 0.
    #[must_use]
    pub fn stream(self) -> Stream {
        Stream {
            key: self,
            index: 0,
        }
    }
}

/// A sequential view over the counter-based generator. `Copy` is deliberate:
/// forking a stream is just copying it (forks that must diverge should use
/// distinct tile/kernel ids instead — divergence by IDENTITY, not by state).
#[derive(Debug, Clone, Copy)]
pub struct Stream {
    key: StreamKey,
    index: u64,
}

impl Stream {
    /// RANDOM ACCESS: the 128 output bits at `index`, independent of any
    /// sequential position. The foundation of replay and shuffle-invariance.
    #[must_use]
    pub fn at(key: StreamKey, index: u64) -> [u32; 4] {
        philox::philox4x32_10(
            [index as u32, (index >> 32) as u32, key.tile, key.kernel],
            [key.seed as u32, (key.seed >> 32) as u32],
        )
    }

    /// Current index (for provenance records / resumable checkpoints).
    #[must_use]
    pub fn index(&self) -> u64 {
        self.index
    }

    /// Resume a stream at a checkpointed index.
    #[must_use]
    pub fn resume(key: StreamKey, index: u64) -> Stream {
        Stream { key, index }
    }

    /// Next 64 uniform bits.
    #[must_use]
    pub fn next_u64(&mut self) -> u64 {
        let block = Self::at(self.key, self.index);
        self.index += 1;
        (u64::from(block[1]) << 32) | u64::from(block[0])
    }

    /// Uniform in [0, 1) with 53 random bits (the standard exact ladder).
    #[must_use]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0) // 2⁻⁵³
    }

    /// Uniform integer in [0, n) via Lemire's widening-multiply method with
    /// the DETERMINISTIC rejection contract: rejected draws advance the
    /// index like any other draw, so the consumed-count is a pure function
    /// of the stream content (replay-safe).
    #[must_use]
    pub fn next_below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "next_below(0) is meaningless");
        loop {
            let x = self.next_u64();
            let m = u128::from(x) * u128::from(n);
            let lo = m as u64;
            if lo >= n.wrapping_neg() % n {
                return (m >> 64) as u64;
            }
            // rejected: index already advanced; try the next block.
        }
    }

    /// Standard normal via Box–Muller on fs-math strict functions —
    /// cross-ISA deterministic sampled values. Consumes exactly 2 draws.
    #[must_use]
    pub fn next_normal(&mut self) -> f64 {
        // u ∈ (0,1]: guard the log; v ∈ [0,1).
        let u = 1.0 - self.next_f64();
        let v = self.next_f64();
        det::sqrt(-2.0 * det::ln(u)) * det::cos(2.0 * std::f64::consts::PI * v)
    }

    /// Exponential(1) via inversion (consumes exactly 1 draw).
    #[must_use]
    pub fn next_exponential(&mut self) -> f64 {
        -det::ln(1.0 - self.next_f64())
    }

    /// Fill a slice with uniform [0,1) values (sequential draws).
    pub fn fill_f64(&mut self, out: &mut [f64]) {
        for v in out {
            *v = self.next_f64();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KEY: StreamKey = StreamKey {
        seed: 0x5EED_0001_DEAD_BEEF,
        kernel: 7,
        tile: 42,
    };

    #[test]
    fn random_access_matches_sequential() {
        let mut s = KEY.stream();
        let seq: Vec<u64> = (0..64).map(|_| s.next_u64()).collect();
        for (i, &want) in seq.iter().enumerate() {
            let block = Stream::at(KEY, i as u64);
            let got = (u64::from(block[1]) << 32) | u64::from(block[0]);
            assert_eq!(got, want, "random access diverged at index {i}");
        }
    }

    #[test]
    fn worker_shuffle_invariance() {
        // Simulate tiles executed in three different worker orders; each
        // tile's draws must be identical regardless (the P2 property that
        // makes results independent of scheduling).
        let tiles: Vec<u32> = (0..16).collect();
        let draw_tile = |tile: u32| -> Vec<f64> {
            let mut s = StreamKey {
                seed: 1234,
                kernel: 3,
                tile,
            }
            .stream();
            (0..32).map(|_| s.next_f64()).collect()
        };
        let baseline: Vec<Vec<f64>> = tiles.iter().map(|&t| draw_tile(t)).collect();
        for order in [
            tiles.iter().rev().copied().collect::<Vec<_>>(),
            tiles
                .iter()
                .step_by(2)
                .chain(tiles.iter().skip(1).step_by(2))
                .copied()
                .collect(),
        ] {
            for &t in &order {
                let redo = draw_tile(t);
                assert!(
                    redo.iter()
                        .zip(&baseline[t as usize])
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "tile {t} draws depended on execution order"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-rand\",\"case\":\"shuffle-invariance\",\"verdict\":\"pass\",\"detail\":\"16 tiles x 3 orders bitwise identical\"}}"
        );
    }

    #[test]
    fn streams_with_different_identities_are_uncorrelated() {
        // Crude cross-correlation check between adjacent tiles/kernels/seeds.
        let corr = |a: StreamKey, b: StreamKey| -> f64 {
            let (mut sa, mut sb) = (a.stream(), b.stream());
            let n = 4096;
            let (mut ma, mut mb, mut cov, mut va, mut vb) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let xs: Vec<(f64, f64)> = (0..n).map(|_| (sa.next_f64(), sb.next_f64())).collect();
            for &(x, y) in &xs {
                ma += x;
                mb += y;
            }
            ma /= f64::from(n);
            mb /= f64::from(n);
            for &(x, y) in &xs {
                cov += (x - ma) * (y - mb);
                va += (x - ma) * (x - ma);
                vb += (y - mb) * (y - mb);
            }
            cov / (va.sqrt() * vb.sqrt())
        };
        for (a, b) in [
            (KEY, StreamKey { tile: 43, ..KEY }),
            (KEY, StreamKey { kernel: 8, ..KEY }),
            (
                KEY,
                StreamKey {
                    seed: KEY.seed ^ 1,
                    ..KEY
                },
            ),
        ] {
            let c = corr(a, b);
            assert!(c.abs() < 0.06, "adjacent identities correlate: {c}");
        }
    }

    #[test]
    fn uniform_chi_square_and_moments() {
        const BINS: usize = 64;
        const N: usize = 64 * 1024;
        let mut s = KEY.stream();
        let mut counts = [0u32; BINS];
        let mut mean = 0.0;
        for _ in 0..N {
            let x = s.next_f64();
            assert!((0.0..1.0).contains(&x));
            counts[(x * BINS as f64) as usize] += 1;
            mean += x;
        }
        mean /= N as f64;
        let expect = (N / BINS) as f64;
        let chi2: f64 = counts
            .iter()
            .map(|&c| (f64::from(c) - expect).powi(2) / expect)
            .sum();
        // 63 dof: mean 63, sd ~11.2; accept within ±5 sd.
        assert!((10.0..=120.0).contains(&chi2), "chi2 {chi2} out of band");
        assert!((mean - 0.5).abs() < 0.005, "mean {mean}");
    }

    #[test]
    fn normal_and_exponential_moments() {
        const N: usize = 200_000;
        let mut s = KEY.stream();
        let (mut m1, mut m2, mut m4) = (0.0, 0.0, 0.0);
        for _ in 0..N {
            let z = s.next_normal();
            m1 += z;
            m2 += z * z;
            m4 += z * z * z * z;
        }
        let n = N as f64;
        assert!((m1 / n).abs() < 0.01, "normal mean {}", m1 / n);
        assert!((m2 / n - 1.0).abs() < 0.02, "normal var {}", m2 / n);
        assert!((m4 / n - 3.0).abs() < 0.12, "normal kurtosis {}", m4 / n);
        let (mut e1, mut e2) = (0.0, 0.0);
        for _ in 0..N {
            let x = s.next_exponential();
            assert!(x >= 0.0);
            e1 += x;
            e2 += x * x;
        }
        assert!((e1 / n - 1.0).abs() < 0.01, "exp mean {}", e1 / n);
        assert!((e2 / n - 2.0).abs() < 0.05, "exp 2nd moment {}", e2 / n);
    }

    #[test]
    fn next_below_is_unbiased_and_replayable() {
        let mut s = KEY.stream();
        let mut counts = [0u32; 7];
        for _ in 0..70_000 {
            counts[s.next_below(7) as usize] += 1;
        }
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (9_500..=10_500).contains(&c),
                "biased bucket {i}: {c} (expect ~10000)"
            );
        }
        // Replay: same key + index range → same values, even through
        // rejection loops (the consumed-count is content-determined).
        let mut a = Stream::resume(KEY, 12345);
        let mut b = Stream::resume(KEY, 12345);
        for _ in 0..1000 {
            assert_eq!(a.next_below(1000), b.next_below(1000));
        }
        assert_eq!(
            a.index(),
            b.index(),
            "rejection consumption must be deterministic"
        );
    }

    #[test]
    fn checkpoint_resume_equality() {
        let mut s = KEY.stream();
        for _ in 0..100 {
            let _ = s.next_normal();
        }
        let ckpt = s.index();
        let tail_a: Vec<u64> = (0..50).map(|_| s.next_u64()).collect();
        let mut resumed = Stream::resume(KEY, ckpt);
        let tail_b: Vec<u64> = (0..50).map(|_| resumed.next_u64()).collect();
        assert_eq!(
            tail_a, tail_b,
            "resume from checkpoint must continue identically"
        );
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}
