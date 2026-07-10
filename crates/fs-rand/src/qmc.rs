//! Quasi–Monte Carlo (plan §6.7): base-2 Sobol' sequences with TRUE Owen
//! nested-uniform scrambling, and rank-1 lattice rules by CBC construction.
//!
//! Why QMC: at FrankenSim's UQ dimensionalities, low-discrepancy points buy
//! 1–2 orders of magnitude over plain MC; Owen scrambling adds unbiasedness
//! and RMSE benefits while PRESERVING the net structure. Both are gated by
//! convergence tests here, not vibes.
//!
//! Owen scrambling done RIGHT via counter-based randomness: nested uniform
//! scrambling assigns an independent flip bit to every node of the binary
//! digit tree. We derive each node's bit from Philox keyed by
//! (seed, dimension) with the counter encoding (bit-depth, prefix-path) —
//! a lazily-materialized random tree with ZERO storage, random-access
//! replayable like everything else in fs-rand. (Hash-based approximations
//! in the Laine–Karras lineage are a recorded PERF follow-up; correctness
//! first.)
//!
//! Direction numbers: the embedded head (dimensions 1..=10) of the Joe–Kuo
//! D6 table; construction preconditions (m_k odd, m_k < 2^k) are ASSERTED
//! at generator creation so a mistyped table fails loudly, and the exact
//! per-dimension stratification tests catch any surviving corruption.
//! The full 21201-dimension table is a recorded data-import follow-up.

use crate::{Stream, StreamKey};

/// Maximum embedded dimension of the v1 table.
pub const MAX_SOBOL_DIM: usize = 10;

/// Joe–Kuo D6 head: (s = degree, a = coefficient bits, m = initial values).
/// Dimension 1 is the van der Corput sequence (handled specially).
const JOE_KUO: [(u32, u32, &[u32]); 9] = [
    (1, 0, &[1]),
    (2, 1, &[1, 3]),
    (3, 1, &[1, 3, 1]),
    (3, 2, &[1, 1, 1]),
    (4, 1, &[1, 1, 3, 3]),
    (4, 4, &[1, 3, 5, 13]),
    (5, 2, &[1, 1, 5, 5, 17]),
    (5, 4, &[1, 1, 5, 5, 5]),
    (5, 7, &[1, 1, 7, 11, 19]),
];

const BITS: u32 = 32;

/// A Sobol' generator over `dim` dimensions with optional Owen scrambling.
#[derive(Debug, Clone)]
pub struct Sobol {
    /// Direction vectors: `v[d][k]` for dimension d, bit k (as 32-bit
    /// binary fractions).
    directions: Vec<[u32; BITS as usize]>,
    /// Owen scrambling seed; `None` = unscrambled net.
    scramble: Option<u64>,
}

impl Sobol {
    /// Unscrambled Sobol' sequence in `dim` dimensions (1..=[`MAX_SOBOL_DIM`]).
    ///
    /// # Panics
    /// If `dim` is 0 or exceeds the embedded table, or if the table violates
    /// the direction-number preconditions (a corrupted table must fail
    /// LOUDLY, not generate a subtly broken net).
    #[must_use]
    pub fn new(dim: usize) -> Sobol {
        assert!(
            (1..=MAX_SOBOL_DIM).contains(&dim),
            "dim {dim} outside 1..={MAX_SOBOL_DIM} (embedded Joe-Kuo head; larger tables are a \
             recorded follow-up)"
        );
        let mut directions = Vec::with_capacity(dim);
        // Dimension 1: van der Corput — v_k = 2^(31-k).
        let mut v0 = [0u32; BITS as usize];
        for (k, v) in v0.iter_mut().enumerate() {
            *v = 1 << (31 - k);
        }
        directions.push(v0);
        for d in 1..dim {
            let (s, a, m) = JOE_KUO[d - 1];
            let s = s as usize;
            for (k, &mk) in m.iter().enumerate() {
                assert!(mk % 2 == 1, "dim {}: m[{k}]={mk} must be odd", d + 1);
                assert!(
                    mk < (2 << k),
                    "dim {}: m[{k}]={mk} must be < 2^{}",
                    d + 1,
                    k + 1
                );
            }
            let mut v = [0u32; BITS as usize];
            for k in 0..BITS as usize {
                if k < s {
                    v[k] = m[k] << (31 - k);
                } else {
                    // Recurrence: v_k = v_{k-s} ^ (v_{k-s} >> s) ^ Σ a_i v_{k-i}.
                    let mut val = v[k - s] ^ (v[k - s] >> s);
                    for i in 1..s {
                        if (a >> (s - 1 - i)) & 1 == 1 {
                            val ^= v[k - i];
                        }
                    }
                    v[k] = val;
                }
            }
            directions.push(v);
        }
        Sobol {
            directions,
            scramble: None,
        }
    }

    /// Owen-scrambled variant (nested uniform scrambling, seed-replayable).
    #[must_use]
    pub fn scrambled(dim: usize, seed: u64) -> Sobol {
        let mut s = Sobol::new(dim);
        s.scramble = Some(seed);
        s
    }

    /// Number of dimensions.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.directions.len()
    }

    /// The n-th point (RANDOM ACCESS; n up to 2³² − 1): Gray-code XOR of
    /// direction vectors, then optional Owen scrambling, then the [0,1)
    /// float ladder.
    pub fn point(&self, n: u32, out: &mut [f64]) {
        assert_eq!(out.len(), self.dim(), "output slice must match dim");
        let gray = n ^ (n >> 1);
        for (d, slot) in out.iter_mut().enumerate() {
            let mut x = 0u32;
            for k in 0..BITS {
                if (gray >> k) & 1 == 1 {
                    x ^= self.directions[d][k as usize];
                }
            }
            if let Some(seed) = self.scramble {
                x = owen_scramble(x, seed, d as u32);
            }
            // x / 2^32: exact ladder into [0, 1).
            *slot = f64::from(x) / 4_294_967_296.0;
        }
    }

    /// Convenience: materialize the first `n` points row-major.
    #[must_use]
    pub fn points(&self, n: u32) -> Vec<f64> {
        let d = self.dim();
        let mut out = vec![0.0; n as usize * d];
        for i in 0..n {
            let start = i as usize * d;
            self.point(i, &mut out[start..start + d]);
        }
        out
    }
}

/// TRUE nested-uniform (Owen) scrambling of a 32-bit digit string: bit b's
/// flip decision is an independent Bernoulli(1/2) determined by the PREFIX
/// (bits above b) — a random binary tree, lazily derived from Philox with
/// counter = (bit index, prefix value), key = (seed, dimension). Determinism
/// and random access come free; no tree is ever stored.
fn owen_scramble(x: u32, seed: u64, dim: u32) -> u32 {
    let key = StreamKey {
        seed,
        kernel: 0x0E11,
        tile: dim,
    };
    let mut y = 0u32;
    for b in 0..BITS {
        // Prefix = the bits ABOVE position b (b=31 is the most significant).
        let bit_pos = 31 - b; // process MSB first
        let prefix = if b == 0 { 0 } else { x >> (bit_pos + 1) };
        // Counter encodes (level, prefix) — unique per tree node.
        let idx = (u64::from(b) << 32) | u64::from(prefix);
        let flip = Stream::at(key, idx)[0] & 1;
        let bit = ((x >> bit_pos) & 1) ^ flip;
        y |= bit << bit_pos;
    }
    y
}

// ---------------------------------------------------------------------------
// Rank-1 lattice rules (CBC construction).
// ---------------------------------------------------------------------------

/// A rank-1 lattice rule: points x_k = frac(k · z / n), k = 0..n.
#[derive(Debug, Clone)]
pub struct Lattice {
    /// Number of points.
    pub n: u32,
    /// Generating vector.
    pub z: Vec<u32>,
}

impl Lattice {
    /// Component-by-component construction in the product-weighted Korobov
    /// space with the Bernoulli-B₂ kernel (weights γ_j = 1): for each
    /// dimension, choose z_j minimizing the worst-case error given the
    /// previously fixed components. O(n²·dim) — fine at tabled sizes.
    ///
    /// # Panics
    /// If `n < 3` or `dim` is 0.
    #[must_use]
    pub fn cbc(n: u32, dim: usize) -> Lattice {
        assert!(n >= 3 && dim >= 1, "lattice needs n >= 3, dim >= 1");
        // B₂(x) = x² − x + 1/6; kernel ω(x) = 1 + γ·B₂({x}).
        let b2 = |x: f64| x * x - x + 1.0 / 6.0;
        let nf = f64::from(n);
        // prod[k] = Π_j (1 + B₂({k z_j / n})) over chosen dims.
        let mut prod = vec![1.0f64; n as usize];
        let mut z = Vec::with_capacity(dim);
        for _ in 0..dim {
            let mut best = (f64::INFINITY, 1u32);
            for cand in 1..n {
                if gcd(cand, n) != 1 {
                    continue;
                }
                // e²(z, cand) ∝ Σ_k prod[k]·(1 + B₂({k·cand/n}))
                let mut err = 0.0;
                for k in 0..n {
                    let frac =
                        f64::from((u64::from(k) * u64::from(cand) % u64::from(n)) as u32) / nf;
                    err += prod[k as usize] * (1.0 + b2(frac));
                }
                // Deterministic tie-breaking: strictly-less keeps lowest cand.
                if err < best.0 {
                    best = (err, cand);
                }
            }
            let chosen = best.1;
            for k in 0..n {
                let frac = f64::from((u64::from(k) * u64::from(chosen) % u64::from(n)) as u32) / nf;
                prod[k as usize] *= 1.0 + b2(frac);
            }
            z.push(chosen);
        }
        Lattice { n, z }
    }

    /// The k-th point.
    pub fn point(&self, k: u32, out: &mut [f64]) {
        assert_eq!(out.len(), self.z.len());
        for (j, slot) in out.iter_mut().enumerate() {
            let prod = u64::from(k) * u64::from(self.z[j]) % u64::from(self.n);
            *slot = f64::from(prod as u32) / f64::from(self.n);
        }
    }

    /// The squared worst-case error in the (γ=1) Korobov space — the CBC
    /// objective, exposed for the convergence-rate diagnostic.
    #[must_use]
    pub fn korobov_error_sq(&self) -> f64 {
        let b2 = |x: f64| x * x - x + 1.0 / 6.0;
        let nf = f64::from(self.n);
        let mut sum = 0.0;
        for k in 0..self.n {
            let mut prod = 1.0;
            for &zj in &self.z {
                let frac =
                    f64::from((u64::from(k) * u64::from(zj) % u64::from(self.n)) as u32) / nf;
                prod *= 1.0 + b2(frac);
            }
            sum += prod;
        }
        sum / nf - 1.0
    }
}

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Baker's transformation (tent map) periodization for non-periodic
/// integrands on lattices: φ(x) = 1 − |2x − 1|.
#[must_use]
pub fn baker(x: f64) -> f64 {
    1.0 - (2.0 * x - 1.0).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// EXACT stratification: for every dimension of a valid base-2 Sobol
    /// net, the first 2^m points put EXACTLY one point in each dyadic bin
    /// [i/2^m, (i+1)/2^m). This mechanically catches direction-table errors.
    #[test]
    fn per_dimension_exact_stratification() {
        let s = Sobol::new(MAX_SOBOL_DIM);
        let mut buf = vec![0.0; MAX_SOBOL_DIM];
        for m in 1..=8u32 {
            let count = 1u32 << m;
            let mut bins = vec![vec![0u32; count as usize]; MAX_SOBOL_DIM];
            for i in 0..count {
                s.point(i, &mut buf);
                for (d, &x) in buf.iter().enumerate() {
                    bins[d][(x * f64::from(count)) as usize] += 1;
                }
            }
            for (d, b) in bins.iter().enumerate() {
                assert!(
                    b.iter().all(|&c| c == 1),
                    "dim {} not (0,m,1)-stratified at m={m}: {b:?}",
                    d + 1
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-rand/qmc\",\"case\":\"stratification\",\"verdict\":\"pass\",\"detail\":\"dims 1..=10, m 1..=8 exact\"}}"
        );
    }

    /// 2D elementary-interval property for the leading pair: the first 2^m
    /// points hit every 2^a × 2^b box (a + b = m) exactly once.
    #[test]
    fn leading_pair_elementary_intervals() {
        let s = Sobol::new(2);
        let mut buf = [0.0; 2];
        for m in 2..=8u32 {
            for a in 0..=m {
                let b = m - a;
                let (na, nb) = (1u32 << a, 1u32 << b);
                let mut boxes = vec![0u32; (na * nb) as usize];
                for i in 0..(1u32 << m) {
                    s.point(i, &mut buf);
                    let ia = (buf[0] * f64::from(na)) as u32;
                    let ib = (buf[1] * f64::from(nb)) as u32;
                    boxes[(ia * nb + ib) as usize] += 1;
                }
                assert!(
                    boxes.iter().all(|&c| c == 1),
                    "(a={a}, b={b}) elementary intervals violated"
                );
            }
        }
    }

    /// Owen scrambling must PRESERVE the net (stratification survives) while
    /// randomizing positions, and must replay bit-identically from its seed.
    #[test]
    fn owen_preserves_net_and_replays() {
        let plain = Sobol::new(4);
        let s1 = Sobol::scrambled(4, 0xA11CE);
        let s2 = Sobol::scrambled(4, 0xA11CE);
        let s3 = Sobol::scrambled(4, 0xB0B);
        let mut buf = vec![0.0; 4];
        // Stratification survives scrambling (nested uniform property).
        for m in 1..=7u32 {
            let count = 1u32 << m;
            let mut bins = vec![vec![0u32; count as usize]; 4];
            for i in 0..count {
                s1.point(i, &mut buf);
                for (d, &x) in buf.iter().enumerate() {
                    bins[d][(x * f64::from(count)) as usize] += 1;
                }
            }
            for b in &bins {
                assert!(
                    b.iter().all(|&c| c == 1),
                    "scrambling broke the net at m={m}"
                );
            }
        }
        // Replayable and seed-sensitive; differs from the plain net.
        let (mut a, mut b, mut c, mut p) = (vec![0.0; 4], vec![0.0; 4], vec![0.0; 4], vec![0.0; 4]);
        let mut differs_seed = false;
        let mut differs_plain = false;
        for i in 0..64 {
            s1.point(i, &mut a);
            s2.point(i, &mut b);
            s3.point(i, &mut c);
            plain.point(i, &mut p);
            assert!(
                a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()),
                "same seed must replay bitwise at point {i}"
            );
            differs_seed |= a.iter().zip(&c).any(|(x, y)| x.to_bits() != y.to_bits());
            differs_plain |= a.iter().zip(&p).any(|(x, y)| x.to_bits() != y.to_bits());
        }
        assert!(
            differs_seed && differs_plain,
            "scrambling must actually scramble"
        );
        println!(
            "{{\"suite\":\"fs-rand/qmc\",\"case\":\"owen\",\"verdict\":\"pass\",\"detail\":\"net preserved, replayable, seed-sensitive\"}}"
        );
    }

    /// The payoff test: on a smooth integrand, scrambled Sobol beats MC
    /// decisively at equal N (this is WHY the crate exists — plan §6.7's
    /// "1-2 orders of magnitude" claim, gated).
    #[test]
    fn qmc_beats_mc_on_genz_product_peak() {
        const DIM: usize = 5;
        // Genz product-peak: f(x) = Π 1/(c² + (x_j − w_j)²), analytic value
        // Π c·(atan(c·(1−w)) + atan(c·w)) ... use c=1, w=0.5 per dim:
        // ∫ 1/(1+(x−.5)²) dx = atan(.5) − atan(−.5) = 2·atan(.5).
        let f = |x: &[f64]| -> f64 {
            x.iter()
                .map(|&v| 1.0 / (1.0 + (v - 0.5) * (v - 0.5)))
                .product()
        };
        let exact = fs_math::det::powi(2.0 * 0.5f64.atan(), i32::try_from(DIM).expect("small"));
        let n = 4096u32;
        // Scrambled-Sobol RMSE over independent randomizations.
        let mut qmc_se = 0.0;
        for rep in 0..8u64 {
            let s = Sobol::scrambled(DIM, 0xC0DE + rep);
            let mut buf = vec![0.0; DIM];
            let mut acc = 0.0;
            for i in 0..n {
                s.point(i, &mut buf);
                acc += f(&buf);
            }
            let err = acc / f64::from(n) - exact;
            qmc_se += err * err;
        }
        let qmc_rmse = (qmc_se / 8.0).sqrt();
        // Plain MC RMSE over the same budget.
        let mut mc_se = 0.0;
        for rep in 0..8u32 {
            let mut st = crate::StreamKey {
                seed: 0xFACE,
                kernel: 9,
                tile: rep,
            }
            .stream();
            let mut acc = 0.0;
            let mut buf = vec![0.0; DIM];
            for _ in 0..n {
                st.fill_f64(&mut buf);
                acc += f(&buf);
            }
            let err = acc / f64::from(n) - exact;
            mc_se += err * err;
        }
        let mc_rmse = (mc_se / 8.0).sqrt();
        assert!(
            qmc_rmse * 5.0 < mc_rmse,
            "scrambled Sobol must beat MC decisively: qmc {qmc_rmse:.2e} vs mc {mc_rmse:.2e}"
        );
        println!(
            "{{\"suite\":\"fs-rand/qmc\",\"case\":\"genz\",\"verdict\":\"pass\",\"detail\":\"rmse qmc={qmc_rmse:.2e} mc={mc_rmse:.2e} at n={n}, dim={DIM}\"}}"
        );
    }

    /// CBC lattices: the Korobov worst-case error must fall near O(n⁻²) in
    /// the tabled range, and beat a bad (non-CBC) generating vector.
    #[test]
    fn cbc_lattice_error_decays_and_beats_naive() {
        let dims = 6;
        let e_small = Lattice::cbc(257, dims).korobov_error_sq();
        let e_big = Lattice::cbc(1031, dims).korobov_error_sq();
        // n grows ×4.01 → error² should drop by ≳ 4² (rate ~ n⁻²⁺ᵋ each in
        // error, i.e. error² ~ n⁻⁴⁺ᵋ; accept a lenient factor 8).
        assert!(
            e_big * 8.0 < e_small,
            "CBC error² must decay: {e_small:.3e} -> {e_big:.3e}"
        );
        // A deliberately poor vector (all components 1) must be worse.
        let naive = Lattice {
            n: 1031,
            z: vec![1; dims],
        };
        assert!(
            e_big * 4.0 < naive.korobov_error_sq(),
            "CBC must beat the naive vector"
        );
        println!(
            "{{\"suite\":\"fs-rand/qmc\",\"case\":\"cbc\",\"verdict\":\"pass\",\"detail\":\"err2 {e_small:.3e}@257 -> {e_big:.3e}@1031\"}}"
        );
    }

    #[test]
    fn baker_periodization_and_input_contracts() {
        assert_eq!(baker(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(baker(0.5).to_bits(), 1.0f64.to_bits());
        assert_eq!(baker(1.0).to_bits(), 0.0f64.to_bits());
        assert!((baker(0.25) - 0.5).abs() < 1e-15);
        // Contract violations refuse loudly.
        assert!(std::panic::catch_unwind(|| Sobol::new(0)).is_err());
        assert!(std::panic::catch_unwind(|| Sobol::new(MAX_SOBOL_DIM + 1)).is_err());
        assert!(std::panic::catch_unwind(|| Lattice::cbc(2, 3)).is_err());
    }
}
