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

/// Unsigned arbitrary-precision integer used only to rank CBC candidates.
///
/// The Bernoulli-B₂ kernel is rational at every lattice point. At a fixed
/// prefix length every candidate has the same denominator, so an exact sum of
/// numerator products determines the declared objective order. Base-2³² limbs
/// keep multiplication by the at-most-67-bit kernel numerator straightforward
/// and make overflow impossible except as an explicit allocation/capacity
/// failure. This avoids importing a general big-integer runtime dependency.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ExactNat {
    /// Little-endian base-2³² limbs; zero has no limbs.
    limbs: Vec<u32>,
}

impl ExactNat {
    const LIMB_BITS: u32 = 32;

    pub(crate) fn zero() -> Self {
        Self { limbs: Vec::new() }
    }

    pub(crate) fn one() -> Self {
        Self { limbs: vec![1] }
    }

    pub(crate) fn normalize(&mut self) {
        while self.limbs.last() == Some(&0) {
            let removed = self
                .limbs
                .pop()
                .expect("a zero limb observed by last() is present");
            debug_assert_eq!(removed, 0);
        }
    }

    pub(crate) fn magnitude_cmp(&self, other: &Self) -> core::cmp::Ordering {
        match self.limbs.len().cmp(&other.limbs.len()) {
            core::cmp::Ordering::Equal => self.limbs.iter().rev().cmp(other.limbs.iter().rev()),
            ordering => ordering,
        }
    }

    /// Add `multiplicand * factor` exactly using base-2³² schoolbook
    /// arithmetic. The score accumulator deliberately retains spare zero
    /// limbs between additions and normalizes once before comparison.
    pub(crate) fn add_mul_factor(&mut self, multiplicand: &Self, factor: u128) {
        if factor == 0 || multiplicand.limbs.is_empty() {
            return;
        }

        let mut factor_limbs = [0_u32; 4];
        let mut remaining = factor;
        let mut factor_len = 0_usize;
        while remaining != 0 {
            factor_limbs[factor_len] = u32::try_from(remaining & u128::from(u32::MAX))
                .expect("a masked base-2^32 limb fits u32");
            remaining >>= Self::LIMB_BITS;
            factor_len += 1;
        }

        let needed = multiplicand
            .limbs
            .len()
            .checked_add(factor_len)
            .and_then(|length| length.checked_add(1))
            .expect("exact CBC accumulator limb capacity overflow");
        if self.limbs.len() < needed {
            self.limbs.resize(needed, 0);
        }

        for (source_index, &source_limb) in multiplicand.limbs.iter().enumerate() {
            let mut carry = 0_u64;
            for (factor_index, &factor_limb) in factor_limbs[..factor_len].iter().enumerate() {
                let destination = source_index + factor_index;
                let wide = u64::from(source_limb)
                    .checked_mul(u64::from(factor_limb))
                    .and_then(|value| value.checked_add(u64::from(self.limbs[destination])))
                    .and_then(|value| value.checked_add(carry))
                    .expect("base-2^32 multiply-add fits u64");
                self.limbs[destination] = u32::try_from(wide & u64::from(u32::MAX))
                    .expect("a masked base-2^32 limb fits u32");
                carry = wide >> Self::LIMB_BITS;
            }

            let mut destination = source_index + factor_len;
            while carry != 0 {
                if destination == self.limbs.len() {
                    self.limbs.push(0);
                }
                let wide = u64::from(self.limbs[destination])
                    .checked_add(carry)
                    .expect("base-2^32 carry propagation fits u64");
                self.limbs[destination] = u32::try_from(wide & u64::from(u32::MAX))
                    .expect("a masked base-2^32 limb fits u32");
                carry = wide >> Self::LIMB_BITS;
                destination += 1;
            }
        }
    }

    /// Add under an admission-owned requested-capacity ceiling.
    ///
    /// Returns the required limb count before touching the accumulator when
    /// the shared storage schema is too small.
    pub(crate) fn add_mul_factor_with_capacity(
        &mut self,
        multiplicand: &Self,
        factor: u128,
        capacity_limbs: usize,
    ) -> Result<(), usize> {
        if factor == 0 || multiplicand.limbs.is_empty() {
            return Ok(());
        }
        let required = multiplicand
            .limbs
            .len()
            .checked_add(factor_limb_count(factor))
            .and_then(|length| length.checked_add(1))
            .expect("exact CBC required limb count overflow");
        if required > capacity_limbs {
            return Err(required);
        }
        self.add_mul_factor(multiplicand, factor);
        debug_assert!(self.limbs.len() <= capacity_limbs);
        Ok(())
    }

    pub(crate) fn mul_assign_factor(&mut self, factor: u128) {
        let multiplicand = Self {
            limbs: core::mem::take(&mut self.limbs),
        };
        self.add_mul_factor(&multiplicand, factor);
        self.normalize();
    }

    /// Multiply while requesting exactly the admission-owned replacement
    /// capacity and retaining the moved old allocation until completion.
    pub(crate) fn mul_assign_factor_with_capacity(
        &mut self,
        factor: u128,
        capacity_limbs: usize,
    ) -> Result<(), usize> {
        if factor == 0 || self.limbs.is_empty() {
            self.limbs.clear();
            return Ok(());
        }
        let required = self
            .limbs
            .len()
            .checked_add(factor_limb_count(factor))
            .and_then(|length| length.checked_add(1))
            .expect("exact CBC required limb count overflow");
        if required > capacity_limbs {
            return Err(required);
        }
        let multiplicand = Self {
            limbs: core::mem::take(&mut self.limbs),
        };
        self.limbs = Vec::with_capacity(capacity_limbs);
        self.add_mul_factor(&multiplicand, factor);
        self.normalize();
        debug_assert!(self.limbs.len() <= capacity_limbs);
        Ok(())
    }

    /// Reserve exact limb capacity so subsequent in-envelope arithmetic never
    /// reallocates (the admission contract's exact-capacity storage
    /// requirement for the execution tranche).
    pub(crate) fn reserve_exact_limbs(&mut self, capacity_limbs: usize) {
        let additional = capacity_limbs.saturating_sub(self.limbs.len());
        self.limbs.reserve_exact(additional);
    }

    /// The normalized little-endian limbs (certificate serialization and
    /// independent checking; callers must normalize first).
    pub(crate) fn limbs(&self) -> &[u32] {
        &self.limbs
    }

    #[cfg(test)]
    fn from_u128(mut value: u128) -> Self {
        let mut result = Self::zero();
        while value != 0 {
            result.limbs.push(
                u32::try_from(value & u128::from(u32::MAX))
                    .expect("a masked base-2^32 limb fits u32"),
            );
            value >>= Self::LIMB_BITS;
        }
        result
    }

    #[cfg(test)]
    fn to_u128(&self) -> u128 {
        assert!(
            self.limbs.len() <= 4,
            "test conversion only covers values representable by u128"
        );
        self.limbs.iter().rev().fold(0_u128, |value, &limb| {
            (value << Self::LIMB_BITS) | u128::from(limb)
        })
    }
}

fn factor_limb_count(mut factor: u128) -> usize {
    let mut limbs = 0_usize;
    while factor != 0 {
        limbs += 1;
        factor >>= ExactNat::LIMB_BITS;
    }
    limbs
}

/// Integer numerator of `1 + B₂(residue/n)` over the candidate-independent
/// denominator `6*n²`.
pub(crate) fn exact_kernel_numerator(n: u32, residue: u32) -> u128 {
    let n = u128::from(n);
    let residue = u128::from(residue);
    let positive = 7 * n * n + 6 * residue * residue;
    positive
        .checked_sub(6 * residue * n)
        .expect("the Bernoulli-B2 kernel numerator is non-negative")
}

pub(crate) fn lattice_residue(point: usize, generator: u32, n: u32) -> u32 {
    let point = u64::try_from(point).expect("a u32-bounded lattice index fits u64");
    u32::try_from(point * u64::from(generator) % u64::from(n))
        .expect("a modular residue below n fits u32")
}

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
    /// previously fixed components. Candidate ordering is exact: the rational
    /// kernel's common denominator is removed and arbitrary-precision integer
    /// numerator sums are compared, with the lower candidate winning exact
    /// ties. There are O(n²·dim) lattice visits; exact-integer limb work grows
    /// with prefix length.
    ///
    /// # Panics
    /// If `n < 3` or `dim` is 0.
    #[must_use]
    pub fn cbc(n: u32, dim: usize) -> Lattice {
        assert!(n >= 3 && dim >= 1, "lattice needs n >= 3, dim >= 1");
        let point_count = usize::try_from(n).expect("lattice point count fits usize");
        // products[k] is the exact numerator product for the chosen prefix at
        // point k. The omitted denominator (6*n²)^prefix_len is common to
        // every candidate in one CBC step.
        let mut products = vec![ExactNat::one(); point_count];
        let mut z = Vec::with_capacity(dim);

        // Every unit modulo n permutes the complete residue set, so every
        // admissible first component has exactly the same one-dimensional
        // score. Candidate 1 is therefore the mathematical lowest-candidate
        // tie winner for every supported n, not a fixture-specific shortcut.
        let first = 1_u32;
        for (point, product) in products.iter_mut().enumerate() {
            let residue = lattice_residue(point, first, n);
            product.mul_assign_factor(exact_kernel_numerator(n, residue));
        }
        z.push(first);

        for _ in 1..dim {
            let mut best: Option<(ExactNat, u32)> = None;
            for candidate in 1..n {
                if gcd(candidate, n) != 1 {
                    continue;
                }
                let mut score = ExactNat::zero();
                for (point, prefix_product) in products.iter().enumerate() {
                    let residue = lattice_residue(point, candidate, n);
                    score.add_mul_factor(prefix_product, exact_kernel_numerator(n, residue));
                }
                score.normalize();

                let replace = match &best {
                    None => true,
                    Some((best_score, best_candidate)) => match score.magnitude_cmp(best_score) {
                        core::cmp::Ordering::Less => true,
                        core::cmp::Ordering::Equal => candidate < *best_candidate,
                        core::cmp::Ordering::Greater => false,
                    },
                };
                if replace {
                    best = Some((score, candidate));
                }
            }

            let (_, chosen) = best.expect("candidate 1 is coprime to every n");
            for (point, product) in products.iter_mut().enumerate() {
                let residue = lattice_residue(point, chosen, n);
                product.mul_assign_factor(exact_kernel_numerator(n, residue));
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

pub(crate) fn gcd(a: u32, b: u32) -> u32 {
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

    fn exact_prefix_score_for_test(n: u32, generators: &[u32]) -> ExactNat {
        let point_count = usize::try_from(n).expect("lattice point count fits usize");
        let mut score = ExactNat::zero();
        for point in 0..point_count {
            let mut product = ExactNat::one();
            for &generator in generators {
                let residue = lattice_residue(point, generator, n);
                product.mul_assign_factor(exact_kernel_numerator(n, residue));
            }
            score.add_mul_factor(&product, 1);
        }
        score.normalize();
        score
    }

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
    fn exact_nat_multiply_add_carries_match_u128() {
        let max_kernel_factor = 7 * u128::from(u32::MAX) * u128::from(u32::MAX);
        let cases = [
            (0_u128, 0_u128, max_kernel_factor),
            (0, u128::from(u32::MAX), max_kernel_factor),
            (u128::from(u64::MAX), (1_u128 << 64) + 1, (1_u128 << 32) + 1),
            (
                (1_u128 << 100) + 17,
                (1_u128 << 63) + 29,
                (1_u128 << 31) + 3,
            ),
        ];

        for (accumulator, multiplicand, factor) in cases {
            let expected = accumulator
                .checked_add(
                    multiplicand
                        .checked_mul(factor)
                        .expect("fixture product fits u128"),
                )
                .expect("fixture multiply-add fits u128");
            let mut computed = ExactNat::from_u128(accumulator);
            computed.add_mul_factor(&ExactNat::from_u128(multiplicand), factor);
            computed.normalize();
            assert_eq!(computed.to_u128(), expected);

            if accumulator == 0 {
                let mut assigned = ExactNat::from_u128(multiplicand);
                assigned.mul_assign_factor(factor);
                assert_eq!(assigned.to_u128(), expected);
            }
        }

        // Independent base-2^32 KAT beyond u128, so the production-sized
        // multi-limb path cannot self-confirm through the bounded conversion.
        let mut wide = ExactNat::from_u128(u128::MAX);
        wide.mul_assign_factor(max_kernel_factor);
        assert_eq!(
            wide.limbs,
            [
                4_294_967_289,
                13,
                4_294_967_289,
                4_294_967_295,
                6,
                4_294_967_282,
                6
            ]
        );
        wide.add_mul_factor(
            &ExactNat::from_u128((1_u128 << 127) + 12_345),
            (1_u128 << 64) + 17,
        );
        wide.normalize();
        assert_eq!(
            wide.limbs,
            [209_858, 14, 12_338, 2_147_483_648, 15, 2_147_483_634, 7]
        );
    }

    /// G0 finite property sweep: multiplication by a unit modulo n merely
    /// permutes the complete residue set, so every admissible first component
    /// has the same exact score and CBC must choose the lowest one.
    #[test]
    fn cbc_first_component_units_are_exact_ties_and_choose_one() {
        for n in 3..=64_u32 {
            let reference = exact_prefix_score_for_test(n, &[1]);
            for candidate in 1..n {
                if gcd(candidate, n) == 1 {
                    assert_eq!(
                        exact_prefix_score_for_test(n, &[candidate]),
                        reference,
                        "unit {candidate} must only permute residues for n={n}"
                    );
                }
            }
            assert_eq!(Lattice::cbc(n, 1).z, [1]);
        }
    }

    #[test]
    fn cbc_resolves_later_exact_ties_and_matches_independent_kat() {
        let candidate_one = exact_prefix_score_for_test(5, &[1, 2, 1, 2, 1]);
        let candidate_two = exact_prefix_score_for_test(5, &[1, 2, 1, 2, 2]);
        assert_eq!(
            candidate_one, candidate_two,
            "the later tie fixture is exact"
        );
        assert_eq!(Lattice::cbc(5, 6).z, [1, 2, 1, 2, 1, 2]);

        // Independently enumerated exact-integer KAT. This size is large
        // enough to exercise multi-limb prefix scores without duplicating the
        // n=257 replay receipt in the integration Casebook.
        assert_eq!(Lattice::cbc(127, 5).z, [1, 29, 24, 56, 35]);
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
